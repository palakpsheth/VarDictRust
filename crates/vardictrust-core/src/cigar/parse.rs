//! CIGAR parse (pileup engine) — ports CigarParser.java.

use crate::bam::{get_mate_reference_name, SamRecord};
use crate::config::Config;
use crate::reference::Reference;
use crate::regions::Region;
use crate::types::VardictError;
use crate::variations::{
    get_variation, get_variation_maybe_mut as get_variation_maybe, Mate, SVStructures, Sclip,
    Variation, VariationData, VariationMap,
};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// Matches `MC_CLIP_CLIP` (`\d+S.*?\d+S`) without regex compilation or backtracking.
#[inline]
fn has_two_soft_clips(mc: &str) -> bool {
    let mut count = 0u8;
    let mut prev_digit = false;
    for &b in mc.as_bytes() {
        if b == b'\n' {
            count = 0;
            prev_digit = false;
        } else if b == b'S' && prev_digit {
            count += 1;
            if count >= 2 {
                return true;
            }
            prev_digit = false;
        } else {
            prev_digit = b.is_ascii_digit();
        }
    }
    false
}

/// CigarParser stores the running state for parsing CIGAR strings and building pileups.
/// Ports `modules/CigarParser.java`.
pub struct CigarParser<'a> {
    pub reference: &'a Reference,
    pub region: &'a Region,
    pub config: &'a Config,
    pub svflag: bool,
    pub(crate) primary_ref_seq: &'a [u8],
    pub(crate) primary_ref_start: u32,
    pub(crate) cov_cap: i64,

    pub non_insertion_variants: crate::javamap::JavaHashMap<i32, VariationMap>,
    pub insertion_variants: crate::javamap::JavaHashMap<i32, VariationMap>,
    // Iteration-order audit (W5.3, 2026-08-23): the String-keyed inner HashMaps
    // below never leak order into output or realignment decisions.
    // position_to_insertion_count/position_to_deletion_count/mnp are consumed
    // ONLY via fill_and_sort_tmp, whose (count,pos,desc) comparator is a total
    // order; splice_count DOES reach stdout (lib.rs -x print) but its VarDictJava
    // counterpart is also a HashMap iterated unordered (CigarParser.java:65,
    // :102-106), so record-set parity is unaffected. All four are `new HashMap<>()`
    // in CigarParser.java:62-65. SVStructures svffus/svrfus (variations/mod.rs)
    // use IndexMap + java_hashmap_key_order because is_best_sv_used (sv.rs)
    // tie-breaks across mate contigs in Java HashMap<String, List<Sclip>> order.
    pub position_to_insertion_count: BTreeMap<i32, rustc_hash::FxHashMap<String, i32>>,
    pub position_to_deletion_count: BTreeMap<i32, rustc_hash::FxHashMap<String, i32>>,
    pub sv_structures: SVStructures,
    pub ref_coverage: rustc_hash::FxHashMap<i32, i32>,
    pub soft_clips_5_end: BTreeMap<i32, Sclip>,
    pub soft_clips_3_end: BTreeMap<i32, Sclip>,
    pub max_read_length: i32,
    pub splice: rustc_hash::FxHashSet<String>,
    pub mnp: BTreeMap<i32, rustc_hash::FxHashMap<String, i32>>,
    pub splice_count: rustc_hash::FxHashMap<String, [i32; 1]>,
    pub discordant_count: i32,
    pub total_reads: i32,
    pub duplicate_reads: i32,
    /// Per-read scratch buffers reused across [`CigarParser::parse_cigar`]
    /// calls (zero per-read allocation after warmup). Never read across
    /// calls: refilled from each record before use.
    pub scratch_seq: Vec<u8>,
    pub scratch_qual: Vec<u8>,
    pub scratch_cigar_str: String,
    pub scratch_elements: Vec<(u32, CigarKind)>,
    /// Per-CIGAR-Match-block string scratch, reused across reads via
    /// per-read take/put-back (same pattern as the sequence buffers).
    pub scratch_s: String,
    pub scratch_ss: String,
    /// Splice-junction key buffer (F3 churn fix).
    pub scratch_splice_key: String,
    /// F4 scratch: per-indel-occurrence desc/quality buffers, take/put-back
    /// per process_insertion/process_deletion call (methods never nest, so
    /// one set suffices).
    pub scratch_f4_desc: String,
    pub scratch_f4_qual: Vec<u8>,
    pub scratch_f4_seg1: String,
    pub scratch_f4_seg2: Vec<u8>,
}

/// Zero-alloc single-base allele key for per-base variation lookups.
/// Returns a `&'static str` for ACGTN (either case — the overwhelmingly common
/// reference bytes); exotic bytes fall back to the caller-provided buffer.
/// Identical keys to `(base as char).to_string()`, without the per-base heap
/// allocation that dominates soft-clip walk-backs on deep data.
#[inline]
fn allele_key<'a>(base: u8, buf: &'a mut [u8; 4]) -> &'a str {
    match base {
        b'A' => "A",
        b'C' => "C",
        b'G' => "G",
        b'T' => "T",
        b'N' => "N",
        b'n' => "n",
        b'a' => "a",
        b'c' => "c",
        b'g' => "g",
        b't' => "t",
        _ => (base as char).encode_utf8(buf),
    }
}

#[inline]
fn push_ascii_base(dst: &mut String, b: u8) {
    dst.push(if b.is_ascii() { b as char } else { 'N' });
}

#[inline]
fn push_ascii_slice(dst: &mut String, bytes: &[u8]) {
    dst.reserve(bytes.len());
    for &b in bytes {
        push_ascii_base(dst, b);
    }
}

/// Static single-base `&str` lookup matching [`push_ascii_base`] + `N`/`n` filtering:
/// ASCII `0x00..=0x7F` map to their 1-byte `&str` except `b'n'` (`0x6E`), which maps
/// to `"N"` alongside `b'N'` and non-ASCII `0x80..=0xFF`.
const SINGLE_BASE_STR: [&str; 256] = [
    "\x00", "\x01", "\x02", "\x03", "\x04", "\x05", "\x06", "\x07", "\x08", "\t", "\n", "\x0B",
    "\x0C", "\r", "\x0E", "\x0F", "\x10", "\x11", "\x12", "\x13", "\x14", "\x15", "\x16", "\x17",
    "\x18", "\x19", "\x1A", "\x1B", "\x1C", "\x1D", "\x1E", "\x1F", " ", "!", "\"", "#", "$", "%",
    "&", "'", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8",
    "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
    "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^",
    "_", "`", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "N", "o", "p", "q",
    "r", "s", "t", "u", "v", "w", "x", "y", "z", "{", "|", "}", "~", "\x7F", "N", "N", "N", "N",
    "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N",
    "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N",
    "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N",
    "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N",
    "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N",
    "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N",
    "N", "N", "N", "N", "N", "N", "N", "N", "N", "N",
];

/// Phred quality at `idx`, defaulting missing/truncated QUAL (`*`) to 0.
/// Ports Java leniency: `getBaseQualityString` may yield `*`/short string where
/// `charAt(i)-33` would throw; Rust must not crash (AGENTS §11) — fail the
/// `base_qual` filter instead. `shortcut: O(1) bounds check, no upgrade needed`.
#[inline]
fn qual_phred<Q: AsRef<[u8]>>(query_qual: &Q, idx: usize) -> f64 {
    query_qual
        .as_ref()
        .get(idx)
        .copied()
        .unwrap_or(33)
        .saturating_sub(33) as f64
}

impl<'a> CigarParser<'a> {
    /// Create a new CigarParser with given reference, region, and config.
    /// Ports `CigarParser.java:61-68`.
    pub fn new(
        reference: &'a Reference,
        region: &'a Region,
        config: &'a Config,
        svflag: bool,
        total_reads: i32,
        duplicate_reads: i32,
        seed_max_read_length: i32,
    ) -> Self {
        let (primary_ref_seq, primary_ref_start) = match (
            reference.loaded_regions.first(),
            reference.region_sequences.first(),
        ) {
            (Some(lr), Some(bases))
                if lr.chr == region.chr && lr.sequence_end >= lr.sequence_start =>
            {
                let max_len = (lr.sequence_end.saturating_sub(lr.sequence_start) as usize)
                    .saturating_add(1)
                    .min(bases.len());
                (&bases[..max_len], lr.sequence_start)
            }
            _ => (&[][..], 1),
        };
        let cov_cap = match reference
            .loaded_regions
            .iter()
            .filter(|lr| lr.chr == region.chr)
            .map(|lr| lr.sequence_end)
            .max()
        {
            Some(max_end) => i64::from(max_end.max(region.end.0)) + 1024,
            None => i64::from(region.end.0) + i64::from(config.max_region_len) + 65_536,
        }
        .min(i64::from(i32::MAX));
        Self {
            reference,
            region,
            config,
            svflag,
            primary_ref_seq,
            primary_ref_start,
            cov_cap,
            non_insertion_variants: crate::javamap::JavaHashMap::new(),
            insertion_variants: crate::javamap::JavaHashMap::new(),
            position_to_insertion_count: BTreeMap::new(),
            position_to_deletion_count: BTreeMap::new(),
            sv_structures: SVStructures::new(),
            ref_coverage: rustc_hash::FxHashMap::default(),
            soft_clips_5_end: BTreeMap::new(),
            soft_clips_3_end: BTreeMap::new(),
            max_read_length: seed_max_read_length,
            splice: rustc_hash::FxHashSet::default(),
            mnp: BTreeMap::new(),
            splice_count: rustc_hash::FxHashMap::default(),
            discordant_count: 0,
            total_reads,
            duplicate_reads,
            scratch_seq: Vec::new(),
            scratch_qual: Vec::new(),
            scratch_cigar_str: String::new(),
            scratch_elements: Vec::new(),
            scratch_s: String::new(),
            scratch_ss: String::new(),
            scratch_splice_key: String::new(),
            scratch_f4_desc: String::new(),
            scratch_f4_qual: Vec::new(),
            scratch_f4_seg1: String::new(),
            scratch_f4_seg2: Vec::new(),
        }
    }

    #[inline]
    fn get_ref_base(&self, pos: u32) -> Option<u8> {
        if let Some(offset) = pos.checked_sub(self.primary_ref_start) {
            if let Some(&b) = self.primary_ref_seq.get(offset as usize) {
                return Some(b);
            }
        }
        self.reference.get_base(&self.region.chr, pos)
    }

    #[inline]
    fn ref_mismatch_not_n(&self, pos: u32, query_base: Option<u8>) -> bool {
        match (self.get_ref_base(pos), query_base) {
            (Some(ref_base), Some(qb)) => ref_base != qb && ref_base != b'N' && ref_base != b'n',
            _ => false,
        }
    }

    /// Parse a single read record and update the pileup counts.
    /// Ports `CigarParser.java:216-310,332-621`.
    pub fn parse_cigar(&mut self, record: &SamRecord) -> Result<(), VardictError> {
        // Drop reads that do not overlap the target region
        let read_start = record.alignment_start.map(|p| p.0 as i32).unwrap_or(0);
        let mut read_end = read_start;
        for &(len, op) in &record.cigar {
            if matches!(
                op,
                CigarKind::Match
                    | CigarKind::Deletion
                    | CigarKind::Skip
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
            ) {
                read_end = read_end.saturating_add(i32::try_from(len).unwrap_or(i32::MAX));
            }
        }
        let region_start = self.region.start.0 as i32;
        let region_end = self.region.end.0 as i32;
        if read_start > region_end || read_end <= region_start {
            return Ok(());
        }

        let direction = (record.flags & 0x10) != 0;
        let insertion_deletion_length = get_insertion_deletion_length(&record.cigar);
        let adjusted_nm = record.nm.map_or(0, |n| {
            n.saturating_sub(i32::try_from(insertion_deletion_length).unwrap_or(i32::MAX))
        });

        // T045 mismatch gate: drop read if NM - (I + D) > max_mismatches (CigarParser.java:244-251)
        if record.nm.is_some() {
            if adjusted_nm > self.config.max_mismatches as i32 {
                return Ok(()); // Drop read
            }
        } else if record.cigar.is_empty() || (record.flags & 0x4) != 0 {
            return Ok(());
        }
        if record.cigar.is_empty() {
            return Ok(());
        }

        let mut start_pos = record.alignment_start.map(|p| p.0 as i32).unwrap_or(0);
        // Amplicon read filter (Java parseCigarWithAmpCase, CigarParser.java:255-258):
        // skip reads that don't map to the amplicon. Runs before scratch
        // fill so dropped reads allocate nothing (cigar_elements still
        // equals record.cigar here — nothing has mutated it yet).
        if let Some(ref amp_param) = self.config.amplicon {
            let (distance_to_amplicon, overlap_fraction) = parse_amp_param(amp_param);
            let is_mate_ref_name_equal =
                record.reference_name.as_deref() == record.mate_reference_name.as_deref();
            if parse_cigar_with_amp_case(
                &record.cigar,
                record,
                is_mate_ref_name_equal,
                self.region,
                distance_to_amplicon,
                overlap_fraction,
            ) {
                return Ok(());
            }
        }

        if !self.svflag && self.ref_coverage.capacity() == 0 {
            let span_cap = (self
                .region
                .end
                .0
                .saturating_sub(self.region.start.0)
                .saturating_add(1) as usize)
                .min(4_096);
            self.ref_coverage.reserve(span_cap);
            self.non_insertion_variants.reserve_nodes(span_cap);
        }

        // Zero-copy read fast path: borrow `&record.sequence` and `&record.quality_scores`
        // directly unless `CigarModifier` fallback actually runs and materializes buffers.
        self.scratch_elements.clear();
        self.scratch_elements.extend_from_slice(&record.cigar);
        let mut cigar_elements = std::mem::take(&mut self.scratch_elements);
        let mut modified_bufs: Option<(Vec<u8>, Vec<u8>)> = None;

        // T046: run CigarModifier if realignment is enabled
        // Ports CigarParser.java:257 — Java's `reference.referenceSequences` is a
        // single per-chromosome Map<Integer, Character> shared across all loaded
        // sub-regions (including secondary regions loaded by SV `partialPipeline`).
        // Select the loaded region on `self.region.chr` whose stored sequence covers
        // the read span `[start_pos, start_pos + query_seq.len()]` (falling back to
        // any overlapping region, then the first region on `self.region.chr`).
        if self.config.realign == 1 {
            let qlen = i32::try_from(record.sequence.len()).unwrap_or(i32::MAX);
            let span_end = read_end.max(read_start.saturating_add(qlen));
            let ext_start = read_start.saturating_sub(qlen);
            let ext_end = span_end.saturating_add(qlen);
            let (ref_slice, ref_start) = self.reference.region_slice_for_span(
                &self.region.chr,
                read_start,
                span_end,
                ext_start,
                ext_end,
            );

            // ports CigarModifier.java:619-791, CigarParser.java:80-95:
            // When modify_cigar returns an error (JavaThrows on malformed CIGAR with non-positive
            // operator lengths), emulate Java's catch block: log the exception and skip the read
            // cleanly without mutating start_pos, query_seq, or pileup structures, restoring all
            // parser-owned scratch buffers for the next read.
            let log_modify_err = |e: crate::types::VardictError| {
                let name = String::from_utf8_lossy(&record.name);
                log::warn!(
                    "There was Exception while processing record on {}. The processing will be continued from the next record. Exception: {}",
                    name, e
                );
                // ports Utils.java:261-277 and CLAUDE.md §5:
                // Count exceptions atomically. While Java terminates with exit(1) after
                // MAX_EXCEPTION_COUNT (10) continued exceptions, VarDictRust never aborts;
                // it logs a warning and continues processing valid reads.
                let count = self
                    .config
                    .exception_counter
                    .fetch_add(1, std::sync::atomic::Ordering::Relaxed)
                    + 1;
                if count > 10 {
                    log::warn!("VarDictRust encountered >10 exceptions. Java would abort here, but we continue (CLAUDE.md §5).");
                }
            };
            let run_fallback_modifier = match *cigar_elements.as_slice() {
                [(mch_u32, CigarKind::Match)] if i32::try_from(mch_u32).is_ok_and(|m| m > 0) => {
                    if let Err(e) = crate::cigar::modify::modify_single_match_in_place(
                        &mut start_pos,
                        mch_u32 as i32,
                        &mut cigar_elements,
                        &record.sequence,
                        ref_slice,
                        ref_start,
                    ) {
                        self.scratch_elements = cigar_elements;
                        log_modify_err(e);
                        return Ok(());
                    }
                    false
                }
                _ => {
                    let chr_len = if ref_slice.is_empty() {
                        0
                    } else {
                        (i64::from(ref_start) + (ref_slice.len() as i64) - 1)
                            .clamp(0, i64::from(u32::MAX)) as u32
                    };
                    match crate::cigar::modify::modify_pure_sm_in_place(
                        &mut start_pos,
                        &mut cigar_elements,
                        &record.sequence,
                        &record.quality_scores,
                        ref_slice,
                        ref_start,
                        chr_len,
                        self.config,
                    ) {
                        Ok(crate::cigar::modify::PureSmOutcome::Handled) => false,
                        Ok(crate::cigar::modify::PureSmOutcome::Fallback) => true,
                        Err(e) => {
                            self.scratch_elements = cigar_elements;
                            log_modify_err(e);
                            return Ok(());
                        }
                    }
                }
            };

            if run_fallback_modifier {
                self.scratch_seq.clear();
                self.scratch_seq.extend_from_slice(&record.sequence);
                self.scratch_qual.clear();
                self.scratch_qual.extend_from_slice(&record.quality_scores);
                let query_seq_buf = std::mem::take(&mut self.scratch_seq);
                let query_qual_buf = std::mem::take(&mut self.scratch_qual);
                let mut cigar_str = std::mem::take(&mut self.scratch_cigar_str);
                cigar_str.clear();
                crate::bam::push_cigar_string(&mut cigar_str, &cigar_elements);
                let mut modifier = crate::cigar::modify::CigarModifier::new(
                    start_pos,
                    cigar_str,
                    query_seq_buf,
                    query_qual_buf,
                    ref_slice,
                    ref_start,
                    &self.reference.seed,
                    i32::try_from(insertion_deletion_length).unwrap_or(i32::MAX),
                    self.max_read_length,
                    self.config,
                );
                match modifier.modify_cigar_preserving_buffers() {
                    Ok(modified) => {
                        start_pos = modified.position;
                        if crate::cigar::parse_cigar_into(&modified.cigar, &mut cigar_elements)
                            .is_none()
                        {
                            // Unreachable in practice (modify emits normalized
                            // CIGARs), but replicate the old keep-previous shape
                            // exactly: refill from the record.
                            cigar_elements.clear();
                            cigar_elements.extend_from_slice(&record.cigar);
                        }
                        modified_bufs = Some((modified.query_seq, modified.query_qual));
                        self.scratch_cigar_str = modified.cigar;
                    }
                    Err(e) => {
                        let (cigar_str, query_seq_buf, query_qual_buf) = modifier.into_buffers();
                        self.scratch_cigar_str = cigar_str;
                        self.scratch_seq = query_seq_buf;
                        self.scratch_qual = query_qual_buf;
                        self.scratch_elements = cigar_elements;
                        log_modify_err(e);
                        return Ok(());
                    }
                }
            }
        }

        // Determine discordant reads
        if !crate::bam::mate_reference_is_equal(record) {
            self.discordant_count += 1;
        }

        // T046: Ignore reads that are softclipped at both ends and both greater than 10 bp
        // Java regex ^\d\dS.*\d\dS$ means: first soft clip exactly 2 digits (10-99),
        // last soft clip at least 2 digits (>= 10)
        if cigar_elements.len() >= 2 {
            if let (
                Some(&(len_first, CigarKind::SoftClip)),
                Some(&(len_last, CigarKind::SoftClip)),
            ) = (cigar_elements.first(), cigar_elements.last())
            {
                if (10..=99).contains(&len_first) && len_last >= 10 {
                    // F8 churn fix: put back scratch before dropping the read.
                    if let Some((s, q)) = modified_bufs {
                        self.scratch_seq = s;
                        self.scratch_qual = q;
                    }
                    self.scratch_elements = cigar_elements;
                    return Ok(()); // Drop read
                }
            }
        }

        let read_len_match_ins = get_match_insertion_length(&cigar_elements);
        if self.config.minmatch != 0 && read_len_match_ins < self.config.minmatch {
            if let Some((s, q)) = modified_bufs {
                self.scratch_seq = s;
                self.scratch_qual = q;
            }
            self.scratch_elements = cigar_elements;
            return Ok(()); // Drop read shorter than minmatch (Java -M)
        }
        let read_len_match_ins_i32 = i32::try_from(read_len_match_ins).unwrap_or(i32::MAX);
        let total_length_including_soft_clipped =
            i32::try_from(get_soft_clipped_length(&cigar_elements)).unwrap_or(i32::MAX);
        if total_length_including_soft_clipped > self.max_read_length {
            self.max_read_length = total_length_including_soft_clipped;
        }

        // If supplementary alignment is present and flag filter is enabled, ignore it
        // Ports CigarParser.java:307-315: supplementary check runs AFTER updating maxReadLength.
        if self.config.skip_secondary_and_supplementary() && (record.flags & 0x800) != 0 {
            if let Some((s, q)) = modified_bufs {
                self.scratch_seq = s;
                self.scratch_qual = q;
            }
            self.scratch_elements = cigar_elements;
            return Ok(());
        }

        // CRISPR mode: skip reads not overlapping the cutting site region
        // Ports CigarParser.java:318, 634-646 (`skipSitesOutRegionOfInterest`).
        if let Some(cut_site) = self.config.crispr_cutting_site {
            if cut_site != 0 {
                if let Some(filter_bp) = self.config.crispr_flanking_size {
                    if filter_bp != 0 {
                        let rlen3 = i64::from(get_aligned_length_amp(&cigar_elements));
                        let cut_site = i64::from(cut_site);
                        let start_pos_64 = i64::from(start_pos);
                        let filter_bp = i64::from(filter_bp);
                        if !(cut_site - start_pos_64 > filter_bp
                            && start_pos_64 + rlen3 - cut_site > filter_bp)
                        {
                            if let Some((s, q)) = modified_bufs {
                                self.scratch_seq = s;
                                self.scratch_qual = q;
                            }
                            self.scratch_elements = cigar_elements;
                            return Ok(());
                        }
                    }
                }
            }
        }

        let (query_seq, query_qual): (&[u8], &[u8]) = match &modified_bufs {
            Some((s, q)) => (s.as_slice(), q.as_slice()),
            None => (record.sequence.as_slice(), record.quality_scores.as_slice()),
        };

        let paired = (record.flags & 0x1) != 0;
        let mate_unmapped = (record.flags & 0x8) != 0;
        let mapping_quality = record.mapping_quality.unwrap_or(0);
        let mapping_quality_i32 = i32::from(mapping_quality);
        if !(paired && mate_unmapped) {
            if mapping_quality > 10 && !self.config.no_sv {
                self.prepare_sv_structures_for_analysis(
                    record,
                    query_qual,
                    adjusted_nm,
                    direction,
                    (record.flags & 0x20) == 0,
                    start_pos,
                    total_length_including_soft_clipped,
                    &cigar_elements,
                );
            }
        }

        // Main CIGAR parser / pileup loop
        let mut read_pos_inc: usize = 0;
        let mut read_pos_exc: usize = 0;
        let mut start = start_pos;
        let mut offset = 0;
        let mut ci = 0;

        let mate_alignment_start = record.mate_alignment_start.map(|p| p.0 as i32).unwrap_or(0);
        let paired_same_chr = (self.config.unique_mode_alignment_enabled
            || self.config.unique_mode_second_in_pair_enabled)
            && is_paired_and_same_chromosome(record);
        let overlap_u = self.config.unique_mode_alignment_enabled && !direction && paired_same_chr;
        let overlap_second = self.config.unique_mode_second_in_pair_enabled
            && (record.flags & 0x80) != 0
            && paired_same_chr;
        let overlap_second_end_i64: Option<i64> = if !overlap_second {
            None
        } else if start_pos >= mate_alignment_start {
            let ref_len = i64::from(get_reference_length(&record.cigar));
            Some(i64::from(mate_alignment_start) + ref_len - 1)
        } else {
            let alignment_end = record.alignment_end.map(|p| p.0 as i32).unwrap_or(0);
            if mate_alignment_start <= alignment_end {
                Some(i64::MAX)
            } else {
                None
            }
        };
        let has_overlap_gate = overlap_u || overlap_second_end_i64.is_some();
        let skip_overlapping = |cur_start: i32| -> bool {
            (overlap_u && cur_start >= mate_alignment_start)
                || overlap_second_end_i64.is_some_and(|max_s| {
                    cur_start >= mate_alignment_start && i64::from(cur_start) <= max_s
                })
        };

        let mut buf_s = std::mem::take(&mut self.scratch_s);
        let mut buf_ss = std::mem::take(&mut self.scratch_ss);
        let mut scratch_splice_key = std::mem::take(&mut self.scratch_splice_key);

        while ci < cigar_elements.len() {
            if has_overlap_gate && skip_overlapping(start) {
                break;
            }
            let Some(&(cigar_element_len, raw_operator)) = cigar_elements.get(ci) else {
                break;
            };
            let operator = if (ci == 0 || ci == cigar_elements.len() - 1)
                && raw_operator == CigarKind::Insertion
            {
                CigarKind::SoftClip
            } else {
                raw_operator
            };

            match operator {
                CigarKind::Skip => {
                    // Reused key buffer + borrowed lookup: zero alloc for
                    // repeat junctions (the common case); two inserts only for
                    // a new junction. Set/count contents identical to the old
                    // clone+entry shape.
                    scratch_splice_key.clear();
                    let _ = write!(
                        scratch_splice_key,
                        "{}-{}",
                        i64::from(start) - 1,
                        i64::from(start) + i64::from(cigar_element_len) - 1
                    );
                    if let Some([cnt]) = self.splice_count.get_mut(scratch_splice_key.as_str()) {
                        *cnt += 1;
                    } else {
                        self.splice.insert(scratch_splice_key.clone());
                        self.splice_count.insert(scratch_splice_key.clone(), [1]);
                    }
                    start =
                        start.saturating_add(i32::try_from(cigar_element_len).unwrap_or(i32::MAX));
                    ci += 1;
                    offset = 0;
                }
                CigarKind::SoftClip => {
                    self.process_soft_clip(
                        record,
                        query_seq,
                        query_qual,
                        &cigar_elements,
                        direction,
                        read_len_match_ins_i32,
                        total_length_including_soft_clipped,
                        adjusted_nm,
                        ci,
                        &mut start,
                        &mut read_pos_inc,
                        &mut read_pos_exc,
                        &mut offset,
                        start_pos,
                    );
                    ci += 1;
                }
                CigarKind::HardClip => {
                    ci += 1;
                    offset = 0;
                }
                CigarKind::Insertion => {
                    offset = 0;
                    ci = self.process_insertion(
                        record,
                        query_seq,
                        query_qual,
                        &cigar_elements,
                        ci,
                        direction,
                        read_len_match_ins_i32,
                        adjusted_nm,
                        &mut start,
                        &mut read_pos_inc,
                        &mut read_pos_exc,
                        &mut offset,
                        start_pos,
                    );
                    ci += 1;
                }
                CigarKind::Deletion => {
                    offset = 0;
                    ci = self.process_deletion(
                        record,
                        query_seq,
                        query_qual,
                        &cigar_elements,
                        ci,
                        direction,
                        read_len_match_ins_i32,
                        adjusted_nm,
                        &mut start,
                        &mut read_pos_inc,
                        &mut read_pos_exc,
                        &mut offset,
                    );
                    ci += 1;
                }
                _ => {
                    // Deal with matching part (Match / SequenceMatch / SequenceMismatch)
                    let mut nmoff = 0;
                    let mut moffset = 0;
                    let mut should_break = false;

                    let vext = self.config.extension_nuc.max(0) as usize;
                    let next_indel_kind = if self.config.realign == 1
                        && vext > 0
                        && !cigar_elements
                            .get(ci + 2)
                            .is_some_and(|&(_, kind)| kind == CigarKind::HardClip)
                    {
                        match cigar_elements.get(ci + 1).map(|&(_, kind)| kind) {
                            Some(CigarKind::Deletion) => Some(CigarKind::Deletion),
                            Some(CigarKind::Insertion) => Some(CigarKind::Insertion),
                            _ => None,
                        }
                    } else {
                        None
                    };
                    let indel_prox_start = if next_indel_kind.is_some() {
                        cigar_element_len.saturating_sub(u32::try_from(vext).unwrap_or(u32::MAX))
                    } else {
                        u32::MAX
                    };

                    let mut i = offset;
                    while i < cigar_element_len {
                        let Some(&ch1) = query_seq.get(read_pos_inc) else {
                            break;
                        };
                        let is_trim_at_opt_t_bases = self.config.trim_bases_after != 0 && {
                            let dist = if direction {
                                i64::from(total_length_including_soft_clipped)
                                    - (read_pos_inc as i64)
                            } else {
                                read_pos_inc as i64
                            };
                            dist > i64::from(self.config.trim_bases_after)
                        };

                        if ch1 == b'N' || ch1 == b'n' {
                            if self.config.include_n {
                                inc_cnt(&mut self.ref_coverage, start, 1);
                            }
                            start = start.saturating_add(1);
                            read_pos_inc += 1;
                            read_pos_exc += 1;
                            i += 1;
                            continue;
                        }

                        let mut q = qual_phred(&query_qual, read_pos_inc);
                        let mut qbases = 1;
                        let mut qibases = 0;
                        buf_ss.clear();
                        let mut start_with_deletion = false;

                        // T048: MNV Growth Look-Ahead
                        // H1: single fused reference lookup (`ref_mismatch_not_n`)
                        // instead of two identical `get_base` calls per base.
                        while start < i32::MAX
                            && (start.saturating_add(1) as u32) >= self.region.start.0
                            && (start.saturating_add(1) as u32) <= self.region.end.0
                            && i + 1 < cigar_element_len
                            && q >= self.config.base_qual
                            && self.ref_mismatch_not_n(
                                start as u32,
                                query_seq.get(read_pos_inc).copied(),
                            )
                        {
                            let next_read_idx = read_pos_inc + 1;
                            let Some(&nuc) = query_seq.get(next_read_idx) else {
                                break;
                            };
                            if (qual_phred(&query_qual, next_read_idx))
                                < self.config.base_qual + 5.0
                            {
                                break;
                            }
                            if nuc == b'N' || nuc == b'n' {
                                break;
                            }
                            if let Some(next_ref_base) =
                                self.get_ref_base(start.saturating_add(1) as u32)
                            {
                                if next_ref_base == b'N' || next_ref_base == b'n' {
                                    break;
                                }
                                if nuc != next_ref_base {
                                    buf_ss.push(nuc as char);
                                    q += qual_phred(&query_qual, next_read_idx);
                                    qbases += 1;
                                    read_pos_inc += 1;
                                    read_pos_exc += 1;
                                    i += 1;
                                    start = start.saturating_add(1);
                                    nmoff += 1;
                                } else {
                                    // CigarParser.java:638 sub-segment lookahead span
                                    let mut subseq_len = 0;
                                    for ssi in 1..=vext {
                                        if i + 1 + u32::try_from(ssi).unwrap_or(u32::MAX)
                                            >= cigar_element_len
                                        {
                                            break;
                                        }
                                        let lookahead_read_idx = read_pos_inc + 1 + ssi;
                                        if let Some(&q_base) = query_seq.get(lookahead_read_idx) {
                                            if let Some(lookahead_ref_base) =
                                                u32::try_from(i64::from(start) + 1 + (ssi as i64))
                                                    .ok()
                                                    .and_then(|p| self.get_ref_base(p))
                                            {
                                                if q_base != lookahead_ref_base {
                                                    subseq_len = ssi + 1;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if subseq_len == 0 {
                                        break;
                                    }
                                    let subseq_read_idx = read_pos_inc + subseq_len;
                                    if subseq_read_idx >= query_seq.len()
                                        || (qual_phred(&query_qual, subseq_read_idx))
                                            < self.config.base_qual + 5.0
                                    {
                                        break;
                                    }
                                    for ssi in 1..=subseq_len {
                                        let idx = read_pos_inc + ssi;
                                        if let Some(&q_base) = query_seq.get(idx) {
                                            push_ascii_base(&mut buf_ss, q_base);
                                        }
                                        q += qual_phred(&query_qual, idx);
                                        qbases += 1;
                                    }
                                    read_pos_inc += subseq_len;
                                    read_pos_exc += subseq_len;
                                    i += u32::try_from(subseq_len).unwrap_or(u32::MAX);
                                    start = start.saturating_add(
                                        i32::try_from(subseq_len).unwrap_or(i32::MAX),
                                    );
                                }
                            } else {
                                break;
                            }
                        }

                        let mut ddlen = 0;
                        let mut grafted = false;
                        if i >= indel_prox_start
                            && next_indel_kind == Some(CigarKind::Deletion)
                            && self.is_closer_than_vext_and_good_base(
                                query_seq,
                                query_qual,
                                &cigar_elements,
                                ci,
                                i,
                                buf_ss.len(),
                                CigarKind::Deletion,
                                start,
                                read_pos_inc,
                            )
                        {
                            grafted = true;
                            buf_s.clear();
                            push_ascii_base(&mut buf_s, ch1);
                            if !buf_ss.is_empty() {
                                buf_s.push('&');
                                buf_s.push_str(&buf_ss);
                            }
                            while i + 1 < cigar_element_len {
                                if let (Some(&qb), Some(&qq)) = (
                                    query_seq.get(read_pos_inc + 1),
                                    query_qual.get(read_pos_inc + 1),
                                ) {
                                    push_ascii_base(&mut buf_s, qb);
                                    q += qq.saturating_sub(33) as f64;
                                    qbases += 1;
                                }
                                i += 1;
                                read_pos_inc += 1;
                                read_pos_exc += 1;
                                start = start.saturating_add(1);
                            }
                            let next_del_len = cigar_elements.get(ci + 1).map_or(0, |&(l, _)| l);
                            let s_clean = buf_s.replacen('&', "", 1);
                            buf_s.clear();
                            let _ = write!(buf_s, "-{}&{}", next_del_len, s_clean);
                            start_with_deletion = true;
                            ddlen = i32::try_from(next_del_len).unwrap_or(i32::MAX);
                            ci += 1;

                            if self.is_two_insertions_ahead(&cigar_elements, ci) {
                                let next_next_len =
                                    cigar_elements.get(ci + 1).map_or(0, |&(l, _)| l as usize);
                                if let Some(sub) = query_seq
                                    .get(read_pos_inc + 1..read_pos_inc + 1 + next_next_len)
                                {
                                    buf_s.push('^');
                                    push_ascii_slice(&mut buf_s, sub);
                                }
                                for qi in 1..=next_next_len {
                                    if let Some(&qq) = query_qual.get(read_pos_inc + 1 + qi) {
                                        q += qq.saturating_sub(33) as f64;
                                        qibases += 1;
                                    }
                                }
                                read_pos_inc += next_next_len;
                                read_pos_exc += next_next_len;
                                ci += 1;
                            }

                            if self.is_next_after_num_matched(&cigar_elements, ci, 1) {
                                let next_next_next_len =
                                    cigar_elements.get(ci + 1).map_or(0, |&(l, _)| l as usize);
                                let offset_read_pos = read_pos_inc + 1;
                                let (toffset, offset_nm) = self.find_offset(
                                    start.saturating_add(ddlen).saturating_add(1),
                                    offset_read_pos,
                                    next_next_next_len,
                                    query_seq,
                                    query_qual,
                                );
                                if toffset != 0 {
                                    moffset = i32::try_from(toffset).unwrap_or(i32::MAX);
                                    nmoff += offset_nm;
                                    buf_s.push('&');
                                    if let Some(sub) =
                                        query_seq.get(offset_read_pos..offset_read_pos + toffset)
                                    {
                                        push_ascii_slice(&mut buf_s, sub);
                                    }
                                    if let Some(q_sub) =
                                        query_qual.get(offset_read_pos..offset_read_pos + toffset)
                                    {
                                        for &b in q_sub {
                                            q += b.saturating_sub(33) as f64;
                                            qibases += 1;
                                        }
                                    }
                                    for osi in 0..moffset {
                                        *self
                                            .ref_coverage
                                            .entry(
                                                start
                                                    .saturating_add(ddlen)
                                                    .saturating_add(1)
                                                    .saturating_add(osi),
                                            )
                                            .or_default() += 1;
                                    }
                                }
                            }
                        } else if i >= indel_prox_start
                            && next_indel_kind == Some(CigarKind::Insertion)
                            && self.is_closer_than_vext_and_good_base(
                                query_seq,
                                query_qual,
                                &cigar_elements,
                                ci,
                                i,
                                buf_ss.len(),
                                CigarKind::Insertion,
                                start,
                                read_pos_inc,
                            )
                        {
                            grafted = true;
                            buf_s.clear();
                            push_ascii_base(&mut buf_s, ch1);
                            if !buf_ss.is_empty() {
                                buf_s.push('&');
                                buf_s.push_str(&buf_ss);
                            }
                            while i + 1 < cigar_element_len {
                                if let (Some(&qb), Some(&qq)) = (
                                    query_seq.get(read_pos_inc + 1),
                                    query_qual.get(read_pos_inc + 1),
                                ) {
                                    push_ascii_base(&mut buf_s, qb);
                                    q += qq.saturating_sub(33) as f64;
                                    qbases += 1;
                                }
                                i += 1;
                                read_pos_inc += 1;
                                read_pos_exc += 1;
                                start = start.saturating_add(1);
                            }
                            let s_clean = buf_s.replacen('&', "", 1);
                            let next_len =
                                cigar_elements.get(ci + 1).map_or(0, |&(l, _)| l as usize);
                            let mut temp_s = s_clean;
                            if let Some(sub) =
                                query_seq.get(read_pos_inc + 1..read_pos_inc + 1 + next_len)
                            {
                                push_ascii_slice(&mut temp_s, sub);
                            }
                            if let (Some(first_part), Some(second_part)) =
                                (temp_s.get(..next_len), temp_s.get(next_len..))
                            {
                                buf_s.clear();
                                let _ = write!(buf_s, "+{}&{}", first_part, second_part);
                            }

                            for qi in 1..=next_len {
                                if let Some(&qq) = query_qual.get(read_pos_inc + 1 + qi) {
                                    q += qq.saturating_sub(33) as f64;
                                    qibases += 1;
                                }
                            }
                            read_pos_inc += next_len;
                            read_pos_exc += next_len;
                            ci += 1;
                            qibases -= 1;
                            qbases += 1;
                        } else if !buf_ss.is_empty() {
                            buf_s.clear();
                            push_ascii_base(&mut buf_s, ch1);
                            buf_s.push('&');
                            buf_s.push_str(&buf_ss);
                        }

                        // Add matching/mismatched Variation (CigarParser.java:582-590)
                        if !is_trim_at_opt_t_bases {
                            let pos = start.saturating_sub(qbases).saturating_add(1);
                            if u32::try_from(pos)
                                .is_ok_and(|p| p >= self.region.start.0 && p <= self.region.end.0)
                            {
                                if !grafted && buf_ss.is_empty() {
                                    let s = SINGLE_BASE_STR
                                        .get(usize::from(ch1))
                                        .copied()
                                        .unwrap_or("N");
                                    if s != "N" {
                                        self.add_variation_for_matching_part(
                                            mapping_quality_i32,
                                            adjusted_nm,
                                            direction,
                                            read_len_match_ins_i32,
                                            nmoff,
                                            s,
                                            false,
                                            q,
                                            qbases,
                                            qibases,
                                            ddlen,
                                            pos,
                                            read_pos_exc as i32,
                                            start,
                                        );
                                    }
                                } else if !buf_s.contains('N') && !buf_s.contains('n') {
                                    self.add_variation_for_matching_part(
                                        mapping_quality_i32,
                                        adjusted_nm,
                                        direction,
                                        read_len_match_ins_i32,
                                        nmoff,
                                        &buf_s,
                                        start_with_deletion,
                                        q,
                                        qbases,
                                        qibases,
                                        ddlen,
                                        pos,
                                        read_pos_exc as i32,
                                        start,
                                    );
                                }
                            }
                        }

                        if start_with_deletion {
                            start = start.saturating_add(ddlen);
                        }

                        if operator != CigarKind::Insertion {
                            start = start.saturating_add(1);
                        }
                        if operator != CigarKind::Deletion {
                            read_pos_inc += 1;
                            read_pos_exc += 1;
                        }

                        if has_overlap_gate && skip_overlapping(start) {
                            should_break = true;
                            break;
                        }

                        i += 1;
                    }

                    if should_break {
                        break;
                    }

                    if moffset != 0 {
                        offset = moffset as u32;
                        read_pos_inc += moffset as usize;
                        start = start.saturating_add(moffset);
                        read_pos_exc += moffset as usize;
                    }
                    ci += 1;
                }
            }

            if start > self.region.end.0 as i32 {
                break;
            }
        }

        // Reclaim per-read buffers into parser scratch for the next read.
        // Output-identical: buffers are refilled from each record before use.
        if let Some((s, q)) = modified_bufs {
            self.scratch_seq = s;
            self.scratch_qual = q;
        }
        self.scratch_elements = cigar_elements;
        self.scratch_s = buf_s;
        self.scratch_ss = buf_ss;
        self.scratch_splice_key = scratch_splice_key;

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn add_variation_for_matching_part(
        &mut self,
        mapping_quality: i32,
        nm: i32,
        dir: bool,
        rlen1: i32,
        nmoff: i32,
        s: &str,
        start_with_deletion: bool,
        mut q: f64,
        qbases: i32,
        qibases: i32,
        ddlen: i32,
        pos: i32,
        read_pos_exc: i32,
        start: i32,
    ) {
        let hv = if s.starts_with('+') {
            increment(&mut self.position_to_insertion_count, pos, s);
            get_variation(&mut self.insertion_variants, pos, s)
        } else {
            if is_mnp(s) {
                increment(&mut self.mnp, pos, s);
            }
            get_variation(&mut self.non_insertion_variants, pos, s)
        };

        hv.inc_dir(dir);
        hv.cnt += 1;

        let tp = if read_pos_exc < rlen1 - read_pos_exc {
            read_pos_exc + 1
        } else {
            rlen1 - read_pos_exc
        };

        q = q / (qbases + qibases) as f64;

        if !hv.pstd && hv.pp != 0 && tp != hv.pp {
            hv.pstd = true;
        }
        if !hv.qstd && hv.pq != 0.0 && q != hv.pq {
            hv.qstd = true;
        }

        hv.pmean += tp as f64;
        hv.qmean += q;
        hv.mqmean += mapping_quality as f64;
        hv.pp = tp;
        hv.pq = q;
        hv.nm += nm.saturating_sub(nmoff) as f64;
        if q >= self.config.base_qual {
            hv.hicnt += 1;
        } else {
            hv.locnt += 1;
        }

        let shift = if s.starts_with('+') && s.contains('&') {
            1
        } else {
            0
        };
        for qi in 1..=(qbases - shift) {
            inc_cnt(
                &mut self.ref_coverage,
                start.saturating_sub(qi).saturating_add(1),
                1,
            );
        }

        if start_with_deletion {
            increment(&mut self.position_to_deletion_count, pos, s);
            let first = i64::from(start) + 1;
            let last = (i64::from(start) + i64::from(ddlen) - 1).min(self.cov_cap);
            for cov_pos in first..=last {
                inc_cnt(&mut self.ref_coverage, cov_pos as i32, 1);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn is_closer_than_vext_and_good_base(
        &self,
        query_seq: &[u8],
        query_qual: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        i: u32,
        ss_len: usize,
        target_kind: CigarKind,
        start: i32,
        read_pos_inc: usize,
    ) -> bool {
        if cigar_elements
            .get(ci + 2)
            .is_some_and(|&(_, kind)| kind == CigarKind::HardClip)
        {
            return false;
        }
        if self.config.realign != 1 {
            return false;
        }
        let Some(&(cigar_element_len, _)) = cigar_elements.get(ci) else {
            return false;
        };
        let Ok(vext) = u32::try_from(self.config.extension_nuc) else {
            return false;
        };

        if cigar_element_len.saturating_sub(i) > vext {
            return false;
        }
        let Some(&(_, next_kind)) = cigar_elements.get(ci + 1) else {
            return false;
        };
        if next_kind != target_kind {
            return false;
        }
        let Ok(pos) = u32::try_from(start) else {
            return false;
        };
        let Some(ref_base) = self.get_ref_base(pos) else {
            return false;
        };
        let Some(&q_base) = query_seq.get(read_pos_inc) else {
            return false;
        };
        let q = query_qual
            .get(read_pos_inc)
            .copied()
            .unwrap_or(33)
            .saturating_sub(33) as f64;
        if q < self.config.base_qual {
            return false;
        }
        ss_len > 0 || q_base != ref_base
    }

    fn is_two_insertions_ahead(&self, cigar_elements: &[(u32, CigarKind)], ci: usize) -> bool {
        cigar_elements
            .get(ci + 1)
            .is_some_and(|&(_, kind)| kind == CigarKind::Insertion)
    }

    fn is_next_after_num_matched(
        &self,
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        number: usize,
    ) -> bool {
        cigar_elements
            .get(ci.saturating_add(number))
            .is_some_and(|&(_, kind)| kind == CigarKind::Match)
    }

    fn is_insertion_or_deletion_with_next_matched(
        &self,
        cigar: &[(u32, CigarKind)],
        ci: usize,
    ) -> bool {
        if self.config.realign != 1 {
            return false;
        }
        let Ok(vext) = u32::try_from(self.config.extension_nuc) else {
            return false;
        };
        let (Some(&(len1, kind1)), Some(&(_, kind2))) = (cigar.get(ci + 1), cigar.get(ci + 2))
        else {
            return false;
        };
        if len1 > vext || kind1 != CigarKind::Match {
            return false;
        }
        if kind2 != CigarKind::Insertion && kind2 != CigarKind::Deletion {
            return false;
        }
        if let Some(&(_, op3)) = cigar.get(ci + 3) {
            if op3 == CigarKind::Insertion || op3 == CigarKind::Deletion {
                return false;
            }
        }
        true
    }

    fn is_next_matched(&self, cigar: &[(u32, CigarKind)], ci: usize) -> bool {
        self.config.realign == 1 && cigar.len() > ci + 1 && cigar[ci + 1].1 == CigarKind::Match
    }

    fn is_next_insertion(&self, cigar: &[(u32, CigarKind)], ci: usize) -> bool {
        self.config.realign == 1 && cigar.len() > ci + 1 && cigar[ci + 1].1 == CigarKind::Insertion
    }

    fn skip_indel_next_to_intron(&self, cigar: &[(u32, CigarKind)], ci: usize) -> bool {
        (cigar.len() > ci + 1 && cigar[ci + 1].1 == CigarKind::Skip)
            || (ci > 0 && cigar[ci - 1].1 == CigarKind::Skip)
    }

    fn append_segments(
        &self,
        query_sequence: &[u8],
        query_quality: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        desc_string: &mut String,
        quality_segment: &mut Vec<u8>,
        m_len: usize,
        indel_len: usize,
        begin: usize,
        is_insertion: bool,
    ) {
        let sub_m = crate::util::subseq(query_sequence, begin, m_len).unwrap_or(&[]);
        desc_string.push('#');
        push_ascii_slice(desc_string, sub_m);

        let qual_m = crate::util::subseq(query_quality, begin, m_len).unwrap_or(&[]);
        quality_segment.extend_from_slice(qual_m);

        desc_string.push('^');
        let is_next_next_ins = cigar_elements[ci + 2].1 == CigarKind::Insertion;
        if is_next_next_ins {
            let sub_ins =
                crate::util::subseq(query_sequence, begin + m_len, indel_len).unwrap_or(&[]);
            push_ascii_slice(desc_string, sub_ins);
        } else {
            let _ = write!(desc_string, "{indel_len}");
        }

        if is_insertion {
            if is_next_next_ins {
                let qual_ins =
                    crate::util::subseq(query_quality, begin + m_len, indel_len).unwrap_or(&[]);
                quality_segment.extend_from_slice(qual_ins);
            } else {
                if let Some(q_char) = crate::util::char_at(query_quality, (begin + m_len) as isize)
                {
                    quality_segment.push(q_char);
                }
            }
        } else {
            if is_next_next_ins {
                let qual_ins =
                    crate::util::subseq(query_quality, begin + m_len, indel_len).unwrap_or(&[]);
                quality_segment.extend_from_slice(qual_ins);
            }
        }
    }

    fn add_variation_for_deletion(
        &mut self,
        mapping_quality: i32,
        nm: i32,
        dir: bool,
        rlen1: i32,
        desc_string: &str,
        quality_segment: &[u8],
        nmoff: i32,
        start: i32,
        read_pos_exc: i32,
        cigar_element_len: u32,
    ) {
        let hv = get_variation(&mut self.non_insertion_variants, start, desc_string);
        increment(&mut self.position_to_deletion_count, start, desc_string);

        hv.inc_dir(dir);
        hv.cnt += 1;

        let tp = if read_pos_exc < rlen1 - read_pos_exc {
            read_pos_exc + 1
        } else {
            rlen1 - read_pos_exc
        };

        let mut tmpq = 0.0;
        if !quality_segment.is_empty() {
            for &b in quality_segment {
                tmpq += b.saturating_sub(33) as f64;
            }
            tmpq /= quality_segment.len() as f64;
        }

        if !hv.pstd && hv.pp != 0 && tp != hv.pp {
            hv.pstd = true;
        }
        if !hv.qstd && hv.pq != 0.0 && tmpq != hv.pq {
            hv.qstd = true;
        }

        hv.pmean += tp as f64;
        hv.qmean += tmpq;
        hv.mqmean += mapping_quality as f64;
        hv.pp = tp;
        hv.pq = tmpq;
        hv.nm += nm.saturating_sub(nmoff) as f64;

        if tmpq >= self.config.base_qual {
            hv.hicnt += 1;
        } else {
            hv.locnt += 1;
        }

        let first = i64::from(start);
        let last = (i64::from(start) + i64::from(cigar_element_len) - 1).min(self.cov_cap);
        for cov_pos in first..=last {
            inc_cnt(&mut self.ref_coverage, cov_pos as i32, 1);
        }
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "per-operator pileup helper borrows split CigarParser fields"
    )]
    fn process_deletion(
        &mut self,
        record: &SamRecord,
        query_seq: &[u8],
        query_qual: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        direction: bool,
        read_len_match_ins: i32,
        adjusted_nm: i32,
        start: &mut i32,
        read_pos_inc: &mut usize,
        read_pos_exc: &mut usize,
        offset: &mut u32,
    ) -> usize {
        let cigar_element_len = cigar_elements[ci].0;
        let cel_i32 = i32::try_from(cigar_element_len).unwrap_or(i32::MAX);
        if self.skip_indel_next_to_intron(cigar_elements, ci) {
            *read_pos_exc += cigar_element_len as usize;
            return ci;
        }

        // F4 scratch: reuse take/put-back buffers instead of allocating per
        // deletion occurrence. Values identical; only allocation sites change.
        let mut desc_string = std::mem::take(&mut self.scratch_f4_desc);
        desc_string.clear();
        let _ = write!(desc_string, "-{cigar_element_len}");
        let mut sequence_to_append = std::mem::take(&mut self.scratch_f4_seg1);
        sequence_to_append.clear();
        let mut quality_of_segment_sb = std::mem::take(&mut self.scratch_f4_seg2);
        quality_of_segment_sb.clear();
        let quality_of_last_segment_before_del = if *read_pos_inc > 0 {
            query_qual.get(*read_pos_inc - 1).copied().unwrap_or(33)
        } else {
            33
        };
        let mut quality_of_segment = std::mem::take(&mut self.scratch_f4_qual);
        quality_of_segment.clear();

        let mut multoffs = 0;
        let mut multoffp = 0;
        let mut nmoff = 0;

        let mut next_ci = ci;

        if self.is_insertion_or_deletion_with_next_matched(cigar_elements, ci) {
            let m_len = cigar_elements[ci + 1].0 as usize;
            let indel_len = cigar_elements[ci + 2].0 as usize;
            let begin = *read_pos_inc;

            quality_of_segment_sb.clear();
            self.append_segments(
                query_seq,
                query_qual,
                cigar_elements,
                ci,
                &mut desc_string,
                &mut quality_of_segment_sb,
                m_len,
                indel_len,
                begin,
                false,
            );
            quality_of_segment.extend_from_slice(&quality_of_segment_sb);

            let is_next_next_del = cigar_elements[ci + 2].1 == CigarKind::Deletion;
            multoffs += m_len + if is_next_next_del { indel_len } else { 0 };
            multoffp += m_len + if !is_next_next_del { indel_len } else { 0 };

            if self.is_next_after_num_matched(cigar_elements, ci, 3) {
                let mut vsn = 0;
                let tn = *read_pos_inc + multoffp;
                let ts = start
                    .saturating_add(i32::try_from(multoffs).unwrap_or(i32::MAX))
                    .saturating_add(cel_i32);
                let vext = self.config.extension_nuc.max(0) as usize;
                let ci3_len = cigar_elements.get(ci + 3).map_or(0, |&(l, _)| l as usize);

                let mut current_offset = 0;
                for vi in 0..ci3_len {
                    if vsn > vext {
                        break;
                    }
                    let tn_idx = tn + vi;
                    if tn_idx >= query_seq.len() {
                        break;
                    }
                    let q_char = query_seq[tn_idx];
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    if (query_qual
                        .get(tn_idx)
                        .copied()
                        .unwrap_or(33)
                        .saturating_sub(33) as f64)
                        < self.config.base_qual
                    {
                        break;
                    }
                    let ts_pos = ts.saturating_add(i32::try_from(vi).unwrap_or(i32::MAX));
                    if let Some(ref_ch) = self.get_ref_base(ts_pos as u32) {
                        if ref_ch == b'N' || ref_ch == b'n' {
                            break;
                        }
                        if q_char != ref_ch {
                            current_offset = vi + 1;
                            nmoff += 1;
                            vsn = 0;
                        } else {
                            vsn += 1;
                        }
                    }
                }
                if current_offset > 0 {
                    *offset = current_offset as u32;
                    if let Some(sub) = query_seq.get(tn..tn + current_offset) {
                        push_ascii_slice(&mut sequence_to_append, sub);
                    }
                    if let Some(sub) = query_qual.get(tn..tn + current_offset) {
                        quality_of_segment.extend_from_slice(sub);
                    }
                }
            }
            next_ci += 2;
        } else if self.is_next_insertion(cigar_elements, ci) {
            let ins_len = cigar_elements[ci + 1].0 as usize;
            if let Some(sub) = query_seq.get(*read_pos_inc..*read_pos_inc + ins_len) {
                desc_string.push('^');
                push_ascii_slice(&mut desc_string, sub);
            }
            if let Some(sub) = query_qual.get(*read_pos_inc..*read_pos_inc + ins_len) {
                quality_of_segment.extend_from_slice(sub);
            }
            multoffp += ins_len;

            if self.is_next_after_num_matched(cigar_elements, ci, 2) {
                let m_len = cigar_elements[ci + 2].0 as usize;
                let mut vsn = 0;
                let tn = *read_pos_inc + multoffp;
                let ts = start.saturating_add(cel_i32);
                let vext = self.config.extension_nuc.max(0) as usize;

                let mut current_offset = 0;
                for vi in 0..m_len {
                    if vsn > vext {
                        break;
                    }
                    let tn_idx = tn + vi;
                    if tn_idx >= query_seq.len() {
                        break;
                    }
                    let q_char = query_seq[tn_idx];
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    if (query_qual
                        .get(tn_idx)
                        .copied()
                        .unwrap_or(33)
                        .saturating_sub(33) as f64)
                        < self.config.base_qual
                    {
                        break;
                    }
                    let ts_pos = ts.saturating_add(i32::try_from(vi).unwrap_or(i32::MAX));
                    if let Some(ref_ch) = self.get_ref_base(ts_pos as u32) {
                        if ref_ch == b'N' || ref_ch == b'n' {
                            break;
                        }
                        if q_char != ref_ch {
                            current_offset = vi + 1;
                            nmoff += 1;
                            vsn = 0;
                        } else {
                            vsn += 1;
                        }
                    }
                }
                if current_offset > 0 {
                    *offset = current_offset as u32;
                    if let Some(sub) = query_seq.get(tn..tn + current_offset) {
                        push_ascii_slice(&mut sequence_to_append, sub);
                    }
                    if let Some(sub) = query_qual.get(tn..tn + current_offset) {
                        quality_of_segment.extend_from_slice(sub);
                    }
                }
            }
            next_ci += 1;
        } else {
            if self.is_next_matched(cigar_elements, ci) {
                let m_len = cigar_elements.get(ci + 1).map_or(0, |&(l, _)| l as usize);
                let mut vsn = 0;
                let vext = self.config.extension_nuc.max(0) as usize;

                let mut current_offset = 0;
                for vi in 0..m_len {
                    if vsn > vext {
                        break;
                    }
                    let idx = *read_pos_inc + vi;
                    if idx >= query_seq.len() {
                        break;
                    }
                    let q_char = query_seq[idx];
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    if (query_qual
                        .get(idx)
                        .copied()
                        .unwrap_or(33)
                        .saturating_sub(33) as f64)
                        < self.config.base_qual
                    {
                        break;
                    }
                    let ts_pos = start
                        .saturating_add(cel_i32)
                        .saturating_add(i32::try_from(vi).unwrap_or(i32::MAX));
                    if let Some(ref_ch) = self.get_ref_base(ts_pos as u32) {
                        if ref_ch == b'N' || ref_ch == b'n' {
                            break;
                        }
                        if q_char != ref_ch {
                            current_offset = vi + 1;
                            nmoff += 1;
                            vsn = 0;
                        } else {
                            vsn += 1;
                        }
                    }
                }
                if current_offset > 0 {
                    *offset = current_offset as u32;
                    if let Some(sub) = query_seq.get(*read_pos_inc..*read_pos_inc + current_offset)
                    {
                        push_ascii_slice(&mut sequence_to_append, sub);
                    }
                    if let Some(sub) = query_qual.get(*read_pos_inc..*read_pos_inc + current_offset)
                    {
                        quality_of_segment.extend_from_slice(sub);
                    }
                }
            }
        }

        if *offset > 0 {
            desc_string.push('&');
            desc_string.push_str(&sequence_to_append);
        }

        let next_q_idx = *read_pos_inc + *offset as usize;
        if next_q_idx >= query_qual.len() {
            quality_of_segment.push(quality_of_last_segment_before_del);
        } else {
            let quality_of_segment_with_offset = query_qual.get(next_q_idx).copied().unwrap_or(33);
            let best_qual = if quality_of_last_segment_before_del > quality_of_segment_with_offset {
                quality_of_last_segment_before_del
            } else {
                quality_of_segment_with_offset
            };
            quality_of_segment.push(best_qual);
        }

        if *start >= self.region.start.0 as i32 && *start <= self.region.end.0 as i32 {
            self.add_variation_for_deletion(
                record.mapping_quality.unwrap_or(0) as i32,
                adjusted_nm,
                direction,
                read_len_match_ins,
                &desc_string,
                &quality_of_segment,
                nmoff,
                *start,
                *read_pos_exc as i32,
                cigar_element_len,
            );
        }

        *start = start
            .saturating_add(cel_i32)
            .saturating_add(i32::try_from(*offset).unwrap_or(i32::MAX))
            .saturating_add(i32::try_from(multoffs).unwrap_or(i32::MAX));
        *read_pos_inc += *offset as usize + multoffp;
        *read_pos_exc += *offset as usize + multoffp;

        // F4 put-back: buffers refilled on next call.
        self.scratch_f4_desc = desc_string;
        self.scratch_f4_seg1 = sequence_to_append;
        self.scratch_f4_qual = quality_of_segment;
        self.scratch_f4_seg2 = quality_of_segment_sb;

        next_ci
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "per-operator pileup helper borrows split CigarParser fields"
    )]
    fn process_insertion(
        &mut self,
        record: &SamRecord,
        query_seq: &[u8],
        query_qual: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        direction: bool,
        read_len_match_ins: i32,
        adjusted_nm: i32,
        start: &mut i32,
        read_pos_inc: &mut usize,
        read_pos_exc: &mut usize,
        offset: &mut u32,
        position: i32,
    ) -> usize {
        let cigar_element_len = cigar_elements[ci].0;
        if self.skip_indel_next_to_intron(cigar_elements, ci) {
            *read_pos_inc += cigar_element_len as usize;
            return ci;
        }

        // F4 scratch: reuse take/put-back buffers instead of allocating per
        // insertion occurrence. Values identical; only allocation sites change.
        let mut desc_string = std::mem::take(&mut self.scratch_f4_desc);
        desc_string.clear();
        if let Some(sub) = query_seq.get(*read_pos_inc..*read_pos_inc + cigar_element_len as usize)
        {
            push_ascii_slice(&mut desc_string, sub);
        }

        let mut quality_string = std::mem::take(&mut self.scratch_f4_qual);
        quality_string.clear();
        if let Some(sub) = query_qual.get(*read_pos_inc..*read_pos_inc + cigar_element_len as usize)
        {
            quality_string.extend_from_slice(sub);
        }

        let mut ss = std::mem::take(&mut self.scratch_f4_seg1);
        ss.clear();
        let mut quality_segment_sb = std::mem::take(&mut self.scratch_f4_seg2);
        quality_segment_sb.clear();
        let mut multoffs = 0;
        let mut multoffp = 0;
        let mut nmoff = 0;

        let mut next_ci = ci;

        if self.is_insertion_or_deletion_with_next_matched(cigar_elements, ci) {
            let m_len = cigar_elements[ci + 1].0 as usize;
            let indel_len = cigar_elements[ci + 2].0 as usize;
            let begin = *read_pos_inc + cigar_element_len as usize;

            quality_segment_sb.clear();
            self.append_segments(
                query_seq,
                query_qual,
                cigar_elements,
                ci,
                &mut desc_string,
                &mut quality_segment_sb,
                m_len,
                indel_len,
                begin,
                true,
            );
            quality_string.extend_from_slice(&quality_segment_sb);

            let is_next_next_del = cigar_elements[ci + 2].1 == CigarKind::Deletion;
            multoffs += m_len + if is_next_next_del { indel_len } else { 0 };
            multoffp += m_len + if !is_next_next_del { indel_len } else { 0 };

            let ci6 = cigar_elements.get(ci + 3).map_or(0, |&(l, _)| l as usize);
            if ci6 != 0
                && cigar_elements
                    .get(ci + 3)
                    .is_some_and(|&(_, k)| k == CigarKind::Match)
            {
                let offset_read_pos = *read_pos_inc + cigar_element_len as usize + multoffp;
                let multoffs_i32 = i32::try_from(multoffs).unwrap_or(i32::MAX);
                let (found_offset, _) = self.find_offset(
                    start.saturating_add(multoffs_i32),
                    offset_read_pos,
                    ci6,
                    query_seq,
                    query_qual,
                );
                *offset = found_offset as u32;
                if found_offset > 0 {
                    if let Some(sub) =
                        query_seq.get(offset_read_pos..offset_read_pos + found_offset)
                    {
                        push_ascii_slice(&mut ss, sub);
                    }
                    if let Some(sub) =
                        query_qual.get(offset_read_pos..offset_read_pos + found_offset)
                    {
                        quality_string.extend_from_slice(sub);
                    }
                    for osi in 0..(found_offset as i32) {
                        inc_cnt(
                            &mut self.ref_coverage,
                            start.saturating_add(multoffs_i32).saturating_add(osi),
                            1,
                        );
                    }
                }
            }
            next_ci += 2;
        } else {
            if self.is_next_matched(cigar_elements, ci) {
                let mut vsn = 0;
                let m_len = cigar_elements.get(ci + 1).map_or(0, |&(l, _)| l as usize);
                let vext = self.config.extension_nuc.max(0) as usize;

                let mut current_offset = 0;
                for vi in 0..m_len {
                    if vsn > vext {
                        break;
                    }
                    let idx = *read_pos_inc + cigar_element_len as usize + vi;
                    let Some(&q_char) = query_seq.get(idx) else {
                        break;
                    };
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    if (query_qual
                        .get(idx)
                        .copied()
                        .unwrap_or(33)
                        .saturating_sub(33) as f64)
                        < self.config.base_qual
                    {
                        break;
                    }
                    let ts_pos = start.saturating_add(i32::try_from(vi).unwrap_or(i32::MAX));
                    if let Some(ref_ch) = self.get_ref_base(ts_pos as u32) {
                        if q_char != ref_ch {
                            current_offset = vi + 1;
                            nmoff += 1;
                            vsn = 0;
                        } else {
                            vsn += 1;
                        }
                    }
                }
                if current_offset > 0 {
                    *offset = current_offset as u32;
                    if let Some(sub) = query_seq.get(
                        *read_pos_inc + cigar_element_len as usize
                            ..*read_pos_inc + cigar_element_len as usize + current_offset,
                    ) {
                        push_ascii_slice(&mut ss, sub);
                    }
                    if let Some(sub) = query_qual.get(
                        *read_pos_inc + cigar_element_len as usize
                            ..*read_pos_inc + cigar_element_len as usize + current_offset,
                    ) {
                        quality_string.extend_from_slice(sub);
                    }
                    for osi in 0..current_offset {
                        inc_cnt(
                            &mut self.ref_coverage,
                            start.saturating_add(i32::try_from(osi).unwrap_or(i32::MAX)),
                            1,
                        );
                    }
                }
            }
        }

        if *offset > 0 {
            desc_string.push('&');
            desc_string.push_str(&ss);
        }

        if *start - 1 >= self.region.start.0 as i32
            && *start - 1 <= self.region.end.0 as i32
            && !desc_string.contains('N')
            && !desc_string.contains('n')
        {
            let mut insertion_position = *start - 1;
            if desc_string
                .bytes()
                .all(|b| matches!(b, b'A' | b'T' | b'G' | b'C'))
                && !desc_string.is_empty()
            {
                let (adjusted_insertion_position, split_idx) = crate::realign::adj_ins_pos_slice(
                    insertion_position,
                    desc_string.as_bytes(),
                    self.reference,
                    &self.region.chr,
                );

                if adjusted_insertion_position != insertion_position
                    && *read_pos_inc as i32 - 1 - (*start - 1 - adjusted_insertion_position) > 0
                {
                    insertion_position = adjusted_insertion_position;
                    if split_idx > 0 {
                        if let (Some(tail), Some(head)) =
                            (desc_string.get(split_idx..), desc_string.get(..split_idx))
                        {
                            ss.clear();
                            ss.push_str(tail);
                            ss.push_str(head);
                            std::mem::swap(&mut desc_string, &mut ss);
                        }
                    }
                }
            }

            ss.clear();
            ss.push('+');
            ss.push_str(&desc_string);
            increment(
                &mut self.position_to_insertion_count,
                insertion_position,
                ss.as_str(),
            );

            let hv = get_variation(
                &mut self.insertion_variants,
                insertion_position,
                ss.as_str(),
            );
            hv.inc_dir(direction);
            hv.cnt += 1;

            let rpe = *read_pos_exc as i32;
            let tp = if rpe < read_len_match_ins - rpe {
                rpe + 1
            } else {
                read_len_match_ins - rpe
            };

            let mut tmpq = 0.0;
            if !quality_string.is_empty() {
                for &b in &quality_string {
                    tmpq += b.saturating_sub(33) as f64;
                }
                tmpq /= quality_string.len() as f64;
            }

            if !hv.pstd && hv.pp != 0 && tp != hv.pp {
                hv.pstd = true;
            }
            if !hv.qstd && hv.pq != 0.0 && tmpq != hv.pq {
                hv.qstd = true;
            }
            hv.pmean += tp as f64;
            hv.qmean += tmpq;
            hv.mqmean += record.mapping_quality.unwrap_or(0) as f64;
            hv.pp = tp;
            hv.pq = tmpq;
            if tmpq >= self.config.base_qual {
                hv.hicnt += 1;
            } else {
                hv.locnt += 1;
            }
            hv.nm += adjusted_nm.saturating_sub(nmoff) as f64;
            let (hv_pstd, hv_qstd) = (hv.pstd, hv.qstd);

            let index_of_insertion_in_query_sequence =
                *read_pos_inc as i32 - 1 - (*start - 1 - insertion_position);
            if insertion_position > position {
                if let Some(sub_base) = query_seq
                    .get(index_of_insertion_in_query_sequence as usize)
                    .copied()
                {
                    if let Some(ref_ch) = self.get_ref_base(insertion_position as u32) {
                        if sub_base == ref_ch {
                            let base_qual_val = query_qual
                                .get(index_of_insertion_in_query_sequence as usize)
                                .copied()
                                .unwrap_or(33)
                                .saturating_sub(33)
                                as f64;
                            if let Some(tv) = get_variation_maybe(
                                &mut self.non_insertion_variants,
                                insertion_position,
                                Some(sub_base),
                            ) {
                                tv.sub_cnt(
                                    direction,
                                    tp,
                                    base_qual_val,
                                    record.mapping_quality.unwrap_or(0) as i32,
                                    adjusted_nm.saturating_sub(nmoff) as f64,
                                    self.config.base_qual,
                                );
                            }
                        }
                    }
                }
            }

            if ci == 1
                && cigar_elements.len() > 0
                && (cigar_elements[0].1 == CigarKind::SoftClip
                    || cigar_elements[0].1 == CigarKind::HardClip)
            {
                if let Some(ref_base_char) = self.get_ref_base(insertion_position as u32) {
                    let mut buf = [0u8; 4];
                    let ref_base_str = allele_key(ref_base_char, &mut buf);
                    let ttref = get_variation(
                        &mut self.non_insertion_variants,
                        insertion_position,
                        ref_base_str,
                    );
                    ttref.inc_dir(direction);
                    ttref.cnt += 1;
                    ttref.pstd = hv_pstd;
                    ttref.qstd = hv_qstd;
                    ttref.pmean += tp as f64;
                    ttref.qmean += tmpq;
                    ttref.mqmean += record.mapping_quality.unwrap_or(0) as f64;
                    ttref.pp = tp;
                    ttref.pq = tmpq;
                    ttref.nm += adjusted_nm.saturating_sub(nmoff) as f64;
                    inc_cnt(&mut self.ref_coverage, insertion_position, 1);
                }
            }
        }

        *read_pos_inc += cigar_element_len as usize + *offset as usize + multoffp;
        *read_pos_exc += cigar_element_len as usize + *offset as usize + multoffp;
        *start = start
            .saturating_add(i32::try_from(*offset).unwrap_or(i32::MAX))
            .saturating_add(i32::try_from(multoffs).unwrap_or(i32::MAX));

        // F4 put-back: buffers refilled on next call.
        self.scratch_f4_desc = desc_string;
        self.scratch_f4_qual = quality_string;
        self.scratch_f4_seg1 = ss;
        self.scratch_f4_seg2 = quality_segment_sb;

        next_ci
    }

    fn is_read_chimeric_with_sa(
        &self,
        record: &SamRecord,
        position: i32,
        sa_tag_string: &str,
        dir: bool,
        is_5_side: bool,
    ) -> bool {
        let mut parts = sa_tag_string.split(',');
        let (Some(sa_chromosome), Some(sa_pos_str), Some(sa_direction_string), Some(sa_cigar)) =
            (parts.next(), parts.next(), parts.next(), parts.next())
        else {
            return false;
        };
        let sa_position: i32 = match sa_pos_str.parse() {
            Ok(p) => p,
            Err(_) => return false,
        };
        let sa_direction_is_forward = sa_direction_string == "+";

        let is_chimeric = if is_5_side {
            crate::cigar::modify::leading_2digit_s_len(sa_cigar).is_some()
        } else {
            crate::cigar::modify::trailing_2digit_s_len(sa_cigar).is_some()
        };

        let sa_chr_matches = if let Some(ref_name) = &record.reference_name {
            sa_chromosome == ref_name
        } else {
            false
        };

        let is_chimeric_flag = ((dir && sa_direction_is_forward)
            || (!dir && !sa_direction_is_forward))
            && sa_chr_matches
            && ((i64::from(sa_position) - i64::from(position)).abs()
                < 2 * i64::from(self.max_read_length))
            && is_chimeric;

        if self.config.verbose && is_chimeric_flag {
            let record_name = String::from_utf8_lossy(&record.name);
            let ref_name = record.reference_name.as_deref().unwrap_or("*");
            crate::verbose_eprintln!(
                "{} {} {} {:?} is ignored as chimeric with SA: {},{},{}",
                record_name,
                ref_name,
                position,
                record.cigar,
                sa_position,
                sa_direction_string,
                sa_cigar
            );
        }

        is_chimeric_flag
    }

    fn sclip5_high_quality_processing(
        &mut self,
        query_sequence: &[u8],
        mapping_quality: i32,
        query_quality: &[u8],
        number_of_mismatches: i32,
        dir: bool,
        sum_of_read_qualities: i32,
        number_of_high_quality_bases: i32,
        number_of_low_quality_bases: i32,
        start: i32,
        cigar_element_len: u32,
    ) {
        if number_of_high_quality_bases >= 1
            && number_of_high_quality_bases > number_of_low_quality_bases
            && start >= self.region.start.0 as i32
            && start <= self.region.end.0 as i32
        {
            let sclip = self
                .soft_clips_5_end
                .entry(start)
                .or_insert_with(Sclip::new);
            let mut si = cigar_element_len as i32 - 1;
            while cigar_element_len as i32 - si <= number_of_high_quality_bases {
                if si < 0 || si >= query_sequence.len() as i32 {
                    break;
                }
                let ch = query_sequence[si as usize];
                let idx = cigar_element_len as i32 - 1 - si;

                let cnts = sclip.nt.entry(idx).or_default();
                *cnts.entry(ch).or_default() += 1;

                let seq_variation = get_variation_from_seq(sclip, idx, ch);
                seq_variation.add_cnt(
                    dir,
                    si - (cigar_element_len as i32 - number_of_high_quality_bases),
                    qual_phred(&query_quality, si as usize),
                    mapping_quality,
                    number_of_mismatches as f64,
                    self.config.base_qual,
                );
                si -= 1;
            }
            sclip.variation.add_cnt(
                dir,
                cigar_element_len as i32,
                sum_of_read_qualities as f64 / number_of_high_quality_bases as f64,
                mapping_quality,
                number_of_mismatches as f64,
                self.config.base_qual,
            );
        }
    }

    fn sclip3_high_quality_processing(
        &mut self,
        query_sequence: &[u8],
        mapping_quality: i32,
        query_quality: &[u8],
        number_of_mismatches: i32,
        dir: bool,
        sum_of_read_qualities: i32,
        number_of_high_quality_bases: i32,
        number_of_low_quality_bases: i32,
        start: i32,
        cigar_element_len: u32,
        read_pos_inc: usize,
    ) {
        if number_of_high_quality_bases >= 1
            && number_of_high_quality_bases > number_of_low_quality_bases
            && start >= self.region.start.0 as i32
            && start <= self.region.end.0 as i32
        {
            let sclip = self
                .soft_clips_3_end
                .entry(start)
                .or_insert_with(Sclip::new);
            for si in 0..number_of_high_quality_bases as usize {
                let idx = read_pos_inc + si;
                if idx >= query_sequence.len() {
                    break;
                }
                let ch = query_sequence[idx];

                let cnts = sclip.nt.entry(si as i32).or_default();
                *cnts.entry(ch).or_default() += 1;

                let seq_variation = get_variation_from_seq(sclip, si as i32, ch);
                seq_variation.add_cnt(
                    dir,
                    number_of_high_quality_bases - si as i32,
                    qual_phred(&query_quality, idx),
                    mapping_quality,
                    number_of_mismatches as f64,
                    self.config.base_qual,
                );
            }
            sclip.variation.add_cnt(
                dir,
                cigar_element_len as i32,
                sum_of_read_qualities as f64 / number_of_high_quality_bases as f64,
                mapping_quality,
                number_of_mismatches as f64,
                self.config.base_qual,
            );
        }
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "per-operator pileup helper borrows split CigarParser fields"
    )]
    fn process_soft_clip(
        &mut self,
        record: &SamRecord,
        query_seq: &[u8],
        query_qual: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        direction: bool,
        _read_len_match_ins: i32,
        total_length_including_soft_clipped: i32,
        adjusted_nm: i32,
        ci: usize,
        start: &mut i32,
        read_pos_inc: &mut usize,
        read_pos_exc: &mut usize,
        offset: &mut u32,
        position: i32,
    ) {
        let cigar_element_len = cigar_elements[ci].0;

        if ci == 0 {
            if !self.config.chimeric {
                // ports CigarParser.java:1110-1128: the SEED_1 branch runs whenever
                // there is no SA tag, including clips >= 20 without one.
                if cigar_element_len >= 20 && record.sa.is_some() {
                    if let Some(sa_tag_string) = &record.sa {
                        if self.is_read_chimeric_with_sa(
                            record,
                            position,
                            sa_tag_string,
                            direction,
                            true,
                        ) {
                            *read_pos_inc += cigar_element_len as usize;
                            *offset = 0;
                            *start = position;
                            return;
                        }
                    }
                } else if cigar_element_len as usize >= 17 {
                    if let Some(reverse_complemented_seed) =
                        reverse_complement_seed17(record, 0, cigar_element_len as usize, true)
                    {
                        if let Some([pos0]) = self
                            .reference
                            .seed
                            .get(reverse_complemented_seed.as_slice())
                        {
                            if (i64::from(*start) - i64::from(*pos0)).abs()
                                < 2 * i64::from(self.max_read_length)
                            {
                                *read_pos_inc += cigar_element_len as usize;
                                *offset = 0;
                                *start = position;
                                if self.config.verbose {
                                    let sequence = get_reverse_complemented_sequence(
                                        record,
                                        0,
                                        cigar_element_len as usize,
                                    );
                                    crate::verbose_eprintln!(
                                        "{} at 5' is a chimeric at {} by SEED 17",
                                        String::from_utf8_lossy(&sequence),
                                        start
                                    );
                                }

                                return;
                            }
                        }
                    }
                }
            }

            let mut element_len = cigar_element_len;
            while element_len > 0 && *start - 1 > 0 {
                // §9: single fetch — the loop condition only needs is_some.
                let Some(ref_base) = self.get_ref_base((*start - 1) as u32) else {
                    break;
                };
                let q_base = match query_seq.get(element_len as usize - 1) {
                    Some(&b) => b,
                    None => break,
                };
                let base_qual = query_qual
                    .get(element_len as usize - 1)
                    .copied()
                    .unwrap_or(33)
                    .saturating_sub(33);
                if q_base == ref_base && base_qual > 10 {
                    let mut buf = [0u8; 4];
                    let ref_base_str = allele_key(ref_base, &mut buf);
                    let variation =
                        get_variation(&mut self.non_insertion_variants, *start - 1, ref_base_str);
                    variation.add_cnt(
                        direction,
                        element_len as i32,
                        base_qual as f64,
                        record.mapping_quality.unwrap_or(0) as i32,
                        adjusted_nm as f64,
                        self.config.base_qual,
                    );
                    inc_cnt(&mut self.ref_coverage, *start - 1, 1);
                    *start -= 1;
                    element_len -= 1;
                } else {
                    break;
                }
            }

            if element_len > 0 {
                let mut sum_of_read_qualities = 0;
                let mut number_of_high_quality_bases = 0;
                let mut number_of_low_quality_bases = 0;
                let mut si = (i32::try_from(element_len).unwrap_or(i32::MAX) - 1)
                    .min(i32::try_from(query_seq.len()).unwrap_or(i32::MAX) - 1);
                while si >= 0 {
                    let Some(&q_char) = query_seq.get(si as usize) else {
                        break;
                    };
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    let base_qual = query_qual
                        .get(si as usize)
                        .copied()
                        .unwrap_or(33)
                        .saturating_sub(33);
                    if base_qual <= 12 {
                        number_of_low_quality_bases += 1;
                    }
                    if number_of_low_quality_bases > 1 {
                        break;
                    }
                    sum_of_read_qualities += base_qual as i32;
                    number_of_high_quality_bases += 1;
                    si -= 1;
                }
                self.sclip5_high_quality_processing(
                    query_seq,
                    record.mapping_quality.unwrap_or(0) as i32,
                    query_qual,
                    adjusted_nm,
                    direction,
                    sum_of_read_qualities,
                    number_of_high_quality_bases,
                    number_of_low_quality_bases,
                    *start,
                    element_len,
                );
            }
        } else if ci == cigar_elements.len() - 1 {
            if !self.config.chimeric {
                // ports CigarParser.java:1204-1216: the SEED_1 branch runs whenever
                // there is no SA tag, including clips >= 20 without one.
                if cigar_element_len >= 20 && record.sa.is_some() {
                    if let Some(sa_tag_string) = &record.sa {
                        if self.is_read_chimeric_with_sa(
                            record,
                            position,
                            sa_tag_string,
                            direction,
                            false,
                        ) {
                            *read_pos_inc += cigar_element_len as usize;
                            *offset = 0;
                            *start = position;
                            return;
                        }
                    }
                } else if cigar_element_len as usize >= 17 {
                    if let Some(reverse_complemented_seed) = reverse_complement_seed17(
                        record,
                        -(cigar_element_len as isize),
                        cigar_element_len as usize,
                        false,
                    ) {
                        if let Some([pos0]) = self
                            .reference
                            .seed
                            .get(reverse_complemented_seed.as_slice())
                        {
                            if (i64::from(*start) - i64::from(*pos0)).abs()
                                < 2 * i64::from(self.max_read_length)
                            {
                                *read_pos_inc += cigar_element_len as usize;
                                *offset = 0;
                                *start = position;
                                if self.config.verbose {
                                    let sequence = get_reverse_complemented_sequence(
                                        record,
                                        -(cigar_element_len as isize),
                                        cigar_element_len as usize,
                                    );
                                    crate::verbose_eprintln!(
                                        "{} at 3' is a chimeric at {} by SEED 17",
                                        String::from_utf8_lossy(&sequence),
                                        start
                                    );
                                }

                                return;
                            }
                        }
                    }
                }
            }

            let mut element_len = i32::try_from(cigar_element_len).unwrap_or(i32::MAX);
            while *read_pos_inc < query_seq.len() {
                // §9: single fetch — the loop condition only needs is_some.
                let Some(ref_base) = self.get_ref_base(*start as u32) else {
                    break;
                };
                let q_base = query_seq[*read_pos_inc];
                let base_qual = query_qual
                    .get(*read_pos_inc)
                    .copied()
                    .unwrap_or(33)
                    .saturating_sub(33);
                if q_base == ref_base && base_qual > 10 {
                    let mut buf = [0u8; 4];
                    let ref_base_str = allele_key(ref_base, &mut buf);
                    let variation =
                        get_variation(&mut self.non_insertion_variants, *start, ref_base_str);
                    variation.add_cnt(
                        direction,
                        total_length_including_soft_clipped
                            .saturating_sub(i32::try_from(*read_pos_exc).unwrap_or(i32::MAX)),
                        base_qual as f64,
                        record.mapping_quality.unwrap_or(0) as i32,
                        adjusted_nm as f64,
                        self.config.base_qual,
                    );
                    inc_cnt(&mut self.ref_coverage, *start, 1);
                    *read_pos_inc += 1;
                    *start = start.saturating_add(1);
                    element_len -= 1;
                    *read_pos_exc += 1;
                } else {
                    break;
                }
            }

            if query_seq.len() > *read_pos_inc {
                let mut sum_of_read_qualities = 0;
                let mut number_of_high_quality_bases = 0;
                let mut number_of_low_quality_bases = 0;
                let mut si = 0i32;
                while si < element_len {
                    let idx = *read_pos_inc + si as usize;
                    if idx >= query_seq.len() {
                        break;
                    }
                    let q_char = query_seq[idx];
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    let base_qual = query_qual
                        .get(idx)
                        .copied()
                        .unwrap_or(33)
                        .saturating_sub(33);
                    if base_qual <= 12 {
                        number_of_low_quality_bases += 1;
                    }
                    if number_of_low_quality_bases > 1 {
                        break;
                    }
                    sum_of_read_qualities += base_qual as i32;
                    number_of_high_quality_bases += 1;
                    si += 1;
                }
                self.sclip3_high_quality_processing(
                    query_seq,
                    record.mapping_quality.unwrap_or(0) as i32,
                    query_qual,
                    adjusted_nm,
                    direction,
                    sum_of_read_qualities,
                    number_of_high_quality_bases,
                    number_of_low_quality_bases,
                    *start,
                    element_len.max(0) as u32,
                    *read_pos_inc,
                );
            }
        }

        *read_pos_inc += cigar_element_len as usize;
        *offset = 0;
        *start = position;
    }

    fn find_offset(
        &self,
        reference_position: i32,
        read_position: usize,
        cigar_length: usize,
        query_sequence: &[u8],
        query_quality: &[u8],
    ) -> (usize, i32) {
        let mut offset = 0;
        let mut tnm = 0;
        let mut vsn = 0;
        let vext = self.config.extension_nuc.max(0) as usize;

        for vi in 0..cigar_length {
            let read_idx = read_position + vi;
            let Some(&ch) = query_sequence.get(read_idx) else {
                break;
            };
            if ch == b'N' || ch == b'n' {
                break;
            }
            let q = qual_phred(&query_quality, read_idx);
            if q < self.config.base_qual {
                break;
            }
            let ref_pos = reference_position.saturating_add(i32::try_from(vi).unwrap_or(i32::MAX));
            if let Some(ref_ch) = self.get_ref_base(ref_pos as u32) {
                if ch != ref_ch {
                    offset = vi + 1;
                    tnm += 1;
                    vsn = 0;
                } else {
                    vsn += 1;
                }
            }
            if vsn > vext {
                break;
            }
        }

        (offset, tnm)
    }

    /// Finalize and return the collected VariationData.
    /// Ports `CigarParser.java:101-150`.
    /// Note (G2): `CigarParser.java:78,162-163` copies `totalReads` and
    /// `duplicateReads` from `InitialData` in `initFromScope(scope)` BEFORE the
    /// `processor.nextRecord()` loop at `:86`, so both fields are 0 here in the
    /// production pipeline and `duprate` always evaluates to `0.0`.
    pub fn into_variation_data(self) -> VariationData {
        let duprate = if self.svflag {
            if (self.discordant_count + 1) as f64
                / (self.total_reads - self.duplicate_reads + 1) as f64
                > 0.5
            {
                0.0
            } else {
                1.0
            }
        } else {
            if self.config.dedup && self.total_reads != 0 {
                format!(
                    "{:.3}",
                    self.duplicate_reads as f64 / self.total_reads as f64
                )
                .parse::<f64>()
                .unwrap_or(0.0)
            } else {
                0.0
            }
        };

        VariationData {
            non_insertion_variants: self.non_insertion_variants,
            insertion_variants: self.insertion_variants,
            position_to_insertion_count: self.position_to_insertion_count,
            position_to_deletion_count: self.position_to_deletion_count,
            sv_structures: self.sv_structures,
            ref_coverage: self.ref_coverage,
            soft_clips_5_end: self.soft_clips_5_end,
            soft_clips_3_end: self.soft_clips_3_end,
            max_read_length: self.max_read_length,
            splice: self.splice,
            mnp: self.mnp,
            splice_count: self.splice_count,
            duprate,
            sv_structures_filtered: false,
            // G3 (mate-window ownership): the pileup region is the owning
            // tile's active range (whole region in whole-region runs). Merge
            // temps built over mate windows are discarded, never emitted
            // through, so stamping unconditionally is safe.
            owning_range: Some((self.region.start, self.region.end)),
            merge_planted: BTreeSet::new(),
        }
    }

    fn prepare_sv_structures_for_analysis(
        &mut self,
        record: &SamRecord,
        query_qual: &[u8],
        number_of_mismatches: i32,
        read_direction: bool,
        mate_direction: bool,
        start: i32,
        total_length_including_soft_clipped: i32,
        cigar_elements: &[(u32, CigarKind)],
    ) {
        let mate_start = record.mate_alignment_start.map(|p| p.0 as i32).unwrap_or(0);
        let mend = mate_start.saturating_add(total_length_including_soft_clipped);

        let mut end = start;
        for &(len, op) in cigar_elements {
            if matches!(op, CigarKind::Match | CigarKind::Deletion | CigarKind::Skip) {
                end = end.saturating_add(i32::try_from(len).unwrap_or(i32::MAX));
            }
        }

        let mut soft5 = 0;
        if let Some(&(len, op)) = cigar_elements.first() {
            if matches!(op, CigarKind::SoftClip | CigarKind::HardClip) {
                let len = len as usize;
                if len > 0 && len <= query_qual.len() {
                    let qual_val = qual_phred(&query_qual, len.saturating_sub(1));
                    if qual_val > self.config.base_qual {
                        soft5 = start;
                    }
                }
            }
        }

        let mut soft3 = 0;
        if let Some(&(len, op)) = cigar_elements.last() {
            if matches!(op, CigarKind::SoftClip | CigarKind::HardClip) {
                let len = len as usize;
                if len > 0 && query_qual.len() >= len {
                    let qual_val = qual_phred(&query_qual, query_qual.len().saturating_sub(len));
                    if qual_val > self.config.base_qual {
                        soft3 = end;
                    }
                }
            }
        }

        let read_dir_num: i32 = if read_direction { -1 } else { 1 };
        let read_dir_num_64 = i64::from(read_dir_num);
        let mate_dir_num: i32 = if mate_direction { 1 } else { -1 };
        let min_d: i64 = 75;
        let start_64 = i64::from(start);

        if crate::bam::mate_reference_is_equal(record) {
            let mut mlen = i64::from(record.template_length);

            if let Some(ref mc) = record.mc {
                if has_two_soft_clips(mc) {
                    return;
                }
            }

            if let Some(mate_mq) = record.mq {
                if mate_mq < 15 {
                    return;
                }
            }

            let mq = record.mapping_quality.unwrap_or(0);

            let min_svc_dist = 1.5;
            let limit = (min_svc_dist * self.max_read_length as f64) as i64;

            if read_dir_num * mate_dir_num == -1
                && (mlen * read_dir_num_64) > 0
                && query_qual.len() > 15
            {
                // Deletion candidate
                mlen = if mate_start > start {
                    i64::from(mend) - start_64
                } else {
                    i64::from(end) - i64::from(mate_start)
                };
                let sv_min_len = i64::from(self.config.sv_read_len)
                    + (self.config.sv_min_af as f64 * self.config.sv_flanking as f64) as i64;
                if mlen.abs() > sv_min_len {
                    let mlen_i32 = mlen.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
                    if read_dir_num == 1 {
                        if self.sv_structures.svfdel.is_empty()
                            || start_64 - i64::from(self.sv_structures.svdelfend) > limit
                        {
                            let mut sclip = Sclip::new();
                            sclip.variation.cnt = 0;
                            self.sv_structures.svfdel.push(sclip);
                        }
                        if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                mlen_i32,
                                soft3,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                        }
                        self.sv_structures.svdelfend = end;
                    } else {
                        if self.sv_structures.svrdel.is_empty()
                            || start_64 - i64::from(self.sv_structures.svdelrend) > limit
                        {
                            let mut sclip = Sclip::new();
                            sclip.variation.cnt = 0;
                            self.sv_structures.svrdel.push(sclip);
                        }
                        if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                mlen_i32,
                                soft5,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                        }
                        self.sv_structures.svdelrend = end;
                    }

                    if !self.sv_structures.svfdel.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svdelfend)).abs() <= limit
                    {
                        if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdel.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svdelrend)).abs() <= limit
                    {
                        if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfdup.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svdupfend)).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdup.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svduprend)).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfinv5.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svinvfend5)).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfinv5.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrinv5.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svinvrend5)).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrinv5.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfinv3.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svinvfend3)).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfinv3.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrinv3.is_empty()
                        && (start_64 - i64::from(self.sv_structures.svinvrend3)).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrinv3.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                }
            } else if read_dir_num * mate_dir_num == -1
                && read_dir_num_64 * mlen < 0
                && query_qual.len() > 15
            {
                // Duplication
                let mlen_i32 = mlen.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
                if read_dir_num == 1 {
                    if self.sv_structures.svfdup.is_empty()
                        || start_64 - i64::from(self.sv_structures.svdupfend) > limit
                    {
                        let mut sclip = Sclip::new();
                        sclip.variation.cnt = 0;
                        self.sv_structures.svfdup.push(sclip);
                    }
                    if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                        add_sv(
                            sdref,
                            start,
                            end,
                            mate_start,
                            mend,
                            read_dir_num,
                            total_length_including_soft_clipped,
                            mlen_i32,
                            soft3,
                            self.max_read_length as f64 / 2.0,
                            qual_phred(&query_qual, 15),
                            mq as f64,
                            number_of_mismatches as f64,
                            self.config.base_qual,
                        );
                    }
                    self.sv_structures.svdupfend = end;
                } else {
                    if self.sv_structures.svrdup.is_empty()
                        || start_64 - i64::from(self.sv_structures.svduprend) > limit
                    {
                        let mut sclip = Sclip::new();
                        sclip.variation.cnt = 0;
                        self.sv_structures.svrdup.push(sclip);
                    }
                    if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                        add_sv(
                            sdref,
                            start,
                            end,
                            mate_start,
                            mend,
                            read_dir_num,
                            total_length_including_soft_clipped,
                            mlen_i32,
                            soft5,
                            self.max_read_length as f64 / 2.0,
                            qual_phred(&query_qual, 15),
                            mq as f64,
                            number_of_mismatches as f64,
                            self.config.base_qual,
                        );
                    }
                    self.sv_structures.svduprend = end;
                }
                if !self.sv_structures.svfdup.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svdupfend)).abs() <= limit
                {
                    if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrdup.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svduprend)).abs() <= limit
                {
                    if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svfdel.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svdelfend)).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrdel.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svdelrend)).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svfinv5.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svinvfend5)).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svfinv5.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrinv5.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svinvrend5)).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svrinv5.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svfinv3.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svinvfend3)).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svfinv3.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrinv3.is_empty()
                    && (start_64 - i64::from(self.sv_structures.svinvrend3)).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svrinv3.last_mut() {
                        sdref.disc += 1;
                    }
                }
            } else if read_dir_num * mate_dir_num == 1 && query_qual.len() > 15 {
                // Inversion
                let mlen_i32 = mlen.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
                let max_rl_64 = i64::from(self.max_read_length);
                if read_dir_num == 1 && mlen != 0 {
                    if mlen < -3 * max_rl_64 {
                        if self.sv_structures.svfinv3.is_empty()
                            || start_64 - i64::from(self.sv_structures.svinvfend3) > limit
                        {
                            let mut sclip = Sclip::new();
                            sclip.variation.cnt = 0;
                            self.sv_structures.svfinv3.push(sclip);
                        }
                        if let Some(sdref) = self.sv_structures.svfinv3.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                mlen_i32,
                                soft3,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                        self.sv_structures.svinvfend3 = end;
                    } else if mlen > 3 * max_rl_64 {
                        if self.sv_structures.svfinv5.is_empty()
                            || start_64 - i64::from(self.sv_structures.svinvfend5) > limit
                        {
                            let mut sclip = Sclip::new();
                            sclip.variation.cnt = 0;
                            self.sv_structures.svfinv5.push(sclip);
                        }
                        if let Some(sdref) = self.sv_structures.svfinv5.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                mlen_i32,
                                soft3,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                        self.sv_structures.svinvfend5 = end;
                    }
                } else if mlen != 0 {
                    if mlen < -3 * max_rl_64 {
                        if self.sv_structures.svrinv3.is_empty()
                            || start_64 - i64::from(self.sv_structures.svinvrend3) > limit
                        {
                            let mut sclip = Sclip::new();
                            sclip.variation.cnt = 0;
                            self.sv_structures.svrinv3.push(sclip);
                        }
                        if let Some(sdref) = self.sv_structures.svrinv3.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                mlen_i32,
                                soft5,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                        self.sv_structures.svinvrend3 = end;
                    } else if mlen > 3 * max_rl_64 {
                        if self.sv_structures.svrinv5.is_empty()
                            || start_64 - i64::from(self.sv_structures.svinvrend5) > limit
                        {
                            let mut sclip = Sclip::new();
                            sclip.variation.cnt = 0;
                            self.sv_structures.svrinv5.push(sclip);
                        }
                        if let Some(sdref) = self.sv_structures.svrinv5.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                mlen_i32,
                                soft5,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                        self.sv_structures.svinvrend5 = end;
                    }
                }

                if mlen != 0 {
                    // ports CigarParser.java:2190-2202: signed comparisons, no abs().
                    if !self.sv_structures.svfdel.is_empty()
                        && start_64 - i64::from(self.sv_structures.svdelfend) <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdel.is_empty()
                        && start_64 - i64::from(self.sv_structures.svdelrend) <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfdup.is_empty()
                        && start_64 - i64::from(self.sv_structures.svdupfend) <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdup.is_empty()
                        && start_64 - i64::from(self.sv_structures.svduprend) <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                }
            }
        } else if query_qual.len() > 15 {
            // Inter-chr translocation (rare: mate mapped elsewhere — the only
            // path that needs the owned name, built here instead of per read).
            let mchr = get_mate_reference_name(record);
            // ports CigarParser.java:2208-2215: the MC-tag and mate-MQ(MQ tag) filters,
            // not the read's own mapping quality.
            let sv_allowed = !(record.mc.as_ref().is_some_and(|mc| has_two_soft_clips(mc))
                || record.mq.is_some_and(|mate_mq| mate_mq < 15));
            let mq = record.mapping_quality.unwrap_or(0);

            let min_svc_dist = 1.5;
            let limit = (min_svc_dist * self.max_read_length as f64) as i64;

            if sv_allowed {
                if read_dir_num == 1 {
                    let last_fus_end = self
                        .sv_structures
                        .svfusfend
                        .get(&mchr)
                        .copied()
                        .unwrap_or(0);
                    if !self.sv_structures.svffus.contains_key(&mchr)
                        || start_64 - i64::from(last_fus_end) > limit
                    {
                        let mut sclip = Sclip::new();
                        sclip.variation.cnt = 0;
                        self.sv_structures
                            .svffus
                            .entry(mchr.clone())
                            .or_default()
                            .push(sclip);
                    }
                    if let Some(sclips) = self.sv_structures.svffus.get_mut(&mchr) {
                        if let Some(sdref) = sclips.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                0,
                                soft3,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                    }
                    self.sv_structures.svfusfend.insert(mchr, end);
                } else {
                    let last_fus_end = self
                        .sv_structures
                        .svfusrend
                        .get(&mchr)
                        .copied()
                        .unwrap_or(0);
                    if !self.sv_structures.svrfus.contains_key(&mchr)
                        || start_64 - i64::from(last_fus_end) > limit
                    {
                        let mut sclip = Sclip::new();
                        sclip.variation.cnt = 0;
                        self.sv_structures
                            .svrfus
                            .entry(mchr.clone())
                            .or_default()
                            .push(sclip);
                    }
                    if let Some(sclips) = self.sv_structures.svrfus.get_mut(&mchr) {
                        if let Some(sdref) = sclips.last_mut() {
                            add_sv(
                                sdref,
                                start,
                                end,
                                mate_start,
                                mend,
                                read_dir_num,
                                total_length_including_soft_clipped,
                                0,
                                soft5,
                                self.max_read_length as f64 / 2.0,
                                qual_phred(&query_qual, 15),
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                    }
                    self.sv_structures.svfusrend.insert(mchr, end);
                }
            }

            let min_sv_pos: i64 = 25;
            if !self.sv_structures.svfdel.is_empty()
                && start_64 - i64::from(self.sv_structures.svdelfend) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrdel.is_empty()
                && start_64 - i64::from(self.sv_structures.svdelrend) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svfdup.is_empty()
                && start_64 - i64::from(self.sv_structures.svdupfend) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrdup.is_empty()
                && start_64 - i64::from(self.sv_structures.svduprend) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svfinv5.is_empty()
                && start_64 - i64::from(self.sv_structures.svinvfend5) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfinv5.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrinv5.is_empty()
                && start_64 - i64::from(self.sv_structures.svinvrend5) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrinv5.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svfinv3.is_empty()
                && start_64 - i64::from(self.sv_structures.svinvfend3) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfinv3.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrinv3.is_empty()
                && start_64 - i64::from(self.sv_structures.svinvrend3) <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrinv3.last_mut() {
                    sdref.disc += 1;
                }
            }
        }
    }
}

fn add_sv(
    sdref: &mut Sclip,
    start_s: i32,
    end_e: i32,
    mate_start_ms: i32,
    mate_end_me: i32,
    dir: i32,
    rlen: i32,
    mlen: i32,
    softp: i32,
    pmean_rp: f64,
    qmean: f64,
    mqmean: f64,
    nm: f64,
    goodq: f64,
) {
    sdref.variation.cnt += 1;
    sdref.variation.inc_dir(dir != 1);

    if qmean >= goodq {
        sdref.variation.hicnt += 1;
    } else {
        sdref.variation.locnt += 1;
    }

    if sdref.start == 0 || sdref.start >= start_s {
        sdref.start = start_s;
    }
    if sdref.end == 0 || sdref.end <= end_e {
        sdref.end = end_e;
    }
    sdref.mates.push(Mate {
        mate_start: mate_start_ms,
        mate_end: mate_end_me,
        mate_len: mlen,
        start: start_s,
        end: end_e,
        pmean: pmean_rp,
        qmean,
        mqmean,
        nm,
    });

    if sdref.mstart == 0 || sdref.mstart >= mate_start_ms {
        sdref.mstart = mate_start_ms;
    }
    if sdref.mend == 0 || sdref.mend <= mate_end_me {
        sdref.mend = mate_start_ms.saturating_add(rlen);
    }
    if softp != 0 {
        if dir == 1 {
            if (i64::from(softp) - i64::from(sdref.end)).abs() < 10 {
                let soft_count = sdref.soft.entry(softp).or_insert(0);
                *soft_count += 1;
            }
        } else {
            if (i64::from(softp) - i64::from(sdref.start)).abs() < 10 {
                let soft_count = sdref.soft.entry(softp).or_insert(0);
                *soft_count += 1;
            }
        }
    }
}

fn is_paired_and_same_chromosome(record: &SamRecord) -> bool {
    (record.flags & 0x1) != 0 && crate::bam::mate_reference_is_equal(record)
}

fn get_reference_length(cigar: &[(u32, CigarKind)]) -> u32 {
    cigar
        .iter()
        .filter(|&&(_, kind)| {
            matches!(
                kind,
                CigarKind::Match
                    | CigarKind::Deletion
                    | CigarKind::Skip
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
            )
        })
        .fold(0u32, |acc, &(len, _)| acc.saturating_add(len))
}

/// Helper function to calculate the length of insertion and deletion CIGAR elements.
/// Ports `CigarParser.java:1636-1644`.
pub fn get_insertion_deletion_length(cigar: &[(u32, CigarKind)]) -> u32 {
    cigar
        .iter()
        .filter(|(_, kind)| matches!(kind, CigarKind::Insertion | CigarKind::Deletion))
        .fold(0u32, |acc, &(len, _)| acc.saturating_add(len))
}

/// Sum of query-consuming CIGAR operators (`M | I | S | = | X`) in `u64` (no `u32` wrap).
/// Used by the BAM record query-consumption consistency gate (SAMv1 §1.4.6).
pub fn query_consumed_len(cigar: &[(u32, CigarKind)]) -> u64 {
    cigar
        .iter()
        .filter(|(_, kind)| {
            matches!(
                kind,
                CigarKind::Match
                    | CigarKind::Insertion
                    | CigarKind::SoftClip
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
            )
        })
        .map(|&(len, _)| u64::from(len))
        .sum()
}

/// Helper to calculate the read length including matched, soft-clipped and insertions segments.
/// Ports `CigarParser.java:1605-1614`.
pub fn get_soft_clipped_length(cigar: &[(u32, CigarKind)]) -> u32 {
    cigar
        .iter()
        .filter(|(_, kind)| {
            matches!(
                kind,
                CigarKind::Match | CigarKind::Insertion | CigarKind::SoftClip
            )
        })
        .fold(0u32, |acc, &(len, _)| acc.saturating_add(len))
}

/// Parse the `-a` amplicon parameter "distance:fraction" into (distance, fraction).
/// Falls back to (10, 0.95) on parse failure, matching Java parseCigarWithAmpCase.
fn parse_amp_param(param: &str) -> (i32, f64) {
    let mut parts = param.split(':');
    let distance = parts.next().and_then(|s| s.parse::<i32>().ok());
    let fraction = parts.next().and_then(|s| s.parse::<f64>().ok());
    match (distance, fraction) {
        (Some(d), Some(f)) => (d, f),
        _ => (10, 0.95),
    }
}

/// Sum of matched and deleted segment lengths (Java getAlignedLength: M + D).
fn get_aligned_length_amp(cigar: &[(u32, CigarKind)]) -> i32 {
    cigar
        .iter()
        .filter(|(_, kind)| matches!(kind, CigarKind::Match | CigarKind::Deletion))
        .fold(0i32, |acc, &(len, _)| {
            acc.saturating_add(i32::try_from(len).unwrap_or(i32::MAX))
        })
}

/// Ports `CigarParser.parseCigarWithAmpCase`. Returns true if the read should be
/// filtered out (does not belong to the amplicon).
fn parse_cigar_with_amp_case(
    cigar_elements: &[(u32, CigarKind)],
    record: &crate::bam::SamRecord,
    is_mate_ref_name_equal: bool,
    region: &crate::regions::Region,
    distance_to_amplicon: i32,
    overlap_fraction: f64,
) -> bool {
    if cigar_elements.is_empty() {
        return false;
    }
    let read_length = i64::from(get_aligned_length_amp(cigar_elements));
    let mut segstart = record.alignment_start.map_or(0i64, |p| i64::from(p.0));
    let mut segend = segstart + read_length - 1;
    let region_start = i64::from(region.start.0);
    let region_end = i64::from(region.end.0);
    let distance_to_amplicon = i64::from(distance_to_amplicon);
    let tlen = i64::from(record.template_length);

    let first_op = cigar_elements.first().map(|&(_, op)| op);
    let last_op = cigar_elements.last().map(|&(_, op)| op);

    // Java `!(ratio > overlapFraction)` is true for NaN too; partial_cmp != Greater matches that.
    let not_gt =
        |ratio: f64| ratio.partial_cmp(&overlap_fraction) != Some(std::cmp::Ordering::Greater);

    if first_op == Some(CigarKind::SoftClip) {
        let ts1 = segstart.max(region_start);
        let te1 = segend.min(region_end);
        not_gt((ts1 - te1).abs() as f64 / (segend - segstart) as f64)
    } else if last_op == Some(CigarKind::SoftClip) {
        let ts1 = segstart.max(region_start);
        let te1 = segend.min(region_end);
        not_gt((te1 - ts1).abs() as f64 / (segend - segstart) as f64)
    } else {
        if is_mate_ref_name_equal && tlen != 0 {
            if tlen > 0 {
                segend = segstart + tlen - 1;
            } else {
                segstart = record.mate_alignment_start.map_or(0i64, |p| i64::from(p.0));
                segend = segstart - tlen - 1;
            }
        }
        let ts1 = segstart.max(region_start);
        let te1 = segend.min(region_end);
        (segstart - region_start).abs() > distance_to_amplicon
            || (segend - region_end).abs() > distance_to_amplicon
            || not_gt((ts1 - te1).abs() as f64 / (segend - segstart) as f64)
    }
}

/// Helper to calculate the read length including matched and insertions segments.
/// Ports `CigarParser.java:1621-1629`.
pub fn get_match_insertion_length(cigar: &[(u32, CigarKind)]) -> u32 {
    cigar
        .iter()
        .filter(|(_, kind)| matches!(kind, CigarKind::Match | CigarKind::Insertion))
        .fold(0u32, |acc, &(len, _)| acc.saturating_add(len))
}

/// Fused single-lookup form of the old MNV guard conjunction
/// `is_has_and_not_equals(ref, chr, pos, qb) && is_not_equals_n(get_base(..))`:
/// one `get_base` (linear region scan + map lookups) per evaluation instead
/// of two identical ones. Truth table identical by construction (unit-tested).
/// `None` query base (read exhausted) yields `false` instead of panicking on a
/// bare index (AGENTS.md §9); the old code only survived that case by never
/// reaching it. Ports `CigarParser.java` MNV-growth mismatch gate.
#[cfg(test)]
#[inline]
fn ref_mismatch_not_n(reference: &Reference, chr: &str, pos: u32, query_base: Option<u8>) -> bool {
    match (reference.get_base(chr, pos), query_base) {
        (Some(ref_base), Some(qb)) => ref_base != qb && ref_base != b'N' && ref_base != b'n',
        _ => false,
    }
}

fn is_mnp(s: &str) -> bool {
    let mut parts = s.split('&');
    if let (Some(first), Some(second), None) = (parts.next(), parts.next(), parts.next()) {
        if first.len() == 1 {
            let is_first_base = first
                .bytes()
                .all(|b| matches!(b, b'A' | b'T' | b'G' | b'C' | b'a' | b't' | b'g' | b'c'));
            let is_second_base = second
                .bytes()
                .all(|b| matches!(b, b'A' | b'T' | b'G' | b'C' | b'a' | b't' | b'g' | b'c'));
            return is_first_base && is_second_base && !second.is_empty();
        }
    }
    false
}

fn increment(hash: &mut BTreeMap<i32, rustc_hash::FxHashMap<String, i32>>, k1: i32, k2: &str) {
    let map = hash.entry(k1).or_default();
    if let Some(cnt) = map.get_mut(k2) {
        *cnt += 1;
    } else {
        map.insert(k2.to_owned(), 1);
    }
}

fn inc_cnt(hash: &mut rustc_hash::FxHashMap<i32, i32>, start: i32, count: i32) {
    *hash.entry(start).or_default() += count;
}

fn get_variation_from_seq(sclip: &mut Sclip, idx: i32, ch: u8) -> &mut Variation {
    let map = sclip.seq.entry(idx).or_default();
    map.entry(ch).or_insert_with(Variation::new)
}

fn reverse_complement_seed17(
    record: &SamRecord,
    start_index: isize,
    length: usize,
    take_rc_head: bool,
) -> Option<[u8; 17]> {
    let start_idx = if start_index < 0 {
        record
            .sequence
            .len()
            .checked_sub(start_index.unsigned_abs())?
    } else {
        usize::try_from(start_index).ok()?
    };
    let end_idx = start_idx.saturating_add(length).min(record.sequence.len());
    let sub = record.sequence.get(start_idx..end_idx)?;
    let raw17 = if take_rc_head {
        sub.get(sub.len().checked_sub(17)?..)?
    } else {
        sub.get(..17)?
    };
    let mut seed = [0u8; 17];
    seed.copy_from_slice(raw17);
    crate::util::reverse_complement(&mut seed);
    Some(seed)
}

fn get_reverse_complemented_sequence(
    record: &SamRecord,
    mut start_index: isize,
    length: usize,
) -> Vec<u8> {
    let read_len = record.sequence.len() as isize;
    if start_index < 0 {
        start_index = read_len + start_index;
    }
    let start_idx = start_index as usize;
    if start_idx < record.sequence.len() {
        let end_idx = (start_idx + length).min(record.sequence.len());
        let mut sub = record.sequence[start_idx..end_idx].to_vec();
        crate::util::reverse_complement(&mut sub);
        sub
    } else {
        Vec::new()
    }
}

/// Main pileup / CIGAR parsing engine interface.
/// Ports `CigarParser.java:75-150`.
pub fn pileup(
    records: &[SamRecord],
    reference: &Reference,
    region: &Region,
    config: &Config,
    svflag: bool,
    total_reads: i32,
    duplicate_reads: i32,
    seed_max_read_length: i32,
) -> Result<VariationData, VardictError> {
    let mut parser = CigarParser::new(
        reference,
        region,
        config,
        svflag,
        total_reads,
        duplicate_reads,
        seed_max_read_length,
    );
    for record in records {
        parser.parse_cigar(record)?;
    }
    Ok(parser.into_variation_data())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Pos1;

    #[test]
    fn allele_key_preserves_case() {
        // `allele_key` must return identical keys to `(base as char).to_string()`
        // (Java char cast): lowercase `n` stays `"n"`, never folds to `"N"`.
        let mut buf = [0u8; 4];
        for (base, want) in [
            (b'A', "A"),
            (b'C', "C"),
            (b'G', "G"),
            (b'T', "T"),
            (b'N', "N"),
            (b'n', "n"),
            (b'a', "a"),
            (b'c', "c"),
            (b'g', "g"),
            (b't', "t"),
        ] {
            assert_eq!(allele_key(base, &mut buf), want);
        }
        assert_eq!(allele_key(b'X', &mut buf), "X");
    }

    #[test]
    fn ref_mismatch_not_n_fuses_double_lookup() {
        // H1 direction: the MNV look-ahead evaluated `is_has_and_not_equals`
        // (one `get_base`: linear region scan + map lookups) AND
        // `is_not_equals_n(get_base(..))` (a second, identical `get_base`) per
        // base. `ref_mismatch_not_n` must equal that conjunction with one lookup.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 20,
            });
        let mut seq = vec![b'A'; 20];
        seq[10] = b'N'; // pos 11
        seq[11] = b'G'; // pos 12
        reference.region_sequences.push(seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), seq);
        // (query_base, pos, expected conjunction)
        for (qb, pos, want) in [
            (Some(b'C'), 10u32, true),  // A != C, ref not N
            (Some(b'A'), 10u32, false), // ref == query
            (Some(b'C'), 11u32, false), // ref is N
            (Some(b'N'), 11u32, false), // ref == query (both N)
            (Some(b'C'), 99u32, false), // pos uncovered -> None
            (Some(b'n'), 10u32, true),  // A != n, ref not N
            (Some(b'A'), 12u32, true),  // G != A, ref not N
            (Some(b'G'), 12u32, false), // ref == query
            (None, 10u32, false),       // read exhausted -> no growth, no panic
        ] {
            assert_eq!(
                ref_mismatch_not_n(&reference, "chr1", pos, qb),
                want,
                "qb={:?} pos={}",
                qb.map(|b| b as char),
                pos
            );
        }
    }

    #[test]
    fn mnv_growth_records_identical_hicnt_locnt_qmean_nm() {
        // Parity anchor for the MNV look-ahead hoist: a two-base MNV must keep
        // producing identical hicnt/locnt/qmean/nm into tovars/realign inputs.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        let seq = vec![b'A'; 100];
        reference.region_sequences.push(seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), seq);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        // Isolate the pileup MNV loop from CigarModifier rewriting (realign=1
        // would soft-clip/shift the mismatching 5' bases before pileup).
        config.realign = 0;
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        let rec = SamRecord {
            name: b"r_mnv".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(17)),
            mapping_quality: Some(60),
            cigar: vec![(8, CigarKind::Match)],
            sequence: b"CGAAAAAA".to_vec(),
            quality_scores: vec![70; 8],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(2),
            sa: None,
            mc: None,
            mq: None,
        };
        parser.parse_cigar(&rec).unwrap();
        let hv = parser
            .non_insertion_variants
            .get(&10)
            .and_then(|entry| entry.map.get("C&G"))
            .expect("MNV variation C&G at pos 10 must exist");
        assert_eq!(hv.cnt, 1);
        assert_eq!(hv.hicnt, 1);
        assert_eq!(hv.locnt, 0);
        assert_eq!(hv.qmean, 37.0);
        assert_eq!(hv.nm, 1.0);
        assert_eq!(hv.pmean, 2.0);
        assert_eq!(hv.mqmean, 60.0);
    }

    #[test]
    fn mismatch_gate() {
        let reference = Reference::new();
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(10),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        config.max_mismatches = 2;

        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        let record_pass = SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(20)),
            mapping_quality: Some(60),
            cigar: vec![(10, CigarKind::Match)],
            sequence: b"AAAAAAAAAA".to_vec(),
            quality_scores: vec![30; 10],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(1),
            sa: None,
            mc: None,
            mq: None,
        };

        let record_indel_pass = SamRecord {
            name: b"read3".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(20)),
            mapping_quality: Some(60),
            cigar: vec![(8, CigarKind::Match), (2, CigarKind::Insertion)],
            sequence: b"AAAAAAAAAA".to_vec(),
            quality_scores: vec![30; 10],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(4),
            sa: None,
            mc: None,
            mq: None,
        };

        assert!(parser.parse_cigar(&record_pass).is_ok());
        assert_eq!(get_insertion_deletion_length(&record_pass.cigar), 0);
        assert_eq!(get_insertion_deletion_length(&record_indel_pass.cigar), 2);
    }

    #[test]
    fn filter_double_softclips() {
        let reference = Reference::new();
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(10),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let config = Config::default();

        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(30)),
            mapping_quality: Some(60),
            cigar: vec![
                (12, CigarKind::SoftClip),
                (10, CigarKind::Match),
                (15, CigarKind::SoftClip),
            ],
            sequence: b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
            quality_scores: vec![30; 37],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: None,
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();
        assert_eq!(parser.max_read_length, 0);
    }

    #[test]
    fn match_loop_test() {
        // Setup Reference with a simple chromosome and sequence
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 20,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(20),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut config = Config::default();
        config.realign = 0;
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        // A read record with matches and mismatches
        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(5)),
            mapping_quality: Some(60),
            cigar: vec![(5, CigarKind::Match)],
            // Ref:   A C G T A
            // Read:  A C T T A   (mismatch G -> T at position 3)
            sequence: b"ACTTA".to_vec(),
            quality_scores: vec![60; 5],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(1),
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        // At position 3 (1-based), we should have a Variation for "T" in non_insertion_variants
        let var_t = parser
            .non_insertion_variants
            .get(&3)
            .unwrap()
            .map
            .get("T")
            .unwrap();
        assert_eq!(var_t.cnt, 1);

        // At position 1 (1-based), we should have a Variation for "A" matching the reference
        let var_a = parser
            .non_insertion_variants
            .get(&1)
            .unwrap()
            .map
            .get("A")
            .unwrap();
        assert_eq!(var_a.cnt, 1);

        // Ref coverage should be incremented at each matched/mismatched position
        assert_eq!(parser.ref_coverage.get(&1), Some(&1));
        assert_eq!(parser.ref_coverage.get(&2), Some(&1));
        assert_eq!(parser.ref_coverage.get(&3), Some(&1));
        assert_eq!(parser.ref_coverage.get(&4), Some(&1));
        assert_eq!(parser.ref_coverage.get(&5), Some(&1));
    }

    #[test]
    fn del_test() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 20,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), b"ACGTACGTACGT".to_vec());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(20),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut config = Config::default();
        config.realign = 0;
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(7)),
            mapping_quality: Some(60),
            cigar: vec![
                (2, CigarKind::Match),
                (3, CigarKind::Deletion),
                (2, CigarKind::Match),
            ],
            sequence: b"ACGT".to_vec(),
            quality_scores: vec![60; 4],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        // At position 3 (1-based), we should have a Deletion Variation for "-3"
        let var_del = parser
            .non_insertion_variants
            .get(&3)
            .unwrap()
            .map
            .get("-3")
            .unwrap();
        assert_eq!(var_del.cnt, 1);
    }

    #[test]
    fn ins_test() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 20,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), b"ACGTACGTACGT".to_vec());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(20),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut config = Config::default();
        config.realign = 0;
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(4)),
            mapping_quality: Some(60),
            cigar: vec![
                (2, CigarKind::Match),
                (2, CigarKind::Insertion),
                (2, CigarKind::Match),
            ],
            sequence: b"ACTTGT".to_vec(),
            quality_scores: vec![60; 6],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        // At position 2 (1-based), we should have an Insertion Variation for "+TT"
        let var_ins = parser
            .insertion_variants
            .get(&2)
            .unwrap()
            .map
            .get("+TT")
            .unwrap();
        assert_eq!(var_ins.cnt, 1);
    }

    #[test]
    fn splice_test() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 20,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), b"ACGTACGTACGT".to_vec());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(20),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut config = Config::default();
        config.realign = 0;
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(9)),
            mapping_quality: Some(60),
            cigar: vec![
                (2, CigarKind::Match),
                (5, CigarKind::Skip),
                (2, CigarKind::Match),
            ],
            sequence: b"ACGT".to_vec(),
            quality_scores: vec![60; 4],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        assert!(parser.splice.contains("2-7"));
        assert_eq!(parser.splice_count.get("2-7"), Some(&[1]));
    }

    #[test]
    fn softclip_test() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 20,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), b"ACGTACGTACGT".to_vec());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(20),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut config = Config::default();
        config.realign = 0;
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(5)),
            mapping_quality: Some(60),
            cigar: vec![(2, CigarKind::SoftClip), (5, CigarKind::Match)],
            sequence: b"ATACGTA".to_vec(),
            quality_scores: vec![60; 7],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        // At position 1 (1-based), we should have a 5' soft clip entry
        assert!(parser.soft_clips_5_end.contains_key(&1));
    }

    fn chimeric_seed_test_record(sa: Option<String>) -> SamRecord {
        // 24bp 3' clip whose reverse-complement seeds uniquely into the reference.
        let clip = *b"CAGTGAATTACAGTTAAAAGATCG";
        let mut seq = b"ACGTA".to_vec();
        seq.extend_from_slice(&clip);
        SamRecord {
            name: b"read1".to_vec(),
            flags: 0,
            reference_name: Some("chrZ".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(5)),
            mapping_quality: Some(60),
            cigar: vec![(5, CigarKind::Match), (24, CigarKind::SoftClip)],
            sequence: seq,
            quality_scores: vec![60; 29],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa,
            mc: None,
            mq: None,
        }
    }

    #[test]
    fn softclip_3p_chimeric_seed_check_runs_without_sa_tag() {
        // ports CigarParser.java:1204-1216: a clip >= SEED_1 without an SA tag must
        // go through the seed-based chimeric check (the SA branch requires an SA tag).
        let mut reference = Reference::new();
        let ref_seq = vec![b'T'; 100];
        reference
            .reference_sequences
            .insert("chrZ".to_string(), ref_seq);
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chrZ".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        let record = chimeric_seed_test_record(None);
        let seq_tail = get_reverse_complemented_sequence(&record, -24, 24);
        let seed = seq_tail[seq_tail.len() - 17..].to_vec();
        assert_eq!(seed.len(), 17);
        reference.seed.insert(seed, vec![30]);

        let region = Region {
            chr: "chrZ".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        parser.parse_cigar(&record).unwrap();

        // The clip is chimeric: no 3' cluster may be registered at the walk position.
        assert!(parser.soft_clips_3_end.is_empty());
    }

    #[test]
    fn softclip_3p_seed_check_skipped_when_sa_tag_present() {
        // ports CigarParser.java:1110-1128: with an SA tag present and len >= 20 only
        // the SA check runs; when it says non-chimeric the clip is registered even if
        // the seed would have matched.
        let mut reference = Reference::new();
        let ref_seq = vec![b'T'; 100];
        reference
            .reference_sequences
            .insert("chrZ".to_string(), ref_seq);
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chrZ".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        // Same unique seed as above; must be ignored because an SA tag exists.
        let record = chimeric_seed_test_record(Some("chrOther,30,+,24M,60,0".to_string()));
        let seq_tail = get_reverse_complemented_sequence(&record, -24, 24);
        let seed = seq_tail[seq_tail.len() - 17..].to_vec();
        reference.seed.insert(seed, vec![30]);

        let region = Region {
            chr: "chrZ".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);

        parser.parse_cigar(&record).unwrap();

        // Non-chimeric per SA: the remaining clip is registered as usual.
        assert!(
            parser.soft_clips_3_end.keys().next().is_some(),
            "expected a 3' soft-clip cluster"
        );
    }

    #[test]
    fn sv_collection_uses_nm_minus_indel_length() {
        // ports CigarParser.java:227-236: SV structures receive NM - (I + D), not raw NM.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 200]);

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        parser.max_read_length = 150;

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0x61, // paired, mate-reverse, read1, forward
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(157)),
            mapping_quality: Some(60),
            cigar: vec![
                (40, CigarKind::Match),
                (3, CigarKind::Insertion),
                (107, CigarKind::Match),
            ],
            sequence: vec![b'A'; 150],
            quality_scores: vec![60; 150],
            mate_reference_name: Some("=".to_string()),
            mate_alignment_start: Some(Pos1(500)),
            template_length: -520, // negative -> duplication branch
            nm: Some(3),           // 3 == the 3I; adjusted NM must be 0
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        let dup = parser.sv_structures.svfdup.last().expect("dup cluster");
        assert_eq!(dup.variation.cnt, 1);
        assert_eq!(dup.variation.nm, 0.0, "nm must be NM minus indel length");
        assert_eq!(dup.mates[0].nm, 0.0);
    }

    #[test]
    fn interchromosomal_sv_ignores_own_mapq_and_mate_mc_tag() {
        // ports CigarParser.java:2208-2215: the fusion branch filters on the MC tag
        // and the mate MQ tag, never on the read's own mapping quality.
        let mut reference = Reference::new();
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 200]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();

        let mk_record = |mc: Option<String>, mq: Option<u8>| SamRecord {
            name: b"read1".to_vec(),
            flags: 0x61, // paired, mate-reverse, read1, forward
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(50)),
            alignment_end: Some(Pos1(73)),
            mapping_quality: Some(13), // above 10 but below 15: must not filter
            cigar: vec![(24, CigarKind::Match)],
            sequence: vec![b'A'; 24],
            quality_scores: vec![60; 24],
            mate_reference_name: Some("chr2".to_string()),
            mate_alignment_start: Some(Pos1(900)),
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc,
            mq,
        };

        // Own MAPQ 13 but no MQ tag: Java collects the pair.
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        parser.max_read_length = 150;
        parser.parse_cigar(&mk_record(None, None)).unwrap();
        let fus = parser
            .sv_structures
            .svffus
            .get("chr2")
            .expect("fusion list");
        assert_eq!(fus.last().unwrap().variation.cnt, 1);

        // Mate MQ tag < 15: Java skips collection.
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        parser.max_read_length = 150;
        parser.parse_cigar(&mk_record(None, Some(10u8))).unwrap();
        assert!(!parser.sv_structures.svffus.contains_key("chr2"));

        // MC tag with clips on both ends: Java skips collection.
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        parser.max_read_length = 150;
        parser
            .parse_cigar(&mk_record(Some("20S50M23S".to_string()), None))
            .unwrap();
        assert!(!parser.sv_structures.svffus.contains_key("chr2"));
    }

    #[test]
    fn del_branch_dup_disc_bump_uses_abs_distance() {
        // ports CigarParser.java:2064: the deletion-tail dup disc bump uses
        // abs(start - svdupfend) <= MIN_D. Signed forms exist only in the
        // inversion-tail block (CigarParser.java:2191-2202).
        let mut reference = Reference::new();
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 300]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(300),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        parser.max_read_length = 100;

        // Dup cluster ending far to the RIGHT of the upcoming read:
        // |start - svdupfend| = 200 > MIN_D -> must NOT bump.
        // Mate-reverse flag (0x20) makes dirnum*mdirnum == -1 with positive
        // template length -> deletion branch. Mate placed far away so the
        // recomputed |mlen| = mateStart+softclips-start = 950 clears Java's
        // INSSIZE+INSSTD gate (:2038) — the disc bumps live inside it.
        // Pre-seeded del cluster guards branch entry: its tail bump uses
        // |start - svdelfend(0)| = 50 <= MINSVCDIST*maxReadLength = 150.
        let mut sc = Sclip::new();
        sc.end = 250;
        parser.sv_structures.svfdup.push(sc);
        parser.sv_structures.svdupfend = 250;
        parser.sv_structures.svfdel.push(Sclip::new());

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0x21, // paired + mate reverse -> deletion branch
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(50)),
            alignment_end: Some(Pos1(177)),
            mapping_quality: Some(60),
            cigar: vec![(128, CigarKind::Match)],
            sequence: vec![b'A'; 128],
            quality_scores: vec![60; 128],
            mate_reference_name: Some("=".to_string()),
            mate_alignment_start: Some(Pos1(1000)),
            template_length: 1000,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        let del = parser.sv_structures.svfdel.first().unwrap();
        assert_eq!(
            del.disc, 1,
            "deletion cluster within MINSVCDIST*maxReadLength must bump"
        );
        let dup = parser.sv_structures.svfdup.first().unwrap();
        assert_eq!(
            dup.disc, 0,
            "abs(start - svdupfend) > MIN_D must not fire on far-left starts"
        );
    }

    #[test]
    fn inversion_tail_disc_bump_uses_signed_distance() {
        // ports CigarParser.java:2191-2202: the INVERSION branch's tail-bump
        // block (reached when mlen != 0 but |mlen| <= 3*maxReadLength, so no
        // inversion cluster is created) uses SIGNED start - tail comparisons
        // with NO abs(). W2.2 erratum restore: the old deleted test exercised
        // this signed site while claiming to cover the abs() deletion-branch
        // site (see del_branch_dup_disc_bump_uses_abs_distance above).
        //
        // svfdup ends at 250, read starts at 50: start - svdupfend = -200,
        // which is <= MIN_D(75) signed -> MUST bump, although |diff| = 200
        // would reject under an abs() form. svrdel ends at -100:
        // start - svdelrend = +150 > 75 -> must NOT bump.
        let mut reference = Reference::new();
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 300]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(300),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        parser.max_read_length = 100;

        let mut dup_sc = Sclip::new();
        dup_sc.end = 250;
        parser.sv_structures.svfdup.push(dup_sc);
        parser.sv_structures.svdupfend = 250;
        let mut rdel_sc = Sclip::new();
        rdel_sc.end = -100;
        parser.sv_structures.svrdel.push(rdel_sc);
        parser.sv_structures.svdelrend = -100;

        let record = SamRecord {
            name: b"read1".to_vec(),
            flags: 0x1, // paired only: read and mate both forward -> dirnum*mdirnum == 1 -> inversion branch
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(50)),
            alignment_end: Some(Pos1(177)),
            mapping_quality: Some(60),
            cigar: vec![(128, CigarKind::Match)],
            sequence: vec![b'A'; 128],
            quality_scores: vec![60; 128],
            mate_reference_name: Some("=".to_string()),
            mate_alignment_start: Some(Pos1(1000)),
            template_length: 100, // != 0 but |mlen| <= 3*max_read_length -> skips cluster creation
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };

        parser.parse_cigar(&record).unwrap();

        let dup = parser.sv_structures.svfdup.first().unwrap();
        assert_eq!(
            dup.disc, 1,
            "signed diff -200 <= MIN_D(75) must bump; an abs() form would reject at 200 > 75"
        );
        let rdel = parser.sv_structures.svrdel.first().unwrap();
        assert_eq!(rdel.disc, 0, "signed diff +150 > MIN_D(75) must not bump");
    }

    #[test]
    fn crispr_filter_uses_1based_start_without_extra_offset() {
        // ports CigarParser.java:288, 292-318, 634-646:
        // 1. `start` is already the 1-based alignment start (`position`), so
        //    `skipSitesOutRegionOfInterest` evaluates
        //    `cutSite - start > filterBp && start + rlen3 - cutSite > filterBp`
        //    without adding an extra +1 to `start`.
        // 2. `rlen3` uses `ALIGNED_LENGTH_MD` (`[MD]` only, excluding `N`/Skip).
        // 3. The filter runs after `maxReadLength` is updated (:309-318).
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "T".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        config.realign = 0;
        config.crispr_cutting_site = Some(50);
        config.crispr_flanking_size = Some(10);

        // Case 1: start = 39, rlen3 = 30.
        // Java: cutSite - start = 50 - 39 = 11 > 10 (true)
        //       start + rlen3 - cutSite = 39 + 30 - 50 = 19 > 10 (true) -> ACCEPTED.
        // Buggy +1 offset gives 50 - 40 = 10 > 10 (false) -> REJECTED.
        let mut parser_accept = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        let read_accept = SamRecord {
            name: b"read_accept".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(39)),
            alignment_end: Some(Pos1(68)),
            mapping_quality: Some(60),
            cigar: vec![(30, CigarKind::Match)],
            sequence: vec![b'A'; 30],
            quality_scores: vec![60; 30],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };
        parser_accept.parse_cigar(&read_accept).unwrap();
        assert_eq!(
            parser_accept.ref_coverage.get(&39),
            Some(&1),
            "read at 1-based start=39 satisfies 50 - 39 > 10 and must be kept"
        );

        // Case 2: start = 30, rlen3 = 30.
        // Java: start + rlen3 - cutSite = 30 + 30 - 50 = 10 > 10 (false) -> SKIPPED,
        //       but maxReadLength is updated before skipSitesOutRegionOfInterest (:309-318).
        // Buggy +1 offset gives 31 + 30 - 50 = 11 > 10 (true) -> ACCEPTED.
        let mut parser_reject = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        let read_reject = SamRecord {
            name: b"read_reject".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(30)),
            alignment_end: Some(Pos1(59)),
            mapping_quality: Some(60),
            cigar: vec![(30, CigarKind::Match)],
            sequence: vec![b'A'; 30],
            quality_scores: vec![60; 30],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };
        parser_reject.parse_cigar(&read_reject).unwrap();
        assert!(
            parser_reject.ref_coverage.is_empty(),
            "read at 1-based start=30 has 30 + 30 - 50 == 10 (not > 10) and must be skipped"
        );
        assert_eq!(
            parser_reject.max_read_length, 30,
            "CigarParser.java:309-318 updates maxReadLength before skipSitesOutRegionOfInterest"
        );

        // Case 3: Spliced CIGAR (15M20N15M at start=30, end=79).
        // Java ALIGNED_LENGTH_MD ([MD] only) gives rlen3 = 30 -> 30 + 30 - 50 = 10 -> SKIPPED.
        // Including N (get_reference_length = 50) gives 30 + 50 - 50 = 30 > 10 -> ACCEPTED.
        let mut parser_spliced = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        let read_spliced = SamRecord {
            alignment_end: Some(Pos1(79)),
            cigar: vec![
                (15, CigarKind::Match),
                (20, CigarKind::Skip),
                (15, CigarKind::Match),
            ],
            ..read_reject
        };
        parser_spliced.parse_cigar(&read_spliced).unwrap();
        assert!(
            parser_spliced.ref_coverage.is_empty(),
            "spliced read 15M20N15M has ALIGNED_LENGTH_MD = 30 (not 50) and must be skipped"
        );
        assert_eq!(parser_spliced.max_read_length, 30);
    }

    #[test]
    fn increment_and_inc_cnt_preserve_counts() {
        let mut counts: BTreeMap<i32, rustc_hash::FxHashMap<String, i32>> = BTreeMap::new();
        increment(&mut counts, 10, "+AT");
        increment(&mut counts, 10, "+AT");
        increment(&mut counts, 10, "+G");
        assert_eq!(counts.get(&10).unwrap().get("+AT"), Some(&2));
        assert_eq!(counts.get(&10).unwrap().get("+G"), Some(&1));

        let reference = Reference::new();
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        inc_cnt(&mut parser.ref_coverage, 42, 1);
        inc_cnt(&mut parser.ref_coverage, 42, 3);
        inc_cnt(&mut parser.ref_coverage, 43, 2);
        assert_eq!(parser.ref_coverage.get(&42), Some(&4));
        assert_eq!(parser.ref_coverage.get(&43), Some(&2));
        assert_eq!(parser.ref_coverage.get(&99), None);
    }

    #[test]
    fn parse_cigar_selects_secondary_loaded_region_for_partial_pipeline() {
        // Issue #42 (Bug A): when SV `partial_pipeline` loads a secondary sub-region
        // into `reference.loaded_regions[1]` / `reference.region_sequences[1]`,
        // `CigarParser::parse_cigar` must pass that secondary region's sequence and
        // `sequence_start` to `CigarModifier` instead of `loaded_regions[0]`.
        let mut reference = Reference::new();
        // Primary region at 1..=200
        let seq0 = vec![b'A'; 200];
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr20".to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
        reference
            .reference_sequences
            .insert("chr20".to_string(), seq0.clone());
        reference.region_sequences.push(seq0);

        // Secondary region at 1000..=1200 where 1050..1075 is "ACGTACGT..." with a 2-bp
        // soft-clip match at the 5' end that CigarModifier rewrites when given the
        // secondary region's reference bases.
        let mut seq1 = vec![b'T'; 201];
        // Positions 1050..=1051 (offsets 50..=51) are 'G', 'C'
        seq1[50] = b'G';
        seq1[51] = b'C';
        for (i, slot) in seq1[52..80].iter_mut().enumerate() {
            *slot = b"ACGT"[i % 4];
        }
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr20".to_string(),
                sequence_start: 1000,
                sequence_end: 1200,
            });
        reference.region_sequences.push(seq1.clone());

        let region = Region {
            chr: "chr20".to_string(),
            start: Pos1(1040),
            end: Pos1(1100),
            gene: "G".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, true, 0, 0, 0);

        // Read aligned at 1052 with 2S28M where the 2S bases ("GC") match 1050..=1051 in seq1.
        let mut read_seq = vec![b'G', b'C'];
        read_seq.extend_from_slice(&seq1[52..80]);
        let record = SamRecord {
            name: b"r_sec".to_vec(),
            flags: 0,
            reference_name: Some("chr20".to_string()),
            alignment_start: Some(Pos1(1052)),
            alignment_end: Some(Pos1(1079)),
            mapping_quality: Some(60),
            cigar: vec![(2, CigarKind::SoftClip), (28, CigarKind::Match)],
            sequence: read_seq,
            quality_scores: vec![60; 30],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };
        parser.parse_cigar(&record).unwrap();
        // Because CigarModifier matched the 2S ("GC") against 1050..=1051 in `region_sequences[1]`,
        // it rewrote 2S28M @ 1052 -> 30M @ 1050, adding reference coverage at 1050 and 1051.
        assert_eq!(
            parser.ref_coverage.get(&1050),
            Some(&1),
            "CigarModifier must use secondary loaded_region[1] reference sequence"
        );
        assert_eq!(parser.ref_coverage.get(&1051), Some(&1));

        // Also verify a spliced read (`15M50N13M2S` at 1052 -> 1130) when an earlier
        // secondary region covers `read_start + qlen` (`1000..=1090`) but only a wider
        // secondary region covers the true reference-consuming `read_end` (`1000..=1200`).
        let mut ref_spliced = Reference::new();
        ref_spliced
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr20".to_string(),
                sequence_start: 1000,
                sequence_end: 1090,
            });
        ref_spliced.region_sequences.push(vec![b'T'; 91]);
        let mut seq2 = vec![b'T'; 201];
        for (i, slot) in seq2[52..67].iter_mut().enumerate() {
            *slot = b"ACGT"[i % 4];
        }
        for (i, slot) in seq2[117..130].iter_mut().enumerate() {
            *slot = b"ACGT"[i % 4];
        }
        seq2[130] = b'G';
        seq2[131] = b'C';
        ref_spliced
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr20".to_string(),
                sequence_start: 1000,
                sequence_end: 1200,
            });
        ref_spliced.region_sequences.push(seq2.clone());

        let region_spliced = Region {
            chr: "chr20".to_string(),
            start: Pos1(1040),
            end: Pos1(1150),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut parser_spliced =
            CigarParser::new(&ref_spliced, &region_spliced, &config, true, 0, 0, 0);
        let mut spliced_seq = Vec::with_capacity(30);
        spliced_seq.extend_from_slice(&seq2[52..67]);
        spliced_seq.extend_from_slice(&seq2[117..130]);
        spliced_seq.extend_from_slice(b"GC"); // 2S matching 1130..=1131 in seq2
        let spliced_record = SamRecord {
            name: b"r_spliced".to_vec(),
            flags: 0,
            reference_name: Some("chr20".to_string()),
            alignment_start: Some(Pos1(1052)),
            alignment_end: Some(Pos1(1129)),
            mapping_quality: Some(60),
            cigar: vec![
                (15, CigarKind::Match),
                (50, CigarKind::Skip),
                (13, CigarKind::Match),
                (2, CigarKind::SoftClip),
            ],
            sequence: spliced_seq,
            quality_scores: vec![60; 30],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };
        parser_spliced.parse_cigar(&spliced_record).unwrap();
        assert_eq!(
            parser_spliced.ref_coverage.get(&1130),
            Some(&1),
            "spliced read spanning 50N past 1090 must select the 1000..=1200 region slice"
        );
        assert_eq!(parser_spliced.ref_coverage.get(&1131), Some(&1));
    }

    #[test]
    fn parse_cigar_handles_empty_and_truncated_quality_scores() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        reference
            .region_sequences
            .push(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec());
        reference.reference_sequences.insert(
            "chr1".to_string(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        );
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        config.base_qual = 0.0;
        let mut parser = CigarParser::new(&reference, &region, &config, true, 0, 0, 0);

        // 1. Completely empty quality_scores (`*` in SAM) with 5' softclip, match/mismatch, indel, and 3' softclip
        let rec_empty_qual = SamRecord {
            name: b"r_empty_qual".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(25)),
            mapping_quality: Some(60),
            cigar: vec![
                (4, CigarKind::SoftClip),
                (5, CigarKind::Match),
                (2, CigarKind::Deletion),
                (1, CigarKind::Match),
                (2, CigarKind::Insertion),
                (4, CigarKind::Match),
                (4, CigarKind::SoftClip),
            ],
            sequence: b"TTTTCGTACATTCGGGGGGG".to_vec(),
            quality_scores: Vec::new(),
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(4),
            sa: None,
            mc: None,
            mq: None,
        };
        parser.parse_cigar(&rec_empty_qual).unwrap();

        // 2. Truncated quality_scores shorter than sequence across 5' and 3' softclips and sub-33 bytes
        let rec_trunc_qual = SamRecord {
            name: b"r_trunc_qual".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(20)),
            alignment_end: Some(Pos1(30)),
            mapping_quality: Some(60),
            cigar: vec![
                (5, CigarKind::SoftClip),
                (6, CigarKind::Match),
                (5, CigarKind::SoftClip),
            ],
            sequence: b"GGGGGCGTACGTTTTT".to_vec(),
            // Only 12 quality bytes for 16-base sequence (so 3' softclip has 1 high-qual + missing bytes)
            // plus a sub-33 byte (10) to verify saturation to Phred 0.
            quality_scores: vec![70, 70, 70, 70, 10, 70, 70, 70, 70, 70, 70, 70],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(0),
            sa: None,
            mc: None,
            mq: None,
        };
        parser.parse_cigar(&rec_trunc_qual).unwrap();
        assert_eq!(qual_phred(&[10u8], 0), 0.0);
        assert_eq!(qual_phred(&Vec::<u8>::new(), 5), 0.0);
    }

    #[test]
    fn chimeric_realignment_debug_log_gated_solely_on_verbose() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chrZ".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        reference
            .reference_sequences
            .insert("chrZ".to_string(), vec![b'T'; 100]);
        let record_sa = chimeric_seed_test_record(Some("chrZ,10,-,5M24S,60,0".to_string()));
        let record_seed = chimeric_seed_test_record(None);
        let seq_tail = get_reverse_complemented_sequence(&record_seed, -24, 24);
        let seed = seq_tail[seq_tail.len() - 17..].to_vec();
        reference.seed.insert(seed, vec![30]);
        let region = Region {
            chr: "chrZ".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "T".to_string(),
            amplicon: None,
        };

        for verbose in [false, true] {
            let mut config = Config::default();
            config.verbose = verbose;
            let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 100, 100);
            assert!(parser.is_read_chimeric_with_sa(
                &record_sa,
                1,
                "chrZ,10,-,5M24S,60,0",
                false,
                false,
            ));
            parser.parse_cigar(&record_sa).unwrap();
            assert!(parser.soft_clips_3_end.is_empty());
            parser.parse_cigar(&record_seed).unwrap();
            assert!(parser.soft_clips_3_end.is_empty());
        }
    }

    #[test]
    fn supplementary_read_updates_max_read_length_before_early_return() {
        // Issue #60 (#18): CigarParser.java:307-315 updates `maxReadLength` from
        // `totalLengthIncludingSoftClipped` before checking `(flags & 0x800) != 0`.
        // Observed on HG008 Ultima UG100 chr20:323309 where a 317 bp supplementary
        // read increases maxReadLength from 287 to 317.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr20".to_string(),
                sequence_start: 1,
                sequence_end: 500,
            });
        let ref_seq = vec![b'A'; 500];
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr20".to_string(), ref_seq);
        let region = Region {
            chr: "chr20".to_string(),
            start: Pos1(1),
            end: Pos1(500),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 500, 287);
        let rec = SamRecord {
            name: b"ug100_supp_317bp".to_vec(),
            flags: 0x800, // supplementary
            reference_name: Some("chr20".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(326)),
            mapping_quality: Some(60),
            cigar: vec![(317, CigarKind::Match)],
            sequence: vec![b'A'; 317],
            quality_scores: vec![73; 317],
            ..SamRecord::default()
        };
        parser.parse_cigar(&rec).unwrap();
        assert_eq!(
            parser.max_read_length, 317,
            "supplementary read must update max_read_length before early return"
        );
    }

    #[test]
    fn m_after_i_n_reference_base_counts_as_mismatch_and_advances_offset() {
        // Issue #60 (#4): CigarParser.java:1680-1695 (`isCloserThenVextAndGoodBase`)
        // compares `q != refCh` directly when `refCh == 'N'`, setting `offset = vi + 1`
        // and incrementing `nmoff` so the `M` after `I` merges into a complex variant
        // instead of breaking early on `N`. Observed on HG002 MGI DNBSEQ chr20:31059956.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr20".to_string(),
                sequence_start: 1,
                sequence_end: 60,
            });
        let mut ref_seq = vec![b'A'; 60];
        // At 1-based pos 11 (index 10), place 'N' right after the 1I at pos 10.
        ref_seq[10] = b'N';
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr20".to_string(), ref_seq);
        let region = Region {
            chr: "chr20".to_string(),
            start: Pos1(1),
            end: Pos1(50),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 1, 60);
        // 10M (pos 1..10) + 1I ('C') + 10M (pos 11..20, where ref[11] == 'N' and read has 'G')
        let mut seq = vec![b'A'; 21];
        seq[10] = b'C'; // 1I
        seq[11] = b'G'; // pos 11 against ref 'N'
        let rec = SamRecord {
            name: b"mgi_n_after_ins".to_vec(),
            flags: 0,
            reference_name: Some("chr20".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(20)),
            mapping_quality: Some(60),
            cigar: vec![
                (10, CigarKind::Match),
                (1, CigarKind::Insertion),
                (10, CigarKind::Match),
            ],
            sequence: seq,
            quality_scores: vec![73; 21],
            nm: Some(2),
            ..SamRecord::default()
        };
        parser.parse_cigar(&rec).unwrap();
        let ins_at_10 = parser
            .insertion_variants
            .get(&10)
            .expect("expected insertion map at pos 10");
        assert!(
            ins_at_10.map.contains_key("+C&G"),
            "expected +C&G complex insertion spanning ref 'N' at pos 11, got keys: {:?}",
            ins_at_10.map.keys().collect::<Vec<_>>()
        );
    }

    #[test]
    fn sclip3_rescue_negative_element_len_does_not_underflow() {
        // Issue #60 (#10, #12): CigarParser.java:611-616, 1238-1258 preserves `offset`
        // when `moffset == 0` on a second `I` and uses signed comparison in the 3'
        // soft-clip rescue loop (`for (int si = 0; si < cigarElementLength; si++)`).
        // Observed on HG002 Ultima UG100 chr20:13697661 and HG008 Roche SBX chr20:72357.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 60,
            });
        let mut ref_seq = vec![b'A'; 60];
        ref_seq[11] = b'T'; // 1-based pos 12
        ref_seq[24] = b'T'; // 1-based pos 25 (i = 14 of 16M, within extension_nuc = 2)
        ref_seq[37] = b'T'; // 1-based pos 38
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(50),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 100, 100);

        let mut seq = vec![b'A'; 42];
        seq[10] = b'G'; // first 2I ("GG")
        seq[11] = b'G';
        seq[13] = b'C'; // mismatch at pos 12 (vi = 1 after first 2I -> offset = 2)
        seq[26] = b'C'; // mismatch at pos 25 (i = 14 of 16M -> absorbed into second 2I)
        seq[28] = b'G'; // second 2I ("GG")
        seq[29] = b'G';
        seq[41] = b'C'; // trailing mismatch at pos 38
        let qual = vec![73; 42];
        let rec = SamRecord {
            name: b"r1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(37)),
            mapping_quality: Some(60),
            cigar: vec![
                (10, CigarKind::Match),
                (2, CigarKind::Insertion),
                (16, CigarKind::Match),
                (2, CigarKind::Insertion),
                (11, CigarKind::Match),
                (1, CigarKind::SoftClip),
            ],
            sequence: seq,
            quality_scores: qual,
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(3),
            sa: None,
            mc: None,
            mq: None,
        };
        parser.parse_cigar(&rec).unwrap();
        assert!(
            parser.soft_clips_3_end.is_empty(),
            "negative remaining 3' soft-clip element_len must not wrap around and populate soft_clips_3_end"
        );
    }

    fn make_non_repeating_ref(len: usize) -> Vec<u8> {
        let bases = b"ACGT";
        let mut s = Vec::with_capacity(len);
        let mut state: u64 = 0x1234_5678_9ABC_DEF0;
        for i in 0..len {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let mut b = bases[(state as usize) % 4];
            if i > 0 && s[i - 1] == b {
                b = bases[((state as usize) + 1) % 4];
            }
            s.push(b);
        }
        s
    }

    #[test]
    fn deletion_tp_after_skipped_intron_del_is_negative() {
        let mut ref_seq = make_non_repeating_ref(1000);
        // Ensure CigarModifier's 3'-end check (which counts 50D in ALIGNED_LENGTH_MND)
        // sees matching bases at 543..553 so 10M is not soft-clipped.
        let tail = ref_seq[493..503].to_vec();
        ref_seq[543..553].copy_from_slice(&tail);
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq.clone());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(99),
            end: Pos1(700),
            gene: "G".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 1, 1000);

        // 20M (101..=120) + 50D (skipped next to 300N) + 300N (121..=420) + 70M (421..=490) + 3D (491..=493) + 10M (494..=503)
        let mut seq = Vec::with_capacity(100);
        seq.extend_from_slice(&ref_seq[100..120]);
        seq.extend_from_slice(&ref_seq[420..490]);
        seq.extend_from_slice(&ref_seq[493..503]);

        let rec = SamRecord {
            name: b"j1_del".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(101)),
            alignment_end: Some(Pos1(553)),
            mapping_quality: Some(60),
            cigar: vec![
                (20, CigarKind::Match),
                (50, CigarKind::Deletion),
                (300, CigarKind::Skip),
                (70, CigarKind::Match),
                (3, CigarKind::Deletion),
                (10, CigarKind::Match),
            ],
            sequence: seq,
            quality_scores: vec![b'I'; 100],
            nm: Some(53),
            ..SamRecord::default()
        };
        parser.parse_cigar(&rec).unwrap();

        let vmap = parser
            .non_insertion_variants
            .get(&491)
            .expect("expected deletion variant at pos 491");
        let del_var = vmap.map.get("-3").expect("expected -3 deletion at pos 491");
        assert_eq!(del_var.pp, -40);
        assert_eq!(del_var.pmean, -40.0);
    }

    #[test]
    fn insertion_tp_after_skipped_intron_del_is_negative() {
        let ref_seq = make_non_repeating_ref(1000);
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq.clone());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(99),
            end: Pos1(700),
            gene: "G".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 1, 1000);

        // 20M (101..=120) + 50D (skipped next to 300N) + 300N (121..=420) + 60M (421..=480) + 2I + 28M (481..=508)
        let mut seq = Vec::with_capacity(110);
        seq.extend_from_slice(&ref_seq[100..120]);
        seq.extend_from_slice(&ref_seq[420..480]);
        // Pick 2I bases different from flanking reference bases so 3' shift does not move it
        let avoid1 = ref_seq[479];
        let avoid2 = ref_seq[480];
        let ins_base = [b'A', b'C', b'G', b'T']
            .into_iter()
            .find(|&b| b != avoid1 && b != avoid2)
            .unwrap();
        seq.push(ins_base);
        seq.push(ins_base);
        seq.extend_from_slice(&ref_seq[480..508]);

        let rec = SamRecord {
            name: b"j1_ins".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(101)),
            alignment_end: Some(Pos1(558)),
            mapping_quality: Some(60),
            cigar: vec![
                (20, CigarKind::Match),
                (50, CigarKind::Deletion),
                (300, CigarKind::Skip),
                (60, CigarKind::Match),
                (2, CigarKind::Insertion),
                (28, CigarKind::Match),
            ],
            sequence: seq,
            quality_scores: vec![b'I'; 110],
            nm: Some(52),
            ..SamRecord::default()
        };
        parser.parse_cigar(&rec).unwrap();

        let ins_key = format!("+{}{}", ins_base as char, ins_base as char);
        let vmap = parser
            .insertion_variants
            .get(&480)
            .expect("expected insertion variant at pos 480");
        let ins_var = vmap
            .map
            .get(&ins_key)
            .expect("expected +2I insertion at pos 480");
        assert_eq!(ins_var.pp, -20);
        assert_eq!(ins_var.pmean, -20.0);
    }

    #[test]
    fn tp_boundary_p_eq_rlen1_and_rlen1_plus_1() {
        let ref_seq = make_non_repeating_ref(1000);
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq.clone());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(99),
            end: Pos1(700),
            gene: "G".to_string(),
            amplicon: None,
        };
        let config = Config::default();

        // Case 1: p == rlen1 -> tp == 0 (20M + 10D + 300N + 70M + 3D + 10M => rlen1=100, p at 3D = 20+10+70 = 100)
        let mut parser0 = CigarParser::new(&reference, &region, &config, false, 0, 1, 1000);
        let mut seq0 = Vec::with_capacity(100);
        seq0.extend_from_slice(&ref_seq[100..120]);
        seq0.extend_from_slice(&ref_seq[420..490]);
        seq0.extend_from_slice(&ref_seq[493..503]);
        let rec0 = SamRecord {
            name: b"j1_p_eq_rlen1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(101)),
            alignment_end: Some(Pos1(513)),
            mapping_quality: Some(60),
            cigar: vec![
                (20, CigarKind::Match),
                (10, CigarKind::Deletion),
                (300, CigarKind::Skip),
                (70, CigarKind::Match),
                (3, CigarKind::Deletion),
                (10, CigarKind::Match),
            ],
            sequence: seq0,
            quality_scores: vec![b'I'; 100],
            nm: Some(13),
            ..SamRecord::default()
        };
        parser0.parse_cigar(&rec0).unwrap();
        let del0 = &parser0.non_insertion_variants.get(&491).unwrap().map["-3"];
        assert_eq!(del0.pp, 0);
        assert_eq!(del0.pmean, 0.0);

        // Case 2: p == rlen1 + 1 -> tp == -1 (20M + 11D + 300N + 70M + 3D + 10M => rlen1=100, p at 3D = 101)
        let mut parser1 = CigarParser::new(&reference, &region, &config, false, 0, 1, 1000);
        let mut seq1 = Vec::with_capacity(100);
        seq1.extend_from_slice(&ref_seq[100..120]);
        seq1.extend_from_slice(&ref_seq[420..490]);
        seq1.extend_from_slice(&ref_seq[493..503]);
        let rec1 = SamRecord {
            name: b"j1_p_eq_rlen1_plus_1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(101)),
            alignment_end: Some(Pos1(514)),
            mapping_quality: Some(60),
            cigar: vec![
                (20, CigarKind::Match),
                (11, CigarKind::Deletion),
                (300, CigarKind::Skip),
                (70, CigarKind::Match),
                (3, CigarKind::Deletion),
                (10, CigarKind::Match),
            ],
            sequence: seq1,
            quality_scores: vec![b'I'; 100],
            nm: Some(14),
            ..SamRecord::default()
        };
        parser1.parse_cigar(&rec1).unwrap();
        let del1 = &parser1.non_insertion_variants.get(&491).unwrap().map["-3"];
        assert_eq!(del1.pp, -1);
        assert_eq!(del1.pmean, -1.0);
    }

    #[test]
    fn eqx_cigar_excluded_from_read_length_helpers_and_sv_end() {
        assert_eq!(
            get_match_insertion_length(&[(100, CigarKind::SequenceMatch)]),
            0
        );
        let sclip_eqx = [
            (10, CigarKind::SoftClip),
            (50, CigarKind::SequenceMatch),
            (1, CigarKind::SequenceMismatch),
            (49, CigarKind::SequenceMatch),
        ];
        assert_eq!(get_soft_clipped_length(&sclip_eqx), 10);
        assert_eq!(query_consumed_len(&sclip_eqx), 110);
        assert_eq!(
            query_consumed_len(&[(268_435_455, CigarKind::Match); 17]),
            17 * 268_435_455u64
        );

        // Without =/X, helpers must match query_consumed_len / M+I sums
        let classic = [
            (10, CigarKind::SoftClip),
            (50, CigarKind::Match),
            (2, CigarKind::Insertion),
            (3, CigarKind::Deletion),
            (48, CigarKind::Match),
        ];
        assert_eq!(get_soft_clipped_length(&classic), 110);
        assert_eq!(get_match_insertion_length(&classic), 100);
        assert_eq!(query_consumed_len(&classic), 110);

        // Test 10=1X10= read through CigarParser:
        // get_reference_length consumes =/X (21 bp), while rlen1 == 0 and max_read_length == 0,
        // giving tp = -10 at 1X (pos 11) and SV end == start.
        let cigar_10eq_1x_10eq = [
            (10, CigarKind::SequenceMatch),
            (1, CigarKind::SequenceMismatch),
            (10, CigarKind::SequenceMatch),
        ];
        assert_eq!(get_reference_length(&cigar_10eq_1x_10eq), 21);
        assert_eq!(get_match_insertion_length(&cigar_10eq_1x_10eq), 0);
        assert_eq!(get_soft_clipped_length(&cigar_10eq_1x_10eq), 0);
        assert_eq!(query_consumed_len(&cigar_10eq_1x_10eq), 21);

        let ref_seq = vec![b'A'; 100];
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "G".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 1, 0);

        let mut seq = vec![b'A'; 21];
        seq[10] = b'T'; // 1X at 1-based pos 11 (read_pos_exc = 10)
        let rec = SamRecord {
            name: b"eqx_read".to_vec(),
            flags: 0x1, // paired with mate on chr2 to also exercise prepare_sv_structures_for_analysis
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(21)),
            mapping_quality: Some(60),
            mate_reference_name: Some("chr2".to_string()),
            mate_alignment_start: Some(Pos1(500)),
            cigar: cigar_10eq_1x_10eq.to_vec(),
            sequence: seq,
            quality_scores: vec![b'I'; 21],
            nm: Some(1),
            ..SamRecord::default()
        };
        parser.parse_cigar(&rec).unwrap();
        assert_eq!(parser.max_read_length, 0);
        assert_eq!(
            parser.sv_structures.svfusfend.get("chr2").copied(),
            Some(1),
            "SV end loop must exclude =/X, keeping end == start (1)"
        );
        let snv = &parser.non_insertion_variants.get(&11).unwrap().map["T"];
        assert_eq!(snv.pp, -10);
        assert_eq!(snv.pmean, -10.0);
    }

    #[test]
    fn high_phred_quality_mean_matches_java() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
        let ref_seq: Vec<u8> = b"ACGTACGTAC".iter().copied().cycle().take(200).collect();
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq.clone());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "T".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        config.realign = 0;

        for phred in [95u8, 100u8] {
            let q_byte = phred + 33; // 128 for Phred 95, 133 for Phred 100 (>= 0x80)
            let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 100, 100);

            // 1. 50M2I48M with raw Phred `phred` on the 2 inserted bases
            let mut ins_seq = Vec::new();
            ins_seq.extend_from_slice(&ref_seq[0..50]);
            ins_seq.extend_from_slice(b"GG");
            ins_seq.extend_from_slice(&ref_seq[50..98]);
            let mut ins_qual = vec![73u8; 100];
            ins_qual[50] = q_byte;
            ins_qual[51] = q_byte;
            let ins_rec = SamRecord {
                name: b"r_ins".to_vec(),
                flags: 0,
                reference_name: Some("chr1".to_string()),
                alignment_start: Some(Pos1(1)),
                alignment_end: Some(Pos1(98)),
                mapping_quality: Some(60),
                cigar: vec![
                    (50, CigarKind::Match),
                    (2, CigarKind::Insertion),
                    (48, CigarKind::Match),
                ],
                sequence: ins_seq,
                quality_scores: ins_qual,
                mate_reference_name: None,
                mate_alignment_start: None,
                template_length: 0,
                nm: Some(2),
                sa: None,
                mc: None,
                mq: None,
            };
            parser.parse_cigar(&ins_rec).unwrap();
            let ins_var = parser
                .insertion_variants
                .iter()
                .find_map(|(_, vm)| vm.map.get("+GG"))
                .expect("expected +GG insertion variant");
            assert_eq!(
                ins_var.qmean,
                f64::from(phred),
                "insertion qmean for Phred {phred} must equal {phred}.0, got {}",
                ins_var.qmean
            );

            // 2. 50M1D50M with raw Phred `phred` at the base preceding the deletion
            let mut del_seq = Vec::new();
            del_seq.extend_from_slice(&ref_seq[0..50]);
            del_seq.extend_from_slice(&ref_seq[51..101]);
            let mut del_qual = vec![73u8; 100];
            del_qual[49] = q_byte;
            let del_rec = SamRecord {
                name: b"r_del".to_vec(),
                flags: 0,
                reference_name: Some("chr1".to_string()),
                alignment_start: Some(Pos1(1)),
                alignment_end: Some(Pos1(101)),
                mapping_quality: Some(60),
                cigar: vec![
                    (50, CigarKind::Match),
                    (1, CigarKind::Deletion),
                    (50, CigarKind::Match),
                ],
                sequence: del_seq,
                quality_scores: del_qual,
                mate_reference_name: None,
                mate_alignment_start: None,
                template_length: 0,
                nm: Some(1),
                sa: None,
                mc: None,
                mq: None,
            };
            parser.parse_cigar(&del_rec).unwrap();
            let del_var = parser
                .non_insertion_variants
                .get(&51)
                .and_then(|vm| vm.map.get("-1"))
                .expect("expected -1 deletion variant at pos 51");
            assert_eq!(
                del_var.qmean,
                f64::from(phred),
                "deletion qmean for Phred {phred} must equal {phred}.0, got {}",
                del_var.qmean
            );
        }
    }

    #[test]
    fn split_at_non_ascii_no_panic() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 60,
            });
        let ref_seq = vec![b'A'; 60];
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(50),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 100, 100);

        let mut seq = vec![b'A'; 22];
        seq[9] = 0xc9; // non-ASCII mismatch within vext before 2I
        seq[10] = 0x80; // non-ASCII inserted base
        seq[11] = b'G';
        let qual = vec![128u8; 22];
        let rec = SamRecord {
            name: b"r_non_ascii".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(20)),
            mapping_quality: Some(60),
            cigar: vec![
                (10, CigarKind::Match),
                (2, CigarKind::Insertion),
                (10, CigarKind::Match),
            ],
            sequence: seq,
            quality_scores: qual,
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: Some(3),
            sa: None,
            mc: None,
            mq: None,
        };
        parser.parse_cigar(&rec).unwrap();
    }

    #[test]
    fn parse_amp_param_iterator_parity() {
        assert_eq!(parse_amp_param("15:0.85"), (15, 0.85));
        assert_eq!(parse_amp_param("15:0.85:extra"), (15, 0.85));
        assert_eq!(parse_amp_param("bad"), (10, 0.95));
        assert_eq!(parse_amp_param(""), (10, 0.95));
        assert_eq!(parse_amp_param("15:"), (10, 0.95));
        assert_eq!(parse_amp_param(":0.85"), (10, 0.95));
    }

    #[test]
    fn is_read_chimeric_with_sa_field_count_parity() {
        let reference = Reference::new();
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "T".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 100, 100);
        parser.max_read_length = 100;

        let rec = SamRecord {
            name: b"r_sa".to_vec(),
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(100)),
            alignment_end: Some(Pos1(199)),
            ..SamRecord::default()
        };

        // < 4 comma-separated fields -> false
        for short_sa in ["", "chr1", "chr1,110", "chr1,110,+"] {
            assert!(
                !parser.is_read_chimeric_with_sa(&rec, 100, short_sa, true, true),
                "expected false for SA tag with < 4 fields: {short_sa:?}"
            );
        }

        // >= 4 comma-separated fields -> parsed without Vec allocation
        assert!(parser.is_read_chimeric_with_sa(&rec, 100, "chr1,110,+,25S75M", true, true));
        assert!(parser.is_read_chimeric_with_sa(&rec, 100, "chr1,110,+,25S75M,60,0;", true, true));
        assert!(parser.is_read_chimeric_with_sa(&rec, 100, "chr1,110,+,75M25S,60,0;", true, false));
        assert!(!parser.is_read_chimeric_with_sa(&rec, 100, "chr1,bad,+,25S75M,60,0;", true, true));
        assert!(!parser.is_read_chimeric_with_sa(&rec, 100, "chr2,110,+,25S75M,60,0;", true, true));
    }

    #[test]
    fn get_ref_base_primary_slice_cache_and_fallbacks_match_reference() {
        let config = Config::default();

        // Case (a) + (b) + (d): Primary region `loaded_regions[0]` covers 10..=59 (50 bp)
        // with `region_sequences[0]` truncated by `SEED1_LEN` (17 bp -> 33 bytes stored,
        // covering 10..=42), and a secondary loaded region `loaded_regions[1]` covering
        // the truncated tail 40..=70.
        let mut ref_with_tail = Reference::new();
        let full_seq: Vec<u8> = (0..100).map(|i| b"ACGT"[i % 4]).collect();
        ref_with_tail
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 10,
                sequence_end: 59,
            });
        // Truncate primary sequence by 17 bp (50 - 17 = 33 bytes -> covers 10..=42)
        ref_with_tail
            .region_sequences
            .push(full_seq[9..9 + 33].to_vec());
        // Secondary loaded region covering 40..=70 (including the 17-bp tail 43..=59)
        ref_with_tail
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 40,
                sequence_end: 70,
            });
        ref_with_tail
            .region_sequences
            .push(full_seq[39..70].to_vec());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(10),
            end: Pos1(59),
            gene: "G".to_string(),
            amplicon: None,
        };
        let parser = CigarParser::new(&ref_with_tail, &region, &config, false, 0, 0, 0);
        assert_eq!(parser.primary_ref_start, 10);
        assert_eq!(parser.primary_ref_seq.len(), 33);

        // Verify exact parity across out-of-bounds (<10, >70), in-bounds primary (10..=42),
        // and the 17-bp SEED1_LEN truncated tail (43..=59) + secondary region extension (60..=70).
        for pos in 0..=80u32 {
            assert_eq!(
                parser.get_ref_base(pos),
                ref_with_tail.get_base(&region.chr, pos),
                "mismatch at pos={pos}"
            );
        }

        // Case (c) + (d): Legacy `reference_sequences`-only test setup (`loaded_regions` empty
        // or without `region_sequences`).
        let mut ref_legacy = Reference::new();
        ref_legacy
            .reference_sequences
            .insert("chr1".to_string(), b"ACGTNacgtn".to_vec());
        let region_legacy = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(10),
            gene: "G".to_string(),
            amplicon: None,
        };
        let parser_legacy = CigarParser::new(&ref_legacy, &region_legacy, &config, false, 0, 0, 0);
        assert!(parser_legacy.primary_ref_seq.is_empty());
        for pos in 0..=15u32 {
            assert_eq!(
                parser_legacy.get_ref_base(pos),
                ref_legacy.get_base(&region_legacy.chr, pos),
                "legacy mismatch at pos={pos}"
            );
        }
    }

    #[test]
    fn reverse_complement_seed17_matches_full_rc_head_and_tail_over_17bp() {
        // 25-bp 5' clip + 10M + 25-bp 3' clip = 60 bp total, asymmetric sequences
        // so the first 17 bp and last 17 bp of each reverse-complemented 25-bp clip differ.
        let clip_5p = b"AAACCCGGGTTTACGTACGTACGTA"; // 25 bp
        let mid_10m = b"TTTTTTTTTT"; // 10 bp
        let clip_3p = b"CGATCGATCGATGGGCCCAAATTTA"; // 25 bp
        let mut seq = Vec::with_capacity(60);
        seq.extend_from_slice(clip_5p);
        seq.extend_from_slice(mid_10m);
        seq.extend_from_slice(clip_3p);

        let record = SamRecord {
            name: b"seed17_25bp".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(26)),
            alignment_end: Some(Pos1(35)),
            mapping_quality: Some(60),
            cigar: vec![
                (25, CigarKind::SoftClip),
                (10, CigarKind::Match),
                (25, CigarKind::SoftClip),
            ],
            sequence: seq,
            quality_scores: vec![60; 60],
            ..SamRecord::default()
        };

        // 5' clip (ci == 0): start_index = 0, length = 25, take_rc_head = true -> full_rc[..17]
        let full_5p = get_reverse_complemented_sequence(&record, 0, 25);
        assert_eq!(full_5p.len(), 25);
        assert_ne!(full_5p.get(..17), full_5p.get(25 - 17..));
        let seed_5p = reverse_complement_seed17(&record, 0, 25, true);
        assert_eq!(seed_5p.as_ref().map(|s| s.as_slice()), full_5p.get(..17));

        // 3' clip (ci == last): start_index = -25, length = 25, take_rc_head = false -> full_rc[len - 17..]
        let full_3p = get_reverse_complemented_sequence(&record, -25, 25);
        assert_eq!(full_3p.len(), 25);
        assert_ne!(full_3p.get(..17), full_3p.get(25 - 17..));
        let seed_3p = reverse_complement_seed17(&record, -25, 25, false);
        assert_eq!(
            seed_3p.as_ref().map(|s| s.as_slice()),
            full_3p.get(full_3p.len() - 17..)
        );

        // Edge cases: < 17 bp clip and out-of-bounds negative start_index return None
        assert_eq!(reverse_complement_seed17(&record, 0, 16, true), None);
        assert_eq!(reverse_complement_seed17(&record, -16, 16, false), None);
        assert_eq!(reverse_complement_seed17(&record, -100, 25, false), None);
    }

    #[test]
    fn has_two_soft_clips_and_2digit_s_len_match_regexes() {
        use crate::cigar::modify::{leading_2digit_s_len, trailing_2digit_s_len};

        let re_mc = regex::Regex::new(r"\d+S.*?\d+S").unwrap();
        let re_5 = regex::Regex::new(r"^(\d\d+)S").unwrap();
        let re_3 = regex::Regex::new(r"(\d\d+)S$").unwrap();

        assert_eq!(trailing_2digit_s_len("80M125S"), Some(125));
        assert_eq!(leading_2digit_s_len("125S80M"), Some(125));
        assert_eq!(leading_2digit_s_len("9S80M"), None);
        assert_eq!(trailing_2digit_s_len("80M9S"), None);

        let samples = [
            "",
            "*",
            "100M",
            "5S95M",
            "95M5S",
            "10S80M",
            "80M10S",
            "80M125S",
            "125S80M",
            "10S80M15S",
            "1S2S",
            "0S100M0S",
            "00S50M",
            "50M00S",
            "10SS",
            "S10S",
            "10H80M10S",
            "10S20M5I30M12S",
            "10S\n20S",
            "10S\n20S30S",
            "99999999999S10M",
            "10M99999999999S",
            "-12S10M",
            "10M-12S",
        ];
        for s in samples {
            assert_eq!(
                has_two_soft_clips(s),
                re_mc.is_match(s),
                "has_two_soft_clips mismatch on {s:?}"
            );
            let exp_5 = re_5
                .captures(s)
                .map(|c| c[1].parse::<i32>().unwrap_or(0) as usize);
            assert_eq!(
                leading_2digit_s_len(s),
                exp_5,
                "leading_2digit_s_len mismatch on {s:?}"
            );
            let exp_3 = re_3
                .captures(s)
                .map(|c| c[1].parse::<i32>().unwrap_or(0) as usize);
            assert_eq!(
                trailing_2digit_s_len(s),
                exp_3,
                "trailing_2digit_s_len mismatch on {s:?}"
            );
        }
    }

    #[test]
    fn parse_cigar_tlen_i32_min_crispr_sa_and_scratch_restores() {
        let mut reference = Reference::new();
        let ref_seq = b"AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTACAAAAAAAAAAAGAAAAAAAAAAGAAAAAAAAAGAAAAAAGAAAAAAAAAGAAAGAAAGAAAAAGAAAAAACAGAAAGAAGAATGAAAATCTGGGTG".to_vec();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: ref_seq.len() as u32,
            });
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq.clone());

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(ref_seq.len() as u32),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        config.minmatch = 25;
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 150);

        // 1. template_length = i32::MIN with query_qual.len() <= 15 (hits :2471 and :2598)
        let short_tlen_min = SamRecord {
            name: b"tlen_min_short".to_vec(),
            flags: 0x1 | 0x10, // paired, reverse (read_dir_num = -1, mate_dir_num = 1)
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(39)),
            mapping_quality: Some(60),
            mate_reference_name: Some("chr1".to_string()),
            mate_alignment_start: Some(Pos1(50)),
            template_length: i32::MIN,
            cigar: vec![(30, CigarKind::Match)],
            sequence: ref_seq[9..39].to_vec(),
            quality_scores: vec![b'I'; 12], // <= 15 -> evaluates :2471 then falls into :2598
            nm: Some(0),
            ..SamRecord::default()
        };
        parser.parse_cigar(&short_tlen_min).unwrap();

        // Also template_length = i32::MIN with query_qual.len() > 15 across deletion/duplication/inversion
        for flags in [0x1 | 0x10, 0x1 | 0x20, 0x1] {
            let rec = SamRecord {
                quality_scores: vec![b'I'; 30],
                flags,
                ..short_tlen_min.clone()
            };
            parser.parse_cigar(&rec).unwrap();
        }

        // Amplicon mode with template_length = i32::MIN hitting :3258 (segstart - template_length - 1)
        let mut amp_cfg = config.clone();
        amp_cfg.amplicon = Some("10:0.95".to_string());
        let mut amp_parser = CigarParser::new(&reference, &region, &amp_cfg, false, 0, 0, 150);
        amp_parser.parse_cigar(&short_tlen_min).unwrap();

        // 2. SA:Z:chr1,2147483647,+,20S80M,60,0; (and reverse extreme position difference)
        let mut sa_seq = vec![b'A'; 100];
        sa_seq[20..100].copy_from_slice(&ref_seq[0..80]);
        let sa_rec = SamRecord {
            name: b"sa_i32_max".to_vec(),
            flags: 0x10,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(80)),
            mapping_quality: Some(60),
            cigar: vec![(20, CigarKind::SoftClip), (80, CigarKind::Match)],
            sequence: sa_seq,
            quality_scores: vec![b'I'; 100],
            nm: Some(0),
            sa: Some("chr1,2147483647,+,20S80M,60,0;".to_string()),
            ..SamRecord::default()
        };
        parser.parse_cigar(&sa_rec).unwrap();
        assert!(!parser.is_read_chimeric_with_sa(
            &sa_rec,
            -10,
            "chr1,2147483647,+,20S80M,60,0;",
            true,
            true,
        ));

        // 3. Warm up all 7 scratch buffers (including scratch_cigar_str, scratch_s, scratch_ss, scratch_splice_key)
        let mut warmup_seq = ref_seq[0..30].to_vec();
        warmup_seq[5] = b'T';
        warmup_seq[6] = b'G';
        let warmup_rec = SamRecord {
            name: b"warmup".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(40)),
            mapping_quality: Some(60),
            cigar: vec![
                (15, CigarKind::Match),
                (10, CigarKind::Skip),
                (15, CigarKind::Match),
            ],
            sequence: warmup_seq,
            quality_scores: vec![b'I'; 30],
            nm: Some(2),
            ..SamRecord::default()
        };
        parser.parse_cigar(&warmup_rec).unwrap();
        let assert_all_7_scratch_populated = |p: &CigarParser<'_>, ctx: &str| {
            assert!(p.scratch_seq.capacity() > 0, "{ctx}: scratch_seq");
            assert!(p.scratch_qual.capacity() > 0, "{ctx}: scratch_qual");
            assert!(
                p.scratch_cigar_str.capacity() > 0,
                "{ctx}: scratch_cigar_str"
            );
            assert!(p.scratch_elements.capacity() > 0, "{ctx}: scratch_elements");
            assert!(p.scratch_s.capacity() > 0, "{ctx}: scratch_s");
            assert!(p.scratch_ss.capacity() > 0, "{ctx}: scratch_ss");
            assert!(
                p.scratch_splice_key.capacity() > 0,
                "{ctx}: scratch_splice_key"
            );
        };
        assert_all_7_scratch_populated(&parser, "after warmup");

        // 3a. Early return when modify_cigar_preserving_buffers returns Err (PR #402 malformed rewrite)
        let pr402_seq = b"AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTACAAAAAAAAAAAAAAAAAAAAGAAAAAAGAAAAAAAAAG".to_vec();
        let pr402_qual = b"BBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF=74??44444444FFFFFFFFFFFFFFFFFFFFFFFFFFFFF".to_vec();
        let err_rec = SamRecord {
            name: b"pr402_err".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(155)),
            mapping_quality: Some(60),
            cigar: crate::cigar::parse_cigar("104M1I2M2D1M1D5M1D2M9D27M").unwrap(),
            sequence: pr402_seq,
            quality_scores: pr402_qual,
            nm: Some(18),
            ..SamRecord::default()
        };
        parser.parse_cigar(&err_rec).unwrap();
        assert_all_7_scratch_populated(&parser, "after modify_cigar Err early return");

        // 3b. Early return on minmatch filter (:425-427)
        let short_match_rec = SamRecord {
            name: b"short_minmatch".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(1)),
            alignment_end: Some(Pos1(10)),
            mapping_quality: Some(60),
            cigar: vec![(10, CigarKind::Match)], // < config.minmatch (25)
            sequence: ref_seq[0..10].to_vec(),
            quality_scores: vec![b'I'; 10],
            nm: Some(0),
            ..SamRecord::default()
        };
        parser.parse_cigar(&short_match_rec).unwrap();
        assert_all_7_scratch_populated(&parser, "after minmatch early return");

        // 4. CRISPR cut site near 2e9 (:452-456) without overflow AND preserving scratch buffers
        let mut crispr_ref = Reference::new();
        let crispr_bases = vec![b'A'; 500];
        crispr_ref
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 2_000_000_001,
                sequence_end: 2_000_000_500,
            });
        crispr_ref.region_sequences.push(crispr_bases);
        let crispr_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(2_000_000_001),
            end: Pos1(2_000_000_500),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut crispr_cfg = Config::default();
        crispr_cfg.crispr_cutting_site = Some(2_000_000_250);
        crispr_cfg.crispr_flanking_size = Some(50);
        let mut crispr_parser =
            CigarParser::new(&crispr_ref, &crispr_region, &crispr_cfg, false, 0, 0, 150);
        crispr_parser.scratch_seq = Vec::with_capacity(128);
        crispr_parser.scratch_qual = Vec::with_capacity(128);
        crispr_parser.scratch_cigar_str = String::with_capacity(32);
        crispr_parser.scratch_elements = Vec::with_capacity(8);
        crispr_parser.scratch_s = String::with_capacity(16);
        crispr_parser.scratch_ss = String::with_capacity(16);
        crispr_parser.scratch_splice_key = String::with_capacity(16);

        // Read with huge Deletion (200_000_000D) starting near 2e9 so start_pos + rlen3 > i32::MAX,
        // and filtered out by CRISPR site check at :452-456.
        let crispr_rec = SamRecord {
            name: b"crispr_2e9".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(2_000_000_240)),
            alignment_end: Some(Pos1(2_147_483_647)),
            mapping_quality: Some(60),
            cigar: vec![
                (20, CigarKind::Match),
                (200_000_000, CigarKind::Deletion),
                (20, CigarKind::Match),
            ],
            sequence: vec![b'A'; 40],
            quality_scores: vec![b'I'; 40],
            nm: None,
            ..SamRecord::default()
        };
        crispr_parser.parse_cigar(&crispr_rec).unwrap();
        assert_all_7_scratch_populated(&crispr_parser, "after CRISPR early return");
    }

    #[test]
    fn add_variation_for_deletion_caps_ref_coverage_on_corrupt_cigar() {
        // 1. Standard 5,000-bp loaded tile
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 5_000,
            });
        reference.region_sequences.push(vec![b'A'; 5_000]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(5_000),
            gene: "G".to_string(),
            amplicon: None,
        };
        let config = Config::default();
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, 150);

        parser.add_variation_for_deletion(
            60,
            0,
            false,
            100,
            "-268435455",
            b"I",
            0,
            100,
            10,
            268_435_455,
        );
        assert!(
            !parser.ref_coverage.is_empty(),
            "ref_coverage should record covered positions"
        );
        let min_pos = *parser.ref_coverage.keys().min().unwrap();
        let max_pos = *parser.ref_coverage.keys().max().unwrap();
        assert_eq!(min_pos, 100);
        assert_eq!(
            max_pos,
            5_000 + 1024,
            "ref_coverage must cap at loaded_end + 1024"
        );
        assert_eq!(parser.ref_coverage.len(), (5_000 + 1024 - 100 + 1) as usize);

        // 2. High-coordinate 5,000-bp loaded tile near 2e9 where start + i as i32 would wrap negative
        let mut hi_ref = Reference::new();
        hi_ref.loaded_regions.push(crate::reference::LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 2_000_000_001,
            sequence_end: 2_000_005_000,
        });
        hi_ref.region_sequences.push(vec![b'A'; 5_000]);
        let hi_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(2_000_000_001),
            end: Pos1(2_000_005_000),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut hi_parser = CigarParser::new(&hi_ref, &hi_region, &config, false, 0, 0, 150);
        hi_parser.add_variation_for_deletion(
            60,
            0,
            false,
            100,
            "-268435455",
            b"I",
            0,
            2_000_000_100,
            10,
            268_435_455,
        );
        assert!(
            hi_parser.ref_coverage.keys().all(|&pos| pos > 0),
            "no negative wrapped position allowed in ref_coverage"
        );
        assert_eq!(
            *hi_parser.ref_coverage.keys().max().unwrap(),
            2_000_005_000 + 1024
        );

        // 3. Exhaustive 256-byte parity of SINGLE_BASE_STR against push_ascii_base + N/n filter
        let mut sb = String::new();
        for b in 0u8..=255u8 {
            sb.clear();
            push_ascii_base(&mut sb, b);
            let mapped = SINGLE_BASE_STR.get(usize::from(b)).copied().unwrap_or("N");
            if sb.contains('N') || sb.contains('n') {
                assert_eq!(mapped, "N", "byte {b:#04x} must map to 'N'");
            } else {
                assert_eq!(
                    mapped,
                    sb.as_str(),
                    "byte {b:#04x} must match push_ascii_base"
                );
            }
        }

        // 4. Bounded lazy pileup pre-sizing (0 B on new / svflag=true; <= 4096 cap on first read)
        let mut wide_ref = Reference::new();
        wide_ref
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 10_000,
            });
        let wide_seq: Vec<u8> = (0..10_000).map(|i| b"ACGT"[i % 4]).collect();
        wide_ref.region_sequences.push(wide_seq.clone());
        let wide_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(10_000),
            gene: "G".to_string(),
            amplicon: None,
        };

        let mut sv_parser = CigarParser::new(&wide_ref, &wide_region, &config, true, 0, 0, 150);
        assert_eq!(sv_parser.ref_coverage.capacity(), 0);
        assert_eq!(sv_parser.non_insertion_variants.capacity(), 0);
        let mapped_30m = SamRecord {
            name: b"r30m".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(39)),
            mapping_quality: Some(60),
            cigar: vec![(30, CigarKind::Match)],
            sequence: wide_seq[9..39].to_vec(),
            quality_scores: vec![b'I'; 30],
            nm: Some(0),
            ..SamRecord::default()
        };
        sv_parser.parse_cigar(&mapped_30m).unwrap();
        assert!(
            sv_parser.ref_coverage.capacity() < 128,
            "svflag=true must NOT pre-reserve 4096-entry ref_coverage (got {})",
            sv_parser.ref_coverage.capacity()
        );

        let mut lazy_parser = CigarParser::new(&wide_ref, &wide_region, &config, false, 0, 0, 150);
        assert_eq!(lazy_parser.ref_coverage.capacity(), 0);
        assert_eq!(lazy_parser.non_insertion_variants.capacity(), 0);
        lazy_parser.parse_cigar(&mapped_30m).unwrap();
        assert!(
            lazy_parser.ref_coverage.capacity() >= 4_096
                && lazy_parser.ref_coverage.capacity() <= 8_192,
            "ref_coverage must lazily reserve bounded 4096 span cap, got {}",
            lazy_parser.ref_coverage.capacity()
        );
        // Zero-copy check on [Match] and pure S/M ([S,M], [M,S], [S,M,S]):
        assert_eq!(lazy_parser.scratch_seq.capacity(), 0);
        assert_eq!(lazy_parser.scratch_qual.capacity(), 0);
        assert_eq!(lazy_parser.scratch_cigar_str.capacity(), 0);

        // 5. Pure S/M in-place modifier zero-copy & exact pileup parity against fallback CigarModifier
        let mut fallback_cfg = Config::default();
        fallback_cfg.adaptor = vec!["ZZZZZZZZZZZZZZZZ".to_string()]; // forces PureSmOutcome::Fallback
        let mut fast_sm_parser =
            CigarParser::new(&wide_ref, &wide_region, &config, false, 0, 0, 150);
        let mut full_sm_parser =
            CigarParser::new(&wide_ref, &wide_region, &fallback_cfg, false, 0, 0, 150);

        let mut sm_seq1 = wide_seq[19..49].to_vec(); // Pos1(24) with 4S26M: 2 rescue + 2 clip + SNV
        sm_seq1[0] = b'T';
        sm_seq1[1] = b'T';
        sm_seq1[10] = if sm_seq1[10] == b'A' { b'C' } else { b'A' };
        let rec_sm = SamRecord {
            name: b"sm".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(24)),
            alignment_end: Some(Pos1(49)),
            mapping_quality: Some(60),
            cigar: vec![(4, CigarKind::SoftClip), (26, CigarKind::Match)],
            sequence: sm_seq1,
            quality_scores: vec![b'I'; 30],
            nm: Some(1),
            ..SamRecord::default()
        };

        let mut ms_seq = wide_seq[59..89].to_vec(); // Pos1(60) with 26M4S: 2 rescue + 2 clip
        ms_seq[28] = b'G';
        ms_seq[29] = b'G';
        let rec_ms = SamRecord {
            name: b"ms".to_vec(),
            flags: 0x10,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(60)),
            alignment_end: Some(Pos1(85)),
            mapping_quality: Some(60),
            cigar: vec![(26, CigarKind::Match), (4, CigarKind::SoftClip)],
            sequence: ms_seq,
            quality_scores: vec![b'I'; 30],
            nm: Some(0),
            ..SamRecord::default()
        };

        let mut sms_seq = wide_seq[99..129].to_vec(); // Pos1(104) with 4S22M4S + terminal mismatch trim
        sms_seq[0] = b'T';
        sms_seq[1] = b'T';
        sms_seq[2] = b'T';
        sms_seq[3] = b'T';
        sms_seq[4] = if sms_seq[4] == b'A' { b'C' } else { b'A' }; // 5' terminal mismatch -> clipped
        sms_seq[25] = if sms_seq[25] == b'A' { b'C' } else { b'A' }; // 3' terminal mismatch -> clipped
        sms_seq[26] = b'G';
        sms_seq[27] = b'G';
        sms_seq[28] = b'G';
        sms_seq[29] = b'G';
        let rec_sms = SamRecord {
            name: b"sms".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(104)),
            alignment_end: Some(Pos1(125)),
            mapping_quality: Some(60),
            cigar: vec![
                (4, CigarKind::SoftClip),
                (22, CigarKind::Match),
                (4, CigarKind::SoftClip),
            ],
            sequence: sms_seq,
            quality_scores: vec![b'I'; 30],
            nm: Some(2),
            ..SamRecord::default()
        };

        for r in [&rec_sm, &rec_ms, &rec_sms] {
            fast_sm_parser.parse_cigar(r).unwrap();
            full_sm_parser.parse_cigar(r).unwrap();
        }
        assert_eq!(
            fast_sm_parser.scratch_seq.capacity(),
            0,
            "pure S/M zero-copy path must never allocate scratch_seq"
        );
        assert_eq!(
            fast_sm_parser.scratch_qual.capacity(),
            0,
            "pure S/M zero-copy path must never allocate scratch_qual"
        );
        assert!(full_sm_parser.scratch_seq.capacity() > 0);
        assert_eq!(fast_sm_parser.ref_coverage, full_sm_parser.ref_coverage);
        let fast_vars: Vec<_> = fast_sm_parser
            .non_insertion_variants
            .iter()
            .map(|(&p, vm)| {
                (
                    p,
                    vm.map
                        .iter()
                        .map(|(k, v)| (k.to_owned(), v.cnt, v.hicnt, v.locnt, v.pmean, v.qmean))
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        let full_vars: Vec<_> = full_sm_parser
            .non_insertion_variants
            .iter()
            .map(|(&p, vm)| {
                (
                    p,
                    vm.map
                        .iter()
                        .map(|(k, v)| (k.to_owned(), v.cnt, v.hicnt, v.locnt, v.pmean, v.qmean))
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        assert_eq!(fast_vars, full_vars);

        // 6. Exact -T (trim_bases_after = 5) forward + reverse trimming with soft-clips
        let mut trim_ref = Reference::new();
        let trim_bases = b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec(); // 30 bp (Pos1(1..=30))
        trim_ref
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        trim_ref.region_sequences.push(trim_bases.clone());
        let trim_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(30),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut trim_cfg = Config::default();
        trim_cfg.trim_bases_after = 5;
        trim_cfg.minmatch = 0;
        trim_cfg.realign = 0; // test CigarParser's exact read_pos_inc offsets without modifier clip changes

        // 6a. Forward strand (!direction) with 2S10M at Pos1(5) (total_len = 12):
        // read_pos_inc across 10M runs 2..=11 -> dist = read_pos_inc <= 5 keeps read_pos_inc 2,3,4,5 (Pos1 5,6,7,8)
        // and trims read_pos_inc 6..=11 (Pos1 9..=14).
        let mut fwd_parser2 =
            CigarParser::new(&trim_ref, &trim_region, &trim_cfg, false, 0, 0, 150);
        let mut fwd_seq2 = vec![b'C', b'A']; // ref Pos1(3)='G', Pos1(4)='T' -> 0 walkback!
        fwd_seq2.extend_from_slice(&trim_bases[4..14]);
        fwd_seq2[5] = if fwd_seq2[5] == b'A' { b'C' } else { b'A' }; // read_pos_inc=5 (Pos1 8, kept)
        fwd_seq2[10] = if fwd_seq2[10] == b'A' { b'C' } else { b'A' }; // read_pos_inc=10 (Pos1 13, trimmed)
        let fwd_rec2 = SamRecord {
            name: b"fwd_t5".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(5)),
            alignment_end: Some(Pos1(14)),
            mapping_quality: Some(60),
            cigar: vec![(2, CigarKind::SoftClip), (10, CigarKind::Match)],
            sequence: fwd_seq2,
            quality_scores: vec![b'I'; 12],
            nm: Some(2),
            ..SamRecord::default()
        };
        fwd_parser2.parse_cigar(&fwd_rec2).unwrap();
        let mut fwd_cov2: Vec<i32> = fwd_parser2.ref_coverage.keys().copied().collect();
        fwd_cov2.sort_unstable();
        assert_eq!(
            fwd_cov2,
            vec![5, 6, 7, 8],
            "forward -T 5 with 2S10M must keep only read_pos_inc 2..=5 (Pos1 5..=8)"
        );
        assert!(fwd_parser2
            .non_insertion_variants
            .get(&8)
            .is_some_and(|vm| vm.map.len() == 1));
        assert!(fwd_parser2.non_insertion_variants.get(&9).is_none());
        assert!(fwd_parser2.non_insertion_variants.get(&13).is_none());

        // Also verify CigarParser.java:376 pre-MNV evaluation: an MNV starting at read_pos_inc=5 (<= 5)
        // extends into read_pos_inc=6 (Pos1 8..=9) because isTrimAtOptTBases is captured BEFORE MNV lookahead.
        let mut fwd_mnv_parser =
            CigarParser::new(&trim_ref, &trim_region, &trim_cfg, false, 0, 0, 150);
        let mut fwd_mnv_seq = vec![b'C', b'A'];
        fwd_mnv_seq.extend_from_slice(&trim_bases[4..14]);
        fwd_mnv_seq[5] = if fwd_mnv_seq[5] == b'A' { b'C' } else { b'A' };
        fwd_mnv_seq[6] = if fwd_mnv_seq[6] == b'A' { b'C' } else { b'A' };
        let fwd_mnv_rec = SamRecord {
            sequence: fwd_mnv_seq,
            ..fwd_rec2
        };
        fwd_mnv_parser.parse_cigar(&fwd_mnv_rec).unwrap();
        let mut fwd_mnv_cov: Vec<i32> = fwd_mnv_parser.ref_coverage.keys().copied().collect();
        fwd_mnv_cov.sort_unstable();
        assert_eq!(fwd_mnv_cov, vec![5, 6, 7, 8, 9]);

        // 6b. Reverse strand (direction) with 10M2S at Pos1(5) (total_len = 12):
        // dist = 12 - read_pos_inc <= 5 keeps read_pos_inc 7,8,9 (Pos1 12,13,14) and trims 0..=6 (Pos1 5..=11).
        let mut rev_parser = CigarParser::new(&trim_ref, &trim_region, &trim_cfg, false, 0, 0, 150);
        let mut rev_seq = trim_bases[4..14].to_vec(); // Pos1(5..=14)
        rev_seq.extend_from_slice(b"AA"); // 2S at Pos1(15)='G', Pos1(16)='T' -> 0 walk-forward!
        rev_seq[2] = if rev_seq[2] == b'A' { b'C' } else { b'A' }; // Pos1(7), dist=10 > 5 -> trimmed
        rev_seq[7] = if rev_seq[7] == b'A' { b'C' } else { b'A' }; // Pos1(12), dist=5 <= 5 -> kept
        let rev_rec = SamRecord {
            name: b"rev_t5".to_vec(),
            flags: 0x10,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(5)),
            alignment_end: Some(Pos1(14)),
            mapping_quality: Some(60),
            cigar: vec![(10, CigarKind::Match), (2, CigarKind::SoftClip)],
            sequence: rev_seq,
            quality_scores: vec![b'I'; 12],
            nm: Some(2),
            ..SamRecord::default()
        };
        rev_parser.parse_cigar(&rev_rec).unwrap();
        let mut rev_cov: Vec<i32> = rev_parser.ref_coverage.keys().copied().collect();
        rev_cov.sort_unstable();
        assert_eq!(
            rev_cov,
            vec![12, 13, 14],
            "reverse -T 5 with 10M2S must keep only read_pos_inc 7..=9 (Pos1 12..=14)"
        );
        assert!(rev_parser.non_insertion_variants.get(&7).is_none());
        assert!(rev_parser.non_insertion_variants.get(&11).is_none());
        assert!(rev_parser.non_insertion_variants.get(&12).is_some());

        // 6c. Extreme trim_bases_after = i32::MIN + negative extension_nuc = -1 + contig end i32::MAX safety
        let mut min_t_cfg = Config::default();
        min_t_cfg.trim_bases_after = i32::MIN;
        min_t_cfg.minmatch = 0;
        let mut min_t_parser =
            CigarParser::new(&trim_ref, &trim_region, &min_t_cfg, false, 0, 0, 150);
        min_t_parser.parse_cigar(&rev_rec).unwrap();
        assert!(min_t_parser.ref_coverage.is_empty());

        let mut end_ref = Reference::new();
        end_ref.loaded_regions.push(crate::reference::LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 2_147_483_640,
            sequence_end: 2_147_483_647,
        });
        end_ref.region_sequences.push(b"ACGTACGT".to_vec());
        let end_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(2_147_483_640),
            end: Pos1(2_147_483_647),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut neg_vext_cfg = Config::default();
        neg_vext_cfg.extension_nuc = -1;
        neg_vext_cfg.minmatch = 0;
        neg_vext_cfg.realign = 0;
        let mut end_parser =
            CigarParser::new(&end_ref, &end_region, &neg_vext_cfg, false, 0, 0, 150);
        let end_rec = SamRecord {
            name: b"contig_end_max".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(2_147_483_644)),
            alignment_end: Some(Pos1(2_147_483_647)),
            mapping_quality: Some(60),
            cigar: vec![
                (2, CigarKind::Match),
                (1, CigarKind::Insertion),
                (2, CigarKind::Match),
            ],
            sequence: b"TTTTT".to_vec(),
            quality_scores: vec![b'I'; 5],
            nm: Some(3),
            ..SamRecord::default()
        };
        end_parser.parse_cigar(&end_rec).unwrap();
        assert!(end_parser.ref_coverage.contains_key(&2_147_483_647));

        let mut neg_vext_realign_cfg = Config::default();
        neg_vext_realign_cfg.extension_nuc = -1;
        neg_vext_realign_cfg.minmatch = 0;
        let mut end_realign_parser = CigarParser::new(
            &end_ref,
            &end_region,
            &neg_vext_realign_cfg,
            false,
            0,
            0,
            150,
        );
        end_realign_parser.parse_cigar(&end_rec).unwrap();
    }
}
