//! CIGAR parse (pileup engine) — ports CigarParser.java.

use crate::bam::{get_cigar_string, get_mate_reference_name, SamRecord};
use crate::config::Config;
use crate::reference::Reference;
use crate::regions::Region;
use crate::sv::is_has_and_not_equals;
use crate::types::VardictError;
use crate::variations::{
    get_variation, get_variation_maybe_mut as get_variation_maybe, Mate, SVStructures, Sclip,
    Variation, VariationData, VariationMap,
};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use regex::Regex;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::LazyLock;

static MC_CLIP_CLIP: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d+S.*?\d+S").unwrap());

struct Offset {
    offset: i32,
    sequence: String,
    quality_sequence: String,
    offset_number_of_mismatches: i32,
}

/// CigarParser stores the running state for parsing CIGAR strings and building pileups.
/// Ports `modules/CigarParser.java`.
pub struct CigarParser<'a> {
    pub reference: &'a Reference,
    pub region: &'a Region,
    pub config: &'a Config,
    pub svflag: bool,

    pub non_insertion_variants: crate::javamap::JavaHashMap<i32, VariationMap>,
    pub insertion_variants: crate::javamap::JavaHashMap<i32, VariationMap>,
    // Iteration-order audit (W5.3, 2026-08-23): the String-keyed inner HashMaps
    // below never leak order into output or realignment decisions.
    // position_to_insertion_count/position_to_deletion_count/mnp are consumed
    // ONLY via fill_and_sort_tmp, whose (count,pos,desc) comparator is a total
    // order; splice_count DOES reach stdout (lib.rs -x print) but its VarDictJava
    // counterpart is also a HashMap iterated unordered (CigarParser.java:65,
    // :102-106), so record-set parity is unaffected. All four are `new HashMap<>()`
    // in CigarParser.java:62-65. The same audit covers SVStructures
    // svffus/svrfus/svfusfend/svfusrend (variations/mod.rs): keyed lookups and
    // per-list operations only (SVStructures.java:38-40).
    pub position_to_insertion_count: BTreeMap<i32, HashMap<String, i32>>,
    pub position_to_deletion_count: BTreeMap<i32, HashMap<String, i32>>,
    pub sv_structures: SVStructures,
    pub ref_coverage: BTreeMap<i32, i32>,
    pub soft_clips_5_end: BTreeMap<i32, Sclip>,
    pub soft_clips_3_end: BTreeMap<i32, Sclip>,
    pub max_read_length: i32,
    pub splice: HashSet<String>,
    pub mnp: BTreeMap<i32, HashMap<String, i32>>,
    pub splice_count: HashMap<String, [i32; 1]>,
    pub discordant_count: i32,
    pub total_reads: i32,
    pub duplicate_reads: i32,
}

fn get_adjusted_nm(record: &SamRecord) -> i32 {
    let insertion_deletion_length = get_insertion_deletion_length(&record.cigar);
    record
        .nm
        .map(|n| n.saturating_sub(insertion_deletion_length as i32))
        .unwrap_or(0)
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
    ) -> Self {
        Self {
            reference,
            region,
            config,
            svflag,
            non_insertion_variants: crate::javamap::JavaHashMap::new(),
            insertion_variants: crate::javamap::JavaHashMap::new(),
            position_to_insertion_count: BTreeMap::new(),
            position_to_deletion_count: BTreeMap::new(),
            sv_structures: SVStructures::new(),
            ref_coverage: BTreeMap::new(),
            soft_clips_5_end: BTreeMap::new(),
            soft_clips_3_end: BTreeMap::new(),
            max_read_length: 0,
            splice: HashSet::new(),
            mnp: BTreeMap::new(),
            splice_count: HashMap::new(),
            discordant_count: 0,
            total_reads,
            duplicate_reads,
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
                read_end += len as i32;
            }
        }
        let region_start = self.region.start.0 as i32;
        let region_end = self.region.end.0 as i32;
        if read_start > region_end || read_end <= region_start {
            return Ok(());
        }

        if self.config.verbose {
            eprintln!(
                "DEBUG_READ: name={} start={} end={} flags={:#x} cigar={}",
                String::from_utf8_lossy(&record.name),
                read_start,
                read_end - 1,
                record.flags,
                get_cigar_string(&record.cigar)
            );
        }

        // If supplementary alignment is present and flag filter is enabled, ignore it
        if self.config.flag_filter != 0 && (record.flags & 0x800) != 0 {
            return Ok(());
        }

        let direction = (record.flags & 0x10) != 0;
        let insertion_deletion_length = get_insertion_deletion_length(&record.cigar);

        // T045 mismatch gate: drop read if NM - (I + D) > max_mismatches
        if let Some(nm_val) = record.nm {
            let total_mismatches = nm_val.saturating_sub(insertion_deletion_length as i32);
            if total_mismatches > self.config.max_mismatches as i32 {
                return Ok(()); // Drop read
            }
        }

        let mut start_pos = record.alignment_start.map(|p| p.0 as i32).unwrap_or(0);
        let mut query_seq = record.sequence.clone();
        let mut query_qual = record.quality_scores.clone();
        let mut cigar_elements = record.cigar.clone();

        // Amplicon read filter (Java parseCigarWithAmpCase, CigarParser.java:255-258):
        // skip reads that don't map to the amplicon.
        if let Some(ref amp_param) = self.config.amplicon {
            let (distance_to_amplicon, overlap_fraction) = parse_amp_param(amp_param);
            let is_mate_ref_name_equal =
                record.reference_name.as_deref() == record.mate_reference_name.as_deref();
            if parse_cigar_with_amp_case(
                &cigar_elements,
                record,
                is_mate_ref_name_equal,
                self.region,
                distance_to_amplicon,
                overlap_fraction,
            ) {
                return Ok(());
            }
        }

        // T046: run CigarModifier if realignment is enabled
        if self.config.realign == 1 {
            let (ref_slice, ref_start) = if let Some(lr) = self
                .reference
                .loaded_regions
                .iter()
                .find(|lr| lr.chr == self.region.chr)
            {
                if let Some(seq) = self.reference.reference_sequences.get(&self.region.chr) {
                    (seq.as_slice(), lr.sequence_start as i32)
                } else {
                    (&[] as &[u8], 1)
                }
            } else {
                (&[] as &[u8], 1)
            };

            let cigar_str = get_cigar_string(&cigar_elements);
            let mut modifier = crate::cigar::modify::CigarModifier::new(
                start_pos,
                &cigar_str,
                &query_seq,
                &query_qual,
                ref_slice,
                ref_start,
                &self.reference.seed,
                insertion_deletion_length as i32,
                self.max_read_length,
                self.config,
            );
            // ports CigarModifier.java:619-791, CigarParser.java:80-95:
            // When modify_cigar returns an error (JavaThrows on malformed CIGAR with non-positive
            // operator lengths), emulate Java's catch block: log the exception and skip the read
            // cleanly without mutating start_pos, query_seq, or pileup structures.
            match modifier.modify_cigar() {
                Ok(modified) => {
                    start_pos = modified.position;
                    query_seq = modified.query_seq;
                    query_qual = modified.query_qual;
                    if let Some(parsed) = crate::cigar::parse_cigar(&modified.cigar) {
                        cigar_elements = parsed;
                    }
                }
                Err(e) => {
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
                    return Ok(());
                }
            }
        }

        // CRISPR mode: skip reads not overlapping the cutting site region
        if let Some(cut_site) = self.config.crispr_cutting_site {
            if cut_site != 0 {
                if let Some(filter_bp) = self.config.crispr_flanking_size {
                    if filter_bp != 0 {
                        let rlen3 = get_reference_length(&cigar_elements) as i32;
                        let start_1based = start_pos + 1;
                        if !(cut_site - start_1based > filter_bp
                            && start_1based + rlen3 - cut_site > filter_bp)
                        {
                            return Ok(());
                        }
                    }
                }
            }
        }

        // T046: cleanupCigar
        let mut sv_cigar_elements = cigar_elements.clone();
        cleanup_cigar(&mut sv_cigar_elements);

        // Determine discordant reads
        if get_mate_reference_name(record) != "=" {
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
                    return Ok(()); // Drop read
                }
            }
        }

        let read_len_match_ins = get_match_insertion_length(&cigar_elements);
        if self.config.minmatch != 0 && read_len_match_ins < self.config.minmatch {
            return Ok(()); // Drop read shorter than minmatch (Java -M)
        }
        let total_length_including_soft_clipped = get_soft_clipped_length(&cigar_elements);
        if total_length_including_soft_clipped > self.max_read_length as u32 {
            self.max_read_length = total_length_including_soft_clipped as i32;
        }

        let paired = (record.flags & 0x1) != 0;
        let mate_unmapped = (record.flags & 0x8) != 0;
        let mapping_quality = record.mapping_quality.unwrap_or(0);
        if !(paired && mate_unmapped) {
            if mapping_quality > 10 && !self.config.no_sv {
                self.prepare_sv_structures_for_analysis(
                    record,
                    &query_qual,
                    get_adjusted_nm(record),
                    direction,
                    (record.flags & 0x20) == 0,
                    start_pos,
                    total_length_including_soft_clipped as i32,
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

        while ci < cigar_elements.len() {
            if self.skip_overlapping_reads(
                record,
                start_pos,
                direction,
                mate_alignment_start,
                start,
            ) {
                break;
            }
            let (cigar_element_len, raw_operator) = cigar_elements[ci];
            let operator = if (ci == 0 || ci == cigar_elements.len() - 1)
                && raw_operator == CigarKind::Insertion
            {
                CigarKind::SoftClip
            } else {
                raw_operator
            };

            match operator {
                CigarKind::Skip => {
                    let key = format!("{}-{}", start - 1, start + cigar_element_len as i32 - 1);
                    self.splice.insert(key.clone());
                    let cnt = self.splice_count.entry(key).or_insert([0]);
                    cnt[0] += 1;
                    start += cigar_element_len as i32;
                    ci += 1;
                    offset = 0;
                }
                CigarKind::SoftClip => {
                    self.process_soft_clip(
                        record,
                        &query_seq,
                        &query_qual,
                        &cigar_elements,
                        direction,
                        read_len_match_ins as i32,
                        total_length_including_soft_clipped as i32,
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
                        &query_seq,
                        &query_qual,
                        &cigar_elements,
                        ci,
                        direction,
                        read_len_match_ins as i32,
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
                        &query_seq,
                        &query_qual,
                        &cigar_elements,
                        ci,
                        direction,
                        read_len_match_ins as i32,
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

                    let mut i = offset;
                    while i < cigar_element_len {
                        if read_pos_inc >= query_seq.len() {
                            break;
                        }
                        let ch1 = query_seq[read_pos_inc];

                        if ch1 == b'N' || ch1 == b'n' {
                            if self.config.include_n {
                                inc_cnt(&mut self.ref_coverage, start, 1);
                            }
                            start += 1;
                            read_pos_inc += 1;
                            read_pos_exc += 1;
                            i += 1;
                            continue;
                        }

                        let mut q = query_qual[read_pos_inc] as f64 - 33.0;
                        let mut qbases = 1;
                        let mut qibases = 0;
                        let mut ss = String::new();
                        let mut start_with_deletion = false;

                        // T048: MNV Growth Look-Ahead
                        while (start + 1) as u32 >= self.region.start.0
                            && (start + 1) as u32 <= self.region.end.0
                            && i + 1 < cigar_element_len
                            && q >= self.config.base_qual
                            && is_has_and_not_equals(
                                self.reference,
                                &self.region.chr,
                                start as u32,
                                query_seq[read_pos_inc],
                            )
                            && is_not_equals_n(
                                self.reference.get_base(&self.region.chr, start as u32),
                            )
                        {
                            let next_read_idx = read_pos_inc + 1;
                            if next_read_idx >= query_seq.len() {
                                break;
                            }
                            if (query_qual[next_read_idx] as f64 - 33.0)
                                < self.config.base_qual + 5.0
                            {
                                break;
                            }
                            let nuc = query_seq[next_read_idx];
                            if nuc == b'N' || nuc == b'n' {
                                break;
                            }
                            if let Some(next_ref_base) = self
                                .reference
                                .get_base(&self.region.chr, (start + 1) as u32)
                            {
                                if next_ref_base == b'N' || next_ref_base == b'n' {
                                    break;
                                }
                                if nuc != next_ref_base {
                                    ss.push(nuc as char);
                                    q += query_qual[next_read_idx] as f64 - 33.0;
                                    qbases += 1;
                                    read_pos_inc += 1;
                                    read_pos_exc += 1;
                                    i += 1;
                                    start += 1;
                                    nmoff += 1;
                                } else {
                                    let mut ssn = 0;
                                    let vext = self.config.extension_nuc as usize;
                                    for ssi in 1..=vext {
                                        if i + 1 + ssi as u32 >= cigar_element_len {
                                            break;
                                        }
                                        let lookahead_read_idx = read_pos_inc + 1 + ssi;
                                        if lookahead_read_idx < query_seq.len() {
                                            if let Some(lookahead_ref_base) =
                                                self.reference.get_base(
                                                    &self.region.chr,
                                                    (start + 1 + ssi as i32) as u32,
                                                )
                                            {
                                                if query_seq[lookahead_read_idx]
                                                    != lookahead_ref_base
                                                {
                                                    ssn = ssi + 1;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    if ssn == 0 {
                                        break;
                                    }
                                    let ssn_read_idx = read_pos_inc + ssn;
                                    if ssn_read_idx >= query_seq.len()
                                        || (query_qual[ssn_read_idx] as f64 - 33.0)
                                            < self.config.base_qual + 5.0
                                    {
                                        break;
                                    }
                                    for ssi in 1..=ssn {
                                        let idx = read_pos_inc + ssi;
                                        ss.push(query_seq[idx] as char);
                                        q += query_qual[idx] as f64 - 33.0;
                                        qbases += 1;
                                    }
                                    read_pos_inc += ssn;
                                    read_pos_exc += ssn;
                                    i += ssn as u32;
                                    start += ssn as i32;
                                }
                            } else {
                                break;
                            }
                        }

                        let mut s = (ch1 as char).to_string();
                        if !ss.is_empty() {
                            s.push('&');
                            s.push_str(&ss);
                        }

                        let mut ddlen = 0;
                        if self.is_closer_than_vext_and_good_base(
                            &query_seq,
                            &query_qual,
                            &cigar_elements,
                            ci,
                            i,
                            ss.len(),
                            CigarKind::Deletion,
                            start,
                            read_pos_inc,
                        ) {
                            while i + 1 < cigar_element_len {
                                s.push(query_seq[read_pos_inc + 1] as char);
                                q += query_qual[read_pos_inc + 1] as f64 - 33.0;
                                qbases += 1;
                                i += 1;
                                read_pos_inc += 1;
                                read_pos_exc += 1;
                                start += 1;
                            }
                            let s_clean = s.replacen('&', "", 1);
                            s = format!("-{}&{}", cigar_elements[ci + 1].0, s_clean);
                            start_with_deletion = true;
                            ddlen = cigar_elements[ci + 1].0 as i32;
                            ci += 1;

                            if self.is_two_insertions_ahead(&cigar_elements, ci) {
                                let next_next_len = cigar_elements[ci + 1].0 as usize;
                                if let Some(sub) = query_seq
                                    .get(read_pos_inc + 1..read_pos_inc + 1 + next_next_len)
                                {
                                    s.push('^');
                                    s.push_str(&String::from_utf8_lossy(sub));
                                }
                                for qi in 1..=next_next_len {
                                    q += query_qual[read_pos_inc + 1 + qi] as f64 - 33.0;
                                    qibases += 1;
                                }
                                read_pos_inc += next_next_len;
                                read_pos_exc += next_next_len;
                                ci += 1;
                            }

                            if self.is_next_after_num_matched(&cigar_elements, ci, 1) {
                                let next_next_next_len = cigar_elements[ci + 1].0 as usize;
                                let tpl = self.find_offset(
                                    start + ddlen + 1,
                                    read_pos_inc + 1,
                                    next_next_next_len,
                                    &query_seq,
                                    &query_qual,
                                );
                                let toffset = tpl.offset;
                                if toffset != 0 {
                                    moffset = toffset;
                                    nmoff += tpl.offset_number_of_mismatches;
                                    s.push('&');
                                    s.push_str(&tpl.sequence);
                                    for b in tpl.quality_sequence.bytes() {
                                        q += b as f64 - 33.0;
                                        qibases += 1;
                                    }
                                    for osi in 0..toffset {
                                        *self
                                            .ref_coverage
                                            .entry(start + ddlen + 1 + osi)
                                            .or_default() += 1;
                                    }
                                }
                            }
                        } else if self.is_closer_than_vext_and_good_base(
                            &query_seq,
                            &query_qual,
                            &cigar_elements,
                            ci,
                            i,
                            ss.len(),
                            CigarKind::Insertion,
                            start,
                            read_pos_inc,
                        ) {
                            while i + 1 < cigar_element_len {
                                s.push(query_seq[read_pos_inc + 1] as char);
                                q += query_qual[read_pos_inc + 1] as f64 - 33.0;
                                qbases += 1;
                                i += 1;
                                read_pos_inc += 1;
                                read_pos_exc += 1;
                                start += 1;
                            }
                            let s_clean = s.replacen('&', "", 1);
                            let next_len = cigar_elements[ci + 1].0 as usize;
                            let mut temp_s = s_clean;
                            if let Some(sub) =
                                query_seq.get(read_pos_inc + 1..read_pos_inc + 1 + next_len)
                            {
                                temp_s.push_str(&String::from_utf8_lossy(sub));
                            }
                            let (first_part, second_part) = temp_s.split_at(next_len);
                            s = format!("+{}&{}", first_part, second_part);

                            for qi in 1..=next_len {
                                q += query_qual[read_pos_inc + 1 + qi] as f64 - 33.0;
                                qibases += 1;
                            }
                            read_pos_inc += next_len;
                            read_pos_exc += next_len;
                            ci += 1;
                            qibases -= 1;
                            qbases += 1;
                        }

                        // Add matching/mismatched Variation
                        let pos = start - qbases as i32 + 1;

                        if pos as u32 >= self.region.start.0
                            && pos as u32 <= self.region.end.0
                            && !s.contains('N')
                            && !s.contains('n')
                        {
                            if self.config.verbose && pos == 11291401 {
                                eprintln!(
                                    "DEBUG_VAR_ADD: pos={} s={} dir={} name={}",
                                    pos,
                                    s,
                                    direction,
                                    String::from_utf8_lossy(&record.name)
                                );
                            }
                            self.add_variation_for_matching_part(
                                record.mapping_quality.unwrap_or(0) as i32,
                                get_adjusted_nm(record),
                                direction,
                                read_len_match_ins as i32,
                                nmoff,
                                &s,
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

                        if start_with_deletion {
                            start += ddlen;
                        }

                        if operator != CigarKind::Insertion {
                            start += 1;
                        }
                        if operator != CigarKind::Deletion {
                            read_pos_inc += 1;
                            read_pos_exc += 1;
                        }

                        if self.skip_overlapping_reads(
                            record,
                            start_pos,
                            direction,
                            mate_alignment_start,
                            start,
                        ) {
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
                        start += moffset;
                        read_pos_exc += moffset as usize;
                    } else {
                        offset = 0;
                    }
                    ci += 1;
                }
            }

            if start > self.region.end.0 as i32 {
                break;
            }
        }

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
        hv.nm += (nm - nmoff) as f64;
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
            inc_cnt(&mut self.ref_coverage, start - qi + 1, 1);
        }

        if start_with_deletion {
            increment(&mut self.position_to_deletion_count, pos, s);
            for qi in 1..ddlen {
                inc_cnt(&mut self.ref_coverage, start + qi, 1);
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
        if cigar_elements.len() > ci + 2 && cigar_elements[ci + 2].1 == CigarKind::HardClip {
            return false;
        }
        if self.config.realign != 1 {
            return false;
        }
        let cigar_element_len = cigar_elements[ci].0;
        let vext = self.config.extension_nuc as u32;

        if cigar_element_len - i > vext {
            return false;
        }
        if cigar_elements.len() <= ci + 1 {
            return false;
        }
        if cigar_elements[ci + 1].1 != target_kind {
            return false;
        }
        if !self.reference.loaded_regions.iter().any(|lr| {
            lr.chr == self.region.chr
                && (start as u32) >= lr.sequence_start
                && (start as u32) <= lr.sequence_end
        }) {
            return false;
        }
        if read_pos_inc >= query_seq.len() {
            return false;
        }
        let q = query_qual[read_pos_inc].saturating_sub(33) as f64;
        if q < self.config.base_qual {
            return false;
        }
        let is_mismatch =
            if let Some(ref_base) = self.reference.get_base(&self.region.chr, start as u32) {
                ss_len > 0 || query_seq[read_pos_inc] != ref_base
            } else {
                false
            };
        is_mismatch
    }

    fn is_two_insertions_ahead(&self, cigar_elements: &[(u32, CigarKind)], ci: usize) -> bool {
        cigar_elements.len() > ci + 1 && cigar_elements[ci + 1].1 == CigarKind::Insertion
    }

    fn is_next_after_num_matched(
        &self,
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        number: usize,
    ) -> bool {
        cigar_elements.len() > ci + number && cigar_elements[ci + number].1 == CigarKind::Match
    }

    fn is_insertion_or_deletion_with_next_matched(
        &self,
        cigar: &[(u32, CigarKind)],
        ci: usize,
    ) -> bool {
        if self.config.realign != 1 {
            return false;
        }
        let vext = self.config.extension_nuc as u32;
        if cigar.len() <= ci + 2 {
            return false;
        }
        if cigar[ci + 1].0 > vext || cigar[ci + 1].1 != CigarKind::Match {
            return false;
        }
        if cigar[ci + 2].1 != CigarKind::Insertion && cigar[ci + 2].1 != CigarKind::Deletion {
            return false;
        }
        if cigar.len() > ci + 3 {
            let op3 = cigar[ci + 3].1;
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
        quality_segment: &mut String,
        m_len: usize,
        indel_len: usize,
        begin: usize,
        is_insertion: bool,
    ) {
        let sub_m = crate::util::subseq(query_sequence, begin, m_len).unwrap_or(&[]);
        desc_string.push('#');
        desc_string.push_str(&String::from_utf8_lossy(sub_m));

        let qual_m = crate::util::subseq(query_quality, begin, m_len).unwrap_or(&[]);
        quality_segment.push_str(&String::from_utf8_lossy(qual_m));

        desc_string.push('^');
        let is_next_next_ins = cigar_elements[ci + 2].1 == CigarKind::Insertion;
        if is_next_next_ins {
            let sub_ins =
                crate::util::subseq(query_sequence, begin + m_len, indel_len).unwrap_or(&[]);
            desc_string.push_str(&String::from_utf8_lossy(sub_ins));
        } else {
            desc_string.push_str(&indel_len.to_string());
        }

        if is_insertion {
            if is_next_next_ins {
                let qual_ins =
                    crate::util::subseq(query_quality, begin + m_len, indel_len).unwrap_or(&[]);
                quality_segment.push_str(&String::from_utf8_lossy(qual_ins));
            } else {
                if let Some(q_char) = crate::util::char_at(query_quality, (begin + m_len) as isize)
                {
                    quality_segment.push(q_char as char);
                }
            }
        } else {
            if is_next_next_ins {
                let qual_ins =
                    crate::util::subseq(query_quality, begin + m_len, indel_len).unwrap_or(&[]);
                quality_segment.push_str(&String::from_utf8_lossy(qual_ins));
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
        quality_segment: &str,
        nmoff: i32,
        start: i32,
        read_pos_exc: usize,
        cigar_element_len: u32,
    ) {
        let hv = get_variation(&mut self.non_insertion_variants, start, desc_string);
        increment(&mut self.position_to_deletion_count, start, desc_string);

        hv.inc_dir(dir);
        hv.cnt += 1;

        let tp = if read_pos_exc < rlen1 as usize - read_pos_exc {
            read_pos_exc as i32 + 1
        } else {
            rlen1 - read_pos_exc as i32
        };

        let mut tmpq = 0.0;
        if !quality_segment.is_empty() {
            for b in quality_segment.bytes() {
                tmpq += (b as f64) - 33.0;
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
        hv.nm += (nm - nmoff) as f64;

        if tmpq >= self.config.base_qual {
            hv.hicnt += 1;
        } else {
            hv.locnt += 1;
        }

        for i in 0..cigar_element_len {
            inc_cnt(&mut self.ref_coverage, start + i as i32, 1);
        }
    }

    fn process_deletion(
        &mut self,
        record: &SamRecord,
        query_seq: &[u8],
        query_qual: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        direction: bool,
        read_len_match_ins: i32,
        start: &mut i32,
        read_pos_inc: &mut usize,
        read_pos_exc: &mut usize,
        offset: &mut u32,
    ) -> usize {
        let cigar_element_len = cigar_elements[ci].0;
        if self.skip_indel_next_to_intron(cigar_elements, ci) {
            *read_pos_exc += cigar_element_len as usize;
            return ci;
        }

        let mut desc_string = format!("-{}", cigar_element_len);
        let mut sequence_to_append = String::new();
        let quality_of_last_segment_before_del = if *read_pos_inc > 0 {
            query_qual[*read_pos_inc - 1]
        } else {
            33
        };
        let mut quality_of_segment = String::new();

        let mut multoffs = 0;
        let mut multoffp = 0;
        let mut nmoff = 0;

        let mut next_ci = ci;

        if self.is_insertion_or_deletion_with_next_matched(cigar_elements, ci) {
            let m_len = cigar_elements[ci + 1].0 as usize;
            let indel_len = cigar_elements[ci + 2].0 as usize;
            let begin = *read_pos_inc;

            let mut quality_of_segment_sb = String::new();
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
            quality_of_segment.push_str(&quality_of_segment_sb);

            let is_next_next_del = cigar_elements[ci + 2].1 == CigarKind::Deletion;
            multoffs += m_len + if is_next_next_del { indel_len } else { 0 };
            multoffp += m_len + if !is_next_next_del { indel_len } else { 0 };

            if self.is_next_after_num_matched(cigar_elements, ci, 3) {
                let mut vsn = 0;
                let tn = *read_pos_inc + multoffp;
                let ts = *start + multoffs as i32 + cigar_element_len as i32;
                let vext = self.config.extension_nuc as usize;
                let ci3_len = cigar_elements[ci + 3].0 as usize;

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
                    if (query_qual[tn_idx].saturating_sub(33) as f64) < self.config.base_qual {
                        break;
                    }
                    let ts_pos = ts + vi as i32;
                    if let Some(ref_ch) = self.reference.get_base(&self.region.chr, ts_pos as u32) {
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
                        sequence_to_append.push_str(&String::from_utf8_lossy(sub));
                    }
                    if let Some(sub) = query_qual.get(tn..tn + current_offset) {
                        quality_of_segment.push_str(&String::from_utf8_lossy(sub));
                    }
                }
            }
            next_ci += 2;
        } else if self.is_next_insertion(cigar_elements, ci) {
            let ins_len = cigar_elements[ci + 1].0 as usize;
            if let Some(sub) = query_seq.get(*read_pos_inc..*read_pos_inc + ins_len) {
                desc_string.push('^');
                desc_string.push_str(&String::from_utf8_lossy(sub));
            }
            if let Some(sub) = query_qual.get(*read_pos_inc..*read_pos_inc + ins_len) {
                quality_of_segment.push_str(&String::from_utf8_lossy(sub));
            }
            multoffp += ins_len;

            if self.is_next_after_num_matched(cigar_elements, ci, 2) {
                let m_len = cigar_elements[ci + 2].0 as usize;
                let mut vsn = 0;
                let tn = *read_pos_inc + multoffp;
                let ts = *start + cigar_element_len as i32;
                let vext = self.config.extension_nuc as usize;

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
                    if (query_qual[tn_idx].saturating_sub(33) as f64) < self.config.base_qual {
                        break;
                    }
                    let ts_pos = ts + vi as i32;
                    if let Some(ref_ch) = self.reference.get_base(&self.region.chr, ts_pos as u32) {
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
                        sequence_to_append.push_str(&String::from_utf8_lossy(sub));
                    }
                    if let Some(sub) = query_qual.get(tn..tn + current_offset) {
                        quality_of_segment.push_str(&String::from_utf8_lossy(sub));
                    }
                }
            }
            next_ci += 1;
        } else {
            if self.is_next_matched(cigar_elements, ci) {
                let m_len = cigar_elements[ci + 1].0 as usize;
                let mut vsn = 0;
                let vext = self.config.extension_nuc as usize;

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
                    if (query_qual[idx].saturating_sub(33) as f64) < self.config.base_qual {
                        break;
                    }
                    let ts_pos = *start + cigar_element_len as i32 + vi as i32;
                    if let Some(ref_ch) = self.reference.get_base(&self.region.chr, ts_pos as u32) {
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
                        sequence_to_append.push_str(&String::from_utf8_lossy(sub));
                    }
                    if let Some(sub) = query_qual.get(*read_pos_inc..*read_pos_inc + current_offset)
                    {
                        quality_of_segment.push_str(&String::from_utf8_lossy(sub));
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
            quality_of_segment.push(quality_of_last_segment_before_del as char);
        } else {
            let quality_of_segment_with_offset = query_qual[next_q_idx];
            let best_qual = if quality_of_last_segment_before_del > quality_of_segment_with_offset {
                quality_of_last_segment_before_del
            } else {
                quality_of_segment_with_offset
            };
            quality_of_segment.push(best_qual as char);
        }

        if *start >= self.region.start.0 as i32 && *start <= self.region.end.0 as i32 {
            self.add_variation_for_deletion(
                record.mapping_quality.unwrap_or(0) as i32,
                get_adjusted_nm(record),
                direction,
                read_len_match_ins,
                &desc_string,
                &quality_of_segment,
                nmoff,
                *start,
                *read_pos_exc,
                cigar_element_len,
            );
        }

        *start += cigar_element_len as i32 + *offset as i32 + multoffs as i32;
        *read_pos_inc += *offset as usize + multoffp;
        *read_pos_exc += *offset as usize + multoffp;

        next_ci
    }

    fn process_insertion(
        &mut self,
        record: &SamRecord,
        query_seq: &[u8],
        query_qual: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        ci: usize,
        direction: bool,
        read_len_match_ins: i32,
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

        let mut desc_string = if let Some(sub) =
            query_seq.get(*read_pos_inc..*read_pos_inc + cigar_element_len as usize)
        {
            String::from_utf8_lossy(sub).into_owned()
        } else {
            String::new()
        };

        let mut quality_string = if let Some(sub) =
            query_qual.get(*read_pos_inc..*read_pos_inc + cigar_element_len as usize)
        {
            String::from_utf8_lossy(sub).into_owned()
        } else {
            String::new()
        };

        let mut ss = String::new();
        let mut multoffs = 0;
        let mut multoffp = 0;
        let mut nmoff = 0;

        let mut next_ci = ci;

        if self.is_insertion_or_deletion_with_next_matched(cigar_elements, ci) {
            let m_len = cigar_elements[ci + 1].0 as usize;
            let indel_len = cigar_elements[ci + 2].0 as usize;
            let begin = *read_pos_inc + cigar_element_len as usize;

            let mut quality_segment_sb = String::new();
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
            quality_string.push_str(&quality_segment_sb);

            let is_next_next_del = cigar_elements[ci + 2].1 == CigarKind::Deletion;
            multoffs += m_len + if is_next_next_del { indel_len } else { 0 };
            multoffp += m_len + if !is_next_next_del { indel_len } else { 0 };

            let ci6 = if cigar_elements.len() > ci + 3 {
                cigar_elements[ci + 3].0 as usize
            } else {
                0
            };
            if ci6 != 0 && cigar_elements[ci + 3].1 == CigarKind::Match {
                let tpl = self.find_offset(
                    *start + multoffs as i32,
                    *read_pos_inc + cigar_element_len as usize + multoffp,
                    ci6,
                    query_seq,
                    query_qual,
                );
                *offset = tpl.offset as u32;
                ss.push_str(&tpl.sequence);
                quality_string.push_str(&tpl.quality_sequence);
                for osi in 0..tpl.offset {
                    inc_cnt(&mut self.ref_coverage, *start + multoffs as i32 + osi, 1);
                }
            }
            next_ci += 2;
        } else {
            if self.is_next_matched(cigar_elements, ci) {
                let mut vsn = 0;
                let m_len = cigar_elements[ci + 1].0 as usize;
                let vext = self.config.extension_nuc as usize;

                let mut current_offset = 0;
                for vi in 0..m_len {
                    if vsn > vext {
                        break;
                    }
                    let idx = *read_pos_inc + cigar_element_len as usize + vi;
                    if idx >= query_seq.len() {
                        break;
                    }
                    let q_char = query_seq[idx];
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    if (query_qual[idx].saturating_sub(33) as f64) < self.config.base_qual {
                        break;
                    }
                    let ts_pos = *start + vi as i32;
                    if let Some(ref_ch) = self.reference.get_base(&self.region.chr, ts_pos as u32) {
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
                    if let Some(sub) = query_seq.get(
                        *read_pos_inc + cigar_element_len as usize
                            ..*read_pos_inc + cigar_element_len as usize + current_offset,
                    ) {
                        ss.push_str(&String::from_utf8_lossy(sub));
                    }
                    if let Some(sub) = query_qual.get(
                        *read_pos_inc + cigar_element_len as usize
                            ..*read_pos_inc + cigar_element_len as usize + current_offset,
                    ) {
                        quality_string.push_str(&String::from_utf8_lossy(sub));
                    }
                    for osi in 0..current_offset {
                        inc_cnt(&mut self.ref_coverage, *start + osi as i32, 1);
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
                let tpl = crate::realign::adj_ins_pos(
                    *start - 1,
                    desc_string.clone(),
                    self.reference,
                    &self.region.chr,
                );
                let adjusted_insertion_position = tpl.base_insert;
                let adjusted_desc_string = tpl.insertion_sequence;

                if *read_pos_inc as i32 - 1 - (*start - 1 - adjusted_insertion_position) > 0 {
                    insertion_position = adjusted_insertion_position;
                    desc_string = adjusted_desc_string;
                }
            }

            let ins_key = format!("+{}", desc_string);
            increment(
                &mut self.position_to_insertion_count,
                insertion_position,
                &ins_key,
            );

            let hv = get_variation(&mut self.insertion_variants, insertion_position, &ins_key);
            hv.inc_dir(direction);
            hv.cnt += 1;

            let tp = if *read_pos_exc < read_len_match_ins as usize - *read_pos_exc {
                *read_pos_exc as i32 + 1
            } else {
                read_len_match_ins - *read_pos_exc as i32
            };

            let mut tmpq = 0.0;
            if !quality_string.is_empty() {
                for b in quality_string.bytes() {
                    tmpq += (b as f64) - 33.0;
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
            hv.nm += (get_adjusted_nm(record) - nmoff) as f64;

            let index_of_insertion_in_query_sequence =
                *read_pos_inc as i32 - 1 - (*start - 1 - insertion_position);
            if insertion_position > position {
                if let Some(sub_base) = query_seq
                    .get(index_of_insertion_in_query_sequence as usize)
                    .copied()
                {
                    if let Some(ref_ch) = self
                        .reference
                        .get_base(&self.region.chr, insertion_position as u32)
                    {
                        if sub_base == ref_ch {
                            let base_qual_val =
                                query_qual[index_of_insertion_in_query_sequence as usize]
                                    .saturating_sub(33) as f64;
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
                                    (get_adjusted_nm(record) - nmoff) as f64,
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
                if let Some(ref_base_char) = self
                    .reference
                    .get_base(&self.region.chr, insertion_position as u32)
                {
                    let ref_base_str = (ref_base_char as char).to_string();
                    let ttref = get_variation(
                        &mut self.non_insertion_variants,
                        insertion_position,
                        &ref_base_str,
                    );
                    ttref.inc_dir(direction);
                    ttref.cnt += 1;
                    ttref.pstd = hv.pstd;
                    ttref.qstd = hv.qstd;
                    ttref.pmean += tp as f64;
                    ttref.qmean += tmpq;
                    ttref.mqmean += record.mapping_quality.unwrap_or(0) as f64;
                    ttref.pp = tp;
                    ttref.pq = tmpq;
                    ttref.nm += (get_adjusted_nm(record) - nmoff) as f64;
                    inc_cnt(&mut self.ref_coverage, insertion_position, 1);
                }
            }
        }

        *read_pos_inc += cigar_element_len as usize + *offset as usize + multoffp;
        *read_pos_exc += cigar_element_len as usize + *offset as usize + multoffp;
        *start += *offset as i32 + multoffs as i32;

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
        let parts: Vec<&str> = sa_tag_string.split(',').collect();
        if parts.len() < 4 {
            return false;
        }
        let sa_chromosome = parts[0];
        let sa_position = parts[1].parse::<i32>().unwrap_or(0);
        let sa_direction_string = parts[2];
        let sa_cigar = parts[3];
        let sa_direction_is_forward = sa_direction_string == "+";

        let sa_5clip = {
            static RE: std::sync::OnceLock<regex::Regex> = std::sync::OnceLock::new();
            RE.get_or_init(|| regex::Regex::new(r"^\d\d+S").unwrap())
        };
        let sa_3clip = {
            static RE: std::sync::OnceLock<regex::Regex> = std::sync::OnceLock::new();
            RE.get_or_init(|| regex::Regex::new(r"\d\d+S$").unwrap())
        };

        let is_chimeric = if is_5_side {
            sa_5clip.is_match(sa_cigar)
        } else {
            sa_3clip.is_match(sa_cigar)
        };

        let sa_chr_matches = if let Some(ref_name) = &record.reference_name {
            sa_chromosome == ref_name
        } else {
            false
        };

        let is_chimeric_flag = ((dir && sa_direction_is_forward)
            || (!dir && !sa_direction_is_forward))
            && sa_chr_matches
            && ((sa_position - position).abs() < 2 * self.max_read_length)
            && is_chimeric;

        if self.config.debug_print.is_some() && is_chimeric_flag {
            let record_name = String::from_utf8_lossy(&record.name);
            let ref_name = record.reference_name.as_deref().unwrap_or("*");
            log::debug!(
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
                    query_quality[si as usize].saturating_sub(33) as f64,
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
                    query_quality[idx].saturating_sub(33) as f64,
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

    fn process_soft_clip(
        &mut self,
        record: &SamRecord,
        query_seq: &[u8],
        query_qual: &[u8],
        cigar_elements: &[(u32, CigarKind)],
        direction: bool,
        _read_len_match_ins: i32,
        total_length_including_soft_clipped: i32,
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
                    let sequence =
                        get_reverse_complemented_sequence(record, 0, cigar_element_len as usize);
                    if sequence.len() >= 17 {
                        let reverse_complemented_seed = &sequence[..17];
                        if let Some(positions) = self.reference.seed.get(reverse_complemented_seed)
                        {
                            if positions.len() == 1
                                && (*start - positions[0] as i32).abs() < 2 * self.max_read_length
                            {
                                *read_pos_inc += cigar_element_len as usize;
                                *offset = 0;
                                *start = position;
                                if self.config.debug_print.is_some() {
                                    log::debug!(
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
                let Some(ref_base) = self
                    .reference
                    .get_base(&self.region.chr, (*start - 1) as u32)
                else {
                    break;
                };
                let q_base = query_seq[element_len as usize - 1];
                let base_qual = query_qual[element_len as usize - 1].saturating_sub(33);
                if q_base == ref_base && base_qual > 10 {
                    let ref_base_str = (ref_base as char).to_string();
                    let variation =
                        get_variation(&mut self.non_insertion_variants, *start - 1, &ref_base_str);
                    variation.add_cnt(
                        direction,
                        element_len as i32,
                        base_qual as f64,
                        record.mapping_quality.unwrap_or(0) as i32,
                        get_adjusted_nm(record) as f64,
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
                let mut si = element_len as i32 - 1;
                while si >= 0 {
                    let q_char = query_seq[si as usize];
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    let base_qual = query_qual[si as usize].saturating_sub(33);
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
                    get_adjusted_nm(record),
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
                    let sequence = get_reverse_complemented_sequence(
                        record,
                        -(cigar_element_len as isize),
                        cigar_element_len as usize,
                    );
                    if sequence.len() >= 17 {
                        let reverse_complemented_seed = &sequence[sequence.len() - 17..];
                        if let Some(positions) = self.reference.seed.get(reverse_complemented_seed)
                        {
                            if positions.len() == 1
                                && (*start - positions[0] as i32).abs() < 2 * self.max_read_length
                            {
                                *read_pos_inc += cigar_element_len as usize;
                                *offset = 0;
                                *start = position;
                                if self.config.debug_print.is_some() {
                                    log::debug!(
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

            let mut element_len = cigar_element_len;
            while *read_pos_inc < query_seq.len() {
                // §9: single fetch — the loop condition only needs is_some.
                let Some(ref_base) = self.reference.get_base(&self.region.chr, *start as u32)
                else {
                    break;
                };
                let q_base = query_seq[*read_pos_inc];
                let base_qual = query_qual[*read_pos_inc].saturating_sub(33);
                if q_base == ref_base && base_qual > 10 {
                    let ref_base_str = (ref_base as char).to_string();
                    let variation =
                        get_variation(&mut self.non_insertion_variants, *start, &ref_base_str);
                    variation.add_cnt(
                        direction,
                        total_length_including_soft_clipped - *read_pos_exc as i32,
                        base_qual as f64,
                        record.mapping_quality.unwrap_or(0) as i32,
                        get_adjusted_nm(record) as f64,
                        self.config.base_qual,
                    );
                    inc_cnt(&mut self.ref_coverage, *start, 1);
                    *read_pos_inc += 1;
                    *start += 1;
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
                for si in 0..element_len as usize {
                    let idx = *read_pos_inc + si;
                    if idx >= query_seq.len() {
                        break;
                    }
                    let q_char = query_seq[idx];
                    if q_char == b'N' || q_char == b'n' {
                        break;
                    }
                    let base_qual = query_qual[idx].saturating_sub(33);
                    if base_qual <= 12 {
                        number_of_low_quality_bases += 1;
                    }
                    if number_of_low_quality_bases > 1 {
                        break;
                    }
                    sum_of_read_qualities += base_qual as i32;
                    number_of_high_quality_bases += 1;
                }
                self.sclip3_high_quality_processing(
                    query_seq,
                    record.mapping_quality.unwrap_or(0) as i32,
                    query_qual,
                    get_adjusted_nm(record),
                    direction,
                    sum_of_read_qualities,
                    number_of_high_quality_bases,
                    number_of_low_quality_bases,
                    *start,
                    element_len,
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
    ) -> Offset {
        let mut offset = 0;
        let mut tnm = 0;
        let mut vsn = 0;
        let vext = self.config.extension_nuc as usize;

        for vi in 0..cigar_length {
            let read_idx = read_position + vi;
            if read_idx >= query_sequence.len() {
                break;
            }
            if query_sequence[read_idx] == b'N' || query_sequence[read_idx] == b'n' {
                break;
            }
            let q = query_quality[read_idx].saturating_sub(33) as f64;
            if q < self.config.base_qual {
                break;
            }
            if let Some(ref_ch) = self
                .reference
                .get_base(&self.region.chr, (reference_position + vi as i32) as u32)
            {
                let ch = query_sequence[read_idx];
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

        let mut ss = String::new();
        let mut q_seq = String::new();
        if offset > 0 {
            if let Some(sub) = query_sequence.get(read_position..read_position + offset) {
                ss = String::from_utf8_lossy(sub).into_owned();
            }
            if let Some(sub) = query_quality.get(read_position..read_position + offset) {
                q_seq = String::from_utf8_lossy(sub).into_owned();
            }
        }

        Offset {
            offset: offset as i32,
            sequence: ss,
            quality_sequence: q_seq,
            offset_number_of_mismatches: tnm,
        }
    }

    /// Finalize and return the collected VariationData.
    /// Ports `CigarParser.java:101-150`.
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
        let mend = mate_start + total_length_including_soft_clipped;

        let mut end = start;
        for &(len, op) in cigar_elements {
            if matches!(
                op,
                CigarKind::Match
                    | CigarKind::Deletion
                    | CigarKind::Skip
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
            ) {
                end += len as i32;
            }
        }

        let mut soft5 = 0;
        if let Some(&(len, op)) = cigar_elements.first() {
            if matches!(op, CigarKind::SoftClip | CigarKind::HardClip) {
                let len = len as usize;
                if len > 0 && len <= query_qual.len() {
                    let qual_val = query_qual[len - 1] as f64 - 33.0;
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
                    let qual_val = query_qual[query_qual.len() - len] as f64 - 33.0;
                    if qual_val > self.config.base_qual {
                        soft3 = end;
                    }
                }
            }
        }

        let read_dir_num = if read_direction { -1 } else { 1 };
        let mate_dir_num = if mate_direction { 1 } else { -1 };
        let min_d = 75;

        let mate_ref_name = get_mate_reference_name(record);
        if mate_ref_name == "=" {
            let mut mlen = record.template_length;

            if let Some(ref mc) = record.mc {
                if MC_CLIP_CLIP.is_match(mc) {
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
            let limit = (min_svc_dist * self.max_read_length as f64) as i32;

            if read_dir_num * mate_dir_num == -1
                && (mlen * read_dir_num) > 0
                && query_qual.len() > 15
            {
                // Deletion candidate
                mlen = if mate_start > start {
                    mend - start
                } else {
                    end - mate_start
                };
                if mlen.abs()
                    > self.config.sv_read_len as i32
                        + (self.config.sv_min_af as f64 * self.config.sv_flanking as f64) as i32
                {
                    if read_dir_num == 1 {
                        if self.sv_structures.svfdel.is_empty()
                            || start - self.sv_structures.svdelfend > limit
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
                                mlen,
                                soft3,
                                self.max_read_length as f64 / 2.0,
                                query_qual[15] as f64 - 33.0,
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                        }
                        self.sv_structures.svdelfend = end;
                    } else {
                        if self.sv_structures.svrdel.is_empty()
                            || start - self.sv_structures.svdelrend > limit
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
                                mlen,
                                soft5,
                                self.max_read_length as f64 / 2.0,
                                query_qual[15] as f64 - 33.0,
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                        }
                        self.sv_structures.svdelrend = end;
                    }

                    if !self.sv_structures.svfdel.is_empty()
                        && (start - self.sv_structures.svdelfend).abs() <= limit
                    {
                        if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdel.is_empty()
                        && (start - self.sv_structures.svdelrend).abs() <= limit
                    {
                        if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfdup.is_empty()
                        && (start - self.sv_structures.svdupfend).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdup.is_empty()
                        && (start - self.sv_structures.svduprend).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfinv5.is_empty()
                        && (start - self.sv_structures.svinvfend5).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfinv5.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrinv5.is_empty()
                        && (start - self.sv_structures.svinvrend5).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrinv5.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfinv3.is_empty()
                        && (start - self.sv_structures.svinvfend3).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfinv3.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrinv3.is_empty()
                        && (start - self.sv_structures.svinvrend3).abs() <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrinv3.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                }
            } else if read_dir_num * mate_dir_num == -1
                && read_dir_num * mlen < 0
                && query_qual.len() > 15
            {
                // Duplication
                if read_dir_num == 1 {
                    if self.sv_structures.svfdup.is_empty()
                        || start - self.sv_structures.svdupfend > limit
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
                            mlen,
                            soft3,
                            self.max_read_length as f64 / 2.0,
                            query_qual[15] as f64 - 33.0,
                            mq as f64,
                            number_of_mismatches as f64,
                            self.config.base_qual,
                        );
                    }
                    self.sv_structures.svdupfend = end;
                } else {
                    if self.sv_structures.svrdup.is_empty()
                        || start - self.sv_structures.svduprend > limit
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
                            mlen,
                            soft5,
                            self.max_read_length as f64 / 2.0,
                            query_qual[15] as f64 - 33.0,
                            mq as f64,
                            number_of_mismatches as f64,
                            self.config.base_qual,
                        );
                    }
                    self.sv_structures.svduprend = end;
                }
                if !self.sv_structures.svfdup.is_empty()
                    && (start - self.sv_structures.svdupfend).abs() <= limit
                {
                    if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrdup.is_empty()
                    && (start - self.sv_structures.svduprend).abs() <= limit
                {
                    if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svfdel.is_empty()
                    && (start - self.sv_structures.svdelfend).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrdel.is_empty()
                    && (start - self.sv_structures.svdelrend).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svfinv5.is_empty()
                    && (start - self.sv_structures.svinvfend5).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svfinv5.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrinv5.is_empty()
                    && (start - self.sv_structures.svinvrend5).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svrinv5.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svfinv3.is_empty()
                    && (start - self.sv_structures.svinvfend3).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svfinv3.last_mut() {
                        sdref.disc += 1;
                    }
                }
                if !self.sv_structures.svrinv3.is_empty()
                    && (start - self.sv_structures.svinvrend3).abs() <= min_d
                {
                    if let Some(sdref) = self.sv_structures.svrinv3.last_mut() {
                        sdref.disc += 1;
                    }
                }
            } else if read_dir_num * mate_dir_num == 1 && query_qual.len() > 15 {
                // Inversion
                if read_dir_num == 1 && mlen != 0 {
                    if mlen < -3 * self.max_read_length {
                        if self.sv_structures.svfinv3.is_empty()
                            || start - self.sv_structures.svinvfend3 > limit
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
                                mlen,
                                soft3,
                                self.max_read_length as f64 / 2.0,
                                query_qual[15] as f64 - 33.0,
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                        self.sv_structures.svinvfend3 = end;
                    } else if mlen > 3 * self.max_read_length {
                        if self.sv_structures.svfinv5.is_empty()
                            || start - self.sv_structures.svinvfend5 > limit
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
                                mlen,
                                soft3,
                                self.max_read_length as f64 / 2.0,
                                query_qual[15] as f64 - 33.0,
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                        self.sv_structures.svinvfend5 = end;
                    }
                } else if mlen != 0 {
                    if mlen < -3 * self.max_read_length {
                        if self.sv_structures.svrinv3.is_empty()
                            || start - self.sv_structures.svinvrend3 > limit
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
                                mlen,
                                soft5,
                                self.max_read_length as f64 / 2.0,
                                query_qual[15] as f64 - 33.0,
                                mq as f64,
                                number_of_mismatches as f64,
                                self.config.base_qual,
                            );
                            sdref.disc += 1;
                        }
                        self.sv_structures.svinvrend3 = end;
                    } else if mlen > 3 * self.max_read_length {
                        if self.sv_structures.svrinv5.is_empty()
                            || start - self.sv_structures.svinvrend5 > limit
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
                                mlen,
                                soft5,
                                self.max_read_length as f64 / 2.0,
                                query_qual[15] as f64 - 33.0,
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
                        && start - self.sv_structures.svdelfend <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdel.is_empty()
                        && start - self.sv_structures.svdelrend <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svfdup.is_empty()
                        && start - self.sv_structures.svdupfend <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                    if !self.sv_structures.svrdup.is_empty()
                        && start - self.sv_structures.svduprend <= min_d
                    {
                        if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                            sdref.disc += 1;
                        }
                    }
                }
            }
        } else if query_qual.len() > 15 {
            // Inter-chr translocation
            let mchr = mate_ref_name;
            // ports CigarParser.java:2208-2215: the MC-tag and mate-MQ(MQ tag) filters,
            // not the read's own mapping quality.
            let sv_allowed = !(record
                .mc
                .as_ref()
                .is_some_and(|mc| MC_CLIP_CLIP.is_match(mc))
                || record.mq.is_some_and(|mate_mq| mate_mq < 15));
            let mq = record.mapping_quality.unwrap_or(0);

            let min_svc_dist = 1.5;
            let limit = (min_svc_dist * self.max_read_length as f64) as i32;

            if sv_allowed {
                if read_dir_num == 1 {
                    let last_fus_end = self
                        .sv_structures
                        .svfusfend
                        .get(&mchr)
                        .copied()
                        .unwrap_or(0);
                    if !self.sv_structures.svffus.contains_key(&mchr)
                        || start - last_fus_end > limit
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
                                query_qual[15] as f64 - 33.0,
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
                        || start - last_fus_end > limit
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
                                query_qual[15] as f64 - 33.0,
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

            let min_sv_pos = 25;
            if !self.sv_structures.svfdel.is_empty()
                && start - self.sv_structures.svdelfend <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfdel.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrdel.is_empty()
                && start - self.sv_structures.svdelrend <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrdel.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svfdup.is_empty()
                && start - self.sv_structures.svdupfend <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfdup.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrdup.is_empty()
                && start - self.sv_structures.svduprend <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrdup.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svfinv5.is_empty()
                && start - self.sv_structures.svinvfend5 <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfinv5.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrinv5.is_empty()
                && start - self.sv_structures.svinvrend5 <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrinv5.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svfinv3.is_empty()
                && start - self.sv_structures.svinvfend3 <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svfinv3.last_mut() {
                    sdref.disc += 1;
                }
            }
            if !self.sv_structures.svrinv3.is_empty()
                && start - self.sv_structures.svinvrend3 <= min_sv_pos
            {
                if let Some(sdref) = self.sv_structures.svrinv3.last_mut() {
                    sdref.disc += 1;
                }
            }
        }
    }

    fn skip_overlapping_reads(
        &self,
        record: &SamRecord,
        position: i32,
        dir: bool,
        mate_alignment_start: i32,
        start: i32,
    ) -> bool {
        if self.config.unique_mode_alignment_enabled
            && is_paired_and_same_chromosome(record)
            && !dir
            && start >= mate_alignment_start
        {
            return true;
        }
        if self.config.unique_mode_second_in_pair_enabled
            && (record.flags & 0x80) != 0
            && is_paired_and_same_chromosome(record)
            && self.is_reads_overlap(record, position, mate_alignment_start, start)
        {
            return true;
        }
        false
    }

    fn is_reads_overlap(
        &self,
        record: &SamRecord,
        position: i32,
        mate_alignment_start: i32,
        start: i32,
    ) -> bool {
        if position >= mate_alignment_start {
            let ref_len = get_reference_length(&record.cigar) as i32;
            start >= mate_alignment_start && start <= (mate_alignment_start + ref_len - 1)
        } else {
            let alignment_end = record.alignment_end.map(|p| p.0 as i32).unwrap_or(0);
            start >= mate_alignment_start && mate_alignment_start <= alignment_end
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
        sdref.mend = mate_start_ms + rlen;
    }
    if softp != 0 {
        if dir == 1 {
            if (softp - sdref.end).abs() < 10 {
                let soft_count = sdref.soft.entry(softp).or_insert(0);
                *soft_count += 1;
            }
        } else {
            if (softp - sdref.start).abs() < 10 {
                let soft_count = sdref.soft.entry(softp).or_insert(0);
                *soft_count += 1;
            }
        }
    }
}

fn is_paired_and_same_chromosome(record: &SamRecord) -> bool {
    (record.flags & 0x1) != 0 && get_mate_reference_name(record) == "="
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
        .map(|&(len, _)| len)
        .sum()
}

/// Helper function to calculate the length of insertion and deletion CIGAR elements.
/// Ports `CigarParser.java:1636-1644`.
pub fn get_insertion_deletion_length(cigar: &[(u32, CigarKind)]) -> u32 {
    cigar
        .iter()
        .filter(|(_, kind)| matches!(kind, CigarKind::Insertion | CigarKind::Deletion))
        .map(|(len, _)| *len)
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
                CigarKind::Match
                    | CigarKind::Insertion
                    | CigarKind::SoftClip
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
            )
        })
        .map(|(len, _)| *len)
        .sum()
}

/// Parse the `-a` amplicon parameter "distance:fraction" into (distance, fraction).
/// Falls back to (10, 0.95) on parse failure, matching Java parseCigarWithAmpCase.
fn parse_amp_param(param: &str) -> (i32, f64) {
    let parts: Vec<&str> = param.split(':').collect();
    let distance = parts.first().and_then(|s| s.parse::<i32>().ok());
    let fraction = parts.get(1).and_then(|s| s.parse::<f64>().ok());
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
        .map(|(len, _)| *len as i32)
        .sum()
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
    let read_length = get_aligned_length_amp(cigar_elements);
    let segstart = record.alignment_start.map(|p| p.0 as i32).unwrap_or(0);
    let mut segstart = segstart;
    let mut segend = segstart + read_length - 1;
    let region_start = region.start.0 as i32;
    let region_end = region.end.0 as i32;

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
        if is_mate_ref_name_equal && record.template_length != 0 {
            if record.template_length > 0 {
                segend = segstart + record.template_length - 1;
            } else {
                segstart = record.mate_alignment_start.map(|p| p.0 as i32).unwrap_or(0);
                segend = segstart - record.template_length - 1;
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
        .filter(|(_, kind)| {
            matches!(
                kind,
                CigarKind::Match
                    | CigarKind::Insertion
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
            )
        })
        .map(|(len, _)| *len)
        .sum()
}

/// Ports `CigarParser.java:1530-1573`.
pub fn cleanup_cigar(elements: &mut Vec<(u32, CigarKind)>) {
    // 1. Process leading elements
    let mut i = 0;
    while i < elements.len() {
        let (len, kind) = elements[i];
        if consumes_both(kind) {
            break;
        }
        if kind == CigarKind::Insertion {
            elements[i] = (len, CigarKind::SoftClip);
            i += 1;
        } else if kind == CigarKind::HardClip {
            elements.remove(i);
        } else {
            i += 1;
        }
    }

    // 2. Process trailing elements
    let mut i = elements.len();
    while i > 0 {
        let idx = i - 1;
        let (len, kind) = elements[idx];
        if consumes_both(kind) {
            break;
        }
        if kind == CigarKind::Insertion {
            elements[idx] = (len, CigarKind::SoftClip);
            i = idx;
        } else if kind == CigarKind::HardClip {
            elements.remove(idx);
            i = idx;
        } else {
            i = idx;
        }
    }
}

fn consumes_both(kind: CigarKind) -> bool {
    matches!(
        kind,
        CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch
    )
}

fn is_not_equals_n(ref_base: Option<u8>) -> bool {
    if let Some(b) = ref_base {
        b != b'N' && b != b'n'
    } else {
        false
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

fn increment(hash: &mut BTreeMap<i32, HashMap<String, i32>>, k1: i32, k2: &str) {
    let map = hash.entry(k1).or_default();
    *map.entry(k2.to_string()).or_default() += 1;
}

fn inc_cnt(hash: &mut BTreeMap<i32, i32>, start: i32, count: i32) {
    *hash.entry(start).or_default() += count;
}

fn get_variation_from_seq(sclip: &mut Sclip, idx: i32, ch: u8) -> &mut Variation {
    let map = sclip.seq.entry(idx).or_default();
    map.entry(ch).or_insert_with(Variation::new)
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
) -> Result<VariationData, VardictError> {
    let mut parser = CigarParser::new(
        reference,
        region,
        config,
        svflag,
        total_reads,
        duplicate_reads,
    );
    for record in records {
        parser.parse_cigar(record)?;
    }
    Ok(parser.into_variation_data())
}

/// Run pileup on a slice of records, updating the existing `VariationData` in-place.
pub fn pileup_into(
    var_data: &mut VariationData,
    records: &[SamRecord],
    reference: &Reference,
    region: &Region,
    config: &Config,
    svflag: bool,
    total_reads: i32,
    duplicate_reads: i32,
) -> Result<(), VardictError> {
    let mut parser = CigarParser {
        reference,
        region,
        config,
        svflag,
        non_insertion_variants: std::mem::take(&mut var_data.non_insertion_variants),
        insertion_variants: std::mem::take(&mut var_data.insertion_variants),
        position_to_insertion_count: std::mem::take(&mut var_data.position_to_insertion_count),
        position_to_deletion_count: std::mem::take(&mut var_data.position_to_deletion_count),
        sv_structures: std::mem::take(&mut var_data.sv_structures),
        ref_coverage: std::mem::take(&mut var_data.ref_coverage),
        soft_clips_5_end: std::mem::take(&mut var_data.soft_clips_5_end),
        soft_clips_3_end: std::mem::take(&mut var_data.soft_clips_3_end),
        max_read_length: var_data.max_read_length,
        splice: std::mem::take(&mut var_data.splice),
        mnp: std::mem::take(&mut var_data.mnp),
        splice_count: std::mem::take(&mut var_data.splice_count),
        discordant_count: 0,
        total_reads,
        duplicate_reads,
    };

    for record in records {
        parser.parse_cigar(record)?;
    }

    let updated = parser.into_variation_data();

    var_data.non_insertion_variants = updated.non_insertion_variants;
    var_data.insertion_variants = updated.insertion_variants;
    var_data.position_to_insertion_count = updated.position_to_insertion_count;
    var_data.position_to_deletion_count = updated.position_to_deletion_count;
    var_data.sv_structures = updated.sv_structures;
    var_data.ref_coverage = updated.ref_coverage;
    var_data.soft_clips_5_end = updated.soft_clips_5_end;
    var_data.soft_clips_3_end = updated.soft_clips_3_end;
    var_data.max_read_length = var_data.max_read_length.max(updated.max_read_length);
    var_data.splice = updated.splice;
    var_data.mnp = updated.mnp;
    var_data.splice_count = updated.splice_count;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Pos1;

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

        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
    fn cleanup_cigar_test() {
        let mut cigar1 = vec![
            (5, CigarKind::HardClip),
            (3, CigarKind::Insertion),
            (10, CigarKind::Match),
            (4, CigarKind::HardClip),
        ];
        cleanup_cigar(&mut cigar1);
        assert_eq!(
            cigar1,
            vec![(3, CigarKind::SoftClip), (10, CigarKind::Match)]
        );

        let mut cigar2 = vec![
            (10, CigarKind::Match),
            (2, CigarKind::Insertion),
            (5, CigarKind::HardClip),
        ];
        cleanup_cigar(&mut cigar2);
        assert_eq!(
            cigar2,
            vec![(10, CigarKind::Match), (2, CigarKind::SoftClip)]
        );
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

        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);

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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);
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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);
        parser.max_read_length = 150;
        parser.parse_cigar(&mk_record(None, None)).unwrap();
        let fus = parser
            .sv_structures
            .svffus
            .get("chr2")
            .expect("fusion list");
        assert_eq!(fus.last().unwrap().variation.cnt, 1);

        // Mate MQ tag < 15: Java skips collection.
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);
        parser.max_read_length = 150;
        parser.parse_cigar(&mk_record(None, Some(10u8))).unwrap();
        assert!(!parser.sv_structures.svffus.contains_key("chr2"));

        // MC tag with clips on both ends: Java skips collection.
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);
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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);
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
        let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0);
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
}
