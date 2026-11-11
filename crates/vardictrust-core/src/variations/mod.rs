//! Variation data structures — intermediate variant accumulators and SV helpers.
//!
//! These port the core data types from VarDictJava's `variations/` package:
//! - `Variation`: intermediate variant accumulator (sums, not means)
//! - `Sclip`: soft-clip data, extends Variation
//! - `Mate`: SV mate data
//! - `Cluster`: SV cluster data, extends Mate
//! - `Variant`: finalized output variant (populated in Phase 3)
//! - `Vars`: per-position variant container

use indexmap::map::{raw_entry_v1::RawEntryMut, RawEntryApiV1};
use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet};

pub type FxBuildHasher = std::hash::BuildHasherDefault<rustc_hash::FxHasher>;

pub mod util;

/// Intermediate variant accumulator — stores running sums (NOT means) of
/// position, quality, and mapping quality. The sums are divided by count
/// later (in ToVarsBuilder) to produce means.
///
/// Ports `variations/Variation.java`.
#[derive(Debug, Clone)]
pub struct Variation {
    /// Variant count.
    pub cnt: i32,
    /// Forward strand count.
    pub fwd: i32,
    /// Reverse strand count.
    pub rev: i32,
    /// Sum of variant positions in read (not yet a mean).
    pub pmean: f64,
    /// Sum of base qualities (not yet a mean).
    pub qmean: f64,
    /// Sum of mapping qualities (not yet a mean).
    pub mqmean: f64,
    /// Sum of number-of-mismatches (not yet a mean).
    pub nm: f64,
    /// Number of high-quality reads with this variant.
    pub hicnt: i32,
    /// Number of low-quality reads with this variant.
    pub locnt: i32,
    /// True when variant found in at least 2 different read positions.
    pub pstd: bool,
    /// True when variant read with at least 2 different base qualities.
    pub qstd: bool,
    /// Position in read for previous instance (used to detect pstd).
    pub pp: i32,
    /// Base quality for previous instance (used to detect qstd).
    pub pq: f64,
    /// Adjusted count for indels due to local realignment.
    pub extracnt: i32,
}

impl Variation {
    /// Create a new zeroed Variation.
    pub fn new() -> Self {
        Self {
            cnt: 0,
            fwd: 0,
            rev: 0,
            pmean: 0.0,
            qmean: 0.0,
            mqmean: 0.0,
            nm: 0.0,
            hicnt: 0,
            locnt: 0,
            pstd: false,
            qstd: false,
            pp: 0,
            pq: 0.0,
            extracnt: 0,
        }
    }

    /// Increment direction counter.
    /// `is_reverse`: true → increment reverse, false → increment forward.
    /// Ports `Variation.incDir(boolean)`.
    pub fn inc_dir(&mut self, is_reverse: bool) {
        if is_reverse {
            self.rev += 1;
        } else {
            self.fwd += 1;
        }
    }

    /// Decrement direction counter.
    /// Ports `Variation.decDir(boolean)`.
    pub fn dec_dir(&mut self, is_reverse: bool) {
        if is_reverse {
            self.rev -= 1;
        } else {
            self.fwd -= 1;
        }
    }

    /// Get direction counter.
    /// Ports `Variation.getDir(boolean)`.
    pub fn get_dir(&self, is_reverse: bool) -> i32 {
        if is_reverse {
            self.rev
        } else {
            self.fwd
        }
    }

    /// Add counts to direction counters.
    /// Ports `Variation.addDir(boolean, int)`.
    pub fn add_dir(&mut self, is_reverse: bool, val: i32) {
        if is_reverse {
            self.rev += val;
        } else {
            self.fwd += val;
        }
    }

    /// Subtract counts from direction counters.
    /// Ports `Variation.subDir(boolean, int)`.
    pub fn sub_dir(&mut self, is_reverse: bool, val: i32) {
        if is_reverse {
            self.rev -= val;
        } else {
            self.fwd -= val;
        }
    }

    /// Add a read's contribution to this variation.
    /// Ports `CigarParser.addCnt` (lines 1407-1424).
    pub fn add_cnt(
        &mut self,
        is_reverse: bool,
        read_pos: i32,
        base_qual: f64,
        map_qual: i32,
        nm: f64,
        goodq: f64,
    ) {
        self.cnt += 1;
        self.inc_dir(is_reverse);
        self.pmean += f64::from(read_pos);
        self.qmean += base_qual;
        self.mqmean += f64::from(map_qual);
        self.nm += nm;
        if base_qual >= goodq {
            self.hicnt += 1;
        } else {
            self.locnt += 1;
        }
    }

    /// Subtract a read's contribution from this variation.
    /// Ports `CigarParser.subCnt` (lines 1378-1396).
    pub fn sub_cnt(
        &mut self,
        is_reverse: bool,
        read_pos: i32,
        base_qual: f64,
        map_qual: i32,
        nm: f64,
        goodq: f64,
    ) {
        self.cnt -= 1;
        self.dec_dir(is_reverse);
        self.pmean -= f64::from(read_pos);
        self.qmean -= base_qual;
        self.mqmean -= f64::from(map_qual);
        self.nm -= nm;
        if base_qual >= goodq {
            self.hicnt -= 1;
        } else {
            self.locnt -= 1;
        }
    }
}

impl Default for Variation {
    fn default() -> Self {
        Self::new()
    }
}

/// Soft-clip data structure. Extends `Variation` with per-position base counts
/// and SV-related fields.
///
/// Ports `variations/Sclip.java`.
#[derive(Debug, Clone)]
pub struct Sclip {
    /// Base variation data (inherited from Variation).
    pub variation: Variation,

    /// Position → (base → count) for high-quality bases in the soft-clipped region.
    /// Uses `BTreeMap` to match Java's `TreeMap` ordering.
    pub nt: BTreeMap<i32, BTreeMap<u8, i32>>,

    /// Position → (base → Variation) for base variations in the soft-clipped region.
    pub seq: BTreeMap<i32, BTreeMap<u8, Variation>>,

    /// Consensus sequence in soft-clipped reads.
    pub sequence: Option<Vec<u8>>,

    /// Whether this sclip has been consumed/used.
    pub used: bool,

    // -- SV fields --
    /// SV: start position.
    pub start: i32,
    /// SV: end position.
    pub end: i32,
    /// SV: mate start.
    pub mstart: i32,
    /// SV: mate end.
    pub mend: i32,
    /// SV: mate length.
    pub mlen: i32,
    /// SV: discordant count.
    pub disc: i32,
    /// SV: soft-clip position.
    pub softp: i32,
    /// SV: softclip positions → counts (insertion-ordered).
    pub soft: IndexMap<i32, i32>,
    /// SV: list of mate records.
    pub mates: Vec<Mate>,
}

impl Sclip {
    /// Create a new Sclip with default fields.
    pub fn new() -> Self {
        Self {
            variation: Variation::new(),
            nt: BTreeMap::new(),
            seq: BTreeMap::new(),
            sequence: None,
            used: false,
            start: 0,
            end: 0,
            mstart: 0,
            mend: 0,
            mlen: 0,
            disc: 0,
            softp: 0,
            soft: IndexMap::new(),
            mates: Vec::new(),
        }
    }
}

impl Default for Sclip {
    fn default() -> Self {
        Self::new()
    }
}

/// SV mate data structure.
///
/// Ports `variations/Mate.java`.
#[derive(Debug, Clone)]
pub struct Mate {
    /// Mate start position.
    pub mate_start: i32,
    /// Mate end position.
    pub mate_end: i32,
    /// Mate length.
    pub mate_len: i32,
    /// Read start position.
    pub start: i32,
    /// Read end position.
    pub end: i32,
    /// Mean read position.
    pub pmean: f64,
    /// Mean base quality.
    pub qmean: f64,
    /// Mean mapping quality.
    pub mqmean: f64,
    /// Number of mismatches.
    pub nm: f64,
}

impl Mate {
    /// Create a new zeroed Mate.
    pub fn new() -> Self {
        Self {
            mate_start: 0,
            mate_end: 0,
            mate_len: 0,
            start: 0,
            end: 0,
            pmean: 0.0,
            qmean: 0.0,
            mqmean: 0.0,
            nm: 0.0,
        }
    }
}

impl Default for Mate {
    fn default() -> Self {
        Self::new()
    }
}

/// SV cluster data structure. Extends `Mate` with a count field.
///
/// Ports `variations/Cluster.java`.
#[derive(Debug, Clone, Default)]
pub struct Cluster {
    /// Mate data (inherited).
    pub mate: Mate,
    /// Count of mates in this cluster.
    pub cnt: i32,
}

impl Cluster {
    /// Create a new cluster with given fields.
    pub fn new(cnt: i32, mate: Mate) -> Self {
        Self { mate, cnt }
    }
}

/// Finalized output variant — separate from `Variation`.
/// Populated in Phase 3 (ToVarsBuilder).
///
/// Ports `variations/Variant.java`.
#[derive(Debug, Clone)]
pub struct Variant {
    /// Variant description string (`$n` — SNP letter, `+seq`, `-len`, etc.).
    pub description: String,
    /// Position coverage (`$cov`).
    pub pos_coverage: i32,
    /// Forward strand reads count.
    pub fwd: i32,
    /// Reverse strand reads count.
    pub rev: i32,
    /// Strand bias flag (0, 1, or 2).
    pub bias: String,
    /// Variant frequency.
    pub frequency: f64,
    /// Mean variant position in read.
    pub mean_position: f64,
    /// Found in ≥2 different positions.
    pub pstd: bool,
    /// Mean base quality.
    pub mean_quality: f64,
    /// Read with ≥2 different qualities.
    pub qstd: bool,
    /// Mean mapping quality.
    pub mean_mapping_quality: f64,
    /// Hi-quality/lo-quality reads ratio.
    pub hq_lq_ratio: f64,
    /// Variant frequency for HQ reads.
    pub hi_freq: f64,
    /// Adjusted AF for indels (realignment).
    pub extra_freq: f64,
    /// Bases to shift 3' for deletions.
    pub shift3: i32,
    /// Microsatellite instability (>1 = MSI).
    pub msi: f64,
    /// Microsatellite unit length in bp.
    pub msint: i32,
    /// Average mismatches for reads with variant.
    pub nm: f64,
    /// Number of high-quality reads with variant.
    pub hicnt: i32,
    /// Position coverage by high-quality reads.
    pub hicov: i32,
    /// Preceding reference sequence.
    pub left_seq: String,
    /// Following reference sequence.
    pub right_seq: String,
    /// Start position.
    pub start_pos: i32,
    /// End position.
    pub end_pos: i32,
    /// Reference reverse strand coverage.
    pub ref_rev_cov: i32,
    /// Reference forward strand coverage.
    pub ref_fwd_cov: i32,
    /// Total position coverage.
    pub total_pos_cov: i32,
    /// Duplication rate.
    pub duprate: f64,
    /// Genotype description.
    pub genotype: Option<String>,
    /// Variant allele (for VCF).
    pub var_allele: String,
    /// Reference allele (for VCF).
    pub ref_allele: String,
    /// Variant type.
    pub var_type: String,
    /// Debug information.
    pub debug: String,
    /// CRISPR information.
    pub crispr: i32,
    /// SV information.
    pub sv: String,
    /// True for Java's "null record": a row printed via `SimpleOutputVariant(null, …)`
    /// whose fields keep the printer defaults — genotype "" (:22), bias "0;0" (:24),
    /// flank seqs "" (OutputVariant.java:22-23) — instead of the normal
    /// null→"0" / empty→"0" mapping (SimpleOutputVariant.java:64, :81-82).
    pub null_record: bool,
}

impl Variant {
    /// Create a new Variant with default values.
    pub fn new() -> Self {
        let mut bias = String::with_capacity(4);
        bias.push('0');
        Self {
            description: String::new(),
            pos_coverage: 0,
            fwd: 0,
            rev: 0,
            bias,
            frequency: 0.0,
            mean_position: 0.0,
            pstd: false,
            mean_quality: 0.0,
            qstd: false,
            mean_mapping_quality: 0.0,
            hq_lq_ratio: 0.0,
            hi_freq: 0.0,
            extra_freq: 0.0,
            shift3: 0,
            msi: 0.0,
            msint: 0,
            nm: 0.0,
            hicnt: 0,
            hicov: 0,
            left_seq: String::new(),
            right_seq: String::new(),
            start_pos: 0,
            end_pos: 0,
            ref_rev_cov: 0,
            ref_fwd_cov: 0,
            total_pos_cov: 0,
            duprate: 0.0,
            genotype: None,
            var_allele: String::new(),
            ref_allele: String::new(),
            var_type: String::new(),
            debug: String::new(),
            crispr: 0,
            sv: String::new(),
            null_record: false,
        }
    }

    /// Per-read debug content string for -D output.
    /// Ports `Variant.debugVariantsContent` (Variant.java:268-282):
    /// `n:cnt:F-fwd:R-rev:freq:bias:meanPos:pstd:meanQual:qstd:hifreq:mapq:qratio`.
    /// Note mapq is printed via raw Java `Double.toString` semantics, not DecimalFormat.
    pub fn debug_variants_content(&self) -> String {
        let freq = format!("{:.4}", self.frequency);
        let mean_pos = format!("{:.1}", self.mean_position);
        let mean_qual = format!("{:.1}", self.mean_quality);
        let hi_freq = format!("{:.4}", self.hi_freq);
        let mapq = java_double_to_string(self.mean_mapping_quality);
        let qratio = format!("{:.3}", self.hq_lq_ratio);
        format!(
            "{}:{}:F-{}:R-{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
            self.description,
            self.fwd + self.rev,
            self.fwd,
            self.rev,
            freq,
            self.bias,
            mean_pos,
            if self.pstd { "1" } else { "0" },
            mean_qual,
            if self.qstd { "1" } else { "0" },
            hi_freq,
            mapq,
            qratio,
        )
    }

    /// Check if this variant is noise.
    pub fn is_noise(&mut self, goodq: f64, lofreq: f64) -> bool {
        let qual = self.mean_quality;
        if ((qual < 4.5 || (qual < 12.0 && !self.qstd)) && self.pos_coverage <= 3)
            || (qual < goodq && self.frequency < 2.0 * lofreq && self.pos_coverage <= 1)
        {
            self.total_pos_cov = self.total_pos_cov.saturating_sub(self.pos_coverage);
            self.pos_coverage = 0;
            self.fwd = 0;
            self.rev = 0;
            self.frequency = 0.0;
            self.hi_freq = 0.0;
            true
        } else {
            false
        }
    }

    /// Adjust coordinates and alleles for a complex variant.
    pub fn adj_complex(&mut self) {
        if self.var_allele.starts_with('<') {
            return;
        }

        // Part 1: prefix trimming
        let mut n = 0usize;
        while self.ref_allele.len().saturating_sub(n) > 1
            && self.var_allele.len().saturating_sub(n) > 1
        {
            if self.ref_allele.as_bytes().get(n) == self.var_allele.as_bytes().get(n) {
                n += 1;
            } else {
                break;
            }
        }
        if n > 0 && self.ref_allele.is_char_boundary(n) && self.var_allele.is_char_boundary(n) {
            self.start_pos = self.start_pos.saturating_add(n as i32);
            if let Some(prefix) = self.ref_allele.get(..n) {
                self.left_seq.push_str(prefix);
            }
            let left_drain = n.min(self.left_seq.len());
            if self.left_seq.is_char_boundary(left_drain) {
                self.left_seq.drain(..left_drain);
            }
            self.ref_allele.drain(..n);
            self.var_allele.drain(..n);
        }

        // Part 2: suffix trimming
        let ref_len = self.ref_allele.len();
        let var_len = self.var_allele.len();
        let mut n = 1usize;
        while ref_len.saturating_sub(n) >= 1 && var_len.saturating_sub(n) >= 1 {
            if self.ref_allele.as_bytes().get(ref_len - n)
                == self.var_allele.as_bytes().get(var_len - n)
            {
                n += 1;
            } else {
                break;
            }
        }
        if n > 1 {
            let trim = n - 1;
            self.end_pos = self.end_pos.saturating_sub(trim as i32);
            let new_ref_len = ref_len.saturating_sub(trim);
            let new_var_len = var_len.saturating_sub(trim);
            if self.ref_allele.is_char_boundary(new_ref_len)
                && self.var_allele.is_char_boundary(new_var_len)
            {
                let right_keep = self.right_seq.len().saturating_sub(trim);
                if self.right_seq.is_char_boundary(right_keep) {
                    self.right_seq.truncate(right_keep);
                }
                if let Some(ref_cut) = self.ref_allele.get(new_ref_len..) {
                    self.right_seq.insert_str(0, ref_cut);
                }
                self.ref_allele.truncate(new_ref_len);
                self.var_allele.truncate(new_var_len);
            }
        }
    }

    /// Find variant type based on ref and var alleles without heap allocation for standard types (`L7`).
    pub fn var_type_cow(&self) -> std::borrow::Cow<'static, str> {
        use std::borrow::Cow;
        if self.ref_allele == self.var_allele && self.ref_allele.len() == 1 {
            return Cow::Borrowed("");
        }
        if self.ref_allele.len() == 1 && self.var_allele.len() == 1 {
            return Cow::Borrowed("SNV");
        }
        // Match ANY_SV: "<(...)>"
        let bytes = self.var_allele.as_bytes();
        if bytes.is_ascii() {
            for w in bytes.windows(5) {
                if let [b'<', c1, c2, c3, b'>'] = *w {
                    if c1 != b'<'
                        && c1 != b'>'
                        && c2 != b'<'
                        && c2 != b'>'
                        && c3 != b'<'
                        && c3 != b'>'
                    {
                        let inner = [c1, c2, c3];
                        return match &inner {
                            b"DEL" => Cow::Borrowed("DEL"),
                            b"DUP" => Cow::Borrowed("DUP"),
                            b"INV" => Cow::Borrowed("INV"),
                            b"INS" => Cow::Borrowed("INS"),
                            b"FUS" => Cow::Borrowed("FUS"),
                            _ => Cow::Owned(String::from_utf8_lossy(&inner).into_owned()),
                        };
                    }
                }
            }
        } else {
            let chars: Vec<char> = self.var_allele.chars().collect();
            for w in chars.windows(5) {
                if let ['<', c1, c2, c3, '>'] = *w {
                    let inner = [c1, c2, c3];
                    if inner.iter().all(|&c| c != '<' && c != '>') {
                        return Cow::Owned(inner.iter().collect());
                    }
                }
            }
        }
        if self.ref_allele.is_empty() || self.var_allele.is_empty() {
            return Cow::Borrowed("Complex");
        }
        if self.ref_allele.as_bytes().first() != self.var_allele.as_bytes().first() {
            return Cow::Borrowed("Complex");
        }
        if self.ref_allele.len() == 1
            && self.var_allele.len() > 1
            && self.var_allele.starts_with(&self.ref_allele)
        {
            return Cow::Borrowed("Insertion");
        }
        if self.ref_allele.len() > 1
            && self.var_allele.len() == 1
            && self.ref_allele.starts_with(&self.var_allele)
        {
            return Cow::Borrowed("Deletion");
        }
        Cow::Borrowed("Complex")
    }

    /// Find variant type based on ref and var alleles.
    pub fn var_type(&self) -> String {
        self.var_type_cow().into_owned()
    }

    /// Check if variant meets standard parameters.
    pub fn is_good_var(
        &self,
        reference_var: Option<&Variant>,
        var_type: &str,
        splice: &rustc_hash::FxHashSet<String>,
        config: &crate::config::Config,
    ) -> bool {
        if self.ref_allele.is_empty() {
            return false;
        }
        let vtype = if var_type.is_empty() {
            self.var_type_cow()
        } else {
            std::borrow::Cow::Borrowed(var_type)
        };
        let vtype = vtype.as_ref();

        if self.frequency < config.freq
            || self.hicnt < config.min_reads as i32
            || self.mean_position < config.read_pos
            || self.mean_quality < config.base_qual
        {
            return false;
        }

        if let Some(ref_var) = reference_var {
            if ref_var.hicnt > config.min_reads as i32 && self.frequency < 0.25 {
                let d = self.mean_mapping_quality
                    + self.ref_allele.len() as f64
                    + self.var_allele.len() as f64;
                let f = (1.0 + d) / (ref_var.mean_mapping_quality + 1.0);
                if (d - 2.0 < 5.0 && ref_var.mean_mapping_quality > 20.0) || f < 0.25 {
                    return false;
                }
            }
        }

        if vtype == "Deletion" && !splice.is_empty() {
            let mut buf = [0u8; 32];
            let len = {
                let mut cursor = std::io::Cursor::new(&mut buf[..]);
                let _ = std::io::Write::write_fmt(
                    &mut cursor,
                    format_args!("{}-{}", self.start_pos, self.end_pos),
                );
                cursor.position() as usize
            };
            if let Some(bytes) = buf.get(..len) {
                if let Ok(key) = std::str::from_utf8(bytes) {
                    if splice.contains(key) {
                        return false;
                    }
                }
            }
        }

        if self.hq_lq_ratio < config.qratio {
            return false;
        }

        if self.frequency > 0.30 {
            return true;
        }

        if self.mean_mapping_quality < config.out_mapq {
            return false;
        }

        if self.msi >= 15.0 && self.frequency <= config.mfreq && self.msint == 1 {
            return false;
        }

        if self.msi >= 12.0 && self.frequency <= config.nmfreq && self.msint > 1 {
            return false;
        }

        if self.bias == "2;1"
            && self.frequency < 0.20
            && (vtype.is_empty()
                || vtype == "SNV"
                || (self.ref_allele.len() < 3 && self.var_allele.len() < 3))
        {
            return false;
        }

        true
    }
}

impl Default for Variant {
    fn default() -> Self {
        Self::new()
    }
}

/// Write an `f64` into `out` matching Java's `Double.toString(double)` (`FloatingDecimal` /
/// Steele & White) without intermediate heap allocations.
pub fn write_java_double_into<W: core::fmt::Write>(out: &mut W, v: f64) -> core::fmt::Result {
    use core::fmt::Write as _;
    if v.is_nan() {
        return out.write_str("NaN");
    }
    if v.is_infinite() {
        return if v.is_sign_positive() {
            out.write_str("Infinity")
        } else {
            out.write_str("-Infinity")
        };
    }
    if v == 0.0 {
        return if v.is_sign_negative() {
            out.write_str("-0.0")
        } else {
            out.write_str("0.0")
        };
    }
    let abs = v.abs();
    let mut tmp = crate::types::StackBuf::<64>::new();
    if (1e-3..1e7).contains(&abs) {
        write!(&mut tmp, "{v}")?;
        let s = tmp.as_str();
        out.write_str(s)?;
        if !s.contains('.') {
            out.write_str(".0")?;
        }
        Ok(())
    } else {
        write!(&mut tmp, "{v:e}")?;
        let s = tmp.as_str();
        if let Some((mantissa, exp)) = s.split_once('e') {
            out.write_str(mantissa)?;
            if !mantissa.contains('.') {
                out.write_str(".0")?;
            }
            out.write_str("E")?;
            out.write_str(exp)
        } else {
            out.write_str(s)
        }
    }
}

/// Format an f64 matching Java's `Double.toString(double)` (`FloatingDecimal` /
/// Steele & White):
/// - `NaN`, `Infinity`, `-Infinity`
/// - `0.0` -> `"0.0"`, `-0.0` -> `"-0.0"`
/// - `1e-3 <= |v| < 1e7`: plain decimal with `.0` suffix on integral values
/// - `|v| < 1e-3` or `|v| >= 1e7`: uppercase `E` scientific notation with at
///   least one digit after the decimal point (`1.0E-4`, `1.0E7`, `1.25E7`).
pub fn java_double_to_string(v: f64) -> String {
    let mut out = String::with_capacity(24);
    let _ = write_java_double_into(&mut out, v);
    out
}

/// Per-position variant container.
///
/// Ports `variations/Vars.java`.
#[derive(Debug, Clone)]
pub struct Vars {
    /// Reference variant.
    pub ref_variant: Option<Variant>,
    /// List of all non-reference variants.
    pub variants: Vec<Variant>,
    /// Variant description string → Variant.
    pub var_desc_to_variant: IndexMap<String, Variant, FxBuildHasher>,
    /// SV information string.
    pub sv: String,
}

impl Vars {
    /// Create a new Vars container.
    pub fn new() -> Self {
        Self {
            ref_variant: None,
            variants: Vec::new(),
            var_desc_to_variant: IndexMap::default(),
            sv: String::new(),
        }
    }

    /// Resolve a non-reference variant by description string, checking
    /// `var_desc_to_variant` when populated and falling back to reverse scan of `variants`.
    #[inline]
    pub fn variant_by_desc(&self, desc: &str) -> Option<&Variant> {
        if !self.var_desc_to_variant.is_empty() {
            if let Some(v) = self.var_desc_to_variant.get(desc) {
                return Some(v);
            }
        }
        self.variants.iter().rfind(|v| v.description == desc)
    }

    /// Update an existing variant matching `desc` in place without inserting into an
    /// empty `var_desc_to_variant` map.
    #[inline]
    pub fn update_variant_by_desc(&mut self, desc: &str, updated: Variant) {
        if let Some(slot) = self.var_desc_to_variant.get_mut(desc) {
            *slot = updated.clone();
        }
        if let Some(item) = self.variants.iter_mut().rfind(|v| v.description == desc) {
            *item = updated;
        }
    }
}

impl Default for Vars {
    fn default() -> Self {
        Self::new()
    }
}

/// Retrieve or insert `Vars` at `position` using `HashMap.get` + `HashMap.put` (`putVal`) semantics.
/// Ports `VariationUtils.getOrPutVars` (`VariationUtils.java:384-391`).
pub fn get_or_put_vars(
    map: &mut crate::javamap::JavaHashMap<i32, Vars>,
    position: i32,
) -> &mut Vars {
    map.get_or_put_with(position, Vars::new)
}

/// Structural variant info stored inside a VariationMap.
/// Ports `collection/VariationMap.java`'s static SV class.
#[derive(Debug, Clone, Default)]
pub struct Sv {
    pub r#type: String,
    pub pairs: i32,
    pub splits: i32,
    pub clusters: i32,
}

/// A wrapper around IndexMap that also supports structural variant info.
/// Ports `collection/VariationMap.java`.
#[derive(Debug, Clone, Default)]
pub struct VariationMap {
    pub map: IndexMap<String, Variation, FxBuildHasher>,
    pub sv: Option<Box<Sv>>,
    /// Set when a mate-window `partial_pipeline` merge plants this map's
    /// position outside the merging tile's active range (G3). The builder
    /// ownership gate (`merge_key_owned`, with edge clamping) then decides
    /// emission: single-region runs own everything so the gate never fires;
    /// tiled runs emit exactly once from the owning tile. Defaults false
    /// (native keys are never planted).
    pub merge_planted: bool,
}

impl VariationMap {
    pub fn new() -> Self {
        Self {
            map: IndexMap::default(),
            sv: None,
            merge_planted: false,
        }
    }

    /// Retrieve or create the SV entry for a given start position.
    /// Ports `VariationMap.getSV`.
    pub fn get_sv(
        hash: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
        start: i32,
    ) -> &mut Sv {
        let entry = hash.get_or_put_with(start, VariationMap::new);
        // §9: get_or_insert_with replaces the just-inserted unwrap; the "SV"
        // variation is ensured separately (remove_sv clears both together, so
        // sv==None implies no "SV" entry — the contains_key check preserves
        // the original insert-on-either-missing behavior).
        let sv = entry.sv.get_or_insert_with(|| Box::new(Sv::default()));
        if !entry.map.contains_key("SV") {
            entry.map.insert("SV".to_string(), Variation::new());
        }
        sv
    }

    /// Remove the SV entry and delete information about SV.
    /// Ports `VariationMap.removeSV`.
    pub fn remove_sv(hash: &mut crate::javamap::JavaHashMap<i32, VariationMap>, start: i32) {
        if let Some(entry) = hash.get_mut(&start) {
            entry.sv = None;
            entry.map.swap_remove("SV");
        }
    }
}

/// SVStructures holds SV-related collections.
/// Ports `data/SVStructures.java`.
#[derive(Debug, Clone, Default)]
pub struct SVStructures {
    pub svdelfend: i32,
    pub svdelrend: i32,
    pub svfdel: Vec<Sclip>,
    pub svrdel: Vec<Sclip>,

    pub svdupfend: i32,
    pub svduprend: i32,
    pub svfdup: Vec<Sclip>,
    pub svrdup: Vec<Sclip>,

    pub svinvfend3: i32,
    pub svinvrend3: i32,
    pub svinvfend5: i32,
    pub svinvrend5: i32,

    pub svfinv3: Vec<Sclip>,
    pub svfinv5: Vec<Sclip>,
    pub svrinv3: Vec<Sclip>,
    pub svrinv5: Vec<Sclip>,

    // svffus/svrfus preserve insertion order via IndexMap so is_best_sv_used
    // (sv.rs) can iterate keys in exact Java HashMap<String, List<Sclip>>
    // bucket + insertion order (SVStructures.java:37-40).
    pub svffus: indexmap::IndexMap<String, Vec<Sclip>>,
    pub svrfus: indexmap::IndexMap<String, Vec<Sclip>>,
    pub svfusfend: indexmap::IndexMap<String, i32>,
    pub svfusrend: indexmap::IndexMap<String, i32>,
}

impl SVStructures {
    pub fn new() -> Self {
        Self::default()
    }
}

/// The complete variation data accumulated during the CIGAR parse / pileup step.
/// Ports `data/scopedata/VariationData.java`.
#[derive(Debug, Clone, Default)]
pub struct VariationData {
    pub non_insertion_variants: crate::javamap::JavaHashMap<i32, VariationMap>,
    pub insertion_variants: crate::javamap::JavaHashMap<i32, VariationMap>,
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
    pub duprate: f64,
    pub sv_structures_filtered: bool,
    /// Owning tile's active range (pileup region), stamped at pileup time.
    /// Mate-window `partial_pipeline` merges (G3) plant keys outside this
    /// range; the ToVarsBuilder ownership gate drops those planted keys so
    /// each position emits from exactly one tile. `None` (manually built or
    /// pre-stamp data) fails open — the gate never fires. Metadata only:
    /// never serialized or emitted.
    pub owning_range: Option<(crate::types::Pos1, crate::types::Pos1)>,
    /// Positions added by mate-window `partial_pipeline` merges (G3).
    /// Keys already present before the merge are native, never planted.
    /// Metadata only.
    pub merge_planted: BTreeSet<i32>,
}

impl VariationData {
    pub fn new() -> Self {
        Self::default()
    }
}

pub fn get_variation<'a>(
    hash: &'a mut crate::javamap::JavaHashMap<i32, VariationMap>,
    start: i32,
    s: &str,
) -> &'a mut Variation {
    let entry = hash.get_or_put_with(start, VariationMap::new);
    match entry.map.raw_entry_mut_v1().from_key(s) {
        RawEntryMut::Occupied(e) => e.into_mut(),
        RawEntryMut::Vacant(v) => v.insert(s.to_owned(), Variation::new()).1,
    }
}

pub fn get_variation_maybe(
    hash: &crate::javamap::JavaHashMap<i32, VariationMap>,
    start: i32,
    ref_base: Option<u8>,
) -> Option<&Variation> {
    let ref_base = ref_base?;
    let entry = hash.get(&start)?;
    let mut buf = [0u8; 4];
    let k: &str = (ref_base as char).encode_utf8(&mut buf);
    entry.map.get(k)
}

pub fn get_variation_maybe_mut(
    hash: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    start: i32,
    ref_base: Option<u8>,
) -> Option<&mut Variation> {
    let ref_base = ref_base?;
    let entry = hash.get_mut(&start)?;
    let mut buf = [0u8; 4];
    let k: &str = (ref_base as char).encode_utf8(&mut buf);
    entry.map.get_mut(k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variation_add_sub_cnt() {
        let mut v = Variation::new();
        let goodq = 25.0;

        // Add a high-quality read
        v.add_cnt(false, 50, 30.0, 60, 2.0, goodq);
        assert_eq!(v.cnt, 1);
        assert_eq!(v.fwd, 1);
        assert_eq!(v.rev, 0);
        assert_eq!(v.pmean, 50.0);
        assert_eq!(v.qmean, 30.0);
        assert_eq!(v.mqmean, 60.0);
        assert_eq!(v.nm, 2.0);
        assert_eq!(v.hicnt, 1);
        assert_eq!(v.locnt, 0);

        // Add a low-quality read (reverse strand)
        v.add_cnt(true, 25, 10.0, 40, 1.0, goodq);
        assert_eq!(v.cnt, 2);
        assert_eq!(v.fwd, 1);
        assert_eq!(v.rev, 1);
        assert_eq!(v.pmean, 75.0);
        assert_eq!(v.qmean, 40.0);
        assert_eq!(v.hicnt, 1);
        assert_eq!(v.locnt, 1);

        // Subtract the low-quality read
        v.sub_cnt(true, 25, 10.0, 40, 1.0, goodq);
        assert_eq!(v.cnt, 1);
        assert_eq!(v.fwd, 1);
        assert_eq!(v.rev, 0);
        assert_eq!(v.pmean, 50.0);
        assert_eq!(v.qmean, 30.0);
        assert_eq!(v.hicnt, 1);
        assert_eq!(v.locnt, 0);
    }

    #[test]
    fn variation_direction() {
        let mut v = Variation::new();
        v.inc_dir(false);
        v.inc_dir(false);
        v.inc_dir(true);
        assert_eq!(v.get_dir(false), 2);
        assert_eq!(v.get_dir(true), 1);
        v.dec_dir(false);
        assert_eq!(v.get_dir(false), 1);
    }

    #[test]
    fn sclip_default() {
        let sc = Sclip::new();
        assert_eq!(sc.variation.cnt, 0);
        assert!(sc.nt.is_empty());
        assert!(sc.seq.is_empty());
        assert!(!sc.used);
    }

    #[test]
    fn mate_cluster_default() {
        let m = Mate::new();
        assert_eq!(m.mate_start, 0);
        let c = Cluster::new(5, m.clone());
        assert_eq!(c.cnt, 5);
        assert_eq!(c.mate.mate_start, 0);
    }

    #[test]
    fn variant_is_good_var() {
        let mut config = crate::config::Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;
        config.mfreq = 0.25;
        config.nmfreq = 0.1;

        let splice = rustc_hash::FxHashSet::default();

        // 1. Basic good variant
        let mut v = Variant::new();
        v.ref_allele = "A".to_string();
        v.var_allele = "G".to_string();
        v.frequency = 0.05;
        v.hicnt = 5;
        v.mean_position = 10.0;
        v.mean_quality = 30.0;
        v.hq_lq_ratio = 2.0;
        v.mean_mapping_quality = 30.0;
        v.bias = "0".to_string();
        assert!(v.is_good_var(None, "SNV", &splice, &config));

        // 2. Fails frequency threshold
        v.frequency = 0.005;
        assert!(!v.is_good_var(None, "SNV", &splice, &config));
        v.frequency = 0.05;

        // 3. Fails hicnt threshold
        v.hicnt = 1;
        assert!(!v.is_good_var(None, "SNV", &splice, &config));
        v.hicnt = 5;

        // 4. Early-accept (frequency > 0.3) overrides mapping quality & other checks
        v.frequency = 0.35;
        v.mean_mapping_quality = 5.0;
        assert!(v.is_good_var(None, "SNV", &splice, &config));

        // 5. Deletion splice check: reject when "{start_pos}-{end_pos}" in splice, accept otherwise
        let mut del = Variant::new();
        del.start_pos = 12345;
        del.end_pos = 12350;
        del.ref_allele = "ACGTAC".to_string();
        del.var_allele = "A".to_string();
        del.frequency = 0.50;
        del.hicnt = 10;
        del.mean_position = 15.0;
        del.mean_quality = 35.0;
        del.hq_lq_ratio = 5.0;
        del.mean_mapping_quality = 60.0;
        del.bias = "0".to_string();

        let mut splice_set = rustc_hash::FxHashSet::default();
        splice_set.insert("99999-99999".to_string());
        assert!(del.is_good_var(None, "Deletion", &splice_set, &config));
        assert!(del.is_good_var(None, "", &splice_set, &config));

        splice_set.insert("12345-12350".to_string());
        assert!(!del.is_good_var(None, "Deletion", &splice_set, &config));
        assert!(!del.is_good_var(None, "", &splice_set, &config));
    }

    #[test]
    fn variant_adj_complex_no_match_5prime_3prime_and_combined() {
        // 1. No bases match (n == 0 at 5', n == 1 at 3')
        let mut v = Variant::new();
        v.start_pos = 100;
        v.end_pos = 102;
        v.ref_allele = "ACG".to_string();
        v.var_allele = "TTA".to_string();
        v.left_seq = "GG".to_string();
        v.right_seq = "CC".to_string();
        v.adj_complex();
        assert_eq!(v.start_pos, 100);
        assert_eq!(v.end_pos, 102);
        assert_eq!(v.ref_allele, "ACG");
        assert_eq!(v.var_allele, "TTA");
        assert_eq!(v.left_seq, "GG");
        assert_eq!(v.right_seq, "CC");

        // 2. 5'-only trim
        let mut v5 = Variant::new();
        v5.start_pos = 100;
        v5.end_pos = 103;
        v5.ref_allele = "ACGT".to_string();
        v5.var_allele = "ACCA".to_string();
        v5.left_seq = "TTGG".to_string();
        v5.right_seq = "CC".to_string();
        v5.adj_complex();
        assert_eq!(v5.start_pos, 102);
        assert_eq!(v5.end_pos, 103);
        assert_eq!(v5.ref_allele, "GT");
        assert_eq!(v5.var_allele, "CA");
        assert_eq!(v5.left_seq, "GGAC");
        assert_eq!(v5.right_seq, "CC");

        // 3. 3'-only trim
        let mut v3 = Variant::new();
        v3.start_pos = 200;
        v3.end_pos = 203;
        v3.ref_allele = "TGAC".to_string();
        v3.var_allele = "CAAC".to_string();
        v3.left_seq = "TT".to_string();
        v3.right_seq = "GGCC".to_string();
        v3.adj_complex();
        assert_eq!(v3.start_pos, 200);
        assert_eq!(v3.end_pos, 201);
        assert_eq!(v3.ref_allele, "TG");
        assert_eq!(v3.var_allele, "CA");
        assert_eq!(v3.left_seq, "TT");
        assert_eq!(v3.right_seq, "ACGG");

        // 4. Combined 5' + 3' trim
        let mut vc = Variant::new();
        vc.start_pos = 300;
        vc.end_pos = 305;
        vc.ref_allele = "AATGGA".to_string();
        vc.var_allele = "AACCGA".to_string();
        vc.left_seq = "CCC".to_string();
        vc.right_seq = "TTT".to_string();
        vc.adj_complex();
        assert_eq!(vc.start_pos, 302);
        assert_eq!(vc.end_pos, 303);
        assert_eq!(vc.ref_allele, "TG");
        assert_eq!(vc.var_allele, "CC");
        assert_eq!(vc.left_seq, "CAA");
        assert_eq!(vc.right_seq, "GAT");

        // 5. Short and empty left_seq / right_seq edge cases (Variant.java:237-238, 252-254)
        let mut ve = Variant::new();
        ve.start_pos = 400;
        ve.end_pos = 404;
        ve.ref_allele = "ACGTA".to_string();
        ve.var_allele = "ACCTA".to_string();
        ve.left_seq = String::new();
        ve.right_seq = String::new();
        ve.adj_complex();
        assert_eq!(ve.start_pos, 402);
        assert_eq!(ve.end_pos, 402);
        assert_eq!(ve.ref_allele, "G");
        assert_eq!(ve.var_allele, "C");
        assert_eq!(ve.left_seq, "");
        assert_eq!(ve.right_seq, "TA");

        let mut vs = Variant::new();
        vs.start_pos = 500;
        vs.end_pos = 504;
        vs.ref_allele = "ACGTA".to_string();
        vs.var_allele = "ACCTA".to_string();
        vs.left_seq = "T".to_string();
        vs.right_seq = "G".to_string();
        vs.adj_complex();
        assert_eq!(vs.start_pos, 502);
        assert_eq!(vs.end_pos, 502);
        assert_eq!(vs.ref_allele, "G");
        assert_eq!(vs.var_allele, "C");
        // In Java: leftSeq = ("T" + "AC").substring(2) -> "C"
        assert_eq!(vs.left_seq, "C");
        // In Java: rightSeq = substr("G", 0, -2) -> "", then "TA" + "" -> "TA"
        assert_eq!(vs.right_seq, "TA");
    }

    #[test]
    fn get_variation_hit_does_not_reorder_or_duplicate_and_latin1_matches() {
        let mut m: crate::javamap::JavaHashMap<i32, VariationMap> =
            crate::javamap::JavaHashMap::new();
        get_variation(&mut m, 100, "A").cnt = 3;
        get_variation(&mut m, 100, "+T").cnt = 1;
        get_variation(&mut m, 100, "C").cnt = 2;
        get_variation(&mut m, 100, "A").cnt += 4;

        let keys: Vec<&str> = m
            .get(&100)
            .unwrap()
            .map
            .keys()
            .map(String::as_str)
            .collect();
        assert_eq!(keys, vec!["A", "+T", "C"]);
        assert_eq!(get_variation_maybe(&m, 100, Some(b'A')).unwrap().cnt, 7);

        let latin1_key = (0xE9u8 as char).to_string();
        get_variation(&mut m, 100, &latin1_key).cnt = 9;
        assert_eq!(get_variation_maybe(&m, 100, Some(0xE9)).unwrap().cnt, 9);
        get_variation_maybe_mut(&mut m, 100, Some(0xE9))
            .unwrap()
            .cnt = 11;
        assert_eq!(get_variation_maybe(&m, 100, Some(0xE9)).unwrap().cnt, 11);
    }

    #[test]
    fn get_variation_two_allele_low_af_keeps_insertion_order() {
        // H1 pin: two-allele low-AF pileup (higher-AF SNV + low-AF insertion);
        // hits must accumulate without reorder or duplicate keys.
        let mut m: crate::javamap::JavaHashMap<i32, VariationMap> =
            crate::javamap::JavaHashMap::new();
        get_variation(&mut m, 50, "G").cnt = 40;
        get_variation(&mut m, 50, "+AT").cnt = 2;
        get_variation(&mut m, 50, "G").cnt += 1;
        get_variation(&mut m, 50, "+AT").cnt += 1;

        let keys: Vec<&str> = m.get(&50).unwrap().map.keys().map(String::as_str).collect();
        assert_eq!(keys, vec!["G", "+AT"]);
        assert_eq!(m.get(&50).unwrap().map.get("G").unwrap().cnt, 41);
        assert_eq!(m.get(&50).unwrap().map.get("+AT").unwrap().cnt, 3);
    }
}

#[cfg(test)]
mod debug_tests {
    use super::*;

    #[test]
    fn debug_variants_content_matches_java_format() {
        // ports Variant.java:268-282 field order and formats
        let mut v = Variant::new();
        v.description = "G".to_string();
        v.fwd = 2;
        v.rev = 4;
        v.frequency = 0.75;
        v.bias = "0;2".to_string();
        v.mean_position = 17.0;
        v.pstd = true;
        v.mean_quality = 36.7;
        v.qstd = true;
        v.hi_freq = 1.0;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 12.0;
        assert_eq!(
            v.debug_variants_content(),
            "G:6:F-2:R-4:0.7500:0;2:17.0:1:36.7:1:1.0000:60.0:12.000"
        );
    }

    #[test]
    fn debug_variants_content_insertion_prefix() {
        // ports Variant.debugVariantsContentInsertion (Variant.java:261-266).
        // hi_freq 1.0/6.0 -> "0.1667" pins Java DecimalFormat("0.0000")
        // (HALF_EVEN default) nearest-rounding; note it exercises NO tie, so
        // HALF_EVEN vs HALF_UP is not distinguished by this case.
        let mut v = Variant::new();
        v.description = "+GGAGG".to_string();
        v.fwd = 1;
        v.rev = 1;
        v.frequency = 0.25;
        v.bias = "2".to_string();
        v.mean_position = 22.0;
        v.mean_quality = 31.5;
        v.hi_freq = 1.0 / 6.0;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 1.0;
        assert_eq!(
            format!("I{}", v.debug_variants_content()),
            "I+GGAGG:2:F-1:R-1:0.2500:2:22.0:0:31.5:0:0.1667:60.0:1.000"
        );
    }

    #[test]
    fn java_double_to_string_appends_zero_decimal() {
        assert_eq!(java_double_to_string(60.0), "60.0");
        assert_eq!(java_double_to_string(35.3), "35.3");
        assert_eq!(java_double_to_string(0.0), "0.0");
        assert_eq!(java_double_to_string(-0.0), "-0.0");
        assert_eq!(java_double_to_string(1.0e-4), "1.0E-4");
        assert_eq!(java_double_to_string(-1.0e-4), "-1.0E-4");
        assert_eq!(java_double_to_string(6.0e-5), "6.0E-5");
        assert_eq!(java_double_to_string(6.2e-4), "6.2E-4");
        assert_eq!(java_double_to_string(1.6e-4), "1.6E-4");
        assert_eq!(java_double_to_string(6.0), "6.0");
        assert_eq!(java_double_to_string(16384.0), "16384.0");
        assert_eq!(java_double_to_string(1.0e7), "1.0E7");
        assert_eq!(java_double_to_string(1.25e7), "1.25E7");
        assert_eq!(java_double_to_string(f64::NAN), "NaN");
        assert_eq!(java_double_to_string(f64::INFINITY), "Infinity");
        assert_eq!(java_double_to_string(f64::NEG_INFINITY), "-Infinity");
    }

    #[test]
    fn write_java_double_into_matches_java_double_to_string() {
        let cases = [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.5, "1.5"),
            (60.0, "60.0"),
            (1e-5, "1.0E-5"),
            (-1e-4, "-1.0E-4"),
            (1e8, "1.0E8"),
            (1.25e7, "1.25E7"),
            (4.503599627370496e15, "4.503599627370496E15"),
            (f64::INFINITY, "Infinity"),
            (f64::NEG_INFINITY, "-Infinity"),
            (f64::NAN, "NaN"),
        ];

        for (v, pinned) in cases {
            let expected = java_double_to_string(v);
            assert_eq!(expected, pinned, "pinned mismatch for {v:?}");

            let mut s = String::new();
            write_java_double_into(&mut s, v).unwrap();
            assert_eq!(s, expected, "String mismatch for {v:?}");

            let mut prefixed = "prefix=".to_string();
            write_java_double_into(&mut prefixed, v).unwrap();
            assert_eq!(prefixed, format!("prefix={expected}"));

            let mut stack = crate::types::StackBuf::<32>::new();
            write_java_double_into(&mut stack, v).unwrap();
            assert_eq!(stack.as_str(), expected, "StackBuf mismatch for {v:?}");
        }
    }

    #[test]
    fn vars_variant_by_desc_and_update_variant_by_desc_coherence() {
        // 1. Both var_desc_to_variant and variants populated
        let mut vars_both = Vars::new();
        let mut v_a1 = Variant::new();
        v_a1.description = "A".to_string();
        v_a1.pos_coverage = 3;
        let mut v_a2 = Variant::new();
        v_a2.description = "A".to_string();
        v_a2.pos_coverage = 7;
        let mut v_only_vec = Variant::new();
        v_only_vec.description = "+TT".to_string();
        v_only_vec.pos_coverage = 5;

        vars_both.variants.push(v_a1);
        vars_both.variants.push(v_a2.clone());
        vars_both.variants.push(v_only_vec);
        vars_both.var_desc_to_variant.insert("A".to_string(), v_a2);

        assert_eq!(vars_both.variant_by_desc("A").unwrap().pos_coverage, 7);
        // Missing from non-empty map falls back to rfind on variants
        assert_eq!(vars_both.variant_by_desc("+TT").unwrap().pos_coverage, 5);
        assert!(vars_both.variant_by_desc("-3").is_none());

        let mut updated_a = Variant::new();
        updated_a.description = "A".to_string();
        updated_a.pos_coverage = 42;
        vars_both.update_variant_by_desc("A", updated_a);

        assert_eq!(
            vars_both.var_desc_to_variant.get("A").unwrap().pos_coverage,
            42
        );
        assert_eq!(vars_both.variants[0].pos_coverage, 3);
        assert_eq!(vars_both.variants[1].pos_coverage, 42);
        assert_eq!(vars_both.variant_by_desc("A").unwrap().pos_coverage, 42);

        // 2. Empty var_desc_to_variant (ToVarsBuilder mode): rfind fallback + no map insertion
        let mut vars_vec_only = Vars::new();
        let mut v1 = Variant::new();
        v1.description = "G".to_string();
        v1.pos_coverage = 2;
        let mut v2 = Variant::new();
        v2.description = "G".to_string();
        v2.pos_coverage = 9;
        vars_vec_only.variants.push(v1);
        vars_vec_only.variants.push(v2);

        assert!(vars_vec_only.var_desc_to_variant.is_empty());
        assert_eq!(vars_vec_only.variant_by_desc("G").unwrap().pos_coverage, 9);
        assert!(vars_vec_only.variant_by_desc("C").is_none());

        let mut updated_g = Variant::new();
        updated_g.description = "G".to_string();
        updated_g.pos_coverage = 99;
        vars_vec_only.update_variant_by_desc("G", updated_g);

        assert!(
            vars_vec_only.var_desc_to_variant.is_empty(),
            "update_variant_by_desc must not populate an empty var_desc_to_variant"
        );
        assert_eq!(vars_vec_only.variants[0].pos_coverage, 2);
        assert_eq!(vars_vec_only.variants[1].pos_coverage, 99);
        assert_eq!(vars_vec_only.variant_by_desc("G").unwrap().pos_coverage, 99);
    }
}
