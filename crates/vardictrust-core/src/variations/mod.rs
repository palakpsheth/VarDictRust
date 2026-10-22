//! Variation data structures — intermediate variant accumulators and SV helpers.
//!
//! These port the core data types from VarDictJava's `variations/` package:
//! - `Variation`: intermediate variant accumulator (sums, not means)
//! - `Sclip`: soft-clip data, extends Variation
//! - `Mate`: SV mate data
//! - `Cluster`: SV cluster data, extends Mate
//! - `Variant`: finalized output variant (populated in Phase 3)
//! - `Vars`: per-position variant container

use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet};

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
        Self {
            description: String::new(),
            pos_coverage: 0,
            fwd: 0,
            rev: 0,
            bias: String::from("0"),
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
        let mut ref_allele = self.ref_allele.clone();
        let mut var_allele = self.var_allele.clone();
        let mut left_seq = self.left_seq.clone();
        let mut right_seq = self.right_seq.clone();
        let mut start_pos = self.start_pos;
        let mut end_pos = self.end_pos;

        // Part 1: prefix trimming
        {
            let ref_bytes = ref_allele.as_bytes();
            let var_bytes = var_allele.as_bytes();
            let mut n = 0;
            while ref_bytes.len().saturating_sub(n) > 1
                && var_bytes.len().saturating_sub(n) > 1
                && ref_bytes.get(n) == var_bytes.get(n)
            {
                n += 1;
            }
            if n > 0 {
                start_pos = start_pos.saturating_add(n as i32);
                if let (Some(ref_suffix_bytes), Some(var_suffix_bytes)) =
                    (ref_bytes.get(n..), var_bytes.get(n..))
                {
                    let ref_suffix = String::from_utf8_lossy(ref_suffix_bytes).into_owned();
                    let var_suffix = String::from_utf8_lossy(var_suffix_bytes).into_owned();
                    if let Some(prefix_bytes) = ref_bytes.get(..n) {
                        let prefix = String::from_utf8_lossy(prefix_bytes);
                        left_seq.push_str(&prefix);
                        if left_seq.len() > n {
                            left_seq = left_seq.get(n..).unwrap_or("").to_string();
                        } else {
                            left_seq.clear();
                        }
                    }
                    ref_allele = ref_suffix;
                    var_allele = var_suffix;
                }
            }
        }

        // Part 2: suffix trimming
        {
            let ref_bytes = ref_allele.as_bytes();
            let var_bytes = var_allele.as_bytes();
            let mut n = 1;
            while ref_bytes.len() > n
                && var_bytes.len() > n
                && ref_bytes.get(ref_bytes.len() - n) == var_bytes.get(var_bytes.len() - n)
            {
                n += 1;
            }
            if n > 1 {
                let trim = n - 1;
                end_pos = end_pos.saturating_sub(trim as i32);
                let ref_trimmed = ref_bytes
                    .get(..ref_bytes.len().saturating_sub(trim))
                    .map(|b| String::from_utf8_lossy(b).into_owned())
                    .unwrap_or_default();
                let var_trimmed = var_bytes
                    .get(..var_bytes.len().saturating_sub(trim))
                    .map(|b| String::from_utf8_lossy(b).into_owned())
                    .unwrap_or_default();
                let ref_cut = ref_bytes
                    .get(ref_bytes.len().saturating_sub(trim)..)
                    .map(|b| String::from_utf8_lossy(b).into_owned())
                    .unwrap_or_default();

                if right_seq.len() >= trim {
                    right_seq = right_seq
                        .get(..right_seq.len().saturating_sub(trim))
                        .unwrap_or("")
                        .to_string();
                } else {
                    right_seq.clear();
                }
                right_seq = format!("{}{}", ref_cut, right_seq);
                ref_allele = ref_trimmed;
                var_allele = var_trimmed;
            }
        }

        self.ref_allele = ref_allele;
        self.var_allele = var_allele;
        self.left_seq = left_seq;
        self.right_seq = right_seq;
        self.start_pos = start_pos;
        self.end_pos = end_pos;
    }

    /// Find variant type based on ref and var alleles.
    pub fn var_type(&self) -> String {
        if self.ref_allele == self.var_allele && self.ref_allele.len() == 1 {
            return String::new();
        }
        if self.ref_allele.len() == 1 && self.var_allele.len() == 1 {
            return "SNV".to_string();
        }
        // Match ANY_SV: "<(...)>"
        let chars: Vec<char> = self.var_allele.chars().collect();
        if chars.len() >= 5 {
            for i in 0..=chars.len() - 5 {
                if chars[i] == '<' && chars[i + 4] == '>' {
                    let inner = &chars[i + 1..i + 4];
                    if inner.iter().all(|&c| c != '<' && c != '>') {
                        return inner.iter().collect();
                    }
                }
            }
        }
        if self.ref_allele.is_empty() || self.var_allele.is_empty() {
            return "Complex".to_string();
        }
        if self.ref_allele.as_bytes()[0] != self.var_allele.as_bytes()[0] {
            return "Complex".to_string();
        }
        if self.ref_allele.len() == 1
            && self.var_allele.len() > 1
            && self.var_allele.starts_with(&self.ref_allele)
        {
            return "Insertion".to_string();
        }
        if self.ref_allele.len() > 1
            && self.var_allele.len() == 1
            && self.ref_allele.starts_with(&self.var_allele)
        {
            return "Deletion".to_string();
        }
        "Complex".to_string()
    }

    /// Check if variant meets standard parameters.
    pub fn is_good_var(
        &self,
        reference_var: Option<&Variant>,
        var_type: &str,
        splice: &std::collections::HashSet<String>,
        config: &crate::config::Config,
    ) -> bool {
        if self.ref_allele.is_empty() {
            return false;
        }
        let vtype = if var_type.is_empty() {
            self.var_type()
        } else {
            var_type.to_string()
        };

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

        if vtype == "Deletion" && splice.contains(&format!("{}-{}", self.start_pos, self.end_pos)) {
            return false;
        }

        if self.hq_lq_ratio < config.qratio {
            return false;
        }

        if self.frequency > 0.30 {
            return true;
        }

        if self.mean_mapping_quality < config.out_mapq as f64 {
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

/// Format an f64 the way Java's `Double.toString` does for the value domain seen
/// in debug output (mapping quality): integral values keep a ".0" suffix
/// ("60.0"), non-integral values use shortest round-trip repr ("35.3").
/// Caveat: exact only for |v| < 1e7 — above that Java switches to scientific
/// notation ("1.0E7") while this still prints plain digits; such magnitudes do
/// not occur in the fields routed through here.
pub fn java_double_to_string(v: f64) -> String {
    if v.is_finite() && v.fract() == 0.0 && v.abs() < 1e7 {
        format!("{v:.1}")
    } else {
        format!("{v}")
    }
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
    pub var_desc_to_variant: IndexMap<String, Variant>,
    /// SV information string.
    pub sv: String,
}

impl Vars {
    /// Create a new Vars container.
    pub fn new() -> Self {
        Self {
            ref_variant: None,
            variants: Vec::new(),
            var_desc_to_variant: IndexMap::new(),
            sv: String::new(),
        }
    }
}

impl Default for Vars {
    fn default() -> Self {
        Self::new()
    }
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
    pub map: IndexMap<String, Variation>,
    pub sv: Option<Sv>,
}

impl VariationMap {
    pub fn new() -> Self {
        Self {
            map: IndexMap::new(),
            sv: None,
        }
    }

    /// Retrieve or create the SV entry for a given start position.
    /// Ports `VariationMap.getSV`.
    pub fn get_sv(
        hash: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
        start: i32,
    ) -> &mut Sv {
        let entry = hash.compute_if_absent(start, |_| VariationMap::new());
        // §9: get_or_insert_with replaces the just-inserted unwrap; the "SV"
        // variation is ensured separately (remove_sv clears both together, so
        // sv==None implies no "SV" entry — the contains_key check preserves
        // the original insert-on-either-missing behavior).
        let sv = entry.sv.get_or_insert_with(Sv::default);
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

    // svffus/svrfus/svfusfend/svfusrend: iteration-order audit (W5.3,
    // 2026-08-23) — see the shared comment on CigarParser's containers in
    // cigar/parse.rs. Keyed lookups + per-list ops only; HashMap in
    // SVStructures.java:38-40.
    pub svffus: std::collections::HashMap<String, Vec<Sclip>>,
    pub svrfus: std::collections::HashMap<String, Vec<Sclip>>,
    pub svfusfend: std::collections::HashMap<String, i32>,
    pub svfusrend: std::collections::HashMap<String, i32>,
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
    pub position_to_insertion_count: BTreeMap<i32, std::collections::HashMap<String, i32>>,
    pub position_to_deletion_count: BTreeMap<i32, std::collections::HashMap<String, i32>>,
    pub sv_structures: SVStructures,
    pub ref_coverage: BTreeMap<i32, i32>,
    pub soft_clips_5_end: BTreeMap<i32, Sclip>,
    pub soft_clips_3_end: BTreeMap<i32, Sclip>,
    pub max_read_length: i32,
    pub splice: std::collections::HashSet<String>,
    pub mnp: BTreeMap<i32, std::collections::HashMap<String, i32>>,
    pub splice_count: std::collections::HashMap<String, [i32; 1]>,
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
    let entry = hash.compute_if_absent(start, |_| VariationMap::new());
    entry
        .map
        .entry(s.to_string())
        .or_insert_with(Variation::new)
}

pub fn get_variation_maybe(
    hash: &crate::javamap::JavaHashMap<i32, VariationMap>,
    start: i32,
    ref_base: Option<u8>,
) -> Option<&Variation> {
    let ref_base = ref_base?;
    let entry = hash.get(&start)?;
    let key = (ref_base as char).to_string();
    entry.map.get(&key)
}

pub fn get_variation_maybe_mut(
    hash: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    start: i32,
    ref_base: Option<u8>,
) -> Option<&mut Variation> {
    let ref_base = ref_base?;
    let entry = hash.get_mut(&start)?;
    let key = (ref_base as char).to_string();
    entry.map.get_mut(&key)
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

        let splice = std::collections::HashSet::new();

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
    }
}
