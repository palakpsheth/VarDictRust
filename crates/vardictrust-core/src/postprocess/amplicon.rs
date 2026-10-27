//! Amplicon post-processing module — processes variants from multiple amplicons
//! at each position and applies amplicon bias filtering.
//!
//! Ports `AmpliconPostProcessModule.java` from VarDictJava v1.8.3.
//!
//! An amplicon BED file has 8 columns: chr, start, end, gene, ampStart, ampEnd, score, strand.
//! Multiple amplicons can overlap at the same position. For each position, variants from all
//! overlapping amplicons are collected. Amplicon bias is detected when different amplicons
//! show different variants at the same position.
//!
//! Output adds 4 extra columns: GDAMP (good amplicon count), TLAMP (total amplicon count),
//! NCAMP (no-coverage amplicon count), AMPFLAG (amplicon bias flag).

use crate::config::Config;
use crate::regions::Region;
use crate::variations::{Variant, Vars};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashSet};

/// Output of amplicon post-processing: a variant with amplicon-specific metrics.
#[derive(Debug, Clone)]
pub struct AmpliconVariant {
    /// The finalized variant.
    pub variant: Variant,
    /// Good amplicon count — number of amplicons where this variant passes isGoodVar.
    pub gdamp: i32,
    /// Total amplicon count — total number of amplicons at this position.
    pub tlamp: i32,
    /// No-coverage amplicon count — amplicons with coverage < maxcov/50.
    pub ncamp: i32,
    /// Amplicon bias flag — true if different amplicons have different good variants.
    pub ampflag: bool,
    /// Output region string (first good variant's amplicon region, or chr:pos-pos).
    pub region: String,
    /// Good amplicon variants + their region strings (-D debug column only; empty when off).
    /// Ports the `goodVariants` list passed to AmpliconOutputVariant (AmpliconOutputVariant.java:48).
    pub good_variants: Vec<(Option<Variant>, String)>,
    /// Bad amplicon variants + their region strings (-D debug column only; empty when off).
    /// Entries may be None where an amplicon has no data at the position
    /// (AmpliconPostProcessModule.java:178 adds a null variant).
    pub bad_variants: Vec<(Option<Variant>, String)>,
}

/// Processes variants from multiple amplicons at each position and applies
/// amplicon bias filtering.
///
/// # Arguments
/// * `config` - Global configuration
/// * `region` - The genomic region being processed
/// * `vars_per_amplicon` - Per-amplicon variant maps (position → Vars)
/// * `amplicons_on_positions` - Map of position → list of (amplicon index, region)
/// * `splice` - Set of splice site strings
///
/// Ports `AmpliconPostProcessModule.process()`.
pub fn amplicon_post_process(
    config: &Config,
    region: &Region,
    vars_per_amplicon: &[crate::javamap::JavaHashMap<i32, Vars>],
    amplicons_on_positions: &BTreeMap<i32, Vec<(usize, Region)>>,
    splice: &HashSet<String>,
) -> Result<Vec<AmpliconVariant>, crate::types::VardictError> {
    // Batch wrapper over the streaming core below: the BTreeMap iterator visits
    // positions ascending with segment-ordered covering lists — the same
    // sequence the lazy -p union stream yields — so rows are identical.
    let mut results = Vec::new();
    {
        let mut collect = |av: AmpliconVariant| -> Result<(), crate::types::VardictError> {
            results.push(av);
            Ok(())
        };
        amplicon_post_process_streaming(
            config,
            region,
            vars_per_amplicon,
            amplicons_on_positions
                .iter()
                .map(|(&position, regions)| (position, regions.clone())),
            splice,
            &mut collect,
        )?;
    }

    if !config.emit_intermediate {
        results = amplicon_select_vcf_rows(results);
    }

    Ok(results)
}

/// Active-interval sweep over a segment's insert intervals.
///
/// Stable-sorts amplicon indices by insert start once, then answers
/// `covering_at(p)` for positions arriving in ascending order: a
/// start-pointer admits newly-covering amplicons while a min-heap of
/// `(end, seg_idx)` evicts expired ones. The live set is a `BTreeSet` keyed
/// by segment index, so each covering list emerges in segment order —
/// exactly the order the old per-position scan produced — at a total cost of
/// O((P + A) log A) heap/set operations plus output-proportional cloning,
/// never O(P x A) predicates.
///
/// Precondition: `covering_at` positions MUST arrive in ascending
/// (non-decreasing) order (debug-asserted). Both callers satisfy this:
/// `amplicon_position_stream` walks merged ranges ascending, and the
/// non-`-p` path in `run_amplicon` iterates a `BTreeSet`.
pub(crate) struct AmpliconCoverSweep<'a> {
    segment: &'a [Region],
    by_start: Vec<usize>,
    starts: Vec<i32>,
    ends: Vec<i32>,
    next_start: usize,
    active: BTreeSet<usize>,
    expiry: BinaryHeap<Reverse<(i32, usize)>>,
    #[cfg(debug_assertions)]
    last_pos: Option<i32>,
}

impl<'a> AmpliconCoverSweep<'a> {
    pub(crate) fn new(segment: &'a [Region]) -> Self {
        let mut starts = Vec::with_capacity(segment.len());
        let mut ends = Vec::with_capacity(segment.len());
        let mut by_start = Vec::new();
        for (idx, region) in segment.iter().enumerate() {
            match region.amplicon.as_ref() {
                Some(amp) => {
                    // Same i32 conversion as the old scan (amplicons_covering_position).
                    starts.push(amp.insert_start.0 as i32);
                    ends.push(amp.insert_end.0 as i32);
                    by_start.push(idx);
                }
                None => {
                    // Covers no position: sentinels the sweep never admits.
                    starts.push(i32::MAX);
                    ends.push(i32::MIN);
                }
            }
        }
        // Stable sort: tie starts keep segment order, matching the old scan.
        by_start.sort_by_key(|idx| match starts.get(*idx) {
            Some(start) => *start,
            None => i32::MAX,
        });
        Self {
            segment,
            by_start,
            starts,
            ends,
            next_start: 0,
            active: BTreeSet::new(),
            expiry: BinaryHeap::new(),
            #[cfg(debug_assertions)]
            last_pos: None,
        }
    }

    /// Covering amplicons at `position`, in segment order. Positions MUST
    /// arrive in ascending (non-decreasing) order.
    pub(crate) fn covering_at(&mut self, position: i32) -> Vec<(usize, Region)> {
        #[cfg(debug_assertions)]
        {
            if let Some(last) = self.last_pos {
                debug_assert!(
                    position >= last,
                    "AmpliconCoverSweep requires ascending positions"
                );
            }
            self.last_pos = Some(position);
        }
        while let Some(idx) = self.by_start.get(self.next_start) {
            let idx = *idx;
            let start = match self.starts.get(idx) {
                Some(start) => *start,
                None => break,
            };
            if start > position {
                break;
            }
            self.active.insert(idx);
            if let Some(end) = self.ends.get(idx).copied() {
                self.expiry.push(Reverse((end, idx)));
            }
            self.next_start = self.next_start.saturating_add(1);
        }
        while let Some(Reverse((end, _))) = self.expiry.peek().copied() {
            if end >= position {
                break;
            }
            if let Some(Reverse((_, idx))) = self.expiry.pop() {
                self.active.remove(&idx);
            }
        }
        let mut covering = Vec::with_capacity(self.active.len());
        for idx in &self.active {
            if let Some(region) = self.segment.get(*idx) {
                covering.push((*idx, region.clone()));
            }
        }
        covering
    }
}

/// Lazy iterator behind `amplicon_position_stream`: walks the merged insert
/// ranges ascending, resolving each position's covering list through one
/// shared `AmpliconCoverSweep` (ascending by construction).
struct AmpliconPositionIter<'a> {
    sweep: AmpliconCoverSweep<'a>,
    merged: Vec<(i32, i32)>,
    range_idx: usize,
    current: i32,
}

impl<'a> Iterator for AmpliconPositionIter<'a> {
    type Item = (i32, Vec<(usize, Region)>);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let (start, end) = *self.merged.get(self.range_idx)?;
            if self.current < start {
                self.current = start;
            }
            if self.current > end {
                // Past this range (or a degenerate start > end, which the
                // old `start..=end` expansion also yielded nothing for).
                self.range_idx = self.range_idx.saturating_add(1);
                continue;
            }
            let position = self.current;
            let covering = self.sweep.covering_at(position);
            if position >= end {
                self.range_idx = self.range_idx.saturating_add(1);
            } else {
                // position < end <= i32::MAX here: cannot overflow.
                self.current = position + 1;
            }
            return Some((position, covering));
        }
    }
}

/// Lazy ascending union of the segment's insert intervals as
/// `(position, covering-amplicons)` pairs.
///
/// Streaming form of the per-base expansion in AmpliconMode.java:57-60
/// (`for p in insertStart..insertEnd`): yields exactly the positions — and the
/// same segment-ordered covering lists — the materialized `BTreeMap` held,
/// without one map entry per base (audit-amplicon-expand). Disjoint ranges are
/// merged first, so gaps between amplicons are never visited. Covering lists
/// resolve through one `AmpliconCoverSweep`: O((P + A) log A) total, never
/// O(P x A) predicates. Working memory is O(segment); positions emerge in
/// ascending order.
pub fn amplicon_position_stream(
    segment: &[Region],
) -> impl Iterator<Item = (i32, Vec<(usize, Region)>)> + '_ {
    let mut ranges: Vec<(i32, i32)> = segment
        .iter()
        .filter_map(|region| {
            region
                .amplicon
                .as_ref()
                .map(|amp| (amp.insert_start.0 as i32, amp.insert_end.0 as i32))
        })
        .collect();
    ranges.sort();
    let mut merged: Vec<(i32, i32)> = Vec::with_capacity(ranges.len());
    for (start, end) in ranges {
        if let Some(last) = merged.last_mut() {
            // Merge overlapping or adjacent ranges: every integer between them
            // is covered by construction, so no position is lost or repeated.
            if start <= last.1.saturating_add(1) {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }
    AmpliconPositionIter {
        sweep: AmpliconCoverSweep::new(segment),
        merged,
        range_idx: 0,
        current: 0,
    }
}

/// Streaming core of `AmpliconPostProcessModule.process()`: one iteration per
/// position, each output row handed to `emit` immediately (the Java module
/// likewise prints per position via `variantPrinter`, without buffering rows).
///
/// `positions` MUST arrive in ascending order with segment-ordered covering
/// lists — both the batch `BTreeMap` and `amplicon_position_stream` satisfy
/// this — so row order matches the batch path exactly.
pub fn amplicon_post_process_streaming(
    config: &Config,
    region: &Region,
    vars_per_amplicon: &[crate::javamap::JavaHashMap<i32, Vars>],
    positions: impl IntoIterator<Item = (i32, Vec<(usize, Region)>)>,
    splice: &HashSet<String>,
    emit: &mut dyn FnMut(AmpliconVariant) -> Result<(), crate::types::VardictError>,
) -> Result<(), crate::types::VardictError> {
    // Scratch collections reused across positions to avoid per-base allocations.
    // good variants: (variant, region_string)
    let mut gvs: Vec<(Variant, String)> = Vec::new();
    // reference variants
    let mut ref_vars: Vec<Variant> = Vec::new();
    // coverage per amplicon
    let mut vcovs: Vec<i32> = Vec::new();
    // good amplicon set: "ampliconNumber-refallele-varallele" (populated lazily)
    let mut goodmap: HashSet<String> = HashSet::new();
    // good variants on each amplicon: amplicon_index → list of good variants
    let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
    // bad variants (-D debug column)
    let mut bad_vars: Vec<(Option<Variant>, String)> = Vec::new();

    // Iterate over positions sorted (BTreeMap already sorted)
    for (position, amplicon_regions) in positions {
        gvs.clear();
        ref_vars.clear();
        vcovs.clear();
        if !goodmap.is_empty() {
            goodmap.clear();
        }
        good_variants_on_amp.clear();

        let mut maxcov: i32 = 0;
        let mut maxaf: f64 = 0.0;

        for &(amplicon_number, ref amp_region) in &amplicon_regions {
            let vtmp = vars_per_amplicon
                .get(amplicon_number)
                .and_then(|m| m.get(&position));

            let variants_on_amplicon = vtmp.map(|v| &v.variants);
            let ref_amplicon = vtmp.and_then(|v| v.ref_variant.as_ref());

            if let Some(variants) = variants_on_amplicon {
                if !variants.is_empty() {
                    let mut good_vars: Vec<Variant> = Vec::new();
                    for tv in variants {
                        vcovs.push(tv.total_pos_cov);
                        if tv.total_pos_cov > maxcov {
                            maxcov = tv.total_pos_cov;
                        }
                        // Compute var_type
                        let vtype = tv.var_type();
                        if tv.is_good_var(ref_amplicon, &vtype, splice, config) {
                            let region_str = format!(
                                "{}:{}-{}",
                                amp_region.chr, amp_region.start.0, amp_region.end.0
                            );
                            let mut tv_clone = tv.clone();
                            tv_clone.var_type = vtype;
                            if tv_clone.frequency > maxaf {
                                maxaf = tv_clone.frequency;
                            }
                            gvs.push((tv_clone.clone(), region_str));
                            good_vars.push(tv_clone);
                        }
                    }
                    if !good_vars.is_empty() {
                        good_variants_on_amp.insert(amplicon_number, good_vars);
                    }
                } else if let Some(ref_var) = ref_amplicon {
                    vcovs.push(ref_var.total_pos_cov);
                } else {
                    vcovs.push(0);
                }
            } else if let Some(ref_var) = ref_amplicon {
                vcovs.push(ref_var.total_pos_cov);
            } else {
                vcovs.push(0);
            }

            if let Some(ref_var) = ref_amplicon {
                ref_vars.push(ref_var.clone());
            }
        }

        // Count nocov: amplicons with coverage < maxcov/50
        let nocov = vcovs
            .iter()
            .filter(|&&t| (t as f64) < (maxcov as f64) / 50.0)
            .count() as i32;

        // Sort gvs by frequency descending
        if gvs.len() > 1 {
            gvs.sort_by(|a, b| {
                b.0.frequency
                    .partial_cmp(&a.0.frequency)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }

        // Sort ref variants by totalPosCoverage descending
        if ref_vars.len() > 1 {
            ref_vars.sort_by_key(|v| std::cmp::Reverse(v.total_pos_cov));
        }

        // Build vref_list (deduplicated good variants)
        let vref_list: Vec<Variant>;
        if gvs.is_empty() {
            // Only reference
            if config.pileup {
                if !ref_vars.is_empty() {
                    vref_list = vec![ref_vars[0].clone()];
                } else {
                    // Output empty variant at this position
                    let mut vref = Variant::new();
                    vref.start_pos = position;
                    vref.end_pos = position;
                    vref.var_type = String::new();
                    emit(AmpliconVariant {
                        variant: vref,
                        gdamp: 0,
                        tlamp: amplicon_regions.len() as i32,
                        ncamp: nocov,
                        ampflag: false,
                        region: format!("{}:{}-{}", region.chr, position, position),
                        good_variants: Vec::new(),
                        bad_variants: Vec::new(),
                    })?;
                    continue;
                }
            } else {
                continue;
            }
        } else {
            vref_list = fill_vref_list(&gvs);
        }

        // Detect amplicon bias
        let flag = is_amp_bias_flag(&good_variants_on_amp);

        let total_amplicons = amplicon_regions.len() as i32;

        for vref in &vref_list {
            let mut current_flag = flag;

            // If flag is set, re-check by looking at top variant's description across all amplicons
            let good_variants = if current_flag {
                if let Some(top_gv) = gvs.first() {
                    let gdnt = &top_gv.0.description;
                    let mut gcnt: Vec<(Variant, String)> = Vec::new();
                    for &(amp_num, ref amp_reg) in &amplicon_regions {
                        let vtmp = vars_per_amplicon
                            .get(amp_num)
                            .and_then(|m| m.get(&position));
                        if let Some(vtmp) = vtmp {
                            if let Some(variant) = vtmp.var_desc_to_variant.get(gdnt) {
                                if variant.is_good_var(
                                    vtmp.ref_variant.as_ref(),
                                    "",
                                    splice,
                                    config,
                                ) {
                                    let region_str = format!(
                                        "{}:{}-{}",
                                        amp_reg.chr, amp_reg.start.0, amp_reg.end.0
                                    );
                                    gcnt.push((variant.clone(), region_str));
                                }
                            }
                        }
                    }
                    if gcnt.len() == gvs.len() {
                        current_flag = false;
                    }
                    // Sort gcnt by frequency descending
                    gcnt.sort_by(|a, b| {
                        b.0.frequency
                            .partial_cmp(&a.0.frequency)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    gcnt
                } else {
                    gvs.clone()
                }
            } else {
                gvs.clone()
            };

            // Count how many amplicons support this variant (initialGvscnt)
            let initial_gvscnt = count_variant_on_amplicons(vref, &good_variants_on_amp);
            let mut current_gvscnt = initial_gvscnt;

            // Compute bad variants: amplicons where this variant is not good but falls
            // within the insert region. Mirrors Java AmpliconPostProcessModule logic.
            let mut bad_count = 0;
            bad_vars.clear();
            if initial_gvscnt != total_amplicons || current_flag {
                if goodmap.is_empty() && !good_variants_on_amp.is_empty() {
                    for (&amp_idx, vars) in &good_variants_on_amp {
                        for tv in vars {
                            goodmap
                                .insert(format!("{}-{}-{}", amp_idx, tv.ref_allele, tv.var_allele));
                        }
                    }
                }
                for &(amp_num, ref amp_reg) in &amplicon_regions {
                    if config.pileup && vref.ref_allele == vref.var_allele {
                        continue;
                    }
                    if !goodmap.is_empty() {
                        let key = format!("{}-{}-{}", amp_num, vref.ref_allele, vref.var_allele);
                        if goodmap.contains(&key) {
                            continue;
                        }
                    }
                    if let Some(ref amp) = amp_reg.amplicon {
                        let insert_start = amp.insert_start.0 as i32;
                        let insert_end = amp.insert_end.0 as i32;
                        if vref.start_pos >= insert_start && vref.end_pos <= insert_end {
                            bad_count += 1;
                            // ports AmpliconPostProcessModule.java:170-179: record the
                            // amplicon's first variant, else its reference variant, else null
                            if config.debug {
                                let reg_str = format!(
                                    "{}:{}-{}",
                                    amp_reg.chr, amp_reg.start.0, amp_reg.end.0
                                );
                                let vtmp = vars_per_amplicon
                                    .get(amp_num)
                                    .and_then(|m| m.get(&position));
                                let bad_var = match vtmp {
                                    Some(v) if !v.variants.is_empty() => {
                                        Some(v.variants[0].clone())
                                    }
                                    Some(v) => v.ref_variant.clone(),
                                    None => None,
                                };
                                bad_vars.push((bad_var, reg_str));
                            }
                        } else if (vref.start_pos < insert_end && insert_end < vref.end_pos)
                            || (vref.start_pos < insert_start && insert_start < vref.end_pos)
                        {
                            // Variant overlaps the amplicon primer
                            if current_gvscnt > 1 {
                                current_gvscnt -= 1;
                            }
                        }
                    }
                }
            }
            let gdamp = current_gvscnt;
            let tlamp = current_gvscnt + bad_count;

            // Apply complex variant adjustment
            let mut out_variant = vref.clone();
            out_variant.var_type = out_variant.var_type();
            if out_variant.var_type == "Complex" {
                out_variant.adj_complex();
            }

            let region_str = if !good_variants.is_empty() {
                good_variants[0].1.clone()
            } else {
                format!("{}:{}-{}", region.chr, position, position)
            };

            // Debug lists are only needed for the -D debug column (zero-cost when off).
            let (debug_good, debug_bad) = if config.debug {
                (
                    good_variants
                        .iter()
                        .map(|(v, r)| (Some(v.clone()), r.clone()))
                        .collect(),
                    std::mem::take(&mut bad_vars),
                )
            } else {
                (Vec::new(), Vec::new())
            };

            emit(AmpliconVariant {
                variant: out_variant,
                gdamp,
                tlamp,
                ncamp: nocov,
                ampflag: current_flag,
                region: region_str,
                good_variants: debug_good,
                bad_variants: debug_bad,
            })?;
        }
    }

    Ok(())
}

/// Global top-1-per-raw-start selection for VCF mode.
///
/// Split out of `amplicon_post_process` so the -p VCF path can reuse it after
/// iterator-fed collection. Bucketing is inherently global (tovars shifts such
/// as the deletion `-1` anchor can move a row's start across positions), so
/// only the intermediate (-p table) path streams rows straight to the writer.
pub fn amplicon_select_vcf_rows(results: Vec<AmpliconVariant>) -> Vec<AmpliconVariant> {
    // VCF mode: same var2vcf_valid.pl selection as simple mode (:41,129-137) —
    // bucket rows by the RAW .var start column, stable af-desc sort per bucket,
    // emit only the top row ($ALL = 1 without -A). See postprocess/simple.rs.
    let mut buckets: BTreeMap<i32, Vec<AmpliconVariant>> = BTreeMap::new();
    for av in results {
        buckets.entry(av.variant.start_pos).or_default().push(av);
    }
    let mut selected = Vec::with_capacity(buckets.len());
    for (_, mut rows) in buckets {
        rows.sort_by(|a, b| {
            b.variant
                .frequency
                .partial_cmp(&a.variant.frequency)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        selected.extend(rows.into_iter().take(1));
    }
    selected
}

/// Fill vref_list from good variants, deduplicating by ref_allele + var_allele.
/// The variant with the highest frequency is kept (gvs is already sorted by frequency desc).
///
/// Ports `AmpliconPostProcessModule.fillVrefList()`.
fn fill_vref_list(gvs: &[(Variant, String)]) -> Vec<Variant> {
    let mut vref_list: Vec<Variant> = Vec::new();
    for (variant, _region_str) in gvs {
        let already_added = vref_list
            .iter()
            .any(|v| v.var_allele == variant.var_allele && v.ref_allele == variant.ref_allele);
        if !already_added {
            vref_list.push(variant.clone());
        }
    }
    vref_list
}

/// Count amplicons where a good variant with matching ref/var allele appears.
///
/// Ports `AmpliconPostProcessModule.countVariantOnAmplicons()`.
fn count_variant_on_amplicons(
    vref: &Variant,
    good_variants_on_amp: &BTreeMap<usize, Vec<Variant>>,
) -> i32 {
    let mut gvscnt: i32 = 0;
    for variants in good_variants_on_amp.values() {
        for variant in variants {
            if variant.ref_allele == vref.ref_allele && variant.var_allele == vref.var_allele {
                gvscnt += 1;
            }
        }
    }
    gvscnt
}

/// Determine if different amplicons contain different variants → set AMPBIAS flag.
///
/// Checks if each amplicon at a position contains identical lists of variants.
/// If some variants are absent between amplicons, or differ by variant description
/// string, returns true.
///
/// Ports `AmpliconPostProcessModule.isAmpBiasFlag()`.
fn is_amp_bias_flag(good_variants_on_amp: &BTreeMap<usize, Vec<Variant>>) -> bool {
    if good_variants_on_amp.is_empty() {
        return false;
    }

    let amplicon_list: Vec<usize> = good_variants_on_amp.keys().copied().collect();
    // BTreeMap keys are already sorted

    for i in 0..amplicon_list.len().saturating_sub(1) {
        let current_amplicon = amplicon_list[i];
        let next_amplicon = amplicon_list[i + 1];

        let current_variants = &good_variants_on_amp[&current_amplicon];
        let next_variants = &good_variants_on_amp[&next_amplicon];

        if current_variants.len() != next_variants.len() {
            return true;
        }

        // Sort both by totalPosCoverage descending for comparison
        let mut current_sorted: Vec<&Variant> = current_variants.iter().collect();
        let mut next_sorted: Vec<&Variant> = next_variants.iter().collect();
        current_sorted.sort_by_key(|v| std::cmp::Reverse(v.total_pos_cov));
        next_sorted.sort_by_key(|v| std::cmp::Reverse(v.total_pos_cov));

        for (var1, var2) in current_sorted.iter().zip(next_sorted.iter()) {
            if var1.description != var2.description {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::regions::AmpliconInfo;
    use crate::types::Pos1;

    fn make_test_config() -> Config {
        let mut config = Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;
        config.mfreq = 0.25;
        config.nmfreq = 0.1;
        config.pileup = false;
        // These tests assert the full per-position table (Tier-1 path); VCF-mode
        // raw-start top-AF selection would collapse same-start rows.
        config.emit_intermediate = true;
        config
    }

    fn make_good_variant(pos: i32, ref_a: &str, var_a: &str, freq: f64) -> Variant {
        let mut v = Variant::new();
        v.ref_allele = ref_a.to_string();
        v.var_allele = var_a.to_string();
        v.description = var_a.to_string();
        v.frequency = freq;
        v.hicnt = 5;
        v.mean_position = 10.0;
        v.mean_quality = 30.0;
        v.hq_lq_ratio = 2.0;
        v.mean_mapping_quality = 30.0;
        v.start_pos = pos;
        v.end_pos = pos;
        v.total_pos_cov = 100;
        v.pos_coverage = 10;
        v
    }

    fn make_amp_region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start: Pos1(start),
            end: Pos1(end),
            gene: "TEST".to_string(),
            amplicon: Some(AmpliconInfo {
                insert_start: Pos1(start),
                insert_end: Pos1(end),
            }),
        }
    }

    #[test]
    fn test_single_amplicon_single_variant() {
        let config = make_test_config();
        let region = make_amp_region("chr1", 100, 200);
        let splice = HashSet::new();

        // One amplicon with one good variant at position 150
        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut vars_150 = Vars::new();
        vars_150
            .variants
            .push(make_good_variant(150, "A", "G", 0.10));
        vars_amp0.insert(150, vars_150);

        let vars_per_amplicon = vec![vars_amp0];

        let mut amplicons_on_positions = BTreeMap::new();
        amplicons_on_positions.insert(150, vec![(0, make_amp_region("chr1", 100, 200))]);

        let results = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .expect("batch wrapper collects infallibly");

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].variant.start_pos, 150);
        assert_eq!(results[0].gdamp, 1);
        assert_eq!(results[0].tlamp, 1);
        assert_eq!(results[0].ncamp, 0);
        assert!(!results[0].ampflag);
    }

    #[test]
    fn test_two_amplicons_same_variant_no_bias() {
        let config = make_test_config();
        let region = make_amp_region("chr1", 100, 300);
        let splice = HashSet::new();

        // Two amplicons, both with the same variant at position 150
        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut vars_150_a = Vars::new();
        vars_150_a
            .variants
            .push(make_good_variant(150, "A", "G", 0.10));
        vars_amp0.insert(150, vars_150_a);

        let mut vars_amp1 = crate::javamap::JavaHashMap::new();
        let mut vars_150_b = Vars::new();
        vars_150_b
            .variants
            .push(make_good_variant(150, "A", "G", 0.12));
        vars_amp1.insert(150, vars_150_b);

        let vars_per_amplicon = vec![vars_amp0, vars_amp1];

        let mut amplicons_on_positions = BTreeMap::new();
        amplicons_on_positions.insert(
            150,
            vec![
                (0, make_amp_region("chr1", 100, 200)),
                (1, make_amp_region("chr1", 150, 300)),
            ],
        );

        let results = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .expect("batch wrapper collects infallibly");

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].gdamp, 2);
        assert_eq!(results[0].tlamp, 2);
        assert!(!results[0].ampflag);
    }

    #[test]
    fn test_two_amplicons_different_variants_bias() {
        let config = make_test_config();
        let region = make_amp_region("chr1", 100, 300);
        let splice = HashSet::new();

        // Two amplicons with different variants at position 150
        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut vars_150_a = Vars::new();
        vars_150_a
            .variants
            .push(make_good_variant(150, "A", "G", 0.10));
        vars_amp0.insert(150, vars_150_a);

        let mut vars_amp1 = crate::javamap::JavaHashMap::new();
        let mut vars_150_b = Vars::new();
        vars_150_b
            .variants
            .push(make_good_variant(150, "A", "T", 0.12));
        vars_amp1.insert(150, vars_150_b);

        let vars_per_amplicon = vec![vars_amp0, vars_amp1];

        let mut amplicons_on_positions = BTreeMap::new();
        amplicons_on_positions.insert(
            150,
            vec![
                (0, make_amp_region("chr1", 100, 200)),
                (1, make_amp_region("chr1", 150, 300)),
            ],
        );

        let results = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .expect("batch wrapper collects infallibly");

        // Both variants should be output, with ampflag = true
        assert_eq!(results.len(), 2);
        // First variant (highest freq) is A>T
        assert_eq!(results[0].variant.var_allele, "T");
        assert!(results[0].ampflag);
        assert_eq!(results[0].gdamp, 1);
        // Second variant is A>G
        assert_eq!(results[1].variant.var_allele, "G");
        assert!(results[1].ampflag);
        assert_eq!(results[1].gdamp, 1);
    }

    #[test]
    fn test_nocov_detection() {
        let config = make_test_config();
        let region = make_amp_region("chr1", 100, 300);
        let splice = HashSet::new();

        // Amplicon 0: high coverage variant
        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut vars_150_a = Vars::new();
        let mut v = make_good_variant(150, "A", "G", 0.10);
        v.total_pos_cov = 1000;
        vars_150_a.variants.push(v);
        vars_amp0.insert(150, vars_150_a);

        // Amplicon 1: very low coverage (below maxcov/50 = 1000/50 = 20)
        let mut vars_amp1 = crate::javamap::JavaHashMap::new();
        let mut vars_150_b = Vars::new();
        let mut v2 = make_good_variant(150, "A", "G", 0.10);
        v2.total_pos_cov = 10; // < 20, so nocov
        vars_150_b.variants.push(v2);
        vars_amp1.insert(150, vars_150_b);

        let vars_per_amplicon = vec![vars_amp0, vars_amp1];

        let mut amplicons_on_positions = BTreeMap::new();
        amplicons_on_positions.insert(
            150,
            vec![
                (0, make_amp_region("chr1", 100, 200)),
                (1, make_amp_region("chr1", 150, 300)),
            ],
        );

        let results = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .expect("batch wrapper collects infallibly");

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].ncamp, 1); // one amplicon has low coverage
    }

    #[test]
    fn test_no_good_variants_no_pileup() {
        let config = make_test_config();
        let region = make_amp_region("chr1", 100, 200);
        let splice = HashSet::new();

        // One amplicon with a variant that fails isGoodVar (low frequency)
        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut vars_150 = Vars::new();
        let mut v = make_good_variant(150, "A", "G", 0.001); // below freq threshold
        v.frequency = 0.001;
        vars_150.variants.push(v);
        vars_amp0.insert(150, vars_150);

        let vars_per_amplicon = vec![vars_amp0];

        let mut amplicons_on_positions = BTreeMap::new();
        amplicons_on_positions.insert(150, vec![(0, make_amp_region("chr1", 100, 200))]);

        let results = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .expect("batch wrapper collects infallibly");

        // No good variants and pileup=false → no output
        assert!(results.is_empty());
    }

    #[test]
    fn test_no_good_variants_with_pileup() {
        let mut config = make_test_config();
        config.pileup = true;
        let region = make_amp_region("chr1", 100, 200);
        let splice = HashSet::new();

        // One amplicon with a reference variant but no good variants
        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut vars_150 = Vars::new();
        let mut ref_v = Variant::new();
        ref_v.ref_allele = "A".to_string();
        ref_v.var_allele = "A".to_string();
        ref_v.start_pos = 150;
        ref_v.end_pos = 150;
        ref_v.total_pos_cov = 50;
        vars_150.ref_variant = Some(ref_v);
        vars_amp0.insert(150, vars_150);

        let vars_per_amplicon = vec![vars_amp0];

        let mut amplicons_on_positions = BTreeMap::new();
        amplicons_on_positions.insert(150, vec![(0, make_amp_region("chr1", 100, 200))]);

        let results = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .expect("batch wrapper collects infallibly");

        // Pileup mode: reference variant is output
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].variant.start_pos, 150);
    }

    #[test]
    fn test_is_amp_bias_flag_empty() {
        let good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        assert!(!is_amp_bias_flag(&good_variants_on_amp));
    }

    #[test]
    fn test_is_amp_bias_flag_same_variants() {
        let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        good_variants_on_amp.insert(0, vec![make_good_variant(150, "A", "G", 0.10)]);
        good_variants_on_amp.insert(1, vec![make_good_variant(150, "A", "G", 0.12)]);
        assert!(!is_amp_bias_flag(&good_variants_on_amp));
    }

    #[test]
    fn test_is_amp_bias_flag_different_variants() {
        let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        good_variants_on_amp.insert(0, vec![make_good_variant(150, "A", "G", 0.10)]);
        good_variants_on_amp.insert(1, vec![make_good_variant(150, "A", "T", 0.12)]);
        assert!(is_amp_bias_flag(&good_variants_on_amp));
    }

    #[test]
    fn test_is_amp_bias_flag_different_count() {
        let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        good_variants_on_amp.insert(
            0,
            vec![
                make_good_variant(150, "A", "G", 0.10),
                make_good_variant(150, "A", "T", 0.05),
            ],
        );
        good_variants_on_amp.insert(1, vec![make_good_variant(150, "A", "G", 0.12)]);
        assert!(is_amp_bias_flag(&good_variants_on_amp));
    }

    #[test]
    fn test_fill_vref_list_dedup() {
        let gvs = vec![
            (
                make_good_variant(150, "A", "G", 0.12),
                "chr1:100-200".to_string(),
            ),
            (
                make_good_variant(150, "A", "G", 0.10),
                "chr1:150-300".to_string(),
            ),
            (
                make_good_variant(150, "A", "T", 0.08),
                "chr1:100-200".to_string(),
            ),
        ];
        let vref_list = fill_vref_list(&gvs);
        // A>G appears twice but should be deduplicated (keep first = highest freq)
        assert_eq!(vref_list.len(), 2);
        assert_eq!(vref_list[0].var_allele, "G");
        assert_eq!(vref_list[1].var_allele, "T");
    }

    #[test]
    fn test_count_variant_on_amplicons() {
        let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        good_variants_on_amp.insert(0, vec![make_good_variant(150, "A", "G", 0.10)]);
        good_variants_on_amp.insert(1, vec![make_good_variant(150, "A", "G", 0.12)]);
        good_variants_on_amp.insert(2, vec![make_good_variant(150, "A", "T", 0.08)]);

        let vref = make_good_variant(150, "A", "G", 0.12);
        assert_eq!(count_variant_on_amplicons(&vref, &good_variants_on_amp), 2);

        let vref_t = make_good_variant(150, "A", "T", 0.08);
        assert_eq!(
            count_variant_on_amplicons(&vref_t, &good_variants_on_amp),
            1
        );
    }

    #[test]
    fn pileup_position_stream_matches_expansion() {
        use crate::regions::AmpliconInfo;
        use crate::types::Pos1;

        fn amp_region(chr: &str, start: u32, end: u32, is: u32, ie: u32) -> Region {
            Region {
                chr: chr.to_string(),
                start: Pos1(start),
                end: Pos1(end),
                gene: String::new(),
                amplicon: Some(AmpliconInfo {
                    insert_start: Pos1(is),
                    insert_end: Pos1(ie),
                }),
            }
        }

        // Overlapping inserts plus a disjoint interval (its gap is never
        // visited) plus one entry without amplicon info (covers nothing).
        let segment = vec![
            amp_region("chr1", 90, 210, 100, 200),
            amp_region("chr1", 140, 260, 150, 250),
            amp_region("chr1", 500, 510, 500, 510),
            Region {
                chr: "chr1".to_string(),
                start: Pos1(300),
                end: Pos1(400),
                gene: String::new(),
                amplicon: None,
            },
        ];

        // Old logic: per-base expansion (AmpliconMode.java:57-60).
        let mut expanded: BTreeMap<i32, Vec<(usize, Region)>> = BTreeMap::new();
        for (n, region) in segment.iter().enumerate() {
            if let Some(ref amp) = region.amplicon {
                for p in amp.insert_start.0 as i32..=amp.insert_end.0 as i32 {
                    expanded.entry(p).or_default().push((n, region.clone()));
                }
            }
        }
        // The stream yields the same pairs, in the same ascending order.
        let streamed: Vec<(i32, Vec<(usize, Region)>)> =
            amplicon_position_stream(&segment).collect();
        assert_eq!(streamed.len(), expanded.len());
        for ((got_p, got_cov), (want_p, want_cov)) in streamed.iter().zip(expanded.iter()) {
            assert_eq!(got_p, want_p);
            assert_eq!(got_cov, want_cov);
        }
    }

    #[test]
    fn sweep_matches_old_scan_incl_tie_starts() {
        use crate::regions::AmpliconInfo;
        use crate::types::Pos1;

        fn amp_region(chr: &str, start: u32, end: u32, is: u32, ie: u32) -> Region {
            Region {
                chr: chr.to_string(),
                start: Pos1(start),
                end: Pos1(end),
                gene: String::new(),
                amplicon: Some(AmpliconInfo {
                    insert_start: Pos1(is),
                    insert_end: Pos1(ie),
                }),
            }
        }

        // Tie insert starts (amps 0 and 3 share start 150) in non-sorted
        // segment order — amp 1 starts earlier despite its later index — plus
        // a nested interval, a disjoint one, and a non-amplicon entry.
        let segment = vec![
            amp_region("chr1", 140, 260, 150, 250),
            amp_region("chr1", 90, 210, 100, 200),
            amp_region("chr1", 500, 510, 500, 510),
            amp_region("chr1", 150, 160, 150, 155),
            Region {
                chr: "chr1".to_string(),
                start: Pos1(300),
                end: Pos1(400),
                gene: String::new(),
                amplicon: None,
            },
            amp_region("chr1", 180, 190, 180, 185),
        ];

        // Old logic: per-base expansion (AmpliconMode.java:57-60).
        let mut expanded: BTreeMap<i32, Vec<(usize, Region)>> = BTreeMap::new();
        for (n, region) in segment.iter().enumerate() {
            if let Some(ref amp) = region.amplicon {
                for p in amp.insert_start.0 as i32..=amp.insert_end.0 as i32 {
                    expanded.entry(p).or_default().push((n, region.clone()));
                }
            }
        }
        // The sweep stream yields the same pairs, pair-for-pair: same
        // positions ascending, same segment-ordered covering lists — tie
        // starts keep segment order via the stable by-start sort.
        let streamed: Vec<(i32, Vec<(usize, Region)>)> =
            amplicon_position_stream(&segment).collect();
        assert_eq!(streamed.len(), expanded.len());
        for ((got_p, got_cov), (want_p, want_cov)) in streamed.iter().zip(expanded.iter()) {
            assert_eq!(got_p, want_p);
            assert_eq!(got_cov, want_cov);
        }
    }

    #[test]
    fn sweep_covering_at_matches_old_scan_sparse() {
        use crate::regions::AmpliconInfo;
        use crate::types::Pos1;

        fn amp_region(chr: &str, start: u32, end: u32, is: u32, ie: u32) -> Region {
            Region {
                chr: chr.to_string(),
                start: Pos1(start),
                end: Pos1(end),
                gene: String::new(),
                amplicon: Some(AmpliconInfo {
                    insert_start: Pos1(is),
                    insert_end: Pos1(ie),
                }),
            }
        }

        let segment = vec![
            amp_region("chr1", 140, 260, 150, 250),
            amp_region("chr1", 90, 210, 100, 200),
            amp_region("chr1", 150, 160, 150, 155),
        ];

        // Sparse ascending positions (the non-`-p` run_amplicon path feeds a
        // BTreeSet in this order): boundaries, tie starts, gaps, outside.
        let mut sweep = AmpliconCoverSweep::new(&segment);
        for p in [99, 100, 149, 150, 151, 155, 156, 180, 200, 201, 250, 251] {
            let got = sweep.covering_at(p);
            let want = crate::amplicons_covering_position(&segment, p);
            assert_eq!(got, want, "covering list differs at {p}");
        }
    }

    #[test]
    fn sweep_dense_chain_bounded() {
        use crate::regions::AmpliconInfo;
        use crate::types::Pos1;

        // Dense synthetic chain: 3000 inserts x 200 bp @100 spacing. The old
        // per-position scan costs P x A ~ 270M predicates here (~13 s in debug);
        // the sweep is output-proportional + O((P + A) log A) (~0.1 s).
        // The 3 s wall bound passes with wide headroom on the sweep and fails
        // on the old approach.
        let amps = 3000usize;
        let segment: Vec<Region> = (0..amps)
            .map(|i| {
                let s = (i as u32) * 100;
                Region {
                    chr: "chr1".to_string(),
                    start: Pos1(s),
                    end: Pos1(s + 199),
                    gene: String::new(),
                    amplicon: Some(AmpliconInfo {
                        insert_start: Pos1(s),
                        insert_end: Pos1(s + 199),
                    }),
                }
            })
            .collect();

        let t0 = std::time::Instant::now();
        let mut npos = 0u64;
        let mut nmem = 0u64;
        for (_, cov) in amplicon_position_stream(&segment) {
            npos += 1;
            nmem += cov.len() as u64;
        }
        let elapsed = t0.elapsed();
        assert_eq!(npos, (amps as u64 - 1) * 100 + 200);
        assert_eq!(nmem, amps as u64 * 200);
        assert!(
            elapsed < std::time::Duration::from_secs(3),
            "dense-chain sweep took too long: {elapsed:?}",
        );
    }

    #[test]
    fn pileup_stream_matches_batch() {
        // Two overlapping amplicons; positions 150-152 cover the good-variant,
        // reference-only, and empty (null row under -p) paths.
        let mut config = make_test_config();
        config.pileup = true;
        let region = make_amp_region("chr1", 100, 300);
        let splice = HashSet::new();

        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut v150 = Vars::new();
        v150.variants.push(make_good_variant(150, "A", "G", 0.10));
        vars_amp0.insert(150, v150);
        let mut v151 = Vars::new();
        let mut ref_v = Variant::new();
        ref_v.ref_allele = "C".to_string();
        ref_v.var_allele = "C".to_string();
        ref_v.start_pos = 151;
        ref_v.end_pos = 151;
        ref_v.total_pos_cov = 40;
        v151.ref_variant = Some(ref_v);
        vars_amp0.insert(151, v151);

        let mut vars_amp1 = crate::javamap::JavaHashMap::new();
        let mut v150b = Vars::new();
        v150b.variants.push(make_good_variant(150, "A", "G", 0.12));
        vars_amp1.insert(150, v150b);

        let vars_per_amplicon = vec![vars_amp0, vars_amp1];
        let segment = vec![
            make_amp_region("chr1", 100, 200),
            make_amp_region("chr1", 150, 300),
        ];

        // Old-style per-base map, restricted to the exercised positions.
        let mut amplicons_on_positions = BTreeMap::new();
        for p in 150..=152 {
            amplicons_on_positions
                .insert(p, vec![(0, segment[0].clone()), (1, segment[1].clone())]);
        }
        let batch = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .expect("batch wrapper collects infallibly");

        let mut streamed = Vec::new();
        {
            let mut collect = |av: AmpliconVariant| -> Result<(), crate::types::VardictError> {
                streamed.push(av);
                Ok(())
            };
            amplicon_post_process_streaming(
                &config,
                &region,
                &vars_per_amplicon,
                amplicon_position_stream(&segment).filter(|(p, _)| (150..=152).contains(p)),
                &splice,
                &mut collect,
            )
            .expect("streaming is infallible here");
        }

        // Variant row + reference row + null row, row-for-row identical.
        assert_eq!(batch.len(), 3);
        assert_eq!(streamed.len(), batch.len());
        assert_eq!(format!("{streamed:?}"), format!("{batch:?}"));
    }

    #[test]
    fn vcf_select_keeps_top_af_per_start() {
        fn row(pos: i32, af: f64, tag: &str) -> AmpliconVariant {
            let mut v = make_good_variant(pos, "A", "G", af);
            v.start_pos = pos;
            v.end_pos = pos;
            AmpliconVariant {
                variant: v,
                gdamp: 1,
                tlamp: 1,
                ncamp: 0,
                ampflag: false,
                region: tag.to_string(),
                good_variants: Vec::new(),
                bad_variants: Vec::new(),
            }
        }

        // Same start collapses to the highest AF; ties keep the first row
        // (stable sort, matching the batch path); starts emerge ascending.
        let selected = amplicon_select_vcf_rows(vec![
            row(150, 0.05, "low"),
            row(151, 0.30, "other"),
            row(150, 0.20, "high"),
            row(150, 0.10, "mid"),
            row(152, 0.10, "first"),
            row(152, 0.10, "second"),
        ]);
        assert_eq!(selected.len(), 3);
        assert_eq!(selected[0].variant.start_pos, 150);
        assert_eq!(selected[0].variant.frequency, 0.20);
        assert_eq!(selected[0].region, "high");
        assert_eq!(selected[1].variant.start_pos, 151);
        assert_eq!(selected[2].variant.start_pos, 152);
        assert_eq!(selected[2].region, "first");
    }

    #[test]
    fn pileup_stream_huge_interval_bounded() {
        // -p over a 2 Mbp insert with zero coverage everywhere: every base
        // yields a null row. Asserts full row-per-base coverage plus a wall
        // bound; the allocator-observed cap lives in
        // tests/amplicon_stream_mem.rs, which runs the same stream in an
        // isolated child process (process-wide VmHWM is polluted by other
        // test threads in-suite — measured 315 MB of unrelated-test noise in
        // one full-suite run vs <64 MB isolated for this exact stream).
        let mut config = make_test_config();
        config.pileup = true;
        let region = make_amp_region("chr1", 1, 2_000_000);
        let segment = vec![region.clone()];
        let vars_per_amplicon = vec![crate::javamap::JavaHashMap::new()];

        let t0 = std::time::Instant::now();
        let mut rows = 0u64;
        {
            let mut count = |_av: AmpliconVariant| -> Result<(), crate::types::VardictError> {
                rows += 1;
                Ok(())
            };
            amplicon_post_process_streaming(
                &config,
                &region,
                &vars_per_amplicon,
                amplicon_position_stream(&segment),
                &HashSet::new(),
                &mut count,
            )
            .expect("streaming is infallible here");
        }
        let elapsed = t0.elapsed();
        assert_eq!(rows, 2_000_000, "every insert base must yield its -p row");
        assert!(
            elapsed < std::time::Duration::from_secs(60),
            "2 Mbp -p stream took too long: {elapsed:?}",
        );
    }
}
