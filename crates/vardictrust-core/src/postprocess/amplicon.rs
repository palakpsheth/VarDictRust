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
use std::borrow::Borrow;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

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
pub fn amplicon_post_process<R: std::borrow::Borrow<Region>>(
    config: &Config,
    region: &Region,
    vars_per_amplicon: &[crate::javamap::JavaHashMap<i32, Vars>],
    amplicons_on_positions: &BTreeMap<i32, Vec<(usize, R)>>,
    splice: &rustc_hash::FxHashSet<String>,
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
                .map(|(&position, regions)| (position, regions.as_slice())),
            splice,
            &mut collect,
        )?;
    }

    if config.compat.strict_region_bounds {
        amplicon_retain_strict_bounds(
            &mut results,
            region,
            amplicons_on_positions,
            config.compat.canonical_sv_vcf_anchoring,
        );
    }

    if !config.emit_intermediate && !config.vcf_all_alleles {
        results = amplicon_select_vcf_rows_compat(results, config.compat.exact_binary_half_even);
    }

    Ok(results)
}

#[inline]
pub(crate) fn region_target_range_i32(region: &Region) -> (i32, i32) {
    match region.amplicon.as_ref() {
        Some(amp) => (
            i32::try_from(amp.insert_start.0.max(1)).unwrap_or(i32::MAX),
            i32::try_from(amp.insert_end.0).unwrap_or(i32::MAX),
        ),
        None => (
            i32::try_from(region.start.0.max(1)).unwrap_or(i32::MAX),
            i32::try_from(region.end.0).unwrap_or(i32::MAX),
        ),
    }
}

/// Filter `results` in-place to the merged 1-based closed amplicon insert bounds of
/// `amplicons_on_positions` (falling back to `fallback_region` when empty,
/// `AstraZeneca-NGS/VarDictJava#369`, `AR-4`), retaining 1-bp-upstream-anchored
/// symbolic SVs when `canonical_sv_anchoring` is active.
pub fn amplicon_retain_strict_bounds<R: std::borrow::Borrow<Region>>(
    results: &mut Vec<AmpliconVariant>,
    fallback_region: &Region,
    amplicons_on_positions: &BTreeMap<i32, Vec<(usize, R)>>,
    canonical_sv_anchoring: bool,
) {
    let mut seen_amp_idxs = rustc_hash::FxHashSet::default();
    let mut ranges = Vec::new();
    for covering in amplicons_on_positions.values() {
        for (amp_idx, reg) in covering {
            if seen_amp_idxs.insert(*amp_idx) {
                ranges.push(region_target_range_i32(reg.borrow()));
            }
        }
    }
    if ranges.is_empty() {
        ranges.push(region_target_range_i32(fallback_region));
    }
    let merged = crate::merge_i32_intervals(ranges);
    results.retain(|av| {
        crate::variant_within_strict_bounds(&merged, &av.variant, canonical_sv_anchoring)
    });
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
    cached_covering: Option<std::sync::Arc<[(usize, &'a Region)]>>,
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
            cached_covering: None,
            #[cfg(debug_assertions)]
            last_pos: None,
        }
    }

    fn advance_to(&mut self, position: i32) {
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
        let mut changed = false;
        while let Some(idx) = self.by_start.get(self.next_start) {
            let idx = *idx;
            let start = match self.starts.get(idx) {
                Some(start) => *start,
                None => break,
            };
            if start > position {
                break;
            }
            if self.active.insert(idx) {
                changed = true;
            }
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
                if self.active.remove(&idx) {
                    changed = true;
                }
            }
        }
        if changed {
            self.cached_covering = None;
        }
    }

    /// Covering amplicon references at `pos`, in segment order, without cloning
    /// `Region` per base. Positions MUST arrive in ascending (non-decreasing) order.
    #[cfg(test)]
    pub(crate) fn covering_refs_at(&mut self, pos: i32) -> Vec<(usize, &'a Region)> {
        self.advance_to(pos);
        let mut covering = Vec::with_capacity(self.active.len());
        for idx in &self.active {
            if let Some(region) = self.segment.get(*idx) {
                covering.push((*idx, region));
            }
        }
        covering
    }

    /// Cached covering amplicon slice at `pos`, in segment order, sharing a single
    /// `Arc` allocation across consecutive positions whose active set is unchanged.
    /// Positions MUST arrive in ascending (non-decreasing) order.
    pub(crate) fn covering_arc_at(&mut self, pos: i32) -> std::sync::Arc<[(usize, &'a Region)]> {
        self.advance_to(pos);
        if let Some(ref cached) = self.cached_covering {
            return std::sync::Arc::clone(cached);
        }
        let mut covering = Vec::with_capacity(self.active.len());
        for idx in &self.active {
            if let Some(region) = self.segment.get(*idx) {
                covering.push((*idx, region));
            }
        }
        let arc: std::sync::Arc<[(usize, &'a Region)]> = covering.into();
        self.cached_covering = Some(std::sync::Arc::clone(&arc));
        arc
    }

    /// Covering amplicons at `position`, in segment order. Positions MUST
    /// arrive in ascending (non-decreasing) order.
    #[cfg(test)]
    pub(crate) fn covering_at(&mut self, position: i32) -> Vec<(usize, Region)> {
        self.advance_to(position);
        let mut covering = Vec::with_capacity(self.active.len());
        for idx in &self.active {
            if let Some(region) = self.segment.get(*idx) {
                covering.push((*idx, region.clone()));
            }
        }
        covering
    }
}

fn merge_insert_ranges(segment: &[Region]) -> Vec<(i32, i32)> {
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
    merged
}

/// Lazy iterator behind [`amplicon_position_ref_stream`] and [`amplicon_position_stream`]:
/// walks the merged insert ranges ascending, resolving each position's covering list
/// through one shared [`AmpliconCoverSweep`] (ascending by construction).
struct AmpliconPositionRefIter<'a> {
    sweep: AmpliconCoverSweep<'a>,
    merged: Vec<(i32, i32)>,
    range_idx: usize,
    current: i32,
}

impl<'a> Iterator for AmpliconPositionRefIter<'a> {
    type Item = (i32, std::sync::Arc<[(usize, &'a Region)]>);

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
            let covering = self.sweep.covering_arc_at(position);
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

/// Zero-clone reference variant of [`amplicon_position_stream`], borrowing
/// `&Region` directly from `segment` instead of cloning `Region` per base.
pub(crate) fn amplicon_position_ref_stream(
    segment: &[Region],
) -> impl Iterator<Item = (i32, std::sync::Arc<[(usize, &Region)]>)> + '_ {
    let merged = merge_insert_ranges(segment);
    AmpliconPositionRefIter {
        sweep: AmpliconCoverSweep::new(segment),
        merged,
        range_idx: 0,
        current: 0,
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
#[deprecated(
    since = "1.0.10",
    note = "use postprocess::amplicon::amplicon_post_process or amplicon_post_process_streaming"
)]
pub fn amplicon_position_stream(
    segment: &[Region],
) -> impl Iterator<Item = (i32, Vec<(usize, Region)>)> + '_ {
    amplicon_position_ref_stream(segment)
        .map(|(pos, refs)| (pos, refs.iter().map(|(i, r)| (*i, (*r).clone())).collect()))
}

#[inline]
fn get_or_init_region_tag<'a>(
    region_tags: &'a mut Vec<String>,
    amp_num: usize,
    amp_reg: &Region,
) -> &'a str {
    if region_tags.len() <= amp_num {
        region_tags.resize_with(amp_num.saturating_add(1), String::new);
    }
    if let Some(tag) = region_tags.get_mut(amp_num) {
        if tag.is_empty() {
            use std::fmt::Write as _;
            let _ = write!(tag, "{}:{}-{}", amp_reg.chr, amp_reg.start.0, amp_reg.end.0);
        }
        return tag.as_str();
    }
    ""
}

/// Streaming core of `AmpliconPostProcessModule.process()`: one iteration per
/// position, each output row handed to `emit` immediately (the Java module
/// likewise prints per position via `variantPrinter`, without buffering rows).
///
/// `positions` MUST arrive in ascending order with segment-ordered covering
/// lists — both the batch `BTreeMap` and `amplicon_position_stream` satisfy
/// this — so row order matches the batch path exactly.
pub fn amplicon_post_process_streaming<R: std::borrow::Borrow<Region>, C: AsRef<[(usize, R)]>>(
    config: &Config,
    region: &Region,
    vars_per_amplicon: &[crate::javamap::JavaHashMap<i32, Vars>],
    positions: impl IntoIterator<Item = (i32, C)>,
    splice: &rustc_hash::FxHashSet<String>,
    emit: &mut dyn FnMut(AmpliconVariant) -> Result<(), crate::types::VardictError>,
) -> Result<(), crate::types::VardictError> {
    amplicon_post_process_streaming_with_ref(
        config,
        region,
        vars_per_amplicon,
        positions,
        splice,
        None,
        emit,
    )
}

pub(crate) fn amplicon_post_process_streaming_with_ref<
    R: std::borrow::Borrow<Region>,
    C: AsRef<[(usize, R)]>,
>(
    config: &Config,
    region: &Region,
    vars_per_amplicon: &[crate::javamap::JavaHashMap<i32, Vars>],
    positions: impl IntoIterator<Item = (i32, C)>,
    splice: &rustc_hash::FxHashSet<String>,
    ref_info: Option<(&crate::reference::Reference, i32)>,
    emit: &mut dyn FnMut(AmpliconVariant) -> Result<(), crate::types::VardictError>,
) -> Result<(), crate::types::VardictError> {
    let skip_eq_non_pileup = !config.pileup && config.compat.sanitize_iupac_alt_alleles;
    let refresh_vt = config.compat.recompute_vartype_after_adj_complex;
    let do_left_align = crate::variations::canonical_left_align_enabled(config);
    let indel_only_shift3 = config.compat.indel_only_shift3;
    let ref_ctx = ref_info.map(|(r, len)| (r, region.chr.as_str(), len));
    let dedup_non_pileup = config.compat.canonical_vcf_complex_left_align && !config.pileup;
    let mut buffered_dedup: Vec<AmpliconVariant> = Vec::new();
    let mut dedup_seen: rustc_hash::FxHashMap<(i32, String, String), usize> =
        rustc_hash::FxHashMap::default();
    // Per-amplicon cached `chr:start-end` strings formatted at most once per segment.
    let mut region_tags: Vec<String> = vec![String::new(); vars_per_amplicon.len()];
    // Scratch collections reused across positions to avoid per-base allocations.
    // good variants: (borrowed variant, amplicon_index into region_tags)
    let mut gvs: Vec<(&Variant, usize)> = Vec::new();
    // reference variants (only populated when config.pileup is active)
    let mut ref_vars: Vec<&Variant> = Vec::new();
    // coverage per amplicon
    let mut vcovs: Vec<i32> = Vec::new();
    // good variants on each amplicon: amplicon_index → list of borrowed good variants
    let mut good_variants_on_amp: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
    // deduplicated good/reference candidate list per position
    let mut vref_list: Vec<&Variant> = Vec::new();
    // re-checked top-allele good variants when AMPBIAS is active
    let mut gcnt_buf: Vec<(&Variant, usize)> = Vec::new();
    // bad variants (-D debug column)
    let mut bad_vars: Vec<(Option<Variant>, String)> = Vec::new();

    // Iterate over positions sorted (BTreeMap already sorted)
    for (position, amplicon_regions) in positions {
        let amplicon_regions = amplicon_regions.as_ref();
        gvs.clear();
        if config.pileup {
            ref_vars.clear();
        }
        vcovs.clear();
        good_variants_on_amp.clear();
        vref_list.clear();

        let mut maxcov: i32 = 0;
        let mut maxaf: f64 = 0.0;

        for (amplicon_number, amp_region_item) in amplicon_regions {
            let amplicon_number = *amplicon_number;
            let amp_region: &Region = amp_region_item.borrow();
            let vtmp = vars_per_amplicon
                .get(amplicon_number)
                .and_then(|m| m.get(&position));

            let variants_on_amplicon = vtmp.map(|v| &v.variants);
            let ref_amplicon = vtmp.and_then(|v| v.ref_variant.as_ref());

            if let Some(variants) = variants_on_amplicon {
                if !variants.is_empty() {
                    let mut good_vars: Vec<&Variant> = Vec::new();
                    let mut any_kept = false;
                    let mut skipped_eq_cov: i32 = 0;
                    for tv in variants {
                        if skip_eq_non_pileup && tv.ref_allele == tv.var_allele {
                            skipped_eq_cov = skipped_eq_cov.max(tv.total_pos_cov);
                            continue;
                        }
                        any_kept = true;
                        vcovs.push(tv.total_pos_cov);
                        if tv.total_pos_cov > maxcov {
                            maxcov = tv.total_pos_cov;
                        }
                        // Compute var_type without allocating
                        let vtype = tv.var_type_cow();
                        if tv.is_good_var(ref_amplicon, &vtype, splice, config) {
                            let _ = get_or_init_region_tag(
                                &mut region_tags,
                                amplicon_number,
                                amp_region,
                            );
                            if tv.frequency > maxaf {
                                maxaf = tv.frequency;
                            }
                            gvs.push((tv, amplicon_number));
                            good_vars.push(tv);
                        }
                    }
                    if !good_vars.is_empty() {
                        good_variants_on_amp.insert(amplicon_number, good_vars);
                    }
                    if !any_kept {
                        // Mirror pure-reference amplicon handling (`AmpliconPostProcessModule.java:94-95`,
                        // lines 521-522 below, and `aggregate_sample_amplicon_positions`): push
                        // `fallback_cov` to `vcovs` without inflating `maxcov` above variant-carrying
                        // amplicons, so a 10,000x wild-type (`0` alt) amplicon never marks a 100x
                        // variant-carrying amplicon as `NCAMP = 1` (`1 < 2 - 1` falsely clearing `AMPBIAS`).
                        let fallback_cov = ref_amplicon
                            .map_or(0, |r| r.total_pos_cov)
                            .max(skipped_eq_cov);
                        vcovs.push(fallback_cov);
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

            if config.pileup {
                if let Some(ref_var) = ref_amplicon {
                    ref_vars.push(ref_var);
                }
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

        // Sort ref variants by totalPosCoverage descending (only when pileup is active)
        if config.pileup && ref_vars.len() > 1 {
            ref_vars.sort_by_key(|v| std::cmp::Reverse(v.total_pos_cov));
        }

        // Build vref_list (deduplicated good variants)
        if gvs.is_empty() {
            // Only reference
            if config.pileup {
                if let Some(&top_ref) = ref_vars.first() {
                    vref_list.push(top_ref);
                } else {
                    if config.compat.suppress_pileup_sentinels {
                        continue;
                    }
                    // Output empty variant at this position (AmpliconPostProcessModule.java:105-106,
                    // AmpliconOutputVariant.java:25,27,54-60: new Variant() has genotype="", bias="0;0",
                    // leftseq="", rightseq="").
                    let mut vref = Variant::new();
                    vref.start_pos = position;
                    vref.end_pos = position;
                    vref.var_type = String::new();
                    vref.bias.clear();
                    vref.bias.push_str("0;0");
                    vref.null_record = true;
                    emit(AmpliconVariant {
                        variant: vref,
                        gdamp: 0,
                        tlamp: 0,
                        ncamp: 0,
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
            fill_vref_list_into(&gvs, &mut vref_list);
        }

        // Detect amplicon bias (AmpliconPostProcessModule.java:135-137)
        let mut current_flag = is_amp_bias_flag(&good_variants_on_amp);
        let mut use_gcnt = false;

        let total_amplicons = amplicon_regions.len() as i32;

        for &max_var in &vref_list {
            // If flag is set, re-check by looking at top variant's description across all amplicons
            if current_flag {
                if let Some(&(top_gv, _)) = gvs.first() {
                    let gdnt = &top_gv.description;
                    gcnt_buf.clear();
                    for (amp_num, amp_reg_item) in amplicon_regions {
                        let amp_num = *amp_num;
                        let amp_reg: &Region = amp_reg_item.borrow();
                        let vtmp = vars_per_amplicon
                            .get(amp_num)
                            .and_then(|m| m.get(&position));
                        if let Some(vtmp) = vtmp {
                            if let Some(variant) = vtmp.variant_by_desc(gdnt) {
                                if (!skip_eq_non_pileup || variant.ref_allele != variant.var_allele)
                                    && variant.is_good_var(
                                        vtmp.ref_variant.as_ref(),
                                        "",
                                        splice,
                                        config,
                                    )
                                {
                                    let _ =
                                        get_or_init_region_tag(&mut region_tags, amp_num, amp_reg);
                                    gcnt_buf.push((variant, amp_num));
                                }
                            }
                        }
                    }
                    if gcnt_buf.len() == gvs.len() {
                        current_flag = false;
                    }
                    // Sort gcnt by frequency descending
                    gcnt_buf.sort_by(|a, b| {
                        b.0.frequency
                            .partial_cmp(&a.0.frequency)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    use_gcnt = true;
                }
            }

            let good_variants: &[(&Variant, usize)] = if use_gcnt { &gcnt_buf } else { &gvs };

            // Count how many amplicons support this variant (initialGvscnt)
            let initial_gvscnt = count_variant_on_amplicons(max_var, &good_variants_on_amp);
            let mut current_gvscnt = initial_gvscnt;

            // Compute bad variants: amplicons where this variant is not good but falls
            // within the insert region. Mirrors Java AmpliconPostProcessModule logic.
            let mut bad_count: i32 = 0;
            bad_vars.clear();
            if initial_gvscnt != total_amplicons || current_flag {
                for (amp_num, amp_reg_item) in amplicon_regions {
                    let amp_num = *amp_num;
                    let amp_reg: &Region = amp_reg_item.borrow();
                    if config.pileup && max_var.ref_allele == max_var.var_allele {
                        continue;
                    }
                    let is_in_goodmap = good_variants_on_amp.get(&amp_num).is_some_and(|vars| {
                        vars.iter().any(|tv| {
                            tv.ref_allele == max_var.ref_allele
                                && tv.var_allele == max_var.var_allele
                        })
                    });
                    if is_in_goodmap {
                        continue;
                    }
                    if let Some(ref amp) = amp_reg.amplicon {
                        let insert_start = amp.insert_start.0 as i32;
                        let insert_end = amp.insert_end.0 as i32;
                        if max_var.start_pos >= insert_start && max_var.end_pos <= insert_end {
                            bad_count = bad_count.saturating_add(1);
                            // ports AmpliconPostProcessModule.java:170-179: record the
                            // amplicon's first variant, else its reference variant, else null
                            if config.debug {
                                let reg_str =
                                    get_or_init_region_tag(&mut region_tags, amp_num, amp_reg)
                                        .to_owned();
                                let vtmp = vars_per_amplicon
                                    .get(amp_num)
                                    .and_then(|m| m.get(&position));
                                let bad_var = match vtmp {
                                    Some(v) if !v.variants.is_empty() => {
                                        v.variants.first().cloned()
                                    }
                                    Some(v) => v.ref_variant.clone(),
                                    None => None,
                                };
                                bad_vars.push((bad_var, reg_str));
                            }
                        } else if (max_var.start_pos < insert_end && insert_end < max_var.end_pos)
                            || (max_var.start_pos < insert_start && insert_start < max_var.end_pos)
                        {
                            // Variant overlaps the amplicon primer
                            if current_gvscnt > 1 {
                                current_gvscnt -= 1;
                            }
                        }
                    }
                }
            }
            if current_flag && current_gvscnt < initial_gvscnt {
                current_flag = false;
            }
            let gdamp = current_gvscnt;
            let tlamp = current_gvscnt.saturating_add(bad_count);

            // Apply complex variant adjustment + Switch 35 / Switch 38 normalization
            let mut out_variant = max_var.clone();
            out_variant.var_type = out_variant.var_type();
            out_variant.normalize_postprocess(
                ref_ctx,
                do_left_align,
                refresh_vt,
                indel_only_shift3,
            );
            if skip_eq_non_pileup && out_variant.ref_allele == out_variant.var_allele {
                continue;
            }

            let region_str = if let Some(&(_, amp_idx)) = good_variants.first() {
                region_tags.get(amp_idx).cloned().unwrap_or_default()
            } else {
                format!("{}:{}-{}", region.chr, position, position)
            };

            // Debug lists are only needed for the -D debug column (zero-cost when off).
            let (debug_good, debug_bad) = if config.debug {
                (
                    good_variants
                        .iter()
                        .map(|&(v, amp_idx)| {
                            let mut gv = v.clone();
                            gv.var_type = gv.var_type();
                            (
                                Some(gv),
                                region_tags.get(amp_idx).cloned().unwrap_or_default(),
                            )
                        })
                        .collect(),
                    std::mem::take(&mut bad_vars),
                )
            } else {
                (Vec::new(), Vec::new())
            };

            let av = AmpliconVariant {
                variant: out_variant,
                gdamp,
                tlamp,
                ncamp: nocov,
                ampflag: current_flag,
                region: region_str,
                good_variants: debug_good,
                bad_variants: debug_bad,
            };
            if dedup_non_pileup {
                if av.variant.null_record {
                    buffered_dedup.push(av);
                } else {
                    let key = (
                        av.variant.start_pos,
                        av.variant.ref_allele.clone(),
                        av.variant.var_allele.clone(),
                    );
                    match dedup_seen.entry(key) {
                        std::collections::hash_map::Entry::Vacant(e) => {
                            e.insert(buffered_dedup.len());
                            buffered_dedup.push(av);
                        }
                        std::collections::hash_map::Entry::Occupied(e) => {
                            let idx = *e.get();
                            if let Some(cur) = buffered_dedup.get_mut(idx) {
                                if crate::variations::post_norm_candidate_beats(
                                    &av.variant,
                                    &cur.variant,
                                ) {
                                    *cur = av;
                                }
                            }
                        }
                    }
                }
            } else {
                emit(av)?;
            }
        }
    }

    if dedup_non_pileup {
        for av in buffered_dedup {
            emit(av)?;
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
    amplicon_select_vcf_rows_compat(results, false)
}

/// Global top-1-per-raw-start selection for VCF mode parameterized by `exact_binary_half_even`.
pub fn amplicon_select_vcf_rows_compat(
    results: Vec<AmpliconVariant>,
    exact_binary_half_even: bool,
) -> Vec<AmpliconVariant> {
    // VCF mode: same var2vcf_valid.pl selection as simple mode (:41,129-137) —
    // bucket rows by the RAW .var start column, stable af-desc sort per bucket,
    // emit only the top row ($ALL = 1 without -A). See postprocess/simple.rs.
    let mut buckets: BTreeMap<i32, AmpliconVariant> = BTreeMap::new();
    for av in results {
        let v_af =
            crate::types::round_half_even_compat(av.variant.frequency, 4, exact_binary_half_even);
        match buckets.entry(av.variant.start_pos) {
            std::collections::btree_map::Entry::Vacant(e) => {
                e.insert(av);
            }
            std::collections::btree_map::Entry::Occupied(mut e) => {
                let cur_af = crate::types::round_half_even_compat(
                    e.get().variant.frequency,
                    4,
                    exact_binary_half_even,
                );
                if v_af > cur_af {
                    e.insert(av);
                }
            }
        }
    }
    buckets.into_values().collect()
}

/// Fill vref_list from good variants, deduplicating by ref_allele + var_allele.
/// The variant with the highest frequency is kept (gvs is already sorted by frequency desc).
///
/// Ports `AmpliconPostProcessModule.fillVrefList()`.
#[cfg(test)]
fn fill_vref_list<'a, V: Borrow<Variant>, T>(gvs: &'a [(V, T)]) -> Vec<&'a Variant> {
    let mut vref_list: Vec<&'a Variant> = Vec::new();
    for (variant, _region_str) in gvs {
        let variant = variant.borrow();
        let already_added = vref_list
            .iter()
            .any(|v| v.var_allele == variant.var_allele && v.ref_allele == variant.ref_allele);
        if !already_added {
            vref_list.push(variant);
        }
    }
    vref_list
}

fn fill_vref_list_into<'v, T>(gvs: &[(&'v Variant, T)], vref_list: &mut Vec<&'v Variant>) {
    vref_list.clear();
    for &(variant, _) in gvs {
        let already_added = vref_list
            .iter()
            .any(|v| v.var_allele == variant.var_allele && v.ref_allele == variant.ref_allele);
        if !already_added {
            vref_list.push(variant);
        }
    }
}

/// Count amplicons where a good variant with matching ref/var allele appears.
///
/// Ports `AmpliconPostProcessModule.countVariantOnAmplicons()`.
fn count_variant_on_amplicons<V: Borrow<Variant>>(
    vref: &Variant,
    good_variants_on_amp: &BTreeMap<usize, Vec<V>>,
) -> i32 {
    let mut gvscnt: i32 = 0;
    for variants in good_variants_on_amp.values() {
        for variant in variants {
            let variant = variant.borrow();
            if variant.ref_allele == vref.ref_allele && variant.var_allele == vref.var_allele {
                gvscnt = gvscnt.saturating_add(1);
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
fn is_amp_bias_flag<V: Borrow<Variant>>(good_variants_on_amp: &BTreeMap<usize, Vec<V>>) -> bool {
    if good_variants_on_amp.len() <= 1 {
        return false;
    }

    let mut current_sorted: Vec<&Variant> = Vec::new();
    let mut next_sorted: Vec<&Variant> = Vec::new();
    let mut values_iter = good_variants_on_amp.values();
    let Some(mut prev_variants) = values_iter.next() else {
        return false;
    };

    // BTreeMap values are already ordered by ascending amplicon index.
    for next_variants in values_iter {
        if prev_variants.len() != next_variants.len() {
            return true;
        }

        if let ([v1], [v2]) = (prev_variants.as_slice(), next_variants.as_slice()) {
            if v1.borrow().description != v2.borrow().description {
                return true;
            }
            prev_variants = next_variants;
            continue;
        }

        if !prev_variants.is_empty() {
            // Sort both by totalPosCoverage descending for comparison
            current_sorted.clear();
            current_sorted.extend(prev_variants.iter().map(Borrow::borrow));
            current_sorted.sort_by_key(|v| Reverse(v.total_pos_cov));

            next_sorted.clear();
            next_sorted.extend(next_variants.iter().map(Borrow::borrow));
            next_sorted.sort_by_key(|v| Reverse(v.total_pos_cov));

            for (var1, var2) in current_sorted.iter().zip(next_sorted.iter()) {
                if var1.description != var2.description {
                    return true;
                }
            }
        }

        prev_variants = next_variants;
    }
    false
}

/// First altered 1-based genomic coordinate for paired somatic amplicon insert-bound checks:
/// advances `start_pos` by 1 bp for 5'-VCF-padded pure deletions (`REF[0] == ALT[0]`, `|ALT| == 1 < |REF|`)
/// and 1-bp-upstream-anchored symbolic SVs (`<DUP>`, `<DEL>`, `<INV>` under `canonical_sv_anchoring`).
#[inline]
fn paired_amp_effective_start(
    start_pos: i32,
    end_pos: i32,
    ref_allele: &str,
    var_allele: &str,
    canonical_sv_anchoring: bool,
) -> i32 {
    let is_vcf_padded_del = end_pos > start_pos
        && ref_allele.len() > 1
        && var_allele.len() == 1
        && ref_allele.as_bytes().first() == var_allele.as_bytes().first();
    let is_anchored_symbolic_sv = canonical_sv_anchoring
        && var_allele.starts_with('<')
        && end_pos > start_pos
        && start_pos >= 1;
    if is_vcf_padded_del || is_anchored_symbolic_sv {
        start_pos.saturating_add(1)
    } else {
        start_pos
    }
}

type CoveringRegionSlice<'a> = std::sync::Arc<[(usize, &'a Region)]>;

/// Single allele entry `(ref_allele, var_allele, concordance)` inside [`SampleAmpliconConcordance`].
pub type AlleleConcordanceEntry = (Box<str>, Box<str>, crate::config::AmpliconConcordance);

fn covering_arc_exact<'a>(segment: &'a [Region], pos: i32) -> CoveringRegionSlice<'a> {
    let mut covering = Vec::new();
    for (idx, region) in segment.iter().enumerate() {
        if let Some(amp) = region.amplicon.as_ref() {
            let start = amp.insert_start.0 as i32;
            let end = amp.insert_end.0 as i32;
            if start <= pos && pos <= end {
                covering.push((idx, region));
            }
        }
    }
    covering.into()
}

/// Per-sample amplicon concordance lookup table produced by [`aggregate_sample_amplicon_vars`]
/// for paired somatic amplicon calling (`AstraZeneca-NGS/VarDictJava#378`,
/// `AstraZeneca-NGS/VarDictJava#155`, `AstraZeneca-NGS/VarDict#109`).
#[derive(Debug, Default, Clone)]
pub struct SampleAmpliconConcordance {
    /// Per-allele concordance keyed by post-`adj_complex` `start_pos` -> `(ref_allele, var_allele, concordance)`.
    pub allele_concordance: rustc_hash::FxHashMap<i32, Vec<AlleleConcordanceEntry>>,
    /// Per-position `(ncamp, ampflag)` baseline across covering amplicons for 0-alt-read / LOH fallback sites.
    pub position_baseline: rustc_hash::FxHashMap<i32, (i32, bool)>,
    /// Whether 1-bp-upstream symbolic SV anchoring (`CompatProfile::canonical_sv_vcf_anchoring`) is active.
    pub canonical_sv_anchoring: bool,
}

impl SampleAmpliconConcordance {
    /// Resolve [`crate::config::AmpliconConcordance`] for `(start_pos, end_pos, ref_allele, var_allele)`
    /// against `segment` with zero heap allocations, falling back to `(gdamp: 0, tlamp, ncamp, ampflag)`
    /// when the allele had zero in-insert candidate reads in this sample.
    pub fn resolve(
        &self,
        start_pos: i32,
        end_pos: i32,
        ref_allele: &str,
        var_allele: &str,
        segment: &[Region],
    ) -> crate::config::AmpliconConcordance {
        if let Some(entries) = self.allele_concordance.get(&start_pos) {
            for (r, a, conc) in entries {
                if r.as_ref() == ref_allele && a.as_ref() == var_allele {
                    return *conc;
                }
            }
        }
        let eff_start = paired_amp_effective_start(
            start_pos,
            end_pos,
            ref_allele,
            var_allele,
            self.canonical_sv_anchoring,
        );
        let tlamp = segment
            .iter()
            .filter(|reg| {
                let (ins_lo, ins_hi) = region_target_range_i32(reg);
                eff_start >= ins_lo && end_pos <= ins_hi
            })
            .count() as i32;
        let (ncamp, ampflag) = self
            .position_baseline
            .get(&start_pos)
            .or_else(|| {
                (eff_start != start_pos)
                    .then(|| self.position_baseline.get(&eff_start))
                    .flatten()
            })
            .copied()
            .unwrap_or((0, false));
        crate::config::AmpliconConcordance {
            gdamp: 0,
            tlamp,
            ncamp,
            ampflag,
        }
    }
}

/// Returns `true` when `v` extends into a 5' or 3' primer outside `amp_reg`'s insert bounds.
#[inline]
pub(crate) fn variant_overlaps_primer(
    v: &Variant,
    amp_reg: &Region,
    canonical_sv_anchoring: bool,
) -> bool {
    amp_reg.amplicon.as_ref().is_some_and(|a| {
        let eff_start = paired_amp_effective_start(
            v.start_pos,
            v.end_pos,
            &v.ref_allele,
            &v.var_allele,
            canonical_sv_anchoring,
        );
        eff_start < a.insert_start.0 as i32 || v.end_pos > a.insert_end.0 as i32
    })
}

/// Aggregate per-amplicon `Vars` for a single sample (tumor or normal) across `segment`
/// (`AstraZeneca-NGS/VarDictJava#378`, `AstraZeneca-NGS/VarDictJava#155`,
/// `AstraZeneca-NGS/VarDict#109`, `palakpsheth/VarDictRust#95`).
pub fn aggregate_sample_amplicon_vars(
    config: &Config,
    segment: &[Region],
    vars_per_amplicon: Vec<crate::javamap::JavaHashMap<i32, Vars>>,
    splice: &rustc_hash::FxHashSet<String>,
) -> (
    crate::javamap::JavaHashMap<i32, Vars>,
    SampleAmpliconConcordance,
) {
    aggregate_sample_amplicon_vars_with_ref(config, segment, vars_per_amplicon, splice, None)
}

pub(crate) fn aggregate_sample_amplicon_vars_with_ref(
    config: &Config,
    segment: &[Region],
    vars_per_amplicon: Vec<crate::javamap::JavaHashMap<i32, Vars>>,
    splice: &rustc_hash::FxHashSet<String>,
    ref_info: Option<(&crate::reference::Reference, i32)>,
) -> (
    crate::javamap::JavaHashMap<i32, Vars>,
    SampleAmpliconConcordance,
) {
    let seg_chr = segment.first().map_or("", |r| r.chr.as_str());
    let ref_ctx = ref_info.map(|(r, len)| (r, seg_chr, len));
    let allow_padded_sv = config.compat.canonical_sv_vcf_anchoring;
    let has_symbolic_sv_at = |p: i32| -> bool {
        vars_per_amplicon.iter().any(|vm| {
            vm.get(&p).is_some_and(|vt| {
                vt.variants.iter().any(|v| {
                    v.var_allele.starts_with('<') && v.end_pos > v.start_pos && v.start_pos >= 1
                })
            })
        })
    };
    let out = if config.pileup {
        let mut extra_sv_positions: BTreeSet<i32> = BTreeSet::new();
        if allow_padded_sv {
            for vm in &vars_per_amplicon {
                for &p in vm.keys() {
                    if covering_arc_exact(segment, p).is_empty() && has_symbolic_sv_at(p) {
                        extra_sv_positions.insert(p);
                    }
                }
            }
        }
        let extra_covs: Vec<(i32, CoveringRegionSlice<'_>)> = extra_sv_positions
            .into_iter()
            .filter_map(|p| {
                let next_cov = covering_arc_exact(segment, p.saturating_add(1));
                (!next_cov.is_empty()).then_some((p, next_cov))
            })
            .collect();
        aggregate_sample_amplicon_positions(
            config,
            &vars_per_amplicon,
            amplicon_position_ref_stream(segment).chain(extra_covs),
            splice,
            ref_ctx,
        )
    } else {
        let mut positions: BTreeSet<i32> = BTreeSet::new();
        for vm in &vars_per_amplicon {
            positions.extend(vm.keys().copied());
        }
        let mut sweep = AmpliconCoverSweep::new(segment);
        aggregate_sample_amplicon_positions(
            config,
            &vars_per_amplicon,
            positions.into_iter().filter_map(|p| {
                let cov = sweep.covering_arc_at(p);
                if !cov.is_empty() {
                    return Some((p, cov));
                }
                if allow_padded_sv && has_symbolic_sv_at(p) {
                    let next_cov = covering_arc_exact(segment, p.saturating_add(1));
                    if !next_cov.is_empty() {
                        return Some((p, next_cov));
                    }
                }
                None
            }),
            splice,
            ref_ctx,
        )
    };
    drop(vars_per_amplicon);
    out
}

fn upsert_concordance_entry(
    concordance: &mut SampleAmpliconConcordance,
    key_start: i32,
    key_ref: &str,
    key_var: &str,
    cand_conc: crate::config::AmpliconConcordance,
    nocov: i32,
    current_flag: bool,
) {
    let entries = concordance.allele_concordance.entry(key_start).or_default();
    if let Some((_, _, existing)) = entries
        .iter_mut()
        .find(|(r, a, _)| r.as_ref() == key_ref && a.as_ref() == key_var)
    {
        if (cand_conc.gdamp, cand_conc.tlamp) > (existing.gdamp, existing.tlamp) {
            *existing = cand_conc;
        }
    } else {
        entries.push((
            key_ref.to_owned().into_boxed_str(),
            key_var.to_owned().into_boxed_str(),
            cand_conc,
        ));
    }
    concordance
        .position_baseline
        .entry(key_start)
        .or_insert((nocov, current_flag));
}

fn aggregate_sample_amplicon_positions<R: Borrow<Region>, C: AsRef<[(usize, R)]>>(
    config: &Config,
    vars_per_amplicon: &[crate::javamap::JavaHashMap<i32, Vars>],
    positions: impl IntoIterator<Item = (i32, C)>,
    splice: &rustc_hash::FxHashSet<String>,
    ref_ctx: Option<(&crate::reference::Reference, &str, i32)>,
) -> (
    crate::javamap::JavaHashMap<i32, Vars>,
    SampleAmpliconConcordance,
) {
    let exact_tie = config.compat.exact_binary_half_even;
    let canonical_sv_anchoring = config.compat.canonical_sv_vcf_anchoring;
    let skip_eq_non_pileup = !config.pileup && config.compat.sanitize_iupac_alt_alleles;
    let refresh_vt = config.compat.recompute_vartype_after_adj_complex;
    let do_left_align = crate::variations::canonical_left_align_enabled(config);
    let indel_only_shift3 = config.compat.indel_only_shift3;
    let min_reads_bias = i32::try_from(config.min_reads_bias).unwrap_or(i32::MAX);
    let mut aligned_vars = crate::javamap::JavaHashMap::new();
    let mut concordance = SampleAmpliconConcordance {
        canonical_sv_anchoring,
        ..SampleAmpliconConcordance::default()
    };

    let mut gvs: Vec<(&Variant, usize)> = Vec::new();
    let mut vcovs: Vec<i32> = Vec::new();
    let mut good_variants_on_amp: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
    let mut best_by_desc: indexmap::IndexMap<
        &str,
        (&Variant, f64, i32, bool),
        crate::variations::FxBuildHasher,
    > = indexmap::IndexMap::default();

    for (position, amplicon_regions) in positions {
        let amplicon_regions = amplicon_regions.as_ref();
        if amplicon_regions.is_empty() {
            continue;
        }
        gvs.clear();
        vcovs.clear();
        good_variants_on_amp.clear();
        best_by_desc.clear();

        let mut maxcov: i32 = 0;
        let mut best_ref: Option<&Variant> = None;
        let mut best_skipped_eq: Option<&Variant> = None;
        let mut fallback_sv: Option<&str> = None;
        let mut top_sv: Option<&str> = None;

        for (amp_num_ref, amp_reg_item) in amplicon_regions {
            let amp_num = *amp_num_ref;
            let amp_reg: &Region = amp_reg_item.borrow();
            let vtmp = vars_per_amplicon
                .get(amp_num)
                .and_then(|m| m.get(&position));
            if let Some(vt) = vtmp {
                if !vt.sv.is_empty() && fallback_sv.is_none() {
                    fallback_sv = Some(vt.sv.as_str());
                }
            }
            let ref_amp = vtmp
                .and_then(|v| v.ref_variant.as_ref())
                .filter(|rv| !variant_overlaps_primer(rv, amp_reg, canonical_sv_anchoring));
            if let Some(rv) = ref_amp {
                if best_ref.is_none_or(|cur| rv.total_pos_cov > cur.total_pos_cov) {
                    best_ref = Some(rv);
                }
            }

            let mut any_in_insert = false;
            let mut skipped_eq_cov: i32 = 0;
            if let Some(vt) = vtmp {
                let mut good_vars: Vec<&Variant> = Vec::new();
                for tv in &vt.variants {
                    if variant_overlaps_primer(tv, amp_reg, canonical_sv_anchoring) {
                        continue;
                    }
                    if skip_eq_non_pileup && tv.ref_allele == tv.var_allele {
                        skipped_eq_cov = skipped_eq_cov.max(tv.total_pos_cov);
                        if tv.total_pos_cov > 0
                            && best_skipped_eq
                                .is_none_or(|cur| tv.total_pos_cov > cur.total_pos_cov)
                        {
                            best_skipped_eq = Some(tv);
                        }
                        continue;
                    }
                    any_in_insert = true;
                    vcovs.push(tv.total_pos_cov);
                    if tv.total_pos_cov > maxcov {
                        maxcov = tv.total_pos_cov;
                    }
                    let vtype = tv.var_type_cow();
                    let is_good = tv.is_good_var(ref_amp, &vtype, splice, config);
                    if is_good {
                        gvs.push((tv, amp_num));
                        good_vars.push(tv);
                    }
                    let af_4dp = crate::types::round_half_even_compat(tv.frequency, 4, exact_tie);
                    match best_by_desc.entry(tv.description.as_str()) {
                        indexmap::map::Entry::Vacant(e) => {
                            e.insert((tv, af_4dp, tv.total_pos_cov, is_good));
                            if top_sv.is_none() && !vt.sv.is_empty() {
                                top_sv = Some(vt.sv.as_str());
                            }
                        }
                        indexmap::map::Entry::Occupied(mut e) => {
                            let cur = e.get_mut();
                            let merged_good = cur.3 || is_good;
                            let better = is_good
                                .cmp(&cur.3)
                                .then_with(|| {
                                    af_4dp
                                        .partial_cmp(&cur.1)
                                        .unwrap_or(std::cmp::Ordering::Equal)
                                })
                                .then_with(|| tv.total_pos_cov.cmp(&cur.2))
                                == std::cmp::Ordering::Greater;
                            if better {
                                *cur = (tv, af_4dp, tv.total_pos_cov, merged_good);
                                if !vt.sv.is_empty() {
                                    top_sv = Some(vt.sv.as_str());
                                }
                            } else {
                                cur.3 = merged_good;
                            }
                        }
                    }
                }
                if !good_vars.is_empty() {
                    good_variants_on_amp.insert(amp_num, good_vars);
                }
            }
            if !any_in_insert {
                let cov = ref_amp.map_or(0, |rv| rv.total_pos_cov).max(skipped_eq_cov);
                vcovs.push(cov);
            }
        }

        let eff_maxcov = if maxcov > 0 {
            maxcov
        } else {
            vcovs.iter().copied().max().unwrap_or(0)
        };
        let nocov = vcovs
            .iter()
            .filter(|&&t| (t as f64) < (eff_maxcov as f64) / 50.0)
            .count() as i32;

        if gvs.len() > 1 {
            gvs.sort_by(|a, b| {
                b.0.frequency
                    .partial_cmp(&a.0.frequency)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }

        let mut current_flag = is_amp_bias_flag(&good_variants_on_amp);
        if current_flag {
            if let Some(&(top_gv, _)) = gvs.first() {
                let gdnt = &top_gv.description;
                let mut gcnt_len = 0usize;
                for (amp_num_ref, amp_reg_item) in amplicon_regions {
                    let amp_num = *amp_num_ref;
                    let amp_reg: &Region = amp_reg_item.borrow();
                    if let Some(vtmp) = vars_per_amplicon
                        .get(amp_num)
                        .and_then(|m| m.get(&position))
                    {
                        if let Some(variant) = vtmp.variant_by_desc(gdnt) {
                            if (!skip_eq_non_pileup || variant.ref_allele != variant.var_allele)
                                && !variant_overlaps_primer(
                                    variant,
                                    amp_reg,
                                    canonical_sv_anchoring,
                                )
                            {
                                let ref_amp = vtmp.ref_variant.as_ref().filter(|rv| {
                                    !variant_overlaps_primer(rv, amp_reg, canonical_sv_anchoring)
                                });
                                if variant.is_good_var(ref_amp, "", splice, config) {
                                    gcnt_len = gcnt_len.saturating_add(1);
                                }
                            }
                        }
                    }
                }
                if gcnt_len == gvs.len() {
                    current_flag = false;
                }
            }
        }

        concordance
            .position_baseline
            .insert(position, (nocov, current_flag));

        for &(rep_var, _, _, _) in best_by_desc.values() {
            let gdamp = count_variant_on_amplicons(rep_var, &good_variants_on_amp);
            let rep_eff_start = paired_amp_effective_start(
                rep_var.start_pos,
                rep_var.end_pos,
                &rep_var.ref_allele,
                &rep_var.var_allele,
                canonical_sv_anchoring,
            );
            let mut bad_count: i32 = 0;
            for (amp_num_ref, amp_reg_item) in amplicon_regions {
                let amp_num = *amp_num_ref;
                let amp_reg: &Region = amp_reg_item.borrow();
                if config.pileup && rep_var.ref_allele == rep_var.var_allele {
                    continue;
                }
                let is_in_goodmap = good_variants_on_amp.get(&amp_num).is_some_and(|vars| {
                    vars.iter().any(|tv| {
                        tv.ref_allele == rep_var.ref_allele && tv.var_allele == rep_var.var_allele
                    })
                });
                if is_in_goodmap {
                    continue;
                }
                let (ins_lo, ins_hi) = region_target_range_i32(amp_reg);
                if rep_eff_start >= ins_lo && rep_var.end_pos <= ins_hi {
                    bad_count = bad_count.saturating_add(1);
                }
            }
            let tlamp = gdamp.saturating_add(bad_count);
            let (key_start, key_ref, key_var) = if rep_var.var_type_cow() == "Complex" {
                let mut tmp = rep_var.clone();
                tmp.adj_complex_compat(refresh_vt);
                if skip_eq_non_pileup && tmp.ref_allele == tmp.var_allele {
                    continue;
                }
                (
                    tmp.start_pos,
                    std::borrow::Cow::Owned(tmp.ref_allele),
                    std::borrow::Cow::Owned(tmp.var_allele),
                )
            } else {
                (
                    rep_var.start_pos,
                    std::borrow::Cow::Borrowed(rep_var.ref_allele.as_str()),
                    std::borrow::Cow::Borrowed(rep_var.var_allele.as_str()),
                )
            };
            let cand_conc = crate::config::AmpliconConcordance {
                gdamp,
                tlamp,
                ncamp: nocov,
                ampflag: current_flag,
            };
            upsert_concordance_entry(
                &mut concordance,
                key_start,
                key_ref.as_ref(),
                key_var.as_ref(),
                cand_conc,
                nocov,
                current_flag,
            );
            if do_left_align && rep_var.ref_allele.len() != rep_var.var_allele.len() {
                let mut la_tmp = rep_var.clone();
                la_tmp.var_type = la_tmp.var_type();
                la_tmp.normalize_postprocess(ref_ctx, true, refresh_vt, indel_only_shift3);
                if !(skip_eq_non_pileup && la_tmp.ref_allele == la_tmp.var_allele)
                    && (la_tmp.start_pos != key_start
                        || la_tmp.ref_allele != key_ref.as_ref()
                        || la_tmp.var_allele != key_var.as_ref())
                {
                    upsert_concordance_entry(
                        &mut concordance,
                        la_tmp.start_pos,
                        &la_tmp.ref_allele,
                        &la_tmp.var_allele,
                        cand_conc,
                        nocov,
                        current_flag,
                    );
                }
            }
        }

        if best_ref.is_some() || best_skipped_eq.is_some() || !best_by_desc.is_empty() {
            let mut out_vars = Vars::new();
            out_vars.ref_variant = match best_ref {
                Some(rv) if rv.total_pos_cov > 0 || best_skipped_eq.is_none() => Some(rv.clone()),
                _ => best_skipped_eq.map(|eq| {
                    let mut rv = Variant::new();
                    rv.start_pos = position;
                    rv.end_pos = position;
                    rv.description.clone_from(&eq.var_allele);
                    rv.ref_allele.clone_from(&eq.ref_allele);
                    rv.var_allele.clone_from(&eq.var_allele);
                    let rf = eq.ref_fwd_cov.saturating_add(eq.fwd);
                    let rr = eq.ref_rev_cov.saturating_add(eq.rev);
                    rv.total_pos_cov = eq.total_pos_cov;
                    rv.ref_fwd_cov = rf;
                    rv.ref_rev_cov = rr;
                    rv.genotype = eq
                        .genotype
                        .clone()
                        .or_else(|| Some(format!("{}/{}", eq.ref_allele, eq.ref_allele)));
                    let ref_bias =
                        crate::variations::util::strand_bias(rf, rr, 0.05, min_reads_bias);
                    rv.bias.clear();
                    if best_by_desc.is_empty() {
                        rv.fwd = 0;
                        rv.rev = 0;
                        let _ =
                            std::fmt::Write::write_fmt(&mut rv.bias, format_args!("{ref_bias};0"));
                    } else {
                        rv.fwd = rf;
                        rv.rev = rr;
                        rv.hicnt = eq.hicnt;
                        rv.mean_mapping_quality = eq.mean_mapping_quality;
                        let _ =
                            std::fmt::Write::write_fmt(&mut rv.bias, format_args!("{ref_bias}"));
                    }
                    rv
                }),
            };
            if !best_by_desc.is_empty() {
                let mut ordered: Vec<(&Variant, f64, i32, bool)> =
                    best_by_desc.values().copied().collect();
                ordered.sort_by(|a, b| {
                    b.3.cmp(&a.3)
                        .then_with(|| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal))
                        .then_with(|| b.2.cmp(&a.2))
                        .then_with(|| a.0.description.cmp(&b.0.description))
                });
                out_vars.variants.reserve(ordered.len());
                for (v, _, _, _) in ordered {
                    out_vars.variants.push(v.clone());
                }
            }
            if let Some(sv_str) = top_sv.or(fallback_sv) {
                out_vars.sv.push_str(sv_str);
            }
            aligned_vars.insert(position, out_vars);
        }
    }

    (aligned_vars, concordance)
}

#[cfg(test)]
#[allow(deprecated)]
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
        let splice = rustc_hash::FxHashSet::default();

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
        let splice = rustc_hash::FxHashSet::default();

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
        let splice = rustc_hash::FxHashSet::default();

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
        let splice = rustc_hash::FxHashSet::default();

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
        let splice = rustc_hash::FxHashSet::default();

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
        let splice = rustc_hash::FxHashSet::default();

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
        let borrowed_empty: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
        assert!(!is_amp_bias_flag(&borrowed_empty));

        // Single amplicon (len == 1) with 0, 1, and >= 2 good variants returns false immediately.
        let v_g = make_good_variant(150, "A", "G", 0.10);
        let v_t = make_good_variant(150, "A", "T", 0.08);
        let mut single_amp0: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
        single_amp0.insert(3, vec![]);
        assert!(!is_amp_bias_flag(&single_amp0));
        single_amp0.insert(3, vec![&v_g]);
        assert!(!is_amp_bias_flag(&single_amp0));
        single_amp0.insert(3, vec![&v_g, &v_t]);
        assert!(!is_amp_bias_flag(&single_amp0));
    }

    #[test]
    fn test_is_amp_bias_flag_same_variants() {
        let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        good_variants_on_amp.insert(0, vec![make_good_variant(150, "A", "G", 0.10)]);
        good_variants_on_amp.insert(1, vec![make_good_variant(150, "A", "G", 0.12)]);
        assert!(!is_amp_bias_flag(&good_variants_on_amp));

        // Borrowed &Variant across 3 non-contiguous amplicon indices with 0, 1, and >= 2 good variants
        let mut empty_multi: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
        empty_multi.insert(0, vec![]);
        empty_multi.insert(4, vec![]);
        assert!(!is_amp_bias_flag(&empty_multi));

        let mut v0_a = make_good_variant(150, "A", "G", 0.10);
        v0_a.total_pos_cov = 80;
        let mut v0_b = make_good_variant(150, "A", "T", 0.06);
        v0_b.total_pos_cov = 120;
        let mut v2_b = make_good_variant(150, "A", "T", 0.09);
        v2_b.total_pos_cov = 200;
        let mut v2_a = make_good_variant(150, "A", "G", 0.14);
        v2_a.total_pos_cov = 90;
        let mut multi_borrowed: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
        // Unsorted insertion order within each Vec; after sorting by total_pos_cov desc both are [T, G]
        multi_borrowed.insert(0, vec![&v0_a, &v0_b]);
        multi_borrowed.insert(2, vec![&v2_b, &v2_a]);
        multi_borrowed.insert(5, vec![&v0_b, &v0_a]);
        assert!(!is_amp_bias_flag(&multi_borrowed));
    }

    #[test]
    fn test_is_amp_bias_flag_different_variants() {
        let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
        good_variants_on_amp.insert(0, vec![make_good_variant(150, "A", "G", 0.10)]);
        good_variants_on_amp.insert(1, vec![make_good_variant(150, "A", "T", 0.12)]);
        assert!(is_amp_bias_flag(&good_variants_on_amp));

        // >= 2 good variants per amplicon where total_pos_cov ordering exposes description mismatch
        let mut v_g_hi = make_good_variant(150, "A", "G", 0.10);
        v_g_hi.total_pos_cov = 150;
        let mut v_t_lo = make_good_variant(150, "A", "T", 0.08);
        v_t_lo.total_pos_cov = 90;
        let mut v_c_lo = make_good_variant(150, "A", "C", 0.07);
        v_c_lo.total_pos_cov = 90;
        let mut borrowed_diff: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
        borrowed_diff.insert(1, vec![&v_t_lo, &v_g_hi]);
        borrowed_diff.insert(4, vec![&v_g_hi, &v_c_lo]);
        assert!(is_amp_bias_flag(&borrowed_diff));
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

        // Multi-allele borrowed (&Variant, usize) keeps top-AF representative per (ref, var) pair
        let v_ag_hi = make_good_variant(150, "A", "G", 0.35);
        let v_ag_lo = make_good_variant(150, "A", "G", 0.15);
        let v_at = make_good_variant(150, "A", "T", 0.22);
        let v_ins = make_good_variant(150, "A", "ACC", 0.11);
        let borrowed_gvs = vec![(&v_ag_hi, 0usize), (&v_at, 1), (&v_ag_lo, 2), (&v_ins, 0)];
        let borrowed_vrefs = fill_vref_list(&borrowed_gvs);
        assert_eq!(borrowed_vrefs.len(), 3);
        assert_eq!(borrowed_vrefs[0].var_allele, "G");
        assert_eq!(borrowed_vrefs[0].frequency, 0.35);
        assert_eq!(borrowed_vrefs[1].var_allele, "T");
        assert_eq!(borrowed_vrefs[2].var_allele, "ACC");
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

        let mut borrowed_on_amp: BTreeMap<usize, Vec<&Variant>> = BTreeMap::new();
        borrowed_on_amp.insert(0, vec![&vref]);
        borrowed_on_amp.insert(2, vec![&vref, &vref_t]);
        assert_eq!(count_variant_on_amplicons(&vref, &borrowed_on_amp), 2);
        assert_eq!(count_variant_on_amplicons(&vref_t, &borrowed_on_amp), 1);
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
        // the sweep is output-proportional + O((P + A) log A) (~0.1 s), rebuilding
        // its cached `Arc` only on `<= 2 * amps` start/end transitions (verified
        // deterministically via `Arc::ptr_eq`, with a 10 s hang backstop).
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
        let mut arc_transitions = 0usize;
        let mut last_arc: Option<std::sync::Arc<[(usize, &Region)]>> = None;
        for (_, cov_arc) in amplicon_position_ref_stream(&segment) {
            if !last_arc
                .as_ref()
                .is_some_and(|prev| std::sync::Arc::ptr_eq(prev, &cov_arc))
            {
                arc_transitions += 1;
                last_arc = Some(cov_arc);
            }
        }
        let elapsed = t0.elapsed();
        assert_eq!(npos, (amps as u64 - 1) * 100 + 200);
        assert_eq!(nmem, amps as u64 * 200);
        assert!(
            arc_transitions <= 2 * amps,
            "amplicon_position_ref_stream must rebuild cached Arc only on start/end transitions, got {arc_transitions}"
        );
        assert!(
            elapsed < std::time::Duration::from_secs(10),
            "dense-chain sweep hung: {elapsed:?}",
        );
    }

    #[test]
    fn pileup_stream_matches_batch() {
        // Two overlapping amplicons; positions 150-152 cover the good-variant,
        // reference-only, and empty (null row under -p) paths.
        let mut config = make_test_config();
        config.pileup = true;
        let region = make_amp_region("chr1", 100, 300);
        let splice = rustc_hash::FxHashSet::default();

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
                &rustc_hash::FxHashSet::default(),
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

    #[test]
    fn primer_overlap_clears_amp_bias_flag_and_pileup_null_row_zeros_tlamp_ncamp() {
        // 1. AmpliconPostProcessModule.java:108-115: pileup null row (no ref_vars) emits tlamp=0, ncamp=0.
        let mut config = make_test_config();
        config.pileup = true;
        let region = make_amp_region("chr1", 100, 200);
        let splice = rustc_hash::FxHashSet::default();
        let vars_per_amplicon = vec![crate::javamap::JavaHashMap::new()];
        let mut amplicons_on_positions = BTreeMap::new();
        amplicons_on_positions.insert(150, vec![(0, make_amp_region("chr1", 100, 200))]);
        let null_rows = amplicon_post_process(
            &config,
            &region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        )
        .unwrap();
        assert_eq!(null_rows.len(), 1);
        assert_eq!(
            (null_rows[0].gdamp, null_rows[0].tlamp, null_rows[0].ncamp),
            (0, 0, 0)
        );

        // 2. AmpliconPostProcessModule.java:186-188: when current_flag is true and a multi-base
        // variant overlaps an amplicon primer so current_gvscnt < initial_gvscnt, current_flag resets to false.
        config.pileup = false;
        let mut v_del = make_good_variant(150, "ACGTACGT", "A", 0.25);
        v_del.end_pos = 157;
        v_del.description = "-7".to_string();
        let mut v_snv = make_good_variant(150, "A", "T", 0.15);
        v_snv.description = "T".to_string();

        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut v0 = Vars::new();
        v0.variants.push(v_del.clone());
        v0.var_desc_to_variant
            .insert("-7".to_string(), v_del.clone());
        vars_amp0.insert(150, v0);

        let mut vars_amp1 = crate::javamap::JavaHashMap::new();
        let mut v1 = Vars::new();
        // Amp 1 has both v_del and v_snv so is_amp_bias_flag is true (different good-variant counts: 1 vs 2)
        v1.variants.push(v_del.clone());
        v1.variants.push(v_snv.clone());
        v1.var_desc_to_variant.insert("-7".to_string(), v_del);
        v1.var_desc_to_variant.insert("T".to_string(), v_snv);
        vars_amp1.insert(150, v1);

        // Amp 2 has insert_end = 153 so v_del (150..157) overlaps its primer and decrements current_gvscnt (2 -> 1)
        let vars_amp2 = crate::javamap::JavaHashMap::new();
        let mut amp2_reg = make_amp_region("chr1", 100, 160);
        if let Some(ref mut amp) = amp2_reg.amplicon {
            amp.insert_start = crate::types::Pos1(110);
            amp.insert_end = crate::types::Pos1(153);
        }

        let mut amps = BTreeMap::new();
        amps.insert(
            150,
            vec![
                (0, make_amp_region("chr1", 100, 200)),
                (1, make_amp_region("chr1", 100, 200)),
                (2, amp2_reg),
            ],
        );
        let res = amplicon_post_process(
            &config,
            &region,
            &[vars_amp0, vars_amp1, vars_amp2],
            &amps,
            &splice,
        )
        .unwrap();
        assert!(
            !res[0].ampflag,
            "primer overlap decrementing gvscnt must clear ampflag"
        );
    }

    #[test]
    fn amplicon_position_ref_stream_matches_owned_without_cloning() {
        use crate::regions::AmpliconInfo;
        use crate::types::Pos1;

        fn amp_reg(chr: &str, gene: &str, start: u32, end: u32, is: u32, ie: u32) -> Region {
            Region {
                chr: chr.to_string(),
                start: Pos1(start),
                end: Pos1(end),
                gene: gene.to_string(),
                amplicon: Some(AmpliconInfo {
                    insert_start: Pos1(is),
                    insert_end: Pos1(ie),
                }),
            }
        }

        let segment = vec![
            amp_reg("chr1", "GENE_A", 140, 260, 150, 250),
            amp_reg("chr1", "GENE_B", 90, 210, 100, 200),
            amp_reg("chr1", "GENE_C", 500, 510, 500, 510),
            amp_reg("chr1", "GENE_D", 150, 160, 150, 155),
            Region {
                chr: "chr1".to_string(),
                start: Pos1(300),
                end: Pos1(400),
                gene: "NON_AMP".to_string(),
                amplicon: None,
            },
            amp_reg("chr1", "GENE_E", 180, 190, 180, 185),
        ];

        let owned: Vec<(i32, Vec<(usize, Region)>)> = amplicon_position_stream(&segment).collect();
        let borrowed: Vec<_> = amplicon_position_ref_stream(&segment).collect();

        assert_eq!(borrowed.len(), owned.len());
        let mut shared_arc_hits = 0usize;
        for w in borrowed.windows(2) {
            let cov0 = &w[0].1;
            let cov1 = &w[1].1;
            let same_indices = cov0.len() == cov1.len()
                && cov0.iter().zip(cov1.iter()).all(|((a, _), (b, _))| a == b);
            if same_indices {
                assert!(
                    std::sync::Arc::ptr_eq(cov0, cov1),
                    "adjacent positions with identical active set must share cached Arc"
                );
                shared_arc_hits += 1;
            }
        }
        assert!(shared_arc_hits > 0);
        for ((bp, bcov), (op, ocov)) in borrowed.iter().zip(owned.iter()) {
            assert_eq!(bp, op);
            assert_eq!(bcov.len(), ocov.len());
            for ((b_idx, b_reg), (o_idx, o_reg)) in bcov.iter().zip(ocov.iter()) {
                assert_eq!(b_idx, o_idx);
                assert_eq!(*b_reg, o_reg);
                assert!(
                    std::ptr::eq(*b_reg, &segment[*b_idx]),
                    "amplicon_position_ref_stream must borrow Region directly from segment"
                );
                assert_eq!(b_reg.gene.as_ptr(), segment[*b_idx].gene.as_ptr());
            }
        }

        // Also verify amplicon_post_process_streaming produces identical output (including -D debug)
        // when driven by amplicon_position_ref_stream vs amplicon_position_stream.
        let mut config = make_test_config();
        config.debug = true;
        let region = segment[0].clone();
        let splice = rustc_hash::FxHashSet::default();

        let mut vars_per_amplicon = vec![crate::javamap::JavaHashMap::new(); segment.len()];
        let mut v0 = Vars::new();
        v0.variants.push(make_good_variant(152, "A", "G", 0.30));
        vars_per_amplicon[0].insert(152, v0);

        let mut v1 = Vars::new();
        v1.variants.push(make_good_variant(152, "A", "T", 0.20));
        vars_per_amplicon[1].insert(152, v1);

        let mut out_owned = Vec::new();
        amplicon_post_process_streaming(
            &config,
            &region,
            &vars_per_amplicon,
            amplicon_position_stream(&segment),
            &splice,
            &mut |av| {
                out_owned.push(av);
                Ok(())
            },
        )
        .unwrap();

        let mut out_ref = Vec::new();
        amplicon_post_process_streaming(
            &config,
            &region,
            &vars_per_amplicon,
            amplicon_position_ref_stream(&segment),
            &splice,
            &mut |av| {
                out_ref.push(av);
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(out_ref.len(), out_owned.len());
        assert_eq!(out_ref.len(), 2);
        for (got, want) in out_ref.iter().zip(out_owned.iter()) {
            assert_eq!(got.variant.var_allele, want.variant.var_allele);
            assert_eq!(got.gdamp, want.gdamp);
            assert_eq!(got.tlamp, want.tlamp);
            assert_eq!(got.ncamp, want.ncamp);
            assert_eq!(got.ampflag, want.ampflag);
            assert_eq!(got.region, want.region);
            assert_eq!(got.good_variants.len(), want.good_variants.len());
            assert_eq!(got.bad_variants.len(), want.bad_variants.len());
            for ((_, gr), (_, wr)) in got.good_variants.iter().zip(want.good_variants.iter()) {
                assert_eq!(gr, wr);
            }
            for ((_, br), (_, wr)) in got.bad_variants.iter().zip(want.bad_variants.iter()) {
                assert_eq!(br, wr);
            }
        }

        // Verify batch `amplicon_post_process` over borrowed `&Region` entries matches owned.
        let mut sweep = AmpliconCoverSweep::new(&segment);
        let mut borrowed_map: BTreeMap<i32, Vec<(usize, &Region)>> = BTreeMap::new();
        let mut owned_map: BTreeMap<i32, Vec<(usize, Region)>> = BTreeMap::new();
        let cov_ref = sweep.covering_refs_at(152);
        let cov_owned: Vec<(usize, Region)> =
            cov_ref.iter().map(|&(i, r)| (i, r.clone())).collect();
        borrowed_map.insert(152, cov_ref);
        owned_map.insert(152, cov_owned);
        let batch_ref =
            amplicon_post_process(&config, &region, &vars_per_amplicon, &borrowed_map, &splice)
                .unwrap();
        let batch_owned =
            amplicon_post_process(&config, &region, &vars_per_amplicon, &owned_map, &splice)
                .unwrap();
        assert_eq!(format!("{batch_ref:?}"), format!("{batch_owned:?}"));

        // Verify non-`-p` streaming `covering_arc_at` + `amplicon_post_process_streaming`
        // (+ `amplicon_select_vcf_rows` in VCF mode) matches `amplicon_post_process` field-for-field
        // across overlapping multi-amplicon `.var` and VCF modes (including out-of-insert keys and
        // multi-allelic same-start AF contests).
        let mut rich_vars = vars_per_amplicon.clone();
        // Out-of-insert position 95 (uncovered by every amplicon insert -> cov.is_empty())
        let mut v_uncov = Vars::new();
        v_uncov.variants.push(make_good_variant(95, "C", "A", 0.50));
        rich_vars[0].insert(95, v_uncov);
        // Adjacent covered position 153 (shares cached Arc with 152 across amp 0, 1, 3) with
        // multi-allelic AF contest and shifted start_pos = 152 colliding with pos 152 in VCF mode
        let mut v_153_a0 = Vars::new();
        let mut shifted_hi = make_good_variant(152, "AT", "A", 0.65);
        shifted_hi.end_pos = 153;
        v_153_a0.variants.push(shifted_hi.clone());
        v_153_a0
            .variants
            .push(make_good_variant(153, "T", "G", 0.25));
        rich_vars[0].insert(153, v_153_a0);
        let mut v_153_a1 = Vars::new();
        v_153_a1.variants.push(shifted_hi);
        rich_vars[1].insert(153, v_153_a1);
        // Covered position 182 (amp 0, 1, 5)
        let mut v_182 = Vars::new();
        v_182.variants.push(make_good_variant(182, "G", "C", 0.40));
        rich_vars[5].insert(182, v_182);

        for emit_intermediate in [true, false] {
            for debug_mode in [false, true] {
                let mut mode_cfg = make_test_config();
                mode_cfg.emit_intermediate = emit_intermediate;
                mode_cfg.debug = debug_mode;

                let mut positions: BTreeSet<i32> = BTreeSet::new();
                for vm in &rich_vars {
                    positions.extend(vm.keys().copied());
                }

                let mut batch_sweep = AmpliconCoverSweep::new(&segment);
                let mut batch_amps: BTreeMap<i32, Vec<(usize, &Region)>> = BTreeMap::new();
                for &p in &positions {
                    let cov = batch_sweep.covering_refs_at(p);
                    if !cov.is_empty() {
                        batch_amps.insert(p, cov);
                    }
                }
                let expected =
                    amplicon_post_process(&mode_cfg, &region, &rich_vars, &batch_amps, &splice)
                        .unwrap();

                let mut stream_sweep = AmpliconCoverSweep::new(&segment);
                let mut streamed = Vec::new();
                amplicon_post_process_streaming(
                    &mode_cfg,
                    &region,
                    &rich_vars,
                    positions.into_iter().filter_map(|p| {
                        let cov = stream_sweep.covering_arc_at(p);
                        (!cov.is_empty()).then_some((p, cov))
                    }),
                    &splice,
                    &mut |av| {
                        streamed.push(av);
                        Ok(())
                    },
                )
                .unwrap();
                if !mode_cfg.emit_intermediate {
                    streamed = amplicon_select_vcf_rows(streamed);
                }

                assert_eq!(streamed.len(), expected.len());
                assert_eq!(format!("{streamed:?}"), format!("{expected:?}"));
                if debug_mode {
                    for av in &expected {
                        assert!(
                            !av.good_variants.is_empty(),
                            "debug_good must be populated for good variant rows under -D"
                        );
                        for (gv_opt, reg_tag) in &av.good_variants {
                            let gv = gv_opt.as_ref().expect("debug_good entry must be Some");
                            assert_eq!(
                                gv.var_type,
                                gv.var_type(),
                                "debug_good var_type must be populated even on AMPBIAS gcnt rows"
                            );
                            assert!(!gv.var_type.is_empty());
                            assert!(!reg_tag.is_empty());
                        }
                    }
                } else {
                    for av in &expected {
                        assert!(av.good_variants.is_empty());
                        assert!(av.bad_variants.is_empty());
                    }
                }
            }
        }

        // Verify pileup == false/true x debug == false/true across ref-only, null-row,
        // Complex variant adjustment, and multi-allele vref_list emission.
        let mut pileup_vars = rich_vars.clone();
        // Add ref_variant on pos 151 across two amplicons with different total_pos_cov (90 vs 140)
        let mut ref_lo = Variant::new();
        ref_lo.ref_allele = "T".to_string();
        ref_lo.var_allele = "T".to_string();
        ref_lo.start_pos = 151;
        ref_lo.end_pos = 151;
        ref_lo.total_pos_cov = 90;
        let mut v151_0 = Vars::new();
        v151_0.ref_variant = Some(ref_lo);
        pileup_vars[0].insert(151, v151_0);

        let mut ref_hi = Variant::new();
        ref_hi.ref_allele = "T".to_string();
        ref_hi.var_allele = "T".to_string();
        ref_hi.start_pos = 151;
        ref_hi.end_pos = 151;
        ref_hi.total_pos_cov = 140;
        let mut v151_1 = Vars::new();
        v151_1.ref_variant = Some(ref_hi);
        pileup_vars[1].insert(151, v151_1);

        // Add Complex variant at pos 154 ("AT" -> "GC" with shared prefix "C" -> adj_complex)
        let mut v_complex = make_good_variant(154, "CAT", "CGC", 0.42);
        v_complex.end_pos = 156;
        v_complex.description = "C#GC^2".to_string();
        let mut v154_0 = Vars::new();
        v154_0.variants.push(v_complex);
        pileup_vars[0].insert(154, v154_0);

        for pileup_mode in [false, true] {
            for debug_mode in [false, true] {
                let mut pcfg = make_test_config();
                pcfg.pileup = pileup_mode;
                pcfg.debug = debug_mode;
                pcfg.emit_intermediate = true;

                let window_owned = amplicon_position_stream(&segment)
                    .filter(|(p, _)| (150..=154).contains(p))
                    .collect::<BTreeMap<_, _>>();
                let batch_rows =
                    amplicon_post_process(&pcfg, &region, &pileup_vars, &window_owned, &splice)
                        .unwrap();

                let mut stream_rows = Vec::new();
                amplicon_post_process_streaming(
                    &pcfg,
                    &region,
                    &pileup_vars,
                    amplicon_position_ref_stream(&segment).filter(|(p, _)| (150..=154).contains(p)),
                    &splice,
                    &mut |av| {
                        stream_rows.push(av);
                        Ok(())
                    },
                )
                .unwrap();

                assert_eq!(format!("{stream_rows:?}"), format!("{batch_rows:?}"));
                let mut rendered_38 = String::new();
                let mut rendered_40 = String::new();
                for av in &stream_rows {
                    crate::output::table::amplicon_38_append(
                        &mut rendered_38,
                        av,
                        &region,
                        "AMP_S1",
                        debug_mode,
                    );
                    rendered_38.push('\n');
                    crate::output::table::amplicon_40_append(
                        &mut rendered_40,
                        av,
                        &region,
                        "AMP_S1",
                        debug_mode,
                    );
                    rendered_40.push('\n');
                }
                assert!(!rendered_38.is_empty());
                assert!(!rendered_40.is_empty());

                if pileup_mode {
                    // Pos 150 is null row, pos 151 picks highest total_pos_cov (140),
                    // pos 152 emits 2 alleles (G, T), pos 153 emits 2 alleles (A, G), pos 154 Complex (start shifted to 155).
                    assert_eq!(stream_rows.len(), 7);
                    assert!(stream_rows[0].variant.null_record);
                    assert_eq!(stream_rows[1].variant.total_pos_cov, 140);
                } else {
                    // Only positions 152 (2 rows), 153 (2 rows), 154 (1 Complex row) emit.
                    assert_eq!(stream_rows.len(), 5);
                }
                let complex_row = stream_rows.last().unwrap();
                assert_eq!(complex_row.variant.var_type, "Complex");
                assert_eq!(complex_row.variant.start_pos, 155);
                if debug_mode {
                    for av in &stream_rows {
                        for (gv_opt, _) in &av.good_variants {
                            let gv = gv_opt.as_ref().unwrap();
                            assert_eq!(gv.var_type, gv.var_type());
                            assert!(!gv.var_type.is_empty());
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_amplicon_vcf_all_alleles_sentinels_and_multi_amplicon_strict_bounds() {
        use crate::regions::AmpliconInfo;
        use crate::types::Pos1;

        let mut base_cfg = Config::default();
        base_cfg.freq = 0.01;
        base_cfg.min_reads = 2;
        base_cfg.read_pos = 5.0;
        base_cfg.base_qual = 22.5;
        base_cfg.qratio = 1.5;
        base_cfg.amplicon = Some("10:0.95".to_string());

        let amp0 = Region {
            chr: "chr1".to_string(),
            start: Pos1(90),
            end: Pos1(160),
            gene: "G1".to_string(),
            amplicon: Some(AmpliconInfo {
                insert_start: Pos1(100),
                insert_end: Pos1(150),
            }),
        };
        let amp1 = Region {
            chr: "chr1".to_string(),
            start: Pos1(170),
            end: Pos1(230),
            gene: "G1".to_string(),
            amplicon: Some(AmpliconInfo {
                insert_start: Pos1(180),
                insert_end: Pos1(220),
            }),
        };
        let segment = vec![amp0.clone(), amp1.clone()];
        let splice = rustc_hash::FxHashSet::default();

        let mk_good = |start: i32, end: i32, r: &str, a: &str, af: f64| -> Variant {
            let mut v = Variant::new();
            v.start_pos = start;
            v.end_pos = end;
            v.ref_allele = r.to_string();
            v.var_allele = a.to_string();
            v.description = a.to_string();
            v.frequency = af;
            v.hicnt = 10;
            v.pos_coverage = 10;
            v.total_pos_cov = 50;
            v.mean_position = 15.0;
            v.pstd = true;
            v.mean_quality = 35.0;
            v.qstd = true;
            v.hq_lq_ratio = 3.0;
            v.mean_mapping_quality = 60.0;
            v
        };

        // 1. vcf_all_alleles = true preserves both alleles at pos 120 in VCF mode.
        let mut vars_amp0 = crate::javamap::JavaHashMap::new();
        let mut v120 = Vars::new();
        v120.variants.push(mk_good(120, 120, "A", "G", 0.45));
        v120.variants.push(mk_good(120, 120, "A", "T", 0.20));
        vars_amp0.insert(120, v120);
        // Left-anchored deletion at 99 (outside insert [100, 150]) and shifted indel at 155 (in gap [151, 179])
        let mut v100 = Vars::new();
        v100.variants.push(mk_good(99, 100, "CA", "C", 0.35));
        vars_amp0.insert(100, v100);
        let mut v150 = Vars::new();
        v150.variants.push(mk_good(155, 155, "T", "TAA", 0.35));
        vars_amp0.insert(150, v150);

        let mut vars_amp1 = crate::javamap::JavaHashMap::new();
        let mut v190 = Vars::new();
        v190.variants.push(mk_good(190, 190, "C", "A", 0.40));
        vars_amp1.insert(190, v190);
        let mut v220 = Vars::new();
        v220.variants.push(mk_good(225, 225, "G", "GTT", 0.30));
        vars_amp1.insert(220, v220);

        let vars_per_amp = vec![vars_amp0, vars_amp1];
        let mut amps_on_pos: BTreeMap<i32, Vec<(usize, &Region)>> = BTreeMap::new();
        amps_on_pos.insert(100, vec![(0, &amp0)]);
        amps_on_pos.insert(120, vec![(0, &amp0)]);
        amps_on_pos.insert(150, vec![(0, &amp0)]);
        amps_on_pos.insert(190, vec![(1, &amp1)]);
        amps_on_pos.insert(220, vec![(1, &amp1)]);

        // Pass `amp1` as `region` (mirroring `segment.last()` in `run_amplicon`, AR-4).
        let single_winner =
            amplicon_post_process(&base_cfg, &amp1, &vars_per_amp, &amps_on_pos, &splice).unwrap();
        assert_eq!(
            single_winner
                .iter()
                .filter(|av| av.variant.start_pos == 120)
                .count(),
            1
        );

        let mut all_cfg = base_cfg.clone();
        all_cfg.vcf_all_alleles = true;
        let all_out =
            amplicon_post_process(&all_cfg, &amp1, &vars_per_amp, &amps_on_pos, &splice).unwrap();
        assert_eq!(
            all_out
                .iter()
                .filter(|av| av.variant.start_pos == 120)
                .count(),
            2
        );

        // 2. strict_region_bounds (#369, AR-4) under COMPAT_OFF:
        //    Even though `region` is `amp1` ([180, 220]), `amp0`'s insert [100, 150] is ALSO kept
        //    (pos 120 survives), while out-of-insert shifted positions 99, 155 (in gap), and 225 are dropped!
        let mut strict_all_cfg = all_cfg.clone();
        strict_all_cfg.compat = crate::config::COMPAT_OFF;
        let strict_out =
            amplicon_post_process(&strict_all_cfg, &amp1, &vars_per_amp, &amps_on_pos, &splice)
                .unwrap();
        let strict_starts: Vec<i32> = strict_out.iter().map(|av| av.variant.start_pos).collect();
        assert_eq!(strict_starts, vec![120, 120, 190]);

        // 3. suppress_pileup_sentinels (UB-2) under COMPAT_OFF:
        let mut p_vars = vec![
            crate::javamap::JavaHashMap::new(),
            crate::javamap::JavaHashMap::new(),
        ];
        let mut ref_101 = Vars::new();
        ref_101.ref_variant = Some(mk_good(101, 101, "G", "G", 1.0));
        p_vars[0].insert(101, ref_101);

        let mut p_cfg = base_cfg.clone();
        p_cfg.pileup = true;
        p_cfg.emit_intermediate = true;
        p_cfg.compat = crate::config::COMPAT_OFF;
        let mut p_streamed = Vec::new();
        amplicon_post_process_streaming(
            &p_cfg,
            &amp1,
            &p_vars,
            amplicon_position_ref_stream(&segment).filter(|(p, _)| *p == 100 || *p == 101),
            &splice,
            &mut |av| {
                p_streamed.push(av);
                Ok(())
            },
        )
        .unwrap();
        // Pos 100 (zero coverage) is suppressed; pos 101 (ref_variant) is emitted.
        assert_eq!(p_streamed.len(), 1);
        assert_eq!(p_streamed[0].variant.start_pos, 101);
        assert!(!p_streamed[0].variant.null_record);

        // 4. exact_binary_half_even (5e-5 AF tie) in amplicon_select_vcf_rows_compat & amplicon_post_process:
        let tie_4dp = f64::from_bits(0x3f0a_36e2_eb1c_432d); // 5e-5
        let mk_av_tie = |alt: &str, af: f64| AmpliconVariant {
            variant: mk_good(120, 120, "A", alt, af),
            gdamp: 1,
            tlamp: 1,
            ncamp: 0,
            ampflag: false,
            region: "chr1:100-150".to_string(),
            good_variants: Vec::new(),
            bad_variants: Vec::new(),
        };
        let vj_sel = amplicon_select_vcf_rows_compat(
            vec![mk_av_tie("G", 0.0), mk_av_tie("T", tie_4dp)],
            false,
        );
        assert_eq!(vj_sel.len(), 1);
        assert_eq!(vj_sel[0].variant.var_allele, "G");

        let off_sel = amplicon_select_vcf_rows_compat(
            vec![mk_av_tie("G", 0.0), mk_av_tie("T", tie_4dp)],
            true,
        );
        assert_eq!(off_sel.len(), 1);
        assert_eq!(off_sel[0].variant.var_allele, "T");
    }

    #[test]
    fn test_aggregate_sample_amplicon_vars_primer_masking_highest_cov_ref_and_concordance() {
        let mut config = make_test_config();
        config.compat = crate::config::COMPAT_OFF;
        let splice = rustc_hash::FxHashSet::default();

        // Two overlapping amplicons:
        // amp0: region 90..=210, insert 100..=200 (both ends primed)
        // amp1: region 130..=250, insert 150..=250 (single-primer 5'-primed: insert_end == region_end == 250)
        let amp0 = Region {
            chr: "chr1".to_string(),
            start: Pos1(90),
            end: Pos1(210),
            gene: "G1".to_string(),
            amplicon: Some(AmpliconInfo {
                insert_start: Pos1(100),
                insert_end: Pos1(200),
            }),
        };
        let amp1 = Region {
            chr: "chr1".to_string(),
            start: Pos1(130),
            end: Pos1(250),
            gene: "G1".to_string(),
            amplicon: Some(AmpliconInfo {
                insert_start: Pos1(150),
                insert_end: Pos1(250),
            }),
        };
        let segment = vec![amp0, amp1];

        let mut vars0 = crate::javamap::JavaHashMap::new();
        let mut vars1 = crate::javamap::JavaHashMap::new();

        // 1. Pos 160 (covered by both amp0 & amp1):
        //    - Concordant SNV A>G in both amp0 (AF=0.25, cov=80) and amp1 (AF=0.25, cov=320) -> picks cov=320 rep, gdamp=2, tlamp=2, ncamp=0, ampflag=false.
        //    - ref_variant: amp0 cov=80, amp1 cov=320 -> picks highest-cov ref_variant (320)!
        let mut v160_0 = Vars::new();
        let mut r160_0 = make_good_variant(160, "A", "A", 0.75);
        r160_0.total_pos_cov = 80;
        v160_0.ref_variant = Some(r160_0);
        let mut g160_0 = make_good_variant(160, "A", "G", 0.25);
        g160_0.total_pos_cov = 80;
        g160_0.pos_coverage = 20;
        v160_0.variants.push(g160_0);
        vars0.insert(160, v160_0);

        let mut v160_1 = Vars::new();
        let mut r160_1 = make_good_variant(160, "A", "A", 0.75);
        r160_1.total_pos_cov = 320;
        v160_1.ref_variant = Some(r160_1);
        let mut g160_1 = make_good_variant(160, "A", "G", 0.25);
        g160_1.total_pos_cov = 320;
        g160_1.pos_coverage = 80;
        v160_1.variants.push(g160_1);
        vars1.insert(160, v160_1);

        // 2. Pos 145 (inside amp0 insert 100..=200, but inside amp1's 5' primer 130..=149!):
        //    amp1 has a primer-dimer artifact C>T at 145; amp0 has ONLY ref_variant (cov=200).
        //    Symmetric primer masking must drop amp1's primer variant!
        let mut v145_0 = Vars::new();
        let mut r145_0 = make_good_variant(145, "C", "C", 1.0);
        r145_0.total_pos_cov = 200;
        v145_0.ref_variant = Some(r145_0);
        vars0.insert(145, v145_0);

        let mut v145_1 = Vars::new();
        v145_1.variants.push(make_good_variant(145, "C", "T", 0.50));
        vars1.insert(145, v145_1);

        // 3. Pos 150 (exact `insert_start` of amp1, also inside amp0):
        //    amp1 also has a deletion 148..=151 overlapping amp1's 5' primer (< 150) which MUST be masked on amp1,
        //    while SNV G>A at exact boundary 150 (`== insert_start`) and SNV T>C at exact SPE boundary 250 (`== insert_end == region_end`) MUST be kept!
        let mut v150_1 = Vars::new();
        let mut del_primer = make_good_variant(148, "TTGA", "T", 0.40);
        del_primer.end_pos = 151;
        del_primer.description = "-3".to_string();
        v150_1.variants.push(del_primer);
        v150_1.variants.push(make_good_variant(150, "G", "A", 0.30));
        vars1.insert(150, v150_1);

        let mut v150_0 = Vars::new();
        v150_0.variants.push(make_good_variant(150, "G", "A", 0.28));
        vars0.insert(150, v150_0);

        let mut v250_1 = Vars::new();
        v250_1.variants.push(make_good_variant(250, "T", "C", 0.35));
        vars1.insert(250, v250_1);

        // 4. Pos 170: pure reference site (simulating Normal on StrongSomatic / StrongLOH fallback)
        //    with amp0 cov=500 and amp1 cov=5 (< 500/50 = 10 -> ncamp=1)!
        let mut v170_0 = Vars::new();
        let mut r170_0 = make_good_variant(170, "G", "G", 1.0);
        r170_0.total_pos_cov = 500;
        v170_0.ref_variant = Some(r170_0);
        vars0.insert(170, v170_0);

        let mut v170_1 = Vars::new();
        let mut r170_1 = make_good_variant(170, "G", "G", 1.0);
        r170_1.total_pos_cov = 5;
        v170_1.ref_variant = Some(r170_1);
        vars1.insert(170, v170_1);

        let (aligned, conc) =
            aggregate_sample_amplicon_vars(&config, &segment, vec![vars0, vars1], &splice);

        // Check pos 160: highest-cov ref (320) + highest-cov tie-broken G rep (320) + (2, 2, 0, false)
        let agg160 = aligned.get(&160).expect("pos 160 present");
        assert_eq!(agg160.ref_variant.as_ref().unwrap().total_pos_cov, 320);
        assert_eq!(agg160.variants.len(), 1);
        assert_eq!(agg160.variants[0].total_pos_cov, 320);
        assert_eq!(
            conc.resolve(160, 160, "A", "G", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 2,
                tlamp: 2,
                ncamp: 0,
                ampflag: false,
            }
        );

        // Check pos 145: primer artifact from amp1 masked out completely; only amp0's ref_variant remains!
        let agg145 = aligned.get(&145).expect("pos 145 ref present");
        assert!(agg145.variants.is_empty());
        assert_eq!(agg145.ref_variant.as_ref().unwrap().total_pos_cov, 200);

        // Check pos 150 & 250: primer-overlapping deletion 148..151 masked out, boundary SNVs at 150 and 250 kept!
        let agg150 = aligned.get(&150).expect("pos 150 present");
        assert_eq!(agg150.variants.len(), 1);
        assert_eq!(agg150.variants[0].var_allele, "A");
        assert_eq!(
            conc.resolve(150, 150, "G", "A", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 2,
                tlamp: 2,
                ncamp: 0,
                ampflag: false,
            }
        );
        let agg250 = aligned.get(&250).expect("pos 250 SPE boundary present");
        assert_eq!(agg250.variants.len(), 1);
        assert_eq!(agg250.variants[0].var_allele, "C");
        assert_eq!(
            conc.resolve(250, 250, "T", "C", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 1,
                tlamp: 1,
                ncamp: 0,
                ampflag: false,
            }
        );

        // Check pos 170 pure-ref fallback concordance: highest-cov ref (500) + gdamp=0, tlamp=2, ncamp=1
        let agg170 = aligned.get(&170).expect("pos 170 ref present");
        assert_eq!(agg170.ref_variant.as_ref().unwrap().total_pos_cov, 500);
        assert_eq!(
            conc.resolve(170, 170, "G", "T", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 0,
                tlamp: 2,
                ncamp: 1,
                ampflag: false,
            }
        );
    }

    #[test]
    fn test_aggregate_sample_amplicon_vars_good_rep_beats_higher_af_bad_rep_tumor_and_normal() {
        let mut config = make_test_config();
        config.compat = crate::config::COMPAT_OFF;
        let splice = rustc_hash::FxHashSet::default();

        let amp0 = make_amp_region("chr1", 100, 220);
        let amp1 = make_amp_region("chr1", 140, 260);
        let seg_fwd = vec![amp0.clone(), amp1.clone()];
        let seg_rev = vec![amp1, amp0];

        let mk_ref = |pos: i32, base: &str, cov: i32| -> Variant {
            let mut rv = make_good_variant(pos, base, base, 1.0);
            rv.total_pos_cov = cov;
            rv.pos_coverage = cov;
            rv.ref_fwd_cov = cov / 2;
            rv.ref_rev_cov = cov - rv.ref_fwd_cov;
            rv
        };

        // Good passing call on deep amplicon (35/200 = 0.1750, hicnt = 35 >= min_reads=2)
        let mk_deep_good = |pos: i32, ref_a: &str, alt_a: &str, vd: i32, dp: i32| -> Variant {
            let mut v = make_good_variant(pos, ref_a, alt_a, vd as f64 / dp as f64);
            v.pos_coverage = vd;
            v.total_pos_cov = dp;
            v.hicnt = vd;
            v.hicov = dp;
            v.fwd = vd / 2;
            v.rev = vd - v.fwd;
            v.ref_fwd_cov = (dp - vd) / 2;
            v.ref_rev_cov = (dp - vd) - v.ref_fwd_cov;
            v.bias = "2;2".to_string();
            v.pstd = true;
            v.qstd = true;
            v
        };

        // Failing 1/1 dropout call with higher raw AF = 1.0000 (`hicnt = 1 < min_reads = 2`, `is_good = false`)
        let mk_dropout_bad = |pos: i32, ref_a: &str, alt_a: &str| -> Variant {
            let mut v = make_good_variant(pos, ref_a, alt_a, 1.0);
            v.pos_coverage = 1;
            v.total_pos_cov = 1;
            v.hicnt = 1;
            v.hicov = 1;
            v.fwd = 1;
            v.rev = 0;
            v.bias = "0;1".to_string();
            v
        };

        // 1. Both `[amp0, amp1]` (good first, bad second) and `[amp1, amp0]` (bad first, good second)
        //    must select the `is_good = true` representative (total_pos_cov = 200, hicnt = 35)!
        for (seg, flip) in [(&seg_fwd, false), (&seg_rev, true)] {
            let mut m_good = crate::javamap::JavaHashMap::new();
            let mut vg = Vars::new();
            vg.ref_variant = Some(mk_ref(160, "A", 200));
            vg.variants.push(mk_deep_good(160, "A", "G", 35, 200));
            m_good.insert(160, vg);

            let mut m_bad = crate::javamap::JavaHashMap::new();
            let mut vb = Vars::new();
            vb.variants.push(mk_dropout_bad(160, "A", "G"));
            m_bad.insert(160, vb);

            let list = if flip {
                vec![m_bad, m_good]
            } else {
                vec![m_good, m_bad]
            };
            let (aligned, conc) = aggregate_sample_amplicon_vars(&config, seg, list, &splice);
            let rep = &aligned.get(&160).expect("pos 160 present").variants[0];
            assert_eq!(
                (rep.total_pos_cov, rep.hicnt),
                (200, 35),
                "passing 35/200 representative must beat failing 1/1 AF=1.0 representative (flip={flip})"
            );
            assert!(rep.is_good_var(
                aligned.get(&160).unwrap().ref_variant.as_ref(),
                &rep.var_type_cow(),
                &splice,
                &config
            ));
            assert_eq!(
                conc.resolve(160, 160, "A", "G", seg),
                crate::config::AmpliconConcordance {
                    gdamp: 1,
                    tlamp: 2,
                    ncamp: 1,
                    ampflag: false,
                }
            );
        }

        // 2. End-to-end somatic classification across Tumor (`StrongSomatic`) and Normal (`Germline` & `StrongLOH`)
        //    when `amp0` is deep `is_good = true` and `amp1` has a `1/1` (`AF = 1.0, hicnt = 1`) failing read:
        let mut t0 = crate::javamap::JavaHashMap::new();
        let mut t1 = crate::javamap::JavaHashMap::new();
        let mut n0 = crate::javamap::JavaHashMap::new();
        let mut n1 = crate::javamap::JavaHashMap::new();

        // Pos 160: StrongSomatic in Tumor (35/200 good + 1/1 bad), pure 200x ref in Normal
        let mut t160_0 = Vars::new();
        t160_0.ref_variant = Some(mk_ref(160, "A", 200));
        t160_0.variants.push(mk_deep_good(160, "A", "G", 35, 200));
        t0.insert(160, t160_0);
        let mut t160_1 = Vars::new();
        t160_1.variants.push(mk_dropout_bad(160, "A", "G"));
        t1.insert(160, t160_1);
        let mut n160_0 = Vars::new();
        n160_0.ref_variant = Some(mk_ref(160, "A", 200));
        n0.insert(160, n160_0);

        // Pos 170: Germline in BOTH Tumor and Normal (70/200 good + 1/1 bad in both -> must stay Germline, NOT StrongSomatic!)
        let mut t170_0 = Vars::new();
        t170_0.ref_variant = Some(mk_ref(170, "C", 200));
        t170_0.variants.push(mk_deep_good(170, "C", "T", 70, 200));
        t0.insert(170, t170_0);
        let mut t170_1 = Vars::new();
        t170_1.variants.push(mk_dropout_bad(170, "C", "T"));
        t1.insert(170, t170_1);
        let mut n170_0 = Vars::new();
        n170_0.ref_variant = Some(mk_ref(170, "C", 200));
        n170_0.variants.push(mk_deep_good(170, "C", "T", 70, 200));
        n0.insert(170, n170_0);
        let mut n170_1 = Vars::new();
        n170_1.variants.push(mk_dropout_bad(170, "C", "T"));
        n1.insert(170, n170_1);

        // Pos 180: StrongLOH in Normal (80/200 good + 1/1 bad), pure 200x ref in Tumor
        let mut t180_0 = Vars::new();
        t180_0.ref_variant = Some(mk_ref(180, "G", 200));
        t0.insert(180, t180_0);
        let mut n180_0 = Vars::new();
        n180_0.ref_variant = Some(mk_ref(180, "G", 200));
        n180_0.variants.push(mk_deep_good(180, "G", "A", 80, 200));
        n0.insert(180, n180_0);
        let mut n180_1 = Vars::new();
        n180_1.variants.push(mk_dropout_bad(180, "G", "A"));
        n1.insert(180, n180_1);

        let (aligned_t, _) =
            aggregate_sample_amplicon_vars(&config, &seg_fwd, vec![t0, t1], &splice);
        let (aligned_n, _) =
            aggregate_sample_amplicon_vars(&config, &seg_fwd, vec![n0, n1], &splice);
        let seg_bounds = [(100i32, 260i32)];
        let seg_reg = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(260),
            gene: "G1".to_string(),
            amplicon: None,
        };
        let svs = crate::postprocess::somatic::somatic_post_process_owned_with_bounds(
            &config,
            &seg_reg,
            aligned_t,
            aligned_n,
            &splice,
            None,
            100,
            Some(&seg_bounds),
        );
        let got_statuses: Vec<(i32, &str)> = svs
            .iter()
            .map(|sv| (sv.variant1.start_pos, sv.status.as_str()))
            .collect();
        assert_eq!(
            got_statuses,
            vec![
                (160, "StrongSomatic"),
                (170, "Germline"),
                (180, "StrongLOH"),
            ]
        );
        assert_eq!(svs[1].variant2.as_ref().unwrap().pos_coverage, 70);
    }

    #[test]
    fn test_vcf_padded_deletion_and_anchored_symbolic_sv_amplicon_bounds() {
        let mut config = make_test_config();
        config.compat = crate::config::COMPAT_OFF;
        let splice = rustc_hash::FxHashSet::default();

        // Amplicon with 5' primer 130..=149, insert 150..=250, 3' primer 251..=270
        let amp = Region {
            chr: "chr1".to_string(),
            start: Pos1(130),
            end: Pos1(270),
            gene: "G1".to_string(),
            amplicon: Some(AmpliconInfo {
                insert_start: Pos1(150),
                insert_end: Pos1(250),
            }),
        };
        let segment = vec![amp.clone()];

        // (a) In-insert 3-bp deletion of 150..=152 (`start_pos = 149`, `ref = "CATG"`, `var = "C"`)
        //     is retained with `gdamp = 1, tlamp = 1`, whereas primer-overlapping 149..=151 (`start_pos = 148`)
        //     is masked!
        let mut t_vars = crate::javamap::JavaHashMap::new();
        let mut v150_t = Vars::new();
        let mut r150 = make_good_variant(150, "A", "A", 1.0);
        r150.total_pos_cov = 200;
        v150_t.ref_variant = Some(r150.clone());

        let mut del_in_insert = make_good_variant(149, "CATG", "C", 0.35);
        del_in_insert.end_pos = 152;
        del_in_insert.description = "-3".to_string();
        v150_t.variants.push(del_in_insert);

        let mut del_in_primer = make_good_variant(148, "TCAT", "T", 0.40);
        del_in_primer.end_pos = 151;
        del_in_primer.description = "-3".to_string();
        v150_t.variants.push(del_in_primer);
        t_vars.insert(150, v150_t);

        // (b) 1-bp-upstream-anchored `<DUP>` at `start_pos = 149` (`bp = 150 = insert_start`,
        //     Origin A map key `149`) under `COMPAT_OFF`
        let mut v149_t = Vars::new();
        let mut dup_sv = make_good_variant(149, "C", "<DUP>", 0.30);
        dup_sv.end_pos = 230;
        dup_sv.description = "+<dup-20>".to_string();
        v149_t.variants.push(dup_sv.clone());
        t_vars.insert(149, v149_t);

        let mut n_vars = crate::javamap::JavaHashMap::new();
        let mut v150_n = Vars::new();
        v150_n.ref_variant = Some(r150);
        n_vars.insert(150, v150_n);

        let (aligned_t, conc_t) =
            aggregate_sample_amplicon_vars(&config, &segment, vec![t_vars], &splice);
        let (aligned_n, conc_n) =
            aggregate_sample_amplicon_vars(&config, &segment, vec![n_vars], &splice);

        let kept_150 = &aligned_t.get(&150).expect("pos 150 retained").variants;
        assert_eq!(kept_150.len(), 1, "only in-insert deletion retained");
        assert_eq!(kept_150[0].start_pos, 149);
        assert_eq!(kept_150[0].ref_allele, "CATG");
        assert_eq!(
            conc_t.resolve(149, 152, "CATG", "C", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 1,
                tlamp: 1,
                ncamp: 0,
                ampflag: false,
            }
        );
        assert_eq!(
            conc_n.resolve(149, 152, "CATG", "C", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 0,
                tlamp: 1,
                ncamp: 0,
                ampflag: false,
            }
        );

        let kept_149 = &aligned_t
            .get(&149)
            .expect("anchored <DUP> at 149 retained")
            .variants;
        assert_eq!(kept_149.len(), 1);
        assert_eq!(kept_149[0].var_allele, "<DUP>");
        assert_eq!(
            conc_t.resolve(149, 230, "C", "<DUP>", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 1,
                tlamp: 1,
                ncamp: 0,
                ampflag: false,
            }
        );
        assert_eq!(
            conc_n.resolve(149, 230, "C", "<DUP>", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 0,
                tlamp: 1,
                ncamp: 0,
                ampflag: false,
            }
        );

        let strict_bounds = [(150i32, 250i32)];
        let svs = crate::postprocess::somatic::somatic_post_process_owned_with_bounds(
            &config,
            &amp,
            aligned_t,
            aligned_n,
            &splice,
            None,
            100,
            Some(&strict_bounds),
        );
        assert!(
            svs.iter()
                .any(|sv| sv.variant1.start_pos == 149 && sv.variant1.var_allele == "<DUP>"),
            "paired somatic amplicon strict_region_bounds must retain anchored <DUP> at 149"
        );

        // Single-sample `amplicon_post_process` (`amplicon_retain_strict_bounds`) also retains
        // 1-bp-upstream-anchored `<DUP>` at `start_pos = 149` (`bp = 150 == insert_start`) under `COMPAT_OFF`
        let mut single_map = crate::javamap::JavaHashMap::new();
        let mut v150_single = Vars::new();
        v150_single.variants.push(dup_sv);
        single_map.insert(150, v150_single);
        let mut amps_on_pos: BTreeMap<i32, Vec<(usize, &Region)>> = BTreeMap::new();
        amps_on_pos.insert(150, vec![(0, &amp)]);
        let single_res =
            amplicon_post_process(&config, &amp, &[single_map], &amps_on_pos, &splice).unwrap();
        assert_eq!(single_res.len(), 1);
        assert_eq!(single_res[0].variant.start_pos, 149);
        assert_eq!(single_res[0].variant.var_allele, "<DUP>");
    }

    #[test]
    fn test_wave5_amplicon_adj_complex_vartype_and_upfront_iupac_eq_concordance_parity() {
        let base_cfg = make_test_config();
        let splice = rustc_hash::FxHashSet::default();
        let amp0 = make_amp_region("chr1", 100, 220);
        let amp1 = make_amp_region("chr1", 120, 240);
        let segment = vec![amp0.clone(), amp1.clone()];

        // 1. Switch 29 (`recompute_vartype_after_adj_complex`, `clintval/vartovcf#109`,
        //    `clintval/vartovcf#147`) in `amplicon_post_process` (`:696`) and `aggregate_sample_amplicon_vars` (`:1298`):
        let mk_complex = |pos: i32, desc: &str, ref_a: &str, var_a: &str| {
            let mut v = make_good_variant(pos, ref_a, var_a, 0.35);
            v.description = desc.to_string();
            v.end_pos = pos + (ref_a.len() as i32) - 1;
            v
        };
        let mut c0 = crate::javamap::JavaHashMap::new();
        let mut c1 = crate::javamap::JavaHashMap::new();
        for (pos, desc, r, a) in [
            (130, "#GT", "GA", "GT"),
            (140, "+TG#AA", "CAA", "CAATG"),
            (150, "-2#AA", "CAATG", "CAA"),
        ] {
            let mut vs0 = Vars::new();
            vs0.variants.push(mk_complex(pos, desc, r, a));
            c0.insert(pos, vs0);
            let mut vs1 = Vars::new();
            vs1.variants.push(mk_complex(pos, desc, r, a));
            c1.insert(pos, vs1);
        }
        let mut amps_on_pos = BTreeMap::new();
        for pos in [130, 140, 150] {
            amps_on_pos.insert(pos, vec![(0, amp0.clone()), (1, amp1.clone())]);
        }
        let vj_amp = amplicon_post_process(
            &base_cfg,
            &amp0,
            &[c0.clone(), c1.clone()],
            &amps_on_pos,
            &splice,
        )
        .expect("vj amplicon_post_process");
        assert!(vj_amp.iter().all(|av| av.variant.var_type == "Complex"));

        let mut off_cfg = base_cfg.clone();
        off_cfg.compat = crate::config::COMPAT_OFF;
        let off_amp = amplicon_post_process(
            &off_cfg,
            &amp0,
            &[c0.clone(), c1.clone()],
            &amps_on_pos,
            &splice,
        )
        .expect("off amplicon_post_process");
        let got_vt: Vec<(i32, &str, &str, &str)> = off_amp
            .iter()
            .map(|av| {
                (
                    av.variant.start_pos,
                    av.variant.ref_allele.as_str(),
                    av.variant.var_allele.as_str(),
                    av.variant.var_type.as_str(),
                )
            })
            .collect();
        assert_eq!(
            got_vt,
            vec![
                (131, "A", "T", "SNV"),
                (142, "A", "ATG", "Insertion"),
                (152, "ATG", "A", "Deletion"),
            ]
        );

        let (agg_c, conc_c) =
            aggregate_sample_amplicon_vars(&off_cfg, &segment, vec![c0, c1], &splice);
        assert_eq!(
            conc_c.resolve(131, 131, "A", "T", &segment),
            crate::config::AmpliconConcordance {
                gdamp: 2,
                tlamp: 2,
                ncamp: 0,
                ampflag: false,
            }
        );
        let n_empty = crate::javamap::JavaHashMap::new();
        let som_c = crate::postprocess::somatic::somatic_post_process_owned(
            &off_cfg, &amp0, agg_c, n_empty, &splice, None, 100,
        );
        assert_eq!(
            som_c
                .iter()
                .map(|s| (s.variant1.start_pos, s.variant1.var_type.as_str()))
                .collect::<Vec<_>>(),
            vec![(131, "SNV"), (142, "Insertion"), (152, "Deletion")]
        );

        // 2. Switch 26 (`sanitize_iupac_alt_alleles`) upfront `REF == ALT` skip in BOTH
        //    `amplicon_post_process` (`:480`) and `aggregate_sample_amplicon_vars` (`:1165`):
        //    - Pos 160: `amp0` has sanitized `"A" -> "A"` (`AF = 0.75`, `DP = 200`) + `"A" -> "G"` (`AF = 0.25`, `DP = 200`);
        //      `amp1` has `"A" -> "G"` (`AF = 0.25`, `DP = 200`).
        //    - Pos 170: `amp0` has ONLY sanitized `"C" -> "C"` (`AF = 1.0`, `DP = 200`, `ref_variant = None`);
        //      `amp1` has `"C" -> "T"` (`AF = 0.30`, `DP = 200`).
        let build_iupac_amp_maps = || {
            let mut m0 = crate::javamap::JavaHashMap::new();
            let mut m1 = crate::javamap::JavaHashMap::new();

            let mut v160_0 = Vars::new();
            let mut eq_160 = make_good_variant(160, "A", "A", 0.75);
            eq_160.total_pos_cov = 200;
            v160_0.variants.push(eq_160);
            let mut g160_0 = make_good_variant(160, "A", "G", 0.25);
            g160_0.total_pos_cov = 200;
            v160_0.variants.push(g160_0);
            m0.insert(160, v160_0);

            let mut v160_1 = Vars::new();
            let mut g160_1 = make_good_variant(160, "A", "G", 0.25);
            g160_1.total_pos_cov = 200;
            v160_1.variants.push(g160_1);
            m1.insert(160, v160_1);

            let mut v170_0 = Vars::new();
            let mut eq_170 = make_good_variant(170, "C", "C", 1.0);
            eq_170.total_pos_cov = 200;
            v170_0.variants.push(eq_170);
            m0.insert(170, v170_0);

            let mut v170_1 = Vars::new();
            let mut t170_1 = make_good_variant(170, "C", "T", 0.30);
            t170_1.total_pos_cov = 200;
            v170_1.variants.push(t170_1);
            m1.insert(170, v170_1);

            (m0, m1)
        };

        let mut iupac_amps_on_pos = BTreeMap::new();
        iupac_amps_on_pos.insert(160, vec![(0, amp0.clone()), (1, amp1.clone())]);
        iupac_amps_on_pos.insert(170, vec![(0, amp0.clone()), (1, amp1.clone())]);

        // Default parity: `"A" -> "A"` and `"C" -> "C"` leak into `gvs` and trigger false `AMPBIAS`!
        let (m0_vj, m1_vj) = build_iupac_amp_maps();
        let vj_iupac = amplicon_post_process(
            &base_cfg,
            &amp0,
            &[m0_vj, m1_vj],
            &iupac_amps_on_pos,
            &splice,
        )
        .expect("vj_iupac");
        assert!(vj_iupac.len() > 2);
        assert!(vj_iupac.iter().any(|av| av.ampflag));

        // With `sanitize_iupac_alt_alleles: true` (isolated and `COMPAT_OFF`):
        let mut iso_iupac_cfg = base_cfg.clone();
        iso_iupac_cfg.compat.sanitize_iupac_alt_alleles = true;
        for cfg in [&iso_iupac_cfg, &off_cfg] {
            let (m0_on, m1_on) = build_iupac_amp_maps();
            let single_on = amplicon_post_process(
                cfg,
                &amp0,
                &[m0_on.clone(), m1_on.clone()],
                &iupac_amps_on_pos,
                &splice,
            )
            .expect("single_on");
            assert_eq!(single_on.len(), 2);
            assert_eq!(single_on[0].variant.start_pos, 160);
            assert_eq!(single_on[0].variant.var_allele, "G");
            assert_eq!(single_on[0].gdamp, 2);
            assert_eq!(single_on[0].tlamp, 2);
            assert_eq!(single_on[0].ncamp, 0);
            assert!(!single_on[0].ampflag);

            assert_eq!(single_on[1].variant.start_pos, 170);
            assert_eq!(single_on[1].variant.var_allele, "T");
            assert_eq!(single_on[1].gdamp, 1);
            assert_eq!(single_on[1].tlamp, 2);
            assert_eq!(single_on[1].ncamp, 0);
            assert!(!single_on[1].ampflag);

            let (agg_on, conc_on) =
                aggregate_sample_amplicon_vars(cfg, &segment, vec![m0_on, m1_on], &splice);
            assert_eq!(agg_on.get(&160).expect("160").variants.len(), 1);
            assert_eq!(agg_on.get(&160).expect("160").variants[0].var_allele, "G");
            assert_eq!(
                conc_on.resolve(160, 160, "A", "G", &segment),
                crate::config::AmpliconConcordance {
                    gdamp: 2,
                    tlamp: 2,
                    ncamp: 0,
                    ampflag: false,
                }
            );
            assert_eq!(agg_on.get(&170).expect("170").variants.len(), 1);
            assert_eq!(agg_on.get(&170).expect("170").variants[0].var_allele, "T");
            assert_eq!(
                conc_on.resolve(170, 170, "C", "T", &segment),
                crate::config::AmpliconConcordance {
                    gdamp: 1,
                    tlamp: 2,
                    ncamp: 0,
                    ampflag: false,
                }
            );

            // Exact `10,000x` pure-ref (`variants = []`) vs `10,000x` skipped `REF == ALT` parity
            // alongside a `100x` variant-carrying amplicon (`AmpliconPostProcessModule.java:94-95`):
            // `amp0` (`0/10,000` variant reads) must NOT inflate `maxcov` to mark variant-carrying
            // `amp1` (`100x`) as `NCAMP = 1` (which would falsely excuse `1 < 2 - 1` `AMPBIAS`).
            let build_skewed_pair = |use_skipped_eq: bool| {
                let mut s0 = crate::javamap::JavaHashMap::new();
                let mut s1 = crate::javamap::JavaHashMap::new();
                let mut vs0 = Vars::new();
                if use_skipped_eq {
                    let mut eq_hi = make_good_variant(170, "C", "C", 1.0);
                    eq_hi.total_pos_cov = 10_000;
                    vs0.variants.push(eq_hi);
                } else {
                    let mut rv_hi = Variant::new();
                    rv_hi.start_pos = 170;
                    rv_hi.end_pos = 170;
                    rv_hi.ref_allele = "C".to_string();
                    rv_hi.var_allele = "C".to_string();
                    rv_hi.total_pos_cov = 10_000;
                    vs0.ref_variant = Some(rv_hi);
                }
                s0.insert(170, vs0);
                let mut vs1 = Vars::new();
                let mut v_lo = make_good_variant(170, "C", "T", 0.50);
                v_lo.total_pos_cov = 100;
                vs1.variants.push(v_lo);
                s1.insert(170, vs1);
                (s0, s1)
            };
            let mut skew_pos = BTreeMap::new();
            skew_pos.insert(170, vec![(0, amp0.clone()), (1, amp1.clone())]);
            let (clean0, clean1) = build_skewed_pair(false);
            let (iupac0, iupac1) = build_skewed_pair(true);
            let single_clean = amplicon_post_process(
                cfg,
                &amp0,
                &[clean0.clone(), clean1.clone()],
                &skew_pos,
                &splice,
            )
            .expect("single_clean");
            let single_iupac = amplicon_post_process(
                cfg,
                &amp0,
                &[iupac0.clone(), iupac1.clone()],
                &skew_pos,
                &splice,
            )
            .expect("single_iupac");
            assert_eq!(single_clean[0].ncamp, 0);
            assert_eq!(single_iupac[0].ncamp, single_clean[0].ncamp);
            let (_, conc_clean) =
                aggregate_sample_amplicon_vars(cfg, &segment, vec![clean0, clean1], &splice);
            let (_, conc_iupac) =
                aggregate_sample_amplicon_vars(cfg, &segment, vec![iupac0, iupac1], &splice);
            assert_eq!(
                conc_iupac.resolve(170, 170, "C", "T", &segment),
                conc_clean.resolve(170, 170, "C", "T", &segment)
            );
            assert_eq!(conc_iupac.resolve(170, 170, "C", "T", &segment).ncamp, 0);

            // End-to-end paired somatic amplicon pure-`eq` reference hydration invariants:
            // (a) `StrongSomatic` at 170 when Normal has ONLY sanitized `"C" -> "C"` (`fwd = 85, rev = 95, ref_fwd = 10, ref_rev = 10, DP = 200`):
            let mut n_pure0 = crate::javamap::JavaHashMap::new();
            let mut vn170 = Vars::new();
            let mut eq_n170 = make_good_variant(170, "C", "C", 0.90);
            eq_n170.total_pos_cov = 200;
            eq_n170.pos_coverage = 180;
            eq_n170.fwd = 85;
            eq_n170.rev = 95;
            eq_n170.ref_fwd_cov = 10;
            eq_n170.ref_rev_cov = 10;
            vn170.variants.push(eq_n170);
            n_pure0.insert(170, vn170);
            let (agg_n_pure, _) = aggregate_sample_amplicon_vars(
                cfg,
                &segment,
                vec![n_pure0, crate::javamap::JavaHashMap::new()],
                &splice,
            );
            let som_ss = crate::postprocess::somatic::somatic_post_process_owned(
                cfg,
                &amp0,
                agg_on.clone(),
                agg_n_pure,
                &splice,
                None,
                100,
            );
            let ss_170 = som_ss
                .iter()
                .find(|s| s.variant1.start_pos == 170)
                .expect("StrongSomatic at 170");
            assert_eq!(ss_170.status, "StrongSomatic");
            let n_var = ss_170.variant2.as_ref().expect("normal ref placeholder");
            assert_eq!(n_var.pos_coverage, 0);
            assert_eq!(n_var.fwd, 0);
            assert_eq!(n_var.rev, 0);
            assert_eq!(n_var.ref_fwd_cov + n_var.ref_rev_cov, 200);
            assert_eq!(n_var.total_pos_cov, 200);
            assert_eq!(n_var.genotype.as_deref(), Some("C/C"));
            assert_eq!(n_var.bias, "2;0");

            // (b) `StrongLOH` at 170 when Tumor has ONLY sanitized `"C" -> "C"` (`fwd = 90, rev = 110, DP = 200`) and Normal has `"C" -> "T"`:
            let mut t_pure0 = crate::javamap::JavaHashMap::new();
            let mut vt170 = Vars::new();
            let mut eq_t170 = make_good_variant(170, "C", "C", 1.0);
            eq_t170.total_pos_cov = 200;
            eq_t170.pos_coverage = 200;
            eq_t170.fwd = 90;
            eq_t170.rev = 110;
            vt170.variants.push(eq_t170);
            t_pure0.insert(170, vt170);
            let (agg_t_pure, _) = aggregate_sample_amplicon_vars(
                cfg,
                &segment,
                vec![t_pure0, crate::javamap::JavaHashMap::new()],
                &splice,
            );
            let som_loh = crate::postprocess::somatic::somatic_post_process_owned(
                cfg, &amp0, agg_t_pure, agg_on, &splice, None, 100,
            );
            let loh_170 = som_loh
                .iter()
                .find(|s| s.variant1.start_pos == 170)
                .expect("StrongLOH at 170");
            assert_eq!(loh_170.status, "StrongLOH");
            assert_eq!(loh_170.variant1.pos_coverage, 0);
            assert_eq!(loh_170.variant1.fwd, 0);
            assert_eq!(loh_170.variant1.rev, 0);
            assert_eq!(
                loh_170.variant1.ref_fwd_cov + loh_170.variant1.ref_rev_cov,
                200
            );
            assert_eq!(loh_170.variant1.total_pos_cov, 200);
            assert_eq!(loh_170.variant1.genotype.as_deref(), Some("C/C"));
            assert_eq!(loh_170.variant1.bias, "2;0");
        }
    }

    #[test]
    fn test_wave6_switch35_switch38_and_ruling3_amplicon_postprocess() {
        let amp_reg = make_amp_region("chr1", 90, 150);
        let splice = rustc_hash::FxHashSet::default();

        let mut m0 = crate::javamap::JavaHashMap::new();
        let mut v100 = Vars::new();
        let mut d100 = make_good_variant(100, "CAT", "C", 0.40);
        d100.end_pos = 102;
        d100.description = "-2AT".to_string();
        d100.pos_coverage = 8;
        d100.total_pos_cov = 20;
        d100.left_seq = "GGGGGG".to_string();
        v100.variants.push(d100);
        m0.insert(100, v100);

        let mut v102 = Vars::new();
        let mut d102 = make_good_variant(102, "TAT", "T", 0.45);
        d102.end_pos = 104;
        d102.description = "-2AT".to_string();
        d102.pos_coverage = 3;
        d102.total_pos_cov = 20;
        d102.left_seq = "GGGGGGCA".to_string();
        v102.variants.push(d102);
        m0.insert(102, v102);

        let mut v110 = Vars::new();
        let mut c110 = make_good_variant(110, "CA", "CG", 0.35);
        c110.end_pos = 111;
        c110.description = "C&G".to_string();
        c110.shift3 = 3;
        c110.msi = 5.0;
        c110.msint = 1;
        v110.variants.push(c110);
        m0.insert(110, v110);

        let vars_per_amp = vec![m0];
        let mut amps_on_pos = BTreeMap::new();
        amps_on_pos.insert(100, vec![(0usize, &amp_reg)]);
        amps_on_pos.insert(102, vec![(0usize, &amp_reg)]);
        amps_on_pos.insert(110, vec![(0usize, &amp_reg)]);

        let mut cfg = Config::default();
        cfg.emit_intermediate = true;
        cfg.freq = 0.01;
        cfg.min_reads = 2;
        cfg.read_pos = 5.0;
        cfg.base_qual = 22.5;
        cfg.qratio = 1.5;
        cfg.compat.canonical_vcf_complex_left_align = true;
        cfg.compat.indel_only_shift3 = true;

        let out = amplicon_post_process(&cfg, &amp_reg, &vars_per_amp, &amps_on_pos, &splice)
            .expect("amplicon_post_process");
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].variant.start_pos, 100);
        assert_eq!(out[0].variant.ref_allele, "CAT");
        assert_eq!(out[0].variant.var_allele, "C");
        assert_eq!(out[0].variant.pos_coverage, 8);
        assert_eq!(out[1].variant.start_pos, 111);
        assert_eq!(out[1].variant.ref_allele, "A");
        assert_eq!(out[1].variant.var_allele, "G");
        assert_eq!(out[1].variant.shift3, 0);
    }
}
