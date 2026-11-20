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
                    for tv in variants {
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
                                if variant.is_good_var(
                                    vtmp.ref_variant.as_ref(),
                                    "",
                                    splice,
                                    config,
                                ) {
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

            // Apply complex variant adjustment — clone max_var once at emission
            let mut out_variant = max_var.clone();
            out_variant.var_type = out_variant.var_type();
            if out_variant.var_type == "Complex" {
                out_variant.adj_complex();
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
    let mut buckets: BTreeMap<i32, AmpliconVariant> = BTreeMap::new();
    for av in results {
        let v_af = crate::types::round_half_even(av.variant.frequency, 4);
        match buckets.entry(av.variant.start_pos) {
            std::collections::btree_map::Entry::Vacant(e) => {
                e.insert(av);
            }
            std::collections::btree_map::Entry::Occupied(mut e) => {
                let cur_af = crate::types::round_half_even(e.get().variant.frequency, 4);
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
}
