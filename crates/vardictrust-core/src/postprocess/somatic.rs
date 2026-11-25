//! Somatic post-processing module — determines somatic STATUS for paired samples.
//!
//! Ports `SomaticPostProcessModule.java` from VarDictJava v1.8.3.
//! Compares variants from two BAM files (tumor=BAM1, normal=BAM2) and
//! determines the somatic STATUS (StrongSomatic, LikelySomatic, Germline,
//! LikelyLOH, AFDiff, SampleSpecific, Deletion, StrongLOH).

use crate::config::Config;
use crate::javamap::JavaHashMap;
use crate::regions::Region;
use crate::types::VardictError;
use crate::variations::{Variant, Vars};

/// Callback that runs the full pipeline (pileup→realign→SV→tovars) on both BAMs
/// combined over a given region, returning aligned variants keyed by position
/// while updating the shared `splice` set and `max_read_length` in place
/// (`SomaticPostProcessModule.java:198-200, 303-311, 398-406`).
#[deprecated(since = "1.0.10", note = "use PositionPipelineFn")]
pub type PipelineFn<'a> = dyn Fn(
        &Region,
        &mut rustc_hash::FxHashSet<String>,
        &mut i32,
    ) -> Result<JavaHashMap<i32, Vars>, VardictError>
    + 'a;

/// Single-position combine callback returning finalized [`Vars`] only at `target_pos`
/// (the outer `aligned_variants` pileup map key, `SomaticPostProcessModule.java:198,307,408`).
pub type PositionPipelineFn<'a> = dyn Fn(
        &Region,
        &mut rustc_hash::FxHashSet<String>,
        &mut i32,
        i32,
    ) -> Result<Option<Vars>, VardictError>
    + 'a;

// STATUS constants matching VarDictJava's SomaticPostProcessModule
const STRONG_SOMATIC: &str = "StrongSomatic";

/// Ports `Patterns.MINUS_NUM_NUM` = `-\d\d` (unanchored find).
fn minus_num_num_find(desc: &str) -> bool {
    desc.as_bytes()
        .windows(3)
        .any(|w| matches!(w, [b'-', d1, d2] if d1.is_ascii_digit() && d2.is_ascii_digit()))
}
const SAMPLE_SPECIFIC: &str = "SampleSpecific";
const DELETION: &str = "Deletion";
const LIKELY_LOH: &str = "LikelyLOH";
const GERMLINE: &str = "Germline";
const STRONG_LOH: &str = "StrongLOH";
const LIKELY_SOMATIC: &str = "LikelySomatic";
const AF_DIFF: &str = "AFDiff";
const COMPLEX: &str = "Complex";
const SNV: &str = "SNV";

/// Output of somatic post-processing: a paired variant with its STATUS.
#[derive(Debug, Clone)]
pub struct SomaticVariant {
    /// Variant from BAM1 (tumor).
    pub variant1: Variant,
    /// Variant from BAM2 (normal). None if sample-specific (no coverage in BAM2).
    pub variant2: Option<Variant>,
    /// Somatic STATUS string.
    pub status: String,
    /// SV info from BAM1.
    pub sv1: String,
    /// SV info from BAM2.
    pub sv2: String,
    /// Whether Java's endVariant (source of left/right seq, shift3, msi, msint)
    /// is variant2 (true) or variant1 (false).
    pub end_is_variant2: bool,
    /// Amplicon concordance for BAM1 (tumor) when paired somatic amplicon mode is active (`AstraZeneca-NGS/VarDictJava#378`).
    pub tumor_amplicon: Option<crate::config::AmpliconConcordance>,
    /// Amplicon concordance for BAM2 (normal) when paired somatic amplicon mode is active (`AstraZeneca-NGS/VarDictJava#378`).
    pub normal_amplicon: Option<crate::config::AmpliconConcordance>,
}

fn collect_sorted_region_positions(
    region: &Region,
    aligned_variants_bam1: &JavaHashMap<i32, Vars>,
    aligned_variants_bam2: &JavaHashMap<i32, Vars>,
) -> Vec<i32> {
    let mut all_positions: Vec<i32> = aligned_variants_bam1
        .keys()
        .chain(aligned_variants_bam2.keys())
        .copied()
        .filter(|&p| u32::try_from(p).is_ok_and(|pu| pu >= region.start.0 && pu <= region.end.0))
        .collect();
    all_positions.sort_unstable();
    all_positions.dedup();
    all_positions
}

fn dispatch_somatic_position(
    config: &Config,
    position: i32,
    v1: Option<std::borrow::Cow<'_, Vars>>,
    v2: Option<std::borrow::Cow<'_, Vars>>,
    chr: &str,
    splice: &mut std::borrow::Cow<'_, rustc_hash::FxHashSet<String>>,
    results: &mut Vec<SomaticVariant>,
    pipeline: Option<&PositionPipelineFn<'_>>,
    max_read_length: &mut i32,
) {
    match (v1, v2) {
        (None, None) => {}
        (None, Some(v2)) => {
            // No coverage for sample 1 → Deletion
            calling_for_one_sample(config, v2, true, DELETION, splice, results);
        }
        (Some(v1), None) => {
            // No coverage for sample 2 → SampleSpecific
            calling_for_one_sample(config, v1, false, SAMPLE_SPECIFIC, splice, results);
        }
        (Some(v1), Some(v2)) => {
            // Both samples have coverage
            calling_for_both_samples(
                config,
                position,
                v1,
                v2,
                chr,
                splice,
                results,
                pipeline,
                max_read_length,
            );
        }
    }
}

#[inline]
fn single_somatic_region_bounds_i32(region: &Region) -> [(i32, i32); 1] {
    [(
        i32::try_from(region.start.0.max(1)).unwrap_or(i32::MAX),
        i32::try_from(region.end.0).unwrap_or(i32::MAX),
    )]
}

/// Performs somatic analysis of variants from both samples (BAM1 and BAM2)
/// using a single-position combine callback [`PositionPipelineFn`].
///
/// Merges all positions from both samples' aligned variants, then for each
/// position in the region determines the somatic STATUS.
///
/// Ports `SomaticPostProcessModule.accept()`.
pub fn somatic_post_process_at(
    config: &Config,
    region: &Region,
    aligned_variants_bam1: &JavaHashMap<i32, Vars>,
    aligned_variants_bam2: &JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
    pipeline: Option<&PositionPipelineFn<'_>>,
    max_read_length: i32,
) -> Vec<SomaticVariant> {
    let mut results = Vec::new();
    let mut splice: std::borrow::Cow<'_, rustc_hash::FxHashSet<String>> =
        std::borrow::Cow::Borrowed(splice);
    let mut max_read_length = max_read_length;

    let all_positions =
        collect_sorted_region_positions(region, aligned_variants_bam1, aligned_variants_bam2);

    for position in all_positions {
        let v1 = aligned_variants_bam1
            .get(&position)
            .map(std::borrow::Cow::Borrowed);
        let v2 = aligned_variants_bam2
            .get(&position)
            .map(std::borrow::Cow::Borrowed);
        dispatch_somatic_position(
            config,
            position,
            v1,
            v2,
            &region.chr,
            &mut splice,
            &mut results,
            pipeline,
            &mut max_read_length,
        );
    }

    if config.compat.strict_region_bounds {
        let single_bounds = single_somatic_region_bounds_i32(region);
        results.retain(|sv| crate::interval_set_contains(&single_bounds, sv.variant1.start_pos));
    }

    results
}

/// Owned fast-path counterpart of [`somatic_post_process_at`], extracting per-position
/// [`Vars`] via `get_mut(&position).map(std::mem::take)` to move [`Variant`] and SV
/// strings without cloning on hot somatic paths.
#[cfg(test)]
pub(crate) fn somatic_post_process_owned(
    config: &Config,
    region: &Region,
    aligned_variants_bam1: JavaHashMap<i32, Vars>,
    aligned_variants_bam2: JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
    pipeline: Option<&PositionPipelineFn<'_>>,
    max_read_length: i32,
) -> Vec<SomaticVariant> {
    let single_bounds = single_somatic_region_bounds_i32(region);
    somatic_post_process_owned_with_bounds(
        config,
        region,
        aligned_variants_bam1,
        aligned_variants_bam2,
        splice,
        pipeline,
        max_read_length,
        config
            .compat
            .strict_region_bounds
            .then_some(single_bounds.as_slice()),
    )
}

/// Owned fast-path counterpart of [`somatic_post_process_at`] using pre-merged chromosome
/// target intervals when `strict_bounds` is supplied (`#369`, `AR-4`).
pub(crate) fn somatic_post_process_owned_with_bounds(
    config: &Config,
    region: &Region,
    mut aligned_variants_bam1: JavaHashMap<i32, Vars>,
    mut aligned_variants_bam2: JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
    pipeline: Option<&PositionPipelineFn<'_>>,
    max_read_length: i32,
    strict_bounds: Option<&[(i32, i32)]>,
) -> Vec<SomaticVariant> {
    let mut results = Vec::new();
    let mut splice: std::borrow::Cow<'_, rustc_hash::FxHashSet<String>> =
        std::borrow::Cow::Borrowed(splice);
    let mut max_read_length = max_read_length;

    let all_positions =
        collect_sorted_region_positions(region, &aligned_variants_bam1, &aligned_variants_bam2);

    for position in all_positions {
        let v1 = aligned_variants_bam1
            .get_mut(&position)
            .map(|v| std::borrow::Cow::Owned(std::mem::take(v)));
        let v2 = aligned_variants_bam2
            .get_mut(&position)
            .map(|v| std::borrow::Cow::Owned(std::mem::take(v)));
        dispatch_somatic_position(
            config,
            position,
            v1,
            v2,
            &region.chr,
            &mut splice,
            &mut results,
            pipeline,
            &mut max_read_length,
        );
    }

    if let Some(merged_bounds) = strict_bounds {
        results.retain(|sv| crate::interval_set_contains(merged_bounds, sv.variant1.start_pos));
    } else if config.compat.strict_region_bounds {
        results.retain(|sv| sv.variant1.start_pos >= 1);
    }

    results
}

/// Performs somatic analysis of variants from both samples (BAM1 and BAM2).
///
/// Wraps [`PipelineFn`] as a [`PositionPipelineFn`] and delegates to
/// [`somatic_post_process_at`].
#[deprecated(
    since = "1.0.10",
    note = "use postprocess::somatic::somatic_post_process_at"
)]
#[allow(deprecated)]
pub fn somatic_post_process(
    config: &Config,
    region: &Region,
    aligned_variants_bam1: &JavaHashMap<i32, Vars>,
    aligned_variants_bam2: &JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
    pipeline: Option<&PipelineFn<'_>>,
    max_read_length: i32,
) -> Vec<SomaticVariant> {
    if let Some(pipe) = pipeline {
        let adapted = |comb_region: &Region,
                       splice_set: &mut rustc_hash::FxHashSet<String>,
                       max_rl: &mut i32,
                       target_pos: i32|
         -> Result<Option<Vars>, VardictError> {
            let mut vars_map = pipe(comb_region, splice_set, max_rl)?;
            Ok(vars_map.remove(&target_pos))
        };
        somatic_post_process_at(
            config,
            region,
            aligned_variants_bam1,
            aligned_variants_bam2,
            splice,
            Some(&adapted),
            max_read_length,
        )
    } else {
        somatic_post_process_at(
            config,
            region,
            aligned_variants_bam1,
            aligned_variants_bam2,
            splice,
            None,
            max_read_length,
        )
    }
}

/// Implements variations analysis from one sample.
///
/// Ports `SomaticPostProcessModule.callingForOneSample()`.
fn calling_for_one_sample(
    config: &Config,
    variants: std::borrow::Cow<'_, Vars>,
    is_first_cover: bool,
    var_label: &str,
    splice: &rustc_hash::FxHashSet<String>,
    results: &mut Vec<SomaticVariant>,
) {
    if variants.variants.is_empty() {
        return;
    }
    let (ref_var_cow, var_iter, sv_cow) = super::unpack_cow_vars(variants);

    let start_idx = results.len();
    for variant_cow in var_iter {
        if variant_cow.ref_allele == variant_cow.var_allele {
            continue;
        }
        let vtype = variant_cow.var_type_cow();
        if !variant_cow.is_good_var(ref_var_cow.as_deref(), &vtype, splice, config) {
            continue;
        }
        let mut v = variant_cow.into_owned();
        if v.var_type != vtype.as_ref() {
            v.var_type.clear();
            v.var_type.push_str(&vtype);
        }
        if v.var_type == COMPLEX {
            v.adj_complex();
        }

        if is_first_cover {
            // Only BAM2 has coverage → Deletion. Java passes tumorVariant=null,
            // so the tumor columns are all zero; header fields come from the variant.
            let mut tumor = crate::variations::Variant::new();
            tumor.start_pos = v.start_pos;
            tumor.end_pos = v.end_pos;
            tumor.ref_allele = v.ref_allele.clone();
            tumor.var_allele = v.var_allele.clone();
            tumor.var_type = v.var_type.clone();
            tumor.shift3 = v.shift3;
            tumor.msi = v.msi;
            tumor.msint = v.msint;
            tumor.left_seq = v.left_seq.clone();
            tumor.right_seq = v.right_seq.clone();
            tumor.debug.clone_from(&v.debug);
            results.push(SomaticVariant {
                variant1: tumor,
                variant2: Some(v),
                status: var_label.to_string(),
                sv1: String::new(),
                sv2: String::new(),
                end_is_variant2: true,
                tumor_amplicon: None,
                normal_amplicon: None,
            });
        } else {
            // Only BAM1 has coverage → SampleSpecific
            results.push(SomaticVariant {
                variant1: v,
                variant2: None,
                status: var_label.to_string(),
                sv1: String::new(),
                sv2: String::new(),
                end_is_variant2: false,
                tumor_amplicon: None,
                normal_amplicon: None,
            });
        }
    }

    if let Some((last, prior)) = results
        .get_mut(start_idx..)
        .and_then(|s| s.split_last_mut())
    {
        if is_first_cover {
            for r in prior {
                r.sv2.clear();
                r.sv2.push_str(&sv_cow);
            }
            last.sv2 = sv_cow.into_owned();
        } else {
            for r in prior {
                r.sv1.clear();
                r.sv1.push_str(&sv_cow);
            }
            last.sv1 = sv_cow.into_owned();
        }
    }
}

/// Run analysis of variations from BAM1 and BAM2.
///
/// Ports `SomaticPostProcessModule.callingForBothSamples()`.
fn calling_for_both_samples(
    config: &Config,
    position: i32,
    v1: std::borrow::Cow<'_, Vars>,
    v2: std::borrow::Cow<'_, Vars>,
    chr: &str,
    splice: &mut std::borrow::Cow<'_, rustc_hash::FxHashSet<String>>,
    results: &mut Vec<SomaticVariant>,
    pipeline: Option<&PositionPipelineFn<'_>>,
    max_read_length: &mut i32,
) {
    if v1.variants.is_empty() && v2.variants.is_empty() {
        return;
    }
    if !v1.variants.is_empty() {
        print_variations_from_first_sample(
            config,
            position,
            v1,
            v2,
            chr,
            splice,
            results,
            pipeline,
            max_read_length,
        );
    } else if !v2.variants.is_empty() {
        // sample 1 has only reference
        print_variations_from_second_sample(
            config,
            v1,
            v2,
            chr,
            position,
            splice,
            results,
            pipeline,
            max_read_length,
        );
    }
}

/// Analyse and output variations from BAM1 based on variations from BAM2.
///
/// Ports `SomaticPostProcessModule.printVariationsFromFirstSample()`.
fn print_variations_from_first_sample(
    config: &Config,
    position: i32,
    v1: std::borrow::Cow<'_, Vars>,
    mut v2: std::borrow::Cow<'_, Vars>,
    chr: &str,
    splice: &mut std::borrow::Cow<'_, rustc_hash::FxHashSet<String>>,
    results: &mut Vec<SomaticVariant>,
    pipeline: Option<&PositionPipelineFn<'_>>,
    max_read_length: &mut i32,
) {
    // Keep `v1` unconsumed until `v1.variants.first()` passes `is_good_var` so the
    // `number_of_processed_variation == 0` LOH fallback preserves `v1.variant_by_desc`,
    // `v1.variants.first()`, `v1.ref_variant`, and `v1.sv`.
    let first_passes = v1.variants.first().is_some_and(|first_v| {
        let vt = first_v.var_type_cow();
        first_v.is_good_var(v1.ref_variant.as_ref(), &vt, splice, config)
    });

    if !first_passes {
        if v2.variants.is_empty() {
            return;
        }
        let (v2_ref_cow, v2_iter, sv2_cow) = super::unpack_cow_vars(v2);
        let start_idx = results.len();
        for v2var_cow in v2_iter {
            let vtype = v2var_cow.var_type_cow();
            if !v2var_cow.is_good_var(v2_ref_cow.as_deref(), &vtype, splice, config) {
                continue;
            }
            let mut v2var = v2var_cow.into_owned();
            if v2var.var_type != vtype.as_ref() {
                v2var.var_type.clear();
                v2var.var_type.push_str(&vtype);
            }
            // potential LOH
            let v1nt = v1.variant_by_desc(&v2var.description);

            if let Some(v1nt) = v1nt {
                if v1nt.ref_allele == v1nt.var_allele {
                    continue;
                }
                let type_str = if v1nt.frequency < config.lofreq {
                    LIKELY_LOH
                } else {
                    GERMLINE
                };
                let mut v1nt_out = v1nt.clone();
                if v2var.var_type == COMPLEX {
                    v1nt_out.adj_complex();
                }
                v1nt_out.var_type = v1nt_out.var_type();
                results.push(SomaticVariant {
                    variant1: v1nt_out,
                    variant2: Some(v2var),
                    status: type_str.to_string(),
                    sv1: v1.sv.clone(),
                    sv2: String::new(),
                    end_is_variant2: true,
                    tumor_amplicon: None,
                    normal_amplicon: None,
                });
            } else {
                if v2var.ref_allele == v2var.var_allele {
                    continue;
                }
                // Build a placeholder variant1 from sample1's reference info
                let mut var_for_print =
                    build_var_for_print_from_v1(&v1, config.compat.reconcile_vd_dp);
                let mut v2var_out = v2var;
                if v2var_out.var_type == COMPLEX {
                    v2var_out.adj_complex();
                }
                // Shared columns (ref/var alleles, type) come from normal variant
                var_for_print.ref_allele = v2var_out.ref_allele.clone();
                var_for_print.var_allele = v2var_out.var_allele.clone();
                var_for_print.var_type = v2var_out.var_type.clone();
                var_for_print.start_pos = v2var_out.start_pos;
                var_for_print.end_pos = v2var_out.end_pos;
                var_for_print.shift3 = v2var_out.shift3;
                var_for_print.msi = v2var_out.msi;
                var_for_print.msint = v2var_out.msint;
                var_for_print.left_seq = v2var_out.left_seq.clone();
                var_for_print.right_seq = v2var_out.right_seq.clone();
                var_for_print.debug.clone_from(&v2var_out.debug);
                results.push(SomaticVariant {
                    variant1: var_for_print,
                    variant2: Some(v2var_out),
                    status: STRONG_LOH.to_string(),
                    sv1: String::new(),
                    sv2: String::new(),
                    end_is_variant2: true,
                    tumor_amplicon: None,
                    normal_amplicon: None,
                });
            }
        }
        if let Some((last, prior)) = results
            .get_mut(start_idx..)
            .and_then(|s| s.split_last_mut())
        {
            for r in prior {
                r.sv2.clear();
                r.sv2.push_str(&sv2_cow);
            }
            last.sv2 = sv2_cow.into_owned();
        }
        return;
    }

    // `first_passes` is true: `number_of_processed_variation` will reach >= 1, so the
    // LOH fallback block never runs and `v1` can be unpacked directly.
    // Ports SomaticPostProcessModule.java:172-186 — in Java, `v2.variants` and
    // `v2.varDescriptionStringToVariants` hold references to the same `Variant`
    // objects. When `determinateType` calls `v2nt.isNoise()`, it mutates `v2nt`
    // in place (`totalPosCoverage -= positionCoverage`, etc.), so a subsequent
    // tumor variant at the same position that falls into the `v2nt == null`
    // branch reads the mutated `v2.variants.get(0).totalPosCoverage`.
    // Keep `v2` via `Cow` (`to_mut()` mutates in place on `Cow::Owned` and clones lazily
    // on `Cow::Borrowed` only when a subsequent tumor variant exists at this locus, #45).
    let (v1_ref_cow, v1_iter, sv1_cow) = super::unpack_cow_vars(v1);
    let v1_len = v1_iter.len();

    let start_idx = results.len();
    for (idx, vref_cow) in v1_iter.enumerate() {
        let vartype_check = vref_cow.var_type_cow();
        if idx > 0 && !vref_cow.is_good_var(v1_ref_cow.as_deref(), &vartype_check, splice, config) {
            break;
        }
        if vref_cow.ref_allele == vref_cow.var_allele {
            continue;
        }

        let mut vref = vref_cow.into_owned();
        if vref.var_type != vartype_check.as_ref() {
            vref.var_type.clear();
            vref.var_type.push_str(&vartype_check);
        }
        if vref.var_type == COMPLEX {
            vref.adj_complex();
        }
        let has_more_v1 = idx + 1 < v1_len;
        let nt = vref.description.as_str();

        // Single-allele Cow::Owned fast-path on v2: when no further v1 variants follow,
        // v2 is Owned, var_desc_to_variant is empty, and v2.variants holds solely `nt`,
        // pop `v2`'s Variant directly without cloning.
        let popped_v2nt = if !has_more_v1 {
            match v2 {
                std::borrow::Cow::Owned(ref mut v2_owned)
                    if v2_owned.var_desc_to_variant.is_empty()
                        && v2_owned.variants.len() == 1
                        && v2_owned
                            .variants
                            .first()
                            .is_some_and(|v| v.description == nt) =>
                {
                    v2_owned.variants.pop()
                }
                _ => None,
            }
        } else {
            None
        };

        let v2nt = match popped_v2nt {
            Some(v) => Some(v),
            None => v2.variant_by_desc(nt).cloned(),
        };

        if let Some(mut v2nt_var) = v2nt {
            let type_str = determinate_type(
                config,
                v2.ref_variant.as_ref(),
                &vref,
                &mut v2nt_var,
                splice,
            );
            if has_more_v1 {
                v2.to_mut().update_variant_by_desc(nt, v2nt_var.clone());
            }
            results.push(SomaticVariant {
                variant1: vref,
                variant2: Some(v2nt_var),
                status: type_str,
                sv1: String::new(),
                sv2: String::new(),
                end_is_variant2: true,
                tumor_amplicon: None,
                normal_amplicon: None,
            });
        } else {
            // sample 1 only, should be strong somatic (SomaticPostProcessModule.java:177)
            let mut type_str: &'static str = STRONG_SOMATIC;
            let mut combined_v2nt: Option<Variant> = None;

            // combineAnalysis gate (:193-197): not SNV, description longer than 10
            // or matching MINUS_NUM_NUM (-\d\d), low pos coverage, and not an SV
            // ('<') → re-run the pipeline on both BAMs combined
            if vref.var_type != SNV
                && (nt.len() > 10 || minus_num_num_find(nt))
                && i64::from(vref.pos_coverage) < i64::from(config.min_reads) + 3
                && !nt.contains('<')
            {
                if let Some(pipeline_fn) = pipeline {
                    match combine_analysis(
                        config,
                        pipeline_fn,
                        &vref,
                        chr,
                        position,
                        splice,
                        max_read_length,
                    ) {
                        CombineVerdict::Skip => {
                            // FALSE → skip this variant (:202-204)
                            continue;
                        }
                        CombineVerdict::Germline(v) => {
                            type_str = GERMLINE;
                            combined_v2nt = Some(*v);
                        }
                        CombineVerdict::Keep => {} // "" keeps StrongSomatic
                    }
                }
            }

            if type_str == STRONG_SOMATIC {
                let var_for_print = match v2 {
                    std::borrow::Cow::Owned(ref mut v2_owned)
                        if !has_more_v1 && v2_owned.variants.is_empty() =>
                    {
                        let mut rv = v2_owned.ref_variant.take();
                        if config.compat.reconcile_vd_dp {
                            if let Some(ref mut vp) = rv {
                                vp.ref_fwd_cov = vp.ref_fwd_cov.max(0);
                                vp.ref_rev_cov = vp.ref_rev_cov.max(0);
                                vp.total_pos_cov = vp
                                    .total_pos_cov
                                    .max(0)
                                    .max(vp.ref_fwd_cov.saturating_add(vp.ref_rev_cov));
                            }
                        }
                        rv
                    }
                    _ => build_var_for_print_from_v2(&v2, config.compat.reconcile_vd_dp),
                };
                // SomaticOutputVariant(vref, vref, vref, varForPrint, …) (:212)
                results.push(SomaticVariant {
                    variant1: vref,
                    variant2: var_for_print,
                    status: type_str.to_string(),
                    sv1: String::new(),
                    sv2: String::new(),
                    end_is_variant2: false,
                    tumor_amplicon: None,
                    normal_amplicon: None,
                });
            } else {
                // SomaticOutputVariant(vref, vref, vref, v2nt, …) (:215): header and
                // end fields stay with vref; the var2 columns carry the
                // combine-adjusted placeholder.
                results.push(SomaticVariant {
                    variant1: vref,
                    variant2: combined_v2nt,
                    status: type_str.to_string(),
                    sv1: String::new(),
                    sv2: String::new(),
                    end_is_variant2: false,
                    tumor_amplicon: None,
                    normal_amplicon: None,
                });
            }
        }
    }

    let sv2_cow: std::borrow::Cow<'_, str> = match v2 {
        std::borrow::Cow::Borrowed(vars) => std::borrow::Cow::Borrowed(vars.sv.as_str()),
        std::borrow::Cow::Owned(vars) => std::borrow::Cow::Owned(vars.sv),
    };
    if let Some((last, prior)) = results
        .get_mut(start_idx..)
        .and_then(|s| s.split_last_mut())
    {
        for r in prior {
            r.sv1.clear();
            r.sv1.push_str(&sv1_cow);
            r.sv2.clear();
            r.sv2.push_str(&sv2_cow);
        }
        last.sv1 = sv1_cow.into_owned();
        last.sv2 = sv2_cow.into_owned();
    }
}

/// Analyse and output variations from BAM2 based on variations from BAM1.
///
/// Ports `SomaticPostProcessModule.printVariationsFromSecondSample()`
/// (SomaticPostProcessModule.java:285-336).
fn print_variations_from_second_sample(
    config: &Config,
    v1: std::borrow::Cow<'_, Vars>,
    v2: std::borrow::Cow<'_, Vars>,
    chr: &str,
    position: i32,
    splice: &mut std::borrow::Cow<'_, rustc_hash::FxHashSet<String>>,
    results: &mut Vec<SomaticVariant>,
    pipeline: Option<&PositionPipelineFn<'_>>,
    max_read_length: &mut i32,
) {
    let mut v1_ref_cow: Option<std::borrow::Cow<'_, Variant>> = match v1 {
        std::borrow::Cow::Borrowed(vars) => {
            vars.ref_variant.as_ref().map(std::borrow::Cow::Borrowed)
        }
        std::borrow::Cow::Owned(Vars { ref_variant, .. }) => {
            ref_variant.map(std::borrow::Cow::Owned)
        }
    };
    let (v2_ref_cow, mut v2_iter, sv2_cow) = super::unpack_cow_vars(v2);

    let start_idx = results.len();
    while let Some(v2var_cow) = v2_iter.next() {
        if v2var_cow.ref_allele == v2var_cow.var_allele {
            continue;
        }
        let vtype = v2var_cow.var_type_cow();
        if !v2var_cow.is_good_var(v2_ref_cow.as_deref(), &vtype, splice, config) {
            continue;
        }
        let mut v2var = v2var_cow.into_owned();
        if v2var.var_type != vtype.as_ref() {
            v2var.var_type.clear();
            v2var.var_type.push_str(&vtype);
        }

        // combineAnalysis gate (:301-302): low normal coverage, not an SV, long
        // description or a -\d\d pattern (MINUS_NUM_NUM, Patterns.java:18)
        let description = &v2var.description;
        let should_combine = i64::from(v2var.pos_coverage) < i64::from(config.min_reads) + 3
            && !description.contains('<')
            && (description.len() > 10 || minus_num_num_find(description));
        let mut type_str: &'static str = STRONG_LOH;
        let mut tumor_from_combine: Option<Variant> = None;
        if should_combine {
            if let Some(pipeline_fn) = pipeline {
                match combine_analysis(
                    config,
                    pipeline_fn,
                    &v2var,
                    chr,
                    position,
                    splice,
                    max_read_length,
                ) {
                    CombineVerdict::Skip => continue, // :313-315
                    CombineVerdict::Germline(v) => {
                        type_str = GERMLINE; // :318-319
                        tumor_from_combine = Some(*v); // varForPrint = v1nt (:320)
                    }
                    CombineVerdict::Keep => {} // "" keeps StrongLOH
                }
            }
        }

        if v2var.var_type == COMPLEX {
            v2var.adj_complex(); // :329-331, before output construction
        }

        let tumor_var = match tumor_from_combine {
            Some(mut t) => {
                // Java's SomaticOutputVariant(v2var, v2var, varForPrint, v2var, …)
                // (:333): header fields come from the ADJUSTED v2var while the
                // var1 columns carry the combine-adjusted placeholder. The Rust
                // printer takes header + var1 columns from variant1, so overlay
                // Java's adjusted-placeholder columns onto the adjusted v2var.
                t.start_pos = v2var.start_pos;
                t.end_pos = v2var.end_pos;
                t.ref_allele = v2var.ref_allele.clone();
                t.var_allele = v2var.var_allele.clone();
                t.var_type = v2var.var_type.clone();
                t.shift3 = v2var.shift3;
                t.msi = v2var.msi;
                t.msint = v2var.msint;
                t.left_seq = v2var.left_seq.clone();
                t.right_seq = v2var.right_seq.clone();
                t.debug.clone_from(&v2var.debug);
                t.duprate = 0.0; // fresh Variant in Java; never set by combineAnalysis
                t
            }
            None => {
                // Fallback: StrongLOH with reference placeholder (:322-327).
                // Header/allele fields are re-sourced from v2var exactly as the
                // SomaticOutputVariant constructor does.
                let mut var_for_print = if v2_iter.is_empty() {
                    match v1_ref_cow.take() {
                        Some(v1ref) => v1ref.into_owned(),
                        None => Variant::new(),
                    }
                } else {
                    match v1_ref_cow.as_deref() {
                        Some(v1ref) => v1ref.clone(),
                        None => Variant::new(),
                    }
                };
                if config.compat.reconcile_vd_dp {
                    var_for_print.ref_fwd_cov = var_for_print.ref_fwd_cov.max(0);
                    var_for_print.ref_rev_cov = var_for_print.ref_rev_cov.max(0);
                    var_for_print.total_pos_cov = var_for_print.total_pos_cov.max(0).max(
                        var_for_print
                            .ref_fwd_cov
                            .saturating_add(var_for_print.ref_rev_cov),
                    );
                }
                var_for_print.ref_allele = v2var.ref_allele.clone();
                var_for_print.var_allele = v2var.var_allele.clone();
                var_for_print.var_type = v2var.var_type.clone();
                var_for_print.start_pos = v2var.start_pos;
                var_for_print.end_pos = v2var.end_pos;
                var_for_print.shift3 = v2var.shift3;
                var_for_print.msi = v2var.msi;
                var_for_print.msint = v2var.msint;
                var_for_print.left_seq = v2var.left_seq.clone();
                var_for_print.right_seq = v2var.right_seq.clone();
                var_for_print.debug.clone_from(&v2var.debug);
                var_for_print
            }
        };

        results.push(SomaticVariant {
            variant1: tumor_var,
            variant2: Some(v2var),
            status: type_str.to_string(),
            sv1: String::new(),
            sv2: String::new(),
            end_is_variant2: true,
            tumor_amplicon: None,
            normal_amplicon: None,
        });
    }

    if let Some((last, prior)) = results
        .get_mut(start_idx..)
        .and_then(|s| s.split_last_mut())
    {
        for r in prior {
            r.sv2.clear();
            r.sv2.push_str(&sv2_cow);
        }
        last.sv2 = sv2_cow.into_owned();
    }
}

/// Verdict of `combineAnalysis` — ports the `CombineAnalysisData.type` string.
enum CombineVerdict {
    /// Java `"FALSE"`: the combined run contradicts the variant — caller skips it.
    Skip,
    /// Java `""`: keep the current status (StrongSomatic / StrongLOH fallback).
    Keep,
    /// Java `"Germline"`: variant2 was filled by subtraction (combined − variant1).
    Germline(Box<Variant>),
}

/// Re-run the pipeline on both BAMs combined to decide whether a sample-specific
/// variant is real: combined evidence ≥ minr over the base sample → Germline,
/// contradictory/absent → FALSE, otherwise keep.
///
/// Ports `SomaticPostProcessModule.combineAnalysis()` (SomaticPostProcessModule.java:385-476).
/// `variant1` is the variant being tested (tumor side in printVariationsFromFirstSample,
/// normal side in printVariationsFromSecondSample); `position` is the pileup iteration
/// position key the caller is processing — the combined result is looked up by that key
/// and the description string (Java getVarMaybe(vars, position, varn, descriptionString),
/// :408), NOT by the variant's own start_pos.
fn combine_analysis(
    config: &Config,
    pipeline: &PositionPipelineFn<'_>,
    variant1: &Variant,
    chr: &str,
    position: i32,
    splice: &mut std::borrow::Cow<'_, rustc_hash::FxHashSet<String>>,
    max_read_length: &mut i32,
) -> CombineVerdict {
    if config.verbose {
        crate::verbose_eprintln!(
            "Start Combine {} {}",
            position,
            crate::util::SanitizedBytes(variant1.description.as_bytes())
        );
    }
    // Don't do it for structural variants (:394)
    if i64::from(variant1.end_pos) - i64::from(variant1.start_pos) > i64::from(config.sv_min_len) {
        return CombineVerdict::Keep;
    }

    // Region spans the variant ± maxReadLength (:398)
    let region = Region {
        chr: chr.to_string(),
        start: crate::types::Pos1(variant1.start_pos.saturating_sub(*max_read_length).max(1) as u32),
        end: crate::types::Pos1(variant1.end_pos.saturating_add(*max_read_length).max(1) as u32),
        gene: String::new(),
        amplicon: None,
    };

    // Java: a pipeline exception escapes combineAnalysis to accept()'s per-position
    // try/catch (:89-91) → printExceptionAndContinue prints it and processing
    // continues from the next position; abort only after >10 continued exceptions
    // (Utils.java:261-278, Configuration.java:300). The port keeps the FALSE-skip
    // control flow but emits the error instead of swallowing it.
    // Always pass the outer `aligned_variants` pileup map key `position` (never
    // `variant1.start_pos`, which `ToVarsBuilder` shifts on deletions/`+AT&GC`/`-3`).
    let combined_vars = match pipeline(&region, splice.to_mut(), max_read_length, position) {
        Ok(vars) => vars,
        Err(e) => {
            if config.verbose {
                crate::verbose_eprintln!(
                    "There was Exception while processing position on {}. The processing will be continued from the next position.: {}",
                    position, e
                );
            }
            return CombineVerdict::Skip;
        }
    };

    // Look up by ITERATION position key and description in the map (:408).
    let Some(combined_var) = combined_vars
        .as_ref()
        .and_then(|pv| pv.variant_by_desc(variant1.description.as_str()))
    else {
        // vref == null → FALSE (:475)
        return CombineVerdict::Skip;
    };

    let base_pc = variant1.pos_coverage;
    let comb_pc = combined_var.pos_coverage;
    let min_reads_bias = i32::try_from(config.min_reads_bias).unwrap_or(i32::MAX);

    if i64::from(comb_pc) - i64::from(base_pc) >= i64::from(config.min_reads) {
        // Enough extra evidence in the other sample: fill variant2 by subtraction
        // (combined − variant1), clamped at zero (:414-436)
        let mut var2 = Variant::new();
        var2.total_pos_cov = combined_var
            .total_pos_cov
            .saturating_sub(variant1.total_pos_cov)
            .max(0);
        var2.pos_coverage = comb_pc.saturating_sub(base_pc).max(0);
        var2.ref_fwd_cov = combined_var
            .ref_fwd_cov
            .saturating_sub(variant1.ref_fwd_cov)
            .max(0);
        var2.ref_rev_cov = combined_var
            .ref_rev_cov
            .saturating_sub(variant1.ref_rev_cov)
            .max(0);
        var2.fwd = combined_var.fwd.saturating_sub(variant1.fwd).max(0);
        var2.rev = combined_var.rev.saturating_sub(variant1.rev).max(0);

        if var2.pos_coverage != 0 {
            // Weighted mean subtraction formula (:439-444)
            let pcv = var2.pos_coverage as f64;
            let ccv = comb_pc as f64;
            let bcv = base_pc as f64;
            var2.mean_position =
                (combined_var.mean_position * ccv - variant1.mean_position * bcv) / pcv;
            var2.mean_quality =
                (combined_var.mean_quality * ccv - variant1.mean_quality * bcv) / pcv;
            var2.mean_mapping_quality = (combined_var.mean_mapping_quality * ccv
                - variant1.mean_mapping_quality * bcv)
                / pcv;
            var2.hi_freq = (combined_var.hi_freq * ccv - variant1.hi_freq * bcv) / pcv;
            var2.extra_freq = (combined_var.extra_freq * ccv - variant1.extra_freq * bcv) / pcv;
            var2.nm = (combined_var.nm * ccv - variant1.nm * bcv) / pcv;
        } else {
            var2.mean_position = 0.0;
            var2.mean_quality = 0.0;
            var2.mean_mapping_quality = 0.0;
            var2.hi_freq = 0.0;
            var2.extra_freq = 0.0;
            var2.nm = 0.0;
        }
        var2.pstd = true; // :453 isAtLeastAt2Positions
        var2.qstd = true; // :454 hasAtLeast2DiffQualities

        if config.compat.reconcile_vd_dp {
            var2.total_pos_cov = var2.total_pos_cov.max(var2.pos_coverage);
            var2.hicnt = combined_var.hicnt.saturating_sub(variant1.hicnt).max(0);
            var2.hicov = combined_var
                .hicov
                .saturating_sub(variant1.hicov)
                .max(var2.hicnt);
            var2.hi_freq = if var2.hi_freq.is_finite() {
                var2.hi_freq.clamp(0.0, 1.0)
            } else {
                0.0
            };
            var2.extra_freq = if var2.extra_freq.is_finite() {
                var2.extra_freq.clamp(0.0, 1.0)
            } else {
                0.0
            };
        }

        if var2.total_pos_cov <= 0 {
            return CombineVerdict::Skip; // :456-458
        }

        var2.frequency = if config.compat.reconcile_vd_dp {
            (var2.pos_coverage as f64 / var2.total_pos_cov as f64).clamp(0.0, 1.0)
        } else {
            var2.pos_coverage as f64 / var2.total_pos_cov as f64
        }; // :460
           // Can't back-calculate; copy from variant1 (:461)
        var2.hq_lq_ratio = variant1.hq_lq_ratio;
        var2.genotype = combined_var.genotype.clone(); // :462
                                                       // Recalculate strand bias from subtracted counts (:463-464); Java strandBias
                                                       // uses conf.bias (default 0.05) and conf.minBiasReads (-B, default 2)
        let ref_bias = crate::variations::util::strand_bias(
            var2.ref_fwd_cov,
            var2.ref_rev_cov,
            0.05,
            min_reads_bias,
        );
        let var_bias =
            crate::variations::util::strand_bias(var2.fwd, var2.rev, 0.05, min_reads_bias);
        var2.bias.clear();
        let _ =
            std::fmt::Write::write_fmt(&mut var2.bias, format_args!("{};{}", ref_bias, var_bias));
        CombineVerdict::Germline(Box::new(var2))
    } else if i64::from(comb_pc) < i64::from(base_pc) - 2 {
        // Combined run produced less than the base sample alone (:466-470)
        if config.verbose {
            // Java config.y (-y); vardict.pl:657 prints the same line.
            crate::verbose_eprintln!(
                "Combine produce less: {} {} {} {} {}",
                crate::util::SanitizedBytes(chr.as_bytes()),
                position,
                crate::util::SanitizedBytes(variant1.description.as_bytes()),
                comb_pc,
                base_pc
            );
        }
        CombineVerdict::Skip
    } else {
        CombineVerdict::Keep
    }
}

/// Analyse two variations and return their somatic STATUS type.
///
/// Ports `SomaticPostProcessModule.determinateType()`.
fn determinate_type(
    config: &Config,
    v2_ref_variant: Option<&Variant>,
    standard_variant: &Variant,
    variant_to_compare: &mut Variant,
    splice: &rustc_hash::FxHashSet<String>,
) -> String {
    let vtc = variant_to_compare;
    let vartype = &standard_variant.var_type;

    let type_str = if vtc.is_good_var(v2_ref_variant, vartype, splice, config) {
        if standard_variant.frequency > (1.0 - config.lofreq)
            && vtc.frequency < 0.8
            && vtc.frequency > 0.2
        {
            LIKELY_LOH
        } else if vtc.frequency < config.lofreq || vtc.pos_coverage <= 1 {
            LIKELY_SOMATIC
        } else {
            GERMLINE
        }
    } else if vtc.frequency < config.lofreq || vtc.pos_coverage <= 1 {
        LIKELY_SOMATIC
    } else {
        AF_DIFF
    };

    // Override: if variant is noise and type is SNV → StrongSomatic.
    // is_noise mutates vtc (zeroes coverage), matching Java's isNoise side effect
    // which propagates to the variant used in the output.
    if vtc.is_noise(config.base_qual, config.lofreq) && *vartype == SNV {
        return STRONG_SOMATIC.to_string();
    }

    type_str.to_string()
}

/// Build a placeholder Variant from v2's first variant or reference for StrongSomatic output.
///
/// Ports the `varForPrint` construction in `printVariationsFromFirstSample` when
/// the variant is not found in sample2.
fn build_var_for_print_from_v2(v2: &Vars, reconcile_vd_dp: bool) -> Option<Variant> {
    let mut out = if let Some(v2r) = v2.variants.first() {
        let mut var_for_print = Variant::new();
        var_for_print.total_pos_cov = if v2r.total_pos_cov != 0 {
            v2r.total_pos_cov
        } else {
            0
        };
        var_for_print.ref_fwd_cov = if v2r.ref_fwd_cov != 0 {
            v2r.ref_fwd_cov
        } else {
            0
        };
        var_for_print.ref_rev_cov = if v2r.ref_rev_cov != 0 {
            v2r.ref_rev_cov
        } else {
            0
        };
        Some(var_for_print)
    } else {
        v2.ref_variant.clone()
    };
    if reconcile_vd_dp {
        if let Some(ref mut vp) = out {
            vp.ref_fwd_cov = vp.ref_fwd_cov.max(0);
            vp.ref_rev_cov = vp.ref_rev_cov.max(0);
            vp.total_pos_cov = vp
                .total_pos_cov
                .max(0)
                .max(vp.ref_fwd_cov.saturating_add(vp.ref_rev_cov));
        }
    }
    out
}

/// Build a placeholder Variant from v1's first variant or reference for StrongLOH output.
///
/// Ports the `varForPrint` construction in `printVariationsFromFirstSample` when
/// numberOfProcessedVariation == 0 and variant is not found in sample1.
fn build_var_for_print_from_v1(v1: &Vars, reconcile_vd_dp: bool) -> Variant {
    let v1_first = v1.variants.first();
    let tcov = match v1_first {
        Some(v1var) if v1var.total_pos_cov != 0 => v1var.total_pos_cov,
        _ => 0,
    };

    let (fwd, rev, genotype) = if let Some(ref v1ref) = v1.ref_variant {
        (
            v1ref.fwd,
            v1ref.rev,
            if let Some(v1var) = v1_first {
                v1var.genotype.clone()
            } else {
                Some(format!("{}/{}", v1ref.description, v1ref.description))
            },
        )
    } else {
        (
            0,
            0,
            if let Some(v1var) = v1_first {
                v1var.genotype.clone()
            } else {
                Some("N/N".to_string())
            },
        )
    };

    let mut var_for_print = Variant::new();
    var_for_print.total_pos_cov = tcov;
    var_for_print.ref_fwd_cov = fwd;
    var_for_print.ref_rev_cov = rev;
    if reconcile_vd_dp {
        var_for_print.ref_fwd_cov = var_for_print.ref_fwd_cov.max(0);
        var_for_print.ref_rev_cov = var_for_print.ref_rev_cov.max(0);
        var_for_print.total_pos_cov = var_for_print.total_pos_cov.max(0).max(
            var_for_print
                .ref_fwd_cov
                .saturating_add(var_for_print.ref_rev_cov),
        );
    }
    var_for_print.genotype = genotype;
    var_for_print
}

#[cfg(test)]
#[allow(deprecated)]
mod tests {
    use super::*;
    use crate::types::Pos1;

    fn test_config() -> Config {
        let mut config = Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;
        config.lofreq = 0.05;
        config
    }

    fn test_region() -> Region {
        Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "TEST".to_string(),
            amplicon: None,
        }
    }

    fn make_good_variant(desc: &str, ref_allele: &str, var_allele: &str, freq: f64) -> Variant {
        let mut v = Variant::new();
        v.description = desc.to_string();
        v.ref_allele = ref_allele.to_string();
        v.var_allele = var_allele.to_string();
        v.frequency = freq;
        v.hicnt = 5;
        v.mean_position = 10.0;
        v.mean_quality = 30.0;
        v.hq_lq_ratio = 2.0;
        v.mean_mapping_quality = 30.0;
        v.pos_coverage = 10;
        v.total_pos_cov = 100;
        v.start_pos = 150;
        v.end_pos = 150;
        v
    }

    #[test]
    fn test_sample_specific() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars = Vars::new();
        vars.variants.push(make_good_variant("G", "A", "G", 0.10));
        vars.sv = String::new();
        bam1.insert(150, vars);

        let bam2 = crate::javamap::JavaHashMap::new();

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, SAMPLE_SPECIFIC);
        assert!(results[0].variant2.is_none());
    }

    #[test]
    fn test_deletion() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let bam1 = crate::javamap::JavaHashMap::new();

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars = Vars::new();
        let mut v2_del = make_good_variant("G", "A", "G", 0.10);
        v2_del.debug = "BAM2_VAR_DBG".to_string();
        vars.variants.push(v2_del);
        vars.sv = String::new();
        bam2.insert(150, vars);

        // 1) Deletion (`callingForOneSample` `:121`): beginVariant = v2var -> DEBUG = v2var.DEBUG
        let results = somatic_post_process_at(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, DELETION);
        assert!(results[0].variant2.is_some());
        assert_eq!(results[0].variant1.debug, "BAM2_VAR_DBG");
        let results_owned =
            somatic_post_process_owned(&config, &region, bam1, bam2, &splice, None, 150);
        assert_eq!(results_owned[0].variant1.debug, "BAM2_VAR_DBG");

        let mk_v1_ref = || {
            let mut r = Variant::new();
            r.description = "A".to_string();
            r.ref_allele = "A".to_string();
            r.var_allele = "A".to_string();
            r.total_pos_cov = 80;
            r.ref_fwd_cov = 40;
            r.ref_rev_cov = 40;
            r.fwd = 40;
            r.rev = 40;
            r.debug = "BAM1_REF_DBG".to_string();
            r
        };

        // 2) `print_variations_from_first_sample` StrongLOH (`:270`): beginVariant = v2var
        let mut bam1_b2 = crate::javamap::JavaHashMap::new();
        let mut vars1_b2 = Vars::new();
        let mut v1_bad = make_good_variant("C", "A", "C", 0.002); // fails freq -> !first_passes
        v1_bad.debug = "BAM1_BAD_DBG".to_string();
        vars1_b2.variants.push(v1_bad);
        vars1_b2.ref_variant = Some(mk_v1_ref());
        bam1_b2.insert(150, vars1_b2);

        let mut bam2_b2 = crate::javamap::JavaHashMap::new();
        let mut vars2_b2 = Vars::new();
        let mut v2_loh1 = make_good_variant("G", "A", "G", 0.40);
        v2_loh1.debug = "BAM2_VAR_DBG".to_string();
        vars2_b2.variants.push(v2_loh1);
        bam2_b2.insert(150, vars2_b2);

        let res_b2 =
            somatic_post_process_at(&config, &region, &bam1_b2, &bam2_b2, &splice, None, 150);
        assert_eq!(res_b2.len(), 1);
        assert_eq!(res_b2[0].status, STRONG_LOH);
        assert_eq!(res_b2[0].variant1.debug, "BAM2_VAR_DBG");
        let res_b2_owned =
            somatic_post_process_owned(&config, &region, bam1_b2, bam2_b2, &splice, None, 150);
        assert_eq!(res_b2_owned[0].variant1.debug, "BAM2_VAR_DBG");

        // 3) `print_variations_from_second_sample` CombineVerdict::Germline (`:333`): beginVariant = v2var
        let mut bam1_b3 = crate::javamap::JavaHashMap::new();
        let mut vars1_b3 = Vars::new();
        vars1_b3.ref_variant = Some(mk_v1_ref());
        bam1_b3.insert(150, vars1_b3);

        let mut bam2_b3 = crate::javamap::JavaHashMap::new();
        let mut vars2_b3 = Vars::new();
        let mut v2_ins = make_insertion_variant();
        v2_ins.debug = "BAM2_VAR_DBG".to_string();
        vars2_b3.variants.push(v2_ins);
        bam2_b3.insert(150, vars2_b3);

        let germ_pipe = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 7, 20);
        let germ_adapted = adapt_pipeline(&germ_pipe);
        let res_b3 = somatic_post_process_at(
            &config,
            &region,
            &bam1_b3,
            &bam2_b3,
            &splice,
            Some(&germ_adapted),
            150,
        );
        assert_eq!(res_b3.len(), 1);
        assert_eq!(res_b3[0].status, GERMLINE);
        assert_eq!(res_b3[0].variant1.debug, "BAM2_VAR_DBG");
        let res_b3_owned = somatic_post_process_owned(
            &config,
            &region,
            bam1_b3,
            bam2_b3,
            &splice,
            Some(&germ_adapted),
            150,
        );
        assert_eq!(res_b3_owned[0].variant1.debug, "BAM2_VAR_DBG");

        // 4) `print_variations_from_second_sample` None StrongLOH (`:333`): beginVariant = v2var
        let mut bam1_b4 = crate::javamap::JavaHashMap::new();
        let mut vars1_b4 = Vars::new();
        vars1_b4.ref_variant = Some(mk_v1_ref());
        bam1_b4.insert(150, vars1_b4);

        let mut bam2_b4 = crate::javamap::JavaHashMap::new();
        let mut vars2_b4 = Vars::new();
        let mut v2_loh2 = make_good_variant("G", "A", "G", 0.40);
        v2_loh2.debug = "BAM2_VAR_DBG".to_string();
        vars2_b4.variants.push(v2_loh2);
        bam2_b4.insert(150, vars2_b4);

        let res_b4 =
            somatic_post_process_at(&config, &region, &bam1_b4, &bam2_b4, &splice, None, 150);
        assert_eq!(res_b4.len(), 1);
        assert_eq!(res_b4[0].status, STRONG_LOH);
        assert_eq!(res_b4[0].variant1.debug, "BAM2_VAR_DBG");
        let res_b4_owned =
            somatic_post_process_owned(&config, &region, bam1_b4, bam2_b4, &splice, None, 150);
        assert_eq!(res_b4_owned[0].variant1.debug, "BAM2_VAR_DBG");
    }

    #[test]
    fn test_germline() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        vars1.variants.push(make_good_variant("G", "A", "G", 0.40));
        vars1
            .var_desc_to_variant
            .insert("G".to_string(), make_good_variant("G", "A", "G", 0.40));
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        let v2 = make_good_variant("G", "A", "G", 0.40);
        vars2.variants.push(v2.clone());
        vars2.var_desc_to_variant.insert("G".to_string(), v2);
        bam2.insert(150, vars2);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, GERMLINE);
    }

    #[test]
    fn test_likely_somatic() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        vars1.variants.push(make_good_variant("G", "A", "G", 0.30));
        vars1
            .var_desc_to_variant
            .insert("G".to_string(), make_good_variant("G", "A", "G", 0.30));
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        // Low frequency in normal → LikelySomatic
        let mut v2 = make_good_variant("G", "A", "G", 0.02);
        v2.frequency = 0.02;
        v2.pos_coverage = 2;
        vars2.variants.push(v2.clone());
        vars2.var_desc_to_variant.insert("G".to_string(), v2);
        bam2.insert(150, vars2);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, LIKELY_SOMATIC);
    }

    #[test]
    fn test_strong_somatic_not_in_normal() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        vars1.variants.push(make_good_variant("G", "A", "G", 0.30));
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        // Different variant in normal — "G" not found in v2
        let v2 = make_good_variant("T", "A", "T", 0.30);
        vars2.variants.push(v2.clone());
        vars2.var_desc_to_variant.insert("T".to_string(), v2);
        bam2.insert(150, vars2);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, STRONG_SOMATIC);
    }

    #[test]
    fn test_strong_loh() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        // BAM1 has no non-ref variants, BAM2 has a good variant
        let mut bam1 = crate::javamap::JavaHashMap::new();
        let vars1 = Vars::new(); // empty variants
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        vars2.variants.push(make_good_variant("G", "A", "G", 0.40));
        bam2.insert(150, vars2);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, STRONG_LOH);
    }

    #[test]
    fn test_strong_loh_complex_adj_before_placeholder() {
        // ports SomaticPostProcessModule.printVariationsFromSecondSample: adjComplex()
        // runs on v2var BEFORE SomaticOutputVariant is constructed from it, so the
        // printed POS/REF/ALT (taken from variant1) are the adjusted ones.
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        bam1.insert(150, Vars::new());

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        let mut v2 = make_good_variant(
            "-13^GCGGGGAGGGGAGAGG",
            "GAGGGCTAGGGGC",
            "GCGGGGAGGGGAGAGG",
            0.25,
        );
        v2.end_pos = 162;
        vars2.variants.push(v2);
        bam2.insert(150, vars2);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, STRONG_LOH);
        let out = &results[0].variant1;
        assert_eq!(out.start_pos, 151);
        assert_eq!(out.ref_allele, "AGGGCTAGGGGC");
        assert_eq!(out.var_allele, "CGGGGAGGGGAGAGG");
        let v2out = results[0].variant2.as_ref().unwrap();
        assert_eq!(v2out.start_pos, 151);
    }

    #[test]
    fn test_likely_loh() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        // BAM1 has a good variant at high frequency, BAM2 has same variant
        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        let v1 = make_good_variant("G", "A", "G", 0.96);
        vars1.variants.push(v1.clone());
        vars1.var_desc_to_variant.insert("G".to_string(), v1);
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        let v2 = make_good_variant("G", "A", "G", 0.50);
        vars2.variants.push(v2.clone());
        vars2.var_desc_to_variant.insert("G".to_string(), v2);
        bam2.insert(150, vars2);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, LIKELY_LOH);
    }

    #[test]
    fn test_position_outside_region_skipped() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars = Vars::new();
        let mut v = make_good_variant("G", "A", "G", 0.10);
        v.start_pos = 50;
        v.end_pos = 50;
        vars.variants.push(v);
        bam1.insert(50, vars); // position 50 < region.start (100)

        let bam2 = crate::javamap::JavaHashMap::new();

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_ref_eq_var_skipped() {
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars = Vars::new();
        let mut v = make_good_variant("A", "A", "A", 0.10);
        v.ref_allele = "A".to_string();
        v.var_allele = "A".to_string();
        vars.variants.push(v);
        bam1.insert(150, vars);

        let bam2 = crate::javamap::JavaHashMap::new();

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 0);
    }

    // ---- Category 11: combineAnalysis (SomaticPostProcessModule.java:385-476) ----

    fn adapt_pipeline<'a, F>(
        pipe: &'a F,
    ) -> impl Fn(
        &Region,
        &mut rustc_hash::FxHashSet<String>,
        &mut i32,
        i32,
    ) -> Result<Option<Vars>, VardictError>
           + 'a
    where
        F: Fn(
                &Region,
                &mut rustc_hash::FxHashSet<String>,
                &mut i32,
            ) -> Result<JavaHashMap<i32, Vars>, VardictError>
            + ?Sized,
    {
        move |comb_region, splice_set, max_rl, target_pos| {
            let mut vars_map = pipe(comb_region, splice_set, max_rl)?;
            Ok(vars_map.remove(&target_pos))
        }
    }

    fn make_insertion_variant() -> Variant {
        // Non-SNV, description longer than 10 chars and low coverage so the
        // combineAnalysis gate opens; Insertion type keeps adj_complex inert.
        let mut v = make_good_variant("+GAACCACTGAGACAGGA", "A", "AGAACCACTGAGACAGGA", 0.03);
        v.pos_coverage = 3;
        v.total_pos_cov = 10;
        v
    }

    fn combined_pipeline(
        pos: i32,
        desc: Option<&str>,
        pc: i32,
        total: i32,
    ) -> impl Fn(
        &Region,
        &mut rustc_hash::FxHashSet<String>,
        &mut i32,
    ) -> Result<crate::javamap::JavaHashMap<i32, Vars>, crate::types::VardictError>
           + '_ {
        move |_region: &Region, _splice: &mut rustc_hash::FxHashSet<String>, _max_rl: &mut i32| {
            let mut vars = Vars::new();
            if let Some(d) = desc {
                let mut v = Variant::new();
                v.description = d.to_string();
                v.pos_coverage = pc;
                v.total_pos_cov = total;
                vars.var_desc_to_variant.insert(d.to_string(), v.clone());
                vars.variants.push(v);
            }
            let mut m = crate::javamap::JavaHashMap::new();
            if pos != 0 {
                m.insert(pos, vars);
            }
            Ok(m)
        }
    }

    #[test]
    fn test_minus_num_num_find() {
        // ports Patterns.MINUS_NUM_NUM = -\d\d (Patterns.java:18)
        assert!(minus_num_num_find("-12^AAAAAATA"));
        assert!(minus_num_num_find("A/-12iAAAAAATA"));
        assert!(!minus_num_num_find("-2"));
        assert!(!minus_num_num_find("+GAACCACTGA"));
        assert!(!minus_num_num_find(""));
    }

    #[test]
    fn test_combine_not_found_is_skip() {
        // Java :475: getVarMaybe miss → FALSE. This is the HCC2218 +30 case.
        let config = test_config();
        let v = make_insertion_variant();
        let pipeline = combined_pipeline(150, None, 0, 0);
        let adapted = adapt_pipeline(&pipeline);
        let empty_splice = rustc_hash::FxHashSet::default();
        let mut splice = std::borrow::Cow::Borrowed(&empty_splice);
        let mut max_rl = 150;
        assert!(matches!(
            combine_analysis(&config, &adapted, &v, "chr1", 150, &mut splice, &mut max_rl),
            CombineVerdict::Skip
        ));
    }

    #[test]
    fn test_combine_lookup_uses_iteration_key_not_start_pos() {
        // Java :408 looks up by the ITERATION position key; the TumorA/Colo829
        // deletions sit at a key that differs from variant start_pos after
        // adjustment. Lookup at the key must hit even when start_pos would miss.
        let config = test_config();
        let mut v = make_insertion_variant();
        v.start_pos = 148; // differs from iteration key 150
        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 3, 10);
        let adapted = adapt_pipeline(&pipeline);
        let empty_splice = rustc_hash::FxHashSet::default();
        let mut splice = std::borrow::Cow::Borrowed(&empty_splice);
        let mut max_rl = 150;
        assert!(matches!(
            combine_analysis(&config, &adapted, &v, "chr1", 150, &mut splice, &mut max_rl),
            CombineVerdict::Keep
        ));
    }

    #[test]
    fn test_combine_keep_when_delta_below_minr() {
        // Java :471-473: found, delta < minr, and not < base-2 → ""
        let config = test_config();
        let v = make_insertion_variant(); // pc=3
        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 4, 20); // delta=1 < minr=2
        let adapted = adapt_pipeline(&pipeline);
        let empty_splice = rustc_hash::FxHashSet::default();
        let mut splice = std::borrow::Cow::Borrowed(&empty_splice);
        let mut max_rl = 150;
        assert!(matches!(
            combine_analysis(&config, &adapted, &v, "chr1", 150, &mut splice, &mut max_rl),
            CombineVerdict::Keep
        ));
    }

    #[test]
    fn test_combine_less_evidence_is_skip() {
        // Java :466-470: comb_pc < base_pc - 2 → FALSE ("Combine produce less")
        let config = test_config();
        let v = make_insertion_variant(); // pc=3
        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 0, 5); // 0 < 1
        let adapted = adapt_pipeline(&pipeline);
        let empty_splice = rustc_hash::FxHashSet::default();
        let mut splice = std::borrow::Cow::Borrowed(&empty_splice);
        let mut max_rl = 150;
        assert!(matches!(
            combine_analysis(&config, &adapted, &v, "chr1", 150, &mut splice, &mut max_rl),
            CombineVerdict::Skip
        ));
    }

    #[test]
    fn test_combine_sv_length_keeps() {
        // Java :394-396: span > SVMINLEN → "" without running the pipeline
        let config = test_config();
        let mut v = make_insertion_variant();
        v.end_pos = v.start_pos + config.sv_min_len as i32 + 1;
        let pipeline = combined_pipeline(150, None, 0, 0);
        let adapted = adapt_pipeline(&pipeline);
        let empty_splice = rustc_hash::FxHashSet::default();
        let mut splice = std::borrow::Cow::Borrowed(&empty_splice);
        let mut max_rl = 150;
        assert!(matches!(
            combine_analysis(&config, &adapted, &v, "chr1", 150, &mut splice, &mut max_rl),
            CombineVerdict::Keep
        ));
        assert!(matches!(splice, std::borrow::Cow::Borrowed(_)));

        // Large sv_min_len (u32::MAX) must not cast negative at somatic.rs:579, and
        // end_pos near i32::MAX / <= 0 must clamp cleanly to [1..=i32::MAX] at somatic.rs:587.
        let mut huge_sv_cfg = test_config();
        huge_sv_cfg.sv_min_len = u32::MAX;
        let seen_regions: std::cell::RefCell<Vec<(u32, u32)>> = std::cell::RefCell::new(Vec::new());
        let spy_pipeline =
            |comb_reg: &Region,
             _splice_set: &mut rustc_hash::FxHashSet<String>,
             _max_rl: &mut i32|
             -> Result<crate::javamap::JavaHashMap<i32, Vars>, VardictError> {
                seen_regions
                    .borrow_mut()
                    .push((comb_reg.start.0, comb_reg.end.0));
                Ok(crate::javamap::JavaHashMap::new())
            };
        let mut v_near_max = make_insertion_variant();
        v_near_max.start_pos = i32::MAX - 10;
        v_near_max.end_pos = i32::MAX - 10;
        let mut rl_150 = 150;
        assert!(matches!(
            combine_analysis(
                &huge_sv_cfg,
                &adapt_pipeline(&spy_pipeline),
                &v_near_max,
                "chr1",
                i32::MAX - 10,
                &mut splice,
                &mut rl_150
            ),
            CombineVerdict::Skip
        ));

        let mut v_nonpos = make_insertion_variant();
        v_nonpos.start_pos = -10;
        v_nonpos.end_pos = -5;
        let mut rl_zero = 0;
        assert!(matches!(
            combine_analysis(
                &huge_sv_cfg,
                &adapt_pipeline(&spy_pipeline),
                &v_nonpos,
                "chr1",
                1,
                &mut splice,
                &mut rl_zero
            ),
            CombineVerdict::Skip
        ));

        assert_eq!(
            seen_regions.into_inner(),
            vec![((i32::MAX - 160) as u32, i32::MAX as u32), (1, 1),]
        );
    }

    #[test]
    fn test_combine_germline_subtracts_other_sample() {
        // Java :413-465: delta >= minr → variant2 filled by subtraction → GERMLINE
        let config = test_config();
        let mut v = make_insertion_variant();
        v.genotype = Some("0/1".to_string());
        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 6, 20);
        let adapted = adapt_pipeline(&pipeline);
        let empty_splice = rustc_hash::FxHashSet::default();
        let mut splice = std::borrow::Cow::Borrowed(&empty_splice);
        let mut max_rl = 150;
        match combine_analysis(&config, &adapted, &v, "chr1", 150, &mut splice, &mut max_rl) {
            CombineVerdict::Germline(var2) => {
                assert_eq!(var2.pos_coverage, 3); // 6 - 3
                assert_eq!(var2.total_pos_cov, 10); // 20 - 10
                assert!((var2.frequency - 0.3).abs() < 1e-9);
                assert!(var2.pstd && var2.qstd); // :453-454
            }
            _ => panic!("expected Germline"),
        }
    }

    #[test]
    fn test_first_sample_false_skips_record() {
        // printVariationsFromFirstSample :202-204 — FALSE drops the tumor-only
        // record entirely (the reverted port printed it as spurious StrongSomatic)
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        vars1.variants.push(make_insertion_variant());
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        bam2.insert(150, Vars::new());

        let pipeline = combined_pipeline(150, None, 0, 0); // combined run finds nothing
        let results = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &splice,
            Some(&pipeline),
            150,
        );
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_first_sample_keep_strong_somatic() {
        // Java "" verdict keeps the StrongSomatic record (:211-213)
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        vars1.variants.push(make_insertion_variant());
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        bam2.insert(150, Vars::new());

        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 3, 10); // delta 0
        let results = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &splice,
            Some(&pipeline),
            150,
        );
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, STRONG_SOMATIC);
        assert!(results[0].variant2.is_none()); // varForPrint null: v2 empty, no ref
    }

    #[test]
    fn test_first_sample_germline_carries_adjusted_placeholder() {
        // Java :215 — Germline verdict prints with the adjusted v2nt as variant2
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        vars1.variants.push(make_insertion_variant());
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        bam2.insert(150, Vars::new());

        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 8, 40); // delta 5 >= minr
        let results = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &splice,
            Some(&pipeline),
            150,
        );
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, GERMLINE);
        let var2 = results[0].variant2.as_ref().unwrap();
        assert_eq!(var2.pos_coverage, 5); // 8 - 3
        assert!((var2.frequency - 5.0 / 30.0).abs() < 1e-9);
        assert!(!results[0].end_is_variant2); // end fields stay on vref (:215)
    }

    #[test]
    fn test_second_sample_false_skips_record() {
        // printVariationsFromSecondSample :313-315 — FALSE skips where the old
        // conflated port fell back to StrongLOH
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        bam1.insert(150, Vars::new());

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        let mut v = make_insertion_variant();
        v.ref_allele = "C".to_string();
        v.var_allele = "CGAACCACTGAGACAGGA".to_string();
        vars2.variants.push(v);
        bam2.insert(150, vars2);

        let pipeline = combined_pipeline(150, None, 0, 0);
        let results = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &splice,
            Some(&pipeline),
            150,
        );
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_second_sample_keep_strong_loh_fallback() {
        // Java "" keeps StrongLOH with sample1's reference as varForPrint
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        let mut vref1 = Variant::new();
        vref1.description = "A".to_string();
        vref1.total_pos_cov = 50;
        vars1.ref_variant = Some(vref1);
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        let mut v = make_insertion_variant();
        v.ref_allele = "C".to_string();
        v.var_allele = "CGAACCACTGAGACAGGA".to_string();
        vars2.variants.push(v);
        bam2.insert(150, vars2);

        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 3, 10);
        let results = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &splice,
            Some(&pipeline),
            150,
        );
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, STRONG_LOH);
        assert_eq!(results[0].variant1.total_pos_cov, 50); // from v1.referenceVariant
    }

    #[test]
    fn test_is_noise_mutation_on_first_v2_variant_propagates_to_second_tumor_allele() {
        // Bug E (`SomaticPostProcessModule.java:172-186`, `Variant.java:205`):
        // When the first tumor variant at a locus matches `v2.variants.get(0)` and
        // `determinateType` calls `v2nt.isNoise()`, `isNoise()` decrements
        // `v2nt.totalPosCoverage -= v2nt.positionCoverage` in place. Because
        // `v2.variants.get(0)` and `v2.varDescriptionStringToVariants.get(nt)` point
        // to the same `Variant` object in Java, a second tumor variant at the same
        // locus that is absent in BAM2 (`v2nt == null`) reads the decremented
        // `v2.variants.get(0).totalPosCoverage`.
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1 = Vars::new();
        let t1 = make_good_variant("G", "A", "G", 0.10);
        let t2 = make_good_variant("T", "A", "T", 0.05);
        vars1.variants.push(t1);
        vars1.variants.push(t2);
        bam1.insert(150, vars1);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars2 = Vars::new();
        // Normal has a 1-read low-quality noise variant "G" at index 0 (total_pos_cov = 212)
        let mut n1 = make_good_variant("G", "A", "G", 1.0 / 212.0);
        n1.pos_coverage = 1;
        n1.total_pos_cov = 212;
        n1.mean_quality = 10.0; // triggers is_noise()
        n1.hicnt = 0;
        vars2.variants.push(n1.clone());
        vars2.var_desc_to_variant.insert("G".to_string(), n1);
        bam2.insert(150, vars2);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].status, STRONG_SOMATIC);
        assert_eq!(results[0].variant2.as_ref().unwrap().total_pos_cov, 211);
        assert_eq!(results[1].status, STRONG_SOMATIC);
        assert_eq!(
            results[1].variant2.as_ref().unwrap().total_pos_cov,
            211,
            "second tumor variant must observe v2.variants[0].total_pos_cov mutated by is_noise()"
        );
    }

    #[test]
    fn combine_analysis_propagates_updated_max_read_length_and_splice_to_subsequent_windows() {
        // SomaticPostProcessModule.java:198-200, 303-311, 398-406:
        // When `combineAnalysis` runs on an indel at position 150 and discovers a longer
        // `maxReadLength` (e.g. 100 -> 250) and new splice junctions, it updates
        // `this.maxReadLength` and the shared `splice` set so the next `combineAnalysis`
        // call at position 180 queries `180 - 250 .. 180 + 250` (`1..430`) and sees the
        // updated `splice` set.
        let config = test_config();
        let region = test_region();
        let splice = rustc_hash::FxHashSet::default();

        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut vars1_150 = Vars::new();
        let mut v150 = make_insertion_variant();
        v150.start_pos = 150;
        v150.end_pos = 150;
        vars1_150.variants.push(v150);
        bam1.insert(150, vars1_150);

        let mut vars1_180 = Vars::new();
        let mut v180 = make_insertion_variant();
        v180.start_pos = 180;
        v180.end_pos = 180;
        vars1_180.variants.push(v180);
        bam1.insert(180, vars1_180);

        let mut bam2 = crate::javamap::JavaHashMap::new();
        bam2.insert(150, Vars::new());
        bam2.insert(180, Vars::new());

        let observed_windows: std::cell::RefCell<Vec<(u32, u32, bool)>> =
            std::cell::RefCell::new(Vec::new());
        let pipeline = |comb_reg: &Region,
                        splice_set: &mut rustc_hash::FxHashSet<String>,
                        max_rl: &mut i32|
         -> Result<crate::javamap::JavaHashMap<i32, Vars>, VardictError> {
            observed_windows.borrow_mut().push((
                comb_reg.start.0,
                comb_reg.end.0,
                splice_set.contains("155-175"),
            ));
            *max_rl = (*max_rl).max(250);
            splice_set.insert("155-175".to_string());
            Ok(crate::javamap::JavaHashMap::new())
        };

        let _ = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &splice,
            Some(&pipeline),
            100,
        );
        let windows = observed_windows.into_inner();
        assert_eq!(windows.len(), 2);
        // First combine_analysis at pos 150 with initial max_read_length = 100 -> [50, 250], splice empty
        assert_eq!(windows[0], (50, 250, false));
        // Second combine_analysis at pos 180 must use updated max_read_length = 250 -> [1, 430] and updated splice
        assert_eq!(windows[1], (1, 430, true));
    }

    #[test]
    fn test_matched_and_sv_variants_keep_splice_cow_borrowed_without_pipeline() {
        let config = test_config();
        let mut region = test_region();
        region.end = Pos1(500);

        let mut initial_splice = rustc_hash::FxHashSet::default();
        initial_splice.insert("105-115".to_string());

        let pipeline_calls = std::cell::Cell::new(0usize);
        let pipeline = |_comb_reg: &Region,
                        _splice_set: &mut rustc_hash::FxHashSet<String>,
                        _max_rl: &mut i32|
         -> Result<crate::javamap::JavaHashMap<i32, Vars>, VardictError> {
            pipeline_calls.set(pipeline_calls.get() + 1);
            Ok(crate::javamap::JavaHashMap::new())
        };

        // 1. Matched Complex variant (triggers adj_complex) + second matched SNV at pos 150
        //    so both `nt` lookup and `v2.to_mut()` + `v2_item.description == nt` execute.
        let mut vars1_150 = Vars::new();
        let mut c1 = make_good_variant(
            "-13^GCGGGGAGGGGAGAGG",
            "GAGGGCTAGGGGC",
            "GCGGGGAGGGGAGAGG",
            0.35,
        );
        c1.end_pos = 162;
        let s1 = make_good_variant("G", "A", "G", 0.30);
        vars1_150.variants.push(c1.clone());
        vars1_150.variants.push(s1.clone());
        vars1_150
            .var_desc_to_variant
            .insert(c1.description.clone(), c1.clone());
        vars1_150
            .var_desc_to_variant
            .insert(s1.description.clone(), s1.clone());

        let mut vars2_150 = Vars::new();
        let mut c2 = c1;
        c2.frequency = 0.30;
        let mut s2 = s1;
        s2.frequency = 0.25;
        vars2_150.variants.push(c2.clone());
        vars2_150.variants.push(s2.clone());
        vars2_150
            .var_desc_to_variant
            .insert(c2.description.clone(), c2);
        vars2_150
            .var_desc_to_variant
            .insert(s2.description.clone(), s2);

        let mut splice_cow: std::borrow::Cow<'_, rustc_hash::FxHashSet<String>> =
            std::borrow::Cow::Borrowed(&initial_splice);
        let mut results = Vec::new();
        let mut max_rl = 150;

        let adapted_pipeline = adapt_pipeline(&pipeline);
        calling_for_both_samples(
            &config,
            150,
            std::borrow::Cow::Borrowed(&vars1_150),
            std::borrow::Cow::Borrowed(&vars2_150),
            "chr1",
            &mut splice_cow,
            &mut results,
            Some(&adapted_pipeline),
            &mut max_rl,
        );
        assert_eq!(pipeline_calls.get(), 0);
        assert!(
            matches!(splice_cow, std::borrow::Cow::Borrowed(_)),
            "matched tumor/normal variants must not clone splice"
        );
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].status, GERMLINE);
        assert_eq!(results[0].variant1.start_pos, 151); // adj_complex trimmed leading 'G'
        assert_eq!(results[1].status, GERMLINE);

        // 2. Unmatched tumor-only SV-span variant entering combine_analysis and returning
        //    early at the sv_min_len guard must also leave splice_cow Borrowed.
        let mut sv_vars1 = Vars::new();
        let mut sv_var = make_insertion_variant();
        sv_var.start_pos = 160;
        sv_var.end_pos = 160 + config.sv_min_len as i32 + 1;
        sv_vars1.variants.push(sv_var);
        let empty_vars2 = Vars::new();

        calling_for_both_samples(
            &config,
            160,
            std::borrow::Cow::Borrowed(&sv_vars1),
            std::borrow::Cow::Borrowed(&empty_vars2),
            "chr1",
            &mut splice_cow,
            &mut results,
            Some(&adapted_pipeline),
            &mut max_rl,
        );
        assert_eq!(pipeline_calls.get(), 0);
        assert!(
            matches!(splice_cow, std::borrow::Cow::Borrowed(_)),
            "SV early return in combine_analysis must not clone splice"
        );
        assert_eq!(results.len(), 3);
        assert_eq!(results[2].status, STRONG_SOMATIC);

        // 3. End-to-end somatic_post_process over both positions never invokes pipeline.
        let mut bam1 = crate::javamap::JavaHashMap::new();
        bam1.insert(150, vars1_150);
        bam1.insert(160, sv_vars1);
        let mut bam2 = crate::javamap::JavaHashMap::new();
        bam2.insert(150, vars2_150);
        bam2.insert(160, empty_vars2);

        let e2e = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &initial_splice,
            Some(&pipeline),
            150,
        );
        assert_eq!(pipeline_calls.get(), 0);
        assert_eq!(e2e.len(), 3);
    }

    #[test]
    fn test_unmatched_combine_analysis_lazily_clones_splice_and_propagates() {
        let config = test_config();
        let region = test_region();

        let mut initial_splice = rustc_hash::FxHashSet::default();
        initial_splice.insert("105-115".to_string());

        // 1. Direct check on calling_for_both_samples (second-sample unmatched path):
        //    splice starts Borrowed, transitions to Owned via Cow::to_mut() inside
        //    combine_analysis, and leaves initial_splice untouched.
        let mut splice_cow: std::borrow::Cow<'_, rustc_hash::FxHashSet<String>> =
            std::borrow::Cow::Borrowed(&initial_splice);
        let empty_v1 = Vars::new();
        let mut unmatched_v2 = Vars::new();
        unmatched_v2.variants.push(make_insertion_variant());
        let mut direct_results = Vec::new();
        let mut max_rl = 150;
        let single_pipeline =
            |_comb_reg: &Region,
             splice_set: &mut rustc_hash::FxHashSet<String>,
             _max_rl: &mut i32|
             -> Result<crate::javamap::JavaHashMap<i32, Vars>, VardictError> {
                splice_set.insert("155-175".to_string());
                Ok(crate::javamap::JavaHashMap::new())
            };
        let adapted_single = adapt_pipeline(&single_pipeline);
        calling_for_both_samples(
            &config,
            150,
            std::borrow::Cow::Borrowed(&empty_v1),
            std::borrow::Cow::Borrowed(&unmatched_v2),
            "chr1",
            &mut splice_cow,
            &mut direct_results,
            Some(&adapted_single),
            &mut max_rl,
        );
        assert!(
            matches!(splice_cow, std::borrow::Cow::Owned(_)),
            "combine_analysis must lazily materialize Owned splice via Cow::to_mut()"
        );
        assert!(splice_cow.contains("105-115"));
        assert!(splice_cow.contains("155-175"));
        assert!(!initial_splice.contains("155-175"));

        // 2. End-to-end somatic_post_process across 3 positions:
        //    - pos 120: matched SNV (does not invoke pipeline)
        //    - pos 150: unmatched tumor-only insertion (clones splice, inserts "155-175")
        //    - pos 180: unmatched normal-only insertion (observes "155-175" in splice)
        let mut bam1 = crate::javamap::JavaHashMap::new();
        let mut bam2 = crate::javamap::JavaHashMap::new();

        let mut m1 = Vars::new();
        let mut mv1 = make_good_variant("G", "A", "G", 0.40);
        mv1.start_pos = 120;
        mv1.end_pos = 120;
        m1.variants.push(mv1.clone());
        m1.var_desc_to_variant.insert("G".to_string(), mv1.clone());
        bam1.insert(120, m1);

        let mut m2 = Vars::new();
        m2.variants.push(mv1.clone());
        m2.var_desc_to_variant.insert("G".to_string(), mv1);
        bam2.insert(120, m2);

        let mut t150 = Vars::new();
        let mut v150 = make_insertion_variant();
        v150.start_pos = 150;
        v150.end_pos = 150;
        t150.variants.push(v150);
        bam1.insert(150, t150);
        bam2.insert(150, Vars::new());

        bam1.insert(180, Vars::new());
        let mut n180 = Vars::new();
        let mut v180 = make_insertion_variant();
        v180.start_pos = 180;
        v180.end_pos = 180;
        n180.variants.push(v180);
        bam2.insert(180, n180);

        let observed: std::cell::RefCell<Vec<(u32, bool, bool)>> =
            std::cell::RefCell::new(Vec::new());
        let e2e_pipeline =
            |comb_reg: &Region,
             splice_set: &mut rustc_hash::FxHashSet<String>,
             _max_rl: &mut i32|
             -> Result<crate::javamap::JavaHashMap<i32, Vars>, VardictError> {
                observed.borrow_mut().push((
                    comb_reg.start.0,
                    splice_set.contains("105-115"),
                    splice_set.contains("155-175"),
                ));
                splice_set.insert("155-175".to_string());
                Ok(crate::javamap::JavaHashMap::new())
            };

        let out = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &initial_splice,
            Some(&e2e_pipeline),
            100,
        );
        assert_eq!(out.len(), 1); // pos 120 Germline; pos 150 & 180 skipped by empty combined map
        assert_eq!(out[0].status, GERMLINE);

        let obs = observed.into_inner();
        assert_eq!(obs.len(), 2);
        assert_eq!(obs[0], (50, true, false)); // pos 150 - 100 = 50; sees initial_splice only
        assert_eq!(obs[1], (80, true, true)); // pos 180 - 100 = 80; sees propagated "155-175"
        assert!(!initial_splice.contains("155-175"));
    }

    #[test]
    fn extreme_somatic_i64_gates_and_bias_buffer_reuse() {
        // Verify pre-allocated capacity-4 bias buffer in Variant::new() is reused without realloc.
        let mut v_buf = Variant::new();
        let initial_ptr = v_buf.bias.as_ptr();
        assert_eq!(v_buf.bias.capacity(), 4);
        v_buf.bias.clear();
        v_buf.bias.push_str("0;0");
        assert_eq!(v_buf.bias, "0;0");
        assert_eq!(v_buf.bias.as_ptr(), initial_ptr);
        v_buf.bias.clear();
        let _ = std::fmt::Write::write_fmt(&mut v_buf.bias, format_args!("{};{}", 2, 1));
        assert_eq!(v_buf.bias, "2;1");
        assert_eq!(v_buf.bias.as_ptr(), initial_ptr);

        // Verify combine_analysis handles min_reads = u32::MAX, min_reads_bias = u32::MAX,
        // and base_pc = i32::MIN without signed i32 subtraction underflow in `comb_pc < base_pc - 2`.
        let mut config = test_config();
        config.min_reads = u32::MAX;
        config.min_reads_bias = u32::MAX;
        let initial_splice = rustc_hash::FxHashSet::default();
        let mut splice_cow = std::borrow::Cow::Borrowed(&initial_splice);
        let mut max_rl = 100;

        let mut base_var = make_good_variant("G", "A", "G", 0.30);
        base_var.start_pos = 150;
        base_var.end_pos = 150;
        base_var.pos_coverage = i32::MIN;

        let pipeline = |_reg: &Region,
                        _sp: &mut rustc_hash::FxHashSet<String>,
                        _mrl: &mut i32|
         -> Result<crate::javamap::JavaHashMap<i32, Vars>, VardictError> {
            let mut map = crate::javamap::JavaHashMap::new();
            let mut pv = Vars::new();
            let mut comb = make_good_variant("G", "A", "G", 0.30);
            comb.pos_coverage = i32::MIN;
            pv.variants.push(comb.clone());
            pv.var_desc_to_variant.insert("G".to_string(), comb);
            map.insert(150, pv);
            Ok(map)
        };

        let adapted = adapt_pipeline(&pipeline);
        let verdict = combine_analysis(
            &config,
            &adapted,
            &base_var,
            "chr1",
            150,
            &mut splice_cow,
            &mut max_rl,
        );
        // i64::from(i32::MIN) < i64::from(i32::MIN) - 2 is FALSE, so Keep (not Skip).
        assert!(matches!(verdict, CombineVerdict::Keep));
    }

    #[test]
    fn somatic_post_process_at_passes_outer_pileup_position_not_shifted_start_pos_and_skips_on_err()
    {
        let config = test_config();
        let mut region = test_region();
        region.start = Pos1(100);
        region.end = Pos1(300);
        let initial_splice = rustc_hash::FxHashSet::default();

        let mut bam1 = JavaHashMap::new();
        let mut bam2 = JavaHashMap::new();

        // Pos 150 (sample 1 tumor-only deletion "-12" whose start_pos was shifted to 149):
        // combine_analysis MUST pass outer key 150 (NOT 149) while region uses start_pos 149!
        let mut t150 = Vars::new();
        let mut v150 = make_good_variant("-12", "ACGTACGTACGTN", "A", 0.30);
        v150.var_type = "Deletion".to_string();
        v150.start_pos = 149;
        v150.end_pos = 161;
        v150.pos_coverage = 3;
        v150.total_pos_cov = 10;
        t150.variants.push(v150);
        bam1.insert(150, t150);
        bam2.insert(150, Vars::new());

        // Pos 200 (sample 2 normal-only +AT&GC complex insertion whose start_pos shifted to 201,
        // and -3 shift3 shifted further to 204): combine_analysis MUST pass outer key 200!
        let mut n200 = Vars::new();
        let mut v200 = make_good_variant("+ACGTACGTAC&GG", "CG", "ACGTACGTACGG", 0.30);
        v200.var_type = "Complex".to_string();
        v200.start_pos = 204;
        v200.end_pos = 206;
        v200.pos_coverage = 3;
        v200.total_pos_cov = 10;
        n200.variants.push(v200);
        bam1.insert(200, Vars::new());
        bam2.insert(200, n200);

        // Pos 250 (sample 1 tumor-only insertion where pipeline returns Err -> Skip, infallible!).
        let mut t250 = Vars::new();
        let mut v250 = make_insertion_variant();
        v250.start_pos = 253; // shifted by shift3
        v250.end_pos = 253;
        t250.variants.push(v250);
        bam1.insert(250, t250);
        bam2.insert(250, Vars::new());

        let seen: std::cell::RefCell<Vec<(u32, u32, i32)>> = std::cell::RefCell::new(Vec::new());
        let pos_pipe = |comb_reg: &Region,
                        _sp: &mut rustc_hash::FxHashSet<String>,
                        _mrl: &mut i32,
                        target_pos: i32|
         -> Result<Option<Vars>, VardictError> {
            seen.borrow_mut()
                .push((comb_reg.start.0, comb_reg.end.0, target_pos));
            match target_pos {
                150 => {
                    let mut pv = Vars::new();
                    let mut comb = Variant::new();
                    comb.description = "-12".to_string();
                    comb.pos_coverage = 6; // delta = 3 >= minr (2) -> Germline
                    comb.total_pos_cov = 20;
                    comb.genotype = Some("A/-12".to_string());
                    pv.variants.push(comb);
                    Ok(Some(pv))
                }
                200 => {
                    let mut pv = Vars::new();
                    let mut comb = Variant::new();
                    comb.description = "+ACGTACGTAC&GG".to_string();
                    comb.pos_coverage = 3; // delta = 0 < minr -> Keep (StrongLOH)
                    comb.total_pos_cov = 20;
                    pv.variants.push(comb);
                    Ok(Some(pv))
                }
                250 => {
                    if comb_reg.start.0 == 203 {
                        Err(VardictError::Unsupported(
                            "synthetic combine failure".into(),
                        ))
                    } else {
                        Ok(None)
                    }
                }
                other => panic!("unexpected target_pos {other}; must be outer pileup key"),
            }
        };

        let at_out = somatic_post_process_at(
            &config,
            &region,
            &bam1,
            &bam2,
            &initial_splice,
            Some(&pos_pipe),
            50,
        );
        assert_eq!(
            seen.borrow().as_slice(),
            &[
                (149 - 50, 161 + 50, 150),
                (204 - 50, 206 + 50, 200),
                (253 - 50, 253 + 50, 250),
            ]
        );
        assert_eq!(at_out.len(), 2);
        assert_eq!(at_out[0].status, GERMLINE);
        assert_eq!(at_out[0].variant1.start_pos, 149);
        assert_eq!(at_out[1].status, STRONG_LOH);

        // Verify somatic_post_process wrapping PipelineFn produces identical output.
        let map_pipe = |comb_reg: &Region,
                        sp: &mut rustc_hash::FxHashSet<String>,
                        mrl: &mut i32|
         -> Result<JavaHashMap<i32, Vars>, VardictError> {
            let mut m = JavaHashMap::new();
            for p in [150, 200, 250] {
                if (comb_reg.start.0..=comb_reg.end.0).contains(&(p as u32)) {
                    if let Some(v) = pos_pipe(comb_reg, sp, mrl, p)? {
                        m.insert(p, v);
                    }
                }
            }
            Ok(m)
        };
        let wrapped_out = somatic_post_process(
            &config,
            &region,
            &bam1,
            &bam2,
            &initial_splice,
            Some(&map_pipe),
            50,
        );
        assert_eq!(wrapped_out.len(), at_out.len());
        assert_eq!(wrapped_out[0].status, at_out[0].status);
        assert_eq!(wrapped_out[1].status, at_out[1].status);
    }

    #[test]
    fn somatic_post_process_rejects_negative_i32_keys_under_large_region_end() {
        let config = Config::default();
        // With `end = u32::MAX`, `(-10_i32 as u32)` (`4294967286`) previously slipped past
        // `(p as u32) >= region.start.0 && (p as u32) <= region.end.0`.
        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(u32::MAX),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut bam1 = JavaHashMap::new();
        let mut vars = Vars::new();
        let mut v = Variant::new();
        v.description = "G".to_string();
        v.ref_allele = "A".to_string();
        v.var_allele = "G".to_string();
        v.frequency = 0.5;
        v.pos_coverage = 10;
        v.total_pos_cov = 20;
        v.hicnt = 10;
        v.hicov = 20;
        v.fwd = 5;
        v.rev = 5;
        v.ref_fwd_cov = 5;
        v.ref_rev_cov = 5;
        v.mean_position = 15.0;
        v.mean_quality = 35.0;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 20.0;
        v.pstd = true;
        v.qstd = true;
        v.start_pos = -10;
        v.end_pos = -10;
        vars.variants.push(v);
        bam1.insert(-10, vars);
        let bam2 = JavaHashMap::new();
        let empty_splice = rustc_hash::FxHashSet::default();

        let out = somatic_post_process_at(&config, &region, &bam1, &bam2, &empty_splice, None, 100);
        assert!(
            out.is_empty(),
            "negative i32 key -10 must be excluded by somatic_post_process_at, got {} records",
            out.len()
        );

        // Exhaustive field-for-field parity between `somatic_post_process_owned` and
        // `somatic_post_process_at` across SampleSpecific, Deletion, StrongSomatic, Germline,
        // LikelyLOH, StrongLOH (number_of_processed_variation == 0 when v1.variants[0] fails
        // is_good_var while v2 passes, plus multi-allelic v2), combineAnalysis, duplicate-
        // description v2 entries, and multi-allelic tumor + is_noise() normal mutation (#45).
        fn assert_variant_eq(a: &Variant, b: &Variant) {
            let Variant {
                description,
                pos_coverage,
                fwd,
                rev,
                bias,
                frequency,
                mean_position,
                pstd,
                mean_quality,
                qstd,
                mean_mapping_quality,
                hq_lq_ratio,
                hi_freq,
                extra_freq,
                shift3,
                msi,
                msint,
                nm,
                hicnt,
                hicov,
                left_seq,
                right_seq,
                start_pos,
                end_pos,
                ref_rev_cov,
                ref_fwd_cov,
                total_pos_cov,
                duprate,
                genotype,
                var_allele,
                ref_allele,
                var_type,
                debug,
                crispr,
                sv,
                null_record,
            } = a;
            assert_eq!(description, &b.description);
            assert_eq!(pos_coverage, &b.pos_coverage);
            assert_eq!(fwd, &b.fwd);
            assert_eq!(rev, &b.rev);
            assert_eq!(bias, &b.bias);
            assert_eq!(frequency.to_bits(), b.frequency.to_bits());
            assert_eq!(mean_position.to_bits(), b.mean_position.to_bits());
            assert_eq!(pstd, &b.pstd);
            assert_eq!(mean_quality.to_bits(), b.mean_quality.to_bits());
            assert_eq!(qstd, &b.qstd);
            assert_eq!(
                mean_mapping_quality.to_bits(),
                b.mean_mapping_quality.to_bits()
            );
            assert_eq!(hq_lq_ratio.to_bits(), b.hq_lq_ratio.to_bits());
            assert_eq!(hi_freq.to_bits(), b.hi_freq.to_bits());
            assert_eq!(extra_freq.to_bits(), b.extra_freq.to_bits());
            assert_eq!(shift3, &b.shift3);
            assert_eq!(msi.to_bits(), b.msi.to_bits());
            assert_eq!(msint, &b.msint);
            assert_eq!(nm.to_bits(), b.nm.to_bits());
            assert_eq!(hicnt, &b.hicnt);
            assert_eq!(hicov, &b.hicov);
            assert_eq!(left_seq, &b.left_seq);
            assert_eq!(right_seq, &b.right_seq);
            assert_eq!(start_pos, &b.start_pos);
            assert_eq!(end_pos, &b.end_pos);
            assert_eq!(ref_rev_cov, &b.ref_rev_cov);
            assert_eq!(ref_fwd_cov, &b.ref_fwd_cov);
            assert_eq!(total_pos_cov, &b.total_pos_cov);
            assert_eq!(duprate.to_bits(), b.duprate.to_bits());
            assert_eq!(genotype, &b.genotype);
            assert_eq!(var_allele, &b.var_allele);
            assert_eq!(ref_allele, &b.ref_allele);
            assert_eq!(var_type, &b.var_type);
            assert_eq!(debug, &b.debug);
            assert_eq!(crispr, &b.crispr);
            assert_eq!(sv, &b.sv);
            assert_eq!(null_record, &b.null_record);
        }

        fn assert_somatic_eq(a: &[SomaticVariant], b: &[SomaticVariant]) {
            assert_eq!(a.len(), b.len());
            for (got, want) in a.iter().zip(b.iter()) {
                let SomaticVariant {
                    variant1,
                    variant2,
                    status,
                    sv1,
                    sv2,
                    end_is_variant2,
                    tumor_amplicon,
                    normal_amplicon,
                } = got;
                assert_variant_eq(variant1, &want.variant1);
                assert_eq!(variant2.is_some(), want.variant2.is_some());
                if let (Some(v2a), Some(v2b)) = (variant2.as_ref(), want.variant2.as_ref()) {
                    assert_variant_eq(v2a, v2b);
                }
                assert_eq!(status, &want.status);
                assert_eq!(sv1, &want.sv1);
                assert_eq!(sv2, &want.sv2);
                assert_eq!(end_is_variant2, &want.end_is_variant2);
                assert_eq!(tumor_amplicon, &want.tumor_amplicon);
                assert_eq!(normal_amplicon, &want.normal_amplicon);
            }
        }

        let cfg = test_config();
        let reg = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(300),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let mut p_bam1 = JavaHashMap::new();
        let mut p_bam2 = JavaHashMap::new();
        p_bam1.insert(-10, Vars::new());

        // 105: SampleSpecific (multi-allelic + Complex + sv1)
        let mut s1_105 = Vars::new();
        s1_105.sv = "105-190-2".to_string();
        let mut v_105a = make_good_variant("G", "A", "G", 0.35);
        v_105a.start_pos = 105;
        v_105a.end_pos = 105;
        let mut v_105b = make_good_variant("#GG", "ACGT", "AGGT", 0.20);
        v_105b.start_pos = 105;
        v_105b.end_pos = 108;
        s1_105.variants.push(v_105a);
        s1_105.variants.push(v_105b);
        p_bam1.insert(105, s1_105);

        // 110: Deletion (only BAM2 covered, multi-allelic + Complex + sv2)
        let mut s2_110 = Vars::new();
        s2_110.sv = "110-195-3".to_string();
        let mut v_110a = make_good_variant("T", "C", "T", 0.40);
        v_110a.start_pos = 110;
        v_110a.end_pos = 110;
        let mut v_110b = make_good_variant("#TT", "CGTA", "CTTA", 0.25);
        v_110b.start_pos = 110;
        v_110b.end_pos = 113;
        s2_110.variants.push(v_110a);
        s2_110.variants.push(v_110b);
        p_bam2.insert(110, s2_110);

        // 115: Single-allele StrongSomatic (not in normal; normal has ref_variant)
        let mut s1_115 = Vars::new();
        s1_115.sv = "115-160-1".to_string();
        let mut v_115 = make_good_variant("G", "A", "G", 0.30);
        v_115.start_pos = 115;
        v_115.end_pos = 115;
        s1_115.variants.push(v_115);
        let mut s2_115 = Vars::new();
        s2_115.sv = "115-160-2".to_string();
        let mut ref_115 = make_good_variant("A", "A", "A", 1.0);
        ref_115.start_pos = 115;
        ref_115.end_pos = 115;
        s2_115.ref_variant = Some(ref_115);
        p_bam1.insert(115, s1_115);
        p_bam2.insert(115, s2_115);

        // 120: Single-allele Germline + single-allele LikelyLOH (122) + single-allele noise SNV StrongSomatic (125)
        let mut s1_120 = Vars::new();
        s1_120.sv = "sv1-120".to_string();
        s1_120.variants.push(make_good_variant("G", "A", "G", 0.45));
        let mut s2_120 = Vars::new();
        s2_120.sv = "sv2-120".to_string();
        s2_120.variants.push(make_good_variant("G", "A", "G", 0.42));
        p_bam1.insert(120, s1_120);
        p_bam2.insert(120, s2_120);

        let mut s1_122 = Vars::new();
        s1_122.variants.push(make_good_variant("T", "C", "T", 0.98));
        let mut s2_122 = Vars::new();
        s2_122.variants.push(make_good_variant("T", "C", "T", 0.50));
        p_bam1.insert(122, s1_122);
        p_bam2.insert(122, s2_122);

        let mut s1_125 = Vars::new();
        s1_125.variants.push(make_good_variant("G", "A", "G", 0.35));
        let mut s2_125 = Vars::new();
        let mut v2_noise_single = make_good_variant("G", "A", "G", 0.02);
        v2_noise_single.mean_quality = 10.0;
        v2_noise_single.qstd = false;
        v2_noise_single.pos_coverage = 2;
        v2_noise_single.total_pos_cov = 100;
        s2_125.variants.push(v2_noise_single);
        p_bam1.insert(125, s1_125);
        p_bam2.insert(125, s2_125);

        // 130: StrongLOH + LikelyLOH + Germline when v1.variants[0] FAILS is_good_var
        // (`number_of_processed_variation == 0`) while multi-allelic v2 passes!
        let mut s1_130 = Vars::new();
        s1_130.sv = "sv1-130".to_string();
        let mut v1_bad = make_good_variant("C", "A", "C", 0.002); // fails freq
        v1_bad.total_pos_cov = 88;
        v1_bad.genotype = Some("A/C".to_string());
        let mut v1_low_loh = make_good_variant("#GG", "ACGT", "AGGT", 0.02); // < lofreq -> LikelyLOH + Complex
        v1_low_loh.start_pos = 130;
        v1_low_loh.end_pos = 133;
        let v1_germ_loh = make_good_variant("T", "A", "T", 0.30); // >= lofreq -> Germline
        s1_130.variants.push(v1_bad);
        s1_130.variants.push(v1_low_loh);
        s1_130.variants.push(v1_germ_loh);
        let mut ref1_130 = Variant::new();
        ref1_130.description = "A".to_string();
        ref1_130.fwd = 40;
        ref1_130.rev = 42;
        s1_130.ref_variant = Some(ref1_130);

        let mut s2_130 = Vars::new();
        s2_130.sv = "sv2-130".to_string();
        let mut v2_loh_only = make_good_variant("#CC", "ACGT", "ACCT", 0.45); // absent in v1 -> StrongLOH
        v2_loh_only.start_pos = 130;
        v2_loh_only.end_pos = 133;
        let mut v2_likely_loh = make_good_variant("#GG", "ACGT", "AGGT", 0.50); // matches v1_low_loh
        v2_likely_loh.start_pos = 130;
        v2_likely_loh.end_pos = 133;
        let v2_germ_loh = make_good_variant("T", "A", "T", 0.48); // matches v1_germ_loh
        s2_130.variants.push(v2_loh_only);
        s2_130.variants.push(v2_likely_loh);
        s2_130.variants.push(v2_germ_loh);
        p_bam1.insert(130, s1_130);
        p_bam2.insert(130, s2_130);

        // 135: Duplicate-description v2 entries (variant_by_desc picks last / var_desc_to_variant)
        let mut s1_135 = Vars::new();
        s1_135.variants.push(make_good_variant("G", "A", "G", 0.40));
        let mut s2_135 = Vars::new();
        let mut v2_dup_first = make_good_variant("G", "A", "G", 0.02);
        v2_dup_first.pos_coverage = 1;
        let v2_dup_second = make_good_variant("G", "A", "G", 0.45);
        s2_135.variants.push(v2_dup_first);
        s2_135.variants.push(v2_dup_second);
        p_bam1.insert(135, s1_135);
        p_bam2.insert(135, s2_135);

        // 137: Single-element v2.variants PLUS distinct var_desc_to_variant entry
        let mut s1_137 = Vars::new();
        s1_137.variants.push(make_good_variant("G", "A", "G", 0.40));
        let mut s2_137 = Vars::new();
        s2_137.variants.push(make_good_variant("T", "A", "T", 0.40));
        s2_137
            .var_desc_to_variant
            .insert("G".to_string(), make_good_variant("G", "A", "G", 0.38));
        p_bam1.insert(137, s1_137);
        p_bam2.insert(137, s2_137);

        // 140: Multi-allelic tumor + is_noise() normal mutation propagation (#45)
        let mut s1_140 = Vars::new();
        s1_140.sv = "sv1-140".to_string();
        s1_140.variants.push(make_good_variant("G", "A", "G", 0.35));
        s1_140.variants.push(make_good_variant("T", "A", "T", 0.25));
        let mut s2_140 = Vars::new();
        s2_140.sv = "sv2-140".to_string();
        let mut v2_noise = make_good_variant("G", "A", "G", 0.02);
        v2_noise.mean_quality = 8.0;
        v2_noise.qstd = false;
        v2_noise.pos_coverage = 3;
        v2_noise.total_pos_cov = 120;
        v2_noise.ref_fwd_cov = 58;
        v2_noise.ref_rev_cov = 59;
        s2_140.variants.push(v2_noise);
        p_bam1.insert(140, s1_140);
        p_bam2.insert(140, s2_140);

        // 150 & 200: combineAnalysis in both print_variations_from_first_sample and
        // print_variations_from_second_sample (Germline, Skip, Keep, splice + max_rl mutation)
        let mut s1_150 = Vars::new();
        let mut v_150a = make_good_variant("-12", "ACGTACGTACGTN", "A", 0.30);
        v_150a.start_pos = 149;
        v_150a.end_pos = 161;
        v_150a.pos_coverage = 3;
        v_150a.total_pos_cov = 10;
        let mut v_150b = make_good_variant("-14", "ACGTACGTACGTACN", "A", 0.25);
        v_150b.start_pos = 149;
        v_150b.end_pos = 163;
        v_150b.pos_coverage = 3;
        v_150b.total_pos_cov = 10;
        s1_150.variants.push(v_150a);
        s1_150.variants.push(v_150b);
        p_bam1.insert(150, s1_150);
        p_bam2.insert(150, Vars::new());

        let mut s1_200 = Vars::new();
        let mut ref1_200 = Variant::new();
        ref1_200.total_pos_cov = 50;
        ref1_200.ref_fwd_cov = 25;
        ref1_200.ref_rev_cov = 25;
        s1_200.ref_variant = Some(ref1_200);
        let mut s2_200 = Vars::new();
        s2_200.sv = "sv2-200".to_string();
        let mut v_200a = make_good_variant("+ACGTACGTAC&GG", "CG", "ACGTACGTACGG", 0.30);
        v_200a.start_pos = 200;
        v_200a.end_pos = 202;
        v_200a.pos_coverage = 3;
        v_200a.total_pos_cov = 10;
        let mut v_200b = make_good_variant("-11", "ACGTACGTACGN", "A", 0.35);
        v_200b.start_pos = 200;
        v_200b.end_pos = 211;
        v_200b.pos_coverage = 3;
        v_200b.total_pos_cov = 10;
        s2_200.variants.push(v_200a);
        s2_200.variants.push(v_200b);
        p_bam1.insert(200, s1_200);
        p_bam2.insert(200, s2_200);

        let combine_pipe = |_comb_reg: &Region,
                            sp: &mut rustc_hash::FxHashSet<String>,
                            mrl: &mut i32,
                            target_pos: i32|
         -> Result<Option<Vars>, VardictError> {
            sp.insert(format!("100-{}", target_pos));
            *mrl = (*mrl).max(175);
            let mut pv = Vars::new();
            match target_pos {
                150 => {
                    let mut comb_a = make_good_variant("-12", "ACGTACGTACGTN", "A", 0.30);
                    comb_a.pos_coverage = 6; // Germline
                    comb_a.total_pos_cov = 20;
                    comb_a.genotype = Some("A/-12".to_string());
                    pv.variants.push(comb_a);
                    // "-14" absent -> Skip
                    Ok(Some(pv))
                }
                200 => {
                    let mut comb_a =
                        make_good_variant("+ACGTACGTAC&GG", "CG", "ACGTACGTACGG", 0.30);
                    comb_a.pos_coverage = 3; // Keep -> StrongLOH
                    comb_a.total_pos_cov = 20;
                    let mut comb_b = make_good_variant("-11", "ACGTACGTACGN", "A", 0.35);
                    comb_b.pos_coverage = 7; // Germline
                    comb_b.total_pos_cov = 25;
                    comb_b.genotype = Some("A/-11".to_string());
                    pv.variants.push(comb_a);
                    pv.variants.push(comb_b);
                    Ok(Some(pv))
                }
                _ => Ok(None),
            }
        };

        let borrowed_som = somatic_post_process_at(
            &cfg,
            &reg,
            &p_bam1,
            &p_bam2,
            &empty_splice,
            Some(&combine_pipe),
            100,
        );
        let owned_som = somatic_post_process_owned(
            &cfg,
            &reg,
            p_bam1,
            p_bam2,
            &empty_splice,
            Some(&combine_pipe),
            100,
        );
        assert_somatic_eq(&owned_som, &borrowed_som);
    }

    #[test]
    fn test_somatic_strict_region_bounds_compat_off() {
        use crate::types::Pos1;

        let mut cfg = Config::default();
        cfg.freq = 0.01;
        cfg.min_reads = 2;
        cfg.read_pos = 5.0;
        cfg.base_qual = 22.5;
        cfg.qratio = 1.5;

        let reg = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let empty_splice = rustc_hash::FxHashSet::default();

        let mut bam1 = JavaHashMap::new();
        let bam2 = JavaHashMap::new();

        // Pileup key 100 (inside reg), but left-anchored deletion start_pos = 99 (< reg.start)
        let mut v_left = Vars::new();
        let mut del_99 = make_good_variant("-1", "CA", "C", 0.40);
        del_99.start_pos = 99;
        del_99.end_pos = 100;
        v_left.variants.push(del_99);
        bam1.insert(100, v_left);

        // Pileup key 150 (inside reg), SNV at 150
        let mut v_in = Vars::new();
        let mut snv_150 = make_good_variant("G", "A", "G", 0.45);
        snv_150.start_pos = 150;
        snv_150.end_pos = 150;
        v_in.variants.push(snv_150);
        bam1.insert(150, v_in);

        // Pileup key 200 (inside reg), right-shifted indel at 204 (> reg.end)
        let mut v_right = Vars::new();
        let mut ins_204 = make_good_variant("+TT", "G", "GTT", 0.35);
        ins_204.start_pos = 204;
        ins_204.end_pos = 204;
        v_right.variants.push(ins_204);
        bam1.insert(200, v_right);

        let compat_default =
            somatic_post_process_at(&cfg, &reg, &bam1, &bam2, &empty_splice, None, 100);
        assert_eq!(compat_default.len(), 3);

        let mut strict_cfg = cfg.clone();
        strict_cfg.compat = crate::config::COMPAT_OFF;
        let strict_borrowed =
            somatic_post_process_at(&strict_cfg, &reg, &bam1, &bam2, &empty_splice, None, 100);
        assert_eq!(strict_borrowed.len(), 1);
        assert_eq!(strict_borrowed[0].variant1.start_pos, 150);

        let strict_owned =
            somatic_post_process_owned(&strict_cfg, &reg, bam1, bam2, &empty_splice, None, 100);
        assert_eq!(strict_owned.len(), 1);
        assert_eq!(strict_owned[0].variant1.start_pos, 150);

        // Verify paired somatic `reconcile_vd_dp` (`COMPAT_VARDICTJAVA_1_8_3` vs `COMPAT_OFF`):
        // 1) `combine_analysis` DP == 0 rescue + HIAF < 0 / ADJAF < 0 clamping + HICNT <= HICOV:
        let mut base_rescue = make_good_variant("-12", "ACGTACGTACGTN", "A", 0.15);
        base_rescue.var_type = "Deletion".to_string();
        base_rescue.pos_coverage = 3;
        base_rescue.total_pos_cov = 20;
        base_rescue.hicnt = 3;
        base_rescue.hicov = 20;
        base_rescue.hi_freq = 0.9;
        base_rescue.extra_freq = 0.8;

        let rescue_pipe = |_r: &Region,
                           _sp: &mut rustc_hash::FxHashSet<String>,
                           _mrl: &mut i32,
                           _p: i32|
         -> Result<Option<Vars>, VardictError> {
            let mut pv = Vars::new();
            let mut comb = base_rescue.clone();
            comb.pos_coverage = 6; // delta = 3 >= min_reads (2)
            comb.total_pos_cov = 20; // combined.DP == variant1.DP -> subtracted DP == 0
            comb.hicnt = 5;
            comb.hicov = 10;
            comb.hi_freq = 0.1; // weighted diff (0.6 - 2.7)/3 = -0.7 < 0
            comb.extra_freq = 0.1; // weighted diff (0.6 - 2.4)/3 = -0.6 < 0
            pv.variants.push(comb);
            Ok(Some(pv))
        };
        let mut sp_cow = std::borrow::Cow::Borrowed(&empty_splice);
        let mut mrl = 100;
        let vj_verdict = combine_analysis(
            &cfg,
            &rescue_pipe,
            &base_rescue,
            "chr1",
            150,
            &mut sp_cow,
            &mut mrl,
        );
        assert!(matches!(vj_verdict, CombineVerdict::Skip));
        let off_verdict = combine_analysis(
            &strict_cfg,
            &rescue_pipe,
            &base_rescue,
            "chr1",
            150,
            &mut sp_cow,
            &mut mrl,
        );
        match off_verdict {
            CombineVerdict::Germline(var2) => {
                assert_eq!(var2.pos_coverage, 3);
                assert_eq!(var2.total_pos_cov, 3);
                assert_eq!(var2.frequency, 1.0);
                assert_eq!(var2.hicnt, 2);
                assert_eq!(var2.hicov, 2);
                assert_eq!(var2.hi_freq, 0.0);
                assert_eq!(var2.extra_freq, 0.0);
            }
            _ => panic!("expected Germline rescue under COMPAT_OFF"),
        }

        // 2) `combine_analysis` VD > DP + HIAF > 1.0 / ADJAF > 1.0 clamping:
        let mut base_over = make_good_variant("-12", "ACGTACGTACGTN", "A", 0.10);
        base_over.var_type = "Deletion".to_string();
        base_over.pos_coverage = 2;
        base_over.total_pos_cov = 18;
        base_over.hicnt = 0;
        base_over.hicov = 18;
        base_over.hi_freq = 0.0;
        base_over.extra_freq = 0.0;

        let over_pipe = |_r: &Region,
                         _sp: &mut rustc_hash::FxHashSet<String>,
                         _mrl: &mut i32,
                         _p: i32|
         -> Result<Option<Vars>, VardictError> {
            let mut pv = Vars::new();
            let mut comb = base_over.clone();
            comb.pos_coverage = 6; // var2.pos_coverage = 4
            comb.total_pos_cov = 20; // var2.total_pos_cov = 2 < 4 (VD > DP)
            comb.hicnt = 6;
            comb.hicov = 20; // 20 - 18 = 2 < 6
            comb.hi_freq = 1.0; // (6.0 - 0.0)/4.0 = 1.5 > 1.0
            comb.extra_freq = 1.0; // 1.5 > 1.0
            pv.variants.push(comb);
            Ok(Some(pv))
        };
        match combine_analysis(
            &cfg,
            &over_pipe,
            &base_over,
            "chr1",
            150,
            &mut sp_cow,
            &mut mrl,
        ) {
            CombineVerdict::Germline(vj_v2) => {
                assert_eq!(vj_v2.pos_coverage, 4);
                assert_eq!(vj_v2.total_pos_cov, 2);
                assert!(vj_v2.frequency > 1.0);
                assert!(vj_v2.hi_freq > 1.0);
                assert!(vj_v2.extra_freq > 1.0);
            }
            _ => panic!("expected Germline under COMPAT_VARDICTJAVA_1_8_3"),
        }
        match combine_analysis(
            &strict_cfg,
            &over_pipe,
            &base_over,
            "chr1",
            150,
            &mut sp_cow,
            &mut mrl,
        ) {
            CombineVerdict::Germline(off_v2) => {
                assert_eq!(off_v2.pos_coverage, 4);
                assert_eq!(off_v2.total_pos_cov, 4);
                assert_eq!(off_v2.frequency, 1.0);
                assert_eq!(off_v2.hicnt, 6);
                assert_eq!(off_v2.hicov, 6);
                assert_eq!(off_v2.hi_freq, 1.0);
                assert_eq!(off_v2.extra_freq, 1.0);
            }
            _ => panic!("expected Germline under COMPAT_OFF"),
        }

        // 3) `build_var_for_print_from_v1` & `build_var_for_print_from_v2` negative ref_fwd clamping:
        let mut neg_v1 = Vars::new();
        let mut bad_first = Variant::new();
        bad_first.total_pos_cov = 5;
        neg_v1.variants.push(bad_first);
        let mut v1ref = Variant::new();
        v1ref.description = "A".to_string();
        v1ref.fwd = -4;
        v1ref.rev = 12;
        neg_v1.ref_variant = Some(v1ref);

        let vp1_vj = build_var_for_print_from_v1(&neg_v1, false);
        assert_eq!(vp1_vj.ref_fwd_cov, -4);
        assert_eq!(vp1_vj.ref_rev_cov, 12);
        assert_eq!(vp1_vj.total_pos_cov, 5);

        let vp1_off = build_var_for_print_from_v1(&neg_v1, true);
        assert_eq!(vp1_off.ref_fwd_cov, 0);
        assert_eq!(vp1_off.ref_rev_cov, 12);
        assert_eq!(vp1_off.total_pos_cov, 12);

        let mut neg_v2 = Vars::new();
        let mut v2_first = Variant::new();
        v2_first.total_pos_cov = 4;
        v2_first.ref_fwd_cov = -3;
        v2_first.ref_rev_cov = 9;
        neg_v2.variants.push(v2_first);

        let vp2_vj = build_var_for_print_from_v2(&neg_v2, false).expect("vp2_vj");
        assert_eq!(vp2_vj.ref_fwd_cov, -3);
        assert_eq!(vp2_vj.ref_rev_cov, 9);
        assert_eq!(vp2_vj.total_pos_cov, 4);

        let vp2_off = build_var_for_print_from_v2(&neg_v2, true).expect("vp2_off");
        assert_eq!(vp2_off.ref_fwd_cov, 0);
        assert_eq!(vp2_off.ref_rev_cov, 9);
        assert_eq!(vp2_off.total_pos_cov, 9);
    }
}
