//! Somatic post-processing module — determines somatic STATUS for paired samples.
//!
//! Ports `SomaticPostProcessModule.java` from VarDictJava v1.8.3.
//! Compares variants from two BAM files (tumor=BAM1, normal=BAM2) and
//! determines the somatic STATUS (StrongSomatic, LikelySomatic, Germline,
//! LikelyLOH, AFDiff, SampleSpecific, Deletion, StrongLOH).

use crate::config::Config;
use crate::regions::Region;
use crate::types::VardictError;
use crate::variations::{Variant, Vars};
use std::collections::HashSet;

/// Callback that runs the full pipeline (pileup→realign→SV→tovars) on both BAMs
/// combined over a given region, returning aligned variants keyed by position.
pub type PipelineFn<'a> =
    dyn Fn(&Region) -> Result<crate::javamap::JavaHashMap<i32, Vars>, VardictError> + 'a;

// STATUS constants matching VarDictJava's SomaticPostProcessModule
const STRONG_SOMATIC: &str = "StrongSomatic";

/// Ports `Patterns.MINUS_NUM_NUM` = `-\d\d` (unanchored find).
fn minus_num_num_find(desc: &str) -> bool {
    let b = desc.as_bytes();
    b.len() >= 3
        && b.windows(3)
            .any(|w| w[0] == b'-' && w[1].is_ascii_digit() && w[2].is_ascii_digit())
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
}

/// Performs somatic analysis of variants from both samples (BAM1 and BAM2).
///
/// Merges all positions from both samples' aligned variants, then for each
/// position in the region determines the somatic STATUS.
///
/// Ports `SomaticPostProcessModule.accept()`.
pub fn somatic_post_process(
    config: &Config,
    region: &Region,
    aligned_variants_bam1: &crate::javamap::JavaHashMap<i32, Vars>,
    aligned_variants_bam2: &crate::javamap::JavaHashMap<i32, Vars>,
    splice: &HashSet<String>,
    pipeline: Option<&PipelineFn<'_>>,
    max_read_length: i32,
) -> Vec<SomaticVariant> {
    let mut results = Vec::new();

    // Merge all positions from both samples and sort
    let mut all_positions: Vec<i32> = aligned_variants_bam1
        .keys()
        .chain(aligned_variants_bam2.keys())
        .copied()
        .collect();
    all_positions.sort_unstable();
    all_positions.dedup();

    for &position in &all_positions {
        if (position as u32) < region.start.0 || (position as u32) > region.end.0 {
            continue;
        }
        let v1 = aligned_variants_bam1.get(&position);
        let v2 = aligned_variants_bam2.get(&position);

        match (v1, v2) {
            (None, None) => continue,
            (None, Some(v2)) => {
                // No coverage for sample 1 → Deletion
                calling_for_one_sample(config, v2, true, DELETION, splice, &mut results);
            }
            (Some(v1), None) => {
                // No coverage for sample 2 → SampleSpecific
                calling_for_one_sample(config, v1, false, SAMPLE_SPECIFIC, splice, &mut results);
            }
            (Some(v1), Some(v2)) => {
                // Both samples have coverage
                calling_for_both_samples(
                    config,
                    position,
                    v1,
                    v2,
                    &region.chr,
                    splice,
                    &mut results,
                    pipeline,
                    max_read_length,
                );
            }
        }
    }

    results
}

/// Implements variations analysis from one sample.
///
/// Ports `SomaticPostProcessModule.callingForOneSample()`.
fn calling_for_one_sample(
    config: &Config,
    variants: &Vars,
    is_first_cover: bool,
    var_label: &str,
    splice: &HashSet<String>,
    results: &mut Vec<SomaticVariant>,
) {
    if variants.variants.is_empty() {
        return;
    }
    for variant in &variants.variants {
        if variant.ref_allele == variant.var_allele {
            continue;
        }
        let mut v = variant.clone();
        v.var_type = v.var_type();
        if !v.is_good_var(variants.ref_variant.as_ref(), &v.var_type, splice, config) {
            continue;
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
            results.push(SomaticVariant {
                variant1: tumor,
                variant2: Some(v),
                status: var_label.to_string(),
                sv1: String::new(),
                sv2: variants.sv.clone(),
                end_is_variant2: true,
            });
        } else {
            // Only BAM1 has coverage → SampleSpecific
            results.push(SomaticVariant {
                variant1: v.clone(),
                variant2: None,
                status: var_label.to_string(),
                sv1: variants.sv.clone(),
                sv2: String::new(),
                end_is_variant2: false,
            });
        }
    }
}

/// Run analysis of variations from BAM1 and BAM2.
///
/// Ports `SomaticPostProcessModule.callingForBothSamples()`.
fn calling_for_both_samples(
    config: &Config,
    position: i32,
    v1: &Vars,
    v2: &Vars,
    chr: &str,
    splice: &HashSet<String>,
    results: &mut Vec<SomaticVariant>,
    pipeline: Option<&PipelineFn<'_>>,
    max_read_length: i32,
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
    v1: &Vars,
    v2: &Vars,
    chr: &str,
    splice: &HashSet<String>,
    results: &mut Vec<SomaticVariant>,
    pipeline: Option<&PipelineFn<'_>>,
    max_read_length: i32,
) {
    let mut number_of_processed_variation: usize = 0;
    // Ports SomaticPostProcessModule.java:172-186 — in Java, `v2.variants` and
    // `v2.varDescriptionStringToVariants` hold references to the same `Variant`
    // objects. When `determinateType` calls `v2nt.isNoise()`, it mutates `v2nt`
    // in place (`totalPosCoverage -= positionCoverage`, etc.), so a subsequent
    // tumor variant at the same position that falls into the `v2nt == null`
    // branch reads the mutated `v2.variants.get(0).totalPosCoverage`.
    let mut v2 = v2.clone();

    // Process variants from sample 1 while they pass isGoodVar
    while number_of_processed_variation < v1.variants.len() {
        let vref_orig = &v1.variants[number_of_processed_variation];
        let vartype_check = vref_orig.var_type();
        if !vref_orig.is_good_var(v1.ref_variant.as_ref(), &vartype_check, splice, config) {
            break;
        }

        let mut vref = vref_orig.clone();
        if vref.ref_allele == vref.var_allele {
            number_of_processed_variation += 1;
            continue;
        }

        let nt = vref.description.clone();
        vref.var_type = vref.var_type();
        if vref.var_type == COMPLEX {
            vref.adj_complex();
        }

        // Look up same descriptionString in sample2
        let v2nt = v2.var_desc_to_variant.get(&nt).cloned();

        if let Some(mut v2nt_clone) = v2nt {
            let type_str = determinate_type(config, &v2, &vref, &mut v2nt_clone, splice);
            if number_of_processed_variation + 1 < v1.variants.len() {
                if let Some(slot) = v2.var_desc_to_variant.get_mut(&nt) {
                    *slot = v2nt_clone.clone();
                }
                for v2_item in &mut v2.variants {
                    if v2_item.description == nt {
                        *v2_item = v2nt_clone.clone();
                        break;
                    }
                }
            }
            results.push(SomaticVariant {
                variant1: vref,
                variant2: Some(v2nt_clone),
                status: type_str,
                sv1: v1.sv.clone(),
                sv2: v2.sv.clone(),
                end_is_variant2: true,
            });
        } else {
            // sample 1 only, should be strong somatic (SomaticPostProcessModule.java:177)
            let var_for_print = build_var_for_print_from_v2(&v2);
            let mut type_str = STRONG_SOMATIC.to_string();
            let mut combined_v2nt: Option<Variant> = None;

            // combineAnalysis gate (:193-197): not SNV, description longer than 10
            // or matching MINUS_NUM_NUM (-\d\d), low pos coverage, and not an SV
            // ('<') → re-run the pipeline on both BAMs combined
            if vref.var_type != SNV
                && (nt.len() > 10 || minus_num_num_find(&nt))
                && vref.pos_coverage < config.min_reads as i32 + 3
                && !nt.contains('<')
            {
                if let Some(pipeline_fn) = pipeline {
                    match combine_analysis(
                        config,
                        pipeline_fn,
                        &vref,
                        chr,
                        position,
                        max_read_length,
                    ) {
                        CombineVerdict::Skip => {
                            // FALSE → skip this variant (:202-204)
                            number_of_processed_variation += 1;
                            continue;
                        }
                        CombineVerdict::Germline(v) => {
                            type_str = GERMLINE.to_string();
                            combined_v2nt = Some(*v);
                        }
                        CombineVerdict::Keep => {} // "" keeps StrongSomatic
                    }
                }
            }

            if type_str == STRONG_SOMATIC {
                // SomaticOutputVariant(vref, vref, vref, varForPrint, …) (:212)
                results.push(SomaticVariant {
                    variant1: vref,
                    variant2: var_for_print,
                    status: STRONG_SOMATIC.to_string(),
                    sv1: v1.sv.clone(),
                    sv2: v2.sv.clone(),
                    end_is_variant2: false,
                });
            } else {
                // SomaticOutputVariant(vref, vref, vref, v2nt, …) (:215): header and
                // end fields stay with vref; the var2 columns carry the
                // combine-adjusted placeholder.
                results.push(SomaticVariant {
                    variant1: vref,
                    variant2: combined_v2nt,
                    status: type_str,
                    sv1: v1.sv.clone(),
                    sv2: v2.sv.clone(),
                    end_is_variant2: false,
                });
            }
        }
        number_of_processed_variation += 1;
    }

    // If no variants from sample1 were processed, check sample2 variants
    if number_of_processed_variation == 0 {
        if v2.variants.is_empty() {
            return;
        }
        for v2var_orig in &v2.variants {
            let mut v2var = v2var_orig.clone();
            v2var.var_type = v2var.var_type();
            if !v2var.is_good_var(v2.ref_variant.as_ref(), &v2var.var_type, splice, config) {
                continue;
            }
            // potential LOH
            let nt = v2var.description.clone();
            let v1nt = v1.var_desc_to_variant.get(&nt);

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
                    sv2: v2.sv.clone(),
                    end_is_variant2: true,
                });
            } else {
                if v2var.ref_allele == v2var.var_allele {
                    continue;
                }
                // Build a placeholder variant1 from sample1's reference info
                let mut var_for_print = build_var_for_print_from_v1(v1);
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
                results.push(SomaticVariant {
                    variant1: var_for_print,
                    variant2: Some(v2var_out),
                    status: STRONG_LOH.to_string(),
                    sv1: String::new(),
                    sv2: v2.sv.clone(),
                    end_is_variant2: true,
                });
            }
        }
    }
}

/// Analyse and output variations from BAM2 based on variations from BAM1.
///
/// Ports `SomaticPostProcessModule.printVariationsFromSecondSample()`
/// (SomaticPostProcessModule.java:285-336).
fn print_variations_from_second_sample(
    config: &Config,
    v1: &Vars,
    v2: &Vars,
    chr: &str,
    position: i32,
    splice: &HashSet<String>,
    results: &mut Vec<SomaticVariant>,
    pipeline: Option<&PipelineFn<'_>>,
    max_read_length: i32,
) {
    for v2var_orig in &v2.variants {
        if v2var_orig.ref_allele == v2var_orig.var_allele {
            continue;
        }
        let mut v2var = v2var_orig.clone();
        v2var.var_type = v2var.var_type();
        if !v2var.is_good_var(v2.ref_variant.as_ref(), &v2var.var_type, splice, config) {
            continue;
        }

        // combineAnalysis gate (:301-302): low normal coverage, not an SV, long
        // description or a -\d\d pattern (MINUS_NUM_NUM, Patterns.java:18)
        let description = &v2var.description;
        let should_combine = v2var.pos_coverage < config.min_reads as i32 + 3
            && !description.contains('<')
            && (description.len() > 10 || minus_num_num_find(description));
        let mut type_str = STRONG_LOH.to_string();
        let mut tumor_from_combine: Option<Variant> = None;
        if should_combine {
            if let Some(pipeline_fn) = pipeline {
                match combine_analysis(config, pipeline_fn, &v2var, chr, position, max_read_length)
                {
                    CombineVerdict::Skip => continue, // :313-315
                    CombineVerdict::Germline(v) => {
                        type_str = GERMLINE.to_string(); // :318-319
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
                t.duprate = 0.0; // fresh Variant in Java; never set by combineAnalysis
                t
            }
            None => {
                // Fallback: StrongLOH with reference placeholder (:322-327).
                // Header/allele fields are re-sourced from v2var exactly as the
                // SomaticOutputVariant constructor does.
                let mut var_for_print = if let Some(ref v1ref) = v1.ref_variant {
                    v1ref.clone()
                } else {
                    Variant::new()
                };
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
                var_for_print
            }
        };

        results.push(SomaticVariant {
            variant1: tumor_var,
            variant2: Some(v2var),
            status: type_str,
            sv1: String::new(),
            sv2: v2.sv.clone(),
            end_is_variant2: true,
        });
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
    pipeline: &PipelineFn<'_>,
    variant1: &Variant,
    chr: &str,
    position: i32,
    max_read_length: i32,
) -> CombineVerdict {
    if config.verbose {
        eprintln!("Start Combine {} {}", position, variant1.description);
    }
    // Don't do it for structural variants (:394)
    if variant1.end_pos - variant1.start_pos > config.sv_min_len as i32 {
        return CombineVerdict::Keep;
    }

    // Region spans the variant ± maxReadLength (:398)
    let region = Region {
        chr: chr.to_string(),
        start: crate::types::Pos1(variant1.start_pos.saturating_sub(max_read_length).max(1) as u32),
        end: crate::types::Pos1(variant1.end_pos.saturating_add(max_read_length) as u32),
        gene: String::new(),
        amplicon: None,
    };

    // Java: a pipeline exception escapes combineAnalysis to accept()'s per-position
    // try/catch (:89-91) → printExceptionAndContinue prints it and processing
    // continues from the next position; abort only after >10 continued exceptions
    // (Utils.java:261-278, Configuration.java:300). The port keeps the FALSE-skip
    // control flow but emits the error instead of swallowing it.
    let combined_vars = match pipeline(&region) {
        Ok(vars) => vars,
        Err(e) => {
            if config.verbose {
                eprintln!(
                    "There was Exception while processing position on {}. The processing will be continued from the next position.: {}",
                    position, e
                );
            }
            return CombineVerdict::Skip;
        }
    };

    // Look up by ITERATION position key and description in the map (:408).
    let Some(combined_var) = combined_vars
        .get(&position)
        .and_then(|pv| pv.var_desc_to_variant.get(variant1.description.as_str()))
    else {
        // vref == null → FALSE (:475)
        return CombineVerdict::Skip;
    };

    let base_pc = variant1.pos_coverage;
    let comb_pc = combined_var.pos_coverage;

    if comb_pc - base_pc >= config.min_reads as i32 {
        // Enough extra evidence in the other sample: fill variant2 by subtraction
        // (combined − variant1), clamped at zero (:414-436)
        let mut var2 = Variant::new();
        var2.total_pos_cov = (combined_var.total_pos_cov - variant1.total_pos_cov).max(0);
        var2.pos_coverage = (comb_pc - base_pc).max(0);
        var2.ref_fwd_cov = (combined_var.ref_fwd_cov - variant1.ref_fwd_cov).max(0);
        var2.ref_rev_cov = (combined_var.ref_rev_cov - variant1.ref_rev_cov).max(0);
        var2.fwd = (combined_var.fwd - variant1.fwd).max(0);
        var2.rev = (combined_var.rev - variant1.rev).max(0);

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

        if var2.total_pos_cov <= 0 {
            return CombineVerdict::Skip; // :456-458
        }

        var2.frequency = var2.pos_coverage as f64 / var2.total_pos_cov as f64; // :460
                                                                               // Can't back-calculate; copy from variant1 (:461)
        var2.hq_lq_ratio = variant1.hq_lq_ratio;
        var2.genotype = combined_var.genotype.clone(); // :462
                                                       // Recalculate strand bias from subtracted counts (:463-464); Java strandBias
                                                       // uses conf.bias (default 0.05) and conf.minBiasReads (-B, default 2)
        let ref_bias = crate::variations::util::strand_bias(
            var2.ref_fwd_cov,
            var2.ref_rev_cov,
            0.05,
            config.min_reads_bias as i32,
        );
        let var_bias = crate::variations::util::strand_bias(
            var2.fwd,
            var2.rev,
            0.05,
            config.min_reads_bias as i32,
        );
        var2.bias = format!("{};{}", ref_bias, var_bias);
        CombineVerdict::Germline(Box::new(var2))
    } else if comb_pc < base_pc - 2 {
        // Combined run produced less than the base sample alone (:466-470)
        if config.verbose {
            // Java config.y (-y); vardict.pl:657 prints the same line.
            eprintln!(
                "Combine produce less: {} {} {} {} {}",
                chr, position, variant1.description, comb_pc, base_pc
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
    variants: &Vars,
    standard_variant: &Variant,
    variant_to_compare: &mut Variant,
    splice: &HashSet<String>,
) -> String {
    let vtc = variant_to_compare;
    let vartype = &standard_variant.var_type;

    let type_str = if vtc.is_good_var(variants.ref_variant.as_ref(), vartype, splice, config) {
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
fn build_var_for_print_from_v2(v2: &Vars) -> Option<Variant> {
    if !v2.variants.is_empty() {
        let v2r = &v2.variants[0];
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
    }
}

/// Build a placeholder Variant from v1's first variant or reference for StrongLOH output.
///
/// Ports the `varForPrint` construction in `printVariationsFromFirstSample` when
/// numberOfProcessedVariation == 0 and variant is not found in sample1.
fn build_var_for_print_from_v1(v1: &Vars) -> Variant {
    let tcov = if !v1.variants.is_empty() {
        let v1var = &v1.variants[0];
        if v1var.total_pos_cov != 0 {
            v1var.total_pos_cov
        } else {
            0
        }
    } else {
        0
    };

    let (fwd, rev, genotype) = if let Some(ref v1ref) = v1.ref_variant {
        (v1ref.fwd, v1ref.rev, {
            if !v1.variants.is_empty() {
                v1.variants[0].genotype.clone()
            } else {
                Some(format!("{}/{}", v1ref.description, v1ref.description))
            }
        })
    } else {
        (0, 0, {
            if !v1.variants.is_empty() {
                v1.variants[0].genotype.clone()
            } else {
                Some("N/N".to_string())
            }
        })
    };

    let mut var_for_print = Variant::new();
    var_for_print.total_pos_cov = tcov;
    var_for_print.ref_fwd_cov = fwd;
    var_for_print.ref_rev_cov = rev;
    var_for_print.genotype = genotype;
    var_for_print
}

#[cfg(test)]
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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

        let bam1 = crate::javamap::JavaHashMap::new();

        let mut bam2 = crate::javamap::JavaHashMap::new();
        let mut vars = Vars::new();
        vars.variants.push(make_good_variant("G", "A", "G", 0.10));
        vars.sv = String::new();
        bam2.insert(150, vars);

        let results = somatic_post_process(&config, &region, &bam1, &bam2, &splice, None, 150);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].status, DELETION);
        assert!(results[0].variant2.is_some());
    }

    #[test]
    fn test_germline() {
        let config = test_config();
        let region = test_region();
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
    ) -> impl Fn(&Region) -> Result<crate::javamap::JavaHashMap<i32, Vars>, crate::types::VardictError>
           + '_ {
        move |_region: &Region| {
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
        assert!(matches!(
            combine_analysis(&config, &pipeline, &v, "chr1", 150, 150),
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
        assert!(matches!(
            combine_analysis(&config, &pipeline, &v, "chr1", 150, 150),
            CombineVerdict::Keep
        ));
    }

    #[test]
    fn test_combine_keep_when_delta_below_minr() {
        // Java :471-473: found, delta < minr, and not < base-2 → ""
        let config = test_config();
        let v = make_insertion_variant(); // pc=3
        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 4, 20); // delta=1 < minr=2
        assert!(matches!(
            combine_analysis(&config, &pipeline, &v, "chr1", 150, 150),
            CombineVerdict::Keep
        ));
    }

    #[test]
    fn test_combine_less_evidence_is_skip() {
        // Java :466-470: comb_pc < base_pc - 2 → FALSE ("Combine produce less")
        let config = test_config();
        let v = make_insertion_variant(); // pc=3
        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 0, 5); // 0 < 1
        assert!(matches!(
            combine_analysis(&config, &pipeline, &v, "chr1", 150, 150),
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
        assert!(matches!(
            combine_analysis(&config, &pipeline, &v, "chr1", 150, 150),
            CombineVerdict::Keep
        ));
    }

    #[test]
    fn test_combine_germline_subtracts_other_sample() {
        // Java :413-465: delta >= minr → variant2 filled by subtraction → GERMLINE
        let config = test_config();
        let mut v = make_insertion_variant();
        v.genotype = Some("0/1".to_string());
        let pipeline = combined_pipeline(150, Some("+GAACCACTGAGACAGGA"), 6, 20);
        match combine_analysis(&config, &pipeline, &v, "chr1", 150, 150) {
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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
        let splice = HashSet::new();

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
}
