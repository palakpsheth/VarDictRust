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
use std::collections::{BTreeMap, HashSet};

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
) -> Vec<AmpliconVariant> {
    let mut results = Vec::new();

    // Iterate over positions sorted (BTreeMap already sorted)
    for (&position, amplicon_regions) in amplicons_on_positions {
        // good variants: (variant, region_string)
        let mut gvs: Vec<(Variant, String)> = Vec::new();
        // reference variants
        let mut ref_vars: Vec<Variant> = Vec::new();
        // coverage per amplicon
        let mut vcovs: Vec<i32> = Vec::new();
        // good amplicon set: "ampliconNumber-refallele-varallele"
        let mut goodmap: HashSet<String> = HashSet::new();
        // good variants on each amplicon: amplicon_index → list of good variants
        let mut good_variants_on_amp: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();

        let mut maxcov: i32 = 0;
        let mut maxaf: f64 = 0.0;

        for &(amplicon_number, ref amp_region) in amplicon_regions {
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
                            goodmap.insert(format!(
                                "{}-{}-{}",
                                amplicon_number, tv_clone.ref_allele, tv_clone.var_allele
                            ));
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
                    results.push(AmpliconVariant {
                        variant: vref,
                        gdamp: 0,
                        tlamp: amplicon_regions.len() as i32,
                        ncamp: nocov,
                        ampflag: false,
                        region: format!("{}:{}-{}", region.chr, position, position),
                        good_variants: Vec::new(),
                        bad_variants: Vec::new(),
                    });
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
                    for &(amp_num, ref amp_reg) in amplicon_regions {
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
            let mut bad_vars: Vec<(Option<Variant>, String)> = Vec::new();
            if initial_gvscnt != total_amplicons || current_flag {
                for &(amp_num, ref amp_reg) in amplicon_regions {
                    let key = format!("{}-{}-{}", amp_num, vref.ref_allele, vref.var_allele);
                    if goodmap.contains(&key) {
                        continue;
                    }
                    if config.pileup && vref.ref_allele == vref.var_allele {
                        continue;
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
                    bad_vars,
                )
            } else {
                (Vec::new(), Vec::new())
            };

            results.push(AmpliconVariant {
                variant: out_variant,
                gdamp,
                tlamp,
                ncamp: nocov,
                ampflag: current_flag,
                region: region_str,
                good_variants: debug_good,
                bad_variants: debug_bad,
            });
        }
    }

    // VCF mode: same var2vcf_valid.pl selection as simple mode (:41,129-137) —
    // bucket rows by the RAW .var start column, stable af-desc sort per bucket,
    // emit only the top row ($ALL = 1 without -A). See postprocess/simple.rs.
    if !config.emit_intermediate {
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
        results = selected;
    }

    results
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
        );

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
        );

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
        );

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
        );

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
        );

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
        );

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
}
