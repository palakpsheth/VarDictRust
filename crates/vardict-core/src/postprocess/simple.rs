//! Simple post-processing module for germline (single-sample) variants.
//! Ports SimplePostProcessModule.java.

use crate::config::Config;
use crate::regions::Region;
use crate::variations::{Variant, Vars};
use std::collections::{BTreeMap, HashSet};

/// Build the null-record row Java prints via `SimpleOutputVariant(null, …)`
/// (SimpleOutputVariant.java:50-52): start/end = the alignedVariants map key,
/// everything else at printer defaults — bias "0;0" (:24), genotype "" (:22),
/// empty flank seqs (OutputVariant.java:22-23) — rendered by table.rs when
/// `null_record` is set.
fn make_null_record(position: i32, sv: &str) -> Variant {
    let mut out = Variant::new();
    out.start_pos = position;
    out.end_pos = position;
    out.bias = "0;0".to_string();
    out.sv = sv.to_string();
    out.null_record = true;
    out
}

/// Single sample mode variant filtering and postprocessing.
///
/// Filters the aligned variants list using is_good_var and pileup rules,
/// then applies complex variant sequence adjustment.
pub fn simple_post_process(
    config: &Config,
    region: &Region,
    aligned_variants: &BTreeMap<i32, Vars>,
    splice: &HashSet<String>,
) -> Vec<Variant> {
    let mut vrefs = Vec::new();

    // Iterate over positions sorted by key
    let mut sorted_positions: Vec<i32> = aligned_variants.keys().copied().collect();
    sorted_positions.sort();

    for &position in &sorted_positions {
        if let Some(variants_on_position) = aligned_variants.get(&position) {
            let sv = &variants_on_position.sv;

            if sv.is_empty() && (position < region.start.0 as i32 || position > region.end.0 as i32)
            {
                continue;
            }

            if variants_on_position.variants.is_empty() {
                if !config.pileup {
                    continue;
                }
                if let Some(ref vref) = variants_on_position.ref_variant {
                    let mut vref = vref.clone();
                    vref.var_type = String::new();
                    vref.sv = sv.clone();
                    vrefs.push(vref);
                } else {
                    // ports SimplePostProcessModule.java:55-58: under -p an entry
                    // with no variants and no reference variant still prints, via
                    // `SimpleOutputVariant(null, region, sv, position)` whose fields
                    // keep the printer defaults (genotype "", bias "0;0", flank seqs
                    // ""). CLAUDE.md §5: live-oracle output is replicated as-is and
                    // tracked upstream, never silently fixed.
                    vrefs.push(make_null_record(position, sv));
                }
            } else {
                let vvar = &variants_on_position.variants;
                let mut valid_vars = Vec::new();
                for vref in vvar {
                    let mut vref = vref.clone();
                    vref.sv = sv.clone();
                    if vref.ref_allele.contains('N') {
                        continue;
                    }
                    if vref.ref_allele == vref.var_allele && !config.pileup {
                        continue;
                    }
                    if vref.start_pos != position && config.pileup && vvar.len() == 1 {
                        if let Some(ref_var) = &variants_on_position.ref_variant {
                            let mut ref_var = ref_var.clone();
                            ref_var.var_type = String::new();
                            ref_var.sv = sv.clone();
                            valid_vars.push(ref_var);
                        } else {
                            // ports SimplePostProcessModule.java:75-79: Java prints a
                            // SimpleOutputVariant wrapping null immediately (the
                            // null-record row at the map key), then continues with a
                            // default `new Variant()` whose startPosition/endPosition
                            // stay 0 — that 0 leaks into .var columns 4-5 (the
                            // upstream pos-0 quirk; replicated per §5, see
                            // DIVERGENCES.md "Known non-bugs").
                            valid_vars.push(make_null_record(position, sv));
                            let mut ref_var = Variant::new();
                            ref_var.var_type = String::new();
                            ref_var.sv = sv.clone();
                            valid_vars.push(ref_var);
                        }
                    }
                    vref.var_type = vref.var_type();

                    let is_good = vref.is_good_var(
                        variants_on_position.ref_variant.as_ref(),
                        &vref.var_type,
                        splice,
                        config,
                    );
                    if !is_good && !config.pileup {
                        continue;
                    }
                    valid_vars.push(vref);
                }

                vrefs.extend(valid_vars);
            }
        }
    }

    // Apply adjustments
    for vref in &mut vrefs {
        if vref.var_type == "Complex" {
            vref.adj_complex();
        }
        if config.crispr_cutting_site.unwrap_or(0) == 0 {
            vref.crispr = 0;
        }
    }

    // VCF mode: ports var2vcf_valid.pl:41,129-137 — bucket rows by the RAW .var
    // start column ($hash{$chr}->{$a[3]}, column 3), walk buckets in ascending
    // numeric position (:129), sort each bucket by column 14 (af) descending
    // (:135), and emit only the top row ($ALL = 1 without -A, :137). The stable
    // sort keeps af ties in .var stream order, like Perl's default mergesort.
    // Bucketing must use the raw start (what simple_36 writes), NOT the
    // ToVarsBuilder position key: indels and SNVs can share one key while
    // having different raw starts, and different keys can share a raw start.
    if !config.emit_intermediate {
        let mut buckets: BTreeMap<i32, Vec<Variant>> = BTreeMap::new();
        for v in vrefs {
            buckets.entry(v.start_pos).or_default().push(v);
        }
        let mut selected = Vec::with_capacity(buckets.len());
        for (_, mut rows) in buckets {
            rows.sort_by(|a, b| {
                b.frequency
                    .partial_cmp(&a.frequency)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            selected.extend(rows.into_iter().take(1));
        }
        vrefs = selected;
    }

    vrefs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Pos1;

    #[test]
    fn test_simple_post_process_filters() {
        let mut config = Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;
        config.pileup = false;

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut aligned_variants = BTreeMap::new();
        let splice = HashSet::new();

        // 1. Position 120: Valid variant
        let mut vars_120 = Vars::new();
        let mut v_120 = Variant::new();
        v_120.ref_allele = "A".to_string();
        v_120.var_allele = "G".to_string();
        v_120.frequency = 0.05;
        v_120.hicnt = 5;
        v_120.mean_position = 10.0;
        v_120.mean_quality = 30.0;
        v_120.hq_lq_ratio = 2.0;
        v_120.mean_mapping_quality = 30.0;
        v_120.start_pos = 120;
        v_120.end_pos = 120;
        vars_120.variants.push(v_120);
        aligned_variants.insert(120, vars_120);

        // 2. Position 130: Low frequency variant (should be filtered out since pileup=false)
        let mut vars_130 = Vars::new();
        let mut v_130 = Variant::new();
        v_130.ref_allele = "T".to_string();
        v_130.var_allele = "C".to_string();
        v_130.frequency = 0.002; // < 0.01
        v_130.hicnt = 5;
        v_130.mean_position = 10.0;
        v_130.mean_quality = 30.0;
        v_130.hq_lq_ratio = 2.0;
        v_130.mean_mapping_quality = 30.0;
        v_130.start_pos = 130;
        v_130.end_pos = 130;
        vars_130.variants.push(v_130);
        aligned_variants.insert(130, vars_130);

        let filtered = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].start_pos, 120);

        // Now run with pileup = true. Low frequency variant should be kept.
        config.pileup = true;
        let filtered_pileup = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(filtered_pileup.len(), 2);
    }

    fn mk_variant(start: i32, end: i32, ref_al: &str, alt: &str, af: f64) -> Variant {
        let mut v = Variant::new();
        v.start_pos = start;
        v.end_pos = end;
        v.ref_allele = ref_al.to_string();
        v.var_allele = alt.to_string();
        v.frequency = af;
        v.hicnt = 5;
        v.mean_position = 10.0;
        v.mean_quality = 30.0;
        v.hq_lq_ratio = 2.0;
        v.mean_mapping_quality = 30.0;
        v
    }

    /// Two rows share the raw .var start column but sit in different
    /// ToVarsBuilder groups: Perl's bucketing (var2vcf_valid.pl:41,137)
    /// suppresses the lower-AF row in VCF mode and keeps both in intermediate.
    #[test]
    fn test_vcf_mode_top_af_per_raw_start() {
        let mut config = Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;

        let region = Region {
            chr: "20".to_string(),
            start: Pos1(20159100),
            end: Pos1(20159300),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        // Different group keys, identical raw start column.
        let mut g1 = Vars::new();
        g1.variants
            .push(mk_variant(20159222, 20159222, "A", "AT", 0.6));
        let mut g2 = Vars::new();
        g2.variants
            .push(mk_variant(20159222, 20159223, "AT", "A", 0.1));
        let mut aligned_variants = BTreeMap::new();
        aligned_variants.insert(20159221, g1);
        aligned_variants.insert(20159222, g2);

        let splice = HashSet::new();
        let kept = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].var_allele, "AT");

        // Intermediate mode emits every row.
        config.emit_intermediate = true;
        let all = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(all.len(), 2);
    }

    /// An indel and an SNV can share one ToVarsBuilder group while having
    /// different raw starts: both are separate Perl buckets and both must be
    /// emitted (regression for the chr20:20179027 CT>C loss, DIVERGENCES P1a).
    #[test]
    fn test_vcf_mode_keeps_adjacent_raw_starts() {
        let mut config = Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;

        let region = Region {
            chr: "20".to_string(),
            start: Pos1(20178900),
            end: Pos1(20179200),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        // Same group key, different raw start columns.
        let mut vars = Vars::new();
        vars.variants
            .push(mk_variant(20179028, 20179028, "T", "C", 0.9));
        vars.variants
            .push(mk_variant(20179027, 20179028, "CT", "C", 0.4));
        let mut aligned_variants = BTreeMap::new();
        aligned_variants.insert(20179028, vars);

        let splice = HashSet::new();
        let kept = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(kept.len(), 2);
        assert_eq!(kept[0].start_pos, 20179027);
        assert_eq!(kept[1].start_pos, 20179028);

        // The top-AF row wins its bucket regardless of insertion order.
        let mut reversed = Vars::new();
        reversed
            .variants
            .push(mk_variant(20179027, 20179028, "CT", "C", 0.4));
        reversed
            .variants
            .push(mk_variant(20179028, 20179028, "T", "C", 0.9));
        let mut aligned_rev = BTreeMap::new();
        aligned_rev.insert(20179028, reversed);
        let kept_rev = simple_post_process(&config, &region, &aligned_rev, &splice);
        assert_eq!(kept_rev[1].var_allele, "C");
    }

    fn pileup_config() -> Config {
        let mut config = Config::default();
        config.freq = 0.0;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;
        // hard_clip case semantics: -p with intermediate emission
        config.pileup = true;
        config.emit_intermediate = true;
        config
    }

    /// ports SimplePostProcessModule.java:73-82: an entry keyed at the map key
    /// whose single variant's raw start differs and which has no reference
    /// variant emits THREE rows — the null-record row at (key,key) with printer
    /// defaults, a default `new Variant()` at (0,0), and the variant itself.
    /// Regression for the live-oracle 6447 + pos-0 rows (P7/P8a).
    #[test]
    fn test_pileup_null_record_fallback_branch() {
        let config = pileup_config();
        let region = Region {
            chr: "test".to_string(),
            start: Pos1(6434),
            end: Pos1(6564),
            gene: "testbed".to_string(),
            amplicon: None,
        };
        let mut vars = Vars::new();
        vars.variants.push(mk_variant(6449, 6449, "G", "AGTT", 1.0));
        let mut aligned_variants = BTreeMap::new();
        aligned_variants.insert(6447, vars);
        let splice = HashSet::new();

        let out = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(out.len(), 3);

        let null_row = &out[0];
        assert!(null_row.null_record);
        assert_eq!(null_row.start_pos, 6447);
        assert_eq!(null_row.end_pos, 6447);
        assert_eq!(null_row.bias, "0;0");
        assert_eq!(null_row.genotype, None);

        let zero_row = &out[1];
        assert!(!zero_row.null_record);
        assert_eq!(zero_row.start_pos, 0);
        assert_eq!(zero_row.end_pos, 0);
        assert_eq!(zero_row.genotype, None);

        assert_eq!(&out[2].start_pos, &6449);
    }

    /// ports SimplePostProcessModule.java:50-58: under -p an entry with no
    /// variants and no reference variant still prints one null-record row.
    #[test]
    fn test_pileup_null_record_empty_entry() {
        let config = pileup_config();
        let region = Region {
            chr: "test".to_string(),
            start: Pos1(6434),
            end: Pos1(6564),
            gene: "testbed".to_string(),
            amplicon: None,
        };
        let mut aligned_variants = BTreeMap::new();
        aligned_variants.insert(6447, Vars::new());
        let splice = HashSet::new();

        let out = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(out.len(), 1);
        assert!(out[0].null_record);
        assert_eq!(out[0].start_pos, 6447);
    }

    /// Without -p nothing changes: empty entries are skipped and the fallback
    /// branch never fires (corpus proves non-pileup modes unaffected).
    #[test]
    fn test_non_pileup_emits_no_null_records() {
        let mut config = pileup_config();
        config.pileup = false;
        let region = Region {
            chr: "test".to_string(),
            start: Pos1(6434),
            end: Pos1(6564),
            gene: "testbed".to_string(),
            amplicon: None,
        };
        let mut vars = Vars::new();
        vars.variants.push(mk_variant(6449, 6449, "G", "AGTT", 1.0));
        let mut aligned_variants = BTreeMap::new();
        aligned_variants.insert(6447, vars);
        aligned_variants.insert(6500, Vars::new());
        let splice = HashSet::new();

        let out = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start_pos, 6449);
        assert!(!out[0].null_record);
    }
}
