//! Simple post-processing module for germline (single-sample) variants.
//! Ports SimplePostProcessModule.java.

use crate::config::Config;
use crate::regions::Region;
use crate::variations::{Variant, Vars};
use std::collections::BTreeMap;

/// Build the null-record row Java prints via `SimpleOutputVariant(null, …)`
/// (SimpleOutputVariant.java:50-52): start/end = the alignedVariants map key,
/// everything else at printer defaults — bias "0;0" (:24), genotype "" (:22),
/// empty flank seqs (OutputVariant.java:22-23) — rendered by table.rs when
/// `null_record` is set (`sv` populated by caller).
fn make_null_record(position: i32) -> Variant {
    let mut out = Variant::new();
    out.start_pos = position;
    out.end_pos = position;
    out.bias.clear();
    out.bias.push_str("0;0");
    out.null_record = true;
    out
}

#[inline]
fn single_region_bounds_i32(region: &Region) -> [(i32, i32); 1] {
    [(
        i32::try_from(region.start.0.max(1)).unwrap_or(i32::MAX),
        i32::try_from(region.end.0).unwrap_or(i32::MAX),
    )]
}

/// Single sample mode variant filtering and postprocessing.
///
/// Filters the aligned variants list using is_good_var and pileup rules,
/// then applies complex variant sequence adjustment.
pub fn simple_post_process(
    config: &Config,
    region: &Region,
    aligned_variants: &crate::javamap::JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
) -> Vec<Variant> {
    let single_bounds = single_region_bounds_i32(region);
    let strict_bounds = config
        .compat
        .strict_region_bounds
        .then_some(single_bounds.as_slice());
    if config.compat.coordinate_sorted_vars {
        let mut keys: Vec<i32> = aligned_variants.keys().copied().collect();
        keys.sort_unstable();
        simple_post_process_core(
            config,
            region,
            keys.into_iter().filter_map(|pos| {
                aligned_variants
                    .get(&pos)
                    .map(|vars| (pos, std::borrow::Cow::Borrowed(vars)))
            }),
            splice,
            strict_bounds,
        )
    } else {
        simple_post_process_core(
            config,
            region,
            aligned_variants
                .iter()
                .map(|(&pos, vars)| (pos, std::borrow::Cow::Borrowed(vars))),
            splice,
            strict_bounds,
        )
    }
}

/// Owned fast-path counterpart of [`simple_post_process`], moving finalized [`Variant`]
/// and SV strings directly out of `aligned_variants` without cloning.
#[cfg(test)]
pub(crate) fn simple_post_process_owned(
    config: &Config,
    region: &Region,
    aligned_variants: crate::javamap::JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
) -> Vec<Variant> {
    let single_bounds = single_region_bounds_i32(region);
    simple_post_process_owned_with_bounds(
        config,
        region,
        aligned_variants,
        splice,
        config
            .compat
            .strict_region_bounds
            .then_some(single_bounds.as_slice()),
    )
}

/// Owned fast-path counterpart of [`simple_post_process`] using pre-merged chromosome
/// target intervals when `strict_bounds` is supplied (`#369`, `AR-4`).
pub(crate) fn simple_post_process_owned_with_bounds(
    config: &Config,
    region: &Region,
    mut aligned_variants: crate::javamap::JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
    strict_bounds: Option<&[(i32, i32)]>,
) -> Vec<Variant> {
    if config.compat.coordinate_sorted_vars {
        let mut keys: Vec<i32> = aligned_variants.keys().copied().collect();
        keys.sort_unstable();
        simple_post_process_core(
            config,
            region,
            keys.into_iter().filter_map(|pos| {
                aligned_variants
                    .get_mut(&pos)
                    .map(|vars| (pos, std::borrow::Cow::Owned(std::mem::take(vars))))
            }),
            splice,
            strict_bounds,
        )
    } else {
        simple_post_process_core(
            config,
            region,
            aligned_variants
                .into_entries()
                .map(|(pos, vars)| (pos, std::borrow::Cow::Owned(vars))),
            splice,
            strict_bounds,
        )
    }
}

fn simple_post_process_core<'a>(
    config: &Config,
    region: &Region,
    entries: impl IntoIterator<Item = (i32, std::borrow::Cow<'a, Vars>)>,
    splice: &rustc_hash::FxHashSet<String>,
    strict_bounds: Option<&[(i32, i32)]>,
) -> Vec<Variant> {
    let mut vrefs = Vec::new();

    let reset_crispr = config.crispr_cutting_site.unwrap_or(0) == 0;
    let refresh_vt = config.compat.recompute_vartype_after_adj_complex;
    let skip_eq_non_pileup = !config.pileup && config.compat.sanitize_iupac_alt_alleles;
    let preserve_shifted_ref = config.pileup && config.compat.pileup_preserve_shifted_ref;

    // Iterate entries in JavaHashMap table-bucket iteration order under default parity
    // (ports SimplePostProcessModule.java:35; see DIVERGENCES.md UB-3), or ascending
    // genomic coordinate order when `CompatProfile::coordinate_sorted_vars` is active.
    for (position, variants_on_position) in entries {
        if variants_on_position.sv.is_empty()
            && (position < region.start.0 as i32 || position > region.end.0 as i32)
        {
            continue;
        }

        let (mut ref_var_cow, vvar_iter, sv_cow) = super::unpack_cow_vars(variants_on_position);
        let vvar_len = vvar_iter.len();

        let pos_start = vrefs.len();
        if vvar_len == 0 {
            if !config.pileup {
                continue;
            }
            if let Some(vref_cow) = ref_var_cow.take() {
                let mut vref = vref_cow.into_owned();
                vref.var_type = String::new();
                if reset_crispr {
                    vref.crispr = 0;
                }
                vrefs.push(vref);
            } else if !config.compat.suppress_pileup_sentinels {
                // ports SimplePostProcessModule.java:55-58: under -p an entry
                // with no variants and no reference variant still prints, via
                // `SimpleOutputVariant(null, region, sv, position)` whose fields
                // keep the printer defaults (genotype "", bias "0;0", flank seqs
                // ""). CLAUDE.md §5: live-oracle output is replicated as-is and
                // tracked upstream, never silently fixed.
                vrefs.push(make_null_record(position));
            }
        } else if preserve_shifted_ref {
            // AstraZeneca-NGS/VarDictJava#330: run per-variant `adj_complex_compat` + `crispr`
            // reset inside the per-position loop BEFORE checking whether >= 1 variant survived
            // at `position` with zero surviving variants keeping `start_pos == position`.
            let had_null_ref = ref_var_cow.is_none();
            let mut single_pre_adj_shifted = false;
            for vref_cow in vvar_iter {
                if vref_cow.ref_allele.contains('N') {
                    continue;
                }
                let vtype_cow = vref_cow.var_type_cow();
                let mut vref = vref_cow.into_owned();
                if vvar_len == 1 && vref.start_pos != position {
                    single_pre_adj_shifted = true;
                }
                if vref.var_type != vtype_cow.as_ref() {
                    vref.var_type.clear();
                    vref.var_type.push_str(&vtype_cow);
                }
                if vref.var_type == "Complex" {
                    vref.adj_complex_compat(refresh_vt);
                }
                if reset_crispr {
                    vref.crispr = 0;
                }
                vrefs.push(vref);
            }
            let has_survived = vrefs.len() > pos_start;
            let any_at_pos = vrefs
                .get(pos_start..)
                .is_some_and(|s| s.iter().any(|v| v.start_pos == position));
            if has_survived && !any_at_pos {
                if let Some(rv_cow) = ref_var_cow.take().filter(|rv| rv.start_pos == position) {
                    let mut ref_var = rv_cow.into_owned();
                    ref_var.var_type = String::new();
                    if reset_crispr {
                        ref_var.crispr = 0;
                    }
                    vrefs.push(ref_var);
                    if let Some(s) = vrefs.get_mut(pos_start..) {
                        if config.compat.coordinate_sorted_vars {
                            s.sort_by_key(|v| v.start_pos);
                        } else if !s.is_empty() {
                            s.rotate_right(1);
                        }
                    }
                } else if vvar_len == 1
                    && !config.compat.suppress_pileup_sentinels
                    && had_null_ref
                    && single_pre_adj_shifted
                {
                    let null_rec = make_null_record(position);
                    let mut zero_ref = Variant::new();
                    zero_ref.var_type = String::new();
                    vrefs.push(null_rec);
                    vrefs.push(zero_ref);
                    if let Some(s) = vrefs.get_mut(pos_start..) {
                        if s.len() >= 2 {
                            s.rotate_right(2);
                        }
                    }
                }
            }
        } else {
            for vref_cow in vvar_iter {
                if vref_cow.ref_allele.contains('N') {
                    continue;
                }
                if vref_cow.ref_allele == vref_cow.var_allele && !config.pileup {
                    continue;
                }
                let vtype_cow = vref_cow.var_type_cow();
                let is_good =
                    vref_cow.is_good_var(ref_var_cow.as_deref(), &vtype_cow, splice, config);
                if !is_good && !config.pileup {
                    continue;
                }
                let mut vref = vref_cow.into_owned();
                if vref.start_pos != position && config.pileup && vvar_len == 1 {
                    if let Some(rv_cow) = ref_var_cow.take() {
                        let mut ref_var = rv_cow.into_owned();
                        ref_var.var_type = String::new();
                        if reset_crispr {
                            ref_var.crispr = 0;
                        }
                        vrefs.push(ref_var);
                    } else if !config.compat.suppress_pileup_sentinels {
                        // ports SimplePostProcessModule.java:75-79: Java prints a
                        // SimpleOutputVariant wrapping null immediately (the
                        // null-record row at the map key), then continues with a
                        // default `new Variant()` whose startPosition/endPosition
                        // stay 0 — that 0 leaks into .var columns 4-5 (the
                        // upstream pos-0 quirk; replicated per §5, see
                        // DIVERGENCES.md "Known non-bugs").
                        vrefs.push(make_null_record(position));
                        let mut ref_var = Variant::new();
                        ref_var.var_type = String::new();
                        vrefs.push(ref_var);
                    }
                }
                if vref.var_type != vtype_cow.as_ref() {
                    vref.var_type.clear();
                    vref.var_type.push_str(&vtype_cow);
                }
                if vref.var_type == "Complex" {
                    vref.adj_complex_compat(refresh_vt);
                    if skip_eq_non_pileup && vref.ref_allele == vref.var_allele {
                        continue;
                    }
                }
                if reset_crispr {
                    vref.crispr = 0;
                }
                vrefs.push(vref);
            }
        }

        if let Some((last, prior)) = vrefs.get_mut(pos_start..).and_then(|s| s.split_last_mut()) {
            for r in prior {
                r.sv.clear();
                r.sv.push_str(&sv_cow);
            }
            last.sv = sv_cow.into_owned();
        }
    }

    if let Some(merged_bounds) = strict_bounds {
        // AstraZeneca-NGS/VarDictJava#369 + AstraZeneca-NGS/VarDictJava#398:
        // retain 1-bp-upstream-padded symbolic SVs anchored at `region.start - 1`.
        let allow_padded_sv = config.compat.canonical_sv_vcf_anchoring;
        vrefs.retain(|v| crate::variant_within_strict_bounds(merged_bounds, v, allow_padded_sv));
    } else if config.compat.strict_region_bounds {
        vrefs.retain(|v| v.start_pos >= 1);
    }

    // VCF mode: ports var2vcf_valid.pl:41,129-137 — bucket rows by the RAW .var
    // start column ($hash{$chr}->{$a[3]}, column 3), walk buckets in ascending
    // numeric position (:129), sort each bucket by column 14 (af) descending
    // (:135), and emit only the top row ($ALL = 1 without -A, :137). The stable
    // sort keeps af ties in .var stream order, like Perl's default mergesort.
    // Bucketing must use the raw start (what simple_36 writes), NOT the
    // ToVarsBuilder position key: indels and SNVs can share one key while
    // having different raw starts, and different keys can share a raw start.
    if !config.emit_intermediate && !config.vcf_all_alleles {
        let exact_tie = config.compat.exact_binary_half_even;
        let mut buckets: BTreeMap<i32, Variant> = BTreeMap::new();
        for v in vrefs {
            let v_af = crate::types::round_half_even_compat(v.frequency, 4, exact_tie);
            match buckets.entry(v.start_pos) {
                std::collections::btree_map::Entry::Vacant(e) => {
                    e.insert(v);
                }
                std::collections::btree_map::Entry::Occupied(mut e) => {
                    let cur_af =
                        crate::types::round_half_even_compat(e.get().frequency, 4, exact_tie);
                    if v_af > cur_af {
                        e.insert(v);
                    }
                }
            }
        }
        vrefs = buckets.into_values().collect();
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

        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        let splice = rustc_hash::FxHashSet::default();

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
        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        aligned_variants.insert(20159221, g1);
        aligned_variants.insert(20159222, g2);

        let splice = rustc_hash::FxHashSet::default();
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
        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        aligned_variants.insert(20179028, vars);

        let splice = rustc_hash::FxHashSet::default();
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
        let mut aligned_rev = crate::javamap::JavaHashMap::new();
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
        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        aligned_variants.insert(6447, vars);
        let splice = rustc_hash::FxHashSet::default();

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
        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        aligned_variants.insert(6447, Vars::new());
        let splice = rustc_hash::FxHashSet::default();

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
        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        aligned_variants.insert(6447, vars);
        aligned_variants.insert(6500, Vars::new());
        let splice = rustc_hash::FxHashSet::default();

        let out = simple_post_process(&config, &region, &aligned_variants, &splice);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].start_pos, 6449);
        assert!(!out[0].null_record);
    }

    #[test]
    fn test_vcf_mode_single_winner_tie_keeps_first_and_higher_af_wins() {
        let mut config = Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let splice = rustc_hash::FxHashSet::default();

        // 1. Two variants at start_pos = 150 with equal af_4dp = 0.2500 keep the first-seen allele.
        let mut vars_tie = Vars::new();
        vars_tie.variants.push(mk_variant(150, 150, "A", "G", 0.25));
        vars_tie.variants.push(mk_variant(150, 150, "A", "T", 0.25));
        let mut aligned_tie = crate::javamap::JavaHashMap::new();
        aligned_tie.insert(150, vars_tie);

        let kept_tie = simple_post_process(&config, &region, &aligned_tie, &splice);
        assert_eq!(kept_tie.len(), 1);
        assert_eq!(kept_tie[0].start_pos, 150);
        assert_eq!(kept_tie[0].var_allele, "G");

        // 2. A second variant at start_pos = 150 with higher af_4dp = 0.5000 (an unsplit <INV>
        // with sv = "0-100-1") replaces the lower-AF 0.2500 variant.
        let mut g_low = Vars::new();
        g_low.variants.push(mk_variant(150, 150, "A", "G", 0.25));
        let mut g_inv = Vars::new();
        g_inv.sv = "0-100-1".to_string();
        g_inv
            .variants
            .push(mk_variant(150, 180, "ACGT", "<INV>", 0.50));
        let mut aligned_inv = crate::javamap::JavaHashMap::new();
        aligned_inv.insert(150, g_low);
        aligned_inv.insert(151, g_inv);

        let kept_inv = simple_post_process(&config, &region, &aligned_inv, &splice);
        assert_eq!(kept_inv.len(), 1);
        assert_eq!(kept_inv[0].start_pos, 150);
        assert_eq!(kept_inv[0].var_allele, "<INV>");
        assert_eq!(kept_inv[0].sv, "0-100-1");

        // 3. Direct library callers passing unrounded `Variant.frequency` do not trigger debug assertions.
        let mut g_unrounded = Vars::new();
        g_unrounded
            .variants
            .push(mk_variant(160, 160, "C", "T", 0.3333333));
        let mut aligned_unrounded = crate::javamap::JavaHashMap::new();
        aligned_unrounded.insert(160, g_unrounded);
        let kept_unrounded = simple_post_process(&config, &region, &aligned_unrounded, &splice);
        assert_eq!(kept_unrounded.len(), 1);
        assert_eq!(kept_unrounded[0].var_allele, "T");

        // 4. Exhaustive field-for-field parity between `simple_post_process_owned` and
        // `simple_post_process` across germline, -p pileup (empty with/without ref_variant,
        // shifted-singleton with/without ref_variant, multi-variant with sv != ""), Complex,
        // CRISPR, and multi-allelic VCF/intermediate buckets.
        fn assert_variant_fields_eq(a: &Variant, b: &Variant) {
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

        let mut parity_map = crate::javamap::JavaHashMap::new();
        // Out-of-region without SV (skipped) vs with SV (kept)
        let mut out_no_sv = Vars::new();
        out_no_sv.variants.push(mk_variant(90, 90, "A", "G", 0.5));
        parity_map.insert(90, out_no_sv);
        let mut out_with_sv = Vars::new();
        out_with_sv.sv = "95-210-4".to_string();
        let mut v_sv_out = mk_variant(95, 95, "A", "<DEL>", 0.4);
        v_sv_out.crispr = 7;
        out_with_sv.variants.push(v_sv_out);
        parity_map.insert(95, out_with_sv);

        // Pos 105: empty variants without ref_variant (null_record under -p)
        let mut empty_no_ref = Vars::new();
        empty_no_ref.sv = "105-180-2".to_string();
        parity_map.insert(105, empty_no_ref);

        // Pos 110: empty variants with ref_variant
        let mut empty_with_ref = Vars::new();
        let mut ref_110 = mk_variant(110, 110, "C", "C", 1.0);
        ref_110.var_type = "SNV".to_string();
        ref_110.crispr = 3;
        empty_with_ref.ref_variant = Some(ref_110);
        empty_with_ref.sv = "110-190-1".to_string();
        parity_map.insert(110, empty_with_ref);

        // Pos 120: shifted-singleton WITHOUT ref_variant + non-empty sv
        let mut shifted_no_ref = Vars::new();
        shifted_no_ref.sv = "120-160-3".to_string();
        shifted_no_ref
            .variants
            .push(mk_variant(122, 122, "G", "GTT", 0.6));
        parity_map.insert(120, shifted_no_ref);

        // Pos 130: shifted-singleton WITH ref_variant + non-empty sv
        let mut shifted_with_ref = Vars::new();
        shifted_with_ref.sv = "130-170-5".to_string();
        let mut ref_130 = mk_variant(130, 130, "T", "T", 0.4);
        ref_130.var_type = "SNV".to_string();
        shifted_with_ref.ref_variant = Some(ref_130);
        shifted_with_ref
            .variants
            .push(mk_variant(132, 134, "TAC", "T", 0.55));
        parity_map.insert(130, shifted_with_ref);

        // Pos 140: vvar.len() > 1 with sv != "", N-ref skip, Complex adj_complex, and
        // shifted start_pos (must NOT trigger vvar.len() == 1 shifted-singleton branch even
        // when one sibling is skipped for 'N' in ref_allele!)
        let mut multi_sv = Vars::new();
        multi_sv.sv = "140-199-8".to_string();
        multi_sv.ref_variant = Some(mk_variant(140, 140, "A", "A", 0.2));
        multi_sv
            .variants
            .push(mk_variant(140, 140, "AN", "AT", 0.5));
        let mut v_complex = mk_variant(141, 144, "ACGT", "AGGT", 0.45);
        v_complex.left_seq = "TT".to_string();
        v_complex.right_seq = "CC".to_string();
        v_complex.crispr = 9;
        multi_sv.variants.push(v_complex);
        multi_sv.variants.push(mk_variant(142, 142, "C", "A", 0.35));
        parity_map.insert(140, multi_sv);

        // Pos 150 & 151: multi-allelic VCF/intermediate bucket tie + higher-AF override + low-AF filter
        let mut bucket_a = Vars::new();
        bucket_a.variants.push(mk_variant(150, 150, "A", "G", 0.25));
        bucket_a.variants.push(mk_variant(150, 150, "A", "T", 0.25));
        bucket_a
            .variants
            .push(mk_variant(150, 150, "A", "C", 0.001));
        bucket_a.variants.push(mk_variant(150, 150, "A", "A", 0.50));
        parity_map.insert(150, bucket_a);
        let mut bucket_b = Vars::new();
        bucket_b.sv = "0-100-1".to_string();
        bucket_b
            .variants
            .push(mk_variant(150, 180, "ACGT", "<INV>", 0.50));
        parity_map.insert(151, bucket_b);

        for pileup in [false, true] {
            for emit_intermediate in [false, true] {
                for vcf_all_alleles in [false, true] {
                    for compat in [
                        crate::config::COMPAT_VARDICTJAVA_1_8_3,
                        crate::config::COMPAT_OFF,
                    ] {
                        for crispr_site in [None, Some(150)] {
                            let mut cfg = config.clone();
                            cfg.pileup = pileup;
                            cfg.emit_intermediate = emit_intermediate;
                            cfg.vcf_all_alleles = vcf_all_alleles;
                            cfg.compat = compat;
                            cfg.crispr_cutting_site = crispr_site;
                            let borrowed = simple_post_process(&cfg, &region, &parity_map, &splice);
                            let owned = simple_post_process_owned(
                                &cfg,
                                &region,
                                parity_map.clone(),
                                &splice,
                            );
                            assert_eq!(owned.len(), borrowed.len());
                            for (got, want) in owned.iter().zip(borrowed.iter()) {
                                assert_variant_fields_eq(got, want);
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_simple_vcf_all_alleles_and_compat_off_sentinels_and_bounds() {
        let mut config = Config::default();
        config.freq = 0.01;
        config.min_reads = 2;
        config.read_pos = 5.0;
        config.base_qual = 22.5;
        config.qratio = 1.5;

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let splice = rustc_hash::FxHashSet::default();

        // 1. vcf_all_alleles = true preserves all passing alleles at the same start_pos in VCF mode.
        let mut g1 = Vars::new();
        g1.variants.push(mk_variant(150, 150, "A", "AT", 0.45));
        g1.variants.push(mk_variant(150, 150, "A", "G", 0.20));
        let mut g2 = Vars::new();
        g2.variants.push(mk_variant(150, 151, "AT", "A", 0.20));
        let mut multi_map = crate::javamap::JavaHashMap::new();
        multi_map.insert(149, g1);
        multi_map.insert(150, g2);

        let single_winner = simple_post_process(&config, &region, &multi_map, &splice);
        assert_eq!(single_winner.len(), 1);
        assert_eq!(single_winner[0].var_allele, "AT");

        let mut all_cfg = config.clone();
        all_cfg.vcf_all_alleles = true;
        let all_alleles = simple_post_process(&all_cfg, &region, &multi_map, &splice);
        assert_eq!(all_alleles.len(), 3);

        // 2. suppress_pileup_sentinels (UB-1) under COMPAT_OFF:
        //    - empty Vars without ref_variant emits 0 rows (no null_record)
        //    - shifted singleton without ref_variant emits ONLY the real shifted variant (no null_record, no POS=0)
        //    - empty Vars WITH ref_variant still emits the real reference row
        let mut p_map = crate::javamap::JavaHashMap::new();
        p_map.insert(110, Vars::new());
        let mut with_ref = Vars::new();
        with_ref.ref_variant = Some(mk_variant(120, 120, "C", "C", 1.0));
        p_map.insert(120, with_ref);
        let mut shifted = Vars::new();
        shifted.variants.push(mk_variant(132, 132, "G", "GTT", 0.8));
        p_map.insert(130, shifted);

        let mut p_compat_off = pileup_config();
        p_compat_off.compat = crate::config::COMPAT_OFF;
        let p_out = simple_post_process(&p_compat_off, &region, &p_map, &splice);
        assert_eq!(p_out.len(), 2);
        assert!(p_out.iter().all(|v| !v.null_record && v.start_pos >= 1));
        let p_starts: Vec<i32> = p_out.iter().map(|v| v.start_pos).collect();
        assert!(p_starts.contains(&120) && p_starts.contains(&132));

        // 3. strict_region_bounds (#369) under COMPAT_OFF:
        //    variants shifted outside [100, 200] (e.g. left-anchored deletion at 99, right-shifted indel at 203,
        //    out-of-bounds SV at 95, or POS <= 0) are dropped while 100 and 200 are retained.
        let mut b_map = crate::javamap::JavaHashMap::new();
        let mut v_left_del = Vars::new();
        v_left_del
            .variants
            .push(mk_variant(99, 100, "CA", "C", 0.50));
        b_map.insert(100, v_left_del);

        let mut v_in_lo = Vars::new();
        v_in_lo.variants.push(mk_variant(100, 100, "A", "G", 0.50));
        b_map.insert(101, v_in_lo);

        let mut v_in_hi = Vars::new();
        v_in_hi.variants.push(mk_variant(200, 200, "T", "C", 0.50));
        b_map.insert(199, v_in_hi);

        let mut v_right_shift = Vars::new();
        v_right_shift
            .variants
            .push(mk_variant(203, 203, "G", "GAA", 0.50));
        b_map.insert(200, v_right_shift);

        let mut v_sv_out = Vars::new();
        v_sv_out.sv = "95-210-4".to_string();
        v_sv_out
            .variants
            .push(mk_variant(95, 210, "A", "<DEL>", 0.50));
        b_map.insert(95, v_sv_out);

        let default_bounds = simple_post_process(&config, &region, &b_map, &splice);
        assert_eq!(default_bounds.len(), 5);

        let mut strict_cfg = config.clone();
        strict_cfg.compat = crate::config::COMPAT_OFF;
        let strict_out = simple_post_process(&strict_cfg, &region, &b_map, &splice);
        let strict_starts: Vec<i32> = strict_out.iter().map(|v| v.start_pos).collect();
        assert_eq!(strict_starts, vec![100, 200]);
        let strict_owned = simple_post_process_owned(&strict_cfg, &region, b_map, &splice);
        assert_eq!(strict_owned.len(), 2);

        // 4. exact_binary_half_even (5e-5 AF tie) in VCF top-1-per-start_pos selection:
        let tie_4dp = f64::from_bits(0x3f0a_36e2_eb1c_432d); // 5e-5
        let mut tie_cfg = config.clone();
        tie_cfg.freq = 0.0;
        let mut tie_vars = Vars::new();
        tie_vars.variants.push(mk_variant(150, 150, "A", "G", 0.0));
        tie_vars
            .variants
            .push(mk_variant(150, 150, "A", "T", tie_4dp));
        let mut tie_map = crate::javamap::JavaHashMap::new();
        tie_map.insert(150, tie_vars);

        let vj_winner = simple_post_process(&tie_cfg, &region, &tie_map, &splice);
        assert_eq!(vj_winner.len(), 1);
        assert_eq!(vj_winner[0].var_allele, "G");

        let mut off_tie_cfg = tie_cfg.clone();
        off_tie_cfg.compat = crate::config::COMPAT_OFF;
        let off_winner = simple_post_process(&off_tie_cfg, &region, &tie_map, &splice);
        assert_eq!(off_winner.len(), 1);
        assert_eq!(off_winner[0].var_allele, "T");
    }

    #[test]
    fn compat_off_coordinate_sorted_vars_orders_var_rows_and_vcf_ties() {
        let mut base_cfg = Config::default();
        base_cfg.freq = 0.01;
        base_cfg.min_reads = 2;
        base_cfg.read_pos = 5.0;
        base_cfg.base_qual = 22.5;
        base_cfg.qratio = 1.5;

        let region = Region {
            chr: "20".to_string(),
            start: Pos1(20_202_800),
            end: Pos1(20_203_000),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let splice = rustc_hash::FxHashSet::default();
        let strict_bounds = single_region_bounds_i32(&region);

        // OpenJDK HashMap spread `(h ^ (h >>> 16)) & 15` at default capacity 16:
        // - 20_202_900 (0x0134_4594) -> bucket 0
        // - 20_202_850 (0x0134_4562) -> bucket 6
        // - 20_202_876 (0x0134_457c, k_hi) -> bucket 8
        // - 20_202_910 (0x0134_459e) -> bucket 10
        // - 20_202_875 (0x0134_457b, k_lo) -> bucket 15
        let k_lo = 20_202_875;
        let k_hi = k_lo + 1;

        let mut map = crate::javamap::JavaHashMap::new();
        let mut vars_lo = Vars::new();
        vars_lo
            .variants
            .push(mk_variant(k_lo, k_lo, "T", "A", 0.05));
        map.insert(k_lo, vars_lo);

        let mut vars_hi = Vars::new();
        vars_hi
            .variants
            .push(mk_variant(k_lo, k_hi, "TA", "T", 0.05));
        map.insert(k_hi, vars_hi);

        for (pos, r, a, af) in [
            (20_202_900, "C", "G", 0.20),
            (20_202_850, "G", "T", 0.15),
            (20_202_910, "A", "C", 0.30),
        ] {
            let mut vars = Vars::new();
            vars.variants.push(mk_variant(pos, pos, r, a, af));
            map.insert(pos, vars);
        }

        let bucket_keys: Vec<i32> = map.keys().copied().collect();
        assert_eq!(
            bucket_keys,
            vec![20_202_900, 20_202_850, k_hi, 20_202_910, k_lo]
        );
        assert!(map.contains_both_and_first_precedes(&k_hi, &k_lo));
        assert!(!bucket_keys.windows(2).all(|w| w[0] <= w[1]));

        let summarize = |rows: &[Variant]| -> Vec<(i32, String, String)> {
            rows.iter()
                .map(|v| (v.start_pos, v.ref_allele.clone(), v.var_allele.clone()))
                .collect()
        };

        // 1. Default parity (`COMPAT_VARDICTJAVA_1_8_3`): JavaHashMap bucket order.
        let mut vj_var_cfg = base_cfg.clone();
        vj_var_cfg.compat = crate::config::COMPAT_VARDICTJAVA_1_8_3;
        vj_var_cfg.emit_intermediate = true;
        let expected_vj_var = vec![
            (20_202_900, "C".to_string(), "G".to_string()),
            (20_202_850, "G".to_string(), "T".to_string()),
            (k_lo, "TA".to_string(), "T".to_string()),
            (20_202_910, "A".to_string(), "C".to_string()),
            (k_lo, "T".to_string(), "A".to_string()),
        ];
        assert_eq!(
            summarize(&simple_post_process(&vj_var_cfg, &region, &map, &splice)),
            expected_vj_var
        );
        assert_eq!(
            summarize(&simple_post_process_owned_with_bounds(
                &vj_var_cfg,
                &region,
                map.clone(),
                &splice,
                None,
            )),
            expected_vj_var
        );

        let mut vj_vcf_cfg = base_cfg.clone();
        vj_vcf_cfg.compat = crate::config::COMPAT_VARDICTJAVA_1_8_3;
        vj_vcf_cfg.emit_intermediate = false;
        let expected_vj_vcf = vec![
            (20_202_850, "G".to_string(), "T".to_string()),
            (k_lo, "TA".to_string(), "T".to_string()),
            (20_202_900, "C".to_string(), "G".to_string()),
            (20_202_910, "A".to_string(), "C".to_string()),
        ];
        assert_eq!(
            summarize(&simple_post_process(&vj_vcf_cfg, &region, &map, &splice)),
            expected_vj_vcf
        );
        assert_eq!(
            summarize(&simple_post_process_owned_with_bounds(
                &vj_vcf_cfg,
                &region,
                map.clone(),
                &splice,
                None,
            )),
            expected_vj_vcf
        );

        // 2. Opt-in `coordinate_sorted_vars = true` (isolated switch and full `COMPAT_OFF`):
        //    ascending map-key order in `.var` and deterministic lowest-key (`k_lo` SNV) winner in VCF.
        let mut isolated_compat = crate::config::COMPAT_VARDICTJAVA_1_8_3;
        isolated_compat.coordinate_sorted_vars = true;
        let expected_sorted_var = vec![
            (20_202_850, "G".to_string(), "T".to_string()),
            (k_lo, "T".to_string(), "A".to_string()),
            (k_lo, "TA".to_string(), "T".to_string()),
            (20_202_900, "C".to_string(), "G".to_string()),
            (20_202_910, "A".to_string(), "C".to_string()),
        ];
        let expected_sorted_vcf = vec![
            (20_202_850, "G".to_string(), "T".to_string()),
            (k_lo, "T".to_string(), "A".to_string()),
            (20_202_900, "C".to_string(), "G".to_string()),
            (20_202_910, "A".to_string(), "C".to_string()),
        ];

        for compat in [isolated_compat, crate::config::COMPAT_OFF] {
            let bounds = compat
                .strict_region_bounds
                .then_some(strict_bounds.as_slice());

            let mut sorted_var_cfg = base_cfg.clone();
            sorted_var_cfg.compat = compat;
            sorted_var_cfg.emit_intermediate = true;
            let got_var_borrowed = simple_post_process(&sorted_var_cfg, &region, &map, &splice);
            let got_var_owned = simple_post_process_owned_with_bounds(
                &sorted_var_cfg,
                &region,
                map.clone(),
                &splice,
                bounds,
            );
            assert!(got_var_borrowed
                .windows(2)
                .all(|w| w[0].start_pos <= w[1].start_pos));
            assert_eq!(summarize(&got_var_borrowed), expected_sorted_var);
            assert_eq!(summarize(&got_var_owned), expected_sorted_var);

            let mut sorted_vcf_cfg = base_cfg.clone();
            sorted_vcf_cfg.compat = compat;
            sorted_vcf_cfg.emit_intermediate = false;
            let got_vcf_borrowed = simple_post_process(&sorted_vcf_cfg, &region, &map, &splice);
            let got_vcf_owned = simple_post_process_owned_with_bounds(
                &sorted_vcf_cfg,
                &region,
                map.clone(),
                &splice,
                bounds,
            );
            assert_eq!(summarize(&got_vcf_borrowed), expected_sorted_vcf);
            assert_eq!(summarize(&got_vcf_owned), expected_sorted_vcf);
        }
    }

    #[test]
    fn test_strict_bounds_retains_1bp_upstream_padded_symbolic_sv_at_region_start() {
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "SV_BOUND".to_string(),
            amplicon: None,
        };
        let strict_bounds = [(100i32, 200i32)];
        let splice = rustc_hash::FxHashSet::default();

        let mut map = crate::javamap::JavaHashMap::new();
        // 1. 1-bp-upstream-padded symbolic <DEL> at start_pos = 99 (99 + 1 == 100 in bounds)
        let mut g_del99 = Vars::new();
        g_del99.sv = "5-2-1".to_string();
        g_del99
            .variants
            .push(mk_variant(99, 180, "A", "<DEL>", 0.45));
        map.insert(99, g_del99);
        // 2. Symbolic <DUP> at start_pos = 98 (98 + 1 == 99 still outside [100, 200] -> dropped)
        let mut g_dup98 = Vars::new();
        g_dup98.sv = "5-2-1".to_string();
        g_dup98
            .variants
            .push(mk_variant(98, 180, "A", "<DUP>", 0.40));
        map.insert(98, g_dup98);
        // 3. Non-symbolic deletion at start_pos = 99 -> dropped even when canonical_sv_vcf_anchoring is true
        let mut g_indel99 = Vars::new();
        g_indel99.sv = "1-0-1".to_string();
        g_indel99
            .variants
            .push(mk_variant(99, 101, "CAT", "C", 0.35));
        map.insert(100, g_indel99);
        // 4. Degenerate symbolic SV at start_pos = 0 (start_pos < 1 -> dropped) and end_pos == start_pos (dropped)
        let mut g_zero = Vars::new();
        g_zero.sv = "3-1-1".to_string();
        g_zero.variants.push(mk_variant(0, 50, "A", "<INV>", 0.30));
        map.insert(1, g_zero);
        let mut g_eq = Vars::new();
        g_eq.sv = "3-1-1".to_string();
        g_eq.variants.push(mk_variant(299, 299, "A", "<DEL>", 0.30));
        map.insert(299, g_eq);

        let multi_bounds = [(100i32, 200i32), (300i32, 400i32)];

        // Without `canonical_sv_vcf_anchoring`, start_pos = 99 is dropped under `strict_region_bounds`.
        let mut cfg_strict_only = Config {
            emit_intermediate: true,
            ..Config::default()
        };
        cfg_strict_only.compat.strict_region_bounds = true;
        cfg_strict_only.compat.canonical_sv_vcf_anchoring = false;
        let out_strict_only = simple_post_process_owned_with_bounds(
            &cfg_strict_only,
            &region,
            map.clone(),
            &splice,
            Some(&multi_bounds),
        );
        assert!(out_strict_only.is_empty());

        // Under COMPAT_OFF (`strict_region_bounds && canonical_sv_vcf_anchoring`), ONLY the valid
        // 1-bp-upstream-padded `<DEL>` at start_pos = 99 (99 + 1 == 100) is retained!
        let cfg_off = Config {
            emit_intermediate: true,
            compat: crate::config::COMPAT_OFF,
            ..Config::default()
        };
        let out_owned = simple_post_process_owned_with_bounds(
            &cfg_off,
            &region,
            map.clone(),
            &splice,
            Some(&multi_bounds),
        );
        assert_eq!(out_owned.len(), 1);
        assert_eq!(out_owned[0].start_pos, 99);
        assert_eq!(out_owned[0].var_allele, "<DEL>");

        let out_borrowed = simple_post_process(&cfg_off, &region, &map, &splice);
        let _ = strict_bounds;
        assert_eq!(out_borrowed.len(), 1);
        assert_eq!(out_borrowed[0].start_pos, 99);
        assert_eq!(out_borrowed[0].var_allele, "<DEL>");
    }

    #[test]
    fn test_wave5_pileup_preserve_shifted_ref_and_adj_complex_vartype_and_iupac_eq_skip() {
        let splice = rustc_hash::FxHashSet::default();
        let region_tp53 = Region {
            chr: "17".to_string(),
            start: Pos1(7_578_700),
            end: Pos1(7_578_730),
            gene: "TP53".to_string(),
            amplicon: None,
        };

        // 1. AstraZeneca-NGS/VarDictJava#330: `-p` at `17:7578712` with 2 shifted deletions
        //    (`7578711 CA>C` and `7578711 CAA>C`) + hydrated `ref_variant` at `7578712`:
        let mut p_del2_map = crate::javamap::JavaHashMap::new();
        let mut vars_712 = Vars::new();
        vars_712.sv = "7578700-7578900-2".to_string();
        let mut ref_712 = mk_variant(7_578_712, 7_578_712, "A", "A", 0.0);
        ref_712.var_type = "SNV".to_string();
        ref_712.crispr = 9;
        vars_712.ref_variant = Some(ref_712);
        vars_712
            .variants
            .push(mk_variant(7_578_711, 7_578_712, "CA", "C", 0.55));
        vars_712
            .variants
            .push(mk_variant(7_578_711, 7_578_713, "CAA", "C", 0.45));
        p_del2_map.insert(7_578_712, vars_712);

        // Default `COMPAT_VARDICTJAVA_1_8_3`: `vvar_len == 2 != 1` drops `ref_variant` at 7578712.
        let vj_p_cfg = pileup_config();
        let vj_del2 = simple_post_process(&vj_p_cfg, &region_tp53, &p_del2_map, &splice);
        assert_eq!(vj_del2.len(), 2);
        assert!(vj_del2.iter().all(|v| v.start_pos == 7_578_711));

        // Isolated `pileup_preserve_shifted_ref` + `coordinate_sorted_vars` (and `COMPAT_OFF`):
        // emits `[del1@7578711, del2@7578711, ref@7578712]` in monotonic order with `sv` stamped and `crispr == 0`.
        let mut iso_p_cfg = pileup_config();
        iso_p_cfg.compat.pileup_preserve_shifted_ref = true;
        iso_p_cfg.compat.coordinate_sorted_vars = true;
        for cfg in [&iso_p_cfg, &{
            let mut c = pileup_config();
            c.compat = crate::config::COMPAT_OFF;
            c
        }] {
            let got_b = simple_post_process(cfg, &region_tp53, &p_del2_map, &splice);
            let got_o = simple_post_process_owned(cfg, &region_tp53, p_del2_map.clone(), &splice);
            assert_eq!(got_b.len(), 3);
            assert_eq!(got_o.len(), 3);
            let starts: Vec<i32> = got_b.iter().map(|v| v.start_pos).collect();
            assert_eq!(starts, vec![7_578_711, 7_578_711, 7_578_712]);
            assert_eq!(got_b[0].ref_allele, "CA");
            assert_eq!(got_b[1].ref_allele, "CAA");
            assert_eq!(got_b[2].var_type, "");
            assert_eq!(got_b[2].crispr, 0);
            assert!(got_b.iter().all(|v| v.sv == "7578700-7578900-2"));
            let starts_o: Vec<i32> = got_o.iter().map(|v| v.start_pos).collect();
            assert_eq!(starts_o, starts);
        }

        // Mixed `[right_shifted > P, left_shifted < P]` under `coordinate_sorted_vars` orders
        // `[left_shifted < P, ref@P, right_shifted > P]` monotonically.
        let mut mixed_map = crate::javamap::JavaHashMap::new();
        let mut mixed_vars = Vars::new();
        mixed_vars.ref_variant = Some(mk_variant(7_578_712, 7_578_712, "A", "A", 0.2));
        mixed_vars
            .variants
            .push(mk_variant(7_578_714, 7_578_714, "G", "GTT", 0.5));
        mixed_vars
            .variants
            .push(mk_variant(7_578_711, 7_578_712, "CA", "C", 0.3));
        mixed_map.insert(7_578_712, mixed_vars);
        let mixed_out = simple_post_process(&iso_p_cfg, &region_tp53, &mixed_map, &splice);
        let mixed_starts: Vec<i32> = mixed_out.iter().map(|v| v.start_pos).collect();
        assert_eq!(mixed_starts, vec![7_578_711, 7_578_712, 7_578_714]);

        // Unhydrated `ref_variant = None` or `start_pos == 0` when `vvar_len >= 2` NEVER emits pos0/null sentinels.
        let mut unhyd_map = crate::javamap::JavaHashMap::new();
        let mut unhyd_vars = Vars::new();
        unhyd_vars.ref_variant = Some(Variant::new()); // unhydrated start_pos == 0
        unhyd_vars
            .variants
            .push(mk_variant(7_578_711, 7_578_712, "CA", "C", 0.5));
        unhyd_vars
            .variants
            .push(mk_variant(7_578_711, 7_578_713, "CAA", "C", 0.4));
        unhyd_map.insert(7_578_712, unhyd_vars);
        let mut iso_no_sup = pileup_config();
        iso_no_sup.compat.pileup_preserve_shifted_ref = true;
        iso_no_sup.compat.suppress_pileup_sentinels = false;
        let unhyd_out = simple_post_process(&iso_no_sup, &region_tp53, &unhyd_map, &splice);
        assert_eq!(unhyd_out.len(), 2);
        assert!(unhyd_out.iter().all(|v| !v.null_record && v.start_pos >= 1));

        // 2. `[del@P-1, snv@P]` non-duplication: surviving `snv@P` keeps `start_pos == P`, so `ref_variant` is NOT added.
        let mut nondup_map = crate::javamap::JavaHashMap::new();
        let mut nondup_vars = Vars::new();
        nondup_vars.ref_variant = Some(mk_variant(7_578_712, 7_578_712, "A", "A", 0.3));
        nondup_vars
            .variants
            .push(mk_variant(7_578_711, 7_578_712, "CA", "C", 0.4));
        nondup_vars
            .variants
            .push(mk_variant(7_578_712, 7_578_712, "A", "G", 0.3));
        nondup_map.insert(7_578_712, nondup_vars);
        let mut off_p_cfg = pileup_config();
        off_p_cfg.compat = crate::config::COMPAT_OFF;
        let nondup_out = simple_post_process(&off_p_cfg, &region_tp53, &nondup_map, &splice);
        assert_eq!(nondup_out.len(), 2);
        assert_eq!(nondup_out[0].start_pos, 7_578_711);
        assert_eq!(nondup_out[1].start_pos, 7_578_712);
        assert_eq!(nondup_out[1].var_allele, "G");

        // 3. `AA>AC` `adj_complex`-shifted `-p` retention + `recompute_vartype_after_adj_complex` (`clintval/vartovcf#109`):
        //    Before `adj_complex`, `start_pos == 7578712 == position`; after `adj_complex_compat(true)`,
        //    it shifts to `7578713 A>C` (`var_type == "SNV"`), so `ref_variant` at `7578712` IS emitted!
        let mut adj_p_map = crate::javamap::JavaHashMap::new();
        let mut adj_p_vars = Vars::new();
        adj_p_vars.ref_variant = Some(mk_variant(7_578_712, 7_578_712, "A", "A", 0.5));
        adj_p_vars
            .variants
            .push(mk_variant(7_578_712, 7_578_713, "AA", "AC", 0.5));
        adj_p_map.insert(7_578_712, adj_p_vars);

        let vj_adj_p = simple_post_process(&vj_p_cfg, &region_tp53, &adj_p_map, &splice);
        assert_eq!(vj_adj_p.len(), 1);
        assert_eq!(vj_adj_p[0].start_pos, 7_578_713);
        assert_eq!(vj_adj_p[0].var_type, "Complex");

        let off_adj_p = simple_post_process(&off_p_cfg, &region_tp53, &adj_p_map, &splice);
        assert_eq!(off_adj_p.len(), 2);
        assert_eq!(off_adj_p[0].start_pos, 7_578_712);
        assert_eq!(off_adj_p[0].var_type, "");
        assert_eq!(off_adj_p[1].start_pos, 7_578_713);
        assert_eq!(off_adj_p[1].ref_allele, "A");
        assert_eq!(off_adj_p[1].var_allele, "C");
        assert_eq!(off_adj_p[1].var_type, "SNV");

        // 4. Switch 29 (`recompute_vartype_after_adj_complex`) & Switch 26 (`sanitize_iupac_alt_alleles`)
        //    in non-pileup `simple_post_process` (`SNV`, `Insertion`, `Deletion`, genuine `Complex`, and `REF == ALT` skip):
        let region_s = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut vt_map = crate::javamap::JavaHashMap::new();
        let mut v_snv = Vars::new();
        v_snv.variants.push(mk_variant(110, 111, "GA", "GT", 0.4));
        vt_map.insert(110, v_snv);
        let mut v_ins = Vars::new();
        v_ins
            .variants
            .push(mk_variant(120, 122, "CAA", "CAATG", 0.4));
        vt_map.insert(120, v_ins);
        let mut v_del = Vars::new();
        v_del
            .variants
            .push(mk_variant(130, 134, "CAATG", "CAA", 0.4));
        vt_map.insert(130, v_del);
        let mut v_mnv = Vars::new();
        v_mnv.variants.push(mk_variant(140, 141, "AC", "GT", 0.4));
        vt_map.insert(140, v_mnv);
        let mut v_eq = Vars::new();
        v_eq.variants.push(mk_variant(150, 151, "AA", "AA", 0.4));
        vt_map.insert(150, v_eq);

        let mut non_p_vj = Config::default();
        non_p_vj.emit_intermediate = true;
        let out_vj_vt = simple_post_process(&non_p_vj, &region_s, &vt_map, &splice);
        assert!(out_vj_vt.iter().all(|v| v.var_type == "Complex"));

        let mut non_p_off = non_p_vj.clone();
        non_p_off.compat = crate::config::COMPAT_OFF;
        let out_off_vt = simple_post_process(&non_p_off, &region_s, &vt_map, &splice);
        let typed: Vec<(i32, &str, &str, &str)> = out_off_vt
            .iter()
            .map(|v| {
                (
                    v.start_pos,
                    v.ref_allele.as_str(),
                    v.var_allele.as_str(),
                    v.var_type.as_str(),
                )
            })
            .collect();
        assert_eq!(
            typed,
            vec![
                (111, "A", "T", "SNV"),
                (122, "A", "ATG", "Insertion"),
                (132, "ATG", "A", "Deletion"),
                (140, "AC", "GT", "Complex"),
            ]
        );
    }
}
