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
    simple_post_process_core(
        config,
        region,
        aligned_variants
            .iter()
            .map(|(&pos, vars)| (pos, std::borrow::Cow::Borrowed(vars))),
        splice,
    )
}

/// Owned fast-path counterpart of [`simple_post_process`], moving finalized [`Variant`]
/// and SV strings directly out of `aligned_variants` without cloning.
pub(crate) fn simple_post_process_owned(
    config: &Config,
    region: &Region,
    aligned_variants: crate::javamap::JavaHashMap<i32, Vars>,
    splice: &rustc_hash::FxHashSet<String>,
) -> Vec<Variant> {
    simple_post_process_core(
        config,
        region,
        aligned_variants
            .into_entries()
            .map(|(pos, vars)| (pos, std::borrow::Cow::Owned(vars))),
        splice,
    )
}

fn simple_post_process_core<'a>(
    config: &Config,
    region: &Region,
    entries: impl IntoIterator<Item = (i32, std::borrow::Cow<'a, Vars>)>,
    splice: &rustc_hash::FxHashSet<String>,
) -> Vec<Variant> {
    let mut vrefs = Vec::new();

    // Iterate entries in JavaHashMap table-bucket iteration order
    // (ports SimplePostProcessModule.java:35; see DIVERGENCES.md UB-3).
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
                vrefs.push(vref);
            } else {
                // ports SimplePostProcessModule.java:55-58: under -p an entry
                // with no variants and no reference variant still prints, via
                // `SimpleOutputVariant(null, region, sv, position)` whose fields
                // keep the printer defaults (genotype "", bias "0;0", flank seqs
                // ""). CLAUDE.md §5: live-oracle output is replicated as-is and
                // tracked upstream, never silently fixed.
                vrefs.push(make_null_record(position));
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
                        vrefs.push(ref_var);
                    } else {
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
        let mut buckets: BTreeMap<i32, Variant> = BTreeMap::new();
        for v in vrefs {
            let v_af = crate::types::round_half_even(v.frequency, 4);
            match buckets.entry(v.start_pos) {
                std::collections::btree_map::Entry::Vacant(e) => {
                    e.insert(v);
                }
                std::collections::btree_map::Entry::Occupied(mut e) => {
                    let cur_af = crate::types::round_half_even(e.get().frequency, 4);
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
                for crispr_site in [None, Some(150)] {
                    let mut cfg = config.clone();
                    cfg.pileup = pileup;
                    cfg.emit_intermediate = emit_intermediate;
                    cfg.crispr_cutting_site = crispr_site;
                    let borrowed = simple_post_process(&cfg, &region, &parity_map, &splice);
                    let owned =
                        simple_post_process_owned(&cfg, &region, parity_map.clone(), &splice);
                    assert_eq!(owned.len(), borrowed.len());
                    for (got, want) in owned.iter().zip(borrowed.iter()) {
                        assert_variant_fields_eq(got, want);
                    }
                }
            }
        }
    }
}
