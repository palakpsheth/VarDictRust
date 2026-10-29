//! Table output serialization formats — intermediate .var output formatters.
//! Ports SimpleOutputVariant.java.

use crate::regions::Region;
use crate::variations::Variant;

/// Print the `-h/--header` column line for the active mode to `w`.
/// Ports CmdParser.java:67 (`printHeader`) and the per-mode headers:
/// SimpleMode.java:114-127, SomaticMode.java:134-150, AmpliconMode.java:135-146.
/// Java prints this unconditionally of output format, before any variant rows.
pub fn print_header<W: std::io::Write>(
    w: &mut W,
    config: &crate::config::Config,
) -> std::io::Result<()> {
    const SIMPLE: &[&str] = &[
        "Sample",
        "Gene",
        "Chr",
        "Start",
        "End",
        "Ref",
        "Alt",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "shift3",
        "MSI",
        "MSI_NT",
        "NM",
        "HiCnt",
        "HiCov",
        "5pFlankSeq",
        "3pFlankSeq",
        "Seg",
        "VarType",
        "Duprate",
        "SV_info",
    ];
    // Paired mode repeats the per-sample block; no CRISPR suffix in v1.8.3.
    const SOMATIC: &[&str] = &[
        "Sample",
        "Gene",
        "Chr",
        "Start",
        "End",
        "Ref",
        "Alt",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "NM",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "NM",
        "shift3",
        "MSI",
        "MSI_NT",
        "5pFlankSeq",
        "3pFlankSeq",
        "Seg",
        "VarLabel",
        "VarType",
        "Duprate1",
        "SV_info1",
        "Duprate2",
        "SV_info2",
    ];
    const AMPLICON: &[&str] = &[
        "Sample",
        "Gene",
        "Chr",
        "Start",
        "End",
        "Ref",
        "Alt",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "shift3",
        "MSI",
        "MSI_NT",
        "NM",
        "HiCnt",
        "HiCov",
        "5pFlankSeq",
        "3pFlankSeq",
        "Seg",
        "VarType",
        "GoodVarCount",
        "TotalVarCount",
        "Nocov",
        "Ampflag",
    ];

    let mut names: Vec<&str> = if config.bam.contains('|') {
        SOMATIC.to_vec()
    } else if config.amplicon.is_some() {
        AMPLICON.to_vec()
    } else {
        SIMPLE.to_vec()
    };
    if !config.bam.contains('|')
        && config.amplicon.is_none()
        && matches!(config.crispr_cutting_site, Some(n) if n != 0)
    {
        // SimpleMode.java:124-126: CRISPR mode appends a column.
        names.push("CRISPR");
    }
    writeln!(w, "{}", names.join("\t"))?;
    Ok(())
}

/// Write simple 36-column representation of finalized Variants to a tab-delimited string.
/// Ports SimpleOutputVariant.java's create_simple_variant_36columns.
pub fn simple_36(variant: &Variant, region: &Region, sample_name: &str, sv: &str) -> String {
    let region_str = format!("{}:{}-{}", region.chr, region.start.0, region.end.0);
    let sv_str = if sv.is_empty() {
        "0".to_string()
    } else {
        sv.to_string()
    };

    let freq_str = fmt_4dp(variant.frequency);
    let pmean_str = fmt_1dp(variant.mean_position);
    let qual_str = fmt_1dp(variant.mean_quality);
    let mapq_str = fmt_1dp(variant.mean_mapping_quality);
    let qratio_str = fmt_3dp(variant.hq_lq_ratio);
    let hifreq_str = fmt_4dp(variant.hi_freq);
    let extrafreq_str = fmt_4dp(variant.extra_freq);
    let msi_str = fmt_3dp(variant.msi);
    let nm_str = fmt_nm(variant.nm);
    let duprate_str = fmt_1dp(variant.duprate);

    // Null records render the printer field defaults, not the normal
    // null→"0" / empty→"0" mapping (SimpleOutputVariant.java:50-52 vs :64,
    // :81-82; OutputVariant.java:22-23).
    let left_seq = if variant.null_record {
        ""
    } else if variant.left_seq.is_empty() {
        "0"
    } else {
        &variant.left_seq
    };
    let right_seq = if variant.null_record {
        ""
    } else if variant.right_seq.is_empty() {
        "0"
    } else {
        &variant.right_seq
    };
    let genotype = if variant.null_record {
        ""
    } else {
        variant.genotype.as_deref().unwrap_or("0")
    };

    let pstd = if variant.pstd { 1 } else { 0 };
    let qstd = if variant.qstd { 1 } else { 0 };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        sample_name,
        region.gene,
        region.chr,
        variant.start_pos,
        variant.end_pos,
        variant.ref_allele,
        variant.var_allele,
        variant.total_pos_cov,
        variant.pos_coverage,
        variant.ref_fwd_cov,
        variant.ref_rev_cov,
        variant.fwd,
        variant.rev,
        genotype,
        freq_str,
        variant.bias,
        pmean_str,
        pstd,
        qual_str,
        qstd,
        mapq_str,
        qratio_str,
        hifreq_str,
        extrafreq_str,
        variant.shift3,
        msi_str,
        variant.msint,
        nm_str,
        variant.hicnt,
        variant.hicov,
        left_seq,
        right_seq,
        region_str,
        variant.var_type,
        duprate_str,
        sv_str
    )
}

/// Write simple 38-column representation of finalized Variants including Fisher Exact test results.
/// Ports SimpleOutputVariant.java's create_simple_variant_38columns.
pub fn simple_38(variant: &Variant, region: &Region, sample_name: &str, sv: &str) -> String {
    let region_str = format!("{}:{}-{}", region.chr, region.start.0, region.end.0);
    let sv_str = if sv.is_empty() {
        "0".to_string()
    } else {
        sv.to_string()
    };

    let fisher = crate::stats::FisherExact::new_java(
        variant.ref_fwd_cov as i32,
        variant.ref_rev_cov as i32,
        variant.fwd as i32,
        variant.rev as i32,
    );
    let pvalue = fisher.get_p_value();
    let oddratio = fisher.get_odds_ratio();

    let hifreq_str = if variant.hi_freq == 0.0 {
        "0".to_string()
    } else {
        format!("{:.4}", crate::types::round_half_even(variant.hi_freq, 4))
    };

    let nm_val = if variant.nm > 0.0 { variant.nm } else { 0.0 };
    let nm_str = if nm_val == 0.0 {
        "0".to_string()
    } else {
        format!("{:.1}", crate::types::round_half_even(nm_val, 1))
    };

    // Same null-record printer-defaults rule as simple_36.
    let left_seq = if variant.null_record {
        ""
    } else if variant.left_seq.is_empty() {
        "0"
    } else {
        &variant.left_seq
    };
    let right_seq = if variant.null_record {
        ""
    } else if variant.right_seq.is_empty() {
        "0"
    } else {
        &variant.right_seq
    };
    let genotype = if variant.null_record {
        ""
    } else {
        variant.genotype.as_deref().unwrap_or("0")
    };

    let pstd = if variant.pstd { 1 } else { 0 };
    let qstd = if variant.qstd { 1 } else { 0 };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        sample_name,
        region.gene,
        region.chr,
        variant.start_pos,
        variant.end_pos,
        variant.ref_allele,
        variant.var_allele,
        variant.total_pos_cov,
        variant.pos_coverage,
        variant.ref_fwd_cov,
        variant.ref_rev_cov,
        variant.fwd,
        variant.rev,
        genotype,
        get_rounded_value_to_print(variant.frequency, 4),
        variant.bias,
        get_rounded_value_to_print(variant.mean_position, 1),
        pstd,
        get_rounded_value_to_print(variant.mean_quality, 1),
        qstd,
        get_rounded_value_to_print(pvalue, 5),
        oddratio,
        get_rounded_value_to_print(variant.mean_mapping_quality, 1),
        get_rounded_value_to_print(variant.hq_lq_ratio, 3),
        hifreq_str,
        get_rounded_value_to_print(variant.extra_freq, 4),
        variant.shift3,
        get_rounded_value_to_print(variant.msi, 3),
        variant.msint,
        nm_str,
        variant.hicnt,
        variant.hicov,
        left_seq,
        right_seq,
        region_str,
        variant.var_type,
        get_rounded_value_to_print(variant.duprate, 2),
        sv_str
    )
}

/// Write amplicon 38-column representation (no fisher).
/// Ports AmpliconOutputVariant.java's create_amplicon_variant_38columns.
/// Columns 1-34 match the simple format; columns 35-38 are
/// GDAMP, TLAMP, NCAMP, AMPFLAG (replacing simple's duprate/sv).
/// When `debug` is set, appends `"\t" + DEBUG` where
/// DEBUG = variant.DEBUG + goodOutput + badOutput (AmpliconOutputVariant.java:113-138).
pub fn amplicon_38(
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) -> String {
    let variant = &av.variant;
    let region_str = &av.region;

    let left_seq = if variant.left_seq.is_empty() {
        "0"
    } else {
        &variant.left_seq
    };
    let right_seq = if variant.right_seq.is_empty() {
        "0"
    } else {
        &variant.right_seq
    };
    let genotype = variant.genotype.as_deref().unwrap_or("0");
    let pstd = if variant.pstd { 1 } else { 0 };
    let qstd = if variant.qstd { 1 } else { 0 };
    let ampflag = if av.ampflag { 1 } else { 0 };

    let line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        sample_name,
        region.gene,
        region.chr,
        variant.start_pos,
        variant.end_pos,
        variant.ref_allele,
        variant.var_allele,
        variant.total_pos_cov,
        variant.pos_coverage,
        variant.ref_fwd_cov,
        variant.ref_rev_cov,
        variant.fwd,
        variant.rev,
        genotype,
        fmt_4dp(variant.frequency),
        variant.bias,
        fmt_1dp(variant.mean_position),
        pstd,
        fmt_1dp(variant.mean_quality),
        qstd,
        fmt_1dp(variant.mean_mapping_quality),
        fmt_3dp(variant.hq_lq_ratio),
        fmt_4dp(variant.hi_freq),
        fmt_4dp(variant.extra_freq),
        variant.shift3,
        fmt_3dp_msi(variant.msi),
        variant.msint,
        fmt_nm(variant.nm),
        variant.hicnt,
        variant.hicov,
        left_seq,
        right_seq,
        region_str,
        variant.var_type,
        av.gdamp,
        av.tlamp,
        av.ncamp,
        ampflag
    );

    if !debug {
        return line;
    }
    // ports AmpliconOutputVariant.java:113-125
    let mut debug_col = variant.debug.clone();
    for (i, (v, r)) in av.good_variants.iter().enumerate() {
        debug_col.push_str(&format!(
            "\tGood{i} {} {}",
            debug_amp_variant(v.as_ref()),
            r
        ));
    }
    for (i, (v, r)) in av.bad_variants.iter().enumerate() {
        debug_col.push_str(&format!("\tBad{i} {} {}", debug_amp_variant(v.as_ref()), r));
    }
    format!("{}\t{}", line, debug_col)
}

/// Write amplicon 40-column representation with Fisher exact test (when config.fisher is set).
/// Ports AmpliconOutputVariant.java's create_amplicon_variant_40columns (lines 145-201).
/// Columns 1-20 match amplicon_38, columns 21-22 are pvalue and oddratio,
/// columns 23-36 match columns 21-34 of amplicon_38, and columns 37-40 are
/// GDAMP, TLAMP, NCAMP, AMPFLAG.
/// When `debug` is set, appends `"\t" + DEBUG` where
/// DEBUG = variant.DEBUG + goodOutput + badOutput.
pub fn amplicon_40(
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) -> String {
    let variant = &av.variant;
    let region_str = &av.region;

    let fisher = crate::stats::FisherExact::new_java(
        variant.ref_fwd_cov as i32,
        variant.ref_rev_cov as i32,
        variant.fwd as i32,
        variant.rev as i32,
    );
    let pvalue = fisher.get_p_value();
    let oddratio = fisher.get_odds_ratio();

    let hifreq_str = if variant.hi_freq == 0.0 {
        "0".to_string()
    } else {
        format!("{:.4}", crate::types::round_half_even(variant.hi_freq, 4))
    };

    let nm_val = if variant.nm > 0.0 { variant.nm } else { 0.0 };
    let nm_str = if nm_val == 0.0 {
        "0".to_string()
    } else {
        format!("{:.1}", crate::types::round_half_even(nm_val, 1))
    };

    let left_seq = if variant.left_seq.is_empty() {
        "0"
    } else {
        &variant.left_seq
    };
    let right_seq = if variant.right_seq.is_empty() {
        "0"
    } else {
        &variant.right_seq
    };
    let genotype = variant.genotype.as_deref().unwrap_or("0");
    let pstd = if variant.pstd { 1 } else { 0 };
    let qstd = if variant.qstd { 1 } else { 0 };
    let ampflag = if av.ampflag { 1 } else { 0 };

    let line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        sample_name,
        region.gene,
        region.chr,
        variant.start_pos,
        variant.end_pos,
        variant.ref_allele,
        variant.var_allele,
        variant.total_pos_cov,
        variant.pos_coverage,
        variant.ref_fwd_cov,
        variant.ref_rev_cov,
        variant.fwd,
        variant.rev,
        genotype,
        get_rounded_value_to_print(variant.frequency, 4),
        variant.bias,
        get_rounded_value_to_print(variant.mean_position, 1),
        pstd,
        get_rounded_value_to_print(variant.mean_quality, 1),
        qstd,
        get_rounded_value_to_print(pvalue, 5),
        oddratio,
        get_rounded_value_to_print(variant.mean_mapping_quality, 1),
        get_rounded_value_to_print(variant.hq_lq_ratio, 3),
        hifreq_str,
        get_rounded_value_to_print(variant.extra_freq, 4),
        variant.shift3,
        get_rounded_value_to_print(variant.msi, 3),
        variant.msint,
        nm_str,
        variant.hicnt,
        variant.hicov,
        left_seq,
        right_seq,
        region_str,
        variant.var_type,
        av.gdamp,
        av.tlamp,
        av.ncamp,
        ampflag
    );

    if !debug {
        return line;
    }
    // ports AmpliconOutputVariant.java:113-125
    let mut debug_col = variant.debug.clone();
    for (i, (v, r)) in av.good_variants.iter().enumerate() {
        debug_col.push_str(&format!(
            "\tGood{i} {} {}",
            debug_amp_variant(v.as_ref()),
            r
        ));
    }
    for (i, (v, r)) in av.bad_variants.iter().enumerate() {
        debug_col.push_str(&format!("\tBad{i} {} {}", debug_amp_variant(v.as_ref()), r));
    }
    format!("{}\t{}", line, debug_col)
}

/// One space-joined entry of the -D amplicon debug column.
/// Ports `AmpliconOutputVariant.debugAmpVariant` (AmpliconOutputVariant.java:258-278):
/// 17 fields; zero-valued freq fields print as "0", position/quality/mapq/qratio are
/// always decimal-formatted. Java dereferences a null variant here (NPE on the
/// no-data bad-variant path); we render an all-zero record instead.
fn debug_amp_variant(v: Option<&Variant>) -> String {
    let default = Variant::new();
    let v = v.unwrap_or(&default);
    let genotype = v.genotype.as_deref().unwrap_or("0");
    let pstd = if v.pstd { 1 } else { 0 };
    let qstd = if v.qstd { 1 } else { 0 };
    let mean_pos = format!("{:.1}", crate::types::round_half_even(v.mean_position, 1));
    let mean_qual = format!("{:.1}", crate::types::round_half_even(v.mean_quality, 1));
    let mapq = format!(
        "{:.1}",
        crate::types::round_half_even(v.mean_mapping_quality, 1)
    );
    let qratio = format!("{:.3}", crate::types::round_half_even(v.hq_lq_ratio, 3));
    format!(
        "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        v.total_pos_cov,
        v.pos_coverage,
        v.ref_fwd_cov,
        v.ref_rev_cov,
        v.fwd,
        v.rev,
        genotype,
        fmt_4dp(v.frequency),
        v.bias,
        mean_pos,
        pstd,
        mean_qual,
        qstd,
        mapq,
        qratio,
        fmt_4dp(v.hi_freq),
        fmt_4dp(v.extra_freq),
    )
}

/// MSI formatter: 3dp but "0" when zero (matches Java msi == 0 ? 0 : format).
fn fmt_3dp_msi(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!("{:.3}", crate::types::round_half_even(v, 3))
    }
}

/// Java-compatible formatter that rounds a float using banker's rounding
/// to a maximum of `digits` decimal places, and strips any trailing zeros.
/// Ports `Utils.getRoundedValueToPrint`. Shared with the VCF writer (vcf.rs).
pub(crate) fn get_rounded_value_to_print(value: f64, digits: i32) -> String {
    if value.is_nan() || value.is_infinite() {
        return value.to_string();
    }
    let rounded = crate::types::round_half_even(value, digits);
    if rounded == 0.0 {
        "0".to_string()
    } else if rounded == rounded.round() {
        format!("{:.0}", rounded)
    } else {
        let mut s = format!("{:.*}", digits as usize, rounded);
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
        }
        s
    }
}

fn fmt_1dp(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!("{:.1}", crate::types::round_half_even(v, 1))
    }
}
fn fmt_3dp(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!("{:.3}", crate::types::round_half_even(v, 3))
    }
}
fn fmt_4dp(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!("{:.4}", crate::types::round_half_even(v, 4))
    }
}
fn fmt_nm(v: f64) -> String {
    if v > 0.0 {
        format!("{:.1}", crate::types::round_half_even(v, 1))
    } else {
        "0".to_string()
    }
}

fn variant_18cols(v: &Variant) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        v.total_pos_cov,
        v.pos_coverage,
        v.ref_fwd_cov,
        v.ref_rev_cov,
        v.fwd,
        v.rev,
        v.genotype.as_deref().unwrap_or("0"),
        fmt_4dp(v.frequency),
        v.bias,
        fmt_1dp(v.mean_position),
        if v.pstd { 1 } else { 0 },
        fmt_1dp(v.mean_quality),
        if v.qstd { 1 } else { 0 },
        fmt_1dp(v.mean_mapping_quality),
        fmt_3dp(v.hq_lq_ratio),
        fmt_4dp(v.hi_freq),
        fmt_4dp(v.extra_freq),
        fmt_nm(v.nm),
    )
}

fn variant_18cols_empty() -> String {
    "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0".to_string()
}

/// Write somatic 55-column representation for paired (tumor/normal) variants.
/// Ports SomaticOutputVariant.java's create_somatic_variant_55columns.
pub fn somatic_55(
    tumor: &Variant,
    normal: Option<&Variant>,
    region: &Region,
    sample_name: &str,
    status: &str,
    sv1: &str,
    sv2: &str,
    end_is_variant2: bool,
) -> String {
    let region_str = format!("{}:{}-{}", region.chr, region.start.0, region.end.0);
    let sv1_str = if sv1.is_empty() {
        "0".to_string()
    } else {
        sv1.to_string()
    };
    let sv2_str = if sv2.is_empty() {
        "0".to_string()
    } else {
        sv2.to_string()
    };
    // endVariant supplies shift3/msi/msint/left/right seq (Java SomaticOutputVariant).
    let end_var = if end_is_variant2 {
        normal.unwrap_or(tumor)
    } else {
        tumor
    };
    let left_seq = if end_var.left_seq.is_empty() {
        "0"
    } else {
        &end_var.left_seq
    };
    let right_seq = if end_var.right_seq.is_empty() {
        "0"
    } else {
        &end_var.right_seq
    };

    let normal_cols = match normal {
        Some(n) => variant_18cols(n),
        None => variant_18cols_empty(),
    };

    let duprate1 = fmt_1dp(tumor.duprate);
    let duprate2 = match normal {
        Some(n) => fmt_1dp(n.duprate),
        None => "0".to_string(),
    };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        sample_name,
        region.gene,
        region.chr,
        tumor.start_pos,
        tumor.end_pos,
        tumor.ref_allele,
        tumor.var_allele,
        variant_18cols(tumor),
        normal_cols,
        end_var.shift3,
        fmt_3dp(end_var.msi),
        end_var.msint,
        left_seq,
        right_seq,
        region_str,
        status,
        tumor.var_type,
    ) + &format!("\t{}\t{}\t{}\t{}", duprate1, sv1_str, duprate2, sv2_str)
}

/// 20-column fisher sample block: the 18 base columns plus pvalue and oddratio,
/// using trailing-zero-stripping formatting (Java getRoundedValueToPrint).
fn variant_20cols_fisher(v: &Variant, pvalue: f64, oddratio: &str) -> String {
    let nm = if v.nm > 0.0 { v.nm } else { 0.0 };
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        v.total_pos_cov,
        v.pos_coverage,
        v.ref_fwd_cov,
        v.ref_rev_cov,
        v.fwd,
        v.rev,
        v.genotype.as_deref().unwrap_or("0"),
        get_rounded_value_to_print(v.frequency, 4),
        v.bias,
        get_rounded_value_to_print(v.mean_position, 1),
        if v.pstd { 1 } else { 0 },
        get_rounded_value_to_print(v.mean_quality, 1),
        if v.qstd { 1 } else { 0 },
        get_rounded_value_to_print(v.mean_mapping_quality, 1),
        get_rounded_value_to_print(v.hq_lq_ratio, 3),
        get_rounded_value_to_print(v.hi_freq, 4),
        get_rounded_value_to_print(v.extra_freq, 4),
        get_rounded_value_to_print(nm, 1),
        get_rounded_value_to_print(pvalue, 5),
        oddratio,
    )
}

fn variant_20cols_fisher_empty(pvalue: f64, oddratio: &str) -> String {
    format!(
        "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t{}\t{}",
        get_rounded_value_to_print(pvalue, 5),
        oddratio,
    )
}

/// Write somatic 61-column representation (with --fisher).
/// Ports SomaticOutputVariant.java's create_somatic_variant_61columns.
pub fn somatic_61(
    tumor: &Variant,
    normal: Option<&Variant>,
    region: &Region,
    sample_name: &str,
    status: &str,
    sv1: &str,
    sv2: &str,
    end_is_variant2: bool,
) -> String {
    let region_str = format!("{}:{}-{}", region.chr, region.start.0, region.end.0);
    let sv1_str = if sv1.is_empty() {
        "0".to_string()
    } else {
        sv1.to_string()
    };
    let sv2_str = if sv2.is_empty() {
        "0".to_string()
    } else {
        sv2.to_string()
    };
    // endVariant supplies shift3/msi/msint/left/right seq (Java SomaticOutputVariant).
    let end_var = if end_is_variant2 {
        normal.unwrap_or(tumor)
    } else {
        tumor
    };
    let left_seq = if end_var.left_seq.is_empty() {
        "0"
    } else {
        &end_var.left_seq
    };
    let right_seq = if end_var.right_seq.is_empty() {
        "0"
    } else {
        &end_var.right_seq
    };

    let fisher1 = crate::stats::FisherExact::new_java(
        tumor.ref_fwd_cov as i32,
        tumor.ref_rev_cov as i32,
        tumor.fwd as i32,
        tumor.rev as i32,
    );
    let pvalue1 = fisher1.get_p_value();
    let oddratio1 = fisher1.get_odds_ratio();

    let normal_cols = match normal {
        Some(n) => {
            let f2 = crate::stats::FisherExact::new_java(
                n.ref_fwd_cov as i32,
                n.ref_rev_cov as i32,
                n.fwd as i32,
                n.rev as i32,
            );
            variant_20cols_fisher(n, f2.get_p_value(), &f2.get_odds_ratio())
        }
        None => {
            let f2 = crate::stats::FisherExact::new_java(0, 0, 0, 0);
            variant_20cols_fisher_empty(f2.get_p_value(), &f2.get_odds_ratio())
        }
    };

    let (tvar_cov, ttotal_cov) = (tumor.pos_coverage, tumor.total_pos_cov);
    let (nvar_cov, ntotal_cov) = match normal {
        Some(n) => (n.pos_coverage, n.total_pos_cov),
        None => (0, 0),
    };
    let tref = (ttotal_cov - tvar_cov).max(0);
    let rref = (ntotal_cov - nvar_cov).max(0);
    let fisher_total = crate::stats::FisherExact::new_java(
        tvar_cov as i32,
        tref as i32,
        nvar_cov as i32,
        rref as i32,
    );
    let pvalue_greater = fisher_total.get_p_value_greater();
    let pvalue_less = fisher_total.get_p_value_less();
    let total_pvalue = if pvalue_less < pvalue_greater {
        pvalue_less
    } else {
        pvalue_greater
    };
    let total_oddratio = fisher_total.get_odds_ratio();

    let duprate1 = get_rounded_value_to_print(tumor.duprate, 1);
    let duprate2 = match normal {
        Some(n) => get_rounded_value_to_print(n.duprate, 1),
        None => "0".to_string(),
    };

    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        sample_name,
        region.gene,
        region.chr,
        tumor.start_pos,
        tumor.end_pos,
        tumor.ref_allele,
        tumor.var_allele,
        variant_20cols_fisher(tumor, pvalue1, &oddratio1),
        normal_cols,
        end_var.shift3,
        fmt_3dp(end_var.msi),
        end_var.msint,
        left_seq,
        right_seq,
        region_str,
    ) + &format!(
        "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        status,
        tumor.var_type,
        duprate1,
        sv1_str,
        duprate2,
        sv2_str,
        get_rounded_value_to_print(total_pvalue, 5),
        total_oddratio,
    )
}

#[cfg(test)]
mod schema_tripwire_tests {
    // Tripwires against schema drift: pin the EXACT intermediate-table column
    // counts (Simple 36/38, Amplicon 38, Somatic 55/61 — the parity contract in
    // CLAUDE.md §5) plus a stable prefix/suffix of field values per row. The
    // .var rows are headerless, so the "columns" are pinned as tab-separated
    // values on one representative record fixture.
    //
    // BY DESIGN these pins cannot catch middle-column swaps: only the column
    // COUNT and the first-7 / last-4 boundary fields are asserted. Value-level
    // drift in the middle columns is covered by the corpus comparators
    // (harness/compare_var.py), which check every cell.
    use super::*;
    use crate::postprocess::amplicon::AmpliconVariant;
    use crate::regions::Region;
    use crate::types::Pos1;

    fn region() -> Region {
        Region {
            chr: "chr1".to_string(),
            start: Pos1(933866),
            end: Pos1(934466),
            gene: "testbed".to_string(),
            amplicon: None,
        }
    }

    fn variant() -> Variant {
        let mut v = Variant::new();
        v.start_pos = 100;
        v.end_pos = 100;
        v.ref_allele = "G".to_string();
        v.var_allele = "C".to_string();
        v.var_type = "SNV".to_string();
        v
    }

    // Columns 1-7 shared by every schema: sample, gene, chr, start, end, ref, alt.
    const PREFIX_7: &str = "SAMPLE\ttestbed\tchr1\t100\t100\tG\tC";

    #[test]
    fn simple_36_pins_column_count_and_boundary_fields() {
        let line = simple_36(&variant(), &region(), "SAMPLE", "");
        assert_eq!(line.split('\t').count(), 36);
        assert!(line.starts_with(PREFIX_7), "prefix drifted: {line}");
        // ... region, var_type, duprate("0"), sv("0")
        assert!(
            line.ends_with("\tchr1:933866-934466\tSNV\t0\t0"),
            "suffix drifted: {line}"
        );
    }

    /// Java's null record keeps the printer field defaults: genotype "" (col 14),
    /// bias "0;0" (col 16), flank seqs "" (cols 31/32) — NOT the normal
    /// null→"0"/empty→"0" mapping (SimpleOutputVariant.java:50-52 vs :64/:81-82).
    /// Byte shape pinned against the P7 live-oracle capture: cols 4/5 = 6447,
    /// cols 6/7 empty, col 14 empty, col 15 "0", col 16 "0;0", cols 31/32
    /// empty, col 34 empty, cols 35/36 "0".
    #[test]
    fn simple_36_null_record_renders_printer_defaults() {
        let mut v = Variant::new();
        v.start_pos = 6447;
        v.end_pos = 6447;
        v.bias = "0;0".to_string();
        v.null_record = true;
        let line = simple_36(&v, &region(), "SAMPLE", "");
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 36);
        assert_eq!(cols[3], "6447");
        assert_eq!(cols[4], "6447");
        assert_eq!(cols[5], ""); // ref allele default
        assert_eq!(cols[6], ""); // var allele default
        assert_eq!(cols[13], ""); // genotype "" — not "0"
        assert_eq!(cols[14], "0"); // af == 0 → "0"
        assert_eq!(cols[15], "0;0"); // bias field default
        assert_eq!(cols[30], ""); // 5pFlankSeq "" — not "0"
        assert_eq!(cols[31], ""); // 3pFlankSeq "" — not "0"
        assert_eq!(cols[33], ""); // varType default
    }

    /// The pos-0 companion row is a plain default Variant through the NORMAL
    /// path: every cell maps to "0" including genotype/bias/seqs.
    #[test]
    fn simple_36_default_variant_renders_all_zeros() {
        let v = Variant::new();
        let line = simple_36(&v, &region(), "SAMPLE", "");
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 36);
        assert_eq!(cols[3], "0");
        assert_eq!(cols[4], "0");
        for i in [13usize, 15, 30, 31] {
            assert_eq!(cols[i], "0", "col {} must be \"0\"", i + 1);
        }
    }

    #[test]
    fn simple_38_pins_column_count_and_boundary_fields() {
        let line = simple_38(&variant(), &region(), "SAMPLE", "");
        assert_eq!(line.split('\t').count(), 38);
        assert!(line.starts_with(PREFIX_7), "prefix drifted: {line}");
        // Same tail shape as 36-col; duprate via getRoundedValueToPrint(0.0, 2) = "0".
        assert!(
            line.ends_with("\tchr1:933866-934466\tSNV\t0\t0"),
            "suffix drifted: {line}"
        );
    }

    #[test]
    fn amplicon_38_pins_column_count_and_boundary_fields() {
        let av = AmpliconVariant {
            variant: variant(),
            gdamp: 1,
            tlamp: 1,
            ncamp: 0,
            ampflag: false,
            region: "chr1:933867-934466".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };
        let line = amplicon_38(&av, &region(), "SAMPLE", false);
        assert_eq!(line.split('\t').count(), 38);
        assert!(line.starts_with(PREFIX_7), "prefix drifted: {line}");
        // Columns 35-38: GDAMP, TLAMP, NCAMP, AMPFLAG (no debug column).
        assert!(line.ends_with("\t1\t1\t0\t0"), "suffix drifted: {line}");
    }

    #[test]
    fn test_amplicon_40_columns() {
        let mut v = variant();
        v.ref_fwd_cov = 10;
        v.ref_rev_cov = 10;
        v.fwd = 5;
        v.rev = 5;
        v.frequency = 0.3333;
        v.mean_position = 25.0;
        v.mean_quality = 30.0;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 1.0;
        v.hi_freq = 0.3333;
        v.extra_freq = 0.0;
        v.shift3 = 0;
        v.msi = 0.0;
        v.msint = 1;
        v.nm = 1.0;
        v.hicnt = 10;
        v.hicov = 30;

        let av = AmpliconVariant {
            variant: v,
            gdamp: 2,
            tlamp: 3,
            ncamp: 1,
            ampflag: true,
            region: "chr1:933867-934466".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };
        let line = amplicon_40(&av, &region(), "SAMPLE", false);
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 40);
        assert!(line.starts_with(PREFIX_7), "prefix drifted: {line}");
        // Column 21 (index 20) is pvalue, Column 22 (index 21) is oddratio
        assert_eq!(cols[20], "1");
        assert_eq!(cols[21], "1");
        // Columns 37-40 (indices 36-39): GDAMP, TLAMP, NCAMP, AMPFLAG
        assert_eq!(cols[36], "2");
        assert_eq!(cols[37], "3");
        assert_eq!(cols[38], "1");
        assert_eq!(cols[39], "1");

        // Verify debug column
        let mut av_dbg = av.clone();
        av_dbg.variant.debug = "C:4:F-2:R-2:1.0000".to_string();
        let dbg_line = amplicon_40(&av_dbg, &region(), "SAMPLE", true);
        assert!(dbg_line.starts_with(&line));
        assert!(dbg_line.contains("\tC:4:F-2:R-2:1.0000"));
        assert_eq!(dbg_line.split('\t').count(), 41);
    }

    #[test]
    fn somatic_55_pins_column_count_and_boundary_fields() {
        let normal = variant();
        let line = somatic_55(
            &variant(),
            Some(&normal),
            &region(),
            "SAMPLE",
            "StrongSomatic",
            "12",
            "",
            true,
        );
        assert_eq!(line.split('\t').count(), 55);
        assert!(line.starts_with(PREFIX_7), "prefix drifted: {line}");
        // Tail: status, var_type, duprate1, sv1, duprate2, sv2("" -> "0").
        assert!(
            line.ends_with("\tStrongSomatic\tSNV\t0\t12\t0\t0"),
            "suffix drifted: {line}"
        );
    }

    #[test]
    fn somatic_61_pins_column_count_and_boundary_fields() {
        let normal = variant();
        let line = somatic_61(
            &variant(),
            Some(&normal),
            &region(),
            "SAMPLE",
            "StrongSomatic",
            "",
            "34",
            true,
        );
        assert_eq!(line.split('\t').count(), 61);
        assert!(line.starts_with(PREFIX_7), "prefix drifted: {line}");
        // Fields 54-59: status, var_type, duprate1, sv1("" -> "0"), duprate2,
        // sv2; 60-61 are the total fisher pvalue/opposite-ratio stats (pinned
        // by the stats tests, not here).
        let tail: Vec<&str> = line.split('\t').skip(53).take(6).collect();
        assert_eq!(
            tail,
            vec!["StrongSomatic", "SNV", "0", "0", "0", "34"],
            "suffix drifted: {line}"
        );
    }
}

#[cfg(test)]
mod debug_column_tests {
    use super::*;
    use crate::postprocess::amplicon::AmpliconVariant;
    use crate::regions::Region;
    use crate::types::Pos1;

    fn region() -> Region {
        Region {
            chr: "chr1".to_string(),
            start: Pos1(933866),
            end: Pos1(934466),
            gene: "testbed".to_string(),
            amplicon: None,
        }
    }

    fn amplicon_variant() -> AmpliconVariant {
        let mut v = Variant::new();
        v.ref_allele = "G".to_string();
        v.var_allele = "C".to_string();
        v.genotype = None; // null genotype must print as "0" in the debug column
        v.frequency = 0.0; // zero freq must print as "0"
        v.total_pos_cov = 4;
        v.pos_coverage = 4;
        AmpliconVariant {
            variant: v.clone(),
            gdamp: 1,
            tlamp: 1,
            ncamp: 0,
            ampflag: false,
            region: "chr1:933867-934466".to_string(),
            good_variants: vec![(Some(v), "chr1:933867-934466".to_string())],
            bad_variants: vec![(None, "chr1:100-200".to_string())],
        }
    }

    #[test]
    fn debug_amp_variant_field_formats() {
        let mut v = Variant::new();
        v.total_pos_cov = 8;
        v.pos_coverage = 3;
        v.fwd = 1;
        v.rev = 2;
        v.genotype = Some("G/+38".to_string());
        v.frequency = 0.375;
        v.bias = "0;2".to_string();
        v.mean_position = 35.25;
        v.pstd = true;
        v.qstd = true;
        v.mean_quality = 0.0; // always formatted → "0.0"
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 6.0;
        v.hi_freq = 0.0; // zero → "0"
        v.extra_freq = 0.375;
        let s = super::debug_amp_variant(Some(&v));
        assert_eq!(
            s,
            "8 3 0 0 1 2 G/+38 0.3750 0;2 35.2 1 0.0 1 60.0 6.000 0 0.3750"
        );
    }

    #[test]
    fn debug_amp_variant_null_renders_zeros() {
        // Java NPEs here (AmpliconOutputVariant.java:258); we render an all-zero record.
        assert_eq!(
            super::debug_amp_variant(None),
            "0 0 0 0 0 0 0 0 0 0.0 0 0.0 0 0.0 0.000 0 0"
        );
    }

    #[test]
    fn amplicon_38_debug_column_only_when_enabled() {
        let mut av = amplicon_variant();
        let rg = region();

        let plain = super::amplicon_38(&av, &rg, "SAMPLE", false);
        assert!(!plain.contains("Good0"));
        assert_eq!(plain.split('\t').count(), 38);

        av.variant.debug = "C:4:F-2:R-2:1.0000".to_string();
        let dbg = super::amplicon_38(&av, &rg, "SAMPLE", true);
        assert!(dbg.contains("\tC:4:F-2:R-2:1.0000\tGood0 "));
        assert!(dbg.contains("\tBad0 0 0 0 0 0 0 0 0 0 0.0 0 0.0 0 0.0 0.000 0 0 chr1:100-200"));
        // genotype null → "0", zero freq → "0"
        assert!(
            dbg.contains("Good0 4 4 0 0 0 0 0 0 0 0.0 0 0.0 0 0.0 0.000 0 0 chr1:933867-934466")
        );
        // 38 base columns + DEBUG column + one tab-separated Good/Bad entry each
        assert_eq!(dbg.split('\t').count(), 41);
    }
}

#[cfg(test)]
mod header_tests {
    use super::print_header;
    use crate::config::Config;

    fn render(config: &Config) -> String {
        let mut buf: Vec<u8> = Vec::new();
        print_header(&mut buf, config).unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn simple_header_matches_java() {
        // ports SimpleMode.java:116-123 (36 columns)
        let cfg = Config::default();
        let h = render(&cfg);
        assert_eq!(h.split('\t').count(), 36);
        assert!(h.starts_with("Sample\tGene\tChr\tStart\tEnd\tRef\tAlt\t"));
        assert!(h.ends_with("Seg\tVarType\tDuprate\tSV_info\n"));
        assert!(!h.contains("CRISPR"));
        // No somatic-only or amplicon-only tokens
        assert!(!h.contains("VarLabel"));
        assert!(!h.contains("Ampflag"));
    }

    #[test]
    fn crispr_appends_column() {
        // ports SimpleMode.java:124-126
        let mut cfg = Config::default();
        cfg.crispr_cutting_site = Some(3);
        let h = render(&cfg);
        assert_eq!(h.split('\t').count(), 37);
        assert!(h.ends_with("\tCRISPR\n"));
    }

    #[test]
    fn somatic_header_matches_java() {
        // ports SomaticMode.java:137-145 (55 columns, repeated block)
        let mut cfg = Config::default();
        cfg.bam = "tumor.bam|normal.bam".to_string();
        let h = render(&cfg);
        let java = [
            "Sample",
            "Gene",
            "Chr",
            "Start",
            "End",
            "Ref",
            "Alt",
            "Depth",
            "AltDepth",
            "RefFwdReads",
            "RefRevReads",
            "AltFwdReads",
            "AltRevReads",
            "Genotype",
            "AF",
            "Bias",
            "PMean",
            "PStd",
            "QMean",
            "QStd",
            "MQ",
            "Sig_Noise",
            "HiAF",
            "ExtraAF",
            "NM",
            "Depth",
            "AltDepth",
            "RefFwdReads",
            "RefRevReads",
            "AltFwdReads",
            "AltRevReads",
            "Genotype",
            "AF",
            "Bias",
            "PMean",
            "PStd",
            "QMean",
            "QStd",
            "MQ",
            "Sig_Noise",
            "HiAF",
            "ExtraAF",
            "NM",
            "shift3",
            "MSI",
            "MSI_NT",
            "5pFlankSeq",
            "3pFlankSeq",
            "Seg",
            "VarLabel",
            "VarType",
            "Duprate1",
            "SV_info1",
            "Duprate2",
            "SV_info2",
        ];
        assert_eq!(h.trim_end(), java.join("\t"));
        assert_eq!(h.split('\t').count(), 55);
        assert!(h.ends_with('\n'));
    }

    #[test]
    fn amplicon_header_matches_java() {
        // ports AmpliconMode.java:136-141 (38 columns)
        let mut cfg = Config::default();
        cfg.amplicon = Some("10:0.95".to_string());
        let h = render(&cfg);
        assert_eq!(h.split('\t').count(), 38);
        assert!(h.ends_with("GoodVarCount\tTotalVarCount\tNocov\tAmpflag\n"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::postprocess::amplicon::AmpliconVariant;
    use crate::regions::Region;
    use crate::types::Pos1;

    #[test]
    fn table_formatters_use_half_even_rounding() {
        assert_eq!(fmt_3dp_msi(0.0), "0");
        assert_eq!(fmt_3dp_msi(1.125), "1.125");
        assert_eq!(fmt_3dp_msi(1.0005), "1.000");

        let mut v = Variant::new();
        v.ref_allele = "A".to_string();
        v.var_allele = "G".to_string();
        v.start_pos = 100;
        v.end_pos = 100;
        v.mean_position = 10.25; // -> 10.2 (half-even)
        v.mean_quality = 30.35; // -> 30.4 (half-even)
        v.mean_mapping_quality = 40.25; // -> 40.2 (half-even)
        v.hq_lq_ratio = 2.125;
        v.hi_freq = 0.125;
        v.nm = 1.25; // -> 1.2 (half-even)
        v.msi = 2.0005;

        let dbg = debug_amp_variant(Some(&v));
        assert!(dbg.contains(" 10.2 0 30.4 0 40.2 2.125 0.1250 "), "{dbg}");

        let rg = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "G".to_string(),
            amplicon: None,
        };
        let s38 = simple_38(&v, &rg, "S", "");
        let cols38: Vec<&str> = s38.split('\t').collect();
        assert_eq!(cols38[24], "0.1250");
        assert_eq!(cols38[29], "1.2");

        let av = AmpliconVariant {
            variant: v,
            gdamp: 1,
            tlamp: 1,
            ncamp: 0,
            ampflag: false,
            region: "chr1:1-200".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };
        let a40 = amplicon_40(&av, &rg, "S", false);
        let cols40: Vec<&str> = a40.split('\t').collect();
        assert_eq!(cols40[24], "0.1250");
        assert_eq!(cols40[29], "1.2");
        assert_eq!(get_rounded_value_to_print(-0.00001, 4), "0");
    }
}
