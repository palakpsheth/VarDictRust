//! Table output serialization formats — intermediate .var output formatters.
//! Ports SimpleOutputVariant.java.

use crate::regions::Region;
use crate::variations::Variant;
use std::fmt::Write as _;

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

/// Write the simple 36-column row into `buf` (cleared first).
/// Byte-identical to [`simple_36`]; the caller's buffer is reused across
/// rows (one reserve, `clear()` per variant) instead of one temp `String`
/// per float plus a final `format!`.
/// Ports SimpleOutputVariant.java's create_simple_variant_36columns.
pub fn simple_36_into(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    sample_name: &str,
    sv: &str,
) {
    buf.clear();
    // Golden p50 row is 213 B, max 505 B: 256 covers the common row, long
    // alleles regrow once. Retained rows keep at most this slack.
    buf.reserve(256);
    buf.push_str(sample_name);
    buf.push('\t');
    buf.push_str(&region.gene);
    buf.push('\t');
    buf.push_str(&region.chr);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.start_pos);
    let _ = write!(buf, "{}\t", variant.end_pos);
    buf.push_str(&variant.ref_allele);
    buf.push('\t');
    buf.push_str(&variant.var_allele);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.total_pos_cov);
    let _ = write!(buf, "{}\t", variant.pos_coverage);
    let _ = write!(buf, "{}\t", variant.ref_fwd_cov);
    let _ = write!(buf, "{}\t", variant.ref_rev_cov);
    let _ = write!(buf, "{}\t", variant.fwd);
    let _ = write!(buf, "{}\t", variant.rev);
    // Null records render the printer field defaults, not the normal
    // null→"0" / empty→"0" mapping (SimpleOutputVariant.java:50-52 vs :64,
    // :81-82; OutputVariant.java:22-23).
    if variant.null_record {
        buf.push('\t');
    } else {
        buf.push_str(variant.genotype.as_deref().unwrap_or("0"));
        buf.push('\t');
    }
    push_keep(buf, variant.frequency, 4);
    buf.push('\t');
    buf.push_str(&variant.bias);
    buf.push('\t');
    push_keep(buf, variant.mean_position, 1);
    buf.push('\t');
    buf.push(if variant.pstd { '1' } else { '0' });
    buf.push('\t');
    push_keep(buf, variant.mean_quality, 1);
    buf.push('\t');
    buf.push(if variant.qstd { '1' } else { '0' });
    buf.push('\t');
    push_keep(buf, variant.mean_mapping_quality, 1);
    buf.push('\t');
    push_keep(buf, variant.hq_lq_ratio, 3);
    buf.push('\t');
    push_keep(buf, variant.hi_freq, 4);
    buf.push('\t');
    push_keep(buf, variant.extra_freq, 4);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.shift3);
    push_keep(buf, variant.msi, 3);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.msint);
    push_nm(buf, variant.nm);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.hicnt);
    let _ = write!(buf, "{}\t", variant.hicov);
    if variant.null_record {
        buf.push('\t');
    } else if variant.left_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.left_seq);
        buf.push('\t');
    }
    if variant.null_record {
        buf.push('\t');
    } else if variant.right_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.right_seq);
        buf.push('\t');
    }
    buf.push_str(&region.chr);
    buf.push(':');
    let _ = write!(buf, "{}-{}\t", region.start.0, region.end.0);
    buf.push_str(&variant.var_type);
    buf.push('\t');
    push_keep(buf, variant.duprate, 1);
    buf.push('\t');
    if sv.is_empty() {
        buf.push('0');
    } else {
        buf.push_str(sv);
    }
}

/// Write simple 36-column representation of finalized Variants to a tab-delimited string.
/// Ports SimpleOutputVariant.java's create_simple_variant_36columns.
pub fn simple_36(variant: &Variant, region: &Region, sample_name: &str, sv: &str) -> String {
    let mut buf = String::new();
    simple_36_into(&mut buf, variant, region, sample_name, sv);
    buf
}

/// Write the simple 38-column row into `buf` (cleared first).
/// Byte-identical to [`simple_38`]; same reuse contract as [`simple_36_into`].
/// Ports SimpleOutputVariant.java's create_simple_variant_38columns.
pub fn simple_38_into(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    sample_name: &str,
    sv: &str,
) {
    buf.clear();
    buf.reserve(320);
    let fisher = crate::stats::FisherExact::new_java(
        variant.ref_fwd_cov as i32,
        variant.ref_rev_cov as i32,
        variant.fwd as i32,
        variant.rev as i32,
    );
    let pvalue = fisher.get_p_value();
    let oddratio = fisher.get_odds_ratio();

    buf.push_str(sample_name);
    buf.push('\t');
    buf.push_str(&region.gene);
    buf.push('\t');
    buf.push_str(&region.chr);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.start_pos);
    let _ = write!(buf, "{}\t", variant.end_pos);
    buf.push_str(&variant.ref_allele);
    buf.push('\t');
    buf.push_str(&variant.var_allele);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.total_pos_cov);
    let _ = write!(buf, "{}\t", variant.pos_coverage);
    let _ = write!(buf, "{}\t", variant.ref_fwd_cov);
    let _ = write!(buf, "{}\t", variant.ref_rev_cov);
    let _ = write!(buf, "{}\t", variant.fwd);
    let _ = write!(buf, "{}\t", variant.rev);
    // Same null-record printer-defaults rule as simple_36.
    if variant.null_record {
        buf.push('\t');
    } else {
        buf.push_str(variant.genotype.as_deref().unwrap_or("0"));
        buf.push('\t');
    }
    push_stripped(buf, variant.frequency, 4);
    buf.push('\t');
    buf.push_str(&variant.bias);
    buf.push('\t');
    push_stripped(buf, variant.mean_position, 1);
    buf.push('\t');
    buf.push(if variant.pstd { '1' } else { '0' });
    buf.push('\t');
    push_stripped(buf, variant.mean_quality, 1);
    buf.push('\t');
    buf.push(if variant.qstd { '1' } else { '0' });
    buf.push('\t');
    push_stripped(buf, pvalue, 5);
    buf.push('\t');
    buf.push_str(&oddratio);
    buf.push('\t');
    push_stripped(buf, variant.mean_mapping_quality, 1);
    buf.push('\t');
    push_stripped(buf, variant.hq_lq_ratio, 3);
    buf.push('\t');
    // hifreq/nm cells keep trailing zeros in BOTH modes (raw-cell rule).
    push_keep(buf, variant.hi_freq, 4);
    buf.push('\t');
    push_stripped(buf, variant.extra_freq, 4);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.shift3);
    push_stripped(buf, variant.msi, 3);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.msint);
    push_nm(buf, variant.nm);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.hicnt);
    let _ = write!(buf, "{}\t", variant.hicov);
    if variant.null_record {
        buf.push('\t');
    } else if variant.left_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.left_seq);
        buf.push('\t');
    }
    if variant.null_record {
        buf.push('\t');
    } else if variant.right_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.right_seq);
        buf.push('\t');
    }
    buf.push_str(&region.chr);
    buf.push(':');
    let _ = write!(buf, "{}-{}\t", region.start.0, region.end.0);
    buf.push_str(&variant.var_type);
    buf.push('\t');
    push_stripped(buf, variant.duprate, 2);
    buf.push('\t');
    if sv.is_empty() {
        buf.push('0');
    } else {
        buf.push_str(sv);
    }
}

/// Write simple 38-column representation of finalized Variants including Fisher Exact test results.
/// Ports SimpleOutputVariant.java's create_simple_variant_38columns.
pub fn simple_38(variant: &Variant, region: &Region, sample_name: &str, sv: &str) -> String {
    let mut buf = String::new();
    simple_38_into(&mut buf, variant, region, sample_name, sv);
    buf
}

/// Write the amplicon 38-column row into `buf` (cleared first).
/// Byte-identical to [`amplicon_38`]; same reuse contract as [`simple_36_into`].
/// Ports AmpliconOutputVariant.java's create_amplicon_variant_38columns.
/// Columns 1-34 match the simple format; columns 35-38 are
/// GDAMP, TLAMP, NCAMP, AMPFLAG (replacing simple's duprate/sv).
/// When `debug` is set, appends `"\t" + DEBUG` where
/// DEBUG = variant.DEBUG + goodOutput + badOutput (AmpliconOutputVariant.java:113-138).
pub fn amplicon_38_into(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) {
    buf.clear();
    buf.reserve(256);
    let variant = &av.variant;

    buf.push_str(sample_name);
    buf.push('\t');
    buf.push_str(&region.gene);
    buf.push('\t');
    buf.push_str(&region.chr);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.start_pos);
    let _ = write!(buf, "{}\t", variant.end_pos);
    buf.push_str(&variant.ref_allele);
    buf.push('\t');
    buf.push_str(&variant.var_allele);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.total_pos_cov);
    let _ = write!(buf, "{}\t", variant.pos_coverage);
    let _ = write!(buf, "{}\t", variant.ref_fwd_cov);
    let _ = write!(buf, "{}\t", variant.ref_rev_cov);
    let _ = write!(buf, "{}\t", variant.fwd);
    let _ = write!(buf, "{}\t", variant.rev);
    buf.push_str(variant.genotype.as_deref().unwrap_or("0"));
    buf.push('\t');
    push_keep(buf, variant.frequency, 4);
    buf.push('\t');
    buf.push_str(&variant.bias);
    buf.push('\t');
    push_keep(buf, variant.mean_position, 1);
    buf.push('\t');
    buf.push(if variant.pstd { '1' } else { '0' });
    buf.push('\t');
    push_keep(buf, variant.mean_quality, 1);
    buf.push('\t');
    buf.push(if variant.qstd { '1' } else { '0' });
    buf.push('\t');
    push_keep(buf, variant.mean_mapping_quality, 1);
    buf.push('\t');
    push_keep(buf, variant.hq_lq_ratio, 3);
    buf.push('\t');
    push_keep(buf, variant.hi_freq, 4);
    buf.push('\t');
    push_keep(buf, variant.extra_freq, 4);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.shift3);
    push_keep(buf, variant.msi, 3);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.msint);
    push_nm(buf, variant.nm);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.hicnt);
    let _ = write!(buf, "{}\t", variant.hicov);
    if variant.left_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.left_seq);
        buf.push('\t');
    }
    if variant.right_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.right_seq);
        buf.push('\t');
    }
    buf.push_str(&av.region);
    buf.push('\t');
    buf.push_str(&variant.var_type);
    buf.push('\t');
    let _ = write!(
        buf,
        "{}\t{}\t{}\t{}",
        av.gdamp,
        av.tlamp,
        av.ncamp,
        if av.ampflag { 1 } else { 0 }
    );

    if !debug {
        return;
    }
    // ports AmpliconOutputVariant.java:113-125
    buf.push('\t');
    buf.push_str(&variant.debug);
    for (i, (v, r)) in av.good_variants.iter().enumerate() {
        let _ = write!(buf, "\tGood{} {} {}", i, debug_amp_variant(v.as_ref()), r);
    }
    for (i, (v, r)) in av.bad_variants.iter().enumerate() {
        let _ = write!(buf, "\tBad{} {} {}", i, debug_amp_variant(v.as_ref()), r);
    }
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
    let mut buf = String::new();
    amplicon_38_into(&mut buf, av, region, sample_name, debug);
    buf
}

/// Write the amplicon 40-column row into `buf` (cleared first).
/// Byte-identical to [`amplicon_40`]; same reuse contract as [`simple_36_into`].
/// Ports AmpliconOutputVariant.java's create_amplicon_variant_40columns (lines 145-201).
/// Columns 1-20 match amplicon_38, columns 21-22 are pvalue and oddratio,
/// columns 23-36 match columns 21-34 of amplicon_38, and columns 37-40 are
/// GDAMP, TLAMP, NCAMP, AMPFLAG.
/// When `debug` is set, appends `"\t" + DEBUG` where
/// DEBUG = variant.DEBUG + goodOutput + badOutput.
pub fn amplicon_40_into(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) {
    buf.clear();
    buf.reserve(320);
    let variant = &av.variant;

    let fisher = crate::stats::FisherExact::new_java(
        variant.ref_fwd_cov as i32,
        variant.ref_rev_cov as i32,
        variant.fwd as i32,
        variant.rev as i32,
    );
    let pvalue = fisher.get_p_value();
    let oddratio = fisher.get_odds_ratio();

    buf.push_str(sample_name);
    buf.push('\t');
    buf.push_str(&region.gene);
    buf.push('\t');
    buf.push_str(&region.chr);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.start_pos);
    let _ = write!(buf, "{}\t", variant.end_pos);
    buf.push_str(&variant.ref_allele);
    buf.push('\t');
    buf.push_str(&variant.var_allele);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.total_pos_cov);
    let _ = write!(buf, "{}\t", variant.pos_coverage);
    let _ = write!(buf, "{}\t", variant.ref_fwd_cov);
    let _ = write!(buf, "{}\t", variant.ref_rev_cov);
    let _ = write!(buf, "{}\t", variant.fwd);
    let _ = write!(buf, "{}\t", variant.rev);
    buf.push_str(variant.genotype.as_deref().unwrap_or("0"));
    buf.push('\t');
    push_stripped(buf, variant.frequency, 4);
    buf.push('\t');
    buf.push_str(&variant.bias);
    buf.push('\t');
    push_stripped(buf, variant.mean_position, 1);
    buf.push('\t');
    buf.push(if variant.pstd { '1' } else { '0' });
    buf.push('\t');
    push_stripped(buf, variant.mean_quality, 1);
    buf.push('\t');
    buf.push(if variant.qstd { '1' } else { '0' });
    buf.push('\t');
    push_stripped(buf, pvalue, 5);
    buf.push('\t');
    buf.push_str(&oddratio);
    buf.push('\t');
    push_stripped(buf, variant.mean_mapping_quality, 1);
    buf.push('\t');
    push_stripped(buf, variant.hq_lq_ratio, 3);
    buf.push('\t');
    push_keep(buf, variant.hi_freq, 4);
    buf.push('\t');
    push_stripped(buf, variant.extra_freq, 4);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.shift3);
    push_stripped(buf, variant.msi, 3);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.msint);
    push_nm(buf, variant.nm);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.hicnt);
    let _ = write!(buf, "{}\t", variant.hicov);
    if variant.left_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.left_seq);
        buf.push('\t');
    }
    if variant.right_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&variant.right_seq);
        buf.push('\t');
    }
    buf.push_str(&av.region);
    buf.push('\t');
    buf.push_str(&variant.var_type);
    buf.push('\t');
    let _ = write!(
        buf,
        "{}\t{}\t{}\t{}",
        av.gdamp,
        av.tlamp,
        av.ncamp,
        if av.ampflag { 1 } else { 0 }
    );

    if !debug {
        return;
    }
    // ports AmpliconOutputVariant.java:113-125
    buf.push('\t');
    buf.push_str(&variant.debug);
    for (i, (v, r)) in av.good_variants.iter().enumerate() {
        let _ = write!(buf, "\tGood{} {} {}", i, debug_amp_variant(v.as_ref()), r);
    }
    for (i, (v, r)) in av.bad_variants.iter().enumerate() {
        let _ = write!(buf, "\tBad{} {} {}", i, debug_amp_variant(v.as_ref()), r);
    }
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
    let mut buf = String::new();
    amplicon_40_into(&mut buf, av, region, sample_name, debug);
    buf
}

/// One space-joined entry of the -D amplicon debug column.
/// Ports `AmpliconOutputVariant.debugAmpVariant` (AmpliconOutputVariant.java:258-278):
/// 17 fields; zero-valued freq fields print as "0", position/quality/mapq/qratio are
/// always decimal-formatted. Java dereferences a null variant here (NPE on the
/// no-data bad-variant path); we render an all-zero record instead.
/// Debug-only path: keeps its own temp-`String` keep-mode helper.
fn fmt_4dp(v: f64) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!("{:.4}", crate::types::round_half_even(v, 4))
    }
}

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

/// Java-compatible formatter that rounds a float using banker's rounding
/// to a maximum of `digits` decimal places, and strips any trailing zeros.
/// Ports `Utils.getRoundedValueToPrint`. Shared with the VCF writer (vcf.rs).
/// Same value mapping as [`push_stripped`]; kept for the VCF rare paths and
/// the formatter-mode pins below.
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

/// Keep-mode float cell: HALF_EVEN at `digits` decimals with trailing zeros
/// kept (`"0"` when zero). Ports the `DecimalFormat("0.0000"/"0.0"/"0.000")`
/// cells of the no-`--fisher` writers (Simple 36-col, Amplicon 38-col,
/// Somatic 55-col) — one of the two `round_half_even` call sites. Pushes
/// into the row buffer without a temp `String`.
pub(crate) fn push_keep(buf: &mut String, v: f64, digits: i32) {
    if v == 0.0 {
        buf.push('0');
    } else {
        let _ = write!(
            buf,
            "{:.*}",
            digits as usize,
            crate::types::round_half_even(v, digits)
        );
    }
}

/// Strip-mode float cell: HALF_EVEN at `digits` decimals with trailing zeros
/// stripped (`"0"` when zero, integers without a dot). Same value mapping as
/// [`get_rounded_value_to_print`] — the other `round_half_even` call site —
/// but pushed into the row buffer in place (no temp `String`): trailing zeros
/// are truncated off the buffer tail.
/// Ports `Utils.getRoundedValueToPrint` (the `--fisher` writers).
pub(crate) fn push_stripped(buf: &mut String, value: f64, digits: i32) {
    if value.is_nan() || value.is_infinite() {
        let _ = write!(buf, "{}", value);
        return;
    }
    let rounded = crate::types::round_half_even(value, digits);
    if rounded == 0.0 {
        buf.push('0');
    } else if rounded == rounded.round() {
        let _ = write!(buf, "{:.0}", rounded);
    } else {
        let start = buf.len();
        let _ = write!(buf, "{:.*}", digits as usize, rounded);
        while buf.len() > start && buf.as_bytes().get(buf.len() - 1) == Some(&b'0') {
            buf.pop();
        }
    }
}

/// NM cell: keep-mode 1dp when positive, `"0"` otherwise (matches `fmt_nm`).
pub(crate) fn push_nm(buf: &mut String, v: f64) {
    if v > 0.0 {
        push_keep(buf, v, 1);
    } else {
        buf.push('0');
    }
}

/// Push the 18-column tumor/normal sample block (keep-mode floats).
/// Replaces `variant_18cols`; same field order and values, no temp `String`.
fn push_18cols(buf: &mut String, v: &Variant) {
    let _ = write!(buf, "{}\t", v.total_pos_cov);
    let _ = write!(buf, "{}\t", v.pos_coverage);
    let _ = write!(buf, "{}\t", v.ref_fwd_cov);
    let _ = write!(buf, "{}\t", v.ref_rev_cov);
    let _ = write!(buf, "{}\t", v.fwd);
    let _ = write!(buf, "{}\t", v.rev);
    buf.push_str(v.genotype.as_deref().unwrap_or("0"));
    buf.push('\t');
    push_keep(buf, v.frequency, 4);
    buf.push('\t');
    buf.push_str(&v.bias);
    buf.push('\t');
    push_keep(buf, v.mean_position, 1);
    buf.push('\t');
    buf.push(if v.pstd { '1' } else { '0' });
    buf.push('\t');
    push_keep(buf, v.mean_quality, 1);
    buf.push('\t');
    buf.push(if v.qstd { '1' } else { '0' });
    buf.push('\t');
    push_keep(buf, v.mean_mapping_quality, 1);
    buf.push('\t');
    push_keep(buf, v.hq_lq_ratio, 3);
    buf.push('\t');
    push_keep(buf, v.hi_freq, 4);
    buf.push('\t');
    push_keep(buf, v.extra_freq, 4);
    buf.push('\t');
    push_nm(buf, v.nm);
}

/// Missing-normal 18-column block: eighteen `"0"` cells (was
/// `variant_18cols_empty`). Inlined so the common true-somatic row (no normal
/// allele) allocates nothing here.
const VARIANT_18COLS_EMPTY: &str = "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";

/// Write the somatic 55-column row into `buf` (cleared first).
/// Byte-identical to [`somatic_55`]; same reuse contract as [`simple_36_into`].
/// Ports SomaticOutputVariant.java's create_somatic_variant_55columns.
pub fn somatic_55_into(
    buf: &mut String,
    tumor: &Variant,
    normal: Option<&Variant>,
    region: &Region,
    sample_name: &str,
    status: &str,
    sv1: &str,
    sv2: &str,
    end_is_variant2: bool,
) {
    buf.clear();
    buf.reserve(384);
    // endVariant supplies shift3/msi/msint/left/right seq (Java SomaticOutputVariant).
    let end_var = if end_is_variant2 {
        normal.unwrap_or(tumor)
    } else {
        tumor
    };

    buf.push_str(sample_name);
    buf.push('\t');
    buf.push_str(&region.gene);
    buf.push('\t');
    buf.push_str(&region.chr);
    buf.push('\t');
    let _ = write!(buf, "{}\t", tumor.start_pos);
    let _ = write!(buf, "{}\t", tumor.end_pos);
    buf.push_str(&tumor.ref_allele);
    buf.push('\t');
    buf.push_str(&tumor.var_allele);
    buf.push('\t');
    push_18cols(buf, tumor);
    buf.push('\t');
    match normal {
        Some(n) => push_18cols(buf, n),
        None => buf.push_str(VARIANT_18COLS_EMPTY),
    }
    buf.push('\t');
    let _ = write!(buf, "{}\t", end_var.shift3);
    push_keep(buf, end_var.msi, 3);
    buf.push('\t');
    let _ = write!(buf, "{}\t", end_var.msint);
    if end_var.left_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&end_var.left_seq);
        buf.push('\t');
    }
    if end_var.right_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&end_var.right_seq);
        buf.push('\t');
    }
    buf.push_str(&region.chr);
    buf.push(':');
    let _ = write!(buf, "{}-{}\t", region.start.0, region.end.0);
    buf.push_str(status);
    buf.push('\t');
    buf.push_str(&tumor.var_type);
    buf.push('\t');
    push_keep(buf, tumor.duprate, 1);
    buf.push('\t');
    if sv1.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(sv1);
        buf.push('\t');
    }
    match normal {
        Some(n) => push_keep(buf, n.duprate, 1),
        None => buf.push('0'),
    }
    buf.push('\t');
    if sv2.is_empty() {
        buf.push('0');
    } else {
        buf.push_str(sv2);
    }
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
    let mut buf = String::new();
    somatic_55_into(
        &mut buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
    );
    buf
}

/// Push the 20-column fisher sample block: the 18 base columns plus pvalue
/// and oddratio, using trailing-zero-stripping formatting (Java
/// getRoundedValueToPrint). Replaces `variant_20cols_fisher` with no temp.
fn push_20cols_fisher(buf: &mut String, v: &Variant, pvalue: f64, oddratio: &str) {
    let _ = write!(buf, "{}\t", v.total_pos_cov);
    let _ = write!(buf, "{}\t", v.pos_coverage);
    let _ = write!(buf, "{}\t", v.ref_fwd_cov);
    let _ = write!(buf, "{}\t", v.ref_rev_cov);
    let _ = write!(buf, "{}\t", v.fwd);
    let _ = write!(buf, "{}\t", v.rev);
    buf.push_str(v.genotype.as_deref().unwrap_or("0"));
    buf.push('\t');
    push_stripped(buf, v.frequency, 4);
    buf.push('\t');
    buf.push_str(&v.bias);
    buf.push('\t');
    push_stripped(buf, v.mean_position, 1);
    buf.push('\t');
    buf.push(if v.pstd { '1' } else { '0' });
    buf.push('\t');
    push_stripped(buf, v.mean_quality, 1);
    buf.push('\t');
    buf.push(if v.qstd { '1' } else { '0' });
    buf.push('\t');
    push_stripped(buf, v.mean_mapping_quality, 1);
    buf.push('\t');
    push_stripped(buf, v.hq_lq_ratio, 3);
    buf.push('\t');
    push_stripped(buf, v.hi_freq, 4);
    buf.push('\t');
    push_stripped(buf, v.extra_freq, 4);
    buf.push('\t');
    if v.nm > 0.0 {
        push_stripped(buf, v.nm, 1);
    } else {
        buf.push('0');
    }
    buf.push('\t');
    push_stripped(buf, pvalue, 5);
    buf.push('\t');
    buf.push_str(oddratio);
}

fn variant_20cols_fisher_empty(pvalue: f64, oddratio: &str) -> String {
    format!(
        "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t{}\t{}",
        get_rounded_value_to_print(pvalue, 5),
        oddratio,
    )
}

/// Write the somatic 61-column row into `buf` (cleared first).
/// Byte-identical to [`somatic_61`]; same reuse contract as [`simple_36_into`].
/// Ports SomaticOutputVariant.java's create_somatic_variant_61columns.
pub fn somatic_61_into(
    buf: &mut String,
    tumor: &Variant,
    normal: Option<&Variant>,
    region: &Region,
    sample_name: &str,
    status: &str,
    sv1: &str,
    sv2: &str,
    end_is_variant2: bool,
) {
    buf.clear();
    buf.reserve(448);
    // endVariant supplies shift3/msi/msint/left/right seq (Java SomaticOutputVariant).
    let end_var = if end_is_variant2 {
        normal.unwrap_or(tumor)
    } else {
        tumor
    };

    let fisher1 = crate::stats::FisherExact::new_java(
        tumor.ref_fwd_cov as i32,
        tumor.ref_rev_cov as i32,
        tumor.fwd as i32,
        tumor.rev as i32,
    );
    let pvalue1 = fisher1.get_p_value();
    let oddratio1 = fisher1.get_odds_ratio();

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

    buf.push_str(sample_name);
    buf.push('\t');
    buf.push_str(&region.gene);
    buf.push('\t');
    buf.push_str(&region.chr);
    buf.push('\t');
    let _ = write!(buf, "{}\t", tumor.start_pos);
    let _ = write!(buf, "{}\t", tumor.end_pos);
    buf.push_str(&tumor.ref_allele);
    buf.push('\t');
    buf.push_str(&tumor.var_allele);
    buf.push('\t');
    push_20cols_fisher(buf, tumor, pvalue1, &oddratio1);
    buf.push('\t');
    match normal {
        Some(n) => {
            let f2 = crate::stats::FisherExact::new_java(
                n.ref_fwd_cov as i32,
                n.ref_rev_cov as i32,
                n.fwd as i32,
                n.rev as i32,
            );
            push_20cols_fisher(buf, n, f2.get_p_value(), &f2.get_odds_ratio());
        }
        None => {
            let f2 = crate::stats::FisherExact::new_java(0, 0, 0, 0);
            buf.push_str(&variant_20cols_fisher_empty(
                f2.get_p_value(),
                &f2.get_odds_ratio(),
            ));
        }
    }
    buf.push('\t');
    let _ = write!(buf, "{}\t", end_var.shift3);
    push_keep(buf, end_var.msi, 3);
    buf.push('\t');
    let _ = write!(buf, "{}\t", end_var.msint);
    if end_var.left_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&end_var.left_seq);
        buf.push('\t');
    }
    if end_var.right_seq.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(&end_var.right_seq);
        buf.push('\t');
    }
    buf.push_str(&region.chr);
    buf.push(':');
    let _ = write!(buf, "{}-{}\t", region.start.0, region.end.0);
    buf.push_str(status);
    buf.push('\t');
    buf.push_str(&tumor.var_type);
    buf.push('\t');
    push_stripped(buf, tumor.duprate, 1);
    buf.push('\t');
    if sv1.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(sv1);
        buf.push('\t');
    }
    match normal {
        Some(n) => push_stripped(buf, n.duprate, 1),
        None => buf.push('0'),
    }
    buf.push('\t');
    if sv2.is_empty() {
        buf.push_str("0\t");
    } else {
        buf.push_str(sv2);
        buf.push('\t');
    }
    push_stripped(buf, total_pvalue, 5);
    buf.push('\t');
    buf.push_str(&total_oddratio);
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
    let mut buf = String::new();
    somatic_61_into(
        &mut buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
    );
    buf
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
mod buffer_reuse_tests {
    // Task 6 pins (behavior only, no alloc/wall asserts):
    // - keep/strip formatter split across simple_38 / amplicon_40 / somatic_61
    // - null-record GENO "" path for 38-col (mirrors the 36-col pin above)
    // - `_into` buffer reuse: byte-equal to one-shot fns, `clear()` per row so
    //   a long row followed by a short row leaves no tail (never elide rows).
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

    fn rich_variant() -> Variant {
        let mut v = Variant::new();
        v.start_pos = 100;
        v.end_pos = 100;
        v.ref_allele = "G".to_string();
        v.var_allele = "C".to_string();
        v.var_type = "SNV".to_string();
        v.total_pos_cov = 40;
        v.pos_coverage = 15;
        v.ref_fwd_cov = 10;
        v.ref_rev_cov = 15;
        v.fwd = 8;
        v.rev = 7;
        v.genotype = Some("G/C".to_string());
        v.frequency = 0.375;
        v.bias = "1;2".to_string();
        v.mean_position = 20.0;
        v.pstd = true;
        v.mean_quality = 30.0;
        v.qstd = true;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 2.5;
        v.hi_freq = 0.375;
        v.extra_freq = 0.125;
        v.shift3 = 0;
        v.msi = 0.5;
        v.msint = 1;
        v.nm = 1.2;
        v.hicnt = 15;
        v.hicov = 40;
        v.left_seq = "ACGT".to_string();
        v.right_seq = "TGCA".to_string();
        v.duprate = 0.25;
        v
    }

    fn null_variant() -> Variant {
        let mut v = Variant::new();
        v.start_pos = 6447;
        v.end_pos = 6447;
        v.bias = "0;0".to_string();
        v.null_record = true;
        v
    }

    #[test]
    fn simple_38_null_record_renders_printer_defaults() {
        // 38-col null path: genotype "" (col 14), bias "0;0" (col 16),
        // flank seqs "" (cols 33/34) — same printer defaults as 36-col.
        let line = simple_38(&null_variant(), &region(), "SAMPLE", "");
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 38);
        assert_eq!(cols[3], "6447");
        assert_eq!(cols[4], "6447");
        assert_eq!(cols[13], "");
        assert_eq!(cols[15], "0;0");
        assert_eq!(cols[32], "");
        assert_eq!(cols[33], "");
    }

    #[test]
    fn simple_keep_vs_strip_split() {
        // Keep mode (36-col, DecimalFormat-style) retains trailing zeros;
        // strip mode (38-col, getRoundedValueToPrint) trims them — except
        // hifreq/nm which stay kept in BOTH modes (Java raw-cell rule).
        let v = rich_variant();
        let l36 = simple_36(&v, &region(), "S", "");
        let c36: Vec<&str> = l36.split('\t').collect();
        assert_eq!(c36[14], "0.3750");
        assert_eq!(c36[21], "2.500");
        assert_eq!(c36[23], "0.1250");
        assert_eq!(c36[25], "0.500");
        assert_eq!(c36[34], "0.2");
        let l38 = simple_38(&v, &region(), "S", "");
        let c38: Vec<&str> = l38.split('\t').collect();
        assert_eq!(c38[14], "0.375");
        assert_eq!(c38[23], "2.5");
        assert_eq!(c38[24], "0.3750");
        assert_eq!(c38[25], "0.125");
        assert_eq!(c38[27], "0.5");
        assert_eq!(c38[29], "1.2");
        assert_eq!(c38[36], "0.25");
    }

    #[test]
    fn amplicon_40_keep_vs_strip_split() {
        let av = AmpliconVariant {
            variant: rich_variant(),
            gdamp: 2,
            tlamp: 3,
            ncamp: 1,
            ampflag: true,
            region: "chr1:933867-934466".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };
        let line40 = amplicon_40(&av, &region(), "S", false);
        let cols: Vec<&str> = line40.split('\t').collect();
        assert_eq!(cols.len(), 40);
        assert_eq!(cols[14], "0.375");
        assert_eq!(cols[24], "0.3750");
        assert_eq!(cols[29], "1.2");
        assert_eq!(cols[36], "2");
        assert_eq!(cols[39], "1");
    }

    #[test]
    fn somatic_61_none_genotype_renders_zero_not_empty() {
        // Somatic/amplicon have no null_record branch: a missing genotype
        // renders "0", NOT the simple-null "". Both 20-col blocks pinned.
        let mut tumor = rich_variant();
        tumor.genotype = None;
        let mut normal = rich_variant();
        normal.genotype = None;
        let line = somatic_61(
            &tumor,
            Some(&normal),
            &region(),
            "SAMPLE",
            "StrongSomatic",
            "",
            "34",
            true,
        );
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 61);
        assert_eq!(cols[13], "0");
        assert_eq!(cols[33], "0");
    }

    #[test]
    fn strip_and_keep_modes_share_half_even() {
        // Both formatter modes round HALF_EVEN; they differ only in whether
        // trailing zeros are kept (no-fisher) or stripped (--fisher).
        assert_eq!(get_rounded_value_to_print(0.375, 4), "0.375");
        assert_eq!(get_rounded_value_to_print(1.0, 4), "1");
        assert_eq!(get_rounded_value_to_print(0.0, 4), "0");
        assert_eq!(get_rounded_value_to_print(-0.00001, 4), "0");
        assert_eq!(get_rounded_value_to_print(f64::NAN, 4), "NaN");
        let mut b = String::new();
        push_keep(&mut b, 0.375, 4);
        assert_eq!(b, "0.3750");
        b.clear();
        push_keep(&mut b, 0.0, 4);
        assert_eq!(b, "0");
        b.clear();
        push_stripped(&mut b, 0.375, 4);
        assert_eq!(b, "0.375");
        b.clear();
        push_stripped(&mut b, f64::NAN, 4);
        assert_eq!(b, "NaN");
    }

    #[test]
    fn into_reuse_matches_one_shot_and_never_elides() {
        // Long row then short null row into the SAME buffer: a missing
        // clear() would concatenate (dup) or drop (elide) rows.
        let rg = region();
        let mut buf = String::new();

        let mut long = rich_variant();
        long.ref_allele = "ACGTACGTACGT".to_string();
        long.var_allele = "TGCATGCATGCA".to_string();

        simple_36_into(&mut buf, &long, &rg, "SAMPLE", "1-2");
        assert_eq!(buf, simple_36(&long, &rg, "SAMPLE", "1-2"));
        simple_36_into(&mut buf, &null_variant(), &rg, "SAMPLE", "");
        assert_eq!(buf, simple_36(&null_variant(), &rg, "SAMPLE", ""));
        assert_eq!(buf.split('\t').count(), 36);

        simple_38_into(&mut buf, &long, &rg, "SAMPLE", "1-2");
        assert_eq!(buf, simple_38(&long, &rg, "SAMPLE", "1-2"));
        simple_38_into(&mut buf, &null_variant(), &rg, "SAMPLE", "");
        assert_eq!(buf, simple_38(&null_variant(), &rg, "SAMPLE", ""));
        assert_eq!(buf.split('\t').count(), 38);

        let mk_av = |v: Variant| AmpliconVariant {
            variant: v,
            gdamp: 2,
            tlamp: 3,
            ncamp: 1,
            ampflag: true,
            region: "chr1:933867-934466".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };
        amplicon_38_into(&mut buf, &mk_av(long.clone()), &rg, "SAMPLE", false);
        assert_eq!(buf, amplicon_38(&mk_av(long.clone()), &rg, "SAMPLE", false));
        amplicon_40_into(&mut buf, &mk_av(long.clone()), &rg, "SAMPLE", false);
        assert_eq!(buf, amplicon_40(&mk_av(long.clone()), &rg, "SAMPLE", false));
        assert_eq!(buf.split('\t').count(), 40);

        somatic_55_into(
            &mut buf,
            &long,
            Some(&long),
            &rg,
            "SAMPLE",
            "StrongSomatic",
            "12",
            "",
            true,
        );
        assert_eq!(
            buf,
            somatic_55(
                &long,
                Some(&long),
                &rg,
                "SAMPLE",
                "StrongSomatic",
                "12",
                "",
                true
            )
        );
        assert_eq!(buf.split('\t').count(), 55);
        somatic_61_into(
            &mut buf,
            &long,
            Some(&long),
            &rg,
            "SAMPLE",
            "StrongSomatic",
            "",
            "34",
            true,
        );
        assert_eq!(
            buf,
            somatic_61(
                &long,
                Some(&long),
                &rg,
                "SAMPLE",
                "StrongSomatic",
                "",
                "34",
                true
            )
        );
        assert_eq!(buf.split('\t').count(), 61);
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
        // Keep-mode pins (were fmt_3dp_msi; same values through push_keep).
        let mut msi = String::new();
        push_keep(&mut msi, 0.0, 3);
        assert_eq!(msi, "0");
        msi.clear();
        push_keep(&mut msi, 1.125, 3);
        assert_eq!(msi, "1.125");
        msi.clear();
        push_keep(&mut msi, 1.0005, 3);
        assert_eq!(msi, "1.000");

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
