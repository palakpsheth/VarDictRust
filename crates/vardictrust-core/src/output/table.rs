//! Table output serialization formats — intermediate .var output formatters.
//! Ports SimpleOutputVariant.java, AmpliconOutputVariant.java, and SomaticOutputVariant.java.

use crate::regions::Region;
use crate::variations::Variant;
use std::fmt::Write as _;

#[inline]
fn push_tsv_clean(buf: &mut String, s: &str) {
    if !s
        .as_bytes()
        .iter()
        .any(|&b| matches!(b, b'\t' | b'\r' | b'\n' | b'\0'))
    {
        buf.push_str(s);
        return;
    }
    for ch in s.chars() {
        if matches!(ch, '\t' | '\r' | '\n' | '\0') {
            buf.push('_');
        } else {
            buf.push(ch);
        }
    }
}

/// Print the `-h/--header` column line for the active mode to `w`.
/// Ports CmdParser.java:67 (`printHeader`) and the per-mode headers:
/// SimpleMode.java:114-127, SomaticMode.java:134-150, AmpliconMode.java:135-146.
/// Java prints this unconditionally of output format, before any variant rows.
pub fn print_header<W: std::io::Write>(
    w: &mut W,
    config: &crate::config::Config,
) -> std::io::Result<()> {
    const PREFIX_7: &[&str] = &["Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt"];
    const SINGLE_MIDDLE_27: &[&str] = &[
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
    ];
    const SOMATIC_SAMPLE_18: &[&str] = &[
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
    ];
    const SOMATIC_TAIL_12: &[&str] = &[
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

    if config.splice {
        writeln!(w, "Sample\tChr\tIntron\tIntron count")?;
        return Ok(());
    }
    let mut names: Vec<&str> = Vec::with_capacity(55);
    names.extend_from_slice(PREFIX_7);
    if config.region.is_none() && config.amplicon.is_some() {
        names.extend_from_slice(SINGLE_MIDDLE_27);
        names.extend_from_slice(&["GoodVarCount", "TotalVarCount", "Nocov", "Ampflag"]);
    } else if config.bam.contains('|') {
        names.extend_from_slice(SOMATIC_SAMPLE_18);
        names.extend_from_slice(SOMATIC_SAMPLE_18);
        names.extend_from_slice(SOMATIC_TAIL_12);
    } else {
        names.extend_from_slice(SINGLE_MIDDLE_27);
        names.extend_from_slice(&["Duprate", "SV_info"]);
        if matches!(config.crispr_cutting_site, Some(n) if n != 0) {
            // SimpleMode.java:124-126: CRISPR mode appends a column.
            names.push("CRISPR");
        }
    }
    writeln!(w, "{}", names.join("\t"))?;
    Ok(())
}

/// Explicit `.var` serialization dialect across VarDictJava v1.8.3's 6 output schemas (`AR-09`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TableModeStyle {
    Simple36,
    Simple38,
    Amplicon38,
    Amplicon40,
    Somatic55,
    Somatic61,
}

impl TableModeStyle {
    #[inline]
    const fn is_fisher(self) -> bool {
        matches!(self, Self::Simple38 | Self::Amplicon40 | Self::Somatic61)
    }
}

/// `.var` row formatting configuration pairing the column schema (`TableModeStyle`)
/// with the rounding tie mode (`exact_binary_half_even`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct TableFmt {
    style: TableModeStyle,
    exact_binary_half_even: bool,
}

impl TableFmt {
    #[inline]
    const fn new(style: TableModeStyle, exact_binary_half_even: bool) -> Self {
        Self {
            style,
            exact_binary_half_even,
        }
    }

    #[inline]
    const fn is_fisher(self) -> bool {
        self.style.is_fisher()
    }

    #[inline]
    fn push_float(self, buf: &mut String, v: f64, digits: i32) {
        if self.is_fisher() {
            push_java_rounded_compat(buf, v, digits, self.exact_binary_half_even);
        } else {
            push_keep_compat(buf, v, digits, self.exact_binary_half_even);
        }
    }

    /// `AR-09` Quirk 1: `Simple38` & `Amplicon40` keep trailing zeros on `hi_freq`,
    /// whereas `Somatic61` strips trailing zeros (`push_java_rounded_compat`).
    #[inline]
    fn push_hi_freq(self, buf: &mut String, hi_freq: f64) {
        match self.style {
            TableModeStyle::Somatic61 => {
                push_java_rounded_compat(buf, hi_freq, 4, self.exact_binary_half_even);
            }
            TableModeStyle::Simple36
            | TableModeStyle::Simple38
            | TableModeStyle::Amplicon38
            | TableModeStyle::Amplicon40
            | TableModeStyle::Somatic55 => {
                push_keep_compat(buf, hi_freq, 4, self.exact_binary_half_even);
            }
        }
    }

    /// `AR-09` Quirk 2: `Simple38` & `Amplicon40` format `nm` via `push_nm_compat` (keeping `.0`),
    /// whereas `Somatic61` formats positive `nm` via `push_java_rounded_compat(..., 1, ...)`.
    #[inline]
    fn push_nm_cell(self, buf: &mut String, nm: f64) {
        match self.style {
            TableModeStyle::Somatic61 => {
                if nm > 0.0 {
                    push_java_rounded_compat(buf, nm, 1, self.exact_binary_half_even);
                } else {
                    buf.push('0');
                }
            }
            TableModeStyle::Simple36
            | TableModeStyle::Simple38
            | TableModeStyle::Amplicon38
            | TableModeStyle::Amplicon40
            | TableModeStyle::Somatic55 => {
                push_nm_compat(buf, nm, self.exact_binary_half_even);
            }
        }
    }

    /// `AR-09` Quirk 3: `Somatic61` keeps trailing zeros on `msi` (`push_keep_compat(..., 3, ...)`),
    /// whereas `Simple38` & `Amplicon40` strip trailing zeros (`push_java_rounded_compat(..., 3, ...)`).
    #[inline]
    fn push_msi_cell(self, buf: &mut String, msi: f64) {
        match self.style {
            TableModeStyle::Simple38 | TableModeStyle::Amplicon40 => {
                push_java_rounded_compat(buf, msi, 3, self.exact_binary_half_even);
            }
            TableModeStyle::Simple36
            | TableModeStyle::Amplicon38
            | TableModeStyle::Somatic55
            | TableModeStyle::Somatic61 => {
                push_keep_compat(buf, msi, 3, self.exact_binary_half_even);
            }
        }
    }

    /// `AR-09` Quirk 4: `Simple38` rounds `duprate` to 2 decimals, whereas `Simple36`,
    /// `Somatic55`, and `Somatic61` format `duprate` to 1 decimal.
    #[inline]
    fn push_duprate_cell(self, buf: &mut String, duprate: f64) {
        match self.style {
            TableModeStyle::Simple38 => {
                push_java_rounded_compat(buf, duprate, 2, self.exact_binary_half_even);
            }
            TableModeStyle::Somatic61 => {
                push_java_rounded_compat(buf, duprate, 1, self.exact_binary_half_even);
            }
            TableModeStyle::Simple36
            | TableModeStyle::Somatic55
            | TableModeStyle::Amplicon38
            | TableModeStyle::Amplicon40 => {
                push_keep_compat(buf, duprate, 1, self.exact_binary_half_even);
            }
        }
    }
}

#[inline]
fn push_row_prefix_7(buf: &mut String, sample_name: &str, region: &Region, v: &Variant) {
    push_tsv_clean(buf, sample_name);
    buf.push('\t');
    push_tsv_clean(buf, &region.gene);
    buf.push('\t');
    push_tsv_clean(buf, &region.chr);
    buf.push('\t');
    crate::util::push_i32(buf, v.start_pos);
    buf.push('\t');
    crate::util::push_i32(buf, v.end_pos);
    buf.push('\t');
    buf.push_str(&v.ref_allele);
    buf.push('\t');
    buf.push_str(&v.var_allele);
    buf.push('\t');
}

#[inline]
fn push_allele_counts_6(buf: &mut String, v: &Variant) {
    crate::util::push_i32(buf, v.total_pos_cov);
    buf.push('\t');
    crate::util::push_i32(buf, v.pos_coverage);
    buf.push('\t');
    crate::util::push_i32(buf, v.ref_fwd_cov);
    buf.push('\t');
    crate::util::push_i32(buf, v.ref_rev_cov);
    buf.push('\t');
    crate::util::push_i32(buf, v.fwd);
    buf.push('\t');
    crate::util::push_i32(buf, v.rev);
    buf.push('\t');
}

#[inline]
fn push_region_span(buf: &mut String, region: &Region) {
    push_tsv_clean(buf, &region.chr);
    buf.push(':');
    crate::util::push_u32(buf, region.start.0);
    buf.push('-');
    crate::util::push_u32(buf, region.end.0);
}

#[inline]
fn push_flank_pair(buf: &mut String, v: &Variant, honor_null_record: bool) {
    let empty_cell = if honor_null_record && v.null_record {
        "\t"
    } else {
        "0\t"
    };
    if v.left_seq.is_empty() {
        buf.push_str(empty_cell);
    } else {
        buf.push_str(&v.left_seq);
        buf.push('\t');
    }
    if v.right_seq.is_empty() {
        buf.push_str(empty_cell);
    } else {
        buf.push_str(&v.right_seq);
        buf.push('\t');
    }
}

#[inline]
fn push_nonempty_or_zero(buf: &mut String, s: &str) {
    if s.is_empty() {
        buf.push('0');
    } else {
        buf.push_str(s);
    }
}

/// Shared columns 1..=34 (no-fisher) / 1..=36 (`--fisher`) for `Simple` and `Amplicon` modes.
fn push_single_core(
    buf: &mut String,
    variant: &Variant,
    sample_name: &str,
    region: &Region,
    amplicon_seg: Option<&str>,
    fmt: TableFmt,
) {
    push_row_prefix_7(buf, sample_name, region, variant);
    push_allele_counts_6(buf, variant);
    // Null records render printer field defaults (`""`), not `"0"` (`SimpleOutputVariant.java:50-52`).
    if !variant.null_record {
        buf.push_str(variant.genotype.as_deref().unwrap_or("0"));
    }
    buf.push('\t');
    fmt.push_float(buf, variant.frequency, 4);
    buf.push('\t');
    buf.push_str(&variant.bias);
    buf.push('\t');
    fmt.push_float(buf, variant.mean_position, 1);
    buf.push('\t');
    buf.push(if variant.pstd { '1' } else { '0' });
    buf.push('\t');
    fmt.push_float(buf, variant.mean_quality, 1);
    buf.push('\t');
    buf.push(if variant.qstd { '1' } else { '0' });
    buf.push('\t');
    if fmt.is_fisher() {
        let fisher = crate::stats::FisherExact::new_java(
            variant.ref_fwd_cov,
            variant.ref_rev_cov,
            variant.fwd,
            variant.rev,
        );
        push_java_rounded_compat(buf, fisher.get_p_value(), 5, fmt.exact_binary_half_even);
        buf.push('\t');
        fisher.write_odds_ratio_into(buf);
        buf.push('\t');
    }
    fmt.push_float(buf, variant.mean_mapping_quality, 1);
    buf.push('\t');
    fmt.push_float(buf, variant.hq_lq_ratio, 3);
    buf.push('\t');
    fmt.push_hi_freq(buf, variant.hi_freq);
    buf.push('\t');
    fmt.push_float(buf, variant.extra_freq, 4);
    buf.push('\t');
    crate::util::push_i32(buf, variant.shift3);
    buf.push('\t');
    fmt.push_msi_cell(buf, variant.msi);
    buf.push('\t');
    crate::util::push_i32(buf, variant.msint);
    buf.push('\t');
    fmt.push_nm_cell(buf, variant.nm);
    buf.push('\t');
    crate::util::push_i32(buf, variant.hicnt);
    buf.push('\t');
    crate::util::push_i32(buf, variant.hicov);
    buf.push('\t');
    push_flank_pair(buf, variant, true);
    match amplicon_seg {
        Some(seg) => push_tsv_clean(buf, seg),
        None => push_region_span(buf, region),
    }
    buf.push('\t');
    buf.push_str(&variant.var_type);
}

fn simple_append_styled(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    sample_name: &str,
    sv: &str,
    fmt: TableFmt,
) {
    buf.reserve(if fmt.is_fisher() { 320 } else { 256 });
    push_single_core(buf, variant, sample_name, region, None, fmt);
    buf.push('\t');
    fmt.push_duprate_cell(buf, variant.duprate);
    buf.push('\t');
    push_nonempty_or_zero(buf, sv);
}

fn amplicon_append_styled(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
    fmt: TableFmt,
) {
    buf.reserve(if fmt.is_fisher() { 320 } else { 256 });
    let variant = &av.variant;
    push_single_core(buf, variant, sample_name, region, Some(&av.region), fmt);
    buf.push('\t');
    crate::util::push_i32(buf, av.gdamp);
    buf.push('\t');
    crate::util::push_i32(buf, av.tlamp);
    buf.push('\t');
    crate::util::push_i32(buf, av.ncamp);
    buf.push('\t');
    buf.push(if av.ampflag { '1' } else { '0' });

    if !debug {
        return;
    }
    // ports AmpliconOutputVariant.java:113-125
    buf.push('\t');
    buf.push_str(&variant.debug);
    for (i, (v, r)) in av.good_variants.iter().enumerate() {
        buf.push_str("\tGood");
        crate::util::push_usize(buf, i);
        buf.push(' ');
        buf.push_str(&debug_amp_variant_compat(
            v.as_ref(),
            fmt.exact_binary_half_even,
        ));
        buf.push(' ');
        buf.push_str(r);
    }
    for (i, (v, r)) in av.bad_variants.iter().enumerate() {
        buf.push_str("\tBad");
        crate::util::push_usize(buf, i);
        buf.push(' ');
        buf.push_str(&debug_amp_variant_compat(
            v.as_ref(),
            fmt.exact_binary_half_even,
        ));
        buf.push(' ');
        buf.push_str(r);
    }
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
    simple_36_append(buf, variant, region, sample_name, sv);
}

/// Append the simple 36-column row to `buf` without clearing it first.
pub(crate) fn simple_36_append(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    sample_name: &str,
    sv: &str,
) {
    simple_append_styled(
        buf,
        variant,
        region,
        sample_name,
        sv,
        TableFmt::new(TableModeStyle::Simple36, false),
    );
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
    simple_38_append(buf, variant, region, sample_name, sv);
}

/// Append the simple 38-column row to `buf` without clearing it first.
pub(crate) fn simple_38_append(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    sample_name: &str,
    sv: &str,
) {
    simple_append_styled(
        buf,
        variant,
        region,
        sample_name,
        sv,
        TableFmt::new(TableModeStyle::Simple38, false),
    );
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
pub fn amplicon_38_into(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) {
    buf.clear();
    amplicon_38_append(buf, av, region, sample_name, debug);
}

/// Append the amplicon 38-column row to `buf` without clearing it first.
pub(crate) fn amplicon_38_append(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) {
    amplicon_append_styled(
        buf,
        av,
        region,
        sample_name,
        debug,
        TableFmt::new(TableModeStyle::Amplicon38, false),
    );
}

/// Write amplicon 38-column representation (no fisher).
/// Ports AmpliconOutputVariant.java's create_amplicon_variant_38columns.
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
pub fn amplicon_40_into(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) {
    buf.clear();
    amplicon_40_append(buf, av, region, sample_name, debug);
}

/// Append the amplicon 40-column row to `buf` without clearing it first.
pub(crate) fn amplicon_40_append(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
) {
    amplicon_append_styled(
        buf,
        av,
        region,
        sample_name,
        debug,
        TableFmt::new(TableModeStyle::Amplicon40, false),
    );
}

/// Write amplicon 40-column representation with Fisher exact test (when config.fisher is set).
/// Ports AmpliconOutputVariant.java's create_amplicon_variant_40columns (lines 145-201).
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
fn fmt_4dp_compat(v: f64, exact_binary_half_even: bool) -> String {
    if v == 0.0 {
        "0".to_string()
    } else {
        format!(
            "{:.4}",
            crate::types::round_half_even_compat(v, 4, exact_binary_half_even)
        )
    }
}

#[cfg(test)]
fn debug_amp_variant(v: Option<&Variant>) -> String {
    debug_amp_variant_compat(v, false)
}

fn debug_amp_variant_compat(v: Option<&Variant>, exact_binary_half_even: bool) -> String {
    let default = Variant::new();
    let v = v.unwrap_or(&default);
    let genotype = v.genotype.as_deref().unwrap_or("0");
    let pstd = if v.pstd { 1 } else { 0 };
    let qstd = if v.qstd { 1 } else { 0 };
    let mean_pos = format!(
        "{:.1}",
        crate::types::round_half_even_compat(v.mean_position, 1, exact_binary_half_even)
    );
    let mean_qual = format!(
        "{:.1}",
        crate::types::round_half_even_compat(v.mean_quality, 1, exact_binary_half_even)
    );
    let mapq = format!(
        "{:.1}",
        crate::types::round_half_even_compat(v.mean_mapping_quality, 1, exact_binary_half_even)
    );
    let qratio = format!(
        "{:.3}",
        crate::types::round_half_even_compat(v.hq_lq_ratio, 3, exact_binary_half_even)
    );
    format!(
        "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        v.total_pos_cov,
        v.pos_coverage,
        v.ref_fwd_cov,
        v.ref_rev_cov,
        v.fwd,
        v.rev,
        genotype,
        fmt_4dp_compat(v.frequency, exact_binary_half_even),
        v.bias,
        mean_pos,
        pstd,
        mean_qual,
        qstd,
        mapq,
        qratio,
        fmt_4dp_compat(v.hi_freq, exact_binary_half_even),
        fmt_4dp_compat(v.extra_freq, exact_binary_half_even),
    )
}

/// Java-compatible formatter that rounds a float using banker's rounding
/// to `digits` decimal places and strips trailing `'0'` characters.
/// Ports `Utils.getRoundedValueToPrint` (`Utils.java:105-109`).
#[cfg(test)]
pub(crate) fn get_rounded_value_to_print(value: f64, digits: i32) -> String {
    let mut buf = String::new();
    push_java_rounded(&mut buf, value, digits);
    buf
}

/// Push `Utils.getRoundedValueToPrint` directly into `buf` without a temp `String`,
/// parameterized by `exact_binary_half_even`.
pub(crate) fn push_java_rounded_compat(
    buf: &mut String,
    value: f64,
    digits: i32,
    exact_binary_half_even: bool,
) {
    if value.is_nan() || value.is_infinite() {
        let _ = write!(buf, "{}", value);
        return;
    }
    if value == value.round() {
        if value == 0.0 {
            if value.is_sign_negative() {
                buf.push_str("-0");
            } else {
                buf.push('0');
            }
        } else {
            let _ = write!(buf, "{:.0}", value);
        }
        return;
    }
    let rounded = crate::types::round_half_even_compat(value, digits, exact_binary_half_even);
    let start = buf.len();
    let _ = write!(buf, "{:.*}", digits.max(0) as usize, rounded);
    while buf.len() > start && buf.as_bytes().get(buf.len() - 1) == Some(&b'0') {
        buf.pop();
    }
}

/// Push `Utils.getRoundedValueToPrint` directly into `buf` without a temp `String`.
/// Used by all `--fisher` `.var` writers (`SimpleOutputVariant` 38-col,
/// `AmpliconOutputVariant` 40-col, `SomaticOutputVariant` 61-col).
#[inline]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "default compat=vardictjava wrapper retained alongside push_java_rounded_compat"
    )
)]
pub(crate) fn push_java_rounded(buf: &mut String, value: f64, digits: i32) {
    push_java_rounded_compat(buf, value, digits, false);
}

/// Keep-mode float cell parameterized by `exact_binary_half_even`: HALF_EVEN at `digits`
/// decimals with trailing zeros kept (`"0"` when zero).
pub(crate) fn push_keep_compat(
    buf: &mut String,
    v: f64,
    digits: i32,
    exact_binary_half_even: bool,
) {
    if v == 0.0 {
        buf.push('0');
    } else {
        let _ = write!(
            buf,
            "{:.*}",
            digits as usize,
            crate::types::round_half_even_compat(v, digits, exact_binary_half_even)
        );
    }
}

/// Keep-mode float cell: HALF_EVEN at `digits` decimals with trailing zeros
/// kept (`"0"` when zero). Ports the `DecimalFormat("0.0000"/"0.0"/"0.000")`
/// cells of the no-`--fisher` writers.
#[inline]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "default compat=vardictjava wrapper retained alongside push_keep_compat"
    )
)]
pub(crate) fn push_keep(buf: &mut String, v: f64, digits: i32) {
    push_keep_compat(buf, v, digits, false);
}

/// Strip-mode float cell for R `write.table` / Perl `var2vcf_*.pl` VCF output,
/// parameterized by `exact_binary_half_even`: HALF_EVEN at `digits` decimals with
/// trailing zeros AND trailing decimal point stripped (`"0"` when rounded == 0,
/// integers without a dot).
pub(crate) fn push_stripped_compat(
    buf: &mut String,
    value: f64,
    digits: i32,
    exact_binary_half_even: bool,
) {
    if value.is_nan() || value.is_infinite() {
        let _ = write!(buf, "{}", value);
        return;
    }
    let rounded = crate::types::round_half_even_compat(value, digits, exact_binary_half_even);
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

/// Strip-mode float cell for R `write.table` / Perl `var2vcf_*.pl` VCF output:
/// HALF_EVEN at `digits` decimals with trailing zeros AND trailing decimal point
/// stripped (`"0"` when rounded == 0, integers without a dot).
#[inline]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "default compat=vardictjava wrapper retained alongside push_stripped_compat"
    )
)]
pub(crate) fn push_stripped(buf: &mut String, value: f64, digits: i32) {
    push_stripped_compat(buf, value, digits, false);
}

/// NM cell parameterized by `exact_binary_half_even`: keep-mode 1dp when positive, `"0"` otherwise.
pub(crate) fn push_nm_compat(buf: &mut String, v: f64, exact_binary_half_even: bool) {
    if v > 0.0 {
        push_keep_compat(buf, v, 1, exact_binary_half_even);
    } else {
        buf.push('0');
    }
}

/// NM cell: keep-mode 1dp when positive, `"0"` otherwise (matches `fmt_nm`).
#[inline]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "default compat=vardictjava wrapper retained alongside push_nm_compat"
    )
)]
pub(crate) fn push_nm(buf: &mut String, v: f64) {
    push_nm_compat(buf, v, false);
}

/// Push the 18-column (`Somatic55`) or 20-column (`Somatic61`) per-sample allele block.
fn push_somatic_sample_cols(buf: &mut String, v: &Variant, fmt: TableFmt) {
    push_allele_counts_6(buf, v);
    buf.push_str(v.genotype.as_deref().unwrap_or("0"));
    buf.push('\t');
    fmt.push_float(buf, v.frequency, 4);
    buf.push('\t');
    buf.push_str(&v.bias);
    buf.push('\t');
    fmt.push_float(buf, v.mean_position, 1);
    buf.push('\t');
    buf.push(if v.pstd { '1' } else { '0' });
    buf.push('\t');
    fmt.push_float(buf, v.mean_quality, 1);
    buf.push('\t');
    buf.push(if v.qstd { '1' } else { '0' });
    buf.push('\t');
    fmt.push_float(buf, v.mean_mapping_quality, 1);
    buf.push('\t');
    fmt.push_float(buf, v.hq_lq_ratio, 3);
    buf.push('\t');
    fmt.push_hi_freq(buf, v.hi_freq);
    buf.push('\t');
    fmt.push_float(buf, v.extra_freq, 4);
    buf.push('\t');
    fmt.push_nm_cell(buf, v.nm);
    if fmt.is_fisher() {
        let fisher =
            crate::stats::FisherExact::new_java(v.ref_fwd_cov, v.ref_rev_cov, v.fwd, v.rev);
        buf.push('\t');
        push_java_rounded_compat(buf, fisher.get_p_value(), 5, fmt.exact_binary_half_even);
        buf.push('\t');
        fisher.write_odds_ratio_into(buf);
    }
}

const VARIANT_18COLS_EMPTY: &str = "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";

#[cfg(test)]
fn variant_20cols_fisher_empty(pvalue: f64, oddratio: &str) -> String {
    format!(
        "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t{}\t{}",
        get_rounded_value_to_print(pvalue, 5),
        oddratio,
    )
}

const VARIANT_20COLS_FISHER_EMPTY: &str =
    "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t0";

#[expect(
    clippy::too_many_arguments,
    reason = "paired somatic row formatter mirrors SomaticOutputVariant 55/61-column arguments"
)]
fn somatic_append_styled(
    buf: &mut String,
    tumor: &Variant,
    normal: Option<&Variant>,
    region: &Region,
    sample_name: &str,
    status: &str,
    sv1: &str,
    sv2: &str,
    end_is_variant2: bool,
    fmt: TableFmt,
) {
    buf.reserve(if fmt.is_fisher() { 448 } else { 384 });
    // endVariant supplies shift3/msi/msint/left/right seq (Java SomaticOutputVariant).
    let end_var = if end_is_variant2 {
        normal.unwrap_or(tumor)
    } else {
        tumor
    };

    push_row_prefix_7(buf, sample_name, region, tumor);
    push_somatic_sample_cols(buf, tumor, fmt);
    buf.push('\t');
    match normal {
        Some(n) => push_somatic_sample_cols(buf, n, fmt),
        None => buf.push_str(if fmt.is_fisher() {
            VARIANT_20COLS_FISHER_EMPTY
        } else {
            VARIANT_18COLS_EMPTY
        }),
    }
    buf.push('\t');
    crate::util::push_i32(buf, end_var.shift3);
    buf.push('\t');
    fmt.push_msi_cell(buf, end_var.msi);
    buf.push('\t');
    crate::util::push_i32(buf, end_var.msint);
    buf.push('\t');
    push_flank_pair(buf, end_var, false);
    push_region_span(buf, region);
    buf.push('\t');
    buf.push_str(status);
    buf.push('\t');
    buf.push_str(&tumor.var_type);
    buf.push('\t');
    fmt.push_duprate_cell(buf, tumor.duprate);
    buf.push('\t');
    push_nonempty_or_zero(buf, sv1);
    buf.push('\t');
    match normal {
        Some(n) => fmt.push_duprate_cell(buf, n.duprate),
        None => buf.push('0'),
    }
    buf.push('\t');
    push_nonempty_or_zero(buf, sv2);

    if fmt.is_fisher() {
        let (tvar_cov, ttotal_cov) = (tumor.pos_coverage, tumor.total_pos_cov);
        let (nvar_cov, ntotal_cov) = match normal {
            Some(n) => (n.pos_coverage, n.total_pos_cov),
            None => (0, 0),
        };
        let tref =
            (i64::from(ttotal_cov) - i64::from(tvar_cov)).clamp(0, i64::from(i32::MAX)) as i32;
        let rref =
            (i64::from(ntotal_cov) - i64::from(nvar_cov)).clamp(0, i64::from(i32::MAX)) as i32;
        let fisher_total = crate::stats::FisherExact::new_java(tvar_cov, tref, nvar_cov, rref);
        let pvalue_greater = fisher_total.get_p_value_greater();
        let pvalue_less = fisher_total.get_p_value_less();
        let total_pvalue = if pvalue_less < pvalue_greater {
            pvalue_less
        } else {
            pvalue_greater
        };
        buf.push('\t');
        push_java_rounded_compat(buf, total_pvalue, 5, fmt.exact_binary_half_even);
        buf.push('\t');
        fisher_total.write_odds_ratio_into(buf);
    }
}

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
    somatic_55_append(
        buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
    );
}

/// Append the somatic 55-column row to `buf` without clearing it first.
pub(crate) fn somatic_55_append(
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
    somatic_append_styled(
        buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
        TableFmt::new(TableModeStyle::Somatic55, false),
    );
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
    somatic_61_append(
        buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
    );
}

/// Append the somatic 61-column row to `buf` without clearing it first.
pub(crate) fn somatic_61_append(
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
    somatic_append_styled(
        buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
        TableFmt::new(TableModeStyle::Somatic61, false),
    );
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

pub(crate) fn simple_row_append_compat(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    sample_name: &str,
    sv: &str,
    fisher: bool,
    exact_binary_half_even: bool,
) {
    let style = if fisher {
        TableModeStyle::Simple38
    } else {
        TableModeStyle::Simple36
    };
    simple_append_styled(
        buf,
        variant,
        region,
        sample_name,
        sv,
        TableFmt::new(style, exact_binary_half_even),
    );
}

#[doc(hidden)]
pub fn simple_row_into_compat(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    sample_name: &str,
    sv: &str,
    fisher: bool,
    exact_binary_half_even: bool,
) {
    buf.clear();
    simple_row_append_compat(
        buf,
        variant,
        region,
        sample_name,
        sv,
        fisher,
        exact_binary_half_even,
    );
}

pub(crate) fn amplicon_row_append_compat(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
    fisher: bool,
    exact_binary_half_even: bool,
) {
    let style = if fisher {
        TableModeStyle::Amplicon40
    } else {
        TableModeStyle::Amplicon38
    };
    amplicon_append_styled(
        buf,
        av,
        region,
        sample_name,
        debug,
        TableFmt::new(style, exact_binary_half_even),
    );
}

#[doc(hidden)]
pub fn amplicon_row_into_compat(
    buf: &mut String,
    av: &crate::postprocess::amplicon::AmpliconVariant,
    region: &Region,
    sample_name: &str,
    debug: bool,
    fisher: bool,
    exact_binary_half_even: bool,
) {
    buf.clear();
    amplicon_row_append_compat(
        buf,
        av,
        region,
        sample_name,
        debug,
        fisher,
        exact_binary_half_even,
    );
}

#[expect(
    clippy::too_many_arguments,
    reason = "paired somatic row formatter mirrors SomaticOutputVariant 55/61-column arguments plus fisher/compat selectors"
)]
pub(crate) fn somatic_row_append_compat(
    buf: &mut String,
    tumor: &Variant,
    normal: Option<&Variant>,
    region: &Region,
    sample_name: &str,
    status: &str,
    sv1: &str,
    sv2: &str,
    end_is_variant2: bool,
    fisher: bool,
    exact_binary_half_even: bool,
) {
    let style = if fisher {
        TableModeStyle::Somatic61
    } else {
        TableModeStyle::Somatic55
    };
    somatic_append_styled(
        buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
        TableFmt::new(style, exact_binary_half_even),
    );
}

#[doc(hidden)]
#[expect(
    clippy::too_many_arguments,
    reason = "paired somatic row formatter mirrors SomaticOutputVariant 55/61-column arguments plus fisher/compat selectors"
)]
pub fn somatic_row_into_compat(
    buf: &mut String,
    tumor: &Variant,
    normal: Option<&Variant>,
    region: &Region,
    sample_name: &str,
    status: &str,
    sv1: &str,
    sv2: &str,
    end_is_variant2: bool,
    fisher: bool,
    exact_binary_half_even: bool,
) {
    buf.clear();
    somatic_row_append_compat(
        buf,
        tumor,
        normal,
        region,
        sample_name,
        status,
        sv1,
        sv2,
        end_is_variant2,
        fisher,
        exact_binary_half_even,
    );
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
        assert_eq!(get_rounded_value_to_print(-0.00001, 4), "-0.");
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
    }

    #[test]
    fn java_rounded_value_to_print_pins() {
        assert_eq!(get_rounded_value_to_print(0.023, 1), "0.");
        assert_eq!(get_rounded_value_to_print(-0.023, 1), "-0.");
        assert_eq!(get_rounded_value_to_print(0.96, 1), "1.");
        assert_eq!(get_rounded_value_to_print(0.003, 2), "0.");
        assert_eq!(get_rounded_value_to_print(0.996, 2), "1.");
        assert_eq!(get_rounded_value_to_print(35.96, 1), "36.");
        assert_eq!(get_rounded_value_to_print(1.0 / 30000.0, 4), "0.");
        assert_eq!(get_rounded_value_to_print(-0.00001, 4), "-0.");
        assert_eq!(get_rounded_value_to_print(2.00001, 4), "2.");
        assert_eq!(get_rounded_value_to_print(2001.0 / 2000.0, 3), "1.");
        assert_eq!(get_rounded_value_to_print(29.94, 1), "29.9");
        assert_eq!(get_rounded_value_to_print(30.0, 1), "30");
        assert_eq!(get_rounded_value_to_print(0.0, 4), "0");
        assert_eq!(get_rounded_value_to_print(-0.0, 4), "-0");

        // VCF R-echo stripped path keeps -0.00001 -> "0" (#60 item 18 guard).
        let mut vcf_buf = String::new();
        push_stripped(&mut vcf_buf, -0.00001, 4);
        assert_eq!(vcf_buf, "0");
    }

    #[test]
    fn format_oracle_table_parity() {
        const ORACLE: &str = include_str!("../../testdata/format/java_fmt_oracle.tsv");
        // The 4 JDK DecimalFormat sub-ulp underflow tie inputs (±5e-(d+1) for d in {3,4,5,8})
        let mut checked = 0usize;
        for (line_idx, line) in ORACLE.lines().enumerate() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 9, "malformed oracle row {}", line_idx + 1);
            let digits: i32 = cols[0].parse().expect("digits");
            let v_bits = u64::from_str_radix(cols[3], 16).expect("value_bits_hex");
            let value = f64::from_bits(v_bits);
            let expected_gvtp = cols[6];
            let expected_dts = cols[8];

            assert_eq!(
                crate::variations::java_double_to_string(value),
                expected_dts,
                "Double.toString mismatch on row {} (value={})",
                line_idx + 1,
                cols[2]
            );

            let actual_gvtp = get_rounded_value_to_print(value, digits);
            assert_eq!(
                actual_gvtp,
                expected_gvtp,
                "getRoundedValueToPrint mismatch on row {} (digits={}, value={})",
                line_idx + 1,
                digits,
                cols[2]
            );
            checked += 1;
        }
        assert_eq!(checked, 432, "expected all 432 oracle rows, got {checked}");
    }

    #[test]
    fn amplicon_null_record_and_header_precedence() {
        let rg = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut v = Variant::new();
        v.start_pos = 100;
        v.end_pos = 100;
        v.bias = "0;0".to_string();
        v.null_record = true;
        let av = AmpliconVariant {
            variant: v,
            gdamp: 0,
            tlamp: 0,
            ncamp: 0,
            ampflag: false,
            region: "chr1:100-100".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };
        let row38 = amplicon_38(&av, &rg, "S", false);
        let c38: Vec<&str> = row38.split('\t').collect();
        assert_eq!(c38.len(), 38);
        assert_eq!(c38[13], "", "null_record genotype must be empty");
        assert_eq!(c38[15], "0;0", "null_record bias must be '0;0'");
        assert_eq!(c38[30], "", "null_record 5pFlankSeq must be empty");
        assert_eq!(c38[31], "", "null_record 3pFlankSeq must be empty");

        let row40 = amplicon_40(&av, &rg, "S", false);
        let c40: Vec<&str> = row40.split('\t').collect();
        assert_eq!(c40.len(), 40);
        assert_eq!(c40[13], "", "null_record genotype must be empty");
        assert_eq!(c40[15], "0;0", "null_record bias must be '0;0'");
        assert_eq!(c40[32], "", "null_record 5pFlankSeq must be empty");
        assert_eq!(c40[33], "", "null_record 3pFlankSeq must be empty");

        // Header precedence: Splicing > (-R overrides Amplicon) > Amplicon > Somatic
        let mut cfg = crate::config::Config::default();
        cfg.splice = true;
        cfg.bam = "t.bam|n.bam".to_string();
        let mut h = Vec::new();
        print_header(&mut h, &cfg).unwrap();
        assert_eq!(
            String::from_utf8(h).unwrap(),
            "Sample\tChr\tIntron\tIntron count\n"
        );

        cfg.splice = false;
        cfg.amplicon = Some("10:0.95".to_string());
        let mut h_amp = Vec::new();
        print_header(&mut h_amp, &cfg).unwrap();
        let h_amp_str = String::from_utf8(h_amp).unwrap();
        assert!(
            h_amp_str.ends_with("GoodVarCount\tTotalVarCount\tNocov\tAmpflag\n"),
            "Amplicon mode takes precedence over paired BAM when -R is absent: {h_amp_str}"
        );
    }

    #[test]
    fn variant_20cols_fisher_empty_constant_matches_java_zero_table() {
        let f2 = crate::stats::FisherExact::new_java(0, 0, 0, 0);
        assert_eq!(
            VARIANT_20COLS_FISHER_EMPTY,
            variant_20cols_fisher_empty(f2.get_p_value(), &f2.get_odds_ratio())
        );
    }

    #[test]
    fn somatic_61_into_extreme_coverage_clamps_without_overflow() {
        let rg = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut tumor = Variant::new();
        tumor.start_pos = 100;
        tumor.end_pos = 100;
        tumor.ref_allele = "A".to_string();
        tumor.var_allele = "G".to_string();
        tumor.var_type = "SNV".to_string();
        tumor.total_pos_cov = i32::MAX;
        tumor.pos_coverage = i32::MIN;

        let mut normal = tumor.clone();
        normal.total_pos_cov = i32::MIN;
        normal.pos_coverage = i32::MAX;

        let line = somatic_61(
            &tumor,
            Some(&normal),
            &rg,
            "SAMPLE",
            "StrongSomatic",
            "",
            "34",
            true,
        );
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 61);
    }

    #[test]
    fn table_writers_sanitize_residual_tsv_delimiters_in_sample_chr_and_gene() {
        let dirty_region = Region {
            chr: "chr1\tINJ\nLINE".to_string(),
            start: Pos1(10),
            end: Pos1(20),
            gene: "TP53\tCOL\r\n\0".to_string(),
            amplicon: None,
        };
        let dirty_sample = "SAMP\t1\n2\r3\0";
        let mut v = Variant::new();
        v.start_pos = 100;
        v.end_pos = 100;
        v.ref_allele = "A".to_string();
        v.var_allele = "G".to_string();
        v.var_type = "SNV".to_string();
        let av = crate::postprocess::amplicon::AmpliconVariant {
            variant: v.clone(),
            gdamp: 1,
            tlamp: 1,
            ncamp: 0,
            ampflag: false,
            region: "chr1\t10-20\n".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };

        for (label, expected_cols, row) in [
            (
                "simple_36",
                36usize,
                simple_36(&v, &dirty_region, dirty_sample, ""),
            ),
            (
                "simple_38",
                38usize,
                simple_38(&v, &dirty_region, dirty_sample, ""),
            ),
            (
                "amplicon_38",
                38usize,
                amplicon_38(&av, &dirty_region, dirty_sample, false),
            ),
            (
                "amplicon_40",
                40usize,
                amplicon_40(&av, &dirty_region, dirty_sample, false),
            ),
            (
                "somatic_55",
                55usize,
                somatic_55(
                    &v,
                    Some(&v),
                    &dirty_region,
                    dirty_sample,
                    "StrongSomatic",
                    "",
                    "",
                    false,
                ),
            ),
            (
                "somatic_61",
                61usize,
                somatic_61(
                    &v,
                    Some(&v),
                    &dirty_region,
                    dirty_sample,
                    "StrongSomatic",
                    "",
                    "",
                    false,
                ),
            ),
        ] {
            assert!(
                !row.contains(['\n', '\r', '\0']),
                "{label} leaked control delimiter: {row:?}"
            );
            assert_eq!(
                row.split('\t').count(),
                expected_cols,
                "{label} column count drifted under dirty inputs: {row:?}"
            );
        }
    }

    #[test]
    fn ar09_table_mode_style_preserves_all_four_vardictjava_quirks() {
        let reg = Region {
            chr: "chr1".to_string(),
            start: Pos1(933866),
            end: Pos1(934466),
            gene: "testbed".to_string(),
            amplicon: None,
        };
        let mut v = Variant::new();
        v.start_pos = 100;
        v.end_pos = 100;
        v.ref_allele = "G".to_string();
        v.var_allele = "C".to_string();
        v.var_type = "SNV".to_string();
        v.total_pos_cov = 40;
        v.pos_coverage = 20;
        v.ref_fwd_cov = 10;
        v.ref_rev_cov = 10;
        v.fwd = 10;
        v.rev = 10;
        v.frequency = 0.5;
        v.hi_freq = 0.5;
        v.msi = 3.5;
        v.msint = 2;
        v.nm = 2.0;
        v.duprate = 0.125;

        let s36_line = simple_36(&v, &reg, "SAMPLE", "");
        let s36: Vec<&str> = s36_line.split('\t').collect();
        let s38_line = simple_38(&v, &reg, "SAMPLE", "");
        let s38: Vec<&str> = s38_line.split('\t').collect();
        let av = AmpliconVariant {
            variant: v.clone(),
            gdamp: 1,
            tlamp: 1,
            ncamp: 0,
            ampflag: false,
            region: "chr1:933867-934466".to_string(),
            good_variants: vec![],
            bad_variants: vec![],
        };
        let a38_line = amplicon_38(&av, &reg, "SAMPLE", false);
        let a38: Vec<&str> = a38_line.split('\t').collect();
        let a40_line = amplicon_40(&av, &reg, "SAMPLE", false);
        let a40: Vec<&str> = a40_line.split('\t').collect();
        let som55_line = somatic_55(&v, Some(&v), &reg, "SAMPLE", "StrongSomatic", "", "", false);
        let som55: Vec<&str> = som55_line.split('\t').collect();
        let som61_line = somatic_61(&v, Some(&v), &reg, "SAMPLE", "StrongSomatic", "", "", false);
        let som61: Vec<&str> = som61_line.split('\t').collect();

        // 1. hi_freq: Simple38 & Amplicon40 keep trailing zeros ("0.5000"); Somatic61 strips ("0.5").
        assert_eq!(s36[22], "0.5000");
        assert_eq!(s38[24], "0.5000");
        assert_eq!(a38[22], "0.5000");
        assert_eq!(a40[24], "0.5000");
        assert_eq!(som55[22], "0.5000");
        assert_eq!(som61[22], "0.5");
        assert_eq!(som61[42], "0.5");

        // 2. nm: Simple38 & Amplicon40 keep ".0" ("2.0"); Somatic61 strips ("2").
        assert_eq!(s36[27], "2.0");
        assert_eq!(s38[29], "2.0");
        assert_eq!(a38[27], "2.0");
        assert_eq!(a40[29], "2.0");
        assert_eq!(som55[24], "2.0");
        assert_eq!(som61[24], "2");
        assert_eq!(som61[44], "2");

        // 3. msi: Somatic61 keeps trailing zeros ("3.500"); Simple38 & Amplicon40 strip ("3.5").
        assert_eq!(s36[25], "3.500");
        assert_eq!(s38[27], "3.5");
        assert_eq!(a38[25], "3.500");
        assert_eq!(a40[27], "3.5");
        assert_eq!(som55[44], "3.500");
        assert_eq!(som61[48], "3.500");

        // 4. duprate: Simple38 rounds to 2dp ("0.12"); Simple36 ("0.1"), Somatic55 ("0.1"), Somatic61 ("0.1") use 1dp.
        assert_eq!(s36[34], "0.1");
        assert_eq!(s38[36], "0.12");
        assert_eq!(som55[51], "0.1");
        assert_eq!(som55[53], "0.1");
        assert_eq!(som61[55], "0.1");
        assert_eq!(som61[57], "0.1");
    }

    #[test]
    fn table_fmt_exact_binary_half_even_ties_and_debug_amp_variant() {
        let tie_3dp = f64::from_bits(0x3f40_624d_d2f1_a9fc); // 5e-4
        let tie_4dp = f64::from_bits(0x3f0a_36e2_eb1c_432d); // 5e-5
        let tie_5dp = f64::from_bits(0x3ed4_f8b5_88e3_68f1); // 5e-6

        let mut b = String::new();
        // Default wrappers (exact_binary_half_even = false)
        push_java_rounded(&mut b, tie_4dp, 4);
        assert_eq!(b, "0.");
        b.clear();
        push_keep(&mut b, tie_4dp, 4);
        assert_eq!(b, "0.0000");
        b.clear();
        push_stripped(&mut b, tie_4dp, 4);
        assert_eq!(b, "0");
        b.clear();
        push_nm(&mut b, 1.25);
        assert_eq!(b, "1.2");
        b.clear();
        push_nm_compat(&mut b, 0.0, true);
        assert_eq!(b, "0");

        for (val, digits, vj_java, vj_keep, vj_strip, off_str) in [
            (tie_3dp, 3, "0.", "0.000", "0", "0.001"),
            (tie_4dp, 4, "0.", "0.0000", "0", "0.0001"),
            (tie_5dp, 5, "0.", "0.00000", "0", "0.00001"),
        ] {
            b.clear();
            push_java_rounded_compat(&mut b, val, digits, false);
            assert_eq!(b, vj_java);
            b.clear();
            push_java_rounded_compat(&mut b, val, digits, true);
            assert_eq!(b, off_str);

            b.clear();
            push_keep_compat(&mut b, val, digits, false);
            assert_eq!(b, vj_keep);
            b.clear();
            push_keep_compat(&mut b, val, digits, true);
            assert_eq!(b, off_str);

            b.clear();
            push_stripped_compat(&mut b, val, digits, false);
            assert_eq!(b, vj_strip);
            b.clear();
            push_stripped_compat(&mut b, val, digits, true);
            assert_eq!(b, off_str);
        }

        let reg = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut v = Variant::new();
        v.start_pos = 100;
        v.end_pos = 100;
        v.ref_allele = "A".to_string();
        v.var_allele = "G".to_string();
        v.var_type = "SNV".to_string();
        v.frequency = tie_4dp;
        v.hi_freq = tie_4dp;
        v.extra_freq = tie_4dp;
        v.hq_lq_ratio = tie_3dp;
        v.msi = tie_3dp;
        v.debug = "DBG".to_string();

        // Simple 36 / 38 default vs COMPAT_OFF
        let mut row = String::new();
        simple_row_into_compat(&mut row, &v, &reg, "S", "", false, false);
        assert_eq!(row, simple_36(&v, &reg, "S", ""));
        let c36_vj: Vec<&str> = row.split('\t').collect();
        assert_eq!(c36_vj[14], "0.0000");
        assert_eq!(c36_vj[21], "0.000");
        assert_eq!(c36_vj[22], "0.0000");
        assert_eq!(c36_vj[23], "0.0000");
        assert_eq!(c36_vj[25], "0.000");

        simple_row_into_compat(&mut row, &v, &reg, "S", "", false, true);
        let c36_off: Vec<&str> = row.split('\t').collect();
        assert_eq!(c36_off[14], "0.0001");
        assert_eq!(c36_off[21], "0.001");
        assert_eq!(c36_off[22], "0.0001");
        assert_eq!(c36_off[23], "0.0001");
        assert_eq!(c36_off[25], "0.001");

        simple_row_into_compat(&mut row, &v, &reg, "S", "", true, false);
        assert_eq!(row, simple_38(&v, &reg, "S", ""));
        let c38_vj: Vec<&str> = row.split('\t').collect();
        assert_eq!(c38_vj[14], "0.");
        assert_eq!(c38_vj[23], "0.");
        assert_eq!(c38_vj[24], "0.0000");
        assert_eq!(c38_vj[25], "0.");
        assert_eq!(c38_vj[27], "0.");

        simple_row_into_compat(&mut row, &v, &reg, "S", "", true, true);
        let c38_off: Vec<&str> = row.split('\t').collect();
        assert_eq!(c38_off[14], "0.0001");
        assert_eq!(c38_off[23], "0.001");
        assert_eq!(c38_off[24], "0.0001");
        assert_eq!(c38_off[25], "0.0001");
        assert_eq!(c38_off[27], "0.001");

        // Amplicon 38 / 40 + -D debug_amp_variant (qratio = 5e-4, freqs = 5e-5)
        let av = AmpliconVariant {
            variant: v.clone(),
            gdamp: 1,
            tlamp: 1,
            ncamp: 0,
            ampflag: false,
            region: "chr1:1-200".to_string(),
            good_variants: vec![(Some(v.clone()), "chr1:1-200".to_string())],
            bad_variants: vec![(Some(v.clone()), "chr1:50-150".to_string())],
        };
        amplicon_row_into_compat(&mut row, &av, &reg, "S", true, false, false);
        assert_eq!(row, amplicon_38(&av, &reg, "S", true));
        assert!(
            row.contains(
                "Good0 0 0 0 0 0 0 0 0.0000 0 0.0 0 0.0 0 0.0 0.000 0.0000 0.0000 chr1:1-200"
            ),
            "unexpected default debug_amp_variant: {row}"
        );

        amplicon_row_into_compat(&mut row, &av, &reg, "S", true, false, true);
        assert!(
            row.contains(
                "Good0 0 0 0 0 0 0 0 0.0001 0 0.0 0 0.0 0 0.0 0.001 0.0001 0.0001 chr1:1-200"
            ),
            "unexpected COMPAT_OFF debug_amp_variant: {row}"
        );
        assert!(
            row.contains(
                "Bad0 0 0 0 0 0 0 0 0.0001 0 0.0 0 0.0 0 0.0 0.001 0.0001 0.0001 chr1:50-150"
            ),
            "unexpected COMPAT_OFF bad debug_amp_variant: {row}"
        );

        // Somatic 55 / 61 default vs COMPAT_OFF
        somatic_row_into_compat(
            &mut row,
            &v,
            Some(&v),
            &reg,
            "S",
            "StrongSomatic",
            "",
            "",
            false,
            false,
            false,
        );
        assert_eq!(
            row,
            somatic_55(&v, Some(&v), &reg, "S", "StrongSomatic", "", "", false)
        );
        somatic_row_into_compat(
            &mut row,
            &v,
            Some(&v),
            &reg,
            "S",
            "StrongSomatic",
            "",
            "",
            false,
            true,
            false,
        );
        assert_eq!(
            row,
            somatic_61(&v, Some(&v), &reg, "S", "StrongSomatic", "", "", false)
        );
        somatic_row_into_compat(
            &mut row,
            &v,
            Some(&v),
            &reg,
            "S",
            "StrongSomatic",
            "",
            "",
            false,
            true,
            true,
        );
        let cs61_off: Vec<&str> = row.split('\t').collect();
        assert_eq!(cs61_off[14], "0.0001");
        assert_eq!(cs61_off[21], "0.001");
        assert_eq!(cs61_off[22], "0.0001");
        assert_eq!(cs61_off[23], "0.0001");
        assert_eq!(cs61_off[48], "0.001");
    }
}
