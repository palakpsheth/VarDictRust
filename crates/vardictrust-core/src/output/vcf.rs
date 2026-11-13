//! VCF format output serialization and var2vcf_valid filtering.
//! Ports var2vcf_valid.pl.

use crate::config::Config;
use crate::output::table::{push_java_rounded, push_keep, push_stripped};
use crate::regions::Region;
use crate::variations::Variant;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::{self, Write};

/// Push a rounded float into `buf` matching the active pipeline:
/// - With `--fisher` (`fisher == true`), `var2vcf_*.pl` consumes Java's
///   `Utils.getRoundedValueToPrint` strings directly without an intermediate R step
///   (preserving trailing `.` on fractional values that round to integers, e.g. `"2."`).
/// - Without `--fisher` (`fisher == false`), `teststrandbias.R` / `testsomatic.R`
///   parses numeric columns and `write.table` strips both trailing zeros and trailing `.`.
#[inline]
fn push_vcf_rounded(buf: &mut String, fisher: bool, value: f64, digits: i32) {
    if fisher {
        push_java_rounded(buf, value, digits);
    } else {
        push_stripped(buf, value, digits);
    }
}

/// Load contig lengths from reference .fai index, falling back to BAM header.
/// Ports print_contigs / print_reference logic and fixes issue #320.
pub fn load_contigs(
    reference_path: &str,
    bam_chrom_lengths: &HashMap<String, usize>,
) -> Vec<(String, usize)> {
    let fai_path = format!("{}.fai", reference_path);
    if let Ok(content) = std::fs::read_to_string(&fai_path) {
        let mut contigs = Vec::new();
        for line in content.lines() {
            if line.is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                let name = parts[0].to_string();
                if let Ok(len) = parts[1].parse::<usize>() {
                    contigs.push((name, len));
                }
            }
        }
        if !contigs.is_empty() {
            return contigs;
        }
    }

    // Fallback: use BAM chromosomes lengths, sorted in standard chromosome order.
    let mut contigs: Vec<(String, usize)> = bam_chrom_lengths
        .iter()
        .map(|(k, v)| (k.clone(), *v))
        .collect();
    contigs.sort_by(|a, b| chrom_compare(&a.0, &b.0));
    contigs
}

/// Standard chromosome ordering for `##contig` fallback headers.
fn chrom_compare(a: &str, b: &str) -> std::cmp::Ordering {
    let clean_a = a.strip_prefix("chr").unwrap_or(a);
    let clean_b = b.strip_prefix("chr").unwrap_or(b);

    let a_is_num = clean_a.chars().all(|c| c.is_ascii_digit());
    let b_is_num = clean_b.chars().all(|c| c.is_ascii_digit());

    match (a_is_num, b_is_num) {
        (true, true) => {
            let val_a: u32 = clean_a.parse().unwrap_or(0);
            let val_b: u32 = clean_b.parse().unwrap_or(0);
            val_a.cmp(&val_b)
        }
        (true, false) => std::cmp::Ordering::Less,
        (false, true) => std::cmp::Ordering::Greater,
        (false, false) => {
            let rank = |s: &str| -> u32 {
                match s {
                    "X" => 1,
                    "Y" => 2,
                    "M" | "MT" => 3,
                    _ => 4,
                }
            };
            let rank_a = rank(clean_a);
            let rank_b = rank(clean_b);
            if rank_a != rank_b {
                rank_a.cmp(&rank_b)
            } else {
                clean_a.cmp(clean_b)
            }
        }
    }
}

/// Ports `var2vcf_valid.pl:254-288` and `var2vcf_paired.pl:260-294` (`sub reorder`):
/// - Deduplicates input chromosomes preserving first-seen order (`unique_chrs`).
/// - `@chrn`: chromosomes containing at least one ASCII digit (`\d`) and NO `'_'`,
///   paired with their digits-only string (`$t =~ s/\D//g`), stable-sorted ascending
///   by exact non-negative integer value of the digits string (without `u64` overflow).
/// - `@nonchrn`: chromosomes without digits or containing `'_'`, excluding
///   `"X"`, `"chrX"`, `"Y"`, `"chrY"`, `"MT"`, `"chrM"` (note `"M"` is NOT excluded).
/// - Sex and mitochondrial chromosomes appended via Perl's `if ... elsif` chain
///   (replicating upstream v1.8.3's `chrX`/`chrY`/`chrM` drop when both unprefixed
///   and `chr`-prefixed names are present; tracked as `--fix var2vcf-chrx-drop`):
///   - `X` else `chrX`
///   - `Y` else `chrY`
///   - `MT` else `chrM`
/// - Followed by `@nonchrn` in first-seen order.
pub fn perl_reorder_chromosomes<S: AsRef<str>>(chrs: &[S]) -> Vec<String> {
    let mut seen = std::collections::HashSet::with_capacity(chrs.len());
    let mut unique_chrs = Vec::with_capacity(chrs.len());
    for c in chrs {
        let s = c.as_ref();
        if seen.insert(s) {
            unique_chrs.push(s);
        }
    }

    let mut chrn: Vec<(String, &str)> = Vec::new();
    let mut nonchrn: Vec<&str> = Vec::new();
    for &c in &unique_chrs {
        if c.bytes().any(|b| b.is_ascii_digit()) && !c.contains('_') {
            let digits: String = c.chars().filter(|ch| ch.is_ascii_digit()).collect();
            chrn.push((digits, c));
        } else if !matches!(c, "X" | "chrX" | "Y" | "chrY" | "MT" | "chrM") {
            nonchrn.push(c);
        }
    }

    chrn.sort_by(|(da, _), (db, _)| {
        let na = da.trim_start_matches('0');
        let nb = db.trim_start_matches('0');
        na.len().cmp(&nb.len()).then_with(|| na.cmp(nb))
    });

    let mut order: Vec<String> = chrn.into_iter().map(|(_, c)| c.to_string()).collect();
    if seen.contains("X") {
        order.push("X".to_string());
    } else if seen.contains("chrX") {
        order.push("chrX".to_string());
    }
    if seen.contains("Y") {
        order.push("Y".to_string());
    } else if seen.contains("chrY") {
        order.push("chrY".to_string());
    }
    if seen.contains("MT") {
        order.push("MT".to_string());
    } else if seen.contains("chrM") {
        order.push("chrM".to_string());
    }
    for c in nonchrn {
        order.push(c.to_string());
    }
    order
}

/// Ports `var2vcf_valid.pl:21` `my $Freq = $opt_f ? $opt_f : 0.02;`:
/// in Perl, `-f 0` (`"0"`) is falsy and falls back to `0.02`, whereas non-`"0"`
/// zero strings like `-f 0.0` are truthy in Perl boolean context and keep `0.0`.
#[inline]
fn valid_freq_threshold(config: &Config) -> f64 {
    if config.af_threshold == 0.0 && !config.af_threshold_truthy_zero {
        0.02
    } else {
        config.af_threshold
    }
}

/// Effective SNV cluster distance `$opt_c` (`var2vcf_valid.pl:29`, `var2vcf_paired.pl:31`),
/// defaulting to `0` when `--vcf-cluster-dist` is unset.
#[inline]
pub fn effective_vcf_cluster_dist(config: &Config) -> i32 {
    config
        .vcf_cluster_dist
        .as_ref()
        .map_or(0, |t| t.value.max(0))
}

#[inline]
fn sanitize_nonneg_i32(val: i32, fallback: i32) -> i32 {
    if val < 0 {
        fallback
    } else {
        val
    }
}

#[inline]
fn effective_gt_freq(config: &Config) -> f64 {
    match config.vcf_gt_freq {
        Some(v) if v.is_finite() && (0.0..=1.0).contains(&v) => {
            if v == 0.0 {
                0.0
            } else {
                v
            }
        }
        _ => 0.2,
    }
}

#[inline]
fn sanitize_nonneg_f64(val: f64, fallback: f64) -> f64 {
    if !val.is_finite() || val < 0.0 {
        fallback
    } else if val == 0.0 {
        0.0
    } else {
        val
    }
}

/// Perl's default numeric stringification when interpolating an NV into a
/// string: C `%.15g` semantics (15 significant digits, trailing zeros and dot
/// trimmed; `%e` form with a signed >=2-digit exponent outside `[1e-4, 1e15)`),
/// written directly to `out` using stack buffers (`StackBuf<64>`).
fn write_perl_num_str<W: std::fmt::Write>(out: &mut W, x: f64) {
    if x.is_nan() {
        let _ = out.write_str("NaN");
        return;
    }
    if x.is_infinite() {
        let _ = out.write_str(if x < 0.0 { "-Inf" } else { "Inf" });
        return;
    }
    let neg = x < 0.0;
    let a = x.abs();
    if a == 0.0 {
        let _ = out.write_str("0");
        return;
    }
    if neg {
        let _ = out.write_char('-');
    }
    let p = 15i32;
    let mut e = crate::types::StackBuf::<64>::new();
    let _ = write!(e, "{:.*e}", (p - 1) as usize, a);
    let e_str = e.as_str();
    let (mant, exp) = e_str.split_once(['e', 'E']).unwrap_or((e_str, "0"));
    let exp10: i32 = exp.parse().unwrap_or(0);
    if exp10 < -4 || exp10 >= p {
        let mant = mant.trim_end_matches('0').trim_end_matches('.');
        let _ = write!(
            out,
            "{}e{}{:02}",
            mant,
            if exp10 < 0 { '-' } else { '+' },
            exp10.abs()
        );
    } else {
        let prec = (p - 1 - exp10).max(0) as usize;
        let mut f = crate::types::StackBuf::<64>::new();
        let _ = write!(f, "{:.*}", prec, a);
        let f_str = f.as_str();
        if f_str.contains('.') {
            let trimmed = f_str.trim_end_matches('0');
            let _ = out.write_str(trimmed.strip_suffix('.').unwrap_or(trimmed));
        } else {
            let _ = out.write_str(f_str);
        }
    }
}

/// Writes a sanitized `VcfThreshold<f64>` label (`<= 64` ASCII bytes) into `out`:
/// re-validates `thr.raw` via [`crate::config::validate_vcf_nonneg_float_token`] and
/// `to_bits()` agreement so direct library `&Config` callers cannot inject newlines,
/// tabs, or unbounded strings into VCF headers/filters, falling back to `write_perl_num_str`.
fn write_vcf_threshold_label<W: std::fmt::Write>(
    out: &mut W,
    thr: &crate::config::VcfThreshold<f64>,
    fallback: f64,
) {
    let val = sanitize_nonneg_f64(thr.value, fallback);
    if let Some(raw) = thr.raw.as_deref() {
        if let Ok((norm_val, norm_raw)) = crate::config::validate_vcf_nonneg_float_token(raw) {
            if norm_val.to_bits() == val.to_bits() {
                let _ = out.write_str(norm_raw);
                return;
            }
        }
    }
    write_perl_num_str(out, val);
}

#[inline]
fn header_float_label<'a>(
    opt: Option<&'a crate::config::VcfThreshold<f64>>,
    default_str: &'static str,
    default_num: f64,
    scratch: &'a mut crate::types::StackBuf<64>,
) -> &'a str {
    match opt {
        None => default_str,
        Some(thr) => {
            *scratch = crate::types::StackBuf::<64>::new();
            write_vcf_threshold_label(scratch, thr, default_num);
            scratch.as_str()
        }
    }
}

/// Writes a sanitized `VcfThreshold<i32>` label (`<= 64` ASCII bytes) into `out`:
/// re-validates `thr.raw` via [`crate::config::validate_vcf_nonneg_i32_token`] and
/// numeric agreement so direct library `&Config` callers cannot inject newlines,
/// tabs, or unbounded strings into VCF headers/filters, falling back to decimal `{val}`.
fn write_vcf_int_threshold_label<W: std::fmt::Write>(
    out: &mut W,
    thr: &crate::config::VcfThreshold<i32>,
    fallback: i32,
) {
    let val = sanitize_nonneg_i32(thr.value, fallback);
    if let Some(raw) = thr.raw.as_deref() {
        if let Ok((norm_val, norm_raw)) = crate::config::validate_vcf_nonneg_i32_token(raw) {
            if norm_val == val {
                let _ = out.write_str(norm_raw);
                return;
            }
        }
    }
    let _ = write!(out, "{val}");
}

#[inline]
fn header_int_label<'a>(
    opt: Option<&'a crate::config::VcfThreshold<i32>>,
    default_str: &'static str,
    default_num: i32,
    scratch: &'a mut crate::types::StackBuf<64>,
) -> &'a str {
    match opt {
        None => default_str,
        Some(thr) => {
            *scratch = crate::types::StackBuf::<64>::new();
            write_vcf_int_threshold_label(scratch, thr, default_num);
            scratch.as_str()
        }
    }
}

/// Print VCF 4.2 compliant header populated with contigs.
pub fn print_vcf_header<W: Write>(
    w: &mut W,
    config: &Config,
    bam_chrom_lengths: &HashMap<String, usize>,
) -> io::Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=VarDict_v1.8.2")?;
    writeln!(w, "##reference={}", config.reference)?;

    let contigs = load_contigs(&config.reference, bam_chrom_lengths);
    for (name, len) in contigs {
        writeln!(w, "##contig=<ID={},length={}>", name, len)?;
    }

    writeln!(w, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">")?;
    writeln!(w, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant Type: SNV Insertion Deletion Complex\">")?;
    writeln!(
        w,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Chr End Position\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=VD,Number=1,Type=Integer,Description=\"Variant Depth\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=BIAS,Number=1,Type=String,Description=\"Strand Bias Info\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=REFBIAS,Number=1,Type=String,Description=\"Reference depth by strand\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=VARBIAS,Number=1,Type=String,Description=\"Variant depth by strand\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=PMEAN,Number=1,Type=Float,Description=\"The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=PSTD,Number=1,Type=Float,Description=\"Position STD in reads\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=QUAL,Number=1,Type=Float,Description=\"Mean quality score in reads\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=QSTD,Number=1,Type=Float,Description=\"Quality score STD in reads\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SBF,Number=1,Type=Float,Description=\"Strand Bias Fisher p-value\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description=\"Strand Bias Odds ratio\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mean Mapping Quality\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SN,Number=1,Type=Float,Description=\"Signal to noise\">"
    )?;
    writeln!(w, "##INFO=<ID=HIAF,Number=1,Type=Float,Description=\"Allele frequency using only high quality bases\">")?;
    writeln!(w, "##INFO=<ID=ADJAF,Number=1,Type=Float,Description=\"Adjusted AF for indels due to local realignment\">")?;
    writeln!(w, "##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description=\"No. of bases to be shifted to 3 prime for deletions due to alternative alignment\">")?;
    writeln!(
        w,
        "##INFO=<ID=MSI,Number=1,Type=Float,Description=\"MicroSatellite. > 1 indicates MSI\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=MSILEN,Number=1,Type=Float,Description=\"MicroSatellite unit length in bp\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=NM,Number=1,Type=Float,Description=\"Mean mismatches in reads\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=LSEQ,Number=1,Type=String,Description=\"5' flanking seq\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=RSEQ,Number=1,Type=String,Description=\"3' flanking seq\">"
    )?;
    writeln!(w, "##INFO=<ID=GDAMP,Number=1,Type=Integer,Description=\"No. of amplicons supporting variant\">")?;
    writeln!(w, "##INFO=<ID=TLAMP,Number=1,Type=Integer,Description=\"Total of amplicons covering variant\">")?;
    writeln!(
        w,
        "##INFO=<ID=NCAMP,Number=1,Type=Integer,Description=\"No. of amplicons don't work\">"
    )?;
    writeln!(w, "##INFO=<ID=AMPFLAG,Number=1,Type=Integer,Description=\"Top variant in amplicons don't match\">")?;
    writeln!(
        w,
        "##INFO=<ID=HICNT,Number=1,Type=Integer,Description=\"High quality variant reads\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=HICOV,Number=1,Type=Integer,Description=\"High quality total reads\">"
    )?;
    writeln!(w, "##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description=\"No. of split reads supporting SV\">")?;
    writeln!(
        w,
        "##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description=\"No. of pairs supporting SV\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type: INV DUP DEL INS FUS\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"The length of SV in bp\">"
    )?;
    writeln!(
        w,
        "##INFO=<ID=DUPRATE,Number=1,Type=Float,Description=\"Duplication rate in fraction\">"
    )?;

    let mut s_q = crate::types::StackBuf::<64>::new();
    let mut s_mapq = crate::types::StackBuf::<64>::new();
    let mut s_p = crate::types::StackBuf::<64>::new();
    let mut s_sn = crate::types::StackBuf::<64>::new();
    let mut s_d = crate::types::StackBuf::<64>::new();
    let mut s_v = crate::types::StackBuf::<64>::new();
    let mut s_msi = crate::types::StackBuf::<64>::new();
    let mut s_nm = crate::types::StackBuf::<64>::new();
    let mut s_c = crate::types::StackBuf::<64>::new();
    let mut s_af = crate::types::StackBuf::<64>::new();

    let qmean = header_float_label(config.vcf_qmean.as_ref(), "22.5", 22.5, &mut s_q);
    let map_qual_val = header_float_label(config.vcf_mapq.as_ref(), "10", 10.0, &mut s_mapq);
    let pmean_val = header_float_label(config.vcf_pmean.as_ref(), "8", 8.0, &mut s_p);
    let sn_val = header_float_label(config.vcf_sn.as_ref(), "1.5", 1.5, &mut s_sn);
    let total_depth = header_int_label(config.vcf_min_depth.as_ref(), "3", 3, &mut s_d);
    let var_depth = header_int_label(config.vcf_min_vd.as_ref(), "2", 2, &mut s_v);
    let opt_i = header_float_label(config.vcf_msi.as_ref(), "12", 12.0, &mut s_msi);
    let opt_m = header_float_label(config.vcf_max_nm.as_ref(), "5.25", 5.25, &mut s_nm);
    let opt_c = header_int_label(config.vcf_cluster_dist.as_ref(), "0", 0, &mut s_c);

    writeln!(
        w,
        "##FILTER=<ID=q{},Description=\"Mean Base Quality Below {}\">",
        qmean, qmean
    )?;
    writeln!(
        w,
        "##FILTER=<ID=Q{},Description=\"Mean Mapping Quality Below {}\">",
        map_qual_val, map_qual_val
    )?;
    writeln!(
        w,
        "##FILTER=<ID=p{},Description=\"Mean Position in Reads Less than {}\">",
        pmean_val, pmean_val
    )?;
    writeln!(
        w,
        "##FILTER=<ID=SN{},Description=\"Signal to Noise Less than {}\">",
        sn_val, sn_val
    )?;
    writeln!(w, "##FILTER=<ID=Bias,Description=\"Strand Bias\">")?;
    writeln!(
        w,
        "##FILTER=<ID=pSTD,Description=\"Position in Reads has STD of 0\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=d{},Description=\"Total Depth < {}\">",
        total_depth, total_depth
    )?;
    writeln!(
        w,
        "##FILTER=<ID=v{},Description=\"Var Depth < {}\">",
        var_depth, var_depth
    )?;
    if let Some(ref thr) = config.vcf_af {
        write_vcf_threshold_label(&mut s_af, thr, valid_freq_threshold(config));
        let af_lbl = s_af.as_str();
        writeln!(
            w,
            "##FILTER=<ID=f{},Description=\"Allele frequency < {}\">",
            af_lbl, af_lbl
        )?;
    } else {
        let freq_val = valid_freq_threshold(config);
        writeln!(
            w,
            "##FILTER=<ID=f{},Description=\"Allele frequency < {}\">",
            freq_val, freq_val
        )?;
    }
    writeln!(w, "##FILTER=<ID=MSI{},Description=\"Variant in MSI region with {} non-monomer MSI or 13 monomer MSI\">", opt_i, opt_i)?;
    writeln!(w, "##FILTER=<ID=NM{},Description=\"Mean mismatches in reads >= {}, thus likely false positive\">", opt_m, opt_m)?;
    writeln!(w, "##FILTER=<ID=InGap,Description=\"The variant is in the deletion gap, thus likely false positive\">")?;
    writeln!(
        w,
        "##FILTER=<ID=InIns,Description=\"The variant is adjacent to an insertion variant\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=Cluster{}bp,Description=\"Two variants are within {} bp\">",
        opt_c, opt_c
    )?;
    writeln!(
        w,
        "##FILTER=<ID=LongMSI,Description=\"The somatic variant is flanked by long A/T (>=14)\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=AMPBIAS,Description=\"Indicate the variant has amplicon bias.\">"
    )?;

    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=VD,Number=1,Type=Integer,Description=\"Variant Depth\">"
    )?;
    writeln!(w, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">")?;
    writeln!(
        w,
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=RD,Number=2,Type=Integer,Description=\"Reference forward, reverse reads\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=ALD,Number=2,Type=Integer,Description=\"Variant forward, reverse reads\">"
    )?;

    let sample_nowhitespace = config.sample.replace(char::is_whitespace, "_");
    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        sample_nowhitespace
    )?;

    Ok(())
}

/// Ports var2vcf_valid.pl:141 / var2vcf_paired.pl:136 — `next unless ($ref)`:
/// a .var row whose REF-allele cell is falsy in Perl (empty, or the literal
/// "0") is skipped and never becomes a VCF record. Reachable via the -p pileup
/// null records (SimplePostProcessModule.java:73-82), which print empty
/// ref/allele cells to the .var table.
pub fn perl_row_skipped(ref_allele: &str) -> bool {
    ref_allele.is_empty() || ref_allele == "0"
}

#[derive(Clone, Copy)]
enum FilterTag<'a> {
    Text(&'a str),
    MinAf(f64),
}

struct VcfFloatTagScratch {
    slots: [crate::types::StackBuf<80>; 8],
    int_slots: [crate::types::StackBuf<80>; 3],
}

impl VcfFloatTagScratch {
    #[inline]
    fn new() -> Self {
        Self {
            slots: std::array::from_fn(|_| crate::types::StackBuf::<80>::new()),
            int_slots: std::array::from_fn(|_| crate::types::StackBuf::<80>::new()),
        }
    }
}

#[inline]
fn format_float_filter_tag<'a>(
    slot: &'a mut crate::types::StackBuf<80>,
    prefix: &str,
    thr: &crate::config::VcfThreshold<f64>,
    fallback: f64,
) -> FilterTag<'a> {
    *slot = crate::types::StackBuf::<80>::new();
    let _ = slot.write_str(prefix);
    write_vcf_threshold_label(slot, thr, fallback);
    FilterTag::Text(slot.as_str())
}

#[inline]
fn format_int_threshold_filter_tag<'a>(
    slot: &'a mut crate::types::StackBuf<80>,
    prefix: &str,
    thr: &crate::config::VcfThreshold<i32>,
    fallback: i32,
    suffix: &str,
) -> FilterTag<'a> {
    *slot = crate::types::StackBuf::<80>::new();
    let _ = slot.write_str(prefix);
    write_vcf_int_threshold_label(slot, thr, fallback);
    let _ = slot.write_str(suffix);
    FilterTag::Text(slot.as_str())
}

struct StackFilterList<'a> {
    tags: [FilterTag<'a>; 16],
    len: usize,
}

impl<'a> StackFilterList<'a> {
    #[inline]
    const fn new() -> Self {
        Self {
            tags: [FilterTag::Text(""); 16],
            len: 0,
        }
    }

    #[inline]
    fn push(&mut self, tag: FilterTag<'a>) {
        if let Some(slot) = self.tags.get_mut(self.len) {
            *slot = tag;
            self.len += 1;
        }
    }

    #[inline]
    fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    fn as_slice(&self) -> &[FilterTag<'a>] {
        self.tags.get(..self.len).unwrap_or(&[])
    }

    fn write_vcf_filter_field(&self, buf: &mut String) {
        let slice = self.as_slice();
        if slice.is_empty() {
            buf.push_str("PASS");
            return;
        }
        for (i, tag) in slice.iter().enumerate() {
            if i > 0 {
                buf.push(';');
            }
            match *tag {
                FilterTag::Text(s) => buf.push_str(s),
                FilterTag::MinAf(af) => {
                    let _ = write!(buf, "f{af}");
                }
            }
        }
    }

    fn to_vec_strings(&self) -> Vec<String> {
        let slice = self.as_slice();
        let mut out = Vec::with_capacity(slice.len());
        for tag in slice {
            match *tag {
                FilterTag::Text(s) => out.push(s.to_string()),
                FilterTag::MinAf(af) => out.push(format!("f{af}")),
            }
        }
        out
    }
}

/// Parse `(splitreads, spanpairs)` from the non-amplicon `.var` SV string
/// (`var2vcf_valid.pl:198`: `($splitreads, $spanpairs, $cluster) = (split(/-/, $tamp)) if (defined($tamp) && $tamp =~ /-/)`).
#[inline]
pub(crate) fn parse_sv_counts(sv: &str) -> (i32, i32) {
    if !sv.contains('-') {
        return (0, 0);
    }
    let mut parts = sv.split('-');
    let splitreads = parts
        .next()
        .and_then(|p| p.parse::<i32>().ok())
        .unwrap_or(0);
    let spanpairs = parts
        .next()
        .and_then(|p| p.parse::<i32>().ok())
        .unwrap_or(0);
    (splitreads, spanpairs)
}

struct ValidThresholds {
    min_depth: i32,
    min_vd: i32,
    fallback_freq: f64,
    freq: f64,
    pmean: f64,
    qmean: f64,
    mapq: f64,
    sn: f64,
    max_nm: f64,
    msi: f64,
    gt_freq: f64,
    opt_c: i32,
}

impl ValidThresholds {
    #[inline]
    fn from_config(config: &Config) -> Self {
        let fallback_freq = valid_freq_threshold(config);
        let freq = config.vcf_af.as_ref().map_or(fallback_freq, |t| {
            sanitize_nonneg_f64(t.value, fallback_freq)
        });
        Self {
            min_depth: config
                .vcf_min_depth
                .as_ref()
                .map_or(3, |t| sanitize_nonneg_i32(t.value, 3)),
            min_vd: config
                .vcf_min_vd
                .as_ref()
                .map_or(2, |t| sanitize_nonneg_i32(t.value, 2)),
            fallback_freq,
            freq,
            pmean: config
                .vcf_pmean
                .as_ref()
                .map_or(8.0, |t| sanitize_nonneg_f64(t.value, 8.0)),
            qmean: config
                .vcf_qmean
                .as_ref()
                .map_or(22.5, |t| sanitize_nonneg_f64(t.value, 22.5)),
            mapq: config
                .vcf_mapq
                .as_ref()
                .map_or(10.0, |t| sanitize_nonneg_f64(t.value, 10.0)),
            sn: config
                .vcf_sn
                .as_ref()
                .map_or(1.5, |t| sanitize_nonneg_f64(t.value, 1.5)),
            max_nm: config
                .vcf_max_nm
                .as_ref()
                .map_or(5.25, |t| sanitize_nonneg_f64(t.value, 5.25)),
            msi: config
                .vcf_msi
                .as_ref()
                .map_or(12.0, |t| sanitize_nonneg_f64(t.value, 12.0)),
            gt_freq: effective_gt_freq(config),
            opt_c: effective_vcf_cluster_dist(config),
        }
    }
}

#[inline]
fn has_non_fisher_valid_filter(
    variant: &Variant,
    config: &Config,
    thr: &ValidThresholds,
    amp: Option<(i32, i32, i32, bool)>,
) -> bool {
    let dp = variant.total_pos_cov;
    let hicnt = variant.hicnt;
    let af = crate::types::round_half_even(variant.frequency, 4);
    let pmean = crate::types::round_half_even(variant.mean_position, 1);
    let qmean = crate::types::round_half_even(variant.mean_quality, 1);
    let mapq = crate::types::round_half_even(variant.mean_mapping_quality, 1);
    let nm = crate::types::round_half_even(variant.nm, 1).max(0.0);
    let hiaf = crate::types::round_half_even(variant.hi_freq, 4);
    let sn = crate::types::round_half_even(variant.hq_lq_ratio, 3);
    let msi = crate::types::round_half_even(variant.msi, 3);
    let hi_check = (hicnt as f64 * hiaf) >= 0.5;

    if (dp < thr.min_depth || hicnt < thr.min_vd) && !hi_check {
        return true;
    }
    if af < thr.freq
        || pmean < thr.pmean
        || (config.vcf_pstd_filter && !variant.pstd && amp.is_none() && af < 0.35)
        || qmean < thr.qmean
        || (mapq < thr.mapq && af < 0.8)
        || sn < thr.sn
        || nm > thr.max_nm
    {
        return true;
    }

    let len_diff = variant.ref_allele.len().abs_diff(variant.var_allele.len());
    if usize::try_from(variant.msint).ok() == Some(len_diff) {
        if (msi > thr.msi && variant.msint > 1 && af < 0.2)
            || (msi >= 13.0 && variant.msint == 1 && af <= 0.275)
            || (hiaf <= 0.275 && msi >= 13.0)
            || (hiaf <= 0.2 && msi >= 8.0 && variant.msint > 1)
        {
            return true;
        }
    }
    false
}

#[inline]
fn bias_non_fisher_precond(variant: &Variant) -> bool {
    let hiaf = crate::types::round_half_even(variant.hi_freq, 4);
    hiaf < 0.25
        && variant.bias == "2;1"
        && (i64::from(variant.end_pos) - i64::from(variant.start_pos)) < 100
}

#[inline]
fn bias_fires_with_fisher(sbf: f64, oddratio_str: &str) -> bool {
    let sbf_r = crate::types::round_half_even(sbf, 5);
    if sbf_r >= 0.01 {
        return false;
    }
    let transformed_or = if oddratio_str == "Inf" {
        0.0
    } else {
        let raw_or = oddratio_str.parse::<f64>().unwrap_or(0.0);
        if raw_or < 1.0 && raw_or > 0.0 {
            1.0 / raw_or
        } else {
            raw_or
        }
    };
    transformed_or > 5.0 || transformed_or == 0.0
}

#[inline]
fn compute_ampbias(amp: Option<(i32, i32, i32, bool)>) -> bool {
    amp.is_some_and(|(gdamp, tlamp, ncamp, ampflag)| {
        i64::from(gdamp) < i64::from(tlamp) - i64::from(ncamp) || ampflag
    })
}

/// Collect `var2vcf_valid.pl` filter tags into `filters` backed by `tag_scratch`.
/// When `pvs` is `Some(pvs_val)`, evaluates the standalone `Cluster{opt_c}bp` check
/// (`var2vcf_valid.pl:132,181` — without any `pvs > 0` guard); when `None`, omits
/// `Cluster` so [`SingleRecordStreamer`] can evaluate it on the chromosome-sorted stream.
/// Returns `(pre_cluster_non_amp_empty, has_ampbias)`.
#[allow(clippy::too_many_arguments)]
fn collect_valid_filters_into<'a>(
    filters: &mut StackFilterList<'a>,
    tag_scratch: &'a mut VcfFloatTagScratch,
    variant: &Variant,
    config: &Config,
    thr: &ValidThresholds,
    pvs: Option<i32>,
    amp: Option<(i32, i32, i32, bool)>,
    sbf: f64,
    oddratio_str: &str,
) -> (bool, bool) {
    let [s_f, s_p, s_q, s_mapq, s_sn, s_nm, s_msi, _s_pval] = &mut tag_scratch.slots;
    let [s_d, s_v, s_cluster] = &mut tag_scratch.int_slots;

    let dp = variant.total_pos_cov;
    let hicnt = variant.hicnt;
    let af = crate::types::round_half_even(variant.frequency, 4);
    let pmean = crate::types::round_half_even(variant.mean_position, 1);
    let qmean = crate::types::round_half_even(variant.mean_quality, 1);
    let mapq = crate::types::round_half_even(variant.mean_mapping_quality, 1);
    let nm = crate::types::round_half_even(variant.nm, 1).max(0.0);
    let hiaf = crate::types::round_half_even(variant.hi_freq, 4);
    let sn = crate::types::round_half_even(variant.hq_lq_ratio, 3);
    let msi = crate::types::round_half_even(variant.msi, 3);
    let hi_check = (hicnt as f64 * hiaf) >= 0.5;

    if dp < thr.min_depth && !hi_check {
        let tag = match config.vcf_min_depth.as_ref() {
            Some(t) => format_int_threshold_filter_tag(s_d, "d", t, 3, ""),
            None => FilterTag::Text("d3"),
        };
        filters.push(tag);
    }

    if hicnt < thr.min_vd && !hi_check {
        let tag = match config.vcf_min_vd.as_ref() {
            Some(t) => format_int_threshold_filter_tag(s_v, "v", t, 2, ""),
            None => FilterTag::Text("v2"),
        };
        filters.push(tag);
    }

    if af < thr.freq {
        if let Some(ref af_thr) = config.vcf_af {
            filters.push(format_float_filter_tag(s_f, "f", af_thr, thr.fallback_freq));
        } else {
            filters.push(FilterTag::MinAf(thr.freq));
        }
    }

    if pmean < thr.pmean {
        let tag = match config.vcf_pmean.as_ref() {
            Some(t) => format_float_filter_tag(s_p, "p", t, 8.0),
            None => FilterTag::Text("p8"),
        };
        filters.push(tag);
    }

    if config.vcf_pstd_filter && !variant.pstd && amp.is_none() && af < 0.35 {
        filters.push(FilterTag::Text("pSTD"));
    }

    if qmean < thr.qmean {
        let tag = match config.vcf_qmean.as_ref() {
            Some(t) => format_float_filter_tag(s_q, "q", t, 22.5),
            None => FilterTag::Text("q22.5"),
        };
        filters.push(tag);
    }

    if mapq < thr.mapq && af < 0.8 {
        let tag = match config.vcf_mapq.as_ref() {
            Some(t) => format_float_filter_tag(s_mapq, "Q", t, 10.0),
            None => FilterTag::Text("Q10"),
        };
        filters.push(tag);
    }

    if sn < thr.sn {
        let tag = match config.vcf_sn.as_ref() {
            Some(t) => format_float_filter_tag(s_sn, "SN", t, 1.5),
            None => FilterTag::Text("SN1.5"),
        };
        filters.push(tag);
    }

    if nm > thr.max_nm {
        let tag = match config.vcf_max_nm.as_ref() {
            Some(t) => format_float_filter_tag(s_nm, "NM", t, 5.25),
            None => FilterTag::Text("NM5.25"),
        };
        filters.push(tag);
    }

    let len_diff = variant.ref_allele.len().abs_diff(variant.var_allele.len());
    let len_matches_msint = usize::try_from(variant.msint).ok() == Some(len_diff);
    let msi_check1 = msi > thr.msi && variant.msint > 1 && af < 0.2 && len_matches_msint;
    let msi_check2 = msi >= 13.0 && variant.msint == 1 && af <= 0.275 && len_matches_msint;
    if msi_check1 || msi_check2 {
        let tag = match config.vcf_msi.as_ref() {
            Some(t) => format_float_filter_tag(s_msi, "MSI", t, 12.0),
            None => FilterTag::Text("MSI12"),
        };
        filters.push(tag);
    }

    if bias_non_fisher_precond(variant) && bias_fires_with_fisher(sbf, oddratio_str) {
        filters.push(FilterTag::Text("Bias"));
    }

    if len_matches_msint
        && ((hiaf <= 0.275 && msi >= 13.0) || (hiaf <= 0.2 && msi >= 8.0 && variant.msint > 1))
    {
        filters.push(FilterTag::Text("LongMSI"));
    }

    let pre_cluster_non_amp_empty = filters.is_empty();

    if let Some(pvs_val) = pvs {
        if variant.var_type == "SNV"
            && pre_cluster_non_amp_empty
            && i64::from(variant.start_pos) - i64::from(pvs_val) < i64::from(thr.opt_c)
        {
            let tag = match config.vcf_cluster_dist.as_ref() {
                Some(t) => format_int_threshold_filter_tag(s_cluster, "Cluster", t, 0, "bp"),
                None => FilterTag::Text("Cluster0bp"),
            };
            filters.push(tag);
        }
    }

    let has_ampbias = compute_ampbias(amp);
    if has_ampbias {
        filters.push(FilterTag::Text("AMPBIAS"));
    }

    (pre_cluster_non_amp_empty, has_ampbias)
}

/// Run var2vcf_valid.pl filter checks on a finalized Variant.
/// Ports var2vcf_valid.pl filter logic (T098-T104).
pub fn filter_variant(
    variant: &Variant,
    config: &Config,
    pvs: i32,
    // Amplicon columns 37-40 of the 40-col table (GDAMP, TLAMP, NCAMP, AMPFLAG);
    // None = simple mode ($isamp false, var2vcf_valid.pl:145).
    amp: Option<(i32, i32, i32, bool)>,
    sbf: f64,
    oddratio_str: &str,
) -> Vec<String> {
    let thr = ValidThresholds::from_config(config);
    let mut tag_scratch = VcfFloatTagScratch::new();
    let mut filters = StackFilterList::new();
    collect_valid_filters_into(
        &mut filters,
        &mut tag_scratch,
        variant,
        config,
        &thr,
        Some(pvs),
        amp,
        sbf,
        oddratio_str,
    );
    filters.to_vec_strings()
}

/// Java's RAW table cell for hifreq (trailing zeros kept): ports
/// SimpleOutputVariant.java:127-129 / :211 — `hifreq == 0 ? "0" :
/// new DecimalFormat("0.0000").format(hifreq)` — NOT getRoundedValueToPrint.
/// The cell survives to Perl verbatim: teststrandbias.R's 22-entry colClasses
/// (:14) is recycled over the input columns, landing "character" on input col 23
/// (hifreq), so R echoes it unchanged and var2vcf_valid.pl:24/:236 interpolates
/// it as-is into `HIAF=`. (The PAIRED writer must stay stripped: testsomatic.R's
/// 48-entry pattern leaves hiaf1/hiaf2 numeric, so R re-formats them.)
#[cfg(test)]
fn hiaf_cell(hi_freq: f64) -> String {
    if hi_freq == 0.0 {
        "0".to_string()
    } else {
        format!("{:.4}", crate::types::round_half_even(hi_freq, 4))
    }
}

#[cfg(test)]
fn perl_num_str(x: f64) -> String {
    let mut out = String::new();
    write_perl_num_str(&mut out, x);
    out
}

#[allow(clippy::too_many_arguments)]
fn format_vcf_record_with_filter<F: FnOnce(&mut String)>(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    config: &Config,
    freq_threshold: f64,
    gt_freq: f64,
    sbf: f64,
    oddratio_str: &str,
    // Amplicon columns 37-40 (GDAMP, TLAMP, NCAMP, AMPFLAG); None = simple mode.
    amp: Option<(i32, i32, i32, bool)>,
    write_filter: F,
) -> (usize, usize) {
    buf.clear();
    // Golden single-sample VCF p50 is 398 B, max 551 B: 512 covers the common
    // record, long alleles regrow once.
    buf.reserve(512);
    let var_type = if variant.var_type.is_empty() {
        "REF"
    } else {
        variant.var_type.as_str()
    };

    // ref==alt (-p pileup reference calls) become ALT ".", GT "0/0".
    let alt_is_dot = variant.ref_allele == variant.var_allele;
    let gt = if alt_is_dot {
        "0/0"
    } else {
        let af = crate::types::round_half_even(variant.frequency, 4);
        if 1.0 - af < gt_freq {
            "1/1"
        } else if af >= 0.5 {
            "1/0"
        } else if af >= freq_threshold {
            "0/1"
        } else {
            "0/0"
        }
    };

    let qual_val = if variant.pos_coverage <= 1 {
        0
    } else {
        // ports var2vcf_valid.pl:190 `int(log($vd)/log(2) * $qual)` — mirrored
        // literally as ln/LN_2 so the computation is bit-exact against Perl's
        // log/log2 division form (f64::log2 can differ by 1 ulp). $qual is
        // table col 19, printed by Java at HALF_EVEN 1dp (SimpleOutputVariant.java:155
        // getRoundedValueToPrint("0.0", …) / :206 DecimalFormat("0.0")), so the
        // multiplicand is the table-rounded value, not the raw mean.
        ((variant.pos_coverage as f64).ln() / std::f64::consts::LN_2
            * crate::types::round_half_even(variant.mean_quality, 1)) as i32
    };

    buf.push_str(&region.chr);
    buf.push('\t');
    let _ = write!(buf, "{}\t", variant.start_pos);
    buf.push_str(".\t");
    buf.push_str(&variant.ref_allele);
    buf.push('\t');
    if alt_is_dot {
        buf.push_str(".\t");
    } else {
        buf.push_str(&variant.var_allele);
        buf.push('\t');
    }
    let _ = write!(buf, "{}\t", qual_val);
    let filter_start = buf.len();
    write_filter(buf);
    let filter_end = buf.len();
    buf.push('\t');

    // INFO (field order mirrors var2vcf_valid.pl:236; values contain no
    // spaces so the perl post-strip is a no-op here).
    buf.push_str("SAMPLE=");
    for c in config.sample.chars() {
        buf.push(if c.is_whitespace() { '_' } else { c });
    }
    buf.push_str(";TYPE=");
    buf.push_str(var_type);
    let _ = write!(buf, ";DP={}", variant.total_pos_cov);
    if !config.omit_end_tag {
        let _ = write!(buf, ";END={}", variant.end_pos);
    }
    let _ = write!(buf, ";VD={}", variant.pos_coverage);
    buf.push_str(";AF=");
    push_vcf_rounded(buf, config.fisher, variant.frequency, 4);
    buf.push_str(";BIAS=");
    for c in variant.bias.chars() {
        buf.push(if c == ';' { ':' } else { c });
    }
    let _ = write!(
        buf,
        ";REFBIAS={}:{};VARBIAS={}:{}",
        variant.ref_fwd_cov, variant.ref_rev_cov, variant.fwd, variant.rev
    );
    buf.push_str(";PMEAN=");
    push_vcf_rounded(buf, config.fisher, variant.mean_position, 1);
    let _ = write!(buf, ";PSTD={}", if variant.pstd { 1 } else { 0 });
    buf.push_str(";QUAL=");
    push_vcf_rounded(buf, config.fisher, variant.mean_quality, 1);
    let _ = write!(buf, ";QSTD={}", if variant.qstd { 1 } else { 0 });
    buf.push_str(";SBF=");
    push_vcf_rounded(buf, config.fisher, sbf, 5);
    // ports var2vcf_valid.pl:147-151 + :236 interpolation: "Inf" → 0; a value in
    // (0,1) is replaced by its reciprocal and stringified with Perl's default
    // numeric format (%.15g); any other table string passes through unchanged.
    buf.push_str(";ODDRATIO=");
    if oddratio_str == "Inf" {
        buf.push('0');
    } else {
        match oddratio_str.parse::<f64>() {
            Ok(x) if x < 1.0 && x > 0.0 => write_perl_num_str(buf, 1.0 / x),
            _ => buf.push_str(oddratio_str),
        }
    }
    buf.push_str(";MQ=");
    push_vcf_rounded(buf, config.fisher, variant.mean_mapping_quality, 1);
    buf.push_str(";SN=");
    push_vcf_rounded(buf, config.fisher, variant.hq_lq_ratio, 3);
    // Java's RAW table cell for hifreq (trailing zeros kept): ports
    // SimpleOutputVariant.java:127-129 / :211 — `hifreq == 0 ? "0" :
    // new DecimalFormat("0.0000").format(hifreq)` — NOT getRoundedValueToPrint.
    // (The PAIRED writer stays stripped: testsomatic.R's 48-entry pattern
    // leaves hiaf1/hiaf2 numeric, so R re-formats them.)
    buf.push_str(";HIAF=");
    push_keep(buf, variant.hi_freq, 4);
    buf.push_str(";ADJAF=");
    push_vcf_rounded(buf, config.fisher, variant.extra_freq, 4);
    let _ = write!(buf, ";SHIFT3={}", variant.shift3);
    buf.push_str(";MSI=");
    push_vcf_rounded(buf, config.fisher, variant.msi, 3);
    let _ = write!(buf, ";MSILEN={}", variant.msint);
    buf.push_str(";NM=");
    if variant.nm > 0.0 {
        push_keep(buf, variant.nm, 1);
    } else {
        buf.push('0');
    }
    let _ = write!(buf, ";HICNT={};HICOV={}", variant.hicnt, variant.hicov);
    buf.push_str(";LSEQ=");
    if variant.left_seq.is_empty() {
        buf.push('0');
    } else {
        buf.push_str(&variant.left_seq);
    }
    buf.push_str(";RSEQ=");
    if variant.right_seq.is_empty() {
        buf.push('0');
    } else {
        buf.push_str(&variant.right_seq);
    }
    // After RSEQ: amplicon mode emits the GDAMP/TLAMP/NCAMP/AMPFLAG tags (:233),
    // simple mode DUPRATE (its col 37 is Duprate) — never both.
    if let Some((gdamp, tlamp, ncamp, ampflag)) = amp {
        let _ = write!(
            buf,
            ";GDAMP={};TLAMP={};NCAMP={};AMPFLAG={}",
            gdamp,
            tlamp,
            ncamp,
            i32::from(ampflag)
        );
    } else {
        // DUPRATE cell format differs by pipeline: with --fisher Perl consumes the
        // 38-col table directly — getRoundedValueToPrint("0.00", duprate)
        // (SimpleOutputVariant.java:175). Without --fisher the 36-col cell is
        // DecimalFormat("0.0")-based (SimpleOutputVariant.java:223) and R echoes it
        // numeric-stripped (teststrandbias.R colClasses leaves input col 35
        // numeric), i.e. a 1dp value.
        buf.push_str(";DUPRATE=");
        if config.fisher {
            push_java_rounded(buf, variant.duprate, 2);
        } else {
            push_stripped(buf, variant.duprate, 1);
        }
        if matches!(config.crispr_cutting_site, Some(n) if n != 0) {
            let _ = write!(buf, ";CRISPR={}", variant.crispr);
        }
    }

    let is_amp = amp.is_some();
    if !is_amp {
        // ports var2vcf_valid.pl:197-212 — the SV block runs only when !$isamp
        let (splitreads, spanpairs) = parse_sv_counts(&variant.sv);

        if !alt_is_dot && variant.var_allele.contains('<') {
            // ports var2vcf_valid.pl:199,206-210 — `$alt =~ /</` (including `#342`
            // mid-allele `<dup-N>` when `--min-splitreads 0`); SVLEN = end - start,
            // incremented only when the ALT contains INV.
            let mut svlen = i64::from(variant.end_pos) - i64::from(variant.start_pos);
            if variant.var_allele.contains("INV") {
                svlen += 1;
            }
            buf.push_str(";SVTYPE=");
            for c in variant.var_allele.chars() {
                if c != '<' && c != '>' {
                    buf.push(c);
                }
            }
            let _ = write!(buf, ";SVLEN={}", svlen);
        }
        let _ = write!(buf, ";SPLITREAD={};SPANPAIR={}", splitreads, spanpairs);
    }

    // FORMAT + sample columns (ports var2vcf_valid.pl:240-241).
    buf.push_str("\tGT:DP:VD:AD:AF:RD:ALD\t");
    buf.push_str(gt);
    buf.push(':');
    let _ = write!(buf, "{}", variant.total_pos_cov);
    buf.push(':');
    let _ = write!(buf, "{}", variant.pos_coverage);
    buf.push(':');
    let rd = i64::from(variant.ref_fwd_cov) + i64::from(variant.ref_rev_cov);
    let vd = variant.pos_coverage;
    if vd != 0 {
        let _ = write!(buf, "{},{}", rd, vd);
    } else {
        let _ = write!(buf, "{}", rd);
    }
    buf.push(':');
    push_vcf_rounded(buf, config.fisher, variant.frequency, 4);
    buf.push(':');
    let _ = write!(
        buf,
        "{},{}:{},{}",
        variant.ref_fwd_cov, variant.ref_rev_cov, variant.fwd, variant.rev
    );
    (filter_start, filter_end)
}

/// Write a variant call as a VCF record line into `buf` (cleared first).
/// Byte-identical to [`format_vcf_record`]; the caller's buffer is reused
/// across rows (one reserve, `clear()` per variant) instead of one temp
/// `String` per INFO field plus a final `format!`.
/// Ports VCF serialization (T105).
#[allow(clippy::too_many_arguments)]
pub fn format_vcf_record_into(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    config: &Config,
    filters: &[String],
    sbf: f64,
    oddratio_str: &str,
    // Amplicon columns 37-40 (GDAMP, TLAMP, NCAMP, AMPFLAG); None = simple mode.
    amp: Option<(i32, i32, i32, bool)>,
) {
    let fallback_freq = valid_freq_threshold(config);
    let freq_threshold = config.vcf_af.as_ref().map_or(fallback_freq, |t| {
        sanitize_nonneg_f64(t.value, fallback_freq)
    });
    let gt_freq = effective_gt_freq(config);
    let _ = format_vcf_record_with_filter(
        buf,
        variant,
        region,
        config,
        freq_threshold,
        gt_freq,
        sbf,
        oddratio_str,
        amp,
        |out| {
            if filters.is_empty() {
                out.push_str("PASS");
            } else {
                for (i, f) in filters.iter().enumerate() {
                    if i > 0 {
                        out.push(';');
                    }
                    out.push_str(f);
                }
            }
        },
    );
}

/// Format a variant call to a VCF record line.
/// Ports VCF serialization (T105).
pub fn format_vcf_record(
    variant: &Variant,
    region: &Region,
    config: &Config,
    filters: &[String],
    sbf: f64,
    oddratio_str: &str,
    // Amplicon columns 37-40 (GDAMP, TLAMP, NCAMP, AMPFLAG); None = simple mode.
    amp: Option<(i32, i32, i32, bool)>,
) -> String {
    let mut buf = String::new();
    format_vcf_record_into(
        &mut buf,
        variant,
        region,
        config,
        filters,
        sbf,
        oddratio_str,
        amp,
    );
    buf
}

/// Filter and write a single-sample/amplicon VCF record line into `buf` without
/// allocating an intermediate `Vec<String>` for FILTER tags.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn format_filtered_vcf_record_into(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    config: &Config,
    pvs: i32,
    sbf: f64,
    oddratio_str: &str,
    // Amplicon columns 37-40 (GDAMP, TLAMP, NCAMP, AMPFLAG); None = simple mode.
    amp: Option<(i32, i32, i32, bool)>,
) {
    let thr = ValidThresholds::from_config(config);
    let mut tag_scratch = VcfFloatTagScratch::new();
    let mut filters = StackFilterList::new();
    collect_valid_filters_into(
        &mut filters,
        &mut tag_scratch,
        variant,
        config,
        &thr,
        Some(pvs),
        amp,
        sbf,
        oddratio_str,
    );
    let _ = format_vcf_record_with_filter(
        buf,
        variant,
        region,
        config,
        thr.freq,
        thr.gt_freq,
        sbf,
        oddratio_str,
        amp,
        |out| filters.write_vcf_filter_field(out),
    );
}

/// Filter and conditionally write a single-sample/amplicon VCF record line into `buf`.
///
/// When `config.vcf_pass_only` (`--pass-only`, `var2vcf_valid.pl:187`) is enabled and
/// any filter fires (`!filters.is_empty()`), clears `buf` and returns `false` without
/// formatting the record body.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn try_format_filtered_vcf_record_into(
    buf: &mut String,
    variant: &Variant,
    region: &Region,
    config: &Config,
    pvs: i32,
    sbf: f64,
    oddratio_str: &str,
    // Amplicon columns 37-40 (GDAMP, TLAMP, NCAMP, AMPFLAG); None = simple mode.
    amp: Option<(i32, i32, i32, bool)>,
) -> bool {
    let thr = ValidThresholds::from_config(config);
    let mut tag_scratch = VcfFloatTagScratch::new();
    let mut filters = StackFilterList::new();
    collect_valid_filters_into(
        &mut filters,
        &mut tag_scratch,
        variant,
        config,
        &thr,
        Some(pvs),
        amp,
        sbf,
        oddratio_str,
    );
    if config.vcf_pass_only && !filters.is_empty() {
        buf.clear();
        return false;
    }
    let _ = format_vcf_record_with_filter(
        buf,
        variant,
        region,
        config,
        thr.freq,
        thr.gt_freq,
        sbf,
        oddratio_str,
        amp,
        |out| filters.write_vcf_filter_field(out),
    );
    true
}

#[inline]
fn is_valid_vcf_filter_span(
    line: &str,
    filter_start: usize,
    filter_end: usize,
    expected_tabs: usize,
) -> bool {
    let Some(filt) = line.get(filter_start..filter_end) else {
        return false;
    };
    let Some(prev_idx) = filter_start.checked_sub(1) else {
        return false;
    };
    let Some(prefix) = line.get(..filter_start) else {
        return false;
    };
    !filt.is_empty()
        && !line.bytes().any(|b| matches!(b, b'\r' | b'\n'))
        && !filt.as_bytes().contains(&b'\t')
        && line.as_bytes().get(prev_idx) == Some(&b'\t')
        && line.as_bytes().get(filter_end) == Some(&b'\t')
        && prefix.bytes().filter(|&b| b == b'\t').count() == 6
        && line.bytes().filter(|&b| b == b'\t').count() == expected_tabs
}

#[inline]
fn replace_vcf_filter_span(
    line: &mut String,
    filter_start: usize,
    filter_end: &mut usize,
    new_filter: &str,
) {
    if line.get(filter_start..*filter_end).is_some() {
        line.replace_range(filter_start..*filter_end, new_filter);
        *filter_end = filter_start.saturating_add(new_filter.len());
    }
}

/// One formatted single-sample / amplicon VCF record awaiting delayed emission
/// (`CHROM\tPOS\t.\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE`, 9 tabs).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SingleVcfRecord {
    line: String,
    filter_start: usize,
    filter_end: usize,
    pub start: i32,
    pub type_snv: bool,
    pub pre_cluster_non_amp_empty: bool,
    pub has_ampbias: bool,
}

impl SingleVcfRecord {
    #[inline]
    pub fn as_str(&self) -> &str {
        &self.line
    }

    #[inline]
    pub fn filter_str(&self) -> &str {
        self.line
            .get(self.filter_start..self.filter_end)
            .unwrap_or("")
    }

    #[inline]
    pub fn filter_offsets(&self) -> (usize, usize) {
        (self.filter_start, self.filter_end)
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        self.line.capacity()
    }

    #[inline]
    pub fn is_pass(&self) -> bool {
        self.filter_str() == "PASS"
    }

    #[inline]
    pub fn is_snv(&self) -> bool {
        self.type_snv
    }

    #[inline]
    pub fn pre_cluster_empty(&self) -> bool {
        self.pre_cluster_non_amp_empty
    }

    #[inline]
    pub(crate) fn set_cluster_filter(&mut self, cluster_tag: &str) {
        replace_vcf_filter_span(
            &mut self.line,
            self.filter_start,
            &mut self.filter_end,
            cluster_tag,
        );
    }

    #[inline]
    pub fn from_parts_checked(
        line: String,
        filter_start: usize,
        filter_end: usize,
        start: i32,
        type_snv: bool,
        pre_cluster_non_amp_empty: bool,
        has_ampbias: bool,
    ) -> Option<Self> {
        if !is_valid_vcf_filter_span(&line, filter_start, filter_end, 9) {
            return None;
        }
        Some(Self {
            line,
            filter_start,
            filter_end,
            start,
            type_snv,
            pre_cluster_non_amp_empty,
            has_ampbias,
        })
    }
}

/// Worker-produced single-sample / amplicon VCF stream event consumed by
/// [`SingleRecordStreamer`] in chromosome-sorted order (`var2vcf_valid.pl:132, 181-205, 236-244`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SingleVcfEvent {
    Skip,
    AmpBiasSnvClusterProbe { start: i32 },
    SvSplitreadFlush { cluster_probe_start: Option<i32> },
    Record(SingleVcfRecord),
}

/// Build a [`SingleVcfEvent`] on a worker thread without evaluating `Cluster{opt_c}bp`
/// (which [`SingleRecordStreamer`] evaluates on the chromosome-sorted stream) and
/// skipping `FisherExact` + `uniroot` MLE under `config.vcf_pass_only` whenever any
/// non-Fisher filter already disqualifies the call.
pub fn build_single_vcf_event(
    scratch: &mut String,
    variant: &Variant,
    region: &Region,
    config: &Config,
    amp: Option<(i32, i32, i32, bool)>,
) -> SingleVcfEvent {
    if perl_row_skipped(&variant.ref_allele) {
        scratch.clear();
        return SingleVcfEvent::Skip;
    }

    let thr = ValidThresholds::from_config(config);
    let var_type = if variant.var_type.is_empty() {
        "REF"
    } else {
        variant.var_type.as_str()
    };
    let type_snv = var_type == "SNV";
    let is_dropped_sv = amp.is_none()
        && variant.var_allele.contains('<')
        && parse_sv_counts(&variant.sv).0 < config.vcf_min_splitreads;

    if config.vcf_pass_only {
        if has_non_fisher_valid_filter(variant, config, &thr, amp) {
            scratch.clear();
            return SingleVcfEvent::Skip;
        }
        let has_ampbias = compute_ampbias(amp);
        if has_ampbias {
            if thr.opt_c <= 0 || !type_snv {
                scratch.clear();
                return SingleVcfEvent::Skip;
            }
            if bias_non_fisher_precond(variant) {
                let fisher = crate::stats::FisherExact::new_with_mode(
                    variant.ref_fwd_cov,
                    variant.ref_rev_cov,
                    variant.fwd,
                    variant.rev,
                    config.fisher,
                );
                let sbf = fisher.get_p_value();
                let oddratio_str = fisher.odds_ratio_buf();
                if bias_fires_with_fisher(sbf, oddratio_str.as_str()) {
                    scratch.clear();
                    return SingleVcfEvent::Skip;
                }
            }
            scratch.clear();
            return SingleVcfEvent::AmpBiasSnvClusterProbe {
                start: variant.start_pos,
            };
        }
        if is_dropped_sv {
            if bias_non_fisher_precond(variant) {
                let fisher = crate::stats::FisherExact::new_with_mode(
                    variant.ref_fwd_cov,
                    variant.ref_rev_cov,
                    variant.fwd,
                    variant.rev,
                    config.fisher,
                );
                let sbf = fisher.get_p_value();
                let oddratio_str = fisher.odds_ratio_buf();
                if bias_fires_with_fisher(sbf, oddratio_str.as_str()) {
                    scratch.clear();
                    return SingleVcfEvent::Skip;
                }
            }
            scratch.clear();
            return SingleVcfEvent::SvSplitreadFlush {
                cluster_probe_start: type_snv.then_some(variant.start_pos),
            };
        }

        let fisher = crate::stats::FisherExact::new_with_mode(
            variant.ref_fwd_cov,
            variant.ref_rev_cov,
            variant.fwd,
            variant.rev,
            config.fisher,
        );
        let sbf = fisher.get_p_value();
        let oddratio_str = fisher.odds_ratio_buf();
        if bias_non_fisher_precond(variant) && bias_fires_with_fisher(sbf, oddratio_str.as_str()) {
            scratch.clear();
            return SingleVcfEvent::Skip;
        }
        let (filter_start, filter_end) = format_vcf_record_with_filter(
            scratch,
            variant,
            region,
            config,
            thr.freq,
            thr.gt_freq,
            sbf,
            oddratio_str.as_str(),
            amp,
            |out| out.push_str("PASS"),
        );
        let line = String::from(scratch.as_str());
        return SingleVcfEvent::Record(SingleVcfRecord {
            line,
            filter_start,
            filter_end,
            start: variant.start_pos,
            type_snv,
            pre_cluster_non_amp_empty: true,
            has_ampbias: false,
        });
    }

    // `!config.vcf_pass_only`
    if is_dropped_sv {
        let cluster_probe_start =
            if type_snv && !has_non_fisher_valid_filter(variant, config, &thr, amp) {
                let bias_fires = if bias_non_fisher_precond(variant) {
                    let fisher = crate::stats::FisherExact::new_with_mode(
                        variant.ref_fwd_cov,
                        variant.ref_rev_cov,
                        variant.fwd,
                        variant.rev,
                        config.fisher,
                    );
                    let sbf = fisher.get_p_value();
                    let oddratio_str = fisher.odds_ratio_buf();
                    bias_fires_with_fisher(sbf, oddratio_str.as_str())
                } else {
                    false
                };
                (!bias_fires).then_some(variant.start_pos)
            } else {
                None
            };
        scratch.clear();
        return SingleVcfEvent::SvSplitreadFlush {
            cluster_probe_start,
        };
    }

    let fisher = crate::stats::FisherExact::new_with_mode(
        variant.ref_fwd_cov,
        variant.ref_rev_cov,
        variant.fwd,
        variant.rev,
        config.fisher,
    );
    let sbf = fisher.get_p_value();
    let oddratio_str = fisher.odds_ratio_buf();
    let mut tag_scratch = VcfFloatTagScratch::new();
    let mut filters = StackFilterList::new();
    let (pre_cluster_non_amp_empty, has_ampbias) = collect_valid_filters_into(
        &mut filters,
        &mut tag_scratch,
        variant,
        config,
        &thr,
        None,
        amp,
        sbf,
        oddratio_str.as_str(),
    );

    let (filter_start, filter_end) = format_vcf_record_with_filter(
        scratch,
        variant,
        region,
        config,
        thr.freq,
        thr.gt_freq,
        sbf,
        oddratio_str.as_str(),
        amp,
        |out| filters.write_vcf_filter_field(out),
    );
    let line = String::from(scratch.as_str());
    SingleVcfEvent::Record(SingleVcfRecord {
        line,
        filter_start,
        filter_end,
        start: variant.start_pos,
        type_snv,
        pre_cluster_non_amp_empty,
        has_ampbias,
    })
}

#[inline]
fn format_cluster_tag_from_int(opt_c: i32) -> crate::types::StackBuf<80> {
    let mut buf = crate::types::StackBuf::<80>::new();
    let _ = write!(buf, "Cluster{}bp", opt_c.max(0));
    buf
}

#[inline]
fn format_cluster_tag_from_config(config: &Config) -> crate::types::StackBuf<80> {
    let mut buf = crate::types::StackBuf::<80>::new();
    match config.vcf_cluster_dist.as_ref() {
        Some(thr) => {
            let _ = buf.write_str("Cluster");
            write_vcf_int_threshold_label(&mut buf, thr, 0);
            let _ = buf.write_str("bp");
        }
        None => {
            let _ = buf.write_str("Cluster0bp");
        }
    }
    buf
}

/// Streaming delayed-print state machine for single-sample / amplicon VCF records
/// (`var2vcf_valid.pl:132, 181-205, 239, 242-244`).
#[derive(Debug)]
pub struct SingleRecordStreamer {
    pending: Option<SingleVcfRecord>,
    pvs: i32,
    opt_c: i32,
    pass_only: bool,
    cluster_tag: crate::types::StackBuf<80>,
}

impl SingleRecordStreamer {
    pub fn new(opt_c: i32) -> Self {
        Self::with_pass_only(opt_c, false)
    }

    pub fn with_pass_only(opt_c: i32, pass_only: bool) -> Self {
        Self {
            pending: None,
            pvs: 0,
            opt_c,
            pass_only,
            cluster_tag: format_cluster_tag_from_int(opt_c),
        }
    }

    pub fn from_config(config: &Config) -> Self {
        Self {
            pending: None,
            pvs: 0,
            opt_c: effective_vcf_cluster_dist(config),
            pass_only: config.vcf_pass_only,
            cluster_tag: format_cluster_tag_from_config(config),
        }
    }

    pub fn push_event<F>(&mut self, event: SingleVcfEvent, mut emit: F) -> io::Result<()>
    where
        F: FnMut(&str) -> io::Result<()>,
    {
        match event {
            SingleVcfEvent::Skip => Ok(()),
            SingleVcfEvent::AmpBiasSnvClusterProbe { start } => {
                if i64::from(start) - i64::from(self.pvs) < i64::from(self.opt_c) {
                    if let Some(prev) = self.pending.as_mut() {
                        prev.set_cluster_filter(self.cluster_tag.as_str());
                    }
                }
                Ok(())
            }
            SingleVcfEvent::SvSplitreadFlush {
                cluster_probe_start,
            } => {
                if let Some(start) = cluster_probe_start {
                    if i64::from(start) - i64::from(self.pvs) < i64::from(self.opt_c) {
                        if let Some(prev) = self.pending.as_mut() {
                            prev.set_cluster_filter(self.cluster_tag.as_str());
                        }
                        if self.pass_only {
                            return Ok(());
                        }
                    }
                }
                if let Some(prev) = self.pending.take() {
                    if !self.pass_only || prev.is_pass() {
                        emit(prev.as_str())?;
                    }
                }
                Ok(())
            }
            SingleVcfEvent::Record(mut cur) => {
                if cur.is_snv()
                    && cur.pre_cluster_empty()
                    && i64::from(cur.start) - i64::from(self.pvs) < i64::from(self.opt_c)
                {
                    if let Some(prev) = self.pending.as_mut() {
                        prev.set_cluster_filter(self.cluster_tag.as_str());
                    }
                    if cur.has_ampbias {
                        let mut cur_tag = crate::types::StackBuf::<96>::new();
                        let _ = cur_tag.write_str(self.cluster_tag.as_str());
                        let _ = cur_tag.write_str(";AMPBIAS");
                        cur.set_cluster_filter(cur_tag.as_str());
                    } else {
                        cur.set_cluster_filter(self.cluster_tag.as_str());
                    }
                }
                if self.pass_only && !cur.is_pass() {
                    return Ok(());
                }
                if let Some(prev) = self.pending.take() {
                    if !self.pass_only || prev.is_pass() {
                        emit(prev.as_str())?;
                    }
                }
                if cur.is_snv() && cur.is_pass() {
                    self.pvs = cur.start;
                }
                self.pending = Some(cur);
                Ok(())
            }
        }
    }

    pub fn finish<F>(mut self, mut emit: F) -> io::Result<()>
    where
        F: FnMut(&str) -> io::Result<()>,
    {
        if let Some(prev) = self.pending.take() {
            if !self.pass_only || prev.is_pass() {
                emit(prev.as_str())?;
            }
        }
        Ok(())
    }
}

// ---- var2vcf_paired.pl port (somatic/paired VCF, W4.5c) ----

/// Paired-writer defaults (var2vcf_paired.pl:19-31).
const PAIRED_MIN_DEPTH: i32 = 5;
const PAIRED_VAR_DEPTH: i32 = 3;
const PAIRED_QMEAN: f64 = 22.5;
const PAIRED_MQMEAN: f64 = 0.0;
const PAIRED_PMEAN: f64 = 8.0;
const PAIRED_SN: f64 = 1.5;
const PAIRED_OPT_I: f64 = 12.0;
const PAIRED_OPT_M: f64 = 5.25;

#[inline]
fn effective_paired_pval_num(config: &Config) -> f64 {
    if !config.vcf_pval.is_finite() || config.vcf_pval < 0.0 {
        0.05
    } else if config.vcf_pval == 0.0 {
        0.0
    } else {
        config.vcf_pval
    }
}

/// Format `$PVAL` (`<= 64` ASCII bytes) for `var2vcf_paired.pl:27,86,200` onto the stack:
/// uses `config.vcf_pval_raw` verbatim when valid and consistent with `config.vcf_pval`,
/// else `write_perl_num_str(out, eff_pval)`.
fn write_vcf_pval_label<W: std::fmt::Write>(out: &mut W, config: &Config) {
    let pval = effective_paired_pval_num(config);
    if let Some(raw) = config.vcf_pval_raw.as_deref() {
        if let Ok((p_norm, norm_raw)) = crate::config::validate_vcf_pval_token(raw) {
            if p_norm.to_bits() == pval.to_bits() {
                let _ = out.write_str(norm_raw);
                return;
            }
        }
    }
    if pval.to_bits() == 0.05_f64.to_bits() {
        let _ = out.write_str("0.05");
    } else {
        write_perl_num_str(out, pval);
    }
}

/// Paired VCF header. Ports var2vcf_paired.pl:51-115 (INFO :60-75,
/// FILTER :76-93, FORMAT :94-112, two-sample #CHROM :115 with `-N` split
/// :45-48 and whitespace→underscore :49).
pub fn print_paired_vcf_header<W: Write>(
    w: &mut W,
    config: &Config,
    bam_chrom_lengths: &HashMap<String, usize>,
) -> io::Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=VarDict_v1.8.2")?;
    writeln!(w, "##reference={}", config.reference)?;

    let contigs = load_contigs(&config.reference, bam_chrom_lengths);
    for (name, len) in contigs {
        writeln!(w, "##contig=<ID={},length={}>", name, len)?;
    }

    // ports var2vcf_paired.pl:60-75
    let info_lines = [
        r#"##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">"#,
        r#"##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">"#,
        r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#,
        r#"##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">"#,
        r#"##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">"#,
        r#"##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"#,
        r#"##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">"#,
        r#"##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSatellite. > 1 indicates MSI">"#,
        r#"##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MSI unit repeat length in bp">"#,
        r#"##INFO=<ID=SSF,Number=1,Type=Float,Description="P-value">"#,
        r#"##INFO=<ID=SOR,Number=1,Type=Float,Description="Odds ratio">"#,
        r#"##INFO=<ID=LSEQ,Number=1,Type=String,Description="5' flanking seq">"#,
        r#"##INFO=<ID=RSEQ,Number=1,Type=String,Description="3' flanking seq">"#,
        r#"##INFO=<ID=STATUS,Number=1,Type=String,Description="Somatic or germline status">"#,
        r#"##INFO=<ID=P0.01Likely,Number=0,Type=Flag,Description="Likely candidate but p-value > 0.01/5**vd2 (means the evidence in tumor sample might be weak, e.g. small diff in AF)">"#,
        r#"##INFO=<ID=InDelLikely,Number=0,Type=Flag,Description="Likely indels more than 2bp are not considered somatic (weak evidence of presence in normal samples)">"#,
    ];
    for line in &info_lines {
        writeln!(w, "{}", line)?;
    }

    // ports var2vcf_paired.pl:76-93 (zero heap String allocations)
    let mut s_q = crate::types::StackBuf::<64>::new();
    let mut s_mapq = crate::types::StackBuf::<64>::new();
    let mut s_p = crate::types::StackBuf::<64>::new();
    let mut s_sn = crate::types::StackBuf::<64>::new();
    let mut s_d = crate::types::StackBuf::<64>::new();
    let mut s_v = crate::types::StackBuf::<64>::new();
    let mut s_msi = crate::types::StackBuf::<64>::new();
    let mut s_nm = crate::types::StackBuf::<64>::new();
    let mut s_c = crate::types::StackBuf::<64>::new();
    let mut s_af = crate::types::StackBuf::<64>::new();
    let mut s_pval = crate::types::StackBuf::<64>::new();

    let qmean = header_float_label(config.vcf_qmean.as_ref(), "22.5", PAIRED_QMEAN, &mut s_q);
    let mqmean = header_float_label(config.vcf_mapq.as_ref(), "0", PAIRED_MQMEAN, &mut s_mapq);
    let pmean = header_float_label(config.vcf_pmean.as_ref(), "8", PAIRED_PMEAN, &mut s_p);
    let sn = header_float_label(config.vcf_sn.as_ref(), "1.5", PAIRED_SN, &mut s_sn);
    let min_depth = header_int_label(
        config.vcf_min_depth.as_ref(),
        "5",
        PAIRED_MIN_DEPTH,
        &mut s_d,
    );
    let min_vd = header_int_label(config.vcf_min_vd.as_ref(), "3", PAIRED_VAR_DEPTH, &mut s_v);
    let opt_i = header_float_label(config.vcf_msi.as_ref(), "12", PAIRED_OPT_I, &mut s_msi);
    let opt_m = header_float_label(config.vcf_max_nm.as_ref(), "5.25", PAIRED_OPT_M, &mut s_nm);
    let opt_c = header_int_label(config.vcf_cluster_dist.as_ref(), "0", 0, &mut s_c);
    write_vcf_pval_label(&mut s_pval, config);
    let pval_lbl = s_pval.as_str();

    writeln!(
        w,
        r#"##FILTER=<ID=q{qmean},Description="Mean Base Quality Below {qmean}">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=Q{mqmean},Description="Mean Mapping Quality Below {mqmean}">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=p{pmean},Description="Mean Position in Reads Less than {pmean}">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=SN{sn},Description="Signal to Noise Less than {sn}">"#
    )?;
    writeln!(w, r#"##FILTER=<ID=Bias,Description="Strand Bias">"#)?;
    writeln!(
        w,
        r#"##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=MAF0.05,Description="Matched sample has AF > 0.05, thus not somatic">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=d{min_depth},Description="Total Depth < {min_depth}">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=v{min_vd},Description="Var Depth < {min_vd}">"#
    )?;
    if let Some(ref thr) = config.vcf_af {
        write_vcf_threshold_label(&mut s_af, thr, config.af_threshold);
        let af_lbl = s_af.as_str();
        writeln!(
            w,
            r#"##FILTER=<ID=f{af_lbl},Description="Allele frequency < {af_lbl}">"#
        )?;
    } else {
        let freq_val = config.af_threshold;
        writeln!(
            w,
            r#"##FILTER=<ID=f{freq_val},Description="Allele frequency < {freq_val}">"#
        )?;
    }
    writeln!(
        w,
        r#"##FILTER=<ID=P{pval_lbl},Description="Not significant with p-value > {pval_lbl}">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=DIFF0.2,Description="Non-somatic or LOH and allele frequency difference < 0.2">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=MSI{opt_i},Description="Variant in MSI region with {opt_i} non-monomer MSI or 12 monomer MSI">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=NM{opt_m},Description="Mean mismatches in reads >= {opt_m}, thus likely false positive">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=InGap,Description="The somatic variant is in the deletion gap, thus likely false positive">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=InIns,Description="The somatic variant is adjacent to an insertion variant">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=Cluster{opt_c}bp,Description="Two somatic variants are within {opt_c} bp">"#
    )?;
    writeln!(
        w,
        r#"##FILTER=<ID=LongAT,Description="The somatic variant is flanked by long A/T (>=14)">"#
    )?;

    // ports var2vcf_paired.pl:94-112
    let format_lines = [
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">"#,
        r#"##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">"#,
        r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"#,
        r#"##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">"#,
        r#"##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">"#,
        r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"#,
        r#"##FORMAT=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">"#,
        r#"##FORMAT=<ID=BIAS,Number=2,Type=String,Description="Strand Bias Info">"#,
        r#"##FORMAT=<ID=PMEAN,Number=1,Type=Float,Description="The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">"#,
        r#"##FORMAT=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">"#,
        r#"##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">"#,
        r#"##FORMAT=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">"#,
        r#"##FORMAT=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">"#,
        r#"##FORMAT=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Odds ratio">"#,
        r#"##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Mean Mapping Quality">"#,
        r#"##FORMAT=<ID=SN,Number=1,Type=Float,Description="Signal to noise">"#,
        r#"##FORMAT=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">"#,
        r#"##FORMAT=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">"#,
    ];
    for line in &format_lines {
        writeln!(w, "{}", line)?;
    }

    // ports var2vcf_paired.pl:43-49,115 — sample names
    let (sample, samplem) = match config.sample.split_once('|') {
        Some((s, m)) if !m.is_empty() => (s.to_string(), m.to_string()),
        Some((s, _)) => (s.to_string(), format!("{}-match", s)),
        None => (config.sample.clone(), format!("{}-match", config.sample)),
    };
    let sample_nowhitespace = sample.replace(char::is_whitespace, "_");
    let samplem_nowhitespace = samplem.replace(char::is_whitespace, "_");
    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\t{}",
        sample_nowhitespace, samplem_nowhitespace
    )?;

    Ok(())
}

/// One formatted paired record awaiting the delayed print.
///
/// Stores the entire VCF row (`CHROM\tPOS\t.\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL`)
/// in a single contiguous [`String`] with byte offsets `[filter_start..filter_end]`
/// pointing to the `FILTER` column.
#[derive(Debug, Clone)]
pub struct PairedVcfRecord {
    line: String,
    filter_start: usize,
    filter_end: usize,
    /// Raw start (for the Cluster rule :229 and pvs tracking :245).
    pub(crate) start: i32,
    pub(crate) type_snv: bool,
    /// Whether STATUS matches /Somatic/ (Cluster-rule precondition).
    pub(crate) status_somatic: bool,
}

impl PairedVcfRecord {
    /// Full tab-delimited paired VCF line (without trailing newline).
    #[inline]
    pub fn as_str(&self) -> &str {
        &self.line
    }

    /// `FILTER` column slice (`"PASS"` or `;`-joined filter tags).
    #[inline]
    pub fn filter_str(&self) -> &str {
        self.line
            .get(self.filter_start..self.filter_end)
            .unwrap_or("")
    }

    /// Returns `true` when the `FILTER` column equals `"PASS"`.
    #[inline]
    pub fn is_pass(&self) -> bool {
        self.filter_str() == "PASS"
    }

    #[inline]
    pub(crate) fn capacity(&self) -> usize {
        self.line.capacity()
    }

    #[inline]
    pub(crate) fn filter_offsets(&self) -> (usize, usize) {
        (self.filter_start, self.filter_end)
    }

    #[inline]
    pub(crate) fn set_cluster_filter(&mut self, cluster_tag: &str) {
        replace_vcf_filter_span(
            &mut self.line,
            self.filter_start,
            &mut self.filter_end,
            cluster_tag,
        );
    }

    #[inline]
    pub(crate) fn from_parts_checked(
        line: String,
        filter_start: usize,
        filter_end: usize,
        start: i32,
        type_snv: bool,
        status_somatic: bool,
    ) -> Option<Self> {
        if !is_valid_vcf_filter_span(&line, filter_start, filter_end, 10) {
            return None;
        }
        Some(Self {
            line,
            filter_start,
            filter_end,
            start,
            type_snv,
            status_somatic,
        })
    }
}

/// Perl numification of non-numeric strings ("NA" → 0, with warning).
fn perl_num(s: &str) -> f64 {
    s.parse().unwrap_or(0.0)
}

/// Odds-ratio transform for per-sample fields. Ports var2vcf_paired.pl:146-155:
/// "Inf" → 0; (0,1) → sprintf("%.2f", 1/x). Returns (string, numeric-after-format);
/// the numeric form feeds the Bias comparisons (:175/:186).
fn paired_or_transform(raw: &str) -> (crate::types::StackBuf<32>, f64) {
    let mut s = crate::types::StackBuf::<32>::new();
    if raw == "Inf" {
        let _ = s.write_str("0");
        return (s, 0.0);
    }
    let x = perl_num(raw);
    if x < 1.0 && x > 0.0 {
        let _ = write!(s, "{:.2}", crate::types::round_half_even(1.0 / x, 2));
        let v = s.as_str().parse::<f64>().unwrap_or(0.0);
        (s, v)
    } else {
        let _ = s.write_str(raw);
        (s, x)
    }
}

/// BIAS field transform. Ports var2vcf_paired.pl:222-225: ";" → ","; bare "0" → "0,0".
#[inline]
fn push_paired_bias(buf: &mut String, bias: &str) {
    if bias == "0" {
        buf.push_str("0,0");
    } else {
        for c in bias.chars() {
            buf.push(if c == ';' { ',' } else { c });
        }
    }
}

/// GT formula shared by both samples. Ports var2vcf_paired.pl:220-221.
fn paired_gt(af: f64, freq: f64, gt_freq: f64) -> &'static str {
    if 1.0 - af < gt_freq {
        "1/1"
    } else if af >= 0.5 {
        "1/0"
    } else if af >= freq {
        "0/1"
    } else {
        "0/0"
    }
}

/// Record QUAL column. Ports var2vcf_paired.pl:228 —
/// `int(log($vd1)/log(2) * $qual1)` on the sample chosen by `$vd1 > $vd2`
/// (ties take the normal side), where qual1/qual2 are the table cells printed
/// at HALF_EVEN 1dp — pass the already-rounded values. Perl lacks a zero guard;
/// both-zero cannot occur for emitted records (placeholders always keep one
/// real side), guarded to 0.
fn paired_qual_col(vd1: i32, qual1: f64, vd2: i32, qual2: f64) -> i32 {
    let (vd, q) = if vd1 > vd2 {
        (vd1, qual1)
    } else {
        (vd2, qual2)
    };
    if vd <= 0 {
        0
    } else {
        ((vd as f64).ln() / std::f64::consts::LN_2 * q) as i32
    }
}

/// Append one per-sample FORMAT column directly into `buf` without clearing it.
/// Ports var2vcf_paired.pl:240-241.
fn push_paired_sample_field(
    buf: &mut String,
    v: &crate::variations::Variant,
    gt: &str,
    sbf: f64,
    oddratio: &(impl std::ops::Deref<Target = str>, f64),
    fisher: bool,
) {
    let mq =
        crate::types::round_half_even(crate::types::round_half_even(v.mean_mapping_quality, 1), 0);
    buf.push_str(gt);
    buf.push(':');
    let _ = write!(buf, "{}", v.total_pos_cov);
    buf.push(':');
    let _ = write!(buf, "{}", v.pos_coverage);
    buf.push(':');
    let _ = write!(buf, "{},{}", v.fwd, v.rev);
    buf.push(':');
    let _ = write!(buf, "{},{}", v.ref_fwd_cov, v.ref_rev_cov);
    buf.push(':');
    let _ = write!(
        buf,
        "{}",
        i64::from(v.ref_fwd_cov) + i64::from(v.ref_rev_cov)
    );
    buf.push(',');
    let _ = write!(buf, "{}", v.pos_coverage);
    buf.push(':');
    push_vcf_rounded(buf, fisher, v.frequency, 4);
    buf.push(':');
    push_paired_bias(buf, &v.bias);
    buf.push(':');
    push_vcf_rounded(buf, fisher, v.mean_position, 1);
    buf.push(':');
    buf.push(if v.pstd { '1' } else { '0' });
    buf.push(':');
    push_vcf_rounded(buf, fisher, v.mean_quality, 1);
    buf.push(':');
    buf.push(if v.qstd { '1' } else { '0' });
    buf.push(':');
    push_vcf_rounded(buf, fisher, sbf, 5);
    buf.push(':');
    buf.push_str(&oddratio.0);
    buf.push(':');
    let _ = write!(buf, "{:.0}", mq);
    buf.push(':');
    push_vcf_rounded(buf, fisher, v.hq_lq_ratio, 3);
    buf.push(':');
    push_vcf_rounded(buf, fisher, v.hi_freq, 4);
    buf.push(':');
    push_vcf_rounded(buf, fisher, v.extra_freq, 4);
    buf.push(':');
    if v.nm > 0.0 {
        push_vcf_rounded(buf, fisher, v.nm, 1);
    } else {
        buf.push('0');
    }
}

/// Per-sample FORMAT string, written into `buf` (cleared first).
/// Ports var2vcf_paired.pl:240-241.
#[cfg(test)]
fn paired_sample_field_into(
    buf: &mut String,
    v: &crate::variations::Variant,
    gt: &str,
    sbf: f64,
    oddratio: &(impl std::ops::Deref<Target = str>, f64),
    fisher: bool,
) {
    buf.clear();
    // 192 B covers the 21-field sample column at typical widths; longer
    // alleles regrow once and the capacity is retained across rows.
    buf.reserve(192);
    push_paired_sample_field(buf, v, gt, sbf, oddratio, fisher);
}

static ZERO_VARIANT: std::sync::OnceLock<crate::variations::Variant> = std::sync::OnceLock::new();

/// Format one somatic variant into a paired VCF record.
/// Ports var2vcf_paired.pl:127-242 field-by-field. The delayed-print /
/// Cluster state machine lives in `finish_paired_records`.
pub fn format_paired_vcf_record(
    sv: &crate::postprocess::somatic::SomaticVariant,
    chr: &str,
    config: &Config,
) -> PairedVcfRecord {
    let v1 = &sv.variant1;
    let v2 = sv
        .variant2
        .as_ref()
        .unwrap_or_else(|| ZERO_VARIANT.get_or_init(crate::variations::Variant::new));

    // Fisher values exactly as somatic_61 computes them (matches the oracle's
    // testsomatic.R / --fisher columns; degenerate all-zero tables give p=1/OR=0
    // in both R and our FisherExact — verified empirically in vardict-oracle).
    let fisher1 = crate::stats::FisherExact::new_with_mode(
        v1.ref_fwd_cov,
        v1.ref_rev_cov,
        v1.fwd,
        v1.rev,
        config.fisher,
    );
    let p1 = fisher1.get_p_value();
    let or1 = paired_or_transform(&fisher1.odds_ratio_buf());

    let fisher2 = crate::stats::FisherExact::new_with_mode(
        v2.ref_fwd_cov,
        v2.ref_rev_cov,
        v2.fwd,
        v2.rev,
        config.fisher,
    );
    let p2 = fisher2.get_p_value();
    let or2 = paired_or_transform(&fisher2.odds_ratio_buf());

    let tref = (i64::from(v1.total_pos_cov) - i64::from(v1.pos_coverage))
        .clamp(0, i64::from(i32::MAX)) as i32;
    let rref = (i64::from(v2.total_pos_cov) - i64::from(v2.pos_coverage))
        .clamp(0, i64::from(i32::MAX)) as i32;
    let fisher_total = crate::stats::FisherExact::new_with_mode(
        v1.pos_coverage,
        tref,
        v2.pos_coverage,
        rref,
        config.fisher,
    );
    // min(greater, less), mirroring testsomatic.R's alternative selection (:31-43).
    let ptot = fisher_total
        .get_p_value_greater()
        .min(fisher_total.get_p_value_less());
    // Total oddratio: Inf→0 only, no reciprocal (:143-145).
    let otot_raw = fisher_total.odds_ratio_buf();
    let sor_str: &str = if otot_raw.as_str() == "Inf" {
        "0"
    } else {
        otot_raw.as_str()
    };

    // Filter comparisons run against table-rounded values: the oracle pipeline
    // formats through the .var table before Perl numifies, so compare at table
    // precision.
    let dp1 = v1.total_pos_cov;
    let vd1 = v1.pos_coverage;
    let dp2 = v2.total_pos_cov;
    let vd2 = v2.pos_coverage;
    let af1 = crate::types::round_half_even(v1.frequency, 4);
    let af2 = crate::types::round_half_even(v2.frequency, 4);
    let qual1 = crate::types::round_half_even(v1.mean_quality, 1);
    let qual2 = crate::types::round_half_even(v2.mean_quality, 1);

    let var_type: &str = if v1.var_type.is_empty() {
        "REF" // :139
    } else {
        &v1.var_type
    };
    let status = &sv.status;

    // endVariant supplies shift3/msi/msint/lseq/rseq (= table columns 48-52).
    let end_var = if sv.end_is_variant2 {
        sv.variant2.as_ref().unwrap_or(v1)
    } else {
        v1
    };

    let eff_min_depth = config.vcf_min_depth.as_ref().map_or(PAIRED_MIN_DEPTH, |t| {
        sanitize_nonneg_i32(t.value, PAIRED_MIN_DEPTH)
    });
    let eff_min_vd = config.vcf_min_vd.as_ref().map_or(PAIRED_VAR_DEPTH, |t| {
        sanitize_nonneg_i32(t.value, PAIRED_VAR_DEPTH)
    });
    let eff_freq = config.vcf_af.as_ref().map_or(config.af_threshold, |t| {
        sanitize_nonneg_f64(t.value, config.af_threshold)
    });
    let eff_pmean = config
        .vcf_pmean
        .as_ref()
        .map_or(PAIRED_PMEAN, |t| sanitize_nonneg_f64(t.value, PAIRED_PMEAN));
    let eff_qmean = config
        .vcf_qmean
        .as_ref()
        .map_or(PAIRED_QMEAN, |t| sanitize_nonneg_f64(t.value, PAIRED_QMEAN));
    let eff_mapq = config.vcf_mapq.as_ref().map_or(PAIRED_MQMEAN, |t| {
        sanitize_nonneg_f64(t.value, PAIRED_MQMEAN)
    });
    let eff_sn = config
        .vcf_sn
        .as_ref()
        .map_or(PAIRED_SN, |t| sanitize_nonneg_f64(t.value, PAIRED_SN));
    let eff_nm = config
        .vcf_max_nm
        .as_ref()
        .map_or(PAIRED_OPT_M, |t| sanitize_nonneg_f64(t.value, PAIRED_OPT_M));
    let eff_msi = config
        .vcf_msi
        .as_ref()
        .map_or(PAIRED_OPT_I, |t| sanitize_nonneg_f64(t.value, PAIRED_OPT_I));
    let eff_gt_freq = effective_gt_freq(config);

    // Normal depth & quality filters (:162-163, :177-183); Q gate dead at MQMEAN=0,
    // kept for fidelity. Without -M (`!config.vcf_paired_strict`), a fully clean NORMAL
    // (`normal_clean == true`) unconditionally overrides any non-PASS FILTER to "PASS"
    // (:217-219), so check `normal_clean` first and short-circuit all tumor/Bias/MSI
    // filter evaluation. With -M (`config.vcf_paired_strict`), germline rescue is disabled.
    let normal_clean = !config.vcf_paired_strict
        && !(dp2 < eff_min_depth
            || vd2 < eff_min_vd
            || af2 < eff_freq
            || crate::types::round_half_even(v2.mean_position, 1) < eff_pmean
            || (!v2.pstd && vd2 < eff_min_depth)
            || qual2 < eff_qmean
            || crate::types::round_half_even(v2.mean_mapping_quality, 1) < eff_mapq
            || crate::types::round_half_even(v2.hq_lq_ratio, 3) < eff_sn
            || crate::types::round_half_even(if v2.nm > 0.0 { v2.nm } else { 0.0 }, 1) >= eff_nm);

    let mut p_likely = false;
    let mut indel_likely = false;
    let mut tag_scratch = VcfFloatTagScratch::new();
    let mut filters = StackFilterList::new();
    if !normal_clean {
        let [s_f, s_p, s_q, s_mapq, s_sn, s_nm, s_msi, s_pval] = &mut tag_scratch.slots;
        let [s_d, s_v, _s_cluster] = &mut tag_scratch.int_slots;
        let pmean1 = crate::types::round_half_even(v1.mean_position, 1);
        let mapq1 = crate::types::round_half_even(v1.mean_mapping_quality, 1);
        let sn1 = crate::types::round_half_even(v1.hq_lq_ratio, 3);
        let nm1 = crate::types::round_half_even(if v1.nm > 0.0 { v1.nm } else { 0.0 }, 1);
        // Perl sees R's round-5 strings from testsomatic.R:31-47, numified before
        // the P:156/:175/:186/:189 threshold checks — compare at table precision.
        let ptot_r = crate::types::round_half_even(ptot, 5);
        let p1_r = crate::types::round_half_even(p1, 5);
        let p2_r = crate::types::round_half_even(p2, 5);

        let mut has_bias = false;
        let mut has_msi = false;

        // Tumor depth filters with StrongSomatic exception (:156-161).
        let ss_exception = status == "StrongSomatic" && ptot_r < 0.15 && af1 * vd1 as f64 >= 0.5;
        if dp1 < eff_min_depth && !ss_exception {
            let tag = match config.vcf_min_depth.as_ref() {
                Some(t) => format_int_threshold_filter_tag(s_d, "d", t, PAIRED_MIN_DEPTH, ""),
                None => FilterTag::Text("d5"),
            };
            filters.push(tag);
        }
        if vd1 < eff_min_vd && !ss_exception {
            let tag = match config.vcf_min_vd.as_ref() {
                Some(t) => format_int_threshold_filter_tag(s_v, "v", t, PAIRED_VAR_DEPTH, ""),
                None => FilterTag::Text("v3"),
            };
            filters.push(tag);
        }

        // Tumor quality filters (:165-175).
        if af1 < eff_freq {
            if let Some(ref af_thr) = config.vcf_af {
                filters.push(format_float_filter_tag(
                    s_f,
                    "f",
                    af_thr,
                    config.af_threshold,
                ));
            } else {
                filters.push(FilterTag::MinAf(config.af_threshold));
            }
        }
        if pmean1 < eff_pmean {
            if let Some(ref p_thr) = config.vcf_pmean {
                filters.push(format_float_filter_tag(s_p, "p", p_thr, PAIRED_PMEAN));
            } else {
                filters.push(FilterTag::Text("p8"));
            }
        }
        if !v1.pstd && vd1 < eff_min_depth {
            filters.push(FilterTag::Text("pSTD"));
        }
        if qual1 < eff_qmean {
            if let Some(ref q_thr) = config.vcf_qmean {
                filters.push(format_float_filter_tag(s_q, "q", q_thr, PAIRED_QMEAN));
            } else {
                filters.push(FilterTag::Text("q22.5"));
            }
        }
        let mapq_low = mapq1 < eff_mapq;
        let mapq_snv_low = mapq1 < 10.0 && var_type == "SNV";
        if mapq_low || mapq_snv_low {
            let q_tag = if let Some(ref mq_thr) = config.vcf_mapq {
                format_float_filter_tag(s_mapq, "Q", mq_thr, PAIRED_MQMEAN)
            } else {
                FilterTag::Text("Q0")
            };
            if mapq_low {
                filters.push(q_tag);
            }
            if mapq_snv_low {
                filters.push(q_tag);
            }
        }
        if sn1 < eff_sn {
            if let Some(ref sn_thr) = config.vcf_sn {
                filters.push(format_float_filter_tag(s_sn, "SN", sn_thr, PAIRED_SN));
            } else {
                filters.push(FilterTag::Text("SN1.5"));
            }
        }
        if nm1 >= eff_nm {
            if let Some(ref nm_thr) = config.vcf_max_nm {
                filters.push(format_float_filter_tag(s_nm, "NM", nm_thr, PAIRED_OPT_M));
            } else {
                filters.push(FilterTag::Text("NM5.25"));
            }
        }
        if v1.bias == "2;1"
            && p1_r < 0.01
            && (or1.1 > 5.0 || or1.1 == 0.0)
            && i64::from(v1.end_pos) - i64::from(v1.start_pos) < 100
        {
            filters.push(FilterTag::Text("Bias"));
            has_bias = true;
        }

        // Normal-side Bias appended to the MAIN list when not already present (:185-186).
        if !has_bias
            && v2.bias == "2;1"
            && p2_r < 0.01
            && (or2.1 > 5.0 || or2.1 == 0.0)
            && i64::from(v1.end_pos) - i64::from(v1.start_pos) < 100
        {
            filters.push(FilterTag::Text("Bias"));
        }

        // MSI gates (:188-194).
        let msilen = end_var.msint as f64;
        let msi = crate::types::round_half_even(end_var.msi, 3);
        let mut push_msi = false;
        if (msi > eff_msi && msilen > 1.0) || (msi > 12.0 && msilen == 1.0) {
            if !(status == "StrongSomatic" && ptot_r < 0.0005) {
                push_msi = true;
                has_msi = true;
            }
        }
        let len_diff = v1.ref_allele.len().abs_diff(v1.var_allele.len());
        if usize::try_from(end_var.msint).ok() == Some(len_diff) && !has_msi {
            if (msi > eff_msi && msilen > 1.0 && af1 < 0.35 && af2 < 0.35)
                || (msi > 12.0 && msilen == 1.0 && af1 < 0.35 && af2 < 0.35)
            {
                push_msi = true;
            }
        }
        if push_msi {
            let msi_tag = if let Some(ref msi_thr) = config.vcf_msi {
                format_float_filter_tag(s_msi, "MSI", msi_thr, PAIRED_OPT_I)
            } else {
                FilterTag::Text("MSI12")
            };
            filters.push(msi_tag);
        }

        // Paired strict `-M` (`--paired-strict`, var2vcf_paired.pl:198-206).
        if config.vcf_paired_strict {
            let eff_pval = effective_paired_pval_num(config);
            if ptot_r > eff_pval {
                let ss_pval_exc = status == "StrongSomatic"
                    && ((ptot_r < 0.25 && af1 > 0.1)
                        || (ptot_r < 0.5 && af1 > 0.20)
                        || (ptot_r < 0.15 && af1 > 0.05));
                if !ss_pval_exc {
                    *s_pval = crate::types::StackBuf::<80>::new();
                    let _ = s_pval.write_str("P");
                    write_vcf_pval_label(s_pval, config);
                    filters.push(FilterTag::Text(s_pval.as_str()));
                }
            } else if status.contains("LikelySomatic") && ptot_r > 0.05 / 5.0_f64.powi(vd2.max(0)) {
                p_likely = true;
            } else if status.contains("Likely")
                && var_type != "SNV"
                && !(v1.ref_allele.len() <= 2 && v1.var_allele.len() <= 2)
            {
                indel_likely = true;
            }
        }
    }

    let gt = paired_gt(af1, eff_freq, eff_gt_freq);
    let gtm = paired_gt(af2, eff_freq, eff_gt_freq);

    let qual_col = paired_qual_col(vd1, qual1, vd2, qual2);

    let sample_tumor = config.sample.split('|').next().unwrap_or(&config.sample);

    // Single contiguous allocation for the entire paired VCF record (`896 B` covers
    // common rows without regrowth).
    let mut line = String::with_capacity(896);
    line.push_str(chr);
    line.push('\t');
    let _ = write!(line, "{}", v1.start_pos);
    line.push_str("\t.\t");
    line.push_str(&v1.ref_allele);
    line.push('\t');
    line.push_str(&v1.var_allele);
    line.push('\t');
    let _ = write!(line, "{}\t", qual_col);

    let filter_start = line.len();
    filters.write_vcf_filter_field(&mut line);
    let filter_end = line.len();
    line.push('\t');

    // INFO (:237-239). No END tag and no SV tags in the paired writer.
    line.push_str("STATUS=");
    line.push_str(status);
    line.push_str(";SAMPLE=");
    for c in sample_tumor.chars() {
        line.push(if c.is_whitespace() { '_' } else { c });
    }
    line.push_str(";TYPE=");
    line.push_str(var_type);
    let _ = write!(line, ";DP={};VD={}", dp1, vd1);
    line.push_str(";AF=");
    push_vcf_rounded(&mut line, config.fisher, v1.frequency, 4);
    let _ = write!(line, ";SHIFT3={}", end_var.shift3);
    line.push_str(";MSI=");
    push_keep(&mut line, end_var.msi, 3);
    let _ = write!(line, ";MSILEN={}", end_var.msint);
    line.push_str(";SSF=");
    push_vcf_rounded(&mut line, config.fisher, ptot, 5);
    line.push_str(";SOR=");
    line.push_str(sor_str);
    line.push_str(";LSEQ=");
    if end_var.left_seq.is_empty() {
        line.push('0');
    } else {
        line.push_str(&end_var.left_seq);
    }
    line.push_str(";RSEQ=");
    if end_var.right_seq.is_empty() {
        line.push('0');
    } else {
        line.push_str(&end_var.right_seq);
    }
    if p_likely {
        line.push_str(";P0.01Likely");
    }
    if indel_likely {
        line.push_str(";InDelLikely");
    }

    line.push_str(
        "\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t",
    );
    push_paired_sample_field(&mut line, v1, gt, p1, &or1, config.fisher);
    line.push('\t');
    push_paired_sample_field(&mut line, v2, gtm, p2, &or2, config.fisher);

    PairedVcfRecord {
        line,
        filter_start,
        filter_end,
        start: v1.start_pos,
        type_snv: var_type == "SNV",
        status_somatic: status.contains("Somatic"),
    }
}

/// Streaming delayed-print state machine for paired VCF records (`var2vcf_paired.pl:126,229-250`).
/// Retains at most 1 lookbehind record (`pending`) per chromosome so whole-chromosome
/// emission uses `O(1)` memory.
#[derive(Debug)]
pub struct PairedRecordStreamer {
    pending: Option<PairedVcfRecord>,
    pvs: i32,
    opt_c: i32,
    pass_only: bool,
    cluster_tag: crate::types::StackBuf<80>,
}

impl PairedRecordStreamer {
    pub fn new(opt_c: i32) -> Self {
        Self::with_pass_only(opt_c, false)
    }

    pub fn with_pass_only(opt_c: i32, pass_only: bool) -> Self {
        Self {
            pending: None,
            pvs: 0,
            opt_c,
            pass_only,
            cluster_tag: format_cluster_tag_from_int(opt_c),
        }
    }

    pub fn from_config(config: &Config) -> Self {
        Self {
            pending: None,
            pvs: 0,
            opt_c: effective_vcf_cluster_dist(config),
            pass_only: config.vcf_pass_only,
            cluster_tag: format_cluster_tag_from_config(config),
        }
    }

    pub fn push<E>(
        &mut self,
        mut cur: PairedVcfRecord,
        mut emit: impl FnMut(&str) -> Result<(), E>,
    ) -> Result<(), E> {
        if let Some(mut prev) = self.pending.take() {
            if prev.is_pass()
                && prev.status_somatic
                && prev.type_snv
                && cur.is_pass()
                && cur.status_somatic
                && cur.type_snv
                && i64::from(cur.start) - i64::from(self.pvs) < i64::from(self.opt_c)
            {
                prev.set_cluster_filter(self.cluster_tag.as_str());
                cur.set_cluster_filter(self.cluster_tag.as_str());
            }
            if !self.pass_only || prev.is_pass() {
                emit(prev.as_str())?;
            }
        }
        if cur.type_snv && cur.is_pass() {
            self.pvs = cur.start; // :245
        }
        self.pending = Some(cur);
        Ok(())
    }

    pub fn finish<E>(mut self, mut emit: impl FnMut(&str) -> Result<(), E>) -> Result<(), E> {
        if let Some(prev) = self.pending.take() {
            if !self.pass_only || prev.is_pass() {
                emit(prev.as_str())?;
            }
        }
        Ok(())
    }
}

/// Delayed-print state machine with optional `-S` (`--pass-only`) post-Cluster suppression
/// (`var2vcf_paired.pl:126,229-236,243-245,248-250`).
pub fn finish_paired_records_with_pass_only(
    records: Vec<PairedVcfRecord>,
    opt_c: i32,
    pass_only: bool,
) -> Vec<String> {
    let mut out = Vec::with_capacity(records.len());
    let mut streamer = PairedRecordStreamer::with_pass_only(opt_c, pass_only);
    for cur in records {
        let _ = streamer.push(cur, |line| -> Result<(), std::convert::Infallible> {
            out.push(line.to_string());
            Ok(())
        });
    }
    let _ = streamer.finish(|line| -> Result<(), std::convert::Infallible> {
        out.push(line.to_string());
        Ok(())
    });
    out
}

/// Delayed-print state machine without `-S` (`--pass-only`) suppression.
pub fn finish_paired_records(records: Vec<PairedVcfRecord>, opt_c: i32) -> Vec<String> {
    finish_paired_records_with_pass_only(records, opt_c, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chrom_compare() {
        assert_eq!(chrom_compare("chr1", "chr2"), std::cmp::Ordering::Less);
        assert_eq!(chrom_compare("chr10", "chr2"), std::cmp::Ordering::Greater);
        assert_eq!(chrom_compare("chrX", "chrY"), std::cmp::Ordering::Less);
        assert_eq!(chrom_compare("chrY", "chrM"), std::cmp::Ordering::Less);
    }

    #[test]
    fn test_perl_row_skipped_empty_ref() {
        // ports var2vcf_valid.pl:141 — `next unless ($ref)` drops the -p null
        // records (empty REF cell) and Perl-falsy "0"; real alleles survive.
        assert!(perl_row_skipped(""));
        assert!(perl_row_skipped("0"));
        assert!(!perl_row_skipped("G"));
        assert!(!perl_row_skipped("ACGT"));
    }

    #[test]
    fn test_header() {
        let config = Config::default();
        let mut lengths = HashMap::new();
        lengths.insert("chr1".to_string(), 1000000);
        lengths.insert("chr2".to_string(), 2000000);
        let mut buf = Vec::new();
        print_vcf_header(&mut buf, &config, &lengths).unwrap();
        let header = String::from_utf8(buf).unwrap();
        assert!(header.contains("##fileformat=VCFv4.2"));
        assert!(header.contains("##source=VarDict_v1.8.2"));
        assert!(header.contains("##contig=<ID=chr1,length=1000000>"));
        assert!(header.contains("##contig=<ID=chr2,length=2000000>"));
        assert!(header.contains("##INFO=<ID=SAMPLE"));
        assert!(header.contains("##INFO=<ID=PMEAN,Number=1,Type=Float,Description=\"The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call\">"));
        assert!(header.contains("##FILTER=<ID=q22.5"));
        assert!(header.contains("##FORMAT=<ID=GT"));
    }

    #[test]
    fn test_contig_lengths_come_from_fai_not_bam_header() {
        // issue #320: when the reference .fai exists, ##contig lengths come from
        // it — never from the BAM header or hardcoded defaults.
        let dir =
            std::env::temp_dir().join(format!("vardict-fai-contig-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("ref.fa");
        std::fs::write(&fa, ">chrFAI\nA\n").unwrap();
        // Distinctive .fai length vs a conflicting BAM-header length.
        // (Trailing cols are not realistic for this single-line fasta; the
        // loader reads only the name/length columns.)
        std::fs::write(fa.with_extension("fa.fai"), "chrFAI\t987654321\t6\t1\t2\n").unwrap();
        let mut config = Config::default();
        config.reference = fa.to_string_lossy().into_owned();
        let mut bam_lengths = HashMap::new();
        bam_lengths.insert("chrFAI".to_string(), 111);
        let mut buf = Vec::new();
        print_vcf_header(&mut buf, &config, &bam_lengths).unwrap();
        let header = String::from_utf8(buf).unwrap();
        assert!(header.contains("##contig=<ID=chrFAI,length=987654321>"));
        assert!(!header.contains("length=111"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_record() {
        let mut variant = Variant::new();
        variant.ref_allele = "A".to_string();
        variant.var_allele = "G".to_string();
        variant.start_pos = 100;
        variant.end_pos = 100;
        variant.total_pos_cov = 100;
        variant.pos_coverage = 10;
        variant.ref_fwd_cov = 40;
        variant.ref_rev_cov = 50;
        variant.fwd = 4;
        variant.rev = 6;
        variant.frequency = 0.1;
        variant.mean_position = 20.0;
        variant.pstd = true;
        variant.mean_quality = 30.0;
        variant.qstd = true;
        variant.mean_mapping_quality = 60.0;
        variant.hq_lq_ratio = 2.0;
        variant.hi_freq = 0.1;
        variant.msi = 0.5;
        variant.msint = 1;
        variant.nm = 1.0;
        variant.hicnt = 10;
        variant.hicov = 100;
        variant.left_seq = "ACGT".to_string();
        variant.right_seq = "TGCA".to_string();
        variant.var_type = "SNV".to_string();
        variant.bias = "2;2".to_string();

        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let config = Config::default();
        let filters = filter_variant(&variant, &config, 0, None, 1.0, "1.5");
        let rec = format_vcf_record(&variant, &region, &config, &filters, 1.0, "1.5", None);
        assert!(rec.contains("chr1\t100\t.\tA\tG"));
        assert!(rec.contains("SAMPLE=SAMPLE"));
    }

    #[test]
    fn format_vcf_record_into_matches_one_shot_and_clears() {
        // Task 6 pin (behavior only): the reused-buffer writer is byte-equal
        // to the one-shot fn, `clear()`s per row, and keeps the -p survivor
        // path (ref==alt → ALT ".", GT 0/0) plus the amp GDAMP path.
        let config = Config::default();
        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let mut survivor = amplicon_passing_variant();
        survivor.var_type = String::new();
        survivor.var_allele = survivor.ref_allele.clone();
        let filters = filter_variant(&survivor, &config, 0, None, 1.0, "1.5");
        let rec = format_vcf_record(&survivor, &region, &config, &filters, 1.0, "1.5", None);
        let mut buf = String::new();
        format_vcf_record_into(
            &mut buf, &survivor, &region, &config, &filters, 1.0, "1.5", None,
        );
        assert_eq!(buf, rec);
        assert!(buf.contains("\t0/0:"), "GT 0/0: {buf}");

        // Reuse the same buffer for a longer SNV record: exact match means
        // no tail leakage from the previous row and no elided row.
        let mut snv = amplicon_passing_variant();
        snv.left_seq = "ACGTACGTACGTACGTACGT".to_string();
        snv.right_seq = "TGCATGCATGCATGCATGCA".to_string();
        let filters2 = filter_variant(&snv, &config, 0, None, 1.0, "1.5");
        let rec2 = format_vcf_record(&snv, &region, &config, &filters2, 1.0, "1.5", None);
        format_vcf_record_into(
            &mut buf, &snv, &region, &config, &filters2, 1.0, "1.5", None,
        );
        assert_eq!(buf, rec2);

        // Amplicon path: GDAMP tags, no DUPRATE/SPLITREAD.
        let amp = Some((4, 4, 0, false));
        let filters3 = filter_variant(&snv, &config, 0, amp, 1.0, "1.5");
        let rec3 = format_vcf_record(&snv, &region, &config, &filters3, 1.0, "1.5", amp);
        format_vcf_record_into(&mut buf, &snv, &region, &config, &filters3, 1.0, "1.5", amp);
        assert_eq!(buf, rec3);
        assert!(buf.contains("GDAMP=4;TLAMP=4;NCAMP=0;AMPFLAG=0"));
    }

    #[test]
    fn test_filters() {
        let config = Config::default();

        // 1. Total depth filter (d3)
        let mut variant = Variant::new();
        variant.ref_allele = "A".to_string();
        variant.var_allele = "G".to_string();
        variant.total_pos_cov = 2; // dp < 3
        variant.hicnt = 1;
        variant.hi_freq = 0.1; // hicnt * hiaf = 0.1 < 0.5
        let filters = filter_variant(&variant, &config, 0, None, 1.0, "1.5");
        assert!(filters.contains(&"d3".to_string()));

        // 2. Frequency filter tag uses config.af_threshold (= var2vcf_valid.pl's
        //    $Freq from -f, default 0.02; harness runs pass it explicitly).
        let mut variant = Variant::new();
        variant.ref_allele = "A".to_string();
        variant.var_allele = "G".to_string();
        variant.total_pos_cov = 100;
        variant.hicnt = 10;
        variant.hi_freq = 0.1;
        variant.frequency = 0.005; // af < 0.01 (default af_threshold)
        let filters = filter_variant(&variant, &config, 0, None, 1.0, "1.5");
        assert!(filters.contains(&"f0.01".to_string()));

        // With -f 0.001 the same AF stays above threshold: annotated never,
        // matching perl :159 (`push @filters,"f$Freq" if $af < $Freq`).
        let mut loose = config.clone();
        loose.af_threshold = 0.001;
        let filters = filter_variant(&variant, &loose, 0, None, 1.0, "1.5");
        assert!(!filters.iter().any(|f| f.starts_with('f')));

        // And a sub-threshold AF gets the exact -f value in the tag.
        variant.frequency = 0.0005;
        let filters = filter_variant(&variant, &loose, 0, None, 1.0, "1.5");
        assert!(filters.contains(&"f0.001".to_string()));
    }

    // ---- amplicon writer deltas (var2vcf_valid.pl $isamp branches) — W4.5d ----

    fn amplicon_passing_variant() -> Variant {
        // Passes every filter: dp>=3, hicnt>=2, af high, pmean>=8, pstd set,
        // qual>=22.5, mapq>=10, sn>=1.5, nm<5.25.
        let mut v = Variant::new();
        v.ref_allele = "A".to_string();
        v.var_allele = "G".to_string();
        v.start_pos = 100;
        v.end_pos = 100;
        v.total_pos_cov = 50;
        v.pos_coverage = 30;
        v.hicnt = 30;
        v.hicov = 50;
        v.ref_fwd_cov = 5;
        v.ref_rev_cov = 15;
        v.fwd = 10;
        v.rev = 20;
        v.frequency = 0.6;
        v.pstd = true;
        v.mean_position = 20.0;
        v.mean_quality = 35.0;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 2.0;
        v.nm = 1.0;
        v.var_type = "SNV".to_string();
        v
    }

    #[test]
    fn test_amp_record_has_gdamp_tags_and_no_sv_or_duprate() {
        let config = Config::default();
        let variant = amplicon_passing_variant();
        let region = Region {
            chr: "11".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let amp = Some((4, 4, 0, false));
        let filters = filter_variant(&variant, &config, 0, amp, 1.0, "1.5");
        assert_eq!(filters, Vec::<String>::new());
        let rec = format_vcf_record(&variant, &region, &config, &filters, 1.0, "1.5", amp);
        assert!(rec.contains("GDAMP=4;TLAMP=4;NCAMP=0;AMPFLAG=0"));
        assert!(!rec.contains("DUPRATE"));
        assert!(!rec.contains("SPLITREAD"));
        assert!(!rec.contains("SPANPAIR"));
    }

    #[test]
    fn test_amp_bias_filter_rules() {
        // ports var2vcf_valid.pl:185 — AMPBIAS iff gamp < tamp-ncamp || ampflag.
        let config = Config::default();
        let variant = amplicon_passing_variant();

        // gamp(2) < tamp(4) - ncamp(1) → AMPBIAS.
        let filters = filter_variant(&variant, &config, 0, Some((2, 4, 1, false)), 1.0, "1.5");
        assert!(filters.contains(&"AMPBIAS".to_string()));

        // gamp(3) == tamp(4) - ncamp(1), no flag → clean.
        let filters = filter_variant(&variant, &config, 0, Some((3, 4, 1, false)), 1.0, "1.5");
        assert!(!filters.contains(&"AMPBIAS".to_string()));

        // Counts equal but ampflag set → AMPBIAS (perl truthy $ampflag).
        let filters = filter_variant(&variant, &config, 0, Some((4, 4, 0, true)), 1.0, "1.5");
        assert!(filters.contains(&"AMPBIAS".to_string()));
    }

    #[test]
    fn test_amp_pstd_suppressed() {
        // ports var2vcf_valid.pl:161 — pSTD filter skipped when $isamp.
        let config = Config::default();
        let mut variant = amplicon_passing_variant();
        variant.pstd = false;
        variant.frequency = 0.1;

        let filters = filter_variant(&variant, &config, 0, Some((1, 1, 0, false)), 1.0, "1.5");
        assert!(!filters.contains(&"pSTD".to_string()));

        let filters = filter_variant(&variant, &config, 0, None, 1.0, "1.5");
        assert!(filters.contains(&"pSTD".to_string()));
    }

    #[test]
    fn test_amp_ref_row_alt_dot_gt_00_single_ad() {
        // ports var2vcf_valid.pl:216-231 — ref==alt rows (-p pileup reference
        // calls) become ALT=".", GT=0/0 with single-value AD in ANY mode; the
        // MRE11a slice emits these.
        let config = Config::default();
        let mut variant = amplicon_passing_variant();
        variant.var_type = String::new(); // → TYPE=REF (:144)
        variant.var_allele = variant.ref_allele.clone(); // ref==alt row (-p)
        let region = Region {
            chr: "11".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let amp = Some((1, 1, 0, false));
        let filters = filter_variant(&variant, &config, 0, amp, 1.0, "1.5");
        let rec = format_vcf_record(&variant, &region, &config, &filters, 1.0, "1.5", amp);
        let fields: Vec<&str> = rec.split('\t').collect();
        assert_eq!(fields[4], "."); // ALT
        let info = fields[7];
        assert!(info.starts_with("SAMPLE=SAMPLE;TYPE=REF;"));
        let sample_col = fields[9];
        let vals: Vec<&str> = sample_col.split(':').collect();
        assert_eq!(vals[0], "0/0"); // GT
                                    // AD = "rd,vd" whenever $vd ne "0" (perl :229) — REF rows still carry
                                    // their variant-depth column; single-value AD only when vd == 0.
        assert_eq!(vals[3], "20,30");

        let mut zero_vd = amplicon_passing_variant();
        zero_vd.ref_allele = "C".to_string();
        zero_vd.var_allele = "C".to_string();
        zero_vd.var_type = String::new();
        zero_vd.pos_coverage = 0;
        let rec0 = format_vcf_record(&zero_vd, &region, &config, &filters, 1.0, "1.5", amp);
        let vals0: Vec<&str> = rec0.split('\t').nth(9).unwrap().split(':').collect();
        assert_eq!(vals0[3], "20"); // rd only (:228)
    }

    #[test]
    fn test_simple_mode_still_emits_sv_and_duprate() {
        // Guard: None keeps the W4.5a simple behavior byte-for-byte.
        let config = Config::default();
        let mut variant = amplicon_passing_variant();
        variant.var_allele = "<DUP130>".to_string();
        variant.end_pos = 140; // SVLEN = end - start (perl :206)
        variant.sv = "3-4".to_string();
        variant.duprate = 0.25;
        let region = Region {
            chr: "chr20".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let filters = filter_variant(&variant, &config, 0, None, 1.0, "1.5");
        let rec = format_vcf_record(&variant, &region, &config, &filters, 1.0, "1.5", None);
        assert!(rec.contains("SVTYPE=DUP130;SVLEN=40"));
        assert!(rec.contains("SPLITREAD=3;SPANPAIR=4"));
        // P2: default (no-fisher) DUPRATE is the 1dp 36-col cell — HALF_EVEN
        // (0.25 → 0.2); the 2dp form is the --fisher path (SimpleOutputVariant
        // .java:223 vs :175). See test_duprate_pipeline_formats.
        assert!(rec.contains("DUPRATE=0.2"));
        assert!(!rec.contains("GDAMP"));
    }

    // ---- var2vcf_paired.pl port tests (W4.5c) ----

    fn somatic_sv(
        v1: crate::variations::Variant,
        v2: Option<crate::variations::Variant>,
        status: &str,
    ) -> crate::postprocess::somatic::SomaticVariant {
        crate::postprocess::somatic::SomaticVariant {
            variant1: v1,
            variant2: v2,
            status: status.to_string(),
            sv1: String::new(),
            sv2: String::new(),
            end_is_variant2: false,
        }
    }

    fn clean_normal() -> crate::variations::Variant {
        // Passes every normal-side gate: dp>=5, vd>=3, af>=0.01, pmean>=8,
        // pstd set, qual>=22.5, sn>=1.5, nm<5.25, bias != 2;1.
        let mut v = crate::variations::Variant::new();
        v.ref_allele = "A".to_string();
        v.var_allele = "G".to_string();
        v.total_pos_cov = 50;
        v.pos_coverage = 10;
        v.ref_fwd_cov = 20;
        v.ref_rev_cov = 20;
        v.fwd = 5;
        v.rev = 5;
        v.frequency = 0.2;
        v.mean_position = 30.0;
        v.pstd = true;
        v.mean_quality = 35.0;
        v.qstd = true;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 2.0;
        v.hi_freq = 0.2;
        v.bias = "0;0".to_string();
        v
    }

    #[test]
    fn test_paired_header_two_samples() {
        let mut config = Config::default();
        config.sample = "SAMPLE".to_string();
        config.af_threshold = 0.02;
        let lengths = HashMap::new();
        let mut buf = Vec::new();
        print_paired_vcf_header(&mut buf, &config, &lengths).unwrap();
        let h = String::from_utf8(buf).unwrap();
        // two-sample #CHROM with -match suffix (:43-49,115)
        assert!(h.contains(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tSAMPLE-match"
        ));
        // paired-specific INFO/FILTER lines (:69-75,:82-87)
        assert!(h.contains("##INFO=<ID=STATUS"));
        assert!(h.contains("##INFO=<ID=SSF"));
        assert!(h.contains(r#"##FILTER=<ID=MAF0.05"#));
        assert!(h.contains(r#"##FILTER=<ID=P0.05"#));
        assert!(h.contains(r#"##FILTER=<ID=DIFF0.2"#));
        assert!(h.contains(r#"##FILTER=<ID=d5,"#));
        assert!(h.contains(r#"##FILTER=<ID=v3,"#));
        // MQ header uses default MQMEAN=0 (:24)
        assert!(h.contains(r#"##FILTER=<ID=Q0,"#));
        // FORMAT block carries the 19 paired keys incl ADJAF/BIAS (:94-112)
        assert!(h.contains(r#"##FORMAT=<ID=ADJAF"#));
        assert!(h.contains(r#"##FORMAT=<ID=ODDRATIO"#));
    }

    #[test]
    fn test_paired_germline_passthrough_and_status_annotation() {
        let config = Config::default(); // freq 0.01
                                        // Tumor fails q22.5 only; normal is clean → germline passthrough (:217-219).
        let mut tumor = clean_normal();
        tumor.var_type = "SNV".to_string();
        tumor.start_pos = 100;
        tumor.end_pos = 100;
        tumor.mean_quality = 10.0;
        let rec = format_paired_vcf_record(
            &somatic_sv(tumor.clone(), Some(clean_normal()), "Germline"),
            "chr1",
            &config,
        );
        assert_eq!(rec.filter_str(), "PASS");
        let full = finish_paired_records(vec![rec], 0);
        let line = &full[0];
        assert!(line.contains("\tPASS\tSTATUS=Germline;SAMPLE=SAMPLE;TYPE=SNV;DP=50;VD=10;AF="));
        // BIAS comma-form (:222-225) and GT from af=0.2 → 0/1 (:220)
        assert!(line.contains(":40,10:0.2:0,0:"), "BIAS field: {}", line);
        assert!(line.contains("0/1:50:10:"), "tumor sample field: {}", line);

        // Same tumor but normal ALSO filtered (low depth) → filters stand.
        let mut bad_normal = clean_normal();
        bad_normal.total_pos_cov = 3; // < d5 → filters2 non-empty
        let rec2 = format_paired_vcf_record(
            &somatic_sv(tumor, Some(bad_normal), "Germline"),
            "chr1",
            &config,
        );
        assert_eq!(rec2.filter_str(), "q22.5");
    }

    #[test]
    fn test_paired_strongsomatic_depth_exception() {
        let config = Config::default();
        // Tumor dp1=4 < d5 and low qual; deep normal without variant makes the
        // total one-sided p tiny → StrongSomatic suppresses d5/v3 (:156-161),
        // leaving only q22.5. Normal vd2=0<3 keeps filters2 non-empty so the
        // germline passthrough cannot mask the result.
        let mut tumor = crate::variations::Variant::new();
        tumor.ref_allele = "C".to_string();
        tumor.var_allele = "T".to_string();
        tumor.total_pos_cov = 4;
        tumor.pos_coverage = 4;
        tumor.ref_fwd_cov = 0;
        tumor.ref_rev_cov = 0;
        tumor.fwd = 4;
        tumor.rev = 0;
        tumor.frequency = 1.0; // af*vd = 4 >= 0.5
        tumor.mean_position = 30.0;
        tumor.pstd = true;
        tumor.mean_quality = 10.0; // → q22.5
        tumor.qstd = true;
        tumor.mean_mapping_quality = 60.0;
        tumor.hq_lq_ratio = 2.0;
        tumor.var_type = "SNV".to_string();
        let mut normal = crate::variations::Variant::new();
        normal.total_pos_cov = 200;
        normal.pos_coverage = 0;
        normal.pstd = true;

        let rec = format_paired_vcf_record(
            &somatic_sv(tumor.clone(), Some(normal.clone()), "StrongSomatic"),
            "chr1",
            &config,
        );
        assert_eq!(rec.filter_str(), "q22.5", "d5/v3 suppressed, q22.5 stands");

        // Control: any other status keeps the depth filters (order :156→:169).
        let rec2 = format_paired_vcf_record(
            &somatic_sv(tumor, Some(normal), "Germline"),
            "chr1",
            &config,
        );
        assert_eq!(rec2.filter_str(), "d5;q22.5");
    }

    #[test]
    fn test_paired_ss_exception_compares_round5_ptot() {
        // Boundary: true total-Fisher p = 0.1499968 (<0.15 raw) but testsomatic.R
        // writes round(p,5) = 0.15 into the table and perl numifies "0.15", so
        // 0.15 < 0.15 is FALSE and the StrongSomatic depth exception must be
        // DENIED (d5 stands). Total table (refFwd=3, refRev=0, altFwd=210,
        // altRev=187); flips if the comparison ever sees unrounded p-values.
        let config = Config::default();
        let mut tumor = crate::variations::Variant::new();
        tumor.ref_allele = "C".to_string();
        tumor.var_allele = "T".to_string();
        tumor.total_pos_cov = 3; // dp1 < 5
        tumor.pos_coverage = 3; // vd1 = 3, tref = 0
        tumor.frequency = 1.0; // af*vd = 3 >= 0.5
        tumor.mean_position = 30.0;
        tumor.pstd = true;
        tumor.mean_quality = 10.0; // q22.5 remains visible when d5/v3 suppressed
        tumor.qstd = true;
        tumor.mean_mapping_quality = 60.0;
        tumor.hq_lq_ratio = 2.0;
        tumor.var_type = "SNV".to_string(); // SNV+mapq<10 pushes Q0
        let mut normal = clean_normal();
        normal.total_pos_cov = 397; // vd2=210 + rref=187
        normal.pos_coverage = 210;
        normal.mean_quality = 10.0; // keeps filters2 non-empty (no passthrough)

        let rec = format_paired_vcf_record(
            &somatic_sv(tumor.clone(), Some(normal.clone()), "StrongSomatic"),
            "chr1",
            &config,
        );
        assert!(
            rec.filter_str().contains("d5"),
            "round-5 ptot=0.15 must NOT grant the exception: {}",
            rec.filter_str()
        );

        // Control: clearly significant p (deep clean normal) still grants it.
        let mut deep = clean_normal();
        deep.total_pos_cov = 300;
        deep.pos_coverage = 0;
        deep.mean_quality = 10.0; // keeps filters2 non-empty (no passthrough)
        let rec2 = format_paired_vcf_record(
            &somatic_sv(tumor, Some(deep), "StrongSomatic"),
            "chr1",
            &config,
        );
        assert_eq!(
            rec2.filter_str(),
            "q22.5",
            "tiny p grants d5/v3 suppression"
        );
    }

    #[test]
    fn test_paired_gt_thresholds_and_bias_zero_form() {
        let config = Config::default();
        let mut tumor = clean_normal();
        tumor.frequency = 0.9; // 1-af < GTFREQ → 1/1
        tumor.bias = "0".to_string(); // bare zero → "0,0" (:224-225)
        let rec = format_paired_vcf_record(
            &somatic_sv(tumor, Some(clean_normal()), "Somatic"),
            "chr1",
            &config,
        );
        assert!(
            rec.as_str().contains("1/1:"),
            "gt 1/1 expected: {}",
            rec.as_str()
        );
        assert!(
            rec.as_str().contains(":0,0:"),
            "bias 0→0,0: {}",
            rec.as_str()
        );
    }

    #[test]
    fn test_finish_paired_records_delays_and_flushes_in_order() {
        let config = Config::default();
        let mk = |start: i32| {
            let mut t = clean_normal();
            t.start_pos = start;
            t.end_pos = start;
            t.var_type = "SNV".to_string();
            format_paired_vcf_record(&somatic_sv(t, None, "Somatic"), "chr1", &config)
        };
        // Two records: first must be buffered until the second arrives, final
        // flush emits both in order (:233-236,248-250). With c=0 no Cluster tag.
        let lines = finish_paired_records(vec![mk(200), mk(300)], 0);
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("chr1\t200\t."));
        assert!(lines[1].starts_with("chr1\t300\t."));
        assert!(!lines.iter().any(|l| l.contains("Cluster")));
    }

    // ---- P2 byte-exactness pass (writer formatting family) ----

    /// Shared record with the four divergence inputs set.
    fn p2_variant() -> Variant {
        let mut v = amplicon_passing_variant();
        v.pos_coverage = 4; // vd=4 → log2 = 2 exactly
        v.mean_quality = 35.47; // table-rounded 35.5
        v.hi_freq = 1.0;
        v.duprate = 0.126;
        v
    }

    fn info_field<'a>(record: &'a str, tag: &str) -> &'a str {
        let info: Vec<&str> = record.split('\t').nth(7).unwrap().split(';').collect();
        info.iter()
            .find(|p| p.starts_with(tag))
            .unwrap_or_else(|| panic!("tag {} missing in {}", tag, record))
    }

    #[test]
    fn test_qual_uses_table_rounded_qual() {
        // ports var2vcf_valid.pl:190 — multiplicand is the HALF_EVEN 1dp table
        // cell: int(2 * 35.5) = 71, where the raw mean gives int(2 * 35.47) = 70.
        let config = Config::default();
        let region = Region {
            chr: "chr20".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "T".to_string(),
            amplicon: None,
        };
        let rec = format_vcf_record(&p2_variant(), &region, &config, &[], 1.0, "0", None);
        assert_eq!(rec.split('\t').nth(5).unwrap(), "71");
        // `($vd le 1)` string-compare → QUAL 0 for vd ≤ 1.
        let mut v1 = p2_variant();
        v1.pos_coverage = 1;
        let rec1 = format_vcf_record(&v1, &region, &config, &[], 1.0, "0", None);
        assert_eq!(rec1.split('\t').nth(5).unwrap(), "0");
    }

    #[test]
    fn test_paired_qual_uses_table_rounded_qual() {
        // ports var2vcf_paired.pl:228 — vd1 > vd2 picks side 1 with its ROUNDED
        // qual: int(2 * 35.5) = 71 (raw mean would print 70).
        let config = Config::default();
        let mut t = clean_normal();
        t.pos_coverage = 4;
        t.mean_quality = 35.47;
        let mut n = clean_normal();
        n.pos_coverage = 3;
        let sv = crate::postprocess::somatic::SomaticVariant {
            variant1: t,
            variant2: Some(n),
            status: "Somatic".to_string(),
            sv1: String::new(),
            sv2: String::new(),
            end_is_variant2: false,
        };
        let pr = format_paired_vcf_record(&sv, "chr1", &config);
        assert_eq!(pr.as_str().split('\t').nth(5).unwrap(), "71");
    }

    #[test]
    fn test_paired_sample_field_into_clears_and_is_deterministic() {
        // Reuse-contract pin: `_into` clears any previous content (no
        // bleed-through between the tumor/normal columns or across rows)
        // and is deterministic for the same inputs.
        let v = clean_normal();
        let oddratio = ("0.5".to_string(), 0.5);
        let mut buf = String::new();
        buf.push_str("LONGER-PREVIOUS-ROW-CONTENT-THAT-MUST-VANISH");
        paired_sample_field_into(&mut buf, &v, "0/1", 0.5, &oddratio, false);
        let first = buf.clone();
        assert!(first.starts_with("0/1:"), "gt leads: {first}");
        assert!(
            !first.contains("VANISH"),
            "previous content cleared: {first}"
        );
        paired_sample_field_into(&mut buf, &v, "0/1", 0.5, &oddratio, false);
        assert_eq!(buf, first, "same inputs, same bytes");
        // Different genotype takes a different path through the same buffer.
        paired_sample_field_into(&mut buf, &v, "1/1", 0.5, &oddratio, false);
        assert!(buf.starts_with("1/1:"), "gt leads: {buf}");
        assert_ne!(buf, first);
    }

    #[test]
    fn test_hiaf_raw_table_cell_format() {
        // Simple writer HIAF must be Java's RAW DecimalFormat("0.0000") cell
        // (zeros kept — R echoes it verbatim via recycled colClasses):
        // oracle bytes e.g. HIAF=1.0000 / HIAF=0.0610.
        let config = Config::default();
        let region = Region {
            chr: "chr20".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "T".to_string(),
            amplicon: None,
        };
        let rec = format_vcf_record(&p2_variant(), &region, &config, &[], 1.0, "0", None);
        assert_eq!(info_field(&rec, "HIAF="), "HIAF=1.0000");
        let mut v = p2_variant();
        v.hi_freq = 0.061;
        let rec = format_vcf_record(&v, &region, &config, &[], 1.0, "0", None);
        assert_eq!(info_field(&rec, "HIAF="), "HIAF=0.0610");
        let mut z = p2_variant();
        z.hi_freq = 0.0;
        let rec = format_vcf_record(&z, &region, &config, &[], 1.0, "0", None);
        assert_eq!(info_field(&rec, "HIAF="), "HIAF=0");
    }

    #[test]
    fn test_oddratio_reciprocal_perl_stringification() {
        // ports var2vcf_valid.pl:147-151/:236 — reciprocal in (0,1) is
        // stringified like Perl's default NV interpolation (%.15g), everything
        // else passes the table string through. Expected values verified with
        // `perl -e 'printf "%.15g", 1/x'`; the first two are live oracle forms.
        assert_eq!(perl_num_str(0.54193_f64.recip()), "1.84525676747919");
        assert_eq!(perl_num_str(0.0992_f64.recip()), "10.0806451612903");
        assert_eq!(perl_num_str(4.0), "4"); // trailing zeros/dot trimmed
        assert_eq!(perl_num_str(-2.5), "-2.5");

        let config = Config::default();
        let region = Region {
            chr: "chr20".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "T".to_string(),
            amplicon: None,
        };
        let rec = format_vcf_record(&p2_variant(), &region, &config, &[], 1.0, "0.54193", None);
        assert_eq!(info_field(&rec, "ODDRATIO="), "ODDRATIO=1.84525676747919");
        // ≥1 passthrough keeps the exact table string; Inf → 0 (:147-148).
        let rec = format_vcf_record(&p2_variant(), &region, &config, &[], 1.0, "1.67216", None);
        assert_eq!(info_field(&rec, "ODDRATIO="), "ODDRATIO=1.67216");
        let rec = format_vcf_record(&p2_variant(), &region, &config, &[], 1.0, "Inf", None);
        assert_eq!(info_field(&rec, "ODDRATIO="), "ODDRATIO=0");
    }

    #[test]
    fn test_duprate_pipeline_formats() {
        // Without --fisher the oracle's DUPRATE is the 36-col cell
        // (DecimalFormat("0.0"), SimpleOutputVariant.java:223) echoed stripped by
        // R → 1dp; with --fisher Perl consumes the 38-col cell directly
        // (getRoundedValueToPrint("0.00"), :175) → 2dp.
        let mut config = Config::default();
        config.fisher = false;
        let region = Region {
            chr: "chr20".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "T".to_string(),
            amplicon: None,
        };
        let rec = format_vcf_record(&p2_variant(), &region, &config, &[], 1.0, "0", None);
        assert_eq!(info_field(&rec, "DUPRATE="), "DUPRATE=0.1");
        config.fisher = true;
        let rec = format_vcf_record(&p2_variant(), &region, &config, &[], 1.0, "0", None);
        assert_eq!(info_field(&rec, "DUPRATE="), "DUPRATE=0.13");
    }

    #[test]
    fn test_var2vcf_valid_zero_af_threshold_falls_back_to_0_02() {
        // Bug D (`var2vcf_valid.pl:21` `my $Freq = $opt_f ? $opt_f : 0.02;`):
        // when `-f 0` (`af_threshold == 0.0`) is passed in single-sample mode,
        // Perl treats `0` as falsy and uses `0.02` for the header `##FILTER=<ID=f0.02>`,
        // the `f0.02` filter check, and the `GT` `0/1 vs 0/0` threshold.
        let mut config = Config::default();
        config.af_threshold = 0.0;
        let mut header_buf = Vec::new();
        print_vcf_header(&mut header_buf, &config, &HashMap::new()).unwrap();
        let header_str = String::from_utf8(header_buf).unwrap();
        assert!(
            header_str.contains("##FILTER=<ID=f0.02,Description=\"Allele frequency < 0.02\">"),
            "var2vcf_valid header with -f 0 must emit f0.02"
        );

        let mut low_af_var = p2_variant();
        low_af_var.frequency = 0.01;
        let filters = filter_variant(&low_af_var, &config, 0, None, 1.0, "1.0");
        assert!(
            filters.contains(&"f0.02".to_string()),
            "variant with AF=0.01 under -f 0 must receive f0.02 filter"
        );
        let region = Region {
            chr: "chr20".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(110),
            gene: "T".to_string(),
            amplicon: None,
        };
        let rec = format_vcf_record(&low_af_var, &region, &config, &filters, 1.0, "1.0", None);
        let sample_col = rec.split('\t').nth(9).unwrap();
        assert!(
            sample_col.starts_with("0/0:"),
            "variant with AF=0.01 < 0.02 under -f 0 must have GT=0/0, got {sample_col}"
        );

        // Conversely, `-f 0.0` is a truthy string in Perl (`"0.0"`), so
        // `af_threshold_truthy_zero = true` preserves `0.0` instead of `0.02`.
        config.af_threshold_truthy_zero = true;
        let filters_truthy = filter_variant(&low_af_var, &config, 0, None, 1.0, "1.0");
        assert!(
            !filters_truthy.iter().any(|f| f.starts_with('f')),
            "variant with AF=0.01 under truthy -f 0.0 must not receive f0.02 filter"
        );
        let rec_truthy = format_vcf_record(
            &low_af_var,
            &region,
            &config,
            &filters_truthy,
            1.0,
            "1.0",
            None,
        );
        let sample_col_truthy = rec_truthy.split('\t').nth(9).unwrap();
        assert!(
            sample_col_truthy.starts_with("0/1:"),
            "variant with AF=0.01 >= 0.0 under truthy -f 0.0 must have GT=0/1, got {sample_col_truthy}"
        );
    }

    #[test]
    fn paired_vcf_uses_half_even_rounding() {
        // Verify paired_or_transform, hiaf_cell, and nm_str use round_half_even
        let inf_or = paired_or_transform("Inf");
        assert_eq!((inf_or.0.as_str(), inf_or.1), ("0", 0.0));
        // 1.0 / 0.8 = 1.25 -> "1.25", 1.25
        let rec_or = paired_or_transform("0.8");
        assert_eq!((rec_or.0.as_str(), rec_or.1), ("1.25", 1.25));
        assert_eq!(hiaf_cell(0.125), "0.1250");

        let config = Config::default();
        let mut t = clean_normal();
        t.mean_mapping_quality = 42.5; // half-even 0dp -> 42
        t.nm = 1.25; // half-even 1dp -> 1.2
        t.msi = 2.125; // half-even 3dp -> 2.125
        let sv = somatic_sv(t.clone(), Some(clean_normal()), "Somatic");
        let rec = format_paired_vcf_record(&sv, "chr1", &config);
        assert!(rec.as_str().contains("MSI=2.125;"), "{}", rec.as_str());
        assert!(rec.as_str().contains(":42:"), "{}", rec.as_str());

        // Two-stage .var (1dp) -> var2vcf_paired.pl sprintf("%.0f") (0dp) rounding:
        // 32.53846 -> 32.5 in .var -> 32 in paired VCF FORMAT:MQ (not 33).
        let mut t_two_stage = clean_normal();
        t_two_stage.mean_mapping_quality = 423.0 / 13.0;
        let sv_two_stage = somatic_sv(t_two_stage, Some(clean_normal()), "Somatic");
        let rec_two_stage = format_paired_vcf_record(&sv_two_stage, "chr1", &config);
        assert!(
            rec_two_stage.as_str().contains(":32:"),
            "{}",
            rec_two_stage.as_str()
        );

        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(200),
            gene: "T".to_string(),
            amplicon: None,
        };
        let single = format_vcf_record(&t, &region, &config, &[], 1.0, "1.0", None);
        assert!(single.contains("NM=1.2;"), "{single}");

        // CRISPR INFO tag when crispr_cutting_site is set (var2vcf_valid.pl:310)
        let mut crispr_cfg = Config::default();
        crispr_cfg.crispr_cutting_site = Some(50);
        let mut t_crispr = clean_normal();
        t_crispr.crispr = 3;
        let single_crispr =
            format_vcf_record(&t_crispr, &region, &crispr_cfg, &[], 1.0, "1.0", None);
        assert!(
            single_crispr.contains(";CRISPR=3;"),
            "expected ;CRISPR=3; in INFO column, got: {single_crispr}"
        );

        // --fisher VCF preserves Java getRoundedValueToPrint trailing '.' on fractional
        // values that round to integers (e.g., SN=2.0001 -> "SN=2." with --fisher vs "SN=2" without).
        let mut t_sn = clean_normal();
        t_sn.hq_lq_ratio = 2.0001;
        let mut fisher_cfg = Config::default();
        fisher_cfg.fisher = true;
        let vcf_no_fisher = format_vcf_record(&t_sn, &region, &config, &[], 1.0, "1.0", None);
        let vcf_fisher = format_vcf_record(&t_sn, &region, &fisher_cfg, &[], 1.0, "1.0", None);
        assert!(
            vcf_no_fisher.contains(";SN=2;"),
            "expected ;SN=2; without --fisher, got: {vcf_no_fisher}"
        );
        assert!(
            vcf_fisher.contains(";SN=2.;"),
            "expected ;SN=2.; with --fisher, got: {vcf_fisher}"
        );

        // omit_end_tag (-E / --no-end, ports var2vcf_valid.pl -E): omits ;END= uniformly on
        // single-sample SNVs and <DEL> SVs while keeping ;SVTYPE=DEL;SVLEN=... on SVs.
        let mut with_end_cfg = Config::default();
        with_end_cfg.omit_end_tag = false;
        let mut no_end_cfg = Config::default();
        no_end_cfg.omit_end_tag = true;
        let snv_var = amplicon_passing_variant();
        let snv_with_end =
            format_vcf_record(&snv_var, &region, &with_end_cfg, &[], 1.0, "1.0", None);
        let snv_no_end = format_vcf_record(&snv_var, &region, &no_end_cfg, &[], 1.0, "1.0", None);
        assert!(
            snv_with_end.contains(";END=100;"),
            "expected ;END=100; when omit_end_tag=false, got: {snv_with_end}"
        );
        assert!(
            !snv_no_end.contains(";END="),
            "expected ;END= omitted when omit_end_tag=true, got: {snv_no_end}"
        );

        let mut del_sv = amplicon_passing_variant();
        del_sv.var_allele = "<DEL>".to_string();
        del_sv.var_type = "Deletion".to_string();
        del_sv.start_pos = 100;
        del_sv.end_pos = 1500;
        del_sv.sv = "5-2-1".to_string();
        let sv_with_end = format_vcf_record(&del_sv, &region, &with_end_cfg, &[], 1.0, "1.0", None);
        let sv_no_end = format_vcf_record(&del_sv, &region, &no_end_cfg, &[], 1.0, "1.0", None);
        assert!(
            sv_with_end.contains(";END=1500;")
                && sv_with_end.contains(";SVTYPE=DEL;SVLEN=1400;"),
            "expected ;END=1500; and ;SVTYPE=DEL;SVLEN=1400; when omit_end_tag=false, got: {sv_with_end}"
        );
        assert!(
            !sv_no_end.contains(";END=") && sv_no_end.contains(";SVTYPE=DEL;SVLEN=1400;"),
            "expected ;END= omitted while keeping ;SVTYPE=DEL;SVLEN=1400; when omit_end_tag=true, got: {sv_no_end}"
        );
    }

    #[test]
    fn test_perl_reorder_chromosomes_matches_var2vcf_sub_reorder() {
        // 1. hg19 lexicographic input: GL000192.1 has digits 1921 and no '_', so it enters
        //    @chrn after chr10 (10) and before chrX; chrUn_gl000211 has '_' so it stays in @nonchrn.
        let hg19 = [
            "chr1",
            "chr10",
            "chr2",
            "chrX",
            "chrY",
            "chrM",
            "GL000192.1",
            "chrUn_gl000211",
        ];
        assert_eq!(
            perl_reorder_chromosomes(&hg19),
            vec![
                "chr1",
                "chr2",
                "chr10",
                "GL000192.1",
                "chrX",
                "chrY",
                "chrM",
                "chrUn_gl000211"
            ]
        );

        // 2. b37 input: hs37d5 (375) and GL000220.1 (2201) precede X/Y/MT; NC_007605 has '_' so tail.
        let b37 = [
            "1",
            "10",
            "2",
            "X",
            "Y",
            "MT",
            "GL000220.1",
            "hs37d5",
            "NC_007605",
        ];
        assert_eq!(
            perl_reorder_chromosomes(&b37),
            vec![
                "1",
                "2",
                "10",
                "hs37d5",
                "GL000220.1",
                "X",
                "Y",
                "MT",
                "NC_007605"
            ]
        );

        // 3. Upstream Perl elsif drop: when both X and chrX (or Y/chrY, MT/chrM) exist,
        //    chrX/chrY/chrM are dropped; 'M' is not skipped by theelsif block and stays in @nonchrn;
        //    ties in @chrn ('chr1' vs '1') preserve first-seen order.
        let drop_case = ["X", "chrX", "chr1", "1", "Y", "chrY", "MT", "chrM", "M"];
        assert_eq!(
            perl_reorder_chromosomes(&drop_case),
            vec!["chr1", "1", "X", "Y", "MT", "M"]
        );

        // 4. 25-digit chromosome name: exact arbitrary-length digit comparison without u64 overflow.
        let huge_digits = ["chr9999999999999999999999999", "chr10", "chr2"];
        assert_eq!(
            perl_reorder_chromosomes(&huge_digits),
            vec!["chr2", "chr10", "chr9999999999999999999999999"]
        );
    }

    #[test]
    fn test_paired_normal_clean_boolean_and_static_filter_tags() {
        let config = Config::default();

        // 1. Both tumor and normal pass filters -> PASS
        let mut clean_tumor = clean_normal();
        clean_tumor.var_type = "SNV".to_string();
        let rec_both_clean = format_paired_vcf_record(
            &somatic_sv(clean_tumor.clone(), Some(clean_normal()), "Somatic"),
            "chr1",
            &config,
        );
        assert_eq!(rec_both_clean.filter_str(), "PASS");

        // 2. Tumor clean (filters.is_empty()), normal fails (!normal_clean) -> PASS (no -M)
        let mut dirty_normal = clean_normal();
        dirty_normal.pos_coverage = 1; // vd2 < 3 -> !normal_clean
        let rec_tumor_clean_normal_dirty = format_paired_vcf_record(
            &somatic_sv(clean_tumor, Some(dirty_normal.clone()), "Somatic"),
            "chr1",
            &config,
        );
        assert_eq!(rec_tumor_clean_normal_dirty.filter_str(), "PASS");

        // 3. Tumor fails d5, v3, Q0 and normal fails (!normal_clean) -> emits d5;v3;Q0
        let mut low_tumor = clean_normal();
        low_tumor.var_type = "SNV".to_string();
        low_tumor.total_pos_cov = 4; // < d5
        low_tumor.pos_coverage = 2; // < v3
        low_tumor.pstd = true;
        low_tumor.mean_mapping_quality = 5.0; // SNV with mapq < 10 -> Q0
        let rec_both_dirty = format_paired_vcf_record(
            &somatic_sv(low_tumor.clone(), Some(dirty_normal), "Germline"),
            "chr1",
            &config,
        );
        assert_eq!(rec_both_dirty.filter_str(), "d5;v3;Q0");

        // 4. Same failing tumor with a clean normal (normal_clean = true) -> germline passthrough PASS
        let rec_passthrough = format_paired_vcf_record(
            &somatic_sv(low_tumor, Some(clean_normal()), "Germline"),
            "chr1",
            &config,
        );
        assert_eq!(rec_passthrough.filter_str(), "PASS");
    }

    #[test]
    fn format_filtered_vcf_record_into_matches_filter_variant_plus_format_vcf_record_into() {
        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(90),
            end: crate::types::Pos1(150),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        config.af_threshold = 0.025;

        // 1. PASS variant
        let pass_var = amplicon_passing_variant();
        let expected_pass_filters = filter_variant(&pass_var, &config, 0, None, 1.0, "1.5");
        assert!(expected_pass_filters.is_empty());
        let mut expected_buf = String::new();
        format_vcf_record_into(
            &mut expected_buf,
            &pass_var,
            &region,
            &config,
            &expected_pass_filters,
            1.0,
            "1.5",
            None,
        );
        let mut actual_buf = String::new();
        format_filtered_vcf_record_into(
            &mut actual_buf,
            &pass_var,
            &region,
            &config,
            0,
            1.0,
            "1.5",
            None,
        );
        assert_eq!(actual_buf, expected_buf);
        assert_eq!(actual_buf.split('\t').nth(6), Some("PASS"));

        // 2. Multi-tag ;-joined filters including f{min_af}
        let mut multi_var = Variant::new();
        multi_var.ref_allele = "AT".to_string();
        multi_var.var_allele = "A".to_string();
        multi_var.start_pos = 100;
        multi_var.end_pos = 101;
        multi_var.total_pos_cov = 2; // d3
        multi_var.pos_coverage = 1;
        multi_var.hicnt = 1; // v2 (1 * 0.01 < 0.5)
        multi_var.hicov = 2;
        multi_var.hi_freq = 0.01;
        multi_var.frequency = 0.01; // f0.025
        multi_var.mean_position = 4.0; // p8
        multi_var.pstd = false; // pSTD
        multi_var.mean_quality = 15.0; // q22.5
        multi_var.mean_mapping_quality = 5.0; // Q10
        multi_var.hq_lq_ratio = 1.0; // SN1.5
        multi_var.nm = 6.0; // NM5.25
        multi_var.msi = 13.5; // MSI12 + LongMSI
        multi_var.msint = 1;
        multi_var.bias = "2;1".to_string(); // Bias
        multi_var.var_type = "Deletion".to_string();

        let expected_multi_filters = filter_variant(&multi_var, &config, 0, None, 0.001, "6.0");
        assert_eq!(
            expected_multi_filters.join(";"),
            "d3;v2;f0.025;p8;pSTD;q22.5;Q10;SN1.5;NM5.25;MSI12;Bias;LongMSI"
        );
        format_vcf_record_into(
            &mut expected_buf,
            &multi_var,
            &region,
            &config,
            &expected_multi_filters,
            0.001,
            "6.0",
            None,
        );
        format_filtered_vcf_record_into(
            &mut actual_buf,
            &multi_var,
            &region,
            &config,
            0,
            0.001,
            "6.0",
            None,
        );
        assert_eq!(actual_buf, expected_buf);

        // 3. Cluster0bp + AMPBIAS ordering (Cluster0bp requires filters.is_empty() before AMPBIAS)
        let amp_bias = Some((1, 3, 0, true));
        let expected_cluster_amp = filter_variant(&pass_var, &config, 105, amp_bias, 1.0, "1.5");
        assert_eq!(expected_cluster_amp, vec!["Cluster0bp", "AMPBIAS"]);
        format_vcf_record_into(
            &mut expected_buf,
            &pass_var,
            &region,
            &config,
            &expected_cluster_amp,
            1.0,
            "1.5",
            amp_bias,
        );
        format_filtered_vcf_record_into(
            &mut actual_buf,
            &pass_var,
            &region,
            &config,
            105,
            1.0,
            "1.5",
            amp_bias,
        );
        assert_eq!(actual_buf, expected_buf);
        assert_eq!(actual_buf.split('\t').nth(6), Some("Cluster0bp;AMPBIAS"));

        // When an earlier filter fires, Cluster0bp is suppressed while AMPBIAS still appends
        let mut low_af_snv = pass_var.clone();
        low_af_snv.frequency = 0.01;
        let expected_no_cluster = filter_variant(&low_af_snv, &config, 105, amp_bias, 1.0, "1.5");
        assert_eq!(expected_no_cluster, vec!["f0.025", "AMPBIAS"]);
        format_vcf_record_into(
            &mut expected_buf,
            &low_af_snv,
            &region,
            &config,
            &expected_no_cluster,
            1.0,
            "1.5",
            amp_bias,
        );
        format_filtered_vcf_record_into(
            &mut actual_buf,
            &low_af_snv,
            &region,
            &config,
            105,
            1.0,
            "1.5",
            amp_bias,
        );
        assert_eq!(actual_buf, expected_buf);

        check_valid_pstd_and_pass_only_semantics(&region);
    }

    fn check_valid_pstd_and_pass_only_semantics(region: &Region) {
        // 1. parse_sv_counts (`var2vcf_valid.pl:198`: `defined($tamp) && $tamp =~ /-/`)
        assert_eq!(parse_sv_counts(""), (0, 0));
        assert_eq!(parse_sv_counts("0"), (0, 0));
        assert_eq!(parse_sv_counts("5"), (0, 0));
        assert_eq!(parse_sv_counts("0-9-1"), (0, 9));
        assert_eq!(parse_sv_counts("5-2-1"), (5, 2));
        assert_eq!(parse_sv_counts("3-4"), (3, 4));

        // 2. `--no-pstd-filter` (`config.vcf_pstd_filter = false`, `var2vcf_valid.pl:30,161`)
        let mut pstd_var = amplicon_passing_variant();
        pstd_var.pstd = false;
        pstd_var.frequency = 0.10; // < 0.35 -> triggers pSTD by default in non-amplicon mode

        let default_cfg = Config::default();
        assert_eq!(
            filter_variant(&pstd_var, &default_cfg, 0, None, 1.0, "1.5"),
            vec!["pSTD"]
        );

        let mut no_pstd_cfg = Config::default();
        no_pstd_cfg.vcf_pstd_filter = false;
        assert!(
            filter_variant(&pstd_var, &no_pstd_cfg, 0, None, 1.0, "1.5").is_empty(),
            "vcf_pstd_filter=false must suppress pSTD in non-amplicon mode"
        );

        // 3. `try_format_filtered_vcf_record_into` with `config.vcf_pass_only` (`var2vcf_valid.pl:187`)
        let mut pass_only_cfg = Config::default();
        pass_only_cfg.vcf_pass_only = true;
        let mut buf = "STALE-CONTENT".to_string();
        let emitted_filtered = try_format_filtered_vcf_record_into(
            &mut buf,
            &pstd_var,
            region,
            &pass_only_cfg,
            0,
            1.0,
            "1.5",
            None,
        );
        assert!(!emitted_filtered, "non-PASS variant must be suppressed");
        assert!(
            buf.is_empty(),
            "suppressed call must clear destination buffer"
        );

        // Amplicon AMPBIAS non-PASS variant is also suppressed under `vcf_pass_only = true`
        let pass_var = amplicon_passing_variant();
        buf.push_str("STALE");
        let emitted_amp_filtered = try_format_filtered_vcf_record_into(
            &mut buf,
            &pass_var,
            region,
            &pass_only_cfg,
            0,
            1.0,
            "1.5",
            Some((1, 4, 0, true)),
        );
        assert!(!emitted_amp_filtered);
        assert!(buf.is_empty());

        // Combining `vcf_pstd_filter = false` + `vcf_pass_only = true` promotes `pstd_var` to PASS
        let mut pass_no_pstd_cfg = pass_only_cfg.clone();
        pass_no_pstd_cfg.vcf_pstd_filter = false;
        let emitted_pass = try_format_filtered_vcf_record_into(
            &mut buf,
            &pstd_var,
            region,
            &pass_no_pstd_cfg,
            0,
            1.0,
            "1.5",
            None,
        );
        assert!(emitted_pass);
        assert_eq!(buf.split('\t').nth(6), Some("PASS"));

        // 4. Mid-allele `<dup-N>` (`#342` family, `var2vcf_valid.pl:199` `$alt =~ /</`):
        //    when `--min-splitreads 0` admits a mid-allele `<dup-12>` variant (`splitreads = 0`),
        //    `format_vcf_record_with_filter` must emit `;SVTYPE=...dup-12...;SVLEN=...` stripping `<`/`>`.
        let mut mid_dup_var = amplicon_passing_variant();
        mid_dup_var.start_pos = 100;
        mid_dup_var.end_pos = 112;
        mid_dup_var.ref_allele = "A".to_string();
        mid_dup_var.var_allele = "ACGT<dup-12>TTAA".to_string();
        mid_dup_var.sv = "0-0-0".to_string();
        let mid_dup_line =
            format_vcf_record(&mid_dup_var, region, &default_cfg, &[], 1.0, "1.5", None);
        assert!(
            mid_dup_line.contains(";SVTYPE=ACGTdup-12TTAA;SVLEN=12;SPLITREAD=0;SPANPAIR=0"),
            "mid-allele <dup-12> must emit SVTYPE/SVLEN matching Perl `$alt =~ /</`: {mid_dup_line}"
        );
    }

    #[test]
    fn format_paired_vcf_record_normal_clean_pass_override_and_multi_filter_parity() {
        let mut config = Config::default();
        config.af_threshold = 0.02;

        // Build tumor failing all tumor gates + tumor Bias + MSI12 (both gate 1 and gate 2)
        let mut tumor = Variant::new();
        tumor.ref_allele = "ATG".to_string();
        tumor.var_allele = "A".to_string(); // len_diff = 2 == msint
        tumor.start_pos = 100;
        tumor.end_pos = 102;
        tumor.total_pos_cov = 4; // d5
        tumor.pos_coverage = 2; // v3
        tumor.ref_fwd_cov = 20;
        tumor.ref_rev_cov = 0;
        tumor.fwd = 0;
        tumor.rev = 2; // p1 < 0.01, OR = Inf -> 0.0 => Bias
        tumor.frequency = 0.01; // f0.02 and < 0.35
        tumor.mean_position = 4.0; // p8
        tumor.pstd = false; // pSTD (vd1 < 5)
        tumor.mean_quality = 15.0; // q22.5
        tumor.mean_mapping_quality = 5.0; // Q0 when SNV
        tumor.hq_lq_ratio = 1.0; // SN1.5
        tumor.nm = 6.0; // NM5.25
        tumor.bias = "2;1".to_string();
        tumor.msi = 15.0; // > 12.0 with msint = 2 => MSI12
        tumor.msint = 2;
        tumor.var_type = "SNV".to_string();

        // Normal passes all normal_clean checks (dp>=5, vd>=3, af>=0.02, pmean>=8, pstd,
        // qual>=22.5, mapq>=0, sn>=1.5, nm<5.25) while ALSO triggering normal-side Bias!
        let mut clean_norm_with_bias = clean_normal();
        clean_norm_with_bias.total_pos_cov = 40;
        clean_norm_with_bias.pos_coverage = 10;
        clean_norm_with_bias.ref_fwd_cov = 30;
        clean_norm_with_bias.ref_rev_cov = 0;
        clean_norm_with_bias.fwd = 0;
        clean_norm_with_bias.rev = 10; // p2 < 0.01, OR = Inf -> 0.0
        clean_norm_with_bias.frequency = 0.25; // >= 0.02 and < 0.35
        clean_norm_with_bias.bias = "2;1".to_string();

        // 1. normal_clean = true overrides tumor filters + Bias + MSI12 to PASS
        let rec_override = format_paired_vcf_record(
            &somatic_sv(
                tumor.clone(),
                Some(clean_norm_with_bias.clone()),
                "Germline",
            ),
            "chr1",
            &config,
        );
        assert_eq!(rec_override.filter_str(), "PASS");

        // 2. normal_clean = false preserves exact ;-joined order without duplicate Bias or MSI12
        let mut dirty_norm_with_bias = clean_norm_with_bias;
        dirty_norm_with_bias.nm = 6.0; // >= 5.25 -> normal_clean = false
        let rec_multi = format_paired_vcf_record(
            &somatic_sv(tumor, Some(dirty_norm_with_bias.clone()), "Germline"),
            "chr1",
            &config,
        );
        assert_eq!(
            rec_multi.filter_str(),
            "d5;v3;f0.02;p8;pSTD;q22.5;Q0;SN1.5;NM5.25;Bias;MSI12"
        );

        // 3. normal-only Bias + MSI gate 2 only (StrongSomatic with tiny ptot suppresses MSI gate 1)
        let mut clean_tumor_msi2 = clean_normal();
        clean_tumor_msi2.ref_allele = "ATG".to_string();
        clean_tumor_msi2.var_allele = "A".to_string(); // len_diff = 2 == msint
        clean_tumor_msi2.start_pos = 100;
        clean_tumor_msi2.end_pos = 102;
        clean_tumor_msi2.total_pos_cov = 100;
        clean_tumor_msi2.pos_coverage = 30;
        clean_tumor_msi2.frequency = 0.30; // < 0.35
        clean_tumor_msi2.msi = 15.0;
        clean_tumor_msi2.msint = 2;
        clean_tumor_msi2.var_type = "Deletion".to_string();

        let mut deep_dirty_norm_bias = dirty_dirty_normal_for_ss(&dirty_norm_with_bias);
        deep_dirty_norm_bias.frequency = 0.01; // < 0.02 (normal_clean = false) and < 0.35 (MSI gate 2)
        let rec_norm_bias_msi2 = format_paired_vcf_record(
            &somatic_sv(
                clean_tumor_msi2,
                Some(deep_dirty_norm_bias),
                "StrongSomatic",
            ),
            "chr1",
            &config,
        );
        assert_eq!(rec_norm_bias_msi2.filter_str(), "Bias;MSI12");

        check_paired_strict_and_pval_semantics();
    }

    fn check_paired_strict_and_pval_semantics() {
        let lengths = HashMap::new();

        // 1. Header `##FILTER=<ID=P$PVAL,...>` formatting (`var2vcf_paired.pl:27,86`)
        let render_hdr = |cfg: &Config| -> String {
            let mut b = Vec::new();
            print_paired_vcf_header(&mut b, cfg, &lengths).unwrap();
            String::from_utf8(b).unwrap()
        };
        assert!(render_hdr(&Config::default())
            .contains(r#"##FILTER=<ID=P0.05,Description="Not significant with p-value > 0.05">"#));

        let mut raw_sci_cfg = Config::default();
        raw_sci_cfg.vcf_pval = 1e-5;
        raw_sci_cfg.vcf_pval_raw = Some("1e-5".to_string());
        assert!(render_hdr(&raw_sci_cfg)
            .contains(r#"##FILTER=<ID=P1e-5,Description="Not significant with p-value > 1e-5">"#));

        let mut raw_trailing_cfg = Config::default();
        raw_trailing_cfg.vcf_pval = 0.05;
        raw_trailing_cfg.vcf_pval_raw = Some("0.050".to_string());
        assert!(render_hdr(&raw_trailing_cfg).contains(
            r#"##FILTER=<ID=P0.050,Description="Not significant with p-value > 0.050">"#
        ));

        let mut prog_sci_cfg = Config::default();
        prog_sci_cfg.vcf_pval = 1e-5;
        prog_sci_cfg.vcf_pval_raw = None;
        assert!(render_hdr(&prog_sci_cfg).contains(
            r#"##FILTER=<ID=P1e-05,Description="Not significant with p-value > 1e-05">"#
        ));

        // 2. Germline rescue disabled under `vcf_paired_strict = true` (`var2vcf_paired.pl:217-219`)
        let mut low_qual_tumor = clean_normal();
        low_qual_tumor.var_type = "SNV".to_string();
        low_qual_tumor.mean_quality = 12.0; // < q22.5; equal tumor/normal counts -> ptot = 1.0 > 0.05
        let germline_sv = somatic_sv(low_qual_tumor, Some(clean_normal()), "Germline");
        assert_eq!(
            format_paired_vcf_record(&germline_sv, "chr1", &Config::default()).filter_str(),
            "PASS"
        );
        let mut strict_cfg = Config::default();
        strict_cfg.vcf_paired_strict = true;
        let strict_germline = format_paired_vcf_record(&germline_sv, "chr1", &strict_cfg);
        assert_eq!(strict_germline.filter_str(), "q22.5;P0.05");

        // Verbatim `P1e-5` vs programmatic `P1e-05` in record FILTER under `vcf_paired_strict = true`
        raw_sci_cfg.vcf_paired_strict = true;
        prog_sci_cfg.vcf_paired_strict = true;
        assert_eq!(
            format_paired_vcf_record(&germline_sv, "chr1", &raw_sci_cfg).filter_str(),
            "q22.5;P1e-5"
        );
        assert_eq!(
            format_paired_vcf_record(&germline_sv, "chr1", &prog_sci_cfg).filter_str(),
            "q22.5;P1e-05"
        );

        // 3. StrongSomatic `P$PVAL` exceptions (`var2vcf_paired.pl:199-200`):
        //    (ptot < 0.25 && af1 > 0.1) || (ptot < 0.5 && af1 > 0.20) || (ptot < 0.15 && af1 > 0.05)
        // Construct counts giving moderate ptot in (0.05, 0.45):
        // tumor: dp=20, vd=6 (tref=14, af1=0.30); normal: dp=20, vd=3 (rref=17) -> Fisher one-sided p ~ 0.22
        let mut ss_tumor = clean_normal();
        ss_tumor.var_type = "SNV".to_string();
        ss_tumor.total_pos_cov = 20;
        ss_tumor.pos_coverage = 6;
        ss_tumor.ref_fwd_cov = 7;
        ss_tumor.ref_rev_cov = 7;
        ss_tumor.fwd = 3;
        ss_tumor.rev = 3;
        ss_tumor.frequency = 0.30; // af1 > 0.20 satisfies `(ptot < 0.5 && af1 > 0.20)`
        let mut ss_normal = clean_normal();
        ss_normal.total_pos_cov = 20;
        ss_normal.pos_coverage = 3;
        ss_normal.ref_fwd_cov = 9;
        ss_normal.ref_rev_cov = 8;
        ss_normal.fwd = 2;
        ss_normal.rev = 1;
        ss_normal.frequency = 0.15;

        let ss_exc_rec = format_paired_vcf_record(
            &somatic_sv(ss_tumor.clone(), Some(ss_normal.clone()), "StrongSomatic"),
            "chr1",
            &strict_cfg,
        );
        assert_eq!(
            ss_exc_rec.filter_str(),
            "PASS",
            "StrongSomatic with ptot<0.5 and af1>0.20 must suppress P0.05: {}",
            ss_exc_rec.as_str()
        );
        assert!(!ss_exc_rec.as_str().contains(";P0.01Likely"));
        assert!(!ss_exc_rec.as_str().contains(";InDelLikely"));

        // Same counts with LikelySomatic do NOT get the StrongSomatic exception -> P0.05 (and NOT P0.01Likely!)
        let ls_no_exc = format_paired_vcf_record(
            &somatic_sv(ss_tumor, Some(ss_normal), "LikelySomatic"),
            "chr1",
            &strict_cfg,
        );
        assert_eq!(ls_no_exc.filter_str(), "P0.05");
        assert!(
            !ls_no_exc.as_str().contains(";P0.01Likely"),
            "ptot > PVAL enters `if` branch so `elsif` P0.01Likely must not fire"
        );

        // 4. `;P0.01Likely` (`var2vcf_paired.pl:201-202, 239`):
        //    `status =~ /LikelySomatic/ && ptot > 0.05 / 5**vd2` when `ptot <= PVAL`.
        //    tumor: dp=40, vd=12 (tref=28); normal: dp=40, vd=2 (rref=38) -> ptot ~ 0.0033 <= 0.05,
        //    while 0.05 / 5^2 = 0.002 < 0.0033 -> sets `;P0.01Likely` right after `RSEQ=...`.
        let mut ls_tumor = clean_normal();
        ls_tumor.var_type = "SNV".to_string();
        ls_tumor.total_pos_cov = 40;
        ls_tumor.pos_coverage = 12;
        ls_tumor.ref_fwd_cov = 14;
        ls_tumor.ref_rev_cov = 14;
        ls_tumor.fwd = 6;
        ls_tumor.rev = 6;
        ls_tumor.frequency = 0.30;
        ls_tumor.right_seq = "ACGT".to_string();
        let mut ls_norm = clean_normal();
        ls_norm.total_pos_cov = 40;
        ls_norm.pos_coverage = 2;
        ls_norm.ref_fwd_cov = 19;
        ls_norm.ref_rev_cov = 19;
        ls_norm.fwd = 1;
        ls_norm.rev = 1;
        ls_norm.frequency = 0.05;

        let p_likely_rec = format_paired_vcf_record(
            &somatic_sv(ls_tumor.clone(), Some(ls_norm), "LikelySomatic"),
            "chr1",
            &strict_cfg,
        );
        assert_eq!(p_likely_rec.filter_str(), "PASS");
        assert!(
            p_likely_rec
                .as_str()
                .contains(";RSEQ=ACGT;P0.01Likely\tGT:"),
            "expected ;P0.01Likely immediately after RSEQ before FORMAT: {}",
            p_likely_rec.as_str()
        );

        // Ultra-deep `vd2 = 10_000`: `5.0_f64.powi(10_000)` saturates cleanly to `+inf`
        // (`0.05 / +inf == 0.0`), never panicking. With `vcf_pval = 1.0`, `ptot_r > 0.0` sets `;P0.01Likely`.
        let mut deep_pval_cfg = strict_cfg.clone();
        deep_pval_cfg.vcf_pval = 1.0;
        let mut ultra_norm = clean_normal();
        ultra_norm.total_pos_cov = 20_000;
        ultra_norm.pos_coverage = 10_000;
        let ultra_rec = format_paired_vcf_record(
            &somatic_sv(ls_tumor, Some(ultra_norm), "LikelySomatic"),
            "chr1",
            &deep_pval_cfg,
        );
        assert!(ultra_rec.as_str().contains(";P0.01Likely\tGT:"));

        // 5. `;InDelLikely` (`var2vcf_paired.pl:203-205, 239`):
        //    `status =~ /Likely/ && type ne "SNV"` unless `length(ref) <= 2 && length(alt) <= 2`.
        let mut indel_tumor = clean_normal();
        indel_tumor.ref_allele = "ACGT".to_string(); // len 4 > 2
        indel_tumor.var_allele = "A".to_string();
        indel_tumor.var_type = "Deletion".to_string();
        indel_tumor.total_pos_cov = 100;
        indel_tumor.pos_coverage = 40;
        indel_tumor.frequency = 0.40;
        indel_tumor.right_seq = "TTAA".to_string();
        let mut clean_zero_norm = clean_normal();
        clean_zero_norm.total_pos_cov = 100;
        clean_zero_norm.pos_coverage = 0; // ptot tiny (0.0 after 5dp) <= 0.05/5^0 = 0.05
        clean_zero_norm.frequency = 0.0;

        let indel_likely_rec = format_paired_vcf_record(
            &somatic_sv(
                indel_tumor.clone(),
                Some(clean_zero_norm.clone()),
                "LikelyLOH",
            ),
            "chr1",
            &strict_cfg,
        );
        assert_eq!(indel_likely_rec.filter_str(), "PASS");
        assert!(
            indel_likely_rec
                .as_str()
                .contains(";RSEQ=TTAA;InDelLikely\tGT:"),
            "expected ;InDelLikely immediately after RSEQ: {}",
            indel_likely_rec.as_str()
        );

        // Short indel (`ref.len() <= 2 && alt.len() <= 2`) must NOT set `;InDelLikely`
        indel_tumor.ref_allele = "AC".to_string();
        indel_tumor.var_allele = "A".to_string();
        let short_indel_rec = format_paired_vcf_record(
            &somatic_sv(indel_tumor, Some(clean_zero_norm), "LikelyLOH"),
            "chr1",
            &strict_cfg,
        );
        assert!(!short_indel_rec.as_str().contains(";InDelLikely"));
    }

    fn dirty_dirty_normal_for_ss(base: &Variant) -> Variant {
        let mut v = base.clone();
        v.total_pos_cov = 200;
        v.pos_coverage = 1; // ptot vs (30, 70) is < 0.0005
        v.ref_fwd_cov = 190;
        v.ref_rev_cov = 0;
        v.fwd = 0;
        v.rev = 9;
        v
    }

    #[test]
    fn paired_streamer_cluster_tag_and_zero_variant_bias() {
        let mut config = Config::default();

        // 1. ZERO_VARIANT singleton & normal = None "0;0" BIAS preservation
        let z1 = ZERO_VARIANT.get_or_init(crate::variations::Variant::new);
        let z2 = ZERO_VARIANT.get_or_init(crate::variations::Variant::new);
        assert!(std::ptr::eq(z1, z2), "ZERO_VARIANT must initialize once");
        assert_eq!(z1.bias, "0");

        let mut explicit_zero = crate::variations::Variant::new();
        explicit_zero.bias = "0;0".to_string();

        for fisher in [false, true] {
            config.fisher = fisher;
            let mut t = clean_normal();
            t.start_pos = 100;
            t.end_pos = 100;
            t.var_type = "SNV".to_string();

            let rec_none = format_paired_vcf_record(
                &somatic_sv(t.clone(), None, "StrongSomatic"),
                "chr1",
                &config,
            );
            let rec_zero_bias = format_paired_vcf_record(
                &somatic_sv(t, Some(explicit_zero.clone()), "StrongSomatic"),
                "chr1",
                &config,
            );
            assert_eq!(rec_none.as_str(), rec_zero_bias.as_str());
            assert_eq!(rec_none.filter_str(), rec_zero_bias.filter_str());
            assert_eq!(rec_none.filter_offsets(), rec_zero_bias.filter_offsets());

            let normal_sample = rec_none
                .as_str()
                .rsplit('\t')
                .next()
                .expect("normal sample column");
            assert_eq!(
                normal_sample.split(':').nth(7),
                Some("0,0"),
                "normal=None must preserve 0;0 -> 0,0 BIAS in FORMAT: {normal_sample}"
            );
        }

        // 2. PairedRecordStreamer Cluster{opt_c}bp formatting in prev.filter and cur.filter
        config.fisher = false;
        let mk = |start: i32, status: &str| {
            let mut t = clean_normal();
            t.start_pos = start;
            t.end_pos = start;
            t.var_type = "SNV".to_string();
            format_paired_vcf_record(&somatic_sv(t, None, status), "chr1", &config)
        };

        let mut pair_lines = Vec::new();
        let mut streamer10 = PairedRecordStreamer::new(10);
        streamer10
            .push(mk(100, "StrongSomatic"), |l| -> Result<(), ()> {
                pair_lines.push(l.to_string());
                Ok(())
            })
            .unwrap();
        streamer10
            .push(mk(106, "LikelySomatic"), |l| -> Result<(), ()> {
                pair_lines.push(l.to_string());
                Ok(())
            })
            .unwrap();
        streamer10
            .finish(|l| -> Result<(), ()> {
                pair_lines.push(l.to_string());
                Ok(())
            })
            .unwrap();
        assert_eq!(pair_lines.len(), 2);
        assert_eq!(pair_lines[0].split('\t').nth(6), Some("Cluster10bp"));
        assert_eq!(pair_lines[1].split('\t').nth(6), Some("Cluster10bp"));

        let lines15 = finish_paired_records(
            vec![
                mk(100, "StrongSomatic"),
                mk(110, "StrongSomatic"),
                mk(120, "StrongSomatic"),
                mk(125, "Germline"),
            ],
            15,
        );
        let filters15: Vec<&str> = lines15
            .iter()
            .map(|l| l.split('\t').nth(6).unwrap())
            .collect();
        assert_eq!(
            filters15,
            vec!["Cluster15bp", "Cluster15bp", "PASS", "PASS"]
        );

        check_single_buffer_paired_record_and_cluster_rewrite();
    }

    fn check_single_buffer_paired_record_and_cluster_rewrite() {
        #[cfg(target_pointer_width = "64")]
        {
            assert_eq!(std::mem::size_of::<PairedVcfRecord>(), 48);
            assert_eq!(std::mem::size_of::<FilterTag<'_>>(), 16);
        }

        let cfg = Config::default();
        let mut t = clean_normal();
        t.start_pos = 100;
        t.end_pos = 100;
        t.var_type = "SNV".to_string();
        let mut rec = format_paired_vcf_record(&somatic_sv(t, None, "StrongSomatic"), "chr1", &cfg);
        assert!(rec.is_pass());
        assert_eq!(rec.filter_str(), "PASS");
        let (fs0, fe0) = rec.filter_offsets();
        assert_eq!(&rec.as_str()[fs0..fe0], "PASS");

        // Length-changing FILTER rewrite ("PASS" 4B -> "Cluster10bp" 11B) preserves prefix & rest
        let before_suffix = rec.as_str()[fe0..].to_string();
        let before_prefix = rec.as_str()[..fs0].to_string();
        rec.set_cluster_filter("Cluster10bp");
        assert!(!rec.is_pass());
        assert_eq!(rec.filter_str(), "Cluster10bp");
        let (fs1, fe1) = rec.filter_offsets();
        assert_eq!((fs1, fe1), (fs0, fs0 + 11));
        assert_eq!(
            rec.as_str(),
            format!("{before_prefix}Cluster10bp{before_suffix}")
        );

        // `from_parts_checked` rejects inverted, out-of-bounds, and mid-UTF8 spans
        let valid_line = rec.as_str().to_string();
        assert!(
            PairedVcfRecord::from_parts_checked(valid_line.clone(), fs1, fe1, 100, true, true)
                .is_some()
        );
        assert!(
            PairedVcfRecord::from_parts_checked(valid_line.clone(), fe1, fs1, 100, true, true)
                .is_none(),
            "inverted span must be rejected"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(
                valid_line.clone(),
                fs1,
                valid_line.len() + 1,
                100,
                true,
                true,
            )
            .is_none(),
            "out-of-bounds span must be rejected"
        );
        let utf8_line = "chr1\t100\t.\tA\tG\t35\tµPASS\tSTATUS=Somatic\tGT\t0/1\t0/0".to_string();
        let mu_pos = utf8_line.find('µ').unwrap();
        assert_eq!('µ'.len_utf8(), 2);
        assert!(
            PairedVcfRecord::from_parts_checked(
                utf8_line.clone(),
                mu_pos,
                mu_pos + 1,
                100,
                true,
                true,
            )
            .is_none(),
            "mid-UTF8 end offset must be rejected"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(
                utf8_line,
                mu_pos + 1,
                mu_pos + 6,
                100,
                true,
                true,
            )
            .is_none(),
            "mid-UTF8 start offset must be rejected"
        );
        // Column-aligned non-empty FILTER span checks (VCF column 7 = index 6, 10 tabs, no \r/\n)
        assert!(
            PairedVcfRecord::from_parts_checked(valid_line.clone(), 0, 4, 100, true, true)
                .is_none(),
            "CHROM column (0..4) must be rejected as FILTER span"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(valid_line.clone(), fs1, fs1, 100, true, true)
                .is_none(),
            "empty FILTER span must be rejected"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(valid_line.clone(), fs1, fe1 - 1, 100, true, true)
                .is_none(),
            "truncated FILTER span without trailing tab must be rejected"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(valid_line.clone(), 5, 8, 100, true, true)
                .is_none(),
            "POS column (col 2, 1 tab before) must be rejected as FILTER span"
        );
        let tab_in_filter_end = valid_line[fe1 + 1..]
            .find('\t')
            .map(|rel| fe1 + 1 + rel)
            .unwrap();
        assert!(
            PairedVcfRecord::from_parts_checked(
                valid_line.clone(),
                fs1,
                tab_in_filter_end,
                100,
                true,
                true,
            )
            .is_none(),
            "span crossing tab into INFO column must be rejected"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(
                "chr1\t100\t.\tA\tG\t35\tPASS\tSTATUS=Somatic".to_string(),
                18,
                22,
                100,
                true,
                true,
            )
            .is_none(),
            "8-column line (< 10 tabs) must be rejected"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(
                format!("{valid_line}\textra"),
                fs1,
                fe1,
                100,
                true,
                true,
            )
            .is_none(),
            "12-column line (> 10 tabs) must be rejected"
        );
        assert!(
            PairedVcfRecord::from_parts_checked(
                format!("{valid_line}\n"),
                fs1,
                fe1,
                100,
                true,
                true,
            )
            .is_none(),
            "line containing trailing newline must be rejected"
        );

        // Direct library-API sanitization (`vcf_pval_raw` injection / NaN `vcf_pval`) &
        // extreme `i32` bounds in `filter_variant` / `format_paired_vcf_record`
        let mut bad_lib_cfg = Config::default();
        bad_lib_cfg.vcf_paired_strict = true;
        bad_lib_cfg.vcf_pval = f64::NAN;
        bad_lib_cfg.vcf_pval_raw = Some("0.01\n##evil".to_string());
        let mut hdr_buf = Vec::new();
        print_paired_vcf_header(&mut hdr_buf, &bad_lib_cfg, &HashMap::new()).unwrap();
        let hdr_str = String::from_utf8(hdr_buf).unwrap();
        assert!(hdr_str.contains("##FILTER=<ID=P0.05,"));
        assert!(!hdr_str.contains("##evil"));

        let mut extreme = clean_normal();
        extreme.start_pos = 0;
        extreme.end_pos = i32::MAX;
        extreme.total_pos_cov = 10;
        extreme.pos_coverage = i32::MIN;
        extreme.ref_fwd_cov = i32::MAX;
        extreme.ref_rev_cov = i32::MAX;
        extreme.msint = i32::MIN;
        extreme.bias = "2;1".to_string();
        extreme.var_allele = "<INV>".to_string();
        let _ = filter_variant(
            &extreme,
            &Config::default(),
            i32::MAX,
            Some((0, i32::MAX, i32::MIN, false)),
            0.0,
            "Inf",
        );
        let mut sv_zero_cfg = Config::default();
        sv_zero_cfg.vcf_min_splitreads = 0;
        let mut vcf_out = String::new();
        let ext_region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        format_filtered_vcf_record_into(
            &mut vcf_out,
            &extreme,
            &ext_region,
            &sv_zero_cfg,
            i32::MAX,
            0.0,
            "Inf",
            None,
        );
        assert!(vcf_out.contains(";SVTYPE=INV;SVLEN=2147483648;"));
        let san_rec = format_paired_vcf_record(
            &somatic_sv(extreme.clone(), Some(extreme), "Germline"),
            "chr1",
            &bad_lib_cfg,
        );
        assert!(san_rec.capacity() >= san_rec.as_str().len());
        assert!(!san_rec.as_str().contains('\n'));
        assert_eq!(san_rec.as_str().bytes().filter(|&b| b == b'\t').count(), 10);

        // `PairedRecordStreamer::with_pass_only(10, true)` suppresses clustered SNVs after
        // `Cluster10bp` mutation while keeping unclustered PASS records (`var2vcf_paired.pl:229-236, 249`)
        let mk_pos = |start: i32, qual: f64| {
            let mut v = clean_normal();
            v.start_pos = start;
            v.end_pos = start;
            v.mean_quality = qual;
            v.var_type = "SNV".to_string();
            format_paired_vcf_record(&somatic_sv(v, None, "StrongSomatic"), "chr1", &cfg)
        };
        let pass_only_lines = finish_paired_records_with_pass_only(
            vec![
                mk_pos(100, 35.0), // clustered with 106 -> Cluster10bp -> suppressed by pass_only
                mk_pos(106, 35.0), // clustered with 100 -> Cluster10bp -> suppressed by pass_only
                mk_pos(150, 10.0), // q22.5 -> suppressed by pass_only
                mk_pos(200, 35.0), // unclustered PASS -> emitted!
            ],
            10,
            true,
        );
        assert_eq!(pass_only_lines.len(), 1);
        assert!(pass_only_lines[0].starts_with("chr1\t200\t."));
        assert_eq!(pass_only_lines[0].split('\t').nth(6), Some("PASS"));
    }

    #[test]
    fn single_and_paired_vcf_headers_respect_all_11_overrides_and_sanitize_raw() {
        use crate::config::VcfThreshold;

        let lengths = HashMap::new();
        let render_single = |cfg: &Config| -> String {
            let mut b = Vec::new();
            print_vcf_header(&mut b, cfg, &lengths).unwrap();
            String::from_utf8(b).unwrap()
        };
        let render_paired = |cfg: &Config| -> String {
            let mut b = Vec::new();
            print_paired_vcf_header(&mut b, cfg, &lengths).unwrap();
            String::from_utf8(b).unwrap()
        };

        // 1. Default (`None`) headers preserve exact upstream text & effective_vcf_cluster_dist = 0
        let def = Config::default();
        assert_eq!(effective_vcf_cluster_dist(&def), 0);
        let sh0 = render_single(&def);
        let ph0 = render_paired(&def);
        assert!(sh0.contains(r#"##FILTER=<ID=Q10,Description="Mean Mapping Quality Below 10">"#));
        assert!(ph0.contains(r#"##FILTER=<ID=Q0,Description="Mean Mapping Quality Below 0">"#));
        assert!(
            sh0.contains(r#"##FILTER=<ID=Cluster0bp,Description="Two variants are within 0 bp">"#)
        );
        assert!(ph0.contains(
            r#"##FILTER=<ID=Cluster0bp,Description="Two somatic variants are within 0 bp">"#
        ));

        // 2. Fractional, zero, raw-token, and programmatic overrides across single & paired headers
        let mut custom = Config::default();
        custom.vcf_qmean = Some(VcfThreshold::new(20.0)); // programmatic -> "20"
        custom.vcf_mapq = Some(VcfThreshold::with_raw(0.0, "0".to_string()));
        custom.vcf_pmean = Some(VcfThreshold::with_raw(7.5, "7.5".to_string()));
        custom.vcf_sn = Some(VcfThreshold::with_raw(2.25, "2.25".to_string()));
        custom.vcf_min_depth = Some(VcfThreshold::new(0));
        custom.vcf_min_vd = Some(VcfThreshold::new(0));
        custom.vcf_af = Some(VcfThreshold::with_raw(0.05, "0.050".to_string()));
        custom.vcf_msi = Some(VcfThreshold::with_raw(10.5, "10.5".to_string()));
        custom.vcf_max_nm = Some(VcfThreshold::with_raw(4.5, "4.5".to_string()));
        custom.vcf_cluster_dist = Some(VcfThreshold::new(15));
        assert_eq!(effective_vcf_cluster_dist(&custom), 15);

        let sh = render_single(&custom);
        let ph = render_paired(&custom);
        for h in [&sh, &ph] {
            assert!(h.contains(r#"##FILTER=<ID=q20,Description="Mean Base Quality Below 20">"#));
            assert!(h.contains(r#"##FILTER=<ID=Q0,Description="Mean Mapping Quality Below 0">"#));
            assert!(h.contains(
                r#"##FILTER=<ID=p7.5,Description="Mean Position in Reads Less than 7.5">"#
            ));
            assert!(
                h.contains(r#"##FILTER=<ID=SN2.25,Description="Signal to Noise Less than 2.25">"#)
            );
            assert!(h.contains(r#"##FILTER=<ID=d0,Description="Total Depth < 0">"#));
            assert!(h.contains(r#"##FILTER=<ID=v0,Description="Var Depth < 0">"#));
            assert!(h.contains(r#"##FILTER=<ID=f0.050,Description="Allele frequency < 0.050">"#));
            assert!(h.contains(
                r#"##FILTER=<ID=NM4.5,Description="Mean mismatches in reads >= 4.5, thus likely false positive">"#
            ));
        }
        // Monomer MSI number stays fixed at 13 in single vs 12 in paired; Cluster text differs
        assert!(sh.contains(
            r#"##FILTER=<ID=MSI10.5,Description="Variant in MSI region with 10.5 non-monomer MSI or 13 monomer MSI">"#
        ));
        assert!(ph.contains(
            r#"##FILTER=<ID=MSI10.5,Description="Variant in MSI region with 10.5 non-monomer MSI or 12 monomer MSI">"#
        ));
        assert!(
            sh.contains(r#"##FILTER=<ID=Cluster15bp,Description="Two variants are within 15 bp">"#)
        );
        assert!(ph.contains(
            r#"##FILTER=<ID=Cluster15bp,Description="Two somatic variants are within 15 bp">"#
        ));

        // 2b. Leading-zero / `+`-prefixed integer tokens (`03`, `+02`, `010`) preserved verbatim
        let mut raw_ints = Config::default();
        raw_ints.vcf_min_depth = Some(VcfThreshold::with_raw(3, "03".to_string()));
        raw_ints.vcf_min_vd = Some(VcfThreshold::with_raw(2, "+02".to_string()));
        raw_ints.vcf_cluster_dist = Some(VcfThreshold::with_raw(10, "010".to_string()));
        let sh_raw = render_single(&raw_ints);
        let ph_raw = render_paired(&raw_ints);
        for h in [&sh_raw, &ph_raw] {
            assert!(h.contains(r#"##FILTER=<ID=d03,Description="Total Depth < 03">"#));
            assert!(h.contains(r#"##FILTER=<ID=v+02,Description="Var Depth < +02">"#));
        }
        assert!(sh_raw.contains(
            r#"##FILTER=<ID=Cluster010bp,Description="Two variants are within 010 bp">"#
        ));
        assert!(ph_raw.contains(
            r#"##FILTER=<ID=Cluster010bp,Description="Two somatic variants are within 010 bp">"#
        ));

        // Explicit `vcf_af = Some(0)` uses `0` even in single-sample mode (`defined()` semantics)
        let mut zero_af = Config::default();
        zero_af.vcf_af = Some(VcfThreshold::with_raw(0.0, "0".to_string()));
        assert!(render_single(&zero_af)
            .contains(r#"##FILTER=<ID=f0,Description="Allele frequency < 0">"#));

        // 3. Direct library `&Config` injection / overlong / mismatched `raw` sanitization
        let mut evil = Config::default();
        evil.vcf_qmean = Some(VcfThreshold::with_raw(
            20.0,
            "20.0\n##evil\tcol".to_string(),
        ));
        evil.vcf_mapq = Some(VcfThreshold::with_raw(12.0, "9".repeat(80)));
        evil.vcf_pmean = Some(VcfThreshold::with_raw(7.5, "99.0".to_string())); // bit mismatch
        evil.vcf_sn = Some(VcfThreshold::new(f64::NAN));
        evil.vcf_min_depth = Some(VcfThreshold::with_raw(3, "03\n##evil".to_string()));
        evil.vcf_min_vd = Some(VcfThreshold::with_raw(2, "99".to_string())); // value mismatch
        evil.vcf_cluster_dist = Some(VcfThreshold::new(-10));
        assert_eq!(effective_vcf_cluster_dist(&evil), 0);

        for h in [render_single(&evil), render_paired(&evil)] {
            assert!(!h.contains("##evil"));
            assert!(!h.contains(&"9".repeat(80)));
            assert!(h.contains("##FILTER=<ID=q20,"));
            assert!(h.contains("##FILTER=<ID=Q12,"));
            assert!(h.contains("##FILTER=<ID=p7.5,"));
            assert!(h.contains("##FILTER=<ID=SN1.5,"));
            assert!(h.contains("##FILTER=<ID=d3,"));
            assert!(h.contains("##FILTER=<ID=v2,"));
            assert!(h.contains("##FILTER=<ID=Cluster0bp,"));
        }
    }

    #[test]
    fn single_and_paired_record_filters_and_gt_respect_all_11_overrides() {
        use crate::config::VcfThreshold;

        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(300),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        // ---- Single-sample record filters & GT (`var2vcf_valid.pl`) ----
        let mut scfg = Config::default();
        scfg.vcf_min_depth = Some(VcfThreshold::new(8));
        scfg.vcf_min_vd = Some(VcfThreshold::new(5));
        scfg.vcf_af = Some(VcfThreshold::with_raw(0.05, "0.050".to_string()));
        scfg.vcf_pmean = Some(VcfThreshold::with_raw(7.5, "7.5".to_string()));
        scfg.vcf_qmean = Some(VcfThreshold::with_raw(25.0, "25.0".to_string()));
        scfg.vcf_mapq = Some(VcfThreshold::with_raw(20.0, "20".to_string()));
        scfg.vcf_sn = Some(VcfThreshold::with_raw(2.5, "2.5".to_string()));
        scfg.vcf_max_nm = Some(VcfThreshold::with_raw(4.0, "4".to_string()));
        scfg.vcf_msi = Some(VcfThreshold::with_raw(10.5, "10.5".to_string()));
        scfg.vcf_gt_freq = Some(0.05);
        scfg.vcf_cluster_dist = Some(VcfThreshold::new(15));

        let mut sv = amplicon_passing_variant();
        sv.start_pos = 10; // < opt_c (15) against initial pvs = 0!
        sv.end_pos = 10;
        // `amplicon_passing_variant()` has hq_lq_ratio = 2.0 < vcf_sn (2.5) -> fires SN2.5 and suppresses Cluster15bp
        assert_eq!(
            filter_variant(&sv, &scfg, 0, None, 1.0, "1.0"),
            vec!["SN2.5"]
        );
        sv.hq_lq_ratio = 3.0;
        // Initial pvs = 0 with start = 10 < 15 fires Cluster15bp on clean SNV (`var2vcf_valid.pl:132,181`)
        assert_eq!(
            filter_variant(&sv, &scfg, 0, None, 1.0, "1.0"),
            vec!["Cluster15bp"]
        );
        let mut try_buf = String::new();
        assert!(try_format_filtered_vcf_record_into(
            &mut try_buf,
            &sv,
            &region,
            &scfg,
            0,
            1.0,
            "1.0",
            None
        ));
        assert_eq!(try_buf.split('\t').nth(6), Some("Cluster15bp"));

        // NM strict `>` in single: nm == 4.0 passes, nm == 4.1 fires NM4
        sv.start_pos = 100;
        sv.end_pos = 100;
        sv.nm = 4.0;
        assert!(filter_variant(&sv, &scfg, 0, None, 1.0, "1.0").is_empty());
        sv.nm = 4.1;
        assert_eq!(filter_variant(&sv, &scfg, 0, None, 1.0, "1.0"), vec!["NM4"]);

        // Depth/VD rescue (`hicnt * hiaf >= 0.5`) vs failure under overridden d8 / v5
        let mut low_dv = amplicon_passing_variant();
        low_dv.hq_lq_ratio = 3.0;
        low_dv.total_pos_cov = 7;
        low_dv.hicnt = 4;
        low_dv.hi_freq = 0.1; // 4 * 0.1 = 0.4 < 0.5 -> d8;v5
        assert_eq!(
            filter_variant(&low_dv, &scfg, 0, None, 1.0, "1.0"),
            vec!["d8", "v5"]
        );
        low_dv.hi_freq = 0.2; // 4 * 0.2 = 0.8 >= 0.5 -> rescued!
        assert!(filter_variant(&low_dv, &scfg, 0, None, 1.0, "1.0").is_empty());

        // Raw integer spellings matching numeric defaults (`03`, `+02`, `010`) emit verbatim filter tags!
        let mut raw_icfg = Config::default();
        raw_icfg.vcf_min_depth = Some(VcfThreshold::with_raw(3, "03".to_string()));
        raw_icfg.vcf_min_vd = Some(VcfThreshold::with_raw(2, "+02".to_string()));
        raw_icfg.vcf_cluster_dist = Some(VcfThreshold::with_raw(10, "010".to_string()));
        let mut raw_low = amplicon_passing_variant();
        raw_low.total_pos_cov = 2;
        raw_low.hicnt = 1;
        raw_low.hi_freq = 0.1; // 1 * 0.1 = 0.1 < 0.5 -> d03;v+02
        assert_eq!(
            filter_variant(&raw_low, &raw_icfg, 0, None, 1.0, "1.0"),
            vec!["d03", "v+02"]
        );
        let mut raw_cl = amplicon_passing_variant();
        raw_cl.start_pos = 5;
        raw_cl.end_pos = 5;
        assert_eq!(
            filter_variant(&raw_cl, &raw_icfg, 0, None, 1.0, "1.0"),
            vec!["Cluster010bp"]
        );

        // `--vcf-mapq 0` disables Q even when mapq < 10 and af < 0.8
        let mut mapq0_cfg = Config::default();
        mapq0_cfg.vcf_mapq = Some(VcfThreshold::with_raw(0.0, "0".to_string()));
        let mut low_mq = amplicon_passing_variant();
        low_mq.hq_lq_ratio = 3.0;
        low_mq.mean_mapping_quality = 5.0;
        low_mq.frequency = 0.2;
        assert!(filter_variant(&low_mq, &mapq0_cfg, 0, None, 1.0, "1.0").is_empty());
        assert_eq!(
            filter_variant(&low_mq, &scfg, 0, None, 1.0, "1.0"),
            vec!["Q20"]
        );
        low_mq.frequency = 0.85; // af >= 0.8 suppresses Q20 in single mode
        assert!(filter_variant(&low_mq, &scfg, 0, None, 1.0, "1.0").is_empty());

        // Non-monomer MSI uses `vcf_msi` (10.5), monomer stays fixed at `>= 13.0`
        let mut msi_var = amplicon_passing_variant();
        msi_var.hq_lq_ratio = 3.0;
        msi_var.ref_allele = "ATG".to_string();
        msi_var.var_allele = "A".to_string();
        msi_var.var_type = "Deletion".to_string();
        msi_var.msint = 2;
        msi_var.msi = 11.0; // > 10.5 (fires MSI10.5) and >= 8 with hiaf=0.5 (>0.2 so no LongMSI)
        msi_var.hi_freq = 0.5;
        msi_var.frequency = 0.15;
        assert_eq!(
            filter_variant(&msi_var, &scfg, 0, None, 1.0, "1.0"),
            vec!["MSI10.5"]
        );
        msi_var.ref_allele = "AT".to_string();
        msi_var.msint = 1; // monomer: msi=11.0 < 13.0 -> no MSI10.5
        assert!(filter_variant(&msi_var, &scfg, 0, None, 1.0, "1.0").is_empty());

        // GT uses `vcf_af` (0.05) and `vcf_gt_freq` (0.05): af=0.90 -> 1/0, af=0.96 -> 1/1, af=0.04 -> 0/0
        let mut gt_var = amplicon_passing_variant();
        gt_var.hq_lq_ratio = 3.0;
        gt_var.frequency = 0.90; // 1 - 0.90 = 0.10 >= 0.05 -> 1/0 (default 0.2 would give 1/1!)
        let rec_90 = format_vcf_record(&gt_var, &region, &scfg, &[], 1.0, "1.0", None);
        assert!(rec_90.split('\t').nth(9).unwrap().starts_with("1/0:"));
        gt_var.frequency = 0.96; // 1 - 0.96 = 0.04 < 0.05 -> 1/1
        let rec_96 = format_vcf_record(&gt_var, &region, &scfg, &[], 1.0, "1.0", None);
        assert!(rec_96.split('\t').nth(9).unwrap().starts_with("1/1:"));
        gt_var.frequency = 0.04; // < vcf_af (0.05) -> f0.050 and GT=0/0
        let f_04 = filter_variant(&gt_var, &scfg, 0, None, 1.0, "1.0");
        assert_eq!(f_04, vec!["f0.050"]);
        let rec_04 = format_vcf_record(&gt_var, &region, &scfg, &f_04, 1.0, "1.0", None);
        assert!(rec_04.split('\t').nth(9).unwrap().starts_with("0/0:"));

        // ---- Paired record filters & GT (`var2vcf_paired.pl`) ----
        let mut pcfg = scfg.clone();
        pcfg.vcf_mapq = Some(VcfThreshold::with_raw(15.0, "15".to_string()));

        // Paired `pSTD` uses `vd < eff_min_depth` (8), dual `Q15;Q15` on SNV with mapq < 10 & < 15,
        // and inclusive `NM >= 4` on both normal_clean and tumor!
        let mut ptumor = clean_normal();
        ptumor.hq_lq_ratio = 3.0;
        ptumor.var_type = "SNV".to_string();
        ptumor.total_pos_cov = 20;
        ptumor.pos_coverage = 6; // >= 5 (default) and >= v5, but < d8 with pstd=false -> pSTD!
        ptumor.pstd = false;
        ptumor.mean_mapping_quality = 5.0; // < 15 AND (< 10 && SNV) -> Q15;Q15
        ptumor.nm = 4.0; // == 4.0 -> fires NM4 in paired (inclusive >=)
        ptumor.frequency = 0.90; // 1 - 0.90 = 0.10 >= 0.05 -> 1/0

        let mut pnorm = clean_normal();
        pnorm.nm = 4.0; // >= 4.0 -> normal_clean = false!
        pnorm.frequency = 0.04; // < vcf_af (0.05) -> gtm = 0/0

        let prec =
            format_paired_vcf_record(&somatic_sv(ptumor, Some(pnorm), "Germline"), "chr1", &pcfg);
        assert_eq!(prec.filter_str(), "pSTD;Q15;Q15;NM4");
        let pcols: Vec<&str> = prec.as_str().split('\t').collect();
        assert!(pcols[9].starts_with("1/0:"), "tumor GT: {}", pcols[9]);
        assert!(pcols[10].starts_with("0/0:"), "normal GT: {}", pcols[10]);

        // Paired raw integer spellings matching defaults (`05`, `+03`, `010`) emit `d05`, `v+03`, and `Cluster010bp`!
        let mut raw_pcfg = Config::default();
        raw_pcfg.vcf_min_depth = Some(VcfThreshold::with_raw(5, "05".to_string()));
        raw_pcfg.vcf_min_vd = Some(VcfThreshold::with_raw(3, "+03".to_string()));
        raw_pcfg.vcf_cluster_dist = Some(VcfThreshold::with_raw(10, "010".to_string()));
        let mut p_low_t = clean_normal();
        p_low_t.total_pos_cov = 4; // < 5 -> d05
        p_low_t.pos_coverage = 2; // < 3 -> v+03
        p_low_t.pstd = true;
        let mut p_low_n = clean_normal();
        p_low_n.pos_coverage = 0; // < 3 -> normal_clean = false
        let prec_raw = format_paired_vcf_record(
            &somatic_sv(p_low_t, Some(p_low_n.clone()), "LikelySomatic"),
            "chr1",
            &raw_pcfg,
        );
        assert_eq!(prec_raw.filter_str(), "d05;v+03");

        let mut ps1 = clean_normal();
        ps1.start_pos = 100;
        ps1.end_pos = 100;
        ps1.var_type = "SNV".to_string();
        let mut ps2 = ps1.clone();
        ps2.start_pos = 105;
        ps2.end_pos = 105;
        let mut pstreamer = PairedRecordStreamer::from_config(&raw_pcfg);
        let mut pout = Vec::new();
        pstreamer
            .push(
                format_paired_vcf_record(
                    &somatic_sv(ps1, Some(p_low_n.clone()), "StrongSomatic"),
                    "chr1",
                    &raw_pcfg,
                ),
                |l| -> Result<(), std::convert::Infallible> {
                    pout.push(l.to_string());
                    Ok(())
                },
            )
            .unwrap();
        pstreamer
            .push(
                format_paired_vcf_record(
                    &somatic_sv(ps2, Some(p_low_n), "StrongSomatic"),
                    "chr1",
                    &raw_pcfg,
                ),
                |l| -> Result<(), std::convert::Infallible> {
                    pout.push(l.to_string());
                    Ok(())
                },
            )
            .unwrap();
        pstreamer
            .finish(|l| -> Result<(), std::convert::Infallible> {
                pout.push(l.to_string());
                Ok(())
            })
            .unwrap();
        assert_eq!(pout.len(), 2);
        assert_eq!(pout[0].split('\t').nth(6), Some("Cluster010bp"));
        assert_eq!(pout[1].split('\t').nth(6), Some("Cluster010bp"));
    }

    #[test]
    fn single_vcf_record_span_checks_and_lazy_event_builder() {
        #[cfg(target_pointer_width = "64")]
        {
            assert_eq!(std::mem::size_of::<SingleVcfRecord>(), 48);
        }

        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(300),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let mut cfg = Config::default();
        let pass_snv = amplicon_passing_variant();

        // 1. Tight capacity on `build_single_vcf_event` even when worker scratch has 4096 B capacity
        let mut scratch = String::with_capacity(4096);
        let ev = build_single_vcf_event(&mut scratch, &pass_snv, &region, &cfg, None);
        let SingleVcfEvent::Record(mut rec) = ev else {
            panic!("expected Record, got {ev:?}");
        };
        assert!(scratch.capacity() >= 4096);
        assert_eq!(rec.capacity(), rec.as_str().len());
        assert!(rec.capacity() < 512);
        assert!(rec.is_pass());
        assert!(rec.is_snv());
        assert!(rec.pre_cluster_empty());

        // 2. `SingleVcfRecord::from_parts_checked` (9 tabs = 10 columns)
        let (fs, fe) = rec.filter_offsets();
        let line = rec.as_str().to_string();
        assert!(
            SingleVcfRecord::from_parts_checked(line.clone(), fs, fe, 100, true, true, false)
                .is_some()
        );
        assert!(
            SingleVcfRecord::from_parts_checked(
                format!("{line}\textra_sample"),
                fs,
                fe,
                100,
                true,
                true,
                false
            )
            .is_none(),
            "11-column paired line must be rejected by SingleVcfRecord"
        );
        assert!(
            SingleVcfRecord::from_parts_checked(line.clone(), fe, fs, 100, true, true, false)
                .is_none()
        );
        assert!(SingleVcfRecord::from_parts_checked(
            format!("{line}\n"),
            fs,
            fe,
            100,
            true,
            true,
            false
        )
        .is_none());

        rec.set_cluster_filter("Cluster10bp;AMPBIAS");
        assert!(!rec.is_pass());
        assert_eq!(rec.filter_str(), "Cluster10bp;AMPBIAS");
        assert_eq!(rec.as_str().bytes().filter(|&b| b == b'\t').count(), 9);

        // 3. `Skip`, `AmpBiasSnvClusterProbe`, and `SvSplitreadFlush` under `--pass-only`
        let mut empty_ref = pass_snv.clone();
        empty_ref.ref_allele.clear();
        assert_eq!(
            build_single_vcf_event(&mut scratch, &empty_ref, &region, &cfg, None),
            SingleVcfEvent::Skip
        );

        cfg.vcf_pass_only = true;
        cfg.vcf_cluster_dist = Some(crate::config::VcfThreshold::new(10));

        let mut low_q = pass_snv.clone();
        low_q.mean_quality = 10.0;
        assert_eq!(
            build_single_vcf_event(&mut scratch, &low_q, &region, &cfg, None),
            SingleVcfEvent::Skip
        );
        assert!(scratch.is_empty());

        // Amplicon SNV with clean pre-cluster filters + AMPBIAS emits probe when opt_c > 0
        assert_eq!(
            build_single_vcf_event(
                &mut scratch,
                &pass_snv,
                &region,
                &cfg,
                Some((1, 4, 0, true))
            ),
            SingleVcfEvent::AmpBiasSnvClusterProbe { start: 100 }
        );
        // When opt_c == 0, AMPBIAS SNV cannot cluster, so worker skips directly
        cfg.vcf_cluster_dist = Some(crate::config::VcfThreshold::new(0));
        assert_eq!(
            build_single_vcf_event(
                &mut scratch,
                &pass_snv,
                &region,
                &cfg,
                Some((1, 4, 0, true))
            ),
            SingleVcfEvent::Skip
        );

        // Dropped `<DEL>` SV (`splitreads = 0 < 1`): flushes iff clean under `vcf_pass_only`
        let mut dropped_sv = pass_snv.clone();
        dropped_sv.var_allele = "<DEL>".to_string();
        dropped_sv.var_type = "Deletion".to_string();
        dropped_sv.sv = "0-4-1".to_string();
        assert_eq!(
            build_single_vcf_event(&mut scratch, &dropped_sv, &region, &cfg, None),
            SingleVcfEvent::SvSplitreadFlush {
                cluster_probe_start: None
            }
        );
        dropped_sv.mean_quality = 10.0; // fails q22.5 -> hits `next` at :187 before :192 flush!
        assert_eq!(
            build_single_vcf_event(&mut scratch, &dropped_sv, &region, &cfg, None),
            SingleVcfEvent::Skip
        );
        cfg.vcf_pass_only = false; // without -S, even filtered <SV> flushes pending at :192-204!
        assert_eq!(
            build_single_vcf_event(&mut scratch, &dropped_sv, &region, &cfg, None),
            SingleVcfEvent::SvSplitreadFlush {
                cluster_probe_start: None
            }
        );
    }

    type TestAmpCols = (i32, i32, i32, bool);

    /// Reference line-by-line simulator of `var2vcf_valid.pl:132-244` for differential testing
    /// against `build_single_vcf_event` + `SingleRecordStreamer`.
    fn simulate_perl_var2vcf_valid(
        rows: &[(Variant, Option<TestAmpCols>)],
        region: &Region,
        config: &Config,
    ) -> Vec<String> {
        let opt_c = effective_vcf_cluster_dist(config);
        let cluster_tag = format_cluster_tag_from_config(config);
        let mut pvs: i32 = 0;
        let mut pending: Option<String> = None;
        let mut out = Vec::new();

        let set_line_filter = |line: &mut String, new_filt: &str| {
            let cols: Vec<&str> = line.split('\t').collect();
            let mut rebuilt = String::new();
            for (idx, col) in cols.iter().enumerate() {
                if idx > 0 {
                    rebuilt.push('\t');
                }
                rebuilt.push_str(if idx == 6 { new_filt } else { col });
            }
            *line = rebuilt;
        };
        let line_is_pass = |line: &str| line.split('\t').nth(6) == Some("PASS");

        for (variant, amp) in rows {
            if perl_row_skipped(&variant.ref_allele) {
                continue;
            }
            let fisher = crate::stats::FisherExact::new_with_mode(
                variant.ref_fwd_cov,
                variant.ref_rev_cov,
                variant.fwd,
                variant.rev,
                config.fisher,
            );
            let sbf = fisher.get_p_value();
            let odd_buf = fisher.odds_ratio_buf();
            let thr = ValidThresholds::from_config(config);
            let mut tag_scratch = VcfFloatTagScratch::new();
            let mut filters = StackFilterList::new();
            let (pre_empty, _has_amp) = collect_valid_filters_into(
                &mut filters,
                &mut tag_scratch,
                variant,
                config,
                &thr,
                Some(pvs),
                *amp,
                sbf,
                odd_buf.as_str(),
            );
            let var_type = if variant.var_type.is_empty() {
                "REF"
            } else {
                variant.var_type.as_str()
            };
            if var_type == "SNV"
                && pre_empty
                && i64::from(variant.start_pos) - i64::from(pvs) < i64::from(opt_c)
            {
                if let Some(prev) = pending.as_mut() {
                    set_line_filter(prev, cluster_tag.as_str());
                }
            }
            if config.vcf_pass_only && !filters.is_empty() {
                continue;
            }
            if let Some(prev) = pending.take() {
                if !config.vcf_pass_only || line_is_pass(&prev) {
                    out.push(prev);
                }
            }
            if amp.is_none()
                && variant.var_allele.contains('<')
                && parse_sv_counts(&variant.sv).0 < config.vcf_min_splitreads
            {
                pending = None;
                continue;
            }
            let filter_vec = filters.to_vec_strings();
            let rec = format_vcf_record(
                variant,
                region,
                config,
                &filter_vec,
                sbf,
                odd_buf.as_str(),
                *amp,
            );
            if var_type == "SNV" && filter_vec.is_empty() {
                pvs = variant.start_pos;
            }
            pending = Some(rec);
        }
        if let Some(prev) = pending.take() {
            if !config.vcf_pass_only || line_is_pass(&prev) {
                out.push(prev);
            }
        }
        out
    }

    #[test]
    fn single_record_streamer_matches_var2vcf_valid_perl_state_machine() {
        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(500),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mk_snv = |pos: i32, qual: f64| {
            let mut v = amplicon_passing_variant();
            v.start_pos = pos;
            v.end_pos = pos;
            v.mean_quality = qual;
            v
        };
        let mut empty_row = mk_snv(25, 35.0);
        empty_row.ref_allele.clear();

        let mut clean_dropped_sv = mk_snv(205, 35.0);
        clean_dropped_sv.var_allele = "<DEL>".to_string();
        clean_dropped_sv.var_type = "Deletion".to_string();
        clean_dropped_sv.end_pos = 260;
        clean_dropped_sv.sv = "0-5-1".to_string();

        let mut dirty_dropped_sv = clean_dropped_sv.clone();
        dirty_dropped_sv.start_pos = 305;
        dirty_dropped_sv.end_pos = 360;
        dirty_dropped_sv.mean_quality = 10.0; // fails q22.5

        let mut bias_amp_snv = mk_snv(408, 35.0);
        bias_amp_snv.hi_freq = 0.10;
        bias_amp_snv.bias = "2;1".to_string();
        bias_amp_snv.ref_fwd_cov = 30;
        bias_amp_snv.ref_rev_cov = 0;
        bias_amp_snv.fwd = 0;
        bias_amp_snv.rev = 12; // triggers Fisher Bias

        let mut snv_typed_dropped_sv = clean_dropped_sv.clone();
        snv_typed_dropped_sv.var_type = "SNV".to_string(); // synthetic library caller edge case
        snv_typed_dropped_sv.start_pos = 458;
        snv_typed_dropped_sv.end_pos = 495;

        let stream_rows: Vec<(Variant, Option<TestAmpCols>)> = vec![
            // 1. start=10 < opt_c(15) against initial pvs=0!
            (mk_snv(10, 35.0), None),
            // 2. Skip empty REF row
            (empty_row, None),
            // 3. PASS SNV at 30 followed by AMPBIAS SNV at 38 (AmpBiasSnvClusterProbe under -S)
            (mk_snv(30, 35.0), Some((4, 4, 0, false))),
            (mk_snv(38, 35.0), Some((1, 4, 0, true))),
            // 4. 3-SNV chain 100, 110 (clusters with 100, pvs stays 100), 120 (120-100>=15 -> PASS!)
            (mk_snv(100, 35.0), None),
            (mk_snv(110, 35.0), None),
            (mk_snv(120, 35.0), None),
            // 5. Clean dropped <DEL> at 205 flushes 200 and clears pending while keeping pvs=200
            (mk_snv(200, 35.0), None),
            (clean_dropped_sv, None),
            (mk_snv(210, 35.0), None),
            // 6. Dirty dropped <DEL> at 305: flushes 300 when !pass_only, skips at :187 when pass_only
            (mk_snv(300, 35.0), None),
            (dirty_dropped_sv, None),
            (mk_snv(310, 35.0), None),
            // 7. PASS SNV at 400 followed by Fisher-Bias + AMPBIAS SNV at 408 (pre_cluster NOT empty!)
            (mk_snv(400, 35.0), Some((4, 4, 0, false))),
            (bias_amp_snv, Some((1, 4, 0, true))),
            // 8. Synthetic `var_type == "SNV"` + `<DEL>` dropped SV at 458 clusters 450 before dropping!
            (mk_snv(450, 35.0), None),
            (snv_typed_dropped_sv, None),
            (mk_snv(480, 35.0), None),
        ];

        for cluster_opt in [
            None,
            Some(crate::config::VcfThreshold::new(0)),
            Some(crate::config::VcfThreshold::new(15)),
            Some(crate::config::VcfThreshold::with_raw(15, "015".to_string())),
            Some(crate::config::VcfThreshold::with_raw(15, "+15".to_string())),
        ] {
            for pass_only in [false, true] {
                let mut cfg = Config::default();
                cfg.vcf_cluster_dist = cluster_opt.clone();
                cfg.vcf_pass_only = pass_only;

                let expected = simulate_perl_var2vcf_valid(&stream_rows, &region, &cfg);

                let mut scratch = String::new();
                let mut streamer = SingleRecordStreamer::from_config(&cfg);
                let mut actual = Vec::new();
                for (var, amp) in &stream_rows {
                    let ev = build_single_vcf_event(&mut scratch, var, &region, &cfg, *amp);
                    streamer
                        .push_event(ev, |line| {
                            actual.push(line.to_string());
                            Ok(())
                        })
                        .unwrap();
                }
                streamer
                    .finish(|line| {
                        actual.push(line.to_string());
                        Ok(())
                    })
                    .unwrap();

                assert_eq!(
                    actual, expected,
                    "mismatch vs var2vcf_valid.pl simulator for cluster_opt={cluster_opt:?}, pass_only={pass_only}"
                );
                if !pass_only {
                    if let Some(ref thr) = cluster_opt {
                        if let Some(raw) = thr.raw.as_deref() {
                            let expected_tag = format!("Cluster{raw}bp");
                            assert!(
                                actual.iter().any(|l| l.contains(&expected_tag)),
                                "expected {expected_tag} in streamed output"
                            );
                        }
                    }
                }
            }
        }

        // Direct library `&Variant` callers with unrounded float fields evaluate filters
        // at `.var` table precision (`round_half_even`), matching formatted VCF fields.
        let mut unrounded_var = mk_snv(500, 22.46); // rounds to QUAL=22.5 >= 22.5
        unrounded_var.frequency = 0.01996; // rounds to AF=0.0200 >= 0.02
        unrounded_var.mean_position = 7.96; // rounds to PMEAN=8.0 >= 8.0
        unrounded_var.mean_mapping_quality = 9.96; // rounds to MQ=10.0 >= 10.0
        unrounded_var.nm = 5.24; // rounds to NM=5.2 <= 5.25
        let mut strict_cfg = Config::default();
        strict_cfg.af_threshold = 0.02;
        assert!(
            filter_variant(&unrounded_var, &strict_cfg, 0, None, 1.0, "0").is_empty(),
            "table-rounded boundary values must PASS in filter_variant"
        );
        strict_cfg.vcf_pass_only = true;
        let mut ev_scratch = String::new();
        let SingleVcfEvent::Record(rec) =
            build_single_vcf_event(&mut ev_scratch, &unrounded_var, &region, &strict_cfg, None)
        else {
            panic!("expected SingleVcfEvent::Record on table-rounded PASS boundary");
        };
        assert!(
            rec.as_str().contains("\t0/1:"),
            "expected table-rounded AF=0.0200 >= 0.02 to emit GT=0/1, got: {}",
            rec.as_str()
        );
    }
}
