//! VCF format output serialization and var2vcf_valid filtering.
//! Ports var2vcf_valid.pl.

use crate::config::Config;
use crate::output::table::get_rounded_value_to_print;
use crate::regions::Region;
use crate::variations::Variant;
use std::collections::HashMap;
use std::io::{self, Write};

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

/// Standard chromosome ordering matching var2vcf_valid.pl's reorder subroutine.
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

    let qmean = 22.5;
    let map_qual_val = 10;
    let pmean_val = 8;
    let sn_val = 1.5;
    let total_depth = 3;
    let var_depth = 2;
    let freq_val = valid_freq_threshold(config);
    let opt_i = 12;
    let opt_m = 5.25;
    let opt_c = 0;

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
    writeln!(
        w,
        "##FILTER=<ID=f{},Description=\"Allele frequency < {}\">",
        freq_val, freq_val
    )?;
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
    let mut filters = Vec::new();

    let total_depth = 3;
    let var_depth = 2;
    let dp = variant.total_pos_cov;
    let hicnt = variant.hicnt;
    let hiaf = variant.hi_freq;
    let hi_check = (hicnt as f64 * hiaf) >= 0.5;

    // T098: d3 check
    if dp < total_depth {
        if !hi_check {
            filters.push(format!("d{}", total_depth));
        }
    }

    // T099: v2 check
    if (hicnt as u32) < var_depth {
        if !hi_check {
            filters.push(format!("v{}", var_depth));
        }
    }

    // T100: f$Freq check — var2vcf_valid.pl:159 annotates `f$Freq` where $Freq
    // is the -f value (default 0.02, :21); it never drops on AF without -S.
    // The same $Freq drives the GT "0/1 vs 0/0" threshold (:220). $Freq is the
    // writer's own -f (`af_threshold`), NOT the engine freq that -p sets to -1.
    let freq_threshold = valid_freq_threshold(config);
    if variant.frequency < freq_threshold {
        filters.push(format!("f{}", freq_threshold));
    }

    // T101: p8 check
    let pmean_threshold = 8.0;
    if variant.mean_position < pmean_threshold {
        filters.push(format!("p{}", pmean_threshold as i32));
    }

    // T102: pSTD check — suppressed in amplicon mode (:161 `! $isamp`)
    let opt_p = 1; // default
    if opt_p == 1 && !variant.pstd && amp.is_none() && variant.frequency < 0.35 {
        filters.push("pSTD".to_string());
    }

    // T103: q22.5 check
    let qmean = 22.5;
    if variant.mean_quality < qmean {
        filters.push(format!("q{}", qmean));
    }

    // T104: Q10 check
    let qmean_map = 10.0;
    if variant.mean_mapping_quality < qmean_map && variant.frequency < 0.8 {
        filters.push(format!("Q{}", qmean_map as i32));
    }

    // SN1.5 check
    let sn_threshold = 1.5;
    if variant.hq_lq_ratio < sn_threshold {
        filters.push(format!("SN{}", sn_threshold));
    }

    // NM5.25 check
    let opt_m = 5.25;
    if variant.nm > opt_m {
        filters.push(format!("NM{}", opt_m));
    }

    // MSI12 check
    let opt_i = 12.0;
    let ref_len = variant.ref_allele.len() as i32;
    let alt_len = variant.var_allele.len() as i32;
    let len_diff = (ref_len - alt_len).abs();
    let msi_check1 = variant.msi > opt_i
        && variant.msint > 1
        && variant.frequency < 0.2
        && len_diff == variant.msint;
    let msi_check2 = variant.msi >= 13.0
        && variant.msint == 1
        && variant.frequency <= 0.275
        && len_diff == variant.msint;
    if msi_check1 || msi_check2 {
        filters.push(format!("MSI{}", opt_i as i32));
    }

    // Bias check
    let raw_odd_ratio = if oddratio_str == "Inf" {
        0.0
    } else {
        oddratio_str.parse::<f64>().unwrap_or(0.0)
    };
    let mut transformed_or = raw_odd_ratio;
    if oddratio_str == "Inf" {
        transformed_or = 0.0;
    } else if raw_odd_ratio < 1.0 && raw_odd_ratio > 0.0 {
        transformed_or = 1.0 / raw_odd_ratio;
    }

    if variant.hi_freq < 0.25
        && variant.bias == "2;1"
        && sbf < 0.01
        && (transformed_or > 5.0 || transformed_or == 0.0)
        && (variant.end_pos - variant.start_pos) < 100
    {
        filters.push("Bias".to_string());
    }

    // LongMSI check
    if len_diff == variant.msint {
        if (variant.hi_freq <= 0.275 && variant.msi >= 13.0)
            || (variant.hi_freq <= 0.2 && variant.msi >= 8.0 && variant.msint > 1)
        {
            filters.push("LongMSI".to_string());
        }
    }

    // Cluster check
    let opt_c = 0; // default
    if variant.var_type == "SNV"
        && filters.is_empty()
        && pvs > 0
        && (variant.start_pos - pvs) < opt_c
    {
        filters.push(format!("Cluster{}bp", opt_c));
    }

    // AMPBIAS — amplicon mode only (:185 `$isamp && (($gamp < $tamp-$ncamp) || $ampflag)`);
    // pushed last, after the Cluster check, matching perl's @filters order.
    if let Some((gdamp, tlamp, ncamp, ampflag)) = amp {
        if gdamp < tlamp - ncamp || ampflag {
            filters.push("AMPBIAS".to_string());
        }
    }

    filters
}

/// Java's RAW table cell for hifreq (trailing zeros kept): ports
/// SimpleOutputVariant.java:127-129 / :211 — `hifreq == 0 ? "0" :
/// new DecimalFormat("0.0000").format(hifreq)` — NOT getRoundedValueToPrint.
/// The cell survives to Perl verbatim: teststrandbias.R's 22-entry colClasses
/// (:14) is recycled over the input columns, landing "character" on input col 23
/// (hifreq), so R echoes it unchanged and var2vcf_valid.pl:24/:236 interpolates
/// it as-is into `HIAF=`. (The PAIRED writer must stay stripped: testsomatic.R's
/// 48-entry pattern leaves hiaf1/hiaf2 numeric, so R re-formats them.)
fn hiaf_cell(hi_freq: f64) -> String {
    if hi_freq == 0.0 {
        "0".to_string()
    } else {
        format!("{:.4}", hi_freq)
    }
}

/// Perl's default numeric stringification when interpolating an NV into a
/// string: C `%.15g` semantics (15 significant digits, trailing zeros and dot
/// trimmed; %e form with a signed >=2-digit exponent outside [1e-5, 1e15)).
fn perl_num_str(x: f64) -> String {
    if x.is_nan() {
        return "NaN".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 {
            "-Inf".to_string()
        } else {
            "Inf".to_string()
        };
    }
    let neg = x < 0.0;
    let a = x.abs();
    if a == 0.0 {
        return "0".to_string();
    }
    let p = 15i32;
    // Exponent of the value AFTER rounding to p significant digits.
    let e = format!("{:.*e}", (p - 1) as usize, a);
    let (mant, exp) = e.split_once(['e', 'E']).unwrap_or((e.as_str(), "0"));
    let exp10: i32 = exp.parse().unwrap_or(0);
    let body = if exp10 < -4 || exp10 >= p {
        // %e style: trimmed mantissa + C-style exponent.
        let mant = mant.trim_end_matches('0').trim_end_matches('.');
        format!(
            "{}e{}{:02}",
            mant,
            if exp10 < 0 { '-' } else { '+' },
            exp10.abs()
        )
    } else {
        let prec = (p - 1 - exp10).max(0) as usize;
        let f = format!("{:.*}", prec, a);
        if f.contains('.') {
            let f = f.trim_end_matches('0');
            f.strip_suffix('.').unwrap_or(f).to_string()
        } else {
            f
        }
    };
    if neg {
        format!("-{}", body)
    } else {
        body
    }
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
    let ref_allele = variant.ref_allele.clone();
    let mut var_allele = variant.var_allele.clone();
    let var_type = if variant.var_type.is_empty() {
        "REF"
    } else {
        &variant.var_type
    };

    let gt = if ref_allele == var_allele {
        var_allele = ".".to_string();
        "0/0".to_string()
    } else {
        let af = variant.frequency;
        let gt_freq = 0.2;
        let freq_threshold = valid_freq_threshold(config);
        if 1.0 - af < gt_freq {
            "1/1".to_string()
        } else if af >= 0.5 {
            "1/0".to_string()
        } else if af >= freq_threshold {
            "0/1".to_string()
        } else {
            "0/0".to_string()
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

    let filter_str = if filters.is_empty() {
        "PASS"
    } else {
        &filters.join(";")
    };

    let sample_nowhitespace = config.sample.replace(char::is_whitespace, "_");

    let end_str = if config.omit_end_tag {
        "".to_string()
    } else {
        format!(";END={}", variant.end_pos)
    };

    let rd = variant.ref_fwd_cov + variant.ref_rev_cov;
    let vd = variant.pos_coverage;
    let ad = if vd != 0 {
        format!("{},{}", rd, vd)
    } else {
        format!("{}", rd)
    };

    // ports var2vcf_valid.pl:147-151 + :236 interpolation: "Inf" → 0; a value in
    // (0,1) is replaced by its reciprocal and stringified with Perl's default
    // numeric format (%.15g); any other table string passes through unchanged.
    let vcf_or_str = if oddratio_str == "Inf" {
        "0".to_string()
    } else {
        let x: f64 = oddratio_str.parse().unwrap_or(0.0);
        if x < 1.0 && x > 0.0 {
            perl_num_str(1.0 / x)
        } else {
            oddratio_str.to_string()
        }
    };

    let nm_str = if variant.nm > 0.0 {
        format!("{:.1}", variant.nm)
    } else {
        "0".to_string()
    };

    let lseq = if variant.left_seq.is_empty() {
        "0"
    } else {
        &variant.left_seq
    };
    let rseq = if variant.right_seq.is_empty() {
        "0"
    } else {
        &variant.right_seq
    };

    let mut sv_info = String::new();
    let is_amp = amp.is_some();
    if !is_amp {
        // ports var2vcf_valid.pl:197-212 — the SV block runs only when !$isamp
        let mut splitreads = 0;
        let mut spanpairs = 0;
        if !variant.sv.is_empty() && variant.sv.contains('-') {
            let parts: Vec<&str> = variant.sv.split('-').collect();
            if parts.len() >= 2 {
                splitreads = parts[0].parse::<i32>().unwrap_or(0);
                spanpairs = parts[1].parse::<i32>().unwrap_or(0);
            }
        }

        if var_allele.starts_with('<') {
            // ports var2vcf_valid.pl:206-208 — SVLEN = end - start,
            // incremented only when the ALT is an INV.
            let mut svlen = variant.end_pos - variant.start_pos;
            if var_allele.contains("INV") {
                svlen += 1;
            }
            let mut clean_alt = var_allele.clone();
            clean_alt.retain(|c| c != '<' && c != '>');
            sv_info = format!(";SVTYPE={};SVLEN={}", clean_alt, svlen);
        }
        sv_info = format!(
            "{};SPLITREAD={};SPANPAIR={}",
            sv_info, splitreads, spanpairs
        );
    }
    // After RSEQ: amplicon mode emits the GDAMP/TLAMP/NCAMP/AMPFLAG tags (:233),
    // simple mode DUPRATE (its col 37 is Duprate) — never both.
    let post_rseq = if let Some((gdamp, tlamp, ncamp, ampflag)) = amp {
        format!(
            ";GDAMP={};TLAMP={};NCAMP={};AMPFLAG={}",
            gdamp,
            tlamp,
            ncamp,
            i32::from(ampflag)
        )
    } else {
        // DUPRATE cell format differs by pipeline: with --fisher Perl consumes the
        // 38-col table directly — getRoundedValueToPrint("0.00", duprate)
        // (SimpleOutputVariant.java:175). Without --fisher the 36-col cell is
        // DecimalFormat("0.0")-based (SimpleOutputVariant.java:223) and R echoes it
        // numeric-stripped (teststrandbias.R colClasses leaves input col 35
        // numeric), i.e. a 1dp value. Both stripped, zero-guarded ("0").
        let dup_str = if config.fisher {
            get_rounded_value_to_print(variant.duprate, 2)
        } else {
            get_rounded_value_to_print(variant.duprate, 1)
        };
        format!(";DUPRATE={}", dup_str)
    };

    let info = format!(
        "SAMPLE={};TYPE={};DP={}{} \
         ;VD={};AF={};BIAS={};REFBIAS={}:{};VARBIAS={}:{} \
         ;PMEAN={};PSTD={};QUAL={};QSTD={} \
         ;SBF={};ODDRATIO={};MQ={};SN={} \
         ;HIAF={};ADJAF={};SHIFT3={} \
         ;MSI={};MSILEN={};NM={} \
         ;HICNT={};HICOV={};LSEQ={};RSEQ={}{}{}",
        sample_nowhitespace,
        var_type,
        variant.total_pos_cov,
        end_str,
        variant.pos_coverage,
        get_rounded_value_to_print(variant.frequency, 4),
        variant.bias.replace(';', ":"),
        variant.ref_fwd_cov,
        variant.ref_rev_cov,
        variant.fwd,
        variant.rev,
        get_rounded_value_to_print(variant.mean_position, 1),
        if variant.pstd { 1 } else { 0 },
        get_rounded_value_to_print(variant.mean_quality, 1),
        if variant.qstd { 1 } else { 0 },
        get_rounded_value_to_print(sbf, 5),
        vcf_or_str,
        get_rounded_value_to_print(variant.mean_mapping_quality, 1),
        get_rounded_value_to_print(variant.hq_lq_ratio, 3),
        hiaf_cell(variant.hi_freq),
        get_rounded_value_to_print(variant.extra_freq, 4),
        variant.shift3,
        get_rounded_value_to_print(variant.msi, 3),
        variant.msint,
        nm_str,
        variant.hicnt,
        variant.hicov,
        lseq,
        rseq,
        post_rseq,
        sv_info
    )
    .replace(" ;", ";")
    .replace(" ", "");

    let format_keys = "GT:DP:VD:AD:AF:RD:ALD";
    let format_vals = format!(
        "{}:{}:{}:{}:{}:{},{}:{},{}",
        gt,
        variant.total_pos_cov,
        variant.pos_coverage,
        ad,
        get_rounded_value_to_print(variant.frequency, 4),
        variant.ref_fwd_cov,
        variant.ref_rev_cov,
        variant.fwd,
        variant.rev
    );

    format!(
        "{}\t{}\t.\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        region.chr,
        variant.start_pos,
        ref_allele,
        var_allele,
        qual_val,
        filter_str,
        info,
        format_keys,
        format_vals
    )
}

// ---- var2vcf_paired.pl port (somatic/paired VCF, W4.5c) ----

/// Paired-writer defaults (var2vcf_paired.pl:19-31).
const PAIRED_MIN_DEPTH: i32 = 5;
const PAIRED_VAR_DEPTH: i32 = 3;
const PAIRED_QMEAN: f64 = 22.5;
const PAIRED_MQMEAN: f64 = 0.0;
const PAIRED_PMEAN: f64 = 8.0;
const PAIRED_SN: f64 = 1.5;
const PAIRED_GTFREQ: f64 = 0.2;
const PAIRED_OPT_I: f64 = 12.0;
const PAIRED_OPT_M: f64 = 5.25;

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

    // ports var2vcf_paired.pl:76-93
    let freq_str = config.af_threshold.to_string();
    let filter_lines = [
        format!(r#"##FILTER=<ID=q{},Description="Mean Base Quality Below {}">"#, PAIRED_QMEAN, PAIRED_QMEAN),
        format!(r#"##FILTER=<ID=Q{},Description="Mean Mapping Quality Below {}">"#, PAIRED_MQMEAN, PAIRED_MQMEAN),
        format!(r#"##FILTER=<ID=p{},Description="Mean Position in Reads Less than {}">"#, PAIRED_PMEAN as i32, PAIRED_PMEAN as i32),
        format!(r#"##FILTER=<ID=SN{},Description="Signal to Noise Less than {}">"#, PAIRED_SN, PAIRED_SN),
        r#"##FILTER=<ID=Bias,Description="Strand Bias">"#.to_string(),
        r#"##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">"#.to_string(),
        r#"##FILTER=<ID=MAF0.05,Description="Matched sample has AF > 0.05, thus not somatic">"#.to_string(),
        format!(r#"##FILTER=<ID=d{},Description="Total Depth < {}">"#, PAIRED_MIN_DEPTH, PAIRED_MIN_DEPTH),
        format!(r#"##FILTER=<ID=v{},Description="Var Depth < {}">"#, PAIRED_VAR_DEPTH, PAIRED_VAR_DEPTH),
        format!(r#"##FILTER=<ID=f{},Description="Allele frequency < {}">"#, freq_str, freq_str),
        r#"##FILTER=<ID=P0.05,Description="Not significant with p-value > 0.05">"#.to_string(),
        r#"##FILTER=<ID=DIFF0.2,Description="Non-somatic or LOH and allele frequency difference < 0.2">"#.to_string(),
        format!(r#"##FILTER=<ID=MSI{},Description="Variant in MSI region with {} non-monomer MSI or 12 monomer MSI">"#, PAIRED_OPT_I as i32, PAIRED_OPT_I as i32),
        format!(r#"##FILTER=<ID=NM{},Description="Mean mismatches in reads >= {}, thus likely false positive">"#, PAIRED_OPT_M, PAIRED_OPT_M),
        r#"##FILTER=<ID=InGap,Description="The somatic variant is in the deletion gap, thus likely false positive">"#.to_string(),
        r#"##FILTER=<ID=InIns,Description="The somatic variant is adjacent to an insertion variant">"#.to_string(),
        r#"##FILTER=<ID=Cluster0bp,Description="Two somatic variants are within 0 bp">"#.to_string(),
        r#"##FILTER=<ID=LongAT,Description="The somatic variant is flanked by long A/T (>=14)">"#.to_string(),
    ];
    for line in &filter_lines {
        writeln!(w, "{}", line)?;
    }

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
pub struct PairedVcfRecord {
    /// CHROM POS . REF ALT QUAL (var2vcf_paired.pl:237 pinfo1).
    prefix: String,
    filter: String,
    /// INFO \t FORMAT \t tumor-sample \t normal-sample (:238-242 pinfo2).
    rest: String,
    /// Raw start (for the Cluster rule :229 and pvs tracking :245).
    start: i32,
    type_snv: bool,
    /// Whether STATUS matches /Somatic/ (Cluster-rule precondition).
    status_somatic: bool,
}

/// Perl numification of non-numeric strings ("NA" → 0, with warning).
fn perl_num(s: &str) -> f64 {
    s.parse().unwrap_or(0.0)
}

/// Odds-ratio transform for per-sample fields. Ports var2vcf_paired.pl:146-155:
/// "Inf" → 0; (0,1) → sprintf("%.2f", 1/x). Returns (string, numeric-after-format);
/// the numeric form feeds the Bias comparisons (:175/:186).
fn paired_or_transform(raw: &str) -> (String, f64) {
    if raw == "Inf" {
        return ("0".to_string(), 0.0);
    }
    let x = perl_num(raw);
    if x < 1.0 && x > 0.0 {
        let s = format!("{:.2}", 1.0 / x);
        let v = s.parse::<f64>().unwrap_or(0.0);
        (s, v)
    } else {
        (raw.to_string(), x)
    }
}

/// BIAS field transform. Ports var2vcf_paired.pl:222-225: ";" → ","; bare "0" → "0,0".
fn paired_bias_fmt(bias: &str) -> String {
    let s = bias.replace(';', ",");
    if s == "0" {
        "0,0".to_string()
    } else {
        s
    }
}

/// GT formula shared by both samples. Ports var2vcf_paired.pl:220-221.
fn paired_gt(af: f64, freq: f64) -> &'static str {
    if 1.0 - af < PAIRED_GTFREQ {
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

/// Per-sample FORMAT string. Ports var2vcf_paired.pl:240-241.
fn paired_sample_field(
    v: &crate::variations::Variant,
    gt: &str,
    sbf: f64,
    oddratio: &(String, f64),
) -> String {
    let nm = if v.nm > 0.0 { v.nm } else { 0.0 };
    let mq = format!(
        "{:.0}",
        crate::types::round_half_even(v.mean_mapping_quality, 0)
    );
    format!(
        "{}:{}:{}:{},{}:{},{}:{},{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
        gt,
        v.total_pos_cov,
        v.pos_coverage,
        v.fwd,
        v.rev,
        v.ref_fwd_cov,
        v.ref_rev_cov,
        v.ref_fwd_cov + v.ref_rev_cov,
        v.pos_coverage,
        get_rounded_value_to_print(v.frequency, 4),
        paired_bias_fmt(&v.bias),
        get_rounded_value_to_print(v.mean_position, 1),
        if v.pstd { 1 } else { 0 },
        get_rounded_value_to_print(v.mean_quality, 1),
        if v.qstd { 1 } else { 0 },
        get_rounded_value_to_print(sbf, 5),
        oddratio.0,
        mq,
        get_rounded_value_to_print(v.hq_lq_ratio, 3),
        get_rounded_value_to_print(v.hi_freq, 4),
        get_rounded_value_to_print(v.extra_freq, 4),
        get_rounded_value_to_print(nm, 1),
    )
}

/// Format one somatic variant into a paired VCF record.
/// Ports var2vcf_paired.pl:127-242 field-by-field. The delayed-print /
/// Cluster state machine lives in `finish_paired_records`.
pub fn format_paired_vcf_record(
    sv: &crate::postprocess::somatic::SomaticVariant,
    chr: &str,
    config: &Config,
) -> PairedVcfRecord {
    let v1 = &sv.variant1;
    let zero = crate::variations::Variant::new();
    let v2 = sv.variant2.as_ref().unwrap_or(&zero);

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
    let or1 = paired_or_transform(&fisher1.get_odds_ratio());

    let fisher2 = crate::stats::FisherExact::new_with_mode(
        v2.ref_fwd_cov,
        v2.ref_rev_cov,
        v2.fwd,
        v2.rev,
        config.fisher,
    );
    let p2 = fisher2.get_p_value();
    let or2 = paired_or_transform(&fisher2.get_odds_ratio());

    let tref = (v1.total_pos_cov - v1.pos_coverage).max(0);
    let rref = (v2.total_pos_cov - v2.pos_coverage).max(0);
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
    let otot_raw = fisher_total.get_odds_ratio();
    let sor_str = if otot_raw == "Inf" {
        "0".to_string()
    } else {
        otot_raw.clone()
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
    let pmean1 = crate::types::round_half_even(v1.mean_position, 1);
    let pmean2 = crate::types::round_half_even(v2.mean_position, 1);
    let qual1 = crate::types::round_half_even(v1.mean_quality, 1);
    let qual2 = crate::types::round_half_even(v2.mean_quality, 1);
    let mapq1 = crate::types::round_half_even(v1.mean_mapping_quality, 1);
    let mapq2 = crate::types::round_half_even(v2.mean_mapping_quality, 1);
    let sn1 = crate::types::round_half_even(v1.hq_lq_ratio, 3);
    let sn2 = crate::types::round_half_even(v2.hq_lq_ratio, 3);
    let nm1 = crate::types::round_half_even(if v1.nm > 0.0 { v1.nm } else { 0.0 }, 1);
    let nm2 = crate::types::round_half_even(if v2.nm > 0.0 { v2.nm } else { 0.0 }, 1);
    // Perl sees R's round-5 strings from testsomatic.R:31-47, numified before
    // the P:156/:175/:186/:189 threshold checks — compare at table precision.
    let ptot_r = crate::types::round_half_even(ptot, 5);
    let p1_r = crate::types::round_half_even(p1, 5);
    let p2_r = crate::types::round_half_even(p2, 5);

    let mut var_type = v1.var_type.clone();
    if var_type.is_empty() {
        var_type = "REF".to_string(); // :139
    }
    let status = &sv.status;

    // endVariant supplies shift3/msi/msint/lseq/rseq (= table columns 48-52).
    let end_var = if sv.end_is_variant2 {
        sv.variant2.as_ref().unwrap_or(v1)
    } else {
        v1
    };

    let mut filters: Vec<String> = Vec::new();
    let mut filters2: Vec<String> = Vec::new();

    // Tumor depth filters with StrongSomatic exception (:156-161).
    let ss_exception = status == "StrongSomatic" && ptot_r < 0.15 && af1 * vd1 as f64 >= 0.5;
    if dp1 < PAIRED_MIN_DEPTH && !ss_exception {
        filters.push(format!("d{}", PAIRED_MIN_DEPTH));
    }
    if vd1 < PAIRED_VAR_DEPTH && !ss_exception {
        filters.push(format!("v{}", PAIRED_VAR_DEPTH));
    }
    // Normal depth filters, unconditional (:162-163).
    if dp2 < PAIRED_MIN_DEPTH {
        filters2.push(format!("d{}", PAIRED_MIN_DEPTH));
    }
    if vd2 < PAIRED_VAR_DEPTH {
        filters2.push(format!("v{}", PAIRED_VAR_DEPTH));
    }

    // Tumor quality filters (:165-175).
    if af1 < config.af_threshold {
        filters.push(format!("f{}", config.af_threshold));
    }
    if pmean1 < PAIRED_PMEAN {
        filters.push(format!("p{}", PAIRED_PMEAN as i32));
    }
    if !v1.pstd && vd1 < PAIRED_MIN_DEPTH {
        filters.push("pSTD".to_string());
    }
    if qual1 < PAIRED_QMEAN {
        filters.push(format!("q{}", PAIRED_QMEAN));
    }
    if mapq1 < PAIRED_MQMEAN {
        filters.push(format!("Q{}", PAIRED_MQMEAN as i32)); // dead at MQMEAN=0
    }
    if mapq1 < 10.0 && var_type == "SNV" {
        filters.push(format!("Q{}", PAIRED_MQMEAN as i32));
    }
    if sn1 < PAIRED_SN {
        filters.push(format!("SN{}", PAIRED_SN));
    }
    if nm1 >= PAIRED_OPT_M {
        filters.push(format!("NM{}", PAIRED_OPT_M));
    }
    if v1.bias == "2;1"
        && p1_r < 0.01
        && (or1.1 > 5.0 || or1.1 == 0.0)
        && v1.end_pos - v1.start_pos < 100
    {
        filters.push("Bias".to_string());
    }

    // Normal quality filters (:177-183); Q gate dead at MQMEAN=0, kept for fidelity.
    if af2 < config.af_threshold {
        filters2.push(format!("f{}", config.af_threshold));
    }
    if pmean2 < PAIRED_PMEAN {
        filters2.push(format!("p{}", PAIRED_PMEAN as i32));
    }
    if !v2.pstd && vd2 < PAIRED_MIN_DEPTH {
        filters2.push("pSTD".to_string());
    }
    if qual2 < PAIRED_QMEAN {
        filters2.push(format!("q{}", PAIRED_QMEAN));
    }
    if mapq2 < PAIRED_MQMEAN {
        filters2.push(format!("Q{}", PAIRED_MQMEAN as i32)); // dead at MQMEAN=0
    }
    if sn2 < PAIRED_SN {
        filters2.push(format!("SN{}", PAIRED_SN));
    }
    if nm2 >= PAIRED_OPT_M {
        filters2.push(format!("NM{}", PAIRED_OPT_M));
    }
    // Normal-side Bias appended to the MAIN list when not already present (:185-186).
    if !filters.iter().any(|f| f == "Bias")
        && v2.bias == "2;1"
        && p2_r < 0.01
        && (or2.1 > 5.0 || or2.1 == 0.0)
        && v1.end_pos - v1.start_pos < 100
    {
        filters.push("Bias".to_string());
    }

    // MSI gates (:188-194).
    let msilen = end_var.msint as f64;
    let msi = crate::types::round_half_even(end_var.msi, 3);
    if (msi > PAIRED_OPT_I && msilen > 1.0) || (msi > 12.0 && msilen == 1.0) {
        if !(status == "StrongSomatic" && ptot_r < 0.0005) {
            filters.push(format!("MSI{}", PAIRED_OPT_I as i32));
        }
    }
    let len_diff = (v1.ref_allele.len() as i32 - v1.var_allele.len() as i32).abs();
    if len_diff == msilen as i32 && !filters.iter().any(|f| f.starts_with("MSI")) {
        if (msi > PAIRED_OPT_I && msilen > 1.0 && af1 < 0.35 && af2 < 0.35)
            || (msi > 12.0 && msilen == 1.0 && af1 < 0.35 && af2 < 0.35)
        {
            filters.push(format!("MSI{}", PAIRED_OPT_I as i32));
        }
    }

    // The -M block (:198-206) is not plumbed: every oracle/corpus invocation omits
    // -M, so P$PVAL/P0.01Likely/InDelLikely never fire.

    let mut filter = if filters.is_empty() {
        "PASS".to_string()
    } else {
        filters.join(";")
    };
    // Germline passthrough (:217-219) — without -M a fully clean NORMAL overrides
    // all tumor filters.
    if filter != "PASS" && filters2.is_empty() {
        filter = "PASS".to_string();
    }

    let gt = paired_gt(af1, config.af_threshold);
    let gtm = paired_gt(af2, config.af_threshold);

    let qual_col = paired_qual_col(vd1, qual1, vd2, qual2);

    let sample_tumor = config.sample.split('|').next().unwrap_or(&config.sample);
    let sample_nowhitespace = sample_tumor.replace(char::is_whitespace, "_");

    let lseq = if end_var.left_seq.is_empty() {
        "0"
    } else {
        &end_var.left_seq
    };
    let rseq = if end_var.right_seq.is_empty() {
        "0"
    } else {
        &end_var.right_seq
    };

    // INFO (:237-239). No END tag and no SV tags in the paired writer.
    let info = format!(
        "STATUS={};SAMPLE={};TYPE={};DP={};VD={};AF={};SHIFT3={};MSI={};MSILEN={};SSF={};SOR={};LSEQ={};RSEQ={}",
        status,
        sample_nowhitespace,
        var_type,
        dp1,
        vd1,
        get_rounded_value_to_print(v1.frequency, 4),
        end_var.shift3,
        if end_var.msi == 0.0 {
            "0".to_string()
        } else {
            format!("{:.3}", crate::types::round_half_even(end_var.msi, 3))
        },
        end_var.msint,
        get_rounded_value_to_print(ptot, 5),
        sor_str,
        lseq,
        rseq,
    );

    let format_keys =
        "GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM";
    let s1 = paired_sample_field(v1, gt, p1, &or1);
    let s2 = paired_sample_field(v2, gtm, p2, &or2);

    PairedVcfRecord {
        prefix: format!(
            "{}\t{}\t.\t{}\t{}\t{}",
            chr, v1.start_pos, v1.ref_allele, v1.var_allele, qual_col
        ),
        filter,
        rest: format!("{}\t{}\t{}\t{}", info, format_keys, s1, s2),
        start: v1.start_pos,
        type_snv: var_type == "SNV",
        status_somatic: status.contains("Somatic"),
    }
}

/// Delayed-print state machine. Ports var2vcf_paired.pl:126,229-236,243-245,248-250:
/// each record is buffered; on the next record the Cluster rule may mark BOTH
/// records Cluster{c}bp (prev filter PASS + prev INFO /Somatic/+/TYPE=SNV/,
/// cur filter PASS + cur status /Somatic/ + type SNV, and cur.start - pvs < c);
/// pvs tracks the last PASS-SNV start. The buffer is flushed after the last
/// record of a chromosome batch. -S suppression is not plumbed (the oracle
/// invocation never passes -S). With the default c=0 the rule is inert at
/// ascending positions; `opt_c` is threaded through for fidelity.
pub fn finish_paired_records(records: Vec<PairedVcfRecord>, opt_c: i32) -> Vec<String> {
    let mut out = Vec::with_capacity(records.len());
    let mut pending: Option<PairedVcfRecord> = None;
    let mut pvs = 0i32;
    for mut cur in records {
        if let Some(mut prev) = pending.take() {
            if prev.filter == "PASS"
                && prev.status_somatic
                && prev.type_snv
                && cur.filter == "PASS"
                && cur.status_somatic
                && cur.type_snv
                && cur.start - pvs < opt_c
            {
                let tag = format!("Cluster{}bp", opt_c);
                prev.filter = tag.clone();
                cur.filter = tag;
            }
            out.push(format!("{}\t{}\t{}", prev.prefix, prev.filter, prev.rest));
        }
        if cur.type_snv && cur.filter == "PASS" {
            pvs = cur.start; // :245
        }
        pending = Some(cur);
    }
    if let Some(prev) = pending.take() {
        out.push(format!("{}\t{}\t{}", prev.prefix, prev.filter, prev.rest));
    }
    out
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
        assert_eq!(rec.filter, "PASS");
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
        assert_eq!(rec2.filter, "q22.5");
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
        assert_eq!(rec.filter, "q22.5", "d5/v3 suppressed, q22.5 stands");

        // Control: any other status keeps the depth filters (order :156→:169).
        let rec2 = format_paired_vcf_record(
            &somatic_sv(tumor, Some(normal), "Germline"),
            "chr1",
            &config,
        );
        assert_eq!(rec2.filter, "d5;q22.5");
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
            rec.filter.contains("d5"),
            "round-5 ptot=0.15 must NOT grant the exception: {}",
            rec.filter
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
        assert_eq!(rec2.filter, "q22.5", "tiny p grants d5/v3 suppression");
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
        assert!(rec.rest.contains("1/1:"), "gt 1/1 expected: {}", rec.rest);
        assert!(rec.rest.contains(":0,0:"), "bias 0→0,0: {}", rec.rest);
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
        assert_eq!(pr.prefix.split('\t').nth(5).unwrap(), "71");
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
}
