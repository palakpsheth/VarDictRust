//! Configuration — mirrors VarDictJava's `Configuration` (§4.2 defaults).

/// Opt-in VCF post-filter numeric override carrying an optional normalized CLI token
/// (<= 64 ASCII bytes) so `##FILTER=<ID=...>` and record `FILTER` tags reproduce Perl
/// dualvar interpolation verbatim (e.g. `--vcf-qmean 20.0` -> `q20.0`, programmatic
/// `VcfThreshold::new(20.0)` -> `q20`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfThreshold<T> {
    pub value: T,
    pub raw: Option<String>,
}

impl<T> VcfThreshold<T> {
    #[inline]
    pub const fn new(value: T) -> Self {
        Self { value, raw: None }
    }

    #[inline]
    pub fn with_raw(value: T, raw: String) -> Self {
        Self {
            value,
            raw: Some(raw),
        }
    }
}

impl<T> From<T> for VcfThreshold<T> {
    #[inline]
    fn from(value: T) -> Self {
        Self::new(value)
    }
}

/// Central configuration struct.
///
/// All VarDictJava defaults from §4.2 are preserved here. Constructed from
/// the CLI args via `From<Cli>`.
#[derive(Debug, Clone)]
pub struct Config {
    pub exception_counter: std::sync::Arc<std::sync::atomic::AtomicUsize>,
    /// Reference FASTA path.
    pub reference: String,
    /// BAM/CRAM path(s). Somatic: "tumor.bam|normal.bam".
    pub bam: String,
    /// Allele frequency threshold (default: 0.01).
    pub freq: f64,
    /// Sample name(s).
    pub sample: String,
    /// True when `-N` was explicitly passed on the CLI (`Configuration.sample` != null before `initResources`).
    pub sample_explicit: bool,
    /// Region string (if -R provided).
    pub region: Option<String>,
    /// BED file path (positional).
    pub bed: Option<String>,
    /// Delimiter for BED files (default: "\t")
    pub delimiter: String,
    /// BED column indices (0-based internally).
    pub col_chr: usize,
    pub col_start: usize,
    pub col_end: usize,
    pub col_seg_start: usize,
    pub col_seg_end: usize,
    pub col_gene: Option<usize>,
    /// Track if column for chromosome was explicitly set
    pub column_for_chromosome: Option<usize>,
    /// Thread count.
    pub threads: usize,
    /// Local realignment (0 or 1).
    pub realign: u8,
    /// Minimum variant reads (default: 2).
    pub min_reads: u32,
    /// Base quality threshold (default: 22.5).
    pub base_qual: f64,
    /// Mapping quality filter.
    pub map_qual: Option<u32>,
    /// Read position filter (default: 5.0).
    pub read_pos: f64,
    /// Max mismatches (default: 8).
    pub max_mismatches: u32,
    /// Minimum match+insertion read length; reads shorter are ignored (Java -M, default: 0).
    pub minmatch: u32,
    /// Low frequency threshold for somatic STATUS determination (default: 0.05, -V flag).
    pub lofreq: f64,
    /// Disable SV calling.
    pub no_sv: bool,
    /// SV minimum length (default: 1000).
    pub sv_min_len: u32,
    /// Emit intermediate .var instead of VCF.
    pub emit_intermediate: bool,
    /// Include Fisher p-value/OR columns.
    pub fisher: bool,
    /// Tile size for bounded memory (None = auto).
    pub tile_size: Option<u32>,
    /// Amplicon mode string (e.g. "10:0.95").
    pub amplicon: Option<String>,

    // -- Additional flags added in T027-T032 --
    /// Extension for BED region (default: 0).
    pub extension: i32,
    /// Number of bases to extend for reads (default: 2).
    pub extension_nuc: i32,
    /// Include 3' UTR.
    pub three_prime: bool,
    /// Minimum # of reads to call with (default: 2).
    pub min_reads_bias: u32,
    /// qratio threshold (default: 1.5).
    pub qratio: f64,
    /// MapQ for output (`Configuration.outputMapq: Double`, default: 0.0).
    pub out_mapq: f64,
    /// Indel size (default: 50).
    pub indel_size: u32,
    /// Zero-based BED coordinates.
    pub zero_based: bool,
    /// True when `-z` was explicitly supplied on the command line (`Configuration.isZeroBasedDefined()`).
    pub zero_based_defined: bool,
    /// SAM flag filter (default: 0x504).
    pub flag_filter: u32,
    /// Raw `-F` CLI string (`Configuration.samfilter: String`, default: `"0x504"`).
    pub samfilter_raw: String,
    /// SV parameters.
    pub sv_read_len: u32,
    pub sv_flanking: u32,
    pub sv_min_af: u32,
    pub max_region_len: u32,
    /// Mono/di frequency thresholds.
    pub mfreq: f64,
    pub nmfreq: f64,
    /// Chimeric read handling.
    pub chimeric: bool,
    /// Delete/dup variant handling.
    pub deldupvar: bool,
    /// Splice mode.
    pub splice: bool,
    /// Pileup-only mode.
    pub pileup: bool,
    /// Debug mode (-D): append per-read debug columns to intermediate output.
    pub debug: bool,
    /// Debug print target.
    pub debug_print: Option<String>,
    /// CRISPR parameters.
    pub crispr_cutting_site: Option<i32>,
    pub crispr_flanking_size: Option<i32>,
    /// Column header output.
    pub column_header: bool,
    /// Dedup flag.
    pub dedup: bool,
    /// List of adaptor sequences.
    pub adaptor: Vec<String>,
    /// Verbose mode.
    pub verbose: bool,
    /// The literal -f value for the VCF writer's $Freq (var2vcf_valid.pl /
    /// var2vcf_paired.pl get it from their own command line, independent of the
    /// engine's pileup override that sets `freq` to -1, CmdParser.java:134-136).
    pub af_threshold: f64,
    /// True when `-f` was given as a non-`"0"` zero string (e.g. `-f 0.0` or `-f 0.00`),
    /// which is truthy in Perl (`var2vcf_valid.pl:21` `my $Freq = $opt_f ? $opt_f : 0.02;`)
    /// whereas `-f 0` (`"0"`) is falsy and falls back to `0.02`.
    pub af_threshold_truthy_zero: bool,
    /// unique mode alignment enabled (-u)
    pub unique_mode_alignment_enabled: bool,
    /// unique mode second in pair enabled (-UN)
    pub unique_mode_second_in_pair_enabled: bool,
    /// Omit END tag in VCF output (matches var2vcf_valid.pl -E)
    pub omit_end_tag: bool,
    /// Include Ns in depth calculation (-K)
    pub include_n: bool,
    /// Opt-in output file path (`--output <PATH>`; plain file or `.gz`/`.bgz` BGZF).
    pub output_path: Option<String>,
    /// Write a Tabix (.tbi) index alongside a BGZF-compressed --output .vcf.gz file (--tbi / --index).
    pub index_output: bool,
    /// Opt-in JSON execution telemetry path (`--stats-json <PATH>`).
    pub stats_json_path: Option<String>,
    /// Trim bases after N bases from 5' end (`-T` / `--trim-bases-after`, default: 0 = disabled).
    pub trim_bases_after: i32,
    /// Sample name regular expression (`-n` / `--sample-regexp`, with outer `/` stripped).
    pub sample_regexp: Option<String>,
    /// Write a CSI (`.csi`) index alongside a BGZF-compressed `--output .vcf.gz` file (`--csi`).
    pub index_csi: bool,
    /// Emit only PASS variants in VCF output (`--pass-only`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -S`).
    pub vcf_pass_only: bool,
    /// Disable germline rescue in paired somatic VCF filtering (`--paired-strict`, ports `var2vcf_paired.pl -M`).
    pub vcf_paired_strict: bool,
    /// Maximum p-value threshold for high-confidence somatic calls (`--vcf-pval`, default: 0.05, ports `var2vcf_paired.pl -P`).
    pub vcf_pval: f64,
    /// Normalized raw CLI token supplied to `--vcf-pval` for verbatim `P$PVAL` formatting (`None` when unset).
    pub vcf_pval_raw: Option<String>,
    /// Minimum split reads required for structural variants in single-sample VCF (`--min-splitreads`, default: 1, ports `var2vcf_valid.pl -T`).
    pub vcf_min_splitreads: i32,
    /// Enable the PSTD=0 single-position read filter in single-sample VCF output (`true` by default; `--no-pstd-filter` sets `false`, ports `var2vcf_valid.pl -P 0`).
    pub vcf_pstd_filter: bool,
    /// Opt-in VCF allele-frequency filter/GT threshold override (`--vcf-af`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -f`).
    pub vcf_af: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum total read depth override (`--vcf-min-depth`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -d`).
    pub vcf_min_depth: Option<i32>,
    /// Opt-in VCF minimum variant depth override (`--vcf-min-vd`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -v`).
    pub vcf_min_vd: Option<i32>,
    /// Opt-in VCF minimum mean base quality override (`--vcf-qmean`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -q`).
    pub vcf_qmean: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum mean mapping quality override (`--vcf-mapq`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -Q`).
    pub vcf_mapq: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum mean read-position override (`--vcf-pmean`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -p`).
    pub vcf_pmean: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum signal-to-noise (`qratio`) override (`--vcf-sn`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -o`).
    pub vcf_sn: Option<VcfThreshold<f64>>,
    /// Opt-in VCF maximum mean mismatches (`NM`) override (`--vcf-max-nm`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -m`).
    pub vcf_max_nm: Option<VcfThreshold<f64>>,
    /// Opt-in VCF homozygous `1/1` allele-frequency margin override (`--vcf-gt-freq`, `1 - AF < F`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -F`).
    pub vcf_gt_freq: Option<f64>,
    /// Opt-in VCF non-monomer MSI repeat-count threshold override (`--vcf-msi`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -I`).
    pub vcf_msi: Option<VcfThreshold<f64>>,
    /// Opt-in VCF SNV cluster window in bp (`--vcf-cluster-dist`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -c`).
    pub vcf_cluster_dist: Option<i32>,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            exception_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            reference: String::new(),
            bam: String::new(),
            freq: 0.01,
            af_threshold: 0.01,
            af_threshold_truthy_zero: false,
            sample: "SAMPLE".to_string(),
            sample_explicit: false,
            region: None,
            bed: None,
            delimiter: "\t".to_string(),
            col_chr: 2,
            col_start: 6,
            col_end: 7,
            col_seg_start: 9,
            col_seg_end: 10,
            col_gene: Some(12),
            column_for_chromosome: None,
            threads: num_cpus(),
            realign: 1,
            min_reads: 2,
            base_qual: 22.5,
            map_qual: None,
            read_pos: 5.0,
            max_mismatches: 8,
            minmatch: 0,
            lofreq: 0.05,
            no_sv: false,
            sv_min_len: 1000,
            emit_intermediate: false,
            fisher: false,
            tile_size: None,
            amplicon: None,
            extension: 0,
            extension_nuc: 2,
            three_prime: false,
            min_reads_bias: 2,
            qratio: 1.5,
            out_mapq: 0.0,
            indel_size: 50,
            zero_based: false,
            zero_based_defined: false,
            flag_filter: 0x504,
            samfilter_raw: "0x504".to_string(),
            sv_read_len: 300,
            sv_flanking: 100,
            sv_min_af: 4,
            max_region_len: 1200,
            mfreq: 0.25,
            nmfreq: 0.1,
            chimeric: false,
            deldupvar: false,
            splice: false,
            pileup: false,
            debug: false,
            debug_print: None,
            crispr_cutting_site: None,
            crispr_flanking_size: None,
            column_header: false,
            dedup: false,
            adaptor: Vec::new(),
            verbose: false,
            unique_mode_alignment_enabled: false,
            unique_mode_second_in_pair_enabled: false,
            omit_end_tag: true,
            include_n: false,
            output_path: None,
            index_output: false,
            stats_json_path: None,
            trim_bases_after: 0,
            sample_regexp: None,
            index_csi: false,
            vcf_pass_only: false,
            vcf_paired_strict: false,
            vcf_pval: 0.05,
            vcf_pval_raw: None,
            vcf_min_splitreads: 1,
            vcf_pstd_filter: true,
            vcf_af: None,
            vcf_min_depth: None,
            vcf_min_vd: None,
            vcf_qmean: None,
            vcf_mapq: None,
            vcf_pmean: None,
            vcf_sn: None,
            vcf_max_nm: None,
            vcf_gt_freq: None,
            vcf_msi: None,
            vcf_cluster_dist: None,
        }
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn validate_opt_vcf_f64_threshold(
    slot: &mut Option<VcfThreshold<f64>>,
    field_name: &str,
) -> Result<(), crate::types::VardictError> {
    if let Some(thr) = slot.as_mut() {
        if !thr.value.is_finite() || thr.value < 0.0 {
            return Err(crate::types::VardictError::Parse(format!(
                "Invalid {field_name} '{}': must be a finite number >= 0.0",
                thr.value
            )));
        }
        if thr.value == 0.0 {
            thr.value = 0.0;
        }
        if let Some(ref mut raw) = thr.raw {
            let (norm_val, norm_str) = validate_vcf_nonneg_float_token(raw).map_err(|reason| {
                crate::types::VardictError::Parse(format!(
                    "Invalid {field_name}.raw '{raw}': {reason}"
                ))
            })?;
            if norm_val.to_bits() != thr.value.to_bits() {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid {field_name}.raw '{raw}': parsed value {norm_val} disagrees with {field_name} '{}'",
                    thr.value
                )));
            }
            if raw.as_str() != norm_str {
                *raw = norm_str.to_string();
            }
        }
    }
    Ok(())
}

impl Config {
    /// Apply pileup mode (`-p`) parameter overrides matching `Configuration.java:162-165`.
    pub fn apply_pileup_overrides(&mut self) {
        if self.pileup {
            self.freq = -1.0;
            self.min_reads = 0;
        }
    }

    /// Ports `!conf.samfilter.equals("0")` in `SamView.java:101` (`0x100` secondary alignment skip)
    /// and `CigarParser.java:367` (`0x800` supplementary alignment skip).
    /// When `samfilter_raw` was updated from default `"0x504"`, checks `samfilter_raw != "0"`;
    /// otherwise falls back to `flag_filter != 0` so programmatic tests setting `config.flag_filter = 0` work.
    pub fn skip_secondary_and_supplementary(&self) -> bool {
        if self.samfilter_raw != "0x504" {
            self.samfilter_raw != "0"
        } else {
            self.flag_filter != 0
        }
    }

    /// Validates and canonicalizes all VCF post-filter thresholds (`vcf_pval`,
    /// `vcf_min_splitreads`, and all 11 `vcf_*` override fields), canonicalizing `-0.0` to
    /// `0.0` and verifying every `raw` CLI token via [`validate_vcf_nonneg_float_token`] +
    /// `f64::to_bits()` agreement so programmatic `&Config` callers share identical guards.
    pub fn validate_vcf_thresholds(&mut self) -> Result<(), crate::types::VardictError> {
        if self.af_threshold == 0.0 {
            self.af_threshold = 0.0;
        }
        if !self.vcf_pval.is_finite() || self.vcf_pval < 0.0 {
            return Err(crate::types::VardictError::Parse(format!(
                "Invalid vcf_pval '{}': must be a finite number >= 0.0",
                self.vcf_pval
            )));
        }
        if self.vcf_pval == 0.0 {
            self.vcf_pval = 0.0;
        }
        if let Some(ref mut raw) = self.vcf_pval_raw {
            let (p_norm, norm_raw) = validate_vcf_nonneg_float_token(raw).map_err(|reason| {
                crate::types::VardictError::Parse(format!("Invalid vcf_pval_raw '{raw}': {reason}"))
            })?;
            if p_norm.to_bits() != self.vcf_pval.to_bits() {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid vcf_pval_raw '{raw}': parsed value {p_norm} disagrees with vcf_pval '{}'",
                    self.vcf_pval
                )));
            }
            if raw.as_str() != norm_raw {
                *raw = norm_raw.to_string();
            }
        }
        if self.vcf_min_splitreads < 0 {
            return Err(crate::types::VardictError::Parse(format!(
                "Invalid vcf_min_splitreads '{}': must be an integer >= 0",
                self.vcf_min_splitreads
            )));
        }

        validate_opt_vcf_f64_threshold(&mut self.vcf_af, "vcf_af")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_qmean, "vcf_qmean")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_mapq, "vcf_mapq")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_pmean, "vcf_pmean")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_sn, "vcf_sn")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_max_nm, "vcf_max_nm")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_msi, "vcf_msi")?;

        for (field_name, val) in [
            ("vcf_min_depth", self.vcf_min_depth),
            ("vcf_min_vd", self.vcf_min_vd),
            ("vcf_cluster_dist", self.vcf_cluster_dist),
        ] {
            if let Some(v) = val {
                if v < 0 {
                    return Err(crate::types::VardictError::Parse(format!(
                        "Invalid {field_name} '{v}': must be an integer >= 0"
                    )));
                }
            }
        }

        if let Some(ref mut gt_freq) = self.vcf_gt_freq {
            if !gt_freq.is_finite() || *gt_freq < 0.0 || *gt_freq > 1.0 {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid vcf_gt_freq '{gt_freq}': must be a finite number in [0.0, 1.0]"
                )));
            }
            if *gt_freq == 0.0 {
                *gt_freq = 0.0;
            }
        }

        Ok(())
    }
}

/// Ports `java.lang.Integer.decode(String)` (`Configuration.java:155` / `SamView.java:53`):
/// - Rejects empty strings and strings with leading/trailing whitespace.
/// - Handles optional leading `+` or `-` sign.
/// - Detects radix from prefix: `0x` / `0X` / `#` -> 16, leading `0` with remaining length > 1 -> 8, else 10.
/// - Rejects a second sign character after the radix prefix (`0x-1`, `0x+1`, `#-1`).
/// - Parses into `i64` and validates `i32::MIN..=i32::MAX` (so `"0xFFFFFFFF"` is rejected).
pub fn java_integer_decode(s: &str) -> Result<i32, String> {
    if s.is_empty() || s.trim() != s {
        return Err(format!("For input string: \"{s}\""));
    }
    let (negative, rest) = if let Some(r) = s.strip_prefix('-') {
        (true, r)
    } else if let Some(r) = s.strip_prefix('+') {
        (false, r)
    } else {
        (false, s)
    };
    let (radix, digits) =
        if let Some(r) = rest.strip_prefix("0x").or_else(|| rest.strip_prefix("0X")) {
            (16u32, r)
        } else if let Some(r) = rest.strip_prefix('#') {
            (16u32, r)
        } else if let Some(r) = rest.strip_prefix('0').filter(|r| !r.is_empty()) {
            (8u32, r)
        } else {
            (10u32, rest)
        };
    if digits.is_empty() || digits.starts_with('-') || digits.starts_with('+') {
        return Err(format!("For input string: \"{s}\""));
    }
    let mag =
        i64::from_str_radix(digits, radix).map_err(|_| format!("For input string: \"{s}\""))?;
    let val = if negative { -mag } else { mag };
    i32::try_from(val).map_err(|_| format!("For input string: \"{s}\""))
}

/// Ports `org.apache.commons.cli.TypeHandler.createNumber(String).doubleValue()` (`CmdParser.java`).
/// - If `s` contains `.`: parses as Java `Double.valueOf(s)` (trims whitespace, strips optional trailing
///   `f`/`F`/`d`/`D`, rejects `NaN`/`Infinity`).
/// - Otherwise: parses as `Long.valueOf(s)` (strict base-10 signed 64-bit integer, no whitespace or exponent)
///   and converts to `f64`.
pub fn parse_java_number_f64(s: &str) -> Result<f64, String> {
    if s.contains('.') {
        let trimmed = s.trim();
        let core = trimmed
            .strip_suffix(['f', 'F', 'd', 'D'])
            .unwrap_or(trimmed);
        if core.is_empty() {
            return Err(format!("Invalid numeric value '{s}'"));
        }
        let v = core
            .parse::<f64>()
            .map_err(|_| format!("Invalid numeric value '{s}'"))?;
        if !v.is_finite() {
            return Err(format!("Invalid numeric value '{s}'"));
        }
        Ok(v)
    } else {
        if s.is_empty() || s.trim() != s {
            return Err(format!("Invalid numeric value '{s}'"));
        }
        let v = s
            .parse::<i64>()
            .map_err(|_| format!("Invalid numeric value '{s}'"))?;
        Ok(v as f64)
    }
}

/// Ports `((Number) cmd.getParsedOptionValue(...)).intValue()` (`CmdParser.java`).
pub fn parse_java_number_i32(s: &str) -> Result<i32, String> {
    if s.contains('.') {
        let v = parse_java_number_f64(s)?;
        Ok(v as i32)
    } else {
        if s.is_empty() || s.trim() != s {
            return Err(format!("Invalid integer value '{s}'"));
        }
        let v = s
            .parse::<i64>()
            .map_err(|_| format!("Invalid integer value '{s}'"))?;
        Ok(v as i32)
    }
}

/// Unsigned `u32` wrapper around [`parse_java_number_i32`] rejecting negative results.
pub fn parse_java_number_u32(s: &str) -> Result<u32, String> {
    let v = parse_java_number_i32(s)?;
    u32::try_from(v).map_err(|_| format!("Invalid non-negative integer value '{s}'"))
}

/// Unsigned `usize` wrapper around [`parse_java_number_i32`] rejecting negative results.
pub fn parse_java_number_usize(s: &str) -> Result<usize, String> {
    if s == "18446744073709551615" {
        return Ok(usize::MAX);
    }
    let v = parse_java_number_i32(s)?;
    usize::try_from(v).map_err(|_| format!("Invalid non-negative integer value '{s}'"))
}

/// Unsigned `u8` wrapper around [`parse_java_number_i32`] rejecting out-of-range results.
pub fn parse_java_number_u8(s: &str) -> Result<u8, String> {
    let v = parse_java_number_i32(s)?;
    u8::try_from(v).map_err(|_| format!("Invalid u8 value '{s}'"))
}

fn strip_outer_slashes(raw: &str) -> &str {
    let pat = raw.strip_prefix('/').unwrap_or(raw);
    pat.strip_suffix('/').unwrap_or(pat)
}

const SAMPLE_REGEX_SIZE_LIMIT: usize = 1 << 20;

fn compile_sample_regex(pat: &str) -> Result<regex::Regex, regex::Error> {
    regex::RegexBuilder::new(pat)
        .size_limit(SAMPLE_REGEX_SIZE_LIMIT)
        .build()
}

fn extract_sample_name_clean(bam_path: &str, clean_pat: Option<&str>) -> String {
    static SAMPLE_RE1: std::sync::LazyLock<Option<regex::Regex>> =
        std::sync::LazyLock::new(|| regex::Regex::new(r"([^/\._]+).sorted[^/]*\.bam").ok());
    static SAMPLE_RE2: std::sync::LazyLock<Option<regex::Regex>> =
        std::sync::LazyLock::new(|| regex::Regex::new(r"([^/]+)[_\.][^/]*bam").ok());

    if let Some(pat) = clean_pat {
        if let Ok(re) = compile_sample_regex(pat) {
            if let Some(caps) = re.captures(bam_path) {
                if let Some(m) = caps.get(1) {
                    return m.as_str().to_string();
                }
            }
        }
    } else if let Some(re1) = SAMPLE_RE1.as_ref() {
        if let Some(caps) = re1.captures(bam_path) {
            if let Some(m) = caps.get(1) {
                return m.as_str().to_string();
            }
        }
    }

    if let Some(re2) = SAMPLE_RE2.as_ref() {
        if let Some(caps) = re2.captures(bam_path) {
            if let Some(m) = caps.get(1) {
                return m.as_str().to_string();
            }
        }
    }

    "null".to_string()
}

/// Extract sample name from a BAM file path (or `bam1|bam2` specification),
/// matching `VarDictLauncher.java:75-102` and `Patterns.java:65-67`:
/// 1. `SAMPLE_PATTERN` (`([^/\\._]+).sorted[^/]*\.bam`) or custom `-n` regex (with leading/trailing `/` stripped)
///    applied to the full `bam_path` (`conf.bam.getBamRaw()`).
/// 2. `SAMPLE_PATTERN2` (`([^/]+)[_\.][^/]*bam`) applied to the full `bam_path`.
/// 3. Fallback `"null"` (`VarDictLauncher.java:76`).
pub fn extract_sample_name(bam_path: &str, sample_regexp: Option<&str>) -> String {
    let clean_pat = sample_regexp
        .map(strip_outer_slashes)
        .filter(|s| !s.is_empty());
    extract_sample_name_clean(bam_path, clean_pat)
}

/// Validates a `-n` / `--sample-regexp` CLI pattern after stripping a single
/// leading `'/'` and single trailing `'/'` (`CmdParser.java:88-92`).
///
/// Uses the same linear-time `regex` DFA compiler (`size_limit = 1 MiB`) as
/// `extract_sample_name` and `resolve_sample_names`; non-regular backtracking
/// constructs (lookaround/backreferences) fail fast per `AGENTS.md §12` rather
/// than silently falling back to `"null"`.
pub fn validate_sample_regexp(raw: &str) -> Result<Option<String>, String> {
    let pat = strip_outer_slashes(raw);
    if pat.is_empty() {
        return Ok(None);
    }
    compile_sample_regex(pat).map_err(|e| format!("Invalid sample name regex '{raw}': {e}"))?;
    Ok(Some(pat.to_string()))
}

/// Validates a non-negative VCF threshold `<FLOAT>` token without allocating.
///
/// Enforces `1..=64` trimmed ASCII float bytes (`[0-9.+-eE]`), rejects non-finite numbers
/// (`NaN`/`inf`/overflow) and negative values (including negative non-zero mantissas that
/// underflow `f64` to `-0.0` such as `"-1e-400"`), and canonicalizes signed zeros
/// (`"-0"` -> `(0.0, "0")`, `"-0.0"` / `"+0.0"` -> `(0.0, "0.0")`, `"+0.00"` -> `(0.0, "0.00")`)
/// while preserving verbatim positive formatting (`"1e-5"`, `"0.050"`, `"20.0"`).
pub fn validate_vcf_nonneg_float_token(s: &str) -> Result<(f64, &str), &'static str> {
    if s.is_empty() || s.trim() != s {
        return Err("must be a non-empty trimmed float token");
    }
    if s.len() > 64 {
        return Err("must be <= 64 ASCII bytes");
    }
    if !s
        .bytes()
        .all(|b| matches!(b, b'0'..=b'9' | b'.' | b'+' | b'-' | b'e' | b'E'))
    {
        return Err("must contain only ASCII float characters [0-9.+-eE]");
    }
    let v = s
        .parse::<f64>()
        .map_err(|_| "invalid floating-point syntax")?;
    let neg_nonzero_mantissa = s.starts_with('-')
        && s.split(['e', 'E'])
            .next()
            .is_some_and(|m| m.bytes().any(|b| b.is_ascii_digit() && b != b'0'));
    if !v.is_finite() || v < 0.0 || neg_nonzero_mantissa {
        return Err("must be a finite number >= 0.0");
    }
    if v == 0.0 {
        Ok((0.0, s.strip_prefix(['-', '+']).unwrap_or(s)))
    } else {
        Ok((v, s))
    }
}

/// Inline alias for [`validate_vcf_nonneg_float_token`] retained for `--vcf-pval` callers.
#[inline]
pub fn validate_vcf_pval_token(s: &str) -> Result<(f64, &str), &'static str> {
    validate_vcf_nonneg_float_token(s)
}

/// Validates and normalizes a `--vcf-pval <FLOAT>` (`var2vcf_paired.pl -P`) token into
/// `(f64, String)` for CLI parsing and [`Config`] validation.
pub fn parse_vcf_pval_token(s: &str) -> Result<(f64, String), String> {
    validate_vcf_nonneg_float_token(s)
        .map(|(v, norm)| (v, norm.to_string()))
        .map_err(|reason| format!("Invalid --vcf-pval value '{s}': {reason}"))
}

/// Parses a non-negative VCF threshold `<FLOAT>` flag (`flag`, e.g. `"--vcf-qmean"`) into
/// [`VcfThreshold<f64>`], preserving the normalized raw token for Perl dualvar formatting.
pub fn parse_vcf_threshold_nonneg_f64(s: &str, flag: &str) -> Result<VcfThreshold<f64>, String> {
    let (value, norm) = validate_vcf_nonneg_float_token(s)
        .map_err(|reason| format!("Invalid {flag} value '{s}': {reason}"))?;
    Ok(VcfThreshold::with_raw(value, norm.to_string()))
}

/// Parses a unit-interval `[0.0, 1.0]` VCF float flag (`flag`, e.g. `"--vcf-gt-freq"`).
/// Canonicalizes signed zero (`-0.0` -> `0.0`) and rejects values `> 1.0` or `< 0.0`.
pub fn parse_vcf_unit_f64(s: &str, flag: &str) -> Result<f64, String> {
    let (value, _) = validate_vcf_nonneg_float_token(s)
        .map_err(|reason| format!("Invalid {flag} value '{s}': {reason}"))?;
    if value > 1.0 {
        return Err(format!(
            "Invalid {flag} value '{s}': must be a finite number in [0.0, 1.0]"
        ));
    }
    Ok(value)
}

/// Parses a trimmed base-10 non-negative `i32` (`>= 0`) for VCF integer threshold flags
/// (`--min-splitreads`, `--vcf-min-depth`, `--vcf-min-vd`, `--vcf-cluster-dist`).
pub fn parse_vcf_nonneg_i32(s: &str, flag: &str) -> Result<i32, String> {
    if s.is_empty() || s.trim() != s {
        return Err(format!("Invalid {flag} value '{s}'"));
    }
    let v = s
        .parse::<i32>()
        .map_err(|e| format!("Invalid {flag} value '{s}': {e}"))?;
    if v < 0 {
        return Err(format!(
            "Invalid {flag} value '{s}': must be an integer >= 0"
        ));
    }
    Ok(v)
}

/// Resolves sample name(s) matching `VarDictLauncher.java:82-86, 206-268` and `var2vcf_paired.pl:45-48`.
/// Splits `bam_raw` ONLY on `'|'` (never on `':'`).
pub fn resolve_sample_names(
    bam_raw: &str,
    explicit_sample: Option<&str>,
    sample_regexp: Option<&str>,
    has_region: bool,
    emit_intermediate_or_splice: bool,
) -> String {
    let clean_pat = sample_regexp
        .map(strip_outer_slashes)
        .filter(|s| !s.is_empty());
    if let Some((bam1, bam2)) = bam_raw.split_once('|') {
        if has_region && emit_intermediate_or_splice {
            let mut base1 = match explicit_sample {
                Some(exp) if clean_pat.is_some() => exp.to_string(),
                Some(exp) => exp.split('|').next().unwrap_or(exp).to_string(),
                None => extract_sample_name_clean(bam_raw, clean_pat),
            };
            if let Some(pat) = clean_pat {
                if let Ok(re) = compile_sample_regex(pat) {
                    if let Some(caps) = re.captures(bam1) {
                        if let Some(m) = caps.get(1) {
                            base1 = m.as_str().to_string();
                        }
                    }
                }
            }
            base1
        } else if let Some(exp) = explicit_sample {
            exp.to_string()
        } else if emit_intermediate_or_splice {
            extract_sample_name_clean(bam_raw, clean_pat)
        } else {
            format!(
                "{}|{}",
                extract_sample_name_clean(bam1, clean_pat),
                extract_sample_name_clean(bam2, clean_pat)
            )
        }
    } else if let Some(exp) = explicit_sample {
        exp.to_string()
    } else {
        extract_sample_name_clean(bam_raw, clean_pat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_vardict() {
        let c = Config::default();
        assert!((c.freq - 0.01).abs() < f64::EPSILON);
        assert_eq!(c.min_reads, 2);
        assert!((c.base_qual - 22.5).abs() < f64::EPSILON);
        assert_eq!(c.realign, 1);
        assert_eq!(c.flag_filter, 0x504);
        assert_eq!(c.samfilter_raw, "0x504");
        assert!((c.mfreq - 0.25).abs() < f64::EPSILON);
        assert!((c.nmfreq - 0.1).abs() < f64::EPSILON);
        assert_eq!(c.col_chr, 2);
        assert_eq!(c.col_start, 6);
        assert_eq!(c.col_end, 7);
        assert_eq!(c.col_seg_start, 9);
        assert_eq!(c.col_seg_end, 10);
        assert_eq!(c.col_gene, Some(12));
        assert_eq!(c.column_for_chromosome, None);
        assert_eq!(c.trim_bases_after, 0);
        assert_eq!(c.sample_regexp, None);
        assert!(!c.index_csi);
    }

    #[test]
    fn java_integer_decode_table() {
        assert_eq!(java_integer_decode("0x504").unwrap(), 0x504);
        assert_eq!(java_integer_decode("0X504").unwrap(), 0x504);
        assert_eq!(java_integer_decode("#504").unwrap(), 0x504);
        assert_eq!(java_integer_decode("1284").unwrap(), 1284);
        assert_eq!(java_integer_decode("+1284").unwrap(), 1284);
        assert_eq!(java_integer_decode("02404").unwrap(), 0o2404);
        assert_eq!(java_integer_decode("0").unwrap(), 0);
        assert_eq!(java_integer_decode("00").unwrap(), 0);
        assert_eq!(java_integer_decode("-0x10").unwrap(), -16);
        assert_eq!(java_integer_decode("-010").unwrap(), -8);
        assert!(java_integer_decode("0x-1").is_err());
        assert!(java_integer_decode("0x+1").is_err());
        assert!(java_integer_decode("08").is_err());
        assert!(java_integer_decode("0xFFFFFFFF").is_err());
        assert!(java_integer_decode(" 10").is_err());
        assert!(java_integer_decode("").is_err());
    }

    #[test]
    fn extract_sample_name_and_pileup_overrides() {
        assert_eq!(
            extract_sample_name("/path/to/MySample.sorted.bam", None),
            "MySample"
        );
        assert_eq!(
            extract_sample_name("sample_ABC_123.bam", Some(r"sample_([A-Z]+)_\d+")),
            "ABC"
        );
        assert_eq!(extract_sample_name("", None), "null");

        let mut cfg = Config::default();
        cfg.apply_pileup_overrides();
        assert!((cfg.freq - 0.01).abs() < f64::EPSILON);
        assert_eq!(cfg.min_reads, 2);

        cfg.pileup = true;
        cfg.apply_pileup_overrides();
        assert!((cfg.freq - (-1.0)).abs() < f64::EPSILON);
        assert_eq!(cfg.min_reads, 0);
    }

    #[test]
    fn extract_sample_name_java_exact() {
        assert_eq!(
            extract_sample_name("/dir/foo_bar.sorted.dup.bam", None),
            "bar"
        );
        assert_eq!(extract_sample_name("/dir/foo.sorted.dup.bam", None), "foo");
        assert_eq!(
            extract_sample_name("/dir/fooXsorted.bam", None),
            "foo",
            "SAMPLE_PATTERN has unescaped '.' before 'sorted'"
        );
        assert_eq!(
            extract_sample_name("aligned_bams/NA12878.bam", None),
            "aligned",
            "SAMPLE_PATTERN2 matches 'aligned_bam' in directory name when applied to raw -b path"
        );
        assert_eq!(extract_sample_name("/dir/foo_bar.bam", None), "foo_bar");
        assert_eq!(extract_sample_name("/dir/foo.cram", None), "null");
        assert_eq!(
            extract_sample_name("t.bam|n.bam", None),
            "t.bam|n",
            "SAMPLE_PATTERN2 is greedy across `|` when applied to full bamRaw"
        );
        assert_eq!(
            extract_sample_name("/dir/sample_42.bam", Some("/sample_(\\d+)/")),
            "42"
        );
        // -n "/" and -n "//" strip outer slashes first to "" and fall through to SAMPLE_RE1
        assert_eq!(
            extract_sample_name("/dir/foo.sorted.dup.bam", Some("/")),
            "foo"
        );
        assert_eq!(
            extract_sample_name("/dir/foo.sorted.dup.bam", Some("//")),
            "foo"
        );

        // validate_sample_regexp checks
        assert_eq!(validate_sample_regexp(""), Ok(None));
        assert_eq!(validate_sample_regexp("/"), Ok(None));
        assert_eq!(validate_sample_regexp("//"), Ok(None));
        assert_eq!(
            validate_sample_regexp("/sample_(\\d+)/"),
            Ok(Some("sample_(\\d+)".to_string()))
        );
        assert!(validate_sample_regexp("([unclosed").is_err());
        assert!(validate_sample_regexp(r"(?<=/)([^/]+)\.bam$").is_err());

        // resolve_sample_names matrix (paired -R .var/-i, paired BED .var/-i, paired VCF, single)
        assert_eq!(
            resolve_sample_names("t.sorted.bam|n.sorted.bam", None, None, true, true),
            "t"
        );
        assert_eq!(
            resolve_sample_names("t.bam|n.bam", None, None, true, true),
            "t.bam|n"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam:t_2.bam|n_1.bam",
                None,
                Some(r":t_(\d+)\.bam"),
                true,
                true
            ),
            "2"
        );
        assert_eq!(
            resolve_sample_names("t_1.bam|n_2.bam", Some("EXP_T|EXP_N"), None, true, true),
            "EXP_T"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam|n_2.bam",
                Some("EXP_T|EXP_N"),
                Some(r"t_(\d+)"),
                true,
                true
            ),
            "1"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam|n_2.bam",
                Some("EXP_T|EXP_N"),
                Some(r"t_(\d+)"),
                false,
                true
            ),
            "EXP_T|EXP_N"
        );
        assert_eq!(
            resolve_sample_names("t.bam|n.bam", None, None, false, true),
            "t.bam|n"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam|n_2.bam",
                Some("EXP_T|EXP_N"),
                Some(r"t_(\d+)"),
                true,
                false
            ),
            "EXP_T|EXP_N"
        );
        assert_eq!(
            resolve_sample_names(
                "dir/T1.sorted.bam|dir/N1.sorted.bam",
                None,
                None,
                false,
                false
            ),
            "T1|N1"
        );
        assert_eq!(
            resolve_sample_names(
                "dir/s_T9.bam|dir/s_N8.bam",
                None,
                Some("s_([A-Z0-9]+)\\.bam"),
                true,
                false
            ),
            "T9|N8"
        );
        assert_eq!(
            resolve_sample_names(
                "dir/s_T9.bam",
                Some("EXPLICIT"),
                Some("s_([A-Z0-9]+)"),
                false,
                false
            ),
            "EXPLICIT"
        );
        assert_eq!(
            resolve_sample_names("a_1.bam:b_2.bam", None, Some(r":b_(\d+)\.bam"), false, true),
            "2"
        );
        // Single outer-slash strip (`CmdParser.java:86-93`): escaped trailing `\/` and `//...//` anchors
        assert_eq!(
            extract_sample_name("/mnt/cohort/SAMP_ESC/aln.bam", Some(r"/cohort\/([^/]+)\//")),
            "SAMP_ESC"
        );
        assert_eq!(
            resolve_sample_names(
                "/mnt/cohort/T_ESC/t.bam|/mnt/cohort/N_ESC/n.bam",
                None,
                Some(r"/cohort\/([^/]+)\//"),
                false,
                false,
            ),
            "T_ESC|N_ESC"
        );
        assert_eq!(
            resolve_sample_names(
                "/mnt/cohort/T_ESC/t.bam|/mnt/cohort/N_ESC/n.bam",
                None,
                Some(r"/cohort\/([^/]+)\//"),
                true,
                true,
            ),
            "T_ESC"
        );
        assert_eq!(
            extract_sample_name("/cohort/ANCHORED/aln.bam", Some(r"//cohort/([^/]+)//")),
            "ANCHORED"
        );
        assert_eq!(
            extract_sample_name("not_cohort/FAKE/aln_1.bam", Some(r"//cohort/([^/]+)//")),
            "aln_1",
            "double-slash pattern strips once to `/cohort/([^/]+)/` and must not match `not_cohort/`"
        );
    }

    #[test]
    fn vcf_threshold_constructors_and_defaults() {
        let c = Config::default();
        assert!(c.vcf_af.is_none());
        assert!(c.vcf_min_depth.is_none());
        assert!(c.vcf_min_vd.is_none());
        assert!(c.vcf_qmean.is_none());
        assert!(c.vcf_mapq.is_none());
        assert!(c.vcf_pmean.is_none());
        assert!(c.vcf_sn.is_none());
        assert!(c.vcf_max_nm.is_none());
        assert!(c.vcf_gt_freq.is_none());
        assert!(c.vcf_msi.is_none());
        assert!(c.vcf_cluster_dist.is_none());

        let t_new = VcfThreshold::new(20.0_f64);
        assert_eq!(t_new.value, 20.0);
        assert_eq!(t_new.raw, None);

        let t_from: VcfThreshold<f64> = 20.0.into();
        assert_eq!(t_from, t_new);

        let t_raw = VcfThreshold::with_raw(20.0, "20.0".to_string());
        assert_eq!(t_raw.value, 20.0);
        assert_eq!(t_raw.raw.as_deref(), Some("20.0"));
    }

    #[test]
    fn validate_vcf_nonneg_float_token_and_parsers() {
        let (v0, r0) = validate_vcf_nonneg_float_token("+0.00").unwrap();
        assert_eq!(v0.to_bits(), 0.0_f64.to_bits());
        assert_eq!(r0, "0.00");

        let (vn0, rn0) = validate_vcf_nonneg_float_token("-0").unwrap();
        assert_eq!(vn0.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rn0, "0");

        let (vn00, rn00) = validate_vcf_nonneg_float_token("-0.0").unwrap();
        assert_eq!(vn00.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rn00, "0.0");

        let (v_sci, r_sci) = validate_vcf_nonneg_float_token("1e-5").unwrap();
        assert_eq!(v_sci.to_bits(), 1e-5_f64.to_bits());
        assert_eq!(r_sci, "1e-5");

        let exact_64 = format!("0.{}", "0".repeat(62));
        assert_eq!(exact_64.len(), 64);
        assert!(validate_vcf_nonneg_float_token(&exact_64).is_ok());

        let over_64 = format!("0.{}", "0".repeat(63));
        assert_eq!(over_64.len(), 65);
        assert!(validate_vcf_nonneg_float_token(&over_64).is_err());

        for bad in [
            "", " 0.05", "0.05 ", "0.05\n", "0.05\t", "0.05;bad", "-0.01", "-1e-400", "-3", "NaN",
            "inf", "Infinity", "1e500", "0.01f", "-", "+", ".",
        ] {
            assert!(
                validate_vcf_nonneg_float_token(bad).is_err(),
                "expected Err for {bad:?}"
            );
        }

        let thr = parse_vcf_threshold_nonneg_f64("+0.00", "--vcf-af").unwrap();
        assert_eq!(thr.value.to_bits(), 0.0_f64.to_bits());
        assert_eq!(thr.raw.as_deref(), Some("0.00"));

        assert_eq!(
            parse_vcf_unit_f64("-0.0", "--vcf-gt-freq")
                .unwrap()
                .to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(parse_vcf_unit_f64("1.0", "--vcf-gt-freq").unwrap(), 1.0);
        assert!(parse_vcf_unit_f64("1.0001", "--vcf-gt-freq").is_err());
        assert!(parse_vcf_unit_f64("-0.1", "--vcf-gt-freq").is_err());

        assert_eq!(parse_vcf_nonneg_i32("0", "--vcf-min-depth").unwrap(), 0);
        assert_eq!(
            parse_vcf_nonneg_i32("15", "--vcf-cluster-dist").unwrap(),
            15
        );
        for bad_int in ["", " 1", "-1", "-3", "2.5", "0x2"] {
            assert!(
                parse_vcf_nonneg_i32(bad_int, "--vcf-min-depth").is_err(),
                "expected Err for {bad_int:?}"
            );
        }
    }

    #[test]
    fn config_validate_vcf_thresholds_checks_all_fields() {
        let mut cfg = Config {
            af_threshold: -0.0,
            vcf_pval: -0.0,
            vcf_pval_raw: Some("-0.0".to_string()),
            vcf_af: Some(VcfThreshold::with_raw(-0.0, "+0.00".to_string())),
            vcf_min_depth: Some(0),
            vcf_min_vd: Some(0),
            vcf_qmean: Some(VcfThreshold::with_raw(20.0, "20.0".to_string())),
            vcf_mapq: Some(VcfThreshold::with_raw(0.0, "-0".to_string())),
            vcf_pmean: Some(VcfThreshold::with_raw(7.5, "7.5".to_string())),
            vcf_sn: Some(VcfThreshold::new(2.0)),
            vcf_max_nm: Some(VcfThreshold::with_raw(4.0, "4.0".to_string())),
            vcf_gt_freq: Some(-0.0),
            vcf_msi: Some(VcfThreshold::with_raw(10.5, "10.5".to_string())),
            vcf_cluster_dist: Some(15),
            ..Config::default()
        };
        cfg.validate_vcf_thresholds().unwrap();
        assert_eq!(cfg.af_threshold.to_bits(), 0.0_f64.to_bits());
        assert_eq!(cfg.vcf_pval.to_bits(), 0.0_f64.to_bits());
        assert_eq!(cfg.vcf_pval_raw.as_deref(), Some("0.0"));
        assert_eq!(
            cfg.vcf_af.as_ref().unwrap().value.to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(cfg.vcf_af.as_ref().unwrap().raw.as_deref(), Some("0.00"));
        assert_eq!(cfg.vcf_mapq.as_ref().unwrap().raw.as_deref(), Some("0"));
        assert_eq!(cfg.vcf_gt_freq.unwrap().to_bits(), 0.0_f64.to_bits());

        // Reject negative integer overrides, >1.0 gt_freq, mismatched/hostile raw tokens
        let mut bad_depth = Config {
            vcf_min_depth: Some(-1),
            ..Config::default()
        };
        assert!(bad_depth.validate_vcf_thresholds().is_err());

        let mut bad_vd = Config {
            vcf_min_vd: Some(-2),
            ..Config::default()
        };
        assert!(bad_vd.validate_vcf_thresholds().is_err());

        let mut bad_cluster = Config {
            vcf_cluster_dist: Some(-5),
            ..Config::default()
        };
        assert!(bad_cluster.validate_vcf_thresholds().is_err());

        let mut bad_gt = Config {
            vcf_gt_freq: Some(1.1),
            ..Config::default()
        };
        assert!(bad_gt.validate_vcf_thresholds().is_err());

        let mut bad_raw_mismatch = Config {
            vcf_qmean: Some(VcfThreshold::with_raw(20.0, "25.0".to_string())),
            ..Config::default()
        };
        assert!(bad_raw_mismatch.validate_vcf_thresholds().is_err());

        let mut bad_raw_inject = Config {
            vcf_qmean: Some(VcfThreshold::with_raw(20.0, "20.0\n##evil".to_string())),
            ..Config::default()
        };
        assert!(bad_raw_inject.validate_vcf_thresholds().is_err());

        let mut bad_raw_over64 = Config {
            vcf_msi: Some(VcfThreshold::with_raw(0.0, format!("0.{}", "0".repeat(63)))),
            ..Config::default()
        };
        assert!(bad_raw_over64.validate_vcf_thresholds().is_err());
    }
}
