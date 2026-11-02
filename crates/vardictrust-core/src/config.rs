//! Configuration — mirrors VarDictJava's `Configuration` (§4.2 defaults).

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
        }
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
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
        } else if rest.starts_with('0') && rest.len() > 1 {
            (8u32, &rest[1..])
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

/// Extract sample name from a BAM file path (or `bam1|bam2` specification),
/// matching `VarDictLauncher.java:75-102` and `Patterns.java:65-67`:
/// 1. `SAMPLE_PATTERN` (`([^/\\._]+).sorted[^/]*\.bam`) or custom `-n` regex (with leading/trailing `/` stripped)
///    applied to the full `bam_path` (`conf.bam.getBamRaw()`).
/// 2. `SAMPLE_PATTERN2` (`([^/]+)[_\.][^/]*bam`) applied to the full `bam_path`.
/// 3. Fallback `"null"` (`VarDictLauncher.java:76`).
pub fn extract_sample_name(bam_path: &str, sample_regexp: Option<&str>) -> String {
    static SAMPLE_RE1: std::sync::LazyLock<Option<regex::Regex>> =
        std::sync::LazyLock::new(|| regex::Regex::new(r"([^/\._]+).sorted[^/]*\.bam").ok());
    static SAMPLE_RE2: std::sync::LazyLock<Option<regex::Regex>> =
        std::sync::LazyLock::new(|| regex::Regex::new(r"([^/]+)[_\.][^/]*bam").ok());

    if let Some(re_str) = sample_regexp.filter(|s| !s.is_empty()) {
        let mut pat = re_str;
        if let Some(s) = pat.strip_prefix('/') {
            pat = s;
        }
        if let Some(s) = pat.strip_suffix('/') {
            pat = s;
        }
        if let Ok(re) = regex::Regex::new(pat) {
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
    }
}
