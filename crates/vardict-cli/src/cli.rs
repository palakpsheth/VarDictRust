//! CLI argument parser — drop-in compatible with VarDictJava flags plus
//! ergonomic extras.
//!
//! ports Main.java / CmdParser.java

use clap::Parser;

/// VarDictJava v1.8.3 flags (CmdParser.java) with no v1.0 implementation.
/// Must be rejected loudly (§12) before Clap sees them.
const UNSUPPORTED_FLAGS: &[&str] = &[
    "-n", // sample-name regular expression
    "-C", // chromosome names are just numbers (deprecated upstream)
    "-T",
    "--trim",       // trim bases from reads
    "--VS",         // SAM validation stringency
    "-Z",           // downsampling fraction (CmdParser.java:125-127, 433-440)
    "--downsample", // long form of -Z
];

/// Rewrite Java-style multi-char single-dash options to their clap long forms.
/// commons-cli accepts `-th` (CmdParser.java:480) and `--th`; clap has no
/// multi-char shorts, so normalize both to `--threads` before parsing.
pub fn normalize_args(args: &[String]) -> Vec<String> {
    args.iter()
        .map(|arg| {
            if arg == "-th" || arg == "--th" {
                "--threads".to_string()
            } else if let Some(v) = arg.strip_prefix("-th=") {
                format!("--threads={v}")
            } else if arg == "-UN" {
                "--UN".to_string()
            } else {
                arg.clone()
            }
        })
        .collect()
}

/// Fail loudly if any known-unimplemented VarDictJava flag was supplied.
/// Runs before Clap parsing, which would otherwise emit a generic
/// "unexpected argument" error instead of the actionable §12 message.
pub fn check_unsupported(args: &[String]) -> Result<(), vardict_core::types::VardictError> {
    for arg in args {
        let flag = arg.split_once('=').map_or(arg.as_str(), |(f, _)| f);
        if UNSUPPORTED_FLAGS.contains(&flag) {
            return Err(vardict_core::types::VardictError::Unsupported(format!(
                "{flag} not implemented in v1.0; use VarDictJava"
            )));
        }
    }
    Ok(())
}

/// VarDictRust: single-binary VarDict variant caller (VarDictJava 1.8.3-compatible).
#[derive(Parser, Debug)]
#[command(
    name = "vardictrust",
    version,
    long_version = concat!(
        env!("CARGO_PKG_VERSION"), "\n",
        "Idiomatic Rust rewrite of VarDictJava v1.8.3 (commit 0886261; Perl/R submodule 009e017).\n",
        "Output validated against VarDictJava v1.8.3 - see README \"Validation status\".\n",
        "Original: VarDict/VarDictJava, AstraZeneca-NGS - https://github.com/AstraZeneca-NGS/VarDictJava\n",
        "Please cite: Lai Z, et al. VarDict: a novel and versatile variant caller for next-generation\n",
        "  sequencing in cancer research. Nucleic Acids Res. 2016;44(11):e108. doi:10.1093/nar/gkw227\n",
        "AI disclosure: written with AI coding assistants under human direction; see README \"AI Disclosure\".\n",
        "License: MIT"
    ),
    about,
    disable_version_flag = true,
    // `-h` belongs to Java's column header (CmdParser.java:223, 67), so clap
    // help moves to `--help`/`-H` (CmdParser.java:222).
    disable_help_flag = true,
    args_override_self = true
)]
pub struct Cli {
    /// Print version.
    #[arg(long = "version", action = clap::ArgAction::Version)]
    _version: (),

    /// Print help.
    #[arg(short = 'H', long = "help", action = clap::ArgAction::Help)]
    _help: (),

    /// Reference FASTA (.fai required). Aliases: --reference. Accepts hg19/hg38/mm10 shortcuts.
    #[arg(short = 'G', long = "reference")]
    pub reference: String,

    /// Indexed BAM/CRAM. Somatic: "tumor.bam|normal.bam". Alias: --bam.
    #[arg(short = 'b', long = "bam")]
    pub bam: String,

    /// Allele-frequency threshold. Alias: --af-threshold. [default: 0.01]
    #[arg(short = 'f', long = "af-threshold", default_value_t = 0.01)]
    pub freq: f64,

    /// Sample name (somatic: "tumor|normal"). Alias: --name.
    #[arg(short = 'N', long = "name")]
    pub sample: Option<String>,

    /// Region chr:start-end[:gene] (instead of BED).
    #[arg(short = 'R', long = "region")]
    pub region: Option<String>,

    /// BED column for chromosome (1-based).
    #[arg(short = 'c')]
    pub col_chr: Option<usize>,

    /// BED column for start (1-based).
    #[arg(short = 'S')]
    pub col_start: Option<usize>,

    /// BED column for end (1-based).
    #[arg(short = 'E')]
    pub col_end: Option<usize>,

    /// BED column for segment starts (1-based).
    #[arg(short = 's')]
    pub col_seg_start: Option<usize>,

    /// BED column for segment ends (1-based).
    #[arg(short = 'e')]
    pub col_seg_end: Option<usize>,

    /// BED column for gene (1-based).
    #[arg(short = 'g')]
    pub col_gene: Option<usize>,

    /// Threads. Bare `-th` = all cores (Java readThreadsCount
    /// CmdParser.java:617-628); explicit `-th 0` clamps to 1 via max(n,1)
    /// (CmdParser.java:162). Absent flag stays all-cores (deliberate perf-only
    /// divergence over Java's 1-thread default — task P1 report).
    /// usize::MAX is the bare-flag sentinel; alias: --threads.
    #[arg(
        long = "threads",
        num_args(0..=1),
        default_missing_value = "18446744073709551615"
    )]
    pub threads: Option<usize>,

    /// Local realignment 0/1. [default: 1]
    #[arg(short = 'k', default_value_t = 1)]
    pub realign: u8,

    /// Amplicon mode "int:float". [default: 10:0.95 when 8-col BED]
    #[arg(short = 'a', long = "amplicon")]
    pub amplicon: Option<String>,

    /// Min variant reads. [default: 2]
    #[arg(short = 'r', default_value_t = 2)]
    pub min_reads: u32,

    /// Min base quality. [default: 22.5]
    #[arg(short = 'q', default_value_t = 22.5)]
    pub base_qual: f64,

    /// Min mapping quality.
    #[arg(short = 'Q')]
    pub map_qual: Option<u32>,

    /// Mean read position filter. [default: 5]
    #[arg(short = 'P', default_value_t = 5.0)]
    pub read_pos: f64,

    /// Max mismatches. [default: 8]
    #[arg(short = 'm', default_value_t = 8)]
    pub max_mismatches: u32,

    /// Minimum match+insertion read length; shorter reads ignored. [default: 0]
    #[arg(short = 'M', default_value_t = 0)]
    pub minmatch: u32,

    /// Min normal AF for somatic. [default: 0.05]
    #[arg(short = 'V', default_value_t = 0.05)]
    pub min_normal_af: f64,

    /// Disable SV calling.
    #[arg(short = 'U', long = "nosv")]
    pub no_sv: bool,

    /// SV minimum length. [default: 1000]
    #[arg(short = 'L', default_value_t = 1000)]
    pub sv_min_len: u32,

    /// Emit the intermediate .var table instead of VCF.
    #[arg(long = "emit-intermediate")]
    pub emit_intermediate: bool,

    /// Use the in-engine Fisher test.
    #[arg(long = "fisher")]
    pub fisher: bool,

    /// Advanced: region tile size (bp) for bounded-memory WGS. [auto]
    #[arg(long = "tile-size")]
    pub tile_size: Option<u32>,

    /// Extension for BED region.
    #[arg(short = 'x', default_value_t = 0)]
    pub extension: i32,

    /// Number of nucleotide extension. [default: 2]
    #[arg(short = 'X', default_value_t = 2)]
    pub extension_nuc: i32,

    /// Shift indel positions to their 3'-most alternative alignment
    /// (ports CmdParser.java:71 moveIndelsTo3; shift3 applied in tovars.rs).
    #[arg(short = '3')]
    pub three_prime: bool,

    /// Min reads for bias. [default: 2]
    #[arg(short = 'B', default_value_t = 2)]
    pub min_reads_bias: u32,

    /// qratio filter. [default: 1.5]
    #[arg(short = 'o', default_value_t = 1.5)]
    pub qratio: f64,

    /// MapQ for output. [default: 0]
    #[arg(short = 'O', default_value_t = 0)]
    pub out_mapq: u32,

    /// Indel size. [default: 50]
    #[arg(short = 'I', default_value_t = 50)]
    pub indel_size: u32,

    /// BED coordinates are 0-based.
    #[arg(short = 'z')]
    pub zero_based: bool,

    /// SAM flag filter. [default: 0x504]
    #[arg(short = 'F', default_value = "0x504")]
    pub flag_filter: String,

    /// SV read length. [default: 300]
    #[arg(short = 'w', default_value_t = 300)]
    pub sv_read_len: u32,

    /// SV flanking. [default: 100]
    #[arg(short = 'W', default_value_t = 100)]
    pub sv_flanking: u32,

    /// SV minimum AF. [default: 4]
    #[arg(short = 'A', default_value_t = 4)]
    pub sv_min_af: u32,

    /// Max region length. [default: 1200]
    #[arg(short = 'Y', default_value_t = 1200)]
    pub max_region_len: u32,

    /// Column header output: print the .var column names before results
    /// (ports CmdParser.java:67 `printHeader` + SimpleMode.java:114-127).
    #[arg(short = 'h', long = "header")]
    pub column_header: bool,

    /// VCF format output. Defined but never read by VarDictJava v1.8.3
    /// (CmdParser.java:224 vs parseCmd 59-186); accepted and ignored for
    /// drop-in compatibility.
    #[arg(short = 'v')]
    pub vcf_output_noop: bool,

    /// Chimeric read handling.
    #[arg(long = "chimeric")]
    pub chimeric: bool,

    /// Delete/dup variant handling.
    #[arg(long = "deldupvar")]
    pub deldupvar: bool,

    /// Mono frequency threshold. [default: 0.25]
    #[arg(long = "mfreq", default_value_t = 0.25)]
    pub mfreq: f64,

    /// Normal mono frequency threshold. [default: 0.1]
    #[arg(long = "nmfreq", default_value_t = 0.1)]
    pub nmfreq: f64,

    /// CRISPR cutting site.
    #[arg(short = 'J')]
    pub crispr_cutting_site: Option<i32>,

    /// CRISPR flanking size.
    #[arg(short = 'j')]
    pub crispr_flanking_size: Option<i32>,

    /// Splice mode.
    #[arg(short = 'i', long = "splice")]
    pub splice: bool,

    /// Pileup-only mode.
    #[arg(short = 'p')]
    pub pileup: bool,

    /// Debug mode: append per-read debug columns to the intermediate output.
    #[arg(short = 'D', long = "debug")]
    pub debug: bool,

    /// Debug print destination: OUT or ERR. Bare `-DP` or an unknown value
    /// normalize to OUT (ports CmdParser.java:171-176).
    #[arg(
        long = "DP",
        num_args(0..=1),
        default_missing_value = "OUT"
    )]
    pub debug_print: Option<String>,

    /// Verbose mode.
    #[arg(short = 'y', long = "verbose")]
    pub verbose: bool,

    /// Dedup: remove duplicated reads, keeping one pair per start position
    /// (ports CmdParser.java:70, 229).
    #[arg(short = 't', long = "dedup")]
    pub dedup: bool,

    /// Filter adaptor sequences so that they are not used in realignment. Multiple adaptors can be supplied by setting them with comma.
    #[arg(long = "adaptor")]
    pub adaptor: Option<String>,

    /// The delimiter for splitting region info. [default: \t]
    #[arg(short = 'd', default_value = "\t")]
    pub delimiter: String,

    /// Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using forward read only (-u).
    #[arg(short = 'u')]
    pub unique_mode_alignment_enabled: bool,

    /// Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using first read only (-UN).
    #[arg(long = "UN")]
    pub unique_mode_second_in_pair_enabled: bool,

    /// Include Ns in depth calculation.
    #[arg(short = 'K')]
    pub include_n: bool,

    /// Positional BED file (when -R absent).
    pub bed: Option<String>,
}

fn all_cores() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

impl Cli {
    /// Convert CLI args to a core Config.
    pub fn to_config(&self) -> vardict_core::config::Config {
        let flag_filter =
            if self.flag_filter.starts_with("0x") || self.flag_filter.starts_with("0X") {
                u32::from_str_radix(&self.flag_filter[2..], 16).unwrap_or(0x504)
            } else {
                self.flag_filter.parse().unwrap_or(0x504)
            };

        let (freq, min_reads) = if self.pileup {
            (-1.0, 0)
        } else {
            (self.freq, self.min_reads)
        };

        vardict_core::config::Config {
            exception_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            reference: self.reference.clone(),
            bam: self.bam.clone(),
            freq,
            // The writer's $Freq is the literal -f value even under -p (perl
            // receives -f from the pipeline, not java's pileup override).
            af_threshold: self.freq,
            sample: self
                .sample
                .clone()
                .unwrap_or_else(|| extract_sample_name(&self.bam)),
            region: self.region.clone(),
            bed: self.bed.clone(),
            delimiter: self.delimiter.clone(),
            col_chr: self.col_chr.map(|c| c.saturating_sub(1)).unwrap_or(2),
            col_start: self.col_start.map(|s| s.saturating_sub(1)).unwrap_or(6),
            col_end: self.col_end.map(|e| e.saturating_sub(1)).unwrap_or(7),
            col_seg_start: self
                .col_seg_start
                .map(|s| s.saturating_sub(1))
                .or_else(|| self.col_start.map(|s| s.saturating_sub(1)))
                .unwrap_or(9),
            col_seg_end: self
                .col_seg_end
                .map(|e| e.saturating_sub(1))
                .or_else(|| self.col_end.map(|e| e.saturating_sub(1)))
                .unwrap_or(10),
            col_gene: self
                .col_gene
                .map(|g| {
                    if g > 0 {
                        Some(g.saturating_sub(1))
                    } else {
                        None
                    }
                })
                .unwrap_or(Some(12)),
            column_for_chromosome: self.col_chr.map(|c| c.saturating_sub(1)),
            // Java: bare -th → availableProcessors (CmdParser.java:621-623),
            // -th N (incl. explicit 0) → max(N,1) (:162). Absent flag stays
            // all-cores (perf-only divergence; no output effect — task P1
            // report).
            threads: match self.threads {
                None => all_cores(),
                Some(n) if n == usize::MAX => all_cores(),
                Some(n) => n.max(1),
            },
            realign: self.realign,
            min_reads,
            base_qual: self.base_qual,
            map_qual: self.map_qual,
            read_pos: self.read_pos,
            max_mismatches: self.max_mismatches,
            minmatch: self.minmatch,
            lofreq: self.min_normal_af,
            no_sv: self.no_sv,
            sv_min_len: self.sv_min_len,
            emit_intermediate: self.emit_intermediate,
            fisher: self.fisher,
            tile_size: self.tile_size,
            amplicon: self.amplicon.clone(),
            extension: self.extension,
            extension_nuc: self.extension_nuc,
            three_prime: self.three_prime,
            min_reads_bias: self.min_reads_bias,
            qratio: self.qratio,
            out_mapq: self.out_mapq,
            indel_size: self.indel_size,
            zero_based: self.zero_based,
            flag_filter,
            sv_read_len: self.sv_read_len,
            sv_flanking: self.sv_flanking,
            sv_min_af: self.sv_min_af,
            max_region_len: self.max_region_len,
            mfreq: self.mfreq,
            nmfreq: self.nmfreq,
            chimeric: self.chimeric,
            deldupvar: self.deldupvar,
            splice: self.splice,
            pileup: self.pileup,
            debug: self.debug,
            // ports CmdParser.java:171-176: bare -DP → "OUT"; unknown
            // values silently fall back to "OUT" (only ERR is honored).
            debug_print: self.debug_print.as_ref().map(|v| {
                if v == "ERR" {
                    v.clone()
                } else {
                    "OUT".to_string()
                }
            }),
            crispr_cutting_site: self.crispr_cutting_site,
            crispr_flanking_size: self.crispr_flanking_size,
            column_header: self.column_header,
            dedup: self.dedup,
            adaptor: self
                .adaptor
                .as_ref()
                .map(|s| s.split(',').map(|x| x.to_string()).collect())
                .unwrap_or_default(),
            verbose: self.verbose,
            unique_mode_alignment_enabled: self.unique_mode_alignment_enabled,
            unique_mode_second_in_pair_enabled: self.unique_mode_second_in_pair_enabled,
            omit_end_tag: self.col_end.is_some(),
            include_n: self.include_n,
        }
    }
}

fn extract_sample_name(bam_path: &str) -> String {
    // Split by '|' to handle somatic BAMs, and take the first BAM
    let first_bam = bam_path.split('|').next().unwrap_or(bam_path);
    // Get the filename
    let filename = if let Some(idx) = first_bam.rfind('/') {
        &first_bam[idx + 1..]
    } else {
        first_bam
    };

    if let Ok(re1) = regex::Regex::new(r"([^/\\._]+).sorted[^/]*\.bam") {
        if let Some(caps) = re1.captures(filename) {
            if let Some(m) = caps.get(1) {
                return m.as_str().to_string();
            }
        }
    }

    if let Ok(re2) = regex::Regex::new(r"([^/\\]+)[_\.][^/\\]*bam") {
        if let Some(caps) = re2.captures(filename) {
            if let Some(m) = caps.get(1) {
                return m.as_str().to_string();
            }
        }
    }

    // Fallback: use the filename itself or "SAMPLE"
    if !filename.is_empty() {
        filename.to_string()
    } else {
        "SAMPLE".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(flags: &[&str]) -> Vec<String> {
        std::iter::once("vardictrust".to_string())
            .chain(flags.iter().map(|s| s.to_string()))
            .collect()
    }

    #[test]
    fn unsupported_flag_errors() {
        assert!(check_unsupported(&args(&["-n", "SAMPLE[0-9]"])).is_err());
        assert!(check_unsupported(&args(&["--downsample=0.5"])).is_err());
        // -Z is Java's downsampling fraction (CmdParser.java:125-127, 433-440).
        assert!(check_unsupported(&args(&["-Z", "0.7"])).is_err());
        // Glued form isn't matched here, but clap rejects it loudly after we return Ok.
        assert!(check_unsupported(&args(&["-Z0.7"])).is_ok());
    }

    #[test]
    fn supported_flag_passes() {
        assert!(check_unsupported(&args(&["-f", "0.05", "-G", "ref.fa"])).is_ok());
        assert!(check_unsupported(&args(&["-D"])).is_ok());
        assert!(check_unsupported(&args(&["--debug"])).is_ok());
        assert!(check_unsupported(&args(&["-h"])).is_ok());
        assert!(check_unsupported(&args(&["-t"])).is_ok());
    }

    #[test]
    fn t_is_dedup_not_threads() {
        // ports CmdParser.java:70 `config.removeDuplicatedReads = cmd.hasOption("t")`
        let cli = Cli::try_parse_from(args(&["-G", "r.fa", "-b", "a.bam", "-t"])).unwrap();
        assert!(cli.dedup);
        assert!(cli.threads.is_none());
        assert!(!cli.column_header);
    }
    #[test]
    fn th_is_threads() {
        // ports CmdParser.java:480 (option th), 617-628 (readThreadsCount)
        let cli = Cli::try_parse_from(super::normalize_args(&args(&[
            "-G", "r.fa", "-b", "a.bam", "-th", "4",
        ])))
        .unwrap();
        assert_eq!(cli.threads, Some(4));
        let cli = Cli::try_parse_from(super::normalize_args(&args(&[
            "-G", "r.fa", "-b", "a.bam", "-th=2", "--th", "8",
        ])))
        .unwrap();
        assert_eq!(cli.threads, Some(8)); // args_override_self
    }

    #[test]
    fn bare_th_means_all_cores() {
        // ports CmdParser.java:621-623: -th without value → availableProcessors
        // (usize::MAX is the bare-flag sentinel, resolved in to_config).
        let cli = Cli::try_parse_from(super::normalize_args(&args(&[
            "-b", "a.bam", "-G", "r.fa", "-th",
        ])))
        .unwrap();
        assert_eq!(cli.threads, Some(usize::MAX));
        assert_eq!(cli.to_config().threads, all_cores());
    }

    #[test]
    fn explicit_zero_threads_means_one() {
        // ports CmdParser.java:162: config.threads = max(readThreadsCount(), 1),
        // so -th 0 clamps to a single thread (NOT all cores — Java's bare -th
        // case is the only one mapping to availableProcessors).
        let cli = Cli::try_parse_from(super::normalize_args(&args(&[
            "-b", "a.bam", "-G", "r.fa", "-th", "0",
        ])))
        .unwrap();
        assert_eq!(cli.threads, Some(0));
        assert_eq!(cli.to_config().threads, 1);
    }

    #[test]
    fn absent_threads_stays_all_cores() {
        // Deliberate divergence over Java's 1-thread default (task P1 report):
        // no -th at all → all cores, perf-only, no output effect.
        let cli = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa"])).unwrap();
        assert!(cli.threads.is_none());
        assert_eq!(cli.to_config().threads, all_cores());
    }

    #[test]
    fn dp_normalizes_bare_and_unknown_to_out() {
        // ports CmdParser.java:171-176 + :505-508: getOptionValue("DP","OUT")
        // + switch default arm — bare --DP → OUT, "ERR" → ERR,
        // anything else → OUT. ("DP" is long-only in commons-cli too.)
        for (argv, expect) in [
            (vec!["--DP"], "OUT"),
            (vec!["--DP", "ERR"], "ERR"),
            (vec!["--DP", "BOGUS"], "OUT"),
        ] {
            let mut a = vec!["-b", "a.bam", "-G", "r.fa"];
            a.extend(&argv);
            let cli = Cli::try_parse_from(super::normalize_args(&args(&a))).unwrap();
            assert_eq!(
                cli.to_config().debug_print.as_deref(),
                Some(expect),
                "argv={argv:?}"
            );
        }
    }

    #[test]
    fn h_prints_header_not_help() {
        // ports CmdParser.java:67 `config.printHeader = cmd.hasOption('h')`
        let cli = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-h"])).unwrap();
        assert!(cli.column_header);
        let cli = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "--header"])).unwrap();
        assert!(cli.column_header);
    }

    #[test]
    fn help_is_dash_h_cap_and_long() {
        // ports CmdParser.java:222 (-H/? print help page)
        use clap::error::ErrorKind;
        let err = Cli::try_parse_from(args(&["-H"])).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::DisplayHelp);
        let err = Cli::try_parse_from(args(&["--help"])).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::DisplayHelp);
    }

    #[test]
    fn v_is_accepted_and_ignored() {
        // Defined but never read by VarDictJava v1.8.3 (CmdParser.java:224 vs parseCmd 59-186).
        let cli = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-v"])).unwrap();
        assert!(cli.vcf_output_noop);
    }
}
