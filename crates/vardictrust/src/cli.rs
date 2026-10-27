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
            } else if let Some(v) = arg
                .strip_prefix("-th=")
                .or_else(|| arg.strip_prefix("--th="))
            {
                format!("--threads={v}")
            } else if let Some(v) = arg
                .strip_prefix("-th")
                .filter(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()))
            {
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
pub fn check_unsupported(args: &[String]) -> Result<(), vardictrust_core::types::VardictError> {
    for arg in args {
        let flag = arg.split_once('=').map_or(arg.as_str(), |(f, _)| f);
        if UNSUPPORTED_FLAGS.contains(&flag) {
            return Err(vardictrust_core::types::VardictError::Unsupported(format!(
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
    #[arg(short = 'f', long = "af-threshold", default_value = "0.01")]
    pub freq: String,

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

    /// Default printer destination: OUT or ERR (ports CmdParser.java:171-176, 505 defaultPrinter).
    /// Bare `-DP` or an unknown value normalize to OUT.
    #[arg(
        long = "DP",
        alias = "defaultPrinter",
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
    pub fn to_config(
        &self,
    ) -> Result<vardictrust_core::config::Config, vardictrust_core::types::VardictError> {
        let flag_filter =
            if self.flag_filter.starts_with("0x") || self.flag_filter.starts_with("0X") {
                u32::from_str_radix(&self.flag_filter[2..], 16).map_err(|e| {
                    vardictrust_core::types::VardictError::Parse(format!(
                        "Invalid -F value '{}': {}",
                        self.flag_filter, e
                    ))
                })?
            } else {
                self.flag_filter.parse().map_err(|e| {
                    vardictrust_core::types::VardictError::Parse(format!(
                        "Invalid -F value '{}': {}",
                        self.flag_filter, e
                    ))
                })?
            };

        if let Some(0) = self.tile_size {
            return Err(vardictrust_core::types::VardictError::Parse(
                "Invalid --tile-size: 0. Tile size must be >= 1".into(),
            ));
        }
        if let Some(n) = self.tile_size {
            if n > vardictrust_core::MAX_TILE_SIZE {
                return Err(vardictrust_core::types::VardictError::Parse(format!(
                    "Invalid --tile-size: {}. Must be <= {} to avoid u32 wrap",
                    n,
                    vardictrust_core::MAX_TILE_SIZE
                )));
            }
        }

        let parsed_freq: f64 = self.freq.parse().map_err(|e| {
            vardictrust_core::types::VardictError::Parse(format!(
                "Invalid -f value '{}': {}",
                self.freq, e
            ))
        })?;
        // In Perl (`var2vcf_valid.pl:21`: `my $Freq = $opt_f ? $opt_f : 0.02;`),
        // only `"0"` (or `""`) is falsy; non-`"0"` zero spellings such as `"0.0"`
        // or `"0.00"` are truthy strings in Perl boolean context and keep `0.0`.
        let af_threshold_truthy_zero = parsed_freq == 0.0 && self.freq.trim() != "0";

        let mut cfg = vardictrust_core::config::Config {
            exception_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            reference: self.reference.clone(),
            bam: self.bam.clone(),
            freq: parsed_freq,
            // The writer's $Freq is the literal -f value even under -p (perl
            // receives -f from the pipeline, not java's pileup override).
            af_threshold: parsed_freq,
            af_threshold_truthy_zero,
            sample: self
                .sample
                .clone()
                .unwrap_or_else(|| vardictrust_core::config::extract_sample_name(&self.bam, None)),
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
            min_reads: self.min_reads,
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
        };
        cfg.apply_pileup_overrides();
        Ok(cfg)
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
        assert_eq!(cli.to_config().unwrap().threads, all_cores());
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
        assert_eq!(cli.to_config().unwrap().threads, 1);
    }

    #[test]
    fn absent_threads_stays_all_cores() {
        // Deliberate divergence over Java's 1-thread default (task P1 report):
        // no -th at all → all cores, perf-only, no output effect.
        let cli = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa"])).unwrap();
        assert!(cli.threads.is_none());
        assert_eq!(cli.to_config().unwrap().threads, all_cores());
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
                cli.to_config().unwrap().debug_print.as_deref(),
                Some(expect),
                "argv={argv:?}"
            );
        }
    }

    #[test]
    fn invalid_flag_filter_returns_parse_error() {
        let cli = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-F", "bogus"])).unwrap();
        let err = cli.to_config().unwrap_err();
        assert!(matches!(
            err,
            vardictrust_core::types::VardictError::Parse(_)
        ));
        assert!(err.to_string().contains("Invalid -F value 'bogus'"));
    }

    #[test]
    fn tile_size_zero_returns_parse_error() {
        let cli =
            Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "--tile-size", "0"])).unwrap();
        let err = cli.to_config().unwrap_err();
        assert!(matches!(
            err,
            vardictrust_core::types::VardictError::Parse(_)
        ));
        assert!(err.to_string().contains("Invalid --tile-size: 0"));
    }

    #[test]
    fn tile_size_max_boundary_500mbp() {
        let ok_cli = Cli::try_parse_from(args(&[
            "-b",
            "a.bam",
            "-G",
            "r.fa",
            "--tile-size",
            "500000000",
        ]))
        .unwrap();
        assert_eq!(
            ok_cli.to_config().unwrap().tile_size,
            Some(vardictrust_core::MAX_TILE_SIZE)
        );

        let err_cli = Cli::try_parse_from(args(&[
            "-b",
            "a.bam",
            "-G",
            "r.fa",
            "--tile-size",
            "500000001",
        ]))
        .unwrap();
        let err = err_cli.to_config().unwrap_err();
        assert!(matches!(
            err,
            vardictrust_core::types::VardictError::Parse(_)
        ));
        assert!(err.to_string().contains("Must be <= 500000000"));
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

    #[test]
    fn freq_distinguishes_perl_falsy_zero_from_truthy_zero_decimal() {
        let cli_zero =
            Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-f", "0"])).unwrap();
        let cfg_zero = cli_zero.to_config().unwrap();
        assert_eq!(cfg_zero.af_threshold, 0.0);
        assert!(!cfg_zero.af_threshold_truthy_zero);

        let cli_zero_dec =
            Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-f", "0.0"])).unwrap();
        let cfg_zero_dec = cli_zero_dec.to_config().unwrap();
        assert_eq!(cfg_zero_dec.af_threshold, 0.0);
        assert!(cfg_zero_dec.af_threshold_truthy_zero);

        let cli_zero_exp =
            Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-f", "0e0"])).unwrap();
        let cfg_zero_exp = cli_zero_exp.to_config().unwrap();
        assert_eq!(cfg_zero_exp.af_threshold, 0.0);
        assert!(cfg_zero_exp.af_threshold_truthy_zero);
    }

    #[test]
    fn unsupported_flags_all_seven_return_clear_error() {
        for flag in ["-n", "-C", "-T", "--trim", "--VS", "-Z", "--downsample"] {
            let err = check_unsupported(&args(&[flag])).unwrap_err();
            assert!(
                matches!(err, vardictrust_core::types::VardictError::Unsupported(_)),
                "expected Unsupported error for {flag}, got {err:?}"
            );
            assert!(
                err.to_string()
                    .contains("not implemented in v1.0; use VarDictJava"),
                "unexpected error message for {flag}: {err}"
            );
            let eq_form = format!("{flag}=val");
            let err_eq =
                check_unsupported(&["vardictrust".to_string(), eq_form.clone()]).unwrap_err();
            assert!(
                matches!(
                    err_eq,
                    vardictrust_core::types::VardictError::Unsupported(_)
                ),
                "expected Unsupported error for {eq_form}, got {err_eq:?}"
            );
        }
    }

    #[test]
    fn normalize_args_all_five_thread_flag_forms() {
        let cases: &[(&[&str], usize)] = &[
            (&["-b", "a.bam", "-G", "r.fa", "-th", "4"], 4),
            (&["-b", "a.bam", "-G", "r.fa", "-th8"], 8),
            (&["-b", "a.bam", "-G", "r.fa", "--th", "4"], 4),
            (&["-b", "a.bam", "-G", "r.fa", "-th=4"], 4),
            (&["-b", "a.bam", "-G", "r.fa", "--th=4"], 4),
        ];
        for &(argv, expected_threads) in cases {
            let norm = normalize_args(&args(argv));
            let cli = Cli::try_parse_from(&norm)
                .unwrap_or_else(|e| panic!("failed to parse {argv:?} (norm={norm:?}): {e}"));
            assert_eq!(cli.threads, Some(expected_threads), "argv={argv:?}");
            assert_eq!(
                cli.to_config().unwrap().threads,
                expected_threads,
                "argv={argv:?}"
            );
        }
    }
}
