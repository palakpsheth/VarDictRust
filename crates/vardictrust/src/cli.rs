//! CLI argument parser — drop-in compatible with VarDictJava flags plus
//! ergonomic extras.
//!
//! ports Main.java / CmdParser.java

use clap::Parser;

#[path = "completions_man.rs"]
pub(crate) mod completions_man;

/// Shell target for `--completions <SHELL>`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, clap::ValueEnum)]
pub enum CompletionShell {
    Bash,
    Zsh,
    Fish,
}

/// VarDictJava v1.8.3 flags (CmdParser.java) with no v1.0 implementation.
/// Must be rejected loudly (§12) before Clap sees them.
const UNSUPPORTED_FLAGS: &[&str] = &[
    "-C",           // chromosome names are just numbers (deprecated upstream)
    "-Z",           // downsampling fraction (CmdParser.java:125-127, 433-440)
    "-downsample",  // single-dash long form of -Z
    "--downsample", // long form of -Z
];

/// GNU-only `--` long flags that have no VarDictJava single-dash equivalent and would
/// otherwise collide with short options or boolean clusters.
const GNU_ONLY_LONG_FLAGS_SINGLE_DASH: &[&str] = &[
    "-pass-only",
    "-paired-strict",
    "-vcf-pval",
    "-min-splitreads",
    "-no-pstd-filter",
    "-vcf-af",
    "-vcf-min-depth",
    "-vcf-min-vd",
    "-vcf-qmean",
    "-vcf-mapq",
    "-vcf-pmean",
    "-vcf-sn",
    "-vcf-max-nm",
    "-vcf-gt-freq",
    "-vcf-msi",
    "-vcf-cluster-dist",
    "-no-end",
    "-output",
    "-tbi",
    "-csi",
    "-index",
    "-stats-json",
    "-parity",
    "-compat",
    "-completions",
    "-manpage",
    "-reference",
    "-bam",
    "-af-threshold",
    "-name",
    "-region",
    "-threads",
    "-emit-intermediate",
    "-tile-size",
    "-version",
    "-help",
];

/// Value parser for signed integer CLI options (`-P`, `-x`, `-X`, `-J`, `-j`, `-T`).
/// Rejects `"-3"` because Commons-CLI treats `-3` as the short option flag `moveIndelsTo3`,
/// leaving the preceding option missing its argument (`MissingArgumentException`).
fn parse_java_signed_i32(s: &str) -> Result<i32, String> {
    if s == "-3" {
        return Err("Missing argument for option (-3 is a flag)".to_string());
    }
    vardictrust_core::config::parse_java_number_i32(s)
}

/// Strict parser for `--vcf-pval <FLOAT>` (`var2vcf_paired.pl -P`).
/// Preserves the normalized raw CLI token for verbatim `P$PVAL` header/FILTER formatting
/// while canonicalizing signed zeros (`-0.0` / `+0.0` -> `(0.0, "0.0")`).
fn parse_vcf_pval_token(s: &str) -> Result<(f64, String), String> {
    if s == "-3" {
        return Err(format!("Invalid --vcf-pval value '{s}'"));
    }
    vardictrust_core::config::parse_vcf_pval_token(s)
}

macro_rules! def_vcf_cli_parsers {
    ($(($fn_name:ident, $flag:literal, $ty:ty, $delegate:path)),* $(,)?) => {
        $(
            fn $fn_name(s: &str) -> Result<$ty, String> {
                $delegate(s, $flag)
            }
        )*
    };
}

def_vcf_cli_parsers! {
    (parse_min_splitreads, "--min-splitreads", i32, vardictrust_core::config::parse_vcf_nonneg_i32),
    (parse_vcf_af, "--vcf-af", vardictrust_core::config::VcfThreshold<f64>, vardictrust_core::config::parse_vcf_threshold_nonneg_f64),
    (parse_vcf_min_depth, "--vcf-min-depth", vardictrust_core::config::VcfThreshold<i32>, vardictrust_core::config::parse_vcf_threshold_nonneg_i32),
    (parse_vcf_min_vd, "--vcf-min-vd", vardictrust_core::config::VcfThreshold<i32>, vardictrust_core::config::parse_vcf_threshold_nonneg_i32),
    (parse_vcf_qmean, "--vcf-qmean", vardictrust_core::config::VcfThreshold<f64>, vardictrust_core::config::parse_vcf_threshold_nonneg_f64),
    (parse_vcf_mapq, "--vcf-mapq", vardictrust_core::config::VcfThreshold<f64>, vardictrust_core::config::parse_vcf_threshold_nonneg_f64),
    (parse_vcf_pmean, "--vcf-pmean", vardictrust_core::config::VcfThreshold<f64>, vardictrust_core::config::parse_vcf_threshold_nonneg_f64),
    (parse_vcf_sn, "--vcf-sn", vardictrust_core::config::VcfThreshold<f64>, vardictrust_core::config::parse_vcf_threshold_nonneg_f64),
    (parse_vcf_max_nm, "--vcf-max-nm", vardictrust_core::config::VcfThreshold<f64>, vardictrust_core::config::parse_vcf_threshold_nonneg_f64),
    (parse_vcf_gt_freq, "--vcf-gt-freq", f64, vardictrust_core::config::parse_vcf_unit_f64),
    (parse_vcf_msi, "--vcf-msi", vardictrust_core::config::VcfThreshold<f64>, vardictrust_core::config::parse_vcf_threshold_nonneg_f64),
    (parse_vcf_cluster_dist, "--vcf-cluster-dist", vardictrust_core::config::VcfThreshold<i32>, vardictrust_core::config::parse_vcf_threshold_nonneg_i32),
}

fn is_known_single_dash_long_flag(flag: &str) -> bool {
    matches!(
        flag,
        "-th"
            | "-UN"
            | "-fisher"
            | "-chimeric"
            | "-deldupvar"
            | "-nosv"
            | "-splice"
            | "-header"
            | "-dedup"
            | "-debug"
            | "-verbose"
            | "-amplicon"
            | "-DP"
            | "-default-printer"
            | "-defaultPrinter"
            | "-adaptor"
            | "-mfreq"
            | "-nmfreq"
            | "-insert-std"
            | "-insert-size"
            | "-ref-extension"
            | "-crispr"
            | "-trim"
            | "-VS"
            | "-downsample"
    ) || flag
        .strip_prefix("-th")
        .is_some_and(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()))
}

fn is_bool_short_flag(c: char) -> bool {
    matches!(
        c,
        'H' | '?' | 'h' | 'v' | 'i' | 'p' | 'D' | 't' | '3' | 'K' | 'u' | 'U' | 'y'
    )
}

/// Inspects a single-dash short option or short-option cluster (`-<chars>`), skipping known
/// multi-char single-dash options (`-nosv`, `-nmfreq`, `-adaptor`, etc.).
/// Returns `(unsupported_short_char, consumes_mandatory_next_token)`.
fn classify_short_cluster(arg: &str) -> (Option<char>, bool) {
    let (head, has_eq) = match arg.split_once('=') {
        Some((h, _)) => (h, true),
        None => (arg, false),
    };
    let Some(rest) = head.strip_prefix('-') else {
        return (None, false);
    };
    if rest.is_empty()
        || rest.starts_with('-')
        || is_known_single_dash_long_flag(head)
        || GNU_ONLY_LONG_FLAGS_SINGLE_DASH.contains(&head)
    {
        return (None, false);
    }
    for (idx, ch) in rest.char_indices() {
        let rem_empty = rest.get(idx + ch.len_utf8()..).is_none_or(str::is_empty);
        let consumes_next = !has_eq && rem_empty;
        match ch {
            'C' => return (Some('C'), false),
            'Z' => return (Some('Z'), consumes_next),
            'G' | 'b' | 'f' | 'N' | 'R' | 'c' | 'S' | 'E' | 's' | 'e' | 'g' | 'r' | 'q' | 'Q'
            | 'P' | 'm' | 'M' | 'V' | 'L' | 'x' | 'X' | 'B' | 'o' | 'O' | 'I' | 'F' | 'w' | 'W'
            | 'A' | 'Y' | 'J' | 'j' | 'd' | 'n' | 'T' => {
                return (None, consumes_next);
            }
            'z' | 'k' | 'a' => return (None, false),
            c if is_bool_short_flag(c) => {}
            _ => return (None, false),
        }
    }
    (None, false)
}

/// Returns `true` when `arg` is an option token (without `=`) that unconditionally consumes
/// `args[i + 1]` as its mandatory argument.
fn takes_mandatory_next_token(arg: &str) -> bool {
    if arg.contains('=') {
        return false;
    }
    if matches!(
        arg,
        "--reference"
            | "--bam"
            | "--af-threshold"
            | "--name"
            | "--region"
            | "--tile-size"
            | "-insert-size"
            | "--insert-size"
            | "-insert-std"
            | "--insert-std"
            | "-ref-extension"
            | "--ref-extension"
            | "-mfreq"
            | "--mfreq"
            | "-nmfreq"
            | "--nmfreq"
            | "-crispr"
            | "--crispr"
            | "-adaptor"
            | "--adaptor"
            | "--output"
            | "--compat"
            | "--completions"
            | "--stats-json"
            | "--vcf-pval"
            | "--min-splitreads"
            | "--vcf-af"
            | "--vcf-min-depth"
            | "--vcf-min-vd"
            | "--vcf-qmean"
            | "--vcf-mapq"
            | "--vcf-pmean"
            | "--vcf-sn"
            | "--vcf-max-nm"
            | "--vcf-gt-freq"
            | "--vcf-msi"
            | "--vcf-cluster-dist"
            | "-trim"
            | "--trim"
            | "-VS"
            | "--VS"
            | "-downsample"
            | "--downsample"
    ) {
        return true;
    }
    classify_short_cluster(arg).1
}

fn arg_has_region_short(arg: &str) -> bool {
    let head = arg.split_once('=').map_or(arg, |(h, _)| h);
    if arg.starts_with("--")
        || is_known_single_dash_long_flag(head)
        || GNU_ONLY_LONG_FLAGS_SINGLE_DASH.contains(&head)
    {
        return false;
    }
    let Some(rest) = head.strip_prefix('-') else {
        return false;
    };
    for ch in rest.chars() {
        if ch == 'R' {
            return true;
        }
        if !is_bool_short_flag(ch) {
            return false;
        }
    }
    false
}

fn scan_has_region(args: &[String]) -> bool {
    let mut i = 1;
    while let Some(arg) = args.get(i) {
        if arg == "--" {
            break;
        }
        if arg == "--region" || arg.starts_with("--region=") || arg_has_region_short(arg) {
            return true;
        }
        i += if takes_mandatory_next_token(arg) {
            2
        } else {
            1
        };
    }
    false
}

/// Splits a boolean-prefixed short option cluster whose non-boolean tail begins with `z`, `k`,
/// `a`, or `R` (e.g. `-pz`, `-pz0`, `-uz`, `-uk`, `-ua`, `-pR`) into `(bool_prefix, tail)`.
fn split_bool_prefix_cluster(arg: &str) -> Option<(&str, &str)> {
    let head = arg.split_once('=').map_or(arg, |(h, _)| h);
    if arg.starts_with("--")
        || is_known_single_dash_long_flag(head)
        || GNU_ONLY_LONG_FLAGS_SINGLE_DASH.contains(&head)
    {
        return None;
    }
    let rest = arg.strip_prefix('-')?;
    let split_idx = rest.char_indices().find_map(|(idx, c)| {
        if is_bool_short_flag(c) {
            None
        } else {
            Some(idx)
        }
    })?;
    if split_idx == 0 {
        return None;
    }
    let bool_prefix = rest.get(..split_idx)?;
    let tail = rest.get(split_idx..)?;
    let mut tail_chars = tail.chars();
    let first = tail_chars.next()?;
    let tail_rest = tail_chars.as_str();
    let should_split = match first {
        'R' => true,
        'z' | 'k' => {
            tail_rest.is_empty()
                || tail_rest.starts_with('=')
                || parse_java_signed_i32(tail_rest).is_ok()
        }
        'a' => {
            tail_rest.is_empty()
                || tail_rest.starts_with('=')
                || tail_rest.contains(':')
                || tail_rest
                    .strip_prefix(['+', '-'])
                    .unwrap_or(tail_rest)
                    .as_bytes()
                    .first()
                    .is_some_and(u8::is_ascii_digit)
        }
        _ => false,
    };
    if should_split {
        Some((bool_prefix, tail))
    } else {
        None
    }
}

fn normalize_optional_or_passthrough(
    effective_arg: &str,
    args: &[String],
    i: usize,
    has_region: bool,
    out: &mut Vec<String>,
) -> usize {
    let head = effective_arg
        .split_once('=')
        .map_or(effective_arg, |(h, _)| h);
    if GNU_ONLY_LONG_FLAGS_SINGLE_DASH.contains(&head) {
        out.push(effective_arg.to_string());
        return i + 1;
    }
    if effective_arg == "-z" {
        let next_is_num = args.get(i + 1).is_some_and(|n| {
            parse_java_signed_i32(n).is_ok() && (has_region || i + 2 < args.len())
        });
        if next_is_num {
            out.push(format!("-z={}", args[i + 1]));
            return i + 2;
        }
        out.push("-z=1".to_string());
        return i + 1;
    } else if let Some(v) = effective_arg
        .strip_prefix("-z")
        .filter(|v| !v.is_empty() && !v.starts_with('=') && parse_java_signed_i32(v).is_ok())
    {
        out.push(format!("-z={v}"));
        return i + 1;
    } else if effective_arg == "-k" {
        if let Some(next) = args.get(i + 1) {
            if parse_java_signed_i32(next).is_ok()
                || (!next.starts_with('-') && (has_region || i + 2 < args.len()))
            {
                out.push(format!("-k={next}"));
                return i + 2;
            }
        }
        out.push("-k=1".to_string());
        return i + 1;
    } else if let Some(v) = effective_arg
        .strip_prefix("-k")
        .filter(|v| !v.is_empty() && !v.starts_with('=') && parse_java_signed_i32(v).is_ok())
    {
        out.push(format!("-k={v}"));
        return i + 1;
    } else if let Some(v) = effective_arg
        .strip_prefix("-a")
        .filter(|v| !v.is_empty() && !v.starts_with('=') && v.contains(':'))
    {
        out.push(format!("-a={v}"));
        return i + 1;
    } else if effective_arg == "-a" || effective_arg == "-amplicon" || effective_arg == "--amplicon"
    {
        match args.get(i + 1) {
            None => {
                out.push("-a=10:0.95".to_string());
                return i + 1;
            }
            Some(next) if next.starts_with('-') => {
                out.push("-a=10:0.95".to_string());
                return i + 1;
            }
            Some(next) if next.is_empty() => {
                out.push("-a=10:0.95".to_string());
                return i + 2;
            }
            Some(next) if next.contains(':') || has_region || i + 2 < args.len() => {
                out.push(format!("-a={next}"));
                return i + 2;
            }
            Some(_) => {
                out.push("-a=10:0.95".to_string());
                return i + 1;
            }
        }
    }

    out.push(effective_arg.to_string());
    if takes_mandatory_next_token(effective_arg) {
        if let Some(next_val) = args.get(i + 1) {
            out.push(next_val.clone());
            return i + 2;
        }
    }
    i + 1
}

/// Rewrite Java-style multi-char single-dash options and optional-arg flags (`-th`, `-UN`, `-z`, `-k`, `-a`)
/// to their clap-compatible forms while preserving mandatory option values and `--` literals verbatim.
pub fn normalize_args(args: &[String]) -> Vec<String> {
    let has_region = scan_has_region(args);
    let mut out = Vec::with_capacity(args.len());
    if let Some(first) = args.first() {
        out.push(first.clone());
    }
    let mut seen_threads = false;
    let mut i = 1;
    while let Some(arg) = args.get(i) {
        if arg == "--" {
            if let Some(tail) = args.get(i..) {
                out.extend_from_slice(tail);
            }
            break;
        }
        if arg == "-th" || arg == "--th" || arg == "--threads" {
            let next_is_val = args.get(i + 1).is_some_and(|n| {
                !n.starts_with('-') && vardictrust_core::config::parse_java_number_usize(n).is_ok()
            });
            if !seen_threads {
                seen_threads = true;
                out.push("--threads".to_string());
                if next_is_val {
                    out.push(args[i + 1].clone());
                }
            }
            i += if next_is_val { 2 } else { 1 };
            continue;
        } else if let Some(v) = arg
            .strip_prefix("-th=")
            .or_else(|| arg.strip_prefix("--th="))
            .or_else(|| arg.strip_prefix("--threads="))
        {
            if !seen_threads {
                seen_threads = true;
                out.push(format!("--threads={v}"));
            }
            i += 1;
            continue;
        } else if let Some(v) = arg
            .strip_prefix("-th")
            .filter(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()))
        {
            if !seen_threads {
                seen_threads = true;
                out.push(format!("--threads={v}"));
            }
            i += 1;
            continue;
        } else if arg == "--?" {
            out.push("-?".to_string());
            i += 1;
            continue;
        } else if matches!(
            arg.as_str(),
            "-UN"
                | "-fisher"
                | "-chimeric"
                | "-deldupvar"
                | "-nosv"
                | "-splice"
                | "-header"
                | "-dedup"
                | "-debug"
                | "-verbose"
        ) {
            out.push(format!("-{arg}"));
            i += 1;
            continue;
        } else if matches!(
            arg.as_str(),
            "-DP"
                | "--DP"
                | "-default-printer"
                | "--default-printer"
                | "-defaultPrinter"
                | "--defaultPrinter"
        ) {
            if let Some(next) = args.get(i + 1).filter(|n| {
                !n.starts_with('-')
                    && (matches!(n.as_str(), "OUT" | "ERR") || has_region || i + 2 < args.len())
            }) {
                out.push(format!("--DP={next}"));
                i += 2;
            } else {
                out.push("--DP=OUT".to_string());
                i += 1;
            }
            continue;
        } else if let Some((flag, val)) = arg.split_once('=').filter(|(f, _)| {
            matches!(
                *f,
                "-DP"
                    | "-default-printer"
                    | "-defaultPrinter"
                    | "-adaptor"
                    | "-mfreq"
                    | "-nmfreq"
                    | "-insert-std"
                    | "-insert-size"
                    | "-ref-extension"
                    | "-crispr"
                    | "-trim"
                    | "-VS"
            )
        }) {
            out.push(format!("-{flag}={val}"));
            i += 1;
            continue;
        } else if let Some(val) = arg.strip_prefix("-amplicon=") {
            out.push(format!("-a={val}"));
            i += 1;
            continue;
        } else if matches!(
            arg.as_str(),
            "-adaptor"
                | "-mfreq"
                | "-nmfreq"
                | "-insert-std"
                | "-insert-size"
                | "-ref-extension"
                | "-crispr"
                | "-trim"
                | "-VS"
        ) {
            out.push(format!("-{arg}"));
            if let Some(next_val) = args.get(i + 1) {
                out.push(next_val.clone());
                i += 2;
            } else {
                i += 1;
            }
            continue;
        }

        if let Some((bool_prefix, tail)) = split_bool_prefix_cluster(arg) {
            out.push(format!("-{bool_prefix}"));
            let tail_flag = format!("-{tail}");
            i = normalize_optional_or_passthrough(&tail_flag, args, i, has_region, &mut out);
            continue;
        }
        i = normalize_optional_or_passthrough(arg, args, i, has_region, &mut out);
    }
    out
}

/// Fail loudly if any known-unimplemented VarDictJava flag was supplied.
/// Runs before Clap parsing, which would otherwise emit a generic
/// "unexpected argument" error instead of the actionable §12 message.
pub fn check_unsupported(args: &[String]) -> Result<(), vardictrust_core::types::VardictError> {
    let mut i = 1;
    while let Some(arg) = args.get(i) {
        if arg == "--" {
            break;
        }
        let head = arg.split_once('=').map_or(arg.as_str(), |(f, _)| f);
        if UNSUPPORTED_FLAGS.contains(&head) {
            return Err(vardictrust_core::types::VardictError::Unsupported(format!(
                "{head} not implemented in v1.0; use VarDictJava"
            )));
        }
        if GNU_ONLY_LONG_FLAGS_SINGLE_DASH.contains(&head) {
            return Err(vardictrust_core::types::VardictError::Parse(format!(
                "Unknown single-dash option '{head}'; did you mean '-{head}'?"
            )));
        }
        let (unsupported_short, consumes_next) = classify_short_cluster(arg);
        if let Some(ch) = unsupported_short {
            return Err(vardictrust_core::types::VardictError::Unsupported(format!(
                "-{ch} not implemented in v1.0; use VarDictJava"
            )));
        }
        i += if consumes_next || takes_mandatory_next_token(arg) {
            2
        } else {
            1
        };
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
        "License: MIT AND GPL-2.0-or-later (see LICENSE-3RD-PARTY)"
    ),
    about,
    disable_version_flag = true,
    // `-h` belongs to Java's column header (CmdParser.java:223, 67), so clap
    // help moves to `--help`/`-H` (CmdParser.java:222).
    disable_help_flag = true,
    args_override_self = true
)]
pub struct Cli {
    // ── 1. Core Inputs & Outputs ─────────────────────────────────────────────
    /// Print version.
    #[arg(
        long = "version",
        action = clap::ArgAction::Version,
        help_heading = "Core Inputs & Outputs"
    )]
    _version: (),

    /// Print help.
    #[arg(
        short = 'H',
        short_alias = '?',
        long = "help",
        action = clap::ArgAction::Help,
        help_heading = "Core Inputs & Outputs"
    )]
    _help: (),

    /// Print shell completion script to stdout and exit (`bash`, `zsh`, or `fish`).
    #[arg(
        long = "completions",
        value_name = "SHELL",
        num_args = 1,
        conflicts_with = "manpage",
        help_heading = "Core Inputs & Outputs"
    )]
    pub completions: Vec<CompletionShell>,

    /// Print standalone `roff` `vardictrust(1)` manpage to stdout and exit.
    #[arg(long = "manpage", help_heading = "Core Inputs & Outputs")]
    pub manpage: bool,

    /// Reference FASTA (.fai required). Aliases: --reference.
    #[arg(
        short = 'G',
        long = "reference",
        required_unless_present_any = ["completions", "manpage"],
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub reference: Vec<String>,

    /// Indexed BAM/CRAM. Somatic: "tumor.bam|normal.bam". Alias: --bam.
    #[arg(
        short = 'b',
        long = "bam",
        required_unless_present_any = ["completions", "manpage"],
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub bam: Vec<String>,

    /// Sample name (somatic: "tumor|normal"). Alias: --name.
    #[arg(
        short = 'N',
        long = "name",
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub sample: Vec<String>,

    /// Regular expression to extract sample name from BAM filenames (ports CmdParser.java:86-93).
    #[arg(
        short = 'n',
        value_name = "REGEXP",
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub sample_regexp: Vec<String>,

    /// Region chr:start-end[:gene] (instead of BED).
    #[arg(
        short = 'R',
        long = "region",
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub region: Vec<String>,

    /// Column header output: print the .var column names before results
    /// (ports CmdParser.java:67 `printHeader` + SimpleMode.java:114-127).
    #[arg(short = 'h', long = "header", help_heading = "Core Inputs & Outputs")]
    pub column_header: bool,

    /// Write output (.vcf or .var) directly to file; .gz/.bgz suffixes use multithreaded BGZF
    #[arg(
        long = "output",
        value_name = "PATH",
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub output: Vec<String>,

    /// Write a Tabix (.tbi) index alongside a BGZF-compressed --output .vcf.gz file.
    #[arg(
        long = "tbi",
        visible_alias = "index",
        help_heading = "Core Inputs & Outputs"
    )]
    pub tbi: bool,

    /// Write a CSI (.csi) index alongside a BGZF-compressed --output .vcf.gz file.
    #[arg(long = "csi", help_heading = "Core Inputs & Outputs")]
    pub index_csi: bool,

    /// Write run telemetry (mode, threads, regions, tiles, elapsed_ms, peak_rss_kb) as JSON to PATH
    #[arg(
        long = "stats-json",
        value_name = "PATH",
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub stats_json: Vec<String>,

    /// Explicitly lock 100% VarDictJava v1.8.3 byte-for-byte parity mode (default in v1.x).
    #[arg(long = "parity", help_heading = "Core Inputs & Outputs")]
    pub parity: bool,

    /// Pin compatibility target (supported in v1.0: 'vardictjava-1.8.3' or '1.8.3').
    #[arg(
        long = "compat",
        value_name = "TARGET",
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub compat: Option<String>,

    /// Emit the intermediate .var table instead of VCF.
    #[arg(long = "emit-intermediate", help_heading = "Core Inputs & Outputs")]
    pub emit_intermediate: bool,

    /// The delimiter for splitting region info. [default: \t]
    #[arg(
        short = 'd',
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub delimiter: Vec<String>,

    /// VCF format output. Defined but never read by VarDictJava v1.8.3
    /// (CmdParser.java:224 vs parseCmd 59-186); accepted and ignored for
    /// drop-in compatibility.
    #[arg(short = 'v', help_heading = "Core Inputs & Outputs")]
    pub vcf_output_noop: bool,

    /// Verbose mode.
    #[arg(short = 'y', long = "verbose", help_heading = "Core Inputs & Outputs")]
    pub verbose: bool,

    /// Debug mode: append per-read debug columns to the intermediate output.
    #[arg(short = 'D', long = "debug", help_heading = "Core Inputs & Outputs")]
    pub debug: bool,

    /// Default printer destination: OUT or ERR (ports CmdParser.java:171-176, 505 defaultPrinter).
    /// Bare `-DP` or an unknown value normalize to OUT.
    #[arg(
        long = "DP",
        visible_alias = "default-printer",
        alias = "defaultPrinter",
        num_args(0..=1),
        default_missing_value = "OUT",
        require_equals = true,
        help_heading = "Core Inputs & Outputs"
    )]
    pub debug_print: Vec<String>,

    // ── 2. VCF Post-Filtering (var2vcf) ──────────────────────────────────────
    /// Omit the END tag in VCF INFO output (ports var2vcf_valid.pl -E).
    #[arg(long = "no-end", help_heading = "VCF Post-Filtering (var2vcf)")]
    pub no_end: bool,

    /// Emit only PASS variants in VCF output (ports var2vcf_valid.pl / var2vcf_paired.pl -S).
    #[arg(long = "pass-only", help_heading = "VCF Post-Filtering (var2vcf)")]
    pub pass_only: bool,

    /// Disable germline rescue in paired somatic VCF filtering (ports var2vcf_paired.pl -M).
    #[arg(long = "paired-strict", help_heading = "VCF Post-Filtering (var2vcf)")]
    pub paired_strict: bool,

    /// Maximum p-value threshold for high-confidence somatic calls in VCF output (ports var2vcf_paired.pl -P). [default: 0.05]
    #[arg(
        long = "vcf-pval",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_pval_token,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_pval: Vec<(f64, String)>,

    /// Minimum split reads required for structural variants in single-sample VCF output (ports var2vcf_valid.pl -T). [default: 1]
    #[arg(
        long = "min-splitreads",
        value_name = "INT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_min_splitreads,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub min_splitreads: Vec<i32>,

    /// Disable the PSTD=0 single-position read filter in single-sample VCF output (ports var2vcf_valid.pl -P 0).
    #[arg(long = "no-pstd-filter", help_heading = "VCF Post-Filtering (var2vcf)")]
    pub no_pstd_filter: bool,

    /// Minimum allele frequency threshold for VCF FILTER/GT (overrides -f in VCF output; ports var2vcf_valid.pl / var2vcf_paired.pl -f).
    #[arg(
        long = "vcf-af",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_af,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_af: Vec<vardictrust_core::config::VcfThreshold<f64>>,

    /// Minimum total read depth for VCF FILTER (ports var2vcf_valid.pl -d [default: 3] / var2vcf_paired.pl -d [default: 5]).
    #[arg(
        long = "vcf-min-depth",
        value_name = "INT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_min_depth,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_min_depth: Vec<vardictrust_core::config::VcfThreshold<i32>>,

    /// Minimum variant depth for VCF FILTER (ports var2vcf_valid.pl -v [default: 2] / var2vcf_paired.pl -v [default: 3]).
    #[arg(
        long = "vcf-min-vd",
        value_name = "INT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_min_vd,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_min_vd: Vec<vardictrust_core::config::VcfThreshold<i32>>,

    /// Minimum mean base quality for VCF FILTER (ports var2vcf_valid.pl / var2vcf_paired.pl -q). [default: 22.5]
    #[arg(
        long = "vcf-qmean",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_qmean,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_qmean: Vec<vardictrust_core::config::VcfThreshold<f64>>,

    /// Minimum mean mapping quality for VCF FILTER (ports var2vcf_valid.pl -Q [default: 10] / var2vcf_paired.pl -Q [default: 0]).
    #[arg(
        long = "vcf-mapq",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_mapq,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_mapq: Vec<vardictrust_core::config::VcfThreshold<f64>>,

    /// Minimum mean read position for VCF FILTER (ports var2vcf_valid.pl / var2vcf_paired.pl -p). [default: 8]
    #[arg(
        long = "vcf-pmean",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_pmean,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_pmean: Vec<vardictrust_core::config::VcfThreshold<f64>>,

    /// Minimum signal-to-noise (qratio) for VCF FILTER (ports var2vcf_valid.pl / var2vcf_paired.pl -o). [default: 1.5]
    #[arg(
        long = "vcf-sn",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_sn,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_sn: Vec<vardictrust_core::config::VcfThreshold<f64>>,

    /// Maximum mean mismatches (NM) for VCF FILTER (ports var2vcf_valid.pl / var2vcf_paired.pl -m). [default: 5.25]
    #[arg(
        long = "vcf-max-nm",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_max_nm,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_max_nm: Vec<vardictrust_core::config::VcfThreshold<f64>>,

    /// Homozygous 1/1 genotype AF margin (1 - AF < FLOAT) in VCF output (ports var2vcf_valid.pl / var2vcf_paired.pl -F). [default: 0.2]
    #[arg(
        long = "vcf-gt-freq",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_gt_freq,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_gt_freq: Vec<f64>,

    /// Maximum non-monomer MSI repeat count for VCF FILTER (ports var2vcf_valid.pl / var2vcf_paired.pl -I). [default: 12]
    #[arg(
        long = "vcf-msi",
        value_name = "FLOAT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_msi,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_msi: Vec<vardictrust_core::config::VcfThreshold<f64>>,

    /// SNV cluster filter window in bp for VCF output (ports var2vcf_valid.pl / var2vcf_paired.pl -c). [default: 0]
    #[arg(
        long = "vcf-cluster-dist",
        value_name = "INT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_vcf_cluster_dist,
        help_heading = "VCF Post-Filtering (var2vcf)"
    )]
    pub vcf_cluster_dist: Vec<vardictrust_core::config::VcfThreshold<i32>>,

    // ── 3. Coordinate & BED Column Mapping ───────────────────────────────────
    /// BED column for chromosome (1-based).
    #[arg(
        short = 'c',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_usize,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub col_chr: Vec<usize>,

    /// BED column for start (1-based).
    #[arg(
        short = 'S',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_usize,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub col_start: Vec<usize>,

    /// BED column for end (1-based).
    #[arg(
        short = 'E',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_usize,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub col_end: Vec<usize>,

    /// BED column for segment starts (1-based).
    #[arg(
        short = 's',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_usize,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub col_seg_start: Vec<usize>,

    /// BED column for segment ends (1-based).
    #[arg(
        short = 'e',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_usize,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub col_seg_end: Vec<usize>,

    /// BED column for gene (1-based).
    #[arg(
        short = 'g',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_usize,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub col_gene: Vec<usize>,

    /// BED coordinates are 0-based (0: 1-based, 1/empty: 0-based; CmdParser.java:101-112).
    #[arg(
        short = 'z',
        num_args(0..=1),
        default_missing_value = "1",
        require_equals = true,
        allow_negative_numbers = true,
        value_parser = vardictrust_core::config::parse_java_number_i32,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub zero_based: Vec<i32>,

    /// Extension for BED region.
    #[arg(
        short = 'x',
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_java_signed_i32,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub extension: Vec<i32>,

    /// Number of nucleotide extension. [default: 2]
    #[arg(
        short = 'X',
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_java_signed_i32,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub extension_nuc: Vec<i32>,

    /// Max region length. [default: 1200]
    #[arg(
        short = 'Y',
        long = "ref-extension",
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Coordinate & BED Column Mapping"
    )]
    pub max_region_len: Vec<u32>,

    // ── 4. Variant Calling & Quality Thresholds ──────────────────────────────
    /// Allele-frequency threshold. Alias: --af-threshold. [default: 0.01]
    #[arg(
        short = 'f',
        long = "af-threshold",
        num_args = 1,
        allow_negative_numbers = true,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub freq: Vec<String>,

    /// Min variant reads. [default: 2]
    #[arg(
        short = 'r',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub min_reads: Vec<u32>,

    /// Min base quality. [default: 22.5]
    #[arg(
        short = 'q',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_f64,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub base_qual: Vec<f64>,

    /// Min mapping quality.
    #[arg(
        short = 'Q',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub map_qual: Vec<u32>,

    /// Mean read position filter. [default: 5]
    #[arg(
        short = 'P',
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_java_signed_i32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub read_pos: Vec<i32>,

    /// Max mismatches. [default: 8]
    #[arg(
        short = 'm',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub max_mismatches: Vec<u32>,

    /// Minimum match+insertion read length; shorter reads ignored. [default: 0]
    #[arg(
        short = 'M',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub minmatch: Vec<u32>,

    /// Min reads for bias. [default: 2]
    #[arg(
        short = 'B',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub min_reads_bias: Vec<u32>,

    /// qratio filter. [default: 1.5]
    #[arg(
        short = 'o',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_f64,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub qratio: Vec<f64>,

    /// MapQ for output. [default: 0]
    #[arg(
        short = 'O',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_f64,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub out_mapq: Vec<f64>,

    /// SAM flag filter. [default: 0x504]
    #[arg(
        short = 'F',
        num_args = 1,
        allow_negative_numbers = true,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub flag_filter: Vec<String>,

    /// Local realignment 0/1. [default: 1]
    #[arg(
        short = 'k',
        num_args(0..=1),
        default_missing_value = "1",
        require_equals = true,
        allow_negative_numbers = true,
        value_parser = vardictrust_core::config::parse_java_number_i32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub realign: Vec<i32>,

    /// Shift indel positions to their 3'-most alternative alignment
    /// (ports CmdParser.java:71 moveIndelsTo3; shift3 applied in tovars.rs).
    #[arg(short = '3', help_heading = "Variant Calling & Quality Thresholds")]
    pub three_prime: bool,

    /// Include Ns in depth calculation.
    #[arg(short = 'K', help_heading = "Variant Calling & Quality Thresholds")]
    pub include_n: bool,

    /// Pileup-only mode.
    #[arg(short = 'p', help_heading = "Variant Calling & Quality Thresholds")]
    pub pileup: bool,

    /// Dedup: remove duplicated reads, keeping one pair per start position
    /// (ports CmdParser.java:70, 229).
    #[arg(
        short = 't',
        long = "dedup",
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub dedup: bool,

    /// Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using forward read only (-u).
    #[arg(short = 'u', help_heading = "Variant Calling & Quality Thresholds")]
    pub unique_mode_alignment_enabled: bool,

    /// Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using first read only (-UN).
    #[arg(long = "UN", help_heading = "Variant Calling & Quality Thresholds")]
    pub unique_mode_second_in_pair_enabled: bool,

    /// Trim bases after [INT] bases in the reads (ports CmdParser.java:122, 411-417). [default: 0]
    #[arg(
        short = 'T',
        long = "trim",
        value_name = "INT",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_java_signed_i32,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub trim_bases_after: Vec<i32>,

    /// SAM/BAM validation stringency: SILENT, LENIENT, or STRICT (ports CmdParser.java:146-148, 491-500).
    #[arg(
        long = "VS",
        value_name = "STRINGENCY",
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub validation_stringency: Vec<String>,

    /// Filter adaptor sequences so that they are not used in realignment. Multiple adaptors can be supplied by setting them with comma.
    #[arg(
        long = "adaptor",
        num_args = 1,
        allow_hyphen_values = true,
        help_heading = "Variant Calling & Quality Thresholds"
    )]
    pub adaptor: Vec<String>,

    // ── 5. Modality Flags: Somatic & Amplicon ────────────────────────────────
    /// Min normal AF for somatic. [default: 0.05]
    #[arg(
        short = 'V',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_f64,
        help_heading = "Modality Flags: Somatic & Amplicon"
    )]
    pub min_normal_af: Vec<f64>,

    /// Use the in-engine Fisher test.
    #[arg(long = "fisher", help_heading = "Modality Flags: Somatic & Amplicon")]
    pub fisher: bool,

    /// Amplicon mode "int:float". [default: 10:0.95 when 8-col BED]
    #[arg(
        short = 'a',
        long = "amplicon",
        num_args(0..=1),
        default_missing_value = "10:0.95",
        require_equals = true,
        help_heading = "Modality Flags: Somatic & Amplicon"
    )]
    pub amplicon: Vec<String>,

    /// Splice mode.
    #[arg(
        short = 'i',
        long = "splice",
        help_heading = "Modality Flags: Somatic & Amplicon"
    )]
    pub splice: bool,

    /// Mono frequency threshold. [default: 0.25]
    #[arg(
        long = "mfreq",
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_f64,
        help_heading = "Modality Flags: Somatic & Amplicon"
    )]
    pub mfreq: Vec<f64>,

    /// Normal mono frequency threshold. [default: 0.1]
    #[arg(
        long = "nmfreq",
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_f64,
        help_heading = "Modality Flags: Somatic & Amplicon"
    )]
    pub nmfreq: Vec<f64>,

    /// CRISPR cutting site.
    #[arg(
        short = 'J',
        long = "crispr",
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_java_signed_i32,
        help_heading = "Modality Flags: Somatic & Amplicon"
    )]
    pub crispr_cutting_site: Vec<i32>,

    /// CRISPR flanking size.
    #[arg(
        short = 'j',
        num_args = 1,
        allow_negative_numbers = true,
        value_parser = parse_java_signed_i32,
        help_heading = "Modality Flags: Somatic & Amplicon"
    )]
    pub crispr_flanking_size: Vec<i32>,

    // ── 6. Structural Variants & Complex Indels ──────────────────────────────
    /// Disable SV calling.
    #[arg(
        short = 'U',
        long = "nosv",
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub no_sv: bool,

    /// SV minimum length. [default: 1000]
    #[arg(
        short = 'L',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub sv_min_len: Vec<u32>,

    /// SV read length. [default: 300]
    #[arg(
        short = 'w',
        long = "insert-size",
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub sv_read_len: Vec<u32>,

    /// SV flanking. [default: 100]
    #[arg(
        short = 'W',
        long = "insert-std",
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub sv_flanking: Vec<u32>,

    /// SV minimum AF. [default: 4]
    #[arg(
        short = 'A',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub sv_min_af: Vec<u32>,

    /// Indel size. [default: 50]
    #[arg(
        short = 'I',
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub indel_size: Vec<u32>,

    /// Chimeric read handling.
    #[arg(
        long = "chimeric",
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub chimeric: bool,

    /// Delete/dup variant handling.
    #[arg(
        long = "deldupvar",
        help_heading = "Structural Variants & Complex Indels"
    )]
    pub deldupvar: bool,

    // ── 7. Performance & Concurrency ─────────────────────────────────────────
    // Dev note: ports CmdParser.java:162, 617-628 (`readThreadsCount`); `usize::MAX` is the
    // bare-flag sentinel, `--threads` is the GNU long alias, and absent `-th` stays all-cores
    // (or `RAYON_NUM_THREADS`) as a deliberate performance-only default over Java's 1 thread
    // (task P1 report).
    /// Threads (`-th [INT]` / `--threads [INT]`). Bare `-th` uses all available cores; explicit `-th 0` clamps to `1`; omitting `-th` defaults to `RAYON_NUM_THREADS` or all cores.
    #[arg(
        long = "threads",
        num_args(0..=1),
        default_missing_value = "18446744073709551615",
        value_parser = vardictrust_core::config::parse_java_number_usize,
        help_heading = "Performance & Concurrency"
    )]
    pub threads: Option<usize>,

    /// Advanced: region tile size (bp) for bounded-memory WGS. [auto]
    #[arg(
        long = "tile-size",
        num_args = 1,
        value_parser = vardictrust_core::config::parse_java_number_u32,
        help_heading = "Performance & Concurrency"
    )]
    pub tile_size: Vec<u32>,

    /// Positional BED file (when -R absent).
    pub bed: Option<String>,
}

fn all_cores() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

impl Cli {
    /// Convert CLI args to a core Config (first-wins on repeated scalar options per Commons-CLI),
    /// isolated from host `RAYON_NUM_THREADS` for deterministic unit tests.
    pub fn to_config(
        &self,
    ) -> Result<vardictrust_core::config::Config, vardictrust_core::types::VardictError> {
        self.to_config_with_rayon_env(None)
    }

    /// Convert CLI args to a core Config (first-wins on repeated scalar options per Commons-CLI),
    /// consulting `rayon_env` (`RAYON_NUM_THREADS`) when `-th` / `--threads` was not passed.
    pub fn to_config_with_rayon_env(
        &self,
        rayon_env: Option<&str>,
    ) -> Result<vardictrust_core::config::Config, vardictrust_core::types::VardictError> {
        let _ = (self.parity, &self.completions, self.manpage);
        if let Some(ref raw_target) = self.compat {
            let target = raw_target.trim();
            if !target.eq_ignore_ascii_case("vardictjava-1.8.3") && target != "1.8.3" {
                return Err(vardictrust_core::types::VardictError::Unsupported(format!(
                    "--compat {raw_target} not supported in v1.0; use --compat vardictjava-1.8.3 or --parity"
                )));
            }
        }

        if let Some(raw) = self.validation_stringency.first() {
            if raw.eq_ignore_ascii_case("SILENT") || raw.eq_ignore_ascii_case("LENIENT") {
                // Supported no-op: matches noodles lenient decoding.
            } else if raw.eq_ignore_ascii_case("STRICT") {
                return Err(vardictrust_core::types::VardictError::Unsupported(
                    "-VS STRICT not implemented in v1.0; use VarDictJava".to_string(),
                ));
            } else {
                return Err(vardictrust_core::types::VardictError::Parse(format!(
                    "Invalid ValidationStringency '{raw}'; expected SILENT, LENIENT, or STRICT"
                )));
            }
        }

        let flag_filter_raw = self
            .flag_filter
            .first()
            .map(|s| s.as_str())
            .unwrap_or("0x504");
        let decoded_flag =
            vardictrust_core::config::java_integer_decode(flag_filter_raw).map_err(|e| {
                vardictrust_core::types::VardictError::Parse(format!(
                    "Invalid -F value '{}': {}",
                    flag_filter_raw, e
                ))
            })?;
        let flag_filter = decoded_flag as u32;

        let tile_size = self.tile_size.first().copied();
        if let Some(0) = tile_size {
            return Err(vardictrust_core::types::VardictError::Parse(
                "Invalid --tile-size: 0. Tile size must be >= 1".into(),
            ));
        }
        if let Some(n) = tile_size {
            if n > vardictrust_core::MAX_TILE_SIZE {
                return Err(vardictrust_core::types::VardictError::Parse(format!(
                    "Invalid --tile-size: {}. Must be <= {} to avoid u32 wrap",
                    n,
                    vardictrust_core::MAX_TILE_SIZE
                )));
            }
        }

        let freq_str = self.freq.first().map(|s| s.as_str()).unwrap_or("0.01");
        let parsed_freq: f64 = vardictrust_core::config::parse_java_number_f64(freq_str)
            .or_else(|e| {
                if freq_str.eq_ignore_ascii_case("0e0") {
                    Ok(0.0)
                } else {
                    Err(e)
                }
            })
            .map_err(|e| {
                vardictrust_core::types::VardictError::Parse(format!(
                    "Invalid -f value '{}': {}",
                    freq_str, e
                ))
            })?;
        // In Perl (`var2vcf_valid.pl:21`: `my $Freq = $opt_f ? $opt_f : 0.02;`),
        // only `"0"` (or `""`) is falsy; non-`"0"` zero spellings such as `"0.0"`
        // or `"0.00"` are truthy strings in Perl boolean context and keep `0.0`.
        let af_threshold_truthy_zero = parsed_freq == 0.0 && freq_str.trim() != "0";
        let af_threshold = if parsed_freq == 0.0 { 0.0 } else { parsed_freq };

        let sample_regexp = match self.sample_regexp.first() {
            Some(raw) => vardictrust_core::config::validate_sample_regexp(raw)
                .map_err(vardictrust_core::types::VardictError::Parse)?
                .map(|_| raw.clone()),
            None => None,
        };
        let bam_str = self.bam.first().cloned().unwrap_or_default();
        let region_opt = self.region.first().cloned();
        let explicit_sample = self.sample.first().map(String::as_str);
        let sample_explicit = explicit_sample.is_some();
        let sample = vardictrust_core::config::resolve_sample_names(
            &bam_str,
            explicit_sample,
            sample_regexp.as_deref(),
            region_opt.is_some(),
            self.emit_intermediate || self.splice,
        );

        let col_chr = self.col_chr.first().copied();
        let col_start = self.col_start.first().copied();
        let col_end = self.col_end.first().copied();
        let col_seg_start = self.col_seg_start.first().copied();
        let col_seg_end = self.col_seg_end.first().copied();
        let col_gene = self.col_gene.first().copied();

        if region_opt.is_none() {
            for (flag, val) in [
                ("-c", col_chr),
                ("-S", col_start),
                ("-E", col_end),
                ("-s", col_seg_start),
                ("-e", col_seg_end),
            ] {
                if val == Some(0) {
                    return Err(vardictrust_core::types::VardictError::Parse(format!(
                        "Invalid {flag} value: 0. BED column indices are 1-based and must be >= 1"
                    )));
                }
            }
        }

        let amplicon = self.amplicon.first().map(|raw| {
            if let Some((d_str, o_str)) = raw.split_once(':') {
                if d_str.parse::<i32>().is_ok() && o_str.parse::<f64>().is_ok() {
                    return raw.clone();
                }
            }
            "10:0.95".to_string()
        });

        let (vcf_pval, vcf_pval_raw) = match self.vcf_pval.first() {
            Some((val, raw)) => (*val, Some(raw.clone())),
            None => (0.05, None),
        };

        let mut cfg = vardictrust_core::config::Config {
            exception_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            reference: self.reference.first().cloned().unwrap_or_default(),
            bam: bam_str,
            freq: parsed_freq,
            // The writer's $Freq is the literal -f value even under -p (perl
            // receives -f from the pipeline, not java's pileup override).
            af_threshold,
            af_threshold_truthy_zero,
            sample,
            sample_explicit,
            region: region_opt,
            bed: self.bed.clone(),
            delimiter: self
                .delimiter
                .first()
                .cloned()
                .unwrap_or_else(|| "\t".to_string()),
            col_chr: col_chr.map(|c| c.saturating_sub(1)).unwrap_or(2),
            col_start: col_start.map(|s| s.saturating_sub(1)).unwrap_or(6),
            col_end: col_end.map(|e| e.saturating_sub(1)).unwrap_or(7),
            col_seg_start: col_seg_start
                .map(|s| s.saturating_sub(1))
                .or_else(|| col_start.map(|s| s.saturating_sub(1)))
                .unwrap_or(9),
            col_seg_end: col_seg_end
                .map(|e| e.saturating_sub(1))
                .or_else(|| col_end.map(|e| e.saturating_sub(1)))
                .unwrap_or(10),
            col_gene: col_gene
                .map(|g| {
                    if g > 0 {
                        Some(g.saturating_sub(1))
                    } else {
                        None
                    }
                })
                .unwrap_or(Some(12)),
            column_for_chromosome: col_chr.map(|c| c.saturating_sub(1)),
            // Java: bare -th → availableProcessors (CmdParser.java:621-623),
            // -th N (incl. explicit 0) → max(N,1) (:162). Absent flag consults
            // valid RAYON_NUM_THREADS (>= 1) when provided, else all-cores.
            threads: match self.threads {
                None => rayon_env
                    .and_then(|s| s.parse::<usize>().ok())
                    .filter(|&n| n >= 1)
                    .unwrap_or_else(all_cores),
                Some(n) if n == usize::MAX => all_cores(),
                Some(n) => n.max(1),
            },
            realign: match self.realign.first().copied().unwrap_or(1) {
                1 => 1,
                _ => 0,
            },
            min_reads: self.min_reads.first().copied().unwrap_or(2),
            base_qual: self.base_qual.first().copied().unwrap_or(22.5),
            map_qual: self.map_qual.first().copied(),
            read_pos: self.read_pos.first().copied().unwrap_or(5) as f64,
            max_mismatches: self.max_mismatches.first().copied().unwrap_or(8),
            minmatch: self.minmatch.first().copied().unwrap_or(0),
            lofreq: self.min_normal_af.first().copied().unwrap_or(0.05),
            no_sv: self.no_sv,
            sv_min_len: self.sv_min_len.first().copied().unwrap_or(1000),
            emit_intermediate: self.emit_intermediate,
            fisher: self.fisher,
            tile_size,
            amplicon,
            extension: self.extension.first().copied().unwrap_or(0),
            extension_nuc: self.extension_nuc.first().copied().unwrap_or(2),
            three_prime: self.three_prime,
            min_reads_bias: self.min_reads_bias.first().copied().unwrap_or(2),
            qratio: self.qratio.first().copied().unwrap_or(1.5),
            out_mapq: self.out_mapq.first().copied().unwrap_or(0.0),
            indel_size: self.indel_size.first().copied().unwrap_or(50),
            zero_based: self.zero_based.first().is_some_and(|&v| v == 1),
            zero_based_defined: !self.zero_based.is_empty(),
            flag_filter,
            samfilter_raw: flag_filter_raw.to_string(),
            sv_read_len: self.sv_read_len.first().copied().unwrap_or(300),
            sv_flanking: self.sv_flanking.first().copied().unwrap_or(100),
            sv_min_af: self.sv_min_af.first().copied().unwrap_or(4),
            max_region_len: self.max_region_len.first().copied().unwrap_or(1200),
            mfreq: self.mfreq.first().copied().unwrap_or(0.25),
            nmfreq: self.nmfreq.first().copied().unwrap_or(0.1),
            chimeric: self.chimeric,
            deldupvar: self.deldupvar,
            splice: self.splice,
            pileup: self.pileup,
            debug: self.debug,
            // ports CmdParser.java:171-176: bare -DP → "OUT"; unknown
            // values silently fall back to "OUT" (only ERR is honored).
            debug_print: self.debug_print.first().map(|v| {
                if v == "ERR" {
                    v.clone()
                } else {
                    "OUT".to_string()
                }
            }),
            crispr_cutting_site: self.crispr_cutting_site.first().copied(),
            crispr_flanking_size: self.crispr_flanking_size.first().copied(),
            column_header: self.column_header,
            dedup: self.dedup,
            adaptor: self
                .adaptor
                .first()
                .map(|s| s.split(',').map(|x| x.to_string()).collect())
                .unwrap_or_default(),
            verbose: self.verbose,
            unique_mode_alignment_enabled: self.unique_mode_alignment_enabled,
            unique_mode_second_in_pair_enabled: self.unique_mode_second_in_pair_enabled,
            omit_end_tag: self.no_end || !self.col_end.is_empty(),
            include_n: self.include_n,
            output_path: self.output.first().cloned(),
            index_output: self.tbi,
            stats_json_path: self.stats_json.first().cloned(),
            trim_bases_after: self.trim_bases_after.first().copied().unwrap_or(0),
            sample_regexp,
            index_csi: self.index_csi,
            vcf_pass_only: self.pass_only,
            vcf_paired_strict: self.paired_strict,
            vcf_pval,
            vcf_pval_raw,
            vcf_min_splitreads: self.min_splitreads.first().copied().unwrap_or(1),
            vcf_pstd_filter: !self.no_pstd_filter,
            vcf_af: self.vcf_af.first().cloned(),
            vcf_min_depth: self.vcf_min_depth.first().cloned(),
            vcf_min_vd: self.vcf_min_vd.first().cloned(),
            vcf_qmean: self.vcf_qmean.first().cloned(),
            vcf_mapq: self.vcf_mapq.first().cloned(),
            vcf_pmean: self.vcf_pmean.first().cloned(),
            vcf_sn: self.vcf_sn.first().cloned(),
            vcf_max_nm: self.vcf_max_nm.first().cloned(),
            vcf_gt_freq: self.vcf_gt_freq.first().copied(),
            vcf_msi: self.vcf_msi.first().cloned(),
            vcf_cluster_dist: self.vcf_cluster_dist.first().cloned(),
        };
        cfg.apply_pileup_overrides();
        cfg.validate_vcf_thresholds()?;
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
        assert!(check_unsupported(&args(&["-C"])).is_err());
        assert!(check_unsupported(&args(&["--downsample=0.5"])).is_err());
        assert!(check_unsupported(&args(&["-downsample=0.5"])).is_err());
        // -Z is Java's downsampling fraction (CmdParser.java:125-127, 433-440).
        assert!(check_unsupported(&args(&["-Z", "0.7"])).is_err());
        assert!(matches!(
            check_unsupported(&args(&["-Z0.7"])),
            Err(vardictrust_core::types::VardictError::Unsupported(_))
        ));
        for clustered_eq in ["-uZ=0.5", "-pC"] {
            assert!(
                matches!(
                    check_unsupported(&args(&[clustered_eq])),
                    Err(vardictrust_core::types::VardictError::Unsupported(_))
                ),
                "{clustered_eq} must fail with VardictError::Unsupported"
            );
        }
        // -n, -T, -trim, --trim, -VS, --VS are now supported in CLI parsing
        for supported in ["-pT=5", "-un=RG", "-T", "--trim", "-VS", "--VS"] {
            assert!(
                check_unsupported(&args(&[supported, "SILENT"])).is_ok(),
                "{supported} must pass check_unsupported"
            );
        }
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
        assert_eq!(cli.threads, Some(2)); // first-wins (CmdParser.java)
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
        use clap::CommandFactory;

        let err = Cli::try_parse_from(args(&["-H"])).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::DisplayHelp);
        let err = Cli::try_parse_from(args(&["--help"])).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::DisplayHelp);

        let mut cmd = Cli::command();
        assert_eq!(cmd.get_arguments().count(), 89);
        let rendered_help = cmd.render_help().to_string();
        for heading in [
            "Core Inputs & Outputs",
            "VCF Post-Filtering (var2vcf)",
            "Coordinate & BED Column Mapping",
            "Variant Calling & Quality Thresholds",
            "Modality Flags: Somatic & Amplicon",
            "Structural Variants & Complex Indels",
            "Performance & Concurrency",
        ] {
            assert!(
                rendered_help.contains(heading),
                "missing help heading '{heading}' in:\n{rendered_help}"
            );
        }

        let comp_cli = Cli::try_parse_from(args(&["--completions", "bash"])).unwrap();
        assert_eq!(comp_cli.completions, vec![CompletionShell::Bash]);
        assert!(!comp_cli.manpage);

        let man_cli = Cli::try_parse_from(args(&["--manpage"])).unwrap();
        assert!(man_cli.manpage);
        assert!(man_cli.completions.is_empty());

        let first_wins =
            Cli::try_parse_from(args(&["--completions", "bash", "--completions", "zsh"])).unwrap();
        assert_eq!(first_wins.completions.first(), Some(&CompletionShell::Bash));

        let missing_req = Cli::try_parse_from(args(&[])).unwrap_err();
        assert_eq!(missing_req.kind(), ErrorKind::MissingRequiredArgument);

        let conflict =
            Cli::try_parse_from(args(&["--completions", "bash", "--manpage"])).unwrap_err();
        assert_eq!(conflict.kind(), ErrorKind::ArgumentConflict);

        assert!(check_unsupported(&args(&["-completions"])).is_err());
        assert!(check_unsupported(&args(&["-manpage"])).is_err());
        assert!(check_unsupported(&args(&["-N", "-completions"])).is_ok());

        for shell in [
            CompletionShell::Bash,
            CompletionShell::Zsh,
            CompletionShell::Fish,
        ] {
            let mut buf = Vec::new();
            completions_man::write_completions(shell, &mut buf).unwrap();
            let text = String::from_utf8(buf).unwrap();
            for needle in [
                "--vcf-cluster-dist",
                "--completions",
                "-H",
                "-h",
                "-UN",
                "--defaultPrinter",
            ] {
                assert!(
                    text.contains(needle),
                    "{shell:?} completions missing '{needle}'"
                );
            }
            if shell == CompletionShell::Zsh {
                for valid_opt in ["-z-[", "-k-[", "-a-[", "--amplicon=-["] {
                    assert!(
                        text.contains(valid_opt),
                        "zsh completions missing valid optional spec '{valid_opt}'"
                    );
                }
                for invalid_opt in ["-z=-", "-k=-", "-a=-"] {
                    assert!(
                        !text.contains(invalid_opt),
                        "zsh completions must not emit invalid short-flag '{invalid_opt}'"
                    );
                }
            }
        }

        let mut man_buf = Vec::new();
        completions_man::write_manpage(&mut man_buf).unwrap();
        let man_text = String::from_utf8(man_buf).unwrap();
        for needle in [
            r"\-\-vcf\-cluster\-dist",
            r"\-\-completions",
            r"\-H",
            r"\-h",
            r"\\t",
            r"VCF 4.2 (VCFv4.2)",
            r"\-c 1 \-S 2 \-E 3 \-g 4 \-z targets.bed",
            r"\-c 1 \-S 2 \-E 3 \-z targets.bed",
            r"\-z [ZERO_BASED]",
            r"\-k [REALIGN]",
            r"\-a [AMPLICON], \-\-amplicon[=AMPLICON]",
        ] {
            assert!(
                man_text.contains(needle),
                "manpage missing escaped '{needle}'"
            );
        }
        for forbidden in ["VCF 4.3", "usize::MAX", "task P1", r"\-z[=", r"\-k[="] {
            assert!(
                !man_text.contains(forbidden),
                "manpage must not contain '{forbidden}'"
            );
        }
        assert_eq!(
            completions_man::escape_roff(".var-row\n'tumor|normal'\\path\u{2013}x\u{2014}y"),
            "\\&.var\\-row\n\\&'tumor|normal'\\\\path\\(enx\\(emy"
        );
        if std::process::Command::new("groff")
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success())
        {
            use std::io::Write as _;
            let mut groff = std::process::Command::new("groff")
                .args(["-Tutf8", "-man", "-ww", "-z"])
                .stdin(std::process::Stdio::piped())
                .stdout(std::process::Stdio::piped())
                .stderr(std::process::Stdio::piped())
                .spawn()
                .expect("spawn groff");
            groff
                .stdin
                .take()
                .expect("groff stdin")
                .write_all(man_text.as_bytes())
                .expect("write man_text to groff");
            let out = groff.wait_with_output().expect("wait groff");
            let stderr = String::from_utf8_lossy(&out.stderr);
            assert!(
                out.status.success() && stderr.is_empty(),
                "groff -Tutf8 -man -ww -z emitted warnings or failed (status={:?}):\n{stderr}",
                out.status
            );
        }
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

        let cli_neg_zero =
            Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-f=-0.0"])).unwrap();
        let cfg_neg_zero = cli_neg_zero.to_config().unwrap();
        assert_eq!(cfg_neg_zero.af_threshold.to_bits(), 0.0_f64.to_bits());
        assert!(cfg_neg_zero.af_threshold_truthy_zero);
    }

    #[test]
    fn unsupported_flags_all_seven_return_clear_error() {
        for flag in ["-C", "-Z", "-downsample", "--downsample"] {
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

        // -VS / --VS: SILENT & LENIENT succeed (case-insensitive, first-wins); STRICT -> Unsupported;
        // anything else -> Parse error.
        for ok_vs in [
            &["-VS", "SILENT"][..],
            &["-VS", "LENIENT"][..],
            &["--VS", "silent"][..],
            &["-VS=lenient"][..],
            &["--VS=SILENT"][..],
            &["-VS", "LENIENT", "-VS", "STRICT"][..],
        ] {
            let mut full = vec!["-G", "r.fa", "-b", "a.bam"];
            full.extend_from_slice(ok_vs);
            let norm = normalize_args(&args(&full));
            check_unsupported(&norm).unwrap();
            assert!(
                Cli::try_parse_from(&norm).unwrap().to_config().is_ok(),
                "expected Ok for {ok_vs:?}"
            );
        }
        for strict_vs in [
            &["-VS", "STRICT"][..],
            &["--VS", "STRICT"][..],
            &["-VS=strict"][..],
            &["--VS=STRICT"][..],
            &["-VS", "STRICT", "-VS", "LENIENT"][..],
        ] {
            let mut full = vec!["-G", "r.fa", "-b", "a.bam"];
            full.extend_from_slice(strict_vs);
            let norm = normalize_args(&args(&full));
            check_unsupported(&norm).unwrap();
            let err = Cli::try_parse_from(&norm).unwrap().to_config().unwrap_err();
            assert!(
                matches!(err, vardictrust_core::types::VardictError::Unsupported(_)),
                "expected Unsupported for {strict_vs:?}, got {err:?}"
            );
            assert!(
                err.to_string()
                    .contains("-VS STRICT not implemented in v1.0; use VarDictJava"),
                "unexpected STRICT error for {strict_vs:?}: {err}"
            );
        }
        let bad_vs = normalize_args(&args(&["-G", "r.fa", "-b", "a.bam", "-VS", "BOGUS"]));
        let bad_vs_err = Cli::try_parse_from(&bad_vs)
            .unwrap()
            .to_config()
            .unwrap_err();
        assert!(matches!(
            bad_vs_err,
            vardictrust_core::types::VardictError::Parse(_)
        ));
        assert!(bad_vs_err
            .to_string()
            .contains("Invalid ValidationStringency 'BOGUS'; expected SILENT, LENIENT, or STRICT"));
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
            assert_eq!(
                cli.to_config_with_rayon_env(Some("2")).unwrap().threads,
                expected_threads,
                "explicit -th must override RAYON_NUM_THREADS for {argv:?}"
            );
        }

        // RAYON_NUM_THREADS precedence: explicit -th N / bare -th / -th 0 > RAYON_NUM_THREADS >= 1 > all_cores()
        let no_th = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa"])).unwrap();
        assert_eq!(
            no_th.to_config_with_rayon_env(Some("3")).unwrap().threads,
            3
        );
        assert_eq!(
            no_th.to_config_with_rayon_env(Some("1")).unwrap().threads,
            1
        );
        for invalid_env in [None, Some(""), Some("0"), Some("bogus"), Some("-2")] {
            assert_eq!(
                no_th.to_config_with_rayon_env(invalid_env).unwrap().threads,
                all_cores(),
                "invalid RAYON_NUM_THREADS={invalid_env:?} must fall back to all_cores()"
            );
        }
        let bare_th =
            Cli::try_parse_from(normalize_args(&args(&["-b", "a.bam", "-G", "r.fa", "-th"])))
                .unwrap();
        assert_eq!(
            bare_th.to_config_with_rayon_env(Some("2")).unwrap().threads,
            all_cores()
        );
        let zero_th = Cli::try_parse_from(normalize_args(&args(&[
            "-b", "a.bam", "-G", "r.fa", "-th", "0",
        ])))
        .unwrap();
        assert_eq!(
            zero_th.to_config_with_rayon_env(Some("6")).unwrap().threads,
            1
        );
    }

    #[test]
    fn java_number_parse_table() {
        use vardictrust_core::config::{parse_java_number_f64, parse_java_number_i32};

        assert_eq!(parse_java_number_i32("2.9").unwrap(), 2);
        assert_eq!(parse_java_number_i32("2.9f").unwrap(), 2);
        assert_eq!(parse_java_number_i32(" 2.9d ").unwrap(), 2);
        assert_eq!(parse_java_number_i32("-100").unwrap(), -100);
        assert_eq!(parse_java_number_f64("22.5f").unwrap(), 22.5);
        assert_eq!(parse_java_number_f64("22.5D").unwrap(), 22.5);
        assert_eq!(parse_java_number_f64("1.5e2").unwrap(), 150.0);

        for bad in ["1e1", "NaN", "inf", "Infinity", " 10", "10 ", ""] {
            assert!(
                parse_java_number_i32(bad).is_err(),
                "expected i32 parse error on '{bad}'"
            );
            assert!(
                parse_java_number_f64(bad).is_err(),
                "expected f64 parse error on '{bad}'"
            );
        }

        // CLI integration: -P truncates float to int, -O accepts float, -x accepts negative,
        // repeated scalar options use first-wins, -a "" falls back to 10:0.95, -z 0 disables zero_based.
        let norm = normalize_args(&args(&[
            "-b",
            "a.bam",
            "-G",
            "r.fa",
            "-P",
            "2.9",
            "-O",
            "12.5",
            "-x",
            "-5",
            "-m",
            "5",
            "-m",
            "10",
            "-a",
            "",
            "-z",
            "0",
            "regions.bed",
        ]));
        let cli = Cli::try_parse_from(&norm).unwrap();
        let cfg = cli.to_config().unwrap();
        assert_eq!(cfg.read_pos, 2.0);
        assert_eq!(cfg.out_mapq, 12.5);
        assert_eq!(cfg.extension, -5);
        assert_eq!(cfg.max_mismatches, 5);
        assert_eq!(cfg.amplicon.as_deref(), Some("10:0.95"));
        assert!(!cfg.zero_based);
        assert!(cfg.zero_based_defined);
        assert_eq!(cfg.bed.as_deref(), Some("regions.bed"));
        assert_eq!(cfg.trim_bases_after, 0);

        // -z 2 disables zero_based (CmdParser.java:75 `1 == getIntValue(cmd, "z", 1)`)
        let norm_z2 = normalize_args(&args(&["-b", "a.bam", "-G", "r.fa", "-z", "2", "r.bed"]));
        let cfg_z2 = Cli::try_parse_from(&norm_z2).unwrap().to_config().unwrap();
        assert!(!cfg_z2.zero_based);
        assert!(cfg_z2.zero_based_defined);

        // -z -3 leaves -z as bare (-z=1) and preserves -3 (three_prime = true)
        let norm_z_3 = normalize_args(&args(&["-b", "a.bam", "-G", "r.fa", "-z", "-3", "r.bed"]));
        let cfg_z_3 = Cli::try_parse_from(&norm_z_3).unwrap().to_config().unwrap();
        assert!(cfg_z_3.zero_based);
        assert!(cfg_z_3.zero_based_defined);
        assert!(cfg_z_3.three_prime);

        // -th 2 -th 8.0 consumes 8.0 via Java number lookahead and keeps first thread count 2
        let norm_th = normalize_args(&args(&[
            "-b", "a.bam", "-G", "r.fa", "-th", "2", "-th", "8.0", "r.bed",
        ]));
        let cfg_th = Cli::try_parse_from(&norm_th).unwrap().to_config().unwrap();
        assert_eq!(cfg_th.threads, 2);

        // -x -3 and -T -3 / -trim -3 / --trim -3 must fail (Commons-CLI treats -3 as flag)
        assert!(Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-x", "-3"])).is_err());
        for bad_trim in [
            &["-T", "-3"][..],
            &["-trim", "-3"][..],
            &["--trim", "-3"][..],
        ] {
            let mut full = vec!["-b", "a.bam", "-G", "r.fa"];
            full.extend_from_slice(bad_trim);
            let norm_bad = normalize_args(&args(&full));
            assert!(
                Cli::try_parse_from(&norm_bad).is_err(),
                "{bad_trim:?} must fail because -3 is moveIndelsTo3"
            );
        }

        // -T / -trim / --trim: first-wins, negative (-5), float truncation (2.9 -> 2), attached -T50
        for (trim_argv, expected_trim) in [
            (&["-T", "50", "-T", "20"][..], 50),
            (&["-trim", "35", "--trim", "10"][..], 35),
            (&["-trim=40"][..], 40),
            (&["--trim", "45"][..], 45),
            (&["--trim=48"][..], 48),
            (&["-T", "-5"][..], -5),
            (&["-T", "2.9"][..], 2),
            (&["-T50"][..], 50),
            (&["-pT=5"][..], 5),
        ] {
            let mut full = vec!["-b", "a.bam", "-G", "r.fa"];
            full.extend_from_slice(trim_argv);
            let norm_t = normalize_args(&args(&full));
            check_unsupported(&norm_t).unwrap();
            let cfg_t = Cli::try_parse_from(&norm_t).unwrap().to_config().unwrap();
            assert_eq!(cfg_t.trim_bases_after, expected_trim, "argv={trim_argv:?}");
        }
    }

    #[test]
    fn output_and_stats_json_flags_parse_with_first_wins() {
        let default_cli = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa"])).unwrap();
        let default_cfg = default_cli.to_config().unwrap();
        assert_eq!(default_cfg.output_path, None);
        assert_eq!(default_cfg.stats_json_path, None);
        assert!(!default_cfg.index_csi);
        assert!(!default_cfg.omit_end_tag);

        let cli = Cli::try_parse_from(args(&[
            "-b",
            "a.bam",
            "-G",
            "r.fa",
            "--output",
            "first.vcf.gz",
            "--output",
            "second.vcf",
            "--stats-json=first_stats.json",
            "--stats-json",
            "second_stats.json",
        ]))
        .unwrap();
        let cfg = cli.to_config().unwrap();
        assert_eq!(cfg.output_path.as_deref(), Some("first.vcf.gz"));
        assert_eq!(cfg.stats_json_path.as_deref(), Some("first_stats.json"));
        assert!(!cfg.index_output);
        assert!(!cfg.index_csi);

        let tbi_cli = Cli::try_parse_from(args(&[
            "-b",
            "a.bam",
            "-G",
            "r.fa",
            "--output",
            "out.vcf.gz",
            "--tbi",
            "--parity",
            "--compat",
            "vardictjava-1.8.3",
        ]))
        .unwrap();
        assert!(tbi_cli.parity);
        assert!(tbi_cli.to_config().unwrap().index_output);

        let idx_cli = Cli::try_parse_from(args(&[
            "-b",
            "a.bam",
            "-G",
            "r.fa",
            "--output",
            "out.vcf.gz",
            "--index",
            "--compat=1.8.3",
        ]))
        .unwrap();
        assert!(idx_cli.to_config().unwrap().index_output);

        // --csi and --no-end / -E omit_end_tag combinations
        let csi_no_end_cli = Cli::try_parse_from(args(&[
            "-b",
            "a.bam",
            "-G",
            "r.fa",
            "--output",
            "out.vcf.gz",
            "--csi",
            "--no-end",
        ]))
        .unwrap();
        let csi_no_end_cfg = csi_no_end_cli.to_config().unwrap();
        assert!(csi_no_end_cfg.index_csi);
        assert!(!csi_no_end_cfg.index_output);
        assert!(csi_no_end_cfg.omit_end_tag);

        let col_end_cfg = Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "-E", "3"]))
            .unwrap()
            .to_config()
            .unwrap();
        assert!(col_end_cfg.omit_end_tag);

        let bad_compat_cli =
            Cli::try_parse_from(args(&["-b", "a.bam", "-G", "r.fa", "--compat", "modern"]))
                .unwrap();
        assert!(matches!(
            bad_compat_cli.to_config().unwrap_err(),
            vardictrust_core::types::VardictError::Unsupported(_)
        ));

        check_var2vcf_postfilter_cli_flags();
    }

    fn check_var2vcf_postfilter_cli_flags() {
        let parse_pipeline = |raw_flags: &[&str]| -> vardictrust_core::config::Config {
            let norm = normalize_args(&args(raw_flags));
            check_unsupported(&norm).unwrap_or_else(|e| {
                panic!("unexpected check_unsupported error on {raw_flags:?} (norm={norm:?}): {e}")
            });
            let cli = Cli::try_parse_from(&norm).unwrap_or_else(|e| {
                panic!("unexpected Cli parse error on {raw_flags:?} (norm={norm:?}): {e}")
            });
            cli.to_config()
                .unwrap_or_else(|e| panic!("unexpected to_config error on {raw_flags:?}: {e}"))
        };

        // 1. Default Config (both Config::default() and minimal CLI args)
        let core_default = vardictrust_core::config::Config::default();
        assert!(!core_default.vcf_pass_only && !core_default.vcf_paired_strict);
        assert!(
            (core_default.vcf_pval - 0.05).abs() < f64::EPSILON
                && core_default.vcf_pval_raw.is_none()
        );
        assert!(core_default.vcf_min_splitreads == 1 && core_default.vcf_pstd_filter);

        let cfg = parse_pipeline(&["-G", "ref.fa", "-b", "in.bam", "-R", "chr1:1-100"]);
        assert!(!cfg.vcf_pass_only && !cfg.vcf_paired_strict);
        assert!((cfg.vcf_pval - 0.05).abs() < f64::EPSILON && cfg.vcf_pval_raw.is_none());
        assert!(cfg.vcf_min_splitreads == 1 && cfg.vcf_pstd_filter);

        // 2. Explicit overrides & first-wins repetition, scientific notation, trailing zeros, -0.0 canonicalization
        let sci_cfg = parse_pipeline(&[
            "-G",
            "ref.fa",
            "-b",
            "in.bam",
            "-R",
            "chr1:1-100",
            "--pass-only",
            "--paired-strict",
            "--no-pstd-filter",
            "--vcf-pval",
            "1e-5",
            "--vcf-pval",
            "0.02",
            "--min-splitreads",
            "3",
            "--min-splitreads",
            "9",
        ]);
        assert!(sci_cfg.vcf_pass_only);
        assert!(sci_cfg.vcf_paired_strict);
        assert!(!sci_cfg.vcf_pstd_filter);
        assert!((sci_cfg.vcf_pval - 1e-5).abs() < 1e-15);
        assert_eq!(sci_cfg.vcf_pval_raw.as_deref(), Some("1e-5"));
        assert_eq!(sci_cfg.vcf_min_splitreads, 3);

        let tz_cfg = parse_pipeline(&[
            "-G",
            "ref.fa",
            "-b",
            "in.bam",
            "-R",
            "chr1:1-100",
            "--vcf-pval=0.050",
            "--min-splitreads=0",
        ]);
        assert!((tz_cfg.vcf_pval - 0.05).abs() < f64::EPSILON);
        assert_eq!(tz_cfg.vcf_pval_raw.as_deref(), Some("0.050"));
        assert_eq!(tz_cfg.vcf_min_splitreads, 0);

        let neg_zero_cfg = parse_pipeline(&[
            "-G",
            "ref.fa",
            "-b",
            "in.bam",
            "-R",
            "chr1:1-100",
            "--vcf-pval",
            "-0.0",
        ]);
        assert_eq!(neg_zero_cfg.vcf_pval, 0.0);
        assert!(!neg_zero_cfg.vcf_pval.is_sign_negative());
        assert_eq!(neg_zero_cfg.vcf_pval_raw.as_deref(), Some("0.0"));

        let pos_zero_cfg = parse_pipeline(&[
            "-G",
            "ref.fa",
            "-b",
            "in.bam",
            "-R",
            "chr1:1-100",
            "--vcf-pval=+0.00",
        ]);
        assert_eq!(pos_zero_cfg.vcf_pval, 0.0);
        assert!(!pos_zero_cfg.vcf_pval.is_sign_negative());
        assert_eq!(pos_zero_cfg.vcf_pval_raw.as_deref(), Some("0.00"));

        // 3. Hostile value & missing-value rejection
        for bad_flags in [
            &["--vcf-pval", "-0.01", "-R", "chr1:1-100"][..],
            &["--vcf-pval", "-3", "-R", "chr1:1-100"][..],
            &["--vcf-pval", "NaN", "-R", "chr1:1-100"][..],
            &["--vcf-pval", "Infinity", "-R", "chr1:1-100"][..],
            &["--vcf-pval", "inf", "-R", "chr1:1-100"][..],
            &["--vcf-pval", "0.01f", "-R", "chr1:1-100"][..],
            &["--vcf-pval", "", "-R", "chr1:1-100"][..],
            &["--vcf-pval", " 0.05", "-R", "chr1:1-100"][..],
            &["--min-splitreads", "-1", "-R", "chr1:1-100"][..],
            &["--min-splitreads", "-3", "-R", "chr1:1-100"][..],
            &["--min-splitreads", "2.5", "-R", "chr1:1-100"][..],
            &["--min-splitreads", "0x2", "-R", "chr1:1-100"][..],
            &["--min-splitreads", "", "-R", "chr1:1-100"][..],
            &["--min-splitreads", " 1", "-R", "chr1:1-100"][..],
            &["--vcf-pval", "--pass-only", "-R", "chr1:1-100"][..],
            &["--min-splitreads", "-R", "chr1:1-100"][..],
        ] {
            let mut full = vec!["-G", "ref.fa", "-b", "in.bam"];
            full.extend_from_slice(bad_flags);
            let norm = normalize_args(&args(&full));
            assert!(
                check_unsupported(&norm).is_err() || Cli::try_parse_from(&norm).is_err(),
                "expected CLI parse failure for {bad_flags:?} (norm={norm:?})"
            );
            assert!(
                Cli::try_parse_from(&norm).is_err(),
                "expected Cli::try_parse_from failure for {bad_flags:?} (norm={norm:?})"
            );
        }

        // 4. Single-dash GNU-flag typo rejection in check_unsupported (bare, =val, =10:0.95, and separate value token)
        for &gnu_typo in GNU_ONLY_LONG_FLAGS_SINGLE_DASH {
            let expected_suggestion = format!("did you mean '-{gnu_typo}'?");
            let norm_bare = normalize_args(&args(&[
                "-G", "r.fa", "-b", "a.bam", gnu_typo, "-R", "c:1-10",
            ]));
            let err_bare = check_unsupported(&norm_bare).expect_err(&format!(
                "expected Parse error for single-dash GNU flag {gnu_typo} (norm={norm_bare:?})"
            ));
            assert!(
                matches!(err_bare, vardictrust_core::types::VardictError::Parse(_)),
                "expected VardictError::Parse for {gnu_typo}, got {err_bare:?}"
            );
            assert!(
                err_bare.to_string().contains(&expected_suggestion),
                "missing suggestion '{expected_suggestion}' in error for {gnu_typo}: {err_bare}"
            );

            for eq_suffix in ["0.01", "10:0.95"] {
                let eq_tok = format!("{gnu_typo}={eq_suffix}");
                let norm_eq = normalize_args(&args(&[
                    "-G",
                    "r.fa",
                    "-b",
                    "a.bam",
                    eq_tok.as_str(),
                    "-R",
                    "c:1-10",
                ]));
                let err_eq = check_unsupported(&norm_eq).expect_err(&format!(
                    "expected Parse error for single-dash GNU flag {eq_tok} (norm={norm_eq:?})"
                ));
                assert!(
                    matches!(err_eq, vardictrust_core::types::VardictError::Parse(_)),
                    "expected VardictError::Parse for {eq_tok}, got {err_eq:?}"
                );
                assert!(
                    err_eq.to_string().contains(&expected_suggestion),
                    "missing suggestion '{expected_suggestion}' in error for {eq_tok}: {err_eq}"
                );
            }
        }

        // Separate value tokens after single-dash GNU flags must also be rejected cleanly
        for (sep_argv, expected_flag) in [
            (
                &[
                    "-G",
                    "r.fa",
                    "-b",
                    "a.bam",
                    "-vcf-pval",
                    "0.01",
                    "-R",
                    "c:1-10",
                ][..],
                "--vcf-pval",
            ),
            (
                &[
                    "-G", "r.fa", "-b", "a.bam", "-output", "x.vcf", "-R", "c:1-10",
                ][..],
                "--output",
            ),
            (
                &[
                    "-G",
                    "r.fa",
                    "-b",
                    "a.bam",
                    "--vcf-pval",
                    "0.01",
                    "-output",
                    "x.vcf",
                    "-R",
                    "c:1-10",
                ][..],
                "--output",
            ),
            (
                &[
                    "-G",
                    "r.fa",
                    "-b",
                    "a.bam",
                    "--output",
                    "x.vcf",
                    "-min-splitreads",
                    "2",
                    "-R",
                    "c:1-10",
                ][..],
                "--min-splitreads",
            ),
        ] {
            let norm_sep = normalize_args(&args(sep_argv));
            let err_sep = check_unsupported(&norm_sep).expect_err(&format!(
                "expected Parse error for {sep_argv:?} (norm={norm_sep:?})"
            ));
            assert!(
                matches!(err_sep, vardictrust_core::types::VardictError::Parse(_)),
                "expected VardictError::Parse for {sep_argv:?}, got {err_sep:?}"
            );
            assert!(
                err_sep
                    .to_string()
                    .contains(&format!("did you mean '{expected_flag}'?")),
                "unexpected error message for {sep_argv:?}: {err_sep}"
            );
        }

        // Valid values passed to `-n` (including hyphenated values like `-pass-only` and normal regexes)
        // must NOT be rejected by `check_unsupported` per `takes_mandatory_next_token`.
        let n_hyphen_cfg = parse_pipeline(&[
            "-G",
            "r.fa",
            "-b",
            "sample-pass-only.bam",
            "-n",
            "-pass-only",
            "-R",
            "c:1-10",
        ]);
        assert_eq!(n_hyphen_cfg.sample_regexp.as_deref(), Some("-pass-only"));
        assert!(!n_hyphen_cfg.vcf_pass_only);

        let n_regex_cfg = parse_pipeline(&[
            "-G",
            "r.fa",
            "-b",
            "cohort_S42.bam",
            "-n",
            "cohort_(S\\d+)",
            "--pass-only",
            "-R",
            "c:1-10",
        ]);
        assert_eq!(n_regex_cfg.sample, "S42");
        assert!(n_regex_cfg.vcf_pass_only);

        // Short boolean clusters (-pz, -pz0, -uk, -ua, -pR) continue to normalize cleanly
        for cluster_ok in [
            &["-G", "r.fa", "-b", "a.bam", "-pz", "regions.bed"][..],
            &["-G", "r.fa", "-b", "a.bam", "-pz0", "regions.bed"][..],
            &["-G", "r.fa", "-b", "a.bam", "-uk", "regions.bed"][..],
            &["-G", "r.fa", "-b", "a.bam", "-ua", "regions.bed"][..],
            &["-G", "r.fa", "-b", "a.bam", "-pR", "chr1:1-100"][..],
        ] {
            let _ = parse_pipeline(cluster_ok);
        }

        // 5. normalize_args token-consumption safety with --vcf-pval and --min-splitreads before -z regions.bed
        let norm_z = normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "a.bam",
            "--vcf-pval",
            "1e-5",
            "--min-splitreads",
            "0",
            "-z",
            "regions.bed",
        ]));
        assert!(
            norm_z.iter().any(|t| t == "-z=1"),
            "expected -z to normalize to -z=1, got {norm_z:?}"
        );
        check_unsupported(&norm_z).unwrap();
        let z_cfg = Cli::try_parse_from(&norm_z).unwrap().to_config().unwrap();
        assert!(z_cfg.zero_based);
        assert!(z_cfg.zero_based_defined);
        assert_eq!(z_cfg.bed.as_deref(), Some("regions.bed"));
        assert!((z_cfg.vcf_pval - 1e-5).abs() < 1e-15);
        assert_eq!(z_cfg.vcf_pval_raw.as_deref(), Some("1e-5"));
        assert_eq!(z_cfg.vcf_min_splitreads, 0);
    }

    #[test]
    fn commons_cli_value_aware_unsupported_and_hyphen_values() {
        for bad_argv in [
            &["-downsample", "0.5"][..],
            &["-downsample=0.5"][..],
            &["-Z0.7"][..],
            &["-C"][..],
            &["-uZ=0.5"][..],
        ] {
            let norm = normalize_args(&args(bad_argv));
            let err = check_unsupported(&norm).expect_err(&format!(
                "expected Unsupported for {bad_argv:?} (norm={norm:?})"
            ));
            assert!(
                matches!(err, vardictrust_core::types::VardictError::Unsupported(_)),
                "expected Unsupported for {bad_argv:?}, got {err:?}"
            );
        }

        // -n (sampleNameRegexp): first-wins, outer `/` stripping, resolve_sample_names wiring, invalid regex -> Parse
        let n_cfg = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "/dir/sample_42.bam",
            "-n",
            "/sample_(\\d+)/",
            "-n",
            "/ignored_(\\w+)/",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(n_cfg.sample_regexp.as_deref(), Some("/sample_(\\d+)/"));
        assert_eq!(n_cfg.sample, "42");

        let empty_slash_n = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "/dir/foo.sorted.bam",
            "-n",
            "//",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(empty_slash_n.sample_regexp, None);
        assert_eq!(empty_slash_n.sample, "foo");

        let escaped_slash_n = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "/mnt/cohort/SAMP_ESC/aln.bam",
            "-n",
            r"/cohort\/([^/]+)\//",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(
            escaped_slash_n.sample_regexp.as_deref(),
            Some(r"/cohort\/([^/]+)\//")
        );
        assert_eq!(escaped_slash_n.sample, "SAMP_ESC");

        let double_slash_n = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "not_cohort/FAKE/aln_1.bam|/cohort/REAL/aln_2.bam",
            "-n",
            r"//cohort/([^/]+)//",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(double_slash_n.sample, "aln_1|REAL");

        let paired_vcf_n = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "dir/s_T9.bam|dir/s_N8.bam",
            "-n",
            "s_([A-Z0-9]+)\\.bam",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(paired_vcf_n.sample, "T9|N8");

        let paired_splice_r_n = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "t_1.bam|n_2.bam",
            "-R",
            "chr1:1-100",
            "-i",
            "-N",
            "EXP_T|EXP_N",
            "-n",
            "t_(\\d+)",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(paired_splice_r_n.sample, "1");

        let bad_n = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "a.bam",
            "-n",
            "([unclosed",
        ])))
        .unwrap()
        .to_config()
        .unwrap_err();
        assert!(matches!(
            bad_n,
            vardictrust_core::types::VardictError::Parse(_)
        ));

        for sample_val in ["-fisher", "-Tumor", "-Control", "-normal", "-trim"] {
            let raw = args(&["-G", "r.fa", "-b", "a.bam", "-N", sample_val]);
            let norm = normalize_args(&raw);
            check_unsupported(&norm)
                .unwrap_or_else(|e| panic!("-N {sample_val} rejected by check_unsupported: {e}"));
            let cli = Cli::try_parse_from(&norm)
                .unwrap_or_else(|e| panic!("-N {sample_val} rejected by Cli: {e}"));
            let cfg = cli.to_config().unwrap();
            assert_eq!(cfg.sample, sample_val);
            assert!(!cfg.fisher, "-N -fisher must not enable --fisher");
        }

        let out_raw = args(&["-G", "r.fa", "-b", "a.bam", "--output", "-Z0.7"]);
        let out_norm = normalize_args(&out_raw);
        check_unsupported(&out_norm).expect("--output -Z0.7 must pass check_unsupported");
        let out_cfg = Cli::try_parse_from(&out_norm)
            .expect("--output -Z0.7 must parse")
            .to_config()
            .unwrap();
        assert_eq!(out_cfg.output_path.as_deref(), Some("-Z0.7"));

        let dd_raw = args(&["-G", "r.fa", "-b", "a.bam", "--", "-downsample.bed"]);
        let dd_norm = normalize_args(&dd_raw);
        check_unsupported(&dd_norm).expect("-- -downsample.bed must pass check_unsupported");
        let dd_cfg = Cli::try_parse_from(&dd_norm)
            .expect("-- -downsample.bed must parse")
            .to_config()
            .unwrap();
        assert_eq!(dd_cfg.bed.as_deref(), Some("-downsample.bed"));
    }

    #[test]
    fn commons_cli_single_dash_long_aliases_help_and_realign() {
        use clap::error::ErrorKind;
        use clap::CommandFactory as _;

        let ad_norm = normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "a.bam",
            "-adaptor",
            "ACGTTGCTC",
        ]));
        check_unsupported(&ad_norm).unwrap();
        let ad_cfg = Cli::try_parse_from(&ad_norm).unwrap().to_config().unwrap();
        assert!(
            ad_cfg.amplicon.is_none(),
            "-adaptor must not collide with -a amplicon"
        );
        assert_eq!(ad_cfg.adaptor, vec!["ACGTTGCTC".to_string()]);

        let amp_norm = normalize_args(&args(&["-G", "r.fa", "-b", "a.bam", "-amplicon", "8:0.9"]));
        let amp_cfg = Cli::try_parse_from(&amp_norm).unwrap().to_config().unwrap();
        assert_eq!(amp_cfg.amplicon.as_deref(), Some("8:0.9"));

        let bool_norm = normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "a.bam",
            "-fisher",
            "-chimeric",
            "-deldupvar",
            "-nosv",
            "-splice",
            "-header",
            "-dedup",
            "-debug",
            "-verbose",
        ]));
        check_unsupported(&bool_norm).unwrap();
        let bool_cfg = Cli::try_parse_from(&bool_norm)
            .unwrap()
            .to_config()
            .unwrap();
        assert!(bool_cfg.fisher);
        assert!(bool_cfg.chimeric);
        assert!(bool_cfg.deldupvar);
        assert!(bool_cfg.no_sv);
        assert!(bool_cfg.splice);
        assert!(bool_cfg.column_header);
        assert!(bool_cfg.dedup);
        assert!(bool_cfg.debug);
        assert!(bool_cfg.verbose);

        for dp_argv in [
            &["-DP", "ERR"][..],
            &["-default-printer", "ERR"][..],
            &["--default-printer", "ERR"][..],
        ] {
            let mut full = vec!["-G", "r.fa", "-b", "a.bam"];
            full.extend_from_slice(dp_argv);
            let norm = normalize_args(&args(&full));
            check_unsupported(&norm).unwrap();
            let cfg = Cli::try_parse_from(&norm).unwrap().to_config().unwrap();
            assert_eq!(cfg.debug_print.as_deref(), Some("ERR"), "argv={dp_argv:?}");
        }

        // Bare optional-value flags (-DP, --DP, -a, -ua) before a positional BED path must NOT
        // swallow the positional BED filename as their optional value.
        for (opt_argv, exp_dp, exp_amp) in [
            (&["-DP", "regions.bed"][..], Some("OUT"), None),
            (&["--DP", "regions.bed"][..], Some("OUT"), None),
            (&["-DP", "ERR", "regions.bed"][..], Some("ERR"), None),
            (&["-a", "regions.bed"][..], None, Some("10:0.95")),
            (&["-ua", "regions.bed"][..], None, Some("10:0.95")),
            (&["-ua", "30:0.5", "regions.bed"][..], None, Some("30:0.5")),
            (&["-ua30:0.5", "regions.bed"][..], None, Some("30:0.5")),
            (&["-a", "30:0.5", "regions.bed"][..], None, Some("30:0.5")),
        ] {
            let mut full = vec!["-G", "r.fa", "-b", "a.bam"];
            full.extend_from_slice(opt_argv);
            let norm = normalize_args(&args(&full));
            check_unsupported(&norm).unwrap();
            let cfg = Cli::try_parse_from(&norm)
                .unwrap_or_else(|e| panic!("failed on {opt_argv:?} (norm={norm:?}): {e}"))
                .to_config()
                .unwrap();
            assert_eq!(cfg.bed.as_deref(), Some("regions.bed"), "argv={opt_argv:?}");
            assert_eq!(cfg.debug_print.as_deref(), exp_dp, "argv={opt_argv:?}");
            assert_eq!(cfg.amplicon.as_deref(), exp_amp, "argv={opt_argv:?}");
        }

        // Safe boolean-prefix short-cluster splitter (-pz, -pz0, -uz 0, -uk, -pR)
        for (cluster_argv, exp_pileup, exp_u, exp_zb, exp_realign) in [
            (&["-pz", "regions.bed"][..], true, false, true, 1u8),
            (&["-pz0", "regions.bed"][..], true, false, false, 1u8),
            (&["-pz=0", "regions.bed"][..], true, false, false, 1u8),
            (&["-uz", "0", "regions.bed"][..], false, true, false, 1u8),
            (&["-uk", "regions.bed"][..], false, true, false, 1u8),
            (&["-uk", "0", "regions.bed"][..], false, true, false, 0u8),
            (&["-uk0", "regions.bed"][..], false, true, false, 0u8),
        ] {
            let mut full = vec!["-G", "r.fa", "-b", "a.bam"];
            full.extend_from_slice(cluster_argv);
            let norm = normalize_args(&args(&full));
            check_unsupported(&norm).unwrap();
            let cfg = Cli::try_parse_from(&norm)
                .unwrap_or_else(|e| panic!("failed on {cluster_argv:?} (norm={norm:?}): {e}"))
                .to_config()
                .unwrap();
            assert_eq!(cfg.pileup, exp_pileup, "argv={cluster_argv:?}");
            assert_eq!(
                cfg.unique_mode_alignment_enabled, exp_u,
                "argv={cluster_argv:?}"
            );
            assert_eq!(cfg.zero_based, exp_zb, "argv={cluster_argv:?}");
            assert_eq!(cfg.realign, exp_realign, "argv={cluster_argv:?}");
            assert_eq!(
                cfg.bed.as_deref(),
                Some("regions.bed"),
                "argv={cluster_argv:?}"
            );
        }

        // Clustered `-a` with unsigned and signed amplicon values (`-ua10:0.95`, `-ua+30:0.5`, `-pa-5:0.95`)
        for (amp_cluster, exp_pileup, exp_u, exp_amp) in [
            (&["-ua10:0.95", "regions.bed"][..], false, true, "10:0.95"),
            (&["-ua+30:0.5", "regions.bed"][..], false, true, "+30:0.5"),
            (&["-pa-5:0.95", "regions.bed"][..], true, false, "-5:0.95"),
        ] {
            let mut full = vec!["-G", "r.fa", "-b", "a.bam"];
            full.extend_from_slice(amp_cluster);
            let norm = normalize_args(&args(&full));
            check_unsupported(&norm).unwrap();
            let cfg = Cli::try_parse_from(&norm)
                .unwrap_or_else(|e| panic!("failed on {amp_cluster:?} (norm={norm:?}): {e}"))
                .to_config()
                .unwrap();
            assert_eq!(cfg.pileup, exp_pileup, "argv={amp_cluster:?}");
            assert_eq!(
                cfg.unique_mode_alignment_enabled, exp_u,
                "argv={amp_cluster:?}"
            );
            assert_eq!(
                cfg.amplicon.as_deref(),
                Some(exp_amp),
                "argv={amp_cluster:?}"
            );
        }

        // -pR recognized by scan_has_region so trailing `-z 0` consumes `"0"` without BED
        let pr_cfg = Cli::try_parse_from(normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "a.bam",
            "-pR",
            "chr1:1-100",
            "-z",
            "0",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert!(pr_cfg.pileup);
        assert_eq!(pr_cfg.region.as_deref(), Some("chr1:1-100"));
        assert!(!pr_cfg.zero_based);
        assert!(pr_cfg.zero_based_defined);

        // Mandatory-value short flags (-Npz, -bpz.bam, -F0x504, -dpz), option values (-N -pz),
        // and `-- -pz` must NEVER be split by the boolean-prefix cluster splitter.
        let non_split_cfg = Cli::try_parse_from(normalize_args(&args(&[
            "-G", "r.fa", "-bpz.bam", "-Npz", "-F0x504", "-dpz", "--", "-pz",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(non_split_cfg.bam, "pz.bam");
        assert_eq!(non_split_cfg.sample, "pz");
        assert_eq!(non_split_cfg.flag_filter, 0x504);
        assert_eq!(non_split_cfg.delimiter, "pz");
        assert_eq!(non_split_cfg.bed.as_deref(), Some("-pz"));
        assert!(!non_split_cfg.pileup);
        assert!(!non_split_cfg.zero_based_defined);

        let opt_val_pz = Cli::try_parse_from(normalize_args(&args(&[
            "-G", "r.fa", "-b", "a.bam", "-N", "-pz",
        ])))
        .unwrap()
        .to_config()
        .unwrap();
        assert_eq!(opt_val_pz.sample, "-pz");
        assert!(!opt_val_pz.pileup);

        let num_norm = normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "a.bam",
            "-mfreq",
            "0.2",
            "-nmfreq",
            "0.15",
            "--insert-std",
            "80",
            "--insert-size",
            "250",
            "--ref-extension",
            "1500",
            "--crispr",
            "1000",
        ]));
        check_unsupported(&num_norm).unwrap();
        let num_cfg = Cli::try_parse_from(&num_norm).unwrap().to_config().unwrap();
        assert_eq!(num_cfg.mfreq, 0.2);
        assert_eq!(num_cfg.nmfreq, 0.15);
        assert_eq!(num_cfg.sv_flanking, 80);
        assert_eq!(num_cfg.sv_read_len, 250);
        assert_eq!(num_cfg.max_region_len, 1500);
        assert_eq!(num_cfg.crispr_cutting_site, Some(1000));

        let q_err = Cli::try_parse_from(normalize_args(&args(&["-?"]))).unwrap_err();
        assert_eq!(q_err.kind(), ErrorKind::DisplayHelp);

        let help_text = Cli::command().render_help().to_string();
        assert!(
            help_text.contains("--index"),
            "--index must be a visible alias in --help: {help_text}"
        );
        assert!(
            !help_text.contains("hg19/hg38/mm10"),
            "-G help must not advertise unimplemented hg19/hg38/mm10 shortcuts: {help_text}"
        );

        for (k_argv, expected_realign) in [
            (&["-k"][..], 1u8),
            (&["-k", "0"][..], 0u8),
            (&["-k", "-1"][..], 0u8),
            (&["-k", "300"][..], 0u8),
        ] {
            let mut full = vec!["-G", "r.fa", "-b", "a.bam"];
            full.extend_from_slice(k_argv);
            let norm = normalize_args(&args(&full));
            check_unsupported(&norm).unwrap();
            let cfg = Cli::try_parse_from(&norm)
                .unwrap_or_else(|e| panic!("failed on {k_argv:?} (norm={norm:?}): {e}"))
                .to_config()
                .unwrap();
            assert_eq!(cfg.realign, expected_realign, "argv={k_argv:?}");
        }
    }

    #[test]
    fn bed_column_zero_rejected_without_region() {
        for flag in ["-c", "-S", "-E", "-s", "-e"] {
            let norm = normalize_args(&args(&["-G", "r.fa", "-b", "a.bam", flag, "0", "r.bed"]));
            let cli = Cli::try_parse_from(&norm).unwrap();
            let err = cli
                .to_config()
                .expect_err(&format!("{flag} 0 without -R must fail"));
            assert!(
                matches!(err, vardictrust_core::types::VardictError::Parse(_)),
                "expected Parse error for {flag} 0, got {err:?}"
            );
        }

        let r_norm = normalize_args(&args(&[
            "-G",
            "r.fa",
            "-b",
            "a.bam",
            "-R",
            "chr1:100-200",
            "-c",
            "0",
        ]));
        assert!(
            Cli::try_parse_from(&r_norm).unwrap().to_config().is_ok(),
            "-R chr1:100-200 -c 0 must succeed"
        );

        let g_norm = normalize_args(&args(&["-G", "r.fa", "-b", "a.bam", "-g", "0", "r.bed"]));
        let g_cfg = Cli::try_parse_from(&g_norm).unwrap().to_config().unwrap();
        assert!(g_cfg.col_gene.is_none());
    }

    #[test]
    fn all_11_var2vcf_threshold_flags_and_typo_rejection() {
        use vardictrust_core::config::VcfThreshold;

        let parse_pipeline = |raw_flags: &[&str]| -> vardictrust_core::config::Config {
            let norm = normalize_args(&args(raw_flags));
            check_unsupported(&norm).expect("check_unsupported must pass");
            let cli = Cli::try_parse_from(&norm)
                .unwrap_or_else(|e| panic!("Cli::try_parse_from failed on {norm:?}: {e}"));
            cli.to_config().expect("to_config must succeed")
        };

        // 1. Defaults: all 11 fields are None when omitted
        let def = parse_pipeline(&["-G", "ref.fa", "-b", "in.bam", "-R", "chr1:1-100"]);
        assert!(def.vcf_af.is_none());
        assert!(def.vcf_min_depth.is_none());
        assert!(def.vcf_min_vd.is_none());
        assert!(def.vcf_qmean.is_none());
        assert!(def.vcf_mapq.is_none());
        assert!(def.vcf_pmean.is_none());
        assert!(def.vcf_sn.is_none());
        assert!(def.vcf_max_nm.is_none());
        assert!(def.vcf_gt_freq.is_none());
        assert!(def.vcf_msi.is_none());
        assert!(def.vcf_cluster_dist.is_none());

        // 2. Explicit values + first-wins + raw token preservation + BED -z safety
        let custom = parse_pipeline(&[
            "-G",
            "ref.fa",
            "-b",
            "in.bam",
            "--vcf-af",
            "0.050",
            "--vcf-af",
            "0.25",
            "--vcf-min-depth",
            "0",
            "--vcf-min-depth",
            "9",
            "--vcf-min-vd=0",
            "--vcf-min-vd=5",
            "--vcf-qmean",
            "20.0",
            "--vcf-mapq",
            "-0",
            "--vcf-pmean",
            "7.5",
            "--vcf-sn",
            "2.00",
            "--vcf-max-nm",
            "4.5",
            "--vcf-gt-freq",
            "+0.00",
            "--vcf-msi",
            "10.5",
            "--vcf-cluster-dist",
            "15",
            "-z",
            "regions.bed",
        ]);
        assert_eq!(
            custom.vcf_af,
            Some(VcfThreshold::with_raw(0.05, "0.050".to_string()))
        );
        assert_eq!(
            custom.vcf_min_depth,
            Some(VcfThreshold::with_raw(0, "0".to_string()))
        );
        assert_eq!(
            custom.vcf_min_vd,
            Some(VcfThreshold::with_raw(0, "0".to_string()))
        );
        assert_eq!(
            custom.vcf_qmean,
            Some(VcfThreshold::with_raw(20.0, "20.0".to_string()))
        );
        assert_eq!(
            custom.vcf_mapq,
            Some(VcfThreshold::with_raw(0.0, "0".to_string()))
        );
        assert_eq!(
            custom.vcf_pmean,
            Some(VcfThreshold::with_raw(7.5, "7.5".to_string()))
        );
        assert_eq!(
            custom.vcf_sn,
            Some(VcfThreshold::with_raw(2.0, "2.00".to_string()))
        );
        assert_eq!(
            custom.vcf_max_nm,
            Some(VcfThreshold::with_raw(4.5, "4.5".to_string()))
        );
        assert_eq!(custom.vcf_gt_freq, Some(0.0));
        assert!(custom.vcf_gt_freq.unwrap().is_sign_positive());
        assert_eq!(
            custom.vcf_msi,
            Some(VcfThreshold::with_raw(10.5, "10.5".to_string()))
        );
        assert_eq!(
            custom.vcf_cluster_dist,
            Some(VcfThreshold::with_raw(15, "15".to_string()))
        );
        assert!(custom.zero_based);
        assert_eq!(custom.bed.as_deref(), Some("regions.bed"));

        // 2b. Leading-zero / `+`-prefixed / `-0` integer CLI tokens preserved verbatim for Perl dualvar labels
        let raw_int_tokens = parse_pipeline(&[
            "-G",
            "ref.fa",
            "-b",
            "in.bam",
            "--vcf-min-depth",
            "03",
            "--vcf-min-vd",
            "+02",
            "--vcf-cluster-dist",
            "010",
            "-R",
            "chr1:1-100",
        ]);
        assert_eq!(
            raw_int_tokens.vcf_min_depth,
            Some(VcfThreshold::with_raw(3, "03".to_string()))
        );
        assert_eq!(
            raw_int_tokens.vcf_min_vd,
            Some(VcfThreshold::with_raw(2, "+02".to_string()))
        );
        assert_eq!(
            raw_int_tokens.vcf_cluster_dist,
            Some(VcfThreshold::with_raw(10, "010".to_string()))
        );

        // 3. Invalid / hostile value rejection across all 11 flags
        let over_64 = format!("0.{}", "0".repeat(63));
        for bad_flags in [
            vec!["--vcf-af", "-0.01", "-R", "chr1:1-100"],
            vec!["--vcf-af", "-1e-400", "-R", "chr1:1-100"],
            vec!["--vcf-af", over_64.as_str(), "-R", "chr1:1-100"],
            vec!["--vcf-min-depth", "-1", "-R", "chr1:1-100"],
            vec!["--vcf-min-depth", "-3", "-R", "chr1:1-100"],
            vec!["--vcf-min-depth", "2.5", "-R", "chr1:1-100"],
            vec!["--vcf-min-vd", "-1", "-R", "chr1:1-100"],
            vec!["--vcf-qmean", "NaN", "-R", "chr1:1-100"],
            vec!["--vcf-mapq", "-1", "-R", "chr1:1-100"],
            vec!["--vcf-pmean", "-3", "-R", "chr1:1-100"],
            vec!["--vcf-sn", "inf", "-R", "chr1:1-100"],
            vec!["--vcf-max-nm", "5.25f", "-R", "chr1:1-100"],
            vec!["--vcf-gt-freq", "1.01", "-R", "chr1:1-100"],
            vec!["--vcf-gt-freq", "-0.1", "-R", "chr1:1-100"],
            vec!["--vcf-msi", "", "-R", "chr1:1-100"],
            vec!["--vcf-cluster-dist", "-1", "-R", "chr1:1-100"],
            vec!["--vcf-cluster-dist", "-R", "chr1:1-100"],
        ] {
            let mut full = vec!["-G", "ref.fa", "-b", "in.bam"];
            full.extend_from_slice(&bad_flags);
            let norm = normalize_args(&args(&full));
            assert!(
                Cli::try_parse_from(&norm).is_err(),
                "expected failure for {bad_flags:?} (norm={norm:?})"
            );
        }

        // 4. Separate-value single-dash typo rejection for all 11 new flags
        for single_dash in [
            "-vcf-af",
            "-vcf-min-depth",
            "-vcf-min-vd",
            "-vcf-qmean",
            "-vcf-mapq",
            "-vcf-pmean",
            "-vcf-sn",
            "-vcf-max-nm",
            "-vcf-gt-freq",
            "-vcf-msi",
            "-vcf-cluster-dist",
        ] {
            let norm = normalize_args(&args(&[
                "-G",
                "r.fa",
                "-b",
                "a.bam",
                single_dash,
                "1",
                "-R",
                "c:1-10",
            ]));
            let err = check_unsupported(&norm).expect_err(&format!(
                "expected Parse error for single-dash typo {single_dash}"
            ));
            let expected = format!("did you mean '-{single_dash}'?");
            assert!(
                err.to_string().contains(&expected),
                "missing '{expected}' in {err}"
            );
        }
    }
}
