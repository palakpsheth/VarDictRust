//! Zero-dependency shell completion (`bash`, `zsh`, `fish`) and `roff` `vardictrust(1)` manpage
//! generators backed directly by `clap::CommandFactory` introspection on [`super::Cli`].

use std::io::{self, Write};

use clap::CommandFactory;

use super::{Cli, CompletionShell};

const HELP_HEADINGS: [&str; 7] = [
    "Core Inputs & Outputs",
    "VCF Post-Filtering (var2vcf)",
    "Coordinate & BED Column Mapping",
    "Variant Calling & Quality Thresholds",
    "Modality Flags: Somatic & Amplicon",
    "Structural Variants & Complex Indels",
    "Performance & Concurrency",
];

/// Escape arbitrary text for `roff` (`man(7)`) output (`AR-18`):
/// - `\` -> `\\`
/// - `-` -> `\-`
/// - `—` (`\u{2014}`) -> `\(em` and `–` (`\u{2013}`) -> `\(en`
/// - line-leading `.` or `'` -> `\&.` or `\&'` (prevents accidental macro invocation)
pub(crate) fn escape_roff(s: &str) -> String {
    let mut out = String::with_capacity(s.len().saturating_add(16));
    let mut at_line_start = true;
    for ch in s.chars() {
        if at_line_start && (ch == '.' || ch == '\'') {
            out.push_str("\\&");
            out.push(ch);
            at_line_start = false;
            continue;
        }
        match ch {
            '\\' => {
                out.push_str("\\\\");
                at_line_start = false;
            }
            '-' => {
                out.push_str("\\-");
                at_line_start = false;
            }
            '\u{2014}' => {
                out.push_str("\\(em");
                at_line_start = false;
            }
            '\u{2013}' => {
                out.push_str("\\(en");
                at_line_start = false;
            }
            '\n' => {
                out.push('\n');
                at_line_start = true;
            }
            _ => {
                out.push(ch);
                at_line_start = false;
            }
        }
    }
    out
}

fn single_line_summary(help: Option<&clap::builder::StyledStr>) -> String {
    let Some(styled) = help else {
        return String::new();
    };
    let raw = styled.to_string();
    let mut out = String::with_capacity(raw.len());
    let mut prev_space = false;
    for ch in raw.chars() {
        if ch.is_whitespace() {
            if !out.is_empty() && !prev_space {
                out.push(' ');
                prev_space = true;
            }
        } else {
            out.push(ch);
            prev_space = false;
        }
    }
    if out.ends_with(' ') {
        out.pop();
    }
    out
}

fn arg_takes_values(arg: &clap::Arg) -> bool {
    arg.get_action().takes_values()
}

fn arg_min_values(arg: &clap::Arg) -> usize {
    if !arg_takes_values(arg) {
        return 0;
    }
    arg.get_num_args().map_or(1, |r| r.min_values())
}

fn arg_value_placeholder(arg: &clap::Arg) -> String {
    if let Some(names) = arg.get_value_names() {
        if let Some(first) = names.first() {
            return first.as_str().to_string();
        }
    }
    arg.get_id().as_str().to_ascii_uppercase()
}

fn collect_all_flag_tokens(cmd: &clap::Command) -> Vec<String> {
    let mut flags = Vec::new();
    for arg in cmd.get_arguments() {
        if arg.is_positional() {
            continue;
        }
        if let Some(shorts) = arg.get_short_and_visible_aliases() {
            for s in shorts {
                if s != '?' {
                    flags.push(format!("-{s}"));
                }
            }
        }
        if let Some(longs) = arg.get_long_and_visible_aliases() {
            for l in longs {
                flags.push(format!("--{l}"));
            }
        }
        if arg.get_id().as_str() == "debug_print" {
            flags.push("--defaultPrinter".to_string());
            flags.push("-defaultPrinter".to_string());
        }
    }
    flags.push("-UN".to_string());
    flags.push("-th".to_string());
    flags
}

fn write_bash_completions<W: Write>(cmd: &clap::Command, out: &mut W) -> io::Result<()> {
    let all_flags = collect_all_flag_tokens(cmd).join(" ");
    let mut mandatory_non_file_flags = Vec::new();
    for arg in cmd.get_arguments() {
        if arg.is_positional() || arg_min_values(arg) == 0 {
            continue;
        }
        let id = arg.get_id().as_str();
        if matches!(
            id,
            "completions"
                | "reference"
                | "bam"
                | "output"
                | "stats_json"
                | "compat"
                | "validation_stringency"
                | "debug_print"
        ) {
            continue;
        }
        if let Some(shorts) = arg.get_short_and_visible_aliases() {
            for s in shorts {
                if s != '?' {
                    mandatory_non_file_flags.push(format!("-{s}"));
                }
            }
        }
        if let Some(longs) = arg.get_long_and_visible_aliases() {
            for l in longs {
                mandatory_non_file_flags.push(format!("--{l}"));
            }
        }
    }
    let scalar_case_pat = mandatory_non_file_flags.join("|");

    writeln!(out, "# Bash completion for vardictrust")?;
    writeln!(out, "_vardictrust() {{")?;
    writeln!(out, "    local cur prev opts")?;
    writeln!(out, "    COMPREPLY=()")?;
    writeln!(out, "    cur=\"${{COMP_WORDS[COMP_CWORD]}}\"")?;
    writeln!(out, "    prev=\"\"")?;
    writeln!(out, "    if [[ ${{COMP_CWORD}} -gt 0 ]]; then")?;
    writeln!(out, "        prev=\"${{COMP_WORDS[COMP_CWORD-1]}}\"")?;
    writeln!(out, "    fi")?;
    writeln!(out, "    opts=\"{all_flags}\"")?;
    writeln!(out)?;
    writeln!(out, "    case \"${{prev}}\" in")?;
    writeln!(out, "        --completions)")?;
    writeln!(
        out,
        "            COMPREPLY=( $(compgen -W \"bash zsh fish\" -- \"${{cur}}\") )"
    )?;
    writeln!(out, "            return 0")?;
    writeln!(out, "            ;;")?;
    writeln!(
        out,
        "        -G|--reference|-b|--bam|--output|--stats-json)"
    )?;
    writeln!(
        out,
        "            COMPREPLY=( $(compgen -f -- \"${{cur}}\") )"
    )?;
    writeln!(out, "            return 0")?;
    writeln!(out, "            ;;")?;
    writeln!(out, "        --compat)")?;
    writeln!(
        out,
        "            COMPREPLY=( $(compgen -W \"vardictjava-1.8.3 1.8.3\" -- \"${{cur}}\") )"
    )?;
    writeln!(out, "            return 0")?;
    writeln!(out, "            ;;")?;
    writeln!(out, "        --VS|-VS)")?;
    writeln!(
        out,
        "            COMPREPLY=( $(compgen -W \"SILENT LENIENT\" -- \"${{cur}}\") )"
    )?;
    writeln!(out, "            return 0")?;
    writeln!(out, "            ;;")?;
    writeln!(
        out,
        "        --DP|-DP|--default-printer|-default-printer|--defaultPrinter|-defaultPrinter)"
    )?;
    writeln!(
        out,
        "            COMPREPLY=( $(compgen -W \"OUT ERR\" -- \"${{cur}}\") )"
    )?;
    writeln!(out, "            return 0")?;
    writeln!(out, "            ;;")?;
    if !scalar_case_pat.is_empty() {
        writeln!(out, "        {scalar_case_pat})")?;
        writeln!(out, "            return 0")?;
        writeln!(out, "            ;;")?;
    }
    writeln!(out, "    esac")?;
    writeln!(out)?;
    writeln!(out, "    if [[ \"${{cur}}\" == -* ]]; then")?;
    writeln!(
        out,
        "        COMPREPLY=( $(compgen -W \"${{opts}}\" -- \"${{cur}}\") )"
    )?;
    writeln!(out, "        return 0")?;
    writeln!(out, "    fi")?;
    writeln!(out)?;
    writeln!(out, "    COMPREPLY=( $(compgen -f -- \"${{cur}}\") )")?;
    writeln!(out, "    return 0")?;
    writeln!(out, "}}")?;
    writeln!(
        out,
        "complete -F _vardictrust -o bashdefault -o default vardictrust"
    )?;
    Ok(())
}

fn escape_zsh_desc(s: &str) -> String {
    let mut out = String::with_capacity(s.len().saturating_add(8));
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '[' => out.push_str("\\["),
            ']' => out.push_str("\\]"),
            ':' => out.push_str("\\:"),
            '\'' => out.push_str("'\\''"),
            _ => out.push(ch),
        }
    }
    out
}

fn zsh_value_spec(arg: &clap::Arg) -> String {
    if !arg_takes_values(arg) {
        return String::new();
    }
    let id = arg.get_id().as_str();
    let placeholder = escape_zsh_desc(&arg_value_placeholder(arg));
    let action = match id {
        "completions" => "(bash zsh fish)",
        "reference" | "bam" | "output" | "stats_json" => "_files",
        "compat" => "(vardictjava-1.8.3 1.8.3)",
        "validation_stringency" => "(SILENT LENIENT)",
        "debug_print" => "(OUT ERR)",
        _ => " ",
    };
    if arg_min_values(arg) == 0 {
        format!("::{placeholder}:{action}")
    } else {
        format!(":{placeholder}:{action}")
    }
}

fn write_zsh_completions<W: Write>(cmd: &clap::Command, out: &mut W) -> io::Result<()> {
    writeln!(out, "#compdef vardictrust")?;
    writeln!(out)?;
    writeln!(out, "_vardictrust() {{")?;
    writeln!(out, "    _arguments -s -S \\")?;

    for arg in cmd.get_arguments() {
        if arg.is_positional() {
            continue;
        }
        let desc = escape_zsh_desc(&single_line_summary(arg.get_help()));
        let val_spec = zsh_value_spec(arg);
        let mut group = Vec::new();
        if let Some(shorts) = arg.get_short_and_visible_aliases() {
            for s in shorts {
                if s != '?' {
                    group.push(format!("-{s}"));
                }
            }
        }
        if let Some(longs) = arg.get_long_and_visible_aliases() {
            for l in longs {
                group.push(format!("--{l}"));
            }
        }
        if arg.get_id().as_str() == "debug_print" {
            group.push("--defaultPrinter".to_string());
        }
        let excl = if group.len() > 1 {
            format!("({})", group.join(" "))
        } else {
            String::new()
        };
        let is_opt_val = arg_takes_values(arg) && arg_min_values(arg) == 0;
        for flag in &group {
            let opt_eq = if is_opt_val {
                if flag.starts_with("--") {
                    "=-"
                } else {
                    "-"
                }
            } else {
                ""
            };
            writeln!(out, "        '{excl}{flag}{opt_eq}[{desc}]{val_spec}' \\")?;
        }
    }
    writeln!(
        out,
        "        '(-UN --UN)-UN[Indicate unique mode (first read in pair; VarDictJava single-dash alias)]' \\"
    )?;
    writeln!(
        out,
        "        '(-th --threads)-th[Threads (VarDictJava single-dash alias)]::THREADS: ' \\"
    )?;
    writeln!(out, "        '::BED file:_files'")?;
    writeln!(out, "}}")?;
    writeln!(out)?;
    writeln!(out, "_vardictrust \"$@\"")?;
    Ok(())
}

fn escape_fish_single_quoted(s: &str) -> String {
    let mut out = String::with_capacity(s.len().saturating_add(8));
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '\'' => out.push_str("\\'"),
            _ => out.push(ch),
        }
    }
    out
}

fn write_fish_completions<W: Write>(cmd: &clap::Command, out: &mut W) -> io::Result<()> {
    let all_flags = collect_all_flag_tokens(cmd).join(" ");
    writeln!(out, "# Fish completions for vardictrust")?;
    writeln!(out, "# Supported options: {all_flags}")?;
    for arg in cmd.get_arguments() {
        if arg.is_positional() {
            continue;
        }
        let mut line = String::from("complete -c vardictrust");
        if let Some(shorts) = arg.get_short_and_visible_aliases() {
            for s in shorts {
                if s != '?' {
                    line.push_str(&format!(" -s {s}"));
                }
            }
        }
        if let Some(longs) = arg.get_long_and_visible_aliases() {
            for l in longs {
                line.push_str(&format!(" -l {l}"));
            }
        }
        if arg.get_id().as_str() == "debug_print" {
            line.push_str(" -l defaultPrinter -o defaultPrinter");
        }
        if arg_min_values(arg) > 0 {
            line.push_str(" -r");
        }
        match arg.get_id().as_str() {
            "completions" => line.push_str(" -f -a 'bash zsh fish'"),
            "compat" => line.push_str(" -f -a 'vardictjava-1.8.3 1.8.3'"),
            "validation_stringency" => line.push_str(" -f -a 'SILENT LENIENT'"),
            "debug_print" => line.push_str(" -f -a 'OUT ERR'"),
            "reference" | "bam" | "output" | "stats_json" => line.push_str(" -F"),
            _ if arg_min_values(arg) > 0 => line.push_str(" -f"),
            _ => {}
        }
        let desc = escape_fish_single_quoted(&single_line_summary(arg.get_help()));
        if !desc.is_empty() {
            line.push_str(&format!(" -d '{desc}'"));
        }
        writeln!(out, "{line}")?;
    }
    writeln!(
        out,
        "complete -c vardictrust -o UN -d 'Indicate unique mode (first read in pair; VarDictJava single-dash alias for --UN)'"
    )?;
    writeln!(
        out,
        "complete -c vardictrust -o th -d 'Threads (VarDictJava single-dash alias for --threads)'"
    )?;
    Ok(())
}

/// Write shell completion script for `shell` (`bash`, `zsh`, or `fish`) to `out`.
pub(crate) fn write_completions<W: Write>(shell: CompletionShell, out: &mut W) -> io::Result<()> {
    let cmd = Cli::command();
    match shell {
        CompletionShell::Bash => write_bash_completions(&cmd, out),
        CompletionShell::Zsh => write_zsh_completions(&cmd, out),
        CompletionShell::Fish => write_fish_completions(&cmd, out),
    }
}

fn format_man_arg_synopsis(arg: &clap::Arg) -> String {
    let mut parts = Vec::new();
    if let Some(shorts) = arg.get_short_and_visible_aliases() {
        for s in shorts {
            parts.push(format!("-{s}"));
        }
    }
    if arg.get_id().as_str() == "_help" {
        parts.push("-?".to_string());
    }
    if arg.get_id().as_str() == "unique_mode_second_in_pair_enabled" {
        parts.push("-UN".to_string());
    }
    if let Some(longs) = arg.get_long_and_visible_aliases() {
        for l in longs {
            parts.push(format!("--{l}"));
        }
    }
    if arg.get_id().as_str() == "threads" {
        parts.insert(0, "-th".to_string());
    }
    let header = if arg_takes_values(arg) {
        let placeholder = arg_value_placeholder(arg);
        if arg_min_values(arg) == 0 {
            let formatted: Vec<String> = parts
                .into_iter()
                .map(|flag| {
                    if flag.starts_with("--") && arg.get_id().as_str() != "threads" {
                        format!("{flag}[={placeholder}]")
                    } else {
                        format!("{flag} [{placeholder}]")
                    }
                })
                .collect();
            formatted.join(", ")
        } else {
            format!("{} <{placeholder}>", parts.join(", "))
        }
    } else {
        parts.join(", ")
    };
    escape_roff(&header)
}

/// Write a standalone `roff` `vardictrust(1)` manpage to `out`.
pub(crate) fn write_manpage<W: Write>(out: &mut W) -> io::Result<()> {
    let cmd = Cli::command();
    let version = escape_roff(env!("CARGO_PKG_VERSION"));

    writeln!(
        out,
        ".TH VARDICTRUST 1 \"\" \"vardictrust {version}\" \"User Commands\""
    )?;
    writeln!(out, ".SH NAME")?;
    writeln!(
        out,
        "vardictrust \\- single\\-binary somatic and germline variant caller (VarDictJava v1.8.3\\-compatible)"
    )?;
    writeln!(out, ".SH SYNOPSIS")?;
    writeln!(
        out,
        ".B vardictrust\n.B \\-G\n.I FASTA\n.B \\-b\n.I BAM\n.RI [ OPTIONS ]\n.RI [ BED ]"
    )?;
    writeln!(out, ".br")?;
    writeln!(
        out,
        ".B vardictrust\n.B \\-G\n.I FASTA\n.B \\-b\n.I \"TUMOR.bam|NORMAL.bam\"\n.B \\-R\n.I CHR:START\\-END[:GENE]\n.RI [ OPTIONS ]"
    )?;
    writeln!(out, ".br")?;
    writeln!(
        out,
        ".B vardictrust\n.B \\-\\-completions\n.I bash|zsh|fish"
    )?;
    writeln!(out, ".br")?;
    writeln!(out, ".B vardictrust\n.B \\-\\-manpage")?;

    writeln!(out, ".SH DESCRIPTION")?;
    writeln!(
        out,
        "{}",
        escape_roff(
            "VarDictRust consolidates the three-stage VarDictJava v1.8.3 + R (teststrandbias.R / testsomatic.R) + Perl (var2vcf_valid.pl / var2vcf_paired.pl) variant-calling pipeline into a single statically-compiled binary with bounded-memory WGS tiling, native BAM/CRAM and BGZF FASTA/BED support, and direct BGZF + Tabix (.tbi) / CSI (.csi) VCF 4.2 (VCFv4.2) emission."
        )
    )?;
    writeln!(out, ".PP")?;
    writeln!(
        out,
        "{}",
        escape_roff(
            "Important CLI note: -h / --header prints the intermediate .var column header row (matching VarDictJava CmdParser.java); use -H, -?, or --help to display CLI help."
        )
    )?;

    writeln!(out, ".SH OPTIONS")?;
    for heading in HELP_HEADINGS {
        writeln!(out, ".SS {}", escape_roff(heading))?;
        for arg in cmd.get_arguments() {
            if arg.is_positional() || arg.get_help_heading() != Some(heading) {
                continue;
            }
            writeln!(out, ".TP")?;
            writeln!(out, ".B {}", format_man_arg_synopsis(arg))?;
            let help = single_line_summary(arg.get_help());
            writeln!(out, "{}", escape_roff(&help))?;
        }
    }

    writeln!(out, ".SS {}", escape_roff("Positional Arguments"))?;
    for arg in cmd.get_arguments() {
        if !arg.is_positional() {
            continue;
        }
        writeln!(out, ".TP")?;
        writeln!(
            out,
            ".B {}",
            escape_roff(&format!("[{}]", arg_value_placeholder(arg)))
        )?;
        let help = single_line_summary(arg.get_help());
        writeln!(out, "{}", escape_roff(&help))?;
    }

    writeln!(out, ".SH EXAMPLES")?;
    writeln!(out, ".TP")?;
    writeln!(
        out,
        "{}",
        escape_roff("Germline calling on a standard 0-based 3/4-column BED file:")
    )?;
    writeln!(
        out,
        ".nf\n{}\n.fi",
        escape_roff(
            "vardictrust -G ref.fa -b sample.bam -N SAMPLE -f 0.01 -c 1 -S 2 -E 3 -g 4 -z targets.bed > variants.vcf"
        )
    )?;
    writeln!(out, ".TP")?;
    writeln!(
        out,
        "{}",
        escape_roff("Paired tumor/normal somatic calling on a 0-based BED file with BGZF output and Tabix index:")
    )?;
    writeln!(
        out,
        ".nf\n{}\n.fi",
        escape_roff(
            "vardictrust -G ref.fa -b \"tumor.bam|normal.bam\" -N \"TUMOR|NORMAL\" -f 0.01 -c 1 -S 2 -E 3 -z targets.bed --output somatic.vcf.gz --tbi"
        )
    )?;

    writeln!(out, ".SH EXIT STATUS")?;
    writeln!(
        out,
        "{}",
        escape_roff(
            "0 on success; non-zero on invalid CLI options, unsupported legacy flags, missing index/input files, or I/O errors."
        )
    )?;

    writeln!(out, ".SH AUTHORS")?;
    writeln!(
        out,
        "{}",
        escape_roff(
            "Original algorithm and VarDictJava v1.8.3 reference: Zhongwu Lai et al., AstraZeneca-NGS (Nucleic Acids Res. 2016;44(11):e108, doi:10.1093/nar/gkw227). VarDictRust port built under the rewrites.bio manifesto."
        )
    )?;
    Ok(())
}
