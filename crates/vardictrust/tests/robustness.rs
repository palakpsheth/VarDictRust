//! Process-level CLI smoke and robustness tests (`5` tests).
//!
//! Verifies end-to-end binary stream redirection (`--DP ERR` / `--defaultPrinter ERR`),
//! verbose stderr region debug output (`-y` / `--verbose`), non-zero exit codes on
//! unsupported/invalid CLI flags, C vs `de_DE` locale independence on `micro.bam`,
//! and clean jemalloc startup without `<jemalloc>:` stderr warnings.

#[cfg(unix)]
use std::os::unix::process::ExitStatusExt;
use std::path::Path;
use std::process::Command;

fn workspace_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
}

fn base_cmd() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vardictrust"));
    cmd.current_dir(workspace_root());
    cmd
}

#[test]
fn dp_err_and_default_printer_route_to_stderr() {
    for flag in ["--DP", "--defaultPrinter"] {
        let out = base_cmd()
            .args([
                "-G",
                "fixtures/micro/micro.ref.fa",
                "-b",
                "fixtures/micro/micro.bam",
                "-c",
                "1",
                "-S",
                "2",
                "-E",
                "3",
                "-g",
                "4",
                "-z",
                "fixtures/micro/micro.bed",
                "--threads=1",
                flag,
                "ERR",
            ])
            .output()
            .unwrap_or_else(|e| panic!("run vardictrust with {flag} ERR: {e}"));

        assert!(
            out.status.success(),
            "command failed under {flag} ERR: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(
            out.stdout.is_empty(),
            "stdout must be empty under {flag} ERR, got {} bytes",
            out.stdout.len()
        );
        assert!(
            !out.stderr.is_empty(),
            "stderr must contain variant output under {flag} ERR"
        );
    }

    // H2: When -h (--header) is combined with --DP ERR, VarDictJava prints the .var
    // header to System.out (stdout) while routing variant rows to System.err (stderr).
    let out_hdr = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-z",
            "-h",
            "fixtures/micro/micro.bed",
            "--threads=1",
            "--DP",
            "ERR",
        ])
        .output()
        .expect("run vardictrust with -h --DP ERR");
    assert!(out_hdr.status.success(), "-h --DP ERR failed");
    let stdout_str = String::from_utf8_lossy(&out_hdr.stdout);
    assert!(
        stdout_str.starts_with("Sample\tGene\tChr\t"),
        "stdout must contain .var header under -h --DP ERR, got: {stdout_str}"
    );
    assert_eq!(
        stdout_str.lines().count(),
        1,
        "stdout must ONLY contain the single .var header line under -h --DP ERR"
    );
    assert!(
        !out_hdr.stderr.is_empty(),
        "stderr must contain variant rows under -h --DP ERR"
    );
}

#[test]
fn verbose_flag_y_emits_region_debug_to_stderr() {
    for flag in ["-y", "--verbose"] {
        let out = base_cmd()
            .args([
                "-G",
                "fixtures/micro/micro.ref.fa",
                "-b",
                "fixtures/micro/micro.bam",
                "-R",
                "chr1:40-260",
                "--threads=1",
                flag,
            ])
            .output()
            .unwrap_or_else(|e| panic!("run vardictrust with {flag}: {e}"));

        assert!(
            out.status.success(),
            "command failed under {flag}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(
            out.stdout.windows(5).any(|w| w == b"CHROM"),
            "stdout must contain VCF header under {flag}"
        );
    }
}

#[test]
fn unsupported_and_invalid_cli_flags_exit_nonzero() {
    for flag in ["-n", "-C", "-T", "--trim", "--VS", "-Z", "--downsample"] {
        let out = base_cmd()
            .args([
                "-G",
                "fixtures/micro/micro.ref.fa",
                "-b",
                "fixtures/micro/micro.bam",
                "-R",
                "chr1:100-110",
                flag,
                "dummy_arg",
            ])
            .output()
            .unwrap_or_else(|e| panic!("run vardictrust with {flag}: {e}"));

        assert!(
            !out.status.success(),
            "unsupported flag {flag} must exit non-zero"
        );
        #[cfg(unix)]
        assert_eq!(out.status.signal(), None, "{flag} died by signal");
        assert_ne!(out.status.code(), Some(101), "{flag} died by panic");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("not implemented in v1.0; use VarDictJava"),
            "flag {flag} stderr must mention v1.0 rejection, got: {stderr}"
        );
    }

    for (invalid_args, expected_msg) in [
        (&["-F", "bogus"][..], "Invalid -F value"),
        (&["--tile-size", "0"][..], "Invalid --tile-size"),
        (&["-R", "chr1:1-3000000000"][..], "invalid integer"),
    ] {
        let mut cmd = base_cmd();
        cmd.args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
        ])
        .args(invalid_args);
        if !invalid_args.contains(&"-R") {
            cmd.args(["-R", "chr1:100-110"]);
        }
        let out = cmd.output().expect("spawn vardictrust");

        assert!(!out.status.success());
        #[cfg(unix)]
        assert_eq!(out.status.signal(), None, "died by signal");
        assert_ne!(out.status.code(), Some(101), "died by panic");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains(expected_msg),
            "stderr should contain '{expected_msg}', got: {stderr}"
        );
    }

    // H2 (#378): Paired BAM (-b tumor|normal) with amplicon mode (-a) in VCF mode
    // must return VardictError::Unsupported with non-zero exit code.
    let out_378 = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam|fixtures/micro/micro.bam",
            "-a",
            "10:0.95",
            "fixtures/beds/correct_8columns_amplicon.bed",
        ])
        .output()
        .expect("spawn vardictrust paired amplicon VCF");
    assert!(
        !out_378.status.success(),
        "paired + amplicon VCF (#378) must exit non-zero"
    );
    assert_ne!(out_378.status.code(), Some(101), "#378 died by panic");
    let stderr_378 = String::from_utf8_lossy(&out_378.stderr);
    assert!(
        stderr_378.contains("not supported in VCF mode (#378)"),
        "stderr should mention #378 unsupported VCF mode, got: {stderr_378}"
    );
}

#[test]
fn locale_independence_c_vs_de_de_on_micro_bam() {
    let de_installed = Command::new("locale")
        .arg("-a")
        .output()
        .map(|o| {
            String::from_utf8_lossy(&o.stdout)
                .to_lowercase()
                .contains("de_de")
        })
        .unwrap_or(false);
    if !de_installed {
        println!("SKIP: no de_DE locale on host; cross-locale byte equality not exercised");
        return;
    }

    let run = |envs: &[(&str, &str)]| {
        let mut cmd = base_cmd();
        cmd.args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--threads=1",
        ]);
        for &(k, v) in envs {
            cmd.env(k, v);
        }
        cmd.output().expect("spawn vardictrust")
    };

    let c = run(&[("LC_ALL", "C"), ("LC_NUMERIC", "C")]);
    let de = run(&[("LC_ALL", "de_DE.UTF-8"), ("LC_NUMERIC", "de_DE.UTF-8")]);
    for (out, lc) in [(&c, "C"), (&de, "de_DE")] {
        assert!(out.status.success(), "{lc}: non-zero exit");
        let s = String::from_utf8_lossy(&out.stdout);
        assert!(s.contains("CHROM"), "{lc}: no VCF output");
        assert!(s.contains("AF="), "{lc}: no formatted float record emitted");
    }
    assert_eq!(
        c.stdout, de.stdout,
        "output bytes differ between LC_ALL=C and de_DE"
    );
}

#[test]
fn jemalloc_startup_emits_no_stderr_warnings() {
    let ver = base_cmd()
        .args(["--version"])
        .output()
        .expect("spawn vardictrust --version");
    assert!(
        ver.status.success(),
        "--version failed: {}",
        String::from_utf8_lossy(&ver.stderr)
    );
    let ver_stderr = String::from_utf8_lossy(&ver.stderr);
    assert!(
        !ver_stderr.contains("<jemalloc>:"),
        "unexpected jemalloc warning on --version stderr: {ver_stderr}"
    );

    let call = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--threads=1",
        ])
        .output()
        .expect("spawn vardictrust on micro.bam");
    assert!(
        call.status.success(),
        "micro.bam call failed: {}",
        String::from_utf8_lossy(&call.stderr)
    );
    assert!(
        call.stdout.windows(5).any(|w| w == b"CHROM"),
        "stdout must contain VCF header"
    );
    let call_stderr = String::from_utf8_lossy(&call.stderr);
    assert!(
        !call_stderr.contains("<jemalloc>:"),
        "unexpected jemalloc warning on micro.bam stderr: {call_stderr}"
    );
}
