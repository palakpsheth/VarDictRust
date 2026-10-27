//! Process-level CLI smoke and robustness tests (`4` tests).
//!
//! Verifies end-to-end binary stream redirection (`--DP ERR` / `--defaultPrinter ERR`),
//! verbose stderr region debug output (`-y` / `--verbose`), non-zero exit codes on
//! unsupported/invalid CLI flags, and C vs `de_DE` locale independence on `micro.bam`.

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
    ] {
        let out = base_cmd()
            .args([
                "-G",
                "fixtures/micro/micro.ref.fa",
                "-b",
                "fixtures/micro/micro.bam",
                "-R",
                "chr1:100-110",
            ])
            .args(invalid_args)
            .output()
            .expect("spawn vardictrust");

        assert!(!out.status.success());
        assert_ne!(out.status.code(), Some(101), "died by panic");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains(expected_msg),
            "stderr should contain '{expected_msg}', got: {stderr}"
        );
    }
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
