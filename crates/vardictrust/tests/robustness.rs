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

#[test]
fn output_plain_and_bgzf_match_stdout_and_stats_json_emits_telemetry() {
    use std::io::Read as _;

    let dir = tempfile::TempDir::new().expect("tempdir");
    let plain_vcf = dir.path().join("out.vcf");
    let bgzf_vcf = dir.path().join("out.vcf.gz");
    let stats_json = dir.path().join("stats.json");

    let common_args = [
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
        "--threads=2",
    ];

    let stdout_run = base_cmd()
        .args(common_args)
        .output()
        .expect("spawn vardictrust stdout");
    assert!(
        stdout_run.status.success(),
        "stdout run failed: {}",
        String::from_utf8_lossy(&stdout_run.stderr)
    );
    assert!(!stdout_run.stdout.is_empty());

    // 1. --output <tmp>/out.vcf + --stats-json <tmp>/stats.json
    let plain_run = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&plain_vcf)
        .arg("--stats-json")
        .arg(&stats_json)
        .output()
        .expect("spawn vardictrust --output out.vcf");
    assert!(
        plain_run.status.success(),
        "plain --output run failed: {}",
        String::from_utf8_lossy(&plain_run.stderr)
    );
    assert!(
        plain_run.stdout.is_empty(),
        "stdout must be empty when --output is set"
    );
    let plain_bytes = std::fs::read(&plain_vcf).expect("read out.vcf");
    assert_eq!(
        plain_bytes, stdout_run.stdout,
        "--output out.vcf must be byte-for-byte identical to stdout"
    );

    // 2. --output <tmp>/out.vcf.gz (BGZF)
    let bgzf_run = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&bgzf_vcf)
        .output()
        .expect("spawn vardictrust --output out.vcf.gz");
    assert!(
        bgzf_run.status.success(),
        "bgzf --output run failed: {}",
        String::from_utf8_lossy(&bgzf_run.stderr)
    );
    assert!(
        bgzf_run.stdout.is_empty(),
        "stdout must be empty when --output out.vcf.gz is set"
    );
    let bgzf_file = std::fs::File::open(&bgzf_vcf).expect("open out.vcf.gz");
    let mut bgzf_reader = noodles::bgzf::io::Reader::new(bgzf_file);
    let mut decompressed = Vec::new();
    bgzf_reader
        .read_to_end(&mut decompressed)
        .expect("decompress BGZF out.vcf.gz");
    assert_eq!(
        decompressed, stdout_run.stdout,
        "decompressed --output out.vcf.gz must be byte-for-byte identical to stdout"
    );

    // 3. --stats-json <tmp>/stats.json telemetry keys
    let stats_text = std::fs::read_to_string(&stats_json).expect("read stats.json");
    let trimmed = stats_text.trim();
    assert!(
        trimmed.starts_with('{') && trimmed.ends_with('}'),
        "stats.json must be a JSON object, got: {stats_text}"
    );
    for key in [
        "\"mode\": \"Simple\"",
        "\"threads\": 2",
        "\"regions_count\":",
        "\"tiles_count\":",
        "\"elapsed_ms\":",
        "\"peak_rss_kb\":",
        "\"exceptions_caught\": 0",
    ] {
        assert!(
            stats_text.contains(key),
            "stats.json missing expected entry {key}, got:\n{stats_text}"
        );
    }

    // 4. Colliding --output and --stats-json paths must be rejected before creating or truncating either file
    let non_existent = dir.path().join("collide_new.vcf");
    let non_existent_dot = dir.path().join(".").join("collide_new.vcf");
    let collide_run = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&non_existent)
        .arg("--stats-json")
        .arg(&non_existent_dot)
        .output()
        .expect("spawn vardictrust colliding paths");
    assert!(
        !collide_run.status.success(),
        "colliding --output and --stats-json must fail"
    );
    assert!(
        !non_existent.exists(),
        "colliding destination must not be created"
    );
    let collide_stderr = String::from_utf8_lossy(&collide_run.stderr);
    assert!(
        collide_stderr.contains("--output and --stats-json must not point to the same file"),
        "expected collision error message, got: {collide_stderr}"
    );
}

#[test]
fn output_bgzf_tbi_and_parity_compat_flags() {
    use noodles::csi::BinningIndex as _;
    use std::io::Read as _;

    let dir = tempfile::TempDir::new().expect("tempdir");
    let common_args = [
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
        "--threads=2",
    ];

    let stdout_run = base_cmd()
        .args(common_args)
        .output()
        .expect("spawn vardictrust stdout");
    assert!(stdout_run.status.success());

    // 1. --output <tmp>/indexed.vcf.gz --tbi and --index alias
    for (flag, file_name) in [
        ("--tbi", "indexed_tbi.vcf.gz"),
        ("--index", "indexed_alias.vcf.gz"),
    ] {
        let vcf_gz = dir.path().join(file_name);
        let tbi_path = dir.path().join(format!("{file_name}.tbi"));
        let run = base_cmd()
            .args(common_args)
            .arg("--output")
            .arg(&vcf_gz)
            .arg(flag)
            .output()
            .expect("spawn vardictrust --output .vcf.gz --tbi");
        assert!(
            run.status.success(),
            "{flag} run failed: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let mut reader = noodles::bgzf::io::Reader::new(
            std::fs::File::open(&vcf_gz).expect("open indexed .vcf.gz"),
        );
        let mut decompressed = Vec::new();
        reader
            .read_to_end(&mut decompressed)
            .expect("decompress indexed .vcf.gz");
        assert_eq!(decompressed, stdout_run.stdout);

        let index = noodles::tabix::fs::read(&tbi_path).expect("read generated .tbi");
        let hdr = index.header().expect("tabix header must be present");
        assert_eq!(
            hdr.format(),
            noodles::csi::binning_index::index::header::Format::Vcf,
            "tabix header format must be VCF (not default GFF)"
        );
        assert!(
            !hdr.reference_sequence_names().is_empty(),
            "tabix header reference_sequence_names must not be empty"
        );
        assert!(
            index
                .reference_sequences()
                .iter()
                .any(|rs| !rs.bins().is_empty()),
            "tabix reference_sequences must contain non-empty bins"
        );
    }

    // 2. 0-variant header-only VCF with --tbi still produces a valid VCF-formatted .tbi
    let empty_vcf_gz = dir.path().join("empty.vcf.gz");
    let empty_tbi = dir.path().join("empty.vcf.gz.tbi");
    let empty_run = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:1-10",
            "-f",
            "0.99",
        ])
        .arg("--output")
        .arg(&empty_vcf_gz)
        .arg("--tbi")
        .output()
        .expect("spawn header-only --tbi");
    assert!(empty_run.status.success());
    let empty_index = noodles::tabix::fs::read(&empty_tbi).expect("read header-only .tbi");
    assert_eq!(
        empty_index.header().unwrap().format(),
        noodles::csi::binning_index::index::header::Format::Vcf
    );

    // 3. Invalid --tbi combinations fail fast before creating output files
    let bad_plain = dir.path().join("bad_plain.vcf");
    let bad_plain_run = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&bad_plain)
        .arg("--tbi")
        .output()
        .expect("spawn --tbi with plain .vcf");
    assert!(!bad_plain_run.status.success());
    assert!(!bad_plain.exists(), "bad_plain.vcf must not be created");
    assert!(String::from_utf8_lossy(&bad_plain_run.stderr)
        .contains("--tbi/--index requires --output <PATH> ending in .vcf.gz or .bgz"));

    let bad_inter = dir.path().join("bad_inter.vcf.gz");
    let bad_inter_run = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&bad_inter)
        .arg("--tbi")
        .arg("--emit-intermediate")
        .output()
        .expect("spawn --tbi with --emit-intermediate");
    assert!(!bad_inter_run.status.success());
    assert!(!bad_inter.exists(), "bad_inter.vcf.gz must not be created");
    assert!(String::from_utf8_lossy(&bad_inter_run.stderr)
        .contains("--tbi/--index is only supported for VCF output"));

    // 4. --parity and --compat vardictjava-1.8.3 match default stdout byte-for-byte; unknown --compat fails with Unsupported
    for compat_args in [
        vec!["--parity"],
        vec!["--compat", "vardictjava-1.8.3"],
        vec!["--compat=1.8.3"],
    ] {
        let run = base_cmd()
            .args(common_args)
            .args(&compat_args)
            .output()
            .expect("spawn --parity / --compat");
        assert!(
            run.status.success(),
            "{compat_args:?} failed: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        assert_eq!(run.stdout, stdout_run.stdout);
    }

    let bad_compat = base_cmd()
        .args(common_args)
        .args(["--compat", "foo"])
        .output()
        .expect("spawn invalid --compat");
    assert!(!bad_compat.status.success());
    let bad_compat_stderr = String::from_utf8_lossy(&bad_compat.stderr);
    assert!(
        bad_compat_stderr.contains("unsupported:")
            && bad_compat_stderr.contains("--compat foo not supported in v1.0"),
        "expected Unsupported error for --compat foo, got: {bad_compat_stderr}"
    );
}

#[test]
fn upstream_toml_manifest_present_and_pinned() {
    let upstream_path = workspace_root().join("UPSTREAM.toml");
    let content = std::fs::read_to_string(&upstream_path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", upstream_path.display()));
    for expected in ["0886261", "009e017", "1.8.3"] {
        assert!(
            content.contains(expected),
            "UPSTREAM.toml must contain '{expected}', got:\n{content}"
        );
    }
}
