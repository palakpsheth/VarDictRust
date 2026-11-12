//! Process-level CLI smoke and robustness tests (`13` tests).
//!
//! Verifies end-to-end binary stream redirection (`--DP ERR` / `--defaultPrinter ERR`),
//! verbose stderr region debug output (`-y` / `--verbose`), non-zero exit codes on
//! unsupported/invalid CLI flags, C vs `de_DE` locale independence on `micro.bam`,
//! clean jemalloc startup without `<jemalloc>:` stderr warnings, `--output` / `--tbi` /
//! `--stats-json` collision & BrokenPipe guards, hostile BAM records, and Commons-CLI
//! single-dash option normalization.

#[cfg(unix)]
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::Once;

fn workspace_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
}

fn data_root() -> PathBuf {
    let root = workspace_root();
    if root
        .join("VarDictJava/testdata/integrationtestcases")
        .is_dir()
    {
        return root.to_path_buf();
    }
    if let Ok(out) = Command::new("git")
        .args([
            "-C",
            root.to_str().unwrap_or("."),
            "rev-parse",
            "--path-format=absolute",
            "--git-common-dir",
        ])
        .output()
    {
        if out.status.success() {
            let common = String::from_utf8_lossy(&out.stdout).trim().to_string();
            if let Some(parent) = Path::new(&common).parent() {
                if parent
                    .join("VarDictJava/testdata/integrationtestcases")
                    .is_dir()
                {
                    return parent.to_path_buf();
                }
            }
        }
    }
    root.to_path_buf()
}

fn ensure_fastas_converted() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let droot = data_root();
        let fasta_dir = droot.join("VarDictJava/testdata/fastas");
        if !fasta_dir.join("GRCh37.fa").is_file()
            || !fasta_dir.join("GRCh37.fa.fai").is_file()
            || !fasta_dir.join("hg19.fa").is_file()
            || !fasta_dir.join("hg19.fa.fai").is_file()
            || !fasta_dir.join("hg38.fa").is_file()
            || !fasta_dir.join("hg38.fa.fai").is_file()
        {
            let script = workspace_root().join("harness/convert_fastas.py");
            let status = Command::new("python3")
                .arg(&script)
                .arg("--data-root")
                .arg(&droot)
                .current_dir(&droot)
                .status()
                .expect("run convert_fastas.py");
            assert!(status.success(), "convert_fastas.py failed");
        }
    });
}

fn write_indexed_bam(
    path: &Path,
    chr: &str,
    len: usize,
    recs: &[noodles::sam::alignment::record_buf::RecordBuf],
) -> std::io::Result<()> {
    use noodles::bam::bai;
    use noodles::core::Position;
    use noodles::csi::binning_index::index::reference_sequence::{bin::Chunk, index::LinearIndex};
    use noodles::csi::binning_index::Indexer;
    use noodles::sam::{
        alignment::io::Write as _,
        header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap},
        Header,
    };

    let rs = HeaderMap::<RefSeqMap>::builder()
        .set_length(len.try_into().expect("nonzero contig length"))
        .build()
        .expect("valid reference sequence header");
    let hd = Header::builder().add_reference_sequence(chr, rs).build();
    let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(path)?);
    let mut w = noodles::bam::io::Writer::from(bgzf);
    w.write_alignment_header(&hd)?;
    let mut indexer = Indexer::<LinearIndex>::new(14, 5).expect("valid BAI parameters");
    for rec in recs {
        let start = w.get_ref().virtual_position();
        w.write_alignment_record(&hd, rec)?;
        let end = w.get_ref().virtual_position();
        let s = rec.alignment_start().expect("aligned record");
        let e = Position::new(s.get() + rec.cigar().alignment_span().saturating_sub(1))
            .unwrap_or(Position::MIN);
        indexer
            .add_record(Some((0, s, e, true)), Chunk::new(start, end))
            .map_err(std::io::Error::other)?;
    }
    w.try_finish()?;
    let mut bai_path = PathBuf::from(path);
    bai_path.set_extension("bam.bai");
    bai::io::Writer::new(std::fs::File::create(bai_path)?).write_index(&indexer.build(1))
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
    for unsupported_args in [
        &["-C", "dummy_arg"][..],
        &["-Z", "dummy_arg"][..],
        &["--downsample", "dummy_arg"][..],
        &["-VS", "STRICT"][..],
        &["--VS", "STRICT"][..],
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
            .args(unsupported_args)
            .output()
            .unwrap_or_else(|e| panic!("run vardictrust with {unsupported_args:?}: {e}"));

        assert!(
            !out.status.success(),
            "unsupported flag {unsupported_args:?} must exit non-zero"
        );
        #[cfg(unix)]
        assert_eq!(
            out.status.signal(),
            None,
            "{unsupported_args:?} died by signal"
        );
        assert_ne!(
            out.status.code(),
            Some(101),
            "{unsupported_args:?} died by panic"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("not implemented in v1.0; use VarDictJava"),
            "flag {unsupported_args:?} stderr must mention v1.0 rejection, got: {stderr}"
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

    let bad_hdr_dp = dir.path().join("bad_hdr_dp.vcf.gz");
    let bad_hdr_dp_run = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&bad_hdr_dp)
        .arg("--tbi")
        .args(["-h", "--DP", "ERR"])
        .output()
        .expect("spawn --tbi with -h --DP ERR");
    assert!(!bad_hdr_dp_run.status.success());
    assert!(
        !bad_hdr_dp.exists(),
        "bad_hdr_dp.vcf.gz must not be created when -h --DP ERR is combined with --tbi"
    );
    assert!(String::from_utf8_lossy(&bad_hdr_dp_run.stderr)
        .contains("--tbi/--index is only supported for VCF output"));

    #[cfg(unix)]
    {
        let alias_vcf_gz = dir.path().join("alias_collision.vcf.gz");
        let alias_tbi = dir.path().join("alias_collision.vcf.gz.tbi");
        std::fs::write(&alias_vcf_gz, b"sentinel").expect("write sentinel");
        std::os::unix::fs::symlink(&alias_vcf_gz, &alias_tbi).expect("symlink .tbi to .vcf.gz");
        let alias_run = base_cmd()
            .args(common_args)
            .arg("--output")
            .arg(&alias_vcf_gz)
            .arg("--tbi")
            .output()
            .expect("spawn --tbi with colliding .tbi symlink");
        assert!(!alias_run.status.success());
        assert_eq!(
            std::fs::read(&alias_vcf_gz).expect("read sentinel after collision"),
            b"sentinel",
            "colliding .tbi symlink must neither truncate nor delete the target file"
        );
        assert!(
            String::from_utf8_lossy(&alias_run.stderr).contains("must not point to the same file")
        );
    }

    // 3b. --csi mutual exclusion with --tbi, .csi collision, /dev/null rejection,
    //     0-variant empty .csi roundtrip, and >512 Mbp (700 Mbp contig, pos = 600_000_100)
    //     --csi success vs --tbi clean error with zero .tmp residue.
    let both_idx_vcf = dir.path().join("both_idx.vcf.gz");
    let both_idx_run = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&both_idx_vcf)
        .args(["--tbi", "--csi"])
        .output()
        .expect("spawn --tbi --csi");
    assert!(!both_idx_run.status.success());
    assert!(
        !both_idx_vcf.exists(),
        "both_idx.vcf.gz must not be created"
    );
    assert!(String::from_utf8_lossy(&both_idx_run.stderr)
        .contains("--tbi/--index and --csi cannot be used together"));

    #[cfg(unix)]
    {
        let alias_csi_vcf = dir.path().join("alias_csi_collision.vcf.gz");
        let alias_csi = dir.path().join("alias_csi_collision.vcf.gz.csi");
        std::fs::write(&alias_csi_vcf, b"sentinel_csi").expect("write csi sentinel");
        std::os::unix::fs::symlink(&alias_csi_vcf, &alias_csi).expect("symlink .csi to .vcf.gz");
        let alias_csi_run = base_cmd()
            .args(common_args)
            .arg("--output")
            .arg(&alias_csi_vcf)
            .arg("--csi")
            .output()
            .expect("spawn --csi with colliding .csi symlink");
        assert!(!alias_csi_run.status.success());
        assert_eq!(
            std::fs::read(&alias_csi_vcf).expect("read csi sentinel"),
            b"sentinel_csi"
        );
        assert!(String::from_utf8_lossy(&alias_csi_run.stderr)
            .contains("must not point to the same file"));

        let dev_null_csi = base_cmd()
            .args(common_args)
            .args(["--output", "/dev/null", "--csi"])
            .output()
            .expect("spawn --output /dev/null --csi");
        assert!(!dev_null_csi.status.success());
        assert!(String::from_utf8_lossy(&dev_null_csi.stderr)
            .contains("--tbi/--csi requires a regular seekable --output file"));
    }

    let empty_csi_vcf = dir.path().join("empty_csi.vcf.gz");
    let empty_csi_path = dir.path().join("empty_csi.vcf.gz.csi");
    let empty_csi_run = base_cmd()
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
        .arg(&empty_csi_vcf)
        .arg("--csi")
        .output()
        .expect("spawn header-only --csi");
    assert!(
        empty_csi_run.status.success(),
        "empty --csi failed: {}",
        String::from_utf8_lossy(&empty_csi_run.stderr)
    );
    let empty_csi_index = noodles::csi::io::Reader::new(
        std::fs::File::open(&empty_csi_path).expect("open empty .csi"),
    )
    .read_index()
    .expect("read header-only .csi");
    let empty_csi_hdr = empty_csi_index.header().expect("empty .csi VCF header");
    assert_eq!(
        empty_csi_hdr.format(),
        noodles::csi::binning_index::index::header::Format::Vcf
    );
    assert!(empty_csi_hdr.reference_sequence_names().is_empty());

    // Synthetic 700,000,000 bp contig (`chr_large`) with a variant at pos = 600,000,100:
    // sparse single-line FASTA (~20 KB real disk blocks) + .bam.csi input index.
    let large_fa = dir.path().join("large_700m.fa");
    let large_fai = dir.path().join("large_700m.fa.fai");
    let win_start_1b: usize = 599_995_001;
    let win_seq: Vec<u8> = (0..10_000usize)
        .map(|i| b"ACGT"[(i.wrapping_mul(7).wrapping_add(i / 4)) % 4])
        .collect();
    {
        use std::io::{Seek as _, SeekFrom, Write as _};
        let mut f = std::fs::File::create(&large_fa).expect("create sparse large_700m.fa");
        f.write_all(b">chr_large\n").expect("write fasta header");
        f.seek(SeekFrom::Start(11 + (win_start_1b as u64 - 1)))
            .expect("seek to 600M window");
        f.write_all(&win_seq).expect("write 600M window bases");
        f.seek(SeekFrom::Start(11 + 700_000_000))
            .expect("seek to 700M contig end");
        f.write_all(b"\n").expect("write trailing newline");
    }
    std::fs::write(
        &large_fai,
        "chr_large\t700000000\t11\t700000000\t700000001\n",
    )
    .expect("write large_700m.fa.fai");

    let large_bam = dir.path().join("large_700m.bam");
    {
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::{
            bin::Chunk, index::BinnedIndex,
        };
        use noodles::csi::binning_index::Indexer;
        use noodles::sam::{
            alignment::{
                io::Write as _,
                record::cigar::{op::Kind, Op},
                record_buf::{Cigar, RecordBuf},
            },
            header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap},
            Header,
        };

        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(std::num::NonZeroUsize::new(700_000_000).unwrap())
            .build()
            .expect("700M contig header");
        let hd = Header::builder()
            .add_reference_sequence("chr_large", rs)
            .build();
        let bgzf = noodles::bgzf::io::Writer::new(
            std::fs::File::create(&large_bam).expect("create large_700m.bam"),
        );
        let mut w = noodles::bam::io::Writer::from(bgzf);
        w.write_alignment_header(&hd).expect("write bam header");
        let mut bam_csi = Indexer::<BinnedIndex>::new(14, 6).expect("valid BAM CSI params");

        let read_start_1b: usize = 600_000_080;
        let win_off = read_start_1b - win_start_1b;
        let ref_50 = win_seq[win_off..win_off + 50].to_vec();
        let alt_b = match ref_50[20] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
        let cigar_50m: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
        for idx in 0..6usize {
            let mut seq = ref_50.clone();
            if idx < 4 {
                seq[20] = alt_b; // pos = 600_000_080 + 20 = 600_000_100
            }
            let s_pos = Position::new(read_start_1b).unwrap();
            let e_pos = Position::new(read_start_1b + 49).unwrap();
            let rec = RecordBuf::builder()
                .set_name(format!("r600m_{idx}"))
                .set_flags(noodles::sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(s_pos)
                .set_cigar(cigar_50m.clone())
                .set_sequence(seq.into())
                .set_quality_scores(vec![40u8; 50].into())
                .build();
            let start_vp = w.get_ref().virtual_position();
            w.write_alignment_record(&hd, &rec)
                .expect("write 600M bam record");
            let end_vp = w.get_ref().virtual_position();
            bam_csi
                .add_record(Some((0, s_pos, e_pos, true)), Chunk::new(start_vp, end_vp))
                .expect("index 600M bam record");
        }
        w.try_finish().expect("finish large_700m.bam");
        let bam_csi_path = dir.path().join("large_700m.bam.csi");
        let mut csi_w = noodles::csi::io::Writer::new(
            std::fs::File::create(&bam_csi_path).expect("create large_700m.bam.csi"),
        );
        csi_w
            .write_index(&bam_csi.build(1))
            .expect("write large_700m.bam.csi");
        csi_w.into_inner().finish().expect("finish bam.csi");
    }

    let large_csi_vcf = dir.path().join("large_out.vcf.gz");
    let large_csi_idx = dir.path().join("large_out.vcf.gz.csi");
    let large_csi_run = base_cmd()
        .arg("-G")
        .arg(&large_fa)
        .arg("-b")
        .arg(&large_bam)
        .args(["-R", "chr_large:600000070-600000140", "--threads=1"])
        .arg("--output")
        .arg(&large_csi_vcf)
        .arg("--csi")
        .output()
        .expect("spawn --csi on 700M contig");
    assert!(
        large_csi_run.status.success(),
        "700M --csi run failed: {}",
        String::from_utf8_lossy(&large_csi_run.stderr)
    );
    let parsed_csi = noodles::csi::io::Reader::new(
        std::fs::File::open(&large_csi_idx).expect("open large_out.vcf.gz.csi"),
    )
    .read_index()
    .expect("read large_out.vcf.gz.csi");
    let parsed_csi_hdr = parsed_csi.header().expect("CSI VCF header present");
    assert_eq!(
        parsed_csi_hdr.format(),
        noodles::csi::binning_index::index::header::Format::Vcf
    );
    let ref_id = parsed_csi_hdr
        .reference_sequence_names()
        .get_index_of(b"chr_large".as_slice())
        .expect("chr_large registered in CSI header");
    let pos_600m = noodles::core::Position::new(600_000_100).unwrap();
    let chunks = parsed_csi
        .query(
            ref_id,
            noodles::core::region::Interval::from(pos_600m..=pos_600m),
        )
        .expect("query CSI at 600,000,100");
    assert!(
        !chunks.is_empty(),
        "CSI query at 600,000,100 returned no chunks"
    );
    {
        use std::io::BufRead as _;
        let mut bgzf_r = noodles::bgzf::io::Reader::new(
            std::fs::File::open(&large_csi_vcf).expect("open large_out.vcf.gz"),
        );
        bgzf_r
            .seek(chunks[0].start())
            .expect("seek BGZF to CSI chunk start");
        let mut queried_line = String::new();
        bgzf_r
            .read_line(&mut queried_line)
            .expect("read queried VCF line");
        assert!(
            queried_line.starts_with("chr_large\t600000100\t"),
            "CSI roundtrip seek must land on chr_large:600000100 record, got: {queried_line}"
        );
    }

    let tbi_600m_dir = dir.path().join("tbi_600m_fail");
    std::fs::create_dir(&tbi_600m_dir).expect("create tbi_600m_fail");
    let large_tbi_vcf = tbi_600m_dir.join("out.vcf.gz");
    let large_tbi_run = base_cmd()
        .arg("-G")
        .arg(&large_fa)
        .arg("-b")
        .arg(&large_bam)
        .args(["-R", "chr_large:600000070-600000140", "--threads=1"])
        .arg("--output")
        .arg(&large_tbi_vcf)
        .arg("--tbi")
        .output()
        .expect("spawn --tbi on 700M contig");
    assert!(!large_tbi_run.status.success());
    assert_ne!(large_tbi_run.status.code(), Some(101));
    let large_tbi_stderr = String::from_utf8_lossy(&large_tbi_run.stderr);
    assert!(
        large_tbi_stderr.contains("536,870,912") && large_tbi_stderr.contains("--csi"),
        "expected actionable >536,870,912 --csi hint, got: {large_tbi_stderr}"
    );
    assert!(
        !large_tbi_vcf.exists(),
        "out.vcf.gz must not exist on >512 Mbp --tbi failure"
    );
    assert!(
        !tbi_600m_dir.join("out.vcf.gz.tbi").exists(),
        ".tbi file must not exist on >512 Mbp failure"
    );
    let tmp_residue: Vec<_> = std::fs::read_dir(&tbi_600m_dir)
        .expect("read tbi_600m_fail")
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains(".tmp"))
        .collect();
    assert!(
        tmp_residue.is_empty(),
        "zero .tbi.tmp residue allowed after >512 Mbp --tbi failure, found: {tmp_residue:?}"
    );

    let pre_vcf_bytes = b"PREEXISTING_VCF_GZ_SENTINEL";
    let pre_tbi_bytes = b"PREEXISTING_TBI_SENTINEL";
    let pre_csi_bytes = b"PREEXISTING_CSI_SENTINEL";
    let large_tbi_idx = tbi_600m_dir.join("out.vcf.gz.tbi");
    let large_csi_sidecar = tbi_600m_dir.join("out.vcf.gz.csi");
    std::fs::write(&large_tbi_vcf, pre_vcf_bytes).expect("seed preexisting out.vcf.gz");
    std::fs::write(&large_tbi_idx, pre_tbi_bytes).expect("seed preexisting out.vcf.gz.tbi");
    std::fs::write(&large_csi_sidecar, pre_csi_bytes).expect("seed preexisting out.vcf.gz.csi");
    let large_tbi_preserve = base_cmd()
        .arg("-G")
        .arg(&large_fa)
        .arg("-b")
        .arg(&large_bam)
        .args(["-R", "chr_large:600000070-600000140", "--threads=1"])
        .arg("--output")
        .arg(&large_tbi_vcf)
        .arg("--tbi")
        .output()
        .expect("spawn --tbi over preexisting files on 700M contig");
    assert!(!large_tbi_preserve.status.success());
    assert_eq!(
        std::fs::read(&large_tbi_vcf).expect("read preexisting out.vcf.gz"),
        pre_vcf_bytes,
        "failing --tbi must leave preexisting out.vcf.gz bytes untouched"
    );
    assert_eq!(
        std::fs::read(&large_tbi_idx).expect("read preexisting out.vcf.gz.tbi"),
        pre_tbi_bytes,
        "failing --tbi must leave preexisting out.vcf.gz.tbi bytes untouched"
    );
    assert_eq!(
        std::fs::read(&large_csi_sidecar).expect("read preexisting out.vcf.gz.csi"),
        pre_csi_bytes,
        "failing pre-commit --tbi must preserve preexisting alternate .csi sidecar"
    );
    let tmp_residue_preserve: Vec<_> = std::fs::read_dir(&tbi_600m_dir)
        .expect("read tbi_600m_fail after preserve")
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains(".tmp"))
        .collect();
    assert!(
        tmp_residue_preserve.is_empty(),
        "zero .tmp residue allowed after preexisting --tbi failure, found: {tmp_residue_preserve:?}"
    );

    // 3c. Directory .tbi / .csi destinations fail before staging/commit and leave preexisting VCF untouched;
    //     rewriting with --tbi after --csi (and --csi after --tbi) invalidates the stale alternate sidecar
    //     while rejecting input / --stats-json alternate-sidecar collisions.
    let dir_target_vcf = dir.path().join("dir_target.vcf.gz");
    let dir_target_sentinel = b"PREEXISTING_VCF_BEFORE_DIR_IDX";
    std::fs::write(&dir_target_vcf, dir_target_sentinel).expect("seed dir_target.vcf.gz");
    for (idx_flag, idx_ext) in [("--tbi", "tbi"), ("--csi", "csi")] {
        let idx_dir = dir.path().join(format!("dir_target.vcf.gz.{idx_ext}"));
        std::fs::create_dir(&idx_dir).expect("create index destination directory");
        let dir_idx_run = base_cmd()
            .args(common_args)
            .arg("--output")
            .arg(&dir_target_vcf)
            .arg(idx_flag)
            .output()
            .expect("spawn with directory index destination");
        assert!(
            !dir_idx_run.status.success(),
            "{idx_flag} with directory destination must fail"
        );
        assert_ne!(dir_idx_run.status.code(), Some(101));
        assert_eq!(
            std::fs::read(&dir_target_vcf).expect("read dir_target.vcf.gz"),
            dir_target_sentinel,
            "directory {idx_flag} destination must leave preexisting VCF unchanged"
        );
        assert!(idx_dir.is_dir(), "{idx_flag} directory must remain intact");
        std::fs::remove_dir(&idx_dir).expect("cleanup index destination directory");
    }

    let swap_vcf = dir.path().join("sidecar_swap.vcf.gz");
    let swap_tbi = dir.path().join("sidecar_swap.vcf.gz.tbi");
    let swap_csi = dir.path().join("sidecar_swap.vcf.gz.csi");
    let swap_csi_first = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&swap_vcf)
        .arg("--csi")
        .output()
        .expect("spawn initial --csi for sidecar swap");
    assert!(
        swap_csi_first.status.success(),
        "initial --csi failed: {}",
        String::from_utf8_lossy(&swap_csi_first.stderr)
    );
    assert!(swap_csi.is_file() && !swap_tbi.exists());

    let swap_to_tbi = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&swap_vcf)
        .arg("--tbi")
        .output()
        .expect("spawn replacement --tbi over prior --csi");
    assert!(
        swap_to_tbi.status.success(),
        "replacement --tbi failed: {}",
        String::from_utf8_lossy(&swap_to_tbi.stderr)
    );
    assert!(swap_tbi.is_file(), "replacement .tbi must exist");
    assert!(
        !swap_csi.exists(),
        "stale .csi sidecar must be invalidated after --tbi replacement"
    );
    let _ = noodles::tabix::fs::read(&swap_tbi).expect("replacement .tbi must parse");

    let swap_to_csi = base_cmd()
        .args(common_args)
        .arg("--output")
        .arg(&swap_vcf)
        .arg("--csi")
        .output()
        .expect("spawn replacement --csi over prior --tbi");
    assert!(
        swap_to_csi.status.success(),
        "replacement --csi failed: {}",
        String::from_utf8_lossy(&swap_to_csi.stderr)
    );
    assert!(swap_csi.is_file(), "replacement .csi must exist");
    assert!(
        !swap_tbi.exists(),
        "stale .tbi sidecar must be invalidated after --csi replacement"
    );
    let _ = noodles::csi::io::Reader::new(std::fs::File::open(&swap_csi).expect("open swap .csi"))
        .read_index()
        .expect("replacement .csi must parse");

    #[cfg(unix)]
    {
        let alt_stats_vcf = dir.path().join("alt_stats_col.vcf.gz");
        let alt_stats_csi = dir.path().join("alt_stats_col.vcf.gz.csi");
        let stats_target = dir.path().join("protected_stats.json");
        std::fs::write(&alt_stats_vcf, b"ALT_STATS_VCF").expect("seed alt_stats_col.vcf.gz");
        std::fs::write(&stats_target, b"{\"keep\":true}").expect("seed protected_stats.json");
        std::os::unix::fs::symlink(&stats_target, &alt_stats_csi)
            .expect("symlink alternate .csi to --stats-json");
        let alt_stats_run = base_cmd()
            .args(common_args)
            .arg("--output")
            .arg(&alt_stats_vcf)
            .arg("--tbi")
            .arg("--stats-json")
            .arg(&stats_target)
            .output()
            .expect("spawn --tbi with alternate .csi colliding with --stats-json");
        assert!(!alt_stats_run.status.success());
        assert!(String::from_utf8_lossy(&alt_stats_run.stderr)
            .contains("must not point to the same file"));
        assert_eq!(
            std::fs::read(&alt_stats_vcf).expect("read alt_stats_col.vcf.gz"),
            b"ALT_STATS_VCF"
        );
        assert_eq!(
            std::fs::read(&stats_target).expect("read protected_stats.json"),
            b"{\"keep\":true}"
        );

        let alt_in_vcf = dir.path().join("alt_input_col.vcf.gz");
        let alt_in_csi = dir.path().join("alt_input_col.vcf.gz.csi");
        let micro_bam_abs = workspace_root().join("fixtures/micro/micro.bam");
        let micro_bam_len = std::fs::metadata(&micro_bam_abs)
            .expect("stat micro.bam")
            .len();
        std::fs::write(&alt_in_vcf, b"ALT_INPUT_VCF").expect("seed alt_input_col.vcf.gz");
        std::os::unix::fs::symlink(&micro_bam_abs, &alt_in_csi)
            .expect("symlink alternate .csi to input BAM");
        let alt_in_run = base_cmd()
            .args(common_args)
            .arg("--output")
            .arg(&alt_in_vcf)
            .arg("--tbi")
            .output()
            .expect("spawn --tbi with alternate .csi colliding with input BAM");
        assert!(!alt_in_run.status.success());
        assert!(
            String::from_utf8_lossy(&alt_in_run.stderr).contains("must not point to the same file")
        );
        assert_eq!(
            std::fs::read(&alt_in_vcf).expect("read alt_input_col.vcf.gz"),
            b"ALT_INPUT_VCF"
        );
        assert_eq!(
            std::fs::metadata(&micro_bam_abs)
                .expect("stat micro.bam after")
                .len(),
            micro_bam_len
        );
    }

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

#[test]
fn broken_pipe_stdout_and_stderr_exit_cleanly_and_preserve_stats_json() {
    use std::io::{BufRead as _, BufReader, Read as _};

    let dir = tempfile::TempDir::new().expect("tempdir");
    // Repeat micro.bed intervals 20x under `-p` so uncompressed output (~1.1 MB)
    // far exceeds the 64 KB BufWriter + 64 KB OS pipe buffer, forcing mid-stream
    // EPIPE (`BrokenPipe`) while `-th 4` workers are active.
    let repeated_bed = dir.path().join("repeated_micro.bed");
    let mut bed_text = String::new();
    for i in 0..20 {
        bed_text.push_str(&format!(
            "chr1\t40\t260\tmicro_chr1_{i}\nchr2\t20\t250\tmicro_chr2_{i}\n"
        ));
    }
    std::fs::write(&repeated_bed, &bed_text).expect("write repeated_micro.bed");

    // 1. Captured stderr + dropped stdout after 1 line -> exit 0, zero panic text, valid stats.json
    let stats_captured = dir.path().join("stats_captured_stderr.json");
    let mut cmd_a = base_cmd();
    cmd_a
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
            "-p",
            "-th",
            "4",
            "-y",
        ])
        .arg("--stats-json")
        .arg(&stats_captured)
        .arg(&repeated_bed)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child_a = cmd_a.spawn().expect("spawn vardictrust for stdout EPIPE");
    drop(cmd_a);

    let mut first_line = String::new();
    {
        let stdout_a = child_a.stdout.take().expect("child stdout");
        let mut reader = BufReader::new(stdout_a);
        let n = reader
            .read_line(&mut first_line)
            .expect("read first VCF header line");
        assert!(n > 0, "expected at least one line on stdout");
        assert!(
            first_line.starts_with("##fileformat=VCF"),
            "unexpected first stdout line: {first_line}"
        );
    } // drop child stdout immediately after 1 line -> downstream EPIPE

    let mut stderr_bytes = Vec::new();
    if let Some(mut err_pipe) = child_a.stderr.take() {
        err_pipe
            .read_to_end(&mut stderr_bytes)
            .expect("drain child stderr");
    }
    let status_a = child_a.wait().expect("wait child_a");
    let stderr_str = String::from_utf8_lossy(&stderr_bytes);
    assert!(
        status_a.success(),
        "BrokenPipe on stdout must exit 0, got {status_a:?}, stderr:\n{stderr_str}"
    );
    #[cfg(unix)]
    assert_eq!(status_a.signal(), None, "child died by signal on EPIPE");
    assert_ne!(status_a.code(), Some(101), "child panicked on EPIPE");
    assert!(
        !stderr_str.contains("panicked"),
        "stderr must not contain panic message: {stderr_str}"
    );

    // 2. Closed stderr under `-y` + dropped stdout after 1 line + `--stats-json` -> exit 0 + valid JSON
    let stats_closed_err = dir.path().join("stats_closed_stderr.json");
    let mut cmd_b = base_cmd();
    cmd_b
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
            "-p",
            "-th",
            "4",
            "-y",
        ])
        .arg("--stats-json")
        .arg(&stats_closed_err)
        .arg(&repeated_bed)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child_b = cmd_b
        .spawn()
        .expect("spawn vardictrust with piped stderr and stdout");
    drop(cmd_b);
    // Close parent's stderr read end immediately so `-y` verbose writes hit EPIPE (`BrokenPipe`).
    drop(child_b.stderr.take());

    {
        let stdout_b = child_b.stdout.take().expect("child_b stdout");
        let mut reader = BufReader::new(stdout_b);
        let mut line = String::new();
        let n = reader.read_line(&mut line).expect("read first line");
        assert!(n > 0 && line.starts_with("##fileformat=VCF"));
    }
    let status_b = child_b.wait().expect("wait child_b");
    assert!(
        status_b.success(),
        "BrokenPipe on both `-y` stderr and stdout must exit 0, got {status_b:?}"
    );
    #[cfg(unix)]
    assert_eq!(status_b.signal(), None, "child_b died by signal");
    assert_ne!(status_b.code(), Some(101), "child_b panicked");

    for stats_path in [&stats_captured, &stats_closed_err] {
        let text = std::fs::read_to_string(stats_path).unwrap_or_else(|e| {
            panic!(
                "stats-json must be preserved on BrokenPipe ({}): {e}",
                stats_path.display()
            )
        });
        let trimmed = text.trim();
        assert!(
            trimmed.starts_with('{') && trimmed.ends_with('}'),
            "stats-json must be a complete JSON object, got:\n{text}"
        );
        for key in [
            "\"mode\": \"Simple\"",
            "\"threads\": 4",
            "\"regions_count\": 40",
            "\"tiles_count\":",
            "\"elapsed_ms\":",
            "\"peak_rss_kb\":",
            "\"exceptions_caught\": 0",
        ] {
            assert!(
                text.contains(key),
                "stats-json ({}) missing '{key}', got:\n{text}",
                stats_path.display()
            );
        }
    }
}

#[test]
fn output_and_stats_json_reject_input_collisions_and_leave_no_empty_residue() {
    let dir = tempfile::TempDir::new().expect("tempdir");
    let root = workspace_root();

    let ref_fa = dir.path().join("micro.ref.fa");
    let ref_fai = dir.path().join("micro.ref.fa.fai");
    let tumor_bam = dir.path().join("tumor.bam");
    let tumor_bai = dir.path().join("tumor.bam.bai");
    let normal_bam = dir.path().join("normal.bam");
    let normal_bai = dir.path().join("normal.bam.bai");
    let bed = dir.path().join("micro.bed");

    std::fs::copy(root.join("fixtures/micro/micro.ref.fa"), &ref_fa).expect("copy ref.fa");
    std::fs::copy(root.join("fixtures/micro/micro.ref.fa.fai"), &ref_fai).expect("copy ref.fa.fai");
    std::fs::copy(root.join("fixtures/micro/micro.bam"), &tumor_bam).expect("copy tumor.bam");
    std::fs::copy(root.join("fixtures/micro/micro.bam.bai"), &tumor_bai)
        .expect("copy tumor.bam.bai");
    std::fs::copy(root.join("fixtures/micro/micro.bam"), &normal_bam).expect("copy normal.bam");
    std::fs::copy(root.join("fixtures/micro/micro.bam.bai"), &normal_bai)
        .expect("copy normal.bam.bai");
    std::fs::copy(root.join("fixtures/micro/micro.bed"), &bed).expect("copy micro.bed");

    let protected = [
        &ref_fa,
        &ref_fai,
        &tumor_bam,
        &tumor_bai,
        &normal_bam,
        &normal_bai,
        &bed,
    ];
    let snapshots: Vec<Vec<u8>> = protected
        .iter()
        .map(|p| std::fs::read(p).expect("snapshot input"))
        .collect();

    let assert_parse_collision = |cmd: &mut Command, label: &str| {
        let out = cmd
            .output()
            .unwrap_or_else(|e| panic!("spawn {label}: {e}"));
        assert!(!out.status.success(), "{label} must fail");
        assert_ne!(out.status.code(), Some(101), "{label} panicked");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("parse:") && stderr.contains("must not point to the same file"),
            "{label} must fail with Parse collision error, got: {stderr}"
        );
    };

    // 1. --output <bam>
    assert_parse_collision(
        base_cmd()
            .arg("-G")
            .arg(&ref_fa)
            .arg("-b")
            .arg(&tumor_bam)
            .args(["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z"])
            .arg(&bed)
            .arg("--output")
            .arg(&tumor_bam),
        "--output <bam>",
    );

    // 2. --output <bam.bai>
    assert_parse_collision(
        base_cmd()
            .arg("-G")
            .arg(&ref_fa)
            .arg("-b")
            .arg(&tumor_bam)
            .args(["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z"])
            .arg(&bed)
            .arg("--output")
            .arg(&tumor_bai),
        "--output <bam.bai>",
    );

    // 3. -i -b "<tumor>|<normal>" --output <normal>
    let paired_arg = format!("{}|{}", tumor_bam.display(), normal_bam.display());
    assert_parse_collision(
        base_cmd()
            .arg("-G")
            .arg(&ref_fa)
            .arg("-i")
            .arg("-b")
            .arg(&paired_arg)
            .args(["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z"])
            .arg(&bed)
            .arg("--output")
            .arg(&normal_bam),
        "-i -b <tumor>|<normal> --output <normal>",
    );

    // 4. --output <symlink_to_ref>
    #[cfg(unix)]
    {
        let symlink_to_ref = dir.path().join("symlink_to_ref.vcf");
        std::os::unix::fs::symlink(&ref_fa, &symlink_to_ref).expect("symlink ref.fa");
        assert_parse_collision(
            base_cmd()
                .arg("-G")
                .arg(&ref_fa)
                .arg("-b")
                .arg(&tumor_bam)
                .args(["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z"])
                .arg(&bed)
                .arg("--output")
                .arg(&symlink_to_ref),
            "--output <symlink_to_ref>",
        );
    }

    // 5. --stats-json <bed>
    assert_parse_collision(
        base_cmd()
            .arg("-G")
            .arg(&ref_fa)
            .arg("-b")
            .arg(&tumor_bam)
            .args(["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z"])
            .arg(&bed)
            .arg("--stats-json")
            .arg(&bed),
        "--stats-json <bed>",
    );

    // Verify zero truncation on any input/companion file
    for (p, expected) in protected.iter().zip(snapshots.iter()) {
        let actual = std::fs::read(p).expect("read protected input after collision checks");
        assert!(
            !actual.is_empty() && &actual == expected,
            "input file {} was truncated or modified!",
            p.display()
        );
    }

    // 6. --output /dev/null --stats-json /dev/null succeeds on Unix
    #[cfg(unix)]
    {
        let dev_null = base_cmd()
            .arg("-G")
            .arg(&ref_fa)
            .arg("-b")
            .arg(&tumor_bam)
            .args([
                "-R",
                "chr1:40-260",
                "--threads=1",
                "--output",
                "/dev/null",
                "--stats-json",
                "/dev/null",
            ])
            .output()
            .expect("spawn --output /dev/null --stats-json /dev/null");
        assert!(
            dev_null.status.success(),
            "--output /dev/null --stats-json /dev/null failed: {}",
            String::from_utf8_lossy(&dev_null.stderr)
        );
    }

    // 7. Malformed `-R chr1:bad --output new.vcf.gz` leaves NO 28-byte `new.vcf.gz` residue
    let new_vcf_gz = dir.path().join("new.vcf.gz");
    let new_stats = dir.path().join("new_stats.json");
    let bad_region = base_cmd()
        .arg("-G")
        .arg(&ref_fa)
        .arg("-b")
        .arg(&tumor_bam)
        .args(["-R", "chr1:bad", "--threads=1"])
        .arg("--output")
        .arg(&new_vcf_gz)
        .arg("--stats-json")
        .arg(&new_stats)
        .output()
        .expect("spawn malformed -R chr1:bad --output new.vcf.gz");
    assert!(!bad_region.status.success(), "-R chr1:bad must fail");
    let bad_stderr = String::from_utf8_lossy(&bad_region.stderr);
    assert!(
        bad_stderr.contains("parse:"),
        "expected Parse error on -R chr1:bad, got: {bad_stderr}"
    );
    assert!(
        !new_vcf_gz.exists(),
        "malformed -R chr1:bad must leave NO 28-byte new.vcf.gz residue"
    );
    assert!(
        !new_stats.exists(),
        "malformed -R chr1:bad must not commit stats.json"
    );
}

#[test]
fn bgzf_tbi_somatic_amplicon_and_sv_end_to_end() {
    use noodles::csi::BinningIndex as _;
    use std::io::Read as _;

    ensure_fastas_converted();
    let droot = data_root();
    let ref_hg19 = droot.join("VarDictJava/testdata/fastas/hg19.fa");
    let bam_dir =
        droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
    let dir = tempfile::TempDir::new().expect("tempdir");

    let query_bgzf_tbi = |vcf_path: &Path, tbi_path: &Path, region_str: &str| -> Vec<String> {
        // Decompress via `bgzf::io::IndexedReader` to verify valid BGZF framing
        let mut bgzf_indexed = noodles::bgzf::io::IndexedReader::new(
            std::fs::File::open(vcf_path).expect("open BGZF"),
            noodles::bgzf::gzi::Index::default(),
        );
        let mut raw_vcf = String::new();
        bgzf_indexed
            .read_to_string(&mut raw_vcf)
            .expect("decompress BGZF stream");
        assert!(
            raw_vcf.contains("#CHROM\tPOS\t"),
            "BGZF output {} missing VCF header",
            vcf_path.display()
        );

        // Read `.tbi` via `noodles::tabix` and verify both binning index chunks and VCF query records
        let index = noodles::tabix::fs::read(tbi_path).expect("read .tbi");
        let hdr = index.header().expect("tabix header");
        assert_eq!(
            hdr.format(),
            noodles::csi::binning_index::index::header::Format::Vcf
        );
        let region: noodles::core::Region = region_str.parse().expect("parse query region");
        let ref_seq_id = hdr
            .reference_sequence_names()
            .get_index_of(region.name())
            .expect("contig present in .tbi header");
        let chunks = index
            .query(ref_seq_id, region.interval())
            .expect("query .tbi chunks");
        assert!(
            !chunks.is_empty(),
            "expected non-empty .tbi chunks for {region_str} in {}",
            tbi_path.display()
        );

        // Seek BGZF virtual-position chunks and collect overlapping VCF lines (honoring INFO/END=)
        let mut vcf_reader = noodles::vcf::io::IndexedReader::new(
            std::fs::File::open(vcf_path).expect("open indexed VCF"),
            index,
        );
        let vcf_header = vcf_reader.read_header().expect("read VCF header");
        let query = vcf_reader
            .query(&vcf_header, &region)
            .expect("query indexed VCF");
        query
            .records()
            .map(|r| {
                let rec = r.expect("valid VCF record");
                let pos = rec
                    .variant_start()
                    .expect("POS present")
                    .expect("valid POS")
                    .get();
                format!(
                    "{}\t{pos}\t.\t{}\t{}\t.\t.\t{}",
                    rec.reference_sequence_name(),
                    rec.reference_bases(),
                    rec.alternate_bases().as_ref(),
                    rec.info().as_ref(),
                )
            })
            .collect()
    };

    // 1. Somatic: `--output som.vcf.gz --tbi` (`-b "tumor|normal"`)
    let som_vcf_gz = dir.path().join("som.vcf.gz");
    let som_tbi = dir.path().join("som.vcf.gz.tbi");
    let som_bams = format!(
        "{}|{}",
        bam_dir.join("ArgumentNA_Tumor.bam").display(),
        bam_dir.join("ArgumentNA_Normal.bam").display()
    );
    let som_run = base_cmd()
        .arg("-G")
        .arg(&ref_hg19)
        .arg("-b")
        .arg(&som_bams)
        .args([
            "-R",
            "chr7:140453010-140453158:testbed",
            "-z",
            "-f",
            "0.001",
            "--threads=2",
        ])
        .arg("--output")
        .arg(&som_vcf_gz)
        .arg("--tbi")
        .output()
        .expect("spawn somatic --tbi");
    assert!(
        som_run.status.success(),
        "somatic --tbi failed: {}",
        String::from_utf8_lossy(&som_run.stderr)
    );
    let som_hits = query_bgzf_tbi(&som_vcf_gz, &som_tbi, "chr7:140453030-140453040");
    assert!(
        som_hits.iter().any(|l| l.contains("STATUS=")),
        "expected somatic STATUS= record in chr7:140453030-140453040, got: {som_hits:?}"
    );

    // 2. Amplicon: `--output amp.vcf.bgz --index` (`-a`)
    let amp_bed = dir.path().join("amp.bed");
    std::fs::write(&amp_bed, "chr1\t81\t140\tmicro_amp\t0\t+\t81\t140\n").expect("write amp.bed");
    let amp_vcf_bgz = dir.path().join("amp.vcf.bgz");
    let amp_tbi = dir.path().join("amp.vcf.bgz.tbi");
    let amp_run = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-a",
            "30:0.5",
            "-f",
            "0.01",
            "-r",
            "1",
            "-z",
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "--threads=2",
        ])
        .arg(&amp_bed)
        .arg("--output")
        .arg(&amp_vcf_bgz)
        .arg("--index")
        .output()
        .expect("spawn amplicon .vcf.bgz --index");
    assert!(
        amp_run.status.success(),
        "amplicon --index failed: {}",
        String::from_utf8_lossy(&amp_run.stderr)
    );
    let amp_hits = query_bgzf_tbi(&amp_vcf_bgz, &amp_tbi, "chr1:95-105");
    assert!(
        amp_hits
            .iter()
            .any(|l| l.starts_with("chr1\t101\t") && (l.contains("GDAMP=") || l.contains("TLAMP="))),
        "expected amplicon GDAMP/TLAMP record at chr1:101, got: {amp_hits:?}"
    );

    // 3. Symbolic SV: `--output sv.vcf.gz --tbi`, mid-span query between POS=89700386 and END=89712332
    let sv_vcf_gz = dir.path().join("sv.vcf.gz");
    let sv_tbi = dir.path().join("sv.vcf.gz.tbi");
    let sv_bam = bam_dir.join("Colo829.chr10.DELDiscReverse.bam");
    let sv_run = base_cmd()
        .arg("-G")
        .arg(&ref_hg19)
        .arg("-b")
        .arg(&sv_bam)
        .args([
            "-R",
            "chr10:89712340-89713976:testbed",
            "-z",
            "-f",
            "0.001",
            "--threads=2",
        ])
        .arg("--output")
        .arg(&sv_vcf_gz)
        .arg("--tbi")
        .output()
        .expect("spawn SV --tbi");
    assert!(
        sv_run.status.success(),
        "SV --tbi failed: {}",
        String::from_utf8_lossy(&sv_run.stderr)
    );
    // Mid-span window chr10:89706000-89706100 sits >5.6 kbp past POS=89700386 and >6.2 kbp
    // before INFO/END=89712332, containing zero point variants; only `END=` indexing returns `<DEL>`.
    let sv_hits = query_bgzf_tbi(&sv_vcf_gz, &sv_tbi, "chr10:89706000-89706100");
    assert!(
        sv_hits
            .iter()
            .any(|l| l.starts_with("chr10\t89700386\t") && l.contains("\t<DEL>\t") && l.contains("END=89712332")),
        "mid-span .tbi query must recover symbolic <DEL> spanning POS=89700386..END=89712332, got: {sv_hits:?}"
    );
}

#[test]
fn hostile_missing_qual_huge_deletion_and_i32_min_tlen_do_not_panic_or_oom() {
    use noodles::core::Position;
    use noodles::sam::alignment::{
        record::{
            cigar::{op::Kind, Op},
            Flags,
        },
        record_buf::{Cigar, QualityScores, RecordBuf},
    };

    let dir = tempfile::TempDir::new().expect("tempdir");
    let ref_fa = dir.path().join("hostile_ref.fa");
    let ref_fai = dir.path().join("hostile_ref.fa.fai");
    let ref_len = 5000usize;
    std::fs::write(&ref_fa, format!(">chr1\n{}\n", "A".repeat(ref_len))).expect("write ref.fa");
    std::fs::write(&ref_fai, format!("chr1\t{ref_len}\t6\t{ref_len}\t5001\n"))
        .expect("write ref.fa.fai");

    let mut snv_seq = vec![b'A'; 50];
    snv_seq[20] = b'C'; // SNV at chr1:100 when alignment_start == 80
    let match_50: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();

    // 1. `*` missing base quality scores (`QualityScores::default()`)
    let missing_qual_read = RecordBuf::builder()
        .set_name("missing_qual_star")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(75usize).unwrap())
        .set_cigar(match_50.clone())
        .set_sequence(snv_seq.clone().into())
        .set_quality_scores(QualityScores::default())
        .build();

    // 2. Valid high-quality SNV support reads at chr1:80 (variant at chr1:100)
    let valid_snv_1 = RecordBuf::builder()
        .set_name("valid_snv_1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(80usize).unwrap())
        .set_cigar(match_50.clone())
        .set_sequence(snv_seq.clone().into())
        .set_quality_scores(vec![40u8; 50].into())
        .build();
    let valid_snv_2 = RecordBuf::builder()
        .set_name("valid_snv_2")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(80usize).unwrap())
        .set_cigar(match_50)
        .set_sequence(snv_seq.into())
        .set_quality_scores(vec![40u8; 50].into())
        .build();

    // 3. Paired read with `TLEN = i32::MIN` (-2147483648) and soft-clip
    let sc_cigar: Cigar = [Op::new(Kind::SoftClip, 10), Op::new(Kind::Match, 40)]
        .into_iter()
        .collect();
    let i32_min_tlen_read = RecordBuf::builder()
        .set_name("tlen_i32_min")
        .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_REVERSE_COMPLEMENTED)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(90usize).unwrap())
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::try_from(500usize).unwrap())
        .set_template_length(i32::MIN)
        .set_cigar(sc_cigar)
        .set_sequence(vec![b'A'; 50].into())
        .set_quality_scores(vec![35u8; 50].into())
        .build();

    // 4. `1M268435455D1M` (28-bit max BAM CIGAR deletion op)
    let huge_del_cigar: Cigar = [
        Op::new(Kind::Match, 1),
        Op::new(Kind::Deletion, 268_435_455),
        Op::new(Kind::Match, 1),
    ]
    .into_iter()
    .collect();
    let huge_del_read = RecordBuf::builder()
        .set_name("huge_268m_del")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(120usize).unwrap())
        .set_cigar(huge_del_cigar)
        .set_sequence(b"AA".to_vec().into())
        .set_quality_scores(vec![40u8; 2].into())
        .build();

    let recs = vec![
        missing_qual_read,
        valid_snv_1,
        valid_snv_2,
        i32_min_tlen_read,
        huge_del_read,
    ];
    let hostile_bam = dir.path().join("hostile.bam");
    write_indexed_bam(&hostile_bam, "chr1", 300_000_000, &recs).expect("write hostile.bam");

    let out = base_cmd()
        .arg("-G")
        .arg(&ref_fa)
        .arg("-b")
        .arg(&hostile_bam)
        .args(["-R", "chr1:1-1000", "-f", "0.01", "--threads=2"])
        .output()
        .expect("spawn vardictrust on hostile.bam");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.success(),
        "hostile BAM must execute cleanly without OOM or panic, status={:?}, stderr:\n{stderr}",
        out.status
    );
    #[cfg(unix)]
    assert_eq!(out.status.signal(), None, "hostile BAM died by signal");
    assert_ne!(out.status.code(), Some(101), "hostile BAM panicked");
    assert!(!stderr.contains("panicked"), "unexpected panic in stderr");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("#CHROM\tPOS\t"));
    assert!(
        stdout.contains("chr1\t100\t"),
        "valid SNV at chr1:100 must still be called alongside hostile records, got:\n{stdout}"
    );
}

#[test]
fn commons_cli_single_dash_flags_and_zero_bed_col_behavior() {
    let base_region_args = [
        "-G",
        "fixtures/micro/micro.ref.fa",
        "-b",
        "fixtures/micro/micro.bam",
        "-R",
        "chr1:40-260",
        "--threads=1",
    ];

    // Baseline default run
    let default_out = base_cmd()
        .args(base_region_args)
        .output()
        .expect("spawn default baseline");
    assert!(default_out.status.success());

    // 1. `-downsample 0.5` -> Unsupported error (exit non-zero, no panic)
    let downsample_out = base_cmd()
        .args(base_region_args)
        .args(["-downsample", "0.5"])
        .output()
        .expect("spawn -downsample 0.5");
    assert!(!downsample_out.status.success());
    assert_ne!(downsample_out.status.code(), Some(101));
    let downsample_stderr = String::from_utf8_lossy(&downsample_out.stderr);
    assert!(
        downsample_stderr.contains("unsupported:")
            && downsample_stderr.contains("-downsample not implemented in v1.0; use VarDictJava"),
        "expected Unsupported error for -downsample 0.5, got: {downsample_stderr}"
    );

    // 2. `-adaptor ACGTACGTACGT` -> accepted Commons-CLI single-dash long flag
    let adaptor_out = base_cmd()
        .args(base_region_args)
        .args(["-adaptor", "ACGTACGTACGT"])
        .output()
        .expect("spawn -adaptor ACGTACGTACGT");
    assert!(
        adaptor_out.status.success(),
        "-adaptor failed: {}",
        String::from_utf8_lossy(&adaptor_out.stderr)
    );
    assert!(String::from_utf8_lossy(&adaptor_out.stdout).contains("chr1\t101\t"));

    // 3. `-fisher` matches `--fisher` byte-for-byte
    let fisher_single = base_cmd()
        .args(base_region_args)
        .arg("-fisher")
        .output()
        .expect("spawn -fisher");
    let fisher_double = base_cmd()
        .args(base_region_args)
        .arg("--fisher")
        .output()
        .expect("spawn --fisher");
    assert!(fisher_single.status.success());
    assert!(fisher_double.status.success());
    assert_eq!(fisher_single.stdout, fisher_double.stdout);
    let fisher_vcf = String::from_utf8_lossy(&fisher_single.stdout);
    assert!(fisher_vcf.contains("SBF=") && fisher_vcf.contains("ODDRATIO="));

    // 4. `-N -fisher` preserves `-fisher` verbatim as the sample name (not rewritten to `--fisher`)
    let sample_dash_fisher = base_cmd()
        .args(base_region_args)
        .args(["-N", "-fisher"])
        .output()
        .expect("spawn -N -fisher");
    assert!(
        sample_dash_fisher.status.success(),
        "-N -fisher failed: {}",
        String::from_utf8_lossy(&sample_dash_fisher.stderr)
    );
    let sample_vcf = String::from_utf8_lossy(&sample_dash_fisher.stdout);
    assert!(
        sample_vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t-fisher\n"),
        "-N -fisher must set VCF sample header to literal '-fisher', got:\n{sample_vcf}"
    );

    // 5. `-k` (bare optional-arg flag) accepted and matches default `realign = 1`
    let bare_k = base_cmd()
        .args(base_region_args)
        .arg("-k")
        .output()
        .expect("spawn bare -k");
    assert!(
        bare_k.status.success(),
        "bare -k failed: {}",
        String::from_utf8_lossy(&bare_k.stderr)
    );
    assert_eq!(bare_k.stdout, default_out.stdout);

    // 6. `-c 0`: Parse error without `-R`, accepted when `-R` is supplied
    let c0_no_region = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-c",
            "0",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-z",
            "fixtures/micro/micro.bed",
        ])
        .output()
        .expect("spawn -c 0 without -R");
    assert!(
        !c0_no_region.status.success(),
        "-c 0 without -R must fail with Parse error"
    );
    assert_ne!(c0_no_region.status.code(), Some(101));
    let c0_stderr = String::from_utf8_lossy(&c0_no_region.stderr);
    assert!(
        c0_stderr.contains("parse:") && c0_stderr.contains("Invalid -c value: 0"),
        "expected Parse error for -c 0 without -R, got: {c0_stderr}"
    );

    let c0_with_region = base_cmd()
        .args(base_region_args)
        .args(["-c", "0"])
        .output()
        .expect("spawn -c 0 with -R");
    assert!(
        c0_with_region.status.success(),
        "-c 0 with -R must succeed, got stderr: {}",
        String::from_utf8_lossy(&c0_with_region.stderr)
    );
    assert_eq!(c0_with_region.stdout, default_out.stdout);

    // 7. Pre-existing --output .vcf.gz is preserved byte-for-byte on failed CLI run
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let pre_vcf = tmp.path().join("existing.vcf.gz");
    std::fs::write(&pre_vcf, b"PREEXISTING_VCF_BYTES").expect("seed existing.vcf.gz");
    let fail_pre = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:not_a_coord",
            "--output",
        ])
        .arg(&pre_vcf)
        .output()
        .expect("spawn failed run with pre-existing --output");
    assert!(!fail_pre.status.success());
    assert_eq!(
        std::fs::read(&pre_vcf).expect("read pre_vcf"),
        b"PREEXISTING_VCF_BYTES",
        "failed CLI run must not truncate pre-existing --output file"
    );

    // 8. Non-UTF-8 OS CLI argument exits cleanly (code != 101, zero panic) on Unix
    #[cfg(unix)]
    {
        use std::ffi::OsStr;
        use std::os::unix::ffi::OsStrExt as _;
        let bad_utf8 = base_cmd()
            .arg(OsStr::from_bytes(b"\xff\xfe"))
            .output()
            .expect("spawn non-UTF-8 arg");
        assert!(!bad_utf8.status.success());
        assert_ne!(bad_utf8.status.code(), Some(101), "non-UTF-8 arg panicked");
        let bad_stderr = String::from_utf8_lossy(&bad_utf8.stderr);
        assert!(
            !bad_stderr.contains("panicked") && bad_stderr.contains("not valid UTF-8"),
            "expected clean UTF-8 Parse error without panic, got: {bad_stderr}"
        );
    }

    // 9. Binary E2E runs for `-n` (sample regex), clustered `-uz 0`, `-VS LENIENT` / `STRICT`,
    //    `--no-end`, and `RAYON_NUM_THREADS=2`.
    let named_bam = tmp.path().join("sample_42_run.bam");
    let named_bai = tmp.path().join("sample_42_run.bam.bai");
    std::fs::copy(
        workspace_root().join("fixtures/micro/micro.bam"),
        &named_bam,
    )
    .expect("copy micro.bam -> sample_42_run.bam");
    std::fs::copy(
        workspace_root().join("fixtures/micro/micro.bam.bai"),
        &named_bai,
    )
    .expect("copy micro.bam.bai -> sample_42_run.bam.bai");

    let n_single = base_cmd()
        .args(["-G", "fixtures/micro/micro.ref.fa", "-b"])
        .arg(&named_bam)
        .args([
            "-R",
            "chr1:40-260",
            "-n",
            r"/sample_(\d+)_run/",
            "--threads=1",
        ])
        .output()
        .expect("spawn -n single-sample VCF");
    assert!(
        n_single.status.success(),
        "-n single-sample failed: {}",
        String::from_utf8_lossy(&n_single.stderr)
    );
    assert!(
        String::from_utf8_lossy(&n_single.stdout)
            .contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t42\n"),
        "-n must extract sample '42' into VCF header"
    );

    let paired_n_arg = format!("{}|fixtures/micro/micro.bam", named_bam.display());
    let n_paired_var = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            &paired_n_arg,
            "-R",
            "chr1:40-260",
            "-n",
            r"/sample_(\d+)_run/",
            "--emit-intermediate",
            "--threads=1",
        ])
        .output()
        .expect("spawn -n paired -R .var");
    assert!(
        n_paired_var.status.success(),
        "-n paired -R .var failed: {}",
        String::from_utf8_lossy(&n_paired_var.stderr)
    );
    let n_paired_stdout = String::from_utf8_lossy(&n_paired_var.stdout);
    assert!(
        n_paired_stdout.lines().any(|l| l.starts_with("42\t")),
        "-n paired -R .var must resolve sample column to '42', got:\n{n_paired_stdout}"
    );

    let cohort_dir = tmp.path().join("cohort").join("SAMP_ESC");
    std::fs::create_dir_all(&cohort_dir).expect("mkdir cohort/SAMP_ESC");
    let cohort_bam = cohort_dir.join("aln.bam");
    let cohort_bai = cohort_dir.join("aln.bam.bai");
    std::fs::copy(
        workspace_root().join("fixtures/micro/micro.bam"),
        &cohort_bam,
    )
    .expect("copy micro.bam -> cohort/SAMP_ESC/aln.bam");
    std::fs::copy(
        workspace_root().join("fixtures/micro/micro.bam.bai"),
        &cohort_bai,
    )
    .expect("copy micro.bam.bai -> cohort/SAMP_ESC/aln.bam.bai");
    for slash_pat in [r"/cohort\/([^/]+)\//", r"//cohort/([^/]+)//"] {
        let slash_single = base_cmd()
            .args(["-G", "fixtures/micro/micro.ref.fa", "-b"])
            .arg(&cohort_bam)
            .args(["-R", "chr1:40-260", "-n", slash_pat, "--threads=1"])
            .output()
            .expect("spawn -n slash pattern VCF");
        assert!(
            slash_single.status.success(),
            "-n {slash_pat} failed: {}",
            String::from_utf8_lossy(&slash_single.stderr)
        );
        assert!(
            String::from_utf8_lossy(&slash_single.stdout)
                .contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP_ESC\n"),
            "-n {slash_pat} must extract SAMP_ESC without double-stripping slashes"
        );
    }

    let uz0_out = base_cmd()
        .args(base_region_args)
        .args(["-uz", "0"])
        .output()
        .expect("spawn clustered -uz 0");
    assert!(
        uz0_out.status.success(),
        "clustered -uz 0 failed: {}",
        String::from_utf8_lossy(&uz0_out.stderr)
    );
    assert!(String::from_utf8_lossy(&uz0_out.stdout).contains("chr1\t101\t"));

    let vs_lenient = base_cmd()
        .args(base_region_args)
        .args(["-VS", "LENIENT"])
        .output()
        .expect("spawn -VS LENIENT");
    assert!(
        vs_lenient.status.success(),
        "-VS LENIENT failed: {}",
        String::from_utf8_lossy(&vs_lenient.stderr)
    );
    assert_eq!(vs_lenient.stdout, default_out.stdout);

    let vs_strict = base_cmd()
        .args(base_region_args)
        .args(["-VS", "STRICT"])
        .output()
        .expect("spawn -VS STRICT");
    assert!(!vs_strict.status.success());
    assert!(String::from_utf8_lossy(&vs_strict.stderr)
        .contains("-VS STRICT not implemented in v1.0; use VarDictJava"));

    assert!(
        String::from_utf8_lossy(&default_out.stdout).contains(";END=101;"),
        "default -R VCF without -E/--no-end must emit ;END=101;"
    );
    let no_end_out = base_cmd()
        .args(base_region_args)
        .arg("--no-end")
        .output()
        .expect("spawn --no-end");
    assert!(
        no_end_out.status.success(),
        "--no-end failed: {}",
        String::from_utf8_lossy(&no_end_out.stderr)
    );
    let no_end_vcf = String::from_utf8_lossy(&no_end_out.stdout);
    assert!(no_end_vcf.contains("chr1\t101\t"));
    assert!(
        !no_end_vcf.contains(";END="),
        "--no-end must omit ;END= from VCF INFO column, got:\n{no_end_vcf}"
    );

    let rayon_stats = tmp.path().join("rayon_threads.json");
    let rayon_run = base_cmd()
        .env("RAYON_NUM_THREADS", "2")
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--stats-json",
        ])
        .arg(&rayon_stats)
        .output()
        .expect("spawn with RAYON_NUM_THREADS=2");
    assert!(
        rayon_run.status.success(),
        "RAYON_NUM_THREADS=2 run failed: {}",
        String::from_utf8_lossy(&rayon_run.stderr)
    );
    let rayon_stats_text = std::fs::read_to_string(&rayon_stats).expect("read rayon_threads.json");
    assert!(
        rayon_stats_text.contains("\"threads\": 2"),
        "expected RAYON_NUM_THREADS=2 to resolve threads=2 in live binary, got:\n{rayon_stats_text}"
    );
}
