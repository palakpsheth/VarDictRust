use std::fs;
use std::path::Path;
use std::process::Command;

fn workspace_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
}

fn vardictrust_bin() -> &'static str {
    env!("CARGO_BIN_EXE_vardictrust")
}

fn base_cmd() -> Command {
    let mut cmd = Command::new(vardictrust_bin());
    cmd.current_dir(workspace_root());
    cmd
}

fn extract_variants(stdout: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(stdout)
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|s| s.to_string())
        .collect()
}

#[test]
fn test_unsupported_flags_rejected() {
    let unsupported = ["-n", "-C", "-T", "--trim", "--VS", "-Z", "--downsample"];

    for flag in unsupported {
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
            "flag {flag} must be rejected loudly with non-zero exit code"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("not implemented in v1.0; use VarDictJava"),
            "flag {flag} stderr must mention v1.0 rejection, got: {stderr}"
        );
    }
}

#[test]
fn test_dp_err_routes_to_stderr() {
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
            "--DP",
            "ERR",
        ])
        .output()
        .expect("run vardictrust with --DP ERR");

    if !out.status.success() {
        panic!(
            "command failed: {}\nstdout: {}\nstderr: {}",
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    assert!(out.status.success());
    assert!(
        out.stdout.is_empty(),
        "stdout must be empty under --DP ERR, but got {} bytes",
        out.stdout.len()
    );
    assert!(
        !out.stderr.is_empty(),
        "stderr must contain variant output under --DP ERR, but was empty"
    );
}

#[test]
fn test_default_printer_long_alias_routes_to_stderr() {
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
            "--defaultPrinter",
            "ERR",
        ])
        .output()
        .expect("run vardictrust with --defaultPrinter ERR");

    if !out.status.success() {
        panic!(
            "command failed: {}\nstdout: {}\nstderr: {}",
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    assert!(out.status.success());
    assert!(
        out.stdout.is_empty(),
        "stdout must be empty under --defaultPrinter ERR, but got {} bytes",
        out.stdout.len()
    );
    assert!(
        !out.stderr.is_empty(),
        "stderr must contain variant output under --defaultPrinter ERR, but was empty"
    );
}

#[test]
fn test_parity_exempt_v_flag() {
    let out_no_v = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
        ])
        .output()
        .expect("run without -v");
    assert!(out_no_v.status.success());

    let out_with_v = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-v",
        ])
        .output()
        .expect("run with -v");
    assert!(out_with_v.status.success());
    assert_eq!(
        out_no_v.stdout, out_with_v.stdout,
        "-v output must be identical to without -v (parity-exempt ignored flag)"
    );
}

#[test]
fn test_minmatch_flag_m() {
    // Reads in micro.bam are ~80bp.
    // -M 200 should filter out all reads.
    let out_filtered = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-M",
            "200",
        ])
        .output()
        .expect("run with -M 200");
    assert!(out_filtered.status.success());
    let variants_filtered = extract_variants(&out_filtered.stdout);
    assert!(
        variants_filtered.is_empty(),
        "-M 200 should filter out reads and produce 0 variants"
    );

    // -M 50 should retain reads and emit variant.
    let out_pass = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-M",
            "50",
        ])
        .output()
        .expect("run with -M 50");
    assert!(out_pass.status.success());
    let variants_pass = extract_variants(&out_pass.stdout);
    assert!(
        !variants_pass.is_empty(),
        "-M 50 should allow reads and produce variants"
    );
}

#[test]
fn test_out_mapq_flag_cap_o() {
    // Reads in micro.bam have MAPQ 60.
    // VarDictJava's isGoodVar checks `meanMappingQuality < conf.mapq` on candidates with AF <= 0.30.
    // Test on NA24385 integration fixture over 800bp window: -O 65 filters out low-MAPQ candidates.
    let ref_fa = "VarDictJava/testdata/fastas/GRCh37.fa";
    let bam = "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_20_20156208_20234556.bam";
    let region = "20:20156208-20157000";

    let out_unfiltered = base_cmd()
        .args(["-G", ref_fa, "-b", bam, "-R", region, "-O", "0"])
        .output()
        .expect("run with -O 0");
    assert!(out_unfiltered.status.success());
    let variants_unfiltered = extract_variants(&out_unfiltered.stdout);

    let out_filtered = base_cmd()
        .args(["-G", ref_fa, "-b", bam, "-R", region, "-O", "65"])
        .output()
        .expect("run with -O 65");
    assert!(out_filtered.status.success());
    let variants_filtered = extract_variants(&out_filtered.stdout);

    assert!(
        variants_filtered.len() < variants_unfiltered.len(),
        "-O 65 should filter candidates relative to -O 0 (got {} vs {})",
        variants_filtered.len(),
        variants_unfiltered.len()
    );
    assert!(
        !variants_filtered.is_empty(),
        "-O 65 should still retain high-MAPQ variants"
    );
}

#[test]
fn test_move_indels_to_3_prime_flag_3() {
    let out = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr2:110-230",
            "-3",
        ])
        .output()
        .expect("run with -3");
    assert!(out.status.success());
    let variants = extract_variants(&out.stdout);
    assert!(
        !variants.is_empty(),
        "-3 should call variants without error"
    );
}

#[test]
fn test_extension_flag_x() {
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
            "-x",
            "50",
        ])
        .output()
        .expect("run with -x 50");
    assert!(out.status.success());
    let variants = extract_variants(&out.stdout);
    assert!(
        !variants.is_empty(),
        "-x 50 should expand regions and call variants"
    );
}

#[test]
fn test_deldupvar_flag() {
    let out = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--deldupvar",
        ])
        .output()
        .expect("run with --deldupvar");
    assert!(out.status.success());
}

#[test]
fn test_dedup_flag_t() {
    let out_short = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-t",
        ])
        .output()
        .expect("run with -t");
    assert!(out_short.status.success());

    let out_long = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--dedup",
        ])
        .output()
        .expect("run with --dedup");
    assert!(out_long.status.success());
    assert_eq!(out_short.stdout, out_long.stdout);
}

#[test]
fn test_verbose_flag_y() {
    let out_short = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-y",
        ])
        .output()
        .expect("run with -y");
    assert!(out_short.status.success());

    let out_long = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--verbose",
        ])
        .output()
        .expect("run with --verbose");
    assert!(out_long.status.success());
}

#[test]
fn test_seg_start_and_end_flags_s_e() {
    let tmp_bed = workspace_root().join("scratch").join("test_flags_se.bed");
    fs::create_dir_all(tmp_bed.parent().unwrap()).unwrap();
    let bed_content = "chr1\t40\t260\tmicro_chr1\t0\t+\t40\t260\t0\t40\t260\tmicro_chr1\n";
    fs::write(&tmp_bed, bed_content).unwrap();

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
            "-s",
            "10",
            "-e",
            "11",
            "-g",
            "4",
            "-z",
            tmp_bed.to_str().unwrap(),
        ])
        .output()
        .expect("run with -s 10 -e 11");

    let _ = fs::remove_file(&tmp_bed);
    assert!(out.status.success());
    let variants = extract_variants(&out.stdout);
    assert!(
        !variants.is_empty(),
        "-s and -e flags should parse segment bounds and emit variants"
    );
}

#[test]
fn test_functional_operational_flags_suite() {
    // Tests --emit-intermediate, --fisher, -u, --UN, -K, -p, -h, -D, -th, etc.
    let out = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--emit-intermediate",
            "--fisher",
            "-h",
            "-D",
            "-u",
            "-K",
            "-f",
            "0.01",
            "-r",
            "2",
            "-B",
            "1",
            "-q",
            "20",
            "-Q",
            "10",
            "-m",
            "5",
            "-P",
            "3",
            "-X",
            "2",
            "-Y",
            "1000",
            "-I",
            "40",
            "-F",
            "0x504",
            "-w",
            "250",
            "-W",
            "80",
            "-A",
            "2",
            "-L",
            "500",
            "-V",
            "0.02",
            "--nosv",
            "--chimeric",
            "--mfreq",
            "0.2",
            "--nmfreq",
            "0.08",
            "-J",
            "100",
            "-j",
            "5",
            "-d",
            "\t",
            "-N",
            "CUSTOM_SAMPLE",
            "--tile-size",
            "500",
            "--adaptor",
            "AGATCGGAAGAGC",
            "-th",
            "2",
        ])
        .output()
        .expect("run complex operational flag combination");

    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert!(!lines.is_empty(), "must emit table lines");
    assert!(
        lines[0].starts_with("Sample\t"),
        "with -h and --emit-intermediate, first line must be column header"
    );
    assert!(
        stdout.contains("CUSTOM_SAMPLE"),
        "output must contain custom sample name passed via -N"
    );
}

#[test]
fn test_unique_un_flag() {
    let out = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "--UN",
        ])
        .output()
        .expect("run with --UN");
    assert!(out.status.success());
}

#[test]
fn test_pileup_mode_p() {
    let out = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:100-105",
            "-p",
        ])
        .output()
        .expect("run with -p");
    assert!(out.status.success());
}

#[test]
fn test_splice_mode_i() {
    let out = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-i",
        ])
        .output()
        .expect("run with -i");
    assert!(out.status.success());
}

#[test]
fn test_threads_normalization_forms() {
    for th_arg in &["-th", "--th", "-th=2", "--threads=2", "--threads", "2"] {
        let mut args = vec![
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:100-110",
        ];
        if th_arg.contains(' ') {
            let parts: Vec<&str> = th_arg.split_whitespace().collect();
            args.extend(parts);
        } else {
            args.push(th_arg);
        }

        let out = base_cmd()
            .args(&args)
            .output()
            .unwrap_or_else(|e| panic!("run with thread arg {th_arg}: {e}"));
        assert!(
            out.status.success(),
            "thread arg form {th_arg} must succeed, stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
}
