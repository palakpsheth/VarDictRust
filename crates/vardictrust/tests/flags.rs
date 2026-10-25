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
    if !workspace_root().join(ref_fa).exists() {
        let status = Command::new("python3")
            .arg("harness/convert_fastas.py")
            .current_dir(workspace_root())
            .status()
            .expect("convert fastas");
        assert!(status.success(), "convert_fastas.py failed");
    }
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
    // 1. On micro.bam, test alias equivalence between -t and --dedup
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

    // 2. Behavioral verification on fixtures/pr402 (which contains 31 duplicate reads at chr7:1):
    // Without -t: all 30 duplicate reads contribute to coverage; variant chr7:30 A>C is called with DP=30, VD=30.
    let out_nodedup = base_cmd()
        .args([
            "-G",
            "fixtures/pr402/ref.fa",
            "-b",
            "fixtures/pr402/pr402.bam",
            "-R",
            "chr7:1-140",
            "-f",
            "0.01",
            "-N",
            "test",
        ])
        .output()
        .expect("run without -t on pr402");
    assert!(out_nodedup.status.success());
    let variants_nodedup = extract_variants(&out_nodedup.stdout);
    assert_eq!(
        variants_nodedup.len(),
        1,
        "without deduplication, variant at chr7:30 should be called with full depth"
    );
    assert!(
        variants_nodedup[0].contains("DP=30;END=30;VD=30"),
        "expected DP=30;VD=30 without -t, got: {}",
        variants_nodedup[0]
    );

    // With -t: all 30 duplicates are filtered by removeDuplicatedReads, leaving 1 read; variant drops below -r 2.
    let out_dedup = base_cmd()
        .args([
            "-G",
            "fixtures/pr402/ref.fa",
            "-b",
            "fixtures/pr402/pr402.bam",
            "-R",
            "chr7:1-140",
            "-f",
            "0.01",
            "-N",
            "test",
            "-t",
        ])
        .output()
        .expect("run with -t on pr402");
    assert!(out_dedup.status.success());
    let variants_dedup = extract_variants(&out_dedup.stdout);
    assert_eq!(
        variants_dedup.len(),
        0,
        "with -t, all duplicate reads must be removed and 0 variants called"
    );

    // With --dedup: must also filter duplicate reads identically
    let out_dedup_long = base_cmd()
        .args([
            "-G",
            "fixtures/pr402/ref.fa",
            "-b",
            "fixtures/pr402/pr402.bam",
            "-R",
            "chr7:1-140",
            "-f",
            "0.01",
            "-N",
            "test",
            "--dedup",
        ])
        .output()
        .expect("run with --dedup on pr402");
    assert!(out_dedup_long.status.success());
    assert_eq!(extract_variants(&out_dedup_long.stdout).len(), 0);
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
    let tmp_bed_a = workspace_root().join("scratch").join("test_flags_se_a.bed");
    let tmp_bed_b = workspace_root().join("scratch").join("test_flags_se_b.bed");
    fs::create_dir_all(tmp_bed_a.parent().unwrap()).unwrap();

    // Bed A: region bounds (cols 2 & 3) = 40..120; segment bounds (cols 10 & 11) = 40..80.
    // chr1:101 is OUTSIDE segment [40, 80], so zero variants must be called.
    let bed_content_a = "chr1\t40\t120\tmicro_chr1\t0\t+\t40\t120\t0\t40\t80\tmicro_chr1\n";
    fs::write(&tmp_bed_a, bed_content_a).unwrap();

    // Bed B: region bounds (cols 2 & 3) = 40..120; segment bounds (cols 10 & 11) = 90..120.
    // chr1:101 is INSIDE segment [90, 120], so variant chr1:101 T>A must be called.
    let bed_content_b = "chr1\t40\t120\tmicro_chr1\t0\t+\t40\t120\t0\t90\t120\tmicro_chr1\n";
    fs::write(&tmp_bed_b, bed_content_b).unwrap();

    let out_a = base_cmd()
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
            "-x",
            "0",
            "-z",
            tmp_bed_a.to_str().unwrap(),
        ])
        .output()
        .expect("run with segment 40-80");

    let out_b = base_cmd()
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
            "-x",
            "0",
            "-z",
            tmp_bed_b.to_str().unwrap(),
        ])
        .output()
        .expect("run with segment 90-120");

    let _ = fs::remove_file(&tmp_bed_a);
    let _ = fs::remove_file(&tmp_bed_b);

    assert!(out_a.status.success());
    assert!(out_b.status.success());

    let vars_a = extract_variants(&out_a.stdout);
    let vars_b = extract_variants(&out_b.stdout);

    assert_eq!(
        vars_a.len(),
        0,
        "segment 40..80 should exclude variant at chr1:101"
    );
    assert_eq!(
        vars_b.len(),
        1,
        "segment 90..120 should include variant at chr1:101"
    );
    assert!(
        vars_b[0].contains("chr1\t101\t"),
        "expected variant at chr1:101, got: {}",
        vars_b[0]
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
            "-k",
            "0",
            "-o",
            "1.5",
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
    for th_arg in &["-th", "--th", "-th=2", "--threads=2", "--threads 2"] {
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

#[test]
fn test_amplicon_flag_a() {
    use std::io::Write;
    let mut amp_bed = tempfile::NamedTempFile::new().expect("temp amplicon bed");
    writeln!(amp_bed, "chr1\t81\t140\tmicro_amp\t0\t+\t81\t140").expect("write bed");
    let bed_path = amp_bed.path().to_str().expect("bed path");

    let out = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-z",
            bed_path,
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-a",
            "30:0.5",
            "-f",
            "0.01",
            "-r",
            "1",
        ])
        .output()
        .expect("run amplicon mode");

    assert!(
        out.status.success(),
        "amplicon run must succeed, stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("GDAMP=1") || stdout.contains("TLAMP=1"),
        "amplicon output must contain amplicon annotations, got:\n{stdout}"
    );
}

#[test]
fn test_qratio_flag_o() {
    let out_pass = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-f",
            "0.01",
            "-o",
            "0.0",
        ])
        .output()
        .expect("run with -o 0.0");
    assert!(out_pass.status.success());
    let variants_pass = extract_variants(&out_pass.stdout);
    assert!(!variants_pass.is_empty(), "-o 0.0 must allow variant calls");

    let out_filter = base_cmd()
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-R",
            "chr1:40-260",
            "-f",
            "0.01",
            "-o",
            "100.0",
        ])
        .output()
        .expect("run with -o 100.0");
    assert!(out_filter.status.success());
    let variants_filter = extract_variants(&out_filter.stdout);
    assert!(
        variants_filter.is_empty(),
        "-o 100.0 must filter out variants with hq_lq_ratio < 100.0"
    );
}

fn ensure_fastas_converted() {
    let ref_fa = workspace_root().join("VarDictJava/testdata/fastas/hg19.fa");
    let grch37_fa = workspace_root().join("VarDictJava/testdata/fastas/GRCh37.fa");
    if !ref_fa.exists() || !grch37_fa.exists() {
        let status = Command::new("python3")
            .arg("harness/convert_fastas.py")
            .current_dir(workspace_root())
            .status()
            .expect("convert fastas");
        assert!(status.success(), "convert_fastas.py failed");
    }
}

#[test]
fn test_cram_soft_masked_lowercase_ref_and_omitted_nm_matches_bam() {
    // Real biological dataset regression test for Bug C + mate-window .bai fallback:
    // Hartwig COLO829 melanoma WGS (`Colo829.chr20.CRISPR.bam`, `chr20:50454281-50455444`)
    // aligned against `VarDictJava/testdata/fastas/hg19.fa` (which contains soft-masked
    // lowercase reference bases `a/c/g/t`). Encodes the real BAM into an external-reference
    // CRAM 3.0 with `NM` tags stripped (`store_nm=0`), and verifies 100% byte-identical
    // VCF output between the BAM and CRAM (including `NM=` INFO values and `NM5.25` filter).
    ensure_fastas_converted();
    let ref_fa = "VarDictJava/testdata/fastas/hg19.fa";
    let bam = "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.CRISPR.bam";
    let region = "chr20:50454281-50455444";

    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let sam_out = Command::new("samtools")
        .args(["view", "-h", bam])
        .current_dir(workspace_root())
        .output()
        .expect("samtools view -h");
    assert!(sam_out.status.success(), "samtools view -h failed");

    // Normalize @SQ SN:chr20 LN to match hg19.fa.fai (59858789) so samtools writes an
    // external-reference CRAM (`embed_ref=0`) and strip all `NM:i:*` tags so noodles-cram
    // must reconstruct `NM` via `compute_cram_nm` against soft-masked lowercase `hg19.fa`.
    let sam_text = String::from_utf8(sam_out.stdout).expect("utf8 sam");
    let mut cleaned_sam = String::with_capacity(sam_text.len());
    for line in sam_text.lines() {
        if line.starts_with("@SQ\tSN:chr20\t") {
            cleaned_sam.push_str("@SQ\tSN:chr20\tLN:59858789\n");
        } else if line.starts_with("@SQ\t") {
            continue;
        } else if line.starts_with('@') {
            cleaned_sam.push_str(line);
            cleaned_sam.push('\n');
        } else {
            let mut cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 6 && cols[6] != "=" && cols[6] != "chr20" {
                cols[6] = "*";
                cols[7] = "0";
            }
            let filtered: Vec<&str> = cols
                .into_iter()
                .filter(|c| !c.starts_with("NM:i:"))
                .collect();
            cleaned_sam.push_str(&filtered.join("\t"));
            cleaned_sam.push('\n');
        }
    }
    let sam_path = tmp_dir.path().join("crispr_no_nm.sam");
    fs::write(&sam_path, cleaned_sam).expect("write sam");

    let cram_path = tmp_dir.path().join("crispr.cram");
    let enc_status = Command::new("samtools")
        .args([
            "view",
            "-C",
            "-T",
            ref_fa,
            "-o",
            cram_path.to_str().unwrap(),
            sam_path.to_str().unwrap(),
        ])
        .current_dir(workspace_root())
        .status()
        .expect("samtools view -C");
    assert!(enc_status.success(), "samtools CRAM encode failed");

    let idx_status = Command::new("samtools")
        .args(["index", cram_path.to_str().unwrap()])
        .current_dir(workspace_root())
        .status()
        .expect("samtools index cram");
    assert!(idx_status.success(), "samtools index cram failed");

    let out_bam = base_cmd()
        .args([
            "-G", ref_fa, "-b", bam, "-N", "COLO829", "-R", region, "-f", "0.001",
        ])
        .output()
        .expect("run on real BAM");
    assert!(out_bam.status.success());

    let out_cram = base_cmd()
        .args([
            "-G",
            ref_fa,
            "-b",
            cram_path.to_str().unwrap(),
            "-N",
            "COLO829",
            "-R",
            region,
            "-f",
            "0.001",
        ])
        .output()
        .expect("run on soft-masked external-ref CRAM");
    assert!(out_cram.status.success());

    let bam_vars = extract_variants(&out_bam.stdout);
    let cram_vars = extract_variants(&out_cram.stdout);
    assert_eq!(bam_vars.len(), 8, "expected 8 real COLO829 variants");
    assert_eq!(
        bam_vars, cram_vars,
        "soft-masked hg19 CRAM with omitted NM must produce byte-identical VCF variants to BAM"
    );

    // Also verify explicit multi-tile (`--tile-size 400`) execution on the same real BAM:
    let out_tiled = base_cmd()
        .args([
            "-G",
            ref_fa,
            "-b",
            bam,
            "-N",
            "COLO829",
            "-R",
            region,
            "-f",
            "0.001",
            "--tile-size",
            "400",
        ])
        .output()
        .expect("run tiled on real BAM");
    assert!(
        out_tiled.status.success(),
        "tiled run failed: {}",
        String::from_utf8_lossy(&out_tiled.stderr)
    );
    assert_eq!(
        bam_vars,
        extract_variants(&out_tiled.stdout),
        "tiled run (--tile-size 400) on real COLO829 BAM must match untiled output"
    );
}

#[test]
fn test_freq_zero_perl_falsy_fallback_vs_truthy_zero_decimal_end_to_end() {
    // Real biological dataset regression test for Bug D (`var2vcf_valid.pl:21`):
    // Runs on real GIAB NA24385 deep WGS (`NA24385_chr_20_20156208_20234556.bam`,
    // `20:20156208-20157500` against `GRCh37.fa`) with `-f 0` vs `-f 0.0`.
    // Under `-f 0`, Perl treats `"0"` as falsy and falls back to `0.02` (`##FILTER=<ID=f0.02>`,
    // tagging sub-2% AF variants with `f0.02` and `GT=0/0`), whereas `-f 0.0` is truthy (`0.0`).
    ensure_fastas_converted();
    let ref_fa = "VarDictJava/testdata/fastas/GRCh37.fa";
    let bam = "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_20_20156208_20234556.bam";
    let region = "20:20156208-20157500";

    let out_f0 = base_cmd()
        .args(["-G", ref_fa, "-b", bam, "-R", region, "-f", "0"])
        .output()
        .expect("run with -f 0");
    assert!(out_f0.status.success());
    let stdout_f0 = String::from_utf8_lossy(&out_f0.stdout);

    let out_f0_dec = base_cmd()
        .args(["-G", ref_fa, "-b", bam, "-R", region, "-f", "0.0"])
        .output()
        .expect("run with -f 0.0");
    assert!(out_f0_dec.status.success());
    let stdout_f0_dec = String::from_utf8_lossy(&out_f0_dec.stdout);

    assert!(
        stdout_f0.contains("##FILTER=<ID=f0.02,Description=\"Allele frequency < 0.02\">"),
        "-f 0 header must emit f0.02 filter description"
    );
    assert!(
        stdout_f0_dec.contains("##FILTER=<ID=f0,Description=\"Allele frequency < 0\">"),
        "-f 0.0 header must emit f0 filter description"
    );

    let vars_f0 = extract_variants(&out_f0.stdout);
    let vars_f0_dec = extract_variants(&out_f0_dec.stdout);
    assert_eq!(vars_f0.len(), vars_f0_dec.len());
    assert!(
        vars_f0
            .iter()
            .any(|v| v.contains("f0.02") && v.contains("\t0/0:")),
        "real sub-2% AF candidate in NA24385 must receive f0.02 filter and GT=0/0 under -f 0"
    );
    assert!(
        !vars_f0_dec.iter().any(|v| v.contains("f0.02")),
        "-f 0.0 must not apply f0.02 filter to sub-2% AF candidates"
    );
}

#[test]
fn test_fisher_flag_java_mode_vcf_and_intermediate_table_end_to_end() {
    // Real biological dataset regression test for Bug H (`FisherExact.java` `--fisher` mode):
    // Runs on real paired tumor/normal BRAF V600E hotspot (`ArgumentNA_Tumor.bam|ArgumentNA_Normal.bam`,
    // `chr7:140453010-140453158` against `hg19.fa`) with `--fisher` in both 61-column `.var`
    // intermediate mode (comparing directly against the checked-in VarDictJava v1.8.3 oracle
    // fixture) and end-to-end VCF mode.
    ensure_fastas_converted();
    let ref_fa = "VarDictJava/testdata/fastas/hg19.fa";
    let tumor = "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/ArgumentNA_Tumor.bam";
    let normal = "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/ArgumentNA_Normal.bam";
    let bam_pair = format!("{tumor}|{normal}");
    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let bed_path = tmp_dir.path().join("argument_na.bed");
    fs::write(&bed_path, "chr7\t140453010\t140453158\ttestbed\n").expect("write bed");

    let out_var = base_cmd()
        .args([
            "-G",
            ref_fa,
            "-b",
            &bam_pair,
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-q",
            "25",
            "-F",
            "0x100",
            "-f",
            "0.001",
            "-k",
            "1",
            "-M",
            "20",
            "-p",
            "--fisher",
            "--emit-intermediate",
            "-z",
            bed_path.to_str().unwrap(),
        ])
        .output()
        .expect("run somatic --fisher --emit-intermediate");
    assert!(out_var.status.success());

    let expected_path = workspace_root().join(
        "VarDictJava/testdata/integrationtestcases/Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p-fisher;.txt",
    );
    let expected_raw = fs::read_to_string(&expected_path).expect("read oracle fixture");
    // Normalize the sample-name column (col 0) just as `run_integration.sh` does.
    let normalize_rows = |text: &str| -> Vec<String> {
        text.lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let mut cols: Vec<&str> = l.split('\t').collect();
                if !cols.is_empty() {
                    cols[0] = "SAMPLE";
                }
                assert_eq!(cols.len(), 61, "somatic --fisher must emit 61 columns");
                cols.join("\t")
            })
            .collect()
    };
    let actual_rows = normalize_rows(&String::from_utf8_lossy(&out_var.stdout));
    let expected_body: String = expected_raw
        .lines()
        .skip(1)
        .collect::<Vec<&str>>()
        .join("\n");
    let expected_rows = normalize_rows(&expected_body);
    assert_eq!(
        actual_rows.len(),
        expected_rows.len(),
        "row count must match VarDictJava v1.8.3 oracle fixture (149 pileup positions)"
    );
    assert_eq!(
        actual_rows, expected_rows,
        "somatic --fisher 61-column intermediate output must match VarDictJava v1.8.3 oracle byte-for-byte"
    );

    // Also verify end-to-end somatic VCF output with `--fisher` on the same real tumor/normal BAMs:
    let out_vcf = base_cmd()
        .args([
            "-G",
            ref_fa,
            "-b",
            &bam_pair,
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-q",
            "25",
            "-F",
            "0x100",
            "-f",
            "0.001",
            "-k",
            "1",
            "-M",
            "20",
            "--fisher",
            "-z",
            bed_path.to_str().unwrap(),
        ])
        .output()
        .expect("run somatic --fisher VCF");
    assert!(out_vcf.status.success());
    let vcf_vars = extract_variants(&out_vcf.stdout);
    assert!(
        !vcf_vars.is_empty()
            && vcf_vars
                .iter()
                .all(|v| v.contains(";SSF=") && v.contains(";SOR=")),
        "somatic --fisher VCF records must emit SSF and SOR INFO fields"
    );
}

#[test]
fn test_real_sv_duplication_and_5prime_deletion_partial_pipeline_end_to_end() {
    // Real biological dataset regression test for Bug A (`CigarParser::parse_cigar` secondary
    // `LoadedRegion` selection during SV `partial_pipeline`), Bug B (`get_reference` `siteEnd = len - 17`
    // truncation), Bug F (`realignins` `<dupN>` 3' flank `substr` offset), and Bug G (`find_sv`
    // 5' deletion `p5_adj + 1` coverage lookup) on real Hartwig COLO829 melanoma WGS BAMs
    // (`Colo829.chr20.MarkingDUPSV.bam` and `Colo829.chr20.CandDel5SV.bam` against `hg19.fa`),
    // verified byte-for-byte against the checked-in VarDictJava v1.8.3 oracle fixtures.
    ensure_fastas_converted();
    let ref_fa = "VarDictJava/testdata/fastas/hg19.fa";
    let tmp_dir = tempfile::tempdir().expect("tempdir");

    let cases = [
        (
            "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.MarkingDUPSV.bam",
            "chr20\t47131954\t47132045\ttestbed\n",
            "VarDictJava/testdata/integrationtestcases/Simple;hg19.fa;Colo829.chr20.MarkingDUPSV.bam;chr20;47131954-47132045;-f 0.001;.txt",
        ),
        (
            "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.CandDel5SV.bam",
            "chr20\t506503\t507295\ttestbed\n",
            "VarDictJava/testdata/integrationtestcases/Simple;hg19.fa;Colo829.chr20.CandDel5SV.bam;chr20;506503-507295;-f 0.001;.txt",
        ),
    ];

    let normalize_rows = |text: &str| -> Vec<String> {
        text.lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let mut cols: Vec<&str> = l.split('\t').collect();
                if !cols.is_empty() {
                    cols[0] = "SAMPLE";
                }
                assert_eq!(cols.len(), 36, "Simple .var must emit 36 columns");
                cols.join("\t")
            })
            .collect()
    };

    for (idx, (bam, bed_line, oracle_rel)) in cases.iter().enumerate() {
        let bed_path = tmp_dir.path().join(format!("sv_case_{idx}.bed"));
        fs::write(&bed_path, bed_line).expect("write bed");

        let out = base_cmd()
            .args([
                "-G",
                ref_fa,
                "-b",
                bam,
                "-c",
                "1",
                "-S",
                "2",
                "-E",
                "3",
                "-g",
                "4",
                "-f",
                "0.001",
                "--emit-intermediate",
                "-z",
                bed_path.to_str().unwrap(),
            ])
            .output()
            .expect("run real SV fixture");
        assert!(out.status.success(), "SV fixture {bam} failed");

        let expected_raw =
            fs::read_to_string(workspace_root().join(oracle_rel)).expect("read oracle");
        let expected_body: String = expected_raw
            .lines()
            .skip(1)
            .collect::<Vec<&str>>()
            .join("\n");
        let mut actual_rows = normalize_rows(&String::from_utf8_lossy(&out.stdout));
        let mut expected_rows = normalize_rows(&expected_body);
        actual_rows.sort();
        expected_rows.sort();
        assert_eq!(
            actual_rows, expected_rows,
            "real SV fixture {bam} must match VarDictJava v1.8.3 oracle byte-for-byte"
        );
    }
}
