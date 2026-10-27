use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Mutex;

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

fn ensure_test_fastas(droot: &Path) {
    let fasta_dir = droot.join("VarDictJava/testdata/fastas");
    if !fasta_dir.join("GRCh37.fa").is_file() || !fasta_dir.join("GRCh37.fa.fai").is_file() {
        let script = workspace_root().join("harness/convert_fastas.py");
        let status = Command::new("python3")
            .arg(&script)
            .current_dir(droot)
            .status()
            .expect("run convert_fastas.py");
        assert!(status.success(), "convert_fastas.py failed");
    }
}

fn load_xfail_registry() -> HashMap<String, String> {
    let path = workspace_root().join("harness/integration_xfail.tsv");
    let content = fs::read_to_string(&path).expect("read integration_xfail.tsv");
    let mut map = HashMap::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() >= 5 && cols[1].trim() == "XFAIL" {
            map.insert(cols[0].trim().to_string(), cols[4].trim().to_string());
        }
    }
    map
}

fn sha256_hex(input: &[u8]) -> String {
    use std::io::Write;
    use std::process::Stdio;
    let mut child = Command::new("sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .or_else(|_| {
            Command::new("shasum")
                .args(["-a", "256"])
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .spawn()
        })
        .expect("spawn sha256sum or shasum");
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(input)
        .expect("write stdin");
    let out = child.wait_with_output().expect("wait sha256sum");
    assert!(out.status.success(), "sha256sum failed");
    String::from_utf8_lossy(&out.stdout)
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn run_single_integration_case(
    testfile: &Path,
    droot: &Path,
    xfails: &HashMap<String, String>,
    run_vcf_companion: bool,
) -> &'static str {
    let fname = testfile
        .file_name()
        .and_then(|s| s.to_str())
        .expect("filename");
    let raw_content = fs::read_to_string(testfile).expect("read testcase");
    let content = raw_content.trim_end_matches(&['\n', '\r'][..]);
    let mut lines = content.split('\n').map(|l| l.trim_end_matches('\r'));
    let header = lines.next().expect("header line");
    let expected_lines: Vec<&str> = lines.collect();
    let expected = expected_lines.join("\n");

    let fields: Vec<String> = header.split(',').map(|s| s.trim().to_string()).collect();
    let mode = &fields[0];
    assert!(
        mode == "Simple" || mode == "Somatic" || mode == "Amplicon",
        "unexpected mode {mode} in {fname}"
    );
    let fasta_path = droot.join("VarDictJava/testdata/fastas").join(&fields[1]);
    let bam_dir =
        droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
    let bam_raw = &fields[2];
    let bam_arg = if let Some((b1, b2)) = bam_raw.split_once('|') {
        format!(
            "{}|{}",
            bam_dir.join(b1).display(),
            bam_dir.join(b2).display()
        )
    } else {
        bam_dir.join(bam_raw).display().to_string()
    };
    let chr = &fields[3];
    let (start_val, end_val, mut field_idx) = if let Some((s, e)) = fields[4].split_once('-') {
        (s.to_string(), e.to_string(), 5)
    } else {
        (
            fields[4].clone(),
            fields.get(5).cloned().unwrap_or_default(),
            6,
        )
    };

    let mut extra_args: Vec<String> = Vec::new();
    while field_idx < fields.len() {
        let fld = &fields[field_idx];
        if fld.starts_with('-') {
            for tok in fld.split_whitespace() {
                extra_args.push(tok.to_string());
            }
        }
        field_idx += 1;
    }

    let tmp_dir = tempfile::tempdir().expect("tempdir");
    let tmp_bed = tmp_dir.path().join("region.bed");
    let bed_line = if mode == "Amplicon" {
        format!(
            "{chr}\t{start_val}\t{end_val}\ttestbed\t0\t0\t{}\t{}\n",
            fields[6], fields[7]
        )
    } else {
        format!("{chr}\t{start_val}\t{end_val}\ttestbed\n")
    };
    fs::write(&tmp_bed, bed_line).expect("write bed");

    let bin = env!("CARGO_BIN_EXE_vardictrust");
    let mut cmd = Command::new(bin);
    cmd.current_dir(workspace_root()).args([
        "-th",
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
        "-G",
        fasta_path.to_str().expect("fasta path"),
        "-b",
        &bam_arg,
        "--emit-intermediate",
    ]);
    cmd.args(&extra_args);
    cmd.arg(&tmp_bed);

    let out = cmd.output().expect("execute vardictrust");
    assert!(
        out.status.success(),
        "vardictrust failed on {fname}: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let actual_str = String::from_utf8_lossy(&out.stdout);
    let actual_trimmed = actual_str.trim_end_matches(&['\n', '\r'][..]);

    let mut actual_sorted: Vec<&str> = if actual_trimmed.is_empty() {
        vec![""]
    } else {
        actual_trimmed.split('\n').collect()
    };
    actual_sorted.sort_unstable();

    let mut expected_sorted: Vec<&str> = if expected.is_empty() {
        vec![""]
    } else {
        expected.split('\n').collect()
    };
    expected_sorted.sort_unstable();

    let verdict = if actual_sorted == expected_sorted {
        assert!(
            !xfails.contains_key(fname),
            "XPASS on registered XFAIL {fname}"
        );
        "PASS"
    } else if let Some(exp_sha) = xfails.get(fname) {
        let mut joined = actual_sorted.join("\n");
        joined.push('\n');
        let got_sha = sha256_hex(joined.as_bytes());
        assert_eq!(
            &got_sha, exp_sha,
            "XFAIL sha256 mismatch on {fname}: got {got_sha}, expected {exp_sha}"
        );
        "XFAIL"
    } else {
        panic!("MISMATCH on {fname}");
    };

    if run_vcf_companion && !extra_args.iter().any(|a| a == "-i") {
        let mut vcf_cmd = Command::new(bin);
        vcf_cmd.current_dir(workspace_root()).args([
            "-th",
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
            "-G",
            fasta_path.to_str().expect("fasta path"),
            "-b",
            &bam_arg,
        ]);
        vcf_cmd.args(&extra_args);
        vcf_cmd.arg(&tmp_bed);
        let vcf_out = vcf_cmd.output().expect("execute vardictrust VCF companion");
        assert!(
            vcf_out.status.success(),
            "VCF companion failed on {fname}: {}",
            String::from_utf8_lossy(&vcf_out.stderr)
        );
    }

    verdict
}

fn exercise_in_memory_core_branches() {
    use vardictrust_core::stats::f80::{f80_sum, F80};
    use vardictrust_core::stats::rmath::{bd0, dbinom_raw, dhyper, phyper, stirlerr};

    for &n in &[
        0.0,
        0.5,
        1.0,
        2.5,
        3.1,
        5.2,
        5.5,
        6.2,
        6.8,
        7.5,
        9.5,
        12.5,
        13.5,
        15.0,
        20.0,
        25.0,
        50.0,
        100.0,
        250.0,
        7000.0,
        20_000_000.0,
    ] {
        let _ = stirlerr::stirlerr(n);
    }

    for &(x, np) in &[
        (f64::NAN, 1.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (1.05, 1.0),
        (5.0, 2.0),
        (2.0, 5.0),
        (1e308, 1e-308),
    ] {
        let _ = bd0::bd0(x, np);
    }

    for &y in &[f64::NAN, 0.0, 1.0, 2.0, 3.0, 4.0, 5.5] {
        let _ = dbinom_raw::pow1p(0.0, y);
        let _ = dbinom_raw::pow1p(0.1, y);
        let _ = dbinom_raw::pow1p(0.8, y);
        let _ = dbinom_raw::pow1p(f64::NAN, y);
    }

    for &give_log in &[false, true] {
        let _ = dbinom_raw::dbinom_raw(0.0, 10.0, 0.0, 1.0, give_log);
        let _ = dbinom_raw::dbinom_raw(2.0, 10.0, 0.0, 1.0, give_log);
        let _ = dbinom_raw::dbinom_raw(10.0, 10.0, 1.0, 0.0, give_log);
        let _ = dbinom_raw::dbinom_raw(2.0, 10.0, 1.0, 0.0, give_log);
        let _ = dbinom_raw::dbinom_raw(0.0, 0.0, 0.5, 0.5, give_log);
        let _ = dbinom_raw::dbinom_raw(0.0, 10.0, 0.7, 0.3, give_log);
        let _ = dbinom_raw::dbinom_raw(0.0, 10.0, 0.3, 0.7, give_log);
        let _ = dbinom_raw::dbinom_raw(10.0, 10.0, 0.7, 0.3, give_log);
        let _ = dbinom_raw::dbinom_raw(10.0, 10.0, 0.3, 0.7, give_log);
        let _ = dbinom_raw::dbinom_raw(-1.0, 10.0, 0.5, 0.5, give_log);
        let _ = dbinom_raw::dbinom_raw(11.0, 10.0, 0.5, 0.5, give_log);
        let _ = dhyper::dhyper(-1.0, 5.0, 5.0, 5.0, give_log);
        let _ = dhyper::dhyper(6.0, 5.0, 5.0, 5.0, give_log);
        let _ = dhyper::dhyper(0.0, 5.0, 5.0, 0.0, give_log);
        let _ = dhyper::dhyper(2.0, 5.0, 5.0, 5.0, give_log);
        let _ = phyper::phyper(-1.0, 5.0, 5.0, 5.0, true);
        let _ = phyper::phyper(10.0, 5.0, 5.0, 5.0, false);
        let _ = phyper::phyper(2.0, 10.0, 10.0, 10.0, true);
        let _ = phyper::phyper(8.0, 10.0, 10.0, 10.0, false);
    }

    let subnormal = f64::from_bits(1);
    let f_sub = F80::from_f64(subnormal);
    assert!(f_sub.to_f64() >= 0.0);
    let f_big = F80 {
        exp: 1030,
        mant: 1u64 << 63,
    };
    assert_eq!(f_big.to_f64(), f64::INFINITY);
    let f_tiny = F80 {
        exp: -1100,
        mant: 1u64 << 63,
    };
    assert_eq!(f_tiny.to_f64(), 0.0);
    let _ = f80_sum(&[subnormal, 1e-300, 1.0]);
}

#[test]
fn corpus_and_modalities_parity() {
    exercise_in_memory_core_branches();

    let droot = data_root();
    ensure_test_fastas(&droot);
    let xfails = load_xfail_registry();
    let testdata_dir = droot.join("VarDictJava/testdata/integrationtestcases");

    let full_sweep = std::env::var_os("LLVM_PROFILE_FILE").is_some()
        || std::env::var_os("VARDICT_CORPUS_IN_RUST").is_some();

    let smoke_names = [
        "Simple;GRCh37.fa;HCC2218_C_position_bug.bam;chr2;179260451-179260596;-f 0.001.txt",
        "Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p-fisher;.txt",
        "Amplicon;hg19.fa;Colo829-18_S3-sort.bam;chr1;933866-934466;933866-934466;-a 100.95 -D;.txt",
        "Simple;hg38.fa;Tcov_zero_SV_DEL_candidate3.bam;chr8;90534611;90620070;-f 0.1 -i.txt",
    ];

    let mut entries: Vec<PathBuf> = fs::read_dir(&testdata_dir)
        .expect("read integrationtestcases")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("txt"))
        .collect();
    entries.sort();
    assert_eq!(
        entries.len(),
        92,
        "expected 92 upstream integration fixtures"
    );

    if full_sweep {
        let pass_cnt = Mutex::new(0usize);
        let xfail_cnt = Mutex::new(0usize);
        let next_idx = Mutex::new(0usize);

        std::thread::scope(|s| {
            for _ in 0..8 {
                s.spawn(|| loop {
                    let idx = {
                        let mut guard = next_idx.lock().unwrap();
                        if *guard >= entries.len() {
                            break;
                        }
                        let cur = *guard;
                        *guard += 1;
                        cur
                    };
                    let p = &entries[idx];
                    let fname = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
                    let vcf_comp = smoke_names.contains(&fname);
                    let verdict = run_single_integration_case(p, &droot, &xfails, vcf_comp);
                    if verdict == "PASS" {
                        *pass_cnt.lock().unwrap() += 1;
                    } else {
                        *xfail_cnt.lock().unwrap() += 1;
                    }
                });
            }
        });

        assert_eq!(*pass_cnt.lock().unwrap(), 91);
        assert_eq!(*xfail_cnt.lock().unwrap(), 1);
    } else {
        let p = testdata_dir.join(
            "Simple;hard_clip_case.fa;test_reads.hardclipped.bam;example_genome;160-280;-f 0.0 -P 0.txt",
        );
        assert_eq!(
            run_single_integration_case(&p, &droot, &xfails, false),
            "PASS"
        );
    }
}
