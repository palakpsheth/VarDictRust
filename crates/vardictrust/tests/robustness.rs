//! Robustness corpus (T158 / Appendix C T158, T159, T09).
//!
//! A variant caller in a clinical pipeline must never panic on adversarial input
//! (BAM with truncated headers, malformed CIGARs, missing FAI, empty BAMs). It may emit
//! graceful non-zero with a stderr message; death by signal (SIGSEGV/SIGABRT) fails.

use std::{
    path::Path,
    process::{Command, Output},
};

use noodles::core::Position;
use noodles::sam::alignment::{
    record::cigar::{op::Kind, Op},
    record_buf::{Cigar, RecordBuf},
};
use tempfile::TempDir;

mod common;
use common::{indexed_bam, snv_record, write_ref};

/// CIGAR shapes that historically crashed Java's modifyCigar (bounds family
/// #341/#243/#253/#283/#7/#38; seeds from cigar/modify.rs crash-corpus unit
/// test + leading/trailing-D and H-before-S). Pure-D / "0M" seeds have zero read
/// length and cannot be encoded in a spec-valid BAM ('*'-SEQ records are dropped
/// by preprocessing), so deletions span past the reference ends instead.
#[rustfmt::skip]
#[allow(clippy::type_complexity)]
const BAD_CIGARS: &[(&str, usize, &[(usize, Kind)])] = &[
    ("del-at-chr-start",   1,   &[(10, Kind::Match), (200, Kind::Deletion), (10, Kind::Match)]),
    ("leading-D",          50,  &[(2, Kind::Deletion), (10, Kind::Match)]),
    ("trailing-D",         70,  &[(10, Kind::Match), (3, Kind::Deletion)]),
    ("giant-ins",          90,  &[(5, Kind::Match), (100, Kind::Insertion), (5, Kind::Match)]),
    ("mixed-mid",          120, &[(1, Kind::Match), (1, Kind::Insertion), (1, Kind::Deletion), (1, Kind::Match)]),
    ("h-before-s",         140, &[(5, Kind::HardClip), (5, Kind::SoftClip), (10, Kind::Match)]),
    ("del-past-ref-end",   985, &[(10, Kind::Match), (500, Kind::Deletion), (10, Kind::Match)]),
];

#[rustfmt::skip]
fn read(pos: usize, ops: &[(usize, Kind)]) -> RecordBuf {
    // SEQ length must equal the CIGAR read length (BAM spec 1.4.10); zero-read-length
    // CIGARs ("0M"/pure-D) are therefore unspecifiable and asserted out here.
    let cigar: Cigar = ops.iter().map(|&(l, k)| Op::new(k, l)).collect();
    let rl = cigar.read_length();
    assert!(rl > 0, "zero-read-length CIGARs are unspecifiable");
    RecordBuf::builder()
        .set_name(format!("r{pos}"))
        .set_flags(noodles::sam::alignment::record::Flags::empty())
        .set_reference_sequence_id(0).set_alignment_start(Position::try_from(pos).unwrap())
        .set_cigar(cigar).set_sequence(vec![b'A'; rl].into()).set_quality_scores(vec![40; rl].into())
        .build()
}

fn run_cli(dir: &TempDir, bed_lines: &str, fasta: &Path, bam: &Path) -> Output {
    let bed = dir.path().join("case.bed");
    std::fs::write(&bed, bed_lines).unwrap();
    Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .args(["-c1", "-S2", "-E3", "-g4", "--threads=1"])
        .arg(format!("-G{}", fasta.display()))
        .arg(format!("-b{}", bam.display()))
        .arg(bed)
        .output()
        .expect("spawn vardictrust")
}

/// No signal death, no Rust-panic leak (exit 101 / "panicked at"); failure only
/// ever as a graceful error message on stderr.
#[rustfmt::skip]
fn assert_graceful(out: &Output, case: &str) {
    assert!(out.status.code().is_some(), "{case}: killed by signal");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.code() != Some(101) && !stderr.contains("panicked at"),
        "{case}: died by panic (exit 101 / 'panicked at' on stderr), not gracefully"
    );
    if !out.status.success() {
        assert!(!out.stderr.is_empty(), "{case}: failed with no stderr message");
    }
}

/// A successful run must print the VCF header — evidence the full pipeline ran.
#[rustfmt::skip]
fn assert_ran(out: &Output, case: &str) {
    assert!(out.status.success() && out.stdout.windows(5).any(|w| w == b"CHROM"),
            "{case}: expected clean exit with VCF header output");
}

#[test]
fn robustness_corpus_never_crashes() {
    let dir = TempDir::new().unwrap();

    // Case 1: empty BAM region query (indexed header-only BAM).
    let empty_bam = dir.path().join("empty.bam");
    indexed_bam(&empty_bam, "chr1", 1000, &[]).unwrap();
    let fa1 = write_ref(&dir, "chr1", 1000, true);
    let out = run_cli(&dir, "chr1\t100\t200\tG1\n", &fa1, &empty_bam);
    assert_graceful(&out, "empty-bam");
    assert_ran(&out, "empty-bam");

    // Case 2: malformed-CIGAR family through the full pipeline.
    let bad_bam = dir.path().join("bad.bam");
    let recs: Vec<_> = BAD_CIGARS
        .iter()
        .map(|(_, pos, ops)| read(*pos, ops))
        .collect();
    indexed_bam(&bad_bam, "chr1", 1000, &recs).unwrap();
    // Region end 999 (not ~300): the del-past-ref-end seed sits at pos 985, so a
    // narrower query drops it pre-parse and the boundary-spanning shape is never
    // exercised.
    let out = run_cli(&dir, "chr1\t50\t999\tG2\n", &fa1, &bad_bam);
    assert_graceful(&out, "bad-cigar");
    assert_ran(&out, "bad-cigar");

    // Case 3: BED end<start and end beyond chromosome length.
    let out = run_cli(
        &dir,
        "chr1\t300\t100\tBAD\nchr1\t900\t5000\tOOB\n",
        &fa1,
        &empty_bam,
    );
    assert_graceful(&out, "bed-oob");

    // Case 4: FASTA missing its .fai index -> graceful loud error, not panic.
    let fa_no_idx = write_ref(&dir, "chr1", 1000, false);
    let out = run_cli(&dir, "chr1\t100\t200\tG4\n", &fa_no_idx, &empty_bam);
    assert_graceful(&out, "missing-fai");
    assert!(!out.status.success(), "missing-fai: expected a loud error");

    // Case 5: region chromosome in the BAM header but absent from the reference FASTA.
    let chrz_bam = dir.path().join("chrz.bam");
    indexed_bam(&chrz_bam, "chrZ", 500, &[read(50, &[(15, Kind::Match)])]).unwrap();
    let out = run_cli(&dir, "chrZ\t50\t200\tG5\n", &fa1, &chrz_bam);
    assert_graceful(&out, "absent-chr");
}

/// T159: float parsing/formatting must be locale-independent end-to-end. Rust
/// std ignores LC_* by design; this pins that assumption against future dep
/// drift (an FFI/C dependency honoring the locale would leak comma-decimal
/// output). Skips with a clear message when no de_DE locale is installed.
#[test]
fn locale_independence_c_vs_de_de() {
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

    let dir = TempDir::new().unwrap();
    let fa = write_ref(&dir, "chr1", 1000, true);
    // Two reads sharing one mismatch → a real variant row (AF/PMEAN/QUAL float
    // cells), so the comparison covers record formatting, not just the header.
    let bam = dir.path().join("loc.bam");
    indexed_bam(
        &bam,
        "chr1",
        1000,
        &[snv_record("r1", 100), snv_record("r2", 100)],
    )
    .unwrap();

    let run = |envs: &[(&str, &str)]| {
        let bed = dir.path().join("loc.bed");
        std::fs::write(&bed, "chr1\t100\t300\tLOC\n").unwrap();
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_vardictrust"));
        cmd.args(["-c1", "-S2", "-E3", "-g4", "--threads=1"])
            .arg(format!("-G{}", fa.display()))
            .arg(format!("-b{}", bam.display()))
            .arg(&bed);
        for (k, v) in envs {
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
        "output bytes differ between LC_ALL=C and de_DE (comma-decimal leakage?)"
    );
}

/// T09 / PR #402: 11 copies of an adversarial read with a homopolymer 3-deletion
/// cluster causing negative-length operator generation. VarDictJava aborts with
/// `Critical exception occurs` after 10 exceptions; VarDictRust must skip each
/// T09 / PR #402: 11 copies of an adversarial read with a homopolymer 3-deletion
/// cluster causing negative-length operator generation. VarDictJava aborts with
/// `Critical exception occurs` after 10 exceptions; VarDictRust must skip each
/// bad read without modifying pileup state, never abort, and emit identical output
/// between 1 and 11 bad reads.
#[test]
fn test_pr402_11_bad_reads_does_not_abort() {
    let samtools_available = Command::new("samtools")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !samtools_available {
        println!("SKIP: samtools not available on host");
        return;
    }

    let workspace_root = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../.."));
    let temp_dir = TempDir::new().unwrap();
    let scratch_dir = temp_dir.path();

    let pr402_bam = workspace_root.join("fixtures/pr402/pr402.bam");
    let out_sam = scratch_dir.join("pr402_11.sam");
    let out_bam = scratch_dir.join("pr402_11.bam");

    let script = format!(
        r#"
        samtools view -h "{0}" > "{1}"
        grep -w "bad_read" "{1}" > "{2}/pr402_bad.sam"
        for i in $(seq 1 10); do
            sed "s/bad_read/bad_read_$i/" "{2}/pr402_bad.sam" >> "{1}"
        done
        samtools view -b "{1}" | samtools sort -o "{3}"
        samtools index "{3}"
    "#,
        pr402_bam.display(),
        out_sam.display(),
        scratch_dir.display(),
        out_bam.display()
    );

    let status = Command::new("bash")
        .arg("-c")
        .arg(script)
        .status()
        .expect("failed to execute bash");
    assert!(status.success());

    let ref_fa = workspace_root.join("fixtures/pr402/ref.fa");
    let pr402_bed = workspace_root.join("fixtures/pr402/pr402.bed");

    let cmd1 = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .args([
            "-G",
            ref_fa.to_str().unwrap(),
            "-b",
            pr402_bam.to_str().unwrap(),
            "-N",
            "SAMPLE",
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            pr402_bed.to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute vardictrust");

    let cmd11 = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .args([
            "-G",
            ref_fa.to_str().unwrap(),
            "-b",
            out_bam.to_str().unwrap(),
            "-N",
            "SAMPLE",
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            pr402_bed.to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute vardictrust");

    assert!(cmd1.status.success());
    assert!(cmd11.status.success());

    let out1 = String::from_utf8_lossy(&cmd1.stdout);
    let out11 = String::from_utf8_lossy(&cmd11.stdout);

    assert_eq!(
        out1, out11,
        "Outputs differ between 1 bad read and 11 bad reads"
    );
}

/// QA row 12 / B3: invalid -F flag syntax must fail with loud, actionable error,
/// not silently fall back to 0x504.
#[test]
fn test_invalid_flag_syntax_fails_with_actionable_error() {
    let dir = TempDir::new().unwrap();
    let fa = write_ref(&dir, "chr1", 1000, true);
    let bam = dir.path().join("empty.bam");
    indexed_bam(&bam, "chr1", 1000, &[]).unwrap();

    let bed = dir.path().join("case.bed");
    std::fs::write(&bed, "chr1\t100\t200\tG1\n").unwrap();

    let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .args(["-G", fa.to_str().unwrap()])
        .args(["-b", bam.to_str().unwrap()])
        .args(["-F", "bogus"])
        .arg(&bed)
        .output()
        .expect("spawn vardictrust");

    assert!(!out.status.success());
    assert_ne!(out.status.code(), Some(101), "died with panic");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Invalid -F value"),
        "stderr should mention invalid -F value: {stderr}"
    );
}

/// QA row 13 / B4: --tile-size 0 must fail immediately with actionable error,
/// not hang in an infinite loop.
#[test]
fn test_zero_tile_size_fails_fast_without_hanging() {
    let dir = TempDir::new().unwrap();
    let fa = write_ref(&dir, "chr1", 1000, true);
    let bam = dir.path().join("empty.bam");
    indexed_bam(&bam, "chr1", 1000, &[]).unwrap();

    let bed = dir.path().join("case.bed");
    std::fs::write(&bed, "chr1\t100\t200\tG1\n").unwrap();

    let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .args(["-G", fa.to_str().unwrap()])
        .args(["-b", bam.to_str().unwrap()])
        .args(["--tile-size", "0"])
        .arg(&bed)
        .output()
        .expect("spawn vardictrust");

    assert!(!out.status.success());
    assert_ne!(out.status.code(), Some(101), "died with panic");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Invalid --tile-size"),
        "stderr should mention invalid --tile-size: {stderr}"
    );
}
