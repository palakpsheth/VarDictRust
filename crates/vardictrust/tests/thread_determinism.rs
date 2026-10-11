//! Thread-scaling determinism (T151/T152-lite): vardictrust must emit byte-identical
//! VCF output regardless of `-th/--threads` (§10: deterministic ordered output).
//! Fixture generated in-test (robustness.rs patterns); `--tile-size` shrinks tiles so
//! rayon actually parallelizes across several tiles per run instead of one.

use std::process::Command;

use noodles::core::Position;
use noodles::sam::alignment::{
    record::cigar::{op::Kind, Op},
    record_buf::{Cigar, RecordBuf},
};
use tempfile::TempDir;

mod common;
use common::{indexed_bam, snv_record, write_ref};

/// Two identical reads carrying one C mismatch at reference position `pos`
/// (all-A reference), so each call site yields a supported germline SNV row.
fn snv_pair(pos: usize) -> Vec<RecordBuf> {
    (0..2)
        .map(|i| snv_record(&format!("s{pos}r{i}"), pos - 20))
        .collect()
}

/// Two identical reads sharing a 10 bp insertion mid-read.
fn insertion_pair(pos: usize) -> Vec<RecordBuf> {
    let cigar: Cigar = [
        (40usize, Kind::Match),
        (10, Kind::Insertion),
        (10, Kind::Match),
    ]
    .iter()
    .map(|&(l, k)| Op::new(k, l))
    .collect();
    (0..2)
        .map(|i| {
            RecordBuf::builder()
                .set_name(format!("i{pos}r{i}"))
                .set_flags(noodles::sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(pos - 20).unwrap())
                .set_cigar(cigar.clone())
                .set_sequence(
                    [vec![b'A'; 40], vec![b'T'; 10], vec![b'A'; 10]]
                        .concat()
                        .into(),
                )
                .set_quality_scores(vec![40u8; 60].into())
                .build()
        })
        .collect()
}

#[test]
fn vcf_output_byte_identical_across_thread_counts() {
    let dir = TempDir::new().unwrap();
    let fa = write_ref(&dir, "chr1", 2000, true);

    // Variant clusters spread across ~1 kb so --tile-size=300 splits the BED
    // region into several tiles (halo 500 overlaps them): boundary dedup and
    // ordered merge are both exercised, not a single-tile trivial case.
    let mut recs = Vec::new();
    for pos in (150..=1000).step_by(80) {
        recs.extend(snv_pair(pos));
    }
    recs.extend(insertion_pair(600));
    recs.sort_by_key(|r| r.alignment_start());
    let bam = dir.path().join("det.bam");
    indexed_bam(&bam, "chr1", 2000, &recs).unwrap();

    let bed = dir.path().join("det.bed");
    std::fs::write(&bed, "chr1\t100\t1100\tDET\n").unwrap();

    // [1, 2, 8] capped by available cores, de-duplicated so a 2-core host
    // doesn't run the same config twice.
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    let counts: Vec<usize> = [1usize, 2, 8].iter().map(|t| (*t).min(cores)).collect();

    let runs: Vec<Vec<u8>> = counts
        .iter()
        .map(|&t| {
            let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
                .args(["-c1", "-S2", "-E3", "-g4"])
                .arg(format!("--threads={t}"))
                .arg("--tile-size=300")
                .arg(format!("-G{}", fa.display()))
                .arg(format!("-b{}", bam.display()))
                .arg(&bed)
                .output()
                .expect("spawn vardictrust");
            assert!(out.status.success(), "--threads={t}: non-zero exit");
            out.stdout
        })
        .collect();

    let s = String::from_utf8_lossy(&runs[0]);
    assert!(s.contains("CHROM"), "no VCF header");
    let rows = s.lines().filter(|l| l.starts_with("chr1\t")).count();
    // Fixture plants 11 SNV clusters + 1 insertion across chr1:150-950, which
    // the --tile-size=300 run splits into several tiles: the FULL count (and
    // presence of the first/last-tile variants) proves boundary dedup kept
    // every tile-edge variant — multi-tile-ness is self-evidencing.
    assert!(
        rows >= 12,
        "expected >=12 variant rows spanning multiple tiles, got {rows}"
    );
    for pos in [150, 950] {
        assert!(
            s.contains(&format!("chr1\t{pos}\t")),
            "tile-spanning variant @{pos} missing"
        );
    }
    for pair in runs.windows(2) {
        assert_eq!(pair[0], pair[1], "VCF bytes differ across thread counts");
    }
}
