//! Engine hot-path benches (T147-lite): pileup/parse_cigar over a synthetic
//! 100-read set with mixed CIGAR ops, and find_bp realignment scanning over
//! synthetic consensus data. Fixtures are built in-code — no committed data.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use vardict_core::bam::SamRecord;
use vardict_core::cigar::parse::pileup;
use vardict_core::config::Config;
use vardict_core::realign::find_bp;
use vardict_core::reference::{LoadedRegion, Reference};
use vardict_core::regions::Region;
use vardict_core::types::Pos1;

const CHR: &str = "chr1";
const REF_LEN: usize = 20_000;

fn make_reference() -> Reference {
    let seq: Vec<u8> = (0..REF_LEN).map(|i| b"ACGT"[i % 4]).collect();
    let mut r = Reference::new();
    r.loaded_regions.push(LoadedRegion {
        chr: CHR.into(),
        sequence_start: 1,
        sequence_end: REF_LEN as u32,
    });
    r.reference_sequences.insert(CHR.into(), seq);
    r
}

/// 100 coordinate-sorted reads, mixed CIGAR ops (M/D/I/S), 150bp each.
#[rustfmt::skip]
fn make_records() -> Vec<SamRecord> {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let mut recs = Vec::new();
    for i in 0..100 {
        let start = 200 + i * 180;
        let (cigar, rl): (Vec<(u32, Kind)>, usize) = match i % 5 {
            0 => (vec![(150, Kind::Match)], 150),
            1 => (vec![(10, Kind::SoftClip), (140, Kind::Match)], 150),
            2 => (vec![(50, Kind::Match), (6, Kind::Deletion), (94, Kind::Match)], 144),
            3 => (vec![(70, Kind::Match), (8, Kind::Insertion), (72, Kind::Match)], 150),
            _ => (vec![(20, Kind::Match), (4, Kind::Deletion), (126, Kind::Match)], 146),
        };
        // Bases follow the reference pattern so most positions match.
        let seq: Vec<u8> = (0..rl).map(|j| b"ACGT"[(start + j) % 4]).collect();
        recs.push(SamRecord {
            name: format!("r{i}").into_bytes(),
            flags: 0,
            reference_name: Some(CHR.into()),
            alignment_start: Some(Pos1(start as u32)),
            alignment_end: None,
            mapping_quality: Some(60),
            cigar,
            sequence: seq,
            quality_scores: vec![b'I'; rl],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: None,
            sa: None,
            mc: None,
            mq: None,
        });
    }
    recs
}

fn bench_parse_cigar(c: &mut Criterion) {
    let reference = make_reference();
    let region = Region {
        chr: CHR.into(),
        start: Pos1(1),
        end: Pos1(REF_LEN as u32),
        gene: "BENCH".into(),
        amplicon: None,
    };
    let records = make_records();
    c.bench_function("parse_cigar/100_reads_mixed_cigars", |b| {
        b.iter(|| {
            black_box(pileup(
                black_box(&records),
                &reference,
                &region,
                &Config::default(),
                false,
                records.len() as i32,
                0,
                0,
            ))
            .expect("pileup")
        })
    });
}

/// find_bp scans an indel-flanking consensus window against the reference
/// in one direction — the inner loop of realign_indels breakpoint search.
fn bench_find_bp(c: &mut Criterion) {
    let reference = make_reference();
    // Synthetic 60bp consensus with ~10% mismatches vs the repeating pattern.
    let consensus: Vec<u8> = (0..60)
        .map(|i| if i % 10 == 7 { b'T' } else { b"ACGT"[i % 4] })
        .collect();
    c.bench_function("find_bp/60bp_consensus_indel40", |b| {
        b.iter(|| {
            black_box(find_bp(
                black_box(&consensus),
                5000,
                &reference,
                -1,
                CHR,
                40,
                REF_LEN as u32,
            ))
        })
    });
}

criterion_group!(benches, bench_parse_cigar, bench_find_bp);
criterion_main!(benches);
