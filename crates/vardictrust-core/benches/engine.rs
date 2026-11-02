//! Engine hot-path benches (T147-lite + Issue N): pileup/parse_cigar, find_bp,
//! find_35_match, javamap/compute_if_absent, round_half_even, modify_cigar, and fisher.
//! Fixtures are built in-code — no committed external genomic data required.

use std::hint::black_box;
use std::sync::Arc;

use criterion::{criterion_group, criterion_main, Criterion};
use vardictrust_core::bam::SamRecord;
use vardictrust_core::cigar::modify::CigarModifier;
use vardictrust_core::cigar::parse::pileup;
use vardictrust_core::config::Config;
use vardictrust_core::javamap::JavaHashMap;
use vardictrust_core::realign::{find_35_match, find_bp};
use vardictrust_core::reference::{LoadedRegion, Reference, SeedMap};
use vardictrust_core::regions::Region;
use vardictrust_core::stats::fisher::FisherExact;
use vardictrust_core::types::{round_half_even, Pos1};

const CHR: &str = "chr1";
const REF_LEN: usize = 20_000;

/// Order-3 de Bruijn 64-bp motif over {A, C, G, T}: every 3-mer occurs once, avoiding period-4 artifacts.
const MOTIF_64: &[u8; 64] = b"AAACAGAATAGACCACGACTAGCAGGAGTATCATGATTCCGCCTCGGCGTCTGCTTGGGTGTTT";

fn make_reference() -> Reference {
    let seq: Vec<u8> = (0..REF_LEN).map(|i| MOTIF_64[i % 64]).collect();
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
        // Bases follow the non-periodic 64-bp reference motif so most positions match.
        let seq: Vec<u8> = (0..rl).map(|j| MOTIF_64[(start + j) % 64]).collect();
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
    let config = Arc::new(Config::default());
    c.bench_function("parse_cigar/100_reads_mixed_cigars", |b| {
        b.iter(|| {
            black_box(pileup(
                black_box(&records),
                &reference,
                &region,
                &config,
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
    // Synthetic 60bp consensus with ~10% mismatches vs the non-periodic 64-bp motif.
    let consensus: Vec<u8> = (0..60)
        .map(|i| if i % 10 == 7 { b'T' } else { MOTIF_64[i % 64] })
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

/// find_35_match scans 5' and 3' soft-clip consensus sequences for reverse-complement
/// overlap during large insertion realignment (realignlgins30).
fn bench_find_35_match(c: &mut Criterion) {
    let seq5: Vec<u8> = (0..128).map(|i| MOTIF_64[(i * 7 + 3) % 64]).collect();
    // Construct seq3 whose reverse-complement Suffix overlaps seq5 prefix over 24 bp with 1 mismatch.
    let mut seq3: Vec<u8> = (0..104).map(|i| MOTIF_64[(i * 11 + 19) % 64]).collect();
    for i in (0..24).rev() {
        let b = seq5[i];
        let comp = match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        };
        seq3.push(if i == 11 { b'A' } else { comp });
    }
    c.bench_function("find_35_match", |b| {
        b.iter(|| black_box(find_35_match(black_box(&seq5), black_box(&seq3), false)))
    });
}

/// JavaHashMap::compute_if_absent and get_or_put_with hot loops across pileup positions with ~10:1 hit-to-miss ratio.
fn bench_javamap_compute_if_absent(c: &mut Criterion) {
    let keys: Vec<i32> = (0..2_000).map(|i| 10_000 + ((i * 37) % 200)).collect();
    c.bench_function("javamap/compute_if_absent", |b| {
        b.iter(|| {
            let mut map: JavaHashMap<i32, u32> = JavaHashMap::new();
            for &k in black_box(&keys) {
                let entry = map.compute_if_absent(k, |_| 0);
                *entry = entry.wrapping_add(1);
            }
            black_box(map.len())
        })
    });
    c.bench_function("javamap/get_or_put", |b| {
        b.iter(|| {
            let mut map: JavaHashMap<i32, u32> = JavaHashMap::new();
            for &k in black_box(&keys) {
                let entry = map.get_or_put_with(k, || 0);
                *entry = entry.wrapping_add(1);
            }
            black_box(map.len())
        })
    });
}

/// Banker's rounding (round_half_even) over realistic AF, quality, and Fisher logdc inputs.
fn bench_round_half_even(c: &mut Criterion) {
    let samples: [(f64, i32); 12] = [
        (0.12345, 4),
        (1.0 / 20_000.0, 4),
        (1.0 / 2_000.0, 3),
        (27.85, 1),
        (37.25, 1),
        (0.015625, 5),
        (0.190625, 5),
        (2.96, 1),
        (0.99999, 4),
        (-3.45678915, 7),
        (14.125, 2),
        (5.0e-9, 8),
    ];
    c.bench_function("round_half_even", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(val, digits) in black_box(&samples) {
                acc += round_half_even(val, digits);
            }
            black_box(acc)
        })
    });
}

/// CigarModifier::modify_cigar over a realistic mix of matched, soft-clipped, and multi-indel reads.
fn bench_modify_cigar(c: &mut Criterion) {
    let config = Arc::new(Config::default());
    let ref_seq: Vec<u8> = (0..2_000).map(|i| MOTIF_64[i % 64]).collect();
    let seed = SeedMap::new();
    let cases: [(&str, usize, usize); 5] = [
        ("150M", 200, 150),
        ("12S138M", 220, 150),
        ("60M2D4M1D86M", 250, 150),
        ("15S50M3I82M", 300, 150),
        ("8M1I141M", 350, 150),
    ];
    let reads: Vec<(i32, String, Vec<u8>, Vec<u8>)> = cases
        .iter()
        .map(|&(cig, pos, len)| {
            let seq: Vec<u8> = (0..len).map(|j| ref_seq[pos + j]).collect();
            let qual: Vec<u8> = (0..len).map(|j| b'!' + ((30 + (j % 10)) as u8)).collect();
            (pos as i32, cig.to_string(), seq, qual)
        })
        .collect();

    c.bench_function("modify_cigar", |b| {
        b.iter(|| {
            let mut total_pos = 0i32;
            for (pos, cig, seq, qual) in black_box(&reads) {
                let modifier = CigarModifier::new(
                    *pos,
                    cig.clone(),
                    seq.clone(),
                    qual.clone(),
                    &ref_seq,
                    1,
                    &seed,
                    10,
                    150,
                    &config,
                );
                if let Ok(res) = modifier.modify_cigar() {
                    total_pos = total_pos.wrapping_add(res.position);
                }
            }
            black_box(total_pos)
        })
    });
}

/// FisherExact contingency tables across R-mode and Java (--fisher) mode.
fn bench_fisher(c: &mut Criterion) {
    let tables: [(i32, i32, i32, i32); 6] = [
        (27, 17, 3, 8),
        (9, 6, 5, 1),
        (10, 8, 1, 3),
        (86, 41, 0, 1),
        (60, 259, 1, 0),
        (145, 494, 12, 19),
    ];
    c.bench_function("fisher", |b| {
        b.iter(|| {
            let mut p_sum = 0.0;
            for &(rf, rr, af, ar) in black_box(&tables) {
                let r_mode = FisherExact::new(rf, rr, af, ar);
                let j_mode = FisherExact::new_java(rf, rr, af, ar);
                p_sum += r_mode.get_p_value() + j_mode.get_p_value();
                black_box(r_mode.get_odds_ratio());
                black_box(j_mode.get_odds_ratio());
            }
            black_box(p_sum)
        })
    });
}

criterion_group!(
    benches,
    bench_parse_cigar,
    bench_find_bp,
    bench_find_35_match,
    bench_javamap_compute_if_absent,
    bench_round_half_even,
    bench_modify_cigar,
    bench_fisher,
);
criterion_main!(benches);
