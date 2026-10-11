//! Shared fixture builders for the vardictrust integration tests (generated
//! in-test: no committed data); consumed by robustness.rs and
//! thread_determinism.rs.

use std::{
    io,
    path::{Path, PathBuf},
};

use noodles::bam::bai;
use noodles::core::Position;
use noodles::csi::binning_index::index::reference_sequence::{bin::Chunk, index::LinearIndex};
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
use tempfile::TempDir;

/// Write a BAM plus its `.bai`.
pub fn indexed_bam(path: &Path, chr: &str, len: usize, recs: &[RecordBuf]) -> io::Result<()> {
    let rs = HeaderMap::<RefSeqMap>::builder()
        .set_length(len.try_into().unwrap())
        .build()
        .unwrap();
    let hd = Header::builder().add_reference_sequence(chr, rs).build();
    let bgzf = noodles::bgzf::Writer::new(std::fs::File::create(path)?);
    let mut w = noodles::bam::io::Writer::from(bgzf);
    w.write_alignment_header(&hd)?;
    let mut indexer = Indexer::<LinearIndex>::new(14, 5);
    for rec in recs {
        let start = w.get_ref().virtual_position();
        w.write_alignment_record(&hd, rec)?;
        let end = w.get_ref().virtual_position();
        let s = rec.alignment_start().expect("aligned record");
        let e = Position::new(s.get() + rec.cigar().alignment_span()).unwrap_or(Position::MIN);
        indexer.add_record(Some((0, s, e, true)), Chunk::new(start, end))?;
    }
    w.try_finish()?;
    let mut bai_path = PathBuf::from(path);
    bai_path.set_extension("bam.bai"); // "x.bam" -> "x.bam.bai"
    bai::io::Writer::new(std::fs::File::create(bai_path)?).write_index(&indexer.build(1))
}

/// Single-line FASTA + hand-written .fai row; unique per call so a stale .fai
/// can't mask the missing-index case.
#[rustfmt::skip]
pub fn write_ref(dir: &TempDir, name: &str, len: usize, with_fai: bool) -> PathBuf {
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let fa = dir.path().join(format!("ref{n}.fa"));
    std::fs::write(&fa, format!(">{name}\n{}\n", "A".repeat(len))).unwrap();
    if with_fai { std::fs::write(fa.with_extension("fa.fai"), format!("{name}\t{len}\t{}\t{len}\t{}\n", name.len() + 2, len + 1)).unwrap(); }
    fa
}

/// One 50M read starting at reference position `start` (all-A except a single
/// C 20 bp in), so two such reads yield a supported germline SNV row against
/// an all-A reference.
pub fn snv_record(name: &str, start: usize) -> RecordBuf {
    let cigar: Cigar = [(50usize, Kind::Match)]
        .iter()
        .map(|&(l, k)| Op::new(k, l))
        .collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(noodles::sam::alignment::record::Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(start).unwrap())
        .set_cigar(cigar)
        .set_sequence([vec![b'A'; 20], vec![b'C'], vec![b'A'; 29]].concat().into())
        .set_quality_scores(vec![40u8; 50].into())
        .build()
}
