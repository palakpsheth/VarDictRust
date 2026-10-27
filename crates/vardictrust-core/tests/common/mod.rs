//! Shared pure-Rust fixture builders and in-process runner helpers for
//! `vardictrust-core` integration tests (`corpus_parity.rs` and `engine_behavior.rs`).
//! Zero `samtools`, `vardictrust`, or `sha256sum` subprocesses.

#![allow(dead_code)]

use std::{
    collections::HashMap,
    fs::File,
    io,
    path::{Path, PathBuf},
    process::Command,
    sync::Once,
};

use noodles::bam::bai;
use noodles::core::Position;
use noodles::csi::binning_index::index::reference_sequence::{bin::Chunk, index::LinearIndex};
use noodles::csi::binning_index::Indexer;
use noodles::sam::{
    alignment::{
        io::Write as _,
        record::{
            cigar::{op::Kind, Op},
            data::field::Tag,
        },
        record_buf::{Cigar, RecordBuf},
    },
    header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap},
    Header,
};
use tempfile::TempDir;
use vardictrust_core::{config::Config, types::VardictError};

pub fn workspace_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
}

pub fn data_root() -> PathBuf {
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

pub fn ensure_fastas_converted() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let droot = data_root();
        let fasta_dir = droot.join("VarDictJava/testdata/fastas");
        if !fasta_dir.join("GRCh37.fa").is_file()
            || !fasta_dir.join("GRCh37.fa.fai").is_file()
            || !fasta_dir.join("hg19.fa").is_file()
            || !fasta_dir.join("hg38.fa").is_file()
        {
            let script = workspace_root().join("harness/convert_fastas.py");
            let status = Command::new("python3")
                .arg(&script)
                .current_dir(&droot)
                .status()
                .expect("run convert_fastas.py");
            assert!(status.success(), "convert_fastas.py failed");
        }
    });
}

/// Write a BAM plus its `.bai` in pure Rust.
pub fn indexed_bam(path: &Path, chr: &str, len: usize, recs: &[RecordBuf]) -> io::Result<()> {
    let rs = HeaderMap::<RefSeqMap>::builder()
        .set_length(len.try_into().unwrap())
        .build()
        .unwrap();
    let hd = Header::builder().add_reference_sequence(chr, rs).build();
    let bgzf = noodles::bgzf::io::Writer::new(File::create(path)?);
    let mut w = noodles::bam::io::Writer::from(bgzf);
    w.write_alignment_header(&hd)?;
    let mut indexer = Indexer::<LinearIndex>::new(14, 5).expect("valid index parameters");
    for rec in recs {
        let start = w.get_ref().virtual_position();
        w.write_alignment_record(&hd, rec)?;
        let end = w.get_ref().virtual_position();
        let s = rec.alignment_start().expect("aligned record");
        let e = Position::new(s.get() + rec.cigar().alignment_span().saturating_sub(1).max(1))
            .unwrap_or(Position::MIN);
        indexer
            .add_record(Some((0, s, e, true)), Chunk::new(start, end))
            .map_err(io::Error::other)?;
    }
    w.try_finish()?;
    let mut bai_path = PathBuf::from(path);
    bai_path.set_extension("bam.bai");
    bai::io::Writer::new(File::create(bai_path)?).write_index(&indexer.build(1))
}

/// Encode an existing BAM into an indexed CRAM 3.0 + `.cram.crai` in pure Rust
/// against `ref_fa` (optionally stripping `NM`/`nM`/`MD` tags so the decoder must
/// reconstruct `NM` from the reference).
pub fn write_indexed_cram(
    bam_path: &Path,
    ref_fa: &Path,
    cram_path: &Path,
    strip_nm: bool,
    only_contig: Option<&str>,
) -> io::Result<()> {
    let fai_path = PathBuf::from(format!("{}.fai", ref_fa.display()));
    let fai = noodles::fasta::fai::fs::read(&fai_path)?;
    let mut fai_lengths: HashMap<String, usize> = HashMap::new();
    for rec in fai.as_ref() {
        let name = String::from_utf8_lossy(rec.name()).into_owned();
        fai_lengths.insert(name, rec.length() as usize);
    }

    let make_repo = || -> io::Result<noodles::fasta::Repository> {
        let fasta_reader = noodles::fasta::io::IndexedReader::new(
            io::BufReader::new(File::open(ref_fa)?),
            fai.clone(),
        );
        Ok(noodles::fasta::Repository::new(
            noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
        ))
    };

    let inspect_repo = make_repo()?;
    let mut bam_reader = noodles::bam::io::reader::Builder.build_from_path(bam_path)?;
    let bam_header = bam_reader.read_header()?;

    // Keep only contigs present in ref_fa.fai (and matching only_contig if specified),
    // setting their lengths to match .fai.
    let mut cram_header_builder = Header::builder();
    let mut old_to_new_id: HashMap<usize, usize> = HashMap::new();
    let mut new_contig_seqs: Vec<Vec<u8>> = Vec::new();

    for (old_id, (name_bytes, _)) in bam_header.reference_sequences().iter().enumerate() {
        let name_str = String::from_utf8_lossy(name_bytes.as_slice()).into_owned();
        if let Some(target_chr) = only_contig {
            if name_str != target_chr {
                continue;
            }
        }
        if let Some(&fai_len) = fai_lengths.get(&name_str) {
            let new_id = new_contig_seqs.len();
            old_to_new_id.insert(old_id, new_id);
            let rs = HeaderMap::<RefSeqMap>::builder()
                .set_length(fai_len.try_into().unwrap())
                .build()
                .unwrap();
            cram_header_builder = cram_header_builder.add_reference_sequence(name_str, rs);
            let seq = inspect_repo
                .get(name_bytes.as_slice())
                .transpose()?
                .map(|s| s.as_ref().as_ref().to_vec())
                .unwrap_or_default();
            new_contig_seqs.push(seq);
        }
    }
    let cram_header = cram_header_builder.build();

    let repo = make_repo()?;
    let mut cram_writer = noodles::cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_path(cram_path)?;
    cram_writer.write_alignment_header(&cram_header)?;

    let rg_tag = Tag::READ_GROUP;
    let nm_tag = Tag::EDIT_DISTANCE;
    let nm_alt_tag = Tag::new(b'n', b'M');
    let md_tag = Tag::new(b'M', b'D');

    for rec_res in bam_reader.record_bufs(&bam_header) {
        let mut rec = rec_res?;
        let Some(old_ref_id) = rec.reference_sequence_id() else {
            continue;
        };
        let Some(&new_ref_id) = old_to_new_id.get(&old_ref_id) else {
            continue;
        };
        *rec.reference_sequence_id_mut() = Some(new_ref_id);

        if let Some(old_mate_id) = rec.mate_reference_sequence_id() {
            if let Some(&new_mate_id) = old_to_new_id.get(&old_mate_id) {
                *rec.mate_reference_sequence_id_mut() = Some(new_mate_id);
            } else {
                *rec.mate_reference_sequence_id_mut() = None;
                *rec.mate_alignment_start_mut() = None;
            }
        }

        rec.data_mut().remove(&rg_tag);
        if strip_nm {
            rec.data_mut().remove(&nm_tag);
            rec.data_mut().remove(&nm_alt_tag);
            rec.data_mut().remove(&md_tag);
        }

        // Align read base casing with soft-masked reference bases where the base
        // matches case-insensitively so CRAM encodes reference-matches without
        // substitution features (reproducing samtools view -C external-ref behavior).
        if let Some(start_pos) = rec.alignment_start() {
            let ref_seq = &new_contig_seqs[new_ref_id];
            let mut ref_idx = start_pos.get().saturating_sub(1);
            let mut read_idx = 0usize;
            let ops: Vec<Op> = rec.cigar().as_ref().to_vec();
            let seq_bytes = rec.sequence_mut().as_mut();
            for op in ops {
                let len = op.len();
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        for _ in 0..len {
                            if let (Some(rb), Some(&fb)) =
                                (seq_bytes.get_mut(read_idx), ref_seq.get(ref_idx))
                            {
                                if rb.eq_ignore_ascii_case(&fb) {
                                    *rb = fb;
                                }
                            }
                            read_idx += 1;
                            ref_idx += 1;
                        }
                    }
                    Kind::Insertion | Kind::SoftClip => {
                        read_idx += len;
                    }
                    Kind::Deletion | Kind::Skip => {
                        ref_idx += len;
                    }
                    Kind::HardClip | Kind::Pad => {}
                }
            }
        }

        cram_writer.write_alignment_record(&cram_header, &rec)?;
    }

    cram_writer.try_finish(&cram_header)?;
    let crai = noodles::cram::fs::index(cram_path)?;
    noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai)
}

/// Single-line FASTA + hand-written `.fai` row; unique per call so a stale `.fai`
/// cannot mask the missing-index case.
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

/// Run the engine in-process in intermediate `.var` table mode (`threads = 1`).
pub fn run_core_var(config: &Config) -> Result<String, VardictError> {
    let mut cfg = config.clone();
    cfg.threads = 1;
    cfg.emit_intermediate = true;
    cfg.apply_pileup_overrides();
    let mut buf = Vec::with_capacity(4096);
    vardictrust_core::run_with_writer(&cfg, &mut buf)?;
    Ok(String::from_utf8(buf).expect("engine emitted valid UTF-8"))
}

/// Run the engine in-process in VCF mode (`threads = 1`).
pub fn run_core_vcf(config: &Config) -> Result<String, VardictError> {
    let mut cfg = config.clone();
    cfg.threads = 1;
    cfg.emit_intermediate = false;
    cfg.apply_pileup_overrides();
    let mut buf = Vec::with_capacity(4096);
    vardictrust_core::run_with_writer(&cfg, &mut buf)?;
    Ok(String::from_utf8(buf).expect("engine emitted valid UTF-8"))
}

/// Extract non-header, non-empty VCF variant lines.
pub fn extract_variants(vcf_text: &str) -> Vec<String> {
    vcf_text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|s| s.to_string())
        .collect()
}

/// Pure-Rust FIPS 180-4 SHA-256 hex digest (0 subprocesses).
pub fn sha256_hex(input: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];

    let bit_len = (input.len() as u64).wrapping_mul(8);
    let mut msg = input.to_vec();
    msg.push(0x80);
    while (msg.len() % 64) != 56 {
        msg.push(0x00);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in msg.chunks_exact(64) {
        let mut w = [0u32; 64];
        for (i, word) in chunk.chunks_exact(4).enumerate() {
            w[i] = u32::from_be_bytes([word[0], word[1], word[2], word[3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let temp1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);

            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }

    let mut out = String::with_capacity(64);
    for word in h {
        use std::fmt::Write as _;
        let _ = write!(out, "{word:08x}");
    }
    out
}
