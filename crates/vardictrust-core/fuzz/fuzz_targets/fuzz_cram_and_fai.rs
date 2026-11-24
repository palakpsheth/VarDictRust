#![no_main]

use std::fs::{self, File};
use std::io::Write as _;
use std::path::Path;
use std::sync::OnceLock;

use libfuzzer_sys::fuzz_target;
use vardictrust_core::bam::{clear_bam_caches, query_region_iter, RecordPreprocessor, SamRecord};
use vardictrust_core::config::Config;
use vardictrust_core::reference::fetch_subsequence;
use vardictrust_core::run_with_writer;
use vardictrust_core::types::{Pos0, Pos1};

/// Maximum accepted libFuzzer input payload (`8 KiB`) to bound per-iteration I/O and runtime.
const MAX_FUZZ_INPUT_BYTES: usize = 8192;
/// Maximum queried or tiled genomic interval span (`4096 bp`) per iteration.
const MAX_REGION_SPAN_BP: u32 = 4096;
/// Maximum records drained per `query_region_iter` probe to prevent mutated index loops.
const MAX_DRAINED_RECORDS: usize = 64;

struct FuzzTemplates {
    plain_fa: Vec<u8>,
    plain_fai: Vec<u8>,
    bgzf_fa: Vec<u8>,
    bgzf_gzi: Vec<u8>,
    valid_bam: Vec<u8>,
    valid_bai: Vec<u8>,
    valid_cram: Vec<u8>,
    valid_crai: Vec<u8>,
}

static TEMPLATES: OnceLock<Option<FuzzTemplates>> = OnceLock::new();

struct CacheClearGuard;

impl Drop for CacheClearGuard {
    fn drop(&mut self) {
        clear_bam_caches();
    }
}

fn create_fuzz_tempdir() -> Option<tempfile::TempDir> {
    let shm_root = std::env::var_os("VDR_FUZZ_SHM_DIR")
        .filter(|v| !v.is_empty())
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::path::PathBuf::from("/dev/shm"));
    tempfile::Builder::new()
        .prefix("vdr-fuzz-")
        .tempdir_in(&shm_root)
        .or_else(|_| tempfile::tempdir())
        .ok()
}

fn init_templates() -> Option<FuzzTemplates> {
    use noodles::csi::binning_index::index::reference_sequence::{bin::Chunk, index::LinearIndex};
    use noodles::csi::binning_index::Indexer;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record::cigar::{op::Kind, Op};
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record::{Flags, MappingQuality};
    use noodles::sam::alignment::record_buf::{
        data::field::Value as BufValue, Cigar, Data, RecordBuf,
    };
    use noodles::sam::header::record::value::{
        map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
    };
    use noodles::sam::Header;

    let temp_dir = create_fuzz_tempdir()?;
    let dir = temp_dir.path();

    // Build a deterministic 512-bp reference sequence (8 lines x 64 bp) on `chr1`.
    let ref_len = 512usize;
    let mut ref_seq = Vec::with_capacity(ref_len);
    for i in 0..ref_len {
        ref_seq.push(b"ACGT"[(i * 7 + 3) % 4]);
    }
    let mut plain_fa = Vec::with_capacity(600);
    plain_fa.extend_from_slice(b">chr1\n");
    let mut split_upos = 0usize;
    for (idx, chunk) in ref_seq.chunks(64).enumerate() {
        if idx == 4 {
            split_upos = plain_fa.len();
        }
        plain_fa.extend_from_slice(chunk);
        plain_fa.push(b'\n');
    }
    let plain_fai = b"chr1\t512\t6\t64\t65\n".to_vec();

    let clean_fa_path = dir.join("clean_ref.fa");
    let clean_fai_path = dir.join("clean_ref.fa.fai");
    fs::write(&clean_fa_path, &plain_fa).ok()?;
    fs::write(&clean_fai_path, &plain_fai).ok()?;

    // Build a 2-block BGZF-compressed FASTA + matching GZI index in memory.
    let mut bgzf_fa = Vec::new();
    let block1_cpos = {
        let mut w = noodles_bgzf::io::Writer::new(&mut bgzf_fa);
        w.write_all(&plain_fa[..split_upos]).ok()?;
        w.flush().ok()?;
        let pos = w.position();
        w.write_all(&plain_fa[split_upos..]).ok()?;
        w.finish().ok()?;
        pos
    };
    let gzi_path = dir.join("init_tmp.gzi");
    let gzi_index = noodles_bgzf::gzi::Index::from(vec![(block1_cpos, split_upos as u64)]);
    noodles_bgzf::gzi::fs::write(&gzi_path, &gzi_index).ok()?;
    let bgzf_gzi = fs::read(&gzi_path).ok()?;
    let _ = fs::remove_file(&gzi_path);

    // Build a minimal valid BAM + BAI and CRAM + CRAI on `chr1`.
    let nz_len = std::num::NonZeroUsize::new(ref_len)?;
    let rs = HeaderMap::<RefSeqMap>::builder()
        .set_length(nz_len)
        .build()
        .ok()?;
    let header = Header::builder().add_reference_sequence("chr1", rs).build();

    let mut records = Vec::new();
    for (idx, start_1) in [(0usize, 40usize), (1usize, 140usize), (2usize, 240usize)] {
        let read_len = 40usize;
        let mut seq = ref_seq[start_1 - 1..start_1 - 1 + read_len].to_vec();
        seq[10] = if seq[10] == b'A' { b'G' } else { b'A' };
        let cigar: Cigar = match idx {
            0 => [Op::new(Kind::Match, read_len)].into_iter().collect(),
            1 => [
                Op::new(Kind::Match, 15),
                Op::new(Kind::Insertion, 2),
                Op::new(Kind::Match, 15),
                Op::new(Kind::Deletion, 2),
                Op::new(Kind::Match, 8),
            ]
            .into_iter()
            .collect(),
            _ => [
                Op::new(Kind::Match, 20),
                Op::new(Kind::Skip, 20),
                Op::new(Kind::Match, 20),
            ]
            .into_iter()
            .collect(),
        };
        let data: Data = [
            (Tag::EDIT_DISTANCE, BufValue::from(1i32)),
            (Tag::new(b'M', b'C'), BufValue::from("40M")),
            (Tag::new(b'M', b'Q'), BufValue::from(60i32)),
        ]
        .into_iter()
        .collect();
        let pos = noodles::core::Position::try_from(start_1).ok()?;
        let mate_pos = noodles::core::Position::try_from(start_1 + 80).ok()?;
        let mq = MappingQuality::try_from(60u8).ok()?;
        let rec = RecordBuf::builder()
            .set_name(format!("r{idx}"))
            .set_flags(Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_alignment_start(pos)
            .set_mapping_quality(mq)
            .set_cigar(cigar)
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(mate_pos)
            .set_template_length(120)
            .set_sequence(seq.into())
            .set_quality_scores(vec![35u8; read_len].into())
            .set_data(data)
            .build();
        records.push(rec);
    }

    let init_bam = dir.join("init_tmp.bam");
    let init_bai = dir.join("init_tmp.bam.bai");
    {
        let bgzf = noodles::bgzf::io::Writer::new(File::create(&init_bam).ok()?);
        let mut bam_writer = noodles::bam::io::Writer::from(bgzf);
        bam_writer.write_alignment_header(&header).ok()?;
        let mut indexer = Indexer::<LinearIndex>::new(14, 5)?;
        for rec in &records {
            let s_vp = bam_writer.get_ref().virtual_position();
            bam_writer.write_alignment_record(&header, rec).ok()?;
            let e_vp = bam_writer.get_ref().virtual_position();
            let s = rec.alignment_start()?;
            let span = rec.cigar().alignment_span().max(1);
            let e = noodles::core::Position::new(s.get() + span - 1)?;
            indexer
                .add_record(Some((0, s, e, true)), Chunk::new(s_vp, e_vp))
                .ok()?;
        }
        bam_writer.try_finish().ok()?;
        noodles::bam::bai::fs::write(&init_bai, &indexer.build(1)).ok()?;
    }
    let valid_bam = fs::read(&init_bam).ok()?;
    let valid_bai = fs::read(&init_bai).ok()?;
    let _ = fs::remove_file(&init_bam);
    let _ = fs::remove_file(&init_bai);

    let init_cram = dir.join("init_tmp.cram");
    let init_crai = dir.join("init_tmp.cram.crai");
    {
        let fasta_reader = noodles::fasta::io::IndexedReader::new(
            std::io::BufReader::new(File::open(&clean_fa_path).ok()?),
            noodles::fasta::fai::fs::read(&clean_fai_path).ok()?,
        );
        let repo = noodles::fasta::Repository::new(
            noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
        );
        let mut cram_writer = noodles::cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repo)
            .build_from_path(&init_cram)
            .ok()?;
        cram_writer.write_alignment_header(&header).ok()?;
        for rec in &records {
            cram_writer.write_alignment_record(&header, rec).ok()?;
        }
        cram_writer.try_finish(&header).ok()?;
        let crai = noodles::cram::fs::index(&init_cram).ok()?;
        noodles::cram::crai::fs::write(&init_crai, &crai).ok()?;
    }
    let valid_cram = fs::read(&init_cram).ok()?;
    let valid_crai = fs::read(&init_crai).ok()?;

    Some(FuzzTemplates {
        plain_fa,
        plain_fai,
        bgzf_fa,
        bgzf_gzi,
        valid_bam,
        valid_bai,
        valid_cram,
        valid_crai,
    })
}

fn mutate_bytes(template: &[u8], patch: &[u8], preserve_prefix: usize) -> Vec<u8> {
    if patch.is_empty() || template.is_empty() {
        return template.to_vec();
    }
    let mut out = template.to_vec();
    let mode = patch[0] % 4;
    match mode {
        0 => {
            // Byte-level XOR / overwrite mutations past `preserve_prefix`.
            let start = preserve_prefix.min(out.len());
            let span = out.len().saturating_sub(start);
            if span > 0 {
                for pair in patch[1..].chunks(2).take(32) {
                    let off = start + (usize::from(pair[0]) % span);
                    let mask = pair.get(1).copied().unwrap_or(0xff).max(1);
                    if let Some(slot) = out.get_mut(off) {
                        *slot ^= mask;
                    }
                }
            }
        }
        1 => {
            // Truncation at arbitrary byte offset.
            let trunc_len = preserve_prefix
                .saturating_add(usize::from(patch.get(1).copied().unwrap_or(0)))
                .min(out.len());
            out.truncate(trunc_len);
        }
        2 => {
            // Splice fuzzer bytes into the middle/tail while capping total size <= MAX_FUZZ_INPUT_BYTES.
            let keep = preserve_prefix.max(out.len() / 2).min(out.len());
            out.truncate(keep);
            let max_extra = MAX_FUZZ_INPUT_BYTES.saturating_sub(out.len());
            out.extend_from_slice(&patch[1..patch.len().min(1 + max_extra)]);
        }
        _ => {
            // Direct overwrite of slice starting at `preserve_prefix`.
            let start = preserve_prefix.min(out.len());
            for (dst, &src) in out[start..].iter_mut().zip(patch[1..].iter()) {
                *dst = src;
            }
        }
    }
    out
}

fn bgzf_compress(raw: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(raw.len().saturating_add(64));
    let mut w = noodles_bgzf::io::Writer::new(&mut out);
    let _ = w.write_all(raw);
    let _ = w.finish();
    out
}

/// Computes `(start_1, end_1)` strictly satisfying `end_1.saturating_sub(start_1) <= MAX_REGION_SPAN_BP`.
fn bounded_pos1_pair(b0: u8, b1: u8, b2: u8, b3: u8, mode: u8) -> (u32, u32) {
    let raw_start = u32::from(u16::from_le_bytes([b0, b1]));
    let span = u32::from(u16::from_le_bytes([b2, b3])) % (MAX_REGION_SPAN_BP + 1); // 0..=4096
    match mode % 4 {
        0 => {
            // Normal low coordinates overlapping chr1:1..=512
            let s = (raw_start % 600).saturating_add(1);
            let e = s.saturating_add(span);
            (s, e)
        }
        1 => {
            // Zero / inverted coordinates (end <= start -> end - start == 0 <= 4096)
            let s = raw_start % 1024;
            let e = s.saturating_sub(u32::from(b2 % 16));
            (s, e)
        }
        2 => {
            // High coordinates near i32::MAX
            let max_u32 = i32::MAX as u32;
            let s = max_u32
                .saturating_sub(MAX_REGION_SPAN_BP)
                .saturating_add(raw_start % 2048);
            let e = s.saturating_add(span).min(max_u32);
            (s, e)
        }
        _ => {
            let s = raw_start.max(1);
            let e = s.saturating_add(span);
            (s, e)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn exercise_alignment_and_runner(
    aln_path: &Path,
    ref_path: &Path,
    bed_path: &Path,
    chr: &str,
    start_1: u32,
    end_1: u32,
    ctrl: u8,
    sub_mode: u8,
) {
    debug_assert!(end_1.saturating_sub(start_1) <= MAX_REGION_SPAN_BP);

    let mut config = Config {
        threads: 1,
        reference: ref_path.to_string_lossy().into_owned(),
        max_region_len: 512,
        realign: u8::from((ctrl & 1) != 0),
        dedup: (ctrl & 2) != 0,
        chimeric: (ctrl & 4) != 0,
        pileup: (ctrl & 8) != 0,
        splice: (ctrl & 16) != 0,
        emit_intermediate: (ctrl & 32) != 0,
        fisher: (ctrl & 64) != 0,
        ..Default::default()
    };
    config.debug = (sub_mode & 0x40) != 0;

    let q_start = Pos1(start_1.max(1)).to_zero_based();
    let q_end_val = q_start.0.saturating_add(end_1.saturating_sub(start_1));
    let q_end = Pos0(q_end_val);
    debug_assert!(q_end.0.saturating_sub(q_start.0) <= MAX_REGION_SPAN_BP);

    if let Ok(mut feed) = query_region_iter(aln_path, chr, q_start, q_end, &config) {
        let mut preproc = RecordPreprocessor::new();
        let mut sam_rec = SamRecord::new();
        for _ in 0..MAX_DRAINED_RECORDS {
            match feed.next_into(&mut sam_rec) {
                Ok(true) => {
                    let _ = preproc.preprocess_record(&sam_rec, &config);
                }
                Ok(false) | Err(_) => break,
            }
        }
    }

    clear_bam_caches();

    let aln_str = aln_path.to_string_lossy().into_owned();
    config.bam = if (ctrl & 128) != 0 {
        format!("{aln_str}|{aln_str}")
    } else {
        aln_str
    };

    // Keep region span strictly <= MAX_REGION_SPAN_BP in both `-R` and amplicon BED modes.
    let reg_s = start_1.max(1);
    let reg_e = reg_s.saturating_add(end_1.saturating_sub(start_1).min(MAX_REGION_SPAN_BP));
    debug_assert!(reg_e.saturating_sub(reg_s) <= MAX_REGION_SPAN_BP);

    if (ctrl & 3) == 3 && !config.bam.contains('|') {
        let ins_s = reg_s.saturating_add(2).min(reg_e);
        let ins_e = reg_e.saturating_sub(2).max(ins_s);
        let bed_line = format!("{chr}\t{reg_s}\t{reg_e}\tGENE1\t0\t+\t{ins_s}\t{ins_e}\n");
        if fs::write(bed_path, bed_line.as_bytes()).is_ok() {
            config.region = None;
            config.bed = Some(bed_path.to_string_lossy().into_owned());
            config.amplicon = Some("10:0.95".to_string());
        } else {
            config.region = Some(format!("{chr}:{reg_s}-{reg_e}"));
        }
    } else {
        config.region = Some(format!("{chr}:{reg_s}-{reg_e}"));
    }

    let _ = run_with_writer(&config, std::io::sink());
}

fuzz_target!(|data: &[u8]| {
    if data.len() < 8 || data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    assert!(
        TEMPLATES.get_or_init(init_templates).is_some(),
        "init_templates failed"
    );
    let Some(tpl) = TEMPLATES.get_or_init(init_templates).as_ref() else {
        return;
    };

    clear_bam_caches();
    // Declare `temp_dir` before `_guard` so reverse drop order runs `CacheClearGuard::drop`
    // (`clear_bam_caches()`, closing pooled BAM/CRAM/FASTA FDs) before unlinking `temp_dir`.
    let Some(temp_dir) = create_fuzz_tempdir() else {
        return;
    };
    let _guard = CacheClearGuard;
    let dir = temp_dir.path();

    let selector = data[0];
    let sub_mode = data[1];
    let ctrl = data[2];
    let (start_1, end_1) = bounded_pos1_pair(data[3], data[4], data[5], data[6], data[7]);
    debug_assert!(end_1.saturating_sub(start_1) <= MAX_REGION_SPAN_BP);

    let payload = &data[8..];
    let mid = payload.len() / 2;
    let (part_a, part_b) = payload.split_at(mid);

    // Use bits [3:2] & [5] of `sub_mode` (independent of `sub_mode % 4` mutation-mode bits [1:0]
    // and bit [4] `use_mutated_ref`) so every index mutation mode queries `chr1` as well as
    // missing/empty contigs.
    let chr = match (sub_mode >> 2) & 3 {
        0..=2 => "chr1",
        _ if (sub_mode & 0x20) == 0 => "chr_missing",
        _ => "",
    };

    let clean_fa = dir.join("clean_ref.fa");
    let clean_fai = dir.join("clean_ref.fa.fai");
    let fuzz_fa = dir.join("fuzz_ref.fa");
    let fuzz_fai = dir.join("fuzz_ref.fa.fai");
    let fuzz_gz_fa = dir.join("fuzz_ref.fa.gz");
    let fuzz_gz_fai = dir.join("fuzz_ref.fa.gz.fai");
    let fuzz_gz_gzi = dir.join("fuzz_ref.fa.gz.gzi");
    let fuzz_bam = dir.join("fuzz_in.bam");
    let fuzz_bai = dir.join("fuzz_in.bam.bai");
    let fuzz_csi = dir.join("fuzz_in.bam.csi");
    let fuzz_cram = dir.join("fuzz_in.cram");
    let fuzz_crai = dir.join("fuzz_in.cram.crai");
    let fuzz_bed = dir.join("fuzz_regions.bed");

    match selector % 4 {
        0 => {
            // Hostile / mutated uncompressed FASTA + .fai
            match sub_mode % 3 {
                0 => {
                    let _ = fs::write(&fuzz_fa, part_a);
                    let _ = fs::write(&fuzz_fai, part_b);
                }
                1 => {
                    let fa_bytes = mutate_bytes(&tpl.plain_fa, part_a, 6);
                    let fuzzed_len = u32::from(u16::from_le_bytes([data[3], data[4]]));
                    let fuzzed_off = u32::from(data[5]);
                    let line_b = u32::from(data[6]).max(1);
                    let line_w = line_b.saturating_add(u32::from(data[7] % 4));
                    let fai_txt = format!("chr1\t{fuzzed_len}\t{fuzzed_off}\t{line_b}\t{line_w}\n");
                    let _ = fs::write(&fuzz_fa, fa_bytes);
                    let _ = fs::write(&fuzz_fai, fai_txt.as_bytes());
                }
                _ => {
                    let _ = fs::write(&fuzz_fa, mutate_bytes(&tpl.plain_fa, part_a, 0));
                    let _ = fs::write(&fuzz_fai, mutate_bytes(&tpl.plain_fai, part_b, 0));
                }
            }
            let _ = fetch_subsequence(&fuzz_fa, chr, Pos1(start_1), Pos1(end_1));
        }
        1 => {
            // Hostile / mutated BGZF-compressed FASTA + .fai + .gzi
            match sub_mode % 4 {
                0 => {
                    // Raw BGZF magic header + hostile payload + raw hostile GZI
                    let mut raw_bgzf = vec![
                        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
                        0x00,
                    ];
                    raw_bgzf.extend_from_slice(part_a);
                    let _ = fs::write(&fuzz_gz_fa, raw_bgzf);
                    let _ = fs::write(&fuzz_gz_fai, &tpl.plain_fai);
                    let _ = fs::write(&fuzz_gz_gzi, part_b);
                }
                1 => {
                    // Valid/mutated BGZF FASTA + structured hostile GZI (count + 16-byte pairs)
                    let fa_gz = mutate_bytes(&tpl.bgzf_fa, part_a, 18);
                    let _ = fs::write(&fuzz_gz_fa, fa_gz);
                    let _ = fs::write(&fuzz_gz_fai, mutate_bytes(&tpl.plain_fai, part_b, 0));
                    let n_entries = (part_b.len() / 16).min(16);
                    let declared = if (ctrl & 1) == 0 {
                        n_entries as u64
                    } else {
                        u64::from(ctrl)
                    };
                    let mut gzi_bytes = Vec::with_capacity(8 + n_entries * 16);
                    gzi_bytes.extend_from_slice(&declared.to_le_bytes());
                    gzi_bytes.extend_from_slice(&part_b[..n_entries * 16]);
                    let _ = fs::write(&fuzz_gz_gzi, gzi_bytes);
                }
                2 => {
                    // Missing GZI index path check (fresh `temp_dir` has no `.gzi` file)
                    let _ = fs::write(&fuzz_gz_fa, &tpl.bgzf_fa);
                    let _ = fs::write(&fuzz_gz_fai, &tpl.plain_fai);
                }
                _ => {
                    let _ = fs::write(&fuzz_gz_fa, mutate_bytes(&tpl.bgzf_fa, part_a, 0));
                    let _ = fs::write(&fuzz_gz_fai, mutate_bytes(&tpl.plain_fai, part_b, 0));
                    let _ = fs::write(&fuzz_gz_gzi, mutate_bytes(&tpl.bgzf_gzi, part_b, 0));
                }
            }
            let _ = fetch_subsequence(&fuzz_gz_fa, chr, Pos1(start_1), Pos1(end_1));
        }
        2 => {
            // Hostile / mutated BAM + .bai / .csi
            if fs::write(&clean_fa, &tpl.plain_fa).is_err()
                || fs::write(&clean_fai, &tpl.plain_fai).is_err()
            {
                return;
            }
            match sub_mode % 4 {
                0 => {
                    // Raw hostile BAM + raw hostile BAI
                    let _ = fs::write(&fuzz_bam, part_a);
                    let _ = fs::write(&fuzz_bai, part_b);
                }
                1 => {
                    // Mutated valid BAM + valid BAI
                    let _ = fs::write(&fuzz_bam, mutate_bytes(&tpl.valid_bam, payload, 4));
                    let _ = fs::write(&fuzz_bai, &tpl.valid_bai);
                }
                2 => {
                    // Valid BAM + mutated BAI
                    let _ = fs::write(&fuzz_bam, &tpl.valid_bam);
                    let _ = fs::write(&fuzz_bai, mutate_bytes(&tpl.valid_bai, payload, 4));
                }
                _ => {
                    // Mutated BAM + BGZF-wrapped hostile CSI index (`CSI\x01` prefix)
                    let _ = fs::write(&fuzz_bam, mutate_bytes(&tpl.valid_bam, part_a, 0));
                    let mut csi_raw = b"CSI\x01".to_vec();
                    csi_raw.extend_from_slice(part_b);
                    let _ = fs::write(&fuzz_csi, bgzf_compress(&csi_raw));
                }
            }
            exercise_alignment_and_runner(
                &fuzz_bam, &clean_fa, &fuzz_bed, chr, start_1, end_1, ctrl, sub_mode,
            );
        }
        _ => {
            // Hostile / mutated CRAM + .crai + reference FASTA
            let use_mutated_ref = (sub_mode & 0x10) != 0;
            let active_ref = if use_mutated_ref {
                let _ = fs::write(&fuzz_fa, mutate_bytes(&tpl.plain_fa, part_b, 6));
                let _ = fs::write(&fuzz_fai, mutate_bytes(&tpl.plain_fai, part_b, 0));
                &fuzz_fa
            } else {
                if fs::write(&clean_fa, &tpl.plain_fa).is_err()
                    || fs::write(&clean_fai, &tpl.plain_fai).is_err()
                {
                    return;
                }
                &clean_fa
            };

            match sub_mode % 4 {
                0 => {
                    // Raw hostile CRAM (with optional CRAM magic) + raw hostile CRAI
                    let mut cram_raw = Vec::with_capacity(part_a.len() + 6);
                    if (ctrl & 1) == 0 {
                        cram_raw.extend_from_slice(b"CRAM\x03\x00");
                    }
                    cram_raw.extend_from_slice(part_a);
                    let _ = fs::write(&fuzz_cram, cram_raw);
                    let _ = fs::write(&fuzz_crai, part_b);
                }
                1 => {
                    // Mutated valid CRAM container/slice bytes + valid CRAI
                    let keep_hdr = if (ctrl & 1) == 0 { 64 } else { 6 };
                    let _ = fs::write(&fuzz_cram, mutate_bytes(&tpl.valid_cram, payload, keep_hdr));
                    let _ = fs::write(&fuzz_crai, &tpl.valid_crai);
                }
                2 => {
                    // Valid/mutated CRAM + gzip-compressed hostile CRAI lines
                    let _ = fs::write(&fuzz_cram, mutate_bytes(&tpl.valid_cram, part_a, 64));
                    let crai_txt = if (ctrl & 1) == 0 {
                        let ref_id = i32::from(data[3] as i8);
                        let aln_start = i32::from(u16::from_le_bytes([data[4], data[5]]));
                        let aln_span = i32::from(u16::from_le_bytes([data[6], data[7]]));
                        let off = u64::from(part_b.first().copied().unwrap_or(0)) * 8;
                        let landmark = u64::from(part_b.get(1).copied().unwrap_or(0));
                        let slice_len = u64::from(part_b.get(2).copied().unwrap_or(64)) * 4;
                        format!(
                            "{ref_id}\t{aln_start}\t{aln_span}\t{off}\t{landmark}\t{slice_len}\n"
                        )
                        .into_bytes()
                    } else {
                        part_b.to_vec()
                    };
                    let _ = fs::write(&fuzz_crai, bgzf_compress(&crai_txt));
                }
                _ => {
                    // Mutated CRAM + byte-mutated CRAI
                    let _ = fs::write(&fuzz_cram, mutate_bytes(&tpl.valid_cram, part_a, 26));
                    let _ = fs::write(&fuzz_crai, mutate_bytes(&tpl.valid_crai, part_b, 2));
                }
            }
            exercise_alignment_and_runner(
                &fuzz_cram, active_ref, &fuzz_bed, chr, start_1, end_1, ctrl, sub_mode,
            );
        }
    }
});
