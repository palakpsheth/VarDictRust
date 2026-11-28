#![no_main]

use std::collections::HashMap;
use std::fs::{self, File};

use libfuzzer_sys::fuzz_target;
use noodles::csi::binning_index::index::reference_sequence::{bin::Chunk, index::LinearIndex};
use noodles::csi::binning_index::Indexer;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::cigar::{op::Kind as CigarKind, Op as CigarOp};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::{
    data::field::Value as BufValue, Cigar as CigarBuf, Data as DataBuf, RecordBuf,
};
use noodles::sam::header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap};
use noodles::sam::Header;
use vardictrust_core::bam::{clear_bam_caches, SamRecord};
use vardictrust_core::cigar::parse::CigarParser;
use vardictrust_core::config::{
    CompatProfile, Config, COMPAT_OFF, COMPAT_SWITCH_COUNT, COMPAT_VARDICTJAVA_1_8_3,
};
use vardictrust_core::postprocess::simple::simple_post_process;
use vardictrust_core::realign::{adj_snv, adjust_mnp, realign_indels};
use vardictrust_core::reference::{LoadedRegion, Reference};
use vardictrust_core::regions::Region;
use vardictrust_core::run_with_writer;
use vardictrust_core::sv::find_all_sv_paths;
use vardictrust_core::tovars::ToVarsBuilder;
use vardictrust_core::types::{Pos1, VardictError};
use vardictrust_core::variations::Variant;

/// Maximum accepted libFuzzer input payload (`4 KiB`) to bound per-iteration runtime.
const MAX_FUZZ_INPUT_BYTES: usize = 4096;

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
        .prefix("vdr-pr-fuzz-")
        .tempdir_in(&shm_root)
        .or_else(|_| tempfile::tempdir())
        .ok()
}

fn decode_compat_profile(sel: u8, mask64: u64) -> CompatProfile {
    match sel % 6 {
        0 => COMPAT_VARDICTJAVA_1_8_3,
        1 => COMPAT_OFF,
        2 => CompatProfile::only_switch((mask64 as usize) % COMPAT_SWITCH_COUNT)
            .unwrap_or(COMPAT_OFF),
        3 => CompatProfile::without_switch((mask64 as usize) % COMPAT_SWITCH_COUNT)
            .unwrap_or(COMPAT_OFF),
        _ => CompatProfile::from_mask(mask64),
    }
}

#[inline]
fn assert_accounting_invariants(
    v: &Variant,
    compat: &CompatProfile,
    is_ref: bool,
    check_postproc_shift3: bool,
) {
    if compat.reconcile_vd_dp || compat.reconcile_rd_vd_dp {
        assert!(
            0 <= v.pos_coverage && v.pos_coverage <= v.total_pos_cov,
            "VD out of [0, DP]: vd={}, dp={}, is_ref={is_ref}",
            v.pos_coverage,
            v.total_pos_cov
        );
        if is_ref {
            assert!(
                0 <= v.ref_fwd_cov
                    && 0 <= v.ref_rev_cov
                    && v.ref_fwd_cov.saturating_add(v.ref_rev_cov) <= v.total_pos_cov,
                "ref_variant RD out of [0, DP]: rf={}, rr={}, dp={}",
                v.ref_fwd_cov,
                v.ref_rev_cov,
                v.total_pos_cov
            );
        }
        assert!(
            v.frequency.is_finite() && (0.0..=1.0).contains(&v.frequency),
            "AF out of [0, 1]: {}",
            v.frequency
        );
        assert!(
            v.extra_freq.is_finite() && (0.0..=1.0).contains(&v.extra_freq),
            "ADJAF out of [0, 1]: {}",
            v.extra_freq
        );
        assert!(
            0 <= v.hicnt && v.hicnt <= v.hicov,
            "HICNT out of [0, HICOV]: hicnt={}, hicov={}",
            v.hicnt,
            v.hicov
        );
        assert!(
            v.hi_freq.is_finite() && (0.0..=1.0).contains(&v.hi_freq),
            "HIAF out of [0, 1]: {}",
            v.hi_freq
        );
    }

    if compat.reconcile_rd_vd_dp {
        let rd = i64::from(v.ref_fwd_cov.max(0)) + i64::from(v.ref_rev_cov.max(0));
        let vd = i64::from(v.pos_coverage.max(0));
        assert!(
            v.ref_fwd_cov >= 0
                && v.ref_rev_cov >= 0
                && v.pos_coverage >= 0
                && rd + vd <= i64::from(v.total_pos_cov),
            "reconcile_rd_vd_dp violated RD ({rd}) + VD ({vd}) <= DP ({})",
            v.total_pos_cov
        );
    }

    if compat.reconcile_strand_and_hq_sums && !v.null_record {
        assert!(
            v.fwd >= 0
                && v.rev >= 0
                && i64::from(v.fwd) + i64::from(v.rev) == i64::from(v.pos_coverage),
            "reconcile_strand_and_hq_sums violated VD ({}) == fwd ({}) + rev ({})",
            v.pos_coverage,
            v.fwd,
            v.rev
        );
        if is_ref {
            assert!(
                0 <= v.hicnt && v.hicnt <= v.hicov.max(0),
                "ref_variant violated 0 <= hicnt ({}) <= hicov ({})",
                v.hicnt,
                v.hicov
            );
        } else {
            assert!(
                0 <= v.hicnt && v.hicnt <= v.pos_coverage,
                "reconcile_strand_and_hq_sums violated 0 <= hicnt ({}) <= VD ({})",
                v.hicnt,
                v.pos_coverage
            );
        }
    }

    if compat.nonnegative_subtracted_means && !v.null_record {
        assert!(
            v.mean_position >= 0.0
                && v.mean_quality >= 0.0
                && v.mean_mapping_quality >= 0.0
                && v.nm >= 0.0,
            "nonnegative_subtracted_means violated non-negative means: pmean={}, qmean={}, mq={}, nm={}",
            v.mean_position,
            v.mean_quality,
            v.mean_mapping_quality,
            v.nm
        );
    }

    if check_postproc_shift3
        && compat.indel_only_shift3
        && !v.null_record
        && !v.ref_allele.is_empty()
        && !v.var_allele.contains('<')
        && v.ref_allele.len() == v.var_allele.len()
    {
        assert_eq!(
            v.shift3, 0,
            "indel_only_shift3 violated SHIFT3 == 0 on equal-length allele ref={} alt={}",
            v.ref_allele, v.var_allele
        );
    }
}

fn write_indexed_bam_bufs(
    bam_path: &std::path::Path,
    chr_len: usize,
    mut raw_recs: Vec<(usize, RecordBuf)>,
) -> Option<()> {
    raw_recs.sort_by_key(|&(s, _)| s);
    let nz_len = std::num::NonZeroUsize::new(chr_len)?;
    let rs = HeaderMap::<RefSeqMap>::builder()
        .set_length(nz_len)
        .build()
        .ok()?;
    let header = Header::builder().add_reference_sequence("chr1", rs).build();
    let bai_path = bam_path.with_extension("bam.bai");

    let bgzf = noodles::bgzf::io::Writer::new(File::create(bam_path).ok()?);
    let mut bam_writer = noodles::bam::io::Writer::from(bgzf);
    bam_writer.write_alignment_header(&header).ok()?;
    let mut indexer = Indexer::<LinearIndex>::new(14, 5)?;
    for (_, rec) in &raw_recs {
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
    noodles::bam::bai::fs::write(&bai_path, &indexer.build(1)).ok()
}

fn exercise_tiled_vs_untiled(base_cfg: &Config, payload: &[u8], ctrl_byte: u8) {
    clear_bam_caches();
    let Some(temp_dir) = create_fuzz_tempdir() else {
        return;
    };
    let _guard = CacheClearGuard;
    let dir = temp_dir.path();

    // Build a 320-bp single-line FASTA + .fai so `tile_size = Some(64)` splits `chr1:16-272` into 4 tiles.
    let ref_len = 320usize;
    let bases = [b'A', b'C', b'G', b'T'];
    let ref_seq: Vec<u8> = (0..ref_len)
        .map(|i| bases[(usize::from(payload[i % payload.len()]) + i * 3) % 4])
        .collect();
    let fa_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let mut fa_bytes = Vec::with_capacity(ref_len + 8);
    fa_bytes.extend_from_slice(b">chr1\n");
    fa_bytes.extend_from_slice(&ref_seq);
    fa_bytes.push(b'\n');
    let fai_txt = format!("chr1\t{ref_len}\t6\t{ref_len}\t{}\n", ref_len + 1);
    if fs::write(&fa_path, &fa_bytes).is_err() || fs::write(&fai_path, fai_txt.as_bytes()).is_err()
    {
        return;
    }

    let build_sample_reads = |salt: usize| -> Option<Vec<(usize, RecordBuf)>> {
        let mut out = Vec::with_capacity(8);
        for (idx, chunk) in payload.chunks(12).take(8).enumerate() {
            let start_1 = 10 + (usize::from(chunk[0]).wrapping_add(salt * 17) % 240); // 10..=249
            let read_len = 36usize;
            if start_1 + read_len > ref_len {
                continue;
            }
            let mut seq = ref_seq[start_1 - 1..start_1 - 1 + read_len].to_vec();
            let mut nm = 0i32;
            let edit_kind = chunk.get(1).copied().unwrap_or(0) % 4;
            let cigar: CigarBuf = match edit_kind {
                0 => {
                    // Clean reference or single SNV read (`36M` or spec-valid `=`/`X`)
                    let mut snv_off = None;
                    if (chunk.get(2).copied().unwrap_or(0) & 1) == 0 {
                        let off = usize::from(chunk.get(3).copied().unwrap_or(10)) % read_len;
                        seq[off] = if seq[off] == b'A' { b'G' } else { b'A' };
                        nm = 1;
                        snv_off = Some(off);
                    }
                    if base_cfg.compat.eqx_cigar_read_length && idx % 2 == 1 {
                        let mut ops = Vec::with_capacity(3);
                        if let Some(off) = snv_off {
                            if off > 0 {
                                ops.push(CigarOp::new(CigarKind::SequenceMatch, off));
                            }
                            ops.push(CigarOp::new(CigarKind::SequenceMismatch, 1));
                            if off + 1 < read_len {
                                ops.push(CigarOp::new(
                                    CigarKind::SequenceMatch,
                                    read_len - 1 - off,
                                ));
                            }
                        } else {
                            ops.push(CigarOp::new(CigarKind::SequenceMatch, read_len));
                        }
                        ops.into_iter().collect()
                    } else {
                        [CigarOp::new(CigarKind::Match, read_len)]
                            .into_iter()
                            .collect()
                    }
                }
                1 => {
                    // 2-bp insertion inside read (`16M2I18M` -> query=36, ref=34)
                    if start_1 + 34 > ref_len {
                        continue;
                    }
                    let ins0 = bases[usize::from(chunk.get(2).copied().unwrap_or(0)) % 4];
                    let ins1 = bases[usize::from(chunk.get(3).copied().unwrap_or(1)) % 4];
                    seq.clear();
                    seq.extend_from_slice(&ref_seq[start_1 - 1..start_1 - 1 + 16]);
                    seq.push(ins0);
                    seq.push(ins1);
                    seq.extend_from_slice(&ref_seq[start_1 - 1 + 16..start_1 - 1 + 34]);
                    nm = 2;
                    [
                        CigarOp::new(CigarKind::Match, 16),
                        CigarOp::new(CigarKind::Insertion, 2),
                        CigarOp::new(CigarKind::Match, 18),
                    ]
                    .into_iter()
                    .collect()
                }
                2 => {
                    // 2-bp deletion inside read (`18M2D18M` -> query=36, ref=38)
                    if start_1 + 38 > ref_len {
                        continue;
                    }
                    seq.clear();
                    seq.extend_from_slice(&ref_seq[start_1 - 1..start_1 - 1 + 18]);
                    seq.extend_from_slice(&ref_seq[start_1 - 1 + 20..start_1 - 1 + 38]);
                    nm = 2;
                    [
                        CigarOp::new(CigarKind::Match, 18),
                        CigarOp::new(CigarKind::Deletion, 2),
                        CigarOp::new(CigarKind::Match, 18),
                    ]
                    .into_iter()
                    .collect()
                }
                _ => {
                    // Outer hard-clip + soft-clip (`4H6S30M2H` -> query=36, ref=30)
                    if start_1 + 30 > ref_len {
                        continue;
                    }
                    let sc_seed = usize::from(chunk.get(2).copied().unwrap_or(0));
                    seq.clear();
                    for k in 0..6usize {
                        seq.push(bases[(sc_seed + k + 1) % 4]);
                    }
                    seq.extend_from_slice(&ref_seq[start_1 - 1..start_1 - 1 + 30]);
                    nm = 0;
                    [
                        CigarOp::new(CigarKind::HardClip, 4),
                        CigarOp::new(CigarKind::SoftClip, 6),
                        CigarOp::new(CigarKind::Match, 30),
                        CigarOp::new(CigarKind::HardClip, 2),
                    ]
                    .into_iter()
                    .collect()
                }
            };
            let pos = noodles::core::Position::try_from(start_1).ok()?;
            let mq = MappingQuality::try_from(60u8).ok()?;
            let flags = if (chunk.get(4).copied().unwrap_or(0) & 1) != 0 {
                Flags::REVERSE_COMPLEMENTED
            } else {
                Flags::empty()
            };
            let data_buf: DataBuf = [(Tag::EDIT_DISTANCE, BufValue::from(nm))]
                .into_iter()
                .collect();
            let rec = RecordBuf::builder()
                .set_name(format!("s{salt}_r{idx}"))
                .set_flags(flags)
                .set_reference_sequence_id(0)
                .set_alignment_start(pos)
                .set_mapping_quality(mq)
                .set_cigar(cigar)
                .set_sequence(seq.into())
                .set_quality_scores(vec![35u8; read_len].into())
                .set_data(data_buf)
                .build();
            out.push((start_1, rec));
        }
        Some(out)
    };

    let Some(t_reads) = build_sample_reads(0) else {
        return;
    };
    if t_reads.is_empty() {
        return;
    }
    let tumor_bam = dir.join("tumor.bam");
    if write_indexed_bam_bufs(&tumor_bam, ref_len, t_reads).is_none() {
        return;
    }

    let paired = (ctrl_byte & 1) != 0;
    let is_amplicon = (ctrl_byte & 2) != 0;
    let emit_inter = (ctrl_byte & 4) != 0;

    let bam_spec = if paired {
        let Some(n_reads) = build_sample_reads(1) else {
            return;
        };
        let normal_bam = dir.join("normal.bam");
        if write_indexed_bam_bufs(&normal_bam, ref_len, n_reads).is_none() {
            return;
        }
        format!("{}|{}", tumor_bam.display(), normal_bam.display())
    } else {
        tumor_bam.to_string_lossy().into_owned()
    };

    let mut untiled_cfg = Config {
        reference: fa_path.to_string_lossy().into_owned(),
        bam: bam_spec,
        sample: if paired {
            "TUMOR|NORMAL".to_string()
        } else {
            "SAMPLE".to_string()
        },
        sample_explicit: true,
        threads: 1,
        realign: base_cfg.realign,
        chimeric: base_cfg.chimeric,
        dedup: base_cfg.dedup,
        three_prime: base_cfg.three_prime,
        unique_mode_alignment_enabled: base_cfg.unique_mode_alignment_enabled,
        pileup: (ctrl_byte & 8) != 0,
        fisher: (ctrl_byte & 16) != 0,
        vcf_all_alleles: (ctrl_byte & 32) != 0,
        emit_intermediate: emit_inter,
        compat: base_cfg.compat,
        tile_size: None,
        ..Config::default()
    };
    untiled_cfg.apply_pileup_overrides();

    if is_amplicon {
        let bed_path = dir.join("amp.bed");
        let bed_txt = "chr1\t20\t140\tAMP1\t0\t+\t28\t132\nchr1\t120\t260\tAMP2\t0\t+\t128\t252\n";
        if fs::write(&bed_path, bed_txt.as_bytes()).is_err() {
            return;
        }
        untiled_cfg.bed = Some(bed_path.to_string_lossy().into_owned());
        untiled_cfg.region = None;
        untiled_cfg.amplicon = Some("10:0.95".to_string());
        untiled_cfg.zero_based = true;
        untiled_cfg.zero_based_defined = true;
        untiled_cfg.col_chr = 0;
        untiled_cfg.col_start = 1;
        untiled_cfg.col_end = 2;
        untiled_cfg.col_seg_start = 6;
        untiled_cfg.col_seg_end = 7;
        untiled_cfg.col_gene = Some(3);
        untiled_cfg.column_for_chromosome = Some(0);
    } else {
        untiled_cfg.region = Some("chr1:16-272:FUZZ".to_string());
    }

    let mut tiled_cfg = untiled_cfg.clone();
    tiled_cfg.tile_size = Some(64);

    let mut untiled_buf = Vec::new();
    let mut tiled_buf = Vec::new();
    let untiled_res = run_with_writer(&untiled_cfg, &mut untiled_buf);
    let tiled_res = run_with_writer(&tiled_cfg, &mut tiled_buf);

    if paired && is_amplicon && !untiled_cfg.compat.paired_somatic_amplicon && !emit_inter {
        assert!(
            matches!(untiled_res, Err(VardictError::Unsupported(_))),
            "expected VardictError::Unsupported on paired -a VCF under !paired_somatic_amplicon, got {untiled_res:?}"
        );
        assert!(
            matches!(tiled_res, Err(VardictError::Unsupported(_))),
            "expected VardictError::Unsupported on tiled paired -a VCF under !paired_somatic_amplicon, got {tiled_res:?}"
        );
        return;
    }

    let () = untiled_res.expect("untiled synthetic run_with_writer must succeed");
    let () = tiled_res.expect("tiled synthetic run_with_writer (tile_size=64) must succeed");

    let untiled_str = String::from_utf8(untiled_buf).expect("valid UTF-8 untiled output");
    let tiled_str = String::from_utf8(tiled_buf).expect("valid UTF-8 tiled output");

    if !emit_inter || untiled_cfg.compat.coordinate_sorted_vars {
        assert_eq!(
            tiled_str, untiled_str,
            "Tiled (tile_size=64) != Untiled (emit_intermediate={emit_inter}, compat={:?})",
            untiled_cfg.compat
        );
    } else {
        let mut u_lines: Vec<&str> = untiled_str.lines().collect();
        let mut t_lines: Vec<&str> = tiled_str.lines().collect();
        u_lines.sort_unstable();
        t_lines.sort_unstable();
        assert_eq!(
            t_lines, u_lines,
            "Tiled (tile_size=64) sorted .var multiset != Untiled"
        );
    }
}

fuzz_target!(|data: &[u8]| {
    if data.len() < 20 || data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    let mode_bits = data[0];
    let high_coord_mode = (mode_bits & 1) != 0;
    let small_max_rl_mode = (mode_bits & 2) != 0;

    // Clamp reference/region span strictly to <= 4096 bp, Pos1 >= 1, start <= end (AR-19, AR-29).
    let span_len = ((usize::from(u16::from_le_bytes([data[1], data[2]])) % 4096) + 1).min(4096);
    let base_start: u32 = if high_coord_mode {
        // Cover coordinates near i32::MAX - 2048 without exceeding i32::MAX.
        let max_u32 = i32::MAX as u32;
        let anchor = max_u32.saturating_sub(2048);
        let jitter = u32::from(data[3]);
        anchor
            .saturating_sub(jitter)
            .min(max_u32.saturating_sub(span_len as u32).saturating_add(1))
            .max(1)
    } else {
        (u32::from(u16::from_le_bytes([data[3], data[4]])) % 2048).saturating_add(1)
    };
    let base_end = base_start
        .saturating_add(span_len as u32 - 1)
        .min(i32::MAX as u32)
        .max(base_start);
    let eff_len = (base_end - base_start + 1) as usize;

    let seed_max_rl: i32 = if small_max_rl_mode {
        i32::from(data[5] % 21) // 0..=20
    } else {
        i32::from(data[5]).saturating_add(21)
    };

    let mask64 = u64::from_le_bytes([
        data[11], data[12], data[13], data[14], data[15], data[16], data[17], data[18],
    ]);
    let compat_sel = (mode_bits >> 6).wrapping_add(data[10] >> 5);

    let mut config = Config::default();
    config.compat = decode_compat_profile(compat_sel, mask64);
    if (mode_bits & 128) != 0 && (data[10] & 1) != 0 {
        config.amplicon = Some("10:0.95".to_string());
    }
    config.realign = u8::from((mode_bits & 4) != 0);
    config.chimeric = (mode_bits & 8) != 0;
    config.dedup = (mode_bits & 16) != 0;
    config.minmatch = u32::from(data[6] % 40);
    config.extension_nuc = i32::from(data[7] as i8).clamp(-2, 15);
    if (mode_bits & 32) != 0 {
        config.trim_bases_after = i32::from(data[8]);
    }

    let extra_bits = data[9];
    let extra_param = data[10];
    let payload = &data[19..];
    if payload.is_empty() {
        return;
    }

    let bases = [b'A', b'C', b'G', b'T', b'N'];
    config.splice = (extra_bits & 1) != 0;
    config.debug = (extra_bits & 64) != 0;
    config.unique_mode_alignment_enabled = (extra_param & 2) != 0;
    config.three_prime = (extra_param & 4) != 0;
    config.vcf_all_alleles = (extra_param & 8) != 0;
    config.pileup = (extra_param & 16) != 0;
    if config.pileup {
        config.apply_pileup_overrides();
    }
    if (extra_bits & 2) != 0 {
        let cut_site = if (extra_bits & 16) != 0 {
            match extra_param & 3 {
                0 => 0,
                1 => -1,
                2 => i32::MIN,
                _ => i32::MAX,
            }
        } else {
            let cut_off = u32::from(extra_param) % (eff_len as u32).max(1);
            base_start.saturating_add(cut_off).min(i32::MAX as u32) as i32
        };
        config.crispr_cutting_site = Some(cut_site);
    }
    if (extra_bits & 4) != 0 {
        config.crispr_flanking_size = Some(i32::from(extra_param as i8));
    }
    if (extra_bits & 8) != 0 {
        let ad_len = (usize::from(extra_param & 0x0f)).saturating_add(4); // 4..=19 bp
        let ad_seq: String = (0..ad_len)
            .map(|i| bases[usize::from(payload[i % payload.len()]) % 4] as char)
            .collect();
        config.adaptor.push(ad_seq);
        if (extra_bits & 32) != 0 {
            config.adaptor.push("AGATCGGAAGAGC".to_string());
        }
    }
    let ref_seq: Vec<u8> = (0..eff_len)
        .map(|i| bases[usize::from(payload[i % payload.len()]) % bases.len()])
        .collect();

    let mut reference = Reference::new();
    reference.loaded_regions.push(LoadedRegion {
        chr: "chr1".to_string(),
        sequence_start: base_start,
        sequence_end: base_end,
    });
    reference.region_sequences.push(ref_seq.clone());
    reference
        .seed
        .index_sequence(&ref_seq, base_start, base_end, base_end);

    let sub_off = usize::from(data[8]) % eff_len;
    let reg_start = base_start.saturating_add(sub_off as u32).min(base_end);
    let reg_end = base_end.max(reg_start);
    let region = Region {
        chr: "chr1".to_string(),
        start: Pos1(reg_start),
        end: Pos1(reg_end),
        gene: "FUZZ".to_string(),
        amplicon: None,
    };

    let mut chr_lengths = HashMap::new();
    chr_lengths.insert("chr1".to_string(), base_end as usize);

    let cigar_kinds = [
        CigarKind::Match,
        CigarKind::Insertion,
        CigarKind::Deletion,
        CigarKind::Skip,
        CigarKind::SoftClip,
        CigarKind::HardClip,
        CigarKind::Pad,
        CigarKind::SequenceMatch,
        CigarKind::SequenceMismatch,
    ];

    let mut pre = vardictrust_core::bam::RecordPreprocessor::new();
    let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, seed_max_rl);
    let mut raw_records: Vec<SamRecord> = Vec::with_capacity(10);
    for (rec_idx, chunk) in payload.chunks(24).take(8).enumerate() {
        if chunk.len() < 6 {
            break;
        }
        let rel_pos = u32::from(u16::from_le_bytes([chunk[0], chunk[1]])) % (eff_len as u32);
        let aln_start = base_start.saturating_add(rel_pos).min(base_end);
        let mut cigar = Vec::new();
        let mut qlen = 0usize;
        let mut rlen = 0u32;
        for pair in chunk[2..].chunks(2).take(6) {
            let kind = cigar_kinds[usize::from(pair[0]) % cigar_kinds.len()];
            let raw_len = u32::from(pair.get(1).copied().unwrap_or(1));
            let op_len = if kind == CigarKind::Deletion && (pair[0] & 0x80) != 0 {
                // Exercise huge corrupt deletion coverage clipping without OOM (Fix #3)
                raw_len.saturating_mul(1_000_000).max(1)
            } else {
                (raw_len % 64).max(1)
            };
            if matches!(
                kind,
                CigarKind::Match
                    | CigarKind::Insertion
                    | CigarKind::SoftClip
                    | CigarKind::SequenceMatch
                    | CigarKind::SequenceMismatch
            ) {
                if qlen.saturating_add(op_len as usize) > 512 {
                    break;
                }
                qlen += op_len as usize;
            }
            if kind.consumes_reference() {
                rlen = rlen.saturating_add(op_len);
            }
            cigar.push((op_len, kind));
        }
        if cigar.is_empty() || qlen == 0 {
            continue;
        }
        let seq: Vec<u8> = (0..qlen)
            .map(|i| bases[usize::from(chunk[(i + 2) % chunk.len()]) % bases.len()])
            .collect();
        let quals: Vec<u8> = (0..qlen)
            .map(|i| b'!' + (chunk[(i + 3) % chunk.len()] % 41))
            .collect();
        let aln_end = aln_start
            .saturating_add(rlen.saturating_sub(1))
            .min(i32::MAX as u32);
        let flags = u16::from(chunk[0] & 0x1b);
        let rec = SamRecord {
            name: format!("r{rec_idx}").into_bytes(),
            flags,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(aln_start)),
            alignment_end: Some(Pos1(aln_end.max(aln_start))),
            mapping_quality: Some(chunk[1] % 61),
            cigar,
            sequence: seq,
            quality_scores: quals,
            nm: Some(i32::from(chunk[2] % 12)),
            mate_reference_name: Some("=".to_string()),
            mate_alignment_start: Some(Pos1(if (extra_bits & 128) != 0 {
                base_start.max(10)
            } else {
                u32::from(chunk[2] % 9).saturating_add(1)
            })),
            template_length: i32::from(chunk[3] as i8) * 10,
            ..SamRecord::default()
        };
        raw_records.push(rec);
    }

    // When `-u` (`unique_mode_alignment_enabled`) is active, inject an agreed overlapping FR
    // mate pair sharing `qname` so deferred `finalize_overlap_strand_rebalance()` is exercised.
    if config.unique_mode_alignment_enabled && eff_len >= 40 {
        let p1 = reg_start;
        let p2 = reg_start.saturating_add(8).min(base_end);
        let m_len = 24u32
            .min(base_end.saturating_sub(p2).saturating_add(1))
            .max(1);
        let mut m_seq = vec![b'G'; m_len as usize];
        if let Some(first) = m_seq.first_mut() {
            *first = b'T';
        }
        let fwd_mate = SamRecord {
            name: b"ovl_pair".to_vec(),
            flags: 0x1 | 0x2 | 0x20 | 0x40,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(p1)),
            alignment_end: Some(Pos1(p1.saturating_add(m_len - 1))),
            mapping_quality: Some(60),
            cigar: vec![(m_len, CigarKind::Match)],
            sequence: m_seq.clone(),
            quality_scores: vec![b'I'; m_len as usize],
            nm: Some(1),
            mate_reference_name: Some("=".to_string()),
            mate_alignment_start: Some(Pos1(p2)),
            template_length: (p2.saturating_sub(p1).saturating_add(m_len)) as i32,
            ..SamRecord::default()
        };
        let rev_mate = SamRecord {
            name: b"ovl_pair".to_vec(),
            flags: 0x1 | 0x2 | 0x10 | 0x80,
            alignment_start: Some(Pos1(p2)),
            alignment_end: Some(Pos1(p2.saturating_add(m_len - 1))),
            mate_alignment_start: Some(Pos1(p1)),
            template_length: -((p2.saturating_sub(p1).saturating_add(m_len)) as i32),
            ..fwd_mate.clone()
        };
        raw_records.push(fwd_mate);
        raw_records.push(rev_mate);
    }

    raw_records.sort_by_key(|r| r.alignment_start.map_or(0, |p| p.0));
    for rec in &raw_records {
        if pre.preprocess_record(rec, &config) {
            let _ = parser.parse_cigar(rec);
        }
    }

    // Production pipeline order (`lib.rs:1305-1339`): extract unrebalanced pileup state first,
    // run `adjust_mnp` -> `realign_indels` -> `find_all_sv_paths` -> `adj_snv`, then finalize
    // `-u` mate-overlap strand rebalancing ONCE after all `sub_dir` subtractions finish.
    let mut var_data = parser.into_variation_data_unrebalanced();
    if small_max_rl_mode {
        var_data.max_read_length = seed_max_rl; // force 0..=20 into realign_indels & find_all_sv_paths
    }
    adjust_mnp(
        &mut reference,
        "chr1",
        base_end,
        &region,
        &config,
        &mut var_data,
    );
    if config.realign == 1 {
        realign_indels(
            None,
            &mut reference,
            "chr1",
            base_end,
            &region,
            &config,
            &mut var_data,
        );
    }
    let _ = find_all_sv_paths(
        &mut var_data,
        &[],
        &mut reference,
        "chr1",
        base_end,
        &region,
        &config,
        &chr_lengths,
    );
    adj_snv(&reference, "chr1", &config, &mut var_data);
    var_data.finalize_overlap_strand_rebalance();

    let builder = ToVarsBuilder::new(&reference, &region, &chr_lengths, &config);
    let aligned = builder.build_aligned_variants(&var_data);

    for (_, vars) in aligned.iter() {
        for v in &vars.variants {
            assert_accounting_invariants(v, &config.compat, false, false);
        }
        if let Some(ref rv) = vars.ref_variant {
            assert_accounting_invariants(rv, &config.compat, true, false);
        }
    }

    if config.amplicon.is_none() {
        let mut var_cfg = config.clone();
        var_cfg.emit_intermediate = true;
        let var_rows = simple_post_process(&var_cfg, &region, &aligned, &var_data.splice);
        for v in &var_rows {
            let is_ref = v.var_type.is_empty() && v.pos_coverage == 0;
            assert_accounting_invariants(v, &var_cfg.compat, is_ref, true);
        }
        if var_cfg.compat.coordinate_sorted_vars {
            let mut prev_snv_pos = i32::MIN;
            for v in &var_rows {
                if v.description.len() == 1
                    && !v.description.starts_with(['+', '-'])
                    && v.var_type == "SNV"
                    && v.sv.is_empty()
                    && !v.null_record
                {
                    assert!(
                        v.start_pos >= prev_snv_pos,
                        "coordinate_sorted_vars SNV .var start_pos regressed: {} < {}",
                        v.start_pos,
                        prev_snv_pos
                    );
                    prev_snv_pos = v.start_pos;
                }
            }
        }

        let mut vcf_cfg = config.clone();
        vcf_cfg.emit_intermediate = false;
        let vcf_rows = simple_post_process(&vcf_cfg, &region, &aligned, &var_data.splice);
        for v in &vcf_rows {
            let is_ref = v.var_type.is_empty() && v.pos_coverage == 0;
            assert_accounting_invariants(v, &vcf_cfg.compat, is_ref, true);
        }
        if vcf_cfg.compat.coordinate_sorted_vars && !vcf_cfg.vcf_all_alleles {
            let mut prev_pos = i32::MIN;
            for v in &vcf_rows {
                assert!(
                    v.start_pos >= prev_pos,
                    "coordinate_sorted_vars VCF start_pos regressed: {} < {}",
                    v.start_pos,
                    prev_pos
                );
                prev_pos = v.start_pos;
            }
        }
    }

    // 1/8-gated end-to-end `run_with_writer` `Tiled (tile_size = Some(64)) == Untiled` check.
    if (mode_bits & 0x07) == 0 && !high_coord_mode && eff_len <= 512 && !config.splice {
        exercise_tiled_vs_untiled(&config, payload, extra_param);
    }
});
