#![no_main]

use libfuzzer_sys::fuzz_target;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use std::collections::HashMap;
use vardictrust_core::bam::SamRecord;
use vardictrust_core::cigar::parse::CigarParser;
use vardictrust_core::config::Config;
use vardictrust_core::realign::{adj_snv, adjust_mnp, realign_indels};
use vardictrust_core::reference::{LoadedRegion, Reference};
use vardictrust_core::regions::Region;
use vardictrust_core::sv::find_all_sv_paths;
use vardictrust_core::tovars::ToVarsBuilder;
use vardictrust_core::types::Pos1;

fuzz_target!(|data: &[u8]| {
    if data.len() < 16 {
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

    let mut config = Config::default();
    config.realign = u8::from((mode_bits & 4) != 0);
    config.chimeric = (mode_bits & 8) != 0;
    config.dedup = (mode_bits & 16) != 0;
    config.minmatch = u32::from(data[6] % 40);
    config.extension_nuc = i32::from(data[7] as i8).clamp(-2, 15);
    if (mode_bits & 32) != 0 {
        config.trim_bases_after = i32::from(data[8]);
    }

    let payload = &data[9..];
    if payload.is_empty() {
        return;
    }

    let bases = [b'A', b'C', b'G', b'T', b'N'];
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
        CigarKind::SequenceMatch,
        CigarKind::SequenceMismatch,
    ];

    let mut parser = CigarParser::new(&reference, &region, &config, false, 0, 0, seed_max_rl);
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
        let flags = u16::from(chunk[0] & 0x13);
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
            mate_alignment_start: Some(Pos1(base_start)),
            template_length: i32::from(chunk[3] as i8) * 10,
            ..SamRecord::default()
        };
        let _ = parser.parse_cigar(&rec);
    }

    let mut var_data = parser.into_variation_data();
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

    let builder = ToVarsBuilder::new(&reference, &region, &chr_lengths, &config);
    let _ = builder.build_aligned_variants(&var_data);
});
