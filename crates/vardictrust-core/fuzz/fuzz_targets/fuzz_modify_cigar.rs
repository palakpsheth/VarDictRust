#![no_main]

use libfuzzer_sys::fuzz_target;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use vardictrust_core::cigar::modify::CigarModifier;
use vardictrust_core::cigar::parse::query_consumed_len;
use vardictrust_core::cigar::parse_cigar;
use vardictrust_core::config::{
    CompatProfile, Config, COMPAT_OFF, COMPAT_SWITCH_COUNT, COMPAT_VARDICTJAVA_1_8_3,
};
use vardictrust_core::reference::SeedMap;

fn decode_compat_profile(mode: u8, mask_bytes: &[u8]) -> CompatProfile {
    let mut raw = [0u8; 8];
    for (i, slot) in raw.iter_mut().enumerate() {
        *slot = mask_bytes
            .get(i % mask_bytes.len().max(1))
            .copied()
            .unwrap_or(0);
    }
    let mask64 = u64::from_le_bytes(raw);
    match mode % 6 {
        0 => COMPAT_VARDICTJAVA_1_8_3,
        1 => COMPAT_OFF,
        2 => CompatProfile::only_switch((mask64 as usize) % COMPAT_SWITCH_COUNT)
            .unwrap_or(COMPAT_OFF),
        3 => CompatProfile::without_switch((mask64 as usize) % COMPAT_SWITCH_COUNT)
            .unwrap_or(COMPAT_OFF),
        _ => CompatProfile::from_mask(mask64),
    }
}

fuzz_target!(|data: &[u8]| {
    if data.len() < 12 {
        return;
    }

    let pos = i32::from_le_bytes([data[0], data[1], 0, 0]).clamp(1, 512);
    let max_rl = if data[2] <= 20 {
        i32::from(data[2])
    } else {
        i32::from(data[2]).clamp(21, 300)
    };
    let compat_mode = data[3] >> 5;
    let offset = i32::from(data[3] & 0x1f).clamp(1, 64);

    let payload = &data[4..];
    let third = payload.len() / 3;
    if third == 0 {
        return;
    }

    let bases = [b'A', b'C', b'G', b'T', b'N'];
    let ref_seq: Vec<u8> = payload[..third]
        .iter()
        .take(4096)
        .map(|&b| bases[(b as usize) % bases.len()])
        .collect();
    if ref_seq.is_empty() {
        return;
    }
    let query_seq: Vec<u8> = payload[third..2 * third]
        .iter()
        .take(4096)
        .map(|&b| bases[(b as usize) % bases.len()])
        .collect();
    let query_qual: Vec<u8> = payload[2 * third..]
        .iter()
        .take(query_seq.len())
        .map(|&b| b'!' + (b % 42))
        .collect();

    let mut config = Config::default();
    config.compat = decode_compat_profile(compat_mode, &data[4..12]);

    let mut seed = SeedMap::new();
    let seq_start = offset as u32;
    let seq_end = seq_start
        .saturating_add(ref_seq.len() as u32)
        .saturating_sub(1);
    seed.index_sequence(&ref_seq, seq_start, seq_end, seq_end);

    // Pass 1: Raw hostile / arbitrary lengths (including `0` and `i32::MAX as u32` ops and `=`/`X`).
    let ops = [b'M', b'I', b'D', b'N', b'S', b'H', b'P', b'=', b'X'];
    let mut cigar = String::new();
    for chunk in payload[2 * third..].chunks(2).take(12) {
        let len: u32 = match chunk[0] {
            0 => 0,
            255 => i32::MAX as u32,
            raw => u32::from(raw % 80),
        };
        let op = ops[(chunk.get(1).copied().unwrap_or(0) as usize) % ops.len()] as char;
        cigar.push_str(&len.to_string());
        cigar.push(op);
    }
    if !cigar.is_empty() {
        let modifier = CigarModifier::new(
            pos, cigar, query_seq, query_qual, &ref_seq, offset, &seed, 10, max_rl, &config,
        );
        if let Ok(modified) = modifier.modify_cigar() {
            let parsed =
                parse_cigar(&modified.cigar).expect("modify_cigar Ok must emit valid CIGAR syntax");
            assert!(
                parsed.iter().all(|&(len, _)| len > 0),
                "modify_cigar emitted non-positive op length in {}",
                modified.cigar
            );
        }
    }

    // Pass 2: Matched-query-length valid SAM CIGARs (`[M, I, D, N, S, H, P]` + `[=, X]` when `eqx_cigar_read_length`, `len >= 1`).
    let sam_ops: &[u8] = if config.compat.eqx_cigar_read_length {
        b"MIDNSHP=X"
    } else {
        b"MIDNSHP"
    };
    let mut valid_cigar = String::new();
    let mut qlen: usize = 0;
    // Only place `H` at outer ends (`SAMv1 §1.4.6`) while keeping `P` and `=X` interior-capable.
    let inner_ops: &[u8] = if config.compat.eqx_cigar_read_length {
        b"MIDNSP=X"
    } else {
        b"MIDNSP"
    };
    let chunks: Vec<&[u8]> = payload[2 * third..].chunks(2).take(12).collect();
    let n_chunks = chunks.len();
    for (c_idx, chunk) in chunks.into_iter().enumerate() {
        let len = usize::from(chunk[0] % 80).max(1);
        let pool = if c_idx == 0 || c_idx + 1 == n_chunks {
            sam_ops
        } else {
            inner_ops
        };
        let op_byte = pool[(chunk.get(1).copied().unwrap_or(0) as usize) % pool.len()];
        let consumes_query = matches!(op_byte, b'M' | b'I' | b'S' | b'=' | b'X');
        if consumes_query {
            if qlen.saturating_add(len) > 4096 {
                break;
            }
            qlen += len;
        }
        valid_cigar.push_str(&len.to_string());
        valid_cigar.push(op_byte as char);
    }
    if !valid_cigar.is_empty() && qlen >= 1 {
        let q_slice = &payload[third..2 * third];
        let qual_slice = &payload[2 * third..];
        let matched_seq: Vec<u8> = (0..qlen)
            .map(|i| bases[usize::from(q_slice[i % q_slice.len()]) % bases.len()])
            .collect();
        let matched_qual: Vec<u8> = (0..qlen)
            .map(|i| b'!' + (qual_slice[i % qual_slice.len()] % 42))
            .collect();
        let modifier = CigarModifier::new(
            pos,
            valid_cigar,
            matched_seq,
            matched_qual,
            &ref_seq,
            offset,
            &seed,
            10,
            max_rl,
            &config,
        );
        if let Ok(modified) = modifier.modify_cigar() {
            let parsed =
                parse_cigar(&modified.cigar).expect("modify_cigar Ok must emit valid CIGAR syntax");
            assert!(
                parsed.iter().all(|&(len, _)| len > 0),
                "modify_cigar emitted non-positive op length in {}",
                modified.cigar
            );
            if config.compat.cigar_strip_outer_hard_clips {
                assert!(
                    !matches!(parsed.first(), Some(&(_, CigarKind::HardClip)))
                        && !matches!(parsed.last(), Some(&(_, CigarKind::HardClip))),
                    "cigar_strip_outer_hard_clips left outer H in {}",
                    modified.cigar
                );
            }
            if config.compat.cigar_homopolymer_merge {
                assert_eq!(
                    query_consumed_len(&parsed),
                    modified.query_seq.len() as u64,
                    "cigar_homopolymer_merge modify_cigar must conserve query consumption: {}",
                    modified.cigar
                );
            }
        }
    }
});
