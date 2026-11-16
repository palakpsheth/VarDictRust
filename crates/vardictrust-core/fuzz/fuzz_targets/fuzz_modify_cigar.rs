#![no_main]

use libfuzzer_sys::fuzz_target;
use vardictrust_core::cigar::modify::CigarModifier;
use vardictrust_core::config::Config;
use vardictrust_core::reference::SeedMap;

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
    let offset = (data[3] as i32).clamp(1, 64);

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

    let ops = [b'M', b'I', b'D', b'N', b'S', b'H', b'P', b'=', b'X'];
    let mut cigar = String::new();
    for chunk in payload[2 * third..].chunks(2).take(12) {
        let len = (chunk[0] % 80) as usize;
        let op = ops[(chunk.get(1).copied().unwrap_or(0) as usize) % ops.len()] as char;
        cigar.push_str(&len.to_string());
        cigar.push(op);
    }
    if cigar.is_empty() {
        return;
    }

    let config = Config::default();
    let mut seed = SeedMap::new();
    let seq_start = offset as u32;
    let seq_end = seq_start
        .saturating_add(ref_seq.len() as u32)
        .saturating_sub(1);
    seed.index_sequence(&ref_seq, seq_start, seq_end, seq_end);
    let modifier = CigarModifier::new(
        pos, cigar, query_seq, query_qual, &ref_seq, offset, &seed, 10, max_rl, &config,
    );
    let _ = modifier.modify_cigar();
});
