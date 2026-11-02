#![no_main]

use libfuzzer_sys::fuzz_target;
use vardictrust_core::realign::{find_35_match, find_bp};
use vardictrust_core::reference::{LoadedRegion, Reference};

fuzz_target!(|data: &[u8]| {
    if data.len() < 8 {
        return;
    }

    let dir = if data[0] & 1 == 0 { 1 } else { -1 };
    let is_len3 = (data[0] & 2) != 0;
    let indel_len = (data[1] as i32).clamp(1, 200);
    let pos = (u16::from_le_bytes([data[2], data[3]]) as i32).clamp(1, 1000);

    let rest = &data[4..];
    let mid = rest.len() / 2;
    let seq5 = &rest[..mid];
    let seq3 = &rest[mid..];

    let _ = find_35_match(seq5, seq3, is_len3);

    let bases = [b'A', b'C', b'G', b'T', b'N'];
    let ref_seq: Vec<u8> = rest
        .iter()
        .take(1024)
        .map(|&b| bases[(b as usize) % bases.len()])
        .collect();
    if !ref_seq.is_empty() {
        let chr_len = ref_seq.len() as u32;
        let mut reference = Reference::new();
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: chr_len,
        });
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);
        let _ = find_bp(seq5, pos, &reference, dir, "chr1", indel_len, chr_len);
    }
});
