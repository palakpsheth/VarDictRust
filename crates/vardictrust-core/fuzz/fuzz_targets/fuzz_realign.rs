#![no_main]

use libfuzzer_sys::fuzz_target;
use vardictrust_core::config::{
    CompatProfile, Config, COMPAT_OFF, COMPAT_SWITCH_COUNT, COMPAT_VARDICTJAVA_1_8_3,
};
use vardictrust_core::realign::{
    adj_snv, find_35_match, find_35_match_compat, find_bp, realign_indels,
};
use vardictrust_core::reference::{LoadedRegion, Reference};
use vardictrust_core::regions::Region;
use vardictrust_core::types::Pos1;
use vardictrust_core::variations::{Sclip, Variation, VariationData, VariationMap};

/// 12-byte header + 2 * 258 bp max overhang.
const MAX_FUZZ_INPUT_BYTES: usize = 528;

fn decode_compat(mode: u8, mask64: u64) -> CompatProfile {
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

fn populate_sclip(sclip: &mut Sclip, overhang: &[u8], tie_probe: bool) {
    let bases = [b'A', b'C', b'G', b'T'];
    sclip.variation.cnt = 8;
    sclip.variation.fwd = 4;
    sclip.variation.rev = 4;
    sclip.variation.hicnt = 8;
    sclip.variation.pmean = 160.0;
    sclip.variation.qmean = 280.0;
    sclip.variation.mqmean = 480.0;
    sclip.variation.pstd = true;
    sclip.variation.qstd = true;
    for (idx, &raw) in overhang.iter().take(24).enumerate() {
        let pos_in_clip = idx as i32;
        let maj_base = bases[usize::from(raw) % 4];
        let min_base = bases[(usize::from(raw) + 1) % 4];
        let nt_entry = sclip.nt.entry(pos_in_clip).or_default();
        let seq_entry = sclip.seq.entry(pos_in_clip).or_default();
        if tie_probe && idx == 0 {
            // Exercise Switch 33 (`sclip_majority_count_first`): high-count lower-Q majority
            // (`cnt=6, qmean=26.0`) vs 1-read higher-Q error (`cnt=1, qmean=40.0`).
            nt_entry.insert(maj_base, 6);
            let mut maj_v = Variation::new();
            maj_v.cnt = 6;
            maj_v.qmean = 26.0;
            seq_entry.insert(maj_base, maj_v);

            nt_entry.insert(min_base, 1);
            let mut min_v = Variation::new();
            min_v.cnt = 1;
            min_v.qmean = 40.0;
            seq_entry.insert(min_base, min_v);
        } else {
            nt_entry.insert(maj_base, 7);
            let mut v = Variation::new();
            v.cnt = 7;
            v.qmean = 35.0;
            seq_entry.insert(maj_base, v);
        }
    }
}

fuzz_target!(|data: &[u8]| {
    if data.len() < 16 || data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    let mode_byte = data[0];
    let dir = if mode_byte & 1 == 0 { 1 } else { -1 };
    let is_len3 = (mode_byte & 2) != 0;
    let tie_probe = (mode_byte & 4) != 0;
    let use_adaptor = (mode_byte & 8) != 0;
    let indel_len = (data[1] as i32).clamp(1, 200);
    let pos = (u16::from_le_bytes([data[2], data[3]]) as i32).clamp(1, 512);
    let mask64 = u64::from_le_bytes([
        data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11],
    ]);

    let mut config = Config::default();
    config.compat = decode_compat(mode_byte >> 4, mask64);
    if use_adaptor {
        config.adaptor.push("AGATCGGAAGAGC".to_string());
    }

    let rest = &data[12..];
    let mid = rest.len() / 2;
    let seq5 = &rest[..mid];
    let seq3 = &rest[mid..];

    let _ = find_35_match(seq5, seq3, is_len3);
    let _ = find_35_match_compat(seq5, seq3, false, config.compat.exact_35_overlap);

    let bases = [b'A', b'C', b'G', b'T', b'N'];
    let mut ref_seq: Vec<u8> = rest
        .iter()
        .take(512)
        .map(|&b| bases[(b as usize) % bases.len()])
        .collect();
    if ref_seq.len() < 64 {
        ref_seq.resize(64, b'A');
    }
    let chr_len = ref_seq.len() as u32;
    let mut reference = Reference::new();
    reference.loaded_regions.push(LoadedRegion {
        chr: "chr1".to_string(),
        sequence_start: 1,
        sequence_end: chr_len,
    });
    reference.region_sequences.push(ref_seq.clone());
    reference.seed.index_sequence(&ref_seq, 1, chr_len, chr_len);
    reference
        .reference_sequences
        .insert("chr1".to_string(), ref_seq.clone());

    let _ = find_bp(seq5, pos, &reference, dir, "chr1", indel_len, chr_len);

    // Exercise `realign_indels` and `adj_snv` under `Switch 33` (`sclip_majority_count_first`),
    // `Switch 37` (`adj_snv_full_sclip_match`), and `Switch 19` (`adaptor_all_consensus_dirs`).
    let anchor = (pos as u32).clamp(4, chr_len.saturating_sub(4).max(4)) as i32;
    let region = Region {
        chr: "chr1".to_string(),
        start: Pos1(1),
        end: Pos1(chr_len),
        gene: "FUZZ_REALIGN".to_string(),
        amplicon: None,
    };
    let mut vd = VariationData::new();
    vd.max_read_length = 100;
    vd.ref_coverage.insert(anchor, 20);

    // Seed SNV candidate at `anchor` + 1-bp / multi-bp 3' soft-clip at `anchor + 1` for `adj_snv` (`Switch 37`).
    let ref_b = ref_seq.get((anchor - 1) as usize).copied().unwrap_or(b'A');
    let alt_b = if ref_b == b'A' { b'G' } else { b'A' };
    let mut vmap = VariationMap::new();
    let mut ref_var = Variation::new();
    ref_var.cnt = 12;
    ref_var.fwd = 6;
    ref_var.rev = 6;
    ref_var.hicnt = 12;
    ref_var.pmean = 240.0;
    ref_var.qmean = 420.0;
    ref_var.mqmean = 720.0;
    vmap.map.insert((ref_b as char).to_string(), ref_var);

    let mut snv_var = Variation::new();
    snv_var.cnt = 4;
    snv_var.fwd = 2;
    snv_var.rev = 2;
    snv_var.hicnt = 4;
    snv_var.pmean = 40.0;
    snv_var.qmean = 140.0;
    snv_var.mqmean = 240.0;
    snv_var.pstd = true;
    snv_var.qstd = true;
    vmap.map.insert((alt_b as char).to_string(), snv_var);
    vd.non_insertion_variants.insert(anchor, vmap);

    // 3' soft-clip at `anchor`: 1-bp vs multi-bp overhang starting with `alt_b`
    let sclip_len = if (mode_byte & 2) != 0 {
        1usize
    } else {
        seq3.len().clamp(2, 16)
    };
    let mut clip3_bytes = Vec::with_capacity(sclip_len);
    clip3_bytes.push(alt_b);
    for i in 1..sclip_len {
        let b = if (mode_byte & 1) != 0 {
            ref_seq
                .get((anchor as usize).saturating_add(i).saturating_sub(1))
                .copied()
                .unwrap_or(b'C')
        } else {
            bases[usize::from(seq3.get(i).copied().unwrap_or(0)) % 4]
        };
        clip3_bytes.push(b);
    }
    let mut sc3 = Sclip::new();
    populate_sclip(&mut sc3, &clip3_bytes, tie_probe);
    vd.soft_clips_3_end.insert(anchor, sc3);

    if !seq5.is_empty() {
        let mut sc5 = Sclip::new();
        populate_sclip(&mut sc5, seq5, tie_probe);
        vd.soft_clips_5_end
            .insert(anchor.saturating_add(6).min(chr_len as i32), sc5);
    }

    realign_indels(
        None,
        &mut reference,
        "chr1",
        chr_len,
        &region,
        &config,
        &mut vd,
    );
    adj_snv(&reference, "chr1", &config, &mut vd);
});
