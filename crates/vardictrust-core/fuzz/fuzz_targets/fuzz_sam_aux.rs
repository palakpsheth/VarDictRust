#![no_main]

use std::io::Cursor;

use libfuzzer_sys::fuzz_target;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use vardictrust_core::bam::{RecordPreprocessor, SamRecord};
use vardictrust_core::config::{
    CompatProfile, Config, COMPAT_OFF, COMPAT_SWITCH_COUNT, COMPAT_VARDICTJAVA_1_8_3,
};
use vardictrust_core::types::Pos1;

/// Maximum accepted libFuzzer input payload (`8 KiB`) to bound per-iteration runtime.
const MAX_FUZZ_INPUT_BYTES: usize = 8192;

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

fn exercise_dup_key_switches(config: &Config, seed_byte: u8) {
    // Exercise Switch 9 (`eqx_cigar_read_length`) and Switch 22 (`cigar_strip_outer_hard_clips`)
    // on paired unmapped-mate `DupKey::UnmappedMateCigar` (`flags = 0x1 | 0x8`, `mate_pos >= 10`).
    let aln_start = u32::from(seed_byte).saturating_add(100);
    let hp_len_l = u32::from(seed_byte % 7).saturating_add(1);
    let hp_len_r = u32::from((seed_byte >> 3) % 7).saturating_add(1);
    let pad_mid = u32::from((seed_byte >> 5) % 4).saturating_add(1);

    let rec_hp_eqx = SamRecord {
        name: b"dup_hp_eqx".to_vec(),
        flags: 0x1 | 0x8,
        reference_name: Some("chr1".to_string()),
        alignment_start: Some(Pos1(aln_start)),
        alignment_end: Some(Pos1(aln_start + 29)),
        mapping_quality: Some(60),
        cigar: vec![
            (hp_len_l, CigarKind::HardClip),
            (20, CigarKind::SequenceMatch),
            (pad_mid, CigarKind::Pad),
            (10, CigarKind::SequenceMismatch),
            (hp_len_r, CigarKind::HardClip),
        ],
        sequence: vec![b'A'; 30],
        quality_scores: vec![b'I'; 30],
        mate_reference_name: Some("chr1".to_string()),
        mate_alignment_start: Some(Pos1(aln_start)),
        ..SamRecord::default()
    };
    let rec_plain = SamRecord {
        name: b"dup_plain".to_vec(),
        cigar: vec![(30, CigarKind::Match)],
        ..rec_hp_eqx.clone()
    };

    let mut dedup_cfg = config.clone();
    dedup_cfg.dedup = true;
    let mut pre = RecordPreprocessor::new();
    assert!(pre.preprocess_record(&rec_hp_eqx, &dedup_cfg));
    let kept_second = pre.preprocess_record(&rec_plain, &dedup_cfg);
    if dedup_cfg.compat.cigar_strip_outer_hard_clips && dedup_cfg.compat.eqx_cigar_read_length {
        assert!(
            !kept_second,
            "outer H + interior P + =/X must deduplicate against 30M when Switch 9 & Switch 22 are active"
        );
        assert_eq!(pre.duplicate_reads, 1);
    }
}

fuzz_target!(|data: &[u8]| {
    if data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    let (ctrl, mask64, sam_payload) = if data.len() >= 10 {
        let m = u64::from_le_bytes([
            data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8],
        ]);
        (data[0], m, &data[9..])
    } else {
        (data.first().copied().unwrap_or(0), 0u64, data)
    };

    let mut config = Config::default();
    config.compat = decode_compat(ctrl >> 4, mask64);
    config.dedup = (ctrl & 1) != 0;
    if (ctrl & 2) != 0 {
        config.flag_filter = 0;
        config.samfilter_raw = "0".to_string();
    }
    if (ctrl & 4) != 0 {
        config.map_qual = Some(u32::from(ctrl & 0x1f));
    }

    exercise_dup_key_switches(&config, ctrl);

    let header_src = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n";
    let mut sam_bytes = Vec::with_capacity(header_src.len() + sam_payload.len() + 1);
    sam_bytes.extend_from_slice(header_src);
    sam_bytes.extend_from_slice(sam_payload);
    sam_bytes.push(b'\n');

    let mut reader = noodles::sam::io::Reader::new(Cursor::new(sam_bytes));
    let Ok(header) = reader.read_header() else {
        return;
    };

    let mut preprocessor = RecordPreprocessor::new();
    let mut sam_rec = SamRecord::new();

    for result in reader.records().take(16) {
        if let Ok(record) = result {
            if sam_rec.fill_from_alignment_record(&record, &header).is_ok() {
                let _ = preprocessor.preprocess_record(&sam_rec, &config);
            }
        }
    }
    assert!(preprocessor.duplicate_reads <= preprocessor.total_reads);
});
