#![no_main]

use std::collections::HashMap;
use std::io::Cursor;

use libfuzzer_sys::fuzz_target;
use vardictrust_core::config::{
    java_integer_decode, parse_compat_target, parse_java_number_f64, parse_java_number_i32,
    parse_workflow_preset, validate_sample_regexp, CompatProfile, Config, VcfThreshold,
    WorkflowPreset, COMPAT_OFF, COMPAT_SWITCH_COUNT, COMPAT_VARDICTJAVA_1_8_3,
};
use vardictrust_core::regions::{parse_bed_from_reader, parse_region_flag};
use vardictrust_core::types::VardictError;

/// Maximum accepted libFuzzer input payload (`8 KiB`) to bound per-iteration runtime.
const MAX_FUZZ_INPUT_BYTES: usize = 8192;

/// Decoded custom BED column mapping (`-c/-S/-E/-g/-s/-e`).
#[derive(Clone, Copy)]
struct ColumnIndices {
    col_chr: usize,
    column_for_chromosome: Option<usize>,
    col_start: usize,
    col_end: usize,
    col_seg_start: usize,
    col_seg_end: usize,
    col_gene: Option<usize>,
}

impl ColumnIndices {
    fn from_selector(selector: u8, b1: u8, b2: u8) -> Self {
        let col_chr = usize::from(b1 & 0x0f);
        let column_for_chromosome = if (selector & 2) != 0 {
            Some(col_chr)
        } else {
            None
        };
        let col_start = usize::from((b1 >> 4) & 0x0f);
        let col_end = usize::from(b2 & 0x0f);
        let seg_mode = (selector >> 2) & 0x03;
        let (col_seg_start, col_seg_end) = match seg_mode {
            0 => (col_start, col_end),
            1 => (usize::from((b2 >> 4) & 0x0f), col_end),
            2 => (col_start, usize::from((b2 >> 4) & 0x0f)),
            _ => (
                usize::from((b2 >> 4) & 0x07),
                usize::from((b2 >> 5) & 0x07).saturating_add(1),
            ),
        };
        let gene_raw = (selector >> 4) & 0x0f;
        let col_gene = if gene_raw == 0 {
            None
        } else {
            Some(usize::from(gene_raw - 1))
        };
        Self {
            col_chr,
            column_for_chromosome,
            col_start,
            col_end,
            col_seg_start,
            col_seg_end,
            col_gene,
        }
    }

    fn apply_to(self, config: &mut Config) {
        config.col_chr = self.col_chr;
        config.column_for_chromosome = self.column_for_chromosome;
        config.col_start = self.col_start;
        config.col_end = self.col_end;
        config.col_seg_start = self.col_seg_start;
        config.col_seg_end = self.col_seg_end;
        config.col_gene = self.col_gene;
    }
}

fn fuzz_cli_and_preset_config(data: &[u8]) {
    let mut raw_mask = [0u8; 8];
    for (i, slot) in raw_mask.iter_mut().enumerate() {
        *slot = data.get(i % data.len()).copied().unwrap_or(0);
    }
    let mask64 = u64::from_le_bytes(raw_mask);
    let _ = CompatProfile::from_mask(mask64);
    let _ = CompatProfile::only_switch((mask64 as usize) % COMPAT_SWITCH_COUNT);
    let _ = CompatProfile::without_switch((mask64 as usize) % COMPAT_SWITCH_COUNT);

    // Exercise `--preset {wgs,wes,amplicon,ctdna}` defaults + explicit `-r` / `-q` floor relaxation.
    let presets = [
        WorkflowPreset::Wgs,
        WorkflowPreset::Wes,
        WorkflowPreset::Amplicon,
        WorkflowPreset::Ctdna,
    ];
    let preset = presets[usize::from(data[0]) % presets.len()];
    let mut preset_cfg = Config::from_preset(preset);
    assert_eq!(preset_cfg.compat, COMPAT_OFF);
    assert_eq!(preset_cfg.preset, Some(preset));

    let orig_vd = preset_cfg.vcf_min_vd.as_ref().map(|t| t.value);
    let explicit_min_reads = u32::from(data.get(1).copied().unwrap_or(1) % 8);
    let explicit_base_qual = f64::from(data.get(2).copied().unwrap_or(10) % 35) * 0.5;
    preset_cfg.min_reads = explicit_min_reads;
    preset_cfg.base_qual = explicit_base_qual;

    let explicit_vd = i32::try_from(explicit_min_reads).unwrap_or(i32::MAX);
    if let Some(cur_vd) = orig_vd {
        if explicit_vd < cur_vd {
            preset_cfg.vcf_min_vd = Some(VcfThreshold::new(explicit_vd));
        }
    }
    let cur_qmean = preset_cfg.vcf_qmean.as_ref().map_or(22.5, |t| t.value);
    let clamped_q = explicit_base_qual.max(0.0);
    if clamped_q < cur_qmean {
        preset_cfg.vcf_qmean = Some(VcfThreshold::new(clamped_q));
    }
    assert!(preset_cfg.validate_vcf_thresholds().is_ok());
    if let (Some(ref vd_thr), Some(cur_vd)) = (&preset_cfg.vcf_min_vd, orig_vd) {
        assert!(vd_thr.value >= 0 && vd_thr.value <= cur_vd && vd_thr.value <= explicit_vd);
    }
    if let Some(ref q_thr) = preset_cfg.vcf_qmean {
        assert!(q_thr.value >= 0.0 && q_thr.value <= cur_qmean);
    }

    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        for &(tok, expected) in &[
            ("wgs", Some(WorkflowPreset::Wgs)),
            ("wes", Some(WorkflowPreset::Wes)),
            ("amplicon", Some(WorkflowPreset::Amplicon)),
            ("ctdna", Some(WorkflowPreset::Ctdna)),
            ("hifi", None),
            ("ont", None),
        ] {
            match (parse_workflow_preset(tok), expected) {
                (Ok(got), Some(exp)) => assert_eq!(got, exp),
                (Err(VardictError::Unsupported(_)), None) => {}
                other => panic!("unexpected parse_workflow_preset({tok:?}): {other:?}"),
            }
        }
        assert_eq!(
            parse_compat_target("vardictjava").ok(),
            Some(COMPAT_VARDICTJAVA_1_8_3.target)
        );
        assert_eq!(parse_compat_target("off").ok(), Some(COMPAT_OFF.target));
    });

    if let Ok(tok) = std::str::from_utf8(&data[..data.len().min(96)]) {
        let first_line = tok.lines().next().unwrap_or(tok);
        let _ = parse_compat_target(first_line);
        let _ = parse_workflow_preset(first_line);
        let _ = java_integer_decode(first_line);
        let _ = parse_java_number_f64(first_line);
        let _ = parse_java_number_i32(first_line);
        let _ = validate_sample_regexp(first_line);
    }
}

fuzz_target!(|data: &[u8]| {
    if data.is_empty() || data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    fuzz_cli_and_preset_config(data);

    let mut config = Config::default();
    config.zero_based = (data[0] & 1) != 0;
    if (data[0] & 2) != 0 {
        config.amplicon = Some("10:0.95".to_string());
    }
    config.extension = ((data[0] >> 2) as i32) - 16;
    if data.len() >= 12 {
        let mask64 = u64::from_le_bytes([
            data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11],
        ]);
        config.compat = match data[0] >> 6 {
            0 => COMPAT_VARDICTJAVA_1_8_3,
            1 => COMPAT_OFF,
            _ => CompatProfile::from_mask(mask64),
        };
    }

    let mut chromosomes_lengths = HashMap::new();
    chromosomes_lengths.insert("chr1".to_string(), 250_000_000usize);
    chromosomes_lengths.insert("1".to_string(), 250_000_000usize);
    chromosomes_lengths.insert("HLA-A*01:01:01:01".to_string(), 3_503usize);
    chromosomes_lengths.insert("chrX".to_string(), 156_040_895usize);

    let (payload, zero_override) = if data.len() >= 5 && (data[1] & 1) != 0 {
        let cols = ColumnIndices::from_selector(data[1], data[2], data[3]);
        cols.apply_to(&mut config);
        config.zero_based_defined = (data[3] & 0x80) != 0;
        let z_ovr = match (data[1] >> 6) & 0x03 {
            0 => None,
            1 => Some(false),
            _ => Some(true),
        };
        (&data[4..], z_ovr)
    } else {
        (&data[1..], None)
    };

    let _ = parse_bed_from_reader(
        Cursor::new(payload),
        &config,
        &chromosomes_lengths,
        zero_override,
    );

    if let Ok(region_str) = std::str::from_utf8(payload) {
        for (idx, line) in region_str.lines().take(8).enumerate() {
            let _ = parse_region_flag(line, &config, &chromosomes_lengths);
            if idx == 0 && line.len() <= 128 {
                let hla_probe = format!("HLA-A*01:01:01:01:{line}");
                let _ = parse_region_flag(&hla_probe, &config, &chromosomes_lengths);
            }
        }
    }
});
