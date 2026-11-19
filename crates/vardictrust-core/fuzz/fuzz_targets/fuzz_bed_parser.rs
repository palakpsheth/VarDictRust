#![no_main]

use std::collections::HashMap;
use std::io::Cursor;

use libfuzzer_sys::fuzz_target;
use vardictrust_core::config::Config;
use vardictrust_core::regions::{parse_bed_from_reader, parse_region_flag};

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

fuzz_target!(|data: &[u8]| {
    if data.is_empty() || data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    let mut config = Config::default();
    config.zero_based = (data[0] & 1) != 0;
    if (data[0] & 2) != 0 {
        config.amplicon = Some("10:0.95".to_string());
    }
    config.extension = ((data[0] >> 2) as i32) - 16;

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
