#![no_main]

use std::collections::HashMap;
use std::io::Cursor;

use libfuzzer_sys::fuzz_target;
use vardictrust_core::config::Config;
use vardictrust_core::regions::{parse_bed_from_reader, parse_region_flag};

fuzz_target!(|data: &[u8]| {
    if data.is_empty() {
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

    let payload = &data[1..];
    let _ = parse_bed_from_reader(Cursor::new(payload), &config, &chromosomes_lengths, None);

    if let Ok(region_str) = std::str::from_utf8(payload) {
        for line in region_str.lines().take(8) {
            let _ = parse_region_flag(line, &config, &chromosomes_lengths);
        }
    }
});
