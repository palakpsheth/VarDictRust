//! Fuzz target: BED text + `chr:start-end` region specs through the
//! vardict-core parsers must never panic; tiled intervals stay in bounds.
#![no_main]

use libfuzzer_sys::fuzz_target;
use std::collections::HashMap;
use vardict_core::config::Config;
use vardict_core::regions;

fn bed_config() -> Config {
    let mut config = Config::default();
    config.delimiter = "\t".to_string();
    config.col_chr = 0;
    config.col_start = 1;
    config.col_end = 2;
    config.col_seg_start = 1;
    config.col_seg_end = 2;
    config.col_gene = Some(3);
    config.zero_based = true;
    config
}

fuzz_target!(|data: &[u8]| {
    let Ok(text) = std::str::from_utf8(data) else {
        return;
    };
    let chr_lengths: HashMap<String, usize> = HashMap::from([
        ("chr1".to_string(), 100_000usize),
        ("chr2".to_string(), 100_000),
    ]);

    // Non-amplicon BED path.
    if let Ok(segs) =
        regions::parse_bed_from_reader(text.as_bytes(), &bed_config(), &chr_lengths, None)
    {
        // Tiling is only meaningful for regions within known contig bounds:
        // the production pipeline validates each region against the reference
        // (fetch_subsequence -> noodles "invalid start bound") before tiling,
        // so absurd spans never reach tile() there. Mirror that domain here;
        // otherwise one 2-billion-bp parsed span expands into ~33M Tiles and
        // OOMs the fuzzer (found during triage 2026-08-24; not an engine bug).
        let max_contig = chr_lengths.values().copied().max().unwrap_or(0);
        for group in segs.iter().flatten() {
            let span = group.end.0.saturating_sub(group.start.0);
            if u64::from(span) <= max_contig as u64 {
                for t in regions::tile(group, 64, 16) {
                    assert!(t.halo_start <= t.start && t.start <= t.end && t.end <= t.halo_end);
                    assert!(t.start <= group.end && t.end <= group.end);
                }
            }
        }
    }

    // Amplicon BED path (grouping/sort logic).
    let mut amp = bed_config();
    amp.amplicon = Some("10:0.95".to_string());
    let _ = regions::parse_bed_from_reader(text.as_bytes(), &amp, &chr_lengths, None);

    // `-R`-style region spec path (first line of the input).
    if let Some(spec) = text.lines().next() {
        let _ = regions::parse_region_flag(spec.trim(), &bed_config(), &chr_lengths);
    }
});
