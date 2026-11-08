//! VarDictRust core engine library.
//!
//! An idiomatic Rust rewrite of VarDictJava v1.8.3 — a sensitive somatic and
//! germline variant caller for next-generation sequencing data.
//!
//! This crate contains all variant-calling logic. The binary crate
//! (`vardictrust`) provides the command-line interface.
//!
//! # Safety
//! This crate forbids unsafe code: `#![forbid(unsafe_code)]`.
//!
//! # Reference
//! Validated against VarDictJava v1.8.3 (commit `0886261`).

#![forbid(unsafe_code)]
#![cfg_attr(
    not(test),
    deny(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::unreachable,
        clippy::todo,
        clippy::string_slice,
    )
)]
#![allow(
    clippy::too_many_arguments,
    clippy::unwrap_or_default,
    clippy::needless_lifetimes,
    clippy::assign_op_pattern,
    clippy::field_reassign_with_default,
    clippy::int_plus_one,
    clippy::unnecessary_cast,
    clippy::let_and_return,
    clippy::len_zero,
    clippy::collapsible_if
)]

pub mod bam;
pub mod cigar;
pub mod config;
pub mod counters;
pub mod javamap;
pub mod output;
pub mod postprocess;
pub mod realign;
pub mod reference;
pub mod regions;
pub mod runner;
pub mod stats;
pub mod sv;
pub mod tovars;
pub mod types;
pub mod util;
pub mod variations;

use std::collections::HashMap;
use std::path::Path;
use types::VardictError;

/// Default region tile size (`100 kbp`) for bounded-memory WGS tiling of multi-megabase intervals.
pub const DEFAULT_TILE_SIZE: u32 = 100_000;

/// Maximum untiled BED interval span (`1 Mbp`): when `--tile-size` is not explicitly specified,
/// regions up to `1 Mbp` execute as a single untiled interval like VarDictJava, while regions
/// exceeding `1 Mbp` (such as `-R chr1:1-10000000` or whole-chromosome sweeps) are tiled into
/// `DEFAULT_TILE_SIZE` (`100 kbp`) windows to preserve `<1.5 GiB` bounded WGS memory.
pub const DEFAULT_UNTILED_SPAN_MAX: u32 = 1_000_000;

/// Maximum `--tile-size` (`500 Mbp`, above human `chr1` `~249 Mbp`) before `u32` wrap guard.
pub const MAX_TILE_SIZE: u32 = 500_000_000;

const MAX_SPLIT_VAR_DEDUP_ENTRIES: usize = 262_144;

#[inline]
fn split_var_line_fingerprint(line: &str) -> (u64, u64) {
    use std::hash::{Hash, Hasher};
    let mut h1 = rustc_hash::FxHasher::default();
    0x9e37_79b9_7f4a_7c15u64.hash(&mut h1);
    line.as_bytes().hash(&mut h1);
    let mut h2 = rustc_hash::FxHasher::default();
    0xc2b2_ae3d_27d4_eb4fu64.hash(&mut h2);
    line.len().hash(&mut h2);
    line.as_bytes().hash(&mut h2);
    (h1.finish(), h2.finish())
}

#[inline]
fn effective_tile_size(reg: &crate::regions::Region, config_tile_size: Option<u32>) -> u32 {
    if let Some(ts) = config_tile_size {
        ts
    } else {
        let span = reg.end.0.saturating_sub(reg.start.0).saturating_add(1);
        if span <= DEFAULT_UNTILED_SPAN_MAX {
            span.max(1)
        } else {
            DEFAULT_TILE_SIZE
        }
    }
}

/// Post-pileup chain selection for [`run_pipeline_on_records`], mirroring the
/// mode-specific pipelines in VarDictJava's SimpleMode/SomaticMode.
enum PipelineTail {
    /// Splice mode: stop after the pileup (SimpleMode.java returns early).
    Stop,
    /// Full post-pileup pipeline (AbstractMode.java:44-56, VariationRealigner.java:73-90):
    /// adjust_mnp unconditionally, realign_indels under `-k 1`, then find_all_sv + adj_snv.
    Full,
}

/// Splits a colon-delimited BAM specification into non-empty file paths,
/// returning `Err(VardictError::Parse(empty_err))` when no path is present.
fn parse_colon_bams(spec: &str, empty_err: &str) -> Result<Vec<String>, VardictError> {
    let bams: Vec<String> = spec
        .split(':')
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();
    if bams.is_empty() {
        return Err(VardictError::Parse(empty_err.to_string()));
    }
    Ok(bams)
}

/// Per-region pipeline shared by simple, amplicon and somatic modes: query the
/// BAM(s), preprocess the records, optionally load the reference, pile up,
/// then run the post-pileup chain selected by `tail`.
///
/// Multiple `bam_paths` are queried in order and merged before a single
/// preprocessing pass (the somatic combineAnalysis shape).
#[allow(clippy::too_many_arguments)]
fn run_pipeline_on_records(
    bam_paths: &[String],
    chr: &str,
    query_start: crate::types::Pos0,
    query_end: crate::types::Pos0,
    ref_region: Option<&crate::regions::Region>,
    reference: &mut crate::reference::Reference,
    pileup_region: &crate::regions::Region,
    tail: PipelineTail,
    config: &config::Config,
    chromosomes_lengths: &HashMap<String, usize>,
    seed_max_read_length: i32,
    initial_splice: Option<std::collections::HashSet<String>>,
) -> Result<crate::variations::VariationData, VardictError> {
    run_pipeline_on_records_ext(
        bam_paths,
        chr,
        query_start,
        query_end,
        ref_region,
        reference,
        pileup_region,
        tail,
        config,
        chromosomes_lengths,
        seed_max_read_length,
        initial_splice,
        &[],
        None,
        None,
    )
}

#[derive(Clone, Copy, Debug)]
struct MaxRlStep {
    bam_idx: usize,
    read_start: u32,
    same_start_overlap_before: usize,
    new_max: i32,
}

#[inline]
fn record_read_end(r: &crate::bam::SamRecord, read_start: i32) -> i32 {
    let mut read_end = read_start;
    use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
    for &(len, op) in &r.cigar {
        if matches!(
            op,
            CigarKind::Match
                | CigarKind::Deletion
                | CigarKind::Skip
                | CigarKind::SequenceMatch
                | CigarKind::SequenceMismatch
        ) {
            read_end += len as i32;
        }
    }
    read_end
}

#[inline]
fn record_overlaps_region(r: &crate::bam::SamRecord, reg_start: u32, reg_end: u32) -> bool {
    let read_start = r.alignment_start.map(|p| p.0 as i32).unwrap_or(0);
    let read_end = record_read_end(r, read_start);
    read_start <= reg_end as i32 && read_end > reg_start as i32
}

#[allow(clippy::too_many_arguments)]
fn run_pipeline_on_records_ext(
    bam_paths: &[String],
    chr: &str,
    query_start: crate::types::Pos0,
    query_end: crate::types::Pos0,
    ref_region: Option<&crate::regions::Region>,
    reference: &mut crate::reference::Reference,
    pileup_region: &crate::regions::Region,
    tail: PipelineTail,
    config: &config::Config,
    chromosomes_lengths: &HashMap<String, usize>,
    seed_max_read_length: i32,
    initial_splice: Option<std::collections::HashSet<String>>,
    max_rl_steps: &[MaxRlStep],
    post_pileup_max_read_length: Option<i32>,
    parent_original_start: Option<crate::types::Pos1>,
) -> Result<crate::variations::VariationData, VardictError> {
    // Task 4 (streaming feed): one lazy pass per BAM, BAM1-then-BAM2, through
    // a single sequential preprocessor into a single CigarParser. Records are
    // dropped promptly after parsing, so peak memory stays flat in depth —
    // the old `raw_records` Vec held the whole tile. Note: `CigarParser.java:78`
    // calls `initFromScope(scope)` BEFORE the `processor.nextRecord()` loop, so
    // `totalReads` and `duplicateReads` remain 0 at finalization (G2 parity).
    // `?` propagates I/O errors with the same Err semantics as the Vec path;
    // nothing is skipped or unwrapped.
    let mut feeds = Vec::with_capacity(bam_paths.len());
    for path in bam_paths {
        feeds.push(crate::bam::query_region_iter(
            Path::new(path),
            chr,
            query_start,
            query_end,
            config,
        )?);
    }

    if let Some(ref_region) = ref_region {
        crate::reference::get_reference(
            ref_region,
            config.max_region_len,
            chromosomes_lengths,
            config,
            reference,
        )?;
    }

    let mut preprocessor = crate::bam::RecordPreprocessor::new();
    let mut parser = crate::cigar::parse::CigarParser::new(
        reference,
        pileup_region,
        config,
        false,
        0,
        0,
        seed_max_read_length,
    );
    if let Some(seed_splice) = initial_splice {
        parser.splice = seed_splice;
    }
    let mut step_cursor = 0usize;
    for (idx, mut feed) in feeds.into_iter().enumerate() {
        preprocessor.reset_for_next_bam();
        while let Some(step) = max_rl_steps.get(step_cursor) {
            if step.bam_idx < idx {
                parser.max_read_length = parser.max_read_length.max(step.new_max);
                step_cursor += 1;
            } else {
                break;
            }
        }
        // Reused record buffer: zero per-read allocation after warmup (the
        // retained-churn slope that OOMs deep tiles).
        let mut r = crate::bam::SamRecord::new();
        let mut any_passed = false;
        let mut cur_pos: Option<u32> = None;
        let mut same_start_ov_idx = 0usize;

        let mut handle_record = |r: &crate::bam::SamRecord,
                                 preprocessor: &mut crate::bam::RecordPreprocessor,
                                 parser: &mut crate::cigar::parse::CigarParser<'_>,
                                 any_passed: &mut bool|
         -> Result<(), VardictError> {
            if !max_rl_steps.is_empty()
                && record_overlaps_region(r, pileup_region.start.0, pileup_region.end.0)
            {
                let r_start = r.alignment_start.map(|p| p.0).unwrap_or(0);
                if cur_pos != Some(r_start) {
                    cur_pos = Some(r_start);
                    same_start_ov_idx = 0;
                }
                while let Some(step) = max_rl_steps.get(step_cursor) {
                    if (
                        step.bam_idx,
                        step.read_start,
                        step.same_start_overlap_before,
                    ) <= (idx, r_start, same_start_ov_idx)
                    {
                        parser.max_read_length = parser.max_read_length.max(step.new_max);
                        step_cursor += 1;
                    } else {
                        break;
                    }
                }
                same_start_ov_idx += 1;
            }
            if preprocessor.preprocess_record(r, config) {
                *any_passed = true;
                parser.parse_cigar(r)?;
            }
            Ok(())
        };

        if feed.next_into(&mut r)? {
            // Issue F Part C (`> 1 Mbp` tile boundary `-t` / `dedup` continuity):
            // `RecordPreprocessor.duplicates` is scoped to the current `alignment_start`
            // (`first_matching_position`). When a non-initial tile starts at `pileup_region.start > orig_start`
            // and its first overlapping read starts at `first_pos < pileup_region.start`, an earlier read
            // at the same `first_pos` whose alignment ended `< pileup_region.start` would have been seen
            // by the untiled parent-region query and seeded `duplicates` at `first_pos`. Re-opening from
            // `first_pos.max(orig_start)` visits the complete same-start group at `first_pos` (and any
            // read ending `<= pileup_region.start` returns immediately from `parse_cigar`).
            let mut switched_feed = false;
            if config.dedup {
                if let (Some(orig_start), Some(first_pos), Some(bam_path)) =
                    (parent_original_start, r.alignment_start, bam_paths.get(idx))
                {
                    if pileup_region.start.0 > orig_start.0 && first_pos.0 < pileup_region.start.0 {
                        let ext_start_1 = first_pos.0.max(orig_start.0);
                        let ext_query_start = crate::types::Pos1(ext_start_1).to_zero_based();
                        let mut ext_feed = crate::bam::query_region_iter(
                            Path::new(bam_path),
                            chr,
                            ext_query_start,
                            query_end,
                            config,
                        )?;
                        while ext_feed.next_into(&mut r)? {
                            if r.alignment_start.is_some_and(|p| p.0 < first_pos.0) {
                                continue;
                            }
                            handle_record(&r, &mut preprocessor, &mut parser, &mut any_passed)?;
                        }
                        switched_feed = true;
                    }
                }
            }
            if !switched_feed {
                handle_record(&r, &mut preprocessor, &mut parser, &mut any_passed)?;
                while feed.next_into(&mut r)? {
                    handle_record(&r, &mut preprocessor, &mut parser, &mut any_passed)?;
                }
            }
        }
        // Ports RecordPreprocessor.java:57-64 (`if (record == null)` advances at most one reader).
        if idx > 0 && !any_passed {
            break;
        }
    }
    let mut var_data = parser.into_variation_data();
    if let Some(post_max) = post_pileup_max_read_length {
        var_data.max_read_length = var_data.max_read_length.max(post_max);
    }
    let chr_len = *chromosomes_lengths.get(chr).unwrap_or(&0) as u32;
    let joined_bams = bam_paths.join(":");

    match tail {
        PipelineTail::Stop => {}
        PipelineTail::Full => {
            crate::realign::adjust_mnp(
                reference,
                chr,
                chr_len,
                pileup_region,
                config,
                &mut var_data,
            );
            if config.realign == 1 {
                crate::realign::realign_indels(
                    Some(bam_paths),
                    reference,
                    chr,
                    chr_len,
                    pileup_region,
                    config,
                    &mut var_data,
                );
            }
            let current_bam = std::path::Path::new(joined_bams.as_str());
            crate::sv::find_all_sv(
                &mut var_data,
                current_bam,
                reference,
                chr,
                chr_len,
                pileup_region,
                config,
                chromosomes_lengths,
            )?;
            crate::realign::adj_snv(reference, chr, config, &mut var_data);
        }
    }
    Ok(var_data)
}

type MultiTileKey = (String, u32, u32, u32, u32);

#[inline]
fn multi_tile_key(t: &crate::regions::Tile) -> MultiTileKey {
    (
        t.chr.clone(),
        t.original_start.0,
        t.original_end.0,
        t.start.0,
        t.end.0,
    )
}

/// Issue F Part C: For a multi-tile parent region (`reg_tiles.len() > 1`), scan the parent
/// region in BAM order to compute:
/// 1. `steps_per_tile[i]`: the exact BAM-ordered `MaxRlStep` trajectory for `reg_tiles[i]`, so
///    every overlapping read on `reg_tiles[i]` observes the exact same `max_read_length` as in an
///    untiled parent-region run.
/// 2. `final_max_rl`: the final `max_read_length` at the end of the parent region's pileup, used by
///    the post-pileup tail (`adjust_mnp`, `realign_indels`, `find_all_sv`, `adj_snv`) and to seed
///    BAM2 in Somatic mode (`SomaticMode.java:104-122`).
fn scan_multi_tile_max_read_lengths(
    bam_paths: &[String],
    reg: &crate::regions::Region,
    reg_tiles: &[crate::regions::Tile],
    config: &config::Config,
    chromosomes_lengths: &HashMap<String, usize>,
    initial_max_rl: i32,
) -> Result<(Vec<Vec<MaxRlStep>>, i32), VardictError> {
    let mut steps_per_tile = vec![Vec::new(); reg_tiles.len()];
    if reg_tiles.is_empty() {
        return Ok((steps_per_tile, initial_max_rl));
    }
    let mut cur_max = initial_max_rl;
    let query_start = reg.start.to_zero_based();
    let query_end = crate::types::Pos0(reg.end.0);

    let mut feeds = Vec::with_capacity(bam_paths.len());
    for path in bam_paths {
        feeds.push(crate::bam::query_region_iter(
            Path::new(path),
            &reg.chr,
            query_start,
            query_end,
            config,
        )?);
    }

    #[derive(Clone, Copy, Debug)]
    struct LoadedIntervalMeta {
        sequence_start: u32,
        sequence_end: u32,
        stored_end: u32,
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Seed12Multiplicity {
        Absent,
        Unique(u32),
        Multiple,
    }

    const REF_SCAN_WINDOW_BP: u32 = 128 * 1024;
    const SEED_SCAN_CHUNK_BP: u32 = 64 * 1024;
    const MAX_SEED12_MEMO_ENTRIES: usize = 65_536;

    /// Bounded-memory reference + 12-mer seed scanner for `scan_multi_tile_max_read_lengths`.
    /// Preserves the exact unpadded parent-region seed scope of an untiled run while keeping
    /// at most a `128 KiB` sliding reference window and `64 KiB` seed-scan chunk in RAM (`AGENTS.md §10`).
    struct BoundedParentRefScanner {
        chr: String,
        chr_len: u32,
        total_ext: u32,
        reg_start: u32,
        reg_end: u32,
        fasta_path: std::path::PathBuf,
        loaded_intervals: Vec<LoadedIntervalMeta>,
        window_interval_idx: Option<usize>,
        window_start: u32,
        window_end: u32,
        window_bases: Vec<u8>,
        seed12_memo: rustc_hash::FxHashMap<[u8; 12], (usize, Seed12Multiplicity)>,
    }

    impl BoundedParentRefScanner {
        fn new(
            reg: &crate::regions::Region,
            chromosomes_lengths: &HashMap<String, usize>,
            config: &config::Config,
        ) -> Self {
            let chr_len = chromosomes_lengths.get(&reg.chr).copied().unwrap_or(0) as u32;
            let total_ext = (config.extension.max(0) as u32).saturating_add(config.max_region_len);
            Self {
                chr: reg.chr.clone(),
                chr_len,
                total_ext,
                reg_start: reg.start.0,
                reg_end: reg.end.0,
                fasta_path: std::path::PathBuf::from(&config.reference),
                loaded_intervals: Vec::new(),
                window_interval_idx: None,
                window_start: 0,
                window_end: 0,
                window_bases: Vec::new(),
                seed12_memo: rustc_hash::FxHashMap::default(),
            }
        }

        fn ensure_interval_loaded(&mut self, start_val: u32, end_val: u32) {
            let sequence_start = start_val.saturating_sub(self.total_ext).max(1);
            let end_ext = end_val.saturating_add(self.total_ext);
            let sequence_end = if self.chr_len > 0 && end_ext > self.chr_len {
                self.chr_len
            } else {
                end_ext
            };
            if self
                .loaded_intervals
                .iter()
                .any(|iv| sequence_start >= iv.sequence_start && sequence_end <= iv.sequence_end)
            {
                return;
            }
            let stored_end = if self.chr_len == sequence_end {
                sequence_end
            } else if sequence_end.saturating_sub(sequence_start) >= 17 {
                sequence_end - 17
            } else {
                sequence_start.saturating_sub(1)
            };
            self.loaded_intervals.push(LoadedIntervalMeta {
                sequence_start,
                sequence_end,
                stored_end,
            });
        }

        fn slice_for_read(
            &mut self,
            read_start: i32,
            span_end: i32,
            ext_start_i32: i32,
            ext_end_i32: i32,
        ) -> Result<(&[u8], i32), VardictError> {
            if self.loaded_intervals.is_empty() {
                self.ensure_interval_loaded(self.reg_start, self.reg_end);
            }
            let ext_start = ext_start_i32.max(1) as u32;
            let ext_end = ext_end_i32.max(1) as u32;
            self.ensure_interval_loaded(ext_start, ext_end);

            let candidates = || {
                self.loaded_intervals.iter().enumerate().map(|(idx, iv)| {
                    let lr_start = iv.sequence_start as i32;
                    let lr_end = iv.stored_end.saturating_add(1) as i32;
                    (idx, *iv, lr_start, lr_end)
                })
            };
            let chosen = candidates()
                .find(|&(_, _, lr_start, lr_end)| {
                    ext_start_i32 >= lr_start && ext_end_i32 <= lr_end
                })
                .or_else(|| {
                    candidates().find(|&(_, _, lr_start, lr_end)| {
                        read_start >= lr_start && span_end <= lr_end
                    })
                })
                .or_else(|| {
                    candidates().find(|&(_, _, lr_start, lr_end)| {
                        span_end > lr_start && read_start < lr_end
                    })
                })
                .or_else(|| candidates().next());

            let Some((iv_idx, iv, _, _)) = chosen else {
                return Ok((&[], 1));
            };
            if iv.stored_end < iv.sequence_start {
                return Ok((&[], iv.sequence_start as i32));
            }

            let req_start = (ext_start_i32.saturating_sub(self.total_ext as i32))
                .max(iv.sequence_start as i32) as u32;
            let req_end = (ext_end_i32.saturating_add(self.total_ext as i32))
                .min(iv.stored_end as i32) as u32;
            let req_end = req_end.max(req_start).min(iv.stored_end);

            if self.window_interval_idx != Some(iv_idx)
                || req_start < self.window_start
                || req_end > self.window_end
            {
                let w_start = req_start.saturating_sub(4096).max(iv.sequence_start);
                let w_end = w_start
                    .saturating_add(REF_SCAN_WINDOW_BP - 1)
                    .max(req_end.saturating_add(4096))
                    .min(iv.stored_end);
                let mut bases = crate::reference::fetch_subsequence(
                    &self.fasta_path,
                    &self.chr,
                    crate::types::Pos1(w_start),
                    crate::types::Pos1(w_end),
                )?;
                for b in &mut bases {
                    b.make_ascii_uppercase();
                }
                self.window_interval_idx = Some(iv_idx);
                self.window_start = w_start;
                self.window_end = w_end;
                self.window_bases = bases;
            }

            Ok((&self.window_bases, self.window_start as i32))
        }

        fn unique_seed12_pos(
            &mut self,
            rc_seed: &[u8],
            read_start: i32,
            span_end: i32,
            cur_max: i32,
        ) -> Result<Option<u32>, VardictError> {
            let Ok(seed12) = <[u8; 12]>::try_from(rc_seed) else {
                return Ok(None);
            };
            let n_intervals = self.loaded_intervals.len();
            if let Some(&(cached_n, mult)) = self.seed12_memo.get(&seed12) {
                if cached_n == n_intervals || mult == Seed12Multiplicity::Multiple {
                    return Ok(match mult {
                        Seed12Multiplicity::Unique(p) => Some(p),
                        Seed12Multiplicity::Absent | Seed12Multiplicity::Multiple => None,
                    });
                }
            }

            // Fast local necessary-condition check: CigarModifier only uses a unique seed position
            // when `(this.position - ref_pos).abs() < 2 * max_read_length`, where `this.position`
            // lies in `[read_start, span_end]`. If the current window covers that entire local
            // band and `seed12` occurs 0 times (or >= 2 times) inside the window, a parent-wide
            // scan can be skipped (or immediately memoized as Multiple).
            let band_lo = read_start.saturating_sub(2 * cur_max).max(1) as u32;
            let band_hi = span_end.saturating_add(2 * cur_max).max(1) as u32;
            if self.window_interval_idx.is_some()
                && band_lo >= self.window_start
                && band_hi.saturating_add(11) <= self.window_end
            {
                let mut local_hits = 0u32;
                for w in self.window_bases.windows(12) {
                    if w == seed12 {
                        local_hits += 1;
                        if local_hits >= 2 {
                            break;
                        }
                    }
                }
                if local_hits == 0 {
                    return Ok(None);
                }
                if local_hits >= 2 {
                    if self.seed12_memo.len() >= MAX_SEED12_MEMO_ENTRIES {
                        self.seed12_memo.clear();
                    }
                    self.seed12_memo
                        .insert(seed12, (n_intervals, Seed12Multiplicity::Multiple));
                    return Ok(None);
                }
            }

            let mut mult = Seed12Multiplicity::Absent;
            'intervals: for (iv_idx, &iv) in self.loaded_intervals.iter().enumerate() {
                let last_indexed_pos = if self.chr_len == iv.sequence_end {
                    iv.sequence_end.saturating_sub(16)
                } else {
                    iv.stored_end
                };
                if last_indexed_pos < iv.sequence_start {
                    continue;
                }
                let mut chunk_start = iv.sequence_start;
                while chunk_start <= last_indexed_pos {
                    let chunk_last_pos = chunk_start
                        .saturating_add(SEED_SCAN_CHUNK_BP - 1)
                        .min(last_indexed_pos);
                    let fetch_end = chunk_last_pos.saturating_add(11).min(iv.sequence_end);
                    let mut chunk_bases = crate::reference::fetch_subsequence(
                        &self.fasta_path,
                        &self.chr,
                        crate::types::Pos1(chunk_start),
                        crate::types::Pos1(fetch_end),
                    )?;
                    for b in &mut chunk_bases {
                        b.make_ascii_uppercase();
                    }
                    let n_pos = (chunk_last_pos - chunk_start + 1) as usize;
                    for i in 0..n_pos {
                        let Some(kmer) = chunk_bases.get(i..i + 12) else {
                            break;
                        };
                        if kmer != seed12 {
                            continue;
                        }
                        let pos = chunk_start + i as u32;
                        let already_covered = self.loaded_intervals[..iv_idx]
                            .iter()
                            .any(|prev| pos >= prev.sequence_start && pos <= prev.stored_end);
                        if already_covered {
                            continue;
                        }
                        mult = match mult {
                            Seed12Multiplicity::Absent => Seed12Multiplicity::Unique(pos),
                            Seed12Multiplicity::Unique(_) | Seed12Multiplicity::Multiple => {
                                Seed12Multiplicity::Multiple
                            }
                        };
                        if mult == Seed12Multiplicity::Multiple {
                            break 'intervals;
                        }
                    }
                    if chunk_last_pos == last_indexed_pos {
                        break;
                    }
                    chunk_start = chunk_last_pos + 1;
                }
            }

            if self.seed12_memo.len() >= MAX_SEED12_MEMO_ENTRIES {
                self.seed12_memo.clear();
            }
            self.seed12_memo.insert(seed12, (n_intervals, mult));
            Ok(match mult {
                Seed12Multiplicity::Unique(p) => Some(p),
                Seed12Multiplicity::Absent | Seed12Multiplicity::Multiple => None,
            })
        }
    }

    let mut preprocessor = crate::bam::RecordPreprocessor::new();
    let mut cigar_elements = Vec::new();
    let mut cigar_str = String::new();
    let mut scanner = BoundedParentRefScanner::new(reg, chromosomes_lengths, config);

    for (idx, mut feed) in feeds.into_iter().enumerate() {
        preprocessor.reset_for_next_bam();
        let mut r = crate::bam::SamRecord::new();
        let mut any_passed = false;
        let mut cur_pos: Option<u32> = None;
        let mut same_start_ends: Vec<i32> = Vec::new();
        while feed.next_into(&mut r)? {
            let r_start_u32 = r.alignment_start.map(|p| p.0).unwrap_or(0);
            let read_start = r_start_u32 as i32;
            let read_end = record_read_end(&r, read_start);
            if cur_pos != Some(r_start_u32) {
                cur_pos = Some(r_start_u32);
                same_start_ends.clear();
            }

            let mut preprocessed = false;
            if config.dedup || (idx > 0 && !any_passed) {
                if !preprocessor.preprocess_record(&r, config) {
                    if read_start <= reg.end.0 as i32 && read_end > reg.start.0 as i32 {
                        same_start_ends.push(read_end);
                    }
                    continue;
                }
                any_passed = true;
                preprocessed = true;
            }
            if crate::cigar::parse::get_soft_clipped_length(&r.cigar) <= cur_max as u32 {
                if read_start <= reg.end.0 as i32 && read_end > reg.start.0 as i32 {
                    same_start_ends.push(read_end);
                }
                continue;
            }
            if !preprocessed && !preprocessor.preprocess_record(&r, config) {
                if read_start <= reg.end.0 as i32 && read_end > reg.start.0 as i32 {
                    same_start_ends.push(read_end);
                }
                continue;
            }
            any_passed = true;

            if read_start > reg.end.0 as i32 || read_end <= reg.start.0 as i32 {
                continue;
            }

            let insertion_deletion_length =
                crate::cigar::parse::get_insertion_deletion_length(&r.cigar);
            if let Some(nm_val) = r.nm {
                let total_mismatches = nm_val.saturating_sub(insertion_deletion_length as i32);
                if total_mismatches > config.max_mismatches as i32 {
                    same_start_ends.push(read_end);
                    continue;
                }
            }

            cigar_elements.clear();
            cigar_elements.extend_from_slice(&r.cigar);

            if config.realign == 1 {
                let qlen = r.sequence.len() as i32;
                let span_end = read_end.max(read_start.saturating_add(qlen));
                let ext_start_i32 = read_start.saturating_sub(qlen);
                let ext_end_i32 = span_end.saturating_add(qlen);
                let (ref_slice_vec, ref_start) = {
                    let (s, st) =
                        scanner.slice_for_read(read_start, span_end, ext_start_i32, ext_end_i32)?;
                    (s.to_vec(), st)
                };
                cigar_str.clear();
                crate::bam::push_cigar_string(&mut cigar_str, &cigar_elements);
                let mut seed_err: Option<VardictError> = None;
                let mut seed_lookup = |rc_seed: &[u8]| -> Option<u32> {
                    match scanner.unique_seed12_pos(rc_seed, read_start, span_end, cur_max) {
                        Ok(pos) => pos,
                        Err(e) => {
                            if seed_err.is_none() {
                                seed_err = Some(e);
                            }
                            None
                        }
                    }
                };
                let modifier = crate::cigar::modify::CigarModifier::with_seed_fn(
                    read_start,
                    std::mem::take(&mut cigar_str),
                    r.sequence.clone(),
                    r.quality_scores.clone(),
                    &ref_slice_vec,
                    ref_start,
                    &mut seed_lookup,
                    insertion_deletion_length as i32,
                    cur_max,
                    config,
                );
                let mod_res = modifier.modify_cigar();
                if let Some(e) = seed_err {
                    return Err(e);
                }
                match mod_res {
                    Ok(modified) => {
                        if crate::cigar::parse_cigar_into(&modified.cigar, &mut cigar_elements)
                            .is_none()
                        {
                            cigar_elements.clear();
                            cigar_elements.extend_from_slice(&r.cigar);
                        }
                        cigar_str = modified.cigar;
                    }
                    Err(_) => {
                        same_start_ends.push(read_end);
                        continue;
                    }
                }
            }

            use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
            if cigar_elements.len() >= 2 {
                if let (
                    Some(&(len_first, CigarKind::SoftClip)),
                    Some(&(len_last, CigarKind::SoftClip)),
                ) = (cigar_elements.first(), cigar_elements.last())
                {
                    if (10..=99).contains(&len_first) && len_last >= 10 {
                        same_start_ends.push(read_end);
                        continue;
                    }
                }
            }

            let read_len_match_ins =
                crate::cigar::parse::get_match_insertion_length(&cigar_elements);
            if config.minmatch != 0 && read_len_match_ins < config.minmatch {
                same_start_ends.push(read_end);
                continue;
            }
            let total_length_including_soft_clipped =
                crate::cigar::parse::get_soft_clipped_length(&cigar_elements);
            if total_length_including_soft_clipped > cur_max as u32 {
                cur_max = total_length_including_soft_clipped as i32;
                for (tile, steps) in reg_tiles.iter().zip(steps_per_tile.iter_mut()) {
                    if r_start_u32 <= tile.end.0 || idx + 1 < bam_paths.len() {
                        let prev_ov = same_start_ends
                            .iter()
                            .filter(|&&e| e > tile.start.0 as i32)
                            .count();
                        let self_ov = usize::from(read_end > tile.start.0 as i32);
                        steps.push(MaxRlStep {
                            bam_idx: idx,
                            read_start: r_start_u32,
                            same_start_overlap_before: prev_ov + self_ov,
                            new_max: cur_max,
                        });
                    }
                }
            }
            same_start_ends.push(read_end);
        }
        if idx > 0 && !any_passed {
            break;
        }
    }

    Ok((steps_per_tile, cur_max))
}

const DEFAULT_VCF_SPILL_BYTES: usize = 8 * 1024 * 1024; // 8 MiB

fn write_spill_u8<W: std::io::Write>(w: &mut W, v: u8) -> std::io::Result<()> {
    w.write_all(&[v])
}

fn write_spill_i32<W: std::io::Write>(w: &mut W, v: i32) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn write_spill_f64<W: std::io::Write>(w: &mut W, v: f64) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn write_spill_str<W: std::io::Write>(w: &mut W, s: &str) -> std::io::Result<()> {
    let len = s.len() as u64;
    w.write_all(&len.to_le_bytes())?;
    w.write_all(s.as_bytes())
}

fn read_spill_i32_opt<R: std::io::Read>(r: &mut R) -> std::io::Result<Option<i32>> {
    let mut buf = [0u8; 4];
    let mut read = 0usize;
    while read < 4 {
        match r.read(buf.get_mut(read..).unwrap_or(&mut []))? {
            0 if read == 0 => return Ok(None),
            0 => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "truncated spill record",
                ));
            }
            n => read += n,
        }
    }
    Ok(Some(i32::from_le_bytes(buf)))
}

fn read_spill_i32<R: std::io::Read>(r: &mut R) -> std::io::Result<i32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_spill_f64<R: std::io::Read>(r: &mut R) -> std::io::Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

fn read_spill_u8<R: std::io::Read>(r: &mut R) -> std::io::Result<u8> {
    let mut buf = [0u8; 1];
    r.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_spill_str<R: std::io::Read>(r: &mut R) -> std::io::Result<String> {
    let mut len_buf = [0u8; 8];
    r.read_exact(&mut len_buf)?;
    let len = u64::from_le_bytes(len_buf) as usize;
    let mut bytes = vec![0u8; len];
    r.read_exact(&mut bytes)?;
    String::from_utf8(bytes).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

trait SpillableCandidate: Sized {
    fn start_pos(&self) -> i32;
    fn af_4dp(&self) -> f64;
    fn approx_heap_bytes(&self) -> usize;
    fn write_to<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()>;
    fn read_from<R: std::io::Read>(r: &mut R) -> std::io::Result<Option<Self>>;
}

struct RawVcfCandidate {
    start_pos: i32,
    af_4dp: f64,
    should_emit: bool,
    vcf_line: String,
}

impl SpillableCandidate for RawVcfCandidate {
    #[inline]
    fn start_pos(&self) -> i32 {
        self.start_pos
    }

    #[inline]
    fn af_4dp(&self) -> f64 {
        self.af_4dp
    }

    #[inline]
    fn approx_heap_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + self.vcf_line.len()
    }

    fn write_to<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        write_spill_i32(w, self.start_pos)?;
        write_spill_f64(w, self.af_4dp)?;
        write_spill_u8(w, u8::from(self.should_emit))?;
        write_spill_str(w, &self.vcf_line)
    }

    fn read_from<R: std::io::Read>(r: &mut R) -> std::io::Result<Option<Self>> {
        let Some(start_pos) = read_spill_i32_opt(r)? else {
            return Ok(None);
        };
        let af_4dp = read_spill_f64(r)?;
        let should_emit = read_spill_u8(r)? != 0;
        let vcf_line = read_spill_str(r)?;
        Ok(Some(Self {
            start_pos,
            af_4dp,
            should_emit,
            vcf_line,
        }))
    }
}

struct SomaticVcfCandidate {
    start_pos: i32,
    af_4dp: f64,
    record: Option<crate::output::vcf::PairedVcfRecord>,
}

impl SpillableCandidate for SomaticVcfCandidate {
    #[inline]
    fn start_pos(&self) -> i32 {
        self.start_pos
    }

    #[inline]
    fn af_4dp(&self) -> f64 {
        self.af_4dp
    }

    #[inline]
    fn approx_heap_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.record.as_ref().map_or(0, |rec| {
                rec.prefix.len() + rec.filter.len() + rec.rest.len()
            })
    }

    fn write_to<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        write_spill_i32(w, self.start_pos)?;
        write_spill_f64(w, self.af_4dp)?;
        match &self.record {
            None => write_spill_u8(w, 0),
            Some(rec) => {
                write_spill_u8(w, 1)?;
                write_spill_str(w, &rec.prefix)?;
                write_spill_str(w, &rec.filter)?;
                write_spill_str(w, &rec.rest)?;
                write_spill_i32(w, rec.start)?;
                write_spill_u8(w, u8::from(rec.type_snv))?;
                write_spill_u8(w, u8::from(rec.status_somatic))
            }
        }
    }

    fn read_from<R: std::io::Read>(r: &mut R) -> std::io::Result<Option<Self>> {
        let Some(start_pos) = read_spill_i32_opt(r)? else {
            return Ok(None);
        };
        let af_4dp = read_spill_f64(r)?;
        let has_rec = read_spill_u8(r)? != 0;
        let record = if has_rec {
            let prefix = read_spill_str(r)?;
            let filter = read_spill_str(r)?;
            let rest = read_spill_str(r)?;
            let start = read_spill_i32(r)?;
            let type_snv = read_spill_u8(r)? != 0;
            let status_somatic = read_spill_u8(r)? != 0;
            Some(crate::output::vcf::PairedVcfRecord {
                prefix,
                filter,
                rest,
                start,
                type_snv,
                status_somatic,
            })
        } else {
            None
        };
        Ok(Some(Self {
            start_pos,
            af_4dp,
            record,
        }))
    }
}

struct SpillTempDir {
    path: std::path::PathBuf,
}

impl SpillTempDir {
    fn new() -> std::io::Result<Self> {
        static NEXT_DIR_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let pid = std::process::id();
        let base = std::env::temp_dir();
        let mut builder = std::fs::DirBuilder::new();
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt as _;
            builder.mode(0o700);
        }
        for _ in 0..1024 {
            let id = NEXT_DIR_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let candidate = base.join(format!("vardictrust_vcf_spill_{pid}_{id}"));
            match builder.create(&candidate) {
                Ok(()) => return Ok(Self { path: candidate }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
        Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "failed to allocate unique VCF spill directory",
        ))
    }

    fn path(&self) -> &std::path::Path {
        &self.path
    }
}

fn create_private_spill_file(path: &std::path::Path) -> std::io::Result<std::fs::File> {
    let mut opts = std::fs::OpenOptions::new();
    opts.write(true).create_new(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt as _;
        opts.mode(0o600);
    }
    opts.open(path)
}

impl Drop for SpillTempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

const MAX_SPILL_MERGE_FANIN: usize = 64;

enum RunSource<T> {
    File(std::io::BufReader<std::fs::File>),
    Mem(std::collections::btree_map::IntoValues<i32, T>),
}

impl<T: SpillableCandidate> RunSource<T> {
    fn next_item(&mut self) -> Result<Option<T>, VardictError> {
        match self {
            Self::File(r) => T::read_from(r).map_err(VardictError::Io),
            Self::Mem(it) => Ok(it.next()),
        }
    }
}

fn merge_run_sources<V, F>(
    mut sources: Vec<RunSource<V>>,
    mut on_winner: F,
) -> Result<(), VardictError>
where
    V: SpillableCandidate,
    F: FnMut(V) -> Result<(), VardictError>,
{
    let mut heads: Vec<Option<V>> = Vec::with_capacity(sources.len());
    let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<(i32, usize)>> =
        std::collections::BinaryHeap::new();
    for (run_idx, src) in sources.iter_mut().enumerate() {
        let first = src.next_item()?;
        if let Some(ref item) = first {
            heap.push(std::cmp::Reverse((item.start_pos(), run_idx)));
        }
        heads.push(first);
    }

    while let Some(std::cmp::Reverse((pos, run_idx))) = heap.pop() {
        let Some(mut winner) = heads.get_mut(run_idx).and_then(Option::take) else {
            continue;
        };
        if let Some(src) = sources.get_mut(run_idx) {
            if let Some(next_item) = src.next_item()? {
                let next_pos = next_item.start_pos();
                if let Some(slot) = heads.get_mut(run_idx) {
                    *slot = Some(next_item);
                }
                heap.push(std::cmp::Reverse((next_pos, run_idx)));
            }
        }
        while let Some(&std::cmp::Reverse((next_pos, next_run_idx))) = heap.peek() {
            if next_pos != pos {
                break;
            }
            heap.pop();
            if let Some(contender) = heads.get_mut(next_run_idx).and_then(Option::take) {
                if contender.af_4dp() > winner.af_4dp() {
                    winner = contender;
                }
            }
            if let Some(src) = sources.get_mut(next_run_idx) {
                if let Some(next_item) = src.next_item()? {
                    let np = next_item.start_pos();
                    if let Some(slot) = heads.get_mut(next_run_idx) {
                        *slot = Some(next_item);
                    }
                    heap.push(std::cmp::Reverse((np, next_run_idx)));
                }
            }
        }
        on_winner(winner)?;
    }
    Ok(())
}

/// Bounded-memory per-chromosome VCF candidate accumulator (`AGENTS.md §10`).
/// Keeps up to `max_bytes` of per-chromosome sorted `BTreeMap<i32, V>` entries in RAM,
/// spilling sorted runs to a temporary directory when exceeded and k-way merging runs
/// in chronological tile order on final chromosome emission.
struct SpillableVcfAccumulator<V: SpillableCandidate> {
    in_memory: indexmap::IndexMap<String, std::collections::BTreeMap<i32, V>>,
    spilled_runs: indexmap::IndexMap<String, Vec<std::path::PathBuf>>,
    current_bytes: usize,
    max_bytes: usize,
    max_fanin: usize,
    temp_dir: Option<SpillTempDir>,
    next_spill_id: usize,
}

impl<V: SpillableCandidate> SpillableVcfAccumulator<V> {
    fn new() -> Self {
        Self::with_max_bytes(DEFAULT_VCF_SPILL_BYTES)
    }

    fn with_max_bytes(max_bytes: usize) -> Self {
        Self::with_limits(max_bytes, MAX_SPILL_MERGE_FANIN)
    }

    fn with_limits(max_bytes: usize, max_fanin: usize) -> Self {
        Self {
            in_memory: indexmap::IndexMap::new(),
            spilled_runs: indexmap::IndexMap::new(),
            current_bytes: 0,
            max_bytes,
            max_fanin: max_fanin.max(2),
            temp_dir: None,
            next_spill_id: 0,
        }
    }

    fn insert(&mut self, chr: &str, cand: V) -> Result<(), VardictError> {
        if !self.spilled_runs.contains_key(chr) {
            self.spilled_runs.insert(chr.to_string(), Vec::new());
        }
        let chr_map = self.in_memory.entry(chr.to_string()).or_default();
        let new_bytes = cand.approx_heap_bytes() + 32;
        match chr_map.entry(cand.start_pos()) {
            std::collections::btree_map::Entry::Vacant(e) => {
                self.current_bytes = self.current_bytes.saturating_add(new_bytes);
                e.insert(cand);
            }
            std::collections::btree_map::Entry::Occupied(mut e) => {
                if cand.af_4dp() > e.get().af_4dp() {
                    let old_bytes = e.get().approx_heap_bytes() + 32;
                    self.current_bytes = self
                        .current_bytes
                        .saturating_sub(old_bytes)
                        .saturating_add(new_bytes);
                    e.insert(cand);
                }
            }
        }
        if self.current_bytes >= self.max_bytes {
            self.spill_to_disk()?;
        }
        Ok(())
    }

    fn spill_to_disk(&mut self) -> Result<(), VardictError> {
        if self.in_memory.is_empty() {
            self.current_bytes = 0;
            return Ok(());
        }
        let dir = match &self.temp_dir {
            Some(d) => d,
            None => {
                let d = SpillTempDir::new().map_err(VardictError::Io)?;
                self.temp_dir.insert(d)
            }
        };
        for (chr, btree) in self.in_memory.drain(..) {
            if btree.is_empty() {
                continue;
            }
            let path = dir
                .path()
                .join(format!("vcf_spill_{}.bin", self.next_spill_id));
            self.next_spill_id += 1;
            let file = create_private_spill_file(&path).map_err(VardictError::Io)?;
            let mut writer = std::io::BufWriter::new(file);
            for (_, cand) in btree {
                cand.write_to(&mut writer).map_err(VardictError::Io)?;
            }
            use std::io::Write as _;
            writer.flush().map_err(VardictError::Io)?;
            self.spilled_runs.entry(chr).or_default().push(path);
        }
        self.current_bytes = 0;
        Ok(())
    }

    fn emit_all<F>(mut self, mut on_chr_winner: F) -> Result<(), VardictError>
    where
        F: FnMut(&str, V) -> Result<(), VardictError>,
    {
        if self.spilled_runs.is_empty() && self.in_memory.is_empty() {
            return Ok(());
        }
        let max_fanin = self.max_fanin.max(2);
        let chr_keys: Vec<String> = self.spilled_runs.keys().cloned().collect();
        for chr in crate::output::vcf::perl_reorder_chromosomes(&chr_keys) {
            let mem_map = self.in_memory.swap_remove(&chr).unwrap_or_default();
            let mut run_paths = self.spilled_runs.swap_remove(&chr).unwrap_or_default();
            if run_paths.is_empty() {
                for (_, cand) in mem_map {
                    on_chr_winner(&chr, cand)?;
                }
                continue;
            }

            // Multi-pass bounded fan-in merge so open file descriptors never exceed `max_fanin`.
            // Merging contiguous chronological chunks preserves first-seen AF tie-breaking.
            while run_paths.len() > max_fanin {
                let mut next_gen = Vec::with_capacity(run_paths.len().div_ceil(max_fanin));
                for chunk in run_paths.chunks(max_fanin) {
                    if chunk.len() == 1 {
                        next_gen.push(chunk[0].clone());
                        continue;
                    }
                    let dir = match &self.temp_dir {
                        Some(d) => d,
                        None => {
                            let d = SpillTempDir::new().map_err(VardictError::Io)?;
                            self.temp_dir.insert(d)
                        }
                    };
                    let out_path = dir
                        .path()
                        .join(format!("vcf_spill_{}.bin", self.next_spill_id));
                    self.next_spill_id += 1;
                    let mut chunk_sources: Vec<RunSource<V>> = Vec::with_capacity(chunk.len());
                    for p in chunk {
                        let f = std::fs::File::open(p).map_err(VardictError::Io)?;
                        chunk_sources.push(RunSource::File(std::io::BufReader::new(f)));
                    }
                    let out_file =
                        create_private_spill_file(&out_path).map_err(VardictError::Io)?;
                    let mut writer = std::io::BufWriter::new(out_file);
                    merge_run_sources(chunk_sources, |winner: V| {
                        winner.write_to(&mut writer).map_err(VardictError::Io)
                    })?;
                    use std::io::Write as _;
                    writer.flush().map_err(VardictError::Io)?;
                    for p in chunk {
                        let _ = std::fs::remove_file(p);
                    }
                    next_gen.push(out_path);
                }
                run_paths = next_gen;
            }

            let mut sources: Vec<RunSource<V>> = Vec::with_capacity(run_paths.len() + 1);
            for p in &run_paths {
                let f = std::fs::File::open(p).map_err(VardictError::Io)?;
                sources.push(RunSource::File(std::io::BufReader::new(f)));
            }
            if !mem_map.is_empty() {
                sources.push(RunSource::Mem(mem_map.into_values()));
            }
            merge_run_sources(sources, |winner| on_chr_winner(&chr, winner))?;
            for p in &run_paths {
                let _ = std::fs::remove_file(p);
            }
        }
        Ok(())
    }
}

/// Output destination wrapper around buffered stdout or stderr.
/// Uses unlocked BufWriter to eliminate mutex contention with worker diagnostics
/// while reducing syscall overhead by ~1000x for deep WGS output.
pub enum OutputDestination {
    Stdout(std::io::BufWriter<std::io::Stdout>),
    Stderr(std::io::BufWriter<std::io::Stderr>),
}

impl std::io::Write for OutputDestination {
    #[inline]
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Stdout(w) => w.write(buf),
            Self::Stderr(w) => w.write(buf),
        }
    }

    #[inline]
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Stdout(w) => w.flush(),
            Self::Stderr(w) => w.flush(),
        }
    }
}

const BUFFER_CAPACITY: usize = 64 * 1024; // 64 KB

/// Creates the persistent buffered output writer according to `--DP` / `--defaultPrinter`.
/// Defaults to stdout; routes to stderr when `--DP ERR` is passed and `-i` (`--splice`) is off
/// (ports CmdParser.java:171-177 and CigarParser.java:105).
pub fn create_output_writer(config: &config::Config) -> OutputDestination {
    if config.debug_print.as_deref() == Some("ERR") && !config.splice {
        OutputDestination::Stderr(std::io::BufWriter::with_capacity(
            BUFFER_CAPACITY,
            std::io::stderr(),
        ))
    } else {
        OutputDestination::Stdout(std::io::BufWriter::with_capacity(
            BUFFER_CAPACITY,
            std::io::stdout(),
        ))
    }
}

/// Internal writer abstraction that preserves `run_with_writer`'s complete-stream contract
/// on injected writers while routing `--DP ERR` column headers to `stdout` when invoked via
/// the CLI [`OutputDestination`].
trait EngineWriter: std::io::Write + Sized {
    fn emit_column_header(&mut self, config: &config::Config) -> Result<(), VardictError> {
        crate::output::table::print_header(self, config).map_err(VardictError::Io)
    }
}

impl EngineWriter for OutputDestination {
    fn emit_column_header(&mut self, config: &config::Config) -> Result<(), VardictError> {
        match self {
            Self::Stdout(w) => {
                crate::output::table::print_header(w, config).map_err(VardictError::Io)
            }
            Self::Stderr(_) => {
                // SimpleMode/SomaticMode/AmpliconMode/SplicingMode printHeader() always writes
                // the header line to System.out even when `--DP ERR` routes variant rows to System.err.
                let mut stdout = std::io::stdout();
                crate::output::table::print_header(&mut stdout, config)
                    .map_err(VardictError::Io)?;
                std::io::Write::flush(&mut stdout).map_err(VardictError::Io)
            }
        }
    }
}

struct InjectedWriter<W>(W);

impl<W: std::io::Write> std::io::Write for InjectedWriter<W> {
    #[inline]
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.write(buf)
    }

    #[inline]
    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

impl<W: std::io::Write> EngineWriter for InjectedWriter<W> {}

/// Backwards-compatible alias for existing call sites.
pub fn get_output_writer(config: &config::Config) -> OutputDestination {
    create_output_writer(config)
}

/// Writes a single line to the output handle, propagating I/O errors (including BrokenPipe)
/// so parallel execution terminates promptly without spinning.
pub fn write_output_line<W: std::io::Write>(
    handle: &mut W,
    line: &str,
) -> Result<(), VardictError> {
    writeln!(handle, "{}", line).map_err(VardictError::Io)
}

/// Calling mode resolved from CLI flags and BED auto-detection (`VarDictLauncher.java:46-66`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallingMode {
    Splicing,
    Amplicon,
    Somatic,
    Simple,
}

/// Ports `VarDictLauncher.start` (`VarDictLauncher.java:46-66`) mode precedence:
/// 1. `SplicingMode` (`-i`)
/// 2. `-R` region flag or non-amplicon -> `SomaticMode` (paired) or `SimpleMode` (single)
/// 3. `AmpliconMode` (`-a` or 8-col amplicon BED auto-detect)
pub fn resolve_mode(
    splice: bool,
    has_region: bool,
    amplicon_effective: bool,
    paired: bool,
) -> CallingMode {
    if splice {
        CallingMode::Splicing
    } else if has_region || !amplicon_effective {
        if paired {
            CallingMode::Somatic
        } else {
            CallingMode::Simple
        }
    } else {
        CallingMode::Amplicon
    }
}

/// Resolves effective [`config::Config`] and [`CallingMode`] matching `VarDictLauncher.initResources`
/// and `VarDictLauncher.start` (`VarDictLauncher.java:46-95, 139-173`).
fn resolve_effective_config(
    config: &config::Config,
) -> Result<(config::Config, CallingMode), VardictError> {
    let mut effective = config.clone();
    let paired = effective.bam.contains('|');
    let has_region = effective.region.is_some();

    if has_region {
        effective.amplicon = None;
        if paired && (effective.emit_intermediate || effective.splice) {
            if effective.sample_explicit {
                if let Some(first_sample) = effective.sample.split('|').next() {
                    effective.sample = first_sample.to_string();
                }
            } else if let Some(first_bam_spec) = effective.bam.split('|').next() {
                let first_bam = first_bam_spec.split(':').next().unwrap_or(first_bam_spec);
                effective.sample = crate::config::extract_sample_name(first_bam, None);
            }
        }
    } else if effective.amplicon.is_none() {
        if let Some(ref bed_path) = effective.bed {
            if crate::regions::detect_amplicon_bed(bed_path, &effective.delimiter)? {
                effective.amplicon = Some("10:0.95".to_string());
                if !effective.zero_based_defined {
                    effective.zero_based = true;
                    effective.zero_based_defined = true;
                }
            }
        }
    }

    let mode = resolve_mode(
        effective.splice,
        has_region,
        effective.amplicon.is_some(),
        paired,
    );

    if mode == CallingMode::Amplicon && paired && !effective.emit_intermediate {
        return Err(VardictError::Unsupported(
            "Paired BAM (-b tumor|normal) with amplicon mode (-a) is not supported in VCF mode (#378); use .var intermediate output or single-sample BAM"
                .to_string(),
        ));
    }

    if paired && matches!(mode, CallingMode::Splicing | CallingMode::Amplicon) {
        if let Some(first_bam_spec) = effective.bam.split('|').next() {
            effective.bam = first_bam_spec.to_string();
        }
    }

    Ok((effective, mode))
}

/// Main entry point for the engine. Called by the CLI binary.
///
/// Takes a fully-resolved configuration and runs the variant-calling pipeline.
pub fn run(config: &config::Config) -> Result<(), VardictError> {
    struct BamCacheGuard;
    impl Drop for BamCacheGuard {
        fn drop(&mut self) {
            crate::bam::clear_bam_caches();
        }
    }
    crate::bam::clear_bam_caches();
    let _bam_cache_guard = BamCacheGuard;

    crate::runner::run_with_worker_stack_for_bam(&config.bam, || {
        let (effective, _) = resolve_effective_config(config)?;
        let writer = create_output_writer(&effective);
        run_with_engine_writer(&effective, writer)
    })
}

/// Runs the variant-calling pipeline writing all output to `writer` without clearing
/// process-wide BAM/CRAM/FASTA caches before/after the call.
pub fn run_with_writer<W: std::io::Write>(
    config: &config::Config,
    writer: W,
) -> Result<(), VardictError> {
    run_with_engine_writer(config, InjectedWriter(writer))
}

fn run_with_engine_writer<W: EngineWriter>(
    config: &config::Config,
    mut writer: W,
) -> Result<(), VardictError> {
    match run_with_writer_inner(config, &mut writer) {
        Ok(()) => writer.flush().map_err(VardictError::Io),
        Err(e) => {
            let _ = writer.flush();
            Err(e)
        }
    }
}

fn run_with_writer_inner<W: EngineWriter>(
    config: &config::Config,
    mut writer: &mut W,
) -> Result<(), VardictError> {
    log::info!("VarDictRust engine starting");

    let (effective, mode) = resolve_effective_config(config)?;
    let config = &effective;

    if mode == CallingMode::Somatic {
        return run_somatic(config, writer);
    }

    let first_bam_spec = config
        .bam
        .split('|')
        .next()
        .ok_or_else(|| VardictError::Parse("No BAM files specified".to_string()))?;
    let bam_list = parse_colon_bams(first_bam_spec, "No BAM files specified")?;
    let first_bam = bam_list
        .first()
        .map(|s| s.as_str())
        .ok_or_else(|| VardictError::Parse("No BAM files specified".to_string()))?;

    // Read reference lengths from cached BAM/CRAM header
    let chromosomes_lengths = crate::bam::read_reference_lengths(Path::new(first_bam), config)?;

    // Parse target regions
    let regions = if let Some(ref r_str) = config.region {
        crate::regions::parse_region_flag(r_str, config, &chromosomes_lengths)?
    } else if let Some(ref bed_path) = config.bed {
        crate::regions::parse_bed(
            bed_path,
            config,
            &chromosomes_lengths,
            config.zero_based_defined.then_some(config.zero_based),
        )?
    } else {
        return Err(VardictError::Parse(
            "Either region (-R) or BED file must be specified".to_string(),
        ));
    };

    // Amplicon mode: process per-amplicon and combine via amplicon post-process
    if mode == CallingMode::Amplicon {
        return run_amplicon(config, &chromosomes_lengths, regions, writer);
    }

    let flat_regions: Vec<crate::regions::Region> = regions.into_iter().flatten().collect();

    // Tile regions into overlapping tiles (regions <= 1 Mbp execute untiled like
    // VarDictJava; multi-megabase intervals > 1 Mbp tile into 100 kbp windows).
    let tile_size = config.tile_size.unwrap_or(DEFAULT_TILE_SIZE);
    let halo = 500;
    let mut tiles = Vec::new();
    let mut multi_tile_seeds: HashMap<MultiTileKey, (Vec<MaxRlStep>, i32)> = HashMap::new();
    for reg in &flat_regions {
        let reg_tiles = crate::regions::tile(reg, effective_tile_size(reg, config.tile_size), halo);
        if reg_tiles.len() > 1 && !config.splice {
            let (steps_per_tile, final_max_rl) = scan_multi_tile_max_read_lengths(
                &bam_list,
                reg,
                &reg_tiles,
                config,
                &chromosomes_lengths,
                0,
            )?;
            for (t, steps) in reg_tiles.iter().zip(steps_per_tile) {
                multi_tile_seeds.insert(multi_tile_key(t), (steps, final_max_rl));
            }
        }
        tiles.extend(reg_tiles);
    }
    let tile_batches = crate::regions::batch_tiles(tiles, tile_size, 16, config.threads);

    let config_clone = config.clone();
    let chromosomes_lengths_clone = chromosomes_lengths.clone();

    // -h/--header column line (CmdParser.java:67; SimpleMode.java:33 calls
    // printHeader() before any variant output).
    if config.column_header {
        writer.emit_column_header(config)?;
    }

    // Print VCF header when not in intermediate or splicing mode
    if !config.emit_intermediate && !config.splice {
        crate::output::vcf::print_vcf_header(&mut writer, config, &chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    enum SimpleTileOutput {
        Lines {
            lines: Vec<String>,
            is_split_tile: bool,
            is_first_split_tile: bool,
        },
        Vcf {
            chr: String,
            records: Vec<RawVcfCandidate>,
        },
    }

    let mut vcf_acc: SpillableVcfAccumulator<RawVcfCandidate> = SpillableVcfAccumulator::new();
    let mut seen_split_var: rustc_hash::FxHashSet<(u64, u64)> = rustc_hash::FxHashSet::default();

    // Process tiles in parallel, print variants in deterministic index order
    crate::runner::run_parallel_batched_for_bam(
        &config.bam,
        tile_batches,
        config.threads,
        move |tile| -> Result<SimpleTileOutput, VardictError> {
            crate::runner::set_current_worker_region(Some(format!(
                "{}:{}-{}",
                tile.chr, tile.start.0, tile.end.0
            )));
            let bed_region = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.halo_start,
                end: tile.halo_end,
                gene: tile.gene.clone(),
                amplicon: tile.amplicon.clone(),
            };

            let region_for_pileup = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.start,
                end: tile.end,
                gene: tile.gene.clone(),
                amplicon: tile.amplicon.clone(),
            };

            let mut reference = crate::reference::Reference::new();
            // Splice pileup runs halo-wide (bed_region) so the owning tile
            // counts junctions anchored at its last active base, then emits
            // owned-only (single-ownership, no cross-tile state).
            let pileup_region = if config_clone.splice {
                &bed_region
            } else {
                &region_for_pileup
            };
            let is_split_tile = tile.start != tile.original_start || tile.end != tile.original_end;
            let (steps_slice, post_max_rl, parent_orig_start) =
                if is_split_tile && !multi_tile_seeds.is_empty() {
                    if let Some((steps, final_max)) = multi_tile_seeds.get(&multi_tile_key(&tile)) {
                        (
                            steps.as_slice(),
                            Some(*final_max),
                            Some(tile.original_start),
                        )
                    } else {
                        (&[] as &[MaxRlStep], None, Some(tile.original_start))
                    }
                } else {
                    (&[] as &[MaxRlStep], None, None)
                };
            let var_data = run_pipeline_on_records_ext(
                &bam_list,
                &tile.chr,
                pileup_region.start.to_zero_based(),
                crate::types::Pos0(pileup_region.end.0),
                Some(&bed_region),
                &mut reference,
                pileup_region,
                if config_clone.splice {
                    PipelineTail::Stop
                } else {
                    PipelineTail::Full
                },
                &config_clone,
                &chromosomes_lengths_clone,
                0,
                None,
                steps_slice,
                post_max_rl,
                parent_orig_start,
            )?;

            let mut lines = Vec::new();

            if config_clone.splice {
                let mut introns: Vec<&String> = var_data.splice_count.keys().collect();
                introns.sort();
                for intron in introns {
                    if !tile_owns_splice_key(&tile, tile.original_start, tile.original_end, intron)
                    {
                        continue;
                    }
                    let count = var_data.splice_count.get(intron).map(|a| a[0]).unwrap_or(0);
                    lines.push(format!(
                        "{}\t{}\t{}\t{}",
                        config_clone.sample, tile.chr, intron, count
                    ));
                }
                return Ok(SimpleTileOutput::Lines {
                    lines,
                    is_split_tile: false,
                    is_first_split_tile: false,
                });
            }

            let builder_region = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.original_start,
                end: tile.original_end,
                gene: tile.gene.clone(),
                amplicon: tile.amplicon.clone(),
            };

            let builder = crate::tovars::ToVarsBuilder::new(
                &reference,
                &builder_region,
                &chromosomes_lengths_clone,
                &config_clone,
            );
            let aligned_variants = builder.build_aligned_variants(&var_data);
            let variants = crate::postprocess::simple::simple_post_process(
                &config_clone,
                &builder_region,
                &aligned_variants,
                &var_data.splice,
            );

            if config_clone.emit_intermediate {
                for v in variants {
                    let mut line = if config_clone.fisher {
                        crate::output::table::simple_38(
                            &v,
                            &builder_region,
                            &config_clone.sample,
                            &v.sv,
                        )
                    } else {
                        crate::output::table::simple_36(
                            &v,
                            &builder_region,
                            &config_clone.sample,
                            &v.sv,
                        )
                    };
                    if config_clone.crispr_cutting_site.unwrap_or(0) != 0 {
                        let _ =
                            std::fmt::Write::write_fmt(&mut line, format_args!("\t{}", v.crispr));
                    }
                    // ports SimpleOutputVariant.java:116-117 (DEBUG appended after crispr)
                    if config_clone.debug {
                        line.push('\t');
                        line.push_str(&v.debug);
                    }
                    lines.push(line);
                }
                Ok(SimpleTileOutput::Lines {
                    lines,
                    is_split_tile,
                    is_first_split_tile: tile.start == tile.original_start,
                })
            } else {
                let mut records = Vec::new();
                for v in variants {
                    let ref_empty = crate::output::vcf::perl_row_skipped(&v.ref_allele);
                    let is_sv = v.var_allele.contains('<');
                    let mut splitreads = 0;
                    if let Some((first, _)) = v.sv.split_once('-') {
                        splitreads = first.parse::<i32>().unwrap_or(0);
                    }
                    let should_emit = !ref_empty && (!is_sv || splitreads >= 1);
                    let vcf_line = if !should_emit {
                        String::new()
                    } else {
                        let fisher = crate::stats::FisherExact::new_with_mode(
                            v.ref_fwd_cov as i32,
                            v.ref_rev_cov as i32,
                            v.fwd as i32,
                            v.rev as i32,
                            config_clone.fisher,
                        );
                        let sbf = fisher.get_p_value();
                        let oddratio_str = fisher.get_odds_ratio();
                        let filters = crate::output::vcf::filter_variant(
                            &v,
                            &config_clone,
                            0,
                            None,
                            sbf,
                            &oddratio_str,
                        );
                        crate::output::vcf::format_vcf_record(
                            &v,
                            &builder_region,
                            &config_clone,
                            &filters,
                            sbf,
                            &oddratio_str,
                            None,
                        )
                    };
                    records.push(RawVcfCandidate {
                        start_pos: v.start_pos,
                        af_4dp: crate::types::round_half_even(v.frequency, 4),
                        should_emit,
                        vcf_line,
                    });
                }
                Ok(SimpleTileOutput::Vcf {
                    chr: tile.chr.clone(),
                    records,
                })
            }
        },
        |output| {
            match output {
                SimpleTileOutput::Lines {
                    lines,
                    is_split_tile,
                    is_first_split_tile,
                } => {
                    if !is_split_tile || is_first_split_tile {
                        seen_split_var.clear();
                    }
                    for line in lines {
                        if is_split_tile {
                            let fp = split_var_line_fingerprint(&line);
                            if seen_split_var.contains(&fp) {
                                continue;
                            }
                            if seen_split_var.len() < MAX_SPLIT_VAR_DEDUP_ENTRIES {
                                seen_split_var.insert(fp);
                            }
                        }
                        write_output_line(&mut writer, &line)?;
                    }
                }
                SimpleTileOutput::Vcf { chr, records } => {
                    for r in records {
                        vcf_acc.insert(&chr, r)?;
                    }
                }
            }
            Ok(())
        },
    )?;

    vcf_acc.emit_all(|_, cand| {
        if cand.should_emit {
            write_output_line(&mut writer, &cand.vcf_line)?;
        }
        Ok(())
    })?;
    writer.flush().map_err(VardictError::Io)?;

    Ok(())
}

/// Process a single amplicon region through the full pipeline, returning the
/// aligned variants map and accumulated splice set (`AmpliconMode.java:79, 98`).
/// The pileup runs over the amplicon region (start..end);
/// the ToVarsBuilder calls variants over the insert region.
fn process_amplicon_region(
    bam: &str,
    region: &crate::regions::Region,
    config: &config::Config,
    chromosomes_lengths: &HashMap<String, usize>,
    seed_splice: std::collections::HashSet<String>,
) -> Result<
    (
        crate::javamap::JavaHashMap<i32, crate::variations::Vars>,
        std::collections::HashSet<String>,
    ),
    VardictError,
> {
    let halo = 500u32;
    let halo_start = crate::types::Pos0(region.start.0.saturating_sub(halo).saturating_sub(1));
    let halo_end = crate::types::Pos0(region.end.0 + halo);

    let mut reference = crate::reference::Reference::new();
    let bam_list = parse_colon_bams(bam, "No BAM files specified")?;
    let var_data = run_pipeline_on_records(
        &bam_list,
        &region.chr,
        halo_start,
        halo_end,
        Some(region),
        &mut reference,
        region,
        PipelineTail::Full,
        config,
        chromosomes_lengths,
        0,
        Some(seed_splice),
    )?;

    // Variant calling over the insert region (amplicon insert bounds)
    let builder_region = if let Some(ref amp) = region.amplicon {
        crate::regions::Region {
            chr: region.chr.clone(),
            start: amp.insert_start,
            end: amp.insert_end,
            gene: region.gene.clone(),
            amplicon: region.amplicon.clone(),
        }
    } else {
        region.clone()
    };

    let builder =
        crate::tovars::ToVarsBuilder::new(&reference, &builder_region, chromosomes_lengths, config);
    let aligned = builder.build_aligned_variants(&var_data);
    Ok((aligned, var_data.splice))
}

/// Covering amplicons at one position: `(amplicon_number, region)` for every
/// segment entry whose insert interval contains `position`, in segment order.
///
/// Interval-overlap form of the per-base expansion in AmpliconMode.java:57-60
/// (`for p in insertStart..insertEnd`): identical `(position → list)` entries
/// without materializing one map entry per base (audit-amplicon-expand).
///
/// Single-query O(A) scan kept for isolated lookups and as the test oracle;
/// the hot paths (`amplicon_position_stream`, non-`-p` `run_amplicon`) sweep
/// via `AmpliconCoverSweep` instead.
#[allow(dead_code)]
pub(crate) fn amplicons_covering_position(
    segment: &[crate::regions::Region],
    position: i32,
) -> Vec<(usize, crate::regions::Region)> {
    segment
        .iter()
        .enumerate()
        .filter(|(_, region)| {
            region.amplicon.as_ref().is_some_and(|amp| {
                position >= amp.insert_start.0 as i32 && position <= amp.insert_end.0 as i32
            })
        })
        .map(|(amplicon_number, region)| (amplicon_number, region.clone()))
        .collect()
}

/// Amplicon mode: process each amplicon region separately, then combine via
/// amplicon post-processing. Ports AmpliconMode.java.
fn run_amplicon<W: EngineWriter>(
    config: &config::Config,
    chromosomes_lengths: &HashMap<String, usize>,
    segments: Vec<Vec<crate::regions::Region>>,
    writer: &mut W,
) -> Result<(), VardictError> {
    let first_bam = config
        .bam
        .split('|')
        .next()
        .ok_or_else(|| VardictError::Parse("No BAM files specified".to_string()))?;

    // -h/--header column line (AmpliconMode.java:39, 135-146).
    if config.column_header {
        writer.emit_column_header(config)?;
    }

    // Print VCF header once when not in intermediate mode (same var2vcf_valid.pl
    // header as simple mode — amplicon rides the same writer, perl :140-145).
    if !config.emit_intermediate {
        crate::output::vcf::print_vcf_header(writer, config, chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    let mut vcf_acc: SpillableVcfAccumulator<RawVcfCandidate> = SpillableVcfAccumulator::new();
    fn format_amp_cand(
        av: &crate::postprocess::amplicon::AmpliconVariant,
        current_region: &crate::regions::Region,
        config: &config::Config,
    ) -> RawVcfCandidate {
        let v = &av.variant;
        let should_emit = !crate::output::vcf::perl_row_skipped(&v.ref_allele);
        let vcf_line = if !should_emit {
            String::new()
        } else {
            let fisher = crate::stats::FisherExact::new_with_mode(
                v.ref_fwd_cov as i32,
                v.ref_rev_cov as i32,
                v.fwd as i32,
                v.rev as i32,
                config.fisher,
            );
            let sbf = fisher.get_p_value();
            let oddratio_str = fisher.get_odds_ratio();
            let amp = Some((av.gdamp, av.tlamp, av.ncamp, av.ampflag));
            let filters = crate::output::vcf::filter_variant(v, config, 0, amp, sbf, &oddratio_str);
            let mut line = String::new();
            crate::output::vcf::format_vcf_record_into(
                &mut line,
                v,
                current_region,
                config,
                &filters,
                sbf,
                &oddratio_str,
                amp,
            );
            line
        };
        RawVcfCandidate {
            start_pos: v.start_pos,
            af_4dp: crate::types::round_half_even(v.frequency, 4),
            should_emit,
            vcf_line,
        }
    }

    type AmpliconVarsList = Vec<crate::javamap::JavaHashMap<i32, crate::variations::Vars>>;

    enum AmpliconSegmentOutput {
        Lines(Vec<String>),
        Vcf {
            chr: String,
            records: Vec<RawVcfCandidate>,
        },
        Pileup {
            segment: Vec<crate::regions::Region>,
            vars_per_amplicon: AmpliconVarsList,
            splice: std::collections::HashSet<String>,
        },
    }

    let segments: Vec<Vec<crate::regions::Region>> =
        segments.into_iter().filter(|s| !s.is_empty()).collect();
    let first_bam_owned = first_bam.to_string();
    let config_clone = config.clone();
    let chromosomes_lengths_clone = chromosomes_lengths.clone();

    // H5: Process independent amplicon segments in parallel via `run_parallel_for_bam`,
    // streaming results into the consumer in exact BED segment order.
    crate::runner::run_parallel_for_bam(
        &config.bam,
        segments,
        config.threads,
        move |segment| -> Result<AmpliconSegmentOutput, VardictError> {
            // §9: let-else merges the empty guard with the last() unwrap.
            let Some(current_region) = segment.last() else {
                return Ok(AmpliconSegmentOutput::Lines(Vec::new()));
            };
            let run_segment_regions = || {
                let mut vars_per_amplicon: AmpliconVarsList = Vec::with_capacity(segment.len());
                let mut splice = std::collections::HashSet::new();
                for region in segment.iter() {
                    crate::runner::set_current_worker_region(Some(format!(
                        "{}:{}-{}",
                        region.chr, region.start.0, region.end.0
                    )));
                    let (aligned, next_splice) = process_amplicon_region(
                        &first_bam_owned,
                        region,
                        &config_clone,
                        &chromosomes_lengths_clone,
                        std::mem::take(&mut splice),
                    )?;
                    vars_per_amplicon.push(aligned);
                    splice = next_splice;
                }
                Ok::<_, VardictError>((vars_per_amplicon, splice))
            };

            let (vars_per_amplicon, splice) = if config_clone.threads <= 1 {
                crate::runner::run_with_worker_stack_for_bam(
                    &config_clone.bam,
                    run_segment_regions,
                )?
            } else {
                run_segment_regions()?
            };

            if config_clone.pileup {
                return Ok(AmpliconSegmentOutput::Pileup {
                    segment,
                    vars_per_amplicon,
                    splice,
                });
            }

            // Bounded-memory (audit-amplicon-expand): visit only positions with
            // variant data and resolve covering amplicons with one ascending
            // sweep (audit-amplicon-sweep). Positions absent from every vars
            // map `continue` without output (AmpliconPostProcessModule.java:129-131),
            // so skipping them is output-identical while bounding the map by
            // data, not interval size.
            let mut amplicons_on_positions: std::collections::BTreeMap<
                i32,
                Vec<(usize, crate::regions::Region)>,
            > = std::collections::BTreeMap::new();
            let mut positions: std::collections::BTreeSet<i32> = std::collections::BTreeSet::new();
            for vars in &vars_per_amplicon {
                positions.extend(vars.keys().copied());
            }
            // BTreeSet iterates ascending: satisfies the sweep precondition,
            // and each covering list comes out in segment order — the same
            // pairs the old per-position scan built, at O((P + A) log A).
            let mut sweep = crate::postprocess::amplicon::AmpliconCoverSweep::new(&segment);
            for p in positions {
                let covering = sweep.covering_at(p);
                if !covering.is_empty() {
                    amplicons_on_positions.insert(p, covering);
                }
            }

            let results = crate::postprocess::amplicon::amplicon_post_process(
                &config_clone,
                current_region,
                &vars_per_amplicon,
                &amplicons_on_positions,
                &splice,
            )?;

            if config_clone.emit_intermediate {
                let mut lines = Vec::with_capacity(results.len());
                for av in results {
                    let mut row_buf = String::new();
                    if config_clone.fisher {
                        crate::output::table::amplicon_40_into(
                            &mut row_buf,
                            &av,
                            current_region,
                            &config_clone.sample,
                            config_clone.debug,
                        );
                    } else {
                        crate::output::table::amplicon_38_into(
                            &mut row_buf,
                            &av,
                            current_region,
                            &config_clone.sample,
                            config_clone.debug,
                        );
                    }
                    lines.push(row_buf);
                }
                Ok(AmpliconSegmentOutput::Lines(lines))
            } else {
                let records = results
                    .iter()
                    .map(|av| format_amp_cand(av, current_region, &config_clone))
                    .collect();
                Ok(AmpliconSegmentOutput::Vcf {
                    chr: current_region.chr.clone(),
                    records,
                })
            }
        },
        |output| {
            match output {
                AmpliconSegmentOutput::Lines(lines) => {
                    for line in lines {
                        write_output_line(writer, &line)?;
                    }
                }
                AmpliconSegmentOutput::Vcf { chr, records } => {
                    for r in records {
                        vcf_acc.insert(&chr, r)?;
                    }
                }
                AmpliconSegmentOutput::Pileup {
                    segment,
                    vars_per_amplicon,
                    splice,
                } => {
                    let Some(current_region) = segment.last() else {
                        return Ok(());
                    };
                    // Position → covering-amplicon pairs (AmpliconMode.java:57-60).
                    //
                    // -p emits a (null/ref) row at EVERY insert position, even with zero
                    // coverage in all amplicons, so visiting every insert position is
                    // load-bearing here — but the per-base map is streamed, never
                    // materialized (audit-amplicon-expand): the lazy union iterator yields
                    // exactly the positions (and segment-ordered covering lists) the old
                    // per-base BTreeMap held, in the same ascending order.
                    if config.emit_intermediate {
                        let mut row_buf = String::new();
                        crate::postprocess::amplicon::amplicon_post_process_streaming(
                            config,
                            current_region,
                            &vars_per_amplicon,
                            crate::postprocess::amplicon::amplicon_position_stream(&segment),
                            &splice,
                            &mut |av| {
                                if config.fisher {
                                    crate::output::table::amplicon_40_into(
                                        &mut row_buf,
                                        &av,
                                        current_region,
                                        &config.sample,
                                        config.debug,
                                    );
                                } else {
                                    crate::output::table::amplicon_38_into(
                                        &mut row_buf,
                                        &av,
                                        current_region,
                                        &config.sample,
                                        config.debug,
                                    );
                                }
                                write_output_line(writer, &row_buf)
                            },
                        )?;
                    } else {
                        crate::postprocess::amplicon::amplicon_post_process_streaming(
                            config,
                            current_region,
                            &vars_per_amplicon,
                            crate::postprocess::amplicon::amplicon_position_stream(&segment),
                            &splice,
                            &mut |av| {
                                vcf_acc.insert(
                                    &current_region.chr,
                                    format_amp_cand(&av, current_region, config),
                                )
                            },
                        )?;
                    }
                }
            }
            Ok(())
        },
    )?;

    vcf_acc.emit_all(|_, cand| {
        if cand.should_emit {
            write_output_line(writer, &cand.vcf_line)?;
        }
        Ok(())
    })?;

    writer.flush().map_err(VardictError::Io)?;

    Ok(())
}

/// Somatic (paired) mode: process tumor BAM, then normal BAM, then compare.
/// Ports SomaticMode.java's processBothBamsInPipeline.
fn run_somatic<W: EngineWriter>(
    config: &config::Config,
    writer: &mut W,
) -> Result<(), VardictError> {
    let parts: Vec<&str> = config.bam.split('|').collect();
    if parts.len() != 2 {
        return Err(VardictError::Parse(
            "Somatic mode requires exactly two BAM files separated by '|'".to_string(),
        ));
    }
    let bam1_bams = parse_colon_bams(parts[0], "No tumor BAM file specified")?;
    let bam2_bams = parse_colon_bams(parts[1], "No normal BAM file specified")?;
    let first_bam1 = bam1_bams
        .first()
        .ok_or_else(|| VardictError::Parse("No tumor BAM file specified".to_string()))?;

    // Read reference lengths from cached BAM/CRAM header
    let chromosomes_lengths = crate::bam::read_reference_lengths(Path::new(first_bam1), config)?;

    // Parse target regions
    let regions = if let Some(ref r_str) = config.region {
        crate::regions::parse_region_flag(r_str, config, &chromosomes_lengths)?
    } else if let Some(ref bed_path) = config.bed {
        crate::regions::parse_bed(
            bed_path,
            config,
            &chromosomes_lengths,
            config.zero_based_defined.then_some(config.zero_based),
        )?
    } else {
        return Err(VardictError::Parse(
            "Either region (-R) or BED file must be specified".to_string(),
        ));
    };

    let flat_regions: Vec<crate::regions::Region> = regions.into_iter().flatten().collect();

    // Tile regions (regions <= 1 Mbp execute untiled like VarDictJava;
    // multi-megabase intervals > 1 Mbp tile into 100 kbp windows).
    let tile_size = config.tile_size.unwrap_or(DEFAULT_TILE_SIZE);
    let halo = 500;
    let mut tiles = Vec::new();
    let mut multi_tile_seeds: HashMap<MultiTileKey, (Vec<MaxRlStep>, i32, Vec<MaxRlStep>, i32)> =
        HashMap::new();
    for reg in &flat_regions {
        let reg_tiles = crate::regions::tile(reg, effective_tile_size(reg, config.tile_size), halo);
        if reg_tiles.len() > 1 {
            let (bam1_steps, bam1_final) = scan_multi_tile_max_read_lengths(
                &bam1_bams,
                reg,
                &reg_tiles,
                config,
                &chromosomes_lengths,
                0,
            )?;
            let (bam2_steps, bam2_final) = scan_multi_tile_max_read_lengths(
                &bam2_bams,
                reg,
                &reg_tiles,
                config,
                &chromosomes_lengths,
                bam1_final,
            )?;
            for ((t, s1), s2) in reg_tiles.iter().zip(bam1_steps).zip(bam2_steps) {
                multi_tile_seeds.insert(multi_tile_key(t), (s1, bam1_final, s2, bam2_final));
            }
        }
        tiles.extend(reg_tiles);
    }
    let tile_batches = crate::regions::batch_tiles(tiles, tile_size, 16, config.threads);

    let config_clone = config.clone();
    let chromosomes_lengths_clone = chromosomes_lengths.clone();
    let bam1_bams_for_combine = bam1_bams.clone();
    let bam2_bams_for_combine = bam2_bams.clone();

    // -h/--header column line (SomaticMode.java:35, 134-150).
    if config.column_header {
        writer.emit_column_header(config)?;
    }

    // Paired VCF header when not in intermediate mode
    if !config.emit_intermediate {
        crate::output::vcf::print_paired_vcf_header(writer, config, &chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    enum SomaticTileOutput {
        Lines {
            lines: Vec<String>,
            is_split_tile: bool,
            is_first_split_tile: bool,
        },
        Vcf {
            chr: String,
            candidates: Vec<SomaticVcfCandidate>,
        },
    }

    let mut vcf_acc: SpillableVcfAccumulator<SomaticVcfCandidate> = SpillableVcfAccumulator::new();
    let mut seen_split_var: rustc_hash::FxHashSet<(u64, u64)> = rustc_hash::FxHashSet::default();

    // Process tiles: run both BAMs through pipeline, then somatic postprocess
    crate::runner::run_parallel_batched_for_bam(
        &config.bam,
        tile_batches,
        config.threads,
        move |tile| -> Result<SomaticTileOutput, VardictError> {
            crate::runner::set_current_worker_region(Some(format!(
                "{}:{}-{}",
                tile.chr, tile.start.0, tile.end.0
            )));
            // Process BAM1 (tumor)
            let bed_region = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.halo_start,
                end: tile.halo_end,
                gene: tile.gene.clone(),
                amplicon: tile.amplicon.clone(),
            };
            let region_for_pileup = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.start,
                end: tile.end,
                gene: tile.gene.clone(),
                amplicon: tile.amplicon.clone(),
            };

            let is_split_tile = tile.start != tile.original_start || tile.end != tile.original_end;
            let (steps1, post1, steps2, post2, parent_orig_start) =
                if is_split_tile && !multi_tile_seeds.is_empty() {
                    if let Some((s1, f1, s2, f2)) = multi_tile_seeds.get(&multi_tile_key(&tile)) {
                        (
                            s1.as_slice(),
                            Some(*f1),
                            s2.as_slice(),
                            Some(*f2),
                            Some(tile.original_start),
                        )
                    } else {
                        (
                            &[] as &[MaxRlStep],
                            None,
                            &[] as &[MaxRlStep],
                            None,
                            Some(tile.original_start),
                        )
                    }
                } else {
                    (&[] as &[MaxRlStep], None, &[] as &[MaxRlStep], None, None)
                };

            let mut reference = crate::reference::Reference::new();
            let mut var_data1 = run_pipeline_on_records_ext(
                &bam1_bams,
                &tile.chr,
                region_for_pileup.start.to_zero_based(),
                crate::types::Pos0(region_for_pileup.end.0),
                Some(&bed_region),
                &mut reference,
                &region_for_pileup,
                PipelineTail::Full,
                &config_clone,
                &chromosomes_lengths_clone,
                0,
                None,
                steps1,
                post1,
                parent_orig_start,
            )?;

            let max_read_length_bam1 = var_data1.max_read_length;

            // Inter-pass trim (ports SomaticMode.java:111-113, AbstractMode.java:44-56):
            // in VarDictJava, ToVarsBuilder for BAM1 runs inside pipeline(initialScope1)
            // and .join()s before BAM2 starts, observing only the initial region reference
            // load plus BAM1's own SV mate loads (never preloading unqueried deletion spans,
            // which would alter rightSeq/RSEQ when a realignlgdel deletion end exceeds
            // region.end + reference_extension). Build aligned1 now and drop pass-1's
            // pileup maps before pass 2 runs so the full VariationData is freed early.
            let builder_region = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.original_start,
                end: tile.original_end,
                gene: tile.gene.clone(),
                amplicon: tile.amplicon.clone(),
            };
            let builder1 = crate::tovars::ToVarsBuilder::new(
                &reference,
                &builder_region,
                &chromosomes_lengths_clone,
                &config_clone,
            );
            let aligned1 = builder1.build_aligned_variants(&var_data1);
            let splice1 = std::mem::take(&mut var_data1.splice);
            drop(var_data1);
            // Process BAM2 (normal), seeded with BAM1's maxReadLength and shared
            // `splice` set (`SomaticMode.java:104-122`); the reference loaded above
            // is reused (no second get_reference).
            let seed_bam2 = post1.map_or(max_read_length_bam1, |f1| f1.max(max_read_length_bam1));
            let post_bam2 = post2.map(|f2| f2.max(max_read_length_bam1));
            let mut var_data2 = run_pipeline_on_records_ext(
                &bam2_bams,
                &tile.chr,
                region_for_pileup.start.to_zero_based(),
                crate::types::Pos0(region_for_pileup.end.0),
                None,
                &mut reference,
                &region_for_pileup,
                PipelineTail::Full,
                &config_clone,
                &chromosomes_lengths_clone,
                seed_bam2,
                Some(splice1),
                steps2,
                post_bam2,
                parent_orig_start,
            )?;

            // Use max of both BAMs' maxReadLength (SomaticPostProcessModule.java:64)
            let combined_max_rl = max_read_length_bam1.max(var_data2.max_read_length);
            let splice = std::mem::take(&mut var_data2.splice);

            // Build aligned variants for BAM2 (aligned1 already built above).
            let builder2 = crate::tovars::ToVarsBuilder::new(
                &reference,
                &builder_region,
                &chromosomes_lengths_clone,
                &config_clone,
            );
            let aligned2 = builder2.build_aligned_variants(&var_data2);
            // Inter-pass trim (cont.): var_data2 is fully consumed here;
            // drop its maps before postprocess.
            drop(var_data2);

            // Run somatic postprocess with combineAnalysis pipeline callback
            let mut bams_for_combine =
                Vec::with_capacity(bam1_bams_for_combine.len() + bam2_bams_for_combine.len());
            bams_for_combine.extend(bam1_bams_for_combine.iter().cloned());
            bams_for_combine.extend(bam2_bams_for_combine.iter().cloned());
            let config_for_combine = config_clone.clone();
            let chr_lens_for_combine = chromosomes_lengths_clone.clone();
            let pipeline_fn = |combine_region: &crate::regions::Region,
                               splice_set: &mut std::collections::HashSet<String>,
                               max_rl: &mut i32|
             -> Result<
                crate::javamap::JavaHashMap<i32, crate::variations::Vars>,
                VardictError,
            > {
                // Query reads from both BAMs over the combine region; the
                // merged record set goes through one preprocessing pass.
                // Ports SomaticPostProcessModule.java:399 (`Reference ref = referenceResource.getReference(region)`):
                // start with a fresh Reference so `run_pipeline_on_records` loads
                // `combine_region` (`reference.seed` and `referenceSequences` bounds
                // scoped to `combine_region`, not the outer tile).
                let mut pipeline_ref = crate::reference::Reference::new();
                let mut var_data = run_pipeline_on_records(
                    &bams_for_combine,
                    &combine_region.chr,
                    crate::types::Pos0(combine_region.start.0.saturating_sub(1)),
                    crate::types::Pos0(combine_region.end.0),
                    Some(combine_region),
                    &mut pipeline_ref,
                    combine_region,
                    PipelineTail::Full,
                    &config_for_combine,
                    &chr_lens_for_combine,
                    *max_rl,
                    Some(splice_set.clone()),
                )?;
                *splice_set = std::mem::take(&mut var_data.splice);
                *max_rl = (*max_rl).max(var_data.max_read_length);

                // ToVarsBuilder
                let builder = crate::tovars::ToVarsBuilder::new(
                    &pipeline_ref,
                    combine_region,
                    &chr_lens_for_combine,
                    &config_for_combine,
                );
                Ok(builder.build_aligned_variants(&var_data))
            };

            let somatic_variants = crate::postprocess::somatic::somatic_post_process(
                &config_clone,
                &builder_region,
                &aligned1,
                &aligned2,
                &splice,
                Some(&pipeline_fn),
                combined_max_rl,
            );

            // Output 55-column format (or 61-column with --fisher), or the
            // paired VCF (var2vcf_paired.pl) when not in intermediate mode.
            if config_clone.emit_intermediate {
                let mut lines = Vec::new();
                for sv in somatic_variants {
                    let mut line = if config_clone.fisher {
                        crate::output::table::somatic_61(
                            &sv.variant1,
                            sv.variant2.as_ref(),
                            &builder_region,
                            &config_clone.sample,
                            &sv.status,
                            &sv.sv1,
                            &sv.sv2,
                            sv.end_is_variant2,
                        )
                    } else {
                        crate::output::table::somatic_55(
                            &sv.variant1,
                            sv.variant2.as_ref(),
                            &builder_region,
                            &config_clone.sample,
                            &sv.status,
                            &sv.sv1,
                            &sv.sv2,
                            sv.end_is_variant2,
                        )
                    };
                    // ports SomaticOutputVariant.java:193-194 (DEBUG = tumor variant's)
                    if config_clone.debug {
                        line.push('\t');
                        line.push_str(&sv.variant1.debug);
                        lines.push(line);
                    } else {
                        lines.push(line);
                    }
                }
                Ok(SomaticTileOutput::Lines {
                    lines,
                    is_split_tile,
                    is_first_split_tile: tile.start == tile.original_start,
                })
            } else {
                // Paired VCF mode: ports var2vcf_paired.pl:41,128-129,136 —
                // bucket rows by the raw .var start column ($hash{$chr}->{$a[3]}),
                // AF-contest per bucket (:128-129, top row only without -A), then
                // the `next unless ($ref)` guard on each SELECTED row.
                let mut buckets: std::collections::BTreeMap<
                    i32,
                    crate::postprocess::somatic::SomaticVariant,
                > = std::collections::BTreeMap::new();
                for sv in somatic_variants {
                    match buckets.entry(sv.variant1.start_pos) {
                        std::collections::btree_map::Entry::Vacant(e) => {
                            e.insert(sv);
                        }
                        std::collections::btree_map::Entry::Occupied(mut e) => {
                            let sv_af = crate::types::round_half_even(sv.variant1.frequency, 4);
                            let cur_af =
                                crate::types::round_half_even(e.get().variant1.frequency, 4);
                            if sv_af > cur_af {
                                e.insert(sv);
                            }
                        }
                    }
                }
                let candidates = buckets
                    .into_values()
                    .map(|sv| {
                        let start_pos = sv.variant1.start_pos;
                        let af_4dp = crate::types::round_half_even(sv.variant1.frequency, 4);
                        let record =
                            if crate::output::vcf::perl_row_skipped(&sv.variant1.ref_allele) {
                                None
                            } else {
                                Some(crate::output::vcf::format_paired_vcf_record(
                                    &sv,
                                    &builder_region.chr,
                                    &config_clone,
                                ))
                            };
                        SomaticVcfCandidate {
                            start_pos,
                            af_4dp,
                            record,
                        }
                    })
                    .collect();
                Ok(SomaticTileOutput::Vcf {
                    chr: builder_region.chr,
                    candidates,
                })
            }
        },
        |output| {
            match output {
                SomaticTileOutput::Lines {
                    lines,
                    is_split_tile,
                    is_first_split_tile,
                } => {
                    if !is_split_tile || is_first_split_tile {
                        seen_split_var.clear();
                    }
                    for line in lines {
                        if is_split_tile {
                            let fp = split_var_line_fingerprint(&line);
                            if seen_split_var.contains(&fp) {
                                continue;
                            }
                            if seen_split_var.len() < MAX_SPLIT_VAR_DEDUP_ENTRIES {
                                seen_split_var.insert(fp);
                            }
                        }
                        write_output_line(writer, &line)?;
                    }
                }
                SomaticTileOutput::Vcf { chr, candidates } => {
                    for cand in candidates {
                        vcf_acc.insert(&chr, cand)?;
                    }
                }
            }
            Ok(())
        },
    )?;

    let mut cur_chr: Option<String> = None;
    let mut streamer = crate::output::vcf::PairedRecordStreamer::new(0);
    vcf_acc.emit_all(|chr, cand| {
        if cur_chr.as_deref() != Some(chr) {
            if cur_chr.is_some() {
                let prev = std::mem::replace(
                    &mut streamer,
                    crate::output::vcf::PairedRecordStreamer::new(0),
                );
                prev.finish(|line| write_output_line(writer, line))?;
            }
            cur_chr = Some(chr.to_string());
        }
        if let Some(rec) = cand.record {
            streamer.push(rec, |line| write_output_line(writer, line))?;
        }
        Ok(())
    })?;
    if cur_chr.is_some() {
        streamer.finish(|line| write_output_line(writer, line))?;
    }

    writer.flush().map_err(VardictError::Io)?;

    Ok(())
}

/// Paired-VCF row selection: bucket by tumor-variant start (var2vcf_paired.pl
/// `$hash{$chr}->{$a[3]}`), keep the highest-AF row per position (`$ALL = 1`
/// without -A; replace-only-if-strictly-greater keeps first-seen ties like the
/// stable mergesort at :128, comparing the 4dp `.var` table AF `$d->[14]`),
/// THEN drop Perl-falsy REF rows (:136 `next unless ($ref)` runs on each
/// selected row inside the $ALL loop — an empty-ref winner means NOTHING is
/// emitted for that position).
#[cfg(test)]
fn select_paired_rows(
    somatic_variants: Vec<crate::postprocess::somatic::SomaticVariant>,
) -> Vec<crate::postprocess::somatic::SomaticVariant> {
    let mut buckets: std::collections::BTreeMap<i32, crate::postprocess::somatic::SomaticVariant> =
        std::collections::BTreeMap::new();
    for sv in somatic_variants {
        match buckets.entry(sv.variant1.start_pos) {
            std::collections::btree_map::Entry::Vacant(e) => {
                e.insert(sv);
            }
            std::collections::btree_map::Entry::Occupied(mut e) => {
                let sv_af = crate::types::round_half_even(sv.variant1.frequency, 4);
                let cur_af = crate::types::round_half_even(e.get().variant1.frequency, 4);
                if sv_af > cur_af {
                    e.insert(sv);
                }
            }
        }
    }
    buckets
        .into_values()
        .filter(|sv| !crate::output::vcf::perl_row_skipped(&sv.variant1.ref_allele))
        .collect()
}

/// ports CigarParser.java:1289-1301
///
/// Single-owner rule for a splice junction key (`A-B`, where `A` is the last
/// pre-`N` reference base): the tile whose active interval contains `A` owns
/// the junction; edge-clamped anchors outside the region go to the first
/// (resp. last) tile. Ownership is over the anchor `A` — aligners never emit
/// a leading `N`, so a junction always has one.
fn tile_owns_splice_key(
    tile: &crate::regions::Tile,
    reg_start: crate::types::Pos1,
    reg_end: crate::types::Pos1,
    key: &str,
) -> bool {
    // Fail-open: keys are internally generated so an unparseable key is
    // unreachable; emit rather than silently drop.
    let Some(anchor) = key.split_once('-').and_then(|(a, _)| a.parse::<u32>().ok()) else {
        return true;
    };
    if tile.start.0 <= anchor && anchor <= tile.end.0 {
        return true;
    }
    if anchor < reg_start.0 && tile.start == reg_start {
        return true;
    }
    if anchor > reg_end.0 && tile.end == reg_end {
        return true;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variations::Variant;

    fn sv(pos: i32, freq: f64, ref_allele: &str) -> crate::postprocess::somatic::SomaticVariant {
        let mut v = Variant::new();
        v.start_pos = pos;
        v.frequency = freq;
        v.ref_allele = ref_allele.to_string();
        crate::postprocess::somatic::SomaticVariant {
            variant1: v,
            variant2: None,
            status: String::new(),
            sv1: String::new(),
            sv2: String::new(),
            end_is_variant2: false,
        }
    }

    #[test]
    fn test_empty_ref_winner_suppresses_position() {
        // ports var2vcf_paired.pl:128-136 — the empty-REF row has the higher AF,
        // so it wins the bucket contest and `next unless ($ref)` then drops it:
        // NOTHING is emitted for that position (the real-allele runner-up must
        // NOT surface, unlike a pre-bucket guard placement).
        let rows = vec![sv(100, 0.5, ""), sv(100, 0.3, "A")];
        assert!(select_paired_rows(rows).is_empty());

        // Real-allele winner still passes through.
        let rows = vec![sv(100, 0.5, "G"), sv(100, 0.3, "A")];
        assert_eq!(select_paired_rows(rows).len(), 1);

        // Perl-falsy "0" REF is dropped too.
        let rows = vec![sv(200, 0.9, "0")];
        assert!(select_paired_rows(rows).is_empty());

        // 4dp .var table AF ties preserve the first-seen row (var2vcf_paired.pl:128
        // sorts `$d->[14]` after 4dp rounding; e.g., combineAnalysis unrounded
        // 2/324 = 0.0061728395... vs ToVarsBuilder 0.0062 both format as 0.0062).
        let rows = vec![
            sv(300, 2.0 / 324.0, "TAGAGAGTCCCAG"),
            sv(300, 0.0062, "TAGAGAG"),
        ];
        let selected = select_paired_rows(rows);
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].variant1.ref_allele, "TAGAGAGTCCCAG");
    }

    #[test]
    fn amplicon_covering_position_matches_expansion_and_scales() {
        use crate::regions::{AmpliconInfo, Region};
        use crate::types::Pos1;

        fn amp_region(chr: &str, start: u32, end: u32, is: u32, ie: u32) -> Region {
            Region {
                chr: chr.to_string(),
                start: Pos1(start),
                end: Pos1(end),
                gene: String::new(),
                amplicon: Some(AmpliconInfo {
                    insert_start: Pos1(is),
                    insert_end: Pos1(ie),
                }),
            }
        }

        // Synthetic overlapping intervals plus one region without amplicon info.
        let segment = vec![
            amp_region("chr1", 90, 210, 100, 200),
            amp_region("chr1", 140, 260, 150, 250),
            Region {
                chr: "chr1".to_string(),
                start: Pos1(300),
                end: Pos1(400),
                gene: String::new(),
                amplicon: None,
            },
        ];

        // Old logic: per-base expansion (AmpliconMode.java:57-60).
        let mut expanded: std::collections::BTreeMap<i32, Vec<usize>> =
            std::collections::BTreeMap::new();
        for (n, region) in segment.iter().enumerate() {
            if let Some(ref amp) = region.amplicon {
                for p in amp.insert_start.0 as i32..=amp.insert_end.0 as i32 {
                    expanded.entry(p).or_default().push(n);
                }
            }
        }
        // New logic agrees at every expanded position, in the same order.
        for (p, want) in &expanded {
            let got: Vec<usize> = amplicons_covering_position(&segment, *p)
                .iter()
                .map(|(n, _)| *n)
                .collect();
            assert_eq!(&got, want, "covering set differs at {p}");
        }
        // ...and reports empty outside all intervals.
        assert!(amplicons_covering_position(&segment, 99).is_empty());
        assert!(amplicons_covering_position(&segment, 251).is_empty());

        // 250 Mbp insert interval: overlap tests answer in milliseconds with no
        // per-base allocation (the old expansion would materialize 250M entries).
        let huge = vec![amp_region("chr1", 1, 250_000_000, 1, 250_000_000)];
        let t0 = std::time::Instant::now();
        assert_eq!(amplicons_covering_position(&huge, 1).len(), 1);
        assert_eq!(amplicons_covering_position(&huge, 125_000_000).len(), 1);
        assert_eq!(amplicons_covering_position(&huge, 250_000_000).len(), 1);
        assert!(amplicons_covering_position(&huge, 250_000_001).is_empty());
        assert!(
            t0.elapsed() < std::time::Duration::from_secs(5),
            "250 Mbp overlap test took too long"
        );
    }

    fn t(start: u32, end: u32) -> crate::regions::Tile {
        crate::regions::Tile {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(start),
            end: crate::types::Pos1(end),
            gene: String::new(),
            halo_start: crate::types::Pos1(start),
            halo_end: crate::types::Pos1(end),
            amplicon: None,
            original_start: crate::types::Pos1(start),
            original_end: crate::types::Pos1(end),
        }
    }

    fn p(v: u32) -> crate::types::Pos1 {
        crate::types::Pos1(v)
    }

    #[test]
    fn splice_key_owned_by_exactly_one_tile() {
        // tiles [100,399] and [400,699] of region [100,699]; junction "399-450"
        assert!(tile_owns_splice_key(
            &t(100, 399),
            p(100),
            p(699),
            "399-450"
        ));
        assert!(!tile_owns_splice_key(
            &t(400, 699),
            p(100),
            p(699),
            "399-450"
        ));
        // boundary: anchor 400 belongs to tile 2 only
        assert!(!tile_owns_splice_key(
            &t(100, 399),
            p(100),
            p(699),
            "400-450"
        ));
        assert!(tile_owns_splice_key(
            &t(400, 699),
            p(100),
            p(699),
            "400-450"
        ));
        // edge clamp: Java counts A = reg_start - 1; first tile owns it
        assert!(tile_owns_splice_key(&t(100, 399), p(100), p(699), "99-150"));
        assert!(!tile_owns_splice_key(
            &t(400, 699),
            p(100),
            p(699),
            "99-150"
        ));
        // fail-open: unparseable key is emitted (never silently dropped)
        assert!(tile_owns_splice_key(&t(100, 399), p(100), p(699), "junk"));
    }

    // --- Task 2 (red): pileup-level splice tile-emission tests ---
    //
    // Per-tile spec mirror of the `run_simple` splice branch plus the Task-1
    // ownership gate: pile up each tile's halo-queried records (via
    // `run_pipeline_on_records`, exactly as production does) and emit one
    // `sample\tchr\tjunction\tcount` line per junction OWNED by the tile.
    // Task 3 implements `tile_owns_splice_key` and applies this same gate in
    // production; these tests must green with no test edits.

    use noodles::sam::alignment::record_buf::RecordBuf;

    fn splice_test_config() -> crate::config::Config {
        let mut config = crate::config::Config::default();
        config.realign = 0;
        config.splice = true;
        config
    }

    fn splice_test_reference() -> crate::reference::Reference {
        let mut reference = crate::reference::Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 1000]);
        reference
    }

    fn splice_test_chr_lens() -> std::collections::HashMap<String, usize> {
        [("chr1".to_string(), 1000usize)].into_iter().collect()
    }

    fn test_record_buf(
        name: &str,
        start: usize,
        cigar_ops: &[(usize, noodles::sam::alignment::record::cigar::op::Kind)],
        seq: Vec<u8>,
    ) -> RecordBuf {
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record_buf::Cigar;
        let cigar: Cigar = cigar_ops
            .iter()
            .map(|&(len, kind)| Op::new(kind, len))
            .collect();
        let quals = vec![40u8; seq.len()];
        RecordBuf::builder()
            .set_name(name)
            .set_flags(noodles::sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(cigar)
            .set_sequence(seq.into())
            .set_quality_scores(quals.into())
            .build()
    }

    /// Mirror of `crates/vardictrust-core/tests/common::indexed_bam`: writes `path`
    /// plus its `.bam.bai` index.
    fn write_test_bam(path: &std::path::Path, chr: &str, len: usize, recs: &[RecordBuf]) {
        use noodles::bam::bai;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::{
            bin::Chunk, index::LinearIndex,
        };
        use noodles::csi::binning_index::Indexer;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;
        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(len.try_into().unwrap())
            .build()
            .unwrap();
        let hd = Header::builder().add_reference_sequence(chr, rs).build();
        let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
        let mut w = noodles::bam::io::Writer::from(bgzf);
        w.write_alignment_header(&hd).unwrap();
        let mut indexer = Indexer::<LinearIndex>::new(14, 5).expect("valid index parameters");
        for rec in recs {
            let start = w.get_ref().virtual_position();
            w.write_alignment_record(&hd, rec).unwrap();
            let end = w.get_ref().virtual_position();
            let s = rec.alignment_start().expect("aligned record");
            let e = Position::new(s.get() + rec.cigar().alignment_span().saturating_sub(1).max(1))
                .unwrap_or(Position::MIN);
            indexer
                .add_record(Some((0, s, e, true)), Chunk::new(start, end))
                .map_err(std::io::Error::other)
                .unwrap();
        }
        w.try_finish().unwrap();
        let mut bai_path = std::path::PathBuf::from(path);
        bai_path.set_extension("bam.bai");
        bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
            .write_index(&indexer.build(1))
            .unwrap();
    }

    /// Desired per-tile splice emission: one line per OWNED junction, in sorted
    /// key order, mirroring the `run_simple` splice-branch line format.
    fn gated_splice_lines(
        tile: &crate::regions::Tile,
        region_start: crate::types::Pos1,
        region_end: crate::types::Pos1,
        sample: &str,
        splice_count: &std::collections::HashMap<String, [i32; 1]>,
    ) -> Vec<String> {
        let mut introns: Vec<&String> = splice_count.keys().collect();
        introns.sort();
        let mut lines = Vec::new();
        for intron in introns {
            if !tile_owns_splice_key(tile, region_start, region_end, intron) {
                continue;
            }
            let count = splice_count.get(intron).map(|a| a[0]).unwrap_or(0);
            lines.push(format!("{}\t{}\t{}\t{}", sample, tile.chr, intron, count));
        }
        lines
    }

    fn pileup_tile(
        bam: &str,
        tile: &crate::regions::Tile,
        reference: &mut crate::reference::Reference,
        config: &crate::config::Config,
        chr_lens: &std::collections::HashMap<String, usize>,
    ) -> crate::variations::VariationData {
        // Pile up over the HALO bounds, not the active bounds: a junction
        // anchored in the tile's active region is recorded when the CIGAR loop
        // reaches its N op, which sits one base past the anchor (`start` after
        // the pre-N match block). Active-region pileup `break`s at
        // `start > region.end` before reaching an N whose anchor is the tile's
        // last base, so the owning tile would never count its own junction.
        // Halo-wide pileup counts every candidate; the ownership gate in
        // `gated_splice_lines` restores single emission.
        let pileup_region = crate::regions::Region {
            chr: tile.chr.clone(),
            start: tile.halo_start,
            end: tile.halo_end,
            gene: tile.gene.clone(),
            amplicon: None,
        };
        let owned = bam.to_string();
        run_pipeline_on_records(
            std::slice::from_ref(&owned),
            &tile.chr,
            tile.halo_start.to_zero_based(),
            crate::types::Pos0(tile.halo_end.0),
            None,
            reference,
            &pileup_region,
            PipelineTail::Stop,
            config,
            chr_lens,
            0,
            None,
        )
        .unwrap()
    }

    #[test]
    fn splice_tile_emission_halo_junction_emitted_once_with_full_count() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // Region chr1:[100,699], tiles of 300 + halo 50:
        //   tile 1 active [100,399] halo [100,449];
        //   tile 2 active [400,699] halo [350,699].
        // One spliced read at 390 (10M 50N 10M): junction "399-449" anchors in
        // tile 1 (399) while the body spans 390..459, so tile 2 counts a
        // junction anchored outside its active region, visible to it only
        // through its halo-queried read. Desired: exactly one line (from the
        // owning tile 1) with the whole-region count.
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(100),
            end: crate::types::Pos1(699),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let tiles = crate::regions::tile(&region, 300, 50);
        assert_eq!(tiles.len(), 2);
        assert_eq!(tiles[0].halo_start, crate::types::Pos1(100));
        assert_eq!(tiles[0].halo_end, crate::types::Pos1(449));
        assert_eq!(tiles[1].halo_start, crate::types::Pos1(350));
        assert_eq!(tiles[1].halo_end, crate::types::Pos1(699));

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("splice.bam");
        write_test_bam(
            &bam_path,
            "chr1",
            1000,
            &[test_record_buf(
                "spliced1",
                390,
                &[(10, Kind::Match), (50, Kind::Skip), (10, Kind::Match)],
                vec![b'A'; 20],
            )],
        );
        let bam = bam_path.to_string_lossy().into_owned();
        let paths = [bam.clone()];

        let config = splice_test_config();
        let mut reference = splice_test_reference();
        let chr_lens = splice_test_chr_lens();

        // Whole-region pileup pins the expected count (oracle cross-check).
        let whole = run_pipeline_on_records(
            &paths,
            "chr1",
            crate::types::Pos1(100).to_zero_based(),
            crate::types::Pos0(699),
            None,
            &mut reference,
            &region,
            PipelineTail::Stop,
            &config,
            &chr_lens,
            0,
            None,
        )
        .unwrap();
        let whole_count = whole.splice_count.get("399-449").map(|a| a[0]).unwrap_or(0);
        assert_eq!(whole_count, 1, "oracle: whole-region junction count");

        let mut all_lines = Vec::new();
        let mut per_tile_counts = Vec::new();
        for tile in &tiles {
            let var_data = pileup_tile(&bam, tile, &mut reference, &config, &chr_lens);
            // Premise: the halo overlap makes the read visible to BOTH tiles'
            // pileups, so without the ownership gate the junction would be
            // emitted twice (unowned duplicate).
            per_tile_counts.push(var_data.splice_count.get("399-449").map(|a| a[0]));
            all_lines.extend(gated_splice_lines(
                tile,
                crate::types::Pos1(100),
                crate::types::Pos1(699),
                &config.sample,
                &var_data.splice_count,
            ));
        }

        // Java-blessed expected splice line (sample\tchr\tjunction\tcount).
        assert_eq!(
            per_tile_counts,
            vec![Some(1), Some(1)],
            "halo overlap: both tiles count the junction (duplicates without the gate)"
        );
        assert_eq!(all_lines, vec!["SAMPLE\tchr1\t399-449\t1".to_string()]);
    }

    #[test]
    fn splice_tile_emission_empty_tile_zero_lines() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // Tile chr1:[400,699] halo [400,699] over a BAM with reads but NO
        // spliced (N-op) reads: pileup yields an empty splice_count, so
        // emission must produce zero lines without panicking on the empty map.
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(400),
            end: crate::types::Pos1(699),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let tiles = crate::regions::tile(&region, 300, 50);
        assert_eq!(tiles.len(), 1);

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("nosplice.bam");
        write_test_bam(
            &bam_path,
            "chr1",
            1000,
            &[
                test_record_buf("plain1", 500, &[(50, Kind::Match)], vec![b'A'; 50]),
                test_record_buf("plain2", 560, &[(50, Kind::Match)], vec![b'A'; 50]),
            ],
        );
        let bam = bam_path.to_string_lossy().into_owned();

        let config = splice_test_config();
        let mut reference = splice_test_reference();
        let chr_lens = splice_test_chr_lens();

        let var_data = pileup_tile(&bam, &tiles[0], &mut reference, &config, &chr_lens);
        assert!(
            var_data.splice_count.is_empty(),
            "no N-op reads so splice_count must be empty"
        );
        let lines = gated_splice_lines(
            &tiles[0],
            crate::types::Pos1(400),
            crate::types::Pos1(699),
            &config.sample,
            &var_data.splice_count,
        );
        assert!(
            lines.is_empty(),
            "empty splice_count must emit zero lines, got: {:?}",
            lines
        );
    }

    #[test]
    fn g3_mate_window_merge_emits_out_of_tile_key_from_both_tiles() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // Task 6 (G3, sibling audit rows 16/46-47): mate-window
        // `partial_pipeline` merges (realign.rs:4012-4063,
        // sv.rs:3464-3478) plant out-of-active-range keys that pass the
        // whole-region builder gate, so the merging tile emits rows owned
        // by another tile -> duplicates in intermediate tables (#1,#3)
        // and somatic VCF (#4).
        //
        // Region chr1:[100,699], tiles of 300 + halo 50:
        //   tile 1 active [100,399]; tile 2 active [400,699].
        // SNV reads at q=200 are owned by tile 1. Tile 2's mate-window
        // merge over [150,250] (the realign.rs:1959-1971 /
        // sv.rs:1374-1383 call shape for a mate breakpoint near 200)
        // plants key 200 in tile 2's maps; the whole-region builder gate
        // (lib.rs:440-446, tovars.rs:116-121) then emits 200 from BOTH
        // tiles while the whole-region run emits it once.
        // RED while G3 is unaddressed: tile 2 emits the tile-1-owned key.
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(100),
            end: crate::types::Pos1(699),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let tiles = crate::regions::tile(&region, 300, 50);
        assert_eq!(tiles.len(), 2);

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("g3.bam");
        // 3 reads, 20M at 191 (1-based), alt C at position 200 (ref A).
        let mut seq = vec![b'A'; 20];
        seq[9] = b'C';
        write_test_bam(
            &bam_path,
            "chr1",
            1000,
            &[
                test_record_buf("snv1", 191, &[(20, Kind::Match)], seq.clone()),
                test_record_buf("snv2", 191, &[(20, Kind::Match)], seq.clone()),
                test_record_buf("snv3", 191, &[(20, Kind::Match)], seq.clone()),
            ],
        );
        let bam = bam_path.to_string_lossy().into_owned();
        let paths = [bam.clone()];

        let config = crate::config::Config::default();
        let chr_lens = splice_test_chr_lens();

        // Per-tile pileup over the halo window with the ACTIVE region (the
        // production shape).
        let pileup_active = |tile: &crate::regions::Tile,
                             reference: &mut crate::reference::Reference|
         -> crate::variations::VariationData {
            let pileup_region = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.start,
                end: tile.end,
                gene: tile.gene.clone(),
                amplicon: None,
            };
            run_pipeline_on_records(
                &paths,
                &tile.chr,
                tile.halo_start.to_zero_based(),
                crate::types::Pos0(tile.halo_end.0),
                None,
                reference,
                &pileup_region,
                PipelineTail::Stop,
                &config,
                &chr_lens,
                0,
                None,
            )
            .unwrap()
        };

        // Emission through the production whole-region builder gate
        // (lib.rs:440-460): builder_region spans the whole BED region.
        let emit = |var_data: &crate::variations::VariationData,
                    reference: &crate::reference::Reference|
         -> Vec<i32> {
            let builder = crate::tovars::ToVarsBuilder::new(reference, &region, &chr_lens, &config);
            let aligned = builder.build_aligned_variants(var_data);
            crate::postprocess::simple::simple_post_process(
                &config,
                &region,
                &aligned,
                &var_data.splice,
            )
            .iter()
            .map(|v| v.start_pos)
            .collect()
        };

        // Tile 1 (owner of q=200).
        let mut ref1 = splice_test_reference();
        let var1 = pileup_active(&tiles[0], &mut ref1);
        assert!(
            var1.non_insertion_variants.contains_key(&200),
            "owner tile 1 must pile up key 200"
        );

        // Tile 2 (merger): its own window holds no reads; the SV/mate path
        // merges the mate window covering q=200 into its maps.
        let mut ref2 = splice_test_reference();
        let mut var2 = pileup_active(&tiles[1], &mut ref2);
        assert!(
            !var2.non_insertion_variants.contains_key(&200),
            "merging tile 2 must not pile up key 200 on its own"
        );
        crate::sv::partial_pipeline(
            &bam_path, &mut var2, &mut ref2, "chr1", 150, 250, &config, &chr_lens,
        )
        .unwrap();
        assert!(
            var2.non_insertion_variants.contains_key(&200),
            "G3 mechanism: mate-window merge planted out-of-active key 200 in tile 2's maps"
        );

        // Whole-region oracle: single pileup + single emission.
        let mut refw = splice_test_reference();
        let varw = run_pipeline_on_records(
            &paths,
            "chr1",
            crate::types::Pos1(100).to_zero_based(),
            crate::types::Pos0(699),
            None,
            &mut refw,
            &region,
            PipelineTail::Stop,
            &config,
            &chr_lens,
            0,
            None,
        )
        .unwrap();

        let mut pos1 = emit(&var1, &ref1);
        let mut pos2 = emit(&var2, &ref2);
        let mut posw = emit(&varw, &refw);
        pos1.sort();
        pos2.sort();
        posw.sort();

        assert_eq!(
            posw.iter().filter(|p| **p == 200).count(),
            1,
            "oracle: whole-region run emits 200 exactly once, got {:?}",
            posw
        );
        assert!(
            !pos2.contains(&200),
            "G3 REPRO: merging tile 2 emits tile-1-owned position 200 \
             (tiled: {pos1:?} + {pos2:?}, whole: {posw:?})",
        );
    }

    #[test]
    fn g3_merge_planted_requires_absence_from_both_maps() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // Thread 3: `merge_planted` attribution is per position but novelty
        // was checked per map, so a native non-insertion entry at `k` plus a
        // mate-window merge adding the FIRST insertion at `k` marked the
        // whole position planted — and the ToVarsBuilder ownership gate then
        // dropped the native entry too at out-of-owning-range positions
        // (native SV outside the active interval is by design: realignlgdel
        // scans soft-clips up to 5000bp beyond the pileup region while the
        // SV gate keeps out-of-region SV rows). A position is planted only
        // when absent from BOTH maps pre-merge.
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(100),
            end: crate::types::Pos1(699),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let tiles = crate::regions::tile(&region, 300, 50);
        assert_eq!(tiles.len(), 2);

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("g3_planted.bam");
        // Insertion reads at 191 (1-based): 10M + 3I + 10M. The mate window
        // [150,250] covers them; tile 2's halo does not.
        let mut seq = vec![b'A'; 10];
        seq.extend_from_slice(b"TTT");
        seq.extend_from_slice(&[b'A'; 10]);
        write_test_bam(
            &bam_path,
            "chr1",
            1000,
            &[
                test_record_buf(
                    "ins1",
                    191,
                    &[(10, Kind::Match), (3, Kind::Insertion), (10, Kind::Match)],
                    seq.clone(),
                ),
                test_record_buf(
                    "ins2",
                    191,
                    &[(10, Kind::Match), (3, Kind::Insertion), (10, Kind::Match)],
                    seq.clone(),
                ),
                test_record_buf(
                    "ins3",
                    191,
                    &[(10, Kind::Match), (3, Kind::Insertion), (10, Kind::Match)],
                    seq.clone(),
                ),
            ],
        );

        let config = crate::config::Config::default();
        let chr_lens = splice_test_chr_lens();
        let merge = |var_data: &mut crate::variations::VariationData,
                     reference: &mut crate::reference::Reference| {
            crate::sv::partial_pipeline(
                &bam_path, var_data, reference, "chr1", 150, 250, &config, &chr_lens,
            )
            .unwrap();
        };

        // Phase 1 (discovery): empty maps -> every merged key is planted.
        let mut ref1 = splice_test_reference();
        let mut var1 = crate::variations::VariationData::new();
        merge(&mut var1, &mut ref1);
        let ins_keys: Vec<i32> = var1.insertion_variants.keys().copied().collect();
        assert!(
            !ins_keys.is_empty(),
            "mate-window merge must plant insertion keys, got non_ins={:?} ins={:?}",
            var1.non_insertion_variants.keys().collect::<Vec<_>>(),
            ins_keys,
        );
        for k in var1
            .non_insertion_variants
            .keys()
            .chain(var1.insertion_variants.keys())
        {
            assert!(
                var1.merge_planted.contains(k),
                "truly-new key {k} must be planted"
            );
        }
        let target = ins_keys[0];

        // Phase 2 (bug scenario): native non-insertion entry at the insertion
        // key (stands in for an SV-path native entry outside the active
        // interval) + mate merge adding the first insertion there.
        let mut ref2 = splice_test_reference();
        let mut var2 = crate::variations::VariationData::new();
        let mut vm = crate::variations::VariationMap::new();
        let mut v = crate::variations::Variation::new();
        v.cnt = 3;
        v.fwd = 2;
        v.rev = 1;
        vm.map.insert("C".to_string(), v);
        var2.non_insertion_variants.insert(target, vm);
        var2.ref_coverage.insert(target, 10);
        merge(&mut var2, &mut ref2);
        assert!(
            var2.insertion_variants.contains_key(&target),
            "mate merge must add the first insertion at {target}"
        );
        assert!(
            !var2.merge_planted.contains(&target),
            "native position {target} must not be planted by a first-insertion merge"
        );

        // Phase 3 (gate): with the owning range excluding the native
        // position, the builder still emits it (pre-fix it was dropped).
        var2.owning_range = Some((tiles[1].start, tiles[1].end));
        let builder = crate::tovars::ToVarsBuilder::new(&ref2, &region, &chr_lens, &config);
        let aligned = builder.build_aligned_variants(&var2);
        assert!(
            aligned.contains_key(&target),
            "native position {target} outside the owning range must survive the builder gate"
        );
    }

    #[test]
    fn test_effective_tile_size_untiled_up_to_1mbp_and_100kbp_above_1mbp() {
        let reg_980k = crate::regions::Region {
            chr: "chr20".to_string(),
            start: crate::types::Pos1(10_020_001),
            end: crate::types::Pos1(11_000_000),
            gene: "G".to_string(),
            amplicon: None,
        };
        assert_eq!(effective_tile_size(&reg_980k, None), 980_000);
        assert_eq!(
            crate::regions::tile(&reg_980k, effective_tile_size(&reg_980k, None), 500).len(),
            1
        );

        let reg_1m = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(1_000_000),
            gene: "G".to_string(),
            amplicon: None,
        };
        assert_eq!(effective_tile_size(&reg_1m, None), 1_000_000);
        assert_eq!(
            crate::regions::tile(&reg_1m, effective_tile_size(&reg_1m, None), 500).len(),
            1
        );

        let reg_10m = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(10_000_000),
            gene: "G".to_string(),
            amplicon: None,
        };
        assert_eq!(effective_tile_size(&reg_10m, None), DEFAULT_TILE_SIZE);
        assert_eq!(
            crate::regions::tile(&reg_10m, effective_tile_size(&reg_10m, None), 500).len(),
            100
        );

        // Explicit --tile-size always overrides auto heuristic.
        assert_eq!(effective_tile_size(&reg_980k, Some(50_000)), 50_000);
        assert_eq!(effective_tile_size(&reg_10m, Some(2_000_000)), 2_000_000);
    }

    #[test]
    fn write_output_line_propagates_broken_pipe_error() {
        struct BrokenPipeWriter;
        impl std::io::Write for BrokenPipeWriter {
            fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::new(
                    std::io::ErrorKind::BrokenPipe,
                    "broken pipe",
                ))
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }

        let mut writer = BrokenPipeWriter;
        let err = write_output_line(&mut writer, "chr1\t100").unwrap_err();
        match err {
            VardictError::Io(e) => assert_eq!(e.kind(), std::io::ErrorKind::BrokenPipe),
            other => panic!("expected VardictError::Io(BrokenPipe), got {other:?}"),
        }

        let mut cfg = crate::config::Config::default();
        assert!(matches!(
            create_output_writer(&cfg),
            OutputDestination::Stdout(_)
        ));
        cfg.debug_print = Some("ERR".to_string());
        assert!(matches!(
            create_output_writer(&cfg),
            OutputDestination::Stderr(_)
        ));
    }

    #[test]
    fn somatic_rejects_empty_normal_bam_spec() {
        for bam_spec in ["tumor.bam|", "tumor.bam|::"] {
            let cfg = crate::config::Config {
                bam: bam_spec.to_string(),
                region: Some("chr1:1-100".to_string()),
                ..Default::default()
            };
            let mut out = Vec::new();
            let err = run_with_writer(&cfg, &mut out).unwrap_err();
            match err {
                VardictError::Parse(msg) => {
                    assert!(
                        msg.contains("No normal BAM file specified"),
                        "expected 'No normal BAM file specified' for {bam_spec:?}, got {msg:?}"
                    );
                }
                other => panic!("expected VardictError::Parse for {bam_spec:?}, got {other:?}"),
            }
        }

        let simple_cfg = crate::config::Config {
            bam: "::".to_string(),
            region: Some("chr1:1-100".to_string()),
            ..Default::default()
        };
        let mut out = Vec::new();
        let err = run_with_writer(&simple_cfg, &mut out).unwrap_err();
        match err {
            VardictError::Parse(msg) => {
                assert!(
                    msg.contains("No BAM files specified"),
                    "expected 'No BAM files specified', got {msg:?}"
                );
            }
            other => panic!("expected VardictError::Parse, got {other:?}"),
        }
    }

    #[test]
    fn somatic_bam2_sv_pipeline_inherits_bam1_splice_set() {
        // SomaticMode.java:104-122 passes the same mutable `Set<String> splice` instance
        // into both `initialScope1` and `initialScope2`, so splice junctions discovered in
        // BAM1 suppress < 250 kbp discordant-pair deletions during BAM2's
        // `StructuralVariantsProcessor.findDELdisc` (`sv.rs:2396`).
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::Cigar;
        use noodles::sam::alignment::RecordBuf;

        let dir = tempfile::TempDir::new().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        let ref_bases = "A".repeat(2000);
        std::fs::write(&fasta_path, format!(">chr1\n{ref_bases}\n")).unwrap();
        std::fs::write(&fai_path, "chr1\t2000\t6\t2000\t2001\n").unwrap();

        let tumor_bam = dir.path().join("tumor.bam");
        let normal_bam = dir.path().join("normal.bam");

        let spliced_cigar: Cigar = [
            Op::new(Kind::Match, 10),
            Op::new(Kind::Skip, 50),
            Op::new(Kind::Match, 10),
        ]
        .into_iter()
        .collect();
        let spliced_read = RecordBuf::builder()
            .set_name("t_splice")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(50usize).unwrap())
            .set_mapping_quality(noodles::sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(spliced_cigar)
            .set_sequence(vec![b'A'; 20].into())
            .set_quality_scores(vec![40u8; 20].into())
            .build();
        write_test_bam(&tumor_bam, "chr1", 2000, &[spliced_read]);

        let mut normal_reads = Vec::new();
        for i in 0..7usize {
            let cigar: Cigar = [Op::new(Kind::Match, 20)].into_iter().collect();
            let rec = RecordBuf::builder()
                .set_name(format!("n_disc_{i}"))
                .set_flags(Flags::SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED)
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(200 + i).unwrap())
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(60).unwrap(),
                )
                .set_cigar(cigar)
                .set_mate_reference_sequence_id(0)
                .set_mate_alignment_start(Position::try_from(1500usize).unwrap())
                .set_template_length(1320)
                .set_sequence(vec![b'A'; 20].into())
                .set_quality_scores(vec![40u8; 20].into())
                .build();
            normal_reads.push(rec);
        }
        write_test_bam(&normal_bam, "chr1", 2000, &normal_reads);

        let cfg = crate::config::Config {
            bam: format!(
                "{}|{}",
                tumor_bam.to_string_lossy(),
                normal_bam.to_string_lossy()
            ),
            reference: fasta_path.to_string_lossy().into_owned(),
            region: Some("chr1:1-900".to_string()),
            emit_intermediate: true,
            min_reads: 2,
            freq: 0.01,
            ..Default::default()
        };
        let mut out = Vec::new();
        run_with_writer(&cfg, &mut out).unwrap();
        let output = String::from_utf8(out).unwrap();
        assert!(
            !output.contains("<DEL>"),
            "BAM2 < 250 kbp discordant DEL must be suppressed when BAM1 observed a splice junction, got:\n{output}"
        );
    }

    fn vcf_record_keys(vcf_text: &str) -> Vec<String> {
        vcf_text
            .lines()
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                format!("{}:{}:{}:{}", cols[0], cols[1], cols[3], cols[4])
            })
            .collect()
    }

    #[test]
    fn test_issue_f_whole_stream_vcf_accumulation_and_reorder() {
        let workspace_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(|p| p.parent())
            .expect("workspace root");
        let micro_ref = workspace_root
            .join("fixtures/micro/micro.ref.fa")
            .to_string_lossy()
            .into_owned();
        let micro_bam = workspace_root
            .join("fixtures/micro/micro.bam")
            .to_string_lossy()
            .into_owned();
        let som_bam = format!("{micro_bam}|{micro_bam}");

        let dir = tempfile::TempDir::new().unwrap();
        let interleaved_bed = dir.path().join("interleaved.bed");
        std::fs::write(
            &interleaved_bed,
            "chr1\t40\t260\tG1\nchr2\t40\t260\tG2\nchr1\t50\t250\tG3\n",
        )
        .unwrap();
        let grouped_bed = dir.path().join("grouped.bed");
        std::fs::write(
            &grouped_bed,
            "chr1\t40\t260\tG1\nchr1\t50\t250\tG3\nchr2\t40\t260\tG2\n",
        )
        .unwrap();
        let chr21_bed = dir.path().join("chr21.bed");
        std::fs::write(&chr21_bed, "chr2\t40\t260\tG2\nchr1\t40\t260\tG1\n").unwrap();

        let run_cfg = |bam: &str, bed: &std::path::Path, emit_intermediate: bool| -> String {
            let cfg = crate::config::Config {
                reference: micro_ref.clone(),
                bam: bam.to_string(),
                sample: "S".to_string(),
                bed: Some(bed.to_string_lossy().into_owned()),
                col_chr: 0,
                col_start: 1,
                col_end: 2,
                col_seg_start: 1,
                col_seg_end: 2,
                col_gene: Some(3),
                column_for_chromosome: Some(0),
                zero_based: true,
                fisher: true,
                emit_intermediate,
                threads: 1,
                ..Default::default()
            };
            let mut out = Vec::new();
            run_with_writer(&cfg, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        };

        let want = vec!["chr1:101:T:A".to_string(), "chr2:150:GCGG:G".to_string()];

        for (name, bed) in [
            ("interleaved", &interleaved_bed),
            ("grouped", &grouped_bed),
            ("chr21", &chr21_bed),
        ] {
            let simple_vcf = run_cfg(&micro_bam, bed, false);
            assert_eq!(
                vcf_record_keys(&simple_vcf),
                want,
                "simple VCF {name} must deduplicate per (chr, start) and emit in perl reorder() order"
            );
            let som_vcf = run_cfg(&som_bam, bed, false);
            assert_eq!(
                vcf_record_keys(&som_vcf),
                want,
                "somatic VCF {name} must deduplicate per (chr, start) and emit in perl reorder() order"
            );
        }

        // Intermediate .var mode must preserve BED stream order without deduplication
        let interleaved_var = run_cfg(&micro_bam, &interleaved_bed, true);
        let var_chrs: Vec<&str> = interleaved_var
            .lines()
            .filter_map(|l| l.split('\t').nth(2))
            .collect();
        assert_eq!(
            var_chrs,
            vec!["chr1", "chr2", "chr1"],
            ".var output must remain in BED region order"
        );

        // Amplicon mode: nested/disjoint segments on chr1 plus chr2 listed before chr1
        let amp_bed = dir.path().join("amp_nested_2chr.bed");
        std::fs::write(
            &amp_bed,
            concat!(
                "chr2\t40\t260\tG2\t0\t+\t50\t250\n",
                "chr1\t40\t260\tA\t0\t+\t79\t125\n",
                "chr1\t70\t100\tB\t0\t+\t81\t90\n",
                "chr1\t84\t130\tC\t0\t+\t94\t120\n",
            ),
        )
        .unwrap();
        let amp_cfg = crate::config::Config {
            reference: micro_ref.clone(),
            bam: micro_bam.clone(),
            sample: "S".to_string(),
            bed: Some(amp_bed.to_string_lossy().into_owned()),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_gene: Some(3),
            col_seg_start: 6,
            col_seg_end: 7,
            amplicon: Some("200:0.1".to_string()),
            zero_based: true,
            fisher: true,
            emit_intermediate: false,
            threads: 1,
            ..Default::default()
        };
        let mut amp_out = Vec::new();
        run_with_writer(&amp_cfg, &mut amp_out).unwrap();
        let amp_vcf = String::from_utf8(amp_out).unwrap();
        assert_eq!(
            vcf_record_keys(&amp_vcf),
            want,
            "amplicon VCF must deduplicate across nested segments and order chr1 before chr2"
        );
    }

    #[test]
    fn test_issue_f_part_c_multi_megabase_tile_continuity() {
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::Cigar;
        use noodles::sam::alignment::RecordBuf;

        let dir = tempfile::TempDir::new().unwrap();
        let ref_len = 1_300_000usize;
        // Build deterministic non-homopolymer reference so realignment k-mers are unambiguous
        let motif = b"ACGTACGTAGCTAGCTTGCAACGT";
        let mut ref_bytes = Vec::with_capacity(ref_len);
        for i in 0..ref_len {
            ref_bytes.push(motif[i % motif.len()]);
        }
        // Anchor at 50,000 (1-based) = 'T'
        ref_bytes[49_999] = b'T';
        // At 100,020 (1-based) = 'T'
        ref_bytes[100_019] = b'T';

        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        let ref_str = String::from_utf8(ref_bytes.clone()).unwrap();
        std::fs::write(&fasta_path, format!(">chr1\n{ref_str}\n")).unwrap();
        std::fs::write(
            &fai_path,
            format!("chr1\t{ref_len}\t6\t{ref_len}\t{}\n", ref_len + 1),
        )
        .unwrap();

        let mk_read = |name: &str,
                       start_1based: usize,
                       flags: Flags,
                       cigar_ops: &[(usize, Kind)],
                       seq: Vec<u8>|
         -> RecordBuf {
            let cigar: Cigar = cigar_ops
                .iter()
                .map(|&(len, kind)| Op::new(kind, len))
                .collect();
            let quals = vec![40u8; seq.len()];
            RecordBuf::builder()
                .set_name(name)
                .set_flags(flags)
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start_1based).unwrap())
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(60).unwrap(),
                )
                .set_cigar(cigar)
                .set_sequence(seq.into())
                .set_quality_scores(quals.into())
                .build()
        };

        // 1. del_short / del_long BAMs:
        //    - 20 reads (10 fwd, 10 rev) with 50M 92D 50M starting at 49,951 (deletion at 50,001..50,092)
        //    - 10 reads (5 fwd, 5 rev) with 100M starting at 49,960 (ending at 50,059 inside the deletion)
        //    - del_long adds one 250M read at 749,900 (Tile 8)
        let mut del_short_reads = Vec::new();
        for i in 0..20usize {
            let flags = if i % 2 == 0 {
                Flags::empty()
            } else {
                Flags::REVERSE_COMPLEMENTED
            };
            let mut seq = Vec::with_capacity(100);
            seq.extend_from_slice(&ref_bytes[49_950..50_000]);
            seq.extend_from_slice(&ref_bytes[50_092..50_142]);
            del_short_reads.push(mk_read(
                &format!("del_{i}"),
                49_951,
                flags,
                &[(50, Kind::Match), (92, Kind::Deletion), (50, Kind::Match)],
                seq,
            ));
        }
        for i in 0..10usize {
            let flags = if i % 2 == 0 {
                Flags::empty()
            } else {
                Flags::REVERSE_COMPLEMENTED
            };
            let seq = ref_bytes[49_959..50_059].to_vec();
            del_short_reads.push(mk_read(
                &format!("ref_{i}"),
                49_960,
                flags,
                &[(100, Kind::Match)],
                seq,
            ));
        }
        del_short_reads.sort_by_key(|r| r.alignment_start().unwrap().get());

        let mut del_long_reads = del_short_reads.clone();
        del_long_reads.push(mk_read(
            "long_250",
            749_900,
            Flags::empty(),
            &[(250, Kind::Match)],
            ref_bytes[749_899..750_149].to_vec(),
        ));

        let del_short_bam = dir.path().join("del_short.bam");
        let del_long_bam = dir.path().join("del_long.bam");
        write_test_bam(&del_short_bam, "chr1", ref_len, &del_short_reads);
        write_test_bam(&del_long_bam, "chr1", ref_len, &del_long_reads);

        // Simple mode del_long: tiled (default 100kb on 1.2Mbp) must equal untiled (250Mbp)
        let run_simple_var =
            |bam: &std::path::Path, tile_size: Option<u32>, dedup: bool| -> String {
                let cfg = crate::config::Config {
                    reference: fasta_path.to_string_lossy().into_owned(),
                    bam: bam.to_string_lossy().into_owned(),
                    sample: "S".to_string(),
                    region: Some("chr1:1-1200000".to_string()),
                    emit_intermediate: true,
                    dedup,
                    tile_size,
                    threads: 2,
                    ..Default::default()
                };
                let mut out = Vec::new();
                run_with_writer(&cfg, &mut out).unwrap();
                String::from_utf8(out).unwrap()
            };

        let del_long_untiled = run_simple_var(&del_long_bam, Some(250_000_000), false);
        let del_long_tiled = run_simple_var(&del_long_bam, None, false);
        assert!(
            del_long_untiled.contains("\t-92/-92\t1.0000\t"),
            "untiled del_long must realign to hom -92/-92 AF=1.0000, got:\n{del_long_untiled}"
        );
        assert_eq!(
            del_long_tiled, del_long_untiled,
            "tiled del_long must match untiled del_long (Tile 1 realigndel sees Tile 8 max_read_length=250)"
        );

        // Control: del_short (max_read_length=100 everywhere) stays A/-92 AF=0.6667
        let del_short_tiled = run_simple_var(&del_short_bam, None, false);
        assert!(
            del_short_tiled.contains("\tA/-92\t0.6667\t"),
            "del_short control must remain het A/-92 AF=0.6667, got:\n{del_short_tiled}"
        );

        // Somatic mode del_long|del_short: BAM1 (del_long) finishes parent region with max_read_length=250
        // and seeds BAM2 (del_short), so tiled somatic must match untiled somatic!
        let run_som_var = |tile_size: Option<u32>| -> String {
            let cfg = crate::config::Config {
                reference: fasta_path.to_string_lossy().into_owned(),
                bam: format!(
                    "{}|{}",
                    del_long_bam.to_string_lossy(),
                    del_short_bam.to_string_lossy()
                ),
                sample: "T|N".to_string(),
                region: Some("chr1:1-1200000".to_string()),
                emit_intermediate: true,
                tile_size,
                threads: 2,
                ..Default::default()
            };
            let mut out = Vec::new();
            run_with_writer(&cfg, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        };
        let som_untiled = run_som_var(Some(250_000_000));
        let som_tiled = run_som_var(None);
        assert_eq!(
            som_tiled, som_untiled,
            "tiled somatic del_long|del_short must match untiled (BAM2 seeded with BAM1 region-final max_read_length=250)"
        );

        // 2. bound BAM (-t duplicate across tile boundary 100,000 -> 100,001):
        //    - Read A at 99,950 (50M, ends at 99,999 in Tile 1)
        //    - Read B at 99,950 (100M, spans 99,950..100,049 into Tile 2, carrying alt 'A' at 100,020)
        //    - 5 ref 'T' reads + 5 alt 'A' reads of 100M starting at 100,001..100,010
        let mut bound_reads = Vec::new();
        bound_reads.push(mk_read(
            "read_A",
            99_950,
            Flags::empty(),
            &[(50, Kind::Match)],
            ref_bytes[99_949..99_999].to_vec(),
        ));
        let mut seq_b = ref_bytes[99_949..100_049].to_vec();
        seq_b[100_019 - 99_949] = b'A';
        bound_reads.push(mk_read(
            "read_B_dup",
            99_950,
            Flags::empty(),
            &[(100, Kind::Match)],
            seq_b,
        ));
        for i in 0..5usize {
            let start = 100_001 + i;
            let seq_ref = ref_bytes[start - 1..start - 1 + 100].to_vec();
            bound_reads.push(mk_read(
                &format!("b_ref_{i}"),
                start,
                Flags::empty(),
                &[(100, Kind::Match)],
                seq_ref,
            ));
            let alt_start = 100_006 + i;
            let mut seq_alt = ref_bytes[alt_start - 1..alt_start - 1 + 100].to_vec();
            seq_alt[100_020 - alt_start] = b'A';
            bound_reads.push(mk_read(
                &format!("b_alt_{i}"),
                alt_start,
                Flags::empty(),
                &[(100, Kind::Match)],
                seq_alt,
            ));
        }
        bound_reads.sort_by_key(|r| r.alignment_start().unwrap().get());
        let bound_bam = dir.path().join("bound.bam");
        write_test_bam(&bound_bam, "chr1", ref_len, &bound_reads);

        let bound_untiled = run_simple_var(&bound_bam, Some(250_000_000), true);
        let bound_tiled = run_simple_var(&bound_bam, None, true);
        assert!(
            bound_untiled.contains("\t100020\t100020\tT\tA\t10\t5\t5\t0\t5\t0\tT/A\t0.5000\t"),
            "untiled bound with -t must have DP=10 VD=5 AF=0.5000, got:\n{bound_untiled}"
        );
        assert_eq!(
            bound_tiled, bound_untiled,
            "tiled bound with -t must drop Read B at 99,950 as duplicate of Read A and match untiled"
        );

        // 3. BAM-order MaxRlStep trajectory when Read A (99,900: 150M, overlaps Tile 2 at 100,001)
        //    precedes Read B (99,901: 90M 310S, ends at 99,991 <= 100,001):
        //    Read A must NOT see Read B's 400bp max_read_length when parsed on Tile 2, whereas a
        //    subsequent Read C (100,010: 100M) on Tile 2 MUST see 400bp.
        let order_reads = vec![
            mk_read(
                "read_A_overlap",
                99_900,
                Flags::empty(),
                &[(150, Kind::Match)],
                ref_bytes[99_899..100_049].to_vec(),
            ),
            mk_read(
                "read_B_non_overlap_400",
                99_901,
                Flags::empty(),
                &[(90, Kind::Match), (310, Kind::SoftClip)],
                ref_bytes[99_900..100_300].to_vec(),
            ),
            mk_read(
                "read_C_after",
                100_010,
                Flags::empty(),
                &[(100, Kind::Match)],
                ref_bytes[100_009..100_109].to_vec(),
            ),
        ];
        let order_bam = dir.path().join("order.bam");
        write_test_bam(&order_bam, "chr1", ref_len, &order_reads);
        let parent_reg = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(1_200_000),
            gene: "G".to_string(),
            amplicon: None,
        };
        let parent_tiles = crate::regions::tile(&parent_reg, 100_000, 500);
        let mut chr_lens = HashMap::new();
        chr_lens.insert("chr1".to_string(), ref_len);
        let order_cfg = crate::config::Config {
            reference: fasta_path.to_string_lossy().into_owned(),
            bam: order_bam.to_string_lossy().into_owned(),
            sample: "S".to_string(),
            ..Default::default()
        };
        let (steps_per_tile, final_max) = scan_multi_tile_max_read_lengths(
            &[order_bam.to_string_lossy().into_owned()],
            &parent_reg,
            &parent_tiles,
            &order_cfg,
            &chr_lens,
            0,
        )
        .unwrap();
        assert_eq!(final_max, 400);
        let tile2_steps = &steps_per_tile[1];
        assert_eq!(tile2_steps.len(), 2);
        // Step 0 is from Read A at 99,900 (same_start_overlap_before = 1, so it applies AFTER Read A)
        assert_eq!(tile2_steps[0].read_start, 99_900);
        assert_eq!(tile2_steps[0].same_start_overlap_before, 1);
        assert_eq!(tile2_steps[0].new_max, 150);
        // Step 1 is from Read B at 99,901 (applies before Read C at 100,010, NOT before Read A at 99,900)
        assert_eq!(tile2_steps[1].read_start, 99_901);
        assert_eq!(tile2_steps[1].same_start_overlap_before, 0);
        assert_eq!(tile2_steps[1].new_max, 400);
    }

    #[test]
    fn test_spillable_vcf_accumulator_k_way_merge_and_tie_break() {
        let mut acc: SpillableVcfAccumulator<RawVcfCandidate> =
            SpillableVcfAccumulator::with_limits(80, 2);

        // Insert chr2 first, then chr1 across multiple spill boundaries (with max_fanin=2
        // so chr1's >=3 spill runs exercise multi-pass intermediate merging).
        acc.insert(
            "chr2",
            RawVcfCandidate {
                start_pos: 500,
                af_4dp: 0.2500,
                should_emit: true,
                vcf_line: "chr2\t500\t.\tA\tG\tfirst_25".to_string(),
            },
        )
        .unwrap();
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 100,
                af_4dp: 0.5000,
                should_emit: true,
                vcf_line: "chr1\t100\t.\tC\tT\tfirst_50".to_string(),
            },
        )
        .unwrap();
        // Equal AF at chr1:100 in a later spill run must NOT overwrite first_50
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 100,
                af_4dp: 0.5000,
                should_emit: true,
                vcf_line: "chr1\t100\t.\tC\tA\tsecond_tie_50".to_string(),
            },
        )
        .unwrap();
        // Third spill run on chr1 with equal AF at 100 so > max_fanin(2) runs merge in passes
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 100,
                af_4dp: 0.5000,
                should_emit: true,
                vcf_line: "chr1\t100\t.\tC\tG\tthird_tie_50".to_string(),
            },
        )
        .unwrap();
        // Higher AF at chr2:500 in a later spill run MUST overwrite first_25
        acc.insert(
            "chr2",
            RawVcfCandidate {
                start_pos: 500,
                af_4dp: 0.7500,
                should_emit: true,
                vcf_line: "chr2\t500\t.\tA\tT\tlater_75".to_string(),
            },
        )
        .unwrap();
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 50,
                af_4dp: 0.3000,
                should_emit: true,
                vcf_line: "chr1\t50\t.\tG\tC\tearly_pos".to_string(),
            },
        )
        .unwrap();

        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt as _;
            let dir = acc.temp_dir.as_ref().expect("spill temp_dir must exist");
            let mode = std::fs::metadata(dir.path()).unwrap().permissions().mode() & 0o777;
            assert_eq!(mode, 0o700, "spill temp_dir must be owner-only (0o700)");
        }

        assert!(
            acc.spilled_runs.get("chr1").map_or(0, Vec::len) >= 3,
            "expected >= 3 disk spill runs on chr1 to exercise max_fanin=2 intermediate merge"
        );

        let mut emitted = Vec::new();
        acc.emit_all(|_chr, cand| {
            if cand.should_emit {
                emitted.push(cand.vcf_line);
            }
            Ok(())
        })
        .unwrap();

        assert_eq!(
            emitted,
            vec![
                "chr1\t50\t.\tG\tC\tearly_pos".to_string(),
                "chr1\t100\t.\tC\tT\tfirst_50".to_string(),
                "chr2\t500\t.\tA\tT\tlater_75".to_string(),
            ]
        );
    }

    #[test]
    fn resolve_mode_truth_table() {
        // 16-row truth table over (splice, has_region, amplicon_effective, paired)
        for splice in [false, true] {
            for has_region in [false, true] {
                for amp in [false, true] {
                    for paired in [false, true] {
                        let expected = if splice {
                            CallingMode::Splicing
                        } else if has_region || !amp {
                            if paired {
                                CallingMode::Somatic
                            } else {
                                CallingMode::Simple
                            }
                        } else {
                            CallingMode::Amplicon
                        };
                        assert_eq!(
                            resolve_mode(splice, has_region, amp, paired),
                            expected,
                            "mismatch for (splice={splice}, has_region={has_region}, amp={amp}, paired={paired})"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn process_amplicon_region_seed_chain() {
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::{Cigar, Data};
        use noodles::sam::alignment::RecordBuf;

        let dir = tempfile::TempDir::new().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        let ref_bases = "A".repeat(400);
        std::fs::write(&fasta_path, format!(">chr1\n{ref_bases}\n")).unwrap();
        std::fs::write(&fai_path, "chr1\t400\t6\t400\t401\n").unwrap();

        let bam_path = dir.path().join("amp.bam");
        let spliced_cigar: Cigar = [
            Op::new(Kind::Match, 10),
            Op::new(Kind::Skip, 20),
            Op::new(Kind::Match, 10),
        ]
        .into_iter()
        .collect();
        let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(0i32))]
            .into_iter()
            .collect();
        let rec = RecordBuf::builder()
            .set_name("r_amp_splice")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(10usize).unwrap())
            .set_mapping_quality(noodles::sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(spliced_cigar)
            .set_sequence(vec![b'A'; 20].into())
            .set_quality_scores(vec![40u8; 20].into())
            .set_data(data)
            .build();
        write_test_bam(&bam_path, "chr1", 400, &[rec]);

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 400usize);
        let cfg = crate::config::Config {
            bam: bam_path.to_string_lossy().into_owned(),
            reference: fasta_path.to_string_lossy().into_owned(),
            amplicon: Some("100:0.1".to_string()),
            emit_intermediate: true,
            ..Default::default()
        };
        let reg1 = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(10),
            end: crate::types::Pos1(49),
            gene: "G".to_string(),
            amplicon: Some(crate::regions::AmpliconInfo {
                insert_start: crate::types::Pos1(12),
                insert_end: crate::types::Pos1(47),
            }),
        };
        let reg2 = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(150),
            end: crate::types::Pos1(250),
            gene: "G".to_string(),
            amplicon: Some(crate::regions::AmpliconInfo {
                insert_start: crate::types::Pos1(160),
                insert_end: crate::types::Pos1(240),
            }),
        };

        let (_, splice1) = process_amplicon_region(
            &cfg.bam,
            &reg1,
            &cfg,
            &chr_lengths,
            std::collections::HashSet::new(),
        )
        .unwrap();
        assert!(splice1.contains("19-39"));

        let (_, splice2) =
            process_amplicon_region(&cfg.bam, &reg2, &cfg, &chr_lengths, splice1).unwrap();
        assert!(
            splice2.contains("19-39"),
            "second amplicon segment in cluster must preserve splice junctions from first segment"
        );
    }

    #[test]
    fn resolve_effective_config_splicing_amplicon_and_zero_based() {
        let dir = tempfile::TempDir::new().unwrap();
        let amp8_bed = dir.path().join("amp8.bed");
        std::fs::write(&amp8_bed, "chr1\t100\t200\tG1\t0\t+\t120\t180\n").unwrap();
        let amp9_bed = dir.path().join("amp9.bed");
        std::fs::write(&amp9_bed, "chr1\t100\t200\tG1\t0\t+\t120\t180\textra\n").unwrap();

        // 1. Splicing mode with -R and paired BAM resolves sample to first sample even without --emit-intermediate
        let splice_cfg = crate::config::Config {
            bam: "tumor.bam|normal.bam".to_string(),
            sample: "TUMOR|NORMAL".to_string(),
            sample_explicit: true,
            region: Some("chr1:100-200".to_string()),
            splice: true,
            emit_intermediate: false,
            ..Default::default()
        };
        let (eff_splice, mode_splice) = resolve_effective_config(&splice_cfg).unwrap();
        assert_eq!(mode_splice, CallingMode::Splicing);
        assert_eq!(eff_splice.sample, "TUMOR");
        assert_eq!(eff_splice.bam, "tumor.bam");

        // 2. Explicit -a on a 9-column BED preserves amplicon parameters
        let amp9_cfg = crate::config::Config {
            bam: "a.bam".to_string(),
            bed: Some(amp9_bed.to_string_lossy().into_owned()),
            amplicon: Some("15:0.90".to_string()),
            ..Default::default()
        };
        let (eff_amp9, mode_amp9) = resolve_effective_config(&amp9_cfg).unwrap();
        assert_eq!(mode_amp9, CallingMode::Amplicon);
        assert_eq!(eff_amp9.amplicon.as_deref(), Some("15:0.90"));

        // 3. 8-column amplicon BED auto-detects even when -c 1 (column_for_chromosome) is set,
        //    and preserves explicit -z 0 (`zero_based_defined = true, zero_based = false`)
        let amp8_z0_cfg = crate::config::Config {
            bam: "a.bam".to_string(),
            bed: Some(amp8_bed.to_string_lossy().into_owned()),
            column_for_chromosome: Some(0),
            zero_based: false,
            zero_based_defined: true,
            ..Default::default()
        };
        let (eff_amp8_z0, mode_amp8_z0) = resolve_effective_config(&amp8_z0_cfg).unwrap();
        assert_eq!(mode_amp8_z0, CallingMode::Amplicon);
        assert_eq!(eff_amp8_z0.amplicon.as_deref(), Some("10:0.95"));
        assert!(
            !eff_amp8_z0.zero_based,
            "explicit -z 0 must not be overwritten by amplicon auto-detect"
        );
    }

    #[test]
    fn test_split_var_line_fingerprint_and_split_tile_dedup() {
        let l1 = "S\tG\tchr1\t1000\t1000\tA\tG\t50\t25";
        let l2 = "S\tG\tchr1\t2000\t2000\tC\tT\t50\t25";
        assert_eq!(
            split_var_line_fingerprint(l1),
            split_var_line_fingerprint(l1)
        );
        assert_ne!(
            split_var_line_fingerprint(l1),
            split_var_line_fingerprint(l2)
        );

        // Simulate the consumer loop across:
        // 1) Split BED region 1 (tile 1: is_split_tile=true, is_first_split_tile=true;
        //    tile 2: is_split_tile=true, is_first_split_tile=false) carrying duplicate `l1`.
        // 2) Separate split BED region 2 (tile 1: is_split_tile=true, is_first_split_tile=true)
        //    carrying an identical line `l1` (must NOT be suppressed across BED regions).
        // 3) Untiled BED region 3 (is_split_tile=false, is_first_split_tile=true) carrying `l1` twice.
        let batches = vec![
            (vec![l1.to_string(), l2.to_string()], true, true),
            (vec![l1.to_string()], true, false),
            (vec![l1.to_string()], true, true),
            (vec![l1.to_string(), l1.to_string()], false, true),
        ];

        let mut seen_split_var: rustc_hash::FxHashSet<(u64, u64)> =
            rustc_hash::FxHashSet::default();
        let mut emitted = Vec::new();
        for (lines, is_split_tile, is_first_split_tile) in batches {
            if !is_split_tile || is_first_split_tile {
                seen_split_var.clear();
            }
            for line in lines {
                if is_split_tile {
                    let fp = split_var_line_fingerprint(&line);
                    if seen_split_var.contains(&fp) {
                        continue;
                    }
                    if seen_split_var.len() < MAX_SPLIT_VAR_DEDUP_ENTRIES {
                        seen_split_var.insert(fp);
                    }
                }
                emitted.push(line);
            }
        }

        assert_eq!(
            emitted,
            vec![
                l1.to_string(), // region 1 tile 1
                l2.to_string(), // region 1 tile 1
                // region 1 tile 2 duplicate l1 suppressed
                l1.to_string(), // region 2 tile 1 (new split region resets seen_split_var)
                l1.to_string(), // region 3 untiled (first copy)
                l1.to_string(), // region 3 untiled (second copy preserved when !is_split_tile)
            ]
        );
    }
}
