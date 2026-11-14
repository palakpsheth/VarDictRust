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

/// Sliding boundary-retention window (`200 kbp` = 2 default `100 kbp` tiles) for
/// cross-tile `.var` intermediate deduplication (`AGENTS.md §10`). Non-SV entries
/// with `anchor_pos < tile.start - SPLIT_VAR_RETENTION_WINDOW_BP` are safely evicted
/// as tiles advance monotonically, keeping deduplication continuously active across
/// whole-chromosome split runs without an arbitrary row-count cutoff.
const SPLIT_VAR_RETENTION_WINDOW_BP: i32 = 200_000;

/// Stable cross-tile variant identity for `.var` deduplication (independent of
/// tile-local depth, AF, status, SV counts, or debug measurements).
#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
struct SplitVarIdentity {
    origin_pos: i32,
    start_pos: i32,
    end_pos: i32,
    ref_allele: String,
    var_allele: String,
    null_record: bool,
    /// `0` = primary row, `1` = `-p` shifted-singleton null companion (`postprocess/simple.rs:103`),
    /// `2` = `-p` shifted-singleton pos-0 companion (`postprocess/simple.rs:104-107`).
    companion_kind: u8,
}

#[derive(Debug, Clone)]
struct SplitVarRow {
    id: SplitVarIdentity,
    is_sv: bool,
    line: String,
}

impl SplitVarRow {
    #[inline]
    fn anchor_pos(&self) -> i32 {
        self.id.origin_pos
    }
}

/// Builds `SplitVarIdentity` entries for a `simple_post_process` variant list,
/// tagging `-p` shifted-singleton companion triples (`[null_record_at_pos, pos0_default, vref]`)
/// with their originating pileup `position` and shifted allele so distinct shifted
/// singletons preserve their position-0 companion rows while cross-tile copies of
/// the same originating record deduplicate together.
fn build_simple_split_var_identities(
    variants: &[crate::variations::Variant],
    is_split_tile: bool,
) -> Vec<SplitVarIdentity> {
    if !is_split_tile {
        return vec![SplitVarIdentity::default(); variants.len()];
    }
    let mut out = Vec::with_capacity(variants.len());
    let mut idx = 0;
    while idx < variants.len() {
        if let (Some(v0), Some(v1), Some(v2)) = (
            variants.get(idx),
            variants.get(idx + 1),
            variants.get(idx + 2),
        ) {
            if v0.null_record && v1.start_pos == 0 && v1.end_pos == 0 && !v1.null_record {
                let origin_pos = v0.start_pos;
                let s_pos = v2.start_pos;
                let e_pos = v2.end_pos;
                let r_all = &v2.ref_allele;
                let v_all = &v2.var_allele;
                for kind in [1u8, 2u8, 0u8] {
                    out.push(SplitVarIdentity {
                        origin_pos,
                        start_pos: s_pos,
                        end_pos: e_pos,
                        ref_allele: r_all.clone(),
                        var_allele: v_all.clone(),
                        null_record: kind == 1,
                        companion_kind: kind,
                    });
                }
                idx += 3;
                continue;
            }
        }
        if let Some(v) = variants.get(idx) {
            out.push(SplitVarIdentity {
                origin_pos: v.start_pos,
                start_pos: v.start_pos,
                end_pos: v.end_pos,
                ref_allele: v.ref_allele.clone(),
                var_allele: v.var_allele.clone(),
                null_record: v.null_record,
                companion_kind: 0,
            });
        }
        idx += 1;
    }
    out
}

/// Cross-tile `.var` intermediate deduplicator using stable variant identity,
/// tile ownership (`tile.start..=tile.end`), within-tile multiplicity preservation,
/// and a sliding boundary-retention window (`SPLIT_VAR_RETENTION_WINDOW_BP`).
#[derive(Default)]
struct SplitVarDeduper {
    emitted_by_pos: std::collections::BTreeMap<i32, Vec<SplitVarIdentity>>,
    emitted_sv: rustc_hash::FxHashSet<SplitVarIdentity>,
    pending_forward: Vec<SplitVarRow>,
}

impl SplitVarDeduper {
    fn is_emitted(&self, row: &SplitVarRow) -> bool {
        if row.is_sv && self.emitted_sv.contains(&row.id) {
            return true;
        }
        self.emitted_by_pos
            .get(&row.anchor_pos())
            .is_some_and(|ids| ids.contains(&row.id))
    }

    fn record_emitted(&mut self, anchor_pos: i32, is_sv: bool, id: SplitVarIdentity) {
        if is_sv {
            self.emitted_sv.insert(id.clone());
        }
        let entry = self.emitted_by_pos.entry(anchor_pos).or_default();
        if !entry.contains(&id) {
            entry.push(id);
        }
    }

    fn flush_all_pending<F>(&mut self, mut emit: F) -> Result<(), VardictError>
    where
        F: FnMut(&str) -> Result<(), VardictError>,
    {
        let pending = std::mem::take(&mut self.pending_forward);
        for p in pending {
            if !self.is_emitted(&p) {
                emit(&p.line)?;
                self.record_emitted(p.anchor_pos(), p.is_sv, p.id);
            }
        }
        Ok(())
    }

    fn flush_before_unsplit_tile<F>(&mut self, emit: F) -> Result<(), VardictError>
    where
        F: FnMut(&str) -> Result<(), VardictError>,
    {
        self.flush_all_pending(emit)?;
        self.emitted_by_pos.clear();
        self.emitted_sv.clear();
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn accept_tile<F>(
        &mut self,
        rows: Vec<SplitVarRow>,
        is_split_tile: bool,
        tile_start: crate::types::Pos1,
        tile_end: crate::types::Pos1,
        original_start: crate::types::Pos1,
        original_end: crate::types::Pos1,
        mut emit: F,
    ) -> Result<(), VardictError>
    where
        F: FnMut(&str) -> Result<(), VardictError>,
    {
        if !is_split_tile {
            self.flush_before_unsplit_tile(&mut emit)?;
            for row in rows {
                emit(&row.line)?;
            }
            return Ok(());
        }

        if tile_start == original_start {
            self.flush_all_pending(&mut emit)?;
            self.emitted_by_pos.clear();
            self.emitted_sv.clear();
        } else {
            let cutoff = (tile_start.0 as i32).saturating_sub(SPLIT_VAR_RETENTION_WINDOW_BP);
            self.emitted_by_pos = self.emitted_by_pos.split_off(&cutoff);
        }

        let is_last_split_tile = tile_end == original_end;

        if !self.pending_forward.is_empty() {
            let old_pending = std::mem::take(&mut self.pending_forward);
            for p in old_pending {
                if is_last_split_tile || p.anchor_pos() <= tile_end.0 as i32 {
                    let owned_by_current = rows.iter().any(|r| {
                        r.id == p.id && (r.anchor_pos() <= tile_end.0 as i32 || is_last_split_tile)
                    });
                    if !owned_by_current && !self.is_emitted(&p) {
                        emit(&p.line)?;
                        self.record_emitted(p.anchor_pos(), p.is_sv, p.id);
                    }
                } else {
                    self.pending_forward.push(p);
                }
            }
        }

        let mut newly_emitted: Vec<(i32, bool, SplitVarIdentity)> = Vec::new();
        let mut newly_pending: Vec<SplitVarRow> = Vec::new();
        for row in rows {
            if self.is_emitted(&row) {
                continue;
            }
            if !is_last_split_tile && row.anchor_pos() > tile_end.0 as i32 {
                self.pending_forward.retain(|p| p.id != row.id);
                newly_pending.push(row);
            } else {
                emit(&row.line)?;
                newly_emitted.push((row.anchor_pos(), row.is_sv, row.id));
            }
        }
        self.pending_forward.extend(newly_pending);
        for (anchor_pos, is_sv, id) in newly_emitted {
            self.record_emitted(anchor_pos, is_sv, id);
        }

        if is_last_split_tile {
            self.flush_all_pending(&mut emit)?;
        }

        Ok(())
    }
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
    initial_splice: Option<rustc_hash::FxHashSet<String>>,
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
            read_end = read_end.saturating_add(i32::try_from(len).unwrap_or(i32::MAX));
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
    initial_splice: Option<rustc_hash::FxHashSet<String>>,
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
        window_bases: std::sync::Arc<[u8]>,
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
                window_bases: std::sync::Arc::from(&b""[..]),
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
        ) -> Result<(std::sync::Arc<[u8]>, i32), VardictError> {
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
                return Ok((std::sync::Arc::from(&b""[..]), 1));
            };
            if iv.stored_end < iv.sequence_start {
                return Ok((std::sync::Arc::from(&b""[..]), iv.sequence_start as i32));
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
                self.window_bases = std::sync::Arc::from(bases);
            }

            Ok((
                std::sync::Arc::clone(&self.window_bases),
                self.window_start as i32,
            ))
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
            let band_lo = read_start.saturating_sub(cur_max.saturating_mul(2)).max(1) as u32;
            let band_hi = span_end.saturating_add(cur_max.saturating_mul(2)).max(1) as u32;
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
                let total_mismatches = nm_val
                    .saturating_sub(i32::try_from(insertion_deletion_length).unwrap_or(i32::MAX));
                if total_mismatches > i32::try_from(config.max_mismatches).unwrap_or(i32::MAX) {
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
                let (ref_slice_arc, ref_start) =
                    scanner.slice_for_read(read_start, span_end, ext_start_i32, ext_end_i32)?;
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
                    std::mem::take(&mut r.sequence),
                    std::mem::take(&mut r.quality_scores),
                    &ref_slice_arc,
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
                        r.sequence = modified.query_seq;
                        r.quality_scores = modified.query_qual;
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

fn write_spill_u64<W: std::io::Write>(w: &mut W, v: u64) -> std::io::Result<()> {
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

fn read_spill_u64<R: std::io::Read>(r: &mut R) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
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

fn read_spill_bool<R: std::io::Read>(r: &mut R) -> std::io::Result<bool> {
    match read_spill_u8(r)? {
        0 => Ok(false),
        1 => Ok(true),
        b => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("corrupt spill boolean byte {b:#x}"),
        )),
    }
}

const MAX_SPILL_STR_BYTES: usize = 256 * 1024 * 1024;
const SPILL_STR_STACK_FAST_BYTES: usize = 4096;

fn read_spill_str<R: std::io::Read>(r: &mut R) -> std::io::Result<String> {
    use std::io::Read as _;
    let mut len_buf = [0u8; 8];
    r.read_exact(&mut len_buf)?;
    let raw_len = u64::from_le_bytes(len_buf);
    let len = usize::try_from(raw_len)
        .ok()
        .filter(|&n| n <= MAX_SPILL_STR_BYTES)
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("spill string length {raw_len} exceeds bound"),
            )
        })?;
    if len == 0 {
        return Ok(String::new());
    }
    if len <= SPILL_STR_STACK_FAST_BYTES {
        let mut stack_buf = [0u8; SPILL_STR_STACK_FAST_BYTES];
        let dst = stack_buf.get_mut(..len).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "spill string stack slice out of bounds",
            )
        })?;
        r.read_exact(dst)?;
        let s = std::str::from_utf8(dst)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        return Ok(s.to_owned());
    }
    let mut bytes = Vec::with_capacity(SPILL_STR_STACK_FAST_BYTES);
    r.take(raw_len).read_to_end(&mut bytes)?;
    if bytes.len() != len {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "truncated spill string",
        ));
    }
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
    event: crate::output::vcf::SingleVcfEvent,
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
        std::mem::size_of::<Self>()
            + match &self.event {
                crate::output::vcf::SingleVcfEvent::Record(rec) => rec.capacity(),
                _ => 0,
            }
    }

    fn write_to<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        write_spill_i32(w, self.start_pos)?;
        write_spill_f64(w, self.af_4dp)?;
        match &self.event {
            crate::output::vcf::SingleVcfEvent::Skip => write_spill_u8(w, 0),
            crate::output::vcf::SingleVcfEvent::AmpBiasSnvClusterProbe { start } => {
                write_spill_u8(w, 1)?;
                write_spill_i32(w, *start)
            }
            crate::output::vcf::SingleVcfEvent::SvSplitreadFlush {
                cluster_probe_start,
            } => {
                write_spill_u8(w, 2)?;
                match cluster_probe_start {
                    None => write_spill_u8(w, 0),
                    Some(start) => {
                        write_spill_u8(w, 1)?;
                        write_spill_i32(w, *start)
                    }
                }
            }
            crate::output::vcf::SingleVcfEvent::Record(rec) => {
                let (filter_start, filter_end) = rec.filter_offsets();
                write_spill_u8(w, 3)?;
                write_spill_str(w, rec.as_str())?;
                write_spill_u64(w, filter_start as u64)?;
                write_spill_u64(w, filter_end as u64)?;
                write_spill_i32(w, rec.start)?;
                write_spill_u8(w, u8::from(rec.type_snv))?;
                write_spill_u8(w, u8::from(rec.pre_cluster_non_amp_empty))?;
                write_spill_u8(w, u8::from(rec.has_ampbias))
            }
        }
    }

    fn read_from<R: std::io::Read>(r: &mut R) -> std::io::Result<Option<Self>> {
        let Some(start_pos) = read_spill_i32_opt(r)? else {
            return Ok(None);
        };
        let af_4dp = read_spill_f64(r)?;
        let event = match read_spill_u8(r)? {
            0 => crate::output::vcf::SingleVcfEvent::Skip,
            1 => {
                let start = read_spill_i32(r)?;
                crate::output::vcf::SingleVcfEvent::AmpBiasSnvClusterProbe { start }
            }
            2 => {
                let cluster_probe_start = if read_spill_bool(r)? {
                    Some(read_spill_i32(r)?)
                } else {
                    None
                };
                crate::output::vcf::SingleVcfEvent::SvSplitreadFlush {
                    cluster_probe_start,
                }
            }
            3 => {
                let line = read_spill_str(r)?;
                let filter_start = usize::try_from(read_spill_u64(r)?)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
                let filter_end = usize::try_from(read_spill_u64(r)?)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
                let start = read_spill_i32(r)?;
                let type_snv = read_spill_bool(r)?;
                let pre_cluster_non_amp_empty = read_spill_bool(r)?;
                let has_ampbias = read_spill_bool(r)?;
                let rec = crate::output::vcf::SingleVcfRecord::from_parts_checked(
                    line,
                    filter_start,
                    filter_end,
                    start,
                    type_snv,
                    pre_cluster_non_amp_empty,
                    has_ampbias,
                )
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "corrupt SingleVcfRecord filter span in spill file",
                    )
                })?;
                crate::output::vcf::SingleVcfEvent::Record(rec)
            }
            tag => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("unknown SingleVcfEvent tag {tag:#x} in spill file"),
                ));
            }
        };
        Ok(Some(Self {
            start_pos,
            af_4dp,
            event,
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
        std::mem::size_of::<Self>() + self.record.as_ref().map_or(0, |rec| rec.capacity())
    }

    fn write_to<W: std::io::Write>(&self, w: &mut W) -> std::io::Result<()> {
        write_spill_i32(w, self.start_pos)?;
        write_spill_f64(w, self.af_4dp)?;
        match &self.record {
            None => write_spill_u8(w, 0),
            Some(rec) => {
                let (filter_start, filter_end) = rec.filter_offsets();
                write_spill_u8(w, 1)?;
                write_spill_str(w, rec.as_str())?;
                write_spill_u64(w, filter_start as u64)?;
                write_spill_u64(w, filter_end as u64)?;
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
        let has_rec = read_spill_bool(r)?;
        let record = if has_rec {
            let line = read_spill_str(r)?;
            let filter_start = usize::try_from(read_spill_u64(r)?)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let filter_end = usize::try_from(read_spill_u64(r)?)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let start = read_spill_i32(r)?;
            let type_snv = read_spill_bool(r)?;
            let status_somatic = read_spill_bool(r)?;
            let rec = crate::output::vcf::PairedVcfRecord::from_parts_checked(
                line,
                filter_start,
                filter_end,
                start,
                type_snv,
                status_somatic,
            )
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "corrupt PairedVcfRecord filter span in spill file",
                )
            })?;
            Some(rec)
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

#[cfg(test)]
thread_local! {
    static TEST_VCF_SPILL_BYTES_OVERRIDE: std::cell::Cell<Option<usize>> =
        const { std::cell::Cell::new(None) };
}

impl<V: SpillableCandidate> SpillableVcfAccumulator<V> {
    fn new() -> Self {
        #[cfg(test)]
        if let Some(override_bytes) = TEST_VCF_SPILL_BYTES_OVERRIDE.with(std::cell::Cell::get) {
            return Self::with_max_bytes(override_bytes);
        }
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
        let chr_map = if let Some(m) = self.in_memory.get_mut(chr) {
            m
        } else {
            if !self.spilled_runs.contains_key(chr) {
                self.spilled_runs.insert(chr.to_string(), Vec::new());
            }
            self.in_memory.entry(chr.to_string()).or_default()
        };
        match chr_map.entry(cand.start_pos()) {
            std::collections::btree_map::Entry::Vacant(e) => {
                let new_bytes = cand.approx_heap_bytes().saturating_add(32);
                self.current_bytes = self.current_bytes.saturating_add(new_bytes);
                e.insert(cand);
            }
            std::collections::btree_map::Entry::Occupied(mut e) => {
                if cand.af_4dp() > e.get().af_4dp() {
                    let new_bytes = cand.approx_heap_bytes().saturating_add(32);
                    let old_bytes = e.get().approx_heap_bytes().saturating_add(32);
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
            let mut writer = std::io::BufWriter::with_capacity(BUFFER_CAPACITY, file);
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
                        chunk_sources.push(RunSource::File(std::io::BufReader::with_capacity(
                            BUFFER_CAPACITY,
                            f,
                        )));
                    }
                    let out_file =
                        create_private_spill_file(&out_path).map_err(VardictError::Io)?;
                    let mut writer = std::io::BufWriter::with_capacity(BUFFER_CAPACITY, out_file);
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
                sources.push(RunSource::File(std::io::BufReader::with_capacity(
                    BUFFER_CAPACITY,
                    f,
                )));
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

fn emit_single_vcf_accumulator<W: std::io::Write>(
    vcf_acc: SpillableVcfAccumulator<RawVcfCandidate>,
    config: &config::Config,
    writer: &mut W,
) -> Result<(), VardictError> {
    let mut cur_chr: Option<String> = None;
    let mut streamer = crate::output::vcf::SingleRecordStreamer::from_config(config);
    vcf_acc.emit_all(|chr, cand| {
        if cur_chr.as_deref() != Some(chr) {
            if cur_chr.is_some() {
                let prev = std::mem::replace(
                    &mut streamer,
                    crate::output::vcf::SingleRecordStreamer::from_config(config),
                );
                prev.finish(|line| writeln!(writer, "{line}"))
                    .map_err(VardictError::Io)?;
            }
            cur_chr = Some(chr.to_string());
        }
        streamer
            .push_event(cand.event, |line| writeln!(writer, "{line}"))
            .map_err(VardictError::Io)
    })?;
    if cur_chr.is_some() {
        streamer
            .finish(|line| writeln!(writer, "{line}"))
            .map_err(VardictError::Io)?;
    }
    Ok(())
}

/// Output destination wrapper around buffered stdout, stderr, plain file, or multithreaded BGZF.
/// Uses unlocked BufWriter to eliminate mutex contention with worker diagnostics
/// while reducing syscall overhead by ~1000x for deep WGS output.
pub enum OutputDestination {
    Stdout(std::io::BufWriter<std::io::Stdout>),
    Stderr(std::io::BufWriter<std::io::Stderr>),
    File(std::io::BufWriter<std::fs::File>),
    Bgzf(Box<noodles::bgzf::io::MultithreadedWriter<std::fs::File>>),
}

impl std::fmt::Debug for OutputDestination {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Stdout(_) => f.debug_tuple("Stdout").finish_non_exhaustive(),
            Self::Stderr(_) => f.debug_tuple("Stderr").finish_non_exhaustive(),
            Self::File(_) => f.debug_tuple("File").finish_non_exhaustive(),
            Self::Bgzf(_) => f.debug_tuple("Bgzf").finish_non_exhaustive(),
        }
    }
}

impl std::io::Write for OutputDestination {
    #[inline]
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Stdout(w) => w.write(buf),
            Self::Stderr(w) => w.write(buf),
            Self::File(w) => w.write(buf),
            Self::Bgzf(w) => w.write(buf),
        }
    }

    #[inline]
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Stdout(w) => w.flush(),
            Self::Stderr(w) => w.flush(),
            Self::File(w) => w.flush(),
            Self::Bgzf(w) => w.flush(),
        }
    }
}

const BUFFER_CAPACITY: usize = 64 * 1024; // 64 KB

fn sync_file_ignoring_special(file: &std::fs::File) -> Result<(), VardictError> {
    // POSIX/Darwin ENOTTY (25): returned by macOS F_FULLFSYNC / fsync on character devices
    // when accessed via an fd alias or race.
    #[cfg(unix)]
    const UNIX_ENOTTY: i32 = 25;

    if !file.metadata().is_ok_and(|m| m.is_file()) {
        return Ok(());
    }
    match file.sync_all() {
        Ok(()) => Ok(()),
        Err(err)
            if matches!(
                err.kind(),
                std::io::ErrorKind::InvalidInput
                    | std::io::ErrorKind::NotSeekable
                    | std::io::ErrorKind::Unsupported
            ) =>
        {
            Ok(())
        }
        #[cfg(unix)]
        Err(err) if err.raw_os_error() == Some(UNIX_ENOTTY) => Ok(()),
        Err(err) => Err(VardictError::Io(err)),
    }
}

static NEXT_TMP_FILE_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);

fn append_os_suffix(path: &Path, suffix: &str) -> std::path::PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(suffix);
    std::path::PathBuf::from(s)
}

fn atomic_tmp_path(target: &Path) -> std::path::PathBuf {
    let seq = NEXT_TMP_FILE_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let pid = std::process::id();
    append_os_suffix(target, &format!(".tmp.{pid}.{seq}"))
}

fn tbi_output_path(vcf_gz_path: &Path) -> std::path::PathBuf {
    append_os_suffix(vcf_gz_path, ".tbi")
}

fn csi_output_path(vcf_gz_path: &Path) -> std::path::PathBuf {
    append_os_suffix(vcf_gz_path, ".csi")
}

fn is_readonly_parent_fallback_eligible(target: &Path, err: &std::io::Error) -> bool {
    !is_special_stream_or_fd_alias(target)
        && matches!(
            err.kind(),
            std::io::ErrorKind::PermissionDenied | std::io::ErrorKind::ReadOnlyFilesystem
        )
        && std::fs::metadata(target).is_ok_and(|m| m.is_file())
}

fn preflight_readonly_index_plan(plan: &OutputIndexPlan) -> std::io::Result<()> {
    let idx_probe = std::fs::OpenOptions::new()
        .write(true)
        .open(&plan.primary_path)
        .map_err(|e| {
            std::io::Error::new(
                e.kind(),
                format!(
                    "cannot write {} index '{}' in read-only parent directory: {e}",
                    plan.flag,
                    plan.primary_path.display()
                ),
            )
        })?;
    if !idx_probe.metadata()?.is_file() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            format!(
                "cannot stage {} index '{}': target is not a writable regular file",
                plan.flag,
                plan.primary_path.display()
            ),
        ));
    }
    drop(idx_probe);

    if plan.alt_path.exists() || std::fs::symlink_metadata(&plan.alt_path).is_ok() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            format!(
                "cannot remove stale {} index '{}' in read-only parent directory",
                plan.alt_flag,
                plan.alt_path.display()
            ),
        ));
    }
    Ok(())
}

fn is_directory_target(path: &Path) -> bool {
    path.is_dir() || std::fs::symlink_metadata(path).is_ok_and(|m| m.is_dir())
}

fn directory_target_io_error(label: &str, path: &Path) -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::IsADirectory,
        format!("{label} '{}' is a directory", path.display()),
    )
}

enum StagedFileTarget {
    AtomicTmp {
        final_path: std::path::PathBuf,
        tmp_path: std::path::PathBuf,
        committed: bool,
    },
    DirectTempStaged {
        final_path: std::path::PathBuf,
        _scratch_dir: SpillTempDir,
        scratch_path: std::path::PathBuf,
    },
    Direct,
}

impl StagedFileTarget {
    fn open(target: &Path) -> std::io::Result<(std::fs::File, Self)> {
        Self::open_with_index_preflight(target, None)
    }

    fn open_with_index_preflight(
        target: &Path,
        index_plan: Option<&OutputIndexPlan>,
    ) -> std::io::Result<(std::fs::File, Self)> {
        if is_directory_target(target) {
            return Err(directory_target_io_error("output target", target));
        }
        if is_special_stream_or_fd_alias(target) {
            let file = open_special_or_fd_alias_for_write(target)?;
            return Ok((file, Self::Direct));
        }
        let tmp_path = atomic_tmp_path(target);
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&tmp_path)
        {
            Ok(file) => Ok((
                file,
                Self::AtomicTmp {
                    final_path: target.to_path_buf(),
                    tmp_path,
                    committed: false,
                },
            )),
            Err(err) if is_readonly_parent_fallback_eligible(target, &err) => {
                // Probe writability WITHOUT `.truncate(true)` so pre-existing target bytes
                // are never clobbered before index preflight and engine execution complete.
                let probe = std::fs::OpenOptions::new().write(true).open(target)?;
                if !probe.metadata()?.is_file() {
                    return Err(err);
                }
                drop(probe);

                if let Some(plan) = index_plan {
                    preflight_readonly_index_plan(plan)?;
                }

                match SpillTempDir::new() {
                    Ok(scratch_dir) => {
                        let scratch_path = scratch_dir.path().join("staged_output.bin");
                        let file = create_private_spill_file(&scratch_path)?;
                        Ok((
                            file,
                            Self::DirectTempStaged {
                                final_path: target.to_path_buf(),
                                _scratch_dir: scratch_dir,
                                scratch_path,
                            },
                        ))
                    }
                    Err(tmp_err)
                        if matches!(
                            tmp_err.kind(),
                            std::io::ErrorKind::PermissionDenied
                                | std::io::ErrorKind::ReadOnlyFilesystem
                        ) =>
                    {
                        let file = std::fs::OpenOptions::new()
                            .write(true)
                            .truncate(true)
                            .open(target)?;
                        Ok((file, Self::Direct))
                    }
                    Err(tmp_err) => Err(tmp_err),
                }
            }
            Err(err) => Err(err),
        }
    }

    fn readable_path<'a>(&'a self, final_path: &'a Path) -> &'a Path {
        match self {
            Self::AtomicTmp { tmp_path, .. } => tmp_path.as_path(),
            Self::DirectTempStaged { scratch_path, .. } => scratch_path.as_path(),
            Self::Direct => final_path,
        }
    }

    fn commit(mut self) -> std::io::Result<()> {
        match &mut self {
            Self::AtomicTmp {
                final_path,
                tmp_path,
                committed,
            } => {
                std::fs::rename(&*tmp_path, &*final_path)?;
                *committed = true;
            }
            Self::DirectTempStaged {
                final_path,
                scratch_path,
                ..
            } => {
                let mut src = std::fs::File::open(&*scratch_path)?;
                let mut dst = std::fs::OpenOptions::new()
                    .write(true)
                    .truncate(true)
                    .open(&*final_path)?;
                std::io::copy(&mut src, &mut dst)?;
                std::io::Write::flush(&mut dst)?;
                sync_file_ignoring_special(&dst).map_err(|e| match e {
                    VardictError::Io(io_e) => io_e,
                    other => std::io::Error::other(other.to_string()),
                })?;
            }
            Self::Direct => {}
        }
        Ok(())
    }
}

impl Drop for StagedFileTarget {
    fn drop(&mut self) {
        if let Self::AtomicTmp {
            tmp_path,
            committed,
            ..
        } = self
        {
            if !*committed {
                let _ = std::fs::remove_file(tmp_path);
            }
        }
    }
}

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

fn is_bgzf_output_path(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz") || ext.eq_ignore_ascii_case("bgz"))
}

fn open_output_writer(
    config: &config::Config,
) -> Result<(OutputDestination, Option<StagedFileTarget>), VardictError> {
    match config.output_path.as_deref() {
        None => Ok((create_output_writer(config), None)),
        Some("-") => Ok((
            OutputDestination::Stdout(std::io::BufWriter::with_capacity(
                BUFFER_CAPACITY,
                std::io::stdout(),
            )),
            None,
        )),
        Some(path_str) => {
            let target = Path::new(path_str);
            let index_plan = output_index_plan(config);
            let (file, staged) =
                StagedFileTarget::open_with_index_preflight(target, index_plan.as_ref())
                    .map_err(VardictError::Io)?;
            if is_bgzf_output_path(target) {
                let worker_count = std::num::NonZeroUsize::new(config.threads.clamp(1, 4))
                    .unwrap_or(std::num::NonZeroUsize::MIN);
                Ok((
                    OutputDestination::Bgzf(Box::new(
                        noodles::bgzf::io::MultithreadedWriter::with_worker_count(
                            worker_count,
                            file,
                        ),
                    )),
                    Some(staged),
                ))
            } else {
                Ok((
                    OutputDestination::File(std::io::BufWriter::with_capacity(
                        BUFFER_CAPACITY,
                        file,
                    )),
                    Some(staged),
                ))
            }
        }
    }
}

/// Internal writer abstraction that preserves `run_with_writer`'s complete-stream contract
/// on injected writers while routing `--DP ERR` column headers to `stdout` when invoked via
/// the CLI [`OutputDestination`].
trait EngineWriter: std::io::Write + Sized {
    fn emit_column_header(&mut self, config: &config::Config) -> Result<(), VardictError> {
        crate::output::table::print_header(self, config).map_err(VardictError::Io)
    }

    fn finish_output(mut self) -> Result<(), VardictError> {
        self.flush().map_err(VardictError::Io)
    }
}

impl EngineWriter for OutputDestination {
    fn emit_column_header(&mut self, config: &config::Config) -> Result<(), VardictError> {
        match self {
            Self::Stdout(_) | Self::File(_) | Self::Bgzf(_) => {
                crate::output::table::print_header(self, config).map_err(VardictError::Io)
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

    fn finish_output(self) -> Result<(), VardictError> {
        use std::io::Write as _;
        match self {
            Self::Stdout(mut w) => w.flush().map_err(VardictError::Io),
            Self::Stderr(mut w) => w.flush().map_err(VardictError::Io),
            Self::File(mut w) => {
                w.flush().map_err(VardictError::Io)?;
                sync_file_ignoring_special(w.get_ref())
            }
            Self::Bgzf(mut w) => {
                let file = w.finish().map_err(VardictError::Io)?;
                sync_file_ignoring_special(&file)
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

const TABIX_MAX_COORD_BP: usize = 536_870_912;

fn scan_bgzf_vcf_index_records<F>(vcf_gz_path: &Path, mut on_record: F) -> Result<(), VardictError>
where
    F: FnMut(
        &str,
        usize,
        usize,
        noodles::core::Position,
        noodles::core::Position,
        noodles::csi::binning_index::index::reference_sequence::bin::Chunk,
    ) -> Result<(), VardictError>,
{
    let mut reader =
        noodles::bgzf::io::Reader::new(std::fs::File::open(vcf_gz_path).map_err(VardictError::Io)?);
    let mut line = String::new();
    loop {
        let start_vp = reader.virtual_position();
        let bytes_read =
            std::io::BufRead::read_line(&mut reader, &mut line).map_err(VardictError::Io)?;
        if bytes_read == 0 {
            break;
        }
        let end_vp = reader.virtual_position();
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.is_empty() || trimmed.starts_with('#') {
            line.clear();
            continue;
        }
        let mut fields = trimmed.split('\t');
        let chrom = fields.next().ok_or_else(|| {
            VardictError::Parse("missing VCF CHROM column while indexing".to_string())
        })?;
        let pos_str = fields.next().ok_or_else(|| {
            VardictError::Parse("missing VCF POS column while indexing".to_string())
        })?;
        let _id = fields.next().ok_or_else(|| {
            VardictError::Parse("missing VCF ID column while indexing".to_string())
        })?;
        let ref_allele = fields.next().ok_or_else(|| {
            VardictError::Parse("missing VCF REF column while indexing".to_string())
        })?;
        let _alt = fields.next();
        let _qual = fields.next();
        let _filter = fields.next();
        let info_end = fields.next().and_then(|info| {
            info.split(';')
                .find_map(|kv| kv.strip_prefix("END=")?.parse::<usize>().ok())
        });
        let pos_usize = pos_str.parse::<usize>().map_err(|e| {
            VardictError::Parse(format!("invalid VCF POS '{pos_str}' while indexing: {e}"))
        })?;
        let start_pos = noodles::core::Position::new(pos_usize).ok_or_else(|| {
            VardictError::Parse(format!(
                "VCF POS must be >= 1 while indexing, got {pos_usize}"
            ))
        })?;
        let end_usize = info_end
            .unwrap_or_else(|| pos_usize.saturating_add(ref_allele.len().saturating_sub(1)))
            .max(pos_usize);
        let end_pos = noodles::core::Position::new(end_usize)
            .ok_or_else(|| VardictError::Parse(format!("invalid VCF END {end_usize}")))?;
        let chunk = noodles::csi::binning_index::index::reference_sequence::bin::Chunk::new(
            start_vp, end_vp,
        );
        on_record(chrom, pos_usize, end_usize, start_pos, end_pos, chunk)?;
        line.clear();
    }
    Ok(())
}

fn stage_index_file<F>(index_path: &Path, write_fn: F) -> Result<StagedFileTarget, VardictError>
where
    F: FnOnce(std::fs::File) -> std::io::Result<std::fs::File>,
{
    let (file, staged) = StagedFileTarget::open(index_path).map_err(VardictError::Io)?;
    let inner_file = write_fn(file).map_err(VardictError::Io)?;
    sync_file_ignoring_special(&inner_file)?;
    drop(inner_file);
    Ok(staged)
}

fn stage_vcf_tabix_index_inner(
    readable_vcf: &Path,
    tbi_path: &Path,
) -> Result<StagedFileTarget, VardictError> {
    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(noodles::csi::binning_index::index::header::Builder::vcf().build());

    scan_bgzf_vcf_index_records(
        readable_vcf,
        |chrom, pos_usize, end_usize, start_pos, end_pos, chunk| {
            if end_usize > TABIX_MAX_COORD_BP {
                return Err(VardictError::Parse(format!(
                    "Tabix (.tbi) index cannot represent coordinates > 536,870,912 bp (saw {chrom}:{pos_usize}-{end_usize}); re-run with --csi instead of --tbi"
                )));
            }
            indexer
                .add_record(chrom, start_pos, end_pos, chunk)
                .map_err(|e| {
                    VardictError::Io(std::io::Error::new(std::io::ErrorKind::InvalidData, e))
                })
        },
    )?;

    let index = indexer.build();
    stage_index_file(tbi_path, |file| {
        let mut tbi_writer = noodles::tabix::io::Writer::new(file);
        tbi_writer.write_index(&index)?;
        tbi_writer.into_inner().finish()
    })
}

fn stage_vcf_tabix_index(
    readable_vcf: &Path,
    final_vcf: &Path,
) -> Result<StagedFileTarget, VardictError> {
    let tbi_path = tbi_output_path(final_vcf);
    if paths_collide(final_vcf, &tbi_path) {
        return Err(VardictError::Parse(format!(
            "--output '{}' and its Tabix index '{}' must not point to the same file",
            final_vcf.display(),
            tbi_path.display()
        )));
    }
    stage_vcf_tabix_index_inner(readable_vcf, &tbi_path)
}

fn stage_vcf_csi_index_inner(
    readable_vcf: &Path,
    csi_path: &Path,
) -> Result<StagedFileTarget, VardictError> {
    let mut indexer = noodles::csi::binning_index::Indexer::<
        noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex,
    >::new(14, 6)
    .ok_or_else(|| VardictError::Parse("invalid CSI binning parameters (14, 6)".to_string()))?;
    let mut names = noodles::csi::binning_index::index::header::ReferenceSequenceNames::new();

    scan_bgzf_vcf_index_records(
        readable_vcf,
        |chrom, _pos_usize, _end_usize, start_pos, end_pos, chunk| {
            let (ref_id, _) = names.insert_full(chrom.as_bytes().into());
            indexer
                .add_record(Some((ref_id, start_pos, end_pos, true)), chunk)
                .map_err(|e| {
                    VardictError::Io(std::io::Error::new(std::io::ErrorKind::InvalidData, e))
                })
        },
    )?;

    let contig_count = names.len();
    let header = noodles::csi::binning_index::index::header::Builder::vcf()
        .set_reference_sequence_names(names)
        .build();
    let index = indexer.set_header(header).build(contig_count);
    stage_index_file(csi_path, |file| {
        let mut csi_writer = noodles::csi::io::Writer::new(file);
        csi_writer.write_index(&index)?;
        csi_writer.into_inner().finish()
    })
}

fn stage_vcf_csi_index(
    readable_vcf: &Path,
    final_vcf: &Path,
) -> Result<StagedFileTarget, VardictError> {
    let csi_path = csi_output_path(final_vcf);
    if paths_collide(final_vcf, &csi_path) {
        return Err(VardictError::Parse(format!(
            "--output '{}' and its CSI index '{}' must not point to the same file",
            final_vcf.display(),
            csi_path.display()
        )));
    }
    stage_vcf_csi_index_inner(readable_vcf, &csi_path)
}

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

impl CallingMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Splicing => "Splicing",
            Self::Amplicon => "Amplicon",
            Self::Somatic => "Somatic",
            Self::Simple => "Simple",
        }
    }
}

#[derive(Default)]
struct RunTelemetryCounters {
    regions_count: std::sync::atomic::AtomicUsize,
    tiles_count: std::sync::atomic::AtomicUsize,
}

fn read_peak_rss_kb() -> u64 {
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return 0;
    };
    let parse_kb = |prefix: &str| -> Option<u64> {
        status.lines().find_map(|line| {
            let rest = line.strip_prefix(prefix)?;
            rest.split_whitespace().next()?.parse::<u64>().ok()
        })
    };
    parse_kb("VmHWM:")
        .or_else(|| parse_kb("VmRSS:"))
        .unwrap_or(0)
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

/// Resolves effective [`config::Config`], [`CallingMode`], and optional preloaded BED lines
/// matching `VarDictLauncher.initResources` and `VarDictLauncher.start`
fn is_regular_seekable_bed(bed_path: &str) -> bool {
    let p = Path::new(bed_path);
    !targets_os_fd_alias(p) && std::fs::metadata(p).is_ok_and(|m| m.is_file())
}

/// Resolves the effective configuration (`amplicon`, `zero_based`, `zero_based_defined`, `sample`)
/// and calling mode without mutating the caller's `Config`. Streams regular BED files in `O(1)`
/// memory and buffers non-regular BED streams (FIFOs/pipes/`/dev/fd/N`) in a single pass (`AR-4`).
fn resolve_effective_config(
    config: &config::Config,
) -> Result<(config::Config, CallingMode, Option<Vec<String>>), VardictError> {
    let mut effective = config.clone();
    if let Some(ref raw) = effective.sample_regexp {
        crate::config::validate_sample_regexp(raw).map_err(VardictError::Parse)?;
    }
    effective.validate_vcf_thresholds()?;
    let paired = effective.bam.contains('|');
    let has_region = effective.region.is_some();
    let mut preloaded_bed_lines: Option<Vec<String>> = None;

    if has_region {
        effective.amplicon = None;
    } else if effective.amplicon.is_none() {
        if let Some(ref bed_path) = effective.bed {
            let is_amplicon = if is_regular_seekable_bed(bed_path) {
                crate::regions::detect_amplicon_bed(bed_path, &effective.delimiter)?
            } else {
                let lines = crate::regions::read_bed_lines(bed_path)?;
                let detected =
                    crate::regions::detect_amplicon_bed_from_lines(&lines, &effective.delimiter)?;
                preloaded_bed_lines = Some(lines);
                detected
            };
            if is_amplicon {
                effective.amplicon = Some("10:0.95".to_string());
                if !effective.zero_based_defined {
                    effective.zero_based = true;
                    effective.zero_based_defined = true;
                }
            }
        }
    }

    if effective.sample_regexp.is_some()
        || (has_region && paired && (effective.emit_intermediate || effective.splice))
    {
        effective.sample = crate::config::resolve_sample_names(
            &effective.bam,
            effective
                .sample_explicit
                .then_some(effective.sample.as_str()),
            effective.sample_regexp.as_deref(),
            has_region,
            effective.emit_intermediate || effective.splice,
        );
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

    validate_index_flags(&effective)?;

    Ok((effective, mode, preloaded_bed_lines))
}

fn validate_index_flags(config: &config::Config) -> Result<(), VardictError> {
    if config.index_output && config.index_csi {
        return Err(VardictError::Parse(
            "--tbi/--index and --csi cannot be used together; choose --tbi (<= 512 Mbp) or --csi (up to 2.14 Gbp)"
                .to_string(),
        ));
    }
    if config.index_output || config.index_csi {
        let flag_label = if config.index_output {
            "--tbi/--index"
        } else {
            "--csi"
        };
        let Some(out_path) = config.output_path.as_deref() else {
            return Err(VardictError::Parse(format!(
                "{flag_label} requires --output <PATH> ending in .vcf.gz or .bgz"
            )));
        };
        let out_path_ref = Path::new(out_path);
        if out_path == "-" || is_special_stream_or_fd_alias(out_path_ref) {
            return Err(VardictError::Parse(
                "--tbi/--csi requires a regular seekable --output file".to_string(),
            ));
        }
        if !is_bgzf_output_path(out_path_ref) {
            return Err(VardictError::Parse(format!(
                "{flag_label} requires --output <PATH> ending in .vcf.gz or .bgz"
            )));
        }
        if config.emit_intermediate || config.splice || config.column_header {
            return Err(VardictError::Parse(format!(
                "{flag_label} is only supported for VCF output (cannot be combined with --emit-intermediate, -i/--splice, or -h/--header)"
            )));
        }
    }
    Ok(())
}

fn resolve_target_path(path: &Path) -> Option<std::path::PathBuf> {
    if let Ok(canon) = std::fs::canonicalize(path) {
        return Some(canon);
    }
    let mut cur = path.to_path_buf();
    for _ in 0..40 {
        match std::fs::read_link(&cur) {
            Ok(target) => {
                cur = if target.is_absolute() {
                    target
                } else {
                    cur.parent()
                        .filter(|p| !p.as_os_str().is_empty())
                        .unwrap_or_else(|| Path::new("."))
                        .join(target)
                };
                if let Ok(canon) = std::fs::canonicalize(&cur) {
                    return Some(canon);
                }
            }
            Err(_) => break,
        }
    }
    normalize_parent_preserving_leaf(&cur)
}

fn normalize_parent_preserving_leaf(p: &Path) -> Option<std::path::PathBuf> {
    let file_name = p.file_name()?;
    let parent = p
        .parent()
        .filter(|par| !par.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let canon_parent = std::fs::canonicalize(parent).ok()?;
    Some(canon_parent.join(file_name))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OsFdAliasKind {
    SelfStdout,
    SelfStderr,
    SelfOtherFd(u32),
    OtherFd,
}

fn is_ascii_digits_nonempty(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

fn classify_self_fd_str(n: &str) -> OsFdAliasKind {
    match n {
        "1" => OsFdAliasKind::SelfStdout,
        "2" => OsFdAliasKind::SelfStderr,
        _ => match n.parse::<u32>() {
            Ok(fd) => OsFdAliasKind::SelfOtherFd(fd),
            Err(_) => OsFdAliasKind::OtherFd,
        },
    }
}

fn classify_os_fd_alias_path(path: &Path) -> Option<OsFdAliasKind> {
    use std::path::Component;

    let mut rooted = false;
    let mut parts: Vec<&str> = Vec::new();
    for comp in path.components() {
        match comp {
            Component::Prefix(_) => return None,
            Component::RootDir => {
                rooted = true;
                parts.clear();
            }
            Component::CurDir => {}
            Component::ParentDir => {
                if !rooted {
                    return None;
                }
                let _ = parts.pop();
            }
            Component::Normal(os) => {
                parts.push(os.to_str()?);
            }
        }
    }
    if !rooted {
        return None;
    }
    match parts.as_slice() {
        ["dev", "stdout"] => Some(OsFdAliasKind::SelfStdout),
        ["dev", "stderr"] => Some(OsFdAliasKind::SelfStderr),
        ["dev", "stdin"] => Some(OsFdAliasKind::OtherFd),
        ["dev", "fd", n] if is_ascii_digits_nonempty(n) => Some(classify_self_fd_str(n)),
        ["proc", pid_seg, "fd", n]
            if is_ascii_digits_nonempty(n)
                && (*pid_seg == "self"
                    || *pid_seg == "thread-self"
                    || is_ascii_digits_nonempty(pid_seg)) =>
        {
            let is_self = *pid_seg == "self"
                || *pid_seg == "thread-self"
                || pid_seg.parse::<u32>().ok() == Some(std::process::id());
            if is_self {
                Some(classify_self_fd_str(n))
            } else {
                Some(OsFdAliasKind::OtherFd)
            }
        }
        _ => None,
    }
}

fn is_os_fd_alias_path(path: &Path) -> bool {
    classify_os_fd_alias_path(path).is_some()
}

fn classify_target_os_fd_alias(path: &Path) -> Option<OsFdAliasKind> {
    if let Some(kind) = classify_os_fd_alias_path(path) {
        return Some(kind);
    }
    let mut cur = path.to_path_buf();
    if let Some(norm) = normalize_parent_preserving_leaf(&cur) {
        cur = norm;
        if let Some(kind) = classify_os_fd_alias_path(&cur) {
            return Some(kind);
        }
    }
    for _ in 0..40 {
        let target = std::fs::read_link(&cur).ok()?;
        cur = if target.is_absolute() {
            target
        } else {
            cur.parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."))
                .join(target)
        };
        if let Some(kind) = classify_os_fd_alias_path(&cur) {
            return Some(kind);
        }
        if let Some(norm) = normalize_parent_preserving_leaf(&cur) {
            cur = norm;
            if let Some(kind) = classify_os_fd_alias_path(&cur) {
                return Some(kind);
            }
        }
    }
    None
}

fn targets_os_fd_alias(path: &Path) -> bool {
    is_os_fd_alias_path(path) || classify_target_os_fd_alias(path).is_some()
}

fn is_special_stream_or_fd_alias(path: &Path) -> bool {
    if targets_os_fd_alias(path) {
        return true;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::FileTypeExt as _;
        if let Ok(m) = std::fs::metadata(path) {
            let ft = m.file_type();
            return ft.is_char_device() || ft.is_fifo() || ft.is_socket();
        }
    }
    #[cfg(not(unix))]
    let _ = path;
    false
}

#[cfg(target_os = "linux")]
fn inspect_linux_self_fdinfo(fd: u32) -> (bool, u64) {
    // Linux `O_APPEND` bit in octal `flags:` of `/proc/self/fdinfo/<fd>`:
    // `0o0010` on Alpha/MIPS/PA-RISC/SPARC, `0o2000` on all other Linux architectures.
    #[cfg(any(
        target_arch = "mips",
        target_arch = "mips64",
        target_arch = "sparc",
        target_arch = "sparc64"
    ))]
    const LINUX_O_APPEND: u64 = 0o0010;
    #[cfg(not(any(
        target_arch = "mips",
        target_arch = "mips64",
        target_arch = "sparc",
        target_arch = "sparc64"
    )))]
    const LINUX_O_APPEND: u64 = 0o2000;

    let Ok(info) = std::fs::read_to_string(format!("/proc/self/fdinfo/{fd}")) else {
        return (false, 0);
    };
    let mut is_append = false;
    let mut pos = 0u64;
    for line in info.lines() {
        if let Some(rest) = line.strip_prefix("flags:") {
            if let Ok(flags) = u64::from_str_radix(rest.trim(), 8) {
                is_append = (flags & LINUX_O_APPEND) != 0;
            }
        } else if let Some(rest) = line.strip_prefix("pos:") {
            if let Ok(parsed_pos) = rest.trim().parse::<u64>() {
                pos = parsed_pos;
            }
        }
    }
    (is_append, pos)
}

fn open_special_or_fd_alias_for_write(target: &Path) -> std::io::Result<std::fs::File> {
    #[cfg(unix)]
    {
        match classify_target_os_fd_alias(target) {
            Some(OsFdAliasKind::SelfStdout) => std::os::fd::AsFd::as_fd(&std::io::stdout())
                .try_clone_to_owned()
                .map(std::fs::File::from),
            Some(OsFdAliasKind::SelfStderr) => std::os::fd::AsFd::as_fd(&std::io::stderr())
                .try_clone_to_owned()
                .map(std::fs::File::from),
            Some(OsFdAliasKind::SelfOtherFd(fd)) => {
                #[cfg(target_os = "linux")]
                {
                    let (is_append, pos) = inspect_linux_self_fdinfo(fd);
                    if is_append {
                        std::fs::OpenOptions::new().append(true).open(target)
                    } else {
                        let mut file = std::fs::OpenOptions::new().write(true).open(target)?;
                        if pos > 0 && file.metadata().is_ok_and(|m| m.is_file()) {
                            let _ = std::io::Seek::seek(&mut file, std::io::SeekFrom::Start(pos));
                        }
                        Ok(file)
                    }
                }
                #[cfg(not(target_os = "linux"))]
                {
                    let _ = fd;
                    std::fs::OpenOptions::new().write(true).open(target)
                }
            }
            Some(OsFdAliasKind::OtherFd) | None => {
                std::fs::OpenOptions::new().write(true).open(target)
            }
        }
    }
    #[cfg(not(unix))]
    {
        std::fs::OpenOptions::new().write(true).open(target)
    }
}

#[cfg(unix)]
fn is_dev_null_metadata(m: &std::fs::Metadata) -> bool {
    use std::os::unix::fs::{FileTypeExt as _, MetadataExt as _};
    if !m.file_type().is_char_device() {
        return false;
    }
    let Ok(dn) = std::fs::metadata("/dev/null") else {
        return false;
    };
    dn.file_type().is_char_device()
        && ((m.dev() == dn.dev() && m.ino() == dn.ino()) || m.rdev() == dn.rdev())
}

fn check_pair_collision(a: &Path, b: &Path) -> Option<bool> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt as _;
        let ma = std::fs::metadata(a).ok();
        if ma.as_ref().is_some_and(is_dev_null_metadata) {
            return Some(false);
        }
        let mb = std::fs::metadata(b).ok();
        if mb.as_ref().is_some_and(is_dev_null_metadata) {
            return Some(false);
        }
        if a == b {
            return Some(true);
        }
        if let (Some(ma), Some(mb)) = (ma, mb) {
            if ma.dev() == mb.dev() && ma.ino() == mb.ino() {
                return Some(true);
            }
        }
        None
    }
    #[cfg(not(unix))]
    {
        if a == b {
            Some(true)
        } else {
            None
        }
    }
}

fn paths_collide(a: &Path, b: &Path) -> bool {
    if let Some(decision) = check_pair_collision(a, b) {
        return decision;
    }
    if let (Some(ra), Some(rb)) = (resolve_target_path(a), resolve_target_path(b)) {
        if let Some(decision) = check_pair_collision(&ra, &rb) {
            return decision;
        }
    }
    false
}

fn push_alignment_and_companion_paths(inputs: &mut Vec<std::path::PathBuf>, bam_path: &Path) {
    inputs.push(bam_path.to_path_buf());
    inputs.push(append_os_suffix(bam_path, ".bai"));
    inputs.push(append_os_suffix(bam_path, ".csi"));
    inputs.push(append_os_suffix(bam_path, ".crai"));
    if bam_path.extension().is_some() {
        inputs.push(bam_path.with_extension("bai"));
        inputs.push(bam_path.with_extension("csi"));
        inputs.push(bam_path.with_extension("crai"));
    }
}

fn collect_protected_input_paths(config: &config::Config) -> Vec<std::path::PathBuf> {
    let mut inputs = Vec::new();
    if !config.reference.is_empty() {
        let ref_path = Path::new(&config.reference);
        inputs.push(ref_path.to_path_buf());
        inputs.push(append_os_suffix(ref_path, ".fai"));
        inputs.push(append_os_suffix(ref_path, ".gzi"));
        inputs.push(crate::reference::resolve_fasta_fai_path(ref_path));
        inputs.push(crate::reference::resolve_fasta_gzi_path(ref_path));
        if let Ok(canon_ref) = std::fs::canonicalize(ref_path) {
            if canon_ref != ref_path {
                inputs.push(append_os_suffix(&canon_ref, ".fai"));
                inputs.push(crate::reference::resolve_fasta_fai_path(&canon_ref));
                inputs.push(crate::reference::resolve_fasta_gzi_path(&canon_ref));
            }
        }
    }
    if config.region.is_none() {
        if let Some(ref bed_str) = config.bed {
            if !bed_str.is_empty() {
                inputs.push(std::path::PathBuf::from(bed_str));
            }
        }
    }
    for bam_spec in config.bam.split('|') {
        if let Ok(bams) = parse_colon_bams(bam_spec, "") {
            for bam_str in bams {
                let bam_path = Path::new(&bam_str);
                push_alignment_and_companion_paths(&mut inputs, bam_path);
                if let Ok(canon_bam) = std::fs::canonicalize(bam_path) {
                    if canon_bam != bam_path {
                        push_alignment_and_companion_paths(&mut inputs, &canon_bam);
                    }
                }
            }
        }
    }
    inputs
}

struct OutputIndexPlan {
    flag: &'static str,
    primary_path: std::path::PathBuf,
    alt_flag: &'static str,
    alt_path: std::path::PathBuf,
}

fn output_index_plan(config: &config::Config) -> Option<OutputIndexPlan> {
    let out_str = config.output_path.as_deref().filter(|&p| p != "-")?;
    let out = Path::new(out_str);
    if config.index_output {
        Some(OutputIndexPlan {
            flag: "--tbi",
            primary_path: tbi_output_path(out),
            alt_flag: "--csi",
            alt_path: csi_output_path(out),
        })
    } else if config.index_csi {
        Some(OutputIndexPlan {
            flag: "--csi",
            primary_path: csi_output_path(out),
            alt_flag: "--tbi",
            alt_path: tbi_output_path(out),
        })
    } else {
        None
    }
}

fn remove_stale_index_sidecar(alt_path: &Path) -> Result<(), VardictError> {
    match std::fs::symlink_metadata(alt_path) {
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(err) => return Err(VardictError::Io(err)),
        Ok(_) => {}
    }
    match std::fs::remove_file(alt_path) {
        Ok(()) => Ok(()),
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(err) => Err(VardictError::Io(err)),
    }
}

#[cfg(unix)]
fn check_active_default_fd_collisions(
    stream_label: &str,
    fd_meta: &std::fs::Metadata,
    protected_inputs: &[std::path::PathBuf],
    stats_opt: Option<&Path>,
) -> Result<(), VardictError> {
    use std::os::unix::fs::MetadataExt as _;

    if is_dev_null_metadata(fd_meta) {
        return Ok(());
    }
    if let Some(stats) = stats_opt {
        if let Ok(st_meta) = std::fs::metadata(stats) {
            if !is_dev_null_metadata(&st_meta)
                && fd_meta.dev() == st_meta.dev()
                && fd_meta.ino() == st_meta.ino()
            {
                return Err(VardictError::Parse(format!(
                    "{stream_label} and --stats-json '{}' must not point to the same file",
                    stats.display()
                )));
            }
        }
    }
    for input in protected_inputs {
        if let Ok(in_meta) = std::fs::metadata(input) {
            if !is_dev_null_metadata(&in_meta)
                && fd_meta.dev() == in_meta.dev()
                && fd_meta.ino() == in_meta.ino()
            {
                return Err(VardictError::Parse(format!(
                    "{stream_label} must not point to the same file as input '{}'",
                    input.display()
                )));
            }
        }
    }
    Ok(())
}

fn validate_io_paths(config: &config::Config) -> Result<(), VardictError> {
    validate_index_flags(config)?;
    if config.stats_json_path.as_deref() == Some("-") {
        return Err(VardictError::Parse(
            "--stats-json '-' (stdout) is not supported; pass a file path or /dev/null".to_string(),
        ));
    }
    let out_opt = config
        .output_path
        .as_deref()
        .filter(|&p| p != "-")
        .map(Path::new);
    let index_plan = output_index_plan(config);
    let stats_opt = config.stats_json_path.as_deref().map(Path::new);

    if let Some(out) = out_opt {
        if let Some(ref plan) = index_plan {
            if paths_collide(out, &plan.primary_path) {
                return Err(VardictError::Parse(format!(
                    "--output '{}' and its {} index '{}' must not point to the same file",
                    out.display(),
                    plan.flag,
                    plan.primary_path.display()
                )));
            }
            if paths_collide(out, &plan.alt_path) {
                return Err(VardictError::Parse(format!(
                    "--output '{}' and its alternate {} index sidecar '{}' must not point to the same file",
                    out.display(),
                    plan.alt_flag,
                    plan.alt_path.display()
                )));
            }
            if paths_collide(&plan.primary_path, &plan.alt_path) {
                return Err(VardictError::Parse(format!(
                    "--output {} index '{}' and alternate {} index sidecar '{}' must not point to the same file",
                    plan.flag,
                    plan.primary_path.display(),
                    plan.alt_flag,
                    plan.alt_path.display()
                )));
            }
        }
        if let Some(stats) = stats_opt {
            if paths_collide(out, stats) {
                return Err(VardictError::Parse(
                    "--output and --stats-json must not point to the same file".to_string(),
                ));
            }
            if let Some(ref plan) = index_plan {
                if paths_collide(&plan.primary_path, stats) {
                    return Err(VardictError::Parse(format!(
                        "--output {} index '{}' and --stats-json '{}' must not point to the same file",
                        plan.flag,
                        plan.primary_path.display(),
                        stats.display()
                    )));
                }
                if paths_collide(&plan.alt_path, stats) {
                    return Err(VardictError::Parse(format!(
                        "--output alternate {} index sidecar '{}' and --stats-json '{}' must not point to the same file",
                        plan.alt_flag,
                        plan.alt_path.display(),
                        stats.display()
                    )));
                }
            }
        }
    }

    let protected_inputs = collect_protected_input_paths(config);
    if out_opt.is_some() || index_plan.is_some() || stats_opt.is_some() {
        for input in &protected_inputs {
            if let Some(out) = out_opt {
                if paths_collide(out, input) {
                    return Err(VardictError::Parse(format!(
                        "--output '{}' must not point to the same file as input '{}'",
                        out.display(),
                        input.display()
                    )));
                }
            }
            if let Some(ref plan) = index_plan {
                if paths_collide(&plan.primary_path, input) {
                    return Err(VardictError::Parse(format!(
                        "--output {} index '{}' must not point to the same file as input '{}'",
                        plan.flag,
                        plan.primary_path.display(),
                        input.display()
                    )));
                }
                if paths_collide(&plan.alt_path, input) {
                    return Err(VardictError::Parse(format!(
                        "--output alternate {} index sidecar '{}' must not point to the same file as input '{}'",
                        plan.alt_flag,
                        plan.alt_path.display(),
                        input.display()
                    )));
                }
            }
            if let Some(stats) = stats_opt {
                if paths_collide(stats, input) {
                    return Err(VardictError::Parse(format!(
                        "--stats-json '{}' must not point to the same file as input '{}'",
                        stats.display(),
                        input.display()
                    )));
                }
            }
        }
    }

    #[cfg(unix)]
    {
        let uses_default_stderr = config.output_path.is_none()
            && config.debug_print.as_deref() == Some("ERR")
            && !config.splice;
        let uses_default_stdout = match config.output_path.as_deref() {
            Some("-") => true,
            Some(_) => false,
            None => !uses_default_stderr || config.column_header,
        };
        if uses_default_stdout {
            if let Ok(stdout_meta) = std::os::fd::AsFd::as_fd(&std::io::stdout())
                .try_clone_to_owned()
                .and_then(|fd| std::fs::File::from(fd).metadata())
            {
                let label = if config.output_path.as_deref() == Some("-") {
                    "--output '-' (stdout)"
                } else {
                    "stdout"
                };
                check_active_default_fd_collisions(
                    label,
                    &stdout_meta,
                    &protected_inputs,
                    stats_opt,
                )?;
            }
        }
        if uses_default_stderr {
            if let Ok(stderr_meta) = std::os::fd::AsFd::as_fd(&std::io::stderr())
                .try_clone_to_owned()
                .and_then(|fd| std::fs::File::from(fd).metadata())
            {
                check_active_default_fd_collisions(
                    "stderr (--DP ERR)",
                    &stderr_meta,
                    &protected_inputs,
                    stats_opt,
                )?;
            }
        }
    }

    if let Some(out) = out_opt {
        if is_directory_target(out) {
            return Err(VardictError::Io(directory_target_io_error("--output", out)));
        }
    }
    if let Some(ref plan) = index_plan {
        if is_directory_target(&plan.primary_path) {
            return Err(VardictError::Io(directory_target_io_error(
                plan.flag,
                &plan.primary_path,
            )));
        }
        if is_directory_target(&plan.alt_path) {
            return Err(VardictError::Io(directory_target_io_error(
                plan.alt_flag,
                &plan.alt_path,
            )));
        }
    }
    if let Some(stats) = stats_opt {
        if is_directory_target(stats) {
            return Err(map_stats_io_err(
                &stats.to_string_lossy(),
                directory_target_io_error("--stats-json", stats),
            ));
        }
    }
    Ok(())
}

fn map_stats_io_err(stats_path: &str, err: std::io::Error) -> VardictError {
    if err.kind() == std::io::ErrorKind::BrokenPipe {
        VardictError::Io(std::io::Error::other(format!(
            "failed to write --stats-json '{stats_path}': {err}"
        )))
    } else {
        VardictError::Io(err)
    }
}

enum PendingStatsJson {
    Staged {
        stats_label: String,
        file: std::fs::File,
        staged: StagedFileTarget,
    },
    DirectDeferred {
        stats_label: String,
        final_path: std::path::PathBuf,
    },
}

impl PendingStatsJson {
    fn open(stats_path: &str) -> Result<Self, VardictError> {
        if stats_path == "-" {
            return Err(VardictError::Parse(
                "--stats-json '-' (stdout) is not supported; pass a file path or /dev/null"
                    .to_string(),
            ));
        }
        let target = Path::new(stats_path);
        if is_directory_target(target) {
            return Err(map_stats_io_err(
                stats_path,
                directory_target_io_error("--stats-json", target),
            ));
        }
        if is_special_stream_or_fd_alias(target) {
            return Ok(Self::DirectDeferred {
                stats_label: stats_path.to_string(),
                final_path: target.to_path_buf(),
            });
        }
        let tmp_path = atomic_tmp_path(target);
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&tmp_path)
        {
            Ok(file) => Ok(Self::Staged {
                stats_label: stats_path.to_string(),
                file,
                staged: StagedFileTarget::AtomicTmp {
                    final_path: target.to_path_buf(),
                    tmp_path,
                    committed: false,
                },
            }),
            Err(err) if is_readonly_parent_fallback_eligible(target, &err) => {
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(target)
                    .map_err(|e| map_stats_io_err(stats_path, e))?;
                Ok(Self::DirectDeferred {
                    stats_label: stats_path.to_string(),
                    final_path: target.to_path_buf(),
                })
            }
            Err(err) => Err(map_stats_io_err(stats_path, err)),
        }
    }

    fn commit(self, json: &[u8]) -> Result<(), VardictError> {
        use std::io::Write as _;
        match self {
            Self::Staged {
                stats_label,
                mut file,
                staged,
            } => {
                file.write_all(json)
                    .map_err(|e| map_stats_io_err(&stats_label, e))?;
                file.flush()
                    .map_err(|e| map_stats_io_err(&stats_label, e))?;
                sync_file_ignoring_special(&file).map_err(|e| match e {
                    VardictError::Io(io_e) => map_stats_io_err(&stats_label, io_e),
                    other => other,
                })?;
                drop(file);
                staged
                    .commit()
                    .map_err(|e| map_stats_io_err(&stats_label, e))
            }
            Self::DirectDeferred {
                stats_label,
                final_path,
            } => {
                let mut file = if is_special_stream_or_fd_alias(&final_path) {
                    open_special_or_fd_alias_for_write(&final_path)
                } else {
                    let is_reg_file = std::fs::metadata(&final_path).is_ok_and(|m| m.is_file());
                    std::fs::OpenOptions::new()
                        .write(true)
                        .truncate(is_reg_file)
                        .open(&final_path)
                }
                .map_err(|e| map_stats_io_err(&stats_label, e))?;
                file.write_all(json)
                    .map_err(|e| map_stats_io_err(&stats_label, e))?;
                file.flush()
                    .map_err(|e| map_stats_io_err(&stats_label, e))?;
                sync_file_ignoring_special(&file).map_err(|e| match e {
                    VardictError::Io(io_e) => map_stats_io_err(&stats_label, io_e),
                    other => other,
                })
            }
        }
    }
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
        let started = std::time::Instant::now();
        validate_io_paths(config)?;
        let (effective, mode, preloaded_bed_lines) = resolve_effective_config(config)?;
        let pending_stats = if let Some(ref stats_path) = effective.stats_json_path {
            Some(PendingStatsJson::open(stats_path)?)
        } else {
            None
        };
        let (writer, staged_out) = open_output_writer(&effective)?;
        let telemetry = RunTelemetryCounters::default();
        let engine_res = run_with_engine_writer_and_telemetry(
            &effective,
            mode,
            preloaded_bed_lines,
            writer,
            Some(&telemetry),
        );

        let emit_stats = |pending: PendingStatsJson| -> Result<(), VardictError> {
            let elapsed_ms = started.elapsed().as_millis();
            let peak_rss_kb = read_peak_rss_kb();
            let regions_count = telemetry
                .regions_count
                .load(std::sync::atomic::Ordering::Relaxed);
            let tiles_count = telemetry
                .tiles_count
                .load(std::sync::atomic::Ordering::Relaxed);
            let exceptions_caught = effective
                .exception_counter
                .load(std::sync::atomic::Ordering::Relaxed);
            let json = format!(
                "{{\n  \"mode\": \"{}\",\n  \"threads\": {},\n  \"regions_count\": {},\n  \"tiles_count\": {},\n  \"elapsed_ms\": {},\n  \"peak_rss_kb\": {},\n  \"exceptions_caught\": {}\n}}\n",
                mode.as_str(),
                effective.threads,
                regions_count,
                tiles_count,
                elapsed_ms,
                peak_rss_kb,
                exceptions_caught,
            );
            pending.commit(json.as_bytes())
        };

        match engine_res {
            Ok(()) => {
                let index_plan = output_index_plan(&effective);
                let staged_idx = if index_plan.is_some() {
                    if let Some(out_path) = effective.output_path.as_deref().filter(|&p| p != "-") {
                        let final_vcf = Path::new(out_path);
                        let readable_vcf = match staged_out.as_ref() {
                            Some(staged) => staged.readable_path(final_vcf),
                            None => final_vcf,
                        };
                        if effective.index_output {
                            Some(stage_vcf_tabix_index(readable_vcf, final_vcf)?)
                        } else {
                            Some(stage_vcf_csi_index(readable_vcf, final_vcf)?)
                        }
                    } else {
                        None
                    }
                } else {
                    None
                };
                if let Some(staged) = staged_out {
                    staged.commit().map_err(VardictError::Io)?;
                }
                if let Some(staged_idx) = staged_idx {
                    staged_idx.commit().map_err(VardictError::Io)?;
                }
                if let Some(ref plan) = index_plan {
                    remove_stale_index_sidecar(&plan.alt_path)?;
                }
                if let Some(pending) = pending_stats {
                    emit_stats(pending)?;
                }
                Ok(())
            }
            Err(err) => {
                drop(staged_out);
                if matches!(effective.output_path.as_deref(), None | Some("-"))
                    && err.is_broken_pipe()
                {
                    if let Some(pending) = pending_stats {
                        emit_stats(pending)?;
                    }
                }
                Err(err)
            }
        }
    })
}

/// Runs the variant-calling pipeline writing all output to `writer` without clearing
/// process-wide BAM/CRAM/FASTA caches before/after the call.
pub fn run_with_writer<W: std::io::Write>(
    config: &config::Config,
    writer: W,
) -> Result<(), VardictError> {
    let (effective, mode, preloaded_bed_lines) = resolve_effective_config(config)?;
    run_with_engine_writer_and_telemetry(
        &effective,
        mode,
        preloaded_bed_lines,
        InjectedWriter(writer),
        None,
    )
}

fn run_with_engine_writer_and_telemetry<W: EngineWriter>(
    config: &config::Config,
    mode: CallingMode,
    preloaded_bed_lines: Option<Vec<String>>,
    mut writer: W,
    telemetry: Option<&RunTelemetryCounters>,
) -> Result<(), VardictError> {
    match run_with_writer_inner(config, mode, preloaded_bed_lines, &mut writer, telemetry) {
        Ok(()) => writer.finish_output(),
        Err(e) => {
            let _ = writer.flush();
            Err(e)
        }
    }
}

fn run_with_writer_inner<W: EngineWriter>(
    config: &config::Config,
    mode: CallingMode,
    preloaded_bed_lines: Option<Vec<String>>,
    mut writer: &mut W,
    telemetry: Option<&RunTelemetryCounters>,
) -> Result<(), VardictError> {
    log::info!("VarDictRust engine starting");

    if mode == CallingMode::Somatic {
        return run_somatic(config, preloaded_bed_lines, writer, telemetry);
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
        drop(preloaded_bed_lines);
        crate::regions::parse_region_flag(r_str, config, &chromosomes_lengths)?
    } else if let Some(ref bed_path) = config.bed {
        match preloaded_bed_lines {
            Some(bed_lines) => crate::regions::parse_bed_from_lines(
                &bed_lines,
                config,
                &chromosomes_lengths,
                config.zero_based_defined.then_some(config.zero_based),
            )?,
            None => crate::regions::parse_bed(
                bed_path,
                config,
                &chromosomes_lengths,
                config.zero_based_defined.then_some(config.zero_based),
            )?,
        }
    } else {
        return Err(VardictError::Parse(
            "Either region (-R) or BED file must be specified".to_string(),
        ));
    };

    // Amplicon mode: process per-amplicon and combine via amplicon post-process
    if mode == CallingMode::Amplicon {
        return run_amplicon(config, &chromosomes_lengths, regions, writer, telemetry);
    }

    let flat_regions: Vec<crate::regions::Region> = regions.into_iter().flatten().collect();

    // Tile regions into overlapping tiles (regions <= 1 Mbp execute untiled like
    // VarDictJava; multi-megabase intervals > 1 Mbp tile into 100 kbp windows).
    let tile_size = config.tile_size.unwrap_or(DEFAULT_TILE_SIZE);
    let halo = 500;
    let mut tiles = Vec::new();
    let mut multi_tile_seeds: HashMap<MultiTileKey, (Vec<MaxRlStep>, i32)> = HashMap::new();
    let reg_tiles_list: Vec<Vec<crate::regions::Tile>> = flat_regions
        .iter()
        .map(|reg| crate::regions::tile(reg, effective_tile_size(reg, config.tile_size), halo))
        .collect();
    let multi_idxs: Vec<usize> = reg_tiles_list
        .iter()
        .enumerate()
        .filter_map(|(i, rt)| (rt.len() > 1 && !config.splice).then_some(i))
        .collect();
    if multi_idxs.len() > 1 && config.threads > 1 {
        let flat_arc = std::sync::Arc::new(flat_regions.clone());
        let reg_tiles_arc = std::sync::Arc::new(reg_tiles_list.clone());
        let bam_list_arc = std::sync::Arc::new(bam_list.clone());
        let cfg_clone = config.clone();
        let chr_lens_arc = std::sync::Arc::new(chromosomes_lengths.clone());
        crate::runner::run_parallel_for_bam(
            &config.bam,
            multi_idxs,
            config.threads,
            move |idx| {
                let reg = flat_arc
                    .get(idx)
                    .ok_or_else(|| VardictError::Parse(format!("Invalid region index {idx}")))?;
                let reg_tiles = reg_tiles_arc.get(idx).ok_or_else(|| {
                    VardictError::Parse(format!("Invalid region tiles index {idx}"))
                })?;
                let (steps_per_tile, final_max_rl) = scan_multi_tile_max_read_lengths(
                    &bam_list_arc,
                    reg,
                    reg_tiles,
                    &cfg_clone,
                    &chr_lens_arc,
                    0,
                )?;
                Ok((idx, steps_per_tile, final_max_rl))
            },
            |(idx, steps_per_tile, final_max_rl)| {
                if let Some(reg_tiles) = reg_tiles_list.get(idx) {
                    for (t, steps) in reg_tiles.iter().zip(steps_per_tile) {
                        multi_tile_seeds.insert(multi_tile_key(t), (steps, final_max_rl));
                    }
                }
                Ok(())
            },
        )?;
    } else {
        for idx in multi_idxs {
            if let (Some(reg), Some(reg_tiles)) = (flat_regions.get(idx), reg_tiles_list.get(idx)) {
                let (steps_per_tile, final_max_rl) = scan_multi_tile_max_read_lengths(
                    &bam_list,
                    reg,
                    reg_tiles,
                    config,
                    &chromosomes_lengths,
                    0,
                )?;
                for (t, steps) in reg_tiles.iter().zip(steps_per_tile) {
                    multi_tile_seeds.insert(multi_tile_key(t), (steps, final_max_rl));
                }
            }
        }
    }
    for reg_tiles in reg_tiles_list {
        tiles.extend(reg_tiles);
    }
    if let Some(tel) = telemetry {
        tel.regions_count
            .store(flat_regions.len(), std::sync::atomic::Ordering::Relaxed);
        tel.tiles_count
            .store(tiles.len(), std::sync::atomic::Ordering::Relaxed);
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
        UnSplitVarBlock(String),
        Lines {
            rows: Vec<SplitVarRow>,
            is_split_tile: bool,
            tile_start: crate::types::Pos1,
            tile_end: crate::types::Pos1,
            original_start: crate::types::Pos1,
            original_end: crate::types::Pos1,
        },
        Vcf {
            chr: String,
            records: Vec<RawVcfCandidate>,
        },
    }

    let mut vcf_acc: SpillableVcfAccumulator<RawVcfCandidate> = SpillableVcfAccumulator::new();
    let mut split_var_deduper = SplitVarDeduper::default();

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

            if config_clone.splice {
                let mut rows = Vec::new();
                let mut introns: Vec<&String> = var_data.splice_count.keys().collect();
                introns.sort();
                for intron in introns {
                    if !tile_owns_splice_key(&tile, tile.original_start, tile.original_end, intron)
                    {
                        continue;
                    }
                    let count = var_data.splice_count.get(intron).map(|a| a[0]).unwrap_or(0);
                    rows.push(SplitVarRow {
                        id: SplitVarIdentity::default(),
                        is_sv: false,
                        line: format!(
                            "{}\t{}\t{}\t{}",
                            config_clone.sample, tile.chr, intron, count
                        ),
                    });
                }
                return Ok(SimpleTileOutput::Lines {
                    rows,
                    is_split_tile: false,
                    tile_start: tile.start,
                    tile_end: tile.end,
                    original_start: tile.original_start,
                    original_end: tile.original_end,
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
                if !is_split_tile {
                    let mut block = String::with_capacity(variants.len().saturating_mul(256));
                    for v in variants {
                        if config_clone.fisher {
                            crate::output::table::simple_38_append(
                                &mut block,
                                &v,
                                &builder_region,
                                &config_clone.sample,
                                &v.sv,
                            );
                        } else {
                            crate::output::table::simple_36_append(
                                &mut block,
                                &v,
                                &builder_region,
                                &config_clone.sample,
                                &v.sv,
                            );
                        }
                        if config_clone.crispr_cutting_site.unwrap_or(0) != 0 {
                            block.push('\t');
                            crate::util::push_i32(&mut block, v.crispr);
                        }
                        // ports SimpleOutputVariant.java:116-117 (DEBUG appended after crispr)
                        if config_clone.debug {
                            block.push('\t');
                            block.push_str(&v.debug);
                        }
                        block.push('\n');
                    }
                    Ok(SimpleTileOutput::UnSplitVarBlock(block))
                } else {
                    let mut rows = Vec::with_capacity(variants.len());
                    let ids = build_simple_split_var_identities(&variants, is_split_tile);
                    for (v, id) in variants.into_iter().zip(ids) {
                        let mut line = String::new();
                        if config_clone.fisher {
                            crate::output::table::simple_38_into(
                                &mut line,
                                &v,
                                &builder_region,
                                &config_clone.sample,
                                &v.sv,
                            );
                        } else {
                            crate::output::table::simple_36_into(
                                &mut line,
                                &v,
                                &builder_region,
                                &config_clone.sample,
                                &v.sv,
                            );
                        }
                        if config_clone.crispr_cutting_site.unwrap_or(0) != 0 {
                            line.push('\t');
                            crate::util::push_i32(&mut line, v.crispr);
                        }
                        // ports SimpleOutputVariant.java:116-117 (DEBUG appended after crispr)
                        if config_clone.debug {
                            line.push('\t');
                            line.push_str(&v.debug);
                        }
                        let is_sv = !v.sv.is_empty() || v.var_allele.contains('<');
                        rows.push(SplitVarRow { id, is_sv, line });
                    }
                    Ok(SimpleTileOutput::Lines {
                        rows,
                        is_split_tile,
                        tile_start: tile.start,
                        tile_end: tile.end,
                        original_start: tile.original_start,
                        original_end: tile.original_end,
                    })
                }
            } else {
                let mut buckets: std::collections::BTreeMap<
                    i32,
                    (f64, crate::variations::Variant),
                > = std::collections::BTreeMap::new();
                for v in variants {
                    let af_4dp = crate::types::round_half_even(v.frequency, 4);
                    match buckets.entry(v.start_pos) {
                        std::collections::btree_map::Entry::Vacant(e) => {
                            e.insert((af_4dp, v));
                        }
                        std::collections::btree_map::Entry::Occupied(mut e) => {
                            if af_4dp > e.get().0 {
                                e.insert((af_4dp, v));
                            }
                        }
                    }
                }
                let mut scratch = String::with_capacity(512);
                let mut records = Vec::with_capacity(buckets.len());
                for (start_pos, (af_4dp, v)) in buckets {
                    let event = crate::output::vcf::build_single_vcf_event(
                        &mut scratch,
                        &v,
                        &builder_region,
                        &config_clone,
                        None,
                    );
                    records.push(RawVcfCandidate {
                        start_pos,
                        af_4dp,
                        event,
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
                SimpleTileOutput::UnSplitVarBlock(block) => {
                    split_var_deduper
                        .flush_before_unsplit_tile(|line| write_output_line(&mut writer, line))?;
                    if !block.is_empty() {
                        std::io::Write::write_all(&mut writer, block.as_bytes())
                            .map_err(VardictError::Io)?;
                    }
                }
                SimpleTileOutput::Lines {
                    rows,
                    is_split_tile,
                    tile_start,
                    tile_end,
                    original_start,
                    original_end,
                } => {
                    split_var_deduper.accept_tile(
                        rows,
                        is_split_tile,
                        tile_start,
                        tile_end,
                        original_start,
                        original_end,
                        |line| write_output_line(&mut writer, line),
                    )?;
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
    split_var_deduper.flush_all_pending(|line| write_output_line(&mut writer, line))?;

    emit_single_vcf_accumulator(vcf_acc, config, &mut writer)?;
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
    seed_splice: rustc_hash::FxHashSet<String>,
) -> Result<
    (
        crate::javamap::JavaHashMap<i32, crate::variations::Vars>,
        rustc_hash::FxHashSet<String>,
    ),
    VardictError,
> {
    let halo = 500u32;
    let halo_start = crate::types::Pos0(region.start.0.saturating_sub(halo).saturating_sub(1));
    let halo_end = crate::types::Pos0(region.end.0.saturating_add(halo));

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
    telemetry: Option<&RunTelemetryCounters>,
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

    type AmpliconVarsList = Vec<crate::javamap::JavaHashMap<i32, crate::variations::Vars>>;

    enum AmpliconSegmentOutput {
        UnSplitVarBlock(String),
        Vcf {
            chr: String,
            records: Vec<RawVcfCandidate>,
        },
        Pileup {
            segment: Vec<crate::regions::Region>,
            vars_per_amplicon: AmpliconVarsList,
            splice: rustc_hash::FxHashSet<String>,
        },
    }

    let segments: Vec<Vec<crate::regions::Region>> =
        segments.into_iter().filter(|s| !s.is_empty()).collect();
    if let Some(tel) = telemetry {
        let total_regions: usize = segments.iter().map(Vec::len).sum();
        tel.regions_count
            .store(total_regions, std::sync::atomic::Ordering::Relaxed);
        tel.tiles_count
            .store(segments.len(), std::sync::atomic::Ordering::Relaxed);
    }
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
                return Ok(AmpliconSegmentOutput::UnSplitVarBlock(String::new()));
            };
            let run_segment_regions = || {
                let mut vars_per_amplicon: AmpliconVarsList = Vec::with_capacity(segment.len());
                let mut splice = rustc_hash::FxHashSet::default();
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
                Vec<(usize, &crate::regions::Region)>,
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
                let covering = sweep.covering_refs_at(p);
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
                let mut block = String::with_capacity(results.len().saturating_mul(256));
                for av in results {
                    if config_clone.fisher {
                        crate::output::table::amplicon_40_append(
                            &mut block,
                            &av,
                            current_region,
                            &config_clone.sample,
                            config_clone.debug,
                        );
                    } else {
                        crate::output::table::amplicon_38_append(
                            &mut block,
                            &av,
                            current_region,
                            &config_clone.sample,
                            config_clone.debug,
                        );
                    }
                    block.push('\n');
                }
                Ok(AmpliconSegmentOutput::UnSplitVarBlock(block))
            } else {
                let mut buckets: std::collections::BTreeMap<
                    i32,
                    (f64, crate::postprocess::amplicon::AmpliconVariant),
                > = std::collections::BTreeMap::new();
                for av in results {
                    let af_4dp = crate::types::round_half_even(av.variant.frequency, 4);
                    match buckets.entry(av.variant.start_pos) {
                        std::collections::btree_map::Entry::Vacant(e) => {
                            e.insert((af_4dp, av));
                        }
                        std::collections::btree_map::Entry::Occupied(mut e) => {
                            if af_4dp > e.get().0 {
                                e.insert((af_4dp, av));
                            }
                        }
                    }
                }
                let mut scratch = String::with_capacity(512);
                let mut records = Vec::with_capacity(buckets.len());
                for (start_pos, (af_4dp, av)) in buckets {
                    let event = crate::output::vcf::build_single_vcf_event(
                        &mut scratch,
                        &av.variant,
                        current_region,
                        &config_clone,
                        Some((av.gdamp, av.tlamp, av.ncamp, av.ampflag)),
                    );
                    records.push(RawVcfCandidate {
                        start_pos,
                        af_4dp,
                        event,
                    });
                }
                Ok(AmpliconSegmentOutput::Vcf {
                    chr: current_region.chr.clone(),
                    records,
                })
            }
        },
        |output| {
            match output {
                AmpliconSegmentOutput::UnSplitVarBlock(block) => {
                    if !block.is_empty() {
                        std::io::Write::write_all(writer, block.as_bytes())
                            .map_err(VardictError::Io)?;
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
                            crate::postprocess::amplicon::amplicon_position_ref_stream(&segment),
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
                        let mut scratch = String::with_capacity(512);
                        crate::postprocess::amplicon::amplicon_post_process_streaming(
                            config,
                            current_region,
                            &vars_per_amplicon,
                            crate::postprocess::amplicon::amplicon_position_ref_stream(&segment),
                            &splice,
                            &mut |av| {
                                let start_pos = av.variant.start_pos;
                                let af_4dp = crate::types::round_half_even(av.variant.frequency, 4);
                                let event = crate::output::vcf::build_single_vcf_event(
                                    &mut scratch,
                                    &av.variant,
                                    current_region,
                                    config,
                                    Some((av.gdamp, av.tlamp, av.ncamp, av.ampflag)),
                                );
                                vcf_acc.insert(
                                    &current_region.chr,
                                    RawVcfCandidate {
                                        start_pos,
                                        af_4dp,
                                        event,
                                    },
                                )
                            },
                        )?;
                    }
                }
            }
            Ok(())
        },
    )?;

    emit_single_vcf_accumulator(vcf_acc, config, writer)?;

    writer.flush().map_err(VardictError::Io)?;

    Ok(())
}

/// Somatic (paired) mode: process tumor BAM, then normal BAM, then compare.
/// Ports SomaticMode.java's processBothBamsInPipeline.
fn run_somatic<W: EngineWriter>(
    config: &config::Config,
    preloaded_bed_lines: Option<Vec<String>>,
    writer: &mut W,
    telemetry: Option<&RunTelemetryCounters>,
) -> Result<(), VardictError> {
    let mut parts = config.bam.split('|');
    let (Some(bam1_spec), Some(bam2_spec), None) = (parts.next(), parts.next(), parts.next())
    else {
        return Err(VardictError::Parse(
            "Somatic mode requires exactly two BAM files separated by '|'".to_string(),
        ));
    };
    let bam1_bams = parse_colon_bams(bam1_spec, "No tumor BAM file specified")?;
    let bam2_bams = parse_colon_bams(bam2_spec, "No normal BAM file specified")?;
    let first_bam1 = bam1_bams
        .first()
        .ok_or_else(|| VardictError::Parse("No tumor BAM file specified".to_string()))?;

    // Read reference lengths from cached BAM/CRAM header
    let chromosomes_lengths = crate::bam::read_reference_lengths(Path::new(first_bam1), config)?;

    // Parse target regions
    let regions = if let Some(ref r_str) = config.region {
        drop(preloaded_bed_lines);
        crate::regions::parse_region_flag(r_str, config, &chromosomes_lengths)?
    } else if let Some(ref bed_path) = config.bed {
        match preloaded_bed_lines {
            Some(bed_lines) => crate::regions::parse_bed_from_lines(
                &bed_lines,
                config,
                &chromosomes_lengths,
                config.zero_based_defined.then_some(config.zero_based),
            )?,
            None => crate::regions::parse_bed(
                bed_path,
                config,
                &chromosomes_lengths,
                config.zero_based_defined.then_some(config.zero_based),
            )?,
        }
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
    let reg_tiles_list: Vec<Vec<crate::regions::Tile>> = flat_regions
        .iter()
        .map(|reg| crate::regions::tile(reg, effective_tile_size(reg, config.tile_size), halo))
        .collect();
    let multi_idxs: Vec<usize> = reg_tiles_list
        .iter()
        .enumerate()
        .filter_map(|(i, rt)| (rt.len() > 1).then_some(i))
        .collect();
    if multi_idxs.len() > 1 && config.threads > 1 {
        let flat_arc = std::sync::Arc::new(flat_regions.clone());
        let reg_tiles_arc = std::sync::Arc::new(reg_tiles_list.clone());
        let bam1_arc = std::sync::Arc::new(bam1_bams.clone());
        let bam2_arc = std::sync::Arc::new(bam2_bams.clone());
        let cfg_clone = config.clone();
        let chr_lens_arc = std::sync::Arc::new(chromosomes_lengths.clone());
        crate::runner::run_parallel_for_bam(
            &config.bam,
            multi_idxs,
            config.threads,
            move |idx| {
                let reg = flat_arc
                    .get(idx)
                    .ok_or_else(|| VardictError::Parse(format!("Invalid region index {idx}")))?;
                let reg_tiles = reg_tiles_arc.get(idx).ok_or_else(|| {
                    VardictError::Parse(format!("Invalid region tiles index {idx}"))
                })?;
                let (bam1_steps, bam1_final) = scan_multi_tile_max_read_lengths(
                    &bam1_arc,
                    reg,
                    reg_tiles,
                    &cfg_clone,
                    &chr_lens_arc,
                    0,
                )?;
                let (bam2_steps, bam2_final) = scan_multi_tile_max_read_lengths(
                    &bam2_arc,
                    reg,
                    reg_tiles,
                    &cfg_clone,
                    &chr_lens_arc,
                    bam1_final,
                )?;
                Ok((idx, bam1_steps, bam1_final, bam2_steps, bam2_final))
            },
            |(idx, bam1_steps, bam1_final, bam2_steps, bam2_final)| {
                if let Some(reg_tiles) = reg_tiles_list.get(idx) {
                    for ((t, s1), s2) in reg_tiles.iter().zip(bam1_steps).zip(bam2_steps) {
                        multi_tile_seeds
                            .insert(multi_tile_key(t), (s1, bam1_final, s2, bam2_final));
                    }
                }
                Ok(())
            },
        )?;
    } else {
        for idx in multi_idxs {
            if let (Some(reg), Some(reg_tiles)) = (flat_regions.get(idx), reg_tiles_list.get(idx)) {
                let (bam1_steps, bam1_final) = scan_multi_tile_max_read_lengths(
                    &bam1_bams,
                    reg,
                    reg_tiles,
                    config,
                    &chromosomes_lengths,
                    0,
                )?;
                let (bam2_steps, bam2_final) = scan_multi_tile_max_read_lengths(
                    &bam2_bams,
                    reg,
                    reg_tiles,
                    config,
                    &chromosomes_lengths,
                    bam1_final,
                )?;
                for ((t, s1), s2) in reg_tiles.iter().zip(bam1_steps).zip(bam2_steps) {
                    multi_tile_seeds.insert(multi_tile_key(t), (s1, bam1_final, s2, bam2_final));
                }
            }
        }
    }
    for reg_tiles in reg_tiles_list {
        tiles.extend(reg_tiles);
    }
    if let Some(tel) = telemetry {
        tel.regions_count
            .store(flat_regions.len(), std::sync::atomic::Ordering::Relaxed);
        tel.tiles_count
            .store(tiles.len(), std::sync::atomic::Ordering::Relaxed);
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
        UnSplitVarBlock(String),
        Lines {
            rows: Vec<SplitVarRow>,
            is_split_tile: bool,
            tile_start: crate::types::Pos1,
            tile_end: crate::types::Pos1,
            original_start: crate::types::Pos1,
            original_end: crate::types::Pos1,
        },
        Vcf {
            chr: String,
            candidates: Vec<SomaticVcfCandidate>,
        },
    }

    let mut vcf_acc: SpillableVcfAccumulator<SomaticVcfCandidate> = SpillableVcfAccumulator::new();
    let mut split_var_deduper = SplitVarDeduper::default();

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
            let use_compact_somatic_ref_mode =
                !config_clone.pileup && !config_clone.debug && config_clone.amplicon.is_none();
            let mut ref_sites1: rustc_hash::FxHashMap<i32, crate::tovars::CompactRefSite> =
                rustc_hash::FxHashMap::default();
            let mut aligned1 = {
                let builder1 = crate::tovars::ToVarsBuilder::new(
                    &reference,
                    &builder_region,
                    &chromosomes_lengths_clone,
                    &config_clone,
                );
                if use_compact_somatic_ref_mode {
                    builder1.build_aligned_variants_with_ref_mode(
                        &var_data1,
                        crate::tovars::RefSiteMode::CompactOut(&mut ref_sites1),
                    )
                } else {
                    builder1.build_aligned_variants(&var_data1)
                }
            };
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
            let aligned2 = if use_compact_somatic_ref_mode {
                let aligned2 = builder2.build_aligned_variants_with_ref_mode(
                    &var_data2,
                    crate::tovars::RefSiteMode::KeepWhere(&|p| {
                        aligned1.get(&p).is_some_and(|v1| !v1.variants.is_empty())
                    }),
                );
                for (&pos, v2) in aligned2.iter() {
                    if !v2.variants.is_empty() && !aligned1.contains_key(&pos) {
                        if let Some(site) = ref_sites1.get(&pos) {
                            aligned1.insert(pos, builder2.hydrate_compact_ref_site(pos, site));
                        }
                    }
                }
                drop(ref_sites1);
                aligned2
            } else {
                drop(ref_sites1);
                builder2.build_aligned_variants(&var_data2)
            };
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
            let pipeline_fn =
                |combine_region: &crate::regions::Region,
                 splice_set: &mut rustc_hash::FxHashSet<String>,
                 max_rl: &mut i32,
                 target_pos: i32|
                 -> Result<Option<crate::variations::Vars>, VardictError> {
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

                    // ToVarsBuilder (single-position finalization at target_pos)
                    let builder = crate::tovars::ToVarsBuilder::new(
                        &pipeline_ref,
                        combine_region,
                        &chr_lens_for_combine,
                        &config_for_combine,
                    );
                    if use_compact_somatic_ref_mode {
                        Ok(builder.build_vars_at_position_with_ref_mode(
                            &var_data,
                            target_pos,
                            crate::tovars::RefSiteMode::DropAll,
                        ))
                    } else {
                        Ok(builder.build_vars_at_position(&var_data, target_pos))
                    }
                };

            let somatic_variants = crate::postprocess::somatic::somatic_post_process_at(
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
                if !is_split_tile {
                    let mut block =
                        String::with_capacity(somatic_variants.len().saturating_mul(384));
                    for sv in somatic_variants {
                        if config_clone.fisher {
                            crate::output::table::somatic_61_append(
                                &mut block,
                                &sv.variant1,
                                sv.variant2.as_ref(),
                                &builder_region,
                                &config_clone.sample,
                                &sv.status,
                                &sv.sv1,
                                &sv.sv2,
                                sv.end_is_variant2,
                            );
                        } else {
                            crate::output::table::somatic_55_append(
                                &mut block,
                                &sv.variant1,
                                sv.variant2.as_ref(),
                                &builder_region,
                                &config_clone.sample,
                                &sv.status,
                                &sv.sv1,
                                &sv.sv2,
                                sv.end_is_variant2,
                            );
                        }
                        // ports SomaticOutputVariant.java:193-194 (DEBUG = tumor variant's)
                        if config_clone.debug {
                            block.push('\t');
                            block.push_str(&sv.variant1.debug);
                        }
                        block.push('\n');
                    }
                    Ok(SomaticTileOutput::UnSplitVarBlock(block))
                } else {
                    let mut rows = Vec::with_capacity(somatic_variants.len());
                    for sv in somatic_variants {
                        let mut line = String::new();
                        if config_clone.fisher {
                            crate::output::table::somatic_61_into(
                                &mut line,
                                &sv.variant1,
                                sv.variant2.as_ref(),
                                &builder_region,
                                &config_clone.sample,
                                &sv.status,
                                &sv.sv1,
                                &sv.sv2,
                                sv.end_is_variant2,
                            );
                        } else {
                            crate::output::table::somatic_55_into(
                                &mut line,
                                &sv.variant1,
                                sv.variant2.as_ref(),
                                &builder_region,
                                &config_clone.sample,
                                &sv.status,
                                &sv.sv1,
                                &sv.sv2,
                                sv.end_is_variant2,
                            );
                        }
                        // ports SomaticOutputVariant.java:193-194 (DEBUG = tumor variant's)
                        if config_clone.debug {
                            line.push('\t');
                            line.push_str(&sv.variant1.debug);
                        }
                        let id = SplitVarIdentity {
                            origin_pos: sv.variant1.start_pos,
                            start_pos: sv.variant1.start_pos,
                            end_pos: sv.variant1.end_pos,
                            ref_allele: sv.variant1.ref_allele.clone(),
                            var_allele: sv.variant1.var_allele.clone(),
                            null_record: sv.variant1.null_record,
                            companion_kind: 0,
                        };
                        let is_sv = !sv.sv1.is_empty()
                            || !sv.sv2.is_empty()
                            || sv.variant1.var_allele.contains('<');
                        rows.push(SplitVarRow { id, is_sv, line });
                    }
                    Ok(SomaticTileOutput::Lines {
                        rows,
                        is_split_tile,
                        tile_start: tile.start,
                        tile_end: tile.end,
                        original_start: tile.original_start,
                        original_end: tile.original_end,
                    })
                }
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
                SomaticTileOutput::UnSplitVarBlock(block) => {
                    split_var_deduper
                        .flush_before_unsplit_tile(|line| write_output_line(writer, line))?;
                    if !block.is_empty() {
                        std::io::Write::write_all(writer, block.as_bytes())
                            .map_err(VardictError::Io)?;
                    }
                }
                SomaticTileOutput::Lines {
                    rows,
                    is_split_tile,
                    tile_start,
                    tile_end,
                    original_start,
                    original_end,
                } => {
                    split_var_deduper.accept_tile(
                        rows,
                        is_split_tile,
                        tile_start,
                        tile_end,
                        original_start,
                        original_end,
                        |line| write_output_line(writer, line),
                    )?;
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
    split_var_deduper.flush_all_pending(|line| write_output_line(writer, line))?;

    let mut cur_chr: Option<String> = None;
    let mut streamer = crate::output::vcf::PairedRecordStreamer::from_config(config);
    vcf_acc.emit_all(|chr, cand| {
        if cur_chr.as_deref() != Some(chr) {
            if cur_chr.is_some() {
                let prev = std::mem::replace(
                    &mut streamer,
                    crate::output::vcf::PairedRecordStreamer::from_config(config),
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
        splice_count: &rustc_hash::FxHashMap<String, [i32; 1]>,
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

    fn mk_test_single_vcf_rec(
        chr: &str,
        pos: i32,
        ref_a: &str,
        alt_a: &str,
        filter: &str,
        info: &str,
        type_snv: bool,
        pre_cluster_non_amp_empty: bool,
        has_ampbias: bool,
    ) -> crate::output::vcf::SingleVcfRecord {
        let prefix = format!("{chr}\t{pos}\t.\t{ref_a}\t{alt_a}\t35\t");
        let filter_start = prefix.len();
        let filter_end = filter_start + filter.len();
        let line = format!("{prefix}{filter}\t{info}\tGT\t0/1");
        crate::output::vcf::SingleVcfRecord::from_parts_checked(
            line,
            filter_start,
            filter_end,
            pos,
            type_snv,
            pre_cluster_non_amp_empty,
            has_ampbias,
        )
        .expect("valid 9-tab SingleVcfRecord")
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
                event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                    "chr2", 500, "A", "G", "PASS", "first_25", true, true, false,
                )),
            },
        )
        .unwrap();
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 100,
                af_4dp: 0.5000,
                event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                    "chr1", 100, "C", "T", "PASS", "first_50", true, true, false,
                )),
            },
        )
        .unwrap();
        // Equal AF at chr1:100 in a later spill run must NOT overwrite first_50
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 100,
                af_4dp: 0.5000,
                event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                    "chr1",
                    100,
                    "C",
                    "A",
                    "PASS",
                    "second_tie_50",
                    true,
                    true,
                    false,
                )),
            },
        )
        .unwrap();
        // Third spill run on chr1 with equal AF at 100 so > max_fanin(2) runs merge in passes
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 100,
                af_4dp: 0.5000,
                event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                    "chr1",
                    100,
                    "C",
                    "G",
                    "PASS",
                    "third_tie_50",
                    true,
                    true,
                    false,
                )),
            },
        )
        .unwrap();
        // Higher AF at chr2:500 in a later spill run MUST overwrite first_25
        acc.insert(
            "chr2",
            RawVcfCandidate {
                start_pos: 500,
                af_4dp: 0.7500,
                event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                    "chr2", 500, "A", "T", "PASS", "later_75", true, true, false,
                )),
            },
        )
        .unwrap();
        acc.insert(
            "chr1",
            RawVcfCandidate {
                start_pos: 50,
                af_4dp: 0.3000,
                event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                    "chr1",
                    50,
                    "G",
                    "C",
                    "PASS",
                    "early_pos",
                    true,
                    true,
                    false,
                )),
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
            if let crate::output::vcf::SingleVcfEvent::Record(rec) = cand.event {
                emitted.push(rec.as_str().to_string());
            }
            Ok(())
        })
        .unwrap();

        assert_eq!(
            emitted,
            vec![
                "chr1\t50\t.\tG\tC\t35\tPASS\tearly_pos\tGT\t0/1".to_string(),
                "chr1\t100\t.\tC\tT\t35\tPASS\tfirst_50\tGT\t0/1".to_string(),
                "chr2\t500\t.\tA\tT\t35\tPASS\tlater_75\tGT\t0/1".to_string(),
            ]
        );

        check_somatic_candidate_spill_and_winner_suppression();

        // Verify exact `current_bytes` accounting across in-memory `Entry::Occupied` transitions:
        // initial `Record`, equal/lower-AF no-op, `Record -> Skip` shrinkage, and `Skip -> Record` growth.
        let mut mem_acc: SpillableVcfAccumulator<RawVcfCandidate> =
            SpillableVcfAccumulator::with_limits(64 * 1024, 8);
        let skip_unit_bytes = std::mem::size_of::<RawVcfCandidate>() + 32;
        let rec1 = mk_test_single_vcf_rec(
            "chr1",
            300,
            "A",
            "G",
            "PASS",
            "long_info_payload_initial_record",
            true,
            true,
            false,
        );
        let rec1_cap = rec1.capacity();
        mem_acc
            .insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 300,
                    af_4dp: 0.4000,
                    event: crate::output::vcf::SingleVcfEvent::Record(rec1),
                },
            )
            .unwrap();
        assert_eq!(mem_acc.current_bytes, skip_unit_bytes + rec1_cap);

        // Equal AF and lower AF are no-ops on `current_bytes`
        for noop_af in [0.4000, 0.2000] {
            mem_acc
                .insert(
                    "chr1",
                    RawVcfCandidate {
                        start_pos: 300,
                        af_4dp: noop_af,
                        event: crate::output::vcf::SingleVcfEvent::Skip,
                    },
                )
                .unwrap();
            assert_eq!(mem_acc.current_bytes, skip_unit_bytes + rec1_cap);
        }

        // Higher-AF `Skip` shrinks `current_bytes` by `rec1_cap`
        mem_acc
            .insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 300,
                    af_4dp: 0.6000,
                    event: crate::output::vcf::SingleVcfEvent::Skip,
                },
            )
            .unwrap();
        assert_eq!(mem_acc.current_bytes, skip_unit_bytes);

        // Even-higher-AF `Record` grows `current_bytes` by `rec2_cap`
        let rec2 = mk_test_single_vcf_rec(
            "chr1",
            300,
            "A",
            "T",
            "PASS",
            "even_longer_info_payload_replacement_record_after_skip",
            true,
            true,
            false,
        );
        let rec2_cap = rec2.capacity();
        mem_acc
            .insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 300,
                    af_4dp: 0.9000,
                    event: crate::output::vcf::SingleVcfEvent::Record(rec2),
                },
            )
            .unwrap();
        assert_eq!(mem_acc.current_bytes, skip_unit_bytes + rec2_cap);
    }

    fn check_somatic_candidate_spill_and_winner_suppression() {
        // 1. Position-winner suppression when highest-AF winner is `SingleVcfEvent::Skip`
        //    (`var2vcf_valid.pl:135-205` — `--pass-only` and `--min-splitreads` run AFTER `$tmp[0]`)
        let mut raw_acc: SpillableVcfAccumulator<RawVcfCandidate> =
            SpillableVcfAccumulator::with_limits(50, 2);
        raw_acc
            .insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 200,
                    af_4dp: 0.8000,
                    event: crate::output::vcf::SingleVcfEvent::Skip,
                },
            )
            .unwrap();
        raw_acc
            .insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 200,
                    af_4dp: 0.3000,
                    event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr1",
                        200,
                        "A",
                        "G",
                        "PASS",
                        "lower_af_pass",
                        true,
                        true,
                        false,
                    )),
                },
            )
            .unwrap();
        raw_acc
            .insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 250,
                    af_4dp: 0.2000,
                    event: crate::output::vcf::SingleVcfEvent::Skip,
                },
            )
            .unwrap();
        raw_acc
            .insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 250,
                    af_4dp: 0.9000,
                    event: crate::output::vcf::SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr1",
                        250,
                        "C",
                        "T",
                        "PASS",
                        "higher_af_pass",
                        true,
                        true,
                        false,
                    )),
                },
            )
            .unwrap();
        let mut raw_emitted = Vec::new();
        raw_acc
            .emit_all(|_chr, cand| {
                if let crate::output::vcf::SingleVcfEvent::Record(rec) = cand.event {
                    raw_emitted.push(rec.as_str().to_string());
                }
                Ok(())
            })
            .unwrap();
        assert_eq!(
            raw_emitted,
            vec!["chr1\t250\t.\tC\tT\t35\tPASS\thigher_af_pass\tGT\t0/1".to_string()],
            "filtered higher-AF Skip tombstone at 200 must suppress lower-AF PASS candidate across spill runs"
        );

        // 2. Single-buffer `SomaticVcfCandidate` u64-offset spill round-trip + k-way merge
        let cfg = crate::config::Config::default();
        let mk_som = |pos: i32, af: f64, alt: &str| -> SomaticVcfCandidate {
            let mut v = crate::variations::Variant::new();
            v.ref_allele = "A".to_string();
            v.var_allele = alt.to_string();
            v.start_pos = pos;
            v.end_pos = pos;
            v.total_pos_cov = 50;
            v.pos_coverage = (50.0 * af).round() as i32;
            v.ref_fwd_cov = 20;
            v.ref_rev_cov = 20;
            v.fwd = 5;
            v.rev = 5;
            v.frequency = af;
            v.mean_position = 25.0;
            v.pstd = true;
            v.mean_quality = 35.0;
            v.qstd = true;
            v.mean_mapping_quality = 60.0;
            v.hq_lq_ratio = 2.0;
            v.var_type = "SNV".to_string();
            let sv = crate::postprocess::somatic::SomaticVariant {
                variant1: v,
                variant2: None,
                status: "StrongSomatic".to_string(),
                sv1: String::new(),
                sv2: String::new(),
                end_is_variant2: false,
            };
            SomaticVcfCandidate {
                start_pos: pos,
                af_4dp: crate::types::round_half_even(af, 4),
                record: Some(crate::output::vcf::format_paired_vcf_record(
                    &sv, "chr1", &cfg,
                )),
            }
        };

        let mut som_acc: SpillableVcfAccumulator<SomaticVcfCandidate> =
            SpillableVcfAccumulator::with_limits(128, 2);
        som_acc.insert("chr1", mk_som(100, 0.40, "G")).unwrap();
        som_acc.insert("chr1", mk_som(100, 0.70, "T")).unwrap();
        som_acc
            .insert(
                "chr1",
                SomaticVcfCandidate {
                    start_pos: 150,
                    af_4dp: 0.9000,
                    record: None,
                },
            )
            .unwrap();
        som_acc.insert("chr1", mk_som(150, 0.30, "C")).unwrap();

        let mut som_lines = Vec::new();
        som_acc
            .emit_all(|_chr, cand| {
                if let Some(rec) = cand.record {
                    assert_eq!(rec.filter_str(), "PASS");
                    som_lines.push(rec.as_str().to_string());
                }
                Ok(())
            })
            .unwrap();
        assert_eq!(som_lines.len(), 1);
        assert!(som_lines[0].starts_with("chr1\t100\t.\tA\tT\t"));

        // 3. Corrupt `u64` filter spans, oversized/truncated string lengths, and corrupt boolean bytes
        //    in `SomaticVcfCandidate::read_from` / `RawVcfCandidate::read_from` fail cleanly with zero panic
        for (bad_start, bad_end) in [(10u64, 5u64), (0u64, 9999u64), (u64::MAX, u64::MAX)] {
            let mut corrupt = Vec::new();
            write_spill_i32(&mut corrupt, 100).unwrap();
            write_spill_f64(&mut corrupt, 0.5).unwrap();
            write_spill_u8(&mut corrupt, 1).unwrap();
            write_spill_str(
                &mut corrupt,
                "chr1\t100\t.\tA\tG\t35\tPASS\tSTATUS=Somatic\tGT\t0/1\t0/0",
            )
            .unwrap();
            write_spill_u64(&mut corrupt, bad_start).unwrap();
            write_spill_u64(&mut corrupt, bad_end).unwrap();
            write_spill_i32(&mut corrupt, 100).unwrap();
            write_spill_u8(&mut corrupt, 1).unwrap();
            write_spill_u8(&mut corrupt, 1).unwrap();
            let err = SomaticVcfCandidate::read_from(&mut corrupt.as_slice())
                .err()
                .expect("corrupt filter span must fail SomaticVcfCandidate::read_from");
            assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        }

        for (bad_len, expected_kind) in [
            (u64::MAX, std::io::ErrorKind::InvalidData),
            (1u64 << 40, std::io::ErrorKind::InvalidData),
            (1024u64 * 1024, std::io::ErrorKind::UnexpectedEof),
        ] {
            let mut corrupt_len = Vec::new();
            write_spill_i32(&mut corrupt_len, 100).unwrap();
            write_spill_f64(&mut corrupt_len, 0.5).unwrap();
            write_spill_u8(&mut corrupt_len, 1).unwrap();
            write_spill_u64(&mut corrupt_len, bad_len).unwrap();
            let err = SomaticVcfCandidate::read_from(&mut corrupt_len.as_slice())
                .err()
                .expect("corrupt/truncated spill string length must fail without panic");
            assert_eq!(err.kind(), expected_kind);

            let mut corrupt_raw_len = Vec::new();
            write_spill_i32(&mut corrupt_raw_len, 100).unwrap();
            write_spill_f64(&mut corrupt_raw_len, 0.5).unwrap();
            write_spill_u8(&mut corrupt_raw_len, 3).unwrap();
            write_spill_u64(&mut corrupt_raw_len, bad_len).unwrap();
            let raw_err = RawVcfCandidate::read_from(&mut corrupt_raw_len.as_slice())
                .err()
                .expect(
                    "corrupt/truncated spill string length must fail RawVcfCandidate::read_from",
                );
            assert_eq!(raw_err.kind(), expected_kind);
        }

        let mut bad_bool = Vec::new();
        write_spill_i32(&mut bad_bool, 100).unwrap();
        write_spill_f64(&mut bad_bool, 0.5).unwrap();
        write_spill_u8(&mut bad_bool, 2).unwrap();
        assert_eq!(
            SomaticVcfCandidate::read_from(&mut bad_bool.as_slice())
                .err()
                .expect("corrupt boolean byte must fail")
                .kind(),
            std::io::ErrorKind::InvalidData
        );
    }

    #[test]
    fn test_raw_vcf_candidate_spill_roundtrip_all_variants_and_corrupt_rejection() {
        use crate::output::vcf::{SingleVcfEvent, SingleVcfRecord};

        let rec_pass_snv =
            mk_test_single_vcf_rec("chr1", 120, "A", "G", "PASS", "DP=50", true, true, false);
        let rec_amp_indel = mk_test_single_vcf_rec(
            "chr1",
            180,
            "AC",
            "A",
            "q22.5;AMPBIAS",
            "DP=40",
            false,
            false,
            true,
        );

        let items = [
            RawVcfCandidate {
                start_pos: 100,
                af_4dp: 0.8500,
                event: SingleVcfEvent::Skip,
            },
            RawVcfCandidate {
                start_pos: 110,
                af_4dp: 0.4500,
                event: SingleVcfEvent::AmpBiasSnvClusterProbe { start: 110 },
            },
            RawVcfCandidate {
                start_pos: 115,
                af_4dp: 0.6000,
                event: SingleVcfEvent::SvSplitreadFlush {
                    cluster_probe_start: Some(115),
                },
            },
            RawVcfCandidate {
                start_pos: 120,
                af_4dp: 0.5000,
                event: SingleVcfEvent::Record(rec_pass_snv),
            },
            RawVcfCandidate {
                start_pos: 180,
                af_4dp: 0.3500,
                event: SingleVcfEvent::Record(rec_amp_indel),
            },
        ];

        assert_eq!(
            items[0].approx_heap_bytes(),
            std::mem::size_of::<RawVcfCandidate>()
        );
        assert_eq!(
            items[1].approx_heap_bytes(),
            std::mem::size_of::<RawVcfCandidate>()
        );
        assert_eq!(
            items[2].approx_heap_bytes(),
            std::mem::size_of::<RawVcfCandidate>()
        );
        let SingleVcfEvent::Record(ref stored_rec) = items[3].event else {
            panic!("expected Record");
        };
        assert_eq!(
            items[3].approx_heap_bytes(),
            std::mem::size_of::<RawVcfCandidate>() + stored_rec.capacity()
        );

        let mut buf = Vec::new();
        for item in &items {
            item.write_to(&mut buf).unwrap();
        }
        let mut slice = buf.as_slice();
        for expected in &items {
            let decoded = RawVcfCandidate::read_from(&mut slice)
                .unwrap()
                .expect("expected decoded RawVcfCandidate");
            assert_eq!(decoded.start_pos, expected.start_pos);
            assert_eq!(decoded.af_4dp.to_bits(), expected.af_4dp.to_bits());
            assert_eq!(decoded.event, expected.event);
        }
        assert!(RawVcfCandidate::read_from(&mut slice).unwrap().is_none());

        // Reject unknown event tags (> 3)
        for bad_tag in [4u8, 7u8, 255u8] {
            let mut bad = Vec::new();
            write_spill_i32(&mut bad, 100).unwrap();
            write_spill_f64(&mut bad, 0.5).unwrap();
            write_spill_u8(&mut bad, bad_tag).unwrap();
            let err = RawVcfCandidate::read_from(&mut bad.as_slice())
                .err()
                .expect("unknown SingleVcfEvent tag must fail");
            assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        }

        // Reject corrupt SingleVcfRecord filter spans and malformed 9-tab invariants
        let valid_line = "chr1\t120\t.\tA\tG\t35\tPASS\tDP=50\tGT\t0/1";
        let valid_fs = "chr1\t120\t.\tA\tG\t35\t".len() as u64;
        let valid_fe = valid_fs + 4;
        for (line, fs, fe) in [
            (valid_line, 25u64, 18u64),
            (valid_line, 0u64, 9999u64),
            (valid_line, u64::MAX, u64::MAX),
            (valid_line, valid_fs, valid_fs), // empty filter
            (valid_line, 5u64, 6u64),         // wrong column (< 6 preceding tabs)
            (
                "chr1\t120\t.\tA\tG\t35\tPASS\tDP=50\tGT",
                valid_fs,
                valid_fe,
            ), // 8 tabs instead of 9
            (
                "chr1\t120\t.\tA\tG\t35\tPASS\tDP=50\tGT\t0/1\textra",
                valid_fs,
                valid_fe,
            ), // 10 tabs instead of 9
            (
                "chr1\t120\t.\tA\tG\t35\tPA\nSS\tDP=50\tGT\t0/1",
                valid_fs,
                valid_fs + 5,
            ), // embedded newline
        ] {
            let mut corrupt = Vec::new();
            write_spill_i32(&mut corrupt, 120).unwrap();
            write_spill_f64(&mut corrupt, 0.5).unwrap();
            write_spill_u8(&mut corrupt, 3).unwrap();
            write_spill_str(&mut corrupt, line).unwrap();
            write_spill_u64(&mut corrupt, fs).unwrap();
            write_spill_u64(&mut corrupt, fe).unwrap();
            write_spill_i32(&mut corrupt, 120).unwrap();
            write_spill_u8(&mut corrupt, 1).unwrap();
            write_spill_u8(&mut corrupt, 1).unwrap();
            write_spill_u8(&mut corrupt, 0).unwrap();
            assert!(SingleVcfRecord::from_parts_checked(
                line.to_string(),
                usize::try_from(fs).unwrap_or(usize::MAX),
                usize::try_from(fe).unwrap_or(usize::MAX),
                120,
                true,
                true,
                false,
            )
            .is_none());
            let err = RawVcfCandidate::read_from(&mut corrupt.as_slice())
                .err()
                .expect("corrupt SingleVcfRecord span/line must fail RawVcfCandidate::read_from");
            assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        }

        // Reject corrupt boolean bytes inside tag 3 (Record)
        for bad_bool_slot in 0..3 {
            let mut corrupt = Vec::new();
            write_spill_i32(&mut corrupt, 120).unwrap();
            write_spill_f64(&mut corrupt, 0.5).unwrap();
            write_spill_u8(&mut corrupt, 3).unwrap();
            write_spill_str(&mut corrupt, valid_line).unwrap();
            write_spill_u64(&mut corrupt, valid_fs).unwrap();
            write_spill_u64(&mut corrupt, valid_fe).unwrap();
            write_spill_i32(&mut corrupt, 120).unwrap();
            for slot in 0..3 {
                write_spill_u8(&mut corrupt, if slot == bad_bool_slot { 2 } else { 1 }).unwrap();
            }
            let err = RawVcfCandidate::read_from(&mut corrupt.as_slice())
                .err()
                .expect("corrupt boolean byte in Record must fail");
            assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        }
    }

    #[test]
    fn test_forced_disk_spill_single_and_somatic_cluster_e2e() {
        use crate::output::vcf::SingleVcfEvent;
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::{Cigar, Data};
        use noodles::sam::alignment::RecordBuf;

        // Part 1: Direct low-budget SpillableVcfAccumulator + emit_single_vcf_accumulator
        // verifying all 4 SingleVcfEvent variants across disk spills and chromosome resets.
        {
            let mut acc: SpillableVcfAccumulator<RawVcfCandidate> =
                SpillableVcfAccumulator::with_limits(64, 2);
            // chr1:30 PASS SNV + chr1:36 AmpBiasSnvClusterProbe (clusters 30 to Cluster10bp; suppressed under pass_only!)
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 30,
                    af_4dp: 0.5000,
                    event: SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr1", 30, "A", "G", "PASS", "DP=40", true, true, false,
                    )),
                },
            )
            .unwrap();
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 36,
                    af_4dp: 0.5000,
                    event: SingleVcfEvent::AmpBiasSnvClusterProbe { start: 36 },
                },
            )
            .unwrap();
            // chr1:60 PASS SNV -> chr1:62 SvSplitreadFlush (flushes 60 as PASS!) -> chr1:65 PASS SNV (gets Cluster10bp against pvs=60!)
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 60,
                    af_4dp: 0.5000,
                    event: SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr1", 60, "C", "T", "PASS", "DP=40", true, true, false,
                    )),
                },
            )
            .unwrap();
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 62,
                    af_4dp: 0.4000,
                    event: SingleVcfEvent::SvSplitreadFlush {
                        cluster_probe_start: None,
                    },
                },
            )
            .unwrap();
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 65,
                    af_4dp: 0.5000,
                    event: SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr1", 65, "G", "A", "PASS", "DP=40", true, true, false,
                    )),
                },
            )
            .unwrap();
            // chr1:90 Skip tombstone at higher AF shadows lower-AF PASS SNV from overlapping tile
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 90,
                    af_4dp: 0.3000,
                    event: SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr1", 90, "T", "C", "PASS", "DP=40", true, true, false,
                    )),
                },
            )
            .unwrap();
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 90,
                    af_4dp: 0.8000,
                    event: SingleVcfEvent::Skip,
                },
            )
            .unwrap();
            // chr1:200 PASS SNV and chr2:203 PASS SNV — chromosome reset prevents cross-chr clustering
            acc.insert(
                "chr1",
                RawVcfCandidate {
                    start_pos: 200,
                    af_4dp: 0.5000,
                    event: SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr1", 200, "A", "C", "PASS", "DP=40", true, true, false,
                    )),
                },
            )
            .unwrap();
            acc.insert(
                "chr2",
                RawVcfCandidate {
                    start_pos: 203,
                    af_4dp: 0.5000,
                    event: SingleVcfEvent::Record(mk_test_single_vcf_rec(
                        "chr2", 203, "G", "T", "PASS", "DP=40", true, true, false,
                    )),
                },
            )
            .unwrap();

            assert!(
                acc.spilled_runs.get("chr1").map_or(0, Vec::len) >= 3,
                "expected disk spills on chr1"
            );

            let cfg = crate::config::Config {
                vcf_cluster_dist: Some(crate::config::VcfThreshold::new(10)),
                vcf_pass_only: true,
                ..Default::default()
            };
            let mut out = Vec::new();
            emit_single_vcf_accumulator(acc, &cfg, &mut out).unwrap();
            let text = String::from_utf8(out).unwrap();
            let keys: Vec<&str> = text
                .lines()
                .map(|l| {
                    let mut cols = l.split('\t');
                    let chr = cols.next().unwrap();
                    let pos = cols.next().unwrap();
                    let _ = (cols.next(), cols.next(), cols.next(), cols.next());
                    let filt = cols.next().unwrap();
                    assert_eq!(filt, "PASS");
                    match (chr, pos) {
                        ("chr1", "60") => "chr1:60",
                        ("chr1", "200") => "chr1:200",
                        ("chr2", "203") => "chr2:203",
                        other => panic!("unexpected emitted record: {other:?} in {l}"),
                    }
                })
                .collect();
            assert_eq!(keys, vec!["chr1:60", "chr1:200", "chr2:203"]);
        }

        // Part 2: Live `run_with_writer` across multi-tile boundary (tile_size = 50) with
        // forced disk spill (`TEST_VCF_SPILL_BYTES_OVERRIDE = Some(1)`) in Simple and Somatic modes.
        let dir = tempfile::TempDir::new().unwrap();
        let ref_len = 200usize;
        let motif = b"ACGTACGTAGCTAGCTTGCAACGT";
        let mut ref_bytes = Vec::with_capacity(ref_len);
        for i in 0..ref_len {
            ref_bytes.push(motif[i % motif.len()]);
        }
        ref_bytes[47] = b'A'; // pos 48 (Tile 1: 1..=50)
        ref_bytes[53] = b'C'; // pos 54 (Tile 2: 51..=100; 54 - 48 = 6 < 10)
        ref_bytes[139] = b'G'; // pos 140 (Tile 3: 101..=150; unclustered PASS)
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        std::fs::write(
            &fasta_path,
            format!(">chr1\n{}\n", String::from_utf8(ref_bytes.clone()).unwrap()),
        )
        .unwrap();
        std::fs::write(&fai_path, format!("chr1\t{ref_len}\t6\t{ref_len}\t201\n")).unwrap();

        let mk_read = |name: &str, start_1: usize, seq: &[u8], rev: bool, nm: i32| -> RecordBuf {
            let flags = if rev {
                Flags::REVERSE_COMPLEMENTED
            } else {
                Flags::empty()
            };
            let cigar: Cigar = [Op::new(Kind::Match, seq.len())].into_iter().collect();
            let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(nm))]
                .into_iter()
                .collect();
            RecordBuf::builder()
                .set_name(name)
                .set_flags(flags)
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start_1).unwrap())
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::new(60).unwrap(),
                )
                .set_cigar(cigar)
                .set_sequence(seq.to_vec().into())
                .set_quality_scores(vec![38u8; seq.len()].into())
                .set_data(data)
                .build()
        };

        let mut tumor_reads = Vec::new();
        let mut normal_reads = Vec::new();
        // Spanning 35..=74 (covers 48 and 54 across tile boundary 50|51)
        for idx in 0..20usize {
            let start_1 = 30 + (idx % 6); // vary start/position so pstd=1
            let mut seq = ref_bytes[(start_1 - 1)..(start_1 - 1 + 40)].to_vec();
            let rev = idx % 2 == 1;
            normal_reads.push(mk_read(&format!("n1_{idx}"), start_1, &seq, rev, 0));
            if idx < 10 {
                seq[48 - start_1] = b'T';
                seq[54 - start_1] = b'T';
                tumor_reads.push(mk_read(&format!("t1_alt_{idx}"), start_1, &seq, rev, 2));
            } else {
                tumor_reads.push(mk_read(&format!("t1_ref_{idx}"), start_1, &seq, rev, 0));
            }
        }
        // Spanning 125..=165 (covers unclustered pos 140)
        for idx in 0..20usize {
            let start_1 = 120 + (idx % 6);
            let mut seq = ref_bytes[(start_1 - 1)..(start_1 - 1 + 40)].to_vec();
            let rev = idx % 2 == 1;
            normal_reads.push(mk_read(&format!("n2_{idx}"), start_1, &seq, rev, 0));
            if idx < 10 {
                seq[140 - start_1] = b'T';
                tumor_reads.push(mk_read(&format!("t2_alt_{idx}"), start_1, &seq, rev, 1));
            } else {
                tumor_reads.push(mk_read(&format!("t2_ref_{idx}"), start_1, &seq, rev, 0));
            }
        }
        tumor_reads.sort_by_key(|r| r.alignment_start().unwrap().get());
        normal_reads.sort_by_key(|r| r.alignment_start().unwrap().get());

        let tumor_bam = dir.path().join("tumor.bam");
        let normal_bam = dir.path().join("normal.bam");
        write_test_bam(&tumor_bam, "chr1", ref_len, &tumor_reads);
        write_test_bam(&normal_bam, "chr1", ref_len, &normal_reads);

        struct ResetSpillGuard;
        impl Drop for ResetSpillGuard {
            fn drop(&mut self) {
                TEST_VCF_SPILL_BYTES_OVERRIDE.with(|c| c.set(None));
            }
        }
        let _guard = ResetSpillGuard;
        // 1-byte limit forces EVERY candidate (including Skip, AmpBiasSnvClusterProbe,
        // and SvSplitreadFlush) to spill to disk runs and merge via BinaryHeap.
        TEST_VCF_SPILL_BYTES_OVERRIDE.with(|c| c.set(Some(1)));

        // Simple mode: multi-tile + forced disk spill tags 48 & 54 with Cluster010bp (raw token), keeps 140 PASS
        let simple_cfg = crate::config::Config {
            reference: fasta_path.to_string_lossy().into_owned(),
            bam: tumor_bam.to_string_lossy().into_owned(),
            region: Some("chr1:1-180".to_string()),
            tile_size: Some(50),
            vcf_cluster_dist: Some(crate::config::VcfThreshold::with_raw(10, "010".to_string())),
            threads: 1,
            ..Default::default()
        };
        let mut simple_buf = Vec::new();
        run_with_writer(&simple_cfg, &mut simple_buf).unwrap();
        let simple_vcf = String::from_utf8(simple_buf).unwrap();
        assert!(simple_vcf.contains(
            r#"##FILTER=<ID=Cluster010bp,Description="Two variants are within 010 bp">"#
        ));
        let simple_filts: Vec<(String, String)> = simple_vcf
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                (cols[1].to_string(), cols[6].to_string())
            })
            .collect();
        assert_eq!(
            simple_filts,
            vec![
                ("48".to_string(), "Cluster010bp".to_string()),
                ("54".to_string(), "Cluster010bp".to_string()),
                ("140".to_string(), "PASS".to_string()),
            ]
        );

        let mut simple_pass_cfg = simple_cfg.clone();
        simple_pass_cfg.vcf_pass_only = true;
        let mut simple_pass_buf = Vec::new();
        run_with_writer(&simple_pass_cfg, &mut simple_pass_buf).unwrap();
        let simple_pass_positions: Vec<String> = String::from_utf8(simple_pass_buf)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap().to_string())
            .collect();
        assert_eq!(simple_pass_positions, vec!["140".to_string()]);

        // Amplicon mode + forced disk spill: AmpBiasSnvClusterProbe at 54 suppresses 48 under --pass-only
        let mut amp_reads = Vec::new();
        // Amplicon 1 (20..=70, len 51): supports SNVs at 48 and 54
        for idx in 0..10usize {
            let mut seq = ref_bytes[19..70].to_vec();
            let rev = idx % 2 == 1;
            if idx < 5 {
                seq[48 - 20] = b'T';
                seq[54 - 20] = b'T';
                amp_reads.push(mk_read(&format!("amp1_alt_{idx}"), 20, &seq, rev, 2));
            } else {
                amp_reads.push(mk_read(&format!("amp1_ref_{idx}"), 20, &seq, rev, 0));
            }
        }
        // Amplicon 2 (35..=85, len 51): supports SNV ONLY at 48 -> 54 has GDAMP=1;TLAMP=2 (AMPBIAS)
        for idx in 0..10usize {
            let mut seq = ref_bytes[34..85].to_vec();
            let rev = idx % 2 == 1;
            if idx < 5 {
                seq[48 - 35] = b'T';
                amp_reads.push(mk_read(&format!("amp2_alt_{idx}"), 35, &seq, rev, 1));
            } else {
                amp_reads.push(mk_read(&format!("amp2_ref_{idx}"), 35, &seq, rev, 0));
            }
        }
        // Amplicon 3 (120..=170, len 51): unclustered PASS SNV at 140
        for idx in 0..10usize {
            let mut seq = ref_bytes[119..170].to_vec();
            let rev = idx % 2 == 1;
            if idx < 5 {
                seq[140 - 120] = b'T';
                amp_reads.push(mk_read(&format!("amp3_alt_{idx}"), 120, &seq, rev, 1));
            } else {
                amp_reads.push(mk_read(&format!("amp3_ref_{idx}"), 120, &seq, rev, 0));
            }
        }
        amp_reads.sort_by_key(|r| r.alignment_start().unwrap().get());
        let amp_bam = dir.path().join("amp_spill.bam");
        write_test_bam(&amp_bam, "chr1", ref_len, &amp_reads);
        let amp_bed = dir.path().join("amp_spill.bed");
        std::fs::write(
            &amp_bed,
            "chr1\t19\t70\tG1\t0\t+\t30\t65\nchr1\t34\t85\tG1\t0\t+\t40\t75\nchr1\t119\t170\tG2\t0\t+\t125\t160\n",
        )
        .unwrap();
        let amp_pass_cfg = crate::config::Config {
            reference: fasta_path.to_string_lossy().into_owned(),
            bam: amp_bam.to_string_lossy().into_owned(),
            bed: Some(amp_bed.to_string_lossy().into_owned()),
            amplicon: Some("10:0.95".to_string()),
            zero_based: true,
            zero_based_defined: true,
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 6,
            col_seg_end: 7,
            col_gene: Some(3),
            column_for_chromosome: Some(0),
            vcf_cluster_dist: Some(crate::config::VcfThreshold::new(10)),
            vcf_pass_only: true,
            threads: 1,
            ..Default::default()
        };
        let mut amp_pass_buf = Vec::new();
        run_with_writer(&amp_pass_cfg, &mut amp_pass_buf).unwrap();
        let amp_pass_positions: Vec<String> = String::from_utf8(amp_pass_buf)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap().to_string())
            .collect();
        assert_eq!(amp_pass_positions, vec!["140".to_string()]);

        // SV splitread flush + forced disk spill across tiles: <DEL> at 51 (splitreads=0 < 1)
        // flushes pending SNV 48 as PASS before clustered SNV 56 arrives.
        let mut sv_reads = Vec::new();
        for idx in 0..20usize {
            let start_1 = 30 + (idx % 6);
            let rev = idx % 2 == 1;
            let flags = if rev {
                Flags::REVERSE_COMPLEMENTED
            } else {
                Flags::empty()
            };
            if idx < 8 {
                let mut seq = ref_bytes[(start_1 - 1)..(start_1 - 1 + 40)].to_vec();
                seq[48 - start_1] = b'T';
                seq[56 - start_1] = b'A';
                sv_reads.push(mk_read(&format!("sv_snv_{idx}"), start_1, &seq, rev, 2));
            } else if idx < 13 {
                // 2bp deletion of 52..=53 -> anchor at pos 51, becomes <DEL> under sv_min_len=2
                let left_len = 51 - start_1 + 1;
                let right_len = 40 - left_len;
                let mut seq = Vec::with_capacity(40);
                seq.extend_from_slice(&ref_bytes[(start_1 - 1)..51]);
                seq.extend_from_slice(&ref_bytes[53..(53 + right_len)]);
                let cigar: Cigar = [
                    Op::new(Kind::Match, left_len),
                    Op::new(Kind::Deletion, 2),
                    Op::new(Kind::Match, right_len),
                ]
                .into_iter()
                .collect();
                let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(2i32))]
                    .into_iter()
                    .collect();
                sv_reads.push(
                    RecordBuf::builder()
                        .set_name(format!("sv_del_{idx}"))
                        .set_flags(flags)
                        .set_reference_sequence_id(0)
                        .set_alignment_start(Position::try_from(start_1).unwrap())
                        .set_mapping_quality(
                            noodles::sam::alignment::record::MappingQuality::new(60).unwrap(),
                        )
                        .set_cigar(cigar)
                        .set_sequence(seq.into())
                        .set_quality_scores(vec![38u8; 40].into())
                        .set_data(data)
                        .build(),
                );
            } else {
                let seq = ref_bytes[(start_1 - 1)..(start_1 - 1 + 40)].to_vec();
                sv_reads.push(mk_read(&format!("sv_ref_{idx}"), start_1, &seq, rev, 0));
            }
        }
        sv_reads.sort_by_key(|r| r.alignment_start().unwrap().get());
        let sv_bam = dir.path().join("sv_spill.bam");
        write_test_bam(&sv_bam, "chr1", ref_len, &sv_reads);
        let sv_flush_cfg = crate::config::Config {
            reference: fasta_path.to_string_lossy().into_owned(),
            bam: sv_bam.to_string_lossy().into_owned(),
            region: Some("chr1:1-100".to_string()),
            tile_size: Some(50),
            sv_min_len: 2,
            vcf_min_splitreads: 1,
            vcf_cluster_dist: Some(crate::config::VcfThreshold::new(10)),
            vcf_pass_only: true,
            threads: 1,
            ..Default::default()
        };
        let mut sv_flush_buf = Vec::new();
        run_with_writer(&sv_flush_cfg, &mut sv_flush_buf).unwrap();
        let sv_flush_positions: Vec<String> = String::from_utf8(sv_flush_buf)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap().to_string())
            .collect();
        assert_eq!(sv_flush_positions, vec!["48".to_string()]);

        // Paired somatic mode: multi-tile + forced disk spill tags StrongSomatic SNVs 48 & 54 with Cluster010bp
        let somatic_cfg = crate::config::Config {
            reference: fasta_path.to_string_lossy().into_owned(),
            bam: format!(
                "{}|{}",
                tumor_bam.to_string_lossy(),
                normal_bam.to_string_lossy()
            ),
            region: Some("chr1:1-180".to_string()),
            tile_size: Some(50),
            vcf_cluster_dist: Some(crate::config::VcfThreshold::with_raw(10, "010".to_string())),
            threads: 1,
            ..Default::default()
        };
        let mut som_buf = Vec::new();
        run_with_writer(&somatic_cfg, &mut som_buf).unwrap();
        let som_vcf = String::from_utf8(som_buf).unwrap();
        assert!(som_vcf.contains(
            r#"##FILTER=<ID=Cluster010bp,Description="Two somatic variants are within 010 bp">"#
        ));
        let som_filts: Vec<(String, String)> = som_vcf
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                (cols[1].to_string(), cols[6].to_string())
            })
            .collect();
        assert_eq!(
            som_filts,
            vec![
                ("48".to_string(), "Cluster010bp".to_string()),
                ("54".to_string(), "Cluster010bp".to_string()),
                ("140".to_string(), "PASS".to_string()),
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
            rustc_hash::FxHashSet::default(),
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
        let (eff_splice, mode_splice, pre_splice) = resolve_effective_config(&splice_cfg).unwrap();
        assert_eq!(mode_splice, CallingMode::Splicing);
        assert_eq!(eff_splice.sample, "TUMOR");
        assert_eq!(eff_splice.bam, "tumor.bam");
        assert!(pre_splice.is_none());

        // 2. Explicit -a on a 9-column BED preserves amplicon parameters without pre-reading BED
        let amp9_cfg = crate::config::Config {
            bam: "a.bam".to_string(),
            bed: Some(amp9_bed.to_string_lossy().into_owned()),
            amplicon: Some("15:0.90".to_string()),
            ..Default::default()
        };
        let (eff_amp9, mode_amp9, pre_amp9) = resolve_effective_config(&amp9_cfg).unwrap();
        assert_eq!(mode_amp9, CallingMode::Amplicon);
        assert_eq!(eff_amp9.amplicon.as_deref(), Some("15:0.90"));
        assert!(pre_amp9.is_none());

        // 3. 8-column amplicon BED auto-detects even when -c 1 (column_for_chromosome) is set,
        //    preserves explicit -z 0 (`zero_based_defined = true, zero_based = false`), streams
        //    regular seekable files without `Vec<String>` slurp (`WF-4`: `pre_amp8_z0.is_none()`),
        //    and single-pass buffers non-regular pipe/FD streams (`Some(len == 1)`).
        let amp8_z0_cfg = crate::config::Config {
            bam: "a.bam".to_string(),
            bed: Some(amp8_bed.to_string_lossy().into_owned()),
            column_for_chromosome: Some(0),
            zero_based: false,
            zero_based_defined: true,
            ..Default::default()
        };
        let (eff_amp8_z0, mode_amp8_z0, pre_amp8_z0) =
            resolve_effective_config(&amp8_z0_cfg).unwrap();
        assert_eq!(mode_amp8_z0, CallingMode::Amplicon);
        assert_eq!(eff_amp8_z0.amplicon.as_deref(), Some("10:0.95"));
        assert!(
            pre_amp8_z0.is_none(),
            "regular BED file must stream in O(1) line memory without preloading Vec<String>"
        );
        assert!(
            !eff_amp8_z0.zero_based,
            "explicit -z 0 must not be overwritten by amplicon auto-detect"
        );
        #[cfg(target_os = "linux")]
        {
            use std::io::Write as _;
            use std::os::fd::AsRawFd as _;
            let (pipe_r, mut pipe_w) = std::io::pipe().unwrap();
            pipe_w
                .write_all(&std::fs::read(&amp8_bed).unwrap())
                .unwrap();
            drop(pipe_w);
            let fifo_cfg = crate::config::Config {
                bam: "a.bam".to_string(),
                bed: Some(format!("/proc/self/fd/{}", pipe_r.as_raw_fd())),
                ..Default::default()
            };
            let (eff_fifo, mode_fifo, pre_fifo) = resolve_effective_config(&fifo_cfg).unwrap();
            assert_eq!(mode_fifo, CallingMode::Amplicon);
            assert_eq!(eff_fifo.amplicon.as_deref(), Some("10:0.95"));
            assert_eq!(
                pre_fifo.as_ref().map(Vec::len),
                Some(1),
                "non-regular pipe/FD BED stream must be preloaded in a single pass"
            );
            drop(pipe_r);
        }

        assert!(is_os_fd_alias_path(Path::new("/dev/stdout")));
        assert!(is_os_fd_alias_path(Path::new("/dev/stderr")));
        assert!(is_os_fd_alias_path(Path::new("/dev/fd/3")));
        assert!(is_os_fd_alias_path(Path::new("/proc/self/fd/2")));
        assert!(!is_os_fd_alias_path(Path::new("/dev/null")));
        #[cfg(unix)]
        {
            let sym_stdout = dir.path().join("sym_stdout.vcf");
            std::os::unix::fs::symlink("/dev/stdout", &sym_stdout).unwrap();
            assert!(targets_os_fd_alias(&sym_stdout));
            assert_eq!(
                classify_target_os_fd_alias(&sym_stdout),
                Some(OsFdAliasKind::SelfStdout)
            );
        }

        check_resolve_effective_config_vcf_flag_guards();
    }

    fn check_resolve_effective_config_vcf_flag_guards() {
        for bad_pval in [f64::NAN, -0.01, f64::INFINITY, f64::NEG_INFINITY] {
            let cfg = crate::config::Config {
                bam: "a.bam".to_string(),
                region: Some("chr1:1-10".to_string()),
                vcf_pval: bad_pval,
                ..Default::default()
            };
            let err = resolve_effective_config(&cfg)
                .expect_err(&format!("expected Parse error for vcf_pval={bad_pval}"));
            assert!(matches!(err, VardictError::Parse(_)));
        }

        let neg_split_cfg = crate::config::Config {
            bam: "a.bam".to_string(),
            region: Some("chr1:1-10".to_string()),
            vcf_min_splitreads: -1,
            ..Default::default()
        };
        let split_err = resolve_effective_config(&neg_split_cfg)
            .expect_err("expected Parse error for vcf_min_splitreads = -1");
        assert!(matches!(split_err, VardictError::Parse(_)));

        let neg_zero_cfg = crate::config::Config {
            bam: "a.bam".to_string(),
            region: Some("chr1:1-10".to_string()),
            vcf_pval: -0.0,
            vcf_pval_raw: Some("-0.0".to_string()),
            vcf_min_splitreads: 0,
            ..Default::default()
        };
        let (eff_zero, _, _) = resolve_effective_config(&neg_zero_cfg).unwrap();
        assert_eq!(eff_zero.vcf_pval, 0.0);
        assert!(eff_zero.vcf_pval.is_sign_positive());
        assert_eq!(eff_zero.vcf_pval_raw.as_deref(), Some("0.0"));
        assert_eq!(eff_zero.vcf_min_splitreads, 0);

        // Programmatic `vcf_pval_raw` syntax, header-injection, underflow, and bit-agreement guards
        for (pval, bad_raw) in [
            (0.05, "0.01"),        // disagrees with vcf_pval = 0.05
            (0.05, "0.05\n##bad"), // header injection / invalid float
            (0.05, " 0.05"),       // untrimmed leading whitespace
            (0.05, "0.05\n"),      // untrimmed trailing newline
            (0.0, ""),             // empty
            (0.0, "-"),            // lone sign
            (0.0, "-1e-400"),      // negative non-zero mantissa underflowing to -0.0
            (0.05, "NaN"),         // non-finite raw
            (0.05, "inf"),         // infinite raw
        ] {
            let bad_raw_cfg = crate::config::Config {
                bam: "a.bam".to_string(),
                region: Some("chr1:1-10".to_string()),
                vcf_pval: pval,
                vcf_pval_raw: Some(bad_raw.to_string()),
                ..Default::default()
            };
            let raw_err = resolve_effective_config(&bad_raw_cfg).expect_err(&format!(
                "expected Parse error for (vcf_pval={pval}, vcf_pval_raw={bad_raw:?})"
            ));
            assert!(matches!(raw_err, VardictError::Parse(_)));
        }

        // Valid matching `vcf_pval = 1e-5` + `vcf_pval_raw = Some("1e-5")` succeeds and preserves raw
        let sci_cfg = crate::config::Config {
            bam: "a.bam".to_string(),
            region: Some("chr1:1-10".to_string()),
            vcf_pval: 1e-5,
            vcf_pval_raw: Some("1e-5".to_string()),
            ..Default::default()
        };
        let (eff_sci, _, _) = resolve_effective_config(&sci_cfg).unwrap();
        assert_eq!(eff_sci.vcf_pval.to_bits(), 1e-5_f64.to_bits());
        assert_eq!(eff_sci.vcf_pval_raw.as_deref(), Some("1e-5"));

        // Programmatic `validate_vcf_thresholds` guard for all new `vcf_*` override fields
        for mut bad_vcf_cfg in [
            crate::config::Config {
                vcf_af: Some(crate::config::VcfThreshold::new(-0.01)),
                ..sci_cfg.clone()
            },
            crate::config::Config {
                vcf_min_depth: Some(crate::config::VcfThreshold::new(-1)),
                ..sci_cfg.clone()
            },
            crate::config::Config {
                vcf_cluster_dist: Some(crate::config::VcfThreshold::new(-2)),
                ..sci_cfg.clone()
            },
            crate::config::Config {
                vcf_min_vd: Some(crate::config::VcfThreshold::with_raw(2, "99".to_string())),
                ..sci_cfg.clone()
            },
            crate::config::Config {
                vcf_gt_freq: Some(1.01),
                ..sci_cfg.clone()
            },
        ] {
            bad_vcf_cfg.vcf_pval_raw = None;
            assert!(matches!(
                resolve_effective_config(&bad_vcf_cfg).unwrap_err(),
                VardictError::Parse(_)
            ));
        }
    }

    #[test]
    fn test_split_var_deduper_ownership_multiplicity_and_retention_window() {
        let mk_row = |pos: i32, r: &str, v: &str, line: &str| SplitVarRow {
            id: SplitVarIdentity {
                origin_pos: pos,
                start_pos: pos,
                end_pos: pos,
                ref_allele: r.to_string(),
                var_allele: v.to_string(),
                null_record: false,
                companion_kind: 0,
            },
            is_sv: false,
            line: line.to_string(),
        };

        let mut deduper = SplitVarDeduper::default();
        let mut emitted = Vec::new();

        // 1) Cross-tile ownership resolution when tile-local depth/AF/status/debug differ:
        //    Parent region [1, 200_000] split into Tile 1 [1, 100_000] and Tile 2 [100_001, 200_000].
        //    - Pos 99_995 A>G is owned by Tile 1; Tile 2 also produces a lookbehind copy with
        //      different depth/AF ("dp=12" vs "dp=50") -> Tile 1's owned copy wins.
        //    - Pos 99_995 A>T (distinct allele at same position) in Tile 1 is also preserved.
        //    - Pos 100_005 C>T is owned by Tile 2; Tile 1 produces a lookahead copy ("dp=8")
        //      and Tile 2 produces its owned copy ("dp=60") -> Tile 2's owned copy wins.
        //    - Pos 100_010 G>A is produced ONLY by Tile 1 (e.g. shifted variant) -> flushed
        //      when Tile 2 arrives without it, so it is never lost.
        deduper
            .accept_tile(
                vec![
                    mk_row(99_995, "A", "G", "pos99995_A_G_owner_tile1_dp50"),
                    mk_row(99_995, "A", "T", "pos99995_A_T_owner_tile1_dp40"),
                    mk_row(100_005, "C", "T", "pos100005_C_T_nonowner_tile1_dp8"),
                    mk_row(100_010, "G", "A", "pos100010_G_A_only_tile1_dp30"),
                ],
                true,
                p(1),
                p(100_000),
                p(1),
                p(200_000),
                |l| {
                    emitted.push(l.to_string());
                    Ok(())
                },
            )
            .unwrap();

        deduper
            .accept_tile(
                vec![
                    mk_row(99_995, "A", "G", "pos99995_A_G_nonowner_tile2_dp12"),
                    mk_row(100_005, "C", "T", "pos100005_C_T_owner_tile2_dp60"),
                ],
                true,
                p(100_001),
                p(200_000),
                p(1),
                p(200_000),
                |l| {
                    emitted.push(l.to_string());
                    Ok(())
                },
            )
            .unwrap();

        assert_eq!(
            emitted,
            vec![
                "pos99995_A_G_owner_tile1_dp50",
                "pos99995_A_T_owner_tile1_dp40",
                "pos100010_G_A_only_tile1_dp30",
                "pos100005_C_T_owner_tile2_dp60",
            ]
        );

        // 2) Within-tile multiplicity and `-p` shifted-singleton position-0 companion rows:
        //    Two shifted singletons without `ref_variant` in the same tile (at pos 100 and pos 200)
        //    each emit [null_record(pos), Variant::new()(0..0), shifted_vref]. Both pos-0 companion
        //    rows must be preserved within the tile, and a third shifted singleton at pos 150_000
        //    in Tile 2 must also preserve its pos-0 companion row, while a cross-tile repeat of
        //    the pos 100 shifted singleton in Tile 2 is suppressed as a unit.
        let mk_shifted_triple = |origin: i32, shifted: i32, r: &str, v: &str| {
            let mut n = crate::variations::Variant::new();
            n.start_pos = origin;
            n.end_pos = origin;
            n.null_record = true;
            let z = crate::variations::Variant::new();
            let mut s = crate::variations::Variant::new();
            s.start_pos = shifted;
            s.end_pos = shifted;
            s.ref_allele = r.to_string();
            s.var_allele = v.to_string();
            vec![n, z, s]
        };

        let mut tile1_vars = mk_shifted_triple(100, 102, "A", "C");
        tile1_vars.extend(mk_shifted_triple(200, 203, "G", "T"));
        let tile1_ids = build_simple_split_var_identities(&tile1_vars, true);
        assert_eq!(tile1_ids.len(), 6);
        assert_ne!(
            tile1_ids[1], tile1_ids[4],
            "distinct shifted singletons must have distinct pos-0 companion identities"
        );
        let tile1_rows: Vec<SplitVarRow> = tile1_ids
            .into_iter()
            .enumerate()
            .map(|(i, id)| SplitVarRow {
                id,
                is_sv: false,
                line: format!("t1_row_{i}"),
            })
            .collect();

        // Tile 2 repeats the (100 -> 102 A>C) triple from Tile 1, plus a new (150_000 -> 150_002 C>A) triple.
        let mut tile2_vars = mk_shifted_triple(100, 102, "A", "C");
        tile2_vars.extend(mk_shifted_triple(150_000, 150_002, "C", "A"));
        let tile2_ids = build_simple_split_var_identities(&tile2_vars, true);
        let tile2_rows: Vec<SplitVarRow> = tile2_ids
            .into_iter()
            .enumerate()
            .map(|(i, id)| SplitVarRow {
                id,
                is_sv: false,
                line: format!("t2_row_{i}"),
            })
            .collect();

        emitted.clear();
        deduper
            .accept_tile(tile1_rows, true, p(1), p(100_000), p(1), p(200_000), |l| {
                emitted.push(l.to_string());
                Ok(())
            })
            .unwrap();
        deduper
            .accept_tile(
                tile2_rows,
                true,
                p(100_001),
                p(200_000),
                p(1),
                p(200_000),
                |l| {
                    emitted.push(l.to_string());
                    Ok(())
                },
            )
            .unwrap();
        assert_eq!(
            emitted,
            vec![
                "t1_row_0", "t1_row_1", "t1_row_2", "t1_row_3", "t1_row_4", "t1_row_5", "t2_row_3",
                "t2_row_4", "t2_row_5",
            ]
        );

        // 3) Sliding boundary-retention window stays active past 262,144 rows while evicting
        //    obsolete positions < tile_start - SPLIT_VAR_RETENTION_WINDOW_BP.
        let mut big_deduper = SplitVarDeduper::default();
        let tile1_big: Vec<SplitVarRow> = (1..=100_000)
            .flat_map(|pos| {
                [
                    mk_row(pos, "A", "C", "ac"),
                    mk_row(pos, "A", "G", "ag"),
                    mk_row(pos, "A", "T", "at"),
                ]
            })
            .collect();
        assert_eq!(tile1_big.len(), 300_000); // > 262,144 rows in Tile 1 alone
        let mut count = 0usize;
        big_deduper
            .accept_tile(tile1_big, true, p(1), p(100_000), p(1), p(500_000), |_| {
                count += 1;
                Ok(())
            })
            .unwrap();
        assert_eq!(count, 300_000);

        // Tile 2 [100_001, 200_000] repeats pos 100_000 (the 300,000th row, well past 262,144):
        // deduplication must still suppress it!
        let mut tile2_emitted = Vec::new();
        big_deduper
            .accept_tile(
                vec![
                    mk_row(100_000, "A", "T", "dup_past_262k"),
                    mk_row(100_001, "C", "G", "new_in_tile2"),
                ],
                true,
                p(100_001),
                p(200_000),
                p(1),
                p(500_000),
                |l| {
                    tile2_emitted.push(l.to_string());
                    Ok(())
                },
            )
            .unwrap();
        assert_eq!(tile2_emitted, vec!["new_in_tile2"]);

        // Tile 5 [400_001, 500_000] advances tile_start to 400_001, so cutoff = 200_001
        // evicts all Tile 1 positions (1..=100_000) from `emitted_by_pos`.
        big_deduper
            .accept_tile(
                vec![mk_row(400_001, "G", "A", "tile5_row")],
                true,
                p(400_001),
                p(500_000),
                p(1),
                p(500_000),
                |_| Ok(()),
            )
            .unwrap();
        assert!(
            !big_deduper.emitted_by_pos.contains_key(&100_000),
            "positions below tile_start - SPLIT_VAR_RETENTION_WINDOW_BP must be evicted"
        );

        // 4) flush_before_unsplit_tile flushes pending forward rows from a preceding split tile
        //    and clears retention sets before raw UnSplitVarBlock emission.
        let mut flush_deduper = SplitVarDeduper::default();
        let mut flushed_lines = Vec::new();
        flush_deduper
            .accept_tile(
                vec![mk_row(100_005, "C", "T", "pending_before_unsplit")],
                true,
                p(1),
                p(100_000),
                p(1),
                p(200_000),
                |l| {
                    flushed_lines.push(l.to_string());
                    Ok(())
                },
            )
            .unwrap();
        assert!(flushed_lines.is_empty());
        flush_deduper
            .flush_before_unsplit_tile(|l| {
                flushed_lines.push(l.to_string());
                Ok(())
            })
            .unwrap();
        assert_eq!(flushed_lines, vec!["pending_before_unsplit"]);
        assert!(flush_deduper.emitted_by_pos.is_empty());
        assert!(flush_deduper.emitted_sv.is_empty());
    }

    #[test]
    fn test_write_vcf_tabix_index_honors_info_end_for_symbolic_sv() {
        use noodles::csi::BinningIndex as _;
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();
        let vcf_gz = dir.path().join("sv.vcf.gz");
        let tbi_path = dir.path().join("sv.vcf.gz.tbi");
        {
            let file = std::fs::File::create(&vcf_gz).unwrap();
            let mut writer = noodles::bgzf::io::Writer::new(file);
            writer
                .write_all(
                    b"##fileformat=VCFv4.3\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                      chr1\t1000\t.\tA\t<DEL>\t100\tPASS\tSVTYPE=DEL;END=50_000_invalid;END=50000\tGT\t0/1\n",
                )
                .unwrap();
            writer.finish().unwrap();
        }
        stage_vcf_tabix_index(&vcf_gz, &vcf_gz)
            .unwrap()
            .commit()
            .unwrap();
        let index = noodles::tabix::fs::read(&tbi_path).unwrap();
        let start = noodles::core::Position::new(30000).unwrap();
        let end = noodles::core::Position::new(30100).unwrap();
        let chunks = index.query(0, (start..=end).into()).unwrap();
        assert!(
            !chunks.is_empty(),
            "symbolic SV at POS=1000 with INFO/END=50000 must overlap query [30000, 30100]"
        );
    }

    #[test]
    fn test_paths_collide_input_guards_and_dev_null() {
        let dir = tempfile::TempDir::new().unwrap();
        let ref_fa = dir.path().join("ref.fa");
        let ref_fai = dir.path().join("ref.fa.fai");
        std::fs::write(
            &ref_fa,
            b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        )
        .unwrap();
        std::fs::write(&ref_fai, b"chr1\t60\t6\t60\t61\n").unwrap();

        let tumor_bam = dir.path().join("tumor.bam");
        let normal_bam = dir.path().join("normal.bam");
        let rec = test_record_buf(
            "r1",
            10,
            &[(10, noodles::sam::alignment::record::cigar::op::Kind::Match)],
            b"ACGTACGTAC".to_vec(),
        );
        write_test_bam(&tumor_bam, "chr1", 60, std::slice::from_ref(&rec));
        write_test_bam(&normal_bam, "chr1", 60, std::slice::from_ref(&rec));
        let tumor_bai = dir.path().join("tumor.bam.bai");
        let normal_bai = dir.path().join("normal.bam.bai");
        let tumor_csi = dir.path().join("tumor.bam.csi");
        std::fs::write(&tumor_csi, b"CSI_SENTINEL").unwrap();

        let extra_cram = dir.path().join("extra.cram");
        let extra_crai = dir.path().join("extra.cram.crai");
        std::fs::write(&extra_cram, b"CRAM_SENTINEL").unwrap();
        std::fs::write(&extra_crai, b"CRAI_SENTINEL").unwrap();

        let bed_file = dir.path().join("targets.bed");
        std::fs::write(&bed_file, b"chr1\t1\t50\tGENE1\n").unwrap();

        let protected = [
            &ref_fa,
            &ref_fai,
            &bed_file,
            &tumor_bam,
            &normal_bam,
            &tumor_bai,
            &normal_bai,
            &tumor_csi,
            &extra_cram,
            &extra_crai,
        ];
        let snapshots: Vec<Vec<u8>> = protected
            .iter()
            .map(|p| std::fs::read(p).unwrap())
            .collect();

        let base_cfg = crate::config::Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: format!(
                "{}:{}|{}",
                tumor_bam.display(),
                extra_cram.display(),
                normal_bam.display()
            ),
            bed: Some(bed_file.to_string_lossy().into_owned()),
            column_for_chromosome: Some(0),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            splice: true, // exercises raw config.bam retaining normal.bam even under -i
            threads: 1,
            ..Default::default()
        };

        for target in protected {
            // 1. --output colliding directly with protected input/index
            let mut out_cfg = base_cfg.clone();
            out_cfg.output_path = Some(target.to_string_lossy().into_owned());
            let err = crate::run(&out_cfg).expect_err(&format!(
                "--output colliding with {} must fail",
                target.display()
            ));
            assert!(
                matches!(err, VardictError::Parse(_)),
                "expected Parse error for --output {}, got {err:?}",
                target.display()
            );

            // 2. --stats-json colliding directly with protected input/index
            let mut stats_cfg = base_cfg.clone();
            stats_cfg.stats_json_path = Some(target.to_string_lossy().into_owned());
            let err = crate::run(&stats_cfg).expect_err(&format!(
                "--stats-json colliding with {} must fail",
                target.display()
            ));
            assert!(
                matches!(err, VardictError::Parse(_)),
                "expected Parse error for --stats-json {}, got {err:?}",
                target.display()
            );

            // 3. --tbi output colliding via symlink with protected input/index
            #[cfg(unix)]
            {
                let safe_vcf_gz = dir.path().join("tbi_guard_check.vcf.gz");
                let safe_tbi = dir.path().join("tbi_guard_check.vcf.gz.tbi");
                let _ = std::fs::remove_file(&safe_tbi);
                std::os::unix::fs::symlink(target, &safe_tbi).unwrap();
                let mut tbi_cfg = base_cfg.clone();
                tbi_cfg.splice = false;
                tbi_cfg.output_path = Some(safe_vcf_gz.to_string_lossy().into_owned());
                tbi_cfg.index_output = true;
                let err = crate::run(&tbi_cfg).expect_err(&format!(
                    "--tbi colliding with {} must fail",
                    target.display()
                ));
                assert!(
                    matches!(err, VardictError::Parse(_)),
                    "expected Parse error for --tbi {}, got {err:?}",
                    target.display()
                );
                assert!(
                    !safe_vcf_gz.exists(),
                    "--tbi collision check must reject before creating {}",
                    safe_vcf_gz.display()
                );
                let _ = std::fs::remove_file(&safe_tbi);
            }
        }

        // Verify zero bytes were truncated or mutated on any input/companion file
        for (p, expected) in protected.iter().zip(snapshots.iter()) {
            assert_eq!(
                &std::fs::read(p).unwrap(),
                expected,
                "protected file {} was truncated or modified!",
                p.display()
            );
        }

        // Verify --tbi error message no longer mislabels -E as --emit-intermediate
        let mut bad_tbi_inter = base_cfg.clone();
        bad_tbi_inter.splice = false;
        bad_tbi_inter.emit_intermediate = true;
        bad_tbi_inter.index_output = true;
        bad_tbi_inter.output_path =
            Some(dir.path().join("ok.vcf.gz").to_string_lossy().into_owned());
        let msg = crate::run(&bad_tbi_inter).unwrap_err().to_string();
        assert!(msg.contains("--emit-intermediate"));
        assert!(!msg.contains("-E/--emit-intermediate"));

        // Verify --output /dev/null --stats-json /dev/null succeeds on Unix
        #[cfg(unix)]
        {
            let dev_null_cfg = crate::config::Config {
                reference: ref_fa.to_string_lossy().into_owned(),
                bam: tumor_bam.to_string_lossy().into_owned(),
                region: Some("chr1:1-50".to_string()),
                output_path: Some("/dev/null".to_string()),
                stats_json_path: Some("/dev/null".to_string()),
                threads: 1,
                ..Default::default()
            };
            crate::run(&dev_null_cfg)
                .expect("--output /dev/null --stats-json /dev/null must succeed on Unix");
        }
    }

    #[test]
    fn test_paths_collide_bgzf_residue_and_atomic_tmp_guards() {
        let dir = tempfile::TempDir::new().unwrap();
        let ref_fa = dir.path().join("ref.fa");
        let ref_fai = dir.path().join("ref.fa.fai");
        std::fs::write(
            &ref_fa,
            b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        )
        .unwrap();
        std::fs::write(&ref_fai, b"chr1\t60\t6\t60\t61\n").unwrap();

        let bam = dir.path().join("sample.bam");
        let rec = test_record_buf(
            "r1",
            10,
            &[(10, noodles::sam::alignment::record::cigar::op::Kind::Match)],
            b"ACGTACGTAC".to_vec(),
        );
        write_test_bam(&bam, "chr1", 60, std::slice::from_ref(&rec));

        // 1. Malformed region (-R chr1:bad) fails before writing any bytes:
        //    - newly created BGZF --output must be unlinked (no 28-byte BGZF_EOF residue)
        //    - pre-existing --stats-json must remain intact with no .tmp.* left behind
        let residue_vcf_gz = dir.path().join("residue.vcf.gz");
        let existing_stats = dir.path().join("existing_stats.json");
        std::fs::write(&existing_stats, b"{\"preserved\":true}\n").unwrap();

        let bad_cfg = crate::config::Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: bam.to_string_lossy().into_owned(),
            region: Some("chr1:bad".to_string()),
            output_path: Some(residue_vcf_gz.to_string_lossy().into_owned()),
            stats_json_path: Some(existing_stats.to_string_lossy().into_owned()),
            threads: 1,
            ..Default::default()
        };
        assert!(crate::run(&bad_cfg).is_err());
        assert!(
            !residue_vcf_gz.exists(),
            "failed run with 0 output bytes must unlink newly created BGZF residue"
        );
        assert_eq!(
            std::fs::read(&existing_stats).unwrap(),
            b"{\"preserved\":true}\n",
            "failed run must not truncate pre-existing --stats-json"
        );

        // 2. Pre-existing --output file is preserved byte-for-byte (neither unlinked nor truncated
        //    into a 28-byte BGZF EOF block) on early error, and dangling symlinks leave zero residue.
        let preexisting_out = dir.path().join("preexisting.vcf.gz");
        std::fs::write(&preexisting_out, b"old").unwrap();
        let mut pre_cfg = bad_cfg.clone();
        pre_cfg.output_path = Some(preexisting_out.to_string_lossy().into_owned());
        pre_cfg.stats_json_path = None;
        assert!(crate::run(&pre_cfg).is_err());
        assert_eq!(
            std::fs::read(&preexisting_out).unwrap(),
            b"old",
            "pre-existing --output file must remain byte-for-byte intact on error"
        );

        #[cfg(unix)]
        {
            let dangling_target = dir.path().join("dangling_target.vcf.gz");
            let dangling_link = dir.path().join("dangling_link.vcf.gz");
            std::os::unix::fs::symlink(&dangling_target, &dangling_link).unwrap();
            let mut sym_cfg = bad_cfg.clone();
            sym_cfg.output_path = Some(dangling_link.to_string_lossy().into_owned());
            sym_cfg.stats_json_path = None;
            assert!(crate::run(&sym_cfg).is_err());
            assert!(
                !dangling_target.exists(),
                "failed run on dangling symlink must not create target residue"
            );
            assert!(
                dangling_link.symlink_metadata().is_ok(),
                "dangling symlink itself must not be removed on failed run"
            );
        }

        // Public OutputDestination implements Debug (M-PUBLIC-DEBUG), and VardictError::is_broken_pipe works.
        let dbg_str = format!(
            "{:?}",
            create_output_writer(&crate::config::Config::default())
        );
        assert!(
            dbg_str.starts_with("Stdout"),
            "unexpected Debug output: {dbg_str}"
        );
        assert!(
            VardictError::Io(std::io::Error::from(std::io::ErrorKind::BrokenPipe)).is_broken_pipe()
        );
        assert!(!VardictError::Parse("bad".to_string()).is_broken_pipe());
        let jm: crate::javamap::JavaHashMap<i32, i32> = Default::default();
        assert!(!format!("{:?}", jm.iter()).is_empty());
        let cached_bgzf = crate::bam::CachedBgzfReader::new(std::fs::File::open(&bam).unwrap(), 42);
        assert!(format!("{cached_bgzf:?}").contains("CachedBgzfReader"));

        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt as _;
            let ro_dir = dir.path().join("readonly_bind_dir");
            std::fs::create_dir(&ro_dir).unwrap();
            let ro_out = ro_dir.join("mounted.vcf.gz");
            let ro_tbi = ro_dir.join("mounted.vcf.gz.tbi");
            let ro_stats = ro_dir.join("mounted_stats.json");
            std::fs::write(&ro_out, b"seed").unwrap();
            std::fs::write(&ro_stats, b"seed_stats").unwrap();
            std::fs::set_permissions(&ro_dir, std::fs::Permissions::from_mode(0o555)).unwrap();

            let ro_enforced = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(ro_dir.join("probe_nonroot"))
                .is_err();
            if ro_enforced {
                // Missing pre-created .tbi inside 0555 dir must fail preflight without truncating ro_out.
                let missing_tbi_cfg = crate::config::Config {
                    reference: ref_fa.to_string_lossy().into_owned(),
                    bam: bam.to_string_lossy().into_owned(),
                    region: Some("chr1:1-50".to_string()),
                    output_path: Some(ro_out.to_string_lossy().into_owned()),
                    index_output: true,
                    stats_json_path: Some(ro_stats.to_string_lossy().into_owned()),
                    threads: 1,
                    ..Default::default()
                };
                assert!(crate::run(&missing_tbi_cfg).is_err());
                assert_eq!(
                    std::fs::read(&ro_out).unwrap(),
                    b"seed",
                    "missing companion .tbi in read-only parent must not truncate existing output"
                );
                assert_eq!(std::fs::read(&ro_stats).unwrap(), b"seed_stats");
            }

            std::fs::set_permissions(&ro_dir, std::fs::Permissions::from_mode(0o755)).unwrap();
            std::fs::write(&ro_tbi, b"seed_tbi").unwrap();
            std::fs::set_permissions(&ro_dir, std::fs::Permissions::from_mode(0o555)).unwrap();

            if ro_enforced {
                // Mid-run / region-parse failure with pre-created .tbi must also leave ro_out untouched.
                let bad_region_cfg = crate::config::Config {
                    reference: ref_fa.to_string_lossy().into_owned(),
                    bam: bam.to_string_lossy().into_owned(),
                    region: Some("chr1:bad".to_string()),
                    output_path: Some(ro_out.to_string_lossy().into_owned()),
                    index_output: true,
                    stats_json_path: Some(ro_stats.to_string_lossy().into_owned()),
                    threads: 1,
                    ..Default::default()
                };
                assert!(crate::run(&bad_region_cfg).is_err());
                assert_eq!(
                    std::fs::read(&ro_out).unwrap(),
                    b"seed",
                    "failed run in read-only parent must not truncate existing output"
                );
                assert_eq!(std::fs::read(&ro_tbi).unwrap(), b"seed_tbi");
                assert_eq!(std::fs::read(&ro_stats).unwrap(), b"seed_stats");
            }

            let ro_cfg = crate::config::Config {
                reference: ref_fa.to_string_lossy().into_owned(),
                bam: bam.to_string_lossy().into_owned(),
                region: Some("chr1:1-50".to_string()),
                output_path: Some(ro_out.to_string_lossy().into_owned()),
                index_output: true,
                stats_json_path: Some(ro_stats.to_string_lossy().into_owned()),
                threads: 1,
                ..Default::default()
            };
            let ro_res = crate::run(&ro_cfg);
            let stats_after = std::fs::read_to_string(&ro_stats).unwrap_or_default();
            let tbi_len = std::fs::metadata(&ro_tbi).map(|m| m.len()).unwrap_or(0);
            std::fs::set_permissions(&ro_dir, std::fs::Permissions::from_mode(0o755)).unwrap();
            ro_res.expect(
                "writable bind-mounted files inside read-only parent directory must succeed via direct fallback",
            );
            assert!(stats_after.contains("\"mode\": \"Simple\""));
            assert!(
                tbi_len > 8,
                "expected valid Tabix index in ro_tbi, len={tbi_len}"
            );
        }

        // 3. stage_vcf_tabix_index uses atomic tmp path and never unlinks existing .tbi on error
        let corrupt_vcf_gz = dir.path().join("corrupt.vcf.gz");
        let existing_tbi = dir.path().join("corrupt.vcf.gz.tbi");
        {
            use std::io::Write as _;
            let mut bgzf_w =
                noodles::bgzf::io::Writer::new(std::fs::File::create(&corrupt_vcf_gz).unwrap());
            bgzf_w.write_all(b"chr1\tNOT_A_POS\t.\tA\tG\n").unwrap();
            bgzf_w.finish().unwrap();
        }
        std::fs::write(&existing_tbi, b"KEEP_EXISTING_TBI").unwrap();
        assert!(stage_vcf_tabix_index(&corrupt_vcf_gz, &corrupt_vcf_gz).is_err());
        assert_eq!(
            std::fs::read(&existing_tbi).unwrap(),
            b"KEEP_EXISTING_TBI",
            "failed stage_vcf_tabix_index must preserve existing .tbi"
        );

        let dir_idx_vcf = dir.path().join("dir_idx.vcf.gz");
        {
            use std::io::Write as _;
            let mut bgzf_w =
                noodles::bgzf::io::Writer::new(std::fs::File::create(&dir_idx_vcf).unwrap());
            bgzf_w
                .write_all(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                .unwrap();
            bgzf_w.finish().unwrap();
        }
        let dir_tbi_target = dir.path().join("dir_idx.vcf.gz.tbi");
        let dir_csi_target = dir.path().join("dir_idx.vcf.gz.csi");
        std::fs::create_dir(&dir_tbi_target).unwrap();
        std::fs::create_dir(&dir_csi_target).unwrap();
        assert!(stage_vcf_tabix_index(&dir_idx_vcf, &dir_idx_vcf).is_err());
        assert!(stage_vcf_csi_index(&dir_idx_vcf, &dir_idx_vcf).is_err());
        assert!(dir_tbi_target.is_dir() && dir_csi_target.is_dir());

        let leftover_tmps: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().contains(".tmp."))
            .collect();
        assert!(
            leftover_tmps.is_empty(),
            "no .tmp.<pid>.<seq> files may remain after errors: {leftover_tmps:?}"
        );

        // 4. Stats BrokenPipe maps to ErrorKind::Other (raw_os_error == None) on Linux
        #[cfg(target_os = "linux")]
        {
            use std::os::fd::AsRawFd as _;
            let (reader, writer) = std::io::pipe().unwrap();
            drop(reader); // close read end so writing stats triggers EPIPE (BrokenPipe)
            let stats_fd_path = format!("/proc/self/fd/{}", writer.as_raw_fd());
            let pipe_cfg = crate::config::Config {
                reference: ref_fa.to_string_lossy().into_owned(),
                bam: bam.to_string_lossy().into_owned(),
                region: Some("chr1:1-50".to_string()),
                output_path: Some("/dev/null".to_string()),
                stats_json_path: Some(stats_fd_path),
                threads: 1,
                ..Default::default()
            };
            let err = crate::run(&pipe_cfg).expect_err("stats write to closed pipe must fail");
            match err {
                VardictError::Io(io_err) => {
                    assert_eq!(io_err.kind(), std::io::ErrorKind::Other);
                    assert!(
                        io_err.raw_os_error().is_none(),
                        "raw_os_error must be None so main.rs does not mask stats EPIPE as exit 0"
                    );
                }
                other => panic!("expected VardictError::Io(Other), got {other:?}"),
            }
            drop(writer);
        }

        // 5. record_read_end saturates instead of overflowing i32
        let mut big_rec = crate::bam::SamRecord::new();
        big_rec.cigar = vec![(
            u32::MAX,
            noodles::sam::alignment::record::cigar::op::Kind::Match,
        )];
        assert_eq!(record_read_end(&big_rec, i32::MAX - 10), i32::MAX);
    }
}
