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
    /// Simple/amplicon: adjust_mnp unconditionally, realign_indels under
    /// `-r 1`, then find_all_sv + adj_snv.
    Full,
    /// Somatic BAM pass (SomaticMode.java): adjust_mnp + realign_indels both
    /// gated on `-r 1`, then find_all_sv + adj_snv unconditionally.
    Somatic,
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
) -> Result<crate::variations::VariationData, VardictError> {
    // Task 4 (streaming feed): one lazy pass per BAM, BAM1-then-BAM2, through
    // a single sequential preprocessor into a single CigarParser. Records are
    // dropped promptly after parsing, so peak memory stays flat in depth —
    // the old `raw_records` Vec held the whole tile. Counts are unknown until
    // the feed is exhausted (parse_cigar never reads them mid-stream), so
    // they are set on the parser just before finalization, exactly as
    // `pileup()` received them. `?` propagates I/O errors with the same Err
    // semantics as the Vec path; nothing is skipped or unwrapped.
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
    for mut feed in feeds {
        // Reused record buffer: zero per-read allocation after warmup (the
        // retained-churn slope that OOMs deep tiles).
        let mut r = crate::bam::SamRecord::new();
        while feed.next_into(&mut r)? {
            if preprocessor.preprocess_record(&r, config) {
                parser.parse_cigar(&r)?;
            }
        }
    }
    parser.total_reads = preprocessor.total_reads as i32;
    parser.duplicate_reads = preprocessor.duplicate_reads as i32;
    let mut var_data = parser.into_variation_data();
    let chr_len = *chromosomes_lengths.get(chr).unwrap_or(&0) as u32;

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
            let current_bam = bam_paths
                .first()
                .map(|p| std::path::Path::new(p.as_str()))
                .unwrap_or(std::path::Path::new(""));
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
        PipelineTail::Somatic => {
            if config.realign == 1 {
                crate::realign::adjust_mnp(
                    reference,
                    chr,
                    chr_len,
                    pileup_region,
                    config,
                    &mut var_data,
                );
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
            let current_bam = bam_paths
                .first()
                .map(|p| std::path::Path::new(p.as_str()))
                .unwrap_or(std::path::Path::new(""));
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
/// Defaults to stdout; routes to stderr when `--DP ERR` is passed (ports CmdParser.java:171-177).
pub fn create_output_writer(config: &config::Config) -> OutputDestination {
    if config.debug_print.as_deref() == Some("ERR") {
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

    let writer = create_output_writer(config);
    run_with_writer(config, writer)
}

/// Runs the variant-calling pipeline writing all output to `writer` without clearing
/// process-wide BAM/CRAM/FASTA caches before/after the call.
pub fn run_with_writer<W: std::io::Write>(
    config: &config::Config,
    mut writer: W,
) -> Result<(), VardictError> {
    log::info!("VarDictRust engine starting");

    // Detect somatic (paired) mode: BAM path contains '|'
    if config.bam.contains('|') {
        return run_somatic(config, &mut writer);
    }

    let first_bam = config
        .bam
        .split('|')
        .next()
        .ok_or_else(|| VardictError::Parse("No BAM files specified".to_string()))?;

    // Read reference lengths from cached BAM/CRAM header
    let chromosomes_lengths = crate::bam::read_reference_lengths(Path::new(first_bam), config)?;

    // Parse target regions
    let regions = if let Some(ref r_str) = config.region {
        crate::regions::parse_region_flag(r_str, config, &chromosomes_lengths)?
    } else if let Some(ref bed_path) = config.bed {
        crate::regions::parse_bed(bed_path, config, &chromosomes_lengths, None)?
    } else {
        return Err(VardictError::Parse(
            "Either region (-R) or BED file must be specified".to_string(),
        ));
    };

    // Amplicon mode: process per-amplicon and combine via amplicon post-process
    if config.amplicon.is_some() {
        return run_amplicon(config, &chromosomes_lengths, regions, &mut writer);
    }

    let flat_regions: Vec<crate::regions::Region> = regions.into_iter().flatten().collect();

    // Tile regions into overlapping tiles (regions <= 1 Mbp execute untiled like
    // VarDictJava; multi-megabase intervals > 1 Mbp tile into 100 kbp windows).
    let tile_size = config.tile_size.unwrap_or(DEFAULT_TILE_SIZE);
    let halo = 500;
    let mut tiles = Vec::new();
    for reg in &flat_regions {
        tiles.extend(crate::regions::tile(
            reg,
            effective_tile_size(reg, config.tile_size),
            halo,
        ));
    }
    let tile_batches = crate::regions::batch_tiles(tiles, tile_size, 16, config.threads);

    let first_bam_owned = first_bam.to_string();
    let config_clone = config.clone();
    let chromosomes_lengths_clone = chromosomes_lengths.clone();

    // -h/--header column line (CmdParser.java:67; SimpleMode.java:33 calls
    // printHeader() before any variant output).
    if config.column_header {
        crate::output::table::print_header(&mut writer, config).map_err(VardictError::Io)?;
    }

    // Print VCF header when not in intermediate mode
    if !config.emit_intermediate {
        crate::output::vcf::print_vcf_header(&mut writer, config, &chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    struct SimpleVcfCandidate {
        start_pos: i32,
        frequency: f64,
        is_sv: bool,
        splitreads: i32,
        ref_empty: bool,
        vcf_line: String,
    }

    enum SimpleTileOutput {
        Lines(Vec<String>),
        Vcf {
            chr: String,
            records: Vec<SimpleVcfCandidate>,
        },
    }

    let mut current_chr = String::new();
    let mut pending: std::collections::BTreeMap<i32, SimpleVcfCandidate> =
        std::collections::BTreeMap::new();

    // Process tiles in parallel, print variants in deterministic index order
    crate::runner::run_parallel_batched(
        tile_batches,
        config.threads,
        move |tile| -> Result<SimpleTileOutput, VardictError> {
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
            let var_data = run_pipeline_on_records(
                std::slice::from_ref(&first_bam_owned),
                &tile.chr,
                tile.halo_start.to_zero_based(),
                crate::types::Pos0(tile.halo_end.0),
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
                return Ok(SimpleTileOutput::Lines(lines));
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
                        line.push_str(&format!("\t{}", v.crispr));
                    }
                    // ports SimpleOutputVariant.java:116-117 (DEBUG appended after crispr)
                    if config_clone.debug {
                        line.push('\t');
                        line.push_str(&v.debug);
                    }
                    lines.push(line);
                }
                Ok(SimpleTileOutput::Lines(lines))
            } else {
                let mut records = Vec::new();
                for v in variants {
                    let ref_empty = crate::output::vcf::perl_row_skipped(&v.ref_allele);
                    let is_sv = v.var_allele.contains('<');
                    let mut splitreads = 0;
                    if !v.sv.is_empty() && v.sv.contains('-') {
                        let parts: Vec<&str> = v.sv.split('-').collect();
                        if !parts.is_empty() {
                            splitreads = parts[0].parse::<i32>().unwrap_or(0);
                        }
                    }
                    let vcf_line = if ref_empty || (is_sv && splitreads < 1) {
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
                    records.push(SimpleVcfCandidate {
                        start_pos: v.start_pos,
                        frequency: v.frequency,
                        is_sv,
                        splitreads,
                        ref_empty,
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
                SimpleTileOutput::Lines(lines) => {
                    for line in lines {
                        write_output_line(&mut writer, &line)?;
                    }
                }
                SimpleTileOutput::Vcf { chr, records } => {
                    if chr != current_chr {
                        for (_, cand) in std::mem::take(&mut pending) {
                            if !cand.ref_empty && (!cand.is_sv || cand.splitreads >= 1) {
                                write_output_line(&mut writer, &cand.vcf_line)?;
                            }
                        }
                        current_chr = chr;
                    }

                    for r in records {
                        match pending.entry(r.start_pos) {
                            std::collections::btree_map::Entry::Vacant(e) => {
                                e.insert(r);
                            }
                            std::collections::btree_map::Entry::Occupied(mut e) => {
                                if r.frequency > e.get().frequency {
                                    e.insert(r);
                                }
                            }
                        }
                    }
                }
            }
            Ok(())
        },
    )?;

    for (_, cand) in pending {
        if !cand.ref_empty && (!cand.is_sv || cand.splitreads >= 1) {
            write_output_line(&mut writer, &cand.vcf_line)?;
        }
    }
    writer.flush().map_err(VardictError::Io)?;

    Ok(())
}

/// Process a single amplicon region through the full pipeline, returning the
/// aligned variants map. The pileup runs over the amplicon region (start..end);
/// the ToVarsBuilder calls variants over the insert region.
fn process_amplicon_region(
    bam: &str,
    region: &crate::regions::Region,
    config: &config::Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<crate::javamap::JavaHashMap<i32, crate::variations::Vars>, VardictError> {
    let halo = 500u32;
    let halo_start = crate::types::Pos0(region.start.0.saturating_sub(halo).saturating_sub(1));
    let halo_end = crate::types::Pos0(region.end.0 + halo);

    // Reference is loaded over the halo bounds; the pileup runs over the
    // amplicon region itself.
    let ref_region = crate::regions::Region {
        chr: region.chr.clone(),
        start: crate::types::Pos1(halo_start.0 + 1),
        end: crate::types::Pos1(halo_end.0),
        gene: region.gene.clone(),
        amplicon: region.amplicon.clone(),
    };

    let mut reference = crate::reference::Reference::new();
    let bam_string = bam.to_string();
    let var_data = run_pipeline_on_records(
        std::slice::from_ref(&bam_string),
        &region.chr,
        halo_start,
        halo_end,
        Some(&ref_region),
        &mut reference,
        region,
        PipelineTail::Full,
        config,
        chromosomes_lengths,
        0,
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
    Ok(builder.build_aligned_variants(&var_data))
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
fn run_amplicon<W: std::io::Write>(
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
        crate::output::table::print_header(writer, config).map_err(VardictError::Io)?;
    }

    // Print VCF header once when not in intermediate mode (same var2vcf_valid.pl
    // header as simple mode — amplicon rides the same writer, perl :140-145).
    if !config.emit_intermediate {
        crate::output::vcf::print_vcf_header(writer, config, chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    for segment in segments {
        // §9: let-else merges the empty guard with the last() unwrap.
        let Some(current_region) = segment.last() else {
            continue;
        };
        let mut vars_per_amplicon: Vec<crate::javamap::JavaHashMap<i32, crate::variations::Vars>> =
            Vec::with_capacity(segment.len());
        let splice = std::collections::HashSet::new();

        for region in segment.iter() {
            let aligned = process_amplicon_region(first_bam, region, config, chromosomes_lengths)?;
            vars_per_amplicon.push(aligned);
        }

        // Position → covering-amplicon pairs (AmpliconMode.java:57-60).
        //
        // -p emits a (null/ref) row at EVERY insert position, even with zero
        // coverage in all amplicons, so visiting every insert position is
        // load-bearing here — but the per-base map is streamed, never
        // materialized (audit-amplicon-expand): the lazy union iterator yields
        // exactly the positions (and segment-ordered covering lists) the old
        // per-base BTreeMap held, in the same ascending order.
        if config.pileup && config.emit_intermediate {
            // Fully streaming -p table path: each row goes straight to the
            // writer, so working memory stays O(segment + one position) no
            // matter how large the insert interval is.
            crate::postprocess::amplicon::amplicon_post_process_streaming(
                config,
                current_region,
                &vars_per_amplicon,
                crate::postprocess::amplicon::amplicon_position_stream(&segment),
                &splice,
                &mut |av| {
                    let line = if config.fisher {
                        crate::output::table::amplicon_40(
                            &av,
                            current_region,
                            &config.sample,
                            config.debug,
                        )
                    } else {
                        crate::output::table::amplicon_38(
                            &av,
                            current_region,
                            &config.sample,
                            config.debug,
                        )
                    };
                    write_output_line(writer, &line)
                },
            )?;
            continue;
        }

        let results = if config.pileup {
            // -p VCF path: same rows as the old per-base expansion, fed by the
            // lazy union iterator (no per-base map). Row collection stays: the
            // var2vcf_valid.pl-equivalent top-1-per-start selection below is
            // global, so these rows cannot go straight to the writer.
            let mut rows = Vec::new();
            {
                let mut collect =
                    |av: crate::postprocess::amplicon::AmpliconVariant| -> Result<(), VardictError> {
                        rows.push(av);
                        Ok(())
                    };
                crate::postprocess::amplicon::amplicon_post_process_streaming(
                    config,
                    current_region,
                    &vars_per_amplicon,
                    crate::postprocess::amplicon::amplicon_position_stream(&segment),
                    &splice,
                    &mut collect,
                )?;
            }
            crate::postprocess::amplicon::amplicon_select_vcf_rows(rows)
        } else {
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

            crate::postprocess::amplicon::amplicon_post_process(
                config,
                current_region,
                &vars_per_amplicon,
                &amplicons_on_positions,
                &splice,
            )?
        };

        for av in results {
            if config.emit_intermediate {
                let line = if config.fisher {
                    crate::output::table::amplicon_40(
                        &av,
                        current_region,
                        &config.sample,
                        config.debug,
                    )
                } else {
                    crate::output::table::amplicon_38(
                        &av,
                        current_region,
                        &config.sample,
                        config.debug,
                    )
                };
                write_output_line(writer, &line)?;
            } else {
                // VCF mode: same writer as simple mode (var2vcf_valid.pl), with the
                // amplicon columns (37-40) selecting $isamp behavior (:145).
                let v = &av.variant;
                if crate::output::vcf::perl_row_skipped(&v.ref_allele) {
                    // ports var2vcf_valid.pl:141 `next unless ($ref)`; rows
                    // arriving here are already AF-selected winners — the
                    // amplicon.rs port of valid.pl:129-137 (bucket by raw
                    // start, keep top row) ran upstream of emission, so a
                    // pre-emit skip equals a post-selection skip.
                    continue;
                }
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
                let filters =
                    crate::output::vcf::filter_variant(v, config, 0, amp, sbf, &oddratio_str);
                let record = crate::output::vcf::format_vcf_record(
                    v,
                    current_region,
                    config,
                    &filters,
                    sbf,
                    &oddratio_str,
                    amp,
                );
                write_output_line(writer, &record)?;
            }
        }
    }

    writer.flush().map_err(VardictError::Io)?;

    Ok(())
}

/// Somatic (paired) mode: process tumor BAM, then normal BAM, then compare.
/// Ports SomaticMode.java's processBothBamsInPipeline.
fn run_somatic<W: std::io::Write>(
    config: &config::Config,
    writer: &mut W,
) -> Result<(), VardictError> {
    let parts: Vec<&str> = config.bam.split('|').collect();
    if parts.len() != 2 {
        return Err(VardictError::Parse(
            "Somatic mode requires exactly two BAM files separated by '|'".to_string(),
        ));
    }
    let bam1_path = parts[0].to_string();
    let bam2_path = parts[1].to_string();

    // Read reference lengths from cached BAM/CRAM header
    let chromosomes_lengths = crate::bam::read_reference_lengths(Path::new(&bam1_path), config)?;

    // Parse target regions
    let regions = if let Some(ref r_str) = config.region {
        crate::regions::parse_region_flag(r_str, config, &chromosomes_lengths)?
    } else if let Some(ref bed_path) = config.bed {
        crate::regions::parse_bed(bed_path, config, &chromosomes_lengths, None)?
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
    for reg in &flat_regions {
        tiles.extend(crate::regions::tile(
            reg,
            effective_tile_size(reg, config.tile_size),
            halo,
        ));
    }
    let tile_batches = crate::regions::batch_tiles(tiles, tile_size, 16, config.threads);

    let config_clone = config.clone();
    let chromosomes_lengths_clone = chromosomes_lengths.clone();
    let bam1_path_for_combine = bam1_path.clone();
    let bam2_path_for_combine = bam2_path.clone();

    // -h/--header column line (SomaticMode.java:35, 134-150).
    if config.column_header {
        crate::output::table::print_header(writer, config).map_err(VardictError::Io)?;
    }

    // Paired VCF header when not in intermediate mode
    if !config.emit_intermediate {
        crate::output::vcf::print_paired_vcf_header(writer, config, &chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    // Process tiles: run both BAMs through pipeline, then somatic postprocess
    crate::runner::run_parallel_batched(
        tile_batches,
        config.threads,
        move |tile| -> Result<Vec<String>, VardictError> {
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

            let mut reference = crate::reference::Reference::new();
            let mut var_data1 = run_pipeline_on_records(
                std::slice::from_ref(&bam1_path),
                &tile.chr,
                tile.halo_start.to_zero_based(),
                crate::types::Pos0(tile.halo_end.0),
                Some(&bed_region),
                &mut reference,
                &region_for_pileup,
                PipelineTail::Somatic,
                &config_clone,
                &chromosomes_lengths_clone,
                0,
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
            // Process BAM2 (normal), seeded with BAM1's maxReadLength; the
            // reference loaded above is reused (no second get_reference).
            let var_data2 = run_pipeline_on_records(
                std::slice::from_ref(&bam2_path),
                &tile.chr,
                tile.halo_start.to_zero_based(),
                crate::types::Pos0(tile.halo_end.0),
                None,
                &mut reference,
                &region_for_pileup,
                PipelineTail::Somatic,
                &config_clone,
                &chromosomes_lengths_clone,
                max_read_length_bam1,
            )?;

            // Use max of both BAMs' maxReadLength
            let _combined_max_rl = max_read_length_bam1.max(var_data2.max_read_length);

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
            let bam1 = bam1_path_for_combine.clone();
            let bam2 = bam2_path_for_combine.clone();
            let config_for_combine = config_clone.clone();
            let chr_lens_for_combine = chromosomes_lengths_clone.clone();
            let pipeline_fn = |combine_region: &crate::regions::Region| -> Result<
                crate::javamap::JavaHashMap<i32, crate::variations::Vars>,
                VardictError,
            > {
                // Query reads from both BAMs over the combine region; the
                // merged record set goes through one preprocessing pass.
                // The working reference is a per-call clone so this callback
                // stays a shared-borrow (Fn) capture, like Java's.
                let bams = vec![bam1.clone(), bam2.clone()];
                let mut pipeline_ref = reference.clone();
                let var_data = run_pipeline_on_records(
                    &bams,
                    &combine_region.chr,
                    crate::types::Pos0(combine_region.start.0.saturating_sub(1)),
                    crate::types::Pos0(combine_region.end.0),
                    None,
                    &mut pipeline_ref,
                    combine_region,
                    PipelineTail::Somatic,
                    &config_for_combine,
                    &chr_lens_for_combine,
                    0,
                )?;

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
                &splice1,
                Some(&pipeline_fn),
                max_read_length_bam1,
            );

            let mut lines = Vec::new();

            // Output 55-column format (or 61-column with --fisher), or the
            // paired VCF (var2vcf_paired.pl) when not in intermediate mode.
            if config_clone.emit_intermediate {
                for sv in somatic_variants {
                    let line = if config_clone.fisher {
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
                        lines.push(format!("{}\t{}", line, sv.variant1.debug));
                    } else {
                        lines.push(line);
                    }
                }
            } else {
                // Paired VCF mode: ports var2vcf_paired.pl:41,128-129,136 —
                // bucket rows by the raw .var start column ($hash{$chr}->{$a[3]}),
                // AF-contest per bucket (:128-129, top row only without -A), then
                // the `next unless ($ref)` guard on each SELECTED row.
                let records: Vec<crate::output::vcf::PairedVcfRecord> =
                    select_paired_rows(somatic_variants)
                        .into_iter()
                        .map(|sv| {
                            crate::output::vcf::format_paired_vcf_record(
                                &sv,
                                &builder_region.chr,
                                &config_clone,
                            )
                        })
                        .collect();
                // Delayed print + Cluster state machine (:229-236,248-250);
                // opt_c default 0 (no -c plumbed to the paired writer) makes the
                // rule inert.
                let finished = crate::output::vcf::finish_paired_records(records, 0);
                lines.extend(finished);
            }
            Ok(lines)
        },
        |lines| {
            for line in lines {
                write_output_line(writer, &line)?;
            }
            Ok(())
        },
    )?;

    writer.flush().map_err(VardictError::Io)?;

    Ok(())
}

/// Paired-VCF row selection: bucket by tumor-variant start (var2vcf_paired.pl
/// `$hash{$chr}->{$a[3]}`), keep the highest-AF row per position (`$ALL = 1`
/// without -A; replace-only-if-strictly-greater keeps first-seen ties like the
/// stable mergesort at :128), THEN drop Perl-falsy REF rows (:136 `next unless
/// ($ref)` runs on each selected row inside the $ALL loop — an empty-ref winner
/// means NOTHING is emitted for that position).
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
                if sv.variant1.frequency > e.get().variant1.frequency {
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
}
