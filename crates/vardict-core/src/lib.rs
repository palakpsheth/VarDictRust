//! VarDictRust core engine library.
//!
//! An idiomatic Rust rewrite of VarDictJava v1.8.3 — a sensitive somatic and
//! germline variant caller for next-generation sequencing data.
//!
//! This crate contains all variant-calling logic. The binary crate
//! (`vardict-cli`) provides the command-line interface.
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
    let mut raw_records = Vec::new();
    for path in bam_paths {
        raw_records.extend(crate::bam::query_region(
            Path::new(path),
            chr,
            query_start,
            query_end,
            config,
        )?);
    }

    let mut preprocessor = crate::bam::RecordPreprocessor::new();
    let mut records = Vec::new();
    for r in raw_records {
        if preprocessor.preprocess_record(&r, config) {
            records.push(r);
        }
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

    let mut var_data = crate::cigar::parse::pileup(
        &records,
        reference,
        pileup_region,
        config,
        false,
        preprocessor.total_reads as i32,
        preprocessor.duplicate_reads as i32,
        seed_max_read_length,
    )?;

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

/// Main entry point for the engine. Called by the CLI binary.
///
/// Takes a fully-resolved configuration and runs the variant-calling pipeline.
pub fn run(config: &config::Config) -> Result<(), VardictError> {
    log::info!("VarDictRust engine starting");

    // Detect somatic (paired) mode: BAM path contains '|'
    if config.bam.contains('|') {
        return run_somatic(config);
    }

    let first_bam = config
        .bam
        .split('|')
        .next()
        .ok_or_else(|| VardictError::Parse("No BAM files specified".to_string()))?;

    // Open first BAM to read references lengths
    let cached_reader = crate::bam::open_reader(Path::new(first_bam), config)?;
    let mut chromosomes_lengths = HashMap::new();
    for (name, ref_seq) in cached_reader.header.reference_sequences().iter() {
        let name_str = String::from_utf8_lossy(name).into_owned();
        let len = ref_seq.length().get();
        chromosomes_lengths.insert(name_str, len);
    }

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
        return run_amplicon(config, &chromosomes_lengths, regions);
    }

    let flat_regions: Vec<crate::regions::Region> = regions.into_iter().flatten().collect();

    // Tile regions into overlapping tiles
    let tile_size = config.tile_size.unwrap_or(100_000);
    let halo = 500;
    let mut tiles = Vec::new();
    for reg in &flat_regions {
        tiles.extend(crate::regions::tile(reg, tile_size, halo));
    }

    let first_bam_owned = first_bam.to_string();
    let config_clone = config.clone();
    let chromosomes_lengths_clone = chromosomes_lengths.clone();

    // -h/--header column line (CmdParser.java:67; SimpleMode.java:33 calls
    // printHeader() before any variant output).
    if config.column_header {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        crate::output::table::print_header(&mut handle, config).map_err(VardictError::Io)?;
    }

    // Print VCF header when not in intermediate mode
    if !config.emit_intermediate {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        crate::output::vcf::print_vcf_header(&mut handle, config, &chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    // Process tiles in parallel, print variants in deterministic index order
    crate::runner::run_parallel(
        tiles,
        config.threads,
        move |tile| -> Result<Vec<String>, VardictError> {
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
            let var_data = run_pipeline_on_records(
                std::slice::from_ref(&first_bam_owned),
                &tile.chr,
                tile.halo_start.to_zero_based(),
                crate::types::Pos0(tile.halo_end.0),
                Some(&bed_region),
                &mut reference,
                &region_for_pileup,
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
                for (intron, count_arr) in &var_data.splice_count {
                    let count = count_arr[0];
                    lines.push(format!(
                        "{}\t{}\t{}\t{}",
                        config_clone.sample, tile.chr, intron, count
                    ));
                }
                return Ok(lines);
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

            for v in variants {
                if config_clone.emit_intermediate {
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
                } else {
                    // ports var2vcf_valid.pl:141 — empty-REF .var rows (the -p
                    // null records) are dropped before VCF conversion. Java's
                    // vars list can hold several rows per position, but rows
                    // arriving here are already AF-selected winners: the
                    // postprocess/simple.rs port of valid.pl:129-137 (bucket
                    // by raw start, keep top row) ran upstream of emission,
                    // so a pre-emit skip equals a post-selection skip.
                    if crate::output::vcf::perl_row_skipped(&v.ref_allele) {
                        continue;
                    }
                    // VCF mode: always compute Fisher for SBF field
                    // (var2vcf_valid.pl never DROPS on AF without -S; it only
                    // annotates the f$Freq filter, :159)
                    let fisher = crate::stats::FisherExact::new(
                        v.ref_fwd_cov as i32,
                        v.ref_rev_cov as i32,
                        v.fwd as i32,
                        v.rev as i32,
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
                    let record = crate::output::vcf::format_vcf_record(
                        &v,
                        &builder_region,
                        &config_clone,
                        &filters,
                        sbf,
                        &oddratio_str,
                        None,
                    );
                    lines.push(record);
                }
            }
            Ok(lines)
        },
        move |lines| {
            let stdout = std::io::stdout();
            let mut handle = stdout.lock();
            use std::io::Write;
            for line in lines {
                writeln!(handle, "{}", line).map_err(VardictError::Io)?;
            }
            Ok(())
        },
    )?;

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

/// Amplicon mode: process each amplicon region separately, then combine via
/// amplicon post-processing. Ports AmpliconMode.java.
fn run_amplicon(
    config: &config::Config,
    chromosomes_lengths: &HashMap<String, usize>,
    segments: Vec<Vec<crate::regions::Region>>,
) -> Result<(), VardictError> {
    let first_bam = config
        .bam
        .split('|')
        .next()
        .ok_or_else(|| VardictError::Parse("No BAM files specified".to_string()))?;

    // -h/--header column line (AmpliconMode.java:39, 135-146).
    if config.column_header {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        crate::output::table::print_header(&mut handle, config).map_err(VardictError::Io)?;
    }

    // Print VCF header once when not in intermediate mode (same var2vcf_valid.pl
    // header as simple mode — amplicon rides the same writer, perl :140-145).
    if !config.emit_intermediate {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        crate::output::vcf::print_vcf_header(&mut handle, config, chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    for segment in segments {
        // §9: let-else merges the empty guard with the last() unwrap.
        let Some(current_region) = segment.last() else {
            continue;
        };
        let mut amplicons_on_positions: std::collections::BTreeMap<
            i32,
            Vec<(usize, crate::regions::Region)>,
        > = std::collections::BTreeMap::new();
        let mut vars_per_amplicon: Vec<crate::javamap::JavaHashMap<i32, crate::variations::Vars>> =
            Vec::new();
        let splice = std::collections::HashSet::new();

        for (amplicon_number, region) in segment.iter().enumerate() {
            if let Some(ref amp) = region.amplicon {
                for p in amp.insert_start.0 as i32..=amp.insert_end.0 as i32 {
                    amplicons_on_positions
                        .entry(p)
                        .or_default()
                        .push((amplicon_number, region.clone()));
                }
            }
            let aligned = process_amplicon_region(first_bam, region, config, chromosomes_lengths)?;
            vars_per_amplicon.push(aligned);
        }

        let results = crate::postprocess::amplicon::amplicon_post_process(
            config,
            current_region,
            &vars_per_amplicon,
            &amplicons_on_positions,
            &splice,
        );

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
                println!("{}", line);
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
                let fisher = crate::stats::FisherExact::new(
                    v.ref_fwd_cov as i32,
                    v.ref_rev_cov as i32,
                    v.fwd as i32,
                    v.rev as i32,
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
                println!("{}", record);
            }
        }
    }

    Ok(())
}

/// Somatic (paired) mode: process tumor BAM, then normal BAM, then compare.
/// Ports SomaticMode.java's processBothBamsInPipeline.
fn run_somatic(config: &config::Config) -> Result<(), VardictError> {
    let parts: Vec<&str> = config.bam.split('|').collect();
    if parts.len() != 2 {
        return Err(VardictError::Parse(
            "Somatic mode requires exactly two BAM files separated by '|'".to_string(),
        ));
    }
    let bam1_path = parts[0].to_string();
    let bam2_path = parts[1].to_string();

    // Open first BAM to read reference lengths
    let cached_reader = crate::bam::open_reader(Path::new(&bam1_path), config)?;
    let mut chromosomes_lengths = HashMap::new();
    for (name, ref_seq) in cached_reader.header.reference_sequences().iter() {
        let name_str = String::from_utf8_lossy(name).into_owned();
        let len = ref_seq.length().get();
        chromosomes_lengths.insert(name_str, len);
    }

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

    // Tile regions
    let tile_size = config.tile_size.unwrap_or(100_000);
    let halo = 500;
    let mut tiles = Vec::new();
    for reg in &flat_regions {
        tiles.extend(crate::regions::tile(reg, tile_size, halo));
    }

    let config_clone = config.clone();
    let chromosomes_lengths_clone = chromosomes_lengths.clone();
    let bam1_path_for_combine = bam1_path.clone();
    let bam2_path_for_combine = bam2_path.clone();

    // -h/--header column line (SomaticMode.java:35, 134-150).
    if config.column_header {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        crate::output::table::print_header(&mut handle, config).map_err(VardictError::Io)?;
    }

    // Paired VCF header when not in intermediate mode
    if !config.emit_intermediate {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        crate::output::vcf::print_paired_vcf_header(&mut handle, config, &chromosomes_lengths)
            .map_err(VardictError::Io)?;
    }

    // Process tiles: run both BAMs through pipeline, then somatic postprocess
    crate::runner::run_parallel(
        tiles,
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
            let var_data1 = run_pipeline_on_records(
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

            let builder_region = crate::regions::Region {
                chr: tile.chr.clone(),
                start: tile.original_start,
                end: tile.original_end,
                gene: tile.gene.clone(),
                amplicon: tile.amplicon.clone(),
            };

            // Build aligned variants for both BAMs
            let builder1 = crate::tovars::ToVarsBuilder::new(
                &reference,
                &builder_region,
                &chromosomes_lengths_clone,
                &config_clone,
            );
            let aligned1 = builder1.build_aligned_variants(&var_data1);

            let builder2 = crate::tovars::ToVarsBuilder::new(
                &reference,
                &builder_region,
                &chromosomes_lengths_clone,
                &config_clone,
            );
            let aligned2 = builder2.build_aligned_variants(&var_data2);

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
                &var_data1.splice,
                Some(&pipeline_fn),
                var_data1.max_read_length,
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
        move |lines| {
            let stdout = std::io::stdout();
            let mut handle = stdout.lock();
            use std::io::Write;
            for line in lines {
                writeln!(handle, "{}", line).map_err(VardictError::Io)?;
            }
            Ok(())
        },
    )?;

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
}
