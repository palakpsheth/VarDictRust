//! BED region parsing and region tiling logic.
//!
//! Ports Region.java and RegionBuilder.java.

use crate::config::Config;
use crate::types::{Pos1, VardictError};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Information for amplicon-based variant calling.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AmpliconInfo {
    pub insert_start: Pos1,
    pub insert_end: Pos1,
}

/// Represents a genomic region parsed from BED or command line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub chr: String,
    pub start: Pos1,
    pub end: Pos1,
    pub gene: String,
    pub amplicon: Option<AmpliconInfo>,
}

/// Represents a tiled sub-interval of a Region with active boundaries and surrounding overlap halo.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tile {
    pub chr: String,
    pub start: Pos1,
    pub end: Pos1,
    pub gene: String,
    pub halo_start: Pos1,
    pub halo_end: Pos1,
    pub amplicon: Option<AmpliconInfo>,
    pub original_start: Pos1,
    pub original_end: Pos1,
}

/// Helper struct for BED row column indexes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BedRowFormat {
    pub chr_column: usize,
    pub start_column: usize,
    pub end_column: usize,
    pub thick_start_column: usize,
    pub thick_end_column: usize,
    pub gene_column: usize,
}

/// Custom format used when 4-column BED file is auto-detected.
/// Note: VarDictJava has a known bug/design where it sets thick_start to column 3 (gene),
/// thick_end to column 1 (start), and gene_column to column 2 (end). We emulate this exactly!
pub const CUSTOM_BED_ROW_FORMAT: BedRowFormat = BedRowFormat {
    chr_column: 0,
    start_column: 1,
    end_column: 2,
    thick_start_column: 3,
    thick_end_column: 1,
    gene_column: 2,
};

/// Format used for amplicon BED files.
pub const AMP_BED_ROW_FORMAT: BedRowFormat = BedRowFormat {
    chr_column: 0,
    start_column: 1,
    end_column: 2,
    thick_start_column: 6,
    thick_end_column: 7,
    gene_column: 3,
};

/// Corrects chromosome name according to active chromosomes in the BAM header.
pub fn correct_chromosome(chromosomes_lengths: &HashMap<String, usize>, chr: &str) -> String {
    if chromosomes_lengths.contains_key(chr) {
        chr.to_string()
    } else if let Some(stripped) = chr.strip_prefix("chr") {
        stripped.to_string()
    } else {
        format!("chr{}", chr)
    }
}

/// Tils a Region into overlapping windows (with a surrounding halo context).
pub fn tile(region: &Region, tile_size: u32, halo: u32) -> Vec<Tile> {
    let tile_size = tile_size.max(1);
    let mut tiles = Vec::new();
    let reg_start = region.start.0;
    let reg_end = region.end.0;

    let mut curr_start = reg_start;
    while curr_start <= reg_end {
        // Subtract first: tile_size >= 1 so this never underflows, and the
        // terminal tile reaches reg_end (saturating add-then-sub would stall
        // one base short at u32::MAX and spin forever).
        let mut curr_end = curr_start.saturating_add(tile_size.saturating_sub(1));
        if curr_end > reg_end {
            curr_end = reg_end;
        }

        let mut halo_start_val = if curr_start > halo {
            curr_start - halo
        } else {
            1
        };
        if halo_start_val < reg_start {
            halo_start_val = reg_start;
        }

        let mut halo_end_val = curr_end.saturating_add(halo);
        if halo_end_val > reg_end {
            halo_end_val = reg_end;
        }

        tiles.push(Tile {
            chr: region.chr.clone(),
            start: Pos1(curr_start),
            end: Pos1(curr_end),
            gene: region.gene.clone(),
            halo_start: Pos1(halo_start_val),
            halo_end: Pos1(halo_end_val),
            amplicon: region.amplicon.clone(),
            original_start: region.start,
            original_end: region.end,
        });

        if curr_end >= reg_end {
            break;
        }

        curr_start = curr_start.saturating_add(tile_size);
    }

    tiles
}

/// Groups consecutive small tiles on the same chromosome into coarse-grained
/// batches (`EBO-3`) so high-region-count Exome BEDs (~300k exons) amortize
/// Rayon/Condvar/channel scheduling and preserve per-thread BGZF block locality (`M-THROUGHPUT`).
///
/// Invariants:
/// - Effective batch size is capped at `(tiles.len() / (num_threads.max(1) * 4)).clamp(1, max_batch_tiles)`
///   so small BEDs (`tiles.len() <= num_threads * 4`) never merge tiles and
///   always keep all worker threads busy.
/// - Split sub-tiles of a large region (`tile.start != tile.original_start || tile.end != tile.original_end`)
///   or tiles whose cumulative base span would exceed `target_batch_bp` are never merged.
/// - Batches never cross chromosome boundaries (`curr.chr != prev.chr`).
pub fn batch_tiles(
    tiles: Vec<Tile>,
    target_batch_bp: u32,
    max_batch_tiles: usize,
    num_threads: usize,
) -> Vec<Vec<Tile>> {
    if tiles.is_empty() {
        return Vec::new();
    }
    let min_batches = num_threads.max(1).saturating_mul(4);
    let effective_max = (tiles.len() / min_batches).clamp(1, max_batch_tiles.max(1));
    if effective_max == 1 {
        return tiles.into_iter().map(|t| vec![t]).collect();
    }

    let mut batches: Vec<Vec<Tile>> = Vec::with_capacity(tiles.len() / effective_max + 1);
    let mut current_batch: Vec<Tile> = Vec::with_capacity(effective_max);
    let mut current_bp: u32 = 0;

    for t in tiles {
        let tile_bp = t
            .halo_end
            .0
            .saturating_sub(t.halo_start.0)
            .saturating_add(1);
        let is_split_tile = t.start != t.original_start || t.end != t.original_end;

        let can_append = !current_batch.is_empty()
            && !is_split_tile
            && current_batch.len() < effective_max
            && current_bp.saturating_add(tile_bp) <= target_batch_bp
            && current_batch.last().is_some_and(|prev| {
                prev.chr == t.chr
                    && prev.start == prev.original_start
                    && prev.end == prev.original_end
            });

        if can_append {
            current_bp = current_bp.saturating_add(tile_bp);
            current_batch.push(t);
        } else {
            if !current_batch.is_empty() {
                batches.push(std::mem::replace(
                    &mut current_batch,
                    Vec::with_capacity(effective_max),
                ));
            }
            current_bp = tile_bp;
            current_batch.push(t);
        }
    }

    if !current_batch.is_empty() {
        batches.push(current_batch);
    }

    batches
}

/// Helper to parse a non-negative `i32`-bounded coordinate from string with comma removal.
///
/// Ports `RegionBuilder.java` (`Integer.parseInt(s.replaceAll(",", ""))`): rejects values
/// above `i32::MAX` (`2_147_483_647`) instead of silently wrapping to negative `i32`.
pub fn parse_int_with_commas(s: &str) -> Result<u32, VardictError> {
    let clean = s.replace(',', "");
    let v = clean
        .trim()
        .parse::<u32>()
        .map_err(|e| VardictError::Parse(format!("invalid integer '{}': {}", s, e)))?;
    if v > i32::MAX as u32 {
        return Err(VardictError::Parse(format!(
            "invalid integer '{}': exceeds 2147483647",
            s
        )));
    }
    Ok(v)
}

/// Ports `String.split(delimiter)` (`java.lang.String.split(regex, 0)`):
/// splits `line` by regex `delim` and discards trailing empty segments when `parts.len() > 1`.
pub fn java_split<'a>(line: &'a str, delim: &str) -> Result<Vec<&'a str>, VardictError> {
    let mut parts: Vec<&'a str> = if delim == "\t" || delim == r"\t" {
        line.split('\t').collect()
    } else if delim.len() == 1
        && !delim.starts_with([
            '.', '$', '^', '{', '[', '(', '|', ')', '*', '+', '?', '\\', '}', ']',
        ])
    {
        line.split(delim).collect()
    } else {
        let re = regex::RegexBuilder::new(delim)
            .unicode(false)
            .build()
            .map_err(|e| VardictError::Parse(format!("invalid delimiter regex '{delim}': {e}")))?;
        re.split(line).collect()
    };
    if parts.len() > 1 {
        while parts.last().is_some_and(|s| s.is_empty()) {
            parts.pop();
        }
    }
    Ok(parts)
}

/// Ports `VarDictLauncher.java:139-176`: scans non-header (`#`, `browser`, `track`) lines of a
/// BED reader and returns `true` iff any line has 8 columns with integer columns `1, 2, 6, 7`
/// satisfying `a[6] >= a[1]` and `a[7] <= a[2]`.
pub fn detect_amplicon_bed_from_reader<R: BufRead>(
    reader: R,
    delimiter: &str,
) -> Result<bool, VardictError> {
    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('#') || line.starts_with("browser") || line.starts_with("track") {
            continue;
        }
        let a = java_split(&line, delimiter)?;
        if let [_, a1, a2, _, _, _, a6, a7] = a.as_slice() {
            if let (Ok(s1), Ok(e1), Ok(s2), Ok(e2)) = (
                parse_int_with_commas(a1),
                parse_int_with_commas(a2),
                parse_int_with_commas(a6),
                parse_int_with_commas(a7),
            ) {
                if s2 >= s1 && e2 <= e1 {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

/// File-path wrapper around [`detect_amplicon_bed_from_reader`].
pub fn detect_amplicon_bed<P: AsRef<Path>>(path: P, delimiter: &str) -> Result<bool, VardictError> {
    let file = File::open(path)?;
    detect_amplicon_bed_from_reader(BufReader::new(file), delimiter)
}

/// Parses BED file in non-amplicon mode.
pub fn parse_non_amp_bed_from_reader<R: BufRead>(
    reader: R,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    zero_based_override: Option<bool>,
) -> Result<Vec<Vec<Region>>, VardictError> {
    let is_zero_based = config.zero_based;
    let format = BedRowFormat {
        chr_column: config.col_chr,
        start_column: config.col_start,
        end_column: config.col_end,
        thick_start_column: config.col_seg_start,
        thick_end_column: config.col_seg_end,
        gene_column: config.col_gene.unwrap_or(usize::MAX),
    };

    let mut segs = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty()
            || line.starts_with('#')
            || line.starts_with("browser")
            || line.starts_with("track")
        {
            continue;
        }

        let cols = java_split(&line, &config.delimiter)?;

        let mut current_format = format;
        let mut current_zero_based = is_zero_based;

        // Auto-detect 4-column BED format (RegionBuilder.java:58-68)
        let effective_override =
            zero_based_override.or_else(|| config.zero_based_defined.then_some(config.zero_based));
        if let (None, [_, c1, c2, _]) = (config.column_for_chromosome, cols.as_slice()) {
            let a1 = parse_int_with_commas(c1)?;
            let a2 = parse_int_with_commas(c2)?;
            if a1 <= a2 {
                current_format = CUSTOM_BED_ROW_FORMAT;
                if effective_override.is_none() {
                    current_zero_based = true;
                }
            }
        }

        let active_zero_based = effective_override.unwrap_or(current_zero_based);

        if cols.len() <= current_format.chr_column
            || cols.len() <= current_format.start_column
            || cols.len() <= current_format.end_column
            || cols.len() <= current_format.thick_start_column
            || cols.len() <= current_format.thick_end_column
        {
            return Err(VardictError::Parse(format!(
                "Region line has {} columns, expected at least {} (for 3-, 5-, or 6-column non-amplicon BED files, pass -c 1 -S 2 -E 3 -z for standard 0-based BED or -c 1 -S 2 -E 3 for 1-based regions)",
                cols.len(),
                current_format.thick_end_column.saturating_add(1)
            )));
        }

        let chr = correct_chromosome(chromosomes_lengths, cols[current_format.chr_column]);
        let cds_start = parse_int_with_commas(cols[current_format.start_column])?;
        let cds_end = parse_int_with_commas(cols[current_format.end_column])?;

        let gene = if current_format.gene_column < cols.len() {
            cols[current_format.gene_column].to_string()
        } else {
            chr.clone()
        };

        let thick_starts: Vec<&str> = cols[current_format.thick_start_column].split(',').collect();
        let thick_ends: Vec<&str> = cols[current_format.thick_end_column].split(',').collect();

        let mut thick_regions = Vec::new();

        for i in 0..thick_starts.len() {
            if i >= thick_ends.len() {
                break;
            }
            if thick_starts[i].trim().is_empty() || thick_ends[i].trim().is_empty() {
                continue;
            }

            let thick_start_val = parse_int_with_commas(thick_starts[i])?;
            let thick_end_val = parse_int_with_commas(thick_ends[i])?;

            if cds_start > thick_end_val {
                continue;
            }
            if cds_end > thick_end_val {
                break;
            }

            let mut thick_start = thick_start_val as i32;
            let mut thick_end = thick_end_val as i32;

            if thick_start < cds_start as i32 {
                thick_start = cds_start as i32;
            }
            if thick_end > cds_end as i32 {
                thick_end = cds_end as i32;
            }

            thick_start = thick_start.saturating_sub(config.extension);
            thick_end = thick_end.saturating_add(config.extension);

            if active_zero_based && thick_start < thick_end {
                thick_start = thick_start.saturating_add(1);
            }

            if thick_start > thick_end {
                continue;
            }

            let final_start = std::cmp::max(1, thick_start) as u32;
            let final_end = std::cmp::max(1, thick_end) as u32;

            thick_regions.push(Region {
                chr: chr.clone(),
                start: Pos1(final_start),
                end: Pos1(final_end),
                gene: gene.clone(),
                amplicon: None,
            });
        }
        segs.push(thick_regions);
    }

    Ok(segs)
}

/// Parses BED file in amplicon mode.
pub fn parse_amp_bed_from_reader<R: BufRead>(
    reader: R,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    zero_based_override: Option<bool>,
) -> Result<Vec<Vec<Region>>, VardictError> {
    let is_zero_based = config.zero_based;
    let mut tsegs: indexmap::IndexMap<String, Vec<Region>> = indexmap::IndexMap::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.trim().is_empty()
            || line.starts_with('#')
            || line.starts_with("browser")
            || line.starts_with("track")
        {
            continue;
        }

        let cols = java_split(&line, &config.delimiter)?;

        if cols.len() <= 7 {
            return Err(VardictError::Parse(format!(
                "Amplicon BED line has fewer than 8 columns (got {}): '{}'",
                cols.len(),
                line.trim()
            )));
        }

        let chr = correct_chromosome(chromosomes_lengths, cols[AMP_BED_ROW_FORMAT.chr_column]);
        let mut start = parse_int_with_commas(cols[AMP_BED_ROW_FORMAT.start_column])? as i32;
        let end = parse_int_with_commas(cols[AMP_BED_ROW_FORMAT.end_column])? as i32;
        let gene = cols[AMP_BED_ROW_FORMAT.gene_column].to_string();
        let mut insert_start =
            parse_int_with_commas(cols[AMP_BED_ROW_FORMAT.thick_start_column])? as i32;
        let insert_end = parse_int_with_commas(cols[AMP_BED_ROW_FORMAT.thick_end_column])? as i32;

        let active_zero_based = zero_based_override.unwrap_or(is_zero_based);

        if active_zero_based && start < end {
            start = start.saturating_add(1);
            insert_start = insert_start.saturating_add(1);
        }

        let final_start = std::cmp::max(1, start) as u32;
        let final_end = std::cmp::max(1, end) as u32;
        let final_insert_start = std::cmp::max(1, insert_start) as u32;
        let final_insert_end = std::cmp::max(1, insert_end) as u32;

        let reg = Region {
            chr: chr.clone(),
            start: Pos1(final_start),
            end: Pos1(final_end),
            gene,
            amplicon: Some(AmpliconInfo {
                insert_start: Pos1(final_insert_start),
                insert_end: Pos1(final_insert_end),
            }),
        };

        tsegs.entry(chr).or_default().push(reg);
    }

    let mut segs = Vec::new();
    if tsegs.is_empty() {
        return Ok(segs);
    }

    let mut current_regions = Vec::new();
    let mut previous_chr = String::new();
    let mut previous_end: i32 = -1;

    let keys: Vec<String> = tsegs.keys().cloned().collect();
    for chr_key in crate::javamap::java_hashmap_key_order(&keys) {
        let Some(chr_regions) = tsegs.swap_remove(chr_key) else {
            continue;
        };
        // Parser invariant: every region pushed above carries amplicon info.
        // §9: propagate as an internal error instead of unwrapping; the
        // sort key is extracted once so the stable order is unchanged.
        let mut keyed = chr_regions
            .into_iter()
            .map(|r| {
                let amp = r.amplicon.as_ref().ok_or_else(|| {
                    VardictError::Other(format!(
                        "internal error: amplicon region {}:{} lacks amplicon info",
                        r.chr, r.start.0
                    ))
                })?;
                Ok((amp.insert_start.0 as i32, amp.insert_end.0 as i32, r))
            })
            .collect::<Result<Vec<_>, VardictError>>()?;

        // Stable sort by insert_start ascending
        keyed.sort_by_key(|(insert_start, _, _)| *insert_start);

        for (insert_start, insert_end, region) in keyed {
            if previous_end != -1 && (region.chr != previous_chr || insert_start > previous_end) {
                segs.push(current_regions);
                current_regions = Vec::new();
            }

            previous_chr = region.chr.clone();
            previous_end = insert_end;
            current_regions.push(region);
        }
    }

    if !current_regions.is_empty() {
        segs.push(current_regions);
    }

    Ok(segs)
}

/// Parses BED file from reader.
pub fn parse_bed_from_reader<R: BufRead>(
    reader: R,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    zero_based_override: Option<bool>,
) -> Result<Vec<Vec<Region>>, VardictError> {
    if config.amplicon.is_some() {
        parse_amp_bed_from_reader(reader, config, chromosomes_lengths, zero_based_override)
    } else {
        parse_non_amp_bed_from_reader(reader, config, chromosomes_lengths, zero_based_override)
    }
}

/// Parses BED file from path.
pub fn parse_bed<P: AsRef<Path>>(
    path: P,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    zero_based_override: Option<bool>,
) -> Result<Vec<Vec<Region>>, VardictError> {
    let file = File::open(path)?;
    parse_bed_from_reader(
        BufReader::new(file),
        config,
        chromosomes_lengths,
        zero_based_override,
    )
}

fn chromosome_in_header(chromosomes_lengths: &HashMap<String, usize>, chr: &str) -> bool {
    if chromosomes_lengths.contains_key(chr) {
        return true;
    }
    if let Some(stripped) = chr.strip_prefix("chr") {
        chromosomes_lengths.contains_key(stripped)
    } else {
        let mut toggled = String::with_capacity(3 + chr.len());
        toggled.push_str("chr");
        toggled.push_str(chr);
        chromosomes_lengths.contains_key(&toggled)
    }
}

fn is_numeric_coord_token(s: &str) -> bool {
    let mut parts = s.split('-');
    let Some(p1) = parts.next() else {
        return false;
    };
    let is_digits_or_commas = |p: &str| {
        let t = p.trim();
        !t.is_empty()
            && t.bytes().all(|b| b.is_ascii_digit() || b == b',')
            && t.bytes().any(|b| b.is_ascii_digit())
    };
    if !is_digits_or_commas(p1) {
        return false;
    }
    match parts.next() {
        None => true,
        Some(p2) => is_digits_or_commas(p2) && parts.next().is_none(),
    }
}

fn parse_coord_segment(range_str: &str) -> Result<(i32, i32), VardictError> {
    let mut parts = range_str.split('-');
    let Some(start_str) = parts.next() else {
        return Err(VardictError::Parse(format!(
            "Invalid coordinate range '{}'",
            range_str
        )));
    };
    let start = parse_int_with_commas(start_str)? as i32;
    let end = match parts.next() {
        None => start,
        Some(end_str) => {
            if parts.next().is_some() {
                return Err(VardictError::Parse(format!(
                    "Invalid coordinate range '{}'",
                    range_str
                )));
            }
            parse_int_with_commas(end_str)? as i32
        }
    };
    Ok((start, end))
}

fn finalize_single_region(
    chr: String,
    mut start: i32,
    end: i32,
    gene: String,
    config: &Config,
) -> Vec<Vec<Region>> {
    start = start.saturating_sub(config.extension);
    let end = end.saturating_add(config.extension);

    let active_zero_based = config.zero_based;
    if active_zero_based && start < end {
        start = start.saturating_add(1);
    }

    if start > end {
        start = end;
    }

    let final_start = std::cmp::max(1, start) as u32;
    let final_end = std::cmp::max(1, end) as u32;

    vec![vec![Region {
        chr,
        start: Pos1(final_start),
        end: Pos1(final_end),
        gene,
        amplicon: None,
    }]]
}

/// Parses a region string (from option -R) into a single nested Region.
pub fn parse_region_flag(
    region_str: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<Vec<Vec<Region>>, VardictError> {
    let split: Vec<&str> = region_str.split(':').collect();
    let Some(&first_seg) = split.first() else {
        return Err(VardictError::Parse(format!(
            "Empty region string: '{}'",
            region_str
        )));
    };

    if split.len() < 2 {
        return Err(VardictError::Parse(format!(
            "Invalid region format (missing coordinates): '{}'",
            region_str
        )));
    }

    if split.len() > 2 && !chromosome_in_header(chromosomes_lengths, first_seg) {
        if !chromosomes_lengths.is_empty() {
            let mut matched_contig_err: Option<VardictError> = None;
            for k in (1..split.len()).rev() {
                let Some(prefix) = split.get(..k) else {
                    continue;
                };
                let candidate = prefix.join(":");
                if chromosome_in_header(chromosomes_lengths, &candidate) {
                    let Some(&coord_str) = split.get(k) else {
                        continue;
                    };
                    match parse_coord_segment(coord_str) {
                        Ok((start, end)) => {
                            let chr = correct_chromosome(chromosomes_lengths, &candidate);
                            let gene = match split.get(k + 1..) {
                                Some(tail) if !tail.is_empty() => tail.join(":"),
                                _ => chr.clone(),
                            };
                            return Ok(finalize_single_region(chr, start, end, gene, config));
                        }
                        Err(err) => {
                            if is_numeric_coord_token(coord_str) {
                                return Err(err);
                            }
                            if matched_contig_err.is_none() {
                                matched_contig_err = Some(err);
                            }
                        }
                    }
                }
            }
            if let Some(err) = matched_contig_err {
                return Err(err);
            }
        } else {
            let min_k = split.len().saturating_sub(2).max(1);
            let mut last_coord_err: Option<VardictError> = None;
            for k in (min_k..split.len()).rev() {
                let Some(&coord_str) = split.get(k) else {
                    continue;
                };
                if !coord_str.contains('-')
                    && k > min_k
                    && split
                        .get(k - 1)
                        .is_some_and(|p| p.contains('-') && parse_coord_segment(p).is_ok())
                {
                    continue;
                }
                match parse_coord_segment(coord_str) {
                    Ok((start, end)) => {
                        let Some(prefix) = split.get(..k) else {
                            continue;
                        };
                        let chr = if k == 1 {
                            correct_chromosome(chromosomes_lengths, first_seg)
                        } else {
                            prefix.join(":")
                        };
                        let gene = match split.get(k + 1) {
                            Some(&g) => g.to_string(),
                            None => chr.clone(),
                        };
                        return Ok(finalize_single_region(chr, start, end, gene, config));
                    }
                    Err(err) => {
                        if is_numeric_coord_token(coord_str) {
                            return Err(err);
                        }
                        last_coord_err = Some(err);
                    }
                }
            }
            if let Some(err) = last_coord_err {
                return Err(err);
            }
        }
    }

    let chr = correct_chromosome(chromosomes_lengths, first_seg);
    let Some(&range_str) = split.get(1) else {
        return Err(VardictError::Parse(format!(
            "Invalid region format (missing coordinates): '{}'",
            region_str
        )));
    };
    let gene = match split.get(2) {
        Some(&g) => g.to_string(),
        None => chr.clone(),
    };

    let range: Vec<&str> = range_str.split('-').collect();
    let Some(&start_str) = range.first() else {
        return Err(VardictError::Parse(format!(
            "Invalid region format (missing coordinates): '{}'",
            region_str
        )));
    };
    let start = parse_int_with_commas(start_str)? as i32;
    let end = match range.get(1) {
        Some(&end_str) => parse_int_with_commas(end_str)? as i32,
        None => start,
    };

    Ok(finalize_single_region(chr, start, end, gene, config))
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn test_correct_chromosome() {
        let mut map = HashMap::new();
        map.insert("chr1".to_string(), 1000);
        map.insert("2".to_string(), 2000);

        assert_eq!(correct_chromosome(&map, "chr1"), "chr1");
        assert_eq!(correct_chromosome(&map, "1"), "chr1");
        assert_eq!(correct_chromosome(&map, "2"), "2");
        assert_eq!(correct_chromosome(&map, "chr2"), "2");
        assert_eq!(correct_chromosome(&map, "chr3"), "3");
    }

    #[test]
    fn test_non_amp_to_regions_raws_1() {
        let mut config = Config::default();
        config.col_chr = 0;
        config.col_start = 1;
        config.col_end = 2;
        config.col_seg_start = 1;
        config.col_seg_end = 2;
        config.col_gene = Some(3);
        config.zero_based = true;
        config.column_for_chromosome = Some(0);
        config.extension = 0;
        config.delimiter = "\t".to_string();

        let raws = "\
Y\t59023460\t59028460\tY_59023460_1
Y\t59028310\t59033310\tY_59028310_2
Y\t59033160\t59038160\tY_59033160_3
";
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("Y".to_string(), 100000000);

        let parsed =
            parse_bed_from_reader(raws.as_bytes(), &config, &chr_lengths, Some(true)).unwrap();

        assert_eq!(parsed.len(), 3);
        assert_eq!(
            parsed[0][0],
            Region {
                chr: "Y".to_string(),
                start: Pos1(59023461),
                end: Pos1(59028460),
                gene: "Y_59023460_1".to_string(),
                amplicon: None
            }
        );
        assert_eq!(
            parsed[1][0],
            Region {
                chr: "Y".to_string(),
                start: Pos1(59028311),
                end: Pos1(59033310),
                gene: "Y_59028310_2".to_string(),
                amplicon: None
            }
        );
    }

    #[test]
    fn test_non_amp_to_regions_raws_2() {
        let mut config = Config::default();
        config.col_chr = 0;
        config.col_start = 1;
        config.col_end = 2;
        config.col_seg_start = 1;
        config.col_seg_end = 2;
        config.col_gene = Some(3);
        config.zero_based = true;
        config.column_for_chromosome = None;
        config.extension = 0;
        config.delimiter = "\t".to_string();

        let raws = "\
Y\t59023460\t59028460\t59028460
";
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("Y".to_string(), 100000000);

        let parsed =
            parse_bed_from_reader(raws.as_bytes(), &config, &chr_lengths, Some(true)).unwrap();

        assert_eq!(parsed.len(), 1);
        assert!(parsed[0].is_empty());
    }

    #[test]
    fn test_amp() {
        // Replicates amp test case in RegionBuilderTest.java
        let mut config = Config::default();
        config.amplicon = Some("10:0.95".to_string());
        config.delimiter = "\t".to_string();
        config.zero_based = true;

        let raws = "\
Y\t59023460\t59028460\tY_59023460_1\t0\t.\t59023460\t59028460
Y\t59028310\t59033310\tY_59028310_2\t0\t.\t59028310\t59033310
Y\t59033160\t59038160\tY_59033160_3\t0\t.\t59033160\t59038160
Y\t59038010\t59043010\tY_59038010_4\t0\t.\t59038010\t59043010
Y\t59042860\t59047860\tY_59042860_5\t0\t.\t59042860\t59047860
";
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("Y".to_string(), 100000000);

        let parsed =
            parse_bed_from_reader(raws.as_bytes(), &config, &chr_lengths, Some(true)).unwrap();

        // 5 overlapping amplicons are grouped into 1 cohort list
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].len(), 5);

        assert_eq!(
            parsed[0][0],
            Region {
                chr: "Y".to_string(),
                start: Pos1(59023461),
                end: Pos1(59028460),
                gene: "Y_59023460_1".to_string(),
                amplicon: Some(AmpliconInfo {
                    insert_start: Pos1(59023461),
                    insert_end: Pos1(59028460),
                }),
            }
        );
        assert_eq!(
            parsed[0][4],
            Region {
                chr: "Y".to_string(),
                start: Pos1(59042861),
                end: Pos1(59047860),
                gene: "Y_59042860_5".to_string(),
                amplicon: Some(AmpliconInfo {
                    insert_start: Pos1(59042861),
                    insert_end: Pos1(59047860),
                }),
            }
        );
    }

    #[test]
    fn test_single_region_flag() {
        // Replicates singleRegion test case in RegionBuilderTest.java
        let mut config = Config::default();
        config.zero_based = true;
        config.extension = 0;

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("Y".to_string(), 100000000);

        let parsed =
            parse_region_flag("Y:59023460-59028460:Y_59023460_1", &config, &chr_lengths).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].len(), 1);
        assert_eq!(
            parsed[0][0],
            Region {
                chr: "Y".to_string(),
                start: Pos1(59023461),
                end: Pos1(59028460),
                gene: "Y_59023460_1".to_string(),
                amplicon: None,
            }
        );
    }

    #[test]
    fn test_tiling_covers_and_overlaps() {
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(1000),
            gene: "GENE".to_string(),
            amplicon: None,
        };

        let tiles = tile(&region, 300, 50);

        // Expecting 4 tiles: 100-399, 400-699, 700-999, 1000-1000
        assert_eq!(tiles.len(), 4);

        assert_eq!(tiles[0].start, Pos1(100));
        assert_eq!(tiles[0].end, Pos1(399));
        assert_eq!(tiles[0].halo_start, Pos1(100)); // capped at reg_start (100)
        assert_eq!(tiles[0].halo_end, Pos1(449)); // 399 + 50

        assert_eq!(tiles[1].start, Pos1(400));
        assert_eq!(tiles[1].end, Pos1(699));
        assert_eq!(tiles[1].halo_start, Pos1(350)); // 400 - 50
        assert_eq!(tiles[1].halo_end, Pos1(749)); // 699 + 50

        assert_eq!(tiles[2].start, Pos1(700));
        assert_eq!(tiles[2].end, Pos1(999));
        assert_eq!(tiles[2].halo_start, Pos1(650)); // 700 - 50
        assert_eq!(tiles[2].halo_end, Pos1(1000)); // capped at reg_end (1000)

        assert_eq!(tiles[3].start, Pos1(1000));
        assert_eq!(tiles[3].end, Pos1(1000));
        assert_eq!(tiles[3].halo_start, Pos1(950)); // 1000 - 50
        assert_eq!(tiles[3].halo_end, Pos1(1000)); // capped at reg_end (1000)

        // Check active union covers the entire 100..1000 region
        let mut curr_expected = 100;
        for t in &tiles {
            assert_eq!(t.start.0, curr_expected);
            curr_expected = t.end.0 + 1;
        }
        assert_eq!(curr_expected, 1001); // Covers up to 1000

        // Check neighboring overlap is 2 * halo (except boundary-limited ones)
        // Tile 0 and 1: tile 1 halo_start is 350, tile 0 halo_end is 449. Overlap is 350..449 (len = 100)
        let overlap_len = tiles[0].halo_end.0 - tiles[1].halo_start.0 + 1;
        assert_eq!(overlap_len, 100);
    }

    #[test]
    fn test_tile_at_u32_max_terminates() {
        // One-base region at the u32 ceiling must yield exactly one tile and
        // terminate (saturating add-then-sub stalls one base short and spins).
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(u32::MAX),
            end: Pos1(u32::MAX),
            gene: "GENE".to_string(),
            amplicon: None,
        };
        let tiles = tile(&region, 1, 50);
        assert_eq!(tiles.len(), 1);
        assert_eq!(tiles[0].start, Pos1(u32::MAX));
        assert_eq!(tiles[0].end, Pos1(u32::MAX));
    }

    #[test]
    fn test_batch_tiles_preserves_wgs_and_small_beds_while_batching_exome() {
        // 1. Small BED (20 tiles on 8 threads <= 8 * 4 = 32): effective max is 1 (no batching)
        let small_tiles: Vec<Tile> = (0..20u32)
            .map(|i| Tile {
                chr: "chr1".to_string(),
                start: Pos1(i * 1000 + 1),
                end: Pos1(i * 1000 + 200),
                gene: "G".to_string(),
                halo_start: Pos1(i * 1000 + 1),
                halo_end: Pos1(i * 1000 + 200),
                amplicon: None,
                original_start: Pos1(i * 1000 + 1),
                original_end: Pos1(i * 1000 + 200),
            })
            .collect();
        let small_batches = batch_tiles(small_tiles, 100_000, 16, 8);
        assert_eq!(small_batches.len(), 20);
        assert!(small_batches.iter().all(|b| b.len() == 1));

        // 2. Large WGS 100 kbp tiles (1000 tiles): each tile span is 100,000 bp, so each stays batch size 1
        let wgs_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100_000_000),
            gene: "WGS".to_string(),
            amplicon: None,
        };
        let wgs_tiles = tile(&wgs_region, 100_000, 500);
        assert_eq!(wgs_tiles.len(), 1000);
        let wgs_batches = batch_tiles(wgs_tiles, 100_000, 16, 8);
        assert_eq!(wgs_batches.len(), 1000);
        assert!(wgs_batches.iter().all(|b| b.len() == 1));

        // 3. Exome BED (320 exons of 200 bp across chr1 and chr2 on 4 threads):
        // effective_max = (320 / 16).clamp(1, 16) = 16; batches break on chromosome change.
        let mut exome_tiles = Vec::with_capacity(320);
        for i in 0..320u32 {
            let chr = if i < 165 { "chr1" } else { "chr2" };
            let start = (i % 165) * 500 + 1;
            let end = start + 199;
            exome_tiles.push(Tile {
                chr: chr.to_string(),
                start: Pos1(start),
                end: Pos1(end),
                gene: format!("G{i}"),
                halo_start: Pos1(start),
                halo_end: Pos1(end),
                amplicon: None,
                original_start: Pos1(start),
                original_end: Pos1(end),
            });
        }
        let exome_batches = batch_tiles(exome_tiles.clone(), 100_000, 16, 4);
        // 165 tiles on chr1 -> 10 batches of 16 + 1 batch of 5 = 11 batches
        // 155 tiles on chr2 -> 9 batches of 16 + 1 batch of 11 = 10 batches
        assert_eq!(exome_batches.len(), 21);
        assert_eq!(exome_batches[10].len(), 5);
        assert!(exome_batches[10].iter().all(|t| t.chr == "chr1"));
        assert!(exome_batches[11].iter().all(|t| t.chr == "chr2"));
        let flattened: Vec<String> = exome_batches
            .into_iter()
            .flatten()
            .map(|t| t.gene)
            .collect();
        let expected: Vec<String> = exome_tiles.into_iter().map(|t| t.gene).collect();
        assert_eq!(flattened, expected);
    }

    /// Task 4 pin: tile active intervals partition the region with no gap
    /// and no overlap (each base owned by exactly one tile, so boundary
    /// variants can never duplicate or drop), halos extend `halo` each side
    /// clamped to the region, and split sub-tiles are never merged by
    /// batching (never merge split sub-tiles across halo).
    #[test]
    fn tile_boundary_partition_dedup_pin() {
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(1_000_000),
            gene: "WGS".to_string(),
            amplicon: None,
        };
        let halo = 500u32;
        let tiles = tile(&region, 100_000, halo);
        assert_eq!(tiles.len(), 10);

        // Active intervals: exact partition, single ownership per base.
        assert_eq!(tiles[0].start, Pos1(1));
        assert_eq!(tiles[tiles.len() - 1].end, Pos1(1_000_000));
        for pair in tiles.windows(2) {
            assert_eq!(
                pair[0].end.0 + 1,
                pair[1].start.0,
                "adjacent active intervals must abut exactly"
            );
            // Interior halos overlap by exactly 2 * halo.
            assert_eq!(
                pair[0].halo_end.0 - pair[1].halo_start.0 + 1,
                2 * halo,
                "adjacent halos must overlap by exactly 2 * halo"
            );
        }

        // Halo clamping at region bounds.
        assert_eq!(tiles[0].halo_start, Pos1(1));
        assert_eq!(tiles[0].halo_end, Pos1(100_000 + halo));
        let last = tiles.len() - 1;
        assert_eq!(tiles[last].halo_start, Pos1(900_001 - halo));
        assert_eq!(tiles[last].halo_end, Pos1(1_000_000));
        // Interior tile halos extend exactly `halo` each side.
        assert_eq!(tiles[1].halo_start, Pos1(100_001 - halo));
        assert_eq!(tiles[1].halo_end, Pos1(200_000 + halo));

        // Split sub-tiles of one region are never merged across the halo.
        let batches = batch_tiles(tiles.clone(), 100_000, 16, 8);
        assert_eq!(batches.len(), 10);
        assert!(batches.iter().all(|b| b.len() == 1));
        let flattened: Vec<(u32, u32)> = batches
            .into_iter()
            .flatten()
            .map(|t| (t.start.0, t.end.0))
            .collect();
        let expected: Vec<(u32, u32)> = tiles.iter().map(|t| (t.start.0, t.end.0)).collect();
        assert_eq!(flattened, expected);
    }

    #[test]
    fn detect_amplicon_bed_truth_table() {
        // 8-col valid amplicon (thick inside outer) -> true
        assert!(detect_amplicon_bed_from_reader(
            "chr1\t100\t200\tGENE\t0\t+\t120\t180\n".as_bytes(),
            "\t"
        )
        .unwrap());
        // Header lines skipped before first data line -> true
        assert!(detect_amplicon_bed_from_reader(
            "#comment\nbrowser pos\ntrack name=foo\nchr1\t100\t200\tGENE\t0\t+\t120\t180\n"
                .as_bytes(),
            "\t"
        )
        .unwrap());
        // 8-col where thick extends outside outer (s2 < s1 or e2 > e1) -> false
        assert!(!detect_amplicon_bed_from_reader(
            "chr1\t100\t200\tGENE\t0\t+\t90\t180\n".as_bytes(),
            "\t"
        )
        .unwrap());
        // 8-col non-numeric col 6 -> false
        assert!(!detect_amplicon_bed_from_reader(
            "chr1\t100\t200\tGENE\t0\t+\tfoo\t180\n".as_bytes(),
            "\t"
        )
        .unwrap());
        // 4-col standard BED -> false
        assert!(
            !detect_amplicon_bed_from_reader("chr1\t100\t200\tGENE\n".as_bytes(), "\t").unwrap()
        );
        // 9-col row followed by a valid 8-col amplicon row -> true (scans all rows)
        assert!(detect_amplicon_bed_from_reader(
            "chr1\t100\t200\tGENE\t0\t+\t120\t180\textra\nchr1\t300\t400\tGENE2\t0\t+\t320\t380\n"
                .as_bytes(),
            "\t"
        )
        .unwrap());
    }

    #[test]
    fn bed_reader_java_semantics() {
        let mut chrs = HashMap::new();
        chrs.insert("chr1".to_string(), 10000usize);
        let cfg = Config::default();

        // 1a. 4-col auto-detect (c4.bed: RegionBuilder.java:16, 48-68) uses CUSTOM_BED_ROW_FORMAT
        //     where gene_column = 2 ("101"), thick_start_column = 3 ("5"), thick_end_column = 1 ("101").
        let c4_bed = "chr1\t101\t101\t5\n";
        let c4_segs =
            parse_non_amp_bed_from_reader(c4_bed.as_bytes(), &cfg, &chrs, Some(false)).unwrap();
        assert_eq!(c4_segs.len(), 1);
        assert_eq!(c4_segs[0].len(), 1);
        assert_eq!(c4_segs[0][0].gene, "101");
        assert_eq!(c4_segs[0][0].start, Pos1(101));
        assert_eq!(c4_segs[0][0].end, Pos1(101));

        // 1b. Explicit column flags (-c 1 -S 2 -E 3 -g 4): preserves trailing space in gene name,
        //     strips trailing empty tabs via java_split, and skips browser/track header lines.
        let mut trail_cfg = Config::default();
        trail_cfg.column_for_chromosome = Some(0);
        trail_cfg.col_chr = 0;
        trail_cfg.col_start = 1;
        trail_cfg.col_end = 2;
        trail_cfg.col_seg_start = 1;
        trail_cfg.col_seg_end = 2;
        trail_cfg.col_gene = Some(3);
        trail_cfg.zero_based = true;
        let bed = "browser position chr1:1-1000\ntrack name=test\nchr1\t100\t200\tGENE \t\t\n";
        let segs = parse_non_amp_bed_from_reader(bed.as_bytes(), &trail_cfg, &chrs, None).unwrap();
        assert_eq!(segs.len(), 1);
        assert_eq!(segs[0].len(), 1);
        assert_eq!(segs[0][0].gene, "GENE ");
        assert_eq!(segs[0][0].start, Pos1(101));
        assert_eq!(segs[0][0].end, Pos1(200));

        // 2. Regex delimiter `\s+` splits runs of whitespace
        let mut regex_cfg = trail_cfg.clone();
        regex_cfg.delimiter = r"\s+".to_string();
        let ws_bed = "chr1   100   200   TP53\n";
        let ws_segs =
            parse_non_amp_bed_from_reader(ws_bed.as_bytes(), &regex_cfg, &chrs, None).unwrap();
        assert_eq!(ws_segs[0][0].gene, "TP53");

        // 3. Inverted region (start > end under 1-based) is skipped rather than clamped
        let mut inv_cfg = Config::default();
        inv_cfg.column_for_chromosome = Some(0);
        inv_cfg.col_chr = 0;
        inv_cfg.col_start = 1;
        inv_cfg.col_end = 2;
        inv_cfg.col_seg_start = 1;
        inv_cfg.col_seg_end = 2;
        inv_cfg.zero_based = false;
        let inv_bed = "chr1\t200\t100\tGENE\n";
        let inv_segs =
            parse_non_amp_bed_from_reader(inv_bed.as_bytes(), &inv_cfg, &chrs, Some(false))
                .unwrap();
        assert_eq!(inv_segs.len(), 1);
        assert!(inv_segs[0].is_empty(), "inverted region must be skipped");

        // 4. 4-col BED with non-numeric col 1/2 errors loudly (RegionBuilderTest.testIncorrectRegions)
        let bad_4col = "chr1\tbad\t200\tGENE\n";
        assert!(parse_non_amp_bed_from_reader(bad_4col.as_bytes(), &cfg, &chrs, None).is_err());
    }

    #[test]
    fn parse_int_with_commas_i32_bound() {
        assert_eq!(parse_int_with_commas("2147483647").unwrap(), 2_147_483_647);
        assert_eq!(
            parse_int_with_commas("2,147,483,647").unwrap(),
            2_147_483_647
        );
        for invalid in ["2147483648", "2,147,483,648", "4,294,967,295", "-5"] {
            assert!(
                parse_int_with_commas(invalid).is_err(),
                "expected Err for '{invalid}'"
            );
        }
    }

    #[test]
    fn parse_region_flag_rejects_wrapping_coords() {
        let mut config = Config::default();
        config.zero_based = false;
        config.extension = 0;
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 1000);

        for invalid in [
            "chr1:2147483648-2147483700",
            "chr1:2147483600-2147483700",
            "chr1:1-3000000000",
        ] {
            let err = parse_region_flag(invalid, &config, &chr_lengths)
                .expect_err("coordinate > i32::MAX must fail");
            assert!(
                err.to_string().contains("invalid integer"),
                "unexpected error: {err}"
            );
        }

        let max_ok =
            parse_region_flag("chr1:2147483647-2147483647", &config, &chr_lengths).unwrap();
        assert_eq!(max_ok[0][0].start, Pos1(2_147_483_647));
        assert_eq!(max_ok[0][0].end, Pos1(2_147_483_647));
    }

    #[test]
    fn bed_rejects_coordinates_above_i32_max() {
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 1000);

        // 3-column BED
        let mut cfg3 = Config::default();
        cfg3.col_chr = 0;
        cfg3.col_start = 1;
        cfg3.col_end = 2;
        cfg3.col_seg_start = 1;
        cfg3.col_seg_end = 2;
        cfg3.delimiter = "\t".to_string();
        assert!(parse_bed_from_reader(
            "chr1\t2147483648\t2147483700\n".as_bytes(),
            &cfg3,
            &chr_lengths,
            Some(true),
        )
        .is_err());

        // 4-column BED
        let mut cfg4 = Config::default();
        cfg4.delimiter = "\t".to_string();
        assert!(parse_bed_from_reader(
            "chr1\t100\t2147483648\tg1\n".as_bytes(),
            &cfg4,
            &chr_lengths,
            Some(true),
        )
        .is_err());

        // 8-column non-amplicon BED with cdsStart / thickEnd = 3,000,000,000
        let mut cfg8 = Config::default();
        cfg8.col_chr = 0;
        cfg8.col_start = 1;
        cfg8.col_end = 2;
        cfg8.col_seg_start = 6;
        cfg8.col_seg_end = 7;
        cfg8.col_gene = Some(3);
        cfg8.delimiter = "\t".to_string();
        assert!(parse_bed_from_reader(
            "chr1\t3000000000\t3000000100\tg1\t0\t.\t100\t200\n".as_bytes(),
            &cfg8,
            &chr_lengths,
            Some(true),
        )
        .is_err());
        assert!(parse_bed_from_reader(
            "chr1\t100\t200\tg1\t0\t.\t100\t3000000000\n".as_bytes(),
            &cfg8,
            &chr_lengths,
            Some(true),
        )
        .is_err());

        // Amplicon BED with insert_start = 2147483648
        let mut cfg_amp = Config::default();
        cfg_amp.amplicon = Some("10:0.95".to_string());
        cfg_amp.delimiter = "\t".to_string();
        assert!(parse_bed_from_reader(
            "chr1\t100\t200\tg1\t0\t.\t2147483648\t200\n".as_bytes(),
            &cfg_amp,
            &chr_lengths,
            Some(true),
        )
        .is_err());
    }

    #[test]
    fn extension_overflow_saturates_and_amplicon_z_wraps() {
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 1000);

        // -x 2147483647 -R chr1:1000-2000 saturates at Pos1(1) / Pos1(i32::MAX as u32) instead of wrapping.
        let mut cfg = Config::default();
        cfg.extension = i32::MAX;
        cfg.zero_based = false;
        let parsed = parse_region_flag("chr1:1000-2000", &cfg, &chr_lengths).unwrap();
        assert_eq!(parsed[0][0].start, Pos1(1));
        assert_eq!(parsed[0][0].end, Pos1(i32::MAX as u32));

        // Non-amplicon BED with -x 2147483647 also saturates at Pos1(1) / Pos1(i32::MAX as u32).
        let mut cfg_bed = Config::default();
        cfg_bed.col_chr = 0;
        cfg_bed.col_start = 1;
        cfg_bed.col_end = 2;
        cfg_bed.col_seg_start = 1;
        cfg_bed.col_seg_end = 2;
        cfg_bed.col_gene = Some(3);
        cfg_bed.column_for_chromosome = Some(0);
        cfg_bed.extension = i32::MAX;
        cfg_bed.delimiter = "\t".to_string();
        let bed_parsed = parse_non_amp_bed_from_reader(
            "chr1\t1000\t2000\tg1\n".as_bytes(),
            &cfg_bed,
            &chr_lengths,
            Some(true),
        )
        .unwrap();
        assert_eq!(bed_parsed[0][0].start, Pos1(1));
        assert_eq!(bed_parsed[0][0].end, Pos1(i32::MAX as u32));

        // Non-amplicon 13-column refGene row with col_gene = None (-g 0) disables gene lookup
        // (falling back to chr) rather than reading default column 13 (index 12).
        let mut cfg_no_gene = Config::default();
        cfg_no_gene.col_gene = None;
        cfg_no_gene.delimiter = "\t".to_string();
        let refgene_13 = "bin\tNM_1\tchr1\t+\t90\t210\t100\t200\t1\t100,\t200,\t0\tTP53\n";
        let no_gene_parsed =
            parse_non_amp_bed_from_reader(refgene_13.as_bytes(), &cfg_no_gene, &chr_lengths, None)
                .unwrap();
        assert_eq!(no_gene_parsed[0][0].gene, "chr1");

        // Amplicon row with insert_start = 2147483647 under -z saturates insert_start + 1 without wrapping.
        let mut cfg_amp = Config::default();
        cfg_amp.amplicon = Some("10:0.95".to_string());
        cfg_amp.delimiter = "\t".to_string();
        cfg_amp.zero_based = true;
        let amp_parsed = parse_bed_from_reader(
            "chr1\t100\t200\tg1\t0\t.\t2147483647\t200\n".as_bytes(),
            &cfg_amp,
            &chr_lengths,
            Some(true),
        )
        .unwrap();
        let amp = amp_parsed[0][0].amplicon.as_ref().unwrap();
        assert_eq!(amp.insert_start, Pos1(i32::MAX as u32));
        assert_eq!(amp.insert_end, Pos1(200));
    }

    #[test]
    fn undersized_5_and_6_column_bed_error_includes_column_flag_hint() {
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 1000);
        let cfg = Config::default();

        for (ncols, line) in [
            (5usize, "chr1\t100\t200\tGENE\t0\n"),
            (6usize, "chr1\t100\t200\tGENE\t0\t+\n"),
        ] {
            let err = parse_non_amp_bed_from_reader(line.as_bytes(), &cfg, &chr_lengths, None)
                .expect_err("5- and 6-column BED without -c/-S/-E must fail");
            let msg = err.to_string();
            assert!(
                msg.contains(&format!(
                    "Region line has {ncols} columns, expected at least"
                )),
                "unexpected error message for {ncols}-col BED: {msg}"
            );
            assert!(
                msg.contains("pass -c 1 -S 2 -E 3 -z for standard 0-based BED"),
                "expected column + -z flag hint in {ncols}-col BED error: {msg}"
            );
        }
    }

    #[test]
    fn parse_region_flag_supports_hla_colon_contigs() {
        let mut config = Config::default();
        config.zero_based = false;
        config.extension = 0;

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("HLA-A*01:01:01:01".to_string(), 3503);
        chr_lengths.insert("HLA-A*01:01".to_string(), 3503);

        // 1. Exact HLA contig without trailing gene (longest prefix wins over shorter prefix)
        let r1 = parse_region_flag("HLA-A*01:01:01:01:100-200", &config, &chr_lengths).unwrap();
        assert_eq!(
            r1[0][0],
            Region {
                chr: "HLA-A*01:01:01:01".to_string(),
                start: Pos1(100),
                end: Pos1(200),
                gene: "HLA-A*01:01:01:01".to_string(),
                amplicon: None,
            }
        );

        // 2. Exact HLA contig with trailing :HLA-A gene and comma-formatted coordinates
        let r2 =
            parse_region_flag("HLA-A*01:01:01:01:100-200:HLA-A", &config, &chr_lengths).unwrap();
        assert_eq!(
            r2[0][0],
            Region {
                chr: "HLA-A*01:01:01:01".to_string(),
                start: Pos1(100),
                end: Pos1(200),
                gene: "HLA-A".to_string(),
                amplicon: None,
            }
        );
        let r2_commas =
            parse_region_flag("HLA-A*01:01:01:01:1,000-2,000:HLA-A", &config, &chr_lengths)
                .unwrap();
        assert_eq!(r2_commas[0][0].start, Pos1(1_000));
        assert_eq!(r2_commas[0][0].end, Pos1(2_000));
        assert_eq!(r2_commas[0][0].gene, "HLA-A");

        // 3. Single-base coordinate HLA-A*01:01:01:01:150
        let r3 = parse_region_flag("HLA-A*01:01:01:01:150", &config, &chr_lengths).unwrap();
        assert_eq!(
            r3[0][0],
            Region {
                chr: "HLA-A*01:01:01:01".to_string(),
                start: Pos1(150),
                end: Pos1(150),
                gene: "HLA-A*01:01:01:01".to_string(),
                amplicon: None,
            }
        );

        // 4. Empty chromosomes_lengths fallback (preserves colon-bearing contig without "chr" prefix)
        let empty = HashMap::new();
        let r_empty = parse_region_flag("HLA-A*01:01:01:01:100-200", &config, &empty).unwrap();
        assert_eq!(r_empty[0][0].chr, "HLA-A*01:01:01:01");
        assert_eq!(r_empty[0][0].start, Pos1(100));
        assert_eq!(r_empty[0][0].end, Pos1(200));
        assert_eq!(r_empty[0][0].gene, "HLA-A*01:01:01:01");

        let r_empty_gene =
            parse_region_flag("HLA-A*01:01:01:01:100-200:HLA-A", &config, &empty).unwrap();
        assert_eq!(r_empty_gene[0][0].chr, "HLA-A*01:01:01:01");
        assert_eq!(r_empty_gene[0][0].start, Pos1(100));
        assert_eq!(r_empty_gene[0][0].end, Pos1(200));
        assert_eq!(r_empty_gene[0][0].gene, "HLA-A");

        // 5. chr-toggled HLA lookup in both directions
        let mut chr_prefixed_header = HashMap::new();
        chr_prefixed_header.insert("chrHLA-A*01:01:01:01".to_string(), 3503);
        let toggled_add = parse_region_flag(
            "HLA-A*01:01:01:01:100-200:HLA-A",
            &config,
            &chr_prefixed_header,
        )
        .unwrap();
        assert_eq!(toggled_add[0][0].chr, "chrHLA-A*01:01:01:01");
        assert_eq!(toggled_add[0][0].gene, "HLA-A");

        let toggled_strip =
            parse_region_flag("chrHLA-A*01:01:01:01:100-200:HLA-A", &config, &chr_lengths).unwrap();
        assert_eq!(toggled_strip[0][0].chr, "HLA-A*01:01:01:01");
        assert_eq!(toggled_strip[0][0].gene, "HLA-A");

        // 6. Out-of-range coordinate on HLA contig errors rather than misparsing allele digits
        assert!(parse_region_flag(
            "HLA-A*01:01:01:01:2147483648-2147483700",
            &config,
            &chr_lengths
        )
        .is_err());
        assert!(
            parse_region_flag("HLA-A*01:01:01:01:2147483648-2147483700", &config, &empty).is_err()
        );

        // 7. All-numeric Entrez gene IDs without header and colon-bearing gene names with header
        let r_entrez_chr = parse_region_flag("chr1:100-200:7157", &config, &empty).unwrap();
        assert_eq!(r_entrez_chr[0][0].chr, "1");
        assert_eq!(r_entrez_chr[0][0].start, Pos1(100));
        assert_eq!(r_entrez_chr[0][0].end, Pos1(200));
        assert_eq!(r_entrez_chr[0][0].gene, "7157");

        let r_entrez_num_chr = parse_region_flag("1:100-200:7157", &config, &empty).unwrap();
        assert_eq!(r_entrez_num_chr[0][0].chr, "chr1");
        assert_eq!(r_entrez_num_chr[0][0].start, Pos1(100));
        assert_eq!(r_entrez_num_chr[0][0].end, Pos1(200));
        assert_eq!(r_entrez_num_chr[0][0].gene, "7157");

        let r_entrez_hla =
            parse_region_flag("HLA-A*01:01:01:01:100-200:3105", &config, &empty).unwrap();
        assert_eq!(r_entrez_hla[0][0].chr, "HLA-A*01:01:01:01");
        assert_eq!(r_entrez_hla[0][0].start, Pos1(100));
        assert_eq!(r_entrez_hla[0][0].end, Pos1(200));
        assert_eq!(r_entrez_hla[0][0].gene, "3105");

        let r_colon_gene = parse_region_flag(
            "HLA-A*01:01:01:01:100-200:HLA-A*01:01",
            &config,
            &chr_lengths,
        )
        .unwrap();
        assert_eq!(r_colon_gene[0][0].chr, "HLA-A*01:01:01:01");
        assert_eq!(r_colon_gene[0][0].start, Pos1(100));
        assert_eq!(r_colon_gene[0][0].end, Pos1(200));
        assert_eq!(r_colon_gene[0][0].gene, "HLA-A*01:01");
    }
}
