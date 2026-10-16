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
        let mut curr_end = curr_start.saturating_add(tile_size).saturating_sub(1);
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

/// Helper to parse a u32 from string with comma removal.
pub fn parse_int_with_commas(s: &str) -> Result<u32, VardictError> {
    let clean = s.replace(",", "");
    clean
        .trim()
        .parse::<u32>()
        .map_err(|e| VardictError::Parse(format!("invalid integer '{}': {}", s, e)))
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
        gene_column: config.col_gene.unwrap_or(12),
    };

    let mut segs = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = trimmed.split(&config.delimiter).collect();

        let mut current_format = format;
        let mut current_zero_based = is_zero_based;

        // Auto-detect 4-column BED format
        if config.column_for_chromosome.is_none() && cols.len() == 4 {
            if let (Ok(a1), Ok(a2)) = (
                parse_int_with_commas(cols[1]),
                parse_int_with_commas(cols[2]),
            ) {
                if a1 <= a2 {
                    current_format = CUSTOM_BED_ROW_FORMAT;
                    // Preserve user's gene column if explicitly set
                    if let Some(user_gene_col) = config.col_gene {
                        current_format.gene_column = user_gene_col;
                    }
                    if zero_based_override.is_none() {
                        current_zero_based = true;
                    }
                }
            }
        }

        let active_zero_based = zero_based_override.unwrap_or(current_zero_based);

        if cols.len() <= current_format.chr_column
            || cols.len() <= current_format.start_column
            || cols.len() <= current_format.end_column
            || cols.len() <= current_format.thick_start_column
            || cols.len() <= current_format.thick_end_column
        {
            return Err(VardictError::Parse(format!(
                "BED line has fewer columns ({}) than required by active format {:?}",
                cols.len(),
                current_format
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

            thick_start -= config.extension;
            thick_end += config.extension;

            if active_zero_based && thick_start < thick_end {
                thick_start += 1;
            }

            if thick_start > thick_end {
                thick_start = thick_end;
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
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = trimmed.split(&config.delimiter).collect();

        if cols.len() <= 7 {
            return Err(VardictError::Parse(format!(
                "Amplicon BED line has fewer than 8 columns (got {}): '{}'",
                cols.len(),
                trimmed
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
            start += 1;
            insert_start += 1;
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

    for (_, chr_regions) in tsegs {
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

/// Parses a region string (from option -R) into a single nested Region.
pub fn parse_region_flag(
    region_str: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<Vec<Vec<Region>>, VardictError> {
    let split: Vec<&str> = region_str.split(':').collect();
    if split.is_empty() {
        return Err(VardictError::Parse(format!(
            "Empty region string: '{}'",
            region_str
        )));
    }
    let chr = correct_chromosome(chromosomes_lengths, split[0]);

    if split.len() < 2 {
        return Err(VardictError::Parse(format!(
            "Invalid region format (missing coordinates): '{}'",
            region_str
        )));
    }

    let gene = if split.len() < 3 {
        chr.clone()
    } else {
        split[2].to_string()
    };

    let range_str = split[1];
    let range: Vec<&str> = range_str.split('-').collect();

    let mut start = parse_int_with_commas(range[0])? as i32;
    let end = if range.len() < 2 {
        start
    } else {
        parse_int_with_commas(range[1])? as i32
    };

    start -= config.extension;
    let end = end + config.extension;

    let active_zero_based = config.zero_based;
    if active_zero_based && start < end {
        start += 1;
    }

    if start > end {
        start = end;
    }

    let final_start = std::cmp::max(1, start) as u32;
    let final_end = std::cmp::max(1, end) as u32;

    Ok(vec![vec![Region {
        chr,
        start: Pos1(final_start),
        end: Pos1(final_end),
        gene,
        amplicon: None,
    }]])
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
}
