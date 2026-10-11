//! Reference FASTA and k-mer seed map loading.
//!
//! Ports Reference.java and ReferenceResource.java.

use crate::config::Config;
use crate::regions::Region as BedRegion;
use crate::types::{Pos1, VardictError};
use noodles::core::Region as NoodlesRegion;
use noodles::fasta::fai;
use noodles::fasta::io::IndexedReader;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

thread_local! {
    static FASTA_READER_CACHE: RefCell<HashMap<PathBuf, IndexedReader<BufReader<File>>>> = RefCell::new(HashMap::new());
}

/// Fetches a subsequence from an indexed FASTA file using a thread-local cached reader.
pub fn fetch_subsequence(
    fasta_path: &Path,
    chr: &str,
    start: Pos1,
    end: Pos1,
) -> Result<Vec<u8>, VardictError> {
    FASTA_READER_CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        // §9: entry API — build the reader only on the vacant path, no
        // just-inserted unwrap on the lookup path.
        let reader = match cache.entry(fasta_path.to_path_buf()) {
            std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
            std::collections::hash_map::Entry::Vacant(entry) => {
                let mut fai_path = fasta_path.to_path_buf();
                let mut ext = fai_path
                    .extension()
                    .map(|e| e.to_string_lossy().into_owned())
                    .unwrap_or_default();
                ext.push_str(".fai");
                fai_path.set_extension(ext);

                let index = fai::fs::read(&fai_path).map_err(VardictError::Io)?;
                let file = File::open(fasta_path).map_err(VardictError::Io)?;
                entry.insert(IndexedReader::new(BufReader::new(file), index))
            }
        };

        // NoodlesRegion parses "chr:start-end"
        let region_str = format!("{}:{}-{}", chr, start.0, end.0);
        let region: NoodlesRegion = region_str.parse().map_err(|e| {
            VardictError::Parse(format!("failed to parse region '{}': {}", region_str, e))
        })?;

        let record = reader.query(&region).map_err(VardictError::Io)?;

        Ok(record.sequence().as_ref().to_vec())
    })
}

/// Represents a loaded contiguous block of the reference sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoadedRegion {
    pub chr: String,
    pub sequence_start: u32,
    pub sequence_end: u32,
}

/// Stores reference bases and k-mer seed maps for sequence realignment.
#[derive(Debug, Clone)]
pub struct Reference {
    pub loaded_regions: Vec<LoadedRegion>,
    /// Maps chromosome name to contiguous bases (upper-case) for the primary loaded region.
    pub reference_sequences: HashMap<String, Vec<u8>>,
    /// Sequences matching loaded_regions by index: region_sequences[i] corresponds to loaded_regions[i].
    pub region_sequences: Vec<Vec<u8>>,
    /// Maps k-mer bytes to 1-based positions.
    pub seed: HashMap<Vec<u8>, Vec<u32>>,
}

impl Reference {
    pub fn new() -> Self {
        Reference {
            loaded_regions: Vec::new(),
            reference_sequences: HashMap::new(),
            region_sequences: Vec::new(),
            seed: HashMap::new(),
        }
    }

    /// Safely gets the base at a 1-based chromosome position.
    pub fn get_base(&self, chr: &str, pos: u32) -> Option<u8> {
        for (idx, lr) in self.loaded_regions.iter().enumerate() {
            if lr.chr == chr && pos >= lr.sequence_start && pos <= lr.sequence_end {
                let offset = (pos - lr.sequence_start) as usize;
                if let Some(bases) = self.region_sequences.get(idx) {
                    if let Some(&b) = bases.get(offset) {
                        return Some(b);
                    }
                }
                if let Some(bases) = self.reference_sequences.get(chr) {
                    if let Some(&b) = bases.get(offset) {
                        return Some(b);
                    }
                }
            }
        }
        None
    }

    /// Populates k-mer mapping positions into the seed map.
    fn add_position_to_seed(&mut self, key_sequence: &[u8], pos: u32) {
        self.seed
            .entry(key_sequence.to_vec())
            .or_default()
            .push(pos);
    }
}

impl Default for Reference {
    fn default() -> Self {
        Self::new()
    }
}

/// Checks if a region is already loaded in Reference.
pub fn is_loaded(chr: &str, sequence_start: u32, sequence_end: u32, reference: &Reference) -> bool {
    for lr in &reference.loaded_regions {
        if lr.chr == chr && sequence_start >= lr.sequence_start && sequence_end <= lr.sequence_end {
            return true;
        }
    }
    false
}

/// Ports ReferenceResource.java:getReference.
pub fn get_reference(
    region: &BedRegion,
    extension: u32,
    chromosomes_lengths: &HashMap<String, usize>,
    config: &Config,
    reference: &mut Reference,
) -> Result<(), VardictError> {
    let start_val = region.start.0;
    let end_val = region.end.0;

    let extend_nuc = config.extension_nuc as u32; // config.numberNucleotideToExtend
    let total_left = extend_nuc + extension;

    let sequence_start = if start_val > total_left {
        start_val - total_left
    } else {
        1
    };

    let len = chromosomes_lengths.get(&region.chr).copied().unwrap_or(0) as u32;
    let total_right = extend_nuc + extension;
    let sequence_end = if len > 0 && end_val + total_right > len {
        len
    } else {
        end_val + total_right
    };

    if is_loaded(&region.chr, sequence_start, sequence_end, reference) {
        return Ok(());
    }

    let fasta_path = Path::new(&config.reference);
    let mut bases = fetch_subsequence(
        fasta_path,
        &region.chr,
        Pos1(sequence_start),
        Pos1(sequence_end),
    )?;

    // Convert bases to uppercase matching exon.toUpperCase()
    for b in &mut bases {
        b.make_ascii_uppercase();
    }

    let loaded_region = LoadedRegion {
        chr: region.chr.clone(),
        sequence_start,
        sequence_end,
    };
    reference.loaded_regions.push(loaded_region);

    // Keep reference_sequences populated for primary chromosome region
    reference
        .reference_sequences
        .entry(region.chr.clone())
        .or_insert_with(|| bases.clone());

    reference.region_sequences.push(bases.clone());

    // k-mer seed indexing logic (ReferenceResource.java:106-125)
    let seed1_len = 17;
    let seed2_len = 12;

    let site_end = if len == sequence_end {
        bases.len()
    } else {
        bases.len().saturating_sub(seed1_len)
    };

    let prev_lrs: Vec<LoadedRegion> =
        reference.loaded_regions[..reference.loaded_regions.len().saturating_sub(1)].to_vec();

    for i in 0..site_end {
        let pos = i as u32 + sequence_start;

        // Don't index positions already covered by a previously loaded region
        if prev_lrs
            .iter()
            .any(|lr| lr.chr == region.chr && pos >= lr.sequence_start && pos <= lr.sequence_end)
        {
            continue;
        }

        // Do not create adaptor sequences for the very end of chromosome
        if len == sequence_end && i > bases.len().saturating_sub(seed1_len) {
            continue;
        }

        if i + seed1_len <= bases.len() {
            let key = &bases[i..i + seed1_len];
            reference.add_position_to_seed(key, pos);
        }

        if i + seed2_len <= bases.len() {
            let key = &bases[i..i + seed2_len];
            reference.add_position_to_seed(key, pos);
        }
    }

    Ok(())
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    #[test]
    fn test_fetch_subsequence_and_get_reference() {
        let temp_dir = std::env::temp_dir();
        let fa_path = temp_dir.join("test_reference.fa");
        let fai_path = temp_dir.join("test_reference.fa.fai");

        // Write a mock FASTA
        let mut fa_file = File::create(&fa_path).unwrap();
        fa_file
            .write_all(b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
            .unwrap();

        // Write a mock FAI
        let mut fai_file = File::create(&fai_path).unwrap();
        fai_file.write_all(b"chr1\t40\t6\t40\t41\n").unwrap();

        // Test fetch_subsequence
        let bases = fetch_subsequence(&fa_path, "chr1", Pos1(1), Pos1(10)).unwrap();
        assert_eq!(bases, b"ACGTACGTAC");

        let bases_mid = fetch_subsequence(&fa_path, "chr1", Pos1(11), Pos1(20)).unwrap();
        assert_eq!(bases_mid, b"GTACGTACGT");

        // Test get_reference
        let bed_region = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(5),
            end: Pos1(25),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut config = Config::default();
        config.reference = fa_path.to_string_lossy().into_owned();
        config.extension_nuc = 2;

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 40);

        let mut reference = Reference::new();
        get_reference(&bed_region, 3, &chr_lengths, &config, &mut reference).unwrap();

        // start=5, extend_nuc=2, extension=3 => sequence_start = 5 - 5 = 1 (since 5 <= 5, clamped at 1)
        // end=25, extend_nuc=2, extension=3 => sequence_end = 25 + 5 = 30
        assert_eq!(reference.loaded_regions.len(), 1);
        assert_eq!(reference.loaded_regions[0].sequence_start, 1);
        assert_eq!(reference.loaded_regions[0].sequence_end, 30);

        // Verify base lookup
        assert_eq!(reference.get_base("chr1", 1), Some(b'A'));
        assert_eq!(reference.get_base("chr1", 2), Some(b'C'));
        assert_eq!(reference.get_base("chr1", 30), Some(b'C'));
        assert_eq!(reference.get_base("chr1", 31), None); // out of bounds

        // Verify seed indexing contains our expected 17-mers and 12-mers
        assert!(!reference.seed.is_empty());

        // Cleanup
        let _ = fs::remove_file(&fa_path);
        let _ = fs::remove_file(&fai_path);
    }

    #[test]
    fn seedmap() {
        let temp_dir = std::env::temp_dir();
        let fa_path = temp_dir.join("test_seedmap.fa");
        let fai_path = temp_dir.join("test_seedmap.fa.fai");

        // Write a mock FASTA with a unique 17-mer seed at the end
        let mut fa_file = File::create(&fa_path).unwrap();
        fa_file
            .write_all(b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNNNNNNNNNNNN\n")
            .unwrap();

        // Write a mock FAI
        let mut fai_file = File::create(&fai_path).unwrap();
        fai_file.write_all(b"chr1\t49\t6\t49\t50\n").unwrap();

        let bed_region = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(49),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut config = Config::default();
        config.reference = fa_path.to_string_lossy().into_owned();
        config.extension_nuc = 0;

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 49);

        let mut reference = Reference::new();
        get_reference(&bed_region, 0, &chr_lengths, &config, &mut reference).unwrap();

        // The seed of 17 N's starts at pos 33 (1-based)
        let unique_seed = b"NNNNNNNNNNNNNNNNN";
        let positions = reference.seed.get(unique_seed.as_slice());
        assert!(positions.is_some(), "Unique seed must be indexed");
        let pos_list = positions.unwrap();
        assert_eq!(
            pos_list.len(),
            1,
            "Unique seed must resolve to exactly one position"
        );
        assert_eq!(pos_list[0], 33, "Unique seed position must be 33");

        // Cleanup
        let _ = fs::remove_file(&fa_path);
        let _ = fs::remove_file(&fai_path);
    }
}
