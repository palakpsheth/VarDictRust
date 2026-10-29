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

struct CachedFastaReader {
    reader: IndexedReader<BufReader<File>>,
}

thread_local! {
    static FASTA_READER_CACHE: RefCell<HashMap<PathBuf, CachedFastaReader>> = RefCell::new(HashMap::new());
}

/// Clears the current thread's cached FASTA readers.
pub(crate) fn clear_fasta_reader_cache() {
    FASTA_READER_CACHE.with(|cache| cache.borrow_mut().clear());
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
        if !cache.contains_key(fasta_path) {
            let mut fai_path = fasta_path.to_path_buf();
            let mut ext = fai_path
                .extension()
                .map(|e| e.to_string_lossy().into_owned())
                .unwrap_or_default();
            ext.push_str(".fai");
            fai_path.set_extension(ext);

            let index = fai::fs::read(&fai_path).map_err(VardictError::Io)?;
            let file = File::open(fasta_path).map_err(VardictError::Io)?;
            cache.insert(
                fasta_path.to_path_buf(),
                CachedFastaReader {
                    reader: IndexedReader::new(BufReader::new(file), index),
                },
            );
        }
        let Some(cached) = cache.get_mut(fasta_path) else {
            return Err(VardictError::Other(
                "FASTA reader cache entry missing after insert".to_string(),
            ));
        };

        // NoodlesRegion parses "chr:start-end"
        let region_str = format!("{}:{}-{}", chr, start.0, end.0);
        let region: NoodlesRegion = region_str.parse().map_err(|e| {
            VardictError::Parse(format!("failed to parse region '{}': {}", region_str, e))
        })?;

        let record = cached.reader.query(&region).map_err(VardictError::Io)?;

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

/// Inline k-mer key up to 17 bytes (`SEED_1`), avoiding per-k-mer heap allocations.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct SeedKey {
    len: u8,
    bytes: [u8; 17], // bytes[len as usize..] ALWAYS zero-filled
}

impl SeedKey {
    #[inline]
    fn make_key(key: &[u8]) -> Option<Self> {
        if key.len() > 17 {
            return None;
        }
        let mut bytes = [0u8; 17];
        bytes.get_mut(..key.len())?.copy_from_slice(key);
        Some(Self {
            len: key.len() as u8,
            bytes,
        })
    }
}

/// Inline single-position storage for unique k-mers, spilling to `Vec<u32>` only on repeats.
#[derive(Clone, Debug, PartialEq, Eq)]
enum SeedPositions {
    One(u32),
    Many(Vec<u32>),
}

/// Deferred k-mer indexing chunk captured by [`get_reference`] so tiles that
/// never query `reference.seed` perform zero k-mer hashing (`EBO-2`).
#[derive(Clone, Debug)]
struct PendingSeedChunk {
    bases: std::sync::Arc<[u8]>,
    sequence_start: u32,
    sequence_end: u32,
    chr_len: u32,
    prev_bounds: Vec<(u32, u32)>,
}

type RawSeedMap = rustc_hash::FxHashMap<SeedKey, SeedPositions>;

/// Lookup-only k-mer seed table mapping 12-mer and 17-mer byte slices to 1-based positions.
///
/// K-mer hashing is deferred until the first query (`get`, `len`, `is_empty`),
/// using `OnceCell` so `get(&self)` returns `Option<&[u32]>` tied to `&self`
/// without runtime borrow guards or `unsafe` (`M-UNSAFE`, `M-MEM-REUSE`).
#[derive(Debug, Default)]
pub struct SeedMap {
    map: std::cell::OnceCell<RawSeedMap>,
    pending: std::cell::RefCell<(RawSeedMap, Vec<PendingSeedChunk>)>,
}

impl Clone for SeedMap {
    fn clone(&self) -> Self {
        let map = std::cell::OnceCell::new();
        if let Some(existing) = self.map.get() {
            let _ = map.set(existing.clone());
        }
        let pending = match self.pending.try_borrow() {
            Ok(guard) => std::cell::RefCell::new(guard.clone()),
            Err(_) => std::cell::RefCell::new((RawSeedMap::default(), Vec::new())),
        };
        Self { map, pending }
    }
}

impl SeedMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns `true` if the k-mer lookup table has been materialized.
    #[inline]
    pub fn is_materialized(&self) -> bool {
        self.map.get().is_some()
    }

    fn insert_into_raw(map: &mut RawSeedMap, key: &[u8], pos: u32) {
        let Some(sk) = SeedKey::make_key(key) else {
            debug_assert!(key.len() <= 17, "seed key exceeds 17 bytes");
            return;
        };
        match map.entry(sk) {
            std::collections::hash_map::Entry::Occupied(mut e) => match e.get_mut() {
                SeedPositions::One(first) => {
                    let first_val = *first;
                    *e.into_mut() = SeedPositions::Many(vec![first_val, pos]);
                }
                SeedPositions::Many(v) => v.push(pos),
            },
            std::collections::hash_map::Entry::Vacant(v) => {
                v.insert(SeedPositions::One(pos));
            }
        }
    }

    /// Ports `ReferenceResource.java:106-125` k-mer indexing for one loaded chunk.
    fn index_chunk(map: &mut RawSeedMap, chunk: &PendingSeedChunk) {
        const SEED1_LEN: usize = 17;
        const SEED2_LEN: usize = 12;

        let bases = chunk.bases.as_ref();
        let site_end = if chunk.chr_len == chunk.sequence_end {
            bases.len()
        } else {
            bases.len().saturating_sub(SEED1_LEN)
        };

        for i in 0..site_end {
            let pos = i as u32 + chunk.sequence_start;

            // Don't index positions already covered by a previously loaded region
            if chunk
                .prev_bounds
                .iter()
                .any(|&(start, end)| pos >= start && pos <= end)
            {
                continue;
            }

            // Do not create adaptor sequences for the very end of chromosome
            if chunk.chr_len == chunk.sequence_end && i.saturating_add(SEED1_LEN) > bases.len() {
                continue;
            }

            if let Some(key) = bases.get(i..i + SEED1_LEN) {
                Self::insert_into_raw(map, key, pos);
            }

            if let Some(key) = bases.get(i..i + SEED2_LEN) {
                Self::insert_into_raw(map, key, pos);
            }
        }
    }

    #[inline]
    fn materialize(&self) -> &RawSeedMap {
        self.map.get_or_init(|| {
            let Ok(mut guard) = self.pending.try_borrow_mut() else {
                return RawSeedMap::default();
            };
            let (mut map, mut chunks) = std::mem::take(&mut *guard);
            for chunk in chunks.drain(..) {
                Self::index_chunk(&mut map, &chunk);
            }
            map
        })
    }

    #[inline]
    pub fn get(&self, key: &[u8]) -> Option<&[u32]> {
        let sk = SeedKey::make_key(key)?;
        match self.materialize().get(&sk)? {
            SeedPositions::One(p) => Some(std::slice::from_ref(p)),
            SeedPositions::Many(v) => Some(v.as_slice()),
        }
    }

    pub fn len(&self) -> usize {
        self.materialize().len()
    }

    pub fn is_empty(&self) -> bool {
        self.materialize().is_empty()
    }

    fn enqueue_chunk(&mut self, chunk: PendingSeedChunk) {
        let staged = self.pending.get_mut();
        if let Some(existing) = self.map.take() {
            staged.0 = existing;
        }
        staged.1.push(chunk);
    }

    #[cfg(test)]
    pub(crate) fn push_pos(&mut self, key: &[u8], pos: u32) {
        let staged = self.pending.get_mut();
        if let Some(existing) = self.map.take() {
            staged.0 = existing;
        }
        for chunk in staged.1.drain(..) {
            Self::index_chunk(&mut staged.0, &chunk);
        }
        Self::insert_into_raw(&mut staged.0, key, pos);
    }

    #[cfg(test)]
    pub fn insert_for_test(&mut self, key: &[u8], positions: Vec<u32>) {
        let sk = SeedKey::make_key(key).expect("test seed key exceeds 17 bytes");
        let staged = self.pending.get_mut();
        if let Some(existing) = self.map.take() {
            staged.0 = existing;
        }
        for chunk in staged.1.drain(..) {
            Self::index_chunk(&mut staged.0, &chunk);
        }
        staged.0.remove(&sk);
        for p in positions {
            self.push_pos(key, p);
        }
    }

    #[cfg(test)]
    pub fn insert(&mut self, key: Vec<u8>, positions: Vec<u32>) {
        self.insert_for_test(&key, positions);
    }
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
    pub seed: SeedMap,
}

impl Reference {
    pub fn new() -> Self {
        Reference {
            loaded_regions: Vec::new(),
            reference_sequences: HashMap::new(),
            region_sequences: Vec::new(),
            seed: SeedMap::new(),
        }
    }

    #[inline]
    fn seq_for_idx<'a>(&'a self, idx: usize, chr: &str) -> Option<&'a [u8]> {
        if let Some(bases) = self.region_sequences.get(idx) {
            Some(bases.as_slice())
        } else if idx == 0 {
            self.reference_sequences.get(chr).map(Vec::as_slice)
        } else {
            None
        }
    }

    /// Safely gets the base at a 1-based chromosome position.
    pub fn get_base(&self, chr: &str, pos: u32) -> Option<u8> {
        for (idx, lr) in self.loaded_regions.iter().enumerate() {
            if lr.chr == chr && pos >= lr.sequence_start && pos <= lr.sequence_end {
                let offset = (pos - lr.sequence_start) as usize;
                if let Some(bases) = self.seq_for_idx(idx, chr) {
                    if let Some(&b) = bases.get(offset) {
                        return Some(b);
                    }
                }
            }
        }
        None
    }

    /// Selects the stored `LoadedRegion` sequence slice on `chr` covering a read's
    /// alignment + soft-clip span (`ext_start..=ext_end`, falling back to
    /// `read_start..=read_end`, any overlapping region, then the first loaded region on `chr`).
    pub fn region_slice_for_span(
        &self,
        chr: &str,
        read_start: i32,
        read_end: i32,
        ext_start: i32,
        ext_end: i32,
    ) -> (&[u8], i32) {
        let candidates = || {
            self.loaded_regions
                .iter()
                .enumerate()
                .filter(|(_, lr)| lr.chr == chr)
                .filter_map(|(idx, lr)| {
                    let seq = self.seq_for_idx(idx, chr)?;
                    let lr_start = lr.sequence_start as i32;
                    let lr_end = lr_start.saturating_add(seq.len() as i32);
                    Some((seq, lr_start, lr_end))
                })
        };
        candidates()
            .find(|&(_, lr_start, lr_end)| ext_start >= lr_start && ext_end <= lr_end)
            .or_else(|| {
                candidates()
                    .find(|&(_, lr_start, lr_end)| read_start >= lr_start && read_end <= lr_end)
            })
            .or_else(|| {
                candidates()
                    .find(|&(_, lr_start, lr_end)| read_end > lr_start && read_start < lr_end)
            })
            .or_else(|| candidates().next())
            .map(|(seq, lr_start, _)| (seq, lr_start))
            .unwrap_or((&[], 1))
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
    let chr_len = chromosomes_lengths.get(&region.chr).copied().unwrap_or(0) as u32;
    get_reference_with_chr_len(region, extension, chr_len, config, reference)
}

/// Ports `ReferenceResource.java:getReference` with a pre-resolved chromosome length (`M-MEM-REUSE`).
pub fn get_reference_with_chr_len(
    region: &BedRegion,
    extension: u32,
    len: u32,
    config: &Config,
    reference: &mut Reference,
) -> Result<(), VardictError> {
    let start_val = region.start.0;
    let end_val = region.end.0;

    // Ports ReferenceResource.java:85-89: uses `instance().conf.numberNucleotideToExtend`
    // (`-x` = `config.extension`, default 0), not `vext` (`-X` = `config.extension_nuc`).
    let total_ext = (config.extension.max(0) as u32).saturating_add(extension);
    let sequence_start = start_val.saturating_sub(total_ext).max(1);

    let end_ext = end_val.saturating_add(total_ext);
    let sequence_end = if len > 0 && end_ext > len {
        len
    } else {
        end_ext
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

    const SEED1_LEN: usize = 17;
    // Ports ReferenceResource.java:109 `if (ref.referenceSequences.containsKey(i + sequenceStart))`:
    // previously loaded regions only populated `ref.referenceSequences` up to their `siteEnd`
    // (their stored length in `region_sequences`), not the trailing 17-bp k-mer lookahead.
    let prev_bounds: Vec<(u32, u32)> = reference
        .loaded_regions
        .iter()
        .enumerate()
        .filter(|(_, lr)| lr.chr == region.chr)
        .filter_map(|(idx, lr)| {
            let stored_len = reference.region_sequences.get(idx).map_or(0, Vec::len) as u32;
            stored_len
                .checked_sub(1)
                .map(|last| (lr.sequence_start, lr.sequence_start.saturating_add(last)))
        })
        .collect();

    let loaded_region = LoadedRegion {
        chr: region.chr.clone(),
        sequence_start,
        sequence_end,
    };
    reference.loaded_regions.push(loaded_region);

    // Defer k-mer seed indexing (`ReferenceResource.java:106-125`) until
    // `reference.seed` is actually queried (`EBO-2`). The full untrimmed `bases`
    // (including the trailing `SEED_1` lookahead) is passed to `PendingSeedChunk`
    // so 17-mer keys at `i < site_end` can slice `bases[i..i + 17]`.
    reference.seed.enqueue_chunk(PendingSeedChunk {
        bases: std::sync::Arc::from(bases.as_slice()),
        sequence_start,
        sequence_end,
        chr_len: len,
        prev_bounds,
    });

    // Ports ReferenceResource.java:106-112:
    // `int siteEnd = len == sequenceEnd ? exon.length() : exon.length() - Configuration.SEED_1;`
    // Only `[0..siteEnd)` is inserted into `ref.referenceSequences`.
    let site_end = if len == sequence_end {
        bases.len()
    } else {
        bases.len().saturating_sub(SEED1_LEN)
    };
    bases.truncate(site_end);

    // Keep reference_sequences populated for primary chromosome region
    reference
        .reference_sequences
        .entry(region.chr.clone())
        .or_insert_with(|| bases.clone());

    reference.region_sequences.push(bases);

    Ok(())
}

#[cfg(test)]
#[expect(
    clippy::field_reassign_with_default,
    reason = "tests configure individual Config fields after Config::default()"
)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    #[derive(Default)]
    struct ReferenceLegacy {
        loaded_regions: Vec<LoadedRegion>,
        reference_sequences: HashMap<String, Vec<u8>>,
        region_sequences: Vec<Vec<u8>>,
        seed: HashMap<Vec<u8>, Vec<u32>>,
    }

    impl ReferenceLegacy {
        fn get_base(&self, chr: &str, pos: u32) -> Option<u8> {
            for (idx, lr) in self.loaded_regions.iter().enumerate() {
                if lr.chr == chr && pos >= lr.sequence_start && pos <= lr.sequence_end {
                    let offset = (pos - lr.sequence_start) as usize;
                    if let Some(bases) = self.region_sequences.get(idx) {
                        if let Some(&b) = bases.get(offset) {
                            return Some(b);
                        }
                    } else if idx == 0 {
                        if let Some(bases) = self.reference_sequences.get(chr) {
                            if let Some(&b) = bases.get(offset) {
                                return Some(b);
                            }
                        }
                    }
                }
            }
            None
        }
    }

    fn get_reference_legacy(
        region: &BedRegion,
        extension: u32,
        chromosomes_lengths: &HashMap<String, usize>,
        config: &Config,
        reference: &mut ReferenceLegacy,
    ) -> Result<(), VardictError> {
        let start_val = region.start.0;
        let end_val = region.end.0;
        let total_ext = (config.extension.max(0) as u32).saturating_add(extension);
        let sequence_start = start_val.saturating_sub(total_ext).max(1);
        let len = chromosomes_lengths.get(&region.chr).copied().unwrap_or(0) as u32;
        let end_ext = end_val.saturating_add(total_ext);
        let sequence_end = if len > 0 && end_ext > len {
            len
        } else {
            end_ext
        };
        if reference.loaded_regions.iter().any(|lr| {
            lr.chr == region.chr
                && sequence_start >= lr.sequence_start
                && sequence_end <= lr.sequence_end
        }) {
            return Ok(());
        }

        let fasta_path = Path::new(&config.reference);
        let mut bases = fetch_subsequence(
            fasta_path,
            &region.chr,
            Pos1(sequence_start),
            Pos1(sequence_end),
        )?;
        for b in &mut bases {
            b.make_ascii_uppercase();
        }
        let seed1_len = 17;
        let seed2_len = 12;
        let site_end = if len == sequence_end {
            bases.len()
        } else {
            bases.len().saturating_sub(seed1_len)
        };
        let prev_bounds: Vec<(u32, u32)> = reference
            .loaded_regions
            .iter()
            .enumerate()
            .filter(|(_, lr)| lr.chr == region.chr)
            .filter_map(|(idx, lr)| {
                let stored_len = reference.region_sequences.get(idx).map_or(0, Vec::len) as u32;
                stored_len
                    .checked_sub(1)
                    .map(|last| (lr.sequence_start, lr.sequence_start.saturating_add(last)))
            })
            .collect();
        for i in 0..site_end {
            let pos = i as u32 + sequence_start;
            if prev_bounds
                .iter()
                .any(|&(start, end)| pos >= start && pos <= end)
            {
                continue;
            }
            if len == sequence_end && i > bases.len().saturating_sub(seed1_len) {
                continue;
            }
            if i + seed1_len <= bases.len() {
                let key = &bases[i..i + seed1_len];
                reference.seed.entry(key.to_vec()).or_default().push(pos);
            }
            if i + seed2_len <= bases.len() {
                let key = &bases[i..i + seed2_len];
                reference.seed.entry(key.to_vec()).or_default().push(pos);
            }
        }
        let mut stored_bases = bases;
        stored_bases.truncate(site_end);
        reference.loaded_regions.push(LoadedRegion {
            chr: region.chr.clone(),
            sequence_start,
            sequence_end,
        });
        reference
            .reference_sequences
            .entry(region.chr.clone())
            .or_insert_with(|| stored_bases.clone());
        reference.region_sequences.push(stored_bases);
        Ok(())
    }

    #[test]
    fn seed_map_and_get_reference_match_legacy() {
        // 1. Direct SeedKey / SeedMap edge cases: length discriminator & >17 bounds
        let mut sm = SeedMap::new();
        assert!(sm.is_empty());
        sm.push_pos(b"AC", 42);
        assert!(!sm.is_empty());
        assert_eq!(sm.len(), 1);
        assert_eq!(sm.get(b"AC"), Some(&[42u32][..]));
        assert_eq!(sm.get(b"AC\0"), None);
        assert_ne!(sm.get(b"AC"), sm.get(b"AC\0"));
        sm.push_pos(b"AC\0", 99);
        assert_eq!(sm.get(b"AC"), Some(&[42u32][..]));
        assert_eq!(sm.get(b"AC\0"), Some(&[99u32][..]));
        assert_eq!(sm.get(b""), None);
        assert_eq!(sm.get(&[b'A'; 11]), None);
        assert_eq!(sm.get(&[b'A'; 13]), None);
        assert_eq!(sm.get(&[b'A'; 18]), None);

        // One -> Many transition
        sm.push_pos(b"AC", 77);
        sm.push_pos(b"AC", 123);
        assert_eq!(sm.get(b"AC"), Some(&[42u32, 77, 123][..]));

        // 2. 5 kb seeded-random FASTA DIFF test vs ReferenceLegacy
        let temp_dir = std::env::temp_dir();
        let pid = std::process::id();
        let fa_path = temp_dir.join(format!("test_seed_diff_{pid}.fa"));
        let fai_path = temp_dir.join(format!("test_seed_diff_{pid}.fa.fai"));

        let alphabet = b"ACGTACGTACGTNRYSWKM";
        let mut rng: u64 = 0x5EED_CAFE_1234_5678;
        let mut next_u64 = || -> u64 {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };

        let total_len = 5000usize;
        let mut seq = Vec::with_capacity(total_len);
        for i in 0..total_len {
            // Inject repeated 17-mers and 12-mers every 400 bp plus lowercase bases
            if i % 400 < 20 {
                let b = b"acgtACGTNNRYSWKMAAAA"[i % 400];
                seq.push(b);
            } else {
                let idx = (next_u64() as usize) % alphabet.len();
                let mut b = alphabet[idx];
                if (next_u64() & 3) == 0 {
                    b = b.to_ascii_lowercase();
                }
                seq.push(b);
            }
        }

        let mut fa_file = File::create(&fa_path).unwrap();
        fa_file.write_all(b">chr1\n").unwrap();
        fa_file.write_all(&seq).unwrap();
        fa_file
            .write_all(b"\n>chr2\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
            .unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(
            fai_file,
            "chr1\t{total_len}\t6\t{total_len}\t{}",
            total_len + 1
        )
        .unwrap();
        writeln!(fai_file, "chr2\t32\t{}\t32\t33", 6 + total_len + 1 + 6).unwrap();

        let mut config = Config::default();
        config.reference = fa_path.to_string_lossy().into_owned();
        config.extension = 50;
        config.extension_nuc = 7; // distinct from `-x` so regression to `extension_nuc` fails immediately

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), total_len);
        chr_lengths.insert("chr2".to_string(), 32);

        let mut new_ref = Reference::new();
        let mut leg_ref = ReferenceLegacy::default();

        // Multi-region overlapping loads + duplicate enclosed load + chromosome-end load + second contig
        let regions = [
            ("chr1", 100u32, 1400u32, 25u32),
            ("chr1", 200, 800, 10),   // already enclosed -> no-op
            ("chr1", 1100, 2800, 30), // overlapping with first region
            ("chr2", 1, 32, 0),       // other chromosome (interleaved)
            ("chr1", 2600, 5000, 40), // reaches sequence_end == len (5000)
        ];

        for &(chr, start, end, ext) in &regions {
            let bed = BedRegion {
                chr: chr.to_string(),
                start: Pos1(start),
                end: Pos1(end),
                gene: "G".to_string(),
                amplicon: None,
            };
            get_reference(&bed, ext, &chr_lengths, &config, &mut new_ref).unwrap();
            get_reference_legacy(&bed, ext, &chr_lengths, &config, &mut leg_ref).unwrap();

            assert_eq!(new_ref.loaded_regions, leg_ref.loaded_regions);
            assert_eq!(
                new_ref.loaded_regions[0].sequence_start, 25,
                "get_reference must extend by config.extension (50) + extension (25) => 100 - 75 = 25, ignoring config.extension_nuc (7)"
            );
            assert_eq!(new_ref.reference_sequences, leg_ref.reference_sequences);
            assert_eq!(new_ref.region_sequences, leg_ref.region_sequences);
            assert_eq!(new_ref.seed.len(), leg_ref.seed.len());

            for (k, leg_positions) in &leg_ref.seed {
                assert_eq!(
                    new_ref.seed.get(k.as_slice()),
                    Some(leg_positions.as_slice()),
                    "seed mismatch for key {:?}",
                    String::from_utf8_lossy(k)
                );
            }
        }

        // Verify all 12-mers and 17-mers across the uppercased 5 kb sequence + negative probes
        let upper_seq: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        for klen in [0usize, 11, 12, 13, 17, 18] {
            if klen == 0 {
                assert_eq!(new_ref.seed.get(b""), None);
                continue;
            }
            for i in 0..=upper_seq.len().saturating_sub(klen) {
                let k = &upper_seq[i..i + klen];
                assert_eq!(
                    new_ref.seed.get(k),
                    leg_ref.seed.get(k).map(Vec::as_slice),
                    "probe mismatch at i={i}, klen={klen}"
                );
            }
        }

        for pos in 1..=(total_len as u32 + 10) {
            assert_eq!(
                new_ref.get_base("chr1", pos),
                leg_ref.get_base("chr1", pos),
                "get_base mismatch at pos={pos}"
            );
        }

        let _ = fs::remove_file(&fa_path);
        let _ = fs::remove_file(&fai_path);
    }

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
        config.extension = 2;

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 40);

        let mut reference = Reference::new();
        get_reference(&bed_region, 3, &chr_lengths, &config, &mut reference).unwrap();

        // start=5, extend_nuc=2, extension=3 => sequence_start = 5 - 5 = 1 (since 5 <= 5, clamped at 1)
        // end=25, extend_nuc=2, extension=3 => sequence_end = 25 + 5 = 30
        // Since sequence_end (30) < chr_len (40), ReferenceResource.java:106-112 sets
        // siteEnd = 30 - SEED_1 (17) = 13, so positions 1..=13 are stored in referenceSequences
        // while 14..=30 serve only as the 17-mer lookahead window.
        assert_eq!(reference.loaded_regions.len(), 1);
        assert_eq!(reference.loaded_regions[0].sequence_start, 1);
        assert_eq!(reference.loaded_regions[0].sequence_end, 30);

        // Verify base lookup matches ReferenceResource.java:106-112 (siteEnd = 13)
        assert_eq!(reference.get_base("chr1", 1), Some(b'A'));
        assert_eq!(reference.get_base("chr1", 2), Some(b'C'));
        assert_eq!(reference.get_base("chr1", 13), Some(b'A'));
        assert_eq!(reference.get_base("chr1", 14), None); // trailing SEED_1 lookahead excluded
        assert_eq!(reference.get_base("chr1", 30), None);
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

    #[test]
    fn test_lazy_seed_map_defers_indexing_until_queried() {
        let temp_dir = std::env::temp_dir();
        let pid = std::process::id();
        let fa_path = temp_dir.join(format!("test_lazy_seed_{pid}.fa"));
        let fai_path = temp_dir.join(format!("test_lazy_seed_{pid}.fa.fai"));

        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNNNNNNNNNNNNACGTACGTACGT";
        let total_len = seq.len();
        let mut fa_file = File::create(&fa_path).unwrap();
        fa_file.write_all(b">chr1\n").unwrap();
        fa_file.write_all(seq).unwrap();
        fa_file.write_all(b"\n").unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(
            fai_file,
            "chr1\t{total_len}\t6\t{total_len}\t{}",
            total_len + 1
        )
        .unwrap();

        let mut config = Config::default();
        config.reference = fa_path.to_string_lossy().into_owned();

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), total_len);

        let mut reference = Reference::new();
        assert!(!reference.seed.is_materialized());

        // Load first sub-region [1, 50]: should NOT materialize SeedMap yet.
        let r1 = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(50),
            gene: "G".to_string(),
            amplicon: None,
        };
        get_reference(&r1, 0, &chr_lengths, &config, &mut reference).unwrap();
        assert!(
            !reference.seed.is_materialized(),
            "get_reference must defer k-mer hashing until seed is queried"
        );
        // Cloning an unmaterialized Reference must preserve laziness.
        let cloned_unmaterialized = reference.clone();
        assert!(!cloned_unmaterialized.seed.is_materialized());

        // First query triggers materialization and returns exact position.
        // Note: r1 [1, 50] has sequence_end=50 < total_len=61, so its site_end is 50 - 17 = 33
        // (positions 1..=33), while 17-mer at position 34 ("NNNNNNNNNNNNNNNNT", inside r1's
        // truncated 17-bp lookahead [34..=50]) is NOT indexed until r2 [30, 61] is loaded.
        let unique_seed = b"NNNNNNNNNNNNNNNNN";
        let lookahead_seed = b"NNNNNNNNNNNNNNNNA";
        assert_eq!(reference.seed.get(unique_seed), Some(&[33u32][..]));
        assert_eq!(
            reference.seed.get(lookahead_seed),
            None,
            "17-mer starting in r1's truncated 17-bp lookahead (pos 34) must not be indexed by r1 alone"
        );
        assert_eq!(
            reference.get_base("chr1", 33),
            Some(b'N'),
            "pos 33 is the last base inside r1's site_end"
        );
        assert_eq!(
            reference.get_base("chr1", 34),
            None,
            "pos 34 is in r1's truncated 17-bp lookahead and must not be stored in region_sequences[0]"
        );
        assert!(reference.seed.is_materialized());
        assert_eq!(
            cloned_unmaterialized.seed.get(unique_seed),
            Some(&[33u32][..])
        );
        assert!(cloned_unmaterialized.seed.is_materialized());

        // Subsequent overlapping get_reference [30, 61] after materialization
        // must re-defer the new chunk and include positions [34..=50] that were in
        // r1's truncated lookahead, without duplicating positions in [1..=33].
        let r2 = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(30),
            end: Pos1(total_len as u32),
            gene: "G".to_string(),
            amplicon: None,
        };
        get_reference(&r2, 0, &chr_lengths, &config, &mut reference).unwrap();
        assert!(
            !reference.seed.is_materialized(),
            "second get_reference must queue chunk lazily"
        );
        assert_eq!(
            reference.seed.get(lookahead_seed),
            Some(&[34u32][..]),
            "overlapping r2 must index pos 34 (which was inside r1's truncated 17-bp lookahead)"
        );
        assert_eq!(
            reference.get_base("chr1", 34),
            Some(b'N'),
            "pos 34 must now be served by r2"
        );

        let mut leg_ref = ReferenceLegacy::default();
        get_reference_legacy(&r1, 0, &chr_lengths, &config, &mut leg_ref).unwrap();
        get_reference_legacy(&r2, 0, &chr_lengths, &config, &mut leg_ref).unwrap();

        assert_eq!(reference.seed.len(), leg_ref.seed.len());
        assert!(reference.seed.is_materialized());
        for (k, leg_positions) in &leg_ref.seed {
            assert_eq!(
                reference.seed.get(k.as_slice()),
                Some(leg_positions.as_slice())
            );
        }

        // Also verify that a non-terminal window of <= 17 bp with zero extensions
        // (which truncates to `stored_len == 0`) does not underflow `stored_len - 1`
        // when a subsequent disjoint region is loaded on the same chromosome.
        let mut short_cfg = Config::default();
        short_cfg.reference = fa_path.to_str().unwrap().to_string();
        short_cfg.extension = 0;
        short_cfg.extension_nuc = 0;
        short_cfg.indel_size = 0;
        let r_tiny = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(10),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut ref_tiny = Reference::new();
        get_reference(&r_tiny, 0, &chr_lengths, &short_cfg, &mut ref_tiny).unwrap();
        assert_eq!(ref_tiny.region_sequences[0].len(), 0);
        get_reference(&r2, 0, &chr_lengths, &short_cfg, &mut ref_tiny).unwrap();

        let mut leg_tiny = ReferenceLegacy::default();
        get_reference_legacy(&r_tiny, 0, &chr_lengths, &short_cfg, &mut leg_tiny).unwrap();
        get_reference_legacy(&r2, 0, &chr_lengths, &short_cfg, &mut leg_tiny).unwrap();
        assert_eq!(ref_tiny.seed.len(), leg_tiny.seed.len());

        let _ = fs::remove_file(&fa_path);
        let _ = fs::remove_file(&fai_path);
    }
}
