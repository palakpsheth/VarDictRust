//! Reference FASTA and k-mer seed map loading.
//!
//! Ports Reference.java and ReferenceResource.java.

use crate::config::Config;
use crate::regions::Region as BedRegion;
use crate::types::{Pos1, VardictError};
use noodles::fasta::fai;
use noodles::fasta::io::IndexedReader;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

const MAX_GZI_BYTES: u64 = 128 * 1024 * 1024;
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const BGZF_CM_FLG: [u8; 2] = [0x08, 0x04];
const BGZF_XLEN_SUBFIELD: [u8; 6] = [0x06, 0x00, b'B', b'C', 0x02, 0x00];

struct CachedFastaReader {
    reader: IndexedReader<noodles::fasta::io::BufReader<File>>,
    fai_lengths: HashMap<Vec<u8>, u64>,
}

thread_local! {
    static FASTA_READER_CACHE: RefCell<HashMap<PathBuf, CachedFastaReader>> = RefCell::new(HashMap::new());
}

/// Clears the current thread's cached FASTA readers.
pub(crate) fn clear_fasta_reader_cache() {
    FASTA_READER_CACHE.with(|cache| cache.borrow_mut().clear());
}

/// Appends `suffix` (e.g. `".fai"`, `".gzi"`) directly to `path`'s OS string without
/// `Path::set_extension` double-dot artifacts on extensionless paths (`genome` -> `genome.fai`).
pub(crate) fn append_os_ext(path: &Path, suffix: &str) -> PathBuf {
    let mut os = path.as_os_str().to_os_string();
    os.push(suffix);
    PathBuf::from(os)
}

fn existing_fai_candidate_for(path: &Path) -> Option<PathBuf> {
    let primary = append_os_ext(path, ".fai");
    if std::fs::symlink_metadata(&primary).is_ok() {
        return Some(primary);
    }
    if let Some(ext) = path.extension() {
        let is_compressed = ext.eq_ignore_ascii_case("gz") || ext.eq_ignore_ascii_case("bgz");
        if !is_compressed {
            let stem_fai = path.with_extension("fai");
            if stem_fai != primary
                && stem_fai != path
                && std::fs::symlink_metadata(&stem_fai).is_ok()
            {
                return Some(stem_fai);
            }
        }
    }
    None
}

/// Resolves the `.fai` index path for `fasta_path`: checks `<fasta_path>.fai` first (using
/// `symlink_metadata` so a dangling primary symlink fails loudly instead of being silently
/// skipped), then `<stem>.fai` when uncompressed, then the same candidates beside
/// `std::fs::canonicalize(fasta_path)` when `fasta_path` is a symlink, falling back to
/// `<fasta_path>.fai`.
pub(crate) fn resolve_fasta_fai_path(fasta_path: &Path) -> PathBuf {
    if let Some(found) = existing_fai_candidate_for(fasta_path) {
        return found;
    }
    if let Ok(canon) = std::fs::canonicalize(fasta_path) {
        if canon != fasta_path {
            if let Some(found) = existing_fai_candidate_for(&canon) {
                return found;
            }
        }
    }
    append_os_ext(fasta_path, ".fai")
}

/// Resolves the `.gzi` index path for a BGZF-compressed FASTA: checks `<fasta_path>.gzi` first
/// (using `symlink_metadata` so a dangling primary symlink fails loudly), then `<canon>.gzi`
/// when `fasta_path` is a symlink, falling back to `<fasta_path>.gzi`.
pub(crate) fn resolve_fasta_gzi_path(fasta_path: &Path) -> PathBuf {
    let primary = append_os_ext(fasta_path, ".gzi");
    if std::fs::symlink_metadata(&primary).is_ok() {
        return primary;
    }
    if let Ok(canon) = std::fs::canonicalize(fasta_path) {
        if canon != fasta_path {
            let canon_gzi = append_os_ext(&canon, ".gzi");
            if std::fs::symlink_metadata(&canon_gzi).is_ok() {
                return canon_gzi;
            }
        }
    }
    primary
}

fn read_up_to<R: Read>(reader: &mut R, buf: &mut [u8]) -> std::io::Result<usize> {
    let mut total = 0usize;
    while let Some(dst) = buf.get_mut(total..) {
        if dst.is_empty() {
            break;
        }
        match reader.read(dst) {
            Ok(0) => break,
            Ok(n) => total = total.saturating_add(n),
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(total)
}

fn read_validated_gzi_index(
    gzi_path: &Path,
    fasta_path: &Path,
    fasta_len: u64,
) -> Result<noodles_bgzf::gzi::Index, VardictError> {
    let mut gzi_file = match File::open(gzi_path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            return Err(VardictError::Parse(format!(
                "BGZF-compressed FASTA '{}' is missing its required GZI index '{}'; run `samtools faidx {}` to generate both .fai and .gzi indexes",
                fasta_path.display(),
                gzi_path.display(),
                fasta_path.display()
            )));
        }
        Err(e) => return Err(VardictError::Io(e)),
    };

    let len = gzi_file.metadata().map_err(VardictError::Io)?.len();
    let Some(payload_len) = len.checked_sub(8) else {
        return Err(VardictError::Parse(format!(
            "Corrupt GZI index '{}' ({len} bytes < 8-byte header) for '{}'; regenerate with `samtools faidx {}`",
            gzi_path.display(),
            fasta_path.display(),
            fasta_path.display()
        )));
    };
    if payload_len % 16 != 0 || len > MAX_GZI_BYTES {
        return Err(VardictError::Parse(format!(
            "Invalid or oversized GZI index '{}' ({len} bytes) for '{}'; regenerate with `samtools faidx {}`",
            gzi_path.display(),
            fasta_path.display(),
            fasta_path.display()
        )));
    }

    let mut count_bytes = [0u8; 8];
    let n_read = read_up_to(&mut gzi_file, &mut count_bytes).map_err(VardictError::Io)?;
    if n_read != 8 {
        return Err(VardictError::Parse(format!(
            "Truncated GZI index '{}' (read {n_read}/8 header bytes); regenerate with `samtools faidx {}`",
            gzi_path.display(),
            fasta_path.display()
        )));
    }
    let declared_entries = u64::from_le_bytes(count_bytes);
    let expected_entries = payload_len / 16;
    if declared_entries != expected_entries {
        return Err(VardictError::Parse(format!(
            "Corrupt GZI index '{}': header declares {declared_entries} entries, but file size ({len} bytes) holds {expected_entries} entries; regenerate with `samtools faidx {}`",
            gzi_path.display(),
            fasta_path.display()
        )));
    }

    gzi_file
        .seek(SeekFrom::Start(0))
        .map_err(VardictError::Io)?;
    let index = noodles_bgzf::gzi::io::Reader::new(BufReader::new(gzi_file))
        .read_index()
        .map_err(|e| {
            VardictError::Parse(format!(
                "Corrupt GZI index '{}' for '{}': {e}; regenerate with `samtools faidx {}`",
                gzi_path.display(),
                fasta_path.display(),
                fasta_path.display()
            ))
        })?;

    let pairs = index.as_ref();
    let explicit_pairs = if pairs.first() == Some(&(0, 0))
        && (pairs.len() as u64) == declared_entries.saturating_add(1)
    {
        pairs.get(1..).unwrap_or(&[])
    } else if declared_entries == 0 && pairs.is_empty() {
        &[][..]
    } else {
        pairs
    };

    let mut prev_c = 0u64;
    let mut prev_u = 0u64;
    for &(c, u) in explicit_pairs {
        if c <= prev_c || u <= prev_u {
            return Err(VardictError::Parse(format!(
                "Corrupt GZI index '{}' for '{}': offsets are not strictly monotonic at ({c}, {u}) after ({prev_c}, {prev_u}); regenerate with `samtools faidx {}`",
                gzi_path.display(),
                fasta_path.display(),
                fasta_path.display()
            )));
        }
        if fasta_len > 0 && c >= fasta_len {
            return Err(VardictError::Parse(format!(
                "Corrupt GZI index '{}' for '{}': compressed offset {c} >= FASTA file length ({fasta_len} bytes); regenerate with `samtools faidx {}`",
                gzi_path.display(),
                fasta_path.display(),
                fasta_path.display()
            )));
        }
        prev_c = c;
        prev_u = u;
    }

    Ok(index)
}

/// Opens an indexed FASTA buffered reader with an explicit `BufReader` capacity for uncompressed
/// FASTA, supporting both uncompressed FASTA (including `< 18`-byte files) and BGZF-compressed
/// FASTA (`.fa.gz` / `.fa.bgz`) paired with a pre-validated `.gzi` index.
pub(crate) fn open_indexed_fasta_buf_reader_with_capacity(
    fasta_path: &Path,
    capacity: usize,
) -> Result<noodles::fasta::io::BufReader<File>, VardictError> {
    let mut file = File::open(fasta_path).map_err(VardictError::Io)?;
    let mut hdr = [0u8; 18];
    let n = read_up_to(&mut file, &mut hdr).map_err(VardictError::Io)?;

    if n < 2 || hdr.get(..2) != Some(GZIP_MAGIC.as_slice()) {
        file.seek(SeekFrom::Start(0)).map_err(VardictError::Io)?;
        return Ok(noodles::fasta::io::BufReader::Uncompressed(
            BufReader::with_capacity(capacity, file),
        ));
    }

    let is_bgzf = n == 18
        && hdr.get(2..4) == Some(BGZF_CM_FLG.as_slice())
        && hdr.get(10..16) == Some(BGZF_XLEN_SUBFIELD.as_slice());
    if !is_bgzf {
        return Err(VardictError::Parse(format!(
            "Reference FASTA '{}' is plain gzip, not BGZF; re-compress with bgzip and index with samtools faidx",
            fasta_path.display()
        )));
    }

    let fasta_len = file.metadata().map_err(VardictError::Io)?.len();
    let gzi_path = resolve_fasta_gzi_path(fasta_path);
    let gzi_index = read_validated_gzi_index(&gzi_path, fasta_path, fasta_len)?;
    file.seek(SeekFrom::Start(0)).map_err(VardictError::Io)?;
    Ok(noodles::fasta::io::BufReader::Bgzf(
        noodles_bgzf::io::IndexedReader::new(file, gzi_index),
    ))
}

/// Opens an indexed FASTA buffered reader with the default 8 KiB uncompressed buffer capacity.
pub(crate) fn open_indexed_fasta_buf_reader(
    fasta_path: &Path,
) -> Result<noodles::fasta::io::BufReader<File>, VardictError> {
    open_indexed_fasta_buf_reader_with_capacity(fasta_path, 8 * 1024)
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
            let fasta_reader = open_indexed_fasta_buf_reader(fasta_path)?;
            let fai_path = resolve_fasta_fai_path(fasta_path);
            let index = fai::fs::read(&fai_path).map_err(VardictError::Io)?;
            let fai_lengths: HashMap<Vec<u8>, u64> = index
                .as_ref()
                .iter()
                .map(|rec| (rec.name().to_vec(), rec.length()))
                .collect();
            cache.insert(
                fasta_path.to_path_buf(),
                CachedFastaReader {
                    reader: IndexedReader::new(fasta_reader, index),
                    fai_lengths,
                },
            );
        }
        let Some(cached) = cache.get_mut(fasta_path) else {
            return Err(VardictError::Other(
                "FASTA reader cache entry missing after insert".to_string(),
            ));
        };

        if let Some(&fai_len) = cached.fai_lengths.get(chr.as_bytes()) {
            if u64::from(start.0) > fai_len || u64::from(end.0) > fai_len {
                return Err(VardictError::Parse(format!(
                    "Query asks for data past end of contig {}: requested {}..={}, contig length {} in '{}'",
                    chr,
                    start.0,
                    end.0,
                    fai_len,
                    fasta_path.display()
                )));
            }
        }

        let region = crate::bam::build_noodles_region(chr, start.0, end.0)?;

        let record = cached.reader.query(&region).map_err(VardictError::Io)?;
        let seq = record.sequence().as_ref();
        let expected_len = if end.0 >= start.0 {
            (end.0 - start.0) as usize + 1
        } else {
            0
        };
        if seq.len() < expected_len {
            return Err(VardictError::Parse(format!(
                "Query asks for data past end of contig {}: requested {}..={} ({} bp), got {} bp from '{}'",
                chr,
                start.0,
                end.0,
                expected_len,
                seq.len(),
                fasta_path.display()
            )));
        }

        Ok(seq.to_vec())
    })
}

/// Represents a loaded contiguous block of the reference sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoadedRegion {
    pub chr: String,
    pub sequence_start: u32,
    pub sequence_end: u32,
}

const MAX_SEED_BYTES: usize = 17;

/// Word-packed k-mer key up to 17 bytes (`SEED_1`), storing `0..=17` payload bytes
/// in bytes `0..17` and the explicit length discriminator in byte `23` for 3-word
/// `[u64; 3]` comparison and `FxHasher` hashing without heap allocation.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SeedKey([u64; 3]);

impl SeedKey {
    #[inline]
    fn from_slice(key: &[u8]) -> Option<Self> {
        if key.len() > MAX_SEED_BYTES {
            return None;
        }
        let mut buf = [0u8; 24];
        buf.get_mut(..key.len())?.copy_from_slice(key);
        let [b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, _b23] =
            buf;
        Some(Self([
            u64::from_ne_bytes([b0, b1, b2, b3, b4, b5, b6, b7]),
            u64::from_ne_bytes([b8, b9, b10, b11, b12, b13, b14, b15]),
            u64::from_ne_bytes([b16, b17, b18, b19, b20, b21, b22, key.len() as u8]),
        ]))
    }
}

impl std::hash::Hash for SeedKey {
    #[inline]
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let [w0, w1, w2] = self.0;
        state.write_u64(w0);
        state.write_u64(w1);
        state.write_u64(w2);
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
        let Some(sk) = SeedKey::from_slice(key) else {
            debug_assert!(key.len() <= MAX_SEED_BYTES, "seed key exceeds 17 bytes");
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

        let mut i = 0;
        while i < site_end {
            let pos = chunk
                .sequence_start
                .saturating_add(u32::try_from(i).unwrap_or(u32::MAX));

            // Don't index positions already covered by a previously loaded region.
            // Jump directly past the covered interval instead of stepping 1 bp at a time (#63 H2).
            if let Some(&(_, end)) = chunk
                .prev_bounds
                .iter()
                .find(|&&(start, end)| pos >= start && pos <= end)
            {
                i = (end.saturating_sub(chunk.sequence_start) as usize).saturating_add(1);
                continue;
            }

            // Do not create adaptor sequences for the very end of chromosome
            if chunk.chr_len == chunk.sequence_end && i.saturating_add(SEED1_LEN) > bases.len() {
                i += 1;
                continue;
            }

            if let Some(key) = bases.get(i..i + SEED1_LEN) {
                Self::insert_into_raw(map, key, pos);
            }

            if let Some(key) = bases.get(i..i + SEED2_LEN) {
                Self::insert_into_raw(map, key, pos);
            }
            i += 1;
        }
    }

    const MAX_SEED_RESERVE: usize = 262_144;

    #[inline]
    fn estimate_seed_reserve<I: IntoIterator<Item = usize>>(non_n_counts: I) -> usize {
        non_n_counts
            .into_iter()
            .fold(0usize, |acc, non_n| {
                acc.saturating_add(non_n.saturating_mul(2))
            })
            .min(Self::MAX_SEED_RESERVE)
    }

    #[inline]
    fn materialize(&self) -> &RawSeedMap {
        self.map.get_or_init(|| {
            let Ok(mut guard) = self.pending.try_borrow_mut() else {
                return RawSeedMap::default();
            };
            let (mut map, mut chunks) = std::mem::take(&mut *guard);
            let est = Self::estimate_seed_reserve(chunks.iter().map(|chunk| {
                chunk
                    .bases
                    .iter()
                    .filter(|&&b| b != b'N' && b != b'n')
                    .count()
            }));
            map.reserve(est);
            for chunk in chunks.drain(..) {
                Self::index_chunk(&mut map, &chunk);
            }
            map
        })
    }

    #[inline]
    pub fn get(&self, key: &[u8]) -> Option<&[u32]> {
        let sk = SeedKey::from_slice(key)?;
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
        let sk = SeedKey::from_slice(key).expect("test seed key exceeds 17 bytes");
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

#[derive(Debug, Clone)]
enum RegionSpanSlot {
    Shared(std::sync::Arc<[u8]>, usize),
    RegionSeq(usize),
    Unpaired,
}

/// Stores reference bases and k-mer seed maps for sequence realignment.
#[derive(Debug, Clone)]
pub struct Reference {
    pub loaded_regions: Vec<LoadedRegion>,
    /// Maps chromosome name to contiguous bases (upper-case) for the primary loaded region.
    pub reference_sequences: HashMap<String, Vec<u8>>,
    /// Sequences matching loaded_regions by index: region_sequences[i] corresponds to loaded_regions[i].
    pub region_sequences: Vec<Vec<u8>>,
    /// Per-`loaded_regions` backing-storage descriptor, sharing the same `Arc<[u8]>` queued in `seed`
    /// for FASTA-loaded regions while tracking synthetic `region_sequences` / `reference_sequences` slots.
    shared_region_spans: Vec<RegionSpanSlot>,
    /// Maps k-mer bytes to 1-based positions.
    pub seed: SeedMap,
}

impl Reference {
    pub fn new() -> Self {
        Reference {
            loaded_regions: Vec::new(),
            reference_sequences: HashMap::new(),
            region_sequences: Vec::new(),
            shared_region_spans: Vec::new(),
            seed: SeedMap::new(),
        }
    }

    fn sync_untracked_slots(&mut self, up_to_loaded_len: usize) {
        let cur_slots = self.shared_region_spans.len();
        if cur_slots >= up_to_loaded_len {
            return;
        }
        let mut next_ri = self
            .shared_region_spans
            .iter()
            .filter(|s| matches!(s, RegionSpanSlot::RegionSeq(_)))
            .count();
        let avail_reg_seqs = self.region_sequences.len().saturating_sub(next_ri);
        let untracked = up_to_loaded_len.saturating_sub(cur_slots);
        let mut idx = cur_slots;
        if cur_slots == 0
            && avail_reg_seqs.saturating_add(1) == untracked
            && self
                .loaded_regions
                .first()
                .is_some_and(|lr0| self.reference_sequences.contains_key(&lr0.chr))
        {
            self.shared_region_spans.push(RegionSpanSlot::Unpaired);
            idx = 1;
        }
        for _ in idx..up_to_loaded_len {
            if next_ri < self.region_sequences.len() {
                self.shared_region_spans
                    .push(RegionSpanSlot::RegionSeq(next_ri));
                next_ri += 1;
            } else {
                self.shared_region_spans.push(RegionSpanSlot::Unpaired);
            }
        }
    }

    #[inline]
    fn seq_for_idx<'a>(&'a self, idx: usize, chr: &str) -> Option<&'a [u8]> {
        if let Some(slot) = self.shared_region_spans.get(idx) {
            return match slot {
                RegionSpanSlot::Shared(arc, len) => arc.get(..*len),
                RegionSpanSlot::RegionSeq(ri) => self.region_sequences.get(*ri).map(Vec::as_slice),
                RegionSpanSlot::Unpaired => (idx == 0)
                    .then(|| self.reference_sequences.get(chr).map(Vec::as_slice))
                    .flatten(),
            };
        }
        if self.shared_region_spans.is_empty()
            && self.region_sequences.len() == self.loaded_regions.len()
        {
            return self.region_sequences.get(idx).map(Vec::as_slice);
        }
        let cur_slots = self.shared_region_spans.len();
        let paired_reg_seqs = self
            .shared_region_spans
            .iter()
            .filter(|s| matches!(s, RegionSpanSlot::RegionSeq(_)))
            .count();
        let avail_reg_seqs = self.region_sequences.len().saturating_sub(paired_reg_seqs);
        let untracked = self.loaded_regions.len().saturating_sub(cur_slots);
        if cur_slots == 0
            && avail_reg_seqs.saturating_add(1) == untracked
            && self
                .loaded_regions
                .first()
                .is_some_and(|lr0| self.reference_sequences.contains_key(&lr0.chr))
        {
            if idx == 0 {
                self.reference_sequences.get(chr).map(Vec::as_slice)
            } else {
                self.region_sequences
                    .get(idx.saturating_sub(1))
                    .map(Vec::as_slice)
            }
        } else {
            let ri = paired_reg_seqs.saturating_add(idx.saturating_sub(cur_slots));
            self.region_sequences
                .get(ri)
                .map(Vec::as_slice)
                .or_else(|| {
                    (idx == 0)
                        .then(|| self.reference_sequences.get(chr).map(Vec::as_slice))
                        .flatten()
                })
        }
    }

    /// Returns the stored sequence slice and 1-based `sequence_start` for the first
    /// `LoadedRegion` matching `chr`.
    #[inline]
    pub fn primary_region_slice(&self, chr: &str) -> Option<(&[u8], u32)> {
        let (idx, lr) = self
            .loaded_regions
            .iter()
            .enumerate()
            .find(|(_, r)| r.chr == chr)?;
        if lr.sequence_end < lr.sequence_start {
            return None;
        }
        let slice = self.seq_for_idx(idx, chr)?;
        let max_len = (lr.sequence_end.saturating_sub(lr.sequence_start) as usize)
            .saturating_add(1)
            .min(slice.len());
        Some((slice.get(..max_len)?, lr.sequence_start))
    }

    /// Safely gets the base at a 1-based chromosome position.
    #[inline]
    pub fn get_base(&self, chr: &str, pos: u32) -> Option<u8> {
        if let Some(lr) = self.loaded_regions.first() {
            if pos >= lr.sequence_start && pos <= lr.sequence_end && lr.chr == chr {
                if let Some(offset) = pos
                    .checked_sub(lr.sequence_start)
                    .and_then(|d| usize::try_from(d).ok())
                {
                    if let Some(bases) = self.seq_for_idx(0, chr) {
                        if let Some(&b) = bases.get(offset) {
                            return Some(b);
                        }
                    }
                }
            }
        }
        for (idx, lr) in self.loaded_regions.iter().enumerate().skip(1) {
            if lr.chr == chr && pos >= lr.sequence_start && pos <= lr.sequence_end {
                if let Some(offset) = pos
                    .checked_sub(lr.sequence_start)
                    .and_then(|d| usize::try_from(d).ok())
                {
                    if let Some(bases) = self.seq_for_idx(idx, chr) {
                        if let Some(&b) = bases.get(offset) {
                            return Some(b);
                        }
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
        if self.loaded_regions.len() == 1 {
            if let Some(lr) = self.loaded_regions.first() {
                if lr.chr == chr {
                    if let Some(seq) = self.seq_for_idx(0, chr) {
                        return (seq, lr.sequence_start as i32);
                    }
                }
            }
            return (&[], 1);
        }
        let mut read_match: Option<(&[u8], i32)> = None;
        let mut overlap_match: Option<(&[u8], i32)> = None;
        let mut first_match: Option<(&[u8], i32)> = None;
        for (idx, lr) in self.loaded_regions.iter().enumerate() {
            if lr.chr != chr {
                continue;
            }
            let Some(seq) = self.seq_for_idx(idx, chr) else {
                continue;
            };
            let lr_start = lr.sequence_start as i32;
            let lr_end = lr_start.saturating_add(seq.len() as i32);
            if ext_start >= lr_start && ext_end <= lr_end {
                return (seq, lr_start);
            }
            if read_match.is_none() && read_start >= lr_start && read_end <= lr_end {
                read_match = Some((seq, lr_start));
            }
            if overlap_match.is_none() && read_end > lr_start && read_start < lr_end {
                overlap_match = Some((seq, lr_start));
            }
            if first_match.is_none() {
                first_match = Some((seq, lr_start));
            }
        }
        read_match
            .or(overlap_match)
            .or(first_match)
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
    let chr_len = u32::try_from(chromosomes_lengths.get(&region.chr).copied().unwrap_or(0))
        .unwrap_or(u32::MAX);
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
    let prev_loaded_len = reference.loaded_regions.len();
    reference.sync_untracked_slots(prev_loaded_len);

    // Ports ReferenceResource.java:109 `if (ref.referenceSequences.containsKey(i + sequenceStart))`:
    // previously loaded regions only populated `ref.referenceSequences` up to their `siteEnd`
    // (their stored length via `seq_for_idx`), not the trailing 17-bp k-mer lookahead.
    let prev_bounds: Vec<(u32, u32)> = reference
        .loaded_regions
        .iter()
        .enumerate()
        .filter(|(_, lr)| lr.chr == region.chr)
        .filter_map(|(idx, lr)| {
            let stored_len = reference.seq_for_idx(idx, &lr.chr).map_or(0, <[u8]>::len) as u32;
            stored_len
                .checked_sub(1)
                .map(|last| (lr.sequence_start, lr.sequence_start.saturating_add(last)))
        })
        .collect();

    // Ports ReferenceResource.java:106-112:
    // `int siteEnd = len == sequenceEnd ? exon.length() : exon.length() - Configuration.SEED_1;`
    // Only `[0..siteEnd)` is exposed as the region sequence, while the full untrimmed
    // `shared_bases` (including the trailing `SEED_1` lookahead) is shared via `Arc<[u8]>`
    // with `PendingSeedChunk` so 17-mer keys at `i < site_end` can slice `bases[i..i + 17]`
    // without allocating duplicate `Vec<u8>` buffers.
    let site_end = if len == sequence_end {
        bases.len()
    } else {
        bases.len().saturating_sub(SEED1_LEN)
    };
    let shared_bases: std::sync::Arc<[u8]> = bases.into();

    let loaded_region = LoadedRegion {
        chr: region.chr.clone(),
        sequence_start,
        sequence_end,
    };
    reference.loaded_regions.push(loaded_region);

    // Defer k-mer seed indexing (`ReferenceResource.java:106-125`) until
    // `reference.seed` is actually queried (`EBO-2`).
    reference.seed.enqueue_chunk(PendingSeedChunk {
        bases: std::sync::Arc::clone(&shared_bases),
        sequence_start,
        sequence_end,
        chr_len: len,
        prev_bounds,
    });

    reference
        .shared_region_spans
        .push(RegionSpanSlot::Shared(shared_bases, site_end));

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
        // 1. Direct SeedKey / SeedMap edge cases: word-packed [u64; 3], length discriminator & >17 bounds
        assert_eq!(std::mem::size_of::<SeedKey>(), 24);
        for len in 0..=MAX_SEED_BYTES {
            let k = vec![b'G'; len];
            assert!(
                SeedKey::from_slice(&k).is_some(),
                "SeedKey::from_slice must accept len={len}"
            );
        }
        assert_eq!(SeedKey::from_slice(&[b'A'; 18]), None);
        assert_eq!(SeedKey::from_slice(&[b'A'; 64]), None);
        {
            use std::hash::{BuildHasher, Hasher};
            let sk_acgt = SeedKey::from_slice(b"ACGT").unwrap();
            let sk_acgt_nul = SeedKey::from_slice(b"ACGT\0").unwrap();
            assert_ne!(sk_acgt, sk_acgt_nul);
            let bh = rustc_hash::FxBuildHasher;
            let mut h1 = bh.build_hasher();
            std::hash::Hash::hash(&sk_acgt, &mut h1);
            let mut h2 = bh.build_hasher();
            std::hash::Hash::hash(&sk_acgt_nul, &mut h2);
            assert_ne!(h1.finish(), h2.finish());
        }

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
        sm.push_pos(b"ACGT", 7);
        sm.push_pos(b"ACGT\0", 8);
        assert_eq!(sm.get(b"ACGT"), Some(&[7u32][..]));
        assert_eq!(sm.get(b"ACGT\0"), Some(&[8u32][..]));
        sm.push_pos(&[b'T'; 17], 555);
        assert_eq!(sm.get(&[b'T'; 17]), Some(&[555u32][..]));
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
            assert!(
                new_ref.reference_sequences.is_empty(),
                "get_reference must not duplicate bases into reference_sequences"
            );
            assert!(
                new_ref.region_sequences.is_empty(),
                "get_reference must not duplicate bases into region_sequences"
            );
            assert_eq!(
                new_ref.shared_region_spans.len(),
                leg_ref.region_sequences.len()
            );
            for (idx, lr) in new_ref.loaded_regions.iter().enumerate() {
                assert_eq!(
                    new_ref.seq_for_idx(idx, &lr.chr),
                    Some(leg_ref.region_sequences[idx].as_slice())
                );
                for pos in lr.sequence_start.saturating_sub(1)..=lr.sequence_end.saturating_add(1) {
                    assert_eq!(
                        new_ref.get_base(&lr.chr, pos),
                        leg_ref.get_base(&lr.chr, pos),
                        "get_base mismatch on {} at pos={pos}",
                        lr.chr
                    );
                }
            }
            let expected_primary = leg_ref
                .loaded_regions
                .iter()
                .enumerate()
                .find(|(_, r)| r.chr == chr)
                .map(|(i, r)| (leg_ref.region_sequences[i].as_slice(), r.sequence_start));
            assert_eq!(new_ref.primary_region_slice(chr), expected_primary);
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
        assert!(reference.reference_sequences.is_empty());
        assert!(reference.region_sequences.is_empty());
        assert_eq!(reference.shared_region_spans.len(), 1);
        let RegionSpanSlot::Shared(ref span0_arc, span0_len) = reference.shared_region_spans[0]
        else {
            panic!("expected RegionSpanSlot::Shared at slot 0");
        };
        assert_eq!(span0_len, 13);
        {
            let pending_guard = reference.seed.pending.borrow();
            assert_eq!(pending_guard.1.len(), 1);
            assert!(
                std::ptr::eq(
                    std::sync::Arc::as_ptr(span0_arc),
                    std::sync::Arc::as_ptr(&pending_guard.1[0].bases),
                ),
                "Reference.shared_region_spans and PendingSeedChunk must share a single Arc<[u8]> allocation"
            );
        }
        let (prim_slice, prim_start) = reference.primary_region_slice("chr1").unwrap();
        assert_eq!(prim_start, 1);
        assert_eq!(prim_slice, b"ACGTACGTACGTA");
        assert!(std::ptr::eq(
            prim_slice.as_ptr(),
            std::sync::Arc::as_ptr(span0_arc).cast::<u8>(),
        ));
        assert_eq!(reference.primary_region_slice("chr_missing"), None);

        // Verify base lookup matches ReferenceResource.java:106-112 (siteEnd = 13)
        assert_eq!(reference.get_base("chr1", 1), Some(b'A'));
        assert_eq!(reference.get_base("chr1", 2), Some(b'C'));
        assert_eq!(reference.get_base("chr1", 13), Some(b'A'));
        assert_eq!(reference.get_base("chr1", 14), None); // trailing SEED_1 lookahead excluded
        assert_eq!(reference.get_base("chr1", 30), None);
        assert_eq!(reference.get_base("chr1", 31), None); // out of bounds

        // Verify seed indexing contains our expected 17-mers and 12-mers
        assert!(!reference.seed.is_empty());

        // Verify mixed-order synchronization between synthetic regions and get_reference_with_chr_len:
        let mut mix_cfg = config.clone();
        mix_cfg.extension = 0;

        // 1. [R0, SH] -> [R0, SH, RS] -> [R0, SH, RS, SH]
        let mut ref_r0_sh_rs = Reference::new();
        ref_r0_sh_rs.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 6,
        });
        ref_r0_sh_rs
            .reference_sequences
            .insert("chr1".to_string(), b"TTTTTT".to_vec());
        let reg_mid = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(9),
            end: Pos1(28),
            gene: "MID".to_string(),
            amplicon: None,
        };
        get_reference_with_chr_len(&reg_mid, 0, 40, &mix_cfg, &mut ref_r0_sh_rs).unwrap();
        assert!(ref_r0_sh_rs.region_sequences.is_empty());
        assert!(matches!(
            ref_r0_sh_rs.shared_region_spans[0],
            RegionSpanSlot::Unpaired
        ));
        assert!(matches!(
            ref_r0_sh_rs.shared_region_spans[1],
            RegionSpanSlot::Shared(_, 3)
        ));

        // Append RS at slot 2 -> [R0, SH, RS]
        ref_r0_sh_rs.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 29,
            sequence_end: 32,
        });
        ref_r0_sh_rs.region_sequences.push(b"GGGG".to_vec());
        ref_r0_sh_rs.loaded_regions.push(LoadedRegion {
            chr: "chr2".to_string(),
            sequence_start: 10,
            sequence_end: 13,
        });
        ref_r0_sh_rs.region_sequences.push(b"CCCC".to_vec());
        assert_eq!(ref_r0_sh_rs.seq_for_idx(0, "chr1"), Some(&b"TTTTTT"[..]));
        assert_eq!(ref_r0_sh_rs.seq_for_idx(1, "chr1"), Some(&b"ACG"[..]));
        assert_eq!(ref_r0_sh_rs.seq_for_idx(2, "chr1"), Some(&b"GGGG"[..]));
        assert_eq!(ref_r0_sh_rs.seq_for_idx(3, "chr2"), Some(&b"CCCC"[..]));
        assert_eq!(
            ref_r0_sh_rs.primary_region_slice("chr1"),
            Some((&b"TTTTTT"[..], 1))
        );
        assert_eq!(
            ref_r0_sh_rs.primary_region_slice("chr2"),
            Some((&b"CCCC"[..], 10))
        );
        for pos in 1..=6 {
            assert_eq!(ref_r0_sh_rs.get_base("chr1", pos), Some(b'T'));
        }
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 9), Some(b'A'));
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 10), Some(b'C'));
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 11), Some(b'G'));
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 12), None);
        for pos in 29..=32 {
            assert_eq!(ref_r0_sh_rs.get_base("chr1", pos), Some(b'G'));
        }
        for pos in 10..=13 {
            assert_eq!(ref_r0_sh_rs.get_base("chr2", pos), Some(b'C'));
        }

        // Append SH at slot 4 -> [R0, SH, RS, RS, SH]
        let reg_end = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(35),
            end: Pos1(40),
            gene: "END".to_string(),
            amplicon: None,
        };
        get_reference_with_chr_len(&reg_end, 0, 40, &mix_cfg, &mut ref_r0_sh_rs).unwrap();
        assert_eq!(ref_r0_sh_rs.region_sequences.len(), 2);
        assert_eq!(ref_r0_sh_rs.shared_region_spans.len(), 5);
        assert_eq!(ref_r0_sh_rs.seq_for_idx(0, "chr1"), Some(&b"TTTTTT"[..]));
        assert_eq!(ref_r0_sh_rs.seq_for_idx(1, "chr1"), Some(&b"ACG"[..]));
        assert_eq!(ref_r0_sh_rs.seq_for_idx(2, "chr1"), Some(&b"GGGG"[..]));
        assert_eq!(ref_r0_sh_rs.seq_for_idx(3, "chr2"), Some(&b"CCCC"[..]));
        assert_eq!(ref_r0_sh_rs.seq_for_idx(4, "chr1"), Some(&b"GTACGT"[..]));
        assert_eq!(
            ref_r0_sh_rs.primary_region_slice("chr1"),
            Some((&b"TTTTTT"[..], 1))
        );
        assert_eq!(
            ref_r0_sh_rs.primary_region_slice("chr2"),
            Some((&b"CCCC"[..], 10))
        );
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 1), Some(b'T'));
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 11), Some(b'G'));
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 30), Some(b'G'));
        assert_eq!(ref_r0_sh_rs.get_base("chr2", 12), Some(b'C'));
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 35), Some(b'G'));
        assert_eq!(ref_r0_sh_rs.get_base("chr1", 40), Some(b'T'));

        // 2. [R0, RS, SH] across single-contig and multi-contig primary_region_slice + get_base:
        let mut ref_r0_rs_sh = Reference::new();
        ref_r0_rs_sh.loaded_regions.push(LoadedRegion {
            chr: "chr0".to_string(),
            sequence_start: 1,
            sequence_end: 6,
        });
        ref_r0_rs_sh
            .reference_sequences
            .insert("chr0".to_string(), b"TTTTTT".to_vec());
        ref_r0_rs_sh.loaded_regions.push(LoadedRegion {
            chr: "chr2".to_string(),
            sequence_start: 10,
            sequence_end: 15,
        });
        ref_r0_rs_sh.region_sequences.push(b"GGGGGG".to_vec());
        assert_eq!(
            ref_r0_rs_sh.primary_region_slice("chr0"),
            Some((&b"TTTTTT"[..], 1))
        );
        assert_eq!(
            ref_r0_rs_sh.primary_region_slice("chr2"),
            Some((&b"GGGGGG"[..], 10))
        );
        let reg_tail = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(21),
            end: Pos1(40),
            gene: "TAIL".to_string(),
            amplicon: None,
        };
        get_reference_with_chr_len(&reg_tail, 0, 40, &mix_cfg, &mut ref_r0_rs_sh).unwrap();
        assert_eq!(ref_r0_rs_sh.region_sequences.len(), 1);
        assert_eq!(ref_r0_rs_sh.shared_region_spans.len(), 3);
        assert_eq!(
            ref_r0_rs_sh.primary_region_slice("chr0"),
            Some((&b"TTTTTT"[..], 1))
        );
        assert_eq!(
            ref_r0_rs_sh.primary_region_slice("chr2"),
            Some((&b"GGGGGG"[..], 10))
        );
        assert_eq!(
            ref_r0_rs_sh.primary_region_slice("chr1"),
            Some((&b"ACGTACGTACGTACGTACGT"[..], 21))
        );
        assert_eq!(ref_r0_rs_sh.get_base("chr0", 1), Some(b'T'));
        assert_eq!(ref_r0_rs_sh.get_base("chr2", 10), Some(b'G'));
        assert_eq!(ref_r0_rs_sh.get_base("chr1", 21), Some(b'A'));
        assert_eq!(ref_r0_rs_sh.get_base("chr1", 40), Some(b'T'));

        let mut ref_same_chr_r0_rs_sh = Reference::new();
        ref_same_chr_r0_rs_sh.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 6,
        });
        ref_same_chr_r0_rs_sh
            .reference_sequences
            .insert("chr1".to_string(), b"TTTTTT".to_vec());
        ref_same_chr_r0_rs_sh.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 10,
            sequence_end: 15,
        });
        ref_same_chr_r0_rs_sh
            .region_sequences
            .push(b"GGGGGG".to_vec());
        assert_eq!(
            ref_same_chr_r0_rs_sh.primary_region_slice("chr1"),
            Some((&b"TTTTTT"[..], 1))
        );
        assert_eq!(ref_same_chr_r0_rs_sh.get_base("chr1", 1), Some(b'T'));
        assert_eq!(ref_same_chr_r0_rs_sh.get_base("chr1", 10), Some(b'G'));
        get_reference_with_chr_len(&reg_tail, 0, 40, &mix_cfg, &mut ref_same_chr_r0_rs_sh).unwrap();
        assert_eq!(
            ref_same_chr_r0_rs_sh.primary_region_slice("chr1"),
            Some((&b"TTTTTT"[..], 1))
        );
        assert_eq!(ref_same_chr_r0_rs_sh.get_base("chr1", 1), Some(b'T'));
        assert_eq!(ref_same_chr_r0_rs_sh.get_base("chr1", 10), Some(b'G'));
        assert_eq!(ref_same_chr_r0_rs_sh.get_base("chr1", 21), Some(b'A'));
        assert_eq!(ref_same_chr_r0_rs_sh.get_base("chr1", 40), Some(b'T'));

        // 3. [SH, RS]
        let mut ref_mix_b = Reference::new();
        let reg_head = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(25),
            gene: "HEAD".to_string(),
            amplicon: None,
        };
        get_reference_with_chr_len(&reg_head, 0, 40, &mix_cfg, &mut ref_mix_b).unwrap();
        ref_mix_b.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 31,
            sequence_end: 36,
        });
        ref_mix_b.region_sequences.push(b"GGGGGG".to_vec());
        assert_eq!(ref_mix_b.region_sequences.len(), 1);
        assert_eq!(ref_mix_b.shared_region_spans.len(), 1);
        assert_eq!(ref_mix_b.seq_for_idx(0, "chr1"), Some(&b"ACGTACGT"[..]));
        assert_eq!(ref_mix_b.seq_for_idx(1, "chr1"), Some(&b"GGGGGG"[..]));
        assert_eq!(
            ref_mix_b.primary_region_slice("chr1"),
            Some((&b"ACGTACGT"[..], 1))
        );
        assert_eq!(ref_mix_b.get_base("chr1", 1), Some(b'A'));
        assert_eq!(ref_mix_b.get_base("chr1", 8), Some(b'T'));
        assert_eq!(ref_mix_b.get_base("chr1", 9), None);
        assert_eq!(ref_mix_b.get_base("chr1", 31), Some(b'G'));
        assert_eq!(ref_mix_b.get_base("chr1", 36), Some(b'G'));

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
        assert!(ref_tiny.region_sequences.is_empty());
        assert!(ref_tiny.reference_sequences.is_empty());
        assert!(matches!(
            ref_tiny.shared_region_spans[0],
            RegionSpanSlot::Shared(_, 0)
        ));
        assert_eq!(ref_tiny.primary_region_slice("chr1"), Some((&[][..], 1)));
        get_reference(&r2, 0, &chr_lengths, &short_cfg, &mut ref_tiny).unwrap();

        let mut leg_tiny = ReferenceLegacy::default();
        get_reference_legacy(&r_tiny, 0, &chr_lengths, &short_cfg, &mut leg_tiny).unwrap();
        get_reference_legacy(&r2, 0, &chr_lengths, &short_cfg, &mut leg_tiny).unwrap();
        assert_eq!(ref_tiny.seed.len(), leg_tiny.seed.len());

        let _ = fs::remove_file(&fa_path);
        let _ = fs::remove_file(&fai_path);
    }

    #[test]
    fn fetch_subsequence_past_end_is_err() {
        let dir = tempfile::tempdir().unwrap();
        let fa_path = dir.path().join("short.fa");
        let fai_path = dir.path().join("short.fa.fai");

        let seq_1000 = b"ACGT".repeat(250);
        let mut fa_file = File::create(&fa_path).unwrap();
        fa_file.write_all(b">chr1\n").unwrap();
        fa_file.write_all(&seq_1000).unwrap();
        fa_file.write_all(b"\n").unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        fai_file.write_all(b"chr1\t1000\t6\t1000\t1001\n").unwrap();

        // Query in-bounds succeeds.
        let in_bounds = fetch_subsequence(&fa_path, "chr1", Pos1(1), Pos1(1000)).unwrap();
        assert_eq!(in_bounds.len(), 1000);

        // Query overlapping past the 1000-bp contig end must return Err, not a truncated slice.
        let err_overlap = fetch_subsequence(&fa_path, "chr1", Pos1(1), Pos1(2000)).unwrap_err();
        assert!(
            err_overlap
                .to_string()
                .contains("Query asks for data past end of contig chr1"),
            "unexpected error: {err_overlap}"
        );

        // Query starting completely past the contig end must return Err, not an empty Vec.
        let err_past = fetch_subsequence(&fa_path, "chr1", Pos1(1500), Pos1(2000)).unwrap_err();
        assert!(
            err_past
                .to_string()
                .contains("Query asks for data past end of contig chr1"),
            "unexpected error: {err_past}"
        );

        // Truncated .fa with stale .fai claiming 2000 bp must also return Err when fewer bytes are read.
        let stale_fa_path = dir.path().join("stale.fa");
        let stale_fai_path = dir.path().join("stale.fa.fai");
        fs::copy(&fa_path, &stale_fa_path).unwrap();
        fs::write(&stale_fai_path, b"chr1\t2000\t6\t2000\t2001\n").unwrap();
        let err_stale =
            fetch_subsequence(&stale_fa_path, "chr1", Pos1(900), Pos1(1200)).unwrap_err();
        assert!(
            err_stale
                .to_string()
                .contains("Query asks for data past end of contig chr1"),
            "unexpected error: {err_stale}"
        );
    }

    #[test]
    fn get_base_primary_fastpath_and_truncated_fallthrough() {
        let mut reference = Reference::new();
        // Primary region claims 1..=50, but region_sequences[0] is truncated to 33 bp
        // (simulating SEED1_LEN = 17 trim in get_reference_with_chr_len).
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 50,
        });
        reference.region_sequences.push(vec![b'A'; 33]);

        // Secondary overlapping region covers 30..=60 with 'G'.
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 30,
            sequence_end: 60,
        });
        reference.region_sequences.push(vec![b'G'; 31]);

        // 1. Primary region fast-path hit (within stored 1..=33).
        assert_eq!(reference.get_base("chr1", 1), Some(b'A'));
        assert_eq!(reference.get_base("chr1", 33), Some(b'A'));

        // 2. Fallthrough to loaded_regions[1] when pos <= loaded_regions[0].sequence_end (50)
        //    but past region_sequences[0].len() (33).
        assert_eq!(reference.get_base("chr1", 34), Some(b'G'));
        assert_eq!(reference.get_base("chr1", 50), Some(b'G'));
        assert_eq!(reference.get_base("chr1", 60), Some(b'G'));

        // 3. Out-of-bounds or different chromosome returns None, and synthetic region_sequences
        //    fallback works through primary_region_slice.
        assert_eq!(reference.get_base("chr1", 0), None);
        assert_eq!(reference.get_base("chr1", 61), None);
        assert_eq!(reference.get_base("chr2", 10), None);
        assert_eq!(
            reference.primary_region_slice("chr1"),
            Some((vec![b'A'; 33].as_slice(), 1))
        );
        assert_eq!(reference.primary_region_slice("chr2"), None);
    }

    #[test]
    fn seed_map_materialize_bounded_reserve_and_parity() {
        // 1. Overflow-safe reserve estimator caps at 262_144 even for usize::MAX inputs
        assert_eq!(SeedMap::estimate_seed_reserve([0usize, 0]), 0);
        assert_eq!(SeedMap::estimate_seed_reserve([500usize, 250]), 1_500);
        assert_eq!(SeedMap::estimate_seed_reserve([200_000usize]), 262_144);
        assert_eq!(
            SeedMap::estimate_seed_reserve([usize::MAX, usize::MAX / 2 + 17]),
            262_144
        );

        // 2. All-'N'/'n' stretches reserve 0 upfront (capacity stays tiny, <= 262_144)
        //    while returning bit-for-bit identical 12-mer/17-mer positions.
        {
            let mut sm_n = SeedMap::new();
            let n_len = 200_000usize;
            let mixed_n: Vec<u8> = (0..n_len)
                .map(|i| if i % 2 == 0 { b'N' } else { b'n' })
                .collect();
            sm_n.enqueue_chunk(PendingSeedChunk {
                bases: std::sync::Arc::from(mixed_n.as_slice()),
                sequence_start: 1,
                sequence_end: n_len as u32,
                chr_len: n_len as u32,
                prev_bounds: Vec::new(),
            });
            let raw = sm_n.materialize();
            assert!(
                raw.capacity() <= 64,
                "all-'N'/'n' stretch must not over-reserve capacity (<= 64 <= 262_144, got {})",
                raw.capacity()
            );

            // End-to-end uppercase 'N' stretch via get_reference vs ReferenceLegacy
            let dir = tempfile::tempdir().unwrap();
            let fa_path = dir.path().join("all_n.fa");
            let fai_path = dir.path().join("all_n.fa.fai");
            let fa_seq = vec![b'N'; n_len];
            let mut fa = File::create(&fa_path).unwrap();
            fa.write_all(b">chr1\n").unwrap();
            fa.write_all(&fa_seq).unwrap();
            fa.write_all(b"\n").unwrap();
            let mut fai = File::create(&fai_path).unwrap();
            writeln!(fai, "chr1\t{n_len}\t6\t{n_len}\t{}", n_len + 1).unwrap();

            let mut cfg = Config::default();
            cfg.reference = fa_path.to_str().unwrap().to_string();
            let mut chr_lengths = HashMap::new();
            chr_lengths.insert("chr1".to_string(), n_len);
            let bed = BedRegion {
                chr: "chr1".to_string(),
                start: Pos1(1),
                end: Pos1(n_len as u32),
                gene: "N".to_string(),
                amplicon: None,
            };
            let mut new_ref = Reference::new();
            let mut leg_ref = ReferenceLegacy::default();
            get_reference(&bed, 0, &chr_lengths, &cfg, &mut new_ref).unwrap();
            get_reference_legacy(&bed, 0, &chr_lengths, &cfg, &mut leg_ref).unwrap();

            let cap = new_ref.seed.materialize().capacity();
            assert!(
                cap <= 16,
                "all-'N' FASTA capacity must remain tiny (<= 16 <= 262_144, got {cap})"
            );
            assert_eq!(new_ref.seed.len(), leg_ref.seed.len());
            for (k, leg_pos) in &leg_ref.seed {
                assert_eq!(new_ref.seed.get(k.as_slice()), Some(leg_pos.as_slice()));
            }
        }

        // 3. Non-'N' stretches pre-reserve 2 * non_n up to the 262_144 cap:
        //    - 1_000 'A's (only 2 distinct keys) must pre-reserve >= 2_000 and <= 262_144.
        //    - Two 200_000 'A' chunks (2 * 400_000 = 800_000 uncapped) must cap reserve at 262_144
        //      (hashbrown bucket table <= 524_288 -> capacity <= 458_752, vs 917_504 uncapped).
        {
            let mut sm_small = SeedMap::new();
            sm_small.enqueue_chunk(PendingSeedChunk {
                bases: std::sync::Arc::from(vec![b'A'; 1_000]),
                sequence_start: 1,
                sequence_end: 1_000,
                chr_len: 1_000,
                prev_bounds: Vec::new(),
            });
            let cap_small = sm_small.materialize().capacity();
            assert!(
                (2_000..=262_144).contains(&cap_small),
                "1000 non-'N' bases must pre-reserve [2000, 262144], got {cap_small}"
            );

            let mut sm_capped = SeedMap::new();
            sm_capped.enqueue_chunk(PendingSeedChunk {
                bases: std::sync::Arc::from(vec![b'A'; 200_000]),
                sequence_start: 1,
                sequence_end: 200_000,
                chr_len: 400_000,
                prev_bounds: Vec::new(),
            });
            sm_capped.enqueue_chunk(PendingSeedChunk {
                bases: std::sync::Arc::from(vec![b'C'; 200_000]),
                sequence_start: 200_001,
                sequence_end: 400_000,
                chr_len: 400_000,
                prev_bounds: vec![(1, 200_000 - 17)],
            });
            let mut baseline_262k = RawSeedMap::default();
            baseline_262k.reserve(262_144);
            let cap_capped = sm_capped.materialize().capacity();
            assert_eq!(
                cap_capped,
                baseline_262k.capacity(),
                "reserve on >262_144 estimated k-mers must be capped at reserve(262_144)"
            );
        }
    }

    #[test]
    fn test_fetch_subsequence_and_get_reference_hla_colon_contig_and_u32_saturation() {
        let dir = tempfile::tempdir().unwrap();
        let fa_path = dir.path().join("hla.fa");
        let fai_path = dir.path().join("hla.fa.fai");
        let contig = "HLA-A*01:01:01:01";
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let header_line = format!(">{contig}\n");
        let mut fa = File::create(&fa_path).unwrap();
        fa.write_all(header_line.as_bytes()).unwrap();
        fa.write_all(seq).unwrap();
        fa.write_all(b"\n").unwrap();

        let mut fai = File::create(&fai_path).unwrap();
        writeln!(
            fai,
            "{contig}\t{}\t{}\t{}\t{}",
            seq.len(),
            header_line.len(),
            seq.len(),
            seq.len() + 1
        )
        .unwrap();

        let sub = fetch_subsequence(&fa_path, contig, Pos1(5), Pos1(16)).unwrap();
        assert_eq!(sub, &seq[4..16]);

        // Inverted range returns clean Parse error rather than panicking inside noodles-fasta
        assert!(fetch_subsequence(&fa_path, contig, Pos1(20), Pos1(10)).is_err());

        // get_reference saturates > u32::MAX chromosomes_lengths entry instead of truncating to 0/small u32
        let mut cfg = Config::default();
        cfg.reference = fa_path.to_str().unwrap().to_string();
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert(contig.to_string(), (u32::MAX as usize).saturating_add(5));
        let bed = BedRegion {
            chr: contig.to_string(),
            start: Pos1(1),
            end: Pos1(20),
            gene: "HLA-A".to_string(),
            amplicon: None,
        };
        let mut reference = Reference::new();
        get_reference(&bed, 0, &chr_lengths, &cfg, &mut reference).unwrap();
        assert_eq!(reference.get_base(contig, 1), Some(b'A'));

        // SeedMap::index_chunk saturates near u32::MAX without debug arithmetic panic
        let mut sm = SeedMap::new();
        sm.enqueue_chunk(PendingSeedChunk {
            bases: std::sync::Arc::from(&seq[..]),
            sequence_start: u32::MAX - 5,
            sequence_end: u32::MAX,
            chr_len: u32::MAX,
            prev_bounds: Vec::new(),
        });
        assert!(!sm.is_empty());
    }

    #[test]
    fn bgzf_fasta_short_uncompressed_hostile_gzi_and_fai_fallbacks() {
        clear_fasta_reader_cache();
        let dir = tempfile::tempdir().unwrap();

        // 1. Multi-block BGZF FASTA + .fai + .gzi matches uncompressed FASTA byte-for-byte
        let plain_fa = dir.path().join("ref.fa");
        let plain_fai = dir.path().join("ref.fa.fai");
        let bgzf_fa = dir.path().join("ref.fa.gz");
        let bgzf_fai = dir.path().join("ref.fa.gz.fai");
        let bgzf_gzi = dir.path().join("ref.fa.gz.gzi");

        let line1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bp
        let line2 = b"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"; // 40 bp
        let mut fasta_bytes = Vec::new();
        fasta_bytes.extend_from_slice(b">chr1\n");
        fasta_bytes.extend_from_slice(line1);
        fasta_bytes.push(b'\n');
        let split_upos = fasta_bytes.len();
        fasta_bytes.extend_from_slice(line2);
        fasta_bytes.push(b'\n');

        fs::write(&plain_fa, &fasta_bytes).unwrap();
        let fai_text = b"chr1\t80\t6\t40\t41\n";
        fs::write(&plain_fai, fai_text).unwrap();
        fs::write(&bgzf_fai, fai_text).unwrap();

        let mut bgzf_writer = noodles_bgzf::io::Writer::new(File::create(&bgzf_fa).unwrap());
        bgzf_writer.write_all(&fasta_bytes[..split_upos]).unwrap();
        bgzf_writer.flush().unwrap();
        let block1_cpos = bgzf_writer.position();
        bgzf_writer.write_all(&fasta_bytes[split_upos..]).unwrap();
        bgzf_writer.finish().unwrap();

        let gzi_index = noodles_bgzf::gzi::Index::from(vec![(block1_cpos, split_upos as u64)]);
        noodles_bgzf::gzi::fs::write(&bgzf_gzi, &gzi_index).unwrap();

        let sub_plain = fetch_subsequence(&plain_fa, "chr1", Pos1(1), Pos1(80)).unwrap();
        let sub_bgzf = fetch_subsequence(&bgzf_fa, "chr1", Pos1(1), Pos1(80)).unwrap();
        assert_eq!(sub_bgzf, sub_plain);
        let sub_cross = fetch_subsequence(&bgzf_fa, "chr1", Pos1(35), Pos1(55)).unwrap();
        assert_eq!(sub_cross, &sub_plain[34..55]);

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 80);
        let bed = BedRegion {
            chr: "chr1".to_string(),
            start: Pos1(10),
            end: Pos1(70),
            gene: "GENE1".to_string(),
            amplicon: None,
        };
        let mut cfg_plain = Config::default();
        cfg_plain.reference = plain_fa.to_str().unwrap().to_string();
        let mut cfg_bgzf = Config::default();
        cfg_bgzf.reference = bgzf_fa.to_str().unwrap().to_string();

        let mut ref_plain = Reference::new();
        let mut ref_bgzf = Reference::new();
        get_reference(&bed, 0, &chr_lengths, &cfg_plain, &mut ref_plain).unwrap();
        get_reference(&bed, 0, &chr_lengths, &cfg_bgzf, &mut ref_bgzf).unwrap();
        for pos in 1..=80u32 {
            assert_eq!(
                ref_bgzf.get_base("chr1", pos),
                ref_plain.get_base("chr1", pos),
                "base mismatch at pos {pos}"
            );
        }

        // 2. Tiny < 18-byte uncompressed FASTA (">1\nACGT\n" = 8 bytes) opens and queries cleanly
        let tiny_fa = dir.path().join("tiny.fa");
        let tiny_fai = dir.path().join("tiny.fa.fai");
        fs::write(&tiny_fa, b">1\nACGT\n").unwrap();
        fs::write(&tiny_fai, b"1\t4\t3\t4\t5\n").unwrap();
        assert!(open_indexed_fasta_buf_reader(&tiny_fa).is_ok());
        assert_eq!(
            fetch_subsequence(&tiny_fa, "1", Pos1(1), Pos1(4)).unwrap(),
            b"ACGT"
        );

        // 3. Plain non-BGZF gzip FASTA returns Err(VardictError::Parse(msg)) mentioning bgzip & samtools faidx
        let plain_gz_fa = dir.path().join("plain_gzip.fa.gz");
        {
            let mut gz = flate2::write::GzEncoder::new(
                File::create(&plain_gz_fa).unwrap(),
                flate2::Compression::default(),
            );
            gz.write_all(&fasta_bytes).unwrap();
            gz.finish().unwrap();
        }
        fs::write(dir.path().join("plain_gzip.fa.gz.fai"), fai_text).unwrap();
        match open_indexed_fasta_buf_reader(&plain_gz_fa) {
            Err(VardictError::Parse(msg)) => {
                assert!(msg.contains("bgzip"), "expected 'bgzip' in '{msg}'");
                assert!(
                    msg.contains("samtools faidx"),
                    "expected 'samtools faidx' in '{msg}'"
                );
            }
            Ok(_) => panic!("expected VardictError::Parse for plain gzip FASTA, got Ok"),
            Err(other) => {
                panic!("expected VardictError::Parse for plain gzip FASTA, got {other:?}")
            }
        }

        // 4. Hostile 8-byte [0xff; 8] .gzi (plus malformed/missing .gzi) returns Err(VardictError::Parse(_)) without panicking
        let hostile_fa = dir.path().join("hostile.fa.gz");
        let hostile_gzi = dir.path().join("hostile.fa.gz.gzi");
        fs::copy(&bgzf_fa, &hostile_fa).unwrap();
        // Missing .gzi
        assert!(matches!(
            open_indexed_fasta_buf_reader(&hostile_fa),
            Err(VardictError::Parse(_))
        ));
        // Hostile [0xff; 8] .gzi (AR-1 capacity-overflow bomb)
        fs::write(&hostile_gzi, [0xffu8; 8]).unwrap();
        assert!(matches!(
            open_indexed_fasta_buf_reader(&hostile_fa),
            Err(VardictError::Parse(_))
        ));
        // Misaligned 12-byte .gzi
        fs::write(&hostile_gzi, [0u8; 12]).unwrap();
        assert!(matches!(
            open_indexed_fasta_buf_reader(&hostile_fa),
            Err(VardictError::Parse(_))
        ));
        // Count mismatch: declares 2 entries in an 8-byte file
        fs::write(&hostile_gzi, 2u64.to_le_bytes()).unwrap();
        assert!(matches!(
            open_indexed_fasta_buf_reader(&hostile_fa),
            Err(VardictError::Parse(_))
        ));

        let hostile_fa_len = fs::metadata(&hostile_fa).unwrap().len();
        let write_gzi_pairs = |pairs: &[(u64, u64)]| -> Vec<u8> {
            let mut buf = Vec::with_capacity(8 + pairs.len() * 16);
            buf.extend_from_slice(&(pairs.len() as u64).to_le_bytes());
            for &(c, u) in pairs {
                buf.extend_from_slice(&c.to_le_bytes());
                buf.extend_from_slice(&u.to_le_bytes());
            }
            buf
        };
        // Unsorted / duplicate compressed or uncompressed offsets in .gzi must fail at open time
        for bad_pairs in [
            &[(0u64, 50u64)][..],
            &[(50u64, 0u64)][..],
            &[(10u64, 20u64), (10u64, 40u64)][..],
            &[(20u64, 20u64), (15u64, 40u64)][..],
            &[(10u64, 40u64), (20u64, 40u64)][..],
            &[(hostile_fa_len, 20u64)][..],
            &[(hostile_fa_len + 64, 20u64)][..],
        ] {
            fs::write(&hostile_gzi, write_gzi_pairs(bad_pairs)).unwrap();
            assert!(
                matches!(
                    open_indexed_fasta_buf_reader(&hostile_fa),
                    Err(VardictError::Parse(_))
                ),
                "expected Parse error for invalid .gzi pairs {bad_pairs:?}"
            );
        }

        // 5. Extensionless `genome` + `genome.fai` (never `genome..fai`) and `panel.fasta` + `panel.fai`
        let extless_fa = dir.path().join("genome");
        let extless_fai = dir.path().join("genome.fai");
        assert_eq!(resolve_fasta_fai_path(&extless_fa), extless_fai);
        fs::write(&extless_fa, b">1\nACGT\n").unwrap();
        fs::write(&extless_fai, b"1\t4\t3\t4\t5\n").unwrap();
        assert_eq!(resolve_fasta_fai_path(&extless_fa), extless_fai);
        assert_eq!(
            fetch_subsequence(&extless_fa, "1", Pos1(1), Pos1(4)).unwrap(),
            b"ACGT"
        );

        let stem_fa = dir.path().join("panel.fasta");
        let stem_fai = dir.path().join("panel.fai");
        fs::write(&stem_fa, b">1\nACGT\n").unwrap();
        fs::write(&stem_fai, b"1\t4\t3\t4\t5\n").unwrap();
        assert_eq!(resolve_fasta_fai_path(&stem_fa), stem_fai);
        assert_eq!(
            fetch_subsequence(&stem_fa, "1", Pos1(1), Pos1(4)).unwrap(),
            b"ACGT"
        );

        #[cfg(unix)]
        {
            // Dangling primary `.fai` symlink (`dangling.fasta.fai -> missing`) must NOT be silently
            // masked by `<stem>.fai` (`dangling.fai`): select primary and fail loudly (`AR-15`).
            let dang_fa = dir.path().join("dangling.fasta");
            let dang_primary_fai = dir.path().join("dangling.fasta.fai");
            let dang_stem_fai = dir.path().join("dangling.fai");
            fs::write(&dang_fa, b">1\nACGT\n").unwrap();
            fs::write(&dang_stem_fai, b"1\t4\t3\t4\t5\n").unwrap();
            std::os::unix::fs::symlink(
                dir.path().join("nonexistent_target.fai"),
                &dang_primary_fai,
            )
            .unwrap();
            assert_eq!(resolve_fasta_fai_path(&dang_fa), dang_primary_fai);
            assert!(
                fetch_subsequence(&dang_fa, "1", Pos1(1), Pos1(4)).is_err(),
                "dangling primary .fai symlink must fail loudly rather than falling back to <stem>.fai"
            );

            // Symlinked uncompressed and BGZF FASTA (`links/link.fa` and `links/link.fa.gz`) resolve
            // `.fai` and `.gzi` beside the canonical target directory when local companions are absent.
            let real_dir = dir.path().join("real_refs");
            let link_dir = dir.path().join("link_refs");
            fs::create_dir_all(&real_dir).unwrap();
            fs::create_dir_all(&link_dir).unwrap();

            let real_plain = real_dir.join("canon_panel.fasta");
            let real_plain_stem_fai = real_dir.join("canon_panel.fai");
            fs::write(&real_plain, b">1\nACGT\n").unwrap();
            fs::write(&real_plain_stem_fai, b"1\t4\t3\t4\t5\n").unwrap();
            let link_plain = link_dir.join("link.fa");
            std::os::unix::fs::symlink(Path::new("../real_refs/canon_panel.fasta"), &link_plain)
                .unwrap();
            assert_eq!(
                resolve_fasta_fai_path(&link_plain),
                fs::canonicalize(&real_plain_stem_fai).unwrap()
            );
            assert_eq!(
                fetch_subsequence(&link_plain, "1", Pos1(1), Pos1(4)).unwrap(),
                b"ACGT"
            );

            let real_bgzf = real_dir.join("canon_ref.fa.gz");
            let real_bgzf_fai = real_dir.join("canon_ref.fa.gz.fai");
            let real_bgzf_gzi = real_dir.join("canon_ref.fa.gz.gzi");
            fs::copy(&bgzf_fa, &real_bgzf).unwrap();
            fs::copy(&bgzf_fai, &real_bgzf_fai).unwrap();
            fs::copy(&bgzf_gzi, &real_bgzf_gzi).unwrap();
            let link_bgzf = link_dir.join("link.fa.gz");
            std::os::unix::fs::symlink(Path::new("../real_refs/canon_ref.fa.gz"), &link_bgzf)
                .unwrap();
            assert_eq!(
                resolve_fasta_fai_path(&link_bgzf),
                fs::canonicalize(&real_bgzf_fai).unwrap()
            );
            assert_eq!(
                resolve_fasta_gzi_path(&link_bgzf),
                fs::canonicalize(&real_bgzf_gzi).unwrap()
            );
            assert_eq!(
                fetch_subsequence(&link_bgzf, "chr1", Pos1(1), Pos1(10)).unwrap(),
                fetch_subsequence(&bgzf_fa, "chr1", Pos1(1), Pos1(10)).unwrap()
            );
        }

        // Compressed `.gz` / `.bgz` must NOT fall back to uncompressed `<stem>.fai`
        let gz_no_fai = dir.path().join("panel.gz");
        assert_eq!(
            resolve_fasta_fai_path(&gz_no_fai),
            dir.path().join("panel.gz.fai")
        );
        assert_eq!(
            resolve_fasta_gzi_path(&bgzf_fa),
            dir.path().join("ref.fa.gz.gzi")
        );
    }
}
