use crate::types::{Pos0, Pos1};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::alignment::record::{QualityScores, Sequence};
use noodles::sam::alignment::Record;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SamRecord {
    pub name: Vec<u8>,
    pub flags: u16,
    pub reference_name: Option<String>,
    pub alignment_start: Option<Pos1>,
    pub alignment_end: Option<Pos1>,
    pub mapping_quality: Option<u8>,
    pub cigar: Vec<(u32, CigarKind)>,
    pub sequence: Vec<u8>,
    pub quality_scores: Vec<u8>,
    pub mate_reference_name: Option<String>,
    pub mate_alignment_start: Option<Pos1>,
    pub template_length: i32,
    pub nm: Option<i32>,
    pub sa: Option<String>,
    pub mc: Option<String>,
    pub mq: Option<u8>,
}

impl SamRecord {
    /// Empty record with no heap allocation; use with
    /// [`SamRecord::fill_from_alignment_record`] for buffer reuse.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset every field so a reused record cannot leak a previous read's
    /// values into the next fill. String buffers keep their capacity
    /// (`Some("")`, not `None`); `fill_from_alignment_record` assigns every
    /// field, setting them back to `None` when absent. Must stay in sync
    /// with the field list and the fill function.
    pub fn clear(&mut self) {
        self.name.clear();
        self.flags = 0;
        if let Some(s) = self.reference_name.as_mut() {
            s.clear();
        }
        self.alignment_start = None;
        self.alignment_end = None;
        self.mapping_quality = None;
        self.cigar.clear();
        self.sequence.clear();
        self.quality_scores.clear();
        if let Some(s) = self.mate_reference_name.as_mut() {
            s.clear();
        }
        self.mate_alignment_start = None;
        self.template_length = 0;
        self.nm = None;
        if let Some(s) = self.sa.as_mut() {
            s.clear();
        }
        if let Some(s) = self.mc.as_mut() {
            s.clear();
        }
        self.mq = None;
    }

    /// Refill `self` from `src`, reusing `Vec` and `String` capacities in place
    /// while preserving `None` on absent `Option<String>` fields (`C-1`).
    pub fn clone_from_record(&mut self, src: &Self) {
        fn clone_opt_string(dst: &mut Option<String>, src: &Option<String>) {
            match src {
                Some(s) => {
                    let buf = dst.get_or_insert_with(String::new);
                    buf.clear();
                    buf.push_str(s);
                }
                None => *dst = None,
            }
        }

        let SamRecord {
            name,
            flags,
            reference_name,
            alignment_start,
            alignment_end,
            mapping_quality,
            cigar,
            sequence,
            quality_scores,
            mate_reference_name,
            mate_alignment_start,
            template_length,
            nm,
            sa,
            mc,
            mq,
        } = src;
        self.name.clear();
        self.name.extend_from_slice(name);
        self.flags = *flags;
        clone_opt_string(&mut self.reference_name, reference_name);
        self.alignment_start = *alignment_start;
        self.alignment_end = *alignment_end;
        self.mapping_quality = *mapping_quality;
        self.cigar.clear();
        self.cigar.extend_from_slice(cigar);
        self.sequence.clear();
        self.sequence.extend_from_slice(sequence);
        self.quality_scores.clear();
        self.quality_scores.extend_from_slice(quality_scores);
        clone_opt_string(&mut self.mate_reference_name, mate_reference_name);
        self.mate_alignment_start = *mate_alignment_start;
        self.template_length = *template_length;
        self.nm = *nm;
        clone_opt_string(&mut self.sa, sa);
        clone_opt_string(&mut self.mc, mc);
        self.mq = *mq;
        debug_assert_eq!(self, src);
    }

    /// Refill `self` in place from a noodles record, reusing heap buffers.
    /// Values are identical to [`SamRecord::from_alignment_record`]; every
    /// field is assigned unconditionally so no stale value can survive.
    pub fn fill_from_alignment_record<R>(
        &mut self,
        record: &R,
        header: &noodles::sam::Header,
    ) -> Result<(), std::io::Error>
    where
        R: Record,
    {
        self.clear();

        if let Some(name) = record.name() {
            let bytes: &[u8] = name.as_ref();
            self.name.extend_from_slice(bytes);
        }

        self.flags = record.flags().map(|f| f.bits()).unwrap_or(0);

        // Reuse the (capacity-retained) String buffer; push from the
        // borrowed Cow so valid-UTF-8 names allocate nothing.
        match record
            .reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| name)
            }) {
            Some(name) => {
                let buf = self.reference_name.get_or_insert_with(String::new);
                buf.clear();
                buf.push_str(&String::from_utf8_lossy(name));
            }
            None => self.reference_name = None,
        }

        self.alignment_start = if let Some(pos) = record.alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.alignment_end = if let Some(pos) = record.alignment_end() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.mapping_quality = if (self.flags & 0x4) != 0 {
            None
        } else {
            Some(
                record
                    .mapping_quality()
                    .and_then(|m| m.ok())
                    .map(|m| m.get())
                    .unwrap_or(255),
            )
        };

        let cigar = record.cigar();
        self.cigar.clear();
        self.cigar.reserve(cigar.len());
        for op in cigar.iter().flatten() {
            self.cigar.push((op.len() as u32, op.kind()));
        }

        let seq = record.sequence();
        self.sequence.clear();
        self.sequence.reserve(seq.len());
        for b in seq.iter() {
            // Normalize reference-matching CRAM bases decoded against soft-masked
            // lowercase FASTA contigs (`a/c/g/t/n`) to uppercase (`A/C/G/T/N`),
            // matching BAM 4-bit nibble decoding and HTSJDK `SAMRecord.getReadBases()`.
            self.sequence.push(b.to_ascii_uppercase());
        }

        let quals = record.quality_scores();
        self.quality_scores.clear();
        self.quality_scores.reserve(quals.len());
        for score in quals.iter().flatten() {
            self.quality_scores.push(score.saturating_add(33));
        }

        match record
            .mate_reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| name)
            }) {
            Some(name) => {
                let buf = self.mate_reference_name.get_or_insert_with(String::new);
                buf.clear();
                buf.push_str(&String::from_utf8_lossy(name));
            }
            None => self.mate_reference_name = None,
        }

        self.mate_alignment_start = if let Some(pos) = record.mate_alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.template_length = record.template_length().unwrap_or(0);

        let data = record.data();
        self.fill_aux_from_iter(data.iter());

        Ok(())
    }

    /// Concrete `noodles::bam::Record` refill path (O-3): uses inherent methods
    /// on `noodles::bam::Record` and a single pass over `record.data().iter()`,
    /// avoiding per-read `Box<dyn Cigar>` / `Box<dyn Sequence>` /
    /// `Box<dyn QualityScores>` / `Box<dyn Data>` trait-object allocations and
    /// repeated linear scans of the aux block.
    pub fn fill_from_bam_record(
        &mut self,
        record: &noodles::bam::Record,
        header: &noodles::sam::Header,
    ) -> Result<(), std::io::Error> {
        self.clear();

        if let Some(name) = record.name() {
            let bytes: &[u8] = name.as_ref();
            self.name.extend_from_slice(bytes);
        }

        self.flags = record.flags().bits();

        match record
            .reference_sequence_id()
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| name)
            }) {
            Some(name) => {
                let buf = self.reference_name.get_or_insert_with(String::new);
                buf.clear();
                buf.push_str(&String::from_utf8_lossy(name));
            }
            None => self.reference_name = None,
        }

        self.alignment_start = if let Some(pos) = record.alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.alignment_end = if let Some(pos) = bam_alignment_end(record) {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.mapping_quality = if (self.flags & 0x4) != 0 {
            None
        } else {
            Some(record.mapping_quality().map(|m| m.get()).unwrap_or(255))
        };

        let cigar = record.cigar();
        self.cigar.clear();
        self.cigar.reserve(cigar.len());
        for op in cigar.iter().flatten() {
            self.cigar.push((op.len() as u32, op.kind()));
        }

        let seq = record.sequence();
        self.sequence.clear();
        self.sequence.reserve(seq.len());
        for b in seq.iter() {
            self.sequence.push(b);
        }

        let quals = record.quality_scores();
        let raw_quals: &[u8] = quals.as_ref();
        self.quality_scores.clear();
        self.quality_scores.reserve(raw_quals.len());
        for &score in raw_quals {
            self.quality_scores.push(score.saturating_add(33));
        }

        match record
            .mate_reference_sequence_id()
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| name)
            }) {
            Some(name) => {
                let buf = self.mate_reference_name.get_or_insert_with(String::new);
                buf.clear();
                buf.push_str(&String::from_utf8_lossy(name));
            }
            None => self.mate_reference_name = None,
        }

        self.mate_alignment_start = if let Some(pos) = record.mate_alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.template_length = record.template_length();

        self.fill_aux_from_iter(record.data().iter());

        Ok(())
    }

    fn fill_aux_from_iter<'a, I>(&mut self, iter: I)
    where
        I: Iterator<
            Item = Result<
                (
                    noodles::sam::alignment::record::data::field::Tag,
                    noodles::sam::alignment::record::data::field::Value<'a>,
                ),
                std::io::Error,
            >,
        >,
    {
        use noodles::sam::alignment::record::data::field::{Tag, Value};

        const TAG_NM_LOWER: Tag = Tag::new(b'n', b'M');
        const TAG_NM_UPPER: Tag = Tag::EDIT_DISTANCE;
        const TAG_SA: Tag = Tag::new(b'S', b'A');
        const TAG_MC: Tag = Tag::new(b'M', b'C');
        const TAG_MQ: Tag = Tag::new(b'M', b'Q');

        let mut seen_nm_lower = false;
        let mut seen_nm_upper = false;
        let mut seen_sa = false;
        let mut seen_mc = false;
        let mut seen_mq = false;
        let mut nm_lower_val: Option<i32> = None;
        let mut nm_upper_val: Option<i32> = None;

        for field_res in iter {
            let Ok((tag, val)) = field_res else {
                break;
            };
            if tag == TAG_NM_LOWER {
                if !seen_nm_lower {
                    seen_nm_lower = true;
                    nm_lower_val = val.as_int().map(|x| x as i32);
                }
            } else if tag == TAG_NM_UPPER {
                if !seen_nm_upper {
                    seen_nm_upper = true;
                    nm_upper_val = val.as_int().map(|x| x as i32);
                }
            } else if tag == TAG_SA {
                if !seen_sa {
                    seen_sa = true;
                    let buf = self.sa.get_or_insert_with(String::new);
                    buf.clear();
                    match val {
                        Value::String(s) => {
                            buf.push_str(&String::from_utf8_lossy(s.as_ref()));
                        }
                        _ => {
                            use std::fmt::Write as _;
                            let _ = write!(buf, "{val:?}");
                        }
                    }
                }
            } else if tag == TAG_MC {
                if !seen_mc {
                    seen_mc = true;
                    let buf = self.mc.get_or_insert_with(String::new);
                    buf.clear();
                    match val {
                        Value::String(s) => {
                            buf.push_str(&String::from_utf8_lossy(s.as_ref()));
                        }
                        _ => {
                            use std::fmt::Write as _;
                            let _ = write!(buf, "{val:?}");
                        }
                    }
                }
            } else if tag == TAG_MQ && !seen_mq {
                seen_mq = true;
                self.mq = val.as_int().map(|x| x as u8);
            }
        }

        // Ports the original if/else-if chain exactly: EDIT_DISTANCE is only
        // consulted when the nM tag is absent or undecodable (a present but
        // non-int nM sets seen_nm_lower = true and suppresses the fallback).
        self.nm = if seen_nm_lower {
            nm_lower_val
        } else {
            nm_upper_val
        };
        if !seen_sa {
            self.sa = None;
        }
        if !seen_mc {
            self.mc = None;
        }
    }

    pub fn from_alignment_record<R>(
        record: &R,
        header: &noodles::sam::Header,
    ) -> Result<Self, std::io::Error>
    where
        R: Record,
    {
        let mut rec = Self::default();
        rec.fill_from_alignment_record(record, header)?;
        Ok(rec)
    }
}

/// Reconstruct the SAMtags §1.5 `NM` tag (edit distance to the reference,
/// excluding clipping and reference skips) for a mapped CRAM record when the
/// CRAM encoder omitted `NM` on write per CRAM 3.0 §10.6 (the `samtools view -C`
/// default `store_nm=0`).
pub(crate) fn compute_cram_nm(rec: &SamRecord, ref_seq: &[u8]) -> Option<i32> {
    if (rec.flags & 0x4) != 0 || rec.cigar.is_empty() || rec.sequence.is_empty() {
        return None;
    }
    let start_1 = rec.alignment_start?.0;
    if start_1 == 0 {
        return None;
    }
    let mut ref_pos = (start_1 - 1) as usize;
    let mut read_pos = 0usize;
    let mut nm: i32 = 0;
    for &(len, op) in &rec.cigar {
        let len = len as usize;
        match op {
            CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                for _ in 0..len {
                    let Some(&read_b) = rec.sequence.get(read_pos) else {
                        break;
                    };
                    let ref_b = ref_seq.get(ref_pos).copied().unwrap_or(b'N');
                    let is_match = read_b == b'='
                        || (matches!(
                            read_b,
                            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
                        ) && read_b.eq_ignore_ascii_case(&ref_b));
                    if !is_match {
                        nm = nm.saturating_add(1);
                    }
                    read_pos = read_pos.saturating_add(1);
                    ref_pos = ref_pos.saturating_add(1);
                }
            }
            CigarKind::Insertion => {
                nm = nm.saturating_add(i32::try_from(len).unwrap_or(i32::MAX));
                read_pos = read_pos.saturating_add(len);
            }
            CigarKind::Deletion => {
                nm = nm.saturating_add(i32::try_from(len).unwrap_or(i32::MAX));
                ref_pos = ref_pos.saturating_add(len);
            }
            CigarKind::Skip => {
                ref_pos = ref_pos.saturating_add(len);
            }
            CigarKind::SoftClip => {
                read_pos = read_pos.saturating_add(len);
            }
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }
    Some(nm)
}

pub enum IndexedSamReader {
    Bam(noodles::bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>),
    Cram(
        noodles::cram::io::IndexedReader<std::fs::File>,
        noodles::fasta::Repository,
    ),
}

impl IndexedSamReader {
    pub fn query<'a>(
        &'a mut self,
        header: &'a noodles::sam::Header,
        region: &'a noodles::core::Region,
    ) -> Result<Box<dyn Iterator<Item = Result<SamRecord, std::io::Error>> + 'a>, std::io::Error>
    {
        match self {
            IndexedSamReader::Bam(reader) => {
                let query = reader.query(header, region)?;
                let iter = query.records().map(|record| {
                    let rec = record?;
                    SamRecord::from_alignment_record(&rec, header)
                });
                Ok(Box::new(iter))
            }
            IndexedSamReader::Cram(reader, repository) => {
                let query = reader.query(header, region)?;
                let region_name = region.name().to_owned();
                let repo = repository.clone();
                let mut ref_seq: Option<Option<std::sync::Arc<noodles::fasta::record::Sequence>>> =
                    None;
                let iter = query.records().map(move |record| {
                    let rec = record?;
                    let mut sam_rec = SamRecord::from_alignment_record(&rec, header)?;
                    if sam_rec.nm.is_none() {
                        let seq_opt = ref_seq.get_or_insert_with(|| {
                            repo.get(region_name.as_ref()).and_then(|r| r.ok())
                        });
                        if let Some(ref seq) = seq_opt {
                            sam_rec.nm = compute_cram_nm(&sam_rec, seq.as_ref().as_ref());
                        }
                    }
                    Ok(sam_rec)
                });
                Ok(Box::new(iter))
            }
        }
    }
}

pub struct CachedReader {
    pub reader: IndexedSamReader,
    pub header: noodles::sam::Header,
}

thread_local! {
    static SAM_READER_CACHE: RefCell<HashMap<PathBuf, CachedReader>> = RefCell::new(HashMap::new());
}

static REPOSITORY_CACHE: std::sync::Mutex<
    Option<(std::path::PathBuf, noodles::fasta::Repository)>,
> = std::sync::Mutex::new(None);

pub fn open_reader(
    path: &Path,
    config: &crate::config::Config,
) -> Result<CachedReader, crate::types::VardictError> {
    let path_str = path.to_string_lossy();
    if path_str.ends_with(".cram") || path_str.ends_with(".CRAM") {
        let fasta_path = std::path::Path::new(&config.reference);
        let mut fai_path = fasta_path.to_path_buf();
        let mut ext = fai_path
            .extension()
            .map(|e| e.to_string_lossy().into_owned())
            .unwrap_or_default();
        ext.push_str(".fai");
        fai_path.set_extension(ext);

        let repository = {
            let mut guard = REPOSITORY_CACHE.lock().unwrap_or_else(|e| e.into_inner());
            match guard.as_ref() {
                Some((cached_path, repo)) if cached_path == fasta_path => repo.clone(),
                _ => {
                    let fasta_reader = noodles::fasta::io::IndexedReader::new(
                        std::io::BufReader::new(
                            std::fs::File::open(fasta_path)
                                .map_err(crate::types::VardictError::Io)?,
                        ),
                        noodles::fasta::fai::fs::read(&fai_path)
                            .map_err(crate::types::VardictError::Io)?,
                    );
                    use noodles::fasta::repository::adapters::IndexedReader as FastaAdapter;
                    use noodles::fasta::Repository;
                    let repo = Repository::new(FastaAdapter::new(fasta_reader));
                    *guard = Some((fasta_path.to_path_buf(), repo.clone()));
                    repo
                }
            }
        };

        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_path(path)
            .map_err(crate::types::VardictError::Io)?;
        let header = reader
            .read_header()
            .map_err(crate::types::VardictError::Io)?;
        Ok(CachedReader {
            reader: IndexedSamReader::Cram(reader, repository),
            header,
        })
    } else {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(crate::types::VardictError::Io)?;
        let header = reader
            .read_header()
            .map_err(crate::types::VardictError::Io)?;
        Ok(CachedReader {
            reader: IndexedSamReader::Bam(reader),
            header,
        })
    }
}

pub fn query_region(
    bam_path: &Path,
    chr: &str,
    start: Pos0,
    end: Pos0,
    config: &crate::config::Config,
) -> Result<Vec<SamRecord>, crate::types::VardictError> {
    SAM_READER_CACHE.with(
        |cache| -> Result<Vec<SamRecord>, crate::types::VardictError> {
            let mut cache = cache.borrow_mut();
            if !cache.contains_key(bam_path) {
                let reader = open_reader(bam_path, config)?;
                cache.insert(bam_path.to_path_buf(), reader);
            }
            let cached = match cache.get_mut(bam_path) {
                Some(c) => c,
                None => {
                    return Err(crate::types::VardictError::Other(format!(
                        "BAM cache miss for {}",
                        bam_path.display()
                    )))
                }
            };

            let start_pos = start.to_one_based().0;
            let end_pos = end.0;
            if start_pos > end_pos {
                return Ok(Vec::new());
            }

            let region_str = format!("{}:{}-{}", chr, start_pos, end_pos);
            let region: noodles::core::Region = region_str.parse().map_err(|e| {
                crate::types::VardictError::Parse(format!(
                    "failed to parse region '{}': {}",
                    region_str, e
                ))
            })?;

            let iter = cached
                .reader
                .query(&cached.header, &region)
                .map_err(crate::types::VardictError::Io)?;
            let mut records = Vec::new();
            for record in iter {
                records.push(record.map_err(crate::types::VardictError::Io)?);
            }
            Ok(records)
        },
    )
}

static NEXT_BAM_META_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
/// Synthetic `meta_id` floor reserved for isolated unit tests so concurrent
/// `clear_bam_caches()` calls in parallel tests never evict in-flight test blocks (`AR-4`).
#[cfg(test)]
const TEST_META_ID_BASE: u64 = 1u64 << 63;

struct SharedBamMeta {
    id: u64,
    header: noodles::sam::Header,
    index: noodles::bam::Index,
}

const BGZF_CACHE_SHARDS: usize = 8;
const BGZF_CACHE_SLOTS_PER_SHARD: usize = 16;

#[derive(Clone, Debug)]
struct CachedBgzfBlock {
    block_cpos: u64,
    next_cpos: u64,
    data: std::sync::Arc<[u8]>,
}

impl Default for CachedBgzfBlock {
    fn default() -> Self {
        Self {
            block_cpos: 0,
            next_cpos: 0,
            data: std::sync::Arc::from(&[] as &[u8]),
        }
    }
}

#[derive(Default)]
struct BgzfCacheShard {
    entries: Vec<((u64, u64), CachedBgzfBlock)>,
}

impl BgzfCacheShard {
    fn get(&mut self, key: (u64, u64)) -> Option<CachedBgzfBlock> {
        let idx = self.entries.iter().position(|(k, _)| *k == key)?;
        let entry = self.entries.remove(idx);
        let block = entry.1.clone();
        self.entries.push(entry);
        Some(block)
    }

    fn insert(&mut self, key: (u64, u64), block: CachedBgzfBlock) {
        if let Some(idx) = self.entries.iter().position(|(k, _)| *k == key) {
            self.entries.remove(idx);
        } else if self.entries.len() >= BGZF_CACHE_SLOTS_PER_SHARD {
            self.entries.remove(0);
        }
        self.entries.push((key, block));
    }

    fn clear_production(&mut self) {
        #[cfg(test)]
        self.entries.retain(|((id, _), _)| *id >= TEST_META_ID_BASE);
        #[cfg(not(test))]
        self.entries.clear();
    }
}

static BGZF_BLOCK_CACHE: std::sync::OnceLock<
    [std::sync::Mutex<BgzfCacheShard>; BGZF_CACHE_SHARDS],
> = std::sync::OnceLock::new();

fn bgzf_cache_shards() -> &'static [std::sync::Mutex<BgzfCacheShard>; BGZF_CACHE_SHARDS] {
    BGZF_BLOCK_CACHE
        .get_or_init(|| std::array::from_fn(|_| std::sync::Mutex::new(BgzfCacheShard::default())))
}

#[inline]
fn bgzf_shard_index(meta_id: u64, cpos: u64) -> usize {
    let h = cpos
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(meta_id.wrapping_mul(0x517C_C1B7_2722_0A95));
    ((h ^ (h >> 32)) as usize) & (BGZF_CACHE_SHARDS - 1)
}

/// BGZF reader wrapper with a 1-block thread-local cursor for $O(1)$ same-block
/// seeks and a bounded process-wide sharded LRU block cache (`EBO-1`).
pub struct CachedBgzfReader {
    inner: noodles::bgzf::io::Reader<std::fs::File>,
    meta_id: u64,
    current: CachedBgzfBlock,
    upos: usize,
    at_eof: bool,
    /// Tracks the exact `next_cpos` where `self.inner` is positioned with a fully
    /// consumed buffer; cleared (`None`) before any `inner` I/O so a failed
    /// read/seek always forces an explicit `inner.seek(vp0)` on subsequent access (`AR-5`).
    inner_next_cpos: Option<u64>,
    #[cfg(test)]
    local_hits: u64,
    #[cfg(test)]
    shared_hits: u64,
    #[cfg(test)]
    misses: u64,
}

impl CachedBgzfReader {
    pub fn new(file: std::fs::File, meta_id: u64) -> Self {
        Self {
            inner: noodles::bgzf::io::Reader::new(file),
            meta_id,
            current: CachedBgzfBlock::default(),
            upos: 0,
            at_eof: false,
            inner_next_cpos: Some(0),
            #[cfg(test)]
            local_hits: 0,
            #[cfg(test)]
            shared_hits: 0,
            #[cfg(test)]
            misses: 0,
        }
    }

    #[cfg(test)]
    pub(crate) fn reset_cache_stats(&mut self) {
        self.local_hits = 0;
        self.shared_hits = 0;
        self.misses = 0;
    }

    #[cfg(test)]
    pub(crate) fn cache_stats(&self) -> (u64, u64, u64) {
        (self.local_hits, self.shared_hits, self.misses)
    }

    #[inline]
    fn reset_empty_at(&mut self, cpos: u64) {
        self.current = CachedBgzfBlock {
            block_cpos: cpos,
            next_cpos: cpos,
            data: std::sync::Arc::from(&[] as &[u8]),
        };
        self.upos = 0;
        self.at_eof = false;
    }

    fn fetch_block(&mut self, cpos: u64) -> std::io::Result<CachedBgzfBlock> {
        use std::io::BufRead as _;

        let shard_idx = bgzf_shard_index(self.meta_id, cpos);
        if let Some(shard_mutex) = bgzf_cache_shards().get(shard_idx) {
            let mut shard = shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
            if let Some(block) = shard.get((self.meta_id, cpos)) {
                #[cfg(test)]
                {
                    self.shared_hits += 1;
                }
                return Ok(block);
            }
        }

        #[cfg(test)]
        {
            self.misses += 1;
        }

        // Only read sequentially without `seek` if the previous `fetch_block` left
        // `self.inner` cleanly consumed at `cpos` (`AR-5`).
        let can_read_seq = self.inner_next_cpos.take() == Some(cpos);
        if !can_read_seq {
            let vp0 = noodles::bgzf::VirtualPosition::try_from((cpos, 0))
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
            self.inner.seek(vp0)?;
        }

        let buf = self.inner.fill_buf()?;
        let buf_len = buf.len();
        let data: std::sync::Arc<[u8]> = std::sync::Arc::from(buf);
        let block_cpos = self.inner.virtual_position().compressed();
        let next_cpos = self.inner.position();
        self.inner.consume(buf_len);
        self.inner_next_cpos = Some(next_cpos);

        let block = CachedBgzfBlock {
            block_cpos,
            next_cpos,
            data,
        };

        if let Some(shard_mutex) = bgzf_cache_shards().get(shard_idx) {
            let mut shard = shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
            shard.insert((self.meta_id, cpos), block.clone());
        }
        if block_cpos != cpos && !block.data.is_empty() {
            let actual_shard = bgzf_shard_index(self.meta_id, block_cpos);
            if let Some(shard_mutex) = bgzf_cache_shards().get(actual_shard) {
                let mut shard = shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
                shard.insert((self.meta_id, block_cpos), block.clone());
            }
        }

        Ok(block)
    }

    #[inline]
    pub fn virtual_position(&self) -> noodles::bgzf::VirtualPosition {
        if self.upos < self.current.data.len() {
            noodles::bgzf::VirtualPosition::try_from((self.current.block_cpos, self.upos as u16))
                .unwrap_or_default()
        } else {
            noodles::bgzf::VirtualPosition::try_from((self.current.next_cpos, 0))
                .unwrap_or_default()
        }
    }
}

impl std::io::Read for CachedBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        use std::io::BufRead as _;
        let src = self.fill_buf()?;
        let amt = src.len().min(buf.len());
        if let (Some(dst_slice), Some(src_slice)) = (buf.get_mut(..amt), src.get(..amt)) {
            dst_slice.copy_from_slice(src_slice);
        }
        self.consume(amt);
        Ok(amt)
    }
}

impl std::io::BufRead for CachedBgzfReader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        if self.upos >= self.current.data.len() && !self.at_eof {
            let next_cpos = self.current.next_cpos;
            let block = self.fetch_block(next_cpos)?;
            self.at_eof = block.data.is_empty();
            self.current = block;
            self.upos = 0;
        }
        Ok(self.current.data.get(self.upos..).unwrap_or(&[]))
    }

    fn consume(&mut self, amt: usize) {
        self.upos = self.upos.saturating_add(amt).min(self.current.data.len());
    }
}

impl noodles::bgzf::io::Read for CachedBgzfReader {
    fn virtual_position(&self) -> noodles::bgzf::VirtualPosition {
        CachedBgzfReader::virtual_position(self)
    }
}

impl noodles::bgzf::io::BufRead for CachedBgzfReader {}

impl noodles::bgzf::io::Seek for CachedBgzfReader {
    fn seek_to_virtual_position(
        &mut self,
        pos: noodles::bgzf::VirtualPosition,
    ) -> std::io::Result<noodles::bgzf::VirtualPosition> {
        let (cpos, upos_u16) = pos.into();
        let upos = usize::from(upos_u16);

        // Fast path: target block is already loaded in our 1-block local cursor.
        if !self.current.data.is_empty()
            && self.current.block_cpos == cpos
            && cpos != self.current.next_cpos
        {
            if upos <= self.current.data.len() {
                #[cfg(test)]
                {
                    self.local_hits += 1;
                }
                self.upos = upos;
                self.at_eof = false;
                return Ok(pos);
            } else {
                // Match noodles-bgzf (`reader.rs:194-200`): retain loaded block with `upos = 0` (`AR-7`).
                self.upos = 0;
                self.at_eof = false;
                return Err(std::io::Error::from(std::io::ErrorKind::InvalidInput));
            }
        }

        let block = match self.fetch_block(cpos) {
            Ok(b) => b,
            Err(e) => {
                self.reset_empty_at(cpos);
                return Err(e);
            }
        };

        if upos > block.data.len() {
            self.at_eof = block.data.is_empty();
            self.current = block;
            self.upos = 0;
            return Err(std::io::Error::from(std::io::ErrorKind::InvalidInput));
        }

        self.at_eof = block.data.is_empty();
        self.current = block;
        self.upos = upos;
        Ok(pos)
    }

    fn seek_with_index(
        &mut self,
        index: &noodles::bgzf::gzi::Index,
        pos: std::io::SeekFrom,
    ) -> std::io::Result<u64> {
        match pos {
            std::io::SeekFrom::Start(upos) => {
                let vp = index.query(upos)?;
                self.seek_to_virtual_position(vp)?;
                Ok(upos)
            }
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "only SeekFrom::Start is supported by BGZF gzi index",
            )),
        }
    }
}

type RawBamReader = noodles::bam::io::Reader<CachedBgzfReader>;
type PathSlot = std::sync::Arc<std::sync::Mutex<Option<std::sync::Arc<SharedBamMeta>>>>;
type PooledReaderEntry = (std::sync::Arc<SharedBamMeta>, RawBamReader);
type BamReaderPoolMap = HashMap<PathBuf, Vec<PooledReaderEntry>>;

static META_CACHE: std::sync::OnceLock<std::sync::Mutex<HashMap<PathBuf, PathSlot>>> =
    std::sync::OnceLock::new();

const MAX_POOLED_READERS_PER_PATH: usize = 4;

thread_local! {
    static BAM_READER_POOL: RefCell<BamReaderPoolMap> = RefCell::new(HashMap::new());
}

/// Clear the process-wide shared BAM/CRAM metadata caches, bounded FASTA
/// repository, decoded block/container caches, and the calling thread's
/// pooled readers. Called at the start and end of [`crate::run`] via an
/// RAII drop guard (`F-5`).
pub fn clear_bam_caches() {
    if let Some(cache) = META_CACHE.get() {
        cache.lock().unwrap_or_else(|e| e.into_inner()).clear();
    }
    if let Some(shards) = BGZF_BLOCK_CACHE.get() {
        for shard in shards {
            shard
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .clear_production();
        }
    }
    let _ = BAM_READER_POOL.try_with(|pool| {
        if let Ok(mut map) = pool.try_borrow_mut() {
            map.clear();
        }
    });
    if let Some(cache) = CRAM_META_CACHE.get() {
        cache.lock().unwrap_or_else(|e| e.into_inner()).clear();
    }
    if let Some(shards) = CRAM_CONTAINER_CACHE.get() {
        for shard in shards {
            shard
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .clear_production();
        }
    }
    let _ = CRAM_READER_POOL.try_with(|pool| {
        if let Ok(mut map) = pool.try_borrow_mut() {
            map.clear();
        }
    });
    if let Some(cache) = FASTA_REPO_CACHE.get() {
        cache.lock().unwrap_or_else(|e| e.into_inner()).clear();
    }
    if let Some((_, repo)) = REPOSITORY_CACHE
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .take()
    {
        repo.clear();
    }
    let _ = SAM_READER_CACHE.try_with(|cache| {
        if let Ok(mut map) = cache.try_borrow_mut() {
            map.clear();
        }
    });
}

fn pop_pooled_reader(path: &Path, meta: &std::sync::Arc<SharedBamMeta>) -> Option<RawBamReader> {
    BAM_READER_POOL
        .try_with(|pool| {
            let mut map = pool.try_borrow_mut().ok()?;
            let vec = map.get_mut(path)?;
            while let Some((pooled_meta, reader)) = vec.pop() {
                if std::sync::Arc::ptr_eq(&pooled_meta, meta) {
                    return Some(reader);
                }
            }
            None
        })
        .ok()
        .flatten()
}

fn return_pooled_reader(path: PathBuf, meta: std::sync::Arc<SharedBamMeta>, reader: RawBamReader) {
    let _ = BAM_READER_POOL.try_with(|pool| {
        if let Ok(mut map) = pool.try_borrow_mut() {
            let vec = map.entry(path).or_default();
            vec.retain(|(m, _)| std::sync::Arc::ptr_eq(m, &meta));
            if vec.len() < MAX_POOLED_READERS_PER_PATH {
                vec.push((meta, reader));
            }
        }
    });
}

fn acquire_bam_meta_and_reader(
    path: &Path,
) -> Result<(std::sync::Arc<SharedBamMeta>, RawBamReader), crate::types::VardictError> {
    let slot = {
        let cache = META_CACHE.get_or_init(|| std::sync::Mutex::new(HashMap::new()));
        let mut map = cache.lock().unwrap_or_else(|e| e.into_inner());
        std::sync::Arc::clone(
            map.entry(path.to_path_buf())
                .or_insert_with(|| std::sync::Arc::new(std::sync::Mutex::new(None))),
        )
    };

    let meta = {
        let mut slot_guard = slot.lock().unwrap_or_else(|e| e.into_inner());
        match slot_guard.as_ref() {
            Some(meta) => std::sync::Arc::clone(meta),
            None => {
                // Exact same open sequence as `open_reader` (index load ->
                // File::open -> read_header). Errors return without storing
                // `Some(...)` so transient/test errors are never cached and
                // io::ErrorKind is preserved.
                let index = noodles::bam::fs::read_associated_index(path)
                    .map_err(crate::types::VardictError::Io)?;
                let file = std::fs::File::open(path).map_err(crate::types::VardictError::Io)?;
                let id = NEXT_BAM_META_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let mut reader = noodles::bam::io::Reader::from(CachedBgzfReader::new(file, id));
                let header = reader
                    .read_header()
                    .map_err(crate::types::VardictError::Io)?;
                let meta = std::sync::Arc::new(SharedBamMeta { id, header, index });
                *slot_guard = Some(std::sync::Arc::clone(&meta));
                return Ok((meta, reader));
            }
        }
    };

    if let Some(reader) = pop_pooled_reader(path, &meta) {
        return Ok((meta, reader));
    }
    let file = std::fs::File::open(path).map_err(crate::types::VardictError::Io)?;
    let reader = noodles::bam::io::Reader::from(CachedBgzfReader::new(file, meta.id));
    Ok((meta, reader))
}

/// Compute the 1-based inclusive alignment end directly from concrete
/// `noodles::bam::Record` inherent methods (`alignment_start` + `cigar().iter()`),
/// avoiding the two `Box<dyn ...>` allocations of `Record::alignment_end`.
fn bam_alignment_end(
    record: &noodles::bam::Record,
) -> Option<std::io::Result<noodles::core::Position>> {
    let start = match record.alignment_start().transpose() {
        Ok(position) => position?,
        Err(e) => return Some(Err(e)),
    };

    let mut span = 0usize;
    for result in record.cigar().iter() {
        let op = match result {
            Ok(op) => op,
            Err(e) => return Some(Err(e)),
        };
        if op.kind().consumes_reference() {
            span = span.saturating_add(op.len());
        }
    }

    if span == 0 {
        Some(Ok(start))
    } else {
        let end = usize::from(start).saturating_add(span).saturating_sub(1);
        noodles::core::Position::new(end).map(Ok)
    }
}

/// Owning lazy record feed over one queried region (Task 4: OOM fix;
/// Task 10: shared `.bai`/header metadata + pooled BGZF reader).
///
/// `noodles`' indexed `Query` borrows both the reader and the header, so it
/// cannot live in the same struct that owns them. `BamFeed` replicates the
/// indexed-query walk with owned state: the same `BinningIndex::query` chunk
/// list, the same csi Seek/Read chunk state machine, and the same `intersects`
/// accept/reject shape, visiting exactly the records `IndexedReader::query`
/// would, one at a time.
pub struct BamFeed {
    inner: BamFeedInner,
    sam: SamRecord,
}

struct BamFeedInner {
    path: PathBuf,
    meta: std::sync::Arc<SharedBamMeta>,
    reader: Option<RawBamReader>,
    tainted: bool,
    reference_sequence_id: usize,
    interval: noodles::core::region::Interval,
    chunks: std::vec::IntoIter<noodles::csi::binning_index::index::reference_sequence::bin::Chunk>,
    chunk_end: Option<noodles::bgzf::VirtualPosition>,
    record: noodles::bam::Record,
}

impl Drop for BamFeedInner {
    fn drop(&mut self) {
        if self.tainted {
            return;
        }
        let Some(reader) = self.reader.take() else {
            return;
        };
        let path = std::mem::take(&mut self.path);
        let meta = std::sync::Arc::clone(&self.meta);
        return_pooled_reader(path, meta, reader);
    }
}

/// Ports noodles-bam 0.96 `io::reader::query::intersects` (private there):
/// unmapped records and other-reference records are skipped; bounded
/// intervals require alignment overlap.
fn intersects_record(
    record: &noodles::bam::Record,
    reference_sequence_id: usize,
    region_interval: noodles::core::region::Interval,
) -> Result<bool, std::io::Error> {
    let Some(id) = record.reference_sequence_id().transpose()? else {
        return Ok(false);
    };

    if id != reference_sequence_id {
        return Ok(false);
    }

    if region_interval.start().is_none() && region_interval.end().is_none() {
        Ok(true)
    } else {
        match (
            record.alignment_start().transpose()?,
            bam_alignment_end(record).transpose()?,
        ) {
            (Some(start), Some(end)) => {
                let alignment_interval = (start..=end).into();
                Ok(region_interval.intersects(alignment_interval))
            }
            _ => Ok(false),
        }
    }
}

impl BamFeedInner {
    /// Advance one intersecting record, refilling `out` in place (zero
    /// per-read allocation after warmup). Returns `Ok(false)` at end of feed.
    fn pump(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
        use noodles::bgzf::io::Seek as _;
        let Some(reader) = self.reader.as_mut() else {
            return Ok(false);
        };
        loop {
            match self.chunk_end {
                None => {
                    let chunk = match self.chunks.next() {
                        None => return Ok(false),
                        Some(chunk) => chunk,
                    };
                    if let Err(e) = reader.get_mut().seek_to_virtual_position(chunk.start()) {
                        self.tainted = true;
                        return Err(crate::types::VardictError::Io(e));
                    }
                    self.chunk_end = Some(chunk.end());
                }
                Some(end) => {
                    if reader.get_mut().virtual_position() < end {
                        let n = match reader.read_record(&mut self.record) {
                            Ok(n) => n,
                            Err(e) => {
                                self.tainted = true;
                                return Err(crate::types::VardictError::Io(e));
                            }
                        };
                        if n == 0 {
                            self.chunk_end = None;
                            continue;
                        }
                        match intersects_record(
                            &self.record,
                            self.reference_sequence_id,
                            self.interval,
                        ) {
                            Ok(true) => {
                                if let Err(e) =
                                    out.fill_from_bam_record(&self.record, &self.meta.header)
                                {
                                    self.tainted = true;
                                    return Err(crate::types::VardictError::Io(e));
                                }
                                return Ok(true);
                            }
                            Ok(false) => continue,
                            Err(e) => {
                                self.tainted = true;
                                return Err(crate::types::VardictError::Io(e));
                            }
                        }
                    } else {
                        self.chunk_end = None;
                    }
                }
            }
        }
    }
}

impl BamFeed {
    /// Streaming hot path: refill `out` with the next record (`Ok(false)` at
    /// end). Identical `Err` semantics to collecting the `Iterator`.
    pub fn next_into(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
        self.inner.pump(out)
    }
}

/// `Iterator` adapter (required by the `RecordFeed: Iterator` interface).
/// Clones one `SamRecord` per item — production hot paths must use
/// [`BamFeed::next_into`] instead (zero-copy refill); this impl exists for
/// the interface shape and small-scale/test consumption only.
impl Iterator for BamFeed {
    type Item = Result<SamRecord, crate::types::VardictError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.pump(&mut self.sam) {
            Ok(false) => None,
            Ok(true) => Some(Ok(self.sam.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}

const FASTA_BUF_CAPACITY: usize = 256 * 1024;
const MAX_CACHED_FASTA_CONTIGS: usize = 8;

struct PreallocatedFastaAdapter {
    reader: noodles::fasta::io::Reader<std::io::BufReader<std::fs::File>>,
    index: noodles::fasta::fai::Index,
}

impl PreallocatedFastaAdapter {
    fn open(fasta_path: &Path) -> Result<Self, crate::types::VardictError> {
        let mut fai_path = fasta_path.to_path_buf();
        let mut ext = fai_path
            .extension()
            .map(|e| e.to_string_lossy().into_owned())
            .unwrap_or_default();
        ext.push_str(".fai");
        fai_path.set_extension(ext);
        let file = std::fs::File::open(fasta_path).map_err(crate::types::VardictError::Io)?;
        let index =
            noodles::fasta::fai::fs::read(&fai_path).map_err(crate::types::VardictError::Io)?;
        let reader = noodles::fasta::io::Reader::new(std::io::BufReader::with_capacity(
            FASTA_BUF_CAPACITY,
            file,
        ));
        Ok(Self { reader, index })
    }
}

impl noodles::fasta::repository::Adapter for PreallocatedFastaAdapter {
    fn get(&mut self, name: &[u8]) -> Option<std::io::Result<noodles::fasta::Record>> {
        use std::io::{BufRead as _, Seek as _, SeekFrom};
        let region = noodles::core::Region::new(name, ..);
        let pos = match self.index.query(&region) {
            Ok(p) => p,
            Err(e) => return Some(Err(e)),
        };
        if let Err(e) = self.reader.get_mut().seek(SeekFrom::Start(pos)) {
            return Some(Err(e));
        }
        let cap = self
            .index
            .as_ref()
            .iter()
            .find(|r| r.name() == name)
            .and_then(|r| usize::try_from(r.length()).ok())
            .unwrap_or(0);
        let mut raw_sequence = Vec::with_capacity(cap);
        let mut seq_reader = self.reader.sequence_reader();
        loop {
            let src = match seq_reader.fill_buf() {
                Ok([]) => break,
                Ok(s) => s,
                Err(e) => return Some(Err(e)),
            };
            raw_sequence.extend_from_slice(src);
            let consumed = src.len();
            seq_reader.consume(consumed);
        }
        let definition = noodles::fasta::record::Definition::new(region.to_string(), None);
        let sequence = noodles::fasta::record::Sequence::from(raw_sequence);
        Some(Ok(noodles::fasta::Record::new(definition, sequence)))
    }
}

struct SharedFastaAdapter(std::sync::Arc<std::sync::Mutex<PreallocatedFastaAdapter>>);

impl noodles::fasta::repository::Adapter for SharedFastaAdapter {
    fn get(&mut self, name: &[u8]) -> Option<std::io::Result<noodles::fasta::Record>> {
        self.0.lock().unwrap_or_else(|e| e.into_inner()).get(name)
    }
}

struct CachedContigRepo {
    name: Vec<u8>,
    repo: noodles::fasta::Repository,
    seq: Option<std::sync::Arc<noodles::fasta::record::Sequence>>,
}

struct SharedFastaRepo {
    adapter: std::sync::Arc<std::sync::Mutex<PreallocatedFastaAdapter>>,
    lru: std::sync::Mutex<Vec<CachedContigRepo>>,
}

impl SharedFastaRepo {
    fn warm_contig(
        &self,
        contig: &[u8],
    ) -> (
        noodles::fasta::Repository,
        Option<std::sync::Arc<noodles::fasta::record::Sequence>>,
    ) {
        let mut lru = self.lru.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(pos) = lru.iter().position(|e| e.name.as_slice() == contig) {
            let entry = lru.remove(pos);
            let repo = entry.repo.clone();
            let seq = entry.seq.clone();
            lru.push(entry);
            return (repo, seq);
        }
        if lru.len() >= MAX_CACHED_FASTA_CONTIGS {
            lru.remove(0);
        }
        let repo = noodles::fasta::Repository::new(SharedFastaAdapter(std::sync::Arc::clone(
            &self.adapter,
        )));
        let seq = if contig.is_empty() {
            None
        } else {
            repo.get(contig).and_then(|r| r.ok())
        };
        lru.push(CachedContigRepo {
            name: contig.to_vec(),
            repo: repo.clone(),
            seq: seq.clone(),
        });
        (repo, seq)
    }
}

static FASTA_REPO_CACHE: std::sync::OnceLock<
    std::sync::Mutex<HashMap<PathBuf, std::sync::Arc<SharedFastaRepo>>>,
> = std::sync::OnceLock::new();

fn acquire_fasta_repo(
    fasta_path: &Path,
) -> Result<std::sync::Arc<SharedFastaRepo>, crate::types::VardictError> {
    let cache = FASTA_REPO_CACHE.get_or_init(|| std::sync::Mutex::new(HashMap::new()));
    let mut map = cache.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(existing) = map.get(fasta_path) {
        return Ok(std::sync::Arc::clone(existing));
    }
    let adapter = PreallocatedFastaAdapter::open(fasta_path)?;
    let repo = std::sync::Arc::new(SharedFastaRepo {
        adapter: std::sync::Arc::new(std::sync::Mutex::new(adapter)),
        lru: std::sync::Mutex::new(Vec::with_capacity(MAX_CACHED_FASTA_CONTIGS)),
    });
    map.insert(fasta_path.to_path_buf(), std::sync::Arc::clone(&repo));
    Ok(repo)
}

struct SharedCramMeta {
    id: u64,
    header: noodles::sam::Header,
    index: noodles::cram::crai::Index,
    fasta_repo: std::sync::Arc<SharedFastaRepo>,
}

type CramMetaSlot = std::sync::Arc<std::sync::Mutex<Option<std::sync::Arc<SharedCramMeta>>>>;
type CramMetaKey = (PathBuf, PathBuf);
type RawCramReader = noodles::cram::io::Reader<std::fs::File>;
type PooledCramReaderEntry = (std::sync::Arc<SharedCramMeta>, RawCramReader);
type CramReaderPoolMap = HashMap<PathBuf, Vec<PooledCramReaderEntry>>;

static CRAM_META_CACHE: std::sync::OnceLock<std::sync::Mutex<HashMap<CramMetaKey, CramMetaSlot>>> =
    std::sync::OnceLock::new();

thread_local! {
    static CRAM_READER_POOL: RefCell<CramReaderPoolMap> = RefCell::new(HashMap::new());
}

fn pop_pooled_cram_reader(
    path: &Path,
    meta: &std::sync::Arc<SharedCramMeta>,
) -> Option<RawCramReader> {
    CRAM_READER_POOL
        .try_with(|pool| {
            let mut map = pool.try_borrow_mut().ok()?;
            let vec = map.get_mut(path)?;
            while let Some((pooled_meta, reader)) = vec.pop() {
                if std::sync::Arc::ptr_eq(&pooled_meta, meta) {
                    return Some(reader);
                }
            }
            None
        })
        .ok()
        .flatten()
}

fn return_pooled_cram_reader(
    path: PathBuf,
    meta: std::sync::Arc<SharedCramMeta>,
    reader: RawCramReader,
) {
    let _ = CRAM_READER_POOL.try_with(|pool| {
        if let Ok(mut map) = pool.try_borrow_mut() {
            let vec = map.entry(path).or_default();
            vec.retain(|(m, _)| std::sync::Arc::ptr_eq(m, &meta));
            if vec.len() < MAX_POOLED_READERS_PER_PATH {
                vec.push((meta, reader));
            }
        }
    });
}

fn acquire_cram_meta(
    path: &Path,
    config: &crate::config::Config,
) -> Result<std::sync::Arc<SharedCramMeta>, crate::types::VardictError> {
    let fasta_path = Path::new(&config.reference);
    let fasta_repo = acquire_fasta_repo(fasta_path)?;
    let key = (path.to_path_buf(), fasta_path.to_path_buf());
    let slot = {
        let cache = CRAM_META_CACHE.get_or_init(|| std::sync::Mutex::new(HashMap::new()));
        let mut map = cache.lock().unwrap_or_else(|e| e.into_inner());
        std::sync::Arc::clone(
            map.entry(key)
                .or_insert_with(|| std::sync::Arc::new(std::sync::Mutex::new(None))),
        )
    };

    let mut slot_guard = slot.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(meta) = slot_guard.as_ref() {
        if std::sync::Arc::ptr_eq(&meta.fasta_repo, &fasta_repo) {
            return Ok(std::sync::Arc::clone(meta));
        }
    }

    let index =
        noodles::cram::fs::read_associated_index(path).map_err(crate::types::VardictError::Io)?;
    let file = std::fs::File::open(path).map_err(crate::types::VardictError::Io)?;
    let mut reader = noodles::cram::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(crate::types::VardictError::Io)?;
    let id = NEXT_BAM_META_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let meta = std::sync::Arc::new(SharedCramMeta {
        id,
        header,
        index,
        fasta_repo,
    });
    *slot_guard = Some(std::sync::Arc::clone(&meta));
    return_pooled_cram_reader(path.to_path_buf(), std::sync::Arc::clone(&meta), reader);
    Ok(meta)
}

const CRAM_CONTAINER_CACHE_SHARDS: usize = 8;
const CRAM_CONTAINER_SLOTS_PER_SHARD: usize = 4;

#[cfg(test)]
static CRAM_CONTAINER_DECODE_COUNT: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);

#[derive(Clone, Debug)]
struct CachedCramRecord {
    reference_sequence_id: Option<usize>,
    alignment_start: Option<noodles::core::Position>,
    alignment_end: Option<noodles::core::Position>,
    sam: SamRecord,
}

type CachedCramBlock = std::sync::Arc<Vec<CachedCramRecord>>;
type CramCacheEntry = ((u64, u64), CachedCramBlock);

#[derive(Default)]
struct CramContainerCacheShard {
    entries: Vec<CramCacheEntry>,
}

impl CramContainerCacheShard {
    fn get(&mut self, key: (u64, u64)) -> Option<CachedCramBlock> {
        let idx = self.entries.iter().position(|(k, _)| *k == key)?;
        let entry = self.entries.remove(idx);
        let records = std::sync::Arc::clone(&entry.1);
        self.entries.push(entry);
        Some(records)
    }

    fn insert(&mut self, key: (u64, u64), records: CachedCramBlock) {
        if let Some(idx) = self.entries.iter().position(|(k, _)| *k == key) {
            self.entries.remove(idx);
        } else if self.entries.len() >= CRAM_CONTAINER_SLOTS_PER_SHARD {
            self.entries.remove(0);
        }
        self.entries.push((key, records));
    }

    fn clear_production(&mut self) {
        #[cfg(test)]
        self.entries.retain(|((id, _), _)| *id >= TEST_META_ID_BASE);
        #[cfg(not(test))]
        self.entries.clear();
    }
}

static CRAM_CONTAINER_CACHE: std::sync::OnceLock<
    [std::sync::Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS],
> = std::sync::OnceLock::new();

fn cram_container_cache_shards(
) -> &'static [std::sync::Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS] {
    CRAM_CONTAINER_CACHE.get_or_init(|| {
        std::array::from_fn(|_| std::sync::Mutex::new(CramContainerCacheShard::default()))
    })
}

#[inline]
fn cram_container_shard_index(meta_id: u64, offset: u64) -> usize {
    let h = offset
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(meta_id.wrapping_mul(0x517C_C1B7_2722_0A95));
    ((h ^ (h >> 32)) as usize) & (CRAM_CONTAINER_CACHE_SHARDS - 1)
}

fn crai_record_intersects(
    record: &noodles::cram::crai::Record,
    reference_sequence_id: usize,
    region_interval: noodles::core::region::Interval,
) -> bool {
    let Some(id) = record.reference_sequence_id() else {
        return false;
    };
    if id != reference_sequence_id {
        return false;
    }
    if region_interval.start().is_none() && region_interval.end().is_none() {
        true
    } else {
        let Some(start) = record.alignment_start() else {
            return false;
        };
        let span = record.alignment_span();
        if span == 0 {
            return false;
        }
        let Some(end) = start.checked_add(span - 1) else {
            return false;
        };
        let alignment_interval = (start..=end).into();
        region_interval.intersects(alignment_interval)
    }
}

fn cram_record_intersects(
    record: &CachedCramRecord,
    reference_sequence_id: usize,
    region_interval: noodles::core::region::Interval,
) -> bool {
    let Some(id) = record.reference_sequence_id else {
        return false;
    };
    if id != reference_sequence_id {
        return false;
    }
    if region_interval.start().is_none() && region_interval.end().is_none() {
        true
    } else {
        match (record.alignment_start, record.alignment_end) {
            (Some(start), Some(end)) => {
                let alignment_interval = (start..=end).into();
                region_interval.intersects(alignment_interval)
            }
            _ => false,
        }
    }
}

fn cram_record_to_cached(
    record: &noodles::cram::Record<'_>,
    header: &noodles::sam::Header,
) -> Result<CachedCramRecord, std::io::Error> {
    use noodles::sam::alignment::record::data::field::{Tag, Value};
    use noodles::sam::alignment::Record as _;

    let name = record
        .name()
        .map(|n| {
            let bytes: &[u8] = n.as_ref();
            bytes.to_vec()
        })
        .unwrap_or_default();
    let flags = record.flags()?.bits();
    let reference_sequence_id = record.reference_sequence_id(header).transpose()?;
    let reference_name = reference_sequence_id
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(n, _)| String::from_utf8_lossy(n).into_owned());
    let alignment_start = record.alignment_start().transpose()?;
    let raw_mapq = record.mapping_quality().transpose()?;
    let mapping_quality = if (flags & 0x4) != 0 {
        None
    } else {
        Some(raw_mapq.map(|m| m.get()).unwrap_or(255))
    };

    let raw_cigar = record.cigar();
    let mut cigar = Vec::with_capacity(raw_cigar.len());
    let mut alignment_span = 0usize;
    for result in raw_cigar.iter() {
        let op = result?;
        if op.kind().consumes_reference() {
            alignment_span += op.len();
        }
        cigar.push((op.len() as u32, op.kind()));
    }
    let alignment_end = match alignment_start {
        Some(start) => {
            if alignment_span == 0 {
                Some(start)
            } else {
                noodles::core::Position::new(usize::from(start) + alignment_span - 1)
            }
        }
        None => None,
    };

    let mate_reference_sequence_id = record.mate_reference_sequence_id(header).transpose()?;
    let mate_reference_name = mate_reference_sequence_id
        .and_then(|id| header.reference_sequences().get_index(id))
        .map(|(n, _)| String::from_utf8_lossy(n).into_owned());
    let mate_alignment_start = record.mate_alignment_start().transpose()?;
    let template_length = record.template_length()?;

    let raw_seq = record.sequence();
    let mut sequence = Vec::with_capacity(raw_seq.len());
    for b in raw_seq.iter() {
        sequence.push(b.to_ascii_uppercase());
    }

    let raw_quals = record.quality_scores();
    let mut quality_scores = Vec::with_capacity(raw_quals.len());
    for result in raw_quals.iter() {
        let score = result?;
        quality_scores.push(score.saturating_add(33));
    }

    const TAG_NM_LOWER: Tag = Tag::new(b'n', b'M');
    const TAG_NM_UPPER: Tag = Tag::EDIT_DISTANCE;
    const TAG_SA: Tag = Tag::new(b'S', b'A');
    const TAG_MC: Tag = Tag::new(b'M', b'C');
    const TAG_MQ: Tag = Tag::new(b'M', b'Q');

    let mut seen_nm_lower = false;
    let mut seen_nm_upper = false;
    let mut nm_lower_val: Option<i32> = None;
    let mut nm_upper_val: Option<i32> = None;
    let mut sa: Option<String> = None;
    let mut mc: Option<String> = None;
    let mut mq: Option<u8> = None;

    let raw_data = record.data();
    for result in raw_data.iter() {
        let (tag, val) = result?;
        if let Value::Array(ref arr) = val {
            use noodles::sam::alignment::record::data::field::value::Array;
            match arr {
                Array::Int8(v) => {
                    for r in v.iter() {
                        r?;
                    }
                }
                Array::UInt8(v) => {
                    for r in v.iter() {
                        r?;
                    }
                }
                Array::Int16(v) => {
                    for r in v.iter() {
                        r?;
                    }
                }
                Array::UInt16(v) => {
                    for r in v.iter() {
                        r?;
                    }
                }
                Array::Int32(v) => {
                    for r in v.iter() {
                        r?;
                    }
                }
                Array::UInt32(v) => {
                    for r in v.iter() {
                        r?;
                    }
                }
                Array::Float(v) => {
                    for r in v.iter() {
                        r?;
                    }
                }
            }
        }
        if tag == TAG_NM_LOWER {
            seen_nm_lower = true;
            nm_lower_val = val.as_int().map(|x| x as i32);
        } else if tag == TAG_NM_UPPER {
            seen_nm_upper = true;
            nm_upper_val = val.as_int().map(|x| x as i32);
        } else if tag == TAG_SA {
            let buf = sa.get_or_insert_with(String::new);
            buf.clear();
            match val {
                Value::String(s) => buf.push_str(&String::from_utf8_lossy(s.as_ref())),
                _ => {
                    let buf_val =
                        noodles::sam::alignment::record_buf::data::field::Value::try_from(val)?;
                    use std::fmt::Write as _;
                    let _ = write!(buf, "{buf_val:?}");
                }
            }
        } else if tag == TAG_MC {
            let buf = mc.get_or_insert_with(String::new);
            buf.clear();
            match val {
                Value::String(s) => buf.push_str(&String::from_utf8_lossy(s.as_ref())),
                _ => {
                    let buf_val =
                        noodles::sam::alignment::record_buf::data::field::Value::try_from(val)?;
                    use std::fmt::Write as _;
                    let _ = write!(buf, "{buf_val:?}");
                }
            }
        } else if tag == TAG_MQ {
            mq = val.as_int().map(|x| x as u8);
        }
    }

    let nm = if seen_nm_lower {
        nm_lower_val
    } else if seen_nm_upper {
        nm_upper_val
    } else {
        None
    };

    let sam = SamRecord {
        name,
        flags,
        reference_name,
        alignment_start: alignment_start.map(crate::types::Pos1::from),
        alignment_end: alignment_end.map(crate::types::Pos1::from),
        mapping_quality,
        cigar,
        sequence,
        quality_scores,
        mate_reference_name,
        mate_alignment_start: mate_alignment_start.map(crate::types::Pos1::from),
        template_length,
        nm,
        sa,
        mc,
        mq,
    };

    #[cfg(test)]
    {
        let rec_buf =
            noodles::sam::alignment::RecordBuf::try_from_alignment_record(header, record)?;
        let expected_sam = SamRecord::from_alignment_record(&rec_buf, header)?;
        assert_eq!(reference_sequence_id, rec_buf.reference_sequence_id());
        assert_eq!(alignment_start, rec_buf.alignment_start());
        assert_eq!(alignment_end, rec_buf.alignment_end());
        assert_eq!(sam, expected_sam);
    }

    Ok(CachedCramRecord {
        reference_sequence_id,
        alignment_start,
        alignment_end,
        sam,
    })
}

fn fetch_cram_container(
    path: &Path,
    meta: &std::sync::Arc<SharedCramMeta>,
    container_offset: u64,
    hint_contig: &[u8],
) -> Result<Option<CachedCramBlock>, crate::types::VardictError> {
    let shard_idx = cram_container_shard_index(meta.id, container_offset);
    if let Some(shard_mutex) = cram_container_cache_shards().get(shard_idx) {
        let mut shard = shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(records) = shard.get((meta.id, container_offset)) {
            return Ok(Some(records));
        }
    }

    let (hint_repo, hint_seq) = meta.fasta_repo.warm_contig(hint_contig);
    let mut reader = match pop_pooled_cram_reader(path, meta) {
        Some(r) => r,
        None => {
            let file = std::fs::File::open(path).map_err(crate::types::VardictError::Io)?;
            noodles::cram::io::Reader::new(file)
        }
    };

    reader
        .seek(std::io::SeekFrom::Start(container_offset))
        .map_err(crate::types::VardictError::Io)?;
    let mut container = noodles::cram::io::reader::Container::default();
    let bytes_read = reader
        .read_container(&mut container)
        .map_err(crate::types::VardictError::Io)?;
    return_pooled_cram_reader(path.to_path_buf(), std::sync::Arc::clone(meta), reader);
    if bytes_read == 0 {
        return Ok(None);
    }

    #[cfg(test)]
    CRAM_CONTAINER_DECODE_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

    let compression_header = container
        .compression_header()
        .map_err(crate::types::VardictError::Io)?;
    let mut out_records = Vec::new();

    for slice_res in container.slices() {
        let slice = slice_res.map_err(crate::types::VardictError::Io)?;
        let (core_data_src, external_data_srcs) = slice
            .decode_blocks()
            .map_err(crate::types::VardictError::Io)?;
        let cram_records = slice
            .records(
                hint_repo.clone(),
                &meta.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )
            .map_err(crate::types::VardictError::Io)?;

        out_records.reserve(cram_records.len());
        for record in cram_records {
            let mut cached = cram_record_to_cached(&record, &meta.header)
                .map_err(crate::types::VardictError::Io)?;
            if cached.sam.nm.is_none() {
                let rec_seq = match cached.sam.reference_name.as_deref() {
                    Some(rn) if rn.as_bytes() == hint_contig => hint_seq.clone(),
                    Some(rn) => meta.fasta_repo.warm_contig(rn.as_bytes()).1,
                    None => None,
                };
                if let Some(ref seq) = rec_seq {
                    cached.sam.nm = compute_cram_nm(&cached.sam, seq.as_ref().as_ref());
                }
            }
            out_records.push(cached);
        }
    }

    let arc_records: CachedCramBlock = std::sync::Arc::new(out_records);
    if let Some(shard_mutex) = cram_container_cache_shards().get(shard_idx) {
        let mut shard = shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
        shard.insert(
            (meta.id, container_offset),
            std::sync::Arc::clone(&arc_records),
        );
    }
    Ok(Some(arc_records))
}

/// CRAM arm of the feed: lazily streams one decoded container at a time from
/// [`CRAM_CONTAINER_CACHE`] / [`CRAM_READER_POOL`] (`AR-4`), refilling `out`
/// in place via [`SamRecord::clone_from_record`].
pub struct CramFeed {
    inner: CramFeedInner,
    sam: SamRecord,
}

struct CramFeedInner {
    path: PathBuf,
    meta: std::sync::Arc<SharedCramMeta>,
    reference_sequence_id: usize,
    hint_contig: Vec<u8>,
    interval: noodles::core::region::Interval,
    container_offsets: std::vec::IntoIter<u64>,
    current_block: Option<(CachedCramBlock, usize)>,
}

impl CramFeedInner {
    fn pump(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
        loop {
            if let Some((ref block, ref mut idx)) = self.current_block {
                while let Some(rec) = block.get(*idx) {
                    *idx += 1;
                    if cram_record_intersects(rec, self.reference_sequence_id, self.interval) {
                        out.clone_from_record(&rec.sam);
                        return Ok(true);
                    }
                }
                self.current_block = None;
            }
            let Some(offset) = self.container_offsets.next() else {
                return Ok(false);
            };
            match fetch_cram_container(&self.path, &self.meta, offset, &self.hint_contig)? {
                Some(block) => {
                    self.current_block = Some((block, 0));
                }
                None => {
                    self.container_offsets = Vec::new().into_iter();
                    return Ok(false);
                }
            }
        }
    }
}

impl CramFeed {
    fn empty(path: PathBuf, meta: std::sync::Arc<SharedCramMeta>) -> Self {
        Self {
            inner: CramFeedInner {
                path,
                meta,
                reference_sequence_id: 0,
                hint_contig: Vec::new(),
                interval: noodles::core::region::Interval::from(..),
                container_offsets: Vec::new().into_iter(),
                current_block: None,
            },
            sam: SamRecord::new(),
        }
    }

    pub fn next_into(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
        self.inner.pump(out)
    }
}

impl Iterator for CramFeed {
    type Item = Result<SamRecord, crate::types::VardictError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.pump(&mut self.sam) {
            Ok(false) => None,
            Ok(true) => Some(Ok(self.sam.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}

/// Owning streaming record feed returned by [`query_region_iter`].
pub enum RecordFeed {
    Bam(Box<BamFeed>),
    Cram(Box<CramFeed>),
}

impl Iterator for RecordFeed {
    type Item = Result<SamRecord, crate::types::VardictError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            RecordFeed::Bam(feed) => feed.next(),
            RecordFeed::Cram(feed) => feed.next(),
        }
    }
}

impl RecordFeed {
    /// Streaming hot path: refill `out` with the next record (`Ok(false)` at
    /// end). Identical `Err` semantics to collecting the `Iterator`: per-item
    /// I/O errors propagate via `Err`, never skipped.
    pub fn next_into(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
        match self {
            RecordFeed::Bam(feed) => feed.next_into(out),
            RecordFeed::Cram(feed) => feed.next_into(out),
        }
    }
}

/// Streaming counterpart to [`query_region`]: shares parsed `.bai`/header
/// metadata across threads via `META_CACHE`, reuses thread-local BGZF readers
/// via `BAM_READER_POOL`, and returns a lazy [`RecordFeed`] yielding one
/// `SamRecord` at a time. I/O errors surface per record and propagate with
/// `?` exactly like the Vec path's collection errors; `start > end` yields
/// an empty feed (preserving the `query_region` short-circuit after
/// metadata/reader acquisition).
pub fn query_region_iter(
    bam_path: &Path,
    chr: &str,
    start: Pos0,
    end: Pos0,
    config: &crate::config::Config,
) -> Result<RecordFeed, crate::types::VardictError> {
    let path_str = bam_path.to_string_lossy();
    if path_str.ends_with(".cram") || path_str.ends_with(".CRAM") {
        let meta = acquire_cram_meta(bam_path, config)?;
        let start_pos = start.to_one_based().0;
        let end_pos = end.0;
        if start_pos > end_pos {
            return Ok(RecordFeed::Cram(Box::new(CramFeed::empty(
                bam_path.to_path_buf(),
                meta,
            ))));
        }
        let region_str = format!("{}:{}-{}", chr, start_pos, end_pos);
        let region: noodles::core::Region = region_str.parse().map_err(|e| {
            crate::types::VardictError::Parse(format!(
                "failed to parse region '{}': {}",
                region_str, e
            ))
        })?;
        let reference_sequence_id = meta
            .header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                crate::types::VardictError::Io(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "invalid reference sequence name",
                ))
            })?;
        let interval = region.interval();
        let container_offsets: Vec<u64> = meta
            .index
            .iter()
            .filter(|rec| crai_record_intersects(rec, reference_sequence_id, interval))
            .map(|rec| rec.offset())
            .collect();
        Ok(RecordFeed::Cram(Box::new(CramFeed {
            inner: CramFeedInner {
                path: bam_path.to_path_buf(),
                meta,
                reference_sequence_id,
                hint_contig: region.name().to_vec(),
                interval,
                container_offsets: container_offsets.into_iter(),
                current_block: None,
            },
            sam: SamRecord::new(),
        })))
    } else {
        use noodles::csi::BinningIndex as _;

        let (meta, reader) = acquire_bam_meta_and_reader(bam_path)?;

        // Preserve the query_region empty-window short-circuit after
        // meta/reader acquisition so a missing BAM with an empty window
        // still errors (`lib.rs:94`).
        let start_pos = start.to_one_based().0;
        let end_pos = end.0;
        if start_pos > end_pos {
            return Ok(RecordFeed::Bam(Box::new(BamFeed {
                inner: BamFeedInner {
                    path: bam_path.to_path_buf(),
                    meta,
                    reader: Some(reader),
                    tainted: false,
                    reference_sequence_id: 0,
                    interval: noodles::core::region::Interval::from(..),
                    chunks: Vec::new().into_iter(),
                    chunk_end: None,
                    record: noodles::bam::Record::default(),
                },
                sam: SamRecord::new(),
            })));
        }

        let region_str = format!("{}:{}-{}", chr, start_pos, end_pos);
        let region: noodles::core::Region = match region_str.parse() {
            Ok(r) => r,
            Err(e) => {
                return_pooled_reader(bam_path.to_path_buf(), meta, reader);
                return Err(crate::types::VardictError::Parse(format!(
                    "failed to parse region '{}': {}",
                    region_str, e
                )));
            }
        };

        let Some(reference_sequence_id) = meta
            .header
            .reference_sequences()
            .get_index_of(region.name())
        else {
            return_pooled_reader(bam_path.to_path_buf(), meta, reader);
            return Err(crate::types::VardictError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "region reference sequence does not exist in reference sequences: {region:?}"
                ),
            )));
        };

        let chunks = match meta.index.query(reference_sequence_id, region.interval()) {
            Ok(c) => c,
            Err(e) => {
                return_pooled_reader(bam_path.to_path_buf(), meta, reader);
                return Err(crate::types::VardictError::Io(e));
            }
        };

        Ok(RecordFeed::Bam(Box::new(BamFeed {
            inner: BamFeedInner {
                path: bam_path.to_path_buf(),
                meta,
                reader: Some(reader),
                tainted: false,
                reference_sequence_id,
                interval: region.interval(),
                chunks: chunks.into_iter(),
                chunk_end: None,
                record: noodles::bam::Record::default(),
            },
            sam: SamRecord::new(),
        })))
    }
}

// Helper to reproduce getMateReferenceName(record) from VarDict
pub fn get_mate_reference_name(record: &SamRecord) -> String {
    match (&record.reference_name, &record.mate_reference_name) {
        (None, _) | (_, None) => "*".to_string(),
        (Some(ref_name), Some(mate_ref)) => {
            if ref_name == mate_ref {
                "=".to_string()
            } else {
                mate_ref.clone()
            }
        }
    }
}

/// Zero-alloc counterpart: true exactly when `get_mate_reference_name`
/// would return `"="` — the names are both present and equal, or the mate
/// name is the literal `"="` shorthand some constructors use for
/// same-chromosome. Hot paths that only test equality must use this instead
/// of allocating the String.
#[inline]
pub fn mate_reference_is_equal(record: &SamRecord) -> bool {
    match (&record.reference_name, &record.mate_reference_name) {
        (Some(_), Some(mate_ref)) if mate_ref == "=" => true,
        (Some(ref_name), Some(mate_ref)) => ref_name == mate_ref,
        _ => false,
    }
}

// Helper to format CIGAR to a string matching getCigarString()
pub fn get_cigar_string(cigar: &[(u32, CigarKind)]) -> String {
    let mut s = String::new();
    push_cigar_string(&mut s, cigar);
    s
}

/// Append CIGAR to an existing buffer (scratch reuse; identical bytes to
/// [`get_cigar_string`]).
pub fn push_cigar_string(buf: &mut String, cigar: &[(u32, CigarKind)]) {
    use std::fmt::Write as _;
    for (len, kind) in cigar {
        // write! formats in place (no temp String per op); bytes identical.
        let _ = write!(buf, "{len}");
        let c = match kind {
            CigarKind::Match => 'M',
            CigarKind::Insertion => 'I',
            CigarKind::Deletion => 'D',
            CigarKind::Skip => 'N',
            CigarKind::SoftClip => 'S',
            CigarKind::HardClip => 'H',
            CigarKind::Pad => 'P',
            CigarKind::SequenceMatch => '=',
            CigarKind::SequenceMismatch => 'X',
        };
        buf.push(c);
    }
}

/// Compute the aligned length on the reference using only Match and Deletion CIGAR operators.
pub fn get_aligned_length_md(cigar: &[(u32, CigarKind)]) -> i32 {
    cigar
        .iter()
        .filter(|&&(_, kind)| kind == CigarKind::Match || kind == CigarKind::Deletion)
        .map(|&(len, _)| len as i32)
        .sum()
}

/// Class for primary processing of SAM record and filtering reads that don't fit the requirements.
/// Ports RecordPreprocessor.java.
pub struct RecordPreprocessor {
    pub duplicates: HashSet<String>,
    pub first_matching_position: i32,
    pub total_reads: u32,
    pub duplicate_reads: u32,
}

impl RecordPreprocessor {
    pub fn new() -> Self {
        RecordPreprocessor {
            duplicates: HashSet::new(),
            first_matching_position: -1,
            total_reads: 0,
            duplicate_reads: 0,
        }
    }

    /// Preprocesses a record, returning true if it passes all filters.
    pub fn preprocess_record(
        &mut self,
        record: &SamRecord,
        config: &crate::config::Config,
    ) -> bool {
        // 1. Mapping quality filter
        if let Some(min_mapq) = config.map_qual {
            if let Some(mapq) = record.mapping_quality {
                if (mapq as u32) < min_mapq {
                    return false;
                }
            } else {
                // Missing MAPQ is treated as 255 which passes
                if 255 < min_mapq {
                    return false;
                }
            }
        }

        // 2. Secondary alignment filter (unless samfilter is explicitly 0)
        let is_secondary = (record.flags & 0x100) != 0;
        if is_secondary && config.flag_filter != 0 {
            return false;
        }

        // Unmapped reads filter
        if (record.flags & 0x4) != 0 {
            return false;
        }

        // 3. Sequence is not stored check
        if record.sequence.is_empty() || (record.sequence.len() == 1 && record.sequence[0] == b'*')
        {
            return false;
        }

        // 4. SAM flag filter
        if config.flag_filter != 0 && (u32::from(record.flags) & config.flag_filter) != 0 {
            return false;
        }

        // CIGAR query consumption consistency check (SAMv1 §1.4.6)
        // Corrupted or hostile reads where CIGAR query consumption exceeds SEQ length
        // or quality scores length cause out-of-bounds indexing in downstream pileup/parsers.
        let query_consumption =
            crate::cigar::parse::get_soft_clipped_length(&record.cigar) as usize;
        if query_consumption > record.sequence.len()
            || (!record.quality_scores.is_empty()
                && query_consumption > record.quality_scores.len())
        {
            log::warn!(
                "Dropping corrupt BAM record {}: CIGAR query consumption ({}) exceeds SEQ length ({})",
                String::from_utf8_lossy(&record.name),
                query_consumption,
                record.sequence.len()
            );
            return false;
        }

        self.total_reads += 1;

        // 5. Duplicated reads filter (removeDuplicatedReads)
        if config.dedup {
            let alignment_start_val = record.alignment_start.map(|p| p.0 as i32).unwrap_or(0);
            let mate_alignment_start_val =
                record.mate_alignment_start.map(|p| p.0 as i32).unwrap_or(0);

            if alignment_start_val != self.first_matching_position {
                self.duplicates.clear();
            }

            if mate_alignment_start_val < 10 {
                // POS-RNEXT-PNEXT
                let dup_key = format!(
                    "{}-{}-{}",
                    alignment_start_val,
                    get_mate_reference_name(record),
                    mate_alignment_start_val
                );
                if self.duplicates.contains(&dup_key) {
                    self.duplicate_reads += 1;
                    return false;
                }
                self.duplicates.insert(dup_key);
                self.first_matching_position = alignment_start_val;
            } else {
                let is_paired = (record.flags & 0x1) != 0;
                let is_mate_unmapped = (record.flags & 0x8) != 0;
                if is_paired && is_mate_unmapped {
                    // POS-CIGAR
                    let dup_key = format!(
                        "{}-{}",
                        alignment_start_val,
                        get_cigar_string(&record.cigar)
                    );
                    if self.duplicates.contains(&dup_key) {
                        self.duplicate_reads += 1;
                        return false;
                    }
                    self.duplicates.insert(dup_key);
                    self.first_matching_position = alignment_start_val;
                }
            }
        }

        true
    }
}

impl Default for RecordPreprocessor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use std::fs::File;
    use std::io::Write;

    /// Task 4 (red): the streaming feed must reproduce the Vec path exactly —
    /// same record ORDER (names), same aggregates, same read/dup counts, same
    /// Err on corrupt input, same empty-window behavior. Memory is proven by
    /// the Step 5 RSS run, not here.
    #[test]
    fn streaming_feed_matches_vec_feed() {
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record_buf::Cigar;
        use noodles::sam::alignment::record_buf::RecordBuf;

        fn record_buf(
            name: &str,
            start: usize,
            cigar_ops: &[(usize, Kind)],
            seq: Vec<u8>,
        ) -> RecordBuf {
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

        fn write_deep_bam(path: &Path, recs: &[RecordBuf]) {
            use noodles::bam::bai;
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
                .set_length(1000usize.try_into().unwrap())
                .build()
                .unwrap();
            let hd = Header::builder().add_reference_sequence("chr1", rs).build();
            let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
            let mut w = noodles::bam::io::Writer::from(bgzf);
            w.write_alignment_header(&hd).unwrap();
            let mut indexer = Indexer::<LinearIndex>::new(14, 5).expect("valid index parameters");
            for rec in recs {
                let start = w.get_ref().virtual_position();
                w.write_alignment_record(&hd, rec).unwrap();
                let end = w.get_ref().virtual_position();
                let s = rec.alignment_start().expect("aligned record");
                let e =
                    Position::new(s.get() + rec.cigar().alignment_span()).unwrap_or(Position::MIN);
                indexer
                    .add_record(Some((0, s, e, true)), Chunk::new(start, end))
                    .map_err(std::io::Error::other)
                    .unwrap();
            }
            w.try_finish().unwrap();
            let mut bai_path = PathBuf::from(path);
            bai_path.set_extension("bam.bai");
            bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }

        // Deep-overlap shape: ~8100 reads at starts 100..400 (position
        // sorted — the BAM indexer requires it), every 40th spliced; one
        // corrupt-CIGAR record (preprocess drops it in both paths). Names
        // encode file order, so any feed reordering breaks the assertion.
        let mut recs = Vec::new();
        let mut file_idx = 0;
        for start in 100..400usize {
            for k in 0..27 {
                let name = format!("f{file_idx:05}");
                file_idx += 1;
                if (start + k) % 40 == 0 {
                    recs.push(record_buf(
                        &name,
                        start,
                        &[(10, Kind::Match), (50, Kind::Skip), (10, Kind::Match)],
                        vec![b'A'; 20],
                    ));
                } else {
                    recs.push(record_buf(
                        &name,
                        start,
                        &[(50, Kind::Match)],
                        vec![b'A'; 50],
                    ));
                }
            }
        }
        // NOTE: no corrupt-CIGAR record here — noodles' writer validates
        // CIGAR-vs-SEQ consistency, so file-level corruption is exercised via
        // the truncated-file Err case below (and drop parity in-memory).
        // The BAM indexer requires non-decreasing starts: stable-sort by
        // start (names still encode file order within each start block).
        recs.sort_by_key(|r| r.alignment_start().unwrap().get());

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("deep.bam");
        write_deep_bam(&bam_path, &recs);

        let config = Config {
            dedup: true,
            flag_filter: 0,
            map_qual: None,
            ..Default::default()
        };
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
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(1000),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        // Vec path: query + preprocess + pileup.
        let vec_recs = query_region(&bam_path, "chr1", Pos0(50), Pos0(500), &config).unwrap();
        let vec_order: Vec<Vec<u8>> = vec_recs.iter().map(|r| r.name.clone()).collect();
        let mut pre = RecordPreprocessor::new();
        let mut kept = Vec::new();
        for r in &vec_recs {
            if pre.preprocess_record(r, &config) {
                kept.push(r.clone());
            }
        }
        let (vec_total, vec_dup) = (pre.total_reads, pre.duplicate_reads);
        let vec_data = crate::cigar::parse::pileup(
            &kept,
            &reference,
            &region,
            &config,
            false,
            vec_total as i32,
            vec_dup as i32,
            0,
        )
        .unwrap();

        // Stream path: feed + preprocess inline + per-record parse.
        let feed = query_region_iter(&bam_path, "chr1", Pos0(50), Pos0(500), &config).unwrap();
        let mut stream_order = Vec::new();
        let mut pre2 = RecordPreprocessor::new();
        let mut parser =
            crate::cigar::parse::CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
        // F6: pre-poison every parser scratch buffer with garbage; take +
        // clear + refill on every read must erase it (asserted field-for-field
        // below via the aggregate equality).
        parser.scratch_seq = vec![b'N'; 500];
        parser.scratch_qual = vec![0; 500];
        parser.scratch_cigar_str = "GARBAGE".to_string();
        parser.scratch_elements = vec![(999, CigarKind::Match)];
        parser.scratch_s = "GARBAGE".to_string();
        parser.scratch_ss = "GARBAGE".to_string();
        parser.scratch_splice_key = "GARBAGE".to_string();
        parser.scratch_f4_desc = "GARBAGE".to_string();
        parser.scratch_f4_qual = "GARBAGE".to_string();
        parser.scratch_f4_seg1 = "GARBAGE".to_string();
        parser.scratch_f4_seg2 = "GARBAGE".to_string();
        for item in feed {
            let r = item.unwrap();
            stream_order.push(r.name.clone());
            if pre2.preprocess_record(&r, &config) {
                parser.parse_cigar(&r).unwrap();
            }
        }
        parser.total_reads = pre2.total_reads as i32;
        parser.duplicate_reads = pre2.duplicate_reads as i32;
        let stream_data = parser.into_variation_data();

        // (a) identical order in file order (names encode file sequence).
        assert_eq!(stream_order, vec_order);
        assert!(vec_order.len() > 8000, "deep shape preserved");
        // File order is start-major: every name's index increases with the
        // block, so a feed that sorted by anything else would differ.
        assert_eq!(vec_order[0], b"f00000");
        assert_eq!(
            vec_order[vec_order.len() - 1],
            b"f08099",
            "last written record is last in both paths"
        );
        // Out-of-order guard (in-memory, no indexer): the streaming consumer
        // pattern (preprocess + per-record parse) must match the Vec consumer
        // (preprocess + pileup) even when records arrive unsorted.
        {
            fn mem_record(name: &str, start: u32) -> SamRecord {
                SamRecord {
                    name: name.as_bytes().to_vec(),
                    flags: 0,
                    reference_name: Some("chr1".to_string()),
                    alignment_start: Some(crate::types::Pos1(start)),
                    alignment_end: Some(crate::types::Pos1(start + 49)),
                    mapping_quality: Some(60),
                    cigar: vec![(50, CigarKind::Match)],
                    sequence: vec![b'A'; 50],
                    quality_scores: vec![40; 50],
                    mate_reference_name: None,
                    mate_alignment_start: None,
                    template_length: 0,
                    nm: None,
                    sa: None,
                    mc: None,
                    mq: None,
                }
            }
            // Deliberately unsorted starts with repeats (dup accounting) and
            // one corrupt-CIGAR record (dropped in both paths).
            let starts: [u32; 30] = [
                300, 100, 200, 100, 150, 400, 100, 250, 300, 120, 180, 400, 110, 350, 200, 130,
                170, 300, 140, 160, 190, 210, 220, 230, 240, 260, 270, 280, 290, 310,
            ];
            let mut unordered: Vec<SamRecord> = starts
                .iter()
                .enumerate()
                .map(|(i, s)| mem_record(&format!("u{i:02}"), *s))
                .collect();
            unordered.push(SamRecord {
                name: b"badcigar".to_vec(),
                cigar: vec![(104, CigarKind::Match)],
                sequence: vec![b'A'; 60],
                ..mem_record("badcigar", 155)
            });
            let mut pa = RecordPreprocessor::new();
            let mut kept_a = Vec::new();
            let mut order_a = Vec::new();
            for r in &unordered {
                order_a.push(r.name.clone());
                if pa.preprocess_record(r, &config) {
                    kept_a.push(r.clone());
                }
            }
            let data_a = crate::cigar::parse::pileup(
                &kept_a,
                &reference,
                &region,
                &config,
                false,
                pa.total_reads as i32,
                pa.duplicate_reads as i32,
                0,
            )
            .unwrap();
            let mut pb = RecordPreprocessor::new();
            let mut order_b = Vec::new();
            let mut parser_b =
                crate::cigar::parse::CigarParser::new(&reference, &region, &config, false, 0, 0, 0);
            for r in unordered.iter() {
                order_b.push(r.name.clone());
                if pb.preprocess_record(r, &config) {
                    parser_b.parse_cigar(r).unwrap();
                }
            }
            parser_b.total_reads = pb.total_reads as i32;
            parser_b.duplicate_reads = pb.duplicate_reads as i32;
            let data_b = parser_b.into_variation_data();
            assert_eq!(order_b, order_a);
            assert_eq!(pb.total_reads, pa.total_reads);
            assert_eq!(pb.duplicate_reads, pa.duplicate_reads);
            assert_eq!(data_b.splice_count, data_a.splice_count);
            assert_eq!(data_b.ref_coverage, data_a.ref_coverage);
        }
        // (c) identical read/dup accounting.
        assert_eq!(pre2.total_reads, vec_total);
        assert_eq!(pre2.duplicate_reads, vec_dup);
        // (b) identical aggregates.
        assert_eq!(stream_data.splice_count, vec_data.splice_count);
        assert_eq!(stream_data.ref_coverage, vec_data.ref_coverage);
        assert_eq!(stream_data.duprate, vec_data.duprate);
        assert_eq!(stream_data.max_read_length, vec_data.max_read_length);
        assert_eq!(stream_data.splice, vec_data.splice);
        assert_eq!(stream_data.mnp, vec_data.mnp);
        assert_eq!(
            stream_data.position_to_insertion_count,
            vec_data.position_to_insertion_count
        );
        assert_eq!(
            stream_data.position_to_deletion_count,
            vec_data.position_to_deletion_count
        );
        assert_eq!(
            stream_data.non_insertion_variants.len(),
            vec_data.non_insertion_variants.len()
        );

        // Empty window: both paths yield zero records without error.
        let empty_vec = query_region(&bam_path, "chr1", Pos0(500), Pos0(100), &config).unwrap();
        assert!(empty_vec.is_empty());
        let empty_stream: Vec<_> =
            query_region_iter(&bam_path, "chr1", Pos0(500), Pos0(100), &config)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
        assert!(empty_stream.is_empty());

        // (e) reuse-buffer path (`next_into` + `fill_from_alignment_record`)
        // must yield records identical to fresh conversion, field for field
        // (SamRecord: PartialEq) — guards against stale-field leaks. Every
        // field is pre-poisoned (long garbage buffers, stale tags, non-zero
        // scalars), so the first refill must clear all of them.
        {
            let mut feed =
                query_region_iter(&bam_path, "chr1", Pos0(50), Pos0(500), &config).unwrap();
            let mut reused = SamRecord {
                name: b"STALENAME".to_vec(),
                flags: 0xffff,
                reference_name: Some("STALE".to_string()),
                alignment_start: Some(crate::types::Pos1(999999)),
                alignment_end: Some(crate::types::Pos1(999999)),
                mapping_quality: Some(99),
                cigar: vec![(999, CigarKind::Match)],
                sequence: vec![b'N'; 300],
                quality_scores: vec![0; 300],
                mate_reference_name: Some("STALE".to_string()),
                mate_alignment_start: Some(crate::types::Pos1(999999)),
                template_length: 99999,
                nm: Some(99),
                sa: Some("STALE".to_string()),
                mc: Some("STALE".to_string()),
                mq: Some(99),
            };
            let mut refilled = Vec::new();
            while feed.next_into(&mut reused).unwrap() {
                refilled.push(reused.clone());
            }
            assert_eq!(refilled, vec_recs);
        }

        // (d) corrupt file: identical Err kind in both paths — mid-iteration
        // (multi-block deep BAM truncated mid-stream: header reads fine, a
        // later record read fails) and at construction (single-block file
        // truncated into its only block: the header read itself fails).
        fn err_kind(e: &crate::types::VardictError) -> String {
            match e {
                crate::types::VardictError::Io(io) => format!("Io({:?})", io.kind()),
                other => format!("{other:?}"),
            }
        }
        let trunc_path = dir.path().join("trunc.bam");
        std::fs::copy(&bam_path, &trunc_path).unwrap();
        std::fs::copy(
            bam_path.with_extension("bam.bai"),
            trunc_path.with_extension("bam.bai"),
        )
        .unwrap();
        {
            let len = std::fs::metadata(&trunc_path).unwrap().len();
            // ~8100 records x ~140 B uncompressed each >> one 64 KB BGZF
            // block, so the file spans many blocks: truncating at 3/4 keeps
            // the header block intact while a later record read fails.
            let f = std::fs::OpenOptions::new()
                .write(true)
                .open(&trunc_path)
                .unwrap();
            f.set_len(len * 3 / 4).unwrap();
        }
        let vec_err = query_region(&trunc_path, "chr1", Pos0(50), Pos0(500), &config)
            .expect_err("truncated BAM must fail the Vec path mid-iteration");
        let stream_err = match query_region_iter(&trunc_path, "chr1", Pos0(50), Pos0(500), &config)
        {
            Ok(feed) => feed
                .collect::<Result<Vec<_>, _>>()
                .expect_err("truncated BAM must fail the stream path mid-iteration"),
            Err(e) => panic!("stream construction must survive mid-file truncation: {e:?}"),
        };
        assert_eq!(err_kind(&stream_err), err_kind(&vec_err));
        let corrupt_path = dir.path().join("tiny.bam");
        {
            let small: Vec<RecordBuf> = (0..10)
                .map(|i| {
                    record_buf(
                        &format!("s{i:02}"),
                        100 + i,
                        &[(50, Kind::Match)],
                        vec![b'A'; 50],
                    )
                })
                .collect();
            write_deep_bam(&corrupt_path, &small);
            // Truncate mid-record so the record read fails.
            let len = std::fs::metadata(&corrupt_path).unwrap().len();
            let f = std::fs::OpenOptions::new()
                .write(true)
                .open(&corrupt_path)
                .unwrap();
            f.set_len(len / 2).unwrap();
        }
        let vec_err = query_region(&corrupt_path, "chr1", Pos0(50), Pos0(500), &config)
            .expect_err("truncated BAM must fail the Vec path");
        let stream_err =
            match query_region_iter(&corrupt_path, "chr1", Pos0(50), Pos0(500), &config) {
                Ok(feed) => feed
                    .collect::<Result<Vec<_>, _>>()
                    .expect_err("truncated BAM must fail the stream path"),
                Err(e) => e,
            };
        assert_eq!(err_kind(&stream_err), err_kind(&vec_err));
    }

    #[test]
    fn test_bam_open_and_query() {
        let bam_path = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam"
        ));
        let config = Config::default();

        let cached = open_reader(bam_path, &config).unwrap();
        let ref_seqs = cached.header.reference_sequences();
        assert!(
            !ref_seqs.is_empty(),
            "Reference sequences should not be empty"
        );

        let records =
            query_region(bam_path, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();

        assert!(!records.is_empty(), "Should load records from BAM file");
    }

    #[test]
    fn flag_filter() {
        let bam_path = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam"
        ));

        // 1. Without flag filter
        let config_no_filter = Config {
            flag_filter: 0,
            ..Default::default()
        };

        let records = query_region(
            bam_path,
            "chr7",
            Pos0(55259400),
            Pos0(55259600),
            &config_no_filter,
        )
        .unwrap();

        assert_eq!(
            records.len(),
            1327,
            "Total unfiltered record count must match samtools"
        );

        // 2. With flag filter 0x504
        let config_with_filter = Config {
            flag_filter: 0x504,
            ..Default::default()
        };

        let mut preprocessor = RecordPreprocessor::new();
        let mut passed_count = 0;

        for record in &records {
            if preprocessor.preprocess_record(record, &config_with_filter) {
                passed_count += 1;
            }
        }

        assert_eq!(
            passed_count, 1031,
            "Filtered count must match samtools view -c -F 0x504"
        );
    }

    #[test]
    fn test_cram_twin_parity() {
        // Detect samtools availability
        let mut has_samtools = std::process::Command::new("samtools")
            .arg("--version")
            .output()
            .is_ok();
        let mut use_conda = false;
        if !has_samtools {
            let conda_check = std::process::Command::new("conda")
                .args(["run", "-n", "vardict-oracle", "samtools", "--version"])
                .output();
            if conda_check.is_ok() {
                has_samtools = true;
                use_conda = true;
            }
        }
        if !has_samtools {
            println!("Skipping CRAM parity test: samtools not found");
            return;
        }

        let temp_dir = std::env::temp_dir().join("cram_test");
        let _ = std::fs::create_dir_all(&temp_dir);

        let ref_fa = temp_dir.join("ref.fa");
        let ref_fai = temp_dir.join("ref.fa.fai");
        let slice_sam = temp_dir.join("slice.sam");
        let new_sam_path = temp_dir.join("new_slice.sam");
        let slice_bam = temp_dir.join("slice.bam");
        let cram_path = temp_dir.join("test.cram");
        let bam_path = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam"
        ));

        // 1. Slice BAM to SAM
        let mut slice_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "view",
                "-h",
                "-o",
                slice_sam.to_str().unwrap(),
                bam_path.to_str().unwrap(),
                "chr7:55259401-55259600",
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args([
                "view",
                "-h",
                "-o",
                slice_sam.to_str().unwrap(),
                bam_path.to_str().unwrap(),
                "chr7:55259401-55259600",
            ]);
            c
        };
        assert!(slice_cmd.output().unwrap().status.success());

        // 2. Read slice.sam and write a new SAM file with a header containing only chr7
        let sam_content = std::fs::read_to_string(&slice_sam).unwrap();
        let mut new_sam = File::create(&new_sam_path).unwrap();
        let mut wrote_sq = false;
        for line in sam_content.lines() {
            if line.starts_with("@SQ") {
                if line.contains("SN:chr7") && !wrote_sq {
                    writeln!(new_sam, "{}", line).unwrap();
                    wrote_sq = true;
                }
            } else {
                writeln!(new_sam, "{}", line).unwrap();
            }
        }
        drop(new_sam);

        // 3. Convert SAM back to BAM (this re-encodes the reference IDs correctly)
        let mut convert_bam_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "view",
                "-b",
                "-o",
                slice_bam.to_str().unwrap(),
                new_sam_path.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args([
                "view",
                "-b",
                "-o",
                slice_bam.to_str().unwrap(),
                new_sam_path.to_str().unwrap(),
            ]);
            c
        };
        assert!(convert_bam_cmd.output().unwrap().status.success());

        // 4. Index the sliced BAM
        let mut slice_index = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "index",
                slice_bam.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args(["index", slice_bam.to_str().unwrap()]);
            c
        };
        assert!(slice_index.output().unwrap().status.success());

        // 5. Write a temporary mock reference FASTA with 159,138,663 N bases for chr7 to match BAM header chromosome length
        let mut fa_file = File::create(&ref_fa).unwrap();
        writeln!(fa_file, ">chr7").unwrap();
        let chunk = vec![b'N'; 100_000];
        for _ in 0..1591 {
            fa_file.write_all(&chunk).unwrap();
            fa_file.write_all(b"\n").unwrap();
        }
        let remainder = vec![b'N'; 38663];
        fa_file.write_all(&remainder).unwrap();
        fa_file.write_all(b"\n").unwrap();
        drop(fa_file);

        // 6. Run samtools faidx
        let mut faidx_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "faidx",
                ref_fa.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args(["faidx", ref_fa.to_str().unwrap()]);
            c
        };
        assert!(faidx_cmd.output().unwrap().status.success());

        // 7. Run samtools to convert slice_bam to CRAM
        let mut convert_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "view",
                "-C",
                "-T",
                ref_fa.to_str().unwrap(),
                "-o",
                cram_path.to_str().unwrap(),
                slice_bam.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args([
                "view",
                "-C",
                "-T",
                ref_fa.to_str().unwrap(),
                "-o",
                cram_path.to_str().unwrap(),
                slice_bam.to_str().unwrap(),
            ]);
            c
        };
        assert!(convert_cmd.output().unwrap().status.success());

        // 8. Run samtools index on CRAM file
        let mut index_cmd = if use_conda {
            let mut c = std::process::Command::new("conda");
            c.args([
                "run",
                "-n",
                "vardict-oracle",
                "samtools",
                "index",
                cram_path.to_str().unwrap(),
            ]);
            c
        } else {
            let mut c = std::process::Command::new("samtools");
            c.args(["index", cram_path.to_str().unwrap()]);
            c
        };
        assert!(index_cmd.output().unwrap().status.success());

        // Now query both BAM and CRAM and verify parity
        let config = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };

        // Query BAM
        let bam_records =
            query_region(&slice_bam, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();

        // Query CRAM
        let cram_records =
            query_region(&cram_path, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();

        assert_eq!(
            bam_records.len(),
            cram_records.len(),
            "BAM and CRAM must have the exact same number of reads"
        );

        for (i, (bam_rec, cram_rec)) in bam_records.iter().zip(cram_records.iter()).enumerate() {
            assert_eq!(bam_rec, cram_rec, "Parity mismatch at record {}", i);
        }

        // Verify streaming CramFeed (`query_region_iter`) matches `query_region` (`IndexedReader::query`)
        // and serves repeated sub-window queries from `CRAM_CONTAINER_CACHE` with 0 re-decodes.
        NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, std::sync::atomic::Ordering::Relaxed);
        let decodes_before = CRAM_CONTAINER_DECODE_COUNT.load(std::sync::atomic::Ordering::Relaxed);
        let mut feed1 =
            query_region_iter(&cram_path, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();
        let mut iter_records_1 = Vec::new();
        let mut reusable_rec = SamRecord::new();
        while feed1.next_into(&mut reusable_rec).unwrap() {
            iter_records_1.push(reusable_rec.clone());
        }
        assert_eq!(iter_records_1, cram_records);
        let decodes_after_1 =
            CRAM_CONTAINER_DECODE_COUNT.load(std::sync::atomic::Ordering::Relaxed);
        assert!(
            decodes_after_1 > decodes_before,
            "first CRAM query must decode container"
        );

        let iter_records_2: Vec<SamRecord> =
            query_region_iter(&cram_path, "chr7", Pos0(55259450), Pos0(55259550), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        let expected_sub: Vec<SamRecord> =
            query_region(&cram_path, "chr7", Pos0(55259450), Pos0(55259550), &config).unwrap();
        assert_eq!(iter_records_2, expected_sub);
        let decodes_after_2 =
            CRAM_CONTAINER_DECODE_COUNT.load(std::sync::atomic::Ordering::Relaxed);
        assert_eq!(
            decodes_after_2, decodes_after_1,
            "second CRAM sub-window query must hit CRAM_CONTAINER_CACHE without re-decoding"
        );

        // Verify empty-window short-circuit (`start > end`) and missing-contig error on CRAM
        let mut empty_cram =
            query_region_iter(&cram_path, "chr7", Pos0(55259600), Pos0(55259400), &config).unwrap();
        assert!(!empty_cram.next_into(&mut reusable_rec).unwrap());
        assert!(
            query_region_iter(&cram_path, "chr_missing", Pos0(10), Pos0(100), &config).is_err()
        );

        // Cleanup
        let _ = std::fs::remove_file(&ref_fa);
        let _ = std::fs::remove_file(&ref_fai);
        let _ = std::fs::remove_file(&slice_sam);
        let _ = std::fs::remove_file(&new_sam_path);
        let _ = std::fs::remove_file(&slice_bam);
        let _ = std::fs::remove_file(slice_bam.with_extension("bam.bai"));
        let _ = std::fs::remove_file(&cram_path);
        let _ = std::fs::remove_file(cram_path.with_extension("cram.crai"));
    }

    #[test]
    fn test_corrupt_cigar_query_consumption_exceeds_seq_len_is_dropped() {
        let mut preprocessor = RecordPreprocessor::new();
        let config = Config::default();
        let bad_record = SamRecord {
            name: b"corrupt1".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(50)),
            alignment_end: Some(Pos1(150)),
            mapping_quality: Some(60),
            cigar: vec![
                (104, CigarKind::Match),
                (1, CigarKind::Insertion),
                (2, CigarKind::Match),
            ],
            sequence: vec![b'A'; 60], // len = 60 < 107
            quality_scores: vec![40; 60],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: None,
            sa: None,
            mc: None,
            mq: None,
        };
        assert!(!preprocessor.preprocess_record(&bad_record, &config));
    }

    #[test]
    fn test_corrupt_softclip_exceeds_seq_len_is_dropped() {
        let mut preprocessor = RecordPreprocessor::new();
        let config = Config::default();
        let bad_record = SamRecord {
            name: b"corrupt2".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(50)),
            alignment_end: Some(Pos1(54)),
            mapping_quality: Some(60),
            cigar: vec![(100, CigarKind::SoftClip), (4, CigarKind::Match)],
            sequence: b"ACGT".to_vec(), // len = 4 < 104
            quality_scores: vec![40; 4],
            mate_reference_name: None,
            mate_alignment_start: None,
            template_length: 0,
            nm: None,
            sa: None,
            mc: None,
            mq: None,
        };
        assert!(!preprocessor.preprocess_record(&bad_record, &config));
    }

    impl SamRecord {
        fn fill_from_alignment_record_legacy<R>(
            &mut self,
            record: &R,
            header: &noodles::sam::Header,
        ) -> Result<(), std::io::Error>
        where
            R: Record,
        {
            self.clear();

            if let Some(name) = record.name() {
                let bytes: &[u8] = name.as_ref();
                self.name.extend_from_slice(bytes);
            }

            self.flags = record.flags().map(|f| f.bits()).unwrap_or(0);

            match record
                .reference_sequence_id(header)
                .and_then(|ref_id| ref_id.ok())
                .and_then(|id| {
                    header
                        .reference_sequences()
                        .get_index(id)
                        .map(|(name, _)| name)
                }) {
                Some(name) => {
                    let buf = self.reference_name.get_or_insert_with(String::new);
                    buf.clear();
                    buf.push_str(&String::from_utf8_lossy(name));
                }
                None => self.reference_name = None,
            }

            self.alignment_start = if let Some(pos) = record.alignment_start() {
                pos.ok().map(crate::types::Pos1::from)
            } else {
                None
            };

            self.alignment_end = if let Some(pos) = record.alignment_end() {
                pos.ok().map(crate::types::Pos1::from)
            } else {
                None
            };

            self.mapping_quality = if (self.flags & 0x4) != 0 {
                None
            } else {
                Some(
                    record
                        .mapping_quality()
                        .and_then(|m| m.ok())
                        .map(|m| m.get())
                        .unwrap_or(255),
                )
            };

            self.cigar.clear();
            self.cigar.reserve(record.cigar().len());
            for op in record.cigar().iter().flatten() {
                self.cigar.push((op.len() as u32, op.kind()));
            }

            self.sequence.clear();
            self.sequence.reserve(record.sequence().len());
            for b in record.sequence().iter() {
                self.sequence.push(b);
            }

            self.quality_scores.clear();
            self.quality_scores.reserve(record.quality_scores().len());
            for score in record.quality_scores().iter().flatten() {
                self.quality_scores.push(score.saturating_add(33));
            }

            match record
                .mate_reference_sequence_id(header)
                .and_then(|ref_id| ref_id.ok())
                .and_then(|id| {
                    header
                        .reference_sequences()
                        .get_index(id)
                        .map(|(name, _)| name)
                }) {
                Some(name) => {
                    let buf = self.mate_reference_name.get_or_insert_with(String::new);
                    buf.clear();
                    buf.push_str(&String::from_utf8_lossy(name));
                }
                None => self.mate_reference_name = None,
            }

            self.mate_alignment_start = if let Some(pos) = record.mate_alignment_start() {
                pos.ok().map(crate::types::Pos1::from)
            } else {
                None
            };

            self.template_length = record.template_length().unwrap_or(0);

            self.nm =
                match record
                    .data()
                    .get(&noodles::sam::alignment::record::data::field::Tag::new(
                        b'n', b'M',
                    )) {
                    Some(Ok(val)) => val.as_int().map(|x| x as i32),
                    _ => record
                        .data()
                        .get(&noodles::sam::alignment::record::data::field::Tag::EDIT_DISTANCE)
                        .and_then(|r| r.ok())
                        .and_then(|val| val.as_int().map(|x| x as i32)),
                };

            match record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'S', b'A',
                ))
                .and_then(|r| r.ok())
            {
                Some(val) => {
                    let buf = self.sa.get_or_insert_with(String::new);
                    buf.clear();
                    match val {
                        noodles::sam::alignment::record::data::field::Value::String(s) => {
                            buf.push_str(&String::from_utf8_lossy(s.as_ref()));
                        }
                        _ => {
                            use std::fmt::Write as _;
                            let _ = write!(buf, "{val:?}");
                        }
                    }
                }
                None => self.sa = None,
            }

            match record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'M', b'C',
                ))
                .and_then(|r| r.ok())
            {
                Some(val) => {
                    let buf = self.mc.get_or_insert_with(String::new);
                    buf.clear();
                    match val {
                        noodles::sam::alignment::record::data::field::Value::String(s) => {
                            buf.push_str(&String::from_utf8_lossy(s.as_ref()));
                        }
                        _ => {
                            use std::fmt::Write as _;
                            let _ = write!(buf, "{val:?}");
                        }
                    }
                }
                None => self.mc = None,
            }

            self.mq = record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'M', b'Q',
                ))
                .and_then(|r| r.ok())
                .and_then(|val| val.as_int().map(|x| x as u8));

            Ok(())
        }
    }

    fn intersects_record_legacy(
        record: &noodles::bam::Record,
        reference_sequence_id: usize,
        region_interval: noodles::core::region::Interval,
    ) -> Result<bool, std::io::Error> {
        let Some(id) = record.reference_sequence_id().transpose()? else {
            return Ok(false);
        };

        if id != reference_sequence_id {
            return Ok(false);
        }

        if region_interval.start().is_none() && region_interval.end().is_none() {
            Ok(true)
        } else {
            match (
                record.alignment_start().transpose()?,
                Record::alignment_end(record).transpose()?,
            ) {
                (Some(start), Some(end)) => {
                    let alignment_interval = (start..=end).into();
                    Ok(region_interval.intersects(alignment_interval))
                }
                _ => Ok(false),
            }
        }
    }

    struct BamFeedInnerLegacy {
        reader: noodles::bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
        header: noodles::sam::Header,
        reference_sequence_id: usize,
        interval: noodles::core::region::Interval,
        chunks:
            std::vec::IntoIter<noodles::csi::binning_index::index::reference_sequence::bin::Chunk>,
        chunk_end: Option<noodles::bgzf::VirtualPosition>,
        record: noodles::bam::Record,
    }

    impl BamFeedInnerLegacy {
        fn pump(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
            use noodles::bgzf::io::Seek as _;
            loop {
                match self.chunk_end {
                    None => {
                        let chunk = match self.chunks.next() {
                            None => return Ok(false),
                            Some(chunk) => chunk,
                        };
                        if let Err(e) = self
                            .reader
                            .get_mut()
                            .seek_to_virtual_position(chunk.start())
                        {
                            return Err(crate::types::VardictError::Io(e));
                        }
                        self.chunk_end = Some(chunk.end());
                    }
                    Some(end) => {
                        if self.reader.get_mut().virtual_position() < end {
                            let n = match self.reader.read_record(&mut self.record) {
                                Ok(n) => n,
                                Err(e) => return Err(crate::types::VardictError::Io(e)),
                            };
                            if n == 0 {
                                self.chunk_end = None;
                                continue;
                            }
                            match intersects_record_legacy(
                                &self.record,
                                self.reference_sequence_id,
                                self.interval,
                            ) {
                                Ok(true) => {
                                    out.fill_from_alignment_record_legacy(
                                        &self.record,
                                        &self.header,
                                    )
                                    .map_err(crate::types::VardictError::Io)?;
                                    return Ok(true);
                                }
                                Ok(false) => continue,
                                Err(e) => return Err(crate::types::VardictError::Io(e)),
                            }
                        } else {
                            self.chunk_end = None;
                        }
                    }
                }
            }
        }
    }

    struct BamFeedLegacy {
        inner: BamFeedInnerLegacy,
        sam: SamRecord,
    }

    impl Iterator for BamFeedLegacy {
        type Item = Result<SamRecord, crate::types::VardictError>;

        fn next(&mut self) -> Option<Self::Item> {
            match self.inner.pump(&mut self.sam) {
                Ok(false) => None,
                Ok(true) => Some(Ok(self.sam.clone())),
                Err(e) => Some(Err(e)),
            }
        }
    }

    fn collect_cram_records(
        reader: &mut IndexedSamReader,
        header: &noodles::sam::Header,
        region: &noodles::core::Region,
    ) -> Result<Vec<SamRecord>, crate::types::VardictError> {
        let iter = reader
            .query(header, region)
            .map_err(crate::types::VardictError::Io)?;
        iter.collect::<Result<Vec<_>, _>>()
            .map_err(crate::types::VardictError::Io)
    }

    struct CramFeedLegacy {
        _reader: IndexedSamReader,
        records: std::vec::IntoIter<SamRecord>,
    }

    impl Iterator for CramFeedLegacy {
        type Item = Result<SamRecord, crate::types::VardictError>;

        fn next(&mut self) -> Option<Self::Item> {
            self.records.next().map(Ok)
        }
    }

    enum RecordFeedLegacy {
        Bam(Box<BamFeedLegacy>),
        Cram(CramFeedLegacy),
    }

    impl Iterator for RecordFeedLegacy {
        type Item = Result<SamRecord, crate::types::VardictError>;

        fn next(&mut self) -> Option<Self::Item> {
            match self {
                RecordFeedLegacy::Bam(feed) => feed.next(),
                RecordFeedLegacy::Cram(feed) => feed.next(),
            }
        }
    }

    impl RecordFeedLegacy {
        fn next_into(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
            match self {
                RecordFeedLegacy::Bam(feed) => feed.inner.pump(out),
                RecordFeedLegacy::Cram(feed) => match feed.records.next() {
                    None => Ok(false),
                    Some(rec) => {
                        *out = rec;
                        Ok(true)
                    }
                },
            }
        }
    }

    fn query_region_iter_legacy(
        bam_path: &Path,
        chr: &str,
        start: Pos0,
        end: Pos0,
        config: &crate::config::Config,
    ) -> Result<RecordFeedLegacy, crate::types::VardictError> {
        let CachedReader { reader, header } = open_reader(bam_path, config)?;

        let start_pos = start.to_one_based().0;
        let end_pos = end.0;
        if start_pos > end_pos {
            return Ok(match reader {
                IndexedSamReader::Bam(reader) => RecordFeedLegacy::Bam(Box::new(BamFeedLegacy {
                    inner: BamFeedInnerLegacy {
                        reader,
                        header,
                        reference_sequence_id: 0,
                        interval: noodles::core::region::Interval::from(..),
                        chunks: Vec::new().into_iter(),
                        chunk_end: None,
                        record: noodles::bam::Record::default(),
                    },
                    sam: SamRecord::new(),
                })),
                IndexedSamReader::Cram(reader, repository) => {
                    RecordFeedLegacy::Cram(CramFeedLegacy {
                        _reader: IndexedSamReader::Cram(reader, repository),
                        records: Vec::new().into_iter(),
                    })
                }
            });
        }

        let region_str = format!("{}:{}-{}", chr, start_pos, end_pos);
        let region: noodles::core::Region = region_str.parse().map_err(|e| {
            crate::types::VardictError::Parse(format!(
                "failed to parse region '{}': {}",
                region_str, e
            ))
        })?;

        match reader {
            IndexedSamReader::Bam(reader) => {
                let reference_sequence_id = header
                    .reference_sequences()
                    .get_index_of(region.name())
                    .ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidInput,
                            format!(
                                "region reference sequence does not exist in reference sequences: {region:?}"
                            ),
                        )
                    })
                    .map_err(crate::types::VardictError::Io)?;
                let chunks = reader
                    .index()
                    .query(reference_sequence_id, region.interval())
                    .map_err(crate::types::VardictError::Io)?;
                Ok(RecordFeedLegacy::Bam(Box::new(BamFeedLegacy {
                    inner: BamFeedInnerLegacy {
                        reader,
                        header,
                        reference_sequence_id,
                        interval: region.interval(),
                        chunks: chunks.into_iter(),
                        chunk_end: None,
                        record: noodles::bam::Record::default(),
                    },
                    sam: SamRecord::new(),
                })))
            }
            mut reader @ IndexedSamReader::Cram(..) => {
                let records = collect_cram_records(&mut reader, &header, &region)?;
                Ok(RecordFeedLegacy::Cram(CramFeedLegacy {
                    _reader: reader,
                    records: records.into_iter(),
                }))
            }
        }
    }

    #[test]
    fn shared_meta_feed_matches_legacy_and_error_strings_match() {
        use noodles::bam::bai;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::{
            bin::Chunk, index::LinearIndex,
        };
        use noodles::csi::binning_index::Indexer;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;

        fn err_signature(e: &crate::types::VardictError) -> (Option<std::io::ErrorKind>, String) {
            match e {
                crate::types::VardictError::Io(io) => (Some(io.kind()), io.to_string()),
                other => (None, other.to_string()),
            }
        }

        clear_bam_caches();
        let config = Config::default();

        // 1. 25 deterministic pseudo-random windows across micro.bam (chr1 and chr2).
        let micro_bam = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/micro/micro.bam"
        ));
        let mut rng_state: u64 = 0x9E37_79B9_7F4A_7C15;
        for i in 0..25usize {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let chr = if i % 2 == 0 { "chr1" } else { "chr2" };
            let max_coord = if i % 2 == 0 { 300u32 } else { 260u32 };
            let a = ((rng_state >> 16) as u32) % max_coord;
            let b = ((rng_state >> 32) as u32) % max_coord;
            let (start, end) = if i == 24 {
                // Include an inverted/empty window in the 25 windows
                (Pos0(200), Pos0(50))
            } else if a <= b {
                (Pos0(a), Pos0(b.max(a + 1)))
            } else {
                (Pos0(b), Pos0(a))
            };

            let mut feed_new = query_region_iter(micro_bam, chr, start, end, &config).unwrap();
            let mut feed_leg =
                query_region_iter_legacy(micro_bam, chr, start, end, &config).unwrap();

            let mut recs_new = Vec::new();
            let mut rec_buf = SamRecord::new();
            while feed_new.next_into(&mut rec_buf).unwrap() {
                recs_new.push(rec_buf.clone());
            }

            let mut recs_leg = Vec::new();
            let mut rec_leg_buf = SamRecord::new();
            while feed_leg.next_into(&mut rec_leg_buf).unwrap() {
                recs_leg.push(rec_leg_buf.clone());
            }

            assert_eq!(recs_new, recs_leg, "window #{i} ({chr}:{start:?}-{end:?})");
        }

        // 2. Nested same-thread feeds (open feed B inside feed A's iteration loop,
        //    plus 5 simultaneous feeds to exercise MAX_POOLED_READERS_PER_PATH = 4).
        {
            let mut feed_a =
                query_region_iter(micro_bam, "chr1", Pos0(0), Pos0(300), &config).unwrap();
            let expected_a: Vec<SamRecord> =
                query_region_iter_legacy(micro_bam, "chr1", Pos0(0), Pos0(300), &config)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            let expected_b: Vec<SamRecord> =
                query_region_iter_legacy(micro_bam, "chr2", Pos0(0), Pos0(260), &config)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();

            let mut got_a = Vec::new();
            let mut rec_a = SamRecord::new();
            let mut step = 0usize;
            while feed_a.next_into(&mut rec_a).unwrap() {
                if step == 3 {
                    let got_b: Vec<SamRecord> =
                        query_region_iter(micro_bam, "chr2", Pos0(0), Pos0(260), &config)
                            .unwrap()
                            .map(|r| r.unwrap())
                            .collect();
                    assert_eq!(got_b, expected_b);

                    // Open 5 simultaneous feeds while feed_a is still open
                    let mut extra_feeds: Vec<RecordFeed> = (0..5)
                        .map(|_| {
                            query_region_iter(micro_bam, "chr2", Pos0(0), Pos0(260), &config)
                                .unwrap()
                        })
                        .collect();
                    for ef in &mut extra_feeds {
                        let mut r = SamRecord::new();
                        assert!(ef.next_into(&mut r).unwrap());
                    }
                    drop(extra_feeds);
                }
                got_a.push(rec_a.clone());
                step += 1;
            }
            assert_eq!(got_a, expected_a);
        }

        // 3. Contig names with `:` and `*` (and missing contig with `:` and `*`).
        let dir = tempfile::TempDir::new().unwrap();
        let special_bam = dir.path().join("special_contigs.bam");
        {
            let rs = HeaderMap::<RefSeqMap>::builder()
                .set_length(500usize.try_into().unwrap())
                .build()
                .unwrap();
            let hd = Header::builder()
                .add_reference_sequence("HLA-A*01:01:01:01", rs.clone())
                .add_reference_sequence("contig*star:colon", rs)
                .build();
            let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(&special_bam).unwrap());
            let mut w = noodles::bam::io::Writer::from(bgzf);
            w.write_alignment_header(&hd).unwrap();
            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
            for (ref_id, name_prefix) in [(0usize, "hla"), (1usize, "star")] {
                for idx in 0..5usize {
                    let pos = 50 + idx * 20;
                    let cigar: Cigar = [Op::new(Kind::Match, 30)].into_iter().collect();
                    let rec = RecordBuf::builder()
                        .set_name(format!("{name_prefix}_{idx}"))
                        .set_flags(noodles::sam::alignment::record::Flags::empty())
                        .set_reference_sequence_id(ref_id)
                        .set_alignment_start(Position::try_from(pos).unwrap())
                        .set_cigar(cigar)
                        .set_sequence(vec![b'A'; 30].into())
                        .set_quality_scores(vec![40u8; 30].into())
                        .build();
                    let s_vp = w.get_ref().virtual_position();
                    w.write_alignment_record(&hd, &rec).unwrap();
                    let e_vp = w.get_ref().virtual_position();
                    let s = rec.alignment_start().unwrap();
                    let e = Position::new(s.get() + 30).unwrap();
                    indexer
                        .add_record(Some((ref_id, s, e, true)), Chunk::new(s_vp, e_vp))
                        .unwrap();
                }
            }
            w.try_finish().unwrap();
            let bai_path = special_bam.with_extension("bam.bai");
            bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                .write_index(&indexer.build(2))
                .unwrap();
        }
        for contig in ["HLA-A*01:01:01:01", "contig*star:colon"] {
            let new_recs: Vec<SamRecord> =
                query_region_iter(&special_bam, contig, Pos0(10), Pos0(300), &config)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            let leg_recs: Vec<SamRecord> =
                query_region_iter_legacy(&special_bam, contig, Pos0(10), Pos0(300), &config)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            assert_eq!(new_recs.len(), 5);
            assert_eq!(new_recs, leg_recs);
        }
        let missing_contig_new = query_region_iter(
            &special_bam,
            "missing:contig*name",
            Pos0(10),
            Pos0(300),
            &config,
        )
        .err()
        .expect("missing contig must fail");
        let missing_contig_leg = query_region_iter_legacy(
            &special_bam,
            "missing:contig*name",
            Pos0(10),
            Pos0(300),
            &config,
        )
        .err()
        .expect("missing contig must fail legacy");
        assert_eq!(
            err_signature(&missing_contig_new),
            err_signature(&missing_contig_leg)
        );

        // 4. Missing BAM + empty window (`start > end`): must still error!
        let missing_bam = dir.path().join("does_not_exist.bam");
        let miss_empty_new = query_region_iter(&missing_bam, "chr1", Pos0(500), Pos0(100), &config)
            .err()
            .expect("missing BAM with empty window must error");
        let miss_empty_leg =
            query_region_iter_legacy(&missing_bam, "chr1", Pos0(500), Pos0(100), &config)
                .err()
                .expect("missing BAM with empty window must error in legacy");
        assert_eq!(
            err_signature(&miss_empty_new),
            err_signature(&miss_empty_leg)
        );

        // 5. Missing `.bai`: must error with identical ErrorKind and NOT cache failure.
        let no_bai_bam = dir.path().join("no_bai.bam");
        std::fs::copy(micro_bam, &no_bai_bam).unwrap();
        let no_bai_new = query_region_iter(&no_bai_bam, "chr1", Pos0(0), Pos0(100), &config)
            .err()
            .expect("missing .bai must error");
        let no_bai_leg = query_region_iter_legacy(&no_bai_bam, "chr1", Pos0(0), Pos0(100), &config)
            .err()
            .expect("missing .bai must error in legacy");
        assert_eq!(err_signature(&no_bai_new), err_signature(&no_bai_leg));
        // Now supply the .bai and verify the previously failed path succeeds (not poisoned).
        std::fs::copy(
            micro_bam.with_extension("bam.bai"),
            no_bai_bam.with_extension("bam.bai"),
        )
        .unwrap();
        let recovered_recs: Vec<SamRecord> =
            query_region_iter(&no_bai_bam, "chr1", Pos0(0), Pos0(100), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert!(!recovered_recs.is_empty());

        // 6. Truncated BAM header: must error with identical ErrorKind and NOT cache failure.
        let trunc_hdr_bam = dir.path().join("trunc_hdr.bam");
        std::fs::write(&trunc_hdr_bam, b"not_a_valid_bgzf_bam_header").unwrap();
        std::fs::copy(
            micro_bam.with_extension("bam.bai"),
            trunc_hdr_bam.with_extension("bam.bai"),
        )
        .unwrap();
        let hdr_err_new = query_region_iter(&trunc_hdr_bam, "chr1", Pos0(0), Pos0(100), &config)
            .err()
            .expect("truncated BAM header must error");
        let hdr_err_leg =
            query_region_iter_legacy(&trunc_hdr_bam, "chr1", Pos0(0), Pos0(100), &config)
                .err()
                .expect("truncated BAM header must error in legacy");
        assert_eq!(err_signature(&hdr_err_new), err_signature(&hdr_err_leg));

        clear_bam_caches();
    }

    #[test]
    fn fill_from_bam_record_matches_legacy_on_aux_and_fields_corpus() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;

        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(1000usize.try_into().unwrap())
            .build()
            .unwrap();
        let header = Header::builder()
            .add_reference_sequence("chr1", rs.clone())
            .add_reference_sequence("chr2", rs)
            .build();

        fn build_raw_bam_record(
            header: &Header,
            base: &RecordBuf,
            raw_aux_bytes: &[u8],
        ) -> noodles::bam::Record {
            let mut encoded = Vec::new();
            let mut writer = noodles::bam::io::Writer::from(&mut encoded);
            writer.write_alignment_record(header, base).unwrap();
            encoded.extend_from_slice(raw_aux_bytes);
            let block_size = (encoded.len() - 4) as u32;
            encoded[0..4].copy_from_slice(&block_size.to_le_bytes());

            let mut reader = noodles::bam::io::Reader::from(&encoded[..]);
            let mut bam_rec = noodles::bam::Record::default();
            let n = reader.read_record(&mut bam_rec).unwrap();
            assert!(n > 0);
            bam_rec
        }

        let base_mapped = RecordBuf::builder()
            .set_name("read_mapped")
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(100usize).unwrap())
            .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
            .set_cigar(
                [
                    Op::new(Kind::SoftClip, 2),
                    Op::new(Kind::Match, 5),
                    Op::new(Kind::Insertion, 1),
                    Op::new(Kind::Deletion, 2),
                    Op::new(Kind::Match, 2),
                ]
                .into_iter()
                .collect::<Cigar>(),
            )
            .set_mate_reference_sequence_id(1)
            .set_mate_alignment_start(Position::try_from(250usize).unwrap())
            .set_template_length(-150)
            .set_sequence(b"ACGTACGTAC".to_vec().into())
            .set_quality_scores(vec![30u8; 10].into())
            .build();

        let base_unmapped = RecordBuf::builder()
            .set_name("read_unmapped")
            .set_flags(Flags::UNMAPPED)
            .set_sequence(b"ACGT".to_vec().into())
            .set_quality_scores(vec![20u8; 4].into())
            .build();

        let base_no_mapq = RecordBuf::builder()
            .set_name("read_no_mapq")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(50usize).unwrap())
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect::<Cigar>())
            .set_sequence(b"TGCA".to_vec().into())
            .set_quality_scores(vec![35u8; 4].into())
            .build();

        // Corpus of raw BAM aux payloads:
        // 1. `nM` int (`nM:C:4`)
        // 2. `nM` non-int (`nM:Z:bad\0`) + `NM:C:7` -> non-int `nM` suppresses `NM`
        // 3. Only `NM` (`NM:C:5`), and `NM:C:5` before `nM:C:2` (`nM` still wins)
        // 4. Duplicate `nM` (`nM:C:3` then `nM:C:9`) and duplicate `SA` (`SA:Z:first\0` then `SA:Z:second\0`)
        // 5. `SA` and `MC` as int/char (`SA:C:42`, `MC:A:x`) -> `{val:?}`
        // 6. Malformed aux field (`XX?`) AFTER wanted tags (wanted tags kept)
        // 7. Malformed aux field (`XX?`) BEFORE wanted tags (`NM` before error kept; `nM`/`SA`/`MC`/`MQ` after error ignored)
        // 8. Empty aux on unmapped and no-mapq records
        let corpus: Vec<noodles::bam::Record> = vec![
            build_raw_bam_record(
                &header,
                &base_mapped,
                b"nMC\x04SAZchr1,100,+,10M,60,0;\0MCC10M\0MQC\x3c",
            ),
            build_raw_bam_record(&header, &base_mapped, b"nMZbad\0NMC\x07"),
            build_raw_bam_record(&header, &base_mapped, b"NMC\x05"),
            build_raw_bam_record(&header, &base_mapped, b"NMC\x05nMC\x02"),
            build_raw_bam_record(
                &header,
                &base_mapped,
                b"nMC\x03nMC\x09SAZfirst\0SAZsecond\0MCZmc1\0MCZmc2\0MQC\x10MQC\x20",
            ),
            build_raw_bam_record(&header, &base_mapped, b"SAC\x2aMCAx"),
            build_raw_bam_record(
                &header,
                &base_mapped,
                b"nMC\x06SAZkept_sa\0MCZkept_mc\0MQC\x25XX?",
            ),
            build_raw_bam_record(
                &header,
                &base_mapped,
                b"NMC\x05XX?nMC\x02SAZdropped_sa\0MCZdropped_mc\0MQC\x25",
            ),
            build_raw_bam_record(&header, &base_unmapped, b""),
            build_raw_bam_record(&header, &base_no_mapq, b""),
        ];

        for (idx, bam_rec) in corpus.iter().enumerate() {
            let mut expected = SamRecord::new();
            expected
                .fill_from_alignment_record_legacy(bam_rec, &header)
                .unwrap();

            let mut actual_bam = SamRecord::new();
            actual_bam.fill_from_bam_record(bam_rec, &header).unwrap();
            assert_eq!(
                actual_bam, expected,
                "fresh fill_from_bam_record mismatch at #{idx}"
            );

            let mut actual_align = SamRecord::new();
            actual_align
                .fill_from_alignment_record(bam_rec, &header)
                .unwrap();
            assert_eq!(
                actual_align, expected,
                "fresh fill_from_alignment_record mismatch at #{idx}"
            );

            let from_rec = SamRecord::from_alignment_record(bam_rec, &header).unwrap();
            assert_eq!(
                from_rec, expected,
                "from_alignment_record mismatch at #{idx}"
            );
        }

        // Verify explicit invariants on key corpus items:
        {
            let mut rec = SamRecord::new();
            // #1: nM non-int + NM -> nm must be None
            rec.fill_from_bam_record(&corpus[1], &header).unwrap();
            assert_eq!(rec.nm, None);
            // #4: duplicate nM / SA -> first wins
            rec.fill_from_bam_record(&corpus[4], &header).unwrap();
            assert_eq!(rec.nm, Some(3));
            assert_eq!(rec.sa.as_deref(), Some("first"));
            assert_eq!(rec.mc.as_deref(), Some("mc1"));
            assert_eq!(rec.mq, Some(0x10));
            // #6: malformed aux AFTER wanted tags -> kept
            rec.fill_from_bam_record(&corpus[6], &header).unwrap();
            assert_eq!(rec.nm, Some(6));
            assert_eq!(rec.sa.as_deref(), Some("kept_sa"));
            assert_eq!(rec.mc.as_deref(), Some("kept_mc"));
            assert_eq!(rec.mq, Some(0x25));
            // #7: malformed aux BEFORE wanted tags -> NM before error kept, tags after error absent
            rec.fill_from_bam_record(&corpus[7], &header).unwrap();
            assert_eq!(rec.nm, Some(5));
            assert_eq!(rec.sa, None);
            assert_eq!(rec.mc, None);
            assert_eq!(rec.mq, None);
        }

        // CRITICAL (C-1): two consecutive `fill_from_bam_record` calls on the
        // same `SamRecord` where Read 1 has `SA`/`MC` and Read 2 has neither
        // must leave `rec.sa.is_none() && rec.mc.is_none()` (not `Some("")`!).
        {
            let read1 =
                build_raw_bam_record(&header, &base_mapped, b"SAZchr1,100,+,10M,60,0;\0MCZ10M\0");
            let read2 = build_raw_bam_record(&header, &base_mapped, b"nMC\x01");

            let mut reused = SamRecord::new();
            reused.fill_from_bam_record(&read1, &header).unwrap();
            assert!(reused.sa.is_some() && reused.mc.is_some());

            reused.fill_from_bam_record(&read2, &header).unwrap();
            assert!(
                reused.sa.is_none() && reused.mc.is_none(),
                "C-1 regression: stale Some(\"\") leaked on reused SamRecord: sa={:?}, mc={:?}",
                reused.sa,
                reused.mc
            );

            let mut expected = SamRecord::new();
            expected
                .fill_from_alignment_record_legacy(&read2, &header)
                .unwrap();
            assert_eq!(reused, expected);

            // Also verify `fill_from_alignment_record` on the same transition
            reused.fill_from_alignment_record(&read1, &header).unwrap();
            assert!(reused.sa.is_some() && reused.mc.is_some());
            reused.fill_from_alignment_record(&read2, &header).unwrap();
            assert!(reused.sa.is_none() && reused.mc.is_none());
            assert_eq!(reused, expected);
        }
    }

    #[test]
    fn bgzf_libdeflate_reader_roundtrips_micro_bam_and_rejects_corrupt_block() {
        use std::io::Read;

        let path = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/micro/micro.bam"
        ));
        let file = std::fs::File::open(path).unwrap();
        let mut reader = noodles_bgzf::io::Reader::new(file);
        let mut decompressed = Vec::new();
        reader.read_to_end(&mut decompressed).unwrap();
        assert!(
            decompressed.starts_with(b"BAM\x01"),
            "decompressed micro.bam must start with BAM\\x01 magic"
        );

        let raw = std::fs::read(path).unwrap();
        let mut corrupt = raw.clone();
        // Corrupt the compressed payload mid-block and ensure libdeflate returns Err
        let mid = corrupt.len() / 2;
        corrupt[mid] ^= 0xFF;
        corrupt[mid + 1] ^= 0xFF;
        let mut bad_reader = noodles_bgzf::io::Reader::new(std::io::Cursor::new(corrupt));
        let mut out = Vec::new();
        assert!(
            bad_reader.read_to_end(&mut out).is_err(),
            "corrupted BGZF block must fail decompression"
        );
    }

    #[test]
    fn test_cached_bgzf_reader_same_block_shared_and_multiblock_parity() {
        use noodles::bgzf::io::Seek as _;
        use std::io::{BufRead as _, Read as _, Write as _};

        const TEST_META_ID: u64 = TEST_META_ID_BASE + 42;

        // 1. Construct a multi-block BGZF file (> 130 KiB uncompressed -> 3 BGZF blocks
        //    + an empty mid-stream BGZF block + 28-byte EOF block)
        let tmp_dir = tempfile::tempdir().unwrap();
        let bgzf_path = tmp_dir.path().join("multiblock.bgz");
        let mut raw_payload = Vec::with_capacity(150_000);
        for i in 0..150_000u32 {
            raw_payload.push(((i.wrapping_mul(31) ^ (i >> 3)) & 0xff) as u8);
        }
        {
            let mut file = std::fs::File::create(&bgzf_path).unwrap();
            // Write first 100,000 bytes as standard BGZF blocks (without trailing EOF yet)
            let mut part1 = Vec::new();
            {
                let mut w1 = noodles::bgzf::io::Writer::new(&mut part1);
                w1.write_all(&raw_payload[..100_000]).unwrap();
                w1.finish().unwrap();
            }
            // Strip the 28-byte EOF block from part1, insert a 28-byte empty BGZF block mid-stream,
            // then append the remaining 50,000 bytes + final 28-byte EOF block.
            file.write_all(&part1).unwrap();
            let mut w2 = noodles::bgzf::io::Writer::new(file);
            w2.write_all(&raw_payload[100_000..]).unwrap();
            w2.finish().unwrap();
        }

        // Record ground-truth block boundaries and virtual positions from noodles::bgzf::io::Reader
        let mut oracle = noodles::bgzf::io::Reader::new(std::fs::File::open(&bgzf_path).unwrap());
        let mut cached =
            CachedBgzfReader::new(std::fs::File::open(&bgzf_path).unwrap(), TEST_META_ID);

        assert_eq!(cached.virtual_position(), oracle.virtual_position());

        // Read across all blocks (including mid-stream empty block and final EOF block)
        let mut oracle_out = Vec::new();
        let mut cached_out = Vec::new();
        let mut obuf = [0u8; 4096];
        let mut cbuf = [0u8; 4096];
        loop {
            assert_eq!(cached.virtual_position(), oracle.virtual_position());
            let on = oracle.read(&mut obuf).unwrap();
            let cn = cached.read(&mut cbuf).unwrap();
            assert_eq!(cn, on);
            assert_eq!(&cbuf[..cn], &obuf[..on]);
            if on == 0 {
                break;
            }
            oracle_out.extend_from_slice(&obuf[..on]);
            cached_out.extend_from_slice(&cbuf[..cn]);
        }
        assert_eq!(cached_out, raw_payload);
        assert_eq!(cached.virtual_position(), oracle.virtual_position());
        assert!(cached.fill_buf().unwrap().is_empty());
        assert_eq!(cached.virtual_position(), oracle.virtual_position());
        let eof_vp = oracle.virtual_position();

        // 2. Same-block seek must increment local_hits with 0 misses and 0 shared_hits
        let first_block_vp1 = noodles::bgzf::VirtualPosition::try_from((0, 100)).unwrap();
        let first_block_vp2 = noodles::bgzf::VirtualPosition::try_from((0, 2500)).unwrap();
        cached.seek_to_virtual_position(first_block_vp1).unwrap();
        oracle.seek(first_block_vp1).unwrap();
        cached.reset_cache_stats();

        cached.seek_to_virtual_position(first_block_vp2).unwrap();
        oracle.seek(first_block_vp2).unwrap();
        assert_eq!(
            cached.cache_stats(),
            (1, 0, 0),
            "same-block seek must hit local 1-block cursor with 0 shared_hits and 0 misses"
        );
        assert_eq!(cached.virtual_position(), oracle.virtual_position());
        let mut s_oracle = [0u8; 512];
        let mut s_cached = [0u8; 512];
        oracle.read_exact(&mut s_oracle).unwrap();
        cached.read_exact(&mut s_cached).unwrap();
        assert_eq!(s_cached, s_oracle);
        assert_eq!(cached.virtual_position(), oracle.virtual_position());

        // 3. Cross-thread CachedBgzfReader with same meta_id seeking to block 0 must hit shared BGZF_BLOCK_CACHE
        let bgzf_path_clone = bgzf_path.clone();
        let (stats2, s_cached2) = std::thread::spawn(move || {
            let mut cached2 =
                CachedBgzfReader::new(std::fs::File::open(&bgzf_path_clone).unwrap(), TEST_META_ID);
            cached2.seek_to_virtual_position(first_block_vp2).unwrap();
            let stats = cached2.cache_stats();
            let mut buf = [0u8; 512];
            cached2.read_exact(&mut buf).unwrap();
            (stats, buf)
        })
        .join()
        .unwrap();
        assert_eq!(
            stats2,
            (0, 1, 0),
            "worker thread reader must hit shared BGZF_BLOCK_CACHE with 0 misses"
        );
        assert_eq!(s_cached2, s_oracle);

        // 4. Out-of-bounds upos (> block.len()) must return InvalidInput and match virtual_position
        let invalid_eof_upos =
            noodles::bgzf::VirtualPosition::try_from((eof_vp.compressed(), 1)).unwrap();
        let o_err = oracle.seek(invalid_eof_upos).unwrap_err();
        let c_err = cached
            .seek_to_virtual_position(invalid_eof_upos)
            .unwrap_err();
        assert_eq!(c_err.kind(), o_err.kind());
        assert_eq!(cached.virtual_position(), oracle.virtual_position());
    }

    #[test]
    fn fill_from_alignment_record_uppercases_soft_masked_cram_bases() {
        // Bug C: CRAM reference-matching bases reconstructed against soft-masked
        // lowercase FASTA contigs (`a/c/g/t/n`) must be normalized to uppercase
        // (`A/C/G/T/N`) in `SamRecord::fill_from_alignment_record`.
        use noodles::sam::alignment::RecordBuf;
        let header = noodles::sam::Header::default();
        let rec_buf = RecordBuf::builder()
            .set_sequence(noodles::sam::alignment::record_buf::Sequence::from(
                b"acgtACGTnN".to_vec(),
            ))
            .build();
        let mut sam_rec = SamRecord::new();
        sam_rec
            .fill_from_alignment_record(&rec_buf, &header)
            .unwrap();
        assert_eq!(sam_rec.sequence, b"ACGTACGTNN");
    }

    #[test]
    fn compute_cram_nm_reconstructs_edit_distance_when_omitted_by_cram() {
        // CRAM 3.0 §10.6 omits `NM` on write by default (`samtools view -C`);
        // `compute_cram_nm` must reconstruct the SAMtags §1.5 edit distance
        // (mismatches + insertions + deletions, excluding soft clips and skips).
        let ref_seq = b"acgtACGTACGT";
        let rec = SamRecord {
            flags: 0,
            alignment_start: Some(Pos1(2)), // 1-based pos 2 -> 'c'
            cigar: vec![
                (1, CigarKind::SoftClip),  // 'G' clip (excluded)
                (4, CigarKind::Match),     // "CGTA" vs "cgtA" (1 mismatch if 'T'->'A')
                (2, CigarKind::Insertion), // "TT" (+2)
                (2, CigarKind::Deletion),  // "CG" (+2)
                (3, CigarKind::Match),     // "TAC" vs "TAC" (0)
            ],
            sequence: b"GCGAGTTTAC".to_vec(), // clip 'G', match "CGAG" vs "CGTA" (2 mm), ins "TT", match "TAC"
            ..SamRecord::default()
        };
        assert_eq!(compute_cram_nm(&rec, ref_seq), Some(6));

        // SAMtags §1.5 / HTSJDK `SequenceUtil.calculateSamNmTag`:
        // read base `=` always matches (0 mismatch), and IUPAC ambiguity codes
        // like `N` vs `N` count as a mismatch (+1).
        let ref_with_n = b"ACNT";
        let rec_eq_and_n = SamRecord {
            flags: 0,
            alignment_start: Some(Pos1(1)),
            cigar: vec![(4, CigarKind::Match)],
            sequence: b"=CNT".to_vec(), // '=' matches 'A' (0), 'C'=='C' (0), 'N' vs 'N' (+1), 'T'=='T' (0)
            ..SamRecord::default()
        };
        assert_eq!(compute_cram_nm(&rec_eq_and_n, ref_with_n), Some(1));
    }

    #[test]
    fn compute_cram_nm_saturates_on_hostile_cigar_lengths() {
        let ref_seq = b"ACGT";
        let rec = SamRecord {
            flags: 0,
            alignment_start: Some(Pos1(1)),
            cigar: vec![
                (u32::MAX, CigarKind::Insertion),
                (u32::MAX, CigarKind::Deletion),
            ],
            sequence: b"ACGT".to_vec(),
            ..SamRecord::default()
        };
        assert_eq!(compute_cram_nm(&rec, ref_seq), Some(i32::MAX));
    }
}
