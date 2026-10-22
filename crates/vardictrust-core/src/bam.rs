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

        // Ports the original if/else-if chain exactly: EDIT_DISTANCE is only
        // consulted when the NM tag is absent or undecodable (a present but
        // non-int NM still suppresses the fallback).
        self.nm = match record
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

    pub fn from_alignment_record<R>(
        record: &R,
        header: &noodles::sam::Header,
    ) -> Result<Self, std::io::Error>
    where
        R: Record,
    {
        let name = if let Some(name) = record.name() {
            let bytes: &[u8] = name.as_ref();
            bytes.to_vec()
        } else {
            Vec::new()
        };

        let flags = record.flags().map(|f| f.bits()).unwrap_or(0);

        let reference_name = record
            .reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            });

        let alignment_start = if let Some(pos) = record.alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        let alignment_end = if let Some(pos) = record.alignment_end() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        let mapping_quality = if (flags & 0x4) != 0 {
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

        let mut cigar = Vec::with_capacity(record.cigar().len());
        for op in record.cigar().iter().flatten() {
            cigar.push((op.len() as u32, op.kind()));
        }

        let mut sequence = Vec::with_capacity(record.sequence().len());
        for b in record.sequence().iter() {
            sequence.push(b);
        }

        let mut quality_scores = Vec::with_capacity(record.quality_scores().len());
        for score in record.quality_scores().iter().flatten() {
            quality_scores.push(score.saturating_add(33));
        }

        let mate_reference_name = record
            .mate_reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            });

        let mate_alignment_start = if let Some(pos) = record.mate_alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        let template_length = record.template_length().unwrap_or(0);

        let mut nm = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'n', b'M',
                ))
        {
            nm = val.as_int().map(|x| x as i32);
        } else if let Some(Ok(val)) = record
            .data()
            .get(&noodles::sam::alignment::record::data::field::Tag::EDIT_DISTANCE)
        {
            nm = val.as_int().map(|x| x as i32);
        }

        let mut sa = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'S', b'A',
                ))
        {
            match val {
                noodles::sam::alignment::record::data::field::Value::String(s) => {
                    sa = Some(String::from_utf8_lossy(s.as_ref()).into_owned());
                }
                _ => {
                    sa = Some(format!("{:?}", val));
                }
            }
        }

        let mut mc = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'M', b'C',
                ))
        {
            match val {
                noodles::sam::alignment::record::data::field::Value::String(s) => {
                    mc = Some(String::from_utf8_lossy(s.as_ref()).into_owned());
                }
                _ => {
                    mc = Some(format!("{:?}", val));
                }
            }
        }

        let mut mq = None;
        if let Some(Ok(val)) =
            record
                .data()
                .get(&noodles::sam::alignment::record::data::field::Tag::new(
                    b'M', b'Q',
                ))
        {
            mq = val.as_int().map(|x| x as u8);
        }

        Ok(SamRecord {
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
        })
    }
}

pub enum IndexedSamReader {
    Bam(noodles::bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>),
    Cram(noodles::cram::io::IndexedReader<std::fs::File>),
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
            IndexedSamReader::Cram(reader) => {
                let query = reader.query(header, region)?;
                let iter = query.records().map(|record| {
                    let rec = record?;
                    SamRecord::from_alignment_record(&rec, header)
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
            .set_reference_sequence_repository(repository)
            .build_from_path(path)
            .map_err(crate::types::VardictError::Io)?;
        let header = reader
            .read_header()
            .map_err(crate::types::VardictError::Io)?;
        Ok(CachedReader {
            reader: IndexedSamReader::Cram(reader),
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

/// Owning lazy record feed over one queried region (Task 4: OOM fix).
///
/// `noodles`' indexed `Query` borrows both the reader and the header, so it
/// cannot live in the same struct that owns them — and the `thread_local!
/// SAM_READER_CACHE` guard is `!Send`, so the per-tile reader is opened
/// fresh (non-cached) inside the rayon worker instead. `BamFeed` therefore
/// replicates the indexed-query walk with owned state: the same
/// `BinningIndex::query` chunk list, the same csi Seek/Read chunk state
/// machine, and the same `intersects` accept/reject shape, visiting exactly
/// the records `IndexedReader::query` would, one at a time.
pub struct BamFeed {
    inner: BamFeedInner,
    sam: SamRecord,
}

struct BamFeedInner {
    reader: noodles::bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: noodles::sam::Header,
    reference_sequence_id: usize,
    interval: noodles::core::region::Interval,
    chunks: std::vec::IntoIter<noodles::csi::binning_index::index::reference_sequence::bin::Chunk>,
    chunk_end: Option<noodles::bgzf::VirtualPosition>,
    record: noodles::bam::Record,
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
            record.alignment_end().transpose()?,
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
                        match intersects_record(
                            &self.record,
                            self.reference_sequence_id,
                            self.interval,
                        ) {
                            Ok(true) => {
                                out.fill_from_alignment_record(&self.record, &self.header)
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

/// CRAM arm of the feed: records are collected once through the dedicated
/// non-cached reader (identical bytes to the cached `query_region` path, so
/// behavior is identical by construction). Only the BAM arm streams lazily —
/// CRAM has no reported OOM case and its container/slice query internals are
/// not replicable without new threads or unsafe.
pub struct CramFeed {
    _reader: IndexedSamReader,
    records: std::vec::IntoIter<SamRecord>,
}

impl Iterator for CramFeed {
    type Item = Result<SamRecord, crate::types::VardictError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next().map(Ok)
    }
}

/// Owning streaming record feed returned by [`query_region_iter`].
pub enum RecordFeed {
    Bam(Box<BamFeed>),
    Cram(CramFeed),
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
            RecordFeed::Cram(feed) => match feed.records.next() {
                None => Ok(false),
                Some(rec) => {
                    *out = rec;
                    Ok(true)
                }
            },
        }
    }
}

/// Streaming counterpart to [`query_region`]: opens a dedicated non-cached
/// reader (never the `SAM_READER_CACHE`, whose guard cannot cross the rayon
/// worker boundary) and returns a lazy [`RecordFeed`] yielding one
/// `SamRecord` at a time. I/O errors surface per record and propagate with
/// `?` exactly like the Vec path's collection errors; `start > end` yields
/// an empty feed (preserving the `query_region` short-circuit).
pub fn query_region_iter(
    bam_path: &Path,
    chr: &str,
    start: Pos0,
    end: Pos0,
    config: &crate::config::Config,
) -> Result<RecordFeed, crate::types::VardictError> {
    let CachedReader { reader, header } = open_reader(bam_path, config)?;

    // Preserve the query_region empty-window short-circuit: no filesystem
    // reads and no region parse, just an empty feed.
    let start_pos = start.to_one_based().0;
    let end_pos = end.0;
    if start_pos > end_pos {
        return Ok(match reader {
            IndexedSamReader::Bam(reader) => RecordFeed::Bam(Box::new(BamFeed {
                inner: BamFeedInner {
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
            IndexedSamReader::Cram(reader) => RecordFeed::Cram(CramFeed {
                _reader: IndexedSamReader::Cram(reader),
                records: Vec::new().into_iter(),
            }),
        });
    }

    let region_str = format!("{}:{}-{}", chr, start_pos, end_pos);
    let region: noodles::core::Region = region_str.parse().map_err(|e| {
        crate::types::VardictError::Parse(format!("failed to parse region '{}': {}", region_str, e))
    })?;

    match reader {
        IndexedSamReader::Bam(reader) => {
            // Replicates IndexedReader::query's setup (resolve_region +
            // BinningIndex::query) with owned results; the borrow ends here.
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
            Ok(RecordFeed::Bam(Box::new(BamFeed {
                inner: BamFeedInner {
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
        IndexedSamReader::Cram(mut reader) => {
            let query = reader
                .query(&header, &region)
                .map_err(crate::types::VardictError::Io)?;
            let mut records = Vec::new();
            for record in query.records() {
                let rec = record.map_err(crate::types::VardictError::Io)?;
                records.push(SamRecord::from_alignment_record(&rec, &header)?);
            }
            Ok(RecordFeed::Cram(CramFeed {
                _reader: IndexedSamReader::Cram(reader),
                records: records.into_iter(),
            }))
        }
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
}
