use crate::types::{Pos0, Pos1};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::alignment::record::{QualityScores, Sequence};
use noodles::sam::alignment::Record;
use std::cell::RefCell;
use std::collections::HashMap;
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

/// Ports HTSJDK `SequenceUtil.bamReadBaseLookup` (127-byte lookup table initialized with `'N'`
/// and populated for `"=ABCDGHKMNRSTVWY"` at `b` and `b + 32`).
///
/// Bytes `>= 0x7F` (`127..=255`) are out of bounds in HTSJDK's `byte[127]` table and return `None`.
const CRAM_BASE_LUT: [u8; 127] = build_cram_base_lut();

const fn build_cram_base_lut() -> [u8; 127] {
    let mut table = [b'N'; 127];
    let mut rem: &mut [u8] = &mut table;
    let mut b: u8 = 0;
    while let [slot, rest @ ..] = rem {
        *slot = match b {
            b'=' | b']' => b'=',
            b'A' | b'a' => b'A',
            b'B' | b'b' => b'B',
            b'C' | b'c' => b'C',
            b'D' | b'd' => b'D',
            b'G' | b'g' => b'G',
            b'H' | b'h' => b'H',
            b'K' | b'k' => b'K',
            b'M' | b'm' => b'M',
            b'N' | b'n' => b'N',
            b'R' | b'r' => b'R',
            b'S' | b's' => b'S',
            b'T' | b't' => b'T',
            b'V' | b'v' => b'V',
            b'W' | b'w' => b'W',
            b'Y' | b'y' => b'Y',
            _ => b'N',
        };
        rem = rest;
        b += 1;
    }
    table
}

#[inline]
pub(crate) fn normalize_cram_base(b: u8) -> Option<u8> {
    CRAM_BASE_LUT.get(usize::from(b)).copied()
}

const BAM_SEQ_PAIR_LUT: [[u8; 2]; 256] = build_bam_seq_pair_lut();

const fn decode_bam_nibble(nibble: u8) -> u8 {
    match nibble & 0x0f {
        0 => b'=',
        1 => b'A',
        2 => b'C',
        3 => b'M',
        4 => b'G',
        5 => b'R',
        6 => b'S',
        7 => b'V',
        8 => b'T',
        9 => b'W',
        10 => b'Y',
        11 => b'H',
        12 => b'K',
        13 => b'D',
        14 => b'B',
        _ => b'N',
    }
}

const fn build_bam_seq_pair_lut() -> [[u8; 2]; 256] {
    let mut table = [[b'N'; 2]; 256];
    let mut rem: &mut [[u8; 2]] = &mut table;
    let mut b: usize = 0;
    while let [slot, rest @ ..] = rem {
        let byte = b as u8;
        *slot = [decode_bam_nibble(byte >> 4), decode_bam_nibble(byte & 0x0f)];
        rem = rest;
        b += 1;
    }
    table
}

#[inline]
fn decode_bam_raw_seq_into(raw_seq: &[u8], n: usize, dst: &mut Vec<u8>) {
    dst.clear();
    dst.reserve(n.min(raw_seq.len().saturating_mul(2)));
    let full_pairs = n / 2;
    if let Some(pair_bytes) = raw_seq.get(..full_pairs) {
        for &byte in pair_bytes {
            let pair = BAM_SEQ_PAIR_LUT.get(usize::from(byte)).unwrap_or(b"NN");
            dst.extend_from_slice(pair);
        }
        if (n & 1) == 1 {
            if let Some(&last_byte) = raw_seq.get(full_pairs) {
                let &[hi, _] = BAM_SEQ_PAIR_LUT
                    .get(usize::from(last_byte))
                    .unwrap_or(b"NN");
                dst.push(hi);
            }
        }
    } else {
        for &byte in raw_seq {
            let &[hi, lo] = BAM_SEQ_PAIR_LUT.get(usize::from(byte)).unwrap_or(b"NN");
            if dst.len() < n {
                dst.push(hi);
            }
            if dst.len() < n {
                dst.push(lo);
            }
        }
    }
}

#[inline]
fn decode_bam_sequence_into(seq: noodles::bam::record::Sequence<'_>, dst: &mut Vec<u8>) {
    decode_bam_raw_seq_into(seq.as_ref(), seq.len(), dst);
}

#[inline]
fn decode_bam_quality_scores_into(
    quals: noodles::bam::record::QualityScores<'_>,
    dst: &mut Vec<u8>,
) {
    dst.clear();
    let raw_qual = quals.as_ref();
    if !raw_qual.is_empty() {
        dst.extend_from_slice(raw_qual);
        for q in dst.iter_mut() {
            *q = q.saturating_add(33);
        }
    }
}

const MAX_FEED_STRING_STASH_CAP: usize = 4096;

#[inline]
fn stash_capped_string(stash: &mut String, mut buf: String) {
    if buf.capacity() <= MAX_FEED_STRING_STASH_CAP && buf.capacity() > stash.capacity() {
        buf.clear();
        *stash = buf;
    }
}

#[inline]
fn set_opt_contig_name(dst: &mut Option<String>, expected: Option<&str>) {
    match expected {
        Some(expected) => {
            if dst.as_deref() != Some(expected) {
                let buf = dst.get_or_insert_with(String::new);
                buf.clear();
                buf.push_str(expected);
            }
        }
        None => *dst = None,
    }
}

#[inline]
fn set_opt_contig_name_with_stash(
    dst: &mut Option<String>,
    expected: Option<&str>,
    stash: &mut String,
) {
    match expected {
        Some(expected) => {
            if dst.as_deref() != Some(expected) {
                let buf = dst.get_or_insert_with(|| std::mem::take(stash));
                buf.clear();
                buf.push_str(expected);
            }
        }
        None => {
            if let Some(buf) = dst.take() {
                stash_capped_string(stash, buf);
            }
        }
    }
}

#[inline]
fn checked_cigar_op_len(len: usize) -> Result<u32, std::io::Error> {
    u32::try_from(len).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("CIGAR operator length {len} exceeds u32::MAX"),
        )
    })
}

#[inline]
fn checked_alignment_end(
    start: noodles::core::Position,
    span: usize,
) -> Result<noodles::core::Position, std::io::Error> {
    if span == 0 {
        return Ok(start);
    }
    usize::from(start)
        .checked_add(span)
        .and_then(|v| v.checked_sub(1))
        .and_then(noodles::core::Position::new)
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "alignment span overflow (start={}, span={span})",
                    usize::from(start)
                ),
            )
        })
}

impl SamRecord {
    /// Empty record with no heap allocation; use with
    /// [`SamRecord::fill_from_alignment_record`] for buffer reuse.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset every field so a reused record cannot leak a previous read's
    /// values into the next fill. `sa`/`mc` keep their `String` capacity
    /// (`Some("")`), while `reference_name` and `mate_reference_name` retain
    /// their previous contig string so same-contig refills elide
    /// `clear() + push_str(...)` and absent contigs assign `None`. Must stay
    /// in sync with the field list and the fill functions.
    pub fn clear(&mut self) {
        self.name.clear();
        self.flags = 0;
        self.alignment_start = None;
        self.alignment_end = None;
        self.mapping_quality = None;
        self.cigar.clear();
        self.sequence.clear();
        self.quality_scores.clear();
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
        // borrowed Cow so valid-UTF-8 names allocate nothing, and elide
        // clear+push_str when the contig name is unchanged.
        let ref_name = record
            .reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name))
            });
        set_opt_contig_name(&mut self.reference_name, ref_name.as_deref());

        let raw_start = record.alignment_start().and_then(Result::ok);
        self.alignment_start = raw_start.map(crate::types::Pos1::from);

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
        let mut span = 0usize;
        let mut cigar_err = false;
        for result in cigar.iter() {
            match result {
                Ok(op) => {
                    let op_len = checked_cigar_op_len(op.len())?;
                    if op.kind().consumes_reference() {
                        span = span.checked_add(op.len()).ok_or_else(|| {
                            std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                "CIGAR alignment_span overflow",
                            )
                        })?;
                    }
                    self.cigar.push((op_len, op.kind()));
                }
                Err(_) => {
                    cigar_err = true;
                    break;
                }
            }
        }

        self.alignment_end = match raw_start {
            Some(start) if !cigar_err => Some(crate::types::Pos1::from(checked_alignment_end(
                start, span,
            )?)),
            _ => None,
        };

        let seq = record.sequence();
        self.sequence.clear();
        self.sequence.reserve(seq.len());
        for (idx, b) in seq.iter().enumerate() {
            // Normalize CRAM/SAM bases through HTSJDK `SequenceUtil.bamReadBaseLookup`
            // (ports `CRAMRecord` / `SAMRecord.getReadBases()`).
            let Some(norm) = normalize_cram_base(b) else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "read {}: non-ASCII base 0x{:02x} at {}",
                        String::from_utf8_lossy(&self.name),
                        b,
                        idx + 1,
                    ),
                ));
            };
            self.sequence.push(norm);
        }

        let quals = record.quality_scores();
        self.quality_scores.clear();
        self.quality_scores.reserve(quals.len());
        for score in quals.iter().flatten() {
            self.quality_scores.push(score.saturating_add(33));
        }

        let mate_ref_name = record
            .mate_reference_sequence_id(header)
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name))
            });
        set_opt_contig_name(&mut self.mate_reference_name, mate_ref_name.as_deref());

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

        let ref_name = record
            .reference_sequence_id()
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name))
            });
        set_opt_contig_name(&mut self.reference_name, ref_name.as_deref());

        self.alignment_start = if let Some(pos) = record.alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.alignment_end = match bam_alignment_end(record) {
            Some(Ok(pos)) => Some(crate::types::Pos1::from(pos)),
            Some(Err(e)) => return Err(e),
            None => None,
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
            self.cigar
                .push((checked_cigar_op_len(op.len())?, op.kind()));
        }

        decode_bam_sequence_into(record.sequence(), &mut self.sequence);
        decode_bam_quality_scores_into(record.quality_scores(), &mut self.quality_scores);

        let mate_ref_name = record
            .mate_reference_sequence_id()
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| {
                header
                    .reference_sequences()
                    .get_index(id)
                    .map(|(name, _)| String::from_utf8_lossy(name))
            });
        set_opt_contig_name(&mut self.mate_reference_name, mate_ref_name.as_deref());

        self.mate_alignment_start = if let Some(pos) = record.mate_alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.template_length = record.template_length();

        self.fill_aux_from_iter(record.data().iter());

        Ok(())
    }

    /// Single-pass `BamFeed` refill path: reuses pre-decoded `cigar_scratch`,
    /// `(alignment_start, alignment_end)`, `O(1)` `contig_names` lookup, and
    /// `<= 4096`-byte-capped `sa`/`mc`/`mate_reference_name` string stashes.
    #[expect(
        clippy::too_many_arguments,
        reason = "hot BamFeed refill borrows individual scratch and stash buffers directly"
    )]
    fn fill_from_bam_feed_record(
        &mut self,
        record: &noodles::bam::Record,
        contig_names: &[Box<str>],
        reference_sequence_id: usize,
        alignment_start: Option<noodles::core::Position>,
        alignment_end: Option<noodles::core::Position>,
        cigar_scratch: &mut Vec<(u32, CigarKind)>,
        sa_buf: &mut String,
        mc_buf: &mut String,
        mate_ref_buf: &mut String,
    ) {
        self.clear();

        if let Some(name) = record.name() {
            let bytes: &[u8] = name.as_ref();
            self.name.extend_from_slice(bytes);
        }

        self.flags = record.flags().bits();

        set_opt_contig_name(
            &mut self.reference_name,
            contig_names.get(reference_sequence_id).map(Box::as_ref),
        );

        self.alignment_start = alignment_start.map(crate::types::Pos1::from);
        self.alignment_end = alignment_end.map(crate::types::Pos1::from);

        self.mapping_quality = if (self.flags & 0x4) != 0 {
            None
        } else {
            Some(record.mapping_quality().map(|m| m.get()).unwrap_or(255))
        };

        std::mem::swap(&mut self.cigar, cigar_scratch);

        decode_bam_sequence_into(record.sequence(), &mut self.sequence);
        decode_bam_quality_scores_into(record.quality_scores(), &mut self.quality_scores);

        let mate_contig = record
            .mate_reference_sequence_id()
            .and_then(|ref_id| ref_id.ok())
            .and_then(|id| contig_names.get(id))
            .map(Box::as_ref);
        set_opt_contig_name_with_stash(&mut self.mate_reference_name, mate_contig, mate_ref_buf);

        self.mate_alignment_start = if let Some(pos) = record.mate_alignment_start() {
            pos.ok().map(crate::types::Pos1::from)
        } else {
            None
        };

        self.template_length = record.template_length();

        self.fill_aux_from_iter_with_stashes(record.data().iter(), sa_buf, mc_buf);
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
        self.fill_aux_from_iter_with_stashes(iter, &mut String::new(), &mut String::new());
    }

    fn fill_aux_from_iter_with_stashes<'a, I>(
        &mut self,
        iter: I,
        sa_stash: &mut String,
        mc_stash: &mut String,
    ) where
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
                    nm_lower_val = val
                        .as_int()
                        .map(|x| x.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32);
                }
            } else if tag == TAG_NM_UPPER {
                if !seen_nm_upper {
                    seen_nm_upper = true;
                    nm_upper_val = val
                        .as_int()
                        .map(|x| x.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32);
                }
            } else if tag == TAG_SA {
                if !seen_sa {
                    seen_sa = true;
                    let buf = self.sa.get_or_insert_with(|| std::mem::take(sa_stash));
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
                    let buf = self.mc.get_or_insert_with(|| std::mem::take(mc_stash));
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
                self.mq = val.as_int().and_then(|x| u8::try_from(x).ok());
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
            if let Some(buf) = self.sa.take() {
                stash_capped_string(sa_stash, buf);
            }
        }
        if !seen_mc {
            if let Some(buf) = self.mc.take() {
                stash_capped_string(mc_stash, buf);
            }
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

#[inline]
fn is_strict_nm_match_base(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'='
    )
}

#[inline]
fn is_sam_nm_base_match(read_b: u8, ref_b: u8) -> bool {
    read_b == b'='
        || (matches!(
            read_b,
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
        ) && read_b.eq_ignore_ascii_case(&ref_b))
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
                if let (Some(read_end), Some(ref_end)) =
                    (read_pos.checked_add(len), ref_pos.checked_add(len))
                {
                    if let (Some(read_slice), Some(ref_slice)) = (
                        rec.sequence.get(read_pos..read_end),
                        ref_seq.get(ref_pos..ref_end),
                    ) {
                        if read_slice == ref_slice
                            && read_slice.iter().copied().all(is_strict_nm_match_base)
                        {
                            read_pos = read_end;
                            ref_pos = ref_end;
                            continue;
                        }
                        let mut r_chunks = read_slice.chunks_exact(8);
                        let mut f_chunks = ref_slice.chunks_exact(8);
                        for (r8, f8) in r_chunks.by_ref().zip(f_chunks.by_ref()) {
                            if r8 == f8 && r8.iter().copied().all(is_strict_nm_match_base) {
                                continue;
                            }
                            for (&read_b, &ref_b) in r8.iter().zip(f8.iter()) {
                                if !is_sam_nm_base_match(read_b, ref_b) {
                                    nm = nm.saturating_add(1);
                                }
                            }
                        }
                        for (&read_b, &ref_b) in
                            r_chunks.remainder().iter().zip(f_chunks.remainder().iter())
                        {
                            if !is_sam_nm_base_match(read_b, ref_b) {
                                nm = nm.saturating_add(1);
                            }
                        }
                        read_pos = read_end;
                        ref_pos = ref_end;
                        continue;
                    }
                }
                for _ in 0..len {
                    let Some(&read_b) = rec.sequence.get(read_pos) else {
                        break;
                    };
                    let ref_b = ref_seq.get(ref_pos).copied().unwrap_or(b'N');
                    if !is_sam_nm_base_match(read_b, ref_b) {
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

#[deprecated(since = "1.0.10", note = "use bam::query_region_iter")]
pub enum IndexedSamReader {
    Bam(noodles::bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>),
    Cram(PathBuf, std::sync::Arc<SharedCramMeta>),
}

#[allow(deprecated)]
impl IndexedSamReader {
    pub fn query<'a>(
        &'a mut self,
        header: &'a noodles::sam::Header,
        region: &'a noodles::core::Region,
    ) -> Result<Box<dyn Iterator<Item = Result<SamRecord, std::io::Error>> + 'a>, std::io::Error>
    {
        match self {
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
                    })?;
                let interval = region.interval();
                let chunks = reader.index().query(reference_sequence_id, interval)?;
                let mut csi_query = noodles::csi::io::Query::new(reader.get_mut(), chunks);
                let mut record = noodles::bam::Record::default();
                let mut cigar_scratch = Vec::new();
                let mut finished = false;
                let iter = std::iter::from_fn(move || {
                    if finished {
                        return None;
                    }
                    loop {
                        match read_bounded_bam_record(&mut csi_query, &mut record) {
                            Ok(0) => {
                                finished = true;
                                return None;
                            }
                            Ok(_) => {
                                match intersects_record(
                                    &record,
                                    reference_sequence_id,
                                    interval,
                                    &mut cigar_scratch,
                                ) {
                                    Ok(Some(_)) => {
                                        return Some(SamRecord::from_alignment_record(
                                            &record, header,
                                        ));
                                    }
                                    Ok(None) => continue,
                                    Err(e) => {
                                        finished = true;
                                        return Some(Err(e));
                                    }
                                }
                            }
                            Err(e) => {
                                finished = true;
                                return Some(Err(e));
                            }
                        }
                    }
                });
                Ok(Box::new(iter))
            }
            IndexedSamReader::Cram(path, meta) => {
                let reference_sequence_id = header
                    .reference_sequences()
                    .get_index_of(region.name())
                    .or_else(|| {
                        meta.header
                            .reference_sequences()
                            .get_index_of(region.name())
                    })
                    .ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidInput,
                            "invalid reference sequence name",
                        )
                    })?;
                let interval = region.interval();
                let container_offsets = meta
                    .indexed_crai
                    .query_offsets(reference_sequence_id, interval);
                let hint_contig = region.name().to_vec();
                let pinned_slot = if container_offsets.is_empty() {
                    None
                } else {
                    Some(meta.fasta_repo.warm_contig_slot(&hint_contig))
                };
                let feed = CramFeed {
                    inner: CramFeedInner {
                        path: path.clone(),
                        meta: std::sync::Arc::clone(meta),
                        reference_sequence_id,
                        hint_contig,
                        pinned_slot,
                        interval,
                        container_offsets: container_offsets.into_iter(),
                        current_block: None,
                    },
                    sam: SamRecord::new(),
                };
                let iter = feed.map(|res| {
                    res.map_err(|e| match e {
                        crate::types::VardictError::Io(io_err) => io_err,
                        other => {
                            std::io::Error::new(std::io::ErrorKind::InvalidData, other.to_string())
                        }
                    })
                });
                Ok(Box::new(iter))
            }
        }
    }
}

#[deprecated(since = "1.0.10", note = "use bam::query_region_iter")]
#[allow(deprecated)]
pub struct CachedReader {
    pub reader: IndexedSamReader,
    pub header: noodles::sam::Header,
}

fn fasta_fai_path(fasta_path: &Path) -> PathBuf {
    crate::reference::resolve_fasta_fai_path(fasta_path)
}

#[inline]
fn is_cram_path(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("cram"))
}

#[inline]
fn is_absent_index_candidate(cand: &Path, err: &std::io::Error) -> bool {
    err.kind() == std::io::ErrorKind::NotFound
        && matches!(
            std::fs::symlink_metadata(cand),
            Err(meta_err) if meta_err.kind() == std::io::ErrorKind::NotFound
        )
}

const MAX_BAI_BYTES: u64 = 512 * 1024 * 1024;
const MAX_CSI_BYTES: u64 = 256 * 1024 * 1024;
const MAX_CSI_DECOMPRESSED_BYTES: u64 = 256 * 1024 * 1024;
const BAI_DEPTH: u8 = 5;
const BAI_MAX_BIN_ID: u32 =
    match noodles::csi::binning_index::index::reference_sequence::Bin::max_id(BAI_DEPTH) {
        Some(id) => id as u32,
        None => 0,
    };
const BAI_METADATA_BIN_ID: u32 =
    match noodles::csi::binning_index::index::reference_sequence::Bin::metadata_id(BAI_DEPTH) {
        Some(id) => id as u32,
        None => 0,
    };
const MAX_BAM_INDEX_REFS: usize = 16_777_216;
const MAX_CSI_TOTAL_BINS: usize = 4_194_304;
const MAX_BAM_HEADER_BYTES: usize = 128 * 1024 * 1024;
const MAX_BAM_REF_NAME_BYTES: usize = 65_536;
const MAX_BAM_RECORD_BYTES: usize = 64 * 1024 * 1024;
const BAM_RECORD_FAST_PATH_BYTES: usize = 65_536;
const MAX_CRAM_HEADER_BYTES: usize = 128 * 1024 * 1024;

fn read_validated_bam_header<R: std::io::Read>(
    r: &mut R,
    path: &Path,
) -> std::io::Result<noodles::sam::Header> {
    use std::io::Read as _;

    let mut magic = [0_u8; 4];
    r.read_exact(&mut magic)?;
    if magic != *b"BAM\x01" {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("invalid BAM magic header in '{}'", path.display()),
        ));
    }

    let mut l_text_buf = [0_u8; 4];
    r.read_exact(&mut l_text_buf)?;
    let raw_l_text = u32::from_le_bytes(l_text_buf);
    let l_text = usize::try_from(raw_l_text)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    if l_text > MAX_BAM_HEADER_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "BAM header l_text ({raw_l_text}) in '{}' exceeds maximum allowed ({MAX_BAM_HEADER_BYTES} bytes)",
                path.display()
            ),
        ));
    }

    let mut budget_used = l_text;
    let mut raw = Vec::with_capacity(12usize.saturating_add(l_text.min(65_536)));
    raw.extend_from_slice(&magic);
    raw.extend_from_slice(&l_text_buf);
    let text_read = r.take(l_text as u64).read_to_end(&mut raw)?;
    if text_read < l_text {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            format!(
                "truncated BAM header text in '{}': expected {l_text} bytes, got {text_read}",
                path.display()
            ),
        ));
    }

    let mut n_ref_buf = [0_u8; 4];
    r.read_exact(&mut n_ref_buf)?;
    let raw_n_ref = u32::from_le_bytes(n_ref_buf);
    let n_ref = usize::try_from(raw_n_ref)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    if n_ref > MAX_BAM_INDEX_REFS {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "BAM header n_ref ({raw_n_ref}) in '{}' exceeds maximum allowed ({MAX_BAM_INDEX_REFS})",
                path.display()
            ),
        ));
    }
    raw.extend_from_slice(&n_ref_buf);

    for _ in 0..n_ref {
        let mut l_name_buf = [0_u8; 4];
        r.read_exact(&mut l_name_buf)?;
        let raw_l_name = u32::from_le_bytes(l_name_buf);
        let l_name = usize::try_from(raw_l_name)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        if l_name == 0 || l_name > MAX_BAM_REF_NAME_BYTES {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "BAM reference l_name ({raw_l_name}) in '{}' outside allowed range 1..={MAX_BAM_REF_NAME_BYTES}",
                    path.display()
                ),
            ));
        }
        let entry_bytes = 8usize.saturating_add(l_name);
        budget_used = budget_used
            .checked_add(entry_bytes)
            .filter(|&b| b <= MAX_BAM_HEADER_BYTES)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "BAM header size in '{}' exceeds maximum allowed ({MAX_BAM_HEADER_BYTES} bytes)",
                        path.display()
                    ),
                )
            })?;
        raw.extend_from_slice(&l_name_buf);
        raw.reserve(l_name.saturating_add(4));
        let name_read = r.take(l_name as u64).read_to_end(&mut raw)?;
        if name_read < l_name {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!("truncated BAM reference name in '{}'", path.display()),
            ));
        }
        if raw.last().copied() != Some(0) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "BAM reference name in '{}' missing NUL terminator",
                    path.display()
                ),
            ));
        }

        let mut l_ref_buf = [0_u8; 4];
        r.read_exact(&mut l_ref_buf)?;
        let raw_l_ref = i32::from_le_bytes(l_ref_buf);
        if raw_l_ref < 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "negative BAM reference length ({raw_l_ref}) in '{}'",
                    path.display()
                ),
            ));
        }
        raw.extend_from_slice(&l_ref_buf);
    }

    noodles::bam::io::Reader::from(raw.as_slice()).read_header()
}

fn read_bounded_bam_record<R: std::io::Read>(
    r: &mut R,
    record: &mut noodles::bam::Record,
) -> std::io::Result<usize> {
    use std::io::Read as _;

    let mut size_buf = [0_u8; 4];
    let mut prefix_read = 0_usize;
    while let Some(rem) = size_buf.get_mut(prefix_read..) {
        if rem.is_empty() {
            break;
        }
        match r.read(rem) {
            Ok(0) => break,
            Ok(n) => prefix_read = prefix_read.saturating_add(n),
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    if prefix_read == 0 {
        return Ok(0);
    }
    if prefix_read < 4 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "truncated BAM record block_size prefix",
        ));
    }

    let raw_block_size = u32::from_le_bytes(size_buf);
    let block_size = usize::try_from(raw_block_size)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    if block_size > MAX_BAM_RECORD_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "BAM record block_size ({raw_block_size}) exceeds maximum allowed ({MAX_BAM_RECORD_BYTES} bytes)"
            ),
        ));
    }

    if block_size <= BAM_RECORD_FAST_PATH_BYTES {
        let chained = std::io::Read::chain(&size_buf[..], r);
        noodles::bam::io::Reader::from(chained).read_record(record)
    } else {
        let mut raw = Vec::with_capacity(4 + BAM_RECORD_FAST_PATH_BYTES);
        raw.extend_from_slice(&size_buf);
        let payload_read = r.take(block_size as u64).read_to_end(&mut raw)?;
        if payload_read < block_size {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!("truncated BAM record: expected {block_size} bytes, got {payload_read}"),
            ));
        }
        noodles::bam::io::Reader::from(raw.as_slice()).read_record(record)
    }
}

fn take_bai_rem(rem: &mut u64, n: u64, path: &Path) -> std::io::Result<()> {
    *rem = rem.checked_sub(n).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            format!("truncated BAI index '{}'", path.display()),
        )
    })?;
    Ok(())
}

fn read_bai_u32_le<R: std::io::Read>(
    reader: &mut R,
    rem: &mut u64,
    path: &Path,
) -> std::io::Result<u32> {
    take_bai_rem(rem, 4, path)?;
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_bai_i32_le<R: std::io::Read>(
    reader: &mut R,
    rem: &mut u64,
    path: &Path,
) -> std::io::Result<i32> {
    take_bai_rem(rem, 4, path)?;
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_validated_bai_index(path: &Path) -> std::io::Result<noodles::bam::bai::Index> {
    use std::io::{Read as _, Seek as _, SeekFrom};

    let file = std::fs::File::open(path)?;
    let file_len = file.metadata()?.len();
    if file_len > MAX_BAI_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "BAI index '{}' ({file_len} bytes) exceeds maximum size of {MAX_BAI_BYTES} bytes",
                path.display()
            ),
        ));
    }

    let mut reader = std::io::BufReader::new(file);
    let mut rem = file_len;

    take_bai_rem(&mut rem, 4, path)?;
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if magic != *b"BAI\x01" {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("invalid BAI header in '{}'", path.display()),
        ));
    }

    let n_ref = read_bai_u32_le(&mut reader, &mut rem, path)?;
    if (n_ref as usize) > MAX_BAM_INDEX_REFS {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "BAI n_ref ({n_ref}) in '{}' exceeds maximum allowed ({MAX_BAM_INDEX_REFS})",
                path.display()
            ),
        ));
    }
    let min_refs_bytes = u64::from(n_ref).checked_mul(8).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("BAI n_ref overflow in '{}'", path.display()),
        )
    })?;
    if min_refs_bytes > rem {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "corrupt BAI index '{}': n_ref ({n_ref}) exceeds remaining {rem} bytes",
                path.display()
            ),
        ));
    }

    for _ in 0..n_ref {
        let n_bin = read_bai_u32_le(&mut reader, &mut rem, path)?;
        if n_bin > BAI_METADATA_BIN_ID {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "BAI n_bin ({n_bin}) in '{}' exceeds maximum allowed ({BAI_METADATA_BIN_ID})",
                    path.display()
                ),
            ));
        }
        let min_bins_bytes = u64::from(n_bin).checked_mul(8).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("BAI n_bin overflow in '{}'", path.display()),
            )
        })?;
        if min_bins_bytes > rem {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "corrupt BAI index '{}': n_bin ({n_bin}) exceeds remaining {rem} bytes",
                    path.display()
                ),
            ));
        }

        let mut seen_bins = [0u64; (BAI_METADATA_BIN_ID as usize / 64) + 1];
        for _ in 0..n_bin {
            let bin_id = read_bai_u32_le(&mut reader, &mut rem, path)?;
            if bin_id > BAI_MAX_BIN_ID && bin_id != BAI_METADATA_BIN_ID {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "invalid BAI bin_id ({bin_id}) in '{}' (max={BAI_MAX_BIN_ID}, metadata={BAI_METADATA_BIN_ID})",
                        path.display()
                    ),
                ));
            }
            let bin_idx = bin_id as usize;
            let mask = 1u64 << (bin_idx % 64);
            let Some(word) = seen_bins.get_mut(bin_idx / 64) else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("BAI bin_id ({bin_id}) out of range in '{}'", path.display()),
                ));
            };
            if (*word & mask) != 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("duplicate BAI bin_id ({bin_id}) in '{}'", path.display()),
                ));
            }
            *word |= mask;
            if bin_id == BAI_METADATA_BIN_ID {
                if rem < 36 {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("truncated BAI metadata bin in '{}'", path.display()),
                    ));
                }
                take_bai_rem(&mut rem, 36, path)?;
                reader.seek_relative(36)?;
            } else {
                let n_chunk = read_bai_i32_le(&mut reader, &mut rem, path)?;
                let n_chunk_u64 = u64::try_from(n_chunk).map_err(|_| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "negative BAI chunk count ({n_chunk}) in '{}'",
                            path.display()
                        ),
                    )
                })?;
                let skip = n_chunk_u64.checked_mul(16).ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("BAI chunk byte count overflow in '{}'", path.display()),
                    )
                })?;
                if skip > rem {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "corrupt BAI index '{}': n_chunk ({n_chunk}) exceeds remaining {rem} bytes",
                            path.display()
                        ),
                    ));
                }
                take_bai_rem(&mut rem, skip, path)?;
                let skip_i64 = i64::try_from(skip).map_err(|_| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("BAI chunk skip overflow in '{}'", path.display()),
                    )
                })?;
                reader.seek_relative(skip_i64)?;
            }
        }

        let n_intv = read_bai_u32_le(&mut reader, &mut rem, path)?;
        let skip = u64::from(n_intv).checked_mul(8).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("BAI interval byte count overflow in '{}'", path.display()),
            )
        })?;
        if skip > rem {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "corrupt BAI index '{}': n_intv ({n_intv}) exceeds remaining {rem} bytes",
                    path.display()
                ),
            ));
        }
        take_bai_rem(&mut rem, skip, path)?;
        let skip_i64 = i64::try_from(skip).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("BAI interval skip overflow in '{}'", path.display()),
            )
        })?;
        reader.seek_relative(skip_i64)?;
    }

    reader.seek(SeekFrom::Start(0))?;
    noodles::bam::bai::io::Reader::new(reader).read_index()
}

fn take_csi_slice<'a>(slice: &mut &'a [u8], len: usize, path: &Path) -> std::io::Result<&'a [u8]> {
    let Some((head, tail)) = slice.split_at_checked(len) else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("truncated CSI index '{}'", path.display()),
        ));
    };
    *slice = tail;
    Ok(head)
}

fn read_csi_u32_le(slice: &mut &[u8], path: &Path) -> std::io::Result<u32> {
    let bytes = take_csi_slice(slice, 4, path)?;
    let arr: [u8; 4] = bytes.try_into().map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("truncated CSI index '{}'", path.display()),
        )
    })?;
    Ok(u32::from_le_bytes(arr))
}

fn read_csi_i32_le(slice: &mut &[u8], path: &Path) -> std::io::Result<i32> {
    let bytes = take_csi_slice(slice, 4, path)?;
    let arr: [u8; 4] = bytes.try_into().map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("truncated CSI index '{}'", path.display()),
        )
    })?;
    Ok(i32::from_le_bytes(arr))
}

fn validate_csi_payload(mut slice: &[u8], path: &Path) -> std::io::Result<()> {
    use noodles::csi::binning_index::index::reference_sequence::Bin;

    let magic = take_csi_slice(&mut slice, 4, path)?;
    if magic != b"CSI\x01" {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("invalid CSI magic header in '{}'", path.display()),
        ));
    }

    let min_shift = read_csi_i32_le(&mut slice, path)?;
    let Ok(min_shift_u8) = u8::try_from(min_shift) else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "invalid CSI min_shift ({min_shift}) in '{}'",
                path.display()
            ),
        ));
    };

    let depth = read_csi_i32_le(&mut slice, path)?;
    let depth_u8 = u8::try_from(depth)
        .ok()
        .filter(|&d| d <= 9)
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("invalid CSI depth ({depth}) in '{}'", path.display()),
            )
        })?;
    if u32::from(min_shift_u8) + u32::from(depth_u8) * 3 >= usize::BITS {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "CSI coordinate shift ({min_shift} + 3 * {depth}) overflows usize::BITS in '{}'",
                path.display()
            ),
        ));
    }
    let (Some(max_bin_id), Some(metadata_id)) = (Bin::max_id(depth_u8), Bin::metadata_id(depth_u8))
    else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("invalid CSI depth ({depth}) in '{}'", path.display()),
        ));
    };

    let l_aux = read_csi_i32_le(&mut slice, path)?;
    let l_aux_usize = usize::try_from(l_aux).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("negative CSI l_aux ({l_aux}) in '{}'", path.display()),
        )
    })?;
    let aux_bytes = take_csi_slice(&mut slice, l_aux_usize, path)?;
    if l_aux_usize > 0 {
        let Some(l_nm_raw) = aux_bytes.get(24..28) else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("truncated CSI aux header in '{}'", path.display()),
            ));
        };
        let l_nm_arr: [u8; 4] = l_nm_raw.try_into().map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("truncated CSI aux header in '{}'", path.display()),
            )
        })?;
        let l_nm = usize::try_from(i32::from_le_bytes(l_nm_arr)).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("negative CSI l_nm in '{}'", path.display()),
            )
        })?;
        if 28usize.checked_add(l_nm) != Some(l_aux_usize) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "corrupt CSI aux header length ({l_aux_usize} != 28 + {l_nm}) in '{}'",
                    path.display()
                ),
            ));
        }
    }

    let n_ref = read_csi_i32_le(&mut slice, path)?;
    let n_ref_usize = usize::try_from(n_ref).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("negative CSI n_ref ({n_ref}) in '{}'", path.display()),
        )
    })?;
    if n_ref_usize > MAX_BAM_INDEX_REFS {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "CSI n_ref ({n_ref}) in '{}' exceeds maximum allowed ({MAX_BAM_INDEX_REFS})",
                path.display()
            ),
        ));
    }
    if n_ref_usize.checked_mul(4).is_none_or(|b| b > slice.len()) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "corrupt CSI index '{}': n_ref ({n_ref}) exceeds remaining {} bytes",
                path.display(),
                slice.len()
            ),
        ));
    }

    let mut total_bins = 0usize;
    for _ in 0..n_ref_usize {
        let n_bin = read_csi_i32_le(&mut slice, path)?;
        let n_bin_usize = usize::try_from(n_bin).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("negative CSI n_bin ({n_bin}) in '{}'", path.display()),
            )
        })?;
        if n_bin_usize > metadata_id {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CSI n_bin ({n_bin}) in '{}' exceeds maximum allowed ({metadata_id})",
                    path.display()
                ),
            ));
        }
        if n_bin_usize.checked_mul(16).is_none_or(|b| b > slice.len()) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "corrupt CSI index '{}': n_bin ({n_bin}) exceeds remaining {} bytes",
                    path.display(),
                    slice.len()
                ),
            ));
        }
        total_bins = total_bins
            .checked_add(n_bin_usize)
            .filter(|&tot| tot <= MAX_CSI_TOTAL_BINS)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CSI total bin count in '{}' exceeds maximum allowed ({MAX_CSI_TOTAL_BINS})",
                        path.display()
                    ),
                )
            })?;

        let mut seen_bins = rustc_hash::FxHashSet::with_capacity_and_hasher(
            n_bin_usize.min(4_096),
            Default::default(),
        );
        for _ in 0..n_bin_usize {
            let bin_id = read_csi_u32_le(&mut slice, path)?;
            let bin_id_usize = usize::try_from(bin_id).map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("invalid CSI bin_id ({bin_id}) in '{}'", path.display()),
                )
            })?;
            if bin_id_usize != metadata_id && bin_id_usize > max_bin_id {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "invalid CSI bin_id ({bin_id}) in '{}' (max={max_bin_id}, metadata={metadata_id})",
                        path.display()
                    ),
                ));
            }
            if !seen_bins.insert(bin_id_usize) {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("duplicate CSI bin_id ({bin_id}) in '{}'", path.display()),
                ));
            }
            let _loffset = take_csi_slice(&mut slice, 8, path)?;
            if bin_id_usize == metadata_id {
                let _metadata = take_csi_slice(&mut slice, 36, path)?;
            } else {
                let n_chunk = read_csi_i32_le(&mut slice, path)?;
                let n_chunk_usize = usize::try_from(n_chunk).map_err(|_| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "negative CSI chunk count ({n_chunk}) in '{}'",
                            path.display()
                        ),
                    )
                })?;
                let chunk_bytes = n_chunk_usize.checked_mul(16).ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("CSI chunk byte count overflow in '{}'", path.display()),
                    )
                })?;
                let _chunks = take_csi_slice(&mut slice, chunk_bytes, path)?;
            }
        }
    }

    Ok(())
}

fn read_validated_csi_index(path: &Path) -> std::io::Result<noodles::csi::Index> {
    use std::io::{Read as _, Seek as _, SeekFrom};

    let mut file = std::fs::File::open(path)?;
    let file_len = file.metadata()?.len();
    if file_len > MAX_CSI_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "CSI index '{}' ({file_len} bytes) exceeds maximum size of {MAX_CSI_BYTES} bytes",
                path.display()
            ),
        ));
    }

    let mut decomp = Vec::new();
    noodles::bgzf::io::Reader::new(&mut file)
        .take(MAX_CSI_DECOMPRESSED_BYTES + 1)
        .read_to_end(&mut decomp)?;
    if (decomp.len() as u64) > MAX_CSI_DECOMPRESSED_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "CSI index '{}' exceeds maximum decompressed size of {MAX_CSI_DECOMPRESSED_BYTES} bytes",
                path.display()
            ),
        ));
    }

    validate_csi_payload(&decomp, path)?;

    file.seek(SeekFrom::Start(0))?;
    noodles::csi::io::Reader::new(file).read_index()
}

fn try_bam_index_for_path(src: &Path) -> Result<noodles::bam::Index, (std::io::Error, bool)> {
    let primary_bai = crate::reference::append_os_ext(src, ".bai");
    let first_not_found = match read_validated_bai_index(&primary_bai) {
        Ok(idx) => return Ok(noodles::bam::Index::Bai(idx)),
        Err(e) if is_absent_index_candidate(&primary_bai, &e) => e,
        Err(e) => return Err((e, false)),
    };

    let primary_csi = crate::reference::append_os_ext(src, ".csi");
    match read_validated_csi_index(&primary_csi) {
        Ok(idx) => return Ok(noodles::bam::Index::Csi(idx)),
        Err(e) if is_absent_index_candidate(&primary_csi, &e) => {}
        Err(e) => return Err((e, false)),
    }

    if src.extension().is_some() {
        let stem_bai = src.with_extension("bai");
        if stem_bai != primary_bai && stem_bai != src {
            match read_validated_bai_index(&stem_bai) {
                Ok(idx) => return Ok(noodles::bam::Index::Bai(idx)),
                Err(e) if is_absent_index_candidate(&stem_bai, &e) => {}
                Err(e) => return Err((e, false)),
            }
        }
        let stem_csi = src.with_extension("csi");
        if stem_csi != primary_csi && stem_csi != src {
            match read_validated_csi_index(&stem_csi) {
                Ok(idx) => return Ok(noodles::bam::Index::Csi(idx)),
                Err(e) if is_absent_index_candidate(&stem_csi, &e) => {}
                Err(e) => return Err((e, false)),
            }
        }
    }

    Err((first_not_found, true))
}

fn read_bam_index_with_fallbacks(path: &Path) -> std::io::Result<noodles::bam::Index> {
    let first_not_found = match try_bam_index_for_path(path) {
        Ok(idx) => return Ok(idx),
        Err((e, false)) => return Err(e),
        Err((e, true)) => e,
    };

    if let Ok(canon) = std::fs::canonicalize(path) {
        if canon != path {
            match try_bam_index_for_path(&canon) {
                Ok(idx) => return Ok(idx),
                Err((e, false)) => return Err(e),
                Err((_, true)) => {}
            }
        }
    }

    Err(first_not_found)
}

const MAX_CRAI_BYTES: u64 = 128 * 1024 * 1024;
const MAX_CRAI_LINE_BYTES: usize = 64 * 1024;
const MAX_CRAM_CONTAINER_BYTES: usize = 256 * 1024 * 1024;
const MAX_CRAM_CONTAINER_RECORDS: usize = 2_000_000;
/// Maximum record count (`65,536`) of a decoded CRAM container eligible for [`CRAM_CONTAINER_CACHE`].
const MAX_CACHED_CRAM_CONTAINER_RECORDS: usize = 65_536;
/// Maximum decoded heap footprint (`16 MiB`) cached per [`CRAM_CONTAINER_CACHE`] container.
///
/// Standard `samtools`/`htslib` (`10,000`-read) and `noodles-cram` (`10,240`-read) `150–250 bp`
/// Illumina containers occupy `~6.3–9.8 MiB` (`HG008` max `6.63 MiB`); values `< 8 MiB` reject
/// standard Illumina containers and breach `perf_budget.py`'s `<= 2.50×` CRAM/BAM wall ratio.
/// Across `4` shards × `4` slots (`16` slots) with a `40 MiB` per-shard byte budget
/// ([`MAX_SHARD_CACHED_CRAM_BYTES`]), decoded-record cache residency stays `<= 160 MiB`
/// (separate from `SharedFastaRepo`'s `512 MiB` CRAM reference cache) while bypassing oversized
/// PacBio/ONT long-read containers.
const MAX_CACHED_CRAM_CONTAINER_BYTES: usize = 16 * 1024 * 1024;

fn parse_crai_line(s: &str) -> std::io::Result<noodles::cram::crai::Record> {
    const UNMAPPED: i32 = -1;
    let mut fields = s.splitn(6, '\t');
    let mut next_field = || -> std::io::Result<&str> {
        fields
            .next()
            .ok_or_else(|| std::io::Error::from(std::io::ErrorKind::UnexpectedEof))
    };
    let ref_seq_raw: i32 = next_field()?
        .parse()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let reference_sequence_id = match ref_seq_raw {
        UNMAPPED => None,
        n => Some(
            usize::try_from(n)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?,
        ),
    };
    let alignment_start_raw: usize = next_field()?
        .parse()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let alignment_start = noodles::core::Position::new(alignment_start_raw);
    let alignment_span: usize = next_field()?
        .parse()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let offset: u64 = next_field()?
        .parse()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let landmark: u64 = next_field()?
        .parse()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let slice_length: u64 = next_field()?
        .parse()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    Ok(noodles::cram::crai::Record::new(
        reference_sequence_id,
        alignment_start,
        alignment_span,
        offset,
        landmark,
        slice_length,
    ))
}

fn read_crai_entries(crai_path: &Path) -> std::io::Result<noodles::cram::crai::Index> {
    use std::io::{BufRead as _, Read as _};

    let file = std::fs::File::open(crai_path)?;
    let gz = flate2::read::GzDecoder::new(file).take(MAX_CRAI_BYTES + 1);
    let mut reader = std::io::BufReader::new(gz);
    let mut total_bytes = 0_u64;
    let mut line_buf: Vec<u8> = Vec::with_capacity(128);
    let mut records = Vec::new();

    loop {
        line_buf.clear();
        loop {
            let available = reader.fill_buf()?;
            if available.is_empty() {
                break;
            }
            if let Some(nl_pos) = available.iter().position(|&b| b == b'\n') {
                let take_len = nl_pos + 1;
                total_bytes = total_bytes.saturating_add(take_len as u64);
                if total_bytes > MAX_CRAI_BYTES {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM index '{}' exceeds maximum decompressed size of {} bytes",
                            crai_path.display(),
                            MAX_CRAI_BYTES
                        ),
                    ));
                }
                if line_buf.len().saturating_add(take_len) > MAX_CRAI_LINE_BYTES {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM index '{}' contains an oversized line (> {} bytes)",
                            crai_path.display(),
                            MAX_CRAI_LINE_BYTES
                        ),
                    ));
                }
                if let Some(prefix) = available.get(..take_len) {
                    line_buf.extend_from_slice(prefix);
                }
                reader.consume(take_len);
                break;
            } else {
                let take_len = available.len();
                total_bytes = total_bytes.saturating_add(take_len as u64);
                if total_bytes > MAX_CRAI_BYTES {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM index '{}' exceeds maximum decompressed size of {} bytes",
                            crai_path.display(),
                            MAX_CRAI_BYTES
                        ),
                    ));
                }
                if line_buf.len().saturating_add(take_len) > MAX_CRAI_LINE_BYTES {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM index '{}' contains an oversized line (> {} bytes)",
                            crai_path.display(),
                            MAX_CRAI_LINE_BYTES
                        ),
                    ));
                }
                line_buf.extend_from_slice(available);
                reader.consume(take_len);
            }
        }

        if line_buf.is_empty() {
            break;
        }
        if line_buf.last() == Some(&b'\n') {
            line_buf.pop();
            if line_buf.last() == Some(&b'\r') {
                line_buf.pop();
            }
        }
        let line_str = std::str::from_utf8(&line_buf)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        records.push(parse_crai_line(line_str)?);
    }

    Ok(records)
}

fn try_cram_index_for_path(
    src: &Path,
) -> Result<noodles::cram::crai::Index, (std::io::Error, bool)> {
    let primary_crai = crate::reference::append_os_ext(src, ".crai");
    let first_not_found = match read_crai_entries(&primary_crai) {
        Ok(idx) => return Ok(idx),
        Err(e) if is_absent_index_candidate(&primary_crai, &e) => e,
        Err(e) => return Err((e, false)),
    };

    if src.extension().is_some() {
        let stem_crai = src.with_extension("crai");
        if stem_crai != primary_crai && stem_crai != src {
            match read_crai_entries(&stem_crai) {
                Ok(idx) => return Ok(idx),
                Err(e) if is_absent_index_candidate(&stem_crai, &e) => {}
                Err(e) => return Err((e, false)),
            }
        }
    }

    Err((first_not_found, true))
}

fn read_cram_index_with_fallbacks(path: &Path) -> std::io::Result<noodles::cram::crai::Index> {
    let first_not_found = match try_cram_index_for_path(path) {
        Ok(idx) => return Ok(idx),
        Err((e, false)) => return Err(e),
        Err((e, true)) => e,
    };

    if let Ok(canon) = std::fs::canonicalize(path) {
        if canon != path {
            match try_cram_index_for_path(&canon) {
                Ok(idx) => return Ok(idx),
                Err((e, false)) => return Err(e),
                Err((_, true)) => {}
            }
        }
    }

    Err(first_not_found)
}

#[deprecated(since = "1.0.10", note = "use bam::query_region_iter")]
#[allow(deprecated)]
pub fn open_reader(
    path: &Path,
    config: &crate::config::Config,
) -> Result<CachedReader, crate::types::VardictError> {
    if is_cram_path(path) {
        let meta = acquire_cram_meta(path, config)?;
        let header = meta.header.clone();
        Ok(CachedReader {
            reader: IndexedSamReader::Cram(path.to_path_buf(), meta),
            header,
        })
    } else {
        let index = read_bam_index_with_fallbacks(path).map_err(crate::types::VardictError::Io)?;
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_path(path)
            .map_err(crate::types::VardictError::Io)?;
        let header = read_validated_bam_header(reader.get_mut(), path)
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
    query_region_iter(bam_path, chr, start, end, config)?.collect()
}

static NEXT_BAM_META_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
static BAM_CACHE_EPOCH: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
/// Synthetic `meta_id` floor reserved for isolated unit tests so concurrent
/// `clear_bam_caches()` calls in parallel tests never evict in-flight test blocks (`AR-4`).
#[cfg(test)]
const TEST_META_ID_BASE: u64 = 1u64 << 63;
#[cfg(test)]
static NEXT_PROD_TEST_META_ID: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(1u64 << 32);
#[cfg(test)]
type OneShotTestHook = Box<dyn FnOnce()>;

#[cfg(test)]
thread_local! {
    static FORCE_PRODUCTION_META_ID: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static TEST_MID_CLEAR_HOOK: RefCell<Option<OneShotTestHook>> = RefCell::new(None);
    static TEST_MID_COLD_OPEN_HOOK: RefCell<Option<OneShotTestHook>> = RefCell::new(None);
    static TEST_MID_CRAM_DECODE_HOOK: RefCell<Option<OneShotTestHook>> = RefCell::new(None);
}

#[inline]
fn next_bam_meta_id() -> u64 {
    #[cfg(test)]
    if FORCE_PRODUCTION_META_ID.with(|c| c.get()) {
        return NEXT_PROD_TEST_META_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
    NEXT_BAM_META_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}

#[inline]
fn tls_epoch_valid(cached_epoch: u64, cur_epoch: u64, meta_id: u64) -> bool {
    #[cfg(test)]
    if meta_id >= TEST_META_ID_BASE {
        return true;
    }
    #[cfg(not(test))]
    let _ = meta_id;
    cached_epoch == cur_epoch
}

struct SharedBamMeta {
    id: u64,
    path: PathBuf,
    header: noodles::sam::Header,
    contig_names: Vec<Box<str>>,
    index: noodles::bam::Index,
}

/// Process-wide BGZF block cache sharding (`32 shards × 16 blocks = 512 blocks` <= 32 MiB total).
const BGZF_CACHE_SHARDS: usize = 32;
const BGZF_CACHE_SLOTS_PER_SHARD: usize = 16;
/// 64-bit golden-ratio fractional constant (`2^64 / phi`, splitmix64) for offset hashing.
const GOLDEN_RATIO_64: u64 = 0x9E37_79B9_7F4A_7C15;
/// 64-bit PCG multiplier for mixing `meta_id` into cache shard indices.
const PCG_MULT_64: u64 = 0x517C_C1B7_2722_0A95;

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

struct LruCacheShard<V, const SLOTS: usize> {
    entries: Vec<((u64, u64), V)>,
}

impl<V, const SLOTS: usize> Default for LruCacheShard<V, SLOTS> {
    fn default() -> Self {
        Self {
            entries: Vec::new(),
        }
    }
}

impl<V: Clone, const SLOTS: usize> LruCacheShard<V, SLOTS> {
    fn get(&mut self, key: (u64, u64)) -> Option<V> {
        let idx = self.entries.iter().position(|(k, _)| *k == key)?;
        let entry = self.entries.remove(idx);
        let val = entry.1.clone();
        self.entries.push(entry);
        Some(val)
    }

    fn insert(&mut self, key: (u64, u64), val: V) {
        if let Some(idx) = self.entries.iter().position(|(k, _)| *k == key) {
            self.entries.remove(idx);
        } else if self.entries.len() >= SLOTS {
            self.entries.remove(0);
        }
        self.entries.push((key, val));
    }

    fn clear_production(&mut self) {
        #[cfg(test)]
        self.entries.retain(|((id, _), _)| *id >= TEST_META_ID_BASE);
        #[cfg(not(test))]
        self.entries.clear();
    }
}

type BgzfCacheShard = LruCacheShard<CachedBgzfBlock, BGZF_CACHE_SLOTS_PER_SHARD>;

static BGZF_BLOCK_CACHE: std::sync::OnceLock<
    [std::sync::Mutex<BgzfCacheShard>; BGZF_CACHE_SHARDS],
> = std::sync::OnceLock::new();

fn bgzf_cache_shards() -> &'static [std::sync::Mutex<BgzfCacheShard>; BGZF_CACHE_SHARDS] {
    BGZF_BLOCK_CACHE
        .get_or_init(|| std::array::from_fn(|_| std::sync::Mutex::new(BgzfCacheShard::default())))
}

#[inline]
fn cache_shard_index<const SHARDS: usize>(meta_id: u64, offset: u64) -> usize {
    const { assert!(SHARDS.is_power_of_two()) };
    let h = offset
        .wrapping_mul(GOLDEN_RATIO_64)
        .wrapping_add(meta_id.wrapping_mul(PCG_MULT_64));
    ((h ^ (h >> 32)) as usize) & (SHARDS - 1)
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

impl std::fmt::Debug for CachedBgzfReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CachedBgzfReader")
            .field("meta_id", &self.meta_id)
            .field("upos", &self.upos)
            .field("at_eof", &self.at_eof)
            .field("inner_next_cpos", &self.inner_next_cpos)
            .finish_non_exhaustive()
    }
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

        let shard_idx = cache_shard_index::<BGZF_CACHE_SHARDS>(self.meta_id, cpos);
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
            let actual_shard = cache_shard_index::<BGZF_CACHE_SHARDS>(self.meta_id, block_cpos);
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
type BamMetaTlsSlots = [Option<(u64, std::sync::Arc<SharedBamMeta>)>; 2];

static META_CACHE: std::sync::OnceLock<std::sync::Mutex<HashMap<PathBuf, PathSlot>>> =
    std::sync::OnceLock::new();

const MAX_POOLED_READERS_PER_PATH: usize = 4;

thread_local! {
    static BAM_READER_POOL: RefCell<BamReaderPoolMap> = RefCell::new(HashMap::new());
    static BAM_META_TLS: RefCell<BamMetaTlsSlots> = const { RefCell::new([None, None]) };
}

/// Clear the process-wide shared BAM/CRAM metadata caches, bounded FASTA
/// repository, decoded block/container caches, and the calling thread's
/// pooled readers. Called at the start and end of [`crate::run`] via an
/// RAII drop guard (`F-5`).
pub fn clear_bam_caches() {
    BAM_CACHE_EPOCH.fetch_add(1, std::sync::atomic::Ordering::Release);
    #[cfg(test)]
    if let Some(hook) = TEST_MID_CLEAR_HOOK
        .try_with(|h| h.try_borrow_mut().ok().and_then(|mut g| g.take()))
        .ok()
        .flatten()
    {
        hook();
    }
    if let Some(cache) = META_CACHE.get() {
        let mut map = cache.lock().unwrap_or_else(|e| e.into_inner());
        #[cfg(test)]
        map.retain(|_, slot| {
            slot.lock()
                .unwrap_or_else(|e| e.into_inner())
                .as_ref()
                .is_some_and(|m| m.id >= TEST_META_ID_BASE)
        });
        #[cfg(not(test))]
        map.clear();
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
    let _ = BAM_META_TLS.try_with(|tls| {
        if let Ok(mut slots) = tls.try_borrow_mut() {
            *slots = [None, None];
        }
    });
    if let Some(cache) = CRAM_META_CACHE.get() {
        let mut map = cache.lock().unwrap_or_else(|e| e.into_inner());
        #[cfg(test)]
        map.retain(|_, slot| {
            slot.lock()
                .unwrap_or_else(|e| e.into_inner())
                .as_ref()
                .is_some_and(|m| m.id >= TEST_META_ID_BASE)
        });
        #[cfg(not(test))]
        map.clear();
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
    let _ = CRAM_META_TLS.try_with(|tls| {
        if let Ok(mut slots) = tls.try_borrow_mut() {
            *slots = [None, None];
        }
    });
    let _ = CRAM_LAST_CONTAINER_TLS.try_with(|tls| {
        if let Ok(mut slots) = tls.try_borrow_mut() {
            *slots = [None, None];
        }
    });
    let _ = CRAM_IO_SCRATCH_TLS.try_with(|tls| {
        if let Ok(mut buf) = tls.try_borrow_mut() {
            *buf = Vec::new();
        }
    });
    if let Some(cache) = FASTA_REPO_CACHE.get() {
        cache.lock().unwrap_or_else(|e| e.into_inner()).clear();
    }
    crate::reference::clear_fasta_reader_cache();
    BAM_CACHE_EPOCH.fetch_add(1, std::sync::atomic::Ordering::Release);
}

type ReaderPoolMap<M, R> = HashMap<PathBuf, Vec<(std::sync::Arc<M>, R)>>;

fn pop_from_pool<M: 'static, R: 'static>(
    pool: &'static std::thread::LocalKey<RefCell<ReaderPoolMap<M, R>>>,
    path: &Path,
    meta: &std::sync::Arc<M>,
) -> Option<R> {
    pool.try_with(|cell| {
        let mut map = cell.try_borrow_mut().ok()?;
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

fn return_to_pool<M: 'static, R: 'static>(
    pool: &'static std::thread::LocalKey<RefCell<ReaderPoolMap<M, R>>>,
    path: &Path,
    meta: std::sync::Arc<M>,
    reader: R,
) {
    let _ = pool.try_with(|cell| {
        if let Ok(mut map) = cell.try_borrow_mut() {
            if let Some(vec) = map.get_mut(path) {
                vec.retain(|(m, _)| std::sync::Arc::ptr_eq(m, &meta));
                if vec.len() < MAX_POOLED_READERS_PER_PATH {
                    vec.push((meta, reader));
                }
            } else if MAX_POOLED_READERS_PER_PATH > 0 {
                map.insert(path.to_path_buf(), vec![(meta, reader)]);
            }
        }
    });
}

fn pop_pooled_reader(path: &Path, meta: &std::sync::Arc<SharedBamMeta>) -> Option<RawBamReader> {
    pop_from_pool(&BAM_READER_POOL, path, meta)
}

fn return_pooled_reader(path: &Path, meta: std::sync::Arc<SharedBamMeta>, reader: RawBamReader) {
    return_to_pool(&BAM_READER_POOL, path, meta, reader);
}

#[inline]
fn lookup_bam_meta_tls(path: &Path, cur_epoch: u64) -> Option<std::sync::Arc<SharedBamMeta>> {
    BAM_META_TLS
        .try_with(|tls| {
            let mut slots = tls.try_borrow_mut().ok()?;
            let [s0, s1] = &mut *slots;
            if let Some((ep, meta)) = s0.as_ref() {
                if tls_epoch_valid(*ep, cur_epoch, meta.id) && meta.path.as_path() == path {
                    return Some(std::sync::Arc::clone(meta));
                }
            }
            if let Some((ep, meta)) = s1.as_ref() {
                if tls_epoch_valid(*ep, cur_epoch, meta.id) && meta.path.as_path() == path {
                    let hit = std::sync::Arc::clone(meta);
                    std::mem::swap(s0, s1);
                    return Some(hit);
                }
            }
            None
        })
        .ok()
        .flatten()
}

#[inline]
fn insert_bam_meta_tls(cur_epoch: u64, meta: &std::sync::Arc<SharedBamMeta>) {
    let _ = BAM_META_TLS.try_with(|tls| {
        if let Ok(mut slots) = tls.try_borrow_mut() {
            let [s0, s1] = &mut *slots;
            if s0.as_ref().is_some_and(|(_, m)| m.path == meta.path) {
                *s0 = Some((cur_epoch, std::sync::Arc::clone(meta)));
            } else {
                *s1 = s0.take();
                *s0 = Some((cur_epoch, std::sync::Arc::clone(meta)));
            }
        }
    });
}

fn acquire_bam_meta_and_reader(
    path: &Path,
) -> Result<(std::sync::Arc<SharedBamMeta>, RawBamReader), crate::types::VardictError> {
    let cur_epoch = BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire);
    if let Some(meta) = lookup_bam_meta_tls(path, cur_epoch) {
        if let Some(reader) = pop_pooled_reader(path, &meta) {
            return Ok((meta, reader));
        }
        let file = std::fs::File::open(path).map_err(crate::types::VardictError::Io)?;
        let reader = noodles::bam::io::Reader::from(CachedBgzfReader::new(file, meta.id));
        return Ok((meta, reader));
    }

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
            Some(meta) => {
                let meta_clone = std::sync::Arc::clone(meta);
                if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) == cur_epoch {
                    insert_bam_meta_tls(cur_epoch, &meta_clone);
                }
                meta_clone
            }
            None => {
                // Exact same open sequence as `open_reader` (index load ->
                // File::open -> read_header). Errors return without storing
                // `Some(...)` so transient/test errors are never cached and
                // io::ErrorKind is preserved.
                let index =
                    read_bam_index_with_fallbacks(path).map_err(crate::types::VardictError::Io)?;
                let file = std::fs::File::open(path).map_err(crate::types::VardictError::Io)?;
                let id = next_bam_meta_id();
                let mut reader = noodles::bam::io::Reader::from(CachedBgzfReader::new(file, id));
                let header = read_validated_bam_header(reader.get_mut(), path)
                    .map_err(crate::types::VardictError::Io)?;
                let contig_names: Vec<Box<str>> = header
                    .reference_sequences()
                    .iter()
                    .map(|(name, _)| Box::<str>::from(String::from_utf8_lossy(name).as_ref()))
                    .collect();
                let meta = std::sync::Arc::new(SharedBamMeta {
                    id,
                    path: path.to_path_buf(),
                    header,
                    contig_names,
                    index,
                });
                #[cfg(test)]
                if let Some(hook) = TEST_MID_COLD_OPEN_HOOK
                    .try_with(|h| h.try_borrow_mut().ok().and_then(|mut g| g.take()))
                    .ok()
                    .flatten()
                {
                    hook();
                }
                if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) == cur_epoch {
                    *slot_guard = Some(std::sync::Arc::clone(&meta));
                    insert_bam_meta_tls(cur_epoch, &meta);
                }
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

/// Reads reference sequence lengths (`SN` -> `LN`) from a BAM/CRAM header,
/// reusing `META_CACHE` / `CRAM_META_CACHE` so repeated calls on the same
/// file do not re-open or re-parse the header.
pub fn read_reference_lengths(
    path: &Path,
    config: &crate::config::Config,
) -> Result<HashMap<String, usize>, crate::types::VardictError> {
    if is_cram_path(path) {
        let meta = acquire_cram_meta(path, config)?;
        let mut lengths = HashMap::with_capacity(meta.header.reference_sequences().len());
        for (name, ref_seq) in meta.header.reference_sequences().iter() {
            lengths.insert(
                String::from_utf8_lossy(name).into_owned(),
                ref_seq.length().get(),
            );
        }
        Ok(lengths)
    } else {
        let (meta, reader) = acquire_bam_meta_and_reader(path)?;
        return_pooled_reader(path, std::sync::Arc::clone(&meta), reader);
        let mut lengths = HashMap::with_capacity(meta.header.reference_sequences().len());
        for (name, ref_seq) in meta.header.reference_sequences().iter() {
            lengths.insert(
                String::from_utf8_lossy(name).into_owned(),
                ref_seq.length().get(),
            );
        }
        Ok(lengths)
    }
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
        if let Err(e) = checked_cigar_op_len(op.len()) {
            return Some(Err(e));
        }
        if op.kind().consumes_reference() {
            span = match span.checked_add(op.len()) {
                Some(s) => s,
                None => {
                    return Some(Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "BAM CIGAR alignment_span overflow",
                    )));
                }
            };
        }
    }

    Some(checked_alignment_end(start, span))
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
    meta: std::sync::Arc<SharedBamMeta>,
    reader: Option<RawBamReader>,
    tainted: bool,
    reference_sequence_id: usize,
    interval: noodles::core::region::Interval,
    chunks: std::vec::IntoIter<noodles::csi::binning_index::index::reference_sequence::bin::Chunk>,
    chunk_end: Option<noodles::bgzf::VirtualPosition>,
    record: noodles::bam::Record,
    cigar_scratch: Vec<(u32, CigarKind)>,
    sa_buf: String,
    mc_buf: String,
    mate_ref_buf: String,
}

impl Drop for BamFeedInner {
    fn drop(&mut self) {
        if self.tainted {
            return;
        }
        let Some(reader) = self.reader.take() else {
            return;
        };
        let meta = std::sync::Arc::clone(&self.meta);
        return_pooled_reader(&self.meta.path, meta, reader);
    }
}

type BamSpanBounds = (
    Option<noodles::core::Position>,
    Option<noodles::core::Position>,
);

/// Ports noodles-bam 0.96 `io::reader::query::intersects` while decoding
/// `record.cigar()` once into `cigar_scratch` and computing `(alignment_start, alignment_end)`
/// (`span == 0 => end = start`).
fn intersects_record(
    record: &noodles::bam::Record,
    reference_sequence_id: usize,
    region_interval: noodles::core::region::Interval,
    cigar_scratch: &mut Vec<(u32, CigarKind)>,
) -> Result<Option<BamSpanBounds>, std::io::Error> {
    let Some(id) = record.reference_sequence_id().transpose()? else {
        return Ok(None);
    };

    if id != reference_sequence_id {
        return Ok(None);
    }

    let unbounded = region_interval.start().is_none() && region_interval.end().is_none();
    let alignment_start = match record.alignment_start().transpose() {
        Ok(pos) => pos,
        Err(e) if !unbounded => return Err(e),
        Err(_) => None,
    };
    if !unbounded && alignment_start.is_none() {
        return Ok(None);
    }

    let raw_cigar = record.cigar();
    cigar_scratch.clear();
    cigar_scratch.reserve(raw_cigar.len());
    let mut span = 0usize;
    let mut cigar_err = false;
    for result in raw_cigar.iter() {
        match result {
            Ok(op) => {
                let op_len = checked_cigar_op_len(op.len())?;
                if op.kind().consumes_reference() {
                    span = span.checked_add(op.len()).ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "BAM CIGAR alignment_span overflow",
                        )
                    })?;
                }
                cigar_scratch.push((op_len, op.kind()));
            }
            Err(e) => {
                if !unbounded {
                    return Err(e);
                }
                cigar_err = true;
                break;
            }
        }
    }

    let alignment_end = match alignment_start {
        Some(start) if !cigar_err => Some(checked_alignment_end(start, span)?),
        _ => None,
    };

    if unbounded {
        Ok(Some((alignment_start, alignment_end)))
    } else {
        match (alignment_start, alignment_end) {
            (Some(start), Some(end)) if region_interval.intersects((start..=end).into()) => {
                Ok(Some((Some(start), Some(end))))
            }
            _ => Ok(None),
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
                        let n = match read_bounded_bam_record(reader.get_mut(), &mut self.record) {
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
                            &mut self.cigar_scratch,
                        ) {
                            Ok(Some((alignment_start, alignment_end))) => {
                                out.fill_from_bam_feed_record(
                                    &self.record,
                                    &self.meta.contig_names,
                                    self.reference_sequence_id,
                                    alignment_start,
                                    alignment_end,
                                    &mut self.cigar_scratch,
                                    &mut self.sa_buf,
                                    &mut self.mc_buf,
                                    &mut self.mate_ref_buf,
                                );
                                return Ok(true);
                            }
                            Ok(None) => continue,
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

impl std::fmt::Debug for BamFeed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BamFeed")
            .field("path", &self.inner.meta.path)
            .field("reference_sequence_id", &self.inner.reference_sequence_id)
            .field("interval", &self.inner.interval)
            .finish_non_exhaustive()
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

const FASTA_BUF_CAPACITY: usize = 1024 * 1024;
const DEFAULT_FASTA_CACHE_BUDGET_BYTES: usize = 512 * 1024 * 1024;
const MAX_POOLED_FASTA_ADAPTERS: usize = 4;
/// Trailing `b'N'` padding appended to every CRAM reference contig loaded into
/// [`PreallocatedFastaAdapter`] so `noodles-cram` can decode reads overhanging
/// the 3' end of a chromosome (`[LN - k, LN + m]`) without panicking inside
/// `get_slice_reference_sequence` or `Record::sequence`, matching HTSJDK/BAM
/// semantics (`Issue K` R1b).
const CRAM_REF_OVERHANG_PAD: usize = 65_536;

struct PreallocatedFastaAdapter {
    reader: noodles::fasta::io::Reader<noodles::fasta::io::BufReader<std::fs::File>>,
    index: noodles::fasta::fai::Index,
    raw_file_len: Option<u64>,
}

impl PreallocatedFastaAdapter {
    fn from_buf_reader_and_index(
        buf_reader: noodles::fasta::io::BufReader<std::fs::File>,
        index: noodles::fasta::fai::Index,
    ) -> Self {
        let raw_file_len = match buf_reader {
            noodles::fasta::io::BufReader::Uncompressed(ref br) => {
                br.get_ref().metadata().ok().map(|m| m.len())
            }
            noodles::fasta::io::BufReader::Bgzf(_) => None,
        };
        Self {
            reader: noodles::fasta::io::Reader::new(buf_reader),
            index,
            raw_file_len,
        }
    }

    fn open(fasta_path: &Path) -> Result<Self, crate::types::VardictError> {
        let buf_reader = crate::reference::open_indexed_fasta_buf_reader_with_capacity(
            fasta_path,
            FASTA_BUF_CAPACITY,
        )?;
        let fai_path = fasta_fai_path(fasta_path);
        let index = crate::reference::read_validated_fai_index(&fai_path)?;
        Ok(Self::from_buf_reader_and_index(buf_reader, index))
    }

    fn open_with_index(
        fasta_path: &Path,
        index: noodles::fasta::fai::Index,
    ) -> Result<Self, crate::types::VardictError> {
        let buf_reader = crate::reference::open_indexed_fasta_buf_reader_with_capacity(
            fasta_path,
            FASTA_BUF_CAPACITY,
        )?;
        Ok(Self::from_buf_reader_and_index(buf_reader, index))
    }
}

impl noodles::fasta::repository::Adapter for PreallocatedFastaAdapter {
    fn get(&mut self, name: &[u8]) -> Option<std::io::Result<noodles::fasta::Record>> {
        use std::io::{BufRead as _, Seek as _, SeekFrom};
        const MAX_PREALLOC_CAP: usize = 64 * 1024 * 1024;

        let Some(fai_rec) = self.index.as_ref().iter().find(|r| r.name() == name) else {
            return Some(Err(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!(
                    "reference sequence '{}' not found in FASTA index",
                    String::from_utf8_lossy(name)
                ),
            )));
        };
        let Ok(cap) = usize::try_from(fai_rec.length()) else {
            return Some(Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "FASTA contig '{}' length ({}) overflows usize",
                    String::from_utf8_lossy(name),
                    fai_rec.length()
                ),
            )));
        };

        let region = noodles::core::Region::new(name, ..);
        let pos = match self.index.query(&region) {
            Ok(p) => p,
            Err(e) => return Some(Err(e)),
        };
        if let Err(e) = self.reader.get_mut().seek(SeekFrom::Start(pos)) {
            return Some(Err(e));
        }

        let max_raw_len = cap.max(1).saturating_add(CRAM_REF_OVERHANG_PAD);
        let prealloc_bases = match self.raw_file_len {
            Some(fl) => {
                let rem_file = usize::try_from(fl.saturating_sub(pos)).unwrap_or(usize::MAX);
                cap.min(rem_file).min(DEFAULT_FASTA_CACHE_BUDGET_BYTES)
            }
            None => cap.min(MAX_PREALLOC_CAP),
        };
        let mut raw_sequence =
            Vec::with_capacity(prealloc_bases.saturating_add(CRAM_REF_OVERHANG_PAD));
        let buf_reader = self.reader.get_mut();
        let mut at_line_start = true;
        let mut pending_cr = false;
        loop {
            let chunk = match buf_reader.fill_buf() {
                Ok([]) => {
                    if pending_cr {
                        if raw_sequence.len().saturating_add(1) > max_raw_len {
                            return Some(Err(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                format!(
                                    "FASTA sequence for '{}' exceeds .fai length ({cap} bp) + overhang bound",
                                    String::from_utf8_lossy(name)
                                ),
                            )));
                        }
                        raw_sequence.push(b'\r');
                    }
                    break;
                }
                Ok(s) => s,
                Err(e) => return Some(Err(e)),
            };
            let mut rem = chunk;
            let mut consumed = 0usize;
            let mut hit_next_record = false;

            if pending_cr {
                pending_cr = false;
                if rem.first() == Some(&b'\n') {
                    rem = rem.get(1..).unwrap_or(&[]);
                    consumed += 1;
                    at_line_start = true;
                } else {
                    if raw_sequence.len().saturating_add(1) > max_raw_len {
                        return Some(Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!(
                                "FASTA sequence for '{}' exceeds .fai length ({cap} bp) + overhang bound",
                                String::from_utf8_lossy(name)
                            ),
                        )));
                    }
                    raw_sequence.push(b'\r');
                }
            }

            while !rem.is_empty() {
                if at_line_start {
                    let ws = rem
                        .iter()
                        .take_while(|&&b| b == b'\n' || b == b'\r')
                        .count();
                    if ws > 0 {
                        rem = rem.get(ws..).unwrap_or(&[]);
                        consumed += ws;
                        if rem.is_empty() {
                            break;
                        }
                    }
                    if rem.first() == Some(&b'>') {
                        hit_next_record = true;
                        break;
                    }
                    at_line_start = false;
                }

                if let Some(nl_pos) = rem.iter().position(|&b| b == b'\n') {
                    let line_raw = rem.get(..nl_pos).unwrap_or(&[]);
                    let line = line_raw.strip_suffix(b"\r").unwrap_or(line_raw);
                    if raw_sequence.len().saturating_add(line.len()) > max_raw_len {
                        return Some(Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!(
                                "FASTA sequence for '{}' exceeds .fai length ({cap} bp) + overhang bound",
                                String::from_utf8_lossy(name)
                            ),
                        )));
                    }
                    raw_sequence.extend_from_slice(line);
                    let step = nl_pos + 1;
                    rem = rem.get(step..).unwrap_or(&[]);
                    consumed += step;
                    at_line_start = true;
                } else {
                    let line = match rem.strip_suffix(b"\r") {
                        Some(stripped) => {
                            pending_cr = true;
                            stripped
                        }
                        None => rem,
                    };
                    if raw_sequence.len().saturating_add(line.len()) > max_raw_len {
                        return Some(Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!(
                                "FASTA sequence for '{}' exceeds .fai length ({cap} bp) + overhang bound",
                                String::from_utf8_lossy(name)
                            ),
                        )));
                    }
                    raw_sequence.extend_from_slice(line);
                    consumed += rem.len();
                    break;
                }
            }

            buf_reader.consume(consumed);
            if hit_next_record {
                break;
            }
        }
        raw_sequence.resize(
            raw_sequence.len().saturating_add(CRAM_REF_OVERHANG_PAD),
            b'N',
        );
        let definition = noodles::fasta::record::Definition::new(region.to_string(), None);
        let sequence = noodles::fasta::record::Sequence::from(raw_sequence);
        Some(Ok(noodles::fasta::Record::new(definition, sequence)))
    }
}

struct LoadedContig {
    repo: std::sync::Mutex<noodles::fasta::Repository>,
    seq: Option<std::sync::Arc<noodles::fasta::record::Sequence>>,
    byte_len: usize,
}

struct ContigSlot {
    name: Vec<u8>,
    loaded: std::sync::OnceLock<Result<LoadedContig, (std::io::ErrorKind, String)>>,
    byte_len: std::sync::atomic::AtomicUsize,
    delegated_bytes: std::sync::atomic::AtomicUsize,
    counted_in_registry: std::sync::atomic::AtomicBool,
}

impl std::fmt::Debug for ContigSlot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ContigSlot")
            .field("name", &String::from_utf8_lossy(&self.name))
            .field(
                "byte_len",
                &self.byte_len.load(std::sync::atomic::Ordering::Relaxed),
            )
            .field(
                "delegated_bytes",
                &self
                    .delegated_bytes
                    .load(std::sync::atomic::Ordering::Relaxed),
            )
            .finish_non_exhaustive()
    }
}

struct RegistryDelegatingAdapter {
    primary_name: Vec<u8>,
    primary_record: Option<std::io::Result<noodles::fasta::Record>>,
    shared: std::sync::Weak<SharedFastaRepo>,
}

impl RegistryDelegatingAdapter {
    fn from_primary_seq(
        primary_name: Vec<u8>,
        seq: Option<&std::sync::Arc<noodles::fasta::record::Sequence>>,
        shared: std::sync::Weak<SharedFastaRepo>,
    ) -> Self {
        let primary_record = seq.map(|s| {
            let def = noodles::fasta::record::Definition::new(
                String::from_utf8_lossy(&primary_name).into_owned(),
                None,
            );
            let cloned_seq = noodles::fasta::record::Sequence::from(s.as_ref().as_ref().to_vec());
            Ok(noodles::fasta::Record::new(def, cloned_seq))
        });
        Self {
            primary_name,
            primary_record,
            shared,
        }
    }
}

impl noodles::fasta::repository::Adapter for RegistryDelegatingAdapter {
    fn get(&mut self, name: &[u8]) -> Option<std::io::Result<noodles::fasta::Record>> {
        if name == self.primary_name.as_slice() {
            if let Some(rec) = self.primary_record.take() {
                return Some(rec);
            }
        }
        let Some(shared) = self.shared.upgrade() else {
            return Some(Err(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!(
                    "FASTA repository dropped while resolving contig '{}'",
                    String::from_utf8_lossy(name)
                ),
            )));
        };
        let slot = shared.warm_contig_slot(name);
        match slot.loaded.get() {
            Some(Ok(loaded)) => {
                let Some(seq_arc) = loaded.seq.as_ref() else {
                    return Some(Err(std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!(
                            "reference sequence '{}' has no loaded bases in FASTA repository",
                            String::from_utf8_lossy(name)
                        ),
                    )));
                };
                let raw_bytes = seq_arc.as_ref().as_ref().to_vec();
                shared.charge_delegated_bytes(&self.primary_name, raw_bytes.len());
                let def = noodles::fasta::record::Definition::new(
                    String::from_utf8_lossy(name).into_owned(),
                    None,
                );
                let seq = noodles::fasta::record::Sequence::from(raw_bytes);
                Some(Ok(noodles::fasta::Record::new(def, seq)))
            }
            Some(Err((kind, msg))) => Some(Err(std::io::Error::new(*kind, msg.clone()))),
            None => Some(Err(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!(
                    "reference sequence '{}' not initialized in FASTA repository",
                    String::from_utf8_lossy(name)
                ),
            ))),
        }
    }
}

struct FastaRegistryState {
    slots: HashMap<Vec<u8>, std::sync::Arc<ContigSlot>>,
    order: Vec<Vec<u8>>,
    resident_bytes: usize,
    inflight_bytes: usize,
    inflight_loads: usize,
    budget_bytes: usize,
}

impl FastaRegistryState {
    fn new(budget_bytes: usize) -> Self {
        Self {
            slots: HashMap::new(),
            order: Vec::new(),
            resident_bytes: 0,
            inflight_bytes: 0,
            inflight_loads: 0,
            budget_bytes,
        }
    }

    fn touch_order(&mut self, contig: &[u8]) {
        if let Some(pos) = self.order.iter().position(|k| k.as_slice() == contig) {
            let key = self.order.remove(pos);
            self.order.push(key);
        } else {
            self.order.push(contig.to_vec());
        }
    }

    fn evict_unpinned_to_budget(&mut self) -> bool {
        self.evict_unpinned_for_incoming(0)
    }

    fn evict_unpinned_for_incoming(&mut self, incoming_bytes: usize) -> bool {
        let mut evicted_any = false;
        while self
            .resident_bytes
            .saturating_add(self.inflight_bytes)
            .saturating_add(incoming_bytes)
            > self.budget_bytes
        {
            let candidate_pos = self.order.iter().position(|key| {
                self.slots
                    .get(key)
                    .is_some_and(|slot_arc| std::sync::Arc::strong_count(slot_arc) == 1)
            });
            let Some(pos) = candidate_pos else {
                break;
            };
            let key = self.order.remove(pos);
            if let Some(evicted) = self.slots.remove(&key) {
                evicted_any = true;
                if evicted
                    .counted_in_registry
                    .swap(false, std::sync::atomic::Ordering::Relaxed)
                {
                    let len = evicted.byte_len.load(std::sync::atomic::Ordering::Relaxed);
                    let delegated = evicted
                        .delegated_bytes
                        .swap(0, std::sync::atomic::Ordering::Relaxed);
                    self.resident_bytes = self
                        .resident_bytes
                        .saturating_sub(len.saturating_add(delegated));
                }
            }
        }
        evicted_any
    }
}

struct InflightFastaLoadGuard<'a> {
    repo: &'a SharedFastaRepo,
    reserved_bytes: usize,
}

impl Drop for InflightFastaLoadGuard<'_> {
    fn drop(&mut self) {
        let mut state = self.repo.registry.lock().unwrap_or_else(|e| e.into_inner());
        state.inflight_loads = state.inflight_loads.saturating_sub(1);
        state.inflight_bytes = state.inflight_bytes.saturating_sub(self.reserved_bytes);
        drop(state);
        self.repo.budget_cvar.notify_all();
    }
}

struct SharedFastaRepo {
    fasta_path: PathBuf,
    fai_index: noodles::fasta::fai::Index,
    clone_gate: std::sync::Mutex<()>,
    budget_cvar: std::sync::Condvar,
    adapter_pool: std::sync::Mutex<Vec<PreallocatedFastaAdapter>>,
    registry: std::sync::Mutex<FastaRegistryState>,
    #[cfg(test)]
    load_count: std::sync::atomic::AtomicUsize,
}

impl SharedFastaRepo {
    #[inline]
    fn unpadded_slice(seq: &noodles::fasta::record::Sequence) -> &[u8] {
        let bytes: &[u8] = seq.as_ref();
        let raw_len = bytes.len().saturating_sub(CRAM_REF_OVERHANG_PAD);
        bytes.get(..raw_len).unwrap_or(bytes)
    }

    fn new(
        fasta_path: PathBuf,
        adapter: PreallocatedFastaAdapter,
        budget_bytes: usize,
    ) -> std::sync::Arc<Self> {
        let fai_index = adapter.index.clone();
        std::sync::Arc::new(Self {
            fasta_path,
            fai_index,
            clone_gate: std::sync::Mutex::new(()),
            budget_cvar: std::sync::Condvar::new(),
            adapter_pool: std::sync::Mutex::new(vec![adapter]),
            registry: std::sync::Mutex::new(FastaRegistryState::new(budget_bytes)),
            #[cfg(test)]
            load_count: std::sync::atomic::AtomicUsize::new(0),
        })
    }

    fn pop_or_open_adapter(&self) -> Result<PreallocatedFastaAdapter, std::io::Error> {
        if let Some(adapter) = self
            .adapter_pool
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .pop()
        {
            return Ok(adapter);
        }
        PreallocatedFastaAdapter::open_with_index(&self.fasta_path, self.fai_index.clone()).map_err(
            |e| match e {
                crate::types::VardictError::Io(io_err) => io_err,
                other => std::io::Error::other(other.to_string()),
            },
        )
    }

    fn return_adapter(&self, adapter: PreallocatedFastaAdapter) {
        let mut pool = self.adapter_pool.lock().unwrap_or_else(|e| e.into_inner());
        if pool.len() < MAX_POOLED_FASTA_ADAPTERS {
            pool.push(adapter);
        }
    }

    fn load_contig_for_slot(
        self: &std::sync::Arc<Self>,
        contig: &[u8],
    ) -> Result<LoadedContig, (std::io::ErrorKind, String)> {
        const MAX_CONCURRENT_FASTA_LOADS: usize = 4;

        if contig.is_empty() {
            let repo = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                primary_name: Vec::new(),
                primary_record: None,
                shared: std::sync::Arc::downgrade(self),
            });
            return Ok(LoadedContig {
                repo: std::sync::Mutex::new(repo),
                seq: None,
                byte_len: 0,
            });
        }
        // Reserve `incoming_bytes` against the 512 MiB FASTA cache budget before reading
        // so multiple workers load non-overlapping chromosomes in parallel without breaching
        // the cache ceiling, while `clone_gate` serializes only `Repository::get`'s short
        // `Sequence::clone` step so two transient `2 * contig_len` copies never overlap.
        let incoming_bytes = self
            .fai_index
            .as_ref()
            .iter()
            .find(|r| r.name() == contig)
            .and_then(|r| usize::try_from(r.length()).ok())
            .map_or(0, |len| len.saturating_add(CRAM_REF_OVERHANG_PAD));
        let _inflight_guard = {
            let mut state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
            loop {
                let _ = state.evict_unpinned_for_incoming(incoming_bytes);
                let fits_budget = state
                    .resident_bytes
                    .saturating_add(state.inflight_bytes)
                    .saturating_add(incoming_bytes)
                    <= state.budget_bytes;
                if state.inflight_loads == 0
                    || (state.inflight_loads < MAX_CONCURRENT_FASTA_LOADS && fits_budget)
                {
                    state.inflight_loads = state.inflight_loads.saturating_add(1);
                    state.inflight_bytes = state.inflight_bytes.saturating_add(incoming_bytes);
                    break InflightFastaLoadGuard {
                        repo: self,
                        reserved_bytes: incoming_bytes,
                    };
                }
                state = self
                    .budget_cvar
                    .wait(state)
                    .unwrap_or_else(|e| e.into_inner());
            }
        };
        let mut adapter = self
            .pop_or_open_adapter()
            .map_err(|e| (e.kind(), e.to_string()))?;
        #[cfg(test)]
        self.load_count
            .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        let record_opt = noodles::fasta::repository::Adapter::get(&mut adapter, contig);
        self.return_adapter(adapter);

        match record_opt {
            Some(Ok(record)) => {
                // Charge full retained sequence length (including CRAM_REF_OVERHANG_PAD).
                let byte_len = record.sequence().len();
                let seq;
                let repo = {
                    let _cg = self.clone_gate.lock().unwrap_or_else(|e| e.into_inner());
                    let r = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                        primary_name: contig.to_vec(),
                        primary_record: Some(Ok(record)),
                        shared: std::sync::Arc::downgrade(self),
                    });
                    seq = r.get(contig).and_then(|res| res.ok());
                    r
                };
                Ok(LoadedContig {
                    repo: std::sync::Mutex::new(repo),
                    seq,
                    byte_len,
                })
            }
            Some(Err(e)) => Err((e.kind(), e.to_string())),
            None => {
                let repo = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                    primary_name: contig.to_vec(),
                    primary_record: None,
                    shared: std::sync::Arc::downgrade(self),
                });
                Ok(LoadedContig {
                    repo: std::sync::Mutex::new(repo),
                    seq: None,
                    byte_len: 0,
                })
            }
        }
    }

    fn charge_delegated_bytes(&self, owner_contig: &[u8], bytes: usize) {
        if bytes == 0 {
            return;
        }
        let mut state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(owner_slot) = state.slots.get(owner_contig).cloned() {
            owner_slot
                .delegated_bytes
                .fetch_add(bytes, std::sync::atomic::Ordering::Relaxed);
            if owner_slot
                .counted_in_registry
                .load(std::sync::atomic::Ordering::Acquire)
            {
                state.resident_bytes = state.resident_bytes.saturating_add(bytes);
                state.evict_unpinned_to_budget();
            }
        }
    }

    fn reclaim_delegated_bytes(self: &std::sync::Arc<Self>, owner_contig: &[u8]) {
        let mut state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
        let Some(owner_slot) = state.slots.get(owner_contig).cloned() else {
            return;
        };
        let delegated = owner_slot
            .delegated_bytes
            .swap(0, std::sync::atomic::Ordering::Relaxed);
        if delegated == 0 {
            return;
        }
        if owner_slot
            .counted_in_registry
            .load(std::sync::atomic::Ordering::Acquire)
        {
            state.resident_bytes = state.resident_bytes.saturating_sub(delegated);
        }
        if let Some(Ok(loaded)) = owner_slot.loaded.get() {
            let fresh_repo =
                noodles::fasta::Repository::new(RegistryDelegatingAdapter::from_primary_seq(
                    owner_slot.name.clone(),
                    loaded.seq.as_ref(),
                    std::sync::Arc::downgrade(self),
                ));
            let _ = fresh_repo.get(&owner_slot.name);
            *loaded.repo.lock().unwrap_or_else(|e| e.into_inner()) = fresh_repo;
        }
    }

    fn warm_contig_slot(self: &std::sync::Arc<Self>, contig: &[u8]) -> std::sync::Arc<ContigSlot> {
        let slot = {
            let mut state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
            if let Some(existing) = state.slots.get(contig).cloned() {
                state.touch_order(contig);
                state.evict_unpinned_to_budget();
                existing
            } else {
                let new_slot = std::sync::Arc::new(ContigSlot {
                    name: contig.to_vec(),
                    loaded: std::sync::OnceLock::new(),
                    byte_len: std::sync::atomic::AtomicUsize::new(0),
                    delegated_bytes: std::sync::atomic::AtomicUsize::new(0),
                    counted_in_registry: std::sync::atomic::AtomicBool::new(false),
                });
                state
                    .slots
                    .insert(contig.to_vec(), std::sync::Arc::clone(&new_slot));
                state.order.push(contig.to_vec());
                new_slot
            }
        };

        let loaded_res = slot
            .loaded
            .get_or_init(|| self.load_contig_for_slot(&slot.name));
        let byte_len = match loaded_res {
            Ok(loaded) => loaded.byte_len,
            Err(_) => 0,
        };
        slot.byte_len
            .store(byte_len, std::sync::atomic::Ordering::Relaxed);

        if !slot
            .counted_in_registry
            .load(std::sync::atomic::Ordering::Acquire)
        {
            let mut state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
            let still_registered = state
                .slots
                .get(contig)
                .is_some_and(|s| std::sync::Arc::ptr_eq(s, &slot));
            if still_registered
                && !slot
                    .counted_in_registry
                    .swap(true, std::sync::atomic::Ordering::AcqRel)
            {
                let delegated = slot
                    .delegated_bytes
                    .load(std::sync::atomic::Ordering::Relaxed);
                state.resident_bytes = state
                    .resident_bytes
                    .saturating_add(byte_len.saturating_add(delegated));
            }
            state.evict_unpinned_to_budget();
        }

        slot
    }

    fn warm_contig(
        self: &std::sync::Arc<Self>,
        contig: &[u8],
    ) -> (
        noodles::fasta::Repository,
        Option<std::sync::Arc<noodles::fasta::record::Sequence>>,
    ) {
        let slot = self.warm_contig_slot(contig);
        Self::extract_from_slot(&slot)
    }

    fn extract_from_slot(
        slot: &std::sync::Arc<ContigSlot>,
    ) -> (
        noodles::fasta::Repository,
        Option<std::sync::Arc<noodles::fasta::record::Sequence>>,
    ) {
        match slot.loaded.get() {
            Some(Ok(loaded)) => {
                let repo = loaded
                    .repo
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .clone();
                (repo, loaded.seq.clone())
            }
            _ => {
                let empty_repo = noodles::fasta::Repository::new(
                    noodles::fasta::repository::adapters::Empty::new(),
                );
                (empty_repo, None)
            }
        }
    }

    fn evict_unpinned(&self) {
        let evicted = {
            let mut state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
            state.evict_unpinned_to_budget()
        };
        if evicted {
            self.budget_cvar.notify_all();
        }
    }

    #[cfg(test)]
    fn set_budget_bytes(&self, budget_bytes: usize) {
        let mut state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
        state.budget_bytes = budget_bytes;
        state.evict_unpinned_to_budget();
    }

    #[cfg(test)]
    fn resident_slot(&self, contig: &[u8]) -> Option<std::sync::Arc<ContigSlot>> {
        let state = self.registry.lock().unwrap_or_else(|e| e.into_inner());
        state.slots.get(contig).cloned()
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
    let repo = SharedFastaRepo::new(
        fasta_path.to_path_buf(),
        adapter,
        DEFAULT_FASTA_CACHE_BUDGET_BYTES,
    );
    map.insert(fasta_path.to_path_buf(), std::sync::Arc::clone(&repo));
    Ok(repo)
}

#[derive(Clone, Copy, Debug)]
struct CraiSpanEntry {
    start: usize,
    end: usize,
    offset: u64,
    orig_idx: usize,
}

#[derive(Clone, Debug, Default)]
struct CraiRefGroup {
    unbounded_offsets: Vec<u64>,
    valid_spans: Vec<CraiSpanEntry>,
    prefix_max_end: Vec<usize>,
    is_sorted_by_start: bool,
}

#[derive(Clone, Debug, Default)]
struct IndexedCrai {
    by_ref: HashMap<usize, CraiRefGroup>,
}

impl IndexedCrai {
    fn from_index(index: &noodles::cram::crai::Index) -> Self {
        let mut by_ref: HashMap<usize, CraiRefGroup> = HashMap::new();
        let mut unbounded_seen: std::collections::HashSet<(usize, u64)> =
            std::collections::HashSet::new();
        for (orig_idx, rec) in index.iter().enumerate() {
            let Some(ref_id) = rec.reference_sequence_id() else {
                continue;
            };
            let group = by_ref.entry(ref_id).or_insert_with(|| CraiRefGroup {
                unbounded_offsets: Vec::new(),
                valid_spans: Vec::new(),
                prefix_max_end: Vec::new(),
                is_sorted_by_start: true,
            });
            let offset = rec.offset();
            if unbounded_seen.insert((ref_id, offset)) {
                group.unbounded_offsets.push(offset);
            }
            let Some(start_pos) = rec.alignment_start() else {
                continue;
            };
            let span = rec.alignment_span();
            if span == 0 {
                continue;
            }
            let Some(end_pos) = start_pos.checked_add(span - 1) else {
                continue;
            };
            let start = usize::from(start_pos);
            let end = usize::from(end_pos);
            if let Some(prev) = group.valid_spans.last() {
                if start < prev.start {
                    group.is_sorted_by_start = false;
                }
            }
            group.valid_spans.push(CraiSpanEntry {
                start,
                end,
                offset,
                orig_idx,
            });
        }

        for group in by_ref.values_mut() {
            if !group.is_sorted_by_start {
                group.valid_spans.sort_by_key(|e| (e.start, e.orig_idx));
            }
            let mut running_max = 0usize;
            group.prefix_max_end.reserve(group.valid_spans.len());
            for entry in &group.valid_spans {
                if entry.end > running_max {
                    running_max = entry.end;
                }
                group.prefix_max_end.push(running_max);
            }
        }

        Self { by_ref }
    }

    #[inline]
    fn has_ref_containers(&self, reference_sequence_id: usize) -> bool {
        self.by_ref
            .get(&reference_sequence_id)
            .is_some_and(|g| !g.unbounded_offsets.is_empty())
    }

    fn query_offsets(
        &self,
        reference_sequence_id: usize,
        region_interval: noodles::core::region::Interval,
    ) -> Vec<u64> {
        let Some(group) = self.by_ref.get(&reference_sequence_id) else {
            return Vec::new();
        };
        if region_interval.start().is_none() && region_interval.end().is_none() {
            crate::count!(CraiRowsExamined, group.unbounded_offsets.len());
            crate::count!(CraiRowsMatched, group.unbounded_offsets.len());
            return group.unbounded_offsets.clone();
        }

        let q_start = region_interval.start().map(usize::from).unwrap_or(1);
        let q_end = region_interval.end().map(usize::from).unwrap_or(usize::MAX);
        if q_start > q_end {
            return Vec::new();
        }

        let upper = group.valid_spans.partition_point(|e| e.start <= q_end);
        let prefix_slice = group.prefix_max_end.get(..upper).unwrap_or(&[]);
        let lower = prefix_slice.partition_point(|&max_end| max_end < q_start);
        let window = group.valid_spans.get(lower..upper).unwrap_or(&[]);
        crate::count!(CraiRowsExamined, window.len());

        let mut container_offsets: Vec<u64> = Vec::new();
        if group.is_sorted_by_start {
            for entry in window {
                if entry.end >= q_start {
                    crate::count!(CraiRowsMatched);
                    if !container_offsets.contains(&entry.offset) {
                        container_offsets.push(entry.offset);
                    }
                }
            }
        } else {
            let mut matched: Vec<CraiSpanEntry> = window
                .iter()
                .copied()
                .filter(|e| e.end >= q_start)
                .collect();
            crate::count!(CraiRowsMatched, matched.len());
            matched.sort_by_key(|e| e.orig_idx);
            for entry in matched {
                if !container_offsets.contains(&entry.offset) {
                    container_offsets.push(entry.offset);
                }
            }
        }
        container_offsets
    }
}

#[doc(hidden)]
pub struct SharedCramMeta {
    id: u64,
    file_len: Option<u64>,
    header: noodles::sam::Header,
    contig_names: Vec<std::sync::Arc<str>>,
    indexed_crai: IndexedCrai,
    spans_by_offset: HashMap<u64, Vec<(usize, usize, usize)>>,
    fasta_repo: std::sync::Arc<SharedFastaRepo>,
    #[cfg(test)]
    decodes: std::sync::atomic::AtomicUsize,
}

impl std::fmt::Debug for SharedCramMeta {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SharedCramMeta")
            .field("id", &self.id)
            .field("file_len", &self.file_len)
            .finish_non_exhaustive()
    }
}

type CramMetaSlot = std::sync::Arc<std::sync::Mutex<Option<std::sync::Arc<SharedCramMeta>>>>;
type CramMetaKey = (PathBuf, PathBuf);
type CramMetaTlsEntry = (u64, PathBuf, PathBuf, std::sync::Arc<SharedCramMeta>);
type RawCramReader = noodles::cram::io::Reader<std::fs::File>;
type CramReaderPoolMap = ReaderPoolMap<SharedCramMeta, RawCramReader>;

static CRAM_META_CACHE: std::sync::OnceLock<std::sync::Mutex<HashMap<CramMetaKey, CramMetaSlot>>> =
    std::sync::OnceLock::new();

thread_local! {
    static CRAM_READER_POOL: RefCell<CramReaderPoolMap> = RefCell::new(HashMap::new());
    static CRAM_META_TLS: RefCell<[Option<CramMetaTlsEntry>; 2]> = const { RefCell::new([None, None]) };
}

#[inline]
fn lookup_cram_meta_tls(
    path: &Path,
    fasta_path: &Path,
    cur_epoch: u64,
) -> Option<std::sync::Arc<SharedCramMeta>> {
    CRAM_META_TLS
        .try_with(|tls| {
            let mut slots = tls.try_borrow_mut().ok()?;
            let [s0, s1] = &mut *slots;
            if let Some((ep, cached_cram, cached_fasta, meta)) = s0.as_ref() {
                if tls_epoch_valid(*ep, cur_epoch, meta.id)
                    && cached_cram.as_path() == path
                    && cached_fasta.as_path() == fasta_path
                {
                    return Some(std::sync::Arc::clone(meta));
                }
            }
            if let Some((ep, cached_cram, cached_fasta, meta)) = s1.as_ref() {
                if tls_epoch_valid(*ep, cur_epoch, meta.id)
                    && cached_cram.as_path() == path
                    && cached_fasta.as_path() == fasta_path
                {
                    let hit = std::sync::Arc::clone(meta);
                    std::mem::swap(s0, s1);
                    return Some(hit);
                }
            }
            None
        })
        .ok()
        .flatten()
}

#[inline]
fn insert_cram_meta_tls(
    cur_epoch: u64,
    path: &Path,
    fasta_path: &Path,
    meta: &std::sync::Arc<SharedCramMeta>,
) {
    let _ = CRAM_META_TLS.try_with(|tls| {
        if let Ok(mut slots) = tls.try_borrow_mut() {
            let [s0, s1] = &mut *slots;
            if s0
                .as_ref()
                .is_some_and(|(_, c, f, _)| c.as_path() == path && f.as_path() == fasta_path)
            {
                *s0 = Some((
                    cur_epoch,
                    path.to_path_buf(),
                    fasta_path.to_path_buf(),
                    std::sync::Arc::clone(meta),
                ));
            } else {
                *s1 = s0.take();
                *s0 = Some((
                    cur_epoch,
                    path.to_path_buf(),
                    fasta_path.to_path_buf(),
                    std::sync::Arc::clone(meta),
                ));
            }
        }
    });
}

fn pop_pooled_cram_reader(
    path: &Path,
    meta: &std::sync::Arc<SharedCramMeta>,
) -> Option<RawCramReader> {
    pop_from_pool(&CRAM_READER_POOL, path, meta)
}

fn return_pooled_cram_reader(
    path: &Path,
    meta: std::sync::Arc<SharedCramMeta>,
    reader: RawCramReader,
) {
    return_to_pool(&CRAM_READER_POOL, path, meta, reader);
}

fn prevalidate_cram_file_header(path: &Path, file: &mut std::fs::File) -> std::io::Result<()> {
    use std::io::{Read as _, Seek as _, SeekFrom};
    const MAX_HEADER_CONTAINER_HDR_BYTES: u64 = 65_536;
    const MAX_HEADER_LANDMARKS: usize = 4096;
    const FILE_HEADER_CONTENT_TYPE: u8 = 0;

    {
        let mut buf_reader = std::io::BufReader::new(&mut *file);
        let mut file_def = [0_u8; 26];
        buf_reader.read_exact(&mut file_def)?;
        if file_def.get(..4) != Some(b"CRAM") {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("invalid CRAM magic header in '{}'", path.display()),
            ));
        }

        let mut len_buf = [0_u8; 4];
        buf_reader.read_exact(&mut len_buf)?;
        let raw_container_len = i32::from_le_bytes(len_buf);
        let Ok(container_len) = usize::try_from(raw_container_len) else {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "negative CRAM header container length ({raw_container_len}) in '{}'",
                    path.display()
                ),
            ));
        };
        if container_len > MAX_CRAM_HEADER_BYTES {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM header container length ({container_len}) in '{}' exceeds maximum allowed ({MAX_CRAM_HEADER_BYTES} bytes)",
                    path.display()
                ),
            ));
        }

        {
            let mut hdr_take = (&mut buf_reader).take(MAX_HEADER_CONTAINER_HDR_BYTES);
            let _ref_seq_id = read_cram_itf8(&mut hdr_take)?;
            let _align_start = read_cram_itf8(&mut hdr_take)?;
            let _align_span = read_cram_itf8(&mut hdr_take)?;
            let _record_count = read_cram_itf8(&mut hdr_take)?;
            let _record_counter = read_cram_ltf8(&mut hdr_take)?;
            let _base_count = read_cram_ltf8(&mut hdr_take)?;
            let _block_count = read_cram_itf8(&mut hdr_take)?;
            let raw_landmarks_len = read_cram_itf8(&mut hdr_take)?;
            let landmarks_len = usize::try_from(raw_landmarks_len)
                .ok()
                .filter(|&n| n <= MAX_HEADER_LANDMARKS)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM header container landmarks count ({raw_landmarks_len}) in '{}' exceeds allowed range 0..={MAX_HEADER_LANDMARKS}",
                            path.display()
                        ),
                    )
                })?;
            for _ in 0..landmarks_len {
                let _ = read_cram_itf8(&mut hdr_take)?;
            }
            let mut crc_buf = [0_u8; 4];
            hdr_take.read_exact(&mut crc_buf)?;
        }

        let mut body_take = (&mut buf_reader).take(container_len as u64);
        let mut blk_prefix = [0_u8; 2];
        body_take.read_exact(&mut blk_prefix)?;
        let [raw_method, content_type] = blk_prefix;
        if raw_method > 1 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "invalid CRAM FileHeader compression method ({raw_method}) in '{}'",
                    path.display()
                ),
            ));
        }
        if content_type != FILE_HEADER_CONTENT_TYPE {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "invalid CRAM FileHeader block content type ({content_type}) in '{}'",
                    path.display()
                ),
            ));
        }
        let _content_id = read_cram_itf8(&mut body_take)?;
        let comp = usize::try_from(read_cram_itf8(&mut body_take)?)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let uncomp = usize::try_from(read_cram_itf8(&mut body_take)?)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        if comp > MAX_CRAM_HEADER_BYTES || uncomp > MAX_CRAM_HEADER_BYTES {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM FileHeader block size (compressed={comp}, uncompressed={uncomp}) in '{}' exceeds maximum allowed ({MAX_CRAM_HEADER_BYTES} bytes)",
                    path.display()
                ),
            ));
        }
        let rem_body = body_take.limit();
        if (comp as u64) > rem_body || (raw_method == 0 && (uncomp as u64) > rem_body) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM FileHeader block size (compressed={comp}, uncompressed={uncomp}) in '{}' exceeds container payload ({container_len} bytes)",
                    path.display()
                ),
            ));
        }

        if raw_method == 0 {
            if uncomp >= 4 {
                let mut sam_len_buf = [0_u8; 4];
                body_take.read_exact(&mut sam_len_buf)?;
                let raw_sam_len = i32::from_le_bytes(sam_len_buf);
                let sam_len = usize::try_from(raw_sam_len)
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
                if sam_len > MAX_CRAM_HEADER_BYTES || sam_len > uncomp.saturating_sub(4) {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM FileHeader SAM text length ({sam_len}) in '{}' exceeds uncompressed block bounds ({uncomp} bytes)",
                            path.display()
                        ),
                    ));
                }
            }
        } else {
            let gz_src = (&mut body_take).take(comp as u64);
            let mut gz =
                flate2::read::GzDecoder::new(gz_src).take((MAX_CRAM_HEADER_BYTES as u64) + 5);
            let mut sam_len_buf = [0_u8; 4];
            gz.read_exact(&mut sam_len_buf)?;
            let raw_sam_len = i32::from_le_bytes(sam_len_buf);
            let sam_len = usize::try_from(raw_sam_len)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            if sam_len > MAX_CRAM_HEADER_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM FileHeader SAM text length ({sam_len}) in '{}' exceeds maximum allowed ({MAX_CRAM_HEADER_BYTES} bytes)",
                        path.display()
                    ),
                ));
            }
            let decomp_rem = std::io::copy(&mut gz, &mut std::io::sink())?;
            if decomp_rem > MAX_CRAM_HEADER_BYTES as u64 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM FileHeader decompressed payload in '{}' exceeds maximum allowed ({MAX_CRAM_HEADER_BYTES} bytes)",
                        path.display()
                    ),
                ));
            }
        }
    }
    file.seek(SeekFrom::Start(0))?;
    Ok(())
}

fn read_validated_cram_header(
    path: &Path,
    mut file: std::fs::File,
) -> std::io::Result<(RawCramReader, noodles::sam::Header)> {
    prevalidate_cram_file_header(path, &mut file)?;
    let mut reader = noodles::cram::io::Reader::new(file);
    let header = reader.read_header()?;
    Ok((reader, header))
}

fn acquire_cram_meta(
    path: &Path,
    config: &crate::config::Config,
) -> Result<std::sync::Arc<SharedCramMeta>, crate::types::VardictError> {
    let fasta_path = Path::new(&config.reference);
    let cur_epoch = BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire);
    if let Some(cached_meta) = lookup_cram_meta_tls(path, fasta_path, cur_epoch) {
        return Ok(cached_meta);
    }

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
            let meta_clone = std::sync::Arc::clone(meta);
            if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) == cur_epoch {
                insert_cram_meta_tls(cur_epoch, path, fasta_path, &meta_clone);
            }
            return Ok(meta_clone);
        }
    }

    let index = read_cram_index_with_fallbacks(path)?;
    let indexed_crai = IndexedCrai::from_index(&index);
    let mut spans_by_offset: HashMap<u64, Vec<(usize, usize, usize)>> = HashMap::new();
    for crai_rec in index.iter() {
        if let (Some(ref_id), Some(start)) =
            (crai_rec.reference_sequence_id(), crai_rec.alignment_start())
        {
            let span = crai_rec.alignment_span();
            if span > 0 {
                let s_idx = usize::from(start);
                let e_idx = s_idx.saturating_add(span - 1);
                spans_by_offset
                    .entry(crai_rec.offset())
                    .or_default()
                    .push((ref_id, s_idx, e_idx));
            }
        }
    }
    let file = std::fs::File::open(path)?;
    let file_len = file
        .metadata()
        .ok()
        .filter(|m| m.is_file())
        .map(|m| m.len());
    let (reader, header) = read_validated_cram_header(path, file)?;
    let contig_names: Vec<std::sync::Arc<str>> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| std::sync::Arc::<str>::from(String::from_utf8_lossy(name).as_ref()))
        .collect();
    let id = next_bam_meta_id();
    let meta = std::sync::Arc::new(SharedCramMeta {
        id,
        file_len,
        header,
        contig_names,
        indexed_crai,
        spans_by_offset,
        fasta_repo,
        #[cfg(test)]
        decodes: std::sync::atomic::AtomicUsize::new(0),
    });
    #[cfg(test)]
    if let Some(hook) = TEST_MID_COLD_OPEN_HOOK
        .try_with(|h| h.try_borrow_mut().ok().and_then(|mut g| g.take()))
        .ok()
        .flatten()
    {
        hook();
    }
    if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) == cur_epoch {
        *slot_guard = Some(std::sync::Arc::clone(&meta));
        return_pooled_cram_reader(path, std::sync::Arc::clone(&meta), reader);
        insert_cram_meta_tls(cur_epoch, path, fasta_path, &meta);
    }
    Ok(meta)
}

#[inline]
fn cram_span_exceeds_fasta(
    span_start: usize,
    span_end: usize,
    sq_ln: usize,
    fai_len: usize,
) -> bool {
    span_start > fai_len
        || (span_end > fai_len
            && (fai_len < sq_ln || span_end > fai_len.saturating_add(CRAM_REF_OVERHANG_PAD)))
}

fn validate_cram_span(
    cram_path: &Path,
    meta: &SharedCramMeta,
    reference_sequence_id: usize,
    span_start: usize,
    span_end: usize,
    hint_bounds: Option<(usize, usize, usize)>,
) -> Result<(), crate::types::VardictError> {
    if let Some((hint_ref_id, sq_ln, fai_len)) = hint_bounds {
        if hint_ref_id == reference_sequence_id {
            if cram_span_exceeds_fasta(span_start, span_end, sq_ln, fai_len) {
                let Some((contig_name, _)) = meta
                    .header
                    .reference_sequences()
                    .get_index(reference_sequence_id)
                else {
                    return Err(crate::types::VardictError::Parse(format!(
                        "CRAM '{}' references invalid reference sequence ID {}",
                        cram_path.display(),
                        reference_sequence_id
                    )));
                };
                return Err(crate::types::VardictError::Parse(format!(
                    "CRAM '{}' slice/record on contig '{}' spans {}..={} (header LN={}), exceeding FASTA '{}' contig length {}",
                    cram_path.display(),
                    String::from_utf8_lossy(contig_name),
                    span_start,
                    span_end,
                    sq_ln,
                    meta.fasta_repo.fasta_path.display(),
                    fai_len,
                )));
            }
            return Ok(());
        }
    }

    let Some((contig_name, ref_seq_map)) = meta
        .header
        .reference_sequences()
        .get_index(reference_sequence_id)
    else {
        return Err(crate::types::VardictError::Parse(format!(
            "CRAM '{}' references invalid reference sequence ID {}",
            cram_path.display(),
            reference_sequence_id
        )));
    };
    let sq_ln = usize::from(ref_seq_map.length());
    let (_, seq_opt) = meta.fasta_repo.warm_contig(contig_name);
    let Some(seq) = seq_opt else {
        return Err(crate::types::VardictError::Parse(format!(
            "CRAM '{}' requires reference contig '{}' (span {}..={}), which is missing from FASTA '{}'",
            cram_path.display(),
            String::from_utf8_lossy(contig_name),
            span_start,
            span_end,
            meta.fasta_repo.fasta_path.display(),
        )));
    };
    let fai_len = SharedFastaRepo::unpadded_slice(&seq).len();
    if cram_span_exceeds_fasta(span_start, span_end, sq_ln, fai_len) {
        return Err(crate::types::VardictError::Parse(format!(
            "CRAM '{}' slice/record on contig '{}' spans {}..={} (header LN={}), exceeding FASTA '{}' contig length {}",
            cram_path.display(),
            String::from_utf8_lossy(contig_name),
            span_start,
            span_end,
            sq_ln,
            meta.fasta_repo.fasta_path.display(),
            fai_len,
        )));
    }
    Ok(())
}

const CRAM_CONTAINER_CACHE_SHARDS: usize = 4;
const CRAM_CONTAINER_SLOTS_PER_SHARD: usize = 4;
/// Per-shard decoded CRAM container byte budget (`40 MiB`; `160 MiB` ceiling across `4` shards).
const MAX_SHARD_CACHED_CRAM_BYTES: usize = 40 * 1024 * 1024;

/// RAII guard pinning a chromosome's [`ContigSlot`] in [`SharedFastaRepo`] across
/// all passes of a somatic tile (`bam1 -> bam2 -> combineAnalysis`) so intermediate
/// [`CramFeedInner::drop`] calls cannot evict `tile.chr` mid-tile.
pub(crate) struct CramContigPinGuard {
    inner: Option<(std::sync::Arc<SharedFastaRepo>, std::sync::Arc<ContigSlot>)>,
}

impl Drop for CramContigPinGuard {
    fn drop(&mut self) {
        if let Some((repo, slot)) = self.inner.take() {
            drop(slot);
            repo.evict_unpinned();
        }
    }
}

/// Warms and pins `chr` in [`SharedFastaRepo`] when `bam_spec` contains at least one
/// `.cram` file with `.crai` containers on `chr`; returns a zero-cost empty guard for BAM specs.
pub(crate) fn pin_cram_contig_for_spec(
    bam_spec: &str,
    chr: &str,
    config: &crate::config::Config,
) -> CramContigPinGuard {
    if !crate::runner::bam_spec_has_cram(bam_spec) || chr.is_empty() || config.reference.is_empty()
    {
        return CramContigPinGuard { inner: None };
    }
    let chr_bytes = chr.as_bytes();
    for raw_part in bam_spec.split(['|', ':']).map(str::trim) {
        if raw_part.is_empty() {
            continue;
        }
        let path = Path::new(raw_part);
        if !is_cram_path(path) {
            continue;
        }
        let Ok(meta) = acquire_cram_meta(path, config) else {
            continue;
        };
        let Some(ref_id) = meta.header.reference_sequences().get_index_of(chr_bytes) else {
            continue;
        };
        if meta.indexed_crai.has_ref_containers(ref_id) {
            let slot = meta.fasta_repo.warm_contig_slot(chr_bytes);
            return CramContigPinGuard {
                inner: Some((std::sync::Arc::clone(&meta.fasta_repo), slot)),
            };
        }
    }
    CramContigPinGuard { inner: None }
}

#[derive(Clone, Debug)]
struct CachedCramRecord {
    reference_sequence_id: Option<usize>,
    alignment_start: Option<noodles::core::Position>,
    alignment_end: Option<noodles::core::Position>,
    reference_name: Option<std::sync::Arc<str>>,
    mate_reference_name: Option<std::sync::Arc<str>>,
    sam: SamRecord,
}

#[derive(Clone, Debug)]
struct CachedCramContainer {
    records: Vec<CachedCramRecord>,
    /// `Some((ref_id, max_span_minus_1))` when every record in `records` belongs to
    /// `Some(ref_id)` with monotonically non-decreasing `alignment_start`, enabling
    /// `O(log N)` binary-search start positioning and `O(1)` post-interval early exit
    /// in [`CramFeedInner::pump`].
    sorted_single_ref: Option<(usize, usize)>,
}

impl CachedCramContainer {
    fn from_records(records: Vec<CachedCramRecord>) -> Self {
        let sorted_single_ref = (|| -> Option<(usize, usize)> {
            let first = records.first()?;
            let ref_id = first.reference_sequence_id?;
            let mut prev_start = usize::from(first.alignment_start?);
            let mut max_span = first
                .alignment_end
                .map_or(0, |e| usize::from(e).saturating_sub(prev_start));
            for rec in records.get(1..).unwrap_or(&[]) {
                if rec.reference_sequence_id != Some(ref_id) {
                    return None;
                }
                let start = usize::from(rec.alignment_start?);
                if start < prev_start {
                    return None;
                }
                prev_start = start;
                if let Some(end) = rec.alignment_end {
                    let span = usize::from(end).saturating_sub(start);
                    if span > max_span {
                        max_span = span;
                    }
                }
            }
            Some((ref_id, max_span))
        })();
        Self {
            records,
            sorted_single_ref,
        }
    }

    #[inline]
    fn start_index_and_upper_bound(
        &self,
        reference_sequence_id: usize,
        interval: noodles::core::region::Interval,
    ) -> (usize, Option<usize>) {
        let Some((block_ref_id, max_span)) = self.sorted_single_ref else {
            return (0, None);
        };
        if block_ref_id != reference_sequence_id {
            return (self.records.len(), Some(0));
        }
        let first_idx = match interval.start() {
            Some(q_start) => {
                let min_start = usize::from(q_start).saturating_sub(max_span);
                self.records.partition_point(|r| {
                    r.alignment_start
                        .is_some_and(|s| usize::from(s) < min_start)
                })
            }
            None => 0,
        };
        let end_bound = interval.end().map(usize::from);
        (first_idx, end_bound)
    }
}

impl std::ops::Deref for CachedCramContainer {
    type Target = [CachedCramRecord];

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.records
    }
}

type CachedCramBlock = std::sync::Arc<CachedCramContainer>;

enum InflightCramState {
    Pending,
    Ready(Option<CachedCramBlock>),
    Failed,
}

struct InflightCramDecode {
    state: std::sync::Mutex<InflightCramState>,
    cvar: std::sync::Condvar,
}

impl InflightCramDecode {
    fn new() -> Self {
        Self {
            state: std::sync::Mutex::new(InflightCramState::Pending),
            cvar: std::sync::Condvar::new(),
        }
    }
}

#[derive(Default)]
struct CramContainerCacheShard {
    lru: LruCacheShard<CachedCramBlock, CRAM_CONTAINER_SLOTS_PER_SHARD>,
    entry_bytes: Vec<((u64, u64), usize)>,
    cached_bytes: usize,
    inflight: Vec<((u64, u64), std::sync::Arc<InflightCramDecode>)>,
    #[cfg(test)]
    decodes: usize,
}

impl CramContainerCacheShard {
    fn take_tracked_bytes(&mut self, key: (u64, u64)) -> usize {
        if let Some(pos) = self.entry_bytes.iter().position(|(k, _)| *k == key) {
            let (_, bytes) = self.entry_bytes.swap_remove(pos);
            self.cached_bytes = self.cached_bytes.saturating_sub(bytes);
            bytes
        } else {
            0
        }
    }

    fn insert_with_budget(
        &mut self,
        key: (u64, u64),
        block: CachedCramBlock,
        block_bytes: usize,
        max_shard_bytes: usize,
    ) {
        if let Some(idx) = self.lru.entries.iter().position(|(k, _)| *k == key) {
            self.lru.entries.remove(idx);
            let _ = self.take_tracked_bytes(key);
        }
        if block_bytes > max_shard_bytes {
            return;
        }
        while !self.lru.entries.is_empty()
            && (self.lru.entries.len() >= CRAM_CONTAINER_SLOTS_PER_SHARD
                || self.cached_bytes.saturating_add(block_bytes) > max_shard_bytes)
        {
            let (evicted_key, _) = self.lru.entries.remove(0);
            let _ = self.take_tracked_bytes(evicted_key);
        }
        self.lru.entries.push((key, block));
        self.entry_bytes.push((key, block_bytes));
        self.cached_bytes = self.cached_bytes.saturating_add(block_bytes);
    }

    #[cfg(test)]
    fn remove_key(&mut self, key: (u64, u64)) {
        if let Some(idx) = self.lru.entries.iter().position(|(k, _)| *k == key) {
            self.lru.entries.remove(idx);
        }
        let _ = self.take_tracked_bytes(key);
        let _ = CRAM_LAST_CONTAINER_TLS.try_with(|tls| {
            if let Ok(mut slots) = tls.try_borrow_mut() {
                for slot in slots.iter_mut() {
                    if slot.as_ref().is_some_and(|(_, k, _)| *k == key) {
                        *slot = None;
                    }
                }
            }
        });
    }

    fn clear_production(&mut self) {
        self.lru.clear_production();
        #[cfg(test)]
        {
            self.entry_bytes
                .retain(|((id, _), _)| *id >= TEST_META_ID_BASE);
            self.cached_bytes = self
                .entry_bytes
                .iter()
                .map(|(_, b)| *b)
                .fold(0usize, usize::saturating_add);
        }
        #[cfg(not(test))]
        {
            self.entry_bytes.clear();
            self.cached_bytes = 0;
        }
    }
}

struct InflightCramGuard<'a> {
    shard_mutex: &'a std::sync::Mutex<CramContainerCacheShard>,
    key: (u64, u64),
    entry: std::sync::Arc<InflightCramDecode>,
    completed: bool,
}

impl InflightCramGuard<'_> {
    fn complete(
        mut self,
        result: Option<CachedCramBlock>,
        max_cached_bytes: usize,
        start_epoch: Option<u64>,
    ) {
        self.completed = true;
        let cache_candidate_bytes = result
            .as_ref()
            .and_then(|arc_records| cacheable_cram_container_bytes(arc_records, max_cached_bytes));
        {
            let mut shard = self.shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
            shard
                .inflight
                .retain(|(k, w)| !(*k == self.key && std::sync::Arc::ptr_eq(w, &self.entry)));
            if let (Some(ref arc_records), Some(block_bytes)) = (&result, cache_candidate_bytes) {
                if start_epoch.is_none_or(|ep| {
                    BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) == ep
                }) {
                    shard.insert_with_budget(
                        self.key,
                        std::sync::Arc::clone(arc_records),
                        block_bytes,
                        MAX_SHARD_CACHED_CRAM_BYTES,
                    );
                    if let Some(ep) = start_epoch {
                        if block_bytes <= MAX_SHARD_CACHED_CRAM_BYTES {
                            store_cram_container_tls(
                                ep,
                                self.key,
                                std::sync::Arc::clone(arc_records),
                            );
                        }
                    }
                }
            }
        }
        let mut st = self.entry.state.lock().unwrap_or_else(|e| e.into_inner());
        *st = InflightCramState::Ready(result);
        self.entry.cvar.notify_all();
    }
}

impl Drop for InflightCramGuard<'_> {
    fn drop(&mut self) {
        if !self.completed {
            {
                let mut shard = self.shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
                shard
                    .inflight
                    .retain(|(k, w)| !(*k == self.key && std::sync::Arc::ptr_eq(w, &self.entry)));
            }
            let mut st = self.entry.state.lock().unwrap_or_else(|e| e.into_inner());
            *st = InflightCramState::Failed;
            self.entry.cvar.notify_all();
        }
    }
}

type CramContainerTlsEntry = (u64, (u64, u64), CachedCramBlock);

thread_local! {
    /// Two-slot thread-local MRU cursor (`tumor` + `normal`) so intra-tile SV mate rescue
    /// and intra-batch adjacent-exon queries reuse the worker's active [`CachedCramBlock`]
    /// lock-free without cross-thread shard collisions at `-th 8`.
    static CRAM_LAST_CONTAINER_TLS: RefCell<[Option<CramContainerTlsEntry>; 2]> =
        const { RefCell::new([None, None]) };
    /// Reusable per-thread byte buffer for reading + prevalidating `<= 4 MiB` CRAM containers
    /// in memory without repeated 1-byte `File` syscalls or `SeekWindowReader` allocations.
    static CRAM_IO_SCRATCH_TLS: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[inline]
fn lookup_cram_container_tls(cur_epoch: u64, key: (u64, u64)) -> Option<CachedCramBlock> {
    CRAM_LAST_CONTAINER_TLS
        .try_with(|tls| {
            let mut slots = tls.try_borrow_mut().ok()?;
            let [s0, s1] = &mut *slots;
            if let Some((ep, k, block)) = s0.as_ref() {
                if *ep == cur_epoch && *k == key {
                    return Some(std::sync::Arc::clone(block));
                }
            }
            if let Some((ep, k, block)) = s1.as_ref() {
                if *ep == cur_epoch && *k == key {
                    let hit = std::sync::Arc::clone(block);
                    std::mem::swap(s0, s1);
                    return Some(hit);
                }
            }
            None
        })
        .ok()
        .flatten()
}

#[inline]
fn store_cram_container_tls(cur_epoch: u64, key: (u64, u64), block: CachedCramBlock) {
    let _ = CRAM_LAST_CONTAINER_TLS.try_with(|tls| {
        if let Ok(mut slots) = tls.try_borrow_mut() {
            let [s0, s1] = &mut *slots;
            if s0.as_ref().is_some_and(|(_, k, _)| *k == key) {
                *s0 = Some((cur_epoch, key, block));
            } else {
                *s1 = s0.take();
                *s0 = Some((cur_epoch, key, block));
            }
        }
    });
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

#[cfg(test)]
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

/// Approximates total heap memory owned by a decoded [`CachedCramBlock`].
fn approx_cram_container_bytes(records: &CachedCramBlock) -> usize {
    let mut total = std::mem::size_of::<CachedCramContainer>().saturating_add(
        records
            .records
            .capacity()
            .saturating_mul(std::mem::size_of::<CachedCramRecord>()),
    );
    for rec in records.iter() {
        let sam = &rec.sam;
        let cigar_bytes = sam
            .cigar
            .capacity()
            .saturating_mul(std::mem::size_of::<(u32, CigarKind)>());
        let rec_heap = sam
            .name
            .capacity()
            .saturating_add(cigar_bytes)
            .saturating_add(sam.sequence.capacity())
            .saturating_add(sam.quality_scores.capacity())
            .saturating_add(sam.reference_name.as_ref().map_or(0, String::capacity))
            .saturating_add(sam.mate_reference_name.as_ref().map_or(0, String::capacity))
            .saturating_add(sam.sa.as_ref().map_or(0, String::capacity))
            .saturating_add(sam.mc.as_ref().map_or(0, String::capacity))
            .saturating_add(rec.reference_name.as_deref().map_or(0, str::len))
            .saturating_add(rec.mate_reference_name.as_deref().map_or(0, str::len));
        total = total.saturating_add(rec_heap);
    }
    total
}

/// Returns `Some(approx_bytes)` when `records` fits within both record-count and heap-byte cache bounds.
#[inline]
fn cacheable_cram_container_bytes(
    records: &CachedCramBlock,
    max_cached_bytes: usize,
) -> Option<usize> {
    if records.len() > MAX_CACHED_CRAM_CONTAINER_RECORDS {
        return None;
    }
    let bytes = approx_cram_container_bytes(records);
    (bytes <= max_cached_bytes).then_some(bytes)
}

#[cfg(test)]
#[inline]
fn should_cache_cram_container(records: &CachedCramBlock, max_cached_bytes: usize) -> bool {
    cacheable_cram_container_bytes(records, max_cached_bytes).is_some()
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

#[inline]
fn copy_opt_arc_str(dst: &mut Option<String>, src: &Option<std::sync::Arc<str>>) {
    match src {
        Some(s) => {
            let buf = dst.get_or_insert_with(String::new);
            buf.clear();
            buf.push_str(s);
        }
        None => *dst = None,
    }
}

fn cram_record_to_cached(
    cram_path: &Path,
    record: &noodles::cram::Record<'_>,
    meta: &SharedCramMeta,
    hint_bounds: Option<(usize, usize, usize)>,
) -> Result<CachedCramRecord, crate::types::VardictError> {
    use noodles::sam::alignment::record::data::field::{Tag, Value};
    use noodles::sam::alignment::Record as _;

    let header = &meta.header;
    let contig_names = &meta.contig_names;
    let flags = record.flags()?.bits();
    let reference_sequence_id = record.reference_sequence_id(header).transpose()?;
    let alignment_start = record.alignment_start().transpose()?;

    // Fast pre-filter for unplaced / unaligned CRAM records: records lacking
    // a reference sequence or alignment start never intersect any indexed
    // region query (`cram_record_intersects` returns `false`), so skip
    // materializing heap buffers (`name`, `cigar`, `sequence`, `quality_scores`,
    // `sa`, `mc`) while still validating raw fields.
    if reference_sequence_id.is_none() || alignment_start.is_none() {
        let _ = record.mapping_quality().transpose()?;
        for result in record.cigar().iter() {
            let op = result?;
            let _ = checked_cigar_op_len(op.len())?;
        }
        let _ = record.mate_reference_sequence_id(header).transpose()?;
        let _ = record.mate_alignment_start().transpose()?;
        let _ = record.template_length()?;
        for (idx, b) in record.sequence().iter().enumerate() {
            if normalize_cram_base(b).is_none() {
                let name = record
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_default();
                return Err(crate::types::VardictError::Parse(format!(
                    "CRAM {} read {}: non-ASCII base 0x{:02x} at {}",
                    cram_path.display(),
                    name,
                    b,
                    idx + 1,
                )));
            }
        }
        for result in record.quality_scores().iter() {
            let _ = result?;
        }
        let raw_data = record.data();
        if !raw_data.is_empty() {
            for result in raw_data.iter() {
                let (_, val) = result?;
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
            }
        }
        return Ok(CachedCramRecord {
            reference_sequence_id,
            alignment_start,
            alignment_end: None,
            reference_name: None,
            mate_reference_name: None,
            sam: SamRecord::new(),
        });
    }
    let mut name = record
        .name()
        .map(|n| {
            let bytes: &[u8] = n.as_ref();
            bytes.to_vec()
        })
        .unwrap_or_default();
    if name.capacity() > name.len() {
        name.shrink_to_fit();
    }
    let reference_name = reference_sequence_id
        .and_then(|id| contig_names.get(id))
        .cloned();
    let raw_mapq = record.mapping_quality().transpose()?;
    let mapping_quality = if (flags & 0x4) != 0 {
        None
    } else {
        Some(raw_mapq.map(|m| m.get()).unwrap_or(255))
    };

    // In `noodles-cram 0.100.0`, `Cigar::len()` calls `self.iter().count()`,
    // boxing and walking all CRAM features just to count ops. Buffer up to 8 ops
    // on the stack so common Illumina reads allocate once at exact capacity without
    // per-record `shrink_to_fit` reallocations.
    let mut cigar_inline = [(0_u32, CigarKind::Match); 8];
    let mut cigar_inline_len = 0usize;
    let mut cigar_spill: Vec<(u32, CigarKind)> = Vec::new();
    let mut alignment_span = 0usize;
    for result in record.cigar().iter() {
        let op = result?;
        let op_len = checked_cigar_op_len(op.len())?;
        if op.kind().consumes_reference() {
            alignment_span = alignment_span.checked_add(op.len()).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "CRAM CIGAR alignment_span overflow",
                )
            })?;
        }
        if cigar_spill.is_empty() {
            if let Some(slot) = cigar_inline.get_mut(cigar_inline_len) {
                *slot = (op_len, op.kind());
                cigar_inline_len += 1;
            } else {
                cigar_spill.reserve_exact(16);
                cigar_spill.extend_from_slice(&cigar_inline);
                cigar_spill.push((op_len, op.kind()));
            }
        } else {
            cigar_spill.push((op_len, op.kind()));
        }
    }
    let cigar = if cigar_spill.is_empty() {
        cigar_inline.get(..cigar_inline_len).unwrap_or(&[]).to_vec()
    } else {
        if cigar_spill.capacity() > cigar_spill.len() {
            cigar_spill.shrink_to_fit();
        }
        cigar_spill
    };
    let alignment_end = match alignment_start {
        Some(start) => Some(checked_alignment_end(start, alignment_span)?),
        None => None,
    };

    if (flags & 0x4) == 0 {
        if let (Some(ref_id), Some(start), Some(end)) =
            (reference_sequence_id, alignment_start, alignment_end)
        {
            validate_cram_span(
                cram_path,
                meta,
                ref_id,
                usize::from(start),
                usize::from(end),
                hint_bounds,
            )?;
        }
    }

    let mate_reference_sequence_id = record.mate_reference_sequence_id(header).transpose()?;
    let mate_reference_name = mate_reference_sequence_id
        .and_then(|id| contig_names.get(id))
        .cloned();
    let mate_alignment_start = record.mate_alignment_start().transpose()?;
    let template_length = record.template_length()?;

    let raw_seq = record.sequence();
    let mut sequence = Vec::with_capacity(raw_seq.len());
    for (idx, b) in raw_seq.iter().enumerate() {
        let Some(norm) = normalize_cram_base(b) else {
            return Err(crate::types::VardictError::Parse(format!(
                "CRAM {} read {}: non-ASCII base 0x{:02x} at {}",
                cram_path.display(),
                String::from_utf8_lossy(&name),
                b,
                idx + 1,
            )));
        };
        sequence.push(norm);
    }
    if sequence.capacity() > sequence.len() {
        sequence.shrink_to_fit();
    }

    let raw_quals = record.quality_scores();
    let mut quality_scores = Vec::with_capacity(raw_quals.len());
    for result in raw_quals.iter() {
        let score = result?;
        quality_scores.push(score.saturating_add(33));
    }
    if quality_scores.capacity() > quality_scores.len() {
        quality_scores.shrink_to_fit();
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
    if !raw_data.is_empty() {
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
                nm_lower_val = val
                    .as_int()
                    .map(|x| x.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32);
            } else if tag == TAG_NM_UPPER {
                seen_nm_upper = true;
                nm_upper_val = val
                    .as_int()
                    .map(|x| x.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32);
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
                mq = val.as_int().and_then(|x| u8::try_from(x).ok());
            }
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
        reference_name: None,
        alignment_start: alignment_start.map(crate::types::Pos1::from),
        alignment_end: alignment_end.map(crate::types::Pos1::from),
        mapping_quality,
        cigar,
        sequence,
        quality_scores,
        mate_reference_name: None,
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
        let mut check_sam = sam.clone();
        check_sam.reference_name = reference_name.as_deref().map(str::to_owned);
        check_sam.mate_reference_name = mate_reference_name.as_deref().map(str::to_owned);
        assert_eq!(reference_sequence_id, rec_buf.reference_sequence_id());
        assert_eq!(alignment_start, rec_buf.alignment_start());
        assert_eq!(alignment_end, rec_buf.alignment_end());
        assert_eq!(check_sam, expected_sam);
    }

    Ok(CachedCramRecord {
        reference_sequence_id,
        alignment_start,
        alignment_end,
        reference_name,
        mate_reference_name,
        sam,
    })
}

fn read_cram_itf8<R: std::io::Read + ?Sized>(reader: &mut R) -> std::io::Result<i32> {
    let mut b0_buf = [0_u8; 1];
    reader.read_exact(&mut b0_buf)?;
    let [b0_raw] = b0_buf;
    let b0 = u32::from(b0_raw);
    let n = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        let mut b1 = [0_u8; 1];
        reader.read_exact(&mut b1)?;
        let [b1_raw] = b1;
        ((b0 & 0x7f) << 8) | u32::from(b1_raw)
    } else if b0 & 0x20 == 0 {
        let mut b12 = [0_u8; 2];
        reader.read_exact(&mut b12)?;
        ((b0 & 0x3f) << 16) | u32::from(u16::from_be_bytes(b12))
    } else if b0 & 0x10 == 0 {
        let mut b123 = [0_u8; 3];
        reader.read_exact(&mut b123)?;
        let [b1, b2, b3] = b123;
        ((b0 & 0x1f) << 24) | (u32::from(b1) << 16) | (u32::from(b2) << 8) | u32::from(b3)
    } else {
        let mut b1234 = [0_u8; 4];
        reader.read_exact(&mut b1234)?;
        let b1_4 = u32::from_be_bytes(b1234);
        ((b0 & 0x0f) << 28) | ((b1_4 & 0xffff_ff0f) >> 4) | (b1_4 & 0x0f)
    };
    Ok(n as i32)
}

fn read_cram_ltf8<R: std::io::Read + ?Sized>(reader: &mut R) -> std::io::Result<i64> {
    let mut b0_buf = [0_u8; 1];
    reader.read_exact(&mut b0_buf)?;
    let [b0_raw] = b0_buf;
    let b0 = u64::from(b0_raw);
    let n = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        let mut b1 = [0_u8; 1];
        reader.read_exact(&mut b1)?;
        let [c1] = b1;
        ((b0 & 0x7f) << 8) | u64::from(c1)
    } else if b0 & 0x20 == 0 {
        let mut b12 = [0_u8; 2];
        reader.read_exact(&mut b12)?;
        ((b0 & 0x3f) << 16) | u64::from(u16::from_be_bytes(b12))
    } else if b0 & 0x10 == 0 {
        let mut b123 = [0_u8; 3];
        reader.read_exact(&mut b123)?;
        let [c1, c2, c3] = b123;
        ((b0 & 0x1f) << 24) | (u64::from(c1) << 16) | (u64::from(c2) << 8) | u64::from(c3)
    } else if b0 & 0x08 == 0 {
        let mut b14 = [0_u8; 4];
        reader.read_exact(&mut b14)?;
        ((b0 & 0x0f) << 32) | u64::from(u32::from_be_bytes(b14))
    } else if b0 & 0x04 == 0 {
        let mut b15 = [0_u8; 5];
        reader.read_exact(&mut b15)?;
        let [c1, c2, c3, c4, c5] = b15;
        ((b0 & 0x07) << 40) | u64::from_be_bytes([0, 0, 0, c1, c2, c3, c4, c5])
    } else if b0 & 0x02 == 0 {
        let mut b16 = [0_u8; 6];
        reader.read_exact(&mut b16)?;
        let [c1, c2, c3, c4, c5, c6] = b16;
        ((b0 & 0x03) << 48) | u64::from_be_bytes([0, 0, c1, c2, c3, c4, c5, c6])
    } else if b0 & 0x01 == 0 {
        let mut b17 = [0_u8; 7];
        reader.read_exact(&mut b17)?;
        let [c1, c2, c3, c4, c5, c6, c7] = b17;
        u64::from_be_bytes([0, c1, c2, c3, c4, c5, c6, c7])
    } else {
        let mut b18 = [0_u8; 8];
        reader.read_exact(&mut b18)?;
        u64::from_be_bytes(b18)
    };
    Ok(n as i64)
}

fn read_cram_uint7<R: std::io::Read + ?Sized>(reader: &mut R) -> std::io::Result<u32> {
    let mut n = 0_u32;
    for _ in 0..5 {
        let mut b_buf = [0_u8; 1];
        reader.read_exact(&mut b_buf)?;
        let [b_raw] = b_buf;
        let b = u32::from(b_raw);
        n = n.wrapping_shl(7) | (b & 0x7f);
        if b & 0x80 == 0 {
            return Ok(n);
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        "CRAM uint7 VLQ exceeds 5 bytes",
    ))
}

/// Windowed `Read + Seek` adapter that retains its in-memory read window across
/// `SeekFrom::Start` jumps (unlike `std::io::BufReader::seek`, which discards its
/// buffer on every `SeekFrom::Start`), collapsing ~250k per-block header/codec
/// seeks and 1-byte VLQ reads in `prevalidate_cram_container_slices` to memory hits.
struct SeekWindowReader<R> {
    inner: R,
    buf: Vec<u8>,
    buf_start: u64,
    buf_len: usize,
    pos: u64,
}

impl<R> SeekWindowReader<R> {
    fn with_capacity(capacity: usize, inner: R) -> Self {
        Self {
            inner,
            buf: vec![0_u8; capacity.max(1)],
            buf_start: 0,
            buf_len: 0,
            pos: 0,
        }
    }
}

impl<R: std::io::Read + std::io::Seek> std::io::Read for SeekWindowReader<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        let buf_end = self.buf_start.saturating_add(self.buf_len as u64);
        if self.pos < self.buf_start || self.pos >= buf_end {
            self.inner.seek(std::io::SeekFrom::Start(self.pos))?;
            self.buf_start = self.pos;
            self.buf_len = 0;
            loop {
                match self.inner.read(&mut self.buf) {
                    Ok(n) => {
                        self.buf_len = n;
                        break;
                    }
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                    Err(e) => return Err(e),
                }
            }
            if self.buf_len == 0 {
                return Ok(0);
            }
        }
        let rel = usize::try_from(self.pos - self.buf_start).unwrap_or(self.buf_len);
        let avail = self.buf.get(rel..self.buf_len).unwrap_or(&[]);
        let n = out.len().min(avail.len());
        if let (Some(dst), Some(src)) = (out.get_mut(..n), avail.get(..n)) {
            dst.copy_from_slice(src);
        }
        self.pos = self.pos.saturating_add(n as u64);
        Ok(n)
    }
}

impl<R: std::io::Seek> std::io::Seek for SeekWindowReader<R> {
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        match pos {
            std::io::SeekFrom::Start(p) => {
                self.pos = p;
                Ok(p)
            }
            std::io::SeekFrom::Current(rel) => {
                let next = if rel >= 0 {
                    self.pos.checked_add(rel as u64)
                } else {
                    self.pos.checked_sub(rel.unsigned_abs())
                }
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "invalid seek to a negative or overflowing position",
                    )
                })?;
                self.pos = next;
                Ok(next)
            }
            std::io::SeekFrom::End(rel) => {
                let next = self.inner.seek(std::io::SeekFrom::End(rel))?;
                self.pos = next;
                Ok(next)
            }
        }
    }
}

fn validate_cram_entropy_stream_at(
    path: &Path,
    file: &mut (impl std::io::Read + std::io::Seek),
    start_offset: u64,
    comp_len: usize,
    block_uncomp: usize,
    is_aac: bool,
    depth: u8,
) -> std::io::Result<usize> {
    use std::io::{Read as _, SeekFrom};

    let end_offset = start_offset.checked_add(comp_len as u64).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "CRAM entropy stream offset overflows u64",
        )
    })?;
    file.seek(SeekFrom::Start(start_offset))?;

    let (flags, uncomp, pos_after_size) = {
        let mut r = (&mut *file).take(comp_len as u64);
        let mut f_buf = [0_u8; 1];
        r.read_exact(&mut f_buf)?;
        let [flags] = f_buf;
        let mut u = block_uncomp;
        if (flags & 0x10) == 0 {
            u = read_cram_uint7(&mut r)? as usize;
            if u > MAX_CRAM_CONTAINER_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM codec inner uncompressed size ({u}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                        path.display()
                    ),
                ));
            }
        }
        let consumed = (comp_len as u64).saturating_sub(r.limit());
        (flags, u, start_offset.saturating_add(consumed))
    };

    if (flags & 0x08) != 0 {
        if depth >= 2 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM striped codec nesting depth in '{}' exceeds limit",
                    path.display()
                ),
            ));
        }
        let (n_chunks, sub_comp_sizes, mut cur_pos) = {
            let rem = end_offset.saturating_sub(pos_after_size);
            let mut r = (&mut *file).take(rem);
            let mut cc_buf = [0_u8; 1];
            r.read_exact(&mut cc_buf)?;
            let [cc] = cc_buf;
            let n = usize::from(cc);
            if n == 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM striped codec in '{}' has zero chunk_count",
                        path.display()
                    ),
                ));
            }
            let mut sizes = Vec::with_capacity(n);
            for _ in 0..n {
                sizes.push(read_cram_uint7(&mut r)? as usize);
            }
            let consumed = rem.saturating_sub(r.limit());
            (n, sizes, pos_after_size.saturating_add(consumed))
        };

        let mut sum_sub_uncomp = 0_usize;
        for (idx, sub_comp) in sub_comp_sizes.into_iter().enumerate() {
            let next_pos = cur_pos
                .checked_add(sub_comp as u64)
                .filter(|&p| p <= end_offset)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM striped sub-chunk compressed size ({sub_comp}) in '{}' exceeds block bounds",
                            path.display()
                        ),
                    )
                })?;
            let (q, rem_u) = (uncomp / n_chunks, uncomp % n_chunks);
            let sub_uncomp = if rem_u > idx { q + 1 } else { q };
            let child_peak = validate_cram_entropy_stream_at(
                path,
                file,
                cur_pos,
                sub_comp,
                sub_uncomp,
                is_aac,
                depth + 1,
            )?;
            sum_sub_uncomp = sum_sub_uncomp.saturating_add(child_peak);
            if sum_sub_uncomp > MAX_CRAM_CONTAINER_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM striped sub-chunks cumulative uncompressed size ({sum_sub_uncomp}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                        path.display()
                    ),
                ));
            }
            cur_pos = next_pos;
        }
        return Ok(uncomp.max(sum_sub_uncomp));
    }

    let mut peak = uncomp;
    let mut cur_pos = pos_after_size;

    if (flags & 0x80) != 0 {
        file.seek(SeekFrom::Start(cur_pos))?;
        let sym_count = {
            let rem = end_offset.saturating_sub(cur_pos);
            let mut r = (&mut *file).take(rem);
            let mut sc_buf = [0_u8; 1];
            r.read_exact(&mut sc_buf)?;
            let [sc] = sc_buf;
            if sc == 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM bit-pack codec in '{}' has zero symbol_count",
                        path.display()
                    ),
                ));
            }
            sc
        };
        let after_map = cur_pos
            .checked_add(1 + u64::from(sym_count))
            .filter(|&p| p <= end_offset)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM bit-pack mapping table in '{}' exceeds block bounds",
                        path.display()
                    ),
                )
            })?;
        file.seek(SeekFrom::Start(after_map))?;
        let (pack_len, next_pos) = {
            let rem = end_offset.saturating_sub(after_map);
            let mut r = (&mut *file).take(rem);
            let len = read_cram_uint7(&mut r)? as usize;
            let consumed = rem.saturating_sub(r.limit());
            (len, after_map.saturating_add(consumed))
        };
        if pack_len > MAX_CRAM_CONTAINER_BYTES {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM bit-pack inner length ({pack_len}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                    path.display()
                ),
            ));
        }
        peak = peak.max(pack_len);
        cur_pos = next_pos;
    }

    if !is_aac && (flags & 0x40) != 0 {
        file.seek(SeekFrom::Start(cur_pos))?;
        let rem = end_offset.saturating_sub(cur_pos);
        let mut r = (&mut *file).take(rem);
        let n = read_cram_uint7(&mut r)?;
        let context_size = (n >> 1) as usize;
        let rle_len = read_cram_uint7(&mut r)? as usize;
        if context_size > MAX_CRAM_CONTAINER_BYTES || rle_len > MAX_CRAM_CONTAINER_BYTES {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM RLE inner size (context={context_size}, len={rle_len}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                    path.display()
                ),
            ));
        }
        peak = peak.max(rle_len).max(context_size);
    }

    Ok(peak)
}

fn validate_cram_block_and_codec(
    path: &Path,
    file: &mut (impl std::io::Read + std::io::Seek),
    block_offset: u64,
    limit_end_offset: u64,
    expected_content_type: u8,
    max_uncompressed_bytes: usize,
    max_name_count: usize,
) -> std::io::Result<(u64, usize)> {
    use std::io::{Read as _, SeekFrom};

    if block_offset >= limit_end_offset {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "CRAM block offset ({block_offset}) in '{}' exceeds slice/container end ({limit_end_offset})",
                path.display()
            ),
        ));
    }
    file.seek(SeekFrom::Start(block_offset))?;
    let max_hdr_window = limit_end_offset.saturating_sub(block_offset);
    let (raw_method, content_type, compressed_size, uncompressed_size, payload_offset) = {
        let mut bounded_hdr = (&mut *file).take(max_hdr_window);
        let mut prefix = [0_u8; 2];
        bounded_hdr.read_exact(&mut prefix)?;
        let [raw_m, c_type] = prefix;
        let _content_id = read_cram_itf8(&mut bounded_hdr)?;
        let comp = usize::try_from(read_cram_itf8(&mut bounded_hdr)?)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let uncomp = usize::try_from(read_cram_itf8(&mut bounded_hdr)?)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let consumed = max_hdr_window.saturating_sub(bounded_hdr.limit());
        (
            raw_m,
            c_type,
            comp,
            uncomp,
            block_offset.saturating_add(consumed),
        )
    };

    if raw_method > 8 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "invalid CRAM block compression method ({raw_method}) in '{}'",
                path.display()
            ),
        ));
    }
    if content_type != expected_content_type {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "invalid CRAM block content type in '{}': expected {expected_content_type}, got {content_type}",
                path.display()
            ),
        ));
    }

    let next_block_offset = payload_offset
        .checked_add(compressed_size as u64)
        .and_then(|end_data| end_data.checked_add(4))
        .filter(|&nb| nb <= limit_end_offset)
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM block compressed size ({compressed_size}) + CRC32 in '{}' exceeds bounds",
                    path.display()
                ),
            )
        })?;

    if uncompressed_size > max_uncompressed_bytes {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "CRAM block uncompressed size ({uncompressed_size}) in '{}' exceeds maximum allowed ({max_uncompressed_bytes})",
                path.display()
            ),
        ));
    }

    if uncompressed_size == 0 {
        if raw_method == 0 && compressed_size != 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM uncompressed raw block in '{}' has compressed_size ({compressed_size}) != 0",
                    path.display()
                ),
            ));
        }
        return Ok((next_block_offset, 0));
    }

    let mut effective_uncomp = uncompressed_size;
    match raw_method {
        0 => {
            if compressed_size != uncompressed_size {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM raw block in '{}' has compressed_size ({compressed_size}) != uncompressed_size ({uncompressed_size})",
                        path.display()
                    ),
                ));
            }
        }
        1..=3 => {}
        4 => {
            file.seek(SeekFrom::Start(payload_offset))?;
            let mut r = (&mut *file).take(compressed_size as u64);
            let mut rans_hdr = [0_u8; 9];
            r.read_exact(&mut rans_hdr)?;
            let [order, _c0, _c1, _c2, _c3, u0, u1, u2, u3] = rans_hdr;
            if order > 1 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "invalid CRAM rANS 4x8 order ({order}) in '{}'",
                        path.display()
                    ),
                ));
            }
            let inner_uncomp = u32::from_le_bytes([u0, u1, u2, u3]) as usize;
            if inner_uncomp > MAX_CRAM_CONTAINER_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM rANS 4x8 inner uncompressed size ({inner_uncomp}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                        path.display()
                    ),
                ));
            }
            effective_uncomp = effective_uncomp.max(inner_uncomp);
        }
        5 | 6 => {
            let inner_peak = validate_cram_entropy_stream_at(
                path,
                file,
                payload_offset,
                compressed_size,
                uncompressed_size,
                raw_method == 6,
                0,
            )?;
            effective_uncomp = effective_uncomp.max(inner_peak);
        }
        7 => {
            file.seek(SeekFrom::Start(payload_offset))?;
            let mut r = (&mut *file).take(compressed_size as u64);
            let inner_uncomp = read_cram_uint7(&mut r)? as usize;
            if inner_uncomp > MAX_CRAM_CONTAINER_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM fqzcomp inner uncompressed size ({inner_uncomp}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                        path.display()
                    ),
                ));
            }
            effective_uncomp = effective_uncomp.max(inner_uncomp);
        }
        8 => {
            const MAX_CRAM_NAME_TOKENIZER_NAMES: usize = 131_072;
            const MAX_CRAM_NAME_TOKENIZER_POSITIONS: usize = 128;

            file.seek(SeekFrom::Start(payload_offset))?;
            let mut r = (&mut *file).take(compressed_size as u64);
            let mut nt_hdr = [0_u8; 9];
            r.read_exact(&mut nt_hdr)?;
            let [u0, u1, u2, u3, n0, n1, n2, n3, sub_method] = nt_hdr;
            let inner_uncomp = u32::from_le_bytes([u0, u1, u2, u3]) as usize;
            let name_count = u32::from_le_bytes([n0, n1, n2, n3]) as usize;
            if inner_uncomp > MAX_CRAM_CONTAINER_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM name_tokenizer uncompressed size ({inner_uncomp}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                        path.display()
                    ),
                ));
            }
            let allowed_names = max_name_count
                .min(MAX_CRAM_NAME_TOKENIZER_NAMES)
                .min(inner_uncomp.saturating_add(1));
            if name_count == 0 || name_count > allowed_names {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM name_tokenizer name_count ({name_count}) in '{}' is zero or exceeds maximum allowed ({allowed_names})",
                        path.display()
                    ),
                ));
            }
            let payload_end = payload_offset.saturating_add(compressed_size as u64);
            let mut cur_pos = payload_offset.saturating_add(9);
            let mut stream_ops = 0_usize;
            let mut tok_new_count = 0_usize;
            let mut tok_uncomp_sum = 0_usize;
            let mut max_stream_uncomp = 0_usize;
            while cur_pos < payload_end {
                stream_ops += 1;
                if stream_ops > 4096 {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM name_tokenizer token stream count in '{}' exceeds maximum allowed (4096)",
                            path.display()
                        ),
                    ));
                }
                file.seek(SeekFrom::Start(cur_pos))?;
                let rem = payload_end.saturating_sub(cur_pos);
                let (ttype, dup_info_opt, sub_comp_opt, consumed) = {
                    let mut tr = (&mut *file).take(rem);
                    let mut tt_buf = [0_u8; 1];
                    tr.read_exact(&mut tt_buf)?;
                    let [tt] = tt_buf;
                    if (tt & 0x40) != 0 {
                        let mut dup_buf = [0_u8; 2];
                        tr.read_exact(&mut dup_buf)?;
                        let [dup_pos, dup_ty_raw] = dup_buf;
                        (tt, Some((dup_pos, dup_ty_raw)), None, 3_u64)
                    } else {
                        let sc = read_cram_uint7(&mut tr)? as usize;
                        let used = rem.saturating_sub(tr.limit());
                        (tt, None, Some(sc), used)
                    }
                };
                let raw_ty = ttype & 0x3f;
                if raw_ty > 9 {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM name_tokenizer stream type ({raw_ty}) in '{}' exceeds maximum allowed (9)",
                            path.display()
                        ),
                    ));
                }
                if (ttype & 0x80) != 0 {
                    tok_new_count += 1;
                    if tok_new_count > MAX_CRAM_NAME_TOKENIZER_POSITIONS {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!(
                                "CRAM name_tokenizer token position count ({tok_new_count}) in '{}' exceeds maximum allowed ({MAX_CRAM_NAME_TOKENIZER_POSITIONS})",
                                path.display()
                            ),
                        ));
                    }
                    if raw_ty != 0 {
                        max_stream_uncomp = max_stream_uncomp.max(name_count);
                        tok_uncomp_sum = tok_uncomp_sum.saturating_add(name_count);
                    }
                } else if tok_new_count == 0 {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM name_tokenizer first token stream in '{}' lacks tok_new flag",
                            path.display()
                        ),
                    ));
                }
                if let Some((dup_pos, dup_ty_raw)) = dup_info_opt {
                    let dup_ty = dup_ty_raw & 0x3f;
                    if usize::from(dup_pos) >= tok_new_count || dup_ty > 9 || dup_ty == 4 {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!(
                                "invalid CRAM name_tokenizer tok_dup reference (pos={dup_pos}, ty={dup_ty}, tok_new_count={tok_new_count}) in '{}'",
                                path.display()
                            ),
                        ));
                    }
                    tok_uncomp_sum = tok_uncomp_sum.saturating_add(max_stream_uncomp);
                    cur_pos = cur_pos.saturating_add(consumed);
                } else if let Some(sub_comp) = sub_comp_opt {
                    let stream_start = cur_pos.saturating_add(consumed);
                    let next_pos = stream_start
                        .checked_add(sub_comp as u64)
                        .filter(|&p| p <= payload_end)
                        .ok_or_else(|| {
                            std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                format!(
                                    "CRAM name_tokenizer stream compressed size ({sub_comp}) in '{}' exceeds block bounds",
                                    path.display()
                                ),
                            )
                        })?;
                    let stream_peak = validate_cram_entropy_stream_at(
                        path,
                        file,
                        stream_start,
                        sub_comp,
                        0,
                        sub_method != 0,
                        0,
                    )?;
                    max_stream_uncomp = max_stream_uncomp.max(stream_peak);
                    tok_uncomp_sum = tok_uncomp_sum.saturating_add(stream_peak);
                    cur_pos = next_pos;
                }
                if tok_uncomp_sum > MAX_CRAM_CONTAINER_BYTES {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM name_tokenizer cumulative stream size ({tok_uncomp_sum}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                            path.display()
                        ),
                    ));
                }
            }
            if tok_new_count == 0 || tok_uncomp_sum < name_count {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM name_tokenizer decoded stream total ({tok_uncomp_sum}, tok_new_count={tok_new_count}) in '{}' fails minimum name_count ({name_count}); name_count exceeds stream capacity",
                        path.display()
                    ),
                ));
            }
            effective_uncomp = effective_uncomp.max(inner_uncomp).max(tok_uncomp_sum);
        }
        _ => {}
    }

    Ok((next_block_offset, effective_uncomp))
}

/// Pre-validates CRAM container landmarks, `CompressionHeader`, `SliceHeader` fields, and
/// every slice's `CoreData`/`ExternalData` outer + inner codec allocation headers directly on
/// `file` before `noodles-cram` slices `container.src` or allocates decoded buffers (`SEC-06`).
fn prevalidate_cram_container_slices(
    path: &Path,
    file: &mut (impl std::io::Read + std::io::Seek),
    body_offset: u64,
    bytes_read: usize,
    landmarks: &[usize],
    meta: Option<&SharedCramMeta>,
    hint_bounds: Option<(usize, usize, usize)>,
) -> std::io::Result<()> {
    use std::io::{Read as _, SeekFrom};
    const MAX_SLICE_HEADER_UNCOMPRESSED_BYTES: usize = 1024 * 1024;
    const COMPRESSION_HEADER_CONTENT_TYPE: u8 = 1;
    const SLICE_HEADER_CONTENT_TYPE: u8 = 2;
    const EXTERNAL_DATA_CONTENT_TYPE: u8 = 4;
    const CORE_DATA_CONTENT_TYPE: u8 = 5;
    const MAX_SLICE_BLOCKS: usize = 4096;

    for &lm in landmarks {
        if lm > bytes_read {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "malformed CRAM container in '{}': landmark ({lm}) exceeds container length ({bytes_read})",
                    path.display()
                ),
            ));
        }
    }
    if landmarks.windows(2).any(|w| matches!(w, [a, b] if *a > *b)) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "malformed CRAM container in '{}': non-monotonic slice landmarks",
                path.display()
            ),
        ));
    }

    // Validate CompressionHeader block and inner codec bounds before `container.compression_header()`
    let comp_hdr_end = landmarks.first().copied().unwrap_or(bytes_read);
    if comp_hdr_end > 0 {
        let comp_hdr_limit = body_offset
            .checked_add(comp_hdr_end as u64)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "CRAM compression header end overflows u64",
                )
            })?;
        let _ = validate_cram_block_and_codec(
            path,
            file,
            body_offset,
            comp_hdr_limit,
            COMPRESSION_HEADER_CONTENT_TYPE,
            MAX_CRAM_CONTAINER_BYTES,
            0,
        )?;
    }

    let mut cumulative_records = 0_usize;
    for (idx, &start) in landmarks.iter().enumerate() {
        let end = landmarks.get(idx + 1).copied().unwrap_or(bytes_read);
        let slice_span = end.saturating_sub(start);
        let slice_offset = body_offset.checked_add(start as u64).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "CRAM slice offset overflows u64",
            )
        })?;
        let slice_end_offset = slice_offset.checked_add(slice_span as u64).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "CRAM slice end offset overflows u64",
            )
        })?;
        file.seek(SeekFrom::Start(slice_offset))?;

        let (mut method, content_type, compressed_size, uncompressed_size, sh_payload_offset) = {
            let mut hdr_bounded = (&mut *file).take(slice_span as u64);
            let mut blk_hdr = [0_u8; 2];
            hdr_bounded.read_exact(&mut blk_hdr)?;
            let [m, ct] = blk_hdr;
            let _content_id = read_cram_itf8(&mut hdr_bounded)?;
            let comp = usize::try_from(read_cram_itf8(&mut hdr_bounded)?)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let uncomp = usize::try_from(read_cram_itf8(&mut hdr_bounded)?)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            let consumed = (slice_span as u64).saturating_sub(hdr_bounded.limit());
            (m, ct, comp, uncomp, slice_offset.saturating_add(consumed))
        };
        if content_type != SLICE_HEADER_CONTENT_TYPE {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "invalid CRAM slice header block content type in '{}': expected 2, got {content_type}",
                    path.display()
                ),
            ));
        }
        let first_data_blk_offset = sh_payload_offset
            .checked_add(compressed_size as u64)
            .and_then(|p| p.checked_add(4))
            .filter(|&nb| nb <= slice_end_offset)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM slice header block size in '{}' exceeds bounds (compressed={compressed_size}, uncompressed={uncompressed_size})",
                        path.display()
                    ),
                )
            })?;
        if uncompressed_size > MAX_SLICE_HEADER_UNCOMPRESSED_BYTES {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM slice header block size in '{}' exceeds bounds (compressed={compressed_size}, uncompressed={uncompressed_size})",
                    path.display()
                ),
            ));
        }
        if uncompressed_size == 0 {
            method = 0;
        }
        file.seek(SeekFrom::Start(sh_payload_offset))?;
        type RawSliceHdrFields = (i32, i32, i32, i32, i64, i32, i32);
        let read_slice_hdr_fields =
            |r: &mut dyn std::io::Read| -> std::io::Result<RawSliceHdrFields> {
                let ref_seq_id = read_cram_itf8(r)?;
                let align_start = read_cram_itf8(r)?;
                let align_span = read_cram_itf8(r)?;
                let rc = read_cram_itf8(r)?;
                let rcounter = read_cram_ltf8(r)?;
                let bcount = read_cram_itf8(r)?;
                let cids_len = read_cram_itf8(r)?;
                Ok((
                    ref_seq_id,
                    align_start,
                    align_span,
                    rc,
                    rcounter,
                    bcount,
                    cids_len,
                ))
            };
        let (
            raw_ref_seq_id,
            raw_align_start,
            raw_align_span,
            raw_record_count,
            raw_record_counter,
            raw_block_count,
            raw_content_ids_len,
        ) = match method {
            0 => {
                if uncompressed_size > 0 && compressed_size != uncompressed_size {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM raw slice header block in '{}' has compressed_size ({compressed_size}) != uncompressed_size ({uncompressed_size})",
                            path.display()
                        ),
                    ));
                }
                let mut bounded = (&mut *file).take(compressed_size as u64);
                read_slice_hdr_fields(&mut bounded)?
            }
            1 => {
                let bounded = (&mut *file).take(compressed_size as u64);
                let mut gz = flate2::read::GzDecoder::new(bounded)
                    .take(MAX_SLICE_HEADER_UNCOMPRESSED_BYTES as u64);
                read_slice_hdr_fields(&mut gz)?
            }
            _ => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "unsupported compression method ({method}) on CRAM slice header block in '{}'",
                        path.display()
                    ),
                ));
            }
        };
        if raw_ref_seq_id < -2 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "invalid CRAM slice reference_sequence_id ({raw_ref_seq_id}) in '{}'",
                    path.display()
                ),
            ));
        }
        if raw_ref_seq_id >= 0 {
            if raw_align_start < 1 || raw_align_span < 1 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "invalid CRAM slice alignment bounds (start={raw_align_start}, span={raw_align_span}) in '{}'",
                        path.display()
                    ),
                ));
            }
            let (Ok(ref_id), Ok(start_pos), Ok(span_len)) = (
                usize::try_from(raw_ref_seq_id),
                usize::try_from(raw_align_start),
                usize::try_from(raw_align_span),
            ) else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("invalid CRAM slice coordinates in '{}'", path.display()),
                ));
            };
            let Some(end_pos) = start_pos.checked_add(span_len - 1) else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM slice alignment span overflow (start={start_pos}, span={span_len}) in '{}'",
                        path.display()
                    ),
                ));
            };
            if let Some(m) = meta {
                validate_cram_span(path, m, ref_id, start_pos, end_pos, hint_bounds).map_err(
                    |err| match err {
                        crate::types::VardictError::Io(io_err) => io_err,
                        other => {
                            std::io::Error::new(std::io::ErrorKind::InvalidData, other.to_string())
                        }
                    },
                )?;
            }
        }
        let slice_records = usize::try_from(raw_record_count).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "negative CRAM slice record_count ({raw_record_count}) in '{}'",
                    path.display()
                ),
            )
        })?;
        cumulative_records = cumulative_records.saturating_add(slice_records);
        if slice_records > MAX_CRAM_CONTAINER_RECORDS
            || cumulative_records > MAX_CRAM_CONTAINER_RECORDS
        {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "CRAM slice/container record count ({cumulative_records}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_RECORDS})",
                    path.display()
                ),
            ));
        }
        if raw_record_counter < 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "negative CRAM slice record_counter ({raw_record_counter}) in '{}'",
                    path.display()
                ),
            ));
        }
        let block_count = usize::try_from(raw_block_count)
            .ok()
            .filter(|&bc| (1..=MAX_SLICE_BLOCKS).contains(&bc))
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM slice block_count ({raw_block_count}) in '{}' exceeds allowed range 1..={MAX_SLICE_BLOCKS}",
                        path.display()
                    ),
                )
            })?;
        let content_ids_len = usize::try_from(raw_content_ids_len)
            .ok()
            .filter(|&n| n <= MAX_SLICE_BLOCKS)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM slice block_content_ids length ({raw_content_ids_len}) in '{}' exceeds maximum allowed ({MAX_SLICE_BLOCKS})",
                        path.display()
                    ),
                )
            })?;
        let _ = content_ids_len;

        let mut cur_blk_offset = first_data_blk_offset;
        let mut cumulative_slice_uncomp = 0_usize;
        for blk_idx in 0..block_count {
            let expected_type = if blk_idx == 0 {
                CORE_DATA_CONTENT_TYPE
            } else {
                EXTERNAL_DATA_CONTENT_TYPE
            };
            let (next_blk_offset, blk_uncomp) = validate_cram_block_and_codec(
                path,
                file,
                cur_blk_offset,
                slice_end_offset,
                expected_type,
                MAX_CRAM_CONTAINER_BYTES,
                slice_records,
            )?;
            cumulative_slice_uncomp = cumulative_slice_uncomp.saturating_add(blk_uncomp);
            if cumulative_slice_uncomp > MAX_CRAM_CONTAINER_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM slice cumulative uncompressed block size ({cumulative_slice_uncomp}) in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_BYTES})",
                        path.display()
                    ),
                ));
            }
            cur_blk_offset = next_blk_offset;
        }
    }

    Ok(())
}

fn fetch_cram_container(
    path: &Path,
    meta: &std::sync::Arc<SharedCramMeta>,
    container_offset: u64,
    hint_contig: &[u8],
    pinned_slot: Option<&std::sync::Arc<ContigSlot>>,
) -> Result<Option<CachedCramBlock>, crate::types::VardictError> {
    fetch_cram_container_in(
        path,
        meta,
        container_offset,
        hint_contig,
        pinned_slot,
        cram_container_cache_shards(),
        MAX_CACHED_CRAM_CONTAINER_BYTES,
    )
}

fn fetch_cram_container_in(
    path: &Path,
    meta: &std::sync::Arc<SharedCramMeta>,
    container_offset: u64,
    hint_contig: &[u8],
    pinned_slot: Option<&std::sync::Arc<ContigSlot>>,
    cache: &[std::sync::Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS],
    max_cached_bytes: usize,
) -> Result<Option<CachedCramBlock>, crate::types::VardictError> {
    use std::io::{Read as _, Seek as _, SeekFrom};

    let key = (meta.id, container_offset);
    let shard_idx = cache_shard_index::<CRAM_CONTAINER_CACHE_SHARDS>(meta.id, container_offset);
    let is_global_cache = CRAM_CONTAINER_CACHE
        .get()
        .is_some_and(|global| std::ptr::eq(cache, global));
    let start_epoch =
        is_global_cache.then(|| BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire));
    if let Some(ep) = start_epoch {
        if let Some(records) = lookup_cram_container_tls(ep, key) {
            return Ok(Some(records));
        }
    }

    let inflight_guard = if let Some(shard_mutex) = cache.get(shard_idx) {
        loop {
            let waiter_opt = {
                let mut shard = shard_mutex.lock().unwrap_or_else(|e| e.into_inner());
                if let Some(records) = shard.lru.get(key) {
                    if let Some(ep) = start_epoch {
                        store_cram_container_tls(ep, key, std::sync::Arc::clone(&records));
                    }
                    return Ok(Some(records));
                }
                if let Some((_, waiter)) = shard.inflight.iter().find(|(k, _)| *k == key) {
                    Some(std::sync::Arc::clone(waiter))
                } else {
                    let inflight = std::sync::Arc::new(InflightCramDecode::new());
                    shard.inflight.push((key, std::sync::Arc::clone(&inflight)));
                    break Some(InflightCramGuard {
                        shard_mutex,
                        key,
                        entry: inflight,
                        completed: false,
                    });
                }
            };
            if let Some(waiter) = waiter_opt {
                let mut st = waiter.state.lock().unwrap_or_else(|e| e.into_inner());
                while matches!(*st, InflightCramState::Pending) {
                    st = waiter.cvar.wait(st).unwrap_or_else(|e| e.into_inner());
                }
                match *st {
                    InflightCramState::Ready(ref opt) => return Ok(opt.clone()),
                    InflightCramState::Failed | InflightCramState::Pending => continue,
                }
            }
        }
    } else {
        None
    };

    let pre_hint_bounds = pinned_slot
        .filter(|slot| slot.name.as_slice() == hint_contig)
        .and_then(|slot| match slot.loaded.get() {
            Some(Ok(loaded)) => {
                let (ref_id, _, mr) = meta.header.reference_sequences().get_full(hint_contig)?;
                let seq = loaded.seq.as_ref()?;
                Some((
                    ref_id,
                    mr.length().get(),
                    SharedFastaRepo::unpadded_slice(seq).len(),
                ))
            }
            _ => None,
        });

    let mut reader = match pop_pooled_cram_reader(path, meta) {
        Some(r) => r,
        None => noodles::cram::io::Reader::new(std::fs::File::open(path)?),
    };

    reader.seek(SeekFrom::Start(container_offset))?;
    let mut len_buf = [0_u8; 4];
    reader.get_mut().read_exact(&mut len_buf)?;
    let raw_container_len = i32::from_le_bytes(len_buf);
    if raw_container_len < 0
        || (raw_container_len as usize) > MAX_CRAM_CONTAINER_BYTES
        || meta
            .file_len
            .is_some_and(|fl| (raw_container_len as u64) > fl.saturating_sub(container_offset))
    {
        return Err(crate::types::VardictError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "CRAM container length ({raw_container_len}) at offset {container_offset} in '{}' exceeds bounds",
                path.display()
            ),
        )));
    }

    const MAX_INMEM_CRAM_CONTAINER_BYTES: usize = 4 * 1024 * 1024;
    const CRAM_HDR_SLACK_BYTES: usize = 4096;
    const MAX_RETAINED_CRAM_SCRATCH_BYTES: usize = 2 * 1024 * 1024;
    let payload_len_est = raw_container_len as usize;

    let mut container = noodles::cram::io::reader::Container::default();
    let inmem_outcome = if payload_len_est <= MAX_INMEM_CRAM_CONTAINER_BYTES {
        CRAM_IO_SCRATCH_TLS
            .try_with(|tls| -> Option<Result<Option<usize>, crate::types::VardictError>> {
                let mut scratch = tls.try_borrow_mut().ok()?;
                let want = payload_len_est.saturating_add(CRAM_HDR_SLACK_BYTES);
                scratch.clear();
                scratch.reserve(want);
                scratch.extend_from_slice(&len_buf);
                let rem_to_read = want.saturating_sub(4);
                let read_res = (&mut *reader.get_mut())
                    .take(rem_to_read as u64)
                    .read_to_end(&mut scratch);
                let Ok(_) = read_res else {
                    return None;
                };
                let mut cur = std::io::Cursor::new(scratch.as_slice());
                let mut mem_reader = noodles::cram::io::Reader::new(&mut cur);
                let bytes_read = match mem_reader.read_container(&mut container) {
                    Ok(n) => n,
                    Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                        container = noodles::cram::io::reader::Container::default();
                        return None;
                    }
                    Err(e) => return Some(Err(crate::types::VardictError::Io(e))),
                };
                if bytes_read == 0 {
                    return Some(Ok(None));
                }
                if container.header().record_count() > MAX_CRAM_CONTAINER_RECORDS {
                    return Some(Err(crate::types::VardictError::Io(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM container record_count ({}) at offset {container_offset} in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_RECORDS})",
                            container.header().record_count(),
                            path.display()
                        ),
                    ))));
                }
                let post_rel = cur.position() as usize;
                let Some(body_rel) = post_rel.checked_sub(bytes_read) else {
                    return Some(Err(crate::types::VardictError::Io(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "CRAM container stream position precedes container payload",
                    ))));
                };
                let mut val_cur =
                    std::io::Cursor::new(scratch.get(..post_rel).unwrap_or(scratch.as_slice()));
                if let Err(e) = prevalidate_cram_container_slices(
                    path,
                    &mut val_cur,
                    body_rel as u64,
                    bytes_read,
                    container.header().landmarks(),
                    Some(meta.as_ref()),
                    pre_hint_bounds,
                ) {
                    return Some(Err(crate::types::VardictError::Io(e)));
                }
                if scratch.capacity() > MAX_RETAINED_CRAM_SCRATCH_BYTES {
                    scratch.clear();
                    scratch.shrink_to(MAX_RETAINED_CRAM_SCRATCH_BYTES);
                }
                Some(Ok(Some(bytes_read)))
            })
            .ok()
            .flatten()
    } else {
        None
    };

    match inmem_outcome {
        Some(Ok(Some(_bytes_read))) => {
            return_pooled_cram_reader(path, std::sync::Arc::clone(meta), reader);
        }
        Some(Ok(None)) => {
            return_pooled_cram_reader(path, std::sync::Arc::clone(meta), reader);
            if let Some(g) = inflight_guard {
                g.complete(None, max_cached_bytes, start_epoch);
            }
            return Ok(None);
        }
        Some(Err(err)) => {
            return Err(err);
        }
        None => {
            reader.seek(SeekFrom::Start(container_offset))?;
            let bytes_read = reader.read_container(&mut container)?;
            if bytes_read == 0 {
                return_pooled_cram_reader(path, std::sync::Arc::clone(meta), reader);
                if let Some(g) = inflight_guard {
                    g.complete(None, max_cached_bytes, start_epoch);
                }
                return Ok(None);
            }
            if container.header().record_count() > MAX_CRAM_CONTAINER_RECORDS {
                return Err(crate::types::VardictError::Io(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "CRAM container record_count ({}) at offset {container_offset} in '{}' exceeds maximum allowed ({MAX_CRAM_CONTAINER_RECORDS})",
                        container.header().record_count(),
                        path.display()
                    ),
                )));
            }
            let post_container_pos = reader.get_mut().stream_position()?;
            let body_offset = post_container_pos
                .checked_sub(bytes_read as u64)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "CRAM container stream position precedes container payload",
                    )
                })?;
            {
                let mut buffered = SeekWindowReader::with_capacity(64 * 1024, reader.get_mut());
                prevalidate_cram_container_slices(
                    path,
                    &mut buffered,
                    body_offset,
                    bytes_read,
                    container.header().landmarks(),
                    Some(meta.as_ref()),
                    pre_hint_bounds,
                )?;
            }
            reader.get_mut().seek(SeekFrom::Start(post_container_pos))?;
            return_pooled_cram_reader(path, std::sync::Arc::clone(meta), reader);
        }
    }

    let (hint_repo, hint_seq) = match pinned_slot {
        Some(slot) if slot.name.as_slice() == hint_contig => {
            SharedFastaRepo::extract_from_slot(slot)
        }
        _ => meta.fasta_repo.warm_contig(hint_contig),
    };
    let hint_bounds = match (
        meta.header.reference_sequences().get_full(hint_contig),
        hint_seq.as_ref(),
    ) {
        (Some((ref_id, _, mr)), Some(seq)) => Some((
            ref_id,
            mr.length().get(),
            SharedFastaRepo::unpadded_slice(seq).len(),
        )),
        _ => None,
    };

    if let Some(spans) = meta.spans_by_offset.get(&container_offset) {
        for &(ref_id, s_idx, e_idx) in spans {
            if let Err(err) = validate_cram_span(path, meta, ref_id, s_idx, e_idx, hint_bounds) {
                meta.fasta_repo.reclaim_delegated_bytes(hint_contig);
                return Err(err);
            }
        }
    }

    crate::count!(CramContainerDecodes);
    #[cfg(test)]
    if is_global_cache {
        meta.decodes
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
    #[cfg(test)]
    if let Some(shard_mutex) = cache.get(shard_idx) {
        shard_mutex
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .decodes += 1;
    }

    // guideline deviation: M-PANIC-CONTINUATION — defensive boundary converting upstream noodles-cram parser panics on malformed untrusted CRAM inputs into io::ErrorKind::InvalidData without shared mutable state
    let decoded_res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(
        || -> Result<Vec<CachedCramRecord>, crate::types::VardictError> {
            let compression_header = container.compression_header()?;
            let record_count = container.header().record_count();
            let mut out_records = Vec::with_capacity(record_count.min(16_384));

            for slice_res in container.slices() {
                let slice = slice_res?;
                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
                let cram_records = slice.records(
                    hint_repo.clone(),
                    &meta.header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                )?;
                if out_records.len().saturating_add(cram_records.len()) > MAX_CRAM_CONTAINER_RECORDS
                {
                    return Err(crate::types::VardictError::Io(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "CRAM container decoded records in '{}' exceed maximum allowed ({MAX_CRAM_CONTAINER_RECORDS})",
                            path.display()
                        ),
                    )));
                }

                crate::count!(CramRecordsVisited, cram_records.len());
                out_records.reserve(cram_records.len());
                for record in cram_records {
                    let mut cached = cram_record_to_cached(path, &record, meta, hint_bounds)?;
                    if cached.sam.nm.is_none()
                        && cached.reference_sequence_id.is_some()
                        && cached.alignment_start.is_some()
                    {
                        match cached.reference_name.as_deref() {
                            Some(rn) if rn.as_bytes() == hint_contig => {
                                if let Some(seq) = hint_seq.as_deref() {
                                    cached.sam.nm = compute_cram_nm(
                                        &cached.sam,
                                        SharedFastaRepo::unpadded_slice(seq),
                                    );
                                }
                            }
                            Some(rn) => {
                                if let Some(ref seq) = meta.fasta_repo.warm_contig(rn.as_bytes()).1
                                {
                                    cached.sam.nm = compute_cram_nm(
                                        &cached.sam,
                                        SharedFastaRepo::unpadded_slice(seq),
                                    );
                                }
                            }
                            None => {}
                        }
                    }
                    out_records.push(cached);
                }
            }
            Ok(out_records)
        },
    ));
    meta.fasta_repo.reclaim_delegated_bytes(hint_contig);

    let mut out_records = match decoded_res {
        Ok(res) => res?,
        Err(_) => {
            return Err(crate::types::VardictError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "malformed CRAM container caused decoder panic",
            )));
        }
    };
    out_records.shrink_to_fit();

    let arc_records: CachedCramBlock =
        std::sync::Arc::new(CachedCramContainer::from_records(out_records));
    #[cfg(test)]
    if let Some(hook) = TEST_MID_CRAM_DECODE_HOOK
        .try_with(|h| h.try_borrow_mut().ok().and_then(|mut g| g.take()))
        .ok()
        .flatten()
    {
        hook();
    }
    if let Some(g) = inflight_guard {
        g.complete(
            Some(std::sync::Arc::clone(&arc_records)),
            max_cached_bytes,
            start_epoch,
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

impl std::fmt::Debug for CramFeed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramFeed")
            .field("path", &self.inner.path)
            .field("reference_sequence_id", &self.inner.reference_sequence_id)
            .field("interval", &self.inner.interval)
            .finish_non_exhaustive()
    }
}

struct CramFeedInner {
    path: PathBuf,
    meta: std::sync::Arc<SharedCramMeta>,
    reference_sequence_id: usize,
    hint_contig: Vec<u8>,
    pinned_slot: Option<std::sync::Arc<ContigSlot>>,
    interval: noodles::core::region::Interval,
    container_offsets: std::vec::IntoIter<u64>,
    current_block: Option<(CachedCramBlock, usize, Option<usize>)>,
}

impl Drop for CramFeedInner {
    fn drop(&mut self) {
        if self.pinned_slot.take().is_some() {
            self.meta.fasta_repo.evict_unpinned();
        }
    }
}

impl CramFeedInner {
    fn pump(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
        loop {
            if let Some((ref block, ref mut idx, end_bound)) = self.current_block {
                while let Some(rec) = block.get(*idx) {
                    *idx += 1;
                    if let (Some(max_end), Some(start)) = (end_bound, rec.alignment_start) {
                        if usize::from(start) > max_end {
                            break;
                        }
                    }
                    if cram_record_intersects(rec, self.reference_sequence_id, self.interval) {
                        crate::count!(CramRecordsYielded);
                        let prev_ref_name = out.reference_name.take();
                        let prev_mate_ref_name = out.mate_reference_name.take();
                        out.clone_from_record(&rec.sam);
                        out.reference_name = prev_ref_name;
                        out.mate_reference_name = prev_mate_ref_name;
                        copy_opt_arc_str(&mut out.reference_name, &rec.reference_name);
                        copy_opt_arc_str(&mut out.mate_reference_name, &rec.mate_reference_name);
                        return Ok(true);
                    }
                }
                self.current_block = None;
            }
            let Some(offset) = self.container_offsets.next() else {
                return Ok(false);
            };
            match fetch_cram_container(
                &self.path,
                &self.meta,
                offset,
                &self.hint_contig,
                self.pinned_slot.as_ref(),
            )? {
                Some(block) => {
                    let (start_idx, end_bound) = block
                        .start_index_and_upper_bound(self.reference_sequence_id, self.interval);
                    self.current_block = Some((block, start_idx, end_bound));
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
                pinned_slot: None,
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
#[derive(Debug)]
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

#[inline]
pub(crate) fn build_noodles_region(
    chr: &str,
    start_pos: u32,
    end_pos: u32,
) -> Result<noodles::core::Region, crate::types::VardictError> {
    if chr.is_empty() {
        return Err(crate::types::VardictError::Parse(format!(
            "Invalid region {chr}:{start_pos}-{end_pos}: empty chromosome name"
        )));
    }
    if start_pos > end_pos {
        return Err(crate::types::VardictError::Parse(format!(
            "Invalid region {chr}:{start_pos}-{end_pos}: start ({start_pos}) exceeds end ({end_pos})"
        )));
    }
    let s_pos = noodles::core::Position::try_from(start_pos as usize)
        .and_then(|s| noodles::core::Position::try_from(end_pos as usize).map(|e| (s, e)))
        .map_err(|e| {
            crate::types::VardictError::Parse(format!(
                "Invalid region {chr}:{start_pos}-{end_pos}: {e}"
            ))
        })?;
    Ok(noodles::core::Region::new(chr, s_pos.0..=s_pos.1))
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
    crate::count!(QueryRegionCalls);
    if is_cram_path(bam_path) {
        let meta = acquire_cram_meta(bam_path, config)?;
        let start_pos = start.to_one_based().0;
        let end_pos = end.0;
        if start_pos > end_pos {
            return Ok(RecordFeed::Cram(Box::new(CramFeed::empty(
                bam_path.to_path_buf(),
                meta,
            ))));
        }
        let region = build_noodles_region(chr, start_pos, end_pos)?;
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
        let container_offsets = meta
            .indexed_crai
            .query_offsets(reference_sequence_id, interval);
        let hint_contig = region.name().to_vec();
        let pinned_slot = if container_offsets.is_empty() {
            None
        } else {
            Some(meta.fasta_repo.warm_contig_slot(&hint_contig))
        };
        Ok(RecordFeed::Cram(Box::new(CramFeed {
            inner: CramFeedInner {
                path: bam_path.to_path_buf(),
                meta,
                reference_sequence_id,
                hint_contig,
                pinned_slot,
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
                    meta,
                    reader: Some(reader),
                    tainted: false,
                    reference_sequence_id: 0,
                    interval: noodles::core::region::Interval::from(..),
                    chunks: Vec::new().into_iter(),
                    chunk_end: None,
                    record: noodles::bam::Record::default(),
                    cigar_scratch: Vec::new(),
                    sa_buf: String::new(),
                    mc_buf: String::new(),
                    mate_ref_buf: String::new(),
                },
                sam: SamRecord::new(),
            })));
        }

        let region = match build_noodles_region(chr, start_pos, end_pos) {
            Ok(r) => r,
            Err(err) => {
                return_pooled_reader(bam_path, meta, reader);
                return Err(err);
            }
        };

        let Some(reference_sequence_id) = meta
            .header
            .reference_sequences()
            .get_index_of(region.name())
        else {
            return_pooled_reader(bam_path, meta, reader);
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
                return_pooled_reader(bam_path, meta, reader);
                return Err(crate::types::VardictError::Io(e));
            }
        };

        Ok(RecordFeed::Bam(Box::new(BamFeed {
            inner: BamFeedInner {
                meta,
                reader: Some(reader),
                tainted: false,
                reference_sequence_id,
                interval: region.interval(),
                chunks: chunks.into_iter(),
                chunk_end: None,
                record: noodles::bam::Record::default(),
                cigar_scratch: Vec::new(),
                sa_buf: String::new(),
                mc_buf: String::new(),
                mate_ref_buf: String::new(),
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

#[inline]
fn cigar_kind_char(kind: CigarKind) -> char {
    match kind {
        CigarKind::Match => 'M',
        CigarKind::Insertion => 'I',
        CigarKind::Deletion => 'D',
        CigarKind::Skip => 'N',
        CigarKind::SoftClip => 'S',
        CigarKind::HardClip => 'H',
        CigarKind::Pad => 'P',
        CigarKind::SequenceMatch => '=',
        CigarKind::SequenceMismatch => 'X',
    }
}

/// Append CIGAR to an existing buffer (scratch reuse; identical bytes to
/// [`get_cigar_string`]).
pub fn push_cigar_string(buf: &mut String, cigar: &[(u32, CigarKind)]) {
    use std::fmt::Write as _;
    for &(len, kind) in cigar {
        // write! formats in place (no temp String per op); bytes identical.
        let _ = write!(buf, "{len}");
        buf.push(cigar_kind_char(kind));
    }
}

/// Compute the aligned length on the reference using only Match and Deletion CIGAR operators.
pub fn get_aligned_length_md(cigar: &[(u32, CigarKind)]) -> i32 {
    cigar
        .iter()
        .filter(|&&(_, kind)| kind == CigarKind::Match || kind == CigarKind::Deletion)
        .fold(0i32, |acc, &(len, _)| {
            acc.saturating_add(i32::try_from(len).unwrap_or(i32::MAX))
        })
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MateRefKey {
    Unmapped,
    Same,
    Other(Box<str>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DupKey {
    MatePos(i32, MateRefKey, i32),
    UnmappedMateCigar(i32, Vec<(u32, CigarKind)>),
}

impl std::hash::Hash for DupKey {
    #[inline]
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            DupKey::MatePos(pos, mate_ref, mate_pos) => {
                0u8.hash(state);
                pos.hash(state);
                mate_ref.hash(state);
                mate_pos.hash(state);
            }
            DupKey::UnmappedMateCigar(pos, ops) => {
                1u8.hash(state);
                pos.hash(state);
                ops.len().hash(state);
                for &(len, kind) in ops {
                    len.hash(state);
                    cigar_kind_char(kind).hash(state);
                }
            }
        }
    }
}

/// Class for primary processing of SAM record and filtering reads that don't fit the requirements.
/// Ports RecordPreprocessor.java.
pub struct RecordPreprocessor {
    pub duplicates: rustc_hash::FxHashSet<DupKey>,
    pub first_matching_position: i32,
    pub total_reads: u32,
    pub duplicate_reads: u32,
}

impl RecordPreprocessor {
    pub fn new() -> Self {
        RecordPreprocessor {
            duplicates: rustc_hash::FxHashSet::default(),
            first_matching_position: -1,
            total_reads: 0,
            duplicate_reads: 0,
        }
    }

    #[inline]
    fn get_mate_ref_key(record: &SamRecord) -> MateRefKey {
        match (&record.reference_name, &record.mate_reference_name) {
            (None, _) | (_, None) => MateRefKey::Unmapped,
            (Some(_), Some(mate_ref)) if mate_ref == "*" => MateRefKey::Unmapped,
            (Some(ref_name), Some(mate_ref)) if ref_name == mate_ref || mate_ref == "=" => {
                MateRefKey::Same
            }
            (Some(_), Some(mate_ref)) => MateRefKey::Other(mate_ref.as_str().into()),
        }
    }

    /// Resets per-BAM duplicate tracking state (`RecordPreprocessor.java:80-82`).
    pub fn reset_for_next_bam(&mut self) {
        self.duplicates.clear();
        self.first_matching_position = -1;
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

        // 2. Secondary alignment filter (unless samfilter is explicitly "0", SamView.java:101)
        let is_secondary = (record.flags & 0x100) != 0;
        if is_secondary && config.skip_secondary_and_supplementary() {
            return false;
        }

        // 3. Sequence is not stored check
        if record.sequence.is_empty() || record.sequence.as_slice() == b"*" {
            return false;
        }

        // 4. SAM flag filter (default 0x504 includes 0x4 unmapped; -F 0 disables it per SamView.java:109)
        if config.flag_filter != 0 && (u32::from(record.flags) & config.flag_filter) != 0 {
            return false;
        }

        // Fused CIGAR query + reference consumption consistency check (SAMv1 §1.4.6)
        // Corrupted or hostile reads where CIGAR query consumption exceeds SEQ length
        // or quality scores length cause out-of-bounds indexing in downstream pileup/parsers.
        let mut query_consumption = 0u64;
        let mut ref_consumption = 0u64;
        for &(n, kind) in &record.cigar {
            let n64 = u64::from(n);
            match kind {
                CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                    query_consumption = query_consumption.saturating_add(n64);
                    ref_consumption = ref_consumption.saturating_add(n64);
                }
                CigarKind::Insertion | CigarKind::SoftClip => {
                    query_consumption = query_consumption.saturating_add(n64);
                }
                CigarKind::Deletion | CigarKind::Skip => {
                    ref_consumption = ref_consumption.saturating_add(n64);
                }
                CigarKind::HardClip | CigarKind::Pad => {}
            }
        }
        if query_consumption > record.sequence.len() as u64
            || (!record.quality_scores.is_empty()
                && query_consumption > record.quality_scores.len() as u64)
        {
            log::warn!(
                "Dropping corrupt BAM record {}: CIGAR query consumption ({}) exceeds SEQ length ({})",
                String::from_utf8_lossy(&record.name),
                query_consumption,
                record.sequence.len()
            );
            return false;
        }

        let start_u64 = record.alignment_start.map_or(0u64, |p| u64::from(p.0));
        let mate_start_u64 = record.mate_alignment_start.map_or(0u64, |p| u64::from(p.0));
        if start_u64 > i32::MAX as u64
            || mate_start_u64 > i32::MAX as u64
            || ref_consumption > i32::MAX as u64
            || start_u64.saturating_sub(1).saturating_add(ref_consumption) > i32::MAX as u64
        {
            log::warn!(
                "Dropping corrupt BAM record {}: reference/mate coordinate or CIGAR span (start={}, mate_start={}, span={}) exceeds i32::MAX",
                String::from_utf8_lossy(&record.name),
                start_u64,
                mate_start_u64,
                ref_consumption
            );
            return false;
        }

        self.total_reads = self.total_reads.saturating_add(1);

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
                let key = DupKey::MatePos(
                    alignment_start_val,
                    Self::get_mate_ref_key(record),
                    mate_alignment_start_val,
                );
                if self.duplicates.contains(&key) {
                    self.duplicate_reads = self.duplicate_reads.saturating_add(1);
                    return false;
                }
                self.duplicates.insert(key);
                self.first_matching_position = alignment_start_val;
            } else {
                let is_paired = (record.flags & 0x1) != 0;
                let is_mate_unmapped = (record.flags & 0x8) != 0;
                if is_paired && is_mate_unmapped {
                    // POS-CIGAR
                    let key = DupKey::UnmappedMateCigar(alignment_start_val, record.cigar.clone());
                    if self.duplicates.contains(&key) {
                        self.duplicate_reads = self.duplicate_reads.saturating_add(1);
                        return false;
                    }
                    self.duplicates.insert(key);
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
        parser.scratch_f4_qual = b"GARBAGE".to_vec();
        parser.scratch_f4_seg1 = "GARBAGE".to_string();
        parser.scratch_f4_seg2 = b"GARBAGE".to_vec();
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

    /// Task 4 pin: a reused `SamRecord` refilled via `next_into` must carry
    /// no prior aux fields — a record WITH `nM`/`NM`/`SA`/`MC`/`MQ` tags
    /// followed by a record WITHOUT them must decode the second with all
    /// aux `None` (per-tile-worker buffer reuse must never leak tags).
    #[test]
    fn bam_reuse_clears_aux() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::{
            data::field::Value as BufValue, Cigar, Data, RecordBuf,
        };
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;

        fn record(name: &str, with_aux: bool) -> RecordBuf {
            let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
            let mut builder = RecordBuf::builder()
                .set_name(name)
                .set_flags(noodles::sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(100usize).unwrap())
                .set_cigar(cigar)
                .set_sequence(vec![b'A'; 50].into())
                .set_quality_scores(vec![40u8; 50].into());
            if with_aux {
                let data: Data = [
                    (Tag::new(b'n', b'M'), BufValue::from(5i32)),
                    (Tag::EDIT_DISTANCE, BufValue::from(7i32)),
                    (Tag::new(b'S', b'A'), BufValue::from("chr1,200,+,50M,60,0;")),
                    (Tag::new(b'M', b'C'), BufValue::from("50M")),
                    (Tag::new(b'M', b'Q'), BufValue::from(60i32)),
                ]
                .into_iter()
                .collect();
                builder = builder.set_data(data);
            }
            builder.build()
        }

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("aux.bam");
        {
            let rs = HeaderMap::<RefSeqMap>::builder()
                .set_length(1000usize.try_into().unwrap())
                .build()
                .unwrap();
            let header = Header::builder().add_reference_sequence("chr1", rs).build();
            let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            let mut writer = noodles::bam::io::Writer::from(bgzf);
            writer.write_alignment_header(&header).unwrap();
            let mut indexer = noodles::csi::binning_index::Indexer::<
                noodles::csi::binning_index::index::reference_sequence::index::LinearIndex,
            >::new(14, 5)
            .expect("valid index parameters");
            for rec in [record("with_aux", true), record("plain", false)] {
                let start = writer.get_ref().virtual_position();
                writer.write_alignment_record(&header, &rec).unwrap();
                let end = writer.get_ref().virtual_position();
                let s = rec.alignment_start().expect("aligned record");
                let e = Position::try_from(100usize + 50 - 1).unwrap();
                indexer
                    .add_record(
                        Some((0, s, e, true)),
                        noodles::csi::binning_index::index::reference_sequence::bin::Chunk::new(
                            start, end,
                        ),
                    )
                    .map_err(std::io::Error::other)
                    .unwrap();
            }
            writer.try_finish().unwrap();
            noodles::bam::bai::io::Writer::new(
                std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap(),
            )
            .write_index(&indexer.build(1))
            .unwrap();
        }

        let config = Config {
            ..Default::default()
        };
        let mut feed = query_region_iter(&bam_path, "chr1", Pos0(50), Pos0(500), &config).unwrap();
        let mut reused = SamRecord::new();
        assert!(feed.next_into(&mut reused).unwrap(), "first record");
        assert_eq!(reused.name, b"with_aux");
        assert_eq!(reused.nm, Some(5), "lowercase nM wins over NM");
        assert_eq!(
            reused.sa.as_deref(),
            Some("chr1,200,+,50M,60,0;"),
            "SA survives BAM roundtrip"
        );
        assert_eq!(reused.mc.as_deref(), Some("50M"));
        assert_eq!(reused.mq, Some(60));
        assert!(feed.next_into(&mut reused).unwrap(), "second record");
        assert_eq!(reused.name, b"plain");
        assert_eq!(reused.sequence.len(), 50);
        assert_eq!(
            reused.nm, None,
            "stale NM must not leak into tagless record"
        );
        assert_eq!(
            reused.sa, None,
            "stale SA must not leak into tagless record"
        );
        assert_eq!(
            reused.mc, None,
            "stale MC must not leak into tagless record"
        );
        assert_eq!(
            reused.mq, None,
            "stale MQ must not leak into tagless record"
        );
        assert!(!feed.next_into(&mut reused).unwrap(), "feed exhausted");
    }

    fn vardict_java_test_bam(rel: &str) -> std::path::PathBuf {
        let local = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(rel);
        if local.exists() {
            return local;
        }
        if let Ok(out) = std::process::Command::new("git")
            .args(["worktree", "list", "--porcelain"])
            .output()
        {
            for line in String::from_utf8_lossy(&out.stdout).lines() {
                if let Some(wt) = line.strip_prefix("worktree ") {
                    let candidate = Path::new(wt.trim()).join(rel);
                    if candidate.exists() {
                        return candidate;
                    }
                }
            }
        }
        local
    }

    #[test]
    fn test_bam_open_and_query() {
        let bam_path_buf = vardict_java_test_bam(
            "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam",
        );
        let bam_path = bam_path_buf.as_path();
        let config = Config::default();

        let ref_seqs = read_reference_lengths(bam_path, &config).unwrap();
        assert!(
            !ref_seqs.is_empty(),
            "Reference sequences should not be empty"
        );

        let iter_records: Vec<SamRecord> =
            query_region_iter(bam_path, "chr7", Pos0(55259400), Pos0(55259600), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        let records =
            query_region(bam_path, "chr7", Pos0(55259400), Pos0(55259600), &config).unwrap();

        assert!(!records.is_empty(), "Should load records from BAM file");
        assert_eq!(records, iter_records);

        #[expect(deprecated, reason = "verify 1.0.x SemVer shim open_reader")]
        {
            let cached = open_reader(bam_path, &config).unwrap();
            assert!(!cached.header.reference_sequences().is_empty());
        }
    }

    #[test]
    fn flag_filter() {
        let bam_path_buf = vardict_java_test_bam(
            "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam",
        );
        let bam_path = bam_path_buf.as_path();

        // 1. Without flag filter
        let config_no_filter = Config {
            flag_filter: 0,
            samfilter_raw: "0".to_string(),
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

        // 3. I1: -F "0x0" has flag_filter == 0 but samfilter_raw != "0", so secondary (0x100) reads are skipped,
        //    whereas -F "0" retains secondary reads; unmapped (0x4) reads increment total_reads under -F 0.
        let mut sec_rec = SamRecord::new();
        sec_rec.flags = 0x100;
        sec_rec.cigar = vec![(10, CigarKind::Match)];
        sec_rec.sequence = vec![b'A'; 10];
        let mut prep_zero = RecordPreprocessor::new();
        assert!(prep_zero.preprocess_record(&sec_rec, &config_no_filter));
        let config_hex_zero = Config {
            flag_filter: 0,
            samfilter_raw: "0x0".to_string(),
            ..Default::default()
        };
        let mut prep_hex_zero = RecordPreprocessor::new();
        assert!(!prep_hex_zero.preprocess_record(&sec_rec, &config_hex_zero));
    }

    #[test]
    fn test_cram_twin_parity() {
        use noodles::bam::bai;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::{
            bin::Chunk, index::LinearIndex,
        };
        use noodles::csi::binning_index::Indexer;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::{
            data::field::Value as BufValue, Cigar, Data, RecordBuf,
        };
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;

        let temp_dir = tempfile::tempdir().unwrap();
        let ref_fa = temp_dir.path().join("ref.fa");
        let ref_fai = temp_dir.path().join("ref.fa.fai");
        let slice_bam = temp_dir.path().join("slice.bam");
        let cram_path = temp_dir.path().join("test.cram");

        // 1. Write a compact 1,000-bp synthetic reference FASTA and .fai index.
        let ref_len = 1000usize;
        let mut ref_seq = Vec::with_capacity(ref_len);
        for i in 0..ref_len {
            ref_seq.push(b"ACGT"[(i * 7 + 3) % 4]);
        }
        {
            let mut fa_file = File::create(&ref_fa).unwrap();
            writeln!(fa_file, ">chr7").unwrap();
            fa_file.write_all(&ref_seq).unwrap();
            fa_file.write_all(b"\n").unwrap();
        }
        {
            let mut fai_file = File::create(&ref_fai).unwrap();
            writeln!(fai_file, "chr7\t{ref_len}\t6\t{ref_len}\t{}", ref_len + 1).unwrap();
        }

        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(ref_len.try_into().unwrap())
            .build()
            .unwrap();
        let header = Header::builder().add_reference_sequence("chr7", rs).build();

        // 2. Build synthetic aligned records spanning chr7:100..=600 (matches, mismatches,
        //    indels, soft-clips, paired tags) with explicit NM tags so BAM and CRAM
        //    decode to identical SamRecords.
        let mut records = Vec::new();
        for idx in 0..24usize {
            let start_1 = 110 + idx * 15;
            let read_len = 40usize;
            let mut seq = ref_seq[start_1 - 1..start_1 - 1 + read_len].to_vec();
            let nm_val = if idx % 3 == 0 {
                // Introduce 1 mismatch at offset 10
                seq[10] = match seq[10] {
                    b'A' => b'C',
                    _ => b'A',
                };
                1i32
            } else {
                0i32
            };
            let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
            let data: Data = [
                (Tag::EDIT_DISTANCE, BufValue::from(nm_val)),
                (Tag::new(b'M', b'C'), BufValue::from("40M")),
                (Tag::new(b'M', b'Q'), BufValue::from(60i32)),
            ]
            .into_iter()
            .collect();
            let rec = RecordBuf::builder()
                .set_name(format!("read_{idx:02}"))
                .set_flags(Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT)
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start_1).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
                .set_cigar(cigar)
                .set_mate_reference_sequence_id(0)
                .set_mate_alignment_start(Position::try_from(start_1 + 100).unwrap())
                .set_template_length(140)
                .set_sequence(seq.into())
                .set_quality_scores(vec![35u8; read_len].into())
                .set_data(data)
                .build();
            records.push((start_1, read_len, rec));
        }

        // 3. Write BAM + .bam.bai in pure Rust.
        {
            let bgzf = noodles::bgzf::io::Writer::new(File::create(&slice_bam).unwrap());
            let mut bam_writer = noodles::bam::io::Writer::from(bgzf);
            bam_writer.write_alignment_header(&header).unwrap();
            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
            for &(start_1, span, ref rec) in &records {
                let s_vp = bam_writer.get_ref().virtual_position();
                bam_writer.write_alignment_record(&header, rec).unwrap();
                let e_vp = bam_writer.get_ref().virtual_position();
                let s = Position::try_from(start_1).unwrap();
                let e = Position::try_from(start_1 + span - 1).unwrap();
                indexer
                    .add_record(Some((0, s, e, true)), Chunk::new(s_vp, e_vp))
                    .unwrap();
            }
            bam_writer.try_finish().unwrap();
            bai::io::Writer::new(File::create(slice_bam.with_extension("bam.bai")).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }

        // 4. Write CRAM + .cram.crai in pure Rust.
        {
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(File::open(&ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut cram_writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            cram_writer.write_alignment_header(&header).unwrap();
            for (_, _, rec) in &records {
                cram_writer.write_alignment_record(&header, rec).unwrap();
            }
            cram_writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&cram_path).unwrap();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        // 5. Query both BAM and CRAM and verify parity.
        let config = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };

        NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, std::sync::atomic::Ordering::Relaxed);
        let bam_records = query_region(&slice_bam, "chr7", Pos0(100), Pos0(600), &config).unwrap();
        assert_eq!(bam_records.len(), records.len());
        let mut reusable_rec = SamRecord::new();

        // Retry across `BAM_CACHE_EPOCH` bumps from sibling parallel tests, explicitly
        // evicting `(meta.id, _)` at the top of each attempt so `TEST_META_ID_BASE`-pinned
        // cache entries from a partial prior attempt cannot mask `decodes_after_1 > decodes_before`.
        loop {
            let ep0 = BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire);
            let meta = acquire_cram_meta(&cram_path, &config).unwrap();
            if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) != ep0 {
                continue;
            }
            for shard_mutex in cram_container_cache_shards() {
                shard_mutex
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .lru
                    .entries
                    .retain(|((mid, _), _)| *mid != meta.id);
            }
            let decodes_before = meta.decodes.load(std::sync::atomic::Ordering::Relaxed);
            let cram_records =
                query_region(&cram_path, "chr7", Pos0(100), Pos0(600), &config).unwrap();
            if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) != ep0 {
                continue;
            }

            assert_eq!(
                bam_records.len(),
                cram_records.len(),
                "BAM and CRAM must have the exact same number of reads"
            );
            for (i, (bam_rec, cram_rec)) in bam_records.iter().zip(cram_records.iter()).enumerate()
            {
                assert_eq!(bam_rec, cram_rec, "Parity mismatch at record {}", i);
            }

            // 6. Verify streaming CramFeed (`query_region_iter`) matches `query_region`
            //    and serves repeated sub-window queries from `CRAM_CONTAINER_CACHE` with 0 re-decodes.
            let decodes_after_1 = meta.decodes.load(std::sync::atomic::Ordering::Relaxed);
            assert!(
                decodes_after_1 > decodes_before,
                "first CRAM query must decode container"
            );
            let mut feed1 =
                query_region_iter(&cram_path, "chr7", Pos0(100), Pos0(600), &config).unwrap();
            let mut iter_records_1 = Vec::new();
            while feed1.next_into(&mut reusable_rec).unwrap() {
                iter_records_1.push(reusable_rec.clone());
            }
            if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) != ep0 {
                continue;
            }
            assert_eq!(iter_records_1, cram_records);

            let iter_records_2: Vec<SamRecord> =
                query_region_iter(&cram_path, "chr7", Pos0(180), Pos0(420), &config)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            let expected_sub: Vec<SamRecord> =
                query_region(&cram_path, "chr7", Pos0(180), Pos0(420), &config).unwrap();
            if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) != ep0 {
                continue;
            }
            assert_eq!(iter_records_2, expected_sub);
            let decodes_after_2 = meta.decodes.load(std::sync::atomic::Ordering::Relaxed);
            assert_eq!(
                decodes_after_2, decodes_after_1,
                "second CRAM sub-window query must hit CRAM_CONTAINER_CACHE without re-decoding"
            );
            break;
        }

        // 7. Verify empty-window short-circuit (`start > end`) and missing-contig error on CRAM.
        let mut empty_cram =
            query_region_iter(&cram_path, "chr7", Pos0(600), Pos0(100), &config).unwrap();
        assert!(!empty_cram.next_into(&mut reusable_rec).unwrap());
        assert!(
            query_region_iter(&cram_path, "chr_missing", Pos0(10), Pos0(100), &config).is_err()
        );

        // 8. Multi-slice CRAM container regression (#46): a `.crai` index with multiple
        //    slice entries pointing to the same container byte offset must deduplicate
        //    the container offset in `query_region_iter` and never duplicate reads.
        let multi_cram_path = temp_dir.path().join("multi_slice.cram");
        std::fs::copy(&cram_path, &multi_cram_path).unwrap();
        let base_crai = noodles::cram::fs::read_associated_index(&cram_path).unwrap();
        assert!(!base_crai.is_empty());
        let first_entry = &base_crai[0];
        let dup_crai = vec![
            first_entry.clone(),
            noodles::cram::crai::Record::new(
                first_entry.reference_sequence_id(),
                Position::try_from(250usize).ok(),
                first_entry.alignment_span(),
                first_entry.offset(),
                first_entry.landmark(),
                first_entry.slice_length(),
            ),
        ];
        noodles::cram::crai::fs::write(multi_cram_path.with_extension("cram.crai"), &dup_crai)
            .unwrap();

        let multi_crai = noodles::cram::fs::read_associated_index(&multi_cram_path).unwrap();
        assert!(
            multi_crai
                .windows(2)
                .any(|w| w[0].offset() == w[1].offset()),
            "multi_slice.cram.crai must contain repeated container offsets across slices"
        );
        let multi_records: Vec<SamRecord> =
            query_region_iter(&multi_cram_path, "chr7", Pos0(100), Pos0(600), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert_eq!(
            multi_records.len(),
            bam_records.len(),
            "multi-slice CRAM container must not duplicate reads"
        );
        assert_eq!(multi_records, bam_records);
    }

    #[test]
    fn flag_filter_preserves_high_u32_bits_and_saturates_qual() {
        // High u32 bit 0x10000 (65536) in `config.flag_filter` must not truncate to u16 (0):
        // `config.flag_filter != 0` stays true so secondary alignments (0x100) are filtered.
        let mut preprocessor = RecordPreprocessor::new();
        let config = Config {
            flag_filter: 0x10000,
            ..Config::default()
        };

        let secondary_rec = SamRecord {
            name: b"sec".to_vec(),
            flags: 0x100,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(10)),
            alignment_end: Some(Pos1(13)),
            mapping_quality: Some(60),
            cigar: vec![(4, CigarKind::Match)],
            sequence: b"ACGT".to_vec(),
            quality_scores: vec![40; 4],
            ..SamRecord::default()
        };
        assert!(
            !preprocessor.preprocess_record(&secondary_rec, &config),
            "flag_filter=0x10000 must remain non-zero (not truncated to u16) and reject secondary reads"
        );

        let primary_rec = SamRecord {
            flags: 0,
            ..secondary_rec
        };
        assert!(preprocessor.preprocess_record(&primary_rec, &config));

        // Verify quality scores saturate at u8::MAX (255) on score.saturating_add(33)
        // and saturate at 0 on byte.saturating_sub(33).
        use noodles::sam::alignment::record_buf::RecordBuf;
        use noodles::sam::Header;
        let header = Header::default();
        let rec_buf = RecordBuf::builder()
            .set_sequence(b"AC".to_vec().into())
            .set_quality_scores(vec![0u8, 250u8].into())
            .build();
        let sam_rec = SamRecord::from_alignment_record(&rec_buf, &header).unwrap();
        assert_eq!(sam_rec.quality_scores, vec![33u8, 255u8]);
        for raw_byte in [0u8, 1, 32, 33] {
            assert_eq!(raw_byte.saturating_sub(33), 0);
        }
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
            for res in record.cigar().iter() {
                let Ok(op) = res else {
                    break;
                };
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

    impl BamFeedLegacy {
        fn next_into(&mut self, out: &mut SamRecord) -> Result<bool, crate::types::VardictError> {
            self.inner.pump(out)
        }
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

    fn query_region_iter_legacy(
        bam_path: &Path,
        chr: &str,
        start: Pos0,
        end: Pos0,
        _config: &crate::config::Config,
    ) -> Result<Box<BamFeedLegacy>, crate::types::VardictError> {
        let index =
            read_bam_index_with_fallbacks(bam_path).map_err(crate::types::VardictError::Io)?;
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_path(bam_path)
            .map_err(crate::types::VardictError::Io)?;
        let header = reader
            .read_header()
            .map_err(crate::types::VardictError::Io)?;

        let start_pos = start.to_one_based().0;
        let end_pos = end.0;
        if start_pos > end_pos {
            return Ok(Box::new(BamFeedLegacy {
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
            }));
        }

        let region = build_noodles_region(chr, start_pos, end_pos)?;

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
        Ok(Box::new(BamFeedLegacy {
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
        }))
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
        .expect_err("missing contig must fail");
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
            .expect_err("missing BAM with empty window must error");
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
            .expect_err("missing .bai must error");
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
            .expect_err("truncated BAM header must error");
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

        const TEST_META_ID: u64 = u64::MAX - 42;

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

    #[test]
    fn record_preprocessor_reset_for_next_bam_clears_duplicates_while_preserving_totals() {
        // Issue #60 (#3): CigarParser.java:108-120 calls `preprocessor.reset()` at the start
        // of each BAM in a multi-BAM colon-separated run (`combineAnalysis` Tumor:Normal),
        // clearing `duplicates` and resetting `firstMatchingPosition = -1` while preserving
        // `totalReads` and `duplicateReads` counters (`RecordPreprocessor.java:80-82`).
        let mut prep = RecordPreprocessor::new();
        let config = Config {
            dedup: true,
            ..Default::default()
        };
        let rec = SamRecord {
            flags: 0x9,
            alignment_start: Some(Pos1(12_214_074)),
            mate_alignment_start: Some(Pos1(12_214_215)),
            template_length: 142,
            cigar: vec![(151, CigarKind::Match)],
            sequence: vec![b'A'; 151],
            ..SamRecord::default()
        };
        assert!(prep.preprocess_record(&rec, &config));
        assert!(!prep.preprocess_record(&rec, &config));
        assert_eq!(prep.total_reads, 2);
        assert_eq!(prep.duplicate_reads, 1);

        prep.reset_for_next_bam();
        assert!(
            prep.preprocess_record(&rec, &config),
            "identical coordinates in second BAM must not be filtered as duplicate of first BAM"
        );
        assert_eq!(prep.total_reads, 3);
        assert_eq!(prep.duplicate_reads, 1);
    }

    #[test]
    fn preprocess_record_eqx_uses_query_consumed_len() {
        let mut prep = RecordPreprocessor::new();
        let config = Config::default();
        let valid_eqx = SamRecord {
            name: b"valid_eqx".to_vec(),
            flags: 0,
            alignment_start: Some(Pos1(10)),
            cigar: vec![
                (10, CigarKind::SequenceMatch),
                (1, CigarKind::SequenceMismatch),
                (10, CigarKind::SequenceMatch),
            ],
            sequence: vec![b'A'; 21],
            quality_scores: vec![b'I'; 21],
            ..SamRecord::default()
        };
        assert!(prep.preprocess_record(&valid_eqx, &config));

        let corrupt_eqx = SamRecord {
            name: b"corrupt_eqx".to_vec(),
            flags: 0,
            alignment_start: Some(Pos1(10)),
            cigar: vec![
                (10, CigarKind::SequenceMatch),
                (1, CigarKind::SequenceMismatch),
                (10, CigarKind::SequenceMatch),
            ],
            sequence: vec![b'A'; 15],
            quality_scores: vec![b'I'; 15],
            ..SamRecord::default()
        };
        assert!(
            !prep.preprocess_record(&corrupt_eqx, &config),
            "corrupt =/X read with query_consumed_len (21) > sequence.len() (15) must be dropped"
        );
    }

    fn write_test_fasta_with_fai(dir: &Path, name: &str, contigs: &[(&str, &[u8])]) -> PathBuf {
        use std::io::Write as _;
        let fa_path = dir.join(name);
        let fai_path = fasta_fai_path(&fa_path);
        let mut fa = std::fs::File::create(&fa_path).unwrap();
        let mut fai = std::fs::File::create(&fai_path).unwrap();
        let mut offset = 0u64;
        for &(chr, seq) in contigs {
            let hdr = format!(">{chr}\n");
            fa.write_all(hdr.as_bytes()).unwrap();
            offset += hdr.len() as u64;
            fa.write_all(seq).unwrap();
            fa.write_all(b"\n").unwrap();
            let len = seq.len() as u64;
            writeln!(fai, "{chr}\t{len}\t{offset}\t{len}\t{}", len + 1).unwrap();
            offset += len + 1;
        }
        fa_path
    }

    #[test]
    fn cram_short_reference_returns_err() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let temp_dir = tempfile::tempdir().unwrap();
        let full_seq = b"ACGT".repeat(250); // 1000 bp
        let short_seq = &full_seq[..500];
        let full_fa = write_test_fasta_with_fai(temp_dir.path(), "full.fa", &[("chr1", &full_seq)]);
        let short_fa =
            write_test_fasta_with_fai(temp_dir.path(), "short.fa", &[("chr1", short_seq)]);

        let sq_map = Map::<ReferenceSequence>::builder()
            .set_length(NonZeroUsize::try_from(1000).unwrap())
            .build()
            .unwrap();
        let header = noodles::sam::Header::builder()
            .add_reference_sequence("chr1", sq_map)
            .build();

        let cram_path = temp_dir.path().join("short_ref_test.cram");
        {
            let fasta_repo = acquire_fasta_repo(&full_fa).unwrap();
            let repo = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                primary_name: Vec::new(),
                primary_record: None,
                shared: std::sync::Arc::downgrade(&fasta_repo),
            });
            let mut writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            writer.write_alignment_header(&header).unwrap();
            let rec = RecordBuf::builder()
                .set_name("r_past")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(600usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60).unwrap())
                .set_cigar(
                    vec![Op::new(Kind::Match, 50)]
                        .into_iter()
                        .collect::<Cigar>(),
                )
                .set_sequence(Sequence::from(full_seq[599..649].to_vec()))
                .set_quality_scores(QualityScores::from(vec![30u8; 50]))
                .build();
            writer.write_alignment_record(&header, &rec).unwrap();
            writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&cram_path).unwrap();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        let short_config = Config {
            reference: short_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        let iter_res: Result<Vec<_>, _> =
            query_region_iter(&cram_path, "chr1", Pos0(550), Pos0(700), &short_config)
                .and_then(|feed| feed.collect());
        let err = iter_res.expect_err("short FASTA must return Err, not panic");
        assert!(
            err.to_string().contains("exceeding FASTA"),
            "unexpected error: {err}"
        );

        let vec_err = query_region(&cram_path, "chr1", Pos0(550), Pos0(700), &short_config)
            .expect_err("query_region with short FASTA must return Err");
        assert!(
            vec_err.to_string().contains("exceeding FASTA"),
            "unexpected error: {vec_err}"
        );
    }

    #[test]
    fn cram_missing_m5_short_ref_still_errors() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let temp_dir = tempfile::tempdir().unwrap();
        let full_seq = b"ACGT".repeat(250); // 1000 bp
        let short_seq = &full_seq[..500];
        let full_fa = write_test_fasta_with_fai(
            temp_dir.path(),
            "full_multi.fa",
            &[("chr1", &full_seq), ("chr2", &full_seq)],
        );
        let short_fa = write_test_fasta_with_fai(
            temp_dir.path(),
            "short_multi.fa",
            &[("chr1", short_seq), ("chr2", &full_seq)],
        );

        let sq1 = Map::<ReferenceSequence>::builder()
            .set_length(NonZeroUsize::try_from(1000).unwrap())
            .build()
            .unwrap();
        let sq2 = Map::<ReferenceSequence>::builder()
            .set_length(NonZeroUsize::try_from(1000).unwrap())
            .build()
            .unwrap();
        let header = noodles::sam::Header::builder()
            .add_reference_sequence("chr1", sq1)
            .add_reference_sequence("chr2", sq2)
            .build();

        // Multi-ref container: noodles sets `ReferenceSequenceContext::Many` and omits slice MD5.
        let cram_path = temp_dir.path().join("multi_ref_no_slice_md5.cram");
        {
            let fasta_repo = acquire_fasta_repo(&full_fa).unwrap();
            let repo = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                primary_name: Vec::new(),
                primary_record: None,
                shared: std::sync::Arc::downgrade(&fasta_repo),
            });
            let mut writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            writer.write_alignment_header(&header).unwrap();
            let rec1 = RecordBuf::builder()
                .set_name("r_chr1_past")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(600usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60).unwrap())
                .set_cigar(
                    vec![Op::new(Kind::Match, 50)]
                        .into_iter()
                        .collect::<Cigar>(),
                )
                .set_sequence(Sequence::from(full_seq[599..649].to_vec()))
                .set_quality_scores(QualityScores::from(vec![30u8; 50]))
                .build();
            let rec2 = RecordBuf::builder()
                .set_name("r_chr2_ok")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(1)
                .set_alignment_start(Position::try_from(100usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60).unwrap())
                .set_cigar(
                    vec![Op::new(Kind::Match, 50)]
                        .into_iter()
                        .collect::<Cigar>(),
                )
                .set_sequence(Sequence::from(full_seq[99..149].to_vec()))
                .set_quality_scores(QualityScores::from(vec![30u8; 50]))
                .build();
            writer.write_alignment_record(&header, &rec1).unwrap();
            writer.write_alignment_record(&header, &rec2).unwrap();
            writer.try_finish(&header).unwrap();

            // Write a .crai index entry for the multi-ref container (noodles::cram::fs::index
            // uses an empty Repository::default() on multi-ref slices).
            let mut rdr = noodles::cram::io::Reader::new(std::fs::File::open(&cram_path).unwrap());
            rdr.read_header().unwrap();
            let offset = rdr.position().unwrap();
            let mut container = noodles::cram::io::reader::Container::default();
            let container_len = rdr.read_container(&mut container).unwrap();
            let landmark = container.header().landmarks()[0] as u64;
            let slice_len = (container_len as u64).saturating_sub(landmark);
            let crai = vec![
                noodles::cram::crai::Record::new(
                    Some(0),
                    Position::try_from(1usize).ok(),
                    1,
                    offset,
                    landmark,
                    slice_len,
                ),
                noodles::cram::crai::Record::new(
                    Some(1),
                    Position::try_from(100usize).ok(),
                    50,
                    offset,
                    landmark,
                    slice_len,
                ),
            ];
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        let short_config = Config {
            reference: short_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        let iter_res: Result<Vec<_>, _> =
            query_region_iter(&cram_path, "chr1", Pos0(0), Pos0(700), &short_config)
                .and_then(|feed| feed.collect());
        let err = iter_res.expect_err("multi-ref slice with short FASTA must return Err");
        assert!(
            err.to_string().contains("exceeding FASTA"),
            "unexpected error: {err}"
        );
        let vec_err = query_region(&cram_path, "chr1", Pos0(0), Pos0(700), &short_config)
            .expect_err("query_region on multi-ref slice with short FASTA must return Err");
        assert!(
            vec_err.to_string().contains("exceeding FASTA"),
            "unexpected error: {vec_err}"
        );
    }

    #[test]
    fn cram_short_reference_in_range_query_ok() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let temp_dir = tempfile::tempdir().unwrap();
        let full_seq = b"ACGT".repeat(250); // 1000 bp
        let short_seq = &full_seq[..500];
        let full_fa =
            write_test_fasta_with_fai(temp_dir.path(), "full_inrange.fa", &[("chr1", &full_seq)]);
        let short_fa =
            write_test_fasta_with_fai(temp_dir.path(), "short_inrange.fa", &[("chr1", short_seq)]);

        let sq_map = Map::<ReferenceSequence>::builder()
            .set_length(NonZeroUsize::try_from(1000).unwrap())
            .build()
            .unwrap();
        let header = noodles::sam::Header::builder()
            .add_reference_sequence("chr1", sq_map)
            .build();

        let cram_path = temp_dir.path().join("in_range.cram");
        {
            let fasta_repo = acquire_fasta_repo(&full_fa).unwrap();
            let repo = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                primary_name: Vec::new(),
                primary_record: None,
                shared: std::sync::Arc::downgrade(&fasta_repo),
            });
            let mut writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            writer.write_alignment_header(&header).unwrap();
            let rec = RecordBuf::builder()
                .set_name("r_in_range")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(100usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60).unwrap())
                .set_cigar(
                    vec![Op::new(Kind::Match, 50)]
                        .into_iter()
                        .collect::<Cigar>(),
                )
                .set_sequence(Sequence::from(full_seq[99..149].to_vec()))
                .set_quality_scores(QualityScores::from(vec![30u8; 50]))
                .build();
            writer.write_alignment_record(&header, &rec).unwrap();
            writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&cram_path).unwrap();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        let short_config = Config {
            reference: short_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        let recs: Vec<SamRecord> =
            query_region_iter(&cram_path, "chr1", Pos0(80), Pos0(200), &short_config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name, b"r_in_range");
    }

    #[test]
    fn cram_overhang_read_matches_bam_twin() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let temp_dir = tempfile::tempdir().unwrap();
        let ref_seq = b"ACGT".repeat(250); // 1000 bp
        let ref_fa =
            write_test_fasta_with_fai(temp_dir.path(), "overhang_ref.fa", &[("chr1", &ref_seq)]);

        let sq_map = Map::<ReferenceSequence>::builder()
            .set_length(NonZeroUsize::try_from(1000).unwrap())
            .build()
            .unwrap();
        let header = noodles::sam::Header::builder()
            .add_reference_sequence("chr1", sq_map)
            .build();

        // Build overhanging reads [981..=1030] (50M, overhanging LN=1000 by 30 bp) with a SNV at pos 990 (idx 9).
        let mut overhang_seq = ref_seq[980..1000].to_vec();
        overhang_seq[9] = if overhang_seq[9] == b'A' { b'T' } else { b'A' };
        overhang_seq.extend_from_slice(&[b'A'; 30]);

        let mut records = Vec::new();
        for i in 0..4 {
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::{data::field::Value as BufValue, Data};
            let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(31i32))]
                .into_iter()
                .collect();
            let rec = RecordBuf::builder()
                .set_name(format!("overhang_{i}"))
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(981usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60).unwrap())
                .set_cigar(
                    vec![Op::new(Kind::Match, 50)]
                        .into_iter()
                        .collect::<Cigar>(),
                )
                .set_sequence(Sequence::from(overhang_seq.clone()))
                .set_quality_scores(QualityScores::from(vec![35u8; 50]))
                .set_data(data)
                .build();
            records.push(rec);
        }

        let bam_path = temp_dir.path().join("overhang.bam");
        {
            let mut bam_writer =
                noodles::bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            bam_writer.write_alignment_header(&header).unwrap();
            for rec in &records {
                bam_writer.write_alignment_record(&header, rec).unwrap();
            }
            bam_writer.try_finish().unwrap();
            let mut reader = noodles::bam::io::Reader::new(std::fs::File::open(&bam_path).unwrap());
            reader.read_header().unwrap();
            let mut record = noodles::bam::Record::default();
            let mut indexer = noodles::csi::binning_index::Indexer::<
                noodles::csi::binning_index::index::reference_sequence::index::LinearIndex,
            >::default();
            let mut start_vp = reader.get_ref().virtual_position();
            while reader.read_record(&mut record).unwrap() != 0 {
                let end_vp = reader.get_ref().virtual_position();
                let chunk = noodles::csi::binning_index::index::reference_sequence::bin::Chunk::new(
                    start_vp, end_vp,
                );
                let align_ctx = (
                    record.reference_sequence_id().transpose().unwrap().unwrap(),
                    record.alignment_start().transpose().unwrap().unwrap(),
                    record.alignment_end().transpose().unwrap().unwrap(),
                    !record.flags().is_unmapped(),
                );
                indexer.add_record(Some(align_ctx), chunk).unwrap();
                start_vp = end_vp;
            }
            noodles::bam::bai::io::Writer::new(
                std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap(),
            )
            .write_index(&indexer.build(1))
            .unwrap();
        }

        let cram_path = temp_dir.path().join("overhang.cram");
        {
            let fasta_repo = acquire_fasta_repo(&ref_fa).unwrap();
            let repo = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                primary_name: Vec::new(),
                primary_record: None,
                shared: std::sync::Arc::downgrade(&fasta_repo),
            });
            let mut cram_writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            cram_writer.write_alignment_header(&header).unwrap();
            for rec in &records {
                cram_writer.write_alignment_record(&header, rec).unwrap();
            }
            cram_writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&cram_path).unwrap();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        let config = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        let bam_recs: Vec<SamRecord> =
            query_region_iter(&bam_path, "chr1", Pos0(950), Pos0(1000), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        let cram_recs: Vec<SamRecord> =
            query_region_iter(&cram_path, "chr1", Pos0(950), Pos0(1000), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert_eq!(bam_recs.len(), 4);
        assert_eq!(bam_recs, cram_recs);

        let mut bam_out = Vec::new();
        let bam_run_cfg = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            bam: bam_path.to_str().unwrap().to_string(),
            region: Some("chr1:970-1000".to_string()),
            ..Default::default()
        };
        crate::run_with_writer(&bam_run_cfg, &mut bam_out).unwrap();

        let mut cram_out = Vec::new();
        let cram_run_cfg = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            bam: cram_path.to_str().unwrap().to_string(),
            region: Some("chr1:970-1000".to_string()),
            ..Default::default()
        };
        crate::run_with_writer(&cram_run_cfg, &mut cram_out).unwrap();

        assert!(!bam_out.is_empty(), "SNV at pos 990 must be called");
        assert_eq!(
            String::from_utf8(bam_out).unwrap(),
            String::from_utf8(cram_out).unwrap()
        );
    }

    #[test]
    fn cram_base_lookup_matches_htsjdk() {
        let mut expected = [b'N'; 127];
        for &b in b"=ABCDGHKMNRSTVWY" {
            expected[b as usize] = b;
            expected[(b + 32) as usize] = b;
        }
        for byte in 0u8..=255u8 {
            if byte < 127 {
                assert_eq!(
                    normalize_cram_base(byte),
                    Some(expected[byte as usize]),
                    "mismatch at byte 0x{byte:02x}"
                );
            } else {
                assert_eq!(
                    normalize_cram_base(byte),
                    None,
                    "byte 0x{byte:02x} >= 0x7f must return None"
                );
            }
        }
    }

    #[test]
    fn cram_non_ascii_base_and_non_iupac_base_normalization() {
        use noodles::core::Position;
        use noodles::cram;
        use noodles::fasta;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{Cigar as CigarBuf, QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use noodles::sam::header::record::value::Map;
        use std::num::NonZeroUsize;

        clear_bam_caches();
        let tmp = tempfile::tempdir().unwrap();
        let ref_seq: Vec<u8> = b"ACGTACGTAC".iter().copied().cycle().take(200).collect();
        let ref_fa = tmp.path().join("ref.fa");
        std::fs::write(
            &ref_fa,
            format!(">chr1\n{}\n", std::str::from_utf8(&ref_seq).unwrap()),
        )
        .unwrap();
        std::fs::write(tmp.path().join("ref.fa.fai"), "chr1\t200\t6\t200\t201\n").unwrap();

        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(200).unwrap()),
            )
            .build();
        let repo = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("chr1", None),
            fasta::record::Sequence::from(ref_seq.clone()),
        )]);

        // 1. Non-ASCII base 0xc9 in CRAM returns Err("non-ASCII base")
        let bad_cram = tmp.path().join("bad_c9.cram");
        {
            let mut seq = ref_seq[19..49].to_vec();
            seq[10] = 0xc9;
            let rec = RecordBuf::builder()
                .set_name("r_c9")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(20).unwrap())
                .set_mapping_quality(sam::alignment::record::MappingQuality::try_from(60).unwrap())
                .set_cigar(
                    [
                        Op::new(Kind::Match, 10),
                        Op::new(Kind::Insertion, 1),
                        Op::new(Kind::Match, 19),
                    ]
                    .into_iter()
                    .collect::<CigarBuf>(),
                )
                .set_sequence(Sequence::from(seq))
                .set_quality_scores(QualityScores::from(vec![35u8; 30]))
                .build();
            let file = std::fs::File::create(&bad_cram).unwrap();
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo.clone())
                .build_from_writer(file);
            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();
            writer.write_alignment_record(&header, &rec).unwrap();
            writer.try_finish(&header).unwrap();
        }
        let bad_index = cram::fs::index(&bad_cram).unwrap();
        cram::crai::fs::write(bad_cram.with_extension("cram.crai"), &bad_index).unwrap();

        let config = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        let err = query_region_iter(&bad_cram, "chr1", Pos0(1), Pos0(200), &config)
            .and_then(|it| it.collect::<Result<Vec<_>, _>>())
            .expect_err("non-ASCII CRAM base 0xc9 must fail with Parse error");
        assert!(
            err.to_string().contains("non-ASCII base 0xc9"),
            "unexpected error: {err}"
        );

        // 2. Non-IUPAC ASCII base 'X' in CRAM normalizes to 'N'
        let x_cram = tmp.path().join("base_x.cram");
        {
            let mut seq = ref_seq[19..49].to_vec();
            seq[10] = b'X';
            let rec = RecordBuf::builder()
                .set_name("r_x")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(20).unwrap())
                .set_mapping_quality(sam::alignment::record::MappingQuality::try_from(60).unwrap())
                .set_cigar(
                    [
                        Op::new(Kind::Match, 10),
                        Op::new(Kind::Insertion, 1),
                        Op::new(Kind::Match, 19),
                    ]
                    .into_iter()
                    .collect::<CigarBuf>(),
                )
                .set_sequence(Sequence::from(seq))
                .set_quality_scores(QualityScores::from(vec![35u8; 30]))
                .build();
            let file = std::fs::File::create(&x_cram).unwrap();
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_writer(file);
            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();
            writer.write_alignment_record(&header, &rec).unwrap();
            writer.try_finish(&header).unwrap();
        }
        let x_index = cram::fs::index(&x_cram).unwrap();
        cram::crai::fs::write(x_cram.with_extension("cram.crai"), &x_index).unwrap();

        let recs: Vec<SamRecord> = query_region_iter(&x_cram, "chr1", Pos0(1), Pos0(200), &config)
            .unwrap()
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].sequence[10], b'N');
    }

    #[test]
    fn crai_offsets_match_linear_oracle() {
        use noodles::core::Position;
        use noodles::cram::crai;

        let mut records: Vec<crai::Record> = Vec::with_capacity(12_000);
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next_u64 = || -> u64 {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state.wrapping_mul(0x2545_F491_4F6C_DD1D)
        };

        // 0..4: sorted reference groups with nested spans, zero-span rows, None-start rows,
        // and multi-slice containers sharing the same offset.
        for ref_id in 0..4usize {
            let mut cur_start = 10usize;
            for i in 0..2500usize {
                if i % 250 == 0 {
                    // Unplaced / zero-span / None-start edge rows
                    records.push(crai::Record::new(
                        None,
                        Position::new(cur_start),
                        100,
                        (ref_id as u64) * 100_000 + (i as u64),
                        0,
                        100,
                    ));
                    records.push(crai::Record::new(
                        Some(ref_id),
                        None,
                        100,
                        (ref_id as u64) * 100_000 + (i as u64) + 1,
                        0,
                        100,
                    ));
                    records.push(crai::Record::new(
                        Some(ref_id),
                        Position::new(cur_start),
                        0,
                        (ref_id as u64) * 100_000 + (i as u64) + 2,
                        0,
                        100,
                    ));
                }
                let step = (next_u64() % 40) as usize;
                cur_start += step;
                // Occasionally emit a very wide span so subsequent narrower records nest inside it
                let span = if i % 37 == 0 {
                    5_000usize
                } else {
                    50 + ((next_u64() % 300) as usize)
                };
                // Share container offsets across adjacent slices every 5 entries
                let offset = ((ref_id as u64) << 32) | ((i as u64 / 5) * 4096);
                records.push(crai::Record::new(
                    Some(ref_id),
                    Position::new(cur_start),
                    span,
                    offset,
                    (i % 5) as u64 * 512,
                    512,
                ));
            }
        }

        // ref_id = 4: intentionally unsorted alignment_start order
        for i in 0..1000usize {
            let start = 1 + ((next_u64() % 50_000) as usize);
            let span = 10 + ((next_u64() % 500) as usize);
            let offset = (4u64 << 32) | ((i as u64 / 3) * 2048);
            records.push(crai::Record::new(
                Some(4),
                Position::new(start),
                span,
                offset,
                0,
                256,
            ));
        }

        assert!(records.len() >= 10_000);
        let indexed = IndexedCrai::from_index(&records);

        let linear_oracle =
            |ref_id: usize, interval: noodles::core::region::Interval| -> Vec<u64> {
                let mut out = Vec::new();
                for rec in records
                    .iter()
                    .filter(|rec| crai_record_intersects(rec, ref_id, interval))
                {
                    let off = rec.offset();
                    if !out.contains(&off) {
                        out.push(off);
                    }
                }
                out
            };

        // Check 100+ random, boundary, point, and unbounded intervals across all ref_ids
        for q in 0..110usize {
            let ref_id = q % 6; // 5 is an empty/missing ref_id
            let interval: noodles::core::region::Interval = if q % 25 == 0 {
                (..).into()
            } else if q % 25 == 1 {
                let p = Position::new(1 + ((next_u64() % 60_000) as usize)).unwrap();
                (p..=p).into()
            } else {
                let a = 1 + ((next_u64() % 60_000) as usize);
                let b = a + ((next_u64() % 2_000) as usize);
                let s = Position::new(a).unwrap();
                let e = Position::new(b).unwrap();
                (s..=e).into()
            };
            let got = indexed.query_offsets(ref_id, interval);
            let expected = linear_oracle(ref_id, interval);
            assert_eq!(
                got, expected,
                "IndexedCrai mismatch at query {q} on ref_id={ref_id}, interval={interval:?}"
            );
        }
    }

    #[test]
    fn warm_contig_loads_once_under_concurrency() {
        use std::io::Write as _;

        let temp_dir = tempfile::tempdir().unwrap();
        let ref_fa = temp_dir.path().join("conc_ref.fa");
        let ref_fai = temp_dir.path().join("conc_ref.fa.fai");
        let seq1 = b"ACGT".repeat(250); // 1000 bp
        let seq2 = b"TGCA".repeat(250); // 1000 bp
        {
            let mut f = std::fs::File::create(&ref_fa).unwrap();
            writeln!(f, ">chr1").unwrap();
            f.write_all(&seq1).unwrap();
            writeln!(f, "\n>chr2").unwrap();
            f.write_all(&seq2).unwrap();
            writeln!(f).unwrap();
        }
        {
            let mut f = std::fs::File::create(&ref_fai).unwrap();
            writeln!(f, "chr1\t1000\t6\t1000\t1001").unwrap();
            writeln!(f, "chr2\t1000\t1013\t1000\t1001").unwrap();
        }

        let adapter = PreallocatedFastaAdapter::open(&ref_fa).unwrap();
        let repo = SharedFastaRepo::new(ref_fa.clone(), adapter, DEFAULT_FASTA_CACHE_BUDGET_BYTES);

        let loads_before = repo.load_count.load(std::sync::atomic::Ordering::SeqCst);
        let barrier = std::sync::Arc::new(std::sync::Barrier::new(8));
        let mut handles = Vec::new();
        for _ in 0..8 {
            let repo_c = std::sync::Arc::clone(&repo);
            let bar_c = std::sync::Arc::clone(&barrier);
            handles.push(std::thread::spawn(move || {
                bar_c.wait();
                let (_r, seq) = repo_c.warm_contig(b"chr1");
                seq.expect("chr1 sequence must load")
            }));
        }
        let mut results = Vec::new();
        for h in handles {
            results.push(h.join().unwrap());
        }
        let loads_after = repo.load_count.load(std::sync::atomic::Ordering::SeqCst);
        assert_eq!(
            loads_after - loads_before,
            1,
            "8 concurrent threads warming chr1 must invoke the underlying loader exactly once"
        );
        for s in &results {
            assert_eq!(SharedFastaRepo::unpadded_slice(s), seq1.as_slice());
            assert!(std::sync::Arc::ptr_eq(s, &results[0]));
        }

        // Warming chr2 must load chr2 once while keeping chr1 resident
        let (_, s2) = repo.warm_contig(b"chr2");
        assert_eq!(
            SharedFastaRepo::unpadded_slice(&s2.unwrap()),
            seq2.as_slice()
        );
        assert_eq!(
            repo.load_count.load(std::sync::atomic::Ordering::SeqCst) - loads_after,
            1
        );
    }

    #[test]
    fn multi_ref_slice_contig_loaded_once() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use noodles::sam::header::record::value::Map;
        use std::io::Write as _;
        use std::num::NonZeroUsize;

        let temp_dir = tempfile::tempdir().unwrap();
        let ref_fa = temp_dir.path().join("multiref.fa");
        let ref_fai = temp_dir.path().join("multiref.fa.fai");
        let cram_path = temp_dir.path().join("multiref.cram");

        let contigs = [("chr1", b'A'), ("chr2", b'C'), ("chr3", b'G')];
        let contig_len = 600usize;
        {
            let mut fa = std::fs::File::create(&ref_fa).unwrap();
            let mut fai = std::fs::File::create(&ref_fai).unwrap();
            let mut offset = 0usize;
            for (name, base) in contigs {
                let header_line = format!(">{name}\n");
                fa.write_all(header_line.as_bytes()).unwrap();
                offset += header_line.len();
                let seq: Vec<u8> = (0..contig_len)
                    .map(|i| if i % 7 == 0 { b'T' } else { base })
                    .collect();
                fa.write_all(&seq).unwrap();
                fa.write_all(b"\n").unwrap();
                writeln!(
                    fai,
                    "{name}\t{contig_len}\t{offset}\t{contig_len}\t{}",
                    contig_len + 1
                )
                .unwrap();
                offset += contig_len + 1;
            }
        }

        let mut header_builder = noodles::sam::Header::builder();
        for (name, _) in contigs {
            header_builder = header_builder.add_reference_sequence(
                name,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(contig_len).unwrap()),
            );
        }
        let header = header_builder.build();

        // Build records across all 3 contigs WITHOUT NM tags so both slice.records
        // and compute_cram_nm require the reference sequence.
        let mut records = Vec::new();
        for (ref_id, (_, base)) in contigs.iter().enumerate() {
            for idx in 0..4usize {
                let start_1 = 50 + idx * 40;
                let read_len = 40usize;
                let mut seq: Vec<u8> = (start_1 - 1..start_1 - 1 + read_len)
                    .map(|i| if i % 7 == 0 { b'T' } else { *base })
                    .collect();
                if idx == 1 {
                    seq[5] = b'A' ^ b'C' ^ seq[5];
                }
                let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
                let rec = RecordBuf::builder()
                    .set_name(format!("r_{ref_id}_{idx}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(ref_id)
                    .set_alignment_start(Position::try_from(start_1).unwrap())
                    .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
                    .set_cigar(cigar)
                    .set_sequence(seq.into())
                    .set_quality_scores(vec![35u8; read_len].into())
                    .build();
                records.push(rec);
            }
        }

        {
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(std::fs::File::open(&ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut cram_writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            cram_writer.write_alignment_header(&header).unwrap();
            for rec in &records {
                cram_writer.write_alignment_record(&header, rec).unwrap();
            }
            cram_writer.try_finish(&header).unwrap();

            // `noodles::cram::fs::index` passes `fasta::Repository::default()` (`// TODO`)
            // on multi-reference slices; read the container offset + landmark directly.
            let mut reader =
                noodles::cram::io::Reader::new(std::fs::File::open(&cram_path).unwrap());
            let _ = reader.read_header().unwrap();
            let container_pos = reader.position().unwrap();
            let mut container = noodles::cram::io::reader::Container::default();
            let container_len = reader.read_container(&mut container).unwrap();
            assert!(container_len > 0);
            let landmark = container.header().landmarks()[0] as u64;
            let slice_len = (container_len - container.header().landmarks()[0]) as u64;
            let crai: Vec<noodles::cram::crai::Record> = (0..contigs.len())
                .map(|ref_id| {
                    noodles::cram::crai::Record::new(
                        Some(ref_id),
                        Position::try_from(50usize).ok(),
                        200,
                        container_pos,
                        landmark,
                        slice_len,
                    )
                })
                .collect();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        let config = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, std::sync::atomic::Ordering::Relaxed);
        let meta = acquire_cram_meta(&cram_path, &config).unwrap();
        let loads_before = meta
            .fasta_repo
            .load_count
            .load(std::sync::atomic::Ordering::SeqCst);

        for (name, _) in contigs {
            let recs: Vec<SamRecord> =
                query_region_iter(&cram_path, name, Pos0(1), Pos0(500), &config)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            assert_eq!(recs.len(), 4);
            assert_eq!(recs[0].nm, Some(0));
            assert_eq!(recs[1].nm, Some(1));
            assert_eq!(recs[0].reference_name.as_deref(), Some(name));
        }

        let loads_after = meta
            .fasta_repo
            .load_count
            .load(std::sync::atomic::Ordering::SeqCst);
        assert_eq!(
            loads_after - loads_before,
            3,
            "each of the 3 contigs must be loaded from disk only once across multi-ref slice decode + NM calculation"
        );
    }

    #[test]
    fn cram_feed_pins_contig_under_512mb_budget() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use noodles::sam::header::record::value::Map;
        use std::io::Write as _;
        use std::num::NonZeroUsize;

        let temp_dir = tempfile::tempdir().unwrap();
        let ref_fa = temp_dir.path().join("pin_budget.fa");
        let ref_fai = temp_dir.path().join("pin_budget.fa.fai");
        let cram_path = temp_dir.path().join("pin_budget.cram");

        let contigs = [("chr1", b'A'), ("chr2", b'C'), ("chr3", b'G')];
        let contig_len = 1000usize;
        {
            let mut fa = std::fs::File::create(&ref_fa).unwrap();
            let mut fai = std::fs::File::create(&ref_fai).unwrap();
            let mut offset = 0usize;
            for (name, base) in contigs {
                let header_line = format!(">{name}\n");
                fa.write_all(header_line.as_bytes()).unwrap();
                offset += header_line.len();
                let seq = vec![base; contig_len];
                fa.write_all(&seq).unwrap();
                fa.write_all(b"\n").unwrap();
                writeln!(
                    fai,
                    "{name}\t{contig_len}\t{offset}\t{contig_len}\t{}",
                    contig_len + 1
                )
                .unwrap();
                offset += contig_len + 1;
            }
        }

        let mut header_builder = noodles::sam::Header::builder();
        for (name, _) in contigs {
            header_builder = header_builder.add_reference_sequence(
                name,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(contig_len).unwrap()),
            );
        }
        let header = header_builder.build();

        let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
        let rec = RecordBuf::builder()
            .set_name("pinned_read_1")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(100usize).unwrap())
            .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
            .set_cigar(cigar)
            .set_sequence(vec![b'A'; 50].into())
            .set_quality_scores(vec![35u8; 50].into())
            .build();

        {
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(std::fs::File::open(&ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut cram_writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            cram_writer.write_alignment_header(&header).unwrap();
            cram_writer.write_alignment_record(&header, &rec).unwrap();
            cram_writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&cram_path).unwrap();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        let config = Config {
            reference: ref_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, std::sync::atomic::Ordering::Relaxed);
        let meta = acquire_cram_meta(&cram_path, &config).unwrap();
        // Set budget to 2 * (contig_len + CRAM_REF_OVERHANG_PAD) so at most 2 padded contigs fit
        let per_contig_bytes = contig_len + CRAM_REF_OVERHANG_PAD;
        meta.fasta_repo.set_budget_bytes(2 * per_contig_bytes);

        // Hold an active CramFeed on chr1 -> pins chr1's ContigSlot
        let mut feed = query_region_iter(&cram_path, "chr1", Pos0(50), Pos0(200), &config).unwrap();
        let slot_chr1 = meta
            .fasta_repo
            .resident_slot(b"chr1")
            .expect("chr1 must be resident");
        // Drop our local clone after checking strong_count: registry (1) + feed.pinned_slot (1) + local (1) = 3
        assert!(std::sync::Arc::strong_count(&slot_chr1) >= 3);
        drop(slot_chr1);

        // Warm chr2 (fits in 2 * contig_len budget alongside chr1)
        let _ = meta.fasta_repo.warm_contig(b"chr2");
        assert!(meta.fasta_repo.resident_slot(b"chr1").is_some());
        assert!(meta.fasta_repo.resident_slot(b"chr2").is_some());

        // Warm chr3 -> exceeds 2 * contig_len budget; chr1 is pinned (strong_count == 2),
        // so unpinned chr2 (strong_count == 1) must be evicted while chr1 stays resident!
        let _ = meta.fasta_repo.warm_contig(b"chr3");
        let pinned_chr1 = meta
            .fasta_repo
            .resident_slot(b"chr1")
            .expect("pinned chr1 must NOT be evicted under budget pressure");
        assert!(std::sync::Arc::strong_count(&pinned_chr1) >= 3);
        drop(pinned_chr1);
        assert!(
            meta.fasta_repo.resident_slot(b"chr2").is_none(),
            "unpinned chr2 must be evicted to stay within budget"
        );
        assert!(
            meta.fasta_repo.resident_slot(b"chr3").is_some(),
            "newly warmed chr3 must remain resident"
        );

        let mut out = SamRecord::new();
        assert!(feed.next_into(&mut out).unwrap());
        assert_eq!(out.reference_name.as_deref(), Some("chr1"));
        assert!(!feed.next_into(&mut out).unwrap());
    }

    #[test]
    fn build_noodles_region_and_typed_dedup_keys() {
        // 1. Validate `build_noodles_region` on valid, special-character, and invalid inputs.
        let reg = build_noodles_region("HLA-A*01:01:01:01", 10, 250).unwrap();
        assert_eq!(reg.name(), "HLA-A*01:01:01:01");
        assert_eq!(
            reg.interval().start().map(usize::from),
            Some(10),
            "1-based start preserved"
        );
        assert_eq!(
            reg.interval().end().map(usize::from),
            Some(250),
            "1-based end preserved"
        );
        assert!(
            build_noodles_region("", 10, 250).is_err(),
            "empty chromosome name must error"
        );
        assert!(
            build_noodles_region("chr1", 0, 250).is_err(),
            "0 start_pos must fail Position::try_from"
        );
        assert!(
            build_noodles_region("chr1", 10, 0).is_err(),
            "0 end_pos must fail Position::try_from"
        );

        // 2. Validate `RecordPreprocessor` typed `DupKey` duplicate detection across all 4 branches.
        let config = Config {
            dedup: true,
            ..Default::default()
        };
        let mut prep = RecordPreprocessor::new();

        let base_rec = SamRecord {
            flags: 0x1,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(100)),
            alignment_end: Some(Pos1(149)),
            mapping_quality: Some(60),
            cigar: vec![(50, CigarKind::Match)],
            sequence: vec![b'A'; 50],
            quality_scores: vec![40; 50],
            ..SamRecord::default()
        };

        // Branch (a): Concordant paired reads (`MateRefKey::Same`) with PNEXT < 10
        let concordant = SamRecord {
            mate_reference_name: Some("chr1".to_string()),
            mate_alignment_start: Some(Pos1(5)),
            ..base_rec.clone()
        };
        assert!(prep.preprocess_record(&concordant, &config));
        assert!(
            !prep.preprocess_record(&concordant, &config),
            "second identical concordant read must be marked duplicate"
        );
        assert!(prep
            .duplicates
            .contains(&DupKey::MatePos(100, MateRefKey::Same, 5)));

        // Branch (b): Unmapped mate with PNEXT < 10 (`MateRefKey::Unmapped`)
        let unmapped_pnext_low = SamRecord {
            flags: 0x1 | 0x8,
            mate_reference_name: None,
            mate_alignment_start: Some(Pos1(5)),
            ..base_rec.clone()
        };
        assert!(
            prep.preprocess_record(&unmapped_pnext_low, &config),
            "unmapped mate at same POS/PNEXT must not collide with MateRefKey::Same"
        );
        assert!(!prep.preprocess_record(&unmapped_pnext_low, &config));
        assert!(prep
            .duplicates
            .contains(&DupKey::MatePos(100, MateRefKey::Unmapped, 5)));

        // Branch (c): Discordant mate contigs (`chr2` vs `chr3`) at the same POS=100 and PNEXT=5
        // (`MateRefKey::Other`), confirming `chr2` and `chr3` do NOT collide!
        let discordant_chr2 = SamRecord {
            mate_reference_name: Some("chr2".to_string()),
            mate_alignment_start: Some(Pos1(5)),
            ..base_rec.clone()
        };
        let discordant_chr3 = SamRecord {
            mate_reference_name: Some("chr3".to_string()),
            mate_alignment_start: Some(Pos1(5)),
            ..base_rec.clone()
        };
        assert!(prep.preprocess_record(&discordant_chr2, &config));
        assert!(
            prep.preprocess_record(&discordant_chr3, &config),
            "discordant mate on chr3 must NOT collide with chr2 at same POS and PNEXT < 10"
        );
        assert!(!prep.preprocess_record(&discordant_chr2, &config));
        assert!(!prep.preprocess_record(&discordant_chr3, &config));
        assert!(prep.duplicates.contains(&DupKey::MatePos(
            100,
            MateRefKey::Other("chr2".into()),
            5
        )));
        assert!(prep.duplicates.contains(&DupKey::MatePos(
            100,
            MateRefKey::Other("chr3".into()),
            5
        )));

        // Branch (d): Paired unmapped mate with PNEXT >= 10 (`DupKey::UnmappedMateCigar`)
        let unmapped_pnext_high_50m = SamRecord {
            flags: 0x1 | 0x8,
            alignment_start: Some(Pos1(200)),
            alignment_end: Some(Pos1(249)),
            mate_reference_name: Some("chr1".to_string()),
            mate_alignment_start: Some(Pos1(200)),
            cigar: vec![(50, CigarKind::Match)],
            ..base_rec.clone()
        };
        let unmapped_pnext_high_clip = SamRecord {
            cigar: vec![(10, CigarKind::SoftClip), (40, CigarKind::Match)],
            ..unmapped_pnext_high_50m.clone()
        };
        assert!(prep.preprocess_record(&unmapped_pnext_high_50m, &config));
        assert!(
            !prep.preprocess_record(&unmapped_pnext_high_50m, &config),
            "duplicate unmapped-mate CIGAR at PNEXT >= 10 must be filtered"
        );
        assert!(
            prep.preprocess_record(&unmapped_pnext_high_clip, &config),
            "distinct CIGAR at same POS with PNEXT >= 10 must not collide"
        );
        assert!(prep.duplicates.contains(&DupKey::UnmappedMateCigar(
            200,
            vec![(50, CigarKind::Match)]
        )));
        assert_eq!(prep.total_reads, 11);
        assert_eq!(prep.duplicate_reads, 5);
    }

    #[test]
    fn validate_cram_span_hint_bounds_matches_uncached_without_locking_fasta_repo() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;
        use std::sync::atomic::Ordering;

        let temp_dir = tempfile::tempdir().unwrap();
        let full_1000 = b"ACGT".repeat(250);
        let full_600 = b"TGCA".repeat(150);
        let full_fa = write_test_fasta_with_fai(
            temp_dir.path(),
            "hint_full.fa",
            &[("chr1", &full_1000), ("chr2", &full_600)],
        );
        // Short chr1 (500 bp vs header LN=1000), full chr2 (600 bp == header LN=600)
        let test_fa = write_test_fasta_with_fai(
            temp_dir.path(),
            "hint_test.fa",
            &[("chr1", &full_1000[..500]), ("chr2", &full_600)],
        );

        let sq1 = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1000).unwrap());
        let sq2 = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(600).unwrap());
        let header = noodles::sam::Header::builder()
            .add_reference_sequence("chr1", sq1)
            .add_reference_sequence("chr2", sq2)
            .build();

        let cram_path = temp_dir.path().join("hint_bounds.cram");
        {
            let fasta_repo = acquire_fasta_repo(&full_fa).unwrap();
            let repo = noodles::fasta::Repository::new(RegistryDelegatingAdapter {
                primary_name: Vec::new(),
                primary_record: None,
                shared: std::sync::Arc::downgrade(&fasta_repo),
            });
            let mut writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            writer.write_alignment_header(&header).unwrap();
            let rec = RecordBuf::builder()
                .set_name("r1")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(10usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60).unwrap())
                .set_cigar([Op::new(Kind::Match, 20)].into_iter().collect::<Cigar>())
                .set_sequence(Sequence::from(full_1000[9..29].to_vec()))
                .set_quality_scores(QualityScores::from(vec![35u8; 20]))
                .build();
            writer.write_alignment_record(&header, &rec).unwrap();
            writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&cram_path).unwrap();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
        }

        NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, Ordering::Relaxed);
        let config = Config {
            reference: test_fa.to_str().unwrap().to_string(),
            ..Default::default()
        };
        let meta = acquire_cram_meta(&cram_path, &config).unwrap();
        assert_eq!(meta.fasta_repo.load_count.load(Ordering::SeqCst), 0);
        assert!(meta.fasta_repo.resident_slot(b"chr1").is_none());
        assert!(meta.fasta_repo.resident_slot(b"chr2").is_none());

        let chr1_cases = [
            (1usize, 500usize, true),
            (500, 500, true),
            (490, 501, false), // fai_len (500) < sq_ln (1000): any span_end > fai_len rejected
            (501, 510, false), // span_start > fai_len rejected
        ];
        let hint_chr1 = Some((0usize, 1000usize, 500usize));
        let mut hinted_chr1_results = Vec::new();
        for &(s, e, expect_ok) in &chr1_cases {
            let res = validate_cram_span(&cram_path, &meta, 0, s, e, hint_chr1)
                .map_err(|err| err.to_string());
            assert_eq!(res.is_ok(), expect_ok, "chr1 span {s}..={e}: {res:?}");
            hinted_chr1_results.push(res);
        }

        let chr2_cases = [
            (1usize, 600usize, true),
            (590, 600 + CRAM_REF_OVERHANG_PAD, true), // fai_len == sq_ln: allowed up to pad
            (590, 600 + CRAM_REF_OVERHANG_PAD + 1, false), // past pad rejected
            (601, 610, false),                        // span_start > fai_len rejected
        ];
        let hint_chr2 = Some((1usize, 600usize, 600usize));
        let mut hinted_chr2_results = Vec::new();
        for &(s, e, expect_ok) in &chr2_cases {
            let res = validate_cram_span(&cram_path, &meta, 1, s, e, hint_chr2)
                .map_err(|err| err.to_string());
            assert_eq!(res.is_ok(), expect_ok, "chr2 span {s}..={e}: {res:?}");
            hinted_chr2_results.push(res);
        }

        // Neither in-bounds Ok(()) nor overhang Err(...) may lock or warm fasta_repo
        assert_eq!(
            meta.fasta_repo.load_count.load(Ordering::SeqCst),
            0,
            "matching hint_bounds must validate without locking or warming fasta_repo"
        );
        assert!(meta.fasta_repo.resident_slot(b"chr1").is_none());
        assert!(meta.fasta_repo.resident_slot(b"chr2").is_none());

        // Compare every outcome and exact error string against uncached fallback (hint_bounds = None)
        for (idx, &(s, e, _)) in chr1_cases.iter().enumerate() {
            let uncached =
                validate_cram_span(&cram_path, &meta, 0, s, e, None).map_err(|err| err.to_string());
            assert_eq!(hinted_chr1_results[idx], uncached);
        }
        for (idx, &(s, e, _)) in chr2_cases.iter().enumerate() {
            let uncached =
                validate_cram_span(&cram_path, &meta, 1, s, e, None).map_err(|err| err.to_string());
            assert_eq!(hinted_chr2_results[idx], uncached);
        }

        // Even after contigs are warm in registry, matching hint_bounds must not lock registry
        // (verified by LRU order [chr1, chr2] remaining untouched when validating chr1).
        assert_eq!(
            meta.fasta_repo.registry.lock().unwrap().order,
            vec![b"chr1".to_vec(), b"chr2".to_vec()]
        );
        validate_cram_span(&cram_path, &meta, 0, 10, 50, hint_chr1).unwrap();
        assert!(validate_cram_span(&cram_path, &meta, 0, 490, 510, hint_chr1).is_err());
        assert_eq!(
            meta.fasta_repo.registry.lock().unwrap().order,
            vec![b"chr1".to_vec(), b"chr2".to_vec()],
            "matching hint_bounds must not touch fasta_repo.registry LRU order"
        );

        // Non-matching reference_sequence_id falls back cleanly, including invalid ref_id
        let bad_id_hinted =
            validate_cram_span(&cram_path, &meta, 99, 1, 10, hint_chr1).unwrap_err();
        let bad_id_uncached = validate_cram_span(&cram_path, &meta, 99, 1, 10, None).unwrap_err();
        assert_eq!(bad_id_hinted.to_string(), bad_id_uncached.to_string());
    }

    #[test]
    fn bam_feed_next_into_preserves_contig_names_and_untouched_eof_cigar() {
        use noodles::bam::bai;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::{
            bin::Chunk, index::LinearIndex,
        };
        use noodles::csi::binning_index::Indexer;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;

        // 1. SamRecord::clear must not wipe reference_name / mate_reference_name to Some("")
        {
            let mut rec = SamRecord {
                reference_name: Some("chr1".to_string()),
                mate_reference_name: Some("chr2".to_string()),
                cigar: vec![(10, CigarKind::Match)],
                ..SamRecord::default()
            };
            rec.clear();
            assert_eq!(
                rec.reference_name.as_deref(),
                Some("chr1"),
                "SamRecord::clear must retain reference_name so same-contig refills elide clear+push_str"
            );
            assert_eq!(
                rec.mate_reference_name.as_deref(),
                Some("chr2"),
                "SamRecord::clear must retain mate_reference_name so same-contig refills elide clear+push_str"
            );
            assert!(rec.cigar.is_empty());
        }

        // 2. Build a multi-contig indexed BAM with:
        //    - chr1 reads covering concordant mate, cross-contig mate, unmapped mate,
        //      zero-ref-span CIGAR (`10I`, span == 0 => end == start), AND a trailing
        //      non-intersecting read in the same BAI leaf bin (`pos=400`, `25M`) to
        //      prove `out.cigar` is untouched when `next_into` returns `Ok(false)` at EOF.
        //    - chr2 reads for alternating multi-contig buffer reuse.
        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("multicontig_feed.bam");
        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(1000usize.try_into().unwrap())
            .build()
            .unwrap();
        let header = Header::builder()
            .add_reference_sequence("chr1", rs.clone())
            .add_reference_sequence("chr2", rs.clone())
            .add_reference_sequence("HLA-B*07:02", rs)
            .build();

        let mk_rec = |name: &str,
                      ref_id: usize,
                      start: usize,
                      ops: &[(usize, Kind)],
                      seq_len: usize,
                      mate: Option<(usize, usize)>|
         -> RecordBuf {
            let cigar: Cigar = ops.iter().map(|&(n, k)| Op::new(k, n)).collect();
            let mut b = RecordBuf::builder()
                .set_name(name)
                .set_flags(if mate.is_some() {
                    Flags::SEGMENTED | Flags::FIRST_SEGMENT
                } else {
                    Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_UNMAPPED
                })
                .set_reference_sequence_id(ref_id)
                .set_alignment_start(Position::try_from(start).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
                .set_cigar(cigar)
                .set_sequence(vec![b'A'; seq_len].into())
                .set_quality_scores(vec![35u8; seq_len].into());
            if let Some((m_ref, m_pos)) = mate {
                b = b
                    .set_mate_reference_sequence_id(m_ref)
                    .set_mate_alignment_start(Position::try_from(m_pos).unwrap());
            }
            b.build()
        };

        let all_recs = vec![
            mk_rec("c1_same", 0, 100, &[(10, Kind::Match)], 10, Some((0, 150))),
            mk_rec(
                "c1_cross",
                0,
                110,
                &[(5, Kind::Match), (2, Kind::Insertion), (3, Kind::Match)],
                10,
                Some((1, 200)),
            ),
            mk_rec(
                "c1_unmapped_mate",
                0,
                120,
                &[(8, Kind::Match), (2, Kind::Deletion), (2, Kind::Match)],
                10,
                None,
            ),
            mk_rec(
                "c1_zero_span",
                0,
                130,
                &[(10, Kind::Insertion)],
                10,
                Some((2, 50)),
            ),
            // Trailing non-intersecting record inside the same BAI bin [0, 16384):
            mk_rec(
                "c1_trailing_outside",
                0,
                400,
                &[(25, Kind::Match)],
                25,
                Some((0, 450)),
            ),
            mk_rec("c2_first", 1, 100, &[(12, Kind::Match)], 12, Some((0, 100))),
            mk_rec(
                "c2_second",
                1,
                115,
                &[(4, Kind::SoftClip), (8, Kind::Match)],
                12,
                None,
            ),
        ];

        {
            let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            let mut w = noodles::bam::io::Writer::from(bgzf);
            w.write_alignment_header(&header).unwrap();
            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
            for rec in &all_recs {
                let s_vp = w.get_ref().virtual_position();
                w.write_alignment_record(&header, rec).unwrap();
                let e_vp = w.get_ref().virtual_position();
                let ref_id = rec.reference_sequence_id().unwrap();
                let s = rec.alignment_start().unwrap();
                let span = rec.cigar().alignment_span();
                let e = if span == 0 {
                    s
                } else {
                    Position::new(s.get() + span - 1).unwrap()
                };
                indexer
                    .add_record(Some((ref_id, s, e, true)), Chunk::new(s_vp, e_vp))
                    .unwrap();
            }
            w.try_finish().unwrap();
            bai::io::Writer::new(
                std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap(),
            )
            .write_index(&indexer.build(3))
            .unwrap();
        }

        let config = Config::default();
        let expected_chr1 = query_region(&bam_path, "chr1", Pos0(99), Pos0(200), &config).unwrap();
        let expected_chr2 = query_region(&bam_path, "chr2", Pos0(99), Pos0(200), &config).unwrap();
        assert_eq!(expected_chr1.len(), 4);
        assert_eq!(expected_chr2.len(), 2);
        // Verify span == 0 => end == start on c1_zero_span
        assert_eq!(expected_chr1[3].alignment_start, Some(Pos1(130)));
        assert_eq!(expected_chr1[3].alignment_end, Some(Pos1(130)));

        // (a) Fresh SamRecord::new() on every call + sentinel cigar preserved on EOF
        {
            let mut feed =
                query_region_iter(&bam_path, "chr1", Pos0(99), Pos0(200), &config).unwrap();
            for exp in &expected_chr1 {
                let mut fresh = SamRecord::new();
                assert!(feed.next_into(&mut fresh).unwrap());
                assert_eq!(&fresh, exp);
            }
            let mut eof_sentinel = SamRecord::new();
            eof_sentinel.cigar = vec![(99, CigarKind::SoftClip)];
            assert!(!feed.next_into(&mut eof_sentinel).unwrap());
            assert_eq!(
                eof_sentinel.cigar,
                vec![(99, CigarKind::SoftClip)],
                "next_into returning Ok(false) after scanning trailing non-intersecting chunk reads must leave out.cigar untouched"
            );
        }

        // (b) Reused SamRecord across chr1: exact names + same-contig String buffer reuse + EOF cigar untouched
        {
            let mut feed =
                query_region_iter(&bam_path, "chr1", Pos0(99), Pos0(200), &config).unwrap();
            let mut reused = SamRecord::new();
            let mut prev_ref_ptr: Option<*const u8> = None;
            for (idx, exp) in expected_chr1.iter().enumerate() {
                assert!(feed.next_into(&mut reused).unwrap());
                assert_eq!(&reused, exp, "reused mismatch at #{idx}");
                let ptr = reused.reference_name.as_deref().unwrap().as_ptr();
                if let Some(prev) = prev_ref_ptr {
                    assert_eq!(
                        ptr, prev,
                        "same-contig refill must retain reference_name allocation in place"
                    );
                }
                prev_ref_ptr = Some(ptr);
            }
            assert_eq!(reused.cigar, vec![(10, CigarKind::Insertion)]);
            assert!(!feed.next_into(&mut reused).unwrap());
            assert_eq!(
                reused.cigar,
                vec![(10, CigarKind::Insertion)],
                "reused out.cigar must remain untouched when next_into returns Ok(false) at EOF"
            );
        }

        // (c) Alternating multi-contig feeds sharing a single SamRecord buffer
        {
            let mut feed1 =
                query_region_iter(&bam_path, "chr1", Pos0(99), Pos0(200), &config).unwrap();
            let mut feed2 =
                query_region_iter(&bam_path, "chr2", Pos0(99), Pos0(200), &config).unwrap();
            let mut shared = SamRecord::new();

            assert!(feed1.next_into(&mut shared).unwrap());
            assert_eq!(&shared, &expected_chr1[0]);
            assert!(feed2.next_into(&mut shared).unwrap());
            assert_eq!(&shared, &expected_chr2[0]);
            assert!(feed1.next_into(&mut shared).unwrap());
            assert_eq!(&shared, &expected_chr1[1]);
            assert!(feed2.next_into(&mut shared).unwrap());
            assert_eq!(&shared, &expected_chr2[1]);
            assert!(!feed2.next_into(&mut shared).unwrap());
            // shared.cigar must still hold expected_chr2[1].cigar after feed2 hits EOF
            assert_eq!(shared.cigar, expected_chr2[1].cigar);
            assert!(feed1.next_into(&mut shared).unwrap());
            assert_eq!(&shared, &expected_chr1[2]);
            assert!(feed1.next_into(&mut shared).unwrap());
            assert_eq!(&shared, &expected_chr1[3]);
            assert!(!feed1.next_into(&mut shared).unwrap());
            assert_eq!(shared.cigar, expected_chr1[3].cigar);
        }

        // (d) BAM_SEQ_PAIR_LUT exhaustive 256-entry verification, even/odd IUPAC nibble decode,
        //     empty QUAL = *, saturating Phred qualities, and truncated raw sequence safety.
        {
            const EXPECTED_NIBBLES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
            for byte in 0u16..=255 {
                let b = byte as u8;
                let expected_pair = [
                    EXPECTED_NIBBLES[(b >> 4) as usize],
                    EXPECTED_NIBBLES[(b & 0x0f) as usize],
                ];
                assert_eq!(BAM_SEQ_PAIR_LUT[usize::from(b)], expected_pair);
            }

            // Truncated raw sequence slices never panic or over-reserve:
            let mut trunc_buf = Vec::new();
            decode_bam_raw_seq_into(&[0x12], 4, &mut trunc_buf);
            assert_eq!(trunc_buf, b"AC");
            decode_bam_raw_seq_into(&[0x12], 3, &mut trunc_buf);
            assert_eq!(trunc_buf, b"AC");
            decode_bam_raw_seq_into(&[], 5, &mut trunc_buf);
            assert!(trunc_buf.is_empty());
            let mut huge_trunc = Vec::new();
            decode_bam_raw_seq_into(&[0x12], usize::MAX, &mut huge_trunc);
            assert_eq!(huge_trunc, b"AC");
            assert!(huge_trunc.capacity() <= 64);

            let iupac_bam_path = dir.path().join("iupac_qual_feed.bam");
            let mk_custom_rec = |name: &str,
                                 start: usize,
                                 seq_bytes: &[u8],
                                 quals_opt: Option<&[u8]>|
             -> RecordBuf {
                let cigar: Cigar = if seq_bytes.is_empty() {
                    Cigar::default()
                } else {
                    [Op::new(Kind::Match, seq_bytes.len())]
                        .into_iter()
                        .collect()
                };
                let mut b = RecordBuf::builder()
                    .set_name(name)
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(start).unwrap())
                    .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
                    .set_cigar(cigar)
                    .set_sequence(seq_bytes.to_vec().into());
                if let Some(q) = quals_opt {
                    b = b.set_quality_scores(q.to_vec().into());
                }
                b.build()
            };

            let mut custom_recs = vec![
                // Even 16-base sequence covering all 16 BAM IUPAC nibbles + valid [0..=93] qualities
                mk_custom_rec(
                    "even_all16",
                    100,
                    b"=ACMGRSVTWYHKDBN",
                    Some(&[0, 1, 10, 20, 30, 33, 40, 50, 60, 70, 80, 85, 90, 91, 92, 93]),
                ),
                // Odd 17-base sequence covering all 16 BAM IUPAC nibbles + trailing '=' AND empty QUAL = *
                mk_custom_rec("odd17_qual_star", 120, b"=ACMGRSVTWYHKDBN=", None),
                // Odd 15-base sequence + Phred qualities
                mk_custom_rec(
                    "odd15_qual",
                    140,
                    b"=ACMGRSVTWYHKDB",
                    Some(&[93, 92, 91, 90, 60, 45, 35, 0, 1, 32, 33, 40, 50, 75, 93]),
                ),
            ];
            // Also test 1-base odd singleton reads for each of the 16 IUPAC nibbles (trailing high nibble path)
            for (idx, &nib_ch) in EXPECTED_NIBBLES.iter().enumerate() {
                custom_recs.push(mk_custom_rec(
                    &format!("odd1_nib_{idx}"),
                    200 + idx * 2,
                    &[nib_ch],
                    if idx % 2 == 0 { None } else { Some(&[93u8]) },
                ));
            }

            {
                let bgzf =
                    noodles::bgzf::io::Writer::new(std::fs::File::create(&iupac_bam_path).unwrap());
                let mut w = noodles::bam::io::Writer::from(bgzf);
                w.write_alignment_header(&header).unwrap();
                let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
                for rec in &custom_recs {
                    let s_vp = w.get_ref().virtual_position();
                    w.write_alignment_record(&header, rec).unwrap();
                    let e_vp = w.get_ref().virtual_position();
                    let ref_id = rec.reference_sequence_id().unwrap();
                    let s = rec.alignment_start().unwrap();
                    let span = rec.cigar().alignment_span().max(1);
                    let e = Position::new(s.get() + span - 1).unwrap();
                    indexer
                        .add_record(Some((ref_id, s, e, true)), Chunk::new(s_vp, e_vp))
                        .unwrap();
                }
                w.try_finish().unwrap();
                bai::io::Writer::new(
                    std::fs::File::create(iupac_bam_path.with_extension("bam.bai")).unwrap(),
                )
                .write_index(&indexer.build(3))
                .unwrap();
            }

            let expected_custom =
                query_region(&iupac_bam_path, "chr1", Pos0(90), Pos0(300), &config).unwrap();
            assert_eq!(expected_custom.len(), custom_recs.len());

            let mut custom_feed =
                query_region_iter(&iupac_bam_path, "chr1", Pos0(90), Pos0(300), &config).unwrap();
            let mut reused_rec = SamRecord::new();
            for (idx, exp) in expected_custom.iter().enumerate() {
                assert!(custom_feed.next_into(&mut reused_rec).unwrap());
                assert_eq!(&reused_rec, exp, "custom IUPAC/QUAL mismatch at #{idx}");
            }
            assert!(!custom_feed.next_into(&mut reused_rec).unwrap());

            // Explicit checks on decoded sequence, empty QUAL = *, and Phred+33
            assert_eq!(expected_custom[0].sequence, b"=ACMGRSVTWYHKDBN");
            assert_eq!(
                expected_custom[0].quality_scores,
                vec![33, 34, 43, 53, 63, 66, 73, 83, 93, 103, 113, 118, 123, 124, 125, 126]
            );
            assert_eq!(expected_custom[1].sequence, b"=ACMGRSVTWYHKDBN=");
            assert!(
                expected_custom[1].quality_scores.is_empty(),
                "QUAL = * must decode to empty quality_scores"
            );
            assert_eq!(expected_custom[2].sequence, b"=ACMGRSVTWYHKDB");
            assert_eq!(
                expected_custom[2].quality_scores,
                vec![126, 125, 124, 123, 93, 78, 68, 33, 34, 65, 66, 73, 83, 108, 126]
            );
            for (idx, &nib_ch) in EXPECTED_NIBBLES.iter().enumerate() {
                assert_eq!(expected_custom[3 + idx].sequence, vec![nib_ch]);
                if idx % 2 == 0 {
                    assert!(expected_custom[3 + idx].quality_scores.is_empty());
                } else {
                    assert_eq!(expected_custom[3 + idx].quality_scores, vec![126]);
                }
            }

            // Raw BAM record with out-of-spec Phred bytes (221..=255, non-uniform so not QUAL=*)
            // must saturate cleanly at 255 in both fill_from_bam_feed_record and fill_from_bam_record.
            let sat_template = mk_custom_rec("sat_raw_qual", 100, b"ACGTACGT", Some(&[10u8; 8]));
            let mut raw_encoded = Vec::new();
            noodles::bam::io::Writer::from(&mut raw_encoded)
                .write_alignment_record(&header, &sat_template)
                .unwrap();
            let q_offset = raw_encoded.len() - 8;
            raw_encoded[q_offset..].copy_from_slice(&[0, 93, 221, 222, 223, 250, 254, 255]);
            let mut raw_bam_rec = noodles::bam::Record::default();
            assert!(
                noodles::bam::io::Reader::from(&raw_encoded[..])
                    .read_record(&mut raw_bam_rec)
                    .unwrap()
                    > 0
            );
            let contig_names: Vec<Box<str>> = vec!["chr1".into()];
            let mut cigar_scratch = vec![(8, CigarKind::Match)];
            let mut sa_buf = String::new();
            let mut mc_buf = String::new();
            let mut mate_ref_buf = String::new();
            let mut sat_feed_out = SamRecord::new();
            sat_feed_out.fill_from_bam_feed_record(
                &raw_bam_rec,
                &contig_names,
                0,
                Position::new(100),
                Position::new(107),
                &mut cigar_scratch,
                &mut sa_buf,
                &mut mc_buf,
                &mut mate_ref_buf,
            );
            let mut sat_bam_out = SamRecord::new();
            sat_bam_out
                .fill_from_bam_record(&raw_bam_rec, &header)
                .unwrap();
            assert_eq!(
                sat_feed_out.quality_scores,
                vec![33, 126, 254, 255, 255, 255, 255, 255]
            );
            assert_eq!(sat_bam_out, sat_feed_out);

            // Transition from non-empty saturating qualities to empty QUAL = * on the same reused SamRecord:
            let star_template = mk_custom_rec("star_raw_qual", 100, b"ACGTACG", None);
            let mut star_encoded = Vec::new();
            noodles::bam::io::Writer::from(&mut star_encoded)
                .write_alignment_record(&header, &star_template)
                .unwrap();
            noodles::bam::io::Reader::from(&star_encoded[..])
                .read_record(&mut raw_bam_rec)
                .unwrap();
            let mut star_cigar = vec![(7, CigarKind::Match)];
            sat_feed_out.fill_from_bam_feed_record(
                &raw_bam_rec,
                &contig_names,
                0,
                Position::new(100),
                Position::new(106),
                &mut star_cigar,
                &mut sa_buf,
                &mut mc_buf,
                &mut mate_ref_buf,
            );
            assert_eq!(sat_feed_out.sequence, b"ACGTACG");
            assert!(sat_feed_out.quality_scores.is_empty());
        }
    }

    #[test]
    fn preprocess_record_drops_i32_overflowing_reference_cigar_span() {
        let config = Config::default();
        let mut prep = RecordPreprocessor::new();

        // 1. Exact hostile record from brief: Pos1(100) + 2_000_000_000N + 500_000_000D (+ 4M)
        let hostile_skip_del = SamRecord {
            name: b"hostile_overflow_nd".to_vec(),
            flags: 0,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(100)),
            alignment_end: Some(Pos1(103)),
            mapping_quality: Some(60),
            cigar: vec![
                (2_000_000_000, CigarKind::Skip),
                (500_000_000, CigarKind::Deletion),
                (4, CigarKind::Match),
            ],
            sequence: b"ACGT".to_vec(),
            quality_scores: vec![40; 4],
            ..SamRecord::default()
        };
        assert!(
            !prep.preprocess_record(&hostile_skip_del, &config),
            "record whose CIGAR reference span overflows i32::MAX must be dropped"
        );
        assert_eq!(prep.total_reads, 0);

        // 2. ref_consumption fits in i32::MAX, but (start - 1) + ref_consumption > i32::MAX
        let hostile_end_overflow = SamRecord {
            name: b"hostile_end_overflow".to_vec(),
            alignment_start: Some(Pos1(1_800_000_000)),
            cigar: vec![(500_000_000, CigarKind::Skip), (4, CigarKind::Match)],
            ..hostile_skip_del.clone()
        };
        assert!(
            !prep.preprocess_record(&hostile_end_overflow, &config),
            "record whose 1-based reference end exceeds i32::MAX must be dropped"
        );
        assert_eq!(prep.total_reads, 0);

        // 3. Exact i32::MAX boundary (100 - 1 + (i32::MAX - 99) == i32::MAX) is accepted,
        //    while +1 base past i32::MAX is rejected.
        let exact_skip = (i32::MAX as u32) - 99 - 4;
        let boundary_ok = SamRecord {
            name: b"boundary_ok".to_vec(),
            alignment_start: Some(Pos1(100)),
            cigar: vec![(exact_skip, CigarKind::Skip), (4, CigarKind::Match)],
            ..hostile_skip_del.clone()
        };
        let boundary_over = SamRecord {
            name: b"boundary_over".to_vec(),
            alignment_start: Some(Pos1(100)),
            cigar: vec![(exact_skip + 1, CigarKind::Skip), (4, CigarKind::Match)],
            ..hostile_skip_del.clone()
        };
        assert!(prep.preprocess_record(&boundary_ok, &config));
        assert_eq!(prep.total_reads, 1);
        assert!(!prep.preprocess_record(&boundary_over, &config));
        assert_eq!(prep.total_reads, 1);

        // Zero-reference-consumption CIGARs (`50S`, `50I`) at `Pos1(i32::MAX as u32 + 1)`
        // must also be rejected before `alignment_start.0 as i32` wraps negative.
        let zero_span_soft = SamRecord {
            name: b"zero_span_soft_past_i32_max".to_vec(),
            alignment_start: Some(Pos1(i32::MAX as u32 + 1)),
            alignment_end: Some(Pos1(i32::MAX as u32 + 1)),
            cigar: vec![(50, CigarKind::SoftClip)],
            sequence: vec![b'A'; 50],
            quality_scores: vec![40; 50],
            ..hostile_skip_del.clone()
        };
        let zero_span_ins = SamRecord {
            name: b"zero_span_ins_past_i32_max".to_vec(),
            cigar: vec![(50, CigarKind::Insertion)],
            ..zero_span_soft.clone()
        };
        assert!(
            !prep.preprocess_record(&zero_span_soft, &config),
            "50S at Pos1(i32::MAX as u32 + 1) must be rejected"
        );
        assert!(
            !prep.preprocess_record(&zero_span_ins, &config),
            "50I at Pos1(i32::MAX as u32 + 1) must be rejected"
        );
        assert_eq!(prep.total_reads, 1);

        // 4. get_aligned_length_md saturates at i32::MAX instead of panicking or wrapping negative
        assert_eq!(
            get_aligned_length_md(&[
                (2_000_000_000, CigarKind::Match),
                (500_000_000, CigarKind::Deletion),
            ]),
            i32::MAX
        );
        assert_eq!(
            get_aligned_length_md(&[(u32::MAX, CigarKind::Deletion)]),
            i32::MAX
        );
        assert_eq!(
            get_aligned_length_md(&[
                (2_000_000_000, CigarKind::Skip),
                (50, CigarKind::Match),
                (10, CigarKind::Deletion),
                (20, CigarKind::Insertion),
            ]),
            60
        );
    }

    #[test]
    fn build_noodles_region_and_query_region_support_hla_colon_contigs() {
        use noodles::bam::bai;
        use noodles::core::Position;
        use noodles::csi::binning_index::{
            index::reference_sequence::bin::Chunk, Indexer as BinningIndexer,
        };
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;

        let contig = "HLA-A*01:01:01:01";
        let reg = build_noodles_region(contig, 100, 200).unwrap();
        assert_eq!(reg.name(), contig.as_bytes());
        assert!(build_noodles_region("", 1, 10).is_err());
        assert!(build_noodles_region(contig, 0, 10).is_err());
        assert!(build_noodles_region(contig, 20, 10).is_err());

        let dir = tempfile::tempdir().unwrap();
        let bam_path = dir.path().join("hla.bam");
        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(3503usize.try_into().unwrap())
            .build()
            .unwrap();
        let header = Header::builder().add_reference_sequence(contig, rs).build();
        let rec = RecordBuf::builder()
            .set_name("hla_read_1")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(120usize).unwrap())
            .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
            .set_cigar([Op::new(Kind::Match, 8)].into_iter().collect::<Cigar>())
            .set_sequence(b"ACGTACGT".to_vec().into())
            .set_quality_scores(vec![35u8; 8].into())
            .build();

        {
            let mut w = noodles::bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            w.write_header(&header).unwrap();
            let s_vp = w.get_ref().virtual_position();
            w.write_alignment_record(&header, &rec).unwrap();
            let e_vp = w.get_ref().virtual_position();
            let mut indexer = BinningIndexer::default();
            indexer
                .add_record(
                    Some((
                        0,
                        Position::try_from(120usize).unwrap(),
                        Position::try_from(127usize).unwrap(),
                        true,
                    )),
                    Chunk::new(s_vp, e_vp),
                )
                .unwrap();
            w.try_finish().unwrap();
            bai::io::Writer::new(
                std::fs::File::create(bam_path.with_extension("bam.bai")).unwrap(),
            )
            .write_index(&indexer.build(1))
            .unwrap();
        }

        let cfg = Config::default();
        let via_vec = query_region(&bam_path, contig, Pos0(99), Pos0(200), &cfg).unwrap();
        assert_eq!(via_vec.len(), 1);
        assert_eq!(via_vec[0].reference_name.as_deref(), Some(contig));
        assert_eq!(via_vec[0].alignment_start, Some(Pos1(120)));

        let via_iter: Vec<SamRecord> =
            query_region_iter(&bam_path, contig, Pos0(99), Pos0(200), &cfg)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
        assert_eq!(via_iter, via_vec);
    }

    #[test]
    fn cigar_span_and_aux_overflow_return_invalid_data_or_clamp() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{
            data::field::Value as BufValue, Cigar, Data, RecordBuf,
        };
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;

        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(1000usize.try_into().unwrap())
            .build()
            .unwrap();
        let header = Header::builder().add_reference_sequence("chr1", rs).build();

        // 1. checked_cigar_op_len & checked_alignment_end return InvalidData on overflow
        assert_eq!(
            checked_alignment_end(Position::new(usize::MAX - 4).unwrap(), 10)
                .unwrap_err()
                .kind(),
            std::io::ErrorKind::InvalidData
        );
        #[cfg(target_pointer_width = "64")]
        {
            let huge_op_rec = RecordBuf::builder()
                .set_name("huge_op")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(10usize).unwrap())
                .set_cigar(
                    [Op::new(Kind::Match, (u32::MAX as usize) + 1)]
                        .into_iter()
                        .collect::<Cigar>(),
                )
                .set_sequence(b"ACGT".to_vec().into())
                .set_quality_scores(vec![30u8; 4].into())
                .build();
            let err = SamRecord::from_alignment_record(&huge_op_rec, &header).unwrap_err();
            assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        }

        // 2. alignment_end overflow with u32-valid CIGAR op (never saturates to usize::MAX or None)
        let end_overflow_rec = RecordBuf::builder()
            .set_name("end_overflow")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(usize::MAX - 2).unwrap())
            .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect::<Cigar>())
            .set_sequence(b"ACGTACGTAC".to_vec().into())
            .set_quality_scores(vec![30u8; 10].into())
            .build();
        let err = SamRecord::from_alignment_record(&end_overflow_rec, &header).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);

        // 3. Multi-op alignment_span usize overflow with u32-valid ops fails via checked_add
        #[cfg(target_pointer_width = "32")]
        {
            let span_overflow_rec = RecordBuf::builder()
                .set_name("span_overflow_32")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(10usize).unwrap())
                .set_cigar(
                    [
                        Op::new(Kind::Skip, usize::MAX - 2),
                        Op::new(Kind::Match, 10),
                    ]
                    .into_iter()
                    .collect::<Cigar>(),
                )
                .set_sequence(b"ACGTACGTAC".to_vec().into())
                .set_quality_scores(vec![30u8; 10].into())
                .build();
            assert_eq!(
                SamRecord::from_alignment_record(&span_overflow_rec, &header)
                    .unwrap_err()
                    .kind(),
                std::io::ErrorKind::InvalidData
            );
        }

        // 4. NM / nM UInt32 > i32::MAX clamps to i32::MAX; MQ out of 0..=255 maps to None
        let aux_data: Data = [
            (Tag::EDIT_DISTANCE, BufValue::UInt32(u32::MAX)),
            (Tag::new(b'M', b'Q'), BufValue::Int32(300)),
        ]
        .into_iter()
        .collect();
        let aux_rec = RecordBuf::builder()
            .set_name("aux_bounds")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(50usize).unwrap())
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect::<Cigar>())
            .set_sequence(b"ACGT".to_vec().into())
            .set_quality_scores(vec![30u8; 4].into())
            .set_data(aux_data)
            .build();
        let parsed_align = SamRecord::from_alignment_record(&aux_rec, &header).unwrap();
        assert_eq!(parsed_align.nm, Some(i32::MAX));
        assert_eq!(parsed_align.mq, None);

        let mut encoded = Vec::new();
        noodles::bam::io::Writer::from(&mut encoded)
            .write_alignment_record(&header, &aux_rec)
            .unwrap();
        let mut raw_bam = noodles::bam::Record::default();
        noodles::bam::io::Reader::from(&encoded[..])
            .read_record(&mut raw_bam)
            .unwrap();
        let mut parsed_bam = SamRecord::new();
        parsed_bam.fill_from_bam_record(&raw_bam, &header).unwrap();
        assert_eq!(parsed_bam.nm, Some(i32::MAX));
        assert_eq!(parsed_bam.mq, None);
    }

    #[test]
    fn preprocess_record_rejects_pnext_above_i32_max_even_without_dedup() {
        let mut config = Config::default();
        config.dedup = false;
        let mut prep = RecordPreprocessor::new();

        let base = SamRecord {
            name: b"mate_pos_test".to_vec(),
            flags: 0x1,
            reference_name: Some("chr1".to_string()),
            alignment_start: Some(Pos1(100)),
            alignment_end: Some(Pos1(103)),
            mapping_quality: Some(60),
            cigar: vec![(4, CigarKind::Match)],
            sequence: b"ACGT".to_vec(),
            quality_scores: vec![40; 4],
            mate_reference_name: Some("chr1".to_string()),
            mate_alignment_start: Some(Pos1(i32::MAX as u32 + 1)),
            ..SamRecord::default()
        };

        // Rejected when dedup == false
        assert!(!prep.preprocess_record(&base, &config));
        assert_eq!(prep.total_reads, 0);

        // Exact i32::MAX mate_alignment_start is accepted
        let ok_rec = SamRecord {
            mate_alignment_start: Some(Pos1(i32::MAX as u32)),
            ..base.clone()
        };
        assert!(prep.preprocess_record(&ok_rec, &config));
        assert_eq!(prep.total_reads, 1);

        // total_reads and duplicate_reads saturate at u32::MAX without debug panic
        config.dedup = true;
        prep.total_reads = u32::MAX;
        prep.duplicate_reads = u32::MAX;
        let dup_rec = SamRecord {
            mate_alignment_start: Some(Pos1(5)),
            ..base
        };
        assert!(prep.preprocess_record(&dup_rec, &config));
        assert_eq!(prep.total_reads, u32::MAX);
        assert!(!prep.preprocess_record(&dup_rec, &config));
        assert_eq!(prep.duplicate_reads, u32::MAX);
    }

    #[test]
    fn companion_and_symlink_index_fallbacks_corrupt_primary_and_cram_adapter_errors() {
        use noodles::fasta::repository::Adapter as _;

        let dir = tempfile::tempdir().unwrap();
        let config = Config::default();

        let micro_bam = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/micro/micro.bam"
        ));
        let micro_bai = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/micro/micro.bam.bai"
        ));
        let baseline: Vec<SamRecord> =
            query_region_iter(micro_bam, "chr1", Pos0(0), Pos0(300), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert!(!baseline.is_empty());

        // 1. Picard companion `<stem>.bai` (`sample.bam` + `sample.bai`) succeeds in both `query_region` and `query_region_iter`
        let comp_dir = dir.path().join("companion");
        std::fs::create_dir_all(&comp_dir).unwrap();
        let sample_bam = comp_dir.join("sample.bam");
        let sample_bai = comp_dir.join("sample.bai");
        std::fs::copy(micro_bam, &sample_bam).unwrap();
        std::fs::copy(micro_bai, &sample_bai).unwrap();

        assert_eq!(
            query_region(&sample_bam, "chr1", Pos0(0), Pos0(300), &config).unwrap(),
            baseline
        );
        let comp_recs: Vec<SamRecord> =
            query_region_iter(&sample_bam, "chr1", Pos0(0), Pos0(300), &config)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert_eq!(comp_recs, baseline);

        // 2. Relative symlink `workdir/link.bam -> ../realdir/micro.bam` (index only in `realdir`)
        #[cfg(unix)]
        {
            let real_dir = dir.path().join("realdir");
            let work_dir = dir.path().join("workdir");
            std::fs::create_dir_all(&real_dir).unwrap();
            std::fs::create_dir_all(&work_dir).unwrap();
            let real_bam = real_dir.join("micro.bam");
            let real_bai = real_dir.join("micro.bam.bai");
            std::fs::copy(micro_bam, &real_bam).unwrap();
            std::fs::copy(micro_bai, &real_bai).unwrap();

            let link_bam = work_dir.join("link.bam");
            std::os::unix::fs::symlink(Path::new("../realdir/micro.bam"), &link_bam).unwrap();
            assert_eq!(
                query_region(&link_bam, "chr1", Pos0(0), Pos0(300), &config).unwrap(),
                baseline
            );
            let link_recs: Vec<SamRecord> =
                query_region_iter(&link_bam, "chr1", Pos0(0), Pos0(300), &config)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            assert_eq!(link_recs, baseline);
        }

        // 3. Corrupt or dangling-symlink sidecars while valid companion/canonical index exists MUST propagate error (never fall back)
        let bad_dir = dir.path().join("corrupt_primary");
        std::fs::create_dir_all(&bad_dir).unwrap();
        let bad_bam = bad_dir.join("sample.bam");
        std::fs::copy(micro_bam, &bad_bam).unwrap();
        std::fs::copy(micro_bai, bad_dir.join("sample.bai")).unwrap();
        std::fs::write(bad_dir.join("sample.bam.bai"), b"NOT_A_VALID_BAI_INDEX").unwrap();
        assert!(query_region(&bad_bam, "chr1", Pos0(0), Pos0(300), &config).is_err());
        assert!(query_region_iter(&bad_bam, "chr1", Pos0(0), Pos0(300), &config).is_err());

        #[cfg(unix)]
        {
            let dang_dir = dir.path().join("dangling_sidecars");
            std::fs::create_dir_all(&dang_dir).unwrap();
            let dang_bam = dang_dir.join("sample.bam");
            std::fs::copy(micro_bam, &dang_bam).unwrap();
            std::fs::copy(micro_bai, dang_dir.join("sample.bai")).unwrap();

            // Dangling primary `.bam.bai` and `.bam.csi` symlinks must NOT fall back to valid `sample.bai`
            let primary_bai_link = dang_dir.join("sample.bam.bai");
            std::os::unix::fs::symlink(Path::new("missing_target.bai"), &primary_bai_link).unwrap();
            assert!(query_region(&dang_bam, "chr1", Pos0(0), Pos0(300), &config).is_err());
            assert!(query_region_iter(&dang_bam, "chr1", Pos0(0), Pos0(300), &config).is_err());
            std::fs::remove_file(&primary_bai_link).unwrap();

            let primary_csi_link = dang_dir.join("sample.bam.csi");
            std::os::unix::fs::symlink(Path::new("missing_target.csi"), &primary_csi_link).unwrap();
            assert!(query_region(&dang_bam, "chr1", Pos0(0), Pos0(300), &config).is_err());
            assert!(query_region_iter(&dang_bam, "chr1", Pos0(0), Pos0(300), &config).is_err());
            std::fs::remove_file(&primary_csi_link).unwrap();

            // Dangling companion `link.bai` on symlinked BAM `link.bam -> ../companion/sample.bam` must NOT fall back to canonical index
            let link_dang_dir = dir.path().join("link_dangling");
            std::fs::create_dir_all(&link_dang_dir).unwrap();
            let link_dang_bam = link_dang_dir.join("link.bam");
            std::os::unix::fs::symlink(Path::new("../companion/sample.bam"), &link_dang_bam)
                .unwrap();
            let stem_bai_link = link_dang_dir.join("link.bai");
            std::os::unix::fs::symlink(Path::new("missing_stem.bai"), &stem_bai_link).unwrap();
            assert!(read_bam_index_with_fallbacks(&link_dang_bam).is_err());
            std::fs::remove_file(&stem_bai_link).unwrap();

            let stem_csi_link = link_dang_dir.join("link.csi");
            std::os::unix::fs::symlink(Path::new("missing_stem.csi"), &stem_csi_link).unwrap();
            assert!(read_bam_index_with_fallbacks(&link_dang_bam).is_err());
            std::fs::remove_file(&stem_csi_link).unwrap();

            // Dangling CRAM `.cram.crai` and `.crai` symlinks must fail with `(err, false)` instead of falling back
            let cram_real_dir = dir.path().join("cram_real");
            let cram_link_dir = dir.path().join("cram_link");
            std::fs::create_dir_all(&cram_real_dir).unwrap();
            std::fs::create_dir_all(&cram_link_dir).unwrap();
            let real_cram = cram_real_dir.join("sample.cram");
            std::fs::write(&real_cram, b"").unwrap();
            noodles::cram::crai::fs::write(
                cram_real_dir.join("sample.cram.crai"),
                &noodles::cram::crai::Index::default(),
            )
            .unwrap();
            assert!(read_cram_index_with_fallbacks(&real_cram).is_ok());

            let link_cram = cram_link_dir.join("link.cram");
            std::os::unix::fs::symlink(Path::new("../cram_real/sample.cram"), &link_cram).unwrap();
            assert!(read_cram_index_with_fallbacks(&link_cram).is_ok());

            let primary_crai_link = cram_link_dir.join("link.cram.crai");
            std::os::unix::fs::symlink(Path::new("missing.crai"), &primary_crai_link).unwrap();
            assert!(read_cram_index_with_fallbacks(&link_cram).is_err());
            std::fs::remove_file(&primary_crai_link).unwrap();

            let stem_crai_link = cram_link_dir.join("link.crai");
            std::os::unix::fs::symlink(Path::new("missing_stem.crai"), &stem_crai_link).unwrap();
            assert!(read_cram_index_with_fallbacks(&link_cram).is_err());
        }

        // 4. `RegistryDelegatingAdapter::get` on missing contig, empty contig, or dropped Weak returns Some(Err(NotFound)), NEVER None
        let ref_fa =
            write_test_fasta_with_fai(dir.path(), "adapter_ref.fa", &[("chr1", b"ACGTACGT")]);
        let adapter = PreallocatedFastaAdapter::open(&ref_fa).unwrap();
        let shared_repo =
            SharedFastaRepo::new(ref_fa.clone(), adapter, DEFAULT_FASTA_CACHE_BUDGET_BYTES);

        let mut live_delegating = RegistryDelegatingAdapter {
            primary_name: Vec::new(),
            primary_record: None,
            shared: std::sync::Arc::downgrade(&shared_repo),
        };
        let missing_res = live_delegating.get(b"chr_missing");
        match missing_res {
            Some(Err(err)) => assert_eq!(err.kind(), std::io::ErrorKind::NotFound),
            other => panic!("expected Some(Err(NotFound)) for missing contig, got {other:?}"),
        }
        let empty_res = live_delegating.get(b"");
        match empty_res {
            Some(Err(err)) => assert_eq!(err.kind(), std::io::ErrorKind::NotFound),
            other => panic!("expected Some(Err(NotFound)) for empty contig, got {other:?}"),
        }

        drop(shared_repo);
        let dropped_res = live_delegating.get(b"chr1");
        match dropped_res {
            Some(Err(err)) => assert_eq!(err.kind(), std::io::ErrorKind::NotFound),
            other => panic!("expected Some(Err(NotFound)) for dropped weak handle, got {other:?}"),
        }

        // 5. `PreallocatedFastaAdapter::get` rejects `.fa.gz` decompression bomb exceeding `.fai` length + CRAM_REF_OVERHANG_PAD
        let bomb_fa = dir.path().join("bomb.fa.gz");
        let bomb_fai = dir.path().join("bomb.fa.gz.fai");
        let bomb_gzi = dir.path().join("bomb.fa.gz.gzi");
        {
            use std::io::Write as _;
            let mut w = noodles_bgzf::io::Writer::new(std::fs::File::create(&bomb_fa).unwrap());
            w.write_all(b">chr1\n").unwrap();
            let big_line = vec![b'A'; CRAM_REF_OVERHANG_PAD + 2048];
            w.write_all(&big_line).unwrap();
            w.write_all(b"\n").unwrap();
            w.finish().unwrap();
            noodles_bgzf::gzi::fs::write(&bomb_gzi, &noodles_bgzf::gzi::Index::default()).unwrap();
            // .fai claims only 16 bp while BGZF payload has > 65_536 + 16 bp
            std::fs::write(&bomb_fai, b"chr1\t16\t6\t16\t17\n").unwrap();
        }
        let mut bomb_adapter = PreallocatedFastaAdapter::open(&bomb_fa).unwrap();
        match bomb_adapter.get(b"chr1") {
            Some(Err(err)) => assert_eq!(err.kind(), std::io::ErrorKind::InvalidData),
            other => panic!("expected Some(Err(InvalidData)) on .fa.gz bomb, got {other:?}"),
        }
    }

    #[test]
    fn batch12_cram_crai_bounds_spoofed_headers_and_landmark_panic_guards() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::MappingQuality;
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use noodles::sam::header::record::value::Map;
        use std::io::Write as _;
        use std::num::NonZeroUsize;

        let dir = tempfile::tempdir().unwrap();

        // 1. PreallocatedFastaAdapter::open rejects malformed .fai (lbc > lw)
        let bad_fa = dir.path().join("bad_adapter.fa");
        std::fs::write(&bad_fa, b">chr1\nACGTACGT\n").unwrap();
        std::fs::write(dir.path().join("bad_adapter.fa.fai"), b"chr1\t8\t6\t9\t8\n").unwrap();
        assert!(
            PreallocatedFastaAdapter::open(&bad_fa).is_err(),
            "PreallocatedFastaAdapter::open must reject malformed .fai via read_validated_fai_index"
        );

        // 2. .crai gzip stream exceeding MAX_CRAI_BYTES (128 MiB) must be rejected
        let bomb_cram = dir.path().join("bomb_idx.cram");
        std::fs::write(&bomb_cram, b"").unwrap();
        let bomb_crai = dir.path().join("bomb_idx.cram.crai");
        {
            let mut gz = flate2::write::GzEncoder::new(
                std::fs::File::create(&bomb_crai).unwrap(),
                flate2::Compression::fast(),
            );
            let line = b"0\t1\t10\t100\t0\t50\n";
            let chunk = line.repeat(4096); // 65_536 bytes per chunk
            let target = 128 * 1024 * 1024 + 65_536;
            let mut written = 0usize;
            while written < target {
                gz.write_all(&chunk).unwrap();
                written += chunk.len();
            }
            gz.finish().unwrap();
        }
        let crai_err = read_cram_index_with_fallbacks(&bomb_cram)
            .expect_err(".crai exceeding 128 MiB decompressed must fail");
        assert_eq!(crai_err.kind(), std::io::ErrorKind::InvalidData);

        // Build a valid base CRAM with 1 read so we can craft hostile binary mutations
        let ref_fa = dir.path().join("cram_sec.fa");
        let ref_fai = dir.path().join("cram_sec.fa.fai");
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        std::fs::write(
            &ref_fa,
            format!(">chr1\n{}\n", std::str::from_utf8(seq).unwrap()),
        )
        .unwrap();
        std::fs::write(
            &ref_fai,
            format!("chr1\t{}\t6\t{}\t{}\n", seq.len(), seq.len(), seq.len() + 1),
        )
        .unwrap();
        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(seq.len()).unwrap()),
            )
            .build();
        let rec = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(5usize).unwrap())
            .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
            .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect::<Cigar>())
            .set_sequence(seq[4..14].to_vec().into())
            .set_quality_scores(vec![35u8; 10].into())
            .build();
        let base_cram = dir.path().join("base.cram");
        {
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(std::fs::File::open(&ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut w = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&base_cram)
                .unwrap();
            w.write_alignment_header(&header).unwrap();
            w.write_alignment_record(&header, &rec).unwrap();
            w.try_finish(&header).unwrap();
        }
        let base_crai = noodles::cram::fs::index(&base_cram).unwrap();
        let container_offset = base_crai[0].offset() as usize;
        let base_bytes = std::fs::read(&base_cram).unwrap();

        let cfg = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            ..Config::default()
        };

        // 3. Spoofed i32::MAX container length in 4-byte LE header must be rejected BEFORE read_container
        let huge_len_cram = dir.path().join("huge_len.cram");
        let mut huge_len_bytes = base_bytes.clone();
        huge_len_bytes[container_offset..container_offset + 4]
            .copy_from_slice(&i32::MAX.to_le_bytes());
        std::fs::write(&huge_len_cram, &huge_len_bytes).unwrap();
        noodles::cram::crai::fs::write(huge_len_cram.with_extension("cram.crai"), &base_crai)
            .unwrap();
        let huge_err: Result<Vec<SamRecord>, _> =
            query_region_iter(&huge_len_cram, "chr1", Pos0(0), Pos0(50), &cfg)
                .and_then(|it| it.collect());
        match huge_err {
            Err(crate::types::VardictError::Io(io_e)) => {
                assert_eq!(io_e.kind(), std::io::ErrorKind::InvalidData);
            }
            other => {
                panic!("expected InvalidData on spoofed i32::MAX CRAM container len, got {other:?}")
            }
        }

        // 4. Out-of-bounds landmark in CRAM container header (valid CRC32) must NOT panic in
        //    either query_region_iter or query_region!
        let oob_lm_cram = dir.path().join("oob_landmark.cram");
        {
            // Craft a minimal valid-CRC container header at `container_offset` with body len = 8
            // and landmarks = [0x3fff] (16_383 > 8), which panics `&self.src[..end]` if unguarded.
            let mut hdr_body = Vec::new();
            hdr_body.extend_from_slice(&8i32.to_le_bytes()); // length = 8
            hdr_body.push(0x00); // ref_seq_id = 0
            hdr_body.push(0x05); // alignment_start = 5
            hdr_body.push(0x0a); // alignment_span = 10
            hdr_body.push(0x01); // record_count = 1
            hdr_body.push(0x00); // record_counter = 0
            hdr_body.push(0x0a); // base_count = 10
            hdr_body.push(0x02); // block_count = 2
            hdr_body.push(0x01); // landmark count = 1
            hdr_body.extend_from_slice(&[0xbf, 0xff]); // ITF-8 landmark = 16_383 (> length 8)
            let mut crc = flate2::Crc::new();
            crc.update(&hdr_body);
            let crc_le = crc.sum().to_le_bytes();

            let mut crafted = base_bytes[..container_offset].to_vec();
            crafted.extend_from_slice(&hdr_body);
            crafted.extend_from_slice(&crc_le);
            crafted.extend_from_slice(&[0u8; 8]); // 8-byte container body
            std::fs::write(&oob_lm_cram, &crafted).unwrap();
            noodles::cram::crai::fs::write(oob_lm_cram.with_extension("cram.crai"), &base_crai)
                .unwrap();
        }
        let oob_feed_res: Result<Vec<SamRecord>, _> =
            query_region_iter(&oob_lm_cram, "chr1", Pos0(0), Pos0(50), &cfg)
                .and_then(|it| it.collect());
        match oob_feed_res {
            Err(crate::types::VardictError::Io(io_e)) => {
                assert_eq!(io_e.kind(), std::io::ErrorKind::InvalidData);
            }
            other => {
                panic!("expected InvalidData on OOB landmark CRAM via CramFeed, got {other:?}")
            }
        }
        let oob_vec_res = query_region(&oob_lm_cram, "chr1", Pos0(0), Pos0(50), &cfg);
        match oob_vec_res {
            Err(crate::types::VardictError::Io(io_e)) => {
                assert_eq!(io_e.kind(), std::io::ErrorKind::InvalidData);
            }
            other => {
                panic!("expected InvalidData on OOB landmark CRAM via query_region, got {other:?}")
            }
        }

        // 5. Spoofed slice header record_count > MAX_CRAM_CONTAINER_RECORDS (with valid container
        //    record_count = 1) must be rejected BEFORE slice.records allocates vec![Record; count]!
        let spoof_slice_cram = dir.path().join("spoof_slice_rc.cram");
        {
            // Build a container with valid CompressionHeader from `base_bytes` and a crafted
            // SliceHeader block claiming record_count = 5_000_000 (ITF-8: [0xe0, 0x4c, 0x4b, 0x40]).
            let mut rdr = noodles::cram::io::Reader::new(std::fs::File::open(&base_cram).unwrap());
            let _ = rdr.read_header().unwrap();
            let c_pos = rdr.position().unwrap() as usize;
            let mut orig_cont = noodles::cram::io::reader::Container::default();
            let c_len = rdr.read_container(&mut orig_cont).unwrap();
            let body_start =
                std::io::Seek::stream_position(rdr.get_mut()).unwrap() as usize - c_len;
            let lm0 = orig_cont.header().landmarks()[0];
            let comp_hdr_bytes = &base_bytes[body_start..body_start + lm0];

            // SliceHeader uncompressed payload:
            let mut sh_payload = vec![
                0x00, // ref_seq_id = 0
                0x05, // alignment_start = 5
                0x0a, // alignment_span = 10
                0xe0, 0x4c, 0x4b, 0x40, // record_count = 5_000_000
                0x00, // record_counter = 0
                0x01, // block_count = 1
                0x00, // block_content_ids len = 0
                0xff, 0xff, 0xff, 0xff, 0x0f, // embedded_ref = -1
            ];
            sh_payload.extend_from_slice(&[0u8; 16]); // MD5 = 0 (None)

            let mut sh_block = vec![
                0x00,                   // method = None
                0x02,                   // content_type = SliceHeader (2)
                0x00,                   // content_id = 0
                sh_payload.len() as u8, // compressed_size (< 128)
                sh_payload.len() as u8, // uncompressed_size (< 128)
            ];
            sh_block.extend_from_slice(&sh_payload);
            let mut blk_crc = flate2::Crc::new();
            blk_crc.update(&sh_block);
            sh_block.extend_from_slice(&blk_crc.sum().to_le_bytes());

            // Add an empty CoreData block (method=0, content_type=1, content_id=0, size=0, raw=0)
            let mut core_blk = vec![0x00, 0x01, 0x00, 0x00, 0x00];
            let mut c_crc = flate2::Crc::new();
            c_crc.update(&core_blk);
            core_blk.extend_from_slice(&c_crc.sum().to_le_bytes());

            let mut cont_body = Vec::new();
            cont_body.extend_from_slice(comp_hdr_bytes);
            let new_lm0 = cont_body.len();
            cont_body.extend_from_slice(&sh_block);
            cont_body.extend_from_slice(&core_blk);

            fn push_itf8(dst: &mut Vec<u8>, v: usize) {
                if v < 0x80 {
                    dst.push(v as u8);
                } else if v < 0x4000 {
                    dst.push(0x80 | ((v >> 8) as u8));
                    dst.push((v & 0xff) as u8);
                } else {
                    panic!("unexpectedly large test ITF-8 value: {v}");
                }
            }

            let mut hdr_body = Vec::new();
            hdr_body.extend_from_slice(&(cont_body.len() as i32).to_le_bytes());
            hdr_body.push(0x00); // ref_seq_id = 0
            hdr_body.push(0x05); // alignment_start = 5
            hdr_body.push(0x0a); // alignment_span = 10
            hdr_body.push(0x01); // container header record_count = 1 (spoofed low!)
            hdr_body.push(0x00); // record_counter = 0
            hdr_body.push(0x0a); // base_count = 10
            hdr_body.push(0x03); // block_count = 3
            hdr_body.push(0x01); // landmark count = 1
            push_itf8(&mut hdr_body, new_lm0);
            let mut h_crc = flate2::Crc::new();
            h_crc.update(&hdr_body);
            hdr_body.extend_from_slice(&h_crc.sum().to_le_bytes());

            let mut crafted = base_bytes[..c_pos].to_vec();
            crafted.extend_from_slice(&hdr_body);
            crafted.extend_from_slice(&cont_body);
            std::fs::write(&spoof_slice_cram, &crafted).unwrap();
            let crai_rec = vec![noodles::cram::crai::Record::new(
                Some(0),
                Position::try_from(5usize).ok(),
                10,
                c_pos as u64,
                new_lm0 as u64,
                (sh_block.len() + core_blk.len()) as u64,
            )];
            noodles::cram::crai::fs::write(spoof_slice_cram.with_extension("cram.crai"), &crai_rec)
                .unwrap();
        }
        let spoof_res: Result<Vec<SamRecord>, _> =
            query_region_iter(&spoof_slice_cram, "chr1", Pos0(0), Pos0(50), &cfg)
                .and_then(|it| it.collect());
        match spoof_res {
            Err(crate::types::VardictError::Io(io_e)) => {
                assert_eq!(
                    io_e.kind(),
                    std::io::ErrorKind::InvalidData,
                    "spoofed slice record_count must fail with InvalidData before UnexpectedEof, got {io_e:?}"
                );
            }
            other => panic!("expected InvalidData on spoofed slice record_count, got {other:?}"),
        }

        // 6. < 400 B crafted containers spoofing 2 GiB (0x7fff_ffff) in:
        //    (a) CoreData Gzip block `uncompressed_size = 0x7fff_ffff` + ExternalData Rans4x8, and
        //    (b) ExternalData Rans4x8 inner codec header `uncompressed_size = 0x7fff_ffff` (while
        //        outer block `uncompressed_size = 16` and CoreData is valid 1 B).
        //    Both `query_region` and `query_region_iter` must reject cleanly
        //    with `InvalidData` ("exceeds") before allocating (`SEC-06`).
        let write_spoof_codec_cram = |dst_cram: &Path,
                                      core_raw_block: Vec<u8>,
                                      ext_raw_block: Vec<u8>| {
            let mut rdr = noodles::cram::io::Reader::new(std::fs::File::open(&base_cram).unwrap());
            let _ = rdr.read_header().unwrap();
            let c_pos = rdr.position().unwrap() as usize;
            let mut orig_cont = noodles::cram::io::reader::Container::default();
            let c_len = rdr.read_container(&mut orig_cont).unwrap();
            let body_start =
                std::io::Seek::stream_position(rdr.get_mut()).unwrap() as usize - c_len;
            let lm0 = orig_cont.header().landmarks()[0];
            let comp_hdr_bytes = &base_bytes[body_start..body_start + lm0];

            let mut sh_payload = vec![
                0x00, // ref_seq_id = 0
                0x05, // alignment_start = 5
                0x0a, // alignment_span = 10
                0x01, // record_count = 1
                0x00, // record_counter = 0
                0x02, // block_count = 2 (CoreData + 1 ExternalData)
                0x01, // block_content_ids len = 1
                0x01, // block_content_ids[0] = 1
                0xff, 0xff, 0xff, 0xff, 0x0f, // embedded_ref = -1
            ];
            sh_payload.extend_from_slice(&[0u8; 16]);

            let mut sh_block = vec![
                0x00,                   // method = None
                0x02,                   // content_type = SliceHeader (2)
                0x00,                   // content_id = 0
                sh_payload.len() as u8, // compressed_size
                sh_payload.len() as u8, // uncompressed_size
            ];
            sh_block.extend_from_slice(&sh_payload);
            let mut sh_crc = flate2::Crc::new();
            sh_crc.update(&sh_block);
            sh_block.extend_from_slice(&sh_crc.sum().to_le_bytes());

            let seal_crc = |mut blk: Vec<u8>| -> Vec<u8> {
                let mut c = flate2::Crc::new();
                c.update(&blk);
                blk.extend_from_slice(&c.sum().to_le_bytes());
                blk
            };
            let core_blk = seal_crc(core_raw_block);
            let ext_blk = seal_crc(ext_raw_block);

            let mut cont_body = Vec::new();
            cont_body.extend_from_slice(comp_hdr_bytes);
            let new_lm0 = cont_body.len();
            cont_body.extend_from_slice(&sh_block);
            cont_body.extend_from_slice(&core_blk);
            cont_body.extend_from_slice(&ext_blk);
            assert!(
                cont_body.len() < 400,
                "crafted CRAM container body must be < 400 B, got {}",
                cont_body.len()
            );

            let mut hdr_body = Vec::new();
            hdr_body.extend_from_slice(&(cont_body.len() as i32).to_le_bytes());
            hdr_body.extend_from_slice(&[0x00, 0x05, 0x0a, 0x01, 0x00, 0x0a, 0x04, 0x01]);
            if new_lm0 < 0x80 {
                hdr_body.push(new_lm0 as u8);
            } else {
                hdr_body.push(0x80 | ((new_lm0 >> 8) as u8));
                hdr_body.push((new_lm0 & 0xff) as u8);
            }
            let mut h_crc = flate2::Crc::new();
            h_crc.update(&hdr_body);
            hdr_body.extend_from_slice(&h_crc.sum().to_le_bytes());

            let mut crafted = base_bytes[..c_pos].to_vec();
            crafted.extend_from_slice(&hdr_body);
            crafted.extend_from_slice(&cont_body);
            std::fs::write(dst_cram, &crafted).unwrap();
            let crai_rec = vec![noodles::cram::crai::Record::new(
                Some(0),
                Position::try_from(5usize).ok(),
                10,
                c_pos as u64,
                new_lm0 as u64,
                (sh_block.len() + core_blk.len() + ext_blk.len()) as u64,
            )];
            noodles::cram::crai::fs::write(dst_cram.with_extension("cram.crai"), &crai_rec)
                .unwrap();
        };

        // Rans4x8 ExternalData block: method=4, content_type=4, content_id=1,
        // compressed_size=9, uncompressed_size=16, payload=[order=0, comp=0u32, uncomp=0x7fff_ffff]
        let rans4x8_2gib_ext = vec![
            0x04, 0x04, 0x01, 0x09, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0x7f,
        ];
        // Gzip CoreData block: method=1, content_type=5, content_id=0,
        // compressed_size=2, uncompressed_size=0x7fff_ffff (ITF-8: f7 ff ff ff 0f), payload=[0x1f, 0x8b]
        let gzip_2gib_core = vec![
            0x01, 0x05, 0x00, 0x02, 0xf7, 0xff, 0xff, 0xff, 0x0f, 0x1f, 0x8b,
        ];
        // Valid 1-byte raw CoreData block so ExternalData validation is reached:
        let valid_1b_core = vec![0x00, 0x05, 0x00, 0x01, 0x01, 0x00];

        // (c) < 400 B NameTokenizer 8.29 GiB + 0-stream bomb:
        // method=8, content_type=4, content_id=1, compressed_size=9, uncompressed_size=16,
        // payload=[inner_uncomp=16u32, name_count=2_000_000u32, sub_method=0]
        let nametok_8gib_ext = vec![
            0x08, 0x04, 0x01, 0x09, 0x10, 0x10, 0x00, 0x00, 0x00, 0x80, 0x84, 0x1e, 0x00, 0x00,
        ];

        // (d) < 400 B NameTokenizer stream-header malformed cases (name_count=1, inner_uncomp=4):
        // - `raw_ty = 12` (`Type::End`, exceeds writable TokenReader types `0..=9`)
        let nametok_bad_ty_ext = vec![
            0x08, 0x04, 0x01, 0x0d, 0x04, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x8c, 0x02, 0x20, 0x00,
        ];
        // - `tok_new | Diff` (synthesizes `type_reader`) followed by `tok_dup` referencing
        //   `DZLen` (`dup_ty = 4`, which panics with `unimplemented!` in `TokenReader::get`)
        let nametok_dup_dzlen_ext = vec![
            0x08, 0x04, 0x01, 0x14, 0x04, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
            0x86, 0x06, 0x20, 0x04, 0x00, 0x00, 0x00, 0x00, 0x46, 0x00, 0x04,
        ];

        let spoof_gzip_core_cram = dir.path().join("spoof_gzip_core_2gib.cram");
        write_spoof_codec_cram(
            &spoof_gzip_core_cram,
            gzip_2gib_core,
            rans4x8_2gib_ext.clone(),
        );

        let spoof_rans_ext_cram = dir.path().join("spoof_rans4x8_ext_2gib.cram");
        write_spoof_codec_cram(
            &spoof_rans_ext_cram,
            valid_1b_core.clone(),
            rans4x8_2gib_ext,
        );

        let spoof_nametok_8gib_cram = dir.path().join("spoof_nametok_8gib.cram");
        write_spoof_codec_cram(
            &spoof_nametok_8gib_cram,
            valid_1b_core.clone(),
            nametok_8gib_ext,
        );

        let spoof_nametok_bad_ty_cram = dir.path().join("spoof_nametok_bad_ty.cram");
        write_spoof_codec_cram(
            &spoof_nametok_bad_ty_cram,
            valid_1b_core.clone(),
            nametok_bad_ty_ext,
        );

        let spoof_nametok_dup_dzlen_cram = dir.path().join("spoof_nametok_dup_dzlen.cram");
        write_spoof_codec_cram(
            &spoof_nametok_dup_dzlen_cram,
            valid_1b_core,
            nametok_dup_dzlen_ext,
        );

        for (case_label, bad_cram, require_exceeds) in [
            ("gzip_core_2gib", &spoof_gzip_core_cram, true),
            ("rans4x8_ext_2gib", &spoof_rans_ext_cram, true),
            ("nametok_ext_8gib", &spoof_nametok_8gib_cram, true),
            ("nametok_bad_ty", &spoof_nametok_bad_ty_cram, true),
            ("nametok_dup_dzlen", &spoof_nametok_dup_dzlen_cram, false),
        ] {
            let feed_res: Result<Vec<SamRecord>, _> =
                query_region_iter(bad_cram, "chr1", Pos0(0), Pos0(50), &cfg)
                    .and_then(|it| it.collect());
            match feed_res {
                Err(crate::types::VardictError::Io(io_e)) => {
                    assert_eq!(io_e.kind(), std::io::ErrorKind::InvalidData, "{case_label}");
                    if require_exceeds {
                        assert!(
                            io_e.to_string().contains("exceeds"),
                            "{case_label} expected 'exceeds' in error, got: {io_e}"
                        );
                    }
                }
                other => panic!(
                    "{case_label}: expected InvalidData via query_region_iter, got {other:?}"
                ),
            }

            let vec_res = query_region(bad_cram, "chr1", Pos0(0), Pos0(50), &cfg);
            match vec_res {
                Err(crate::types::VardictError::Io(io_e)) => {
                    assert_eq!(io_e.kind(), std::io::ErrorKind::InvalidData, "{case_label}");
                    if require_exceeds {
                        assert!(
                            io_e.to_string().contains("exceeds"),
                            "{case_label} expected 'exceeds' in query_region error, got: {io_e}"
                        );
                    }
                }
                other => {
                    panic!("{case_label}: expected InvalidData via query_region, got {other:?}")
                }
            }
        }

        // Verify valid `noodles-cram` `NameTokenizer` ExternalData block passes prevalidation cleanly:
        let valid_nt_payload: &[u8] = &[
            0x58, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x06,
            0x00, 0x04, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00, 0x06, 0x18, 0x00, 0x0c, 0x00, 0x01, 0x00, 0x00, 0x0e, 0x02,
            0x00, 0x22, 0x25, 0x00, 0x00, 0xbc, 0x00, 0x00, 0x00, 0xbc, 0x00, 0x00, 0x00, 0xbc,
            0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x01, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02,
            0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x01,
            0x1b, 0x00, 0x04, 0x00, 0x31, 0x37, 0x49, 0x00, 0x01, 0x01, 0x01, 0x01, 0x00, 0x0c,
            0x02, 0x00, 0x00, 0x04, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00,
            0x80, 0x17, 0x00, 0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00,
            0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00,
            0x01, 0x5f, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x03, 0x0a, 0x00, 0x01,
            0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x03, 0x19, 0x00, 0x04, 0x00, 0x22, 0x3d, 0x00, 0x02, 0x01, 0x01,
            0x00, 0x0c, 0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
            0x01, 0x00, 0x04, 0x15, 0x00, 0x01, 0x05, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x17, 0x00,
            0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00,
            0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00,
            0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x07, 0x17, 0x00, 0x04, 0x00, 0x02, 0x00, 0x03, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00,
            0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x80, 0x17, 0x00,
            0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00,
            0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00,
            0x80, 0x00, 0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x03, 0x01, 0x00, 0xa8,
            0x00, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x07, 0x1a, 0x00, 0x08, 0x00, 0x7b, 0x7c, 0x00, 0x00, 0x06, 0x01, 0x01, 0x00, 0x7c,
            0x20, 0x00, 0x00, 0xe0, 0x00, 0x00, 0x00, 0xe0, 0x00, 0x00, 0x00, 0xe0, 0x00, 0x00,
            0x80, 0x17, 0x00, 0x03, 0x02, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00,
            0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00,
            0x01, 0x3a, 0x00, 0x01, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x07, 0x00, 0x04, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x07, 0x22, 0x00, 0x0c, 0x00, 0x06, 0x2d, 0x64, 0x65, 0x00, 0xb2, 0xf0, 0x00,
            0x0a, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0xcd, 0x0b, 0x08, 0x00, 0xaf, 0x0e,
            0x08, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x17, 0x00, 0x03, 0x02,
            0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x3a, 0x00, 0x01, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x80, 0x17, 0x00, 0x03, 0x03, 0x07, 0x00, 0x03, 0x01, 0x00, 0xa8, 0x00, 0x00,
            0x00, 0xa8, 0x00, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x80, 0x00, 0x00, 0x03, 0x1d,
            0x00, 0x08, 0x00, 0x06, 0x21, 0xa3, 0xe3, 0x00, 0x04, 0x01, 0x01, 0x01, 0x01, 0x00,
            0x6e, 0x20, 0x00, 0x00, 0x58, 0x20, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x02,
            0x00, 0x04, 0x15, 0x00, 0x02, 0x05, 0x00, 0x02, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x07, 0x19, 0x00, 0x04,
            0x00, 0x21, 0x3f, 0x00, 0x02, 0x01, 0x01, 0x00, 0x08, 0x02, 0x00, 0x00, 0x0c, 0x02,
            0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x17, 0x00, 0x03, 0x02,
            0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0xac,
            0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x02, 0x15, 0x00, 0x01, 0x23, 0x00, 0x01, 0x00,
            0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00,
            0x00, 0x80, 0x17, 0x00, 0x03, 0x07, 0x0a, 0x00, 0x01, 0x03, 0x00, 0x00, 0x02, 0x00,
            0x00, 0xac, 0x00, 0x00, 0x00, 0xac, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x07, 0x17,
            0x00, 0x04, 0x00, 0x09, 0x00, 0x03, 0x01, 0x00, 0x0c, 0x02, 0x00, 0x00, 0xa8, 0x00,
            0x00, 0x00, 0xa8, 0x00, 0x00, 0x00, 0xa8, 0x00, 0x00, 0x80, 0x15, 0x00, 0x03, 0x0c,
            0x00, 0x04, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00,
            0x00, 0x80, 0x00, 0x00,
        ];
        let clen = valid_nt_payload.len();
        let mut valid_nt_blk = vec![
            0x08, // method = NameTokenizer
            0x04, // content_type = ExternalData
            0x01, // content_id = 1
            0x80 | ((clen >> 8) as u8),
            (clen & 0xff) as u8, // 2-byte ITF-8 compressed_size
            0x58,                // uncompressed_size = 88
        ];
        valid_nt_blk.extend_from_slice(valid_nt_payload);
        let mut nt_crc = flate2::Crc::new();
        nt_crc.update(&valid_nt_blk);
        valid_nt_blk.extend_from_slice(&nt_crc.sum().to_le_bytes());
        let valid_nt_path = dir.path().join("valid_nametok_blk.bin");
        std::fs::write(&valid_nt_path, &valid_nt_blk).unwrap();
        let mut valid_nt_file = std::fs::File::open(&valid_nt_path).unwrap();
        let (end_off, peak_uncomp) = validate_cram_block_and_codec(
            &valid_nt_path,
            &mut valid_nt_file,
            0,
            valid_nt_blk.len() as u64,
            4,
            MAX_CRAM_CONTAINER_BYTES,
            3,
        )
        .expect("valid NameTokenizer block must pass prevalidation");
        assert_eq!(end_off, valid_nt_blk.len() as u64);
        assert!(peak_uncomp >= 88);
    }

    #[test]
    fn cram_container_cache_caches_standard_10k_illumina_and_bypasses_over_16mib() {
        use noodles::core::Position;
        use noodles::sam::{
            alignment::{
                io::Write as _,
                record::{
                    cigar::{op::Kind, Op},
                    Flags, MappingQuality,
                },
                record_buf::{Cigar, RecordBuf},
            },
            header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap},
            Header,
        };
        use std::fs::File;
        use std::num::NonZeroUsize;

        let chr1_arc: std::sync::Arc<str> = std::sync::Arc::from("chr1");
        let build_cached_block = |count: usize, read_len: usize| -> CachedCramBlock {
            let mut out = Vec::with_capacity(count);
            for idx in 0..count {
                let mut sam = SamRecord::new();
                sam.name = format!("read_{idx:05}").into_bytes();
                let mut cigar = Vec::with_capacity(4);
                cigar.push((read_len as u32, CigarKind::Match));
                sam.cigar = cigar;
                sam.sequence = vec![b'A'; read_len];
                sam.quality_scores = vec![35u8; read_len];
                out.push(CachedCramRecord {
                    reference_sequence_id: Some(0),
                    alignment_start: Position::new(100),
                    alignment_end: Position::new(100 + read_len - 1),
                    reference_name: Some(std::sync::Arc::clone(&chr1_arc)),
                    mate_reference_name: None,
                    sam,
                });
            }
            out.shrink_to_fit();
            std::sync::Arc::new(CachedCramContainer::from_records(out))
        };

        // Part A: Standard 10,000 × 150 bp Illumina CRAM container (~6.28 MiB > 4 MiB)
        // must be eligible under MAX_CACHED_CRAM_CONTAINER_BYTES (16 MiB), whereas an
        // oversized 240 × 40 kbp long-read container (~18.39 MiB > 16 MiB) must bypass.
        let illumina_10k = build_cached_block(10_000, 150);
        let illumina_bytes = approx_cram_container_bytes(&illumina_10k);
        assert!(
            illumina_bytes > 4 * 1024 * 1024,
            "10k × 150 bp Illumina container ({illumina_bytes} B) exceeds 4 MiB"
        );
        assert!(
            should_cache_cram_container(&illumina_10k, MAX_CACHED_CRAM_CONTAINER_BYTES),
            "10k × 150 bp Illumina container ({illumina_bytes} B) must be cached under MAX_CACHED_CRAM_CONTAINER_BYTES ({MAX_CACHED_CRAM_CONTAINER_BYTES} B)"
        );
        drop(illumina_10k);

        let long_read_240 = build_cached_block(240, 40_000);
        let long_bytes = approx_cram_container_bytes(&long_read_240);
        assert!(
            long_bytes > 16 * 1024 * 1024,
            "240 × 40 kbp long-read container ({long_bytes} B) must exceed 16 MiB"
        );
        assert!(
            !should_cache_cram_container(&long_read_240, MAX_CACHED_CRAM_CONTAINER_BYTES),
            "oversized long-read container ({long_bytes} B) must bypass CRAM_CONTAINER_CACHE"
        );
        drop(long_read_240);

        // Part B: Isolated end-to-end fetch_cram_container_in hit vs byte-cap bypass
        // on a compact 32 × 150 bp CRAM (zero shared global CRAM_CONTAINER_CACHE state).
        let temp_dir = tempfile::tempdir().unwrap();
        let ref_fa = temp_dir.path().join("short_ref.fa");
        let ref_fai = temp_dir.path().join("short_ref.fa.fai");
        let cram_path = temp_dir.path().join("compact_illumina.cram");

        let ref_len = 2_000usize;
        let ref_bases: Vec<u8> = (0..ref_len)
            .map(|i| b"ACGT"[(i.wrapping_mul(13).wrapping_add(i / 7)) % 4])
            .collect();
        std::fs::write(
            &ref_fa,
            format!(">chr1\n{}\n", std::str::from_utf8(&ref_bases).unwrap()),
        )
        .unwrap();
        std::fs::write(
            &ref_fai,
            format!("chr1\t{ref_len}\t6\t{ref_len}\t{}\n", ref_len + 1),
        )
        .unwrap();

        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(NonZeroUsize::new(ref_len).unwrap())
            .build()
            .unwrap();
        let header = Header::builder().add_reference_sequence("chr1", rs).build();

        let read_count = 32usize;
        let read_len = 150usize;
        let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
        let container_offset = {
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(File::open(&ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut cram_writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&cram_path)
                .unwrap();
            cram_writer.write_alignment_header(&header).unwrap();
            for idx in 0..read_count {
                let start_1 = 50 + idx * 10;
                let seq = ref_bases[(start_1 - 1)..(start_1 - 1 + read_len)].to_vec();
                let rec = RecordBuf::builder()
                    .set_name(format!("short_read_{idx:02}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(start_1).unwrap())
                    .set_mapping_quality(MappingQuality::new(60).unwrap())
                    .set_cigar(cigar.clone())
                    .set_sequence(seq.into())
                    .set_quality_scores(vec![35u8; read_len].into())
                    .build();
                cram_writer.write_alignment_record(&header, &rec).unwrap();
            }
            cram_writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&cram_path).unwrap();
            assert_eq!(
                crai.len(),
                1,
                "expected 32 reads in a single CRAM container"
            );
            let off = crai[0].offset();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
            off
        };

        let config = crate::config::Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            ..Default::default()
        };
        let meta = acquire_cram_meta(&cram_path, &config).unwrap();
        let pinned_slot = meta.fasta_repo.warm_contig_slot(b"chr1");

        // (i) Cache-eligible cap inserts into isolated shard (even when global BAM_CACHE_EPOCH
        //     advances mid-decode) and serves repeat call via Arc::ptr_eq
        let hit_cache: [std::sync::Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS] =
            std::array::from_fn(|_| std::sync::Mutex::new(CramContainerCacheShard::default()));
        TEST_MID_CRAM_DECODE_HOOK.with(|h| {
            *h.borrow_mut() = Some(Box::new(|| {
                BAM_CACHE_EPOCH.fetch_add(1, std::sync::atomic::Ordering::Release);
            }));
        });
        let block_1 = fetch_cram_container_in(
            &cram_path,
            &meta,
            container_offset,
            b"chr1",
            Some(&pinned_slot),
            &hit_cache,
            MAX_CACHED_CRAM_CONTAINER_BYTES,
        )
        .unwrap()
        .expect("non-empty CRAM container");
        assert_eq!(block_1.len(), read_count);
        let hit_entries: usize = hit_cache
            .iter()
            .map(|s| s.lock().unwrap().lru.entries.len())
            .sum();
        assert_eq!(
            hit_entries, 1,
            "eligible container must be cached in isolated shard despite global epoch bump"
        );

        let block_2 = fetch_cram_container_in(
            &cram_path,
            &meta,
            container_offset,
            b"chr1",
            Some(&pinned_slot),
            &hit_cache,
            MAX_CACHED_CRAM_CONTAINER_BYTES,
        )
        .unwrap()
        .expect("cached CRAM container");
        assert!(
            std::sync::Arc::ptr_eq(&block_1, &block_2),
            "second fetch must hit isolated CRAM container cache"
        );

        // (ii) Tight 1 KiB cap (< ~20 KiB container) bypasses cache insertion while decoding identically
        let bypass_cache: [std::sync::Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS] =
            std::array::from_fn(|_| std::sync::Mutex::new(CramContainerCacheShard::default()));
        let bypassed = fetch_cram_container_in(
            &cram_path,
            &meta,
            container_offset,
            b"chr1",
            Some(&pinned_slot),
            &bypass_cache,
            1024,
        )
        .unwrap()
        .expect("decoded CRAM container on bypass");
        let bypass_entries: usize = bypass_cache
            .iter()
            .map(|s| s.lock().unwrap().lru.entries.len())
            .sum();
        assert_eq!(bypass_entries, 0, "oversized container must bypass cache");
        assert_eq!(bypassed.len(), block_1.len());
        for (a, b) in block_1.iter().zip(bypassed.iter()) {
            assert_eq!(a.sam, b.sam);
        }

        // (iii) Global `CRAM_CONTAINER_CACHE` remains epoch-protected when `BAM_CACHE_EPOCH`
        //       advances mid-decode, and caches normally when `BAM_CACHE_EPOCH` is unchanged.
        //       Uses isolated `meta_global.id >= TEST_META_ID_BASE` so concurrent sibling
        //       `clear_bam_caches()` (`LruCacheShard::clear_production`) cannot evict `global_key`.
        {
            NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, std::sync::atomic::Ordering::Relaxed);
            let global_cram = temp_dir.path().join("global_epoch_guard.cram");
            std::fs::copy(&cram_path, &global_cram).unwrap();
            std::fs::copy(
                cram_path.with_extension("cram.crai"),
                global_cram.with_extension("cram.crai"),
            )
            .unwrap();
            let meta_global = acquire_cram_meta(&global_cram, &config).unwrap();
            assert!(meta_global.id >= TEST_META_ID_BASE);
            let global_key = (meta_global.id, container_offset);
            let global_shard_idx =
                cache_shard_index::<CRAM_CONTAINER_CACHE_SHARDS>(meta_global.id, container_offset);

            TEST_MID_CRAM_DECODE_HOOK.with(|h| {
                *h.borrow_mut() = Some(Box::new(|| {
                    BAM_CACHE_EPOCH.fetch_add(1, std::sync::atomic::Ordering::Release);
                }));
            });
            let uncached_on_race = fetch_cram_container(
                &global_cram,
                &meta_global,
                container_offset,
                b"chr1",
                Some(&pinned_slot),
            )
            .unwrap()
            .expect("decoded CRAM container across global epoch bump");
            assert_eq!(uncached_on_race.len(), read_count);
            let cached_after_race = cram_container_cache_shards()[global_shard_idx]
                .lock()
                .unwrap()
                .lru
                .get(global_key)
                .is_some();
            assert!(
                !cached_after_race,
                "global CRAM_CONTAINER_CACHE must skip insertion when BAM_CACHE_EPOCH advances mid-decode"
            );

            loop {
                let ep0 = BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire);
                let clean_global = fetch_cram_container(
                    &global_cram,
                    &meta_global,
                    container_offset,
                    b"chr1",
                    Some(&pinned_slot),
                )
                .unwrap()
                .expect("clean global CRAM container decode");
                if BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire) != ep0 {
                    continue;
                }
                assert_eq!(clean_global.len(), read_count);
                let hit = cram_container_cache_shards()[global_shard_idx]
                    .lock()
                    .unwrap()
                    .lru
                    .get(global_key);
                assert!(
                    hit.is_some_and(|b| std::sync::Arc::ptr_eq(&b, &clean_global)),
                    "global CRAM_CONTAINER_CACHE must insert when BAM_CACHE_EPOCH is stable"
                );
                break;
            }
            cram_container_cache_shards()[global_shard_idx]
                .lock()
                .unwrap()
                .remove_key(global_key);
        }
    }

    fn write_indexed_bam(bam_path: &Path, chr_len: usize, read_names: &[&str]) {
        use noodles::bam::bai;
        use noodles::core::Position;
        use noodles::csi::binning_index::{
            index::reference_sequence::{bin::Chunk, index::LinearIndex},
            Indexer,
        };
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;
        use std::fs::File;

        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(chr_len.try_into().unwrap())
            .build()
            .unwrap();
        let header = Header::builder().add_reference_sequence("chr1", rs).build();
        let tmp_bam = bam_path.with_extension("bam.tmp");
        let tmp_bai = bam_path.with_extension("bam.bai.tmp");
        let bgzf = noodles::bgzf::io::Writer::new(File::create(&tmp_bam).unwrap());
        let mut writer = noodles::bam::io::Writer::from(bgzf);
        writer.write_alignment_header(&header).unwrap();
        let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
        let read_len = 20usize;
        let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
        for (idx, &qname) in read_names.iter().enumerate() {
            let start_1 = 10 + idx * 30;
            let rec = RecordBuf::builder()
                .set_name(qname)
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start_1).unwrap())
                .set_mapping_quality(MappingQuality::new(60).unwrap())
                .set_cigar(cigar.clone())
                .set_sequence(b"ACGTACGTACGTACGTACGT".to_vec().into())
                .set_quality_scores(vec![35u8; read_len].into())
                .build();
            let s_vp = writer.get_ref().virtual_position();
            writer.write_alignment_record(&header, &rec).unwrap();
            let e_vp = writer.get_ref().virtual_position();
            let s = Position::try_from(start_1).unwrap();
            let e = Position::try_from(start_1 + read_len - 1).unwrap();
            indexer
                .add_record(Some((0, s, e, true)), Chunk::new(s_vp, e_vp))
                .unwrap();
        }
        writer.try_finish().unwrap();
        bai::io::Writer::new(File::create(&tmp_bai).unwrap())
            .write_index(&indexer.build(1))
            .unwrap();
        std::fs::rename(&tmp_bam, bam_path).unwrap();
        std::fs::rename(&tmp_bai, bam_path.with_extension("bam.bai")).unwrap();
    }

    struct ForceProdIdGuard;
    impl ForceProdIdGuard {
        fn enter() -> Self {
            FORCE_PRODUCTION_META_ID.with(|c| c.set(true));
            Self
        }
    }
    impl Drop for ForceProdIdGuard {
        fn drop(&mut self) {
            FORCE_PRODUCTION_META_ID.with(|c| c.set(false));
        }
    }

    #[test]
    fn test_bam_meta_tls_zero_stat_and_epoch_invalidation() {
        use std::sync::atomic::Ordering;

        let _prod_guard = ForceProdIdGuard::enter();

        let dir = tempfile::tempdir().unwrap();
        let bam_a = dir.path().join("sample_a.bam");
        let bam_b = dir.path().join("sample_b.bam");
        let config = Config::default();

        // Retry the TLS hit verification if a parallel unit test calls `clear_bam_caches()`
        // and bumps `BAM_CACHE_EPOCH` mid-window.
        let (meta_a_v1, meta_b_v1) = loop {
            write_indexed_bam(&bam_a, 1000, &["a_v1_r1"]);
            write_indexed_bam(&bam_b, 1200, &["b_v1_r1", "b_v1_r2"]);
            let epoch_start = BAM_CACHE_EPOCH.load(Ordering::Acquire);

            let recs_a1 = query_region(&bam_a, "chr1", Pos0(0), Pos0(500), &config).unwrap();
            let recs_b1 = query_region(&bam_b, "chr1", Pos0(0), Pos0(500), &config).unwrap();
            if BAM_CACHE_EPOCH.load(Ordering::Acquire) != epoch_start {
                continue;
            }
            assert_eq!(recs_a1.len(), 1);
            assert_eq!(recs_b1.len(), 2);

            let (tls_a, tls_b) = BAM_META_TLS.with(|tls| {
                let slots = tls.borrow();
                let s0 = slots[0].as_ref().map(|(_, m)| std::sync::Arc::clone(m));
                let s1 = slots[1].as_ref().map(|(_, m)| std::sync::Arc::clone(m));
                (s0, s1)
            });
            let (Some(m0), Some(m1)) = (tls_a, tls_b) else {
                panic!("both BAM_META_TLS slots must be populated after querying A and B");
            };
            // Most-recently-used slot[0] is bam_b, slot[1] is bam_a.
            assert_eq!(m0.path, bam_b);
            assert_eq!(m1.path, bam_a);

            // Evict A and B from META_CACHE and unlink `.bai` files: alternating queries on
            // the same thread must still succeed via 2-slot `BAM_META_TLS` + `BAM_READER_POOL`
            // without locking `META_CACHE` or calling `stat()` / opening `.bai`.
            if let Some(cache) = META_CACHE.get() {
                let mut map = cache.lock().unwrap_or_else(|e| e.into_inner());
                map.remove(&bam_a);
                map.remove(&bam_b);
            }
            std::fs::remove_file(bam_a.with_extension("bam.bai")).unwrap();
            std::fs::remove_file(bam_b.with_extension("bam.bai")).unwrap();

            let r_a2 = query_region(&bam_a, "chr1", Pos0(0), Pos0(500), &config);
            let r_b2 = query_region(&bam_b, "chr1", Pos0(0), Pos0(500), &config);
            let r_a3 = query_region(&bam_a, "chr1", Pos0(0), Pos0(500), &config);
            let r_b3 = query_region(&bam_b, "chr1", Pos0(0), Pos0(500), &config);
            if BAM_CACHE_EPOCH.load(Ordering::Acquire) != epoch_start {
                continue;
            }
            assert_eq!(r_a2.unwrap().len(), 1);
            assert_eq!(r_b2.unwrap().len(), 2);
            assert_eq!(r_a3.unwrap().len(), 1);
            assert_eq!(r_b3.unwrap().len(), 2);
            break (m1, m0);
        };
        assert_eq!(
            meta_a_v1
                .header
                .reference_sequences()
                .get(b"chr1".as_slice())
                .unwrap()
                .length()
                .get(),
            1000
        );
        assert_eq!(
            meta_b_v1
                .header
                .reference_sequences()
                .get(b"chr1".as_slice())
                .unwrap()
                .length()
                .get(),
            1200
        );

        // Also warm a worker thread's `BAM_META_TLS` before overwriting `bam_a` in place
        // so `clear_bam_caches()` proves cross-thread `BAM_CACHE_EPOCH` invalidation.
        write_indexed_bam(&bam_a, 1000, &["a_v1_r1"]);
        write_indexed_bam(&bam_b, 1200, &["b_v1_r1", "b_v1_r2"]);
        let (tx_cmd, rx_cmd) = std::sync::mpsc::channel::<()>();
        let (tx_ack, rx_ack) = std::sync::mpsc::channel::<(usize, Vec<Vec<u8>>)>();
        let bam_a_worker = bam_a.clone();
        let bam_b_worker = bam_b.clone();
        let worker = std::thread::spawn(move || {
            let _worker_prod_guard = ForceProdIdGuard::enter();
            let cfg = Config::default();
            // Warm worker's 2-slot TLS with A and B (v1)
            let w_a1 = query_region(&bam_a_worker, "chr1", Pos0(0), Pos0(500), &cfg).unwrap();
            let w_b1 = query_region(&bam_b_worker, "chr1", Pos0(0), Pos0(500), &cfg).unwrap();
            assert_eq!(w_a1.len(), 1);
            assert_eq!(w_b1.len(), 2);
            tx_ack
                .send((
                    read_reference_lengths(&bam_a_worker, &cfg).unwrap()["chr1"],
                    w_a1.into_iter().map(|r| r.name).collect(),
                ))
                .unwrap();

            // Wait until main thread overwrites `bam_a` and calls `clear_bam_caches()`
            rx_cmd.recv().unwrap();
            let lens = read_reference_lengths(&bam_a_worker, &cfg).unwrap();
            let w_a2 = query_region(&bam_a_worker, "chr1", Pos0(0), Pos0(500), &cfg).unwrap();
            tx_ack
                .send((lens["chr1"], w_a2.into_iter().map(|r| r.name).collect()))
                .unwrap();
        });

        let (worker_len_v1, worker_names_v1) = rx_ack.recv().unwrap();
        assert_eq!(worker_len_v1, 1000);
        assert_eq!(worker_names_v1, vec![b"a_v1_r1".to_vec()]);

        // Overwrite `bam_a` + `.bai` in place with updated header LN=2400 and 3 new reads,
        // then call `clear_bam_caches()` (bumping `BAM_CACHE_EPOCH`).
        let epoch_before_clear = BAM_CACHE_EPOCH.load(Ordering::Acquire);
        write_indexed_bam(&bam_a, 2400, &["a_v2_r1", "a_v2_r2", "a_v2_r3"]);
        clear_bam_caches();
        let epoch_after_clear = BAM_CACHE_EPOCH.load(Ordering::Acquire);
        assert!(
            epoch_after_clear > epoch_before_clear,
            "clear_bam_caches() must bump BAM_CACHE_EPOCH"
        );

        let updated_lens = read_reference_lengths(&bam_a, &config).unwrap();
        assert_eq!(updated_lens.get("chr1"), Some(&2400));
        let updated_recs = query_region(&bam_a, "chr1", Pos0(0), Pos0(500), &config).unwrap();
        let updated_names: Vec<Vec<u8>> = updated_recs.into_iter().map(|r| r.name).collect();
        assert_eq!(
            updated_names,
            vec![
                b"a_v2_r1".to_vec(),
                b"a_v2_r2".to_vec(),
                b"a_v2_r3".to_vec()
            ]
        );

        tx_cmd.send(()).unwrap();
        let (worker_len_v2, worker_names_v2) = rx_ack.recv().unwrap();
        worker.join().unwrap();
        assert_eq!(worker_len_v2, 2400);
        assert_eq!(
            worker_names_v2,
            vec![
                b"a_v2_r1".to_vec(),
                b"a_v2_r2".to_vec(),
                b"a_v2_r3".to_vec()
            ]
        );

        // Verify `<= 4096`-byte cap on `BamFeedInner` string stashes (`AR-14`).
        let mut stash = String::new();
        let mut small = String::with_capacity(256);
        small.push_str("chr1,100,+,50M,60,0;");
        stash_capped_string(&mut stash, small);
        assert!(stash.is_empty() && stash.capacity() >= 256 && stash.capacity() <= 4096);
        let oversized = String::with_capacity(8192);
        stash_capped_string(&mut stash, oversized);
        assert!(stash.capacity() <= 4096);
    }

    #[test]
    fn test_compute_cram_nm_counts_matching_n_and_iupac_as_mismatch() {
        for b in 0u8..=u8::MAX {
            let expected = matches!(
                b,
                b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'='
            );
            assert_eq!(
                is_strict_nm_match_base(b),
                expected,
                "byte 0x{b:02x} ({:?}) strict NM match classification mismatch",
                char::from(b)
            );
        }

        // 1. Matching ACGT / soft-masked acgt / '=' runs across multiple 8-byte chunks + remainder -> NM = 0
        let ref_mixed = b"ACGTACGTacgtacgtACGTACGTacgtacgtACGTa";
        let rec_exact = SamRecord {
            flags: 0,
            alignment_start: Some(Pos1(1)),
            cigar: vec![(ref_mixed.len() as u32, CigarKind::Match)],
            sequence: ref_mixed.to_vec(),
            ..SamRecord::default()
        };
        assert_eq!(compute_cram_nm(&rec_exact, ref_mixed), Some(0));

        let rec_upper_vs_softmask = SamRecord {
            sequence: ref_mixed.to_ascii_uppercase(),
            ..rec_exact.clone()
        };
        assert_eq!(compute_cram_nm(&rec_upper_vs_softmask, ref_mixed), Some(0));

        let rec_eq = SamRecord {
            sequence: vec![b'='; ref_mixed.len()],
            ..rec_exact.clone()
        };
        assert_eq!(
            compute_cram_nm(&rec_eq, b"NNNNRRRRACGTNNNNRRRRACGTNNNNRRRRACGTN"),
            Some(0)
        );

        // 2. SAMtags §1.5 (`AR-02`): matching 'N' == 'N', 'n' == 'n', and IUPAC 'R' == 'R' / 'Y' == 'Y'
        //    both inside 8-byte aligned chunks and in the <8-byte tail MUST each increment NM by +1!
        let ref_ambig = b"ACGTNACGRACGTACGTnACGTYACGTN";
        // Ambiguity positions (0-based): 4 ('N'), 8 ('R'), 17 ('n'), 22 ('Y'), 27 ('N') => 5 total
        let rec_ambig = SamRecord {
            flags: 0,
            alignment_start: Some(Pos1(1)),
            cigar: vec![(ref_ambig.len() as u32, CigarKind::Match)],
            sequence: ref_ambig.to_vec(),
            ..SamRecord::default()
        };
        assert_eq!(compute_cram_nm(&rec_ambig, ref_ambig), Some(5));

        // 3. 3' reference overhang: positions past `ref_seq.len()` compare against implicit 'N'
        //    ('=' matches -> +0, 'A' vs implicit 'N' -> +1, 'N' vs implicit 'N' -> +1).
        let short_ref = b"ACGTACGTACGT"; // 12 bp
        let rec_overhang = SamRecord {
            flags: 0,
            alignment_start: Some(Pos1(1)),
            cigar: vec![(15, CigarKind::Match)],
            sequence: b"ACGTACGTACGT=AN".to_vec(), // 12 exact + '=' (0) + 'A' (+1) + 'N' (+1) = 2
            ..SamRecord::default()
        };
        assert_eq!(compute_cram_nm(&rec_overhang, short_ref), Some(2));
    }

    #[test]
    fn test_bam_cache_epoch_post_clear_release_and_cold_miss_guard() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::header::record::value::{
            map::ReferenceSequence as RefSeqMap, Map as HeaderMap,
        };
        use noodles::sam::Header;
        use std::fs::File;
        use std::sync::atomic::Ordering;

        let _prod_guard = ForceProdIdGuard::enter();

        let dir = tempfile::tempdir().unwrap();
        let bam_race = dir.path().join("toctou_race.bam");
        let cfg = Config::default();

        // 1. TOCTOU mid-clear race: populate `META_CACHE` with v1 (`LN = 1000`), overwrite
        //    `bam_race` on disk with v2 (`LN = 3000`), then have a worker thread query `bam_race`
        //    mid-`clear_bam_caches()` BEFORE `META_CACHE` is cleared. Once `clear_bam_caches()`
        //    completes, its post-clear `Release` epoch bump must invalidate the worker's TLS.
        write_indexed_bam(&bam_race, 1000, &["v1_only"]);
        assert_eq!(
            read_reference_lengths(&bam_race, &cfg).unwrap().get("chr1"),
            Some(&1000)
        );
        write_indexed_bam(&bam_race, 3000, &["v2_r1", "v2_r2"]);

        let (tx_continue, rx_continue) = std::sync::mpsc::channel::<()>();
        let (tx_done, rx_done) = std::sync::mpsc::channel::<(usize, Vec<Vec<u8>>)>();
        let (tx_handle, rx_handle) = std::sync::mpsc::channel::<std::thread::JoinHandle<()>>();
        let bam_for_hook = bam_race.clone();
        TEST_MID_CLEAR_HOOK.with(|h| {
            let (tx_warmed, rx_warmed) = std::sync::mpsc::channel::<()>();
            *h.borrow_mut() = Some(Box::new(move || {
                let handle = std::thread::spawn(move || {
                    let _worker_prod_guard = ForceProdIdGuard::enter();
                    let wcfg = Config::default();
                    // Reads stale `v1` from `META_CACHE` mid-sweep and populates worker TLS.
                    let _ = read_reference_lengths(&bam_for_hook, &wcfg).unwrap();
                    let _ = query_region(&bam_for_hook, "chr1", Pos0(0), Pos0(500), &wcfg).unwrap();
                    tx_warmed.send(()).unwrap();

                    rx_continue.recv().unwrap();
                    let post_len = read_reference_lengths(&bam_for_hook, &wcfg).unwrap()["chr1"];
                    let post_names: Vec<Vec<u8>> =
                        query_region(&bam_for_hook, "chr1", Pos0(0), Pos0(500), &wcfg)
                            .unwrap()
                            .into_iter()
                            .map(|r| r.name)
                            .collect();
                    tx_done.send((post_len, post_names)).unwrap();
                });
                rx_warmed.recv().unwrap();
                tx_handle.send(handle).unwrap();
            }));
        });

        clear_bam_caches();
        tx_continue.send(()).unwrap();
        let (post_clear_len, post_clear_names) = rx_done.recv().unwrap();
        rx_handle.recv().unwrap().join().unwrap();
        assert_eq!(
            post_clear_len, 3000,
            "post-clear Release bump must invalidate worker TLS populated mid-clear"
        );
        assert_eq!(post_clear_names, vec![b"v2_r1".to_vec(), b"v2_r2".to_vec()]);

        // 2. Cold-miss epoch mismatch guard (BAM + CRAM): if `BAM_CACHE_EPOCH` advances while
        //    a cold miss has already opened pre-clear disk files, neither `META_CACHE` /
        //    `CRAM_META_CACHE` nor TLS may cache the pre-clear metadata.
        let bam_cold = dir.path().join("cold_race.bam");
        write_indexed_bam(&bam_cold, 1100, &["cold_v1"]);
        let bam_cold_hook = bam_cold.clone();
        TEST_MID_COLD_OPEN_HOOK.with(|h| {
            *h.borrow_mut() = Some(Box::new(move || {
                write_indexed_bam(&bam_cold_hook, 4200, &["cold_v2_a", "cold_v2_b"]);
                BAM_CACHE_EPOCH.fetch_add(1, Ordering::Release);
            }));
        });
        let first_cold_len = read_reference_lengths(&bam_cold, &cfg).unwrap()["chr1"];
        assert_eq!(first_cold_len, 1100);
        // Next query must NOT hit stale 1100 metadata in `META_CACHE` or `BAM_META_TLS`:
        let second_cold_len = read_reference_lengths(&bam_cold, &cfg).unwrap()["chr1"];
        assert_eq!(
            second_cold_len, 4200,
            "cold-miss BAM open spanning epoch bump must not cache stale SharedBamMeta"
        );

        // Same cold-miss epoch mismatch guard on CRAM:
        let ref_fa = dir.path().join("cold_ref.fa");
        let ref_fai = dir.path().join("cold_ref.fa.fai");
        std::fs::write(&ref_fa, b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
        std::fs::write(&ref_fai, b"chr1\t32\t6\t32\t33\n").unwrap();
        let cram_cold = dir.path().join("cold_race.cram");
        fn write_cram(cram_path: &Path, ref_fa: &Path, ref_fai: &Path, sq_len: usize) {
            let rs = HeaderMap::<RefSeqMap>::builder()
                .set_length(sq_len.try_into().unwrap())
                .build()
                .unwrap();
            let header = Header::builder().add_reference_sequence("chr1", rs).build();
            let tmp_cram = cram_path.with_extension("cram.tmp");
            let tmp_crai = cram_path.with_extension("cram.crai.tmp");
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(File::open(ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut cram_writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&tmp_cram)
                .unwrap();
            cram_writer.write_alignment_header(&header).unwrap();
            cram_writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(&tmp_cram).unwrap();
            noodles::cram::crai::fs::write(&tmp_crai, &crai).unwrap();
            std::fs::rename(&tmp_cram, cram_path).unwrap();
            std::fs::rename(&tmp_crai, cram_path.with_extension("cram.crai")).unwrap();
        }
        write_cram(&cram_cold, &ref_fa, &ref_fai, 32);
        let cram_cfg = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            ..Default::default()
        };
        let cram_cold_hook = cram_cold.clone();
        let ref_fa_hook = ref_fa.clone();
        let ref_fai_hook = ref_fai.clone();
        TEST_MID_COLD_OPEN_HOOK.with(|h| {
            *h.borrow_mut() = Some(Box::new(move || {
                write_cram(&cram_cold_hook, &ref_fa_hook, &ref_fai_hook, 64);
                BAM_CACHE_EPOCH.fetch_add(1, Ordering::Release);
            }));
        });
        assert_eq!(
            read_reference_lengths(&cram_cold, &cram_cfg).unwrap()["chr1"],
            32
        );
        assert_eq!(
            read_reference_lengths(&cram_cold, &cram_cfg).unwrap()["chr1"],
            64,
            "cold-miss CRAM open spanning epoch bump must not cache stale SharedCramMeta"
        );
    }

    #[test]
    fn fill_from_alignment_record_terminates_on_malformed_text_sam_cigar() {
        let header_src = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n";
        // Exact 41-byte `fuzz_sam_aux` `$`-CIGAR reproducer + prefix-valid malformed CIGAR line.
        let fuzz_repro: &[u8] =
            b"\t\t\xbd\x1a1\t1\t1\t$\t1\x1a\t\xbfO\x00\x00\xa8O\x00FHDI\x00\xa8+\t1\t1\t1\t1\x08\x01\t+";
        let mut sam_bytes = Vec::new();
        sam_bytes.extend_from_slice(header_src);
        sam_bytes.extend_from_slice(fuzz_repro);
        sam_bytes.push(b'\n');
        sam_bytes.extend_from_slice(b"r2\t0\tchr1\t10\t60\t2M$2M\t*\t0\t0\tACGT\tIIII\n");

        let mut reader = noodles::sam::io::Reader::new(std::io::Cursor::new(sam_bytes));
        let header = reader.read_header().unwrap();
        let cfg = Config::default();
        let mut preproc = RecordPreprocessor::new();

        let records: Vec<noodles::sam::Record> =
            reader.records().take(16).filter_map(Result::ok).collect();
        assert_eq!(records.len(), 2);

        let mut rec0 = SamRecord::new();
        rec0.fill_from_alignment_record(&records[0], &header)
            .unwrap();
        assert!(rec0.cigar.is_empty());
        assert!(rec0.alignment_end.is_none());
        let _ = preproc.preprocess_record(&rec0, &cfg);

        let mut rec0_legacy = SamRecord::new();
        rec0_legacy
            .fill_from_alignment_record_legacy(&records[0], &header)
            .unwrap();
        assert!(rec0_legacy.cigar.is_empty());
        assert!(rec0_legacy.alignment_end.is_none());

        let mut rec1 = SamRecord::new();
        rec1.fill_from_alignment_record(&records[1], &header)
            .unwrap();
        assert_eq!(rec1.cigar, vec![(2, CigarKind::Match)]);
        assert!(rec1.alignment_end.is_none());
    }

    #[test]
    fn bam_index_loader_rejects_oversized_bai_and_csi_counts() {
        use noodles::bam::bai;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::{
            bin::Chunk,
            index::{BinnedIndex, LinearIndex},
        };
        use noodles::csi::binning_index::Indexer;
        use noodles::sam::{
            alignment::{
                io::Write as _,
                record::{
                    cigar::{op::Kind, Op},
                    Flags, MappingQuality,
                },
                record_buf::{Cigar, QualityScores as BufQual, RecordBuf, Sequence as BufSeq},
            },
            header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap},
            Header,
        };
        use std::io::Write as _;

        fn bgzf_wrap(raw: &[u8]) -> Vec<u8> {
            let mut out = Vec::new();
            let mut w = noodles_bgzf::io::Writer::new(&mut out);
            w.write_all(raw).unwrap();
            w.finish().unwrap();
            out
        }

        assert_eq!(BAI_DEPTH, 5);
        assert_eq!(BAI_MAX_BIN_ID, 37_448);
        assert_eq!(BAI_METADATA_BIN_ID, 37_450);
        assert_eq!(MAX_CSI_TOTAL_BINS, 4_194_304);

        let dir = tempfile::tempdir().unwrap();
        let bam_path = dir.path().join("hostile_idx.bam");
        let bai_path = dir.path().join("hostile_idx.bam.bai");
        let csi_path = dir.path().join("hostile_idx.bam.csi");

        let rs = HeaderMap::<RefSeqMap>::builder()
            .set_length(1000.try_into().unwrap())
            .build()
            .unwrap();
        let header = Header::builder().add_reference_sequence("chr1", rs).build();
        let rec = RecordBuf::builder()
            .set_name("read1")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(50usize).unwrap())
            .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
            .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect::<Cigar>())
            .set_sequence(BufSeq::from(b"ACGTACGTAC".to_vec()))
            .set_quality_scores(BufQual::from(vec![35u8; 10]))
            .build();

        let mut bai_indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
        let mut csi_indexer = Indexer::<BinnedIndex>::new(14, 5).unwrap();
        let mut csi_hdr_indexer = Indexer::<BinnedIndex>::new(14, 5).unwrap();
        let (s_vp, e_vp) = {
            let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            let mut w = noodles::bam::io::Writer::from(bgzf);
            w.write_alignment_header(&header).unwrap();
            let s_vp = w.get_ref().virtual_position();
            w.write_alignment_record(&header, &rec).unwrap();
            let e_vp = w.get_ref().virtual_position();
            let s = Position::try_from(50usize).unwrap();
            let e = Position::try_from(59usize).unwrap();
            let chunk = Chunk::new(s_vp, e_vp);
            bai_indexer
                .add_record(Some((0, s, e, true)), chunk)
                .unwrap();
            csi_indexer
                .add_record(Some((0, s, e, true)), chunk)
                .unwrap();
            csi_hdr_indexer
                .add_record(Some((0, s, e, true)), chunk)
                .unwrap();
            w.try_finish().unwrap();
            (s_vp, e_vp)
        };

        let cfg = Config::default();
        let clone_bam_case = |stem: &str| -> PathBuf {
            let dst = dir.path().join(format!("{stem}.bam"));
            std::fs::copy(&bam_path, &dst).unwrap();
            dst
        };

        // 1. Valid `.bai` loads and queries cleanly (including max valid normal bin `37_448`).
        let ok_bai_bam = clone_bam_case("ok_bai");
        bai::fs::write(ok_bai_bam.with_extension("bam.bai"), &bai_indexer.build(1)).unwrap();
        let valid_bai_bytes = std::fs::read(ok_bai_bam.with_extension("bam.bai")).unwrap();
        let bai_reads: Vec<SamRecord> =
            query_region_iter(&ok_bai_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert_eq!(bai_reads.len(), 1);
        assert_eq!(bai_reads[0].name, b"read1");

        let max_bin_bam = clone_bam_case("max_bin_37448");
        let mut bai_max_valid_bin = b"BAI\x01".to_vec();
        bai_max_valid_bin.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bai_max_valid_bin.extend_from_slice(&1u32.to_le_bytes()); // n_bin = 1
        bai_max_valid_bin.extend_from_slice(&37_448u32.to_le_bytes()); // bin_id = 37_448
        bai_max_valid_bin.extend_from_slice(&0i32.to_le_bytes()); // n_chunk = 0
        bai_max_valid_bin.extend_from_slice(&0u32.to_le_bytes()); // n_intv = 0
        std::fs::write(max_bin_bam.with_extension("bam.bai"), &bai_max_valid_bin).unwrap();
        assert!(read_bam_index_with_fallbacks(&max_bin_bam).is_ok());
        let empty_max_bin_reads: Vec<SamRecord> =
            query_region_iter(&max_bin_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
        assert!(empty_max_bin_reads.is_empty());

        // 2. Valid `.csi` (headerless, headered, `min_shift = 0`, and `min_shift = 4`) loads identically.
        let ok_csi_bam = clone_bam_case("ok_csi");
        noodles::csi::fs::write(ok_csi_bam.with_extension("bam.csi"), &csi_indexer.build(1))
            .unwrap();
        let csi_reads: Vec<SamRecord> =
            query_region_iter(&ok_csi_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert_eq!(csi_reads, bai_reads);

        for valid_min_shift in [0u8, 4u8] {
            let shift_bam = clone_bam_case(&format!("ok_csi_shift_{valid_min_shift}"));
            let mut shift_indexer = Indexer::<BinnedIndex>::new(valid_min_shift, 5).unwrap();
            let s = Position::try_from(50usize).unwrap();
            let e = Position::try_from(59usize).unwrap();
            shift_indexer
                .add_record(Some((0, s, e, true)), Chunk::new(s_vp, e_vp))
                .unwrap();
            noodles::csi::fs::write(shift_bam.with_extension("bam.csi"), &shift_indexer.build(1))
                .unwrap();
            let shift_reads: Vec<SamRecord> =
                query_region_iter(&shift_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            assert_eq!(
                shift_reads, bai_reads,
                "valid CSI min_shift={valid_min_shift} depth=5 must load and query cleanly"
            );
        }

        {
            let d9_bam = clone_bam_case("ok_csi_depth_9");
            let mut d9_indexer = Indexer::<BinnedIndex>::new(14, 9).unwrap();
            let s = Position::try_from(50usize).unwrap();
            let e = Position::try_from(59usize).unwrap();
            d9_indexer
                .add_record(Some((0, s, e, true)), Chunk::new(s_vp, e_vp))
                .unwrap();
            noodles::csi::fs::write(d9_bam.with_extension("bam.csi"), &d9_indexer.build(1))
                .unwrap();
            let d9_reads: Vec<SamRecord> =
                query_region_iter(&d9_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                    .unwrap()
                    .map(|r| r.unwrap())
                    .collect();
            assert_eq!(
                d9_reads, bai_reads,
                "valid CSI min_shift=14 depth=9 must load and query cleanly"
            );
        }

        let hdr_csi_bam = clone_bam_case("ok_csi_hdr");
        let mut names = noodles::csi::binning_index::index::header::ReferenceSequenceNames::new();
        names.insert(b"chr1".as_slice().into());
        csi_hdr_indexer = csi_hdr_indexer.set_header(
            noodles::csi::binning_index::index::header::Builder::sam()
                .set_reference_sequence_names(names)
                .build(),
        );
        noodles::csi::fs::write(
            hdr_csi_bam.with_extension("bam.csi"),
            &csi_hdr_indexer.build(1),
        )
        .unwrap();
        let headered_csi_reads: Vec<SamRecord> =
            query_region_iter(&hdr_csi_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
        assert_eq!(headered_csi_reads, bai_reads);

        // Missing index still preserves `ErrorKind::NotFound` (`bam_path` has no index on disk).
        assert_eq!(
            read_bam_index_with_fallbacks(&bam_path).unwrap_err().kind(),
            std::io::ErrorKind::NotFound
        );

        let assert_io_invalid_data = |err: crate::types::VardictError, label: &str| match err {
            crate::types::VardictError::Io(io_err) => assert_eq!(
                io_err.kind(),
                std::io::ErrorKind::InvalidData,
                "{label}: expected InvalidData, got {io_err:?}"
            ),
            other => panic!("{label}: expected VardictError::Io(InvalidData), got {other:?}"),
        };

        // 3. Corrupt `.bai` cases (with complete trailing fields so truncation never masks validation):
        //    `bin_id = 37_449`, `bin_id = 50_000`, `n_bin = 37_451`, decoded fuzz seed `gprBLSgoLSgoCg==`,
        //    fuzzer `n_ref = 0x6573776f`, oversized `n_bin`, oversized/negative `n_chunk`,
        //    truncated metadata bin, oversized `n_intv`.
        let build_single_bin_bai = |bin_id: u32| -> Vec<u8> {
            let mut v = b"BAI\x01".to_vec();
            v.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
            v.extend_from_slice(&1u32.to_le_bytes()); // n_bin = 1
            v.extend_from_slice(&bin_id.to_le_bytes());
            v.extend_from_slice(&0i32.to_le_bytes()); // n_chunk = 0
            v.extend_from_slice(&0u32.to_le_bytes()); // n_intv = 0
            v
        };
        let bai_bin_37449 = build_single_bin_bai(37_449);
        let bai_bin_50000 = build_single_bin_bai(50_000);

        let mut bai_nbin_37451 = Vec::with_capacity(4 + 8 + 37_451 * 8 + 36);
        bai_nbin_37451.extend_from_slice(b"BAI\x01");
        bai_nbin_37451.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bai_nbin_37451.extend_from_slice(&37_451u32.to_le_bytes()); // n_bin = 37_451 (> 37_450)
        for bid in 0..37_450u32 {
            bai_nbin_37451.extend_from_slice(&bid.to_le_bytes());
            bai_nbin_37451.extend_from_slice(&0i32.to_le_bytes()); // n_chunk = 0
        }
        bai_nbin_37451.extend_from_slice(&37_450u32.to_le_bytes()); // metadata bin
        bai_nbin_37451.extend_from_slice(&2u32.to_le_bytes()); // metadata n_chunk = 2
        bai_nbin_37451.extend_from_slice(&[0u8; 32]); // ref_beg, ref_end, n_mapped, n_unmapped
        bai_nbin_37451.extend_from_slice(&0u32.to_le_bytes()); // n_intv = 0

        // Decoded base64 fuzz seed `gprBLSgoLSgoCg==` (`[130, 154, 193, 45, 40, 40, 45, 40, 40, 10]`)
        // both raw and applied via `fuzz_cram_and_fai`'s XOR mutation past `BAI\x01`.
        const FUZZ_SEED_GPR: &[u8] = b"\x82\x9a\xc1-((-((\n";
        let mut bai_fuzz_gpr_mutated = valid_bai_bytes.clone();
        let gpr_off = 4 + (usize::from(FUZZ_SEED_GPR[9]) % (bai_fuzz_gpr_mutated.len() - 4));
        bai_fuzz_gpr_mutated[gpr_off] ^= 0xff;

        let mut bai_oversized_nbin = b"BAI\x01".to_vec();
        bai_oversized_nbin.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bai_oversized_nbin.extend_from_slice(&0x4000_0000u32.to_le_bytes()); // n_bin
        bai_oversized_nbin.extend_from_slice(&0u32.to_le_bytes()); // n_intv

        let mut bai_oversized_nchunk = b"BAI\x01".to_vec();
        bai_oversized_nchunk.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bai_oversized_nchunk.extend_from_slice(&1u32.to_le_bytes()); // n_bin = 1
        bai_oversized_nchunk.extend_from_slice(&0u32.to_le_bytes()); // bin_id = 0
        bai_oversized_nchunk.extend_from_slice(&0x4000_0000i32.to_le_bytes()); // n_chunk
        bai_oversized_nchunk.extend_from_slice(&0u32.to_le_bytes()); // n_intv

        let mut bai_neg_nchunk = b"BAI\x01".to_vec();
        bai_neg_nchunk.extend_from_slice(&1u32.to_le_bytes());
        bai_neg_nchunk.extend_from_slice(&1u32.to_le_bytes());
        bai_neg_nchunk.extend_from_slice(&0u32.to_le_bytes());
        bai_neg_nchunk.extend_from_slice(&(-1i32).to_le_bytes());
        bai_neg_nchunk.extend_from_slice(&0u32.to_le_bytes());

        let mut bai_trunc_meta = b"BAI\x01".to_vec();
        bai_trunc_meta.extend_from_slice(&1u32.to_le_bytes());
        bai_trunc_meta.extend_from_slice(&1u32.to_le_bytes());
        bai_trunc_meta.extend_from_slice(&37_450u32.to_le_bytes());
        bai_trunc_meta.extend_from_slice(&[0u8; 12]);

        let mut bai_oversized_nintv = b"BAI\x01".to_vec();
        bai_oversized_nintv.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bai_oversized_nintv.extend_from_slice(&0u32.to_le_bytes()); // n_bin = 0
        bai_oversized_nintv.extend_from_slice(&0x4000_0000u32.to_le_bytes()); // n_intv

        let mut bai_dup_bin_4681 = b"BAI\x01".to_vec();
        bai_dup_bin_4681.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bai_dup_bin_4681.extend_from_slice(&2u32.to_le_bytes()); // n_bin = 2
        for _ in 0..2 {
            bai_dup_bin_4681.extend_from_slice(&4681u32.to_le_bytes());
            bai_dup_bin_4681.extend_from_slice(&0i32.to_le_bytes()); // n_chunk = 0
        }
        bai_dup_bin_4681.extend_from_slice(&0u32.to_le_bytes()); // n_intv = 0

        for (label, bad_bai) in [
            ("bai_bin_37449", bai_bin_37449),
            ("bai_bin_50000", bai_bin_50000),
            ("bai_nbin_37451", bai_nbin_37451),
            ("bai_dup_bin_4681", bai_dup_bin_4681),
            ("bai_fuzz_seed_gpr_raw", FUZZ_SEED_GPR.to_vec()),
            ("bai_fuzz_seed_gpr_mutated", bai_fuzz_gpr_mutated),
            ("bai_owser_nref", b"BAI\x01owser-".to_vec()),
            ("bai_oversized_nbin", bai_oversized_nbin),
            ("bai_oversized_nchunk", bai_oversized_nchunk),
            ("bai_neg_nchunk", bai_neg_nchunk),
            ("bai_trunc_meta", bai_trunc_meta),
            ("bai_oversized_nintv", bai_oversized_nintv),
        ] {
            std::fs::write(&bai_path, &bad_bai).unwrap();
            let idx_err = read_bam_index_with_fallbacks(&bam_path)
                .expect_err(&format!("{label}: expected BAI index error"));
            assert_eq!(
                idx_err.kind(),
                std::io::ErrorKind::InvalidData,
                "{label}: expected InvalidData, got {idx_err:?}"
            );
            if label == "bai_dup_bin_4681" {
                assert!(
                    idx_err.to_string().contains("duplicate BAI bin_id"),
                    "{label}: expected prevalidator duplicate BAI bin_id error, got {idx_err}"
                );
            }
            let q_err = query_region_iter(&bam_path, "chr1", Pos0(0), Pos0(200), &cfg)
                .expect_err(&format!("{label}: expected query_region_iter error"));
            assert_io_invalid_data(q_err, label);
        }
        std::fs::remove_file(&bai_path).unwrap();

        // 4. Corrupt BGZF `.csi` cases (with complete trailing fields):
        //    `bin_id = 37_449` (`depth = 5`), `n_bin = 37_451`, duplicate `bin_id` at `depth = 9`,
        //    total bins > `MAX_CSI_TOTAL_BINS` (`4_194_304`) with complete payload bytes,
        //    `(min_shift = 64, depth = 0)`, `(min_shift = 50, depth = 5)`, oversized `n_ref`,
        //    `n_bin`, `n_chunk`, `l_aux`, and misaligned `l_aux != 28 + l_nm`.
        let csi_hdr_custom = |min_shift: i32, depth: i32, l_aux: i32| -> Vec<u8> {
            let mut v = b"CSI\x01".to_vec();
            v.extend_from_slice(&min_shift.to_le_bytes());
            v.extend_from_slice(&depth.to_le_bytes());
            v.extend_from_slice(&l_aux.to_le_bytes());
            v
        };
        let csi_hdr = |l_aux: i32| -> Vec<u8> { csi_hdr_custom(14, 5, l_aux) };

        let mut csi_bin_37449 = csi_hdr(0);
        csi_bin_37449.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_bin_37449.extend_from_slice(&1i32.to_le_bytes()); // n_bin = 1
        csi_bin_37449.extend_from_slice(&37_449u32.to_le_bytes()); // bin_id = 37_449
        csi_bin_37449.extend_from_slice(&0u64.to_le_bytes()); // loffset = 0
        csi_bin_37449.extend_from_slice(&0i32.to_le_bytes()); // n_chunk = 0

        let mut csi_nbin_37451 = Vec::with_capacity(24 + 37_451 * 16 + 36);
        csi_nbin_37451.extend_from_slice(&csi_hdr(0));
        csi_nbin_37451.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_nbin_37451.extend_from_slice(&37_451i32.to_le_bytes()); // n_bin = 37_451
        for bid in 0..37_450u32 {
            csi_nbin_37451.extend_from_slice(&bid.to_le_bytes());
            csi_nbin_37451.extend_from_slice(&0u64.to_le_bytes()); // loffset
            csi_nbin_37451.extend_from_slice(&0i32.to_le_bytes()); // n_chunk = 0
        }
        csi_nbin_37451.extend_from_slice(&37_450u32.to_le_bytes()); // metadata bin
        csi_nbin_37451.extend_from_slice(&0u64.to_le_bytes()); // loffset
        csi_nbin_37451.extend_from_slice(&2u32.to_le_bytes()); // metadata n_chunk = 2
        csi_nbin_37451.extend_from_slice(&[0u8; 32]);

        let mut csi_d9_dup_bin = csi_hdr_custom(14, 9, 0);
        csi_d9_dup_bin.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_d9_dup_bin.extend_from_slice(&2i32.to_le_bytes()); // n_bin = 2
        for _ in 0..2 {
            csi_d9_dup_bin.extend_from_slice(&0u32.to_le_bytes()); // repeated bin_id = 0
            csi_d9_dup_bin.extend_from_slice(&0u64.to_le_bytes()); // loffset = 0
            csi_d9_dup_bin.extend_from_slice(&0i32.to_le_bytes()); // n_chunk = 0
        }
        let dup_pre_err = validate_csi_payload(&csi_d9_dup_bin, Path::new("csi_d9_dup_bin.csi"))
            .expect_err(
                "validate_csi_payload must reject duplicate CSI bin_id before noodles allocates",
            );
        assert_eq!(dup_pre_err.kind(), std::io::ErrorKind::InvalidData);
        assert!(dup_pre_err.to_string().contains("duplicate CSI bin_id"));

        const OVER_TOTAL_CSI_BINS: usize = MAX_CSI_TOTAL_BINS + 1;
        let mut csi_d9_over_total_bins = vec![0u8; 24 + OVER_TOTAL_CSI_BINS * 16];
        csi_d9_over_total_bins[..16].copy_from_slice(&csi_hdr_custom(14, 9, 0));
        csi_d9_over_total_bins[16..20].copy_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_d9_over_total_bins[20..24].copy_from_slice(&(OVER_TOTAL_CSI_BINS as i32).to_le_bytes());
        for (bid, slot) in csi_d9_over_total_bins[24..]
            .chunks_exact_mut(16)
            .enumerate()
        {
            slot[..4].copy_from_slice(&(bid as u32).to_le_bytes());
        }

        let mut csi_d7_multiref_over_total_bins = vec![0u8; 20 + 8 + OVER_TOTAL_CSI_BINS * 16];
        csi_d7_multiref_over_total_bins[..16].copy_from_slice(&csi_hdr_custom(14, 7, 0));
        csi_d7_multiref_over_total_bins[16..20].copy_from_slice(&2i32.to_le_bytes()); // n_ref = 2
        let ref0_bins = 2_097_153usize;
        let ref1_bins = OVER_TOTAL_CSI_BINS - ref0_bins; // 2_097_152 (both <= metadata_id(7) == 2_396_746)
        csi_d7_multiref_over_total_bins[20..24].copy_from_slice(&(ref0_bins as i32).to_le_bytes());
        let ref0_end = 24 + ref0_bins * 16;
        for (bid, slot) in csi_d7_multiref_over_total_bins[24..ref0_end]
            .chunks_exact_mut(16)
            .enumerate()
        {
            slot[..4].copy_from_slice(&(bid as u32).to_le_bytes());
        }
        csi_d7_multiref_over_total_bins[ref0_end..ref0_end + 4]
            .copy_from_slice(&(ref1_bins as i32).to_le_bytes());
        for (bid, slot) in csi_d7_multiref_over_total_bins[ref0_end + 4..]
            .chunks_exact_mut(16)
            .enumerate()
        {
            slot[..4].copy_from_slice(&(bid as u32).to_le_bytes());
        }

        let mut csi_shift_64_d0 = csi_hdr_custom(64, 0, 0);
        csi_shift_64_d0.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_shift_64_d0.extend_from_slice(&0i32.to_le_bytes()); // n_bin = 0

        let mut csi_shift_50_d5 = csi_hdr_custom(50, 5, 0);
        csi_shift_50_d5.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_shift_50_d5.extend_from_slice(&0i32.to_le_bytes()); // n_bin = 0

        let mut csi_bad_nref = csi_hdr(0);
        csi_bad_nref.extend_from_slice(&0x4000_0000i32.to_le_bytes());

        let mut csi_bad_nbin = csi_hdr(0);
        csi_bad_nbin.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_bad_nbin.extend_from_slice(&0x4000_0000i32.to_le_bytes()); // n_bin

        let mut csi_bad_nchunk = csi_hdr(0);
        csi_bad_nchunk.extend_from_slice(&1i32.to_le_bytes()); // n_ref = 1
        csi_bad_nchunk.extend_from_slice(&1i32.to_le_bytes()); // n_bin = 1
        csi_bad_nchunk.extend_from_slice(&0u32.to_le_bytes()); // bin_id = 0
        csi_bad_nchunk.extend_from_slice(&0u64.to_le_bytes()); // loffset
        csi_bad_nchunk.extend_from_slice(&0x4000_0000i32.to_le_bytes()); // n_chunk

        let csi_bad_laux = csi_hdr(0x4000_0000);

        let mut csi_desync_laux = csi_hdr(40);
        csi_desync_laux.extend_from_slice(&2i32.to_le_bytes()); // format = VCF
        csi_desync_laux.extend_from_slice(&1i32.to_le_bytes()); // col_seq
        csi_desync_laux.extend_from_slice(&2i32.to_le_bytes()); // col_beg
        csi_desync_laux.extend_from_slice(&0i32.to_le_bytes()); // col_end
        csi_desync_laux.extend_from_slice(&35i32.to_le_bytes()); // meta
        csi_desync_laux.extend_from_slice(&0i32.to_le_bytes()); // skip
        csi_desync_laux.extend_from_slice(&0i32.to_le_bytes()); // l_nm = 0 (28 != 40!)
        csi_desync_laux.extend_from_slice(&0x7fff_ffffi32.to_le_bytes()); // hidden n_ref bomb
        csi_desync_laux.extend_from_slice(&[0u8; 8]);
        csi_desync_laux.extend_from_slice(&0i32.to_le_bytes());

        for (label, bad_csi) in [
            ("csi_bin_37449", csi_bin_37449),
            ("csi_nbin_37451", csi_nbin_37451),
            ("csi_d9_dup_bin", csi_d9_dup_bin),
            ("csi_d9_over_total_bins", csi_d9_over_total_bins),
            (
                "csi_d7_multiref_over_total_bins",
                csi_d7_multiref_over_total_bins,
            ),
            ("csi_shift_64_d0", csi_shift_64_d0),
            ("csi_shift_50_d5", csi_shift_50_d5),
            ("csi_bad_nref", csi_bad_nref),
            ("csi_bad_nbin", csi_bad_nbin),
            ("csi_bad_nchunk", csi_bad_nchunk),
            ("csi_bad_laux", csi_bad_laux),
            ("csi_desync_laux", csi_desync_laux),
        ] {
            std::fs::write(&csi_path, bgzf_wrap(&bad_csi)).unwrap();
            let idx_err = read_bam_index_with_fallbacks(&bam_path)
                .expect_err(&format!("{label}: expected CSI index error"));
            assert_eq!(
                idx_err.kind(),
                std::io::ErrorKind::InvalidData,
                "{label}: expected InvalidData, got {idx_err:?}"
            );
            if label == "csi_d9_dup_bin" {
                assert!(
                    idx_err.to_string().contains("duplicate CSI bin_id"),
                    "{label}: expected prevalidator duplicate CSI bin_id error, got {idx_err}"
                );
            }
            if matches!(
                label,
                "csi_d9_over_total_bins" | "csi_d7_multiref_over_total_bins"
            ) {
                assert!(
                    idx_err.to_string().contains("total bin count"),
                    "{label}: expected total bin count error, got {idx_err}"
                );
            }
            let q_err = query_region_iter(&bam_path, "chr1", Pos0(0), Pos0(200), &cfg)
                .expect_err(&format!("{label}: expected query_region_iter error"));
            assert_io_invalid_data(q_err, label);
        }
        std::fs::remove_file(&csi_path).unwrap();

        // 5. Hostile BAM header (`l_name`, `l_text`) & record (`block_size`) bounds:
        //    4 GiB (`u32::MAX`), oversized `l_name = 65_537` (with complete trailing fields),
        //    oversized `block_size > 64 MiB`, and truncated `block_size` payloads/prefixes.
        let mut bam_lname_65537 = b"BAM\x01".to_vec();
        bam_lname_65537.extend_from_slice(&0u32.to_le_bytes()); // l_text = 0
        bam_lname_65537.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bam_lname_65537.extend_from_slice(&65_537u32.to_le_bytes()); // l_name = 65_537 (> 65_536)
        bam_lname_65537.extend(std::iter::repeat_n(b'c', 65_536));
        bam_lname_65537.push(0); // NUL terminator
        bam_lname_65537.extend_from_slice(&1000u32.to_le_bytes()); // l_ref = 1000

        let mut bam_lname_4gib = b"BAM\x01".to_vec();
        bam_lname_4gib.extend_from_slice(&0u32.to_le_bytes()); // l_text = 0
        bam_lname_4gib.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
        bam_lname_4gib.extend_from_slice(&u32::MAX.to_le_bytes()); // l_name = 4 GiB

        let mut bam_ltext_4gib = b"BAM\x01".to_vec();
        bam_ltext_4gib.extend_from_slice(&u32::MAX.to_le_bytes()); // l_text = 4 GiB

        for (label, raw_hdr) in [
            ("bam_lname_65537", bam_lname_65537),
            ("bam_lname_4gib", bam_lname_4gib),
            ("bam_ltext_4gib", bam_ltext_4gib),
        ] {
            let hdr_bam = dir.path().join(format!("{label}.bam"));
            std::fs::write(hdr_bam.with_extension("bam.bai"), &valid_bai_bytes).unwrap();
            std::fs::write(&hdr_bam, bgzf_wrap(&raw_hdr)).unwrap();
            let q_err = query_region_iter(&hdr_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                .expect_err(&format!("{label}: expected BAM header rejection"));
            assert_io_invalid_data(q_err, label);
            #[expect(deprecated, reason = "verify open_reader BAM header hardening")]
            {
                let open_err = open_reader(&hdr_bam, &cfg)
                    .err()
                    .unwrap_or_else(|| panic!("{label}: expected open_reader rejection"));
                assert_io_invalid_data(open_err, label);
            }
        }

        // Valid tiny BAM header (`chr1`, LN=1000) + hostile record `block_size` cases.
        let valid_raw_bam_hdr = {
            let mut h = b"BAM\x01".to_vec();
            let sam_txt = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
            h.extend_from_slice(&(sam_txt.len() as u32).to_le_bytes());
            h.extend_from_slice(sam_txt);
            h.extend_from_slice(&1u32.to_le_bytes()); // n_ref = 1
            h.extend_from_slice(&5u32.to_le_bytes()); // l_name = 5
            h.extend_from_slice(b"chr1\0");
            h.extend_from_slice(&1000u32.to_le_bytes()); // l_ref = 1000
            h
        };
        let write_single_block_bam_and_bai = |stem: &str, rec_payload: &[u8]| -> PathBuf {
            let case_bam = dir.path().join(format!("{stem}.bam"));
            let mut raw = valid_raw_bam_hdr.clone();
            let rec_upos = raw.len() as u16;
            raw.extend_from_slice(rec_payload);
            let end_upos = raw.len().min(65_535) as u16;
            std::fs::write(&case_bam, bgzf_wrap(&raw)).unwrap();
            let c_start = noodles::bgzf::VirtualPosition::try_from((0u64, rec_upos)).unwrap();
            let c_end = if raw.len() <= 65_535 {
                noodles::bgzf::VirtualPosition::try_from((0u64, end_upos)).unwrap()
            } else {
                noodles::bgzf::VirtualPosition::try_from((1u64 << 20, 0u16)).unwrap()
            };
            let mut idxr = Indexer::<LinearIndex>::new(14, 5).unwrap();
            idxr.add_record(
                Some((
                    0,
                    Position::try_from(50usize).unwrap(),
                    Position::try_from(59usize).unwrap(),
                    true,
                )),
                Chunk::new(c_start, c_end),
            )
            .unwrap();
            bai::fs::write(case_bam.with_extension("bam.bai"), &idxr.build(1)).unwrap();
            case_bam
        };

        // Positive check for valid `65_536 < block_size <= 64 MiB` slow-path in `read_bounded_bam_record`.
        {
            let long_bam = dir.path().join("long_70k.bam");
            let long_bai = dir.path().join("long_70k.bam.bai");
            let long_len = 50_000usize; // 4-bit seq (25 KB) + qual (50 KB) => block_size > 75 KB > 65_536
            let long_rec = RecordBuf::builder()
                .set_name("long_read_70k")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(50usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
                .set_cigar(
                    [
                        Op::new(Kind::SoftClip, long_len - 10),
                        Op::new(Kind::Match, 10),
                    ]
                    .into_iter()
                    .collect::<Cigar>(),
                )
                .set_sequence(BufSeq::from(vec![b'A'; long_len]))
                .set_quality_scores(BufQual::from(vec![35u8; long_len]))
                .build();
            let bgzf = noodles::bgzf::io::Writer::new(std::fs::File::create(&long_bam).unwrap());
            let mut w = noodles::bam::io::Writer::from(bgzf);
            w.write_alignment_header(&header).unwrap();
            let ls_vp = w.get_ref().virtual_position();
            w.write_alignment_record(&header, &long_rec).unwrap();
            let le_vp = w.get_ref().virtual_position();
            w.try_finish().unwrap();
            let mut lidx = Indexer::<LinearIndex>::new(14, 5).unwrap();
            lidx.add_record(
                Some((
                    0,
                    Position::try_from(50usize).unwrap(),
                    Position::try_from(59usize).unwrap(),
                    true,
                )),
                Chunk::new(ls_vp, le_vp),
            )
            .unwrap();
            bai::fs::write(&long_bai, &lidx.build(1)).unwrap();
            let long_reads: Vec<SamRecord> =
                query_region_iter(&long_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap();
            assert_eq!(long_reads.len(), 1);
            assert_eq!(long_reads[0].sequence.len(), long_len);
        }

        for (label, bad_bs) in [
            ("bam_block_size_4gib", u32::MAX),
            ("bam_block_size_over_64mib", (64 * 1024 * 1024 + 1) as u32),
        ] {
            let mut payload = bad_bs.to_le_bytes().to_vec();
            payload.extend_from_slice(&[0u8; 64]);
            let case_bam = write_single_block_bam_and_bai(label, &payload);
            let feed_err = query_region_iter(&case_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .expect_err(&format!("{label}: expected BamFeed block_size rejection"));
            assert_io_invalid_data(feed_err, label);
            #[expect(deprecated, reason = "verify IndexedSamReader::query block_size bound")]
            {
                let mut cached = open_reader(&case_bam, &cfg).unwrap();
                let reg = build_noodles_region("chr1", 1, 200).unwrap();
                let q_err = cached
                    .reader
                    .query(&cached.header, &reg)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .expect_err(&format!("{label}: expected IndexedSamReader rejection"));
                assert_eq!(
                    q_err.kind(),
                    std::io::ErrorKind::InvalidData,
                    "{label}: expected InvalidData, got {q_err:?}"
                );
            }
        }

        // Truncated 64 MiB `block_size` record (16 trailing bytes) + 2-byte truncated prefix -> UnexpectedEof.
        for (label, trunc_payload) in [
            ("bam_trunc_prefix_2b", vec![0x24u8, 0x00]),
            ("bam_trunc_64mib_record", {
                let mut p = (64u32 * 1024 * 1024).to_le_bytes().to_vec();
                p.extend_from_slice(&[0u8; 16]);
                p
            }),
        ] {
            let case_bam = write_single_block_bam_and_bai(label, &trunc_payload);
            let feed_err = query_region_iter(&case_bam, "chr1", Pos0(0), Pos0(200), &cfg)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .expect_err(&format!("{label}: expected UnexpectedEof"));
            match feed_err {
                crate::types::VardictError::Io(io_err) => assert_eq!(
                    io_err.kind(),
                    std::io::ErrorKind::UnexpectedEof,
                    "{label}: expected UnexpectedEof, got {io_err:?}"
                ),
                other => panic!("{label}: expected Io(UnexpectedEof), got {other:?}"),
            }
        }

        // 6. Hostile CRAM `FileHeader` container (> 128 MiB) & `SliceHeader` coordinate bounds.
        let cram_ref_fa = dir.path().join("cram_bounds.fa");
        let cram_ref_fai = dir.path().join("cram_bounds.fa.fai");
        let cram_ref_seq = b"ACGT".repeat(250); // 1000 bp
        std::fs::write(
            &cram_ref_fa,
            format!(">chr1\n{}\n", std::str::from_utf8(&cram_ref_seq).unwrap()),
        )
        .unwrap();
        std::fs::write(&cram_ref_fai, "chr1\t1000\t6\t1000\t1001\n").unwrap();
        let base_cram = dir.path().join("base_bounds.cram");
        {
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(std::fs::File::open(&cram_ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&cram_ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut cw = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(&base_cram)
                .unwrap();
            cw.write_alignment_header(&header).unwrap();
            cw.write_alignment_record(&header, &rec).unwrap();
            cw.try_finish(&header).unwrap();
        }
        let base_cram_bytes = std::fs::read(&base_cram).unwrap();
        let base_crai = noodles::cram::fs::index(&base_cram).unwrap();
        let cram_cfg = Config {
            reference: cram_ref_fa.to_string_lossy().into_owned(),
            ..Default::default()
        };

        // Spoof initial CRAM header container `length = 128 MiB + 1` at offset 26..30 (recomputing header CRC32).
        let huge_hdr_cram = dir.path().join("huge_hdr.cram");
        {
            let mut rdr = noodles::cram::io::Reader::new(std::io::Cursor::new(&base_cram_bytes));
            rdr.read_file_definition().unwrap();
            let mut hc_rdr = rdr.header_reader();
            let _ = hc_rdr.container_reader().unwrap();
            let hdr_end = rdr.get_ref().position() as usize; // right after 4-byte container header CRC32
            let mut crafted = base_cram_bytes.clone();
            crafted[26..30].copy_from_slice(&((128 * 1024 * 1024 + 1) as i32).to_le_bytes());
            let mut crc = flate2::Crc::new();
            crc.update(&crafted[26..hdr_end - 4]);
            crafted[hdr_end - 4..hdr_end].copy_from_slice(&crc.sum().to_le_bytes());
            std::fs::write(&huge_hdr_cram, &crafted).unwrap();
            noodles::cram::crai::fs::write(huge_hdr_cram.with_extension("cram.crai"), &base_crai)
                .unwrap();
        }
        assert_io_invalid_data(
            acquire_cram_meta(&huge_hdr_cram, &cram_cfg)
                .expect_err("CRAM FileHeader container > 128 MiB must error"),
            "cram_huge_file_header",
        );

        // Spoof `SliceHeader` coordinates (`ref_seq_id = -3`, `align_start = 0`, `align_span = 0`,
        // and out-of-FASTA `align_start = 2000, align_span = 10`) on a 0-record slice while `.crai` claims `50..=59`.
        {
            let mut rdr = noodles::cram::io::Reader::new(std::fs::File::open(&base_cram).unwrap());
            let _ = rdr.read_header().unwrap();
            let c_pos = rdr.position().unwrap() as usize;
            let mut orig_cont = noodles::cram::io::reader::Container::default();
            let c_len = rdr.read_container(&mut orig_cont).unwrap();
            let body_start =
                std::io::Seek::stream_position(rdr.get_mut()).unwrap() as usize - c_len;
            let lm0 = orig_cont.header().landmarks()[0];
            let comp_hdr_bytes = &base_cram_bytes[body_start..body_start + lm0];

            let encode_itf8 = |v: i32, out: &mut Vec<u8>| {
                let u = v as u32;
                if (0..0x80).contains(&v) {
                    out.push(u as u8);
                } else if (0..0x4000).contains(&v) {
                    out.push(0x80 | ((u >> 8) as u8));
                    out.push((u & 0xff) as u8);
                } else {
                    out.push(0xf0 | ((u >> 28) as u8));
                    out.push(((u >> 20) & 0xff) as u8);
                    out.push(((u >> 12) & 0xff) as u8);
                    out.push(((u >> 4) & 0xff) as u8);
                    out.push((u & 0x0f) as u8);
                }
            };

            for (label, raw_ref_id, raw_start, raw_span) in [
                ("cram_slice_ref_neg3", -3i32, 1i32, 10i32),
                ("cram_slice_start_zero", 0, 0, 10),
                ("cram_slice_span_zero", 0, 1, 0),
                ("cram_slice_past_fasta", 0, 2000, 10),
            ] {
                let bad_cram = dir.path().join(format!("{label}.cram"));
                let mut sh_payload = Vec::new();
                encode_itf8(raw_ref_id, &mut sh_payload);
                encode_itf8(raw_start, &mut sh_payload);
                encode_itf8(raw_span, &mut sh_payload);
                sh_payload.extend_from_slice(&[
                    0x00, // record_count = 0
                    0x00, // record_counter = 0
                    0x01, // block_count = 1
                    0x00, // block_content_ids len = 0
                    0xff, 0xff, 0xff, 0xff, 0x0f, // embedded_ref = -1
                ]);
                sh_payload.extend_from_slice(&[0u8; 16]);

                let mut sh_block = vec![
                    0x00,
                    0x02,
                    0x00,
                    sh_payload.len() as u8,
                    sh_payload.len() as u8,
                ];
                sh_block.extend_from_slice(&sh_payload);
                let mut sh_crc = flate2::Crc::new();
                sh_crc.update(&sh_block);
                sh_block.extend_from_slice(&sh_crc.sum().to_le_bytes());

                let mut core_blk = vec![0x00, 0x05, 0x00, 0x00, 0x00];
                let mut c_crc = flate2::Crc::new();
                c_crc.update(&core_blk);
                core_blk.extend_from_slice(&c_crc.sum().to_le_bytes());

                let mut cont_body = Vec::new();
                cont_body.extend_from_slice(comp_hdr_bytes);
                let new_lm0 = cont_body.len();
                cont_body.extend_from_slice(&sh_block);
                cont_body.extend_from_slice(&core_blk);

                let mut hdr_body = Vec::new();
                hdr_body.extend_from_slice(&(cont_body.len() as i32).to_le_bytes());
                hdr_body.extend_from_slice(&[0x00, 0x32, 0x0a, 0x00, 0x00, 0x00, 0x03, 0x01]);
                encode_itf8(new_lm0 as i32, &mut hdr_body);
                let mut h_crc = flate2::Crc::new();
                h_crc.update(&hdr_body);
                hdr_body.extend_from_slice(&h_crc.sum().to_le_bytes());

                let mut crafted = base_cram_bytes[..c_pos].to_vec();
                crafted.extend_from_slice(&hdr_body);
                crafted.extend_from_slice(&cont_body);
                std::fs::write(&bad_cram, &crafted).unwrap();
                let crai_rec = vec![noodles::cram::crai::Record::new(
                    Some(0),
                    Position::try_from(50usize).ok(),
                    10,
                    c_pos as u64,
                    new_lm0 as u64,
                    (sh_block.len() + core_blk.len()) as u64,
                )];
                noodles::cram::crai::fs::write(bad_cram.with_extension("cram.crai"), &crai_rec)
                    .unwrap();

                let slice_err = query_region_iter(&bad_cram, "chr1", Pos0(0), Pos0(200), &cram_cfg)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .expect_err(&format!("{label}: expected CRAM SliceHeader rejection"));
                assert_io_invalid_data(slice_err, label);
            }
        }
    }

    #[test]
    fn cram_container_cache_coalesces_concurrent_decodes() {
        use noodles::core::Position;
        use noodles::sam::{
            alignment::{
                io::Write as _,
                record::{
                    cigar::{op::Kind, Op},
                    Flags, MappingQuality,
                },
                record_buf::{Cigar, QualityScores as BufQual, RecordBuf, Sequence as BufSeq},
            },
            header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap},
            Header,
        };
        use std::sync::{Arc, Barrier, Mutex};

        let dir = tempfile::tempdir().unwrap();
        let ref_fa = dir.path().join("sf_ref.fa");
        let ref_fai = dir.path().join("sf_ref.fa.fai");
        let cram_a = dir.path().join("sf_a.cram");
        let cram_b = dir.path().join("sf_b.cram");

        let ref_len = 2048usize;
        let ref_seq: Vec<u8> = (0..ref_len).map(|i| b"ACGT"[(i * 5 + 1) % 4]).collect();
        let mut fa_bytes = b">chr1\n".to_vec();
        fa_bytes.extend_from_slice(&ref_seq);
        fa_bytes.push(b'\n');
        std::fs::write(&ref_fa, &fa_bytes).unwrap();
        std::fs::write(
            &ref_fai,
            format!("chr1\t{ref_len}\t6\t{ref_len}\t{}\n", ref_len + 1),
        )
        .unwrap();

        let write_test_cram = |cram_path: &Path, read_count: usize| -> u64 {
            let rs = HeaderMap::<RefSeqMap>::builder()
                .set_length(ref_len.try_into().unwrap())
                .build()
                .unwrap();
            let header = Header::builder().add_reference_sequence("chr1", rs).build();
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(std::fs::File::open(&ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut writer = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(cram_path)
                .unwrap();
            writer.write_alignment_header(&header).unwrap();
            for i in 0..read_count {
                let start_1 = 10 + (i * 7) % 1000;
                let rec = RecordBuf::builder()
                    .set_name(format!("sf_r{i}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(start_1).unwrap())
                    .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
                    .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect::<Cigar>())
                    .set_sequence(BufSeq::from(
                        ref_seq[start_1 - 1..start_1 - 1 + 100].to_vec(),
                    ))
                    .set_quality_scores(BufQual::from(vec![35u8; 100]))
                    .build();
                writer.write_alignment_record(&header, &rec).unwrap();
            }
            writer.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(cram_path).unwrap();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
            crai[0].offset()
        };

        let container_offset = write_test_cram(&cram_a, 128);
        let _ = write_test_cram(&cram_b, 8);

        let cfg = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            ..Default::default()
        };
        clear_bam_caches();
        NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, std::sync::atomic::Ordering::Relaxed);

        // 1. Verify 2-slot `CRAM_META_TLS` retains both paired somatic CRAM handles (`T -> N -> T -> N`),
        //    retrying if a sibling parallel test bumps `BAM_CACHE_EPOCH` mid-cold-open.
        let (meta_a, meta_b, cur_epoch) = loop {
            let epoch_start = BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire);
            let ma = acquire_cram_meta(&cram_a, &cfg).unwrap();
            let mb = acquire_cram_meta(&cram_b, &cfg).unwrap();
            let cur = BAM_CACHE_EPOCH.load(std::sync::atomic::Ordering::Acquire);
            if cur == epoch_start {
                break (ma, mb, cur);
            }
        };
        assert!(Arc::ptr_eq(
            &lookup_cram_meta_tls(&cram_a, &ref_fa, cur_epoch).expect("slot hit A"),
            &meta_a
        ));
        assert!(Arc::ptr_eq(
            &lookup_cram_meta_tls(&cram_b, &ref_fa, cur_epoch).expect("slot hit B"),
            &meta_b
        ));
        assert!(Arc::ptr_eq(
            &lookup_cram_meta_tls(&cram_a, &ref_fa, cur_epoch).expect("slot hit A again"),
            &meta_a
        ));

        // 2. Verify RAII `InflightCramGuard::drop` (`Failed` state) wakes blocked waiters cleanly,
        //    allowing a waiter to retry and decode once without hanging or caching errors.
        let retry_cache: Arc<[Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS]> =
            Arc::new(std::array::from_fn(|_| {
                Mutex::new(CramContainerCacheShard::default())
            }));
        let shard_idx =
            cache_shard_index::<CRAM_CONTAINER_CACHE_SHARDS>(meta_a.id, container_offset);
        let key = (meta_a.id, container_offset);
        let fake_inflight = Arc::new(InflightCramDecode::new());
        retry_cache[shard_idx]
            .lock()
            .unwrap()
            .inflight
            .push((key, Arc::clone(&fake_inflight)));
        let cram_a_retry = cram_a.clone();
        let meta_a_retry = Arc::clone(&meta_a);
        let retry_cache_w = Arc::clone(&retry_cache);
        let entered = Arc::new(Barrier::new(2));
        let entered_w = Arc::clone(&entered);
        let waiter_handle = std::thread::spawn(move || {
            entered_w.wait();
            fetch_cram_container_in(
                &cram_a_retry,
                &meta_a_retry,
                container_offset,
                b"chr1",
                None,
                &retry_cache_w,
                MAX_CACHED_CRAM_CONTAINER_BYTES,
            )
            .unwrap()
            .expect("decoded after Failed wakeup")
        });
        entered.wait();
        std::thread::sleep(std::time::Duration::from_millis(10));
        drop(InflightCramGuard {
            shard_mutex: &retry_cache[shard_idx],
            key,
            entry: fake_inflight,
            completed: false,
        });
        let recovered = waiter_handle.join().unwrap();
        assert_eq!(recovered.len(), 128);

        // 3. Spawn 8 concurrent worker threads synced on an `Arc<Barrier>` hitting the same
        //    uncached `(meta_a, container_offset)`: all 8 receive `Arc::ptr_eq` blocks with 1 decode.
        let shared_cache: Arc<[Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS]> =
            Arc::new(std::array::from_fn(|_| {
                Mutex::new(CramContainerCacheShard::default())
            }));
        let pinned_slot = meta_a.fasta_repo.warm_contig_slot(b"chr1");
        let barrier = Arc::new(Barrier::new(8));
        let mut handles = Vec::with_capacity(8);
        for _ in 0..8 {
            let b = Arc::clone(&barrier);
            let cache_cl = Arc::clone(&shared_cache);
            let meta_cl = Arc::clone(&meta_a);
            let path_cl = cram_a.clone();
            let slot_cl = pinned_slot.clone();
            handles.push(std::thread::spawn(move || {
                b.wait();
                fetch_cram_container_in(
                    &path_cl,
                    &meta_cl,
                    container_offset,
                    b"chr1",
                    Some(&slot_cl),
                    &cache_cl,
                    MAX_CACHED_CRAM_CONTAINER_BYTES,
                )
                .unwrap()
                .expect("concurrent CRAM container decode")
            }));
        }
        let results: Vec<CachedCramBlock> =
            handles.into_iter().map(|h| h.join().unwrap()).collect();
        assert_eq!(results.len(), 8);
        assert_eq!(results[0].len(), 128);
        for block in &results[1..] {
            assert!(
                Arc::ptr_eq(&results[0], block),
                "all 8 concurrent threads must share the exact singleflight CachedCramBlock Arc"
            );
        }
        let total_decodes: usize = shared_cache.iter().map(|s| s.lock().unwrap().decodes).sum();
        assert_eq!(
            total_decodes, 1,
            "8 concurrent fetches on the same uncached CRAM container must perform exactly 1 decode"
        );
    }

    #[test]
    fn cram_4x4_shard_byte_budget_lru_eviction_and_somatic_contig_pin_guard() {
        use noodles::core::Position;
        use noodles::sam::{
            alignment::{
                io::Write as _,
                record::{
                    cigar::{op::Kind, Op},
                    Flags, MappingQuality,
                },
                record_buf::{Cigar, QualityScores as BufQual, RecordBuf, Sequence as BufSeq},
            },
            header::record::value::{map::ReferenceSequence as RefSeqMap, Map as HeaderMap},
            Header,
        };

        assert_eq!(CRAM_CONTAINER_CACHE_SHARDS, 4);
        assert_eq!(CRAM_CONTAINER_SLOTS_PER_SHARD, 4);
        assert_eq!(MAX_SHARD_CACHED_CRAM_BYTES, 40 * 1024 * 1024);
        assert_eq!(DEFAULT_FASTA_CACHE_BUDGET_BYTES, 512 * 1024 * 1024);

        // 1. Per-shard byte budget (`40 MiB` across `4` slots) evicts LRU entries until
        //    total cached bytes in the shard `<= MAX_SHARD_CACHED_CRAM_BYTES`.
        let chr1_arc: std::sync::Arc<str> = std::sync::Arc::from("chr1");
        let build_block = |read_count: usize, read_len: usize| -> CachedCramBlock {
            let mut out = Vec::with_capacity(read_count);
            for idx in 0..read_count {
                let mut sam = SamRecord::new();
                sam.name = format!("r{idx:05}").into_bytes();
                sam.name.shrink_to_fit();
                sam.cigar = vec![(read_len as u32, CigarKind::Match)];
                sam.sequence = vec![b'A'; read_len];
                sam.quality_scores = vec![35u8; read_len];
                out.push(CachedCramRecord {
                    reference_sequence_id: Some(0),
                    alignment_start: Position::new(10),
                    alignment_end: Position::new(10 + read_len - 1),
                    reference_name: Some(std::sync::Arc::clone(&chr1_arc)),
                    mate_reference_name: None,
                    sam,
                });
            }
            out.shrink_to_fit();
            std::sync::Arc::new(CachedCramContainer::from_records(out))
        };

        // Choose block sizes each `<= MAX_CACHED_CRAM_CONTAINER_BYTES` (16 MiB):
        // b1 ~ 10.5 MiB, b2 ~ 14.5 MiB, b3 ~ 14.5 MiB (sum ~ 39.5 MiB <= 40 MiB),
        // then touch b1 (MRU) and insert b4 ~ 15.5 MiB: evicting only b2 leaves ~ 40.5 MiB > 40 MiB,
        // so BOTH b2 and b3 must be evicted while promoted b1 and new b4 remain!
        let b1 = build_block(1_024, 5_200); // ~10.41 MiB
        let b2 = build_block(1_024, 7_250); // ~14.41 MiB
        let b3 = build_block(1_024, 7_250); // ~14.41 MiB
        let b4 = build_block(1_024, 7_900); // ~15.68 MiB
        for b in [&b1, &b2, &b3, &b4] {
            assert!(should_cache_cram_container(
                b,
                MAX_CACHED_CRAM_CONTAINER_BYTES
            ));
        }
        let s1 = approx_cram_container_bytes(&b1);
        let s2 = approx_cram_container_bytes(&b2);
        let s3 = approx_cram_container_bytes(&b3);
        let s4 = approx_cram_container_bytes(&b4);
        assert!(s1 + s2 + s3 <= MAX_SHARD_CACHED_CRAM_BYTES);
        assert!(s1 + s3 + s4 > MAX_SHARD_CACHED_CRAM_BYTES);
        assert!(s1 + s4 <= MAX_SHARD_CACHED_CRAM_BYTES);

        let mut shard = CramContainerCacheShard::default();
        shard.insert_with_budget(
            (1, 100),
            std::sync::Arc::clone(&b1),
            s1,
            MAX_SHARD_CACHED_CRAM_BYTES,
        );
        shard.insert_with_budget(
            (1, 200),
            std::sync::Arc::clone(&b2),
            s2,
            MAX_SHARD_CACHED_CRAM_BYTES,
        );
        shard.insert_with_budget(
            (1, 300),
            std::sync::Arc::clone(&b3),
            s3,
            MAX_SHARD_CACHED_CRAM_BYTES,
        );
        assert_eq!(shard.lru.entries.len(), 3);
        assert_eq!(shard.cached_bytes, s1 + s2 + s3);
        assert!(shard.cached_bytes <= MAX_SHARD_CACHED_CRAM_BYTES);

        // Promote (1, 100) to MRU so LRU order is [(1, 200), (1, 300), (1, 100)].
        assert!(shard.lru.get((1, 100)).is_some());
        shard.insert_with_budget(
            (1, 400),
            std::sync::Arc::clone(&b4),
            s4,
            MAX_SHARD_CACHED_CRAM_BYTES,
        );
        assert!(
            shard.cached_bytes <= MAX_SHARD_CACHED_CRAM_BYTES,
            "shard cached_bytes ({}) must stay <= MAX_SHARD_CACHED_CRAM_BYTES ({MAX_SHARD_CACHED_CRAM_BYTES})",
            shard.cached_bytes
        );
        assert_eq!(shard.cached_bytes, s1 + s4);
        assert_eq!(shard.lru.entries.len(), 2);
        assert!(
            shard.lru.get((1, 200)).is_none(),
            "oldest LRU entry b2 must be evicted"
        );
        assert!(
            shard.lru.get((1, 300)).is_none(),
            "second LRU entry b3 must also be evicted to meet 40 MiB cap"
        );
        assert!(
            shard.lru.get((1, 100)).is_some(),
            "MRU-promoted b1 must remain cached"
        );
        assert!(
            shard.lru.get((1, 400)).is_some(),
            "newly inserted b4 must be cached"
        );
        drop((b1, b2, b3, b4, shard));

        // 2. Verify `pin_cram_contig_for_spec` keeps `tile.chr` pinned across consecutive
        //    `CramFeed` drops (`bam1 -> bam2 -> combineAnalysis`) until the guard drops,
        //    and `cram_record_to_cached` shrinks `SamRecord` buffers to exact length.
        let dir = tempfile::tempdir().unwrap();
        let ref_fa = dir.path().join("som_pin.fa");
        let ref_fai = dir.path().join("som_pin.fa.fai");
        let tumor_cram = dir.path().join("tumor_pin.cram");
        let normal_cram = dir.path().join("normal_pin.cram");

        let contig_len = 600usize;
        let seq1 = vec![b'A'; contig_len];
        let seq2 = vec![b'C'; contig_len];
        std::fs::write(
            &ref_fa,
            format!(
                ">chr1\n{}\n>chr2\n{}\n",
                std::str::from_utf8(&seq1).unwrap(),
                std::str::from_utf8(&seq2).unwrap()
            ),
        )
        .unwrap();
        std::fs::write(
            &ref_fai,
            format!(
                "chr1\t{contig_len}\t6\t{contig_len}\t{}\nchr2\t{contig_len}\t{}\t{contig_len}\t{}\n",
                contig_len + 1,
                6 + contig_len + 1 + 6,
                contig_len + 1
            ),
        )
        .unwrap();

        let write_one_cram = |cram_path: &Path, qname: &str| -> u64 {
            let mut hb = Header::builder();
            for name in ["chr1", "chr2"] {
                let rs = HeaderMap::<RefSeqMap>::builder()
                    .set_length(contig_len.try_into().unwrap())
                    .build()
                    .unwrap();
                hb = hb.add_reference_sequence(name, rs);
            }
            let header = hb.build();
            let fasta_reader = noodles::fasta::io::IndexedReader::new(
                std::io::BufReader::new(std::fs::File::open(&ref_fa).unwrap()),
                noodles::fasta::fai::fs::read(&ref_fai).unwrap(),
            );
            let repo = noodles::fasta::Repository::new(
                noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader),
            );
            let mut cw = noodles::cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(cram_path)
                .unwrap();
            cw.write_alignment_header(&header).unwrap();
            let rec = RecordBuf::builder()
                .set_name(qname)
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(100usize).unwrap())
                .set_mapping_quality(MappingQuality::try_from(60u8).unwrap())
                .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect::<Cigar>())
                .set_sequence(BufSeq::from(vec![b'A'; 100]))
                .set_quality_scores(BufQual::from(vec![35u8; 100]))
                .build();
            cw.write_alignment_record(&header, &rec).unwrap();
            cw.try_finish(&header).unwrap();
            let crai = noodles::cram::fs::index(cram_path).unwrap();
            let off = crai[0].offset();
            noodles::cram::crai::fs::write(cram_path.with_extension("cram.crai"), &crai).unwrap();
            off
        };

        let t_off = write_one_cram(&tumor_cram, "tumor_r1");
        let _ = write_one_cram(&normal_cram, "normal_r1");
        let bam_spec = format!("{}|{}", tumor_cram.display(), normal_cram.display());
        let cfg = Config {
            bam: bam_spec.clone(),
            reference: ref_fa.to_string_lossy().into_owned(),
            ..Default::default()
        };

        NEXT_BAM_META_ID.fetch_max(TEST_META_ID_BASE, std::sync::atomic::Ordering::Relaxed);
        let meta_t = acquire_cram_meta(&tumor_cram, &cfg).unwrap();
        let per_contig_bytes = contig_len + CRAM_REF_OVERHANG_PAD;
        // Budget holds at most 1 contig (`chr1` OR `chr2`).
        meta_t.fasta_repo.set_budget_bytes(per_contig_bytes);

        // Verify `cram_record_to_cached` shrinks `name`/`cigar`/`sequence`/`quality_scores` to exact length.
        let shrink_cache: [std::sync::Mutex<CramContainerCacheShard>; CRAM_CONTAINER_CACHE_SHARDS] =
            std::array::from_fn(|_| std::sync::Mutex::new(CramContainerCacheShard::default()));
        let decoded_block = fetch_cram_container_in(
            &tumor_cram,
            &meta_t,
            t_off,
            b"chr1",
            None,
            &shrink_cache,
            MAX_CACHED_CRAM_CONTAINER_BYTES,
        )
        .unwrap()
        .unwrap();
        assert_eq!(decoded_block.len(), 1);
        let cached_sam = &decoded_block[0].sam;
        assert_eq!(cached_sam.name.capacity(), cached_sam.name.len());
        assert_eq!(cached_sam.cigar.capacity(), cached_sam.cigar.len());
        assert_eq!(cached_sam.sequence.capacity(), cached_sam.sequence.len());
        assert_eq!(
            cached_sam.quality_scores.capacity(),
            cached_sam.quality_scores.len()
        );

        // Simulate a concurrent worker holding `chr2` pinned (`resident_bytes > budget`):
        // 1) Baseline without tile guard: dropping `feed1` unpins and evicts `chr1`.
        let chr2_worker_pin = meta_t.fasta_repo.warm_contig_slot(b"chr2");
        {
            let mut feed1 =
                query_region_iter(&tumor_cram, "chr1", Pos0(50), Pos0(250), &cfg).unwrap();
            let mut rec = SamRecord::new();
            assert!(feed1.next_into(&mut rec).unwrap());
            drop(feed1);
        }
        assert!(
            meta_t.fasta_repo.resident_slot(b"chr1").is_none(),
            "without tile guard, CramFeedInner::drop must evict chr1 when over budget"
        );

        // 2) With `pin_cram_contig_for_spec`: `chr1` stays pinned across `tumor` feed drop AND
        //    `normal` feed drop (exact same `ContigSlot` Arc, zero chr1 reload) and evicts cleanly
        //    once `tile_pin` drops.
        let tile_pin = pin_cram_contig_for_spec(&bam_spec, "chr1", &cfg);
        let initial_chr1_slot = meta_t
            .fasta_repo
            .resident_slot(b"chr1")
            .expect("chr1 must be warmed and pinned by pin_cram_contig_for_spec");

        {
            let mut feed_t =
                query_region_iter(&tumor_cram, "chr1", Pos0(50), Pos0(250), &cfg).unwrap();
            let mut rec = SamRecord::new();
            assert!(feed_t.next_into(&mut rec).unwrap());
            drop(feed_t);
        }
        let after_t_slot = meta_t
            .fasta_repo
            .resident_slot(b"chr1")
            .expect("tile guard must keep chr1 pinned after BAM1 CramFeed drops");
        assert!(std::sync::Arc::ptr_eq(&initial_chr1_slot, &after_t_slot));
        drop(after_t_slot);

        let loads_before_n = meta_t
            .fasta_repo
            .load_count
            .load(std::sync::atomic::Ordering::SeqCst);
        {
            let mut feed_n =
                query_region_iter(&normal_cram, "chr1", Pos0(50), Pos0(250), &cfg).unwrap();
            let mut rec = SamRecord::new();
            assert!(feed_n.next_into(&mut rec).unwrap());
            drop(feed_n);
        }
        assert_eq!(
            meta_t
                .fasta_repo
                .load_count
                .load(std::sync::atomic::Ordering::SeqCst),
            loads_before_n,
            "BAM2 CramFeed on pinned chr1 must incur 0 FASTA loads"
        );
        let after_n_slot = meta_t
            .fasta_repo
            .resident_slot(b"chr1")
            .expect("tile guard must keep chr1 pinned after BAM2 CramFeed drops");
        assert!(std::sync::Arc::ptr_eq(&initial_chr1_slot, &after_n_slot));
        drop((initial_chr1_slot, after_n_slot));

        assert!(meta_t.fasta_repo.resident_slot(b"chr1").is_some());
        assert!(meta_t.fasta_repo.resident_slot(b"chr2").is_some());
        drop(tile_pin);
        assert!(
            meta_t.fasta_repo.resident_slot(b"chr1").is_none(),
            "dropping CramContigPinGuard must unpin and evict chr1 under budget pressure"
        );
        drop(chr2_worker_pin);
    }
}
