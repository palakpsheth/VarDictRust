//! CIGAR processing: read-level CIGAR rewrite and pileup engine.
//!
//! - `modify`: ports CigarModifier.java (ordered fixpoint, bounds-checked)
//! - `parse`: ports CigarParser.java (pileup engine)

pub mod modify;
pub mod parse;

use noodles::sam::alignment::record::cigar::op::Kind;

/// Parse a CIGAR string (e.g. "10M2I5M") into `(len, op)` pairs.
///
/// The op kind is noodles' [`Kind`], the same representation used by decoded
/// BAM records (`SamRecord::cigar`). Mapping to chars is SAM spec Table 1,
/// 1:1 with noodles' own (crate-private) text parser:
/// `M I D N S H P = X`.
pub fn parse_cigar(s: &str) -> Option<Vec<(u32, Kind)>> {
    let mut elements = Vec::new();
    parse_cigar_into(s, &mut elements)?;
    Some(elements)
}

/// Parse a CIGAR string into a reused buffer (scratch reuse; identical
/// elements to [`parse_cigar`], `None` on the same malformed inputs).
pub fn parse_cigar_into(s: &str, elements: &mut Vec<(u32, Kind)>) -> Option<()> {
    elements.clear();
    let bytes = s.as_bytes();
    let mut num_start = 0;

    for (i, &b) in bytes.iter().enumerate() {
        if b.is_ascii_digit() {
            continue;
        }
        // b should be a CIGAR operation character
        let kind = match b {
            b'M' => Kind::Match,
            b'I' => Kind::Insertion,
            b'D' => Kind::Deletion,
            b'N' => Kind::Skip,
            b'S' => Kind::SoftClip,
            b'H' => Kind::HardClip,
            b'P' => Kind::Pad,
            b'=' => Kind::SequenceMatch,
            b'X' => Kind::SequenceMismatch,
            _ => return None,
        };
        if i == num_start {
            // No length prefix — invalid
            return None;
        }
        let len_str = std::str::from_utf8(&bytes[num_start..i]).ok()?;
        let len: u32 = len_str.parse().ok()?;
        elements.push((len, kind));
        num_start = i + 1;
    }

    // If there are trailing digits with no operation, the CIGAR is malformed
    if num_start != bytes.len() {
        return None;
    }

    Some(())
}

/// A modified cigar record, mirroring VarDictJava's `ModifiedCigar` class.
///
/// The CigarModifier may adjust the alignment start position and trim the
/// query sequence/quality alongside the CIGAR string.
#[derive(Debug, Clone)]
pub struct ModifiedCigar {
    /// Adjusted alignment start position (0-based).
    pub position: i32,
    /// Modified CIGAR string.
    pub cigar: String,
    /// Modified query sequence (as bytes).
    pub query_seq: Vec<u8>,
    /// Modified query quality scores (as bytes, Phred+33).
    pub query_qual: Vec<u8>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_valid() {
        assert_eq!(
            parse_cigar("50M2I49M"),
            Some(vec![
                (50, Kind::Match),
                (2, Kind::Insertion),
                (49, Kind::Match)
            ])
        );
        assert_eq!(
            parse_cigar("5H10S30M"),
            Some(vec![
                (5, Kind::HardClip),
                (10, Kind::SoftClip),
                (30, Kind::Match)
            ])
        );
        assert_eq!(
            parse_cigar("50=2X49M"),
            Some(vec![
                (50, Kind::SequenceMatch),
                (2, Kind::SequenceMismatch),
                (49, Kind::Match)
            ])
        );
        assert_eq!(parse_cigar("").unwrap(), Vec::new()); // empty is valid (no elements)
    }

    #[test]
    fn parse_invalid() {
        assert!(parse_cigar("M").is_none()); // no length
        assert!(parse_cigar("10").is_none()); // trailing digits
        assert!(parse_cigar("10Z").is_none()); // invalid op
    }
}
