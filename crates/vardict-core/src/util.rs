//! Utility functions — bounds-checked sequence helpers.
//!
//! These replace Java's `String.charAt(i)` / `String.substring(i, j)` with
//! safe alternatives that return `Option` instead of panicking on out-of-bounds.

/// Bounds-checked byte access, supporting negative indices (Python-style).
///
/// Negative indices count from the end: `char_at(b"ACGT", -1)` → `Some(b'T')`.
/// Returns `None` if the index is out of bounds.
pub fn char_at(seq: &[u8], idx: isize) -> Option<u8> {
    let len = seq.len() as isize;
    let resolved = if idx < 0 { len + idx } else { idx };
    if resolved < 0 || resolved >= len {
        None
    } else {
        Some(seq[resolved as usize])
    }
}

/// Bounds-checked subsequence extraction.
///
/// Returns `None` if `start` is out of bounds. If `start + len` exceeds the
/// slice length, returns the available suffix (matching Java's lenient behavior
/// in many VarDict paths).
pub fn subseq(seq: &[u8], start: usize, len: usize) -> Option<&[u8]> {
    if start > seq.len() {
        return None;
    }
    let end = (start + len).min(seq.len());
    Some(&seq[start..end])
}

/// Reverse-complement a DNA byte sequence in place.
pub fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = complement(*base);
    }
}

/// Complement a single DNA base.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        _ => b'N',
    }
}

/// Collect reference bases over the 1-based inclusive range `from..=to`.
///
/// Positions < 1 are skipped; positions missing from the reference render as
/// `missing` (`None` skips them entirely).
pub fn ref_bases_with(
    reference: &crate::reference::Reference,
    chr: &str,
    from: i32,
    to: i32,
    missing: Option<u8>,
) -> Vec<u8> {
    let mut sb = Vec::new();
    for pos in from..=to {
        if pos > 0 {
            match reference.get_base(chr, pos as u32) {
                Some(b) => sb.push(b),
                None => {
                    if let Some(fill) = missing {
                        sb.push(fill);
                    }
                }
            }
        }
    }
    sb
}

/// Like [`ref_bases_with`] but skipping positions missing from the reference.
pub fn ref_bases(
    reference: &crate::reference::Reference,
    chr: &str,
    from: i32,
    to: i32,
) -> Vec<u8> {
    ref_bases_with(reference, chr, from, to, None)
}

/// Helper to construct a reference sequence from start to end position
/// (inclusive), as a lossy UTF-8 string.
///
/// Under 1-based coordinates.
pub fn join_ref(
    reference: &crate::reference::Reference,
    chr: &str,
    start: i32,
    end: i32,
) -> String {
    String::from_utf8_lossy(&ref_bases(reference, chr, start, end)).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_char_at_positive() {
        let seq = b"ACGT";
        assert_eq!(char_at(seq, 0), Some(b'A'));
        assert_eq!(char_at(seq, 3), Some(b'T'));
        assert_eq!(char_at(seq, 4), None);
    }

    #[test]
    fn test_char_at_negative() {
        let seq = b"ACGT";
        assert_eq!(char_at(seq, -1), Some(b'T'));
        assert_eq!(char_at(seq, -4), Some(b'A'));
        assert_eq!(char_at(seq, -5), None);
    }

    #[test]
    fn test_char_at_empty() {
        let seq: &[u8] = b"";
        assert_eq!(char_at(seq, 0), None);
        assert_eq!(char_at(seq, -1), None);
    }

    #[test]
    fn test_subseq_normal() {
        let seq = b"ACGTACGT";
        assert_eq!(subseq(seq, 2, 3), Some(b"GTA".as_ref()));
    }

    #[test]
    fn test_subseq_overflow() {
        let seq = b"ACGT";
        // Request goes past end — returns available suffix
        assert_eq!(subseq(seq, 2, 10), Some(b"GT".as_ref()));
    }

    #[test]
    fn test_subseq_out_of_bounds() {
        let seq = b"ACGT";
        assert_eq!(subseq(seq, 5, 1), None);
    }

    #[test]
    fn test_reverse_complement() {
        let mut seq = b"ACGT".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(&seq, b"ACGT"); // palindrome
        let mut seq2 = b"AACG".to_vec();
        reverse_complement(&mut seq2);
        assert_eq!(&seq2, b"CGTT");
    }
}
