//! Utility functions — bounds-checked sequence helpers.
//!
//! These replace Java's `String.charAt(i)` / `String.substring(i, j)` with
//! safe alternatives that return `Option` instead of panicking on out-of-bounds.

/// Compile a static regex literal once (`LazyLock` / `OnceLock` initializers).
///
/// Centralizes the single sanctioned `.unwrap()` on compile-time regex patterns
/// so `vardictrust-core` enforces the no-panic lint ratchet across all runtime
/// library code (`AGENTS.md` §9).
#[allow(
    clippy::unwrap_used,
    reason = "static regex literal validated by unit tests"
)]
#[inline]
pub(crate) fn static_regex(pattern: &'static str) -> regex::Regex {
    regex::Regex::new(pattern).unwrap()
}

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

/// Best-effort diagnostic output to `stderr` under `-y` / `--verbose`,
/// ignoring write errors instead of panicking (`AGENTS.md` §9).
#[macro_export]
macro_rules! verbose_eprintln {
    ($($arg:tt)*) => {{
        use std::io::Write as _;
        let _ = writeln!(std::io::stderr(), $($arg)*);
    }};
}

/// Reverse-complement a DNA byte sequence in place.
pub fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = complement(*base);
    }
}

/// Reverse-complement an exact `N`-byte DNA slice into a stack array without heap allocation,
/// returning `None` when `slice.len() != N`. Matches `cigar/modify.rs::rc_seq` byte-for-byte.
#[allow(
    dead_code,
    reason = "shared Wave 0 primitive for cigar/lib seed callers"
)]
#[inline]
pub(crate) fn reverse_complement_bytes<const N: usize>(slice: &[u8]) -> Option<[u8; N]> {
    if slice.len() != N {
        return None;
    }
    let mut out = [0u8; N];
    out.copy_from_slice(slice);
    reverse_complement(&mut out);
    Some(out)
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
    if to < from || to < 1 {
        return Vec::new();
    }
    let eff_from = from.max(1);
    let span = ((to as i64) - (eff_from as i64) + 1).clamp(0, 1 << 20) as usize;
    let mut sb = Vec::with_capacity(span);
    for pos in eff_from..=to {
        match reference.get_base(chr, pos as u32) {
            Some(b) => sb.push(b),
            None => {
                if let Some(fill) = missing {
                    sb.push(fill);
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

/// Writes reference bases over the 1-based inclusive range `start..=end` into `out`,
/// clearing `out` first and reserving the span capacity without intermediate `Vec<u8>`
/// or `String` allocations (`L2`/`L8`).
pub fn ref_bases_into(
    out: &mut String,
    reference: &crate::reference::Reference,
    chr: &str,
    start: i32,
    end: i32,
) {
    out.clear();
    if end < start || end < 1 {
        return;
    }
    let eff_start = start.max(1);
    let span = ((end as i64) - (eff_start as i64) + 1).clamp(0, 1 << 20) as usize;
    if out.capacity() == 0 {
        out.reserve_exact(span);
    } else {
        out.reserve(span);
    }
    for pos in eff_start..=end {
        if let Some(b) = reference.get_base(chr, pos as u32) {
            if b.is_ascii() {
                out.push(b as char);
            } else {
                out.clear();
                out.push_str(&String::from_utf8_lossy(&ref_bases(
                    reference, chr, start, end,
                )));
                return;
            }
        }
    }
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
    let mut out = String::new();
    ref_bases_into(&mut out, reference, chr, start, end);
    out
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

    fn make_test_reference() -> crate::reference::Reference {
        let mut reference = crate::reference::Reference::new();
        let bases = b"ACGTACGTNNACGT".to_vec();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: bases.len() as u32,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), bases.clone());
        reference.region_sequences.push(bases);
        reference
    }

    #[test]
    fn ref_bases_into_matches_legacy() {
        let reference = make_test_reference();
        let mut scratch = String::new();
        for start in -2..=16 {
            for end in -2..=16 {
                ref_bases_into(&mut scratch, &reference, "chr1", start, end);
                let legacy = String::from_utf8_lossy(&ref_bases(&reference, "chr1", start, end))
                    .into_owned();
                assert_eq!(
                    scratch, legacy,
                    "mismatch for range {start}..={end}: got {scratch:?}, expected {legacy:?}"
                );
            }
        }
    }

    #[test]
    fn join_ref_capacity_exact() {
        let reference = make_test_reference();
        let s = join_ref(&reference, "chr1", 2, 10);
        assert_eq!(s, "CGTACGTNN");
        assert_eq!(s.len(), 9);
        assert_eq!(s.capacity(), 9);
    }

    #[test]
    fn reverse_complement_bytes_12_matches_rc_seq_and_checks_length() {
        assert_eq!(
            reverse_complement_bytes::<12>(b"AAACCCGGGTTN"),
            Some(*b"NAACCCGGGTTT")
        );

        for sample in [
            &b"AAACCCGGGTTN"[..],
            &b"AaCcGgTtNnXx"[..],
            &b"acgtACGTnN12"[..],
            &b"ttttggggcccc"[..],
        ] {
            let mut expected = sample.to_vec();
            reverse_complement(&mut expected);
            let got = reverse_complement_bytes::<12>(sample);
            assert_eq!(
                got.as_ref().map(|a| a.as_slice()),
                Some(expected.as_slice()),
                "mismatch for {sample:?}"
            );
        }

        assert_eq!(reverse_complement_bytes::<12>(b""), None);
        assert_eq!(reverse_complement_bytes::<12>(b"AAACCCGGGTT"), None);
        assert_eq!(reverse_complement_bytes::<12>(b"AAACCCGGGTTNN"), None);
    }

    #[test]
    fn verbose_eprintln_macro_formats_without_panic() {
        crate::verbose_eprintln!("verbose_eprintln smoke check: {} {}", "ok", 123);
    }
}
