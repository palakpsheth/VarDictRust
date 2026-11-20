//! CIGAR modify — ports CigarModifier.java (ordered fixpoint, bounds-checked)
//!
//! normalizes indels near soft-clips, combines adjacent small indels into
//! complex variants, and cleans up end-of-read alignment mismatches.

use crate::cigar::ModifiedCigar;
use crate::config::Config;
use crate::reference::SeedMap;
use crate::util::{char_at, reverse_complement, reverse_complement_bytes, subseq};
#[cfg(test)]
use regex::Regex;
use std::fmt::Write as _;
#[cfg(test)]
use std::sync::LazyLock;

// Define all the regexes from Patterns.java (retained solely for `#[cfg(test)]` legacy oracle parity)
#[cfg(test)]
static BEGIN_NUMBER_D: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)D"));
#[cfg(test)]
static END_NUMBER_D: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)D$"));
#[cfg(test)]
static BEGIN_NUMBER_I: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)I"));
#[cfg(test)]
static END_NUMBER_I: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)I$"));

#[cfg(test)]
static SA_CIGAR_D_S_5CLIP_GROUP: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d\d+)S"));
#[cfg(test)]
static SA_CIGAR_D_S_3CLIP_GROUP: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d\d+)S$"));

#[cfg(test)]
static BEGIN_NUMBER_S_NUMBER_IOR_D: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d+)S(\d+)([ID])"));
#[cfg(test)]
static NUMBER_IOR_D_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d+)S$"));
#[cfg(test)]
static BEGIN_NUMBER_S_NUMBER_M_NUMBER_IOR_D: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d+)S(\d+)M(\d+)([ID])"));
#[cfg(test)]
static NUMBER_IOR_D_NUMBER_M_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d+)M(\d+)S$"));

#[cfg(test)]
static BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d)M(\d+)([ID])(\d+)M"));
#[cfg(test)]
static BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M_OPT: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^\dM\d+[ID]\d+M"));
#[cfg(test)]
static NUMBER_IOR_D_DIGIT_M_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d)M$"));
#[cfg(test)]
static NUMBER_IOR_D_NUMBER_M_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d+)M$"));

#[cfg(test)]
static D_M_D_DD_M_D_I_D_M_D_DD: LazyLock<Regex> = LazyLock::new(|| {
    crate::util::static_regex(r"^(.*?)(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M")
});
#[cfg(test)]
static D_M_D_DD_M_D_I_D_M_D_DD_PRIM: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M"));

#[cfg(test)]
static THREE_DELETIONS_PATTERN: LazyLock<Regex> = LazyLock::new(|| {
    crate::util::static_regex(r"^(.*?)(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M")
});
#[cfg(test)]
static THREE_INDELS_PATTERN: LazyLock<Regex> = LazyLock::new(|| {
    crate::util::static_regex(r"^(.*?)(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M")
});

#[cfg(test)]
static DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"\d+M\d+[DI]\d+M\d+[DI]\d+M\d+[DI]\d+M"));
#[cfg(test)]
static DM_DD_DM_DD_DM_DD_DM: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"\d+M\d+D\d+M\d+D\d+M\d+D\d+M"));

#[cfg(test)]
static DIG_D_DIG_M_DIG_DI_DIGI: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)D(\d+)M(\d+)([DI])(\d+I)?"));
#[cfg(test)]
static NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?"));

#[cfg(test)]
static DIG_D_DIG_D: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)D(\d+)D"));
#[cfg(test)]
static DIG_I_DIG_I: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)I(\d+)I"));

#[cfg(test)]
static ANY_NUMBER_M_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(.*?)(\d+)M(\d+)S$"));
#[cfg(test)]
static BEGIN_ANY_DIG_M_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(.*?)(\d+)M$"));
#[cfg(test)]
static DIG_M_END: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"\d+M$"));
#[cfg(test)]
static BEGIN_DIG_M: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)M"));
#[cfg(test)]
static DIG_S_DIG_M: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)S(\d+)M"));
#[cfg(test)]
static DIG_M_DIG_S_END: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"\d+M\d+S$"));

#[cfg(test)]
static ALIGNED_LENGTH_MND: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)[MND]"));
#[cfg(test)]
static SOFT_CLIPPED: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)[MIS]"));
#[cfg(test)]
static RE_START_S: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^\d+S"));
#[cfg(test)]
static RE_END_2S: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"\d\d+S$"));
#[cfg(test)]
static CIGAR_ELEMENT_PATTERN: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\-?\d+)[MIDNSHPX=]"));

const LOWQUAL: u8 = 10;
const SEED_2: usize = 12;
/// Maximum intervening `M` span (`15 bp`) between adjacent indels eligible to be merged into a
/// single complex `D`/`I` CIGAR event (`CigarModifier.java:530-690`, `M-DOCUMENTED-MAGIC`).
/// Changing this threshold would alter multi-indel combination vs VarDictJava v1.8.3.
const MAX_MULTI_INDEL_MID_SPAN: i32 = 15;

/// Zero-allocation 256-bit set for counting distinct `u8` bases in inner loops.
#[derive(Default, Clone, Copy)]
pub(crate) struct ByteSet([u64; 4]);

impl ByteSet {
    #[inline]
    pub(crate) const fn new() -> Self {
        Self([0; 4])
    }

    #[inline]
    pub(crate) fn insert(&mut self, b: u8) {
        if let Some(word) = self.0.get_mut(usize::from(b >> 6)) {
            *word |= 1u64 << (b & 63);
        }
    }

    #[inline]
    pub(crate) const fn len(&self) -> usize {
        let [w0, w1, w2, w3] = self.0;
        (w0.count_ones() + w1.count_ones() + w2.count_ones() + w3.count_ones()) as usize
    }
}

/// Helper to parse matched group as integer with default 0.
#[inline]
fn to_int(s: &str) -> i32 {
    s.parse::<i32>().unwrap_or(0)
}

/// Matches `regex`'s Unicode `\d` (`\p{Decimal_Number}`, Unicode 16.0.0) so `\d` and `\D`
/// cursor checks match `Patterns.java` / `regex::Regex` identically on all UTF-8 inputs.
#[inline]
fn is_regex_decimal_digit(ch: char) -> bool {
    if ch.is_ascii() {
        return ch.is_ascii_digit();
    }
    matches!(
        ch,
        '\u{0660}'..='\u{0669}'
            | '\u{06F0}'..='\u{06F9}'
            | '\u{07C0}'..='\u{07C9}'
            | '\u{0966}'..='\u{096F}'
            | '\u{09E6}'..='\u{09EF}'
            | '\u{0A66}'..='\u{0A6F}'
            | '\u{0AE6}'..='\u{0AEF}'
            | '\u{0B66}'..='\u{0B6F}'
            | '\u{0BE6}'..='\u{0BEF}'
            | '\u{0C66}'..='\u{0C6F}'
            | '\u{0CE6}'..='\u{0CEF}'
            | '\u{0D66}'..='\u{0D6F}'
            | '\u{0DE6}'..='\u{0DEF}'
            | '\u{0E50}'..='\u{0E59}'
            | '\u{0ED0}'..='\u{0ED9}'
            | '\u{0F20}'..='\u{0F29}'
            | '\u{1040}'..='\u{1049}'
            | '\u{1090}'..='\u{1099}'
            | '\u{17E0}'..='\u{17E9}'
            | '\u{1810}'..='\u{1819}'
            | '\u{1946}'..='\u{194F}'
            | '\u{19D0}'..='\u{19D9}'
            | '\u{1A80}'..='\u{1A89}'
            | '\u{1A90}'..='\u{1A99}'
            | '\u{1B50}'..='\u{1B59}'
            | '\u{1BB0}'..='\u{1BB9}'
            | '\u{1C40}'..='\u{1C49}'
            | '\u{1C50}'..='\u{1C59}'
            | '\u{A620}'..='\u{A629}'
            | '\u{A8D0}'..='\u{A8D9}'
            | '\u{A900}'..='\u{A909}'
            | '\u{A9D0}'..='\u{A9D9}'
            | '\u{A9F0}'..='\u{A9F9}'
            | '\u{AA50}'..='\u{AA59}'
            | '\u{ABF0}'..='\u{ABF9}'
            | '\u{FF10}'..='\u{FF19}'
            | '\u{104A0}'..='\u{104A9}'
            | '\u{10D30}'..='\u{10D39}'
            | '\u{10D40}'..='\u{10D49}'
            | '\u{11066}'..='\u{1106F}'
            | '\u{110F0}'..='\u{110F9}'
            | '\u{11136}'..='\u{1113F}'
            | '\u{111D0}'..='\u{111D9}'
            | '\u{112F0}'..='\u{112F9}'
            | '\u{11450}'..='\u{11459}'
            | '\u{114D0}'..='\u{114D9}'
            | '\u{11650}'..='\u{11659}'
            | '\u{116C0}'..='\u{116C9}'
            | '\u{116D0}'..='\u{116E3}'
            | '\u{11730}'..='\u{11739}'
            | '\u{118E0}'..='\u{118E9}'
            | '\u{11950}'..='\u{11959}'
            | '\u{11BF0}'..='\u{11BF9}'
            | '\u{11C50}'..='\u{11C59}'
            | '\u{11D50}'..='\u{11D59}'
            | '\u{11DA0}'..='\u{11DA9}'
            | '\u{11F50}'..='\u{11F59}'
            | '\u{16130}'..='\u{16139}'
            | '\u{16A60}'..='\u{16A69}'
            | '\u{16AC0}'..='\u{16AC9}'
            | '\u{16B50}'..='\u{16B59}'
            | '\u{16D70}'..='\u{16D79}'
            | '\u{1CCF0}'..='\u{1CCF9}'
            | '\u{1D7CE}'..='\u{1D7FF}'
            | '\u{1E140}'..='\u{1E149}'
            | '\u{1E2F0}'..='\u{1E2F9}'
            | '\u{1E4F0}'..='\u{1E4F9}'
            | '\u{1E5F1}'..='\u{1E5FA}'
            | '\u{1E950}'..='\u{1E959}'
            | '\u{1FBF0}'..='\u{1FBF9}'
    )
}

#[cfg(test)]
thread_local! {
    static SCAN_DIGITS_FWD_CALLS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

#[cfg(test)]
fn reset_scan_digits_fwd_calls() {
    SCAN_DIGITS_FWD_CALLS.with(|c| c.set(0));
}

#[cfg(test)]
fn scan_digits_fwd_calls() -> u64 {
    SCAN_DIGITS_FWD_CALLS.with(|c| c.get())
}

#[inline]
fn scan_digits_fwd(s: &str, idx: usize) -> (usize, usize) {
    #[cfg(test)]
    SCAN_DIGITS_FWD_CALLS.with(|c| c.set(c.get().saturating_add(1)));
    let bytes = s.as_bytes();
    let mut end = idx;
    let mut char_count = 0usize;
    while let Some(&b) = bytes.get(end) {
        if b.is_ascii_digit() {
            end += 1;
            char_count += 1;
        } else if b >= 0x80 {
            if let Some(ch) = s.get(end..).and_then(|tail| tail.chars().next()) {
                if is_regex_decimal_digit(ch) {
                    end += ch.len_utf8();
                    char_count += 1;
                    continue;
                }
            }
            break;
        } else {
            break;
        }
    }
    (end, char_count)
}

#[inline]
fn scan_digits_rev(s: &str, end: usize) -> (usize, usize) {
    let bytes = s.as_bytes();
    let mut start = end;
    let mut char_count = 0usize;
    while start > 0 {
        if let Some(&b) = bytes.get(start - 1) {
            if b.is_ascii_digit() {
                start -= 1;
                char_count += 1;
                continue;
            } else if b >= 0x80 {
                if let Some(ch) = s.get(..start).and_then(|head| head.chars().next_back()) {
                    if is_regex_decimal_digit(ch) {
                        start -= ch.len_utf8();
                        char_count += 1;
                        continue;
                    }
                }
            }
        }
        break;
    }
    (start, char_count)
}

/// Forward cursor parsing `<digits><ascii_op>` starting at `idx`, returning `(val, op, digit_chars, next_idx)`.
#[inline]
fn next_cigar_op(s: &str, idx: usize) -> Option<(i32, u8, usize, usize)> {
    let (d_end, digit_chars) = scan_digits_fwd(s, idx);
    if digit_chars == 0 {
        return None;
    }
    let &op = s.as_bytes().get(d_end)?;
    if op >= 0x80 {
        return None;
    }
    let val = to_int(s.get(idx..d_end)?);
    Some((val, op, digit_chars, d_end + 1))
}

/// Reverse cursor parsing `<digits><ascii_op>` ending at `end`, returning `(op_start, val, op, digit_chars)`.
#[inline]
fn prev_cigar_op(s: &str, end: usize) -> Option<(usize, i32, u8, usize)> {
    let op_pos = end.checked_sub(1)?;
    let &op = s.as_bytes().get(op_pos)?;
    if op >= 0x80 {
        return None;
    }
    let (d_start, digit_chars) = scan_digits_rev(s, op_pos);
    if digit_chars == 0 {
        return None;
    }
    let val = to_int(s.get(d_start..op_pos)?);
    Some((d_start, val, op, digit_chars))
}

/// Advances `idx` across a maximal digit run without an ASCII op or across contiguous `\D` characters
/// so leftmost scanners stay strictly $O(N)$ on arbitrary / hostile digit runs (`AR-05`).
#[inline]
fn step_past_current_run(s: &str, idx: usize) -> usize {
    let bytes = s.as_bytes();
    let (d_end, digit_chars) = scan_digits_fwd(s, idx);
    if digit_chars > 0 {
        if bytes.get(d_end).is_some_and(|&b| b < 0x80) {
            d_end + 1
        } else {
            d_end
        }
    } else {
        let mut cur = idx;
        while let Some(&b) = bytes.get(cur) {
            if b < 0x80 {
                if b.is_ascii_digit() {
                    break;
                }
                cur += 1;
            } else if let Some(ch) = s.get(cur..).and_then(|t| t.chars().next()) {
                if is_regex_decimal_digit(ch) {
                    break;
                }
                cur += ch.len_utf8();
            } else {
                cur += 1;
            }
        }
        cur
    }
}

/// Matches `SA_CIGAR_D_S_5CLIP_GROUP` (`^(\d\d+)S`), returning `to_int(group1) as usize`.
#[inline]
pub(crate) fn leading_2digit_s_len(s: &str) -> Option<usize> {
    let (val, op, digit_chars, _) = next_cigar_op(s, 0)?;
    if op != b'S' || digit_chars < 2 {
        return None;
    }
    Some(val as usize)
}

/// Matches `SA_CIGAR_D_S_3CLIP_GROUP` (`(\d\d+)S$`), returning `to_int(group1) as usize`.
#[inline]
pub(crate) fn trailing_2digit_s_len(s: &str) -> Option<usize> {
    let (_, val, op, digit_chars) = prev_cigar_op(s, s.len())?;
    if op != b'S' || digit_chars < 2 {
        return None;
    }
    Some(val as usize)
}

#[inline]
fn match_begin_s_ior_d(s: &str) -> Option<(std::ops::Range<usize>, i32, i32, u8)> {
    let (g1, b'S', _, idx1) = next_cigar_op(s, 0)? else {
        return None;
    };
    let (g2, op, _, end) = next_cigar_op(s, idx1)?;
    if !matches!(op, b'I' | b'D') {
        return None;
    }
    Some((0..end, g1, g2, op))
}

#[inline]
fn match_ior_d_s_end(s: &str) -> Option<(std::ops::Range<usize>, i32, u8, i32)> {
    let (s_start, g3, b'S', _) = prev_cigar_op(s, s.len())? else {
        return None;
    };
    let (start, g1, op, _) = prev_cigar_op(s, s_start)?;
    if !matches!(op, b'I' | b'D') {
        return None;
    }
    Some((start..s.len(), g1, op, g3))
}

#[inline]
fn match_begin_s_m_ior_d(s: &str) -> Option<(std::ops::Range<usize>, i32, i32, i32, u8)> {
    let (g1, b'S', _, idx1) = next_cigar_op(s, 0)? else {
        return None;
    };
    let (tmid, b'M', _, idx2) = next_cigar_op(s, idx1)? else {
        return None;
    };
    let (g3, op, _, end) = next_cigar_op(s, idx2)?;
    if !matches!(op, b'I' | b'D') {
        return None;
    }
    Some((0..end, g1, tmid, g3, op))
}

#[inline]
fn match_ior_d_m_s_end(s: &str) -> Option<(std::ops::Range<usize>, i32, u8, i32, i32)> {
    let (s_start, g4, b'S', _) = prev_cigar_op(s, s.len())? else {
        return None;
    };
    let (m_start, tmid, b'M', _) = prev_cigar_op(s, s_start)? else {
        return None;
    };
    let (start, g1, op, _) = prev_cigar_op(s, m_start)?;
    if !matches!(op, b'I' | b'D') {
        return None;
    }
    Some((start..s.len(), g1, op, tmid, g4))
}

#[inline]
fn match_begin_digit_m_ior_d_m(s: &str) -> Option<(std::ops::Range<usize>, i32, i32, u8, i32)> {
    let (tmid, b'M', 1, idx1) = next_cigar_op(s, 0)? else {
        return None;
    };
    let (g2, op, _, idx2) = next_cigar_op(s, idx1)?;
    if !matches!(op, b'I' | b'D') {
        return None;
    }
    let (mlen, b'M', _, end) = next_cigar_op(s, idx2)? else {
        return None;
    };
    Some((0..end, tmid, g2, op, mlen))
}

#[inline]
fn match_ior_d_digit_m_end(s: &str) -> Option<(std::ops::Range<usize>, i32, u8, i32)> {
    let (m_start, tmid, b'M', 1) = prev_cigar_op(s, s.len())? else {
        return None;
    };
    let (start, g1, op, _) = prev_cigar_op(s, m_start)?;
    if !matches!(op, b'I' | b'D') {
        return None;
    }
    Some((start..s.len(), g1, op, tmid))
}

#[inline]
fn find_two_del_ins(s: &str) -> Option<(std::ops::Range<usize>, usize, [i32; 7])> {
    let line = s.split('\n').next().unwrap_or("");
    let mut idx = 0;
    while idx < line.len() {
        if let Some((g2, op1, _, i1)) = next_cigar_op(line, idx) {
            if op1 == b'M' {
                if let Some((g3, b'D', _, i2)) = next_cigar_op(line, i1) {
                    if let Some((g4, b'M', _, i3)) = next_cigar_op(line, i2) {
                        if let Some((g5, b'I', _, i4)) = next_cigar_op(line, i3) {
                            if let Some((g6, b'M', _, i5)) = next_cigar_op(line, i4) {
                                if let Some((g7, b'D', _, i6)) = next_cigar_op(line, i5) {
                                    if let Some((g8, b'M', _, end)) = next_cigar_op(line, i6) {
                                        return Some((idx..end, idx, [g2, g3, g4, g5, g6, g7, g8]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            idx = i1;
        } else {
            idx = step_past_current_run(line, idx);
        }
    }
    None
}

#[inline]
fn find_three_deletions(s: &str) -> Option<(std::ops::Range<usize>, usize, [i32; 7])> {
    let line = s.split('\n').next().unwrap_or("");
    let mut idx = 0;
    while idx < line.len() {
        if let Some((g2, op1, _, i1)) = next_cigar_op(line, idx) {
            if op1 == b'M' {
                if let Some((g3, b'D', _, i2)) = next_cigar_op(line, i1) {
                    if let Some((g4, b'M', _, i3)) = next_cigar_op(line, i2) {
                        if let Some((g5, b'D', _, i4)) = next_cigar_op(line, i3) {
                            if let Some((g6, b'M', _, i5)) = next_cigar_op(line, i4) {
                                if let Some((g7, b'D', _, i6)) = next_cigar_op(line, i5) {
                                    if let Some((g8, b'M', _, end)) = next_cigar_op(line, i6) {
                                        return Some((idx..end, idx, [g2, g3, g4, g5, g6, g7, g8]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            idx = i1;
        } else {
            idx = step_past_current_run(line, idx);
        }
    }
    None
}

struct ThreeIndelsMatch {
    range: std::ops::Range<usize>,
    ov5_len: usize,
    g2: i32,
    g3: i32,
    op1: u8,
    g5: i32,
    g6: i32,
    op2: u8,
    g8: i32,
    g9: i32,
    op3: u8,
    g11: i32,
}

struct CloseIndelMatch {
    range: std::ops::Range<usize>,
    first_len: i32,
    mid_len: i32,
    second_len: i32,
    second_op: u8,
    trailing_i: Option<i32>,
}

struct PrefixedCloseIndelMatch {
    prefix_ch: char,
    range_after_prefix: std::ops::Range<usize>,
    first_i_len: i32,
    mid_len: i32,
    second_len: i32,
    second_op: u8,
    trailing_i: Option<i32>,
}

#[inline]
fn find_three_indels(s: &str) -> Option<ThreeIndelsMatch> {
    let line = s.split('\n').next().unwrap_or("");
    let mut idx = 0;
    while idx < line.len() {
        if let Some((g2, op_m1, _, i1)) = next_cigar_op(line, idx) {
            if op_m1 == b'M' {
                if let Some((g3, op1, _, i2)) = next_cigar_op(line, i1) {
                    if matches!(op1, b'D' | b'I') {
                        if let Some((g5, b'M', _, i3)) = next_cigar_op(line, i2) {
                            if let Some((g6, op2, _, i4)) = next_cigar_op(line, i3) {
                                if matches!(op2, b'D' | b'I') {
                                    if let Some((g8, b'M', _, i5)) = next_cigar_op(line, i4) {
                                        if let Some((g9, op3, _, i6)) = next_cigar_op(line, i5) {
                                            if matches!(op3, b'D' | b'I') {
                                                if let Some((g11, b'M', _, end)) =
                                                    next_cigar_op(line, i6)
                                                {
                                                    return Some(ThreeIndelsMatch {
                                                        range: idx..end,
                                                        ov5_len: idx,
                                                        g2,
                                                        g3,
                                                        op1,
                                                        g5,
                                                        g6,
                                                        op2,
                                                        g8,
                                                        g9,
                                                        op3,
                                                        g11,
                                                    });
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            idx = i1;
        } else {
            idx = step_past_current_run(line, idx);
        }
    }
    None
}

#[inline]
fn find_dig_d_dig_m_dig_di_digi(s: &str) -> Option<CloseIndelMatch> {
    let mut idx = 0;
    while idx < s.len() {
        if let Some((g1, op1, _, i1)) = next_cigar_op(s, idx) {
            if op1 == b'D' {
                if let Some((g2, b'M', _, i2)) = next_cigar_op(s, i1) {
                    if let Some((g3, op4, _, i3)) = next_cigar_op(s, i2) {
                        if matches!(op4, b'D' | b'I') {
                            let (end, trailing_i) = match next_cigar_op(s, i3) {
                                Some((ti, b'I', _, i4)) => (i4, Some(ti)),
                                _ => (i3, None),
                            };
                            return Some(CloseIndelMatch {
                                range: idx..end,
                                first_len: g1,
                                mid_len: g2,
                                second_len: g3,
                                second_op: op4,
                                trailing_i,
                            });
                        }
                    }
                }
            }
            idx = i1;
        } else {
            idx = step_past_current_run(s, idx);
        }
    }
    None
}

#[inline]
fn find_notdig_dig_i_dig_m_dig_di_digi(s: &str) -> Option<PrefixedCloseIndelMatch> {
    let mut idx = 0;
    while idx < s.len() {
        if let Some((g2, op1, _, i1)) = next_cigar_op(s, idx) {
            if op1 == b'I' && idx > 0 {
                if let Some(prefix_ch) = s.get(..idx).and_then(|h| h.chars().next_back()) {
                    if !is_regex_decimal_digit(prefix_ch) {
                        if let Some((g3, b'M', _, i2)) = next_cigar_op(s, i1) {
                            if let Some((g4, op5, _, i3)) = next_cigar_op(s, i2) {
                                if matches!(op5, b'D' | b'I') {
                                    let (end, trailing_i) = match next_cigar_op(s, i3) {
                                        Some((ti, b'I', _, i4)) => (i4, Some(ti)),
                                        _ => (i3, None),
                                    };
                                    return Some(PrefixedCloseIndelMatch {
                                        prefix_ch,
                                        range_after_prefix: idx..end,
                                        first_i_len: g2,
                                        mid_len: g3,
                                        second_len: g4,
                                        second_op: op5,
                                        trailing_i,
                                    });
                                }
                            }
                        }
                    }
                }
            }
            idx = i1;
        } else {
            idx = step_past_current_run(s, idx);
        }
    }
    None
}

#[inline]
fn find_adjacent_same_op(s: &str, target: u8) -> Option<(std::ops::Range<usize>, i32, i32)> {
    let mut idx = 0;
    while idx < s.len() {
        if let Some((g1, op1, _, i1)) = next_cigar_op(s, idx) {
            if op1 == target {
                if let Some((g2, op2, _, end)) = next_cigar_op(s, i1) {
                    if op2 == target {
                        return Some((idx..end, g1, g2));
                    }
                }
            }
            idx = i1;
        } else {
            idx = step_past_current_run(s, idx);
        }
    }
    None
}

/// Single-pass replacement for `(global_find_sum(&ALIGNED_LENGTH_MND, s), global_find_sum(&SOFT_CLIPPED, s))`,
/// matching `(\d+)[MND]` and `(\d+)[MIS]` byte-for-byte without regex or heap allocation.
#[inline]
fn sum_mnd_and_mis_prefix(s: &str) -> (i32, i32) {
    let mut mnd: i32 = 0;
    let mut mis: i32 = 0;
    let mut idx = 0;
    while idx < s.len() {
        if let Some((val, op, _, next_idx)) = next_cigar_op(s, idx) {
            let is_mnd = matches!(op, b'M' | b'N' | b'D');
            let is_mis = matches!(op, b'M' | b'I' | b'S');
            if is_mnd {
                mnd = mnd.saturating_add(val);
            }
            if is_mis {
                mis = mis.saturating_add(val);
            }
            idx = next_idx;
        } else {
            idx = step_past_current_run(s, idx);
        }
    }
    (mnd, mis)
}

#[inline]
fn cigar_flags(s: &str) -> (bool, bool, bool) {
    let mut has_d = false;
    let mut has_i = false;
    let mut has_s = false;
    for &b in s.as_bytes() {
        match b {
            b'D' => has_d = true,
            b'I' => has_i = true,
            b'S' => has_s = true,
            _ => {}
        }
    }
    (has_d, has_i, has_s)
}

/// Matches `ANY_NUMBER_M_NUMBER_S_END` (`^(.*?)(\d+)M(\d+)S$`) without regex or heap allocation.
#[inline]
fn match_any_number_m_number_s_end(s: &str) -> Option<(&str, i32, i32)> {
    let (s_start, d2, b'S', _) = prev_cigar_op(s, s.len())? else {
        return None;
    };
    let (d1_start, d1, b'M', _) = prev_cigar_op(s, s_start)? else {
        return None;
    };
    let prefix = s.get(..d1_start)?;
    if prefix.as_bytes().contains(&b'\n') {
        return None;
    }
    Some((prefix, d1, d2))
}

/// Matches `BEGIN_ANY_DIG_M_END` (`^(.*?)(\d+)M$`) without regex or heap allocation.
#[inline]
fn match_begin_any_dig_m_end(s: &str) -> Option<(&str, i32)> {
    let (d1_start, d1, b'M', _) = prev_cigar_op(s, s.len())? else {
        return None;
    };
    let prefix = s.get(..d1_start)?;
    if prefix.as_bytes().contains(&b'\n') {
        return None;
    }
    Some((prefix, d1))
}

/// Matches `DIG_S_DIG_M` (`^(\d+)S(\d+)M`) without regex or heap allocation.
#[inline]
fn match_dig_s_dig_m(s: &str) -> Option<(i32, i32)> {
    let (d1, b'S', _, idx1) = next_cigar_op(s, 0)? else {
        return None;
    };
    let (d2, b'M', _, _) = next_cigar_op(s, idx1)? else {
        return None;
    };
    Some((d1, d2))
}

/// Matches `BEGIN_DIG_M` (`^(\d+)M`) without regex or heap allocation.
#[inline]
fn match_begin_dig_m(s: &str) -> Option<i32> {
    let (d1, b'M', _, _) = next_cigar_op(s, 0)? else {
        return None;
    };
    Some(d1)
}

/// Scans `(\-?\d+)[MIDNSHPX=]` across `s` and validates that all matched lengths are `> 0`.
#[inline]
fn validate_cigar_element_lengths(s: &str) -> Result<(), crate::types::VardictError> {
    let bytes = s.as_bytes();
    let mut idx = 0;
    while idx < bytes.len() {
        let start = idx;
        if bytes.get(idx) == Some(&b'-') {
            idx += 1;
        }
        let (d_end, digit_chars) = scan_digits_fwd(s, idx);
        if digit_chars > 0 {
            if let Some(&op) = bytes.get(d_end) {
                if matches!(
                    op,
                    b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'='
                ) {
                    let num_str = s.get(start..d_end).unwrap_or("");
                    let len: i32 = num_str.parse().unwrap_or(0);
                    if len <= 0 {
                        return Err(crate::types::VardictError::JavaThrows(format!(
                            "invalid CIGAR length {} in {}",
                            len, s
                        )));
                    }
                    idx = d_end + 1;
                    continue;
                }
            }
            idx = d_end;
        } else if start == idx {
            idx += if bytes.get(idx).is_some_and(u8::is_ascii) {
                1
            } else {
                s.get(idx..)
                    .and_then(|t| t.chars().next())
                    .map_or(1, char::len_utf8)
            };
        }
    }
    Ok(())
}

/// Sum all matches of a regex group 1 parsed as integers.
#[cfg(test)]
fn global_find_sum(re: &Regex, s: &str) -> i32 {
    let mut sum: i32 = 0;
    for cap in re.captures_iter(s) {
        if let Some(m) = cap.get(1) {
            sum = sum.saturating_add(to_int(m.as_str()));
        }
    }
    sum
}

// ---------------------------------------------------------------------------
// Reference helper functions matching VariationUtils
// ---------------------------------------------------------------------------

#[inline]
fn is_has_and_equals_char(ch1: u8, reference: &[u8], index: i32, ref_start: i32) -> bool {
    usize::try_from(i64::from(index) - i64::from(ref_start))
        .ok()
        .and_then(|local_idx| reference.get(local_idx))
        .is_some_and(|&b| b == ch1)
}

#[inline]
fn is_has_and_equals_str(
    reference: &[u8],
    index1: i64,
    ref_start: i32,
    str_bytes: &[u8],
    index2: i32,
) -> bool {
    let Some(&ref_b) = usize::try_from(index1 - i64::from(ref_start))
        .ok()
        .and_then(|local_idx| reference.get(local_idx))
    else {
        return false;
    };
    let Some(&str_b) = usize::try_from(index2).ok().and_then(|i| str_bytes.get(i)) else {
        return false;
    };
    ref_b == str_b
}

#[inline]
fn is_has_and_not_equals_char(ch1: u8, reference: &[u8], index: i32, ref_start: i32) -> bool {
    usize::try_from(i64::from(index) - i64::from(ref_start))
        .ok()
        .and_then(|local_idx| reference.get(local_idx))
        .is_some_and(|&b| b != ch1)
}

/// Scans up to `mch` bases from the 5' start of a match block (at `position` on reference
/// and `query_offset` in `query_seq`), resetting on any mismatch and stopping after 3
/// consecutive matches (`rmch >= 3`) or reference out-of-bounds.
/// Ports the 5' mismatch scan in `CigarModifier.java:683-702` (`combineBeginDigM`) and
/// `CigarModifier.java:656-672` (`combineDigSDigM`).
#[inline]
pub(crate) fn count_start_mismatch_clip(
    position: i32,
    mch: i32,
    query_offset: i32,
    query_seq: &[u8],
    reference: &[u8],
    ref_start: i32,
) -> i32 {
    let mut rn = 0i32;
    let mut rrn = 0i32;
    let mut rmch = 0i32;
    while rrn < mch && rn < mch {
        let ref_pos = position.saturating_add(rrn);
        let local_ref = ref_pos.saturating_sub(ref_start);
        let Some(&ref_base) = usize::try_from(local_ref)
            .ok()
            .and_then(|i| reference.get(i))
        else {
            break;
        };
        let q_idx = query_offset.saturating_add(rrn);
        if let Some(&q_base) = usize::try_from(q_idx).ok().and_then(|i| query_seq.get(i)) {
            if ref_base != q_base {
                rn = rrn.saturating_add(1);
                rmch = 0;
            } else {
                rmch = rmch.saturating_add(1);
            }
        }
        rrn = rrn.saturating_add(1);
        if rmch >= 3 {
            break;
        }
    }
    rn
}

/// Scans up to `mch` bases backwards from the 3' end of a match block (ending at `refoff`
/// on reference and `rdoff` in `query_seq`), resetting on any mismatch and stopping after
/// 3 consecutive matches (`rmch >= 3`) or reference out-of-bounds.
/// Ports the 3' mismatch scan in `CigarModifier.java:704-728` (`captureMisSoftly3Mismatches`)
/// and `CigarModifier.java:765-785` (`captureMisSoftlyMS`).
#[inline]
pub(crate) fn count_end_mismatch_clip(
    refoff: i32,
    rdoff: i32,
    mch: i32,
    query_seq: &[u8],
    reference: &[u8],
    ref_start: i32,
) -> i32 {
    let mut rn = 0i32;
    let mut rrn = 0i32;
    let mut rmch = 0i32;
    while rrn < mch && rn < mch {
        let ref_pos = refoff.saturating_sub(rrn).saturating_sub(1);
        let local_ref = ref_pos.saturating_sub(ref_start);
        let Some(&ref_base) = usize::try_from(local_ref)
            .ok()
            .and_then(|i| reference.get(i))
        else {
            break;
        };
        if rrn < rdoff {
            let q_idx = rdoff.saturating_sub(rrn).saturating_sub(1);
            if let Some(&q_base) = usize::try_from(q_idx).ok().and_then(|i| query_seq.get(i)) {
                if ref_base != q_base {
                    rn = rrn.saturating_add(1);
                    rmch = 0;
                } else {
                    rmch = rmch.saturating_add(1);
                }
            }
        }
        rrn = rrn.saturating_add(1);
        if rmch >= 3 {
            break;
        }
    }
    rn
}

/// Fast-path for reads whose initial CIGAR is a single `Match(orig_mch)` (`orig_mch > 0`).
/// Ports `CigarModifier.java:61-368` for a single-`M` CIGAR where `has_d == false`,
/// `has_i == false`, and `has_s == false`, so only `capture_mis_softly_3mismatches` (3' end)
/// and `combine_begin_dig_m` (5' end) can fire.
pub fn modify_single_match_in_place(
    position: &mut i32,
    orig_mch: i32,
    cigar_elements: &mut Vec<(u32, noodles::sam::alignment::record::cigar::op::Kind)>,
    query_seq: &[u8],
    reference: &[u8],
    ref_start: i32,
) -> Result<(), crate::types::VardictError> {
    use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

    let pos = *position;

    // Step A: 3' endpoint mismatch -> soft-clip (`capture_mis_softly_3mismatches("", orig_mch)`)
    let refoff = pos.saturating_add(orig_mch);
    let rn3 = count_end_mismatch_clip(refoff, orig_mch, orig_mch, query_seq, reference, ref_start);
    let (mch_after_3, clip3) = if (1..=3).contains(&rn3) {
        (orig_mch.saturating_sub(rn3), rn3 as u32)
    } else {
        (orig_mch, 0u32)
    };

    // Step B: 5' endpoint mismatch -> soft-clip (`combine_begin_dig_m(mch_after_3)`)
    let rn5 = count_start_mismatch_clip(pos, mch_after_3, 0, query_seq, reference, ref_start);
    let (final_mch, clip5) = if (1..=3).contains(&rn5) {
        (mch_after_3.saturating_sub(rn5), rn5 as u32)
    } else {
        (mch_after_3, 0u32)
    };

    // Validate positive CIGAR element lengths (`validate_cigar_element_lengths`)
    if final_mch <= 0 {
        let bad_cigar = match (clip5 > 0, clip3 > 0) {
            (true, true) => format!("{clip5}S{final_mch}M{clip3}S"),
            (true, false) => format!("{clip5}S{final_mch}M"),
            (false, true) => format!("{final_mch}M{clip3}S"),
            (false, false) => format!("{final_mch}M"),
        };
        return Err(crate::types::VardictError::JavaThrows(format!(
            "invalid CIGAR length {final_mch} in {bad_cigar}"
        )));
    }

    if clip5 == 0 && clip3 == 0 {
        return Ok(());
    }

    *position = pos.saturating_add(clip5 as i32);
    cigar_elements.clear();
    if clip5 > 0 {
        cigar_elements.push((clip5, CigarKind::SoftClip));
    }
    cigar_elements.push((final_mch as u32, CigarKind::Match));
    if clip3 > 0 {
        cigar_elements.push((clip3, CigarKind::SoftClip));
    }
    Ok(())
}

/// Outcome of [`modify_pure_sm_in_place`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PureSmOutcome {
    /// Read CIGAR and start position were updated in place with `query_seq` / `query_qual` untouched.
    Handled,
    /// Preconditions not met; caller must fall back to [`CigarModifier::modify_cigar_preserving_buffers`].
    Fallback,
}

/// Shared numeric core for `CigarModifier.java:398-465` (`captureMisSoftlyMS`).
/// Note that `refoff` and `rdoff` intentionally do NOT advance across Phase 1 rescue into
/// Phase 2 rescue or terminal mismatch clipping, matching `CigarModifier.java:422-462`.
#[inline]
pub(crate) fn step_capture_mis_softly_ms(
    refoff: i32,
    rdoff: i32,
    mut mch: i32,
    mut soft: i32,
    query_seq: &[u8],
    query_qual: &[u8],
    reference: &[u8],
    ref_start: i32,
) -> (i32, i32, bool) {
    let mut rn = 0i32;
    let mut rn_set = ByteSet::new();
    let mut changed = false;

    while rn < soft
        && is_has_and_equals_str(
            reference,
            i64::from(refoff) + i64::from(rn),
            ref_start,
            query_seq,
            rdoff.saturating_add(rn),
        )
        && char_at(query_qual, rdoff.saturating_add(rn) as isize)
            .unwrap_or(0)
            .saturating_sub(33)
            > LOWQUAL
    {
        rn += 1;
    }

    if rn > 0 {
        mch = mch.saturating_add(rn);
        soft -= rn;
        changed = true;
        rn = 0;
    }

    if soft > 0 {
        while rn + 1 < soft
            && is_has_and_equals_str(
                reference,
                i64::from(refoff) + i64::from(rn) + 1,
                ref_start,
                query_seq,
                rdoff.saturating_add(rn).saturating_add(1),
            )
            && char_at(
                query_qual,
                rdoff.saturating_add(rn).saturating_add(1) as isize,
            )
            .unwrap_or(0)
            .saturating_sub(33)
                > LOWQUAL
        {
            rn += 1;
            let ref_pos = i64::from(refoff) + i64::from(rn) + 1;
            if let Some(&b) = usize::try_from(ref_pos - i64::from(ref_start))
                .ok()
                .and_then(|idx| reference.get(idx))
            {
                rn_set.insert(b);
            }
        }
        let rn_nt = rn_set.len();
        if rn > 4 && rn_nt > 1 {
            mch = mch.saturating_add(rn + 1);
            soft -= rn + 1;
            changed = true;
        }
        if rn == 0 {
            rn = count_end_mismatch_clip(refoff, rdoff, mch, query_seq, reference, ref_start);
            if rn > 0 && rn < mch {
                soft = soft.saturating_add(rn);
                mch -= rn;
                changed = true;
            }
        }
    }

    (mch, soft, changed)
}

/// Shared numeric core for `CigarModifier.java:285-352` (`combineDigSDigM`).
/// Returns `(final_pos, cigar_soft, cigar_mch, cigar_changed)`, faithfully reproducing
/// `DIG_S_DIG_M` no-op behavior when Phase 2 rescues `soft` down to `0` (`cigar_soft == 0`)
/// and Phase 3 (`rn == 0`) still advances `position` (`CigarModifier.java:318-349`).
#[inline]
pub(crate) fn step_combine_dig_s_dig_m(
    mut pos: i32,
    mut soft: i32,
    mut mch: i32,
    query_seq: &[u8],
    query_qual: &[u8],
    reference: &[u8],
    ref_start: i32,
) -> (i32, i32, i32, bool) {
    let mut rn = 0i32;
    let mut rn_set = ByteSet::new();
    let mut cigar_soft = soft;
    let mut cigar_mch = mch;
    let mut cigar_changed = false;

    while rn < soft
        && is_has_and_equals_str(
            reference,
            i64::from(pos) - i64::from(rn) - 1,
            ref_start,
            query_seq,
            soft - rn - 1,
        )
        && char_at(query_qual, (soft - rn - 1) as isize)
            .unwrap_or(0)
            .saturating_sub(33)
            > LOWQUAL
    {
        rn += 1;
    }

    if rn > 0 {
        mch = mch.saturating_add(rn);
        soft -= rn;
        cigar_soft = soft;
        cigar_mch = mch;
        cigar_changed = true;
        pos = pos.saturating_sub(rn);
        rn = 0;
    }

    if soft > 0 {
        while rn + 1 < soft
            && is_has_and_equals_str(
                reference,
                i64::from(pos) - i64::from(rn) - 2,
                ref_start,
                query_seq,
                soft - rn - 2,
            )
            && char_at(query_qual, (soft - rn - 2) as isize)
                .unwrap_or(0)
                .saturating_sub(33)
                > LOWQUAL
        {
            rn += 1;
            let ref_pos = i64::from(pos) - i64::from(rn) - 2;
            if let Some(&b) = usize::try_from(ref_pos - i64::from(ref_start))
                .ok()
                .and_then(|idx| reference.get(idx))
            {
                rn_set.insert(b);
            }
        }
        let rn_nt = rn_set.len();
        if (rn > 4 && rn_nt > 1)
            || is_has_and_equals_str(
                reference,
                i64::from(pos) - 1,
                ref_start,
                query_seq,
                soft - 1,
            )
        {
            mch = mch.saturating_add(rn + 1);
            soft -= rn + 1;
            cigar_soft = soft;
            cigar_mch = mch;
            cigar_changed = true;
            pos = pos.saturating_sub(rn + 1);
        }
        if rn == 0 {
            rn = count_start_mismatch_clip(pos, mch, soft, query_seq, reference, ref_start);
            if rn > 0 && rn < mch {
                soft = soft.saturating_add(rn);
                mch -= rn;
                if cigar_soft > 0 {
                    cigar_soft = soft;
                    cigar_mch = mch;
                    cigar_changed = true;
                }
                pos = pos.saturating_add(rn);
            }
        }
    }

    (pos, cigar_soft, cigar_mch, cigar_changed)
}

#[inline]
fn unique_seed_map_pos(seed: &SeedMap, rc_seed: &[u8]) -> Option<u32> {
    seed.get_unique(rc_seed)
}

#[inline]
fn is_chimeric_seed_hit(seed_pos: Option<u32>, position: i32, max_read_length: i32) -> bool {
    seed_pos.is_some_and(|ref_pos_u32| {
        let ref_pos = ref_pos_u32 as i32;
        (i64::from(position) - i64::from(ref_pos)).abs() < 2 * i64::from(max_read_length)
    })
}

#[inline]
fn replace_cigar_range(
    cigar_str: &mut String,
    range: std::ops::Range<usize>,
    args: std::fmt::Arguments<'_>,
) {
    if range.start <= range.end
        && cigar_str.is_char_boundary(range.start)
        && cigar_str.is_char_boundary(range.end)
    {
        let mut buf = crate::types::StackBuf::<64>::new();
        if buf.write_fmt(args).is_ok() {
            cigar_str.replace_range(range, buf.as_str());
        } else {
            let fallback = std::fmt::format(args);
            cigar_str.replace_range(range, &fallback);
        }
    }
}

/// Zero-copy in-place CIGAR modifier for 2-op (`[S, M]`, `[M, S]`) and 3-op (`[S, M, S]`) reads
/// whose `query_seq` and `query_qual` are guaranteed untouched by `CigarModifier.java:61-243`.
#[expect(
    clippy::too_many_arguments,
    reason = "in-place pure-SM CIGAR modifier borrows disjoint CigarParser fields without allocation"
)]
pub(crate) fn modify_pure_sm_in_place(
    position: &mut i32,
    cigar_elements: &mut Vec<(u32, noodles::sam::alignment::record::cigar::op::Kind)>,
    query_seq: &[u8],
    query_qual: &[u8],
    reference: &[u8],
    ref_start: i32,
    seed: &SeedMap,
    max_read_length: i32,
    config: &Config,
) -> Result<PureSmOutcome, crate::types::VardictError> {
    use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

    let (s5, mch, s3, has_s5, has_s3) = match *cigar_elements.as_slice() {
        [(s5_u32, CigarKind::SoftClip), (m_u32, CigarKind::Match)] => {
            let (Ok(s5), Ok(mch)) = (i32::try_from(s5_u32), i32::try_from(m_u32)) else {
                return Ok(PureSmOutcome::Fallback);
            };
            if s5 <= 0 || mch <= 0 {
                return Ok(PureSmOutcome::Fallback);
            }
            (s5, mch, 0i32, true, false)
        }
        [(m_u32, CigarKind::Match), (s3_u32, CigarKind::SoftClip)] => {
            let (Ok(mch), Ok(s3)) = (i32::try_from(m_u32), i32::try_from(s3_u32)) else {
                return Ok(PureSmOutcome::Fallback);
            };
            if mch <= 0 || s3 <= 0 {
                return Ok(PureSmOutcome::Fallback);
            }
            (0i32, mch, s3, false, true)
        }
        [(s5_u32, CigarKind::SoftClip), (m_u32, CigarKind::Match), (s3_u32, CigarKind::SoftClip)] =>
        {
            let (Ok(s5), Ok(mch), Ok(s3)) = (
                i32::try_from(s5_u32),
                i32::try_from(m_u32),
                i32::try_from(s3_u32),
            ) else {
                return Ok(PureSmOutcome::Fallback);
            };
            if s5 <= 0 || mch <= 0 || s3 <= 0 {
                return Ok(PureSmOutcome::Fallback);
            }
            (s5, mch, s3, true, true)
        }
        _ => return Ok(PureSmOutcome::Fallback),
    };

    // 1. Sequence/quality length invariants (`!query_seq.is_empty()` and op sum == len).
    if query_seq.is_empty() || query_seq.len() != query_qual.len() {
        return Ok(PureSmOutcome::Fallback);
    }
    let Some(expected_qlen) = (s5 as usize)
        .checked_add(mch as usize)
        .and_then(|v| v.checked_add(s3 as usize))
    else {
        return Ok(PureSmOutcome::Fallback);
    };
    if expected_qlen != query_seq.len() {
        return Ok(PureSmOutcome::Fallback);
    }

    // 2. Contig/loaded-reference bounds (no overhang past contig boundaries).
    let pos0 = *position;
    if pos0 <= 0 || reference.is_empty() {
        return Ok(PureSmOutcome::Fallback);
    }
    let pos0_i64 = i64::from(pos0);
    let m_end_i64 = pos0_i64 + i64::from(mch) - 1;
    let ref_start_i64 = i64::from(ref_start);
    let ref_end_i64 = ref_start_i64 + (reference.len() as i64) - 1;
    if pos0_i64 < ref_start_i64 || m_end_i64 > ref_end_i64 {
        return Ok(PureSmOutcome::Fallback);
    }

    // 3. No adaptor trimming.
    if !config.adaptor.is_empty() {
        return Ok(PureSmOutcome::Fallback);
    }

    // 4. Mirror Patterns.java:156,159 (`^(\d\d+)S` / `(\d\d+)S$`) chimeric SEED_2 check (`AR-01`):
    //    5' `>= 10S` enters Java's `if` branch (masking `else if` 3' check even for 10S/11S or 5' miss).
    if !config.chimeric {
        if s5 >= 10 {
            if (s5 as usize) >= SEED_2 {
                let s5_usize = s5 as usize;
                let seed_pos = query_seq
                    .get((s5_usize - SEED_2)..s5_usize)
                    .and_then(reverse_complement_bytes::<SEED_2>)
                    .and_then(|rc_seed| unique_seed_map_pos(seed, &rc_seed));
                if is_chimeric_seed_hit(seed_pos, pos0, max_read_length) {
                    return Ok(PureSmOutcome::Fallback);
                }
            }
        } else if (s3 as usize) >= SEED_2 {
            let seed_pos = query_seq
                .len()
                .checked_sub(s3 as usize)
                .and_then(|idx| query_seq.get(idx..))
                .and_then(|s| s.get(..SEED_2))
                .and_then(reverse_complement_bytes::<SEED_2>)
                .and_then(|rc_seed| unique_seed_map_pos(seed, &rc_seed));
            if is_chimeric_seed_hit(seed_pos, pos0, max_read_length) {
                return Ok(PureSmOutcome::Fallback);
            }
        }
    }

    let mut pos = pos0;
    let mut cur_s5 = s5;
    let mut cur_mch = mch;
    let mut cur_s3 = s3;

    // Step 1 (3' end FIRST — CigarModifier.java:230-235):
    let refoff = pos.saturating_add(cur_mch);
    let rdoff = cur_s5.saturating_add(cur_mch);
    if has_s3 {
        let (m3, s3_out, _) = step_capture_mis_softly_ms(
            refoff, rdoff, cur_mch, cur_s3, query_seq, query_qual, reference, ref_start,
        );
        cur_mch = m3;
        cur_s3 = s3_out;
    } else {
        let rn3 = count_end_mismatch_clip(refoff, rdoff, cur_mch, query_seq, reference, ref_start);
        if (1..=3).contains(&rn3) {
            cur_mch -= rn3;
            cur_s3 = rn3;
        }
    }

    // Step 2 (5' start SECOND — CigarModifier.java:237-242):
    if has_s5 {
        let (new_pos, s5_out, m5_out, _) = step_combine_dig_s_dig_m(
            pos, cur_s5, cur_mch, query_seq, query_qual, reference, ref_start,
        );
        pos = new_pos;
        cur_s5 = s5_out;
        cur_mch = m5_out;
    } else {
        let rn5 = count_start_mismatch_clip(pos, cur_mch, 0, query_seq, reference, ref_start);
        if (1..=3).contains(&rn5) {
            cur_mch = cur_mch.saturating_sub(rn5);
            cur_s5 = rn5;
            pos = pos.saturating_add(rn5);
        }
    }

    // Validate strictly positive element lengths (`validate_cigar_element_lengths`).
    if cur_mch <= 0 {
        let bad_cigar = match (cur_s5 > 0, cur_s3 > 0) {
            (true, true) => format!("{cur_s5}S{cur_mch}M{cur_s3}S"),
            (true, false) => format!("{cur_s5}S{cur_mch}M"),
            (false, true) => format!("{cur_mch}M{cur_s3}S"),
            (false, false) => format!("{cur_mch}M"),
        };
        return Err(crate::types::VardictError::JavaThrows(format!(
            "invalid CIGAR length {cur_mch} in {bad_cigar}"
        )));
    }

    *position = pos;
    cigar_elements.clear();
    if cur_s5 > 0 {
        cigar_elements.push((cur_s5 as u32, CigarKind::SoftClip));
    }
    cigar_elements.push((cur_mch as u32, CigarKind::Match));
    if cur_s3 > 0 {
        cigar_elements.push((cur_s3 as u32, CigarKind::SoftClip));
    }
    Ok(PureSmOutcome::Handled)
}

// Helper to reverse complement sequence string
fn rc_seq(seq: &[u8]) -> Vec<u8> {
    let mut v = seq.to_vec();
    reverse_complement(&mut v);
    v
}

enum SeedBackend<'a> {
    Map(&'a SeedMap),
    Closure(&'a mut dyn FnMut(&[u8]) -> Option<u32>),
}

/// The main CIGAR Modifier engine porting CigarModifier.java.
pub struct CigarModifier<'a> {
    position: i32,
    cigar_str: String,
    query_seq: Vec<u8>,
    query_qual: Vec<u8>,
    reference: &'a [u8],
    ref_start: i32,
    seed: SeedBackend<'a>,
    indel: i32,
    max_read_length: i32,
    config: &'a Config,
}

impl<'a> CigarModifier<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        position: i32,
        cigar_str: impl Into<String>,
        query_seq: impl Into<Vec<u8>>,
        query_qual: impl Into<Vec<u8>>,
        reference: &'a [u8],
        ref_start: i32,
        seed: &'a SeedMap,
        indel: i32,
        max_read_length: i32,
        config: &'a Config,
    ) -> Self {
        CigarModifier {
            position,
            cigar_str: cigar_str.into(),
            query_seq: query_seq.into(),
            query_qual: query_qual.into(),
            reference,
            ref_start,
            seed: SeedBackend::Map(seed),
            indel,
            max_read_length,
            config,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn with_seed_fn(
        position: i32,
        cigar_str: impl Into<String>,
        query_seq: impl Into<Vec<u8>>,
        query_qual: impl Into<Vec<u8>>,
        reference: &'a [u8],
        ref_start: i32,
        seed_fn: &'a mut dyn FnMut(&[u8]) -> Option<u32>,
        indel: i32,
        max_read_length: i32,
        config: &'a Config,
    ) -> Self {
        CigarModifier {
            position,
            cigar_str: cigar_str.into(),
            query_seq: query_seq.into(),
            query_qual: query_qual.into(),
            reference,
            ref_start,
            seed: SeedBackend::Closure(seed_fn),
            indel,
            max_read_length,
            config,
        }
    }

    #[inline]
    fn unique_seed_pos(&mut self, rc_seed: &[u8]) -> Option<u32> {
        match &mut self.seed {
            SeedBackend::Map(map) => unique_seed_map_pos(map, rc_seed),
            SeedBackend::Closure(f) => f(rc_seed),
        }
    }

    /// Reclaims the internal CIGAR, sequence, and quality buffers (e.g. when
    /// [`Self::modify_cigar_preserving_buffers`] returns `Err` so `CigarParser`
    /// can restore its scratch allocations).
    pub fn into_buffers(self) -> (String, Vec<u8>, Vec<u8>) {
        (self.cigar_str, self.query_seq, self.query_qual)
    }

    /// Perform all CIGAR modifications. No panic allowed!
    /// Takes `self` by value: the modifier is single-use per read, and the
    /// caller reclaims the buffers from the returned [`ModifiedCigar`].
    pub fn modify_cigar(mut self) -> Result<ModifiedCigar, crate::types::VardictError> {
        self.modify_cigar_preserving_buffers()
    }

    /// Perform all CIGAR modifications in place, moving `cigar_str`, `query_seq`,
    /// and `query_qual` into [`ModifiedCigar`] on `Ok` and retaining them inside
    /// `self` on `Err` for recovery via [`Self::into_buffers`].
    pub fn modify_cigar_preserving_buffers(
        &mut self,
    ) -> Result<ModifiedCigar, crate::types::VardictError> {
        //flag is set to true if CIGAR string is modified and should be looked at again
        let flag = true;
        let this = self;
        let (mut has_d, mut has_i, mut has_s) = cigar_flags(&this.cigar_str);

        if has_d {
            // 1. Strip leading deletion
            if let Some((dlen, b'D', _, end)) = next_cigar_op(&this.cigar_str, 0) {
                this.position = this.position.saturating_add(dlen);
                this.cigar_str.replace_range(0..end, "");
                (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
            }
            // Strip trailing deletion
            if has_d {
                if let Some((start, _, b'D', _)) =
                    prev_cigar_op(&this.cigar_str, this.cigar_str.len())
                {
                    let end = this.cigar_str.len();
                    this.cigar_str.replace_range(start..end, "");
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                }
            }
        }

        if has_i {
            // 2. Replace leading insertion with soft clip
            if let Some((ilen, b'I', _, end)) = next_cigar_op(&this.cigar_str, 0) {
                replace_cigar_range(&mut this.cigar_str, 0..end, format_args!("{ilen}S"));
                (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
            }
            // Replace trailing insertion with soft clip
            if has_i {
                if let Some((start, ilen, b'I', _)) =
                    prev_cigar_op(&this.cigar_str, this.cigar_str.len())
                {
                    let end = this.cigar_str.len();
                    replace_cigar_range(&mut this.cigar_str, start..end, format_args!("{ilen}S"));
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                }
            }
        }

        // 3. Chimeric soft-clip detection
        if has_s {
            if let Some(cigar_element) = leading_2digit_s_len(&this.cigar_str) {
                if !this.config.chimeric && cigar_element >= SEED_2 {
                    let rc_seed_opt = subseq(&this.query_seq, 0, cigar_element)
                        .and_then(|sseq| {
                            sseq.len()
                                .checked_sub(SEED_2)
                                .and_then(|idx| sseq.get(idx..))
                        })
                        .and_then(reverse_complement_bytes::<SEED_2>);
                    if let Some(rc_seed) = rc_seed_opt {
                        let seed_pos = this.unique_seed_pos(&rc_seed);
                        if is_chimeric_seed_hit(seed_pos, this.position, this.max_read_length) {
                            let verbose_seq = if this.config.verbose {
                                subseq(&this.query_seq, 0, cigar_element).map(rc_seq)
                            } else {
                                None
                            };
                            if let Some((_, b'S', _, end)) = next_cigar_op(&this.cigar_str, 0) {
                                this.cigar_str.replace_range(0..end, "");
                            }
                            (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                            if this.query_seq.len() >= cigar_element {
                                this.query_seq.drain(..cigar_element);
                                if this.query_qual.len() >= cigar_element {
                                    this.query_qual.drain(..cigar_element);
                                } else {
                                    this.query_qual.clear();
                                }
                            }
                            if let Some(sequence) = verbose_seq {
                                crate::verbose_eprintln!(
                                    "{} at 5' is a chimeric at {} by SEED {}",
                                    crate::util::SanitizedBytes(&sequence),
                                    this.position,
                                    SEED_2
                                );
                            }
                        }
                    }
                }
            } else if let Some(cigar_element) = trailing_2digit_s_len(&this.cigar_str) {
                if !this.config.chimeric && cigar_element >= SEED_2 {
                    let seq_len = this.query_seq.len();
                    if seq_len >= cigar_element {
                        let start_idx = seq_len - cigar_element;
                        let rc_seed_opt = subseq(&this.query_seq, start_idx, cigar_element)
                            .and_then(|sseq| sseq.get(..SEED_2))
                            .and_then(reverse_complement_bytes::<SEED_2>);
                        if let Some(rc_seed) = rc_seed_opt {
                            let seed_pos = this.unique_seed_pos(&rc_seed);
                            if is_chimeric_seed_hit(seed_pos, this.position, this.max_read_length) {
                                let verbose_seq = if this.config.verbose {
                                    subseq(&this.query_seq, start_idx, cigar_element).map(rc_seq)
                                } else {
                                    None
                                };
                                if let Some((start, _, b'S', dc)) =
                                    prev_cigar_op(&this.cigar_str, this.cigar_str.len())
                                {
                                    if dc >= 2 {
                                        let end = this.cigar_str.len();
                                        this.cigar_str.replace_range(start..end, "");
                                    }
                                }
                                (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                                this.query_seq.truncate(start_idx);
                                if this.query_qual.len() >= start_idx {
                                    this.query_qual.truncate(start_idx);
                                } else {
                                    this.query_qual.clear();
                                }
                                if let Some(sequence) = verbose_seq {
                                    crate::verbose_eprintln!(
                                        "{} at 3' is a chimeric at {} by SEED {}",
                                        crate::util::SanitizedBytes(&sequence),
                                        this.position,
                                        SEED_2
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        // 4. Main indel simplification loop
        let mut flag_loop = flag;
        let mut iter_count = 0;
        let max_iter = this
            .cigar_str
            .bytes()
            .filter(u8::is_ascii_alphabetic)
            .count()
            + 1;
        while flag_loop && this.indel > 0 && (has_d || has_i) {
            iter_count += 1;
            if iter_count > max_iter {
                log::warn!(
                    "Infinite loop check triggered in modify_cigar! CIGAR: {}",
                    crate::util::SanitizedBytes(this.cigar_str.as_bytes())
                );
                break;
            }
            flag_loop = false;

            // S followed by I/D
            if has_s && (has_d || has_i) {
                if let Some((range, g1, g2, op)) = match_begin_s_ior_d(&this.cigar_str) {
                    let slen = g1.saturating_add(if op == b'I' { g2 } else { 0 });
                    this.position = this
                        .position
                        .saturating_add(if op == b'D' { g2 } else { 0 });
                    replace_cigar_range(&mut this.cigar_str, range, format_args!("{slen}S"));
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            // I/D followed by S
            if has_s && (has_d || has_i) {
                if let Some((range, g1, op, g3)) = match_ior_d_s_end(&this.cigar_str) {
                    let slen = g3.saturating_add(if op == b'I' { g1 } else { 0 });
                    replace_cigar_range(&mut this.cigar_str, range, format_args!("{slen}S"));
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            // S-M(<=10)-I/D
            if has_s && (has_d || has_i) {
                if let Some((range, g1, tmid, g3, op)) = match_begin_s_m_ior_d(&this.cigar_str) {
                    if tmid <= 10 {
                        let slen =
                            g1.saturating_add(tmid)
                                .saturating_add(if op == b'I' { g3 } else { 0 });
                        this.position = this
                            .position
                            .saturating_add(tmid.saturating_add(if op == b'D' { g3 } else { 0 }));
                        replace_cigar_range(&mut this.cigar_str, range, format_args!("{slen}S"));
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                        flag_loop = true;
                    }
                }
            }

            // I/D-M(<=10)-S
            if has_s && (has_d || has_i) {
                if let Some((range, g1, op, tmid, g4)) = match_ior_d_m_s_end(&this.cigar_str) {
                    if tmid <= 10 {
                        let slen =
                            g4.saturating_add(tmid)
                                .saturating_add(if op == b'I' { g1 } else { 0 });
                        replace_cigar_range(&mut this.cigar_str, range, format_args!("{slen}S"));
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                        flag_loop = true;
                    }
                }
            }

            // Single-digit M at start with I/D
            if has_d || has_i {
                if let Some((range, tmid, g2, op, mlen)) =
                    match_begin_digit_m_ior_d_m(&this.cigar_str)
                {
                    flag_loop = this.begin_digit_m_number_ior_d_number_m_range(
                        range,
                        tmid,
                        mlen,
                        op == b'I',
                        op == b'D',
                        g2,
                    );
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                }
            }

            // Single-digit M at end with I/D
            if has_d || has_i {
                if let Some((range, g1, op, tmid)) = match_ior_d_digit_m_end(&this.cigar_str) {
                    let slen = tmid.saturating_add(if op == b'I' { g1 } else { 0 });
                    replace_cigar_range(&mut this.cigar_str, range, format_args!("{slen}S"));
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            // Combine close deletions/insertions into complex
            // Keep outer if/else-if precedence intact so a leftmost match with `mid > 15`
            // still masks the `else` branches (`AR-05`, `AR-07`).
            let two_del_ins = if has_d && has_i {
                find_two_del_ins(&this.cigar_str)
            } else {
                None
            };
            if let Some((range, ov5_len, [g2, g3, g4, g5, g6, g7, g8])) = two_del_ins {
                if g4.saturating_add(g6) <= MAX_MULTI_INDEL_MID_SPAN {
                    flag_loop = this.two_deletions_insertion_to_complex_range(
                        range, ov5_len, g2, g3, g4, g5, g6, g7, g8, flag_loop,
                    );
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                }
            } else {
                let three_del = if has_d {
                    find_three_deletions(&this.cigar_str)
                } else {
                    None
                };
                if let Some((range, ov5_len, [g2, g3, g4, g5, g6, g7, g8])) = three_del {
                    if g4.saturating_add(g6) <= MAX_MULTI_INDEL_MID_SPAN {
                        flag_loop = this.three_deletions_range(
                            range, ov5_len, g2, g3, g4, g5, g6, g7, g8, flag_loop,
                        );
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    }
                } else if has_d || has_i {
                    if let Some(m) = find_three_indels(&this.cigar_str) {
                        if m.g5.saturating_add(m.g8) <= MAX_MULTI_INDEL_MID_SPAN {
                            flag_loop = this.three_indels_range(m, flag_loop);
                            (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                        }
                    }
                }
            }

            if has_d {
                if let Some(m) = find_dig_d_dig_m_dig_di_digi(&this.cigar_str) {
                    if m.mid_len <= MAX_MULTI_INDEL_MID_SPAN {
                        let mut dlen = m.first_len.saturating_add(m.mid_len);
                        let mut ilen = m.mid_len;
                        if m.second_op == b'I' {
                            ilen = ilen.saturating_add(m.second_len);
                        } else if m.second_op == b'D' {
                            dlen = dlen.saturating_add(m.second_len);
                            if let Some(ti) = m.trailing_i {
                                ilen = ilen.saturating_add(ti);
                            }
                        }
                        replace_cigar_range(
                            &mut this.cigar_str,
                            m.range,
                            format_args!("{dlen}D{ilen}I"),
                        );
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                        flag_loop = true;
                    }
                }
            }

            if has_i {
                if let Some(m) = find_notdig_dig_i_dig_m_dig_di_digi(&this.cigar_str) {
                    if m.prefix_ch != 'D'
                        && m.prefix_ch != 'H'
                        && m.mid_len <= MAX_MULTI_INDEL_MID_SPAN
                    {
                        let mut dlen = m.mid_len;
                        let mut ilen = m.first_i_len.saturating_add(m.mid_len);
                        if m.second_op == b'I' {
                            ilen = ilen.saturating_add(m.second_len);
                        } else if m.second_op == b'D' {
                            dlen = dlen.saturating_add(m.second_len);
                            if let Some(ti) = m.trailing_i {
                                ilen = ilen.saturating_add(ti);
                            }
                        }
                        replace_cigar_range(
                            &mut this.cigar_str,
                            m.range_after_prefix,
                            format_args!("{dlen}D{ilen}I"),
                        );
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                        flag_loop = true;
                    }
                }
            }

            if has_d {
                if let Some((range, d1, d2)) = find_adjacent_same_op(&this.cigar_str, b'D') {
                    let dlen = d1.saturating_add(d2);
                    replace_cigar_range(&mut this.cigar_str, range, format_args!("{dlen}D"));
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            if has_i {
                if let Some((range, i1, i2)) = find_adjacent_same_op(&this.cigar_str, b'I') {
                    let ilen = i1.saturating_add(i2);
                    replace_cigar_range(&mut this.cigar_str, range, format_args!("{ilen}I"));
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }
        }

        // End-of-read mismatch -> soft-clip
        if let Some((ov5, mch, soft)) = match_any_number_m_number_s_end(&this.cigar_str) {
            let ov5_len = ov5.len();
            this.capture_mis_softly_ms(ov5_len, mch, soft);
        } else if let Some((ov5, mch)) = match_begin_any_dig_m_end(&this.cigar_str) {
            let ov5_len = ov5.len();
            this.capture_mis_softly_3mismatches(ov5_len, mch);
        }

        // Start-of-read soft-clip recovery
        if let Some((soft, mch)) = match_dig_s_dig_m(&this.cigar_str) {
            this.combine_dig_s_dig_m(soft, mch);
        } else if let Some(mch) = match_begin_dig_m(&this.cigar_str) {
            this.combine_begin_dig_m(mch);
        }

        // ports CigarModifier.java:619-791, CigarParser.java:278
        // Validate that all CIGAR element lengths are strictly positive (> 0).
        // If homopolymer overrun produces non-positive operator lengths (e.g. -1M),
        // Java's TextCigarCodec.decode() throws IllegalArgumentException: Malformed CIGAR string.
        validate_cigar_element_lengths(&this.cigar_str)?;

        Ok(ModifiedCigar {
            position: this.position,
            cigar: std::mem::take(&mut this.cigar_str),
            query_seq: std::mem::take(&mut this.query_seq),
            query_qual: std::mem::take(&mut this.query_qual),
        })
    }

    #[cfg(test)]
    pub fn modify_cigar_legacy(mut self) -> Result<ModifiedCigar, crate::types::VardictError> {
        //flag is set to true if CIGAR string is modified and should be looked at again
        let flag = true;
        let this = &mut self;
        // 1. Strip leading deletion
        if let Some(caps) = BEGIN_NUMBER_D.captures(&this.cigar_str) {
            let dlen = to_int(caps.get(1).unwrap().as_str());
            this.position += dlen;
            this.cigar_str = BEGIN_NUMBER_D.replace(&this.cigar_str, "").into_owned();
        }
        // Strip trailing deletion
        if END_NUMBER_D.is_match(&this.cigar_str) {
            this.cigar_str = END_NUMBER_D.replace(&this.cigar_str, "").into_owned();
        }

        // 2. Replace leading insertion with soft clip
        if let Some(caps) = BEGIN_NUMBER_I.captures(&this.cigar_str) {
            let ilen = to_int(caps.get(1).unwrap().as_str());
            this.cigar_str = BEGIN_NUMBER_I
                .replace(&this.cigar_str, &format!("{}S", ilen))
                .into_owned();
        }
        // Replace trailing insertion with soft clip
        if let Some(caps) = END_NUMBER_I.captures(&this.cigar_str) {
            let ilen = to_int(caps.get(1).unwrap().as_str());
            this.cigar_str = END_NUMBER_I
                .replace(&this.cigar_str, &format!("{}S", ilen))
                .into_owned();
        }

        // 3. Chimeric soft-clip detection
        let sc5_opt = SA_CIGAR_D_S_5CLIP_GROUP.captures(&this.cigar_str);
        let sc3_opt = SA_CIGAR_D_S_3CLIP_GROUP.captures(&this.cigar_str);

        if let Some(caps) = sc5_opt {
            let cigar_element = to_int(caps.get(1).unwrap().as_str()) as usize;
            if !this.config.chimeric && cigar_element >= SEED_2 {
                if let Some(sseq) = subseq(&this.query_seq, 0, cigar_element) {
                    let sequence = rc_seq(sseq);
                    if sequence.len() >= SEED_2 {
                        let rc_seed = &sequence[0..SEED_2];
                        if let Some(ref_pos_u32) = this.unique_seed_pos(rc_seed) {
                            let ref_pos = ref_pos_u32 as i32;
                            if (this.position - ref_pos).abs() < 2 * this.max_read_length {
                                this.cigar_str =
                                    RE_START_S.replace(&this.cigar_str, "").into_owned();
                                if this.query_seq.len() >= cigar_element {
                                    this.query_seq = this.query_seq[cigar_element..].to_vec();
                                    this.query_qual = this.query_qual[cigar_element..].to_vec();
                                }
                                if this.config.verbose {
                                    eprintln!(
                                        "{} at 5' is a chimeric at {} by SEED {}",
                                        String::from_utf8_lossy(&sequence),
                                        this.position,
                                        SEED_2
                                    );
                                }
                            }
                        }
                    }
                }
            }
        } else if let Some(caps) = sc3_opt {
            let cigar_element = to_int(caps.get(1).unwrap().as_str()) as usize;
            if !this.config.chimeric && cigar_element >= SEED_2 {
                let seq_len = this.query_seq.len();
                if seq_len >= cigar_element {
                    let start_idx = seq_len - cigar_element;
                    if let Some(sseq) = subseq(&this.query_seq, start_idx, cigar_element) {
                        let sequence = rc_seq(sseq);
                        if sequence.len() >= SEED_2 {
                            let start_seed = sequence.len() - SEED_2;
                            let rc_seed = &sequence[start_seed..];
                            if let Some(ref_pos_u32) = this.unique_seed_pos(rc_seed) {
                                let ref_pos = ref_pos_u32 as i32;
                                if (this.position - ref_pos).abs() < 2 * this.max_read_length {
                                    this.cigar_str =
                                        RE_END_2S.replace(&this.cigar_str, "").into_owned();
                                    this.query_seq = this.query_seq[0..start_idx].to_vec();
                                    this.query_qual = this.query_qual[0..start_idx].to_vec();
                                    if this.config.verbose {
                                        eprintln!(
                                            "{} at 3' is a chimeric at {} by SEED {}",
                                            String::from_utf8_lossy(&sequence),
                                            this.position,
                                            SEED_2
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // 4. Main indel simplification loop
        let mut flag_loop = flag;
        let mut iter_count = 0;
        let max_iter = this
            .cigar_str
            .bytes()
            .filter(u8::is_ascii_alphabetic)
            .count()
            + 1;
        while flag_loop && this.indel > 0 {
            iter_count += 1;
            if iter_count > max_iter {
                log::warn!(
                    "Infinite loop check triggered in modify_cigar! CIGAR: {}",
                    this.cigar_str
                );
                break;
            }
            flag_loop = false;

            // S followed by I/D
            if let Some(caps) = BEGIN_NUMBER_S_NUMBER_IOR_D.captures(&this.cigar_str) {
                let g1 = to_int(caps.get(1).unwrap().as_str());
                let g2 = to_int(caps.get(2).unwrap().as_str());
                let op = caps.get(3).unwrap().as_str();
                let tslen = format!("{}S", g1 + if op == "I" { g2 } else { 0 });
                this.position += if op == "D" { g2 } else { 0 };
                this.cigar_str = BEGIN_NUMBER_S_NUMBER_IOR_D
                    .replace(&this.cigar_str, &tslen)
                    .into_owned();
                flag_loop = true;
            }

            // I/D followed by S
            if let Some(caps) = NUMBER_IOR_D_NUMBER_S_END.captures(&this.cigar_str) {
                let g1 = to_int(caps.get(1).unwrap().as_str());
                let op = caps.get(2).unwrap().as_str();
                let g3 = to_int(caps.get(3).unwrap().as_str());
                let tslen = format!("{}S", g3 + if op == "I" { g1 } else { 0 });
                this.cigar_str = NUMBER_IOR_D_NUMBER_S_END
                    .replace(&this.cigar_str, &tslen)
                    .into_owned();
                flag_loop = true;
            }

            // S-M(<=10)-I/D
            if let Some(caps) = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IOR_D.captures(&this.cigar_str) {
                let tmid = to_int(caps.get(2).unwrap().as_str());
                if tmid <= 10 {
                    let g1 = to_int(caps.get(1).unwrap().as_str());
                    let g3 = to_int(caps.get(3).unwrap().as_str());
                    let op = caps.get(4).unwrap().as_str();
                    let tslen = format!("{}S", g1 + tmid + if op == "I" { g3 } else { 0 });
                    this.position += tmid + if op == "D" { g3 } else { 0 };
                    this.cigar_str = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IOR_D
                        .replace(&this.cigar_str, &tslen)
                        .into_owned();
                    flag_loop = true;
                }
            }

            // I/D-M(<=10)-S
            if let Some(caps) = NUMBER_IOR_D_NUMBER_M_NUMBER_S_END.captures(&this.cigar_str) {
                let tmid = to_int(caps.get(3).unwrap().as_str());
                if tmid <= 10 {
                    let g1 = to_int(caps.get(1).unwrap().as_str());
                    let op = caps.get(2).unwrap().as_str();
                    let g4 = to_int(caps.get(4).unwrap().as_str());
                    let tslen = format!("{}S", g4 + tmid + if op == "I" { g1 } else { 0 });
                    this.cigar_str = NUMBER_IOR_D_NUMBER_M_NUMBER_S_END
                        .replace(&this.cigar_str, &tslen)
                        .into_owned();
                    flag_loop = true;
                }
            }

            // Single-digit M at start with I/D
            if let Some(caps) = BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M.captures(&this.cigar_str) {
                let tmid = to_int(caps.get(1).unwrap().as_str());
                let mlen = to_int(caps.get(4).unwrap().as_str());
                let op = caps.get(3).unwrap().as_str().to_string();
                let g2 = to_int(caps.get(2).unwrap().as_str());
                flag_loop = this.begin_digit_m_number_ior_d_number_m(tmid, mlen, &op, g2);
            }

            // Single-digit M at end with I/D
            if let Some(caps) = NUMBER_IOR_D_DIGIT_M_END.captures(&this.cigar_str) {
                let tmid = to_int(caps.get(3).unwrap().as_str());
                let op = caps.get(2).unwrap().as_str();
                let g1 = to_int(caps.get(1).unwrap().as_str());
                let tslen = format!("{}S", tmid + if op == "I" { g1 } else { 0 });
                this.cigar_str = NUMBER_IOR_D_NUMBER_M_END
                    .replace(&this.cigar_str, &tslen)
                    .into_owned();
                flag_loop = true;
            }

            // Combine close deletions/insertions into complex
            if let Some(caps) = D_M_D_DD_M_D_I_D_M_D_DD.captures(&this.cigar_str) {
                let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
                let g2 = to_int(caps.get(2).unwrap().as_str());
                let g3 = to_int(caps.get(3).unwrap().as_str());
                let g4 = to_int(caps.get(4).unwrap().as_str());
                let g5 = to_int(caps.get(5).unwrap().as_str());
                let g6 = to_int(caps.get(6).unwrap().as_str());
                let g7 = to_int(caps.get(7).unwrap().as_str());
                let g8 = to_int(caps.get(8).unwrap().as_str());
                flag_loop = this.two_deletions_insertion_to_complex(
                    &ov5, g2, g3, g4, g5, g6, g7, g8, flag_loop,
                );
            } else if let Some(caps) = THREE_DELETIONS_PATTERN.captures(&this.cigar_str) {
                let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
                let g2 = to_int(caps.get(2).unwrap().as_str());
                let g3 = to_int(caps.get(3).unwrap().as_str());
                let g4 = to_int(caps.get(4).unwrap().as_str());
                let g5 = to_int(caps.get(5).unwrap().as_str());
                let g6 = to_int(caps.get(6).unwrap().as_str());
                let g7 = to_int(caps.get(7).unwrap().as_str());
                let g8 = to_int(caps.get(8).unwrap().as_str());
                flag_loop = this.three_deletions(&ov5, g2, g3, g4, g5, g6, g7, g8, flag_loop);
            } else if let Some(caps) = THREE_INDELS_PATTERN.captures(&this.cigar_str) {
                let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
                let g2 = to_int(caps.get(2).unwrap().as_str());
                let g3 = to_int(caps.get(3).unwrap().as_str());
                let op1 = caps.get(4).unwrap().as_str().to_string();
                let g5 = to_int(caps.get(5).unwrap().as_str());
                let g6 = to_int(caps.get(6).unwrap().as_str());
                let op2 = caps.get(7).unwrap().as_str().to_string();
                let g8 = to_int(caps.get(8).unwrap().as_str());
                let g9 = to_int(caps.get(9).unwrap().as_str());
                let op3 = caps.get(10).unwrap().as_str().to_string();
                let g11 = to_int(caps.get(11).unwrap().as_str());
                flag_loop = this.three_indels(
                    &ov5, g2, g3, &op1, g5, g6, &op2, g8, g9, &op3, g11, flag_loop,
                );
            }

            if let Some(caps) = DIG_D_DIG_M_DIG_DI_DIGI.captures(&this.cigar_str) {
                let g1 = to_int(caps.get(1).unwrap().as_str());
                let g2 = to_int(caps.get(2).unwrap().as_str());
                let g3 = to_int(caps.get(3).unwrap().as_str());
                let op = caps.get(4).unwrap().as_str().to_string();
                let istr = caps.get(5).map(|m| m.as_str().to_string());
                flag_loop =
                    this.combine_to_close_to_correct(g1, g2, g3, &op, istr.as_deref(), flag_loop);
            }

            if let Some(caps) = NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI.captures(&this.cigar_str) {
                let op_prefix = caps.get(1).unwrap().as_str().to_string();
                if op_prefix != "D" && op_prefix != "H" {
                    let g2 = to_int(caps.get(2).unwrap().as_str());
                    let g3 = to_int(caps.get(3).unwrap().as_str());
                    let g4 = to_int(caps.get(4).unwrap().as_str());
                    let op = caps.get(5).unwrap().as_str().to_string();
                    let istr = caps.get(6).map(|m| m.as_str().to_string());
                    flag_loop = this.combine_to_close_to_one(
                        &op_prefix,
                        g2,
                        g3,
                        g4,
                        &op,
                        istr.as_deref(),
                        flag_loop,
                    );
                }
            }

            if let Some(caps) = DIG_D_DIG_D.captures(&this.cigar_str) {
                let dlen =
                    to_int(caps.get(1).unwrap().as_str()) + to_int(caps.get(2).unwrap().as_str());
                this.cigar_str = DIG_D_DIG_D
                    .replace(&this.cigar_str, &format!("{}D", dlen))
                    .into_owned();
                flag_loop = true;
            }

            if let Some(caps) = DIG_I_DIG_I.captures(&this.cigar_str) {
                let ilen =
                    to_int(caps.get(1).unwrap().as_str()) + to_int(caps.get(2).unwrap().as_str());
                this.cigar_str = DIG_I_DIG_I
                    .replace(&this.cigar_str, &format!("{}I", ilen))
                    .into_owned();
                flag_loop = true;
            }
        }

        // End-of-read mismatch -> soft-clip
        if let Some(caps) = ANY_NUMBER_M_NUMBER_S_END.captures(&this.cigar_str) {
            let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
            let mch = to_int(caps.get(2).unwrap().as_str());
            let soft = to_int(caps.get(3).unwrap().as_str());
            this.capture_mis_softly_ms_legacy(&ov5, mch, soft);
        } else if let Some(caps) = BEGIN_ANY_DIG_M_END.captures(&this.cigar_str) {
            let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
            let mch = to_int(caps.get(2).unwrap().as_str());
            this.capture_mis_softly_3mismatches_legacy(&ov5, mch);
        }

        // Start-of-read soft-clip recovery
        if let Some(caps) = DIG_S_DIG_M.captures(&this.cigar_str) {
            let soft = to_int(caps.get(1).unwrap().as_str());
            let mch = to_int(caps.get(2).unwrap().as_str());
            this.combine_dig_s_dig_m_legacy(soft, mch);
        } else if let Some(caps) = BEGIN_DIG_M.captures(&this.cigar_str) {
            let mch = to_int(caps.get(1).unwrap().as_str());
            this.combine_begin_dig_m_legacy(mch);
        }

        // ports CigarModifier.java:619-791, CigarParser.java:278
        // Validate that all CIGAR element lengths are strictly positive (> 0).
        // If homopolymer overrun produces non-positive operator lengths (e.g. -1M),
        // Java's TextCigarCodec.decode() throws IllegalArgumentException: Malformed CIGAR string.
        for cap in CIGAR_ELEMENT_PATTERN.captures_iter(&self.cigar_str) {
            let len: i32 = cap[1].parse().unwrap_or(0);
            if len <= 0 {
                return Err(crate::types::VardictError::JavaThrows(format!(
                    "invalid CIGAR length {} in {}",
                    len, self.cigar_str
                )));
            }
        }

        Ok(ModifiedCigar {
            position: self.position,
            cigar: self.cigar_str,
            query_seq: self.query_seq,
            query_qual: self.query_qual,
        })
    }

    fn combine_begin_dig_m(&mut self, mut mch: i32) {
        let rn = count_start_mismatch_clip(
            self.position,
            mch,
            0,
            &self.query_seq,
            self.reference,
            self.ref_start,
        );
        if rn > 0 && rn <= 3 {
            mch = mch.saturating_sub(rn);
            if let Some((_, b'M', _, end)) = next_cigar_op(&self.cigar_str, 0) {
                replace_cigar_range(&mut self.cigar_str, 0..end, format_args!("{rn}S{mch}M"));
            }
            self.position = self.position.saturating_add(rn);
        }
    }

    #[cfg(test)]
    fn combine_begin_dig_m_legacy(&mut self, mut mch: i32) {
        let rn = count_start_mismatch_clip(
            self.position,
            mch,
            0,
            &self.query_seq,
            self.reference,
            self.ref_start,
        );
        if rn > 0 && rn <= 3 {
            mch = mch.saturating_sub(rn);
            self.cigar_str = BEGIN_DIG_M
                .replace(&self.cigar_str, &format!("{}S{}M", rn, mch))
                .into_owned();
            self.position = self.position.saturating_add(rn);
        }
    }

    fn combine_dig_s_dig_m(&mut self, soft: i32, mch: i32) {
        let (new_pos, cigar_soft, cigar_mch, changed) = step_combine_dig_s_dig_m(
            self.position,
            soft,
            mch,
            &self.query_seq,
            &self.query_qual,
            self.reference,
            self.ref_start,
        );
        self.position = new_pos;
        if changed {
            if let Some((_, b'S', _, idx1)) = next_cigar_op(&self.cigar_str, 0) {
                if let Some((_, b'M', _, end)) = next_cigar_op(&self.cigar_str, idx1) {
                    let span = 0..end;
                    if cigar_soft > 0 {
                        replace_cigar_range(
                            &mut self.cigar_str,
                            span,
                            format_args!("{cigar_soft}S{cigar_mch}M"),
                        );
                    } else {
                        replace_cigar_range(
                            &mut self.cigar_str,
                            span,
                            format_args!("{cigar_mch}M"),
                        );
                    }
                }
            }
        }
    }

    #[cfg(test)]
    fn combine_dig_s_dig_m_legacy(&mut self, soft: i32, mch: i32) {
        let (new_pos, cigar_soft, cigar_mch, changed) = step_combine_dig_s_dig_m(
            self.position,
            soft,
            mch,
            &self.query_seq,
            &self.query_qual,
            self.reference,
            self.ref_start,
        );
        self.position = new_pos;
        if changed {
            if cigar_soft > 0 {
                self.cigar_str = DIG_S_DIG_M
                    .replace(&self.cigar_str, &format!("{cigar_soft}S{cigar_mch}M"))
                    .into_owned();
            } else {
                self.cigar_str = DIG_S_DIG_M
                    .replace(&self.cigar_str, &format!("{cigar_mch}M"))
                    .into_owned();
            }
        }
    }

    fn capture_mis_softly_3mismatches(&mut self, ov5_len: usize, mut mch: i32) {
        let mut refoff = self.position.saturating_add(mch);
        let mut rdoff = mch;
        if ov5_len > 0 {
            let (mnd, mis) = sum_mnd_and_mis_prefix(self.cigar_str.get(..ov5_len).unwrap_or(""));
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }
        let rn = count_end_mismatch_clip(
            refoff,
            rdoff,
            mch,
            &self.query_seq,
            self.reference,
            self.ref_start,
        );
        mch -= rn;
        if rn > 0 && rn <= 3 {
            let end = self.cigar_str.len();
            replace_cigar_range(
                &mut self.cigar_str,
                ov5_len..end,
                format_args!("{mch}M{rn}S"),
            );
        }
    }

    #[cfg(test)]
    fn capture_mis_softly_3mismatches_legacy(&mut self, ov5: &str, mut mch: i32) {
        let mut refoff = self.position + mch;
        let mut rdoff = mch;
        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }
        let rn = count_end_mismatch_clip(
            refoff,
            rdoff,
            mch,
            &self.query_seq,
            self.reference,
            self.ref_start,
        );
        mch -= rn;
        if rn > 0 && rn <= 3 {
            self.cigar_str = DIG_M_END
                .replace(&self.cigar_str, &format!("{}M{}S", mch, rn))
                .into_owned();
        }
    }

    fn capture_mis_softly_ms(&mut self, ov5_len: usize, mch: i32, soft: i32) {
        let mut refoff = self.position.saturating_add(mch);
        let mut rdoff = mch;
        if ov5_len > 0 {
            let (mnd, mis) = sum_mnd_and_mis_prefix(self.cigar_str.get(..ov5_len).unwrap_or(""));
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }
        let (new_mch, new_soft, changed) = step_capture_mis_softly_ms(
            refoff,
            rdoff,
            mch,
            soft,
            &self.query_seq,
            &self.query_qual,
            self.reference,
            self.ref_start,
        );
        if changed {
            let end = self.cigar_str.len();
            let span = ov5_len..end;
            if new_soft > 0 {
                replace_cigar_range(
                    &mut self.cigar_str,
                    span,
                    format_args!("{new_mch}M{new_soft}S"),
                );
            } else {
                replace_cigar_range(&mut self.cigar_str, span, format_args!("{new_mch}M"));
            }
        }
    }

    #[cfg(test)]
    fn capture_mis_softly_ms_legacy(&mut self, ov5: &str, mch: i32, soft: i32) {
        let mut refoff = self.position + mch;
        let mut rdoff = mch;
        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }
        let (new_mch, new_soft, changed) = step_capture_mis_softly_ms(
            refoff,
            rdoff,
            mch,
            soft,
            &self.query_seq,
            &self.query_qual,
            self.reference,
            self.ref_start,
        );
        if changed {
            if new_soft > 0 {
                self.cigar_str = DIG_M_DIG_S_END
                    .replace(&self.cigar_str, &format!("{new_mch}M{new_soft}S"))
                    .into_owned();
            } else {
                self.cigar_str = DIG_M_DIG_S_END
                    .replace(&self.cigar_str, &format!("{new_mch}M"))
                    .into_owned();
            }
        }
    }

    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    fn combine_to_close_to_one(
        &mut self,
        op_prefix: &str,
        g2: i32,
        g3: i32,
        g4: i32,
        op: &str,
        istr: Option<&str>,
        mut flag: bool,
    ) -> bool {
        if g3 <= MAX_MULTI_INDEL_MID_SPAN {
            let mut dlen = g3;
            let mut ilen = g2.saturating_add(g3);
            if op == "I" {
                ilen = ilen.saturating_add(g4);
            } else if op == "D" {
                dlen = dlen.saturating_add(g4);
                if let Some(istr_str) = istr {
                    if !istr_str.is_empty() {
                        if let Some(prefix) = istr_str.get(..istr_str.len() - 1) {
                            ilen = ilen.saturating_add(to_int(prefix));
                        }
                    }
                }
            }
            // Replace matching pattern
            let replace_str = format!("{}{:?}D{:?}I", op_prefix, dlen, ilen);
            self.cigar_str = NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI
                .replace(&self.cigar_str, &replace_str)
                .into_owned();
            flag = true;
        }
        flag
    }

    #[cfg(test)]
    fn combine_to_close_to_correct(
        &mut self,
        g1: i32,
        g2: i32,
        g3: i32,
        op: &str,
        istr: Option<&str>,
        mut flag: bool,
    ) -> bool {
        if g2 <= MAX_MULTI_INDEL_MID_SPAN {
            let mut dlen = g1.saturating_add(g2);
            let mut ilen = g2;
            if op == "I" {
                ilen = ilen.saturating_add(g3);
            } else if op == "D" {
                dlen = dlen.saturating_add(g3);
                if let Some(istr_str) = istr {
                    if !istr_str.is_empty() {
                        if let Some(prefix) = istr_str.get(..istr_str.len() - 1) {
                            ilen = ilen.saturating_add(to_int(prefix));
                        }
                    }
                }
            }
            self.cigar_str = DIG_D_DIG_M_DIG_DI_DIGI
                .replace(&self.cigar_str, &format!("{}D{}I", dlen, ilen))
                .into_owned();
            flag = true;
        }
        flag
    }

    fn three_indels_range(&mut self, m: ThreeIndelsMatch, flag: bool) -> bool {
        let ThreeIndelsMatch {
            range,
            ov5_len,
            g2,
            g3,
            op1,
            g5,
            g6,
            op2,
            g8,
            g9,
            op3,
            g11,
        } = m;
        let mid = g5.saturating_add(g8);
        if mid > MAX_MULTI_INDEL_MID_SPAN {
            return flag;
        }
        let mut tslen = mid;
        if op1 == b'I' {
            tslen = tslen.saturating_add(g3);
        }
        if op2 == b'I' {
            tslen = tslen.saturating_add(g6);
        }
        if op3 == b'I' {
            tslen = tslen.saturating_add(g9);
        }

        let mut dlen = mid;
        if op1 == b'D' {
            dlen = dlen.saturating_add(g3);
        }
        if op2 == b'D' {
            dlen = dlen.saturating_add(g6);
        }
        if op3 == b'D' {
            dlen = dlen.saturating_add(g9);
        }

        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g11;

        if ov5_len > 0 {
            let (mnd, mis) = sum_mnd_and_mis_prefix(self.cigar_str.get(..ov5_len).unwrap_or(""));
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }

        let mut rn = 0;
        while rdoff.saturating_add(rn) < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get(rdoff.saturating_add(rn) as usize) {
                if is_has_and_equals_char(
                    q_base,
                    self.reference,
                    refoff.saturating_add(rn),
                    self.ref_start,
                ) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr = rdoff_corr.saturating_add(rn);
        dlen = dlen.saturating_sub(rn);
        tslen = tslen.saturating_sub(rn);

        if tslen <= 0 {
            dlen = dlen.saturating_sub(tslen);
            rm = rm.saturating_add(tslen);
            if dlen == 0 {
                rdoff_corr = rdoff_corr.saturating_add(rm);
                replace_cigar_range(&mut self.cigar_str, range, format_args!("{rdoff_corr}M"));
            } else if dlen < 0 {
                tslen = dlen.saturating_neg();
                rm = rm.saturating_add(dlen);
                if rm < 0 {
                    rdoff_corr = rdoff_corr.saturating_add(rm);
                    replace_cigar_range(
                        &mut self.cigar_str,
                        range,
                        format_args!("{rdoff_corr}M{tslen}I"),
                    );
                } else {
                    replace_cigar_range(
                        &mut self.cigar_str,
                        range,
                        format_args!("{rdoff_corr}M{tslen}I{rm}M"),
                    );
                }
            } else {
                replace_cigar_range(
                    &mut self.cigar_str,
                    range,
                    format_args!("{rdoff_corr}M{dlen}D{rm}M"),
                );
            }
        } else if dlen == 0 {
            replace_cigar_range(
                &mut self.cigar_str,
                range,
                format_args!("{rdoff_corr}M{tslen}I{rm}M"),
            );
        } else if dlen < 0 {
            rm = rm.saturating_add(dlen);
            replace_cigar_range(
                &mut self.cigar_str,
                range,
                format_args!("{rdoff_corr}M{tslen}I{rm}M"),
            );
        } else {
            replace_cigar_range(
                &mut self.cigar_str,
                range,
                format_args!("{rdoff_corr}M{dlen}D{tslen}I{rm}M"),
            );
        }
        true
    }

    fn three_deletions_range(
        &mut self,
        range: std::ops::Range<usize>,
        ov5_len: usize,
        g2: i32,
        g3: i32,
        g4: i32,
        g5: i32,
        g6: i32,
        g7: i32,
        g8: i32,
        flag: bool,
    ) -> bool {
        let mid = g4.saturating_add(g6);
        if mid > MAX_MULTI_INDEL_MID_SPAN {
            return flag;
        }
        let mut tslen = mid;
        let mut dlen = g3
            .saturating_add(g4)
            .saturating_add(g5)
            .saturating_add(g6)
            .saturating_add(g7);
        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if ov5_len > 0 {
            let (mnd, mis) = sum_mnd_and_mis_prefix(self.cigar_str.get(..ov5_len).unwrap_or(""));
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }
        let mut rn = 0;
        while rdoff.saturating_add(rn) < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get(rdoff.saturating_add(rn) as usize) {
                if is_has_and_equals_char(
                    q_base,
                    self.reference,
                    refoff.saturating_add(rn),
                    self.ref_start,
                ) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr = rdoff_corr.saturating_add(rn);
        dlen = dlen.saturating_sub(rn);
        tslen = tslen.saturating_sub(rn);
        if tslen <= 0 {
            dlen = dlen.saturating_sub(tslen);
            rm = rm.saturating_add(tslen);
            replace_cigar_range(
                &mut self.cigar_str,
                range,
                format_args!("{rdoff_corr}M{dlen}D{rm}M"),
            );
        } else {
            replace_cigar_range(
                &mut self.cigar_str,
                range,
                format_args!("{rdoff_corr}M{dlen}D{tslen}I{rm}M"),
            );
        }
        true
    }

    fn two_deletions_insertion_to_complex_range(
        &mut self,
        range: std::ops::Range<usize>,
        ov5_len: usize,
        g2: i32,
        g3: i32,
        g4: i32,
        g5: i32,
        g6: i32,
        g7: i32,
        g8: i32,
        flag: bool,
    ) -> bool {
        let mid = g4.saturating_add(g6);
        if mid > MAX_MULTI_INDEL_MID_SPAN {
            return flag;
        }
        let mut tslen = mid.saturating_add(g5);
        let mut dlen = g3.saturating_add(mid).saturating_add(g7);
        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if ov5_len > 0 {
            let (mnd, mis) = sum_mnd_and_mis_prefix(self.cigar_str.get(..ov5_len).unwrap_or(""));
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }

        let mut rn = 0;
        while rdoff.saturating_add(rn) < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get(rdoff.saturating_add(rn) as usize) {
                if is_has_and_equals_char(
                    q_base,
                    self.reference,
                    refoff.saturating_add(rn),
                    self.ref_start,
                ) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr = rdoff_corr.saturating_add(rn);
        dlen = dlen.saturating_sub(rn);
        tslen = tslen.saturating_sub(rn);
        if tslen <= 0 {
            dlen = dlen.saturating_sub(tslen);
            rm = rm.saturating_add(tslen);
            replace_cigar_range(
                &mut self.cigar_str,
                range,
                format_args!("{rdoff_corr}M{dlen}D{rm}M"),
            );
        } else {
            replace_cigar_range(
                &mut self.cigar_str,
                range,
                format_args!("{rdoff_corr}M{dlen}D{tslen}I{rm}M"),
            );
        }
        true
    }

    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    fn three_indels(
        &mut self,
        ov5: &str,
        g2: i32,
        g3: i32,
        op1: &str,
        g5: i32,
        g6: i32,
        op2: &str,
        g8: i32,
        g9: i32,
        op3: &str,
        g11: i32,
        flag: bool,
    ) -> bool {
        let mid = g5.saturating_add(g8);
        if mid > MAX_MULTI_INDEL_MID_SPAN {
            return flag;
        }
        let mut tslen = mid;
        if op1 == "I" {
            tslen = tslen.saturating_add(g3);
        }
        if op2 == "I" {
            tslen = tslen.saturating_add(g6);
        }
        if op3 == "I" {
            tslen = tslen.saturating_add(g9);
        }

        let mut dlen = mid;
        if op1 == "D" {
            dlen = dlen.saturating_add(g3);
        }
        if op2 == "D" {
            dlen = dlen.saturating_add(g6);
        }
        if op3 == "D" {
            dlen = dlen.saturating_add(g9);
        }

        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g11;

        if !ov5.is_empty() {
            let (mnd, mis) = sum_mnd_and_mis_prefix(ov5);
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }

        let mut rn = 0;
        while rdoff.saturating_add(rn) < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get(rdoff.saturating_add(rn) as usize) {
                if is_has_and_equals_char(
                    q_base,
                    self.reference,
                    refoff.saturating_add(rn),
                    self.ref_start,
                ) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr = rdoff_corr.saturating_add(rn);
        dlen = dlen.saturating_sub(rn);
        tslen = tslen.saturating_sub(rn);

        let mut new_cigar = format!("{}M", rdoff_corr);
        if tslen <= 0 {
            dlen = dlen.saturating_sub(tslen);
            rm = rm.saturating_add(tslen);
            if dlen == 0 {
                rdoff_corr = rdoff_corr.saturating_add(rm);
                new_cigar = format!("{}M", rdoff_corr);
            } else if dlen < 0 {
                tslen = dlen.saturating_neg();
                rm = rm.saturating_add(dlen);
                if rm < 0 {
                    rdoff_corr = rdoff_corr.saturating_add(rm);
                    new_cigar = format!("{}M{}I", rdoff_corr, tslen);
                } else {
                    new_cigar = format!("{}{}I{}M", new_cigar, tslen, rm);
                }
            } else {
                new_cigar = format!("{}{}D{}M", new_cigar, dlen, rm);
            }
        } else if dlen == 0 {
            new_cigar = format!("{}{}I{}M", new_cigar, tslen, rm);
        } else if dlen < 0 {
            rm = rm.saturating_add(dlen);
            new_cigar = format!("{}{}I{}M", new_cigar, tslen, rm);
        } else {
            new_cigar = format!("{}{}D{}I{}M", new_cigar, dlen, tslen, rm);
        }

        self.cigar_str = DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM
            .replace(&self.cigar_str, &new_cigar)
            .into_owned();
        true
    }

    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    fn three_deletions(
        &mut self,
        ov5: &str,
        g2: i32,
        g3: i32,
        g4: i32,
        g5: i32,
        g6: i32,
        g7: i32,
        g8: i32,
        flag: bool,
    ) -> bool {
        let mid = g4.saturating_add(g6);
        if mid > MAX_MULTI_INDEL_MID_SPAN {
            return flag;
        }
        let mut tslen = mid;
        let mut dlen = g3
            .saturating_add(g4)
            .saturating_add(g5)
            .saturating_add(g6)
            .saturating_add(g7);
        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if !ov5.is_empty() {
            let (mnd, mis) = sum_mnd_and_mis_prefix(ov5);
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }
        let mut rn = 0;
        while rdoff.saturating_add(rn) < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get(rdoff.saturating_add(rn) as usize) {
                if is_has_and_equals_char(
                    q_base,
                    self.reference,
                    refoff.saturating_add(rn),
                    self.ref_start,
                ) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr = rdoff_corr.saturating_add(rn);
        dlen = dlen.saturating_sub(rn);
        tslen = tslen.saturating_sub(rn);
        let mut new_cigar = format!("{}M", rdoff_corr);
        if tslen <= 0 {
            dlen = dlen.saturating_sub(tslen);
            rm = rm.saturating_add(tslen);
            new_cigar = format!("{}M{}D{}M", rdoff_corr, dlen, rm);
        } else {
            new_cigar = format!("{}{}D{}I{}M", new_cigar, dlen, tslen, rm);
        }

        self.cigar_str = DM_DD_DM_DD_DM_DD_DM
            .replace(&self.cigar_str, &new_cigar)
            .into_owned();
        true
    }

    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    fn two_deletions_insertion_to_complex(
        &mut self,
        ov5: &str,
        g2: i32,
        g3: i32,
        g4: i32,
        g5: i32,
        g6: i32,
        g7: i32,
        g8: i32,
        flag: bool,
    ) -> bool {
        let mid = g4.saturating_add(g6);
        if mid > MAX_MULTI_INDEL_MID_SPAN {
            return flag;
        }
        let mut tslen = mid.saturating_add(g5);
        let mut dlen = g3.saturating_add(mid).saturating_add(g7);
        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if !ov5.is_empty() {
            let (mnd, mis) = sum_mnd_and_mis_prefix(ov5);
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }

        let mut rn = 0;
        while rdoff.saturating_add(rn) < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get(rdoff.saturating_add(rn) as usize) {
                if is_has_and_equals_char(
                    q_base,
                    self.reference,
                    refoff.saturating_add(rn),
                    self.ref_start,
                ) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr = rdoff_corr.saturating_add(rn);
        dlen = dlen.saturating_sub(rn);
        tslen = tslen.saturating_sub(rn);
        let mut new_cigar = format!("{}M", rdoff_corr);
        if tslen <= 0 {
            dlen = dlen.saturating_sub(tslen);
            rm = rm.saturating_add(tslen);
            new_cigar = format!("{}M{}D{}M", rdoff_corr, dlen, rm);
        } else {
            new_cigar = format!("{}{}D{}I{}M", new_cigar, dlen, tslen, rm);
        }

        self.cigar_str = D_M_D_DD_M_D_I_D_M_D_DD_PRIM
            .replace(&self.cigar_str, &new_cigar)
            .into_owned();
        true
    }

    fn begin_digit_m_number_ior_d_number_m_range(
        &mut self,
        range: std::ops::Range<usize>,
        tmid: i32,
        mut mlen: i32,
        is_ins: bool,
        is_del: bool,
        g2: i32,
    ) -> bool {
        let mut tn = 0;

        let mut tslen = tmid.saturating_add(if is_ins { g2 } else { 0 });
        self.position =
            self.position
                .saturating_add(tmid.saturating_add(if is_del { g2 } else { 0 }));

        while tn < mlen {
            if let Some(&q_base) = self.query_seq.get(tslen.saturating_add(tn) as usize) {
                if is_has_and_not_equals_char(
                    q_base,
                    self.reference,
                    self.position.saturating_add(tn),
                    self.ref_start,
                ) {
                    tn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        tslen = tslen.saturating_add(tn);
        mlen -= tn;
        self.position = self.position.saturating_add(tn);

        replace_cigar_range(&mut self.cigar_str, range, format_args!("{tslen}S{mlen}M"));
        true
    }

    #[cfg(test)]
    fn begin_digit_m_number_ior_d_number_m(
        &mut self,
        tmid: i32,
        mut mlen: i32,
        op: &str,
        g2: i32,
    ) -> bool {
        let mut tn = 0;

        let mut tslen = tmid.saturating_add(if op == "I" { g2 } else { 0 });
        self.position = self
            .position
            .saturating_add(tmid.saturating_add(if op == "D" { g2 } else { 0 }));

        while tn < mlen {
            if let Some(&q_base) = self.query_seq.get(tslen.saturating_add(tn) as usize) {
                if is_has_and_not_equals_char(
                    q_base,
                    self.reference,
                    self.position.saturating_add(tn),
                    self.ref_start,
                ) {
                    tn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        tslen = tslen.saturating_add(tn);
        mlen -= tn;
        self.position = self.position.saturating_add(tn);

        self.cigar_str = BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M_OPT
            .replace(&self.cigar_str, &format!("{}S{}M", tslen, mlen))
            .into_owned();
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;

    fn get_default_config() -> Config {
        Config::default()
    }

    #[test]
    fn strip_del() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let reference = b"";
        // 1. Starts with D: should adjust position and strip leading D.
        let cm = CigarModifier::new(
            10,
            "3D10M5D",
            b"ACGTACGTAC",
            b"IIIIIIIIII",
            reference,
            1,
            &seed,
            1,
            150,
            &config,
        );
        let res = cm.modify_cigar();
        assert_eq!(res.as_ref().unwrap().position, 13);
        assert_eq!(res.as_ref().unwrap().cigar, "10M");
    }

    #[test]
    fn ins_to_sclip() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let reference = b"";
        // 2. Starts or ends with I: should replace with S.
        let cm = CigarModifier::new(
            10,
            "3I10M2I",
            b"ACGTACGTACGTACG",
            b"IIIIIIIIIIIIIII",
            reference,
            1,
            &seed,
            1,
            150,
            &config,
        );
        let res = cm.modify_cigar();
        assert_eq!(res.as_ref().unwrap().position, 10);
        assert_eq!(res.as_ref().unwrap().cigar, "3S10M2S");
    }

    #[test]
    fn lead_trail_softclip_indel() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let reference = b"";
        // Starts with softclip followed by insertion or deletion: should combine them.
        let cm = CigarModifier::new(
            10,
            "3S2I10M",
            b"ACGTACGTACGTACG",
            b"IIIIIIIIIIIIIII",
            reference,
            1,
            &seed,
            1,
            150,
            &config,
        );
        let res = cm.modify_cigar();
        assert_eq!(res.as_ref().unwrap().position, 10);
        assert_eq!(res.as_ref().unwrap().cigar, "5S10M");

        let cm = CigarModifier::new(
            10,
            "3S2D10M",
            b"ACGTACGTACGTACG",
            b"IIIIIIIIIIIIIII",
            reference,
            1,
            &seed,
            1,
            150,
            &config,
        );
        let res = cm.modify_cigar();
        assert_eq!(res.as_ref().unwrap().position, 12);
        assert_eq!(res.as_ref().unwrap().cigar, "3S10M");
    }

    #[test]
    fn double_dels_ins() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let reference = b"";
        // Combine adjacent deletions
        let cm = CigarModifier::new(
            10,
            "10M2D3D10M",
            b"ACGTACGTACGTACG",
            b"IIIIIIIIIIIIIII",
            reference,
            1,
            &seed,
            1,
            150,
            &config,
        );
        let res = cm.modify_cigar();
        assert_eq!(res.as_ref().unwrap().cigar, "10M5D10M");

        // Combine adjacent insertions
        let cm = CigarModifier::new(
            10,
            "10M2I3I10M",
            b"ACGTACGTACGTACG",
            b"IIIIIIIIIIIIIII",
            reference,
            1,
            &seed,
            1,
            150,
            &config,
        );
        let res = cm.modify_cigar();
        assert_eq!(res.as_ref().unwrap().cigar, "10M5I10M");
    }

    #[test]
    fn crash_corpus_fixes() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let reference = b"";
        // Feed various weird CIGAR strings that used to crash VarDictJava due to bounds check errors
        let inputs = vec!["0M", "100D", "100I", "5M100D5M", "5M100I5M", "1M1I1D1M"];
        for cig in inputs {
            let cm = CigarModifier::new(
                10,
                cig,
                b"ACGTACGTACGTACG",
                b"IIIIIIIIIIIIIII",
                reference,
                1,
                &seed,
                1,
                150,
                &config,
            );
            let _res = cm.modify_cigar();
            // Verify that it doesn't crash/panic!
        }
    }

    #[test]
    fn test_pr402_cigar_java_throws() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let query_seq = b"AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTACAAAAAAAAAAAAAAAAAAAAGAAAAAAGAAAAAAAAAG";
        let query_qual = b"BBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF=74??44444444FFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        let ref_seq = b"AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTACAAAAAAAAAAAGAAAAAAAAAAGAAAAAAAAAGAAAAAAGAAAAAAAAAGAAAGAAAGAAAAAGAAAAAACAGAAAGAAGAATGAAAATCTGGGTG";

        let cm = CigarModifier::new(
            1,
            "104M1I2M2D1M1D5M1D2M9D27M",
            query_seq,
            query_qual,
            ref_seq,
            1,
            &seed,
            3,
            150,
            &config,
        );

        let res = cm.modify_cigar();
        assert!(
            res.is_err(),
            "PR #402 read should trigger JavaThrows error: {:?}",
            res
        );
        match res.unwrap_err() {
            crate::types::VardictError::JavaThrows(msg) => {
                assert!(msg.contains("invalid"), "Unexpected message: {msg}");
            }
            other => panic!("Expected JavaThrows, got: {:?}", other),
        }
    }

    #[test]
    fn test_siblings_non_positive_operator_lengths_trigger_java_throws() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let ref_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let query_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let query_qual = vec![b'I'; query_seq.len()];

        // Adversarial inputs for three_deletions, two_deletions_insertion_to_complex, three_indels
        // with short M components that cause non-positive operator lengths.
        let adversarial_cases = ["10M1D1M1D1M1D1M", "10M1I1M1D1M1D1M", "10M1D1M1I1M1D1M"];

        for cigar in adversarial_cases {
            let cm = CigarModifier::new(
                1,
                cigar,
                query_seq,
                query_qual.clone(),
                ref_seq,
                1,
                &seed,
                3,
                150,
                &config,
            );
            let res = cm.modify_cigar();
            // If the modification produces an invalid CIGAR, it must return Err(JavaThrows).
            // If it succeeds, every CIGAR operator length must be strictly positive.
            match res {
                Ok(modified) => {
                    let parsed = crate::cigar::parse_cigar(&modified.cigar)
                        .unwrap_or_else(|| panic!("Failed to parse CIGAR {}", modified.cigar));
                    for (len, _op) in parsed {
                        assert!(
                            len > 0,
                            "CIGAR {} produced non-positive length in {}",
                            cigar,
                            modified.cigar
                        );
                    }
                }
                Err(crate::types::VardictError::JavaThrows(msg)) => {
                    assert!(
                        msg.contains("invalid CIGAR length"),
                        "Expected invalid CIGAR length error for {cigar}, got: {msg}"
                    );
                }
                Err(other) => panic!("Unexpected error for {cigar}: {:?}", other),
            }
        }
    }

    fn assert_byte_matchers_equal_regex(s: &str) {
        let expected_ms = ANY_NUMBER_M_NUMBER_S_END.captures(s).map(|caps| {
            (
                caps.get(1).map_or("", |m| m.as_str()),
                to_int(caps.get(2).unwrap().as_str()),
                to_int(caps.get(3).unwrap().as_str()),
            )
        });
        assert_eq!(
            match_any_number_m_number_s_end(s),
            expected_ms,
            "ANY_NUMBER_M_NUMBER_S_END mismatch on {s:?}"
        );

        let expected_m_end = BEGIN_ANY_DIG_M_END.captures(s).map(|caps| {
            (
                caps.get(1).map_or("", |m| m.as_str()),
                to_int(caps.get(2).unwrap().as_str()),
            )
        });
        assert_eq!(
            match_begin_any_dig_m_end(s),
            expected_m_end,
            "BEGIN_ANY_DIG_M_END mismatch on {s:?}"
        );

        let expected_sm = DIG_S_DIG_M.captures(s).map(|caps| {
            (
                to_int(caps.get(1).unwrap().as_str()),
                to_int(caps.get(2).unwrap().as_str()),
            )
        });
        assert_eq!(
            match_dig_s_dig_m(s),
            expected_sm,
            "DIG_S_DIG_M mismatch on {s:?}"
        );

        let expected_bm = BEGIN_DIG_M
            .captures(s)
            .map(|caps| to_int(caps.get(1).unwrap().as_str()));
        assert_eq!(
            match_begin_dig_m(s),
            expected_bm,
            "BEGIN_DIG_M mismatch on {s:?}"
        );

        let expected_valid: Result<(), String> = (|| {
            for cap in CIGAR_ELEMENT_PATTERN.captures_iter(s) {
                let len: i32 = cap[1].parse().unwrap_or(0);
                if len <= 0 {
                    return Err(format!("invalid CIGAR length {} in {}", len, s));
                }
            }
            Ok(())
        })();
        let actual_valid: Result<(), String> =
            validate_cigar_element_lengths(s).map_err(|e| match e {
                crate::types::VardictError::JavaThrows(msg) => msg,
                other => format!("{other:?}"),
            });
        assert_eq!(
            actual_valid, expected_valid,
            "CIGAR_ELEMENT_PATTERN mismatch on {s:?}"
        );

        let expected_lead2s = SA_CIGAR_D_S_5CLIP_GROUP
            .captures(s)
            .map(|caps| to_int(caps.get(1).unwrap().as_str()) as usize);
        assert_eq!(
            leading_2digit_s_len(s),
            expected_lead2s,
            "SA_CIGAR_D_S_5CLIP_GROUP mismatch on {s:?}"
        );

        let expected_trail2s = SA_CIGAR_D_S_3CLIP_GROUP
            .captures(s)
            .map(|caps| to_int(caps.get(1).unwrap().as_str()) as usize);
        assert_eq!(
            trailing_2digit_s_len(s),
            expected_trail2s,
            "SA_CIGAR_D_S_3CLIP_GROUP mismatch on {s:?}"
        );

        let expected_sums = (
            global_find_sum(&ALIGNED_LENGTH_MND, s),
            global_find_sum(&SOFT_CLIPPED, s),
        );
        assert_eq!(
            sum_mnd_and_mis_prefix(s),
            expected_sums,
            "sum_mnd_and_mis_prefix mismatch on {s:?}"
        );
    }

    #[allow(clippy::too_many_arguments)]
    fn assert_modify_cigar_matches_legacy(
        position: i32,
        cigar: &str,
        query_seq: &[u8],
        query_qual: &[u8],
        reference: &[u8],
        ref_start: i32,
        seed: &SeedMap,
        indel: i32,
        max_read_length: i32,
        config: &Config,
    ) {
        let cm_new = CigarModifier::new(
            position,
            cigar,
            query_seq,
            query_qual,
            reference,
            ref_start,
            seed,
            indel,
            max_read_length,
            config,
        );
        let cm_leg = CigarModifier::new(
            position,
            cigar,
            query_seq,
            query_qual,
            reference,
            ref_start,
            seed,
            indel,
            max_read_length,
            config,
        );
        let res_new = cm_new.modify_cigar();
        let res_leg = cm_leg.modify_cigar_legacy();
        match (&res_new, &res_leg) {
            (Ok(a), Ok(b)) => {
                assert_eq!(a.position, b.position, "position mismatch for {cigar}");
                assert_eq!(a.cigar, b.cigar, "cigar mismatch for {cigar}");
                assert_eq!(a.query_seq, b.query_seq, "query_seq mismatch for {cigar}");
                assert_eq!(
                    a.query_qual, b.query_qual,
                    "query_qual mismatch for {cigar}"
                );
            }
            (
                Err(crate::types::VardictError::JavaThrows(msg_a)),
                Err(crate::types::VardictError::JavaThrows(msg_b)),
            ) => {
                assert_eq!(msg_a, msg_b, "error message mismatch for {cigar}");
            }
            _ => panic!("result mismatch for {cigar}: new={res_new:?} vs legacy={res_leg:?}"),
        }
    }

    #[test]
    fn modify_cigar_byte_matchers_and_full_pipeline_match_legacy() {
        // 1. Exhaustive check over all 72^1 + 72^2 + 72^3 (378,504) 1-, 2-, and 3-element
        // CIGAR strings from ops {M, I, D, N, S, H, P, =, X} x lengths {0, 1, 2, 9, 10, 11, 100, 99999999999}
        // plus negative lengths (-1, -10, "10M-5M", "-").
        let ops = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];
        let lens = ["0", "1", "2", "9", "10", "11", "100", "99999999999"];
        let mut elems = Vec::with_capacity(ops.len() * lens.len());
        for &len in &lens {
            for &op in &ops {
                elems.push(format!("{len}{op}"));
            }
        }
        assert_eq!(elems.len(), 72);

        let mut buf = String::with_capacity(48);
        for e1 in &elems {
            assert_byte_matchers_equal_regex(e1);
            for e2 in &elems {
                buf.clear();
                buf.push_str(e1);
                buf.push_str(e2);
                assert_byte_matchers_equal_regex(&buf);
                for e3 in &elems {
                    buf.truncate(e1.len() + e2.len());
                    buf.push_str(e3);
                    assert_byte_matchers_equal_regex(&buf);
                }
            }
        }

        // Negative & malformed edge cases
        for neg in ["-1", "-10"] {
            for &op in &ops {
                assert_byte_matchers_equal_regex(&format!("{neg}{op}"));
                assert_byte_matchers_equal_regex(&format!("10M{neg}{op}"));
                assert_byte_matchers_equal_regex(&format!("{neg}{op}10M"));
                assert_byte_matchers_equal_regex(&format!("5S{neg}{op}10M2S"));
            }
        }
        for special in [
            "",
            "-",
            "--5M",
            "-10-5M",
            "10M-5M",
            "10M-1M5S",
            "0S10M",
            "A-5M",
            "5S\n-5M",
            "10Z-5M",
            "10ä-5M",
            "5S\n10M",
            "5S\n10M5S",
            "10٠M",
            "10٠S5M",
            "10²5M",
            "10½5M5S",
        ] {
            assert_byte_matchers_equal_regex(special);
        }

        // Verify `is_regex_decimal_digit` matches `regex::Regex` `\d` across every Unicode scalar value.
        let re_dig = crate::util::static_regex(r"^\d$");
        let mut ch_buf = [0u8; 4];
        for ch in (0..=0x10FFFFu32).filter_map(char::from_u32) {
            let s = ch.encode_utf8(&mut ch_buf);
            assert_eq!(
                is_regex_decimal_digit(ch),
                re_dig.is_match(s),
                "is_regex_decimal_digit mismatch on {ch:?} (U+{:04X})",
                ch as u32
            );
        }

        // 2. Full CigarModifier::modify_cigar DIFF test against modify_cigar_legacy
        let config = get_default_config();
        let ref_hetero = b"ACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTTACGTACGTACGT";
        let ref_homo = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut seed = SeedMap::new();
        // Populate seed with 12-mer at start of rc(query_seq) so chimeric 5'/3' soft-clip triggers
        let sclip_15 = b"ACGTACGTACGTACG";
        let rc_15 = rc_seq(sclip_15);
        seed.insert(rc_15[0..SEED_2].to_vec(), vec![15]);
        seed.insert(rc_15[rc_15.len() - SEED_2..].to_vec(), vec![15]);

        let explicit_cigars = [
            "2I100M",
            "100M2I",
            "2I100M3I",
            "15I50M",
            "50M15I",
            "15S50M",
            "50M15S",
            "10M1I5M1I5M1I10M",
            "10M2D5M3D5M1D10M",
            "10M2D5M3I5M1D10M",
            "5M2D5M3D10M",
            "5M2D5M3D2I10M",
            "5M2I5M3I10M",
            "5M2I5M3D2I10M",
            "3M2I20M",
            "3M2D20M",
            "20M2I3M",
            "20M2D3M",
            "104M1I2M2D1M1D5M1D2M9D27M",
            "10M1D1M1D1M1D1M",
            "10M1I1M1D1M1D1M",
            "10M1D1M1I1M1D1M",
        ];
        for cig in explicit_cigars {
            for reference in [ref_hetero.as_slice(), ref_homo.as_slice()] {
                let qseq = reference[..120].to_vec();
                let qqual = vec![b'I'; qseq.len()];
                assert_modify_cigar_matches_legacy(
                    10, cig, &qseq, &qqual, reference, 1, &seed, 3, 150, &config,
                );
            }
        }

        // 20,000 seeded-random 1..6-op reads
        let mut rng_state: u64 = 0xDEAD_BEEF_CAFE_BABEu64;
        let mut next_u64 = || -> u64 {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            rng_state
        };

        let rand_ops = b"MMMIIDDSNH";
        let rand_lens = [1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 35];
        let bases = b"ACGT";

        let mut cigar_buf = String::with_capacity(64);
        let mut qseq_buf = Vec::with_capacity(160);
        let mut qqual_buf = Vec::with_capacity(160);

        for i in 0..20_000 {
            let num_ops = (next_u64() % 6 + 1) as usize;
            cigar_buf.clear();
            let mut read_len = 0usize;
            for _ in 0..num_ops {
                let op = rand_ops[(next_u64() as usize) % rand_ops.len()] as char;
                let len = rand_lens[(next_u64() as usize) % rand_lens.len()];
                use std::fmt::Write;
                let _ = write!(&mut cigar_buf, "{len}{op}");
                if matches!(op, 'M' | 'I' | 'S') {
                    read_len += len as usize;
                }
            }
            if read_len == 0 {
                read_len = 20;
            }
            let reference = if i % 4 == 0 {
                ref_homo.as_slice()
            } else {
                ref_hetero.as_slice()
            };
            let position = ((next_u64() % 40) + 1) as i32;
            let indel = (next_u64() % 4) as i32;

            qseq_buf.clear();
            qqual_buf.clear();
            for j in 0..read_len {
                let ref_idx = ((position - 1) as usize + j) % reference.len();
                let base = if (next_u64() % 8) == 0 {
                    bases[(next_u64() as usize) % 4]
                } else {
                    reference[ref_idx]
                };
                qseq_buf.push(base);
                let qual = if (next_u64() % 6) == 0 {
                    b'#' // low qual
                } else {
                    b'I' // high qual
                };
                qqual_buf.push(qual);
            }
            // Occasionally inject the exact 15-bp chimeric soft-clip sequence at 5' or 3'
            if read_len >= 15 && (i % 10 == 0) {
                qseq_buf[..15].copy_from_slice(sclip_15);
            } else if read_len >= 15 && (i % 10 == 5) {
                let start = read_len - 15;
                qseq_buf[start..].copy_from_slice(sclip_15);
            }

            assert_modify_cigar_matches_legacy(
                position, &cigar_buf, &qseq_buf, &qqual_buf, reference, 1, &seed, indel, 150,
                &config,
            );
        }
    }

    #[test]
    fn modify_cigar_runs_past_100_merges() {
        let mut ref_seq = Vec::with_capacity(4500);
        let bases = b"ACGT";
        let mut state: u64 = 0xCAFE_BABE_1234_5678;
        for i in 0..4500 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let mut b = bases[(state as usize) % 4];
            if i > 0 && ref_seq[i - 1] == b {
                b = bases[((state as usize) + 1) % 4];
            }
            ref_seq.push(b);
        }
        let seed = SeedMap::new();
        let config = get_default_config();

        for n in [199usize, 201, 249] {
            let input_cigar = format!("20M{}1D20M", "1D10M".repeat(n));
            let num_d = n + 1;
            let num_pairs = num_d / 2;
            let expected_cigar = format!(
                "20M12D10I{}20M",
                "10M12D10I".repeat(num_pairs.saturating_sub(1))
            );

            // Construct read matching reference at pos 100 with each 1D skipped
            let mut qseq = Vec::with_capacity(40 + n * 10);
            let mut ref_cursor = 99usize; // 0-based for pos 100
            qseq.extend_from_slice(&ref_seq[ref_cursor..ref_cursor + 20]);
            ref_cursor += 20;
            for _ in 0..n {
                ref_cursor += 1; // 1D
                qseq.extend_from_slice(&ref_seq[ref_cursor..ref_cursor + 10]);
                ref_cursor += 10;
            }
            ref_cursor += 1; // trailing 1D
            qseq.extend_from_slice(&ref_seq[ref_cursor..ref_cursor + 20]);
            let qqual = vec![b'I'; qseq.len()];

            let mod_new = CigarModifier::new(
                100,
                input_cigar.clone(),
                qseq.clone(),
                qqual.clone(),
                &ref_seq,
                1,
                &seed,
                num_d as i32,
                qseq.len() as i32,
                &config,
            )
            .modify_cigar()
            .unwrap();
            let mod_leg = CigarModifier::new(
                100,
                input_cigar,
                qseq,
                qqual,
                &ref_seq,
                1,
                &seed,
                num_d as i32,
                3000,
                &config,
            )
            .modify_cigar_legacy()
            .unwrap();

            assert_eq!(mod_new.cigar, expected_cigar, "failed for N={n}");
            assert_eq!(mod_leg.cigar, expected_cigar, "legacy failed for N={n}");
        }

        // Negative control: first D-M-D has M=20 (> 15), blocking first-match combine_to_close_to_correct
        let neg_cigar = format!("20M1D20M{}1D20M", "1D10M".repeat(249));
        let mut qseq_neg = Vec::with_capacity(60 + 249 * 10);
        let mut ref_cursor = 99usize;
        qseq_neg.extend_from_slice(&ref_seq[ref_cursor..ref_cursor + 20]);
        ref_cursor += 21;
        qseq_neg.extend_from_slice(&ref_seq[ref_cursor..ref_cursor + 20]);
        ref_cursor += 20;
        for _ in 0..249 {
            ref_cursor += 1;
            qseq_neg.extend_from_slice(&ref_seq[ref_cursor..ref_cursor + 10]);
            ref_cursor += 10;
        }
        ref_cursor += 1;
        qseq_neg.extend_from_slice(&ref_seq[ref_cursor..ref_cursor + 20]);
        let qqual_neg = vec![b'I'; qseq_neg.len()];
        let mod_neg = CigarModifier::new(
            100,
            neg_cigar.clone(),
            qseq_neg,
            qqual_neg,
            &ref_seq,
            1,
            &seed,
            251,
            3000,
            &config,
        )
        .modify_cigar()
        .unwrap();
        assert_eq!(mod_neg.cigar, neg_cigar);
    }

    fn assert_single_match_fastpath_matches_full(
        position: i32,
        mch: i32,
        query_seq: &[u8],
        reference: &[u8],
        ref_start: i32,
        seed: &SeedMap,
        config: &Config,
    ) {
        use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

        let mut fast_pos = position;
        let mut fast_elements = vec![(mch as u32, CigarKind::Match)];
        let fast_res = modify_single_match_in_place(
            &mut fast_pos,
            mch,
            &mut fast_elements,
            query_seq,
            reference,
            ref_start,
        );

        let cigar_str = format!("{mch}M");
        let query_qual = vec![b'I'; query_seq.len()];
        let full_res = CigarModifier::new(
            position,
            cigar_str,
            query_seq.to_vec(),
            query_qual,
            reference,
            ref_start,
            seed,
            0,
            150,
            config,
        )
        .modify_cigar();

        match (&fast_res, &full_res) {
            (Ok(()), Ok(modified)) => {
                let mut full_elements = Vec::new();
                assert!(
                    crate::cigar::parse_cigar_into(&modified.cigar, &mut full_elements).is_some(),
                    "Failed to parse modified CIGAR {}",
                    modified.cigar
                );
                assert_eq!(
                    fast_pos, modified.position,
                    "position mismatch for pos={position} mch={mch} ref_start={ref_start}"
                );
                assert_eq!(
                    fast_elements, full_elements,
                    "cigar_elements mismatch for pos={position} mch={mch} ref_start={ref_start}"
                );
            }
            (
                Err(crate::types::VardictError::JavaThrows(fast_msg)),
                Err(crate::types::VardictError::JavaThrows(full_msg)),
            ) => {
                assert_eq!(
                    fast_msg, full_msg,
                    "JavaThrows message mismatch for pos={position} mch={mch}"
                );
                assert_eq!(
                    fast_pos, position,
                    "fastpath mutated position on error for pos={position} mch={mch}"
                );
            }
            _ => panic!(
                "result mismatch for pos={position} mch={mch}: fast={fast_res:?} vs full={full_res:?}"
            ),
        }
    }

    #[test]
    fn single_match_fastpath_matches_full_modify_cigar() {
        let config = get_default_config();
        let seed = SeedMap::new();
        let reference = b"ACGTACGTACGTAAACCCGGGTTTACGTACGT";
        let ref_start = 10i32;

        // 1. Sweep mch in 1..=12, 5' mismatches 0..=5, 3' mismatches 0..=5, and positions
        // including left overhang (pos < ref_start), in-bounds, and right overhang.
        for mch in 1..=12i32 {
            for pos in [
                ref_start - 2,
                ref_start,
                ref_start + 5,
                ref_start + (reference.len() as i32) - mch,
                ref_start + (reference.len() as i32) - mch + 2,
            ] {
                for mis5 in 0..=5i32 {
                    for mis3 in 0..=5i32 {
                        let mut qseq = Vec::with_capacity(mch as usize);
                        for i in 0..mch {
                            let rpos = pos + i - ref_start;
                            let rbase = if rpos >= 0 && (rpos as usize) < reference.len() {
                                reference[rpos as usize]
                            } else {
                                b'A'
                            };
                            let is_mis = i < mis5 || i >= (mch - mis3);
                            let qbase = if is_mis {
                                if rbase == b'A' {
                                    b'T'
                                } else {
                                    b'A'
                                }
                            } else {
                                rbase
                            };
                            qseq.push(qbase);
                        }
                        assert_single_match_fastpath_matches_full(
                            pos, mch, &qseq, reference, ref_start, &seed, &config,
                        );
                    }
                }
            }
        }

        // 2. Exhaustive bitmask sweep over all match/mismatch patterns for mch in 1..=10
        // (covers non-contiguous mismatches, mch <= 3 all-mismatch JavaThrows, etc.)
        for mch in 1..=10i32 {
            let pos = ref_start + 2;
            for mask in 0..(1u32 << mch) {
                let mut qseq = Vec::with_capacity(mch as usize);
                for i in 0..mch {
                    let rbase = reference[(pos + i - ref_start) as usize];
                    let is_mis = ((mask >> i) & 1) != 0;
                    qseq.push(if is_mis {
                        if rbase == b'A' {
                            b'T'
                        } else {
                            b'A'
                        }
                    } else {
                        rbase
                    });
                }
                assert_single_match_fastpath_matches_full(
                    pos, mch, &qseq, reference, ref_start, &seed, &config,
                );
            }
        }

        // 3. 'N' bases in query and/or reference, plus truncated query_seq edge cases
        let ref_with_n = b"ACNTACGTACGTANNC";
        for mch in 1..=12i32 {
            let mut qseq_n = ref_with_n[..mch as usize].to_vec();
            qseq_n[0] = b'N';
            if mch >= 2 {
                qseq_n[(mch - 1) as usize] = b'N';
            }
            assert_single_match_fastpath_matches_full(
                1, mch, &qseq_n, ref_with_n, 1, &seed, &config,
            );
            // Truncated query_seq shorter than mch
            if mch >= 2 {
                assert_single_match_fastpath_matches_full(
                    1,
                    mch,
                    &qseq_n[..(mch as usize - 1)],
                    ref_with_n,
                    1,
                    &seed,
                    &config,
                );
            }
        }

        // 4. Exhaustive + targeted parity for `modify_pure_sm_in_place` across `[S, M]`,
        //    `[M, S]`, and `[S, M, S]` against `CigarModifier::modify_cigar`.
        use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
        let check_pure_sm = |pos: i32,
                             elems: &[(u32, CigarKind)],
                             qseq: &[u8],
                             qqual: &[u8],
                             ref_bytes: &[u8],
                             rstart: i32,
                             cfg: &Config|
         -> PureSmOutcome {
            let mut fast_pos = pos;
            let mut fast_elems = elems.to_vec();
            let fast_res = modify_pure_sm_in_place(
                &mut fast_pos,
                &mut fast_elems,
                qseq,
                qqual,
                ref_bytes,
                rstart,
                &seed,
                150,
                cfg,
            );
            let mut cigar_str = String::new();
            crate::bam::push_cigar_string(&mut cigar_str, elems);
            let full_res = CigarModifier::new(
                pos,
                cigar_str,
                qseq.to_vec(),
                qqual.to_vec(),
                ref_bytes,
                rstart,
                &seed,
                0,
                150,
                cfg,
            )
            .modify_cigar();
            match (fast_res, full_res) {
                (Ok(PureSmOutcome::Fallback), _) => {
                    assert_eq!(fast_pos, pos, "Fallback mutated start_pos");
                    assert_eq!(
                        fast_elems.as_slice(),
                        elems,
                        "Fallback mutated cigar_elements"
                    );
                    PureSmOutcome::Fallback
                }
                (Ok(PureSmOutcome::Handled), Ok(modified)) => {
                    let mut full_elems = Vec::new();
                    assert!(
                        crate::cigar::parse_cigar_into(&modified.cigar, &mut full_elems).is_some(),
                        "Failed to parse modified CIGAR {}",
                        modified.cigar
                    );
                    assert_eq!(
                        fast_pos, modified.position,
                        "position mismatch for elems={elems:?} pos={pos}"
                    );
                    assert_eq!(
                        fast_elems, full_elems,
                        "cigar mismatch for elems={elems:?} pos={pos}: fast={fast_elems:?} vs full={}",
                        modified.cigar
                    );
                    assert_eq!(modified.query_seq.as_slice(), qseq);
                    assert_eq!(modified.query_qual.as_slice(), qqual);
                    PureSmOutcome::Handled
                }
                (
                    Err(crate::types::VardictError::JavaThrows(fast_msg)),
                    Err(crate::types::VardictError::JavaThrows(full_msg)),
                ) => {
                    assert_eq!(fast_msg, full_msg);
                    assert_eq!(fast_pos, pos, "Err mutated start_pos");
                    assert_eq!(fast_elems.as_slice(), elems, "Err mutated cigar_elements");
                    PureSmOutcome::Handled
                }
                (f, m) => panic!("pure S/M parity mismatch for {elems:?} at {pos}: {f:?} vs {m:?}"),
            }
        };

        let chr_end = (ref_start as usize + reference.len() - 1) as u32;
        for shape in [
            vec![(1, CigarKind::SoftClip), (6, CigarKind::Match)],
            vec![(3, CigarKind::SoftClip), (5, CigarKind::Match)],
            vec![(7, CigarKind::SoftClip), (6, CigarKind::Match)],
            vec![(6, CigarKind::Match), (3, CigarKind::SoftClip)],
            vec![(6, CigarKind::Match), (7, CigarKind::SoftClip)],
            vec![
                (2, CigarKind::SoftClip),
                (5, CigarKind::Match),
                (3, CigarKind::SoftClip),
            ],
            vec![
                (6, CigarKind::SoftClip),
                (2, CigarKind::Match),
                (6, CigarKind::SoftClip),
            ],
            vec![(4, CigarKind::SoftClip), (2, CigarKind::Match)],
        ] {
            let s5 = if shape.first().is_some_and(|e| e.1 == CigarKind::SoftClip) {
                shape[0].0 as i32
            } else {
                0
            };
            let qlen: usize = shape.iter().map(|e| e.0 as usize).sum();
            for pos in [
                ref_start,
                ref_start + 3,
                ref_start + 10,
                chr_end as i32 - 5,
                chr_end as i32 - 1,
            ] {
                for qual_pat in [0u8, 1, 2] {
                    for mask in [
                        0u32,
                        0b1,
                        0b11,
                        0b111,
                        0b10101,
                        (1u32 << qlen.min(16)) - 1,
                        0b1100011,
                    ] {
                        let mut qseq = Vec::with_capacity(qlen);
                        let mut qqual = Vec::with_capacity(qlen);
                        for idx in 0..qlen {
                            let r_idx = pos - s5 + (idx as i32) - ref_start;
                            let rbase = if r_idx >= 0 && (r_idx as usize) < reference.len() {
                                reference[r_idx as usize]
                            } else {
                                b'C'
                            };
                            let is_mis = ((mask >> (idx % 16)) & 1) != 0;
                            qseq.push(if is_mis {
                                if rbase == b'A' {
                                    b'T'
                                } else {
                                    b'A'
                                }
                            } else {
                                rbase
                            });
                            let qb = match qual_pat {
                                0 => b'I',
                                1 => {
                                    if idx == 0 || idx + 1 == qlen {
                                        b'#' // Phred 2 <= LOWQUAL (10)
                                    } else {
                                        b'I'
                                    }
                                }
                                _ => b'#',
                            };
                            qqual.push(qb);
                        }
                        let _ = check_pure_sm(
                            pos, &shape, &qseq, &qqual, reference, ref_start, &config,
                        );
                    }
                }
            }
        }

        // Targeted 1: `[M, S]` rescues entire 3' clip (`s3 -> 0`), becoming `[M]`, then trims 5' mismatches (`[S, M]`).
        let pos_ms = ref_start + 4;
        let mut qseq_ms =
            reference[(pos_ms - ref_start) as usize..(pos_ms - ref_start + 9) as usize].to_vec();
        qseq_ms[0] = if qseq_ms[0] == b'A' { b'T' } else { b'A' };
        qseq_ms[1] = if qseq_ms[1] == b'A' { b'T' } else { b'A' };
        let qqual_ms = vec![b'I'; 9];
        assert_eq!(
            check_pure_sm(
                pos_ms,
                &[(6, CigarKind::Match), (3, CigarKind::SoftClip)],
                &qseq_ms,
                &qqual_ms,
                reference,
                ref_start,
                &config,
            ),
            PureSmOutcome::Handled
        );

        // Targeted 2: `final_mch <= 0` (`2M3S` with all mismatches -> `0M5S`? No, `rn < mch` in MS,
        // so test `2S2M` where 3' trims 2 mismatches -> `2S0M2S` and 5' fails to rescue -> JavaThrows!).
        let pos_err = ref_start + 5;
        let mut qseq_err = vec![b'N'; 4]; // mismatches reference everywhere
        let qqual_err = vec![b'I'; 4];
        let mut check_err_pos = pos_err;
        let orig_err_elems = vec![(2, CigarKind::SoftClip), (2, CigarKind::Match)];
        let mut check_err_elems = orig_err_elems.clone();
        assert!(modify_pure_sm_in_place(
            &mut check_err_pos,
            &mut check_err_elems,
            &qseq_err,
            &qqual_err,
            reference,
            ref_start,
            &seed,
            150,
            &config,
        )
        .is_err());
        assert_eq!(check_err_pos, pos_err);
        assert_eq!(check_err_elems, orig_err_elems);
        qseq_err[0] = reference[(pos_err - 2 - ref_start) as usize];
        qseq_err[1] = reference[(pos_err - 1 - ref_start) as usize];
        // Now `2S2M` drops to `2S0M2S` in step 1, then rescues `2S` in step 2 -> `2M2S` (`Handled`, NOT Err!).
        assert_eq!(
            check_pure_sm(
                pos_err,
                &orig_err_elems,
                &qseq_err,
                &qqual_err,
                reference,
                ref_start,
                &config,
            ),
            PureSmOutcome::Handled
        );

        // Targeted 3: Precondition Fallback checks (chimeric clip >= 12 with SEED_2 hit, contig overhang, adaptor, qual mismatch).
        let long_clip_elems = [(12, CigarKind::SoftClip), (10, CigarKind::Match)];
        let long_seq = vec![b'A'; 22];
        let long_qual = vec![b'I'; 22];
        assert_eq!(
            check_pure_sm(
                ref_start + 2,
                &long_clip_elems,
                &long_seq,
                &long_qual,
                reference,
                ref_start,
                &config,
            ),
            PureSmOutcome::Handled
        );
        let mut seed_chim = SeedMap::new();
        seed_chim.insert(vec![b'T'; SEED_2], vec![(ref_start + 4) as u32]);
        let mut chim_pos = ref_start + 2;
        let mut chim_elems = long_clip_elems.to_vec();
        assert_eq!(
            modify_pure_sm_in_place(
                &mut chim_pos,
                &mut chim_elems,
                &long_seq,
                &long_qual,
                reference,
                ref_start,
                &seed_chim,
                150,
                &config,
            )
            .unwrap(),
            PureSmOutcome::Fallback
        );
        let mut chim_cfg = config.clone();
        chim_cfg.chimeric = true;
        assert_eq!(
            modify_pure_sm_in_place(
                &mut chim_pos,
                &mut chim_elems,
                &long_seq,
                &long_qual,
                reference,
                ref_start,
                &seed_chim,
                150,
                &chim_cfg,
            )
            .unwrap(),
            PureSmOutcome::Handled
        );
        let mut adaptor_cfg = config.clone();
        adaptor_cfg.adaptor.push("AGATCGGAAGAG".to_string());
        assert_eq!(
            check_pure_sm(
                ref_start + 2,
                &[(3, CigarKind::SoftClip), (6, CigarKind::Match)],
                &[b'A'; 9],
                &[b'I'; 9],
                reference,
                ref_start,
                &adaptor_cfg,
            ),
            PureSmOutcome::Fallback
        );
        // Contig overhang (`start_pos + m_len - 1 > ref_end_i64`) and truncated quality (`len != seq.len()`)
        assert_eq!(
            check_pure_sm(
                ref_start + 10,
                &[(2, CigarKind::SoftClip), (6, CigarKind::Match)],
                &[b'A'; 8],
                &[b'I'; 8],
                &reference[..13],
                ref_start,
                &config,
            ),
            PureSmOutcome::Fallback
        );
        assert_eq!(
            check_pure_sm(
                ref_start + 2,
                &[(2, CigarKind::SoftClip), (6, CigarKind::Match)],
                &[b'A'; 8],
                &[b'I'; 5],
                reference,
                ref_start,
                &config,
            ),
            PureSmOutcome::Fallback
        );
    }

    #[test]
    fn byte_set_matches_hash_set_across_all_byte_values() {
        let mut bs = ByteSet::new();
        let mut hs = std::collections::HashSet::<u8>::new();
        assert_eq!(bs.len(), 0);
        assert_eq!(bs.len(), hs.len());
        assert_eq!(ByteSet::default().len(), 0);

        // Single-byte and repeated-byte insertions across all 4 u64 words
        for &b in &[
            0u8, 63, 64, 127, 128, 191, 192, 255, b'A', b'C', b'G', b'T', b'N',
        ] {
            bs.insert(b);
            hs.insert(b);
            assert_eq!(bs.len(), hs.len(), "mismatch after inserting {b}");
            // Repeat insertion must be idempotent
            bs.insert(b);
            hs.insert(b);
            assert_eq!(bs.len(), hs.len(), "mismatch after repeating {b}");
        }

        // Full 0u8..=255u8 sweep from fresh sets
        let mut full_bs = ByteSet::new();
        let mut full_hs = std::collections::HashSet::<u8>::new();
        for b in 0u8..=255u8 {
            assert_eq!(full_bs.len(), full_hs.len());
            full_bs.insert(b);
            full_hs.insert(b);
            assert_eq!(full_bs.len(), full_hs.len());
            full_bs.insert(b);
            assert_eq!(full_bs.len(), full_hs.len());
        }
        assert_eq!(full_bs.len(), 256);
    }

    #[test]
    fn chimeric_25bp_asymmetric_soft_clips_and_sum_mnd_mis_match_legacy() {
        let config = get_default_config();
        let reference = b"ACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTT";

        // 1. Asymmetric 25-bp 5' soft-clip (> SEED_2 = 12):
        // rc_seq(sseq)[0..12] corresponds to reverse_complement(sseq[25 - 12..25]), NOT sseq[0..12].
        let sclip5 = b"AAACCCGGGTTTGCTAGCTAGCTAG"; // 25 bp
        assert_eq!(sclip5.len(), 25);
        let rc5_full = rc_seq(sclip5);
        assert_ne!(&rc5_full[..SEED_2], &rc5_full[25 - SEED_2..]);

        let mut qseq5 = Vec::with_capacity(75);
        qseq5.extend_from_slice(sclip5);
        qseq5.extend_from_slice(&reference[9..59]);
        let qqual5 = vec![b'I'; qseq5.len()];

        // Positive match: seed holds rc_seq(sseq)[0..12] == rc(sseq[13..25])
        let mut seed5_hit = SeedMap::new();
        seed5_hit.insert(rc5_full[..SEED_2].to_vec(), vec![20]);
        assert_modify_cigar_matches_legacy(
            10, "25S50M", &qseq5, &qqual5, reference, 1, &seed5_hit, 0, 150, &config,
        );
        let res5 = CigarModifier::new(
            10,
            "25S50M",
            qseq5.clone(),
            qqual5.clone(),
            reference,
            1,
            &seed5_hit,
            0,
            150,
            &config,
        )
        .modify_cigar()
        .unwrap();
        assert_eq!(res5.cigar, "50M");
        assert_eq!(res5.query_seq.len(), 50);

        // Negative control: seed holds the WRONG end rc(sseq[0..12]) == rc5_full[13..25]
        let mut seed5_miss = SeedMap::new();
        seed5_miss.insert(rc5_full[25 - SEED_2..].to_vec(), vec![20]);
        assert_modify_cigar_matches_legacy(
            10,
            "25S50M",
            &qseq5,
            &qqual5,
            reference,
            1,
            &seed5_miss,
            0,
            150,
            &config,
        );
        let res5_miss = CigarModifier::new(
            10,
            "25S50M",
            qseq5.clone(),
            qqual5.clone(),
            reference,
            1,
            &seed5_miss,
            0,
            150,
            &config,
        )
        .modify_cigar()
        .unwrap();
        assert!(res5_miss.cigar.ends_with("M") && res5_miss.query_seq.len() == 75);

        // 2. Asymmetric 25-bp 3' soft-clip (> SEED_2 = 12):
        // rc_seq(sseq)[len - 12..] corresponds to reverse_complement(sseq[0..12]), NOT sseq[13..25].
        let sclip3 = b"TTTGGGCCCAAACGATCGATCGATC"; // 25 bp
        assert_eq!(sclip3.len(), 25);
        let rc3_full = rc_seq(sclip3);
        assert_ne!(&rc3_full[..SEED_2], &rc3_full[25 - SEED_2..]);

        let mut qseq3 = Vec::with_capacity(75);
        qseq3.extend_from_slice(&reference[9..59]);
        qseq3.extend_from_slice(sclip3);
        let qqual3 = vec![b'I'; qseq3.len()];

        // Positive match: seed holds rc_seq(sseq)[13..25] == rc(sseq[0..12])
        let mut seed3_hit = SeedMap::new();
        seed3_hit.insert(rc3_full[25 - SEED_2..].to_vec(), vec![20]);
        assert_modify_cigar_matches_legacy(
            10, "50M25S", &qseq3, &qqual3, reference, 1, &seed3_hit, 0, 150, &config,
        );
        let res3 = CigarModifier::new(
            10,
            "50M25S",
            qseq3.clone(),
            qqual3.clone(),
            reference,
            1,
            &seed3_hit,
            0,
            150,
            &config,
        )
        .modify_cigar()
        .unwrap();
        assert_eq!(res3.cigar, "50M");
        assert_eq!(res3.query_seq.len(), 50);

        // Negative control: seed holds the WRONG end rc(sseq[13..25]) == rc3_full[0..12]
        let mut seed3_miss = SeedMap::new();
        seed3_miss.insert(rc3_full[..SEED_2].to_vec(), vec![20]);
        assert_modify_cigar_matches_legacy(
            10,
            "50M25S",
            &qseq3,
            &qqual3,
            reference,
            1,
            &seed3_miss,
            0,
            150,
            &config,
        );
        let res3_miss = CigarModifier::new(
            10,
            "50M25S",
            qseq3.clone(),
            qqual3.clone(),
            reference,
            1,
            &seed3_miss,
            0,
            150,
            &config,
        )
        .modify_cigar()
        .unwrap();
        assert_eq!(res3_miss.query_seq.len(), 75);

        // Extreme position near i32::MIN must not overflow (this.position - ref_pos).abs()
        let _ = CigarModifier::new(
            i32::MIN,
            "25S50M",
            qseq5.clone(),
            qqual5.clone(),
            reference,
            1,
            &seed5_hit,
            0,
            150,
            &config,
        )
        .modify_cigar();
        let _ = CigarModifier::new(
            i32::MIN,
            "50M25S",
            qseq3,
            qqual3,
            reference,
            1,
            &seed3_hit,
            0,
            150,
            &config,
        )
        .modify_cigar();

        // 3. sum_mnd_and_mis_prefix matches (ALIGNED_LENGTH_MND, SOFT_CLIPPED)
        for sample in ["12S45M3I20M5D100N", "15I30M", "-1M10S", ""] {
            let expected = (
                global_find_sum(&ALIGNED_LENGTH_MND, sample),
                global_find_sum(&SOFT_CLIPPED, sample),
            );
            assert_eq!(
                sum_mnd_and_mis_prefix(sample),
                expected,
                "sum_mnd_and_mis_prefix mismatch for {sample:?}"
            );
        }

        // 4. modify_cigar_preserving_buffers + into_buffers retains capacities on Err
        let mut err_mod = CigarModifier::new(
            1,
            String::with_capacity(64),
            Vec::with_capacity(256),
            Vec::with_capacity(256),
            reference,
            1,
            &seed5_hit,
            0,
            150,
            &config,
        );
        err_mod.cigar_str.push_str("0M");
        err_mod.query_seq.extend_from_slice(b"ACGT");
        err_mod.query_qual.extend_from_slice(b"IIII");
        assert!(err_mod.modify_cigar_preserving_buffers().is_err());
        let (c_buf, s_buf, q_buf) = err_mod.into_buffers();
        assert!(c_buf.capacity() >= 64);
        assert!(s_buf.capacity() >= 256);
        assert!(q_buf.capacity() >= 256);

        // 5. Huge multi-insertion / multi-deletion CIGAR lengths inside indel simplification
        //    loop (`indel > 0`) saturate at i32::MAX without overflowing or panicking.
        let empty_seed = SeedMap::new();
        for (pos, hostile_cigar, expected_sub) in [
            (10, "10M1500000000I1500000000I10M", "2147483647I"),
            (10, "10M1500000000D1500000000D10M", "2147483647D"),
            (10, "1500000000S1500000000I10M", "2147483647S"),
            (10, "10M1500000000I1500000000S", "2147483647S"),
            (10, "1500000000S5M1500000000I10M", "2147483647S"),
            (10, "10M1500000000I5M1500000000S", "2147483647S"),
            (10, "10M2147483647I5M", "2147483647S"),
            (10, "5M2147483647I10M", "2147483647S"),
            (10, "10M1500000000I5M1500000000I10M", "2147483647I"),
            (1_500_000_000, "1500000000S1500000000D10M", "1500000000S"),
        ] {
            let got = CigarModifier::new(
                pos,
                hostile_cigar,
                [b'A'; 20],
                [b'I'; 20],
                b"",
                1,
                &empty_seed,
                1,
                150,
                &config,
            )
            .modify_cigar()
            .unwrap_or_else(|e| panic!("unexpected Err on {hostile_cigar}: {e:?}"));
            assert!(
                got.cigar.contains(expected_sub),
                "expected {expected_sub:?} in modified CIGAR for {hostile_cigar}, got {:?}",
                got.cigar
            );
            // Also verify non-empty reference never panics (returns Ok or JavaThrows cleanly).
            let _ = CigarModifier::new(
                pos,
                hostile_cigar,
                &reference[..20],
                [b'I'; 20],
                reference,
                1,
                &empty_seed,
                1,
                150,
                &config,
            )
            .modify_cigar();
        }
    }

    #[test]
    fn test_modify_pure_sm_chimeric_2digit_precedence_and_100s_parity() {
        use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

        let config = get_default_config();
        let mut reference = Vec::with_capacity(400);
        for idx in 0..400usize {
            reference.push(b"ACGTACGTGCTACGTAAACCGGTT"[(idx * 7 + idx / 3) % 24]);
        }
        let ref_start = 1i32;
        let max_read_length = 150i32;

        let compare_pure_sm_and_legacy = |pos: i32,
                                          elems: &[(u32, CigarKind)],
                                          qseq: &[u8],
                                          qqual: &[u8],
                                          seed: &SeedMap,
                                          cfg: &Config|
         -> (PureSmOutcome, ModifiedCigar) {
            let mut fast_pos = pos;
            let mut fast_elems = elems.to_vec();
            let fast_res = modify_pure_sm_in_place(
                &mut fast_pos,
                &mut fast_elems,
                qseq,
                qqual,
                &reference,
                ref_start,
                seed,
                max_read_length,
                cfg,
            )
            .expect("pure S/M must not error on valid read");

            let mut cigar_str = String::new();
            crate::bam::push_cigar_string(&mut cigar_str, elems);
            let full = CigarModifier::new(
                pos,
                cigar_str.clone(),
                qseq.to_vec(),
                qqual.to_vec(),
                &reference,
                ref_start,
                seed,
                0,
                max_read_length,
                cfg,
            )
            .modify_cigar()
            .expect("modify_cigar must succeed");
            let legacy = CigarModifier::new(
                pos,
                cigar_str,
                qseq.to_vec(),
                qqual.to_vec(),
                &reference,
                ref_start,
                seed,
                0,
                max_read_length,
                cfg,
            )
            .modify_cigar_legacy()
            .expect("modify_cigar_legacy must succeed");

            assert_eq!(full.position, legacy.position);
            assert_eq!(full.cigar, legacy.cigar);
            assert_eq!(full.query_seq, legacy.query_seq);
            assert_eq!(full.query_qual, legacy.query_qual);

            match fast_res {
                PureSmOutcome::Handled => {
                    let mut fast_cigar = String::new();
                    crate::bam::push_cigar_string(&mut fast_cigar, &fast_elems);
                    assert_eq!(fast_pos, legacy.position);
                    assert_eq!(fast_cigar, legacy.cigar);
                    assert_eq!(qseq, legacy.query_seq.as_slice());
                    assert_eq!(qqual, legacy.query_qual.as_slice());
                }
                PureSmOutcome::Fallback => {
                    assert_eq!(fast_pos, pos);
                    assert_eq!(fast_elems.as_slice(), elems);
                }
            }
            (fast_res, full)
        };

        // Helper to build a read matching `reference` across `mch` at `pos` with mismatching soft-clips.
        let build_sm_read = |pos: i32, s5: usize, mch: usize, s3: usize| -> (Vec<u8>, Vec<u8>) {
            let mut qseq = Vec::with_capacity(s5 + mch + s3);
            let clip_pat = b"TTGGAACCTTGGCCAATTGGCCAA";
            for i in 0..s5 {
                let mut b = clip_pat[i % clip_pat.len()];
                let r_idx = (pos - (s5 as i32) + (i as i32) - ref_start) as isize;
                if r_idx >= 0
                    && (r_idx as usize) < reference.len()
                    && reference[r_idx as usize] == b
                {
                    b = if b == b'A' { b'C' } else { b'A' };
                }
                qseq.push(b);
            }
            let m_start = (pos - ref_start) as usize;
            qseq.extend_from_slice(&reference[m_start..m_start + mch]);
            for i in 0..s3 {
                let mut b = clip_pat[(i + 5) % clip_pat.len()];
                let r_idx = (pos + (mch as i32) + (i as i32) - ref_start) as isize;
                if r_idx >= 0
                    && (r_idx as usize) < reference.len()
                    && reference[r_idx as usize] == b
                {
                    b = if b == b'G' { b'T' } else { b'G' };
                }
                qseq.push(b);
            }
            let qqual = vec![b'I'; qseq.len()];
            (qseq, qqual)
        };

        let empty_seed = SeedMap::new();

        // 1. [15S, 100M], [100M, 15S], [100S, 50M], [50M, 100S] without and with unique SEED_2 hits (`AR-01`).
        for (pos, s5, mch, s3) in [
            (50i32, 15usize, 100usize, 0usize),
            (50, 0, 100, 15),
            (120, 100, 50, 0),
            (50, 0, 50, 100),
        ] {
            let elems: Vec<(u32, CigarKind)> = match (s5 > 0, s3 > 0) {
                (true, false) => vec![
                    (s5 as u32, CigarKind::SoftClip),
                    (mch as u32, CigarKind::Match),
                ],
                (false, true) => vec![
                    (mch as u32, CigarKind::Match),
                    (s3 as u32, CigarKind::SoftClip),
                ],
                _ => unreachable!(),
            };
            let (qseq, qqual) = build_sm_read(pos, s5, mch, s3);

            // (a) No seed hit -> Handled in-place (even for >= 15S and >= 100S!)
            let (out_empty, mod_empty) =
                compare_pure_sm_and_legacy(pos, &elems, &qseq, &qqual, &empty_seed, &config);
            assert_eq!(out_empty, PureSmOutcome::Handled);
            assert_eq!(mod_empty.query_seq.len(), s5 + mch + s3);

            let rc_seed = if s5 > 0 {
                crate::util::reverse_complement_bytes::<SEED_2>(&qseq[(s5 - SEED_2)..s5]).unwrap()
            } else {
                let start3 = qseq.len() - s3;
                crate::util::reverse_complement_bytes::<SEED_2>(&qseq[start3..start3 + SEED_2])
                    .unwrap()
            };

            // (b) Distant unique seed (>= 2 * max_read_length = 300 away) -> Handled in-place
            let mut distant_seed = SeedMap::new();
            distant_seed.insert(rc_seed.to_vec(), vec![(pos + 2 * max_read_length) as u32]);
            let (out_dist, _) =
                compare_pure_sm_and_legacy(pos, &elems, &qseq, &qqual, &distant_seed, &config);
            assert_eq!(out_dist, PureSmOutcome::Handled);

            // (c) Multi-position seed (non-unique) -> Handled in-place
            let mut multi_seed = SeedMap::new();
            multi_seed.insert(rc_seed.to_vec(), vec![pos as u32 + 5, pos as u32 + 20]);
            let (out_multi, _) =
                compare_pure_sm_and_legacy(pos, &elems, &qseq, &qqual, &multi_seed, &config);
            assert_eq!(out_multi, PureSmOutcome::Handled);

            // (d) Unique SEED_2 hit within < 2 * max_read_length -> Fallback & truncates (including 100S! AR-01)
            let mut hit_seed = SeedMap::new();
            hit_seed.insert(rc_seed.to_vec(), vec![(pos + 10) as u32]);
            let (out_hit, mod_hit) =
                compare_pure_sm_and_legacy(pos, &elems, &qseq, &qqual, &hit_seed, &config);
            assert_eq!(out_hit, PureSmOutcome::Fallback);
            assert_eq!(mod_hit.cigar, format!("{mch}M"));
            assert_eq!(mod_hit.query_seq.len(), mch);
        }

        // 2. [10S, 100M, 15S] where 3' 15S has a chimeric SEED_2 hit:
        //    5' `10S` matches `^(\d\d+)S`, entering Java's `if` branch (`10 < 12` no-op) and masking `else if`!
        let pos_sms = 60i32;
        let elems_10s_15s = [
            (10, CigarKind::SoftClip),
            (100, CigarKind::Match),
            (15, CigarKind::SoftClip),
        ];
        let (qseq_10s_15s, qqual_10s_15s) = build_sm_read(pos_sms, 10, 100, 15);
        let start3_sms = qseq_10s_15s.len() - 15;
        let rc_seed_3 = crate::util::reverse_complement_bytes::<SEED_2>(
            &qseq_10s_15s[start3_sms..start3_sms + SEED_2],
        )
        .unwrap();
        let mut seed_3_hit = SeedMap::new();
        seed_3_hit.insert(rc_seed_3.to_vec(), vec![(pos_sms + 15) as u32]);

        let (out_masked, mod_masked) = compare_pure_sm_and_legacy(
            pos_sms,
            &elems_10s_15s,
            &qseq_10s_15s,
            &qqual_10s_15s,
            &seed_3_hit,
            &config,
        );
        assert_eq!(out_masked, PureSmOutcome::Handled);
        assert_eq!(mod_masked.cigar, "10S100M15S");
        assert_eq!(mod_masked.query_seq.len(), 125);

        // Contrast: [9S, 100M, 15S] (1-digit 5' clip) DOES enter `else if` and truncate 3' 15S!
        let elems_9s_15s = [
            (9, CigarKind::SoftClip),
            (100, CigarKind::Match),
            (15, CigarKind::SoftClip),
        ];
        let (qseq_9s_15s, qqual_9s_15s) = build_sm_read(pos_sms, 9, 100, 15);
        let start3_9 = qseq_9s_15s.len() - 15;
        let rc_seed_3_9 = crate::util::reverse_complement_bytes::<SEED_2>(
            &qseq_9s_15s[start3_9..start3_9 + SEED_2],
        )
        .unwrap();
        let mut seed_3_9_hit = SeedMap::new();
        seed_3_9_hit.insert(rc_seed_3_9.to_vec(), vec![(pos_sms + 15) as u32]);
        let (out_9s, mod_9s) = compare_pure_sm_and_legacy(
            pos_sms,
            &elems_9s_15s,
            &qseq_9s_15s,
            &qqual_9s_15s,
            &seed_3_9_hit,
            &config,
        );
        assert_eq!(out_9s, PureSmOutcome::Fallback);
        assert_eq!(mod_9s.cigar, "9S100M");
        assert_eq!(mod_9s.query_seq.len(), 109);

        // 3. Multi-indel rewrites, mid = 20 > 15, 3-indel ov5-prefixed CIGARs, endpoint helpers, and non-ASCII `\D` (`AR-07`).
        let q_long = reference[49..250].to_vec();
        let qual_long = vec![b'I'; q_long.len()];
        for cigar_case in [
            "10M1I5M1D20M",
            "10M1D5M1I20M",
            "10M1D5M1D2I20M",
            "10M1I5M1D2I20M",
            "10M1I20M1D10M",
            "10M1D20M1I10M",
            "5S10M1D4M1I4M1D15M",
            "5S10M1D10M1I10M1D15M",
            "10M1D10M1I10M1D5M1D5M1D10M",
            "5S10M1D4M1D4M1D15M",
            "5S10M1I4M1I4M1D15M",
            "3D20M2D",
            "3I20M2I",
            "5S2I20M",
            "20M2D5S",
            "5S6M2I20M",
            "20M2I6M5S",
            "4M2I20M",
            "20M2I4M",
            "12M5D20M",
            "9M5D20M",
            "104M5D12M",
            "104M5D2M",
            "25D",
            "25I",
            "5S\n10M1D4M1I4M1D15M",
            "5S\n10M1D4M1I4M1D15M5S",
            "10²1I5M1D20M",
            "10½1I5M1D20M",
            "10٠1I5M1D20M",
            "5S10٠1I5M1D20M",
            "10\u{10D40}1I5M1D20M",
            "A-5M",
            "5S\n-5M",
            "10Z-5M",
            "10ä-5M",
            "10H1I5M1D20M",
            "10D1I5M1D20M",
            "10H1I5M1I4M1D20M",
            "10D1I5M1I4M1D20M",
            "10M2D5M3D4I20M",
            "10M2D5M3I4I20M",
            "10M2I5M3I4I20M",
            "10M2D3D10M",
            "10M2I3I10M",
            "10ä1I5M1D20M",
            "10€1I5M1I20M",
            "10🧬1I5M1D2I20M",
            "10ä1I20M1D10M",
        ] {
            let got = CigarModifier::new(
                50,
                cigar_case,
                q_long.clone(),
                qual_long.clone(),
                &reference,
                ref_start,
                &empty_seed,
                3,
                max_read_length,
                &config,
            )
            .modify_cigar();
            let want = CigarModifier::new(
                50,
                cigar_case,
                q_long.clone(),
                qual_long.clone(),
                &reference,
                ref_start,
                &empty_seed,
                3,
                max_read_length,
                &config,
            )
            .modify_cigar_legacy();

            match (got, want) {
                (Ok(g), Ok(w)) => {
                    assert_eq!(g.position, w.position, "pos mismatch on {cigar_case}");
                    assert_eq!(g.cigar, w.cigar, "cigar mismatch on {cigar_case}");
                    assert_eq!(g.query_seq, w.query_seq, "seq mismatch on {cigar_case}");
                    assert_eq!(g.query_qual, w.query_qual, "qual mismatch on {cigar_case}");
                }
                (
                    Err(crate::types::VardictError::JavaThrows(gm)),
                    Err(crate::types::VardictError::JavaThrows(wm)),
                ) => {
                    assert_eq!(gm, wm, "err mismatch on {cigar_case}");
                }
                (g, w) => panic!("result mismatch on {cigar_case}: {g:?} vs {w:?}"),
            }
        }

        // 4. O(N) maximal-digit-run guard (`AR-05`): 10,000-digit runs must match legacy in <= 64 scan_digits_fwd calls.
        for long_cigar in [
            format!("{}D1I", "9".repeat(10_000)),
            format!("5M{}N1D20M", "9".repeat(10_000)),
        ] {
            let want_long = CigarModifier::new(
                50,
                long_cigar.as_str(),
                q_long.clone(),
                qual_long.clone(),
                &reference,
                ref_start,
                &empty_seed,
                3,
                max_read_length,
                &config,
            )
            .modify_cigar_legacy();

            let cm_timed = CigarModifier::new(
                50,
                long_cigar.clone(),
                q_long.clone(),
                qual_long.clone(),
                &reference,
                ref_start,
                &empty_seed,
                3,
                max_read_length,
                &config,
            );
            super::reset_scan_digits_fwd_calls();
            let t0 = std::time::Instant::now();
            let got_long = cm_timed.modify_cigar();
            let elapsed = t0.elapsed();
            let calls = super::scan_digits_fwd_calls();
            match (got_long, want_long) {
                (Ok(g), Ok(w)) => {
                    assert_eq!(g.position, w.position);
                    assert_eq!(g.cigar, w.cigar);
                    assert_eq!(g.query_seq, w.query_seq);
                    assert_eq!(g.query_qual, w.query_qual);
                }
                (
                    Err(crate::types::VardictError::JavaThrows(gm)),
                    Err(crate::types::VardictError::JavaThrows(wm)),
                ) => {
                    assert_eq!(gm, wm);
                }
                (g, w) => panic!("long_cigar result mismatch: {g:?} vs {w:?}"),
            }
            assert!(
                calls <= 64,
                "10,000-digit CIGAR must scan maximal digit runs in <= 64 scan_digits_fwd calls, got {calls}"
            );
            assert!(
                elapsed < std::time::Duration::from_secs(10),
                "10,000-digit CIGAR hung: {elapsed:?}"
            );
        }
    }
}
