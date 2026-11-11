//! CIGAR modify — ports CigarModifier.java (ordered fixpoint, bounds-checked)
//!
//! normalizes indels near soft-clips, combines adjacent small indels into
//! complex variants, and cleans up end-of-read alignment mismatches.

use crate::cigar::ModifiedCigar;
use crate::config::Config;
use crate::reference::SeedMap;
use crate::util::{char_at, reverse_complement, subseq};
use regex::Regex;
use std::sync::LazyLock;

// Define all the regexes from Patterns.java
static BEGIN_NUMBER_D: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)D"));
static END_NUMBER_D: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)D$"));
static BEGIN_NUMBER_I: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)I"));
static END_NUMBER_I: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)I$"));

#[cfg(test)]
static SA_CIGAR_D_S_5CLIP_GROUP: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d\d+)S"));
#[cfg(test)]
static SA_CIGAR_D_S_3CLIP_GROUP: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d\d+)S$"));

static BEGIN_NUMBER_S_NUMBER_IOR_D: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d+)S(\d+)([ID])"));
static NUMBER_IOR_D_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d+)S$"));
static BEGIN_NUMBER_S_NUMBER_M_NUMBER_IOR_D: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d+)S(\d+)M(\d+)([ID])"));
static NUMBER_IOR_D_NUMBER_M_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d+)M(\d+)S$"));

static BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(\d)M(\d+)([ID])(\d+)M"));
static BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M_OPT: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^\dM\d+[ID]\d+M"));
static NUMBER_IOR_D_DIGIT_M_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d)M$"));
static NUMBER_IOR_D_NUMBER_M_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)([ID])(\d+)M$"));

static D_M_D_DD_M_D_I_D_M_D_DD: LazyLock<Regex> = LazyLock::new(|| {
    crate::util::static_regex(r"^(.*?)(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M")
});
static D_M_D_DD_M_D_I_D_M_D_DD_PRIM: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M"));

static THREE_DELETIONS_PATTERN: LazyLock<Regex> = LazyLock::new(|| {
    crate::util::static_regex(r"^(.*?)(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M")
});
static THREE_INDELS_PATTERN: LazyLock<Regex> = LazyLock::new(|| {
    crate::util::static_regex(r"^(.*?)(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M")
});

static DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"\d+M\d+[DI]\d+M\d+[DI]\d+M\d+[DI]\d+M"));
static DM_DD_DM_DD_DM_DD_DM: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"\d+M\d+D\d+M\d+D\d+M\d+D\d+M"));

static DIG_D_DIG_M_DIG_DI_DIGI: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)D(\d+)M(\d+)([DI])(\d+I)?"));
static NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?"));

static DIG_D_DIG_D: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)D(\d+)D"));
static DIG_I_DIG_I: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)I(\d+)I"));

#[cfg(test)]
static ANY_NUMBER_M_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(.*?)(\d+)M(\d+)S$"));
#[cfg(test)]
static BEGIN_ANY_DIG_M_END: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"^(.*?)(\d+)M$"));
static DIG_M_END: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"\d+M$"));
static BEGIN_DIG_M: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)M"));
static DIG_S_DIG_M: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^(\d+)S(\d+)M"));
static DIG_M_DIG_S_END: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"\d+M\d+S$"));

static ALIGNED_LENGTH_MND: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\d+)[MND]"));
static SOFT_CLIPPED: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"(\d+)[MIS]"));
static RE_START_S: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"^\d+S"));
static RE_END_2S: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"\d\d+S$"));
#[cfg(test)]
static CIGAR_ELEMENT_PATTERN: LazyLock<Regex> =
    LazyLock::new(|| crate::util::static_regex(r"(\-?\d+)[MIDNSHPX=]"));

const LOWQUAL: u8 = 10;
const SEED_2: usize = 12;

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
        (self.0[0].count_ones()
            + self.0[1].count_ones()
            + self.0[2].count_ones()
            + self.0[3].count_ones()) as usize
    }
}

/// Helper to parse matched group as integer with default 0.
fn to_int(s: &str) -> i32 {
    s.parse::<i32>().unwrap_or(0)
}

/// Matches `SA_CIGAR_D_S_5CLIP_GROUP` (`^(\d\d+)S`), returning `to_int(group1) as usize`.
#[inline]
pub(crate) fn leading_2digit_s_len(s: &str) -> Option<usize> {
    let bytes = s.as_bytes();
    let mut s_pos = 0;
    while bytes.get(s_pos).is_some_and(u8::is_ascii_digit) {
        s_pos += 1;
    }
    if s_pos < 2 || bytes.get(s_pos).copied() != Some(b'S') {
        return None;
    }
    Some(to_int(s.get(..s_pos)?) as usize)
}

/// Matches `SA_CIGAR_D_S_3CLIP_GROUP` (`(\d\d+)S$`), returning `to_int(group1) as usize`.
#[inline]
pub(crate) fn trailing_2digit_s_len(s: &str) -> Option<usize> {
    let bytes = s.as_bytes();
    let s_pos = bytes.len().checked_sub(1)?;
    if *bytes.get(s_pos)? != b'S' {
        return None;
    }
    let mut d_start = s_pos;
    while d_start > 0 && bytes.get(d_start - 1).is_some_and(u8::is_ascii_digit) {
        d_start -= 1;
    }
    if s_pos - d_start < 2 {
        return None;
    }
    Some(to_int(s.get(d_start..s_pos)?) as usize)
}

/// Single-pass replacement for `(global_find_sum(&ALIGNED_LENGTH_MND, s), global_find_sum(&SOFT_CLIPPED, s))`,
/// matching `(\d+)[MND]` and `(\d+)[MIS]` byte-for-byte without regex or heap allocation.
#[inline]
fn sum_mnd_and_mis_prefix(s: &str) -> (i32, i32) {
    let bytes = s.as_bytes();
    let mut mnd: i32 = 0;
    let mut mis: i32 = 0;
    let mut idx = 0;
    while idx < bytes.len() {
        if bytes.get(idx).is_some_and(u8::is_ascii_digit) {
            let start = idx;
            while bytes.get(idx).is_some_and(u8::is_ascii_digit) {
                idx += 1;
            }
            if let Some(&op) = bytes.get(idx) {
                let is_mnd = matches!(op, b'M' | b'N' | b'D');
                let is_mis = matches!(op, b'M' | b'I' | b'S');
                if is_mnd || is_mis {
                    let val = s.get(start..idx).map_or(0, to_int);
                    if is_mnd {
                        mnd = mnd.saturating_add(val);
                    }
                    if is_mis {
                        mis = mis.saturating_add(val);
                    }
                    idx += 1;
                }
            }
        } else {
            idx += 1;
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
    let bytes = s.as_bytes();
    let s_pos = bytes.len().checked_sub(1)?;
    if *bytes.get(s_pos)? != b'S' {
        return None;
    }
    let mut d2_start = s_pos;
    while d2_start > 0 && bytes.get(d2_start - 1).is_some_and(|b| b.is_ascii_digit()) {
        d2_start -= 1;
    }
    if d2_start == s_pos {
        return None;
    }
    let m_pos = d2_start.checked_sub(1)?;
    if *bytes.get(m_pos)? != b'M' {
        return None;
    }
    let mut d1_start = m_pos;
    while d1_start > 0 && bytes.get(d1_start - 1).is_some_and(|b| b.is_ascii_digit()) {
        d1_start -= 1;
    }
    if d1_start == m_pos {
        return None;
    }
    let prefix = s.get(..d1_start)?;
    let d1 = to_int(s.get(d1_start..m_pos)?);
    let d2 = to_int(s.get(d2_start..s_pos)?);
    Some((prefix, d1, d2))
}

/// Matches `BEGIN_ANY_DIG_M_END` (`^(.*?)(\d+)M$`) without regex or heap allocation.
#[inline]
fn match_begin_any_dig_m_end(s: &str) -> Option<(&str, i32)> {
    let bytes = s.as_bytes();
    let m_pos = bytes.len().checked_sub(1)?;
    if *bytes.get(m_pos)? != b'M' {
        return None;
    }
    let mut d1_start = m_pos;
    while d1_start > 0 && bytes.get(d1_start - 1).is_some_and(|b| b.is_ascii_digit()) {
        d1_start -= 1;
    }
    if d1_start == m_pos {
        return None;
    }
    let prefix = s.get(..d1_start)?;
    let d1 = to_int(s.get(d1_start..m_pos)?);
    Some((prefix, d1))
}

/// Matches `DIG_S_DIG_M` (`^(\d+)S(\d+)M`) without regex or heap allocation.
#[inline]
fn match_dig_s_dig_m(s: &str) -> Option<(i32, i32)> {
    let bytes = s.as_bytes();
    let mut s_pos = 0;
    while bytes.get(s_pos).is_some_and(|b| b.is_ascii_digit()) {
        s_pos += 1;
    }
    if s_pos == 0 || bytes.get(s_pos).copied() != Some(b'S') {
        return None;
    }
    let d2_start = s_pos + 1;
    let mut m_pos = d2_start;
    while bytes.get(m_pos).is_some_and(|b| b.is_ascii_digit()) {
        m_pos += 1;
    }
    if m_pos == d2_start || bytes.get(m_pos).copied() != Some(b'M') {
        return None;
    }
    let d1 = to_int(s.get(..s_pos)?);
    let d2 = to_int(s.get(d2_start..m_pos)?);
    Some((d1, d2))
}

/// Matches `BEGIN_DIG_M` (`^(\d+)M`) without regex or heap allocation.
#[inline]
fn match_begin_dig_m(s: &str) -> Option<i32> {
    let bytes = s.as_bytes();
    let mut m_pos = 0;
    while bytes.get(m_pos).is_some_and(|b| b.is_ascii_digit()) {
        m_pos += 1;
    }
    if m_pos == 0 || bytes.get(m_pos).copied() != Some(b'M') {
        return None;
    }
    Some(to_int(s.get(..m_pos)?))
}

/// Scans `(\-?\d+)[MIDNSHPX=]` across `s` and validates that all matched lengths are `> 0`.
#[inline]
fn validate_cigar_element_lengths(s: &str) -> Result<(), crate::types::VardictError> {
    let bytes = s.as_bytes();
    let mut idx = 0;
    while idx < bytes.len() {
        let start = idx;
        if bytes.get(idx) == Some(&b'-') && bytes.get(idx + 1).is_some_and(|b| b.is_ascii_digit()) {
            idx += 1;
        }
        let dig_start = idx;
        while bytes.get(idx).is_some_and(|b| b.is_ascii_digit()) {
            idx += 1;
        }
        if idx > dig_start {
            if let Some(&op) = bytes.get(idx) {
                if matches!(
                    op,
                    b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'='
                ) {
                    let num_str = s.get(start..idx).unwrap_or("");
                    let len: i32 = num_str.parse().unwrap_or(0);
                    if len <= 0 {
                        return Err(crate::types::VardictError::JavaThrows(format!(
                            "invalid CIGAR length {} in {}",
                            len, s
                        )));
                    }
                    idx += 1;
                }
            }
        } else {
            idx += 1;
        }
    }
    Ok(())
}

/// Sum all matches of a regex group 1 parsed as integers.
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

/// Zero-copy in-place CIGAR modifier for 2-op (`[S, M]`, `[M, S]`) and 3-op (`[S, M, S]`) reads
/// whose `query_seq` and `query_qual` are guaranteed untouched by `CigarModifier.java:61-243`.
pub(crate) fn modify_pure_sm_in_place(
    position: &mut i32,
    cigar_elements: &mut Vec<(u32, noodles::sam::alignment::record::cigar::op::Kind)>,
    query_seq: &[u8],
    query_qual: &[u8],
    reference: &[u8],
    ref_start: i32,
    chr_len: u32,
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
    if pos0_i64 < ref_start_i64 || m_end_i64 > i64::from(chr_len) || m_end_i64 > ref_end_i64 {
        return Ok(PureSmOutcome::Fallback);
    }

    // 3. No adaptor trimming.
    if !config.adaptor.is_empty() {
        return Ok(PureSmOutcome::Fallback);
    }

    // 4. Reject clips that could trigger chimeric SEED_2 sequence/quality truncation.
    if !config.chimeric && ((s5 as usize) >= SEED_2 || (s3 as usize) >= SEED_2) {
        return Ok(PureSmOutcome::Fallback);
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
            SeedBackend::Map(map) => match map.get(rc_seed) {
                Some([pos]) => Some(*pos),
                _ => None,
            },
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
            if let Some(caps) = BEGIN_NUMBER_D.captures(&this.cigar_str) {
                let dlen = to_int(caps.get(1).map_or("", |m| m.as_str()));
                this.position = this.position.saturating_add(dlen);
                this.cigar_str = BEGIN_NUMBER_D.replace(&this.cigar_str, "").into_owned();
                (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
            }
            // Strip trailing deletion
            if has_d && END_NUMBER_D.is_match(&this.cigar_str) {
                this.cigar_str = END_NUMBER_D.replace(&this.cigar_str, "").into_owned();
                (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
            }
        }

        if has_i {
            // 2. Replace leading insertion with soft clip
            if let Some(caps) = BEGIN_NUMBER_I.captures(&this.cigar_str) {
                let ilen = to_int(caps.get(1).map_or("", |m| m.as_str()));
                this.cigar_str = BEGIN_NUMBER_I
                    .replace(&this.cigar_str, &format!("{}S", ilen))
                    .into_owned();
                (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
            }
            // Replace trailing insertion with soft clip
            if has_i {
                if let Some(caps) = END_NUMBER_I.captures(&this.cigar_str) {
                    let ilen = to_int(caps.get(1).map_or("", |m| m.as_str()));
                    this.cigar_str = END_NUMBER_I
                        .replace(&this.cigar_str, &format!("{}S", ilen))
                        .into_owned();
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
                        .and_then(crate::util::reverse_complement_bytes::<SEED_2>);
                    if let Some(rc_seed) = rc_seed_opt {
                        if let Some(ref_pos_u32) = this.unique_seed_pos(&rc_seed) {
                            let ref_pos = ref_pos_u32 as i32;
                            if (i64::from(this.position) - i64::from(ref_pos)).abs()
                                < 2 * i64::from(this.max_read_length)
                            {
                                let verbose_seq = if this.config.verbose {
                                    subseq(&this.query_seq, 0, cigar_element).map(rc_seq)
                                } else {
                                    None
                                };
                                this.cigar_str =
                                    RE_START_S.replace(&this.cigar_str, "").into_owned();
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
                                        String::from_utf8_lossy(&sequence),
                                        this.position,
                                        SEED_2
                                    );
                                }
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
                            .and_then(crate::util::reverse_complement_bytes::<SEED_2>);
                        if let Some(rc_seed) = rc_seed_opt {
                            if let Some(ref_pos_u32) = this.unique_seed_pos(&rc_seed) {
                                let ref_pos = ref_pos_u32 as i32;
                                if (i64::from(this.position) - i64::from(ref_pos)).abs()
                                    < 2 * i64::from(this.max_read_length)
                                {
                                    let verbose_seq = if this.config.verbose {
                                        subseq(&this.query_seq, start_idx, cigar_element)
                                            .map(rc_seq)
                                    } else {
                                        None
                                    };
                                    this.cigar_str =
                                        RE_END_2S.replace(&this.cigar_str, "").into_owned();
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
        while flag_loop && this.indel > 0 && (has_d || has_i) {
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
            if has_s && (has_d || has_i) {
                if let Some(caps) = BEGIN_NUMBER_S_NUMBER_IOR_D.captures(&this.cigar_str) {
                    let g1 = to_int(caps.get(1).map_or("", |m| m.as_str()));
                    let g2 = to_int(caps.get(2).map_or("", |m| m.as_str()));
                    let op = caps.get(3).map_or("", |m| m.as_str());
                    let tslen = format!("{}S", g1.saturating_add(if op == "I" { g2 } else { 0 }));
                    this.position = this.position.saturating_add(if op == "D" { g2 } else { 0 });
                    this.cigar_str = BEGIN_NUMBER_S_NUMBER_IOR_D
                        .replace(&this.cigar_str, &tslen)
                        .into_owned();
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            // I/D followed by S
            if has_s && (has_d || has_i) {
                if let Some(caps) = NUMBER_IOR_D_NUMBER_S_END.captures(&this.cigar_str) {
                    let g1 = to_int(caps.get(1).map_or("", |m| m.as_str()));
                    let op = caps.get(2).map_or("", |m| m.as_str());
                    let g3 = to_int(caps.get(3).map_or("", |m| m.as_str()));
                    let tslen = format!("{}S", g3.saturating_add(if op == "I" { g1 } else { 0 }));
                    this.cigar_str = NUMBER_IOR_D_NUMBER_S_END
                        .replace(&this.cigar_str, &tslen)
                        .into_owned();
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            // S-M(<=10)-I/D
            if has_s && (has_d || has_i) {
                if let Some(caps) = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IOR_D.captures(&this.cigar_str) {
                    let tmid = to_int(caps.get(2).map_or("", |m| m.as_str()));
                    if tmid <= 10 {
                        let g1 = to_int(caps.get(1).map_or("", |m| m.as_str()));
                        let g3 = to_int(caps.get(3).map_or("", |m| m.as_str()));
                        let op = caps.get(4).map_or("", |m| m.as_str());
                        let tslen = format!(
                            "{}S",
                            g1.saturating_add(tmid)
                                .saturating_add(if op == "I" { g3 } else { 0 })
                        );
                        this.position = this
                            .position
                            .saturating_add(tmid.saturating_add(if op == "D" { g3 } else { 0 }));
                        this.cigar_str = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IOR_D
                            .replace(&this.cigar_str, &tslen)
                            .into_owned();
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                        flag_loop = true;
                    }
                }
            }

            // I/D-M(<=10)-S
            if has_s && (has_d || has_i) {
                if let Some(caps) = NUMBER_IOR_D_NUMBER_M_NUMBER_S_END.captures(&this.cigar_str) {
                    let tmid = to_int(caps.get(3).map_or("", |m| m.as_str()));
                    if tmid <= 10 {
                        let g1 = to_int(caps.get(1).map_or("", |m| m.as_str()));
                        let op = caps.get(2).map_or("", |m| m.as_str());
                        let g4 = to_int(caps.get(4).map_or("", |m| m.as_str()));
                        let tslen = format!(
                            "{}S",
                            g4.saturating_add(tmid)
                                .saturating_add(if op == "I" { g1 } else { 0 })
                        );
                        this.cigar_str = NUMBER_IOR_D_NUMBER_M_NUMBER_S_END
                            .replace(&this.cigar_str, &tslen)
                            .into_owned();
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                        flag_loop = true;
                    }
                }
            }

            // Single-digit M at start with I/D
            if has_d || has_i {
                if let Some(caps) = BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M.captures(&this.cigar_str) {
                    let tmid = to_int(caps.get(1).map_or("", |m| m.as_str()));
                    let mlen = to_int(caps.get(4).map_or("", |m| m.as_str()));
                    let op = caps.get(3).map_or("", |m| m.as_str()).to_string();
                    let g2 = to_int(caps.get(2).map_or("", |m| m.as_str()));
                    flag_loop = this.begin_digit_m_number_ior_d_number_m(tmid, mlen, &op, g2);
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                }
            }

            // Single-digit M at end with I/D
            if has_d || has_i {
                if let Some(caps) = NUMBER_IOR_D_DIGIT_M_END.captures(&this.cigar_str) {
                    let tmid = to_int(caps.get(3).map_or("", |m| m.as_str()));
                    let op = caps.get(2).map_or("", |m| m.as_str());
                    let g1 = to_int(caps.get(1).map_or("", |m| m.as_str()));
                    let tslen = format!("{}S", tmid.saturating_add(if op == "I" { g1 } else { 0 }));
                    this.cigar_str = NUMBER_IOR_D_NUMBER_M_END
                        .replace(&this.cigar_str, &tslen)
                        .into_owned();
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            // Combine close deletions/insertions into complex
            let caps_d_m_d = if has_d && has_i {
                D_M_D_DD_M_D_I_D_M_D_DD.captures(&this.cigar_str)
            } else {
                None
            };
            if let Some(caps) = caps_d_m_d {
                let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
                let g2 = to_int(caps.get(2).map_or("", |m| m.as_str()));
                let g3 = to_int(caps.get(3).map_or("", |m| m.as_str()));
                let g4 = to_int(caps.get(4).map_or("", |m| m.as_str()));
                let g5 = to_int(caps.get(5).map_or("", |m| m.as_str()));
                let g6 = to_int(caps.get(6).map_or("", |m| m.as_str()));
                let g7 = to_int(caps.get(7).map_or("", |m| m.as_str()));
                let g8 = to_int(caps.get(8).map_or("", |m| m.as_str()));
                flag_loop = this.two_deletions_insertion_to_complex(
                    &ov5, g2, g3, g4, g5, g6, g7, g8, flag_loop,
                );
                (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
            } else {
                let caps_three_del = if has_d {
                    THREE_DELETIONS_PATTERN.captures(&this.cigar_str)
                } else {
                    None
                };
                if let Some(caps) = caps_three_del {
                    let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
                    let g2 = to_int(caps.get(2).map_or("", |m| m.as_str()));
                    let g3 = to_int(caps.get(3).map_or("", |m| m.as_str()));
                    let g4 = to_int(caps.get(4).map_or("", |m| m.as_str()));
                    let g5 = to_int(caps.get(5).map_or("", |m| m.as_str()));
                    let g6 = to_int(caps.get(6).map_or("", |m| m.as_str()));
                    let g7 = to_int(caps.get(7).map_or("", |m| m.as_str()));
                    let g8 = to_int(caps.get(8).map_or("", |m| m.as_str()));
                    flag_loop = this.three_deletions(&ov5, g2, g3, g4, g5, g6, g7, g8, flag_loop);
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                } else if has_d || has_i {
                    if let Some(caps) = THREE_INDELS_PATTERN.captures(&this.cigar_str) {
                        let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
                        let g2 = to_int(caps.get(2).map_or("", |m| m.as_str()));
                        let g3 = to_int(caps.get(3).map_or("", |m| m.as_str()));
                        let op1 = caps.get(4).map_or("", |m| m.as_str()).to_string();
                        let g5 = to_int(caps.get(5).map_or("", |m| m.as_str()));
                        let g6 = to_int(caps.get(6).map_or("", |m| m.as_str()));
                        let op2 = caps.get(7).map_or("", |m| m.as_str()).to_string();
                        let g8 = to_int(caps.get(8).map_or("", |m| m.as_str()));
                        let g9 = to_int(caps.get(9).map_or("", |m| m.as_str()));
                        let op3 = caps.get(10).map_or("", |m| m.as_str()).to_string();
                        let g11 = to_int(caps.get(11).map_or("", |m| m.as_str()));
                        flag_loop = this.three_indels(
                            &ov5, g2, g3, &op1, g5, g6, &op2, g8, g9, &op3, g11, flag_loop,
                        );
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    }
                }
            }

            if has_d {
                if let Some(caps) = DIG_D_DIG_M_DIG_DI_DIGI.captures(&this.cigar_str) {
                    let g1 = to_int(caps.get(1).map_or("", |m| m.as_str()));
                    let g2 = to_int(caps.get(2).map_or("", |m| m.as_str()));
                    let g3 = to_int(caps.get(3).map_or("", |m| m.as_str()));
                    let op = caps.get(4).map_or("", |m| m.as_str()).to_string();
                    let istr = caps.get(5).map(|m| m.as_str().to_string());
                    flag_loop = this.combine_to_close_to_correct(
                        g1,
                        g2,
                        g3,
                        &op,
                        istr.as_deref(),
                        flag_loop,
                    );
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                }
            }

            if has_i {
                if let Some(caps) = NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI.captures(&this.cigar_str) {
                    let op_prefix = caps.get(1).map_or("", |m| m.as_str()).to_string();
                    if op_prefix != "D" && op_prefix != "H" {
                        let g2 = to_int(caps.get(2).map_or("", |m| m.as_str()));
                        let g3 = to_int(caps.get(3).map_or("", |m| m.as_str()));
                        let g4 = to_int(caps.get(4).map_or("", |m| m.as_str()));
                        let op = caps.get(5).map_or("", |m| m.as_str()).to_string();
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
                        (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    }
                }
            }

            if has_d {
                if let Some(caps) = DIG_D_DIG_D.captures(&this.cigar_str) {
                    let dlen = to_int(caps.get(1).map_or("", |m| m.as_str()))
                        .saturating_add(to_int(caps.get(2).map_or("", |m| m.as_str())));
                    this.cigar_str = DIG_D_DIG_D
                        .replace(&this.cigar_str, &format!("{}D", dlen))
                        .into_owned();
                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                    flag_loop = true;
                }
            }

            if has_i {
                if let Some(caps) = DIG_I_DIG_I.captures(&this.cigar_str) {
                    let ilen = to_int(caps.get(1).map_or("", |m| m.as_str()))
                        .saturating_add(to_int(caps.get(2).map_or("", |m| m.as_str())));
                    this.cigar_str = DIG_I_DIG_I
                        .replace(&this.cigar_str, &format!("{}I", ilen))
                        .into_owned();
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
            this.combine_dig_s_dig_m(soft, mch);
        } else if let Some(caps) = BEGIN_DIG_M.captures(&this.cigar_str) {
            let mch = to_int(caps.get(1).unwrap().as_str());
            this.combine_begin_dig_m(mch);
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

    fn capture_mis_softly_3mismatches(&mut self, ov5_len: usize, mch: i32) {
        let mut refoff = self.position.saturating_add(mch);
        let mut rdoff = mch;
        if ov5_len > 0 {
            let (mnd, mis) = sum_mnd_and_mis_prefix(self.cigar_str.get(..ov5_len).unwrap_or(""));
            refoff = refoff.saturating_add(mnd);
            rdoff = rdoff.saturating_add(mis);
        }
        self.apply_capture_mis_softly_3mismatches(refoff, rdoff, mch);
    }

    #[cfg(test)]
    fn capture_mis_softly_3mismatches_legacy(&mut self, ov5: &str, mch: i32) {
        let mut refoff = self.position + mch;
        let mut rdoff = mch;
        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }
        self.apply_capture_mis_softly_3mismatches(refoff, rdoff, mch);
    }

    fn apply_capture_mis_softly_3mismatches(&mut self, refoff: i32, rdoff: i32, mut mch: i32) {
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
        self.apply_capture_mis_softly_ms(refoff, rdoff, mch, soft);
    }

    #[cfg(test)]
    fn capture_mis_softly_ms_legacy(&mut self, ov5: &str, mch: i32, soft: i32) {
        let mut refoff = self.position + mch;
        let mut rdoff = mch;
        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }
        self.apply_capture_mis_softly_ms(refoff, rdoff, mch, soft);
    }

    fn apply_capture_mis_softly_ms(&mut self, refoff: i32, rdoff: i32, mch: i32, soft: i32) {
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
        if g3 <= 15 {
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

    fn combine_to_close_to_correct(
        &mut self,
        g1: i32,
        g2: i32,
        g3: i32,
        op: &str,
        istr: Option<&str>,
        mut flag: bool,
    ) -> bool {
        if g2 <= 15 {
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
        mut flag: bool,
    ) -> bool {
        let mut tslen = g5.saturating_add(g8);
        if op1 == "I" {
            tslen = tslen.saturating_add(g3);
        }
        if op2 == "I" {
            tslen = tslen.saturating_add(g6);
        }
        if op3 == "I" {
            tslen = tslen.saturating_add(g9);
        }

        let mut dlen = g5.saturating_add(g8);
        if op1 == "D" {
            dlen = dlen.saturating_add(g3);
        }
        if op2 == "D" {
            dlen = dlen.saturating_add(g6);
        }
        if op3 == "D" {
            dlen = dlen.saturating_add(g9);
        }

        let mid = g5.saturating_add(g8);
        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g11;

        if !ov5.is_empty() {
            refoff = refoff.saturating_add(global_find_sum(&ALIGNED_LENGTH_MND, ov5));
            rdoff = rdoff.saturating_add(global_find_sum(&SOFT_CLIPPED, ov5));
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

        if mid <= 15 {
            self.cigar_str = DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM
                .replace(&self.cigar_str, &new_cigar)
                .into_owned();
            flag = true;
        }
        flag
    }

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
        mut flag: bool,
    ) -> bool {
        let mut tslen = g4.saturating_add(g6);
        let mut dlen = g3
            .saturating_add(g4)
            .saturating_add(g5)
            .saturating_add(g6)
            .saturating_add(g7);
        let mid = g4.saturating_add(g6);
        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if !ov5.is_empty() {
            refoff = refoff.saturating_add(global_find_sum(&ALIGNED_LENGTH_MND, ov5));
            rdoff = rdoff.saturating_add(global_find_sum(&SOFT_CLIPPED, ov5));
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

        if mid <= 15 {
            self.cigar_str = DM_DD_DM_DD_DM_DD_DM
                .replace(&self.cigar_str, &new_cigar)
                .into_owned();
            flag = true;
        }
        flag
    }

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
        mut flag: bool,
    ) -> bool {
        let mut tslen = g4.saturating_add(g5).saturating_add(g6);
        let mut dlen = g3.saturating_add(g4).saturating_add(g6).saturating_add(g7);
        let mid = g4.saturating_add(g6);
        let mut refoff = self.position.saturating_add(g2);
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if !ov5.is_empty() {
            refoff = refoff.saturating_add(global_find_sum(&ALIGNED_LENGTH_MND, ov5));
            rdoff = rdoff.saturating_add(global_find_sum(&SOFT_CLIPPED, ov5));
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

        if mid <= 15 {
            self.cigar_str = D_M_D_DD_M_D_I_D_M_D_DD_PRIM
                .replace(&self.cigar_str, &new_cigar)
                .into_owned();
            flag = true;
        }
        flag
    }

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
        for special in ["", "-", "--5M", "-10-5M", "10M-5M", "10M-1M5S", "0S10M"] {
            assert_byte_matchers_equal_regex(special);
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
                             chr_len: u32,
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
                chr_len,
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
                            pos, &shape, &qseq, &qqual, reference, ref_start, chr_end, &config,
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
                chr_end,
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
            chr_end,
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
                chr_end,
                &config,
            ),
            PureSmOutcome::Handled
        );

        // Targeted 3: Precondition Fallback checks (chimeric clip >= 12, contig overhang, adaptor, qual mismatch).
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
                chr_end,
                &config,
            ),
            PureSmOutcome::Fallback
        );
        let mut chim_cfg = config.clone();
        chim_cfg.chimeric = true;
        assert_eq!(
            check_pure_sm(
                ref_start + 2,
                &long_clip_elems,
                &long_seq,
                &long_qual,
                reference,
                ref_start,
                chr_end,
                &chim_cfg,
            ),
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
                chr_end,
                &adaptor_cfg,
            ),
            PureSmOutcome::Fallback
        );
        // Contig overhang (`start_pos + m_len - 1 > chr_len`) and truncated quality (`len != seq.len()`)
        assert_eq!(
            check_pure_sm(
                ref_start + 10,
                &[(2, CigarKind::SoftClip), (6, CigarKind::Match)],
                &[b'A'; 8],
                &[b'I'; 8],
                reference,
                ref_start,
                (ref_start + 12) as u32,
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
                chr_end,
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
}
