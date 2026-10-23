//! CIGAR modify — ports CigarModifier.java (ordered fixpoint, bounds-checked)
//!
//! normalizes indels near soft-clips, combines adjacent small indels into
//! complex variants, and cleans up end-of-read alignment mismatches.

use crate::cigar::ModifiedCigar;
use crate::config::Config;
use crate::reference::SeedMap;
use crate::util::{char_at, reverse_complement, subseq};
use regex::Regex;
use std::collections::HashSet;
use std::sync::LazyLock;

// Define all the regexes from Patterns.java
static BEGIN_NUMBER_D: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\d+)D").unwrap());
static END_NUMBER_D: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(\d+)D$").unwrap());
static BEGIN_NUMBER_I: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\d+)I").unwrap());
static END_NUMBER_I: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(\d+)I$").unwrap());

static SA_CIGAR_D_S_5CLIP_GROUP: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(\d\d+)S").unwrap());
static SA_CIGAR_D_S_3CLIP_GROUP: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\d\d+)S$").unwrap());

static BEGIN_NUMBER_S_NUMBER_IOR_D: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(\d+)S(\d+)([ID])").unwrap());
static NUMBER_IOR_D_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\d+)([ID])(\d+)S$").unwrap());
static BEGIN_NUMBER_S_NUMBER_M_NUMBER_IOR_D: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(\d+)S(\d+)M(\d+)([ID])").unwrap());
static NUMBER_IOR_D_NUMBER_M_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\d+)([ID])(\d+)M(\d+)S$").unwrap());

static BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(\d)M(\d+)([ID])(\d+)M").unwrap());
static BEGIN_DIGIT_M_NUMBER_IOR_D_NUMBER_M_OPT: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^\dM\d+[ID]\d+M").unwrap());
static NUMBER_IOR_D_DIGIT_M_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\d+)([ID])(\d)M$").unwrap());
static NUMBER_IOR_D_NUMBER_M_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\d+)([ID])(\d+)M$").unwrap());

static D_M_D_DD_M_D_I_D_M_D_DD: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(.*?)(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M").unwrap());
static D_M_D_DD_M_D_I_D_M_D_DD_PRIM: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M").unwrap());

static THREE_DELETIONS_PATTERN: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(.*?)(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M").unwrap());
static THREE_INDELS_PATTERN: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^(.*?)(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M").unwrap()
});

static DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\d+M\d+[DI]\d+M\d+[DI]\d+M\d+[DI]\d+M").unwrap());
static DM_DD_DM_DD_DM_DD_DM: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\d+M\d+D\d+M\d+D\d+M\d+D\d+M").unwrap());

static DIG_D_DIG_M_DIG_DI_DIGI: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\d+)D(\d+)M(\d+)([DI])(\d+I)?").unwrap());
static NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?").unwrap());

static DIG_D_DIG_D: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(\d+)D(\d+)D").unwrap());
static DIG_I_DIG_I: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(\d+)I(\d+)I").unwrap());

#[cfg(test)]
static ANY_NUMBER_M_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(.*?)(\d+)M(\d+)S$").unwrap());
#[cfg(test)]
static BEGIN_ANY_DIG_M_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(.*?)(\d+)M$").unwrap());
static DIG_M_END: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d+M$").unwrap());
static BEGIN_DIG_M: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\d+)M").unwrap());
static DIG_S_DIG_M: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\d+)S(\d+)M").unwrap());
static DIG_M_DIG_S_END: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d+M\d+S$").unwrap());

static ALIGNED_LENGTH_MND: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(\d+)[MND]").unwrap());
static SOFT_CLIPPED: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(\d+)[MIS]").unwrap());
static RE_START_S: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\d+S").unwrap());
static RE_END_2S: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d\d+S$").unwrap());
static RE_S_M: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\d+)S(\d+)M").unwrap());
#[cfg(test)]
static CIGAR_ELEMENT_PATTERN: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\-?\d+)[MIDNSHPX=]").unwrap());

const LOWQUAL: u8 = 10;
const SEED_2: usize = 12;

/// Helper to parse matched group as integer with default 0.
fn to_int(s: &str) -> i32 {
    s.parse::<i32>().unwrap_or(0)
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
    let mut sum = 0;
    for cap in re.captures_iter(s) {
        if let Some(m) = cap.get(1) {
            sum += to_int(m.as_str());
        }
    }
    sum
}

// ---------------------------------------------------------------------------
// Reference helper functions matching VariationUtils
// ---------------------------------------------------------------------------

#[inline]
fn is_has_and_equals_char(ch1: u8, reference: &[u8], index: i32, ref_start: i32) -> bool {
    let local_idx = index - ref_start;
    if local_idx < 0 || local_idx >= reference.len() as i32 {
        false
    } else {
        reference[local_idx as usize] == ch1
    }
}

#[inline]
fn is_has_and_equals_str(
    reference: &[u8],
    index1: i32,
    ref_start: i32,
    str_bytes: &[u8],
    index2: i32,
) -> bool {
    let local_idx = index1 - ref_start;
    if local_idx < 0 || local_idx >= reference.len() as i32 {
        return false;
    }
    if index2 < 0 || index2 >= str_bytes.len() as i32 {
        return false;
    }
    reference[local_idx as usize] == str_bytes[index2 as usize]
}

#[inline]
fn is_has_and_not_equals_char(ch1: u8, reference: &[u8], index: i32, ref_start: i32) -> bool {
    let local_idx = index - ref_start;
    if local_idx < 0 || local_idx >= reference.len() as i32 {
        false
    } else {
        reference[local_idx as usize] != ch1
    }
}

#[inline]
fn is_has_and_not_equals_str(
    reference: &[u8],
    index1: i32,
    ref_start: i32,
    str_bytes: &[u8],
    index2: i32,
) -> bool {
    let local_idx = index1 - ref_start;
    if local_idx < 0 || local_idx >= reference.len() as i32 {
        return false;
    }
    if index2 < 0 || index2 >= str_bytes.len() as i32 {
        return false;
    }
    reference[local_idx as usize] != str_bytes[index2 as usize]
}

// Helper to reverse complement sequence string
fn rc_seq(seq: &[u8]) -> Vec<u8> {
    let mut v = seq.to_vec();
    reverse_complement(&mut v);
    v
}

/// The main CIGAR Modifier engine porting CigarModifier.java.
pub struct CigarModifier<'a> {
    position: i32,
    cigar_str: String,
    query_seq: Vec<u8>,
    query_qual: Vec<u8>,
    reference: &'a [u8],
    ref_start: i32,
    seed: &'a SeedMap,
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
            seed,
            indel,
            max_read_length,
            config,
        }
    }

    /// Perform all CIGAR modifications. No panic allowed!
    /// Takes `self` by value: the modifier is single-use per read, and the
    /// caller reclaims the buffers from the returned [`ModifiedCigar`].
    pub fn modify_cigar(mut self) -> Result<ModifiedCigar, crate::types::VardictError> {
        //flag is set to true if CIGAR string is modified and should be looked at again
        let flag = true;
        let this = &mut self;
        let (mut has_d, mut has_i, mut has_s) = cigar_flags(&this.cigar_str);

        if has_d {
            // 1. Strip leading deletion
            if let Some(caps) = BEGIN_NUMBER_D.captures(&this.cigar_str) {
                let dlen = to_int(caps.get(1).map_or("", |m| m.as_str()));
                this.position += dlen;
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
            let sc5_opt = SA_CIGAR_D_S_5CLIP_GROUP.captures(&this.cigar_str);
            let sc3_opt = SA_CIGAR_D_S_3CLIP_GROUP.captures(&this.cigar_str);

            if let Some(caps) = sc5_opt {
                let cigar_element = to_int(caps.get(1).map_or("", |m| m.as_str())) as usize;
                if !this.config.chimeric && cigar_element >= SEED_2 {
                    if let Some(sseq) = subseq(&this.query_seq, 0, cigar_element) {
                        let sequence = rc_seq(sseq);
                        if let Some(rc_seed) = sequence.get(0..SEED_2) {
                            if let Some([ref_pos_u32]) = this.seed.get(rc_seed) {
                                let ref_pos = *ref_pos_u32 as i32;
                                if (this.position - ref_pos).abs() < 2 * this.max_read_length {
                                    this.cigar_str =
                                        RE_START_S.replace(&this.cigar_str, "").into_owned();
                                    (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                                    if this.query_seq.len() >= cigar_element {
                                        this.query_seq = this
                                            .query_seq
                                            .get(cigar_element..)
                                            .unwrap_or_default()
                                            .to_vec();
                                        this.query_qual = this
                                            .query_qual
                                            .get(cigar_element..)
                                            .unwrap_or_default()
                                            .to_vec();
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
                let cigar_element = to_int(caps.get(1).map_or("", |m| m.as_str())) as usize;
                if !this.config.chimeric && cigar_element >= SEED_2 {
                    let seq_len = this.query_seq.len();
                    if seq_len >= cigar_element {
                        let start_idx = seq_len - cigar_element;
                        if let Some(sseq) = subseq(&this.query_seq, start_idx, cigar_element) {
                            let sequence = rc_seq(sseq);
                            if sequence.len() >= SEED_2 {
                                let start_seed = sequence.len() - SEED_2;
                                if let Some(rc_seed) = sequence.get(start_seed..) {
                                    if let Some([ref_pos_u32]) = this.seed.get(rc_seed) {
                                        let ref_pos = *ref_pos_u32 as i32;
                                        if (this.position - ref_pos).abs()
                                            < 2 * this.max_read_length
                                        {
                                            this.cigar_str =
                                                RE_END_2S.replace(&this.cigar_str, "").into_owned();
                                            (has_d, has_i, has_s) = cigar_flags(&this.cigar_str);
                                            this.query_seq = this
                                                .query_seq
                                                .get(0..start_idx)
                                                .unwrap_or_default()
                                                .to_vec();
                                            this.query_qual = this
                                                .query_qual
                                                .get(0..start_idx)
                                                .unwrap_or_default()
                                                .to_vec();
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
            }
        }

        // 4. Main indel simplification loop
        let mut flag_loop = flag;
        let mut iter_count = 0;
        while flag_loop && this.indel > 0 && (has_d || has_i) {
            iter_count += 1;
            if iter_count > 100 {
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
                    let tslen = format!("{}S", g1 + if op == "I" { g2 } else { 0 });
                    this.position += if op == "D" { g2 } else { 0 };
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
                    let tslen = format!("{}S", g3 + if op == "I" { g1 } else { 0 });
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
                        let tslen = format!("{}S", g1 + tmid + if op == "I" { g3 } else { 0 });
                        this.position += tmid + if op == "D" { g3 } else { 0 };
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
                        let tslen = format!("{}S", g4 + tmid + if op == "I" { g1 } else { 0 });
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
                    let tslen = format!("{}S", tmid + if op == "I" { g1 } else { 0 });
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
                        + to_int(caps.get(2).map_or("", |m| m.as_str()));
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
                        + to_int(caps.get(2).map_or("", |m| m.as_str()));
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
            let ov5 = ov5.to_string();
            this.capture_mis_softly_ms(&ov5, mch, soft);
        } else if let Some((ov5, mch)) = match_begin_any_dig_m_end(&this.cigar_str) {
            let ov5 = ov5.to_string();
            this.capture_mis_softly_3mismatches(&ov5, mch);
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
        validate_cigar_element_lengths(&self.cigar_str)?;

        Ok(ModifiedCigar {
            position: self.position,
            cigar: self.cigar_str,
            query_seq: self.query_seq,
            query_qual: self.query_qual,
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
                        if let Some(positions) = this.seed.get(rc_seed) {
                            if positions.len() == 1 {
                                let ref_pos = positions[0] as i32;
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
                            if let Some(positions) = this.seed.get(rc_seed) {
                                if positions.len() == 1 {
                                    let ref_pos = positions[0] as i32;
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
        }

        // 4. Main indel simplification loop
        let mut flag_loop = flag;
        let mut iter_count = 0;
        while flag_loop && this.indel > 0 {
            iter_count += 1;
            if iter_count > 100 {
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
            this.capture_mis_softly_ms(&ov5, mch, soft);
        } else if let Some(caps) = BEGIN_ANY_DIG_M_END.captures(&this.cigar_str) {
            let ov5 = caps.get(1).map_or("", |m| m.as_str()).to_string();
            let mch = to_int(caps.get(2).unwrap().as_str());
            this.capture_mis_softly_3mismatches(&ov5, mch);
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
        let mut rn = 0;
        let mut rrn = 0;
        let mut rmch = 0;
        while rrn < mch && rn < mch {
            let ref_pos = self.position + rrn;
            if ref_pos - self.ref_start < 0
                || ref_pos - self.ref_start >= self.reference.len() as i32
            {
                break;
            }
            if is_has_and_not_equals_str(
                self.reference,
                ref_pos,
                self.ref_start,
                &self.query_seq,
                rrn,
            ) {
                rn = rrn + 1;
                rmch = 0;
            } else if is_has_and_equals_str(
                self.reference,
                ref_pos,
                self.ref_start,
                &self.query_seq,
                rrn,
            ) {
                rmch += 1;
            }
            rrn += 1;
            if rmch >= 3 {
                break;
            }
        }
        if rn > 0 && rn <= 3 {
            mch -= rn;
            self.cigar_str = BEGIN_DIG_M
                .replace(&self.cigar_str, &format!("{}S{}M", rn, mch))
                .into_owned();
            self.position += rn;
        }
    }

    fn combine_dig_s_dig_m(&mut self, mut soft: i32, mut mch: i32) {
        let mut rn = 0;
        let mut rn_set = HashSet::new();

        while rn < soft
            && is_has_and_equals_str(
                self.reference,
                self.position - rn - 1,
                self.ref_start,
                &self.query_seq,
                soft - rn - 1,
            )
            && char_at(&self.query_qual, (soft - rn - 1) as isize)
                .unwrap_or(0)
                .saturating_sub(33)
                > LOWQUAL
        {
            rn += 1;
        }

        if rn > 0 {
            mch += rn;
            soft -= rn;
            if soft > 0 {
                self.cigar_str = DIG_S_DIG_M
                    .replace(&self.cigar_str, &format!("{}S{}M", soft, mch))
                    .into_owned();
            } else {
                self.cigar_str = DIG_S_DIG_M
                    .replace(&self.cigar_str, &format!("{}M", mch))
                    .into_owned();
            }
            self.position -= rn;
            rn = 0;
        }

        if soft > 0 {
            while rn + 1 < soft
                && is_has_and_equals_str(
                    self.reference,
                    self.position - rn - 2,
                    self.ref_start,
                    &self.query_seq,
                    soft - rn - 2,
                )
                && char_at(&self.query_qual, (soft - rn - 2) as isize)
                    .unwrap_or(0)
                    .saturating_sub(33)
                    > LOWQUAL
            {
                rn += 1;
                let ref_pos = self.position - rn - 2;
                let local_idx = ref_pos - self.ref_start;
                if local_idx >= 0 && local_idx < self.reference.len() as i32 {
                    rn_set.insert(self.reference[local_idx as usize]);
                }
            }
            let rn_nt = rn_set.len();
            if (rn > 4 && rn_nt > 1)
                || is_has_and_equals_str(
                    self.reference,
                    self.position - 1,
                    self.ref_start,
                    &self.query_seq,
                    soft - 1,
                )
            {
                mch += rn + 1;
                soft -= rn + 1;
                if soft > 0 {
                    self.cigar_str = RE_S_M
                        .replace(&self.cigar_str, &format!("{}S{}M", soft, mch))
                        .into_owned();
                } else {
                    self.cigar_str = RE_S_M
                        .replace(&self.cigar_str, &format!("{}M", mch))
                        .into_owned();
                }
                self.position -= rn + 1;
            }
            if rn == 0 {
                let mut rrn = 0;
                let mut rmch = 0;
                while rrn < mch && rn < mch {
                    let ref_pos = self.position + rrn;
                    if ref_pos - self.ref_start < 0
                        || ref_pos - self.ref_start >= self.reference.len() as i32
                    {
                        break;
                    }
                    if is_has_and_not_equals_str(
                        self.reference,
                        ref_pos,
                        self.ref_start,
                        &self.query_seq,
                        soft + rrn,
                    ) {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if is_has_and_equals_str(
                        self.reference,
                        ref_pos,
                        self.ref_start,
                        &self.query_seq,
                        soft + rrn,
                    ) {
                        rmch += 1;
                    }
                    rrn += 1;
                    if rmch >= 3 {
                        break;
                    }
                }
                if rn > 0 && rn < mch {
                    soft += rn;
                    mch -= rn;
                    self.cigar_str = RE_S_M
                        .replace(&self.cigar_str, &format!("{}S{}M", soft, mch))
                        .into_owned();
                    self.position += rn;
                }
            }
        }
    }

    fn capture_mis_softly_3mismatches(&mut self, ov5: &str, mut mch: i32) {
        let mut refoff = self.position + mch;
        let mut rdoff = mch;
        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }
        let mut rn = 0;
        let mut rrn = 0;
        let mut rmch = 0;
        while rrn < mch && rn < mch {
            let ref_pos = refoff - rrn - 1;
            if ref_pos - self.ref_start < 0
                || ref_pos - self.ref_start >= self.reference.len() as i32
            {
                break;
            }
            if rrn < rdoff
                && is_has_and_not_equals_str(
                    self.reference,
                    ref_pos,
                    self.ref_start,
                    &self.query_seq,
                    rdoff - rrn - 1,
                )
            {
                rn = rrn + 1;
                rmch = 0;
            } else if rrn < rdoff
                && is_has_and_equals_str(
                    self.reference,
                    ref_pos,
                    self.ref_start,
                    &self.query_seq,
                    rdoff - rrn - 1,
                )
            {
                rmch += 1;
            }
            rrn += 1;
            if rmch >= 3 {
                break;
            }
        }
        mch -= rn;
        if rn > 0 && rn <= 3 {
            self.cigar_str = DIG_M_END
                .replace(&self.cigar_str, &format!("{}M{}S", mch, rn))
                .into_owned();
        }
    }

    fn capture_mis_softly_ms(&mut self, ov5: &str, mut mch: i32, mut soft: i32) {
        let mut refoff = self.position + mch;
        let mut rdoff = mch;
        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }
        let mut rn = 0;
        let mut rn_set = HashSet::new();

        while rn < soft
            && is_has_and_equals_str(
                self.reference,
                refoff + rn,
                self.ref_start,
                &self.query_seq,
                rdoff + rn,
            )
            && char_at(&self.query_qual, (rdoff + rn) as isize)
                .unwrap_or(0)
                .saturating_sub(33)
                > LOWQUAL
        {
            rn += 1;
        }

        if rn > 0 {
            mch += rn;
            soft -= rn;
            if soft > 0 {
                self.cigar_str = DIG_M_DIG_S_END
                    .replace(&self.cigar_str, &format!("{}M{}S", mch, soft))
                    .into_owned();
            } else {
                self.cigar_str = DIG_M_DIG_S_END
                    .replace(&self.cigar_str, &format!("{}M", mch))
                    .into_owned();
            }
            rn = 0;
        }

        if soft > 0 {
            while rn + 1 < soft
                && is_has_and_equals_str(
                    self.reference,
                    refoff + rn + 1,
                    self.ref_start,
                    &self.query_seq,
                    rdoff + rn + 1,
                )
                && char_at(&self.query_qual, (rdoff + rn + 1) as isize)
                    .unwrap_or(0)
                    .saturating_sub(33)
                    > LOWQUAL
            {
                rn += 1;
                let ref_pos = refoff + rn + 1;
                let local_idx = ref_pos - self.ref_start;
                if local_idx >= 0 && local_idx < self.reference.len() as i32 {
                    rn_set.insert(self.reference[local_idx as usize]);
                }
            }
            let rn_nt = rn_set.len();
            if rn > 4 && rn_nt > 1 {
                mch += rn + 1;
                soft -= rn + 1;
                if soft > 0 {
                    self.cigar_str = DIG_M_DIG_S_END
                        .replace(&self.cigar_str, &format!("{}M{}S", mch, soft))
                        .into_owned();
                } else {
                    self.cigar_str = DIG_M_DIG_S_END
                        .replace(&self.cigar_str, &format!("{}M", mch))
                        .into_owned();
                }
            }
            if rn == 0 {
                let mut rrn = 0;
                let mut rmch = 0;
                while rrn < mch && rn < mch {
                    let ref_pos = refoff - rrn - 1;
                    if ref_pos - self.ref_start < 0
                        || ref_pos - self.ref_start >= self.reference.len() as i32
                    {
                        break;
                    }
                    if rrn < rdoff
                        && is_has_and_not_equals_str(
                            self.reference,
                            ref_pos,
                            self.ref_start,
                            &self.query_seq,
                            rdoff - rrn - 1,
                        )
                    {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if rrn < rdoff
                        && is_has_and_equals_str(
                            self.reference,
                            ref_pos,
                            self.ref_start,
                            &self.query_seq,
                            rdoff - rrn - 1,
                        )
                    {
                        rmch += 1;
                    }
                    rrn += 1;
                    if rmch >= 3 {
                        break;
                    }
                }
                if rn > 0 && rn < mch {
                    soft += rn;
                    mch -= rn;
                    self.cigar_str = DIG_M_DIG_S_END
                        .replace(&self.cigar_str, &format!("{}M{}S", mch, soft))
                        .into_owned();
                }
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
            let mut ilen = g2 + g3;
            if op == "I" {
                ilen += g4;
            } else if op == "D" {
                dlen += g4;
                if let Some(istr_str) = istr {
                    if !istr_str.is_empty() {
                        ilen += to_int(&istr_str[0..istr_str.len() - 1]);
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
            let mut dlen = g1 + g2;
            let mut ilen = g2;
            if op == "I" {
                ilen += g3;
            } else if op == "D" {
                dlen += g3;
                if let Some(istr_str) = istr {
                    if !istr_str.is_empty() {
                        ilen += to_int(&istr_str[0..istr_str.len() - 1]);
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
        let mut tslen = g5 + g8;
        if op1 == "I" {
            tslen += g3;
        }
        if op2 == "I" {
            tslen += g6;
        }
        if op3 == "I" {
            tslen += g9;
        }

        let mut dlen = g5 + g8;
        if op1 == "D" {
            dlen += g3;
        }
        if op2 == "D" {
            dlen += g6;
        }
        if op3 == "D" {
            dlen += g9;
        }

        let mid = g5 + g8;
        let mut refoff = self.position + g2;
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g11;

        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }

        let mut rn = 0;
        while rdoff + rn < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get((rdoff + rn) as usize) {
                if is_has_and_equals_char(q_base, self.reference, refoff + rn, self.ref_start) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr += rn;
        dlen -= rn;
        tslen -= rn;

        let mut new_cigar = format!("{}M", rdoff_corr);
        if tslen <= 0 {
            dlen -= tslen;
            rm += tslen;
            if dlen == 0 {
                rdoff_corr += rm;
                new_cigar = format!("{}M", rdoff_corr);
            } else if dlen < 0 {
                tslen = -dlen;
                rm += dlen;
                if rm < 0 {
                    rdoff_corr += rm;
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
            rm += dlen;
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
        let mut tslen = g4 + g6;
        let mut dlen = g3 + g4 + g5 + g6 + g7;
        let mid = g4 + g6;
        let mut refoff = self.position + g2;
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }
        let mut rn = 0;
        while rdoff + rn < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get((rdoff + rn) as usize) {
                if is_has_and_equals_char(q_base, self.reference, refoff + rn, self.ref_start) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr += rn;
        dlen -= rn;
        tslen -= rn;
        let mut new_cigar = format!("{}M", rdoff_corr);
        if tslen <= 0 {
            dlen -= tslen;
            rm += tslen;
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
        let mut tslen = g4 + g5 + g6;
        let mut dlen = g3 + g4 + g6 + g7;
        let mid = g4 + g6;
        let mut refoff = self.position + g2;
        let mut rdoff = g2;
        let mut rdoff_corr = g2;
        let mut rm = g8;

        if !ov5.is_empty() {
            refoff += global_find_sum(&ALIGNED_LENGTH_MND, ov5);
            rdoff += global_find_sum(&SOFT_CLIPPED, ov5);
        }

        let mut rn = 0;
        while rdoff + rn < self.query_seq.len() as i32 {
            if let Some(&q_base) = self.query_seq.get((rdoff + rn) as usize) {
                if is_has_and_equals_char(q_base, self.reference, refoff + rn, self.ref_start) {
                    rn += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        rdoff_corr += rn;
        dlen -= rn;
        tslen -= rn;
        let mut new_cigar = format!("{}M", rdoff_corr);
        if tslen <= 0 {
            dlen -= tslen;
            rm += tslen;
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

        let mut tslen = tmid + if op == "I" { g2 } else { 0 };
        self.position += tmid + if op == "D" { g2 } else { 0 };

        while tn < mlen {
            if let Some(&q_base) = self.query_seq.get((tslen + tn) as usize) {
                if is_has_and_not_equals_char(
                    q_base,
                    self.reference,
                    self.position + tn,
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
        tslen += tn;
        mlen -= tn;
        self.position += tn;

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

        let rand_ops = [b'M', b'M', b'M', b'I', b'I', b'D', b'D', b'S', b'N', b'H'];
        let rand_lens = [1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 35];
        let bases = [b'A', b'C', b'G', b'T'];

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
}
