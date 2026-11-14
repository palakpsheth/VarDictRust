//! Local realignment and insertion adjustment — ports VariationRealigner.java.

use crate::reference::Reference;
use crate::types::{Pos0, Pos1};
use crate::variations::{
    get_variation, get_variation_maybe, get_variation_maybe_mut, Sclip, Variation, VariationMap,
};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, OnceLock};

/// Result of local insertion adjustment.
pub struct BaseInsertion {
    pub base_insert: i32,
    pub insertion_sequence: String,
    pub base_insert2: i32,
}

/// Adjusts the insertion position on a byte slice without allocating, returning
/// `(adjusted_bi, split_idx)` where `split_idx > 0` indicates rotating `ins` at `split_idx`.
/// Ports `VariationRealigner.java:2187-2201`.
pub(crate) fn adj_ins_pos_slice(
    mut bi: i32,
    ins: &[u8],
    reference: &Reference,
    chr: &str,
) -> (i32, usize) {
    let len = ins.len();
    if len == 0 {
        return (bi, 0);
    }
    let mut n = 1usize;
    while let Some(&ins_char) = ins.get(len.saturating_sub(n)) {
        if let Some(ref_char) = u32::try_from(bi)
            .ok()
            .and_then(|p| reference.get_base(chr, p))
        {
            if ref_char == ins_char {
                n += 1;
                if n > len {
                    n = 1;
                }
                bi = bi.saturating_sub(1);
                continue;
            }
        }
        break;
    }
    let split_idx = if n > 1 { len + 1 - n } else { 0 };
    (bi, split_idx)
}

/// Adjusts the insertion position if it falls in a repeat / homopolymer region.
/// Ports `VariationRealigner.java:2187-2201`.
pub fn adj_ins_pos(bi: i32, mut ins: String, reference: &Reference, chr: &str) -> BaseInsertion {
    let (adjusted_bi, split_idx) = adj_ins_pos_slice(bi, ins.as_bytes(), reference, chr);
    if split_idx > 0 && split_idx < ins.len() {
        if ins.is_char_boundary(split_idx) {
            // M-MEM-REUSE: rotating valid UTF-8 at a char boundary stays valid UTF-8 in place.
            let mut bytes = std::mem::take(&mut ins).into_bytes();
            bytes.rotate_left(split_idx);
            if let Ok(rotated) = String::from_utf8(bytes) {
                ins = rotated;
            }
        } else if let (Some(tail), Some(head)) = (
            ins.as_bytes().get(split_idx..),
            ins.as_bytes().get(..split_idx),
        ) {
            let mut new_ins = String::with_capacity(ins.len());
            new_ins.push_str(&String::from_utf8_lossy(tail));
            new_ins.push_str(&String::from_utf8_lossy(head));
            ins = new_ins;
        }
    }
    BaseInsertion {
        base_insert: adjusted_bi,
        insertion_sequence: ins,
        base_insert2: adjusted_bi,
    }
}

#[inline]
fn strip_markers(seq: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if !seq.contains(&b'#') && !seq.contains(&b'^') {
        std::borrow::Cow::Borrowed(seq)
    } else {
        std::borrow::Cow::Owned(
            seq.iter()
                .copied()
                .filter(|&b| b != b'#' && b != b'^')
                .collect(),
        )
    }
}

/// Check if seq1 matches seq2 with no more than MM mismatches and mismatches ratio < 0.15.
///
/// Ports `VariationRealigner.ismatch`.
pub fn is_match(seq1: &[u8], seq2: &[u8], dir: i32, mm_limit: i32) -> bool {
    let clean_seq2 = strip_markers(seq2);

    let mut mm = 0;
    let limit = seq1.len().min(clean_seq2.len());
    for n in 0..limit {
        let idx2 = if dir == -1 {
            clean_seq2.len().saturating_sub(1).saturating_sub(n)
        } else {
            n
        };
        if let (Some(&base1), Some(&base2)) = (seq1.get(n), clean_seq2.get(idx2)) {
            if base1 != base2 {
                mm += 1;
            }
        }
    }

    if seq1.is_empty() {
        return false;
    }
    mm <= mm_limit && (mm as f64) / (seq1.len() as f64) < 0.15
}

/// Check if sequence matches reference starting at position (1-based), in direction dir, with no more than mm_limit mismatches.
///
/// Ports `VariationRealigner.ismatchref`.
pub fn is_match_ref(
    sequence: &[u8],
    reference: &Reference,
    chr: &str,
    position: i32,
    dir: i32,
    mm_limit: i32,
) -> bool {
    let mut mm = 0;
    for n in 0..sequence.len() {
        let ref_ch = match i32::try_from(n)
            .ok()
            .and_then(|ni| dir.checked_mul(ni))
            .and_then(|off| position.checked_add(off))
            .and_then(|p| u32::try_from(p).ok())
            .and_then(|p| reference.get_base(chr, p))
        {
            Some(c) => c,
            None => return false,
        };

        let seq_ch = if dir == -1 {
            sequence[sequence.len() - 1 - n]
        } else {
            sequence[n]
        };

        if seq_ch != ref_ch {
            mm += 1;
        }
    }

    if sequence.is_empty() {
        return false;
    }
    mm <= mm_limit && (mm as f64) / (sequence.len() as f64) < 0.15
}

/// Find breakpoint position.
///
/// Ports `VariationRealigner.findbp`.
pub fn find_bp(
    sequence: &[u8],
    start_position: i32,
    reference: &Reference,
    direction: i32,
    chr: &str,
    indel_size: i32,
    chr_len: u32,
) -> i32 {
    let maxmm = 3;
    let chr_len_i32 = i32::try_from(chr_len).unwrap_or(i32::MAX);
    let mut bp = 0;
    let mut score = 0;

    for n in 0..indel_size {
        let mut mm = 0;
        let mut i = 0;
        let mut m = crate::cigar::modify::ByteSet::new();

        while i < sequence.len() {
            let Some(ref_pos) = i32::try_from(i)
                .ok()
                .and_then(|ii| direction.checked_mul(ii))
                .and_then(|di| {
                    direction
                        .checked_mul(n)
                        .and_then(|dn| start_position.checked_add(dn))
                        .and_then(|p| p.checked_add(di))
                })
            else {
                break;
            };
            if ref_pos < 1 || ref_pos > chr_len_i32 {
                break;
            }

            let ref_ch = match reference.get_base(chr, ref_pos as u32) {
                Some(c) => c,
                None => break,
            };

            if sequence[i] == ref_ch {
                m.insert(sequence[i]);
            } else {
                mm += 1;
            }

            if mm > maxmm - n / 100 {
                break;
            }
            i += 1;
        }

        if m.len() < 3 {
            continue;
        }

        if mm <= maxmm - n / 100
            && i as i32 >= (sequence.len() as i32) - 2
            && i as i32 >= 8 + n / 10
            && i > 0
            && (mm as f64) / (i as f64) < 0.12
        {
            let Ok(lbp) = i32::try_from(
                i64::from(start_position) + i64::from(direction) * i64::from(n)
                    - i64::from(if direction < 0 { direction } else { 0 }),
            ) else {
                continue;
            };
            if mm == 0 && i == sequence.len() {
                return lbp;
            } else if (i as i32) - mm > score {
                bp = lbp;
                score = (i as i32) - mm;
            }
        }
    }

    bp
}

/// Find breakpoint insertion.
///
/// Ports `VariationRealigner.findbi`.
pub fn find_bi(
    seq: &[u8],
    position: i32,
    reference: &Reference,
    dir: i32,
    chr: &str,
    chr_len: u32,
) -> BaseInsertion {
    let maxmm = 3;
    let chr_len_i32 = i32::try_from(chr_len).unwrap_or(i32::MAX);
    let dir_ext = if dir == -1 { 1 } else { 0 };
    let mut score = 0;
    let mut bi = 0;
    let mut ins = String::new();
    let mut bi2 = 0;

    for n in 6..seq.len() {
        if position.saturating_add(6) >= chr_len_i32 {
            break;
        }

        let mut mm = 0;
        let mut i = 0;
        let mut m = crate::cigar::modify::ByteSet::new();

        while i + n < seq.len() {
            let Some(ref_pos) = i32::try_from(i)
                .ok()
                .and_then(|ii| dir.checked_mul(ii))
                .and_then(|di| position.checked_add(di))
                .and_then(|p| p.checked_sub(dir_ext))
            else {
                break;
            };
            if ref_pos < 1 || ref_pos > chr_len_i32 {
                break;
            }

            let ref_ch = match reference.get_base(chr, ref_pos as u32) {
                Some(c) => c,
                None => break,
            };

            if seq[i + n] != ref_ch {
                mm += 1;
            } else {
                m.insert(seq[i + n]);
            }

            if mm > maxmm {
                break;
            }
            i += 1;
        }

        let mnt = m.len() as i32;
        if mnt < 2 {
            continue;
        }

        let ratio = if i > 0 { (mm as f64) / (i as f64) } else { 0.0 };
        if (mnt >= 3 && (i + n) >= seq.len() - 1 && i >= 8 && ratio < 0.15)
            || (mnt >= 2 && mm == 0 && (i + n) == seq.len() && n >= 20 && i >= 8)
        {
            let mut insert = seq[0..n].to_vec();
            let mut extra = Vec::new();
            let mut ept = 0;

            while n + ept + 1 < seq.len() {
                let ref_ch1 = i32::try_from(ept)
                    .ok()
                    .and_then(|ei| ei.checked_mul(dir))
                    .and_then(|off| position.checked_add(off))
                    .and_then(|p| p.checked_sub(dir_ext))
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|p| reference.get_base(chr, p));
                let ref_ch2 = i32::try_from(ept)
                    .ok()
                    .and_then(|ei| ei.checked_add(1))
                    .and_then(|ei1| ei1.checked_mul(dir))
                    .and_then(|off| position.checked_add(off))
                    .and_then(|p| p.checked_sub(dir_ext))
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|p| reference.get_base(chr, p));

                let matches1 = ref_ch1.map(|c| c == seq[n + ept]).unwrap_or(false);
                let matches2 = ref_ch2.map(|c| c == seq[n + ept + 1]).unwrap_or(false);

                if !matches1 || !matches2 {
                    extra.push(seq[n + ept]);
                    ept += 1;
                } else {
                    break;
                }
            }

            if dir == -1 {
                let mut ins_vec = [insert, extra.clone()].concat();
                ins_vec.reverse();
                if !extra.is_empty() {
                    let insert_idx = ins_vec.len() - extra.len();
                    ins_vec.insert(insert_idx, b'&');
                }
                let ins_str = String::from_utf8_lossy(&ins_vec).into_owned();

                if mm == 0 && (i + n) == seq.len() {
                    let mut bi_val = position
                        .saturating_sub(1)
                        .saturating_sub(extra.len() as i32);
                    let mut ins_val = ins_str;
                    let mut bi2_val = position.saturating_sub(1);
                    if extra.is_empty() {
                        let tpl = adj_ins_pos(bi_val, ins_val, reference, chr);
                        bi_val = tpl.base_insert;
                        ins_val = tpl.insertion_sequence;
                        bi2_val = tpl.base_insert2;
                    }
                    return BaseInsertion {
                        base_insert: bi_val,
                        insertion_sequence: ins_val,
                        base_insert2: bi2_val,
                    };
                } else if (i as i32) - mm > score {
                    bi = position
                        .saturating_sub(1)
                        .saturating_sub(extra.len() as i32);
                    ins = ins_str;
                    bi2 = position.saturating_sub(1);
                    score = (i as i32) - mm;
                }
            } else {
                let mut s = -1;
                if !extra.is_empty() {
                    let mut tmp = insert;
                    tmp.push(b'&');
                    tmp.extend(&extra);
                    let ins_str = String::from_utf8_lossy(&tmp).into_owned();

                    if mm == 0 && (i + n) == seq.len() {
                        return BaseInsertion {
                            base_insert: position.saturating_add(s),
                            insertion_sequence: ins_str,
                            base_insert2: position
                                .saturating_add(s)
                                .saturating_add(extra.len() as i32),
                        };
                    } else if (i as i32) - mm > score {
                        bi = position.saturating_add(s);
                        ins = ins_str;
                        bi2 = position
                            .saturating_add(s)
                            .saturating_add(extra.len() as i32);
                        score = (i as i32) - mm;
                    }
                } else {
                    while s >= -(n as i32) {
                        let ins_ch = crate::util::char_at(&insert, s as isize);
                        let ref_ch = position
                            .checked_add(s)
                            .and_then(|p| u32::try_from(p).ok())
                            .and_then(|p| reference.get_base(chr, p));
                        if ins_ch.is_some() && ins_ch == ref_ch {
                            s -= 1;
                        } else {
                            break;
                        }
                    }

                    if s < -1 {
                        let shift = -(s + 1);
                        if shift > 0 {
                            let len_ins = insert.len();
                            let start_idx = len_ins - (shift as usize);
                            let suffix = insert[start_idx..].to_vec();
                            insert.truncate(start_idx);
                            insert = [suffix, insert].concat();
                        }
                    }

                    let ins_str = String::from_utf8_lossy(&insert).into_owned();

                    if mm == 0 && (i + n) == seq.len() {
                        let mut bi_val = position.saturating_add(s);
                        let mut ins_val = ins_str;
                        let mut bi2_val = position
                            .saturating_add(s)
                            .saturating_add(extra.len() as i32);
                        if extra.is_empty() {
                            let tpl = adj_ins_pos(bi_val, ins_val, reference, chr);
                            bi_val = tpl.base_insert;
                            ins_val = tpl.insertion_sequence;
                            bi2_val = tpl.base_insert2;
                        }
                        return BaseInsertion {
                            base_insert: bi_val,
                            insertion_sequence: ins_val,
                            base_insert2: bi2_val,
                        };
                    } else if (i as i32) - mm > score {
                        bi = position.saturating_add(s);
                        ins = ins_str;
                        bi2 = position
                            .saturating_add(s)
                            .saturating_add(extra.len() as i32);
                        score = (i as i32) - mm;
                    }
                }
            }
        }
    }

    if bi2 == bi && !ins.is_empty() && bi != 0 {
        let tpl = adj_ins_pos(bi, ins, reference, chr);
        bi = tpl.base_insert;
        ins = tpl.insertion_sequence;
    }

    BaseInsertion {
        base_insert: bi,
        insertion_sequence: ins,
        base_insert2: bi2,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mismatch {
    pub mismatch_sequence: String,
    pub mismatch_position: i32,
    pub end: i32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MismatchResult {
    pub mismatches: Vec<Mismatch>,
    pub scp: Vec<i32>,
    pub nm: i32,
    pub misp: i32,
    pub misnt: String,
}

/// Find mismatches and softclipping positions at the 5' end.
///
/// Ports `VariationRealigner.findMM5`.
pub fn find_mm5(
    reference: &Reference,
    chr: &str,
    position: i32,
    wupseq: &[u8],
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
) -> MismatchResult {
    let clean_seq2 = strip_markers(wupseq);

    let longmm = 3;
    let mut mismatches = Vec::new();
    let mut n = 0;
    let mut mn: i32 = 0;
    let mut mcnt = 0;
    let mut mm_buf = [0u8; 3];
    let mut mm_len = 0usize;
    let mut sc5p = Vec::new();

    loop {
        let seq_idx = (clean_seq2.len() as i32)
            .saturating_sub(1)
            .saturating_sub(n);
        let Some(&c_seq) = usize::try_from(seq_idx)
            .ok()
            .and_then(|idx| clean_seq2.get(idx))
        else {
            break;
        };
        let Some(ref_pos) = position.checked_sub(n) else {
            break;
        };
        let ref_ch = u32::try_from(ref_pos)
            .ok()
            .and_then(|p| reference.get_base(chr, p));

        if let Some(c_ref) = ref_ch {
            if c_seq != c_ref && mcnt < longmm {
                mm_len += 1;
                let start = mm_buf.len().saturating_sub(mm_len);
                if let Some(slot) = mm_buf.get_mut(start) {
                    *slot = c_seq;
                }
                let mismatch_str = mm_buf
                    .get(start..)
                    .map(|s| String::from_utf8_lossy(s).into_owned())
                    .unwrap_or_default();
                mismatches.push(Mismatch {
                    mismatch_sequence: mismatch_str,
                    mismatch_position: ref_pos,
                    end: 5,
                });
                n += 1;
                mcnt += 1;
                continue;
            }
        }
        break;
    }
    if let Some(p_next) = position.checked_add(1) {
        sc5p.push(p_next);
    }

    let misp = 0;
    let misnt: Option<u8> = None;
    if mm_len == 1 {
        loop {
            let seq_idx = (clean_seq2.len() as i32)
                .saturating_sub(1)
                .saturating_sub(n);
            let Some(&c_seq) = usize::try_from(seq_idx)
                .ok()
                .and_then(|idx| clean_seq2.get(idx))
            else {
                break;
            };
            let ref_ch = position
                .checked_sub(n)
                .and_then(|p| u32::try_from(p).ok())
                .and_then(|p| reference.get_base(chr, p));
            if let Some(c_ref) = ref_ch {
                if c_seq == c_ref {
                    n += 1;
                    if n != 0 {
                        mn += 1;
                    }
                    continue;
                }
            }
            break;
        }

        if mn > 1 {
            // Java's n2-scanning loop guard is `-1 - n - 1 - n2 >= 0`, which is never
            // satisfied for n,n2 >= 0 (VarDictJava v1.8.3 dropped `seq.length()` from the
            // index bound). The loop therefore never executes, n2 stays 0, and the
            // `n2 > 2` branch is never taken — misp/misnt remain unset. Replicate that
            // exactly: always take the else branch (mark the position-n soft clip used).
            if let Some(target_pos) = position.checked_sub(n) {
                sc5p.push(target_pos);
                if let Some(sclip) = soft_clips_5_end.get_mut(&target_pos) {
                    sclip.used = true;
                }
            }
        }
    }

    MismatchResult {
        mismatches,
        scp: sc5p,
        nm: mn,
        misp,
        misnt: misnt.map(|b| (b as char).to_string()).unwrap_or_default(),
    }
}

/// Find mismatches and softclipping positions at the 3' end.
///
/// Ports `VariationRealigner.findMM3`.
pub fn find_mm3(
    reference: &Reference,
    chr: &str,
    p: i32,
    sanpseq: &[u8],
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
) -> MismatchResult {
    let clean_seq2 = strip_markers(sanpseq);

    let longmm = 3;
    let mut mismatches = Vec::new();
    let mut n = 0;
    let mut mn: i32 = 0;
    let mut mcnt = 0;
    let mut sc3p = Vec::new();
    let mut mm_buf = [0u8; 4];
    let mut mm_len = 0usize;

    while let Some(&c_seq) = clean_seq2.get(n) {
        let ref_ch = i32::try_from(n)
            .ok()
            .and_then(|ni| p.checked_add(ni))
            .and_then(|pos| u32::try_from(pos).ok())
            .and_then(|pos| reference.get_base(chr, pos));
        if let Some(c_ref) = ref_ch {
            if c_ref == c_seq {
                n += 1;
                continue;
            }
        }
        break;
    }
    let tbp_opt = i32::try_from(n).ok().and_then(|ni| p.checked_add(ni));
    if let Some(tbp) = tbp_opt {
        sc3p.push(tbp);
    }

    while mcnt <= longmm {
        let (Some(&c_seq), Some(tbp)) = (clean_seq2.get(n), tbp_opt) else {
            break;
        };
        let ref_ch = i32::try_from(n)
            .ok()
            .and_then(|ni| p.checked_add(ni))
            .and_then(|pos| u32::try_from(pos).ok())
            .and_then(|pos| reference.get_base(chr, pos));
        if let Some(c_ref) = ref_ch {
            if c_ref != c_seq {
                if let Some(slot) = mm_buf.get_mut(mm_len) {
                    *slot = c_seq;
                }
                mm_len += 1;
                let mismatch_str = mm_buf
                    .get(..mm_len)
                    .map(|s| String::from_utf8_lossy(s).into_owned())
                    .unwrap_or_default();
                mismatches.push(Mismatch {
                    mismatch_sequence: mismatch_str,
                    mismatch_position: tbp,
                    end: 3,
                });
                n += 1;
                mcnt += 1;
                continue;
            }
        }
        break;
    }

    let mut misp = 0;
    let mut misnt = None;
    if mm_len == 1 {
        while let Some(&c_seq) = clean_seq2.get(n) {
            let ref_ch = i32::try_from(n)
                .ok()
                .and_then(|ni| p.checked_add(ni))
                .and_then(|pos| u32::try_from(pos).ok())
                .and_then(|pos| reference.get_base(chr, pos));
            if let Some(c_ref) = ref_ch {
                if c_seq == c_ref {
                    n += 1;
                    if n != 0 {
                        mn += 1;
                    }
                    continue;
                }
            }
            break;
        }

        if mn > 1 {
            let mut n2 = 0i32;
            while let Some(&c_seq) = usize::try_from(n2)
                .ok()
                .and_then(|n2u| n.checked_add(n2u))
                .and_then(|idx| idx.checked_add(1))
                .and_then(|idx| clean_seq2.get(idx))
            {
                let ref_ch = i32::try_from(n)
                    .ok()
                    .and_then(|ni| p.checked_add(ni))
                    .and_then(|pos| pos.checked_add(1))
                    .and_then(|pos| pos.checked_add(n2))
                    .and_then(|pos| u32::try_from(pos).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                if let Some(c_ref) = ref_ch {
                    if c_seq == c_ref {
                        n2 += 1;
                        continue;
                    }
                }
                break;
            }

            if n2 > 2
                && usize::try_from(n2)
                    .ok()
                    .and_then(|n2u| n.checked_add(n2u))
                    .and_then(|idx| idx.checked_add(1))
                    .is_some_and(|idx| idx < clean_seq2.len())
            {
                if let Some(base_pos) = i32::try_from(n).ok().and_then(|ni| p.checked_add(ni)) {
                    misp = base_pos;
                    misnt = clean_seq2.get(n).copied();
                    if let Some(target_pos) = base_pos.checked_add(n2) {
                        sc3p.push(target_pos);
                        if let Some(sclip) = soft_clips_3_end.get_mut(&target_pos) {
                            sclip.used = true;
                        }
                    }
                }
                mn = mn.saturating_add(n2);
            } else if let Some(target_pos) = i32::try_from(n).ok().and_then(|ni| p.checked_add(ni))
            {
                sc3p.push(target_pos);
                if let Some(sclip) = soft_clips_3_end.get_mut(&target_pos) {
                    sclip.used = true;
                }
            }
        }
    }

    MismatchResult {
        mismatches,
        scp: sc3p,
        nm: mn,
        misp,
        misnt: misnt.map(|b| (b as char).to_string()).unwrap_or_default(),
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SortPositionDescription {
    pub position: i32,
    pub description_string: String,
    pub count: i32,
}

/// Ports `VariationRealigner.fillAndSortTmp`.
pub fn fill_and_sort_tmp(
    changes: &BTreeMap<i32, rustc_hash::FxHashMap<String, i32>>,
) -> Vec<SortPositionDescription> {
    let mut tmp = Vec::new();
    for (&position, map) in changes {
        for (description_string, &count) in map {
            tmp.push(SortPositionDescription {
                position,
                description_string: description_string.clone(),
                count,
            });
        }
    }
    tmp.sort_by(|o1, o2| {
        let f = o2.count.cmp(&o1.count);
        if f != std::cmp::Ordering::Equal {
            return f;
        }
        let f = o1.position.cmp(&o2.position);
        if f != std::cmp::Ordering::Equal {
            return f;
        }
        o2.description_string.cmp(&o1.description_string)
    });
    tmp
}

/// Fill adaptor maps from sequence string list.
///
/// Results are cached keyed by the exact input list: `config.adaptor` is fixed
/// for a run, so per-candidate `realigndel` calls from `realignlgdel` reuse
/// one build. Only the input-keyed maps are cached — `dels5` is always the
/// caller's own and `realigndel` is never skipped.
// M-FAST-HASHER: FxHash for trusted internal keys (`config.adaptor` is
// fixed per run, never attacker-controlled) instead of SipHash.
// Values are `Arc` so a cache hit can take a cheap handle under the lock and
// deep-clone only after the guard drops — otherwise every hit serializes two
// BTreeMap copies across tile workers.
static ADAPTOR_MAP_CACHE: OnceLock<Mutex<rustc_hash::FxHashMap<Vec<String>, Arc<AdaptorMapPair>>>> =
    OnceLock::new();

type AdaptorMapPair = (BTreeMap<Vec<u8>, i32>, BTreeMap<Vec<u8>, i32>);

fn adaptor_map_cache() -> &'static Mutex<rustc_hash::FxHashMap<Vec<String>, Arc<AdaptorMapPair>>> {
    ADAPTOR_MAP_CACHE.get_or_init(|| Mutex::new(rustc_hash::FxHashMap::default()))
}

pub fn fill_adaptor_maps(adaptor_seqs: &[String]) -> AdaptorMapPair {
    // Fast path: the default config carries no adaptors (`Config::default`
    // leaves `adaptor` empty), so every default-config candidate would take
    // the shared-cache lock and hash an empty key for a build that is a
    // no-op. Return the two empty maps directly instead.
    if adaptor_seqs.is_empty() {
        return (BTreeMap::new(), BTreeMap::new());
    }
    // Borrowed lookup: hits (the common case — `config.adaptor` is fixed per
    // run) skip the `to_vec` + per-`String` clones entirely. `Vec<String>:
    // Borrow<[String]>` lets the map be probed with the caller's slice.
    // The guard is scoped so the deep clone below runs unlocked.
    let hit = match adaptor_map_cache().lock() {
        Ok(cache) => cache.get(adaptor_seqs).map(Arc::clone),
        Err(_) => None,
    };
    if let Some(pair) = hit {
        return (*pair).clone();
    }
    let mut adaptor_forward = BTreeMap::new();
    let mut adaptor_reverse = BTreeMap::new();
    let adseed = 6;
    for sequence in adaptor_seqs {
        let seq_bytes = sequence.as_bytes();
        for i in 0..=6 {
            // Ports VarDictLauncher.java:108:
            //   for (int i = 0; i <= 6 && i + Configuration.ADSEED < sequence.length(); i++)
            if i + adseed >= seq_bytes.len() {
                break;
            }
            let forward_seed = seq_bytes[i..i + adseed].to_vec();
            // Ports VarDictLauncher.java:110 -> Utils.java:234-244 ->
            // htsjdk.samtools.util.SequenceUtil.complement(byte), which complements
            // only ACGTacgt and returns any other byte unchanged.
            let comp_reverse_seed: Vec<u8> = forward_seed
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    b'a' => b't',
                    b'c' => b'g',
                    b'g' => b'c',
                    b't' => b'a',
                    other => other,
                })
                .collect();
            adaptor_forward.insert(forward_seed, (i + 1) as i32);
            adaptor_reverse.insert(comp_reverse_seed, (i + 1) as i32);
        }
    }
    let built = (adaptor_forward, adaptor_reverse);
    // shortcut: fixed small cap; per-run inputs are ~1 distinct list, so this
    // only guards pathological caller variation. Upgrade to LRU if ever >cap.
    // The owned key is built only on a miss (hits return above, keyless).
    if let Ok(mut cache) = adaptor_map_cache().lock() {
        if cache.len() < 16 {
            cache.insert(adaptor_seqs.to_vec(), Arc::new(built.clone()));
        }
    }
    built
}

#[derive(Clone, Debug, Default)]
pub(crate) struct AdaptorSet(Option<Arc<rustc_hash::FxHashSet<[u8; 6]>>>);

impl AdaptorSet {
    #[inline]
    pub(crate) fn contains_seed(&self, seed: &[u8; 6]) -> bool {
        self.0.as_ref().is_some_and(|set| set.contains(seed))
    }
}

type AdaptorSetPair = (AdaptorSet, AdaptorSet);

static ADAPTOR_SET_CACHE: OnceLock<Mutex<rustc_hash::FxHashMap<Vec<String>, AdaptorSetPair>>> =
    OnceLock::new();

fn adaptor_set_cache() -> &'static Mutex<rustc_hash::FxHashMap<Vec<String>, AdaptorSetPair>> {
    ADAPTOR_SET_CACHE.get_or_init(|| Mutex::new(rustc_hash::FxHashMap::default()))
}

pub(crate) fn fill_adaptor_sets(adaptor_seqs: &[String]) -> AdaptorSetPair {
    if adaptor_seqs.is_empty() {
        return (AdaptorSet::default(), AdaptorSet::default());
    }
    let hit = match adaptor_set_cache().lock() {
        Ok(cache) => cache.get(adaptor_seqs).cloned(),
        Err(_) => None,
    };
    if let Some(pair) = hit {
        return pair;
    }
    let mut adaptor_forward = rustc_hash::FxHashSet::default();
    let mut adaptor_reverse = rustc_hash::FxHashSet::default();
    let adseed = 6;
    for sequence in adaptor_seqs {
        let seq_bytes = sequence.as_bytes();
        for i in 0..=6 {
            if i + adseed >= seq_bytes.len() {
                break;
            }
            let Some(slice) = seq_bytes.get(i..i + adseed) else {
                break;
            };
            let mut forward_seed = [0u8; 6];
            forward_seed.copy_from_slice(slice);
            let mut comp_reverse_seed = [0u8; 6];
            for (dst, &b) in comp_reverse_seed.iter_mut().zip(forward_seed.iter().rev()) {
                *dst = match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    b'a' => b't',
                    b'c' => b'g',
                    b'g' => b'c',
                    b't' => b'a',
                    other => other,
                };
            }
            adaptor_forward.insert(forward_seed);
            adaptor_reverse.insert(comp_reverse_seed);
        }
    }
    let built = (
        if adaptor_forward.is_empty() {
            AdaptorSet(None)
        } else {
            AdaptorSet(Some(Arc::new(adaptor_forward)))
        },
        if adaptor_reverse.is_empty() {
            AdaptorSet(None)
        } else {
            AdaptorSet(Some(Arc::new(adaptor_reverse)))
        },
    );
    if let Ok(mut cache) = adaptor_set_cache().lock() {
        if cache.len() < 16 {
            cache.insert(adaptor_seqs.to_vec(), built.clone());
        }
    }
    built
}

#[cfg(test)]
mod task3_tests {
    use super::*;
    use rustc_hash::FxHashMap;

    /// s8 INV-SV guard: `fill_adaptor_maps` may only cache keyed by its exact
    /// input. A global (non-keyed) cache would return `ma` for input `b`.
    #[test]
    fn fill_adaptor_maps_is_keyed_by_input_contents() {
        let empty = fill_adaptor_maps(&[]);
        assert!(empty.0.is_empty() && empty.1.is_empty());

        let a = vec!["ACGTACGTACGT".to_string()];
        let b = vec!["TTTTTTTTTTTT".to_string()];
        let ma1 = fill_adaptor_maps(&a);
        let ma2 = fill_adaptor_maps(&a);
        assert_eq!(
            ma1, ma2,
            "same adaptor input must yield identical maps (cache hit must be exact)"
        );
        assert!(ma1.0.keys().all(|k| k.len() == 6));
        let mb = fill_adaptor_maps(&b);
        assert_ne!(
            ma1.0, mb.0,
            "different adaptor input must not reuse a cached forward map"
        );
        assert_ne!(
            ma1.1, mb.1,
            "different adaptor input must not reuse a cached reverse map"
        );
        // Empty input after fills must still return empty maps (fast path
        // skips the cache; the empty key must never collide with a fill).
        let empty_after = fill_adaptor_maps(&[]);
        assert!(empty_after.0.is_empty() && empty_after.1.is_empty());

        let empty_sets = fill_adaptor_sets(&[]);
        assert!(!empty_sets.0.contains_seed(b"ACGTAC"));
        let sa1 = fill_adaptor_sets(&a);
        let sa2 = fill_adaptor_sets(&a);
        assert!(sa1.0.contains_seed(b"ACGTAC") && sa2.0.contains_seed(b"ACGTAC"));
        assert!(sa1.1.contains_seed(b"GTACGT") && sa2.1.contains_seed(b"GTACGT"));
        let sb = fill_adaptor_sets(&b);
        assert!(sb.0.contains_seed(b"TTTTTT") && !sb.0.contains_seed(b"ACGTAC"));
        assert!(sb.1.contains_seed(b"AAAAAA") && !sb.1.contains_seed(b"GTACGT"));
    }

    /// Ports `VarDictLauncher.java:106-117` (`i <= 6 && i + ADSEED < sequence.length()`)
    /// and `Utils.java:228-244` (`htsjdk.samtools.util.SequenceUtil.complement(byte)`).
    #[test]
    fn fill_adaptor_maps_java_seed_count() {
        let base_alphabet = b"ACGTGHKMRSTVWY";
        for len in 0..=14usize {
            let seq: String = (0..len)
                .map(|idx| base_alphabet[(idx * 3 + len) % base_alphabet.len()] as char)
                .collect();
            let (fwd, rev) = fill_adaptor_maps(std::slice::from_ref(&seq));
            let expected_count = (0..=6usize).filter(|&i| i + 6 < len).count();
            assert_eq!(
                fwd.len(),
                expected_count,
                "forward seed count mismatch for adaptor length {len} ({seq})"
            );
            assert_eq!(
                rev.len(),
                expected_count,
                "reverse seed count mismatch for adaptor length {len} ({seq})"
            );
        }

        for (len, expected) in [(6usize, 0usize), (9, 3), (10, 4), (12, 6), (13, 7), (33, 7)] {
            let seq: String = (0..len)
                .map(|i| {
                    let b = match i % 4 {
                        0 => b'A',
                        1 => b'C',
                        2 => b'G',
                        _ => b'T',
                    };
                    // Make every 6-mer distinct
                    (if i == 6 {
                        b'C'
                    } else if i == 7 {
                        b'A'
                    } else {
                        b
                    }) as char
                })
                .collect();
            let (fwd, _) = fill_adaptor_maps(&[seq]);
            let max_i = fwd.values().copied().max().unwrap_or(0) as usize;
            assert_eq!(max_i, expected, "max seed index 1-based for len={len}");
        }

        // 12-bp adaptor `AAAAAATTGCTC`: Java loop stops at i=5 (`5 + 6 < 12`),
        // so trailing 6-mer at i=6 (`TTGCTC`) is never inserted.
        let (fwd12, _) = fill_adaptor_maps(&["AAAAAATTGCTC".to_string()]);
        assert_eq!(fwd12.len(), 6);
        assert!(
            !fwd12.contains_key(b"TTGCTC".as_slice()),
            "12-bp adaptor must not include trailing 6-mer at i=6"
        );

        // Non-ACGT IUPAC byte ('R') in forward seed `ACGTRG` (reversed `GRTGCA`):
        // htsjdk SequenceUtil.complement(byte) complements only ACGTacgt and leaves
        // 'R' unchanged -> `CRACGT` (not `CNACGT`).
        let (fwd_iupac, rev_iupac) = fill_adaptor_maps(&["ACGTRGCTCAAAA".to_string()]);
        assert_eq!(fwd_iupac.get(b"ACGTRG".as_slice()), Some(&1));
        assert_eq!(rev_iupac.get(b"CRACGT".as_slice()), Some(&1));
        assert!(!rev_iupac.contains_key(b"CNACGT".as_slice()));

        // Trailing empty element (as produced when CLI splits on `,` cleanly or skips empty)
        let (fwd_clean, rev_clean) = fill_adaptor_maps(&["ACGTTGCTC".to_string()]);
        let (fwd_empty_tail, rev_empty_tail) =
            fill_adaptor_maps(&["ACGTTGCTC".to_string(), String::new()]);
        assert_eq!(fwd_clean, fwd_empty_tail);
        assert_eq!(rev_clean, rev_empty_tail);
    }

    fn task3_reference() -> Reference {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 2000,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 2000]);
        reference
    }

    fn task3_sclip(consensus: u8, cnt: i32) -> Sclip {
        let mut sc = Sclip::new();
        sc.variation.cnt = cnt;
        for idx in 0..8 {
            let mut nt = BTreeMap::new();
            nt.insert(consensus, 10);
            sc.nt.insert(idx, nt);
        }
        sc
    }

    /// s8 INV-SV branch pin (`-N^SEQ<INVn>SEQ` shape): `realigndel` must take
    /// the `inv5`/`inv3` override path and must consume each call's own `dels5`
    /// input. Regresses if `dels5` is cached across calls or the call is
    /// bypassed (soft clips stay unused, counts never move).
    #[test]
    fn realigndel_inv_sv_branch_consumes_per_call_dels5() {
        let reference = task3_reference();
        let config = crate::config::Config::default();
        let mut non_ins = crate::javamap::JavaHashMap::new();
        let mut ref_cov = rustc_hash::FxHashMap::default();
        let mut sc5 = BTreeMap::new();
        let mut sc3 = BTreeMap::new();

        // Call 1: INV deletion at p1; 5' soft clip consensus matches inv3.
        let p1 = 1000;
        let vn1 = "-50^AAAAAAAA<INV123>TTTTTTTT".to_string();
        let mut dels5_a: BTreeMap<i32, FxHashMap<String, i32>> = BTreeMap::new();
        let mut inner_a = FxHashMap::default();
        inner_a.insert(vn1.clone(), 7);
        dels5_a.insert(p1, inner_a);
        sc5.insert(1050, task3_sclip(b'T', 5));

        realigndel(
            None,
            &dels5_a,
            &reference,
            "chr1",
            2000,
            &config,
            &mut non_ins,
            &mut ref_cov,
            &mut sc5,
            &mut sc3,
            150,
        );

        assert!(
            sc5.get(&1050).map(|s| s.used).unwrap_or(false),
            "realigndel must absorb the 5' soft clip via the inv3 override (call skipped?)"
        );
        let cnt1 = non_ins
            .get(&p1)
            .and_then(|e| e.map.get(vn1.as_str()))
            .map(|v| v.cnt)
            .unwrap_or(-999);
        assert_eq!(cnt1, 5, "vn1 must absorb the soft-clip count");
        assert_eq!(ref_cov.get(&p1).copied().unwrap_or(-999), 5);

        // Call 2: a *different* dels5 at p2. If dels5 were cached, vn2 would
        // never apply and the fresh soft clip would stay unused.
        let p2 = 1200;
        let vn2 = "-30^CCCCCCCC<INV7>GGGGGGGG".to_string();
        let mut dels5_b: BTreeMap<i32, FxHashMap<String, i32>> = BTreeMap::new();
        let mut inner_b = FxHashMap::default();
        inner_b.insert(vn2.clone(), 7);
        dels5_b.insert(p2, inner_b);
        sc5.insert(1230, task3_sclip(b'G', 3));

        realigndel(
            None,
            &dels5_b,
            &reference,
            "chr1",
            2000,
            &config,
            &mut non_ins,
            &mut ref_cov,
            &mut sc5,
            &mut sc3,
            150,
        );

        assert!(
            sc5.get(&1230).map(|s| s.used).unwrap_or(false),
            "second dels5 input must be honored (dels5 must not be cached)"
        );
        let cnt2 = non_ins
            .get(&p2)
            .and_then(|e| e.map.get(vn2.as_str()))
            .map(|v| v.cnt)
            .unwrap_or(-999);
        assert_eq!(cnt2, 3, "vn2 must absorb the second soft-clip count");
        assert_eq!(ref_cov.get(&p2).copied().unwrap_or(-999), 3);
        // First call's effects intact and not double-applied.
        let cnt1_again = non_ins
            .get(&p1)
            .and_then(|e| e.map.get(vn1.as_str()))
            .map(|v| v.cnt)
            .unwrap_or(-999);
        assert_eq!(cnt1_again, 5);
    }
}

#[inline]
fn cigar_contains_dlen_d(
    cigar: &[(u32, noodles::sam::alignment::record::cigar::op::Kind)],
    dlen: i32,
) -> bool {
    use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
    if dlen < 0 {
        return false;
    }
    let d = dlen as u64;
    let m = 10u64.pow(d.checked_ilog10().unwrap_or(0) + 1);
    cigar
        .iter()
        .any(|&(len, k)| k == CigarKind::Deletion && (len as u64) >= d && (len as u64) % m == d)
}

#[inline]
fn accum_no_passing_read_record(
    record: &crate::bam::SamRecord,
    dlen: i32,
    start: i32,
    end: i32,
    cnt: &mut i32,
    midcnt: &mut i32,
) {
    if cigar_contains_dlen_d(&record.cigar, dlen) {
        return;
    }
    let read_start = match record.alignment_start {
        Some(Pos1(p)) => i32::try_from(p).unwrap_or(i32::MAX),
        None => return,
    };
    let read_len = crate::bam::get_aligned_length_md(&record.cigar);
    let read_end = read_start.saturating_add(read_len);

    if read_end > end.saturating_add(2) && read_start < start.saturating_sub(2) {
        *cnt = cnt.saturating_add(1);
    }
    if read_start < start.saturating_sub(2) && read_end > start && read_end < end {
        *midcnt = midcnt.saturating_add(1);
    }
}

/// Check if there are no passing reads spanning the start-end coordinates.
///
/// Ports `VariationRealigner.noPassingReads`.
pub fn no_passing_reads(
    chr: &str,
    start: i32,
    end: i32,
    bams: &[String],
    config: &crate::config::Config,
) -> bool {
    let mut cnt = 0;
    let mut midcnt = 0;
    let dlen = end.saturating_sub(start);

    'bam: for bam in bams {
        let bam_path = std::path::Path::new(bam);
        let start_pos0 = Pos0(start.saturating_sub(1).max(0) as u32);
        let end_pos0 = Pos0(end.max(0) as u32);
        // Streaming counterpart of the old query_region-Vec shape: the whole
        // BAM is still skipped on any query/read error, exactly as the old
        // collecting query did. The old path failed atomically (no records
        // from a failed BAM were ever counted), so snapshot both counters
        // and restore them when the feed fails mid-query.
        let saved_cnt = cnt;
        let saved_midcnt = midcnt;
        let mut feed =
            match crate::bam::query_region_iter(bam_path, chr, start_pos0, end_pos0, config) {
                Ok(feed) => feed,
                Err(e) => {
                    log::warn!("no_passing_reads: failed to query bam '{}': {}", bam, e);
                    continue;
                }
            };

        let mut record = crate::bam::SamRecord::new();
        loop {
            match feed.next_into(&mut record) {
                Ok(false) => break,
                Ok(true) => {}
                Err(e) => {
                    log::warn!("no_passing_reads: failed to query bam '{}': {}", bam, e);
                    cnt = saved_cnt;
                    midcnt = saved_midcnt;
                    continue 'bam;
                }
            }
            accum_no_passing_read_record(&record, dlen, start, end, &mut cnt, &mut midcnt);
        }
    }

    if config.verbose {
        crate::verbose_eprintln!(
            "    Passing Read CNT: {} {} {} {} {}",
            cnt,
            chr,
            start,
            end,
            midcnt
        );
    }

    cnt <= 0 && midcnt.saturating_add(1) > 0
}

#[cfg(test)]
mod task13_tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

    #[test]
    fn no_passing_reads_saturates_huge_cigar_without_i32_overflow() {
        let mut record = crate::bam::SamRecord::new();
        record.alignment_start = Some(Pos1(10));
        record.cigar = vec![
            (2_000_000_000, CigarKind::Match),
            (500_000_000, CigarKind::Deletion),
        ];
        let mut cnt = 0;
        let mut midcnt = 0;
        accum_no_passing_read_record(&record, 50, 100, 150, &mut cnt, &mut midcnt);
        assert_eq!(cnt, 1);
        assert_eq!(midcnt, 0);

        // Extreme negative coordinates (`start = i32::MIN`, `end = -5`) must not
        // overflow `start - 1`, wrap `end as u32`, or overflow `midcnt + 1`.
        let mut neg_cnt = 0;
        let mut neg_midcnt = i32::MAX;
        accum_no_passing_read_record(
            &record,
            (-5i32).saturating_sub(i32::MIN),
            i32::MIN,
            -5,
            &mut neg_cnt,
            &mut neg_midcnt,
        );
        assert_eq!(neg_cnt, 0);
        assert_eq!(neg_midcnt, i32::MAX);
        assert!(no_passing_reads(
            "chr1",
            i32::MIN,
            -5,
            &["__nonexistent_bam_for_i32_min_guard__.bam".to_string()],
            &crate::config::Config::default(),
        ));
    }

    fn cigar_contains_dlen_d_legacy(cigar: &[(u32, CigarKind)], dlen: i32) -> bool {
        let dlenqr = format!("{dlen}D");
        let cigar_str = crate::bam::get_cigar_string(cigar);
        cigar_str.contains(&dlenqr)
    }

    #[test]
    fn cigar_contains_dlen_d_matches_string_contains() {
        let kinds = [
            CigarKind::Match,
            CigarKind::Insertion,
            CigarKind::Deletion,
            CigarKind::Skip,
            CigarKind::SoftClip,
            CigarKind::HardClip,
            CigarKind::Pad,
            CigarKind::SequenceMatch,
            CigarKind::SequenceMismatch,
        ];

        // 1. Exhaustive single-op: len, dlen in 0..2000 across Deletion and non-Deletion ops
        for len in 0..2000u32 {
            let del_cigar = [(len, CigarKind::Deletion)];
            let match_cigar = [(len, CigarKind::Match)];
            for dlen in 0..2000i32 {
                assert_eq!(
                    cigar_contains_dlen_d(&del_cigar, dlen),
                    cigar_contains_dlen_d_legacy(&del_cigar, dlen),
                    "single-op del mismatch for len={len}, dlen={dlen}"
                );
                assert!(!cigar_contains_dlen_d(&match_cigar, dlen));
            }
        }

        // 2. Explicit substring parity checks (`15D` vs `5D`, `11D` vs `1D`, `21D` vs `1D`, start/mid/multi-D)
        assert!(cigar_contains_dlen_d(&[(15, CigarKind::Deletion)], 5));
        assert!(cigar_contains_dlen_d(&[(11, CigarKind::Deletion)], 1));
        assert!(cigar_contains_dlen_d(&[(21, CigarKind::Deletion)], 1));
        assert!(!cigar_contains_dlen_d(&[(5, CigarKind::Deletion)], 15));
        assert!(cigar_contains_dlen_d(
            &[(10, CigarKind::Deletion), (50, CigarKind::Match)],
            10
        ));
        assert!(cigar_contains_dlen_d(
            &[
                (50, CigarKind::Match),
                (10, CigarKind::Deletion),
                (50, CigarKind::Match)
            ],
            10
        ));
        assert!(cigar_contains_dlen_d(
            &[
                (10, CigarKind::Deletion),
                (20, CigarKind::Match),
                (100, CigarKind::Deletion)
            ],
            100
        ));

        // 3. Extremes: len in {0, 1, 99999, 100000, u32::MAX}, dlen in {0, -1, -3, i32::MIN, i32::MAX}
        let extreme_lens = [0u32, 1, 5, 15, 99_999, 100_000, i32::MAX as u32, u32::MAX];
        let extreme_dlens = [
            i32::MIN,
            -100,
            -3,
            -1,
            0,
            1,
            5,
            15,
            95,
            99_999,
            100_000,
            i32::MAX,
        ];
        assert!(!cigar_contains_dlen_d(&[], 5));
        assert!(!cigar_contains_dlen_d(&[], 0));
        for &len in &extreme_lens {
            for &k in &kinds {
                let c = [(len, k)];
                for &dlen in &extreme_dlens {
                    assert_eq!(
                        cigar_contains_dlen_d(&c, dlen),
                        cigar_contains_dlen_d_legacy(&c, dlen),
                        "extreme mismatch for len={len}, kind={k:?}, dlen={dlen}"
                    );
                }
            }
        }

        // 4. 50,000 seeded multi-op CIGARs
        let mut rng: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut cigar = Vec::with_capacity(6);
        for case_idx in 0..50_000usize {
            cigar.clear();
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let n_ops = ((rng >> 32) as usize % 6) + 1;
            for _ in 0..n_ops {
                rng = rng
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let kind = kinds[(rng >> 16) as usize % kinds.len()];
                let len = match (rng >> 24) & 3 {
                    0 => (rng as u32) % 25,
                    1 => (rng as u32) % 500,
                    2 => (rng as u32) % 100_000,
                    _ => rng as u32,
                };
                cigar.push((len, kind));
            }
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let dlen = match (rng >> 40) & 3 {
                0 => -(((rng as i32).abs() % 50) + 1),
                1 => (rng as u32 % 25) as i32,
                2 => (rng as u32 % 500) as i32,
                _ => (rng as u32 % 100_000) as i32,
            };
            assert_eq!(
                cigar_contains_dlen_d(&cigar, dlen),
                cigar_contains_dlen_d_legacy(&cigar, dlen),
                "random multi-op mismatch at #{case_idx}: cigar={cigar:?}, dlen={dlen}"
            );
        }
    }
}

/// Realign deletions.
///
/// Ports `VariationRealigner.realigndel`.
pub fn realigndel(
    bams: Option<&[String]>,
    position_to_deletion_count: &BTreeMap<i32, rustc_hash::FxHashMap<String, i32>>,
    reference: &Reference,
    chr: &str,
    chr_len: u32,
    config: &crate::config::Config,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut rustc_hash::FxHashMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_sets(&config.adaptor);
    let tmp = fill_and_sort_tmp(position_to_deletion_count);
    let mut extra_buf = String::new();
    let mut wupseq = String::new();
    let mut sanpseq = String::new();

    for tpl in &tmp {
        let p = tpl.position;
        let vn = &tpl.description_string;
        let dcnt = tpl.count;

        if config.verbose {
            let cov = ref_coverage.get(&p).copied().unwrap_or(0);
            crate::verbose_eprintln!("  Realigndel for: {} {} {} cov: {}", p, vn, dcnt, cov);
        }

        let mut dellen = 0;
        if let Some(m) = parse_begin_minus_number(vn) {
            dellen = m.parse::<i32>().unwrap_or(0);
        }
        if let Some(m) = parse_up_number_end(vn) {
            dellen = dellen.saturating_add(m.parse::<i32>().unwrap_or(0));
        }

        let mut extrains = "";
        extra_buf.clear();
        let mut inv5 = "";
        let mut inv3 = "";

        if let Some((i5, i3)) = parse_minus_number_atgnc_sv_atgnc_end(vn) {
            inv5 = i5;
            inv3 = i3;
        } else if let Some((_digits, rest)) = parse_begin_minus_number_any(vn) {
            extra_buf.reserve(rest.len());
            for ch in rest.chars() {
                if !matches!(ch, '^' | '&' | '#') {
                    extra_buf.push(ch);
                }
            }
            extrains = parse_caret_atgnc(vn).unwrap_or("");
        }

        let chr_len_i32 = i32::try_from(chr_len).unwrap_or(i32::MAX);
        let wustart = p.saturating_sub(200).max(1);
        if !inv3.is_empty() {
            wupseq.clear();
            wupseq.push_str(inv3);
        } else {
            crate::util::ref_bases_into(&mut wupseq, reference, chr, wustart, p.saturating_sub(1));
            wupseq.push_str(&extra_buf);
        }

        let sanend = p.saturating_add(200).min(chr_len_i32);
        let start_pos_ref64 =
            i64::from(p) + i64::from(dellen) + (extra_buf.len() as i64) - (extrains.len() as i64);
        sanpseq.clear();
        if !inv5.is_empty() {
            sanpseq.push_str(inv5);
        } else {
            sanpseq.push_str(&extra_buf);
            if start_pos_ref64 <= i64::from(sanend) {
                if let Ok(start_pos_ref) = i32::try_from(start_pos_ref64) {
                    crate::util::push_ref_bases(
                        &mut sanpseq,
                        reference,
                        chr,
                        start_pos_ref,
                        sanend,
                    );
                }
            }
        }

        let r3 = find_mm3(reference, chr, p, sanpseq.as_bytes(), soft_clips_3_end);
        let r5 = if let Ok(r5_pos) = i32::try_from(start_pos_ref64 - 1) {
            find_mm5(reference, chr, r5_pos, wupseq.as_bytes(), soft_clips_5_end)
        } else {
            MismatchResult::default()
        };

        let MismatchResult {
            mismatches: mm3,
            scp: ref sc3p,
            nm: nm3,
            misp: misp3,
            misnt: ref misnt3,
        } = r3;
        let MismatchResult {
            mismatches: mm5,
            scp: ref sc5p,
            nm: nm5,
            misp: misp5,
            misnt: ref misnt5,
        } = r5;
        let mm3_len_is_one = mm3.len() == 1;
        let mm5_len_is_one = mm5.len() == 1;

        if config.verbose {
            crate::verbose_eprintln!(
                "  Mismatches: misp3: {}-{} misp5: {}-{} sclip3: {:?} sclip5: {:?}",
                misp3,
                misnt3,
                misp5,
                misnt5,
                sc3p,
                sc5p
            );
        }

        let mut mmm = mm3;
        mmm.extend(mm5);

        for mismatch in mmm {
            let mut mm = mismatch.mismatch_sequence;
            let mp = mismatch.mismatch_position;
            let me = mismatch.end;

            if let (Some(head), Some(tail)) = (mm.get(..1), mm.get(1..)) {
                if !tail.is_empty() {
                    mm = format!("{head}&{tail}");
                }
            }

            let mut should_adjust = false;
            if let Some(entry) = non_insertion_variants.get(&mp) {
                if let Some(tv) = entry.map.get(&mm) {
                    if tv.cnt > 0 {
                        let mean_qual = tv.qmean / (tv.cnt as f64);
                        let mean_pos = tv.pmean / (tv.cnt as f64);
                        let threshold_pos = if me == 3 {
                            nm3.saturating_add(4)
                        } else {
                            nm5.saturating_add(4)
                        };
                        if mean_qual >= config.base_qual
                            && mean_pos <= (threshold_pos as f64)
                            && (tv.cnt < dcnt.saturating_add(dellen)
                                && (tv.cnt as f64) / (dcnt as f64) < 8.0)
                        {
                            should_adjust = true;
                        }
                    }
                }
            }

            if !should_adjust {
                continue;
            }

            let tv = match non_insertion_variants
                .get_mut(&mp)
                .and_then(|entry| entry.map.swap_remove(&mm))
            {
                Some(tv) => tv,
                None => continue,
            };

            if config.verbose {
                let cov = ref_coverage.get(&p).copied().unwrap_or(0);
                crate::verbose_eprintln!(
                    "  Realigndel Adj: {} {} {} {} {} {} {} {} cov: {}",
                    mm,
                    mp,
                    me,
                    nm3,
                    nm5,
                    p,
                    tv.cnt,
                    tv.qmean,
                    cov
                );
            }

            if mp > p && me == 5 {
                let mean_pos = if tv.cnt != 0 {
                    tv.pmean / (tv.cnt as f64)
                } else {
                    0.0
                };
                let f = if mean_pos != 0.0 {
                    ((mp - p) as f64) / mean_pos
                } else {
                    1.0
                };
                let mut f = f;
                if f > 1.0 {
                    f = 1.0;
                }
                let cov_entry = ref_coverage.entry(p).or_default();
                *cov_entry = cov_entry.saturating_add((tv.cnt as f64 * f) as i32);

                if let Some(ref_base) = reference.get_base(chr, p as u32) {
                    if let Some(ref_var) =
                        get_variation_maybe_mut(non_insertion_variants, p, Some(ref_base))
                    {
                        crate::variations::util::adj_ref_cnt(&tv, ref_var, dellen);
                    }
                }
            }

            let vref = get_variation(non_insertion_variants, p, vn);
            crate::variations::util::adj_cnt(vref, &tv, None);

            if mp > p && me == 3 {
                if let Some(ref_base) = reference.get_base(chr, p as u32) {
                    if let Some(lref) =
                        get_variation_maybe_mut(non_insertion_variants, p, Some(ref_base))
                    {
                        lref.cnt = lref.cnt.saturating_sub(tv.cnt);
                        lref.hicnt = lref.hicnt.saturating_sub(tv.hicnt);
                        lref.locnt = lref.locnt.saturating_sub(tv.locnt);
                        lref.pmean -= tv.pmean;
                        lref.qmean -= tv.qmean;
                        lref.mqmean -= tv.mqmean;
                        lref.nm -= tv.nm;
                        lref.sub_dir(true, tv.get_dir(true));
                        lref.sub_dir(false, tv.get_dir(false));
                        crate::variations::util::correct_cnt(lref);
                    }
                }
            }

            if let Some(entry) = non_insertion_variants.get(&mp) {
                if entry.map.is_empty() {
                    non_insertion_variants.remove(&mp);
                }
            }

            if config.verbose {
                let cov = ref_coverage.get(&p).copied().unwrap_or(0);
                let final_vref = get_variation(non_insertion_variants, p, vn);
                crate::verbose_eprintln!(
                    "  Realigndel AdjA: {} {} {} {} {} {} {} {} cov: {}",
                    mm,
                    mp,
                    me,
                    nm3,
                    nm5,
                    p,
                    final_vref.cnt,
                    final_vref.qmean,
                    cov
                );
            }
        }

        if misp3 != 0 && mm3_len_is_one {
            let mut remove_misp3 = false;
            if let Some(entry) = non_insertion_variants.get(&misp3) {
                if let Some(var) = entry.map.get(misnt3) {
                    if var.cnt < dcnt {
                        remove_misp3 = true;
                    }
                }
            }
            if remove_misp3 {
                if let Some(entry) = non_insertion_variants.get_mut(&misp3) {
                    entry.map.swap_remove(misnt3);
                    if entry.map.is_empty() {
                        non_insertion_variants.remove(&misp3);
                    }
                }
            }
        }

        if misp5 != 0 && mm5_len_is_one {
            let mut remove_misp5 = false;
            if let Some(entry) = non_insertion_variants.get(&misp5) {
                if let Some(var) = entry.map.get(misnt5) {
                    if var.cnt < dcnt {
                        remove_misp5 = true;
                    }
                }
            }
            if remove_misp5 {
                if let Some(entry) = non_insertion_variants.get_mut(&misp5) {
                    entry.map.swap_remove(misnt5);
                    if entry.map.is_empty() {
                        non_insertion_variants.remove(&misp5);
                    }
                }
            }
        }

        for &sc5pp in sc5p {
            if let Some(tv) = soft_clips_5_end.get_mut(&sc5pp) {
                if !tv.used {
                    let seq = crate::variations::util::find_conseq_sets(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if dcnt <= 2 && tv.variation.cnt.checked_div(dcnt).is_none_or(|r| r > 5) {
                        continue;
                    }
                    if config.verbose {
                        let cov = ref_coverage.get(&p).copied().unwrap_or(0);
                        let wupseq_rev: String = wupseq.chars().rev().collect();
                        crate::verbose_eprintln!(
                            "  Realigndel 5: {} {} seq: '{}' Wuseq: {} cnt: {} {} {} {} cov: {}",
                            p,
                            sc5pp,
                            String::from_utf8_lossy(&seq),
                            wupseq_rev,
                            tv.variation.cnt,
                            dcnt,
                            vn,
                            p,
                            cov
                        );
                    }
                    if !seq.is_empty() && is_match(&seq, wupseq.as_bytes(), -1, 3) {
                        if sc5pp > p {
                            let cov_entry = ref_coverage.entry(p).or_default();
                            *cov_entry = cov_entry.saturating_add(tv.variation.cnt);
                        }
                        let vref = get_variation(non_insertion_variants, p, vn);
                        crate::variations::util::adj_cnt(vref, &tv.variation, None);

                        if let Some(tv_mut) = soft_clips_5_end.get_mut(&sc5pp) {
                            tv_mut.used = true;
                        }

                        if config.verbose {
                            let cov = ref_coverage.get(&p).copied().unwrap_or(0);
                            let wupseq_rev: String = wupseq.chars().rev().collect();
                            let final_tv_cnt = soft_clips_5_end
                                .get(&sc5pp)
                                .map_or(0, |sc| sc.variation.cnt);
                            crate::verbose_eprintln!(
                                "  Realigndel 5: {} {} {} {} {} {} {} {} used cov: {}",
                                p,
                                sc5pp,
                                String::from_utf8_lossy(&seq),
                                wupseq_rev,
                                final_tv_cnt,
                                dcnt,
                                vn,
                                p,
                                cov
                            );
                        }
                    }
                }
            }
        }

        for &sc3pp in sc3p {
            if let Some(tv) = soft_clips_3_end.get_mut(&sc3pp) {
                if !tv.used {
                    let seq = crate::variations::util::find_conseq_sets(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if dcnt <= 2 && tv.variation.cnt.checked_div(dcnt).is_none_or(|r| r > 5) {
                        continue;
                    }

                    let slice_start = sc3pp.saturating_sub(p).max(0) as usize;
                    let sanpseq_slice = sanpseq.get(slice_start..).unwrap_or("");

                    if config.verbose {
                        crate::verbose_eprintln!(
                            "  Realigndel 3: {} {} seq '{}' Sanseq: {} cnt: {} {} {} {} {} {}",
                            p,
                            sc3pp,
                            String::from_utf8_lossy(&seq),
                            sanpseq,
                            tv.variation.cnt,
                            dcnt,
                            vn,
                            p,
                            dellen,
                            sanpseq_slice
                        );
                    }

                    if !seq.is_empty() && is_match(&seq, sanpseq_slice.as_bytes(), 1, 3) {
                        if config.verbose {
                            crate::verbose_eprintln!(
                                "  Realigndel 3: {} {} {} {} {} {} {} {} used",
                                p,
                                sc3pp,
                                String::from_utf8_lossy(&seq),
                                sanpseq,
                                tv.variation.cnt,
                                dcnt,
                                vn,
                                p
                            );
                        }
                        if sc3pp <= p {
                            let cov_entry = ref_coverage.entry(p).or_default();
                            *cov_entry = cov_entry.saturating_add(tv.variation.cnt);
                        }

                        let mut lref_owned = None;
                        if sc3pp > p {
                            if let Some(ref_base) = reference.get_base(chr, p as u32) {
                                if let Some(lref) = get_variation_maybe_mut(
                                    non_insertion_variants,
                                    p,
                                    Some(ref_base),
                                ) {
                                    lref_owned = Some(lref.clone());
                                }
                            }
                        }

                        let vref = get_variation(non_insertion_variants, p, vn);
                        crate::variations::util::adj_cnt(vref, &tv.variation, None);

                        if let Some(mut lref) = lref_owned {
                            lref.cnt = lref.cnt.saturating_sub(tv.variation.cnt);
                            lref.hicnt = lref.hicnt.saturating_sub(tv.variation.hicnt);
                            lref.locnt = lref.locnt.saturating_sub(tv.variation.locnt);
                            lref.pmean -= tv.variation.pmean;
                            lref.qmean -= tv.variation.qmean;
                            lref.mqmean -= tv.variation.mqmean;
                            lref.nm -= tv.variation.nm;
                            lref.sub_dir(true, tv.variation.get_dir(true));
                            lref.sub_dir(false, tv.variation.get_dir(false));
                            crate::variations::util::correct_cnt(&mut lref);
                            if let Some(ref_base) = reference.get_base(chr, p as u32) {
                                if let Some(lref_mut) = get_variation_maybe_mut(
                                    non_insertion_variants,
                                    p,
                                    Some(ref_base),
                                ) {
                                    *lref_mut = lref;
                                }
                            }
                        }

                        if let Some(tv_mut) = soft_clips_3_end.get_mut(&sc3pp) {
                            tv_mut.used = true;
                        }
                    }
                }
            }
        }

        let pe64 = start_pos_ref64;
        let pe = pe64.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;

        let mut h_owned = None;
        if let Some(ref_base) = reference.get_base(chr, p as u32) {
            if let Some(h) = get_variation_maybe(non_insertion_variants, p, Some(ref_base)) {
                h_owned = Some(h.clone());
            }
        }

        if let Some(bams_slice) = bams {
            let pe_diff64 = pe64 - i64::from(p);
            if !bams_slice.is_empty()
                && pe_diff64 >= 5
                && pe_diff64 < i64::from(max_read_length.saturating_sub(10))
            {
                if let Some(ref h) = h_owned {
                    if h.cnt != 0 {
                        let vref = get_variation(non_insertion_variants, p, vn);
                        let threshold_vref_cnt = 2.0
                            * (h.cnt as f64)
                            * (1.0 - (pe_diff64 as f64) / (max_read_length as f64));
                        if (vref.cnt as f64) > threshold_vref_cnt
                            && no_passing_reads(chr, p, pe, bams_slice, config)
                        {
                            crate::variations::util::adj_cnt(vref, h, None);
                            if let Some(ref_base) = reference.get_base(chr, p as u32) {
                                if let Some(h_mut) = get_variation_maybe_mut(
                                    non_insertion_variants,
                                    p,
                                    Some(ref_base),
                                ) {
                                    h_mut.cnt = 0;
                                    h_mut.hicnt = 0;
                                    h_mut.locnt = 0;
                                    h_mut.pmean = 0.0;
                                    h_mut.qmean = 0.0;
                                    h_mut.mqmean = 0.0;
                                    h_mut.nm = 0.0;
                                    h_mut.rev = 0;
                                    h_mut.fwd = 0;
                                    crate::variations::util::correct_cnt(h_mut);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let tmp_size = tmp.len();
    if tmp_size > 1 {
        for i in (1..tmp_size).rev() {
            let tpl = &tmp[i];
            let p = tpl.position;
            let vn = &tpl.description_string;
            if !non_insertion_variants.contains_key(&p) {
                continue;
            }

            if let Some(tn) = parse_minus_number_amp_atgcs_end(vn) {
                let mut vref_cnt = 0;
                let mut vref_exists = false;
                if let Some(entry) = non_insertion_variants.get(&p) {
                    if let Some(vref) = entry.map.get(vn) {
                        vref_cnt = vref.cnt;
                        vref_exists = true;
                    }
                }

                if !vref_exists {
                    continue;
                }

                let mut tref_cnt = 0;
                let mut tref_exists = false;
                if let Some(entry) = non_insertion_variants.get(&p) {
                    if let Some(tref) = entry.map.get(tn) {
                        tref_cnt = tref.cnt;
                        tref_exists = true;
                    }
                }

                if tref_exists && vref_cnt < tref_cnt {
                    let tv = match non_insertion_variants
                        .get_mut(&p)
                        .and_then(|entry| entry.map.swap_remove(vn))
                    {
                        Some(tv) => tv,
                        None => continue,
                    };

                    let tref = get_variation(non_insertion_variants, p, tn);
                    crate::variations::util::adj_cnt(tref, &tv, None);

                    if let Some(entry) = non_insertion_variants.get(&p) {
                        if entry.map.is_empty() {
                            non_insertion_variants.remove(&p);
                        }
                    }
                }
            }
        }
    }
}

/// Realign insertions.
///
/// Ports `VariationRealigner.realignins`.
pub fn realignins(
    _bams: Option<&[String]>,
    position_to_insertion_count: &BTreeMap<i32, rustc_hash::FxHashMap<String, i32>>,
    reference: &Reference,
    chr: &str,
    chr_len: u32,
    config: &crate::config::Config,
    insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut rustc_hash::FxHashMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
) -> String {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_sets(&config.adaptor);
    let tmp = fill_and_sort_tmp(position_to_insertion_count);
    let mut new_ins = String::new();
    let mut tn_buf = String::new();
    let mut wupseq = String::new();
    let mut sanpseq = String::new();

    for tpl in &tmp {
        let p = tpl.position;
        let vn = &tpl.description_string;
        let insertion_count = tpl.count;

        if config.verbose {
            crate::verbose_eprintln!("  Realign Ins: {} {} {}", p, vn, insertion_count);
        }

        let Some(insert) = parse_begin_plus_atgc(vn) else {
            continue;
        };
        let insert_len_i32 = i32::try_from(insert.len()).unwrap_or(i32::MAX);

        let mut ins3 = "";
        let mut inslen = insert_len_i32;
        if let Some((dup_num_str, ins3_match)) = parse_dup_num_atgc(vn) {
            ins3 = ins3_match;
            let dup_len = dup_num_str.parse::<i32>().unwrap_or(0);
            inslen = inslen
                .saturating_add(dup_len)
                .saturating_add(i32::try_from(ins3.len()).unwrap_or(i32::MAX));
        }

        let extra = parse_amp_atgc(vn).unwrap_or("");
        let compm = parse_hash_atgc(vn).unwrap_or("");
        let newdel = parse_up_number_end(vn).map_or(0, |s| s.parse::<i32>().unwrap_or(0));

        normalize_realignins_tn_into(&mut tn_buf, vn);

        let chr_len_i32 = i32::try_from(chr_len).unwrap_or(i32::MAX);
        let wustart = p.saturating_sub(150).max(1);
        crate::util::ref_bases_into(&mut wupseq, reference, chr, wustart, p);
        wupseq.push_str(&tn_buf);

        let sanend = p
            .saturating_add(i32::try_from(vn.len()).unwrap_or(i32::MAX))
            .saturating_add(100)
            .min(chr_len_i32);

        sanpseq.clear();
        let find_mm3_res = if !ins3.is_empty() {
            let sv_flanking = 50i32;
            let p3_64 =
                i64::from(p) + i64::from(inslen) - (ins3.len() as i64) + i64::from(sv_flanking);
            if ins3.len() > sv_flanking as usize {
                // Ports VariationRealigner.java:744-746:
                // `sanpseq = substr(ins3, Configuration.SVFLANK - ins3.length());`
                // Because `SVFLANK - ins3.length() < 0`, `Utils.substr(ins3, neg)`
                // starts at `ins3.length() + (SVFLANK - ins3.length()) == SVFLANK`.
                let start_idx = sv_flanking as usize;
                sanpseq.push_str(ins3.get(start_idx..).unwrap_or(""));
            }
            if let Some(p_plus_1) = p.checked_add(1) {
                crate::util::push_ref_bases(
                    &mut sanpseq,
                    reference,
                    chr,
                    p_plus_1,
                    p.saturating_add(101),
                );
            }
            if let Ok(p3_plus_1) = i32::try_from(p3_64 + 1) {
                find_mm3(
                    reference,
                    chr,
                    p3_plus_1,
                    sanpseq.as_bytes(),
                    soft_clips_3_end,
                )
            } else {
                MismatchResult::default()
            }
        } else {
            let start_pos_ref64 =
                i64::from(p) + (extra.len() as i64) + 1 + (compm.len() as i64) + i64::from(newdel);
            sanpseq.push_str(&tn_buf);
            if start_pos_ref64 <= i64::from(sanend) {
                if let Ok(start_pos_ref) = i32::try_from(start_pos_ref64) {
                    crate::util::push_ref_bases(
                        &mut sanpseq,
                        reference,
                        chr,
                        start_pos_ref,
                        sanend,
                    );
                }
            }
            if let Some(p_plus_1) = p.checked_add(1) {
                find_mm3(
                    reference,
                    chr,
                    p_plus_1,
                    sanpseq.as_bytes(),
                    soft_clips_3_end,
                )
            } else {
                MismatchResult::default()
            }
        };

        let p5_pos64 =
            i64::from(p) + (extra.len() as i64) + (compm.len() as i64) + i64::from(newdel);
        let find_mm5_res = if let Ok(p5_pos) = i32::try_from(p5_pos64) {
            find_mm5(reference, chr, p5_pos, wupseq.as_bytes(), soft_clips_5_end)
        } else {
            MismatchResult::default()
        };

        let MismatchResult {
            mismatches: mm3,
            scp: ref sc3p,
            nm: nm3,
            misp: misp3,
            misnt: ref misnt3,
        } = find_mm3_res;
        let MismatchResult {
            mismatches: mm5,
            scp: ref sc5p,
            nm: nm5,
            misp: misp5,
            misnt: ref misnt5,
        } = find_mm5_res;
        let mm3_len_is_one = mm3.len() == 1;
        let mm5_len_is_one = mm5.len() == 1;

        let mut mmm = mm3;
        mmm.extend(mm5);

        for mismatch in mmm {
            let mut mismatch_bases = mismatch.mismatch_sequence;
            let mismatch_position = mismatch.mismatch_position;
            let mismatch_end = mismatch.end;

            if let (Some(head), Some(tail)) = (mismatch_bases.get(..1), mismatch_bases.get(1..)) {
                if !tail.is_empty() {
                    mismatch_bases = format!("{head}&{tail}");
                }
            }

            let mut should_adjust = false;
            if let Some(entry) = non_insertion_variants.get(&mismatch_position) {
                if let Some(variation) = entry.map.get(&mismatch_bases) {
                    if variation.cnt > 0 {
                        let mean_qual = variation.qmean / (variation.cnt as f64);
                        let mean_pos = variation.pmean / (variation.cnt as f64);
                        let threshold_pos = if mismatch_end == 3 {
                            nm3.saturating_add(4)
                        } else {
                            nm5.saturating_add(4)
                        };
                        if mean_qual >= config.base_qual
                            && mean_pos <= (threshold_pos as f64)
                            && (variation.cnt < insertion_count.saturating_add(insert_len_i32)
                                && (variation.cnt as f64) / (insertion_count as f64) < 8.0)
                        {
                            should_adjust = true;
                        }
                    }
                }
            }

            if !should_adjust {
                continue;
            }

            let variation = match non_insertion_variants
                .get_mut(&mismatch_position)
                .and_then(|entry| entry.map.swap_remove(&mismatch_bases))
            {
                Some(variation) => variation,
                None => continue,
            };

            if config.verbose {
                let cov = ref_coverage.get(&p).copied().unwrap_or(0);
                crate::verbose_eprintln!(
                    "    insMM: {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    mismatch_bases,
                    mismatch_position,
                    mismatch_end,
                    nm3,
                    nm5,
                    vn,
                    insertion_count,
                    variation.cnt,
                    variation.qmean,
                    variation.pmean,
                    cov
                );
            }

            if mismatch_position > p && mismatch_end == 5 {
                let cov_entry = ref_coverage.entry(p).or_default();
                *cov_entry = cov_entry.saturating_add(variation.cnt);
            }

            let mut lref_val = None;
            if mismatch_position > p && mismatch_end == 3 {
                if let Some(ref_base) = reference.get_base(chr, p as u32) {
                    let mut ch_buf = [0u8; 4];
                    let ref_base_str: &str = (ref_base as char).encode_utf8(&mut ch_buf);
                    if let Some(entry) = non_insertion_variants.get_mut(&p) {
                        if let Some(lref) = entry.map.get_mut(ref_base_str) {
                            lref_val = Some(lref);
                        }
                    }
                }
            }

            let vref_mut = get_variation(insertion_variants, p, vn);
            crate::variations::util::adj_cnt(vref_mut, &variation, lref_val);

            if let Some(entry) = non_insertion_variants.get_mut(&mismatch_position) {
                if entry.map.is_empty() {
                    non_insertion_variants.remove(&mismatch_position);
                }
            }
        }

        if misp3 != 0 && mm3_len_is_one {
            let mut remove_misp3 = false;
            if let Some(entry) = non_insertion_variants.get(&misp3) {
                if let Some(var) = entry.map.get(misnt3) {
                    if var.cnt < insertion_count {
                        remove_misp3 = true;
                    }
                }
            }
            if remove_misp3 {
                if let Some(entry) = non_insertion_variants.get_mut(&misp3) {
                    entry.map.swap_remove(misnt3);
                    if entry.map.is_empty() {
                        non_insertion_variants.remove(&misp3);
                    }
                }
            }
        }

        if misp5 != 0 && mm5_len_is_one {
            let mut remove_misp5 = false;
            if let Some(entry) = non_insertion_variants.get(&misp5) {
                if let Some(var) = entry.map.get(misnt5) {
                    if var.cnt < insertion_count {
                        remove_misp5 = true;
                    }
                }
            }
            if remove_misp5 {
                if let Some(entry) = non_insertion_variants.get_mut(&misp5) {
                    entry.map.swap_remove(misnt5);
                    if entry.map.is_empty() {
                        non_insertion_variants.remove(&misp5);
                    }
                }
            }
        }

        for &sc5pp in sc5p {
            if let Some(tv) = soft_clips_5_end.get_mut(&sc5pp) {
                if config.verbose {
                    crate::verbose_eprintln!(
                        "    55: {} {} VN: '{}'  5' seq: ^{}^\n",
                        p,
                        sc5pp,
                        vn,
                        wupseq
                    );
                }
                if !tv.used {
                    let seq = crate::variations::util::find_conseq_sets(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if config.verbose {
                        crate::verbose_eprintln!(
                            "    ins5: {} {} {} {} VN: {} iCnt: {} vCnt: {}",
                            p,
                            sc5pp,
                            String::from_utf8_lossy(&seq),
                            wupseq,
                            vn,
                            insertion_count,
                            tv.variation.cnt
                        );
                    }
                    if !seq.is_empty() && is_match(&seq, wupseq.as_bytes(), -1, 3) {
                        if config.verbose {
                            crate::verbose_eprintln!(
                                "      ins5: {} {} {} {} VN: {} iCnt: {} cCnt: {} used",
                                p,
                                sc5pp,
                                String::from_utf8_lossy(&seq),
                                wupseq,
                                vn,
                                insertion_count,
                                tv.variation.cnt
                            );
                        }
                        if sc5pp > p {
                            let cov_entry = ref_coverage.entry(p).or_default();
                            *cov_entry = cov_entry.saturating_add(tv.variation.cnt);
                        }
                        let vref_mut = get_variation(insertion_variants, p, vn);
                        crate::variations::util::adj_cnt(vref_mut, &tv.variation, None);

                        if let Some(tv_mut) = soft_clips_5_end.get_mut(&sc5pp) {
                            tv_mut.used = true;
                        }
                    }
                }
            }
        }

        for &sc3pp in sc3p {
            if let Some(tv) = soft_clips_3_end.get_mut(&sc3pp) {
                if config.verbose {
                    crate::verbose_eprintln!(
                        "    33: {} {} VN: '{}'  3' seq: ^{}^\n",
                        p,
                        sc3pp,
                        vn,
                        sanpseq
                    );
                }
                if !tv.used {
                    let seq = crate::variations::util::find_conseq_sets(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if config.verbose {
                        crate::verbose_eprintln!(
                            "    ins3: {} {} {} {} VN: {} iCnt: {} vCnt: {}",
                            p,
                            sc3pp,
                            String::from_utf8_lossy(&seq),
                            sanpseq,
                            vn,
                            insertion_count,
                            tv.variation.cnt
                        );
                    }
                    let mseq = if !ins3.is_empty() {
                        sanpseq.as_bytes()
                    } else {
                        substr_from(
                            sanpseq.as_bytes(),
                            sc3pp.saturating_sub(p).saturating_sub(1) as isize,
                        )
                    };

                    if !seq.is_empty() && is_match(&seq, mseq, 1, 3) {
                        if config.verbose {
                            crate::verbose_eprintln!(
                                "      ins3: {} {} {} VN: {} iCnt: {} vCnt: {} used",
                                p,
                                sc3pp,
                                String::from_utf8_lossy(&seq),
                                vn,
                                insertion_count,
                                tv.variation.cnt
                            );
                        }
                        if sc3pp <= p
                            || (insert.len() as f64)
                                > (tv.variation.pmean / (tv.variation.cnt as f64))
                        {
                            let cov_entry = ref_coverage.entry(p).or_default();
                            *cov_entry = cov_entry.saturating_add(tv.variation.cnt);
                        }

                        let mut lref_val = None;
                        let mut insert_gt_mean_pos = false;
                        if tv.variation.cnt > 0
                            && (insert.len() as f64)
                                > (tv.variation.pmean / (tv.variation.cnt as f64))
                        {
                            insert_gt_mean_pos = true;
                        }

                        if !insert_gt_mean_pos && sc3pp > p {
                            if let Some(ref_base) = reference.get_base(chr, p as u32) {
                                let mut ch_buf = [0u8; 4];
                                let ref_base_str: &str =
                                    (ref_base as char).encode_utf8(&mut ch_buf);
                                if let Some(entry) = non_insertion_variants.get_mut(&p) {
                                    if let Some(lref) = entry.map.get_mut(ref_base_str) {
                                        lref_val = Some(lref);
                                    }
                                }
                            }
                        }

                        let vref_mut = get_variation(insertion_variants, p, vn);
                        crate::variations::util::adj_cnt(vref_mut, &tv.variation, lref_val);

                        if let Some(tv_mut) = soft_clips_3_end.get_mut(&sc3pp) {
                            tv_mut.used = true;
                        }

                        if insert.len() + 1 == vn.len()
                            && insert.len() > max_read_length as usize
                            && i64::from(sc3pp) >= i64::from(p) + 1 + (insert.len() as i64)
                        {
                            let mut flag = 0;
                            let offset =
                                sc3pp.saturating_sub(p).saturating_sub(1) as usize % insert.len();
                            let mut tvn = vn.clone();
                            let insert_bytes = insert.as_bytes();

                            for seqi in 0..seq.len() {
                                if seqi + offset < insert_bytes.len() {
                                    let seq_char = seq[seqi] as char;
                                    let ins_char = insert_bytes[seqi + offset] as char;
                                    if seq_char != ins_char {
                                        flag += 1;
                                        let shift = seqi + offset + 1;
                                        if let (Some(prefix), Some(suffix)) =
                                            (tvn.get(..shift), tvn.get(shift + 1..))
                                        {
                                            tvn = format!("{prefix}{seq_char}{suffix}");
                                        }
                                    }
                                }
                            }
                            if flag > 0 {
                                if let Some(entry) = insertion_variants.get_mut(&p) {
                                    if let Some(variation) = entry.map.swap_remove(vn) {
                                        entry.map.insert(tvn.clone(), variation);
                                    }
                                }
                                new_ins = tvn;
                            }
                        }
                    }
                }
            }
        }

        if !sc3p.is_empty() && !sc5p.is_empty() {
            let first3 = sc3p[0];
            let first5 = sc5p[0];
            let diff35 = i64::from(first3) - i64::from(first5);
            // ports VariationRealigner.java:921: the bound stays a double (no truncation).
            if first3 > first5.saturating_add(3)
                && (diff35 as f64) < (max_read_length as f64) * 0.75
            {
                if let Some(ref_base) = reference.get_base(chr, p as u32) {
                    let mut ch_buf = [0u8; 4];
                    let ref_base_str: &str = (ref_base as char).encode_utf8(&mut ch_buf);
                    if let Some(entry) = non_insertion_variants.get_mut(&p) {
                        if let Some(ref_var) = entry.map.get_mut(ref_base_str) {
                            crate::variations::util::adj_ref_factor(
                                ref_var,
                                ((diff35 - 1) as f64) / (max_read_length as f64),
                            );
                        }
                    }
                }
                let vref_mut = get_variation(insertion_variants, p, vn);
                crate::variations::util::adj_ref_factor(
                    vref_mut,
                    -((diff35 - 1) as f64) / (max_read_length as f64),
                );
            }
        }
    }

    for i in (1..tmp.len()).rev() {
        let tpl = &tmp[i];
        let p = tpl.position;
        let vn = &tpl.description_string;

        if !insertion_variants.contains_key(&p) {
            continue;
        }

        let mut vref_cnt = 0;
        let mut vref_exists = false;
        if let Some(entry) = insertion_variants.get(&p) {
            if let Some(vref) = entry.map.get(vn) {
                vref_cnt = vref.cnt;
                vref_exists = true;
            }
        }

        if !vref_exists {
            continue;
        }

        if let Some(tn) = parse_atgs_amp_atgs_end(vn) {
            let mut tref_cnt = 0;
            let mut tref_exists = false;
            if let Some(entry) = insertion_variants.get(&p) {
                if let Some(tref) = entry.map.get(tn) {
                    tref_cnt = tref.cnt;
                    tref_exists = true;
                }
            }

            if tref_exists && vref_cnt < tref_cnt {
                let tv = match insertion_variants
                    .get_mut(&p)
                    .and_then(|entry| entry.map.swap_remove(vn))
                {
                    Some(tv) => tv,
                    None => continue,
                };

                let mut ref_base_var = None;
                if let Some(ref_base) = reference.get_base(chr, p as u32) {
                    let mut ch_buf = [0u8; 4];
                    let ref_base_str: &str = (ref_base as char).encode_utf8(&mut ch_buf);
                    if let Some(entry) = non_insertion_variants.get_mut(&p) {
                        if let Some(lref) = entry.map.get_mut(ref_base_str) {
                            ref_base_var = Some(lref);
                        }
                    }
                }

                let tref = get_variation(insertion_variants, p, tn);
                crate::variations::util::adj_cnt(tref, &tv, ref_base_var);

                if let Some(entry) = insertion_variants.get_mut(&p) {
                    if entry.map.is_empty() {
                        insertion_variants.remove(&p);
                    }
                }
            }
        }
    }

    new_ins
}

/// Ports `VariationRealigner.realignlgdel`.
pub fn realignlgdel(
    bams: Option<&[String]>,
    reference: &mut Reference,
    chr: &str,
    chr_len: u32,
    region: &crate::regions::Region,
    config: &crate::config::Config,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut rustc_hash::FxHashMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
    svfdel: &mut Vec<Sclip>,
    svrdel: &mut Vec<Sclip>,
    splice: &mut rustc_hash::FxHashSet<String>,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_sets(&config.adaptor);
    let min_reads = i32::try_from(config.min_reads).unwrap_or(i32::MAX);
    let indel_size = i32::try_from(config.indel_size).unwrap_or(i32::MAX);
    let extension = 5000;
    let longmm = 3;
    let region_start_i32 = i32::try_from(region.start.0).unwrap_or(i32::MAX);
    let region_end_i32 = i32::try_from(region.end.0).unwrap_or(i32::MAX);
    let min_clip_pos = region_start_i32.saturating_sub(extension);
    let max_clip_pos = region_end_i32.saturating_add(extension);

    // Work on 5' clipped reads

    let mut tmp5 = Vec::new();
    for (&p, sc) in soft_clips_5_end.iter() {
        if p < min_clip_pos || p > max_clip_pos {
            continue;
        }
        tmp5.push((p, sc.variation.cnt));
    }
    tmp5.sort_by(|a, b| {
        let cmp_cnt = b.1.cmp(&a.1);
        if cmp_cnt != std::cmp::Ordering::Equal {
            cmp_cnt
        } else {
            a.0.cmp(&b.0)
        }
    });

    let mut svcov = 0;

    for (p, _) in tmp5 {
        let (cnt, seq, sc5_variation) = {
            let Some(sc) = soft_clips_5_end.get_mut(&p) else {
                continue;
            };
            let cnt = sc.variation.cnt;
            if cnt < min_reads {
                break;
            }
            if sc.used {
                continue;
            }
            let seq = crate::variations::util::find_conseq_sets(
                sc,
                5,
                &adaptor_forward,
                &adaptor_reverse,
            );
            if seq.len() < 7 {
                continue;
            }
            (cnt, seq, sc.variation.clone())
        };

        if config.verbose {
            crate::verbose_eprintln!(
                "  Working Realignlgdel: 5' {} '{}' {}",
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }

        let mut bp = find_bp(
            &seq,
            p.saturating_sub(5),
            reference,
            -1,
            chr,
            indel_size,
            chr_len,
        );
        let mut extra = Vec::new();
        let mut extra_str = Vec::new();

        if bp == 0 {
            if crate::variations::util::is_low_complex_seq(&seq) {
                continue;
            }
            let match_res =
                crate::sv::find_match(&seq, reference, p, -1, 17, 1, chr, config.verbose);
            bp = match_res.base_position;
            extra_str = match_res.matched_sequence;
            let p_minus_bp = i64::from(p) - i64::from(bp);
            if !(bp != 0 && p_minus_bp > 15 && p_minus_bp < 150000) {
                continue;
            }
            bp = bp.saturating_add(1);

            let mut lists = [&mut *svfdel, &mut *svrdel];
            let sv_mark = crate::sv::mark_sv(bp, p, &mut lists, max_read_length, config.verbose);
            svcov = sv_mark.0;
            let clusters = sv_mark.1;
            let pairs = sv_mark.2;

            if svcov == 0 && cnt <= min_reads {
                continue;
            }

            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, bp);
            sv.r#type = "DEL".to_string();
            sv.pairs = sv.pairs.saturating_add(pairs);
            sv.splits = sv.splits.saturating_add(cnt);
            sv.clusters = sv.clusters.saturating_add(clusters);

            if bp < region_start_i32 {
                let tts = bp.saturating_sub(max_read_length);
                let mut tte = bp.saturating_add(max_read_length);
                if tte >= region_start_i32 {
                    tte = region_start_i32.saturating_sub(1);
                }
                if let Err(e) = ensure_reference_loaded(
                    reference,
                    chr,
                    tts,
                    tte,
                    max_read_length,
                    chr_len,
                    config,
                ) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    region,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                    splice,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }
        }

        let mut dellen = p.saturating_sub(bp);
        let mut en = 0;
        let gt;

        if extra_str.is_empty() {
            while en < seq.len() {
                let ref_ch = i32::try_from(en)
                    .ok()
                    .and_then(|ei| bp.checked_sub(ei))
                    .and_then(|p| p.checked_sub(1))
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                if ref_ch.is_some() && ref_ch != Some(seq[en]) {
                    extra.push(seq[en]);
                    en += 1;
                } else {
                    break;
                }
            }
            if !extra.is_empty() {
                extra.reverse();
                let extra_str_val = String::from_utf8_lossy(&extra).into_owned();
                gt = format!("-{}&{}", dellen, extra_str_val).into_bytes();
                bp = bp.saturating_sub(extra.len() as i32);
            } else {
                gt = format!("-{}", dellen).into_bytes();
            }
        } else {
            dellen = dellen.saturating_sub(extra_str.len() as i32);
            if dellen == 0 {
                gt = format!(
                    "-{}^{}",
                    extra_str.len(),
                    String::from_utf8_lossy(&extra_str)
                )
                .into_bytes();
            } else {
                gt = format!("-{}&{}", dellen, String::from_utf8_lossy(&extra_str)).into_bytes();
            }
        }

        if config.verbose {
            crate::verbose_eprintln!(
                "  Found Realignlgdel: {} {} 5' {} '{}' {}",
                bp,
                String::from_utf8_lossy(&gt),
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }

        let mut n: i32 = 0;
        if extra.is_empty() && extra_str.is_empty() {
            loop {
                let base_bp = bp
                    .checked_add(n)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                let base_bp_del = bp
                    .checked_add(dellen)
                    .and_then(|p| p.checked_add(n))
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                if base_bp.is_some() && base_bp_del.is_some() && base_bp == base_bp_del {
                    if n == i32::MAX {
                        break;
                    }
                    n = n.saturating_add(1);
                } else {
                    break;
                }
            }
        }

        let mut sc3p = bp.saturating_add(n);
        let mut str_seq = Vec::new();
        let mut mcnt = 0;
        while mcnt <= longmm {
            let base_bp = bp
                .checked_add(n)
                .and_then(|p| u32::try_from(p).ok())
                .and_then(|pos| reference.get_base(chr, pos));
            let base_bp_del = bp
                .checked_add(dellen)
                .and_then(|p| p.checked_add(n))
                .and_then(|p| u32::try_from(p).ok())
                .and_then(|pos| reference.get_base(chr, pos));
            if let (Some(b_bp), Some(b_bp_del)) = (base_bp, base_bp_del) {
                if b_bp != b_bp_del {
                    str_seq.push(b_bp_del);
                    mcnt += 1;
                    if n == i32::MAX {
                        break;
                    }
                    n = n.saturating_add(1);
                    continue;
                }
            }
            break;
        }

        if str_seq.len() == 1 {
            let mut nm: i32 = 0;
            loop {
                let base_bp = bp
                    .checked_add(n)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                let base_bp_del = bp
                    .checked_add(dellen)
                    .and_then(|p| p.checked_add(n))
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                if base_bp.is_some() && base_bp_del.is_some() && base_bp == base_bp_del {
                    if n == i32::MAX {
                        break;
                    }
                    n = n.saturating_add(1);
                    if n != 0 {
                        nm = nm.saturating_add(1);
                    }
                } else {
                    break;
                }
            }
            if nm >= 3 && !soft_clips_3_end.contains_key(&sc3p) {
                sc3p = bp.saturating_add(n);
            }
        }

        let mut should_continue = false;
        if let Some(entry) = non_insertion_variants.get(&bp) {
            if entry.sv.is_some()
                && !soft_clips_3_end.contains_key(&sc3p)
                && svcov == 0
                && cnt <= min_reads
            {
                crate::variations::VariationMap::remove_sv(non_insertion_variants, bp);
                should_continue = true;
            }
        }
        if should_continue {
            continue;
        }

        let gt_str = String::from_utf8_lossy(&gt).into_owned();
        {
            let tv = crate::variations::get_variation(non_insertion_variants, bp, &gt_str);
            tv.qstd = true;
            tv.pstd = true;
            crate::variations::util::adj_cnt(tv, &sc5_variation, None);
        }

        if let Some(sc) = soft_clips_5_end.get_mut(&p) {
            sc.used = bp != 0;
        }

        if !ref_coverage.contains_key(&bp) {
            if let Some(&cov_p) = ref_coverage.get(&p) {
                ref_coverage.insert(bp, cov_p);
            }
        }

        if dellen < indel_size {
            for tp in bp..bp.saturating_add(dellen) {
                inc_cnt(ref_coverage, tp, sc5_variation.cnt);
            }
        }

        if let Some(sclip) = soft_clips_3_end.get_mut(&sc3p) {
            if !sclip.used {
                {
                    let tv = crate::variations::get_variation(non_insertion_variants, bp, &gt_str);
                    crate::variations::util::adj_cnt(tv, &sclip.variation, None);
                }
                if sc3p > bp {
                    let ref_base_ch = reference.get_base(chr, bp as u32);
                    if let Some(ref_var) = crate::variations::get_variation_maybe_mut(
                        non_insertion_variants,
                        bp,
                        ref_base_ch,
                    ) {
                        ref_var.cnt = ref_var.cnt.saturating_sub(sclip.variation.cnt);
                        ref_var.hicnt = ref_var.hicnt.saturating_sub(sclip.variation.hicnt);
                        ref_var.locnt = ref_var.locnt.saturating_sub(sclip.variation.locnt);
                        ref_var.pmean -= sclip.variation.pmean;
                        ref_var.qmean -= sclip.variation.qmean;
                        ref_var.mqmean -= sclip.variation.mqmean;
                        ref_var.nm -= sclip.variation.nm;
                        ref_var.sub_dir(true, sclip.variation.get_dir(true));
                        ref_var.sub_dir(false, sclip.variation.get_dir(false));
                        crate::variations::util::correct_cnt(ref_var);
                    }
                }

                if sc3p == bp && dellen < indel_size {
                    for tp in bp..bp.saturating_add(dellen) {
                        inc_cnt(ref_coverage, tp, sclip.variation.cnt);
                    }
                }

                for ip in bp.saturating_add(1)..sc3p {
                    let ref_ch_opt = dellen
                        .checked_add(ip)
                        .and_then(|p| u32::try_from(p).ok())
                        .and_then(|pos| reference.get_base(chr, pos));
                    let mut cnt_is_zero = false;
                    if let Some(vv) = crate::variations::get_variation_maybe_mut(
                        non_insertion_variants,
                        ip,
                        ref_ch_opt,
                    ) {
                        rm_cnt(vv, &sclip.variation);
                        if vv.cnt == 0 {
                            cnt_is_zero = true;
                        }
                    }
                    if cnt_is_zero {
                        if let Some(ref_ch) = ref_ch_opt {
                            let mut ch_buf = [0u8; 4];
                            let ref_ch_str: &str = (ref_ch as char).encode_utf8(&mut ch_buf);
                            if let Some(vm) = non_insertion_variants.get_mut(&ip) {
                                vm.map.shift_remove(ref_ch_str);
                            }
                        }
                    }
                    if let Some(vm) = non_insertion_variants.get(&ip) {
                        if vm.map.is_empty() {
                            non_insertion_variants.remove(&ip);
                        }
                    }
                }

                sclip.used = bp != 0;
            }
        }

        let mut dels5 = BTreeMap::new();
        let mut inner_map = rustc_hash::FxHashMap::default();
        let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
        inner_map.insert(gt_str.clone(), tv_cnt);
        dels5.insert(bp, inner_map);

        realigndel(
            bams,
            &dels5,
            reference,
            chr,
            chr_len,
            config,
            non_insertion_variants,
            ref_coverage,
            soft_clips_5_end,
            soft_clips_3_end,
            max_read_length,
        );

        if let Some(entry) = non_insertion_variants.get_mut(&bp) {
            if let Some(ref mut sv) = entry.sv {
                if let Some(inner) = dels5.get(&bp) {
                    if let Some(&dels5_cnt) = inner.get(&gt_str) {
                        let tv_cnt = entry.map.get(&gt_str).map(|v| v.cnt).unwrap_or(0);
                        sv.splits = sv.splits.saturating_add(tv_cnt.saturating_sub(dels5_cnt));
                    }
                }
            }
        }

        let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
        if svcov > tv_cnt {
            let tv = crate::variations::get_variation(non_insertion_variants, bp, &gt_str);
            add_var_factor(
                tv,
                (i64::from(svcov) - i64::from(tv_cnt)) as f64 / tv_cnt as f64,
            );
        }

        if config.verbose {
            let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
            crate::verbose_eprintln!(
                "  Found lgdel done: {} {} {} 5' {} {}\n",
                bp,
                gt_str,
                p,
                String::from_utf8_lossy(&seq),
                tv_cnt
            );
        }
    }

    // Work on 3' clipped reads
    let mut tmp3 = Vec::new();
    for (&p, sc) in soft_clips_3_end.iter() {
        if p < min_clip_pos || p > max_clip_pos {
            continue;
        }
        tmp3.push((p, sc.variation.cnt));
    }
    tmp3.sort_by(|a, b| {
        let cmp_cnt = b.1.cmp(&a.1);
        if cmp_cnt != std::cmp::Ordering::Equal {
            cmp_cnt
        } else {
            a.0.cmp(&b.0)
        }
    });

    for (p, _) in tmp3 {
        let (cnt, seq) = {
            let Some(sc) = soft_clips_3_end.get_mut(&p) else {
                continue;
            };
            let cnt = sc.variation.cnt;
            if cnt < min_reads {
                break;
            }
            if sc.used {
                continue;
            }
            let seq = crate::variations::util::find_conseq_sets(
                sc,
                3,
                &adaptor_forward,
                &adaptor_reverse,
            );
            if seq.len() < 7 {
                continue;
            }
            (cnt, seq)
        };
        let Some(mut sc3v) = soft_clips_3_end.get(&p).cloned() else {
            continue;
        };

        if config.verbose {
            crate::verbose_eprintln!(
                "  Working Realignlgdel: 3' {} '{}' {}",
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }

        let mut bp = find_bp(
            &seq,
            p.saturating_add(5),
            reference,
            1,
            chr,
            indel_size,
            chr_len,
        );
        let mut extra = Vec::new();
        let mut extra_str = Vec::new();

        if bp == 0 {
            if crate::variations::util::is_low_complex_seq(&seq) {
                continue;
            }
            let match_res =
                crate::sv::find_match(&seq, reference, p, 1, 17, 1, chr, config.verbose);
            bp = match_res.base_position;
            extra_str = match_res.matched_sequence;
            let bp_minus_p = i64::from(bp) - i64::from(p);
            let p_minus_bp = i64::from(p) - i64::from(bp);
            if !(bp != 0 && bp_minus_p > 15 && p_minus_bp < 150000) {
                continue;
            }

            let mut lists = [&mut *svfdel, &mut *svrdel];
            let sv_mark = crate::sv::mark_sv(p, bp, &mut lists, max_read_length, config.verbose);
            svcov = sv_mark.0;
            let clusters = sv_mark.1;
            let pairs = sv_mark.2;

            if svcov == 0 && cnt <= min_reads {
                continue;
            }

            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, p);
            sv.r#type = "DEL".to_string();
            sv.pairs = sv.pairs.saturating_add(pairs);
            sv.splits = sv.splits.saturating_add(cnt);
            sv.clusters = sv.clusters.saturating_add(clusters);

            if bp > region_end_i32 {
                let mut tts = bp.saturating_sub(max_read_length);
                let tte = bp.saturating_add(max_read_length);
                if tts <= region_end_i32 {
                    tts = region_end_i32.saturating_add(1);
                }
                if let Err(e) = ensure_reference_loaded(
                    reference,
                    chr,
                    tts,
                    tte,
                    max_read_length,
                    chr_len,
                    config,
                ) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    region,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                    splice,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }
        }

        let mut dellen = bp.saturating_sub(p);
        let mut en = 0;
        if !extra_str.is_empty() {
            dellen = dellen.saturating_sub(extra_str.len() as i32);
        } else {
            while en < seq.len() {
                let ref_ch = i32::try_from(en)
                    .ok()
                    .and_then(|ei| bp.checked_add(ei))
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                if ref_ch.is_some() && ref_ch != Some(seq[en]) {
                    extra.push(seq[en]);
                    en += 1;
                } else {
                    break;
                }
            }
        }

        let mut gt = format!("-{}", dellen).into_bytes();
        let mut sc5p = bp;
        bp = p;

        if !extra.is_empty() {
            let extra_str_val = String::from_utf8_lossy(&extra).into_owned();
            gt = format!("-{}&{}", dellen, extra_str_val).into_bytes();
            sc5p = sc5p.saturating_add(extra.len() as i32);
        } else if !extra_str.is_empty() {
            let extra_str_val = String::from_utf8_lossy(&extra_str).into_owned();
            gt = format!("-{}&{}", dellen, extra_str_val).into_bytes();
        } else {
            loop {
                let base_bp_prev = bp
                    .checked_sub(1)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                let base_bp_del_prev = bp
                    .checked_add(dellen)
                    .and_then(|p| p.checked_sub(1))
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                if base_bp_prev.is_some()
                    && base_bp_del_prev.is_some()
                    && base_bp_prev == base_bp_del_prev
                {
                    let Some(next_bp) = bp.checked_sub(1) else {
                        break;
                    };
                    bp = next_bp;
                    if bp != 0 {
                        sc5p = sc5p.saturating_sub(1);
                    }
                } else {
                    break;
                }
            }

            if bp != p {
                let mut move_sv = None;
                if let Some(entry) = non_insertion_variants.get_mut(&p) {
                    if let Some(sv) = entry.sv.take() {
                        move_sv = Some(sv);
                        entry.map.shift_remove("SV");
                    }
                }
                if let Some(sv) = move_sv {
                    let dest = non_insertion_variants.get_or_put_with(bp, VariationMap::new);
                    dest.sv = Some(sv);
                    dest.map
                        .insert("SV".to_string(), crate::variations::Variation::new());
                }
            }
        }

        let mut should_continue = false;
        if let Some(entry) = non_insertion_variants.get(&bp) {
            if entry.sv.is_some()
                && !soft_clips_5_end.contains_key(&sc5p)
                && svcov == 0
                && cnt <= min_reads
            {
                crate::variations::VariationMap::remove_sv(non_insertion_variants, bp);
                should_continue = true;
            }
        }
        if should_continue {
            continue;
        }

        let gt_str = String::from_utf8_lossy(&gt).into_owned();
        if config.verbose {
            crate::verbose_eprintln!(
                "  Found Realignlgdel: bp: {} {} 3' {} 5'clip: {} '{}' {}",
                bp,
                gt_str,
                p,
                sc5p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }

        {
            let tv = crate::variations::get_variation(non_insertion_variants, bp, &gt_str);
            tv.qstd = true;
            tv.pstd = true;
        }

        if dellen < indel_size {
            let tp_end = bp
                .saturating_add(dellen)
                .saturating_add(extra.len() as i32)
                .saturating_add(extra_str.len() as i32);
            for tp in bp..tp_end {
                inc_cnt(ref_coverage, tp, sc3v.variation.cnt);
            }
        }

        if !ref_coverage.contains_key(&bp) {
            if let Some(&cov_p_prev) = ref_coverage.get(&(p.saturating_sub(1))) {
                ref_coverage.insert(bp, cov_p_prev);
            } else {
                ref_coverage.insert(bp, sc3v.variation.cnt);
            }
        }

        sc3v.variation.pmean += (dellen as f64) * (sc3v.variation.cnt as f64);
        {
            let tv = crate::variations::get_variation(non_insertion_variants, bp, &gt_str);
            crate::variations::util::adj_cnt(tv, &sc3v.variation, None);
        }
        sc3v.used = p.saturating_add(dellen) != 0;

        if let Some(sc) = soft_clips_3_end.get_mut(&p) {
            sc.used = p.saturating_add(dellen) != 0;
            sc.sequence = sc3v.sequence;
            sc.variation.pmean = sc3v.variation.pmean;
        }

        let mut dels5 = BTreeMap::new();
        let mut map = rustc_hash::FxHashMap::default();
        let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
        map.insert(gt_str.clone(), tv_cnt);
        dels5.insert(bp, map);

        realigndel(
            bams,
            &dels5,
            reference,
            chr,
            chr_len,
            config,
            non_insertion_variants,
            ref_coverage,
            soft_clips_5_end,
            soft_clips_3_end,
            max_read_length,
        );

        if let Some(entry) = non_insertion_variants.get_mut(&bp) {
            if let Some(ref mut sv) = entry.sv {
                if let Some(inner) = dels5.get(&bp) {
                    if let Some(&dels5_cnt) = inner.get(&gt_str) {
                        let tv_cnt = entry.map.get(&gt_str).map(|v| v.cnt).unwrap_or(0);
                        sv.splits = sv.splits.saturating_add(tv_cnt.saturating_sub(dels5_cnt));
                    }
                }
            }
        }

        if config.verbose {
            let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
            crate::verbose_eprintln!(
                "  Found lgdel: {} {} {} 3' '{}' {}\n",
                bp,
                gt_str,
                p,
                String::from_utf8_lossy(&seq),
                tv_cnt
            );
        }

        let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
        if svcov > tv_cnt {
            let tv = crate::variations::get_variation(non_insertion_variants, bp, &gt_str);
            add_var_factor(
                tv,
                (i64::from(svcov) - i64::from(tv_cnt)) as f64 / tv_cnt as f64,
            );
        }
    }
}

/// Ports `VariationRealigner.realignlgins30`.
pub fn realignlgins30(
    bams: Option<&[String]>,
    reference: &mut Reference,
    chr: &str,
    chr_len: u32,
    region: &crate::regions::Region,
    config: &crate::config::Config,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut rustc_hash::FxHashMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_sets(&config.adaptor);
    let extension = 5000;
    let region_start_i32 = i32::try_from(region.start.0).unwrap_or(i32::MAX);
    let region_end_i32 = i32::try_from(region.end.0).unwrap_or(i32::MAX);
    let min_clip_pos = region_start_i32.saturating_sub(extension);
    let max_clip_pos = region_end_i32.saturating_add(extension);

    let mut tmp5 = Vec::new();
    for (&p, sc) in soft_clips_5_end.iter() {
        if p < min_clip_pos || p > max_clip_pos {
            continue;
        }
        tmp5.push((p, sc.variation.cnt));
    }
    tmp5.sort_by(|a, b| {
        let cmp_cnt = b.1.cmp(&a.1);
        if cmp_cnt != std::cmp::Ordering::Equal {
            cmp_cnt
        } else {
            a.0.cmp(&b.0)
        }
    });

    let mut tmp3 = Vec::new();
    for (&p, sc) in soft_clips_3_end.iter() {
        if p < min_clip_pos || p > max_clip_pos {
            continue;
        }
        tmp3.push((p, sc.variation.cnt));
    }
    tmp3.sort_by(|a, b| {
        let cmp_cnt = b.1.cmp(&a.1);
        if cmp_cnt != std::cmp::Ordering::Equal {
            cmp_cnt
        } else {
            a.0.cmp(&b.0)
        }
    });

    let max_dist_53 = (max_read_length as f64 * 2.5) as i32;
    let max_dist_35 = max_read_length.saturating_sub(10);
    for (p5, cnt5) in tmp5 {
        if soft_clips_5_end.get(&p5).is_none_or(|sc| sc.used) {
            continue;
        }
        let mut sc5_used = false;
        for (p3, cnt3) in tmp3.iter().copied() {
            if sc5_used {
                break;
            }
            if i64::from(p5) - i64::from(p3) > i64::from(max_dist_53)
                || i64::from(p3) - i64::from(p5) > i64::from(max_dist_35)
            {
                continue;
            }
            if soft_clips_3_end.get(&p3).is_none_or(|sc| sc.used) {
                continue;
            }
            let seq5 = {
                let Some(sc5v) = soft_clips_5_end.get_mut(&p5) else {
                    break;
                };
                let s = crate::variations::util::find_conseq_sets(
                    sc5v,
                    5,
                    &adaptor_forward,
                    &adaptor_reverse,
                );
                sc5_used = sc5v.used;
                s
            };
            let seq3 = {
                let Some(sc3v) = soft_clips_3_end.get_mut(&p3) else {
                    continue;
                };
                crate::variations::util::find_conseq_sets(
                    sc3v,
                    3,
                    &adaptor_forward,
                    &adaptor_reverse,
                )
            };
            if seq5.len() <= 10 || seq3.len() <= 10 {
                continue;
            }
            if config.verbose {
                crate::verbose_eprintln!(
                    "  Working lgins30: {} {} 3: {} {} 5: {} {}",
                    p3,
                    p5,
                    String::from_utf8_lossy(&seq3),
                    cnt3,
                    String::from_utf8_lossy(&reverse(&seq5)),
                    cnt5
                );
            }
            if !((cnt5 as f64) / (cnt3 as f64) >= 0.08 && (cnt5 as f64) / (cnt3 as f64) <= 12.0) {
                continue;
            }
            let match35 = find_35_match(&seq5, &seq3, config.verbose);
            let bp5 = match35.matched_5_end as isize;
            let bp3 = match35.matched_3_end as isize;
            let score = match35.max_matched_length;
            if score == 0 {
                continue;
            }
            let smscore = (score / 2) as isize;
            let bp3_plus_sm = bp3.saturating_add(smscore);
            let mut ins = if bp3_plus_sm > 1 {
                substr(&seq3, 0, 1_isize.saturating_sub(bp3_plus_sm)).to_vec()
            } else {
                seq3.clone()
            };
            let bp5_plus_sm = bp5.saturating_add(smscore);
            if bp5_plus_sm > 0 {
                ins.extend_from_slice(&reverse(substr(&seq5, 0, bp5_plus_sm)));
            }
            if crate::variations::util::is_low_complex_seq(&ins) {
                if config.verbose {
                    crate::verbose_eprintln!(
                        "  Discard low complexity insertion found {}.",
                        String::from_utf8_lossy(&ins)
                    );
                }
                continue;
            }
            let bi;
            let ins_str;
            if config.verbose {
                crate::verbose_eprintln!(
                    "  Found candidate lgins30: {} {} {}",
                    p3,
                    p5,
                    String::from_utf8_lossy(&ins)
                );
            }
            if p5 > p3 {
                if seq3.len() > ins.len() {
                    let sub = substr_from(&seq3, ins.len() as isize);
                    let jref = join_ref_int(
                        reference,
                        chr,
                        p5,
                        p5.saturating_add((seq3.len() - ins.len()) as i32)
                            .saturating_add(2),
                    );
                    if !is_match(sub, &jref, 1, 3) {
                        continue;
                    }
                }
                if seq5.len() > ins.len() {
                    let sub = substr_from(&seq5, ins.len() as isize);
                    let jref = join_ref_int(
                        reference,
                        chr,
                        p3.saturating_sub((seq5.len() - ins.len()) as i32)
                            .saturating_sub(2),
                        p3.saturating_sub(1),
                    );
                    if !is_match(sub, &jref, -1, 3) {
                        continue;
                    }
                }
                if config.verbose {
                    crate::verbose_eprintln!(
                        "  Found lgins30 complex: {} {} {} {}",
                        p3,
                        p5,
                        ins.len(),
                        String::from_utf8_lossy(&ins)
                    );
                }
                let tmp = join_ref_int(reference, chr, p3, p5.saturating_sub(1));
                let p3_minus_p5 = p3.saturating_sub(p5);
                if tmp.len() > ins.len() {
                    ins_str = format!("{}^{}", p3_minus_p5, String::from_utf8_lossy(&ins));
                    bi = p3;
                } else if tmp.len() < ins.len() {
                    let part1 = substr(&ins, 0, ins.len() as isize - tmp.len() as isize);
                    let part2 = substr_from(&ins, p3_minus_p5 as isize);
                    ins_str = format!(
                        "+{}&{}",
                        String::from_utf8_lossy(part1),
                        String::from_utf8_lossy(part2)
                    );
                    bi = p3.saturating_sub(1);
                } else {
                    ins_str = format!("-{}^{}", ins.len(), String::from_utf8_lossy(&ins));
                    bi = p3;
                }
            } else {
                if seq3.len() > ins.len() {
                    let sub = substr_from(&seq3, ins.len() as isize);
                    let jref = join_ref_int(
                        reference,
                        chr,
                        p5,
                        p5.saturating_add((seq3.len() - ins.len()) as i32)
                            .saturating_add(2),
                    );
                    if !is_match(sub, &jref, 1, 3) {
                        continue;
                    }
                }
                if seq5.len() > ins.len() {
                    let sub = substr_from(&seq5, ins.len() as isize);
                    let jref = join_ref_int(
                        reference,
                        chr,
                        p3.saturating_sub((seq5.len() - ins.len()) as i32)
                            .saturating_sub(2),
                        p3.saturating_sub(1),
                    );
                    if !is_match(sub, &jref, -1, 3) {
                        continue;
                    }
                }
                let mut tmp = Vec::new();
                let p3_minus_p5 = p3.saturating_sub(p5);
                if ins.len() <= p3_minus_p5.max(0) as usize {
                    let span = p3_minus_p5.saturating_add(ins.len() as i32);
                    let mut rpt = 2;
                    let mut tnr = 3;
                    while ((span as f64 / tnr as f64) / ins.len() as f64) > 1.0 {
                        if span % tnr == 0 {
                            rpt += 1;
                        }
                        tnr += 1;
                    }
                    let limit = p5 as f64 + span as f64 / rpt as f64 - ins.len() as f64;
                    tmp.extend_from_slice(&join_ref_double(reference, chr, p5, limit));
                    ins_str = format!(
                        "+{}{}",
                        String::from_utf8_lossy(&tmp),
                        String::from_utf8_lossy(&ins)
                    );
                } else {
                    tmp.extend_from_slice(&join_ref_int(reference, chr, p5, p3.saturating_sub(1)));
                    if (ins.len() - tmp.len()) % 2 == 0 {
                        let tex = (ins.len() - tmp.len()) / 2;
                        let part1 = substr(&ins, 0, tex as isize);
                        let part2 = substr_from(&ins, tex as isize);
                        let mut check_concat = tmp.clone();
                        check_concat.extend_from_slice(part1);
                        if check_concat == part2 {
                            ins_str = format!("+{}", String::from_utf8_lossy(part2));
                        } else {
                            ins_str = format!(
                                "+{}{}",
                                String::from_utf8_lossy(&tmp),
                                String::from_utf8_lossy(&ins)
                            );
                        }
                    } else {
                        ins_str = format!(
                            "+{}{}",
                            String::from_utf8_lossy(&tmp),
                            String::from_utf8_lossy(&ins)
                        );
                    }
                }
                if config.verbose {
                    crate::verbose_eprintln!(
                        "Found lgins30: {} {} {} {} + {}",
                        p3,
                        p5,
                        ins.len(),
                        String::from_utf8_lossy(&tmp),
                        ins_str
                    );
                }
                bi = p5.saturating_sub(1);
            }

            if let Some(sc3v) = soft_clips_3_end.get_mut(&p3) {
                sc3v.used = true;
            }
            if let Some(sc5v) = soft_clips_5_end.get_mut(&p5) {
                sc5v.used = true;
            }
            let is_ins = ins_str.starts_with('+');
            let is_del = ins_str.starts_with('-');

            let vref = if is_ins {
                get_variation(insertion_variants, bi, &ins_str)
            } else {
                get_variation(non_insertion_variants, bi, &ins_str)
            };
            vref.pstd = true;
            vref.qstd = true;
            inc_cnt(ref_coverage, bi, cnt5);

            if config.verbose {
                crate::verbose_eprintln!(" lgins30 Found: '{}' {} {} {}", ins_str, bi, bp3, bp5);
            }

            let sc3v_var = soft_clips_3_end
                .get(&p3)
                .map_or(Variation::new(), |sc| sc.variation.clone());
            let sc5v_var = soft_clips_5_end
                .get(&p5)
                .map_or(Variation::new(), |sc| sc.variation.clone());

            if is_ins {
                let ref_base = reference.get_base(chr, bi as u32);
                let mvref_mut = get_variation_maybe_mut(non_insertion_variants, bi, ref_base);
                // Mutate mvref in adj_cnt for sc3v:
                let vref = get_variation(insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc3v_var, mvref_mut);

                // Mutate for sc5v:
                let vref = get_variation(insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc5v_var, None);

                let has_bams = bams.is_some_and(|b| !b.is_empty());
                let mvref_cnt =
                    get_variation_maybe(non_insertion_variants, bi, ref_base).map_or(0, |m| m.cnt);
                let bams_slice = match bams {
                    Some(b) => b,
                    None => &[],
                };
                let vref = get_variation(insertion_variants, bi, &ins_str);
                let p3_minus_p5 = p3.saturating_sub(p5);
                if has_bams
                    && p3_minus_p5 >= 5
                    && p3_minus_p5 < max_read_length.saturating_sub(10)
                    && mvref_cnt != 0
                    && no_passing_reads(chr, p5, p3, bams_slice, config)
                    && vref.cnt > mvref_cnt.saturating_mul(2)
                {
                    let mvref_mut = get_variation_maybe_mut(non_insertion_variants, bi, ref_base);
                    if let Some(m_mut) = mvref_mut {
                        let m_clone = m_mut.clone();
                        let vref = get_variation(insertion_variants, bi, &ins_str);
                        crate::variations::util::adj_cnt(vref, &m_clone, Some(m_mut));
                    }
                }

                let vref = get_variation(insertion_variants, bi, &ins_str);
                let mut tins = BTreeMap::new();
                let mut map = rustc_hash::FxHashMap::default();
                map.insert(ins_str.clone(), vref.cnt);
                tins.insert(bi, map);
                realignins(
                    bams,
                    &tins,
                    reference,
                    chr,
                    chr_len,
                    config,
                    insertion_variants,
                    non_insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                    max_read_length,
                );
            } else if is_del {
                let ref_base = reference.get_base(chr, bi as u32);

                let vref = get_variation(non_insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc3v_var, None);

                let vref = get_variation(non_insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc5v_var, None);

                if let Some(mvref) = get_variation_maybe_mut(non_insertion_variants, bi, ref_base) {
                    mvref.cnt = mvref.cnt.saturating_sub(sc3v_var.cnt);
                    mvref.hicnt = mvref.hicnt.saturating_sub(sc3v_var.hicnt);
                    mvref.locnt = mvref.locnt.saturating_sub(sc3v_var.locnt);
                    mvref.pmean -= sc3v_var.pmean;
                    mvref.qmean -= sc3v_var.qmean;
                    mvref.mqmean -= sc3v_var.mqmean;
                    mvref.nm -= sc3v_var.nm;
                    mvref.sub_dir(true, sc3v_var.get_dir(true));
                    mvref.sub_dir(false, sc3v_var.get_dir(false));
                    crate::variations::util::correct_cnt(mvref);
                }

                let vref = get_variation(non_insertion_variants, bi, &ins_str);
                let mut tdel = BTreeMap::new();
                let mut map = rustc_hash::FxHashMap::default();
                map.insert(ins_str.clone(), vref.cnt);
                tdel.insert(bi, map);
                realigndel(
                    None,
                    &tdel,
                    reference,
                    chr,
                    chr_len,
                    config,
                    non_insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                    max_read_length,
                );
            } else {
                let vref = get_variation(non_insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc3v_var, None);
                let vref = get_variation(non_insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc5v_var, None);
            }
            break;
        }
    }
    if config.verbose {
        crate::verbose_eprintln!("Done: lgins30\n");
    }
}

/// Ports `VariationRealigner.realignlgins`.
pub fn realignlgins(
    bams: Option<&[String]>,
    reference: &mut Reference,
    chr: &str,
    chr_len: u32,
    region: &crate::regions::Region,
    config: &crate::config::Config,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut rustc_hash::FxHashMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
    svfdup: &mut Vec<Sclip>,
    svrdup: &mut Vec<Sclip>,
    splice: &mut rustc_hash::FxHashSet<String>,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_sets(&config.adaptor);
    let min_reads = i32::try_from(config.min_reads).unwrap_or(i32::MAX);
    let sv_min_len = i32::try_from(config.sv_min_len).unwrap_or(i32::MAX);
    let extension = 5000;
    let region_start_i32 = i32::try_from(region.start.0).unwrap_or(i32::MAX);
    let region_end_i32 = i32::try_from(region.end.0).unwrap_or(i32::MAX);
    let min_clip_pos = region_start_i32.saturating_sub(extension);
    let max_clip_pos = region_end_i32.saturating_add(extension);

    let mut tmp = Vec::new();
    for (&p, sc5v) in soft_clips_5_end.iter() {
        if p < min_clip_pos || p > max_clip_pos {
            continue;
        }
        tmp.push((p, sc5v.variation.cnt));
    }
    tmp.sort_by(|a, b| {
        let cmp_cnt = b.1.cmp(&a.1);
        if cmp_cnt != std::cmp::Ordering::Equal {
            cmp_cnt
        } else {
            a.0.cmp(&b.0)
        }
    });

    for (p, cnt) in tmp {
        if cnt < min_reads {
            break;
        }
        if soft_clips_5_end.get(&p).is_none_or(|sc| sc.used) {
            continue;
        }
        let seq = {
            let sc5v = match soft_clips_5_end.get_mut(&p) {
                Some(sc) => sc,
                None => continue,
            };
            crate::variations::util::find_conseq_sets(sc5v, 0, &adaptor_forward, &adaptor_reverse)
        };
        if seq.is_empty() {
            continue;
        }
        if config.verbose {
            crate::verbose_eprintln!(
                "  Working lgins: 5: {} {} cnt: {}",
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }
        if seq.len() < 12 {
            continue;
        }

        let tpl = find_bi(&seq, p, reference, -1, chr, chr_len);
        let mut bi = tpl.base_insert;
        let mut ins = tpl.insertion_sequence.into_bytes();

        if bi == 0 {
            if crate::variations::util::is_low_complex_seq(&seq) {
                continue;
            }
            let match_res =
                crate::sv::find_match(&seq, reference, p, -1, 17, 1, chr, config.verbose);
            bi = match_res.base_position;
            let extra = match_res.matched_sequence;
            let bi_minus_p = i64::from(bi) - i64::from(p);
            if !(bi != 0 && bi_minus_p > 15 && bi_minus_p < 150000) {
                continue;
            }

            if bi > region_end_i32 {
                let mut tts = bi.saturating_sub(max_read_length);
                let tte = bi.saturating_add(max_read_length);
                if tts <= region_end_i32 {
                    tts = region_end_i32.saturating_add(1);
                }
                if let Err(e) = ensure_reference_loaded(
                    reference,
                    chr,
                    tts,
                    tte,
                    max_read_length,
                    chr_len,
                    config,
                ) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    region,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                    splice,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }

            const SVFLANK: i32 = 50;
            if bi.saturating_sub(p) > sv_min_len.saturating_add(2 * SVFLANK) {
                let mut tmp_ins = join_ref_int(reference, chr, p, p.saturating_add(SVFLANK - 1));
                let dup_str = format!(
                    "<dup{}>",
                    bi.saturating_sub(p).saturating_sub(2 * SVFLANK - 1)
                );
                tmp_ins.extend_from_slice(dup_str.as_bytes());
                let jref = join_ref_for_5_lgins(
                    reference,
                    chr,
                    bi.saturating_sub(SVFLANK - 1),
                    bi,
                    &seq,
                    &extra,
                );
                tmp_ins.extend_from_slice(&jref);
                ins = tmp_ins;
            } else {
                ins = join_ref_for_5_lgins(reference, chr, p, bi, &seq, &extra);
            }
            ins.extend_from_slice(&extra);

            let mut sv_lists = [&mut *svfdup, &mut *svrdup];
            let (clusters, pairs) =
                crate::sv::mark_dup_sv(p, bi, &mut sv_lists, max_read_length, config.verbose);

            let sc5v_vars_count = soft_clips_5_end.get(&p).map_or(0, |sc| sc.variation.cnt);
            let cov_bi = ref_coverage.get(&bi).copied();
            let p_minus_1 = p.saturating_sub(1);
            let cov_p_minus_1 = ref_coverage.get(&p_minus_1).copied();

            if cov_p_minus_1.is_none()
                || cov_p_minus_1.is_some_and(|cp| cov_bi.is_some_and(|cb| cp < cb))
            {
                if let Some(c_bi) = cov_bi {
                    ref_coverage.insert(p_minus_1, c_bi);
                } else {
                    ref_coverage.insert(p_minus_1, sc5v_vars_count);
                }
            } else {
                if sc5v_vars_count > cov_p_minus_1.unwrap_or(0) {
                    inc_cnt(ref_coverage, p_minus_1, sc5v_vars_count);
                }
            }

            bi = p_minus_1;
            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, bi);
            sv.r#type = "DUP".to_string();
            sv.pairs = sv.pairs.saturating_add(pairs);
            sv.splits = sv.splits.saturating_add(cnt);
            sv.clusters = sv.clusters.saturating_add(clusters);
        }

        let plus_ins = format!("+{}", String::from_utf8_lossy(&ins));
        if config.verbose {
            crate::verbose_eprintln!(
                "  Found candidate lgins from 5: {} {} {} {}",
                bi,
                plus_ins,
                p,
                String::from_utf8_lossy(&seq)
            );
        }
        let iref = get_variation(insertion_variants, bi, &plus_ins);
        iref.pstd = true;
        iref.qstd = true;

        let sc5v_var = soft_clips_5_end
            .get(&p)
            .map_or(Variation::new(), |sc| sc.variation.clone());
        let iref = get_variation(insertion_variants, bi, &plus_ins);
        crate::variations::util::adj_cnt(iref, &sc5v_var, None);

        let mut rpflag = true;
        for (i, &ch) in ins.iter().enumerate() {
            let ref_ch = i32::try_from(i)
                .ok()
                .and_then(|ii| bi.checked_add(1).and_then(|p| p.checked_add(ii)))
                .and_then(|p| u32::try_from(p).ok())
                .and_then(|pos| reference.get_base(chr, pos));
            if ref_ch != Some(ch) {
                rpflag = false;
                break;
            }
        }

        if non_insertion_variants
            .get(&bi)
            .is_some_and(|entry| entry.sv.is_none())
        {
            inc_cnt(ref_coverage, bi, sc5v_var.cnt);
        }

        let mut len = ins.len() as i32;
        if ins.contains(&b'&') {
            len -= 1;
        }
        let seq_len = soft_clips_5_end
            .get(&p)
            .map_or(0, |sc| sc.seq.keys().max().copied().unwrap_or(0))
            .saturating_add(1);

        let entries_to_process: Vec<(i32, u8, Variation)> = {
            let sc5v = soft_clips_5_end.get(&p);
            sc5v.map_or_else(Vec::new, |sc5v| {
                let mut res = Vec::new();
                for ii in len.saturating_add(1)..seq_len {
                    let Some(pii) = bi.checked_sub(ii.saturating_sub(len)) else {
                        break;
                    };
                    if let Some(map) = sc5v.seq.get(&ii) {
                        for (&tnt, tv) in map {
                            res.push((pii, tnt, tv.clone()));
                        }
                    }
                }
                res
            })
        };
        for (pii, tnt, tv) in entries_to_process {
            let mut ch_buf = [0u8; 4];
            let tnt_str: &str = (tnt as char).encode_utf8(&mut ch_buf);
            let tvr = get_variation(non_insertion_variants, pii, tnt_str);
            crate::variations::util::adj_cnt(tvr, &tv, None);
            tvr.pstd = true;
            tvr.qstd = true;
            inc_cnt(ref_coverage, pii, tv.cnt);
        }

        if let Some(sc5v) = soft_clips_5_end.get_mut(&p) {
            sc5v.used = bi.saturating_add(len) != 0;
        }

        let tins_vars_cnt = {
            let iref = get_variation(insertion_variants, bi, &plus_ins);
            iref.cnt
        };
        let mut tins = BTreeMap::new();
        let mut map = rustc_hash::FxHashMap::default();
        map.insert(plus_ins.clone(), tins_vars_cnt);
        tins.insert(bi, map);

        let newins = realignins(
            bams,
            &tins,
            reference,
            chr,
            chr_len,
            config,
            insertion_variants,
            non_insertion_variants,
            ref_coverage,
            soft_clips_5_end,
            soft_clips_3_end,
            max_read_length,
        );
        let newins = if newins.is_empty() { plus_ins } else { newins };

        if let Some(entry) = non_insertion_variants.get_mut(&bi) {
            if entry.sv.is_some() {
                let kref = get_variation(insertion_variants, bi, &newins);
                let diff = kref.cnt.saturating_sub(tins_vars_cnt);
                if let Some(sv) = entry.sv.as_mut() {
                    sv.splits = sv.splits.saturating_add(diff);
                }
            }
        }

        let ref_base = reference.get_base(chr, bi as u32);
        let mref_cloned = get_variation_maybe(non_insertion_variants, bi, ref_base).cloned();
        let has_bams = bams.is_some_and(|b| !b.is_empty());
        let bams_slice = match bams {
            Some(b) => b,
            None => &[],
        };
        if let Some(mref) = mref_cloned {
            if rpflag
                && has_bams
                && ins.len() >= 5
                && (ins.len() as i32) < max_read_length.saturating_sub(10)
                && mref.cnt != 0
                && no_passing_reads(
                    chr,
                    bi,
                    bi.saturating_add(ins.len() as i32),
                    bams_slice,
                    config,
                )
            {
                let m_cnt = mref.cnt;
                let kref = get_variation(insertion_variants, bi, &newins);
                if kref.cnt > m_cnt.saturating_mul(2) {
                    let mref_mut = get_variation_maybe_mut(non_insertion_variants, bi, ref_base);
                    if let Some(m_mut) = mref_mut {
                        let m_clone = m_mut.clone();
                        let kref = get_variation(insertion_variants, bi, &newins);
                        crate::variations::util::adj_cnt(kref, &m_clone, Some(m_mut));
                    }
                }
            }
        }
    }

    let mut tmp3 = Vec::new();
    for (&p, sc3v) in soft_clips_3_end.iter() {
        if p < min_clip_pos || p > max_clip_pos {
            continue;
        }
        tmp3.push((p, sc3v.variation.cnt));
    }
    tmp3.sort_by(|a, b| {
        let cmp_cnt = b.1.cmp(&a.1);
        if cmp_cnt != std::cmp::Ordering::Equal {
            cmp_cnt
        } else {
            a.0.cmp(&b.0)
        }
    });

    for (p, cnt) in tmp3 {
        if cnt < min_reads {
            break;
        }
        if soft_clips_3_end.get(&p).is_none_or(|sc| sc.used) {
            continue;
        }
        let seq = {
            let sc3v = match soft_clips_3_end.get_mut(&p) {
                Some(sc) => sc,
                None => continue,
            };
            crate::variations::util::find_conseq_sets(sc3v, 0, &adaptor_forward, &adaptor_reverse)
        };
        if seq.is_empty() {
            continue;
        }
        if config.verbose {
            crate::verbose_eprintln!(
                "  Working lgins 3: {} {} cnt: {}",
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }
        if seq.len() < 12 {
            continue;
        }

        let tpl = find_bi(&seq, p, reference, 1, chr, chr_len);
        let mut bi = tpl.base_insert;
        let mut ins = tpl.insertion_sequence.into_bytes();

        let orig_p = p;
        let mut p = p;

        if bi == 0 {
            if crate::variations::util::is_low_complex_seq(&seq) {
                continue;
            }
            let match_res =
                crate::sv::find_match(&seq, reference, p, 1, 17, 1, chr, config.verbose);
            bi = match_res.base_position;
            let extra = match_res.matched_sequence;
            let p_minus_bi = i64::from(p) - i64::from(bi);
            if !(bi != 0 && p_minus_bi > 15 && p_minus_bi < 150000) {
                continue;
            }

            if bi < region_start_i32 {
                let tts = bi.saturating_sub(max_read_length);
                let mut tte = bi.saturating_add(max_read_length);
                if tte >= region_start_i32 {
                    tte = region_start_i32.saturating_sub(1);
                }
                if let Err(e) = ensure_reference_loaded(
                    reference,
                    chr,
                    tts,
                    tte,
                    max_read_length,
                    chr_len,
                    config,
                ) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    region,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                    splice,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }

            let mut shift5: i32 = 0;
            loop {
                let ref_p_minus_1 = p
                    .checked_sub(1)
                    .and_then(|pos| u32::try_from(pos).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                let ref_bi_minus_1 = bi
                    .checked_sub(1)
                    .and_then(|pos| u32::try_from(pos).ok())
                    .and_then(|pos| reference.get_base(chr, pos));
                if ref_p_minus_1.is_some()
                    && ref_bi_minus_1.is_some()
                    && ref_p_minus_1 == ref_bi_minus_1
                {
                    let (Some(next_p), Some(next_bi)) = (p.checked_sub(1), bi.checked_sub(1))
                    else {
                        break;
                    };
                    p = next_p;
                    bi = next_bi;
                    shift5 = shift5.saturating_add(1);
                } else {
                    break;
                }
            }

            const SVFLANK: i32 = 50;
            if p.saturating_sub(bi) > sv_min_len.saturating_add(2 * SVFLANK) {
                ins = join_ref_for_3_lgins(
                    reference,
                    chr,
                    bi,
                    bi.saturating_add(SVFLANK - 1),
                    shift5,
                    &seq,
                    &extra,
                );
                let dup_str = format!("<dup{}>", p.saturating_sub(bi).saturating_sub(2 * SVFLANK));
                ins.extend_from_slice(dup_str.as_bytes());
                let jref = join_ref_int(
                    reference,
                    chr,
                    p.saturating_sub(SVFLANK),
                    p.saturating_sub(1),
                );
                ins.extend_from_slice(&jref);
            } else {
                ins = join_ref_for_3_lgins(
                    reference,
                    chr,
                    bi,
                    p.saturating_sub(1),
                    shift5,
                    &seq,
                    &extra,
                );
            }
            ins.extend_from_slice(&extra);

            let mut sv_lists = [&mut *svfdup, &mut *svrdup];
            let (clusters, pairs) = crate::sv::mark_dup_sv(
                bi,
                p.saturating_sub(1),
                &mut sv_lists,
                max_read_length,
                config.verbose,
            );
            bi = bi.saturating_sub(1);

            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, bi);
            sv.r#type = "DUP".to_string();
            sv.pairs = sv.pairs.saturating_add(pairs);
            sv.splits = sv.splits.saturating_add(cnt);
            sv.clusters = sv.clusters.saturating_add(clusters);

            let sc3v_vars_count = soft_clips_3_end
                .get(&orig_p)
                .map_or(0, |sc| sc.variation.cnt);
            let cov_p = ref_coverage.get(&p).copied();
            let cov_bi = ref_coverage.get(&bi).copied();
            if cov_bi.is_none() || cov_bi.is_some_and(|cb| cov_p.is_some_and(|cp| cb < cp)) {
                if let Some(c_p) = cov_p {
                    ref_coverage.insert(bi, c_p);
                } else {
                    ref_coverage.insert(bi, sc3v_vars_count);
                }
            } else {
                if sc3v_vars_count > cov_bi.unwrap_or(0) {
                    inc_cnt(ref_coverage, bi, sc3v_vars_count);
                }
            }
        }

        let plus_ins = format!("+{}", String::from_utf8_lossy(&ins));
        if config.verbose {
            crate::verbose_eprintln!(
                "  Found candidate lgins from 3: {} {} {} {}",
                bi,
                plus_ins,
                p,
                String::from_utf8_lossy(&seq)
            );
        }
        let iref = get_variation(insertion_variants, bi, &plus_ins);
        iref.pstd = true;
        iref.qstd = true;

        let ref_base = reference.get_base(chr, bi as u32);
        let mut lref_mut = get_variation_maybe_mut(non_insertion_variants, bi, ref_base);
        let sc3v_mean_pos = soft_clips_3_end
            .get(&orig_p)
            .map_or(0.0, |sc| sc.variation.pmean);
        if p as f64 - bi as f64 > sc3v_mean_pos / (cnt as f64) {
            lref_mut = None;
        }

        let sc3v_var = soft_clips_3_end
            .get(&orig_p)
            .map_or(Variation::new(), |sc| sc.variation.clone());
        let iref = get_variation(insertion_variants, bi, &plus_ins);
        if let Some(ref_mut) = lref_mut {
            crate::variations::util::adj_cnt(iref, &sc3v_var, Some(ref_mut));
        } else {
            crate::variations::util::adj_cnt(iref, &sc3v_var, None);
        }

        let mut rpflag = true;
        for (i, &ch) in ins.iter().enumerate() {
            let ref_ch = i32::try_from(i)
                .ok()
                .and_then(|ii| bi.checked_add(1).and_then(|p| p.checked_add(ii)))
                .and_then(|p| u32::try_from(p).ok())
                .and_then(|pos| reference.get_base(chr, pos));
            if ref_ch != Some(ch) {
                rpflag = false;
                break;
            }
        }

        let mut len = ins.len() as i32;
        if ins.contains(&b'&') {
            len -= 1;
        }
        let len_seq = soft_clips_3_end
            .get(&orig_p)
            .map_or(0, |sc| sc.seq.keys().max().copied().unwrap_or(0))
            .saturating_add(1);

        let entries_to_process: Vec<(i32, u8, Variation)> = {
            let sc3v = soft_clips_3_end.get(&orig_p);
            sc3v.map_or_else(Vec::new, |sc3v| {
                let mut res = Vec::new();
                for ii in len..len_seq {
                    let Some(pii) = p.checked_add(ii.saturating_sub(len)) else {
                        break;
                    };
                    if let Some(map) = sc3v.seq.get(&ii) {
                        for (&tnt, tv) in map {
                            res.push((pii, tnt, tv.clone()));
                        }
                    }
                }
                res
            })
        };
        for (pii, tnt, tv) in entries_to_process {
            let mut ch_buf = [0u8; 4];
            let tnt_str: &str = (tnt as char).encode_utf8(&mut ch_buf);
            let tvr = get_variation(non_insertion_variants, pii, tnt_str);
            crate::variations::util::adj_cnt(tvr, &tv, None);
            tvr.pstd = true;
            tvr.qstd = true;
            inc_cnt(ref_coverage, pii, tv.cnt);
        }

        if let Some(sc3v) = soft_clips_3_end.get_mut(&orig_p) {
            sc3v.used = true;
        }

        let tins_vars_cnt = {
            let iref = get_variation(insertion_variants, bi, &plus_ins);
            iref.cnt
        };
        let mut tins = BTreeMap::new();
        let mut map = rustc_hash::FxHashMap::default();
        map.insert(plus_ins.clone(), tins_vars_cnt);
        tins.insert(bi, map);

        let newins = realignins(
            bams,
            &tins,
            reference,
            chr,
            chr_len,
            config,
            insertion_variants,
            non_insertion_variants,
            ref_coverage,
            soft_clips_5_end,
            soft_clips_3_end,
            max_read_length,
        );
        let newins = if newins.is_empty() { plus_ins } else { newins };

        if let Some(entry) = non_insertion_variants.get_mut(&bi) {
            if entry.sv.is_some() {
                let iref = get_variation(insertion_variants, bi, &newins);
                let diff = iref.cnt.saturating_sub(tins_vars_cnt);
                if let Some(sv) = entry.sv.as_mut() {
                    sv.splits = sv.splits.saturating_add(diff);
                }
            }
        }

        let mref_cloned = get_variation_maybe(non_insertion_variants, bi, ref_base).cloned();
        let has_bams = bams.is_some_and(|b| !b.is_empty());
        let bams_slice = match bams {
            Some(b) => b,
            None => &[],
        };
        if let Some(mref) = mref_cloned {
            if rpflag
                && has_bams
                && ins.len() >= 5
                && (ins.len() as i32) < max_read_length.saturating_sub(10)
                && mref.cnt != 0
                && no_passing_reads(
                    chr,
                    bi,
                    bi.saturating_add(ins.len() as i32),
                    bams_slice,
                    config,
                )
            {
                let iref = get_variation(insertion_variants, bi, &newins);
                let m_cnt = mref.cnt;
                if iref.cnt > m_cnt.saturating_mul(2) {
                    let mref_mut = get_variation_maybe_mut(non_insertion_variants, bi, ref_base);
                    if let Some(m_mut) = mref_mut {
                        let m_clone = m_mut.clone();
                        let iref = get_variation(insertion_variants, bi, &newins);
                        crate::variations::util::adj_cnt(iref, &m_clone, Some(m_mut));
                    }
                }
            }
        }
    }
}

/// Main orchestrator of local realignment.
/// Ports `VariationRealigner.realignIndels`.
pub fn realign_indels(
    bams: Option<&[String]>,
    reference: &mut Reference,
    chr: &str,
    chr_len: u32,
    region: &crate::regions::Region,
    config: &crate::config::Config,
    var_data: &mut crate::variations::VariationData,
) {
    if !config.no_sv {
        crate::sv::filter_all_sv_structures(var_data);
    }

    if config.verbose {
        crate::verbose_eprintln!("Start local realignment...");
    }

    // 1. realigndel
    realigndel(
        bams,
        &var_data.position_to_deletion_count,
        reference,
        chr,
        chr_len,
        config,
        &mut var_data.non_insertion_variants,
        &mut var_data.ref_coverage,
        &mut var_data.soft_clips_5_end,
        &mut var_data.soft_clips_3_end,
        var_data.max_read_length,
    );

    // 2. realignins
    realignins(
        bams,
        &var_data.position_to_insertion_count,
        reference,
        chr,
        chr_len,
        config,
        &mut var_data.insertion_variants,
        &mut var_data.non_insertion_variants,
        &mut var_data.ref_coverage,
        &mut var_data.soft_clips_5_end,
        &mut var_data.soft_clips_3_end,
        var_data.max_read_length,
    );

    // 3. realignlgdel
    realignlgdel(
        bams,
        reference,
        chr,
        chr_len,
        region,
        config,
        &mut var_data.non_insertion_variants,
        &mut var_data.insertion_variants,
        &mut var_data.ref_coverage,
        &mut var_data.soft_clips_5_end,
        &mut var_data.soft_clips_3_end,
        var_data.max_read_length,
        &mut var_data.sv_structures.svfdel,
        &mut var_data.sv_structures.svrdel,
        &mut var_data.splice,
    );

    // 4. realignlgins30
    realignlgins30(
        bams,
        reference,
        chr,
        chr_len,
        region,
        config,
        &mut var_data.non_insertion_variants,
        &mut var_data.insertion_variants,
        &mut var_data.ref_coverage,
        &mut var_data.soft_clips_5_end,
        &mut var_data.soft_clips_3_end,
        var_data.max_read_length,
    );

    // 5. realignlgins
    realignlgins(
        bams,
        reference,
        chr,
        chr_len,
        region,
        config,
        &mut var_data.non_insertion_variants,
        &mut var_data.insertion_variants,
        &mut var_data.ref_coverage,
        &mut var_data.soft_clips_5_end,
        &mut var_data.soft_clips_3_end,
        var_data.max_read_length,
        &mut var_data.sv_structures.svfdup,
        &mut var_data.sv_structures.svrdup,
        &mut var_data.splice,
    );

    if config.verbose {
        crate::verbose_eprintln!("Finished local realignment.");
    }
}

fn substr(seq: &[u8], mut begin: isize, len: isize) -> &[u8] {
    let s_len = isize::try_from(seq.len()).unwrap_or(isize::MAX);
    if begin < 0 {
        begin = s_len.saturating_add(begin);
    }
    if begin < 0 {
        begin = 0;
    }
    if begin >= s_len || len == 0 {
        return b"";
    }
    let Ok(start) = usize::try_from(begin) else {
        return b"";
    };

    if len > 0 {
        let end_isize = begin.saturating_add(len).min(s_len);
        let Ok(end) = usize::try_from(end_isize) else {
            return b"";
        };
        seq.get(start..end).unwrap_or(b"")
    } else {
        let end_isize = s_len.saturating_add(len);
        if end_isize <= begin {
            b""
        } else {
            let Ok(end) = usize::try_from(end_isize) else {
                return b"";
            };
            seq.get(start..end).unwrap_or(b"")
        }
    }
}

fn substr_from(seq: &[u8], idx: isize) -> &[u8] {
    let s_len = isize::try_from(seq.len()).unwrap_or(isize::MAX);
    let start_isize = if idx >= 0 {
        idx.min(s_len)
    } else {
        s_len.saturating_add(idx).max(0)
    };
    let Ok(start) = usize::try_from(start_isize) else {
        return b"";
    };
    seq.get(start..).unwrap_or(b"")
}

fn reverse(seq: &[u8]) -> Vec<u8> {
    let mut v = seq.to_vec();
    v.reverse();
    v
}

fn join_ref_int(reference: &Reference, chr: &str, from: i32, to: i32) -> Vec<u8> {
    crate::util::ref_bases(reference, chr, from, to)
}

fn join_ref_double(reference: &Reference, chr: &str, from: i32, to: f64) -> Vec<u8> {
    // The float bound iterates while (i as f64) < to; that is the inclusive
    // integer range from..=ceil(to)-1 (empty when NaN/negative).
    crate::util::ref_bases(reference, chr, from, (to.ceil() as i32).saturating_sub(1))
}

/// Ports `VariationUtils.java:207-224` (`joinRefFor5Lgins`).
fn join_ref_for_5_lgins(
    reference: &Reference,
    chr: &str,
    from: i32,
    to: i32,
    seq: &[u8],
    extra: &[u8],
) -> Vec<u8> {
    let mut sb = Vec::new();
    let keep_len = seq.len().saturating_sub(extra.len());
    if from <= to {
        for pos64 in i64::from(from)..=i64::from(to) {
            let diff64 = i64::from(to) - pos64;
            if diff64 >= 0 && (diff64 as usize) < keep_len {
                if let Some(&b) = (diff64 as usize)
                    .checked_add(extra.len())
                    .and_then(|idx| seq.get(idx))
                {
                    sb.push(b);
                }
            } else if pos64 >= 1 && pos64 <= i64::from(i32::MAX) {
                if let Some(ch) = reference.get_base(chr, pos64 as u32) {
                    sb.push(ch);
                }
            }
        }
    }
    sb
}

/// Ports `VariationUtils.java:236-255` (`joinRefFor3Lgins`).
fn join_ref_for_3_lgins(
    reference: &Reference,
    chr: &str,
    from: i32,
    to: i32,
    shift5: i32,
    seq: &[u8],
    extra: &[u8],
) -> Vec<u8> {
    let mut sb = Vec::new();
    let keep_len = seq.len().saturating_sub(extra.len());
    if from <= to {
        for pos64 in i64::from(from)..=i64::from(to) {
            let diff = pos64 - i64::from(from);
            let idx64 = diff - i64::from(shift5);
            if idx64 >= 0 && (idx64 as usize) < keep_len {
                if let Some(&b) = (idx64 as usize)
                    .checked_add(extra.len())
                    .and_then(|idx| seq.get(idx))
                {
                    sb.push(b);
                }
            } else if pos64 >= 1 && pos64 <= i64::from(i32::MAX) {
                if let Some(ch) = reference.get_base(chr, pos64 as u32) {
                    sb.push(ch);
                }
            }
        }
    }
    sb
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Match35 {
    pub matched_5_end: i32,
    pub matched_3_end: i32,
    pub max_matched_length: i32,
}

#[inline]
fn swar_mismatch_count_u64(a: u64, b: u64) -> i32 {
    let diff = a ^ b;
    if diff == 0 {
        return 0;
    }
    let nz = (diff | ((diff & 0x7F7F_7F7F_7F7F_7F7F).wrapping_add(0x7F7F_7F7F_7F7F_7F7F)))
        & 0x8080_8080_8080_8080;
    nz.count_ones() as i32
}

#[inline]
fn load_u64_le(slice: &[u8], offset: usize) -> u64 {
    let mut buf = [0u8; 8];
    if let Some(src) = slice.get(offset..offset + 8) {
        buf.copy_from_slice(src);
    }
    u64::from_le_bytes(buf)
}

pub fn find_35_match(seq5: &[u8], seq3: &[u8], verbose: bool) -> Match35 {
    let long_mismatch = 2;
    let mut max_matched_length = 0;
    let mut b3 = 0;
    let mut b5 = 0;

    let len5 = seq5.len();
    let len3 = seq3.len();

    if len5 <= 8 || len3 <= 9 {
        return Match35 {
            matched_5_end: 0,
            matched_3_end: 0,
            max_matched_length: 0,
        };
    }

    /// Maximum stack elements (`512` bytes for `rev3_stack`, `4 KiB` for `rev3_w0_stack`,
    /// covering `len3 <= 512 bp`) before falling back to heap allocation on 2 MiB rayon stacks.
    const MAX_STACK_REV3: usize = 512;

    let mut rev3_stack = [0u8; MAX_STACK_REV3];
    let mut rev3_heap = Vec::new();
    let rev3: &[u8] = if let Some(dst) = rev3_stack.get_mut(..len3) {
        for (slot, &b) in dst.iter_mut().zip(seq3.iter().rev()) {
            *slot = b;
        }
        dst
    } else {
        rev3_heap.extend(seq3.iter().rev().copied());
        &rev3_heap
    };

    // Precompute the first 8-byte word of `rev3` at each starting offset `j - 1`
    // for `j in 1..(len3 - 8)` (`num_j = len3 - 9` words; since `j <= len3 - 9`,
    // `(j - 1) + 8 <= len3 - 2 < len3` is always in bounds).
    let num_j = len3 - 9;
    let mut rev3_w0_stack = [0u64; MAX_STACK_REV3];
    let mut rev3_w0_heap = Vec::new();
    let rev3_w0: &[u64] = if let Some(dst) = rev3_w0_stack.get_mut(..num_j) {
        for (idx, slot) in dst.iter_mut().enumerate() {
            *slot = load_u64_le(rev3, idx);
        }
        dst
    } else {
        rev3_w0_heap.reserve(num_j);
        for idx in 0..num_j {
            rev3_w0_heap.push(load_u64_le(rev3, idx));
        }
        &rev3_w0_heap
    };

    for i in 0..(len5 - 8) {
        let rem5 = len5 - i;
        let w5_0 = load_u64_le(seq5, i);
        for (j_idx, &w3_0) in rev3_w0.iter().enumerate() {
            // Fast SWAR reject: since `i < len5 - 8` (`rem5 >= 9`) and
            // `j < len3 - 8` (`rem3 >= 10`), `valid_len = rem3.min(rem5) >= 9`.
            // If the first 8 bytes alone have > 2 mismatches, the scalar loop is
            // guaranteed to break at `total_length <= 7` with `number_of_mismatch == 3`,
            // yielding `score <= 4 < 8`.
            let mm0 = swar_mismatch_count_u64(w3_0, w5_0);
            if mm0 > long_mismatch {
                continue;
            }

            let j = j_idx + 1;
            let rem3 = len3 - j + 1;
            let valid_len = rem3.min(rem5);
            let s3 = rev3.get(j_idx..(j_idx + valid_len)).unwrap_or(&[]);
            let s5 = seq5.get(i..(i + valid_len)).unwrap_or(&[]);
            let mut number_of_mismatch = mm0;
            let mut total_length = 8usize;

            while total_length + 8 <= valid_len {
                let w3 = load_u64_le(s3, total_length);
                let w5 = load_u64_le(s5, total_length);
                let chunk_mm = swar_mismatch_count_u64(w3, w5);
                if number_of_mismatch + chunk_mm <= long_mismatch {
                    number_of_mismatch += chunk_mm;
                    total_length += 8;
                } else {
                    break;
                }
            }

            if let (Some(rem_s3), Some(rem_s5)) = (s3.get(total_length..), s5.get(total_length..)) {
                for (&b3_ch, &b5_ch) in rem_s3.iter().zip(rem_s5.iter()) {
                    if b3_ch != b5_ch {
                        number_of_mismatch += 1;
                        if number_of_mismatch > long_mismatch {
                            break;
                        }
                    }
                    total_length += 1;
                }
            }
            if number_of_mismatch <= long_mismatch && rem5 < rem3 {
                // Legacy while-loop sentinel iteration at i + total_length == len5
                number_of_mismatch += 1;
                if number_of_mismatch <= long_mismatch {
                    total_length += 1;
                }
            }
            let score = total_length as i32 - number_of_mismatch;
            if score > max_matched_length
                && score > 8
                && (number_of_mismatch as f64) / (total_length as f64) < 0.1
                && (total_length + j >= len3 || i + total_length >= len5)
            {
                max_matched_length = score;
                b3 = j as i32;
                b5 = i as i32;
                if verbose {
                    crate::verbose_eprintln!(
                        "      Found 35 Match, seq5_len={}, seq3_len={}, total_length={}, mismatch={}, i={}, j={}",
                        len5, len3, total_length, number_of_mismatch, i, j
                    );
                }
                return Match35 {
                    matched_5_end: b5,
                    matched_3_end: b3,
                    max_matched_length,
                };
            }
        }
    }
    Match35 {
        matched_5_end: b5,
        matched_3_end: b3,
        max_matched_length,
    }
}

#[cfg(test)]
mod task7_tests {
    #[test]
    fn find_35_match_matches_legacy_on_random_and_adversarial_corpus() {
        crate::variations::util::tests::run_find_35_match_diff_suite();
    }
}

fn inc_cnt(hash: &mut rustc_hash::FxHashMap<i32, i32>, start: i32, count: i32) {
    let entry = hash.entry(start).or_insert(0);
    *entry = entry.saturating_add(count);
}

fn rm_cnt(vref: &mut crate::variations::Variation, tv: &crate::variations::Variation) {
    vref.cnt = vref.cnt.saturating_sub(tv.cnt);
    vref.hicnt = vref.hicnt.saturating_sub(tv.hicnt);
    vref.locnt = vref.locnt.saturating_sub(tv.locnt);
    vref.pmean -= tv.pmean;
    vref.qmean -= tv.qmean;
    vref.mqmean -= tv.mqmean;
    // ports VariationRealigner.java:2732 (rmCnt does not subtract nm)
    vref.sub_dir(true, tv.get_dir(true));
    vref.sub_dir(false, tv.get_dir(false));
    crate::variations::util::correct_cnt(vref);
}

fn add_var_factor(vref: &mut crate::variations::Variation, factor_f: f64) {
    if factor_f < -1.0 {
        return;
    }
    vref.cnt = vref.cnt.saturating_add((factor_f * vref.cnt as f64) as i32);
    vref.hicnt = vref
        .hicnt
        .saturating_add((factor_f * vref.hicnt as f64) as i32);
    vref.locnt = vref
        .locnt
        .saturating_add((factor_f * vref.locnt as f64) as i32);
    vref.pmean += factor_f * vref.pmean;
    vref.qmean += factor_f * vref.qmean;
    vref.mqmean += factor_f * vref.mqmean;
    vref.nm += factor_f * vref.nm;
    vref.fwd = vref.fwd.saturating_add((factor_f * vref.fwd as f64) as i32);
    vref.rev = vref.rev.saturating_add((factor_f * vref.rev as f64) as i32);
}

fn ensure_reference_loaded(
    reference: &mut Reference,
    chr: &str,
    start: i32,
    end: i32,
    max_read_length: i32,
    chr_len: u32,
    config: &crate::config::Config,
) -> Result<(), crate::types::VardictError> {
    if crate::reference::is_loaded(chr, start.max(1) as u32, end.max(1) as u32, reference) {
        return Ok(());
    }
    let mod_region = crate::regions::Region {
        chr: chr.to_string(),
        start: Pos1(start.max(1) as u32),
        end: Pos1(end.max(1) as u32),
        gene: String::new(),
        amplicon: None,
    };
    crate::reference::get_reference_with_chr_len(
        &mod_region,
        max_read_length.max(0) as u32,
        chr_len,
        config,
        reference,
    )
}

fn merge_sclip(dest: &mut Sclip, src: Sclip) {
    let d = &mut dest.variation;
    d.cnt = d.cnt.saturating_add(src.variation.cnt);
    d.fwd = d.fwd.saturating_add(src.variation.fwd);
    d.rev = d.rev.saturating_add(src.variation.rev);
    d.pmean += src.variation.pmean;
    d.qmean += src.variation.qmean;
    d.mqmean += src.variation.mqmean;
    d.nm += src.variation.nm;
    d.hicnt = d.hicnt.saturating_add(src.variation.hicnt);
    d.locnt = d.locnt.saturating_add(src.variation.locnt);
    d.extracnt = d.extracnt.saturating_add(src.variation.extracnt);
    if src.variation.pstd {
        d.pstd = true;
    }
    if src.variation.qstd {
        d.qstd = true;
    }
    // pp/pq are last-seen trackers (not sums); no consumer reads them on a
    // merged sclip, so accumulating them here is inert.
    d.pp = d.pp.saturating_add(src.variation.pp);
    d.pq += src.variation.pq;

    // SV fields: Java's partialPipeline re-parses into the SAME live Sclip
    // entries via CigarParser.addSV/addadddisccnt (CigarParser.java:2295-2341),
    // so replicate its min/max accumulation for start/end/mstart/mend. A zero
    // source field means that entry saw no SV reads and must not clobber dest.
    if src.start != 0 && (dest.start == 0 || dest.start >= src.start) {
        dest.start = src.start;
    }
    if src.end != 0 && (dest.end == 0 || dest.end <= src.end) {
        dest.end = src.end;
    }
    if src.mstart != 0 && (dest.mstart == 0 || dest.mstart >= src.mstart) {
        dest.mstart = src.mstart;
    }
    if src.mend != 0 && (dest.mend == 0 || dest.mend <= src.mend) {
        // Java stores mateStart_ms + rlen of the triggering read; the closest
        // accumulated equivalent is the source entry's own mend.
        dest.mend = src.mend;
    }
    dest.disc = dest.disc.saturating_add(src.disc);
    for (pos, cnt) in src.soft {
        let entry = dest.soft.entry(pos).or_insert(0);
        *entry = entry.saturating_add(cnt);
    }
    dest.mates.extend(src.mates);

    for (idx, cnts) in src.nt {
        let dest_cnts = dest.nt.entry(idx).or_default();
        for (base, c) in cnts {
            let entry = dest_cnts.entry(base).or_insert(0);
            *entry = entry.saturating_add(c);
        }
    }
    for (idx, bases) in src.seq {
        let dest_bases = dest.seq.entry(idx).or_default();
        for (base, var) in bases {
            let dv = dest_bases.entry(base).or_insert_with(Variation::new);
            dv.cnt = dv.cnt.saturating_add(var.cnt);
            dv.fwd = dv.fwd.saturating_add(var.fwd);
            dv.rev = dv.rev.saturating_add(var.rev);
            dv.pmean += var.pmean;
            dv.qmean += var.qmean;
            dv.mqmean += var.mqmean;
            dv.nm += var.nm;
            dv.hicnt = dv.hicnt.saturating_add(var.hicnt);
            dv.locnt = dv.locnt.saturating_add(var.locnt);
            dv.extracnt = dv.extracnt.saturating_add(var.extracnt);
            if var.pstd {
                dv.pstd = true;
            }
            if var.qstd {
                dv.qstd = true;
            }
            dv.pp = dv.pp.saturating_add(var.pp);
            dv.pq += var.pq;
        }
    }
    // VariationUtils.findconseq caches the consensus on the Sclip
    // (VariationUtils.java:69-71). The merged counts make both sides' caches
    // stale, so clear ours to force recomputation from the merged counts.
    dest.sequence = None;
}

fn partial_pipeline(
    bams: Option<&[String]>,
    chr: &str,
    tts: i32,
    tte: i32,
    reference: &Reference,
    config: &crate::config::Config,
    region: &crate::regions::Region,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut rustc_hash::FxHashMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    splice: &mut rustc_hash::FxHashSet<String>,
) -> Result<(), crate::types::VardictError> {
    let bams = match bams {
        Some(b) => b,
        None => return Ok(()),
    };
    if bams.is_empty() {
        return Ok(());
    }

    let mut preprocessor = crate::bam::RecordPreprocessor::new();
    let pileup_region = crate::regions::Region {
        chr: chr.to_string(),
        start: Pos1(tts.max(1) as u32),
        end: Pos1(tte.max(1) as u32),
        gene: String::new(),
        amplicon: None,
    };
    // Streaming counterpart of the old query_region-Vec + pileup shape: one
    // lazy pass per BAM, BAM1-then-BAM2, through a single preprocessor into a
    // single parser. Peak memory stays flat in candidate-window depth.
    let mut parser =
        crate::cigar::parse::CigarParser::new(reference, &pileup_region, config, true, 0, 0, 0);
    for (idx, bam) in bams.iter().enumerate() {
        preprocessor.reset_for_next_bam();
        let bam_path = std::path::Path::new(bam);
        let mut feed = crate::bam::query_region_iter(
            bam_path,
            chr,
            Pos0(tts.max(0) as u32),
            Pos0(tte.max(0) as u32),
            config,
        )?;
        let mut r = crate::bam::SamRecord::new();
        let mut any_passed = false;
        while feed.next_into(&mut r)? {
            if preprocessor.preprocess_record(&r, config) {
                any_passed = true;
                parser.parse_cigar(&r)?;
            }
        }
        if idx > 0 && !any_passed {
            break;
        }
    }

    let var_data = parser.into_variation_data();
    splice.extend(var_data.splice);

    // G3 (mate-window ownership): attribute — don't drop — merged
    // mate-window VARIANT keys outside the owning tile's active range
    // (`region` is the tile's active range in tiled runs, the whole BED
    // region otherwise). The builder ownership gate (post-find_all_sv, with
    // SV data + both ranges + edge clamping) adjudicates emission, exactly
    // like the sv-path `merge_planted` set. Dropping here zeroes mate
    // ref-support for out-of-BED SV rows (corpus CandDel5SV, MarkingDUPSV,
    // DelOutOfBound): at merge time temp maps never carry SV data (the
    // parser only uses svflag for duprate), so no merge-time predicate can
    // spare them. Only variant maps are attributed: ref_coverage/soft-clips
    // are trigger state, not emission (the builder iterates variant maps),
    // and gating them would perturb downstream sv-trigger firing vs Java.
    let region_start_i32 = i32::try_from(region.start.0).unwrap_or(i32::MAX);
    let region_end_i32 = i32::try_from(region.end.0).unwrap_or(i32::MAX);
    let mate_owned = |pos: i32| pos >= region_start_i32 && pos <= region_end_i32;

    for (pos, cov) in var_data.ref_coverage {
        let entry = ref_coverage.entry(pos).or_insert(0);
        *entry = entry.saturating_add(cov);
    }
    // ports VariationRealigner.java:1663-1671 / CigarParser.java:1922-1927: Java's
    // partialPipeline mutates the SAME live Sclip maps, so re-parsed reads
    // ACCUMULATE into existing entries. Merge by summing all counts; never
    // discard fresh counts when the position already exists.
    for (pos, sc) in var_data.soft_clips_5_end {
        merge_sclip(soft_clips_5_end.entry(pos).or_insert_with(Sclip::new), sc);
    }
    for (pos, sc) in var_data.soft_clips_3_end {
        merge_sclip(soft_clips_3_end.entry(pos).or_insert_with(Sclip::new), sc);
    }
    for (pos, var_map) in var_data.non_insertion_variants {
        let dest = non_insertion_variants.get_or_put_with(pos, VariationMap::new);
        if !mate_owned(pos) {
            dest.merge_planted = true;
        }
        for (allele, var) in var_map.map {
            if let Some(v) = dest.map.get_mut(&allele) {
                v.cnt = v.cnt.saturating_add(var.cnt);
                v.fwd = v.fwd.saturating_add(var.fwd);
                v.rev = v.rev.saturating_add(var.rev);
                v.pmean += var.pmean;
                v.qmean += var.qmean;
                v.mqmean += var.mqmean;
                v.nm += var.nm;
                v.hicnt = v.hicnt.saturating_add(var.hicnt);
                v.locnt = v.locnt.saturating_add(var.locnt);
                v.extracnt = v.extracnt.saturating_add(var.extracnt);
                if var.pstd {
                    v.pstd = true;
                }
                if var.qstd {
                    v.qstd = true;
                }
                v.pp = v.pp.saturating_add(var.pp);
                v.pq += var.pq;
            } else {
                dest.map.insert(allele, var);
            }
        }
        if let Some(sv) = var_map.sv {
            if let Some(dest_sv) = &mut dest.sv {
                dest_sv.pairs = dest_sv.pairs.saturating_add(sv.pairs);
                dest_sv.splits = dest_sv.splits.saturating_add(sv.splits);
                dest_sv.clusters = dest_sv.clusters.saturating_add(sv.clusters);
            } else {
                dest.sv = Some(sv);
            }
        }
    }
    for (pos, var_map) in var_data.insertion_variants {
        let dest = insertion_variants.get_or_put_with(pos, VariationMap::new);
        if !mate_owned(pos) {
            dest.merge_planted = true;
        }
        for (allele, var) in var_map.map {
            if let Some(v) = dest.map.get_mut(&allele) {
                v.cnt = v.cnt.saturating_add(var.cnt);
                v.fwd = v.fwd.saturating_add(var.fwd);
                v.rev = v.rev.saturating_add(var.rev);
                v.pmean += var.pmean;
                v.qmean += var.qmean;
                v.mqmean += var.mqmean;
                v.nm += var.nm;
                v.hicnt = v.hicnt.saturating_add(var.hicnt);
                v.locnt = v.locnt.saturating_add(var.locnt);
                v.extracnt = v.extracnt.saturating_add(var.extracnt);
                if var.pstd {
                    v.pstd = true;
                }
                if var.qstd {
                    v.qstd = true;
                }
                v.pp = v.pp.saturating_add(var.pp);
                v.pq += var.pq;
            } else {
                dest.map.insert(allele, var);
            }
        }
        if let Some(sv) = var_map.sv {
            if let Some(dest_sv) = &mut dest.sv {
                dest_sv.pairs = dest_sv.pairs.saturating_add(sv.pairs);
                dest_sv.splits = dest_sv.splits.saturating_add(sv.splits);
                dest_sv.clusters = dest_sv.clusters.saturating_add(sv.clusters);
            } else {
                dest.sv = Some(sv);
            }
        }
    }

    Ok(())
}

/// Adjust MNP when there are breakpoints within MNP.
/// Ports `VariationRealigner.java:277-382`.
fn format_mnp_key_into(out: &mut String, raw_buf: &mut Vec<u8>, bytes: &[u8]) {
    out.clear();
    if let (Some(&first), Some(rest)) = (bytes.first(), bytes.get(1..)) {
        if !rest.is_empty() {
            if bytes.is_ascii() {
                out.reserve(bytes.len() + 1);
                out.push(first as char);
                out.push('&');
                if let Ok(s) = std::str::from_utf8(rest) {
                    out.push_str(s);
                }
                return;
            }
            raw_buf.clear();
            raw_buf.reserve(bytes.len() + 1);
            raw_buf.push(first);
            raw_buf.push(b'&');
            raw_buf.extend_from_slice(rest);
            out.push_str(&String::from_utf8_lossy(raw_buf));
            return;
        }
    }
    if bytes.is_ascii() {
        if let Ok(s) = std::str::from_utf8(bytes) {
            out.push_str(s);
        }
    } else {
        out.push_str(&String::from_utf8_lossy(bytes));
    }
}

pub fn adjust_mnp(
    reference: &Reference,
    chr: &str,
    _chr_len: u32,
    _region: &crate::regions::Region,
    config: &crate::config::Config,
    var_data: &mut crate::variations::VariationData,
) {
    let tmp = fill_and_sort_tmp(&var_data.mnp);
    let adaptor_forward = AdaptorSet::default();
    let adaptor_reverse = AdaptorSet::default();
    let mut left = String::new();
    let mut right = String::new();
    let mut raw_buf = Vec::new();

    for tpl in tmp {
        let position = tpl.position;
        let vn = tpl.description_string;

        let mut vref = {
            let vars_on_pos = match var_data.non_insertion_variants.get(&position) {
                Some(m) => m,
                None => continue,
            };
            match vars_on_pos.map.get(&vn) {
                Some(v) => v.clone(),
                None => continue,
            }
        };

        if config.verbose {
            crate::verbose_eprintln!("  AdjMnt: {} {} {}", position, vn, vref.cnt);
        }

        let mnt = vn.replacen('&', "", 1);
        let mnt_bytes = mnt.as_bytes();
        let mnt_len = mnt_bytes.len();
        for i in 0..mnt_len.saturating_sub(1) {
            let left_len = i + 1;
            let left_bytes = substr(
                mnt_bytes,
                0,
                isize::try_from(left_len).unwrap_or(isize::MAX),
            );
            format_mnp_key_into(&mut left, &mut raw_buf, left_bytes);

            let right_len = mnt_len - i - 1;
            let right_len_isize = isize::try_from(right_len).unwrap_or(isize::MAX);
            let right_bytes = substr(mnt_bytes, right_len_isize.saturating_neg(), right_len_isize);
            format_mnp_key_into(&mut right, &mut raw_buf, right_bytes);

            // --- LEFT CHECK ---
            let mut left_tref = None;
            {
                if let Some(vars_on_pos) = var_data.non_insertion_variants.get(&position) {
                    if let Some(tref) = vars_on_pos.map.get(&left) {
                        if tref.cnt > 0
                            && tref.cnt < vref.cnt
                            && (tref.pmean / tref.cnt as f64) <= (i + 1) as f64
                        {
                            left_tref = Some(tref.clone());
                        }
                    }
                }
            }
            if let Some(tref_clone) = left_tref {
                if config.verbose {
                    crate::verbose_eprintln!(
                        "    AdjMnt Left: {} {} Left: {} Cnt: {}",
                        position,
                        vn,
                        left,
                        tref_clone.cnt
                    );
                }
                if let Some(vars_on_pos) = var_data.non_insertion_variants.get_mut(&position) {
                    if let Some(vref_mut) = vars_on_pos.map.get_mut(&vn) {
                        crate::variations::util::adj_cnt(vref_mut, &tref_clone, None);
                        vref = vref_mut.clone();
                    }
                    vars_on_pos.map.shift_remove(&left);
                }
            }

            // --- RIGHT CHECK ---
            let right_pos = position
                .saturating_add(i32::try_from(i).unwrap_or(i32::MAX))
                .saturating_add(1);
            let mut right_tref = None;
            {
                if let Some(right_vars) = var_data.non_insertion_variants.get(&right_pos) {
                    if let Some(tref) = right_vars.map.get(&right) {
                        if tref.cnt >= 0 && tref.cnt < vref.cnt {
                            right_tref = Some(tref.clone());
                        }
                    }
                }
            }
            if let Some(tref_clone) = right_tref {
                if config.verbose {
                    crate::verbose_eprintln!(
                        "    AdjMnt Right: {} {} Right: {} Cnt: {}",
                        position,
                        vn,
                        right,
                        tref_clone.cnt
                    );
                }
                if let Some(vars_on_pos) = var_data.non_insertion_variants.get_mut(&position) {
                    if let Some(vref_mut) = vars_on_pos.map.get_mut(&vn) {
                        crate::variations::util::adj_cnt(vref_mut, &tref_clone, None);
                        vref = vref_mut.clone();
                    }
                }
                inc_cnt(&mut var_data.ref_coverage, position, tref_clone.cnt);
                if let Some(right_vars) = var_data.non_insertion_variants.get_mut(&right_pos) {
                    right_vars.map.shift_remove(&right);
                }
            }
        }

        let mnt_len_i32 = i32::try_from(mnt_len).unwrap_or(i32::MAX);

        // --- SOFTCLIPS3END CHECK ---
        let mut sc3_vars_to_add = None;
        {
            if let Some(sc3v) = var_data.soft_clips_3_end.get_mut(&position) {
                if !sc3v.used {
                    let seq = crate::variations::util::find_conseq_sets(
                        sc3v,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if seq.starts_with(mnt_bytes) {
                        let tail_matches = match seq.get(mnt_len..) {
                            Some([]) => true,
                            Some(tail) => is_match_ref(
                                tail,
                                reference,
                                chr,
                                position.saturating_add(mnt_len_i32),
                                1,
                                3,
                            ),
                            None => false,
                        };
                        if tail_matches {
                            sc3_vars_to_add = Some(sc3v.variation.clone());
                        }
                    }
                }
            }
        }
        if let Some(sc3_vars) = sc3_vars_to_add {
            if let Some(vars_on_pos) = var_data.non_insertion_variants.get_mut(&position) {
                if let Some(vref_mut) = vars_on_pos.map.get_mut(&vn) {
                    crate::variations::util::adj_cnt(vref_mut, &sc3_vars, None);
                    vref = vref_mut.clone();
                }
            }
            inc_cnt(&mut var_data.ref_coverage, position, sc3_vars.cnt);
            if let Some(sc3v) = var_data.soft_clips_3_end.get_mut(&position) {
                sc3v.used = true;
            }
        }

        // --- SOFTCLIPS5END CHECK ---
        let right_pos = position.saturating_add(mnt_len_i32);
        let mut sc5_vars_to_add = None;
        {
            if let Some(sc5v) = var_data.soft_clips_5_end.get_mut(&right_pos) {
                if !sc5v.used {
                    let seq = crate::variations::util::find_conseq_sets(
                        sc5v,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if !seq.is_empty() && seq.len() >= mnt_len {
                        let mut seq_rev = seq;
                        seq_rev.reverse();
                        if seq_rev.ends_with(mnt_bytes) {
                            let prefix_len = seq_rev.len().saturating_sub(mnt_len);
                            let head_matches = match seq_rev.get(..prefix_len) {
                                Some([]) => true,
                                Some(head) => is_match_ref(
                                    head,
                                    reference,
                                    chr,
                                    position.saturating_sub(1),
                                    -1,
                                    3,
                                ),
                                None => false,
                            };
                            if head_matches {
                                sc5_vars_to_add = Some(sc5v.variation.clone());
                            }
                        }
                    }
                }
            }
        }
        if let Some(sc5_vars) = sc5_vars_to_add {
            if let Some(vars_on_pos) = var_data.non_insertion_variants.get_mut(&position) {
                if let Some(vref_mut) = vars_on_pos.map.get_mut(&vn) {
                    crate::variations::util::adj_cnt(vref_mut, &sc5_vars, None);
                    vref = vref_mut.clone();
                }
            }
            inc_cnt(&mut var_data.ref_coverage, position, sc5_vars.cnt);
            if let Some(sc5v) = var_data.soft_clips_5_end.get_mut(&right_pos) {
                sc5v.used = true;
            }
        }
    }
}

/// Rescue some SNVs that are softly clipped due to being at the end of the read.
/// Ports `StructuralVariantsProcessor.java:1977-2027`.
pub fn adj_snv(
    reference: &Reference,
    chr: &str,
    _config: &crate::config::Config,
    var_data: &mut crate::variations::VariationData,
) {
    let adaptor_forward = AdaptorSet::default();
    let adaptor_reverse = AdaptorSet::default();

    let keys_5: Vec<i32> = var_data.soft_clips_5_end.keys().copied().collect();
    for position in keys_5 {
        let sc5v = match var_data.soft_clips_5_end.get_mut(&position) {
            Some(sc) => sc,
            None => continue,
        };
        if sc5v.used {
            continue;
        }
        let seq =
            crate::variations::util::find_conseq_sets(sc5v, 0, &adaptor_forward, &adaptor_reverse);
        if seq.len() > 5 {
            continue;
        }
        let bp_bytes = substr(&seq, 0, 1);
        let bp_lossy;
        let bp: &str = match std::str::from_utf8(bp_bytes) {
            Ok(s) => s,
            Err(_) => {
                bp_lossy = String::from_utf8_lossy(bp_bytes);
                bp_lossy.as_ref()
            }
        };

        let previous_position = position.saturating_sub(1);
        let mut should_adj = false;
        let sc5v_vars_count = sc5v.variation.cnt;

        if let Some(vars_on_prev) = var_data.non_insertion_variants.get_mut(&previous_position) {
            if vars_on_prev.map.contains_key(bp) {
                should_adj = true;
                if seq.len() > 1 {
                    if let Some(pos_sub_2) = position.checked_sub(2) {
                        let ref_base = reference.get_base(chr, pos_sub_2 as u32);
                        if ref_base.is_none() || ref_base != seq.get(1).copied() {
                            should_adj = false;
                        }
                    } else {
                        should_adj = false;
                    }
                }
                if should_adj {
                    let sc5v_var = sc5v.variation.clone();
                    if let Some(vref_mut) = vars_on_prev.map.get_mut(bp) {
                        crate::variations::util::adj_cnt(vref_mut, &sc5v_var, None);
                    }
                    inc_cnt(
                        &mut var_data.ref_coverage,
                        previous_position,
                        sc5v_vars_count,
                    );
                }
            }
        }
        if should_adj {
            if let Some(sc5v_mut) = var_data.soft_clips_5_end.get_mut(&position) {
                sc5v_mut.used = true;
            }
        }
    }

    let keys_3: Vec<i32> = var_data.soft_clips_3_end.keys().copied().collect();
    for position in keys_3 {
        let sc3v = match var_data.soft_clips_3_end.get_mut(&position) {
            Some(sc) => sc,
            None => continue,
        };
        if sc3v.used {
            continue;
        }
        let seq =
            crate::variations::util::find_conseq_sets(sc3v, 0, &adaptor_forward, &adaptor_reverse);
        if seq.len() > 5 {
            continue;
        }
        let bp_bytes = substr(&seq, 0, 1);
        let bp_lossy;
        let bp: &str = match std::str::from_utf8(bp_bytes) {
            Ok(s) => s,
            Err(_) => {
                bp_lossy = String::from_utf8_lossy(bp_bytes);
                bp_lossy.as_ref()
            }
        };

        let mut should_adj = false;
        let sc3v_vars_count = sc3v.variation.cnt;

        if let Some(vars_on_pos) = var_data.non_insertion_variants.get_mut(&position) {
            if vars_on_pos.map.contains_key(bp) {
                should_adj = true;
                if seq.len() > 1 {
                    if let Some(pos_add_1) = position.checked_add(1) {
                        let ref_base = reference.get_base(chr, pos_add_1 as u32);
                        if ref_base.is_none() || ref_base != seq.get(1).copied() {
                            should_adj = false;
                        }
                    } else {
                        should_adj = false;
                    }
                }
                if should_adj {
                    let sc3v_var = sc3v.variation.clone();
                    if let Some(vref_mut) = vars_on_pos.map.get_mut(bp) {
                        crate::variations::util::adj_cnt(vref_mut, &sc3v_var, None);
                    }
                    inc_cnt(&mut var_data.ref_coverage, position, sc3v_vars_count);
                }
            }
        }
        if should_adj {
            if let Some(sc3v_mut) = var_data.soft_clips_3_end.get_mut(&position) {
                sc3v_mut.used = true;
            }
        }
    }
}

#[inline]
fn is_atgc(b: u8) -> bool {
    matches!(b, b'A' | b'T' | b'G' | b'C')
}

#[inline]
fn is_atgnc(b: u8) -> bool {
    matches!(b, b'A' | b'T' | b'G' | b'N' | b'C')
}

fn parse_begin_minus_number(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    if bytes.first().copied() != Some(b'-') {
        return None;
    }
    let mut end = 1;
    while bytes.get(end).is_some_and(u8::is_ascii_digit) {
        end += 1;
    }
    if end == 1 {
        return None;
    }
    vn.get(1..end)
}

fn parse_up_number_end(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    let mut i = bytes.len();
    while i > 0 && bytes.get(i - 1).is_some_and(u8::is_ascii_digit) {
        i -= 1;
    }
    if i == bytes.len() || i == 0 || bytes.get(i - 1).copied() != Some(b'^') {
        return None;
    }
    vn.get(i..)
}

fn parse_minus_number_atgnc_sv_atgnc_end(vn: &str) -> Option<(&str, &str)> {
    let bytes = vn.as_bytes();
    let mut pos = if bytes.first().copied() == Some(b'-') {
        1
    } else {
        0
    };
    let digits_start = pos;
    while bytes.get(pos).is_some_and(u8::is_ascii_digit) {
        pos += 1;
    }
    if pos == digits_start || bytes.get(pos).copied() != Some(b'^') {
        return None;
    }
    pos += 1;
    let g1_start = pos;
    while bytes.get(pos).copied().is_some_and(is_atgnc) {
        pos += 1;
    }
    if pos == g1_start || bytes.get(pos).copied() != Some(b'<') {
        return None;
    }
    let g1_end = pos;
    pos += 1;

    let tail = vn.get(pos..)?;
    let mut chars = tail.chars();
    for _ in 0..3 {
        let ch = chars.next()?;
        if ch == '\n' {
            return None;
        }
        pos += ch.len_utf8();
    }

    let sv_digits_start = pos;
    while bytes.get(pos).is_some_and(u8::is_ascii_digit) {
        pos += 1;
    }
    if pos == sv_digits_start || bytes.get(pos).copied() != Some(b'>') {
        return None;
    }
    pos += 1;

    let g2_start = pos;
    while bytes.get(pos).copied().is_some_and(is_atgnc) {
        pos += 1;
    }
    if pos == g2_start || pos != bytes.len() {
        return None;
    }
    Some((vn.get(g1_start..g1_end)?, vn.get(g2_start..pos)?))
}

fn parse_begin_minus_number_any(vn: &str) -> Option<(&str, &str)> {
    let bytes = vn.as_bytes();
    if bytes.first().copied() != Some(b'-') {
        return None;
    }
    let mut pos = 1;
    while bytes.get(pos).is_some_and(u8::is_ascii_digit) {
        pos += 1;
    }
    if pos == 1 {
        return None;
    }
    let rest_end = bytes
        .get(pos..)
        .and_then(|s| s.iter().position(|&b| b == b'\n'))
        .map_or(bytes.len(), |off| pos + off);
    Some((vn.get(1..pos)?, vn.get(pos..rest_end)?))
}

fn parse_caret_atgnc(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    let mut i = 0;
    while let Some(&b) = bytes.get(i) {
        if b == b'^' && bytes.get(i + 1).copied().is_some_and(is_atgnc) {
            let start = i + 1;
            let mut end = start + 1;
            while bytes.get(end).copied().is_some_and(is_atgnc) {
                end += 1;
            }
            return vn.get(start..end);
        }
        i += 1;
    }
    None
}

fn parse_minus_number_amp_atgcs_end(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    let mut i = bytes.len();
    while i > 0 && bytes.get(i - 1).copied().is_some_and(is_atgc) {
        i -= 1;
    }
    if i == bytes.len() || i == 0 || bytes.get(i - 1).copied() != Some(b'&') {
        return None;
    }
    let amp_idx = i - 1;
    let mut d = amp_idx;
    while d > 0 && bytes.get(d - 1).is_some_and(u8::is_ascii_digit) {
        d -= 1;
    }
    if d == amp_idx || d == 0 || bytes.get(d - 1).copied() != Some(b'-') {
        return None;
    }
    vn.get(d - 1..amp_idx)
}

fn parse_begin_plus_atgc(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    if bytes.first().copied() != Some(b'+') {
        return None;
    }
    let mut end = 1;
    while bytes.get(end).copied().is_some_and(is_atgc) {
        end += 1;
    }
    if end == 1 {
        return None;
    }
    vn.get(1..end)
}

fn parse_dup_num_atgc(vn: &str) -> Option<(&str, &str)> {
    let bytes = vn.as_bytes();
    let mut i = bytes.len();
    while i > 0 && bytes.get(i - 1).copied().is_some_and(is_atgc) {
        i -= 1;
    }
    if i == bytes.len() || i == 0 || bytes.get(i - 1).copied() != Some(b'>') {
        return None;
    }
    let gt_idx = i - 1;
    let mut d = gt_idx;
    while d > 0 && bytes.get(d - 1).is_some_and(u8::is_ascii_digit) {
        d -= 1;
    }
    if d == gt_idx || d < 4 || bytes.get(d - 4..d) != Some(b"<dup") {
        return None;
    }
    Some((vn.get(d..gt_idx)?, vn.get(i..)?))
}

fn parse_amp_atgc(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    let mut i = 0;
    while let Some(&b) = bytes.get(i) {
        if b == b'&' && bytes.get(i + 1).copied().is_some_and(is_atgc) {
            let start = i + 1;
            let mut end = start + 1;
            while bytes.get(end).copied().is_some_and(is_atgc) {
                end += 1;
            }
            return vn.get(start..end);
        }
        i += 1;
    }
    None
}

fn parse_hash_atgc(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    let mut i = 0;
    while let Some(&b) = bytes.get(i) {
        if b == b'#' && bytes.get(i + 1).copied().is_some_and(is_atgc) {
            let start = i + 1;
            let mut end = start + 1;
            while bytes.get(end).copied().is_some_and(is_atgc) {
                end += 1;
            }
            return vn.get(start..end);
        }
        i += 1;
    }
    None
}

fn parse_atgs_amp_atgs_end(vn: &str) -> Option<&str> {
    let bytes = vn.as_bytes();
    let mut i = bytes.len();
    while i > 0 && bytes.get(i - 1).copied().is_some_and(is_atgc) {
        i -= 1;
    }
    if i == bytes.len() || i == 0 || bytes.get(i - 1).copied() != Some(b'&') {
        return None;
    }
    let amp_idx = i - 1;
    let mut d = amp_idx;
    while d > 0 && bytes.get(d - 1).copied().is_some_and(is_atgc) {
        d -= 1;
    }
    if d == amp_idx || d == 0 || bytes.get(d - 1).copied() != Some(b'+') {
        return None;
    }
    vn.get(d - 1..amp_idx)
}

fn strip_caret_digits_end(s: &str) -> &str {
    let bytes = s.as_bytes();
    let mut i = bytes.len();
    while i > 0 && bytes.get(i - 1).is_some_and(u8::is_ascii_digit) {
        i -= 1;
    }
    if i < bytes.len() && i > 0 && bytes.get(i - 1).copied() == Some(b'^') {
        s.get(..i - 1).unwrap_or(s)
    } else {
        s
    }
}

fn normalize_realignins_tn_into(tn_buf: &mut String, vn: &str) {
    tn_buf.clear();
    let rest = vn.strip_prefix('+').unwrap_or(vn);
    tn_buf.reserve(rest.len());
    for ch in rest.chars() {
        if ch != '&' && ch != '#' {
            tn_buf.push(ch);
        }
    }
    let trimmed_len = strip_caret_digits_end(tn_buf.as_str()).len();
    tn_buf.truncate(trimmed_len);
    tn_buf.retain(|ch| ch != '^');
}

#[cfg(test)]
mod tests {
    use crate::reference::LoadedRegion;
    use regex::Regex;
    use std::sync::LazyLock;

    static BEGIN_MINUS_NUMBER: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"^-(\d+)"));
    static UP_NUMBER_END: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"\^(\d+)$"));
    static MINUS_NUMBER_ATGNC_SV_ATGNC_END: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"^-?\d+\^([ATGNC]+)<...\d+>([ATGNC]+)$"));
    static BEGIN_MINUS_NUMBER_ANY: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"^-(\d+)(.*)"));
    static CARET_ATGNC: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"\^([ATGNC]+)"));
    static MINUS_NUMBER_AMP_ATGCS_END: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"(-\d+)&[ATGC]+$"));

    static BEGIN_PLUS_ATGC: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"^\+([ATGC]+)"));
    static DUP_NUM_ATGC: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"<dup(\d+)>([ATGC]+)$"));
    static AMP_ATGC: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"&([ATGC]+)"));
    static HASH_ATGC: LazyLock<Regex> = LazyLock::new(|| crate::util::static_regex(r"#([ATGC]+)"));
    static ATGS_AMP_ATGS_END: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"(\+[ATGC]+)&[ATGC]+$"));
    static RE_CARET_DIGITS_END: LazyLock<Regex> =
        LazyLock::new(|| crate::util::static_regex(r"\^\d+$"));

    /// W1.2 fix round 1: merge_sclip must union the SV state (mates, soft
    /// position counts, disc, start/end/mstart/mend min/maxes) exactly like
    /// Java's shared-entry accumulation in CigarParser.addSV (2295-2341).
    #[test]
    fn test_merge_sclip_sv_field_union() {
        use super::{merge_sclip, Sclip};
        use crate::variations::Mate;

        let mate = |ms: i32, me: i32| Mate {
            mate_start: ms,
            mate_end: me,
            ..Mate::new()
        };

        let mut dest = Sclip::new();
        dest.variation.cnt = 3;
        dest.start = 100;
        dest.end = 200;
        dest.mstart = 50;
        dest.mend = 250;
        dest.disc = 1;
        dest.soft.insert(95, 2);
        dest.mates.push(mate(50, 250));
        dest.sequence = Some(b"OLD".to_vec());

        let mut src = Sclip::new();
        src.variation.cnt = 2;
        src.start = 80; // smaller -> wins
        src.end = 300; // larger -> wins
        src.mstart = 40; // smaller -> wins
        src.mend = 400; // larger -> wins
        src.disc = 4;
        src.soft.insert(95, 1); // overlapping key sums
        src.soft.insert(105, 3); // disjoint key unions
        src.mates.push(mate(40, 400));
        src.mates.push(mate(60, 260));

        merge_sclip(&mut dest, src);

        assert_eq!(dest.variation.cnt, 5);
        assert_eq!(dest.start, 80);
        assert_eq!(dest.end, 300);
        assert_eq!(dest.mstart, 40);
        assert_eq!(dest.mend, 400);
        assert_eq!(dest.disc, 5);
        assert_eq!(dest.soft.get(&95), Some(&3));
        assert_eq!(dest.soft.get(&105), Some(&3));
        assert_eq!(dest.mates.len(), 3);
        assert_eq!(dest.mates[0].mate_start, 50); // dest mates stay first
                                                  // merged counts invalidate any cached consensus
        assert!(dest.sequence.is_none());

        // A source entry with no SV data must not clobber accumulated values.
        let mut dest2 = Sclip::new();
        dest2.start = 100;
        dest2.end = 200;
        merge_sclip(&mut dest2, Sclip::new());
        assert_eq!(dest2.start, 100);
        assert_eq!(dest2.end, 200);
    }

    #[test]
    fn test_adjust_mnp() {
        let mut reference = Reference::default();
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 100,
        });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);

        let mut var_data = crate::variations::VariationData::default();

        // Add candidate MNP
        var_data
            .mnp
            .entry(50)
            .or_default()
            .insert("T&A".to_string(), 5);

        // Add the MNP to non_insertion_variants
        let mut mnp_var = Variation::new();
        mnp_var.cnt = 5;
        mnp_var.pmean = 250.0; // mean = 50.0
        mnp_var.hicnt = 5;
        mnp_var.fwd = 5;
        var_data
            .non_insertion_variants
            .get_or_put_with(50, crate::variations::VariationMap::new)
            .map
            .insert("T&A".to_string(), mnp_var);

        // Add Left constituent (position 50, "T")
        let mut left_var = Variation::new();
        left_var.cnt = 2;
        left_var.pmean = 2.0; // mean position = 1.0 <= i + 1 (1.0)
        left_var.hicnt = 2;
        left_var.fwd = 2;
        var_data
            .non_insertion_variants
            .get_or_put_with(50, crate::variations::VariationMap::new)
            .map
            .insert("T".to_string(), left_var);

        // Add Right constituent (position 51, "A")
        let mut right_var = Variation::new();
        right_var.cnt = 2;
        right_var.pmean = 2.0;
        right_var.hicnt = 2;
        right_var.fwd = 2;
        var_data
            .non_insertion_variants
            .get_or_put_with(51, crate::variations::VariationMap::new)
            .map
            .insert("A".to_string(), right_var);

        let config = crate::config::Config::default();
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: String::new(),
            amplicon: None,
        };

        adjust_mnp(&reference, "chr1", 100, &region, &config, &mut var_data);

        // Verify MNP count was adjusted (5 + 2 + 2 = 9)
        let updated_mnp = var_data
            .non_insertion_variants
            .get(&50)
            .unwrap()
            .map
            .get("T&A")
            .unwrap();
        assert_eq!(updated_mnp.cnt, 9);
        assert_eq!(updated_mnp.fwd, 9);

        // Verify constituents were removed
        assert!(!var_data
            .non_insertion_variants
            .get(&50)
            .unwrap()
            .map
            .contains_key("T"));
        assert!(!var_data
            .non_insertion_variants
            .get(&51)
            .unwrap()
            .map
            .contains_key("A"));
    }

    /// Build a minimal fixture: MNP "C&T" (mnt = "CT") at position 50 plus an
    /// unused 3' soft clip at the same position whose consensus is `conseq`.
    fn mnp_sclip3_fixture(conseq: Option<&[u8]>) -> (Reference, crate::variations::VariationData) {
        let mut reference = Reference::default();
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 100,
        });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);

        let mut var_data = crate::variations::VariationData::default();

        var_data
            .mnp
            .entry(50)
            .or_default()
            .insert("C&T".to_string(), 1);

        let mut mnp_var = Variation::new();
        mnp_var.cnt = 1;
        mnp_var.hicnt = 1;
        var_data
            .non_insertion_variants
            .get_or_put_with(50, crate::variations::VariationMap::new)
            .map
            .insert("C&T".to_string(), mnp_var);

        if let Some(seq) = conseq {
            let mut sc = Sclip::new();
            sc.variation.cnt = 2;
            sc.sequence = Some(seq.to_vec());
            var_data.soft_clips_3_end.insert(50, sc);
        }

        (reference, var_data)
    }

    fn adjust_mnp_on(var_data: &mut crate::variations::VariationData, reference: &Reference) {
        let config = crate::config::Config::default();
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: String::new(),
            amplicon: None,
        };
        adjust_mnp(reference, "chr1", 100, &region, &config, var_data);
    }

    /// Ports the marking rule of `VariationRealigner.java:346-357`: when a 3'
    /// soft clip's consensus starts with the MNP description (`mnt`) and either
    /// has no remainder or the remainder matches the reference (`ismatchref`),
    /// adjustMNP folds the clip into the variant and marks it `used`, so that
    /// realigndel skips it. Regression for DIVERGENCES.md Category 5
    /// (threeIndels chr20:62358591).
    #[test]
    fn test_adjust_mnp_marks_3p_sclip_used() {
        // Exact-length consensus ("CT" == mnt): ismatchref arm short-circuits.
        let (reference, mut var_data) = mnp_sclip3_fixture(Some(b"CT"));
        adjust_mnp_on(&mut var_data, &reference);
        let sc = var_data.soft_clips_3_end.get(&50).unwrap();
        assert!(sc.used, "exact-length matching consensus must mark used");
        let vref = var_data
            .non_insertion_variants
            .get(&50)
            .unwrap()
            .map
            .get("C&T")
            .unwrap();
        assert_eq!(vref.cnt, 3, "clip reads must be folded into the variant");

        // Longer consensus whose remainder matches reference (all 'A').
        let (reference, mut var_data) = mnp_sclip3_fixture(Some(b"CTAAAA"));
        adjust_mnp_on(&mut var_data, &reference);
        assert!(var_data.soft_clips_3_end.get(&50).unwrap().used);
        let vref = var_data
            .non_insertion_variants
            .get(&50)
            .unwrap()
            .map
            .get("C&T")
            .unwrap();
        assert_eq!(vref.cnt, 3, "remainder-match must also fold clip counts");

        // Consensus not starting with mnt must stay unused.
        let (reference, mut var_data) = mnp_sclip3_fixture(Some(b"AA"));
        adjust_mnp_on(&mut var_data, &reference);
        let sc = var_data.soft_clips_3_end.get(&50).unwrap();
        assert!(!sc.used);
        let vref = var_data
            .non_insertion_variants
            .get(&50)
            .unwrap()
            .map
            .get("C&T")
            .unwrap();
        assert_eq!(vref.cnt, 1, "no folding when consensus does not match");

        // Consensus starting with mnt but whose remainder mismatches the
        // all-'A' reference must stay unused too (ismatchref fails).
        let (reference, mut var_data) = mnp_sclip3_fixture(Some(b"CTTTTT"));
        adjust_mnp_on(&mut var_data, &reference);
        let sc = var_data.soft_clips_3_end.get(&50).unwrap();
        assert!(!sc.used, "non-matching remainder must not mark used");
        let vref = var_data
            .non_insertion_variants
            .get(&50)
            .unwrap()
            .map
            .get("C&T")
            .unwrap();
        assert_eq!(vref.cnt, 1, "no folding when remainder mismatches");
    }

    #[test]
    fn test_adj_snv() {
        let mut reference = Reference::default();
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 100,
        });
        let mut ref_seq = vec![b'A'; 100];
        ref_seq[48] = b'C'; // position 49 (1-based index) is C
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);

        let mut var_data = crate::variations::VariationData::default();

        // Add soft clip 5' at position 50 (1-based index)
        let mut sc = Sclip::new();
        sc.variation.cnt = 3;
        sc.variation.fwd = 3;
        sc.sequence = Some(b"C".to_vec()); // length <= 5
        var_data.soft_clips_5_end.insert(50, sc);

        // Add SNV candidate bp ("C") at position 49 (position - 1)
        let mut snv = Variation::new();
        snv.cnt = 2;
        snv.fwd = 2;
        var_data
            .non_insertion_variants
            .get_or_put_with(49, crate::variations::VariationMap::new)
            .map
            .insert("C".to_string(), snv);

        let config = crate::config::Config::default();

        adj_snv(&reference, "chr1", &config, &mut var_data);

        // Verify count adjusted (2 + 3 = 5)
        let updated_snv = var_data
            .non_insertion_variants
            .get(&49)
            .unwrap()
            .map
            .get("C")
            .unwrap();
        assert_eq!(updated_snv.cnt, 5);

        // Verify soft clip is marked used
        assert!(var_data.soft_clips_5_end.get(&50).unwrap().used);
    }
    use super::*;

    #[test]
    fn test_ismatch() {
        assert!(is_match(b"ATGC", b"ATGC", 1, 3));
        assert!(is_match(b"ATGCATGCAT", b"ATGCATGCAC", 1, 3));
        assert!(!is_match(b"ATGCATGCAT", b"ATGCATGCCC", 1, 3));
        assert!(is_match(b"ATGC", b"CGTA", -1, 3));
        assert!(is_match(b"ATGC", b"C#G^TA", -1, 3));
    }

    #[test]
    fn test_ismatchref() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 10,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), b"ATGCATGCAT".to_vec());

        assert!(is_match_ref(b"ATG", &reference, "chr1", 1, 1, 3));
        assert!(is_match_ref(b"ATGCATGCAC", &reference, "chr1", 1, 1, 3));
        assert!(!is_match_ref(b"ATGCATGCCC", &reference, "chr1", 1, 1, 3));
        assert!(is_match_ref(b"TGC", &reference, "chr1", 8, -1, 3));
    }

    #[test]
    fn test_findbp() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        reference.reference_sequences.insert(
            "chr1".to_string(),
            b"ATGCATGCATGCATGCATGCATGCATGCAT".to_vec(),
        );

        let bp = find_bp(b"ATGCATGCAT", 1, &reference, 1, "chr1", 50, 30);
        assert_eq!(bp, 1);
    }

    #[test]
    fn test_findbi() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        reference.reference_sequences.insert(
            "chr1".to_string(),
            b"ATGCATGCATGCATGCATGCATGCATGCAT".to_vec(),
        );

        let bi = find_bi(b"AAAAAAATGCATGCATGC", 5, &reference, 1, "chr1", 30);
        assert_eq!(bi.base_insert, 4);
        assert_eq!(bi.insertion_sequence, "AAAAAA");
    }

    #[test]
    fn test_findmm() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        reference.reference_sequences.insert(
            "chr1".to_string(),
            b"ATGCATGCATGCATGCATGCATGCATGCAT".to_vec(),
        );

        let mut soft_clips_5_end = BTreeMap::new();
        let res5 = find_mm5(&reference, "chr1", 4, b"ATGG", &mut soft_clips_5_end);
        assert_eq!(res5.mismatches.len(), 1);
        assert_eq!(res5.mismatches[0].mismatch_sequence, "G");
        assert_eq!(res5.mismatches[0].mismatch_position, 4);

        let mut soft_clips_3_end = BTreeMap::new();
        let res3 = find_mm3(&reference, "chr1", 1, b"ATG", &mut soft_clips_3_end);
        assert_eq!(res3.mismatches.len(), 0);
    }

    #[test]
    fn test_realigndel() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        reference.reference_sequences.insert(
            "chr1".to_string(),
            b"AAAAAAAAATGCTTAAAAAAAAAAAAAAA".to_vec(),
        );

        let mut position_to_deletion_count = BTreeMap::new();
        let mut del_map = rustc_hash::FxHashMap::default();
        del_map.insert("-3".to_string(), 5);
        position_to_deletion_count.insert(10, del_map);

        let mut non_insertion_variants = crate::javamap::JavaHashMap::new();

        let mut map11 = VariationMap::new();
        map11.map.insert(
            "T".to_string(),
            crate::variations::Variation {
                cnt: 4,
                qmean: 100.0,
                pmean: 12.0,
                ..Default::default()
            },
        );
        non_insertion_variants.insert(11, map11);

        let mut map10 = VariationMap::new();
        map10.map.insert(
            "-3".to_string(),
            crate::variations::Variation {
                cnt: 5,
                qmean: 120.0,
                pmean: 50.0,
                ..Default::default()
            },
        );
        non_insertion_variants.insert(10, map10);

        let mut ref_coverage = rustc_hash::FxHashMap::default();
        ref_coverage.insert(10, 10);
        ref_coverage.insert(11, 10);

        let mut soft_clips_5_end = BTreeMap::new();
        let mut soft_clips_3_end = BTreeMap::new();

        let config = crate::config::Config::default();

        realigndel(
            None,
            &position_to_deletion_count,
            &reference,
            "chr1",
            30,
            &config,
            &mut non_insertion_variants,
            &mut ref_coverage,
            &mut soft_clips_5_end,
            &mut soft_clips_3_end,
            100,
        );

        // Verification:
        // Variant at 11 should be removed (since swap_remove is used and it was the only variant, entry is removed from map)
        assert!(!non_insertion_variants.contains_key(&11));

        // Variant "-3" at 10 should have its count increased from 5 to 9 (5 + 4)
        let entry10 = non_insertion_variants.get(&10).unwrap();
        let var10 = entry10.map.get("-3").unwrap();
        assert_eq!(var10.cnt, 9);
    }

    #[test]
    fn test_realignins() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        // Position 10: T, others A
        reference.reference_sequences.insert(
            "chr1".to_string(),
            b"AAAAAAAAATAAAAAAAAAAAAAAAAAAAA".to_vec(),
        );

        let mut position_to_insertion_count = BTreeMap::new();
        let mut ins_map = rustc_hash::FxHashMap::default();
        ins_map.insert("+C".to_string(), 5);
        position_to_insertion_count.insert(10, ins_map);

        let mut insertion_variants = crate::javamap::JavaHashMap::new();
        let mut map10_ins = VariationMap::new();
        map10_ins.map.insert(
            "+C".to_string(),
            crate::variations::Variation {
                cnt: 5,
                qmean: 120.0,
                pmean: 50.0,
                ..Default::default()
            },
        );
        insertion_variants.insert(10, map10_ins);

        let mut non_insertion_variants = crate::javamap::JavaHashMap::new();
        let mut map10_non = VariationMap::new();
        map10_non.map.insert(
            "C".to_string(),
            crate::variations::Variation {
                cnt: 4,
                qmean: 100.0,
                pmean: 12.0,
                ..Default::default()
            },
        );
        non_insertion_variants.insert(10, map10_non);

        let mut ref_coverage = rustc_hash::FxHashMap::default();
        ref_coverage.insert(10, 10);

        let mut soft_clips_5_end = BTreeMap::new();
        let mut soft_clips_3_end = BTreeMap::new();

        let config = crate::config::Config::default();

        realignins(
            None,
            &position_to_insertion_count,
            &reference,
            "chr1",
            30,
            &config,
            &mut insertion_variants,
            &mut non_insertion_variants,
            &mut ref_coverage,
            &mut soft_clips_5_end,
            &mut soft_clips_3_end,
            100,
        );

        // Verification:
        // Variant at 10 (mismatch "C") should be removed
        assert!(!non_insertion_variants.contains_key(&10));

        // Variant "+C" at 10 should have its count increased from 5 to 9 (5 + 4)
        let entry10 = insertion_variants.get(&10).unwrap();
        let var10 = entry10.map.get("+C").unwrap();
        assert_eq!(var10.cnt, 9);
    }

    #[test]
    fn test_lgdel() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        let mut ref_seq = vec![b'A'; 100];
        let seq = b"ACGTACGTACGTACGTA";
        ref_seq[9..26].copy_from_slice(seq);
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);
        reference.seed.insert(seq.to_vec(), vec![10]);

        let mut soft_clips_5_end = BTreeMap::new();
        let mut sc = Sclip::new();
        sc.variation.cnt = 5;
        let mut seq_rev = seq.to_vec();
        seq_rev.reverse();
        sc.sequence = Some(seq_rev);
        soft_clips_5_end.insert(50, sc);

        let mut non_insertion_variants = crate::javamap::JavaHashMap::new();
        let mut insertion_variants = crate::javamap::JavaHashMap::new();
        let mut ref_coverage = rustc_hash::FxHashMap::default();
        ref_coverage.insert(50, 10);
        let mut soft_clips_3_end = BTreeMap::new();

        let mut svfdel = Vec::new();
        let mut svrdel = Vec::new();
        let mut splice = rustc_hash::FxHashSet::default();

        let config = crate::config::Config::default();
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: String::new(),
            amplicon: None,
        };

        realignlgdel(
            None,
            &mut reference,
            "chr1",
            100,
            &region,
            &config,
            &mut non_insertion_variants,
            &mut insertion_variants,
            &mut ref_coverage,
            &mut soft_clips_5_end,
            &mut soft_clips_3_end,
            100,
            &mut svfdel,
            &mut svrdel,
            &mut splice,
        );

        // Verification:
        // Variant "-23" should be created at position 27
        let entry = non_insertion_variants.get(&27);
        assert!(entry.is_some(), "Variant map at 27 must be created");
        let var = entry.unwrap().map.get("-23");
        assert!(var.is_some(), "Deletion -23 must be created");
        assert_eq!(var.unwrap().cnt, 5);

        // Softclip 50 should be marked as used
        assert!(soft_clips_5_end.get(&50).unwrap().used);
    }

    #[test]
    fn test_lgins() {
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        let mut ref_seq = vec![b'A'; 100];
        let seq = b"ACGTACGTACGTACGTA";
        ref_seq[69..86].copy_from_slice(seq);
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);
        reference.seed.insert(seq.to_vec(), vec![70]);

        let mut soft_clips_5_end = BTreeMap::new();
        let mut sc = Sclip::new();
        sc.variation.cnt = 5;
        let mut seq_rev = seq.to_vec();
        seq_rev.reverse();
        sc.sequence = Some(seq_rev);
        soft_clips_5_end.insert(50, sc);

        let mut non_insertion_variants = crate::javamap::JavaHashMap::new();
        let mut insertion_variants = crate::javamap::JavaHashMap::new();
        let mut ref_coverage = rustc_hash::FxHashMap::default();
        ref_coverage.insert(50, 10);
        let mut soft_clips_3_end = BTreeMap::new();

        let mut svfdup = Vec::new();
        let mut svrdup = Vec::new();
        let mut splice = rustc_hash::FxHashSet::default();

        let config = crate::config::Config::default();
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: String::new(),
            amplicon: None,
        };

        realignlgins(
            None,
            &mut reference,
            "chr1",
            100,
            &region,
            &config,
            &mut non_insertion_variants,
            &mut insertion_variants,
            &mut ref_coverage,
            &mut soft_clips_5_end,
            &mut soft_clips_3_end,
            100,
            &mut svfdup,
            &mut svrdup,
            &mut splice,
        );

        // Verification: Softclip should be used
        assert!(soft_clips_5_end.get(&50).unwrap().used);
        assert!(
            !insertion_variants.is_empty(),
            "Insertion variants must be created"
        );
    }

    #[test]
    fn test_realigndel_softclip_ratio_uses_java_integer_division() {
        // ports VariationRealigner.java:561, 587:
        // `if (dcnt <= 2 && tv.varsCount / dcnt > 5) continue;`
        // 1. For dcnt = 2 and varsCount = 11, 11 / 2 == 5 (not > 5) -> rescued.
        // 2. For dcnt = 2 and varsCount = 12, 12 / 2 == 6 > 5 -> skipped.
        // 3. For dcnt = 0, checked_div returns None -> skipped without panic (AGENTS.md §9).
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        // 1..=9: 'A', 10..=12: 'T' (deleted -3 at p=10), 13..=30: 'C'
        let mut ref_seq = vec![b'A'; 9];
        ref_seq.extend_from_slice(b"TTT");
        ref_seq.extend_from_slice(&[b'C'; 18]);
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);
        let config = crate::config::Config::default();

        // Sub-case 1: cnt = 11, dcnt = 2 -> 11 / 2 == 5 (not > 5) -> rescued
        {
            let mut position_to_deletion_count = BTreeMap::new();
            let mut del_map = rustc_hash::FxHashMap::default();
            del_map.insert("-3".to_string(), 2);
            position_to_deletion_count.insert(10, del_map);

            let mut non_insertion_variants = crate::javamap::JavaHashMap::new();
            let mut map10 = VariationMap::new();
            map10.map.insert(
                "-3".to_string(),
                crate::variations::Variation {
                    cnt: 2,
                    qmean: 60.0,
                    pmean: 20.0,
                    ..Default::default()
                },
            );
            non_insertion_variants.insert(10, map10);

            let mut ref_coverage = rustc_hash::FxHashMap::default();
            ref_coverage.insert(10, 20);

            let mut soft_clips_5_end = BTreeMap::new();
            let mut sc5 = Sclip::new();
            sc5.variation.cnt = 11;
            sc5.sequence = Some(b"AAAAAA".to_vec());
            soft_clips_5_end.insert(13, sc5);

            let mut soft_clips_3_end = BTreeMap::new();
            let mut sc3 = Sclip::new();
            sc3.variation.cnt = 11;
            sc3.sequence = Some(b"CCCCCC".to_vec());
            soft_clips_3_end.insert(10, sc3);

            realigndel(
                None,
                &position_to_deletion_count,
                &reference,
                "chr1",
                30,
                &config,
                &mut non_insertion_variants,
                &mut ref_coverage,
                &mut soft_clips_5_end,
                &mut soft_clips_3_end,
                100,
            );

            assert!(
                soft_clips_5_end.get(&13).unwrap().used,
                "5' softclip with cnt=11, dcnt=2 has 11/2 == 5 (not > 5) and must be rescued"
            );
            assert!(
                soft_clips_3_end.get(&10).unwrap().used,
                "3' softclip with cnt=11, dcnt=2 has 11/2 == 5 (not > 5) and must be rescued"
            );
            let var10 = non_insertion_variants
                .get(&10)
                .unwrap()
                .map
                .get("-3")
                .unwrap();
            assert_eq!(var10.cnt, 2 + 11 + 11);
        }

        // Sub-case 2: cnt = 12, dcnt = 2 -> 12 / 2 == 6 > 5 -> skipped
        // Sub-case 3: cnt = 4, dcnt = 0 -> checked_div is None -> skipped without panic
        for (dcnt, clip_cnt) in [(2, 12), (0, 4)] {
            let mut position_to_deletion_count = BTreeMap::new();
            let mut del_map = rustc_hash::FxHashMap::default();
            del_map.insert("-3".to_string(), dcnt);
            position_to_deletion_count.insert(10, del_map);

            let mut non_insertion_variants = crate::javamap::JavaHashMap::new();
            let mut map10 = VariationMap::new();
            map10.map.insert(
                "-3".to_string(),
                crate::variations::Variation {
                    cnt: 2,
                    qmean: 60.0,
                    pmean: 20.0,
                    ..Default::default()
                },
            );
            non_insertion_variants.insert(10, map10);

            let mut ref_coverage = rustc_hash::FxHashMap::default();
            ref_coverage.insert(10, 20);

            let mut soft_clips_5_end = BTreeMap::new();
            let mut sc5 = Sclip::new();
            sc5.variation.cnt = clip_cnt;
            sc5.sequence = Some(b"AAAAAA".to_vec());
            soft_clips_5_end.insert(13, sc5);

            let mut soft_clips_3_end = BTreeMap::new();
            let mut sc3 = Sclip::new();
            sc3.variation.cnt = clip_cnt;
            sc3.sequence = Some(b"CCCCCC".to_vec());
            soft_clips_3_end.insert(10, sc3);

            realigndel(
                None,
                &position_to_deletion_count,
                &reference,
                "chr1",
                30,
                &config,
                &mut non_insertion_variants,
                &mut ref_coverage,
                &mut soft_clips_5_end,
                &mut soft_clips_3_end,
                100,
            );

            assert!(!soft_clips_5_end.get(&13).unwrap().used);
            assert!(!soft_clips_3_end.get(&10).unwrap().used);
            let var10 = non_insertion_variants
                .get(&10)
                .unwrap()
                .map
                .get("-3")
                .unwrap();
            assert_eq!(var10.cnt, 2);
        }
    }

    #[test]
    fn test_lgdel_dellen_zero_includes_extra_len_in_variant_key() {
        // ports VariationRealigner.java:1078-1079:
        // `dellen -= EXTRA.length();`
        // `gt = dellen == 0 ? "-" + EXTRA.length() + "^" + EXTRA : "-" + dellen + "&" + EXTRA;`
        // Reachable when `p - bp_initial == 16` (passes `p - bp > 15` at :1028),
        // `bp++` makes `dellen = p - bp = 15`, and `findMatch` returns a 15-bp
        // `EXTRA` (`ii = 15`), yielding `dellen == 0`.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        let mut ref_seq = vec![b'A'; 100];
        let seed_seq = b"ACGTACGTACGTACGTA"; // 17 bp placed at 1-based 28..=44
        ref_seq[27..44].copy_from_slice(seed_seq);
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);
        reference.seed.insert(seed_seq.to_vec(), vec![28]);

        // 5' softclip at p = 60 with cnt = 5.
        // Forward orientation (seq_vec after reverse in find_match with dir=-1):
        // 17-bp seed_seq (matching 28..=44) + 15-bp EXTRA of 'C's (mismatching 'A's at 45..=59).
        // Initial bp in find_match = 28 + 32 - 0 - 1 = 59.
        // At ii = 15: bp_mut = 59 - 15 = 44, sseq = seed_seq (matches 28..=44),
        // EXTRA = "CCCCCCCCCCCCCCC" (len 15).
        // In realignlgdel: p - bp = 60 - 44 = 16 > 15; bp++ -> 45;
        // dellen = 60 - 45 = 15; dellen -= 15 -> 0.
        let extra_15 = b"CCCCCCCCCCCCCCC";
        let mut forward_clip = seed_seq.to_vec();
        forward_clip.extend_from_slice(extra_15);
        let mut sc_seq = forward_clip;
        sc_seq.reverse(); // stored reversed in 5' Sclip

        let mut soft_clips_5_end = BTreeMap::new();
        let mut sc = Sclip::new();
        sc.variation.cnt = 5;
        sc.sequence = Some(sc_seq);
        soft_clips_5_end.insert(60, sc);

        let mut non_insertion_variants = crate::javamap::JavaHashMap::new();
        let mut insertion_variants = crate::javamap::JavaHashMap::new();
        let mut ref_coverage = rustc_hash::FxHashMap::default();
        ref_coverage.insert(60, 10);
        let mut soft_clips_3_end = BTreeMap::new();
        let mut svfdel = Vec::new();
        let mut svrdel = Vec::new();
        let mut splice = rustc_hash::FxHashSet::default();

        let mut config = crate::config::Config::default();
        config.indel_size = 5; // ensure find_bp returns 0 so find_match runs
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: String::new(),
            amplicon: None,
        };

        realignlgdel(
            None,
            &mut reference,
            "chr1",
            100,
            &region,
            &config,
            &mut non_insertion_variants,
            &mut insertion_variants,
            &mut ref_coverage,
            &mut soft_clips_5_end,
            &mut soft_clips_3_end,
            100,
            &mut svfdel,
            &mut svrdel,
            &mut splice,
        );

        let entry45 = non_insertion_variants
            .get(&45)
            .expect("variant at bp=45 must be created");
        assert!(
            entry45.map.contains_key("-15^CCCCCCCCCCCCCCC"),
            "expected Java key '-15^CCCCCCCCCCCCCCC', got keys: {:?}",
            entry45.map.keys().collect::<Vec<_>>()
        );
    }

    /// Mid-query failure must skip the whole BAM, exactly like the old
    /// `query_region`-Vec path (which failed atomically before any record
    /// from that BAM was counted). Without the counter snapshot, the
    /// spanning reads counted before the failure flip the result from true
    /// to false and change downstream deletion realignment.
    #[test]
    fn no_passing_reads_midquery_failure_skips_whole_bam() {
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

        fn write_indexed_bam(path: &std::path::Path, recs: &[RecordBuf]) {
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
                .set_length(10000usize.try_into().unwrap())
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
            let mut bai_path = std::path::PathBuf::from(path);
            bai_path.set_extension("bam.bai");
            bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("span.bam");
        // 3000 spanning reads (start 50, 100M) over [100,110]: ~450 KB
        // uncompressed across several BGZF blocks, so truncating at 1/2
        // keeps the header block intact while a later record read fails.
        let recs: Vec<RecordBuf> = (0..3000)
            .map(|i| {
                record_buf(
                    &format!("s{i:05}"),
                    50,
                    &[(100, Kind::Match)],
                    vec![b'A'; 100],
                )
            })
            .collect();
        write_indexed_bam(&bam_path, &recs);
        let bam = bam_path.to_string_lossy().into_owned();
        let config = crate::config::Config::default();

        // Intact: spanning reads exist -> passing reads found.
        assert!(
            !super::no_passing_reads("chr1", 100, 110, std::slice::from_ref(&bam), &config),
            "intact BAM has spanning reads"
        );

        // Truncate mid-stream (index copied over, chunks point past EOF).
        let trunc_path = dir.path().join("trunc.bam");
        std::fs::copy(&bam_path, &trunc_path).unwrap();
        std::fs::copy(
            bam_path.with_extension("bam.bai"),
            trunc_path.with_extension("bam.bai"),
        )
        .unwrap();
        {
            let len = std::fs::metadata(&trunc_path).unwrap().len();
            let f = std::fs::OpenOptions::new()
                .write(true)
                .open(&trunc_path)
                .unwrap();
            f.set_len(len / 2).unwrap();
        }
        let trunc = trunc_path.to_string_lossy().into_owned();

        // Old Vec-path premise: the truncated query fails atomically.
        assert!(
            crate::bam::query_region(
                std::path::Path::new(&trunc),
                "chr1",
                crate::types::Pos0(99),
                crate::types::Pos0(110),
                &config,
            )
            .is_err(),
            "truncated BAM must fail the Vec path"
        );

        // Fixed streaming path: the failed BAM is skipped whole, matching
        // the old atomic behavior (pre-fix this returned false from the
        // retained spanning-read prefix).
        assert!(
            super::no_passing_reads("chr1", 100, 110, &[trunc], &config),
            "mid-query failure must skip the whole BAM"
        );
    }

    #[test]
    fn test_dup_ins3_flank_substring_offset_matches_java_substr() {
        // Bug F (`VariationRealigner.java:744-746`):
        // `if (ins3.length() > Configuration.SVFLANK) { sanpseq = substr(ins3, Configuration.SVFLANK - ins3.length()); }`
        // Because `SVFLANK - ins3.length() < 0`, `Utils.substr(ins3, neg)` starts at
        // `ins3.length() + (SVFLANK - ins3.length()) == SVFLANK` (50), matching
        // `super::substr(ins3, 50 - ins3.len() as isize, None)`.
        let sv_flank = 50usize;
        let ins3 = "A".repeat(50) + &"C".repeat(35) + &"G".repeat(40);
        let java_substr = super::substr(
            ins3.as_bytes(),
            (sv_flank as isize) - (ins3.len() as isize),
            ins3.len() as isize,
        );
        assert_eq!(java_substr, &ins3.as_bytes()[sv_flank..]);

        // Also exercise `realignins` directly on a `<dup100>` variant with `ins3.len() = 55 > 50`:
        // `ins3` has 50 'A's followed by "CAAAA", so `sanpseq = ins3[50..]` starts with 'C' at
        // `p3 + 1 = 253` (where `p = 100`, `insert = "AT"`, `dup_len = 100`, `ins3.len() = 55`,
        // `inslen = 157`, `p3 = 100 + 157 - 55 + 50 = 252`). With `start_idx = 50`, `find_mm3`
        // detects the 'C' mismatch at pos 253 and `realignins` absorbs `non_insertion_variants[253]["C"]`
        // into `insertion_variants[100][vn]`.
        let mut reference = crate::reference::Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 500,
            });
        reference.region_sequences.push(vec![b'A'; 500]);

        let ins3_dup = format!("{}CAAAA", "A".repeat(50));
        let vn = format!("+AT<dup100>{}", ins3_dup);
        let mut pos_to_ins = std::collections::BTreeMap::new();
        let mut inner = rustc_hash::FxHashMap::default();
        inner.insert(vn.clone(), 3);
        pos_to_ins.insert(100, inner);

        let mut insertion_variants = crate::javamap::JavaHashMap::new();
        let ins_v = super::get_variation(&mut insertion_variants, 100, &vn);
        ins_v.cnt = 3;
        ins_v.qmean = 90.0;
        ins_v.pmean = 30.0;

        let mut non_insertion_variants = crate::javamap::JavaHashMap::new();
        let snv = super::get_variation(&mut non_insertion_variants, 253, "C");
        snv.cnt = 1;
        snv.qmean = 35.0;
        snv.pmean = 2.0;

        let mut ref_coverage = rustc_hash::FxHashMap::default();
        let mut sc5 = std::collections::BTreeMap::new();
        let mut sc3 = std::collections::BTreeMap::new();
        let config = crate::config::Config::default();

        super::realignins(
            None,
            &pos_to_ins,
            &reference,
            "chr1",
            500,
            &config,
            &mut insertion_variants,
            &mut non_insertion_variants,
            &mut ref_coverage,
            &mut sc5,
            &mut sc3,
            100,
        );

        let updated_cnt = insertion_variants
            .get(&100)
            .and_then(|m| m.map.get(&vn))
            .map(|v| v.cnt)
            .unwrap_or(0);
        assert_eq!(
            updated_cnt, 4,
            "realignins must absorb the 3' flank mismatch at pos 253 when sanpseq starts at offset 50"
        );
        assert!(
            !non_insertion_variants.contains_key(&253),
            "absorbed mismatch at pos 253 must be removed from non_insertion_variants"
        );
    }

    #[test]
    fn test_ensure_reference_loaded_appends_secondary_region_and_populates_seeds() {
        // Ports VariationRealigner.java:149, 250, 916, 1020:
        // When `realignlgdel` or `realignlgins` loads an out-of-window mate region,
        // Java calls `ReferenceResource.getReference(modifiedRegion, maxReadLength, ...)`
        // which appends a second `LoadedRegion` (preserving the primary region's bounds)
        // and populates `reference.seed` for the newly loaded interval.
        let pid = std::process::id();
        let fasta_path = std::env::temp_dir().join(format!("vr_ensure_ref_{pid}.fa"));
        let fai_path = std::env::temp_dir().join(format!("vr_ensure_ref_{pid}.fa.fai"));
        let seq_line = format!("{}{}\n", "A".repeat(350), "CGTACGTACGTACGTAC".repeat(15));
        let chr_len = seq_line.trim().len();
        std::fs::write(&fasta_path, format!(">chr1\n{}", seq_line)).unwrap();
        std::fs::write(
            &fai_path,
            format!("chr1\t{}\t6\t{}\t{}\n", chr_len, chr_len, chr_len + 1),
        )
        .unwrap();

        let mut config = crate::config::Config::default();
        config.reference = fasta_path.to_string_lossy().into_owned();

        let mut reference = crate::reference::Reference::new();
        let primary_region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(100),
            gene: String::new(),
            amplicon: None,
        };
        let mut chr_lengths = std::collections::HashMap::new();
        chr_lengths.insert("chr1".to_string(), chr_len);
        crate::reference::get_reference(&primary_region, 50, &chr_lengths, &config, &mut reference)
            .unwrap();
        assert_eq!(reference.loaded_regions.len(), 1);
        let primary_bounds = (
            reference.loaded_regions[0].sequence_start,
            reference.loaded_regions[0].sequence_end,
        );

        super::ensure_reference_loaded(
            &mut reference,
            "chr1",
            400,
            450,
            50,
            chr_len as u32,
            &config,
        )
        .unwrap();

        let _ = std::fs::remove_file(&fasta_path);
        let _ = std::fs::remove_file(&fai_path);

        assert_eq!(
            reference.loaded_regions.len(),
            2,
            "ensure_reference_loaded must append a secondary LoadedRegion via get_reference rather than overwriting loaded_regions[0]"
        );
        assert_eq!(
            (
                reference.loaded_regions[0].sequence_start,
                reference.loaded_regions[0].sequence_end
            ),
            primary_bounds,
            "primary LoadedRegion bounds must be preserved"
        );
    }

    #[test]
    fn realignlgins_small_max_read_length_does_not_underflow() {
        // J1: When max_read_length < 10 (e.g. 0 on pure =/X reads),
        // `(ins.len() as i32) < max_read_length - 10` must evaluate to false in signed i32
        // without wrapping `usize` in release or panicking in debug.
        use crate::config::Config;
        use crate::javamap::JavaHashMap;
        use crate::reference::{LoadedRegion, Reference};
        use crate::regions::Region;
        use crate::types::Pos1;
        use crate::variations::{Sclip, Variation, VariationMap};

        let chr = "chr1";
        let mut bases: Vec<u8> = (0..300)
            .map(|i| b"ACGT"[(i * 7 + (i / 3) * 5 + (i / 11)) % 4])
            .collect();
        // Make 101..=120 non-homopolymer repeating unit so rpflag can be true
        let motif = b"ACGTTGCAACGTTGCAACGT";
        bases[100..120].copy_from_slice(motif);
        bases[120..140].copy_from_slice(motif);

        let mut reference = Reference::new();
        reference.loaded_regions.push(LoadedRegion {
            chr: chr.to_string(),
            sequence_start: 1,
            sequence_end: 300,
        });
        reference
            .reference_sequences
            .insert(chr.to_string(), bases.clone());
        reference.region_sequences.push(bases.clone());

        let region = Region {
            chr: chr.to_string(),
            start: Pos1(1),
            end: Pos1(300),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut config = Config::default();
        config.min_reads = 1;

        for (max_rl, expected_ref_cnt) in [(30i32, 0i32), (0, 1), (5, 1), (9, 1)] {
            let mut non_ins: JavaHashMap<i32, VariationMap> = JavaHashMap::new();
            let mut ins_vars: JavaHashMap<i32, VariationMap> = JavaHashMap::new();
            let mut ref_cov = rustc_hash::FxHashMap::default();
            let mut sc5 = std::collections::BTreeMap::new();
            let mut sc3 = std::collections::BTreeMap::new();
            let mut svfdup = Vec::new();
            let mut svrdup = Vec::new();
            let mut splice = rustc_hash::FxHashSet::default();

            // Populate a 5' soft-clip at p = 101 whose consensus matches motif + diverse flanking reference
            let mut sclip = Sclip::new();
            sclip.variation.cnt = 5;
            sclip.variation.hicnt = 5;
            sclip.variation.pmean = 100.0;
            sclip.variation.qmean = 150.0;
            // 5' sclip seq is stored reversed from breakpoint p=101 going left:
            // construct a 21-bp clip where find_bi (n >= 6) finds a 6-bp insertion at bi=100
            let mut clip_seq = Vec::new();
            clip_seq.extend_from_slice(&bases[85..100]); // 15 bp diverse flank matching 86..=100
            clip_seq.extend_from_slice(&motif[..6]); // 6 bp inserted repeat matching 101..=106
            for (idx, &b) in clip_seq.iter().rev().enumerate() {
                let mut v = Variation::new();
                v.cnt = 5;
                v.hicnt = 5;
                v.qmean = 150.0;
                sclip.seq.entry(idx as i32).or_default().insert(b, v);
                *sclip
                    .nt
                    .entry(idx as i32)
                    .or_default()
                    .entry(b)
                    .or_insert(0) += 5;
            }
            sc5.insert(101, sclip);

            // Reference variant at bi=100 with cnt=1
            let mref = crate::variations::get_variation(
                &mut non_ins,
                100,
                &(bases[99] as char).to_string(),
            );
            mref.cnt = 1;
            mref.hicnt = 1;

            let dummy_bams = vec!["nonexistent.bam".to_string()];
            super::realignlgins(
                Some(&dummy_bams),
                &mut reference,
                chr,
                300,
                &region,
                &config,
                &mut non_ins,
                &mut ins_vars,
                &mut ref_cov,
                &mut sc5,
                &mut sc3,
                max_rl,
                &mut svfdup,
                &mut svrdup,
                &mut splice,
            );
            assert!(
                ins_vars
                    .get(&100)
                    .is_some_and(|m| m.map.contains_key("+ACGTTG")),
                "expected +ACGTTG insertion candidate at bi=100 for max_read_length={max_rl}"
            );
            // When max_rl=30, (ins.len() as i32) < max_rl - 10 (6 < 20) is true so mref is absorbed (cnt=0).
            // When max_rl in [0, 5, 9], signed (ins.len() as i32) < max_rl - 10 is false so mref stays 1.
            let ref_cnt = non_ins
                .get(&100)
                .and_then(|m| m.map.get(&(bases[99] as char).to_string()))
                .map(|v| v.cnt)
                .unwrap_or(0);
            assert_eq!(
                ref_cnt, expected_ref_cnt,
                "unexpected mref cnt at bi=100 for max_read_length={max_rl}"
            );
        }
    }

    #[test]
    fn fill_adaptor_maps_java_seed_count() {
        // I2: VarDictLauncher.java:108 uses strict `< sequence.length()` (`i + 6 < len`),
        // and SequenceUtil.complement complements only ACGTacgt (leaving N unchanged).
        let (fwd6, rev6) = super::fill_adaptor_maps(&["ACGTAC".to_string()]);
        assert!(fwd6.is_empty());
        assert!(rev6.is_empty());

        let (fwd7, rev7) = super::fill_adaptor_maps(&["ACGNTAC".to_string()]);
        assert_eq!(fwd7.len(), 1);
        assert_eq!(fwd7.get(b"ACGNTA".as_slice()), Some(&1));
        assert_eq!(rev7.get(b"TANCGT".as_slice()), Some(&1));
    }

    #[test]
    fn partial_pipeline_merges_splice() {
        use noodles::core::Position;
        use noodles::sam::alignment::io::Write as AlignmentWrite;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::{
            Cigar, Data, QualityScores, RecordBuf, Sequence,
        };
        use std::num::NonZeroUsize;

        let pid = std::process::id();
        let bam_path = std::env::temp_dir().join(format!("vr_partial_splice_{pid}.bam"));
        let header = noodles::sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                noodles::sam::header::record::value::map::Map::<
                    noodles::sam::header::record::value::map::ReferenceSequence,
                >::new(NonZeroUsize::new(200).unwrap()),
            )
            .build();
        {
            use noodles::csi::binning_index::index::reference_sequence::{
                bin::Chunk, index::LinearIndex,
            };
            use noodles::csi::binning_index::Indexer;
            let file = std::fs::File::create(&bam_path).unwrap();
            let mut writer = noodles::bam::io::Writer::new(file);
            writer.write_header(&header).unwrap();
            // 10M20N10M starting at 1-based pos 10 -> splice "19-40"
            let cigar: Cigar = [
                Op::new(Kind::Match, 10),
                Op::new(Kind::Skip, 20),
                Op::new(Kind::Match, 10),
            ]
            .into_iter()
            .collect();
            let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(0i32))]
                .into_iter()
                .collect();
            let rec = RecordBuf::builder()
                .set_name("r_splice")
                .set_flags(noodles::sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(10usize).unwrap())
                .set_mapping_quality(
                    noodles::sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                )
                .set_cigar(cigar)
                .set_sequence(Sequence::from(b"AAAAAAAAAAAAAAAAAAAA".to_vec()))
                .set_quality_scores(QualityScores::from(vec![35u8; 20]))
                .set_data(data)
                .build();
            let start_vp = writer.get_ref().virtual_position();
            writer.write_alignment_record(&header, &rec).unwrap();
            let end_vp = writer.get_ref().virtual_position();
            writer.try_finish().unwrap();

            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
            let s = Position::try_from(10usize).unwrap();
            let e = Position::try_from(49usize).unwrap();
            indexer
                .add_record(Some((0, s, e, true)), Chunk::new(start_vp, end_vp))
                .unwrap();
            let bai_path = bam_path.with_extension("bam.bai");
            noodles::bam::bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }

        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 200]);

        let config = crate::config::Config::default();
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: String::new(),
            amplicon: None,
        };
        let mut non_ins = crate::javamap::JavaHashMap::new();
        let mut ins = crate::javamap::JavaHashMap::new();
        let mut ref_cov = rustc_hash::FxHashMap::default();
        let mut sc5 = BTreeMap::new();
        let mut sc3 = BTreeMap::new();
        let mut splice = rustc_hash::FxHashSet::default();
        let bams = vec![bam_path.to_string_lossy().into_owned()];

        super::partial_pipeline(
            Some(&bams),
            "chr1",
            1,
            100,
            &reference,
            &config,
            &region,
            &mut non_ins,
            &mut ins,
            &mut ref_cov,
            &mut sc5,
            &mut sc3,
            &mut splice,
        )
        .unwrap();

        let _ = std::fs::remove_file(&bam_path);
        let _ = std::fs::remove_file(bam_path.with_extension("bam.bai"));

        assert!(
            splice.contains("19-39"),
            "partial_pipeline must merge parser.splice into caller's splice set, got: {:?}",
            splice
        );
    }

    #[test]
    fn find_bi_non_ascii_consensus_no_panic() {
        let mut reference = crate::reference::Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 60,
            });
        let ref_seq: Vec<u8> = b"ACGTACGTAC".iter().copied().cycle().take(60).collect();
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq.clone());

        // Build a consensus `seq` where `n = 6` has non-ASCII bytes `[0xc9, 0x80, ...]`
        // and the next 10 bases match `ref_pos = 30 - i - 1` (`dir = -1`) with 1 mismatch
        // at `i = 0` (`mm = 1`, `i = 10`, `mnt = 4 >= 3`, `ratio = 0.10 < 0.15`),
        // exercising `dir == -1` and `!extra.is_empty()` at `realign.rs:284-290`.
        let mut consensus = vec![0xc9, 0x80, b'A', b'C', b'G', b'T'];
        for i in 0..10usize {
            let ref_pos = (30 - (i as i32) - 1) as usize;
            let ref_base = ref_seq[ref_pos - 1];
            if i == 0 {
                consensus.push(if ref_base == b'A' { b'T' } else { b'A' });
            } else {
                consensus.push(ref_base);
            }
        }
        let res = super::find_bi(&consensus, 30, &reference, -1, "chr1", 60);
        assert!(
            res.insertion_sequence.contains('&'),
            "expected '&' in insertion_sequence, got: {}",
            res.insertion_sequence
        );

        // `join_ref_for_3_lgins` (`VariationUtils.java:236-255`):
        // 1) When `extra.len() >= seq.len()`, must NOT underflow `usize` and must NOT return early:
        //    the `for (int i = from; i <= to; i++)` loop falls through to `else` and emits reference bases.
        let fallback_gt =
            super::join_ref_for_3_lgins(&reference, "chr1", 11, 16, 0, b"GG", b"TTTTT");
        assert_eq!(fallback_gt, &ref_seq[10..16]);
        let fallback_eq = super::join_ref_for_3_lgins(&reference, "chr1", 11, 16, 2, b"GG", b"TT");
        assert_eq!(fallback_eq, &ref_seq[10..16]);

        // 2) When `0 < extra.len() < seq.len()` and `shift5 > 0`, emits ref for `diff < shift5`,
        //    `seq[diff - shift5 + extra.len()]` while `diff - shift5 < seq.len() - extra.len()`,
        //    and ref afterward.
        let mixed = super::join_ref_for_3_lgins(&reference, "chr1", 11, 16, 1, b"XXGGG", b"XX");
        let mut expected_mixed = vec![ref_seq[10]];
        expected_mixed.extend_from_slice(b"GGG");
        expected_mixed.extend_from_slice(&ref_seq[14..16]);
        assert_eq!(mixed, expected_mixed);

        // 3) Near-`i32::MAX` bounds (`to == i32::MAX`, `chr_len == u32::MAX`, `region.end == i32::MAX`)
        //    must not overflow `i32` loop counters or window arithmetic.
        let mut max_ref = crate::reference::Reference::new();
        max_ref.loaded_regions.push(crate::reference::LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: (i32::MAX - 5) as u32,
            sequence_end: i32::MAX as u32,
        });
        max_ref.region_sequences.push(b"TGCAAC".to_vec());
        assert_eq!(
            super::join_ref_for_3_lgins(&max_ref, "chr1", i32::MAX - 3, i32::MAX, 0, b"G", b"TTTT",),
            b"CAAC"
        );
        assert_eq!(
            super::join_ref_for_5_lgins(&max_ref, "chr1", i32::MAX - 3, i32::MAX, b"G", b"TTTT"),
            b"CAAC"
        );
        let _ = super::find_bi(
            b"ACGTACGTACGTACGT",
            i32::MAX - 2,
            &max_ref,
            1,
            "chr1",
            u32::MAX,
        );
        let _ = super::find_bp(
            b"ACGTACGTACGT",
            i32::MAX - 2,
            &max_ref,
            1,
            "chr1",
            50,
            u32::MAX,
        );
        let max_region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1((i32::MAX - 5) as u32),
            end: Pos1(i32::MAX as u32),
            gene: String::new(),
            amplicon: None,
        };
        let mut max_var_data = crate::variations::VariationData::default();
        max_var_data.max_read_length = 150;
        max_var_data
            .position_to_deletion_count
            .entry(i32::MAX - 2)
            .or_default()
            .insert("-10".to_string(), 5);
        max_var_data
            .position_to_insertion_count
            .entry(i32::MAX - 2)
            .or_default()
            .insert("+ACGTACGT".to_string(), 5);
        let mut sc_max = super::Sclip::new();
        sc_max.variation.cnt = 5;
        for idx in 0..14 {
            let mut nt = BTreeMap::new();
            nt.insert(b'A', 10);
            sc_max.nt.insert(idx, nt);
        }
        max_var_data
            .soft_clips_5_end
            .insert(i32::MAX - 1, sc_max.clone());
        max_var_data.soft_clips_3_end.insert(i32::MAX - 1, sc_max);
        super::realign_indels(
            None,
            &mut max_ref,
            "chr1",
            u32::MAX,
            &max_region,
            &crate::config::Config::default(),
            &mut max_var_data,
        );

        // `realignlgdel` 5'-lgdel homopolymer shift reaching `bp + dellen + n == i32::MAX`
        // must terminate cleanly (`checked_add` -> `None`) at `n == 2` instead of clamping
        // onto `i32::MAX` and looping indefinitely.
        {
            let mut lgdel_ref = crate::reference::Reference::new();
            lgdel_ref
                .loaded_regions
                .push(crate::reference::LoadedRegion {
                    chr: "chr1".to_string(),
                    sequence_start: (i32::MAX - 25) as u32,
                    sequence_end: i32::MAX as u32,
                });
            // `[i32::MAX - 25 ..= i32::MAX - 16]` = `TCAGTCAGTC` (10 bp),
            // `[i32::MAX - 15 ..= i32::MAX]` = 16 `'A'`s.
            lgdel_ref
                .region_sequences
                .push(b"TCAGTCAGTCAAAAAAAAAAAAAAAA".to_vec());

            let lgdel_region = crate::regions::Region {
                chr: "chr1".to_string(),
                start: Pos1((i32::MAX - 25) as u32),
                end: Pos1(i32::MAX as u32),
                gene: String::new(),
                amplicon: None,
            };
            let mut lgdel_sc5 = BTreeMap::new();
            let mut sc5 = super::Sclip::new();
            sc5.variation.cnt = 5;
            for (idx, &b) in b"CTGACTGACT".iter().enumerate() {
                let mut nt = BTreeMap::new();
                nt.insert(b, 10);
                sc5.nt.insert(idx as i32, nt);
            }
            lgdel_sc5.insert(i32::MAX - 1, sc5);

            let mut lgdel_sc3 = BTreeMap::new();
            let mut sc3 = super::Sclip::new();
            sc3.variation.cnt = 4;
            lgdel_sc3.insert(i32::MAX - 13, sc3);

            let mut lgdel_non_ins = crate::javamap::JavaHashMap::new();
            let mut lgdel_ins = crate::javamap::JavaHashMap::new();
            let mut lgdel_cov = rustc_hash::FxHashMap::default();
            let mut lgdel_splice = rustc_hash::FxHashSet::default();
            let mut svfdel = Vec::new();
            let mut svrdel = Vec::new();
            super::realignlgdel(
                None,
                &mut lgdel_ref,
                "chr1",
                i32::MAX as u32,
                &lgdel_region,
                &crate::config::Config::default(),
                &mut lgdel_non_ins,
                &mut lgdel_ins,
                &mut lgdel_cov,
                &mut lgdel_sc5,
                &mut lgdel_sc3,
                150,
                &mut svfdel,
                &mut svrdel,
                &mut lgdel_splice,
            );
            let called = lgdel_non_ins
                .get(&(i32::MAX - 15))
                .and_then(|vm| vm.map.get("-14"))
                .unwrap();
            assert_eq!(called.cnt, 9);

            // `find_mm5` and `find_mm3` at `i32::MAX` must not push clamped `i32::MAX`
            // into `scp` when the true soft-clip target position exceeds `i32::MAX`.
            let mut sc5_at_max = BTreeMap::new();
            sc5_at_max.insert(i32::MAX, super::Sclip::new());
            let mm5_at_max = super::find_mm5(&lgdel_ref, "chr1", i32::MAX, b"A", &mut sc5_at_max);
            assert!(mm5_at_max.scp.is_empty());
            assert!(!sc5_at_max[&i32::MAX].used);

            let mut sc3_at_max = BTreeMap::new();
            sc3_at_max.insert(i32::MAX, super::Sclip::new());
            let mm3_at_max = super::find_mm3(&lgdel_ref, "chr1", i32::MAX, b"A", &mut sc3_at_max);
            assert!(mm3_at_max.scp.is_empty());
            assert!(!sc3_at_max[&i32::MAX].used);

            // `realigndel` near `i32::MAX` with `p + dellen + extra.len() > i32::MAX`
            // but `p + dellen + extra.len() - extrains.len() <= i32::MAX`:
            // `p = i32::MAX - 14, vn = "-8^CCCCCC&TGTG"` (`dellen = 8, extrains = "CCCCCC" (6), extra = "CCCCCCTGTG" (10)`)
            // -> `start_pos_ref = (i32::MAX - 14) + 8 + 10 - 6 = i32::MAX - 2`, so `r5` starts at `i32::MAX - 3`
            // and `find_mm5` probes `sc5` at `(i32::MAX - 3) + 1 = i32::MAX - 2` (not clamped `i32::MAX - 6`).
            let mut del_counts = BTreeMap::new();
            let mut vn_map = rustc_hash::FxHashMap::default();
            vn_map.insert("-8^CCCCCC&TGTG".to_string(), 5);
            del_counts.insert(i32::MAX - 14, vn_map);
            let mut rdel_non_ins = crate::javamap::JavaHashMap::new();
            super::get_variation(&mut rdel_non_ins, i32::MAX - 14, "-8^CCCCCC&TGTG").cnt = 5;
            let mut rdel_cov = rustc_hash::FxHashMap::default();
            let mut rdel_sc5 = BTreeMap::new();
            let mut rdel_sc5_clip = super::Sclip::new();
            rdel_sc5_clip.variation.cnt = 3;
            for (idx, &b) in b"GTGTCCCCCCACTCAGTGACT".iter().enumerate() {
                let mut nt = BTreeMap::new();
                nt.insert(b, 10);
                rdel_sc5_clip.nt.insert(idx as i32, nt);
            }
            rdel_sc5.insert(i32::MAX - 2, rdel_sc5_clip);
            let mut rdel_sc3 = BTreeMap::new();
            super::realigndel(
                None,
                &del_counts,
                &lgdel_ref,
                "chr1",
                i32::MAX as u32,
                &crate::config::Config::default(),
                &mut rdel_non_ins,
                &mut rdel_cov,
                &mut rdel_sc5,
                &mut rdel_sc3,
                150,
            );
            assert!(
                rdel_sc5[&(i32::MAX - 2)].used,
                "realigndel must compute start_pos_ref in i64 without clamping before - extrains.len()"
            );
        }
    }

    #[test]
    fn strip_markers_and_mm_stack_buffers_match_with_and_without_markers() {
        // 1. strip_markers returns Borrowed when no '#'/'^' present, Owned when present.
        let clean = super::strip_markers(b"ACGTACGT");
        assert!(matches!(clean, std::borrow::Cow::Borrowed(b"ACGTACGT")));
        let stripped = super::strip_markers(b"AC#GT^AC#GT");
        assert!(matches!(stripped, std::borrow::Cow::Owned(_)));
        assert_eq!(stripped.as_ref(), b"ACGTACGT");

        // 2. is_match forward (dir = 1) and reverse (dir = -1) with and without markers.
        let ref_fwd = b"ACGTACGTAC";
        assert!(super::is_match(ref_fwd, b"ACGTACGTAC", 1, 1));
        assert!(super::is_match(ref_fwd, b"AC#GT^ACGTAC", 1, 1));
        assert!(super::is_match(ref_fwd, b"CATGCATGCA", -1, 1));
        assert!(super::is_match(ref_fwd, b"CA#TG^CATGCA", -1, 1));

        // 3. find_mm5 & find_mm3 stack buffers and marker stripping.
        let mut reference = crate::reference::Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 30,
            });
        let ref_bases = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        reference.region_sequences.push(ref_bases);

        let mut sc5_plain = BTreeMap::new();
        let mut sc5_marked = BTreeMap::new();
        // 3 mismatches at 5' end (reading right-to-left from end of wupseq: 'T', 'G', 'C')
        let mm5_plain = super::find_mm5(&reference, "chr1", 10, b"AAAACGT", &mut sc5_plain);
        let mm5_marked = super::find_mm5(&reference, "chr1", 10, b"AA#AA^CG#T", &mut sc5_marked);
        assert_eq!(mm5_plain, mm5_marked);
        assert_eq!(mm5_plain.mismatches.len(), 3);
        assert_eq!(mm5_plain.mismatches[0].mismatch_sequence, "T");
        assert_eq!(mm5_plain.mismatches[1].mismatch_sequence, "GT");
        assert_eq!(mm5_plain.mismatches[2].mismatch_sequence, "CGT");

        let mut sc3_plain = BTreeMap::new();
        let mut sc3_marked = BTreeMap::new();
        // 2 matching 'A's followed by 4 mismatches ('C', 'G', 'T', 'C') at 3' end (mcnt <= 3)
        let mm3_plain = super::find_mm3(&reference, "chr1", 10, b"AACGTC", &mut sc3_plain);
        let mm3_marked = super::find_mm3(&reference, "chr1", 10, b"A#A^CG#TC", &mut sc3_marked);
        assert_eq!(mm3_plain, mm3_marked);
        assert_eq!(mm3_plain.mismatches.len(), 4);
        assert_eq!(mm3_plain.mismatches[0].mismatch_sequence, "C");
        assert_eq!(mm3_plain.mismatches[1].mismatch_sequence, "CG");
        assert_eq!(mm3_plain.mismatches[2].mismatch_sequence, "CGT");
        assert_eq!(mm3_plain.mismatches[3].mismatch_sequence, "CGTC");
    }

    #[test]
    fn adj_ins_pos_slice_homopolymer_1bp_and_dinuc_rotation() {
        let mut reference = crate::reference::Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 20,
            });
        // 1..=4: GGGG, 5..=9: AAAAA, 10: T, 11..=16: ACACAC, 17..=20: GGGG
        let ref_bases = b"GGGGAAAAATACACACGGGG".to_vec();
        reference.region_sequences.push(ref_bases);

        let rotate = |s: &str, split_idx: usize| -> String {
            if split_idx > 0 {
                format!(
                    "{}{}",
                    s.get(split_idx..).unwrap(),
                    s.get(..split_idx).unwrap()
                )
            } else {
                s.to_string()
            }
        };

        // 1. 1-bp homopolymer shift: bi shifts 9 -> 4 while split_idx == 0.
        let (bi_homo, split_homo) = super::adj_ins_pos_slice(9, b"A", &reference, "chr1");
        let expected_homo = super::adj_ins_pos(9, "A".to_string(), &reference, "chr1");
        assert_eq!(bi_homo, 4);
        assert_eq!(split_homo, 0);
        assert_eq!(bi_homo, expected_homo.base_insert);
        assert_eq!(rotate("A", split_homo), expected_homo.insertion_sequence);

        // 2. Dinucleotide odd-step rotation: "CA" at bi = 15 shifts 15 -> 10 with split_idx == 1 ("AC").
        let (bi_dinuc, split_dinuc) = super::adj_ins_pos_slice(15, b"CA", &reference, "chr1");
        let mut dinuc_in = String::with_capacity(16);
        dinuc_in.push_str("CA");
        let dinuc_ptr = dinuc_in.as_ptr();
        let expected_dinuc = super::adj_ins_pos(15, dinuc_in, &reference, "chr1");
        assert_eq!(bi_dinuc, 10);
        assert_eq!(split_dinuc, 1);
        assert_eq!(bi_dinuc, expected_dinuc.base_insert);
        assert_eq!(rotate("CA", split_dinuc), "AC");
        assert_eq!(rotate("CA", split_dinuc), expected_dinuc.insertion_sequence);
        assert_eq!(expected_dinuc.insertion_sequence.as_ptr(), dinuc_ptr);

        // 3. Dinucleotide full-period wrap: "AC" at bi = 16 shifts 16 -> 10 with split_idx == 0 ("AC").
        let (bi_wrap, split_wrap) = super::adj_ins_pos_slice(16, b"AC", &reference, "chr1");
        let expected_wrap = super::adj_ins_pos(16, "AC".to_string(), &reference, "chr1");
        assert_eq!(bi_wrap, 10);
        assert_eq!(split_wrap, 0);
        assert_eq!(bi_wrap, expected_wrap.base_insert);
        assert_eq!(rotate("AC", split_wrap), expected_wrap.insertion_sequence);
    }

    #[test]
    fn extreme_realign_saturating_helpers_and_u32_max_config_clamps() {
        let mut v = Variation::new();
        v.cnt = i32::MIN;
        v.hicnt = i32::MIN;
        v.locnt = i32::MIN;
        v.fwd = i32::MIN;
        v.rev = i32::MIN;
        let mut tv = Variation::new();
        tv.cnt = 10;
        tv.hicnt = 10;
        tv.locnt = 10;
        tv.fwd = 5;
        tv.rev = 5;
        super::rm_cnt(&mut v, &tv);
        assert_eq!(v.cnt, 0);
        assert_eq!(v.fwd, 0);
        assert_eq!(v.rev, 0);

        let mut vmax = Variation::new();
        vmax.cnt = i32::MAX - 5;
        vmax.hicnt = i32::MAX - 5;
        vmax.fwd = i32::MAX - 5;
        super::add_var_factor(&mut vmax, 0.5);
        assert_eq!(vmax.cnt, i32::MAX);
        assert_eq!(vmax.hicnt, i32::MAX);
        assert_eq!(vmax.fwd, i32::MAX);

        let mut dest = Sclip::new();
        dest.variation.cnt = i32::MAX - 1;
        dest.disc = i32::MAX - 1;
        let mut src = Sclip::new();
        src.variation.cnt = 10;
        src.disc = 10;
        super::merge_sclip(&mut dest, src);
        assert_eq!(dest.variation.cnt, i32::MAX);
        assert_eq!(dest.disc, i32::MAX);

        // Verify realignlgdel and realignlgins clamp u32::MAX thresholds cleanly without wrapping to -1.
        let mut reference = Reference::default();
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 200,
        });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 200]);
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: String::new(),
            amplicon: None,
        };
        let mut cfg = crate::config::Config::default();
        cfg.min_reads = u32::MAX;
        cfg.indel_size = u32::MAX;
        cfg.sv_min_len = u32::MAX;

        let mut sc5 = Sclip::new();
        sc5.variation.cnt = 5;
        sc5.sequence = Some(vec![b'C'; 20]);
        let mut soft5 = std::collections::BTreeMap::new();
        soft5.insert(50, sc5);
        let mut soft3 = std::collections::BTreeMap::new();
        let mut non_ins = crate::javamap::JavaHashMap::new();
        let mut ins_vars = crate::javamap::JavaHashMap::new();
        let mut ref_cov = rustc_hash::FxHashMap::default();
        let mut svfdel = Vec::new();
        let mut svrdel = Vec::new();
        let mut svfdup = Vec::new();
        let mut svrdup = Vec::new();
        let mut splice = rustc_hash::FxHashSet::default();

        super::realignlgdel(
            None,
            &mut reference,
            "chr1",
            200,
            &region,
            &cfg,
            &mut non_ins,
            &mut ins_vars,
            &mut ref_cov,
            &mut soft5,
            &mut soft3,
            100,
            &mut svfdel,
            &mut svrdel,
            &mut splice,
        );
        super::realignlgins(
            None,
            &mut reference,
            "chr1",
            200,
            &region,
            &cfg,
            &mut non_ins,
            &mut ins_vars,
            &mut ref_cov,
            &mut soft5,
            &mut soft3,
            100,
            &mut svfdup,
            &mut svrdup,
            &mut splice,
        );
        assert!(!soft5.get(&50).unwrap().used);
    }

    #[test]
    fn adj_mnt_non_utf8_left_right_and_substr_extreme_isize_do_not_panic() {
        use super::{adjust_mnp, substr, substr_from, Reference, Sclip};
        use crate::types::Pos1;
        use crate::variations::{Variation, VariationData, VariationMap};

        let mut reference = Reference::default();
        reference.loaded_regions.push(LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 100,
        });
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);

        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: String::new(),
            amplicon: None,
        };
        let config = crate::config::Config::default();

        // Case 1: `vn` containing `"\u{FFFD}A"` (`\xEF\xBF\xBD\x41` after stripping `&`)
        // exercises `left_bytes`/`right_bytes` `&` insertion without UTF-8 boundary panic.
        let mut var_data = VariationData::new();
        let vn1 = "\u{FFFD}&A".to_string();
        var_data.mnp.entry(10).or_default().insert(vn1.clone(), 5);
        let mut vm1 = VariationMap::new();
        let mut v1 = Variation::new();
        v1.cnt = 5;
        vm1.map.insert(vn1, v1);
        var_data.non_insertion_variants.insert(10, vm1);

        // Case 2: `vn` whose lossy form matches 2 invalid UTF-8 bytes `[0xff, 0xfe]`
        // (`"\u{FFFD}\u{FFFD}AC"`, 8 UTF-8 bytes vs 4 raw bytes in `seq`) exercises
        // 3' and 5' soft-clip prefix/suffix matching without out-of-bounds slicing.
        let vn2 = "\u{FFFD}&\u{FFFD}AC".to_string();
        var_data.mnp.entry(20).or_default().insert(vn2.clone(), 5);
        let mut vm2 = VariationMap::new();
        let mut v2 = Variation::new();
        v2.cnt = 5;
        vm2.map.insert(vn2, v2);
        var_data.non_insertion_variants.insert(20, vm2);

        let mut sc3 = Sclip::new();
        sc3.variation.cnt = 3;
        sc3.sequence = Some(vec![0xff, 0xfe, b'A', b'C']);
        var_data.soft_clips_3_end.insert(20, sc3);

        let mut sc5 = Sclip::new();
        sc5.variation.cnt = 3;
        // Reversed: [0xff, 0xfe, b'A', b'C'] at the end, length 8 >= mnt.len() (8)
        sc5.sequence = Some(vec![b'C', b'A', 0xfe, 0xff, b'A', b'A', b'A', b'A']);
        var_data.soft_clips_5_end.insert(28, sc5);

        adjust_mnp(&reference, "chr1", 100, &region, &config, &mut var_data);

        // Extreme `isize` bounds on `substr` and `substr_from`.
        let seq = b"ACGTACGT";
        assert_eq!(substr(seq, 1, isize::MAX), b"CGTACGT");
        assert_eq!(substr(seq, 0, isize::MAX), b"ACGTACGT");
        assert_eq!(substr(seq, isize::MIN, 3), b"ACG");
        assert_eq!(substr(seq, -isize::MAX, 3), b"ACG");
        assert_eq!(substr(seq, isize::MAX, 3), b"");
        assert_eq!(substr(seq, 2, isize::MIN), b"");
        assert_eq!(substr(seq, 2, -isize::MAX), b"");
        assert_eq!(substr(b"", isize::MIN, isize::MAX), b"");
        assert_eq!(substr(b"", isize::MAX, isize::MIN), b"");

        assert_eq!(substr_from(seq, isize::MIN), b"ACGTACGT");
        assert_eq!(substr_from(seq, -isize::MAX), b"ACGTACGT");
        assert_eq!(substr_from(seq, isize::MAX), b"");
        assert_eq!(substr_from(seq, 3), b"TACGT");
        assert_eq!(substr_from(seq, -3), b"CGT");
        assert_eq!(substr_from(b"", isize::MIN), b"");
        assert_eq!(substr_from(b"", isize::MAX), b"");
    }

    #[test]
    fn indel_desc_byte_parsers_match_all_12_legacy_regexes() {
        use super::{
            normalize_realignins_tn_into, parse_amp_atgc, parse_atgs_amp_atgs_end,
            parse_begin_minus_number, parse_begin_minus_number_any, parse_begin_plus_atgc,
            parse_caret_atgnc, parse_dup_num_atgc, parse_hash_atgc,
            parse_minus_number_amp_atgcs_end, parse_minus_number_atgnc_sv_atgnc_end,
            parse_up_number_end, strip_caret_digits_end,
        };

        let mut cases: Vec<String> = vec![
            "-12",
            "-12^ATG<dup34>CGT",
            "-12^ATN<inv5>CGN",
            "12^ATN<inv5>CGN",
            "-12^ATG<du>5>CGT",
            "-12^ATG<1234>CGT",
            "-12^ATG<a\nb5>CGT",
            "+ATN",
            "+NAT",
            "+AT<dup5>ATN",
            "+ATN&GC",
            "+AT&GCN",
            "-5^3^ACGT",
            "-5^^ACGT",
            "+AT^1#2",
            "+AT^1&2",
            "+AT^1&",
            "+AT&GC&TT",
            "-12^ATG&GC",
            "-12&GC&TT",
            "prefix-12&GC",
            "prefix+AT&GC",
            "+AT<dup5>GC&TT",
            "+AT<dup5>GC",
            "<dup0>A",
            "<dup12>ACGT",
            "&12&ACGT",
            "#12#ACGT",
            "-",
            "+",
            "^12",
            "^^12",
            "-12^34",
            "-12\n^ATG",
            "",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let tokens = [
            "", "+", "-", "1", "12", "^", "&", "#", "<dup5>", "<inv9>", "A", "AT", "ATN", "N", "GC",
        ];
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        for _ in 0..300 {
            let mut s = String::new();
            for _ in 0..5 {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                let idx = ((state >> 33) as usize) % tokens.len();
                s.push_str(tokens[idx]);
            }
            cases.push(s);
        }

        let mut tn_buf = String::new();
        for vn in &cases {
            let vn = vn.as_str();

            let exp1 = BEGIN_MINUS_NUMBER
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(
                parse_begin_minus_number(vn),
                exp1,
                "BEGIN_MINUS_NUMBER mismatch on {vn:?}"
            );

            let exp2 = UP_NUMBER_END
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(
                parse_up_number_end(vn),
                exp2,
                "UP_NUMBER_END mismatch on {vn:?}"
            );

            let exp3 = MINUS_NUMBER_ATGNC_SV_ATGNC_END
                .captures(vn)
                .and_then(|c| Some((c.get(1)?.as_str(), c.get(2)?.as_str())));
            assert_eq!(
                parse_minus_number_atgnc_sv_atgnc_end(vn),
                exp3,
                "MINUS_NUMBER_ATGNC_SV_ATGNC_END mismatch on {vn:?}"
            );

            let exp4 = BEGIN_MINUS_NUMBER_ANY
                .captures(vn)
                .and_then(|c| Some((c.get(1)?.as_str(), c.get(2)?.as_str())));
            assert_eq!(
                parse_begin_minus_number_any(vn),
                exp4,
                "BEGIN_MINUS_NUMBER_ANY mismatch on {vn:?}"
            );

            let exp5 = CARET_ATGNC
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(
                parse_caret_atgnc(vn),
                exp5,
                "CARET_ATGNC mismatch on {vn:?}"
            );

            let exp6 = MINUS_NUMBER_AMP_ATGCS_END
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(
                parse_minus_number_amp_atgcs_end(vn),
                exp6,
                "MINUS_NUMBER_AMP_ATGCS_END mismatch on {vn:?}"
            );

            let exp7 = BEGIN_PLUS_ATGC
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(
                parse_begin_plus_atgc(vn),
                exp7,
                "BEGIN_PLUS_ATGC mismatch on {vn:?}"
            );

            let exp8 = DUP_NUM_ATGC
                .captures(vn)
                .and_then(|c| Some((c.get(1)?.as_str(), c.get(2)?.as_str())));
            assert_eq!(
                parse_dup_num_atgc(vn),
                exp8,
                "DUP_NUM_ATGC mismatch on {vn:?}"
            );

            let exp9 = AMP_ATGC
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(parse_amp_atgc(vn), exp9, "AMP_ATGC mismatch on {vn:?}");

            let exp10 = HASH_ATGC
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(parse_hash_atgc(vn), exp10, "HASH_ATGC mismatch on {vn:?}");

            let exp11 = ATGS_AMP_ATGS_END
                .captures(vn)
                .and_then(|c| c.get(1).map(|m| m.as_str()));
            assert_eq!(
                parse_atgs_amp_atgs_end(vn),
                exp11,
                "ATGS_AMP_ATGS_END mismatch on {vn:?}"
            );

            let exp12 = RE_CARET_DIGITS_END.replace(vn, "");
            assert_eq!(
                strip_caret_digits_end(vn),
                exp12.as_ref(),
                "RE_CARET_DIGITS_END mismatch on {vn:?}"
            );

            let legacy_tn = vn.strip_prefix('+').unwrap_or(vn).to_string();
            let legacy_tn = legacy_tn.replace(['&', '#'], "");
            let legacy_tn = RE_CARET_DIGITS_END.replace(&legacy_tn, "").to_string();
            let legacy_tn = legacy_tn.replace('^', "");
            normalize_realignins_tn_into(&mut tn_buf, vn);
            assert_eq!(
                tn_buf, legacy_tn,
                "normalize_realignins_tn_into mismatch on {vn:?}"
            );
        }
    }
}
