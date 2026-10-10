//! Local realignment and insertion adjustment — ports VariationRealigner.java.

use crate::reference::Reference;
use crate::types::{Pos0, Pos1};
use crate::variations::{
    get_variation, get_variation_maybe, get_variation_maybe_mut, Sclip, Variation, VariationMap,
};
use regex::Regex;
use std::collections::BTreeMap;
use std::sync::LazyLock;

/// Result of local insertion adjustment.
pub struct BaseInsertion {
    pub base_insert: i32,
    pub insertion_sequence: String,
    pub base_insert2: i32,
}

/// Adjusts the insertion position if it falls in a repeat / homopolymer region.
/// Ports `VariationRealigner.java:2187-2201`.
pub fn adj_ins_pos(
    mut bi: i32,
    mut ins: String,
    reference: &Reference,
    chr: &str,
) -> BaseInsertion {
    let len = ins.len();
    if len == 0 {
        return BaseInsertion {
            base_insert: bi,
            insertion_sequence: ins,
            base_insert2: bi,
        };
    }
    let mut n = 1;
    loop {
        let ins_char = ins.as_bytes()[len - n];
        if let Some(ref_char) = reference.get_base(chr, bi as u32) {
            if ref_char == ins_char {
                n += 1;
                if n > len {
                    n = 1;
                }
                bi -= 1;
                continue;
            }
        }
        break;
    }
    if n > 1 {
        let shift = 1 - n as isize;
        let bytes = ins.as_bytes();
        let len_isize = len as isize;
        let start1 = (len_isize + shift) as usize;
        let mut new_ins = String::new();
        if start1 < len {
            new_ins.push_str(&String::from_utf8_lossy(&bytes[start1..]));
        }
        if start1 > 0 {
            new_ins.push_str(&String::from_utf8_lossy(&bytes[..start1]));
        }
        ins = new_ins;
    }
    BaseInsertion {
        base_insert: bi,
        insertion_sequence: ins,
        base_insert2: bi,
    }
}

/// Check if seq1 matches seq2 with no more than MM mismatches and mismatches ratio < 0.15.
///
/// Ports `VariationRealigner.ismatch`.
pub fn is_match(seq1: &[u8], seq2: &[u8], dir: i32, mm_limit: i32) -> bool {
    let clean_seq2: Vec<u8> = seq2
        .iter()
        .filter(|&&b| b != b'#' && b != b'^')
        .copied()
        .collect();

    let mut mm = 0;
    let limit = seq1.len().min(clean_seq2.len());
    for n in 0..limit {
        let base2 = if dir == -1 {
            clean_seq2[clean_seq2.len() - 1 - n]
        } else {
            clean_seq2[n]
        };
        if seq1[n] != base2 {
            mm += 1;
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
        let ref_pos = position + dir * (n as i32);
        let ref_ch = match reference.get_base(chr, ref_pos as u32) {
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
    let mut bp = 0;
    let mut score = 0;

    for n in 0..indel_size {
        let mut mm = 0;
        let mut i = 0;
        let mut m = std::collections::HashSet::new();

        while i < sequence.len() {
            let ref_pos = start_position + direction * n + direction * (i as i32);
            if ref_pos < 1 || ref_pos > chr_len as i32 {
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
            let lbp = start_position + direction * n - if direction < 0 { direction } else { 0 };
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
    let dir_ext = if dir == -1 { 1 } else { 0 };
    let mut score = 0;
    let mut bi = 0;
    let mut ins = String::new();
    let mut bi2 = 0;

    for n in 6..seq.len() {
        if position + 6 >= chr_len as i32 {
            break;
        }

        let mut mm = 0;
        let mut i = 0;
        let mut m = std::collections::HashSet::new();

        while i + n < seq.len() {
            let ref_pos = position + dir * (i as i32) - dir_ext;
            if ref_pos < 1 || ref_pos > chr_len as i32 {
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
                let ref_pos1 = position + (ept as i32) * dir - dir_ext;
                let ref_pos2 = position + ((ept as i32) + 1) * dir - dir_ext;

                let ref_ch1 = reference.get_base(chr, ref_pos1 as u32);
                let ref_ch2 = reference.get_base(chr, ref_pos2 as u32);

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
                let mut ins_str = String::from_utf8_lossy(&ins_vec).into_owned();
                if !extra.is_empty() {
                    let insert_idx = ins_str.len() - extra.len();
                    ins_str.insert(insert_idx, '&');
                }

                if mm == 0 && (i + n) == seq.len() {
                    let mut bi_val = position - 1 - (extra.len() as i32);
                    let mut ins_val = ins_str;
                    let mut bi2_val = position - 1;
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
                    bi = position - 1 - (extra.len() as i32);
                    ins = ins_str;
                    bi2 = position - 1;
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
                            base_insert: position + s,
                            insertion_sequence: ins_str,
                            base_insert2: position + s + (extra.len() as i32),
                        };
                    } else if (i as i32) - mm > score {
                        bi = position + s;
                        ins = ins_str;
                        bi2 = position + s + (extra.len() as i32);
                        score = (i as i32) - mm;
                    }
                } else {
                    while s >= -(n as i32) {
                        let ins_ch = crate::util::char_at(&insert, s as isize);
                        let ref_pos = position + s;
                        let ref_ch = reference.get_base(chr, ref_pos as u32);
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
                        let mut bi_val = position + s;
                        let mut ins_val = ins_str;
                        let mut bi2_val = position + s + (extra.len() as i32);
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
                        bi = position + s;
                        ins = ins_str;
                        bi2 = position + s + (extra.len() as i32);
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

#[derive(Debug, Clone, PartialEq, Eq)]
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
    let seq: Vec<u8> = wupseq
        .iter()
        .filter(|&&b| b != b'#' && b != b'^')
        .copied()
        .collect();

    let longmm = 3;
    let mut mismatches = Vec::new();
    let mut n = 0;
    let mut mn = 0;
    let mut mcnt = 0;
    let mut str_bytes = Vec::new();
    let mut sc5p = Vec::new();

    loop {
        let seq_idx = (seq.len() as i32) - 1 - n;
        if seq_idx < 0 || seq_idx >= seq.len() as i32 {
            break;
        }
        let ref_pos = position - n;
        let ref_ch = reference.get_base(chr, ref_pos as u32);

        if let Some(c_ref) = ref_ch {
            let c_seq = seq[seq_idx as usize];
            if c_seq != c_ref && mcnt < longmm {
                str_bytes.insert(0, c_seq);
                let mismatch_str = String::from_utf8_lossy(&str_bytes).into_owned();
                mismatches.push(Mismatch {
                    mismatch_sequence: mismatch_str,
                    mismatch_position: position - n,
                    end: 5,
                });
                n += 1;
                mcnt += 1;
                continue;
            }
        }
        break;
    }
    sc5p.push(position + 1);

    let misp = 0;
    let misnt: Option<u8> = None;
    if str_bytes.len() == 1 {
        loop {
            let seq_idx = (seq.len() as i32) - 1 - n;
            if seq_idx < 0 || seq_idx >= seq.len() as i32 {
                break;
            }
            let ref_pos = position - n;
            let ref_ch = reference.get_base(chr, ref_pos as u32);
            if let Some(c_ref) = ref_ch {
                if seq[seq_idx as usize] == c_ref {
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
            sc5p.push(position - n);
            let target_pos = position - n;
            if let Some(sclip) = soft_clips_5_end.get_mut(&target_pos) {
                sclip.used = true;
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
    let seq: Vec<u8> = sanpseq
        .iter()
        .filter(|&&b| b != b'#' && b != b'^')
        .copied()
        .collect();

    let longmm = 3;
    let mut mismatches = Vec::new();
    let mut n = 0;
    let mut mn = 0;
    let mut mcnt = 0;
    let mut sc3p = Vec::new();
    let mut str_bytes = Vec::new();

    while n < seq.len() {
        let ref_pos = p + (n as i32);
        let ref_ch = reference.get_base(chr, ref_pos as u32);
        if let Some(c_ref) = ref_ch {
            if c_ref == seq[n] {
                n += 1;
                continue;
            }
        }
        break;
    }
    sc3p.push(p + (n as i32));
    let tbp = p + (n as i32);

    while mcnt <= longmm && n < seq.len() {
        let ref_pos = p + (n as i32);
        let ref_ch = reference.get_base(chr, ref_pos as u32);
        if let Some(c_ref) = ref_ch {
            if c_ref != seq[n] {
                str_bytes.push(seq[n]);
                let mismatch_str = String::from_utf8_lossy(&str_bytes).into_owned();
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
    if str_bytes.len() == 1 {
        while n < seq.len() {
            let ref_pos = p + (n as i32);
            let ref_ch = reference.get_base(chr, ref_pos as u32);
            if let Some(c_ref) = ref_ch {
                if seq[n] == c_ref {
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
            while n + n2 as usize + 1 < seq.len() {
                let ref_pos = p + (n as i32) + 1 + n2;
                let ref_ch = reference.get_base(chr, ref_pos as u32);
                if let Some(c_ref) = ref_ch {
                    if seq[n + n2 as usize + 1] == c_ref {
                        n2 += 1;
                        continue;
                    }
                }
                break;
            }

            if n2 > 2 && n + (n2 as usize) + 1 < seq.len() {
                sc3p.push(p + (n as i32) + n2);
                misp = p + (n as i32);
                misnt = Some(seq[n]);
                let target_pos = p + (n as i32) + n2;
                if let Some(sclip) = soft_clips_3_end.get_mut(&target_pos) {
                    sclip.used = true;
                }
                mn += n2;
            } else {
                sc3p.push(p + (n as i32));
                let target_pos = p + (n as i32);
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
    changes: &BTreeMap<i32, std::collections::HashMap<String, i32>>,
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
pub fn fill_adaptor_maps(
    adaptor_seqs: &[String],
) -> (BTreeMap<Vec<u8>, i32>, BTreeMap<Vec<u8>, i32>) {
    let mut adaptor_forward = BTreeMap::new();
    let mut adaptor_reverse = BTreeMap::new();
    let adseed = 6;
    for sequence in adaptor_seqs {
        let seq_bytes = sequence.as_bytes();
        for i in 0..=6 {
            if i + adseed <= seq_bytes.len() {
                let forward_seed = seq_bytes[i..i + adseed].to_vec();
                let mut comp_reverse_seed = forward_seed.clone();
                crate::util::reverse_complement(&mut comp_reverse_seed);
                adaptor_forward.insert(forward_seed, (i + 1) as i32);
                adaptor_reverse.insert(comp_reverse_seed, (i + 1) as i32);
            }
        }
    }
    (adaptor_forward, adaptor_reverse)
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
    let dlen = end - start;
    let dlenqr = format!("{}D", dlen);

    for bam in bams {
        let bam_path = std::path::Path::new(bam);
        let start_pos0 = Pos0((start - 1).max(0) as u32);
        let end_pos0 = Pos0(end as u32);
        let records = match crate::bam::query_region(bam_path, chr, start_pos0, end_pos0, config) {
            Ok(recs) => recs,
            Err(e) => {
                log::warn!("no_passing_reads: failed to query bam '{}': {}", bam, e);
                continue;
            }
        };

        for record in records {
            let cigar_str = crate::bam::get_cigar_string(&record.cigar);
            if cigar_str.contains(&dlenqr) {
                continue;
            }
            let read_start = match record.alignment_start {
                Some(Pos1(p)) => p as i32,
                None => continue,
            };
            let read_len = crate::bam::get_aligned_length_md(&record.cigar);
            let read_end = read_start + read_len;

            if read_end > end + 2 && read_start < start - 2 {
                cnt += 1;
            }
            if read_start < start - 2 && read_end > start && read_end < end {
                midcnt += 1;
            }
        }
    }

    if config.verbose {
        eprintln!(
            "    Passing Read CNT: {} {} {} {} {}",
            cnt, chr, start, end, midcnt
        );
    }

    cnt <= 0 && midcnt + 1 > 0
}

static BEGIN_MINUS_NUMBER: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^-(\d+)").unwrap());
static UP_NUMBER_END: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\^(\d+)$").unwrap());
static MINUS_NUMBER_ATGNC_SV_ATGNC_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^-?\d+\^([ATGNC]+)<...\d+>([ATGNC]+)$").unwrap());
static BEGIN_MINUS_NUMBER_ANY: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^-(\d+)(.*)").unwrap());
static CARET_ATGNC: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\^([ATGNC]+)").unwrap());
static MINUS_NUMBER_AMP_ATGCS_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(-\d+)&[ATGC]+$").unwrap());

static BEGIN_PLUS_ATGC: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\+([ATGC]+)").unwrap());
static DUP_NUM_ATGC: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"<dup(\d+)>([ATGC]+)$").unwrap());
static AMP_ATGC: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"&([ATGC]+)").unwrap());
static HASH_ATGC: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"#([ATGC]+)").unwrap());
static ATGS_AMP_ATGS_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\+[ATGC]+)&[ATGC]+$").unwrap());
static RE_LEADING_PLUS: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\+").unwrap());
static RE_CARET_DIGITS_END: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\^\d+$").unwrap());

/// Realign deletions.
///
/// Ports `VariationRealigner.realigndel`.
pub fn realigndel(
    bams: Option<&[String]>,
    position_to_deletion_count: &BTreeMap<i32, std::collections::HashMap<String, i32>>,
    reference: &Reference,
    chr: &str,
    chr_len: u32,
    config: &crate::config::Config,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut BTreeMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_maps(&config.adaptor);
    let tmp = fill_and_sort_tmp(position_to_deletion_count);

    for tpl in &tmp {
        let p = tpl.position;
        let vn = &tpl.description_string;
        let dcnt = tpl.count;

        if config.verbose {
            let cov = ref_coverage.get(&p).copied().unwrap_or(0);
            eprintln!("  Realigndel for: {} {} {} cov: {}", p, vn, dcnt, cov);
        }

        let mut dellen = 0;
        if let Some(caps) = BEGIN_MINUS_NUMBER.captures(vn) {
            if let Some(m) = caps.get(1) {
                dellen = m.as_str().parse::<i32>().unwrap_or(0);
            }
        }
        if let Some(caps) = UP_NUMBER_END.captures(vn) {
            if let Some(m) = caps.get(1) {
                dellen += m.as_str().parse::<i32>().unwrap_or(0);
            }
        }

        let mut extrains = String::new();
        let mut extra = String::new();
        let mut inv5 = String::new();
        let mut inv3 = String::new();

        if let Some(caps) = MINUS_NUMBER_ATGNC_SV_ATGNC_END.captures(vn) {
            inv5 = caps
                .get(1)
                .map(|m| m.as_str().to_string())
                .unwrap_or_default();
            inv3 = caps
                .get(2)
                .map(|m| m.as_str().to_string())
                .unwrap_or_default();
        } else if let Some(caps) = BEGIN_MINUS_NUMBER_ANY.captures(vn) {
            if let Some(m) = caps.get(2) {
                extra = m.as_str().replace(&['^', '&', '#'][..], "");
            }
            if let Some(caps2) = CARET_ATGNC.captures(vn) {
                extrains = caps2
                    .get(1)
                    .map(|m| m.as_str().to_string())
                    .unwrap_or_default();
            }
        }

        let wustart = if p - 200 > 1 { p - 200 } else { 1 };
        let mut wupseq = crate::util::join_ref(reference, chr, wustart, p - 1);
        wupseq.push_str(&extra);
        if !inv3.is_empty() {
            wupseq = inv3.clone();
        }

        let sanend = if p + 200 > chr_len as i32 {
            chr_len as i32
        } else {
            p + 200
        };
        let start_pos_ref = p + dellen + (extra.len() as i32) - (extrains.len() as i32);
        let mut sanpseq = extra.clone();
        sanpseq.push_str(&crate::util::join_ref(
            reference,
            chr,
            start_pos_ref,
            sanend,
        ));
        if !inv5.is_empty() {
            sanpseq = inv5.clone();
        }

        let r3 = find_mm3(reference, chr, p, sanpseq.as_bytes(), soft_clips_3_end);
        let r5 = find_mm5(
            reference,
            chr,
            p + dellen + (extra.len() as i32) - (extrains.len() as i32) - 1,
            wupseq.as_bytes(),
            soft_clips_5_end,
        );

        let mm3 = &r3.mismatches;
        let sc3p = &r3.scp;
        let nm3 = r3.nm;
        let misp3 = r3.misp;
        let misnt3 = &r3.misnt;

        let mm5 = &r5.mismatches;
        let sc5p = &r5.scp;
        let nm5 = r5.nm;
        let misp5 = r5.misp;
        let misnt5 = &r5.misnt;

        if config.verbose {
            eprintln!(
                "  Mismatches: misp3: {}-{} misp5: {}-{} sclip3: {:?} sclip5: {:?}",
                misp3, misnt3, misp5, misnt5, sc3p, sc5p
            );
        }

        let mut mmm = mm3.clone();
        mmm.extend(mm5.clone());

        for mismatch in mmm {
            let mut mm = mismatch.mismatch_sequence.clone();
            let mp = mismatch.mismatch_position;
            let me = mismatch.end;

            if mm.len() > 1 {
                mm = format!("{}&{}", &mm[0..1], &mm[1..]);
            }

            let mut should_adjust = false;
            if let Some(entry) = non_insertion_variants.get(&mp) {
                if let Some(tv) = entry.map.get(&mm) {
                    if tv.cnt > 0 {
                        let mean_qual = tv.qmean / (tv.cnt as f64);
                        let mean_pos = tv.pmean / (tv.cnt as f64);
                        let threshold_pos = if me == 3 { nm3 + 4 } else { nm5 + 4 };
                        if mean_qual >= config.base_qual
                            && mean_pos <= (threshold_pos as f64)
                            && (tv.cnt < dcnt + dellen && (tv.cnt as f64) / (dcnt as f64) < 8.0)
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
                eprintln!(
                    "  Realigndel Adj: {} {} {} {} {} {} {} {} cov: {}",
                    mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean, cov
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
                *ref_coverage.entry(p).or_default() += (tv.cnt as f64 * f) as i32;

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
                        lref.cnt -= tv.cnt;
                        lref.hicnt -= tv.hicnt;
                        lref.locnt -= tv.locnt;
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
                eprintln!(
                    "  Realigndel AdjA: {} {} {} {} {} {} {} {} cov: {}",
                    mm, mp, me, nm3, nm5, p, final_vref.cnt, final_vref.qmean, cov
                );
            }
        }

        if misp3 != 0 && mm3.len() == 1 {
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

        if misp5 != 0 && mm5.len() == 1 {
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
                    let seq = crate::variations::util::find_conseq(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if dcnt <= 2 && (tv.variation.cnt as f64) / (dcnt as f64) > 5.0 {
                        continue;
                    }
                    if config.verbose {
                        let cov = ref_coverage.get(&p).copied().unwrap_or(0);
                        let wupseq_rev: String = wupseq.chars().rev().collect();
                        eprintln!(
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
                            *ref_coverage.entry(p).or_default() += tv.variation.cnt;
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
                            eprintln!(
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
                    let seq = crate::variations::util::find_conseq(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if dcnt <= 2 && (tv.variation.cnt as f64) / (dcnt as f64) > 5.0 {
                        continue;
                    }

                    let slice_start = (sc3pp - p).max(0) as usize;
                    let sanpseq_slice = if slice_start < sanpseq.len() {
                        &sanpseq[slice_start..]
                    } else {
                        ""
                    };

                    if config.verbose {
                        eprintln!(
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
                            eprintln!(
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
                            *ref_coverage.entry(p).or_default() += tv.variation.cnt;
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
                            lref.cnt -= tv.variation.cnt;
                            lref.hicnt -= tv.variation.hicnt;
                            lref.locnt -= tv.variation.locnt;
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

        let pe = p + dellen + (extra.len() as i32) - (extrains.len() as i32);

        let mut h_owned = None;
        if let Some(ref_base) = reference.get_base(chr, p as u32) {
            if let Some(h) = get_variation_maybe(non_insertion_variants, p, Some(ref_base)) {
                h_owned = Some(h.clone());
            }
        }

        if let Some(bams_slice) = bams {
            if !bams_slice.is_empty() && pe - p >= 5 && pe - p < max_read_length - 10 {
                if let Some(ref h) = h_owned {
                    if h.cnt != 0 {
                        let vref = get_variation(non_insertion_variants, p, vn);
                        let threshold_vref_cnt = 2.0
                            * (h.cnt as f64)
                            * (1.0 - ((pe - p) as f64) / (max_read_length as f64));
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
                                    let dir_t = h_mut.get_dir(true);
                                    let dir_f = h_mut.get_dir(false);
                                    h_mut.add_dir(true, -dir_t);
                                    h_mut.add_dir(false, -dir_f);
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

            if let Some(caps) = MINUS_NUMBER_AMP_ATGCS_END.captures(vn) {
                if let Some(tn_match) = caps.get(1) {
                    let tn = tn_match.as_str().to_string();

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
                        if let Some(tref) = entry.map.get(&tn) {
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

                        let tref = get_variation(non_insertion_variants, p, &tn);
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
}

/// Realign insertions.
///
/// Ports `VariationRealigner.realignins`.
pub fn realignins(
    _bams: Option<&[String]>,
    position_to_insertion_count: &BTreeMap<i32, std::collections::HashMap<String, i32>>,
    reference: &Reference,
    chr: &str,
    chr_len: u32,
    config: &crate::config::Config,
    insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut BTreeMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
) -> String {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_maps(&config.adaptor);
    let tmp = fill_and_sort_tmp(position_to_insertion_count);
    let mut new_ins = String::new();

    for tpl in &tmp {
        let p = tpl.position;
        let vn = &tpl.description_string;
        let insertion_count = tpl.count;

        if config.verbose {
            eprintln!("  Realign Ins: {} {} {}", p, vn, insertion_count);
        }

        let insert = if let Some(caps) = BEGIN_PLUS_ATGC.captures(vn) {
            caps.get(1)
                .map_or(String::new(), |m| m.as_str().to_string())
        } else {
            continue;
        };

        let mut ins3 = String::new();
        let mut inslen = insert.len() as i32;
        if let Some(caps) = DUP_NUM_ATGC.captures(vn) {
            ins3 = caps
                .get(2)
                .map_or(String::new(), |m| m.as_str().to_string());
            let dup_len = caps
                .get(1)
                .map_or(0, |m| m.as_str().parse::<i32>().unwrap_or(0));
            inslen += dup_len + ins3.len() as i32;
        }

        let mut extra = String::new();
        if let Some(caps) = AMP_ATGC.captures(vn) {
            extra = caps
                .get(1)
                .map_or(String::new(), |m| m.as_str().to_string());
        }

        let mut compm = String::new();
        if let Some(caps) = HASH_ATGC.captures(vn) {
            compm = caps
                .get(1)
                .map_or(String::new(), |m| m.as_str().to_string());
        }

        let mut newdel = 0;
        if let Some(caps) = UP_NUMBER_END.captures(vn) {
            newdel = caps
                .get(1)
                .map_or(0, |m| m.as_str().parse::<i32>().unwrap_or(0));
        }

        let tn = RE_LEADING_PLUS.replace(vn, "").to_string();
        let tn = tn.replace(['&', '#'], "");
        let tn = RE_CARET_DIGITS_END.replace(&tn, "").to_string();
        let tn = tn.replace('^', "");

        let wustart = if p - 150 > 1 { p - 150 } else { 1 };
        let mut wupseq = crate::util::join_ref(reference, chr, wustart, p);
        wupseq.push_str(&tn);

        let sanend = {
            let mut end_pos = p + (vn.len() as i32) + 100;
            if end_pos > chr_len as i32 {
                end_pos = chr_len as i32;
            }
            end_pos
        };

        let mut sanpseq = String::new();
        let find_mm3_res;

        if !ins3.is_empty() {
            let sv_flanking = 50;
            let p3 = p + inslen - (ins3.len() as i32) + sv_flanking;
            if ins3.len() > sv_flanking as usize {
                let start_idx = ins3.len() - sv_flanking as usize;
                sanpseq = ins3[start_idx..].to_string();
            }
            sanpseq.push_str(&crate::util::join_ref(reference, chr, p + 1, p + 101));
            find_mm3_res = find_mm3(reference, chr, p3 + 1, sanpseq.as_bytes(), soft_clips_3_end);
        } else {
            let start_pos_ref = p + (extra.len() as i32) + 1 + (compm.len() as i32) + newdel;
            sanpseq = tn.clone();
            sanpseq.push_str(&crate::util::join_ref(
                reference,
                chr,
                start_pos_ref,
                sanend,
            ));
            find_mm3_res = find_mm3(reference, chr, p + 1, sanpseq.as_bytes(), soft_clips_3_end);
        }

        let p5_pos = p + (extra.len() as i32) + (compm.len() as i32) + newdel;
        let find_mm5_res = find_mm5(reference, chr, p5_pos, wupseq.as_bytes(), soft_clips_5_end);

        let mm3 = &find_mm3_res.mismatches;
        let sc3p = &find_mm3_res.scp;
        let nm3 = find_mm3_res.nm;
        let misp3 = find_mm3_res.misp;
        let misnt3 = &find_mm3_res.misnt;

        let mm5 = &find_mm5_res.mismatches;
        let sc5p = &find_mm5_res.scp;
        let nm5 = find_mm5_res.nm;
        let misp5 = find_mm5_res.misp;
        let misnt5 = &find_mm5_res.misnt;

        let mut mmm = mm3.clone();
        mmm.extend(mm5.clone());

        for mismatch in mmm {
            let mut mismatch_bases = mismatch.mismatch_sequence.clone();
            let mismatch_position = mismatch.mismatch_position;
            let mismatch_end = mismatch.end;

            if mismatch_bases.len() > 1 {
                mismatch_bases = format!("{}&{}", &mismatch_bases[0..1], &mismatch_bases[1..]);
            }

            let mut should_adjust = false;
            if let Some(entry) = non_insertion_variants.get(&mismatch_position) {
                if let Some(variation) = entry.map.get(&mismatch_bases) {
                    if variation.cnt > 0 {
                        let mean_qual = variation.qmean / (variation.cnt as f64);
                        let mean_pos = variation.pmean / (variation.cnt as f64);
                        let threshold_pos = if mismatch_end == 3 { nm3 + 4 } else { nm5 + 4 };
                        if mean_qual >= config.base_qual
                            && mean_pos <= (threshold_pos as f64)
                            && (variation.cnt < insertion_count + (insert.len() as i32)
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
                eprintln!(
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
                *ref_coverage.entry(p).or_default() += variation.cnt;
            }

            let mut lref_val = None;
            if mismatch_position > p && mismatch_end == 3 {
                if let Some(ref_base) = reference.get_base(chr, p as u32) {
                    let ref_base_str = (ref_base as char).to_string();
                    if let Some(entry) = non_insertion_variants.get_mut(&p) {
                        if let Some(lref) = entry.map.get_mut(&ref_base_str) {
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

        if misp3 != 0 && mm3.len() == 1 {
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

        if misp5 != 0 && mm5.len() == 1 {
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
                    eprintln!(
                        "    55: {} {} VN: '{}'  5' seq: ^{}^\n",
                        p, sc5pp, vn, wupseq
                    );
                }
                if !tv.used {
                    let seq = crate::variations::util::find_conseq(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if config.verbose {
                        eprintln!(
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
                            eprintln!(
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
                            *ref_coverage.entry(p).or_default() += tv.variation.cnt;
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
                    eprintln!(
                        "    33: {} {} VN: '{}'  3' seq: ^{}^\n",
                        p, sc3pp, vn, sanpseq
                    );
                }
                if !tv.used {
                    let seq = crate::variations::util::find_conseq(
                        tv,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if config.verbose {
                        eprintln!(
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
                        substr_from(sanpseq.as_bytes(), (sc3pp - p - 1) as isize)
                    };

                    if !seq.is_empty() && is_match(&seq, mseq, 1, 3) {
                        if config.verbose {
                            eprintln!(
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
                            *ref_coverage.entry(p).or_default() += tv.variation.cnt;
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
                                let ref_base_str = (ref_base as char).to_string();
                                if let Some(entry) = non_insertion_variants.get_mut(&p) {
                                    if let Some(lref) = entry.map.get_mut(&ref_base_str) {
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
                            && sc3pp >= p + 1 + (insert.len() as i32)
                        {
                            let mut flag = 0;
                            let offset = (sc3pp - p - 1) as usize % insert.len();
                            let mut tvn = vn.clone();
                            let insert_bytes = insert.as_bytes();

                            for seqi in 0..seq.len() {
                                if seqi + offset < insert_bytes.len() {
                                    let seq_char = seq[seqi] as char;
                                    let ins_char = insert_bytes[seqi + offset] as char;
                                    if seq_char != ins_char {
                                        flag += 1;
                                        let shift = seqi + offset + 1;
                                        if shift < tvn.len() {
                                            tvn = format!(
                                                "{}{}{}",
                                                &tvn[..shift],
                                                seq_char,
                                                &tvn[shift + 1..]
                                            );
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
            // ports VariationRealigner.java:921: the bound stays a double (no truncation).
            if first3 > first5 + 3 && ((first3 - first5) as f64) < (max_read_length as f64) * 0.75 {
                if let Some(ref_base) = reference.get_base(chr, p as u32) {
                    let ref_base_str = (ref_base as char).to_string();
                    if let Some(entry) = non_insertion_variants.get_mut(&p) {
                        if let Some(ref_var) = entry.map.get_mut(&ref_base_str) {
                            crate::variations::util::adj_ref_factor(
                                ref_var,
                                ((first3 - first5 - 1) as f64) / (max_read_length as f64),
                            );
                        }
                    }
                }
                let vref_mut = get_variation(insertion_variants, p, vn);
                crate::variations::util::adj_ref_factor(
                    vref_mut,
                    -((first3 - first5 - 1) as f64) / (max_read_length as f64),
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

        if let Some(caps) = ATGS_AMP_ATGS_END.captures(vn) {
            let tn = caps
                .get(1)
                .map_or(String::new(), |m| m.as_str().to_string());
            let mut tref_cnt = 0;
            let mut tref_exists = false;
            if let Some(entry) = insertion_variants.get(&p) {
                if let Some(tref) = entry.map.get(&tn) {
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
                    let ref_base_str = (ref_base as char).to_string();
                    if let Some(entry) = non_insertion_variants.get_mut(&p) {
                        if let Some(lref) = entry.map.get_mut(&ref_base_str) {
                            ref_base_var = Some(lref);
                        }
                    }
                }

                let tref = get_variation(insertion_variants, p, &tn);
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
    ref_coverage: &mut BTreeMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
    svfdel: &mut Vec<Sclip>,
    svrdel: &mut Vec<Sclip>,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_maps(&config.adaptor);
    let extension = 5000;
    let longmm = 3;

    // Work on 5' clipped reads

    let mut tmp5 = Vec::new();
    for (&p, sc) in soft_clips_5_end.iter() {
        if p < region.start.0 as i32 - extension || p > region.end.0 as i32 + extension {
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
        let mut sc5v = match soft_clips_5_end.get(&p).cloned() {
            Some(sc) => sc,
            None => continue,
        };
        let cnt = sc5v.variation.cnt;
        if cnt < config.min_reads as i32 {
            break;
        }
        if sc5v.used {
            continue;
        }
        let seq =
            crate::variations::util::find_conseq(&mut sc5v, 5, &adaptor_forward, &adaptor_reverse);
        if let Some(entry) = soft_clips_5_end.get_mut(&p) {
            entry.used = sc5v.used;
            entry.sequence = sc5v.sequence.clone();
        }
        if seq.is_empty() {
            continue;
        }
        if seq.len() < 7 {
            continue;
        }

        if config.verbose {
            eprintln!(
                "  Working Realignlgdel: 5' {} '{}' {}",
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }

        let mut bp = find_bp(
            &seq,
            p - 5,
            reference,
            -1,
            chr,
            config.indel_size as i32,
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
            if !(bp != 0 && p - bp > 15 && p - bp < 150000) {
                continue;
            }
            bp += 1;

            let mut lists = [&mut *svfdel, &mut *svrdel];
            let sv_mark = crate::sv::mark_sv(bp, p, &mut lists, max_read_length, config.verbose);
            svcov = sv_mark.0;
            let clusters = sv_mark.1;
            let pairs = sv_mark.2;

            if svcov == 0 && cnt <= config.min_reads as i32 {
                continue;
            }

            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, bp);
            sv.r#type = "DEL".to_string();
            sv.pairs += pairs;
            sv.splits += cnt;
            sv.clusters += clusters;

            if bp < region.start.0 as i32 {
                let tts = bp - max_read_length;
                let mut tte = bp + max_read_length;
                if bp + max_read_length >= region.start.0 as i32 {
                    tte = region.start.0 as i32 - 1;
                }
                if let Err(e) = ensure_reference_loaded(reference, chr, tts, tte, config) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }
        }

        let mut dellen = p - bp;
        let mut en = 0;
        let gt;

        if extra_str.is_empty() {
            while en < seq.len() {
                let ref_pos = bp - en as i32 - 1;
                let ref_ch = reference.get_base(chr, ref_pos as u32);
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
                bp -= extra.len() as i32;
            } else {
                gt = format!("-{}", dellen).into_bytes();
            }
        } else {
            dellen -= extra_str.len() as i32;
            if dellen == 0 {
                gt = format!("-^{}", String::from_utf8_lossy(&extra_str)).into_bytes();
            } else {
                gt = format!("-{}&{}", dellen, String::from_utf8_lossy(&extra_str)).into_bytes();
            }
        }

        if config.verbose {
            eprintln!(
                "  Found Realignlgdel: {} {} 5' {} '{}' {}",
                bp,
                String::from_utf8_lossy(&gt),
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }

        let mut n = 0;
        if extra.is_empty() && extra_str.is_empty() {
            loop {
                let base_bp = reference.get_base(chr, (bp + n) as u32);
                let base_bp_del = reference.get_base(chr, (bp + dellen + n) as u32);
                if base_bp.is_some() && base_bp_del.is_some() && base_bp == base_bp_del {
                    n += 1;
                } else {
                    break;
                }
            }
        }

        let mut sc3p = bp + n;
        let mut str_seq = Vec::new();
        let mut mcnt = 0;
        while mcnt <= longmm {
            let base_bp = reference.get_base(chr, (bp + n) as u32);
            let base_bp_del = reference.get_base(chr, (bp + dellen + n) as u32);
            if let (Some(b_bp), Some(b_bp_del)) = (base_bp, base_bp_del) {
                if b_bp != b_bp_del {
                    str_seq.push(b_bp_del);
                    n += 1;
                    mcnt += 1;
                    continue;
                }
            }
            break;
        }

        if str_seq.len() == 1 {
            let mut nm = 0;
            loop {
                let base_bp = reference.get_base(chr, (bp + n) as u32);
                let base_bp_del = reference.get_base(chr, (bp + dellen + n) as u32);
                if base_bp.is_some() && base_bp_del.is_some() && base_bp == base_bp_del {
                    n += 1;
                    if n != 0 {
                        nm += 1;
                    }
                } else {
                    break;
                }
            }
            if nm >= 3 && !soft_clips_3_end.contains_key(&sc3p) {
                sc3p = bp + n;
            }
        }

        let mut should_continue = false;
        if let Some(entry) = non_insertion_variants.get(&bp) {
            if entry.sv.is_some()
                && !soft_clips_3_end.contains_key(&sc3p)
                && svcov == 0
                && cnt <= config.min_reads as i32
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
            crate::variations::util::adj_cnt(tv, &sc5v.variation, None);
        }
        sc5v.used = bp != 0;

        if let Some(sc) = soft_clips_5_end.get_mut(&p) {
            sc.used = bp != 0;
            sc.sequence = sc5v.sequence;
        }

        if !ref_coverage.contains_key(&bp) {
            if let Some(&cov_p) = ref_coverage.get(&p) {
                ref_coverage.insert(bp, cov_p);
            }
        }

        if dellen < config.indel_size as i32 {
            for tp in bp..bp + dellen {
                inc_cnt(ref_coverage, tp, sc5v.variation.cnt);
            }
        }

        if let Some(mut sclip) = soft_clips_3_end.get(&sc3p).cloned() {
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
                        ref_var.cnt -= sclip.variation.cnt;
                        ref_var.hicnt -= sclip.variation.hicnt;
                        ref_var.locnt -= sclip.variation.locnt;
                        ref_var.pmean -= sclip.variation.pmean;
                        ref_var.qmean -= sclip.variation.qmean;
                        ref_var.mqmean -= sclip.variation.mqmean;
                        ref_var.nm -= sclip.variation.nm;
                        ref_var.sub_dir(true, sclip.variation.get_dir(true));
                        ref_var.sub_dir(false, sclip.variation.get_dir(false));
                        crate::variations::util::correct_cnt(ref_var);
                    }
                }

                if sc3p == bp && dellen < config.indel_size as i32 {
                    for tp in bp..bp + dellen {
                        inc_cnt(ref_coverage, tp, sclip.variation.cnt);
                    }
                }

                for ip in bp + 1..sc3p {
                    let ref_ch_opt = reference.get_base(chr, (dellen + ip) as u32);
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
                            let ref_ch_str = (ref_ch as char).to_string();
                            if let Some(vm) = non_insertion_variants.get_mut(&ip) {
                                vm.map.shift_remove(&ref_ch_str);
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
                if let Some(sc) = soft_clips_3_end.get_mut(&sc3p) {
                    sc.used = bp != 0;
                }
            }
        }

        let mut dels5 = BTreeMap::new();
        let mut inner_map = std::collections::HashMap::new();
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
                        sv.splits += tv_cnt - dels5_cnt;
                    }
                }
            }
        }

        let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
        if svcov > tv_cnt {
            let tv = crate::variations::get_variation(non_insertion_variants, bp, &gt_str);
            add_var_factor(tv, (svcov - tv_cnt) as f64 / tv_cnt as f64);
        }

        if config.verbose {
            let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
            eprintln!(
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
        if p < region.start.0 as i32 - extension || p > region.end.0 as i32 + extension {
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
        let mut sc3v = match soft_clips_3_end.get(&p).cloned() {
            Some(sc) => sc,
            None => continue,
        };
        let cnt = sc3v.variation.cnt;
        if cnt < config.min_reads as i32 {
            break;
        }
        if sc3v.used {
            continue;
        }
        let seq =
            crate::variations::util::find_conseq(&mut sc3v, 3, &adaptor_forward, &adaptor_reverse);
        if let Some(entry) = soft_clips_3_end.get_mut(&p) {
            entry.used = sc3v.used;
            entry.sequence = sc3v.sequence.clone();
        }
        if seq.is_empty() {
            continue;
        }
        if seq.len() < 7 {
            continue;
        }

        if config.verbose {
            eprintln!(
                "  Working Realignlgdel: 3' {} '{}' {}",
                p,
                String::from_utf8_lossy(&seq),
                cnt
            );
        }

        let mut bp = find_bp(
            &seq,
            p + 5,
            reference,
            1,
            chr,
            config.indel_size as i32,
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
            if !(bp != 0 && bp - p > 15 && p - bp < 150000) {
                continue;
            }

            let mut lists = [&mut *svfdel, &mut *svrdel];
            let sv_mark = crate::sv::mark_sv(p, bp, &mut lists, max_read_length, config.verbose);
            svcov = sv_mark.0;
            let clusters = sv_mark.1;
            let pairs = sv_mark.2;

            if svcov == 0 && cnt <= config.min_reads as i32 {
                continue;
            }

            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, p);
            sv.r#type = "DEL".to_string();
            sv.pairs += pairs;
            sv.splits += cnt;
            sv.clusters += clusters;

            if bp > region.end.0 as i32 {
                let mut tts = bp - max_read_length;
                let tte = bp + max_read_length;
                if bp - max_read_length <= region.end.0 as i32 {
                    tts = region.end.0 as i32 + 1;
                }
                if let Err(e) = ensure_reference_loaded(reference, chr, tts, tte, config) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }
        }

        let mut dellen = bp - p;
        let mut en = 0;
        if !extra_str.is_empty() {
            dellen -= extra_str.len() as i32;
        } else {
            while en < seq.len() {
                let ref_pos = bp + en as i32;
                let ref_ch = reference.get_base(chr, ref_pos as u32);
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
            sc5p += extra.len() as i32;
        } else if !extra_str.is_empty() {
            let extra_str_val = String::from_utf8_lossy(&extra_str).into_owned();
            gt = format!("-{}&{}", dellen, extra_str_val).into_bytes();
        } else {
            loop {
                let base_bp_prev = reference.get_base(chr, (bp - 1) as u32);
                let base_bp_del_prev = reference.get_base(chr, (bp + dellen - 1) as u32);
                if base_bp_prev.is_some()
                    && base_bp_del_prev.is_some()
                    && base_bp_prev == base_bp_del_prev
                {
                    bp -= 1;
                    if bp != 0 {
                        sc5p -= 1;
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
                    let dest =
                        non_insertion_variants.compute_if_absent(bp, |_| VariationMap::new());
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
                && cnt <= config.min_reads as i32
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
            eprintln!(
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

        if dellen < config.indel_size as i32 {
            for tp in bp..bp + dellen + extra.len() as i32 + extra_str.len() as i32 {
                inc_cnt(ref_coverage, tp, sc3v.variation.cnt);
            }
        }

        if !ref_coverage.contains_key(&bp) {
            if let Some(&cov_p_prev) = ref_coverage.get(&(p - 1)) {
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
        sc3v.used = p + dellen != 0;

        if let Some(sc) = soft_clips_3_end.get_mut(&p) {
            sc.used = p + dellen != 0;
            sc.sequence = sc3v.sequence;
            sc.variation.pmean = sc3v.variation.pmean;
        }

        let mut dels5 = BTreeMap::new();
        let mut map = std::collections::HashMap::new();
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
                        sv.splits += tv_cnt - dels5_cnt;
                    }
                }
            }
        }

        if config.verbose {
            let tv_cnt = crate::variations::get_variation(non_insertion_variants, bp, &gt_str).cnt;
            eprintln!(
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
            add_var_factor(tv, (svcov - tv_cnt) as f64 / tv_cnt as f64);
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
    ref_coverage: &mut BTreeMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_maps(&config.adaptor);
    let extension = 5000;

    let mut tmp5 = Vec::new();
    for (&p, sc) in soft_clips_5_end.iter() {
        if p < region.start.0 as i32 - extension || p > region.end.0 as i32 + extension {
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
        if p < region.start.0 as i32 - extension || p > region.end.0 as i32 + extension {
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

    for (p5, cnt5) in tmp5 {
        if soft_clips_5_end.get(&p5).is_none_or(|sc| sc.used) {
            continue;
        }
        for (p3, cnt3) in tmp3.iter().copied() {
            if soft_clips_5_end.get(&p5).is_none_or(|sc| sc.used) {
                break;
            }
            if soft_clips_3_end.get(&p3).is_none_or(|sc| sc.used) {
                continue;
            }
            if p5 - p3 > (max_read_length as f64 * 2.5) as i32 {
                continue;
            }
            if p3 - p5 > max_read_length - 10 {
                continue;
            }
            let seq5 = {
                let sc5v = match soft_clips_5_end.get_mut(&p5) {
                    Some(sc) => sc,
                    None => continue,
                };
                crate::variations::util::find_conseq(sc5v, 5, &adaptor_forward, &adaptor_reverse)
            };
            let seq3 = {
                let sc3v = match soft_clips_3_end.get_mut(&p3) {
                    Some(sc) => sc,
                    None => continue,
                };
                crate::variations::util::find_conseq(sc3v, 3, &adaptor_forward, &adaptor_reverse)
            };
            if seq5.len() <= 10 || seq3.len() <= 10 {
                continue;
            }
            if config.verbose {
                eprintln!(
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
            let bp5 = match35.matched_5_end;
            let bp3 = match35.matched_3_end;
            let score = match35.max_matched_length;
            if score == 0 {
                continue;
            }
            let smscore = score / 2;
            let mut ins = if bp3 + smscore > 1 {
                substr(&seq3, 0, (-(bp3 + smscore) + 1) as isize).to_vec()
            } else {
                seq3.clone()
            };
            if bp5 + smscore > 0 {
                ins.extend_from_slice(&reverse(substr(&seq5, 0, (bp5 + smscore) as isize)));
            }
            if crate::variations::util::is_low_complex_seq(&ins) {
                if config.verbose {
                    eprintln!(
                        "  Discard low complexity insertion found {}.",
                        String::from_utf8_lossy(&ins)
                    );
                }
                continue;
            }
            let bi;
            let ins_str;
            if config.verbose {
                eprintln!(
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
                        p5 + seq3.len() as i32 - ins.len() as i32 + 2,
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
                        p3 - (seq5.len() as i32 - ins.len() as i32) - 2,
                        p3 - 1,
                    );
                    if !is_match(sub, &jref, -1, 3) {
                        continue;
                    }
                }
                if config.verbose {
                    eprintln!(
                        "  Found lgins30 complex: {} {} {} {}",
                        p3,
                        p5,
                        ins.len(),
                        String::from_utf8_lossy(&ins)
                    );
                }
                let tmp = join_ref_int(reference, chr, p3, p5 - 1);
                if tmp.len() > ins.len() {
                    ins_str = format!("{}^{}", p3 - p5, String::from_utf8_lossy(&ins));
                    bi = p3;
                } else if tmp.len() < ins.len() {
                    let part1 = substr(&ins, 0, ins.len() as isize - tmp.len() as isize);
                    let part2 = substr_from(&ins, (p3 - p5) as isize);
                    ins_str = format!(
                        "+{}&{}",
                        String::from_utf8_lossy(part1),
                        String::from_utf8_lossy(part2)
                    );
                    bi = p3 - 1;
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
                        p5 + seq3.len() as i32 - ins.len() as i32 + 2,
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
                        p3 - (seq5.len() as i32 - ins.len() as i32) - 2,
                        p3 - 1,
                    );
                    if !is_match(sub, &jref, -1, 3) {
                        continue;
                    }
                }
                let mut tmp = Vec::new();
                if ins.len() <= (p3 - p5) as usize {
                    let mut rpt = 2;
                    let mut tnr = 3;
                    while (((p3 - p5 + ins.len() as i32) as f64 / tnr as f64) / ins.len() as f64)
                        > 1.0
                    {
                        if (p3 - p5 + ins.len() as i32) % tnr == 0 {
                            rpt += 1;
                        }
                        tnr += 1;
                    }
                    let limit = p5 as f64 + (p3 - p5 + ins.len() as i32) as f64 / rpt as f64
                        - ins.len() as f64;
                    tmp.extend_from_slice(&join_ref_double(reference, chr, p5, limit));
                    ins_str = format!(
                        "+{}{}",
                        String::from_utf8_lossy(&tmp),
                        String::from_utf8_lossy(&ins)
                    );
                } else {
                    tmp.extend_from_slice(&join_ref_int(reference, chr, p5, p3 - 1));
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
                    eprintln!(
                        "Found lgins30: {} {} {} {} + {}",
                        p3,
                        p5,
                        ins.len(),
                        String::from_utf8_lossy(&tmp),
                        ins_str
                    );
                }
                bi = p5 - 1;
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
                eprintln!(" lgins30 Found: '{}' {} {} {}", ins_str, bi, bp3, bp5);
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
                let mvref_cloned = mvref_mut.as_deref().cloned();

                // Mutate mvref in adj_cnt for sc3v:
                let vref = get_variation(insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc3v_var, mvref_mut);

                // Mutate for sc5v:
                let vref = get_variation(insertion_variants, bi, &ins_str);
                crate::variations::util::adj_cnt(vref, &sc5v_var, None);

                let has_bams = bams.is_some_and(|b| !b.is_empty());
                let mvref_cnt = mvref_cloned.as_ref().map_or(0, |m| m.cnt);
                let bams_slice = match bams {
                    Some(b) => b,
                    None => &[],
                };
                let vref = get_variation(insertion_variants, bi, &ins_str);
                if has_bams
                    && p3 - p5 >= 5
                    && p3 - p5 < max_read_length - 10
                    && mvref_cnt != 0
                    && no_passing_reads(chr, p5, p3, bams_slice, config)
                    && vref.cnt > 2 * mvref_cnt
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
                let mut map = std::collections::HashMap::new();
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
                    mvref.cnt -= sc3v_var.cnt;
                    mvref.hicnt -= sc3v_var.hicnt;
                    mvref.locnt -= sc3v_var.locnt;
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
                let mut map = std::collections::HashMap::new();
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
        eprintln!("Done: lgins30\n");
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
    ref_coverage: &mut BTreeMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
    max_read_length: i32,
    svfdup: &mut Vec<Sclip>,
    svrdup: &mut Vec<Sclip>,
) {
    let (adaptor_forward, adaptor_reverse) = fill_adaptor_maps(&config.adaptor);
    let extension = 5000;

    let mut tmp = Vec::new();
    for (&p, sc5v) in soft_clips_5_end.iter() {
        if p < region.start.0 as i32 - extension || p > region.end.0 as i32 + extension {
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
        if cnt < config.min_reads as i32 {
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
            crate::variations::util::find_conseq(sc5v, 0, &adaptor_forward, &adaptor_reverse)
        };
        if seq.is_empty() {
            continue;
        }
        if config.verbose {
            eprintln!(
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
            if !(bi != 0 && bi - p > 15 && bi - p < 150000) {
                continue;
            }

            if bi > region.end.0 as i32 {
                let mut tts = bi - max_read_length;
                let tte = bi + max_read_length;
                if bi - max_read_length <= region.end.0 as i32 {
                    tts = region.end.0 as i32 + 1;
                }
                if let Err(e) = ensure_reference_loaded(reference, chr, tts, tte, config) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }

            const SVFLANK: i32 = 50;
            if bi - p > config.sv_min_len as i32 + 2 * SVFLANK {
                let mut tmp_ins = join_ref_int(reference, chr, p, p + SVFLANK - 1);
                let dup_str = format!("<dup{}>", bi - p - 2 * SVFLANK + 1);
                tmp_ins.extend_from_slice(dup_str.as_bytes());
                let jref = join_ref_for_5_lgins(reference, chr, bi - SVFLANK + 1, bi, &seq, &extra);
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
            let cov_p_minus_1 = ref_coverage.get(&(p - 1)).copied();

            if cov_p_minus_1.is_none()
                || cov_p_minus_1.is_some_and(|cp| cov_bi.is_some_and(|cb| cp < cb))
            {
                if let Some(c_bi) = cov_bi {
                    ref_coverage.insert(p - 1, c_bi);
                } else {
                    ref_coverage.insert(p - 1, sc5v_vars_count);
                }
            } else {
                if sc5v_vars_count > cov_p_minus_1.unwrap_or(0) {
                    inc_cnt(ref_coverage, p - 1, sc5v_vars_count);
                }
            }

            bi = p - 1;
            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, bi);
            sv.r#type = "DUP".to_string();
            sv.pairs += pairs;
            sv.splits += cnt;
            sv.clusters += clusters;
        }

        let ins_str = String::from_utf8_lossy(&ins).into_owned();
        if config.verbose {
            eprintln!(
                "  Found candidate lgins from 5: {} +{} {} {}",
                bi,
                ins_str,
                p,
                String::from_utf8_lossy(&seq)
            );
        }
        let iref = get_variation(insertion_variants, bi, &format!("+{}", ins_str));
        iref.pstd = true;
        iref.qstd = true;

        let sc5v_var = soft_clips_5_end
            .get(&p)
            .map_or(Variation::new(), |sc| sc.variation.clone());
        let iref = get_variation(insertion_variants, bi, &format!("+{}", ins_str));
        crate::variations::util::adj_cnt(iref, &sc5v_var, None);

        let mut rpflag = true;
        for (i, &ch) in ins.iter().enumerate() {
            let ref_ch = reference.get_base(chr, (bi + 1 + i as i32) as u32);
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
            + 1;

        let entries_to_process: Vec<(i32, u8, Variation)> = {
            let sc5v = soft_clips_5_end.get(&p);
            sc5v.map_or_else(Vec::new, |sc5v| {
                let mut res = Vec::new();
                for ii in (len + 1)..seq_len {
                    let pii = bi - ii + len;
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
            let tvr = get_variation(non_insertion_variants, pii, &format!("{}", tnt as char));
            crate::variations::util::adj_cnt(tvr, &tv, None);
            tvr.pstd = true;
            tvr.qstd = true;
            inc_cnt(ref_coverage, pii, tv.cnt);
        }

        if let Some(sc5v) = soft_clips_5_end.get_mut(&p) {
            sc5v.used = bi + len != 0;
        }

        let tins_vars_cnt = {
            let iref = get_variation(insertion_variants, bi, &format!("+{}", ins_str));
            iref.cnt
        };
        let mut tins = BTreeMap::new();
        let mut map = std::collections::HashMap::new();
        map.insert(format!("+{}", ins_str), tins_vars_cnt);
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
        let newins = if newins.is_empty() {
            format!("+{}", ins_str)
        } else {
            newins
        };

        if let Some(entry) = non_insertion_variants.get_mut(&bi) {
            if entry.sv.is_some() {
                let kref = get_variation(insertion_variants, bi, &newins);
                let diff = kref.cnt - tins_vars_cnt;
                if let Some(sv) = entry.sv.as_mut() {
                    sv.splits += diff;
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
                && ins.len() < max_read_length as usize - 10
                && mref.cnt != 0
                && no_passing_reads(chr, bi, bi + ins.len() as i32, bams_slice, config)
            {
                let m_cnt = mref.cnt;
                let kref = get_variation(insertion_variants, bi, &newins);
                if kref.cnt > 2 * m_cnt {
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
        if p < region.start.0 as i32 - extension || p > region.end.0 as i32 + extension {
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
        if cnt < config.min_reads as i32 {
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
            crate::variations::util::find_conseq(sc3v, 0, &adaptor_forward, &adaptor_reverse)
        };
        if seq.is_empty() {
            continue;
        }
        if config.verbose {
            eprintln!(
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

        if bi == 0 {
            if crate::variations::util::is_low_complex_seq(&seq) {
                continue;
            }
            let match_res =
                crate::sv::find_match(&seq, reference, p, 1, 17, 1, chr, config.verbose);
            bi = match_res.base_position;
            let extra = match_res.matched_sequence;
            if !(bi != 0 && p - bi > 15 && p - bi < 150000) {
                continue;
            }

            if bi < region.start.0 as i32 {
                let tts = bi - max_read_length;
                let mut tte = bi + max_read_length;
                if bi + max_read_length >= region.start.0 as i32 {
                    tte = region.start.0 as i32 - 1;
                }
                if let Err(e) = ensure_reference_loaded(reference, chr, tts, tte, config) {
                    log::warn!("ensure_reference_loaded failed: {:?}", e);
                }
                if let Err(e) = partial_pipeline(
                    bams,
                    chr,
                    tts,
                    tte,
                    reference,
                    config,
                    non_insertion_variants,
                    insertion_variants,
                    ref_coverage,
                    soft_clips_5_end,
                    soft_clips_3_end,
                ) {
                    log::warn!("partial_pipeline failed: {:?}", e);
                }
            }

            let mut p_mut = p;
            let mut bi_mut = bi;
            let mut shift5 = 0;
            loop {
                let ref_p_minus_1 = reference.get_base(chr, (p_mut - 1) as u32);
                let ref_bi_minus_1 = reference.get_base(chr, (bi_mut - 1) as u32);
                if ref_p_minus_1.is_some()
                    && ref_bi_minus_1.is_some()
                    && ref_p_minus_1 == ref_bi_minus_1
                {
                    p_mut -= 1;
                    bi_mut -= 1;
                    shift5 += 1;
                } else {
                    break;
                }
            }

            const SVFLANK: i32 = 50;
            if p_mut - bi_mut > config.sv_min_len as i32 + 2 * SVFLANK {
                ins = join_ref_for_3_lgins(
                    reference,
                    chr,
                    bi_mut,
                    bi_mut + SVFLANK - 1,
                    shift5,
                    &seq,
                    &extra,
                );
                let dup_str = format!("<dup{}>", p_mut - bi_mut - 2 * SVFLANK);
                ins.extend_from_slice(dup_str.as_bytes());
                let jref = join_ref_int(reference, chr, p_mut - SVFLANK, p_mut - 1);
                ins.extend_from_slice(&jref);
            } else {
                ins = join_ref_for_3_lgins(reference, chr, bi_mut, p_mut - 1, shift5, &seq, &extra);
            }
            ins.extend_from_slice(&extra);

            let mut sv_lists = [&mut *svfdup, &mut *svrdup];
            let (clusters, pairs) = crate::sv::mark_dup_sv(
                bi_mut,
                p_mut - 1,
                &mut sv_lists,
                max_read_length,
                config.verbose,
            );
            bi = bi_mut - 1;

            let sv = crate::variations::VariationMap::get_sv(non_insertion_variants, bi);
            sv.r#type = "DUP".to_string();
            sv.pairs += pairs;
            sv.splits += cnt;
            sv.clusters += clusters;

            let sc3v_vars_count = soft_clips_3_end.get(&p).map_or(0, |sc| sc.variation.cnt);
            let cov_p = ref_coverage.get(&p_mut).copied();
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

        let ins_str = String::from_utf8_lossy(&ins).into_owned();
        if config.verbose {
            eprintln!(
                "  Found candidate lgins from 3: {} +{} {} {}",
                bi,
                ins_str,
                p,
                String::from_utf8_lossy(&seq)
            );
        }
        let iref = get_variation(insertion_variants, bi, &format!("+{}", ins_str));
        iref.pstd = true;
        iref.qstd = true;

        let ref_base = reference.get_base(chr, bi as u32);
        let mut lref_mut = get_variation_maybe_mut(non_insertion_variants, bi, ref_base);
        let sc3v_mean_pos = soft_clips_3_end
            .get(&p)
            .map_or(0.0, |sc| sc.variation.pmean);
        if p as f64 - bi as f64 > sc3v_mean_pos / (cnt as f64) {
            lref_mut = None;
        }

        let sc3v_var = soft_clips_3_end
            .get(&p)
            .map_or(Variation::new(), |sc| sc.variation.clone());
        let iref = get_variation(insertion_variants, bi, &format!("+{}", ins_str));
        if let Some(ref_mut) = lref_mut {
            crate::variations::util::adj_cnt(iref, &sc3v_var, Some(ref_mut));
        } else {
            crate::variations::util::adj_cnt(iref, &sc3v_var, None);
        }

        let mut rpflag = true;
        for (i, &ch) in ins.iter().enumerate() {
            let ref_ch = reference.get_base(chr, (bi + 1 + i as i32) as u32);
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
            .get(&p)
            .map_or(0, |sc| sc.seq.keys().max().copied().unwrap_or(0))
            + 1;

        let entries_to_process: Vec<(i32, u8, Variation)> = {
            let sc3v = soft_clips_3_end.get(&p);
            sc3v.map_or_else(Vec::new, |sc3v| {
                let mut res = Vec::new();
                for ii in len..len_seq {
                    let pii = p + ii - len;
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
            let tvr = get_variation(non_insertion_variants, pii, &format!("{}", tnt as char));
            crate::variations::util::adj_cnt(tvr, &tv, None);
            tvr.pstd = true;
            tvr.qstd = true;
            inc_cnt(ref_coverage, pii, tv.cnt);
        }

        if let Some(sc3v) = soft_clips_3_end.get_mut(&p) {
            sc3v.used = true;
        }

        let tins_vars_cnt = {
            let iref = get_variation(insertion_variants, bi, &format!("+{}", ins_str));
            iref.cnt
        };
        let mut tins = BTreeMap::new();
        let mut map = std::collections::HashMap::new();
        map.insert(format!("+{}", ins_str), tins_vars_cnt);
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
        let newins = if newins.is_empty() {
            format!("+{}", ins_str)
        } else {
            newins
        };

        if let Some(entry) = non_insertion_variants.get_mut(&bi) {
            if entry.sv.is_some() {
                let iref = get_variation(insertion_variants, bi, &newins);
                let diff = iref.cnt - tins_vars_cnt;
                if let Some(sv) = entry.sv.as_mut() {
                    sv.splits += diff;
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
                && ins.len() < max_read_length as usize - 10
                && mref.cnt != 0
                && no_passing_reads(chr, bi, bi + ins.len() as i32, bams_slice, config)
            {
                let iref = get_variation(insertion_variants, bi, &newins);
                let m_cnt = mref.cnt;
                if iref.cnt > 2 * m_cnt {
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
        eprintln!("Start local realignment...");
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
    );

    if config.verbose {
        eprintln!("Finished local realignment.");
    }
}

fn substr(seq: &[u8], mut begin: isize, len: isize) -> &[u8] {
    let s_len = seq.len() as isize;
    if begin < 0 {
        begin = s_len + begin;
    }
    if begin < 0 {
        begin = 0;
    }
    if begin >= s_len {
        return b"";
    }

    if len > 0 {
        let end = (begin + len).min(s_len) as usize;
        &seq[begin as usize..end]
    } else if len == 0 {
        b""
    } else {
        let end = s_len + len;
        if end < begin {
            b""
        } else {
            &seq[begin as usize..end as usize]
        }
    }
}

fn substr_from(seq: &[u8], idx: isize) -> &[u8] {
    let s_len = seq.len() as isize;
    if idx >= 0 {
        let start = idx.min(s_len) as usize;
        &seq[start..]
    } else {
        let start = (s_len + idx).max(0) as usize;
        &seq[start..]
    }
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

fn join_ref_for_5_lgins(
    reference: &Reference,
    chr: &str,
    from: i32,
    to: i32,
    seq: &[u8],
    extra: &[u8],
) -> Vec<u8> {
    let mut sb = Vec::new();
    for i in from..=to {
        let idx = (to - i) as usize + extra.len();
        if idx < seq.len() {
            sb.push(seq[idx]);
        } else {
            if i >= 1 {
                if let Some(ch) = reference.get_base(chr, i as u32) {
                    sb.push(ch);
                }
            }
        }
    }
    sb
}

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
    for i in from..=to {
        let diff = i - from;
        if diff >= shift5 && ((diff - shift5) as usize) < seq.len() - extra.len() {
            let idx = (diff - shift5) as usize + extra.len();
            if idx < seq.len() {
                sb.push(seq[idx]);
            }
        } else {
            if i >= 1 {
                if let Some(ch) = reference.get_base(chr, i as u32) {
                    sb.push(ch);
                }
            }
        }
    }
    sb
}

pub struct Match35 {
    pub matched_5_end: i32,
    pub matched_3_end: i32,
    pub max_matched_length: i32,
}

pub fn find_35_match(seq5: &[u8], seq3: &[u8], verbose: bool) -> Match35 {
    let long_mismatch = 2;
    let mut max_matched_length = 0;
    let mut b3 = 0;
    let mut b5 = 0;

    let len5 = seq5.len();
    let len3 = seq3.len();

    if len5 < 8 || len3 < 8 {
        return Match35 {
            matched_5_end: 0,
            matched_3_end: 0,
            max_matched_length: 0,
        };
    }

    for i in 0..(len5 - 8) {
        for j in 1..(len3 - 8) {
            let mut number_of_mismatch = 0;
            let mut total_length = 0;
            while total_length + j <= len3 && i + total_length <= len5 {
                let match_char = if i + total_length < len5 {
                    let idx3 = len3 - j - total_length;
                    seq3[idx3] == seq5[i + total_length]
                } else {
                    false
                };
                if !match_char {
                    number_of_mismatch += 1;
                }
                if number_of_mismatch > long_mismatch {
                    break;
                }
                total_length += 1;
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
                    eprintln!(
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

fn inc_cnt(hash: &mut BTreeMap<i32, i32>, start: i32, count: i32) {
    *hash.entry(start).or_insert(0) += count;
}

fn rm_cnt(vref: &mut crate::variations::Variation, tv: &crate::variations::Variation) {
    vref.cnt -= tv.cnt;
    vref.hicnt -= tv.hicnt;
    vref.locnt -= tv.locnt;
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
    vref.cnt += (factor_f * vref.cnt as f64) as i32;
    vref.hicnt += (factor_f * vref.hicnt as f64) as i32;
    vref.locnt += (factor_f * vref.locnt as f64) as i32;
    vref.pmean += factor_f * vref.pmean;
    vref.qmean += factor_f * vref.qmean;
    vref.mqmean += factor_f * vref.mqmean;
    vref.nm += factor_f * vref.nm;
    vref.fwd += (factor_f * vref.fwd as f64) as i32;
    vref.rev += (factor_f * vref.rev as f64) as i32;
}

fn ensure_reference_loaded(
    reference: &mut Reference,
    chr: &str,
    start: i32,
    end: i32,
    config: &crate::config::Config,
) -> Result<(), crate::types::VardictError> {
    if crate::reference::is_loaded(chr, start as u32, end as u32, reference) {
        return Ok(());
    }
    let mut min_start = start.max(1) as u32;
    let mut max_end = end.max(1) as u32;
    for lr in &reference.loaded_regions {
        if lr.chr == chr {
            min_start = min_start.min(lr.sequence_start);
            max_end = max_end.max(lr.sequence_end);
        }
    }
    let fasta_path = std::path::Path::new(&config.reference);
    let mut bases =
        crate::reference::fetch_subsequence(fasta_path, chr, Pos1(min_start), Pos1(max_end))?;
    for b in &mut bases {
        b.make_ascii_uppercase();
    }
    if let Some(idx) = reference
        .loaded_regions
        .iter()
        .position(|lr| lr.chr == *chr)
    {
        reference.loaded_regions[idx].sequence_start = min_start;
        reference.loaded_regions[idx].sequence_end = max_end;
        if idx < reference.region_sequences.len() {
            reference.region_sequences[idx] = bases.clone();
        }
    } else {
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: min_start,
                sequence_end: max_end,
            });
        reference.region_sequences.push(bases.clone());
    }
    reference.reference_sequences.insert(chr.to_string(), bases);
    Ok(())
}

fn merge_sclip(dest: &mut Sclip, src: Sclip) {
    let d = &mut dest.variation;
    d.cnt += src.variation.cnt;
    d.fwd += src.variation.fwd;
    d.rev += src.variation.rev;
    d.pmean += src.variation.pmean;
    d.qmean += src.variation.qmean;
    d.mqmean += src.variation.mqmean;
    d.nm += src.variation.nm;
    d.hicnt += src.variation.hicnt;
    d.locnt += src.variation.locnt;
    d.extracnt += src.variation.extracnt;
    if src.variation.pstd {
        d.pstd = true;
    }
    if src.variation.qstd {
        d.qstd = true;
    }
    // pp/pq are last-seen trackers (not sums); no consumer reads them on a
    // merged sclip, so accumulating them here is inert.
    d.pp += src.variation.pp;
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
    dest.disc += src.disc;
    for (pos, cnt) in src.soft {
        *dest.soft.entry(pos).or_insert(0) += cnt;
    }
    dest.mates.extend(src.mates);

    for (idx, cnts) in src.nt {
        let dest_cnts = dest.nt.entry(idx).or_default();
        for (base, c) in cnts {
            *dest_cnts.entry(base).or_insert(0) += c;
        }
    }
    for (idx, bases) in src.seq {
        let dest_bases = dest.seq.entry(idx).or_default();
        for (base, var) in bases {
            let dv = dest_bases.entry(base).or_insert_with(Variation::new);
            dv.cnt += var.cnt;
            dv.fwd += var.fwd;
            dv.rev += var.rev;
            dv.pmean += var.pmean;
            dv.qmean += var.qmean;
            dv.mqmean += var.mqmean;
            dv.nm += var.nm;
            dv.hicnt += var.hicnt;
            dv.locnt += var.locnt;
            dv.extracnt += var.extracnt;
            if var.pstd {
                dv.pstd = true;
            }
            if var.qstd {
                dv.qstd = true;
            }
            dv.pp += var.pp;
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
    non_insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    insertion_variants: &mut crate::javamap::JavaHashMap<i32, VariationMap>,
    ref_coverage: &mut BTreeMap<i32, i32>,
    soft_clips_5_end: &mut BTreeMap<i32, Sclip>,
    soft_clips_3_end: &mut BTreeMap<i32, Sclip>,
) -> Result<(), crate::types::VardictError> {
    let bams = match bams {
        Some(b) => b,
        None => return Ok(()),
    };
    if bams.is_empty() {
        return Ok(());
    }

    let mut records = Vec::new();
    let mut preprocessor = crate::bam::RecordPreprocessor::new();
    for bam in bams {
        let bam_path = std::path::Path::new(bam);
        let raw_records = crate::bam::query_region(
            bam_path,
            chr,
            Pos0(tts.max(0) as u32),
            Pos0(tte as u32),
            config,
        )?;
        for r in raw_records {
            if preprocessor.preprocess_record(&r, config) {
                records.push(r);
            }
        }
    }

    let pileup_region = crate::regions::Region {
        chr: chr.to_string(),
        start: Pos1(tts as u32),
        end: Pos1(tte as u32),
        gene: String::new(),
        amplicon: None,
    };

    let var_data = crate::cigar::parse::pileup(
        &records,
        reference,
        &pileup_region,
        config,
        true,
        preprocessor.total_reads as i32,
        preprocessor.duplicate_reads as i32,
    )?;

    for (pos, cov) in var_data.ref_coverage {
        *ref_coverage.entry(pos).or_insert(0) += cov;
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
        let dest = non_insertion_variants.compute_if_absent(pos, |_| VariationMap::new());
        for (allele, var) in var_map.map {
            dest.map
                .entry(allele)
                .and_modify(|v| {
                    v.cnt += var.cnt;
                    v.fwd += var.fwd;
                    v.rev += var.rev;
                    v.pmean += var.pmean;
                    v.qmean += var.qmean;
                    v.mqmean += var.mqmean;
                    v.nm += var.nm;
                    v.hicnt += var.hicnt;
                    v.locnt += var.locnt;
                    v.extracnt += var.extracnt;
                    if var.pstd {
                        v.pstd = true;
                    }
                    if var.qstd {
                        v.qstd = true;
                    }
                    v.pp += var.pp;
                    v.pq += var.pq;
                })
                .or_insert(var);
        }
        if let Some(sv) = var_map.sv {
            if let Some(dest_sv) = &mut dest.sv {
                dest_sv.pairs += sv.pairs;
                dest_sv.splits += sv.splits;
                dest_sv.clusters += sv.clusters;
            } else {
                dest.sv = Some(sv);
            }
        }
    }
    for (pos, var_map) in var_data.insertion_variants {
        let dest = insertion_variants.compute_if_absent(pos, |_| VariationMap::new());
        for (allele, var) in var_map.map {
            dest.map
                .entry(allele)
                .and_modify(|v| {
                    v.cnt += var.cnt;
                    v.fwd += var.fwd;
                    v.rev += var.rev;
                    v.pmean += var.pmean;
                    v.qmean += var.qmean;
                    v.mqmean += var.mqmean;
                    v.nm += var.nm;
                    v.hicnt += var.hicnt;
                    v.locnt += var.locnt;
                    v.extracnt += var.extracnt;
                    if var.pstd {
                        v.pstd = true;
                    }
                    if var.qstd {
                        v.qstd = true;
                    }
                    v.pp += var.pp;
                    v.pq += var.pq;
                })
                .or_insert(var);
        }
        if let Some(sv) = var_map.sv {
            if let Some(dest_sv) = &mut dest.sv {
                dest_sv.pairs += sv.pairs;
                dest_sv.splits += sv.splits;
                dest_sv.clusters += sv.clusters;
            } else {
                dest.sv = Some(sv);
            }
        }
    }

    Ok(())
}

/// Adjust MNP when there are breakpoints within MNP.
/// Ports `VariationRealigner.java:277-382`.
pub fn adjust_mnp(
    reference: &Reference,
    chr: &str,
    _chr_len: u32,
    _region: &crate::regions::Region,
    config: &crate::config::Config,
    var_data: &mut crate::variations::VariationData,
) {
    let tmp = fill_and_sort_tmp(&var_data.mnp);
    let adaptor_forward = BTreeMap::new();
    let adaptor_reverse = BTreeMap::new();

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
            eprintln!("  AdjMnt: {} {} {}", position, vn, vref.cnt);
        }

        let mnt = vn.replacen('&', "", 1);
        let mnt_len = mnt.len();
        for i in 0..mnt_len.saturating_sub(1) {
            let left_len = i + 1;
            let left_bytes = substr(mnt.as_bytes(), 0, left_len as isize);
            let mut left = String::from_utf8_lossy(left_bytes).into_owned();
            if left.len() > 1 {
                left.insert(1, '&');
            }

            let right_len = mnt_len - i - 1;
            let right_bytes = substr(mnt.as_bytes(), -(right_len as isize), right_len as isize);
            let mut right = String::from_utf8_lossy(right_bytes).into_owned();
            if right.len() > 1 {
                right.insert(1, '&');
            }

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
                    eprintln!(
                        "    AdjMnt Left: {} {} Left: {} Cnt: {}",
                        position, vn, left, tref_clone.cnt
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
            let right_pos = position + i as i32 + 1;
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
                    eprintln!(
                        "    AdjMnt Right: {} {} Right: {} Cnt: {}",
                        position, vn, right, tref_clone.cnt
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

        // --- SOFTCLIPS3END CHECK ---
        let mut sc3_vars_to_add = None;
        {
            if let Some(sc3v) = var_data.soft_clips_3_end.get_mut(&position) {
                if !sc3v.used {
                    let seq = crate::variations::util::find_conseq(
                        sc3v,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    let seq_str = String::from_utf8_lossy(&seq);
                    if seq_str.starts_with(&mnt) {
                        if seq.len() == mnt.len()
                            || is_match_ref(
                                &seq[mnt.len()..],
                                reference,
                                chr,
                                position + mnt.len() as i32,
                                1,
                                3,
                            )
                        {
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
        let right_pos = position + mnt_len as i32;
        let mut sc5_vars_to_add = None;
        {
            if let Some(sc5v) = var_data.soft_clips_5_end.get_mut(&right_pos) {
                if !sc5v.used {
                    let seq = crate::variations::util::find_conseq(
                        sc5v,
                        0,
                        &adaptor_forward,
                        &adaptor_reverse,
                    );
                    if !seq.is_empty() && seq.len() >= mnt_len {
                        let mut seq_rev = seq.clone();
                        seq_rev.reverse();
                        let seq_rev_str = String::from_utf8_lossy(&seq_rev);
                        if seq_rev_str.ends_with(&mnt) {
                            if seq.len() == mnt_len
                                || is_match_ref(
                                    &seq_rev[0..seq_rev.len() - mnt_len],
                                    reference,
                                    chr,
                                    position - 1,
                                    -1,
                                    3,
                                )
                            {
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
    let adaptor_forward = BTreeMap::new();
    let adaptor_reverse = BTreeMap::new();

    let keys_5: Vec<i32> = var_data.soft_clips_5_end.keys().copied().collect();
    for position in keys_5 {
        let sc5v = match var_data.soft_clips_5_end.get_mut(&position) {
            Some(sc) => sc,
            None => continue,
        };
        if sc5v.used {
            continue;
        }
        let seq = crate::variations::util::find_conseq(sc5v, 0, &adaptor_forward, &adaptor_reverse);
        if seq.len() > 5 {
            continue;
        }
        let bp_bytes = substr(&seq, 0, 1);
        let bp = String::from_utf8_lossy(bp_bytes).into_owned();

        let previous_position = position - 1;
        let mut should_adj = false;
        let sc5v_vars_count = sc5v.variation.cnt;

        if let Some(vars_on_prev) = var_data.non_insertion_variants.get_mut(&previous_position) {
            if vars_on_prev.map.contains_key(&bp) {
                should_adj = true;
                if seq.len() > 1 {
                    if let Some(pos_sub_2) = position.checked_sub(2) {
                        let ref_base = reference.get_base(chr, pos_sub_2 as u32);
                        if ref_base.is_none_or(|rb| rb != seq[1]) {
                            should_adj = false;
                        }
                    } else {
                        should_adj = false;
                    }
                }
                if should_adj {
                    let sc5v_clone = sc5v.clone();
                    if let Some(vref_mut) = vars_on_prev.map.get_mut(&bp) {
                        crate::variations::util::adj_cnt(vref_mut, &sc5v_clone.variation, None);
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
        let seq = crate::variations::util::find_conseq(sc3v, 0, &adaptor_forward, &adaptor_reverse);
        if seq.len() > 5 {
            continue;
        }
        let bp_bytes = substr(&seq, 0, 1);
        let bp = String::from_utf8_lossy(bp_bytes).into_owned();

        let mut should_adj = false;
        let sc3v_vars_count = sc3v.variation.cnt;

        if let Some(vars_on_pos) = var_data.non_insertion_variants.get_mut(&position) {
            if vars_on_pos.map.contains_key(&bp) {
                should_adj = true;
                if seq.len() > 1 {
                    if let Some(pos_add_1) = position.checked_add(1) {
                        let ref_base = reference.get_base(chr, pos_add_1 as u32);
                        if ref_base.is_none_or(|rb| rb != seq[1]) {
                            should_adj = false;
                        }
                    } else {
                        should_adj = false;
                    }
                }
                if should_adj {
                    let sc3v_clone = sc3v.clone();
                    if let Some(vref_mut) = vars_on_pos.map.get_mut(&bp) {
                        crate::variations::util::adj_cnt(vref_mut, &sc3v_clone.variation, None);
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

#[cfg(test)]
mod tests {
    use crate::reference::LoadedRegion;

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
            .compute_if_absent(50, |_| crate::variations::VariationMap::new())
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
            .compute_if_absent(50, |_| crate::variations::VariationMap::new())
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
            .compute_if_absent(51, |_| crate::variations::VariationMap::new())
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
            .compute_if_absent(50, |_| crate::variations::VariationMap::new())
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
            .compute_if_absent(49, |_| crate::variations::VariationMap::new())
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
        let mut del_map = std::collections::HashMap::new();
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

        let mut ref_coverage = BTreeMap::new();
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
        let mut ins_map = std::collections::HashMap::new();
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

        let mut ref_coverage = BTreeMap::new();
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
        let mut ref_coverage = BTreeMap::new();
        ref_coverage.insert(50, 10);
        let mut soft_clips_3_end = BTreeMap::new();

        let mut svfdel = Vec::new();
        let mut svrdel = Vec::new();

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
        let mut ref_coverage = BTreeMap::new();
        ref_coverage.insert(50, 10);
        let mut soft_clips_3_end = BTreeMap::new();

        let mut svfdup = Vec::new();
        let mut svrdup = Vec::new();

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
        );

        // Verification: Softclip should be used
        assert!(soft_clips_5_end.get(&50).unwrap().used);
        assert!(
            !insertion_variants.is_empty(),
            "Insertion variants must be created"
        );
    }
}
