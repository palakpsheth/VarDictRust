//! CIGAR modify — ports CigarModifier.java (ordered fixpoint, bounds-checked)
//!
//! normalizes indels near soft-clips, combines adjacent small indels into
//! complex variants, and cleans up end-of-read alignment mismatches.

use crate::cigar::ModifiedCigar;
use crate::config::Config;
use crate::util::{char_at, reverse_complement, subseq};
use regex::Regex;
use std::collections::{HashMap, HashSet};
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

static ANY_NUMBER_M_NUMBER_S_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^(.*?)(\d+)M(\d+)S$").unwrap());
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
static CIGAR_ELEMENT_PATTERN: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\-?\d+)[MIDNSHPX=]").unwrap());

const LOWQUAL: u8 = 10;
const SEED_2: usize = 12;

/// Helper to parse matched group as integer with default 0.
fn to_int(s: &str) -> i32 {
    s.parse::<i32>().unwrap_or(0)
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
    seed: &'a HashMap<Vec<u8>, Vec<u32>>,
    indel: i32,
    max_read_length: i32,
    config: &'a Config,
}

impl<'a> CigarModifier<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        position: i32,
        cigar_str: &str,
        query_seq: &[u8],
        query_qual: &[u8],
        reference: &'a [u8],
        ref_start: i32,
        seed: &'a HashMap<Vec<u8>, Vec<u32>>,
        indel: i32,
        max_read_length: i32,
        config: &'a Config,
    ) -> Self {
        CigarModifier {
            position,
            cigar_str: cigar_str.to_string(),
            query_seq: query_seq.to_vec(),
            query_qual: query_qual.to_vec(),
            reference,
            ref_start,
            seed,
            indel,
            max_read_length,
            config,
        }
    }

    /// Perform all CIGAR modifications. No panic allowed!
    pub fn modify_cigar(&mut self) -> Result<ModifiedCigar, crate::types::VardictError> {
        //flag is set to true if CIGAR string is modified and should be looked at again
        let flag = true;
        // Wrap implementation in safe checks matching Java try/catch.
        // If anything fails, we fall back gracefully to the original state.
        let orig_pos = self.position;
        let orig_cigar = self.cigar_str.clone();
        let orig_seq = self.query_seq.clone();
        let orig_qual = self.query_qual.clone();

        let run_modification = move |this: &mut Self| -> Result<(), crate::types::VardictError> {
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
                                            this.query_seq =
                                                this.query_seq[cigar_element..].to_vec();
                                            this.query_qual =
                                                this.query_qual[cigar_element..].to_vec();
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
                                        if (this.position - ref_pos).abs()
                                            < 2 * this.max_read_length
                                        {
                                            this.cigar_str =
                                                RE_END_2S.replace(&this.cigar_str, "").into_owned();
                                            this.query_seq = this.query_seq[0..start_idx].to_vec();
                                            this.query_qual =
                                                this.query_qual[0..start_idx].to_vec();
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
                    flag_loop = this.combine_to_close_to_correct(
                        g1,
                        g2,
                        g3,
                        &op,
                        istr.as_deref(),
                        flag_loop,
                    );
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
                    let dlen = to_int(caps.get(1).unwrap().as_str())
                        + to_int(caps.get(2).unwrap().as_str());
                    this.cigar_str = DIG_D_DIG_D
                        .replace(&this.cigar_str, &format!("{}D", dlen))
                        .into_owned();
                    flag_loop = true;
                }

                if let Some(caps) = DIG_I_DIG_I.captures(&this.cigar_str) {
                    let ilen = to_int(caps.get(1).unwrap().as_str())
                        + to_int(caps.get(2).unwrap().as_str());
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

            Ok(())
        };

        if let Err(e) = run_modification(self) {
            // Restore original state
            self.position = orig_pos;
            self.cigar_str = orig_cigar;
            self.query_seq = orig_seq;
            self.query_qual = orig_qual;
            log::warn!("CigarModifier exception (restored original state): {}", e);
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
            cigar: self.cigar_str.clone(),
            query_seq: self.query_seq.clone(),
            query_qual: self.query_qual.clone(),
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
    use std::collections::HashMap;

    fn get_default_config() -> Config {
        Config::default()
    }

    #[test]
    fn strip_del() {
        let config = get_default_config();
        let seed = HashMap::new();
        let reference = b"";
        // 1. Starts with D: should adjust position and strip leading D.
        let mut cm = CigarModifier::new(
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
        let seed = HashMap::new();
        let reference = b"";
        // 2. Starts or ends with I: should replace with S.
        let mut cm = CigarModifier::new(
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
        let seed = HashMap::new();
        let reference = b"";
        // Starts with softclip followed by insertion or deletion: should combine them.
        let mut cm = CigarModifier::new(
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

        let mut cm = CigarModifier::new(
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
        let seed = HashMap::new();
        let reference = b"";
        // Combine adjacent deletions
        let mut cm = CigarModifier::new(
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
        let mut cm = CigarModifier::new(
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
        let seed = HashMap::new();
        let reference = b"";
        // Feed various weird CIGAR strings that used to crash VarDictJava due to bounds check errors
        let inputs = vec!["0M", "100D", "100I", "5M100D5M", "5M100I5M", "1M1I1D1M"];
        for cig in inputs {
            let mut cm = CigarModifier::new(
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
        let seed = HashMap::new();
        let query_seq = b"AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTACAAAAAAAAAAAAAAAAAAAAGAAAAAAGAAAAAAAAAG";
        let query_qual = b"BBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF=74??44444444FFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        let ref_seq = b"AATCCTCTGTTTGGAAACCAGCCCGATTCAAGGAGGGTTCTGATGCACTGCGGTGAATTTTTGGCAATGAGCGGGCCAGCAGCTCAATAGAGGCGAGAATCTACAAAAAAAAAAAGAAAAAAAAAAGAAAAAAAAAGAAAAAAGAAAAAAAAAGAAAGAAAGAAAAAGAAAAAACAGAAAGAAGAATGAAAATCTGGGTG";

        let mut cm = CigarModifier::new(
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
        let seed = HashMap::new();
        let ref_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let query_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let query_qual = vec![b'I'; query_seq.len()];

        // Adversarial inputs for three_deletions, two_deletions_insertion_to_complex, three_indels
        // with short M components that cause non-positive operator lengths.
        let adversarial_cases = ["10M1D1M1D1M1D1M", "10M1I1M1D1M1D1M", "10M1D1M1I1M1D1M"];

        for cigar in adversarial_cases {
            let mut cm = CigarModifier::new(
                1,
                cigar,
                query_seq,
                &query_qual,
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
}
