use crate::variations::{Sclip, Variation};
use std::collections::BTreeMap;

/// Calculate strand bias flag.
///
/// Ports `VariationUtils.strandBias`.
pub fn strand_bias(
    forward_count: i32,
    reverse_count: i32,
    bias_cutoff: f64,
    min_bias_reads: i32,
) -> i32 {
    // using p=0.01, because prop.test(1,12) = 0.01
    if forward_count + reverse_count <= 12 {
        return if forward_count * reverse_count > 0 {
            2
        } else {
            0
        };
    }

    let total = (forward_count + reverse_count) as f64;
    let forward_ratio = (forward_count as f64) / total;
    let reverse_ratio = (reverse_count as f64) / total;

    if forward_ratio >= bias_cutoff
        && reverse_ratio >= bias_cutoff
        && forward_count >= min_bias_reads
        && reverse_count >= min_bias_reads
    {
        2
    } else {
        1
    }
}

/// Check if the sequence is low complexity.
/// Low-complexity sequences are simple repeats or regions highly enriched for just one letter.
///
/// Ports `VariationRealigner.islowcomplexseq`.
pub fn is_low_complex_seq(seq: &[u8]) -> bool {
    let len = seq.len();
    if len == 0 {
        return true;
    }
    let mut nt_cnt = 0;

    for base in [b'A', b'T', b'G', b'C'] {
        let cnt = seq
            .iter()
            .filter(|&&b| b == base || b == base.to_ascii_lowercase())
            .count();
        if cnt > 0 {
            nt_cnt += 1;
        }
        if (cnt as f64) / (len as f64) > 0.75 {
            return true;
        }
    }

    nt_cnt < 3
}

/// Find the consensus sequence in soft-clipped reads.
/// Consensus is called if the matched nucleotides are >90% of all softly clipped nucleotides.
///
/// Ports `VariationUtils.findconseq`.
pub fn find_conseq(
    soft_clip: &mut Sclip,
    dir: i32,
    adaptor_forward: &BTreeMap<Vec<u8>, i32>,
    adaptor_reverse: &BTreeMap<Vec<u8>, i32>,
) -> Vec<u8> {
    if let Some(ref seq) = soft_clip.sequence {
        return seq.clone();
    }

    let mut total = 0;
    let mut match_cnt = 0;
    let mut seq = Vec::new();
    let mut flag = false;

    for (&position_in_sclip, nv) in &soft_clip.nt {
        let mut max_count = 0;
        let mut max_quality = 0.0;
        let mut chosen_base: Option<u8> = None;
        let mut total_count = 0;

        for (&current_base, &current_count) in nv {
            total_count += current_count;

            let has_seq_qual = soft_clip
                .seq
                .get(&position_in_sclip)
                .and_then(|m| m.get(&current_base));
            let mean_qual = has_seq_qual.map(|v| v.qmean).unwrap_or(0.0);

            if current_count > max_count || (has_seq_qual.is_some() && mean_qual > max_quality) {
                max_count = current_count;
                chosen_base = Some(current_base);
                max_quality = mean_qual;
            }
        }

        let vars_count = soft_clip.variation.cnt;
        if position_in_sclip == 3
            && soft_clip.nt.len() >= 6
            && vars_count > 0
            && (total_count as f64) / (vars_count as f64) < 0.2
            && total_count <= 2
        {
            break;
        }

        if (total_count - max_count > 2 || max_count <= total_count - max_count)
            && total_count > 0
            && (max_count as f64) / (total_count as f64) < 0.8
        {
            if flag {
                break;
            }
            flag = true;
        }

        total += total_count;
        match_cnt += max_count;
        if let Some(b) = chosen_base {
            seq.push(b);
        }
    }

    let nt_size = soft_clip.nt.len() as i32;
    let seq_len = seq.len() as i32;
    let mut final_seq = Vec::new();

    if total != 0
        && (match_cnt as f64) / (total as f64) > 0.9
        && (seq_len as f64) / 1.5 > (nt_size - seq_len) as f64
        && ((seq_len as f64) / (nt_size as f64) > 0.8 || nt_size - seq_len < 10 || seq_len > 25)
    {
        final_seq = seq;
    }

    if !final_seq.is_empty() && final_seq.len() > 12 {
        let matches_a7 = final_seq
            .get(1..8)
            .is_some_and(|s| s.iter().all(|&b| b == b'A' || b == b'a'));
        let matches_t7 = final_seq
            .get(1..8)
            .is_some_and(|s| s.iter().all(|&b| b == b'T' || b == b't'));
        if matches_a7 || matches_t7 || is_low_complex_seq(&final_seq) {
            soft_clip.used = true;
        }
    }

    let mut cleared_by_adaptor = false;
    if !final_seq.is_empty() && final_seq.len() >= 6 {
        if dir == 3 {
            if let Some(seed) = final_seq.get(0..6) {
                if adaptor_forward.contains_key(seed) {
                    final_seq.clear();
                    cleared_by_adaptor = true;
                }
            }
        } else if dir == 5 {
            if let Some(seed) = final_seq.get(0..6) {
                let mut rev_seed = seed.to_vec();
                rev_seed.reverse();
                if adaptor_reverse.contains_key(&rev_seed) {
                    final_seq.clear();
                    cleared_by_adaptor = true;
                }
            }
        }
    }

    if !cleared_by_adaptor {
        soft_clip.sequence = Some(final_seq.clone());
    }
    final_seq
}

/// Correct counts for negative values.
///
/// Ports `VariationUtils.correctCnt`.
pub fn correct_cnt(var_to_correct: &mut Variation) {
    if var_to_correct.cnt < 0 {
        var_to_correct.cnt = 0;
    }
    if var_to_correct.hicnt < 0 {
        var_to_correct.hicnt = 0;
    }
    if var_to_correct.locnt < 0 {
        var_to_correct.locnt = 0;
    }
    if var_to_correct.pmean < 0.0 {
        var_to_correct.pmean = 0.0;
    }
    if var_to_correct.qmean < 0.0 {
        var_to_correct.qmean = 0.0;
    }
    if var_to_correct.mqmean < 0.0 {
        var_to_correct.mqmean = 0.0;
    }
    if var_to_correct.get_dir(true) < 0 {
        let val = var_to_correct.get_dir(true);
        var_to_correct.add_dir(true, -val);
    }
    if var_to_correct.get_dir(false) < 0 {
        let val = var_to_correct.get_dir(false);
        var_to_correct.add_dir(false, -val);
    }
}

/// Adjust the count.
/// If `reference_var` is not None, the count for reference is also adjusted.
/// Variant counts of `variant` are added to `var_to_add` and removed from `reference_var`.
///
/// Ports `VariationUtils.adjCnt`.
pub fn adj_cnt(
    var_to_add: &mut Variation,
    variant: &Variation,
    reference_var: Option<&mut Variation>,
) {
    var_to_add.cnt += variant.cnt;
    var_to_add.extracnt += variant.cnt;
    var_to_add.hicnt += variant.hicnt;
    var_to_add.locnt += variant.locnt;
    var_to_add.pmean += variant.pmean;
    var_to_add.qmean += variant.qmean;
    var_to_add.mqmean += variant.mqmean;
    var_to_add.nm += variant.nm;
    var_to_add.pstd = true;
    var_to_add.qstd = true;
    var_to_add.add_dir(true, variant.get_dir(true));
    var_to_add.add_dir(false, variant.get_dir(false));

    if let Some(ref_var) = reference_var {
        ref_var.cnt -= variant.cnt;
        ref_var.hicnt -= variant.hicnt;
        ref_var.locnt -= variant.locnt;
        ref_var.pmean -= variant.pmean;
        ref_var.qmean -= variant.qmean;
        ref_var.mqmean -= variant.mqmean;
        ref_var.nm -= variant.nm;
        ref_var.sub_dir(true, variant.get_dir(true));
        ref_var.sub_dir(false, variant.get_dir(false));
        correct_cnt(ref_var);
    }
}

/// Adjust a variation by a given factor.
///
/// Ports `VariationRealigner.adjRefFactor`.
pub fn adj_ref_factor(ref_var: &mut Variation, mut factor_f: f64) {
    if factor_f > 1.0 {
        factor_f = 1.0;
    }
    if factor_f < -1.0 {
        return;
    }

    let old_vars_count = ref_var.cnt;
    ref_var.cnt -= (factor_f * (old_vars_count as f64)) as i32;
    ref_var.hicnt -= (factor_f * (ref_var.hicnt as f64)) as i32;
    ref_var.locnt -= (factor_f * (ref_var.locnt as f64)) as i32;

    let mut factor_cnt = if old_vars_count != 0 {
        ((ref_var.cnt - old_vars_count).abs() as f64) / (old_vars_count as f64)
    } else {
        1.0
    };

    if (factor_f < 0.0 && factor_cnt > 0.0) || (factor_f > 0.0 && factor_cnt < 0.0) {
        factor_cnt = -factor_cnt;
    }

    ref_var.pmean -= ref_var.pmean * factor_cnt;
    ref_var.qmean -= ref_var.qmean * factor_cnt;
    ref_var.mqmean -= ref_var.mqmean * factor_cnt;

    ref_var.nm -= factor_f * ref_var.nm;
    ref_var.sub_dir(true, (factor_f * (ref_var.get_dir(true) as f64)) as i32);
    ref_var.sub_dir(false, (factor_f * (ref_var.get_dir(false) as f64)) as i32);

    correct_cnt(ref_var);
}

/// Adjust reference count by factor.
///
/// Ports `VariationRealigner.adjRefCnt`.
pub fn adj_ref_cnt(tv: &Variation, ref_var: &mut Variation, len: i32) {
    let f = if tv.pmean != 0.0 && tv.cnt != 0 {
        ((tv.pmean / (tv.cnt as f64)) - (len as f64) + 1.0) / (tv.pmean / (tv.cnt as f64))
    } else {
        0.0
    };
    let mut f = f;
    if f < 0.0 {
        return;
    }
    if f > 1.0 {
        f = 1.0;
    }

    ref_var.cnt -= (f * (tv.cnt as f64)) as i32;
    ref_var.hicnt -= (f * (tv.hicnt as f64)) as i32;
    ref_var.locnt -= (f * (tv.locnt as f64)) as i32;
    ref_var.pmean -= f * tv.pmean;
    ref_var.qmean -= f * tv.qmean;
    ref_var.mqmean -= f * tv.mqmean;
    ref_var.nm -= f * tv.nm;
    ref_var.sub_dir(true, (f * (tv.get_dir(true) as f64)) as i32);
    ref_var.sub_dir(false, (f * (tv.get_dir(false) as f64)) as i32);
    correct_cnt(ref_var);
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::realign::{find_35_match, Match35};
    use crate::variations::{Sclip, Variation};
    use std::collections::BTreeMap;

    fn find_35_match_legacy(seq5: &[u8], seq3: &[u8], verbose: bool) -> Match35 {
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

    pub(crate) fn run_find_35_match_diff_suite() {
        // 1. Short boundary cases (< 8, == 8, == 9)
        for l5 in 0..=10 {
            for l3 in 0..=10 {
                let s5 = vec![b'A'; l5];
                let s3 = vec![b'A'; l3];
                assert_eq!(
                    find_35_match(&s5, &s3, false),
                    find_35_match_legacy(&s5, &s3, false),
                    "short boundary mismatch at l5={l5}, l3={l3}"
                );
            }
        }

        // 2. Adversarial boundary cases for rem5 < rem3, rem5 == rem3, rem3 < rem5
        //    specifically exercising the sentinel iteration at i + total_length == len5
        //    with 0, 1, 2, 3, 4 mismatches before the sentinel.
        let base_motif = b"ACGTACGTACGTACGTACGTACGTACGT";
        for rem5_len in 9..=22 {
            for rem3_len in 9..=22 {
                for mm_count in 0..=4 {
                    let s5_prefix = &base_motif[..rem5_len];
                    let mut s5 = s5_prefix.to_vec();
                    let s5_len = s5.len();
                    for k in 0..mm_count.min(s5_len) {
                        s5[k * 2 % s5_len] = b'N';
                    }
                    // Build s3 so that reversed prefix of length rem3_len aligns with s5_prefix
                    let mut s3_rev = vec![b'T'; rem3_len];
                    let copy_len = rem5_len.min(rem3_len);
                    s3_rev[..copy_len].copy_from_slice(&base_motif[..copy_len]);
                    let s3: Vec<u8> = s3_rev.into_iter().rev().collect();
                    assert_eq!(
                        find_35_match(&s5, &s3, false),
                        find_35_match_legacy(&s5, &s3, false),
                        "adversarial mismatch at rem5={rem5_len}, rem3={rem3_len}, mm={mm_count}"
                    );
                }
            }
        }

        // 3. 50,000 seeded random pairs (random & synthetic reverse overlaps with 0..=4 mismatches)
        const ALPHABET: &[u8; 5] = b"ACGTN";
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next_u64 = || -> u64 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        for case_idx in 0..50_000 {
            let len5 = (next_u64() as usize) % 60;
            let len3 = (next_u64() as usize) % 60;
            let mut s5: Vec<u8> = (0..len5)
                .map(|_| ALPHABET[(next_u64() as usize) % 5])
                .collect();
            let mut s3: Vec<u8> = (0..len3)
                .map(|_| ALPHABET[(next_u64() as usize) % 5])
                .collect();

            if case_idx % 2 == 1 && len5 >= 10 && len3 >= 10 {
                // Inject a reverse-matching segment with 0..=4 mismatches
                let max_ov = (len5 - 1).min(len3 - 1);
                let ov = 8 + ((next_u64() as usize) % (max_ov - 7).max(1));
                let i_off = (next_u64() as usize) % (len5 - ov).max(1);
                let j_off = 1 + ((next_u64() as usize) % (len3 - ov).max(1));
                for t in 0..ov {
                    if i_off + t < len5 && j_off + t <= len3 {
                        s3[len3 - j_off - t] = s5[i_off + t];
                    }
                }
                let mm = (next_u64() as usize) % 5;
                for _ in 0..mm {
                    let t = (next_u64() as usize) % ov;
                    if i_off + t < len5 {
                        s5[i_off + t] ^= 0x01;
                    }
                }
            }

            let actual = find_35_match(&s5, &s3, false);
            let expected = find_35_match_legacy(&s5, &s3, false);
            assert_eq!(actual, expected, "random case {case_idx} failed");
        }

        // 4. Stack-to-heap spill boundary cases (510..=600) and high-bit byte values (0x80..=0xFF)
        for &l3 in &[510usize, 512, 513, 521, 522, 600] {
            for &l5 in &[16usize, 512, 522] {
                let mut s5: Vec<u8> = (0..l5)
                    .map(|idx| 0x80u8.wrapping_add((idx as u8).wrapping_mul(31)))
                    .collect();
                let mut s3: Vec<u8> = (0..l3)
                    .map(|idx| 0x80u8.wrapping_add((idx as u8).wrapping_mul(17)))
                    .collect();
                // Plant a 14-bp reverse match with 1 high-bit mismatch near the tail
                for t in 0..14 {
                    s3[t] = s5[l5 - 14 + t];
                }
                s5[l5 - 7] ^= 0x81;
                assert_eq!(
                    find_35_match(&s5, &s3, false),
                    find_35_match_legacy(&s5, &s3, false),
                    "spill/high-bit mismatch at l5={l5}, l3={l3}"
                );
            }
        }
    }

    #[test]
    fn find_35_match_matches_legacy_on_random_and_adversarial_corpus() {
        run_find_35_match_diff_suite();
    }

    #[test]
    fn realignlgins30_hoisted_distance_and_sc5_used_match_behavior() {
        let mut reference = crate::reference::Reference::new();
        reference
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 5000]);
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(100),
            end: crate::types::Pos1(4000),
            gene: String::new(),
            amplicon: None,
        };
        let config = crate::config::Config::default();
        let mut non_ins = crate::javamap::JavaHashMap::new();
        let mut ins = crate::javamap::JavaHashMap::new();
        let mut ref_cov = rustc_hash::FxHashMap::default();
        let mut sclip5 = BTreeMap::new();
        let mut sclip3 = BTreeMap::new();

        // Populate a 5' softclip at p5 = 1000 whose consensus is >12bp of 'A' (low complexity),
        // which causes `find_conseq` to mark `sc5.used = true` on the first in-range `p3`.
        let mut sc5 = Sclip::new();
        sc5.variation.cnt = 5;
        for idx in 0..15 {
            let mut nt = BTreeMap::new();
            nt.insert(b'A', 10);
            sc5.nt.insert(idx, nt);
        }
        sclip5.insert(1000, sc5);

        // Add an out-of-range 3' softclip (p3 = 200, distance 800 > 150 * 2.5 = 375) with higher cnt
        // so it sorts first in `tmp3`; distance hoisting must skip it without touching `sc5` or `sc3`.
        let mut sc3_far = Sclip::new();
        sc3_far.variation.cnt = 10;
        for idx in 0..15 {
            let mut nt = BTreeMap::new();
            nt.insert(b'T', 10);
            sc3_far.nt.insert(idx, nt);
        }
        sclip3.insert(200, sc3_far);

        // Add two in-range 3' softclips (p3 = 1010 and p3 = 1020).
        // The first in-range p3 (1010) triggers `find_conseq` on `sc5` (setting `sc5.used = true`)
        // and on `sc3[1010]`, then continues; the second in-range p3 (1020) must immediately hit
        // `if sc5_used { break; }` and leave `sc3[1020].sequence` as `None`.
        for (pos, cnt) in [(1010, 8), (1020, 4)] {
            let mut sc3 = Sclip::new();
            sc3.variation.cnt = cnt;
            for idx in 0..15 {
                let mut nt = BTreeMap::new();
                nt.insert(b'C', 10);
                sc3.nt.insert(idx, nt);
            }
            sclip3.insert(pos, sc3);
        }

        crate::realign::realignlgins30(
            None,
            &mut reference,
            "chr1",
            5000,
            &region,
            &config,
            &mut non_ins,
            &mut ins,
            &mut ref_cov,
            &mut sclip5,
            &mut sclip3,
            150,
        );

        assert!(sclip5.get(&1000).unwrap().used);
        assert!(
            sclip3.get(&200).unwrap().sequence.is_none(),
            "far p3=200 should be skipped by distance filter before find_conseq"
        );
        assert!(
            sclip3.get(&1010).unwrap().sequence.is_some(),
            "first in-range p3=1010 should run find_conseq"
        );
        assert!(
            sclip3.get(&1020).unwrap().sequence.is_none(),
            "second in-range p3=1020 should be skipped because sc5_used broke out of the loop"
        );
    }

    #[test]
    fn test_findconseq() {
        let mut sc = Sclip::new();
        // nv map at position 0: 'A' -> 10, 'T' -> 1
        let mut nt0 = BTreeMap::new();
        nt0.insert(b'A', 10);
        nt0.insert(b'T', 1);
        sc.nt.insert(0, nt0);

        // nv map at position 1: 'A' -> 12
        let mut nt1 = BTreeMap::new();
        nt1.insert(b'A', 12);
        sc.nt.insert(1, nt1);

        // nv map at position 2: 'G' -> 11
        let mut nt2 = BTreeMap::new();
        nt2.insert(b'G', 11);
        sc.nt.insert(2, nt2);

        let adaptor_forward = BTreeMap::new();
        let adaptor_reverse = BTreeMap::new();

        let consensus = find_conseq(&mut sc, 3, &adaptor_forward, &adaptor_reverse);
        // total: 11 + 12 + 11 = 34
        // match: 10 + 12 + 11 = 33
        // match / total = 33 / 34 = 0.97 > 0.9
        // seq len = 3. nt_size = 3.
        // seq_len / 1.5 = 2. nt_size - seq_len = 0.
        // 2 > 0.
        // seq_len / nt_size = 1.0 > 0.8
        // So consensus is accepted: "AAG"
        assert_eq!(consensus, b"AAG");
        assert_eq!(sc.sequence, Some(b"AAG".to_vec()));
    }

    #[test]
    fn test_findconseq_rejected() {
        let mut sc = Sclip::new();
        let mut nt0 = BTreeMap::new();
        nt0.insert(b'A', 5);
        nt0.insert(b'T', 5);
        sc.nt.insert(0, nt0);

        let adaptor_forward = BTreeMap::new();
        let adaptor_reverse = BTreeMap::new();

        let consensus = find_conseq(&mut sc, 3, &adaptor_forward, &adaptor_reverse);
        assert_eq!(consensus, b"");
    }

    #[test]
    fn find_conseq_caches_empty_consensus_without_recomputing() {
        let mut sc = Sclip::new();
        let mut nt0 = BTreeMap::new();
        nt0.insert(b'A', 5);
        nt0.insert(b'T', 5);
        sc.nt.insert(0, nt0);

        let adaptor_forward = BTreeMap::new();
        let adaptor_reverse = BTreeMap::new();

        let first = find_conseq(&mut sc, 3, &adaptor_forward, &adaptor_reverse);
        assert!(first.is_empty());
        assert_eq!(sc.sequence, Some(Vec::new()));

        // Populate a valid consensus into `sc.nt` without clearing `sc.sequence`.
        // Because `sc.sequence` is already `Some(vec![])`, a second `find_conseq`
        // call must return the cached empty consensus without recomputing from `sc.nt`.
        sc.nt.clear();
        for pos in 0..3 {
            let mut nt = BTreeMap::new();
            nt.insert(b'G', 12);
            sc.nt.insert(pos, nt);
        }
        let second = find_conseq(&mut sc, 3, &adaptor_forward, &adaptor_reverse);
        assert!(
            second.is_empty(),
            "expected cached empty consensus Some(vec![]) to be returned without recomputing"
        );

        // However, if a valid consensus (>=6bp) was cleared solely by direction-specific
        // adaptor matching (`dir=3` or `dir=5`), `sc.sequence` must remain `None` so a
        // subsequent `dir=0` call recomputes and returns the unfiltered consensus.
        let mut sc_adapt = Sclip::new();
        for pos in 0..8 {
            let mut nt = BTreeMap::new();
            nt.insert(b'C', 12);
            sc_adapt.nt.insert(pos, nt);
        }
        let mut fwd_adapt = BTreeMap::new();
        fwd_adapt.insert(b"CCCCCC".to_vec(), 1);
        let dir3 = find_conseq(&mut sc_adapt, 3, &fwd_adapt, &adaptor_reverse);
        assert!(
            dir3.is_empty(),
            "dir=3 should be cleared by forward adaptor"
        );
        assert!(
            sc_adapt.sequence.is_none(),
            "adaptor-cleared consensus must not lock sc.sequence to empty"
        );
        let dir0 = find_conseq(&mut sc_adapt, 0, &fwd_adapt, &adaptor_reverse);
        assert_eq!(
            dir0, b"CCCCCCCC",
            "subsequent dir=0 call must recompute unfiltered consensus"
        );
    }

    #[test]
    fn test_is_low_complex_seq() {
        assert!(is_low_complex_seq(b"AAAAAAAAAAAAA"));
        assert!(is_low_complex_seq(b"A"));
        assert!(is_low_complex_seq(b"AT"));
        assert!(!is_low_complex_seq(b"ATGC"));
    }

    #[test]
    fn test_strand_bias() {
        assert_eq!(strand_bias(5, 5, 0.05, 2), 2);
        assert_eq!(strand_bias(5, 0, 0.05, 2), 0);
        assert_eq!(strand_bias(10, 10, 0.05, 2), 2);
        assert_eq!(strand_bias(20, 0, 0.05, 2), 1);
    }

    #[test]
    fn test_correct_cnt() {
        let mut v = Variation::new();
        v.cnt = -5;
        v.pmean = -1.0;
        v.add_dir(true, -3);
        correct_cnt(&mut v);
        assert_eq!(v.cnt, 0);
        assert_eq!(v.pmean, 0.0);
        assert_eq!(v.get_dir(true), 0);
    }

    #[test]
    fn test_adj_cnt() {
        let mut v1 = Variation::new();
        let mut v2 = Variation::new();
        v2.cnt = 5;
        v2.hicnt = 3;
        v2.add_dir(true, 4);

        adj_cnt(&mut v1, &v2, None);
        assert_eq!(v1.cnt, 5);
        assert_eq!(v1.hicnt, 3);
        assert_eq!(v1.get_dir(true), 4);
    }

    #[test]
    fn test_adj_ref_cnt() {
        let mut tv = Variation::new();
        tv.cnt = 10;
        tv.hicnt = 10;
        tv.pmean = 100.0; // mean pos is 100.0 / 10 = 10
        tv.add_dir(true, 10);

        let mut ref_var = Variation::new();
        ref_var.cnt = 20;
        ref_var.hicnt = 20;
        ref_var.pmean = 200.0;
        ref_var.add_dir(true, 20);

        // len = 5
        // tv.pmean / tv.cnt = 10
        // f = (10 - 5 + 1) / 10 = 6 / 10 = 0.6
        // ref.cnt -= 0.6 * 10 = 6 => 14
        adj_ref_cnt(&tv, &mut ref_var, 5);
        assert_eq!(ref_var.cnt, 14);
        assert_eq!(ref_var.hicnt, 14);
    }
}
