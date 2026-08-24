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
        if !seq.is_empty() {
            return seq.clone();
        }
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
        let matches_a7 =
            final_seq.len() >= 8 && final_seq[1..8].iter().all(|&b| b == b'A' || b == b'a');
        let matches_t7 =
            final_seq.len() >= 8 && final_seq[1..8].iter().all(|&b| b == b'T' || b == b't');
        if matches_a7 || matches_t7 || is_low_complex_seq(&final_seq) {
            soft_clip.used = true;
        }
    }

    if !final_seq.is_empty() && final_seq.len() >= 6 {
        if dir == 3 {
            let seed = &final_seq[0..6];
            if adaptor_forward.contains_key(seed) {
                final_seq.clear();
            }
        } else if dir == 5 {
            let seed = &final_seq[0..6];
            let mut rev_seed = seed.to_vec();
            rev_seed.reverse();
            if adaptor_reverse.contains_key(&rev_seed) {
                final_seq.clear();
            }
        }
    }

    soft_clip.sequence = Some(final_seq.clone());
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
mod tests {
    use super::*;
    use crate::variations::{Sclip, Variation};
    use std::collections::BTreeMap;

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
