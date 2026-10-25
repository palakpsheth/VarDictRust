//! Structural Variant (SV) helpers and calling algorithms.
//!
//! Ports `StructuralVariantsProcessor.java`.
//!
//! Module-level allows kept post-W5.4 (blanket list audited & shrunk):
//! - `unused_variables`: parameters kept for Java-signature parity.
//! - `non_snake_case`: Java-verbatim identifiers (`rQmean`, …) kept for port traceability.
//! - `unused_mut`: bindings stay `mut` where Java mutates them even when this port's flow doesn't.
//! - `unused_assignments`: Java assigns before overwrite; dead writes kept for structural parity.
//! - `clippy::neg_multiply`: `-1 *` forms mirror the Java expressions.
//! - `clippy::map_entry`: contains_key+insert mirrors CigarParser.java flow.
#![allow(
    unused_variables,
    non_snake_case,
    unused_mut,
    clippy::neg_multiply,
    clippy::map_entry,
    unused_assignments
)]

use crate::config::Config;
use crate::reference::{get_reference, is_loaded, Reference};
use crate::types::{Pos1, VardictError};
use crate::variations::{Cluster, Mate, Sclip, Variation};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Check cluster for start and end of SV mates.
/// Ports `VariationRealigner.checkCluster`.
pub fn check_cluster(mates: &mut [Mate], rlen: i32) -> Cluster {
    if mates.is_empty() {
        return Cluster::default();
    }
    mates.sort_by_key(|mate| mate.mate_start);

    let mut clusters: Vec<Cluster> = Vec::new();
    let first_mate = &mates[0];
    clusters.push(Cluster::new(
        0,
        Mate {
            mate_start: first_mate.mate_start,
            mate_end: first_mate.mate_end,
            mate_len: 0,
            start: first_mate.start,
            end: first_mate.end,
            pmean: 0.0,
            qmean: 0.0,
            mqmean: 0.0,
            nm: 0.0,
        },
    ));

    let mut cur = 0;
    let min_svc_dist = 1.5;
    let limit = (min_svc_dist * rlen as f64) as i32;

    for mate_m in mates.iter() {
        let current_cluster_end = clusters[cur].mate.mate_end;
        if mate_m.mate_start - current_cluster_end > limit {
            cur += 1;
            clusters.push(Cluster::new(
                0,
                Mate {
                    mate_start: mate_m.mate_start,
                    mate_end: mate_m.mate_end,
                    mate_len: 0,
                    start: mate_m.start,
                    end: mate_m.end,
                    pmean: 0.0,
                    qmean: 0.0,
                    mqmean: 0.0,
                    nm: 0.0,
                },
            ));
        }

        let current_cluster = &mut clusters[cur];
        current_cluster.cnt += 1;
        current_cluster.mate.mate_len += mate_m.mate_len;

        if mate_m.mate_end > current_cluster.mate.mate_end {
            current_cluster.mate.mate_end = mate_m.mate_end;
        }
        if mate_m.start < current_cluster.mate.start {
            current_cluster.mate.start = mate_m.start;
        }
        if mate_m.end > current_cluster.mate.end {
            current_cluster.mate.end = mate_m.end;
        }
        current_cluster.mate.pmean += mate_m.pmean;
        current_cluster.mate.qmean += mate_m.qmean;
        current_cluster.mate.mqmean += mate_m.mqmean;
        current_cluster.mate.nm += mate_m.nm;
    }

    clusters.sort_by_key(|c| std::cmp::Reverse(c.cnt));

    let first_cluster = &clusters[0];
    if (first_cluster.cnt as f64) / (mates.len() as f64) >= 0.60 {
        Cluster::new(
            first_cluster.cnt,
            Mate {
                mate_start: first_cluster.mate.mate_start,
                mate_end: first_cluster.mate.mate_end,
                mate_len: first_cluster.mate.mate_len / first_cluster.cnt,
                start: first_cluster.mate.start,
                end: first_cluster.mate.end,
                pmean: first_cluster.mate.pmean,
                qmean: first_cluster.mate.qmean,
                mqmean: first_cluster.mate.mqmean,
                nm: first_cluster.mate.nm,
            },
        )
    } else {
        Cluster::default()
    }
}

/// Filter possible SVs by checking created clusters, updates data of possible SV from cluster data.
/// Ports `VariationRealigner.filterSV`.
pub fn filter_sv(sv_list: &mut [Sclip], max_read_length: i32) {
    for sv in sv_list.iter_mut() {
        if sv.mates.is_empty() {
            continue;
        }
        let cluster = check_cluster(&mut sv.mates, max_read_length);

        if cluster.mate.mate_start != 0 {
            sv.mstart = cluster.mate.mate_start;
            sv.mend = cluster.mate.mate_end;
            sv.variation.cnt = cluster.cnt;
            sv.mlen = cluster.mate.mate_len;
            sv.start = cluster.mate.start;
            sv.end = cluster.mate.end;
            sv.variation.pmean = cluster.mate.pmean;
            sv.variation.qmean = cluster.mate.qmean;
            sv.variation.mqmean = cluster.mate.mqmean;
            sv.variation.nm = cluster.mate.nm;
        } else {
            sv.used = true;
        }

        // Too many unhappy mates are false positive
        if sv.disc != 0 && (sv.variation.cnt as f64) / (sv.disc as f64) < 0.5 {
            if !((sv.variation.cnt as f64) / (sv.disc as f64) >= 0.35 && sv.variation.cnt >= 5) {
                sv.used = true;
            }
        }

        let mut soft_sorted: Vec<(i32, i32)> = sv.soft.iter().map(|(&k, &v)| (k, v)).collect();
        soft_sorted.sort_by_key(|&(_, count)| std::cmp::Reverse(count));

        sv.softp = if soft_sorted.is_empty() {
            0
        } else {
            soft_sorted[0].0
        };
    }
}

/// Filter all possible structural variants structures.
/// Ports `VariationRealigner.filterAllSVStructures`.
pub fn filter_all_sv_structures(var_data: &mut crate::variations::VariationData) {
    if var_data.sv_structures_filtered {
        return;
    }
    let max_read_length = var_data.max_read_length;

    filter_sv(&mut var_data.sv_structures.svfinv3, max_read_length);
    filter_sv(&mut var_data.sv_structures.svrinv3, max_read_length);
    filter_sv(&mut var_data.sv_structures.svfinv5, max_read_length);
    filter_sv(&mut var_data.sv_structures.svrinv5, max_read_length);
    filter_sv(&mut var_data.sv_structures.svfdel, max_read_length);
    filter_sv(&mut var_data.sv_structures.svrdel, max_read_length);
    filter_sv(&mut var_data.sv_structures.svfdup, max_read_length);
    filter_sv(&mut var_data.sv_structures.svrdup, max_read_length);

    for sv_list in var_data.sv_structures.svffus.values_mut() {
        filter_sv(sv_list, max_read_length);
    }
    for sv_list in var_data.sv_structures.svrfus.values_mut() {
        filter_sv(sv_list, max_read_length);
    }

    var_data.sv_structures_filtered = true;
}

/// Check if the SV candidate with the highest count at `softp` is used.
/// Ports `SOFTP2SV.containsKey(p5) && SOFTP2SV.get(p5).get(0).used`.
pub fn is_best_sv_used(var_data: &crate::variations::VariationData, softp: i32) -> bool {
    let mut best_cnt = -1;
    let mut best_used = false;

    let mut check_sv = |sv: &Sclip| {
        if sv.softp == softp && sv.variation.cnt > best_cnt {
            best_cnt = sv.variation.cnt;
            best_used = sv.used;
        }
    };

    for sv in &var_data.sv_structures.svfinv3 {
        check_sv(sv);
    }
    for sv in &var_data.sv_structures.svrinv3 {
        check_sv(sv);
    }
    for sv in &var_data.sv_structures.svfinv5 {
        check_sv(sv);
    }
    for sv in &var_data.sv_structures.svrinv5 {
        check_sv(sv);
    }
    for sv in &var_data.sv_structures.svfdel {
        check_sv(sv);
    }
    for sv in &var_data.sv_structures.svrdel {
        check_sv(sv);
    }
    for sv in &var_data.sv_structures.svfdup {
        check_sv(sv);
    }
    for sv in &var_data.sv_structures.svrdup {
        check_sv(sv);
    }

    // Deterministic forward-before-reverse traversal: each map's keys are
    // sorted, and svffus is fully scanned before svrfus, preserving Java's
    // precedence class (Sclip lists are per-key, so ties keep the first
    // highest-count candidate in forward-then-reverse order).
    let mut ffus_keys: Vec<&String> = var_data.sv_structures.svffus.keys().collect();
    ffus_keys.sort();
    for key in ffus_keys {
        if let Some(sv_list) = var_data.sv_structures.svffus.get(key) {
            for sv in sv_list {
                check_sv(sv);
            }
        }
    }
    let mut rfus_keys: Vec<&String> = var_data.sv_structures.svrfus.keys().collect();
    rfus_keys.sort();
    for key in rfus_keys {
        if let Some(sv_list) = var_data.sv_structures.svrfus.get(key) {
            for sv in sv_list {
                check_sv(sv);
            }
        }
    }

    best_cnt >= 0 && best_used
}

/// Represents a structural variant match returned by `find_match`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Match {
    pub base_position: i32,
    pub matched_sequence: Vec<u8>,
}

/// Helper function to check if the reference base at `pos` differs from `query_base`.
/// Ports `VariationUtils.isHasAndNotEquals`.
pub fn is_has_and_not_equals(reference: &Reference, chr: &str, pos: u32, query_base: u8) -> bool {
    if let Some(ref_base) = reference.get_base(chr, pos) {
        ref_base != query_base
    } else {
        false
    }
}

/// Find matching position on reference with given seed length and mismatch limit.
/// Ports `StructuralVariantsProcessor.findMatch` (both overloads, defaults MM=3).
pub fn find_match(
    seq: &[u8],
    reference: &Reference,
    position: i32,
    dir: i32,
    seed_len: usize,
    mm_limit: i32,
    chr: &str,
    verbose: bool,
) -> Match {
    let mut seq_vec = seq.to_vec();
    if dir == -1 {
        seq_vec.reverse();
    }

    if verbose {
        eprintln!(
            "    Working Match {} {} {} SEED: {}",
            position,
            String::from_utf8_lossy(&seq_vec),
            dir,
            seed_len
        );
    }

    let mut extra = Vec::new();
    if seq_vec.len() < seed_len {
        return Match {
            base_position: 0,
            matched_sequence: Vec::new(),
        };
    }

    for i in (0..=seq_vec.len() - seed_len).rev() {
        let Some(seed) = seq_vec.get(i..i + seed_len) else {
            continue;
        };
        if let Some([first_seed_u32]) = reference.seed.get(seed) {
            let first_seed = *first_seed_u32 as i32;
            let mut bp = if dir == 1 {
                first_seed - i as i32
            } else {
                first_seed + seq_vec.len() as i32 - i as i32 - 1
            };

            if crate::realign::is_match_ref(&seq_vec, reference, chr, bp, dir, mm_limit) {
                let mut mm = if dir == -1 { -1 } else { 0 };
                loop {
                    if let Some(ch) = crate::util::char_at(&seq_vec, mm) {
                        if is_has_and_not_equals(reference, chr, bp as u32, ch) {
                            extra.push(ch);
                            bp += dir;
                            mm += dir as isize;
                            continue;
                        }
                    }
                    break;
                }
                if !extra.is_empty() && dir == -1 {
                    extra.reverse();
                }
                if verbose {
                    eprintln!(
                        "      Found SV BP: {} BP: {} SEEDpos {} {} {} {} {} extra: {}",
                        dir,
                        bp,
                        first_seed,
                        position,
                        String::from_utf8_lossy(seed),
                        i,
                        String::from_utf8_lossy(&seq_vec),
                        String::from_utf8_lossy(&extra)
                    );
                }
                return Match {
                    base_position: bp,
                    matched_sequence: extra,
                };
            } else {
                // For complex indels, allowing some mismatches at the end up to 15bp or 20% length
                let mut sseq = seq_vec.clone();
                let mut eqcnt = 0;
                let mut bp_mut = bp;
                for ii in 1..=15 {
                    bp_mut += dir;
                    if dir == 1 {
                        if sseq.len() > 1 {
                            sseq = sseq[1..].to_vec();
                        } else {
                            sseq.clear();
                        }
                        if sseq.is_empty() {
                            break;
                        }
                        if is_has_and_not_equals(reference, chr, bp_mut as u32, sseq[0]) {
                            continue;
                        }
                        eqcnt += 1;
                        if sseq.len() > 1
                            && is_has_and_not_equals(reference, chr, (bp_mut + 1) as u32, sseq[1])
                        {
                            continue;
                        }
                        extra = seq_vec[0..ii].to_vec();
                    } else {
                        if sseq.len() > 1 {
                            sseq = sseq[0..sseq.len() - 1].to_vec();
                        } else {
                            sseq.clear();
                        }
                        if sseq.is_empty() {
                            break;
                        }
                        if is_has_and_not_equals(
                            reference,
                            chr,
                            bp_mut as u32,
                            sseq[sseq.len() - 1],
                        ) {
                            continue;
                        }
                        eqcnt += 1;
                        if sseq.len() > 1
                            && is_has_and_not_equals(
                                reference,
                                chr,
                                (bp_mut - 1) as u32,
                                sseq[sseq.len() - 2],
                            )
                        {
                            continue;
                        }
                        extra = seq_vec[seq_vec.len() - ii..].to_vec();
                    }

                    if eqcnt >= 3 && (eqcnt as f64) / (ii as f64) > 0.5 {
                        break;
                    }
                    if verbose {
                        eprintln!(
                            "      FoundSEED SV BP: {} BP: {} SEEDpos {} {} {} {} {} EXTRA: {}",
                            dir,
                            bp_mut,
                            first_seed,
                            position,
                            String::from_utf8_lossy(seed),
                            i,
                            String::from_utf8_lossy(&seq_vec),
                            String::from_utf8_lossy(&extra)
                        );
                    }
                    if crate::realign::is_match_ref(&sseq, reference, chr, bp_mut, dir, 1) {
                        return Match {
                            base_position: bp_mut,
                            matched_sequence: extra,
                        };
                    }
                }
            }
        }
    }

    Match {
        base_position: 0,
        matched_sequence: Vec::new(),
    }
}

/// Find matching position on reverse strand (for inversions).
pub fn find_match_rev(
    seq: &[u8],
    reference: &Reference,
    position: i32,
    dir: i32,
    seed_len: usize,
    mm_limit: i32,
    chr: &str,
    verbose: bool,
) -> Match {
    let mut seq_vec = seq.to_vec();
    if dir == 1 {
        seq_vec.reverse();
    }
    for base in seq_vec.iter_mut() {
        *base = crate::util::complement(*base);
    }

    if verbose {
        eprintln!(
            "    Working MatchRev {} {} {} SEED: {}",
            position,
            String::from_utf8_lossy(&seq_vec),
            dir,
            seed_len
        );
    }

    let mut extra = Vec::new();
    if seq_vec.len() < seed_len {
        return Match {
            base_position: 0,
            matched_sequence: Vec::new(),
        };
    }

    for i in (0..=seq_vec.len() - seed_len).rev() {
        let Some(seed) = seq_vec.get(i..i + seed_len) else {
            continue;
        };
        if let Some([first_seed_u32]) = reference.seed.get(seed) {
            let first_seed = *first_seed_u32 as i32;
            let mut bp = if dir == 1 {
                first_seed + seq_vec.len() as i32 - i as i32 - 1
            } else {
                first_seed - i as i32
            };

            if crate::realign::is_match_ref(&seq_vec, reference, chr, bp, -1 * dir, mm_limit) {
                if verbose {
                    eprintln!(
                        "      Found SV BP (reverse): {} BP: {} SEEDpos: {} {} {} {} {}",
                        dir,
                        bp,
                        first_seed,
                        position,
                        String::from_utf8_lossy(seed),
                        i,
                        String::from_utf8_lossy(&seq_vec)
                    );
                }
                return Match {
                    base_position: bp,
                    matched_sequence: Vec::new(),
                };
            } else {
                // For complex indels, allowing some mismatches at the end up to 15bp or 20% length
                let mut sseq = seq_vec.clone();
                let mut eqcnt = 0;
                let mut bp_mut = bp;
                for ii in 1..=15 {
                    bp_mut -= dir;
                    if dir == -1 {
                        if sseq.len() > 1 {
                            sseq = sseq[1..].to_vec();
                        } else {
                            sseq.clear();
                        }
                        if sseq.is_empty() {
                            break;
                        }
                        if is_has_and_not_equals(reference, chr, bp_mut as u32, sseq[0]) {
                            continue;
                        }
                        eqcnt += 1;
                        if sseq.len() > 1
                            && is_has_and_not_equals(reference, chr, (bp_mut + 1) as u32, sseq[1])
                        {
                            continue;
                        }
                        extra = seq_vec[0..ii].to_vec();
                    } else {
                        if sseq.len() > 1 {
                            sseq = sseq[0..sseq.len() - 1].to_vec();
                        } else {
                            sseq.clear();
                        }
                        if sseq.is_empty() {
                            break;
                        }
                        if is_has_and_not_equals(
                            reference,
                            chr,
                            bp_mut as u32,
                            sseq[sseq.len() - 1],
                        ) {
                            continue;
                        }
                        eqcnt += 1;
                        if sseq.len() > 1
                            && is_has_and_not_equals(
                                reference,
                                chr,
                                (bp_mut - 1) as u32,
                                sseq[sseq.len() - 2],
                            )
                        {
                            continue;
                        }
                        extra = seq_vec[seq_vec.len() - ii..].to_vec();
                    }

                    if eqcnt >= 3 && (eqcnt as f64) / (ii as f64) > 0.5 {
                        break;
                    }
                    if verbose {
                        eprintln!(
                            "      FoundSEED SV BP (reverse): {} BP: {} SEEDpos {} {} {} {} {} EXTRA: {}",
                            dir,
                            bp_mut,
                            first_seed,
                            position,
                            String::from_utf8_lossy(seed),
                            i,
                            String::from_utf8_lossy(&seq_vec),
                            String::from_utf8_lossy(&extra)
                        );
                    }
                    if crate::realign::is_match_ref(&sseq, reference, chr, bp_mut, -1 * dir, 1) {
                        return Match {
                            base_position: bp_mut,
                            matched_sequence: extra,
                        };
                    }
                }
            }
        }
    }

    Match {
        base_position: 0,
        matched_sequence: Vec::new(),
    }
}

/// Ports `StructuralVariantsProcessor.isOverlap`.
pub fn is_overlap(start1: i32, end1: i32, start2: i32, end2: i32, rlen: i32) -> bool {
    if start1 >= end2 || start2 >= end1 {
        return false;
    }
    let mut positions = [start1, end1, start2, end2];
    positions.sort_unstable();

    let ins = positions[2] - positions[1];
    if (end1 != start1) && (end2 != start2) {
        let ratio1 = (ins as f64) / ((end1 - start1) as f64);
        let ratio2 = (ins as f64) / ((end2 - start2) as f64);
        if ratio1 > 0.75 && ratio2 > 0.75 {
            return true;
        }
    }
    if (positions[1] - positions[0]) + (positions[3] - positions[2]) < 3 * rlen {
        return true;
    }
    false
}

/// Ports `StructuralVariantsProcessor.markSV`.
pub fn mark_sv(
    start: i32,
    end: i32,
    structural_variants_sv: &mut [&mut Vec<Sclip>],
    rlen: i32,
    verbose: bool,
) -> (i32, i32, i32) {
    let mut cov = 0;
    let mut pairs = 0;
    let mut cnt = 0;

    for current_sclips_sr in structural_variants_sv {
        for sv_r in current_sclips_sr.iter_mut() {
            let (start2, end2) = if sv_r.start < sv_r.mstart {
                (sv_r.end, sv_r.mstart)
            } else {
                (sv_r.mend, sv_r.start)
            };

            if verbose {
                eprintln!(
                    "   Marking SV {} {} {} {} cnt: {}",
                    start, end, start2, end2, sv_r.variation.cnt
                );
            }

            if is_overlap(start, end, start2, end2, rlen) {
                if verbose {
                    eprintln!(
                        "       SV {} {} {} {} cnt: {} marked",
                        start, end, start2, end2, sv_r.variation.cnt
                    );
                }
                sv_r.used = true;
                cnt += 1;
                pairs += sv_r.variation.cnt;
                if sv_r.end != sv_r.start {
                    cov += ((sv_r.variation.cnt * rlen) as f64 / (sv_r.end - sv_r.start) as f64)
                        as i32
                        + 1;
                }
            }
        }
    }
    (cov, cnt, pairs)
}

/// Ports `StructuralVariantsProcessor.markDUPSV`.
pub fn mark_dup_sv(
    start: i32,
    end: i32,
    structural_variants_sv: &mut [&mut Vec<Sclip>],
    rlen: i32,
    verbose: bool,
) -> (i32, i32) {
    let mut pairs = 0;
    let mut cnt = 0;

    for current_sclips_sr in structural_variants_sv {
        for sv_r in current_sclips_sr.iter_mut() {
            let (start2, end2) = if sv_r.start < sv_r.mstart {
                (sv_r.start, sv_r.mend)
            } else {
                (sv_r.mstart, sv_r.end)
            };

            if verbose {
                eprintln!(
                    "   Marking DUP SV {} {} {} {} cnt: {}",
                    start, end, start2, end2, sv_r.variation.cnt
                );
            }

            if is_overlap(start, end, start2, end2, rlen) {
                if verbose {
                    eprintln!(
                        "       DUP SV {} {} {} {} cnt: {} marked",
                        start, end, start2, end2, sv_r.variation.cnt
                    );
                }
                sv_r.used = true;
                cnt += 1;
                pairs += sv_r.variation.cnt;
            }
        }
    }
    (cnt, pairs)
}

#[derive(Debug, Clone, Copy)]
pub struct PairsData {
    pub pairs: i32,
    pub pmean: f64,
    pub qmean: f64,
    pub qmean_mapping: f64,
    pub nm: f64,
}

/// Check discordant pair mapping evidence for structural variants.
/// Ports `StructuralVariantsProcessor.checkPairs`.
pub fn check_pairs(
    start: i32,
    end: i32,
    sv: &mut [&mut Vec<Sclip>],
    max_read_length: i32,
    verbose: bool,
) -> PairsData {
    let mut pairs = 0;
    let mut pairs_data = PairsData {
        pairs: 0,
        pmean: 0.0,
        qmean: 0.0,
        qmean_mapping: 0.0,
        nm: 0.0,
    };

    for sv_cluster in sv {
        for svr in sv_cluster.iter_mut() {
            if svr.used {
                continue;
            }
            let mut s = (svr.start + svr.end) / 2;
            let mut e = (svr.mstart + svr.mend) / 2;
            if s > e {
                std::mem::swap(&mut s, &mut e);
            }
            if !is_overlap(start, end, s, e, max_read_length) {
                continue;
            }
            if svr.variation.cnt > pairs {
                pairs_data = PairsData {
                    pairs: svr.variation.cnt,
                    pmean: svr.variation.pmean,
                    qmean: svr.variation.qmean,
                    qmean_mapping: svr.variation.mqmean,
                    nm: svr.variation.nm,
                };
                pairs = svr.variation.cnt;
            }
            svr.used = true;
            if verbose {
                eprintln!("      Pair [{}:{}] found and marked.", s, e);
            }
        }
    }
    pairs_data
}

struct SortPositionSclip {
    position: i32,
    sclip: Sclip,
    count: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    Five,
    Three,
}

/// Helper function to perform structural variant calling for inversions.
/// Ports `StructuralVariantsProcessor.findINVsub`.
pub fn find_inv_sub(
    svref_idx: usize,
    dir: i32,
    side: Side,
    var_data: &mut crate::variations::VariationData,
    bam_path: &std::path::Path,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let minr = config.min_reads as i32;
    let verbose = config.verbose;
    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_maps(&config.adaptor);
    let max_read_length = var_data.max_read_length;

    let n = match svref_idx {
        0 => var_data.sv_structures.svfinv5.len(),
        1 => var_data.sv_structures.svrinv5.len(),
        2 => var_data.sv_structures.svfinv3.len(),
        3 => var_data.sv_structures.svrinv3.len(),
        _ => {
            return Err(VardictError::Other(format!(
                "internal error: invalid inversion candidate index {svref_idx}"
            )));
        }
    };
    for idx in 0..n {
        let (used, cnt, mstart, mend, start, end) = {
            let inv = match svref_idx {
                0 => &var_data.sv_structures.svfinv5[idx],
                1 => &var_data.sv_structures.svrinv5[idx],
                2 => &var_data.sv_structures.svfinv3[idx],
                3 => &var_data.sv_structures.svrinv3[idx],
                _ => {
                    return Err(VardictError::Other(format!(
                        "internal error: invalid inversion candidate index {svref_idx}"
                    )));
                }
            };
            (
                inv.used,
                inv.variation.cnt,
                inv.mstart,
                inv.mend,
                inv.start,
                inv.end,
            )
        };
        if used || cnt < minr {
            continue;
        }

        let mut soft_sorted: Vec<(i32, i32)> = {
            let inv = match svref_idx {
                0 => &var_data.sv_structures.svfinv5[idx],
                1 => &var_data.sv_structures.svrinv5[idx],
                2 => &var_data.sv_structures.svfinv3[idx],
                3 => &var_data.sv_structures.svrinv3[idx],
                _ => {
                    return Err(VardictError::Other(format!(
                        "internal error: invalid inversion candidate index {svref_idx}"
                    )));
                }
            };
            inv.soft.iter().map(|(&k, &v)| (k, v)).collect()
        };
        soft_sorted.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
        let mut softp = if soft_sorted.is_empty() {
            0
        } else {
            soft_sorted[0].0
        };

        // W3.3b: removed the `softp != 0 && is_best_sv_used(softp)` skip.
        // Java findINVsub has NO SOFTP2SV-style gate here — its only skip is
        // scv.used on the fetched clip entry (:534, or :561 in the scan
        // branch); the SOFTP2SV best-cluster check exists solely in findsv
        // (StructuralVariantsProcessor.java:705/:835), which is ported where
        // it belongs (find_sv). The extra gate skipped clusters whose
        // highest-cnt anchor cluster was used (e.g. marked by findDEL's
        // markSV or filterSV) while their clip entry was still unused —
        // records Java produces went missing.

        if !is_loaded(chr, mstart as u32, mend as u32, reference) {
            let mod_region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(mstart.max(1) as u32),
                end: Pos1(mend as u32),
                gene: String::new(),
                amplicon: None,
            };
            get_reference(&mod_region, 500, chromosomes_lengths, config, reference)?;
            partial_pipeline(
                bam_path,
                var_data,
                reference,
                chr,
                mstart - 200,
                mend + 200,
                config,
                chromosomes_lengths,
            )?;
        }

        let mut bp = 0;
        let mut scv_clone = None;
        let mut scv_used_pos: i32 = softp;
        let mut seq = Vec::new();
        let mut extra: Vec<u8> = Vec::new();

        if softp != 0 {
            if dir == 1 {
                if let Some(scv) = var_data.soft_clips_3_end.get(&softp) {
                    if !scv.used {
                        scv_clone = Some(scv.clone());
                    }
                }
            } else {
                if let Some(scv) = var_data.soft_clips_5_end.get(&softp) {
                    if !scv.used {
                        scv_clone = Some(scv.clone());
                    }
                }
            }

            let mut scv = match scv_clone.as_ref() {
                Some(s) => s.clone(),
                None => continue,
            };

            seq = crate::variations::util::find_conseq(
                &mut scv,
                0,
                &adaptor_forward,
                &adaptor_reverse,
            );
            if seq.is_empty() {
                continue;
            }

            let match_res = find_match_rev(&seq, reference, softp, dir, 17, 3, chr, verbose);
            bp = match_res.base_position;
            extra = match_res.matched_sequence;
            if bp == 0 {
                let match_res_lenient =
                    find_match_rev(&seq, reference, softp, dir, 12, 0, chr, verbose);
                bp = match_res_lenient.base_position;
                extra = match_res_lenient.matched_sequence;
            }
            if bp == 0 {
                continue;
            }
        } else {
            // Look within 2*maxReadLength to find soft clips not associated with discordant pairs
            let sp = if dir == 1 { end } else { start };
            for i in 1..=(2 * max_read_length) {
                let cp = sp + i * dir;
                let mut scv_check = None;
                if dir == 1 {
                    if let Some(scv) = var_data.soft_clips_3_end.get(&cp) {
                        if !scv.used {
                            scv_check = Some(scv.clone());
                        }
                    }
                } else {
                    if let Some(scv) = var_data.soft_clips_5_end.get(&cp) {
                        if !scv.used {
                            scv_check = Some(scv.clone());
                        }
                    }
                }

                let mut scv = match scv_check {
                    Some(s) => s,
                    None => continue,
                };

                let seq_check = crate::variations::util::find_conseq(
                    &mut scv,
                    0,
                    &adaptor_forward,
                    &adaptor_reverse,
                );
                if seq_check.is_empty() {
                    continue;
                }

                let match_res = find_match_rev(&seq_check, reference, cp, dir, 17, 3, chr, verbose);
                bp = match_res.base_position;
                extra = match_res.matched_sequence;
                if bp == 0 {
                    let match_res_lenient =
                        find_match_rev(&seq_check, reference, cp, dir, 12, 0, chr, verbose);
                    bp = match_res_lenient.base_position;
                    extra = match_res_lenient.matched_sequence;
                }
                if bp == 0 {
                    continue;
                }
                softp = cp;
                scv_clone = Some(scv);
                scv_used_pos = cp;
                seq = seq_check;
                let bp_diff = if dir == 1 {
                    (bp - mend).abs()
                } else {
                    (bp - mstart).abs()
                };
                if (bp_diff as f64) < 1.5 * (max_read_length as f64) {
                    // MINSVCDIST = 1.5
                    break;
                }
            }
            if bp == 0 {
                continue;
            }
        }

        let mut scv = match scv_clone {
            Some(s) => s,
            None => continue,
        };

        if side == Side::Five {
            if dir == -1 {
                bp -= 1;
            }
        } else {
            if dir == 1 {
                bp += 1;
                if bp != 0 {
                    softp -= 1;
                }
            } else {
                softp -= 1;
            }
        }

        if side == Side::Three {
            std::mem::swap(&mut bp, &mut softp);
        }

        // Refine bp and softp by shifting while complement matches
        if (dir == -1 && side == Side::Five) || (dir == 1 && side == Side::Three) {
            while let (Some(ch_s), Some(ch_b)) = (
                reference.get_base(chr, softp as u32),
                reference.get_base(chr, bp as u32),
            ) {
                if ch_s == crate::util::complement(ch_b) {
                    softp += 1;
                    if softp != 0 {
                        bp -= 1;
                    }
                } else {
                    break;
                }
            }
        }
        while let (Some(ch_s), Some(ch_b)) = (
            reference.get_base(chr, (softp - 1) as u32),
            reference.get_base(chr, (bp + 1) as u32),
        ) {
            if ch_s == crate::util::complement(ch_b) {
                softp -= 1;
                if softp != 0 {
                    bp += 1;
                }
            } else {
                break;
            }
        }

        let inv_mlen = {
            let inv = match svref_idx {
                0 => &var_data.sv_structures.svfinv5[idx],
                1 => &var_data.sv_structures.svrinv5[idx],
                2 => &var_data.sv_structures.svfinv3[idx],
                3 => &var_data.sv_structures.svrinv3[idx],
                _ => {
                    return Err(VardictError::Other(format!(
                        "internal error: invalid inversion candidate index {svref_idx}"
                    )));
                }
            };
            inv.mlen
        };

        // ports StructuralVariantsProcessor.java:627-628: the gate compares
        // bp - softp (> 150, and /|mlen| < 1.5) BEFORE len = bp - softp + 1
        // is computed inside the block. Gating on len instead accepted a gap
        // of exactly 150 and divided (gap+1) by |mlen| — e.g. gap=299,
        // mlen=200 is accepted by Java (1.495 < 1.5) but was rejected here.
        let gap = bp - softp;
        let abs_mlen = (inv_mlen as f64).abs();
        if bp > softp && gap > 150 && (gap as f64) / abs_mlen < 1.5 {
            let len = gap + 1;
            let sv_flank = 50; // Configuration.SVFLANK
            let ins5 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                bp - sv_flank + 1,
                bp,
            ));
            let ins3 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                softp,
                softp + sv_flank - 1,
            ));
            let ins = if len - 2 * sv_flank <= 0 {
                complement_seq(&crate::util::join_ref(reference, chr, softp, bp))
            } else {
                let mut tmp_ins = ins5;
                tmp_ins.push_str(&format!("<inv{}>", len - 2 * sv_flank));
                tmp_ins.push_str(&ins3);
                tmp_ins
            };
            let ins = apply_inv_extra(ins, dir, &extra);

            let gt_full = format!("-{}^{}", len, ins);

            {
                let variation = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    softp,
                    &gt_full,
                );
                variation.pstd = true;
                variation.qstd = true;
                crate::variations::util::adj_cnt(variation, &scv.variation, None);
            }

            {
                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    softp,
                );
                sv.r#type = "INV".to_string();
                sv.splits += scv.variation.cnt;
                sv.pairs += cnt;
                sv.clusters += 1;
            }

            if dir == -1 {
                if let Some(ref_base) = reference.get_base(chr, softp as u32) {
                    let ref_base_str = (ref_base as char).to_string();
                    if let Some(entry) = var_data.non_insertion_variants.get_mut(&softp) {
                        if let Some(ref_var) = entry.map.get_mut(&ref_base_str) {
                            ref_var.cnt -= scv.variation.cnt;
                            ref_var.hicnt -= scv.variation.hicnt;
                            ref_var.locnt -= scv.variation.locnt;
                            ref_var.pmean -= scv.variation.pmean;
                            ref_var.qmean -= scv.variation.qmean;
                            ref_var.mqmean -= scv.variation.mqmean;
                            ref_var.nm -= scv.variation.nm;
                            ref_var.sub_dir(true, scv.variation.get_dir(true));
                            ref_var.sub_dir(false, scv.variation.get_dir(false));
                            crate::variations::util::correct_cnt(ref_var);
                        }
                    }
                }
            }

            let mut dels5 = BTreeMap::new();
            let mut map = HashMap::new();
            map.insert(gt_full.clone(), cnt);
            dels5.insert(softp, map);

            let ref_cov_val = if let Some(&prev_cov) = var_data.ref_coverage.get(&(softp - 1)) {
                prev_cov
            } else {
                cnt
            };
            var_data.ref_coverage.insert(softp, ref_cov_val);
            // ports StructuralVariantsProcessor.java:663: mark the LIVE soft clip
            // map entry used (dir==1 -> 3' end map), under the position it was
            // fetched from — Java mutates the object itself, which survives the
            // subsequent softp/bp shifts.
            if dir == 1 {
                if let Some(scv_live) = var_data.soft_clips_3_end.get_mut(&scv_used_pos) {
                    scv_live.used = true;
                }
            } else if let Some(scv_live) = var_data.soft_clips_5_end.get_mut(&scv_used_pos) {
                scv_live.used = true;
            }

            let bams_vec: Vec<String> = config.bam.split('|').map(|s| s.to_string()).collect();
            let bams_opt = if bams_vec.is_empty() {
                None
            } else {
                Some(&bams_vec[..])
            };
            let chr_len = *chromosomes_lengths.get(chr).unwrap_or(&0) as u32;

            crate::realign::realigndel(
                bams_opt,
                &dels5,
                reference,
                chr,
                chr_len,
                config,
                &mut var_data.non_insertion_variants,
                &mut var_data.ref_coverage,
                &mut var_data.soft_clips_5_end,
                &mut var_data.soft_clips_3_end,
                max_read_length,
            );

            // Update the used status of the processed candidate inversion
            match svref_idx {
                0 => var_data.sv_structures.svfinv5[idx].used = true,
                1 => var_data.sv_structures.svrinv5[idx].used = true,
                2 => var_data.sv_structures.svfinv3[idx].used = true,
                3 => var_data.sv_structures.svrinv3[idx].used = true,
                _ => {
                    return Err(VardictError::Other(format!(
                        "internal error: invalid inversion candidate index {svref_idx}"
                    )));
                }
            };

            // Java findINVsub's success block (StructuralVariantsProcessor.java
            // 643-675) performs exactly three mutations before `return vref;`:
            // inv.used = true (:644), scv.used = true (:663, the LIVE map entry
            // fetched at the original position — handled above via scv_used_pos)
            // and the embedded realigndel (:665-667). It does NOT mark the soft
            // clip under the post-shift/swapped softp, and it does NOT run
            // markSV over any cluster lists: VariationRealigner.realigndel
            // contains no markSV call (the markSV sites live in realignlgdel
            // VR:1032/VR:1244 which ran earlier in realignIndels, and in
            // findDEL/findDELdisc/findINVdisc which run before/after findINV),
            // and the earlier W1.2 trace proved Java mutates the fetched entry
            // object rather than re-marking by position. Extra used-marking
            // here suppressed detections Java performs (missing-record class).
            // W3.3a: removed the mark_sv x4 sweep + shifted-softp marking.

            // Java findINVsub returns after the FIRST successful INV detection
            // (StructuralVariantsProcessor.java:675 `return vref;`) — later
            // clusters in this list are left for nothing else (the call ends).
            // Continuing the loop here made Rust build extra INVs from clusters
            // Java never reached (NA24385 chr20:20201754 duplicate INV source).
            return Ok(());
        }
    }

    Ok(())
}

/// Find inversion structural variants using softclip mappings.
/// Ports `StructuralVariantsProcessor.findINV`.
pub fn find_inv(
    var_data: &mut crate::variations::VariationData,
    bam_path: &std::path::Path,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    find_inv_sub(
        0,
        1,
        Side::Five,
        var_data,
        bam_path,
        reference,
        chr,
        config,
        chromosomes_lengths,
    )?;
    find_inv_sub(
        1,
        -1,
        Side::Five,
        var_data,
        bam_path,
        reference,
        chr,
        config,
        chromosomes_lengths,
    )?;
    find_inv_sub(
        2,
        1,
        Side::Three,
        var_data,
        bam_path,
        reference,
        chr,
        config,
        chromosomes_lengths,
    )?;
    find_inv_sub(
        3,
        -1,
        Side::Three,
        var_data,
        bam_path,
        reference,
        chr,
        config,
        chromosomes_lengths,
    )?;
    Ok(())
}

/// Find deletion structural variants using softclip mappings.
/// Ports `StructuralVariantsProcessor.findDEL`.
pub fn find_del(
    var_data: &mut crate::variations::VariationData,
    bam_path: &std::path::Path,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let minr = config.min_reads as i32;
    let verbose = config.verbose;
    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_maps(&config.adaptor);

    // 5' deletion softclips
    let n_svfdel = var_data.sv_structures.svfdel.len();
    for idx in 0..n_svfdel {
        let (used, vars_count, start, end, mstart, mend, softp_cand, del_variation) = {
            let del = &var_data.sv_structures.svfdel[idx];
            if del.used || del.variation.cnt < minr {
                continue;
            }
            let mut soft_sorted: Vec<(i32, i32)> = del.soft.iter().map(|(&k, &v)| (k, v)).collect();
            soft_sorted.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
            let softp = if soft_sorted.is_empty() {
                0
            } else {
                soft_sorted[0].0
            };
            (
                del.used,
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                softp,
                del.variation.clone(),
            )
        };

        if used || vars_count < minr {
            continue;
        }

        if softp_cand != 0 {
            let mut scv_clone = None;
            if let Some(scv) = var_data.soft_clips_3_end.get(&softp_cand) {
                if !scv.used {
                    scv_clone = Some(scv.clone());
                }
            }
            let mut scv = match scv_clone {
                Some(s) => s,
                None => {
                    continue;
                }
            };

            let seq = crate::variations::util::find_conseq(
                &mut scv,
                0,
                &adaptor_forward,
                &adaptor_reverse,
            );
            if seq.is_empty() || seq.len() < 12 {
                continue;
            }

            if !is_loaded(chr, mstart as u32, mend as u32, reference) {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(mstart.max(1) as u32),
                    end: Pos1(mend as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
                partial_pipeline(
                    bam_path,
                    var_data,
                    reference,
                    chr,
                    mstart - 200,
                    mend + 200,
                    config,
                    chromosomes_lengths,
                )?;
            }

            let match_res = find_match(&seq, reference, softp_cand, 1, 17, 3, chr, verbose);
            let mut bp = match_res.base_position;
            if bp == 0 {
                continue;
            }
            if !(bp - softp_cand > 30
                && is_overlap(softp_cand, bp, end, mstart, var_data.max_read_length))
            {
                continue;
            }
            bp -= 1;
            let mut softp = softp_cand;
            let mut dellen = bp - softp + 1;
            while reference.get_base(chr, bp as u32).is_some()
                && reference.get_base(chr, (softp - 1) as u32).is_some()
                && reference.get_base(chr, bp as u32) == reference.get_base(chr, (softp - 1) as u32)
            {
                bp -= 1;
                if bp != 0 {
                    softp -= 1;
                }
            }
            dellen = bp - softp + 1;

            {
                let variation = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    softp,
                    &format!("-{}", dellen),
                );
                variation.cnt = 0;
            }

            {
                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    softp,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs += vars_count;
                sv.splits += scv.variation.cnt;
                sv.clusters += 1;
            }

            let ref_cov_val = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
            if ref_cov_val < vars_count {
                var_data.ref_coverage.insert(softp, vars_count);
            }
            let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
            let softp_cov = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
            if softp_cov < bp_cov {
                var_data.ref_coverage.insert(softp, bp_cov);
            }

            {
                let variation_mut = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    softp,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(variation_mut, &scv.variation, None);

                if let Some(ref_base) = reference.get_base(chr, softp as u32) {
                    let ref_base_str = (ref_base as char).to_string();
                    if let Some(entry) = var_data.non_insertion_variants.get_mut(&softp) {
                        if let Some(ref_var) = entry.map.get_mut(&ref_base_str) {
                            ref_var.cnt -= scv.variation.cnt;
                            ref_var.hicnt -= scv.variation.hicnt;
                            ref_var.locnt -= scv.variation.locnt;
                            ref_var.pmean -= scv.variation.pmean;
                            ref_var.qmean -= scv.variation.qmean;
                            ref_var.mqmean -= scv.variation.mqmean;
                            ref_var.nm -= scv.variation.nm;
                            ref_var.sub_dir(true, scv.variation.get_dir(true));
                            ref_var.sub_dir(false, scv.variation.get_dir(false));
                            crate::variations::util::correct_cnt(ref_var);
                        }
                    }
                }
            }

            let mcnt = vars_count;
            let mut tv = Variation::new();
            tv.cnt = mcnt;
            tv.hicnt = mcnt;
            tv.fwd = mcnt / 2;
            tv.rev = mcnt - mcnt / 2;
            tv.qmean = del_variation.qmean * (mcnt as f64) / (vars_count as f64);
            tv.pmean = del_variation.pmean * (mcnt as f64) / (vars_count as f64);
            tv.mqmean = del_variation.mqmean * (mcnt as f64) / (vars_count as f64);
            tv.nm = del_variation.nm * (mcnt as f64) / (vars_count as f64);

            {
                let variation_mut2 = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    softp,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(variation_mut2, &tv, None);
            }

            var_data.sv_structures.svfdel[idx].used = true;
            if let Some(scv_mut) = var_data.soft_clips_3_end.get_mut(&softp_cand) {
                scv_mut.used = true;
            }

            let mut svrdel_mut_refs: Vec<&mut Vec<Sclip>> =
                vec![&mut var_data.sv_structures.svrdel];
            mark_sv(
                softp,
                bp,
                &mut svrdel_mut_refs,
                var_data.max_read_length,
                verbose,
            );
        } else {
            if !is_loaded(chr, mstart as u32, mend as u32, reference) {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(mstart.max(1) as u32),
                    end: Pos1(mend as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
            }
            let mut matched_j = None;
            let max_read_length = var_data.max_read_length;
            for (&i, scv) in &mut var_data.soft_clips_3_end {
                if scv.used {
                    continue;
                }
                if !(i >= end - 3 && i - end < 3 * max_read_length) {
                    continue;
                }
                let seq = crate::variations::util::find_conseq(
                    scv,
                    0,
                    &adaptor_forward,
                    &adaptor_reverse,
                );
                if seq.is_empty() || seq.len() < 12 {
                    continue;
                }
                let match_res = find_match(&seq, reference, i, 1, 17, 3, chr, verbose);
                let mut bp = match_res.base_position;
                if bp == 0 {
                    let match_res_lenient = find_match(&seq, reference, i, 1, 12, 0, chr, verbose);
                    bp = match_res_lenient.base_position;
                }
                if bp == 0 {
                    continue;
                }
                if !(bp - i > 30 && is_overlap(i, bp, end, mstart, var_data.max_read_length)) {
                    continue;
                }
                bp -= 1;
                matched_j = Some((i, bp, scv.clone()));
                break;
            }

            if let Some((softp, bp, scv)) = matched_j {
                let dellen = bp - softp + 1;
                {
                    let variation = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        softp,
                        &format!("-{}", dellen),
                    );
                    variation.cnt = 0;
                }

                {
                    let sv = crate::variations::VariationMap::get_sv(
                        &mut var_data.non_insertion_variants,
                        softp,
                    );
                    sv.r#type = "DEL".to_string();
                    sv.pairs += vars_count;
                    sv.splits += scv.variation.cnt;
                    sv.clusters += 1;
                }

                let ref_cov_val = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
                if ref_cov_val < vars_count {
                    var_data.ref_coverage.insert(softp, vars_count);
                }
                let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
                let softp_cov = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
                if softp_cov < bp_cov {
                    var_data.ref_coverage.insert(softp, bp_cov);
                }

                {
                    let variation_mut = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        softp,
                        &format!("-{}", dellen),
                    );
                    crate::variations::util::adj_cnt(variation_mut, &scv.variation, None);
                }

                let mcnt = vars_count;
                let mut tv = Variation::new();
                tv.cnt = mcnt;
                tv.hicnt = mcnt;
                tv.fwd = mcnt / 2;
                tv.rev = mcnt - mcnt / 2;
                tv.qmean = del_variation.qmean * (mcnt as f64) / (vars_count as f64);
                tv.pmean = del_variation.pmean * (mcnt as f64) / (vars_count as f64);
                tv.mqmean = del_variation.mqmean * (mcnt as f64) / (vars_count as f64);
                tv.nm = del_variation.nm * (mcnt as f64) / (vars_count as f64);

                {
                    let variation_mut2 = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        softp,
                        &format!("-{}", dellen),
                    );
                    crate::variations::util::adj_cnt(variation_mut2, &tv, None);
                }

                var_data.sv_structures.svfdel[idx].used = true;
                if let Some(scv_mut) = var_data.soft_clips_3_end.get_mut(&softp) {
                    scv_mut.used = true;
                }

                let mut svrdel_mut_refs: Vec<&mut Vec<Sclip>> =
                    vec![&mut var_data.sv_structures.svrdel];
                mark_sv(
                    softp,
                    bp,
                    &mut svrdel_mut_refs,
                    var_data.max_read_length,
                    verbose,
                );
            }
        }
    }

    // 3' deletion softclips
    let n_svrdel = var_data.sv_structures.svrdel.len();
    for idx in 0..n_svrdel {
        let (used, vars_count, start, end, mstart, mend, softp_cand, del_variation) = {
            let del = &var_data.sv_structures.svrdel[idx];
            if del.used || del.variation.cnt < minr {
                continue;
            }
            let mut soft_sorted: Vec<(i32, i32)> = del.soft.iter().map(|(&k, &v)| (k, v)).collect();
            soft_sorted.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
            let softp = if soft_sorted.is_empty() {
                0
            } else {
                soft_sorted[0].0
            };
            (
                del.used,
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                softp,
                del.variation.clone(),
            )
        };

        if used || vars_count < minr {
            continue;
        }

        if softp_cand != 0 {
            let mut scv_clone = None;
            if let Some(scv) = var_data.soft_clips_5_end.get(&softp_cand) {
                if !scv.used {
                    scv_clone = Some(scv.clone());
                }
            }
            let mut scv = match scv_clone {
                Some(s) => s,
                None => continue,
            };

            let seq = crate::variations::util::find_conseq(
                &mut scv,
                0,
                &adaptor_forward,
                &adaptor_reverse,
            );
            if seq.is_empty() || seq.len() < 12 {
                continue;
            }

            let is_loaded_before = is_loaded(chr, mstart as u32, mend as u32, reference);
            if !is_loaded_before {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(mstart.max(1) as u32),
                    end: Pos1(mend as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
                partial_pipeline(
                    bam_path,
                    var_data,
                    reference,
                    chr,
                    mstart - 200,
                    mend + 200,
                    config,
                    chromosomes_lengths,
                )?;
            }

            let match_res = find_match(&seq, reference, softp_cand, -1, 17, 3, chr, verbose);
            let mut bp = match_res.base_position;
            if bp == 0 {
                let match_res_lenient =
                    find_match(&seq, reference, softp_cand, -1, 12, 0, chr, verbose);
                bp = match_res_lenient.base_position;
            }
            if bp == 0 {
                continue;
            }
            let overlap_ok = softp_cand - bp > 30
                && is_overlap(bp, softp_cand, mend, start, var_data.max_read_length);
            if !overlap_ok {
                continue;
            }
            bp += 1;
            let mut softp = softp_cand - 1;
            let mut dellen = softp - bp + 1;

            {
                let variation = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    &format!("-{}", dellen),
                );
                variation.cnt = 0;
            }

            {
                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    bp,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs += vars_count;
                sv.splits += scv.variation.cnt;
                sv.clusters += 1;
            }

            let ref_cov_val = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
            if ref_cov_val < vars_count {
                var_data.ref_coverage.insert(bp, vars_count);
            }
            let softp_cov = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
            let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
            if bp_cov < softp_cov {
                var_data.ref_coverage.insert(bp, softp_cov);
            }

            {
                let variation_mut = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(variation_mut, &scv.variation, None);
            }

            let mcnt = vars_count;
            let mut tv = Variation::new();
            tv.cnt = mcnt;
            tv.hicnt = mcnt;
            tv.fwd = mcnt / 2;
            tv.rev = mcnt - mcnt / 2;
            tv.qmean = del_variation.qmean * (mcnt as f64) / (vars_count as f64);
            tv.pmean = del_variation.pmean * (mcnt as f64) / (vars_count as f64);
            tv.mqmean = del_variation.mqmean * (mcnt as f64) / (vars_count as f64);
            tv.nm = del_variation.nm * (mcnt as f64) / (vars_count as f64);

            {
                let variation_mut2 = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(variation_mut2, &tv, None);
            }

            var_data.sv_structures.svrdel[idx].used = true;
            if let Some(scv_mut) = var_data.soft_clips_5_end.get_mut(&softp_cand) {
                scv_mut.used = true;
            }

            let mut svfdel_mut_refs: Vec<&mut Vec<Sclip>> =
                vec![&mut var_data.sv_structures.svfdel];
            mark_sv(
                bp,
                softp_cand,
                &mut svfdel_mut_refs,
                var_data.max_read_length,
                verbose,
            );
        } else {
            if !is_loaded(chr, mstart as u32, mend as u32, reference) {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(mstart.max(1) as u32),
                    end: Pos1(mend as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
            }
            let mut matched_j = None;
            let max_read_length = var_data.max_read_length;
            for (&i, scv) in &mut var_data.soft_clips_5_end {
                if scv.used {
                    continue;
                }
                if !(i <= start + 3 && start - i < 3 * max_read_length) {
                    continue;
                }
                let seq = crate::variations::util::find_conseq(
                    scv,
                    0,
                    &adaptor_forward,
                    &adaptor_reverse,
                );
                if seq.is_empty() || seq.len() < 12 {
                    continue;
                }
                let match_res = find_match(&seq, reference, i, -1, 17, 3, chr, verbose);
                let mut bp = match_res.base_position;
                if bp == 0 {
                    let match_res_lenient = find_match(&seq, reference, i, -1, 12, 0, chr, verbose);
                    bp = match_res_lenient.base_position;
                }
                if bp == 0 {
                    continue;
                }
                if !(i - bp > 30 && is_overlap(bp, i, mend, start, var_data.max_read_length)) {
                    continue;
                }
                bp += 1;
                matched_j = Some((i, bp, scv.clone()));
                break;
            }

            if let Some((softp_cand, bp, scv)) = matched_j {
                let mut softp = softp_cand - 1;
                let dellen = softp - bp + 1;

                {
                    let variation = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        bp,
                        &format!("-{}", dellen),
                    );
                    variation.cnt = 0;
                }

                {
                    let sv = crate::variations::VariationMap::get_sv(
                        &mut var_data.non_insertion_variants,
                        bp,
                    );
                    sv.r#type = "DEL".to_string();
                    sv.pairs += vars_count;
                    sv.splits += scv.variation.cnt;
                    sv.clusters += 1;
                }

                let ref_cov_val = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
                if ref_cov_val < vars_count {
                    var_data.ref_coverage.insert(bp, vars_count);
                }
                let softp_cov = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
                let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
                if bp_cov < softp_cov {
                    var_data.ref_coverage.insert(bp, softp_cov);
                }
                *var_data.ref_coverage.entry(bp).or_insert(0) += scv.variation.cnt;

                {
                    let variation_mut = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        bp,
                        &format!("-{}", dellen),
                    );
                    crate::variations::util::adj_cnt(variation_mut, &scv.variation, None);
                }

                let mcnt = vars_count;
                let mut tv = Variation::new();
                tv.cnt = mcnt;
                tv.hicnt = mcnt;
                tv.fwd = mcnt / 2;
                tv.rev = mcnt - mcnt / 2;
                tv.qmean = del_variation.qmean * (mcnt as f64) / (vars_count as f64);
                tv.pmean = del_variation.pmean * (mcnt as f64) / (vars_count as f64);
                tv.mqmean = del_variation.mqmean * (mcnt as f64) / (vars_count as f64);
                tv.nm = del_variation.nm * (mcnt as f64) / (vars_count as f64);

                {
                    let variation_mut2 = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        bp,
                        &format!("-{}", dellen),
                    );
                    crate::variations::util::adj_cnt(variation_mut2, &tv, None);
                }

                var_data.sv_structures.svrdel[idx].used = true;
                if let Some(scv_mut) = var_data.soft_clips_5_end.get_mut(&softp_cand) {
                    scv_mut.used = true;
                }

                let mut svfdel_mut_refs: Vec<&mut Vec<Sclip>> =
                    vec![&mut var_data.sv_structures.svfdel];
                mark_sv(
                    bp,
                    softp_cand,
                    &mut svfdel_mut_refs,
                    var_data.max_read_length,
                    verbose,
                );
            }
        }
    }

    Ok(())
}

/// Find structural variants (DEL, INV) using remaining softclip evidence.
/// Ports `StructuralVariantsProcessor.findsv` and `fillAndSortTmpSV`.
pub fn find_sv(
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    chr: &str,
    curseg_start: i32,
    curseg_end: i32,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let minr = config.min_reads as i32;
    let verbose = config.verbose;
    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_maps(&config.adaptor);

    // Ports StructuralVariantsProcessor.java:973 fillAndSortTmpSV:
    // if (position < CURSEG.start || position > CURSEG.end) continue;
    let mut tmp5: Vec<SortPositionSclip> = var_data
        .soft_clips_5_end
        .iter()
        .filter_map(|(&pos, scv)| {
            if scv.used || scv.variation.cnt < minr || pos < curseg_start || pos > curseg_end {
                None
            } else {
                Some(SortPositionSclip {
                    position: pos,
                    sclip: scv.clone(),
                    count: scv.variation.cnt,
                })
            }
        })
        .collect();
    tmp5.sort_by_key(|item| std::cmp::Reverse(item.count));

    let mut tmp3: Vec<SortPositionSclip> = var_data
        .soft_clips_3_end
        .iter()
        .filter_map(|(&pos, scv)| {
            if scv.used || scv.variation.cnt < minr || pos < curseg_start || pos > curseg_end {
                None
            } else {
                Some(SortPositionSclip {
                    position: pos,
                    sclip: scv.clone(),
                    count: scv.variation.cnt,
                })
            }
        })
        .collect();
    tmp3.sort_by_key(|item| std::cmp::Reverse(item.count));

    for tuple5 in tmp5 {
        let p5 = tuple5.position;
        let mut sc5v = tuple5.sclip;
        let cnt5 = tuple5.count;
        if cnt5 < minr {
            break;
        }

        // Check if soft clip actually got used in the meantime
        if let Some(sc_mut) = var_data.soft_clips_5_end.get(&p5) {
            if sc_mut.used {
                continue;
            }
        }

        // SOFTP2SV check: skip if the SV at this position was already used
        if is_best_sv_used(var_data, p5) {
            continue;
        }

        let seq =
            crate::variations::util::find_conseq(&mut sc5v, 0, &adaptor_forward, &adaptor_reverse);
        if seq.is_empty() || seq.len() < 12 {
            continue;
        }

        let match_res = find_match(&seq, reference, p5, -1, 17, 3, chr, verbose);
        let mut bp = match_res.base_position;
        if bp != 0 {
            if bp < p5 {
                let mut svfdel = std::mem::take(&mut var_data.sv_structures.svfdel);
                let mut svrdel = std::mem::take(&mut var_data.sv_structures.svrdel);
                let pairs_data = {
                    let mut refs = vec![&mut svfdel, &mut svrdel];
                    check_pairs(bp, p5, &mut refs, var_data.max_read_length, verbose)
                };
                var_data.sv_structures.svfdel = svfdel;
                var_data.sv_structures.svrdel = svrdel;

                let pairs = pairs_data.pairs;
                if pairs == 0 {
                    continue;
                }

                let p5_adj = p5 - 1;
                let bp_adj = bp + 1;
                let dellen = p5_adj - bp_adj + 1;

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp_adj,
                    &format!("-{}", dellen),
                );
                vref.cnt = 0;

                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    bp_adj,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs += pairs;
                sv.splits += cnt5;
                sv.clusters += 1;

                let bp_cov = var_data.ref_coverage.get(&bp_adj).copied().unwrap_or(0);
                if bp_cov < pairs + sc5v.variation.cnt {
                    var_data
                        .ref_coverage
                        .insert(bp_adj, pairs + sc5v.variation.cnt);
                }
                // Ports StructuralVariantsProcessor.java:732,748: Java mutates `p5--`
                // at line 732 before `refCoverage.get(p5 + 1)` at line 748, so `p5 + 1`
                // in Java is `p5_adj + 1` (the original `p5` before decrement).
                let p5_cov = var_data
                    .ref_coverage
                    .get(&(p5_adj + 1))
                    .copied()
                    .unwrap_or(0);
                let bp_cov_updated = var_data.ref_coverage.get(&bp_adj).copied().unwrap_or(0);
                if bp_cov_updated < p5_cov {
                    var_data.ref_coverage.insert(bp_adj, p5_cov);
                }

                let vref_mut = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp_adj,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(vref_mut, &sc5v.variation, None);

                let mut tmp_var = Variation::new();
                tmp_var.cnt = pairs;
                tmp_var.hicnt = pairs;
                tmp_var.fwd = pairs / 2;
                tv_adj_vals(&mut tmp_var, pairs_data);

                let vref_mut2 = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp_adj,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(vref_mut2, &tmp_var, None);

                if let Some(sc_mut) = var_data.soft_clips_5_end.get_mut(&p5) {
                    sc_mut.used = true;
                }
            }
        } else {
            let match_res_rev = find_match_rev(&seq, reference, p5, -1, 17, 3, chr, verbose);
            let mut bp = match_res_rev.base_position;
            if bp == 0 {
                continue;
            }
            if (bp - p5).abs() <= 50 {
                continue;
            }
            let mut p5_adj = p5;
            if bp < p5_adj {
                // Keep bp at 3' side and p5_adj at 5' side
                std::mem::swap(&mut bp, &mut p5_adj);
            }
            bp -= 1;
            while let (Some(b1), Some(b2)) = (
                reference.get_base(chr, (bp + 1) as u32),
                reference.get_base(chr, (p5_adj - 1) as u32),
            ) {
                if crate::util::complement(b1) != b2 {
                    break;
                }
                p5_adj -= 1;
                if p5_adj != 0 {
                    bp += 1;
                }
            }

            let ins5 = complement_seq(&crate::util::join_ref(reference, chr, bp - 50 + 1, bp));
            let ins3 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                p5_adj,
                p5_adj + 50 - 1,
            ));
            let mid = bp - p5_adj - (ins5.len() as i32) - (ins3.len() as i32) + 1;

            let extra = String::from_utf8_lossy(&match_res_rev.matched_sequence);
            let mut vn = format!("-{}^{}<inv{}>{}{}", bp - p5_adj + 1, ins5, mid, ins3, extra);
            if mid <= 0 {
                let tins = complement_seq(&crate::util::join_ref(reference, chr, p5_adj, bp));
                vn = format!("-{}^{}{}", bp - p5_adj + 1, tins, extra);
            }

            let vref =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p5_adj, &vn);
            vref.cnt = 0;

            let sv = crate::variations::VariationMap::get_sv(
                &mut var_data.non_insertion_variants,
                p5_adj,
            );
            sv.r#type = "INV".to_string();
            sv.pairs += 0;
            sv.splits += cnt5;
            sv.clusters += 0;

            let vref_mut =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p5_adj, &vn);
            crate::variations::util::adj_cnt(vref_mut, &sc5v.variation, None);

            let p5_cov = var_data.ref_coverage.get(&p5_adj).copied().unwrap_or(0);
            var_data.ref_coverage.insert(p5_adj, p5_cov + cnt5);
            let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
            let p5_cov_updated = var_data.ref_coverage.get(&p5_adj).copied().unwrap_or(0);
            if p5_cov_updated < bp_cov {
                var_data.ref_coverage.insert(p5_adj, bp_cov);
            }

            if let Some(sc_mut) = var_data.soft_clips_5_end.get_mut(&p5) {
                sc_mut.used = true;
            }
        }
    }

    for tuple3 in tmp3 {
        let p3 = tuple3.position;
        let mut sc3v = tuple3.sclip;
        let cnt3 = tuple3.count;
        if cnt3 < minr {
            break;
        }

        if let Some(sc_mut) = var_data.soft_clips_3_end.get(&p3) {
            if sc_mut.used {
                continue;
            }
        }

        // SOFTP2SV check: skip if the SV at this position was already used
        if is_best_sv_used(var_data, p3) {
            continue;
        }

        let seq =
            crate::variations::util::find_conseq(&mut sc3v, 0, &adaptor_forward, &adaptor_reverse);

        if seq.is_empty() || seq.len() < 12 {
            continue;
        }

        let match_res = find_match(&seq, reference, p3, 1, 17, 3, chr, verbose);
        let mut bp = match_res.base_position;

        if bp != 0 {
            if bp > p3 {
                let mut svfdel = std::mem::take(&mut var_data.sv_structures.svfdel);
                let mut svrdel = std::mem::take(&mut var_data.sv_structures.svrdel);
                let pairs_data = {
                    let mut refs = vec![&mut svfdel, &mut svrdel];
                    check_pairs(p3, bp, &mut refs, var_data.max_read_length, verbose)
                };
                var_data.sv_structures.svfdel = svfdel;
                var_data.sv_structures.svrdel = svrdel;

                let pairs = pairs_data.pairs;
                if pairs == 0 {
                    continue;
                }

                let dellen = bp - p3;
                let mut bp_adj = bp - 1;
                let mut p3_adj = p3;
                while reference.get_base(chr, bp_adj as u32).is_some()
                    && reference.get_base(chr, (p3_adj - 1) as u32).is_some()
                    && reference.get_base(chr, bp_adj as u32)
                        == reference.get_base(chr, (p3_adj - 1) as u32)
                {
                    bp_adj -= 1;
                    if bp_adj != 0 {
                        p3_adj -= 1;
                    }
                }

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    p3_adj,
                    &format!("-{}", dellen),
                );
                vref.cnt = 0;

                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    p3_adj,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs += pairs;
                sv.splits += cnt3;
                sv.clusters += 1;

                if !var_data.ref_coverage.contains_key(&p3_adj) {
                    var_data
                        .ref_coverage
                        .insert(p3_adj, pairs + sc3v.variation.cnt);
                }
                if var_data.ref_coverage.contains_key(&bp_adj)
                    && var_data.ref_coverage.get(&bp_adj).copied().unwrap_or(0)
                        < var_data.ref_coverage.get(&p3_adj).copied().unwrap_or(0)
                {
                    var_data.ref_coverage.insert(
                        bp_adj,
                        var_data.ref_coverage.get(&p3_adj).copied().unwrap_or(0),
                    );
                }

                let vref_mut = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    p3_adj,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(vref_mut, &sc3v.variation, None);

                let mut tmp_var = Variation::new();
                tmp_var.cnt = pairs;
                tmp_var.hicnt = pairs;
                tmp_var.fwd = pairs / 2;
                tv_adj_vals(&mut tmp_var, pairs_data);

                let vref_mut2 = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    p3_adj,
                    &format!("-{}", dellen),
                );
                crate::variations::util::adj_cnt(vref_mut2, &tmp_var, None);

                if let Some(sc_mut) = var_data.soft_clips_3_end.get_mut(&p3) {
                    sc_mut.used = true;
                }
            }
        } else {
            let match_res_rev = find_match_rev(&seq, reference, p3, 1, 17, 3, chr, verbose);
            let mut bp = match_res_rev.base_position;
            if bp == 0 {
                continue;
            }
            if (bp - p3).abs() <= 50 {
                continue;
            }
            let mut p3_adj = p3;
            if bp < p3_adj {
                std::mem::swap(&mut bp, &mut p3_adj);
                p3_adj += 1;
                bp -= 1;
            }
            while let (Some(b1), Some(b2)) = (
                reference.get_base(chr, (bp + 1) as u32),
                reference.get_base(chr, (p3_adj - 1) as u32),
            ) {
                if crate::util::complement(b1) != b2 {
                    break;
                }
                p3_adj -= 1;
                if p3_adj != 0 {
                    bp += 1;
                }
            }

            let ins5 = complement_seq(&crate::util::join_ref(reference, chr, bp - 50 + 1, bp));
            let ins3 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                p3_adj,
                p3_adj + 50 - 1,
            ));
            let mid = bp - p3_adj - 100 + 1;

            let extra = String::from_utf8_lossy(&match_res_rev.matched_sequence);
            let mut vn = format!("-{}^{}{}<inv{}>{}", bp - p3_adj + 1, extra, ins5, mid, ins3);
            if mid <= 0 {
                let tins = complement_seq(&crate::util::join_ref(reference, chr, p3_adj, bp));
                vn = format!("-{}^{}{}", bp - p3_adj + 1, extra, tins);
            }

            let vref =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p3_adj, &vn);
            vref.cnt = 0;

            let sv = crate::variations::VariationMap::get_sv(
                &mut var_data.non_insertion_variants,
                p3_adj,
            );
            sv.r#type = "INV".to_string();
            sv.pairs += 0;
            sv.splits += cnt3;
            sv.clusters += 0;

            let vref_mut =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p3_adj, &vn);
            crate::variations::util::adj_cnt(vref_mut, &sc3v.variation, None);

            let p3_cov = var_data.ref_coverage.get(&p3_adj).copied().unwrap_or(0);
            var_data.ref_coverage.insert(p3_adj, p3_cov + cnt3);
            let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
            let p3_cov_updated = var_data.ref_coverage.get(&p3_adj).copied().unwrap_or(0);
            if p3_cov_updated < bp_cov {
                var_data.ref_coverage.insert(p3_adj, bp_cov);
            }

            if let Some(sc_mut) = var_data.soft_clips_3_end.get_mut(&p3) {
                sc_mut.used = true;
            }
        }
    }

    Ok(())
}

/// ports StructuralVariantsProcessor.java:625-632: dir==1 prepends
/// revcomp(extra); dir==-1 appends extra as-is. `extra` is findMatchRev's
/// matched suffix, not the consensus sequence.
fn apply_inv_extra(ins: String, dir: i32, extra: &[u8]) -> String {
    if dir == 1 && !extra.is_empty() {
        let mut combined = complement_seq(&String::from_utf8_lossy(extra));
        combined.push_str(&ins);
        combined
    } else if dir == -1 && !extra.is_empty() {
        let mut out = ins;
        out.push_str(&String::from_utf8_lossy(extra));
        out
    } else {
        ins
    }
}

fn complement_seq(s: &str) -> String {
    let mut out = String::new();
    for b in s.as_bytes() {
        out.push(crate::util::complement(*b) as char);
    }
    // Reverse complement
    let mut out_bytes = out.into_bytes();
    out_bytes.reverse();
    String::from_utf8_lossy(&out_bytes).into_owned()
}

fn tv_adj_vals(tv: &mut Variation, pairs_data: PairsData) {
    tv.hicnt = pairs_data.pairs;
    tv.rev = pairs_data.pairs - tv.fwd;
    tv.pmean = pairs_data.pmean;
    tv.qmean = pairs_data.qmean;
    tv.mqmean = pairs_data.qmean_mapping;
    tv.nm = pairs_data.nm;
}

/// Find deletion structural variants from discordant pairs only.
/// Ports `StructuralVariantsProcessor.findDELdisc`.
pub fn find_del_disc(
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let minr = config.min_reads as i32;
    let verbose = config.verbose;
    let mindist = 8 * var_data.max_read_length;

    // svfdel (5' candidate deletions)
    let n_svfdel = var_data.sv_structures.svfdel.len();
    for idx in 0..n_svfdel {
        let (used, vars_count, start, end, mstart, mend, softp) = {
            let del = &var_data.sv_structures.svfdel[idx];
            if del.used || del.variation.cnt < minr + 5 {
                continue;
            }
            if !var_data.splice.is_empty() && del.mlen.abs() < 250000 {
                continue;
            }
            if del.mstart <= del.end + mindist {
                continue;
            }
            if del.variation.mqmean / (del.variation.cnt as f64) <= 35.0 {
                // DISCPAIRQUAL = 35
                continue;
            }
            let mlen = del.mstart - del.end - var_data.max_read_length / (del.variation.cnt + 1);
            if !(mlen > 0 && mlen > mindist) {
                continue;
            }
            let bp = if del.softp != 0 {
                del.softp
            } else {
                del.end + (var_data.max_read_length / (del.variation.cnt + 1)) / 2
            };
            (
                del.used,
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                bp,
            )
        };

        if used {
            continue;
        }

        let mlen = mstart - end - var_data.max_read_length / (vars_count + 1);
        if reference.get_base(chr, softp as u32).is_none() {
            let mod_region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1((softp - 150).max(1) as u32),
                end: Pos1((softp + 150) as u32),
                gene: String::new(),
                amplicon: None,
            };
            get_reference(
                &mod_region,
                if mlen < 1000 { mlen as u32 } else { 1000 },
                chromosomes_lengths,
                config,
                reference,
            )?;
        }

        let vref = crate::variations::get_variation(
            &mut var_data.non_insertion_variants,
            softp,
            &format!("-{}", mlen),
        );
        vref.cnt = 0;

        let sv =
            crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, softp);
        sv.r#type = "DEL".to_string();
        sv.splits += var_data
            .soft_clips_3_end
            .get(&(end + 1))
            .map(|s| s.variation.cnt)
            .unwrap_or(0);
        sv.splits += var_data
            .soft_clips_5_end
            .get(&mstart)
            .map(|s| s.variation.cnt)
            .unwrap_or(0);
        sv.pairs += vars_count;
        sv.clusters += 1;

        let mut tv = Variation::new();
        tv.cnt = 2 * vars_count;
        tv.hicnt = 2 * vars_count;
        tv.fwd = vars_count;
        tv.rev = vars_count;
        {
            let del = &var_data.sv_structures.svfdel[idx];
            tv.qmean = 2.0 * del.variation.qmean;
            tv.pmean = 2.0 * del.variation.pmean;
            tv.mqmean = 2.0 * del.variation.mqmean;
            tv.nm = 2.0 * del.variation.nm;
        }
        let vref_mut = crate::variations::get_variation(
            &mut var_data.non_insertion_variants,
            softp,
            &format!("-{}", mlen),
        );
        crate::variations::util::adj_cnt(vref_mut, &tv, None);

        if !var_data.ref_coverage.contains_key(&softp) {
            var_data.ref_coverage.insert(softp, 2 * vars_count);
        }
        var_data.sv_structures.svfdel[idx].used = true;

        let mut svrdel_mut_refs: Vec<&mut Vec<Sclip>> = vec![&mut var_data.sv_structures.svrdel];
        mark_sv(
            end,
            mstart,
            &mut svrdel_mut_refs,
            var_data.max_read_length,
            verbose,
        );
    }

    // svrdel (3' candidate deletions)
    let n_svrdel = var_data.sv_structures.svrdel.len();
    for idx in 0..n_svrdel {
        let (used, vars_count, start, end, mstart, mend, softp) = {
            let del = &var_data.sv_structures.svrdel[idx];
            if del.used || del.variation.cnt < minr + 5 {
                continue;
            }
            if !var_data.splice.is_empty() && del.mlen.abs() < 250000 {
                continue;
            }
            if del.start <= del.mend + mindist {
                continue;
            }
            if del.variation.mqmean / (del.variation.cnt as f64) <= 35.0 {
                continue;
            }
            let mlen = del.start - del.mend - var_data.max_read_length / (del.variation.cnt + 1);
            if !(mlen > 0 && mlen > mindist) {
                continue;
            }
            let bp = del.mend + ((var_data.max_read_length / (del.variation.cnt + 1)) / 2);
            (
                del.used,
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                bp,
            )
        };

        if used {
            continue;
        }

        let mlen = start - mend - var_data.max_read_length / (vars_count + 1);
        if reference.get_base(chr, softp as u32).is_none() {
            let mod_region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1((softp - 150).max(1) as u32),
                end: Pos1((softp + 150) as u32),
                gene: String::new(),
                amplicon: None,
            };
            get_reference(
                &mod_region,
                if mlen < 1000 { mlen as u32 } else { 1000 },
                chromosomes_lengths,
                config,
                reference,
            )?;
        }

        let vref = crate::variations::get_variation(
            &mut var_data.non_insertion_variants,
            softp,
            &format!("-{}", mlen),
        );
        vref.cnt = 0;

        let sv =
            crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, softp);
        sv.r#type = "DEL".to_string();
        sv.splits += var_data
            .soft_clips_3_end
            .get(&(mend + 1))
            .map(|s| s.variation.cnt)
            .unwrap_or(0);
        sv.splits += var_data
            .soft_clips_5_end
            .get(&start)
            .map(|s| s.variation.cnt)
            .unwrap_or(0);
        sv.pairs += vars_count;
        sv.clusters += 1;

        let mut tv = Variation::new();
        tv.cnt = 2 * vars_count;
        tv.hicnt = 2 * vars_count;
        tv.fwd = vars_count;
        tv.rev = vars_count;
        {
            let del = &var_data.sv_structures.svrdel[idx];
            tv.qmean = 2.0 * del.variation.qmean;
            tv.pmean = 2.0 * del.variation.pmean;
            tv.mqmean = 2.0 * del.variation.mqmean;
            tv.nm = 2.0 * del.variation.nm;
        }
        let vref_mut = crate::variations::get_variation(
            &mut var_data.non_insertion_variants,
            softp,
            &format!("-{}", mlen),
        );
        crate::variations::util::adj_cnt(vref_mut, &tv, None);

        if !var_data.ref_coverage.contains_key(&softp) {
            var_data.ref_coverage.insert(softp, 2 * vars_count);
        }
        let start_cov = var_data.ref_coverage.get(&start).copied().unwrap_or(0);
        let bp_cov = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
        if start_cov > 0 && bp_cov < start_cov {
            var_data.ref_coverage.insert(softp, start_cov);
        }
        var_data.sv_structures.svrdel[idx].used = true;

        let mod_region_mate = crate::regions::Region {
            chr: chr.to_string(),
            start: Pos1((mstart - 100).max(1) as u32),
            end: Pos1((mend + 100) as u32),
            gene: String::new(),
            amplicon: None,
        };
        get_reference(
            &mod_region_mate,
            200,
            chromosomes_lengths,
            config,
            reference,
        )?;

        let mut svfdel_mut_refs: Vec<&mut Vec<Sclip>> = vec![&mut var_data.sv_structures.svfdel];
        mark_sv(
            mend,
            start,
            &mut svfdel_mut_refs,
            var_data.max_read_length,
            verbose,
        );
    }

    Ok(())
}

/// Find inversion structural variants from discordant pairs only.
/// Ports `StructuralVariantsProcessor.findINVdisc`.
// Java-faithful `!(a < b)` / negated partial-ord comparisons — keep as-is.
#[allow(clippy::neg_cmp_op_on_partial_ord, clippy::nonminimal_bool)]
pub fn find_inv_disc(
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let minr = config.min_reads as i32;
    let verbose = config.verbose;
    let max_read_length = var_data.max_read_length;
    let sv_flank = 50; // Configuration.SVFLANK

    // Part 1: svfinv5 and svrinv5
    let n_svfinv5 = var_data.sv_structures.svfinv5.len();
    let n_svrinv5 = var_data.sv_structures.svrinv5.len();
    for i in 0..n_svfinv5 {
        if var_data.sv_structures.svfinv5[i].used {
            continue;
        }
        let cnt = var_data.sv_structures.svfinv5[i].variation.cnt;
        let me = var_data.sv_structures.svfinv5[i].mend;
        let ms = var_data.sv_structures.svfinv5[i].mstart;
        let end = var_data.sv_structures.svfinv5[i].end;
        let start = var_data.sv_structures.svfinv5[i].start;
        let nm = var_data.sv_structures.svfinv5[i].variation.nm;
        let pmean = var_data.sv_structures.svfinv5[i].variation.pmean;
        let qmean = var_data.sv_structures.svfinv5[i].variation.qmean;
        let Qmean = var_data.sv_structures.svfinv5[i].variation.mqmean;

        if !(Qmean / (cnt as f64) > 35.0) {
            continue;
        }

        for j in 0..n_svrinv5 {
            if var_data.sv_structures.svrinv5[j].used {
                continue;
            }
            let rcnt = var_data.sv_structures.svrinv5[j].variation.cnt;
            let rstart = var_data.sv_structures.svrinv5[j].start;
            let rms = var_data.sv_structures.svrinv5[j].mstart;
            let rnm = var_data.sv_structures.svrinv5[j].variation.nm;
            let rpmean = var_data.sv_structures.svrinv5[j].variation.pmean;
            let rqmean = var_data.sv_structures.svrinv5[j].variation.qmean;
            let rQmean = var_data.sv_structures.svrinv5[j].variation.mqmean;

            if !(rQmean / (rcnt as f64) > 35.0) {
                continue;
            }
            if !(cnt + rcnt > minr + 5) {
                continue;
            }

            if is_overlap(end, me, rstart, rms, max_read_length) {
                let bp = ((end + rstart) / 2).abs();
                let pe = ((me + rms) / 2).abs();

                if reference.get_base(chr, pe as u32).is_none() {
                    let modified_region = crate::regions::Region {
                        chr: chr.to_string(),
                        start: Pos1((pe - 150).max(1) as u32),
                        end: Pos1((pe + 150) as u32),
                        gene: String::new(),
                        amplicon: None,
                    };
                    get_reference(
                        &modified_region,
                        300,
                        chromosomes_lengths,
                        config,
                        reference,
                    )?;
                }

                let len = pe - bp + 1;
                let ins5 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    bp,
                    bp + sv_flank - 1,
                ));
                let ins3 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    pe - sv_flank + 1,
                    pe,
                ));

                let mut ins = format!("{}<inv{}>{}", ins3, len - 2 * sv_flank, ins5);
                if len - 2 * sv_flank <= 0 {
                    ins = complement_seq(&crate::util::join_ref(reference, chr, bp, pe));
                }

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    &format!("-{}^{}", len, ins),
                );

                var_data.sv_structures.svfinv5[i].used = true;
                var_data.sv_structures.svrinv5[j].used = true;
                vref.pstd = true;
                vref.qstd = true;

                let mut tmp = Variation::new();
                tmp.cnt = cnt + rcnt;
                tmp.hicnt = cnt + rcnt;
                tmp.fwd = cnt;
                tmp.rev = rcnt;
                tmp.qmean = qmean + rqmean;
                tmp.pmean = pmean + rpmean;
                tmp.mqmean = Qmean + rQmean;
                tmp.nm = nm + rnm;
                crate::variations::util::adj_cnt(vref, &tmp, None);

                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    bp,
                );
                sv.r#type = "INV".to_string();
                sv.pairs += cnt;
                sv.splits += var_data
                    .soft_clips_5_end
                    .get(&start)
                    .map_or(0, |s| s.variation.cnt);
                sv.splits += var_data
                    .soft_clips_5_end
                    .get(&ms)
                    .map_or(0, |s| s.variation.cnt);
                sv.clusters += 1;

                if !var_data.ref_coverage.contains_key(&bp) {
                    var_data.ref_coverage.insert(bp, 2 * cnt);
                }

                let mut svfinv3_mut_refs: Vec<&mut Vec<Sclip>> = vec![
                    &mut var_data.sv_structures.svfinv3,
                    &mut var_data.sv_structures.svrinv3,
                ];
                mark_sv(bp, pe, &mut svfinv3_mut_refs, max_read_length, verbose);
            }
        }
    }

    // Part 2: svfinv3 and svrinv3
    let n_svfinv3 = var_data.sv_structures.svfinv3.len();
    let n_svrinv3 = var_data.sv_structures.svrinv3.len();
    for i in 0..n_svfinv3 {
        if var_data.sv_structures.svfinv3[i].used {
            continue;
        }
        let cnt = var_data.sv_structures.svfinv3[i].variation.cnt;
        let me = var_data.sv_structures.svfinv3[i].mend;
        let end = var_data.sv_structures.svfinv3[i].end;
        let nm = var_data.sv_structures.svfinv3[i].variation.nm;
        let pmean = var_data.sv_structures.svfinv3[i].variation.pmean;
        let qmean = var_data.sv_structures.svfinv3[i].variation.qmean;
        let Qmean = var_data.sv_structures.svfinv3[i].variation.mqmean;

        for j in 0..n_svrinv3 {
            if var_data.sv_structures.svrinv3[j].used {
                continue;
            }
            let rcnt = var_data.sv_structures.svrinv3[j].variation.cnt;
            let rstart = var_data.sv_structures.svrinv3[j].start;
            let rms = var_data.sv_structures.svrinv3[j].mstart;
            let rnm = var_data.sv_structures.svrinv3[j].variation.nm;
            let rpmean = var_data.sv_structures.svrinv3[j].variation.pmean;
            let rqmean = var_data.sv_structures.svrinv3[j].variation.qmean;
            let rQmean = var_data.sv_structures.svrinv3[j].variation.mqmean;

            if !(rQmean / (rcnt as f64) > 35.0) {
                continue;
            }
            if !(cnt + rcnt > minr + 5) {
                continue;
            }

            if is_overlap(me, end, rms, rstart, max_read_length) {
                let pe = ((end + rstart) / 2).abs();
                let bp = ((me + rms) / 2).abs();

                if reference.get_base(chr, bp as u32).is_none() {
                    let modified_region = crate::regions::Region {
                        chr: chr.to_string(),
                        start: Pos1((bp - 150).max(1) as u32),
                        end: Pos1((bp + 150) as u32),
                        gene: String::new(),
                        amplicon: None,
                    };
                    get_reference(
                        &modified_region,
                        300,
                        chromosomes_lengths,
                        config,
                        reference,
                    )?;
                }

                let len = pe - bp + 1;
                let ins5 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    bp,
                    bp + sv_flank - 1,
                ));
                let ins3 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    pe - sv_flank + 1,
                    pe,
                ));

                let mut ins = format!("{}<inv{}>{}", ins3, len - 2 * sv_flank, ins5);
                if len - 2 * sv_flank <= 0 {
                    ins = complement_seq(&crate::util::join_ref(reference, chr, bp, pe));
                }

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    &format!("-{}^{}", len, ins),
                );

                var_data.sv_structures.svfinv3[i].used = true;
                var_data.sv_structures.svrinv3[j].used = true;
                vref.pstd = true;
                vref.qstd = true;

                let mut tmp = Variation::new();
                tmp.cnt = cnt + rcnt;
                tmp.hicnt = cnt + rcnt;
                tmp.fwd = cnt;
                tmp.rev = rcnt;
                tmp.qmean = qmean + rqmean;
                tmp.pmean = pmean + rpmean;
                tmp.mqmean = Qmean + rQmean;
                tmp.nm = nm + rnm;
                crate::variations::util::adj_cnt(vref, &tmp, None);

                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    bp,
                );
                sv.r#type = "INV".to_string();
                sv.pairs += cnt;
                sv.splits += var_data
                    .soft_clips_3_end
                    .get(&(end + 1))
                    .map_or(0, |s| s.variation.cnt);
                sv.splits += var_data
                    .soft_clips_3_end
                    .get(&(me + 1))
                    .map_or(0, |s| s.variation.cnt);
                sv.clusters += 1;

                if !var_data.ref_coverage.contains_key(&bp) {
                    var_data.ref_coverage.insert(bp, 2 * cnt);
                }

                let mut svfinv5_mut_refs: Vec<&mut Vec<Sclip>> = vec![
                    &mut var_data.sv_structures.svfinv5,
                    &mut var_data.sv_structures.svrinv5,
                ];
                mark_sv(bp, pe, &mut svfinv5_mut_refs, max_read_length, verbose);
            }
        }
    }

    Ok(())
}

/// Find duplication structural variants from discordant pairs only.
/// Ports `StructuralVariantsProcessor.findDUPdisc`.
pub fn find_dup_disc(
    var_data: &mut crate::variations::VariationData,
    bam_path: &std::path::Path,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let minr = config.min_reads as i32;
    let verbose = config.verbose;
    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_maps(&config.adaptor);

    let n_svfdup = var_data.sv_structures.svfdup.len();
    for idx in 0..n_svfdup {
        let (used, vars_count, start, mut end, mstart, mend, qmean, Qmean, pmean, nm) = {
            let dup = &var_data.sv_structures.svfdup[idx];
            if dup.used || dup.variation.cnt < minr + 5 {
                continue;
            }
            if dup.variation.mqmean / (dup.variation.cnt as f64) <= 35.0 {
                continue;
            }
            (
                dup.used,
                dup.variation.cnt,
                dup.start,
                dup.end,
                dup.mstart,
                dup.mend,
                dup.variation.qmean,
                dup.variation.mqmean,
                dup.variation.pmean,
                dup.variation.nm,
            )
        };

        if used {
            continue;
        }

        let mut mlen = end - mstart + var_data.max_read_length / vars_count;
        let mut bp = mstart - (var_data.max_read_length / vars_count) / 2;
        let mut pe = end;

        if !is_loaded(chr, mstart as u32, mend as u32, reference) {
            if reference.get_base(chr, bp as u32).is_none() {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1((bp - 150).max(1) as u32),
                    end: Pos1((bp + 150) as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
            }
            partial_pipeline(
                bam_path,
                var_data,
                reference,
                chr,
                mstart - 200,
                mend + 200,
                config,
                chromosomes_lengths,
            )?;
        }

        let mut cntf = vars_count;
        let mut cntr = vars_count;
        let mut qmeanf = qmean;
        let mut qmeanr = qmean;
        let mut Qmeanf = Qmean;
        let mut Qmeanr = Qmean;
        let mut pmeanf = pmean;
        let mut pmeanr = pmean;
        let mut nmf = nm;
        let mut nmr = nm;

        let has_soft = {
            let dup = &var_data.sv_structures.svfdup[idx];
            !dup.soft.is_empty()
        };

        if has_soft {
            let mut soft_sorted: Vec<(i32, i32)> = {
                let dup = &var_data.sv_structures.svfdup[idx];
                dup.soft.iter().map(|(&k, &v)| (k, v)).collect()
            };
            soft_sorted.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
            if !soft_sorted.is_empty() {
                pe = soft_sorted[0].0;
            }

            let scv_3_clone = match var_data.soft_clips_3_end.get(&pe) {
                Some(s) if !s.used => s.clone(),
                _ => continue,
            };

            let mut current_sclip3 = scv_3_clone;
            cntf = current_sclip3.variation.cnt;
            qmeanf = current_sclip3.variation.qmean;
            Qmeanf = current_sclip3.variation.mqmean;
            pmeanf = current_sclip3.variation.pmean;
            nmf = current_sclip3.variation.nm;

            let seq = crate::variations::util::find_conseq(
                &mut current_sclip3,
                0,
                &adaptor_forward,
                &adaptor_reverse,
            );
            let match_res = find_match(&seq, reference, bp, 1, 17, 3, chr, verbose);
            let mut tbp = match_res.base_position;

            if tbp != 0 && tbp < pe {
                if let Some(scv_mut) = var_data.soft_clips_3_end.get_mut(&pe) {
                    scv_mut.used = true;
                }
                while reference.get_base(chr, (pe - 1) as u32).is_some()
                    && reference.get_base(chr, (tbp - 1) as u32).is_some()
                    && reference.get_base(chr, (pe - 1) as u32)
                        == reference.get_base(chr, (tbp - 1) as u32)
                {
                    tbp -= 1;
                    if tbp != 0 {
                        pe -= 1;
                    }
                }
                mlen = pe - tbp;
                bp = tbp;
                pe -= 1;
                end = pe;

                if let Some(current_sclip5) = var_data.soft_clips_5_end.get(&bp) {
                    cntr = current_sclip5.variation.cnt;
                    qmeanr = current_sclip5.variation.qmean;
                    Qmeanr = current_sclip5.variation.mqmean;
                    pmeanr = current_sclip5.variation.pmean;
                    nmr = current_sclip5.variation.nm;
                }
            }
        }

        let ins5 = crate::util::join_ref(reference, chr, bp, bp + 50 - 1);
        let ins3 = crate::util::join_ref(reference, chr, pe - 50 + 1, pe);
        let ins = format!("{}<dup{}>{}", ins5, mlen - 100, ins3);

        let vref = crate::variations::get_variation(
            &mut var_data.insertion_variants,
            bp,
            &format!("+{}", ins),
        );
        vref.cnt = 0;

        let sv = crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, bp);
        sv.r#type = "DUP".to_string();
        sv.pairs += vars_count;
        let dup_softp = var_data.sv_structures.svfdup[idx].softp;
        sv.splits += var_data
            .soft_clips_3_end
            .get(&dup_softp)
            .map(|s| s.variation.cnt)
            .unwrap_or(0);
        sv.clusters += 1;

        let tcnt = cntr + cntf;
        let mut tmp = Variation::new();
        tmp.cnt = tcnt;
        tmp.hicnt = tcnt;
        tmp.extracnt = tcnt;
        tmp.fwd = cntf;
        tmp.rev = cntr;
        tmp.qmean = qmeanf + qmeanr;
        tmp.pmean = pmeanf + pmeanr;
        tmp.mqmean = Qmeanf + Qmeanr;
        tmp.nm = nmf + nmr;

        let vref_mut = crate::variations::get_variation(
            &mut var_data.insertion_variants,
            bp,
            &format!("+{}", ins),
        );
        crate::variations::util::adj_cnt(vref_mut, &tmp, None);

        var_data.sv_structures.svfdup[idx].used = true;

        if !var_data.ref_coverage.contains_key(&bp) {
            var_data.ref_coverage.insert(bp, tcnt);
        }
        let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
        let end_cov = var_data.ref_coverage.get(&end).copied().unwrap_or(0);
        if bp_cov < end_cov {
            var_data.ref_coverage.insert(bp, end_cov);
        }

        let mut svrdup_mut_refs: Vec<&mut Vec<Sclip>> = vec![&mut var_data.sv_structures.svrdup];
        let (clusters, _) = mark_dup_sv(
            bp,
            pe,
            &mut svrdup_mut_refs,
            var_data.max_read_length,
            verbose,
        );
        {
            let sv =
                crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, bp);
            sv.clusters += clusters;
        }
    }

    // svrdup loop
    let n_svrdup = var_data.sv_structures.svrdup.len();
    for idx in 0..n_svrdup {
        let (used, vars_count, start, end, mstart, mend) = {
            let dup = &var_data.sv_structures.svrdup[idx];
            if dup.used || dup.variation.cnt < minr + 5 {
                continue;
            }
            if dup.variation.mqmean / (dup.variation.cnt as f64) <= 35.0 {
                continue;
            }
            (
                dup.used,
                dup.variation.cnt,
                dup.start,
                dup.end,
                dup.mstart,
                dup.mend,
            )
        };

        if used {
            continue;
        }

        let mut mlen = mend - start + var_data.max_read_length / vars_count;
        let mut bp = start - (var_data.max_read_length / vars_count) / 2;
        let mut pe = mlen + bp - 1;
        let mut tpe = pe;

        if !is_loaded(chr, mstart as u32, mend as u32, reference) {
            if reference.get_base(chr, pe as u32).is_none() {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1((pe - 150).max(1) as u32),
                    end: Pos1((pe + 150) as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
            }
            partial_pipeline(
                bam_path,
                var_data,
                reference,
                chr,
                mstart - 200,
                mend + 200,
                config,
                chromosomes_lengths,
            )?;
        }

        let mut cntf = vars_count;
        let mut cntr = vars_count;
        let mut qmeanf = {
            let dup = &var_data.sv_structures.svrdup[idx];
            dup.variation.qmean
        };
        let mut qmeanr = qmeanf;
        let mut Qmeanf = {
            let dup = &var_data.sv_structures.svrdup[idx];
            dup.variation.mqmean
        };
        let mut Qmeanr = Qmeanf;
        let mut pmeanf = {
            let dup = &var_data.sv_structures.svrdup[idx];
            dup.variation.pmean
        };
        let mut pmeanr = pmeanf;
        let mut nmf = {
            let dup = &var_data.sv_structures.svrdup[idx];
            dup.variation.nm
        };
        let mut nmr = nmf;

        let has_soft = {
            let dup = &var_data.sv_structures.svrdup[idx];
            !dup.soft.is_empty()
        };

        if has_soft {
            let mut soft_sorted: Vec<(i32, i32)> = {
                let dup = &var_data.sv_structures.svrdup[idx];
                dup.soft.iter().map(|(&k, &v)| (k, v)).collect()
            };
            soft_sorted.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
            if !soft_sorted.is_empty() {
                bp = soft_sorted[0].0;
            }

            let scv_5_clone = match var_data.soft_clips_5_end.get(&bp) {
                Some(s) if !s.used => s.clone(),
                _ => continue,
            };

            let mut current_sclip5 = scv_5_clone;
            cntr = current_sclip5.variation.cnt;
            qmeanr = current_sclip5.variation.qmean;
            Qmeanr = current_sclip5.variation.mqmean;
            pmeanr = current_sclip5.variation.pmean;
            nmr = current_sclip5.variation.nm;

            let seq = crate::variations::util::find_conseq(
                &mut current_sclip5,
                0,
                &adaptor_forward,
                &adaptor_reverse,
            );
            let match_res = find_match(&seq, reference, pe, -1, 17, 3, chr, verbose);
            let tbp = match_res.base_position;

            if tbp != 0 && tbp > bp {
                if let Some(scv_mut) = var_data.soft_clips_5_end.get_mut(&bp) {
                    scv_mut.used = true;
                }
                pe = tbp;
                mlen = pe - bp + 1;
                tpe = pe + 1;
                while reference.get_base(chr, tpe as u32).is_some()
                    && reference
                        .get_base(chr, (bp + (tpe - pe - 1)) as u32)
                        .is_some()
                    && reference.get_base(chr, tpe as u32)
                        == reference.get_base(chr, (bp + (tpe - pe - 1)) as u32)
                {
                    tpe += 1;
                }
                if let Some(current_sclip3) = var_data.soft_clips_3_end.get(&tpe) {
                    cntf = current_sclip3.variation.cnt;
                    qmeanf = current_sclip3.variation.qmean;
                    Qmeanf = current_sclip3.variation.mqmean;
                    pmeanf = current_sclip3.variation.pmean;
                    nmf = current_sclip3.variation.nm;
                }
            }
        }

        let ins5 = crate::util::join_ref(reference, chr, bp, bp + 50 - 1);
        let ins3 = crate::util::join_ref(reference, chr, pe - 50 + 1, pe);
        let ins = format!("{}<dup{}>{}", ins5, mlen - 100, ins3);

        let vref = crate::variations::get_variation(
            &mut var_data.insertion_variants,
            bp,
            &format!("+{}", ins),
        );
        vref.cnt = 0;

        {
            let sv =
                crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, bp);
            sv.r#type = "DUP".to_string();
            sv.pairs += vars_count;
            sv.splits += var_data
                .soft_clips_5_end
                .get(&bp)
                .map(|s| s.variation.cnt)
                .unwrap_or(0);
            sv.splits += var_data
                .soft_clips_3_end
                .get(&tpe)
                .map(|s| s.variation.cnt)
                .unwrap_or(0);
            sv.clusters += 1;
        }

        let tcnt = cntr + cntf;
        let mut tmp = Variation::new();
        tmp.cnt = tcnt;
        tmp.hicnt = tcnt;
        tmp.extracnt = tcnt;
        tmp.fwd = cntf;
        tmp.rev = cntr;
        tmp.qmean = qmeanf + qmeanr;
        tmp.pmean = pmeanf + pmeanr;
        tmp.mqmean = Qmeanf + Qmeanr;
        tmp.nm = nmf + nmr;

        let vref_mut = crate::variations::get_variation(
            &mut var_data.insertion_variants,
            bp,
            &format!("+{}", ins),
        );
        crate::variations::util::adj_cnt(vref_mut, &tmp, None);

        var_data.sv_structures.svrdup[idx].used = true;

        if !var_data.ref_coverage.contains_key(&bp) {
            var_data.ref_coverage.insert(bp, tcnt);
        }
        let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
        let mend_cov = var_data.ref_coverage.get(&mend).copied().unwrap_or(0);
        if bp_cov < mend_cov {
            var_data.ref_coverage.insert(bp, mend_cov);
        }

        let mut svfdup_mut_refs: Vec<&mut Vec<Sclip>> = vec![&mut var_data.sv_structures.svfdup];
        let (clusters, _) = mark_dup_sv(
            bp,
            pe,
            &mut svfdup_mut_refs,
            var_data.max_read_length,
            verbose,
        );
        {
            let sv =
                crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, bp);
            sv.clusters += clusters;
        }
    }

    Ok(())
}

/// Run structural variant calling on the accumulated variation data of a tile.
/// Ports `StructuralVariantsProcessor.process`.
pub fn find_all_sv(
    var_data: &mut crate::variations::VariationData,
    bam_path: &std::path::Path,
    reference: &mut Reference,
    chr: &str,
    chr_len: u32,
    pileup_region: &crate::regions::Region,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    if config.no_sv {
        return Ok(());
    }
    filter_all_sv_structures(var_data);

    if verbose_log(config) {
        eprintln!("Start Structural Variants: DEL\n");
    }
    find_del(
        var_data,
        bam_path,
        reference,
        chr,
        config,
        chromosomes_lengths,
    )?;

    if verbose_log(config) {
        eprintln!("Start Structural Variants: INV\n");
    }
    find_inv(
        var_data,
        bam_path,
        reference,
        chr,
        config,
        chromosomes_lengths,
    )?;

    if verbose_log(config) {
        eprintln!("Start Structural Variants\n");
    }
    find_sv(
        var_data,
        reference,
        chr,
        pileup_region.start.0 as i32,
        pileup_region.end.0 as i32,
        config,
        chromosomes_lengths,
    )?;

    if verbose_log(config) {
        eprintln!("Start Structural Variants: DEL discordant pairs only\n");
    }
    find_del_disc(var_data, reference, chr, config, chromosomes_lengths)?;

    if verbose_log(config) {
        eprintln!("Start Structural Variants: INV discordant pairs only\n");
    }
    find_inv_disc(var_data, reference, chr, config, chromosomes_lengths)?;

    if verbose_log(config) {
        eprintln!("Start Structural Variants: DUP discordant pairs only\n");
    }
    find_dup_disc(
        var_data,
        bam_path,
        reference,
        chr,
        config,
        chromosomes_lengths,
    )?;

    Ok(())
}

fn verbose_log(config: &Config) -> bool {
    config.verbose // matches instance().conf.y in Java
}

/// Dynamically query BAM and run pileup/CIGAR parse on a region that is not currently loaded.
/// Ports getMode().partialPipeline(...) in Java.
pub fn partial_pipeline(
    bam_path: &std::path::Path,
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    chr: &str,
    start: i32,
    end: i32,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    // Query BAM records in the modified region
    let query_start = crate::types::Pos0((start.max(1) - 1) as u32);
    let query_end = crate::types::Pos0(end.max(1) as u32);

    // Streaming counterpart of the old query_region-Vec + pileup_into shape:
    // one lazy pass through a single preprocessor into a parser seeded with
    // the taken maps, written back exactly as pileup_into does (fresh
    // SVStructures discarded, discordant seeded 0 — duprate never propagates
    // out of this path either way). Peak memory stays flat in depth.
    let mut feed = crate::bam::query_region_iter(bam_path, chr, query_start, query_end, config)?;

    // Preprocess records
    let mut preprocessor = crate::bam::RecordPreprocessor::new();

    let region = crate::regions::Region {
        chr: chr.to_string(),
        start: Pos1(start.max(1) as u32),
        end: Pos1(end.max(1) as u32),
        gene: String::new(),
        amplicon: None,
    };

    // Java creates a fresh SVStructures in CigarParser(true) for partialPipeline
    // (CigarParser.java:66) and discards it when partialPipeline completes, so the
    // primary region's sv_structures is never mutated by mate-region pileup.
    //
    // G3 (mate-window ownership): snapshot variant-map keys so keys planted
    // by this mate-window re-pileup are attributed in `merge_planted` for the
    // ToVarsBuilder ownership gate. Pre-existing keys are native (the merge
    // only accumulates into them), never planted.
    let before_non: BTreeSet<i32> = var_data.non_insertion_variants.keys().copied().collect();
    let before_ins: BTreeSet<i32> = var_data.insertion_variants.keys().copied().collect();
    let orig_sv_structures = std::mem::take(&mut var_data.sv_structures);
    {
        let mut parser = crate::cigar::parse::CigarParser {
            reference,
            region: &region,
            config,
            svflag: true, // svflag is true for partial pipeline
            non_insertion_variants: std::mem::take(&mut var_data.non_insertion_variants),
            insertion_variants: std::mem::take(&mut var_data.insertion_variants),
            position_to_insertion_count: std::mem::take(&mut var_data.position_to_insertion_count),
            position_to_deletion_count: std::mem::take(&mut var_data.position_to_deletion_count),
            sv_structures: std::mem::take(&mut var_data.sv_structures),
            ref_coverage: std::mem::take(&mut var_data.ref_coverage),
            soft_clips_5_end: std::mem::take(&mut var_data.soft_clips_5_end),
            soft_clips_3_end: std::mem::take(&mut var_data.soft_clips_3_end),
            max_read_length: var_data.max_read_length,
            splice: std::mem::take(&mut var_data.splice),
            mnp: std::mem::take(&mut var_data.mnp),
            splice_count: std::mem::take(&mut var_data.splice_count),
            discordant_count: 0,
            total_reads: 0,
            duplicate_reads: 0,
            scratch_seq: Vec::new(),
            scratch_qual: Vec::new(),
            scratch_cigar_str: String::new(),
            scratch_elements: Vec::new(),
            scratch_s: String::new(),
            scratch_ss: String::new(),
            scratch_splice_key: String::new(),
            scratch_f4_desc: String::new(),
            scratch_f4_qual: String::new(),
            scratch_f4_seg1: String::new(),
            scratch_f4_seg2: String::new(),
        };
        let mut r = crate::bam::SamRecord::new();
        while feed.next_into(&mut r)? {
            if preprocessor.preprocess_record(&r, config) {
                parser.parse_cigar(&r)?;
            }
        }
        parser.total_reads = preprocessor.total_reads as i32;
        parser.duplicate_reads = preprocessor.duplicate_reads as i32;
        var_data.non_insertion_variants = parser.non_insertion_variants;
        var_data.insertion_variants = parser.insertion_variants;
        var_data.position_to_insertion_count = parser.position_to_insertion_count;
        var_data.position_to_deletion_count = parser.position_to_deletion_count;
        var_data.sv_structures = parser.sv_structures;
        var_data.ref_coverage = parser.ref_coverage;
        var_data.soft_clips_5_end = parser.soft_clips_5_end;
        var_data.soft_clips_3_end = parser.soft_clips_3_end;
        var_data.max_read_length = var_data.max_read_length.max(parser.max_read_length);
        var_data.splice = parser.splice;
        var_data.mnp = parser.mnp;
        var_data.splice_count = parser.splice_count;
    }
    var_data.sv_structures = orig_sv_structures;
    // Attribute mate-window-planted keys: positions present now but
    // absent from BOTH maps before the merge. Pre-existing keys are native
    // (the merge only accumulates into them), never planted — even when the
    // merge adds the FIRST entry of the other map at that position, since
    // the ToVarsBuilder ownership gate drops the whole position.
    for &k in var_data.non_insertion_variants.keys() {
        if !before_non.contains(&k) && !before_ins.contains(&k) {
            var_data.merge_planted.insert(k);
        }
    }
    for &k in var_data.insertion_variants.keys() {
        if !before_non.contains(&k) && !before_ins.contains(&k) {
            var_data.merge_planted.insert(k);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression test for W1.4 (NA24385 chr20:20201754 duplicate INV): Java's
    /// findINVsub returns after the FIRST successful INV per call (`return
    /// vref;`, StructuralVariantsProcessor.java:675). Rust kept scanning the
    /// cluster list after a success and built a second INV from a later
    /// cluster Java never reached — the record then carried SV info `3-2-1` /
    /// cov 92 instead of the nested-loop `3-0-0` / cov 109.
    #[test]
    fn test_find_inv_sub_returns_after_first_success() {
        let chr = "chrT";
        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });

        // Cluster A consensus: revcomp (all-G) matches reference ending at 400.
        let mut ref_seq = vec![b'A'; 1000];
        for g in &mut ref_seq[380..400] {
            *g = b'G';
        }
        // Seed window (i=3, 17bp of all-G) -> unique position; bp = P+len-i-1.
        reference.seed.insert(vec![b'G'; 17], vec![384]); // bp = 384+20-3-1 = 400
        reference
            .reference_sequences
            .insert(chr.to_string(), ref_seq);

        let mut var_data = crate::variations::VariationData::default();
        var_data.max_read_length = 100;
        var_data.ref_coverage.insert(159, 77);

        // Cluster A: disc-pair cluster with empty soft map -> scan branch from
        // end=150 finds SC3@160 first; |bp-mend|=50 < 1.5*100 so it breaks.
        let mut inv_a = Sclip::default();
        inv_a.variation.cnt = 2;
        inv_a.start = 100;
        inv_a.end = 150;
        inv_a.mstart = 300;
        inv_a.mend = 350;
        inv_a.mlen = -200;
        var_data.sv_structures.svfinv5.push(inv_a);

        // Cluster B: would ALSO succeed if the loop kept going after A
        // (scan from 350 hits SC3@520, bp=900, len=381, 381/260 < 1.5).
        let mut inv_b = Sclip::default();
        inv_b.variation.cnt = 4;
        inv_b.start = 300;
        inv_b.end = 350;
        inv_b.mstart = 850;
        inv_b.mend = 950;
        inv_b.mlen = -260;
        var_data.sv_structures.svfinv5.push(inv_b);

        // SC3@160: consensus all-C -> find_match_rev revcomps it to all-G,
        // which matches the engineered reference region ending at 400.
        let mut scv_a = Sclip::default();
        scv_a.variation.cnt = 3;
        scv_a.sequence = Some(vec![b'C'; 20]);
        var_data.soft_clips_3_end.insert(160, scv_a);

        // SC3@520: consensus all-T revcomp is all-A == background ref;
        // seed 17xA unique at 884 -> bp = 884+20-3-1 = 900.
        let mut scv_b = Sclip::default();
        scv_b.variation.cnt = 5;
        scv_b.sequence = Some(vec![b'T'; 20]);
        var_data.soft_clips_3_end.insert(520, scv_b);
        reference.seed.insert(vec![b'A'; 17], vec![884]);

        let mut chromosomes_lengths = HashMap::new();
        chromosomes_lengths.insert(chr.to_string(), 1000usize);
        let config = crate::config::Config::default();

        find_inv_sub(
            0,
            1,
            Side::Five,
            &mut var_data,
            std::path::Path::new(""),
            &mut reference,
            chr,
            &config,
            &chromosomes_lengths,
        )
        .expect("find_inv_sub failed");

        // First cluster processed and consumed its soft clip.
        assert!(var_data.sv_structures.svfinv5[0].used);
        assert!(var_data.soft_clips_3_end.get(&160).unwrap().used);

        // The INV was registered at softp=160 with Java's field assignments:
        // splits += scv.varsCount(3), pairs += inv.varsCount(2), clusters++.
        let entry = var_data.non_insertion_variants.get(&160).unwrap();
        assert!(entry.map.keys().any(|k| k.starts_with("-241^")));
        let sv = entry.sv.as_ref().unwrap();
        assert_eq!(sv.r#type, "INV");
        assert_eq!((sv.splits, sv.pairs, sv.clusters), (3, 2, 1));
        // Coverage propagated from softp-1 (Java refCoverage.put).
        assert_eq!(var_data.ref_coverage.get(&160), Some(&77));

        // THE REGRESSION: the second cluster must be untouched — Java's early
        // return leaves it for nothing else. Before the fix it was consumed
        // here (used=true) and emitted a second INV at softp=520.
        assert!(!var_data.sv_structures.svfinv5[1].used);
        assert!(!var_data.soft_clips_3_end.get(&520).unwrap().used);
        assert!(!var_data.non_insertion_variants.contains_key(&520));
    }

    /// Regression test for W3.3a (Cat 8 latent item 1): Java findINVsub's
    /// success block performs exactly three mutations — inv.used
    /// (StructuralVariantsProcessor.java:644), scv.used on the LIVE fetched
    /// entry (:663) and the embedded realigndel (:665-667) — then returns.
    /// It never marks the soft clip under the post-shift/swapped softp and
    /// never sweeps markSV over the four INV lists (no markSV call exists in
    /// realigndel; the realignlgdel/findDEL/findDELdisc/findINVdisc sites run
    /// outside findINV or use their own coordinates). Rust used to do both,
    /// suppressing records Java keeps. Scenario: side=3/dir=1 detection swaps
    /// bp and softp (Java :605-608), so the old code marked a decoy clip at
    /// the swapped-in coordinate and swept four decoy clusters.
    #[test]
    fn test_find_inv_sub_success_marks_only_java_equivalents() {
        let chr = "chrT";
        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        // Matched region ending at bp=200: seed [G;17] -> [184],
        // bp = 184 + 20 - 3 - 1 = 200.
        let mut ref_seq = vec![b'A'; 1000];
        for g in &mut ref_seq[180..200] {
            *g = b'G';
        }
        reference.seed.insert(vec![b'G'; 17], vec![184]);
        reference
            .reference_sequences
            .insert(chr.to_string(), ref_seq);

        let mut var_data = crate::variations::VariationData::default();
        var_data.max_read_length = 100;
        var_data.ref_coverage.insert(200, 77);

        // INV disc-pair cluster (svfinv3, dir=1, side=3) with empty soft map
        // -> scan branch from end=500 finds SC3@520 (i=20); |bp-mend|=140 <
        // 150 breaks. Side::Three/dir=1: bp 200->201, softp 520->519, swap ->
        // bp=519, softp=201; gap=318, |mlen|=260 -> 318/260 < 1.5 accepted.
        let mut inv = Sclip::default();
        inv.variation.cnt = 3;
        inv.start = 100;
        inv.end = 500;
        inv.mstart = 300;
        inv.mend = 340;
        inv.mlen = -260;
        var_data.sv_structures.svfinv3.push(inv);

        // SC3@520: consensus all-C -> revcomp all-G matches the engineered
        // region ending at 200.
        let mut scv = Sclip::default();
        scv.variation.cnt = 5;
        scv.sequence = Some(vec![b'C'; 20]);
        var_data.soft_clips_3_end.insert(520, scv);

        // Decoy clip at the POST-SWAP softp (201): the removed marking
        // flagged it used; Java leaves it usable.
        let mut decoy_clip = Sclip::default();
        decoy_clip.variation.cnt = 9;
        var_data.soft_clips_3_end.insert(201, decoy_clip);

        // Decoy INV clusters overlapping (bp=519, softp=201): the removed
        // mark_sv x4 sweep flagged them used; Java leaves them available.
        for list_idx in 0..4 {
            let mut decoy = Sclip::default();
            decoy.variation.cnt = 7;
            decoy.start = 120;
            decoy.end = 160;
            decoy.mstart = 600;
            decoy.mend = 640;
            decoy.mlen = -200;
            match list_idx {
                0 => var_data.sv_structures.svfinv5.push(decoy),
                1 => var_data.sv_structures.svrinv5.push(decoy),
                2 => var_data.sv_structures.svfinv3.push(decoy),
                _ => var_data.sv_structures.svrinv3.push(decoy),
            }
        }

        let mut chromosomes_lengths = HashMap::new();
        chromosomes_lengths.insert(chr.to_string(), 1000usize);
        let config = crate::config::Config::default();

        find_inv_sub(
            2,
            1,
            Side::Three,
            &mut var_data,
            std::path::Path::new(""),
            &mut reference,
            chr,
            &config,
            &chromosomes_lengths,
        )
        .expect("find_inv_sub failed");

        // The INV was detected and registered under softp=201 (len=319).
        let entry = var_data.non_insertion_variants.get(&201).unwrap();
        assert!(entry.map.keys().any(|k| k.starts_with("-319^")));

        // Java-equivalent markings happened: the fetched live clip (@520)
        // and the processed cluster itself.
        assert!(var_data.soft_clips_3_end.get(&520).unwrap().used);
        assert!(var_data.sv_structures.svfinv3[0].used);

        // THE REMOVED MARKINGS: decoy clip at the swapped coordinate and the
        // four decoy clusters stay untouched (available for later detection,
        // as in Java).
        assert!(!var_data.soft_clips_3_end.get(&201).unwrap().used);
        assert!(!var_data.sv_structures.svfinv5[0].used);
        assert!(!var_data.sv_structures.svrinv5[0].used);
        // svfinv3[1] is the decoy; svfinv3[0] was legitimately consumed.
        assert!(!var_data.sv_structures.svfinv3[1].used);
        assert!(!var_data.sv_structures.svrinv3[0].used);
    }

    /// Regression test for W3.3b (Cat 8 latent item 2): Java findINVsub has no
    /// SOFTP2SV-style gate — it only skips on `scv.used`
    /// (StructuralVariantsProcessor.java:534, scan branch :561). The
    /// best-cluster check exists solely in findsv (:705/:835). A cluster whose
    /// highest-cnt anchor at `softp` is used (e.g. marked by findDEL's markSV)
    /// while its own clip entry is unused MUST still be processed. The removed
    /// Rust gate dropped such records entirely.
    #[test]
    fn test_find_inv_sub_ignores_used_anchor_cluster() {
        let chr = "chrT";
        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        // Matched region ending at bp=400: seed [G;17] -> [384],
        // bp = 384 + 20 - 3 - 1 = 400.
        let mut ref_seq = vec![b'A'; 1000];
        for g in &mut ref_seq[380..400] {
            *g = b'G';
        }
        reference.seed.insert(vec![b'G'; 17], vec![384]);
        reference
            .reference_sequences
            .insert(chr.to_string(), ref_seq);

        let mut var_data = crate::variations::VariationData::default();
        var_data.max_read_length = 100;
        var_data.ref_coverage.insert(159, 77);

        // INV cluster anchored via its soft map at softp=160.
        let mut inv = Sclip::default();
        inv.variation.cnt = 3;
        inv.start = 100;
        inv.end = 150;
        inv.mstart = 300;
        inv.mend = 350;
        inv.mlen = -200;
        inv.soft.insert(160, 2);
        var_data.sv_structures.svfinv5.push(inv);

        // Used anchor cluster at the same softp (as if findDEL's markSV had
        // flagged it): best cnt at 160 is 10 > 3 and used=true.
        let mut anchor = Sclip::default();
        anchor.variation.cnt = 10;
        anchor.softp = 160;
        anchor.used = true;
        var_data.sv_structures.svfdel.push(anchor);

        // SC3@160 is NOT used: Java's scv.used gate (:534) lets processing
        // proceed. Consensus all-C -> revcomp all-G matches region @400.
        let mut scv = Sclip::default();
        scv.variation.cnt = 3;
        scv.sequence = Some(vec![b'C'; 20]);
        var_data.soft_clips_3_end.insert(160, scv);

        let mut chromosomes_lengths = HashMap::new();
        chromosomes_lengths.insert(chr.to_string(), 1000usize);
        let config = crate::config::Config::default();

        find_inv_sub(
            0,
            1,
            Side::Five,
            &mut var_data,
            std::path::Path::new(""),
            &mut reference,
            chr,
            &config,
            &chromosomes_lengths,
        )
        .expect("find_inv_sub failed");

        // The INV was detected despite the used anchor cluster (gap=240,
        // |mlen|=200 -> 240/200 < 1.5 accepted): record under softp=160,
        // len=241, live clip + cluster consumed exactly as Java does.
        let entry = var_data.non_insertion_variants.get(&160).unwrap();
        assert!(entry.map.keys().any(|k| k.starts_with("-241^")));
        assert!(var_data.soft_clips_3_end.get(&160).unwrap().used);
        assert!(var_data.sv_structures.svfinv5[0].used);
    }

    /// Builds a side=5/dir=1 scan-branch scenario whose matched position is
    /// always bp=400 (seed [G;17] -> [384], revcomp all-C consensus), with the
    /// soft clip at `clip_pos` so that gap = bp - softp = 400 - clip_pos.
    fn build_gate_scenario(
        clip_pos: i32,
        mlen: i32,
        mend: i32,
    ) -> (
        crate::variations::VariationData,
        Reference,
        HashMap<String, usize>,
    ) {
        let chr = "chrT".to_string();
        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.clone(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        let mut ref_seq = vec![b'A'; 1000];
        for g in &mut ref_seq[380..400] {
            *g = b'G';
        }
        reference.seed.insert(vec![b'G'; 17], vec![384]);
        reference.reference_sequences.insert(chr.clone(), ref_seq);

        let mut var_data = crate::variations::VariationData::default();
        var_data.max_read_length = 100;
        let mut inv = Sclip::default();
        inv.variation.cnt = 3;
        inv.start = clip_pos - 60;
        inv.end = clip_pos - 1; // scan finds SC3@clip_pos at i=1
        inv.mstart = 300;
        inv.mend = mend; // |400-mend| < 150 breaks the scan
        inv.mlen = mlen;
        var_data.sv_structures.svfinv5.push(inv);
        let mut scv = Sclip::default();
        scv.variation.cnt = 3;
        scv.sequence = Some(vec![b'C'; 20]);
        var_data.soft_clips_3_end.insert(clip_pos, scv);
        let mut chromosomes_lengths = HashMap::new();
        chromosomes_lengths.insert(chr, 1000usize);
        (var_data, reference, chromosomes_lengths)
    }

    /// Regression test for W3.3c (Cat 8 latent item 3): Java gates on
    /// `bp - softp > 150 && (bp-softp)/|mlen| < 1.5`
    /// (StructuralVariantsProcessor.java:627) BEFORE computing
    /// `len = bp - softp + 1` (:628). Rust gated on `len`, which (a) accepted
    /// a gap of exactly 150 and (b) divided (gap+1) by |mlen| — rejecting e.g.
    /// gap=299/mlen=200 where Java computes 1.495 < 1.5 and accepts.
    #[test]
    fn test_find_inv_sub_success_gate_uses_java_gap_semantics() {
        let config = crate::config::Config::default();
        // (clip_pos => gap=400-clip_pos, mlen, mend, java_accepts, gt_prefix)
        let cases: [(i32, i32, i32, bool, &str); 4] = [
            (250, -200, 350, false, ""),     // gap=150: Java rejects (>150 strict)
            (249, -200, 350, true, "-152^"), // gap=151: both accept
            (101, -200, 450, true, "-300^"), // gap=299/200=1.495: len-gate rejected, Java accepts
            (100, -200, 450, false, ""),     // gap=300/200=1.5: both reject (<1.5 strict)
        ];
        for (clip_pos, mlen, mend, accepted, prefix) in cases {
            let (mut var_data, mut reference, chromosomes_lengths) =
                build_gate_scenario(clip_pos, mlen, mend);
            find_inv_sub(
                0,
                1,
                Side::Five,
                &mut var_data,
                std::path::Path::new(""),
                &mut reference,
                "chrT",
                &config,
                &chromosomes_lengths,
            )
            .expect("find_inv_sub failed");
            if accepted {
                let entry = var_data.non_insertion_variants.get(&clip_pos).unwrap();
                assert!(
                    entry.map.keys().any(|k| k.starts_with(prefix)),
                    "clip {clip_pos}: expected {prefix} record"
                );
                assert!(var_data.soft_clips_3_end.get(&clip_pos).unwrap().used);
            } else {
                assert!(
                    !var_data.non_insertion_variants.contains_key(&clip_pos),
                    "clip {clip_pos}: record must be gated off"
                );
                assert!(!var_data.soft_clips_3_end.get(&clip_pos).unwrap().used);
                assert!(!var_data.sv_structures.svfinv5[0].used);
            }
        }
    }

    /// Regression test for W1.2 (Colo829 INVSV/LGDEL3): the INV description tail
    /// must be built from findMatchRev's matched suffix (`extra`) — appended
    /// as-is when dir==-1, revcomp'd+prepended when dir==1 — never from
    /// revcomp(consensus). Values from Colo829.chr20.INVSV @46843164 where the
    /// wrong tail suppressed realigndel's T@47214669 absorption (cnt 2 vs 7).
    #[test]
    fn test_apply_inv_extra_uses_matched_suffix_not_consensus() {
        // Java trace: consensus ATACTGTTCCGGAC, EXTRA: T, dir=-1
        let ins = "ins5<inv371406>ins3".to_string();
        let out = apply_inv_extra(ins.clone(), -1, b"T");
        assert_eq!(out, "ins5<inv371406>ins3T");
        // The old bug: appending revcomp(consensus) instead of extra
        assert_ne!(out, format!("{}GTCCGGAACAGTAT", ins));

        // dir==-1 with empty extra leaves ins unchanged (Java !extra.isEmpty())
        assert_eq!(apply_inv_extra(ins.clone(), -1, b""), ins);

        // dir==1 prepends revcomp(extra)
        let out1 = apply_inv_extra("MID".to_string(), 1, b"AT");
        assert_eq!(out1, "ATMID");
        assert_eq!(apply_inv_extra("MID".to_string(), -1, b"AT"), "MIDAT");
    }

    #[test]
    fn test_is_overlap() {
        assert!(is_overlap(10, 20, 15, 25, 100));
        assert!(!is_overlap(10, 20, 30, 40, 100)); // too far
        assert!(!is_overlap(10, 20, 21, 30, 5));
    }

    #[test]
    fn test_find_match() {
        let mut reference = Reference::default();
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

        // When dir == 1
        let match_res = find_match(seq, &reference, 5, 1, 17, 1, "chr1", false);
        assert_eq!(match_res.base_position, 10);
        assert!(match_res.matched_sequence.is_empty());

        // When dir == -1
        let mut seq_rev = seq.to_vec();
        seq_rev.reverse();
        let match_res_rev = find_match(&seq_rev, &reference, 5, -1, 17, 1, "chr1", false);
        assert_eq!(match_res_rev.base_position, 26);
        assert!(match_res_rev.matched_sequence.is_empty());
    }

    #[test]
    fn test_check_cluster() {
        let mut mates = vec![
            Mate {
                mate_start: 100,
                mate_end: 250,
                mate_len: 150,
                start: 10,
                end: 110,
                pmean: 50.0,
                qmean: 30.0,
                mqmean: 40.0,
                nm: 0.0,
            },
            Mate {
                mate_start: 120,
                mate_end: 270,
                mate_len: 150,
                start: 12,
                end: 112,
                pmean: 52.0,
                qmean: 32.0,
                mqmean: 42.0,
                nm: 1.0,
            },
        ];

        let cluster = check_cluster(&mut mates, 100);
        assert_eq!(cluster.cnt, 2);
        assert_eq!(cluster.mate.mate_start, 100);
        assert_eq!(cluster.mate.mate_end, 270);
        assert_eq!(cluster.mate.mate_len, 150);
        assert_eq!(cluster.mate.start, 10);
        assert_eq!(cluster.mate.end, 112);
        assert_eq!(cluster.mate.pmean, 102.0);
        assert_eq!(cluster.mate.qmean, 62.0);
        assert_eq!(cluster.mate.mqmean, 82.0);
        assert_eq!(cluster.mate.nm, 1.0);
    }

    #[test]
    fn test_find_del_retries_with_seed2_on_short_softclips() {
        // ports StructuralVariantsProcessor.java:251-262, 337-354, 412-423:
        // findDEL gates on `seq.length() < Configuration.SEED_2` (12) and
        // retries `findMatch(seq, reference, softp, dir, Configuration.SEED_2, 0)`
        // across all three retry branches when the 17-mer (SEED_1) lookup returns bp == 0.
        let chr = "chrT";
        let seed12 = b"ACGTACGTACGT"; // 12-mer (SEED_2), < 17 bp (SEED_1)
        let clip14 = b"ACGTACGTACGTCC"; // 14 bp: >= SEED_2 (12) but < SEED_1 (17)

        // 1. svfdel scan branch (softp == 0, scans soft_clips_3_end, dir = 1, :251-262)
        {
            let mut reference = Reference::default();
            reference
                .loaded_regions
                .push(crate::reference::LoadedRegion {
                    chr: chr.to_string(),
                    sequence_start: 1,
                    sequence_end: 1000,
                });
            let mut ref_seq = vec![b'A'; 1000];
            ref_seq[399..413].copy_from_slice(clip14); // 1-based 400..=413
            reference.seed.insert(seed12.to_vec(), vec![400]);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq);

            let mut var_data = crate::variations::VariationData::default();
            var_data.max_read_length = 100;
            let mut del = Sclip::default();
            del.variation.cnt = 4;
            del.start = 100;
            del.end = 150;
            del.mstart = 380;
            del.mend = 450;
            var_data.sv_structures.svfdel.push(del);

            let mut sc3 = Sclip::default();
            sc3.variation.cnt = 3;
            sc3.sequence = Some(clip14.to_vec());
            var_data.soft_clips_3_end.insert(155, sc3); // i >= end - 3 && i - end < 300

            let mut chromosomes_lengths = HashMap::new();
            chromosomes_lengths.insert(chr.to_string(), 1000usize);
            let config = crate::config::Config::default();

            find_del(
                &mut var_data,
                std::path::Path::new(""),
                &mut reference,
                chr,
                &config,
                &chromosomes_lengths,
            )
            .unwrap();

            // bp = 400 -> bp -= 1 (399), dellen = 399 - 155 + 1 = 245 at softp = 155
            let entry = var_data
                .non_insertion_variants
                .get(&155)
                .expect("svfdel scan branch must call DEL via SEED_2=12 retry");
            assert!(entry.map.contains_key("-245"));
        }

        // 2. svrdel anchored branch (softp != 0, soft_clips_5_end, dir = -1, :337-354)
        {
            let mut reference = Reference::default();
            reference
                .loaded_regions
                .push(crate::reference::LoadedRegion {
                    chr: chr.to_string(),
                    sequence_start: 1,
                    sequence_end: 1000,
                });
            let mut ref_seq = vec![b'A'; 1000];
            ref_seq[186..200].copy_from_slice(clip14); // 1-based 187..=200
            reference.seed.insert(seed12.to_vec(), vec![187]);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq);

            let mut var_data = crate::variations::VariationData::default();
            var_data.max_read_length = 100;
            let mut del = Sclip::default();
            del.variation.cnt = 4;
            del.start = 400;
            del.end = 450;
            del.mstart = 150;
            del.mend = 210;
            del.soft.insert(400, 3); // softp_cand = 400
            var_data.sv_structures.svrdel.push(del);

            let mut rev_clip = clip14.to_vec();
            rev_clip.reverse(); // 5' softclip stores reversed sequence
            let mut sc5 = Sclip::default();
            sc5.variation.cnt = 3;
            sc5.sequence = Some(rev_clip);
            var_data.soft_clips_5_end.insert(400, sc5);

            let mut chromosomes_lengths = HashMap::new();
            chromosomes_lengths.insert(chr.to_string(), 1000usize);
            let config = crate::config::Config::default();

            find_del(
                &mut var_data,
                std::path::Path::new(""),
                &mut reference,
                chr,
                &config,
                &chromosomes_lengths,
            )
            .unwrap();

            // find_match(dir=-1): first_seed=187, i=0, len=14 -> bp = 187 + 14 - 0 - 1 = 200
            // bp += 1 -> 201, softp = 399, dellen = 399 - 201 + 1 = 199 at bp = 201
            let entry = var_data
                .non_insertion_variants
                .get(&201)
                .expect("svrdel anchored branch must call DEL via SEED_2=12 retry");
            assert!(entry.map.contains_key("-199"));
        }

        // 3. svrdel scan branch (softp == 0, scans soft_clips_5_end, dir = -1, :412-423)
        {
            let mut reference = Reference::default();
            reference
                .loaded_regions
                .push(crate::reference::LoadedRegion {
                    chr: chr.to_string(),
                    sequence_start: 1,
                    sequence_end: 1000,
                });
            let mut ref_seq = vec![b'A'; 1000];
            ref_seq[186..200].copy_from_slice(clip14); // 1-based 187..=200
            reference.seed.insert(seed12.to_vec(), vec![187]);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq);

            let mut var_data = crate::variations::VariationData::default();
            var_data.max_read_length = 100;
            let mut del = Sclip::default();
            del.variation.cnt = 4;
            del.start = 400;
            del.end = 450;
            del.mstart = 150;
            del.mend = 210;
            // del.soft left empty so softp_cand == 0 -> enters svrdel scan branch
            var_data.sv_structures.svrdel.push(del);

            let mut rev_clip = clip14.to_vec();
            rev_clip.reverse(); // 5' softclip stores reversed sequence
            let mut sc5 = Sclip::default();
            sc5.variation.cnt = 3;
            sc5.sequence = Some(rev_clip);
            var_data.soft_clips_5_end.insert(400, sc5); // i <= start + 3 && start - i < 300

            let mut chromosomes_lengths = HashMap::new();
            chromosomes_lengths.insert(chr.to_string(), 1000usize);
            let config = crate::config::Config::default();

            find_del(
                &mut var_data,
                std::path::Path::new(""),
                &mut reference,
                chr,
                &config,
                &chromosomes_lengths,
            )
            .unwrap();

            let entry = var_data
                .non_insertion_variants
                .get(&201)
                .expect("svrdel scan branch must call DEL via SEED_2=12 retry");
            assert!(entry.map.contains_key("-199"));
        }
    }

    #[test]
    fn test_is_best_sv_used_keeps_forward_before_reverse_tie() {
        // Equal-count tie at the same softp: the forward candidate is scanned
        // first, so an unused forward entry beats a used reverse entry even
        // when the reverse key sorts earlier.
        use crate::variations::VariationData;
        let mut var_data = VariationData::default();
        let mut fwd = crate::variations::Sclip::new();
        fwd.softp = 100;
        fwd.variation.cnt = 5;
        fwd.used = false;
        let mut rev = crate::variations::Sclip::new();
        rev.softp = 100;
        rev.variation.cnt = 5;
        rev.used = true;
        var_data
            .sv_structures
            .svffus
            .insert("chr2".to_string(), vec![fwd]);
        var_data
            .sv_structures
            .svrfus
            .insert("chr1".to_string(), vec![rev]);
        assert!(!is_best_sv_used(&var_data, 100));
    }

    #[test]
    fn test_find_sv_5prime_deletion_uses_post_decrement_p5_plus_1_ref_coverage() {
        // Bug G (`StructuralVariantsProcessor.java:732,748`):
        // In the 5' candidate deletion branch of `find_sv`, Java decrements `p5--` (`p5_adj = p5 - 1`)
        // at line 732 before querying `refCoverage.get(p5 + 1)` at line 748, so `p5 + 1` in Java
        // is `p5_adj + 1` (the original `p5` before decrement).
        use crate::config::Config;
        use crate::reference::{LoadedRegion, Reference};
        use crate::variations::{Sclip, VariationData};
        use std::collections::HashMap;

        let chr = "chr1";
        let mut bases = vec![b'A'; 600];
        let seed_5p = b"CGTACGTACGTACGTAC"; // 17 bp
                                            // Place seed_5p at 1-based positions 184..=200 so a 5' soft-clip at p5 = 400
                                            // matches backwards to bp = 200 (`bp_adj = 201`, `p5_adj = 399`, `dellen = 199`).
        bases[183..200].copy_from_slice(seed_5p);

        let mut reference = Reference::default();
        reference.loaded_regions.push(LoadedRegion {
            chr: chr.to_string(),
            sequence_start: 1,
            sequence_end: 600,
        });
        reference.seed.insert(seed_5p.to_vec(), vec![184]);
        reference
            .reference_sequences
            .insert(chr.to_string(), bases.clone());
        reference.region_sequences.push(bases);

        let mut var_data = VariationData::default();
        var_data.max_read_length = 100;

        let mut rev_clip = seed_5p.to_vec();
        rev_clip.reverse();
        let mut sc5 = Sclip::default();
        sc5.variation.cnt = 3;
        sc5.variation.pmean = 60.0;
        sc5.variation.qmean = 90.0;
        sc5.sequence = Some(rev_clip);
        var_data.soft_clips_5_end.insert(400, sc5);

        let mut cluster = Sclip::default();
        cluster.start = 150;
        cluster.end = 200;
        cluster.mstart = 400;
        cluster.mend = 450;
        cluster.softp = 400;
        cluster.variation.cnt = 2;
        cluster.variation.qmean = 60.0;
        cluster.variation.pmean = 40.0;
        var_data.sv_structures.svfdel.push(cluster);

        // Set ref_coverage at original `p5` (400 == `p5_adj + 1`) to 77, while leaving `p5 + 1` (401) at 0.
        var_data.ref_coverage.insert(400, 77);

        let mut config = Config::default();
        config.min_reads = 2;
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert(chr.to_string(), 600usize);

        find_sv(
            &mut var_data,
            &mut reference,
            chr,
            1,
            600,
            &config,
            &chr_lengths,
        )
        .unwrap();

        assert_eq!(
            var_data.ref_coverage.get(&201).copied(),
            Some(77),
            "find_sv must promote ref_coverage[bp_adj] (201) to ref_coverage[p5_adj + 1] (400)"
        );
    }
}
