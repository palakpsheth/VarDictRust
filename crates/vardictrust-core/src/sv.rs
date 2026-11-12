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
use std::collections::{BTreeMap, HashMap};

#[inline]
fn midpoint_i32(a: i32, b: i32) -> i32 {
    ((i64::from(a) + i64::from(b)) / 2) as i32
}

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
        if i64::from(mate_m.mate_start) - i64::from(current_cluster_end) > i64::from(limit) {
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
        current_cluster.cnt = current_cluster.cnt.wrapping_add(1);
        // VariationRealigner.java:241 — Java 32-bit signed `int +=` wraps on overflow.
        current_cluster.mate.mate_len = current_cluster.mate.mate_len.wrapping_add(mate_m.mate_len);

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
                mate_len: first_cluster.mate.mate_len.wrapping_div(first_cluster.cnt),
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

        sv.softp = best_soft_pos(&sv.soft).unwrap_or(0);
    }
}

#[inline]
fn best_soft_pos(soft: &indexmap::IndexMap<i32, i32>) -> Option<i32> {
    let mut best_pos = None;
    let mut best_count = i32::MIN;
    for (&k, &v) in soft {
        if v > best_count {
            best_count = v;
            best_pos = Some(k);
        }
    }
    best_pos
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

    // Deterministic forward-before-reverse traversal in Java HashMap<String, List<Sclip>>
    // key order (bucket order by String.hashCode, preserving first-insertion order within
    // a bucket), with svffus fully scanned before svrfus (VariationRealigner.java:143-164).
    let scan_fus = |fus_map: &indexmap::IndexMap<String, Vec<Sclip>>,
                    best_cnt: &mut i32,
                    best_used: &mut bool| {
        if fus_map.len() <= 1 {
            for sv_list in fus_map.values() {
                for sv in sv_list {
                    if sv.softp == softp && sv.variation.cnt > *best_cnt {
                        *best_cnt = sv.variation.cnt;
                        *best_used = sv.used;
                    }
                }
            }
            return;
        }
        if !fus_map.values().any(|sv_list| {
            sv_list
                .iter()
                .any(|sv| sv.softp == softp && sv.variation.cnt > *best_cnt)
        }) {
            return;
        }
        for key in crate::javamap::java_hashmap_key_order(fus_map.keys()) {
            if let Some(sv_list) = fus_map.get(key) {
                for sv in sv_list {
                    if sv.softp == softp && sv.variation.cnt > *best_cnt {
                        *best_cnt = sv.variation.cnt;
                        *best_used = sv.used;
                    }
                }
            }
        }
    };
    scan_fus(
        &var_data.sv_structures.svffus,
        &mut best_cnt,
        &mut best_used,
    );
    scan_fus(
        &var_data.sv_structures.svrfus,
        &mut best_cnt,
        &mut best_used,
    );

    best_cnt >= 0 && best_used
}

/// Syncs `Sclip.used` (set on poly-A/T or low-complexity sequences in
/// `VariationUtils.java:125-134`) and cached `Sclip.sequence` (`:148`) from a
/// temporary `Sclip` clone back to the live `soft_clips_{5,3}_end` map entry.
#[inline]
fn sync_sclip_after_conseq(map: &mut BTreeMap<i32, Sclip>, pos: i32, scv: &Sclip) {
    if let Some(live) = map.get_mut(&pos) {
        live.used = scv.used;
        live.sequence.clone_from(&scv.sequence);
    }
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

#[inline]
fn is_has_and_not_equals_i64(reference: &Reference, chr: &str, pos64: i64, query_base: u8) -> bool {
    i32::try_from(pos64)
        .ok()
        .and_then(|p| u32::try_from(p).ok())
        .is_some_and(|pos| is_has_and_not_equals(reference, chr, pos, query_base))
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
        crate::verbose_eprintln!(
            "    Working Match {} {} {} SEED: {}",
            position,
            String::from_utf8_lossy(&seq_vec),
            dir,
            seed_len
        );
    }

    let mut extra_slice: &[u8] = &[];
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
            let first_seed64 = i64::from(*first_seed_u32);
            let mut bp64 = if dir == 1 {
                first_seed64 - (i as i64)
            } else {
                first_seed64 + (seq_vec.len() as i64) - (i as i64) - 1
            };

            if i32::try_from(bp64).ok().is_some_and(|bp| {
                crate::realign::is_match_ref(&seq_vec, reference, chr, bp, dir, mm_limit)
            }) {
                let mut extra = extra_slice.to_vec();
                let mut mm = if dir == -1 { -1 } else { 0 };
                loop {
                    if let Some(ch) = crate::util::char_at(&seq_vec, mm) {
                        if is_has_and_not_equals_i64(reference, chr, bp64, ch) {
                            extra.push(ch);
                            bp64 += i64::from(dir);
                            mm += dir as isize;
                            continue;
                        }
                    }
                    break;
                }
                let Ok(bp_out) = i32::try_from(bp64) else {
                    continue;
                };
                if !extra.is_empty() && dir == -1 {
                    extra.reverse();
                }
                if verbose {
                    crate::verbose_eprintln!(
                        "      Found SV BP: {} BP: {} SEEDpos {} {} {} {} {} extra: {}",
                        dir,
                        bp_out,
                        first_seed64,
                        position,
                        String::from_utf8_lossy(seed),
                        i,
                        String::from_utf8_lossy(&seq_vec),
                        String::from_utf8_lossy(&extra)
                    );
                }
                return Match {
                    base_position: bp_out,
                    matched_sequence: extra,
                };
            } else {
                // For complex indels, allowing some mismatches at the end up to 15bp or 20% length
                let mut sseq: &[u8] = &seq_vec;
                let mut eqcnt = 0;
                let mut bp_mut64 = bp64;
                for ii in 1..=15 {
                    bp_mut64 += i64::from(dir);
                    if dir == 1 {
                        sseq = sseq.get(1..).unwrap_or(&[]);
                        let Some(b0) = crate::util::char_at(sseq, 0) else {
                            break;
                        };
                        if is_has_and_not_equals_i64(reference, chr, bp_mut64, b0) {
                            continue;
                        }
                        eqcnt += 1;
                        if crate::util::char_at(sseq, 1).is_some_and(|b1| {
                            is_has_and_not_equals_i64(reference, chr, bp_mut64 + 1, b1)
                        }) {
                            continue;
                        }
                        extra_slice = seq_vec.get(..ii).unwrap_or(&seq_vec);
                    } else {
                        sseq = sseq.get(..sseq.len().saturating_sub(1)).unwrap_or(&[]);
                        let Some(last) = crate::util::char_at(sseq, -1) else {
                            break;
                        };
                        if is_has_and_not_equals_i64(reference, chr, bp_mut64, last) {
                            continue;
                        }
                        eqcnt += 1;
                        if crate::util::char_at(sseq, -2).is_some_and(|prev| {
                            is_has_and_not_equals_i64(reference, chr, bp_mut64 - 1, prev)
                        }) {
                            continue;
                        }
                        extra_slice = seq_vec
                            .get(seq_vec.len().saturating_sub(ii)..)
                            .unwrap_or(&seq_vec);
                    }

                    if eqcnt >= 3 && (eqcnt as f64) / (ii as f64) > 0.5 {
                        break;
                    }
                    if verbose {
                        crate::verbose_eprintln!(
                            "      FoundSEED SV BP: {} BP: {} SEEDpos {} {} {} {} {} EXTRA: {}",
                            dir,
                            bp_mut64,
                            first_seed64,
                            position,
                            String::from_utf8_lossy(seed),
                            i,
                            String::from_utf8_lossy(&seq_vec),
                            String::from_utf8_lossy(extra_slice)
                        );
                    }
                    if let Ok(bp_mut) = i32::try_from(bp_mut64) {
                        if crate::realign::is_match_ref(sseq, reference, chr, bp_mut, dir, 1) {
                            return Match {
                                base_position: bp_mut,
                                matched_sequence: extra_slice.to_vec(),
                            };
                        }
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
        crate::verbose_eprintln!(
            "    Working MatchRev {} {} {} SEED: {}",
            position,
            String::from_utf8_lossy(&seq_vec),
            dir,
            seed_len
        );
    }

    let mut extra_slice: &[u8] = &[];
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
            let first_seed64 = i64::from(*first_seed_u32);
            let bp64 = if dir == 1 {
                first_seed64 + (seq_vec.len() as i64) - (i as i64) - 1
            } else {
                first_seed64 - (i as i64)
            };

            if let Ok(bp) = i32::try_from(bp64) {
                if crate::realign::is_match_ref(&seq_vec, reference, chr, bp, -dir, mm_limit) {
                    if verbose {
                        crate::verbose_eprintln!(
                            "      Found SV BP (reverse): {} BP: {} SEEDpos: {} {} {} {} {}",
                            dir,
                            bp,
                            first_seed64,
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
                }
            }
            // For complex indels, allowing some mismatches at the end up to 15bp or 20% length
            let mut sseq: &[u8] = &seq_vec;
            let mut eqcnt = 0;
            let mut bp_mut64 = bp64;
            for ii in 1..=15 {
                bp_mut64 -= i64::from(dir);
                if dir == -1 {
                    sseq = sseq.get(1..).unwrap_or(&[]);
                    let Some(b0) = crate::util::char_at(sseq, 0) else {
                        break;
                    };
                    if is_has_and_not_equals_i64(reference, chr, bp_mut64, b0) {
                        continue;
                    }
                    eqcnt += 1;
                    if crate::util::char_at(sseq, 1).is_some_and(|b1| {
                        is_has_and_not_equals_i64(reference, chr, bp_mut64 + 1, b1)
                    }) {
                        continue;
                    }
                    extra_slice = seq_vec.get(..ii).unwrap_or(&seq_vec);
                } else {
                    sseq = sseq.get(..sseq.len().saturating_sub(1)).unwrap_or(&[]);
                    let Some(last) = crate::util::char_at(sseq, -1) else {
                        break;
                    };
                    if is_has_and_not_equals_i64(reference, chr, bp_mut64, last) {
                        continue;
                    }
                    eqcnt += 1;
                    if crate::util::char_at(sseq, -2).is_some_and(|prev| {
                        is_has_and_not_equals_i64(reference, chr, bp_mut64 - 1, prev)
                    }) {
                        continue;
                    }
                    extra_slice = seq_vec
                        .get(seq_vec.len().saturating_sub(ii)..)
                        .unwrap_or(&seq_vec);
                }

                if eqcnt >= 3 && (eqcnt as f64) / (ii as f64) > 0.5 {
                    break;
                }
                if verbose {
                    crate::verbose_eprintln!(
                        "      FoundSEED SV BP (reverse): {} BP: {} SEEDpos {} {} {} {} {} EXTRA: {}",
                        dir,
                        bp_mut64,
                        first_seed64,
                        position,
                        String::from_utf8_lossy(seed),
                        i,
                        String::from_utf8_lossy(&seq_vec),
                        String::from_utf8_lossy(extra_slice)
                    );
                }
                if let Ok(bp_mut) = i32::try_from(bp_mut64) {
                    if crate::realign::is_match_ref(sseq, reference, chr, bp_mut, -dir, 1) {
                        return Match {
                            base_position: bp_mut,
                            matched_sequence: extra_slice.to_vec(),
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

    let ins = positions[2].saturating_sub(positions[1]);
    if (end1 != start1) && (end2 != start2) {
        let ratio1 = (ins as f64) / ((i64::from(end1) - i64::from(start1)) as f64);
        let ratio2 = (ins as f64) / ((i64::from(end2) - i64::from(start2)) as f64);
        if ratio1 > 0.75 && ratio2 > 0.75 {
            return true;
        }
    }
    if positions[1]
        .saturating_sub(positions[0])
        .saturating_add(positions[3].saturating_sub(positions[2]))
        < 3i32.saturating_mul(rlen)
    {
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
    let mut cov: i32 = 0;
    let mut pairs: i32 = 0;
    let mut cnt: i32 = 0;

    for current_sclips_sr in structural_variants_sv {
        for sv_r in current_sclips_sr.iter_mut() {
            let (start2, end2) = if sv_r.start < sv_r.mstart {
                (sv_r.end, sv_r.mstart)
            } else {
                (sv_r.mend, sv_r.start)
            };

            if verbose {
                crate::verbose_eprintln!(
                    "   Marking SV {} {} {} {} cnt: {}",
                    start,
                    end,
                    start2,
                    end2,
                    sv_r.variation.cnt
                );
            }

            if is_overlap(start, end, start2, end2, rlen) {
                if verbose {
                    crate::verbose_eprintln!(
                        "       SV {} {} {} {} cnt: {} marked",
                        start,
                        end,
                        start2,
                        end2,
                        sv_r.variation.cnt
                    );
                }
                sv_r.used = true;
                cnt += 1;
                pairs += sv_r.variation.cnt;
                if sv_r.end != sv_r.start {
                    cov = cov.saturating_add(
                        (((i64::from(sv_r.variation.cnt) * i64::from(rlen)) as f64
                            / (i64::from(sv_r.end) - i64::from(sv_r.start)) as f64)
                            as i32)
                            .saturating_add(1),
                    );
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
                crate::verbose_eprintln!(
                    "   Marking DUP SV {} {} {} {} cnt: {}",
                    start,
                    end,
                    start2,
                    end2,
                    sv_r.variation.cnt
                );
            }

            if is_overlap(start, end, start2, end2, rlen) {
                if verbose {
                    crate::verbose_eprintln!(
                        "       DUP SV {} {} {} {} cnt: {} marked",
                        start,
                        end,
                        start2,
                        end2,
                        sv_r.variation.cnt
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
            let mut s = midpoint_i32(svr.start, svr.end);
            let mut e = midpoint_i32(svr.mstart, svr.mend);
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
                crate::verbose_eprintln!("      Pair [{}:{}] found and marked.", s, e);
            }
        }
    }
    pairs_data
}

struct SortPositionCount {
    position: i32,
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

        let mut softp = {
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
            best_soft_pos(&inv.soft).unwrap_or(0)
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
                end: Pos1(mend.max(1) as u32),
                gene: String::new(),
                amplicon: None,
            };
            get_reference(&mod_region, 500, chromosomes_lengths, config, reference)?;
            partial_pipeline(
                bam_path,
                var_data,
                reference,
                chr,
                mstart.saturating_sub(200).max(1),
                mend.saturating_add(200),
                config,
                chromosomes_lengths,
            )?;
        }

        let mut bp = 0;
        let mut scv_variation_opt = None;
        let mut scv_used_pos: i32 = softp;
        let mut seq = Vec::new();
        let mut extra: Vec<u8> = Vec::new();

        if softp != 0 {
            let sclip_map = if dir == 1 {
                &mut var_data.soft_clips_3_end
            } else {
                &mut var_data.soft_clips_5_end
            };
            let Some(scv_mut) = sclip_map.get_mut(&softp) else {
                continue;
            };
            if scv_mut.used {
                continue;
            }
            seq = crate::variations::util::find_conseq(
                scv_mut,
                0,
                &adaptor_forward,
                &adaptor_reverse,
            );
            if seq.is_empty() {
                continue;
            }
            let scv_var = scv_mut.variation.clone();

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
            scv_variation_opt = Some(scv_var);
        } else {
            // Look within 2*maxReadLength to find soft clips not associated with discordant pairs
            let sp = if dir == 1 { end } else { start };
            for i in 1..=2i32.saturating_mul(max_read_length) {
                let cp = sp.saturating_add(i.saturating_mul(dir));
                let sclip_map = if dir == 1 {
                    &mut var_data.soft_clips_3_end
                } else {
                    &mut var_data.soft_clips_5_end
                };
                let Some(scv_mut) = sclip_map.get_mut(&cp) else {
                    continue;
                };
                if scv_mut.used {
                    continue;
                }
                let seq_check = crate::variations::util::find_conseq(
                    scv_mut,
                    0,
                    &adaptor_forward,
                    &adaptor_reverse,
                );
                if seq_check.is_empty() {
                    continue;
                }
                let scv_var = scv_mut.variation.clone();

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
                scv_variation_opt = Some(scv_var);
                scv_used_pos = cp;
                seq = seq_check;
                let bp_diff = if dir == 1 {
                    bp.saturating_sub(mend).saturating_abs()
                } else {
                    bp.saturating_sub(mstart).saturating_abs()
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

        let scv_variation = match scv_variation_opt {
            Some(v) => v,
            None => continue,
        };

        if side == Side::Five {
            if dir == -1 {
                bp = bp.saturating_sub(1);
            }
        } else {
            if dir == 1 {
                bp = bp.saturating_add(1);
                if bp != 0 {
                    softp = softp.saturating_sub(1);
                }
            } else {
                softp = softp.saturating_sub(1);
            }
        }

        if side == Side::Three {
            std::mem::swap(&mut bp, &mut softp);
        }

        // Refine bp and softp by shifting while complement matches
        if (dir == -1 && side == Side::Five) || (dir == 1 && side == Side::Three) {
            while let (Some(ch_s), Some(ch_b)) = (
                u32::try_from(softp)
                    .ok()
                    .and_then(|pos| reference.get_base(chr, pos)),
                u32::try_from(bp)
                    .ok()
                    .and_then(|pos| reference.get_base(chr, pos)),
            ) {
                if ch_s == crate::util::complement(ch_b) {
                    let Some(next_softp) = softp.checked_add(1) else {
                        break;
                    };
                    softp = next_softp;
                    if softp != 0 {
                        let Some(next_bp) = bp.checked_sub(1) else {
                            break;
                        };
                        bp = next_bp;
                    }
                } else {
                    break;
                }
            }
        }
        while let (Some(ch_s), Some(ch_b)) = (
            softp
                .checked_sub(1)
                .and_then(|p| u32::try_from(p).ok())
                .and_then(|pos| reference.get_base(chr, pos)),
            bp.checked_add(1)
                .and_then(|p| u32::try_from(p).ok())
                .and_then(|pos| reference.get_base(chr, pos)),
        ) {
            if ch_s == crate::util::complement(ch_b) {
                let Some(next_softp) = softp.checked_sub(1) else {
                    break;
                };
                softp = next_softp;
                if softp != 0 {
                    let Some(next_bp) = bp.checked_add(1) else {
                        break;
                    };
                    bp = next_bp;
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
        let gap = bp.saturating_sub(softp);
        let abs_mlen = (inv_mlen as f64).abs();
        if bp > softp && gap > 150 && (gap as f64) / abs_mlen < 1.5 {
            let len = gap.saturating_add(1);
            let sv_flank = 50; // Configuration.SVFLANK
            let ins5 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                bp.saturating_sub(sv_flank - 1),
                bp,
            ));
            let ins3 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                softp,
                softp.saturating_add(sv_flank - 1),
            ));
            let ins = if len.saturating_sub(2 * sv_flank) <= 0 {
                complement_seq(&crate::util::join_ref(reference, chr, softp, bp))
            } else {
                let mut tmp_ins = ins5;
                tmp_ins.push_str(&format!("<inv{}>", len.saturating_sub(2 * sv_flank)));
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
                crate::variations::util::adj_cnt(variation, &scv_variation, None);
            }

            {
                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    softp,
                );
                sv.r#type = "INV".to_string();
                sv.splits += scv_variation.cnt;
                sv.pairs += cnt;
                sv.clusters += 1;
            }

            if dir == -1 {
                if let Some(ref_base) = reference.get_base(chr, softp as u32) {
                    let ref_base_str = (ref_base as char).to_string();
                    if let Some(entry) = var_data.non_insertion_variants.get_mut(&softp) {
                        if let Some(ref_var) = entry.map.get_mut(&ref_base_str) {
                            ref_var.cnt -= scv_variation.cnt;
                            ref_var.hicnt -= scv_variation.hicnt;
                            ref_var.locnt -= scv_variation.locnt;
                            ref_var.pmean -= scv_variation.pmean;
                            ref_var.qmean -= scv_variation.qmean;
                            ref_var.mqmean -= scv_variation.mqmean;
                            ref_var.nm -= scv_variation.nm;
                            ref_var.sub_dir(true, scv_variation.get_dir(true));
                            ref_var.sub_dir(false, scv_variation.get_dir(false));
                            crate::variations::util::correct_cnt(ref_var);
                        }
                    }
                }
            }

            let mut dels5: BTreeMap<i32, rustc_hash::FxHashMap<String, i32>> = BTreeMap::new();
            let mut map = rustc_hash::FxHashMap::default();
            map.insert(gt_full.clone(), cnt);
            dels5.insert(softp, map);

            let ref_cov_val =
                if let Some(&prev_cov) = var_data.ref_coverage.get(&softp.saturating_sub(1)) {
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

            let bams_vec: Vec<String> = bam_path
                .to_string_lossy()
                .split(':')
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .collect();
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
            let softp = best_soft_pos(&del.soft).unwrap_or(0);
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
            sync_sclip_after_conseq(&mut var_data.soft_clips_3_end, softp_cand, &scv);
            if seq.is_empty() || seq.len() < 12 {
                continue;
            }

            if !is_loaded(chr, mstart as u32, mend as u32, reference) {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(mstart.max(1) as u32),
                    end: Pos1(mend.max(1) as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
                partial_pipeline(
                    bam_path,
                    var_data,
                    reference,
                    chr,
                    mstart.saturating_sub(200).max(1),
                    mend.saturating_add(200),
                    config,
                    chromosomes_lengths,
                )?;
            }

            let match_res = find_match(&seq, reference, softp_cand, 1, 17, 3, chr, verbose);
            let mut bp = match_res.base_position;
            if bp == 0 {
                continue;
            }
            if !(bp.saturating_sub(softp_cand) > 30
                && is_overlap(softp_cand, bp, end, mstart, var_data.max_read_length))
            {
                continue;
            }
            bp = bp.saturating_sub(1);
            let mut softp = softp_cand;
            let mut dellen = bp.saturating_sub(softp).saturating_add(1);
            while let (Some(b_bp), Some(b_softp_prev)) = (
                u32::try_from(bp)
                    .ok()
                    .and_then(|pos| reference.get_base(chr, pos)),
                softp
                    .checked_sub(1)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos)),
            ) {
                if b_bp == b_softp_prev {
                    let Some(next_bp) = bp.checked_sub(1) else {
                        break;
                    };
                    bp = next_bp;
                    if bp != 0 {
                        softp = softp.saturating_sub(1);
                    }
                } else {
                    break;
                }
            }
            dellen = bp.saturating_sub(softp).saturating_add(1);

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
                    end: Pos1(mend.max(1) as u32),
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
                if !(i >= end.saturating_sub(3)
                    && i.saturating_sub(end) < 3i32.saturating_mul(max_read_length))
                {
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
                if !(bp.saturating_sub(i) > 30
                    && is_overlap(i, bp, end, mstart, var_data.max_read_length))
                {
                    continue;
                }
                bp = bp.saturating_sub(1);
                matched_j = Some((i, bp, scv.clone()));
                break;
            }

            if let Some((softp, bp, scv)) = matched_j {
                let dellen = bp.saturating_sub(softp).saturating_add(1);
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
            let softp = best_soft_pos(&del.soft).unwrap_or(0);
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
            sync_sclip_after_conseq(&mut var_data.soft_clips_5_end, softp_cand, &scv);
            if seq.is_empty() || seq.len() < 12 {
                continue;
            }

            let is_loaded_before = is_loaded(chr, mstart as u32, mend as u32, reference);
            if !is_loaded_before {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(mstart.max(1) as u32),
                    end: Pos1(mend.max(1) as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(&mod_region, 300, chromosomes_lengths, config, reference)?;
                partial_pipeline(
                    bam_path,
                    var_data,
                    reference,
                    chr,
                    mstart.saturating_sub(200).max(1),
                    mend.saturating_add(200),
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
            let overlap_ok = softp_cand.saturating_sub(bp) > 30
                && is_overlap(bp, softp_cand, mend, start, var_data.max_read_length);
            if !overlap_ok {
                continue;
            }
            bp = bp.saturating_add(1);
            let mut softp = softp_cand.saturating_sub(1);
            let mut dellen = softp.saturating_sub(bp).saturating_add(1);

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
                    end: Pos1(mend.max(1) as u32),
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
                if !(i <= start.saturating_add(3)
                    && start.saturating_sub(i) < 3i32.saturating_mul(max_read_length))
                {
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
                if !(i.saturating_sub(bp) > 30
                    && is_overlap(bp, i, mend, start, var_data.max_read_length))
                {
                    continue;
                }
                bp = bp.saturating_add(1);
                matched_j = Some((i, bp, scv.clone()));
                break;
            }

            if let Some((softp_cand, bp, scv)) = matched_j {
                let mut softp = softp_cand.saturating_sub(1);
                let dellen = softp.saturating_sub(bp).saturating_add(1);

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
    let mut tmp5: Vec<SortPositionCount> = var_data
        .soft_clips_5_end
        .iter()
        .filter_map(|(&pos, scv)| {
            if scv.used || scv.variation.cnt < minr || pos < curseg_start || pos > curseg_end {
                None
            } else {
                Some(SortPositionCount {
                    position: pos,
                    count: scv.variation.cnt,
                })
            }
        })
        .collect();
    tmp5.sort_by_key(|item| std::cmp::Reverse(item.count));

    let mut tmp3: Vec<SortPositionCount> = var_data
        .soft_clips_3_end
        .iter()
        .filter_map(|(&pos, scv)| {
            if scv.used || scv.variation.cnt < minr || pos < curseg_start || pos > curseg_end {
                None
            } else {
                Some(SortPositionCount {
                    position: pos,
                    count: scv.variation.cnt,
                })
            }
        })
        .collect();
    tmp3.sort_by_key(|item| std::cmp::Reverse(item.count));

    for tuple5 in tmp5 {
        let p5 = tuple5.position;
        let cnt5 = tuple5.count;
        if cnt5 < minr {
            break;
        }

        // Check if soft clip actually got used in the meantime
        if let Some(sc_mut) = var_data.soft_clips_5_end.get(&p5) {
            if sc_mut.used {
                continue;
            }
        } else {
            continue;
        }

        // SOFTP2SV check: skip if the SV at this position was already used
        if is_best_sv_used(var_data, p5) {
            continue;
        }

        let (seq, sc5v_variation) = {
            let Some(sc5v) = var_data.soft_clips_5_end.get_mut(&p5) else {
                continue;
            };
            let seq =
                crate::variations::util::find_conseq(sc5v, 0, &adaptor_forward, &adaptor_reverse);
            (seq, sc5v.variation.clone())
        };
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

                let p5_adj = p5.saturating_sub(1);
                let bp_adj = bp.saturating_add(1);
                let dellen = p5_adj.saturating_sub(bp_adj).saturating_add(1);

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
                if bp_cov < pairs + sc5v_variation.cnt {
                    var_data
                        .ref_coverage
                        .insert(bp_adj, pairs + sc5v_variation.cnt);
                }
                // Ports StructuralVariantsProcessor.java:732,748: Java mutates `p5--`
                // at line 732 before `refCoverage.get(p5 + 1)` at line 748, so `p5 + 1`
                // in Java is `p5_adj + 1` (the original `p5` before decrement).
                let p5_cov = var_data
                    .ref_coverage
                    .get(&p5_adj.saturating_add(1))
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
                crate::variations::util::adj_cnt(vref_mut, &sc5v_variation, None);

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
            }
        } else {
            let match_res_rev = find_match_rev(&seq, reference, p5, -1, 17, 3, chr, verbose);
            let mut bp = match_res_rev.base_position;
            if bp == 0 {
                continue;
            }
            if bp.saturating_sub(p5).saturating_abs() <= 50 {
                continue;
            }
            let mut p5_adj = p5;
            if bp < p5_adj {
                // Keep bp at 3' side and p5_adj at 5' side
                std::mem::swap(&mut bp, &mut p5_adj);
            }
            bp = bp.saturating_sub(1);
            while let (Some(b1), Some(b2)) = (
                bp.checked_add(1)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos)),
                p5_adj
                    .checked_sub(1)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos)),
            ) {
                if crate::util::complement(b1) != b2 {
                    break;
                }
                let Some(next_p5) = p5_adj.checked_sub(1) else {
                    break;
                };
                p5_adj = next_p5;
                if p5_adj != 0 {
                    let Some(next_bp) = bp.checked_add(1) else {
                        break;
                    };
                    bp = next_bp;
                }
            }

            let ins5 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                bp.saturating_sub(49),
                bp,
            ));
            let ins3 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                p5_adj,
                p5_adj.saturating_add(49),
            ));
            let span = bp.saturating_sub(p5_adj).saturating_add(1);
            let mid = span
                .saturating_sub(ins5.len() as i32)
                .saturating_sub(ins3.len() as i32);

            let extra = String::from_utf8_lossy(&match_res_rev.matched_sequence);
            let mut vn = format!("-{}^{}<inv{}>{}{}", span, ins5, mid, ins3, extra);
            if mid <= 0 {
                let tins = complement_seq(&crate::util::join_ref(reference, chr, p5_adj, bp));
                vn = format!("-{}^{}{}", span, tins, extra);
            }

            let _vref =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p5_adj, &vn);

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
            crate::variations::util::adj_cnt(vref_mut, &sc5v_variation, None);

            let p5_cov = var_data.ref_coverage.get(&p5_adj).copied().unwrap_or(0);
            var_data.ref_coverage.insert(p5_adj, p5_cov + cnt5);
            let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
            let p5_cov_updated = var_data.ref_coverage.get(&p5_adj).copied().unwrap_or(0);
            if p5_cov_updated < bp_cov {
                var_data.ref_coverage.insert(p5_adj, bp_cov);
            }
        }
    }

    for tuple3 in tmp3 {
        let p3 = tuple3.position;
        let cnt3 = tuple3.count;
        if cnt3 < minr {
            break;
        }

        if let Some(sc_mut) = var_data.soft_clips_3_end.get(&p3) {
            if sc_mut.used {
                continue;
            }
        } else {
            continue;
        }

        // SOFTP2SV check: skip if the SV at this position was already used
        if is_best_sv_used(var_data, p3) {
            continue;
        }

        let (seq, sc3v_variation) = {
            let Some(sc3v) = var_data.soft_clips_3_end.get_mut(&p3) else {
                continue;
            };
            let seq =
                crate::variations::util::find_conseq(sc3v, 0, &adaptor_forward, &adaptor_reverse);
            (seq, sc3v.variation.clone())
        };

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

                let dellen = bp.saturating_sub(p3);
                let mut bp_adj = bp.saturating_sub(1);
                let mut p3_adj = p3;
                while let (Some(b_bp), Some(b_p3_prev)) = (
                    u32::try_from(bp_adj)
                        .ok()
                        .and_then(|pos| reference.get_base(chr, pos)),
                    p3_adj
                        .checked_sub(1)
                        .and_then(|p| u32::try_from(p).ok())
                        .and_then(|pos| reference.get_base(chr, pos)),
                ) {
                    if b_bp == b_p3_prev {
                        let Some(next_bp) = bp_adj.checked_sub(1) else {
                            break;
                        };
                        bp_adj = next_bp;
                        if bp_adj != 0 {
                            p3_adj = p3_adj.saturating_sub(1);
                        }
                    } else {
                        break;
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
                        .insert(p3_adj, pairs + sc3v_variation.cnt);
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
                crate::variations::util::adj_cnt(vref_mut, &sc3v_variation, None);

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
            }
        } else {
            let match_res_rev = find_match_rev(&seq, reference, p3, 1, 17, 3, chr, verbose);
            let mut bp = match_res_rev.base_position;
            if bp == 0 {
                continue;
            }
            if bp.saturating_sub(p3).saturating_abs() <= 50 {
                continue;
            }
            let mut p3_adj = p3;
            if bp < p3_adj {
                std::mem::swap(&mut bp, &mut p3_adj);
                p3_adj = p3_adj.saturating_add(1);
                bp = bp.saturating_sub(1);
            }
            while let (Some(b1), Some(b2)) = (
                bp.checked_add(1)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos)),
                p3_adj
                    .checked_sub(1)
                    .and_then(|p| u32::try_from(p).ok())
                    .and_then(|pos| reference.get_base(chr, pos)),
            ) {
                if crate::util::complement(b1) != b2 {
                    break;
                }
                let Some(next_p3) = p3_adj.checked_sub(1) else {
                    break;
                };
                p3_adj = next_p3;
                if p3_adj != 0 {
                    let Some(next_bp) = bp.checked_add(1) else {
                        break;
                    };
                    bp = next_bp;
                }
            }

            let ins5 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                bp.saturating_sub(49),
                bp,
            ));
            let ins3 = complement_seq(&crate::util::join_ref(
                reference,
                chr,
                p3_adj,
                p3_adj.saturating_add(49),
            ));
            let span = bp.saturating_sub(p3_adj).saturating_add(1);
            let mid = span.saturating_sub(100);

            let extra = String::from_utf8_lossy(&match_res_rev.matched_sequence);
            let mut vn = format!("-{}^{}{}<inv{}>{}", span, extra, ins5, mid, ins3);
            if mid <= 0 {
                let tins = complement_seq(&crate::util::join_ref(reference, chr, p3_adj, bp));
                vn = format!("-{}^{}{}", span, extra, tins);
            }

            let _vref =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p3_adj, &vn);

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
            crate::variations::util::adj_cnt(vref_mut, &sc3v_variation, None);

            let p3_cov = var_data.ref_coverage.get(&p3_adj).copied().unwrap_or(0);
            var_data.ref_coverage.insert(p3_adj, p3_cov + cnt3);
            let bp_cov = var_data.ref_coverage.get(&bp).copied().unwrap_or(0);
            let p3_cov_updated = var_data.ref_coverage.get(&p3_adj).copied().unwrap_or(0);
            if p3_cov_updated < bp_cov {
                var_data.ref_coverage.insert(p3_adj, bp_cov);
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
    let mindist = 8i32.saturating_mul(var_data.max_read_length);

    // svfdel (5' candidate deletions)
    let n_svfdel = var_data.sv_structures.svfdel.len();
    for idx in 0..n_svfdel {
        let (used, vars_count, start, end, mstart, mend, softp) = {
            let del = &var_data.sv_structures.svfdel[idx];
            if del.used || del.variation.cnt < minr.saturating_add(5) {
                continue;
            }
            if !var_data.splice.is_empty() && i64::from(del.mlen).abs() < 250_000 {
                continue;
            }
            if del.mstart <= del.end.saturating_add(mindist) {
                continue;
            }
            if del.variation.mqmean / (del.variation.cnt as f64) <= 35.0 {
                // DISCPAIRQUAL = 35
                continue;
            }
            let mlen = del
                .mstart
                .saturating_sub(del.end)
                .saturating_sub(var_data.max_read_length / del.variation.cnt.saturating_add(1));
            if !(mlen > 0 && mlen > mindist) {
                continue;
            }
            let bp = if del.softp != 0 {
                del.softp
            } else {
                del.end.saturating_add(
                    (var_data.max_read_length / del.variation.cnt.saturating_add(1)) / 2,
                )
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

        let mlen = mstart
            .saturating_sub(end)
            .saturating_sub(var_data.max_read_length / vars_count.saturating_add(1));
        if reference.get_base(chr, softp as u32).is_none() {
            let mod_region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(softp.saturating_sub(150).max(1) as u32),
                end: Pos1(softp.saturating_add(150).max(1) as u32),
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
            .get(&end.saturating_add(1))
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
        let (used, vars_count, start, end, mstart, mend, softp, del_softp) = {
            let del = &var_data.sv_structures.svrdel[idx];
            if del.used || del.variation.cnt < minr.saturating_add(5) {
                continue;
            }
            if !var_data.splice.is_empty() && i64::from(del.mlen).abs() < 250_000 {
                continue;
            }
            if del.start <= del.mend.saturating_add(mindist) {
                continue;
            }
            if del.variation.mqmean / (del.variation.cnt as f64) <= 35.0 {
                continue;
            }
            let mlen = del
                .start
                .saturating_sub(del.mend)
                .saturating_sub(var_data.max_read_length / del.variation.cnt.saturating_add(1));
            if !(mlen > 0 && mlen > mindist) {
                continue;
            }
            let bp = del.mend.saturating_add(
                (var_data.max_read_length / del.variation.cnt.saturating_add(1)) / 2,
            );
            (
                del.used,
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                bp,
                del.softp,
            )
        };

        if used {
            continue;
        }

        let mlen = start
            .saturating_sub(mend)
            .saturating_sub(var_data.max_read_length / vars_count.saturating_add(1));
        if reference.get_base(chr, softp as u32).is_none() {
            let mod_region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(softp.saturating_sub(150).max(1) as u32),
                end: Pos1(softp.saturating_add(150).max(1) as u32),
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
            .get(&mend.saturating_add(1))
            .map(|s| s.variation.cnt)
            .unwrap_or(0);
        sv.splits += var_data
            .soft_clips_5_end
            .get(&start)
            .map(|s| s.variation.cnt)
            .unwrap_or(0);
        sv.pairs += vars_count;
        sv.clusters += 1;

        if del_softp != 0 {
            if let Some(sc5) = var_data.soft_clips_5_end.get_mut(&del_softp) {
                sc5.used = true;
            }
        }

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
            start: Pos1(mstart.saturating_sub(100).max(1) as u32),
            end: Pos1(mend.saturating_add(100).max(1) as u32),
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
                let bp = midpoint_i32(end, rstart).saturating_abs();
                let pe = midpoint_i32(me, rms).saturating_abs();

                if reference.get_base(chr, pe as u32).is_none() {
                    let modified_region = crate::regions::Region {
                        chr: chr.to_string(),
                        start: Pos1(pe.saturating_sub(150).max(1) as u32),
                        end: Pos1(pe.saturating_add(150).max(1) as u32),
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

                let len = pe.saturating_sub(bp).saturating_add(1);
                let ins5 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    bp,
                    bp.saturating_add(sv_flank - 1),
                ));
                let ins3 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    pe.saturating_sub(sv_flank - 1),
                    pe,
                ));

                let mut ins = format!("{}<inv{}>{}", ins3, len.saturating_sub(2 * sv_flank), ins5);
                if len.saturating_sub(2 * sv_flank) <= 0 {
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
                let pe = midpoint_i32(end, rstart).saturating_abs();
                let bp = midpoint_i32(me, rms).saturating_abs();

                if reference.get_base(chr, bp as u32).is_none() {
                    let modified_region = crate::regions::Region {
                        chr: chr.to_string(),
                        start: Pos1(bp.saturating_sub(150).max(1) as u32),
                        end: Pos1(bp.saturating_add(150).max(1) as u32),
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

                let len = pe.saturating_sub(bp).saturating_add(1);
                let ins5 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    bp,
                    bp.saturating_add(sv_flank - 1),
                ));
                let ins3 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    pe.saturating_sub(sv_flank - 1),
                    pe,
                ));

                let mut ins = format!("{}<inv{}>{}", ins3, len.saturating_sub(2 * sv_flank), ins5);
                if len.saturating_sub(2 * sv_flank) <= 0 {
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
                    .get(&end.saturating_add(1))
                    .map_or(0, |s| s.variation.cnt);
                sv.splits += var_data
                    .soft_clips_3_end
                    .get(&me.saturating_add(1))
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

        let mut mlen = end
            .saturating_sub(mstart)
            .saturating_add(var_data.max_read_length / vars_count);
        let mut bp = mstart.saturating_sub((var_data.max_read_length / vars_count) / 2);
        let mut pe = end;

        if !is_loaded(chr, mstart as u32, mend as u32, reference) {
            if reference.get_base(chr, bp as u32).is_none() {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(bp.saturating_sub(150).max(1) as u32),
                    end: Pos1(bp.saturating_add(150).max(1) as u32),
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
                mstart.saturating_sub(200).max(1),
                mend.saturating_add(200),
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
            if let Some(best_pe) = best_soft_pos(&var_data.sv_structures.svfdup[idx].soft) {
                pe = best_pe;
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
            sync_sclip_after_conseq(&mut var_data.soft_clips_3_end, pe, &current_sclip3);
            let match_res = find_match(&seq, reference, bp, 1, 17, 3, chr, verbose);
            let mut tbp = match_res.base_position;

            if tbp != 0 && tbp < pe {
                if let Some(scv_mut) = var_data.soft_clips_3_end.get_mut(&pe) {
                    scv_mut.used = true;
                }
                while let (Some(b_pe_prev), Some(b_tbp_prev)) = (
                    pe.checked_sub(1)
                        .and_then(|p| u32::try_from(p).ok())
                        .and_then(|pos| reference.get_base(chr, pos)),
                    tbp.checked_sub(1)
                        .and_then(|p| u32::try_from(p).ok())
                        .and_then(|pos| reference.get_base(chr, pos)),
                ) {
                    if b_pe_prev == b_tbp_prev {
                        let Some(next_tbp) = tbp.checked_sub(1) else {
                            break;
                        };
                        tbp = next_tbp;
                        if tbp != 0 {
                            pe = pe.saturating_sub(1);
                        }
                    } else {
                        break;
                    }
                }
                mlen = pe.saturating_sub(tbp);
                bp = tbp;
                pe = pe.saturating_sub(1);
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

        let ins5 = crate::util::join_ref(reference, chr, bp, bp.saturating_add(49));
        let ins3 = crate::util::join_ref(reference, chr, pe.saturating_sub(49), pe);
        let ins = format!("{}<dup{}>{}", ins5, mlen.saturating_sub(100), ins3);

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

        let mut mlen = mend
            .saturating_sub(start)
            .saturating_add(var_data.max_read_length / vars_count);
        let mut bp = start.saturating_sub((var_data.max_read_length / vars_count) / 2);
        let mut pe = (i64::from(mlen) + i64::from(bp) - 1)
            .clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
        let mut tpe64 = i64::from(pe);

        if !is_loaded(chr, mstart as u32, mend as u32, reference) {
            if reference.get_base(chr, pe as u32).is_none() {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(pe.saturating_sub(150).max(1) as u32),
                    end: Pos1(pe.saturating_add(150).max(1) as u32),
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
                mstart.saturating_sub(200).max(1),
                mend.saturating_add(200),
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
            if let Some(best_bp) = best_soft_pos(&var_data.sv_structures.svrdup[idx].soft) {
                bp = best_bp;
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
            sync_sclip_after_conseq(&mut var_data.soft_clips_5_end, bp, &current_sclip5);
            let match_res = find_match(&seq, reference, pe, -1, 17, 3, chr, verbose);
            let tbp = match_res.base_position;

            if tbp != 0 && tbp > bp {
                if let Some(scv_mut) = var_data.soft_clips_5_end.get_mut(&bp) {
                    scv_mut.used = true;
                }
                pe = tbp;
                mlen = pe.saturating_sub(bp).saturating_add(1);
                tpe64 = i64::from(pe) + 1;
                while let (Some(b_tpe), Some(b_bp)) = (
                    i32::try_from(tpe64)
                        .ok()
                        .and_then(|p| u32::try_from(p).ok())
                        .and_then(|pos| reference.get_base(chr, pos)),
                    i32::try_from(i64::from(bp) + (tpe64 - i64::from(pe) - 1))
                        .ok()
                        .and_then(|p| u32::try_from(p).ok())
                        .and_then(|pos| reference.get_base(chr, pos)),
                ) {
                    if b_tpe == b_bp {
                        tpe64 += 1;
                    } else {
                        break;
                    }
                }
                if let Ok(tpe_i32) = i32::try_from(tpe64) {
                    if let Some(current_sclip3) = var_data.soft_clips_3_end.get(&tpe_i32) {
                        cntf = current_sclip3.variation.cnt;
                        qmeanf = current_sclip3.variation.qmean;
                        Qmeanf = current_sclip3.variation.mqmean;
                        pmeanf = current_sclip3.variation.pmean;
                        nmf = current_sclip3.variation.nm;
                    }
                }
            }
        }

        let ins5 = crate::util::join_ref(reference, chr, bp, bp.saturating_add(49));
        let ins3 = crate::util::join_ref(reference, chr, pe.saturating_sub(49), pe);
        let ins = format!("{}<dup{}>{}", ins5, mlen.saturating_sub(100), ins3);

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
            if let Ok(tpe_i32) = i32::try_from(tpe64) {
                sv.splits += var_data
                    .soft_clips_3_end
                    .get(&tpe_i32)
                    .map(|s| s.variation.cnt)
                    .unwrap_or(0);
            }
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
        crate::verbose_eprintln!("Start Structural Variants: DEL\n");
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
        crate::verbose_eprintln!("Start Structural Variants: INV\n");
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
        crate::verbose_eprintln!("Start Structural Variants\n");
    }
    find_sv(
        var_data,
        reference,
        chr,
        i32::try_from(pileup_region.start.0).unwrap_or(i32::MAX),
        i32::try_from(pileup_region.end.0).unwrap_or(i32::MAX),
        config,
        chromosomes_lengths,
    )?;

    if verbose_log(config) {
        crate::verbose_eprintln!("Start Structural Variants: DEL discordant pairs only\n");
    }
    find_del_disc(var_data, reference, chr, config, chromosomes_lengths)?;

    if verbose_log(config) {
        crate::verbose_eprintln!("Start Structural Variants: INV discordant pairs only\n");
    }
    find_inv_disc(var_data, reference, chr, config, chromosomes_lengths)?;

    if verbose_log(config) {
        crate::verbose_eprintln!("Start Structural Variants: DUP discordant pairs only\n");
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
    crate::count!(SvPartialPipeline);
    // Query BAM records in the modified region
    let query_start = crate::types::Pos0((start.max(1) - 1) as u32);
    let query_end = crate::types::Pos0(end.max(1) as u32);

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
    // G3 (mate-window ownership): record new keys planted by this mate-window
    // re-pileup via JavaHashMap's O(|Delta|) new_key_log so they are attributed in
    // `merge_planted` for the ToVarsBuilder ownership gate without cloning full-tile
    // key sets. Pre-existing keys are native (the merge only accumulates into them),
    // never planted.
    var_data.non_insertion_variants.enable_new_key_log();
    var_data.insertion_variants.enable_new_key_log();
    let orig_sv_structures = std::mem::take(&mut var_data.sv_structures);

    // Streaming counterpart of the old query_region-Vec + pileup_into shape:
    // one lazy pass through a single preprocessor into a parser seeded with
    // the taken maps, written back unconditionally (fresh SVStructures
    // discarded, discordant seeded 0 — duprate never propagates out of this
    // path either way) so `var_data` is never wiped on `Err` and `new_key_log`
    // is never left armed for the remainder of the tile. Peak memory stays flat
    // in depth.
    let run_res: Result<(), VardictError> = {
        let mut parser = crate::cigar::parse::CigarParser::new(
            reference,
            &region,
            config,
            true, // svflag is true for partial pipeline
            0,
            0,
            var_data.max_read_length,
        );
        parser.non_insertion_variants = std::mem::take(&mut var_data.non_insertion_variants);
        parser.insertion_variants = std::mem::take(&mut var_data.insertion_variants);
        parser.position_to_insertion_count =
            std::mem::take(&mut var_data.position_to_insertion_count);
        parser.position_to_deletion_count =
            std::mem::take(&mut var_data.position_to_deletion_count);
        parser.sv_structures = std::mem::take(&mut var_data.sv_structures);
        parser.ref_coverage = std::mem::take(&mut var_data.ref_coverage);
        parser.soft_clips_5_end = std::mem::take(&mut var_data.soft_clips_5_end);
        parser.soft_clips_3_end = std::mem::take(&mut var_data.soft_clips_3_end);
        parser.splice = std::mem::take(&mut var_data.splice);
        parser.mnp = std::mem::take(&mut var_data.mnp);
        parser.splice_count = std::mem::take(&mut var_data.splice_count);

        let run_res: Result<(), VardictError> = (|| {
            let bam_str = bam_path.to_string_lossy();
            let mut feeds = Vec::with_capacity(2);
            for single_bam in bam_str.split(':').filter(|s| !s.is_empty()) {
                feeds.push(crate::bam::query_region_iter(
                    std::path::Path::new(single_bam),
                    chr,
                    query_start,
                    query_end,
                    config,
                )?);
            }

            let mut preprocessor = crate::bam::RecordPreprocessor::new();
            let mut r = crate::bam::SamRecord::new();
            for (idx, mut feed) in feeds.into_iter().enumerate() {
                preprocessor.reset_for_next_bam();
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
            Ok(())
        })();

        var_data.non_insertion_variants = parser.non_insertion_variants;
        var_data.insertion_variants = parser.insertion_variants;
        var_data.position_to_insertion_count = parser.position_to_insertion_count;
        var_data.position_to_deletion_count = parser.position_to_deletion_count;
        var_data.ref_coverage = parser.ref_coverage;
        var_data.soft_clips_5_end = parser.soft_clips_5_end;
        var_data.soft_clips_3_end = parser.soft_clips_3_end;
        var_data.splice = parser.splice;
        var_data.mnp = parser.mnp;
        var_data.splice_count = parser.splice_count;
        if run_res.is_ok() {
            var_data.max_read_length = var_data.max_read_length.max(parser.max_read_length);
        }
        run_res
    };
    var_data.sv_structures = orig_sv_structures;
    let delta_non = var_data.non_insertion_variants.take_new_key_log();
    let mut delta_ins = var_data.insertion_variants.take_new_key_log();
    crate::count!(SvSnapshotKeys, delta_non.len() + delta_ins.len());

    run_res?;

    // Attribute mate-window-planted keys: positions present now but
    // absent from BOTH maps before the merge (`(N_1 ∪ I_1) \ (N_0 ∪ I_0)`).
    // Pre-existing keys are native (the merge only accumulates into them),
    // never planted — even when the merge adds the FIRST entry of the other
    // map at that position, since the ToVarsBuilder ownership gate drops the
    // whole position.
    delta_ins.sort_unstable();
    for pos in delta_non {
        if !var_data.insertion_variants.contains_key(&pos) || delta_ins.binary_search(&pos).is_ok()
        {
            var_data.merge_planted.insert(pos);
        }
    }
    for &pos in &delta_ins {
        if !var_data.non_insertion_variants.contains_key(&pos) {
            var_data.merge_planted.insert(pos);
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
    fn is_best_sv_used_follows_java_string_hash_order() {
        // chrX hashes to bucket 5 at cap=16 (hashCode=3052875), while chr1 hashes to bucket 10 (hashCode=3052836).
        // Java HashMap iterates [chrX, chr1] regardless of insertion order, whereas lexicographic sort iterates [chr1, chrX].
        use crate::variations::{Sclip, VariationData};

        // Case a: chr1 unused (inserted first), chrX used -> Java scans chrX first -> true
        let mut var_data_a = VariationData::default();
        let mut sc_chr1 = Sclip::new();
        sc_chr1.softp = 100;
        sc_chr1.variation.cnt = 4;
        sc_chr1.used = false;
        let mut sc_chrx = Sclip::new();
        sc_chrx.softp = 100;
        sc_chrx.variation.cnt = 4;
        sc_chrx.used = true;
        var_data_a
            .sv_structures
            .svffus
            .insert("chr1".to_string(), vec![sc_chr1]);
        var_data_a
            .sv_structures
            .svffus
            .insert("chrX".to_string(), vec![sc_chrx]);
        assert!(is_best_sv_used(&var_data_a, 100));

        // Case b: chr1 used (inserted first), chrX unused -> Java scans chrX first -> false
        let mut var_data_b = VariationData::default();
        let mut sc_chr1_b = Sclip::new();
        sc_chr1_b.softp = 100;
        sc_chr1_b.variation.cnt = 4;
        sc_chr1_b.used = true;
        let mut sc_chrx_b = Sclip::new();
        sc_chrx_b.softp = 100;
        sc_chrx_b.variation.cnt = 4;
        sc_chrx_b.used = false;
        var_data_b
            .sv_structures
            .svffus
            .insert("chr1".to_string(), vec![sc_chr1_b]);
        var_data_b
            .sv_structures
            .svffus
            .insert("chrX".to_string(), vec![sc_chrx_b]);
        assert!(!is_best_sv_used(&var_data_b, 100));
    }

    #[test]
    fn is_best_sv_used_same_bucket_keeps_insertion_order() {
        // chr1 (hashCode=3052836) and chr12 (hashCode=94637966) both map to bucket 10 at cap=16.
        // Java HashMap keeps insertion order within a bucket.
        use crate::variations::{Sclip, VariationData};

        // Insert chr12 (unused) then chr1 (used) -> [chr12, chr1] -> false
        let mut var_data_1 = VariationData::default();
        let mut sc_chr12 = Sclip::new();
        sc_chr12.softp = 100;
        sc_chr12.variation.cnt = 4;
        sc_chr12.used = false;
        let mut sc_chr1 = Sclip::new();
        sc_chr1.softp = 100;
        sc_chr1.variation.cnt = 4;
        sc_chr1.used = true;
        var_data_1
            .sv_structures
            .svffus
            .insert("chr12".to_string(), vec![sc_chr12]);
        var_data_1
            .sv_structures
            .svffus
            .insert("chr1".to_string(), vec![sc_chr1]);
        assert!(!is_best_sv_used(&var_data_1, 100));

        // Insert chr1 (used) then chr12 (unused) -> [chr1, chr12] -> true
        let mut var_data_2 = VariationData::default();
        let mut sc_chr1_2 = Sclip::new();
        sc_chr1_2.softp = 100;
        sc_chr1_2.variation.cnt = 4;
        sc_chr1_2.used = true;
        let mut sc_chr12_2 = Sclip::new();
        sc_chr12_2.softp = 100;
        sc_chr12_2.variation.cnt = 4;
        sc_chr12_2.used = false;
        var_data_2
            .sv_structures
            .svffus
            .insert("chr1".to_string(), vec![sc_chr1_2]);
        var_data_2
            .sv_structures
            .svffus
            .insert("chr12".to_string(), vec![sc_chr12_2]);
        assert!(is_best_sv_used(&var_data_2, 100));
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
        assert!(
            !var_data.soft_clips_5_end.get(&400).unwrap().used,
            "StructuralVariantsProcessor.java:687-965 findsv() never marks sc5v.used = true"
        );

        // Conversely, when `findconseq(sc5v, 0)` encounters a poly-A / low-complexity
        // consensus (> 12 bp with >= 7 'A's, VariationUtils.java:125-134), it sets
        // `softClip.used = true` on the live `Sclip` reference (`sync_sclip_after_conseq`).
        let mut sc5_polya = Sclip::default();
        sc5_polya.variation.cnt = 3;
        sc5_polya.variation.hicnt = 3;
        for pos_idx in 0..15 {
            sc5_polya.nt.entry(pos_idx).or_default().insert(b'A', 3);
        }
        var_data.soft_clips_5_end.insert(450, sc5_polya);

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

        let live_polya = var_data.soft_clips_5_end.get(&450).unwrap();
        assert!(
            live_polya.used,
            "find_sv must propagate find_conseq's poly-A/low-complexity used=true to the live soft_clips_5_end entry"
        );
        assert_eq!(
            live_polya.sequence.as_deref(),
            Some(&b"AAAAAAAAAAAAAAA"[..])
        );
    }

    #[test]
    fn test_find_sv_inversion_accumulates_existing_variation_count() {
        // Ports StructuralVariantsProcessor.java:781-795, 922-936:
        // Unlike the deletion branches in `findsv` (which set `vref.cnt = 0`),
        // the 5' and 3' inversion branches call `getVariation(...)` without zeroing
        // `vref.cnt` and then accumulate via `adjCnt(vref, scv.vars)`.
        use crate::config::Config;
        use crate::reference::{LoadedRegion, Reference};
        use crate::variations::{Sclip, VariationData};
        use std::collections::HashMap;

        let chr = "chr1";
        let mut bases = vec![b'A'; 600];
        let seed_g = [b'G'; 17];
        bases[383..400].copy_from_slice(&seed_g); // 1-based 384..=400

        let mut reference = Reference::default();
        reference.loaded_regions.push(LoadedRegion {
            chr: chr.to_string(),
            sequence_start: 1,
            sequence_end: 600,
        });
        reference.seed.insert(seed_g.to_vec(), vec![384]);
        reference
            .reference_sequences
            .insert(chr.to_string(), bases.clone());
        reference.region_sequences.push(bases);

        let mut var_data = VariationData::default();
        var_data.max_read_length = 100;

        let mut sc3 = Sclip::default();
        sc3.variation.cnt = 3;
        sc3.variation.hicnt = 3;
        sc3.variation.pmean = 60.0;
        sc3.variation.qmean = 90.0;
        sc3.sequence = Some(vec![b'C'; 17]);
        var_data.soft_clips_3_end.insert(100, sc3);

        // Pre-populate the exact 3' inversion key at pos 100 with cnt = 2.
        let ins5 = complement_seq(&crate::util::join_ref(&reference, chr, 351, 400));
        let ins3 = complement_seq(&crate::util::join_ref(&reference, chr, 100, 149));
        let vn = format!("-301^{}<inv201>{}", ins5, ins3);
        let pre_vref =
            crate::variations::get_variation(&mut var_data.non_insertion_variants, 100, &vn);
        pre_vref.cnt = 2;
        pre_vref.hicnt = 2;

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

        let post_cnt = var_data
            .non_insertion_variants
            .get(&100)
            .and_then(|m| m.map.get(&vn))
            .map(|v| v.cnt)
            .unwrap_or(0);
        assert_eq!(
            post_cnt, 5,
            "find_sv inversion branch must accumulate onto existing vref.cnt (2 + 3 = 5) instead of zeroing vref.cnt"
        );
    }

    #[test]
    fn test_find_match_and_rev_short_seq_below_seed_1() {
        let reference = Reference::default();
        let short_seq = b"ACGTACGTACGT"; // 12 bp < SEED_1 (17)
        for dir in [1, -1] {
            let m = find_match(short_seq, &reference, 100, dir, 17, 3, "chr1", false);
            assert_eq!(m.base_position, 0);
            assert!(m.matched_sequence.is_empty());

            let mr = find_match_rev(short_seq, &reference, 100, dir, 17, 3, "chr1", false);
            assert_eq!(mr.base_position, 0);
            assert!(mr.matched_sequence.is_empty());
        }
    }

    #[test]
    fn test_find_match_and_rev_sseq_shrinks_to_len_le_2_without_panic() {
        // Seed length 12 (SEED_2) with a 12-bp query: during `for ii in 1..=15`,
        // `sseq` shrinks from 12 down through len == 2, len == 1, and len == 0 (break).
        let chr = "chr1";
        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
        // Reference is all 'A', while query is all 'C' (complement 'G'), so every
        // base mismatches the reference and `sseq` trims all the way to empty.
        reference
            .reference_sequences
            .insert(chr.to_string(), vec![b'A'; 200]);
        reference.seed.insert(vec![b'C'; 12], vec![50]);
        reference.seed.insert(vec![b'G'; 12], vec![50]);

        let seq12 = [b'C'; 12];
        for dir in [1, -1] {
            let m = find_match(&seq12, &reference, 50, dir, 12, 0, chr, false);
            assert_eq!(m.base_position, 0);
            assert!(m.matched_sequence.is_empty());

            let mr = find_match_rev(&seq12, &reference, 50, dir, 12, 0, chr, false);
            assert_eq!(mr.base_position, 0);
            assert!(mr.matched_sequence.is_empty());
        }
    }

    #[test]
    fn test_find_match_and_rev_trimmed_seed_recovery_both_dirs() {
        let chr = "chr1";
        let seed17 = b"ACGTACGTACGTACGTA";
        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
        let mut ref_seq = vec![b'C'; 200];
        // Place seed17 at 1-based positions 102..=118 (0-based 101..118).
        // Surrounding positions 100..=101 and 119..=120 remain 'C'.
        ref_seq[101..118].copy_from_slice(seed17);
        reference
            .reference_sequences
            .insert(chr.to_string(), ref_seq);
        reference.seed.insert(seed17.to_vec(), vec![102]);

        // 1. find_match, dir = 1:
        // seq = "TT" + seed17 (len 19, seed at i = 2 -> initial bp = 102 - 2 = 100).
        // Primary is_match_ref(mm_limit = 1) fails due to 2 leading 'T' vs 'C' mismatches.
        // At ii = 2, sseq trims "TT", bp_mut = 102, extra = "TT".
        let mut seq_fwd = b"TT".to_vec();
        seq_fwd.extend_from_slice(seed17);
        let m_fwd = find_match(&seq_fwd, &reference, 50, 1, 17, 1, chr, false);
        assert_eq!(m_fwd.base_position, 102);
        assert_eq!(m_fwd.matched_sequence, b"TT");

        // 2. find_match, dir = -1:
        // After internal reverse, seq_vec = seed17 + "TT" (seed at i = 0 -> initial bp = 102 + 19 - 0 - 1 = 120).
        // Primary is_match_ref(mm_limit = 1) fails due to 2 trailing 'T' vs 'C' mismatches at 119..=120.
        // At ii = 2, sseq trims trailing "TT", bp_mut = 118, extra = "TT".
        let mut seq_rev_input = seed17.to_vec();
        seq_rev_input.extend_from_slice(b"TT");
        seq_rev_input.reverse();
        let m_rev = find_match(&seq_rev_input, &reference, 150, -1, 17, 1, chr, false);
        assert_eq!(m_rev.base_position, 118);
        assert_eq!(m_rev.matched_sequence, b"TT");

        // 3. find_match_rev, dir = -1:
        // find_match_rev(dir = -1) complements seq without reversing.
        // Pre-complement input so complemented seq_vec = "TT" + seed17.
        let seq_inv_neg: Vec<u8> = seq_fwd
            .iter()
            .map(|&b| crate::util::complement(b))
            .collect();
        let mr_neg = find_match_rev(&seq_inv_neg, &reference, 50, -1, 17, 1, chr, false);
        assert_eq!(mr_neg.base_position, 102);
        assert_eq!(mr_neg.matched_sequence, b"TT");

        // 4. find_match_rev, dir = 1:
        // find_match_rev(dir = 1) reverses and complements seq.
        // Pre-revcomp input so revcomp'd seq_vec = seed17 + "TT".
        let mut seq_inv_pos: Vec<u8> = seed17
            .iter()
            .chain(b"TT".iter())
            .map(|&b| crate::util::complement(b))
            .collect();
        seq_inv_pos.reverse();
        let mr_pos = find_match_rev(&seq_inv_pos, &reference, 150, 1, 17, 1, chr, false);
        assert_eq!(mr_pos.base_position, 118);
        assert_eq!(mr_pos.matched_sequence, b"TT");
    }

    #[test]
    fn test_find_match_rev_failed_trimmed_seed_before_full_match_returns_empty_extra() {
        let chr = "chr1";
        // Complemented query `seq_vec` (len 21, seed_len = 17):
        // - At i = 4 (checked first in reverse loop): seed_b matches at ref pos 50,
        //   where ref has 4 matching bases at 47..=50 so the trim loop sets `extra_slice`,
        //   but the rest of ref (51..=66) is 'T' so `is_match_ref(sseq, ..., 1)` fails.
        // - At i = 0 (checked later): seed_a matches at ref pos 150, where the entire
        //   21-bp `seq_vec` matches the reference at 150..=170.
        // The full-match branch at i = 0 must return an empty `matched_sequence`.
        let seq_vec = b"AAAACCCCCCCCCCCCCCCCC";
        let seed_a = &seq_vec[0..17];
        let seed_b = &seq_vec[4..21];

        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 250,
            });
        let mut ref_seq = vec![b'T'; 250];
        // At 1-based 46..=50 (0-based 45..50), place 'A' so at ii = 1 (bp_mut = 47),
        // b0 = 'A' and b1 = 'A' match ref at 47 and 48, setting `extra_slice = "A"`,
        // while the trailing 'C's mismatch 'T' so `is_match_ref(..., 1)` fails.
        ref_seq[45..50].copy_from_slice(b"AAAAA");
        // At 1-based 150..=170 (0-based 149..170), place full `seq_vec`.
        ref_seq[149..170].copy_from_slice(seq_vec);
        reference
            .reference_sequences
            .insert(chr.to_string(), ref_seq);
        reference.seed.insert(seed_b.to_vec(), vec![50]);
        reference.seed.insert(seed_a.to_vec(), vec![150]);

        let input_for_neg: Vec<u8> = seq_vec
            .iter()
            .map(|&b| crate::util::complement(b))
            .collect();
        let m = find_match_rev(&input_for_neg, &reference, 100, -1, 17, 1, chr, false);
        assert_eq!(m.base_position, 150);
        assert!(
            m.matched_sequence.is_empty(),
            "full match after a failed trimmed seed must not leak extra_slice: {:?}",
            String::from_utf8_lossy(&m.matched_sequence)
        );
    }

    #[test]
    fn midpoint_i32_and_sv_overflow_guards() {
        assert_eq!(midpoint_i32(1_500_000_000, 1_500_000_200), 1_500_000_100);
        assert_eq!(midpoint_i32(-1_500_000_000, -1_500_000_200), -1_500_000_100);
        assert_eq!(midpoint_i32(i32::MAX, i32::MAX), i32::MAX);
        assert_eq!(midpoint_i32(i32::MIN, i32::MIN), i32::MIN);

        // VariationRealigner.java:241 — Java 32-bit signed `int +=` wraps on overflow without panicking.
        let make_mate = |mate_start: i32, mate_end: i32, mate_len: i32| Mate {
            mate_start,
            mate_end,
            mate_len,
            start: 100,
            end: 200,
            pmean: 10.0,
            qmean: 30.0,
            mqmean: 60.0,
            nm: 0.0,
        };
        let mut mates = vec![make_mate(500, 600, i32::MIN), make_mate(510, 610, i32::MIN)];
        let cluster = check_cluster(&mut mates, 100);
        assert_eq!(cluster.cnt, 2);
        assert_eq!(cluster.mate.mate_len, 0);

        // Deep cluster (8,000 mates × 1,000,000 bp) wraps to -589,934,592 / 8,000 == -73,741
        // (|mlen| < 250,000 splice threshold at StructuralVariantsProcessor.java:996), matching Java.
        let mut deep_mates: Vec<Mate> = (0..8_000)
            .map(|i| make_mate(500 + (i % 10), 600 + (i % 10), 1_000_000))
            .collect();
        let deep_cluster = check_cluster(&mut deep_mates, 100);
        assert_eq!(deep_cluster.cnt, 8_000);
        assert_eq!(deep_cluster.mate.mate_len, -73_741);
        assert!(i64::from(deep_cluster.mate.mate_len).abs() < 250_000);

        // check_pairs at > 1.07 Gbp uses midpoint_i32 without i32 sum overflow.
        let mut svr = Sclip::default();
        svr.start = 1_500_000_000;
        svr.end = 1_500_000_200;
        svr.mstart = 1_500_000_400;
        svr.mend = 1_500_000_600;
        svr.variation.cnt = 4;
        let mut cluster_vec = vec![svr];
        let mut refs: Vec<&mut Vec<Sclip>> = vec![&mut cluster_vec];
        let pd = check_pairs(1_500_000_100, 1_500_000_500, &mut refs, 100, false);
        assert_eq!(pd.pairs, 4);
        assert!(cluster_vec[0].used);

        // find_inv_disc at > 1.07 Gbp (both 5' and 3' loops) uses midpoint_i32 without overflow.
        // Separate 3' clusters by > 300 bp from 5' clusters so mark_sv in Part 1 does not
        // pre-mark svfinv3[0]/svrinv3[0] before Part 2 executes.
        let chr = "chr1";
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1_500_000_000,
                sequence_end: 1_500_003_000,
            });
        reference.region_sequences.push(vec![b'A'; 3_001]);

        let mut var_data = crate::variations::VariationData::default();
        var_data.max_read_length = 100;

        let mut invf5 = Sclip::default();
        invf5.variation.cnt = 5;
        invf5.variation.mqmean = 300.0;
        invf5.start = 1_500_000_000;
        invf5.end = 1_500_000_100;
        invf5.mstart = 1_500_000_450;
        invf5.mend = 1_500_000_500;
        var_data.sv_structures.svfinv5.push(invf5);

        let mut invr5 = Sclip::default();
        invr5.variation.cnt = 5;
        invr5.variation.mqmean = 300.0;
        invr5.start = 1_500_000_120;
        invr5.end = 1_500_000_200;
        invr5.mstart = 1_500_000_520;
        invr5.mend = 1_500_000_600;
        var_data.sv_structures.svrinv5.push(invr5);

        let mut invf3 = Sclip::default();
        invf3.variation.cnt = 5;
        invf3.variation.mqmean = 300.0;
        invf3.start = 1_500_002_450;
        invf3.end = 1_500_002_550;
        invf3.mstart = 1_500_002_150;
        invf3.mend = 1_500_002_250;
        var_data.sv_structures.svfinv3.push(invf3);

        let mut invr3 = Sclip::default();
        invr3.variation.cnt = 5;
        invr3.variation.mqmean = 300.0;
        invr3.start = 1_500_002_560;
        invr3.end = 1_500_002_640;
        invr3.mstart = 1_500_002_260;
        invr3.mend = 1_500_002_340;
        var_data.sv_structures.svrinv3.push(invr3);

        let config = Config::default();
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert(chr.to_string(), 2_000_000_000);
        find_inv_disc(&mut var_data, &mut reference, chr, &config, &chr_lengths).unwrap();
        assert!(var_data.sv_structures.svfinv5[0].used);
        assert!(var_data.sv_structures.svrinv5[0].used);
        assert!(var_data.sv_structures.svfinv3[0].used);
        assert!(var_data.sv_structures.svrinv3[0].used);
        assert!(var_data.non_insertion_variants.contains_key(&1_500_000_110));
        assert!(var_data.non_insertion_variants.contains_key(&1_500_002_255));

        // find_del_disc with non-empty splice and del.mlen == i32::MIN does not panic on .abs().
        let mut del_data = crate::variations::VariationData::default();
        del_data.max_read_length = 100;
        del_data.splice.insert("100-200".to_string());
        let mut fdel = Sclip::default();
        fdel.variation.cnt = 10;
        fdel.mlen = i32::MIN;
        fdel.start = 100;
        fdel.end = 200;
        fdel.mstart = 300; // <= end + mindist (1000) -> exits cleanly after splice check
        fdel.mend = 400;
        del_data.sv_structures.svfdel.push(fdel);

        let mut rdel = Sclip::default();
        rdel.variation.cnt = 10;
        rdel.mlen = i32::MIN;
        rdel.start = 300; // <= mend + mindist (1000) -> exits cleanly after splice check
        rdel.end = 400;
        rdel.mstart = 100;
        rdel.mend = 200;
        del_data.sv_structures.svrdel.push(rdel);

        find_del_disc(&mut del_data, &mut reference, chr, &config, &chr_lengths).unwrap();

        // Near-`i32::MAX` window/overlap saturation guards (`3 * rlen + 200`, `end + wsize`, `end + rlen`).
        assert!(is_overlap(
            i32::MAX - 10,
            i32::MAX,
            i32::MAX - 5,
            i32::MAX,
            1_000,
        ));
        let mut max_sclip_list = vec![Sclip::default()];
        let mut max_lists: Vec<&mut Vec<Sclip>> = vec![&mut max_sclip_list];
        let _ = mark_sv(
            i32::MAX - 500,
            i32::MAX - 10,
            &mut max_lists,
            i32::MAX / 2,
            false,
        );

        // `find_match` (`dir = -1`) and `find_match_rev` (`dir = 1`) must compute
        // `first_seed + seq.len() - i - 1` in `i64` without clamping `first_seed + seq.len()`
        // at `i32::MAX` (`first_seed = i32::MAX - 40, seq.len() = 50, i = 33 -> bp = i32::MAX - 24`).
        {
            let start_fm = (i32::MAX - 80) as u32;
            let mut fm_seq = vec![b'T'; 81];
            // Place a 50-bp distinct region at `i32::MAX - 73 ..= i32::MAX - 24` (offsets 7..57),
            // where the last 17 bp (`i = 33..50`, ref `i32::MAX - 40 ..= i32::MAX - 24`) is the unique seed.
            let target50 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
            assert_eq!(target50.len(), 50);
            fm_seq[7..57].copy_from_slice(target50);
            let seed17 = &target50[33..50];
            let mut fm_ref = Reference::new();
            fm_ref.loaded_regions.push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: start_fm,
                sequence_end: i32::MAX as u32,
            });
            fm_ref.region_sequences.push(fm_seq);
            fm_ref
                .seed
                .insert(seed17.to_vec(), vec![(i32::MAX - 40) as u32]);

            let query_dir_neg1: Vec<u8> = target50.iter().rev().copied().collect();
            let m_fwd = find_match(
                &query_dir_neg1,
                &fm_ref,
                i32::MAX - 10,
                -1,
                17,
                3,
                chr,
                false,
            );
            assert_eq!(m_fwd.base_position, i32::MAX - 24);

            // `find_match_rev` with `dir = 1` reverses and complements `seq` into `target50`.
            let query_rev_dir1: Vec<u8> = target50
                .iter()
                .rev()
                .map(|&b| crate::util::complement(b))
                .collect();
            let m_rev = find_match_rev(
                &query_rev_dir1,
                &fm_ref,
                i32::MAX - 10,
                1,
                17,
                3,
                chr,
                false,
            );
            assert_eq!(m_rev.base_position, i32::MAX - 24);
        }

        // `find_dup_disc` (`svrdup` branch) homopolymer shift reaching `tpe == i32::MAX`
        // must advance past `i32::MAX` via `tpe64` and terminate instead of looping forever.
        {
            let start_max = (i32::MAX - 300) as u32;
            let mut dup_seq = vec![b'A'; 301];
            let unique_seed = b"CGATCGATCGATCGATC";
            dup_seq[280..297].copy_from_slice(unique_seed);
            let mut dup_ref = Reference::new();
            dup_ref.loaded_regions.push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: start_max,
                sequence_end: i32::MAX as u32,
            });
            dup_ref.region_sequences.push(dup_seq);
            dup_ref
                .seed
                .insert(unique_seed.to_vec(), vec![(i32::MAX - 20) as u32]);

            let bp_dup = i32::MAX - 240;
            let mut dup_cluster = Sclip::default();
            dup_cluster.start = i32::MAX - 250;
            dup_cluster.end = i32::MAX - 200;
            dup_cluster.mstart = i32::MAX - 60;
            dup_cluster.mend = i32::MAX - 1;
            dup_cluster.variation.cnt = 10;
            dup_cluster.variation.mqmean = 600.0;
            dup_cluster.soft.insert(bp_dup, 10);

            let fwd_clip = b"CGATCGATCGATCGATCAA";
            let mut sclip5 = Sclip::new();
            sclip5.variation.cnt = 10;
            sclip5.variation.mqmean = 600.0;
            for (idx, &b) in fwd_clip.iter().rev().enumerate() {
                let mut nt = BTreeMap::new();
                nt.insert(b, 10);
                sclip5.nt.insert(idx as i32, nt);
            }

            let mut dup_data = crate::variations::VariationData::default();
            dup_data.max_read_length = 150;
            dup_data.sv_structures.svrdup.push(dup_cluster);
            dup_data.soft_clips_5_end.insert(bp_dup, sclip5);
            // Also seed a 3'-end clip at `i32::MAX` to confirm `tpe64 == i32::MAX + 1`
            // does NOT clamp back onto `i32::MAX` when probing `soft_clips_3_end`.
            let mut bogus_sc3 = Sclip::new();
            bogus_sc3.variation.cnt = 999;
            dup_data.soft_clips_3_end.insert(i32::MAX, bogus_sc3);

            let mut dup_chr_lens = HashMap::new();
            dup_chr_lens.insert(chr.to_string(), i32::MAX as usize);
            find_dup_disc(
                &mut dup_data,
                std::path::Path::new(""),
                &mut dup_ref,
                chr,
                &config,
                &dup_chr_lens,
            )
            .unwrap();
            assert!(dup_data.insertion_variants.contains_key(&bp_dup));
            let dup_sv = dup_data
                .non_insertion_variants
                .get(&bp_dup)
                .and_then(|vm| vm.sv.as_ref())
                .unwrap();
            assert_eq!(dup_sv.r#type, "DUP");
            assert_eq!(dup_sv.splits, 10);
        }

        // `sv::partial_pipeline` O(|Delta| log |Delta|) `new_key_log` exact set-difference
        // (`(N_1 ∪ I_1) \ (N_0 ∪ I_0)`) with a `.remove()` freelist hole in `N_0`,
        // overlapping/disjoint keys across `N_0` and `I_0`, and exception-safe restoration.
        {
            use noodles::core::Position;
            use noodles::csi::binning_index::index::reference_sequence::{
                bin::Chunk, index::LinearIndex,
            };
            use noodles::csi::binning_index::Indexer;
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
            let bam_path = std::env::temp_dir().join(format!("vr_sv_delta_log_{pid}.bam"));
            let bai_path = bam_path.with_extension("bam.bai");
            let header = noodles::sam::Header::builder()
                .add_reference_sequence(
                    chr,
                    noodles::sam::header::record::value::map::Map::<
                        noodles::sam::header::record::value::map::ReferenceSequence,
                    >::new(NonZeroUsize::new(200).unwrap()),
                )
                .build();
            {
                let file = std::fs::File::create(&bam_path).unwrap();
                let mut writer = noodles::bam::io::Writer::new(file);
                writer.write_header(&header).unwrap();
                // Two clean `11M2I15M` reads (28 bp each, single indel >= 11 bp from both ends):
                // - r1 @ 10 (ref 10..=35): plants `+CC` at 20, SNV `T` at 21, ref bases across 22..=23, SNV `T` at 26
                // - r2 @ 14 (ref 14..=39): plants `+GG` at 24 and ref bases across 14..=39
                let cigar: Cigar = [
                    Op::new(Kind::Match, 11),
                    Op::new(Kind::Insertion, 2),
                    Op::new(Kind::Match, 15),
                ]
                .into_iter()
                .collect();
                let data1: Data = [(Tag::EDIT_DISTANCE, BufValue::from(4i32))]
                    .into_iter()
                    .collect();
                let rec1 = RecordBuf::builder()
                    .set_name("r_mate_1")
                    .set_flags(noodles::sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(10usize).unwrap())
                    .set_mapping_quality(
                        noodles::sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                    )
                    .set_cigar(cigar.clone())
                    .set_sequence(Sequence::from(b"AAAAAAAAAAACCTAAAATAAAAAAAAA".to_vec()))
                    .set_quality_scores(QualityScores::from(vec![35u8; 28]))
                    .set_data(data1)
                    .build();
                let data2: Data = [(Tag::EDIT_DISTANCE, BufValue::from(2i32))]
                    .into_iter()
                    .collect();
                let rec2 = RecordBuf::builder()
                    .set_name("r_mate_2")
                    .set_flags(noodles::sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(14usize).unwrap())
                    .set_mapping_quality(
                        noodles::sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                    )
                    .set_cigar(cigar)
                    .set_sequence(Sequence::from(b"AAAAAAAAAAAGGAAAAAAAAAAAAAAA".to_vec()))
                    .set_quality_scores(QualityScores::from(vec![35u8; 28]))
                    .set_data(data2)
                    .build();
                let start_vp = writer.get_ref().virtual_position();
                writer.write_alignment_record(&header, &rec1).unwrap();
                writer.write_alignment_record(&header, &rec2).unwrap();
                let end_vp = writer.get_ref().virtual_position();
                writer.try_finish().unwrap();

                let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
                let s = Position::try_from(10usize).unwrap();
                let e = Position::try_from(39usize).unwrap();
                indexer
                    .add_record(Some((0, s, e, true)), Chunk::new(start_vp, end_vp))
                    .unwrap();
                noodles::bam::bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                    .write_index(&indexer.build(1))
                    .unwrap();
            }

            let mut pp_ref = Reference::new();
            pp_ref.loaded_regions.push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
            pp_ref.region_sequences.push(vec![b'A'; 200]);
            pp_ref
                .reference_sequences
                .insert(chr.to_string(), vec![b'A'; 200]);
            let mut pp_chr_lens = HashMap::new();
            pp_chr_lens.insert(chr.to_string(), 200usize);

            let mut pp_data = crate::variations::VariationData::default();
            // Seed N_0:
            // - 5: disjoint outside window
            // - 20: in N_0 only, read adds 20 to Delta_I -> MUST NOT be planted
            // - 21: inserted then `.remove(&21)`d to leave a freelist hole in N_0 -> read re-inserts 21 -> MUST be planted!
            // - 22: in both N_0 and I_0 -> MUST NOT be planted
            let _ = crate::variations::get_variation(&mut pp_data.non_insertion_variants, 5, "T");
            let _ = crate::variations::get_variation(&mut pp_data.non_insertion_variants, 20, "T");
            let _ = crate::variations::get_variation(&mut pp_data.non_insertion_variants, 21, "T");
            assert!(pp_data.non_insertion_variants.remove(&21).is_some());
            let _ = crate::variations::get_variation(&mut pp_data.non_insertion_variants, 22, "T");
            // Seed I_0:
            // - 6: disjoint outside window
            // - 22: shared with N_0
            // - 23: in I_0 only, read adds 23 to Delta_N -> MUST NOT be planted!
            let _ = crate::variations::get_variation(&mut pp_data.insertion_variants, 6, "+TT");
            let _ = crate::variations::get_variation(&mut pp_data.insertion_variants, 22, "+TT");
            let _ = crate::variations::get_variation(&mut pp_data.insertion_variants, 23, "+TT");

            let before_union: std::collections::BTreeSet<i32> = pp_data
                .non_insertion_variants
                .keys()
                .chain(pp_data.insertion_variants.keys())
                .copied()
                .collect();

            partial_pipeline(
                &bam_path,
                &mut pp_data,
                &mut pp_ref,
                chr,
                8,
                45,
                &config,
                &pp_chr_lens,
            )
            .unwrap();

            let after_union: std::collections::BTreeSet<i32> = pp_data
                .non_insertion_variants
                .keys()
                .chain(pp_data.insertion_variants.keys())
                .copied()
                .collect();
            let expected_planted: std::collections::BTreeSet<i32> =
                after_union.difference(&before_union).copied().collect();

            assert_eq!(pp_data.merge_planted, expected_planted);
            assert!(pp_data.insertion_variants.contains_key(&20));
            assert!(!pp_data.merge_planted.contains(&20)); // N_0 only + Delta_I
            assert!(pp_data.merge_planted.contains(&21)); // N_0 freelist hole reused by Delta_N
            assert!(!pp_data.merge_planted.contains(&22)); // N_0 ∩ I_0
            assert!(pp_data.non_insertion_variants.contains_key(&23));
            assert!(!pp_data.merge_planted.contains(&23)); // I_0 only + Delta_N
            assert!(pp_data.insertion_variants.contains_key(&24));
            assert!(pp_data.merge_planted.contains(&24)); // Delta_N ∩ Delta_I
            assert!(pp_data.merge_planted.contains(&26)); // Delta_N only

            // Verify new_key_log was disarmed after Ok return:
            let _ = crate::variations::get_variation(&mut pp_data.non_insertion_variants, 180, "C");
            let _ = crate::variations::get_variation(&mut pp_data.insertion_variants, 181, "+C");
            assert!(pp_data.non_insertion_variants.take_new_key_log().is_empty());
            assert!(pp_data.insertion_variants.take_new_key_log().is_empty());

            // Verify error-path safety: second BAM in colon-separated path fails to open AFTER
            // `std::mem::take` moved all maps into `parser`. Maps and `sv_structures` must be
            // restored unconditionally, `merge_planted` untouched, and `new_key_log` disarmed.
            let mut err_data = crate::variations::VariationData::default();
            let _ = crate::variations::get_variation(&mut err_data.non_insertion_variants, 42, "G");
            let _ = crate::variations::get_variation(&mut err_data.insertion_variants, 43, "+GG");
            err_data.ref_coverage.insert(42, 9);
            err_data.splice.insert("10-20".to_string());
            err_data.sv_structures.svfdel.push(Sclip::default());
            let failing_bam_spec = format!(
                "{}:/nonexistent/missing_vr_task4_{pid}.bam",
                bam_path.display()
            );
            let err_res = partial_pipeline(
                std::path::Path::new(&failing_bam_spec),
                &mut err_data,
                &mut pp_ref,
                chr,
                18,
                35,
                &config,
                &pp_chr_lens,
            );
            assert!(err_res.is_err());
            assert!(err_data.non_insertion_variants.contains_key(&42));
            assert!(err_data.insertion_variants.contains_key(&43));
            assert_eq!(err_data.ref_coverage.get(&42), Some(&9));
            assert!(err_data.splice.contains("10-20"));
            assert_eq!(err_data.sv_structures.svfdel.len(), 1);
            assert!(err_data.merge_planted.is_empty());

            // Confirm `new_key_log` is disarmed after `Err` too:
            let _ =
                crate::variations::get_variation(&mut err_data.non_insertion_variants, 190, "T");
            let _ = crate::variations::get_variation(&mut err_data.insertion_variants, 191, "+T");
            assert!(err_data
                .non_insertion_variants
                .take_new_key_log()
                .is_empty());
            assert!(err_data.insertion_variants.take_new_key_log().is_empty());

            let _ = std::fs::remove_file(&bam_path);
            let _ = std::fs::remove_file(&bai_path);
        }
    }
}
