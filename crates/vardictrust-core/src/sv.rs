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
use crate::realign::{btree_range_inclusive, btree_range_mut_inclusive};
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
    let Some(first_mate) = mates.first() else {
        return Cluster::default();
    };
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
        let current_cluster_end = clusters.get(cur).map_or(0, |c| c.mate.mate_end);
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

        let Some(current_cluster) = clusters.get_mut(cur) else {
            continue;
        };
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

    let Some(first_cluster) = clusters.first() else {
        return Cluster::default();
    };
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

#[derive(Clone, Copy, Debug)]
enum SvListSlot {
    Finv3(usize),
    Rinv3(usize),
    Finv5(usize),
    Rinv5(usize),
    Fdel(usize),
    Rdel(usize),
    Fdup(usize),
    Rdup(usize),
    Ffus(usize, usize),
    Rfus(usize, usize),
}

#[derive(Clone, Copy, Debug)]
struct SclipSlot {
    best_cnt: i32,
    slot: SvListSlot,
}

/// Builds an upfront `SOFTP2SV` index mapping each non-zero `softp` to the `SclipSlot`
/// of the first highest-`variation.cnt` candidate across `SVStructures`, matching
/// `VariationRealigner.java:143-205` (`putToSoft2SV` + stable descending sort on `varsCount`).
fn build_softp2sv_map(
    svs: &crate::variations::SVStructures,
) -> rustc_hash::FxHashMap<i32, SclipSlot> {
    let mut map = rustc_hash::FxHashMap::default();
    if svs.svfinv3.is_empty()
        && svs.svrinv3.is_empty()
        && svs.svfinv5.is_empty()
        && svs.svrinv5.is_empty()
        && svs.svfdel.is_empty()
        && svs.svrdel.is_empty()
        && svs.svfdup.is_empty()
        && svs.svrdup.is_empty()
        && svs.svffus.is_empty()
        && svs.svrfus.is_empty()
    {
        return map;
    }

    let put = |map: &mut rustc_hash::FxHashMap<i32, SclipSlot>, sv: &Sclip, slot: SvListSlot| {
        if sv.softp == 0 || sv.variation.cnt < 0 {
            return;
        }
        match map.entry(sv.softp) {
            std::collections::hash_map::Entry::Occupied(mut e) => {
                if sv.variation.cnt > e.get().best_cnt {
                    e.insert(SclipSlot {
                        best_cnt: sv.variation.cnt,
                        slot,
                    });
                }
            }
            std::collections::hash_map::Entry::Vacant(v) => {
                v.insert(SclipSlot {
                    best_cnt: sv.variation.cnt,
                    slot,
                });
            }
        }
    };

    for (i, sv) in svs.svfinv3.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Finv3(i));
    }
    for (i, sv) in svs.svrinv3.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Rinv3(i));
    }
    for (i, sv) in svs.svfinv5.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Finv5(i));
    }
    for (i, sv) in svs.svrinv5.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Rinv5(i));
    }
    for (i, sv) in svs.svfdel.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Fdel(i));
    }
    for (i, sv) in svs.svrdel.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Rdel(i));
    }
    for (i, sv) in svs.svfdup.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Fdup(i));
    }
    for (i, sv) in svs.svrdup.iter().enumerate() {
        put(&mut map, sv, SvListSlot::Rdup(i));
    }

    let scan_fus = |map: &mut rustc_hash::FxHashMap<i32, SclipSlot>,
                    fus_map: &indexmap::IndexMap<String, Vec<Sclip>>,
                    is_fwd: bool| {
        let mk_slot = |ei: usize, si: usize| {
            if is_fwd {
                SvListSlot::Ffus(ei, si)
            } else {
                SvListSlot::Rfus(ei, si)
            }
        };
        if fus_map.len() <= 1 {
            for (ei, (_, sv_list)) in fus_map.iter().enumerate() {
                for (si, sv) in sv_list.iter().enumerate() {
                    put(map, sv, mk_slot(ei, si));
                }
            }
            return;
        }
        if !fus_map.values().any(|sv_list| {
            sv_list.iter().any(|sv| {
                sv.softp != 0
                    && sv.variation.cnt >= 0
                    && map
                        .get(&sv.softp)
                        .is_none_or(|cur| sv.variation.cnt > cur.best_cnt)
            })
        }) {
            return;
        }
        for key in crate::javamap::java_hashmap_key_order(fus_map.keys()) {
            if let Some((ei, _, sv_list)) = fus_map.get_full(key) {
                for (si, sv) in sv_list.iter().enumerate() {
                    put(map, sv, mk_slot(ei, si));
                }
            }
        }
    };
    scan_fus(&mut map, &svs.svffus, true);
    scan_fus(&mut map, &svs.svrfus, false);
    map
}

/// Dereferences live `.used` from `&SVStructures` via `SclipSlot` so mid-loop
/// `check_pairs` `.used = true` mutations (`AR-08`) are observed immediately.
#[inline]
fn is_softp_sv_used(
    softp2sv: &rustc_hash::FxHashMap<i32, SclipSlot>,
    svs: &crate::variations::SVStructures,
    softp: i32,
) -> bool {
    if softp == 0 {
        return false;
    }
    let Some(entry) = softp2sv.get(&softp) else {
        return false;
    };
    match entry.slot {
        SvListSlot::Finv3(i) => svs.svfinv3.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Rinv3(i) => svs.svrinv3.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Finv5(i) => svs.svfinv5.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Rinv5(i) => svs.svrinv5.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Fdel(i) => svs.svfdel.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Rdel(i) => svs.svrdel.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Fdup(i) => svs.svfdup.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Rdup(i) => svs.svrdup.get(i).is_some_and(|sv| sv.used),
        SvListSlot::Ffus(ei, si) => svs
            .svffus
            .get_index(ei)
            .and_then(|(_, list)| list.get(si))
            .is_some_and(|sv| sv.used),
        SvListSlot::Rfus(ei, si) => svs
            .svrfus
            .get_index(ei)
            .and_then(|(_, list)| list.get(si))
            .is_some_and(|sv| sv.used),
    }
}

/// Check if the SV candidate with the highest count at `softp` is used.
/// Ports `SOFTP2SV.containsKey(p5) && SOFTP2SV.get(p5).get(0).used`.
pub fn is_best_sv_used(var_data: &crate::variations::VariationData, softp: i32) -> bool {
    if softp == 0 {
        return false;
    }
    let softp2sv = build_softp2sv_map(&var_data.sv_structures);
    is_softp_sv_used(&softp2sv, &var_data.sv_structures, softp)
}

#[inline]
fn fmt_del_key<'a>(
    dellen: i32,
    key_buf: &'a mut crate::types::StackBuf<16>,
    key_heap: &'a mut String,
) -> &'a str {
    use std::fmt::Write as _;
    *key_buf = crate::types::StackBuf::new();
    if write!(key_buf, "-{dellen}").is_ok() {
        key_buf.as_str()
    } else {
        key_heap.clear();
        let _ = write!(key_heap, "-{dellen}");
        key_heap.as_str()
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
    if !verbose && seq.len() < seed_len {
        return Match {
            base_position: 0,
            matched_sequence: Vec::new(),
        };
    }

    let mut stack_buf = [0u8; 512];
    let mut heap_buf: Vec<u8>;
    let seq_slice: &[u8] = if dir == -1 {
        if let Some(dst) = stack_buf.get_mut(..seq.len()) {
            for (d, &s) in dst.iter_mut().zip(seq.iter().rev()) {
                *d = s;
            }
            dst
        } else {
            heap_buf = seq.to_vec();
            heap_buf.reverse();
            heap_buf.as_slice()
        }
    } else {
        seq
    };

    if verbose {
        crate::verbose_eprintln!(
            "    Working Match {} {} {} SEED: {}",
            position,
            String::from_utf8_lossy(seq_slice),
            dir,
            seed_len
        );
    }

    let mut extra_slice: &[u8] = &[];
    if seq_slice.len() < seed_len {
        return Match {
            base_position: 0,
            matched_sequence: Vec::new(),
        };
    }

    for i in (0..=seq_slice.len() - seed_len).rev() {
        let Some(seed) = seq_slice.get(i..i + seed_len) else {
            continue;
        };
        if let Some([first_seed_u32]) = reference.seed.get(seed) {
            let first_seed64 = i64::from(*first_seed_u32);
            let mut bp64 = if dir == 1 {
                first_seed64 - (i as i64)
            } else {
                first_seed64 + (seq_slice.len() as i64) - (i as i64) - 1
            };

            if i32::try_from(bp64).ok().is_some_and(|bp| {
                crate::realign::is_match_ref(seq_slice, reference, chr, bp, dir, mm_limit)
            }) {
                let mut extra = extra_slice.to_vec();
                let mut mm = if dir == -1 { -1 } else { 0 };
                loop {
                    if let Some(ch) = crate::util::char_at(seq_slice, mm) {
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
                        String::from_utf8_lossy(seq_slice),
                        String::from_utf8_lossy(&extra)
                    );
                }
                return Match {
                    base_position: bp_out,
                    matched_sequence: extra,
                };
            } else {
                // For complex indels, allowing some mismatches at the end up to 15bp or 20% length
                let mut sseq: &[u8] = seq_slice;
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
                        extra_slice = seq_slice.get(..ii).unwrap_or(seq_slice);
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
                        extra_slice = seq_slice
                            .get(seq_slice.len().saturating_sub(ii)..)
                            .unwrap_or(seq_slice);
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
                            String::from_utf8_lossy(seq_slice),
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
    if !verbose && seq.len() < seed_len {
        return Match {
            base_position: 0,
            matched_sequence: Vec::new(),
        };
    }

    let mut stack_buf = [0u8; 512];
    let mut heap_buf: Vec<u8>;
    let seq_slice: &[u8] = if let Some(dst) = stack_buf.get_mut(..seq.len()) {
        if dir == 1 {
            for (d, &s) in dst.iter_mut().zip(seq.iter().rev()) {
                *d = crate::util::complement(s);
            }
        } else {
            for (d, &s) in dst.iter_mut().zip(seq.iter()) {
                *d = crate::util::complement(s);
            }
        }
        dst
    } else {
        heap_buf = seq.to_vec();
        if dir == 1 {
            heap_buf.reverse();
        }
        for base in heap_buf.iter_mut() {
            *base = crate::util::complement(*base);
        }
        heap_buf.as_slice()
    };

    if verbose {
        crate::verbose_eprintln!(
            "    Working MatchRev {} {} {} SEED: {}",
            position,
            String::from_utf8_lossy(seq_slice),
            dir,
            seed_len
        );
    }

    let mut extra_slice: &[u8] = &[];
    if seq_slice.len() < seed_len {
        return Match {
            base_position: 0,
            matched_sequence: Vec::new(),
        };
    }

    for i in (0..=seq_slice.len() - seed_len).rev() {
        let Some(seed) = seq_slice.get(i..i + seed_len) else {
            continue;
        };
        if let Some([first_seed_u32]) = reference.seed.get(seed) {
            let first_seed64 = i64::from(*first_seed_u32);
            let bp64 = if dir == 1 {
                first_seed64 + (seq_slice.len() as i64) - (i as i64) - 1
            } else {
                first_seed64 - (i as i64)
            };

            if let Ok(bp) = i32::try_from(bp64) {
                if crate::realign::is_match_ref(
                    seq_slice,
                    reference,
                    chr,
                    bp,
                    dir.saturating_neg(),
                    mm_limit,
                ) {
                    if verbose {
                        crate::verbose_eprintln!(
                            "      Found SV BP (reverse): {} BP: {} SEEDpos: {} {} {} {} {}",
                            dir,
                            bp,
                            first_seed64,
                            position,
                            String::from_utf8_lossy(seed),
                            i,
                            String::from_utf8_lossy(seq_slice)
                        );
                    }
                    return Match {
                        base_position: bp,
                        matched_sequence: Vec::new(),
                    };
                }
            }
            // For complex indels, allowing some mismatches at the end up to 15bp or 20% length
            let mut sseq: &[u8] = seq_slice;
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
                    extra_slice = seq_slice.get(..ii).unwrap_or(seq_slice);
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
                    extra_slice = seq_slice
                        .get(seq_slice.len().saturating_sub(ii)..)
                        .unwrap_or(seq_slice);
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
                        String::from_utf8_lossy(seq_slice),
                        String::from_utf8_lossy(extra_slice)
                    );
                }
                if let Ok(bp_mut) = i32::try_from(bp_mut64) {
                    if crate::realign::is_match_ref(
                        sseq,
                        reference,
                        chr,
                        bp_mut,
                        dir.saturating_neg(),
                        1,
                    ) {
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
                cnt = cnt.saturating_add(1);
                pairs = pairs.saturating_add(sv_r.variation.cnt);
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
    let mut pairs: i32 = 0;
    let mut cnt: i32 = 0;

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
                cnt = cnt.saturating_add(1);
                pairs = pairs.saturating_add(sv_r.variation.cnt);
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

#[inline]
fn single_bam_paths(bam_path: &std::path::Path) -> Vec<String> {
    let s = bam_path.to_string_lossy();
    if s.is_empty() {
        Vec::new()
    } else {
        vec![s.into_owned()]
    }
}

#[inline]
fn saturated_chr_len(chromosomes_lengths: &HashMap<String, usize>, chr: &str) -> u32 {
    u32::try_from(*chromosomes_lengths.get(chr).unwrap_or(&0)).unwrap_or(u32::MAX)
}

fn seed_bam_piled_spans(chr: &str, reference: &Reference) -> Vec<(i32, i32)> {
    let raw: Vec<(i32, i32)> = reference
        .loaded_regions
        .iter()
        .filter(|lr| lr.chr == chr)
        .map(|lr| {
            (
                i32::try_from(lr.sequence_start).unwrap_or(i32::MAX),
                i32::try_from(lr.sequence_end).unwrap_or(i32::MAX),
            )
        })
        .collect();
    crate::merge_i32_intervals(raw)
}

#[inline]
fn spans_contain(spans: &[(i32, i32)], start: i32, end: i32) -> bool {
    let idx = spans.partition_point(|&(lo, _)| lo <= start);
    idx.checked_sub(1)
        .and_then(|cand| spans.get(cand))
        .is_some_and(|&(lo, hi)| lo <= start && end <= hi)
}

fn record_bam_piled_span(spans: &mut Vec<(i32, i32)>, start: i32, end: i32) {
    spans.push((start, end));
    *spans = crate::merge_i32_intervals(std::mem::take(spans));
}

#[expect(
    clippy::too_many_arguments,
    reason = "mate FASTA check threads region bounds, read-reach window, and reference context"
)]
fn ensure_mate_fasta_loaded(
    chr: &str,
    reg_start: i32,
    reg_end: i32,
    extension: u32,
    need_lo: i32,
    need_hi: i32,
    chromosomes_lengths: &HashMap<String, usize>,
    config: &Config,
    reference: &mut Reference,
) -> Result<(), VardictError> {
    let s = reg_start.max(1);
    let e = reg_end.max(s);
    let chr_len_u32 = saturated_chr_len(chromosomes_lengths, chr);
    let chr_len = i32::try_from(chr_len_u32).unwrap_or(i32::MAX);
    let (check_lo, check_hi) = if chr_len > 0 && s <= chr_len {
        let lo = need_lo.min(s).clamp(1, chr_len);
        let hi = need_hi.max(e).clamp(lo, chr_len);
        (lo, hi)
    } else {
        let lo = need_lo.min(s).max(1);
        let hi = need_hi.max(e).max(lo);
        (lo, hi)
    };
    if !is_loaded(chr, check_lo as u32, check_hi as u32, reference) {
        let req_pad = u32::try_from(s.saturating_sub(check_lo).max(check_hi.saturating_sub(e)))
            .unwrap_or(extension);
        let mod_region = crate::regions::Region {
            chr: chr.to_string(),
            start: Pos1(s as u32),
            end: Pos1(e as u32),
            gene: String::new(),
            amplicon: None,
        };
        get_reference(
            &mod_region,
            extension.max(req_pad),
            chromosomes_lengths,
            config,
            reference,
        )?;
    }
    Ok(())
}

/// Minimum soft-clip consensus length (`12 bp`, matching `Configuration.SEED_2` in
/// `StructuralVariantsProcessor.findMatchRev`) required before `sv_lazy_mate_pileup`
/// loads or probes an inversion mate window (`M-DOCUMENTED-MAGIC`).
const MIN_INV_SCLIP_LEN: usize = 12;

/// Upper ceiling (`32 KiB`) on explicit sequence materialization under
/// `CompatProfile::sv_respect_min_len` (`AstraZeneca-NGS/VarDictJava#328`,
/// `AstraZeneca-NGS/VarDict#157`) so huge `-L` settings (e.g. `-L 1000000`) keep
/// bounded per-tile heap usage (`AGENTS.md §10`).
pub(crate) const MAX_EXPLICIT_SV_SEQ_BP: i32 = 32_768;

fn ensure_explicit_sv_span_loaded(
    chr: &str,
    lo: i32,
    hi: i32,
    chromosomes_lengths: &HashMap<String, usize>,
    config: &Config,
    reference: &mut Reference,
) -> Result<bool, VardictError> {
    if lo < 1 || hi < lo {
        return Ok(false);
    }
    let lo_u32 = lo as u32;
    let hi_u32 = hi as u32;
    if (!is_loaded(chr, lo_u32, hi_u32, reference)
        || reference.get_base(chr, lo_u32).is_none()
        || reference.get_base(chr, hi_u32).is_none())
        && !config.reference.is_empty()
    {
        let mod_region = crate::regions::Region {
            chr: chr.to_string(),
            start: Pos1(lo_u32),
            end: Pos1(hi_u32),
            gene: String::new(),
            amplicon: None,
        };
        get_reference(&mod_region, 100, chromosomes_lengths, config, reference)?;
    }
    Ok(is_loaded(chr, lo_u32, hi_u32, reference)
        && reference.get_base(chr, lo_u32).is_some()
        && reference.get_base(chr, hi_u32).is_some())
}

fn ensure_explicit_del_sv_span_loaded(
    chr: &str,
    del_pos: i32,
    dellen: i32,
    sv_min: i32,
    chromosomes_lengths: &HashMap<String, usize>,
    config: &Config,
    reference: &mut Reference,
) -> Result<(), VardictError> {
    if config.compat.sv_respect_min_len
        && 0 < dellen
        && dellen <= sv_min.min(MAX_EXPLICIT_SV_SEQ_BP)
    {
        let _ = ensure_explicit_sv_span_loaded(
            chr,
            del_pos.max(1),
            del_pos.saturating_add(dellen),
            chromosomes_lengths,
            config,
            reference,
        )?;
    }
    Ok(())
}

fn max_splice_and_del_ref_span(var_data: &crate::variations::VariationData) -> i32 {
    let mut max_splice: i32 = 0;
    for sp in &var_data.splice {
        if let Some((lhs, rhs)) = sp.rsplit_once('-') {
            if let (Ok(s), Ok(e)) = (lhs.parse::<i64>(), rhs.parse::<i64>()) {
                let span = i32::try_from(e.saturating_sub(s).max(0)).unwrap_or(i32::MAX);
                max_splice = max_splice.max(span);
            }
        }
    }
    let mut max_del: i32 = 0;
    for del_map in var_data.position_to_deletion_count.values() {
        for key in del_map.keys() {
            if let Some(rest) = key.strip_prefix('-') {
                let digits_len = rest.bytes().take_while(u8::is_ascii_digit).count();
                if let Some(num_str) = rest.get(..digits_len) {
                    if let Ok(dlen) = num_str.parse::<i32>() {
                        max_del = max_del.max(dlen);
                    }
                }
            }
        }
    }
    max_splice.saturating_add(max_del)
}

#[derive(Clone, Copy)]
enum MateFastaMode {
    /// `find_del_paths`: ensure `[mstart, mend]` FASTA is loaded in both lazy fallback and legacy paths.
    MateSpan,
    /// `find_dup_disc_paths`: ensure `flank_pos ± 150` is loaded before mate-window pileup.
    DupFlank(i32),
    /// `find_inv_disc_paths`: `sv_lazy_mate_pileup`-only mate pileup (`!bam_paths.is_empty()`).
    LazyInvDiscOnly,
}

fn ensure_flank_base_loaded(
    chr: &str,
    flank_pos: i32,
    extension: u32,
    chromosomes_lengths: &HashMap<String, usize>,
    config: &Config,
    reference: &mut Reference,
) -> Result<(), VardictError> {
    if reference.get_base(chr, flank_pos as u32).is_none() {
        let mod_region = crate::regions::Region {
            chr: chr.to_string(),
            start: Pos1(flank_pos.saturating_sub(150).max(1) as u32),
            end: Pos1(flank_pos.saturating_add(150).max(1) as u32),
            gene: String::new(),
            amplicon: None,
        };
        get_reference(
            &mod_region,
            extension,
            chromosomes_lengths,
            config,
            reference,
        )?;
    }
    Ok(())
}

#[expect(
    clippy::too_many_arguments,
    reason = "shared SV mate-window helper threads region, BAM paths, pileup state, and FASTA context"
)]
fn ensure_sv_mate_window_piled(
    chr: &str,
    mstart: i32,
    mend: i32,
    extension: u32,
    mode: MateFastaMode,
    bam_paths: &[String],
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    bam_piled_spans: &mut Vec<(i32, i32)>,
) -> Result<(), VardictError> {
    let mate_lo = mstart.saturating_sub(200).max(1);
    let mate_hi = mend.saturating_add(200);
    if config.compat.sv_lazy_mate_pileup {
        if let MateFastaMode::DupFlank(flank_pos) = mode {
            ensure_flank_base_loaded(
                chr,
                flank_pos,
                extension,
                chromosomes_lengths,
                config,
                reference,
            )?;
        }
        if matches!(mode, MateFastaMode::LazyInvDiscOnly) && bam_paths.is_empty() {
            return Ok(());
        }
        if !spans_contain(bam_piled_spans, mstart.max(1), mend.max(1)) {
            let eff_max_rl = var_data.max_read_length.max(0);
            ensure_mate_fasta_loaded(
                chr,
                mstart,
                mend,
                extension,
                mate_lo.saturating_sub(eff_max_rl),
                mate_hi.saturating_add(eff_max_rl),
                chromosomes_lengths,
                config,
                reference,
            )?;
            partial_pipeline_paths(
                bam_paths,
                var_data,
                reference,
                chr,
                mate_lo,
                mate_hi,
                config,
                chromosomes_lengths,
            )?;
            record_bam_piled_span(bam_piled_spans, mate_lo, mate_hi);
        } else if matches!(
            mode,
            MateFastaMode::MateSpan | MateFastaMode::LazyInvDiscOnly
        ) {
            ensure_mate_fasta_loaded(
                chr,
                mstart,
                mend,
                extension,
                mstart.max(1),
                mend.max(1),
                chromosomes_lengths,
                config,
                reference,
            )?;
        }
    } else if !matches!(mode, MateFastaMode::LazyInvDiscOnly)
        && !is_loaded(chr, mstart as u32, mend as u32, reference)
    {
        match mode {
            MateFastaMode::MateSpan => {
                let mod_region = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(mstart.max(1) as u32),
                    end: Pos1(mend.max(1) as u32),
                    gene: String::new(),
                    amplicon: None,
                };
                get_reference(
                    &mod_region,
                    extension,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
            }
            MateFastaMode::DupFlank(flank_pos) => {
                ensure_flank_base_loaded(
                    chr,
                    flank_pos,
                    extension,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
            }
            MateFastaMode::LazyInvDiscOnly => {}
        }
        partial_pipeline_paths(
            bam_paths,
            var_data,
            reference,
            chr,
            mate_lo,
            mate_hi,
            config,
            chromosomes_lengths,
        )?;
    }
    Ok(())
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
    let paths = single_bam_paths(bam_path);
    let mut bam_piled_spans = if config.compat.sv_lazy_mate_pileup {
        seed_bam_piled_spans(chr, reference)
    } else {
        Vec::new()
    };
    find_inv_sub_paths(
        svref_idx,
        dir,
        side,
        var_data,
        &paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )
}

fn find_inv_sub_paths(
    svref_idx: usize,
    dir: i32,
    side: Side,
    var_data: &mut crate::variations::VariationData,
    bam_paths: &[String],
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    bam_piled_spans: &mut Vec<(i32, i32)>,
) -> Result<(), VardictError> {
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
    if n == 0 {
        return Ok(());
    }
    let minr = i32::try_from(config.min_reads).unwrap_or(i32::MAX);
    let verbose = config.verbose;
    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_sets(&config.adaptor);
    let conseq_dir = if config.compat.adaptor_all_consensus_dirs {
        if dir == 1 {
            3
        } else {
            5
        }
    } else {
        0
    };
    let max_read_length = var_data.max_read_length;
    let mut max_extra_ref_span = if config.compat.sv_lazy_mate_pileup {
        max_splice_and_del_ref_span(var_data)
    } else {
        0
    };
    for idx in 0..n {
        let (cnt, mstart, mend, start, end, mut softp) = {
            let Some(inv) = (match svref_idx {
                0 => var_data.sv_structures.svfinv5.get(idx),
                1 => var_data.sv_structures.svrinv5.get(idx),
                2 => var_data.sv_structures.svfinv3.get(idx),
                3 => var_data.sv_structures.svrinv3.get(idx),
                _ => {
                    return Err(VardictError::Other(format!(
                        "internal error: invalid inversion candidate index {svref_idx}"
                    )));
                }
            }) else {
                continue;
            };
            if inv.used || inv.variation.cnt < minr {
                continue;
            }
            (
                inv.variation.cnt,
                inv.mstart,
                inv.mend,
                inv.start,
                inv.end,
                best_soft_pos(&inv.soft).unwrap_or(0),
            )
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

        let mate_lo = mstart.saturating_sub(200).max(1);
        let mate_hi = mend.saturating_add(200);
        if config.compat.sv_lazy_mate_pileup {
            let eff_max_rl = max_read_length.max(0);
            let reach = eff_max_rl.saturating_add(max_extra_ref_span);
            let widened_lo = mate_lo.saturating_sub(reach);
            let widened_hi = mate_hi.saturating_add(reach);
            let probe_overlaps = if softp != 0 {
                widened_lo <= softp && softp <= widened_hi
            } else if max_read_length > 0 {
                let sp = if dir == 1 { end } else { start };
                let a = sp.saturating_add(dir);
                let b = sp.saturating_add(2i32.saturating_mul(max_read_length).saturating_mul(dir));
                let probe_lo = a.min(b);
                let probe_hi = a.max(b);
                widened_lo <= probe_hi && probe_lo <= widened_hi
            } else {
                false
            };
            if !spans_contain(bam_piled_spans, mstart.max(1), mend.max(1)) && probe_overlaps {
                ensure_mate_fasta_loaded(
                    chr,
                    mstart,
                    mend,
                    500,
                    widened_lo,
                    widened_hi,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
                partial_pipeline_paths(
                    bam_paths,
                    var_data,
                    reference,
                    chr,
                    mate_lo,
                    mate_hi,
                    config,
                    chromosomes_lengths,
                )?;
                record_bam_piled_span(bam_piled_spans, mate_lo, mate_hi);
                max_extra_ref_span = max_splice_and_del_ref_span(var_data);
            }
        } else if !is_loaded(chr, mstart as u32, mend as u32, reference) {
            let mod_region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(mstart.max(1) as u32),
                end: Pos1(mend.max(1) as u32),
                gene: String::new(),
                amplicon: None,
            };
            get_reference(&mod_region, 500, chromosomes_lengths, config, reference)?;
            partial_pipeline_paths(
                bam_paths,
                var_data,
                reference,
                chr,
                mate_lo,
                mate_hi,
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
            seq = crate::variations::util::find_conseq_sets(
                scv_mut,
                conseq_dir,
                &adaptor_forward,
                &adaptor_reverse,
                config.compat.sclip_majority_count_first,
            );
            if seq.is_empty()
                || (config.compat.sv_lazy_mate_pileup && seq.len() < MIN_INV_SCLIP_LEN)
            {
                continue;
            }
            let scv_var = scv_mut.variation.clone();

            if config.compat.sv_lazy_mate_pileup {
                ensure_mate_fasta_loaded(
                    chr,
                    mstart,
                    mend,
                    500,
                    mstart.max(1),
                    mend.max(1),
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
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
                let seq_check = crate::variations::util::find_conseq_sets(
                    scv_mut,
                    conseq_dir,
                    &adaptor_forward,
                    &adaptor_reverse,
                    config.compat.sclip_majority_count_first,
                );
                if seq_check.is_empty()
                    || (config.compat.sv_lazy_mate_pileup && seq_check.len() < MIN_INV_SCLIP_LEN)
                {
                    continue;
                }
                let scv_var = scv_mut.variation.clone();

                if config.compat.sv_lazy_mate_pileup {
                    ensure_mate_fasta_loaded(
                        chr,
                        mstart,
                        mend,
                        500,
                        mstart.max(1),
                        mend.max(1),
                        chromosomes_lengths,
                        config,
                        reference,
                    )?;
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
            let Some(inv) = (match svref_idx {
                0 => var_data.sv_structures.svfinv5.get(idx),
                1 => var_data.sv_structures.svrinv5.get(idx),
                2 => var_data.sv_structures.svfinv3.get(idx),
                3 => var_data.sv_structures.svrinv3.get(idx),
                _ => {
                    return Err(VardictError::Other(format!(
                        "internal error: invalid inversion candidate index {svref_idx}"
                    )));
                }
            }) else {
                continue;
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
            if config.compat.sv_lazy_mate_pileup
                && !spans_contain(bam_piled_spans, mstart.max(1), mend.max(1))
            {
                let eff_max_rl = max_read_length.max(0);
                ensure_mate_fasta_loaded(
                    chr,
                    mstart,
                    mend,
                    500,
                    mate_lo.saturating_sub(eff_max_rl),
                    mate_hi.saturating_add(eff_max_rl),
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
                partial_pipeline_paths(
                    bam_paths,
                    var_data,
                    reference,
                    chr,
                    mate_lo,
                    mate_hi,
                    config,
                    chromosomes_lengths,
                )?;
                record_bam_piled_span(bam_piled_spans, mate_lo, mate_hi);
            }
            let len = gap.saturating_add(1);
            let sv_flank = 50; // Configuration.SVFLANK
            let sv_min = i32::try_from(config.sv_min_len).unwrap_or(i32::MAX);
            let explicit_inv = if config.compat.sv_respect_min_len {
                len < sv_min.min(MAX_EXPLICIT_SV_SEQ_BP)
                    && ensure_explicit_sv_span_loaded(
                        chr,
                        softp,
                        bp,
                        chromosomes_lengths,
                        config,
                        reference,
                    )?
            } else {
                len.saturating_sub(2 * sv_flank) <= 0
            };
            let ins = if explicit_inv {
                complement_seq(&crate::util::join_ref(reference, chr, softp, bp))
            } else {
                let (f5_start, f3_end) = if config.compat.sv_respect_min_len {
                    (
                        bp.saturating_sub(sv_flank - 1).max(softp),
                        softp.saturating_add(sv_flank - 1).min(bp),
                    )
                } else {
                    (
                        bp.saturating_sub(sv_flank - 1),
                        softp.saturating_add(sv_flank - 1),
                    )
                };
                let ins5 = complement_seq(&crate::util::join_ref(reference, chr, f5_start, bp));
                let ins3 = complement_seq(&crate::util::join_ref(reference, chr, softp, f3_end));
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
                variation.tile_planted = true;
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
                sv.splits = sv.splits.saturating_add(scv_variation.cnt);
                sv.pairs = sv.pairs.saturating_add(cnt);
                sv.clusters = sv.clusters.saturating_add(1);
            }

            if dir == -1 {
                if let Some(ref_base) = reference.get_base(chr, softp as u32) {
                    let mut ch_buf = [0u8; 4];
                    let ref_base_str: &str = (ref_base as char).encode_utf8(&mut ch_buf);
                    if let Some(entry) = var_data.non_insertion_variants.get_mut(&softp) {
                        if let Some(ref_var) = entry.map.get_mut(ref_base_str) {
                            ref_var.cnt = ref_var.cnt.saturating_sub(scv_variation.cnt);
                            ref_var.hicnt = ref_var.hicnt.saturating_sub(scv_variation.hicnt);
                            ref_var.locnt = ref_var.locnt.saturating_sub(scv_variation.locnt);
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

            let filtered_owned: Vec<String>;
            let bams_opt: Option<&[String]> = if bam_paths.is_empty() {
                None
            } else if bam_paths.iter().all(|s| !s.is_empty()) {
                Some(bam_paths)
            } else {
                filtered_owned = bam_paths
                    .iter()
                    .filter(|s| !s.is_empty())
                    .cloned()
                    .collect();
                if filtered_owned.is_empty() {
                    None
                } else {
                    Some(&filtered_owned[..])
                }
            };
            let chr_len = saturated_chr_len(chromosomes_lengths, chr);

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
            if let Some(inv_mut) = match svref_idx {
                0 => var_data.sv_structures.svfinv5.get_mut(idx),
                1 => var_data.sv_structures.svrinv5.get_mut(idx),
                2 => var_data.sv_structures.svfinv3.get_mut(idx),
                3 => var_data.sv_structures.svrinv3.get_mut(idx),
                _ => {
                    return Err(VardictError::Other(format!(
                        "internal error: invalid inversion candidate index {svref_idx}"
                    )));
                }
            } {
                inv_mut.used = true;
            }

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
    let paths = single_bam_paths(bam_path);
    let mut bam_piled_spans = if config.compat.sv_lazy_mate_pileup {
        seed_bam_piled_spans(chr, reference)
    } else {
        Vec::new()
    };
    find_inv_paths(
        var_data,
        &paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )
}

fn find_inv_paths(
    var_data: &mut crate::variations::VariationData,
    bam_paths: &[String],
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    bam_piled_spans: &mut Vec<(i32, i32)>,
) -> Result<(), VardictError> {
    find_inv_sub_paths(
        0,
        1,
        Side::Five,
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        bam_piled_spans,
    )?;
    find_inv_sub_paths(
        1,
        -1,
        Side::Five,
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        bam_piled_spans,
    )?;
    find_inv_sub_paths(
        2,
        1,
        Side::Three,
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        bam_piled_spans,
    )?;
    find_inv_sub_paths(
        3,
        -1,
        Side::Three,
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        bam_piled_spans,
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
    let paths = single_bam_paths(bam_path);
    let mut bam_piled_spans = if config.compat.sv_lazy_mate_pileup {
        seed_bam_piled_spans(chr, reference)
    } else {
        Vec::new()
    };
    find_del_paths(
        var_data,
        &paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )
}

fn find_del_paths(
    var_data: &mut crate::variations::VariationData,
    bam_paths: &[String],
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    bam_piled_spans: &mut Vec<(i32, i32)>,
) -> Result<(), VardictError> {
    if var_data.sv_structures.svfdel.is_empty() && var_data.sv_structures.svrdel.is_empty() {
        return Ok(());
    }
    let minr = i32::try_from(config.min_reads).unwrap_or(i32::MAX);
    let sv_min = i32::try_from(config.sv_min_len).unwrap_or(i32::MAX);
    let verbose = config.verbose;
    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_sets(&config.adaptor);
    let (conseq_dir5, conseq_dir3) = if config.compat.adaptor_all_consensus_dirs {
        (5, 3)
    } else {
        (0, 0)
    };
    let mut seq_stack = [0u8; 512];
    let mut seq_heap = Vec::new();
    let mut key_buf = crate::types::StackBuf::<16>::new();
    let mut key_heap = String::new();

    // 5' deletion softclips
    let n_svfdel = var_data.sv_structures.svfdel.len();
    for idx in 0..n_svfdel {
        let (vars_count, start, end, mstart, mend, softp_cand, del_variation) = {
            let Some(del) = var_data.sv_structures.svfdel.get(idx) else {
                continue;
            };
            if del.used || del.variation.cnt < minr {
                continue;
            }
            let softp = best_soft_pos(&del.soft).unwrap_or(0);
            (
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                softp,
                del.variation.clone(),
            )
        };

        if softp_cand != 0 {
            let (seq_len, scv_var) = {
                let Some(scv) = var_data.soft_clips_3_end.get_mut(&softp_cand) else {
                    continue;
                };
                if scv.used {
                    continue;
                }
                let s = crate::variations::util::find_conseq_slice(
                    scv,
                    conseq_dir3,
                    &adaptor_forward,
                    &adaptor_reverse,
                    config.compat.sclip_majority_count_first,
                );
                if s.len() < 12 {
                    continue;
                }
                if let Some(dst) = seq_stack.get_mut(..s.len()) {
                    dst.copy_from_slice(s);
                } else {
                    seq_heap.clear();
                    seq_heap.extend_from_slice(s);
                }
                (s.len(), scv.variation.clone())
            };
            let seq: &[u8] = seq_stack.get(..seq_len).unwrap_or(&seq_heap);

            ensure_sv_mate_window_piled(
                chr,
                mstart,
                mend,
                300,
                MateFastaMode::MateSpan,
                bam_paths,
                var_data,
                reference,
                config,
                chromosomes_lengths,
                bam_piled_spans,
            )?;

            let match_res = find_match(seq, reference, softp_cand, 1, 17, 3, chr, verbose);
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
            ensure_explicit_del_sv_span_loaded(
                chr,
                softp,
                dellen,
                sv_min,
                chromosomes_lengths,
                config,
                reference,
            )?;
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
            ensure_explicit_del_sv_span_loaded(
                chr,
                softp,
                dellen,
                sv_min,
                chromosomes_lengths,
                config,
                reference,
            )?;
            let del_key = fmt_del_key(dellen, &mut key_buf, &mut key_heap);

            {
                let variation = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    softp,
                    del_key,
                );
                variation.reset_counts_compat(config.compat.reconcile_strand_and_hq_sums);
                variation.tile_planted = true;
            }

            {
                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    softp,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs = sv.pairs.saturating_add(vars_count);
                sv.splits = sv.splits.saturating_add(scv_var.cnt);
                sv.clusters = sv.clusters.saturating_add(1);
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
                    del_key,
                );
                crate::variations::util::adj_cnt(variation_mut, &scv_var, None);

                if let Some(ref_base) = reference.get_base(chr, softp as u32) {
                    let mut ch_buf = [0u8; 4];
                    let ref_base_str: &str = (ref_base as char).encode_utf8(&mut ch_buf);
                    if let Some(entry) = var_data.non_insertion_variants.get_mut(&softp) {
                        if let Some(ref_var) = entry.map.get_mut(ref_base_str) {
                            ref_var.cnt = ref_var.cnt.saturating_sub(scv_var.cnt);
                            ref_var.hicnt = ref_var.hicnt.saturating_sub(scv_var.hicnt);
                            ref_var.locnt = ref_var.locnt.saturating_sub(scv_var.locnt);
                            ref_var.pmean -= scv_var.pmean;
                            ref_var.qmean -= scv_var.qmean;
                            ref_var.mqmean -= scv_var.mqmean;
                            ref_var.nm -= scv_var.nm;
                            ref_var.sub_dir(true, scv_var.get_dir(true));
                            ref_var.sub_dir(false, scv_var.get_dir(false));
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
                    del_key,
                );
                crate::variations::util::adj_cnt(variation_mut2, &tv, None);
            }

            if let Some(del_mut) = var_data.sv_structures.svfdel.get_mut(idx) {
                del_mut.used = true;
            }

            let mut svrdel_mut_refs = [&mut var_data.sv_structures.svrdel];
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
            let limit = 3i32.saturating_mul(max_read_length);
            let range_lo = end.saturating_sub(3);
            let hi64 = i64::from(end) + i64::from(limit) - 1;
            let (eff_lo, eff_hi) = if hi64 < i64::from(range_lo) {
                (1i32, 0i32)
            } else {
                (range_lo, hi64.min(i64::from(i32::MAX)) as i32)
            };
            for (&i, scv) in
                btree_range_mut_inclusive(&mut var_data.soft_clips_3_end, eff_lo, eff_hi)
            {
                if scv.used {
                    continue;
                }
                if !(i >= end.saturating_sub(3)
                    && i.saturating_sub(end) < 3i32.saturating_mul(max_read_length))
                {
                    continue;
                }
                let seq = crate::variations::util::find_conseq_slice(
                    scv,
                    conseq_dir3,
                    &adaptor_forward,
                    &adaptor_reverse,
                    config.compat.sclip_majority_count_first,
                );
                if seq.len() < 12 {
                    continue;
                }
                let match_res = find_match(seq, reference, i, 1, 17, 3, chr, verbose);
                let mut bp = match_res.base_position;
                if bp == 0 {
                    let match_res_lenient = find_match(seq, reference, i, 1, 12, 0, chr, verbose);
                    bp = match_res_lenient.base_position;
                }
                if bp == 0 {
                    continue;
                }
                if !(bp.saturating_sub(i) > 30 && is_overlap(i, bp, end, mstart, max_read_length)) {
                    continue;
                }
                bp = bp.saturating_sub(1);
                matched_j = Some((i, bp, scv.variation.clone()));
                break;
            }

            if let Some((softp, bp, scv_var)) = matched_j {
                let dellen = bp.saturating_sub(softp).saturating_add(1);
                ensure_explicit_del_sv_span_loaded(
                    chr,
                    softp,
                    dellen,
                    sv_min,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
                let del_key = fmt_del_key(dellen, &mut key_buf, &mut key_heap);
                {
                    let variation = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        softp,
                        del_key,
                    );
                    variation.reset_counts_compat(config.compat.reconcile_strand_and_hq_sums);
                    variation.tile_planted = true;
                }

                {
                    let sv = crate::variations::VariationMap::get_sv(
                        &mut var_data.non_insertion_variants,
                        softp,
                    );
                    sv.r#type = "DEL".to_string();
                    sv.pairs = sv.pairs.saturating_add(vars_count);
                    sv.splits = sv.splits.saturating_add(scv_var.cnt);
                    sv.clusters = sv.clusters.saturating_add(1);
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
                        del_key,
                    );
                    crate::variations::util::adj_cnt(variation_mut, &scv_var, None);
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
                        del_key,
                    );
                    crate::variations::util::adj_cnt(variation_mut2, &tv, None);
                }

                if let Some(del_mut) = var_data.sv_structures.svfdel.get_mut(idx) {
                    del_mut.used = true;
                }

                let mut svrdel_mut_refs = [&mut var_data.sv_structures.svrdel];
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
        let (vars_count, start, end, mstart, mend, softp_cand, del_variation) = {
            let Some(del) = var_data.sv_structures.svrdel.get(idx) else {
                continue;
            };
            if del.used || del.variation.cnt < minr {
                continue;
            }
            let softp = best_soft_pos(&del.soft).unwrap_or(0);
            (
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                softp,
                del.variation.clone(),
            )
        };

        if softp_cand != 0 {
            let (seq_len, scv_var) = {
                let Some(scv) = var_data.soft_clips_5_end.get_mut(&softp_cand) else {
                    continue;
                };
                if scv.used {
                    continue;
                }
                let s = crate::variations::util::find_conseq_slice(
                    scv,
                    conseq_dir5,
                    &adaptor_forward,
                    &adaptor_reverse,
                    config.compat.sclip_majority_count_first,
                );
                if s.len() < 12 {
                    continue;
                }
                if let Some(dst) = seq_stack.get_mut(..s.len()) {
                    dst.copy_from_slice(s);
                } else {
                    seq_heap.clear();
                    seq_heap.extend_from_slice(s);
                }
                (s.len(), scv.variation.clone())
            };
            let seq: &[u8] = seq_stack.get(..seq_len).unwrap_or(&seq_heap);

            ensure_sv_mate_window_piled(
                chr,
                mstart,
                mend,
                300,
                MateFastaMode::MateSpan,
                bam_paths,
                var_data,
                reference,
                config,
                chromosomes_lengths,
                bam_piled_spans,
            )?;

            let match_res = find_match(seq, reference, softp_cand, -1, 17, 3, chr, verbose);
            let mut bp = match_res.base_position;
            if bp == 0 {
                let match_res_lenient =
                    find_match(seq, reference, softp_cand, -1, 12, 0, chr, verbose);
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
            ensure_explicit_del_sv_span_loaded(
                chr,
                bp,
                dellen,
                sv_min,
                chromosomes_lengths,
                config,
                reference,
            )?;
            let del_key = fmt_del_key(dellen, &mut key_buf, &mut key_heap);

            {
                let variation = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    del_key,
                );
                variation.reset_counts_compat(config.compat.reconcile_strand_and_hq_sums);
                variation.tile_planted = true;
            }

            {
                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    bp,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs = sv.pairs.saturating_add(vars_count);
                sv.splits = sv.splits.saturating_add(scv_var.cnt);
                sv.clusters = sv.clusters.saturating_add(1);
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
                    del_key,
                );
                crate::variations::util::adj_cnt(variation_mut, &scv_var, None);
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
                    del_key,
                );
                crate::variations::util::adj_cnt(variation_mut2, &tv, None);
            }

            if let Some(del_mut) = var_data.sv_structures.svrdel.get_mut(idx) {
                del_mut.used = true;
            }

            let mut svfdel_mut_refs = [&mut var_data.sv_structures.svfdel];
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
            let limit = 3i32.saturating_mul(max_read_length);
            let range_hi = start.saturating_add(3);
            let lo64 = i64::from(start) - i64::from(limit) + 1;
            let (eff_lo, eff_hi) = if lo64 > i64::from(range_hi) {
                (1i32, 0i32)
            } else {
                (lo64.max(i64::from(i32::MIN)) as i32, range_hi)
            };
            for (&i, scv) in
                btree_range_mut_inclusive(&mut var_data.soft_clips_5_end, eff_lo, eff_hi)
            {
                if scv.used {
                    continue;
                }
                if !(i <= start.saturating_add(3)
                    && start.saturating_sub(i) < 3i32.saturating_mul(max_read_length))
                {
                    continue;
                }
                let seq = crate::variations::util::find_conseq_slice(
                    scv,
                    conseq_dir5,
                    &adaptor_forward,
                    &adaptor_reverse,
                    config.compat.sclip_majority_count_first,
                );
                if seq.len() < 12 {
                    continue;
                }
                let match_res = find_match(seq, reference, i, -1, 17, 3, chr, verbose);
                let mut bp = match_res.base_position;
                if bp == 0 {
                    let match_res_lenient = find_match(seq, reference, i, -1, 12, 0, chr, verbose);
                    bp = match_res_lenient.base_position;
                }
                if bp == 0 {
                    continue;
                }
                if !(i.saturating_sub(bp) > 30 && is_overlap(bp, i, mend, start, max_read_length)) {
                    continue;
                }
                bp = bp.saturating_add(1);
                matched_j = Some((i, bp, scv.variation.clone()));
                break;
            }

            if let Some((softp_cand, bp, scv_var)) = matched_j {
                let mut softp = softp_cand.saturating_sub(1);
                let dellen = softp.saturating_sub(bp).saturating_add(1);
                ensure_explicit_del_sv_span_loaded(
                    chr,
                    bp,
                    dellen,
                    sv_min,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
                let del_key = fmt_del_key(dellen, &mut key_buf, &mut key_heap);

                {
                    let variation = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        bp,
                        del_key,
                    );
                    variation.reset_counts_compat(config.compat.reconcile_strand_and_hq_sums);
                    variation.tile_planted = true;
                }

                {
                    let sv = crate::variations::VariationMap::get_sv(
                        &mut var_data.non_insertion_variants,
                        bp,
                    );
                    sv.r#type = "DEL".to_string();
                    sv.pairs = sv.pairs.saturating_add(vars_count);
                    sv.splits = sv.splits.saturating_add(scv_var.cnt);
                    sv.clusters = sv.clusters.saturating_add(1);
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
                let bp_cov_entry = var_data.ref_coverage.entry(bp).or_insert(0);
                *bp_cov_entry = bp_cov_entry.saturating_add(scv_var.cnt);

                {
                    let variation_mut = crate::variations::get_variation(
                        &mut var_data.non_insertion_variants,
                        bp,
                        del_key,
                    );
                    crate::variations::util::adj_cnt(variation_mut, &scv_var, None);
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
                        del_key,
                    );
                    crate::variations::util::adj_cnt(variation_mut2, &tv, None);
                }

                if let Some(del_mut) = var_data.sv_structures.svrdel.get_mut(idx) {
                    del_mut.used = true;
                }

                let mut svfdel_mut_refs = [&mut var_data.sv_structures.svfdel];
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
    if var_data.soft_clips_5_end.is_empty() && var_data.soft_clips_3_end.is_empty() {
        return Ok(());
    }
    let minr = i32::try_from(config.min_reads).unwrap_or(i32::MAX);
    let sv_min = i32::try_from(config.sv_min_len).unwrap_or(i32::MAX);
    let reconcile = config.compat.reconcile_strand_and_hq_sums;
    let verbose = config.verbose;

    // Ports StructuralVariantsProcessor.java:973 fillAndSortTmpSV:
    // if (position < CURSEG.start || position > CURSEG.end) continue;
    let mut tmp5: Vec<SortPositionCount> =
        btree_range_inclusive(&var_data.soft_clips_5_end, curseg_start, curseg_end)
            .filter_map(|(&pos, scv)| {
                if scv.used || scv.variation.cnt < minr {
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

    let mut tmp3: Vec<SortPositionCount> =
        btree_range_inclusive(&var_data.soft_clips_3_end, curseg_start, curseg_end)
            .filter_map(|(&pos, scv)| {
                if scv.used || scv.variation.cnt < minr {
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

    if tmp5.is_empty() && tmp3.is_empty() {
        return Ok(());
    }

    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_sets(&config.adaptor);
    let (conseq_dir5, conseq_dir3) = if config.compat.adaptor_all_consensus_dirs {
        (5, 3)
    } else {
        (0, 0)
    };
    let softp2sv = build_softp2sv_map(&var_data.sv_structures);
    let mut key_buf = crate::types::StackBuf::<16>::new();
    let mut key_heap = String::new();

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
        if is_softp_sv_used(&softp2sv, &var_data.sv_structures, p5) {
            continue;
        }

        let (seq, sc5v_variation) = {
            let Some(sc5v) = var_data.soft_clips_5_end.get_mut(&p5) else {
                continue;
            };
            let seq = crate::variations::util::find_conseq_sets(
                sc5v,
                conseq_dir5,
                &adaptor_forward,
                &adaptor_reverse,
                config.compat.sclip_majority_count_first,
            );
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
                    let mut refs = [&mut svfdel, &mut svrdel];
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
                ensure_explicit_del_sv_span_loaded(
                    chr,
                    bp_adj,
                    dellen,
                    sv_min,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
                let del_key = fmt_del_key(dellen, &mut key_buf, &mut key_heap);

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp_adj,
                    del_key,
                );
                vref.reset_counts_compat(reconcile);
                vref.tile_planted = true;

                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    bp_adj,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs = sv.pairs.saturating_add(pairs);
                sv.splits = sv.splits.saturating_add(cnt5);
                sv.clusters = sv.clusters.saturating_add(1);

                let bp_cov = var_data.ref_coverage.get(&bp_adj).copied().unwrap_or(0);
                let pair_split_cov = pairs.saturating_add(sc5v_variation.cnt);
                if bp_cov < pair_split_cov {
                    var_data.ref_coverage.insert(bp_adj, pair_split_cov);
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
                    del_key,
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
                    del_key,
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

            let span = bp.saturating_sub(p5_adj).saturating_add(1);
            let extra = String::from_utf8_lossy(&match_res_rev.matched_sequence);
            let explicit_inv = config.compat.sv_respect_min_len
                && span < sv_min.min(MAX_EXPLICIT_SV_SEQ_BP)
                && ensure_explicit_sv_span_loaded(
                    chr,
                    p5_adj,
                    bp,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
            let vn = if explicit_inv {
                let tins = complement_seq(&crate::util::join_ref(reference, chr, p5_adj, bp));
                format!("-{}^{}{}", span, tins, extra)
            } else if config.compat.sv_respect_min_len {
                let ins5 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    bp.saturating_sub(49).max(p5_adj),
                    bp,
                ));
                let ins3 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    p5_adj,
                    p5_adj.saturating_add(49).min(bp),
                ));
                let mid = span
                    .saturating_sub(ins5.len() as i32)
                    .saturating_sub(ins3.len() as i32);
                format!("-{}^{}<inv{}>{}{}", span, ins5, mid, ins3, extra)
            } else {
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
                let mid = span
                    .saturating_sub(ins5.len() as i32)
                    .saturating_sub(ins3.len() as i32);
                if mid <= 0 {
                    let tins = complement_seq(&crate::util::join_ref(reference, chr, p5_adj, bp));
                    format!("-{}^{}{}", span, tins, extra)
                } else {
                    format!("-{}^{}<inv{}>{}{}", span, ins5, mid, ins3, extra)
                }
            };

            let vref =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p5_adj, &vn);
            vref.tile_planted = true;

            let sv = crate::variations::VariationMap::get_sv(
                &mut var_data.non_insertion_variants,
                p5_adj,
            );
            sv.r#type = "INV".to_string();
            sv.splits = sv.splits.saturating_add(cnt5);

            let vref_mut =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p5_adj, &vn);
            crate::variations::util::adj_cnt(vref_mut, &sc5v_variation, None);

            let p5_cov = var_data.ref_coverage.get(&p5_adj).copied().unwrap_or(0);
            var_data
                .ref_coverage
                .insert(p5_adj, p5_cov.saturating_add(cnt5));
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
        if is_softp_sv_used(&softp2sv, &var_data.sv_structures, p3) {
            continue;
        }

        let (seq, sc3v_variation) = {
            let Some(sc3v) = var_data.soft_clips_3_end.get_mut(&p3) else {
                continue;
            };
            let seq = crate::variations::util::find_conseq_sets(
                sc3v,
                conseq_dir3,
                &adaptor_forward,
                &adaptor_reverse,
                config.compat.sclip_majority_count_first,
            );
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
                    let mut refs = [&mut svfdel, &mut svrdel];
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
                ensure_explicit_del_sv_span_loaded(
                    chr,
                    p3_adj,
                    dellen,
                    sv_min,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
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
                ensure_explicit_del_sv_span_loaded(
                    chr,
                    p3_adj,
                    dellen,
                    sv_min,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
                let del_key = fmt_del_key(dellen, &mut key_buf, &mut key_heap);

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    p3_adj,
                    del_key,
                );
                vref.reset_counts_compat(reconcile);
                vref.tile_planted = true;

                let sv = crate::variations::VariationMap::get_sv(
                    &mut var_data.non_insertion_variants,
                    p3_adj,
                );
                sv.r#type = "DEL".to_string();
                sv.pairs = sv.pairs.saturating_add(pairs);
                sv.splits = sv.splits.saturating_add(cnt3);
                sv.clusters = sv.clusters.saturating_add(1);

                if !var_data.ref_coverage.contains_key(&p3_adj) {
                    var_data
                        .ref_coverage
                        .insert(p3_adj, pairs.saturating_add(sc3v_variation.cnt));
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
                    del_key,
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
                    del_key,
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

            let span = bp.saturating_sub(p3_adj).saturating_add(1);
            let extra = String::from_utf8_lossy(&match_res_rev.matched_sequence);
            let explicit_inv = config.compat.sv_respect_min_len
                && span < sv_min.min(MAX_EXPLICIT_SV_SEQ_BP)
                && ensure_explicit_sv_span_loaded(
                    chr,
                    p3_adj,
                    bp,
                    chromosomes_lengths,
                    config,
                    reference,
                )?;
            let vn = if explicit_inv {
                let tins = complement_seq(&crate::util::join_ref(reference, chr, p3_adj, bp));
                format!("-{}^{}{}", span, extra, tins)
            } else if config.compat.sv_respect_min_len {
                let ins5 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    bp.saturating_sub(49).max(p3_adj),
                    bp,
                ));
                let ins3 = complement_seq(&crate::util::join_ref(
                    reference,
                    chr,
                    p3_adj,
                    p3_adj.saturating_add(49).min(bp),
                ));
                let mid = span
                    .saturating_sub(ins5.len() as i32)
                    .saturating_sub(ins3.len() as i32);
                format!("-{}^{}{}<inv{}>{}", span, extra, ins5, mid, ins3)
            } else {
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
                let mid = span.saturating_sub(100);
                if mid <= 0 {
                    let tins = complement_seq(&crate::util::join_ref(reference, chr, p3_adj, bp));
                    format!("-{}^{}{}", span, extra, tins)
                } else {
                    format!("-{}^{}{}<inv{}>{}", span, extra, ins5, mid, ins3)
                }
            };

            let vref =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p3_adj, &vn);
            vref.tile_planted = true;

            let sv = crate::variations::VariationMap::get_sv(
                &mut var_data.non_insertion_variants,
                p3_adj,
            );
            sv.r#type = "INV".to_string();
            sv.splits = sv.splits.saturating_add(cnt3);

            let vref_mut =
                crate::variations::get_variation(&mut var_data.non_insertion_variants, p3_adj, &vn);
            crate::variations::util::adj_cnt(vref_mut, &sc3v_variation, None);

            let p3_cov = var_data.ref_coverage.get(&p3_adj).copied().unwrap_or(0);
            var_data
                .ref_coverage
                .insert(p3_adj, p3_cov.saturating_add(cnt3));
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
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        out.push(crate::util::complement(b));
    }
    String::from_utf8(out).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

fn tv_adj_vals(tv: &mut Variation, pairs_data: PairsData) {
    tv.hicnt = pairs_data.pairs;
    tv.rev = pairs_data.pairs.saturating_sub(tv.fwd);
    tv.pmean = pairs_data.pmean;
    tv.qmean = pairs_data.qmean;
    tv.mqmean = pairs_data.qmean_mapping;
    tv.nm = pairs_data.nm;
}

/// Find deletion structural variants from discordant pairs only.
///
/// Ports `StructuralVariantsProcessor.findDELdisc`.
///
/// # Errors
///
/// Returns [`VardictError`] if reference FASTA loading fails for a candidate breakpoint.
pub fn find_del_disc(
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let mut bam_piled_spans = Vec::new();
    find_del_disc_paths(
        var_data,
        &[],
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )
}

#[expect(
    clippy::neg_cmp_op_on_partial_ord,
    reason = "Java-faithful negated partial-ord comparisons"
)]
fn find_del_disc_paths(
    var_data: &mut crate::variations::VariationData,
    bam_paths: &[String],
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    bam_piled_spans: &mut Vec<(i32, i32)>,
) -> Result<(), VardictError> {
    let verbose = config.verbose;
    let mindist = 8i32.saturating_mul(var_data.max_read_length);
    let min_reads_plus_5 = i64::from(config.min_reads) + 5;
    let sv_min = i32::try_from(config.sv_min_len).unwrap_or(i32::MAX);
    let reconcile = config.compat.reconcile_strand_and_hq_sums;
    let mut key_buf = crate::types::StackBuf::<16>::new();
    let mut key_heap = String::new();

    // svfdel (5' candidate deletions)
    let n_svfdel = var_data.sv_structures.svfdel.len();
    for idx in 0..n_svfdel {
        let (vars_count, _start, end, mstart, mend, softp, mlen) = {
            let Some(del) = var_data.sv_structures.svfdel.get(idx) else {
                continue;
            };
            if del.used || i64::from(del.variation.cnt) < min_reads_plus_5 {
                continue;
            }
            if !var_data.splice.is_empty() && i64::from(del.mlen).abs() < 250_000 {
                continue;
            }
            if del.mstart <= del.end.saturating_add(mindist) {
                continue;
            }
            if del.variation.cnt <= 0 || !(del.variation.mqmean / (del.variation.cnt as f64) > 35.0)
            {
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
            if config.compat.sv_respect_min_len && mlen < sv_min {
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
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                bp,
                mlen,
            )
        };

        if config.compat.sv_lazy_mate_pileup && !bam_paths.is_empty() {
            ensure_sv_mate_window_piled(
                chr,
                mstart,
                mend,
                300,
                MateFastaMode::LazyInvDiscOnly,
                bam_paths,
                var_data,
                reference,
                config,
                chromosomes_lengths,
                bam_piled_spans,
            )?;
        }

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

        let del_key = fmt_del_key(mlen, &mut key_buf, &mut key_heap);
        let vref =
            crate::variations::get_variation(&mut var_data.non_insertion_variants, softp, del_key);
        vref.reset_counts_compat(reconcile);
        vref.tile_planted = true;

        let sv =
            crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, softp);
        sv.r#type = "DEL".to_string();
        sv.splits = sv.splits.saturating_add(
            var_data
                .soft_clips_3_end
                .get(&end.saturating_add(1))
                .map(|s| s.variation.cnt)
                .unwrap_or(0),
        );
        sv.splits = sv.splits.saturating_add(
            var_data
                .soft_clips_5_end
                .get(&mstart)
                .map(|s| s.variation.cnt)
                .unwrap_or(0),
        );
        sv.pairs = sv.pairs.saturating_add(vars_count);
        sv.clusters = sv.clusters.saturating_add(1);

        let twice_vars = 2i32.saturating_mul(vars_count);
        let mut tv = Variation::new();
        tv.cnt = twice_vars;
        tv.hicnt = twice_vars;
        tv.fwd = vars_count;
        tv.rev = vars_count;
        if let Some(del) = var_data.sv_structures.svfdel.get(idx) {
            tv.qmean = 2.0 * del.variation.qmean;
            tv.pmean = 2.0 * del.variation.pmean;
            tv.mqmean = 2.0 * del.variation.mqmean;
            tv.nm = 2.0 * del.variation.nm;
        }
        let vref_mut =
            crate::variations::get_variation(&mut var_data.non_insertion_variants, softp, del_key);
        crate::variations::util::adj_cnt(vref_mut, &tv, None);

        if !var_data.ref_coverage.contains_key(&softp) {
            var_data.ref_coverage.insert(softp, twice_vars);
        }
        if let Some(del_mut) = var_data.sv_structures.svfdel.get_mut(idx) {
            del_mut.used = true;
        }

        let mut svrdel_mut_refs = [&mut var_data.sv_structures.svrdel];
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
        let (vars_count, start, end, mstart, mend, softp, del_softp, mlen) = {
            let Some(del) = var_data.sv_structures.svrdel.get(idx) else {
                continue;
            };
            if del.used || i64::from(del.variation.cnt) < min_reads_plus_5 {
                continue;
            }
            if !var_data.splice.is_empty() && i64::from(del.mlen).abs() < 250_000 {
                continue;
            }
            if del.start <= del.mend.saturating_add(mindist) {
                continue;
            }
            if del.variation.cnt <= 0 || !(del.variation.mqmean / (del.variation.cnt as f64) > 35.0)
            {
                continue;
            }
            let mlen = del
                .start
                .saturating_sub(del.mend)
                .saturating_sub(var_data.max_read_length / del.variation.cnt.saturating_add(1));
            if !(mlen > 0 && mlen > mindist) {
                continue;
            }
            if config.compat.sv_respect_min_len && mlen < sv_min {
                continue;
            }
            let bp = del.mend.saturating_add(
                (var_data.max_read_length / del.variation.cnt.saturating_add(1)) / 2,
            );
            (
                del.variation.cnt,
                del.start,
                del.end,
                del.mstart,
                del.mend,
                bp,
                del.softp,
                mlen,
            )
        };

        if config.compat.sv_lazy_mate_pileup && !bam_paths.is_empty() {
            ensure_sv_mate_window_piled(
                chr,
                mstart,
                mend,
                300,
                MateFastaMode::LazyInvDiscOnly,
                bam_paths,
                var_data,
                reference,
                config,
                chromosomes_lengths,
                bam_piled_spans,
            )?;
        }

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

        let del_key = fmt_del_key(mlen, &mut key_buf, &mut key_heap);
        let vref =
            crate::variations::get_variation(&mut var_data.non_insertion_variants, softp, del_key);
        vref.reset_counts_compat(reconcile);
        vref.tile_planted = true;

        let sv =
            crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, softp);
        sv.r#type = "DEL".to_string();
        sv.splits = sv.splits.saturating_add(
            var_data
                .soft_clips_3_end
                .get(&mend.saturating_add(1))
                .map(|s| s.variation.cnt)
                .unwrap_or(0),
        );
        sv.splits = sv.splits.saturating_add(
            var_data
                .soft_clips_5_end
                .get(&start)
                .map(|s| s.variation.cnt)
                .unwrap_or(0),
        );
        sv.pairs = sv.pairs.saturating_add(vars_count);
        sv.clusters = sv.clusters.saturating_add(1);

        if del_softp != 0 {
            if let Some(sc5) = var_data.soft_clips_5_end.get_mut(&del_softp) {
                sc5.used = true;
            }
        }

        let twice_vars = 2i32.saturating_mul(vars_count);
        let mut tv = Variation::new();
        tv.cnt = twice_vars;
        tv.hicnt = twice_vars;
        tv.fwd = vars_count;
        tv.rev = vars_count;
        if let Some(del) = var_data.sv_structures.svrdel.get(idx) {
            tv.qmean = 2.0 * del.variation.qmean;
            tv.pmean = 2.0 * del.variation.pmean;
            tv.mqmean = 2.0 * del.variation.mqmean;
            tv.nm = 2.0 * del.variation.nm;
        }
        let vref_mut =
            crate::variations::get_variation(&mut var_data.non_insertion_variants, softp, del_key);
        crate::variations::util::adj_cnt(vref_mut, &tv, None);

        if !var_data.ref_coverage.contains_key(&softp) {
            var_data.ref_coverage.insert(softp, twice_vars);
        }
        let start_cov = var_data.ref_coverage.get(&start).copied().unwrap_or(0);
        let bp_cov = var_data.ref_coverage.get(&softp).copied().unwrap_or(0);
        if start_cov > 0 && bp_cov < start_cov {
            var_data.ref_coverage.insert(softp, start_cov);
        }
        if let Some(del_mut) = var_data.sv_structures.svrdel.get_mut(idx) {
            del_mut.used = true;
        }

        let mod_region_mate = crate::regions::Region {
            chr: chr.to_string(),
            start: Pos1(mstart.saturating_sub(100).max(1) as u32),
            end: Pos1(mend.saturating_add(100).max(1) as u32),
            gene: String::new(),
            amplicon: None,
        };
        if !config.compat.sv_lazy_mate_pileup
            || !is_loaded(
                chr,
                mod_region_mate.start.0,
                mod_region_mate.end.0,
                reference,
            )
        {
            get_reference(
                &mod_region_mate,
                200,
                chromosomes_lengths,
                config,
                reference,
            )?;
        }

        let mut svfdel_mut_refs = [&mut var_data.sv_structures.svfdel];
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
pub fn find_inv_disc(
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    let mut bam_piled_spans = Vec::new();
    find_inv_disc_paths(
        var_data,
        &[],
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )
}

#[expect(
    clippy::neg_cmp_op_on_partial_ord,
    clippy::nonminimal_bool,
    reason = "Java-faithful negated partial-ord comparisons"
)]
fn find_inv_disc_paths(
    var_data: &mut crate::variations::VariationData,
    bam_paths: &[String],
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    bam_piled_spans: &mut Vec<(i32, i32)>,
) -> Result<(), VardictError> {
    let min_reads_plus_5 = i64::from(config.min_reads) + 5;
    let sv_min = i32::try_from(config.sv_min_len).unwrap_or(i32::MAX);
    let verbose = config.verbose;
    let max_read_length = var_data.max_read_length;
    let sv_flank = 50; // Configuration.SVFLANK

    // Part 1: svfinv5 and svrinv5
    let n_svfinv5 = var_data.sv_structures.svfinv5.len();
    let n_svrinv5 = var_data.sv_structures.svrinv5.len();
    for i in 0..n_svfinv5 {
        let (cnt, me, ms, end, start, nm, pmean, qmean, Qmean) = {
            let Some(invf) = var_data.sv_structures.svfinv5.get(i) else {
                continue;
            };
            if invf.used {
                continue;
            }
            (
                invf.variation.cnt,
                invf.mend,
                invf.mstart,
                invf.end,
                invf.start,
                invf.variation.nm,
                invf.variation.pmean,
                invf.variation.qmean,
                invf.variation.mqmean,
            )
        };

        if cnt <= 0 || !(Qmean / (cnt as f64) > 35.0) {
            continue;
        }

        for j in 0..n_svrinv5 {
            let (rcnt, rstart, rms, rme, rnm, rpmean, rqmean, rQmean) = {
                let Some(invr) = var_data.sv_structures.svrinv5.get(j) else {
                    continue;
                };
                if invr.used {
                    continue;
                }
                (
                    invr.variation.cnt,
                    invr.start,
                    invr.mstart,
                    invr.mend,
                    invr.variation.nm,
                    invr.variation.pmean,
                    invr.variation.qmean,
                    invr.variation.mqmean,
                )
            };

            if rcnt <= 0 || !(rQmean / (rcnt as f64) > 35.0) {
                continue;
            }
            if !(i64::from(cnt) + i64::from(rcnt) > min_reads_plus_5) {
                continue;
            }

            if is_overlap(end, me, rstart, rms, max_read_length) {
                let bp = midpoint_i32(end, rstart).saturating_abs();
                let pe = midpoint_i32(me, rms).saturating_abs();
                let len = pe.saturating_sub(bp).saturating_add(1);
                if config.compat.sv_respect_min_len && len < sv_min {
                    continue;
                }

                if config.compat.sv_lazy_mate_pileup && !bam_paths.is_empty() {
                    ensure_sv_mate_window_piled(
                        chr,
                        ms,
                        me,
                        500,
                        MateFastaMode::LazyInvDiscOnly,
                        bam_paths,
                        var_data,
                        reference,
                        config,
                        chromosomes_lengths,
                        bam_piled_spans,
                    )?;
                    ensure_sv_mate_window_piled(
                        chr,
                        rms,
                        rme,
                        500,
                        MateFastaMode::LazyInvDiscOnly,
                        bam_paths,
                        var_data,
                        reference,
                        config,
                        chromosomes_lengths,
                        bam_piled_spans,
                    )?;
                }

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
                if config.compat.sv_respect_min_len && reference.get_base(chr, bp as u32).is_none()
                {
                    ensure_flank_base_loaded(chr, bp, 300, chromosomes_lengths, config, reference)?;
                }

                let explicit_inv = if config.compat.sv_respect_min_len {
                    len < sv_min.min(MAX_EXPLICIT_SV_SEQ_BP)
                        && ensure_explicit_sv_span_loaded(
                            chr,
                            bp,
                            pe,
                            chromosomes_lengths,
                            config,
                            reference,
                        )?
                } else {
                    len.saturating_sub(2 * sv_flank) <= 0
                };
                let ins = if explicit_inv {
                    complement_seq(&crate::util::join_ref(reference, chr, bp, pe))
                } else {
                    let (f5_end, f3_start) = if config.compat.sv_respect_min_len {
                        (
                            bp.saturating_add(sv_flank - 1).min(pe),
                            pe.saturating_sub(sv_flank - 1).max(bp),
                        )
                    } else {
                        (
                            bp.saturating_add(sv_flank - 1),
                            pe.saturating_sub(sv_flank - 1),
                        )
                    };
                    let ins5 = complement_seq(&crate::util::join_ref(reference, chr, bp, f5_end));
                    let ins3 = complement_seq(&crate::util::join_ref(reference, chr, f3_start, pe));
                    format!("{}<inv{}>{}", ins3, len.saturating_sub(2 * sv_flank), ins5)
                };

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    &format!("-{}^{}", len, ins),
                );

                if let Some(invf_mut) = var_data.sv_structures.svfinv5.get_mut(i) {
                    invf_mut.used = true;
                }
                if let Some(invr_mut) = var_data.sv_structures.svrinv5.get_mut(j) {
                    invr_mut.used = true;
                }
                vref.tile_planted = true;
                vref.pstd = true;
                vref.qstd = true;

                let mut tmp = Variation::new();
                let sum_cnt = cnt.saturating_add(rcnt);
                tmp.cnt = sum_cnt;
                tmp.hicnt = sum_cnt;
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
                sv.pairs = sv.pairs.saturating_add(cnt);
                sv.splits = sv.splits.saturating_add(
                    var_data
                        .soft_clips_5_end
                        .get(&start)
                        .map_or(0, |s| s.variation.cnt),
                );
                sv.splits = sv.splits.saturating_add(
                    var_data
                        .soft_clips_5_end
                        .get(&ms)
                        .map_or(0, |s| s.variation.cnt),
                );
                sv.clusters = sv.clusters.saturating_add(1);

                if !var_data.ref_coverage.contains_key(&bp) {
                    var_data.ref_coverage.insert(bp, 2i32.saturating_mul(cnt));
                }

                let mut svfinv3_mut_refs = [
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
        let (cnt, ms, me, end, nm, pmean, qmean, Qmean) = {
            let Some(invf) = var_data.sv_structures.svfinv3.get(i) else {
                continue;
            };
            if invf.used {
                continue;
            }
            // Note: StructuralVariantsProcessor.java:1239-1246 does NOT gate invf3 on
            // Qmean/cnt > DISCPAIRQUAL (only invr3 at :1259), so preserve exact valid-input
            // parity (`cnt >= 1`) while rejecting corrupted `cnt <= 0`.
            if invf.variation.cnt <= 0 {
                continue;
            }
            (
                invf.variation.cnt,
                invf.mstart,
                invf.mend,
                invf.end,
                invf.variation.nm,
                invf.variation.pmean,
                invf.variation.qmean,
                invf.variation.mqmean,
            )
        };

        for j in 0..n_svrinv3 {
            let (rcnt, rstart, rms, rme, rnm, rpmean, rqmean, rQmean) = {
                let Some(invr) = var_data.sv_structures.svrinv3.get(j) else {
                    continue;
                };
                if invr.used {
                    continue;
                }
                (
                    invr.variation.cnt,
                    invr.start,
                    invr.mstart,
                    invr.mend,
                    invr.variation.nm,
                    invr.variation.pmean,
                    invr.variation.qmean,
                    invr.variation.mqmean,
                )
            };

            if rcnt <= 0 || !(rQmean / (rcnt as f64) > 35.0) {
                continue;
            }
            if !(i64::from(cnt) + i64::from(rcnt) > min_reads_plus_5) {
                continue;
            }

            if is_overlap(me, end, rms, rstart, max_read_length) {
                let pe = midpoint_i32(end, rstart).saturating_abs();
                let bp = midpoint_i32(me, rms).saturating_abs();
                let len = pe.saturating_sub(bp).saturating_add(1);
                if config.compat.sv_respect_min_len && len < sv_min {
                    continue;
                }

                if config.compat.sv_lazy_mate_pileup && !bam_paths.is_empty() {
                    ensure_sv_mate_window_piled(
                        chr,
                        ms,
                        me,
                        500,
                        MateFastaMode::LazyInvDiscOnly,
                        bam_paths,
                        var_data,
                        reference,
                        config,
                        chromosomes_lengths,
                        bam_piled_spans,
                    )?;
                    ensure_sv_mate_window_piled(
                        chr,
                        rms,
                        rme,
                        500,
                        MateFastaMode::LazyInvDiscOnly,
                        bam_paths,
                        var_data,
                        reference,
                        config,
                        chromosomes_lengths,
                        bam_piled_spans,
                    )?;
                }

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
                if config.compat.sv_respect_min_len && reference.get_base(chr, pe as u32).is_none()
                {
                    ensure_flank_base_loaded(chr, pe, 300, chromosomes_lengths, config, reference)?;
                }

                let explicit_inv = if config.compat.sv_respect_min_len {
                    len < sv_min.min(MAX_EXPLICIT_SV_SEQ_BP)
                        && ensure_explicit_sv_span_loaded(
                            chr,
                            bp,
                            pe,
                            chromosomes_lengths,
                            config,
                            reference,
                        )?
                } else {
                    len.saturating_sub(2 * sv_flank) <= 0
                };
                let ins = if explicit_inv {
                    complement_seq(&crate::util::join_ref(reference, chr, bp, pe))
                } else {
                    let (f5_end, f3_start) = if config.compat.sv_respect_min_len {
                        (
                            bp.saturating_add(sv_flank - 1).min(pe),
                            pe.saturating_sub(sv_flank - 1).max(bp),
                        )
                    } else {
                        (
                            bp.saturating_add(sv_flank - 1),
                            pe.saturating_sub(sv_flank - 1),
                        )
                    };
                    let ins5 = complement_seq(&crate::util::join_ref(reference, chr, bp, f5_end));
                    let ins3 = complement_seq(&crate::util::join_ref(reference, chr, f3_start, pe));
                    format!("{}<inv{}>{}", ins3, len.saturating_sub(2 * sv_flank), ins5)
                };

                let vref = crate::variations::get_variation(
                    &mut var_data.non_insertion_variants,
                    bp,
                    &format!("-{}^{}", len, ins),
                );

                if let Some(invf_mut) = var_data.sv_structures.svfinv3.get_mut(i) {
                    invf_mut.used = true;
                }
                if let Some(invr_mut) = var_data.sv_structures.svrinv3.get_mut(j) {
                    invr_mut.used = true;
                }
                vref.tile_planted = true;
                vref.pstd = true;
                vref.qstd = true;

                let mut tmp = Variation::new();
                let sum_cnt = cnt.saturating_add(rcnt);
                tmp.cnt = sum_cnt;
                tmp.hicnt = sum_cnt;
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
                sv.pairs = sv.pairs.saturating_add(cnt);
                sv.splits = sv.splits.saturating_add(
                    var_data
                        .soft_clips_3_end
                        .get(&end.saturating_add(1))
                        .map_or(0, |s| s.variation.cnt),
                );
                sv.splits = sv.splits.saturating_add(
                    var_data
                        .soft_clips_3_end
                        .get(&me.saturating_add(1))
                        .map_or(0, |s| s.variation.cnt),
                );
                sv.clusters = sv.clusters.saturating_add(1);

                if !var_data.ref_coverage.contains_key(&bp) {
                    var_data.ref_coverage.insert(bp, 2i32.saturating_mul(cnt));
                }

                let mut svfinv5_mut_refs = [
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
    let paths = single_bam_paths(bam_path);
    let mut bam_piled_spans = if config.compat.sv_lazy_mate_pileup {
        seed_bam_piled_spans(chr, reference)
    } else {
        Vec::new()
    };
    find_dup_disc_paths(
        var_data,
        &paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )
}

// Java-faithful `!(a > b)` negated partial-ord comparisons — keep as-is.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn find_dup_disc_paths(
    var_data: &mut crate::variations::VariationData,
    bam_paths: &[String],
    reference: &mut Reference,
    chr: &str,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
    bam_piled_spans: &mut Vec<(i32, i32)>,
) -> Result<(), VardictError> {
    if var_data.sv_structures.svfdup.is_empty() && var_data.sv_structures.svrdup.is_empty() {
        return Ok(());
    }
    let min_reads_plus_5 = i64::from(config.min_reads) + 5;
    let sv_min = i32::try_from(config.sv_min_len).unwrap_or(i32::MAX);
    let reconcile = config.compat.reconcile_strand_and_hq_sums;
    let verbose = config.verbose;
    let (adaptor_forward, adaptor_reverse) = crate::realign::fill_adaptor_sets(&config.adaptor);
    let (conseq_dir5, conseq_dir3) = if config.compat.adaptor_all_consensus_dirs {
        (5, 3)
    } else {
        (0, 0)
    };

    /// Under `CompatProfile::canonical_sv_vcf_anchoring` (`clintval/vartovcf#192`), shift
    /// discordant/clip `<dupN>` map keys from `bp` to `bp - 1` (`pre_anchored_dup = true`)
    /// whenever `bp > 1` and `reference` holds a loaded non-`N` base at `bp - 1`, aligning
    /// `sv.rs` `<DUP>` keys with `realignlgins` (`bi = bp - 1`) for paired `SomaticMode`.
    #[inline]
    fn resolve_dup_anchor(
        reference: &Reference,
        chr: &str,
        bp: i32,
        canonical_sv_vcf_anchoring: bool,
    ) -> (i32, bool) {
        if canonical_sv_vcf_anchoring
            && bp > 1
            && u32::try_from(bp - 1)
                .ok()
                .and_then(|pos| reference.get_base(chr, pos))
                .is_some_and(|b| b != b'N' && b != b'n')
        {
            (bp - 1, true)
        } else {
            (bp, false)
        }
    }

    let n_svfdup = var_data.sv_structures.svfdup.len();
    for idx in 0..n_svfdup {
        let (vars_count, start, mut end, mstart, mend, qmean, Qmean, pmean, nm) = {
            let Some(dup) = var_data.sv_structures.svfdup.get(idx) else {
                continue;
            };
            if dup.used || i64::from(dup.variation.cnt) < min_reads_plus_5 {
                continue;
            }
            if dup.variation.cnt <= 0 || !(dup.variation.mqmean / (dup.variation.cnt as f64) > 35.0)
            {
                continue;
            }
            (
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

        let mut mlen = end
            .saturating_sub(mstart)
            .saturating_add(var_data.max_read_length / vars_count);
        let mut bp = mstart.saturating_sub((var_data.max_read_length / vars_count) / 2);
        let mut pe = end;
        let mut mark_pe = pe;

        ensure_sv_mate_window_piled(
            chr,
            mstart,
            mend,
            300,
            MateFastaMode::DupFlank(bp),
            bam_paths,
            var_data,
            reference,
            config,
            chromosomes_lengths,
            bam_piled_spans,
        )?;

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

        let best_soft = var_data
            .sv_structures
            .svfdup
            .get(idx)
            .filter(|dup| !dup.soft.is_empty())
            .map(|dup| best_soft_pos(&dup.soft));
        let mut split_refined = false;

        if let Some(best_pe_opt) = best_soft {
            if let Some(best_pe) = best_pe_opt {
                pe = best_pe;
            }
            mark_pe = pe;

            let mut tbp = {
                let Some(current_sclip3) = var_data.soft_clips_3_end.get_mut(&pe) else {
                    continue;
                };
                if current_sclip3.used {
                    continue;
                }
                cntf = current_sclip3.variation.cnt;
                qmeanf = current_sclip3.variation.qmean;
                Qmeanf = current_sclip3.variation.mqmean;
                pmeanf = current_sclip3.variation.pmean;
                nmf = current_sclip3.variation.nm;

                let seq = crate::variations::util::find_conseq_slice(
                    current_sclip3,
                    conseq_dir3,
                    &adaptor_forward,
                    &adaptor_reverse,
                    config.compat.sclip_majority_count_first,
                );
                find_match(seq, reference, bp, 1, 17, 3, chr, verbose).base_position
            };

            if tbp != 0 && tbp < pe {
                split_refined = true;
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
                mark_pe = pe;
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

        let mut emit_explicit_dup = false;
        if config.compat.sv_respect_min_len && mlen < sv_min {
            if !split_refined {
                continue;
            }
            if mlen > 0
                && mlen < MAX_EXPLICIT_SV_SEQ_BP
                && ensure_explicit_sv_span_loaded(
                    chr,
                    bp.saturating_sub(1).max(1),
                    pe,
                    chromosomes_lengths,
                    config,
                    reference,
                )?
            {
                emit_explicit_dup = true;
            }
        }

        let (plus_ins, anchor, pre_anchored_dup) = if emit_explicit_dup {
            let dup_seq = crate::util::join_ref(reference, chr, bp, pe);
            let mut s = String::with_capacity(1 + dup_seq.len());
            s.push('+');
            if bp > 1 {
                s.push_str(&dup_seq);
                (s, bp - 1, false)
            } else {
                if let Some(tail) = dup_seq.get(1..) {
                    s.push_str(tail);
                }
                if let Some(b1) = reference
                    .get_base(chr, 1)
                    .or_else(|| dup_seq.as_bytes().first().copied())
                {
                    s.push(b1 as char);
                }
                (s, 1, false)
            }
        } else {
            let (f5_end, f3_start) = if config.compat.sv_respect_min_len {
                (bp.saturating_add(49).min(pe), pe.saturating_sub(49).max(bp))
            } else {
                (bp.saturating_add(49), pe.saturating_sub(49))
            };
            let ins5 = crate::util::join_ref(reference, chr, bp, f5_end);
            let ins3 = crate::util::join_ref(reference, chr, f3_start, pe);
            let ins = format!("{}<dup{}>{}", ins5, mlen.saturating_sub(100), ins3);
            let (anc, pre_anc) =
                resolve_dup_anchor(reference, chr, bp, config.compat.canonical_sv_vcf_anchoring);
            (format!("+{ins}"), anc, pre_anc)
        };

        let vref =
            crate::variations::get_variation(&mut var_data.insertion_variants, anchor, &plus_ins);
        if !(emit_explicit_dup && vref.cnt > 0) {
            vref.reset_counts_compat(reconcile);
        }
        vref.pre_anchored_dup = pre_anchored_dup;
        vref.tile_planted = true;

        let sv =
            crate::variations::VariationMap::get_sv(&mut var_data.non_insertion_variants, anchor);
        sv.r#type = "DUP".to_string();
        sv.pairs = sv.pairs.saturating_add(vars_count);
        let dup_softp = var_data
            .sv_structures
            .svfdup
            .get(idx)
            .map_or(0, |d| d.softp);
        sv.splits = sv.splits.saturating_add(
            var_data
                .soft_clips_3_end
                .get(&dup_softp)
                .map(|s| s.variation.cnt)
                .unwrap_or(0),
        );
        sv.clusters = sv.clusters.saturating_add(1);

        let tcnt = cntr.saturating_add(cntf);
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

        let vref_mut =
            crate::variations::get_variation(&mut var_data.insertion_variants, anchor, &plus_ins);
        crate::variations::util::adj_cnt(vref_mut, &tmp, None);
        vref_mut.pre_anchored_dup = pre_anchored_dup;

        if let Some(dup_mut) = var_data.sv_structures.svfdup.get_mut(idx) {
            dup_mut.used = true;
        }

        if !var_data.ref_coverage.contains_key(&anchor) {
            var_data.ref_coverage.insert(anchor, tcnt);
        }
        let anchor_cov = var_data.ref_coverage.get(&anchor).copied().unwrap_or(0);
        let end_cov = var_data.ref_coverage.get(&end).copied().unwrap_or(0);
        if anchor_cov < end_cov {
            var_data.ref_coverage.insert(anchor, end_cov);
        }

        let mut svrdup_mut_refs = [&mut var_data.sv_structures.svrdup];
        let mark_end = if config.compat.mark_dup_sv_inclusive_end {
            mark_pe
        } else {
            pe
        };
        let (clusters, _) = mark_dup_sv(
            bp,
            mark_end,
            &mut svrdup_mut_refs,
            var_data.max_read_length,
            verbose,
        );
        {
            let sv = crate::variations::VariationMap::get_sv(
                &mut var_data.non_insertion_variants,
                anchor,
            );
            sv.clusters = sv.clusters.saturating_add(clusters);
        }
    }

    // svrdup loop
    let n_svrdup = var_data.sv_structures.svrdup.len();
    for idx in 0..n_svrdup {
        let (vars_count, start, end, mstart, mend) = {
            let Some(dup) = var_data.sv_structures.svrdup.get(idx) else {
                continue;
            };
            if dup.used || i64::from(dup.variation.cnt) < min_reads_plus_5 {
                continue;
            }
            if dup.variation.cnt <= 0 || !(dup.variation.mqmean / (dup.variation.cnt as f64) > 35.0)
            {
                continue;
            }
            (dup.variation.cnt, dup.start, dup.end, dup.mstart, dup.mend)
        };

        let mut mlen = mend
            .saturating_sub(start)
            .saturating_add(var_data.max_read_length / vars_count);
        let mut bp = start.saturating_sub((var_data.max_read_length / vars_count) / 2);
        let mut pe = (i64::from(mlen) + i64::from(bp) - 1)
            .clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
        let mut tpe64 = i64::from(pe);

        ensure_sv_mate_window_piled(
            chr,
            mstart,
            mend,
            300,
            MateFastaMode::DupFlank(pe),
            bam_paths,
            var_data,
            reference,
            config,
            chromosomes_lengths,
            bam_piled_spans,
        )?;

        let mut cntf = vars_count;
        let mut cntr = vars_count;
        let (mut qmeanf, mut Qmeanf, mut pmeanf, mut nmf, best_soft) = {
            let Some(dup) = var_data.sv_structures.svrdup.get(idx) else {
                continue;
            };
            let best = (!dup.soft.is_empty()).then(|| best_soft_pos(&dup.soft));
            (
                dup.variation.qmean,
                dup.variation.mqmean,
                dup.variation.pmean,
                dup.variation.nm,
                best,
            )
        };
        let mut qmeanr = qmeanf;
        let mut Qmeanr = Qmeanf;
        let mut pmeanr = pmeanf;
        let mut nmr = nmf;

        let mut split_refined = false;

        if let Some(best_bp_opt) = best_soft {
            if let Some(best_bp) = best_bp_opt {
                bp = best_bp;
            }

            let tbp = {
                let Some(current_sclip5) = var_data.soft_clips_5_end.get_mut(&bp) else {
                    continue;
                };
                if current_sclip5.used {
                    continue;
                }
                cntr = current_sclip5.variation.cnt;
                qmeanr = current_sclip5.variation.qmean;
                Qmeanr = current_sclip5.variation.mqmean;
                pmeanr = current_sclip5.variation.pmean;
                nmr = current_sclip5.variation.nm;

                let seq = crate::variations::util::find_conseq_slice(
                    current_sclip5,
                    conseq_dir5,
                    &adaptor_forward,
                    &adaptor_reverse,
                    config.compat.sclip_majority_count_first,
                );
                find_match(seq, reference, pe, -1, 17, 3, chr, verbose).base_position
            };

            if tbp != 0 && tbp > bp {
                split_refined = true;
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

        let mut emit_explicit_dup = false;
        if config.compat.sv_respect_min_len && mlen < sv_min {
            if !split_refined {
                continue;
            }
            if mlen > 0
                && mlen < MAX_EXPLICIT_SV_SEQ_BP
                && ensure_explicit_sv_span_loaded(
                    chr,
                    bp.saturating_sub(1).max(1),
                    pe,
                    chromosomes_lengths,
                    config,
                    reference,
                )?
            {
                emit_explicit_dup = true;
            }
        }

        let (plus_ins, anchor, pre_anchored_dup) = if emit_explicit_dup {
            let dup_seq = crate::util::join_ref(reference, chr, bp, pe);
            let mut s = String::with_capacity(1 + dup_seq.len());
            s.push('+');
            if bp > 1 {
                s.push_str(&dup_seq);
                (s, bp - 1, false)
            } else {
                if let Some(tail) = dup_seq.get(1..) {
                    s.push_str(tail);
                }
                if let Some(b1) = reference
                    .get_base(chr, 1)
                    .or_else(|| dup_seq.as_bytes().first().copied())
                {
                    s.push(b1 as char);
                }
                (s, 1, false)
            }
        } else {
            let (f5_end, f3_start) = if config.compat.sv_respect_min_len {
                (bp.saturating_add(49).min(pe), pe.saturating_sub(49).max(bp))
            } else {
                (bp.saturating_add(49), pe.saturating_sub(49))
            };
            let ins5 = crate::util::join_ref(reference, chr, bp, f5_end);
            let ins3 = crate::util::join_ref(reference, chr, f3_start, pe);
            let ins = format!("{}<dup{}>{}", ins5, mlen.saturating_sub(100), ins3);
            let (anc, pre_anc) =
                resolve_dup_anchor(reference, chr, bp, config.compat.canonical_sv_vcf_anchoring);
            (format!("+{ins}"), anc, pre_anc)
        };

        let vref =
            crate::variations::get_variation(&mut var_data.insertion_variants, anchor, &plus_ins);
        if !(emit_explicit_dup && vref.cnt > 0) {
            vref.reset_counts_compat(reconcile);
        }
        vref.pre_anchored_dup = pre_anchored_dup;
        vref.tile_planted = true;

        {
            let sv = crate::variations::VariationMap::get_sv(
                &mut var_data.non_insertion_variants,
                anchor,
            );
            sv.r#type = "DUP".to_string();
            sv.pairs = sv.pairs.saturating_add(vars_count);
            sv.splits = sv.splits.saturating_add(
                var_data
                    .soft_clips_5_end
                    .get(&bp)
                    .map(|s| s.variation.cnt)
                    .unwrap_or(0),
            );
            if let Ok(tpe_i32) = i32::try_from(tpe64) {
                sv.splits = sv.splits.saturating_add(
                    var_data
                        .soft_clips_3_end
                        .get(&tpe_i32)
                        .map(|s| s.variation.cnt)
                        .unwrap_or(0),
                );
            }
            sv.clusters = sv.clusters.saturating_add(1);
        }

        let tcnt = cntr.saturating_add(cntf);
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

        let vref_mut =
            crate::variations::get_variation(&mut var_data.insertion_variants, anchor, &plus_ins);
        crate::variations::util::adj_cnt(vref_mut, &tmp, None);
        vref_mut.pre_anchored_dup = pre_anchored_dup;

        if let Some(dup_mut) = var_data.sv_structures.svrdup.get_mut(idx) {
            dup_mut.used = true;
        }

        if !var_data.ref_coverage.contains_key(&anchor) {
            var_data.ref_coverage.insert(anchor, tcnt);
        }
        let anchor_cov = var_data.ref_coverage.get(&anchor).copied().unwrap_or(0);
        let mend_cov = var_data.ref_coverage.get(&mend).copied().unwrap_or(0);
        if anchor_cov < mend_cov {
            var_data.ref_coverage.insert(anchor, mend_cov);
        }

        let mut svfdup_mut_refs = [&mut var_data.sv_structures.svfdup];
        let (clusters, _) = mark_dup_sv(
            bp,
            pe,
            &mut svfdup_mut_refs,
            var_data.max_read_length,
            verbose,
        );
        {
            let sv = crate::variations::VariationMap::get_sv(
                &mut var_data.non_insertion_variants,
                anchor,
            );
            sv.clusters = sv.clusters.saturating_add(clusters);
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
    let paths = single_bam_paths(bam_path);
    find_all_sv_paths(
        var_data,
        &paths,
        reference,
        chr,
        chr_len,
        pileup_region,
        config,
        chromosomes_lengths,
    )
}

/// Colon-safe multi-BAM entry point for structural variant calling (`SEC-09b`).
pub fn find_all_sv_paths(
    var_data: &mut crate::variations::VariationData,
    bam_paths: &[String],
    reference: &mut Reference,
    chr: &str,
    _chr_len: u32,
    pileup_region: &crate::regions::Region,
    config: &Config,
    chromosomes_lengths: &HashMap<String, usize>,
) -> Result<(), VardictError> {
    if config.no_sv {
        return Ok(());
    }
    filter_all_sv_structures(var_data);

    let mut bam_piled_spans = if config.compat.sv_lazy_mate_pileup {
        let s = i32::try_from(pileup_region.start.0).unwrap_or(i32::MAX);
        let e = i32::try_from(pileup_region.end.0).unwrap_or(i32::MAX);
        if s <= e {
            vec![(s, e)]
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    if verbose_log(config) {
        crate::verbose_eprintln!("Start Structural Variants: DEL\n");
    }
    find_del_paths(
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )?;

    if verbose_log(config) {
        crate::verbose_eprintln!("Start Structural Variants: INV\n");
    }
    find_inv_paths(
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
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
    find_del_disc_paths(
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )?;

    if verbose_log(config) {
        crate::verbose_eprintln!("Start Structural Variants: INV discordant pairs only\n");
    }
    find_inv_disc_paths(
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
    )?;

    if verbose_log(config) {
        crate::verbose_eprintln!("Start Structural Variants: DUP discordant pairs only\n");
    }
    find_dup_disc_paths(
        var_data,
        bam_paths,
        reference,
        chr,
        config,
        chromosomes_lengths,
        &mut bam_piled_spans,
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
    let paths = single_bam_paths(bam_path);
    partial_pipeline_paths(
        &paths,
        var_data,
        reference,
        chr,
        start,
        end,
        config,
        chromosomes_lengths,
    )
}

fn partial_pipeline_paths(
    bam_paths: &[String],
    var_data: &mut crate::variations::VariationData,
    reference: &mut Reference,
    chr: &str,
    start: i32,
    end: i32,
    config: &Config,
    _chromosomes_lengths: &HashMap<String, usize>,
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
    let skip_owned_tile_range = var_data.split_parent_start.and(var_data.owning_range);
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
            let mut feeds = Vec::with_capacity(bam_paths.len());
            for single_bam in bam_paths.iter().filter(|s| !s.is_empty()) {
                feeds.push(crate::bam::query_region_iter(
                    std::path::Path::new(single_bam.as_str()),
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
                parser.reset_for_next_bam();
                let mut any_passed = false;
                while feed.next_into(&mut r)? {
                    if preprocessor.preprocess_record(&r, config) {
                        any_passed = true;
                        if let Some((own_s, own_e)) = skip_owned_tile_range {
                            if crate::record_overlaps_region(&r, own_s.0, own_e.0) {
                                continue;
                            }
                        }
                        parser.parse_cigar(&r)?;
                    }
                }
                if idx > 0 && !any_passed {
                    break;
                }
            }
            Ok(())
        })();

        parser.finish_ref_coverage();
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

    #[test]
    fn extreme_sv_min_reads_zero_neg_cnt_and_min_dir_guards() {
        let chr = "chr1";
        let mut reference = Reference::default();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 2_000,
            });
        reference
            .reference_sequences
            .insert(chr.to_string(), vec![b'A'; 2_000]);
        reference.region_sequences.push(vec![b'A'; 2_000]);

        // 1. find_match_rev with dir = i32::MIN uses saturating_neg without overflow.
        let seq20 = [b'C'; 20];
        let mr = find_match_rev(&seq20, &reference, 100, i32::MIN, 17, 1, chr, false);
        assert_eq!(mr.base_position, 0);

        // 2. mark_sv and mark_dup_sv saturate pairs on i32::MAX counts.
        let mk_sclip = |cnt: i32| {
            let mut sc = Sclip::default();
            sc.start = 100;
            sc.end = 200;
            sc.mstart = 500;
            sc.mend = 600;
            sc.variation.cnt = cnt;
            sc
        };
        let mut del_a = vec![mk_sclip(i32::MAX)];
        let mut del_b = vec![mk_sclip(i32::MAX)];
        let mut del_lists = [&mut del_a, &mut del_b];
        let (_cov, cnt_marked, pairs_marked) = mark_sv(200, 500, &mut del_lists, 100, false);
        assert_eq!(cnt_marked, 2);
        assert_eq!(pairs_marked, i32::MAX);

        let mut dup_a = vec![mk_sclip(i32::MAX)];
        let mut dup_b = vec![mk_sclip(i32::MAX)];
        let mut dup_lists = [&mut dup_a, &mut dup_b];
        let (dup_cnt, dup_pairs) = mark_dup_sv(100, 600, &mut dup_lists, 100, false);
        assert_eq!(dup_cnt, 2);
        assert_eq!(dup_pairs, i32::MAX);

        // 3. Discordant SV scanners (find_del_disc, find_inv_disc Part 1 & Part 2, find_dup_disc)
        //    must widen `min_reads + 5` to i64 (tested at 2_147_483_645 and u32::MAX) AND guard
        //    `cnt <= 0` so `cnt = -1` (with negative mqmean yielding positive ratio > 35.0) or
        //    `cnt = 0` can never reach `(cnt + 1) / 2 == 0` division-by-zero.
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert(chr.to_string(), 2_000usize);

        for min_reads in [0u32, 2_147_483_645u32, u32::MAX] {
            for bad_cnt in [-1i32, 0i32] {
                let mut cfg = Config::default();
                cfg.min_reads = min_reads;

                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 100;

                let mk_bad = |s: i32, e: i32, ms: i32, me: i32| {
                    let mut sc = Sclip::default();
                    sc.start = s;
                    sc.end = e;
                    sc.mstart = ms;
                    sc.mend = me;
                    sc.variation.cnt = bad_cnt;
                    sc.variation.mqmean = if bad_cnt < 0 { -100.0 } else { 0.0 };
                    sc
                };
                let mut partner = Sclip::default();
                partner.start = 150;
                partner.end = 250;
                partner.mstart = 1_150;
                partner.mend = 1_250;
                partner.variation.cnt = 10;
                partner.variation.mqmean = 500.0;

                vd.sv_structures.svfdel.push(mk_bad(100, 200, 1_100, 1_200));
                vd.sv_structures.svrdel.push(mk_bad(1_100, 1_200, 100, 200));

                vd.sv_structures
                    .svfinv5
                    .push(mk_bad(100, 200, 1_100, 1_200));
                vd.sv_structures.svrinv5.push(partner.clone());

                // Part 2: invf3 has cnt <= 0 while invr3 has rcnt = 10 so cnt + rcnt > 5 when min_reads == 0.
                vd.sv_structures
                    .svfinv3
                    .push(mk_bad(100, 200, 1_100, 1_200));
                vd.sv_structures.svrinv3.push(partner);

                vd.sv_structures.svfdup.push(mk_bad(100, 200, 1_100, 1_200));
                vd.sv_structures.svrdup.push(mk_bad(1_100, 1_200, 100, 200));

                find_del_disc(&mut vd, &mut reference, chr, &cfg, &chr_lengths).unwrap();
                find_inv_disc(&mut vd, &mut reference, chr, &cfg, &chr_lengths).unwrap();
                find_dup_disc(
                    &mut vd,
                    std::path::Path::new("/nonexistent/unused.bam"),
                    &mut reference,
                    chr,
                    &cfg,
                    &chr_lengths,
                )
                .unwrap();

                assert!(
                    vd.non_insertion_variants.is_empty() && vd.insertion_variants.is_empty(),
                    "cnt <= 0 or huge min_reads must not emit discordant SV candidates"
                );
            }
        }
    }

    #[test]
    fn find_all_sv_paths_preserves_colon_in_bam_filename_and_saturates_chr_len() {
        use noodles::bam;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles::csi::binning_index::Indexer;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::{
            Cigar, Data, QualityScores, RecordBuf, Sequence,
        };
        use std::num::NonZeroUsize;

        let dir = tempfile::TempDir::new().unwrap();
        let chr = "chr1";

        // 1. Write 200 bp FASTA + .fai so `find_dup_disc` can dynamically load mate reference.
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        let ref_seq = b"A".repeat(200);
        std::fs::write(
            &fasta_path,
            format!(">{chr}\n{}\n", String::from_utf8_lossy(&ref_seq)),
        )
        .unwrap();
        std::fs::write(&fai_path, format!("{chr}\t200\t6\t200\t201\n")).unwrap();

        // 2. Write indexed BAM + .bai whose filename contains literal `:` characters.
        let bam_path = dir.path().join("colon:sample:tumor.bam");
        let bai_path = dir.path().join("colon:sample:tumor.bam.bai");
        {
            let header = sam::Header::builder()
                .add_reference_sequence(
                    chr,
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(NonZeroUsize::new(200).unwrap()),
                )
                .build();
            let mut writer = bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            writer.write_header(&header).unwrap();

            let cigar: Cigar = [Op::new(Kind::Match, 20)].into_iter().collect();
            let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(1i32))]
                .into_iter()
                .collect();
            let mut seq = vec![b'A'; 20];
            seq[5] = b'C'; // 1-based pos 125 when alignment_start = 120
            let rec = RecordBuf::builder()
                .set_name("mate_colon_1")
                .set_flags(sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(120usize).unwrap())
                .set_mapping_quality(sam::alignment::record::MappingQuality::try_from(60).unwrap())
                .set_cigar(cigar)
                .set_sequence(Sequence::from(seq))
                .set_quality_scores(QualityScores::from(vec![35u8; 20]))
                .set_data(data)
                .build();

            let start_vp = writer.get_ref().virtual_position();
            writer.write_alignment_record(&header, &rec).unwrap();
            let end_vp = writer.get_ref().virtual_position();
            writer.try_finish().unwrap();

            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
            indexer
                .add_record(
                    Some((
                        0,
                        Position::try_from(120usize).unwrap(),
                        Position::try_from(139usize).unwrap(),
                        true,
                    )),
                    Chunk::new(start_vp, end_vp),
                )
                .unwrap();
            bam::bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }

        // Only load 1..=60 initially so `find_dup_disc` must run `partial_pipeline` on 120..=150.
        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: 60,
            });
        reference.region_sequences.push(vec![b'A'; 60]);
        reference
            .reference_sequences
            .insert(chr.to_string(), vec![b'A'; 60]);

        let mut var_data = crate::variations::VariationData::default();
        var_data.max_read_length = 50;
        let mut dup = Sclip::default();
        dup.start = 10;
        dup.end = 30;
        dup.mstart = 120;
        dup.mend = 150;
        dup.variation.cnt = 10;
        dup.variation.mqmean = 600.0;
        dup.variation.qmean = 350.0;
        dup.variation.pmean = 100.0;
        var_data.sv_structures.svfdup.push(dup);

        let mut config = Config::default();
        config.min_reads = 2;
        config.reference = fasta_path.to_string_lossy().into_owned();

        let mut chr_lengths = HashMap::new();
        chr_lengths.insert(chr.to_string(), 200usize);
        let pileup_region = crate::regions::Region {
            chr: chr.to_string(),
            start: Pos1(1),
            end: Pos1(60),
            gene: String::new(),
            amplicon: None,
        };

        let bam_paths = vec![bam_path.to_string_lossy().into_owned()];
        find_all_sv_paths(
            &mut var_data,
            &bam_paths,
            &mut reference,
            chr,
            200,
            &pileup_region,
            &config,
            &chr_lengths,
        )
        .expect("find_all_sv_paths must open BAM filenames containing ':' without splitting");

        assert!(
            var_data.non_insertion_variants.contains_key(&125),
            "mate-window partial_pipeline must have read colon-named BAM and planted SNV at 125"
        );

        // Also verify `saturated_chr_len` saturates >4 GiB `usize` lengths to `u32::MAX`
        // instead of wrapping via `as u32`.
        let mut huge_lens = HashMap::new();
        huge_lens.insert(chr.to_string(), (u32::MAX as usize).saturating_add(500));
        assert_eq!(saturated_chr_len(&huge_lens, chr), u32::MAX);
        assert_eq!(saturated_chr_len(&huge_lens, "missing_chr"), 0);

        // Ports StructuralVariantsProcessor.java:1508 — `partial_pipeline_paths` seeds
        // `CigarParser` with `var_data.max_read_length` and discards the returned scope,
        // so a longer 151 bp mate-window read must NOT write back to `var_data.max_read_length`.
        let bam_151 = dir.path().join("mate_151bp.bam");
        let bai_151 = dir.path().join("mate_151bp.bam.bai");
        {
            let header = sam::Header::builder()
                .add_reference_sequence(
                    chr,
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(NonZeroUsize::new(200).unwrap()),
                )
                .build();
            let mut writer = bam::io::Writer::new(std::fs::File::create(&bam_151).unwrap());
            writer.write_header(&header).unwrap();

            let cigar: Cigar = [Op::new(Kind::Match, 151)].into_iter().collect();
            let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(0i32))]
                .into_iter()
                .collect();
            let rec = RecordBuf::builder()
                .set_name("mate_151")
                .set_flags(sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(10usize).unwrap())
                .set_mapping_quality(sam::alignment::record::MappingQuality::try_from(60).unwrap())
                .set_cigar(cigar)
                .set_sequence(Sequence::from(vec![b'A'; 151]))
                .set_quality_scores(QualityScores::from(vec![35u8; 151]))
                .set_data(data)
                .build();
            let start_vp = writer.get_ref().virtual_position();
            writer.write_alignment_record(&header, &rec).unwrap();
            let end_vp = writer.get_ref().virtual_position();
            writer.try_finish().unwrap();

            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();
            indexer
                .add_record(
                    Some((
                        0,
                        Position::try_from(10usize).unwrap(),
                        Position::try_from(160usize).unwrap(),
                        true,
                    )),
                    Chunk::new(start_vp, end_vp),
                )
                .unwrap();
            bam::bai::io::Writer::new(std::fs::File::create(&bai_151).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }
        let mut ref_200 = Reference::new();
        ref_200.loaded_regions.push(crate::reference::LoadedRegion {
            chr: chr.to_string(),
            sequence_start: 1,
            sequence_end: 200,
        });
        ref_200.region_sequences.push(vec![b'A'; 200]);
        ref_200
            .reference_sequences
            .insert(chr.to_string(), vec![b'A'; 200]);
        let mut vd_mrl = crate::variations::VariationData::default();
        vd_mrl.max_read_length = 100;
        let bam_151_paths = vec![bam_151.to_string_lossy().into_owned()];
        partial_pipeline_paths(
            &bam_151_paths,
            &mut vd_mrl,
            &mut ref_200,
            chr,
            1,
            200,
            &config,
            &chr_lengths,
        )
        .unwrap();
        assert_eq!(vd_mrl.max_read_length, 100);
    }

    #[test]
    fn test_softp2sv_slot_map_live_used_mutation_and_max_rl_zero_ranges() {
        // 1. Direct build_softp2sv_map + is_softp_sv_used invariants (softp == 0 skipped, live .used observed)
        let mut vd = crate::variations::VariationData::default();
        let mut zero_sv = Sclip::new();
        zero_sv.softp = 0;
        zero_sv.variation.cnt = 99;
        zero_sv.used = true;
        vd.sv_structures.svfdel.push(zero_sv);

        let mut del_sv = Sclip::new();
        del_sv.softp = 250;
        del_sv.start = 95;
        del_sv.end = 105;
        del_sv.mstart = 295;
        del_sv.mend = 305;
        del_sv.variation.cnt = 8;
        del_sv.used = false;
        vd.sv_structures.svfdel.push(del_sv);

        let softp2sv = build_softp2sv_map(&vd.sv_structures);
        assert!(
            !is_softp_sv_used(&softp2sv, &vd.sv_structures, 0),
            "softp == 0 must never be indexed in SOFTP2SV"
        );
        assert!(!is_softp_sv_used(&softp2sv, &vd.sv_structures, 250));
        vd.sv_structures.svfdel[1].used = true;
        assert!(
            is_softp_sv_used(&softp2sv, &vd.sv_structures, 250),
            "SclipSlot must dereference live .used mutations without rebuilding map"
        );
        vd.sv_structures.svfdel[1].used = false;

        // 2. End-to-end find_sv mid-loop check_pairs .used = true visibility (AR-08):
        //    - p5_first = 300 (cnt = 10) matches forward at bp = 100 -> check_pairs(100, 300) marks svfdel[1] (softp=250) .used = true!
        //    - p5_second = 250 (cnt = 5, processed second) would otherwise match reverse-complement at bp = 400 and plant an INV at 250,
        //      UNLESS find_sv's SOFTP2SV lookup sees svfdel[1].used == true mutated during p5_first's iteration!
        let ref_len = 600usize;
        let motif = b"ACGTACGTAGCTAGCTTGCAACGTTAGCATGC";
        let mut ref_seq: Vec<u8> = (0..ref_len).map(|i| motif[i % motif.len()]).collect();
        let anchor_100 = b"TTCGAAGGTTCCAAGGTTCC"; // 20 bp unique forward target ending at 1-based 100 (indices 80..100)
        ref_seq[80..100].copy_from_slice(anchor_100);
        let inv_target_400 = b"GGGAAATTTCCCGGGAAATT"; // 20 bp at 1-based 401..=420 (indices 400..420)
        ref_seq[400..420].copy_from_slice(inv_target_400);

        let mut reference = Reference::new();
        reference
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: ref_len as u32,
            });
        reference.region_sequences.push(ref_seq.clone());
        reference
            .reference_sequences
            .insert("chr1".to_string(), ref_seq);
        reference.seed.insert_for_test(&anchor_100[3..20], vec![84]); // 17-mer at 1-based 84 so bp = 84 + 16 = 100!
        let mut rc_inv_17 = inv_target_400[0..17].to_vec();
        crate::util::reverse_complement(&mut rc_inv_17);
        reference
            .seed
            .insert_for_test(&inv_target_400[0..17], vec![401]);

        vd.max_read_length = 100;
        // Build 5' soft-clip at 300 (cnt = 10): consensus in 5' Sclip is stored 3'->5' from clip boundary,
        // and find_conseq_sets(sc5v, 0) returns the consensus string (reversed back by find_match dir=-1).
        // Pre-populating Sclip.sequence bypasses nt consensus building!
        let mut sc5_300 = Sclip::new();
        sc5_300.variation.cnt = 10;
        sc5_300.variation.hicnt = 10;
        sc5_300.sequence = Some(anchor_100.iter().rev().copied().collect());
        vd.soft_clips_5_end.insert(300, sc5_300);

        let mut sc5_250 = Sclip::new();
        sc5_250.variation.cnt = 5;
        sc5_250.variation.hicnt = 5;
        let mut inv_clip_seq = inv_target_400.to_vec();
        crate::util::reverse_complement(&mut inv_clip_seq);
        sc5_250.sequence = Some(inv_clip_seq.into_iter().rev().collect());
        vd.soft_clips_5_end.insert(250, sc5_250);

        let mut cfg = Config::default();
        cfg.min_reads = 2;
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), ref_len);

        find_sv(&mut vd, &mut reference, "chr1", 1, 500, &cfg, &chr_lengths).unwrap();
        assert!(
            vd.sv_structures.svfdel[1].used,
            "p5=300 check_pairs must have marked svfdel[1] used"
        );
        assert!(
            vd.non_insertion_variants.contains_key(&101),
            "p5=300 DEL must be called at bp_adj = 101"
        );
        assert!(
            vd.non_insertion_variants
                .iter()
                .all(|(_, vm)| vm.sv.as_ref().is_none_or(|sv| sv.r#type != "INV")),
            "p5=250 MUST be skipped because svfdel[1] (softp=250) became used mid-loop!"
        );

        // 3. find_all_sv_paths with max_read_length == 0 and -1 must not panic on BTreeMap range bounds (AR-07)
        for max_rl in [-1, 0, 1] {
            let mut vd_z = crate::variations::VariationData::default();
            vd_z.max_read_length = max_rl;
            vd_z.sv_structures_filtered = true;
            let mut fdel = Sclip::new();
            fdel.start = 90;
            fdel.end = 100;
            fdel.mstart = 200;
            fdel.mend = 210;
            fdel.variation.cnt = 5;
            vd_z.sv_structures.svfdel.push(fdel);
            let mut rdel = Sclip::new();
            rdel.start = 200;
            rdel.end = 210;
            rdel.mstart = 90;
            rdel.mend = 100;
            rdel.variation.cnt = 5;
            vd_z.sv_structures.svrdel.push(rdel);
            for pos in 95..=105 {
                let mut s3 = Sclip::new();
                s3.variation.cnt = 5;
                s3.sequence = Some(b"ACGT".to_vec()); // < 12 bp after window check
                vd_z.soft_clips_3_end.insert(pos, s3);
            }
            for pos in 195..=205 {
                let mut s5 = Sclip::new();
                s5.variation.cnt = 5;
                s5.sequence = Some(b"ACGT".to_vec());
                vd_z.soft_clips_5_end.insert(pos, s5);
            }
            let reg_z = crate::regions::Region {
                chr: "chr1".to_string(),
                start: Pos1(1),
                end: Pos1(500),
                gene: String::new(),
                amplicon: None,
            };
            find_all_sv_paths(
                &mut vd_z,
                &[],
                &mut reference,
                "chr1",
                ref_len as u32,
                &reg_z,
                &cfg,
                &chr_lengths,
            )
            .unwrap();
        }
    }

    #[test]
    fn ub6_mark_dup_sv_inclusive_end_marks_pre_decrement_boundary_only_when_split_decremented() {
        let chr = "chr1";
        let ref_len = 1_200usize;
        let motif = b"ACGTACGTAGCTAGCTTGCAACGTTAGCATGC";
        let mut ref_seq: Vec<u8> = (0..ref_len).map(|i| motif[i % motif.len()]).collect();
        // Place unique 20-bp anchor at 1-based 200..=219 (0-based 199..219), ensuring
        // 1-based 199 ('G') != 1-based 499 ('A') so the homology shift while-loop stops immediately.
        ref_seq[198] = b'G';
        ref_seq[498] = b'A';
        let anchor_200 = b"TTCGAAGGTTCCAAGGTTCC";
        ref_seq[199..219].copy_from_slice(anchor_200);

        let build_ref = || {
            let mut r = Reference::new();
            r.loaded_regions.push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: ref_len as u32,
            });
            r.region_sequences.push(ref_seq.clone());
            r.reference_sequences
                .insert(chr.to_string(), ref_seq.clone());
            // 17-mer at i=3 (anchor_200[3..20]) placed at 1-based 203 -> bp = 203 - 3 = 200.
            r.seed.insert_for_test(&anchor_200[3..20], vec![203]);
            r
        };

        let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);

        // Part 1: Split-refined svfdup (pe_orig = 500, tbp = 200 -> pe decremented to 499).
        // svrdup starts at 499 (start2 = 499, end2 = 600):
        // - COMPAT_VARDICTJAVA_1_8_3 passes pe = 499 -> 499 >= 499 -> NOT marked (clusters == 1).
        // - COMPAT_OFF passes mark_pe = 500 -> 499 < 500 -> marked (clusters == 2), while ins3 & end keep pe = 499.
        let run_split_case = |compat: crate::config::CompatProfile| {
            let mut reference = build_ref();
            let mut vd = crate::variations::VariationData::default();
            vd.max_read_length = 150;

            let mut fdup = Sclip::new();
            fdup.start = 450;
            fdup.end = 500;
            fdup.mstart = 200;
            fdup.mend = 260;
            fdup.softp = 500;
            fdup.soft.insert(500, 5);
            fdup.variation.cnt = 10;
            fdup.variation.mqmean = 600.0;
            fdup.variation.qmean = 350.0;
            fdup.variation.pmean = 200.0;
            vd.sv_structures.svfdup.push(fdup);

            let mut sc3 = Sclip::new();
            sc3.variation.cnt = 5;
            sc3.variation.mqmean = 300.0;
            sc3.variation.qmean = 175.0;
            sc3.variation.pmean = 100.0;
            sc3.sequence = Some(anchor_200.to_vec());
            vd.soft_clips_3_end.insert(500, sc3);

            let mut rdup = Sclip::new();
            rdup.mstart = 499;
            rdup.mend = 540;
            rdup.start = 560;
            rdup.end = 600;
            rdup.variation.cnt = 10;
            rdup.variation.mqmean = 600.0;
            rdup.variation.qmean = 350.0;
            rdup.variation.pmean = 200.0;
            vd.sv_structures.svrdup.push(rdup);

            let mut cfg = Config::default();
            cfg.min_reads = 2;
            cfg.sv_min_len = 100;
            cfg.compat = compat;

            find_dup_disc(
                &mut vd,
                std::path::Path::new(""),
                &mut reference,
                chr,
                &cfg,
                &chr_lengths,
            )
            .unwrap();
            vd
        };

        let vd_leg = run_split_case(crate::config::COMPAT_VARDICTJAVA_1_8_3);
        let sv_leg = vd_leg
            .non_insertion_variants
            .get(&200)
            .and_then(|vm| vm.sv.as_ref())
            .unwrap();
        assert_eq!(sv_leg.clusters, 1);

        let vd_ub6_only = run_split_case(crate::config::CompatProfile {
            mark_dup_sv_inclusive_end: true,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        });
        assert_eq!(
            vd_ub6_only
                .non_insertion_variants
                .get(&200)
                .and_then(|vm| vm.sv.as_ref())
                .unwrap()
                .clusters,
            2
        );

        let vd_off = run_split_case(crate::config::COMPAT_OFF);
        // Under full `COMPAT_OFF` (`canonical_sv_vcf_anchoring = true`), `<dupN>` shifts from
        // `bp = 200` to `anchor = bp - 1 = 199` while keeping identical `ins5`/`ins3` flanks.
        let sv_off = vd_off
            .non_insertion_variants
            .get(&199)
            .and_then(|vm| vm.sv.as_ref())
            .unwrap();
        assert_eq!(sv_off.clusters, 2);
        assert!(vd_off.sv_structures.svrdup[0].used);
        let keys_leg: Vec<&str> = vd_leg
            .insertion_variants
            .get(&200)
            .unwrap()
            .map
            .keys()
            .collect();
        let keys_off: Vec<&str> = vd_off
            .insertion_variants
            .get(&199)
            .unwrap()
            .map
            .keys()
            .collect();
        assert_eq!(keys_leg, keys_off);

        // Part 2: Pure discordant svfdup (empty soft map AND tbp == 0) must NOT add +1 under COMPAT_OFF!
        for with_unmatched_soft in [false, true] {
            let mut reference = build_ref();
            let mut vd = crate::variations::VariationData::default();
            vd.max_read_length = 150;

            let mut fdup = Sclip::new();
            fdup.start = 450;
            fdup.end = 500;
            fdup.mstart = 210;
            fdup.mend = 260;
            fdup.variation.cnt = 10;
            fdup.variation.mqmean = 600.0;
            fdup.variation.qmean = 350.0;
            fdup.variation.pmean = 200.0;
            if with_unmatched_soft {
                fdup.softp = 500;
                fdup.soft.insert(500, 4);
                let mut sc3 = Sclip::new();
                sc3.variation.cnt = 4;
                sc3.variation.mqmean = 240.0;
                sc3.variation.qmean = 140.0;
                sc3.variation.pmean = 80.0;
                sc3.sequence = Some(vec![b'T'; 20]); // tbp == 0 (unmatched) -> pe stays 500
                vd.soft_clips_3_end.insert(500, sc3);
            }
            vd.sv_structures.svfdup.push(fdup);

            // svrdup starting at 500 (start2 = 500, end2 = 600): would be falsely marked if pe + 1 (501) were passed!
            let mut rdup = Sclip::new();
            rdup.mstart = 500;
            rdup.mend = 540;
            rdup.start = 560;
            rdup.end = 600;
            rdup.variation.cnt = 4; // < min_reads + 5 so svrdup loop won't self-mark it
            rdup.variation.mqmean = 240.0;
            vd.sv_structures.svrdup.push(rdup);

            let mut cfg = Config::default();
            cfg.min_reads = 2;
            cfg.sv_min_len = 100;
            cfg.compat = crate::config::COMPAT_OFF;

            find_dup_disc(
                &mut vd,
                std::path::Path::new(""),
                &mut reference,
                chr,
                &cfg,
                &chr_lengths,
            )
            .unwrap();

            assert!(
                !vd.sv_structures.svrdup[0].used,
                "pure discordant / unmatched svfdup (with_unmatched_soft={with_unmatched_soft}) must keep mark_pe == pe == 500, not pe + 1"
            );
        }
    }

    #[test]
    fn sv_lazy_mate_pileup_avoids_double_pile_skips_unmatched_bam_and_preserves_shared_and_side3_mates(
    ) {
        use noodles::bam;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles::csi::binning_index::Indexer;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::{
            Cigar, Data, QualityScores, RecordBuf, Sequence,
        };
        use std::num::NonZeroUsize;

        let dir = tempfile::TempDir::new().unwrap();
        let chr = "chr1";
        let ref_len = 3_000usize;
        let motif = b"ACGTACGTAGCTAGCTTGCAACGTTAGCATGC";
        let mut ref_seq: Vec<u8> = (0..ref_len).map(|i| motif[i % motif.len()]).collect();
        // Distant Side::Three target ending at 1-based 200 (0-based 180..200)
        let inv3_target_200 = b"GGGAAATTTCCCGGGAAATT";
        ref_seq[180..200].copy_from_slice(inv3_target_200);
        // Ensure no complement homology shift at softp=201 / bp=1899 or softp=200 / bp=1900:
        // At 1-based 201 (idx 200) put 'A', at 1-based 1899 (idx 1898) put 'A' (complement('A')='T' != 'A');
        // At 1-based 200 (idx 199) is 'T', at 1-based 1900 (idx 1899) put 'T' (complement('T')='A' != 'T').
        ref_seq[200] = b'A';
        ref_seq[1898] = b'A';
        ref_seq[1899] = b'T';

        // Distant Side::Five shared-window target ending at 1-based 1580 (0-based 1560..1580)
        let inv5_target_1580 = b"CCCTTTAAAGGGCCCTTTAA";
        ref_seq[1560..1580].copy_from_slice(inv5_target_1580);
        ref_seq[1580] = b'C';
        ref_seq[199] = b'T'; // keep 200 intact

        let fasta_path = dir.path().join("lazy_sv.fa");
        let fai_path = dir.path().join("lazy_sv.fa.fai");
        std::fs::write(
            &fasta_path,
            format!(">{chr}\n{}\n", String::from_utf8_lossy(&ref_seq)),
        )
        .unwrap();
        std::fs::write(
            &fai_path,
            format!("{chr}\t{ref_len}\t6\t{ref_len}\t{}\n", ref_len + 1),
        )
        .unwrap();

        // Build indexed BAM with sorted reads at:
        // - 195..=214 (12 reads matching ref_seq, giving ref_coverage[200] = 12, ref_coverage[201] = 12, and planting 201 in non_insertion_variants)
        // - 310..=329 (4 reads inside primary region 1..=600 with SNV at 315 to detect double-piling)
        // - 1545..=1564 (6 reads with SNV at 1550 to detect shared-window partial_pipeline execution)
        let bam_path = dir.path().join("lazy_sv.bam");
        let bai_path = dir.path().join("lazy_sv.bam.bai");
        {
            let header = sam::Header::builder()
                .add_reference_sequence(
                    chr,
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(NonZeroUsize::new(ref_len).unwrap()),
                )
                .build();
            let mut writer = bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            writer.write_header(&header).unwrap();
            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();

            let mut write_rec = |name: &str, start_1based: usize, seq_bytes: Vec<u8>, nm: i32| {
                let len = seq_bytes.len();
                let cigar: Cigar = [Op::new(Kind::Match, len)].into_iter().collect();
                let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(nm))]
                    .into_iter()
                    .collect();
                let rec = RecordBuf::builder()
                    .set_name(name)
                    .set_flags(sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(start_1based).unwrap())
                    .set_mapping_quality(
                        sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                    )
                    .set_cigar(cigar)
                    .set_sequence(Sequence::from(seq_bytes))
                    .set_quality_scores(QualityScores::from(vec![35u8; len]))
                    .set_data(data)
                    .build();
                let start_vp = writer.get_ref().virtual_position();
                writer.write_alignment_record(&header, &rec).unwrap();
                let end_vp = writer.get_ref().virtual_position();
                indexer
                    .add_record(
                        Some((
                            0,
                            Position::try_from(start_1based).unwrap(),
                            Position::try_from(start_1based + len - 1).unwrap(),
                            true,
                        )),
                        Chunk::new(start_vp, end_vp),
                    )
                    .unwrap();
            };

            for k in 0..12 {
                write_rec(&format!("mate200_{k}"), 195, ref_seq[194..214].to_vec(), 0);
            }
            for k in 0..4 {
                let mut s = ref_seq[309..329].to_vec();
                s[5] = if s[5] == b'A' { b'G' } else { b'A' }; // 1-based 315
                write_rec(&format!("intra310_{k}"), 310, s, 1);
            }
            for k in 0..6 {
                let mut s = ref_seq[1544..1564].to_vec();
                s[5] = if s[5] == b'A' { b'G' } else { b'A' }; // 1-based 1550
                write_rec(&format!("mate1550_{k}"), 1545, s, 1);
            }
            writer.try_finish().unwrap();
            bam::bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }

        let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);

        // (a) Intra-tile SV [mstart=305, mend=340] inside already-loaded primary region 1..=600
        //     must NOT double-pile reads at 315 (VD stays 4, ref_coverage[315] stays 10).
        {
            let mut cfg = Config::default();
            cfg.min_reads = 2;
            cfg.sv_min_len = 100;
            cfg.reference = fasta_path.to_string_lossy().into_owned();
            cfg.compat = crate::config::COMPAT_OFF;

            let primary_reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1),
                end: Pos1(600),
                gene: String::new(),
                amplicon: None,
            };
            let mut reference = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.max_read_length = 100;
            vd.ref_coverage.insert(315, 10);
            let snv_base = if ref_seq[314] == b'A' { "G" } else { "A" };
            crate::variations::get_variation(&mut vd.non_insertion_variants, 315, snv_base).cnt = 4;

            let mut inv = Sclip::new();
            inv.start = 100;
            inv.end = 140;
            inv.mstart = 305;
            inv.mend = 340;
            inv.mlen = -220;
            inv.variation.cnt = 4;
            vd.sv_structures.svfinv5.push(inv);

            let mut fdup = Sclip::new();
            fdup.start = 100;
            fdup.end = 140;
            fdup.mstart = 305;
            fdup.mend = 340;
            fdup.variation.cnt = 10;
            fdup.variation.mqmean = 600.0;
            vd.sv_structures.svfdup.push(fdup);

            let bam_paths = vec![bam_path.to_string_lossy().into_owned()];
            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg,
                &chr_lengths,
            )
            .unwrap();

            assert_eq!(vd.ref_coverage.get(&315), Some(&10));
            assert_eq!(
                vd.non_insertion_variants
                    .get(&315)
                    .unwrap()
                    .map
                    .get(snv_base)
                    .unwrap()
                    .cnt,
                4
            );
        }

        // (b) Poly-G (used=true) / short (<12 bp) distant cluster runs 0 FASTA fetches and 0 BAM queries,
        //     and unmatched >=12 bp distant cluster (bp == 0) runs 0 BAM queries under COMPAT_OFF!
        {
            let mut reference = Reference::new();
            reference
                .loaded_regions
                .push(crate::reference::LoadedRegion {
                    chr: chr.to_string(),
                    sequence_start: 1,
                    sequence_end: 500,
                });
            reference.region_sequences.push(ref_seq[0..500].to_vec());

            let mut vd = crate::variations::VariationData::default();
            vd.max_read_length = 100;
            let mut inv_polyg = Sclip::new();
            inv_polyg.start = 100;
            inv_polyg.end = 150;
            inv_polyg.mstart = 2_000;
            inv_polyg.mend = 2_100;
            inv_polyg.mlen = -1_900;
            inv_polyg.soft.insert(150, 5);
            inv_polyg.variation.cnt = 5;
            vd.sv_structures.svfinv5.push(inv_polyg);

            // Poly-G clip already marked used=true by find_conseq low-complexity filter
            let mut sc_polyg = Sclip::new();
            sc_polyg.used = true;
            sc_polyg.sequence = Some(vec![b'G'; 30]);
            sc_polyg.variation.cnt = 5;
            vd.soft_clips_3_end.insert(150, sc_polyg);

            // Also short <12bp clip on second cluster with softp == 0
            let mut inv_short = Sclip::new();
            inv_short.start = 200;
            inv_short.end = 250;
            inv_short.mstart = 2_200;
            inv_short.mend = 2_300;
            inv_short.mlen = -2_000;
            inv_short.variation.cnt = 5;
            vd.sv_structures.svfinv5.push(inv_short);
            let mut sc_short = Sclip::new();
            sc_short.sequence = Some(b"ACGTACGTACG".to_vec()); // 11 bp < 12 bp
            sc_short.variation.cnt = 5;
            vd.soft_clips_3_end.insert(260, sc_short);

            let mut cfg_off = Config::default();
            cfg_off.min_reads = 2;
            cfg_off.reference = "/nonexistent/must_not_fetch_fasta.fa".to_string();
            cfg_off.compat = crate::config::COMPAT_OFF;

            // Both FASTA and BAM paths are nonexistent: succeeding proves 0 FASTA fetches and 0 BAM queries!
            find_inv_sub(
                0,
                1,
                Side::Five,
                &mut vd,
                std::path::Path::new("/nonexistent/must_not_query_bam.bam"),
                &mut reference,
                chr,
                &cfg_off,
                &chr_lengths,
            )
            .expect(
                "poly-G / <12bp distant clusters must perform 0 FASTA fetches and 0 BAM queries",
            );
            assert_eq!(reference.loaded_regions.len(), 1);

            // Unmatched >=12 bp clip loads FASTA once, hits bp == 0, and performs 0 BAM queries!
            let mut vd_unmatched = crate::variations::VariationData::default();
            vd_unmatched.max_read_length = 100;
            let mut inv_u = Sclip::new();
            inv_u.start = 100;
            inv_u.end = 150;
            inv_u.mstart = 1_500;
            inv_u.mend = 1_600;
            inv_u.mlen = -1_400;
            inv_u.soft.insert(150, 5);
            inv_u.variation.cnt = 5;
            vd_unmatched.sv_structures.svfinv5.push(inv_u);
            let mut sc_u = Sclip::new();
            sc_u.sequence = Some(vec![b'T'; 20]); // 20 bp >= 12, does not match FASTA -> bp == 0
            sc_u.variation.cnt = 5;
            vd_unmatched.soft_clips_3_end.insert(150, sc_u);

            cfg_off.reference = fasta_path.to_string_lossy().into_owned();
            find_inv_sub(
                0,
                1,
                Side::Five,
                &mut vd_unmatched,
                std::path::Path::new("/nonexistent/must_not_query_bam.bam"),
                &mut reference,
                chr,
                &cfg_off,
                &chr_lengths,
            )
            .expect("unmatched >=12bp distant cluster must perform 0 BAM queries");
            assert!(is_loaded(chr, 1_500, 1_600, &reference));
        }

        // (c) 2-cluster shared mate window [1500, 1600]:
        //     Cluster 0 loads FASTA + fails find_match_rev (bp == 0), then Cluster 1 (and find_dup_disc)
        //     MUST STILL run partial_pipeline_paths and see mate-window BAM reads at 1550!
        {
            let mut cfg_off = Config::default();
            cfg_off.min_reads = 2;
            cfg_off.reference = fasta_path.to_string_lossy().into_owned();
            cfg_off.compat = crate::config::COMPAT_OFF;

            let primary_reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1),
                end: Pos1(500),
                gene: String::new(),
                amplicon: None,
            };

            // Subcase c1: Cluster 0 (bp == 0) followed by Cluster 1 (matching INV at [1500, 1600])
            let mut ref_c1 = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut ref_c1).unwrap();
            let mut vd_c1 = crate::variations::VariationData::default();
            vd_c1.max_read_length = 100;

            let mut c0 = Sclip::new();
            c0.start = 100;
            c0.end = 150;
            c0.mstart = 1_500;
            c0.mend = 1_600;
            c0.mlen = -1_400;
            c0.soft.insert(150, 4);
            c0.variation.cnt = 4;
            vd_c1.sv_structures.svfinv5.push(c0.clone());
            let mut sc_c0 = Sclip::new();
            sc_c0.sequence = Some(vec![b'T'; 20]); // unmatched -> bp == 0 after loading FASTA
            sc_c0.variation.cnt = 4;
            vd_c1.soft_clips_3_end.insert(150, sc_c0.clone());

            let mut c1 = Sclip::new();
            c1.start = 160;
            c1.end = 200;
            c1.mstart = 1_500;
            c1.mend = 1_600;
            c1.mlen = -1_380;
            c1.soft.insert(200, 5);
            c1.variation.cnt = 5;
            vd_c1.sv_structures.svfinv5.push(c1);
            let mut rc_1580 = inv5_target_1580.to_vec();
            crate::util::reverse_complement(&mut rc_1580);
            let mut sc_c1 = Sclip::new();
            sc_c1.sequence = Some(rc_1580);
            sc_c1.variation.cnt = 5;
            vd_c1.soft_clips_3_end.insert(200, sc_c1);

            find_inv_sub(
                0,
                1,
                Side::Five,
                &mut vd_c1,
                &bam_path,
                &mut ref_c1,
                chr,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();
            assert!(
                vd_c1.non_insertion_variants.contains_key(&1_550),
                "Cluster 1 must still run partial_pipeline_paths after Cluster 0 loaded FASTA only"
            );

            // Subcase c2: Cluster 0 in find_inv_paths (bp == 0, FASTA-only) followed by find_dup_disc on [1500, 1600]
            let mut ref_c2 = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut ref_c2).unwrap();
            let mut vd_c2 = crate::variations::VariationData::default();
            vd_c2.sv_structures_filtered = true;
            vd_c2.max_read_length = 100;
            vd_c2.sv_structures.svfinv5.push(c0);
            vd_c2.soft_clips_3_end.insert(150, sc_c0);

            let mut fdup = Sclip::new();
            fdup.start = 100;
            fdup.end = 150;
            fdup.mstart = 1_500;
            fdup.mend = 1_600;
            fdup.variation.cnt = 10;
            fdup.variation.mqmean = 600.0;
            fdup.variation.qmean = 350.0;
            fdup.variation.pmean = 200.0;
            vd_c2.sv_structures.svfdup.push(fdup);

            let bam_paths = vec![bam_path.to_string_lossy().into_owned()];
            find_all_sv_paths(
                &mut vd_c2,
                &bam_paths,
                &mut ref_c2,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();
            assert!(
                vd_c2.non_insertion_variants.contains_key(&1_550),
                "find_dup_disc must still run partial_pipeline_paths when find_inv loaded FASTA only"
            );
        }

        // (d) Side::Three distant inversion (primary tile 1800..=2200, distant mate window 190..=220):
        //     swapped softp = 201 MUST be recorded in merge_planted AND read mate-window ref_coverage[200] == 12!
        {
            let primary_reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1_800),
                end: Pos1(2_200),
                gene: String::new(),
                amplicon: None,
            };
            let run_side3 = |compat: crate::config::CompatProfile| {
                let mut cfg = Config::default();
                cfg.min_reads = 2;
                cfg.reference = fasta_path.to_string_lossy().into_owned();
                cfg.compat = compat;

                let mut reference = Reference::new();
                get_reference(&primary_reg, 0, &chr_lengths, &cfg, &mut reference).unwrap();

                let mut vd = crate::variations::VariationData::default();
                vd.owning_range = Some((Pos1(1_800), Pos1(2_200)));
                vd.max_read_length = 100;
                vd.ref_coverage.insert(1_900, 30);

                let mut inv3 = Sclip::new();
                inv3.start = 1_850;
                inv3.end = 1_900;
                inv3.mstart = 190;
                inv3.mend = 220;
                inv3.mlen = -1_300;
                inv3.soft.insert(1_900, 5);
                inv3.variation.cnt = 3;
                vd.sv_structures.svfinv3.push(inv3);

                let mut rc_200 = inv3_target_200.to_vec();
                crate::util::reverse_complement(&mut rc_200);
                let mut sc3_1900 = Sclip::new();
                sc3_1900.sequence = Some(rc_200);
                sc3_1900.variation.cnt = 5;
                sc3_1900.variation.hicnt = 5;
                sc3_1900.variation.qmean = 175.0;
                sc3_1900.variation.mqmean = 300.0;
                sc3_1900.variation.pmean = 100.0;
                vd.soft_clips_3_end.insert(1_900, sc3_1900);

                find_inv_sub(
                    2,
                    1,
                    Side::Three,
                    &mut vd,
                    &bam_path,
                    &mut reference,
                    chr,
                    &cfg,
                    &chr_lengths,
                )
                .unwrap();
                vd
            };

            let vd_leg = run_side3(crate::config::COMPAT_VARDICTJAVA_1_8_3);
            let vd_off = run_side3(crate::config::COMPAT_OFF);

            for (label, vd) in [("legacy", &vd_leg), ("compat_off", &vd_off)] {
                let entry_201 = vd
                    .non_insertion_variants
                    .get(&201)
                    .unwrap_or_else(|| panic!("{label}: missing INV at swapped softp=201"));
                assert!(
                    entry_201.map.keys().any(|k| k.starts_with("-1699^")),
                    "{label}: expected -1699^ INV at 201"
                );
                assert!(
                    vd.merge_planted.contains(&201),
                    "{label}: swapped softp=201 in distant mate window MUST be recorded in merge_planted"
                );
                assert_eq!(
                    vd.ref_coverage.get(&201),
                    Some(&12),
                    "{label}: swapped softp=201 MUST read mate-window ref_coverage[200] == 12"
                );
            }
        }
    }

    #[test]
    fn find_del_fasta_only_load_does_not_starve_svrdel_side3_inv_or_dup_disc_mate_pileup() {
        use noodles::bam;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
        use noodles::csi::binning_index::Indexer;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::{
            Cigar, Data, QualityScores, RecordBuf, Sequence,
        };
        use std::num::NonZeroUsize;

        let dir = tempfile::TempDir::new().unwrap();
        let chr = "chr1";
        let ref_len = 4_500usize;
        let motif = b"ACGTACGTAGCTAGCTTGCAACGTTAGCATGC";
        let mut ref_seq: Vec<u8> = (0..ref_len).map(|i| motif[i % motif.len()]).collect();
        let inv3_target_200 = b"GGGAAATTTCCCGGGAAATT";
        ref_seq[180..200].copy_from_slice(inv3_target_200);
        ref_seq[200] = b'A';
        ref_seq[1898] = b'A';
        ref_seq[1899] = b'T';
        ref_seq[199] = b'T';
        ref_seq[2598] = b'A'; // 1-based 2599 != 'G' so 15S20M at 2600 stays anchored at 2600
        let chr_end_inv_target = b"TTTAAACCCGGGTTTAAACC";
        ref_seq[4470..4490].copy_from_slice(chr_end_inv_target); // 1-based 4471..=4490
        ref_seq[4490] = b'A';
        ref_seq[3901] = b'T';
        ref_seq[4489] = b'C';
        ref_seq[3900] = b'C';

        let fasta_path = dir.path().join("del_fasta_only.fa");
        let fai_path = dir.path().join("del_fasta_only.fa.fai");
        std::fs::write(
            &fasta_path,
            format!(">{chr}\n{}\n", String::from_utf8_lossy(&ref_seq)),
        )
        .unwrap();
        std::fs::write(
            &fai_path,
            format!("{chr}\t{ref_len}\t6\t{ref_len}\t{}\n", ref_len + 1),
        )
        .unwrap();

        let bam_path = dir.path().join("del_fasta_only.bam");
        let bai_path = dir.path().join("del_fasta_only.bam.bai");
        {
            let header = sam::Header::builder()
                .add_reference_sequence(
                    chr,
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(NonZeroUsize::new(ref_len).unwrap()),
                )
                .build();
            let mut writer = bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            writer.write_header(&header).unwrap();
            let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();

            let mut write_rec = |name: &str, start_1based: usize, seq_bytes: Vec<u8>, nm: i32| {
                let len = seq_bytes.len();
                let cigar: Cigar = [Op::new(Kind::Match, len)].into_iter().collect();
                let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(nm))]
                    .into_iter()
                    .collect();
                let rec = RecordBuf::builder()
                    .set_name(name)
                    .set_flags(sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(start_1based).unwrap())
                    .set_mapping_quality(
                        sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                    )
                    .set_cigar(cigar)
                    .set_sequence(Sequence::from(seq_bytes))
                    .set_quality_scores(QualityScores::from(vec![35u8; len]))
                    .set_data(data)
                    .build();
                let start_vp = writer.get_ref().virtual_position();
                writer.write_alignment_record(&header, &rec).unwrap();
                let end_vp = writer.get_ref().virtual_position();
                indexer
                    .add_record(
                        Some((
                            0,
                            Position::try_from(start_1based).unwrap(),
                            Position::try_from(start_1based + len - 1).unwrap(),
                            true,
                        )),
                        Chunk::new(start_vp, end_vp),
                    )
                    .unwrap();
            };

            for k in 0..12 {
                write_rec(&format!("mate200_{k}"), 195, ref_seq[194..214].to_vec(), 0);
            }
            for k in 0..7 {
                let mut s = ref_seq[1079..1099].to_vec();
                s[5] = if s[5] == b'A' { b'G' } else { b'A' }; // 1-based 1085
                write_rec(&format!("mate1085_{k}"), 1080, s, 1);
            }
            for k in 0..8 {
                let mut s = ref_seq[1544..1564].to_vec();
                s[5] = if s[5] == b'A' { b'G' } else { b'A' }; // 1-based 1550
                write_rec(&format!("mate1550_{k}"), 1545, s, 1);
            }
            for k in 0..6 {
                let mut s = ref_seq[2009..2029].to_vec();
                s[5] = if s[5] == b'A' { b'G' } else { b'A' }; // 1-based 2015
                write_rec(&format!("mate2015_{k}"), 2010, s, 1);
            }
            for k in 0..3 {
                let mut seq_bytes = b"TTCGAAGGTTCCAAG".to_vec();
                seq_bytes.extend_from_slice(&ref_seq[2599..2619]);
                let cigar: Cigar = [Op::new(Kind::SoftClip, 15), Op::new(Kind::Match, 20)]
                    .into_iter()
                    .collect();
                let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(0i32))]
                    .into_iter()
                    .collect();
                let rec = RecordBuf::builder()
                    .set_name(format!("mate2600_sc5_{k}"))
                    .set_flags(sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(2600).unwrap())
                    .set_mapping_quality(
                        sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                    )
                    .set_cigar(cigar)
                    .set_sequence(Sequence::from(seq_bytes))
                    .set_quality_scores(QualityScores::from(vec![35u8; 35]))
                    .set_data(data)
                    .build();
                let start_vp = writer.get_ref().virtual_position();
                writer.write_alignment_record(&header, &rec).unwrap();
                let end_vp = writer.get_ref().virtual_position();
                indexer
                    .add_record(
                        Some((
                            0,
                            Position::try_from(2600).unwrap(),
                            Position::try_from(2619).unwrap(),
                            true,
                        )),
                        Chunk::new(start_vp, end_vp),
                    )
                    .unwrap();
            }
            for k in 0..4 {
                let mut seq_bytes = ref_seq[2869..2889].to_vec(); // exon 1: 2870..=2889
                let mut exon2 = ref_seq[3409..3429].to_vec(); // exon 2 across 520N: 3410..=3429
                exon2[5] = if exon2[5] == b'A' { b'G' } else { b'A' }; // SNV at 1-based 3415
                seq_bytes.extend_from_slice(&exon2);
                let cigar: Cigar = [
                    Op::new(Kind::Match, 20),
                    Op::new(Kind::Skip, 520),
                    Op::new(Kind::Match, 20),
                ]
                .into_iter()
                .collect();
                let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(1i32))]
                    .into_iter()
                    .collect();
                let rec = RecordBuf::builder()
                    .set_name(format!("spliced_520n_{k}"))
                    .set_flags(sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(2870).unwrap())
                    .set_mapping_quality(
                        sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                    )
                    .set_cigar(cigar)
                    .set_sequence(Sequence::from(seq_bytes))
                    .set_quality_scores(QualityScores::from(vec![35u8; 40]))
                    .set_data(data)
                    .build();
                let start_vp = writer.get_ref().virtual_position();
                writer.write_alignment_record(&header, &rec).unwrap();
                let end_vp = writer.get_ref().virtual_position();
                indexer
                    .add_record(
                        Some((
                            0,
                            Position::try_from(2870).unwrap(),
                            Position::try_from(3429).unwrap(),
                            true,
                        )),
                        Chunk::new(start_vp, end_vp),
                    )
                    .unwrap();
            }
            writer.try_finish().unwrap();
            bam::bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                .write_index(&indexer.build(1))
                .unwrap();
        }

        let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);
        let primary_reg = crate::regions::Region {
            chr: chr.to_string(),
            start: Pos1(1_800),
            end: Pos1(2_200),
            gene: String::new(),
            amplicon: None,
        };
        let bam_paths = vec![bam_path.to_string_lossy().into_owned()];

        let mut cfg_off = Config::default();
        cfg_off.min_reads = 2;
        cfg_off.sv_min_len = 100;
        cfg_off.reference = fasta_path.to_string_lossy().into_owned();
        cfg_off.compat = crate::config::COMPAT_OFF;

        // 1. svfdel (softp_cand == 0, FASTA-only get_reference on [190, 220]) -> Side::Three find_inv_paths:
        //    must STILL run partial_pipeline_paths, populate ref_coverage[201] == 12, and record 201 in merge_planted!
        {
            let mut reference = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.owning_range = Some((Pos1(1_800), Pos1(2_200)));
            vd.max_read_length = 100;
            vd.ref_coverage.insert(1_900, 30);

            let mut fdel_no_soft = Sclip::new();
            fdel_no_soft.start = 1_820;
            fdel_no_soft.end = 1_860;
            fdel_no_soft.mstart = 190;
            fdel_no_soft.mend = 220;
            fdel_no_soft.variation.cnt = 4;
            vd.sv_structures.svfdel.push(fdel_no_soft);

            let mut inv3 = Sclip::new();
            inv3.start = 1_850;
            inv3.end = 1_900;
            inv3.mstart = 190;
            inv3.mend = 220;
            inv3.mlen = -1_300;
            inv3.soft.insert(1_900, 5);
            inv3.variation.cnt = 3;
            vd.sv_structures.svfinv3.push(inv3);

            let mut rc_200 = inv3_target_200.to_vec();
            crate::util::reverse_complement(&mut rc_200);
            let mut sc3_1900 = Sclip::new();
            sc3_1900.sequence = Some(rc_200);
            sc3_1900.variation.cnt = 5;
            sc3_1900.variation.hicnt = 5;
            sc3_1900.variation.qmean = 175.0;
            sc3_1900.variation.mqmean = 300.0;
            sc3_1900.variation.pmean = 100.0;
            vd.soft_clips_3_end.insert(1_900, sc3_1900);

            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            assert!(
                vd.merge_planted.contains(&201),
                "Side::Three swapped softp=201 must still be recorded in merge_planted after svfdel softp_cand==0 FASTA load"
            );
            assert_eq!(
                vd.ref_coverage.get(&201),
                Some(&12),
                "Side::Three swapped softp=201 must read mate-window ref_coverage[200] == 12 after svfdel softp_cand==0 FASTA load"
            );
            assert_eq!(
                reference.loaded_regions.len(),
                2,
                "svfdel softp_cand==0 + Side::Three INV must not append duplicate LoadedRegions"
            );
        }

        // 2. svfdel (softp_cand == 0, FASTA-only get_reference on [1500, 1600]) -> svrdel (softp_cand != 0):
        //    svrdel must STILL run partial_pipeline_paths and populate mate-window ref_coverage[1550] == 8!
        {
            let mut reference = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.max_read_length = 100;

            let mut fdel_no_soft = Sclip::new();
            fdel_no_soft.start = 1_820;
            fdel_no_soft.end = 1_860;
            fdel_no_soft.mstart = 1_500;
            fdel_no_soft.mend = 1_600;
            fdel_no_soft.variation.cnt = 4;
            vd.sv_structures.svfdel.push(fdel_no_soft);

            let mut rdel_with_soft = Sclip::new();
            rdel_with_soft.start = 1_950;
            rdel_with_soft.end = 1_990;
            rdel_with_soft.mstart = 1_500;
            rdel_with_soft.mend = 1_600;
            rdel_with_soft.soft.insert(1_950, 5);
            rdel_with_soft.variation.cnt = 4;
            vd.sv_structures.svrdel.push(rdel_with_soft);

            let mut sc5_1950 = Sclip::new();
            sc5_1950.sequence = Some(vec![b'T'; 20]); // >= 12 bp so svrdel enters mate-load block
            sc5_1950.variation.cnt = 5;
            vd.soft_clips_5_end.insert(1_950, sc5_1950);

            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            assert_eq!(
                vd.ref_coverage.get(&1_550),
                Some(&8),
                "svrdel softp_cand!=0 must still run partial_pipeline_paths after svfdel softp_cand==0 FASTA load"
            );
            assert!(vd.non_insertion_variants.contains_key(&1_550));
        }

        // 3. svfdel (softp_cand == 0, FASTA-only get_reference on [1500, 1600]) -> find_dup_disc_paths:
        //    find_dup_disc_paths must STILL run partial_pipeline_paths and populate ref_coverage[1550] == 8!
        {
            let mut reference = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.max_read_length = 100;

            let mut fdel_no_soft = Sclip::new();
            fdel_no_soft.start = 1_820;
            fdel_no_soft.end = 1_860;
            fdel_no_soft.mstart = 1_500;
            fdel_no_soft.mend = 1_600;
            fdel_no_soft.variation.cnt = 4;
            vd.sv_structures.svfdel.push(fdel_no_soft);

            let mut fdup = Sclip::new();
            fdup.start = 1_820;
            fdup.end = 1_860;
            fdup.mstart = 1_500;
            fdup.mend = 1_600;
            fdup.variation.cnt = 10;
            fdup.variation.mqmean = 600.0;
            fdup.variation.qmean = 350.0;
            fdup.variation.pmean = 200.0;
            vd.sv_structures.svfdup.push(fdup);

            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            assert_eq!(
                vd.ref_coverage.get(&1_550),
                Some(&8),
                "find_dup_disc_paths must still run partial_pipeline_paths after svfdel softp_cand==0 FASTA load"
            );
            assert!(vd.non_insertion_variants.contains_key(&1_550));
        }

        // 4. svfdel (softp_cand == 0, ext=300 FASTA load on [1520, 1580] -> LoadedRegion 1220..=1880)
        //    followed by Case-1 probe-overlapping svfinv5 (softp = 1820 inside widened [1220, 1880]):
        //    must NOT append a duplicate ext=500 LoadedRegion AND must still run partial_pipeline_paths!
        {
            let mut reference = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.max_read_length = 100;

            let mut fdel_no_soft = Sclip::new();
            fdel_no_soft.start = 1_820;
            fdel_no_soft.end = 1_860;
            fdel_no_soft.mstart = 1_520;
            fdel_no_soft.mend = 1_580;
            fdel_no_soft.variation.cnt = 4;
            vd.sv_structures.svfdel.push(fdel_no_soft);

            let mut inv_case1 = Sclip::new();
            inv_case1.start = 1_800;
            inv_case1.end = 1_820;
            inv_case1.mstart = 1_520;
            inv_case1.mend = 1_580;
            inv_case1.mlen = -300;
            inv_case1.soft.insert(1_820, 4);
            inv_case1.variation.cnt = 4;
            vd.sv_structures.svfinv5.push(inv_case1);

            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            assert_eq!(
                reference.loaded_regions.len(),
                2,
                "Case-1 find_inv_sub_paths must guard get_reference with !is_loaded after find_del ext=300 FASTA load"
            );
            assert_eq!(vd.ref_coverage.get(&1_550), Some(&8));
        }

        // 5. Reproducer (a): After svfdel (softp_cand == 0, +300 -> [1220, 1880]) + find_inv_paths on [1520, 1580],
        //    bam_piled_spans must record ONLY [1320, 1780] (NOT phantom [1020, 2080]!), so a downstream cluster on
        //    [1050, 1150] runs partial_pipeline_paths (ref_coverage[1085] == 7) AND find_dup_disc_paths emits a
        //    non-empty 50-bp reference flank around <dup...>!
        {
            let mut reference = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.max_read_length = 100;

            let mut fdel_no_soft = Sclip::new();
            fdel_no_soft.start = 1_820;
            fdel_no_soft.end = 1_860;
            fdel_no_soft.mstart = 1_520;
            fdel_no_soft.mend = 1_580;
            fdel_no_soft.variation.cnt = 4;
            vd.sv_structures.svfdel.push(fdel_no_soft);

            let mut inv_case1 = Sclip::new();
            inv_case1.start = 1_800;
            inv_case1.end = 1_820;
            inv_case1.mstart = 1_520;
            inv_case1.mend = 1_580;
            inv_case1.mlen = -300;
            inv_case1.soft.insert(1_820, 4);
            inv_case1.variation.cnt = 4;
            vd.sv_structures.svfinv5.push(inv_case1);

            let mut fdup_1050 = Sclip::new();
            fdup_1050.start = 1_820;
            fdup_1050.end = 1_860;
            fdup_1050.mstart = 1_050;
            fdup_1050.mend = 1_150;
            fdup_1050.variation.cnt = 10;
            fdup_1050.variation.mqmean = 600.0;
            fdup_1050.variation.qmean = 350.0;
            fdup_1050.variation.pmean = 200.0;
            vd.sv_structures.svfdup.push(fdup_1050);

            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            assert_eq!(
                vd.ref_coverage.get(&1_085),
                Some(&7),
                "downstream [1050, 1150] cluster must NOT be starved by phantom +500 bam_piled_spans"
            );
            // bp = 1050 - (100 / 10) / 2 = 1045 (`COMPAT_OFF` shifts `<dupN>` key to `bp - 1 = 1044`);
            // verify 50-bp ins5 flank before "<dup"
            let dup_entry = vd
                .insertion_variants
                .get(&1_044)
                .expect("expected DUP insertion at anchor = bp - 1 = 1044 under COMPAT_OFF");
            let dup_key = dup_entry.map.keys().next().unwrap();
            let prefix_before_dup = dup_key.split("<dup").next().unwrap();
            assert_eq!(
                prefix_before_dup.len(),
                1 + 50,
                "expected '+' + 50-bp loaded reference flank before <dup, got {dup_key}"
            );
        }

        // 6. Reproducer (b): Primary tile 100..=300; svfdel (softp_cand == 0 on [1500, 1600] -> FASTA [1200, 1900])
        //    followed by shifted svrdel (softp_cand != 0) on [1750, 1850] (mate_hi = 2050):
        //    must load FASTA over 1901..=2150 BEFORE partial_pipeline_paths so reads at 2010..=2029 pile up cleanly!
        {
            let low_primary_reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(100),
                end: Pos1(300),
                gene: String::new(),
                amplicon: None,
            };
            let mut reference = Reference::new();
            get_reference(&low_primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.max_read_length = 100;

            let mut fdel_no_soft = Sclip::new();
            fdel_no_soft.start = 150;
            fdel_no_soft.end = 200;
            fdel_no_soft.mstart = 1_500;
            fdel_no_soft.mend = 1_600;
            fdel_no_soft.variation.cnt = 4;
            vd.sv_structures.svfdel.push(fdel_no_soft);

            let mut rdel_shifted = Sclip::new();
            rdel_shifted.start = 220;
            rdel_shifted.end = 260;
            rdel_shifted.mstart = 1_750;
            rdel_shifted.mend = 1_850;
            rdel_shifted.soft.insert(220, 5);
            rdel_shifted.variation.cnt = 4;
            vd.sv_structures.svrdel.push(rdel_shifted);

            let mut sc5_220 = Sclip::new();
            sc5_220.sequence = Some(vec![b'T'; 20]);
            sc5_220.variation.cnt = 5;
            vd.soft_clips_5_end.insert(220, sc5_220);

            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &low_primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            assert!(
                is_loaded(chr, 1_901, 2_150, &reference),
                "shifted [1750, 1850] mate pileup must load FASTA across [mate_lo - max_rl, mate_hi + max_rl]"
            );
            assert_eq!(
                vd.ref_coverage.get(&2_015),
                Some(&6),
                "reads at 2010..=2029 in mate_hi window must see loaded FASTA bases during partial_pipeline_paths"
            );
        }

        // 7. `find_inv_disc` under `sv_lazy_mate_pileup` (both Part 2 `svfinv3+svrinv3` and Part 1 `svfinv5+svrinv5`):
        //    Part 2 must pile distant mate window [195, 215] BEFORE get_variation/ref_coverage so bp=205 gets
        //    real mate-window ref_coverage[205]==12 (not fabricated 2*cnt=10) AND bp=205 is in merge_planted;
        //    Part 1 must pile distant mate window [2600, 2620] so soft_clips_5_end[ms=2600] contributes splits=3!
        {
            let mut reference = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.owning_range = Some((Pos1(1_800), Pos1(2_200)));
            vd.max_read_length = 100;

            let mut invf3 = Sclip::new();
            invf3.start = 1_850;
            invf3.end = 1_900;
            invf3.mstart = 195;
            invf3.mend = 205;
            invf3.variation.cnt = 5;
            invf3.variation.mqmean = 300.0;
            invf3.variation.qmean = 175.0;
            invf3.variation.pmean = 100.0;
            vd.sv_structures.svfinv3.push(invf3);

            let mut invr3 = Sclip::new();
            invr3.start = 1_860;
            invr3.end = 1_910;
            invr3.mstart = 205;
            invr3.mend = 215;
            invr3.variation.cnt = 4;
            invr3.variation.mqmean = 240.0;
            invr3.variation.qmean = 140.0;
            invr3.variation.pmean = 80.0;
            vd.sv_structures.svrinv3.push(invr3);

            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            assert_eq!(
                vd.ref_coverage.get(&205),
                Some(&12),
                "find_inv_disc Part 2 must read real mate-window ref_coverage[205] == 12 instead of fabricating 2*cnt == 10"
            );
            assert!(
                vd.merge_planted.contains(&205),
                "find_inv_disc Part 2 distant mate breakpoint bp=205 must be recorded in merge_planted"
            );

            let mut reference_p1 = Reference::new();
            get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference_p1).unwrap();
            let mut vd_p1 = crate::variations::VariationData::default();
            vd_p1.sv_structures_filtered = true;
            vd_p1.max_read_length = 100;
            vd_p1.ref_coverage.insert(1_880, 25);

            let mut invf5 = Sclip::new();
            invf5.start = 1_850;
            invf5.end = 1_900;
            invf5.mstart = 2_600;
            invf5.mend = 2_620;
            invf5.variation.cnt = 5;
            invf5.variation.mqmean = 300.0;
            invf5.variation.qmean = 175.0;
            invf5.variation.pmean = 100.0;
            vd_p1.sv_structures.svfinv5.push(invf5);

            let mut invr5 = Sclip::new();
            invr5.start = 1_860;
            invr5.end = 1_910;
            invr5.mstart = 2_605;
            invr5.mend = 2_625;
            invr5.variation.cnt = 4;
            invr5.variation.mqmean = 240.0;
            invr5.variation.qmean = 140.0;
            invr5.variation.pmean = 80.0;
            vd_p1.sv_structures.svrinv5.push(invr5);

            find_all_sv_paths(
                &mut vd_p1,
                &bam_paths,
                &mut reference_p1,
                chr,
                ref_len as u32,
                &primary_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();

            let sv_p1 = vd_p1
                .non_insertion_variants
                .get(&1_880)
                .and_then(|e| e.sv.as_ref())
                .expect("expected Part 1 INV at bp=1880");
            assert_eq!(
                sv_p1.splits, 3,
                "find_inv_disc Part 1 must pile mate window [2600, 2620] so soft_clips_5_end[ms=2600] contributes splits=3"
            );
        }

        // 8. Spliced RNA-seq (`N`) / deletion (`D`) reference reach in `find_inv_sub_paths` `probe_overlaps`:
        //    Spliced read at 2870 (`20M520N20M`) starts inside `[mate_lo=2460, mate_hi=2880]` and spans a
        //    520-bp intron (`2889-3409`) to exon 2 at `3410..=3429` (past `mend + 500 = 3180` and outside
        //    unwidened `mate_hi + max_rl = 2940`, while probe from `end=3400` is `[3401, 3520]`)!
        {
            let spliced_reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(3_550),
                end: Pos1(3_750),
                gene: String::new(),
                amplicon: None,
            };
            let make_inv3_spliced = || {
                let mut inv = Sclip::new();
                inv.start = 3_360;
                inv.end = 3_400; // probe [3401, 3520] overlaps widened_hi = 2880 + 60 + 520 = 3460!
                inv.mstart = 2_660;
                inv.mend = 2_680;
                inv.mlen = -700;
                inv.variation.cnt = 4;
                inv
            };

            // Negative control without splice/deletion span: widened_hi = 2940 < 3401 -> no mate pileup
            let mut ref_plain = Reference::new();
            get_reference(&spliced_reg, 0, &chr_lengths, &cfg_off, &mut ref_plain).unwrap();
            let mut vd_plain = crate::variations::VariationData::default();
            vd_plain.sv_structures_filtered = true;
            vd_plain.max_read_length = 60;
            vd_plain.sv_structures.svfinv3.push(make_inv3_spliced());
            find_all_sv_paths(
                &mut vd_plain,
                &bam_paths,
                &mut ref_plain,
                chr,
                ref_len as u32,
                &spliced_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();
            assert_eq!(vd_plain.ref_coverage.get(&2_875), None);

            // Positive via `splice` ("2889-3409" -> 520N)
            let mut reference = Reference::new();
            get_reference(&spliced_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();
            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.max_read_length = 60;
            vd.splice.insert("2889-3409".to_string());
            vd.sv_structures.svfinv3.push(make_inv3_spliced());
            find_all_sv_paths(
                &mut vd,
                &bam_paths,
                &mut reference,
                chr,
                ref_len as u32,
                &spliced_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();
            assert!(
                is_loaded(chr, 3_181, 3_460, &reference),
                "spliced 520N reach must load FASTA across 3181..=3460 (past mend + 500 = 3180)"
            );
            assert_eq!(
                vd.ref_coverage.get(&2_875),
                Some(&4),
                "spliced 520N span must widen probe_overlaps and trigger Case-1 mate pileup on [2460, 2880]"
            );

            // Positive via `position_to_deletion_count` ("-520" -> 520D)
            let mut ref_del = Reference::new();
            get_reference(&spliced_reg, 0, &chr_lengths, &cfg_off, &mut ref_del).unwrap();
            let mut vd_del = crate::variations::VariationData::default();
            vd_del.sv_structures_filtered = true;
            vd_del.max_read_length = 60;
            vd_del
                .position_to_deletion_count
                .entry(3_560)
                .or_default()
                .insert("-520".to_string(), 2);
            vd_del.sv_structures.svfinv3.push(make_inv3_spliced());
            find_all_sv_paths(
                &mut vd_del,
                &bam_paths,
                &mut ref_del,
                chr,
                ref_len as u32,
                &spliced_reg,
                &cfg_off,
                &chr_lengths,
            )
            .unwrap();
            assert_eq!(
                vd_del.ref_coverage.get(&2_875),
                Some(&4),
                "long -520 deletion span must also widen probe_overlaps and trigger Case-1 mate pileup"
            );
        }

        // 9. `mend > chr_len` FASTA-only clamp in `find_inv_sub_paths` (`mend = 4520 > ref_len = 4500`):
        //    When full contig 1..=4500 is already loaded in `reference`, clamped `ensure_mate_fasta_loaded`
        //    must recognize `[4460, 4500]` as loaded without attempting an unclamped disk reload!
        {
            let full_reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1),
                end: Pos1(ref_len as u32),
                gene: String::new(),
                amplicon: None,
            };
            let mut reference = Reference::new();
            get_reference(&full_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();

            let mut cfg_inmem = cfg_off.clone();
            cfg_inmem.reference = dir
                .path()
                .join("must_not_reopen_fasta.fa")
                .to_string_lossy()
                .into_owned();

            let mut vd = crate::variations::VariationData::default();
            vd.sv_structures_filtered = true;
            vd.max_read_length = 100;
            vd.ref_coverage.insert(3_901, 20);

            let mut inv_end = Sclip::new();
            inv_end.start = 3_850;
            inv_end.end = 3_900;
            inv_end.mstart = 4_460;
            inv_end.mend = 4_520; // > ref_len (4500)
            inv_end.mlen = 500;
            inv_end.variation.cnt = 3;
            vd.sv_structures.svfinv5.push(inv_end);

            // First cp=3901 is non-matching >=12bp soft-clip; second cp=3902 matches 4490
            let mut sc_miss = Sclip::new();
            sc_miss.sequence = Some(vec![b'G'; 20]);
            sc_miss.variation.cnt = 2;
            vd.soft_clips_3_end.insert(3_901, sc_miss);

            let mut rc_end = chr_end_inv_target.to_vec();
            crate::util::reverse_complement(&mut rc_end);
            let mut sc_hit = Sclip::new();
            sc_hit.sequence = Some(rc_end);
            sc_hit.variation.cnt = 4;
            vd.soft_clips_3_end.insert(3_902, sc_hit);

            let mut bam_piled = vec![(1, ref_len as i32)];
            find_inv_sub_paths(
                0,
                1,
                Side::Five,
                &mut vd,
                &[],
                &mut reference,
                chr,
                &cfg_inmem,
                &chr_lengths,
                &mut bam_piled,
            )
            .expect(
                "mend > chr_len must clamp to chr_len and see preloaded contig without disk I/O",
            );
            assert!(
                vd.non_insertion_variants.contains_key(&3_902),
                "expected INV at softp=3902 when mend > chr_len clamps cleanly"
            );
        }
    }

    mod wave4_task_2d_tests {
        use super::*;
        use std::collections::BTreeMap;

        const TRUSEQ_ADAPTER_3P: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        const TRUSEQ_ADAPTER_5P: &[u8] = b"TCTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTG";

        fn fresh_sclip_with_nt(seq: &[u8], cnt: i32) -> Sclip {
            let mut sc = Sclip::new();
            sc.variation.cnt = cnt;
            sc.variation.hicnt = cnt;
            sc.variation.fwd = cnt / 2;
            sc.variation.rev = cnt - cnt / 2;
            sc.variation.qmean = 35.0 * (cnt as f64);
            sc.variation.mqmean = 60.0 * (cnt as f64);
            sc.variation.pmean = 25.0 * (cnt as f64);
            for (idx, &b) in seq.iter().enumerate() {
                sc.nt.insert(idx as i32, BTreeMap::from([(b, cnt)]));
            }
            assert!(sc.sequence.is_none());
            sc
        }

        fn build_ref(chr: &str, ref_len: usize) -> (Reference, Vec<u8>) {
            let motif = b"ACGTTGCAAGCTTAGCGATCGATCCGTAAGCT";
            let seq: Vec<u8> = (0..ref_len).map(|i| motif[i % motif.len()]).collect();
            let mut r = Reference::new();
            r.loaded_regions.push(crate::reference::LoadedRegion {
                chr: chr.to_string(),
                sequence_start: 1,
                sequence_end: ref_len as u32,
            });
            r.region_sequences.push(seq.clone());
            r.reference_sequences.insert(chr.to_string(), seq.clone());
            (r, seq)
        }

        /// Exercises all 10 flipped `sv.rs` `find_conseq*` call sites (`find_inv_sub_paths`
        /// `:1458` & `:1514` across both `dir == 1` and `dir == -1`, `find_del_paths` `:2031`,
        /// `:2223`, `:2357`, `:2513`, `find_sv` `:2707`, `:2909`, and `find_dup_disc_paths`
        /// `:3970`, `:4176`) under `adaptor_all_consensus_dirs` (`AstraZeneca-NGS/VarDictJava#240`,
        /// `AstraZeneca-NGS/VarDictJava#241`, `AstraZeneca-NGS/VarDict#119`).
        #[test]
        fn sv_all_10_consensus_sites_strip_3p_and_5p_truseq_adapters() {
            let chr = "chr1";
            let ref_len = 2_000usize;
            let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);
            let truseq_vec = vec![String::from_utf8_lossy(TRUSEQ_ADAPTER_3P).into_owned()];

            let run_all_10_sv_sites = |compat: crate::config::CompatProfile| {
                let (mut reference, _) = build_ref(chr, ref_len);
                let mut cfg = Config::default();
                cfg.min_reads = 2;
                cfg.adaptor = truseq_vec.clone();
                cfg.compat = compat;

                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 100;
                vd.sv_structures_filtered = true;
                let mut bam_piled = vec![(1, ref_len as i32)];

                // Sites 1 & 2 (`find_inv_sub_paths` `:1458` softp != 0 & `:1514` softp == 0 scan,
                // across both `dir = 1` -> `soft_clips_3_end` and `dir = -1` -> `soft_clips_5_end`):
                let mut inv_f5_soft = Sclip::new();
                inv_f5_soft.start = 100;
                inv_f5_soft.end = 150;
                inv_f5_soft.mstart = 400;
                inv_f5_soft.mend = 450;
                inv_f5_soft.variation.cnt = 4;
                inv_f5_soft.soft.insert(150, 4);
                vd.sv_structures.svfinv5.push(inv_f5_soft);
                vd.soft_clips_3_end
                    .insert(150, fresh_sclip_with_nt(TRUSEQ_ADAPTER_3P, 4));

                let mut inv_f5_scan = Sclip::new();
                inv_f5_scan.start = 200;
                inv_f5_scan.end = 250;
                inv_f5_scan.mstart = 500;
                inv_f5_scan.mend = 550;
                inv_f5_scan.variation.cnt = 4;
                vd.sv_structures.svfinv5.push(inv_f5_scan);
                vd.soft_clips_3_end
                    .insert(260, fresh_sclip_with_nt(TRUSEQ_ADAPTER_3P, 4));

                let mut inv_r5_soft = Sclip::new();
                inv_r5_soft.start = 300;
                inv_r5_soft.end = 350;
                inv_r5_soft.mstart = 600;
                inv_r5_soft.mend = 650;
                inv_r5_soft.variation.cnt = 4;
                inv_r5_soft.soft.insert(300, 4);
                vd.sv_structures.svrinv5.push(inv_r5_soft);
                vd.soft_clips_5_end
                    .insert(300, fresh_sclip_with_nt(TRUSEQ_ADAPTER_5P, 4));

                let mut inv_r5_scan = Sclip::new();
                inv_r5_scan.start = 380;
                inv_r5_scan.end = 420;
                inv_r5_scan.mstart = 700;
                inv_r5_scan.mend = 750;
                inv_r5_scan.variation.cnt = 4;
                vd.sv_structures.svrinv5.push(inv_r5_scan);
                vd.soft_clips_5_end
                    .insert(370, fresh_sclip_with_nt(TRUSEQ_ADAPTER_5P, 4));

                find_inv_sub_paths(
                    0,
                    1,
                    Side::Five,
                    &mut vd,
                    &[],
                    &mut reference,
                    chr,
                    &cfg,
                    &chr_lengths,
                    &mut bam_piled,
                )
                .unwrap();
                find_inv_sub_paths(
                    1,
                    -1,
                    Side::Five,
                    &mut vd,
                    &[],
                    &mut reference,
                    chr,
                    &cfg,
                    &chr_lengths,
                    &mut bam_piled,
                )
                .unwrap();

                // Sites 3, 4, 5, 6 (`find_del_paths`: `:2031`, `:2223`, `:2357`, `:2513`):
                let mut fdel_soft = Sclip::new();
                fdel_soft.start = 800;
                fdel_soft.end = 850;
                fdel_soft.mstart = 1_100;
                fdel_soft.mend = 1_150;
                fdel_soft.variation.cnt = 4;
                fdel_soft.soft.insert(850, 4);
                vd.sv_structures.svfdel.push(fdel_soft);
                vd.soft_clips_3_end
                    .insert(850, fresh_sclip_with_nt(TRUSEQ_ADAPTER_3P, 4));

                let mut fdel_scan = Sclip::new();
                fdel_scan.start = 900;
                fdel_scan.end = 950;
                fdel_scan.mstart = 1_200;
                fdel_scan.mend = 1_250;
                fdel_scan.variation.cnt = 4;
                vd.sv_structures.svfdel.push(fdel_scan);
                vd.soft_clips_3_end
                    .insert(960, fresh_sclip_with_nt(TRUSEQ_ADAPTER_3P, 4));

                let mut rdel_soft = Sclip::new();
                rdel_soft.start = 1_300;
                rdel_soft.end = 1_350;
                rdel_soft.mstart = 1_000;
                rdel_soft.mend = 1_050;
                rdel_soft.variation.cnt = 4;
                rdel_soft.soft.insert(1_300, 4);
                vd.sv_structures.svrdel.push(rdel_soft);
                vd.soft_clips_5_end
                    .insert(1_300, fresh_sclip_with_nt(TRUSEQ_ADAPTER_5P, 4));

                let mut rdel_scan = Sclip::new();
                rdel_scan.start = 1_400;
                rdel_scan.end = 1_450;
                rdel_scan.mstart = 1_060;
                rdel_scan.mend = 1_090;
                rdel_scan.variation.cnt = 4;
                vd.sv_structures.svrdel.push(rdel_scan);
                vd.soft_clips_5_end
                    .insert(1_390, fresh_sclip_with_nt(TRUSEQ_ADAPTER_5P, 4));

                find_del_paths(
                    &mut vd,
                    &[],
                    &mut reference,
                    chr,
                    &cfg,
                    &chr_lengths,
                    &mut bam_piled,
                )
                .unwrap();

                // Sites 7 & 8 (`find_sv`: `:2707` tmp5 & `:2909` tmp3):
                let mut vd_sv = crate::variations::VariationData::default();
                vd_sv.max_read_length = 100;
                vd_sv
                    .soft_clips_5_end
                    .insert(500, fresh_sclip_with_nt(TRUSEQ_ADAPTER_5P, 5));
                vd_sv
                    .soft_clips_3_end
                    .insert(550, fresh_sclip_with_nt(TRUSEQ_ADAPTER_3P, 5));
                find_sv(
                    &mut vd_sv,
                    &mut reference,
                    chr,
                    1,
                    ref_len as i32,
                    &cfg,
                    &chr_lengths,
                )
                .unwrap();

                // Sites 9 & 10 (`find_dup_disc_paths`: `:3970` svfdup & `:4176` svrdup):
                let mut vd_dup = crate::variations::VariationData::default();
                vd_dup.max_read_length = 100;
                let mut fdup = Sclip::new();
                fdup.start = 1_550;
                fdup.end = 1_600;
                fdup.mstart = 1_200;
                fdup.mend = 1_250;
                fdup.variation.cnt = 10;
                fdup.variation.mqmean = 600.0;
                fdup.variation.qmean = 350.0;
                fdup.variation.pmean = 200.0;
                fdup.softp = 1_600;
                fdup.soft.insert(1_600, 5);
                vd_dup.sv_structures.svfdup.push(fdup);
                vd_dup
                    .soft_clips_3_end
                    .insert(1_600, fresh_sclip_with_nt(TRUSEQ_ADAPTER_3P, 5));

                let mut rdup = Sclip::new();
                rdup.start = 1_650;
                rdup.end = 1_700;
                rdup.mstart = 1_850;
                rdup.mend = 1_900;
                rdup.variation.cnt = 10;
                rdup.variation.mqmean = 600.0;
                rdup.variation.qmean = 350.0;
                rdup.variation.pmean = 200.0;
                rdup.softp = 1_650;
                rdup.soft.insert(1_650, 5);
                vd_dup.sv_structures.svrdup.push(rdup);
                vd_dup
                    .soft_clips_5_end
                    .insert(1_650, fresh_sclip_with_nt(TRUSEQ_ADAPTER_5P, 5));

                find_dup_disc_paths(
                    &mut vd_dup,
                    &[],
                    &mut reference,
                    chr,
                    &cfg,
                    &chr_lengths,
                    &mut bam_piled,
                )
                .unwrap();

                (vd, vd_sv, vd_dup)
            };

            let (vd_leg, vd_sv_leg, vd_dup_leg) =
                run_all_10_sv_sites(crate::config::COMPAT_VARDICTJAVA_1_8_3);
            for pos3 in [150, 260, 850, 960] {
                assert_eq!(
                    vd_leg.soft_clips_3_end[&pos3].sequence.as_deref(),
                    Some(TRUSEQ_ADAPTER_3P),
                    "legacy pos3={pos3}"
                );
            }
            for pos5 in [300, 370, 1_300, 1_390] {
                assert_eq!(
                    vd_leg.soft_clips_5_end[&pos5].sequence.as_deref(),
                    Some(TRUSEQ_ADAPTER_5P),
                    "legacy pos5={pos5}"
                );
            }
            assert_eq!(
                vd_sv_leg.soft_clips_5_end[&500].sequence.as_deref(),
                Some(TRUSEQ_ADAPTER_5P)
            );
            assert_eq!(
                vd_sv_leg.soft_clips_3_end[&550].sequence.as_deref(),
                Some(TRUSEQ_ADAPTER_3P)
            );
            assert_eq!(
                vd_dup_leg.soft_clips_3_end[&1_600].sequence.as_deref(),
                Some(TRUSEQ_ADAPTER_3P)
            );
            assert_eq!(
                vd_dup_leg.soft_clips_5_end[&1_650].sequence.as_deref(),
                Some(TRUSEQ_ADAPTER_5P)
            );

            let (vd_off, vd_sv_off, vd_dup_off) = run_all_10_sv_sites(crate::config::COMPAT_OFF);
            for pos3 in [150, 260, 850, 960] {
                assert_eq!(
                    vd_off.soft_clips_3_end[&pos3].sequence.as_deref(),
                    Some(&b""[..]),
                    "COMPAT_OFF must strip 3' TruSeq adapter at pos3={pos3}"
                );
            }
            for pos5 in [300, 370, 1_300, 1_390] {
                assert_eq!(
                    vd_off.soft_clips_5_end[&pos5].sequence.as_deref(),
                    Some(&b""[..]),
                    "COMPAT_OFF must strip 5' TruSeq adapter at pos5={pos5}"
                );
            }
            assert_eq!(
                vd_sv_off.soft_clips_5_end[&500].sequence.as_deref(),
                Some(&b""[..]),
                "find_sv tmp5 must strip 5' TruSeq adapter"
            );
            assert_eq!(
                vd_sv_off.soft_clips_3_end[&550].sequence.as_deref(),
                Some(&b""[..]),
                "find_sv tmp3 must strip 3' TruSeq adapter"
            );
            assert_eq!(
                vd_dup_off.soft_clips_3_end[&1_600].sequence.as_deref(),
                Some(&b""[..]),
                "find_dup_disc svfdup must strip 3' TruSeq adapter"
            );
            assert_eq!(
                vd_dup_off.soft_clips_5_end[&1_650].sequence.as_deref(),
                Some(&b""[..]),
                "find_dup_disc svrdup must strip 5' TruSeq adapter"
            );
        }

        /// Verifies Tumor split-read `<DUP>` (`realignlgins` at `bp - 1`) and Normal discordant
        /// `<DUP>` (`sv.rs` `svfdup` and `svrdup` shifted to `bp - 1`) land at the identical map
        /// key `bp - 1` with `pre_anchored_dup == true` under `canonical_sv_vcf_anchoring`
        /// (`clintval/vartovcf#192`), while falling back to `bp` (`pre_anchored_dup == false`)
        /// under `COMPAT_VARDICTJAVA_1_8_3`, `bp == 1`, or unloaded/`'N'` `bp - 1`.
        #[test]
        fn dup_key_alignment_pairs_tumor_realignlgins_and_normal_sv_at_bp_minus_1_with_fallbacks() {
            let chr = "chr1";
            let ref_len = 1_200usize;
            let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);
            let region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1),
                end: Pos1(ref_len as u32),
                gene: String::new(),
                amplicon: None,
            };

            let (mut reference, mut ref_seq) = build_ref(chr, ref_len);
            ref_seq[298] = b'G'; // 1-based 299 (anchor = bp - 1)
            ref_seq[599] = b'T'; // 1-based 600 (pe)
            let anchor_581_600 = b"CGATCGTTAACCGGTTAACC"; // 1-based 581..=600
            ref_seq[580..600].copy_from_slice(anchor_581_600);
            reference.region_sequences[0].clone_from(&ref_seq);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq);
            reference
                .seed
                .insert_for_test(&anchor_581_600[3..20], vec![584]);

            let mut cfg_off = Config::default();
            cfg_off.bam = "tumor.bam|normal.bam".to_string();
            cfg_off.min_reads = 2;
            cfg_off.sv_min_len = 100;
            cfg_off.compat = crate::config::COMPAT_OFF;

            // 1. Tumor sample: split-read `<dup201>` on [300..=600] via `realignlgins` -> key 299
            let mut vd_tumor = crate::variations::VariationData::default();
            vd_tumor.max_read_length = 100;
            vd_tumor.ref_coverage.insert(299, 30);
            vd_tumor.ref_coverage.insert(300, 30);
            vd_tumor.ref_coverage.insert(600, 30);
            let mut rev_581_600 = anchor_581_600.to_vec();
            rev_581_600.reverse();
            vd_tumor
                .soft_clips_5_end
                .insert(300, fresh_sclip_with_nt(&rev_581_600, 8));
            let mut splice = rustc_hash::FxHashSet::default();
            crate::realign::realignlgins(
                None,
                &mut reference,
                chr,
                ref_len as u32,
                &region,
                &cfg_off,
                &mut vd_tumor.non_insertion_variants,
                &mut vd_tumor.insertion_variants,
                &mut vd_tumor.ref_coverage,
                &mut vd_tumor.soft_clips_5_end,
                &mut vd_tumor.soft_clips_3_end,
                vd_tumor.max_read_length,
                &mut vd_tumor.sv_structures.svfdup,
                &mut vd_tumor.sv_structures.svrdup,
                &mut splice,
            );

            // 2. Normal sample (`svfdup` AND `svrdup` discordant-only `<dup201>` on [300..=600]):
            //    `max_read_length = 100, cnt = 10 -> max_read_length / cnt = 10`:
            //    - `svfdup`: `mstart = 305, end = 591` -> `bp = 305 - 5 = 300`, `pe = 591`? Wait:
            //      for exact `[300..=600]` (`pe = 600, mlen = 301`): `end = 600` -> `pe = 600`,
            //      `mstart = 305` -> `bp = 305 - 5 = 300`, `mlen = 600 - 305 + 10 = 305`? Let's set
            //      `max_read_length = 16, cnt = 10 -> 16 / 10 = 1, (16 / 10) / 2 = 0`:
            //      - `svfdup`: `mstart = 300, end = 600` -> `bp = 300, pe = 600, mlen = 600 - 300 + 1 = 301`!
            //      - `svrdup`: `start = 300, mend = 600` -> `bp = 300, mlen = 600 - 300 + 1 = 301, pe = 301 + 300 - 1 = 600`!
            for use_svfdup in [true, false] {
                let mut vd_normal = crate::variations::VariationData::default();
                vd_normal.max_read_length = 16;
                let mut dup = Sclip::new();
                dup.variation.cnt = 10;
                dup.variation.mqmean = 600.0;
                dup.variation.qmean = 350.0;
                dup.variation.pmean = 200.0;
                if use_svfdup {
                    dup.start = 560;
                    dup.end = 600;
                    dup.mstart = 300;
                    dup.mend = 340;
                    vd_normal.sv_structures.svfdup.push(dup);
                } else {
                    dup.start = 300;
                    dup.end = 340;
                    dup.mstart = 560;
                    dup.mend = 600;
                    vd_normal.sv_structures.svrdup.push(dup);
                }

                find_dup_disc(
                    &mut vd_normal,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_off,
                    &chr_lengths,
                )
                .unwrap();

                assert!(
                    !vd_normal.insertion_variants.contains_key(&300),
                    "Normal discordant <DUP> must shift from bp=300 to bp-1=299 (use_svfdup={use_svfdup})"
                );
                let tumor_ins = vd_tumor.insertion_variants.get(&299).unwrap();
                let normal_ins = vd_normal
                    .insertion_variants
                    .get(&299)
                    .expect("Normal discordant <DUP> must land at anchor=299");
                let (t_desc, t_var) = tumor_ins.map.iter().next().unwrap();
                let (n_desc, n_var) = normal_ins.map.iter().next().unwrap();
                assert_eq!(
                    t_desc, n_desc,
                    "Tumor realignlgins and Normal sv.rs <DUP> must share exact allele description at 299 (use_svfdup={use_svfdup})"
                );
                assert!(t_var.pre_anchored_dup && n_var.pre_anchored_dup);
                assert!(
                    vd_normal
                        .non_insertion_variants
                        .get(&299)
                        .and_then(|vm| vm.sv.as_ref())
                        .is_some_and(|sv| sv.r#type == "DUP"),
                    "Normal SV metadata must also anchor at 299 so ToVarsBuilder visits 299"
                );
                assert!(vd_normal.ref_coverage.get(&299).copied().unwrap_or(0) > 0);

                // Verify `ToVarsBuilder` materializes BOTH Tumor and Normal `<DUP>` at key `299`!
                let builder =
                    crate::tovars::ToVarsBuilder::new(&reference, &region, &chr_lengths, &cfg_off);
                let aligned_t = builder.build_aligned_variants(&vd_tumor);
                let aligned_n = builder.build_aligned_variants(&vd_normal);
                assert_eq!(
                    aligned_t.get(&299).unwrap().variants[0].description,
                    aligned_n.get(&299).unwrap().variants[0].description
                );
            }

            // 3. Fallback guards:
            // (a) `COMPAT_VARDICTJAVA_1_8_3`: stays at `bp = 300` with `pre_anchored_dup == false`.
            {
                let mut cfg_leg = cfg_off.clone();
                cfg_leg.compat = crate::config::COMPAT_VARDICTJAVA_1_8_3;
                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 16;
                let mut dup = Sclip::new();
                dup.start = 560;
                dup.end = 600;
                dup.mstart = 300;
                dup.mend = 340;
                dup.variation.cnt = 10;
                dup.variation.mqmean = 600.0;
                dup.variation.qmean = 350.0;
                dup.variation.pmean = 200.0;
                vd.sv_structures.svfdup.push(dup);
                find_dup_disc(
                    &mut vd,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_leg,
                    &chr_lengths,
                )
                .unwrap();
                let (_, v) = vd
                    .insertion_variants
                    .get(&300)
                    .unwrap()
                    .map
                    .iter()
                    .next()
                    .unwrap();
                assert!(!v.pre_anchored_dup);
                assert!(!vd.insertion_variants.contains_key(&299));
            }

            // (b) Contig start `bp == 1`: stays at `bp = 1` with `pre_anchored_dup == false`.
            {
                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 16;
                let mut dup = Sclip::new();
                dup.start = 260;
                dup.end = 300;
                dup.mstart = 1;
                dup.mend = 40;
                dup.variation.cnt = 10;
                dup.variation.mqmean = 600.0;
                dup.variation.qmean = 350.0;
                dup.variation.pmean = 200.0;
                vd.sv_structures.svfdup.push(dup);
                find_dup_disc(
                    &mut vd,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_off,
                    &chr_lengths,
                )
                .unwrap();
                let (_, v) = vd
                    .insertion_variants
                    .get(&1)
                    .unwrap()
                    .map
                    .iter()
                    .next()
                    .unwrap();
                assert!(!v.pre_anchored_dup);
                assert!(!vd.insertion_variants.contains_key(&0));
            }

            // (c) `'N'` / `'n'` or unloaded reference base at `bp - 1`: stays at `bp` with `pre_anchored_dup == false`.
            for n_byte in [Some(b'N'), Some(b'n'), None] {
                let (mut ref_n, mut seq_n) = build_ref(chr, ref_len);
                if let Some(nb) = n_byte {
                    seq_n[298] = nb; // 1-based 299 (`bp - 1`) is 'N' or 'n'
                    ref_n.region_sequences[0].clone_from(&seq_n);
                    ref_n.reference_sequences.insert(chr.to_string(), seq_n);
                } else {
                    // Unloaded at 299: restrict loaded span to 300..=1200!
                    ref_n.loaded_regions[0].sequence_start = 300;
                    ref_n.region_sequences[0] = seq_n[299..].to_vec();
                    ref_n.reference_sequences.clear();
                }

                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 16;
                let mut dup = Sclip::new();
                dup.start = 560;
                dup.end = 600;
                dup.mstart = 300;
                dup.mend = 340;
                dup.variation.cnt = 10;
                dup.variation.mqmean = 600.0;
                dup.variation.qmean = 350.0;
                dup.variation.pmean = 200.0;
                vd.sv_structures.svfdup.push(dup);
                let mut bam_piled = vec![(1, ref_len as i32)];
                find_dup_disc_paths(
                    &mut vd,
                    &[],
                    &mut ref_n,
                    chr,
                    &cfg_off,
                    &chr_lengths,
                    &mut bam_piled,
                )
                .unwrap();
                let (_, v) = vd
                    .insertion_variants
                    .get(&300)
                    .unwrap_or_else(|| panic!("expected fallback at bp=300 for {n_byte:?}"))
                    .map
                    .iter()
                    .next()
                    .unwrap();
                assert!(!v.pre_anchored_dup);
                assert!(!vd.insertion_variants.contains_key(&299));
            }
        }

        #[test]
        fn wave5_3c_switch24_sv_resets_clear_aux_while_preserving_explicit_dup_accumulation() {
            // AstraZeneca-NGS/VarDictJava#316: sv.rs cnt=0 reset sites must clear all 9 auxiliary
            // Variation accumulators when reconcile_strand_and_hq_sums is active, EXCEPT when
            // find_dup_disc emits an explicit +<seq> DUP (`emit_explicit_dup && vref.cnt > 0`),
            // where pre-existing CIGAR/realignlgins +<seq> counts must accumulate cleanly!
            let chr = "chr1";
            let ref_len = 3_000usize;
            let (mut reference, mut ref_seq) = build_ref(chr, ref_len);
            let anchor_200_219 = b"CGATCGTTAACCGGTTAACC"; // 1-based 200..=219
            ref_seq[199..219].copy_from_slice(anchor_200_219);
            reference.region_sequences[0].clone_from(&ref_seq);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq.clone());
            reference
                .seed
                .insert_for_test(&anchor_200_219[0..17], vec![200]);
            let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);

            // 1. Discordant DEL reset site (`find_del_disc` svfdel): pre-seed stale Variation at softp=455 for "-1500"
            for reconcile in [false, true] {
                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 100; // mindist = 800
                let mut del = Sclip::new();
                del.start = 400;
                del.end = 450;
                del.mstart = 1_959;
                del.mend = 2_010;
                del.softp = 455;
                del.variation.cnt = 10; // mlen = 1959 - 450 - (100 / 11) = 1500 > 1000
                del.variation.mqmean = 600.0;
                del.variation.qmean = 350.0;
                del.variation.pmean = 200.0;
                del.variation.nm = 5.0;
                vd.sv_structures.svfdel.push(del);

                let stale =
                    crate::variations::get_variation(&mut vd.non_insertion_variants, 455, "-1500");
                stale.cnt = 3;
                stale.fwd = 9;
                stale.rev = 8;
                stale.hicnt = 11;
                stale.locnt = 4;
                stale.extracnt = 6;
                stale.qmean = 999.0;

                let cfg = Config {
                    min_reads: 2,
                    compat: crate::config::CompatProfile {
                        reconcile_strand_and_hq_sums: reconcile,
                        ..crate::config::COMPAT_VARDICTJAVA_1_8_3
                    },
                    ..Config::default()
                };
                find_del_disc(&mut vd, &mut reference, chr, &cfg, &chr_lengths).unwrap();
                let v = vd
                    .non_insertion_variants
                    .get(&455)
                    .unwrap()
                    .map
                    .get("-1500")
                    .unwrap();
                assert_eq!(v.cnt, 20);
                if reconcile {
                    assert_eq!(
                        (v.fwd, v.rev, v.hicnt, v.locnt, v.extracnt),
                        (10, 10, 20, 0, 20)
                    );
                    assert_eq!(v.qmean, 700.0);
                } else {
                    assert_eq!(
                        (v.fwd, v.rev, v.hicnt, v.locnt, v.extracnt),
                        (19, 18, 31, 4, 26)
                    );
                }
            }

            // 2. Explicit 285-bp DUP (`[200..=484]`, anchor = 199) in `find_dup_disc` with pre-existing
            //    CIGAR/realignlgins `+<285bp>` count = 6 MUST accumulate (`6 + (5 + 10) = 21`) instead of resetting!
            let dup_285 = String::from_utf8_lossy(&ref_seq[199..484]).into_owned();
            assert_eq!(dup_285.len(), 285);
            let plus_285 = format!("+{dup_285}");
            let mut vd_acc = crate::variations::VariationData::default();
            vd_acc.max_read_length = 100;
            let pre =
                crate::variations::get_variation(&mut vd_acc.insertion_variants, 199, &plus_285);
            pre.cnt = 6;
            pre.hicnt = 6;
            pre.fwd = 3;
            pre.rev = 3;
            pre.extracnt = 6;
            pre.qmean = 210.0;
            pre.mqmean = 360.0;
            pre.pmean = 120.0;

            let mut fdup = Sclip::new();
            fdup.start = 435;
            fdup.end = 485;
            fdup.mstart = 205;
            fdup.mend = 255;
            fdup.softp = 485;
            fdup.soft.insert(485, 5);
            fdup.variation.cnt = 10;
            fdup.variation.mqmean = 600.0;
            fdup.variation.qmean = 350.0;
            fdup.variation.pmean = 200.0;
            vd_acc.sv_structures.svfdup.push(fdup);

            let mut sc3 = Sclip::new();
            sc3.variation.cnt = 5;
            sc3.variation.mqmean = 300.0;
            sc3.variation.qmean = 175.0;
            sc3.variation.pmean = 100.0;
            sc3.sequence = Some(anchor_200_219.to_vec());
            vd_acc.soft_clips_3_end.insert(485, sc3);

            // Ensure non-homology stop at 199 vs 484: ref[198] != ref[483]
            let b_484 = ref_seq[483];
            ref_seq[198] = if b_484 == b'A' { b'G' } else { b'A' };
            reference.region_sequences[0].clone_from(&ref_seq);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq);

            let cfg_both = Config {
                min_reads: 2,
                sv_min_len: 1000,
                compat: crate::config::COMPAT_OFF,
                ..Config::default()
            };
            find_dup_disc(
                &mut vd_acc,
                std::path::Path::new(""),
                &mut reference,
                chr,
                &cfg_both,
                &chr_lengths,
            )
            .unwrap();
            let merged = vd_acc
                .insertion_variants
                .get(&199)
                .expect("expected explicit +<285bp> at anchor=199")
                .map
                .get(&plus_285)
                .expect("expected exact +<285bp> key");
            assert_eq!(
                merged.cnt, 21,
                "pre-existing realignlgins/CIGAR +<seq> count (6) must accumulate with split+pair count (15)"
            );
            assert_eq!((merged.fwd, merged.rev, merged.hicnt), (8, 13, 21));
        }

        #[test]
        fn wave5_3c_switch27_sv_respect_min_len_flt3_itd_285bp_split_inv_ceiling_and_bp1() {
            // AstraZeneca-NGS/VarDictJava#328, AstraZeneca-NGS/VarDict#157:
            // Comprehensive Switch 27 test covering:
            // (a) 285 bp FLT3-ITD split-refined svfdup & svrdup (< sv_min_len=1000) -> explicit +<285bp> at bp - 1,
            //     while discordant-only < 1000 bp DUP, DEL, and INV are suppressed!
            // (b) Split-read exact find_inv_sub & find_sv (tmp5 & tmp3) (< sv_min_len=1000) -> explicit complement_seq
            //     when loaded, falling back safely to symbolic <invN> when unloaded!
            // (c) -L 1000000 ceiling (MAX_EXPLICIT_SV_SEQ_BP = 32_768): 40,000 bp DUP & INV stay symbolic!
            // (d) Contig start bp == 1 tandem duplication [1..=pe]: right-rotated 1 bp onto bp = 1
            //     (+<dup_seq[1..]><ref_base_1>) so VCF at POS=1 has REF=b1 and mutant haplotype dup_seq + dup_seq!
            let chr = "chr1";
            let ref_len = 45_000usize;
            let (mut reference, mut ref_seq) = build_ref(chr, ref_len);

            // Unique 20-bp anchors for 285-bp FLT3-ITD [300..=584], bp==1 DUP [1..=200], and 250-bp INV [700..=949]
            ref_seq[0] = b'G'; // pos 1 = 'G'
            ref_seq[199] = b'T'; // pos 200 = 'T', ref[200] (pos 201) != ref[0] ('G') to stop homology walk
            ref_seq[200] = b'A';
            let bp1_anchor_1_20 = b"GCTTAGCTTAGCTTAGCTTA"; // 1-based 1..=20 (starts with 'G')
            ref_seq[0..20].copy_from_slice(bp1_anchor_1_20);

            ref_seq[298] = b'C'; // pos 299 = 'C' != pos 584 ('T')
            let flt3_300_319 = b"ATCGATCGTTAACCGGTTAA"; // 1-based 300..=319
            ref_seq[299..319].copy_from_slice(flt3_300_319);
            let flt3_565_584 = b"TTAACCGGTTAACCGGACGT"; // 1-based 565..=584 (ends with 'T' != 'C' at 299)
            ref_seq[564..584].copy_from_slice(flt3_565_584);
            ref_seq[584] = b'G'; // pos 585 = 'G' != pos 300 ('A')

            // 250-bp INV [700..=949]: ensure no complement homology beyond [700, 949]
            // pos 699 = 'A', pos 950 = 'A' (complement('A') = 'T' != 'A')
            ref_seq[698] = b'A';
            ref_seq[949] = b'A';
            let inv_930_949 = b"GGCCAATTGGCCAATTACGT"; // 1-based 930..=949 ('T' at 949; pos 700 = 'T' -> comp('T')='A'!='T')
            ref_seq[929..949].copy_from_slice(inv_930_949);
            ref_seq[699] = b'T'; // pos 700 = 'T'

            // 40,000-bp ceiling anchor at [1000..=40999]
            ref_seq[998] = b'C';
            let ceil_1000_1019 = b"AACCTTGGCCAAGGTTCCAA"; // 1-based 1000..=1019
            ref_seq[999..1019].copy_from_slice(ceil_1000_1019);
            ref_seq[40_998] = b'T'; // pos 40999 = 'T' != 'C' at 999

            reference.region_sequences[0].clone_from(&ref_seq);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq.clone());
            reference
                .seed
                .insert_for_test(&bp1_anchor_1_20[0..17], vec![1]);
            reference
                .seed
                .insert_for_test(&flt3_300_319[0..17], vec![300]);
            reference
                .seed
                .insert_for_test(&flt3_565_584[3..20], vec![568]);
            reference
                .seed
                .insert_for_test(&inv_930_949[3..20], vec![933]);
            reference
                .seed
                .insert_for_test(&ceil_1000_1019[0..17], vec![1000]);

            let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);
            let region = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1),
                end: Pos1(ref_len as u32),
                gene: "FLT3".to_string(),
                amplicon: None,
            };

            let sw27 = crate::config::CompatProfile {
                sv_respect_min_len: true,
                ..crate::config::COMPAT_VARDICTJAVA_1_8_3
            };
            let cfg_1000 = Config {
                min_reads: 2,
                sv_min_len: 1000,
                compat: sw27,
                ..Config::default()
            };

            // (a1) 285-bp FLT3-ITD split-refined via BOTH svfdup and svrdup -> explicit +<285bp> at 299!
            let expected_flt3_seq = String::from_utf8_lossy(&ref_seq[299..584]).into_owned();
            assert_eq!(expected_flt3_seq.len(), 285);
            for use_svfdup in [true, false] {
                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 100;
                vd.ref_coverage.insert(299, 40);
                let mut dup = Sclip::new();
                dup.variation.cnt = 10;
                dup.variation.mqmean = 600.0;
                dup.variation.qmean = 350.0;
                dup.variation.pmean = 200.0;
                if use_svfdup {
                    dup.start = 535;
                    dup.end = 585;
                    dup.mstart = 305;
                    dup.mend = 355;
                    dup.softp = 585;
                    dup.soft.insert(585, 6);
                    vd.sv_structures.svfdup.push(dup);
                    let mut sc3 = Sclip::new();
                    sc3.variation.cnt = 6;
                    sc3.variation.mqmean = 360.0;
                    sc3.variation.qmean = 210.0;
                    sc3.variation.pmean = 120.0;
                    sc3.sequence = Some(flt3_300_319.to_vec());
                    vd.soft_clips_3_end.insert(585, sc3);
                } else {
                    dup.start = 300;
                    dup.end = 350;
                    dup.mstart = 535;
                    dup.mend = 584;
                    dup.softp = 300;
                    dup.soft.insert(300, 6);
                    vd.sv_structures.svrdup.push(dup);
                    let mut rev_565_584 = flt3_565_584.to_vec();
                    rev_565_584.reverse();
                    let mut sc5 = Sclip::new();
                    sc5.variation.cnt = 6;
                    sc5.variation.mqmean = 360.0;
                    sc5.variation.qmean = 210.0;
                    sc5.variation.pmean = 120.0;
                    sc5.sequence = Some(rev_565_584);
                    vd.soft_clips_5_end.insert(300, sc5);
                }
                find_dup_disc(
                    &mut vd,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_1000,
                    &chr_lengths,
                )
                .unwrap();
                let ins_299 = vd
                    .insertion_variants
                    .get(&299)
                    .unwrap_or_else(|| panic!("missing 299 for use_svfdup={use_svfdup}"));
                let (desc, _) = ins_299.map.iter().next().unwrap();
                assert_eq!(desc, &format!("+{expected_flt3_seq}"));
                assert!(!desc.contains("<dup"));

                let aligned =
                    crate::tovars::ToVarsBuilder::new(&reference, &region, &chr_lengths, &cfg_1000)
                        .build_aligned_variants(&vd);
                let v = &aligned.get(&299).unwrap().variants[0];
                assert_eq!(v.start_pos, 299);
                assert_eq!(v.ref_allele, "C");
                assert_eq!(v.var_allele, format!("C{expected_flt3_seq}"));
            }

            // (a2) Discordant-only < 1000 bp DUP, DEL, and INV are ALL suppressed under sv_respect_min_len!
            {
                let mut vd_disc = crate::variations::VariationData::default();
                vd_disc.max_read_length = 50; // mindist = 400 so mlen = 700 > mindist but < 1000!
                let mut fdup = Sclip::new();
                fdup.start = 540;
                fdup.end = 584;
                fdup.mstart = 300;
                fdup.mend = 340;
                fdup.variation.cnt = 10;
                fdup.variation.mqmean = 600.0;
                vd_disc.sv_structures.svfdup.push(fdup);

                let mut fdel = Sclip::new();
                fdel.start = 100;
                fdel.end = 150;
                fdel.mstart = 860;
                fdel.mend = 910;
                fdel.variation.cnt = 10; // mlen = 860 - 150 - 4 = 706 < 1000
                fdel.variation.mqmean = 600.0;
                vd_disc.sv_structures.svfdel.push(fdel);

                let mut invf5 = Sclip::new();
                invf5.start = 100;
                invf5.end = 150;
                invf5.mstart = 700;
                invf5.mend = 750;
                invf5.variation.cnt = 5;
                invf5.variation.mqmean = 300.0;
                vd_disc.sv_structures.svfinv5.push(invf5);
                let mut invr5 = Sclip::new();
                invr5.start = 110;
                invr5.end = 160;
                invr5.mstart = 705;
                invr5.mend = 755;
                invr5.variation.cnt = 5; // bp = 130, pe = 727 -> len = 598 < 1000
                invr5.variation.mqmean = 300.0;
                vd_disc.sv_structures.svrinv5.push(invr5);

                find_dup_disc(
                    &mut vd_disc,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_1000,
                    &chr_lengths,
                )
                .unwrap();
                find_del_disc(&mut vd_disc, &mut reference, chr, &cfg_1000, &chr_lengths).unwrap();
                find_inv_disc(&mut vd_disc, &mut reference, chr, &cfg_1000, &chr_lengths).unwrap();
                assert!(vd_disc.insertion_variants.is_empty());
                assert!(vd_disc.non_insertion_variants.is_empty());
            }

            // (b) Split-read exact find_inv_sub & find_sv (`[700..=949]`, len = 250 < 1000) -> explicit complement_seq!
            let expected_inv_rc = complement_seq(&String::from_utf8_lossy(&ref_seq[699..949]));
            assert_eq!(expected_inv_rc.len(), 250);
            {
                let mut rc_clip = inv_930_949.to_vec();
                crate::util::reverse_complement(&mut rc_clip);

                // find_inv_sub (Side::Five, dir = 1): soft clip at 700 matching revcomp(930..=949) -> bp = 949, len = 250
                let mut vd_inv = crate::variations::VariationData::default();
                vd_inv.max_read_length = 100;
                vd_inv.ref_coverage.insert(700, 30);
                let mut inv5 = Sclip::new();
                inv5.start = 650;
                inv5.end = 700;
                inv5.mstart = 910;
                inv5.mend = 950;
                inv5.mlen = -250;
                inv5.soft.insert(700, 6);
                inv5.variation.cnt = 4;
                vd_inv.sv_structures.svfinv5.push(inv5);
                let mut sc3 = Sclip::new();
                sc3.sequence = Some(rc_clip.clone());
                sc3.variation.cnt = 6;
                sc3.variation.hicnt = 6;
                sc3.variation.qmean = 210.0;
                sc3.variation.mqmean = 360.0;
                sc3.variation.pmean = 120.0;
                vd_inv.soft_clips_3_end.insert(700, sc3);

                find_inv_sub(
                    0,
                    1,
                    Side::Five,
                    &mut vd_inv,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_1000,
                    &chr_lengths,
                )
                .unwrap();
                let inv_key = vd_inv
                    .non_insertion_variants
                    .get(&700)
                    .unwrap()
                    .map
                    .keys()
                    .find(|k| k.starts_with("-250^"))
                    .expect("expected -250^ explicit INV");
                assert_eq!(inv_key, &format!("-250^{expected_inv_rc}"));
                assert!(!inv_key.contains("<inv"));

                // find_sv (tmp3): same 3' soft clip at 700 without discordant cluster -> explicit -250^<rc>
                let mut vd_fsv = crate::variations::VariationData::default();
                vd_fsv.max_read_length = 100;
                let mut sc3_fsv = Sclip::new();
                sc3_fsv.sequence = Some(rc_clip);
                sc3_fsv.variation.cnt = 6;
                sc3_fsv.variation.hicnt = 6;
                sc3_fsv.variation.qmean = 210.0;
                sc3_fsv.variation.mqmean = 360.0;
                sc3_fsv.variation.pmean = 120.0;
                vd_fsv.soft_clips_3_end.insert(700, sc3_fsv);
                find_sv(
                    &mut vd_fsv,
                    &mut reference,
                    chr,
                    650,
                    1000,
                    &cfg_1000,
                    &chr_lengths,
                )
                .unwrap();
                let fsv_key = vd_fsv
                    .non_insertion_variants
                    .get(&700)
                    .unwrap()
                    .map
                    .keys()
                    .find(|k| k.starts_with("-250^"))
                    .expect("expected find_sv -250^ explicit INV");
                assert_eq!(fsv_key, &format!("-250^{expected_inv_rc}"));
                assert!(!fsv_key.contains("<inv"));

                // Unloaded interior gap [750..=850] (two disjoint LoadedRegions [1..=749] and [851..=1000],
                // empty config.reference) must safely fall back to symbolic <inv150> without panicking!
                let mut disjoint_ref = Reference::new();
                disjoint_ref
                    .loaded_regions
                    .push(crate::reference::LoadedRegion {
                        chr: chr.to_string(),
                        sequence_start: 1,
                        sequence_end: 749,
                    });
                disjoint_ref.region_sequences.push(ref_seq[0..749].to_vec());
                disjoint_ref
                    .loaded_regions
                    .push(crate::reference::LoadedRegion {
                        chr: chr.to_string(),
                        sequence_start: 851,
                        sequence_end: 1_000,
                    });
                disjoint_ref
                    .region_sequences
                    .push(ref_seq[850..1_000].to_vec());
                disjoint_ref
                    .seed
                    .insert_for_test(&inv_930_949[3..20], vec![933]);
                let mut vd_unloaded = crate::variations::VariationData::default();
                vd_unloaded.max_read_length = 100;
                let mut rc_clip2 = inv_930_949.to_vec();
                crate::util::reverse_complement(&mut rc_clip2);
                let mut sc3_u = Sclip::new();
                sc3_u.sequence = Some(rc_clip2);
                sc3_u.variation.cnt = 6;
                sc3_u.variation.hicnt = 6;
                vd_unloaded.soft_clips_3_end.insert(700, sc3_u);
                find_sv(
                    &mut vd_unloaded,
                    &mut disjoint_ref,
                    chr,
                    650,
                    1000,
                    &cfg_1000,
                    &chr_lengths,
                )
                .unwrap();
                let u_key = vd_unloaded
                    .non_insertion_variants
                    .get(&700)
                    .unwrap()
                    .map
                    .keys()
                    .find(|k| k.starts_with("-250^"))
                    .expect("expected symbolic fallback INV when interior span is unloaded");
                assert!(
                    u_key.contains("<inv150>"),
                    "disjoint unloaded span must fall back safely to symbolic <inv150>, got {u_key}"
                );

                // Disjoint unloaded interior gap [2401..=5549] on a 5,000-bp split-read DEL
                // ([1200..=6199], -L 10000):
                // 1) Without FASTA path (`config.reference` empty), ToVarsBuilder must fall back
                //    cleanly to symbolic `<DEL>` (never splicing a truncated REF!).
                // 2) With `config.reference` pointing to FASTA, `find_sv` loads the full deletion
                //    span and ToVarsBuilder materializes the exact 5,001-bp REF (anchor + 5,000 bp).
                let mut del_disjoint_ref = Reference::new();
                del_disjoint_ref
                    .loaded_regions
                    .push(crate::reference::LoadedRegion {
                        chr: chr.to_string(),
                        sequence_start: 1,
                        sequence_end: 2_400,
                    });
                del_disjoint_ref
                    .region_sequences
                    .push(ref_seq[0..2_400].to_vec());
                del_disjoint_ref
                    .loaded_regions
                    .push(crate::reference::LoadedRegion {
                        chr: chr.to_string(),
                        sequence_start: 5_550,
                        sequence_end: 6_999,
                    });
                del_disjoint_ref
                    .region_sequences
                    .push(ref_seq[5_549..6_999].to_vec());
                del_disjoint_ref
                    .seed
                    .insert_for_test(&ref_seq[6_199..6_216], vec![6_200]);

                let mk_split_del_vd = || {
                    let mut vd_del = crate::variations::VariationData::default();
                    vd_del.max_read_length = 100;
                    vd_del.ref_coverage.insert(1_200, 30);
                    let mut sc3_d = Sclip::new();
                    sc3_d.sequence = Some(ref_seq[6_199..6_220].to_vec());
                    sc3_d.variation.cnt = 6;
                    sc3_d.variation.hicnt = 6;
                    sc3_d.variation.qmean = 210.0;
                    sc3_d.variation.mqmean = 360.0;
                    sc3_d.variation.pmean = 120.0;
                    vd_del.soft_clips_3_end.insert(1_200, sc3_d);
                    let mut fdel = Sclip::new();
                    fdel.start = 1_110;
                    fdel.end = 1_199;
                    fdel.mstart = 6_200;
                    fdel.mend = 6_290;
                    fdel.variation.cnt = 6;
                    fdel.variation.hicnt = 6;
                    fdel.variation.qmean = 210.0;
                    fdel.variation.mqmean = 360.0;
                    fdel.variation.pmean = 120.0;
                    vd_del.sv_structures.svfdel.push(fdel);
                    vd_del
                };

                let cfg_unloaded_10k = Config {
                    min_reads: 2,
                    sv_min_len: 10_000,
                    compat: crate::config::COMPAT_OFF,
                    ..Config::default()
                };
                let mut vd_del_unloaded = mk_split_del_vd();
                find_sv(
                    &mut vd_del_unloaded,
                    &mut del_disjoint_ref,
                    chr,
                    1_000,
                    1_200,
                    &cfg_unloaded_10k,
                    &chr_lengths,
                )
                .unwrap();
                let del_reg = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(1_000),
                    end: Pos1(1_200),
                    gene: "DEL5K".to_string(),
                    amplicon: None,
                };
                let al_del_unloaded = crate::tovars::ToVarsBuilder::new(
                    &del_disjoint_ref,
                    &del_reg,
                    &chr_lengths,
                    &cfg_unloaded_10k,
                )
                .build_aligned_variants(&vd_del_unloaded);
                let v_unloaded = &al_del_unloaded
                    .get(&1_200)
                    .expect("expected DEL entry at 1200")
                    .variants[0];
                assert_eq!(
                    v_unloaded.var_allele, "<DEL>",
                    "unloaded interior gap must fall back to <DEL> rather than truncating REF"
                );
                assert_eq!(v_unloaded.ref_allele.len(), 1);
                assert_eq!(v_unloaded.start_pos, 1_199);
                assert_eq!(v_unloaded.end_pos, 6_199);

                let tmp_fa_dir = tempfile::tempdir().unwrap();
                let tmp_fa_path = tmp_fa_dir.path().join("del5k.fa");
                let mut fa_bytes = b">chr1\n".to_vec();
                fa_bytes.extend_from_slice(&ref_seq[0..8_000]);
                fa_bytes.push(b'\n');
                std::fs::write(&tmp_fa_path, &fa_bytes).unwrap();
                std::fs::write(
                    tmp_fa_dir.path().join("del5k.fa.fai"),
                    "chr1\t8000\t6\t8000\t8001\n",
                )
                .unwrap();
                let cfg_loaded_10k = Config {
                    reference: tmp_fa_path.to_string_lossy().into_owned(),
                    min_reads: 2,
                    sv_min_len: 10_000,
                    compat: crate::config::COMPAT_OFF,
                    ..Config::default()
                };
                let chr_len_8k = HashMap::from([(chr.to_string(), 8_000usize)]);
                let mut vd_del_loaded = mk_split_del_vd();
                find_sv(
                    &mut vd_del_loaded,
                    &mut del_disjoint_ref,
                    chr,
                    1_000,
                    1_200,
                    &cfg_loaded_10k,
                    &chr_len_8k,
                )
                .unwrap();
                let al_del_loaded = crate::tovars::ToVarsBuilder::new(
                    &del_disjoint_ref,
                    &del_reg,
                    &chr_len_8k,
                    &cfg_loaded_10k,
                )
                .build_aligned_variants(&vd_del_loaded);
                let v_loaded = &al_del_loaded
                    .get(&1_200)
                    .expect("expected loaded DEL entry at 1200")
                    .variants[0];
                assert_eq!(v_loaded.start_pos, 1_199);
                assert_eq!(v_loaded.end_pos, 6_199);
                assert_eq!(v_loaded.ref_allele.len(), 5_001);
                assert_eq!(
                    v_loaded.ref_allele.as_bytes(),
                    &ref_seq[1_198..6_199],
                    "loaded split-read DEL under -L 10000 must materialize all 5,001 REF bases"
                );
                assert_eq!(v_loaded.var_allele.as_bytes(), &ref_seq[1_198..1_199]);
            }

            // (c) -L 1000000 ceiling: 40,000-bp split-refined DUP ([1000..=40999], mlen=40000 >= 32768)
            //     must NOT materialize a 40 KB allele — falls back safely to symbolic <dup39900>!
            {
                let cfg_1m = Config {
                    min_reads: 2,
                    sv_min_len: 1_000_000,
                    compat: sw27,
                    ..Config::default()
                };
                let mut vd_ceil = crate::variations::VariationData::default();
                vd_ceil.max_read_length = 100;
                let mut fdup = Sclip::new();
                fdup.start = 40_950;
                fdup.end = 41_000;
                fdup.mstart = 1_005;
                fdup.mend = 1_055;
                fdup.softp = 41_000;
                fdup.soft.insert(41_000, 6);
                fdup.variation.cnt = 10;
                fdup.variation.mqmean = 600.0;
                vd_ceil.sv_structures.svfdup.push(fdup);
                let mut sc3 = Sclip::new();
                sc3.variation.cnt = 6;
                sc3.sequence = Some(ceil_1000_1019.to_vec());
                vd_ceil.soft_clips_3_end.insert(41_000, sc3);

                find_dup_disc(
                    &mut vd_ceil,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_1m,
                    &chr_lengths,
                )
                .unwrap();
                let (ceil_desc, _) = vd_ceil
                    .insertion_variants
                    .get(&1000)
                    .unwrap()
                    .map
                    .iter()
                    .next()
                    .unwrap();
                assert!(
                    ceil_desc.contains("<dup39900>"),
                    "40kb DUP >= MAX_EXPLICIT_SV_SEQ_BP must stay symbolic, got len={}",
                    ceil_desc.len()
                );
            }

            // (d) Contig start bp == 1 ([1..=200], mlen = 200 < 1000):
            //     anchored at bp = 1 with right-rotated 1 bp `+<dup_seq[1..]><ref_base_1>` so VCF at POS=1
            //     has REF = "G", ALT = "G" + dup_seq[1..] + "G", and ALT + ref[2..=200] == dup_seq + dup_seq!
            {
                let dup_1_200 = String::from_utf8_lossy(&ref_seq[0..200]).into_owned();
                assert_eq!(dup_1_200.len(), 200);
                let mut vd_bp1 = crate::variations::VariationData::default();
                vd_bp1.max_read_length = 100;
                vd_bp1.ref_coverage.insert(1, 30);
                let mut fdup = Sclip::new();
                fdup.start = 155;
                fdup.end = 201;
                fdup.mstart = 5;
                fdup.mend = 55;
                fdup.softp = 201;
                fdup.soft.insert(201, 6);
                fdup.variation.cnt = 10;
                fdup.variation.mqmean = 600.0;
                fdup.variation.qmean = 350.0;
                fdup.variation.pmean = 200.0;
                vd_bp1.sv_structures.svfdup.push(fdup);
                let mut sc3 = Sclip::new();
                sc3.variation.cnt = 6;
                sc3.variation.mqmean = 360.0;
                sc3.variation.qmean = 210.0;
                sc3.variation.pmean = 120.0;
                sc3.sequence = Some(bp1_anchor_1_20.to_vec());
                vd_bp1.soft_clips_3_end.insert(201, sc3);

                find_dup_disc(
                    &mut vd_bp1,
                    std::path::Path::new(""),
                    &mut reference,
                    chr,
                    &cfg_1000,
                    &chr_lengths,
                )
                .unwrap();
                assert!(!vd_bp1.insertion_variants.contains_key(&0));
                let (bp1_desc, _) = vd_bp1
                    .insertion_variants
                    .get(&1)
                    .expect("expected explicit DUP at anchor = 1")
                    .map
                    .iter()
                    .next()
                    .unwrap();
                assert_eq!(bp1_desc, &format!("+{}G", &dup_1_200[1..]));

                let aligned =
                    crate::tovars::ToVarsBuilder::new(&reference, &region, &chr_lengths, &cfg_1000)
                        .build_aligned_variants(&vd_bp1);
                let v1 = &aligned.get(&1).unwrap().variants[0];
                assert_eq!(v1.start_pos, 1);
                assert_eq!(v1.ref_allele, "G");
                // Mutant haplotype = ALT + ref[2..=200] == dup_seq + dup_seq!
                let mutant_hap = format!("{}{}", v1.var_allele, &dup_1_200[1..]);
                assert_eq!(mutant_hap, format!("{dup_1_200}{dup_1_200}"));
            }
        }

        #[test]
        fn wave6_lane1_switch33_sv_consensus_respects_majority_count_first() {
            let chr = "chr1";
            let ref_len = 2_000usize;
            let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);
            let (mut base_ref, mut ref_seq) = build_ref(chr, ref_len);
            // Replace 'T' with 'C' inside 1_100..1_125 so every true motif base < 'T'
            for b in &mut ref_seq[1_099..1_125] {
                if *b == b'T' {
                    *b = b'C';
                }
            }
            ref_seq[798] = b'A'; // 1-based 799
            ref_seq[1_098] = b'G'; // 1-based 1099 (prevents left-shift past 800..=1099)
            base_ref.region_sequences[0].clone_from(&ref_seq);
            base_ref
                .reference_sequences
                .insert(chr.to_string(), ref_seq.clone());
            let target_20bp = &ref_seq[1_099..1_119]; // 1-based 1100..=1119
            base_ref
                .seed
                .insert_for_test(&target_20bp[3..20], vec![1_103]);

            let run_sv_del = |majority_first: bool| -> bool {
                let mut reference = base_ref.clone();
                let cfg = Config {
                    min_reads: 2,
                    compat: crate::config::CompatProfile {
                        sclip_majority_count_first: majority_first,
                        ..crate::config::COMPAT_VARDICTJAVA_1_8_3
                    },
                    ..Config::default()
                };
                let mut vd = crate::variations::VariationData::default();
                vd.max_read_length = 100;
                vd.sv_structures_filtered = true;
                let mut bam_piled = vec![(1, ref_len as i32)];

                let mut fdel = Sclip::new();
                fdel.start = 750;
                fdel.end = 800;
                fdel.mstart = 1_100;
                fdel.mend = 1_150;
                fdel.variation.cnt = 5;
                fdel.soft.insert(800, 21);
                vd.sv_structures.svfdel.push(fdel);

                let mut sc3 = Sclip::new();
                sc3.variation.cnt = 21;
                sc3.variation.hicnt = 21;
                sc3.variation.qmean = 33.0 * 21.0;
                sc3.variation.mqmean = 60.0 * 21.0;
                sc3.variation.pmean = 25.0 * 21.0;
                for (idx, &true_b) in target_20bp.iter().enumerate() {
                    if idx == 4 || idx == 5 {
                        sc3.nt
                            .insert(idx as i32, BTreeMap::from([(true_b, 3), (b'T', 2)]));
                        let mut v_true = Variation::new();
                        v_true.cnt = 3;
                        v_true.qmean = 3.0 * 18.0;
                        let mut v_err = Variation::new();
                        v_err.cnt = 2;
                        v_err.qmean = 2.0 * 38.0;
                        sc3.seq.insert(
                            idx as i32,
                            BTreeMap::from([(true_b, v_true), (b'T', v_err)]),
                        );
                    } else {
                        sc3.nt
                            .insert(idx as i32, BTreeMap::from([(true_b, 20), (b'T', 1)]));
                        let mut v_true = Variation::new();
                        v_true.cnt = 20;
                        v_true.qmean = 20.0 * 33.0;
                        let mut v_err = Variation::new();
                        v_err.cnt = 1;
                        v_err.qmean = 1.0 * 38.0;
                        sc3.seq.insert(
                            idx as i32,
                            BTreeMap::from([(true_b, v_true), (b'T', v_err)]),
                        );
                    }
                }
                vd.soft_clips_3_end.insert(800, sc3);

                find_del_paths(
                    &mut vd,
                    &[],
                    &mut reference,
                    chr,
                    &cfg,
                    &chr_lengths,
                    &mut bam_piled,
                )
                .unwrap();
                vd.non_insertion_variants
                    .get(&800)
                    .is_some_and(|vm| vm.map.contains_key("-300"))
            };

            assert!(
                !run_sv_del(false),
                "legacy sclip_majority_count_first=false drops split-refined SV deletion due to single Q38 error read"
            );
            assert!(
                run_sv_del(true),
                "sclip_majority_count_first=true must preserve 20-read majority consensus in SV split-read refinement"
            );
        }

        #[test]
        fn pr1_lane2_f06_find_del_disc_lazy_mate_pileup_across_split_tiles() {
            use noodles::bam;
            use noodles::core::Position;
            use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
            use noodles::csi::binning_index::index::reference_sequence::index::LinearIndex;
            use noodles::csi::binning_index::Indexer;
            use noodles::sam;
            use noodles::sam::alignment::io::Write as _;
            use noodles::sam::alignment::record::cigar::op::Kind;
            use noodles::sam::alignment::record::cigar::Op;
            use noodles::sam::alignment::record::data::field::Tag;
            use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
            use noodles::sam::alignment::record_buf::{
                Cigar, Data, QualityScores, RecordBuf, Sequence,
            };
            use std::num::NonZeroUsize;

            // 1. Direct split-tile lazy mate-window pileup test for BOTH `svfdel` (distant mate
            //    `[2600, 2620]` contributing `soft_clips_5_end[2600] -> splits=4`) and `svrdel`
            //    (distant mate `[195, 200]` populating real mate-window `ref_coverage[205]==12`,
            //    wild-type `Variation` at `205`, and `merge_planted.contains(&205)` without
            //    duplicate `LoadedRegion` entries).
            let tmp_dir = tempfile::tempdir().unwrap();
            let chr = "chr1";
            let ref_len = 4_000usize;
            let motif = b"ACGTACGTAGCTAGCTTGCAACGTTAGCATGC";
            let mut ref_seq: Vec<u8> = (0..ref_len).map(|i| motif[i % motif.len()]).collect();
            ref_seq[2598] = b'A'; // 1-based 2599 != 'G' so 15S20M at 2600 stays anchored at 2600
            let fasta_path = tmp_dir.path().join("f06_del_disc.fa");
            let fai_path = tmp_dir.path().join("f06_del_disc.fa.fai");
            std::fs::write(
                &fasta_path,
                format!(">{chr}\n{}\n", String::from_utf8_lossy(&ref_seq)),
            )
            .unwrap();
            std::fs::write(
                &fai_path,
                format!("{chr}\t{ref_len}\t6\t{ref_len}\t{}\n", ref_len + 1),
            )
            .unwrap();

            let bam_path = tmp_dir.path().join("f06_del_disc.bam");
            let bai_path = tmp_dir.path().join("f06_del_disc.bam.bai");
            {
                let header = sam::Header::builder()
                    .add_reference_sequence(
                        chr,
                        sam::header::record::value::Map::<
                            sam::header::record::value::map::ReferenceSequence,
                        >::new(NonZeroUsize::new(ref_len).unwrap()),
                    )
                    .build();
                let mut writer = bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
                writer.write_header(&header).unwrap();
                let mut indexer = Indexer::<LinearIndex>::new(14, 5).unwrap();

                // 12 wild-type reads covering 195..=214 (covers svrdel breakpoint bp=205)
                // plus 4 SNV reads at 120..=139 (pos 130, inside [100, 200] mate window for E2E
                // multi-tile test without overlapping [195, 200] in subtest 1b)
                let alt_130 = if ref_seq[129] == b'A' { b'T' } else { b'A' };
                for k in 0..4 {
                    let mut seq_bytes = ref_seq[119..139].to_vec();
                    seq_bytes[10] = alt_130; // 1-based pos 130
                    let cigar: Cigar = [Op::new(Kind::Match, 20)].into_iter().collect();
                    let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(1i32))]
                        .into_iter()
                        .collect();
                    let rec = RecordBuf::builder()
                        .set_name(format!("snv130_{k}"))
                        .set_flags(sam::alignment::record::Flags::empty())
                        .set_reference_sequence_id(0)
                        .set_alignment_start(Position::try_from(120).unwrap())
                        .set_mapping_quality(
                            sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                        )
                        .set_cigar(cigar)
                        .set_sequence(Sequence::from(seq_bytes))
                        .set_quality_scores(QualityScores::from(vec![35u8; 20]))
                        .set_data(data)
                        .build();
                    let start_vp = writer.get_ref().virtual_position();
                    writer.write_alignment_record(&header, &rec).unwrap();
                    let end_vp = writer.get_ref().virtual_position();
                    indexer
                        .add_record(
                            Some((
                                0,
                                Position::try_from(120).unwrap(),
                                Position::try_from(139).unwrap(),
                                true,
                            )),
                            Chunk::new(start_vp, end_vp),
                        )
                        .unwrap();
                }

                for k in 0..12 {
                    let seq_bytes = ref_seq[194..214].to_vec();
                    let cigar: Cigar = [Op::new(Kind::Match, 20)].into_iter().collect();
                    let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(0i32))]
                        .into_iter()
                        .collect();
                    let rec = RecordBuf::builder()
                        .set_name(format!("wt205_{k}"))
                        .set_flags(sam::alignment::record::Flags::empty())
                        .set_reference_sequence_id(0)
                        .set_alignment_start(Position::try_from(195).unwrap())
                        .set_mapping_quality(
                            sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                        )
                        .set_cigar(cigar)
                        .set_sequence(Sequence::from(seq_bytes))
                        .set_quality_scores(QualityScores::from(vec![35u8; 20]))
                        .set_data(data)
                        .build();
                    let start_vp = writer.get_ref().virtual_position();
                    writer.write_alignment_record(&header, &rec).unwrap();
                    let end_vp = writer.get_ref().virtual_position();
                    indexer
                        .add_record(
                            Some((
                                0,
                                Position::try_from(195).unwrap(),
                                Position::try_from(214).unwrap(),
                                true,
                            )),
                            Chunk::new(start_vp, end_vp),
                        )
                        .unwrap();
                }

                // 8 discordant reverse reads covering 1450..=1549 with forward mates at 100..=199
                // (forms native `svrdel` in Tile 3 [1001, 1500] targeting `bp = 205`, `mlen = 1239`)
                for k in 0..8 {
                    let seq_bytes = ref_seq[1449..1549].to_vec();
                    let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
                    let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(0i32))]
                        .into_iter()
                        .collect();
                    let flags = sam::alignment::record::Flags::SEGMENTED
                        | sam::alignment::record::Flags::FIRST_SEGMENT
                        | sam::alignment::record::Flags::REVERSE_COMPLEMENTED;
                    let rec = RecordBuf::builder()
                        .set_name(format!("rdel1450_{k}"))
                        .set_flags(flags)
                        .set_reference_sequence_id(0)
                        .set_alignment_start(Position::try_from(1450).unwrap())
                        .set_mapping_quality(
                            sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                        )
                        .set_cigar(cigar)
                        .set_mate_reference_sequence_id(0)
                        .set_mate_alignment_start(Position::try_from(100).unwrap())
                        .set_template_length(-1450)
                        .set_sequence(Sequence::from(seq_bytes))
                        .set_quality_scores(QualityScores::from(vec![35u8; 100]))
                        .set_data(data)
                        .build();
                    let start_vp = writer.get_ref().virtual_position();
                    writer.write_alignment_record(&header, &rec).unwrap();
                    let end_vp = writer.get_ref().virtual_position();
                    indexer
                        .add_record(
                            Some((
                                0,
                                Position::try_from(1450).unwrap(),
                                Position::try_from(1549).unwrap(),
                                true,
                            )),
                            Chunk::new(start_vp, end_vp),
                        )
                        .unwrap();
                }

                // 4 5'-soft-clipped reads at 2600 (covers svfdel mstart=2600 -> splits=4)
                for k in 0..4 {
                    let mut seq_bytes = b"GGGTTCCAAAGGTTC".to_vec();
                    seq_bytes.extend_from_slice(&ref_seq[2599..2619]);
                    let cigar: Cigar = [Op::new(Kind::SoftClip, 15), Op::new(Kind::Match, 20)]
                        .into_iter()
                        .collect();
                    let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(0i32))]
                        .into_iter()
                        .collect();
                    let rec = RecordBuf::builder()
                        .set_name(format!("sc2600_{k}"))
                        .set_flags(sam::alignment::record::Flags::empty())
                        .set_reference_sequence_id(0)
                        .set_alignment_start(Position::try_from(2600).unwrap())
                        .set_mapping_quality(
                            sam::alignment::record::MappingQuality::try_from(60).unwrap(),
                        )
                        .set_cigar(cigar)
                        .set_sequence(Sequence::from(seq_bytes))
                        .set_quality_scores(QualityScores::from(vec![35u8; 35]))
                        .set_data(data)
                        .build();
                    let start_vp = writer.get_ref().virtual_position();
                    writer.write_alignment_record(&header, &rec).unwrap();
                    let end_vp = writer.get_ref().virtual_position();
                    indexer
                        .add_record(
                            Some((
                                0,
                                Position::try_from(2600).unwrap(),
                                Position::try_from(2619).unwrap(),
                                true,
                            )),
                            Chunk::new(start_vp, end_vp),
                        )
                        .unwrap();
                }

                writer.try_finish().unwrap();
                bam::bai::io::Writer::new(std::fs::File::create(&bai_path).unwrap())
                    .write_index(&indexer.build(1))
                    .unwrap();
            }

            let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);
            let primary_reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1_400),
                end: Pos1(1_700),
                gene: String::new(),
                amplicon: None,
            };
            let bam_paths = vec![bam_path.to_string_lossy().into_owned()];
            let cfg_off = Config {
                min_reads: 2,
                sv_min_len: 100,
                reference: fasta_path.to_string_lossy().into_owned(),
                compat: crate::config::COMPAT_OFF,
                ..Config::default()
            };

            // 1a. Pure discordant `svfdel` (`soft` empty) in [1450, 1500] with mate [2600, 2620]:
            //     `find_all_sv_paths` -> `find_del_disc_paths` must pile [2600, 2620] and recover `splits == 4`.
            {
                let mut reference = Reference::new();
                get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();
                let mut vd = crate::variations::VariationData::default();
                vd.sv_structures_filtered = true;
                vd.max_read_length = 100;
                vd.ref_coverage.insert(1_505, 20);

                let mut fdel = Sclip::new();
                fdel.start = 1_450;
                fdel.end = 1_500;
                fdel.mstart = 2_600;
                fdel.mend = 2_620;
                fdel.variation.cnt = 8; // >= min_reads + 5 (7), mlen = 2600 - 1500 - 11 = 1089 > mindist (800)
                fdel.variation.mqmean = 480.0;
                fdel.variation.qmean = 280.0;
                fdel.variation.pmean = 200.0;
                vd.sv_structures.svfdel.push(fdel);

                // Backward-compat wrapper `find_del_disc` without BAM paths gets splits == 0:
                let mut vd_no_paths = vd.clone();
                let mut ref_no_paths = Reference::new();
                get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut ref_no_paths).unwrap();
                find_del_disc(
                    &mut vd_no_paths,
                    &mut ref_no_paths,
                    chr,
                    &cfg_off,
                    &chr_lengths,
                )
                .unwrap();
                assert_eq!(
                    vd_no_paths
                        .non_insertion_variants
                        .get(&1_505)
                        .unwrap()
                        .sv
                        .as_ref()
                        .unwrap()
                        .splits,
                    0
                );

                // `find_all_sv_paths` -> `find_del_disc_paths` piles [2600, 2620] and recovers splits == 4:
                find_all_sv_paths(
                    &mut vd,
                    &bam_paths,
                    &mut reference,
                    chr,
                    ref_len as u32,
                    &primary_reg,
                    &cfg_off,
                    &chr_lengths,
                )
                .unwrap();
                let sv_fdel = vd
                    .non_insertion_variants
                    .get(&1_505)
                    .unwrap()
                    .sv
                    .as_ref()
                    .unwrap();
                assert_eq!(
                    sv_fdel.splits, 4,
                    "find_del_disc_paths svfdel must pile distant mate window [2600, 2620] and recover soft_clips_5_end[2600] splits=4"
                );
            }

            // 1b. Pure discordant `svrdel` (`soft` empty) in [1450, 1500] with mate [195, 200]:
            //     `bp = 200 + (100 / 9) / 2 = 205`; `find_del_disc_paths` must pile [195, 200] BEFORE
            //     `ref_coverage` check so `ref_coverage[205] == 12` (real wild-type coverage, not fabricated `2*8=16`),
            //     record `205` in `merge_planted`, and avoid duplicate `LoadedRegion` entries.
            {
                let mut reference = Reference::new();
                get_reference(&primary_reg, 0, &chr_lengths, &cfg_off, &mut reference).unwrap();
                let mut vd = crate::variations::VariationData::default();
                vd.sv_structures_filtered = true;
                vd.owning_range = Some((Pos1(1_400), Pos1(1_700)));
                vd.max_read_length = 100;

                let mut rdel = Sclip::new();
                rdel.start = 1_450;
                rdel.end = 1_500;
                rdel.mstart = 195;
                rdel.mend = 200;
                rdel.variation.cnt = 8; // mlen = 1450 - 200 - 11 = 1239 > mindist (800)
                rdel.variation.mqmean = 480.0;
                rdel.variation.qmean = 280.0;
                rdel.variation.pmean = 200.0;
                vd.sv_structures.svrdel.push(rdel);

                find_all_sv_paths(
                    &mut vd,
                    &bam_paths,
                    &mut reference,
                    chr,
                    ref_len as u32,
                    &primary_reg,
                    &cfg_off,
                    &chr_lengths,
                )
                .unwrap();

                assert_eq!(
                    vd.ref_coverage.get(&205),
                    Some(&12),
                    "find_del_disc_paths svrdel must read real mate-window ref_coverage[205] == 12 instead of fabricating 2*cnt == 16"
                );
                let ref_ch = (ref_seq[204] as char).to_string();
                assert_eq!(
                    vd.non_insertion_variants
                        .get(&205)
                        .unwrap()
                        .map
                        .get(&ref_ch)
                        .map(|v| v.cnt),
                    Some(12),
                    "svrdel breakpoint 205 must retain 12 piled wild-type reads for RD"
                );
                assert!(
                    vd.merge_planted.contains(&205),
                    "svrdel distant mate breakpoint 205 must be recorded in merge_planted"
                );
                assert_eq!(
                    reference.loaded_regions.len(),
                    2,
                    "svrdel lazy mate pileup must not append duplicate LoadedRegions"
                );

                // Also plant incidental unowned mate-window SNV + insertion at `205` (`tile_planted == false`)
                // and verify `ToVarsBuilder::build_aligned_variants` preserves `<DEL>` + 12 WT reads (`RD = 12`)
                // while dropping incidental unowned mate-window alleles!
                let sp_alt = if ref_seq[204] == b'A' { "T" } else { "A" };
                {
                    let sp_snv = crate::variations::get_variation(
                        &mut vd.non_insertion_variants,
                        205,
                        sp_alt,
                    );
                    sp_snv.cnt = 4;
                    sp_snv.hicnt = 4;
                    sp_snv.fwd = 2;
                    sp_snv.rev = 2;
                    sp_snv.qmean = 140.0;
                    sp_snv.mqmean = 240.0;
                    sp_snv.pmean = 40.0;
                    let sp_ins =
                        crate::variations::get_variation(&mut vd.insertion_variants, 205, "+ACGT");
                    sp_ins.cnt = 3;
                    sp_ins.hicnt = 3;
                    sp_ins.fwd = 2;
                    sp_ins.rev = 1;
                    sp_ins.qmean = 105.0;
                    sp_ins.mqmean = 180.0;
                    sp_ins.pmean = 30.0;
                }
                let full_reg = crate::regions::Region {
                    chr: chr.to_string(),
                    start: Pos1(1),
                    end: Pos1(2_200),
                    gene: String::new(),
                    amplicon: None,
                };
                let builder = crate::tovars::ToVarsBuilder::new(
                    &reference,
                    &full_reg,
                    &chr_lengths,
                    &cfg_off,
                );
                let aligned = builder.build_aligned_variants(&vd);
                let vars_205 = aligned.get(&205).expect(
                    "ToVarsBuilder must NOT drop tile-planted svrdel <DEL> at merge_planted pos 205",
                );
                assert_eq!(
                    vars_205.variants.len(),
                    1,
                    "ToVarsBuilder must emit only tile-planted <DEL> at 205 and drop unowned mate-window SNV/indel: {:?}",
                    vars_205
                        .variants
                        .iter()
                        .map(|v| (&v.description, &v.var_allele))
                        .collect::<Vec<_>>()
                );
                assert_eq!(vars_205.variants[0].var_allele, "<DEL>");
                assert_eq!(vars_205.variants[0].pos_coverage, 16);
                assert_eq!(
                    vars_205.variants[0].ref_fwd_cov + vars_205.variants[0].ref_rev_cov,
                    12,
                    "ToVarsBuilder must retain exact 12 wild-type reads (RD=12) on heterozygous <DEL> at 205"
                );

                // Full multi-tile `crate::run_with_writer` (`tile_size = Some(500)` vs untiled `Some(50_000_000)`)
                // on `f06_del_disc.bam` (`chr1:1-2200`):
                let run_f06 = |tile_sz: u32, out_vcf: bool| -> String {
                    let cfg = Config {
                        bam: bam_path.to_string_lossy().into_owned(),
                        reference: fasta_path.to_string_lossy().into_owned(),
                        region: Some("chr1:1-2200".to_string()),
                        min_reads: 2,
                        sv_min_len: 100,
                        freq: 0.01,
                        emit_intermediate: !out_vcf,
                        tile_size: Some(tile_sz),
                        threads: 1,
                        compat: crate::config::COMPAT_OFF,
                        ..Config::default()
                    };
                    let mut out = Vec::new();
                    crate::run_with_writer(&cfg, &mut out).unwrap();
                    String::from_utf8(out)
                        .unwrap()
                        .lines()
                        .filter(|l| !l.starts_with('#'))
                        .collect::<Vec<_>>()
                        .join("\n")
                };
                let untiled_vcf = run_f06(50_000_000, true);
                let tiled_vcf = run_f06(500, true);
                assert!(
                    untiled_vcf.contains("\t204\t.\tT\t<DEL>\t")
                        && untiled_vcf.contains(";DP=28;END=1443;VD=16;AF=0.5714;")
                        && untiled_vcf.contains("\t0/1:28:16:12,16:0.5714:12,0:8,8")
                        && untiled_vcf.contains("\t130\t.\t"),
                    "untiled run must emit SNV at 130 and heterozygous <DEL> at 204 with RD=12, DP=28, AF=0.5714:\n{untiled_vcf}"
                );
                assert_eq!(
                    tiled_vcf, untiled_vcf,
                    "multi-tile crate::run (tile_size=500) VCF must match untiled on heterozygous distant-mate <DEL> without dropping <DEL> or duplicating mate-window SNVs"
                );
                assert_eq!(
                    run_f06(500, false),
                    run_f06(50_000_000, false),
                    "multi-tile crate::run (tile_size=500) .var must match untiled on heterozygous distant-mate <DEL>"
                );
            }

            // 2. Real-BAM split-tile parity on `Colo829.chr10.DELDiscReverse.bam` (`chr10:89700300-89713976`,
            //    `--tile-size 200` vs `--tile-size 50000000`):
            let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .parent()
                .unwrap()
                .to_path_buf();
            let data_root = if root.join("VarDictJava/testdata/fastas/hg19.fa").exists() {
                root
            } else {
                std::path::PathBuf::from("/mnt/media_backups/VarDictRust")
            };
            let fa = data_root.join("VarDictJava/testdata/fastas/hg19.fa");
            let bam = data_root.join(
                "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr10.DELDiscReverse.bam",
            );
            if fa.exists() && bam.exists() {
                let run_with_tile = |tile_size: u32| -> String {
                    let out_file = tmp_dir.path().join(format!("deldisc_{tile_size}.vcf"));
                    let cfg = Config {
                        reference: fa.to_string_lossy().into_owned(),
                        bam: bam.to_string_lossy().into_owned(),
                        region: Some("chr10:89700300-89713976".to_string()),
                        freq: 0.001,
                        af_threshold: 0.001,
                        emit_intermediate: false,
                        tile_size: Some(tile_size),
                        threads: 1,
                        output_path: Some(out_file.to_string_lossy().into_owned()),
                        compat: crate::config::COMPAT_OFF,
                        ..Config::default()
                    };
                    crate::run(&cfg).unwrap();
                    std::fs::read_to_string(&out_file).unwrap()
                };
                let untiled = run_with_tile(50_000_000);
                let tiled = run_with_tile(200);
                let untiled_del = untiled
                    .lines()
                    .find(|l| !l.starts_with('#') && l.contains("<DEL>"))
                    .expect("untiled DELDiscReverse must emit <DEL>");
                let tiled_del = tiled
                    .lines()
                    .find(|l| !l.starts_with('#') && l.contains("<DEL>"))
                    .expect("tiled DELDiscReverse must emit <DEL>");
                assert_eq!(tiled_del, untiled_del);
            }
        }

        #[test]
        fn pr1_lane2_f14b_find_inv_disc_cutoff_unification_under_l50() {
            let chr = "chr1";
            let ref_len = 2_000usize;
            let chr_lengths = HashMap::from([(chr.to_string(), ref_len)]);
            let (mut reference, mut ref_seq) = build_ref(chr, ref_len);

            let mut vd = crate::variations::VariationData::default();
            vd.max_read_length = 50;
            vd.ref_coverage.insert(500, 30);
            // Construct overlapping svfinv5 and svrinv5 with 80-bp inversion (bp=500, pe=579, len=80):
            let mut invf = Sclip::new();
            invf.start = 460;
            invf.end = 500;
            invf.mstart = 540;
            invf.mend = 579;
            invf.variation.cnt = 5;
            invf.variation.mqmean = 60.0 * 5.0;
            invf.variation.qmean = 35.0 * 5.0;
            invf.variation.pmean = 20.0 * 5.0;
            vd.sv_structures.svfinv5.push(invf);

            let mut invr = Sclip::new();
            invr.start = 500;
            invr.end = 540;
            invr.mstart = 579;
            invr.mend = 619;
            invr.variation.cnt = 5;
            invr.variation.mqmean = 60.0 * 5.0;
            invr.variation.qmean = 35.0 * 5.0;
            invr.variation.pmean = 20.0 * 5.0;
            vd.sv_structures.svrinv5.push(invr);

            let cfg = Config {
                min_reads: 2,
                sv_min_len: 50,
                compat: crate::config::CompatProfile {
                    sv_respect_min_len: true,
                    ..crate::config::COMPAT_VARDICTJAVA_1_8_3
                },
                ..Config::default()
            };
            let mut piled = Vec::new();
            find_inv_disc_paths(
                &mut vd,
                &[],
                &mut reference,
                chr,
                &cfg,
                &chr_lengths,
                &mut piled,
            )
            .unwrap();

            let keys: Vec<String> = vd
                .non_insertion_variants
                .get(&500)
                .map(|vm| vm.map.keys().map(str::to_string).collect())
                .unwrap_or_default();
            assert!(
                keys.iter().any(|k| k.contains("<inv-20>")),
                "80-bp INV under -L 50 with sv_respect_min_len=true must stay symbolic <inv-20>, got {keys:?}"
            );

            // Run ToVarsBuilder under BOTH isolated `sv_respect_min_len=true` (`canonical_sv_vcf_anchoring=false`)
            // and full `COMPAT_OFF` to prove `<inv-20>` yields `var_allele == "<INV>"` and non-empty 1-bp `ref_allele`:
            let reg = crate::regions::Region {
                chr: chr.to_string(),
                start: Pos1(1),
                end: Pos1(ref_len as u32),
                gene: "INV80".to_string(),
                amplicon: None,
            };
            let al_iso = crate::tovars::ToVarsBuilder::new(&reference, &reg, &chr_lengths, &cfg)
                .build_aligned_variants(&vd);
            let v_iso = &al_iso.get(&500).expect("expected INV at 500").variants[0];
            assert_eq!(v_iso.var_allele, "<INV>");
            assert_eq!(
                v_iso.ref_allele.len(),
                1,
                "isolated sv_respect_min_len must populate 1-bp ref_allele on <inv-20>, got {:?}",
                v_iso.ref_allele
            );

            let cfg_off_l50 = Config {
                min_reads: 2,
                sv_min_len: 50,
                compat: crate::config::COMPAT_OFF,
                ..Config::default()
            };
            let al_off =
                crate::tovars::ToVarsBuilder::new(&reference, &reg, &chr_lengths, &cfg_off_l50)
                    .build_aligned_variants(&vd);
            let v_off = &al_off.get(&500).expect("expected INV at 500").variants[0];
            assert_eq!(v_off.var_allele, "<INV>");
            assert_eq!(v_off.ref_allele.len(), 1);

            // Also verify `find_sv` soft-clip-only 80-bp inversions (`51..=100 bp`, both 3' and 5' paths):
            // 1-based 700..=779 (span = 80), distinct boundary bases at 699 ('A') and 780 ('A'):
            let motif_760_779 = b"CCTTAAGGCCTTAAGGCCAT";
            ref_seq[698] = b'A'; // pos 699
            ref_seq[779] = b'A'; // pos 780
            ref_seq[759..779].copy_from_slice(motif_760_779); // 1-based 760..=779
            let motif_700_719 = b"GGAATTCCGGAATTCCGGTA";
            ref_seq[699..719].copy_from_slice(motif_700_719); // 1-based 700..=719
            reference.region_sequences[0].clone_from(&ref_seq);
            reference
                .reference_sequences
                .insert(chr.to_string(), ref_seq.clone());
            reference
                .seed
                .insert_for_test(&motif_760_779[3..20], vec![763]);
            reference
                .seed
                .insert_for_test(&motif_700_719[0..17], vec![700]);

            // 3' soft-clip at 700 matching RC(760..=779) -> bp=779, span=80:
            let mut vd_fsv3 = crate::variations::VariationData::default();
            vd_fsv3.max_read_length = 100;
            vd_fsv3.ref_coverage.insert(700, 25);
            let mut rc_3 = motif_760_779.to_vec();
            crate::util::reverse_complement(&mut rc_3);
            let mut sc3_700 = Sclip::new();
            sc3_700.sequence = Some(rc_3);
            sc3_700.variation.cnt = 6;
            sc3_700.variation.hicnt = 6;
            sc3_700.variation.qmean = 210.0;
            sc3_700.variation.mqmean = 360.0;
            sc3_700.variation.pmean = 120.0;
            vd_fsv3.soft_clips_3_end.insert(700, sc3_700);
            find_sv(
                &mut vd_fsv3,
                &mut reference,
                chr,
                650,
                850,
                &cfg,
                &chr_lengths,
            )
            .unwrap();
            let fsv3_keys: Vec<String> = vd_fsv3
                .non_insertion_variants
                .get(&700)
                .map(|vm| vm.map.keys().map(str::to_string).collect())
                .unwrap_or_default();
            assert!(
                fsv3_keys.iter().any(|k| k.contains("<inv-20>")),
                "find_sv 3' 80-bp soft-clip INV under -L 50 must emit symbolic <inv-20>, got {fsv3_keys:?}"
            );
            let al_fsv3 = crate::tovars::ToVarsBuilder::new(&reference, &reg, &chr_lengths, &cfg)
                .build_aligned_variants(&vd_fsv3);
            let v_fsv3 = &al_fsv3.get(&700).unwrap().variants[0];
            assert_eq!(v_fsv3.var_allele, "<INV>");
            assert_eq!(v_fsv3.ref_allele.len(), 1);

            // 5' soft-clip at 780 (`p5 = 780`, `dir = -1`) matching RC(700..=719) -> bp=700, swapped p5_adj=700, bp=779, span=80:
            let mut vd_fsv5 = crate::variations::VariationData::default();
            vd_fsv5.max_read_length = 100;
            vd_fsv5.ref_coverage.insert(700, 25);
            vd_fsv5.ref_coverage.insert(780, 25);
            // In find_match_rev(..., dir=-1), seq is complemented in-place (not reversed!), so `complement(seq)` == `motif_700_719`:
            let comp_5: Vec<u8> = motif_700_719
                .iter()
                .map(|&b| crate::util::complement(b))
                .collect();
            let mut sc5_780 = Sclip::new();
            sc5_780.sequence = Some(comp_5);
            sc5_780.variation.cnt = 6;
            sc5_780.variation.hicnt = 6;
            sc5_780.variation.qmean = 210.0;
            sc5_780.variation.mqmean = 360.0;
            sc5_780.variation.pmean = 120.0;
            vd_fsv5.soft_clips_5_end.insert(780, sc5_780);
            find_sv(
                &mut vd_fsv5,
                &mut reference,
                chr,
                650,
                850,
                &cfg,
                &chr_lengths,
            )
            .unwrap();
            let fsv5_keys: Vec<String> = vd_fsv5
                .non_insertion_variants
                .get(&700)
                .map(|vm| vm.map.keys().map(str::to_string).collect())
                .unwrap_or_default();
            assert!(
                fsv5_keys.iter().any(|k| k.contains("<inv-20>")),
                "find_sv 5' 80-bp soft-clip INV under -L 50 must emit symbolic <inv-20>, got {fsv5_keys:?}"
            );
            let al_fsv5 = crate::tovars::ToVarsBuilder::new(&reference, &reg, &chr_lengths, &cfg)
                .build_aligned_variants(&vd_fsv5);
            let v_fsv5 = &al_fsv5.get(&700).unwrap().variants[0];
            assert_eq!(v_fsv5.var_allele, "<INV>");
            assert_eq!(v_fsv5.ref_allele.len(), 1);
        }
    }
}
