//! ToVarsBuilder parses intermediate VariationData and finalizes into Variant records.
//! Ports ToVarsBuilder.java.

use crate::config::Config;
use crate::reference::Reference;
use crate::regions::Region;
use crate::types::round_half_even_compat;
use crate::variations::{Variant, VariationData};
use std::collections::HashMap;

/// Flanking bases captured per side around `<dup±N>` tokens
/// (`StructuralVariantsProcessor.java:1430, 1565`, `ToVarsBuilder.java:619, 699-700`:
/// `Configuration.SVFLANK = 50`).
const SV_FLANK_BP: i64 = 50;

#[derive(Debug, Clone)]
pub struct MsiResult {
    pub msi: f64,
    pub shift3: i32,
    pub msint: String,
}

/// Compact summary of a pure reference-only site for deferred Somatic hydration.
#[derive(Debug, Clone)]
pub struct CompactRefSite {
    pub total_cov: i32,
    pub duprate: f64,
    pub cnt: crate::variations::Variation,
}

/// Controls how pure reference-only sites are handled during
/// [`ToVarsBuilder::build_aligned_variants_with_ref_mode`].
pub enum RefSiteMode<'a> {
    KeepAll,
    DropAll,
    KeepWhere(&'a dyn Fn(i32) -> bool),
    CompactOut(&'a mut rustc_hash::FxHashMap<i32, CompactRefSite>),
}

impl std::fmt::Debug for RefSiteMode<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::KeepAll => f.write_str("KeepAll"),
            Self::DropAll => f.write_str("DropAll"),
            Self::KeepWhere(_) => f.debug_tuple("KeepWhere").field(&"<fn>").finish(),
            Self::CompactOut(map) => f.debug_tuple("CompactOut").field(map).finish(),
        }
    }
}

/// ToVarsBuilder converts intermediate VariationData into finalized Variant structures.
pub struct ToVarsBuilder<'a> {
    pub reference: &'a Reference,
    pub region: &'a Region,
    pub chromosomes_lengths: &'a HashMap<String, usize>,
    pub config: &'a Config,
    primary_ref_seq: &'a [u8],
    primary_ref_start: u32,
    chr_len: i32,
}

/// G3 (mate-window ownership) single-owner rule for a merge-planted key.
///
/// The tile whose active interval contains the position owns it; planted
/// keys outside the called region clamp to the first (resp. last) tile —
/// the same edge-clamp principle as `tile_owns_splice_key` (Task 3).
/// `None` (manually built data) fails open. Single-region runs stamp the
/// whole called region, so every planted key is owned and the gate never
/// fires (Java parity by construction).
fn merge_key_owned(
    owning_range: Option<(crate::types::Pos1, crate::types::Pos1)>,
    builder_region: &Region,
    position: i32,
) -> bool {
    let Some((start, end)) = owning_range else {
        return true;
    };
    if start.0 as i32 <= position && position <= end.0 as i32 {
        return true;
    }
    if position < builder_region.start.0 as i32 && start == builder_region.start {
        return true;
    }
    if position > builder_region.end.0 as i32 && end == builder_region.end {
        return true;
    }
    false
}

#[inline]
fn clamp_nonneg_pos_zero(v: f64) -> f64 {
    if v > 0.0 {
        v
    } else {
        0.0
    }
}

/// Reconciles `(cnt, fwd, rev, hicnt)` under `CompatProfile::reconcile_strand_and_hq_sums`
/// (`AstraZeneca-NGS/VarDictJava#316`, `clintval/vartovcf#120`, `clintval/vartovcf#142`)
/// into `(vd, fwd, rev, hicnt, locnt)` without signed `.clamp(0, cnt)` panics or `0 / 0`.
#[inline]
pub(crate) fn reconcile_strand_and_hq(
    cnt: i32,
    fwd: i32,
    rev: i32,
    hicnt: i32,
) -> (i32, i32, i32, i32, i32) {
    let vd = cnt.max(0);
    let hicnt_eff = hicnt.max(0).min(vd);
    let locnt_eff = vd.saturating_sub(hicnt_eff);
    let mut fwd_eff = fwd.max(0);
    let mut rev_eff = rev.max(0);
    if vd == 0 {
        fwd_eff = 0;
        rev_eff = 0;
    } else if fwd_eff == 0 && rev_eff == 0 {
        fwd_eff = vd;
        rev_eff = 0;
    } else {
        let sum = i64::from(fwd_eff) + i64::from(rev_eff);
        let vd64 = i64::from(vd);
        if sum < vd64 {
            let add = (vd64 - sum) as i32;
            if fwd_eff >= rev_eff {
                fwd_eff = fwd_eff.saturating_add(add);
            } else {
                rev_eff = rev_eff.saturating_add(add);
            }
        } else if sum > vd64 {
            let excess = sum - vd64;
            if fwd_eff >= rev_eff {
                let sub_major = excess.min(i64::from(fwd_eff));
                fwd_eff -= sub_major as i32;
                rev_eff -= (excess - sub_major) as i32;
            } else {
                let sub_major = excess.min(i64::from(rev_eff));
                rev_eff -= sub_major as i32;
                fwd_eff -= (excess - sub_major) as i32;
            }
        }
    }
    (vd, fwd_eff, rev_eff, hicnt_eff, locnt_eff)
}

#[derive(Clone, Copy)]
struct CommonFieldFlags {
    exact_binary_half_even: bool,
    reconcile_dp: bool,
    reconcile_strand_hq: bool,
    nonnegative_subtracted_means: bool,
}

/// Fill the fields shared by every create_variant/create_insertion record
/// (ToVarsBuilder.java's common tvref initialization).
#[allow(clippy::too_many_arguments)]
fn fill_common_fields(
    tvref: &mut Variant,
    description: &str,
    cnt: &crate::variations::Variation,
    fwd: i32,
    rev: i32,
    bias: i32,
    ttcov: i32,
    base_quality: f64,
    mapping_quality: f64,
    hicnt: i32,
    locnt: i32,
    hicov: i32,
    duprate: f64,
    flags: CommonFieldFlags,
) {
    tvref.description = description.to_string();
    tvref.fwd = fwd;
    tvref.rev = rev;
    tvref.bias.clear();
    let _ = std::fmt::Write::write_fmt(&mut tvref.bias, format_args!("{}", bias));
    let mut pmean =
        round_half_even_compat(cnt.pmean / cnt.cnt as f64, 1, flags.exact_binary_half_even);
    let mut qmean = base_quality;
    let mut mqmean = mapping_quality;
    let mut nm = round_half_even_compat(cnt.nm / cnt.cnt as f64, 1, flags.exact_binary_half_even);
    if flags.nonnegative_subtracted_means {
        pmean = clamp_nonneg_pos_zero(pmean);
        qmean = clamp_nonneg_pos_zero(qmean);
        mqmean = clamp_nonneg_pos_zero(mqmean);
        nm = clamp_nonneg_pos_zero(nm);
    }
    tvref.mean_position = pmean;
    tvref.pstd = cnt.pstd;
    tvref.mean_quality = qmean;
    tvref.qstd = cnt.qstd;
    tvref.mean_mapping_quality = mqmean;
    tvref.hq_lq_ratio = hicnt as f64 / if locnt != 0 { locnt as f64 } else { 0.5 };
    tvref.shift3 = 0;
    tvref.msi = 0.0;
    tvref.nm = nm;
    tvref.hicnt = hicnt;
    tvref.hicov = hicov;
    tvref.duprate = duprate;
    if flags.reconcile_dp {
        let pos_cov_eff = cnt.cnt.max(0);
        let ttcov_eff = ttcov.max(pos_cov_eff).max(0);
        tvref.pos_coverage = pos_cov_eff;
        tvref.frequency = if ttcov_eff > 0 && pos_cov_eff > 0 {
            round_half_even_compat(
                pos_cov_eff as f64 / ttcov_eff as f64,
                4,
                flags.exact_binary_half_even,
            )
        } else {
            0.0
        };
        tvref.hi_freq = if hicov > 0 && hicnt > 0 {
            hicnt as f64 / hicov as f64
        } else {
            0.0
        };
        tvref.extra_freq = if cnt.extracnt > 0 && ttcov_eff > 0 {
            (cnt.extracnt.min(ttcov_eff) as f64) / (ttcov_eff as f64)
        } else {
            0.0
        };
        tvref.total_pos_cov = ttcov_eff;
    } else if flags.reconcile_strand_hq {
        let pos_cov_eff = cnt.cnt.max(0);
        tvref.pos_coverage = pos_cov_eff;
        tvref.frequency = if pos_cov_eff > 0 && ttcov > 0 {
            round_half_even_compat(
                pos_cov_eff as f64 / ttcov as f64,
                4,
                flags.exact_binary_half_even,
            )
        } else {
            0.0
        };
        tvref.hi_freq = if hicov > 0 && hicnt > 0 {
            hicnt as f64 / hicov as f64
        } else {
            0.0
        };
        tvref.extra_freq = if cnt.extracnt != 0 {
            cnt.extracnt as f64 / ttcov as f64
        } else {
            0.0
        };
    } else {
        tvref.pos_coverage = cnt.cnt;
        tvref.frequency = round_half_even_compat(
            cnt.cnt as f64 / ttcov as f64,
            4,
            flags.exact_binary_half_even,
        );
        tvref.hi_freq = if hicov > 0 {
            hicnt as f64 / hicov as f64
        } else {
            0.0
        };
        tvref.extra_freq = if cnt.extracnt != 0 {
            cnt.extracnt as f64 / ttcov as f64
        } else {
            0.0
        };
    }
}

enum SortedKeySlice<'a> {
    Stack([&'a str; 8], usize),
    Heap(Vec<&'a str>),
}

impl<'a> SortedKeySlice<'a> {
    #[inline]
    fn new(iter: impl ExactSizeIterator<Item = &'a str>, extra: Option<&'a str>) -> Self {
        let iter_len = iter.len();
        let total_len = iter_len.saturating_add(usize::from(extra.is_some()));
        if total_len <= 8 {
            let mut arr = [""; 8];
            for (slot, item) in arr.iter_mut().zip(iter) {
                *slot = item;
            }
            if let Some(extra_key) = extra {
                if let Some(slot) = arr.get_mut(iter_len) {
                    *slot = extra_key;
                }
            }
            if let Some(s) = arr.get_mut(..total_len) {
                s.sort_unstable();
            }
            Self::Stack(arr, total_len)
        } else {
            let mut v = Vec::with_capacity(total_len);
            v.extend(iter);
            if let Some(extra_key) = extra {
                v.push(extra_key);
            }
            v.sort_unstable();
            Self::Heap(v)
        }
    }

    #[inline]
    fn as_slice(&self) -> &[&'a str] {
        match self {
            Self::Stack(arr, len) => arr.get(..*len).unwrap_or(&[]),
            Self::Heap(v) => v.as_slice(),
        }
    }
}

enum PositionIter<I> {
    Bucket(I),
    Sorted(std::vec::IntoIter<i32>),
}

impl<'m, I: Iterator<Item = &'m i32>> Iterator for PositionIter<I> {
    type Item = i32;

    #[inline]
    fn next(&mut self) -> Option<i32> {
        match self {
            Self::Bucket(it) => it.next().copied(),
            Self::Sorted(it) => it.next(),
        }
    }
}

#[derive(Clone, Copy)]
struct NonInsBuildCtx {
    position: i32,
    total_pos_coverage: i32,
    hicov: i32,
    duprate: f64,
}

#[derive(Clone, Copy)]
struct RefSiteCounts {
    fwd: i32,
    rev: i32,
    extracnt: i32,
}

#[inline]
fn lookup_raw_extracnt(
    position: i32,
    desc: &str,
    insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    orig_non_ins: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    overlay: &rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
) -> i32 {
    insertion_variants
        .get(&position)
        .and_then(|m| m.map.get(desc))
        .or_else(|| {
            overlay
                .get(&position)
                .or_else(|| orig_non_ins.get(&position))
                .and_then(|m| m.map.get(desc))
        })
        .map_or(0, |c| c.extracnt)
}

#[inline]
fn is_sw26_iupac_ref_key(key: &str, eff_ref_base: u8) -> bool {
    if key == "I" || key == "SV" {
        return false;
    }
    let [b] = *key.as_bytes() else {
        return false;
    };
    reduce_ref_byte_sw26(b) == eff_ref_base
}

#[inline]
fn merge_variation_into_ref(
    dst: &mut crate::variations::Variation,
    src: &crate::variations::Variation,
) {
    if dst.cnt == 0 && dst.extracnt == 0 && dst.fwd == 0 && dst.rev == 0 {
        *dst = src.clone();
        return;
    }
    if src.cnt > 0 {
        dst.pstd = dst.pstd || src.pstd || dst.pp != src.pp;
        dst.qstd = dst.qstd || src.qstd || dst.pq != src.pq;
    }
    dst.cnt = dst.cnt.saturating_add(src.cnt);
    dst.fwd = dst.fwd.saturating_add(src.fwd);
    dst.rev = dst.rev.saturating_add(src.rev);
    dst.hicnt = dst.hicnt.saturating_add(src.hicnt);
    dst.locnt = dst.locnt.saturating_add(src.locnt);
    dst.extracnt = dst.extracnt.saturating_add(src.extracnt);
    dst.pmean += src.pmean;
    dst.qmean += src.qmean;
    dst.mqmean += src.mqmean;
    dst.nm += src.nm;
    dst.overlap_agreed_pairs = dst
        .overlap_agreed_pairs
        .saturating_add(src.overlap_agreed_pairs);
}

fn merge_iupac_ref_alleles_in_place(
    vm: &mut crate::variations::VariationMap,
    eff_ref_base: u8,
    eff_ref_str: &str,
) {
    let alias_keys: Vec<String> = vm
        .map
        .keys()
        .filter(|k| *k != eff_ref_str && is_sw26_iupac_ref_key(k, eff_ref_base))
        .map(str::to_string)
        .collect();
    if alias_keys.is_empty() {
        return;
    }
    for k in alias_keys {
        if let Some(alias_var) = vm.map.shift_remove(&k) {
            let dst = vm.map.get_or_insert_variation(eff_ref_str);
            merge_variation_into_ref(dst, &alias_var);
        }
    }
}

#[inline]
fn lookup_raw_ref_extracnt(
    position: i32,
    desc: &str,
    eff_ref_base: Option<u8>,
    sanitize_iupac: bool,
    insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    orig_non_ins: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    overlay: &rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
) -> i32 {
    if sanitize_iupac {
        if let Some(eff_ref) = eff_ref_base {
            if is_sw26_iupac_ref_key(desc, eff_ref) {
                if let Some(ov) = overlay.get(&position) {
                    return ov.map.get(desc).map_or(0, |c| c.extracnt);
                }
                if let Some(orig) = orig_non_ins.get(&position) {
                    let mut sum: i32 = 0;
                    for (k, v) in &orig.map {
                        if is_sw26_iupac_ref_key(k, eff_ref) {
                            sum = sum.saturating_add(v.extracnt);
                        }
                    }
                    return sum;
                }
                return 0;
            }
        }
    }
    lookup_raw_extracnt(position, desc, insertion_variants, orig_non_ins, overlay)
}

#[inline]
fn lookup_pre_anchored_dup(
    position: i32,
    desc: &str,
    insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    orig_non_ins: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    overlay: &rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
) -> bool {
    insertion_variants
        .get(&position)
        .and_then(|m| m.map.get(desc))
        .or_else(|| {
            overlay
                .get(&position)
                .or_else(|| orig_non_ins.get(&position))
                .and_then(|m| m.map.get(desc))
        })
        .is_some_and(|c| c.pre_anchored_dup)
}

#[inline]
fn long_ins_inserted_bases(stripped_desc: &str, description_string: &str) -> i32 {
    if !stripped_desc
        .bytes()
        .any(|b| matches!(b, b'&' | b'#' | b'^'))
    {
        return (stripped_desc.len().min(i32::MAX as usize) as i32).max(1);
    }
    let head_len = stripped_desc
        .bytes()
        .take_while(|b| !matches!(b, b'&' | b'#' | b'^'))
        .count();
    let tail_ins_len = if let Some((_, tail)) = get_hash_group_caret_group(description_string) {
        if get_begin_digits(tail).is_some() {
            0
        } else {
            tail.bytes().take_while(|&b| b != b'&').count()
        }
    } else if let Some(caret_ins) = get_caret_atgnc(description_string) {
        caret_ins.len()
    } else {
        0
    };
    let total = head_len.saturating_add(tail_ins_len).min(i32::MAX as usize) as i32;
    total.max(1)
}

impl<'a> ToVarsBuilder<'a> {
    /// Create a new ToVarsBuilder with reference, region, chromosomes_lengths and config.
    pub fn new(
        reference: &'a Reference,
        region: &'a Region,
        chromosomes_lengths: &'a HashMap<String, usize>,
        config: &'a Config,
    ) -> Self {
        let (primary_ref_seq, primary_ref_start) = reference
            .primary_region_slice(&region.chr)
            .unwrap_or((&[], 1));
        let chr_len = chromosomes_lengths
            .get(&region.chr)
            .map_or(0, |len| (*len).min(i32::MAX as usize) as i32);
        Self {
            reference,
            region,
            chromosomes_lengths,
            config,
            primary_ref_seq,
            primary_ref_start,
            chr_len,
        }
    }

    #[inline]
    fn get_ref_base(&self, pos: u32) -> Option<u8> {
        if let Some(offset) = pos.checked_sub(self.primary_ref_start) {
            if let Some(&b) = self.primary_ref_seq.get(offset as usize) {
                return Some(b);
            }
        }
        self.reference.get_base(&self.region.chr, pos)
    }

    #[inline]
    fn get_ref_base_i32(&self, pos: i32) -> Option<u8> {
        u32::try_from(pos).ok().and_then(|pu| self.get_ref_base(pu))
    }

    #[inline]
    fn effective_ref_base_i32(&self, pos: i32) -> Option<u8> {
        self.get_ref_base_i32(pos).map(|b| {
            if self.config.compat.sanitize_iupac_alt_alleles {
                reduce_ref_byte_sw26(b)
            } else {
                b
            }
        })
    }

    #[inline]
    fn normalize_var_map_iupac_ref<'m>(
        &self,
        position: i32,
        var_map: &'m crate::variations::VariationMap,
    ) -> std::borrow::Cow<'m, crate::variations::VariationMap> {
        if !self.config.compat.sanitize_iupac_alt_alleles {
            return std::borrow::Cow::Borrowed(var_map);
        }
        let Some(eff_ref_base) = self.effective_ref_base_i32(position) else {
            return std::borrow::Cow::Borrowed(var_map);
        };
        let mut buf = [0u8; 4];
        let eff_ref_str: &str = (eff_ref_base as char).encode_utf8(&mut buf);
        let has_alias = var_map
            .map
            .keys()
            .any(|k| k != eff_ref_str && is_sw26_iupac_ref_key(k, eff_ref_base));
        if !has_alias {
            return std::borrow::Cow::Borrowed(var_map);
        }
        let mut owned = var_map.clone();
        merge_iupac_ref_alleles_in_place(&mut owned, eff_ref_base, eff_ref_str);
        std::borrow::Cow::Owned(owned)
    }

    /// Returns the 1-bp-upstream reference base at `anchor_pos - 1` iff `anchor_pos > 1`
    /// and `anchor_pos - 1` is loaded in `self.reference` with a non-`N`/`n` base
    /// (`CompatProfile::canonical_sv_vcf_anchoring`; `AstraZeneca-NGS/VarDictJava#398`,
    /// `clintval/vartovcf#192`).
    #[inline]
    fn upstream_non_n_anchor_base(&self, anchor_pos: i32) -> Option<u8> {
        if anchor_pos <= 1 {
            return None;
        }
        self.effective_ref_base_i32(anchor_pos - 1)
            .filter(|&b| b != b'N' && b != b'n')
    }

    #[inline]
    fn is_explicit_del_span_loaded(&self, lo: i32, hi: i32) -> bool {
        if lo < 1 || hi < lo {
            return false;
        }
        let lo_u32 = lo as u32;
        let hi_u32 = hi as u32;
        crate::reference::is_loaded(&self.region.chr, lo_u32, hi_u32, self.reference)
            && self.get_ref_base(lo_u32).is_some()
            && self.get_ref_base(hi_u32).is_some()
    }

    /// Convert VariationData into a map of position -> Vars container.
    pub fn build_aligned_variants(
        &self,
        var_data: &VariationData,
    ) -> crate::javamap::JavaHashMap<i32, crate::variations::Vars> {
        let ref_mode = if !self.config.bam.contains('|')
            && !self.config.pileup
            && self.config.amplicon.is_none()
        {
            RefSiteMode::DropAll
        } else {
            RefSiteMode::KeepAll
        };
        self.build_aligned_variants_with_ref_mode(var_data, ref_mode)
    }

    /// Evaluate and finalize a single `target_pos` from `var_data`.
    pub fn build_vars_at_position(
        &self,
        var_data: &VariationData,
        target_pos: i32,
    ) -> Option<crate::variations::Vars> {
        self.build_vars_at_position_with_ref_mode(var_data, target_pos, RefSiteMode::KeepAll)
    }

    /// Evaluate and finalize a single `target_pos` from `var_data` under `ref_mode`.
    ///
    /// Matches `self.build_aligned_variants_with_ref_mode(var_data, ref_mode).get(&target_pos).cloned()`
    /// while evaluating only `target_pos` (and `target_pos - 1`'s `overlay[target_pos]` insertion
    /// contribution when `target_pos - 1` precedes `target_pos` in `JavaHashMap` iteration order under
    /// default `--compat=vardictjava`, or ascending coordinate order under `--compat=off`).
    pub(crate) fn build_vars_at_position_with_ref_mode(
        &self,
        var_data: &VariationData,
        target_pos: i32,
        mut ref_mode: RefSiteMode<'_>,
    ) -> Option<crate::variations::Vars> {
        if !var_data.non_insertion_variants.contains_key(&target_pos) {
            return None;
        }
        let duprate = var_data.duprate;
        let mut overlay: rustc_hash::FxHashMap<i32, crate::variations::VariationMap> =
            rustc_hash::FxHashMap::default();

        if let Some(prev_pos) = target_pos.checked_sub(1) {
            let prev_precedes = if self.config.compat.coordinate_sorted_vars {
                var_data.non_insertion_variants.contains_key(&prev_pos)
            } else {
                var_data
                    .non_insertion_variants
                    .contains_both_and_first_precedes(&prev_pos, &target_pos)
            };
            if prev_pos >= 1 && var_data.insertion_variants.contains_key(&prev_pos) && prev_precedes
            {
                if let Some((_, prev_total_cov)) =
                    self.check_position_guards(var_data, &overlay, prev_pos, &mut ref_mode)
                {
                    self.apply_prev_insertion_overlay(
                        &var_data.insertion_variants,
                        &var_data.non_insertion_variants,
                        &mut overlay,
                        &var_data.ref_coverage,
                        prev_pos,
                        prev_total_cov,
                    );
                }
            }
        }

        let (raw_var_map, total_cov) =
            self.check_position_guards(var_data, &overlay, target_pos, &mut ref_mode)?;
        let norm_var_map = self.normalize_var_map_iupac_ref(target_pos, raw_var_map);
        let var_map: &crate::variations::VariationMap = &norm_var_map;

        let ins_at_pos = var_data.insertion_variants.get(&target_pos);
        let hicov = self.calc_hicov(ins_at_pos, var_map);

        let map_len = var_map.map.len() + usize::from(var_map.sv.is_some());
        let ins_len = ins_at_pos.map_or(0, |m| m.map.len());
        let mut var = Vec::with_capacity(map_len + ins_len);
        let mut debug_lines = Vec::new();
        let mut target_vars = crate::variations::Vars::new();

        {
            let sorted_keys = SortedKeySlice::new(
                var_map.map.keys(),
                if var_map.sv.is_some() {
                    Some("SV")
                } else {
                    None
                },
            );
            self.push_non_insertion_variants(
                NonInsBuildCtx {
                    position: target_pos,
                    total_pos_coverage: total_cov,
                    hicov,
                    duprate,
                },
                var_map,
                &mut var,
                &mut debug_lines,
                sorted_keys.as_slice(),
                |sv| {
                    target_vars.sv = format!("{}-{}-{}", sv.splits, sv.pairs, sv.clusters);
                },
            );
        }
        drop(norm_var_map);

        let current_total_cov = self.create_insertion(
            &var_data.insertion_variants,
            &var_data.non_insertion_variants,
            &mut overlay,
            &var_data.ref_coverage,
            duprate,
            target_pos,
            total_cov,
            &mut var,
            &mut debug_lines,
            hicov,
        );

        self.sort_variants(&mut var);
        let maxfreq = self.collect_vars_into(&mut target_vars, target_pos, &mut var);

        if self.should_drop_low_freq(maxfreq) {
            return None;
        }

        let mut scratch1 = String::new();
        let mut scratch2 = String::new();
        self.collect_reference_variants(
            target_pos,
            current_total_cov,
            &mut target_vars,
            &var_data.ref_coverage,
            &var_data.insertion_variants,
            &var_data.non_insertion_variants,
            &overlay,
            duprate,
            &debug_lines,
            &mut scratch1,
            &mut scratch2,
        );

        if var_data.non_insertion_variants.treeified() || var_data.insertion_variants.treeified() {
            log::debug!(
                "JavaHashMap OpenJDK TreeNode bin active for chr {}:{}-{}",
                crate::util::SanitizedBytes(self.region.chr.as_bytes()),
                self.region.start,
                self.region.end
            );
        }

        Some(target_vars)
    }

    /// Convert VariationData into a map of position -> Vars container using `ref_mode`
    /// to elide or compact pure reference-only sites.
    pub fn build_aligned_variants_with_ref_mode(
        &self,
        var_data: &VariationData,
        ref_mode: RefSiteMode<'_>,
    ) -> crate::javamap::JavaHashMap<i32, crate::variations::Vars> {
        self.build_aligned_variants_inner(var_data, ref_mode, None)
    }

    /// Convert `VariationData` into a map of `position -> Vars` under `ref_mode`, clearing
    /// doomed sub-threshold non-reference allele lists at positions selected by `prune_doomed_at`
    /// before Phase-2 flank/MSI/genotype finalization while preserving `JavaHashMap` key order.
    pub(crate) fn build_aligned_variants_pruned(
        &self,
        var_data: &VariationData,
        ref_mode: RefSiteMode<'_>,
        prune_doomed_at: &dyn Fn(i32) -> bool,
    ) -> crate::javamap::JavaHashMap<i32, crate::variations::Vars> {
        self.build_aligned_variants_inner(var_data, ref_mode, Some(prune_doomed_at))
    }

    fn build_aligned_variants_inner(
        &self,
        var_data: &VariationData,
        mut ref_mode: RefSiteMode<'_>,
        prune_doomed_at: Option<&dyn Fn(i32) -> bool>,
    ) -> crate::javamap::JavaHashMap<i32, crate::variations::Vars> {
        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        let duprate = var_data.duprate;
        let can_prune_doomed = !self.config.pileup
            && !self.config.debug
            && self.config.amplicon.is_none()
            && prune_doomed_at.is_some();
        let min_reads = i32::try_from(self.config.min_reads).unwrap_or(i32::MAX);

        // Sparse copy-on-write overlay for positions mutated by create_insertion (position + 1);
        // lookup-only, never iterated.
        let mut overlay: rustc_hash::FxHashMap<i32, crate::variations::VariationMap> =
            rustc_hash::FxHashMap::default();

        let num_positions = var_data.non_insertion_variants.len();
        // Default --compat=vardictjava uses JavaHashMap iteration order; --compat=off
        // (CompatProfile::coordinate_sorted_vars) iterates positions in ascending coordinate order.
        match ref_mode {
            RefSiteMode::KeepAll => {
                aligned_variants.reserve_nodes(num_positions);
            }
            RefSiteMode::CompactOut(ref mut out) => {
                let region_span = (self.region.end.0.saturating_sub(self.region.start.0) as usize)
                    .saturating_add(1);
                out.reserve(num_positions.min(region_span));
            }
            RefSiteMode::DropAll | RefSiteMode::KeepWhere(_) => {}
        }

        // Task 6: per-position row buffers, allocated once per build and
        // cleared per position below (`clear()` retains capacity). Bound 4
        // covers the common 1-2 allele shape (MNP ties included); larger
        // positions regrow as before. `debug_lines` stays `new()`: it is
        // only pushed when `config.debug`, so the common path allocates
        // nothing. Order/skip semantics untouched — capacity only.
        let mut var = Vec::with_capacity(4);
        let mut debug_lines = Vec::new();
        let mut scratch1 = String::new();
        let mut scratch2 = String::new();

        // Process non-insertion variants (JavaHashMap bucket order by default; ascending
        // coordinate order under `CompatProfile::coordinate_sorted_vars`).
        let positions = if self.config.compat.coordinate_sorted_vars {
            let mut sorted: Vec<i32> = var_data.non_insertion_variants.keys().copied().collect();
            sorted.sort_unstable();
            PositionIter::Sorted(sorted.into_iter())
        } else {
            PositionIter::Bucket(var_data.non_insertion_variants.keys())
        };
        for position in positions {
            // `raw_var_map`'s shared NLL borrow of `overlay` and `var_data.non_insertion_variants`
            // ends at `drop(norm_var_map)` below, before `&mut overlay` is passed to `create_insertion`.
            let Some((raw_var_map, total_cov)) =
                self.check_position_guards(var_data, &overlay, position, &mut ref_mode)
            else {
                continue;
            };
            let norm_var_map = self.normalize_var_map_iupac_ref(position, raw_var_map);
            let var_map: &crate::variations::VariationMap = &norm_var_map;

            let ins_at_pos = var_data.insertion_variants.get(&position);
            let hicov = self.calc_hicov(ins_at_pos, var_map);

            let map_len = var_map.map.len() + usize::from(var_map.sv.is_some());
            let ins_len = ins_at_pos.map_or(0, |m| m.map.len());
            var.clear();
            var.reserve(map_len + ins_len);
            debug_lines.clear();

            {
                let sorted_keys = SortedKeySlice::new(
                    var_map.map.keys(),
                    if var_map.sv.is_some() {
                        Some("SV")
                    } else {
                        None
                    },
                );
                self.create_variant(
                    duprate,
                    &mut aligned_variants,
                    position,
                    var_map,
                    total_cov,
                    &mut var,
                    &mut debug_lines,
                    sorted_keys.as_slice(),
                    hicov,
                );
            }
            drop(norm_var_map);

            let current_total_cov = self.create_insertion(
                &var_data.insertion_variants,
                &var_data.non_insertion_variants,
                &mut overlay,
                &var_data.ref_coverage,
                duprate,
                position,
                total_cov,
                &mut var,
                &mut debug_lines,
                hicov,
            );

            let ref_base = self.effective_ref_base_i32(position).unwrap_or(b'N');
            let mut ref_buf = [0u8; 4];
            let ref_base_str: &str = (ref_base as char).encode_utf8(&mut ref_buf);

            let mut maxfreq = 0.0;
            let mut non_ref_rows = 0usize;
            let mut all_non_ref_doomed = true;
            for tvar in &var {
                if tvar.description != ref_base_str {
                    non_ref_rows = non_ref_rows.saturating_add(1);
                    if tvar.frequency > maxfreq {
                        maxfreq = tvar.frequency;
                    }
                    if all_non_ref_doomed
                        && !(tvar.hicnt < min_reads
                            || tvar.mean_position < self.config.read_pos
                            || tvar.mean_quality < self.config.base_qual)
                    {
                        all_non_ref_doomed = false;
                    }
                }
            }

            if self.should_drop_low_freq(maxfreq) {
                if !var.is_empty() {
                    let _ = crate::variations::get_or_put_vars(&mut aligned_variants, position);
                }
                aligned_variants.remove(&position);
                var.clear();
                continue;
            }

            if can_prune_doomed
                && non_ref_rows > 0
                && prune_doomed_at.is_some_and(|pred| pred(position))
                && all_non_ref_doomed
            {
                let _ = crate::variations::get_or_put_vars(&mut aligned_variants, position);
                var.clear();
                continue;
            }

            self.sort_variants(&mut var);
            let entry = crate::variations::get_or_put_vars(&mut aligned_variants, position);
            let _ = Self::drain_vars_into(entry, ref_base_str, non_ref_rows, &mut var);
            self.collect_reference_variants(
                position,
                current_total_cov,
                entry,
                &var_data.ref_coverage,
                &var_data.insertion_variants,
                &var_data.non_insertion_variants,
                &overlay,
                duprate,
                &debug_lines,
                &mut scratch1,
                &mut scratch2,
            );
        }

        if var_data.non_insertion_variants.treeified()
            || var_data.insertion_variants.treeified()
            || aligned_variants.treeified()
        {
            log::debug!(
                "JavaHashMap OpenJDK TreeNode bin active for chr {}:{}-{}",
                crate::util::SanitizedBytes(self.region.chr.as_bytes()),
                self.region.start,
                self.region.end
            );
        }
        aligned_variants
    }

    #[inline]
    fn should_drop_low_freq(&self, maxfreq: f64) -> bool {
        !self.config.pileup
            && maxfreq <= self.config.freq
            && self.config.amplicon.is_none()
            && !self.config.bam.contains('|')
    }

    /// Evaluates per-position skip guards (`ToVarsBuilder.java:103-128`, G3 `merge_planted`,
    /// and [`RefSiteMode`] elision/compaction), returning `(var_map, total_cov)` when kept.
    fn check_position_guards<'m>(
        &self,
        var_data: &'m VariationData,
        overlay: &'m rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
        position: i32,
        ref_mode: &mut RefSiteMode<'_>,
    ) -> Option<(&'m crate::variations::VariationMap, i32)> {
        let var_map = overlay
            .get(&position)
            .or_else(|| var_data.non_insertion_variants.get(&position))?;

        // Skip position if there are no variants on position
        if var_map.map.is_empty() && !var_data.insertion_variants.contains_key(&position) {
            return None;
        }

        // G3 (mate-window ownership): a key planted by a mate-window
        // merge emits only from its owning tile. Provenance comes from
        // the sv-path `merge_planted` set or the realign-path per-map
        // flag (set when the merge plants outside the merging tile's
        // range); the clamp rule owns every key in single-region runs,
        // so out-of-BED SV rows keep full mate ref-support there.
        if (var_data.merge_planted.contains(&position) || var_map.merge_planted)
            && !merge_key_owned(var_data.owning_range, self.region, position)
        {
            return None;
        }

        if (var_map.sv.is_none() || self.config.deldupvar)
            && !u32::try_from(position)
                .is_ok_and(|pu| pu >= self.region.start.0 && pu <= self.region.end.0)
        {
            return None;
        }

        if var_map.sv.is_none() && !var_data.ref_coverage.contains_key(&position) {
            return None;
        }

        let total_cov = var_data.ref_coverage.get(&position).copied().unwrap_or(0);
        match ref_mode {
            RefSiteMode::KeepAll => {}
            RefSiteMode::DropAll => {
                if self.is_pure_ref_variation(position, var_map, &var_data.insertion_variants) {
                    return None;
                }
            }
            RefSiteMode::KeepWhere(pred) => {
                if !self.config.debug
                    && var_map.sv.is_none()
                    && !var_data.insertion_variants.contains_key(&position)
                    && self.is_pure_ref_variation(position, var_map, &var_data.insertion_variants)
                    && !pred(position)
                {
                    return None;
                }
            }
            RefSiteMode::CompactOut(out) => {
                if !self.config.debug
                    && var_map.sv.is_none()
                    && !var_data.insertion_variants.contains_key(&position)
                    && !var_data
                        .insertion_variants
                        .contains_key(&position.saturating_sub(1))
                    && total_cov > 0
                {
                    let norm_map = self.normalize_var_map_iupac_ref(position, var_map);
                    if norm_map.map.len() == 1 {
                        if let (Some(sole_key), Some(ref_base)) = (
                            norm_map.map.keys().next(),
                            self.effective_ref_base_i32(position),
                        ) {
                            let mut buf = [0u8; 4];
                            let ref_base_str = (ref_base as char).encode_utf8(&mut buf);
                            if sole_key != "I" && sole_key != "SV" && sole_key == ref_base_str {
                                if let Some(cnt) = norm_map.map.get(sole_key) {
                                    if cnt.cnt > 0 {
                                        out.insert(
                                            position,
                                            CompactRefSite {
                                                total_cov,
                                                duprate: var_data.duprate,
                                                cnt: cnt.clone(),
                                            },
                                        );
                                        return None;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if total_cov == 0 {
            return None;
        }

        Some((var_map, total_cov))
    }

    /// Materialize a `Vars` container with its finalized `ref_variant` from a [`CompactRefSite`].
    pub fn hydrate_compact_ref_site(
        &self,
        position: i32,
        site: &CompactRefSite,
    ) -> crate::variations::Vars {
        let ref_base = self.effective_ref_base_i32(position).unwrap_or(b'N');
        let mut buf = [0u8; 4];
        let sole_key = (ref_base as char).encode_utf8(&mut buf);
        let mut vref = self.build_ref_variant_from_cnt(
            position,
            site.total_cov,
            site.duprate,
            &site.cnt,
            sole_key,
            &[],
        );
        if self.config.compat.nonnegative_subtracted_means {
            vref.mean_position = clamp_nonneg_pos_zero(vref.mean_position);
            vref.mean_quality = clamp_nonneg_pos_zero(vref.mean_quality);
            vref.mean_mapping_quality = clamp_nonneg_pos_zero(vref.mean_mapping_quality);
            vref.nm = clamp_nonneg_pos_zero(vref.nm);
        }
        let mut vars = crate::variations::Vars::new();
        vars.ref_variant = Some(vref);
        vars
    }

    #[cfg(test)]
    fn is_the_same_variation_on_ref(
        &self,
        position: i32,
        vars_at_cur_position: &crate::variations::VariationMap,
        insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    ) -> bool {
        if self.config.pileup || self.config.bam.contains('|') || self.config.amplicon.is_some() {
            return false;
        }
        self.is_pure_ref_variation(position, vars_at_cur_position, insertion_variants)
    }

    fn is_pure_ref_variation(
        &self,
        position: i32,
        vars_at_cur_position: &crate::variations::VariationMap,
        insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    ) -> bool {
        let has_ins = insertion_variants.contains_key(&position);
        if self.config.compat.sanitize_iupac_alt_alleles {
            if has_ins || vars_at_cur_position.map.is_empty() {
                return false;
            }
            let Some(eff_ref_base) = self.effective_ref_base_i32(position) else {
                return false;
            };
            return vars_at_cur_position
                .map
                .keys()
                .all(|k| is_sw26_iupac_ref_key(k, eff_ref_base));
        }
        let sole_key: &str = match (vars_at_cur_position.map.len(), has_ins) {
            (1, false) => match vars_at_cur_position.map.keys().next() {
                Some(k) => k,
                None => return false,
            },
            (0, true) => "I",
            (1, true) => match vars_at_cur_position.map.keys().next() {
                Some("I") => "I",
                _ => return false,
            },
            _ => return false,
        };
        if let Some(ref_base) = self.get_ref_base_i32(position) {
            let mut buf = [0u8; 4];
            let ref_base_str = (ref_base as char).encode_utf8(&mut buf);
            return sole_key == ref_base_str;
        }
        false
    }

    fn calc_hicov(
        &self,
        _insertion_variations: Option<&crate::variations::VariationMap>,
        non_insertion_variations: &crate::variations::VariationMap,
    ) -> i32 {
        let reconcile_strand_hq = self.config.compat.reconcile_strand_and_hq_sums;
        let clamp_nonneg =
            self.config.compat.reconcile_vd_dp || self.config.compat.reconcile_rd_vd_dp;
        let mut hicov: i32 = 0;
        for (key, val) in &non_insertion_variations.map {
            if key == "SV" || key.starts_with('+') {
                continue;
            }
            let h = if reconcile_strand_hq {
                let vd = val.cnt.max(0);
                val.hicnt.max(0).min(vd)
            } else if clamp_nonneg {
                val.hicnt.max(0)
            } else {
                val.hicnt
            };
            hicov = hicov.saturating_add(h);
        }
        hicov
    }

    /// ports VariationUtils.java:163-175
    #[inline]
    fn join_ref_into(&self, out: &mut String, start: i32, end: i32) {
        out.clear();
        if end < start || end < 1 {
            return;
        }
        let Some((min_s, max_e)) = self.reference.loaded_span_bounds(&self.region.chr) else {
            return;
        };
        if u64::from(min_s) > (end.max(0) as u64) {
            return;
        }
        let min_i32 = i32::try_from(min_s).unwrap_or(i32::MAX);
        let max_i32 = i32::try_from(max_e).unwrap_or(i32::MAX);
        let eff_start = start.max(1).max(min_i32);
        let eff_end = end.min(max_i32);
        if eff_end < eff_start {
            return;
        }
        if (eff_end as i64) - (eff_start as i64) > crate::util::MAX_CONTIGUOUS_REF_SCAN
            && self.reference.loaded_regions.len() > 1
        {
            let (subspans, span) = crate::util::merged_loaded_subspans(
                self.reference,
                &self.region.chr,
                eff_start,
                eff_end,
            );
            if out.capacity() == 0 {
                out.reserve_exact(span);
            } else {
                out.reserve(span);
            }
            for (s, e) in subspans {
                for p in s..=e {
                    if let Some(b) = self.get_ref_base_i32(p) {
                        if b.is_ascii() {
                            out.push(b as char);
                        } else {
                            crate::util::ref_bases_into(
                                out,
                                self.reference,
                                &self.region.chr,
                                eff_start,
                                eff_end,
                            );
                            return;
                        }
                    }
                }
            }
            return;
        }
        let span = ((eff_end as i64) - (eff_start as i64) + 1).clamp(0, 1 << 20) as usize;
        if out.capacity() == 0 {
            out.reserve_exact(span);
        } else {
            out.reserve(span);
        }
        for p in eff_start..=eff_end {
            if let Some(b) = self.get_ref_base_i32(p) {
                if b.is_ascii() {
                    out.push(b as char);
                } else {
                    crate::util::ref_bases_into(
                        out,
                        self.reference,
                        &self.region.chr,
                        eff_start,
                        eff_end,
                    );
                    return;
                }
            }
        }
    }

    fn proceed_vref_is_deletion(
        &self,
        position: i32,
        dellen: i32,
        leftseq: &mut String,
        tseq: &mut String,
    ) -> (f64, i32, i32) {
        self.join_ref_into(
            leftseq,
            std::cmp::max(position.saturating_sub(70), 1),
            position.saturating_sub(1),
        );
        let chr_len = self.chr_len;
        self.join_ref_into(
            tseq,
            position,
            std::cmp::min(position.saturating_add(dellen).saturating_add(70), chr_len),
        );

        let tseq_bytes = tseq.as_bytes();
        let dellen_usize = dellen.max(0) as usize;
        let tseq1 = tseq_bytes.get(..dellen_usize).unwrap_or(tseq_bytes);
        let tseq2 = tseq_bytes.get(dellen_usize..).unwrap_or(&[]);

        let (mut msi, shift3, mut msint_len) = find_msi_metrics(tseq1, tseq2, leftseq.as_bytes());
        let (tmsi, _, tmsint_len) = find_msi_metrics(leftseq.as_bytes(), tseq2, &[]);
        if msi < tmsi {
            msi = tmsi;
            msint_len = tmsint_len;
        }
        if dellen > 0 && msi <= (shift3 as f64) / (dellen as f64) {
            msi = (shift3 as f64) / (dellen as f64);
        }

        (msi, shift3, msint_len)
    }

    fn proceed_vref_is_insertion(
        &self,
        position: i32,
        vn: &str,
        leftseq: &mut String,
        tseq2: &mut String,
    ) -> (f64, i32, i32) {
        let tseq1 = vn.strip_prefix('+').unwrap_or(vn);
        self.join_ref_into(
            leftseq,
            std::cmp::max(position.saturating_sub(50), 1),
            position,
        );
        let chr_len = self.chr_len;
        if let Some(right_start) = position.checked_add(1) {
            self.join_ref_into(
                tseq2,
                right_start,
                std::cmp::min(position.saturating_add(70), chr_len),
            );
        } else {
            tseq2.clear();
        }

        let (mut msi, shift3, mut msint_len) =
            find_msi_metrics(tseq1.as_bytes(), tseq2.as_bytes(), leftseq.as_bytes());
        let (tmsi, _, tmsint_len) = find_msi_metrics(leftseq.as_bytes(), tseq2.as_bytes(), &[]);
        if msi < tmsi {
            msi = tmsi;
            msint_len = tmsint_len;
        }
        if !tseq1.is_empty() && msi <= (shift3 as f64) / (tseq1.len() as f64) {
            msi = (shift3 as f64) / (tseq1.len() as f64);
        }

        (msi, shift3, msint_len)
    }

    fn build_variant_from_cnt(
        &self,
        description_string: &str,
        cnt: &crate::variations::Variation,
        total_pos_coverage: i32,
        hicov: i32,
        duprate: f64,
        is_ref_allele: bool,
    ) -> Variant {
        let exact_binary_half_even = self.config.compat.exact_binary_half_even;
        let reconcile_strand_hq = self.config.compat.reconcile_strand_and_hq_sums;
        let reconcile_dp =
            self.config.compat.reconcile_vd_dp || self.config.compat.reconcile_rd_vd_dp;
        let (_vd, fwd, rev, hicnt_eff, locnt) = if reconcile_strand_hq {
            if is_ref_allele {
                let vd = cnt.cnt.max(0);
                let hc = cnt.hicnt.max(0).min(vd);
                (
                    vd,
                    cnt.fwd.max(0),
                    cnt.rev.max(0),
                    hc,
                    vd.saturating_sub(hc),
                )
            } else {
                reconcile_strand_and_hq(cnt.cnt, cnt.fwd, cnt.rev, cnt.hicnt)
            }
        } else {
            let hc = if reconcile_dp {
                cnt.hicnt.max(0)
            } else {
                cnt.hicnt
            };
            (cnt.cnt, cnt.fwd, cnt.rev, hc, cnt.locnt)
        };
        let min_reads_bias = i32::try_from(self.config.min_reads_bias).unwrap_or(i32::MAX);
        let bias = strand_bias(fwd, rev, min_reads_bias);
        let base_quality =
            round_half_even_compat(cnt.qmean / cnt.cnt as f64, 1, exact_binary_half_even);
        let mapping_quality =
            round_half_even_compat(cnt.mqmean / cnt.cnt as f64, 1, exact_binary_half_even);

        let (ttcov, hicov_eff) = if reconcile_dp {
            let pc = cnt.cnt.max(0);
            (
                total_pos_coverage.max(pc).max(0),
                hicov.max(hicnt_eff).max(0),
            )
        } else {
            let mut ttcov = total_pos_coverage;
            if cnt.cnt > total_pos_coverage
                && cnt.extracnt > 0
                && i64::from(cnt.cnt) - i64::from(total_pos_coverage) < i64::from(cnt.extracnt)
            {
                ttcov = cnt.cnt;
            }
            let hc_cov = if reconcile_strand_hq {
                hicov.max(hicnt_eff).max(0)
            } else {
                hicov
            };
            (ttcov, hc_cov)
        };

        let mut tvref = Variant::empty_internal();
        fill_common_fields(
            &mut tvref,
            description_string,
            cnt,
            fwd,
            rev,
            bias,
            ttcov,
            base_quality,
            mapping_quality,
            hicnt_eff,
            locnt,
            hicov_eff,
            duprate,
            CommonFieldFlags {
                exact_binary_half_even,
                reconcile_dp,
                reconcile_strand_hq,
                nonnegative_subtracted_means: self.config.compat.nonnegative_subtracted_means,
            },
        );
        tvref
    }

    fn build_ref_variant_from_cnt(
        &self,
        position: i32,
        total_cov: i32,
        duprate: f64,
        cnt: &crate::variations::Variation,
        sole_key: &str,
        debug_lines: &[String],
    ) -> Variant {
        let hicov = if sole_key == "SV" || sole_key.starts_with('+') {
            0
        } else if self.config.compat.reconcile_strand_and_hq_sums {
            let vd = cnt.cnt.max(0);
            cnt.hicnt.max(0).min(vd)
        } else if self.config.compat.reconcile_vd_dp || self.config.compat.reconcile_rd_vd_dp {
            cnt.hicnt.max(0)
        } else {
            cnt.hicnt
        };
        let mut tvref = self.build_variant_from_cnt(sole_key, cnt, total_cov, hicov, duprate, true);
        self.update_ref_variant(
            position,
            total_cov,
            &mut tvref,
            RefSiteCounts {
                fwd: cnt.fwd,
                rev: cnt.rev,
                extracnt: cnt.extracnt,
            },
            duprate,
            debug_lines,
        );
        tvref
    }

    fn create_variant(
        &self,
        duprate: f64,
        aligned_vars: &mut crate::javamap::JavaHashMap<i32, crate::variations::Vars>,
        position: i32,
        non_insertion_variations: &crate::variations::VariationMap,
        total_pos_coverage: i32,
        var: &mut Vec<Variant>,
        debug_lines: &mut Vec<String>,
        keys: &[&str],
        hicov: i32,
    ) {
        self.push_non_insertion_variants(
            NonInsBuildCtx {
                position,
                total_pos_coverage,
                hicov,
                duprate,
            },
            non_insertion_variations,
            var,
            debug_lines,
            keys,
            |sv| {
                let entry = crate::variations::get_or_put_vars(aligned_vars, position);
                entry.sv = format!("{}-{}-{}", sv.splits, sv.pairs, sv.clusters);
            },
        );
    }

    fn push_non_insertion_variants<F>(
        &self,
        ctx: NonInsBuildCtx,
        non_insertion_variations: &crate::variations::VariationMap,
        var: &mut Vec<Variant>,
        debug_lines: &mut Vec<String>,
        keys: &[&str],
        mut on_sv: F,
    ) where
        F: FnMut(&crate::variations::Sv),
    {
        let ref_base = self.effective_ref_base_i32(ctx.position).unwrap_or(b'N');
        let mut ref_buf = [0u8; 4];
        let ref_base_str: &str = (ref_base as char).encode_utf8(&mut ref_buf);

        for &description_string in keys {
            if description_string == "SV" {
                if let Some(ref sv) = non_insertion_variations.sv {
                    on_sv(sv);
                }
                continue;
            }
            if let Some(cnt) = non_insertion_variations.map.get(description_string) {
                if cnt.cnt == 0 {
                    continue;
                }
                let tvref = self.build_variant_from_cnt(
                    description_string,
                    cnt,
                    ctx.total_pos_coverage,
                    ctx.hicov,
                    ctx.duprate,
                    description_string == ref_base_str,
                );

                if self.config.debug {
                    debug_lines.push(
                        tvref.debug_variants_content_compat(
                            self.config.compat.exact_binary_half_even,
                        ),
                    );
                }

                var.push(tvref);
            }
        }
    }

    #[inline]
    fn adjust_insertion_coverage_and_next_ref(
        &self,
        orig_non_ins: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
        overlay: &mut rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
        ref_coverage: &rustc_hash::FxHashMap<i32, i32>,
        position: i32,
        cnt: &crate::variations::Variation,
        total_pos_coverage: &mut i32,
        mut ttcov: i32,
    ) -> i32 {
        if cnt.cnt > *total_pos_coverage
            && cnt.extracnt != 0
            && i64::from(cnt.cnt) - i64::from(*total_pos_coverage) < i64::from(cnt.extracnt)
        {
            ttcov = cnt.cnt;
        }

        if ttcov < cnt.cnt {
            ttcov = cnt.cnt;
            if let Some(next_pos) = position.checked_add(1) {
                if let Some(&cov_next) = ref_coverage.get(&next_pos) {
                    if ttcov < cov_next.saturating_sub(cnt.cnt) {
                        ttcov = cov_next;
                        if let Some(ref_base) = self.effective_ref_base_i32(next_pos) {
                            let mut buf = [0u8; 4];
                            let ref_base_str: &str = (ref_base as char).encode_utf8(&mut buf);
                            if let std::collections::hash_map::Entry::Vacant(e) =
                                overlay.entry(next_pos)
                            {
                                if let Some(orig) = orig_non_ins.get(&next_pos) {
                                    let mut cloned = orig.clone();
                                    if self.config.compat.sanitize_iupac_alt_alleles {
                                        merge_iupac_ref_alleles_in_place(
                                            &mut cloned,
                                            ref_base,
                                            ref_base_str,
                                        );
                                    }
                                    e.insert(cloned);
                                }
                            } else if self.config.compat.sanitize_iupac_alt_alleles {
                                if let Some(existing) = overlay.get_mut(&next_pos) {
                                    merge_iupac_ref_alleles_in_place(
                                        existing,
                                        ref_base,
                                        ref_base_str,
                                    );
                                }
                            }
                            if let Some(next_non_ins) = overlay.get_mut(&next_pos) {
                                if let Some(next_var) = next_non_ins.map.get_mut(ref_base_str) {
                                    if self.config.unique_mode_alignment_enabled
                                        && self.config.compat.mate_overlap_both_strands
                                    {
                                        let avail_f = i64::from(next_var.fwd.max(0));
                                        let avail_r = i64::from(next_var.rev.max(0));
                                        let sub_f = i64::from(cnt.fwd.max(0));
                                        let sub_r = i64::from(cnt.rev.max(0));
                                        if avail_f + avail_r >= sub_f + sub_r {
                                            let rem_f = (sub_f - avail_f).max(0);
                                            let rem_r = (sub_r - avail_r).max(0);
                                            next_var.fwd =
                                                ((avail_f - sub_f).max(0) - rem_r).max(0) as i32;
                                            next_var.rev =
                                                ((avail_r - sub_r).max(0) - rem_f).max(0) as i32;
                                        } else if self.config.compat.reconcile_strand_and_hq_sums {
                                            next_var.fwd = 0;
                                            next_var.rev = 0;
                                        } else {
                                            next_var.fwd = next_var.fwd.saturating_sub(cnt.fwd);
                                            next_var.rev = next_var.rev.saturating_sub(cnt.rev);
                                        }
                                    } else {
                                        next_var.fwd = next_var.fwd.saturating_sub(cnt.fwd);
                                        next_var.rev = next_var.rev.saturating_sub(cnt.rev);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            *total_pos_coverage = ttcov;
        }
        ttcov
    }

    fn apply_prev_insertion_overlay(
        &self,
        insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
        orig_non_ins: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
        overlay: &mut rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
        ref_coverage: &rustc_hash::FxHashMap<i32, i32>,
        position: i32,
        mut total_pos_coverage: i32,
    ) {
        let Some(insertion_variations) = insertion_variants.get(&position) else {
            return;
        };
        let insertion_description_strings =
            SortedKeySlice::new(insertion_variations.map.keys(), None);

        for &description_string in insertion_description_strings.as_slice() {
            if description_string
                .strip_prefix('+')
                .unwrap_or(description_string)
                .is_empty()
            {
                continue;
            }
            let mut ttcov = total_pos_coverage;
            if description_string.contains('&') {
                if let Some(&cov_next) = position
                    .checked_add(1)
                    .and_then(|next_pos| ref_coverage.get(&next_pos))
                {
                    total_pos_coverage = cov_next;
                    ttcov = cov_next;
                }
            }
            if let Some(cnt) = insertion_variations.map.get(description_string) {
                if cnt.cnt == 0 {
                    continue;
                }
                let _ = self.adjust_insertion_coverage_and_next_ref(
                    orig_non_ins,
                    overlay,
                    ref_coverage,
                    position,
                    cnt,
                    &mut total_pos_coverage,
                    ttcov,
                );
            }
        }
    }

    fn create_insertion(
        &self,
        insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
        orig_non_ins: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
        overlay: &mut rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
        ref_coverage: &rustc_hash::FxHashMap<i32, i32>,
        duprate: f64,
        position: i32,
        mut total_pos_coverage: i32,
        var: &mut Vec<Variant>,
        debug_lines: &mut Vec<String>,
        mut hicov: i32,
    ) -> i32 {
        if let Some(insertion_variations) = insertion_variants.get(&position) {
            let insertion_description_strings =
                SortedKeySlice::new(insertion_variations.map.keys(), None);

            for &description_string in insertion_description_strings.as_slice() {
                if description_string
                    .strip_prefix('+')
                    .unwrap_or(description_string)
                    .is_empty()
                {
                    continue;
                }
                let mut ttcov = total_pos_coverage;
                if description_string.contains('&') {
                    if let Some(&cov_next) = position
                        .checked_add(1)
                        .and_then(|next_pos| ref_coverage.get(&next_pos))
                    {
                        total_pos_coverage = cov_next;
                        ttcov = cov_next;
                    }
                }
                if let Some(cnt) = insertion_variations.map.get(description_string) {
                    if cnt.cnt == 0 {
                        continue;
                    }
                    let exact_binary_half_even = self.config.compat.exact_binary_half_even;
                    let reconcile_strand_hq = self.config.compat.reconcile_strand_and_hq_sums;
                    let reconcile_dp =
                        self.config.compat.reconcile_vd_dp || self.config.compat.reconcile_rd_vd_dp;
                    let (_vd, fwd, rev, hicnt, locnt) = if reconcile_strand_hq {
                        reconcile_strand_and_hq(cnt.cnt, cnt.fwd, cnt.rev, cnt.hicnt)
                    } else {
                        let hc = if reconcile_dp {
                            cnt.hicnt.max(0)
                        } else {
                            cnt.hicnt
                        };
                        (cnt.cnt, cnt.fwd, cnt.rev, hc, cnt.locnt)
                    };
                    let min_reads_bias =
                        i32::try_from(self.config.min_reads_bias).unwrap_or(i32::MAX);
                    let bias = strand_bias(fwd, rev, min_reads_bias);
                    let base_quality = round_half_even_compat(
                        cnt.qmean / cnt.cnt as f64,
                        1,
                        exact_binary_half_even,
                    );
                    let mapping_quality = round_half_even_compat(
                        cnt.mqmean / cnt.cnt as f64,
                        1,
                        exact_binary_half_even,
                    );

                    ttcov = self.adjust_insertion_coverage_and_next_ref(
                        orig_non_ins,
                        overlay,
                        ref_coverage,
                        position,
                        cnt,
                        &mut total_pos_coverage,
                        ttcov,
                    );

                    if hicov < hicnt {
                        hicov = hicnt;
                    }
                    if reconcile_dp || reconcile_strand_hq {
                        hicov = hicov.max(0);
                    }

                    let mut tvref = Variant::empty_internal();
                    fill_common_fields(
                        &mut tvref,
                        description_string,
                        cnt,
                        fwd,
                        rev,
                        bias,
                        ttcov,
                        base_quality,
                        mapping_quality,
                        hicnt,
                        locnt,
                        hicov,
                        duprate,
                        CommonFieldFlags {
                            exact_binary_half_even,
                            reconcile_dp,
                            reconcile_strand_hq,
                            nonnegative_subtracted_means: self
                                .config
                                .compat
                                .nonnegative_subtracted_means,
                        },
                    );

                    if self.config.debug {
                        // ports Variant.debugVariantsContentInsertion (Variant.java:261-266)
                        let mut s = String::from("I");
                        s.push_str(&tvref.debug_variants_content_compat(exact_binary_half_even));
                        debug_lines.push(s);
                    }

                    var.push(tvref);
                }
            }
        }
        total_pos_coverage
    }

    #[cfg(test)]
    fn collect_vars_at_position(
        &self,
        aligned_variants: &mut crate::javamap::JavaHashMap<i32, crate::variations::Vars>,
        position: i32,
        var: &mut Vec<Variant>,
    ) -> f64 {
        // H1: bound per-element Vec::push regrowth (dhat: ~8% on deep
        // pileup). Guarded by `!var.is_empty()` so an empty `var` creates
        // nothing, preserving exact JavaHashMap entry-creation timing.
        if var.is_empty() {
            return 0.0;
        }
        let entry = crate::variations::get_or_put_vars(aligned_variants, position);
        self.collect_vars_into(entry, position, var)
    }

    fn collect_vars_into(
        &self,
        entry: &mut crate::variations::Vars,
        position: i32,
        var: &mut Vec<Variant>,
    ) -> f64 {
        let ref_base = self.effective_ref_base_i32(position).unwrap_or(b'N');
        let mut buf = [0u8; 4];
        let ref_base_str: &str = (ref_base as char).encode_utf8(&mut buf);

        let non_ref_rows = var
            .iter()
            .filter(|t| t.description.as_str() != ref_base_str)
            .count();
        Self::drain_vars_into(entry, ref_base_str, non_ref_rows, var)
    }

    fn drain_vars_into(
        entry: &mut crate::variations::Vars,
        ref_base_str: &str,
        non_ref_rows: usize,
        var: &mut Vec<Variant>,
    ) -> f64 {
        let mut maxfreq = 0.0;
        if non_ref_rows > 0 {
            if entry.variants.is_empty() {
                entry.variants.reserve_exact(non_ref_rows);
            } else {
                entry.variants.reserve(non_ref_rows);
            }
        }

        for tvar in var.drain(..) {
            if tvar.description == ref_base_str {
                entry.ref_variant = Some(tvar);
            } else {
                if tvar.frequency > maxfreq {
                    maxfreq = tvar.frequency;
                }
                entry.variants.push(tvar);
            }
        }
        maxfreq
    }

    fn sort_variants(&self, var: &mut [Variant]) {
        var.sort_by(|o1, o2| {
            let metric1 = o1.mean_quality * o1.pos_coverage as f64;
            let metric2 = o2.mean_quality * o2.pos_coverage as f64;
            metric2
                .partial_cmp(&metric1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| o1.description.cmp(&o2.description))
        });
    }

    #[allow(clippy::too_many_arguments)]
    fn collect_reference_variants(
        &self,
        position: i32,
        total_pos_coverage: i32,
        variations_at_pos: &mut crate::variations::Vars,
        ref_coverage: &rustc_hash::FxHashMap<i32, i32>,
        insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
        orig_non_ins: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
        overlay: &rustc_hash::FxHashMap<i32, crate::variations::VariationMap>,
        duprate: f64,
        debug_lines: &[String],
        scratch1: &mut String,
        scratch2: &mut String,
    ) {
        if variations_at_pos.ref_variant.is_none() && variations_at_pos.variants.is_empty() {
            // Ports ToVarsBuilder.java:609 + :163: when both referenceVariant == null
            // and variants.isEmpty(), Java throws NullPointerException at line 609
            // (caught at line 163), aborting collectReferenceVariants after
            // getOrPutVars(alignedVariants, position) inserted the empty Vars.
            return;
        }

        let mut reference_forward_coverage = 0;
        let mut reference_reverse_coverage = 0;

        if let Some(ref ref_var) = variations_at_pos.ref_variant {
            reference_forward_coverage = ref_var.fwd;
            reference_reverse_coverage = ref_var.rev;
        }

        let mut reference_forward_coverage_adj = reference_forward_coverage;
        let mut reference_reverse_coverage_adj = reference_reverse_coverage;
        let mut total_pos_coverage_adj = total_pos_coverage;

        if total_pos_coverage_adj > ref_coverage.get(&position).copied().unwrap_or(0) {
            if let Some(next_pos) = position.checked_add(1) {
                if let Some(ref_base_next) = self.effective_ref_base_i32(next_pos) {
                    let mut buf = [0u8; 4];
                    let ref_base_next_str: &str = (ref_base_next as char).encode_utf8(&mut buf);
                    if let Some(next_ov) = overlay.get(&next_pos) {
                        if let Some(tpref) = next_ov.map.get(ref_base_next_str) {
                            reference_forward_coverage_adj = tpref.fwd;
                            reference_reverse_coverage_adj = tpref.rev;
                        }
                    } else if let Some(orig_next) = orig_non_ins.get(&next_pos) {
                        if self.config.compat.sanitize_iupac_alt_alleles {
                            let mut sum_f: i32 = 0;
                            let mut sum_r: i32 = 0;
                            let mut any = false;
                            for (k, v) in &orig_next.map {
                                if is_sw26_iupac_ref_key(k, ref_base_next) {
                                    sum_f = sum_f.saturating_add(v.fwd);
                                    sum_r = sum_r.saturating_add(v.rev);
                                    any = true;
                                }
                            }
                            if any {
                                reference_forward_coverage_adj = sum_f;
                                reference_reverse_coverage_adj = sum_r;
                            }
                        } else if let Some(tpref) = orig_next.map.get(ref_base_next_str) {
                            reference_forward_coverage_adj = tpref.fwd;
                            reference_reverse_coverage_adj = tpref.rev;
                        }
                    }
                }
            }
        }

        let mut changed_ref_variant_at_pos = false;

        if !variations_at_pos.variants.is_empty() {
            let standard_dup_symbolic_alt = self.config.compat.standard_dup_symbolic_alt;
            let canonical_sv_vcf_anchoring = self.config.compat.canonical_sv_vcf_anchoring;
            let sv_respect_min_len = self.config.compat.sv_respect_min_len;
            let sv_min_len = if sv_respect_min_len {
                self.config
                    .sv_min_len
                    .min(crate::sv::MAX_EXPLICIT_SV_SEQ_BP as u32)
            } else {
                self.config.sv_min_len
            };
            let mut gt1_buf = crate::types::StackBuf::<24>::new();
            let genotype1_current: std::borrow::Cow<'_, str> =
                if let Some(ref ref_var) = variations_at_pos.ref_variant {
                    if ref_var.frequency >= self.config.freq {
                        if ref_var.description.starts_with('+') {
                            format_plus_genotype1(
                                &ref_var.description,
                                &mut gt1_buf,
                                standard_dup_symbolic_alt,
                                canonical_sv_vcf_anchoring,
                                sv_respect_min_len,
                                sv_min_len,
                            )
                        } else {
                            std::borrow::Cow::Borrowed(ref_var.description.as_str())
                        }
                    } else if let Some(first_var) = variations_at_pos.variants.first() {
                        if first_var.description.starts_with('+') {
                            format_plus_genotype1(
                                &first_var.description,
                                &mut gt1_buf,
                                standard_dup_symbolic_alt,
                                canonical_sv_vcf_anchoring,
                                sv_respect_min_len,
                                sv_min_len,
                            )
                        } else {
                            std::borrow::Cow::Owned(first_var.description.clone())
                        }
                    } else if ref_var.description.starts_with('+') {
                        format_plus_genotype1(
                            &ref_var.description,
                            &mut gt1_buf,
                            standard_dup_symbolic_alt,
                            canonical_sv_vcf_anchoring,
                            sv_respect_min_len,
                            sv_min_len,
                        )
                    } else {
                        std::borrow::Cow::Borrowed(ref_var.description.as_str())
                    }
                } else if let Some(first_var) = variations_at_pos.variants.first() {
                    if first_var.description.starts_with('+') {
                        format_plus_genotype1(
                            &first_var.description,
                            &mut gt1_buf,
                            standard_dup_symbolic_alt,
                            canonical_sv_vcf_anchoring,
                            sv_respect_min_len,
                            sv_min_len,
                        )
                    } else {
                        std::borrow::Cow::Owned(first_var.description.clone())
                    }
                } else {
                    return;
                };

            let chr_len = self.chr_len;

            for vref in &mut variations_at_pos.variants {
                let mut genotype1_current_inner: std::borrow::Cow<'_, str> =
                    std::borrow::Cow::Borrowed(genotype1_current.as_ref());
                let mut gt2_buf = crate::types::StackBuf::<24>::new();
                let mut genotype2: std::borrow::Cow<'_, str> = if vref.description.starts_with('+')
                {
                    format_plus_genotype(
                        &mut gt2_buf,
                        (vref.description.len() as i32).saturating_sub(1),
                    )
                } else {
                    std::borrow::Cow::Borrowed(vref.description.as_str())
                };

                let description_string = &vref.description;
                let deletion_length = get_begin_minus_number(description_string).unwrap_or(0);

                let mut end_position = position;
                if description_string.starts_with('-') {
                    end_position = (i64::from(position) + i64::from(deletion_length) - 1)
                        .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                        as i32;
                }

                let mut refallele = String::new();
                let mut varallele;

                let mut shift3 = 0;
                let mut msi = 0.0;
                let mut msint_len = 0;

                let mut start_position = position;
                let mut canonical_symbolic_sv_locked = false;
                let mut canonical_symbolic_del = false;

                if let Some(stripped_desc) = description_string.strip_prefix('+') {
                    if !description_string.contains('&')
                        && !description_string.contains('#')
                        && !description_string.contains("<dup")
                    {
                        (msi, shift3, msint_len) = self.proceed_vref_is_insertion(
                            position,
                            description_string,
                            scratch1,
                            scratch2,
                        );
                    }
                    let shift_ins_alleles = self.config.three_prime
                        && self.config.compat.shift3_3prime_alleles
                        && shift3 > 0;
                    if self.config.three_prime {
                        start_position = start_position.saturating_add(shift3);
                        end_position = end_position.saturating_add(shift3);
                    }
                    let ins_anchor_pos = if shift_ins_alleles {
                        start_position
                    } else {
                        position
                    };
                    let ref_base_byte = self.effective_ref_base_i32(ins_anchor_pos).unwrap_or(b'N');
                    refallele.push(ref_base_byte as char);
                    varallele = String::with_capacity(1 + stripped_desc.len());
                    varallele.push(ref_base_byte as char);
                    if shift_ins_alleles
                        && !stripped_desc.is_empty()
                        && !stripped_desc.contains('<')
                        && stripped_desc.is_ascii()
                    {
                        let rot = (shift3 as usize) % stripped_desc.len();
                        if let (Some(tail), Some(head)) =
                            (stripped_desc.get(rot..), stripped_desc.get(..rot))
                        {
                            varallele.push_str(tail);
                            varallele.push_str(head);
                        } else {
                            varallele.push_str(stripped_desc);
                        }
                    } else {
                        varallele.push_str(stripped_desc);
                    }

                    let parse_neg_dup = standard_dup_symbolic_alt
                        || canonical_sv_vcf_anchoring
                        || sv_respect_min_len;
                    if parse_neg_dup {
                        if let Some(dup_val) = get_dup_num_compat(&varallele, parse_neg_dup) {
                            let dup_span = ((2 * SV_FLANK_BP) + i64::from(dup_val)).max(1);
                            if canonical_sv_vcf_anchoring {
                                // clintval/vartovcf#192 & clintval/vartovcf#107:
                                // - Origin A (pre_anchored_dup == true, already at bp - 1):
                                //   keep start_position and set end_position = start_position + dup_span.
                                // - Origin B (pre_anchored_dup == false, unshifted at bp):
                                //   shift start_position to bp - 1 iff loaded & non-N, setting
                                //   end_position = start_position + dup_span (= bp + dup_span - 1).
                                let pre_anchored = lookup_pre_anchored_dup(
                                    position,
                                    description_string,
                                    insertion_variants,
                                    orig_non_ins,
                                    overlay,
                                );
                                if pre_anchored {
                                    end_position = (i64::from(start_position) + dup_span)
                                        .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                                        as i32;
                                } else if let Some(b_prev) =
                                    self.upstream_non_n_anchor_base(start_position)
                                {
                                    start_position = start_position.saturating_sub(1);
                                    refallele.clear();
                                    refallele.push(b_prev as char);
                                    end_position = (i64::from(start_position) + dup_span)
                                        .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                                        as i32;
                                } else {
                                    end_position = (i64::from(start_position) + dup_span - 1)
                                        .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                                        as i32;
                                }
                                canonical_symbolic_sv_locked = true;
                            } else {
                                end_position = (i64::from(start_position) + dup_span - 1)
                                    .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                                    as i32;
                            }
                            genotype2 = format_plus_genotype(&mut gt2_buf, dup_span);
                            varallele = "<DUP>".to_string();
                        } else {
                            let exceeds_sv_min = if sv_respect_min_len {
                                let ins_len =
                                    long_ins_inserted_bases(stripped_desc, description_string);
                                u32::try_from(ins_len).unwrap_or(u32::MAX) >= sv_min_len
                            } else {
                                varallele.len() > sv_min_len as usize
                            };
                            if exceeds_sv_min {
                                if canonical_sv_vcf_anchoring {
                                    // Origin C (clintval/vartovcf#107, AstraZeneca-NGS/VarDictJava#398,
                                    // ToVarsBuilder.java:692-694): varallele carries the leading 1-bp
                                    // reference base, so use exact inserted base count ins_len so
                                    // END - POS == ins_len == SVLEN and clamp refallele to 1 bp.
                                    let anchor_base =
                                        self.effective_ref_base_i32(start_position).unwrap_or(b'N');
                                    refallele.clear();
                                    refallele.push(anchor_base as char);
                                    let ins_len =
                                        long_ins_inserted_bases(stripped_desc, description_string);
                                    end_position = start_position.saturating_add(ins_len);
                                    genotype2 =
                                        format_plus_genotype(&mut gt2_buf, i64::from(ins_len));
                                    canonical_symbolic_sv_locked = true;
                                } else {
                                    end_position = end_position.saturating_add(
                                        varallele.len().min(i32::MAX as usize) as i32,
                                    );
                                }
                                varallele = "<DUP>".to_string();
                            }
                        }
                    } else {
                        if varallele.len() > sv_min_len as usize {
                            end_position = end_position
                                .saturating_add(varallele.len().min(i32::MAX as usize) as i32);
                            varallele = "<DUP>".to_string();
                        }
                        if let Some(dup_val) = get_dup_num(&varallele) {
                            let dup_span = (2 * SV_FLANK_BP) + i64::from(dup_val);
                            end_position = (i64::from(start_position) + dup_span - 1)
                                .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                                as i32;
                            genotype2 = format_plus_genotype(&mut gt2_buf, dup_span);
                            varallele = "<DUP>".to_string();
                        }
                    }
                } else if description_string.starts_with('-') {
                    let can_materialize_del = i64::from(deletion_length) < i64::from(sv_min_len)
                        && (!self.config.compat.sv_respect_min_len
                            || self.is_explicit_del_span_loaded(position, end_position));
                    if can_materialize_del {
                        // replace first -digits with empty string
                        varallele = if let Some(rest) = description_string.strip_prefix('-') {
                            let digit_len = rest.bytes().take_while(|b| b.is_ascii_digit()).count();
                            if digit_len > 0 {
                                rest.get(digit_len..).unwrap_or("").to_string()
                            } else {
                                description_string.clone()
                            }
                        } else {
                            description_string.clone()
                        };

                        (msi, shift3, msint_len) = self.proceed_vref_is_deletion(
                            position,
                            deletion_length,
                            scratch1,
                            scratch2,
                        );

                        let parse_neg_inv = canonical_sv_vcf_anchoring || sv_respect_min_len;
                        if get_inv_num_compat(description_string, parse_neg_inv).is_some() {
                            varallele = "<INV>".to_string();
                            genotype2 =
                                std::borrow::Cow::Owned(format!("<INV{}>", deletion_length));
                        }
                    } else if matches_begin_minus_number_caret(description_string)
                        || ((canonical_sv_vcf_anchoring || sv_respect_min_len)
                            && get_inv_num_compat(description_string, true).is_some())
                    {
                        varallele = "<INV>".to_string();
                        genotype2 = std::borrow::Cow::Owned(format!("<INV{}>", deletion_length));
                    } else {
                        varallele = "<DEL>".to_string();
                    }

                    let simple_del = !description_string.contains('&')
                        && !description_string.contains('#')
                        && !description_string.contains('^');
                    let shift_del_alleles = simple_del
                        && self.config.three_prime
                        && self.config.compat.shift3_3prime_alleles
                        && shift3 > 0;
                    if canonical_sv_vcf_anchoring && varallele == "<INV>" {
                        // clintval/vartovcf#192: anchor <INV> (both -len^... and -len<invN>)
                        // at del_seq_pos - 1 when loaded & non-N so END - POS == deletion_length.
                        let del_seq_pos = if shift_del_alleles {
                            position.saturating_add(shift3)
                        } else {
                            position
                        };
                        end_position = (i64::from(del_seq_pos) + i64::from(deletion_length) - 1)
                            .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                            as i32;
                        if let Some(b_prev) = self.upstream_non_n_anchor_base(del_seq_pos) {
                            start_position = del_seq_pos - 1;
                            refallele.push(b_prev as char);
                        } else {
                            start_position = del_seq_pos;
                            let base_cur = self.effective_ref_base_i32(del_seq_pos).unwrap_or(b'N');
                            refallele.push(base_cur as char);
                        }
                        canonical_symbolic_sv_locked = true;
                    } else if canonical_sv_vcf_anchoring && varallele == "<DEL>" {
                        // AstraZeneca-NGS/VarDictJava#398, AstraZeneca-NGS/VarDictJava#225,
                        // AstraZeneca-NGS/VarDictJava#229, AstraZeneca-NGS/VarDict#111,
                        // AstraZeneca-NGS/VarDictJava#382, clintval/vartovcf#192:
                        // anchor both simple and complex (&/#/^) symbolic <DEL> at position - 1
                        // when loaded & non-N (else unshifted >= 1), clamping refallele to 1 bp.
                        if let Some(b_prev) = self.upstream_non_n_anchor_base(position) {
                            start_position = position - 1;
                            refallele.push(b_prev as char);
                        } else {
                            start_position = position;
                            let base_cur = self.effective_ref_base_i32(position).unwrap_or(b'N');
                            refallele.push(base_cur as char);
                        }
                        canonical_symbolic_del = true;
                    } else {
                        if simple_del {
                            if self.config.three_prime {
                                start_position = start_position.saturating_add(shift3);
                                if shift_del_alleles {
                                    end_position = end_position.saturating_add(shift3);
                                }
                            }
                            let anchor_lookup_pos = if shift_del_alleles {
                                position.saturating_add(shift3)
                            } else {
                                position
                            };
                            let base_prev = anchor_lookup_pos
                                .checked_sub(1)
                                .and_then(|p| self.effective_ref_base_i32(p))
                                .unwrap_or(b'N');
                            if varallele != "<DEL>" {
                                varallele.clear();
                                varallele.push(base_prev as char);
                            }
                            refallele.push(base_prev as char);
                            start_position = start_position.saturating_sub(1);
                        }

                        let del_seq_pos = if shift_del_alleles {
                            position.saturating_add(shift3)
                        } else {
                            position
                        };
                        if matches_some_sv_numbers_compat(description_string, sv_respect_min_len) {
                            let base_cur = self.effective_ref_base_i32(del_seq_pos).unwrap_or(b'N');
                            refallele.clear();
                            refallele.push(base_cur as char);
                        } else if can_materialize_del {
                            let del_end = (i64::from(del_seq_pos) + i64::from(deletion_length) - 1)
                                .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                                as i32;
                            if !self.config.compat.sv_respect_min_len
                                || self.is_explicit_del_span_loaded(del_seq_pos, del_end)
                            {
                                self.join_ref_into(scratch1, del_seq_pos, del_end);
                                refallele.push_str(scratch1);
                            } else {
                                varallele = "<DEL>".to_string();
                            }
                        }
                    }
                } else {
                    self.join_ref_into(
                        scratch1,
                        std::cmp::max(position.saturating_sub(30), 1),
                        position.saturating_add(1),
                    );
                    if let Some(tseq2_start) = position.checked_add(2) {
                        self.join_ref_into(
                            scratch2,
                            tseq2_start,
                            std::cmp::min(position.saturating_add(70), chr_len),
                        );
                    } else {
                        scratch2.clear();
                    }

                    (msi, shift3, msint_len) =
                        find_msi_metrics(scratch1.as_bytes(), scratch2.as_bytes(), &[]);

                    let base_cur = self.effective_ref_base_i32(position).unwrap_or(b'N');
                    refallele.push(base_cur as char);
                    varallele = description_string.clone();
                }

                if let Some(extra) = get_amp_atgc(description_string) {
                    if canonical_symbolic_sv_locked {
                        // Canonical symbolic <DUP> / <INV> already fixed start_position,
                        // end_position, 1-bp refallele, and varallele.
                    } else if canonical_symbolic_del {
                        let extra_len = extra.len().min(i32::MAX as usize) as i32;
                        if let Some(ext_start) = end_position.checked_add(1) {
                            self.join_ref_into(
                                scratch1,
                                ext_start,
                                end_position.saturating_add(extra_len),
                            );
                        } else {
                            scratch1.clear();
                        }
                        if !scratch1.is_empty() {
                            genotype1_current_inner.to_mut().push_str(scratch1);
                        }
                        end_position = end_position.saturating_add(extra_len);
                        if let Some(&cov_prev) =
                            position.checked_sub(1).and_then(|p| ref_coverage.get(&p))
                        {
                            total_pos_coverage_adj = cov_prev;
                        }
                        if vref.pos_coverage > total_pos_coverage_adj {
                            total_pos_coverage_adj = vref.pos_coverage;
                        }
                        if total_pos_coverage_adj > 0 {
                            vref.frequency =
                                vref.pos_coverage as f64 / total_pos_coverage_adj as f64;
                        }
                    } else {
                        let extra_len = extra.len().min(i32::MAX as usize) as i32;
                        if let Some(idx) = varallele.find('&') {
                            varallele.remove(idx);
                        }
                        if let Some(ext_start) = end_position.checked_add(1) {
                            self.join_ref_into(
                                scratch1,
                                ext_start,
                                end_position.saturating_add(extra_len),
                            );
                        } else {
                            scratch1.clear();
                        }
                        refallele.push_str(scratch1);
                        if !scratch1.is_empty() {
                            genotype1_current_inner.to_mut().push_str(scratch1);
                        }
                        end_position = end_position.saturating_add(extra_len);

                        if let Some(vextra) = get_amp_atgc(&varallele) {
                            let vextra_len = vextra.len().min(i32::MAX as usize) as i32;
                            if let Some(idx) = varallele.find('&') {
                                varallele.remove(idx);
                            }
                            if let Some(vext_start) = end_position.checked_add(1) {
                                self.join_ref_into(
                                    scratch1,
                                    vext_start,
                                    end_position.saturating_add(vextra_len),
                                );
                            } else {
                                scratch1.clear();
                            }
                            refallele.push_str(scratch1);
                            if !scratch1.is_empty() {
                                genotype1_current_inner.to_mut().push_str(scratch1);
                            }
                            end_position = end_position.saturating_add(vextra_len);
                        }

                        if description_string.starts_with('+') {
                            if let Some(ch) = refallele.chars().next() {
                                refallele.drain(..ch.len_utf8());
                            }
                            if let Some(ch) = varallele.chars().next() {
                                varallele.drain(..ch.len_utf8());
                            }
                            start_position = start_position.saturating_add(1);
                        }

                        if varallele == "<DEL>" && !refallele.is_empty() {
                            let ref_base =
                                self.effective_ref_base_i32(start_position).unwrap_or(b'N');
                            refallele.clear();
                            refallele.push(ref_base as char);
                            if let Some(&cov_prev) = start_position
                                .checked_sub(1)
                                .and_then(|p| ref_coverage.get(&p))
                            {
                                total_pos_coverage_adj = cov_prev;
                            }
                            if vref.pos_coverage > total_pos_coverage_adj {
                                total_pos_coverage_adj = vref.pos_coverage;
                            }
                            vref.frequency =
                                vref.pos_coverage as f64 / total_pos_coverage_adj as f64;
                        }
                    }
                }

                if let Some((matched_sequence, tail)) =
                    get_hash_group_caret_group(description_string)
                {
                    if canonical_symbolic_sv_locked {
                        // Keep canonical <DUP> / <INV> coordinates and 1-bp refallele intact.
                    } else if canonical_symbolic_del {
                        let matched_len = matched_sequence.len().min(i32::MAX as usize) as i32;
                        end_position = end_position.saturating_add(matched_len);
                        if let Some(deletion) = get_begin_digits(tail) {
                            end_position = end_position.saturating_add(deletion);
                        }
                    } else {
                        let matched_len = matched_sequence.len().min(i32::MAX as usize) as i32;
                        let prev_end = end_position;
                        end_position = end_position.saturating_add(matched_len);
                        if let Some(m_start) = prev_end.checked_add(1) {
                            self.join_ref_into(scratch1, m_start, end_position);
                        } else {
                            scratch1.clear();
                        }
                        refallele.push_str(scratch1);

                        if let Some(deletion) = get_begin_digits(tail) {
                            if let Some(del_start) = end_position.checked_add(1) {
                                self.join_ref_into(
                                    scratch1,
                                    del_start,
                                    end_position.saturating_add(deletion),
                                );
                            } else {
                                scratch1.clear();
                            }
                            refallele.push_str(scratch1);
                            end_position = end_position.saturating_add(deletion);
                        }

                        // replaceFirst("#", "").replaceFirst("\\^(\\d+)?", "") in Java:
                        let cleaned_varallele = varallele.replacen('#', "", 1);
                        let cleaned_varallele =
                            if let Some((prefix, rest)) = cleaned_varallele.split_once('^') {
                                let tail = rest.trim_start_matches(|c: char| c.is_ascii_digit());
                                format!("{prefix}{tail}")
                            } else {
                                cleaned_varallele
                            };
                        varallele = cleaned_varallele;
                    }

                    if genotype1_current_inner.contains('#') {
                        genotype1_current_inner =
                            std::borrow::Cow::Owned(genotype1_current_inner.replacen('#', "m", 1));
                    }
                    if genotype1_current_inner.contains('^') {
                        genotype1_current_inner =
                            std::borrow::Cow::Owned(genotype1_current_inner.replacen('^', "i", 1));
                    }
                    if genotype2.contains('#') {
                        genotype2 = std::borrow::Cow::Owned(genotype2.replacen('#', "m", 1));
                    }
                    if genotype2.contains('^') {
                        genotype2 = std::borrow::Cow::Owned(genotype2.replacen('^', "i", 1));
                    }
                }

                if get_caret_atgnc(description_string).is_some() {
                    if !canonical_symbolic_sv_locked && !canonical_symbolic_del {
                        varallele = varallele.replacen('^', "", 1);
                    }
                    if genotype1_current_inner.contains('^') {
                        genotype1_current_inner =
                            std::borrow::Cow::Owned(genotype1_current_inner.replacen('^', "i", 1));
                    }
                    if genotype2.contains('^') {
                        genotype2 = std::borrow::Cow::Owned(genotype2.replacen('^', "i", 1));
                    }
                }

                // CRISPR Site adjustment
                if let Some(cut_site) = self.config.crispr_cutting_site {
                    if cut_site != 0 && refallele.len() > 1 && varallele.len() > 1 {
                        let mut n = 0usize;
                        while refallele.len() > n.saturating_add(1)
                            && varallele.len() > n.saturating_add(1)
                            && refallele.as_bytes().get(n) == varallele.as_bytes().get(n)
                        {
                            n += 1;
                        }
                        if n != 0 {
                            start_position =
                                start_position.saturating_add(n.min(i32::MAX as usize) as i32);
                            refallele = refallele.get(n..).unwrap_or("").to_string();
                            varallele = varallele.get(n..).unwrap_or("").to_string();
                        }
                    }

                    if cut_site != 0
                        && refallele.len() != varallele.len()
                        && refallele.starts_with(varallele.chars().next().unwrap_or('\0'))
                        && !(start_position == cut_site || end_position == cut_site)
                    {
                        let cut = i64::from(cut_site);
                        let start64 = i64::from(start_position);
                        let end64 = i64::from(end_position);
                        let already_shifted_3prime = self.config.three_prime
                            && self.config.compat.shift3_3prime_alleles
                            && shift3 > 0
                            && (description_string.starts_with('+')
                                || (!description_string.contains('&')
                                    && !description_string.contains('#')
                                    && !description_string.contains('^')));
                        let shift3_64 = if already_shifted_3prime {
                            0
                        } else {
                            i64::from(shift3)
                        };
                        let mut n64 = 0_i64;
                        let dis = std::cmp::min((cut - start64).abs(), (cut - end64).abs());
                        if start64 < cut {
                            while start64 + n64 < cut && n64 < shift3_64 && end64 + n64 != cut {
                                n64 += 1;
                            }
                            if (start64 + n64 - cut).abs() > dis && (end64 + n64 - cut).abs() > dis
                            {
                                n64 = 0;
                            }
                        }
                        if end64 < cut && n64 == 0 && (end64 - cut).abs() <= (start64 - cut).abs() {
                            while end64 + n64 < cut && n64 < shift3_64 {
                                n64 += 1;
                            }
                        }
                        let n = n64.clamp(0, i64::from(i32::MAX)) as i32;
                        if n > 0 {
                            start_position = start_position.saturating_add(n);
                            end_position = end_position.saturating_add(n);
                            self.join_ref_into(&mut refallele, start_position, end_position);
                            let mut tva = String::new();
                            if refallele.len() < varallele.len() {
                                tva = varallele.get(1..).unwrap_or("").to_string();
                                if tva.len() > 1 {
                                    let ttn = (n as usize) % tva.len();
                                    if ttn != 0 {
                                        if let (Some(suffix), Some(prefix)) =
                                            (tva.get(ttn..), tva.get(..ttn))
                                        {
                                            tva = format!("{suffix}{prefix}");
                                        }
                                    }
                                }
                            }
                            let first_char =
                                self.effective_ref_base_i32(start_position).unwrap_or(b'N') as char;
                            varallele = format!("{}{}", first_char, tva);
                        }
                        vref.crispr = n;
                    }
                }

                self.join_ref_into(
                    &mut vref.left_seq,
                    std::cmp::max(start_position.saturating_sub(20), 1),
                    start_position.saturating_sub(1),
                );
                if let Some(right_start) = end_position.checked_add(1) {
                    self.join_ref_into(
                        &mut vref.right_seq,
                        right_start,
                        std::cmp::min(end_position.saturating_add(20), chr_len),
                    );
                } else {
                    vref.right_seq.clear();
                }

                let genotype = {
                    let mut gt =
                        String::with_capacity(genotype1_current_inner.len() + 1 + genotype2.len());
                    if genotype1_current_inner
                        .bytes()
                        .any(|b| matches!(b, b'&' | b'#' | b'^'))
                    {
                        for ch in genotype1_current_inner.chars() {
                            match ch {
                                '&' | '#' => {}
                                '^' => gt.push('i'),
                                _ => gt.push(ch),
                            }
                        }
                    } else {
                        gt.push_str(&genotype1_current_inner);
                    }
                    gt.push('/');
                    if genotype2.bytes().any(|b| matches!(b, b'&' | b'#' | b'^')) {
                        for ch in genotype2.chars() {
                            match ch {
                                '&' | '#' => {}
                                '^' => gt.push('i'),
                                _ => gt.push(ch),
                            }
                        }
                    } else {
                        gt.push_str(&genotype2);
                    }
                    gt
                };
                drop(genotype1_current_inner);
                drop(genotype2);

                let exact_binary_half_even = self.config.compat.exact_binary_half_even;
                if self.config.compat.reconcile_vd_dp || self.config.compat.reconcile_rd_vd_dp {
                    if vref.pos_coverage < 0 {
                        vref.pos_coverage = 0;
                    }
                    if vref.hicnt < 0 {
                        vref.hicnt = 0;
                    }
                    let mut final_dp = total_pos_coverage_adj
                        .max(vref.total_pos_cov)
                        .max(vref.pos_coverage)
                        .max(0);
                    if self.config.compat.reconcile_rd_vd_dp {
                        let rd_vd_sum = reference_forward_coverage_adj
                            .max(0)
                            .saturating_add(reference_reverse_coverage_adj.max(0))
                            .saturating_add(vref.pos_coverage.max(0));
                        final_dp = final_dp.max(rd_vd_sum);
                    }
                    if final_dp > vref.total_pos_cov && final_dp > 0 {
                        let raw_extracnt = lookup_raw_extracnt(
                            position,
                            vref.description.as_str(),
                            insertion_variants,
                            orig_non_ins,
                            overlay,
                        );
                        vref.frequency = round_half_even_compat(
                            vref.pos_coverage as f64 / final_dp as f64,
                            4,
                            exact_binary_half_even,
                        );
                        vref.extra_freq = if raw_extracnt > 0 {
                            round_half_even_compat(
                                raw_extracnt.min(final_dp) as f64 / final_dp as f64,
                                4,
                                exact_binary_half_even,
                            )
                        } else {
                            0.0
                        };
                    } else if final_dp > 0 {
                        vref.frequency = round_half_even_compat(
                            vref.pos_coverage as f64 / final_dp as f64,
                            4,
                            exact_binary_half_even,
                        );
                        vref.extra_freq =
                            round_half_even_compat(vref.extra_freq, 4, exact_binary_half_even);
                    } else {
                        vref.frequency = 0.0;
                        vref.extra_freq = 0.0;
                    }
                    vref.hicov = vref.hicov.max(vref.hicnt);
                    vref.hi_freq = if vref.hicov > 0 {
                        round_half_even_compat(
                            vref.hicnt as f64 / vref.hicov as f64,
                            4,
                            exact_binary_half_even,
                        )
                    } else {
                        0.0
                    };
                    vref.total_pos_cov = final_dp;
                } else {
                    vref.extra_freq =
                        round_half_even_compat(vref.extra_freq, 4, exact_binary_half_even);
                    vref.frequency =
                        round_half_even_compat(vref.frequency, 4, exact_binary_half_even);
                    vref.hi_freq = round_half_even_compat(vref.hi_freq, 4, exact_binary_half_even);
                    vref.total_pos_cov = total_pos_coverage_adj;
                }
                validate_ref_allele_compat_in_place(
                    &mut refallele,
                    self.config.compat.sanitize_iupac_alt_alleles,
                );
                if self.config.compat.sanitize_iupac_alt_alleles {
                    validate_alt_allele_in_place(&mut varallele);
                }
                if self.config.compat.indel_only_shift3
                    && !varallele.starts_with('<')
                    && refallele.len() == varallele.len()
                {
                    shift3 = 0;
                }
                vref.msi = round_half_even_compat(msi, 3, exact_binary_half_even);
                vref.msint = msint_len;
                vref.shift3 = shift3;
                vref.start_pos = start_position;
                vref.end_pos = end_position;
                vref.ref_allele = refallele;
                vref.var_allele = varallele;
                vref.genotype = Some(genotype);
                vref.ref_fwd_cov = reference_forward_coverage_adj;
                vref.ref_rev_cov = reference_reverse_coverage_adj;

                let ref_bias = variations_at_pos
                    .ref_variant
                    .as_ref()
                    .map_or("0", |ref_var| ref_var.bias.as_str());
                vref.bias.insert(0, ';');
                vref.bias.insert_str(0, ref_bias);

                vref.var_type = if vref.var_allele.starts_with('+')
                    || vref.description.starts_with('+')
                {
                    "Insertion".to_string()
                } else if vref.var_allele.starts_with('-') || vref.description.starts_with('-') {
                    "Deletion".to_string()
                } else if vref.description.contains('&') {
                    "Complex".to_string()
                } else if vref.description.len() > 1 {
                    "MNV".to_string()
                } else {
                    "SNV".to_string()
                };

                self.adjust_variant_counts(position, vref);

                if start_position != position && self.config.pileup {
                    changed_ref_variant_at_pos = true;
                }

                // ports ToVarsBuilder.java:947 (constructDebugLines inside the
                // non-reference-variant loop — every variant at the position shares
                // the same joined debug string)
                construct_debug_lines(debug_lines, vref, self.config.debug);
            }

            if self.config.no_sv {
                variations_at_pos
                    .variants
                    .retain(|vref| !matches_any_sv(&vref.var_allele));
            }
        } else if let Some(ref mut vref) = variations_at_pos.ref_variant {
            let eff_ref_base = self.effective_ref_base_i32(position);
            let raw_extracnt =
                if self.config.compat.reconcile_vd_dp || self.config.compat.reconcile_rd_vd_dp {
                    lookup_raw_ref_extracnt(
                        position,
                        vref.description.as_str(),
                        eff_ref_base,
                        self.config.compat.sanitize_iupac_alt_alleles,
                        insertion_variants,
                        orig_non_ins,
                        overlay,
                    )
                } else {
                    0
                };
            self.update_ref_variant(
                position,
                total_pos_coverage_adj,
                vref,
                RefSiteCounts {
                    fwd: reference_forward_coverage_adj,
                    rev: reference_reverse_coverage_adj,
                    extracnt: raw_extracnt,
                },
                duprate,
                debug_lines,
            );
        } else {
            variations_at_pos.ref_variant = Some(Variant::new());
        }

        if self.config.pileup && self.config.compat.pileup_preserve_shifted_ref {
            variations_at_pos
                .ref_variant
                .get_or_insert_with(Variant::new);
        }

        if variations_at_pos.ref_variant.is_some()
            && self.config.pileup
            && (changed_ref_variant_at_pos
                || self.config.amplicon.is_some()
                || self.config.compat.pileup_preserve_shifted_ref)
        {
            let eff_ref_base = self.effective_ref_base_i32(position);
            if let Some(ref mut vref) = variations_at_pos.ref_variant {
                let raw_extracnt = if self.config.compat.reconcile_vd_dp
                    || self.config.compat.reconcile_rd_vd_dp
                {
                    lookup_raw_ref_extracnt(
                        position,
                        vref.description.as_str(),
                        eff_ref_base,
                        self.config.compat.sanitize_iupac_alt_alleles,
                        insertion_variants,
                        orig_non_ins,
                        overlay,
                    )
                } else {
                    0
                };
                self.update_ref_variant(
                    position,
                    total_pos_coverage_adj,
                    vref,
                    RefSiteCounts {
                        fwd: reference_forward_coverage_adj,
                        rev: reference_reverse_coverage_adj,
                        extracnt: raw_extracnt,
                    },
                    duprate,
                    debug_lines,
                );
            }
        }
    }

    fn update_ref_variant(
        &self,
        position: i32,
        total_pos_coverage: i32,
        vref: &mut Variant,
        ref_counts: RefSiteCounts,
        duprate: f64,
        debug_lines: &[String],
    ) {
        let exact_binary_half_even = self.config.compat.exact_binary_half_even;
        if self.config.compat.reconcile_vd_dp || self.config.compat.reconcile_rd_vd_dp {
            let fwd = ref_counts.fwd.max(0);
            let rev = ref_counts.rev.max(0);
            let ref_cov_eff = fwd.saturating_add(rev);
            let final_dp = total_pos_coverage.max(vref.total_pos_cov).max(ref_cov_eff);
            vref.total_pos_cov = final_dp;
            vref.pos_coverage = 0;
            vref.frequency = 0.0;
            vref.ref_fwd_cov = fwd;
            vref.ref_rev_cov = rev;
            vref.extra_freq = if ref_counts.extracnt > 0 && final_dp > 0 {
                round_half_even_compat(
                    ref_counts.extracnt.min(final_dp) as f64 / final_dp as f64,
                    4,
                    exact_binary_half_even,
                )
            } else {
                0.0
            };
            vref.hicnt = vref.hicnt.max(0);
            vref.hicov = vref.hicov.max(vref.hicnt);
            vref.hi_freq = if vref.hicov > 0 {
                round_half_even_compat(
                    vref.hicnt as f64 / vref.hicov as f64,
                    4,
                    exact_binary_half_even,
                )
            } else {
                0.0
            };
        } else {
            vref.total_pos_cov = total_pos_coverage;
            vref.pos_coverage = 0;
            vref.frequency = 0.0;
            vref.ref_fwd_cov = ref_counts.fwd;
            vref.ref_rev_cov = ref_counts.rev;
            vref.hi_freq = round_half_even_compat(vref.hi_freq, 4, exact_binary_half_even);
        }
        vref.fwd = 0;
        vref.rev = 0;
        vref.msi = 0.0;
        vref.msint = 0;
        if !vref.bias.contains(';') {
            vref.bias.push_str(";0");
        }
        vref.shift3 = 0;
        vref.start_pos = position;
        vref.end_pos = position;

        let ref_base = self.get_ref_base_i32(position).unwrap_or(b'N');
        let (valid_base, base_ch) = if self.config.compat.sanitize_iupac_alt_alleles {
            let reduced = reduce_ref_byte_sw26(ref_base) as char;
            (reduced, reduced)
        } else {
            (reduce_iupac_byte(ref_base) as char, ref_base as char)
        };

        vref.ref_allele.clear();
        vref.ref_allele.push(valid_base);
        vref.var_allele.clear();
        vref.var_allele.push(valid_base);
        let mut gt = String::with_capacity(3);
        gt.push(base_ch);
        gt.push('/');
        gt.push(base_ch);
        vref.genotype = Some(gt);
        vref.left_seq.clear();
        vref.right_seq.clear();
        vref.duprate = duprate;

        construct_debug_lines(debug_lines, vref, self.config.debug);
    }

    fn adjust_variant_counts(&self, _position: i32, vref: &mut Variant) {
        if vref.ref_fwd_cov < 0 {
            vref.ref_fwd_cov = 0;
        }
        if vref.ref_rev_cov < 0 {
            vref.ref_rev_cov = 0;
        }
        if vref.fwd < 0 {
            vref.fwd = 0;
        }
        if vref.rev < 0 {
            vref.rev = 0;
        }
        if self.config.compat.reconcile_strand_and_hq_sums {
            let (vd, fwd, rev, hicnt, locnt) =
                reconcile_strand_and_hq(vref.pos_coverage, vref.fwd, vref.rev, vref.hicnt);
            vref.pos_coverage = vd;
            vref.fwd = fwd;
            vref.rev = rev;
            vref.hicnt = hicnt;
            vref.hq_lq_ratio = hicnt as f64 / if locnt != 0 { locnt as f64 } else { 0.5 };
            vref.hicov = vref.hicov.max(hicnt).max(0);
            vref.hi_freq = if vref.hicov > 0 && hicnt > 0 {
                round_half_even_compat(
                    hicnt as f64 / vref.hicov as f64,
                    4,
                    self.config.compat.exact_binary_half_even,
                )
            } else {
                0.0
            };
            let min_reads_bias = i32::try_from(self.config.min_reads_bias).unwrap_or(i32::MAX);
            let alt_bias = strand_bias(fwd, rev, min_reads_bias);
            if let Some(semi_idx) = vref.bias.find(';') {
                vref.bias.truncate(semi_idx + 1);
            } else {
                vref.bias.clear();
            }
            let _ = std::fmt::Write::write_fmt(&mut vref.bias, format_args!("{alt_bias}"));
        }
        if self.config.compat.reconcile_rd_vd_dp {
            let rd_vd_sum = vref
                .ref_fwd_cov
                .max(0)
                .saturating_add(vref.ref_rev_cov.max(0))
                .saturating_add(vref.pos_coverage.max(0));
            vref.total_pos_cov = vref.total_pos_cov.max(rd_vd_sum);
        }
    }
}

/// Join the per-position debug content strings into `vref.debug`.
/// Ports `ToVarsBuilder.constructDebugLines` (ToVarsBuilder.java:1002-1014):
/// only when debug is enabled; empty list yields an empty string.
fn construct_debug_lines(debug_lines: &[String], vref: &mut Variant, enabled: bool) {
    if enabled {
        vref.debug = debug_lines.join(" & ");
    }
}

#[inline]
fn format_plus_genotype<D: std::fmt::Display>(
    buf: &mut crate::types::StackBuf<24>,
    val: D,
) -> std::borrow::Cow<'_, str> {
    use std::fmt::Write as _;
    *buf = crate::types::StackBuf::new();
    if write!(buf, "+{val}").is_ok() {
        std::borrow::Cow::Borrowed(buf.as_str())
    } else {
        std::borrow::Cow::Owned(format!("+{val}"))
    }
}

#[inline]
fn format_plus_genotype1<'b>(
    desc: &str,
    buf: &'b mut crate::types::StackBuf<24>,
    standard_dup_symbolic_alt: bool,
    canonical_sv_vcf_anchoring: bool,
    sv_respect_min_len: bool,
    sv_min_len: u32,
) -> std::borrow::Cow<'b, str> {
    let parse_neg_dup =
        standard_dup_symbolic_alt || canonical_sv_vcf_anchoring || sv_respect_min_len;
    if parse_neg_dup {
        if let Some(dup_val) = get_dup_num_compat(desc, true) {
            let dup_span = ((2 * SV_FLANK_BP) + i64::from(dup_val)).max(1);
            return format_plus_genotype(buf, dup_span);
        }
    } else if let Some(dup_val) = get_dup_num(desc) {
        return format_plus_genotype(buf, (dup_val as u32).saturating_add(SV_FLANK_BP as u32));
    }
    if canonical_sv_vcf_anchoring {
        if let Some(stripped_desc) = desc.strip_prefix('+') {
            let ins_len = long_ins_inserted_bases(stripped_desc, desc);
            let exceeds = if sv_respect_min_len {
                u32::try_from(ins_len).unwrap_or(u32::MAX) >= sv_min_len
            } else {
                stripped_desc.len().saturating_add(1) > sv_min_len as usize
            };
            if exceeds {
                return format_plus_genotype(buf, ins_len);
            }
        }
    }
    format_plus_genotype(buf, (desc.len() as i32).saturating_sub(1))
}

/// Strand bias calculation matching VarDict
pub fn strand_bias(fwd: i32, rev: i32, min_reads_bias: i32) -> i32 {
    crate::variations::util::strand_bias(fwd, rev, 0.05, min_reads_bias)
}

#[inline]
fn reduce_iupac_byte(b: u8) -> u8 {
    match b {
        b'M' | b'R' | b'W' | b'V' | b'H' | b'D' => b'A',
        b'S' | b'Y' | b'B' => b'C',
        b'K' => b'G',
        b'm' | b'r' | b'w' | b'v' | b'h' | b'd' => b'a',
        b's' | b'y' | b'b' => b'c',
        b'k' => b'g',
        other => other,
    }
}

#[inline]
fn reduce_iupac_byte_sw26(b: u8) -> u8 {
    match b {
        b'M' | b'R' | b'W' | b'V' | b'H' | b'D' | b'm' | b'r' | b'w' | b'v' | b'h' | b'd' => b'A',
        b'S' | b'Y' | b'B' | b's' | b'y' | b'b' => b'C',
        b'K' | b'k' => b'G',
        other => other,
    }
}

#[inline]
fn reduce_ref_byte_sw26(b: u8) -> u8 {
    reduce_iupac_byte_sw26(b).to_ascii_uppercase()
}

fn validate_ref_allele_in_place(ref_allele: &mut String) {
    if !ref_allele.bytes().any(|b| reduce_iupac_byte(b) != b) {
        return;
    }
    let mut bytes = std::mem::take(ref_allele).into_bytes();
    for b in &mut bytes {
        *b = reduce_iupac_byte(*b);
    }
    if let Ok(s) = String::from_utf8(bytes) {
        *ref_allele = s;
    }
}

fn validate_ref_allele_compat_in_place(ref_allele: &mut String, sanitize_iupac: bool) {
    if !sanitize_iupac {
        validate_ref_allele_in_place(ref_allele);
        return;
    }
    if !ref_allele.bytes().any(|b| reduce_ref_byte_sw26(b) != b) {
        return;
    }
    let mut bytes = std::mem::take(ref_allele).into_bytes();
    for b in &mut bytes {
        *b = reduce_ref_byte_sw26(*b);
    }
    if let Ok(s) = String::from_utf8(bytes) {
        *ref_allele = s;
    }
}

/// Sanitize `ALT` allele IUPAC ambiguity bases (`AstraZeneca-NGS/VarDictJava#267`,
/// `AstraZeneca-NGS/VarDictJava#287`, `AstraZeneca-NGS/VarDict#137`, `bcbio/bcbio-nextgen#3321`),
/// skipping closed `<...>` symbolic tokens (`<DUP>`, `<DEL>`, `<INV>`, `<dup-5>`, etc.) while
/// still reducing malformed unclosed `'<'` prefixes.
fn validate_alt_allele_in_place(alt_allele: &mut String) {
    if matches!(alt_allele.as_str(), "<DUP>" | "<DEL>" | "<INV>") {
        return;
    }
    if !alt_allele.bytes().any(|b| reduce_iupac_byte_sw26(b) != b) {
        return;
    }
    let mut bytes = std::mem::take(alt_allele).into_bytes();
    let mut i = 0usize;
    while let Some(&b) = bytes.get(i) {
        if b == b'<' {
            if let Some(rel_close) = bytes
                .get(i + 1..)
                .and_then(|tail| tail.iter().position(|&c| c == b'>'))
            {
                i = i.saturating_add(rel_close).saturating_add(2);
                continue;
            }
        }
        if let Some(slot) = bytes.get_mut(i) {
            *slot = reduce_iupac_byte_sw26(b);
        }
        i += 1;
    }
    if let Ok(s) = String::from_utf8(bytes) {
        *alt_allele = s;
    }
}

/// Validate reference allele sequence according to VCF 4.3 IUPAC ambiguity rules.
#[deprecated(
    since = "1.0.10",
    note = "internal allele normalization runs in place via ToVarsBuilder"
)]
pub fn validate_ref_allele(ref_allele: &str) -> String {
    let mut out = ref_allele.to_string();
    validate_ref_allele_in_place(&mut out);
    out
}

fn find_msi_slice<'a>(tseq1: &'a [u8], tseq2: &[u8], left: &[u8]) -> (f64, i32, &'a [u8]) {
    let mut msicnt = 0.0;
    let mut best_nmsi: usize = 0;

    let tseq1_len = tseq1.len();
    let limit = std::cmp::min(tseq1_len, 6);

    let mut stack_combined = [0u8; 256];
    let mut heap_combined = Vec::new();
    let combined_slice: &[u8] = if left.is_empty() || limit == 0 {
        tseq1
    } else {
        let total_len = left.len().saturating_add(tseq1_len);
        if let Some(dst) = stack_combined.get_mut(..total_len) {
            if let Some(left_dst) = dst.get_mut(..left.len()) {
                left_dst.copy_from_slice(left);
            }
            if let Some(right_dst) = dst.get_mut(left.len()..total_len) {
                right_dst.copy_from_slice(tseq1);
            }
            dst
        } else {
            heap_combined.reserve_exact(total_len);
            heap_combined.extend_from_slice(left);
            heap_combined.extend_from_slice(tseq1);
            heap_combined.as_slice()
        }
    };

    for nmsi in 1..=limit {
        let Some(msint) = tseq1.get(tseq1_len.saturating_sub(nmsi)..) else {
            continue;
        };

        let count = count_repetitions_at_end(combined_slice, msint);
        let mut curmsi = count as f64;
        curmsi += count_repetitions_at_start(tseq2, msint) as f64;

        if curmsi > msicnt {
            best_nmsi = nmsi;
            msicnt = curmsi;
        }
    }

    let mut shift3: usize = 0;
    if tseq1_len > 0 {
        while let Some(&t2_char) = tseq2.get(shift3) {
            let Some(&tseq_char) = (if shift3 < tseq1_len {
                tseq1.get(shift3)
            } else {
                tseq2.get(shift3 - tseq1_len)
            }) else {
                break;
            };
            if tseq_char == t2_char {
                shift3 += 1;
            } else {
                break;
            }
        }
    }

    let msint_bytes = if best_nmsi > 0 {
        tseq1
            .get(tseq1_len.saturating_sub(best_nmsi)..)
            .unwrap_or(&[])
    } else {
        &[]
    };

    (msicnt, shift3 as i32, msint_bytes)
}

/// Find microsatellite instability (MSI) motif and coordinate shift bounds.
#[deprecated(
    since = "1.0.10",
    note = "internal MSI computation runs allocation-free via ToVarsBuilder"
)]
pub fn find_msi(tseq1: &[u8], tseq2: &[u8], left: &[u8]) -> MsiResult {
    let (msi, shift3, msint_bytes) = find_msi_slice(tseq1, tseq2, left);
    MsiResult {
        msi,
        shift3,
        msint: String::from_utf8_lossy(msint_bytes).into_owned(),
    }
}

#[inline]
fn find_msi_metrics(tseq1: &[u8], tseq2: &[u8], left: &[u8]) -> (f64, i32, i32) {
    let (msi, shift3, msint_bytes) = find_msi_slice(tseq1, tseq2, left);
    let msint_len = if msint_bytes.is_ascii() {
        msint_bytes.len() as i32
    } else {
        String::from_utf8_lossy(msint_bytes).len() as i32
    };
    (msi, shift3, msint_len)
}

fn count_repetitions_at_end(seq: &[u8], pattern: &[u8]) -> usize {
    if pattern.is_empty() {
        return 0;
    }
    let mut count = 0;
    let mut current = seq;
    while current.ends_with(pattern) {
        count += 1;
        current = current.get(..current.len() - pattern.len()).unwrap_or(&[]);
    }
    count
}

fn count_repetitions_at_start(seq: &[u8], pattern: &[u8]) -> usize {
    if pattern.is_empty() {
        return 0;
    }
    let mut count = 0;
    let mut current = seq;
    while current.starts_with(pattern) {
        count += 1;
        current = current.get(pattern.len()..).unwrap_or(&[]);
    }
    count
}

fn get_dup_num_compat(s: &str, allow_negative: bool) -> Option<i32> {
    if let Some(idx) = s.find("<dup") {
        let rest = s.get(idx + 4..)?;
        if allow_negative {
            if let Some(after_minus) = rest.strip_prefix('-') {
                let digits_len = after_minus
                    .bytes()
                    .take_while(|b| b.is_ascii_digit())
                    .count();
                if digits_len > 0 {
                    let full_len = digits_len.saturating_add(1);
                    return rest.get(..full_len)?.parse().ok();
                }
                return None;
            }
        }
        let len = rest.bytes().take_while(|b| b.is_ascii_digit()).count();
        if len > 0 {
            return rest.get(..len)?.parse().ok();
        }
    }
    None
}

fn get_dup_num(s: &str) -> Option<i32> {
    get_dup_num_compat(s, false)
}

fn get_begin_minus_number(s: &str) -> Option<i32> {
    if let Some(rest) = s.strip_prefix('-') {
        let len = rest.bytes().take_while(|b| b.is_ascii_digit()).count();
        if len > 0 {
            return rest.get(..len)?.parse().ok();
        }
    }
    None
}

fn get_inv_num_compat(s: &str, allow_negative: bool) -> Option<i32> {
    if let Some(idx) = s.find("<inv") {
        let rest = s.get(idx + 4..)?;
        if allow_negative {
            if let Some(after_minus) = rest.strip_prefix('-') {
                let digits_len = after_minus
                    .bytes()
                    .take_while(|b| b.is_ascii_digit())
                    .count();
                if digits_len > 0 {
                    let full_len = digits_len.saturating_add(1);
                    return rest.get(..full_len)?.parse().ok();
                }
                return None;
            }
        }
        let len = rest.bytes().take_while(|b| b.is_ascii_digit()).count();
        if len > 0 {
            return rest.get(..len)?.parse().ok();
        }
    }
    None
}

fn matches_begin_minus_number_caret(s: &str) -> bool {
    if let Some(rest) = s.strip_prefix('-') {
        let digits_count = rest.bytes().take_while(|b| b.is_ascii_digit()).count();
        digits_count > 0 && rest.get(digits_count..).is_some_and(|r| r.starts_with('^'))
    } else {
        false
    }
}

fn matches_some_sv_numbers_compat(s: &str, allow_neg: bool) -> bool {
    if s.is_ascii() {
        let bytes = s.as_bytes();
        if bytes.len() < 6 {
            return false;
        }
        for i in 0..=bytes.len() - 6 {
            if bytes.get(i) == Some(&b'<') {
                if let Some(sv_type) = bytes.get(i + 1..i + 4) {
                    if sv_type.iter().all(|&c| c != b'<' && c != b'>') {
                        let mut j = i + 4;
                        if allow_neg && bytes.get(j) == Some(&b'-') {
                            j += 1;
                        }
                        let mut has_digits = false;
                        while let Some(&b) = bytes.get(j) {
                            if b.is_ascii_digit() {
                                has_digits = true;
                                j += 1;
                            } else {
                                break;
                            }
                        }
                        if has_digits && bytes.get(j) == Some(&b'>') {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
    let chars: Vec<char> = s.chars().collect();
    if chars.len() < 6 {
        return false;
    }
    for i in 0..=chars.len() - 6 {
        if chars.get(i) == Some(&'<') {
            if let Some(sv_type) = chars.get(i + 1..i + 4) {
                if sv_type.iter().all(|&c| c != '<' && c != '>') {
                    let mut j = i + 4;
                    if allow_neg && chars.get(j) == Some(&'-') {
                        j += 1;
                    }
                    let mut has_digits = false;
                    while let Some(&c) = chars.get(j) {
                        if c.is_ascii_digit() {
                            has_digits = true;
                            j += 1;
                        } else {
                            break;
                        }
                    }
                    if has_digits && chars.get(j) == Some(&'>') {
                        return true;
                    }
                }
            }
        }
    }
    false
}

fn get_amp_atgc(s: &str) -> Option<&str> {
    let (_, rest) = s.split_once('&')?;
    let len = rest
        .bytes()
        .take_while(|&b| matches!(b, b'A' | b'T' | b'G' | b'C'))
        .count();
    if len == 0 {
        None
    } else {
        rest.get(..len)
    }
}

fn get_hash_group_caret_group(s: &str) -> Option<(&str, &str)> {
    let (_, rest) = s.split_once('#')?;
    let (g1, g2) = rest.split_once('^')?;
    if !g1.is_empty() && !g2.is_empty() {
        Some((g1, g2))
    } else {
        None
    }
}

fn get_begin_digits(s: &str) -> Option<i32> {
    let len = s.bytes().take_while(|b| b.is_ascii_digit()).count();
    if len > 0 {
        s.get(..len)?.parse().ok()
    } else {
        None
    }
}

fn get_caret_atgnc(s: &str) -> Option<&str> {
    let (_, rest) = s.split_once('^')?;
    let len = rest
        .bytes()
        .take_while(|&b| matches!(b, b'A' | b'T' | b'G' | b'C' | b'N'))
        .count();
    if len == 0 {
        None
    } else {
        rest.get(..len)
    }
}

fn matches_any_sv(s: &str) -> bool {
    if s.is_ascii() {
        for w in s.as_bytes().windows(5) {
            if let [b'<', c1, c2, c3, b'>'] = *w {
                if c1 != b'<' && c1 != b'>' && c2 != b'<' && c2 != b'>' && c3 != b'<' && c3 != b'>'
                {
                    return true;
                }
            }
        }
        return false;
    }
    let chars: Vec<char> = s.chars().collect();
    for w in chars.windows(5) {
        if let ['<', c1, c2, c3, '>'] = *w {
            if c1 != '<' && c1 != '>' && c2 != '<' && c2 != '>' && c3 != '<' && c3 != '>' {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::regions::Region;
    use crate::types::Pos1;
    use crate::variations::{Sv, Variation, VariationData, VariationMap};

    #[test]
    fn test_deldupvar_boundary_filtering() {
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 300,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 300]);
        ref_genome.region_sequences.push(vec![b'A'; 300]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 300)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(200),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut var_data = VariationData::new();
        let pos_outside = 250;
        var_data.ref_coverage.insert(pos_outside, 20);

        let mut vm = VariationMap::new();
        vm.sv = Some(Box::new(Sv {
            r#type: "INV".to_string(),
            pairs: 2,
            splits: 3,
            clusters: 1,
        }));
        let mut v = Variation::new();
        v.cnt = 10;
        v.fwd = 5;
        v.rev = 5;
        vm.map.insert("G".to_string(), v);
        var_data.non_insertion_variants.insert(pos_outside, vm);

        // When deldupvar = false (default), SV outside region boundaries is kept
        let mut cfg_no_del = Config::default();
        cfg_no_del.deldupvar = false;
        cfg_no_del.pileup = true;
        let builder_no_del = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_no_del);
        let aligned_no_del = builder_no_del.build_aligned_variants(&var_data);
        assert!(
            aligned_no_del.contains_key(&pos_outside),
            "Expected SV outside region to be retained when deldupvar=false"
        );

        // When deldupvar = true, SV outside region boundaries is discarded
        let mut cfg_del = Config::default();
        cfg_del.deldupvar = true;
        cfg_del.pileup = true;
        let builder_del = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_del);
        let aligned_del = builder_del.build_aligned_variants(&var_data);
        assert!(
            !aligned_del.contains_key(&pos_outside),
            "Expected SV outside region to be filtered out when deldupvar=true"
        );
    }

    /// Ports SomaticMode.java:111-113 / ToVarsBuilder.java: when a deletion
    /// extends past the loaded reference window (e.g. a realignlgdel call
    /// whose end exceeds region.end + reference_extension), ToVarsBuilder
    /// observes unloaded bases after `end_position` and leaves `right_seq`
    /// empty (`""`, serialized as `0` / `RSEQ=0`).
    #[test]
    fn unloaded_deletion_end_leaves_right_seq_empty() {
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 1000,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 1000]);
        ref_genome.region_sequences.push(vec![b'A'; 1000]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 2000)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(1000),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut var_data = VariationData::new();
        var_data.ref_coverage.insert(990, 20);
        let mut vm = VariationMap::new();
        let mut v = Variation::new();
        v.cnt = 10;
        v.fwd = 5;
        v.rev = 5;
        v.qmean = 300.0;
        v.mqmean = 600.0;
        v.pmean = 250.0;
        v.hicnt = 10;
        vm.map.insert("-100".to_string(), v);
        var_data.non_insertion_variants.insert(990, vm);

        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);
        let aligned = builder.build_aligned_variants(&var_data);
        let vars = aligned.get(&990).expect("variant at 990 should be built");
        assert_eq!(vars.variants.len(), 1);
        assert_eq!(vars.variants[0].end_pos, 1089);
        assert_eq!(vars.variants[0].right_seq, "");
    }

    fn is_the_same_variation_on_ref_legacy(
        builder: &ToVarsBuilder<'_>,
        position: i32,
        vars_at_cur_position: &crate::variations::VariationMap,
        insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    ) -> bool {
        let mut vk: std::collections::HashSet<String> =
            vars_at_cur_position.map.keys().map(str::to_owned).collect();
        if insertion_variants.contains_key(&position) {
            vk.insert("I".to_string());
        }
        if vk.len() == 1 {
            if let Some(ref_base) = builder
                .reference
                .get_base(&builder.region.chr, position as u32)
            {
                let ref_base_str = (ref_base as char).to_string();
                if vk.contains(&ref_base_str) {
                    let has_bam2 = builder.config.bam.contains('|');
                    let do_pileup = builder.config.pileup;
                    let is_amplicon = builder.config.amplicon.is_some();
                    if !do_pileup && !has_bam2 && !is_amplicon {
                        return true;
                    }
                }
            }
        }
        false
    }

    #[test]
    fn is_the_same_variation_on_ref_matches_legacy_hashset() {
        let key_sets: &[&[&str]] = &[
            &[],
            &["A"],
            &["I"],
            &["C"],
            &["A", "G"],
            &["A", "I"],
            &["I", "G"],
        ];
        let ref_bases: &[Option<u8>] = &[Some(b'A'), Some(b'I'), Some(b'C'), Some(0xE9), None];
        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        for &ref_base in ref_bases {
            let mut ref_genome = Reference::new();
            if let Some(b) = ref_base {
                ref_genome
                    .loaded_regions
                    .push(crate::reference::LoadedRegion {
                        chr: "chr1".to_string(),
                        sequence_start: 1,
                        sequence_end: 10,
                    });
                ref_genome
                    .reference_sequences
                    .insert("chr1".to_string(), vec![b; 10]);
                ref_genome.region_sequences.push(vec![b; 10]);
            }

            for &keys in key_sets {
                let mut vm = VariationMap::new();
                for &k in keys {
                    vm.map.insert(k.to_string(), Variation::new());
                }

                for has_ins in [false, true] {
                    let mut ins_map = crate::javamap::JavaHashMap::new();
                    if has_ins {
                        ins_map.insert(5, VariationMap::new());
                    }

                    for (pileup, bam, amplicon) in [
                        (false, "sample.bam", None),
                        (true, "sample.bam", None),
                        (false, "tumor.bam|normal.bam", None),
                        (false, "sample.bam", Some("10:0.95".to_string())),
                    ] {
                        let mut cfg = Config::default();
                        cfg.pileup = pileup;
                        cfg.bam = bam.to_string();
                        cfg.amplicon = amplicon;

                        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);
                        let actual = builder.is_the_same_variation_on_ref(5, &vm, &ins_map);
                        let expected =
                            is_the_same_variation_on_ref_legacy(&builder, 5, &vm, &ins_map);
                        assert_eq!(
                            actual, expected,
                            "mismatch for keys={keys:?}, has_ins={has_ins}, ref_base={ref_base:?}, pileup={pileup}, bam={bam}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn build_aligned_variants_cow_overlay_decrements_next_pos_and_preserves_input() {
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);
        ref_genome.region_sequences.push(vec![b'A'; 100]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut var_data = VariationData::new();
        // Position 10: total_cov = 10, insertion "+T" has cnt = 15 (> total_cov, extracnt == 0),
        // fwd = 9, rev = 6.
        // Position 11: ref_coverage = 40 (so ttcov=15 < 40 - 15 = 25), triggering the P+1
        // ref base ('A') fwd/rev decrement by (9, 6).
        var_data.ref_coverage.insert(10, 10);
        var_data.ref_coverage.insert(11, 40);

        let mut vm10 = VariationMap::new();
        let mut ref10 = Variation::new();
        ref10.cnt = 10;
        ref10.fwd = 5;
        ref10.rev = 5;
        ref10.qmean = 300.0;
        ref10.mqmean = 600.0;
        ref10.pmean = 200.0;
        ref10.hicnt = 10;
        vm10.map.insert("A".to_string(), ref10);
        var_data.non_insertion_variants.insert(10, vm10);

        let mut ins10_map = VariationMap::new();
        let mut ins10 = Variation::new();
        ins10.cnt = 15;
        ins10.fwd = 9;
        ins10.rev = 6;
        ins10.extracnt = 0;
        ins10.qmean = 450.0;
        ins10.mqmean = 900.0;
        ins10.pmean = 300.0;
        ins10.hicnt = 15;
        ins10_map.map.insert("+T".to_string(), ins10);
        var_data.insertion_variants.insert(10, ins10_map);

        let mut vm11 = VariationMap::new();
        let mut ref11 = Variation::new();
        ref11.cnt = 30;
        ref11.fwd = 20;
        ref11.rev = 10;
        ref11.qmean = 900.0;
        ref11.mqmean = 1800.0;
        ref11.pmean = 600.0;
        ref11.hicnt = 30;
        vm11.map.insert("A".to_string(), ref11);
        let mut alt11 = Variation::new();
        alt11.cnt = 10;
        alt11.fwd = 5;
        alt11.rev = 5;
        alt11.qmean = 300.0;
        alt11.mqmean = 600.0;
        alt11.pmean = 200.0;
        alt11.hicnt = 10;
        vm11.map.insert("G".to_string(), alt11);
        var_data.non_insertion_variants.insert(11, vm11);

        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);
        let aligned = builder.build_aligned_variants(&var_data);

        // Input var_data at position 11 must remain untouched (fwd=20, rev=10)
        let orig_ref11 = var_data
            .non_insertion_variants
            .get(&11)
            .and_then(|m| m.map.get("A"))
            .expect("orig ref11");
        assert_eq!((orig_ref11.fwd, orig_ref11.rev), (20, 10));

        // Position 10's insertion variant should see adjusted ref cov from P+1: (20-9, 10-6) = (11, 4)
        let vars10 = aligned.get(&10).expect("vars at 10");
        assert_eq!(vars10.variants.len(), 1);
        assert_eq!(
            (
                vars10.variants[0].ref_fwd_cov,
                vars10.variants[0].ref_rev_cov
            ),
            (11, 4)
        );
        assert!(vars10.variant_by_desc("+T").is_some());

        // Position 11's finalized variants must also see the decremented ref_variant (11, 4)
        let vars11 = aligned.get(&11).expect("vars at 11");
        let ref_var11 = vars11.ref_variant.as_ref().expect("ref_variant at 11");
        assert_eq!((ref_var11.fwd, ref_var11.rev), (11, 4));
        assert_eq!(
            (
                vars11.variants[0].ref_fwd_cov,
                vars11.variants[0].ref_rev_cov
            ),
            (11, 4)
        );
        assert!(vars11.variant_by_desc("G").is_some());

        // Issue #60: ToVarsBuilder.java:337-340 performs signed subtraction
        // (`nextVar.fwd -= fwd; nextVar.rev -= rev;`) without clamping at 0.
        let mut var_data_signed = var_data.clone();
        var_data_signed
            .insertion_variants
            .get_mut(&10)
            .unwrap()
            .map
            .get_mut("+T")
            .unwrap()
            .fwd = 12;
        var_data_signed
            .insertion_variants
            .get_mut(&10)
            .unwrap()
            .map
            .get_mut("+T")
            .unwrap()
            .rev = 3;
        var_data_signed
            .non_insertion_variants
            .get_mut(&11)
            .unwrap()
            .map
            .get_mut("A")
            .unwrap()
            .fwd = 10;
        var_data_signed
            .non_insertion_variants
            .get_mut(&11)
            .unwrap()
            .map
            .get_mut("A")
            .unwrap()
            .rev = 20;
        let aligned_signed = builder.build_aligned_variants(&var_data_signed);
        let ref_var11_signed = aligned_signed
            .get(&11)
            .and_then(|v| v.ref_variant.as_ref())
            .expect("ref_variant at 11");
        assert_eq!((ref_var11_signed.fwd, ref_var11_signed.rev), (-2, 17));
    }

    #[test]
    fn build_aligned_variants_degenerate_plus_allele_does_not_panic() {
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);
        ref_genome.region_sequences.push(vec![b'A'; 100]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut var_data = VariationData::new();
        var_data.ref_coverage.insert(10, 20);

        let mut vm10 = VariationMap::new();
        let mut ref10 = Variation::new();
        ref10.cnt = 10;
        ref10.fwd = 5;
        ref10.rev = 5;
        ref10.qmean = 300.0;
        ref10.mqmean = 600.0;
        ref10.pmean = 200.0;
        ref10.hicnt = 10;
        vm10.map.insert("A".to_string(), ref10.clone());
        // Also test a degenerate "+" in non_insertion_variants and insertion_variants
        vm10.map.insert("+".to_string(), ref10.clone());
        var_data.non_insertion_variants.insert(10, vm10);

        let mut ins10_map = VariationMap::new();
        ins10_map.map.insert("+".to_string(), ref10);
        var_data.insertion_variants.insert(10, ins10_map);

        let mut cfg = Config::default();
        cfg.crispr_cutting_site = Some(12);
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);
        let aligned = builder.build_aligned_variants(&var_data);
        for (_, vars) in aligned.iter() {
            for v in &vars.variants {
                assert!(v.msi.is_finite(), "MSI must remain finite, got {}", v.msi);
            }
        }
    }

    #[test]
    fn h1_two_allele_low_af_drain_row_set_and_mnp_tiebreak_order() {
        // H1 pin: two-allele low-AF pileup drain. `var` arrives sorted
        // (MNP tie-break shape: tied mean_quality*pos_coverage broken by
        // description, CPP-9); the drain must split ref vs non-ref rows,
        // keep the sorted order in `entry.variants`, and report maxfreq.
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);
        ref_genome.region_sequences.push(vec![b'A'; 100]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);

        let mk = |description: &str, frequency: f64, mean_quality: f64, pos_coverage: i32| {
            let mut v = Variant::new();
            v.description = description.to_string();
            v.frequency = frequency;
            v.mean_quality = mean_quality;
            v.pos_coverage = pos_coverage;
            v
        };
        // Pileup order: MNP alt, ref, SNV alt; tied metric 60 on the alts,
        // ref metric 50 so it sorts last (metric value is synthetic here).
        let mut var = vec![
            mk("TA", 0.02, 30.0, 2),
            mk("A", 0.90, 10.0, 5),
            mk("G", 0.03, 20.0, 3),
        ];
        builder.sort_variants(&mut var);
        assert_eq!(
            var.iter()
                .map(|v| v.description.as_str())
                .collect::<Vec<_>>(),
            vec!["G", "TA", "A"],
            "MNP tie-break: tied metric must sort by description"
        );

        let mut aligned = crate::javamap::JavaHashMap::new();
        let maxfreq = builder.collect_vars_at_position(&mut aligned, 10, &mut var);

        assert!(var.is_empty(), "drain must consume every row");
        assert_eq!(maxfreq, 0.03);
        let entry = aligned.get(&10).expect("entry at 10");
        assert_eq!(
            entry.ref_variant.as_ref().expect("ref row").description,
            "A"
        );
        assert_eq!(
            entry
                .variants
                .iter()
                .map(|v| v.description.as_str())
                .collect::<Vec<_>>(),
            vec!["G", "TA"]
        );
    }

    #[test]
    fn h1_early_reserve_preserves_hashmap_resize_order_across_threshold() {
        // Regression for the H1 pre-reserve `compute_if_absent`: it must not
        // shift any key insertion relative to the drain-only path. 13 keys
        // cross the 16-bucket / 0.75-load-factor resize boundary (threshold
        // 12), where a moved insert would change JavaHashMap iteration order
        // and therefore `.var` row order (UB-3). Same key, same position in
        // the sequence -> identical table layout.
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);
        ref_genome.region_sequences.push(vec![b'A'; 100]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);

        let positions: Vec<i32> = (1..=12).chain(std::iter::once(17)).collect();
        let mk_row = || {
            let mut v = Variant::new();
            v.description = "G".to_string(); // non-ref (ref base is A)
            v.frequency = 0.1;
            v.mean_quality = 30.0;
            v.pos_coverage = 2;
            v
        };

        // New path: early reserve + drain (one non-ref row per position).
        let mut aligned = crate::javamap::JavaHashMap::new();
        for &p in &positions {
            let mut var = vec![mk_row()];
            builder.collect_vars_at_position(&mut aligned, p, &mut var);
        }

        // Legacy path: drain-only insertion (first row's get_or_put_with).
        let mut legacy = crate::javamap::JavaHashMap::new();
        for &p in &positions {
            legacy
                .get_or_put_with(p, crate::variations::Vars::new)
                .variants
                .push(mk_row());
        }

        assert_eq!(aligned.len(), 13);
        assert_eq!(legacy.len(), 13);
        assert_eq!(
            aligned.keys().copied().collect::<Vec<i32>>(),
            legacy.keys().copied().collect::<Vec<i32>>(),
            "pre-reserve must not change JavaHashMap iteration order across the resize"
        );
    }

    #[test]
    fn zero_alloc_pattern_helpers_match_empty_and_non_matching() {
        // get_amp_atgc
        assert_eq!(get_amp_atgc("-3&ATGC"), Some("ATGC"));
        assert_eq!(get_amp_atgc("+AT&GC&TT"), Some("GC"));
        assert_eq!(get_amp_atgc("&A#AC^2"), Some("A"));
        assert_eq!(get_amp_atgc(""), None);
        assert_eq!(get_amp_atgc("ATGC"), None);
        assert_eq!(get_amp_atgc("&"), None);
        assert_eq!(get_amp_atgc("&N"), None);
        assert_eq!(get_amp_atgc("&123"), None);

        // get_hash_group_caret_group (ports Patterns.HASH_GROUP_CARET_GROUP `#(.+)\^(.+)`)
        assert_eq!(get_hash_group_caret_group("-12#AC^3"), Some(("AC", "3")));
        assert_eq!(
            get_hash_group_caret_group("-5#ATGCN^12ACN"),
            Some(("ATGCN", "12ACN"))
        );
        assert_eq!(get_hash_group_caret_group("-2#AR^3"), Some(("AR", "3")));
        assert_eq!(get_hash_group_caret_group("#A^G&TT"), Some(("A", "G&TT")));
        assert_eq!(get_hash_group_caret_group(""), None);
        assert_eq!(get_hash_group_caret_group("#^3"), None);
        assert_eq!(get_hash_group_caret_group("#AC^"), None);
        assert_eq!(get_hash_group_caret_group("AC^3"), None);
        assert_eq!(get_hash_group_caret_group("#AC"), None);

        // get_caret_atgnc
        assert_eq!(get_caret_atgnc("-3^ATGCN"), Some("ATGCN"));
        assert_eq!(get_caret_atgnc("^AC12"), Some("AC"));
        assert_eq!(get_caret_atgnc(""), None);
        assert_eq!(get_caret_atgnc("ATGCN"), None);
        assert_eq!(get_caret_atgnc("^"), None);
        assert_eq!(get_caret_atgnc("^-3"), None);
        assert_eq!(get_caret_atgnc("^2"), None);
    }

    #[test]
    fn deletion_prefix_strip_handles_digits_and_non_digits() {
        let strip_del = |description_string: &str| -> String {
            if let Some(rest) = description_string.strip_prefix('-') {
                let digit_len = rest.bytes().take_while(|b| b.is_ascii_digit()).count();
                if digit_len > 0 {
                    rest.get(digit_len..).unwrap_or("").to_string()
                } else {
                    description_string.to_string()
                }
            } else {
                description_string.to_string()
            }
        };

        assert_eq!(strip_del("-3"), "");
        assert_eq!(strip_del("-12&AT"), "&AT");
        assert_eq!(strip_del("-3^2"), "^2");
        assert_eq!(strip_del("-12#AC^3"), "#AC^3");
        assert_eq!(strip_del("-A"), "-A");
        assert_eq!(strip_del("-"), "-");
    }

    #[test]
    fn primary_ref_slice_cache_and_join_ref_into_match_reference_across_all_layouts() {
        let chr_lengths = HashMap::from([("chr1".to_string(), 200)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(10),
            end: Pos1(120),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let cfg = Config::default();

        // (a) & (b): Primary loaded region [10..=100] with 17-bp SEED1_LEN truncated tail
        // (stored region_sequences[0] has length 74 covering [10..=83], while [84..=100] is
        // truncated in primary and covered by secondary loaded_regions[1] covering [80..=130]).
        let mut ref_with_tail = Reference::new();
        ref_with_tail
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 10,
                sequence_end: 100,
            });
        let primary_bases: Vec<u8> = (0..74)
            .map(|i| b"ACGT"[i % 4]) // 74 bp instead of 91 bp (17-bp SEED1_LEN truncated tail)
            .collect();
        ref_with_tail
            .reference_sequences
            .insert("chr1".to_string(), primary_bases.clone());
        ref_with_tail.region_sequences.push(primary_bases);

        ref_with_tail
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 80,
                sequence_end: 130,
            });
        let secondary_bases: Vec<u8> = (0..=50).map(|i| b"TGCA"[i % 4]).collect();
        ref_with_tail.region_sequences.push(secondary_bases);

        // (c): reference_sequences-only setup (region_sequences is empty)
        let mut ref_map_only = Reference::new();
        ref_map_only
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 10,
                sequence_end: 60,
            });
        let map_only_bases: Vec<u8> = (0..=50).map(|i| b"CGAT"[i % 4]).collect();
        ref_map_only
            .reference_sequences
            .insert("chr1".to_string(), map_only_bases);

        for reference in [&ref_with_tail, &ref_map_only] {
            let builder = ToVarsBuilder::new(reference, &region, &chr_lengths, &cfg);

            // Check single-base lookup across OOB before start, primary slice, 17-bp truncated tail,
            // secondary region, and OOB past end.
            for pos in 0..=150u32 {
                assert_eq!(
                    builder.get_ref_base(pos),
                    reference.get_base(&region.chr, pos),
                    "get_ref_base mismatch at pos={pos}"
                );
            }

            // (d): Check join_ref_into across negative/OOB, reversed, primary, tail-spanning, and secondary ranges.
            let ranges = [
                (-20, -5),
                (-10, 0),
                (-5, 15),
                (1, 9),
                (10, 50),
                (70, 95), // spans primary into 17-bp truncated tail covered by loaded_regions[1]
                (84, 100), // strictly inside the 17-bp truncated tail of loaded_regions[0]
                (95, 140), // spans secondary into OOB past end
                (50, 20), // end < start
                (135, 160),
            ];
            let mut actual = String::new();
            let mut expected = String::new();
            for (start, end) in ranges {
                builder.join_ref_into(&mut actual, start, end);
                crate::util::ref_bases_into(&mut expected, reference, &region.chr, start, end);
                assert_eq!(
                    actual, expected,
                    "join_ref_into mismatch for range {start}..={end}"
                );
            }
        }
    }

    #[test]
    fn sorted_key_slice_stack_and_heap_parity() {
        let pool: [&str; 12] = [
            "T", "G", "C", "A", "-3", "-1", "&AT", "#AC^2", "^GC", "GA", "CT", "TT",
        ];
        for n in 0..=pool.len() {
            let sub = pool.get(..n).unwrap_or(&[]);
            for extra in [None, Some("SV")] {
                let sorted = SortedKeySlice::new(sub.iter().copied(), extra);
                let total_len = n + usize::from(extra.is_some());
                if total_len <= 8 {
                    assert!(
                        matches!(sorted, SortedKeySlice::Stack(_, len) if len == total_len),
                        "expected Stack variant for total_len={total_len}"
                    );
                } else {
                    assert!(
                        matches!(sorted, SortedKeySlice::Heap(ref v) if v.len() == total_len),
                        "expected Heap variant for total_len={total_len}"
                    );
                }

                let mut expected: Vec<String> = sub.iter().map(|s| (*s).to_string()).collect();
                if let Some(ex) = extra {
                    expected.push(ex.to_string());
                }
                expected.sort();
                let expected_refs: Vec<&str> = expected.iter().map(String::as_str).collect();
                assert_eq!(sorted.as_slice(), expected_refs.as_slice());
            }
        }

        // End-to-end build_aligned_variants check across both <=8 (Stack) and >8 (Heap) + SV
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);
        ref_genome.region_sequences.push(vec![b'A'; 100]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };
        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);

        let mk_var = |cnt: i32, q: f64| {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.fwd = cnt;
            v.rev = 0;
            v.qmean = q * cnt as f64;
            v.mqmean = 60.0 * cnt as f64;
            v.pmean = 20.0 * cnt as f64;
            v.hicnt = cnt;
            v
        };

        let mut var_data = VariationData::new();
        var_data.ref_coverage.insert(20, 100);
        let mut vm20 = VariationMap::new();
        vm20.sv = Some(Box::new(Sv {
            r#type: "INV".to_string(),
            splits: 2,
            pairs: 4,
            clusters: 1,
        }));
        for (idx, key) in ["A", "C", "G", "T", "-1", "-2", "&GT", "CA", "GA", "TA"]
            .iter()
            .enumerate()
        {
            vm20.map
                .insert((*key).to_string(), mk_var(10, 30.0 + idx as f64));
        }
        var_data.non_insertion_variants.insert(20, vm20);

        let mut ins20 = VariationMap::new();
        for (idx, key) in [
            "+A", "+C", "+G", "+T", "+AA", "+CC", "+GG", "+TT", "+AC", "+GT",
        ]
        .iter()
        .enumerate()
        {
            ins20
                .map
                .insert((*key).to_string(), mk_var(5, 25.0 + idx as f64));
        }
        var_data.insertion_variants.insert(20, ins20);

        let aligned = builder.build_aligned_variants(&var_data);
        let vars20 = aligned.get(&20).expect("vars at 20");
        assert_eq!(vars20.sv, "2-4-1");
        assert_eq!(vars20.variants.len(), 19);
        assert_eq!(
            vars20
                .variants
                .iter()
                .find(|v| v.description == "&GT")
                .and_then(|v| v.genotype.as_deref()),
            Some("AAA/GT")
        );
    }

    #[test]
    #[expect(
        deprecated,
        reason = "verifies 1.0.x SemVer shims find_msi and validate_ref_allele"
    )]
    fn find_msi_metrics_and_validate_ref_allele_in_place_parity() {
        // 1. Clean ACGTN fast path preserves buffer pointer and capacity
        let mut clean = String::with_capacity(64);
        clean.push_str("ACGTNACGTN");
        let ptr_before = clean.as_ptr();
        let cap_before = clean.capacity();
        validate_ref_allele_in_place(&mut clean);
        assert_eq!(clean, "ACGTNACGTN");
        assert_eq!(clean.as_ptr(), ptr_before);
        assert_eq!(clean.capacity(), cap_before);

        let mut clean_lower = String::with_capacity(32);
        clean_lower.push_str("acgtn");
        let lower_ptr = clean_lower.as_ptr();
        let lower_cap = clean_lower.capacity();
        validate_ref_allele_in_place(&mut clean_lower);
        assert_eq!(clean_lower, "acgtn");
        assert_eq!(clean_lower.as_ptr(), lower_ptr);
        assert_eq!(clean_lower.capacity(), lower_cap);

        // 2. IUPAC ambiguity remapping preserves capacity and matches validate_ref_allele
        let iupac_input = "ACGTRYSWKMBDHVNacgtryswkmbdhvn";
        let expected_iupac = "ACGTACCAGACAAANacgtaccagacaaan";
        let mut iupac = String::with_capacity(128);
        iupac.push_str(iupac_input);
        let iupac_cap = iupac.capacity();
        validate_ref_allele_in_place(&mut iupac);
        assert_eq!(iupac, expected_iupac);
        assert_eq!(iupac, validate_ref_allele(iupac_input));
        assert_eq!(iupac.capacity(), iupac_cap);

        // Non-ASCII UTF-8 + IUPAC parity
        let mut utf8_iupac = String::from("éMßkACGT");
        validate_ref_allele_in_place(&mut utf8_iupac);
        assert_eq!(utf8_iupac, "éAßgACGT");
        assert_eq!(utf8_iupac, validate_ref_allele("éMßkACGT"));

        // Exhaustive ASCII byte parity check (0..=127)
        for b in 0u8..=127u8 {
            let ch = b as char;
            let legacy_mapped = match ch {
                'M' | 'R' | 'W' | 'V' | 'H' | 'D' => 'A',
                'S' | 'Y' | 'B' => 'C',
                'K' => 'G',
                'm' | 'r' | 'w' | 'v' | 'h' | 'd' => 'a',
                's' | 'y' | 'b' => 'c',
                'k' => 'g',
                other => other,
            };
            let mut s = ch.to_string();
            validate_ref_allele_in_place(&mut s);
            assert_eq!(s, legacy_mapped.to_string(), "byte {b} ({ch:?}) mismatch");
        }

        // 3. find_msi_metrics vs find_msi parity across empty, repeats, >256B heap, and lossy UTF-8
        let long_left = vec![b'A'; 300];
        let cases: &[(&[u8], &[u8], &[u8])] = &[
            (b"", b"ACGT", b"ACGT"),
            (b"A", b"AAAAACGT", b"TTTAAA"),
            (b"AT", b"ATATATGC", b"CCATAT"),
            (b"CAG", b"CAGCAGCAGTT", b"GGCAGCAG"),
            (b"ACGTAC", b"ACGTACACGTAC", b"ACGTAC"),
            (b"AA", b"AAAA", &long_left),
            (&[0xFF, b'A'], &[0xFF, b'A', 0xFF, b'A'], &[0xFF, b'A']),
        ];
        for &(tseq1, tseq2, left) in cases {
            let full = find_msi(tseq1, tseq2, left);
            let (msi, shift3, msint_len) = find_msi_metrics(tseq1, tseq2, left);
            assert_eq!(msi, full.msi);
            assert_eq!(shift3, full.shift3);
            assert_eq!(msint_len, full.msint.len() as i32);
        }
    }

    #[test]
    fn tovars_omits_desc_map_clones() {
        let mut ref_genome = Reference::new();
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let chr_len_usize = seq.len();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: chr_len_usize as u32,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);

        let chr_lengths = HashMap::from([("chr1".to_string(), chr_len_usize)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(chr_len_usize as u32),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut var_data = VariationData::new();
        // Position 10 (ref 'C'): ref + SNV 'T' + insertion '+GG'
        var_data.ref_coverage.insert(10, 35);
        let mut vm10 = VariationMap::new();
        let mut ref10 = Variation::new();
        ref10.cnt = 20;
        ref10.fwd = 12;
        ref10.rev = 8;
        ref10.qmean = 600.0;
        ref10.mqmean = 1200.0;
        ref10.pmean = 400.0;
        ref10.hicnt = 20;
        vm10.map.insert("C".to_string(), ref10);

        let mut snv10 = Variation::new();
        snv10.cnt = 15;
        snv10.fwd = 15;
        snv10.rev = 0; // total > 12, rev == 0 => strand bias = 1
        snv10.qmean = 480.0;
        snv10.mqmean = 900.0;
        snv10.pmean = 300.0;
        snv10.hicnt = 15;
        vm10.map.insert("T".to_string(), snv10);
        var_data.non_insertion_variants.insert(10, vm10);

        let mut ins10_map = VariationMap::new();
        let mut ins10 = Variation::new();
        ins10.cnt = 6;
        ins10.fwd = 3;
        ins10.rev = 3;
        ins10.qmean = 180.0;
        ins10.mqmean = 360.0;
        ins10.pmean = 120.0;
        ins10.hicnt = 6;
        ins10_map.map.insert("+GG".to_string(), ins10);
        var_data.insertion_variants.insert(10, ins10_map);

        // Position 20 (no ref variant, deletion '-2'): exercises 0;{bias} and chr_len clamping
        var_data.ref_coverage.insert(20, 8);
        let mut vm20 = VariationMap::new();
        let mut del20 = Variation::new();
        del20.cnt = 8;
        del20.fwd = 4;
        del20.rev = 4;
        del20.qmean = 240.0;
        del20.mqmean = 480.0;
        del20.pmean = 160.0;
        del20.hicnt = 8;
        vm20.map.insert("-2".to_string(), del20);
        var_data.non_insertion_variants.insert(20, vm20);

        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);
        let aligned = builder.build_aligned_variants(&var_data);

        let vars10 = aligned.get(&10).expect("vars at 10");
        assert!(
            vars10.var_desc_to_variant.is_empty(),
            "ToVarsBuilder must keep var_desc_to_variant empty"
        );
        assert_eq!(vars10.variants.len(), 2);
        for v in &vars10.variants {
            let resolved = vars10
                .variant_by_desc(&v.description)
                .expect("variant_by_desc must resolve emitted variant");
            assert_eq!(resolved.description, v.description);
            assert_eq!(resolved.start_pos, v.start_pos);
            assert_eq!(resolved.end_pos, v.end_pos);
            assert_eq!(resolved.ref_allele, v.ref_allele);
            assert_eq!(resolved.var_allele, v.var_allele);
            assert_eq!(resolved.bias, v.bias);
            assert_eq!(resolved.pos_coverage, v.pos_coverage);
            assert_eq!(resolved.total_pos_cov, v.total_pos_cov);
        }
        assert_eq!(vars10.variant_by_desc("T").unwrap().bias, "2;1");
        assert_eq!(vars10.variant_by_desc("+GG").unwrap().bias, "2;2");
        assert!(vars10.variant_by_desc("A").is_none());

        let vars20 = aligned.get(&20).expect("vars at 20");
        assert!(vars20.var_desc_to_variant.is_empty());
        assert_eq!(vars20.variant_by_desc("-2").unwrap().bias, "0;2");

        // --- Somatic RefSiteMode / CompactRefSite + hydrate_compact_ref_site parity ---
        let assert_variant_eq = |a: &Variant, b: &Variant, ctx: &str| {
            assert_eq!(a.description, b.description, "{ctx}: description");
            assert_eq!(a.pos_coverage, b.pos_coverage, "{ctx}: pos_coverage");
            assert_eq!(a.fwd, b.fwd, "{ctx}: fwd");
            assert_eq!(a.rev, b.rev, "{ctx}: rev");
            assert_eq!(a.bias, b.bias, "{ctx}: bias");
            assert_eq!(
                a.frequency.to_bits(),
                b.frequency.to_bits(),
                "{ctx}: frequency"
            );
            assert_eq!(
                a.mean_position.to_bits(),
                b.mean_position.to_bits(),
                "{ctx}: mean_position"
            );
            assert_eq!(a.pstd, b.pstd, "{ctx}: pstd");
            assert_eq!(
                a.mean_quality.to_bits(),
                b.mean_quality.to_bits(),
                "{ctx}: mean_quality"
            );
            assert_eq!(a.qstd, b.qstd, "{ctx}: qstd");
            assert_eq!(
                a.mean_mapping_quality.to_bits(),
                b.mean_mapping_quality.to_bits(),
                "{ctx}: mean_mapping_quality"
            );
            assert_eq!(
                a.hq_lq_ratio.to_bits(),
                b.hq_lq_ratio.to_bits(),
                "{ctx}: hq_lq_ratio"
            );
            assert_eq!(a.hi_freq.to_bits(), b.hi_freq.to_bits(), "{ctx}: hi_freq");
            assert_eq!(
                a.extra_freq.to_bits(),
                b.extra_freq.to_bits(),
                "{ctx}: extra_freq"
            );
            assert_eq!(a.shift3, b.shift3, "{ctx}: shift3");
            assert_eq!(a.msi.to_bits(), b.msi.to_bits(), "{ctx}: msi");
            assert_eq!(a.msint, b.msint, "{ctx}: msint");
            assert_eq!(a.nm.to_bits(), b.nm.to_bits(), "{ctx}: nm");
            assert_eq!(a.hicnt, b.hicnt, "{ctx}: hicnt");
            assert_eq!(a.hicov, b.hicov, "{ctx}: hicov");
            assert_eq!(a.left_seq, b.left_seq, "{ctx}: left_seq");
            assert_eq!(a.right_seq, b.right_seq, "{ctx}: right_seq");
            assert_eq!(a.start_pos, b.start_pos, "{ctx}: start_pos");
            assert_eq!(a.end_pos, b.end_pos, "{ctx}: end_pos");
            assert_eq!(a.ref_rev_cov, b.ref_rev_cov, "{ctx}: ref_rev_cov");
            assert_eq!(a.ref_fwd_cov, b.ref_fwd_cov, "{ctx}: ref_fwd_cov");
            assert_eq!(a.total_pos_cov, b.total_pos_cov, "{ctx}: total_pos_cov");
            assert_eq!(a.duprate.to_bits(), b.duprate.to_bits(), "{ctx}: duprate");
            assert_eq!(a.genotype, b.genotype, "{ctx}: genotype");
            assert_eq!(a.var_allele, b.var_allele, "{ctx}: var_allele");
            assert_eq!(a.ref_allele, b.ref_allele, "{ctx}: ref_allele");
            assert_eq!(a.var_type, b.var_type, "{ctx}: var_type");
            assert_eq!(a.debug, b.debug, "{ctx}: debug");
            assert_eq!(a.crispr, b.crispr, "{ctx}: crispr");
            assert_eq!(a.sv, b.sv, "{ctx}: sv");
            assert_eq!(a.null_record, b.null_record, "{ctx}: null_record");
        };

        // Build a somatic fixture with IUPAC R (pos 6) and Y (pos 7), extracnt > 0 ttcov
        // override (pos 5), p - 1 insertion overlay[p] (pos 10 -> pos 11), "I" key (pos 14),
        // cnt == 0 (pos 16), sv.is_some() (pos 18), sub-threshold somatic allele (pos 22),
        // and pure-ref sites (pos 25, 26).
        let mut som_seq = b"ACGTARYTACGTAIATAAAGTACGTACGTA".to_vec();
        som_seq.resize(40, b'A');
        som_seq[29] = b'I'; // 1-based pos 30: hostile FASTA 'I' + BAM2 insertion
                            // 1-based positions: 5='A', 6='R', 7='Y', 10='C', 11='G', 14='I', 16='T', 18='A', 22='A', 25='T', 26='A', 30='I'
        let mut som_ref = Reference::new();
        som_ref.loaded_regions.push(crate::reference::LoadedRegion {
            chr: "chr1".to_string(),
            sequence_start: 1,
            sequence_end: 40,
        });
        som_ref
            .reference_sequences
            .insert("chr1".to_string(), som_seq.clone());
        som_ref.region_sequences.push(som_seq);
        let som_chr_lens = HashMap::from([("chr1".to_string(), 40)]);
        let som_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(40),
            gene: "SOM".to_string(),
            amplicon: None,
        };

        let mut som_vd = VariationData::new();
        som_vd.duprate = 0.125;

        // pos 5 ('A'): cnt.cnt (25) > total_cov (20) with extracnt = 10 -> ttcov = 25
        som_vd.ref_coverage.insert(5, 20);
        let mut vm5 = VariationMap::new();
        let mut r5 = Variation::new();
        r5.cnt = 25;
        r5.fwd = 15;
        r5.rev = 10;
        r5.extracnt = 10;
        r5.hicnt = 22;
        r5.locnt = 3;
        r5.qmean = 775.0;
        r5.mqmean = 1480.0;
        r5.pmean = 525.0;
        r5.nm = 12.5;
        r5.pstd = true;
        r5.qstd = true;
        vm5.map.insert("A".to_string(), r5);
        som_vd.non_insertion_variants.insert(5, vm5);

        // pos 6 ('R'): IUPAC ambiguity -> valid_base 'A', genotype "R/R"
        som_vd.ref_coverage.insert(6, 18);
        let mut vm6 = VariationMap::new();
        let mut r6 = Variation::new();
        r6.cnt = 18;
        r6.fwd = 18;
        r6.rev = 0;
        r6.hicnt = 18;
        r6.qmean = 540.0;
        r6.mqmean = 1080.0;
        r6.pmean = 360.0;
        vm6.map.insert("R".to_string(), r6);
        som_vd.non_insertion_variants.insert(6, vm6);

        // pos 7 ('Y'): IUPAC ambiguity -> valid_base 'C', genotype "Y/Y"
        som_vd.ref_coverage.insert(7, 14);
        let mut vm7 = VariationMap::new();
        let mut r7 = Variation::new();
        r7.cnt = 14;
        r7.fwd = 7;
        r7.rev = 7;
        r7.hicnt = 12;
        r7.locnt = 2;
        r7.qmean = 420.0;
        r7.mqmean = 840.0;
        r7.pmean = 280.0;
        vm7.map.insert("Y".to_string(), r7);
        som_vd.non_insertion_variants.insert(7, vm7);

        // pos 10 ('C') + insertion "+AT" that decrements overlay[11] ('G')
        som_vd.ref_coverage.insert(10, 10);
        som_vd.ref_coverage.insert(11, 40);
        let mut vm10s = VariationMap::new();
        let mut r10s = Variation::new();
        r10s.cnt = 10;
        r10s.fwd = 5;
        r10s.rev = 5;
        r10s.hicnt = 10;
        r10s.qmean = 300.0;
        r10s.mqmean = 600.0;
        r10s.pmean = 200.0;
        vm10s.map.insert("C".to_string(), r10s);
        som_vd.non_insertion_variants.insert(10, vm10s);

        let mut ins10s = VariationMap::new();
        let mut i10s = Variation::new();
        i10s.cnt = 15;
        i10s.fwd = 9;
        i10s.rev = 6;
        i10s.hicnt = 15;
        i10s.qmean = 450.0;
        i10s.mqmean = 900.0;
        i10s.pmean = 300.0;
        ins10s.map.insert("+AT".to_string(), i10s);
        som_vd.insertion_variants.insert(10, ins10s);

        // pos 11 ('G'): pure ref-only, immediately after pos 10 insertion -> MUST stay off CompactOut
        let mut vm11s = VariationMap::new();
        let mut r11s = Variation::new();
        r11s.cnt = 35;
        r11s.fwd = 20;
        r11s.rev = 15;
        r11s.hicnt = 35;
        r11s.qmean = 1050.0;
        r11s.mqmean = 2100.0;
        r11s.pmean = 700.0;
        vm11s.map.insert("G".to_string(), r11s);
        som_vd.non_insertion_variants.insert(11, vm11s);

        // pos 14 ('I'): sole key "I" -> excluded from CompactOut
        som_vd.ref_coverage.insert(14, 8);
        let mut vm14 = VariationMap::new();
        let mut r14 = Variation::new();
        r14.cnt = 8;
        r14.fwd = 4;
        r14.rev = 4;
        r14.hicnt = 8;
        r14.qmean = 240.0;
        r14.mqmean = 480.0;
        r14.pmean = 160.0;
        vm14.map.insert("I".to_string(), r14);
        som_vd.non_insertion_variants.insert(14, vm14);

        // pos 16 ('T'): cnt == 0 -> excluded from CompactOut, yields Vars { ref_variant: None, variants: [] }
        som_vd.ref_coverage.insert(16, 12);
        let mut vm16 = VariationMap::new();
        vm16.map.insert("T".to_string(), Variation::new());
        som_vd.non_insertion_variants.insert(16, vm16);

        // pos 18 ('A'): sv.is_some() -> excluded from CompactOut
        som_vd.ref_coverage.insert(18, 10);
        let mut vm18 = VariationMap::new();
        vm18.sv = Some(Box::new(Sv {
            r#type: "INV".to_string(),
            splits: 2,
            pairs: 3,
            clusters: 1,
        }));
        let mut r18 = Variation::new();
        r18.cnt = 10;
        r18.fwd = 6;
        r18.rev = 4;
        r18.hicnt = 10;
        r18.qmean = 300.0;
        r18.mqmean = 600.0;
        r18.pmean = 200.0;
        vm18.map.insert("A".to_string(), r18);
        som_vd.non_insertion_variants.insert(18, vm18);

        // pos 22 ('A'): sub-threshold non-ref SNV "G" (1 / 100 = 0.01 <= freq 0.05)
        som_vd.ref_coverage.insert(22, 100);
        let mut vm22 = VariationMap::new();
        let mut r22 = Variation::new();
        r22.cnt = 99;
        r22.fwd = 50;
        r22.rev = 49;
        r22.hicnt = 99;
        r22.qmean = 2970.0;
        r22.mqmean = 5940.0;
        r22.pmean = 1980.0;
        vm22.map.insert("A".to_string(), r22);
        let mut alt22 = Variation::new();
        alt22.cnt = 1;
        alt22.fwd = 1;
        alt22.rev = 0;
        alt22.hicnt = 1;
        alt22.qmean = 32.0;
        alt22.mqmean = 60.0;
        alt22.pmean = 25.0;
        vm22.map.insert("G".to_string(), alt22);
        som_vd.non_insertion_variants.insert(22, vm22);

        // pos 25 ('T') & pos 26 ('A'): pure ref-only sites
        for (p, key) in [(25, "T"), (26, "A")] {
            som_vd.ref_coverage.insert(p, 20);
            let mut vmp = VariationMap::new();
            let mut rp = Variation::new();
            rp.cnt = 20;
            rp.fwd = 10;
            rp.rev = 10;
            rp.hicnt = 20;
            rp.qmean = 600.0;
            rp.mqmean = 1200.0;
            rp.pmean = 400.0;
            vmp.map.insert(key.to_string(), rp);
            som_vd.non_insertion_variants.insert(p, vmp);
        }

        // pos 30 (hostile FASTA 'I'): BAM2-only insertion "+GG" with sole non-insertion key "I"
        som_vd.ref_coverage.insert(30, 20);
        let mut vm30 = VariationMap::new();
        let mut r30 = Variation::new();
        r30.cnt = 12;
        r30.fwd = 6;
        r30.rev = 6;
        r30.hicnt = 12;
        r30.qmean = 360.0;
        r30.mqmean = 720.0;
        r30.pmean = 240.0;
        vm30.map.insert("I".to_string(), r30);
        som_vd.non_insertion_variants.insert(30, vm30);
        let mut ins30 = VariationMap::new();
        let mut i30 = Variation::new();
        i30.cnt = 8;
        i30.fwd = 4;
        i30.rev = 4;
        i30.hicnt = 8;
        i30.qmean = 240.0;
        i30.mqmean = 480.0;
        i30.pmean = 160.0;
        ins30.map.insert("+GG".to_string(), i30);
        som_vd.insertion_variants.insert(30, ins30);

        let mut som_cfg = Config::default();
        som_cfg.bam = "tumor.bam|normal.bam".to_string();
        som_cfg.freq = 0.05;
        let som_builder = ToVarsBuilder::new(&som_ref, &som_region, &som_chr_lens, &som_cfg);

        assert_eq!(format!("{:?}", RefSiteMode::KeepAll), "KeepAll");
        assert_eq!(format!("{:?}", RefSiteMode::DropAll), "DropAll");
        assert!(format!("{:?}", RefSiteMode::KeepWhere(&|_| false)).starts_with("KeepWhere"));

        let keep_all =
            som_builder.build_aligned_variants_with_ref_mode(&som_vd, RefSiteMode::KeepAll);

        let mut compact_sites: rustc_hash::FxHashMap<i32, CompactRefSite> =
            rustc_hash::FxHashMap::default();
        assert!(
            format!("{:?}", RefSiteMode::CompactOut(&mut compact_sites)).starts_with("CompactOut")
        );
        let mut compact_aligned = som_builder.build_aligned_variants_with_ref_mode(
            &som_vd,
            RefSiteMode::CompactOut(&mut compact_sites),
        );

        // Pure ref sites (5, 6, 7, 25, 26) must be in compact_sites and omitted from compact_aligned;
        // special sites (10, 11 after p-1 ins, 14 "I", 16 cnt==0, 18 sv, 22 non-ref, 30 FASTA 'I' + ins) stay on standard path.
        for p in [5, 6, 7, 25, 26] {
            assert!(
                compact_sites.contains_key(&p),
                "expected pos {p} in compact_sites"
            );
            assert!(
                !compact_aligned.contains_key(&p),
                "expected pos {p} omitted from compact_aligned"
            );
        }
        for p in [10, 11, 14, 16, 18, 22, 30] {
            assert!(
                !compact_sites.contains_key(&p),
                "pos {p} must be excluded from compact_sites"
            );
            assert!(
                compact_aligned.contains_key(&p),
                "pos {p} must remain in compact_aligned"
            );
        }
        assert!(compact_aligned.get(&16).unwrap().ref_variant.is_none());
        assert!(compact_aligned.get(&16).unwrap().variants.is_empty());

        for (&p, site) in &compact_sites {
            let hydrated = som_builder.hydrate_compact_ref_site(p, site);
            compact_aligned.insert(p, hydrated);
        }
        assert_eq!(compact_aligned.len(), keep_all.len());
        for (&p, expected_vars) in keep_all.iter() {
            let actual_vars = compact_aligned.get(&p).expect("missing hydrated pos");
            assert_eq!(actual_vars.sv, expected_vars.sv, "pos {p} sv");
            assert_eq!(
                actual_vars.variants.len(),
                expected_vars.variants.len(),
                "pos {p} variants len"
            );
            match (&actual_vars.ref_variant, &expected_vars.ref_variant) {
                (Some(a), Some(b)) => assert_variant_eq(a, b, &format!("pos {p} ref_variant")),
                (None, None) => {}
                _ => panic!("pos {p} ref_variant Option mismatch"),
            }
            for (av, bv) in actual_vars
                .variants
                .iter()
                .zip(expected_vars.variants.iter())
            {
                assert_variant_eq(av, bv, &format!("pos {p} alt"));
            }
        }
        // Verify ttcov = cnt.cnt (25) extra_freq (10/25 = 0.4) vs total_pos_cov = 20 at pos 5,
        // and IUPAC valid_base vs raw genotype at pos 6 ('R') & pos 7 ('Y').
        let rv5 = compact_aligned
            .get(&5)
            .unwrap()
            .ref_variant
            .as_ref()
            .unwrap();
        assert_eq!(rv5.total_pos_cov, 20);
        assert_eq!(rv5.extra_freq.to_bits(), (10.0_f64 / 25.0_f64).to_bits());
        let rv6 = compact_aligned
            .get(&6)
            .unwrap()
            .ref_variant
            .as_ref()
            .unwrap();
        assert_eq!(
            (rv6.ref_allele.as_str(), rv6.var_allele.as_str()),
            ("A", "A")
        );
        assert_eq!(rv6.genotype.as_deref(), Some("R/R"));
        let rv7 = compact_aligned
            .get(&7)
            .unwrap()
            .ref_variant
            .as_ref()
            .unwrap();
        assert_eq!(
            (rv7.ref_allele.as_str(), rv7.var_allele.as_str()),
            ("C", "C")
        );
        assert_eq!(rv7.genotype.as_deref(), Some("Y/Y"));

        // KeepWhere & DropAll in somatic mode retain sub-threshold non-ref allele at pos 22.
        let keep_where = som_builder
            .build_aligned_variants_with_ref_mode(&som_vd, RefSiteMode::KeepWhere(&|p| p == 25));
        assert!(keep_where.contains_key(&22));
        assert_eq!(keep_where.get(&22).unwrap().variants.len(), 1);
        assert!(keep_where.contains_key(&25));
        assert!(!keep_where.contains_key(&26));
        assert_variant_eq(
            keep_where.get(&25).unwrap().ref_variant.as_ref().unwrap(),
            keep_all.get(&25).unwrap().ref_variant.as_ref().unwrap(),
            "KeepWhere pos 25",
        );
        assert!(
            keep_where.contains_key(&30),
            "KeepWhere must retain BAM2 insertion at hostile FASTA 'I' site even when pred(30) == false"
        );
        assert_variant_eq(
            keep_where.get(&30).unwrap().variant_by_desc("+GG").unwrap(),
            keep_all.get(&30).unwrap().variant_by_desc("+GG").unwrap(),
            "KeepWhere pos 30 +GG insertion",
        );

        let som_drop_all =
            som_builder.build_aligned_variants_with_ref_mode(&som_vd, RefSiteMode::DropAll);
        assert!(som_drop_all.contains_key(&22));
        assert_eq!(som_drop_all.get(&22).unwrap().variants.len(), 1);
        assert!(!som_drop_all.contains_key(&25));
        assert!(!som_drop_all.contains_key(&26));

        // Single-sample with debug=true still uses DropAll so JavaHashMap order matches Java.
        let mut single_dbg_cfg = Config::default();
        single_dbg_cfg.bam = "sample.bam".to_string();
        single_dbg_cfg.debug = true;
        single_dbg_cfg.freq = 0.005;
        let single_dbg_builder =
            ToVarsBuilder::new(&som_ref, &som_region, &som_chr_lens, &single_dbg_cfg);
        let single_dbg_aligned = single_dbg_builder.build_aligned_variants(&som_vd);
        assert!(!single_dbg_aligned.contains_key(&25));
        assert!(single_dbg_aligned.contains_key(&22));

        // --- Large-contig (i32::MAX) & negative CRISPR (-J i32::MIN) overflow hardening ---
        let max_pos = i32::MAX;
        let near_max = i32::MAX - 5;
        let win_start = (i32::MAX - 100) as u32;
        let win_bases = vec![b'A'; 101]; // covers [i32::MAX - 100 ..= i32::MAX]
        let mut big_ref = Reference::new();
        big_ref.loaded_regions.push(crate::reference::LoadedRegion {
            chr: "chrBig".to_string(),
            sequence_start: win_start,
            sequence_end: i32::MAX as u32,
        });
        big_ref.region_sequences.push(win_bases);
        let huge_chr_len = (i32::MAX as usize).saturating_add(50_000);
        let big_chr_lens = HashMap::from([("chrBig".to_string(), huge_chr_len)]);
        let big_region = Region {
            chr: "chrBig".to_string(),
            start: Pos1(win_start),
            end: Pos1(i32::MAX as u32),
            gene: "BIG".to_string(),
            amplicon: None,
        };

        let mut big_vd = VariationData::new();
        big_vd.ref_coverage.insert(near_max, 30);
        big_vd.ref_coverage.insert(max_pos, 30);

        let mk_cnt = |cnt: i32| {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - v.fwd;
            v.hicnt = cnt;
            v.qmean = 30.0 * cnt as f64;
            v.mqmean = 60.0 * cnt as f64;
            v.pmean = 20.0 * cnt as f64;
            v
        };

        let mut vm_near = VariationMap::new();
        vm_near.map.insert("A".to_string(), mk_cnt(10));
        vm_near.map.insert("-6".to_string(), mk_cnt(8)); // spans [i32::MAX - 5 ..= i32::MAX] -> end_pos == i32::MAX
        vm_near.map.insert("-20".to_string(), mk_cnt(6)); // clamps end_pos == i32::MAX
        vm_near.map.insert("G".to_string(), mk_cnt(6));
        big_vd.non_insertion_variants.insert(near_max, vm_near);

        let mut ins_near = VariationMap::new();
        ins_near.map.insert("+AA".to_string(), mk_cnt(5));
        big_vd.insertion_variants.insert(near_max, ins_near);

        let mut vm_end = VariationMap::new();
        vm_end.map.insert("A".to_string(), mk_cnt(15));
        vm_end.map.insert("C".to_string(), mk_cnt(8));
        big_vd.non_insertion_variants.insert(max_pos, vm_end);

        let mut ins_end = VariationMap::new();
        ins_end.map.insert("+TT".to_string(), mk_cnt(7));
        big_vd.insertion_variants.insert(max_pos, ins_end);

        let mut big_cfg = Config::default();
        big_cfg.three_prime = true;
        big_cfg.crispr_cutting_site = Some(i32::MIN);
        big_cfg.sv_min_len = u32::MAX;
        let big_builder = ToVarsBuilder::new(&big_ref, &big_region, &big_chr_lens, &big_cfg);
        assert_eq!(big_builder.chr_len, i32::MAX);

        let big_aligned = big_builder.build_aligned_variants(&big_vd);
        let vars_near = big_aligned.get(&near_max).expect("vars at i32::MAX - 5");
        let del6 = vars_near.variant_by_desc("-6").expect("-6 at near_max");
        assert_eq!(del6.end_pos, i32::MAX);
        assert_eq!(
            del6.right_seq, "",
            "end_position == i32::MAX must leave right_seq empty via checked_add(1)"
        );
        let del20 = vars_near.variant_by_desc("-20").expect("-20 at near_max");
        assert_eq!(del20.end_pos, i32::MAX);
        assert_eq!(del20.right_seq, "");

        let vars_max = big_aligned.get(&max_pos).expect("vars at i32::MAX");
        let ins_tt = vars_max.variant_by_desc("+TT").expect("+TT at i32::MAX");
        assert_eq!(ins_tt.end_pos, i32::MAX);
        assert_eq!(ins_tt.right_seq, "");
        let snv_c = vars_max.variant_by_desc("C").expect("C at i32::MAX");
        assert_eq!(snv_c.end_pos, i32::MAX);
        assert_eq!(snv_c.right_seq, "");
    }

    #[test]
    fn extracnt_i64_gate_and_utf8_safe_allele_trimming() {
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), vec![b'A'; 100]);
        ref_genome.region_sequences.push(vec![b'A'; 100]);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mut vd = VariationData::new();
        // Negative ref_coverage (i32::MIN) + positive cnt.cnt exercises
        // `i64::from(cnt.cnt) - i64::from(total_pos_coverage) < i64::from(cnt.extracnt)`
        // without i32 subtraction overflow in both build_variant_from_cnt and create_insertion.
        vd.ref_coverage.insert(10, i32::MIN);

        let mk_var = |cnt: i32, extracnt: i32| {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.extracnt = extracnt;
            v.fwd = cnt;
            v.hicnt = cnt;
            v.qmean = 30.0 * (cnt as f64);
            v.mqmean = 60.0 * (cnt as f64);
            v.pmean = 20.0 * (cnt as f64);
            v
        };

        let mut vm10 = VariationMap::new();
        // Complex allele with '&' and non-ASCII leading UTF-8 scalar ('ä' = 2 bytes)
        // exercises in-place '&' removal and UTF-8-safe leading scalar drain.
        vm10.map.insert("ä&CG".to_string(), mk_var(5, 10));
        vm10.map.insert("-2^ä&TT".to_string(), mk_var(4, 10));
        vd.non_insertion_variants.insert(10, vm10);

        let mut ins10 = VariationMap::new();
        ins10.map.insert("+AC&GT".to_string(), mk_var(6, 10));
        vd.insertion_variants.insert(10, ins10);

        let mut cfg = Config::default();
        cfg.min_reads_bias = u32::MAX;
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);
        let aligned = builder.build_aligned_variants(&vd);
        let vars10 = aligned.get(&10).expect("vars at 10");
        assert!(
            !vars10.variants.is_empty(),
            "variants must finalize without UTF-8 boundary or i32 subtraction panic"
        );
        for v in &vars10.variants {
            assert!(!v.var_allele.contains('&'));
        }
    }

    #[test]
    fn build_vars_at_position_matches_full_build_aligned_variants_across_bucket_orders_and_skip_guards(
    ) {
        use crate::variations::Sv;

        let assert_var_fields_eq = |a: &Variant, b: &Variant, ctx: &str| {
            assert_eq!(a.description, b.description, "{ctx}: description");
            assert_eq!(a.pos_coverage, b.pos_coverage, "{ctx}: pos_coverage");
            assert_eq!(a.fwd, b.fwd, "{ctx}: fwd");
            assert_eq!(a.rev, b.rev, "{ctx}: rev");
            assert_eq!(a.bias, b.bias, "{ctx}: bias");
            assert_eq!(
                a.frequency.to_bits(),
                b.frequency.to_bits(),
                "{ctx}: frequency"
            );
            assert_eq!(
                a.mean_position.to_bits(),
                b.mean_position.to_bits(),
                "{ctx}: mean_position"
            );
            assert_eq!(a.pstd, b.pstd, "{ctx}: pstd");
            assert_eq!(
                a.mean_quality.to_bits(),
                b.mean_quality.to_bits(),
                "{ctx}: mean_quality"
            );
            assert_eq!(a.qstd, b.qstd, "{ctx}: qstd");
            assert_eq!(
                a.mean_mapping_quality.to_bits(),
                b.mean_mapping_quality.to_bits(),
                "{ctx}: mean_mapping_quality"
            );
            assert_eq!(
                a.hq_lq_ratio.to_bits(),
                b.hq_lq_ratio.to_bits(),
                "{ctx}: hq_lq_ratio"
            );
            assert_eq!(a.hi_freq.to_bits(), b.hi_freq.to_bits(), "{ctx}: hi_freq");
            assert_eq!(
                a.extra_freq.to_bits(),
                b.extra_freq.to_bits(),
                "{ctx}: extra_freq"
            );
            assert_eq!(a.shift3, b.shift3, "{ctx}: shift3");
            assert_eq!(a.msi.to_bits(), b.msi.to_bits(), "{ctx}: msi");
            assert_eq!(a.msint, b.msint, "{ctx}: msint");
            assert_eq!(a.nm.to_bits(), b.nm.to_bits(), "{ctx}: nm");
            assert_eq!(a.hicnt, b.hicnt, "{ctx}: hicnt");
            assert_eq!(a.hicov, b.hicov, "{ctx}: hicov");
            assert_eq!(a.left_seq, b.left_seq, "{ctx}: left_seq");
            assert_eq!(a.right_seq, b.right_seq, "{ctx}: right_seq");
            assert_eq!(a.start_pos, b.start_pos, "{ctx}: start_pos");
            assert_eq!(a.end_pos, b.end_pos, "{ctx}: end_pos");
            assert_eq!(a.ref_rev_cov, b.ref_rev_cov, "{ctx}: ref_rev_cov");
            assert_eq!(a.ref_fwd_cov, b.ref_fwd_cov, "{ctx}: ref_fwd_cov");
            assert_eq!(a.total_pos_cov, b.total_pos_cov, "{ctx}: total_pos_cov");
            assert_eq!(a.duprate.to_bits(), b.duprate.to_bits(), "{ctx}: duprate");
            assert_eq!(a.genotype, b.genotype, "{ctx}: genotype");
            assert_eq!(a.var_allele, b.var_allele, "{ctx}: var_allele");
            assert_eq!(a.ref_allele, b.ref_allele, "{ctx}: ref_allele");
            assert_eq!(a.var_type, b.var_type, "{ctx}: var_type");
            assert_eq!(a.debug, b.debug, "{ctx}: debug");
            assert_eq!(a.crispr, b.crispr, "{ctx}: crispr");
            assert_eq!(a.sv, b.sv, "{ctx}: sv");
            assert_eq!(a.null_record, b.null_record, "{ctx}: null_record");
        };

        let assert_vars_opt_eq = |actual: &Option<crate::variations::Vars>,
                                  expected: &Option<crate::variations::Vars>,
                                  ctx: &str| {
            match (actual, expected) {
                (None, None) => {}
                (Some(a), Some(b)) => {
                    assert_eq!(a.sv, b.sv, "{ctx}: Vars.sv");
                    assert_eq!(
                        a.var_desc_to_variant.len(),
                        b.var_desc_to_variant.len(),
                        "{ctx}: Vars.var_desc_to_variant.len()"
                    );
                    match (&a.ref_variant, &b.ref_variant) {
                        (None, None) => {}
                        (Some(ra), Some(rb)) => {
                            assert_var_fields_eq(ra, rb, &format!("{ctx}: ref_variant"));
                        }
                        (ra, rb) => {
                            panic!("{ctx}: ref_variant Option mismatch: {ra:?} vs {rb:?}")
                        }
                    }
                    assert_eq!(
                        a.variants.len(),
                        b.variants.len(),
                        "{ctx}: Vars.variants.len()"
                    );
                    for (idx, (va, vb)) in a.variants.iter().zip(b.variants.iter()).enumerate() {
                        assert_var_fields_eq(va, vb, &format!("{ctx}: variants[{idx}]"));
                    }
                }
                (a, b) => panic!("{ctx}: Option<Vars> mismatch: {a:?} vs {b:?}"),
            }
        };

        let mk_cnt = |cnt: i32, fwd: i32, rev: i32| {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.fwd = fwd;
            v.rev = rev;
            v.hicnt = cnt;
            v.pstd = true;
            v.qstd = true;
            v.qmean = 32.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 24.0 * f64::from(cnt);
            v.nm = 1.0 * f64::from(cnt);
            v
        };

        // Reference 1..=120: repeating ACGT, with hostile 'I' at 1-based pos 40.
        let mut seq: Vec<u8> = (0..120).map(|i| b"ACGT"[i % 4]).collect();
        seq[39] = b'I';
        let ref_base_str = |pos: i32| -> String {
            let ch = seq[(pos - 1) as usize] as char;
            ch.to_string()
        };

        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 120,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq.clone());
        let chr_lengths = HashMap::from([("chr1".to_string(), 120)]);

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(10),
            end: Pos1(80),
            gene: "SOM".to_string(),
            amplicon: None,
        };

        let mut vd = VariationData::new();
        vd.duprate = 0.0625;
        // Owning range [12, 80]: in-region pos 11 can test unowned merge_planted skip guard.
        vd.owning_range = Some((Pos1(12), Pos1(80)));

        // 1) Out-of-region prev_pos = region.start - 1 (9) -> region.start (10):
        //    pos 9 has insertion "+AT" with ttcov < cnt.cnt < cov_next - cnt.cnt,
        //    which MUST be skipped because 9 < region.start (10) and sv is None.
        vd.ref_coverage.insert(9, 8);
        let mut vm9 = VariationMap::new();
        vm9.map.insert(ref_base_str(9), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(9, vm9);
        let mut ins9 = VariationMap::new();
        ins9.map.insert("+AT".to_string(), mk_cnt(12, 7, 5));
        vd.insertion_variants.insert(9, ins9);

        vd.ref_coverage.insert(10, 40);
        let mut vm10 = VariationMap::new();
        vm10.map.insert(ref_base_str(10), mk_cnt(32, 16, 16));
        vm10.map.insert("-2".to_string(), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(10, vm10);

        // 2) Unowned merge_planted prev_pos (11 -> 12):
        //    11 >= region.start (10) so it passes region bounds, but 11 < owning_range.start (12).
        vd.merge_planted.insert(11);
        vd.ref_coverage.insert(11, 8);
        let mut vm11 = VariationMap::new();
        vm11.map.insert(ref_base_str(11), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(11, vm11);
        let mut ins11 = VariationMap::new();
        ins11.map.insert("+GG".to_string(), mk_cnt(12, 6, 6));
        vd.insertion_variants.insert(11, ins11);

        vd.ref_coverage.insert(12, 40);
        let mut vm12 = VariationMap::new();
        vm12.map.insert(ref_base_str(12), mk_cnt(30, 15, 15));
        vm12.map.insert("-3".to_string(), mk_cnt(10, 5, 5));
        vd.non_insertion_variants.insert(12, vm12);

        // 3) Adjacent insertions & deletions across (p - 1, p, p + 1) = (13, 14, 15) where
        //    bucket(13) < bucket(14) < bucket(15) at capacity 32/64:
        //    - pos 13 insertion decrements overlay[14] ref fwd/rev;
        //    - pos 14 has deletion "-2" AND insertion "+CC" that decrements overlay[15] (read back
        //      immediately by pos 14's collect_reference_variants AND seen when evaluating pos 15);
        //    - pos 15 has SNV + deletion.
        vd.ref_coverage.insert(13, 8);
        let mut vm13 = VariationMap::new();
        vm13.map.insert(ref_base_str(13), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(13, vm13);
        let mut ins13 = VariationMap::new();
        ins13.map.insert("+ACGT".to_string(), mk_cnt(12, 7, 5));
        vd.insertion_variants.insert(13, ins13);

        vd.ref_coverage.insert(14, 36);
        let mut vm14 = VariationMap::new();
        vm14.map.insert(ref_base_str(14), mk_cnt(28, 15, 13));
        vm14.map.insert("-2".to_string(), mk_cnt(6, 3, 3));
        vd.non_insertion_variants.insert(14, vm14);
        let mut ins14 = VariationMap::new();
        ins14.map.insert("+CC".to_string(), mk_cnt(40, 22, 18));
        vd.insertion_variants.insert(14, ins14);

        vd.ref_coverage.insert(15, 95);
        let mut vm15 = VariationMap::new();
        vm15.map.insert(ref_base_str(15), mk_cnt(80, 45, 35));
        vm15.map.insert("-4".to_string(), mk_cnt(15, 8, 7));
        vd.non_insertion_variants.insert(15, vm15);

        // 4) Reverse bucket order p < (p - 1) at capacity 32:
        //    With <= 24 total keys, JavaHashMap capacity is 32, so bucket(32) = 0 < bucket(31) = 31!
        //    pos 32 is visited BEFORE pos 31 in JavaHashMap::iter(), so pos 31's insertion MUST NOT
        //    decrement overlay[32]!
        vd.ref_coverage.insert(31, 8);
        let mut vm31 = VariationMap::new();
        vm31.map.insert(ref_base_str(31), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(31, vm31);
        let mut ins31 = VariationMap::new();
        ins31.map.insert("+TTAA".to_string(), mk_cnt(12, 7, 5));
        vd.insertion_variants.insert(31, ins31);

        vd.ref_coverage.insert(32, 40);
        let mut vm32 = VariationMap::new();
        vm32.map.insert(ref_base_str(32), mk_cnt(30, 17, 13));
        vm32.map.insert("-2".to_string(), mk_cnt(10, 5, 5));
        vd.non_insertion_variants.insert(32, vm32);

        // 5) total_cov == 0 guard on prev_pos (20 -> 21) even when sv.is_some(), and at target_pos (24):
        vd.ref_coverage.insert(20, 0);
        let mut vm20 = VariationMap::new();
        vm20.sv = Some(Box::new(Sv {
            r#type: "DEL".to_string(),
            splits: 1,
            pairs: 1,
            clusters: 1,
        }));
        vm20.map.insert(ref_base_str(20), mk_cnt(5, 3, 2));
        vd.non_insertion_variants.insert(20, vm20);
        let mut ins20 = VariationMap::new();
        ins20.map.insert("+AA".to_string(), mk_cnt(10, 5, 5));
        vd.insertion_variants.insert(20, ins20);

        vd.ref_coverage.insert(21, 35);
        let mut vm21 = VariationMap::new();
        vm21.map.insert(ref_base_str(21), mk_cnt(28, 14, 14));
        vm21.map.insert("-2".to_string(), mk_cnt(7, 4, 3));
        vd.non_insertion_variants.insert(21, vm21);

        vd.ref_coverage.insert(24, 0);
        let mut vm24 = VariationMap::new();
        vm24.map.insert("-2".to_string(), mk_cnt(5, 3, 2));
        vd.non_insertion_variants.insert(24, vm24);

        // 6) Missing ref_coverage without SV on prev_pos (26 -> 27):
        let mut vm26 = VariationMap::new();
        vm26.map.insert(ref_base_str(26), mk_cnt(6, 3, 3));
        vd.non_insertion_variants.insert(26, vm26);
        let mut ins26 = VariationMap::new();
        ins26.map.insert("+GG".to_string(), mk_cnt(10, 5, 5));
        vd.insertion_variants.insert(26, ins26);

        vd.ref_coverage.insert(27, 35);
        let mut vm27 = VariationMap::new();
        vm27.map.insert(ref_base_str(27), mk_cnt(27, 14, 13));
        vm27.map.insert("-2".to_string(), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(27, vm27);

        // 7) Pure-ref DropAll sites: ordinary pure-ref at pos 35, plus hostile FASTA 'I' at
        //    prev_pos 40 (empty non_insertion_variants[40] + insertion_variants[40], which
        //    is_pure_ref_variation treats as pure-ref under DropAll -> skips overlay[41]!).
        vd.ref_coverage.insert(35, 25);
        let mut vm35 = VariationMap::new();
        vm35.map.insert(ref_base_str(35), mk_cnt(25, 13, 12));
        vd.non_insertion_variants.insert(35, vm35);

        vd.ref_coverage.insert(40, 8);
        vd.non_insertion_variants.insert(40, VariationMap::new());
        let mut ins40 = VariationMap::new();
        ins40.map.insert("+GGGG".to_string(), mk_cnt(12, 7, 5));
        vd.insertion_variants.insert(40, ins40);

        vd.ref_coverage.insert(41, 40);
        let mut vm41 = VariationMap::new();
        vm41.map.insert(ref_base_str(41), mk_cnt(32, 16, 16));
        vm41.map.insert("-2".to_string(), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(41, vm41);

        // 8) Sub-threshold non-ref site at pos 50 (1 / 100 = 0.01 <= freq 0.05) to test
        //    somatic Keep/DropAll retention vs single-sample maxfreq <= config.freq removal.
        vd.ref_coverage.insert(50, 100);
        let mut vm50 = VariationMap::new();
        vm50.map.insert(ref_base_str(50), mk_cnt(99, 50, 49));
        vm50.map.insert("-2".to_string(), mk_cnt(1, 1, 0));
        vd.non_insertion_variants.insert(50, vm50);

        // Confirm capacity is 32 so (13 < 14 < 15) holds AND (32 < 31) reverses order!
        assert_eq!(vd.non_insertion_variants.capacity(), 32);
        let key_order: Vec<i32> = vd.non_insertion_variants.keys().copied().collect();
        let idx_of = |k: i32| key_order.iter().position(|&x| x == k).unwrap();
        assert!(idx_of(13) < idx_of(14) && idx_of(14) < idx_of(15));
        assert!(
            idx_of(32) < idx_of(31),
            "expected pos 32 before pos 31 in cap-32 JavaHashMap"
        );

        let mut som_cfg = Config::default();
        som_cfg.bam = "tumor.bam|normal.bam".to_string();
        som_cfg.freq = 0.05;
        som_cfg.debug = true;
        let som_builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &som_cfg);

        let full_drop = som_builder.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::DropAll);
        let full_keep = som_builder.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::KeepAll);

        for pos in (0..=85).chain([i32::MIN, -1, i32::MAX]) {
            let single_drop =
                som_builder.build_vars_at_position_with_ref_mode(&vd, pos, RefSiteMode::DropAll);
            let expected_drop = full_drop.get(&pos).cloned();
            assert_vars_opt_eq(
                &single_drop,
                &expected_drop,
                &format!("som DropAll pos={pos}"),
            );

            let single_keep = som_builder.build_vars_at_position(&vd, pos);
            let expected_keep = full_keep.get(&pos).cloned();
            assert_vars_opt_eq(
                &single_keep,
                &expected_keep,
                &format!("som KeepAll pos={pos}"),
            );
        }

        // Single-sample mode (no '|' in bam, debug = false): exercises maxfreq <= config.freq removal at pos 50.
        let mut single_cfg = Config::default();
        single_cfg.bam = "combined.bam".to_string();
        single_cfg.freq = 0.05;
        let single_builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &single_cfg);
        let single_full =
            single_builder.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::DropAll);
        assert!(!single_full.contains_key(&50));
        for pos in (0..=85).chain([i32::MIN, -1, i32::MAX]) {
            let got =
                single_builder.build_vars_at_position_with_ref_mode(&vd, pos, RefSiteMode::DropAll);
            let exp = single_full.get(&pos).cloned();
            assert_vars_opt_eq(&got, &exp, &format!("single DropAll pos={pos}"));
        }
    }

    #[test]
    fn tovars_rejects_negative_positions_and_clamps_join_ref_into_to_loaded_span() {
        let mut ref_genome = Reference::new();
        ref_genome.region_sequences.push(b"ACGTACGTAC".to_vec());
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 10,
            });
        // Wide region with `end = u32::MAX` so `(-5_i32 as u32)` (`4294967291`) would have
        // satisfied `(pos as u32) >= region.start.0 && (pos as u32) <= region.end.0` before AR-4!
        let wide_region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(u32::MAX),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut chr_lengths = std::collections::HashMap::new();
        chr_lengths.insert("chr1".to_string(), 10usize);
        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &wide_region, &chr_lengths, &cfg);

        let mut vd = VariationData::default();
        vd.ref_coverage.insert(-5, 20);
        let mut vmap = crate::variations::VariationMap::default();
        let mut v = crate::variations::Variation::default();
        v.cnt = 10;
        v.hicnt = 10;
        v.pmean = 100.0;
        v.qmean = 350.0;
        v.qstd = true;
        v.pstd = true;
        vmap.map.insert("G".to_string(), v);
        vd.non_insertion_variants.insert(-5, vmap);

        assert!(
            builder.build_vars_at_position(&vd, -5).is_none(),
            "negative i32 position -5 must be rejected even when region.end == u32::MAX"
        );

        // join_ref_into over a 200 Mbp range when only 1..=10 is loaded must clamp to 1..=10
        // in O(loaded_span) without reserving 1 MiB or looping 200,000,000 times.
        crate::reference::reset_get_base_calls();
        let t0 = std::time::Instant::now();
        let mut out = String::new();
        builder.join_ref_into(&mut out, 1, 200_000_000);
        let calls = crate::reference::get_base_calls();
        assert_eq!(out, "ACGTACGTAC");
        assert!(
            out.capacity() <= 64,
            "join_ref_into must reserve clamped loaded span capacity, got {}",
            out.capacity()
        );
        assert!(
            calls <= 16,
            "join_ref_into over 200 Mbp with 1..=10 loaded must visit <= 16 bases, got {calls}"
        );
        assert!(
            t0.elapsed() < std::time::Duration::from_secs(10),
            "join_ref_into over 200 Mbp hung: {:?}",
            t0.elapsed()
        );
    }

    #[test]
    fn join_ref_into_and_complex_del_skip_disjoint_200mbp_unloaded_gap() {
        let mut ref_genome = Reference::new();
        let seq1: Vec<u8> = b"ACGT".repeat(325); // 1..=1_300 (1,300 bp)
        let seq2: Vec<u8> = b"TGCA".repeat(200).into_iter().chain([b'G']).collect(); // 199_999_800..=200_000_600 (801 bp)
        ref_genome.region_sequences.push(seq1.clone());
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 1_300,
            });
        ref_genome.region_sequences.push(seq2.clone());
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 199_999_800,
                sequence_end: 200_000_600,
            });

        let region = Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(1_000),
            gene: "G".to_string(),
            amplicon: None,
        };
        let mut chr_lengths = std::collections::HashMap::new();
        chr_lengths.insert("chr1".to_string(), 250_000_000usize);
        let cfg = Config::default();
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);

        // 1. Direct `join_ref_into` across a 200 Mbp unloaded gap between two `LoadedRegion`s.
        let mut expected = String::from_utf8(seq1).unwrap();
        expected.push_str(std::str::from_utf8(&seq2[..201]).unwrap());
        assert_eq!(expected.len(), 1_501);

        crate::reference::reset_get_base_calls();
        let t0 = std::time::Instant::now();
        let mut out = String::new();
        builder.join_ref_into(&mut out, 1, 200_000_000);
        let elapsed = t0.elapsed();
        let join_calls = crate::reference::get_base_calls();

        assert_eq!(out, expected);
        assert!(
            out.capacity() <= 4_096,
            "join_ref_into across 200 Mbp unloaded gap must reserve merged loaded span, got {}",
            out.capacity()
        );
        assert!(
            join_calls <= 4_096,
            "join_ref_into across 200 Mbp unloaded gap must query <= 4,096 loaded bases, got {join_calls}"
        );
        assert!(
            elapsed < std::time::Duration::from_secs(10),
            "join_ref_into across 200 Mbp unloaded gap hung: {elapsed:?}"
        );

        // 2. End-to-end reachability under default CLI flags (`-L 1000 -I 50 -k 1`):
        //    `+A#C^200000000` at position 10 invokes `tovars.rs:1524` `join_ref_into` up to
        //    `end_position + 200_000_000` without being gated by `sv_min_len`.
        let mut vd = VariationData::default();
        vd.ref_coverage.insert(10, 20);
        let mut vmap = crate::variations::VariationMap::default();
        let mut v = crate::variations::Variation::default();
        v.cnt = 10;
        v.hicnt = 10;
        v.pmean = 200.0;
        v.qmean = 350.0;
        v.mqmean = 600.0;
        v.qstd = true;
        v.pstd = true;
        vmap.map.insert("+A#C^200000000".to_string(), v);
        vd.non_insertion_variants.insert(10, vmap);

        crate::reference::reset_get_base_calls();
        let t1 = std::time::Instant::now();
        let built = builder
            .build_vars_at_position(&vd, 10)
            .expect("expected variant at pos 10");
        let elapsed_e2e = t1.elapsed();
        let e2e_calls = crate::reference::get_base_calls();
        assert_eq!(built.variants.len(), 1);
        assert!(
            e2e_calls < 10_000,
            "build_vars_at_position for +A#C^200000000 across 200 Mbp gap must query < 10,000 bases, got {e2e_calls}"
        );
        assert!(
            elapsed_e2e < std::time::Duration::from_secs(10),
            "build_vars_at_position for +A#C^200000000 across 200 Mbp gap hung: {elapsed_e2e:?}"
        );
    }

    #[test]
    fn test_build_aligned_variants_pruned_matches_postprocess_with_del_amp_and_somatic_bam2() {
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let ref_base_str = |pos: i32| -> String { (seq[(pos - 1) as usize] as char).to_string() };

        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 200,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq.clone());
        let chr_lengths = HashMap::from([("chr1".to_string(), 200usize)]);

        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(200),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32, hicnt: i32, fwd: i32, rev: i32, pos_val: f64, qual_val: f64| {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = hicnt;
            v.locnt = cnt - hicnt;
            v.fwd = fwd;
            v.rev = rev;
            v.pmean = pos_val * f64::from(cnt);
            v.qmean = qual_val * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.nm = 1.0 * f64::from(cnt);
            v.pstd = true;
            v.qstd = true;
            v
        };

        // --- 1. Single-sample (!pileup, !debug) ---
        let mut vd = VariationData::new();
        vd.duprate = 0.05;

        // pos 10: doomed by hicnt = 1 < min_reads = 2 (freq = 1/20 = 0.05 > config.freq = 0.01)
        vd.ref_coverage.insert(10, 20);
        let mut vm10 = VariationMap::new();
        vm10.map
            .insert(ref_base_str(10), mk_var(19, 19, 10, 9, 20.0, 35.0));
        vm10.map
            .insert("A".to_string(), mk_var(1, 1, 1, 0, 20.0, 35.0));
        vd.non_insertion_variants.insert(10, vm10);

        // pos 20: doomed by mean_position = 3.0 < config.read_pos = 5.0
        vd.ref_coverage.insert(20, 20);
        let mut vm20 = VariationMap::new();
        vm20.map
            .insert(ref_base_str(20), mk_var(17, 17, 9, 8, 20.0, 35.0));
        vm20.map
            .insert("A".to_string(), mk_var(3, 3, 2, 1, 3.0, 35.0));
        vd.non_insertion_variants.insert(20, vm20);

        // pos 30: doomed by mean_quality = 20.0 < config.base_qual = 22.5
        vd.ref_coverage.insert(30, 20);
        let mut vm30 = VariationMap::new();
        vm30.map
            .insert(ref_base_str(30), mk_var(17, 17, 9, 8, 20.0, 35.0));
        vm30.map
            .insert("A".to_string(), mk_var(3, 3, 2, 1, 20.0, 20.0));
        vd.non_insertion_variants.insert(30, vm30);

        // pos 50: long >30-byte MNP/Complex description with &ATGC extension (AR-12) + +len insertion
        vd.ref_coverage.insert(50, 20);
        let mut vm50 = VariationMap::new();
        let long_complex = "G#ACGTACGTACGTACGTACGTACGTACGTACGT&TT";
        assert!(long_complex.len() > 30);
        vm50.map
            .insert(long_complex.to_string(), mk_var(12, 12, 6, 6, 22.0, 35.0));
        vd.non_insertion_variants.insert(50, vm50);
        let mut ins50 = VariationMap::new();
        ins50
            .map
            .insert("+ACGTACGTACGT".to_string(), mk_var(8, 8, 4, 4, 22.0, 35.0));
        vd.insertion_variants.insert(50, ins50);

        // pos 55: ref_var with strand bias ("1") + +len insertion ("+12" / "+4") (AR-04, AR-12)
        vd.ref_coverage.insert(55, 32);
        let mut vm55 = VariationMap::new();
        vm55.map
            .insert(ref_base_str(55), mk_var(20, 20, 20, 0, 22.0, 35.0));
        vd.non_insertion_variants.insert(55, vm55);
        let mut ins55 = VariationMap::new();
        ins55
            .map
            .insert("+TTAA".to_string(), mk_var(12, 12, 6, 6, 22.0, 35.0));
        vd.insertion_variants.insert(55, ins55);

        // pos 58: pure insertion site (no ref_var) -> exercises gt1_buf ("+12") + gt2_buf ("+12")
        vd.ref_coverage.insert(58, 12);
        vd.non_insertion_variants.insert(58, VariationMap::new());
        let mut ins58 = VariationMap::new();
        ins58.map.insert(
            "+ACGTACGTACGT".to_string(),
            mk_var(12, 12, 6, 6, 22.0, 35.0),
        );
        vd.insertion_variants.insert(58, ins58);

        // pos 80: <DEL>&ATGC candidate (hicnt = 3 >= min_reads) whose initial freq = 3/400 = 0.0075
        // is below config.freq (0.01) before tovars.rs:1540 recomputes it against ref_coverage[79] = 15
        // to 3/15 = 0.20 >= config.freq, paired with a doomed hicnt = 1 allele so maxfreq = 5/400 > 0.01.
        vd.ref_coverage.insert(79, 15);
        vd.ref_coverage.insert(80, 400);
        let mut vm80 = VariationMap::new();
        vm80.map
            .insert(ref_base_str(80), mk_var(392, 392, 196, 196, 20.0, 35.0));
        vm80.map
            .insert("A".to_string(), mk_var(5, 1, 3, 2, 20.0, 25.0));
        vm80.map
            .insert("-60&ATGC".to_string(), mk_var(3, 3, 2, 1, 20.0, 35.0));
        vd.non_insertion_variants.insert(80, vm80);

        let mut single_cfg = Config::default();
        single_cfg.bam = "sample.bam".to_string();
        single_cfg.sample = "S1".to_string();
        single_cfg.freq = 0.01;
        single_cfg.min_reads = 2;
        single_cfg.read_pos = 5.0;
        single_cfg.base_qual = 22.5;
        single_cfg.sv_min_len = 50;
        single_cfg.pileup = false;
        single_cfg.debug = false;
        single_cfg.emit_intermediate = true;

        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &single_cfg);
        let unpruned = builder.build_aligned_variants(&vd);
        let pruned = builder.build_aligned_variants_pruned(&vd, RefSiteMode::DropAll, &|_| true);

        let keys_unpruned: Vec<i32> = unpruned.keys().copied().collect();
        let keys_pruned: Vec<i32> = pruned.keys().copied().collect();
        assert_eq!(
            keys_pruned, keys_unpruned,
            "JavaHashMap key iteration order must be identical"
        );

        for doomed_pos in [10, 20, 30] {
            assert!(
                !unpruned.get(&doomed_pos).unwrap().variants.is_empty(),
                "unpruned pos {doomed_pos} should carry unpruned variant"
            );
            let doomed_entry = pruned.get(&doomed_pos).unwrap();
            assert!(
                doomed_entry.variants.is_empty(),
                "pruned pos {doomed_pos} must have empty variants while keeping map key"
            );
            assert_eq!(
                doomed_entry.variants.capacity(),
                0,
                "doomed-pruned pos {doomed_pos} must elide entry.variants allocation"
            );
            assert!(
                doomed_entry.var_desc_to_variant.is_empty(),
                "doomed-pruned pos {doomed_pos} must keep var_desc_to_variant empty"
            );
        }

        let vars50 = pruned.get(&50).expect("pos 50");
        let complex50 = vars50
            .variant_by_desc(long_complex)
            .expect("long complex at pos 50");
        assert!(
            complex50.genotype.as_ref().unwrap().len() > 60,
            "long >30-byte MNP/Complex genotype must not truncate: {:?}",
            complex50.genotype
        );
        let vars55 = pruned.get(&55).expect("pos 55");
        let ins55_var = vars55.variant_by_desc("+TTAA").expect("+TTAA at pos 55");
        assert_eq!(ins55_var.bias, "1;2");
        assert_eq!(ins55_var.genotype.as_deref(), Some("G/+4"));
        let vars58 = pruned.get(&58).expect("pos 58");
        let ins58_var = vars58
            .variant_by_desc("+ACGTACGTACGT")
            .expect("+ACGTACGTACGT at pos 58");
        assert_eq!(ins58_var.genotype.as_deref(), Some("+12/+12"));

        let vars80 = pruned.get(&80).expect("pos 80");
        let del80 = vars80
            .variant_by_desc("-60&ATGC")
            .expect("<DEL>&ATGC at pos 80 must survive pruning");
        assert_eq!(del80.var_allele, "<DEL>");
        assert_eq!(del80.frequency, 0.2);

        let splice = rustc_hash::FxHashSet::default();
        for emit_inter in [true, false] {
            let mut post_cfg = single_cfg.clone();
            post_cfg.emit_intermediate = emit_inter;
            let out_unpruned = crate::postprocess::simple::simple_post_process(
                &post_cfg, &region, &unpruned, &splice,
            );
            let out_pruned = crate::postprocess::simple::simple_post_process(
                &post_cfg, &region, &pruned, &splice,
            );
            assert!(!out_pruned.is_empty());
            let mut buf_unpruned = String::new();
            let mut buf_pruned = String::new();
            for v in &out_unpruned {
                crate::output::table::simple_36_append(
                    &mut buf_unpruned,
                    v,
                    &region,
                    &post_cfg.sample,
                    &v.sv,
                );
                buf_unpruned.push('\n');
            }
            for v in &out_pruned {
                crate::output::table::simple_36_append(
                    &mut buf_pruned,
                    v,
                    &region,
                    &post_cfg.sample,
                    &v.sv,
                );
                buf_pruned.push('\n');
            }
            assert_eq!(
                buf_pruned, buf_unpruned,
                "simple_post_process byte parity failed for emit_intermediate={emit_inter}"
            );
        }

        // Guard check: pileup, debug, or amplicon disables doomed pruning even if pred returns true.
        for mode_idx in 0..3 {
            let mut guard_cfg = single_cfg.clone();
            match mode_idx {
                0 => guard_cfg.pileup = true,
                1 => guard_cfg.debug = true,
                _ => guard_cfg.amplicon = Some("10:0.95".to_string()),
            }
            let guard_builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &guard_cfg);
            let guard_aligned =
                guard_builder.build_aligned_variants_pruned(&vd, RefSiteMode::DropAll, &|_| true);
            assert!(
                !guard_aligned.get(&10).unwrap().variants.is_empty(),
                "doomed pruning must be disabled for mode_idx={mode_idx}"
            );
        }

        // --- 2. Somatic BAM1 (CompactOut) + BAM2 (KeepWhere + pruned) ---
        let mut som_cfg = Config::default();
        som_cfg.bam = "tumor.bam|normal.bam".to_string();
        som_cfg.sample = "TUMOR|NORMAL".to_string();
        som_cfg.freq = 0.01;
        som_cfg.min_reads = 2;
        som_cfg.read_pos = 5.0;
        som_cfg.base_qual = 22.5;
        som_cfg.pileup = false;
        som_cfg.debug = false;
        som_cfg.emit_intermediate = true;

        let som_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "SOM".to_string(),
            amplicon: None,
        };

        let mut vd1 = VariationData::new();
        vd1.duprate = 0.02;
        let mut vd2 = VariationData::new();
        vd2.duprate = 0.03;

        // pos 25 (ref 'A'): Tumor has cnt = 8 variant "G", Normal has cnt = 1 (hicnt = 1 < min_reads) "G"
        // -> Tumor has !variants.is_empty(), so Normal at pos 25 MUST NOT be pruned!
        vd1.ref_coverage.insert(25, 20);
        let mut vm1_25 = VariationMap::new();
        vm1_25
            .map
            .insert(ref_base_str(25), mk_var(12, 12, 6, 6, 20.0, 35.0));
        vm1_25
            .map
            .insert("G".to_string(), mk_var(8, 8, 4, 4, 20.0, 35.0));
        vd1.non_insertion_variants.insert(25, vm1_25);

        vd2.ref_coverage.insert(25, 20);
        let mut vm2_25 = VariationMap::new();
        vm2_25
            .map
            .insert(ref_base_str(25), mk_var(19, 19, 10, 9, 20.0, 35.0));
        vm2_25
            .map
            .insert("G".to_string(), mk_var(1, 1, 1, 0, 20.0, 35.0));
        vd2.non_insertion_variants.insert(25, vm2_25);

        // pos 40 (ref 'T'): Tumor is pure-ref (cnt = 20), Normal has cnt = 1 error "A" (hicnt = 1 < min_reads)
        // -> BAM2 prunes pos 40 so pos 40 does NOT hydrate ref_sites1!
        vd1.ref_coverage.insert(40, 20);
        let mut vm1_40 = VariationMap::new();
        vm1_40
            .map
            .insert(ref_base_str(40), mk_var(20, 20, 10, 10, 20.0, 35.0));
        vd1.non_insertion_variants.insert(40, vm1_40);

        vd2.ref_coverage.insert(40, 20);
        let mut vm2_40 = VariationMap::new();
        vm2_40
            .map
            .insert(ref_base_str(40), mk_var(19, 19, 10, 9, 20.0, 35.0));
        vm2_40
            .map
            .insert("A".to_string(), mk_var(1, 1, 1, 0, 20.0, 35.0));
        vd2.non_insertion_variants.insert(40, vm2_40);

        // pos 60 (ref 'T'): Tumor is pure-ref (cnt = 20), Normal has cnt = 3 StrongLOH "G" (hicnt = 3 >= min_reads)
        // -> BAM2 keeps pos 60, hydrates ref_sites1[60], and emits StrongLOH!
        vd1.ref_coverage.insert(60, 20);
        let mut vm1_60 = VariationMap::new();
        vm1_60
            .map
            .insert(ref_base_str(60), mk_var(20, 20, 10, 10, 20.0, 35.0));
        vd1.non_insertion_variants.insert(60, vm1_60);

        vd2.ref_coverage.insert(60, 20);
        let mut vm2_60 = VariationMap::new();
        vm2_60
            .map
            .insert(ref_base_str(60), mk_var(17, 17, 9, 8, 20.0, 35.0));
        vm2_60
            .map
            .insert("G".to_string(), mk_var(3, 3, 2, 1, 20.0, 35.0));
        vd2.non_insertion_variants.insert(60, vm2_60);

        let som_builder = ToVarsBuilder::new(&ref_genome, &som_region, &chr_lengths, &som_cfg);
        let keep_all_1 =
            som_builder.build_aligned_variants_with_ref_mode(&vd1, RefSiteMode::KeepAll);
        let keep_all_2 =
            som_builder.build_aligned_variants_with_ref_mode(&vd2, RefSiteMode::KeepAll);

        let mut ref_sites1: rustc_hash::FxHashMap<i32, CompactRefSite> =
            rustc_hash::FxHashMap::default();
        let mut aligned1 = som_builder
            .build_aligned_variants_with_ref_mode(&vd1, RefSiteMode::CompactOut(&mut ref_sites1));
        assert!(ref_sites1.contains_key(&40));
        assert!(ref_sites1.contains_key(&60));

        let aligned2 = som_builder.build_aligned_variants_pruned(
            &vd2,
            RefSiteMode::KeepWhere(&|p| aligned1.get(&p).is_some_and(|v1| !v1.variants.is_empty())),
            &|p| aligned1.get(&p).is_none_or(|v1| v1.variants.is_empty()),
        );
        assert!(
            aligned2.get(&40).unwrap().variants.is_empty(),
            "BAM2 doomed error at pos 40 must be pruned"
        );
        assert!(
            !aligned2.get(&25).unwrap().variants.is_empty(),
            "BAM2 sub-min_reads allele at Tumor variant pos 25 must NOT be pruned"
        );
        assert!(
            !aligned2.get(&60).unwrap().variants.is_empty(),
            "BAM2 StrongLOH at pos 60 must survive"
        );

        for (&pos, v2) in aligned2.iter() {
            if !v2.variants.is_empty() && !aligned1.contains_key(&pos) {
                if let Some(site) = ref_sites1.get(&pos) {
                    aligned1.insert(pos, som_builder.hydrate_compact_ref_site(pos, site));
                }
            }
        }
        assert!(
            !aligned1.contains_key(&40),
            "pos = 40 must NOT hydrate ref_sites1"
        );
        assert!(
            aligned1.contains_key(&60),
            "pos = 60 StrongLOH MUST hydrate ref_sites1"
        );

        let som_keep_all = crate::postprocess::somatic::somatic_post_process_at(
            &som_cfg,
            &som_region,
            &keep_all_1,
            &keep_all_2,
            &splice,
            None,
            100,
        );
        let som_pruned = crate::postprocess::somatic::somatic_post_process_at(
            &som_cfg,
            &som_region,
            &aligned1,
            &aligned2,
            &splice,
            None,
            100,
        );
        assert_eq!(som_pruned.len(), 2);

        let render_som = |rows: &[crate::postprocess::somatic::SomaticVariant]| -> String {
            let mut out = String::new();
            for sv in rows {
                crate::output::table::somatic_55_append(
                    &mut out,
                    &sv.variant1,
                    sv.variant2.as_ref(),
                    &som_region,
                    &som_cfg.sample,
                    &sv.status,
                    &sv.sv1,
                    &sv.sv2,
                    sv.end_is_variant2,
                );
                out.push('\n');
            }
            out
        };
        assert_eq!(
            render_som(&som_pruned),
            render_som(&som_keep_all),
            "somatic_post_process_at output must be byte-identical to KeepAll"
        );
        assert_eq!(
            aligned2.get(&40).unwrap().variants.capacity(),
            0,
            "somatic BAM2 doomed-pruned pos 40 must have capacity() == 0"
        );
        assert!(aligned2.get(&40).unwrap().var_desc_to_variant.is_empty());
        assert!(aligned2.get(&25).unwrap().var_desc_to_variant.is_empty());
        assert!(aligned2.get(&60).unwrap().var_desc_to_variant.is_empty());

        // --- 3. Variant::empty_internal() + UB-3 12..16 resize threshold sweep (AR-02) ---
        assert_eq!(Variant::new().bias, "0");
        let empty_internal = Variant::empty_internal();
        assert!(empty_internal.bias.is_empty());
        assert_eq!(empty_internal.bias.capacity(), 0);

        // Build a 16-position fixture crossing the 12-entry (16 * 0.75) resize threshold
        // where capacity 16 vs 32 permutes bucket order (e.g. pos 32 vs pos 16, pos 17 vs pos 1).
        let mut ub3_vd = VariationData::new();
        ub3_vd.duprate = 0.01;
        // 11 doomed non-ref positions (hicnt = 1 < min_reads = 2, freq = 1/20 = 0.05 > 0.01):
        // keys chosen to collide differently at cap 16 vs cap 32.
        let doomed_keys = [1, 2, 3, 5, 6, 7, 17, 18, 19, 32, 33];
        for &p in &doomed_keys {
            ub3_vd.ref_coverage.insert(p, 20);
            let mut vm = VariationMap::new();
            vm.map
                .insert(ref_base_str(p), mk_var(19, 19, 10, 9, 20.0, 35.0));
            let alt = if ref_base_str(p) == "A" { "G" } else { "A" };
            vm.map
                .insert(alt.to_string(), mk_var(1, 1, 1, 0, 20.0, 35.0));
            ub3_vd.non_insertion_variants.insert(p, vm);
        }
        // 12th entry (pos 16): surviving variant -> brings map size to 12 (right at cap-16 threshold!).
        ub3_vd.ref_coverage.insert(16, 20);
        let mut vm_16 = VariationMap::new();
        vm_16
            .map
            .insert(ref_base_str(16), mk_var(14, 14, 7, 7, 22.0, 35.0));
        let alt16 = if ref_base_str(16) == "C" { "T" } else { "C" };
        vm_16
            .map
            .insert(alt16.to_string(), mk_var(6, 6, 3, 3, 22.0, 35.0));
        ub3_vd.non_insertion_variants.insert(16, vm_16);

        // Pos 24 (bucket 24 > buckets 0..=19 of the 12 kept entries at cap 32): cnt == 0 -> empty `var`!
        // Visited when `aligned_variants.len() == 12` (exact cap-16 load-factor threshold):
        // in single-sample `should_drop_low_freq(0.0) == true`, an empty `var` MUST NOT
        // call `get_or_put_vars` and therefore MUST NOT trigger the 12->13 resize!
        ub3_vd.ref_coverage.insert(24, 20);
        let mut vm_24 = VariationMap::new();
        vm_24.map.insert("G".to_string(), Variation::new()); // cnt == 0 -> var.is_empty()
        ub3_vd.non_insertion_variants.insert(24, vm_24);

        // First test BEFORE adding the 13th non-empty low-freq entry: with 12 non-empty entries
        // + empty-var pos 24 (13 keys in ub3_vd -> cap 32, pos 24 visited last at len == 12),
        // single-sample capacity MUST stay 16 (no spurious resize!).
        let ub3_12_unpruned = builder.build_aligned_variants(&ub3_vd);
        let ub3_12_pruned =
            builder.build_aligned_variants_pruned(&ub3_vd, RefSiteMode::DropAll, &|_| true);
        assert_eq!(ub3_12_unpruned.len(), 12);
        assert_eq!(ub3_12_unpruned.capacity(), 16);
        assert_eq!(ub3_12_pruned.capacity(), 16);
        assert_eq!(
            ub3_12_pruned.keys().copied().collect::<Vec<_>>(),
            ub3_12_unpruned.keys().copied().collect::<Vec<_>>(),
            "empty var at threshold=12 must not trigger resize"
        );

        // Now add pos 25 (bucket 25 > bucket 24) with `!var.is_empty()` (cnt = 1 / 200 -> freq = 0.005 <= 0.01)
        // so `should_drop_low_freq` drops it AFTER `get_or_put_vars` inserts the 13th entry and
        // resizes `aligned_variants` from 16 -> 32 buckets (keeping len == 12, capacity == 32!).
        ub3_vd.ref_coverage.insert(25, 200);
        let mut vm_25 = VariationMap::new();
        vm_25
            .map
            .insert(ref_base_str(25), mk_var(199, 199, 100, 99, 20.0, 35.0));
        let alt25 = if ref_base_str(25) == "A" { "G" } else { "A" };
        vm_25
            .map
            .insert(alt25.to_string(), mk_var(1, 1, 1, 0, 20.0, 35.0));
        ub3_vd.non_insertion_variants.insert(25, vm_25);

        // Also add pure-ref (`non_ref_rows == 0`) sites at pos 49 and pos 64 for KeepAll / somatic BAM2.
        for &rp in &[49, 64] {
            ub3_vd.ref_coverage.insert(rp, 24);
            let mut vmr = VariationMap::new();
            vmr.map
                .insert(ref_base_str(rp), mk_var(24, 24, 12, 12, 21.0, 35.0));
            ub3_vd.non_insertion_variants.insert(rp, vmr);
        }

        let ub3_unpruned = builder.build_aligned_variants(&ub3_vd);
        let ub3_pruned =
            builder.build_aligned_variants_pruned(&ub3_vd, RefSiteMode::DropAll, &|_| true);
        assert_eq!(ub3_unpruned.len(), 12);
        assert_eq!(
            ub3_unpruned.capacity(),
            32,
            "13th transient low-freq !var.is_empty() entry must resize 16 -> 32 before remove"
        );
        assert_eq!(ub3_pruned.len(), 12);
        assert_eq!(ub3_pruned.capacity(), 32);
        assert_ne!(
            ub3_pruned.keys().copied().collect::<Vec<_>>(),
            ub3_12_pruned.keys().copied().collect::<Vec<_>>(),
            "16 -> 32 resize must permute bucket order vs pre-resize cap-16 map"
        );
        assert_eq!(
            ub3_pruned.keys().copied().collect::<Vec<_>>(),
            ub3_unpruned.keys().copied().collect::<Vec<_>>(),
            "UB-3 key order across 12->13 transient resize must match unpruned exactly"
        );
        for &dk in &doomed_keys {
            let e = ub3_pruned.get(&dk).unwrap();
            assert!(e.variants.is_empty());
            assert_eq!(e.variants.capacity(), 0);
            assert!(e.var_desc_to_variant.is_empty());
        }
        let surv16 = ub3_pruned.get(&16).unwrap();
        let exp16 = ub3_unpruned.get(&16).unwrap();
        assert_eq!(surv16.variants.len(), 1);
        assert!(surv16.var_desc_to_variant.is_empty());
        assert_eq!(surv16.variants[0].start_pos, 16);
        assert_eq!(surv16.variants[0].genotype, exp16.variants[0].genotype);
        assert_eq!(
            surv16.ref_variant.as_ref().unwrap().pos_coverage,
            exp16.ref_variant.as_ref().unwrap().pos_coverage,
            "surviving pos 16 ref_variant must match unpruned"
        );

        // Somatic BAM2 `KeepAll` + `build_aligned_variants_pruned` with `prune_doomed_at`:
        // - `non_ref_rows == 0` ref-only sites (49, 64) MUST still execute `collect_reference_variants`
        //   even when `prune_doomed_at` returns `true`!
        // - `prune_doomed_at(1) == false` (Tumor has variant) MUST finalize both `variants` & `ref_variant`!
        // - `prune_doomed_at(p) == true` on other doomed keys MUST elide `entry.variants` (`capacity() == 0`)!
        let ub3_som_ref =
            som_builder.build_aligned_variants_with_ref_mode(&ub3_vd, RefSiteMode::KeepAll);
        let ub3_som_pruned =
            som_builder.build_aligned_variants_pruned(&ub3_vd, RefSiteMode::KeepAll, &|p| p != 1);
        assert_eq!(
            ub3_som_pruned.keys().copied().collect::<Vec<_>>(),
            ub3_som_ref.keys().copied().collect::<Vec<_>>(),
            "somatic UB-3 key order across 16-entry map must match KeepAll exactly"
        );
        for &rp in &[49, 64] {
            let got_r = ub3_som_pruned.get(&rp).unwrap();
            let exp_r = ub3_som_ref.get(&rp).unwrap();
            assert!(got_r.variants.is_empty());
            assert_eq!(got_r.variants.capacity(), 0);
            assert!(got_r.var_desc_to_variant.is_empty());
            let rv = got_r
                .ref_variant
                .as_ref()
                .expect("non_ref_rows == 0 must finalize ref_variant");
            assert_eq!(rv.start_pos, rp);
            assert_eq!(rv.genotype, exp_r.ref_variant.as_ref().unwrap().genotype);
            assert_eq!(rv.bias, exp_r.ref_variant.as_ref().unwrap().bias);
        }
        // Pos 24 (empty `var` in somatic mode): present in map with `ref_variant: None`, `variants: []`.
        let empty24 = ub3_som_pruned.get(&24).unwrap();
        assert!(empty24.ref_variant.is_none());
        assert!(empty24.variants.is_empty());
        assert_eq!(empty24.variants.capacity(), 0);
        assert!(empty24.var_desc_to_variant.is_empty());
        // Pos 1 (`prune_doomed_at(1) == false`): NOT pruned, `variants[0]` finalized & `ref_variant` retained.
        let unpruned_1 = ub3_som_pruned.get(&1).unwrap();
        let exp_1 = ub3_som_ref.get(&1).unwrap();
        assert_eq!(unpruned_1.variants.len(), 1);
        assert_eq!(unpruned_1.variants[0].start_pos, 1);
        assert!(unpruned_1.variants[0].bias.contains(';'));
        assert_eq!(
            unpruned_1.ref_variant.as_ref().unwrap().pos_coverage,
            exp_1.ref_variant.as_ref().unwrap().pos_coverage
        );
        assert_eq!(
            unpruned_1.ref_variant.as_ref().unwrap().bias,
            exp_1.ref_variant.as_ref().unwrap().bias
        );
        assert!(unpruned_1.var_desc_to_variant.is_empty());
        // Other doomed keys (`prune_doomed_at(p) == true`, plus low-freq pos 25 in somatic mode): pruned with cap 0.
        for &dk in doomed_keys.iter().filter(|&&p| p != 1).chain([&25]) {
            let de = ub3_som_pruned.get(&dk).unwrap();
            assert!(de.variants.is_empty());
            assert_eq!(de.variants.capacity(), 0);
            assert!(de.var_desc_to_variant.is_empty());
        }
    }

    #[test]
    fn drain_vars_into_reserves_exact_capacity_for_fresh_vars() {
        let mut ref_var = Variant::empty_internal();
        ref_var.description = "A".to_string();
        ref_var.frequency = 0.8;

        let mut alt_var = Variant::empty_internal();
        alt_var.description = "G".to_string();
        alt_var.frequency = 0.2;

        let mut var = vec![ref_var, alt_var];
        let mut entry = crate::variations::Vars::new();
        let maxfreq = ToVarsBuilder::drain_vars_into(&mut entry, "A", 1, &mut var);

        assert!((maxfreq - 0.2).abs() < 1e-12);
        assert!(var.is_empty());
        assert!(entry.ref_variant.is_some());
        assert_eq!(entry.variants.len(), 1);
        assert_eq!(
            entry.variants.capacity(),
            1,
            "fresh 1-allele Vars populated via drain_vars_into must allocate exact capacity 1 instead of 4"
        );

        // Non-empty `entry.variants` falls back to amortised `reserve(non_ref_rows)`.
        let mut second_alt = Variant::empty_internal();
        second_alt.description = "T".to_string();
        second_alt.frequency = 0.25;
        let mut extra = vec![second_alt];
        let maxfreq2 = ToVarsBuilder::drain_vars_into(&mut entry, "A", 1, &mut extra);
        assert!((maxfreq2 - 0.25).abs() < 1e-12);
        assert_eq!(entry.variants.len(), 2);
        assert!(entry.variants.capacity() >= 2);
    }

    #[test]
    fn compat_off_standard_dup_symbolic_alt_normalizes_negative_dup_and_preserves_legacy() {
        assert_eq!(get_dup_num("<dup-42>"), None);
        assert_eq!(get_dup_num("<dup80>"), Some(80));
        assert_eq!(get_dup_num_compat("<dup-42>", false), None);
        assert_eq!(get_dup_num_compat("<dup-42>", true), Some(-42));
        assert_eq!(get_dup_num_compat("<dup-120>", true), Some(-120));
        assert_eq!(get_dup_num_compat("<dup80>", true), Some(80));
        assert_eq!(get_dup_num_compat("<dup->", true), None);
        assert_eq!(get_dup_num_compat("<dup>", true), None);

        let mut ref_genome = Reference::new();
        let mut seq = vec![b'A'; 400];
        seq[99] = b'G'; // pos 100
        seq[119] = b'G'; // pos 120
        seq[139] = b'T'; // pos 140
        seq[159] = b'C'; // pos 160
        seq[179] = b'G'; // pos 180
        seq[199] = b'T'; // pos 200
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 400,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);

        let chr_lengths = HashMap::from([("chr1".to_string(), 400)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(400),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mk_dup_var = |cnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - (cnt / 2);
            v.qmean = 35.0 * cnt as f64;
            v.mqmean = 60.0 * cnt as f64;
            v.pmean = 25.0 * cnt as f64;
            v
        };

        let flank50 = "A".repeat(50);
        let desc_neg42 = format!("+{flank50}<dup-42>{flank50}");
        let desc_neg120 = format!("+{flank50}<dup-120>{flank50}");
        let desc_pos80 = format!("+{flank50}<dup80>{flank50}");
        let desc_neg20 = format!("+{flank50}<dup-20>{flank50}");

        let mut var_data = VariationData::new();

        // Pos 100: homozygous <dup-42> (mlen = 58 -> dup_span = 58, end_pos = 157, genotype = +58/+58)
        var_data.ref_coverage.insert(100, 10);
        let mut vm100 = VariationMap::new();
        vm100.sv = Some(Box::new(Sv {
            r#type: "DUP".to_string(),
            pairs: 4,
            splits: 6,
            clusters: 1,
        }));
        var_data.non_insertion_variants.insert(100, vm100);
        let mut ins100 = VariationMap::new();
        ins100.map.insert(desc_neg42.clone(), mk_dup_var(10));
        var_data.insertion_variants.insert(100, ins100);

        // Pos 120: heterozygous <dup-42> (ref 'G' freq 0.5 >= 0.01 -> genotype = G/+58)
        var_data.ref_coverage.insert(120, 20);
        let mut vm120 = VariationMap::new();
        vm120.sv = Some(Box::new(Sv {
            r#type: "DUP".to_string(),
            pairs: 4,
            splits: 6,
            clusters: 1,
        }));
        vm120.map.insert("G".to_string(), mk_dup_var(10));
        var_data.non_insertion_variants.insert(120, vm120);
        let mut ins120 = VariationMap::new();
        ins120.map.insert(desc_neg42.clone(), mk_dup_var(10));
        var_data.insertion_variants.insert(120, ins120);

        // Pos 140: homozygous <dup-120> (extreme clamp -> dup_span = 1, end_pos = 140, genotype = +1/+1)
        var_data.ref_coverage.insert(140, 10);
        let mut vm140 = VariationMap::new();
        vm140.sv = Some(Box::new(Sv {
            r#type: "DUP".to_string(),
            pairs: 4,
            splits: 6,
            clusters: 1,
        }));
        var_data.non_insertion_variants.insert(140, vm140);
        let mut ins140 = VariationMap::new();
        ins140.map.insert(desc_neg120.clone(), mk_dup_var(10));
        var_data.insertion_variants.insert(140, ins140);

        // Pos 160: malformed "+<dup->" and "+<dup>" must not match <dup-N>
        var_data.ref_coverage.insert(160, 20);
        let mut vm160 = VariationMap::new();
        vm160.map.insert("C".to_string(), mk_dup_var(10));
        var_data.non_insertion_variants.insert(160, vm160);
        let mut ins160 = VariationMap::new();
        ins160.map.insert("+<dup->".to_string(), mk_dup_var(5));
        ins160.map.insert("+<dup>".to_string(), mk_dup_var(5));
        var_data.insertion_variants.insert(160, ins160);

        // Pos 180: homozygous positive <dup80> (180 bp, END = 359; legacy genotype1 = +130 vs COMPAT_OFF = +180)
        var_data.ref_coverage.insert(180, 10);
        let mut vm180 = VariationMap::new();
        vm180.sv = Some(Box::new(Sv {
            r#type: "DUP".to_string(),
            pairs: 4,
            splits: 6,
            clusters: 1,
        }));
        var_data.non_insertion_variants.insert(180, vm180);
        let mut ins180 = VariationMap::new();
        ins180.map.insert(desc_pos80.clone(), mk_dup_var(10));
        var_data.insertion_variants.insert(180, ins180);

        // Pos 200: homozygous <dup-20> (80 bp, END = 279, genotype = +80/+80 even under -L 50)
        var_data.ref_coverage.insert(200, 10);
        let mut vm200 = VariationMap::new();
        vm200.sv = Some(Box::new(Sv {
            r#type: "DUP".to_string(),
            pairs: 4,
            splits: 6,
            clusters: 1,
        }));
        var_data.non_insertion_variants.insert(200, vm200);
        let mut ins200 = VariationMap::new();
        ins200.map.insert(desc_neg20.clone(), mk_dup_var(10));
        var_data.insertion_variants.insert(200, ins200);

        // 1. Legacy COMPAT_VARDICTJAVA_1_8_3 preserves upstream #342 & ToVarsBuilder.java:619 behavior verbatim
        let cfg_legacy = Config::default();
        let builder_legacy = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_legacy);
        let aligned_legacy = builder_legacy.build_aligned_variants(&var_data);

        let leg100 = &aligned_legacy.get(&100).unwrap().variants[0];
        assert!(leg100.var_allele.contains("<dup-42>"));
        assert_eq!(leg100.start_pos, 100);
        assert_eq!(leg100.end_pos, 100);
        assert_eq!(leg100.genotype.as_deref(), Some("+108/+108"));

        let leg140 = &aligned_legacy.get(&140).unwrap().variants[0];
        assert!(leg140.var_allele.contains("<dup-120>"));
        assert_eq!(leg140.end_pos, 140);

        let leg180 = &aligned_legacy.get(&180).unwrap().variants[0];
        assert_eq!(leg180.var_allele, "<DUP>");
        assert_eq!(leg180.start_pos, 180);
        assert_eq!(leg180.end_pos, 359);
        assert_eq!(leg180.genotype.as_deref(), Some("+130/+180"));

        // Legacy -L 50 keeps varallele.len() > sv_min_len BEFORE get_dup_num
        let cfg_legacy_l50 = Config {
            sv_min_len: 50,
            ..Config::default()
        };
        let aligned_legacy_l50 =
            ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_legacy_l50)
                .build_aligned_variants(&var_data);
        let leg_l50_180 = &aligned_legacy_l50.get(&180).unwrap().variants[0];
        assert_eq!(leg_l50_180.var_allele, "<DUP>");
        assert_eq!(leg_l50_180.end_pos, 288);
        assert_eq!(leg_l50_180.genotype.as_deref(), Some("+130/+107"));

        // 2a. Isolated standard_dup_symbolic_alt (without canonical_sv_vcf_anchoring) normalizes
        //     <dup-N> to <DUP> with positive i64 span (>= 1) and unifies genotype1 at unshifted bp
        let mut cfg_std_dup_only = Config::default();
        cfg_std_dup_only.compat = crate::config::CompatProfile {
            standard_dup_symbolic_alt: true,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        };
        let aligned_std_only =
            ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_std_dup_only)
                .build_aligned_variants(&var_data);

        let std100 = &aligned_std_only.get(&100).unwrap().variants[0];
        assert_eq!(std100.ref_allele, "G");
        assert_eq!(std100.var_allele, "<DUP>");
        assert_eq!(std100.start_pos, 100);
        assert_eq!(std100.end_pos, 157);
        assert_eq!(std100.genotype.as_deref(), Some("+58/+58"));

        let std180 = &aligned_std_only.get(&180).unwrap().variants[0];
        assert_eq!(std180.ref_allele, "G");
        assert_eq!(std180.var_allele, "<DUP>");
        assert_eq!(std180.start_pos, 180);
        assert_eq!(std180.end_pos, 359);
        assert_eq!(std180.genotype.as_deref(), Some("+180/+180"));

        // 2b. Full COMPAT_OFF also activates canonical_sv_vcf_anchoring, shifting unanchored
        //     (pre_anchored_dup == false) <dup±N> start_pos to bp - 1 with 1-bp padding ref_allele
        //     so END - POS == dup_span
        let mut cfg_off = Config::default();
        cfg_off.compat = crate::config::COMPAT_OFF;
        let builder_off = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off);
        let aligned_off = builder_off.build_aligned_variants(&var_data);

        let off100 = &aligned_off.get(&100).unwrap().variants[0];
        assert_eq!(off100.ref_allele, "A");
        assert_eq!(off100.var_allele, "<DUP>");
        assert_eq!(off100.start_pos, 99);
        assert_eq!(off100.end_pos, 157);
        assert_eq!(off100.genotype.as_deref(), Some("+58/+58"));

        let off120 = &aligned_off.get(&120).unwrap().variants[0];
        assert_eq!(off120.ref_allele, "A");
        assert_eq!(off120.var_allele, "<DUP>");
        assert_eq!(off120.start_pos, 119);
        assert_eq!(off120.end_pos, 177);
        assert_eq!(off120.genotype.as_deref(), Some("G/+58"));

        let off140 = &aligned_off.get(&140).unwrap().variants[0];
        assert_eq!(off140.ref_allele, "A");
        assert_eq!(off140.var_allele, "<DUP>");
        assert_eq!(off140.start_pos, 139);
        assert_eq!(off140.end_pos, 140);
        assert_eq!(off140.genotype.as_deref(), Some("+1/+1"));

        let off160 = &aligned_off.get(&160).unwrap().variants;
        assert_eq!(off160.len(), 2);
        assert!(off160.iter().all(|v| v.var_allele != "<DUP>"));

        let off180 = &aligned_off.get(&180).unwrap().variants[0];
        assert_eq!(off180.ref_allele, "A");
        assert_eq!(off180.var_allele, "<DUP>");
        assert_eq!(off180.start_pos, 179);
        assert_eq!(off180.end_pos, 359);
        assert_eq!(off180.genotype.as_deref(), Some("+180/+180"));

        let off200 = &aligned_off.get(&200).unwrap().variants[0];
        assert_eq!(off200.ref_allele, "A");
        assert_eq!(off200.var_allele, "<DUP>");
        assert_eq!(off200.start_pos, 199);
        assert_eq!(off200.end_pos, 279);
        assert_eq!(off200.genotype.as_deref(), Some("+80/+80"));

        // 3. COMPAT_OFF under -L 50 (sv_min_len = 50 <= 110) evaluates get_dup_num_compat BEFORE sv_min_len
        let cfg_off_l50 = Config {
            sv_min_len: 50,
            compat: crate::config::COMPAT_OFF,
            ..Config::default()
        };
        let aligned_off_l50 = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off_l50)
            .build_aligned_variants(&var_data);

        let l50_180 = &aligned_off_l50.get(&180).unwrap().variants[0];
        assert_eq!(l50_180.ref_allele, "A");
        assert_eq!(l50_180.var_allele, "<DUP>");
        assert_eq!(l50_180.start_pos, 179);
        assert_eq!(l50_180.end_pos, 359);
        assert_eq!(l50_180.genotype.as_deref(), Some("+180/+180"));

        let l50_200 = &aligned_off_l50.get(&200).unwrap().variants[0];
        assert_eq!(l50_200.ref_allele, "A");
        assert_eq!(l50_200.var_allele, "<DUP>");
        assert_eq!(l50_200.start_pos, 199);
        assert_eq!(l50_200.end_pos, 279);
        assert_eq!(l50_200.genotype.as_deref(), Some("+80/+80"));
    }

    #[test]
    fn compat_off_shift3_3prime_alleles_shifts_deletion_and_insertion_and_guards_crispr() {
        // 1-based positions 10..=15 = G A T A C G, flanked by 'C' everywhere else
        let mut seq = vec![b'C'; 100];
        for (offset, base) in [b'G', b'A', b'T', b'A', b'C', b'G'].into_iter().enumerate() {
            seq[9 + offset] = base;
        }
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - (cnt / 2);
            v.qmean = 35.0 * cnt as f64;
            v.mqmean = 60.0 * cnt as f64;
            v.pmean = 20.0 * cnt as f64;
            v
        };

        // Case A: Dinucleotide insertion "+AT" at pos 10 (10..=14 = G A T A C -> shift3 = 3, 3 % 2 = 1)
        let mut ins_vd = VariationData::new();
        ins_vd.ref_coverage.insert(10, 20);
        let mut vm10 = VariationMap::new();
        vm10.map.insert("G".to_string(), mk_var(10));
        ins_vd.non_insertion_variants.insert(10, vm10);
        let mut ins10 = VariationMap::new();
        ins10.map.insert("+AT".to_string(), mk_var(10));
        ins_vd.insertion_variants.insert(10, ins10);

        // Legacy -3 keeps unshifted REF='G' and unrotated ALT="GAT" at shifted start_pos=13 (#292)
        let mut cfg_leg_3 = Config::default();
        cfg_leg_3.three_prime = true;
        let leg_ins_b = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg_3);
        let leg_ins = leg_ins_b
            .build_aligned_variants(&ins_vd)
            .get(&10)
            .unwrap()
            .variants[0]
            .clone();
        assert_eq!(leg_ins.shift3, 3);
        assert_eq!(leg_ins.start_pos, 13);
        assert_eq!(leg_ins.end_pos, 13);
        assert_eq!(leg_ins.ref_allele, "G");
        assert_eq!(leg_ins.var_allele, "GAT");

        // COMPAT_OFF -3 shifts REF to 13:A and rotates "+AT" by 3 % 2 = 1 -> "TA" => ALT="ATA"
        let mut cfg_off_3 = Config::default();
        cfg_off_3.three_prime = true;
        cfg_off_3.compat = crate::config::COMPAT_OFF;
        let off_ins_b = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off_3);
        let off_ins = off_ins_b
            .build_aligned_variants(&ins_vd)
            .get(&10)
            .unwrap()
            .variants[0]
            .clone();
        assert_eq!(off_ins.shift3, 3);
        assert_eq!(off_ins.start_pos, 13);
        assert_eq!(off_ins.end_pos, 13);
        assert_eq!(off_ins.ref_allele, "A");
        assert_eq!(off_ins.var_allele, "ATA");

        // COMPAT_OFF -3 + CRISPR cut_site = 15 must NOT double-shift past pos 13
        let mut cfg_off_3_crispr = cfg_off_3.clone();
        cfg_off_3_crispr.crispr_cutting_site = Some(15);
        let off_ins_c_b = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off_3_crispr);
        let off_ins_c = off_ins_c_b
            .build_aligned_variants(&ins_vd)
            .get(&10)
            .unwrap()
            .variants[0]
            .clone();
        assert_eq!(off_ins_c.start_pos, 13);
        assert_eq!(off_ins_c.end_pos, 13);
        assert_eq!(off_ins_c.ref_allele, "A");
        assert_eq!(off_ins_c.var_allele, "ATA");
        assert_eq!(off_ins_c.crispr, 0);

        // Case B: Dinucleotide deletion "-2" at pos 11 (10..=15 = G A T A C G -> shift3 = 1)
        let mut del_vd = VariationData::new();
        del_vd.ref_coverage.insert(11, 20);
        let mut vm11 = VariationMap::new();
        vm11.map.insert("A".to_string(), mk_var(10));
        vm11.map.insert("-2".to_string(), mk_var(10));
        del_vd.non_insertion_variants.insert(11, vm11);

        // Legacy -3 shifts start_pos to 11, leaves end_pos=12, REF="GAT", ALT="G" (#292)
        let leg_del = leg_ins_b
            .build_aligned_variants(&del_vd)
            .get(&11)
            .unwrap()
            .variants[0]
            .clone();
        assert_eq!(leg_del.shift3, 1);
        assert_eq!(leg_del.start_pos, 11);
        assert_eq!(leg_del.end_pos, 12);
        assert_eq!(leg_del.ref_allele, "GAT");
        assert_eq!(leg_del.var_allele, "G");

        // COMPAT_OFF -3 shifts start_pos=11, end_pos=13, REF="ATA" (matches FASTA 11..=13), ALT="A"
        let off_del = off_ins_b
            .build_aligned_variants(&del_vd)
            .get(&11)
            .unwrap()
            .variants[0]
            .clone();
        assert_eq!(off_del.shift3, 1);
        assert_eq!(off_del.start_pos, 11);
        assert_eq!(off_del.end_pos, 13);
        assert_eq!(off_del.ref_allele, "ATA");
        assert_eq!(off_del.var_allele, "A");

        // COMPAT_OFF -3 + CRISPR cut_site = 15 must NOT double-shift past pos 11..=13
        let off_del_c = off_ins_c_b
            .build_aligned_variants(&del_vd)
            .get(&11)
            .unwrap()
            .variants[0]
            .clone();
        assert_eq!(off_del_c.start_pos, 11);
        assert_eq!(off_del_c.end_pos, 13);
        assert_eq!(off_del_c.ref_allele, "ATA");
        assert_eq!(off_del_c.var_allele, "A");
        assert_eq!(off_del_c.crispr, 0);

        // Without -3, COMPAT_OFF + CRISPR cut_site = 15 still performs normal CRISPR shift
        let mut cfg_off_no3_crispr = Config::default();
        cfg_off_no3_crispr.compat = crate::config::COMPAT_OFF;
        cfg_off_no3_crispr.crispr_cutting_site = Some(15);
        let off_no3_c_b =
            ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off_no3_crispr);
        let off_del_no3_c = off_no3_c_b
            .build_aligned_variants(&del_vd)
            .get(&11)
            .unwrap()
            .variants[0]
            .clone();
        assert_eq!(off_del_no3_c.start_pos, 11);
        assert_eq!(off_del_no3_c.end_pos, 13);
        assert_eq!(off_del_no3_c.ref_allele, "ATA");
        assert_eq!(off_del_no3_c.var_allele, "A");
        assert_eq!(off_del_no3_c.crispr, 1);
    }

    #[test]
    fn compat_off_exact_binary_half_even_rounds_5e_minus_5_ties_in_tovars_and_debug() {
        let seq = vec![b'A'; 50];
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 50,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);

        let chr_lengths = HashMap::from([("chr1".to_string(), 50)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(50),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32, hicnt: i32, extracnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = hicnt;
            v.extracnt = extracnt;
            v.fwd = cnt;
            v.rev = 0;
            v.qmean = 35.0 * cnt as f64;
            v.mqmean = 60.0 * cnt as f64;
            v.pmean = 20.0 * cnt as f64;
            v
        };

        // Pos 10: 1 / 20_000 == 5e-5 for SNV "G" and insertion "+T" (frequency, hi_freq, extra_freq, debug)
        let mut vd = VariationData::new();
        vd.ref_coverage.insert(10, 20_000);
        let mut vm10 = VariationMap::new();
        vm10.map.insert("A".to_string(), mk_var(19_999, 19_999, 0));
        vm10.map.insert("G".to_string(), mk_var(1, 1, 1));
        vd.non_insertion_variants.insert(10, vm10);
        let mut ins10 = VariationMap::new();
        ins10.map.insert("+T".to_string(), mk_var(1, 1, 1));
        vd.insertion_variants.insert(10, ins10);

        // Pos 20: ref "A" has hicnt = 1 out of hicov = 20_000 (tests update_ref_variant hi_freq)
        vd.ref_coverage.insert(20, 20_000);
        let mut vm20 = VariationMap::new();
        vm20.map.insert("A".to_string(), mk_var(1, 1, 0));
        vm20.map.insert("C".to_string(), mk_var(19_999, 19_999, 0));
        vd.non_insertion_variants.insert(20, vm20);

        let mut cfg_leg = Config::default();
        cfg_leg.pileup = true;
        cfg_leg.amplicon = Some("10:0.95".to_string());
        cfg_leg.freq = -1.0;
        cfg_leg.debug = true;
        let aligned_leg = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg)
            .build_aligned_variants(&vd);
        let leg_vars10 = aligned_leg.get(&10).unwrap();
        let leg_snv = leg_vars10
            .variants
            .iter()
            .find(|v| v.description == "G")
            .unwrap();
        let leg_ins = leg_vars10
            .variants
            .iter()
            .find(|v| v.description == "+T")
            .unwrap();
        assert_eq!(leg_snv.frequency, 0.0);
        assert_eq!(leg_snv.hi_freq, 0.0);
        assert_eq!(leg_snv.extra_freq, 0.0);
        assert_eq!(leg_ins.frequency, 0.0);
        assert!(leg_snv.debug.contains("G:1:F-1:R-0:0.0000:"));
        assert!(leg_snv.debug.contains("I+T:1:F-1:R-0:0.0000:"));
        assert_eq!(
            aligned_leg
                .get(&20)
                .unwrap()
                .ref_variant
                .as_ref()
                .unwrap()
                .hi_freq,
            0.0
        );

        let mut cfg_off = cfg_leg.clone();
        cfg_off.compat = crate::config::COMPAT_OFF;
        let aligned_off = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off)
            .build_aligned_variants(&vd);
        let off_vars10 = aligned_off.get(&10).unwrap();
        let off_snv = off_vars10
            .variants
            .iter()
            .find(|v| v.description == "G")
            .unwrap();
        let off_ins = off_vars10
            .variants
            .iter()
            .find(|v| v.description == "+T")
            .unwrap();
        assert_eq!(off_snv.frequency, 0.0001);
        assert_eq!(off_snv.hi_freq, 0.0001);
        assert_eq!(off_snv.extra_freq, 0.0001);
        assert_eq!(off_ins.frequency, 0.0001);
        assert!(off_snv.debug.contains("G:1:F-1:R-0:0.0001:"));
        assert!(off_snv.debug.contains("I+T:1:F-1:R-0:0.0001:"));
        assert_eq!(
            aligned_off
                .get(&20)
                .unwrap()
                .ref_variant
                .as_ref()
                .unwrap()
                .hi_freq,
            0.0001
        );
    }

    #[test]
    fn compat_off_reconcile_vd_dp_clamps_dp_af_adjaf_and_hicov() {
        let seq = vec![b'A'; 100];
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);

        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "TEST".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32, hicnt: i32, extracnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = hicnt;
            v.extracnt = extracnt;
            v.fwd = cnt / 2;
            v.rev = cnt - (cnt / 2);
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 20.0 * f64::from(cnt);
            v
        };

        let mut vd = VariationData::new();

        // (a) #92: cnt.cnt (15) > total_pos_coverage (10) with extracnt == 0
        vd.ref_coverage.insert(10, 10);
        let mut vm10 = VariationMap::new();
        vm10.map.insert("G".to_string(), mk_var(15, 15, 0));
        vd.non_insertion_variants.insert(10, vm10);

        // (b) #386 case 1: DP = 2, VD = 4, extracnt = 2 (4 - 2 < 2 is false in Java)
        vd.ref_coverage.insert(20, 2);
        let mut vm20 = VariationMap::new();
        vm20.map.insert("T".to_string(), mk_var(4, 4, 2));
        vd.non_insertion_variants.insert(20, vm20);

        // (c) #386 case 2: DP = 2, VD = 204, extracnt = 5 (204 - 2 < 5 is false -> AF=102, ADJAF=2.5)
        vd.ref_coverage.insert(30, 2);
        let mut vm30 = VariationMap::new();
        vm30.map.insert("C".to_string(), mk_var(204, 204, 5));
        vd.non_insertion_variants.insert(30, vm30);

        // (c2) extracnt > ttcov on both SNV and insertion (ADJAF > 1.0 in legacy)
        vd.ref_coverage.insert(35, 8);
        let mut vm35 = VariationMap::new();
        vm35.map.insert("G".to_string(), mk_var(6, 6, 12));
        vd.non_insertion_variants.insert(35, vm35);
        let mut ins35 = VariationMap::new();
        ins35.map.insert("+AA".to_string(), mk_var(5, 5, 15));
        vd.insertion_variants.insert(35, ins35);

        // (d) #400 / VarDict #113: '+' key in non_insertion_variants skipped by calc_hicov -> hicnt (10) > hicov (4)
        vd.ref_coverage.insert(40, 14);
        let mut vm40 = VariationMap::new();
        vm40.map.insert("A".to_string(), mk_var(4, 4, 0));
        vm40.map.insert("+AC#GT^2".to_string(), mk_var(10, 10, 0));
        vd.non_insertion_variants.insert(40, vm40);

        // (e) Co-occurring insertion raising total_pos_coverage_adj (final_dp = 20 > vref.total_pos_cov = 6)
        //     with raw_extracnt = 8 > 6 (Binding Ruling 1 integer-lossless lookup on both non-ins and ins).
        vd.ref_coverage.insert(50, 5);
        vd.ref_coverage.insert(51, 20);
        let mut vm50 = VariationMap::new();
        vm50.map.insert("G".to_string(), mk_var(6, 6, 8));
        vd.non_insertion_variants.insert(50, vm50);
        let mut ins50 = VariationMap::new();
        ins50.map.insert("+AT".to_string(), mk_var(7, 7, 1));
        vd.insertion_variants.insert(50, ins50);

        vd.ref_coverage.insert(60, 5);
        vd.ref_coverage.insert(61, 20);
        let mut vm60 = VariationMap::new();
        vm60.map.insert("A".to_string(), mk_var(5, 5, 0));
        vd.non_insertion_variants.insert(60, vm60);
        let mut ins60 = VariationMap::new();
        // "+CG" is sorted before "+TT": 6 - 5 < 8 sets ttcov = 6 without raising total_pos_coverage = 5;
        // then "+TT" (cnt = 7, extracnt = 0) raises total_pos_coverage to cov_next = 20!
        ins60.map.insert("+CG".to_string(), mk_var(6, 6, 8));
        ins60.map.insert("+TT".to_string(), mk_var(7, 7, 0));
        vd.insertion_variants.insert(60, ins60);

        let mut cfg_leg = Config::default();
        cfg_leg.debug = true;
        let leg_builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg);
        let aligned_leg = leg_builder.build_aligned_variants(&vd);

        let leg10 = &aligned_leg.get(&10).unwrap().variants[0];
        assert_eq!(leg10.pos_coverage, 15);
        assert_eq!(leg10.total_pos_cov, 10);
        assert_eq!(leg10.frequency, 1.5);

        let leg20 = &aligned_leg.get(&20).unwrap().variants[0];
        assert_eq!(leg20.pos_coverage, 4);
        assert_eq!(leg20.total_pos_cov, 2);
        assert_eq!(leg20.frequency, 2.0);
        assert_eq!(leg20.extra_freq, 1.0);

        let leg30 = &aligned_leg.get(&30).unwrap().variants[0];
        assert_eq!(leg30.pos_coverage, 204);
        assert_eq!(leg30.total_pos_cov, 2);
        assert_eq!(leg30.frequency, 102.0);
        assert_eq!(leg30.extra_freq, 2.5);

        let leg40 = &aligned_leg.get(&40).unwrap().variants[0];
        assert_eq!(leg40.hicnt, 10);
        assert_eq!(leg40.hicov, 4);
        assert_eq!(leg40.hi_freq, 2.5);

        let mut cfg_off = cfg_leg.clone();
        cfg_off.compat = crate::config::COMPAT_OFF;
        let off_builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off);
        let aligned_off = off_builder.build_aligned_variants(&vd);

        let off10 = &aligned_off.get(&10).unwrap().variants[0];
        assert_eq!(off10.pos_coverage, 15);
        assert_eq!(off10.total_pos_cov, 15);
        assert_eq!(off10.frequency, 1.0);
        assert_eq!(off10.extra_freq, 0.0);

        let off20 = &aligned_off.get(&20).unwrap().variants[0];
        assert_eq!(off20.pos_coverage, 4);
        assert_eq!(off20.total_pos_cov, 4);
        assert_eq!(off20.frequency, 1.0);
        assert_eq!(off20.extra_freq, 0.5);

        let off30 = &aligned_off.get(&30).unwrap().variants[0];
        assert_eq!(off30.pos_coverage, 204);
        assert_eq!(off30.total_pos_cov, 204);
        assert_eq!(off30.frequency, 1.0);
        assert_eq!(off30.extra_freq, 0.0245);

        let off35_vars = &aligned_off.get(&35).unwrap().variants;
        let off35_snv = off35_vars.iter().find(|v| v.description == "G").unwrap();
        let off35_ins = off35_vars.iter().find(|v| v.description == "+AA").unwrap();
        assert_eq!(off35_snv.total_pos_cov, 8);
        assert_eq!(off35_snv.extra_freq, 1.0);
        assert_eq!(off35_ins.total_pos_cov, 8);
        assert_eq!(off35_ins.extra_freq, 1.0);

        let off40 = &aligned_off.get(&40).unwrap().variants[0];
        assert_eq!(off40.hicnt, 10);
        assert_eq!(off40.hicov, 10);
        assert_eq!(off40.hi_freq, 1.0);
        assert!(off40.debug.contains(":1.0000:"));

        let off50_vars = &aligned_off.get(&50).unwrap().variants;
        let off50_snv = off50_vars.iter().find(|v| v.description == "G").unwrap();
        let off50_ins = off50_vars.iter().find(|v| v.description == "+AT").unwrap();
        assert_eq!(off50_snv.total_pos_cov, 20);
        assert_eq!(off50_snv.frequency, 0.3);
        assert_eq!(off50_snv.extra_freq, 0.4);
        assert_eq!(off50_ins.total_pos_cov, 20);
        assert_eq!(off50_ins.frequency, 0.35);
        assert_eq!(off50_ins.extra_freq, 0.05);

        let off60_vars = &aligned_off.get(&60).unwrap().variants;
        let off60_cg = off60_vars.iter().find(|v| v.description == "+CG").unwrap();
        assert_eq!(off60_cg.total_pos_cov, 20);
        assert_eq!(off60_cg.frequency, 0.3);
        assert_eq!(off60_cg.extra_freq, 0.4);

        for &pos in aligned_off.keys() {
            let vars = aligned_off.get(&pos).unwrap();
            let single = off_builder
                .build_vars_at_position_with_ref_mode(&vd, pos, RefSiteMode::DropAll)
                .unwrap();
            assert_eq!(single.variants.len(), vars.variants.len());
            for (v, sv) in vars.variants.iter().zip(single.variants.iter()) {
                assert!(0 <= v.pos_coverage && v.pos_coverage <= v.total_pos_cov);
                assert!(0.0 <= v.frequency && v.frequency <= 1.0);
                assert!(0.0 <= v.extra_freq && v.extra_freq <= 1.0);
                assert!(0 <= v.hicnt && v.hicnt <= v.hicov);
                assert!(0.0 <= v.hi_freq && v.hi_freq <= 1.0);
                assert_eq!(sv.total_pos_cov, v.total_pos_cov);
                assert_eq!(sv.frequency.to_bits(), v.frequency.to_bits());
                assert_eq!(sv.extra_freq.to_bits(), v.extra_freq.to_bits());
                assert_eq!(sv.hicov, v.hicov);
                assert_eq!(sv.hi_freq.to_bits(), v.hi_freq.to_bits());
            }
            if let Some(ref rv) = vars.ref_variant {
                assert!(0 <= rv.pos_coverage && rv.pos_coverage <= rv.total_pos_cov);
                assert!(0 <= rv.hicnt && rv.hicnt <= rv.hicov);
                assert!(0.0 <= rv.hi_freq && rv.hi_freq <= 1.0);
                assert!(0.0 <= rv.extra_freq && rv.extra_freq <= 1.0);
            }
        }

        // (f) Finalized reference variants (`Vars.ref_variant` via `hydrate_compact_ref_site` /
        //     `build_ref_variant_from_cnt` and `collect_reference_variants` -> `update_ref_variant`):
        //     - pos 70: pure-reference site with `RD = 15 > ref_coverage = 10`, `cnt.hicnt = -1`
        //       (from `CigarParser` `Variation::sub_cnt`), and `cnt.extracnt = 5`.
        //     - pos 75: variant + co-occurring insertion raising `total_pos_coverage_adj` from 10
        //       to `ref_coverage[76] = 25` above ref's initial `ttcov = 12`, with ref `hicnt = -1`
        //       and `extracnt = 15` (> 12, testing lossless `raw_extracnt` in `update_ref_variant`).
        let mut ref_vd = VariationData::new();
        ref_vd.ref_coverage.insert(70, 10);
        let mut vm70 = VariationMap::new();
        let mut r70 = mk_var(15, -1, 7);
        r70.fwd = 8;
        r70.rev = 7;
        r70.locnt = 2;
        vm70.map.insert("A".to_string(), r70);
        ref_vd.non_insertion_variants.insert(70, vm70);

        ref_vd.ref_coverage.insert(75, 10);
        ref_vd.ref_coverage.insert(76, 25);
        let mut vm75 = VariationMap::new();
        let mut r75 = mk_var(12, -1, 15);
        r75.fwd = 6;
        r75.rev = 6;
        r75.locnt = 2;
        vm75.map.insert("A".to_string(), r75);
        vm75.map.insert("G".to_string(), mk_var(4, 4, 0));
        ref_vd.non_insertion_variants.insert(75, vm75);
        let mut ins75 = VariationMap::new();
        ins75.map.insert("+AT".to_string(), mk_var(11, 11, 0));
        ref_vd.insertion_variants.insert(75, ins75);

        let mut vm76 = VariationMap::new();
        vm76.map.insert("A".to_string(), mk_var(25, 25, 0));
        ref_vd.non_insertion_variants.insert(76, vm76);

        // Legacy preserves negative hicnt/hicov (-1), clobbered DP=10 (< RD=15), and unrounded extra_freq (7/15)
        let mut leg_som_cfg = Config::default();
        leg_som_cfg.bam = "tumor.bam|normal.bam".to_string();
        let leg_som_b = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &leg_som_cfg);
        let leg_keep =
            leg_som_b.build_aligned_variants_with_ref_mode(&ref_vd, RefSiteMode::KeepAll);
        let leg_rv70 = leg_keep.get(&70).unwrap().ref_variant.as_ref().unwrap();
        assert_eq!(leg_rv70.total_pos_cov, 10);
        assert_eq!(leg_rv70.hicnt, -1);
        assert_eq!(leg_rv70.hicov, -1);
        assert_eq!(
            leg_rv70.extra_freq.to_bits(),
            (7.0_f64 / 15.0_f64).to_bits()
        );

        // COMPAT_OFF somatic CompactOut + hydrate_compact_ref_site vs KeepAll
        let mut off_som_cfg = leg_som_cfg.clone();
        off_som_cfg.compat = crate::config::COMPAT_OFF;
        let off_som_b = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &off_som_cfg);
        let off_keep =
            off_som_b.build_aligned_variants_with_ref_mode(&ref_vd, RefSiteMode::KeepAll);
        let mut compact_sites = rustc_hash::FxHashMap::default();
        let _ = off_som_b.build_aligned_variants_with_ref_mode(
            &ref_vd,
            RefSiteMode::CompactOut(&mut compact_sites),
        );
        let hydrated_70 = off_som_b
            .hydrate_compact_ref_site(70, compact_sites.get(&70).expect("pos 70 compact site"));
        let hyd_rv70 = hydrated_70.ref_variant.as_ref().unwrap();
        let keep_rv70 = off_keep.get(&70).unwrap().ref_variant.as_ref().unwrap();

        for rv in [hyd_rv70, keep_rv70] {
            assert_eq!(rv.total_pos_cov, 15);
            assert_eq!(rv.pos_coverage, 0);
            assert_eq!(rv.frequency, 0.0);
            assert_eq!(rv.ref_fwd_cov, 8);
            assert_eq!(rv.ref_rev_cov, 7);
            assert_eq!(rv.hicnt, 0);
            assert_eq!(rv.hicov, 0);
            assert_eq!(rv.hi_freq, 0.0);
            assert_eq!(rv.extra_freq, 0.4667);
        }

        // Un-updated ref_variant on pos 75 (when !pileup) also clamps hicnt >= 0 and hicov >= hicnt,
        // and negative ref hicnt (-1) does not deflate SNV "G"'s hicov below 4!
        let off_rv75_unupdated = off_keep.get(&75).unwrap().ref_variant.as_ref().unwrap();
        assert_eq!(off_rv75_unupdated.hicnt, 0);
        assert_eq!(off_rv75_unupdated.hicov, 4);
        assert_eq!(off_rv75_unupdated.hi_freq, 0.0);
        let off_snv75 = off_keep
            .get(&75)
            .unwrap()
            .variants
            .iter()
            .find(|v| v.description == "G")
            .unwrap();
        assert_eq!(off_snv75.hicnt, 4);
        assert_eq!(off_snv75.hicov, 4);
        assert_eq!(off_snv75.hi_freq, 1.0);

        // -p pileup + Amplicon (`collect_reference_variants:2107`) reconciles ref_variant DP to 25
        // and recomputes extra_freq losslessly from raw_extracnt = 15 -> 15 / 25 = 0.6!
        let mut off_pileup_amp_cfg = Config::default();
        off_pileup_amp_cfg.pileup = true;
        off_pileup_amp_cfg.amplicon = Some("10:0.95".to_string());
        off_pileup_amp_cfg.compat = crate::config::COMPAT_OFF;
        let off_pa_b = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &off_pileup_amp_cfg);
        let off_pa = off_pa_b.build_aligned_variants(&ref_vd);
        let pa_rv75 = off_pa.get(&75).unwrap().ref_variant.as_ref().unwrap();
        assert_eq!(pa_rv75.total_pos_cov, 25);
        assert_eq!(pa_rv75.pos_coverage, 0);
        assert_eq!(pa_rv75.frequency, 0.0);
        assert_eq!(pa_rv75.ref_fwd_cov, 7);
        assert_eq!(pa_rv75.ref_rev_cov, 7);
        assert_eq!(pa_rv75.hicnt, 0);
        assert_eq!(pa_rv75.hicov, 4);
        assert_eq!(pa_rv75.hi_freq, 0.0);
        assert_eq!(pa_rv75.extra_freq, 0.6);

        // Verify `adjust_insertion_coverage_and_next_ref` cross-strand shortfall transfer under `-u`
        // (`Switch 34` `mate_overlap_both_strands` + `unique_mode_alignment_enabled`):
        // At pos 76 (`cov_next = 25`), skew `vm76["A"]` to `fwd = 2, rev = 23`, while `ins75["+AT"]`
        // has `cnt = 11` (`fwd = 5, rev = 6` after `-u` rebalancing; vs `fwd = 0, rev = 11` before):
        // - With `-u` + `Switch 34`: `rem_f = 5 - 2 = 3` transfers onto `rev` (`23 - 6 - 3 = 14`),
        //   preserving exact `ref_fwd_cov + ref_rev_cov == 0 + 14 == 14` (`total_pos_cov = 25`, NOT `0 + 17 + 11 = 28`!).
        // - Without `-u` (`unique_mode_alignment_enabled = false`), `Switch 34` is a 100% no-op!
        let mut skew_vd = ref_vd.clone();
        if let Some(v76) = skew_vd
            .non_insertion_variants
            .get_mut(&76)
            .and_then(|m| m.map.get_mut("A"))
        {
            v76.fwd = 2;
            v76.rev = 23;
        }
        let mut u_off_cfg = off_pileup_amp_cfg.clone();
        u_off_cfg.unique_mode_alignment_enabled = true;
        let u_off_res = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &u_off_cfg)
            .build_aligned_variants(&skew_vd);
        let u_rv75 = u_off_res.get(&75).unwrap().ref_variant.as_ref().unwrap();
        assert_eq!(
            (u_rv75.ref_fwd_cov, u_rv75.ref_rev_cov, u_rv75.total_pos_cov),
            (0, 14, 25)
        );

        let mut nou_loo34_cfg = off_pileup_amp_cfg.clone();
        nou_loo34_cfg.compat.mate_overlap_both_strands = false;
        let nou_off_res =
            ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &off_pileup_amp_cfg)
                .build_aligned_variants(&skew_vd);
        let nou_loo_res = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &nou_loo34_cfg)
            .build_aligned_variants(&skew_vd);
        let nou_off_rv = nou_off_res.get(&75).unwrap().ref_variant.as_ref().unwrap();
        let nou_loo_rv = nou_loo_res.get(&75).unwrap().ref_variant.as_ref().unwrap();
        assert_eq!(
            (
                nou_off_rv.ref_fwd_cov,
                nou_off_rv.ref_rev_cov,
                nou_off_rv.total_pos_cov
            ),
            (
                nou_loo_rv.ref_fwd_cov,
                nou_loo_rv.ref_rev_cov,
                nou_loo_rv.total_pos_cov
            )
        );
    }

    #[test]
    fn compat_off_coordinate_sorted_vars_tovars_and_single_pos_agree() {
        let seq: Vec<u8> = (0..80).map(|i| b"ACGT"[i % 4]).collect();
        let ref_base_str = |pos: i32| -> String { (seq[(pos - 1) as usize] as char).to_string() };

        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 80,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq.clone());
        let chr_lengths = HashMap::from([("chr1".to_string(), 80)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(80),
            gene: "SOM".to_string(),
            amplicon: None,
        };

        let mk_cnt = |cnt: i32, fwd: i32, rev: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.fwd = fwd;
            v.rev = rev;
            v.hicnt = cnt;
            v.qmean = 32.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 24.0 * f64::from(cnt);
            v
        };

        let mut vd = VariationData::new();
        // Populate >= 13 keys (<= 24) so JavaHashMap capacity is 32, placing pos 32 (bucket 0)
        // BEFORE pos 31 (bucket 31) in JavaHashMap iteration order!
        for p in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24] {
            vd.ref_coverage.insert(p, 20);
            let mut vm = VariationMap::new();
            vm.map.insert(ref_base_str(p), mk_cnt(20, 10, 10));
            vd.non_insertion_variants.insert(p, vm);
        }

        vd.ref_coverage.insert(31, 8);
        let mut vm31 = VariationMap::new();
        vm31.map.insert(ref_base_str(31), mk_cnt(8, 4, 4));
        vd.non_insertion_variants.insert(31, vm31);
        let mut ins31 = VariationMap::new();
        ins31.map.insert("+TTAA".to_string(), mk_cnt(12, 7, 5));
        vd.insertion_variants.insert(31, ins31);

        vd.ref_coverage.insert(32, 40);
        let mut vm32 = VariationMap::new();
        vm32.map.insert(ref_base_str(32), mk_cnt(30, 17, 13));
        vm32.map.insert("-2".to_string(), mk_cnt(10, 5, 5));
        vd.non_insertion_variants.insert(32, vm32);

        assert_eq!(vd.non_insertion_variants.capacity(), 32);
        assert!(
            !vd.non_insertion_variants
                .contains_both_and_first_precedes(&31, &32),
            "pos 32 (bucket 0) must precede pos 31 (bucket 31) in cap-32 JavaHashMap"
        );

        // 1. Legacy COMPAT_VARDICTJAVA_1_8_3: pos 32 visited before pos 31 -> un-decremented ref fwd/rev (17, 13)
        let mut cfg_leg = Config::default();
        cfg_leg.bam = "tumor.bam|normal.bam".to_string();
        let leg_builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg);
        let leg_full = leg_builder.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::KeepAll);
        let leg_32 = &leg_full.get(&32).unwrap().variants[0];
        assert_eq!((leg_32.ref_fwd_cov, leg_32.ref_rev_cov), (17, 13));
        let leg_single_32 = leg_builder.build_vars_at_position(&vd, 32).unwrap();
        assert_eq!(
            (
                leg_single_32.variants[0].ref_fwd_cov,
                leg_single_32.variants[0].ref_rev_cov
            ),
            (17, 13)
        );

        // 2. COMPAT_OFF: ascending coordinate order visits pos 31 before pos 32 -> decrements overlay[32] to (10, 8),
        //    and single-position build_vars_at_position_with_ref_mode agrees 100% across all positions!
        let mut cfg_off = cfg_leg.clone();
        cfg_off.compat = crate::config::COMPAT_OFF;
        let off_builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off);
        let off_full_keep =
            off_builder.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::KeepAll);
        let off_full_drop =
            off_builder.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::DropAll);

        let off_32 = &off_full_keep.get(&32).unwrap().variants[0];
        assert_eq!((off_32.ref_fwd_cov, off_32.ref_rev_cov), (10, 8));

        for pos in 1..=40 {
            let single_keep = off_builder.build_vars_at_position(&vd, pos);
            let exp_keep = off_full_keep.get(&pos).cloned();
            assert_eq!(
                single_keep.as_ref().map(|v| v.variants.len()),
                exp_keep.as_ref().map(|v| v.variants.len()),
                "KeepAll variants len at pos={pos}"
            );
            if let (Some(sk), Some(ek)) = (&single_keep, &exp_keep) {
                for (sv, ev) in sk.variants.iter().zip(ek.variants.iter()) {
                    assert_eq!(
                        (sv.ref_fwd_cov, sv.ref_rev_cov, sv.total_pos_cov),
                        (ev.ref_fwd_cov, ev.ref_rev_cov, ev.total_pos_cov),
                        "KeepAll variant cov at pos={pos}"
                    );
                }
            }

            let single_drop =
                off_builder.build_vars_at_position_with_ref_mode(&vd, pos, RefSiteMode::DropAll);
            let exp_drop = off_full_drop.get(&pos).cloned();
            assert_eq!(
                single_drop.as_ref().map(|v| v.variants.len()),
                exp_drop.as_ref().map(|v| v.variants.len()),
                "DropAll variants len at pos={pos}"
            );
            if let (Some(sd), Some(ed)) = (&single_drop, &exp_drop) {
                for (sv, ev) in sd.variants.iter().zip(ed.variants.iter()) {
                    assert_eq!(
                        (sv.ref_fwd_cov, sv.ref_rev_cov, sv.total_pos_cov),
                        (ev.ref_fwd_cov, ev.ref_rev_cov, ev.total_pos_cov),
                        "DropAll variant cov at pos={pos}"
                    );
                }
            }
        }
    }

    #[test]
    fn compat_off_canonical_sv_vcf_anchoring_three_dup_origins_inv_complex_del_and_fallbacks() {
        let mut seq: Vec<u8> = (0..600).map(|i| b"ACGT"[i % 4]).collect();
        // Explicit deterministic bases (1-based pos -> index pos - 1):
        seq[0] = b'G'; // pos 1 (contig start)
        seq[48] = b'N'; // pos 49 ('N' upstream of pos 50)
        seq[49] = b'C'; // pos 50
        seq[58] = b'n'; // pos 59 ('n' upstream of pos 60)
        seq[59] = b'T'; // pos 60
        seq[97] = b'A'; // pos 98
        seq[98] = b'C'; // pos 99 (pre-anchored <dup80> key = bp - 1)
        seq[99] = b'G'; // pos 100
        seq[118] = b'T'; // pos 119 (upstream of unanchored <dup80> at 120)
        seq[119] = b'G'; // pos 120
        seq[138] = b'A'; // pos 139
        seq[139] = b'C'; // pos 140 (Origin C1 pure long insertion)
        seq[149] = b'G'; // pos 150 (Origin C2 +150bp#A^GC)
        seq[159] = b'T'; // pos 160 (Origin C3 +150bp&AT)
        seq[198] = b'C'; // pos 199 (upstream of -200^ACGT <INV> at 200)
        seq[199] = b'G'; // pos 200
        seq[218] = b'T'; // pos 219 (upstream of -80^...<inv40>... <INV> at 220)
        seq[219] = b'A'; // pos 220
        seq[238] = b'G'; // pos 239 (upstream of -50<inv50> <INV> at 240)
        seq[239] = b'C'; // pos 240
        seq[258] = b'A'; // pos 259 (upstream of -150#A^GC <DEL> at 260)
        seq[259] = b'T'; // pos 260
        seq[268] = b'G'; // pos 269 (upstream of -150#AC^10 <DEL> at 270)
        seq[269] = b'C'; // pos 270
        seq[278] = b'T'; // pos 279 (upstream of -150&AT <DEL> at 280)
        seq[279] = b'A'; // pos 280

        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 600,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq.clone());

        let chr_lengths = HashMap::from([("chr1".to_string(), 600)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(600),
            gene: "SV_ANCHOR".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32, pre_anchored_dup: bool| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - (cnt / 2);
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 25.0 * f64::from(cnt);
            v.pstd = true;
            v.qstd = true;
            v.pre_anchored_dup = pre_anchored_dup;
            v
        };
        let sv_box = |ty: &str| {
            Some(Box::new(Sv {
                r#type: ty.to_string(),
                pairs: 4,
                splits: 6,
                clusters: 1,
            }))
        };

        let flank50 = "A".repeat(50);
        let desc_dup80 = format!("+{flank50}<dup80>{flank50}"); // dup_span = 180
        let ins150 = "C".repeat(150);
        let desc_long_ins_pure = format!("+{ins150}");
        let desc_long_ins_hash = format!("+{ins150}#A^GC");
        let desc_long_ins_amp = format!("+{ins150}&AT");
        let flank20 = "T".repeat(20);
        let desc_inv_caret_sv = format!("-80^{flank20}<inv40>{flank20}");

        let mut vd = VariationData::new();

        // Origin A: pre_anchored_dup = true at bi = 99 (already bp - 1)
        vd.ref_coverage.insert(99, 10);
        let mut vm99 = VariationMap::new();
        vm99.sv = sv_box("DUP");
        vd.non_insertion_variants.insert(99, vm99);
        let mut ins99 = VariationMap::new();
        ins99.map.insert(desc_dup80.clone(), mk_var(10, true));
        vd.insertion_variants.insert(99, ins99);

        // Origin B: pre_anchored_dup = false at bp = 120
        vd.ref_coverage.insert(120, 20);
        let mut vm120 = VariationMap::new();
        vm120.sv = sv_box("DUP");
        vm120.map.insert("G".to_string(), mk_var(10, false));
        vd.non_insertion_variants.insert(120, vm120);
        let mut ins120 = VariationMap::new();
        ins120.map.insert(desc_dup80.clone(), mk_var(10, false));
        vd.insertion_variants.insert(120, ins120);

        // Origin C1: pure +150bp (> sv_min_len = 100) at pos 140
        vd.ref_coverage.insert(140, 10);
        vd.non_insertion_variants.insert(140, VariationMap::new());
        let mut ins140 = VariationMap::new();
        ins140
            .map
            .insert(desc_long_ins_pure.clone(), mk_var(10, false));
        vd.insertion_variants.insert(140, ins140);

        // Origin C2: +150bp#A^GC (> sv_min_len = 100) at pos 150
        vd.ref_coverage.insert(150, 10);
        vd.non_insertion_variants.insert(150, VariationMap::new());
        let mut ins150_map = VariationMap::new();
        ins150_map
            .map
            .insert(desc_long_ins_hash.clone(), mk_var(10, false));
        vd.insertion_variants.insert(150, ins150_map);

        // Origin C3: +150bp&AT (> sv_min_len = 100) at pos 160
        vd.ref_coverage.insert(160, 10);
        vd.non_insertion_variants.insert(160, VariationMap::new());
        let mut ins160 = VariationMap::new();
        ins160
            .map
            .insert(desc_long_ins_amp.clone(), mk_var(10, false));
        vd.insertion_variants.insert(160, ins160);

        // <INV> Form 1: -200^ACGT at pos 200 (>= sv_min_len = 100, BEGIN_MINUS_NUMBER_CARET)
        vd.ref_coverage.insert(200, 10);
        let mut vm200 = VariationMap::new();
        vm200.sv = sv_box("INV");
        vm200.map.insert("-200^ACGT".to_string(), mk_var(10, false));
        vd.non_insertion_variants.insert(200, vm200);

        // <INV> Form 2: -80^...<inv40>... at pos 220 (< sv_min_len = 100, INV_NUM + SOME_SV_NUMBERS)
        vd.ref_coverage.insert(220, 10);
        let mut vm220 = VariationMap::new();
        vm220.sv = sv_box("INV");
        vm220
            .map
            .insert(desc_inv_caret_sv.clone(), mk_var(10, false));
        vd.non_insertion_variants.insert(220, vm220);

        // <INV> Form 3: -50<inv50> without '^' at pos 240 (simple_del == true!)
        vd.ref_coverage.insert(240, 10);
        let mut vm240 = VariationMap::new();
        vm240.sv = sv_box("INV");
        vm240
            .map
            .insert("-50<inv50>".to_string(), mk_var(10, false));
        vd.non_insertion_variants.insert(240, vm240);

        // Complex <DEL> 1: -150#A^GC at pos 260
        vd.ref_coverage.insert(260, 10);
        let mut vm260 = VariationMap::new();
        vm260.sv = sv_box("DEL");
        vm260.map.insert("-150#A^GC".to_string(), mk_var(10, false));
        vd.non_insertion_variants.insert(260, vm260);

        // Complex <DEL> 2: -150#AC^10 at pos 270
        vd.ref_coverage.insert(270, 10);
        let mut vm270 = VariationMap::new();
        vm270.sv = sv_box("DEL");
        vm270
            .map
            .insert("-150#AC^10".to_string(), mk_var(10, false));
        vd.non_insertion_variants.insert(270, vm270);

        // Complex <DEL> 3: -150&AT at pos 280
        vd.ref_coverage.insert(280, 10);
        let mut vm280 = VariationMap::new();
        vm280.sv = sv_box("DEL");
        vm280.map.insert("-150&AT".to_string(), mk_var(10, false));
        vd.non_insertion_variants.insert(280, vm280);

        // Fallback 1: position == 1 (contig start — must NEVER produce POS = 0 or REF = "N")
        vd.ref_coverage.insert(1, 40);
        let mut vm1 = VariationMap::new();
        vm1.sv = sv_box("DEL");
        vm1.map.insert("-150".to_string(), mk_var(10, false));
        vm1.map.insert("-150#A^GC".to_string(), mk_var(10, false));
        vm1.map.insert("-200^ACGT".to_string(), mk_var(10, false));
        vd.non_insertion_variants.insert(1, vm1);
        let mut ins1 = VariationMap::new();
        ins1.map.insert(desc_dup80.clone(), mk_var(10, false));
        vd.insertion_variants.insert(1, ins1);

        // Fallback 2: 'N' at pos 49 (upstream of pos 50) and 'n' at pos 59 (upstream of pos 60)
        for p in [50, 60] {
            vd.ref_coverage.insert(p, 40);
            let mut vmp = VariationMap::new();
            vmp.sv = sv_box("DEL");
            vmp.map.insert("-150".to_string(), mk_var(10, false));
            vmp.map.insert("-150#A^GC".to_string(), mk_var(10, false));
            vmp.map.insert("-200^ACGT".to_string(), mk_var(10, false));
            vd.non_insertion_variants.insert(p, vmp);
            let mut insp = VariationMap::new();
            insp.map.insert(desc_dup80.clone(), mk_var(10, false));
            vd.insertion_variants.insert(p, insp);
        }

        // 1. Default COMPAT_VARDICTJAVA_1_8_3 preserves legacy coordinates & alleles verbatim
        let cfg_leg = Config {
            sv_min_len: 100,
            compat: crate::config::COMPAT_VARDICTJAVA_1_8_3,
            ..Config::default()
        };
        let aligned_leg = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg)
            .build_aligned_variants(&vd);
        let leg140 = &aligned_leg.get(&140).unwrap().variants[0];
        assert_eq!(
            (leg140.start_pos, leg140.end_pos, leg140.var_allele.as_str()),
            (140, 291, "<DUP>")
        );
        let leg150 = &aligned_leg.get(&150).unwrap().variants[0];
        assert_eq!(leg150.ref_allele.len(), 2);
        let leg160 = &aligned_leg.get(&160).unwrap().variants[0];
        assert_eq!(leg160.var_allele, "DUP>");
        let leg200 = &aligned_leg.get(&200).unwrap().variants[0];
        assert_eq!(
            (
                leg200.start_pos,
                leg200.end_pos,
                leg200.ref_allele.as_str(),
                leg200.var_allele.as_str()
            ),
            (200, 399, "", "<INV>")
        );
        let leg240 = &aligned_leg.get(&240).unwrap().variants[0];
        assert_eq!(leg240.var_allele, "G");
        let leg270 = &aligned_leg.get(&270).unwrap().variants[0];
        assert_eq!(leg270.ref_allele.len(), 12);
        let leg1_simple_del = aligned_leg
            .get(&1)
            .unwrap()
            .variants
            .iter()
            .find(|v| v.description == "-150")
            .unwrap();
        assert_eq!(
            (
                leg1_simple_del.start_pos,
                leg1_simple_del.ref_allele.as_str()
            ),
            (0, "N")
        );

        // 2. Isolated canonical_sv_vcf_anchoring AND full COMPAT_OFF satisfy all canonical SV invariants
        let isolated_anchor = crate::config::CompatProfile {
            canonical_sv_vcf_anchoring: true,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        };
        for compat in [isolated_anchor, crate::config::COMPAT_OFF] {
            let cfg = Config {
                sv_min_len: 100,
                compat,
                ..Config::default()
            };
            let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg);
            let aligned = builder.build_aligned_variants(&vd);

            // Origin A (pre_anchored_dup = true at 99): stays at 99, END = 279 (END - POS == 180)
            let v99 = &aligned.get(&99).unwrap().variants[0];
            assert_eq!(
                (
                    v99.start_pos,
                    v99.end_pos,
                    v99.ref_allele.as_str(),
                    v99.var_allele.as_str()
                ),
                (99, 279, "C", "<DUP>")
            );
            assert_eq!(v99.genotype.as_deref(), Some("+180/+180"));

            // Origin B (pre_anchored_dup = false at 120): shifts to 119, END = 299 (END - POS == 180)
            let v120 = &aligned.get(&120).unwrap().variants[0];
            assert_eq!(
                (
                    v120.start_pos,
                    v120.end_pos,
                    v120.ref_allele.as_str(),
                    v120.var_allele.as_str()
                ),
                (119, 299, "T", "<DUP>")
            );
            assert_eq!(v120.genotype.as_deref(), Some("G/+180"));

            // Origin C1: pure +150bp at 140 -> POS = 140, END = 290, genotype = +150/+150
            let v140 = &aligned.get(&140).unwrap().variants[0];
            assert_eq!(
                (
                    v140.start_pos,
                    v140.end_pos,
                    v140.ref_allele.as_str(),
                    v140.var_allele.as_str()
                ),
                (140, 290, "C", "<DUP>")
            );
            assert_eq!(v140.genotype.as_deref(), Some("+150/+150"));

            // Origin C2: +150bp#A^GC at 150 -> POS = 150, END = 302 (152 inserted bases), 1-bp REF = "G"
            let v150 = &aligned.get(&150).unwrap().variants[0];
            assert_eq!(
                (
                    v150.start_pos,
                    v150.end_pos,
                    v150.ref_allele.as_str(),
                    v150.var_allele.as_str()
                ),
                (150, 302, "G", "<DUP>")
            );
            assert_eq!(v150.genotype.as_deref(), Some("+152/+152"));

            // Origin C3: +150bp&AT at 160 -> POS = 160, END = 310 (150 inserted bases), 1-bp REF = "T", ALT = "<DUP>"
            let v160 = &aligned.get(&160).unwrap().variants[0];
            assert_eq!(
                (
                    v160.start_pos,
                    v160.end_pos,
                    v160.ref_allele.as_str(),
                    v160.var_allele.as_str()
                ),
                (160, 310, "T", "<DUP>")
            );
            assert_eq!(v160.genotype.as_deref(), Some("+150/+150"));

            // <INV> Form 1 (-200^ACGT at 200): POS = 199, END = 399 (END - POS == 200), REF = "C"
            let v200 = &aligned.get(&200).unwrap().variants[0];
            assert_eq!(
                (
                    v200.start_pos,
                    v200.end_pos,
                    v200.ref_allele.as_str(),
                    v200.var_allele.as_str()
                ),
                (199, 399, "C", "<INV>")
            );
            assert_eq!(v200.genotype.as_deref(), Some("-200iACGT/<INV200>"));

            // <INV> Form 2 (-80^...<inv40>... at 220): POS = 219, END = 299 (END - POS == 80), REF = "T"
            let v220 = &aligned.get(&220).unwrap().variants[0];
            assert_eq!(
                (
                    v220.start_pos,
                    v220.end_pos,
                    v220.ref_allele.as_str(),
                    v220.var_allele.as_str()
                ),
                (219, 299, "T", "<INV>")
            );

            // <INV> Form 3 (-50<inv50> without '^' at 240): POS = 239, END = 289 (END - POS == 50), REF = "G", ALT = "<INV>"
            let v240 = &aligned.get(&240).unwrap().variants[0];
            assert_eq!(
                (
                    v240.start_pos,
                    v240.end_pos,
                    v240.ref_allele.as_str(),
                    v240.var_allele.as_str()
                ),
                (239, 289, "G", "<INV>")
            );
            assert_eq!(v240.genotype.as_deref(), Some("-50<inv50>/<INV50>"));

            // Complex <DEL> 1 (-150#A^GC at 260): POS = 259, END = 410, 1-bp REF = "A", ALT = "<DEL>"
            let v260 = &aligned.get(&260).unwrap().variants[0];
            assert_eq!(
                (
                    v260.start_pos,
                    v260.end_pos,
                    v260.ref_allele.as_str(),
                    v260.var_allele.as_str()
                ),
                (259, 410, "A", "<DEL>")
            );

            // Complex <DEL> 2 (-150#AC^10 at 270): POS = 269, END = 431, 1-bp REF = "G", ALT = "<DEL>"
            let v270 = &aligned.get(&270).unwrap().variants[0];
            assert_eq!(
                (
                    v270.start_pos,
                    v270.end_pos,
                    v270.ref_allele.as_str(),
                    v270.var_allele.as_str()
                ),
                (269, 431, "G", "<DEL>")
            );

            // Complex <DEL> 3 (-150&AT at 280): POS = 279, END = 431, 1-bp REF = "T", ALT = "<DEL>"
            let v280 = &aligned.get(&280).unwrap().variants[0];
            assert_eq!(
                (
                    v280.start_pos,
                    v280.end_pos,
                    v280.ref_allele.as_str(),
                    v280.var_allele.as_str()
                ),
                (279, 431, "T", "<DEL>")
            );

            // Fallbacks at pos 1 (contig start), pos 50 ('N' at 49), and pos 60 ('n' at 59)
            for (p, exp_ref) in [(1, "G"), (50, "C"), (60, "T")] {
                let vars_p = &aligned.get(&p).unwrap().variants;
                assert_eq!(vars_p.len(), 4, "expected 4 variants at fallback pos {p}");
                for v in vars_p {
                    assert_eq!(v.start_pos, p, "unshifted start_pos at fallback pos {p}");
                    assert_eq!(v.ref_allele, exp_ref, "1-bp non-N ref_allele at pos {p}");
                }
                let del_s = vars_p.iter().find(|v| v.description == "-150").unwrap();
                assert_eq!(
                    (del_s.end_pos, del_s.var_allele.as_str()),
                    (p + 149, "<DEL>")
                );
                let del_c = vars_p
                    .iter()
                    .find(|v| v.description == "-150#A^GC")
                    .unwrap();
                assert_eq!(
                    (del_c.end_pos, del_c.var_allele.as_str()),
                    (p + 150, "<DEL>")
                );
                let inv_p = vars_p
                    .iter()
                    .find(|v| v.description == "-200^ACGT")
                    .unwrap();
                assert_eq!(
                    (inv_p.end_pos, inv_p.var_allele.as_str()),
                    (p + 199, "<INV>")
                );
                let dup_p = vars_p.iter().find(|v| v.description == desc_dup80).unwrap();
                assert_eq!(
                    (dup_p.end_pos, dup_p.var_allele.as_str()),
                    (p + 179, "<DUP>")
                );
            }

            // Single-position build_vars_at_position agrees 100% with build_aligned_variants
            for &pos in aligned.keys() {
                let single = builder
                    .build_vars_at_position_with_ref_mode(&vd, pos, RefSiteMode::DropAll)
                    .unwrap();
                let full = aligned.get(&pos).unwrap();
                assert_eq!(single.variants.len(), full.variants.len());
                for (sv, fv) in single.variants.iter().zip(full.variants.iter()) {
                    assert_eq!(
                        (
                            &sv.start_pos,
                            &sv.end_pos,
                            &sv.ref_allele,
                            &sv.var_allele,
                            &sv.genotype
                        ),
                        (
                            &fv.start_pos,
                            &fv.end_pos,
                            &fv.ref_allele,
                            &fv.var_allele,
                            &fv.genotype
                        ),
                        "single vs batch mismatch at pos {pos}"
                    );
                }
            }
        }

        // 3. Isolated canonical_sv_vcf_anchoring at default sv_min_len = 1000 also shifts Origin B <dup80>
        let cfg_iso_1000 = Config {
            sv_min_len: 1000,
            compat: isolated_anchor,
            ..Config::default()
        };
        let aligned_iso_1000 =
            ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_iso_1000)
                .build_aligned_variants(&vd);
        let iso99 = &aligned_iso_1000.get(&99).unwrap().variants[0];
        assert_eq!(
            (
                iso99.start_pos,
                iso99.end_pos,
                iso99.ref_allele.as_str(),
                iso99.var_allele.as_str()
            ),
            (99, 279, "C", "<DUP>")
        );
        let iso120 = &aligned_iso_1000.get(&120).unwrap().variants[0];
        assert_eq!(
            (
                iso120.start_pos,
                iso120.end_pos,
                iso120.ref_allele.as_str(),
                iso120.var_allele.as_str()
            ),
            (119, 299, "T", "<DUP>")
        );

        // 4. Fallback 3: unloaded position - 1 (sequence_start = 100, variant at pos = 100)
        let mut partial_ref = Reference::new();
        partial_ref
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 100,
                sequence_end: 500,
            });
        partial_ref.region_sequences.push(seq[99..500].to_vec());
        let partial_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(100),
            end: Pos1(500),
            gene: "PARTIAL".to_string(),
            amplicon: None,
        };
        let mut partial_vd = VariationData::new();
        partial_vd.ref_coverage.insert(100, 40);
        let mut pvm = VariationMap::new();
        pvm.sv = sv_box("DEL");
        pvm.map.insert("-150".to_string(), mk_var(10, false));
        pvm.map.insert("-150#A^GC".to_string(), mk_var(10, false));
        pvm.map.insert("-200^ACGT".to_string(), mk_var(10, false));
        partial_vd.non_insertion_variants.insert(100, pvm);
        let mut pins = VariationMap::new();
        pins.map.insert(desc_dup80, mk_var(10, false));
        partial_vd.insertion_variants.insert(100, pins);

        let cfg_partial = Config {
            sv_min_len: 100,
            compat: crate::config::COMPAT_OFF,
            ..Config::default()
        };
        let partial_aligned =
            ToVarsBuilder::new(&partial_ref, &partial_region, &chr_lengths, &cfg_partial)
                .build_aligned_variants(&partial_vd);
        for v in &partial_aligned.get(&100).unwrap().variants {
            assert_eq!(
                (v.start_pos, v.ref_allele.as_str()),
                (100, "G"),
                "unloaded pos 99 must keep pos 100 and 'G' for {}",
                v.description
            );
        }
    }

    #[test]
    fn test_isolated_canonical_sv_anchoring_normalizes_negative_dup_span_and_genotype1() {
        let mut seq = b"GTCG".repeat(100);
        assert_eq!(seq.len(), 400);
        seq[0] = b'G'; // pos 1
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 400,
            });
        ref_genome.region_sequences.push(seq);
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), 400);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(400),
            gene: "G1".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32, pre_anchored_dup: bool| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - (cnt / 2);
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 25.0 * f64::from(cnt);
            v.pstd = true;
            v.qstd = true;
            v.pre_anchored_dup = pre_anchored_dup;
            v
        };
        let sv_box = |ty: &str| {
            Some(Box::new(Sv {
                r#type: ty.to_string(),
                pairs: 4,
                splits: 6,
                clusters: 1,
            }))
        };

        let flank50 = "ACGT".repeat(12) + "AC";
        let desc_neg30 = format!("+{flank50}<dup-30>{flank50}");
        let desc_pos20 = format!("+{flank50}<dup20>{flank50}");

        let mut vd = VariationData::new();
        // Pos 99: homozygous `pre_anchored_dup = true` `<dup-30>` (already shifted to bp - 1 = 99; ref[98] = 'C')
        vd.ref_coverage.insert(99, 12);
        let mut vm99 = VariationMap::new();
        vm99.sv = sv_box("DUP");
        vd.non_insertion_variants.insert(99, vm99);
        let mut ins99 = VariationMap::new();
        ins99.map.insert(desc_neg30.clone(), mk_var(12, true));
        vd.insertion_variants.insert(99, ins99);

        // Pos 120: heterozygous `pre_anchored_dup = false` `<dup-30>` (shifts 120 -> 119; ref[118] = 'C', ref[119] = 'G')
        vd.ref_coverage.insert(120, 40);
        let mut vm120 = VariationMap::new();
        vm120.sv = sv_box("DUP");
        vm120.map.insert("G".to_string(), mk_var(30, false));
        vd.non_insertion_variants.insert(120, vm120);
        let mut ins120 = VariationMap::new();
        ins120.map.insert(desc_neg30, mk_var(10, false));
        vd.insertion_variants.insert(120, ins120);

        // Pos 160: homozygous `pre_anchored_dup = false` `<dup20>` (`dup_span = 2 * 50 + 20 = 120`)
        vd.ref_coverage.insert(160, 12);
        let mut vm160 = VariationMap::new();
        vm160.sv = sv_box("DUP");
        vd.non_insertion_variants.insert(160, vm160);
        let mut ins160 = VariationMap::new();
        ins160.map.insert(desc_pos20, mk_var(12, false));
        vd.insertion_variants.insert(160, ins160);

        let isolated_anchor = crate::config::CompatProfile {
            canonical_sv_vcf_anchoring: true,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        };
        assert!(!isolated_anchor.standard_dup_symbolic_alt);
        let cfg = Config {
            sv_min_len: 1000,
            compat: isolated_anchor,
            ..Config::default()
        };
        let aligned = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg)
            .build_aligned_variants(&vd);

        let v99 = &aligned.get(&99).unwrap().variants[0];
        assert_eq!(
            (
                v99.start_pos,
                v99.end_pos,
                v99.ref_allele.as_str(),
                v99.var_allele.as_str(),
                v99.genotype.as_deref()
            ),
            (99, 169, "C", "<DUP>", Some("+70/+70"))
        );

        let v120 = &aligned.get(&120).unwrap().variants[0];
        assert_eq!(
            (
                v120.start_pos,
                v120.end_pos,
                v120.ref_allele.as_str(),
                v120.var_allele.as_str(),
                v120.genotype.as_deref()
            ),
            (119, 189, "C", "<DUP>", Some("G/+70"))
        );

        let v160 = &aligned.get(&160).unwrap().variants[0];
        assert_eq!(
            (
                v160.start_pos,
                v160.end_pos,
                v160.ref_allele.as_str(),
                v160.var_allele.as_str(),
                v160.genotype.as_deref()
            ),
            (159, 279, "C", "<DUP>", Some("+120/+120"))
        );
    }

    #[test]
    fn wave5_3c_switch23_nonnegative_subtracted_means_floors_negatives_neg_zero_and_nan() {
        // bcbio/bcbio-nextgen#1704: subtracted sums (< 0.0), -0.0, and NaN must floor to +0.0
        // inside fill_common_fields (before can_prune_doomed & sort_variants) and hydrate_compact_ref_site.
        let seq = vec![b'A'; 100];
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);
        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "G1".to_string(),
            amplicon: None,
        };

        let mut vd = VariationData::new();
        // Pos 10: SNV "G" with negative subtracted sums and SNV "C" with tiny negative sums (-0.01 -> -0.0)
        // plus insertion "+TT" with negative sums. Also check sort_variants order when "G" has higher cnt (10)
        // than "T" (cnt = 2, qmean = 0.0 -> metric 0.0): in legacy, G's negative mean_quality (-5.0 * 10 = -50)
        // sorts AFTER T (0.0); under Switch 23, G's mean_quality floors to 0.0 -> metric tie 0.0 broken
        // lexicographically ("G" before "T")!
        vd.ref_coverage.insert(10, 30);
        let mut vm10 = VariationMap::new();
        let mut v_neg = Variation::new();
        v_neg.cnt = 10;
        v_neg.hicnt = 10;
        v_neg.fwd = 5;
        v_neg.rev = 5;
        v_neg.pmean = -25.0;
        v_neg.qmean = -50.0;
        v_neg.mqmean = -120.0;
        v_neg.nm = -4.0;
        vm10.map.insert("G".to_string(), v_neg);

        let mut v_neg_zero = Variation::new();
        v_neg_zero.cnt = 10;
        v_neg_zero.hicnt = 10;
        v_neg_zero.fwd = 5;
        v_neg_zero.rev = 5;
        v_neg_zero.pmean = -0.1; // -0.1 / 10 = -0.01 -> rounds to -0.0 at 1 decimal!
        v_neg_zero.qmean = -0.1;
        v_neg_zero.mqmean = -0.1;
        v_neg_zero.nm = -0.1;
        vm10.map.insert("C".to_string(), v_neg_zero);

        let mut v_tie = Variation::new();
        v_tie.cnt = 2;
        v_tie.hicnt = 2;
        v_tie.fwd = 1;
        v_tie.rev = 1;
        vm10.map.insert("T".to_string(), v_tie);
        vd.non_insertion_variants.insert(10, vm10);

        let mut ins10 = VariationMap::new();
        let mut v_ins_neg = Variation::new();
        v_ins_neg.cnt = 4;
        v_ins_neg.hicnt = 4;
        v_ins_neg.fwd = 2;
        v_ins_neg.rev = 2;
        v_ins_neg.pmean = -12.0;
        v_ins_neg.qmean = -40.0;
        v_ins_neg.mqmean = -60.0;
        v_ins_neg.nm = -2.0;
        ins10.map.insert("+TT".to_string(), v_ins_neg);
        vd.insertion_variants.insert(10, ins10);

        // Pos 20: doomed pruning check when read_pos = 0.0 and base_qual = 0.0:
        // In legacy, negative mean_position (-2.0 < 0.0) prunes pos 20 as doomed;
        // under Switch 23, mean_position floors to 0.0 (>= 0.0), so pos 20 is NOT pruned!
        vd.ref_coverage.insert(20, 20);
        let mut vm20 = VariationMap::new();
        let mut v20 = Variation::new();
        v20.cnt = 5;
        v20.hicnt = 5;
        v20.fwd = 3;
        v20.rev = 2;
        v20.pmean = -10.0;
        v20.qmean = -15.0;
        v20.mqmean = 200.0;
        vm20.map.insert("G".to_string(), v20);
        vd.non_insertion_variants.insert(20, vm20);

        let mut cfg_leg = Config::default();
        cfg_leg.read_pos = 0.0;
        cfg_leg.base_qual = 0.0;
        let b_leg = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg);
        let al_leg = b_leg.build_aligned_variants_pruned(&vd, RefSiteMode::DropAll, &|p| p == 20);
        let leg_vars10 = &al_leg.get(&10).unwrap().variants;
        let leg_g = leg_vars10.iter().find(|v| v.description == "G").unwrap();
        assert_eq!(leg_g.mean_position, -2.5);
        assert_eq!(leg_g.mean_quality, -5.0);
        assert_eq!(leg_g.mean_mapping_quality, -12.0);
        assert_eq!(leg_g.nm, -0.4);
        let leg_c = leg_vars10.iter().find(|v| v.description == "C").unwrap();
        assert!(leg_c.mean_position.is_sign_negative());
        // Legacy pruned pos 20 because -2.0 < 0.0
        assert!(al_leg.get(&20).unwrap().variants.is_empty());

        let mut cfg_sw23 = cfg_leg.clone();
        cfg_sw23.compat = crate::config::CompatProfile {
            nonnegative_subtracted_means: true,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        };
        let b_sw23 = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_sw23);
        let al_sw23 = b_sw23.build_aligned_variants_pruned(&vd, RefSiteMode::DropAll, &|p| p == 20);
        let sw23_vars10 = &al_sw23.get(&10).unwrap().variants;
        // "+TT", "C", "G", "T" all have mean_quality == 0.0, so lexicographic tie-break puts "G" before "T"!
        let g_idx = sw23_vars10
            .iter()
            .position(|v| v.description == "G")
            .unwrap();
        let t_idx = sw23_vars10
            .iter()
            .position(|v| v.description == "T")
            .unwrap();
        assert!(
            g_idx < t_idx,
            "floored mean_quality must run before sort_variants"
        );
        for desc in ["G", "C", "+TT"] {
            let v = sw23_vars10.iter().find(|v| v.description == desc).unwrap();
            assert_eq!(
                v.mean_position.to_bits(),
                0.0_f64.to_bits(),
                "{desc} mean_position"
            );
            assert_eq!(
                v.mean_quality.to_bits(),
                0.0_f64.to_bits(),
                "{desc} mean_quality"
            );
            assert_eq!(
                v.mean_mapping_quality.to_bits(),
                0.0_f64.to_bits(),
                "{desc} mean_mapping_quality"
            );
            assert_eq!(v.nm.to_bits(), 0.0_f64.to_bits(), "{desc} nm");
        }
        // Pos 20 survived doomed-pruning because floored mean_position (0.0) >= read_pos (0.0)
        assert_eq!(al_sw23.get(&20).unwrap().variants.len(), 1);

        // hydrate_compact_ref_site with negative sums AND cnt == 0 (0.0 / 0.0 = NaN -> +0.0)
        let mut neg_ref_cnt = Variation::new();
        neg_ref_cnt.cnt = 4;
        neg_ref_cnt.pmean = -8.0;
        neg_ref_cnt.qmean = -20.0;
        neg_ref_cnt.mqmean = -40.0;
        neg_ref_cnt.nm = -1.2;
        let hyd_neg = b_sw23.hydrate_compact_ref_site(
            30,
            &CompactRefSite {
                total_cov: 10,
                duprate: 0.0,
                cnt: neg_ref_cnt,
            },
        );
        let rv_neg = hyd_neg.ref_variant.unwrap();
        assert_eq!(rv_neg.mean_position.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rv_neg.mean_quality.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rv_neg.mean_mapping_quality.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rv_neg.nm.to_bits(), 0.0_f64.to_bits());

        let hyd_nan = b_sw23.hydrate_compact_ref_site(
            31,
            &CompactRefSite {
                total_cov: 10,
                duprate: 0.0,
                cnt: Variation::new(), // 0.0 / 0.0 = NaN
            },
        );
        let rv_nan = hyd_nan.ref_variant.unwrap();
        assert_eq!(rv_nan.mean_position.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rv_nan.mean_quality.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rv_nan.mean_mapping_quality.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rv_nan.nm.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn wave5_3c_switch24_reconcile_strand_and_hq_sums_handles_neg_vd_spill_zero_strands_and_bias() {
        // AstraZeneca-NGS/VarDictJava#316, clintval/vartovcf#120, clintval/vartovcf#142:
        // Isolated Switch 24 test (with reconcile_vd_dp = false, reconcile_rd_vd_dp = false).
        let seq = vec![b'A'; 100];
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);
        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "G1".to_string(),
            amplicon: None,
        };

        let mut vd = VariationData::new();
        // Pos 10: ref "A" has fwd=8, rev=8 -> ref_bias = 2.
        // - SNV "G": cnt = 2 (vd = 2), fwd = 6, rev = 6, hicnt = 5, locnt = 9 ->
        //   fwd+rev=12 > 2 subtracts 6 from fwd (0) and spills 4 onto rev (2) -> (0, 2),
        //   hicnt = 2, locnt = 0 -> hq_lq_ratio = 2 / 0.5 = 4.0, bias = "2;0" (not "2;2"!).
        // - SNV "T": cnt = 3 (vd = 3), fwd = 0, rev = 0, hicnt = -2 ->
        //   fwd+rev == 0 && vd > 0 -> (3, 0), hicnt = 0, locnt = 3 -> hq_lq_ratio = 0.0, bias = "2;0".
        // - SNV "C": cnt = -3 (pos_coverage = -3!), fwd = 4, rev = 2, hicnt = 5 ->
        //   must NOT panic on .clamp(0, -3)! Reconciles vd = 0, fwd = 0, rev = 0, hicnt = 0, locnt = 0, bias = "2;0".
        // - Insertion "+AT": cnt = 4, fwd = 1, rev = 1, hicnt = 7 ->
        //   sum = 2 < 4 adds +2 to larger/fwd -> (3, 1), hicnt = 4, locnt = 0 -> bias = "2;2".
        vd.ref_coverage.insert(10, 20);
        let mut vm10 = VariationMap::new();
        let mut v_ref = Variation::new();
        v_ref.cnt = 16;
        v_ref.hicnt = 16;
        v_ref.fwd = 8;
        v_ref.rev = 8;
        v_ref.qmean = 560.0;
        v_ref.mqmean = 960.0;
        v_ref.pmean = 320.0;
        vm10.map.insert("A".to_string(), v_ref);

        let mut v_spill = Variation::new();
        v_spill.cnt = 2;
        v_spill.fwd = 6;
        v_spill.rev = 6;
        v_spill.hicnt = 5;
        v_spill.locnt = 9;
        v_spill.qmean = 70.0;
        v_spill.mqmean = 120.0;
        v_spill.pmean = 40.0;
        vm10.map.insert("G".to_string(), v_spill);

        let mut v_zero_strands = Variation::new();
        v_zero_strands.cnt = 3;
        v_zero_strands.fwd = 0;
        v_zero_strands.rev = 0;
        v_zero_strands.hicnt = -2;
        v_zero_strands.locnt = 10;
        v_zero_strands.qmean = 105.0;
        v_zero_strands.mqmean = 180.0;
        v_zero_strands.pmean = 60.0;
        vm10.map.insert("T".to_string(), v_zero_strands);

        let mut v_neg_cnt = Variation::new();
        v_neg_cnt.cnt = -3;
        v_neg_cnt.fwd = 4;
        v_neg_cnt.rev = 2;
        v_neg_cnt.hicnt = 5;
        v_neg_cnt.locnt = 1;
        v_neg_cnt.qmean = -90.0;
        v_neg_cnt.mqmean = -180.0;
        v_neg_cnt.pmean = -60.0;
        vm10.map.insert("C".to_string(), v_neg_cnt);
        vd.non_insertion_variants.insert(10, vm10);

        let mut ins10 = VariationMap::new();
        let mut v_ins = Variation::new();
        v_ins.cnt = 4;
        v_ins.fwd = 1;
        v_ins.rev = 1;
        v_ins.hicnt = 7;
        v_ins.locnt = 3;
        v_ins.qmean = 140.0;
        v_ins.mqmean = 240.0;
        v_ins.pmean = 80.0;
        ins10.map.insert("+AT".to_string(), v_ins);
        vd.insertion_variants.insert(10, ins10);

        let mut cfg_sw24 = Config::default();
        cfg_sw24.pileup = true; // keep all variants regardless of freq
        cfg_sw24.compat = crate::config::CompatProfile {
            reconcile_strand_and_hq_sums: true,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        };
        assert!(!cfg_sw24.compat.reconcile_vd_dp);
        let builder = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_sw24);
        let aligned = builder.build_aligned_variants(&vd);
        let vars10 = &aligned.get(&10).unwrap().variants;

        let g = vars10.iter().find(|v| v.description == "G").unwrap();
        assert_eq!((g.pos_coverage, g.fwd, g.rev, g.hicnt), (2, 0, 2, 2));
        assert_eq!(g.hq_lq_ratio, 4.0);
        // calc_hicov clamps v_ref (16) + v_spill (min(5,2)=2) + v_zero (0) + v_neg (0) = 18
        assert_eq!(g.hicov, 18);
        assert_eq!(g.bias, "2;0");

        let t = vars10.iter().find(|v| v.description == "T").unwrap();
        assert_eq!((t.pos_coverage, t.fwd, t.rev, t.hicnt), (3, 3, 0, 0));
        assert_eq!(t.hq_lq_ratio, 0.0);
        assert_eq!(t.bias, "2;0");

        let c = vars10.iter().find(|v| v.description == "C").unwrap();
        assert_eq!((c.pos_coverage, c.fwd, c.rev, c.hicnt), (0, 0, 0, 0));
        assert_eq!(c.frequency, 0.0);
        assert_eq!(c.bias, "2;0");

        let ins = vars10.iter().find(|v| v.description == "+AT").unwrap();
        assert_eq!(
            (ins.pos_coverage, ins.fwd, ins.rev, ins.hicnt),
            (4, 3, 1, 4)
        );
        assert_eq!(ins.hq_lq_ratio, 8.0);
        assert_eq!(ins.bias, "2;2");

        // Also test direct adjust_variant_counts on a Variant with pos_coverage = -3 and "2;1" bias prefix!
        let mut raw_v = Variant::new();
        raw_v.pos_coverage = -3;
        raw_v.fwd = 5;
        raw_v.rev = 4;
        raw_v.hicnt = 6;
        raw_v.bias = "2;1".to_string();
        builder.adjust_variant_counts(10, &mut raw_v);
        assert_eq!(
            (raw_v.pos_coverage, raw_v.fwd, raw_v.rev, raw_v.hicnt),
            (0, 0, 0, 0)
        );
        assert_eq!(raw_v.bias, "2;0");

        let mut raw_spill = Variant::new();
        raw_spill.pos_coverage = 2;
        raw_spill.fwd = 6;
        raw_spill.rev = 6;
        raw_spill.hicnt = 5;
        raw_spill.bias = "1;2".to_string();
        builder.adjust_variant_counts(10, &mut raw_spill);
        assert_eq!(
            (
                raw_spill.pos_coverage,
                raw_spill.fwd,
                raw_spill.rev,
                raw_spill.hicnt
            ),
            (2, 0, 2, 2)
        );
        assert_eq!(raw_spill.bias, "1;0");
    }

    #[test]
    fn wave5_3c_switch25_reconcile_rd_vd_dp_isolated_enforces_dp_ge_rd_plus_vd() {
        // AstraZeneca-NGS/VarDict#168, AstraZeneca-NGS/VarDictJava#323:
        // Isolated reconcile_rd_vd_dp = true with reconcile_vd_dp = false!
        let seq = vec![b'A'; 100];
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);
        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "G1".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32, fwd: i32, rev: i32, extracnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = fwd;
            v.rev = rev;
            v.extracnt = extracnt;
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 20.0 * f64::from(cnt);
            v
        };

        let mut vd = VariationData::new();
        // Pos 10: ref_coverage = 12, but RD = 8 + 7 = 15, SNV "G" VD = 6 (extracnt = 9), ins "+T" VD = 5!
        // Legacy keeps DP = 12 (< RD + VD = 21); isolated Switch 25 raises SNV DP to 15 + 6 = 21
        // (AF = 6/21 = 0.2857, extra_freq = 9/21 = 0.4286) and insertion DP to 15 + 5 = 20!
        vd.ref_coverage.insert(10, 12);
        let mut vm10 = VariationMap::new();
        vm10.map.insert("A".to_string(), mk_var(15, 8, 7, 0));
        vm10.map.insert("G".to_string(), mk_var(6, 3, 3, 9));
        vd.non_insertion_variants.insert(10, vm10);
        let mut ins10 = VariationMap::new();
        ins10.map.insert("+T".to_string(), mk_var(5, 3, 2, 0));
        vd.insertion_variants.insert(10, ins10);

        // Pos 20: pure reference site with ref_coverage = 10, RD = 9 + 6 = 15, extracnt = 3.
        vd.ref_coverage.insert(20, 10);
        let mut vm20 = VariationMap::new();
        vm20.map.insert("A".to_string(), mk_var(15, 9, 6, 3));
        vd.non_insertion_variants.insert(20, vm20);

        let mut cfg_leg = Config::default();
        cfg_leg.bam = "tumor.bam|normal.bam".to_string();
        let al_leg = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg)
            .build_aligned_variants_with_ref_mode(&vd, RefSiteMode::KeepAll);
        assert_eq!(al_leg.get(&10).unwrap().variants[0].total_pos_cov, 12);
        assert_eq!(
            al_leg
                .get(&20)
                .unwrap()
                .ref_variant
                .as_ref()
                .unwrap()
                .total_pos_cov,
            10
        );

        let mut cfg_sw25 = cfg_leg.clone();
        cfg_sw25.compat = crate::config::CompatProfile {
            reconcile_rd_vd_dp: true,
            reconcile_vd_dp: false,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        };
        assert!(!cfg_sw25.compat.reconcile_vd_dp);
        let b_sw25 = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_sw25);
        let al_sw25 = b_sw25.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::KeepAll);
        let vars10 = &al_sw25.get(&10).unwrap().variants;
        let g = vars10.iter().find(|v| v.description == "G").unwrap();
        assert_eq!(g.total_pos_cov, 21);
        assert_eq!(g.frequency, 0.2857);
        assert_eq!(g.extra_freq, 0.4286);
        let ins = vars10.iter().find(|v| v.description == "+T").unwrap();
        assert_eq!(ins.total_pos_cov, 20);
        assert_eq!(ins.frequency, 0.25);

        let rv20 = al_sw25.get(&20).unwrap().ref_variant.as_ref().unwrap();
        assert_eq!(rv20.total_pos_cov, 15);
        assert_eq!(rv20.extra_freq, 0.2);

        let mut compact = rustc_hash::FxHashMap::default();
        let _ =
            b_sw25.build_aligned_variants_with_ref_mode(&vd, RefSiteMode::CompactOut(&mut compact));
        let hyd20 = b_sw25.hydrate_compact_ref_site(20, compact.get(&20).unwrap());
        assert_eq!(hyd20.ref_variant.unwrap().total_pos_cov, 15);
    }

    #[test]
    fn wave5_3c_switch26_sanitize_iupac_alt_alleles_preserves_closed_tags_and_reduces_unclosed() {
        // AstraZeneca-NGS/VarDictJava#267, AstraZeneca-NGS/VarDictJava#287,
        // AstraZeneca-NGS/VarDict#137, bcbio/bcbio-nextgen#3321:
        // Verify symbolic <INV>/<DUP>/<DEL>, mixed "RC<dup-5>GT" -> "AC<dup-5>GT",
        // malformed unclosed "<dup-5RT" -> "<aup-5AT", and ToVarsBuilder end-to-end gating.
        for sym in ["<INV>", "<DUP>", "<DEL>"] {
            let mut s = sym.to_string();
            validate_alt_allele_in_place(&mut s);
            assert_eq!(s, sym, "symbolic SV tag {sym} must be preserved verbatim");
        }

        let mut mixed = "RC<dup-5>GT".to_string();
        validate_alt_allele_in_place(&mut mixed);
        assert_eq!(mixed, "AC<dup-5>GT");

        let mut unclosed = "<dup-5RT".to_string();
        validate_alt_allele_in_place(&mut unclosed);
        assert_eq!(unclosed, "<Aup-5AT");

        let mut all_iupac = "MRWVHDSYBKmrwvhdsybkACGTN".to_string();
        validate_alt_allele_in_place(&mut all_iupac);
        assert_eq!(all_iupac, "AAAAAACCCGAAAAAACCCGACGTN");

        // End-to-end ToVarsBuilder check: legacy preserves IUPAC ALT "R" and "+RGT",
        // while Switch 26 reduces them to "A" and "AAAGT" while keeping "<INV>" untouched!
        let seq = vec![b'A'; 300];
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 300,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);
        let chr_lengths = HashMap::from([("chr1".to_string(), 300)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(300),
            gene: "G1".to_string(),
            amplicon: None,
        };

        let mk = |cnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - cnt / 2;
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 20.0 * f64::from(cnt);
            v
        };

        let mut vd = VariationData::new();
        vd.ref_coverage.insert(50, 20);
        let mut vm50 = VariationMap::new();
        vm50.map.insert("Y".to_string(), mk(8));
        vm50.map.insert("-180^AAAA<inv80>TTTT".to_string(), mk(6));
        vd.non_insertion_variants.insert(50, vm50);
        let mut ins50 = VariationMap::new();
        ins50.map.insert("+RKG".to_string(), mk(5));
        vd.insertion_variants.insert(50, ins50);

        let cfg_leg = Config::default();
        let al_leg = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg)
            .build_aligned_variants(&vd);
        let v_leg = &al_leg.get(&50).unwrap().variants;
        assert_eq!(
            v_leg
                .iter()
                .find(|v| v.description == "Y")
                .unwrap()
                .var_allele,
            "Y"
        );
        assert_eq!(
            v_leg
                .iter()
                .find(|v| v.description == "+RKG")
                .unwrap()
                .var_allele,
            "ARKG"
        );

        let cfg_sw26 = Config {
            compat: crate::config::CompatProfile {
                sanitize_iupac_alt_alleles: true,
                ..crate::config::COMPAT_VARDICTJAVA_1_8_3
            },
            ..Config::default()
        };
        let al_sw26 = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_sw26)
            .build_aligned_variants(&vd);
        let v_sw26 = &al_sw26.get(&50).unwrap().variants;
        assert_eq!(
            v_sw26
                .iter()
                .find(|v| v.description == "Y")
                .unwrap()
                .var_allele,
            "C"
        );
        assert_eq!(
            v_sw26
                .iter()
                .find(|v| v.description == "+RKG")
                .unwrap()
                .var_allele,
            "AAGG"
        );
        assert_eq!(
            v_sw26
                .iter()
                .find(|v| v.description.starts_with("-180^"))
                .unwrap()
                .var_allele,
            "<INV>"
        );
    }

    #[test]
    fn wave5_3c_switch28_pileup_preserve_shifted_ref_hydrates_100pct_indel_and_complex_sites() {
        // AstraZeneca-NGS/VarDictJava#330: under -p (pileup = true) and pileup_preserve_shifted_ref,
        // both 100%-indel positions (ref_variant == None, start_position != position) and
        // adj_complex-candidate positions (start_position == position, changed_ref_variant_at_pos == false)
        // must carry a FASTA-hydrated ref_variant!
        let mut seq = vec![b'A'; 100];
        seq[49] = b'G'; // pos 50 = 'G'
        seq[59] = b'C'; // pos 60 = 'C'
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 100,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);
        let chr_lengths = HashMap::from([("chr1".to_string(), 100)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(100),
            gene: "G1".to_string(),
            amplicon: None,
        };

        let mk = |cnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - cnt / 2;
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 20.0 * f64::from(cnt);
            v
        };

        let mut vd = VariationData::new();
        // Pos 50: 100% "-2" deletion (no "G" ref entry -> ref_variant initially None; shifts to 49)
        vd.ref_coverage.insert(50, 20);
        let mut vm50 = VariationMap::new();
        vm50.map.insert("-2".to_string(), mk(20));
        vd.non_insertion_variants.insert(50, vm50);

        // Pos 60: 100% MNV "CA" (later trimmed by adj_complex from 60 -> 61; start_position == 60 in ToVarsBuilder!)
        vd.ref_coverage.insert(60, 18);
        let mut vm60 = VariationMap::new();
        vm60.map.insert("CA".to_string(), mk(18));
        vd.non_insertion_variants.insert(60, vm60);

        let mut cfg_leg = Config::default();
        cfg_leg.pileup = true;
        let al_leg = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg)
            .build_aligned_variants(&vd);
        assert!(
            al_leg.get(&50).unwrap().ref_variant.is_none(),
            "legacy -p drops ref_variant on 100%-indel site"
        );
        assert!(
            al_leg.get(&60).unwrap().ref_variant.is_none(),
            "legacy -p drops ref_variant on 100%-MNV site"
        );

        let mut cfg_sw28 = cfg_leg.clone();
        cfg_sw28.compat = crate::config::CompatProfile {
            pileup_preserve_shifted_ref: true,
            ..crate::config::COMPAT_VARDICTJAVA_1_8_3
        };
        let b_sw28 = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_sw28);
        let al_sw28 = b_sw28.build_aligned_variants(&vd);

        let rv50 = al_sw28
            .get(&50)
            .unwrap()
            .ref_variant
            .as_ref()
            .expect("Switch 28 must hydrate ref_variant at 100%-indel pos 50");
        assert_eq!(
            (
                rv50.start_pos,
                rv50.end_pos,
                rv50.ref_allele.as_str(),
                rv50.var_allele.as_str(),
                rv50.genotype.as_deref(),
                rv50.total_pos_cov,
                rv50.bias.as_str()
            ),
            (50, 50, "G", "G", Some("G/G"), 20, "0;0")
        );

        let rv60 = al_sw28
            .get(&60)
            .unwrap()
            .ref_variant
            .as_ref()
            .expect("Switch 28 must hydrate ref_variant at adj_complex-candidate pos 60");
        assert_eq!(
            (
                rv60.start_pos,
                rv60.end_pos,
                rv60.ref_allele.as_str(),
                rv60.var_allele.as_str(),
                rv60.genotype.as_deref(),
                rv60.total_pos_cov,
                rv60.bias.as_str()
            ),
            (60, 60, "C", "C", Some("C/C"), 18, "0;0")
        );

        // Single-position builder must match full builder
        let single50 = b_sw28.build_vars_at_position(&vd, 50).unwrap();
        assert_eq!(single50.ref_variant.unwrap().ref_allele, "G");
    }

    #[test]
    fn test_wave6_switch38_indel_only_shift3_clamps_snv_and_mnv_preserving_msi() {
        use crate::types::Pos1;
        use crate::variations::{Variation, VariationData, VariationMap};

        // Build a homopolymer/dinucleotide repeat locus so SNV & MNV trigger legacy shift3 > 0:
        // 1..=40: "A" * 40 (so SNV at pos 20 has tseq1 = ref[1..=21] = "A"*21, tseq2 = ref[22..=90] = "A"*29 -> shift3 > 0!)
        // 41..=80: "AT" * 20 (so 2-bp MNV "GC" at pos 50 has tseq1 ending in "AT", tseq2 starting with "ATAT..." -> shift3 > 0!)
        let mut seq = Vec::new();
        seq.extend_from_slice(&[b'A'; 40]); // 1..=40
        for _ in 0..20 {
            seq.extend_from_slice(b"AT"); // 41..=80
        }
        seq.extend_from_slice(b"CAAAGTTCGA"); // 81..=90: 81=C, 82=A, 83=A, 84=A, 85=G
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: seq.len() as u32,
            });
        ref_genome.region_sequences.push(seq.clone());
        let mut chr_lengths = HashMap::new();
        chr_lengths.insert("chr1".to_string(), seq.len());
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(seq.len() as u32),
            gene: "G".to_string(),
            amplicon: None,
        };

        let mk = |cnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - cnt / 2;
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 20.0 * f64::from(cnt);
            v
        };

        let mut vd = VariationData::new();
        // 1) SNV A->G at pos 20 in poly-A run:
        vd.ref_coverage.insert(20, 30);
        let mut vm20 = VariationMap::new();
        vm20.map.insert("A".to_string(), mk(20));
        vm20.map.insert("G".to_string(), mk(10));
        vd.non_insertion_variants.insert(20, vm20);

        // 2) 2-bp MNV AT->GC at pos 51 (description "G&C") in AT-repeat:
        vd.ref_coverage.insert(51, 30);
        let mut vm51 = VariationMap::new();
        vm51.map.insert("A".to_string(), mk(20));
        vm51.map.insert("G&C".to_string(), mk(10));
        vd.non_insertion_variants.insert(51, vm51);

        // 3) Real 1-bp deletion "-1" in poly-A ("AAAG" -> pos 82 deletes 'A', shift3 = 2):
        vd.ref_coverage.insert(82, 30);
        let mut vm82 = VariationMap::new();
        vm82.map.insert("A".to_string(), mk(20));
        vm82.map.insert("-1".to_string(), mk(10));
        vd.non_insertion_variants.insert(82, vm82);

        let cfg_leg = Config::default();
        let al_leg = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_leg)
            .build_aligned_variants(&vd);
        let snv_leg = &al_leg.get(&20).unwrap().variants[0];
        let mnv_leg = &al_leg.get(&51).unwrap().variants[0];
        let del_leg = &al_leg.get(&82).unwrap().variants[0];
        assert!(snv_leg.shift3 > 0, "legacy SNV in repeat has shift3 > 0");
        assert!(mnv_leg.shift3 > 0, "legacy MNV in repeat has shift3 > 0");
        assert_eq!(del_leg.shift3, 2);

        let mut cfg_sw38 = Config::default();
        cfg_sw38.compat.indel_only_shift3 = true;
        let al_sw38 = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_sw38)
            .build_aligned_variants(&vd);
        let snv_on = &al_sw38.get(&20).unwrap().variants[0];
        let mnv_on = &al_sw38.get(&51).unwrap().variants[0];
        let del_on = &al_sw38.get(&82).unwrap().variants[0];

        assert_eq!(snv_on.shift3, 0, "Switch 38 must clamp SNV shift3 to 0");
        assert_eq!(snv_on.msi.to_bits(), snv_leg.msi.to_bits());
        assert_eq!(snv_on.msint, snv_leg.msint);

        assert_eq!(mnv_on.shift3, 0, "Switch 38 must clamp MNV shift3 to 0");
        assert_eq!(mnv_on.msi.to_bits(), mnv_leg.msi.to_bits());
        assert_eq!(mnv_on.msint, mnv_leg.msint);

        assert_eq!(del_on.shift3, 2, "Switch 38 must preserve indel shift3");
        assert_eq!(del_on.msi.to_bits(), del_leg.msi.to_bits());
        assert_eq!(del_on.msint, del_leg.msint);
    }

    #[test]
    fn pr1_lane2_f12_switch26_uppercase_iupac_and_wildtype_ref_routing() {
        // 1. Lowercase IUPAC ambiguity bases reduce to UPPERCASE A/C/G under Switch 26
        let mut all_iupac = "MRWVHDSYBKmrwvhdsybkACGTN".to_string();
        validate_alt_allele_in_place(&mut all_iupac);
        assert_eq!(all_iupac, "AAAAAACCCGAAAAAACCCGACGTN");

        // 2. FASTA 'Y' at pos 27 with 50 wild-type 'C' reads (24 fwd, 26 rev) + 40 ALT 'T' reads (20 fwd, 20 rev),
        //    and FASTA 'R' at pos 60 with 15 'A' reads + 5 'M' reads (both reducing to 'A') in -p pileup mode!
        let mut seq = vec![b'A'; 120];
        seq[26] = b'Y'; // 1-based pos 27 ('Y' -> 'C')
        seq[59] = b'R'; // 1-based pos 60 ('R' -> 'A')
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr6".to_string(),
                sequence_start: 1,
                sequence_end: 120,
            });
        ref_genome
            .reference_sequences
            .insert("chr6".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);
        let chr_lengths = HashMap::from([("chr6".to_string(), 120)]);
        let region = Region {
            chr: "chr6".to_string(),
            start: Pos1(1),
            end: Pos1(120),
            gene: "IUPAC_TEST".to_string(),
            amplicon: None,
        };

        let mk_var = |cnt: i32, fwd: i32, rev: i32, q: f64, p: f64| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = fwd;
            v.rev = rev;
            v.qmean = q * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = p * f64::from(cnt);
            v.nm = 1.0 * f64::from(cnt);
            v.pstd = true;
            v.qstd = true;
            v
        };

        let mut vd = VariationData::new();
        vd.ref_coverage.insert(27, 90);
        let mut vm27 = VariationMap::new();
        vm27.map
            .insert("C".to_string(), mk_var(50, 24, 26, 37.0, 65.0));
        vm27.map
            .insert("T".to_string(), mk_var(40, 20, 20, 36.0, 55.0));
        vd.non_insertion_variants.insert(27, vm27);

        vd.ref_coverage.insert(60, 20);
        let mut vm60 = VariationMap::new();
        vm60.map
            .insert("A".to_string(), mk_var(15, 8, 7, 38.0, 70.0));
        vm60.map
            .insert("M".to_string(), mk_var(5, 2, 3, 34.0, 50.0));
        vd.non_insertion_variants.insert(60, vm60);

        let cfg_off = Config {
            compat: crate::config::COMPAT_OFF,
            ..Config::default()
        };
        let builder_off = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off);
        let aligned_off = builder_off.build_aligned_variants(&vd);

        // Pos 60 ('A' + 'M' on 'R' FASTA) is pure wild-type under Switch 26 -> dropped in DropAll!
        assert!(
            !aligned_off.contains_key(&60),
            "pure wild-type IUPAC site must be dropped in non-pileup DropAll mode"
        );

        // Pos 27 ('C' + 'T' on 'Y' FASTA): 'C' routed to ref_variant, ONLY 'T' in variants!
        let vars27 = aligned_off.get(&27).expect("pos 27 must emit 'T' SNV");
        assert_eq!(
            vars27.variants.len(),
            1,
            "wild-type 'C' must be omitted from variants"
        );
        let alt_t = &vars27.variants[0];
        assert_eq!(alt_t.ref_allele, "C");
        assert_eq!(alt_t.var_allele, "T");
        assert_eq!((alt_t.ref_fwd_cov, alt_t.ref_rev_cov), (24, 26));
        assert_eq!((alt_t.fwd, alt_t.rev), (20, 20));
        assert_eq!(alt_t.hicov, 90);
        assert_eq!(alt_t.genotype.as_deref(), Some("C/T"));
        assert_eq!(alt_t.bias, "2;2");

        // Pileup (-p) at pos 60: single merged reference row with genotype "A/A" (not "R/R") and weighted means!
        let mut cfg_pile = cfg_off.clone();
        cfg_pile.pileup = true;
        let aligned_pile = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_pile)
            .build_aligned_variants(&vd);
        let vars60_pile = aligned_pile.get(&60).expect("pos 60 must emit in -p mode");
        assert!(
            vars60_pile.variants.is_empty(),
            "-p must not double-emit wild-type alleles in entry.variants"
        );
        let rv60 = vars60_pile
            .ref_variant
            .as_ref()
            .expect("pos 60 ref_variant");
        assert_eq!(rv60.ref_allele, "A");
        assert_eq!(rv60.var_allele, "A");
        assert_eq!(rv60.genotype.as_deref(), Some("A/A"));
        assert_eq!((rv60.ref_fwd_cov, rv60.ref_rev_cov), (10, 10));
        assert_eq!(rv60.mean_quality, 37.0); // (15*38 + 5*34) / 20 = 740 / 20 = 37.0
        assert_eq!(rv60.mean_position, 65.0); // (15*70 + 5*50) / 20 = 1300 / 20 = 65.0
    }

    #[test]
    fn pr1_lane2_f14b_tovars_complex_ins_len_cutoff_and_neg_dup() {
        let seq = vec![b'A'; 300];
        let mut ref_genome = Reference::new();
        ref_genome
            .loaded_regions
            .push(crate::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: 1,
                sequence_end: 300,
            });
        ref_genome
            .reference_sequences
            .insert("chr1".to_string(), seq.clone());
        ref_genome.region_sequences.push(seq);
        let chr_lengths = HashMap::from([("chr1".to_string(), 300)]);
        let region = Region {
            chr: "chr1".to_string(),
            start: Pos1(1),
            end: Pos1(300),
            gene: "SV_CUTOFF".to_string(),
            amplicon: None,
        };

        let mk = |cnt: i32| -> Variation {
            let mut v = Variation::new();
            v.cnt = cnt;
            v.hicnt = cnt;
            v.fwd = cnt / 2;
            v.rev = cnt - cnt / 2;
            v.qmean = 35.0 * f64::from(cnt);
            v.mqmean = 60.0 * f64::from(cnt);
            v.pmean = 20.0 * f64::from(cnt);
            v
        };

        // 48-bp inserted sequence with 2-bp '&CC' matched flank: ins_len == 48 < 50 (`varallele.len() == 52 > 50`)
        let ins48_amp = format!("+{}&CC", "C".repeat(48));
        // 80-bp `<dup-20>` symbolic token with only `sv_respect_min_len: true`
        let dup80 = format!("+{}<dup-20>{}", "C".repeat(40), "G".repeat(40));

        let mut vd = VariationData::new();
        vd.ref_coverage.insert(50, 30);
        vd.non_insertion_variants.insert(50, VariationMap::new());
        let mut ins50 = VariationMap::new();
        ins50.map.insert(ins48_amp, mk(10));
        vd.insertion_variants.insert(50, ins50);

        vd.ref_coverage.insert(120, 30);
        vd.non_insertion_variants.insert(120, VariationMap::new());
        let mut ins120 = VariationMap::new();
        ins120.map.insert(dup80, mk(10));
        vd.insertion_variants.insert(120, ins120);

        // 80-bp `<inv-20>` symbolic inversion with ONLY `sv_respect_min_len: true` (`canonical_sv_vcf_anchoring: false`):
        let inv80 = format!("-80^{}<inv-20>{}", "T".repeat(50), "T".repeat(50));
        vd.ref_coverage.insert(200, 30);
        let mut non_ins200 = VariationMap::new();
        non_ins200.map.insert(inv80, mk(10));
        let sv200 = crate::variations::VariationMap::get_sv(&mut vd.non_insertion_variants, 200);
        sv200.r#type = "INV".to_string();
        sv200.splits = 10;
        if let Some(vm200) = vd.non_insertion_variants.get_mut(&200) {
            vm200.map.extend(non_ins200.map);
        }

        let cfg = Config {
            sv_min_len: 50,
            compat: crate::config::CompatProfile {
                sv_respect_min_len: true,
                ..crate::config::COMPAT_VARDICTJAVA_1_8_3
            },
            ..Config::default()
        };
        let aligned = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg)
            .build_aligned_variants(&vd);
        let v50 = &aligned.get(&50).unwrap().variants[0];
        assert_ne!(
            v50.var_allele, "<DUP>",
            "48-bp complex insertion (< -L 50) must remain explicit nucleotide sequence"
        );
        let v120 = &aligned.get(&120).unwrap().variants[0];
        assert_eq!(v120.var_allele, "<DUP>");
        assert_eq!(v120.end_pos - v120.start_pos + 1, 80);
        let v200 = &aligned.get(&200).unwrap().variants[0];
        assert_eq!(v200.var_allele, "<INV>");
        assert_eq!(
            v200.ref_allele, "A",
            "isolated sv_respect_min_len must populate 1-bp REF on <inv-20> instead of empty string"
        );

        let cfg_off = Config {
            sv_min_len: 50,
            compat: crate::config::COMPAT_OFF,
            ..Config::default()
        };
        let aligned_off = ToVarsBuilder::new(&ref_genome, &region, &chr_lengths, &cfg_off)
            .build_aligned_variants(&vd);
        let v50_off = &aligned_off.get(&50).unwrap().variants[0];
        assert_ne!(
            v50_off.var_allele, "<DUP>",
            "48-bp complex insertion (< -L 50) under --compat=off must remain explicit nucleotide sequence"
        );
        let v120_off = &aligned_off.get(&120).unwrap().variants[0];
        assert_eq!(v120_off.var_allele, "<DUP>");
        assert_eq!(v120_off.end_pos - v120_off.start_pos, 80);
    }
}
