//! ToVarsBuilder parses intermediate VariationData and finalizes into Variant records.
//! Ports ToVarsBuilder.java.

use crate::config::Config;
use crate::reference::Reference;
use crate::regions::Region;
use crate::types::round_half_even;
use crate::variations::{Variant, VariationData};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct MsiResult {
    pub msi: f64,
    pub shift3: i32,
    pub msint: String,
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
) {
    tvref.description = description.to_string();
    tvref.pos_coverage = cnt.cnt;
    tvref.fwd = fwd;
    tvref.rev = rev;
    tvref.bias.clear();
    let _ = std::fmt::Write::write_fmt(&mut tvref.bias, format_args!("{}", bias));
    tvref.frequency = round_half_even(cnt.cnt as f64 / ttcov as f64, 4);
    tvref.mean_position = round_half_even(cnt.pmean / cnt.cnt as f64, 1);
    tvref.pstd = cnt.pstd;
    tvref.mean_quality = base_quality;
    tvref.qstd = cnt.qstd;
    tvref.mean_mapping_quality = mapping_quality;
    tvref.hq_lq_ratio = hicnt as f64 / if locnt != 0 { locnt as f64 } else { 0.5 };
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
    tvref.shift3 = 0;
    tvref.msi = 0.0;
    tvref.nm = round_half_even(cnt.nm / cnt.cnt as f64, 1);
    tvref.hicnt = hicnt;
    tvref.hicov = hicov;
    tvref.duprate = duprate;
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

impl<'a> ToVarsBuilder<'a> {
    /// Create a new ToVarsBuilder with reference, region, chromosomes_lengths and config.
    pub fn new(
        reference: &'a Reference,
        region: &'a Region,
        chromosomes_lengths: &'a HashMap<String, usize>,
        config: &'a Config,
    ) -> Self {
        let (primary_ref_seq, primary_ref_start) = match (
            reference.loaded_regions.first(),
            reference.region_sequences.first(),
        ) {
            (Some(lr), Some(bases))
                if lr.chr == region.chr && lr.sequence_end >= lr.sequence_start =>
            {
                let max_len = (lr.sequence_end.saturating_sub(lr.sequence_start) as usize)
                    .saturating_add(1)
                    .min(bases.len());
                (bases.get(..max_len).unwrap_or(&[]), lr.sequence_start)
            }
            _ => (&[][..], 1),
        };
        let chr_len = chromosomes_lengths.get(&region.chr).copied().unwrap_or(0) as i32;
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

    /// Convert VariationData into a map of position -> Vars container.
    pub fn build_aligned_variants(
        &self,
        var_data: &VariationData,
    ) -> crate::javamap::JavaHashMap<i32, crate::variations::Vars> {
        let mut aligned_variants = crate::javamap::JavaHashMap::new();
        let duprate = var_data.duprate;

        // Sparse copy-on-write overlay for positions mutated by create_insertion (position + 1);
        // lookup-only, never iterated.
        let mut overlay: rustc_hash::FxHashMap<i32, crate::variations::VariationMap> =
            rustc_hash::FxHashMap::default();

        let iter_positions: Vec<i32> = var_data.non_insertion_variants.keys().copied().collect();
        // Do not sort positions, use the JavaHashMap iteration order!
        if self.config.pileup || self.config.amplicon.is_some() || self.config.bam.contains('|') {
            aligned_variants.reserve_nodes(iter_positions.len());
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

        // Process non-insertion variants
        for position in iter_positions {
            // `var_map`'s shared NLL borrow of `overlay` and `var_data.non_insertion_variants`
            // ends at `create_variant` below, before `&mut overlay` is passed to `create_insertion`.
            let var_map = match overlay
                .get(&position)
                .or_else(|| var_data.non_insertion_variants.get(&position))
            {
                Some(m) => m,
                None => continue,
            };

            // Skip position if there are no variants on position
            if var_map.map.is_empty() && !var_data.insertion_variants.contains_key(&position) {
                continue;
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
                continue;
            }

            if (var_map.sv.is_none() || self.config.deldupvar)
                && ((position as u32) < self.region.start.0
                    || (position as u32) > self.region.end.0)
            {
                continue;
            }

            if var_map.sv.is_none() && !var_data.ref_coverage.contains_key(&position) {
                continue;
            }

            if self.is_the_same_variation_on_ref(position, var_map, &var_data.insertion_variants) {
                continue;
            }

            let total_cov = var_data.ref_coverage.get(&position).copied().unwrap_or(0);
            if total_cov == 0 {
                continue;
            }

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

            let mut current_total_cov = total_cov;
            current_total_cov = self.create_insertion(
                &var_data.insertion_variants,
                &var_data.non_insertion_variants,
                &mut overlay,
                &var_data.ref_coverage,
                duprate,
                position,
                current_total_cov,
                &mut var,
                &mut debug_lines,
                hicov,
            );

            self.sort_variants(&mut var);

            let maxfreq = self.collect_vars_at_position(&mut aligned_variants, position, &mut var);

            if !self.config.pileup
                && maxfreq <= self.config.freq
                && self.config.amplicon.is_none()
                && !self.config.bam.contains('|')
            {
                aligned_variants.remove(&position);
                continue;
            }

            let entry = crate::variations::get_or_put_vars(&mut aligned_variants, position);
            self.collect_reference_variants(
                position,
                current_total_cov,
                entry,
                &var_data.ref_coverage,
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
            log::warn!(
                "JavaHashMap treeified bin encountered for chr {}:{}-{}; key iteration order may diverge from Java LinkedList order",
                self.region.chr, self.region.start, self.region.end
            );
        }
        aligned_variants
    }

    fn is_the_same_variation_on_ref(
        &self,
        position: i32,
        vars_at_cur_position: &crate::variations::VariationMap,
        insertion_variants: &crate::javamap::JavaHashMap<i32, crate::variations::VariationMap>,
    ) -> bool {
        if self.config.pileup || self.config.bam.contains('|') || self.config.amplicon.is_some() {
            return false;
        }
        let has_ins = insertion_variants.contains_key(&position);
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
        if let Some(ref_base) = self.get_ref_base(position as u32) {
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
        let mut hicov = 0;
        for (key, val) in &non_insertion_variations.map {
            if key == "SV" || key.starts_with('+') {
                continue;
            }
            hicov += val.hicnt;
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
        let eff_start = start.max(1);
        let span = ((end as i64) - (eff_start as i64) + 1).clamp(0, 1 << 20) as usize;
        if out.capacity() == 0 {
            out.reserve_exact(span);
        } else {
            out.reserve(span);
        }
        for p in eff_start..=end {
            if let Some(b) = self.get_ref_base(p as u32) {
                if b.is_ascii() {
                    out.push(b as char);
                } else {
                    crate::util::ref_bases_into(out, self.reference, &self.region.chr, start, end);
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
        self.join_ref_into(leftseq, std::cmp::max(position - 70, 1), position - 1);
        let chr_len = self.chr_len;
        self.join_ref_into(
            tseq,
            position,
            std::cmp::min(position + dellen + 70, chr_len),
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
        self.join_ref_into(leftseq, std::cmp::max(position - 50, 1), position);
        let chr_len = self.chr_len;
        self.join_ref_into(tseq2, position + 1, std::cmp::min(position + 70, chr_len));

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
        for &description_string in keys {
            if description_string == "SV" {
                if let Some(ref sv) = non_insertion_variations.sv {
                    let entry = crate::variations::get_or_put_vars(aligned_vars, position);

                    entry.sv = format!("{}-{}-{}", sv.splits, sv.pairs, sv.clusters);
                }
                continue;
            }
            if let Some(cnt) = non_insertion_variations.map.get(description_string) {
                if cnt.cnt == 0 {
                    continue;
                }
                let fwd = cnt.fwd;
                let rev = cnt.rev;
                let bias = strand_bias(fwd, rev, self.config.min_reads_bias as i32);
                let base_quality = round_half_even(cnt.qmean / cnt.cnt as f64, 1);
                let mapping_quality = round_half_even(cnt.mqmean / cnt.cnt as f64, 1);
                let hicnt = cnt.hicnt;
                let locnt = cnt.locnt;

                let mut ttcov = total_pos_coverage;
                if cnt.cnt > total_pos_coverage
                    && cnt.extracnt > 0
                    && cnt.cnt - total_pos_coverage < cnt.extracnt
                {
                    ttcov = cnt.cnt;
                }

                let mut tvref = Variant::new();
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
                );

                if self.config.debug {
                    debug_lines.push(tvref.debug_variants_content());
                }

                var.push(tvref);
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
                    if let Some(&cov_next) = ref_coverage.get(&(position + 1)) {
                        total_pos_coverage = cov_next;
                        ttcov = cov_next;
                    }
                }
                if let Some(cnt) = insertion_variations.map.get(description_string) {
                    if cnt.cnt == 0 {
                        continue;
                    }
                    let fwd = cnt.fwd;
                    let rev = cnt.rev;
                    let bias = strand_bias(fwd, rev, self.config.min_reads_bias as i32);
                    let base_quality = round_half_even(cnt.qmean / cnt.cnt as f64, 1);
                    let mapping_quality = round_half_even(cnt.mqmean / cnt.cnt as f64, 1);
                    let hicnt = cnt.hicnt;
                    let locnt = cnt.locnt;

                    if cnt.cnt > total_pos_coverage
                        && cnt.extracnt != 0
                        && cnt.cnt - total_pos_coverage < cnt.extracnt
                    {
                        ttcov = cnt.cnt;
                    }

                    if ttcov < cnt.cnt {
                        ttcov = cnt.cnt;
                        if let Some(&cov_next) = ref_coverage.get(&(position + 1)) {
                            if ttcov < cov_next - cnt.cnt {
                                ttcov = cov_next;
                                if let Some(ref_base) = self.get_ref_base((position + 1) as u32) {
                                    let mut buf = [0u8; 4];
                                    let ref_base_str: &str =
                                        (ref_base as char).encode_utf8(&mut buf);
                                    let next_pos = position + 1;
                                    if let std::collections::hash_map::Entry::Vacant(e) =
                                        overlay.entry(next_pos)
                                    {
                                        if let Some(orig) = orig_non_ins.get(&next_pos) {
                                            e.insert(orig.clone());
                                        }
                                    }
                                    if let Some(next_non_ins) = overlay.get_mut(&next_pos) {
                                        if let Some(next_var) =
                                            next_non_ins.map.get_mut(ref_base_str)
                                        {
                                            next_var.fwd -= fwd;
                                            next_var.rev -= rev;
                                        }
                                    }
                                }
                            }
                        }
                        total_pos_coverage = ttcov;
                    }

                    if hicov < hicnt {
                        hicov = hicnt;
                    }

                    let mut tvref = Variant::new();
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
                    );

                    if self.config.debug {
                        // ports Variant.debugVariantsContentInsertion (Variant.java:261-266)
                        let mut s = String::from("I");
                        s.push_str(&tvref.debug_variants_content());
                        debug_lines.push(s);
                    }

                    var.push(tvref);
                }
            }
        }
        total_pos_coverage
    }

    fn collect_vars_at_position(
        &self,
        aligned_variants: &mut crate::javamap::JavaHashMap<i32, crate::variations::Vars>,
        position: i32,
        var: &mut Vec<Variant>,
    ) -> f64 {
        let mut maxfreq = 0.0;
        let ref_base = self.get_ref_base(position as u32).unwrap_or(b'N');
        let mut buf = [0u8; 4];
        let ref_base_str: &str = (ref_base as char).encode_utf8(&mut buf);

        // H1: bound per-element Vec::push regrowth (dhat: ~8% on deep
        // pileup). The entry starts empty on this path (one drain per
        // position per build), so a single reserve of the exact non-ref row
        // count covers every push below. Guarded by non_ref_rows > 0, so
        // entry creation timing is unchanged: the drain loop creates the
        // entry on its first row exactly as before (a hit here performs no
        // structural change), and an empty `var` creates nothing.
        let non_ref_rows = var
            .iter()
            .filter(|t| t.description.as_str() != ref_base_str)
            .count();
        let mut reserved = false;

        for tvar in var.drain(..) {
            if tvar.description == ref_base_str {
                let entry = crate::variations::get_or_put_vars(aligned_variants, position);

                entry.ref_variant = Some(tvar);
            } else {
                let entry = crate::variations::get_or_put_vars(aligned_variants, position);
                if !reserved {
                    entry.variants.reserve(non_ref_rows);
                    reserved = true;
                }

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
            if let Some(next_non_ins) = overlay
                .get(&(position + 1))
                .or_else(|| orig_non_ins.get(&(position + 1)))
            {
                if let Some(ref_base_next) = self.get_ref_base((position + 1) as u32) {
                    let mut buf = [0u8; 4];
                    let ref_base_next_str: &str = (ref_base_next as char).encode_utf8(&mut buf);
                    if let Some(tpref) = next_non_ins.map.get(ref_base_next_str) {
                        reference_forward_coverage_adj = tpref.fwd;
                        reference_reverse_coverage_adj = tpref.rev;
                    }
                }
            }
        }

        let mut positions_for_changed_ref_variant = Vec::new();

        if !variations_at_pos.variants.is_empty() {
            let format_plus_genotype1 = |desc: &str| -> String {
                if let Some(dup_val) = get_dup_num(desc) {
                    format!("+{}", 50 + dup_val as u32)
                } else {
                    format!("+{}", desc.len() as i32 - 1)
                }
            };

            let genotype1_current: std::borrow::Cow<'_, str> =
                if let Some(ref ref_var) = variations_at_pos.ref_variant {
                    if ref_var.frequency >= self.config.freq {
                        if ref_var.description.starts_with('+') {
                            std::borrow::Cow::Owned(format_plus_genotype1(&ref_var.description))
                        } else {
                            std::borrow::Cow::Borrowed(ref_var.description.as_str())
                        }
                    } else if let Some(first_var) = variations_at_pos.variants.first() {
                        if first_var.description.starts_with('+') {
                            std::borrow::Cow::Owned(format_plus_genotype1(&first_var.description))
                        } else {
                            std::borrow::Cow::Owned(first_var.description.clone())
                        }
                    } else if ref_var.description.starts_with('+') {
                        std::borrow::Cow::Owned(format_plus_genotype1(&ref_var.description))
                    } else {
                        std::borrow::Cow::Borrowed(ref_var.description.as_str())
                    }
                } else if let Some(first_var) = variations_at_pos.variants.first() {
                    if first_var.description.starts_with('+') {
                        std::borrow::Cow::Owned(format_plus_genotype1(&first_var.description))
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
                let mut genotype2: std::borrow::Cow<'_, str> = if vref.description.starts_with('+')
                {
                    std::borrow::Cow::Owned(format!("+{}", vref.description.len() as i32 - 1))
                } else {
                    std::borrow::Cow::Borrowed(vref.description.as_str())
                };

                let description_string = &vref.description;
                let deletion_length = get_begin_minus_number(description_string).unwrap_or(0);

                let mut end_position = position;
                if description_string.starts_with('-') {
                    end_position = position + deletion_length - 1;
                }

                let mut refallele = String::new();
                let mut varallele;

                let mut shift3 = 0;
                let mut msi = 0.0;
                let mut msint_len = 0;

                let mut start_position = position;

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
                    if self.config.three_prime {
                        start_position += shift3;
                        end_position += shift3;
                    }
                    let ref_base_byte = self.get_ref_base(position as u32).unwrap_or(b'N');
                    refallele.push(ref_base_byte as char);
                    varallele = String::with_capacity(1 + stripped_desc.len());
                    varallele.push(ref_base_byte as char);
                    varallele.push_str(stripped_desc);

                    if varallele.len() > self.config.sv_min_len as usize {
                        end_position += varallele.len() as i32;
                        varallele = "<DUP>".to_string();
                    }
                    if let Some(dup_val) = get_dup_num(&varallele) {
                        end_position = start_position + (2 * 50 + dup_val) - 1;
                        genotype2 = std::borrow::Cow::Owned(format!("+{}", 2 * 50 + dup_val));
                        varallele = "<DUP>".to_string();
                    }
                } else if description_string.starts_with('-') {
                    if deletion_length < self.config.sv_min_len as i32 {
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

                        if get_inv_num(description_string).is_some() {
                            varallele = "<INV>".to_string();
                            genotype2 =
                                std::borrow::Cow::Owned(format!("<INV{}>", deletion_length));
                        }
                    } else if matches_begin_minus_number_caret(description_string) {
                        varallele = "<INV>".to_string();
                        genotype2 = std::borrow::Cow::Owned(format!("<INV{}>", deletion_length));
                    } else {
                        varallele = "<DEL>".to_string();
                    }

                    if !description_string.contains('&')
                        && !description_string.contains('#')
                        && !description_string.contains('^')
                    {
                        if self.config.three_prime {
                            start_position += shift3;
                        }
                        if varallele != "<DEL>" {
                            let base_prev =
                                self.get_ref_base((position - 1) as u32).unwrap_or(b'N');
                            varallele.clear();
                            varallele.push(base_prev as char);
                        }
                        let base_prev = self.get_ref_base((position - 1) as u32).unwrap_or(b'N');
                        refallele.push(base_prev as char);
                        start_position -= 1;
                    }

                    if matches_some_sv_numbers(description_string) {
                        let base_cur = self.get_ref_base(position as u32).unwrap_or(b'N');
                        refallele.clear();
                        refallele.push(base_cur as char);
                    } else if deletion_length < self.config.sv_min_len as i32 {
                        self.join_ref_into(scratch1, position, position + deletion_length - 1);
                        refallele.push_str(scratch1);
                    }
                } else {
                    self.join_ref_into(scratch1, std::cmp::max(position - 30, 1), position + 1);
                    self.join_ref_into(
                        scratch2,
                        position + 2,
                        std::cmp::min(position + 70, chr_len),
                    );

                    (msi, shift3, msint_len) =
                        find_msi_metrics(scratch1.as_bytes(), scratch2.as_bytes(), &[]);

                    let base_cur = self.get_ref_base(position as u32).unwrap_or(b'N');
                    refallele.push(base_cur as char);
                    varallele = description_string.clone();
                }

                if let Some(extra) = get_amp_atgc(description_string) {
                    varallele = varallele.replacen('&', "", 1);
                    self.join_ref_into(
                        scratch1,
                        end_position + 1,
                        end_position + extra.len() as i32,
                    );
                    refallele.push_str(scratch1);
                    if !scratch1.is_empty() {
                        genotype1_current_inner.to_mut().push_str(scratch1);
                    }
                    end_position += extra.len() as i32;

                    if let Some(vextra) = get_amp_atgc(&varallele) {
                        let vextra_len = vextra.len() as i32;
                        varallele = varallele.replacen('&', "", 1);
                        self.join_ref_into(scratch1, end_position + 1, end_position + vextra_len);
                        refallele.push_str(scratch1);
                        if !scratch1.is_empty() {
                            genotype1_current_inner.to_mut().push_str(scratch1);
                        }
                        end_position += vextra_len;
                    }

                    if description_string.starts_with('+') {
                        refallele = refallele.get(1..).unwrap_or("").to_string();
                        varallele = varallele.get(1..).unwrap_or("").to_string();
                        start_position += 1;
                    }

                    if varallele == "<DEL>" && !refallele.is_empty() {
                        let ref_base = self.get_ref_base(start_position as u32).unwrap_or(b'N');
                        refallele.clear();
                        refallele.push(ref_base as char);
                        if let Some(&cov_prev) = ref_coverage.get(&(start_position - 1)) {
                            total_pos_coverage_adj = cov_prev;
                        }
                        if vref.pos_coverage > total_pos_coverage_adj {
                            total_pos_coverage_adj = vref.pos_coverage;
                        }
                        vref.frequency = vref.pos_coverage as f64 / total_pos_coverage_adj as f64;
                    }
                }

                if let Some((matched_sequence, tail)) =
                    get_hash_group_caret_group(description_string)
                {
                    end_position += matched_sequence.len() as i32;
                    self.join_ref_into(
                        scratch1,
                        end_position - matched_sequence.len() as i32 + 1,
                        end_position,
                    );
                    refallele.push_str(scratch1);

                    if let Some(deletion) = get_begin_digits(tail) {
                        self.join_ref_into(scratch1, end_position + 1, end_position + deletion);
                        refallele.push_str(scratch1);
                        end_position += deletion;
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
                    varallele = varallele.replacen('^', "", 1);
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
                        let mut n = 0;
                        while refallele.len() > n + 1
                            && varallele.len() > n + 1
                            && refallele.as_bytes().get(n) == varallele.as_bytes().get(n)
                        {
                            n += 1;
                        }
                        if n != 0 {
                            start_position += n as i32;
                            refallele = refallele.get(n..).unwrap_or("").to_string();
                            varallele = varallele.get(n..).unwrap_or("").to_string();
                        }
                    }

                    if cut_site != 0
                        && refallele.len() != varallele.len()
                        && refallele.starts_with(varallele.chars().next().unwrap_or('\0'))
                        && !(start_position == cut_site || end_position == cut_site)
                    {
                        let mut n = 0;
                        let dis = std::cmp::min(
                            (cut_site - start_position).abs(),
                            (cut_site - end_position).abs(),
                        );
                        if start_position < cut_site {
                            while start_position + n < cut_site
                                && n < shift3
                                && end_position + n != cut_site
                            {
                                n += 1;
                            }
                            if (start_position + n - cut_site).abs() > dis
                                && (end_position + n - cut_site).abs() > dis
                            {
                                n = 0;
                            }
                        }
                        if end_position < cut_site
                            && n == 0
                            && (end_position - cut_site).abs() <= (start_position - cut_site).abs()
                        {
                            while end_position + n < cut_site && n < shift3 {
                                n += 1;
                            }
                        }
                        if n > 0 {
                            start_position += n;
                            end_position += n;
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
                                self.get_ref_base(start_position as u32).unwrap_or(b'N') as char;
                            varallele = format!("{}{}", first_char, tva);
                        }
                        vref.crispr = n;
                    }
                }

                self.join_ref_into(
                    &mut vref.left_seq,
                    std::cmp::max(start_position - 20, 1),
                    start_position - 1,
                );
                self.join_ref_into(
                    &mut vref.right_seq,
                    end_position + 1,
                    std::cmp::min(end_position + 20, chr_len),
                );

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

                vref.extra_freq = round_half_even(vref.extra_freq, 4);
                vref.frequency = round_half_even(vref.frequency, 4);
                vref.hi_freq = round_half_even(vref.hi_freq, 4);
                vref.msi = round_half_even(msi, 3);
                vref.msint = msint_len;
                vref.shift3 = shift3;
                vref.start_pos = start_position;
                vref.end_pos = end_position;
                validate_ref_allele_in_place(&mut refallele);
                vref.ref_allele = refallele;
                vref.var_allele = varallele;
                vref.genotype = Some(genotype);
                vref.total_pos_cov = total_pos_coverage_adj;
                vref.ref_fwd_cov = reference_forward_coverage_adj;
                vref.ref_rev_cov = reference_reverse_coverage_adj;

                let ref_bias = variations_at_pos
                    .ref_variant
                    .as_ref()
                    .map_or("0", |ref_var| ref_var.bias.as_str());
                vref.bias.reserve(ref_bias.len().saturating_add(1));
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
                    positions_for_changed_ref_variant.push(position);
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
            self.update_ref_variant(
                position,
                total_pos_coverage_adj,
                vref,
                reference_forward_coverage_adj,
                reference_reverse_coverage_adj,
                duprate,
                debug_lines,
            );
        } else {
            variations_at_pos.ref_variant = Some(Variant::new());
        }

        if variations_at_pos.ref_variant.is_some()
            && self.config.pileup
            && (positions_for_changed_ref_variant.contains(&position)
                || self.config.amplicon.is_some())
        {
            if let Some(ref mut vref) = variations_at_pos.ref_variant {
                self.update_ref_variant(
                    position,
                    total_pos_coverage_adj,
                    vref,
                    reference_forward_coverage_adj,
                    reference_reverse_coverage_adj,
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
        reference_forward_coverage: i32,
        reference_reverse_coverage: i32,
        duprate: f64,
        debug_lines: &[String],
    ) {
        vref.total_pos_cov = total_pos_coverage;
        vref.pos_coverage = 0;
        vref.frequency = 0.0;
        vref.ref_fwd_cov = reference_forward_coverage;
        vref.ref_rev_cov = reference_reverse_coverage;
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
        vref.hi_freq = round_half_even(vref.hi_freq, 4);

        let ref_base = self.get_ref_base(position as u32).unwrap_or(b'N');
        let valid_base = match ref_base {
            b'M' | b'R' | b'W' | b'V' | b'H' | b'D' => b'A',
            b'S' | b'Y' | b'B' => b'C',
            b'K' => b'G',
            b'm' | b'r' | b'w' | b'v' | b'h' | b'd' => b'a',
            b's' | b'y' | b'b' => b'c',
            b'k' => b'g',
            _ => ref_base,
        } as char;

        vref.ref_allele.clear();
        vref.ref_allele.push(valid_base);
        vref.var_allele.clear();
        vref.var_allele.push(valid_base);
        let base_ch = ref_base as char;
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

/// Strand bias calculation matching VarDict
pub fn strand_bias(fwd: i32, rev: i32, min_reads_bias: i32) -> i32 {
    crate::variations::util::strand_bias(fwd, rev, 0.05, min_reads_bias)
}

fn validate_ref_allele_in_place(ref_allele: &mut String) {
    if !ref_allele.bytes().any(|b| {
        matches!(
            b,
            b'M' | b'R'
                | b'W'
                | b'V'
                | b'H'
                | b'D'
                | b'S'
                | b'Y'
                | b'B'
                | b'K'
                | b'm'
                | b'r'
                | b'w'
                | b'v'
                | b'h'
                | b'd'
                | b's'
                | b'y'
                | b'b'
                | b'k'
        )
    }) {
        return;
    }
    let mut bytes = std::mem::take(ref_allele).into_bytes();
    for b in &mut bytes {
        *b = match *b {
            b'M' | b'R' | b'W' | b'V' | b'H' | b'D' => b'A',
            b'S' | b'Y' | b'B' => b'C',
            b'K' => b'G',
            b'm' | b'r' | b'w' | b'v' | b'h' | b'd' => b'a',
            b's' | b'y' | b'b' => b'c',
            b'k' => b'g',
            other => other,
        };
    }
    if let Ok(s) = String::from_utf8(bytes) {
        *ref_allele = s;
    }
}

/// Validate reference allele sequence according to VCF 4.3 IUPAC ambiguity rules.
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

fn get_dup_num(s: &str) -> Option<i32> {
    if let Some(idx) = s.find("<dup") {
        let rest = s.get(idx + 4..)?;
        let len = rest.bytes().take_while(|b| b.is_ascii_digit()).count();
        if len > 0 {
            return rest.get(..len)?.parse().ok();
        }
    }
    None
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

fn get_inv_num(s: &str) -> Option<i32> {
    if let Some(idx) = s.find("<inv") {
        let rest = s.get(idx + 4..)?;
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

fn matches_some_sv_numbers(s: &str) -> bool {
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
    }
}
