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
    tvref.bias = format!("{}", bias);
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

impl<'a> ToVarsBuilder<'a> {
    /// Create a new ToVarsBuilder with reference, region, chromosomes_lengths and config.
    pub fn new(
        reference: &'a Reference,
        region: &'a Region,
        chromosomes_lengths: &'a HashMap<String, usize>,
        config: &'a Config,
    ) -> Self {
        Self {
            reference,
            region,
            chromosomes_lengths,
            config,
        }
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

        let mut var = Vec::new();
        let mut keys = Vec::new();
        let mut debug_lines = Vec::new();

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

            let hicov = self.calc_hicov(var_data.insertion_variants.get(&position), var_map);

            var.clear();
            keys.clear();
            keys.extend(var_map.map.keys().cloned());
            if var_map.sv.is_some() {
                keys.push("SV".to_string());
            }
            keys.sort();

            debug_lines.clear();

            self.create_variant(
                duprate,
                &mut aligned_variants,
                position,
                var_map,
                total_cov,
                &mut var,
                &mut debug_lines,
                &keys,
                hicov,
            );

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

            if let Some(entry) = aligned_variants.get_mut(&position) {
                self.collect_reference_variants(
                    position,
                    current_total_cov,
                    entry,
                    &var_data.ref_coverage,
                    &var_data.non_insertion_variants,
                    &overlay,
                    duprate,
                    &debug_lines,
                );
            }
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
                Some(k) => k.as_str(),
                None => return false,
            },
            (0, true) => "I",
            (1, true) => match vars_at_cur_position.map.keys().next() {
                Some(k) if k == "I" => "I",
                _ => return false,
            },
            _ => return false,
        };
        if let Some(ref_base) = self.reference.get_base(&self.region.chr, position as u32) {
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
    fn join_ref(&self, start: i32, end: i32) -> String {
        String::from_utf8_lossy(&crate::util::ref_bases_with(
            self.reference,
            &self.region.chr,
            start,
            end,
            None,
        ))
        .into_owned()
    }

    fn proceed_vref_is_deletion(&self, position: i32, dellen: i32) -> MsiResult {
        let leftseq = self.join_ref(std::cmp::max(position - 70, 1), position - 1);
        let chr_len = self
            .chromosomes_lengths
            .get(&self.region.chr)
            .copied()
            .unwrap_or(0) as i32;
        let tseq = self.join_ref(position, std::cmp::min(position + dellen + 70, chr_len));

        let tseq_bytes = tseq.as_bytes();
        let tseq1 = if tseq_bytes.len() >= dellen as usize {
            &tseq_bytes[..dellen as usize]
        } else {
            tseq_bytes
        };
        let tseq2 = if tseq_bytes.len() >= dellen as usize {
            &tseq_bytes[dellen as usize..]
        } else {
            &[]
        };

        let msi_res = find_msi(tseq1, tseq2, leftseq.as_bytes());
        let mut msi = msi_res.msi;
        let mut msint = msi_res.msint;
        let shift3 = msi_res.shift3;

        let msi_without_left = find_msi(leftseq.as_bytes(), tseq2, &[]);
        let tmsi = msi_without_left.msi;
        let tmsint = msi_without_left.msint;
        if msi < tmsi {
            msi = tmsi;
            msint = tmsint;
        }
        if dellen > 0 && msi <= (shift3 as f64) / (dellen as f64) {
            msi = (shift3 as f64) / (dellen as f64);
        }

        MsiResult { msi, shift3, msint }
    }

    fn proceed_vref_is_insertion(&self, position: i32, vn: &str) -> MsiResult {
        let tseq1 = vn.strip_prefix('+').unwrap_or(vn);
        let leftseq = self.join_ref(std::cmp::max(position - 50, 1), position);
        let chr_len = self
            .chromosomes_lengths
            .get(&self.region.chr)
            .copied()
            .unwrap_or(0) as i32;
        let tseq2 = self.join_ref(position + 1, std::cmp::min(position + 70, chr_len));

        let msi_res = find_msi(tseq1.as_bytes(), tseq2.as_bytes(), leftseq.as_bytes());
        let mut msi = msi_res.msi;
        let mut msint = msi_res.msint;
        let shift3 = msi_res.shift3;

        let msi_without_left = find_msi(leftseq.as_bytes(), tseq2.as_bytes(), &[]);
        let tmsi = msi_without_left.msi;
        let tmsint = msi_without_left.msint;
        if msi < tmsi {
            msi = tmsi;
            msint = tmsint;
        }
        if !tseq1.is_empty() && msi <= (shift3 as f64) / (tseq1.len() as f64) {
            msi = (shift3 as f64) / (tseq1.len() as f64);
        }

        MsiResult { msi, shift3, msint }
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
        keys: &[String],
        hicov: i32,
    ) {
        for description_string in keys {
            if description_string == "SV" {
                if let Some(ref sv) = non_insertion_variations.sv {
                    let entry = aligned_vars
                        .compute_if_absent(position, |_| crate::variations::Vars::new());

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
            let mut insertion_description_strings: Vec<String> =
                insertion_variations.map.keys().cloned().collect();
            insertion_description_strings.sort();

            for description_string in &insertion_description_strings {
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
                                if let Some(ref_base) = self
                                    .reference
                                    .get_base(&self.region.chr, (position + 1) as u32)
                                {
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
                                            next_var.fwd = (next_var.fwd - fwd).max(0);
                                            next_var.rev = (next_var.rev - rev).max(0);
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
        let ref_base = self
            .reference
            .get_base(&self.region.chr, position as u32)
            .unwrap_or(b'N');
        let mut buf = [0u8; 4];
        let ref_base_str: &str = (ref_base as char).encode_utf8(&mut buf);

        for tvar in var.drain(..) {
            if tvar.description == ref_base_str {
                let entry = aligned_variants
                    .compute_if_absent(position, |_| crate::variations::Vars::new());

                entry.ref_variant = Some(tvar);
            } else {
                let entry = aligned_variants
                    .compute_if_absent(position, |_| crate::variations::Vars::new());

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
    ) {
        let mut reference_forward_coverage = 0;
        let mut reference_reverse_coverage = 0;

        let genotype1 = if let Some(ref ref_var) = variations_at_pos.ref_variant {
            if ref_var.frequency >= self.config.freq {
                ref_var.description.clone()
            } else if let Some(first_var) = variations_at_pos.variants.first() {
                first_var.description.clone()
            } else {
                ref_var.description.clone()
            }
        } else if let Some(first_var) = variations_at_pos.variants.first() {
            first_var.description.clone()
        } else {
            String::new()
        };

        if let Some(ref ref_var) = variations_at_pos.ref_variant {
            reference_forward_coverage = ref_var.fwd;
            reference_reverse_coverage = ref_var.rev;
        }

        let mut genotype1_current = genotype1;
        if genotype1_current.starts_with('+') {
            if let Some(dup_val) = get_dup_num(&genotype1_current) {
                genotype1_current = format!("+{}", 50 + dup_val as u32);
            } else {
                genotype1_current = format!("+{}", genotype1_current.len() as i32 - 1);
            }
        }

        let mut reference_forward_coverage_adj = reference_forward_coverage;
        let mut reference_reverse_coverage_adj = reference_reverse_coverage;
        let mut total_pos_coverage_adj = total_pos_coverage;

        if total_pos_coverage_adj > ref_coverage.get(&position).copied().unwrap_or(0) {
            if let Some(next_non_ins) = overlay
                .get(&(position + 1))
                .or_else(|| orig_non_ins.get(&(position + 1)))
            {
                if let Some(ref_base_next) = self
                    .reference
                    .get_base(&self.region.chr, (position + 1) as u32)
                {
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
            for vref in &mut variations_at_pos.variants {
                let mut genotype1_current_inner = genotype1_current.clone();
                let mut genotype2 = vref.description.clone();
                if genotype2.starts_with('+') {
                    genotype2 = format!("+{}", genotype2.len() as i32 - 1);
                }

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
                let mut msint = String::new();

                let mut start_position = position;

                if let Some(stripped_desc) = description_string.strip_prefix('+') {
                    if !description_string.contains('&')
                        && !description_string.contains('#')
                        && !description_string.contains("<dup")
                    {
                        let ref_variant_msi =
                            self.proceed_vref_is_insertion(position, description_string);
                        msi = ref_variant_msi.msi;
                        shift3 = ref_variant_msi.shift3;
                        msint = ref_variant_msi.msint;
                    }
                    if self.config.three_prime {
                        start_position += shift3;
                        end_position += shift3;
                    }
                    let ref_base_byte = self
                        .reference
                        .get_base(&self.region.chr, position as u32)
                        .unwrap_or(b'N');
                    refallele = (ref_base_byte as char).to_string();
                    varallele = format!("{}{}", refallele, stripped_desc);

                    if varallele.len() > self.config.sv_min_len as usize {
                        end_position += varallele.len() as i32;
                        varallele = "<DUP>".to_string();
                    }
                    if let Some(dup_val) = get_dup_num(&varallele) {
                        end_position = start_position + (2 * 50 + dup_val) - 1;
                        genotype2 = format!("+{}", 2 * 50 + dup_val);
                        varallele = "<DUP>".to_string();
                    }
                } else if description_string.starts_with('-') {
                    if deletion_length < self.config.sv_min_len as i32 {
                        // replace first -digits with empty string
                        let re = format!("-{}", deletion_length);
                        varallele = description_string.replacen(&re, "", 1);

                        let ref_variant_msi =
                            self.proceed_vref_is_deletion(position, deletion_length);
                        msi = ref_variant_msi.msi;
                        shift3 = ref_variant_msi.shift3;
                        msint = ref_variant_msi.msint;

                        if get_inv_num(description_string).is_some() {
                            varallele = "<INV>".to_string();
                            genotype2 = format!("<INV{}>", deletion_length);
                        }
                    } else if matches_begin_minus_number_caret(description_string) {
                        varallele = "<INV>".to_string();
                        genotype2 = format!("<INV{}>", deletion_length);
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
                            let base_prev = self
                                .reference
                                .get_base(&self.region.chr, (position - 1) as u32)
                                .unwrap_or(b'N');
                            varallele = (base_prev as char).to_string();
                        }
                        let base_prev = self
                            .reference
                            .get_base(&self.region.chr, (position - 1) as u32)
                            .unwrap_or(b'N');
                        refallele = (base_prev as char).to_string();
                        start_position -= 1;
                    }

                    if matches_some_sv_numbers(description_string) {
                        let base_cur = self
                            .reference
                            .get_base(&self.region.chr, position as u32)
                            .unwrap_or(b'N');
                        refallele = (base_cur as char).to_string();
                    } else if deletion_length < self.config.sv_min_len as i32 {
                        refallele += &self.join_ref(position, position + deletion_length - 1);
                    }
                } else {
                    let tseq1 = self.join_ref(std::cmp::max(position - 30, 1), position + 1);
                    let chr_len = self
                        .chromosomes_lengths
                        .get(&self.region.chr)
                        .copied()
                        .unwrap_or(0) as i32;
                    let tseq2 = self.join_ref(position + 2, std::cmp::min(position + 70, chr_len));

                    let msi_res = find_msi(tseq1.as_bytes(), tseq2.as_bytes(), &[]);
                    msi = msi_res.msi;
                    shift3 = msi_res.shift3;
                    msint = msi_res.msint;

                    let base_cur = self
                        .reference
                        .get_base(&self.region.chr, position as u32)
                        .unwrap_or(b'N');
                    refallele = (base_cur as char).to_string();
                    varallele = description_string.clone();
                }

                if let Some(extra) = get_amp_atgc(description_string) {
                    varallele = varallele.replacen('&', "", 1);
                    let tch = self.join_ref(end_position + 1, end_position + extra.len() as i32);
                    refallele += &tch;
                    genotype1_current_inner += &tch;
                    end_position += extra.len() as i32;

                    if let Some(vextra) = get_amp_atgc(&varallele) {
                        varallele = varallele.replacen('&', "", 1);
                        let tch =
                            self.join_ref(end_position + 1, end_position + vextra.len() as i32);
                        refallele += &tch;
                        genotype1_current_inner += &tch;
                        end_position += vextra.len() as i32;
                    }

                    if description_string.starts_with('+') {
                        refallele = refallele.get(1..).unwrap_or("").to_string();
                        varallele = varallele.get(1..).unwrap_or("").to_string();
                        start_position += 1;
                    }

                    if varallele == "<DEL>" && refallele.len() >= 1 {
                        let ref_base = self
                            .reference
                            .get_base(&self.region.chr, start_position as u32)
                            .unwrap_or(b'N');
                        refallele = (ref_base as char).to_string();
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
                    refallele += &self.join_ref(
                        end_position - matched_sequence.len() as i32 + 1,
                        end_position,
                    );

                    if let Some(deletion) = get_begin_digits(tail) {
                        refallele += &self.join_ref(end_position + 1, end_position + deletion);
                        end_position += deletion;
                    }

                    // replaceFirst("#", "").replaceFirst("\\^(\\d+)?", "") in Java:
                    let cleaned_varallele = varallele.replacen('#', "", 1);
                    let cleaned_varallele = if let Some(idx) = cleaned_varallele.find('^') {
                        let rest = &cleaned_varallele[idx + 1..];
                        let digits_len = rest.chars().take_while(|c| c.is_ascii_digit()).count();
                        let mut new_s = cleaned_varallele[..idx].to_string();
                        new_s.push_str(&rest[digits_len..]);
                        new_s
                    } else {
                        cleaned_varallele
                    };
                    varallele = cleaned_varallele;

                    genotype1_current_inner = genotype1_current_inner
                        .replacen('#', "m", 1)
                        .replacen('^', "i", 1);
                    genotype2 = genotype2.replacen('#', "m", 1).replacen('^', "i", 1);
                }

                if get_caret_atgnc(description_string).is_some() {
                    varallele = varallele.replacen('^', "", 1);
                    genotype1_current_inner = genotype1_current_inner.replacen('^', "i", 1);
                    genotype2 = genotype2.replacen('^', "i", 1);
                }

                // CRISPR Site adjustment
                if let Some(cut_site) = self.config.crispr_cutting_site {
                    if cut_site != 0 && refallele.len() > 1 && varallele.len() > 1 {
                        let mut n = 0;
                        while refallele.len() > n + 1
                            && varallele.len() > n + 1
                            && refallele.as_bytes()[n] == varallele.as_bytes()[n]
                        {
                            n += 1;
                        }
                        if n != 0 {
                            start_position += n as i32;
                            refallele = refallele[n..].to_string();
                            varallele = varallele[n..].to_string();
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
                            refallele = self.join_ref(start_position, end_position);
                            let mut tva = String::new();
                            if refallele.len() < varallele.len() {
                                tva = varallele.get(1..).unwrap_or("").to_string();
                                if tva.len() > 1 {
                                    let ttn = (n as usize) % tva.len();
                                    if ttn != 0 {
                                        tva = format!("{}{}", &tva[ttn..], &tva[..ttn]);
                                    }
                                }
                            }
                            let first_char =
                                self.reference
                                    .get_base(&self.region.chr, start_position as u32)
                                    .unwrap_or(b'N') as char;
                            varallele = format!("{}{}", first_char, tva);
                        }
                        vref.crispr = n;
                    }
                }

                vref.left_seq =
                    self.join_ref(std::cmp::max(start_position - 20, 1), start_position - 1);
                let chr_len = self
                    .chromosomes_lengths
                    .get(&self.region.chr)
                    .copied()
                    .unwrap_or(0) as i32;
                vref.right_seq =
                    self.join_ref(end_position + 1, std::cmp::min(end_position + 20, chr_len));

                let mut genotype = format!("{}/{}", genotype1_current_inner, genotype2);
                genotype = genotype.replace(['&', '#'], "").replace('^', "i");

                vref.extra_freq = round_half_even(vref.extra_freq, 4);
                vref.frequency = round_half_even(vref.frequency, 4);
                vref.hi_freq = round_half_even(vref.hi_freq, 4);
                vref.msi = round_half_even(msi, 3);
                vref.msint = msint.len() as i32;
                vref.shift3 = shift3;
                vref.start_pos = start_position;
                vref.end_pos = end_position;
                vref.ref_allele = validate_ref_allele(&refallele);
                vref.var_allele = varallele;
                vref.genotype = Some(genotype);
                vref.total_pos_cov = total_pos_coverage_adj;
                vref.ref_fwd_cov = reference_forward_coverage_adj;
                vref.ref_rev_cov = reference_reverse_coverage_adj;

                let bias_final = if let Some(ref ref_var) = variations_at_pos.ref_variant {
                    format!("{};{}", ref_var.bias, vref.bias)
                } else {
                    format!("0;{}", vref.bias)
                };
                vref.bias = bias_final;

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

            // Rebuild var_desc_to_variant from finalized variants (total_pos_cov now set)
            variations_at_pos.var_desc_to_variant.clear();
            for v in &variations_at_pos.variants {
                variations_at_pos
                    .var_desc_to_variant
                    .insert(v.description.clone(), v.clone());
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
            vref.bias = format!("{};0", vref.bias);
        }
        vref.shift3 = 0;
        vref.start_pos = position;
        vref.end_pos = position;
        vref.hi_freq = round_half_even(vref.hi_freq, 4);

        let ref_base = self
            .reference
            .get_base(&self.region.chr, position as u32)
            .unwrap_or(b'N');
        let ref_base_str = (ref_base as char).to_string();

        vref.ref_allele = validate_ref_allele(&ref_base_str);
        vref.var_allele = validate_ref_allele(&ref_base_str);
        vref.genotype = Some(format!("{}/{}", ref_base_str, ref_base_str));
        vref.left_seq = String::new();
        vref.right_seq = String::new();
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

/// Validate reference allele sequence according to VCF 4.3 IUPAC ambiguity rules.
pub fn validate_ref_allele(ref_allele: &str) -> String {
    ref_allele
        .chars()
        .map(|c| match c {
            'M' | 'R' | 'W' | 'V' | 'H' | 'D' => 'A',
            'S' | 'Y' | 'B' => 'C',
            'K' => 'G',
            'm' | 'r' | 'w' | 'v' | 'h' | 'd' => 'a',
            's' | 'y' | 'b' => 'c',
            'k' => 'g',
            _ => c,
        })
        .collect()
}

/// Find microsatellite instability (MSI) motif and coordinate shift bounds.
pub fn find_msi(tseq1: &[u8], tseq2: &[u8], left: &[u8]) -> MsiResult {
    let mut msicnt = 0.0;
    let mut maxmsi = Vec::new();

    let tseq1_len = tseq1.len();
    let limit = std::cmp::min(tseq1_len, 6);

    for nmsi in 1..=limit {
        let msint = &tseq1[tseq1_len - nmsi..];

        let count = if !left.is_empty() {
            let left_len = left.len();
            if left_len + tseq1_len <= 256 {
                let mut combined = [0u8; 256];
                combined[..left_len].copy_from_slice(left);
                combined[left_len..left_len + tseq1_len].copy_from_slice(tseq1);
                count_repetitions_at_end(&combined[..left_len + tseq1_len], msint)
            } else {
                let mut combined = Vec::with_capacity(left_len + tseq1_len);
                combined.extend_from_slice(left);
                combined.extend_from_slice(tseq1);
                count_repetitions_at_end(&combined, msint)
            }
        } else {
            count_repetitions_at_end(tseq1, msint)
        };

        let mut curmsi = count as f64;
        curmsi += count_repetitions_at_start(tseq2, msint) as f64;

        if curmsi > msicnt {
            maxmsi = msint.to_vec();
            msicnt = curmsi;
        }
    }

    let mut shift3 = 0;
    if tseq1_len > 0 {
        while shift3 < tseq2.len() {
            let tseq_char = if shift3 < tseq1_len {
                tseq1[shift3]
            } else {
                tseq2[shift3 - tseq1_len]
            };
            if tseq_char == tseq2[shift3] {
                shift3 += 1;
            } else {
                break;
            }
        }
    }

    MsiResult {
        msi: msicnt,
        shift3: shift3 as i32,
        msint: String::from_utf8_lossy(&maxmsi).into_owned(),
    }
}

fn count_repetitions_at_end(seq: &[u8], pattern: &[u8]) -> usize {
    if pattern.is_empty() {
        return 0;
    }
    let mut count = 0;
    let mut current = seq;
    while current.ends_with(pattern) {
        count += 1;
        current = &current[..current.len() - pattern.len()];
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
        current = &current[pattern.len()..];
    }
    count
}

fn get_dup_num(s: &str) -> Option<i32> {
    if let Some(idx) = s.find("<dup") {
        let start = idx + 4;
        let digits: String = s
            .get(start..)?
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        if !digits.is_empty() {
            return digits.parse().ok();
        }
    }
    None
}

fn get_begin_minus_number(s: &str) -> Option<i32> {
    if s.starts_with('-') {
        let digits: String = s
            .get(1..)?
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        if !digits.is_empty() {
            return digits.parse().ok();
        }
    }
    None
}

fn get_inv_num(s: &str) -> Option<i32> {
    if let Some(idx) = s.find("<inv") {
        let start = idx + 4;
        let digits: String = s
            .get(start..)?
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        if !digits.is_empty() {
            return digits.parse().ok();
        }
    }
    None
}

fn matches_begin_minus_number_caret(s: &str) -> bool {
    if s.starts_with('-') {
        if let Some(rest) = s.get(1..) {
            let digits_count = rest.chars().take_while(|c| c.is_ascii_digit()).count();
            digits_count > 0 && rest.get(digits_count..).is_some_and(|r| r.starts_with('^'))
        } else {
            false
        }
    } else {
        false
    }
}

fn matches_some_sv_numbers(s: &str) -> bool {
    let chars: Vec<char> = s.chars().collect();
    if chars.len() < 6 {
        return false;
    }
    for i in 0..=chars.len() - 6 {
        if chars[i] == '<' {
            let sv_type = &chars[i + 1..i + 4];
            if sv_type.iter().all(|&c| c != '<' && c != '>') {
                let mut j = i + 4;
                let mut has_digits = false;
                while j < chars.len() && chars[j].is_ascii_digit() {
                    has_digits = true;
                    j += 1;
                }
                if has_digits && j < chars.len() && chars[j] == '>' {
                    return true;
                }
            }
        }
    }
    false
}

fn get_amp_atgc(s: &str) -> Option<String> {
    if let Some(idx) = s.find('&') {
        let extra: String = s
            .get(idx + 1..)?
            .chars()
            .take_while(|c| matches!(c, 'A' | 'T' | 'G' | 'C'))
            .collect();
        if !extra.is_empty() {
            return Some(extra);
        }
    }
    None
}

fn get_hash_group_caret_group(s: &str) -> Option<(String, &str)> {
    if let Some(hash_idx) = s.find('#') {
        let rest = s.get(hash_idx + 1..)?;
        if let Some(caret_idx) = rest.find('^') {
            let group1 = rest.get(..caret_idx)?.to_string();
            let group2 = rest.get(caret_idx + 1..)?;
            if !group1.is_empty() && !group2.is_empty() {
                return Some((group1, group2));
            }
        }
    }
    None
}

fn get_begin_digits(s: &str) -> Option<i32> {
    let digits: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
    if !digits.is_empty() {
        digits.parse().ok()
    } else {
        None
    }
}

fn get_caret_atgnc(s: &str) -> Option<String> {
    if let Some(idx) = s.find('^') {
        let extra: String = s
            .get(idx + 1..)?
            .chars()
            .take_while(|c| matches!(c, 'A' | 'T' | 'G' | 'N' | 'C'))
            .collect();
        if !extra.is_empty() {
            return Some(extra);
        }
    }
    None
}

fn matches_any_sv(s: &str) -> bool {
    let chars: Vec<char> = s.chars().collect();
    if chars.len() < 5 {
        return false;
    }
    for i in 0..=chars.len() - 5 {
        if chars[i] == '<' && chars[i + 4] == '>' {
            let inner = &chars[i + 1..i + 4];
            if inner.iter().all(|&c| c != '<' && c != '>') {
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
        vm.sv = Some(Sv {
            r#type: "INV".to_string(),
            pairs: 2,
            splits: 3,
            clusters: 1,
        });
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
            vars_at_cur_position.map.keys().cloned().collect();
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
        assert!(vars10.var_desc_to_variant.contains_key("+T"));

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
        assert!(vars11.var_desc_to_variant.contains_key("G"));
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
}
