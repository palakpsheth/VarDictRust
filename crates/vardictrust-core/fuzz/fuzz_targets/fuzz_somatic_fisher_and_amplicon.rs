#![no_main]

use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet, HashSet};

use libfuzzer_sys::fuzz_target;
use vardictrust_core::config::{Config, COMPAT_OFF, COMPAT_VARDICTJAVA_1_8_3};
use vardictrust_core::javamap::JavaHashMap;
use vardictrust_core::postprocess::amplicon::{
    amplicon_post_process, amplicon_post_process_streaming, amplicon_retain_strict_bounds,
    amplicon_select_vcf_rows_compat,
};
use vardictrust_core::postprocess::simple::simple_post_process;
use vardictrust_core::postprocess::somatic::{somatic_post_process_at, PositionPipelineFn};
use vardictrust_core::regions::{AmpliconInfo, Region};
use vardictrust_core::stats::fisher::FisherExact;
use vardictrust_core::types::{Pos1, VardictError};
use vardictrust_core::variations::{Variant, Vars};

/// Maximum accepted libFuzzer input payload (`4 KiB`) to bound per-iteration runtime.
const MAX_FUZZ_INPUT_BYTES: usize = 4096;
/// Maximum number of distinct genomic positions per `JavaHashMap<i32, Vars>` or `BTreeMap`.
const MAX_POSITIONS: usize = 32;
/// Maximum `Region` and `AmpliconInfo` interval span (`256 bp`, `AR-06`).
const MAX_SPAN_BP: u32 = 256;
/// Floor (`> 16_777_216 = MAX_DENSE_FISHER_SUPPORT`) triggering unimodal `FisherExact`.
const UNIMODAL_FLOOR: i32 = 16_777_220;

struct ByteCursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> ByteCursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    fn u8(&mut self) -> u8 {
        if self.data.is_empty() {
            return 0;
        }
        let idx = self.pos % self.data.len();
        self.pos = self.pos.wrapping_add(1);
        self.data.get(idx).copied().unwrap_or(0)
    }

    fn bool(&mut self) -> bool {
        (self.u8() & 1) != 0
    }

    fn u16(&mut self) -> u16 {
        u16::from_le_bytes([self.u8(), self.u8()])
    }

    fn i32(&mut self) -> i32 {
        i32::from_le_bytes([self.u8(), self.u8(), self.u8(), self.u8()])
    }

    fn f64_bounded(&mut self, max_val: f64) -> f64 {
        match self.u8() & 0x0f {
            0 => 0.0,
            1 => -0.0,
            2 => 0.01,
            3 => 0.05,
            4 => 0.25,
            5 => 0.5,
            6 => 1.0,
            7 => f64::from_bits(0x3f0a_36e2_eb1c_432d), // 5e-5 HALF_EVEN tie
            _ => (f64::from(self.u16()) / 65535.0) * max_val,
        }
    }

    fn dna_string(&mut self, max_len: usize) -> String {
        const BASES: &[u8; 5] = b"ACGTN";
        let len = (usize::from(self.u8()) % max_len.max(1)).saturating_add(1);
        let mut out = String::with_capacity(len);
        for _ in 0..len {
            let b = BASES[usize::from(self.u8()) % BASES.len()] as char;
            out.push(b);
        }
        out
    }
}

fn check_fisher_table(ref_fwd: i32, ref_rev: i32, alt_fwd: i32, alt_rev: i32) {
    debug_assert!(ref_fwd >= 0 && ref_rev >= 0 && alt_fwd >= 0 && alt_rev >= 0);
    let mut or_buf = String::with_capacity(32);
    for java_fisher in [false, true] {
        let fe = FisherExact::new_with_mode(ref_fwd, ref_rev, alt_fwd, alt_rev, java_fisher);
        for p in [
            fe.get_p_value(),
            fe.get_p_value_greater(),
            fe.get_p_value_less(),
        ] {
            assert!(
                !p.is_nan(),
                "FisherExact p-value NaN for ({ref_fwd}, {ref_rev}, {alt_fwd}, {alt_rev}, java_fisher={java_fisher})"
            );
            assert!(
                (0.0..=1.0).contains(&p),
                "FisherExact p-value {p} outside [0, 1] for ({ref_fwd}, {ref_rev}, {alt_fwd}, {alt_rev}, java_fisher={java_fisher})"
            );
        }
        or_buf.clear();
        fe.write_odds_ratio_into(&mut or_buf);
        assert!(
            !or_buf.is_empty(),
            "write_odds_ratio_into empty for ({ref_fwd}, {ref_rev}, {alt_fwd}, {alt_rev}, java_fisher={java_fisher})"
        );
    }
}

fn fuzz_fisher_exact(cur: &mut ByteCursor<'_>) {
    // 1. Standard dense counts in 0..=2_048 (including degenerate zero marginals).
    let (rf, rr, af, ar) = match cur.u8() & 0x0f {
        0 => (0, 0, 0, 0),
        1 => (
            0,
            0,
            i32::from(cur.u16() % 2049),
            i32::from(cur.u16() % 2049),
        ),
        2 => (
            i32::from(cur.u16() % 2049),
            0,
            i32::from(cur.u16() % 2049),
            0,
        ),
        _ => (
            i32::from(cur.u16() % 2049),
            i32::from(cur.u16() % 2049),
            i32::from(cur.u16() % 2049),
            i32::from(cur.u16() % 2049),
        ),
    };
    check_fisher_table(rf, rr, af, ar);

    // 2. 1/64-gated medium dense counts in 0..=8_192.
    let gate = cur.u8() & 0x3f;
    if gate == 0 {
        let m_rf = i32::from(cur.u16() % 8193);
        let m_rr = i32::from(cur.u16() % 8193);
        let m_af = i32::from(cur.u16() % 8193);
        let m_ar = i32::from(cur.u16() % 8193);
        check_fisher_table(m_rf, m_rr, m_af, m_ar);
    } else if gate == 1 {
        // 3. 1/64-gated bounded-sigma unimodal ultra-deep shapes (AR-03).
        let j0 = i32::from(cur.u16() % 512);
        let j1 = i32::from(cur.u16() % 512);
        if cur.bool() {
            let c1 = i32::from(cur.u16() % 2049);
            let c2 = i32::from(cur.u16() % 2049);
            check_fisher_table(UNIMODAL_FLOOR + j0, c1, c2, UNIMODAL_FLOOR + j1);
        } else {
            let a0 = i32::from(cur.u16() % 513);
            let a1 = i32::from(cur.u16() % 513);
            check_fisher_table(1_100_000_000 + j0, 1_100_000_000 + j1, a0, a1);
        }
    }
}

fn synth_variant_at(cur: &mut ByteCursor<'_>, anchor_pos: i32, reconcile_vd_dp: bool) -> Variant {
    let mut v = Variant::new();
    let kind_sel = cur.u8() % 10;
    let shift = i32::from(cur.u8() % 5) - 1;
    v.start_pos = anchor_pos.saturating_add(shift).max(1);

    match kind_sel {
        0 => {
            // Reference-matching row
            let b = b"ACGT"[usize::from(cur.u8()) % 4] as char;
            v.ref_allele = b.to_string();
            v.var_allele = b.to_string();
            v.description = b.to_string();
            v.end_pos = v.start_pos;
        }
        1..=3 => {
            // SNV
            let r_idx = usize::from(cur.u8()) % 4;
            let a_idx = (r_idx + 1 + usize::from(cur.u8() % 3)) % 4;
            let rb = b"ACGT"[r_idx] as char;
            let ab = b"ACGT"[a_idx] as char;
            v.ref_allele = rb.to_string();
            v.var_allele = ab.to_string();
            v.description = ab.to_string();
            v.end_pos = v.start_pos;
        }
        4 => {
            // Deletion matching MINUS_NUM_NUM (-\d\d) to exercise combineAnalysis & splice filter
            let del_len = usize::from(cur.u8() % 14) + 10; // 10..=23
            let deleted = cur.dna_string(del_len);
            v.ref_allele = format!("A{deleted}");
            v.var_allele = "A".to_string();
            v.description = format!("-{}{deleted}", deleted.len());
            v.end_pos = v.start_pos.saturating_add(deleted.len() as i32);
        }
        5 => {
            // Long insertion (nt.len() > 10) to exercise combineAnalysis
            let ins = cur.dna_string(14);
            v.ref_allele = "G".to_string();
            v.var_allele = format!("G{ins}");
            v.description = format!("+{ins}");
            v.end_pos = v.start_pos;
        }
        6 => {
            // Complex MNP with shared prefix/suffix to exercise adj_complex()
            let mid_r = cur.dna_string(6);
            let mid_a = cur.dna_string(6);
            v.ref_allele = format!("AC{mid_r}GT");
            v.var_allele = format!("AC{mid_a}GT");
            v.description = format!("-{}&{}", v.ref_allele.len(), v.var_allele);
            v.end_pos = v.start_pos.saturating_add(v.ref_allele.len() as i32 - 1);
        }
        7 => {
            // Symbolic SV allele
            let sv_tok = ["<DEL>", "<DUP>", "<INV>", "<INS>"][usize::from(cur.u8()) % 4];
            v.ref_allele = "C".to_string();
            v.var_allele = sv_tok.to_string();
            v.description = sv_tok.to_string();
            v.end_pos = v.start_pos.saturating_add(i32::from(cur.u8() % 120) + 1);
        }
        8 => {
            // N-containing reference allele
            v.ref_allele = "AN".to_string();
            v.var_allele = "AT".to_string();
            v.description = "#AT".to_string();
            v.end_pos = v.start_pos.saturating_add(1);
        }
        _ => {
            v.ref_allele = cur.dna_string(16);
            v.var_allele = cur.dna_string(16);
            v.description = format!("^{}", v.var_allele);
            v.end_pos = v
                .start_pos
                .saturating_add(i32::from(cur.u8() % 16))
                .max(v.start_pos);
        }
    }

    let fwd = i32::from(cur.u8() % 64);
    let rev = i32::from(cur.u8() % 64);
    v.fwd = fwd;
    v.rev = rev;
    v.pos_coverage = if cur.bool() {
        // Low coverage window (< min_reads + 3) so combineAnalysis triggers
        i32::from(cur.u8() % 5) + 1
    } else {
        fwd.saturating_add(rev).max(1)
    };
    v.ref_fwd_cov = i32::from(cur.u8() % 80);
    v.ref_rev_cov = i32::from(cur.u8() % 80);
    v.total_pos_cov = v
        .pos_coverage
        .saturating_add(v.ref_fwd_cov)
        .saturating_add(v.ref_rev_cov)
        .max(1);
    v.hicnt = if cur.bool() {
        v.pos_coverage.max(2)
    } else {
        i32::from(cur.u8() % 16)
    };
    v.hicov = v.total_pos_cov.max(v.hicnt);

    // Bias towards passing is_good_var while occasionally injecting noise-level normal rows (#45).
    if (cur.u8() & 7) == 0 {
        // Noise candidate for Variant::is_noise()
        v.mean_quality = f64::from(cur.u8() % 11);
        v.qstd = false;
        v.pos_coverage = i32::from(cur.u8() % 3) + 1;
        v.frequency = 0.02;
    } else {
        v.mean_quality = 24.0 + f64::from(cur.u8() % 18);
        v.qstd = cur.bool();
        v.frequency = cur.f64_bounded(1.0).clamp(0.0, 1.0);
    }
    v.hi_freq = v.frequency;
    v.extra_freq = v.frequency;
    v.mean_position = 6.0 + f64::from(cur.u8() % 30);
    v.pstd = cur.bool();
    v.mean_mapping_quality = f64::from(cur.u8() % 61);
    v.hq_lq_ratio = 1.6 + f64::from(cur.u8() % 20);
    v.msi = f64::from(cur.u8() % 18);
    v.msint = i32::from(cur.u8() % 5);
    v.shift3 = i32::from(cur.u8() % 8);
    v.nm = f64::from(cur.u8() % 6);
    v.bias = ["0;0", "2;2", "2;1", "0;2", "1;2"][usize::from(cur.u8()) % 5].to_string();
    v.left_seq = cur.dna_string(8);
    v.right_seq = cur.dna_string(8);
    v.crispr = i32::from(cur.u8() as i8);
    v.var_type = v.var_type();
    if reconcile_vd_dp {
        v.pos_coverage = v.pos_coverage.max(0);
        v.ref_fwd_cov = v.ref_fwd_cov.max(0);
        v.ref_rev_cov = v.ref_rev_cov.max(0);
        v.total_pos_cov = v
            .total_pos_cov
            .max(v.pos_coverage)
            .max(v.ref_fwd_cov.saturating_add(v.ref_rev_cov))
            .max(0);
        v.hicnt = v.hicnt.max(0);
        v.hicov = v.hicov.max(v.hicnt).max(0);
        v.frequency = if v.total_pos_cov > 0 {
            (v.pos_coverage as f64 / v.total_pos_cov as f64).clamp(0.0, 1.0)
        } else {
            0.0
        };
        v.hi_freq = if v.hicov > 0 {
            (v.hicnt as f64 / v.hicov as f64).clamp(0.0, 1.0)
        } else {
            0.0
        };
        v.extra_freq = if v.extra_freq.is_finite() {
            v.extra_freq.clamp(0.0, 1.0)
        } else {
            0.0
        };
    }
    v
}

fn synth_vars_at(cur: &mut ByteCursor<'_>, pos: i32, reconcile_vd_dp: bool) -> Vars {
    let mut vars = Vars::new();
    if cur.bool() {
        let mut rv = synth_variant_at(cur, pos, reconcile_vd_dp);
        rv.ref_allele = "A".to_string();
        rv.var_allele = "A".to_string();
        rv.description = "A".to_string();
        rv.var_type.clear();
        vars.ref_variant = Some(rv);
    }
    if (cur.u8() & 7) == 0 {
        vars.sv = format!("{}-{}-{}", cur.u8() % 8, cur.u8() % 32 + 1, cur.u8() % 4);
    }
    let n_alts = usize::from(cur.u8() % 4); // 0..=3 variants per position
    let populate_desc_map = cur.bool();
    for _ in 0..n_alts {
        let v = synth_variant_at(cur, pos, reconcile_vd_dp);
        if populate_desc_map {
            vars.var_desc_to_variant
                .insert(v.description.clone(), v.clone());
        }
        vars.variants.push(v);
    }
    vars
}

fn synth_vars_map(
    cur: &mut ByteCursor<'_>,
    reg_start: u32,
    reg_span: u32,
    reconcile_vd_dp: bool,
) -> JavaHashMap<i32, Vars> {
    let mut map = JavaHashMap::new();
    let count = (usize::from(cur.u8()) % MAX_POSITIONS).saturating_add(1);
    let span_with_halo = reg_span.saturating_add(8).max(1);
    let origin = reg_start.saturating_sub(4).max(1);
    for _ in 0..count {
        if map.len() >= MAX_POSITIONS {
            break;
        }
        let off = u32::from(cur.u16()) % span_with_halo;
        let pos = origin.saturating_add(off).min(i32::MAX as u32) as i32;
        let entry = synth_vars_at(cur, pos, reconcile_vd_dp);
        map.insert(pos, entry);
    }
    map
}

fuzz_target!(|data: &[u8]| {
    if data.len() < 16 || data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    let mut cur = ByteCursor::new(data);

    // Part A: FisherExact R vs Java (--fisher) statistical kernel.
    fuzz_fisher_exact(&mut cur);

    // Bounded genomic interval (<= 256 bp, Pos1 >= 1, AR-06).
    let mode = cur.u8();
    let reg_start = if (mode & 1) != 0 {
        let top = (i32::MAX as u32).saturating_sub(MAX_SPAN_BP + 16);
        top.saturating_sub(u32::from(cur.u8())).max(1)
    } else {
        (u32::from(cur.u16()) % 50_000).saturating_add(1)
    };
    let reg_span = (u32::from(cur.u8()) % MAX_SPAN_BP).saturating_add(1); // 1..=256
    let reg_end = reg_start
        .saturating_add(reg_span - 1)
        .min(i32::MAX as u32)
        .max(reg_start);

    let region = Region {
        chr: "chr1".to_string(),
        start: Pos1(reg_start),
        end: Pos1(reg_end),
        gene: "FUZZ".to_string(),
        amplicon: None,
    };

    let mut base_cfg = Config::default();
    base_cfg.pileup = (mode & 2) != 0;
    base_cfg.emit_intermediate = (mode & 4) != 0;
    base_cfg.vcf_pass_only = (mode & 8) != 0;
    base_cfg.debug = (mode & 16) != 0;
    base_cfg.fisher = (mode & 32) != 0;
    base_cfg.crispr_cutting_site = match cur.u8() % 4 {
        0 => None,
        1 => Some(0),
        2 => Some(reg_start.saturating_add(reg_span / 2) as i32),
        _ => Some(cur.i32()),
    };
    base_cfg.crispr_flanking_size = match cur.u8() % 3 {
        0 => None,
        1 => Some(0),
        _ => Some(i32::from(cur.u8() as i8)),
    };
    base_cfg.freq = if base_cfg.pileup && (mode & 64) != 0 {
        -1.0
    } else {
        cur.f64_bounded(0.25)
    };
    base_cfg.min_reads = u32::from(cur.u8() % 5);
    base_cfg.lofreq = cur.f64_bounded(0.2).max(0.001);
    base_cfg.vcf_all_alleles = cur.bool();
    base_cfg.compat = if cur.bool() {
        COMPAT_OFF
    } else {
        COMPAT_VARDICTJAVA_1_8_3
    };

    let reconcile_vd_dp = base_cfg.compat.reconcile_vd_dp;
    let aligned_bam1 = synth_vars_map(&mut cur, reg_start, reg_span, reconcile_vd_dp);
    let aligned_bam2 = synth_vars_map(&mut cur, reg_start, reg_span, reconcile_vd_dp);

    // Populate optional splice junction strings matching synthesized deletions.
    let mut splice = HashSet::default();
    if (mode & 128) != 0 {
        for (_, vars) in aligned_bam1.iter().take(4) {
            if let Some(v) = vars.variants.first() {
                splice.insert(format!("{}-{}", v.start_pos, v.end_pos));
            }
        }
    }

    // Part B: simple_post_process (single-sample germline / pileup / CRISPR / VCF bucketing).
    let simple_rows = simple_post_process(&base_cfg, &region, &aligned_bam1, &splice);
    if base_cfg.crispr_cutting_site.unwrap_or(0) == 0 {
        for row in &simple_rows {
            assert_eq!(row.crispr, 0);
        }
    }
    if !base_cfg.emit_intermediate && !base_cfg.vcf_all_alleles {
        for pair in simple_rows.windows(2) {
            assert!(pair[0].start_pos < pair[1].start_pos);
        }
    } else if !base_cfg.emit_intermediate && base_cfg.vcf_all_alleles {
        let mut inter_cfg = base_cfg.clone();
        inter_cfg.emit_intermediate = true;
        let inter_rows = simple_post_process(&inter_cfg, &region, &aligned_bam1, &splice);
        assert_eq!(simple_rows.len(), inter_rows.len());
    }
    if base_cfg.compat.strict_region_bounds {
        for row in &simple_rows {
            assert!(row.start_pos >= region.start.0 as i32 && row.start_pos <= region.end.0 as i32);
        }
    }

    // Part C: somatic_post_process_at with mock PositionPipelineFn returning Ok(Some), Ok(None), or Err.
    let max_rl = i32::from(cur.u8() % 150) + 25;
    let use_pipeline = cur.bool();
    let pipe_seed = cur.u8();
    let call_counter = Cell::new(0u8);
    let comb_template = synth_vars_at(&mut cur, reg_start as i32, reconcile_vd_dp);
    let mock_pipe = |comb_reg: &Region,
                     splice_set: &mut HashSet<String, _>,
                     rl: &mut i32,
                     target_pos: i32|
     -> Result<Option<Vars>, VardictError> {
        let step = call_counter.get();
        call_counter.set(step.wrapping_add(1));
        let selector = pipe_seed.wrapping_add(step).wrapping_add(target_pos as u8);
        match selector % 4 {
            0 => Err(VardictError::Parse("mock combine error".to_string())),
            1 => Ok(None),
            _ => {
                *rl = (*rl).max(50);
                if (selector & 4) != 0 {
                    splice_set.insert(format!("{}-{}", comb_reg.start.0, comb_reg.end.0));
                }
                let mut out = comb_template.clone();
                if let Some(v1_entry) = aligned_bam1.get(&target_pos) {
                    if let Some(first_v1) = v1_entry.variants.first() {
                        if (selector & 8) != 0 && !out.variants.is_empty() {
                            out.variants[0].description = first_v1.description.clone();
                            out.variants[0].start_pos = first_v1.start_pos;
                            out.variants[0].end_pos = first_v1.end_pos;
                        }
                    }
                }
                Ok(Some(out))
            }
        }
    };
    let pipe_ref: Option<&PositionPipelineFn<'_>> =
        if use_pipeline { Some(&mock_pipe) } else { None };
    let somatic_rows = somatic_post_process_at(
        &base_cfg,
        &region,
        &aligned_bam1,
        &aligned_bam2,
        &splice,
        pipe_ref,
        max_rl,
    );
    for sv in &somatic_rows {
        assert!(!sv.status.is_empty());
        if reconcile_vd_dp {
            for v in std::iter::once(&sv.variant1).chain(sv.variant2.as_ref()) {
                assert!(
                    0 <= v.pos_coverage && v.pos_coverage <= v.total_pos_cov,
                    "somatic reconcile_vd_dp violated 0 <= VD ({}) <= DP ({}) on status={}",
                    v.pos_coverage,
                    v.total_pos_cov,
                    sv.status
                );
                assert!(
                    v.frequency.is_finite() && (0.0..=1.0).contains(&v.frequency),
                    "somatic reconcile_vd_dp violated AF in [0, 1]: {}",
                    v.frequency
                );
                assert!(
                    v.hi_freq.is_finite() && (0.0..=1.0).contains(&v.hi_freq),
                    "somatic reconcile_vd_dp violated HIAF in [0, 1]: {}",
                    v.hi_freq
                );
                assert!(
                    v.extra_freq.is_finite() && (0.0..=1.0).contains(&v.extra_freq),
                    "somatic reconcile_vd_dp violated ADJAF in [0, 1]: {}",
                    v.extra_freq
                );
                assert!(
                    0 <= v.hicnt && v.hicnt <= v.hicov,
                    "somatic reconcile_vd_dp violated 0 <= HICNT ({}) <= HICOV ({})",
                    v.hicnt,
                    v.hicov
                );
            }
        }
    }

    // Part D: amplicon_post_process + amplicon_post_process_streaming + amplicon_select_vcf_rows
    // over 1..=6 overlapping Regions clamped to <= 256 bp and sorted BTreeMap<i32, Vec<(usize, &Region)>>.
    let mut amp_cfg = base_cfg.clone();
    amp_cfg.amplicon = Some("10:0.95".to_string());

    let num_amps = (usize::from(cur.u8()) % 6).saturating_add(1); // 1..=6
    let mut amp_regions = Vec::with_capacity(num_amps);
    let mut vars_per_amplicon = Vec::with_capacity(num_amps);
    let mut candidate_positions = BTreeSet::new();

    for idx in 0..num_amps {
        let rel_off = u32::from(cur.u8()) % reg_span;
        let a_start = reg_start.saturating_add(rel_off).min(reg_end);
        let a_span = (u32::from(cur.u8()) % MAX_SPAN_BP).saturating_add(1);
        let a_end = a_start
            .saturating_add(a_span - 1)
            .min(reg_start.saturating_add(MAX_SPAN_BP - 1))
            .min(i32::MAX as u32)
            .max(a_start);
        let ins_trim_l = u32::from(cur.u8() % 8);
        let ins_trim_r = u32::from(cur.u8() % 8);
        let ins_start = a_start.saturating_add(ins_trim_l).min(a_end);
        let ins_end = a_end.saturating_sub(ins_trim_r).max(ins_start);

        let amp_info = if idx > 0 && (cur.u8() & 0x0f) == 0 {
            None
        } else {
            Some(AmpliconInfo {
                insert_start: Pos1(ins_start),
                insert_end: Pos1(ins_end),
            })
        };
        amp_regions.push(Region {
            chr: "chr1".to_string(),
            start: Pos1(a_start),
            end: Pos1(a_end),
            gene: format!("AMP{idx}"),
            amplicon: amp_info,
        });

        let amp_span = a_end.saturating_sub(a_start).saturating_add(1);
        let amp_map = synth_vars_map(&mut cur, a_start, amp_span, reconcile_vd_dp);
        for &p in amp_map.keys() {
            if candidate_positions.len() < MAX_POSITIONS {
                candidate_positions.insert(p);
            }
        }
        if candidate_positions.len() < MAX_POSITIONS {
            candidate_positions.insert(ins_start as i32);
        }
        vars_per_amplicon.push(amp_map);
    }

    let mut amplicons_on_positions: BTreeMap<i32, Vec<(usize, &Region)>> = BTreeMap::new();
    for pos in candidate_positions {
        let mut covering = Vec::new();
        for (idx, amp_reg) in amp_regions.iter().enumerate() {
            if let Some(ref amp) = amp_reg.amplicon {
                if pos >= amp.insert_start.0 as i32 && pos <= amp.insert_end.0 as i32 {
                    covering.push((idx, amp_reg));
                }
            }
        }
        if !covering.is_empty() {
            amplicons_on_positions.insert(pos, covering);
        }
    }

    let batch_amp = amplicon_post_process(
        &amp_cfg,
        &region,
        &vars_per_amplicon,
        &amplicons_on_positions,
        &splice,
    )
    .expect("infallible collector in amplicon_post_process");

    let mut streamed_amp = Vec::new();
    amplicon_post_process_streaming(
        &amp_cfg,
        &region,
        &vars_per_amplicon,
        amplicons_on_positions
            .iter()
            .map(|(&pos, regs)| (pos, regs.as_slice())),
        &splice,
        &mut |av| {
            streamed_amp.push(av);
            Ok(())
        },
    )
    .expect("infallible collector in amplicon_post_process_streaming");

    if amp_cfg.compat.strict_region_bounds {
        amplicon_retain_strict_bounds(&mut streamed_amp, &region, &amplicons_on_positions);
    }
    let finalized_streamed = if !amp_cfg.emit_intermediate && !amp_cfg.vcf_all_alleles {
        amplicon_select_vcf_rows_compat(streamed_amp, amp_cfg.compat.exact_binary_half_even)
    } else {
        streamed_amp
    };
    assert_eq!(batch_amp.len(), finalized_streamed.len());
    for (lhs, rhs) in batch_amp.iter().zip(finalized_streamed.iter()) {
        assert_eq!(lhs.variant.start_pos, rhs.variant.start_pos);
        assert_eq!(lhs.variant.var_allele, rhs.variant.var_allele);
        assert_eq!(lhs.gdamp, rhs.gdamp);
        assert_eq!(lhs.tlamp, rhs.tlamp);
        assert_eq!(lhs.ncamp, rhs.ncamp);
        assert_eq!(lhs.ampflag, rhs.ampflag);
    }
});
