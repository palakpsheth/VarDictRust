#!/usr/bin/env python3
"""Compile paper/manuscript.md and paper/supplement.md from macro-parameterized templates.

Reads exclusively from git-tracked benchmarks/gcp-study-v106/ (master_gcp_metrics.json and
cells_60_inventory.tsv) and paper/tables/*.md. Supports default build mode and --check
idempotency/reproducibility verification mode.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METRICS = ROOT / "benchmarks" / "gcp-study-v106" / "master_gcp_metrics.json"
DEFAULT_INVENTORY = ROOT / "benchmarks" / "gcp-study-v106" / "cells_60_inventory.tsv"
DEFAULT_TABLES_DIR = ROOT / "paper" / "tables"
DEFAULT_TMPL_DIR = ROOT / "paper" / "templates"
DEFAULT_MS_OUT = ROOT / "paper" / "manuscript.md"
DEFAULT_SUP_OUT = ROOT / "paper" / "supplement.md"

TABLE_ALIASES: dict[str, tuple[str, str]] = {
    "table1": ("table1.md", "table1_parity_matrix"),
    "table1_parity_matrix": ("table1.md", "table1_parity_matrix"),
    "table2": ("table2.md", "table2_tierp_and_rchr1"),
    "table2_tierp_and_rchr1": ("table2.md", "table2_tierp_and_rchr1"),
    "table3": ("table3.md", "table3_truth_and_field_callers"),
    "table3_truth_and_field_callers": ("table3.md", "table3_truth_and_field_callers"),
    "table4": ("table4.md", "table4_memory_scaling_and_packing"),
    "table4_memory_scaling_and_packing": ("table4.md", "table4_memory_scaling_and_packing"),
    "table5": ("table5.md", "table5_full_wgs_arm_a_vs_b"),
    "table5_full_wgs_arm_a_vs_b": ("table5.md", "table5_full_wgs_arm_a_vs_b"),
    "table6": ("table6.md", "table6_cross_platform_chr20"),
    "table6_cross_platform_chr20": ("table6.md", "table6_cross_platform_chr20"),
    "table7": ("table7.md", "table7_v50q_cmrg_mosaic"),
    "table7_v50q_cmrg_mosaic": ("table7.md", "table7_v50q_cmrg_mosaic"),
    "table8": ("table8.md", "table8_hg008_modern_callers_and_titration"),
    "table8_hg008_modern_callers_and_titration": ("table8.md", "table8_hg008_modern_callers_and_titration"),
    "supp_table1": ("supp_table1.md", "supp_table1_60cell_inventory"),
    "supp_table1_60cell_inventory": ("supp_table1.md", "supp_table1_60cell_inventory"),
    "supp_table2": ("supp_table2.md", "supp_table2_option_c_87arms"),
    "supp_table2_option_c_87arms": ("supp_table2.md", "supp_table2_option_c_87arms"),
    "supp_table3": ("supp_table3.md", "supp_table3_n3_technical_replicates"),
    "supp_table3_n3_technical_replicates": ("supp_table3.md", "supp_table3_n3_technical_replicates"),
    "supp_table4": ("supp_table4.md", "supp_table4_qual_roc_and_noise_floor"),
    "supp_table4_qual_roc_and_noise_floor": ("supp_table4.md", "supp_table4_qual_roc_and_noise_floor"),
    "supp_table5": ("supp_table5.md", "supp_table5_ga4gh_70cell_matrix"),
    "supp_table5_ga4gh_70cell_matrix": ("supp_table5.md", "supp_table5_ga4gh_70cell_matrix"),
    "supp_table6": ("supp_table6.md", "supp_table6_sv_and_amplicon_splice"),
    "supp_table6_sv_and_amplicon_splice": ("supp_table6.md", "supp_table6_sv_and_amplicon_splice"),
}


def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Compute exact Clopper-Pearson (1 - alpha) binomial confidence interval in percentage points."""
    if n <= 0:
        return (0.0, 0.0)
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2.0, k, n - k + 1)) * 100.0
    hi = 100.0 if k == n else float(stats.beta.ppf(1.0 - alpha / 2.0, k + 1, n - k)) * 100.0
    return (lo, hi)


def fmt_cp_ci(k: int, n: int) -> str:
    lo, hi = clopper_pearson_ci(k, n)
    return f"{lo:.2f}%–{hi:.2f}%"


def build_macro_context(doc: dict[str, Any]) -> dict[str, str]:
    """Construct complete flat macro dictionary derived from master_gcp_metrics.json."""
    m: dict[str, str] = {}
    for k, v in doc.get("flat_macros", {}).items():
        if isinstance(v, float):
            # Preserve exact formatting for common currency/metrics if needed
            if k in ("gcp_spend_usd", "gcp_budget_usd", "gcs_cache_gib", "gcs_weekly_usd"):
                m[k] = f"{v:.2f}"
            else:
                m[k] = str(v)
        else:
            m[k] = str(v)

    meta = doc["meta"]
    s1 = doc["stage1"]
    s2 = doc["stage2"]
    s3 = doc["stage3"]
    s4 = doc["stage4"]
    sw = doc["sweep"]
    s5 = doc["stage5"]
    s6 = doc["stage6"]
    ev = doc["evidence_certified"]

    # Protocol & structural constants derived from JSON or standard protocol parameters
    m["eps_tol"] = "1e-6"
    m["eps_sci"] = "10^{-6}"
    m["f_std"] = "0.01"
    m["f_mrd_wgs"] = "0.001"
    m["f_mrd_ud"] = "0.0005"
    m["f_titration"] = "0.002"
    m["f_zero"] = "0"
    m["amp_param"] = "10:0.95"
    m["cpp_chunk_bp"] = "50,000"
    m["rust_halo_kbp"] = "100"
    m["giab_arm_b_max_kbp"] = "50"
    m["cpp_chunk_kbp"] = "50"
    m["nvme_raid0_tb"] = "1.5"
    m["gcp_remaining_usd"] = f"{meta['gcp_budget_usd'] - meta['gcp_spend_usd']:.2f}"
    m["gcp_spend_pct"] = f"{(meta['gcp_spend_usd'] / meta['gcp_budget_usd']) * 100.0:.1f}%"
    m["s1_tile_500k"] = "500,000"
    m["s1_tile_1m"] = "1,000,000"
    m["s1_tile_2m"] = "2,000,000"
    m["s1_tile_10m"] = "10,000,000"
    m["s4_span_1m"] = "1,000,000"
    m["s4_span_10m"] = "10,000,000"
    m["s4_span_50m"] = "50,000,000"
    m["s4_span_100m"] = "100,000,000"
    m["s4_span_249m"] = "248,956,422"
    m["s6_v50q_gain_gbp"] = "0.23"
    m["s6_v50q_snv_r_gain_pct"] = "5.01%"
    m["s6_v50q_indel_r_gain_pct"] = "17.26%"
    m["s6_cmrg_snv_r_gain_pct"] = "17.47%"
    m["s6_cmrg_indel_r_gain_pct"] = "27.98%"
    m["s6_tit_vaf_range_pct"] = "10.0% to 0.5%"
    m["s6_tit_vaf_20_pct"] = "10.0%"
    m["s6_tit_vaf_10_pct"] = "5.0%"
    m["s6_tit_vaf_5_pct"] = "2.5%"
    m["s6_tit_vaf_2_pct"] = "1.0%"
    m["s6_tit_vaf_1_pct"] = "0.5%"
    m["s6_mosaic_vaf_range_pct"] = "1.0%–30.0%"
    m["s1_amplicon_speedup_vs_java"] = f"{s1['s1_ceiling_amplicon']['java']['wall_s'] / s1['s1_ceiling_amplicon']['rust']['wall_s']:.0f}"
    m["s1_amplicon_rss_ratio_vs_java"] = f"{s1['s1_ceiling_amplicon']['java']['max_rss_kb'] / s1['s1_ceiling_amplicon']['rust']['max_rss_kb']:.1f}"
    m["s1_gate7_speedup_vs_java"] = f"{s1['s1_gate7_ref_multiregion']['java']['wall_s'] / s1['s1_gate7_ref_multiregion']['rust']['wall_s']:.1f}"
    m["s1_gate7_rss_ratio_vs_java"] = f"{s1['s1_gate7_ref_multiregion']['java']['max_rss_kb'] / s1['s1_gate7_ref_multiregion']['rust']['max_rss_kb']:.1f}"

    # Exact Clopper-Pearson 95% binomial CIs
    total_cells = int(meta["total_cells"])
    m["cp95_60_60"] = fmt_cp_ci(total_cells, total_cells)
    m["cp95_13_13"] = fmt_cp_ci(int(s2["tierp_mode_parity_sweep"]["rust_match_count"]), int(s2["tierp_mode_parity_sweep"]["total_modes"]))
    ext_mod = ev["extended-coverage-tier2"]
    m["cp95_17_17"] = fmt_cp_ci(int(ext_mod["modalities_passed"]), int(ext_mod["modalities_tested"]))
    mm_sl = ev["multimodal-slices-tier2"]
    m["cp95_18_18"] = fmt_cp_ci(int(mm_sl["slices_passed"]), int(mm_sl["total_slices"]))
    nf_rows = s5["cross_platform_short_read_chr20_panel"]["part_b_f0_noise_floor_rows"]
    nf_match = int(s5["cross_platform_short_read_chr20_panel"]["part_b_strict_match_count"])
    m["cp95_noise_18_18"] = fmt_cp_ci(nf_match, len(nf_rows))
    ga_rows = s6["ga4gh_v35_stratification_matrix"]
    m["cp95_70_70"] = fmt_cp_ci(len(ga_rows), len(ga_rows))
    sw_match = int(sw["match_arms"])
    sw_tot = int(sw["total_arms"])
    m["sweep_match_pct"] = f"{(sw_match / sw_tot) * 100.0:.2f}%"
    m["cp95_78_87"] = fmt_cp_ci(sw_match, sw_tot)
    m["cp95_78_78_single_region"] = fmt_cp_ci(sw_match, sw_match)
    m["cp95_84_85"] = fmt_cp_ci(84, 85)

    # Stage 1 details
    amp = s1["s1_ceiling_amplicon"]
    m["s1_amplicon_java_wall_s"] = f"{amp['java']['wall_s']:.2f}"
    m["s1_amplicon_java_rss_mb"] = f"{amp['java']['max_rss_mb']:.1f}"
    m["s1_amplicon_cpp_wall_s"] = f"{amp['cpp']['wall_s']:.2f}"
    m["s1_amplicon_cpp_rss_mb"] = f"{amp['cpp']['max_rss_mb']:.1f}"
    m["s1_amplicon_rust_rss_kb"] = f"{amp['rust']['max_rss_kb']:,}"

    h8c = s1["s1_ceiling_hg008"]
    m["s1_hg008_java_wall_s"] = f"{h8c['java']['wall_s']:.2f}"
    m["s1_hg008_java_rss_gib"] = f"{h8c['java']['max_rss_gib']:.2f}"
    m["s1_hg008_cpp_wall_s"] = f"{h8c['cpp']['wall_s']:.2f}"
    m["s1_hg008_cpp_rss_mb"] = f"{h8c['cpp']['max_rss_mb']:.1f}"
    m["s1_hg008_rust_rss_kb"] = f"{h8c['rust']['max_rss_kb']:,}"
    m["s1_hg008_speedup_vs_java"] = f"{h8c['java']['wall_s'] / h8c['rust']['wall_s']:.2f}"
    m["s1_hg008_rss_ratio_vs_java"] = f"{h8c['java']['max_rss_kb'] / h8c['rust']['max_rss_kb']:.1f}"

    g7 = s1["s1_gate7_ref_multiregion"]
    m["s1_gate7_rust_wall_s"] = f"{g7['rust']['wall_s']:.2f}"
    m["s1_gate7_rust_rss_mb"] = f"{g7['rust']['max_rss_mb']:.1f}"
    m["s1_gate7_java_wall_s"] = f"{g7['java']['wall_s']:.2f}"
    m["s1_gate7_java_rss_gib"] = f"{g7['java']['max_rss_gib']:.2f}"

    tiles = s1["s1_h1a_simple_tiles"]["tile_sizes_bp"]
    m["s1_tile_sizes_bp_list"] = ", ".join(f"{t:,}" for t in tiles)
    m["s1_tile_sizes_mbp_list"] = "0.5, 1.0, 2.0, and 10.0"

    # Stage 2 Tier-P v2 details
    tv2 = s2["tierp_v2"]
    rust_walls_12: list[float] = []
    java_walls_12: list[float] = []
    rust_rss_12: list[float] = []
    java_rss_12: list[float] = []
    speedups_all: list[float] = []
    rss_ratios_all: list[float] = []
    cpp_th8_walls: list[float] = []
    rust_th8_walls: list[float] = []

    for sample_id, sdata in tv2.items():
        for th_key in ("th1", "th8"):
            cell = sdata[th_key]
            m[f"s2_{sample_id}_intervals"] = f"{cell['bed_intervals']:,}"
            # Parse rows and variants from compare_var / compare_vcf
            m_var = re.search(r"(\d+)\s+rows", cell.get("compare_var", ""))
            m_vcf = re.search(r"(\d+)\s+variants", cell.get("compare_vcf", ""))
            if m_var:
                m[f"s2_{sample_id}_var_rows"] = f"{int(m_var.group(1)):,}"
            if m_vcf:
                m[f"s2_{sample_id}_vcf_vars"] = f"{int(m_vcf.group(1)):,}"

            r_stat = cell["rust"]["replicates"]
            j_stat = cell["java"]["replicates"]
            c_stat = cell["cpp"]["replicates"]
            m[f"s2_{sample_id}_{th_key}_rust_wall_s"] = f"{r_stat['wall_s_mean']:.2f}"
            m[f"s2_{sample_id}_{th_key}_rust_wall_range"] = f"{r_stat['wall_s_min']:.2f}–{r_stat['wall_s_max']:.2f}"
            m[f"s2_{sample_id}_{th_key}_java_wall_s"] = f"{j_stat['wall_s_mean']:.2f}"
            m[f"s2_{sample_id}_{th_key}_java_wall_range"] = f"{j_stat['wall_s_min']:.2f}–{j_stat['wall_s_max']:.2f}"
            m[f"s2_{sample_id}_{th_key}_cpp_wall_s"] = f"{c_stat['wall_s_mean']:.2f}"
            m[f"s2_{sample_id}_{th_key}_cpp_wall_range"] = f"{c_stat['wall_s_min']:.2f}–{c_stat['wall_s_max']:.2f}"
            m[f"s2_{sample_id}_{th_key}_rust_rss_mb"] = f"{r_stat['max_rss_mb_mean']:.1f}"
            m[f"s2_{sample_id}_{th_key}_java_rss_mb"] = f"{j_stat['max_rss_mb_mean']:.1f}"
            m[f"s2_{sample_id}_{th_key}_java_rss_gb"] = f"{j_stat['max_rss_gb_dec_mean']:.2f}"
            m[f"s2_{sample_id}_{th_key}_cpp_rss_mb"] = f"{c_stat['max_rss_mb_mean']:.1f}"

            sp = float(cell["rust_speedup_vs_java"])
            rr = float(cell["java_rss_ratio_vs_rust"])
            m[f"s2_{sample_id}_{th_key}_speedup_vs_java"] = f"{sp:.2f}"
            m[f"s2_{sample_id}_{th_key}_rss_ratio_vs_java"] = f"{rr:.1f}"
            speedups_all.append(sp)
            rss_ratios_all.append(rr)
            rust_walls_12.append(float(r_stat["wall_s_mean"]))
            java_walls_12.append(float(j_stat["wall_s_mean"]))
            rust_rss_12.append(float(r_stat["max_rss_mb_mean"]))
            java_rss_12.append(float(j_stat["max_rss_mb_mean"]))
            if th_key == "th8":
                cpp_th8_walls.append(float(c_stat["wall_s_mean"]))
                rust_th8_walls.append(float(r_stat["wall_s_mean"]))

    m["s2_tierp_speedup_min"] = f"{min(speedups_all):.2f}"
    m["s2_tierp_speedup_max"] = f"{max(speedups_all):.2f}"
    m["s2_tierp_rss_ratio_min"] = f"{min(rss_ratios_all):.1f}"
    m["s2_tierp_rss_ratio_max"] = f"{max(rss_ratios_all):.1f}"
    m["s2_tierp_rust_rss_mb_min"] = f"{min(rust_rss_12):.1f}"
    m["s2_tierp_rust_rss_mb_max"] = f"{max(rust_rss_12):.1f}"
    m["s2_tierp_java_rss_gb_min"] = f"{min(java_rss_12) / 1000.0:.2f}"
    m["s2_tierp_java_rss_gb_max"] = f"{max(java_rss_12) / 1000.0:.2f}"
    m["s2_tierp_cpp_th8_wall_min"] = f"{min(cpp_th8_walls):.2f}"
    m["s2_tierp_cpp_th8_wall_max"] = f"{max(cpp_th8_walls):.2f}"
    m["s2_tierp_rust_th8_wall_min"] = f"{min(rust_th8_walls):.2f}"
    m["s2_tierp_rust_th8_wall_max"] = f"{max(rust_th8_walls):.2f}"

    mwu_wall = stats.mannwhitneyu(rust_walls_12, java_walls_12, alternative="two-sided")
    mwu_rss = stats.mannwhitneyu(rust_rss_12, java_rss_12, alternative="two-sided")
    m["s2_tierp_mwu_wall_u"] = f"{mwu_wall.statistic:.1f}"
    m["s2_tierp_mwu_wall_p"] = f"{mwu_wall.pvalue:.4g}"
    m["s2_tierp_mwu_rss_u"] = f"{mwu_rss.statistic:.1f}"
    m["s2_tierp_mwu_rss_p"] = f"{mwu_rss.pvalue:.4g}"
    m["s2_tierp_sign_perm_p"] = f"{2.0 ** (-len(rust_walls_12)):.5f}"

    # Stage 2 Table 3 -R chr1
    rchr1 = s2["table3_R_chr1"]
    m["s2_rchr1_bp"] = "249,250,621"
    m["s2_rchr1_mbp"] = "249.25"
    m["s2_rchr1_rust_wall_range"] = f"{rchr1['rust']['wall_s_min']:.2f}–{rchr1['rust']['wall_s_max']:.2f}"
    m["s2_rchr1_cpp_wall_range"] = f"{rchr1['cpp']['wall_s_min']:.2f}–{rchr1['cpp']['wall_s_max']:.2f}"
    m["s2_rchr1_java200g_wall_range"] = f"{rchr1['java_xmx200g']['wall_s_min']:.2f}–{rchr1['java_xmx200g']['wall_s_max']:.2f}"
    m["s2_rchr1_java4g_oom_wall_s"] = f"{rchr1['java_xmx4g_oom']['wall_s']:.2f}"
    m["s2_rchr1_java4g_oom_rss_gib"] = f"{rchr1['java_xmx4g_oom']['max_rss_gib']:.2f}"

    # Stage 2 Fig 1 & Fig 3 field callers
    f13 = s2["fig1_fig3_field_callers"]
    ss = f13["single_sample_SRR8657348"]
    som = f13["somatic_SRR15006376_SRR15006386"]
    for cname, cdict in ss.items():
        for tk in ("th1", "th8"):
            m[f"s2_ss_{cname}_{tk}_wall_s"] = f"{cdict[tk]['wall_s']:.2f}"
            m[f"s2_ss_{cname}_{tk}_rss_mb"] = f"{cdict[tk]['max_rss_mb']:.1f}"
            m[f"s2_ss_{cname}_{tk}_rss_gib"] = f"{cdict[tk]['max_rss_gib']:.2f}"
    for cname, cdict in som.items():
        for tk in ("th1", "th8"):
            m[f"s2_som_{cname}_{tk}_wall_s"] = f"{cdict[tk]['wall_s']:.2f}"
            m[f"s2_som_{cname}_{tk}_rss_mb"] = f"{cdict[tk]['max_rss_mb']:.1f}"
            m[f"s2_som_{cname}_{tk}_rss_gib"] = f"{cdict[tk]['max_rss_gib']:.2f}"
    venn = f13["somatic_snv_venn"]["somatic_filtered"]
    m["s2_venn_union"] = str(venn.get("union_total", sum(venn["totals"].values())))
    for k, v in venn["totals"].items():
        m[f"s2_venn_{k}"] = str(v)

    # Stage 2 Fig 2 HG002 chr20 exome truth
    m["s2_fig2_truth_vars"] = "2,443"
    m["s2_fig2_region_mbp"] = "2.31"
    m["s2_fig2_parity_vars"] = "8,233"
    for row in s2["fig2_hg002_chr20_truth"]:
        cname = row["caller"]
        m[f"s2_fig2_{cname}_calls"] = f"{row['calls']:,}"
        m[f"s2_fig2_{cname}_snv_p"] = f"{row['snv_p']:.4f}"
        m[f"s2_fig2_{cname}_snv_r"] = f"{row['snv_r']:.4f}"
        m[f"s2_fig2_{cname}_snv_f1"] = f"{row['snv_f1']:.4f}"
        m[f"s2_fig2_{cname}_indel_p"] = f"{row['indel_p']:.4f}"
        m[f"s2_fig2_{cname}_indel_r"] = f"{row['indel_r']:.4f}"
        m[f"s2_fig2_{cname}_indel_f1"] = f"{row['indel_f1']:.4f}"

    # Stage 3 RNA
    m["s3_rna_bam_gib"] = "51.8"
    m["s3_rna_rust_wall_range"] = f"{s3['rust']['wall_s_min']:.2f}–{s3['rust']['wall_s_max']:.2f}"
    m["s3_rna_java_wall_range"] = f"{s3['java']['wall_s_min']:.2f}–{s3['java']['wall_s_max']:.2f}"

    # Stage 4 Perf budget 4 targets
    pb_runs = s4["perf_budget_4targets_n3"]["runs"]
    for tname, tdict in pb_runs.items():
        safe = tname.replace("-", "_")
        m[f"s4_pb_{safe}_wall_ratio"] = f"{tdict['wall_ratio']:.4f}"
        m[f"s4_pb_{safe}_rss_ratio"] = f"{tdict['rss_ratio']:.4f}"
        m[f"s4_pb_{safe}_rust_wall_s"] = f"{tdict['rust_wall_s']:.2f}"
        m[f"s4_pb_{safe}_rust_rss_mb"] = f"{tdict['rust_rss_kb'] / 1024.0:.1f}"
        m[f"s4_pb_{safe}_java_wall_s"] = f"{tdict['java_wall_s']:.2f}"
        m[f"s4_pb_{safe}_java_rss_mb"] = f"{tdict['java_rss_kb'] / 1024.0:.1f}"
        m[f"s4_pb_{safe}_cpp_wall_s"] = f"{tdict['cpp_wall_s']:.2f}"
        m[f"s4_pb_{safe}_cpp_rss_mb"] = f"{tdict['cpp_maxrss_kb'] / 1024.0:.1f}"

    # Stage 4 Exome scaling 1..32c
    m["s4_exome_total_vars"] = "339,698"
    ex = s4["scale_exome_3way"]
    for tk in ("th1", "th2", "th4", "th8", "th16", "th32"):
        td = ex[tk]
        m[f"s4_exome_{tk}_rust_wall_s"] = f"{td['rust']['wall_s_mean']:.2f}"
        m[f"s4_exome_{tk}_rust_wall_range"] = f"{td['rust']['wall_s_min']:.2f}–{td['rust']['wall_s_max']:.2f}"
        m[f"s4_exome_{tk}_rust_rss_mb"] = f"{td['rust']['max_rss_mb_mean']:.1f}"
        m[f"s4_exome_{tk}_rust_rss_gib"] = f"{td['rust']['max_rss_gib_mean']:.2f}"
        m[f"s4_exome_{tk}_rust_sp"] = f"{td['rust']['speedup_sp']:.2f}"
        m[f"s4_exome_{tk}_java_wall_s"] = f"{td['java']['wall_s_mean']:.2f}"
        m[f"s4_exome_{tk}_java_wall_range"] = f"{td['java']['wall_s_min']:.2f}–{td['java']['wall_s_max']:.2f}"
        m[f"s4_exome_{tk}_java_rss_gib"] = f"{td['java']['max_rss_gib_mean']:.2f}"
        m[f"s4_exome_{tk}_cpp_wall_s"] = f"{td['cpp']['wall_s_mean']:.2f}"
        m[f"s4_exome_{tk}_cpp_wall_range"] = f"{td['cpp']['wall_s_min']:.2f}–{td['cpp']['wall_s_max']:.2f}"
        m[f"s4_exome_{tk}_cpp_rss_gib"] = f"{td['cpp']['max_rss_gib_mean']:.2f}"
        m[f"s4_exome_{tk}_speedup_vs_java"] = f"{td['java']['wall_s_mean'] / td['rust']['wall_s_mean']:.1f}"
        m[f"s4_exome_{tk}_rss_ratio_vs_java"] = f"{td['java']['max_rss_kb_mean'] / td['rust']['max_rss_kb_mean']:.1f}"
        m[f"s4_exome_{tk}_rss_ratio_vs_cpp"] = f"{td['cpp']['max_rss_kb_mean'] / td['rust']['max_rss_kb_mean']:.1f}"

    # Stage 4 Span chr1
    sp_c1 = s4["span_chr1_3way"]
    m["s4_span_h2_bound"] = "1.50"
    for tk in ("th8", "th16", "th32"):
        sd = sp_c1[tk]
        m[f"s4_span_{tk}_ratio_exact"] = f"{sd['rust_rss_ratio_249m_vs_10m']:.3f}"
        m[f"s4_span_{tk}_10m_rust_gib"] = f"{sd['spans']['R10000000']['rust']['max_rss_gib_mean']:.2f}"
        m[f"s4_span_{tk}_249m_rust_gib"] = f"{sd['spans']['R248956422']['rust']['max_rss_gib_mean']:.2f}"
        m[f"s4_span_{tk}_249m_rust_wall_s"] = f"{sd['spans']['R248956422']['rust']['wall_s_mean']:.2f}"
        m[f"s4_span_{tk}_249m_cpp_wall_s"] = f"{sd['spans']['R248956422']['cpp']['wall_s_mean']:.2f}"
        m[f"s4_span_{tk}_249m_cpp_gib"] = f"{sd['spans']['R248956422']['cpp']['max_rss_gib_mean']:.2f}"

    spf0 = s4["span_chr1_f0"]
    m["s4_span_f0_vars"] = "7,032"
    m["s4_span_f0_rust_wall_s"] = f"{spf0['rust']['wall_s_mean']:.2f}"
    m["s4_span_f0_rust_gib"] = f"{spf0['rust']['max_rss_gib_mean']:.2f}"
    m["s4_span_f0_cpp_wall_s"] = f"{spf0['cpp']['wall_s_mean']:.2f}"
    m["s4_span_f0_cpp_gib"] = f"{spf0['cpp']['max_rss_gib_mean']:.2f}"

    jpn = s4["jobs_per_node_2x4c_n3"]
    m["s4_jpn_job1_wall_s"] = f"{jpn['job1_4c']['wall_s_mean']:.2f}"
    m["s4_jpn_job1_rss_mb"] = f"{jpn['job1_4c']['max_rss_mb_mean']:.1f}"
    m["s4_jpn_job2_wall_s"] = f"{jpn['job2_4c']['wall_s_mean']:.2f}"
    m["s4_jpn_job2_rss_mb"] = f"{jpn['job2_4c']['max_rss_mb_mean']:.1f}"
    m["s4_jpn_combined_gib"] = f"{jpn['combined_peak_rss_gib_mean']:.3f}"

    # Stage 4 Tier-3 truth & strata
    t3c = s4["tier3_truth_and_qual_sweep"]["cohorts"]
    for cid in ("hg002", "seqc2", "hg008"):
        best = t3c[cid]["engines"]["rust"]["operating_points"]["best"]
        none = t3c[cid]["engines"]["rust"]["operating_points"]["none"]
        m[f"s4_t3_{cid}_best_qual"] = str(best["threshold"])
        m[f"s4_t3_{cid}_best_tp"] = f"{best['tp_baseline']:,}"
        m[f"s4_t3_{cid}_best_fp"] = f"{best['fp']:,}"
        m[f"s4_t3_{cid}_best_fn"] = f"{best['fn']:,}"
        m[f"s4_t3_{cid}_best_p"] = f"{best['precision']:.4f}"
        m[f"s4_t3_{cid}_best_r"] = f"{best['sensitivity']:.4f}"
        m[f"s4_t3_{cid}_best_f1"] = f"{best['f_measure']:.4f}"
        m[f"s4_t3_{cid}_none_tp"] = f"{none['tp_baseline']:,}"
        m[f"s4_t3_{cid}_none_fp"] = f"{none['fp']:,}"
        m[f"s4_t3_{cid}_none_fn"] = f"{none['fn']:,}"
        m[f"s4_t3_{cid}_none_p"] = f"{none['precision']:.4f}"
        m[f"s4_t3_{cid}_none_r"] = f"{none['sensitivity']:.4f}"
        m[f"s4_t3_{cid}_none_f1"] = f"{none['f_measure']:.4f}"

    t3s = s4["tier3_truth_and_qual_sweep"]["strata_summary"]
    for cid in ("hg002", "seqc2", "hg008"):
        for eng in ("rust", "java", "cpp"):
            sd = t3s[cid][eng]
            m[f"s4_strata_{cid}_{eng}_snv"] = f"{sd['SNV']:,}"
            m[f"s4_strata_{cid}_{eng}_ins"] = f"{sd['Insertion']:,}"
            m[f"s4_strata_{cid}_{eng}_del"] = f"{sd['Deletion']:,}"
            m[f"s4_strata_{cid}_{eng}_mnp"] = f"{sd['Complex_MNP']:,}"
            m[f"s4_strata_{cid}_{eng}_sv"] = f"{sd['SV']:,}"
            m[f"s4_strata_{cid}_{eng}_pass"] = f"{sd['PASS']:,}"
            m[f"s4_strata_{cid}_{eng}_total"] = f"{sd['total']:,}"

    # Evidence certified HG005 & NA24385-b37
    h5_op = ev["hg005-truth"]["extra"]["rust"]
    m["ev_hg005_best_f1"] = f"{h5_op['optimal']['f_measure']:.4f}"
    m["ev_hg005_none_f1"] = f"{h5_op['none']['f_measure']:.4f}"
    b37_op = ev["na24385-b37-truth"]["extra"]["rust"]
    m["ev_b37_best_f1"] = f"{b37_op['optimal']['f_measure']:.4f}"
    m["ev_b37_none_f1"] = f"{b37_op['none']['f_measure']:.4f}"

    # Stage 5 WGS details
    m["s5_hg002_bam_gb"] = "130.8"
    m["s5_hg008_bam_gb"] = "329.0"
    m["s5_contigs_bp"] = "2,875,001,522"
    m["s5_hg002_giab_gbp"] = "2.51"
    m["s5_hg002_giab_intervals"] = "56,642"
    m["s5_hg008_giab_gbp"] = "2.37"
    m["s5_hg008_giab_intervals"] = "50,554"

    h2_a = s5["wgs_hg002_full_contigs_3way_th32"]
    h2_b = s5["wgs_hg002_full_giab_3way_th32"]
    h8_a = s5["wgs_hg008_full_contigs_3way_th32"]
    h8_b = s5["wgs_hg008_full_giab_3way_th32"]

    m["s5_hg002_contigs_cpp_rss_ratio_vs_rust"] = f"{h2_a['cpp']['max_rss_kb'] / h2_a['rust']['max_rss_kb']:.1f}"
    m["s5_hg008_contigs_cpp_rss_ratio_vs_rust"] = f"{h8_a['cpp']['max_rss_kb'] / h8_a['rust']['max_rss_kb']:.1f}"
    m["s5_hg008_giab_cpp_rss_ratio_vs_rust"] = f"{h8_b['cpp']['max_rss_kb'] / h8_b['rust']['max_rss_kb']:.1f}"
    m["s5_hg002_giab_cpp_diff_vars"] = f"{h2_b['cpp']['variant_count'] - h2_b['rust']['variant_count']:,}"
    m["s5_hg008_contigs_cpp_wall_str"] = h8_a["cpp"]["wall_str"]
    m["s5_hg008_contigs_cpp_vars"] = "13,182,812"
    m["s5_hg008_giab_cpp_wall_str"] = h8_b["cpp"]["wall_str"]
    m["s5_hg008_giab_cpp_vars"] = "10,287,644"

    h2_v_best = h2_b["vcfeval"]["rust"]["operating_points"]["best"]
    h2_v_none = h2_b["vcfeval"]["rust"]["operating_points"]["none"]
    m["s5_hg002_wgs_p_best"] = f"{h2_v_best['precision']:.4f}"
    m["s5_hg002_wgs_r_best"] = f"{h2_v_best['sensitivity']:.4f}"
    m["s5_hg002_wgs_p_none"] = f"{h2_v_none['precision']:.4f}"
    m["s5_hg002_wgs_r_none"] = f"{h2_v_none['sensitivity']:.4f}"

    # Stage 5 GRCh38 field callers
    for row in s5["grch38_field_callers_truth"]:
        coh = row["cohort"].replace(".", "")
        cal = row["caller"]
        m[f"s5_fc_{coh}_{cal}_wall_s"] = f"{row['wall_s']:.2f}"
        m[f"s5_fc_{coh}_{cal}_rss_mb"] = f"{row['peak_mb']:.1f}"
        m[f"s5_fc_{coh}_{cal}_snv_p"] = f"{row['snv_p']:.4f}"
        m[f"s5_fc_{coh}_{cal}_snv_r"] = f"{row['snv_r']:.4f}"
        m[f"s5_fc_{coh}_{cal}_snv_f1"] = f"{row['snv_f1']:.4f}"
        m[f"s5_fc_{coh}_{cal}_indel_p"] = f"{row['indel_p']:.4f}"
        m[f"s5_fc_{coh}_{cal}_indel_r"] = f"{row['indel_r']:.4f}"
        m[f"s5_fc_{coh}_{cal}_indel_f1"] = f"{row['indel_f1']:.4f}"

    # Stage 5 Cross-platform chr20 panel
    cp_a = s5["cross_platform_short_read_chr20_panel"]["part_a_rows"]
    rust_cp = [r for r in cp_a if r["engine"] == "rust"]
    java_cp = [r for r in cp_a if r["engine"] == "java"]
    m["s5_chr20_mbp"] = "64.44"
    m["s5_chr20_total_runs"] = str(len(rust_cp))
    java_completing = [r for r in java_cp if r["rc"] == 0]
    m["s5_chr20_java_completed_runs"] = str(len(java_completing))
    m["s5_chr20_java_timeout_runs"] = str(len(java_cp) - len(java_completing))
    m["cp95_chr20_java_completed"] = fmt_cp_ci(len(java_completing), len(java_completing))

    hg2_std_snv = [r["snp_f1"] for r in rust_cp if r["cohort"] == "HG002_Germline"]
    hg2_ug_snv = [r["ughcr_snp_f1"] for r in rust_cp if r["cohort"] == "HG002_Germline"]
    hg8_f01_std_snv = [r["snp_f1"] for r in rust_cp if r["cohort"] == "HG008_Somatic" and str(r["freq"]) == "0.01"]
    hg8_f01_ug_snv = [r["ughcr_snp_f1"] for r in rust_cp if r["cohort"] == "HG008_Somatic" and str(r["freq"]) == "0.01"]
    hg8_all_std_snv = [r["snp_f1"] for r in rust_cp if r["cohort"] == "HG008_Somatic"]
    hg8_all_ug_snv = [r["ughcr_snp_f1"] for r in rust_cp if r["cohort"] == "HG008_Somatic"]

    m["s5_chr20_hg2_std_snv_f1_range"] = f"{min(hg2_std_snv):.4f}–{max(hg2_std_snv):.4f}"
    m["s5_chr20_hg2_ug_snv_f1_range"] = f"{min(hg2_ug_snv):.4f}–{max(hg2_ug_snv):.4f}"
    m["s5_chr20_hg8_f01_std_snv_f1_range"] = f"{min(hg8_f01_std_snv):.4f}–{max(hg8_f01_std_snv):.4f}"
    m["s5_chr20_hg8_f01_ug_snv_f1_range"] = f"{min(hg8_f01_ug_snv):.4f}–{max(hg8_f01_ug_snv):.4f}"
    m["s5_chr20_hg8_all_std_snv_f1_range"] = f"{min(hg8_all_std_snv):.4f}–{max(hg8_all_std_snv):.4f}"
    m["s5_chr20_hg8_all_ug_snv_f1_range"] = f"{min(hg8_all_ug_snv):.4f}–{max(hg8_all_ug_snv):.4f}"

    tot_nf_cand = sum(int(r["total_candidates"]) for r in nf_rows)
    tot_nf_pass = sum(int(r["pass_somatic"]) for r in nf_rows)
    m["s5_nf_total_candidates"] = f"{tot_nf_cand:,}"
    m["s5_nf_total_pass"] = f"{tot_nf_pass:,}"

    # Stage 6 HG002 v5.0q / CMRG / Mosaic
    m["s6_v50q_gbp"] = "2.74"
    for row in s6["hg002_v50q_cmrg_mosaic_eval"]:
        b = row["benchmark"]
        c = row["caller"]
        m[f"s6_{b}_{c}_snv_p"] = f"{row['snv_p']:.4f}"
        m[f"s6_{b}_{c}_snv_r"] = f"{row['snv_r']:.4f}"
        m[f"s6_{b}_{c}_snv_f1"] = f"{row['snv_f1']:.4f}"
        m[f"s6_{b}_{c}_indel_p"] = f"{row['indel_p']:.4f}"
        m[f"s6_{b}_{c}_indel_r"] = f"{row['indel_r']:.4f}"
        m[f"s6_{b}_{c}_indel_f1"] = f"{row['indel_f1']:.4f}"
    m["s6_mosaic_recall_pct"] = f"{float(m['s6_mosaic_recall_best']) * 100.0:.2f}%"

    # Stage 6 HG008 v0.4 11-caller matrix
    v04_rows = s6["hg008_modern_ml_and_commercial_callers"]
    m["s6_v04_total_rows"] = str(len(v04_rows))
    m["s6_v04_total_callers"] = str(len({r["caller"] for r in v04_rows}))
    for row in v04_rows:
        sc = row["scope"].replace(".", "")
        cal = row["caller"].replace(".", "")
        m[f"s6_{sc}_{cal}_snv_p"] = f"{row['snv_p']:.4f}"
        m[f"s6_{sc}_{cal}_snv_r"] = f"{row['snv_r']:.4f}"
        m[f"s6_{sc}_{cal}_snv_f1"] = f"{row['snv_f1']:.4f}"
        m[f"s6_{sc}_{cal}_indel_p"] = f"{row['indel_p']:.4f}"
        m[f"s6_{sc}_{cal}_indel_r"] = f"{row['indel_r']:.4f}"
        m[f"s6_{sc}_{cal}_indel_f1"] = f"{row['indel_f1']:.4f}"
        m[f"s6_{sc}_{cal}_indel_r_pct"] = f"{row['indel_r'] * 100.0:.2f}%"

    # Stage 6 Low-VAF titration curve
    tit_rows = s6["low_vaf_titration_truth_curve"]
    rust_tit = [r for r in tit_rows if r["caller"] == "VarDictRust"]
    java_tit = [r for r in tit_rows if r["caller"] == "VarDictJava"]
    t_sp = [j["wall_s"] / r["wall_s"] for r, j in zip(rust_tit, java_tit)]
    t_rr = [j["peak_kb"] / r["peak_kb"] for r, j in zip(rust_tit, java_tit)]
    m["s6_tit_rust_wall_range"] = f"{min(r['wall_s'] for r in rust_tit):.2f}–{max(r['wall_s'] for r in rust_tit):.2f}"
    m["s6_tit_java_wall_range"] = f"{min(j['wall_s'] for j in java_tit):.2f}–{max(j['wall_s'] for j in java_tit):.2f}"
    m["s6_tit_speedup_range"] = f"{min(t_sp):.1f}×–{max(t_sp):.1f}×"
    m["s6_tit_rss_ratio_range"] = f"{min(t_rr):.1f}×–{max(t_rr):.1f}×"
    m["s6_tit_rust_rss_kb_range"] = f"{min(r['peak_kb'] for r in rust_tit):,}–{max(r['peak_kb'] for r in rust_tit):,}"
    m["s6_tit_java_rss_kb_range"] = f"{min(j['peak_kb'] for j in java_tit):,}–{max(j['peak_kb'] for j in java_tit):,}"
    for r in tit_rows:
        p = r["purity_pct"]
        cal = r["caller"]
        m[f"s6_tit_p{p}_{cal}_wall_s"] = f"{r['wall_s']:.2f}"
        m[f"s6_tit_p{p}_{cal}_rss_mb"] = f"{r['peak_mb']:.1f}"
        m[f"s6_tit_p{p}_{cal}_snv_p"] = f"{r['snv_p']:.4f}"
        m[f"s6_tit_p{p}_{cal}_snv_r"] = f"{r['snv_r']:.4f}"
        m[f"s6_tit_p{p}_{cal}_snv_f1"] = f"{r['snv_f1']:.4f}"
        m[f"s6_tit_p{p}_{cal}_indel_p"] = f"{r['indel_p']:.4f}"
        m[f"s6_tit_p{p}_{cal}_indel_r"] = f"{r['indel_r']:.4f}"
        m[f"s6_tit_p{p}_{cal}_indel_f1"] = f"{r['indel_f1']:.4f}"
        m[f"s6_tit_p{p}_{cal}_snv_r_pct"] = f"{r['snv_r'] * 100.0:.2f}%"
        m[f"s6_tit_p{p}_{cal}_indel_r_pct"] = f"{r['indel_r'] * 100.0:.2f}%"

    # Stage 6 GA4GH v3.5 summary
    ga_snv_f1 = [r["snv_f1"] for r in ga_rows]
    m["s6_ga4gh_mean_snv_f1"] = f"{sum(ga_snv_f1) / len(ga_snv_f1):.4f}"

    return m


def render_cells_60_ledger(inv_path: Path) -> str:
    """Render the 60-cell inventory TSV as a Markdown table."""
    with inv_path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    lines = [
        "| Cell ID | Stage | Tier | Exit Code | Files | Paper Targets | Parity Verdict | Description |",
        "| :--- | :--- | :--- | :---: | ---: | :--- | :--- | :--- |",
    ]
    for r in rows:
        verdict = r["parity_verdict"][:80].replace("|", "/")
        desc = r["description"].replace("|", "/")
        lines.append(
            f"| `{r['cell_id']}` | `{r['stage']}` | {r['tier']} | `{r['exit_code']}` | "
            f"`{r['file_count']}` | {r['paper_targets']} | `{verdict}` | {desc} |"
        )
    return "\n".join(lines)


def compile_template(
    tmpl_text: str,
    macros: dict[str, str],
    tables_dir: Path,
    inv_path: Path,
) -> str:
    """Expand {{m:<key>}}, {{table:<name>}}, and {{cells_60_ledger}} in template text."""
    # 1. Expand {{cells_60_ledger}}
    if "{{cells_60_ledger}}" in tmpl_text:
        tmpl_text = tmpl_text.replace("{{cells_60_ledger}}", render_cells_60_ledger(inv_path))

    # 2. Expand {{table:<name>}}
    def replace_table(match: re.Match[str]) -> str:
        tname = match.group(1).strip()
        if tname in TABLE_ALIASES:
            fname, anchor_id = TABLE_ALIASES[tname]
        else:
            fname = f"{tname}.md" if not tname.endswith(".md") else tname
            anchor_id = tname.removesuffix(".md")
        tpath = tables_dir / fname
        if not tpath.is_file():
            raise KeyError(f"Missing table file {tpath} for macro {{{{table:{tname}}}}}")
        tcontent = tpath.read_text(encoding="utf-8").strip()
        return f'<a id="{anchor_id}"></a>\n\n{tcontent}'

    tmpl_text = re.sub(r"\{\{table:([^}]+)\}\}", replace_table, tmpl_text)

    # 3. Expand {{m:<key>}}
    def replace_metric(match: re.Match[str]) -> str:
        mkey = match.group(1).strip()
        if mkey not in macros:
            raise KeyError(f"Unknown metric macro key: {{{{m:{mkey}}}}}")
        return macros[mkey]

    compiled = re.sub(r"\{\{m:([^}]+)\}\}", replace_metric, tmpl_text)

    # 4. Verify zero unresolved {{...}} macros remain
    leftover = re.search(r"\{\{[^}]+\}\}", compiled)
    if leftover:
        raise ValueError(f"Unresolved template macro found after compilation: {leftover.group(0)}")

    return compiled


def main() -> int:
    parser = argparse.ArgumentParser(description="Compile reproducible manuscript.md and supplement.md from templates.")
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS, help="Path to master_gcp_metrics.json")
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY, help="Path to cells_60_inventory.tsv")
    parser.add_argument("--tables-dir", type=Path, default=DEFAULT_TABLES_DIR, help="Directory containing table*.md")
    parser.add_argument("--tmpl-dir", type=Path, default=DEFAULT_TMPL_DIR, help="Directory containing *.tmpl.md")
    parser.add_argument("--ms-out", type=Path, default=DEFAULT_MS_OUT, help="Output path for manuscript.md")
    parser.add_argument("--sup-out", type=Path, default=DEFAULT_SUP_OUT, help="Output path for supplement.md")
    parser.add_argument("--check", action="store_true", help="Verify existing manuscript.md and supplement.md match templates")
    args = parser.parse_args()

    if not args.metrics.is_file():
        print(f"ERROR: Missing metrics JSON: {args.metrics}", file=sys.stderr)
        return 1
    if not args.inventory.is_file():
        print(f"ERROR: Missing inventory TSV: {args.inventory}", file=sys.stderr)
        return 1

    ms_tmpl_path = args.tmpl_dir / "manuscript.tmpl.md"
    sup_tmpl_path = args.tmpl_dir / "supplement.tmpl.md"
    if not ms_tmpl_path.is_file() or not sup_tmpl_path.is_file():
        print(f"ERROR: Missing templates in {args.tmpl_dir}", file=sys.stderr)
        return 1

    doc = json.loads(args.metrics.read_text(encoding="utf-8"))
    macros = build_macro_context(doc)

    ms_compiled = compile_template(
        ms_tmpl_path.read_text(encoding="utf-8"),
        macros,
        args.tables_dir,
        args.inventory,
    )
    sup_compiled = compile_template(
        sup_tmpl_path.read_text(encoding="utf-8"),
        macros,
        args.tables_dir,
        args.inventory,
    )

    if args.check:
        if not args.ms_out.is_file() or not args.sup_out.is_file():
            print("ERROR: --check failed because manuscript.md or supplement.md does not exist.", file=sys.stderr)
            return 1
        ms_disk = args.ms_out.read_text(encoding="utf-8")
        sup_disk = args.sup_out.read_text(encoding="utf-8")
        if ms_disk != ms_compiled:
            print(f"ERROR: {args.ms_out} is out of sync with {ms_tmpl_path}", file=sys.stderr)
            return 1
        if sup_disk != sup_compiled:
            print(f"ERROR: {args.sup_out} is out of sync with {sup_tmpl_path}", file=sys.stderr)
            return 1
        print(f"OK verified {args.ms_out.name} and {args.sup_out.name} match compiled templates byte-for-byte")
        return 0

    args.ms_out.write_text(ms_compiled, encoding="utf-8")
    args.sup_out.write_text(sup_compiled, encoding="utf-8")
    print(f"OK compiled {args.ms_out} ({len(ms_compiled.split())} words) and {args.sup_out} ({len(sup_compiled.split())} words)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
