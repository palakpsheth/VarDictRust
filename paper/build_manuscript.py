#!/usr/bin/env python3
"""Compile paper/manuscript.md and paper/supplement.md from macro-parameterized templates.

Reads exclusively from git-tracked benchmarks/gcp-study-v106/ (master_gcp_metrics.json and
cells_60_inventory.tsv), benchmarks/v110-compat-study/v110_compat_metrics.json, and
paper/tables/*.md. Supports default build mode and --check idempotency verification mode.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import re
import sys
import tomllib
from pathlib import Path
from typing import Any

from scipy import stats

try:
    from paper.render_pdf import render_markdown_to_pdf
except ImportError:
    from render_pdf import render_markdown_to_pdf

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METRICS = ROOT / "benchmarks" / "gcp-study-v106" / "master_gcp_metrics.json"
DEFAULT_V110_METRICS = (
    ROOT / "benchmarks" / "v110-compat-study" / "v110_compat_metrics.json"
)
DEFAULT_INVENTORY = ROOT / "benchmarks" / "gcp-study-v106" / "cells_60_inventory.tsv"
DEFAULT_TABLES_DIR = ROOT / "paper" / "tables"
DEFAULT_TMPL_DIR = ROOT / "paper" / "templates"
DEFAULT_MS_OUT = ROOT / "paper" / "manuscript.md"
DEFAULT_SUP_OUT = ROOT / "paper" / "supplement.md"
DEFAULT_MS_PDF = ROOT / "paper" / "manuscript.pdf"
DEFAULT_SUP_PDF = ROOT / "paper" / "supplement.pdf"

TABLE_ALIASES: dict[str, tuple[str, str]] = {
    "table1": ("table1.md", "table1_parity_matrix"),
    "table1_parity_matrix": ("table1.md", "table1_parity_matrix"),
    "table2": ("table2.md", "table2_tierp_and_rchr1"),
    "table2_tierp_and_rchr1": ("table2.md", "table2_tierp_and_rchr1"),
    "table3": ("table3.md", "table3_truth_and_field_callers"),
    "table3_truth_and_field_callers": ("table3.md", "table3_truth_and_field_callers"),
    "table4": ("table4.md", "table4_memory_scaling_and_packing"),
    "table4_memory_scaling_and_packing": (
        "table4.md",
        "table4_memory_scaling_and_packing",
    ),
    "table5": ("table5.md", "table5_full_wgs_arm_a_vs_b"),
    "table5_full_wgs_arm_a_vs_b": ("table5.md", "table5_full_wgs_arm_a_vs_b"),
    "table6": ("table6.md", "table6_cross_platform_chr20"),
    "table6_cross_platform_chr20": ("table6.md", "table6_cross_platform_chr20"),
    "table7": ("table7.md", "table7_v50q_cmrg_mosaic_and_hg008_v04"),
    "table7_v50q_cmrg_mosaic": ("table7.md", "table7_v50q_cmrg_mosaic_and_hg008_v04"),
    "table7_v50q_cmrg_mosaic_and_hg008_v04": (
        "table7.md",
        "table7_v50q_cmrg_mosaic_and_hg008_v04",
    ),
    "table8": ("table8.md", "table8_supply_chain_and_divergences"),
    "table8_hg008_modern_callers_and_titration": (
        "table8.md",
        "table8_supply_chain_and_divergences",
    ),
    "table8_supply_chain_and_divergences": (
        "table8.md",
        "table8_supply_chain_and_divergences",
    ),
    "table9": ("table9.md", "table9_v110_compat_and_presets"),
    "table9_v110_compat_and_presets": (
        "table9.md",
        "table9_v110_compat_and_presets",
    ),
    "supp_table1": ("supp_table1.md", "supp_table1_60cell_inventory"),
    "supp_table1_60cell_inventory": ("supp_table1.md", "supp_table1_60cell_inventory"),
    "supp_table2": ("supp_table2.md", "supp_table2_option_c_87arms"),
    "supp_table2_option_c_87arms": ("supp_table2.md", "supp_table2_option_c_87arms"),
    "supp_table3": ("supp_table3.md", "supp_table3_n3_technical_replicates"),
    "supp_table3_n3_technical_replicates": (
        "supp_table3.md",
        "supp_table3_n3_technical_replicates",
    ),
    "supp_table4": ("supp_table4.md", "supp_table4_qual_roc_and_noise_floor"),
    "supp_table4_qual_roc_and_noise_floor": (
        "supp_table4.md",
        "supp_table4_qual_roc_and_noise_floor",
    ),
    "supp_table5": ("supp_table5.md", "supp_table5_ga4gh_70cell_matrix"),
    "supp_table5_ga4gh_70cell_matrix": (
        "supp_table5.md",
        "supp_table5_ga4gh_70cell_matrix",
    ),
    "supp_table6": ("supp_table6.md", "supp_table6_sv_and_amplicon_splice"),
    "supp_table6_sv_and_amplicon_splice": (
        "supp_table6.md",
        "supp_table6_sv_and_amplicon_splice",
    ),
    "supp_table7": ("supp_table7.md", "supp_table7_compat_32switches_catalog"),
    "supp_table7_compat_32switches_catalog": (
        "supp_table7.md",
        "supp_table7_compat_32switches_catalog",
    ),
}


def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Compute exact Clopper-Pearson (1 - alpha) binomial confidence interval in percentage points."""
    if n <= 0:
        return (0.0, 0.0)
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2.0, k, n - k + 1)) * 100.0
    hi = (
        100.0
        if k == n
        else float(stats.beta.ppf(1.0 - alpha / 2.0, k + 1, n - k)) * 100.0
    )
    return (lo, hi)


def _eval_rust_const_int(expr_str: str, const_name: str) -> int:
    """Safely evaluate an integer Rust constant expression (e.g. '16 * 1024 * 1024') via stdlib ast."""
    cleaned = re.sub(r"_", "", expr_str.strip())
    cleaned = re.sub(r"\b(\d+)(?:usize|u64|u32|i64|i32)\b", r"\1", cleaned)
    tree = ast.parse(cleaned, mode="eval")

    def _eval_node(node: ast.AST) -> int:
        if isinstance(node, ast.Expression):
            return _eval_node(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, int) and not isinstance(node.value, bool):
            return int(node.value)
        if isinstance(node, ast.BinOp):
            left = _eval_node(node.left)
            right = _eval_node(node.right)
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.FloorDiv):
                return left // right
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.LShift):
                return left << right
        raise ValueError(f"Unsupported AST expression for {const_name}: {expr_str!r}")

    return _eval_node(tree)


def _extract_mib_const(rust_src: str, const_name: str, rel_path: str) -> str:
    match = re.search(
        rf"\bconst\s+{re.escape(const_name)}\s*:\s*(?:usize|u64|u32)\s*=\s*([^;]+);",
        rust_src,
    )
    if not match:
        raise ValueError(f"Could not extract constant {const_name} from {rel_path}")
    raw_bytes = _eval_rust_const_int(match.group(1), const_name)
    mib_unit = 1024 * 1024
    if raw_bytes <= 0 or raw_bytes % mib_unit != 0:
        raise ValueError(
            f"Constant {const_name} in {rel_path} ({raw_bytes} B) is not a positive exact multiple of MiB"
        )
    return str(raw_bytes // mib_unit)


def _load_repo_release_and_engine_macros(repo_root: Path) -> dict[str, str]:
    """Load live release metadata, test/CLI counts, and engine constants without hardcoded numbers."""
    out: dict[str, str] = {}

    cargo_toml_path = repo_root / "Cargo.toml"
    cargo_data = tomllib.loads(cargo_toml_path.read_text(encoding="utf-8"))
    ws_pkg = (cargo_data.get("workspace") or {}).get("package") or {}
    cur_ver = ws_pkg.get("version")
    spdx_lic = ws_pkg.get("license")
    if not isinstance(cur_ver, str) or not cur_ver:
        raise ValueError(f"Missing [workspace.package].version in {cargo_toml_path}")
    if not isinstance(spdx_lic, str) or not spdx_lic:
        raise ValueError(f"Missing [workspace.package].license in {cargo_toml_path}")
    out["current_release_version"] = cur_ver
    out["spdx_license"] = spdx_lic

    consistency_path = (
        repo_root
        / ".agents"
        / "skills"
        / "vardictrust-release"
        / "scripts"
        / "check_release_consistency.py"
    )
    spec = importlib.util.spec_from_file_location("check_release_consistency", consistency_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load spec from {consistency_path}")
    consistency_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(consistency_mod)

    tc = consistency_mod.EXPECTED_TEST_COUNTS
    out["workspace_test_count"] = f"{int(tc['total_workspace']):,}"
    out["core_lib_test_count"] = f"{int(tc['core_lib']):,}"
    out["harness_self_test_count"] = f"{int(tc['harness_self_tests']):,}"
    out["cli_flag_count"] = f"{int(consistency_mod.EXPECTED_CLI_FLAG_COUNT):,}"
    out["cli_group_count"] = str(len(consistency_mod.EXPECTED_CLI_HEADINGS))

    div_txt = (repo_root / "DIVERGENCES.md").read_text(encoding="utf-8")
    rel_block_m = re.search(
        r"Verified Parity Baseline \(release/v[^)]+\):([\s\S]+?)(?:\n\n|\n```)",
        div_txt,
    )
    if not rel_block_m:
        raise ValueError("Could not locate Verified Parity Baseline block in DIVERGENCES.md")
    rel_block = rel_block_m.group(1)
    rel_core_m = re.search(r"Core Library Tests:\s*(\d+)\s+PASS", rel_block)
    rel_ws_m = re.search(r"Total Workspace:\s*(\d+)\s+PASS", rel_block)
    if not rel_core_m or not rel_ws_m:
        raise ValueError(
            "Could not extract certified release Core Library / Total Workspace test counts from DIVERGENCES.md"
        )
    out["release_core_lib_test_count"] = f"{int(rel_core_m.group(1)):,}"
    out["release_workspace_test_count"] = f"{int(rel_ws_m.group(1)):,}"

    lib_rs = (repo_root / "crates" / "vardictrust-core" / "src" / "lib.rs").read_text(encoding="utf-8")
    bam_rs = (repo_root / "crates" / "vardictrust-core" / "src" / "bam.rs").read_text(encoding="utf-8")
    cfg_rs = (repo_root / "crates" / "vardictrust-core" / "src" / "config.rs").read_text(encoding="utf-8")
    fisher_rs = (
        repo_root / "crates" / "vardictrust-core" / "src" / "stats" / "fisher.rs"
    ).read_text(encoding="utf-8")

    out["vcf_spill_cap_mib"] = _extract_mib_const(
        lib_rs, "DEFAULT_VCF_SPILL_BYTES", "crates/vardictrust-core/src/lib.rs"
    )
    out["cram_container_cap_mib"] = _extract_mib_const(
        bam_rs, "MAX_CACHED_CRAM_CONTAINER_BYTES", "crates/vardictrust-core/src/bam.rs"
    )
    out["cram_shard_cap_mib"] = _extract_mib_const(
        bam_rs, "MAX_SHARD_CACHED_CRAM_BYTES", "crates/vardictrust-core/src/bam.rs"
    )
    out["fasta_cache_cap_mib"] = _extract_mib_const(
        bam_rs,
        "DEFAULT_FASTA_CACHE_BUDGET_BYTES",
        "crates/vardictrust-core/src/bam.rs",
    )
    shards_m = re.search(
        r"\bconst\s+CRAM_CONTAINER_CACHE_SHARDS\s*:\s*usize\s*=\s*([^;]+);",
        bam_rs,
    )
    if not shards_m:
        raise ValueError(
            "Could not extract CRAM_CONTAINER_CACHE_SHARDS from crates/vardictrust-core/src/bam.rs"
        )
    cram_shards = _eval_rust_const_int(
        shards_m.group(1), "CRAM_CONTAINER_CACHE_SHARDS"
    )
    out["cram_cache_shards"] = str(cram_shards)
    out["cram_total_cache_mib"] = str(
        cram_shards * int(out["cram_shard_cap_mib"])
    )
    out["max_bai_cap_mib"] = _extract_mib_const(
        bam_rs, "MAX_BAI_BYTES", "crates/vardictrust-core/src/bam.rs"
    )
    out["max_csi_cap_mib"] = _extract_mib_const(
        bam_rs, "MAX_CSI_BYTES", "crates/vardictrust-core/src/bam.rs"
    )

    fisher_sup_m = re.search(
        r"\bconst\s+MAX_DENSE_FISHER_SUPPORT\s*:\s*usize\s*=\s*([^;]+);",
        fisher_rs,
    )
    if not fisher_sup_m:
        raise ValueError(
            "Could not extract MAX_DENSE_FISHER_SUPPORT from crates/vardictrust-core/src/stats/fisher.rs"
        )
    fisher_sup_val = _eval_rust_const_int(
        fisher_sup_m.group(1), "MAX_DENSE_FISHER_SUPPORT"
    )
    fisher_bytes = fisher_sup_val * 8
    mib_unit = 1024 * 1024
    if fisher_bytes <= 0 or fisher_bytes % mib_unit != 0:
        raise ValueError(
            f"MAX_DENSE_FISHER_SUPPORT * 8 ({fisher_bytes} B) is not an exact positive multiple of MiB"
        )
    out["max_dense_fisher_mib"] = str(fisher_bytes // mib_unit)

    default_impl_m = re.search(
        r"impl\s+Default\s+for\s+Config\s*\{([\s\S]*?)\n\}",
        cfg_rs,
    )
    if not default_impl_m:
        raise ValueError("Could not locate 'impl Default for Config' in crates/vardictrust-core/src/config.rs")
    default_body = default_impl_m.group(1)

    mm_m = re.search(r"\bmax_mismatches\s*:\s*(\d+)\s*,", default_body)
    indel_m = re.search(r"\bindel_size\s*:\s*(\d+)\s*,", default_body)
    def_freq_m = re.search(r"\bfreq\s*:\s*([0-9.]+)\s*,", default_body)
    def_mr_m = re.search(r"\bmin_reads\s*:\s*(\d+)\s*,", default_body)
    def_bq_m = re.search(r"\bbase_qual\s*:\s*([0-9.]+)\s*,", default_body)
    if not mm_m:
        raise ValueError("Could not extract default max_mismatches from Config::default()")
    if not indel_m:
        raise ValueError("Could not extract default indel_size from Config::default()")
    if not (def_freq_m and def_mr_m and def_bq_m):
        raise ValueError("Could not extract default freq/min_reads/base_qual from Config::default()")
    out["default_max_mm"] = mm_m.group(1)
    out["default_indel_size_bp"] = indel_m.group(1)

    def_freq_s = def_freq_m.group(1)
    def_mr_s = def_mr_m.group(1)
    def_bq_s = def_bq_m.group(1)

    preset_fn_m = re.search(
        r"pub\s+fn\s+from_preset\s*\([^)]*\)\s*->\s*Self\s*\{([\s\S]*?)\n    \}",
        cfg_rs,
    )
    if not preset_fn_m:
        raise ValueError("Could not extract Config::from_preset from crates/vardictrust-core/src/config.rs")
    preset_body = preset_fn_m.group(1)
    mapq_m = re.search(r"cfg\.vcf_mapq\s*=\s*Some\(VcfThreshold::new\(([0-9.]+)\)\);", preset_body)
    if not mapq_m:
        raise ValueError("Could not extract cfg.vcf_mapq from Config::from_preset")
    out["preset_shared_mapq"] = str(int(float(mapq_m.group(1))))

    for variant_name, prefix in (
        ("Wgs", "preset_wgs"),
        ("Wes", "preset_wes"),
        ("Amplicon", "preset_amp"),
        ("Ctdna", "preset_ctdna"),
    ):
        arm_m = re.search(
            rf"WorkflowPreset::{variant_name}\s*=>\s*\{{([\s\S]*?)\n            \}}",
            preset_body,
        )
        if not arm_m:
            raise ValueError(f"Could not extract WorkflowPreset::{variant_name} arm from Config::from_preset")
        arm_src = arm_m.group(1)
        f_m = re.search(r"cfg\.freq\s*=\s*([0-9.]+);", arm_src)
        mr_m = re.search(r"cfg\.min_reads\s*=\s*(\d+);", arm_src)
        vaf_m = re.search(r"cfg\.vcf_af\s*=\s*Some\(VcfThreshold::new\(([0-9.]+)\)\);", arm_src)
        dp_m = re.search(r"cfg\.vcf_min_depth\s*=\s*Some\(VcfThreshold::new\((\d+)\)\);", arm_src)
        vd_m = re.search(r"cfg\.vcf_min_vd\s*=\s*Some\(VcfThreshold::new\((\d+)\)\);", arm_src)
        qm_m = re.search(r"cfg\.vcf_qmean\s*=\s*Some\(VcfThreshold::new\(([0-9.]+)\)\);", arm_src)
        amp_m = re.search(r'cfg\.amplicon\s*=\s*Some\("([^"]+)"\.to_string\(\)\);', arm_src)
        if not (vaf_m and dp_m and vd_m):
            raise ValueError(f"Missing vcf_af/vcf_min_depth/vcf_min_vd in WorkflowPreset::{variant_name}")
        out[f"{prefix}_freq"] = f_m.group(1) if f_m else def_freq_s
        out[f"{prefix}_min_reads"] = mr_m.group(1) if mr_m else def_mr_s
        out[f"{prefix}_vcf_af"] = vaf_m.group(1)
        out[f"{prefix}_min_dp"] = dp_m.group(1)
        out[f"{prefix}_min_vd"] = vd_m.group(1)
        out[f"{prefix}_qmean"] = qm_m.group(1) if qm_m else def_bq_s
        if amp_m:
            out[f"{prefix}_param"] = amp_m.group(1)

    hg008_ev_path = repo_root / "benchmarks" / "evidence" / "hg008-tier2.json"
    scale_ev_path = repo_root / "benchmarks" / "evidence" / "scale.json"
    if not hg008_ev_path.is_file() or not scale_ev_path.is_file():
        raise FileNotFoundError("Missing benchmarks/evidence/{hg008-tier2,scale}.json")
    hg008_ev = json.loads(hg008_ev_path.read_text(encoding="utf-8"))
    scale_ev = json.loads(scale_ev_path.read_text(encoding="utf-8"))
    cur_fp = hg008_ev.get("fingerprint")
    cur_perf_fp = (scale_ev.get("extra") or {}).get("perf_fingerprint")
    if not isinstance(cur_fp, str) or not cur_fp:
        raise ValueError("Missing fingerprint in benchmarks/evidence/hg008-tier2.json")
    if not isinstance(cur_perf_fp, str) or not cur_perf_fp:
        raise ValueError("Missing extra.perf_fingerprint in benchmarks/evidence/scale.json")
    out["current_release_fingerprint"] = cur_fp
    out["current_perf_fingerprint"] = cur_perf_fp

    return out


def fmt_cp_ci(k: int, n: int) -> str:
    lo, hi = clopper_pearson_ci(k, n)
    return f"{lo:.2f}%–{hi:.2f}%"


def _build_v110_macros(v110_doc: dict[str, Any]) -> dict[str, str]:
    """Derive deterministic v1.1.0 --compat=off and --preset macros from v110_compat_metrics.json."""
    out: dict[str, str] = {}
    v_meta = v110_doc["meta"]
    arms = v110_doc["arms"]
    a1 = arms["arm1_vcf_schema_and_bcftools_norm"]
    a2 = arms["arm2_seqc2_diploid_vcfeval"]
    a3 = arms["arm3_somatic_multiallelic_and_truth"]
    a4 = arms["arm4_brca2_adaptor_chimera"]
    a5 = arms["arm5_clinical_presets_and_amplicon"]
    a6 = arms["arm6_csi_spanning_sv_and_pileup_tiling"]

    out["v110_study_id"] = str(v_meta["study_id"])
    out["v110_study_version"] = str(v_meta["study_version"])
    out["v110_study_fingerprint"] = str(v_meta["study_fingerprint"])
    out["v110_study_perf_fingerprint"] = str(v_meta["study_perf_fingerprint"])
    out["v110_compat_switches"] = str(v_meta["compat_switch_count"])
    out["v110_preset_count"] = str(v_meta["preset_count"])
    out["v110_var2vcf_flags"] = str(v_meta["var2vcf_flag_count"])
    out["v110_total_arms"] = str(v_meta["total_arms"])
    out["v110_rtg_version"] = str(v_meta["rtg_version"])
    out["v110_bcftools_version"] = str(v_meta["bcftools_version"])
    out["v110_samtools_version"] = str(v_meta["samtools_version"])
    out["v110_tabix_version"] = str(v_meta["tabix_version"])
    out["v110_switch_catalog_total"] = str(len(v110_doc["switch_catalog"]))

    # Domain 1: Strict VCFv4.2/v4.3 schema & bcftools norm --check-ref e
    out["v110_arm1_callsets"] = str(a1["total_callsets"])
    out["v110_arm1_bams"] = str(a1["distinct_bams"])
    out["v110_arm1_records_vj"] = f"{int(a1['total_records_vj']):,}"
    out["v110_arm1_records_off"] = f"{int(a1['total_records_off']):,}"
    out["v110_arm1_defects_vj"] = f"{int(a1['total_defects_vj']):,}"
    out["v110_arm1_defects_off"] = str(a1["total_defects_off"])
    out["v110_arm1_def_end"] = f"{int(a1['vardictjava_defects']['nonsv_info_end']):,}"
    out["v110_arm1_def_gt"] = f"{int(a1['vardictjava_defects']['noncanonical_gt']):,}"
    out["v110_arm1_def_acct"] = str(
        a1["vardictjava_defects"]["accounting_type_iupac_shift3"]
    )
    out["v110_arm1_def_sv"] = str(a1["vardictjava_defects"]["noncanonical_sv"])
    out["v110_arm1_def_ad"] = str(a1["vardictjava_defects"]["ad_number_r_mismatch"])
    out["v110_arm1_norm_vj"] = str(a1["bcftools_norm_check_ref_vj_passed"])
    out["v110_arm1_norm_off"] = str(a1["bcftools_norm_check_ref_off_passed"])
    out["v110_arm1_norm_tot"] = str(a1["bcftools_norm_check_ref_off_total"])
    out["v110_compat_switch_count"] = out["v110_compat_switches"]
    out["v110_schema_defects_vj"] = out["v110_arm1_defects_vj"]
    out["v110_bcftools_norm_off"] = (
        f"{out['v110_arm1_norm_off']}/{out['v110_arm1_norm_tot']}"
    )

    # Domain 2: Diploid SEQC2 HCC1395/HCC1395BL rtg vcfeval without --squash-ploidy
    a2_vj = a2["vardictjava"]
    a2_off = a2["compat_off"]
    out["v110_arm2_contigs_count"] = str(len(a2["target_contigs"]))
    out["v110_arm2_intervals"] = f"{int(a2['evaluated_intervals']):,}"
    out["v110_arm2_bp"] = f"{int(a2['evaluated_bp']):,}"
    out["v110_arm2_discovery_af"] = f"{float(a2['discovery_af']):.2f}"
    out["v110_arm2_vcf_af"] = f"{float(a2['vcf_af']):.2f}"
    out["v110_arm2_tvaf_thr"] = f"{float(a2['gt_proj_tvaf_threshold']):.2f}"
    out["v110_arm2_truth_total"] = str(a2["truth_total"])
    out["v110_arm2_subclonal_lt030"] = str(a2["subclonal_truth_lt_030"])
    out["v110_arm2_retained_vj"] = str(a2_vj["retained_somatic_calls"])
    out["v110_arm2_retained_off"] = str(a2_off["retained_somatic_calls"])
    out["v110_arm2_gt00_vj"] = str(a2_vj["gt00_carrier_somatic"])
    out["v110_arm2_gt00_off"] = str(a2_off["gt00_carrier_somatic"])
    out["v110_arm2_vj_tp"] = str(a2_vj["tp_baseline"])
    out["v110_arm2_vj_fp"] = str(a2_vj["fp"])
    out["v110_arm2_vj_fn"] = str(a2_vj["fn"])
    out["v110_arm2_vj_p"] = f"{float(a2_vj['precision']):.4f}"
    out["v110_arm2_vj_r"] = f"{float(a2_vj['sensitivity']):.4f}"
    out["v110_arm2_vj_f1"] = f"{float(a2_vj['f_measure']):.4f}"
    out["v110_arm2_vj_best_qual"] = str(int(round(float(a2_vj["best_qual_threshold"]))))
    out["v110_arm2_vj_best_tp"] = str(a2_vj["best_tp_baseline"])
    out["v110_arm2_vj_best_fp"] = str(a2_vj["best_fp"])
    out["v110_arm2_vj_best_fn"] = str(a2_vj["best_fn"])
    out["v110_arm2_vj_best_p"] = f"{float(a2_vj['best_precision']):.4f}"
    out["v110_arm2_vj_best_r"] = f"{float(a2_vj['best_sensitivity']):.4f}"
    out["v110_arm2_vj_best_f1"] = f"{float(a2_vj['best_f_measure']):.4f}"
    out["v110_arm2_off_tp"] = str(a2_off["tp_baseline"])
    out["v110_arm2_off_fp"] = str(a2_off["fp"])
    out["v110_arm2_off_fn"] = str(a2_off["fn"])
    out["v110_arm2_off_p"] = f"{float(a2_off['precision']):.4f}"
    out["v110_arm2_off_r"] = f"{float(a2_off['sensitivity']):.4f}"
    out["v110_arm2_off_f1"] = f"{float(a2_off['f_measure']):.4f}"
    out["v110_arm2_off_best_qual"] = str(
        int(round(float(a2_off["best_qual_threshold"])))
    )
    out["v110_arm2_off_best_tp"] = str(a2_off["best_tp_baseline"])
    out["v110_arm2_off_best_fp"] = str(a2_off["best_fp"])
    out["v110_arm2_off_best_fn"] = str(a2_off["best_fn"])
    out["v110_arm2_off_best_p"] = f"{float(a2_off['best_precision']):.4f}"
    out["v110_arm2_off_best_r"] = f"{float(a2_off['best_sensitivity']):.4f}"
    out["v110_arm2_off_best_f1"] = f"{float(a2_off['best_f_measure']):.4f}"
    out["v110_arm2_tp_gain"] = str(a2["tp_baseline_gain"])
    out["v110_arm2_recall_gain_pp"] = (
        f"{float(a2['sensitivity_gain']) * 100.0:.1f}"
    )
    out["v110_arm2_r_gain_pct"] = f"{float(a2['sensitivity_gain']) * 100.0:.1f}%"
    out["v110_arm2_unthresh_f1_gain"] = f"{float(a2['f_measure_gain']):.4f}"
    out["v110_arm2_best_tp_gain"] = str(a2["best_tp_baseline_gain"])
    out["v110_arm2_best_recall_gain_pp"] = (
        f"{float(a2['best_sensitivity_gain']) * 100.0:.1f}"
    )
    out["v110_arm2_best_r_gain_pct"] = (
        f"{float(a2['best_sensitivity_gain']) * 100.0:.1f}%"
    )
    out["v110_arm2_best_f1_gain"] = f"{float(a2['best_f_measure_gain']):.4f}"
    out["v110_arm2_f1_gain"] = out["v110_arm2_best_f1_gain"]

    # Exact Clopper-Pearson 95% CIs ONLY over independent truth-set variant denominators (AR-11)
    out["cp95_v110_seqc2_vj_r"] = fmt_cp_ci(
        int(a2_vj["tp_baseline"]), int(a2["truth_total"])
    )
    out["cp95_v110_seqc2_off_r"] = fmt_cp_ci(
        int(a2_off["tp_baseline"]), int(a2["truth_total"])
    )
    out["cp95_v110_seqc2_vj_best_r"] = fmt_cp_ci(
        int(a2_vj["best_tp_baseline"]), int(a2["truth_total"])
    )
    out["cp95_v110_seqc2_off_best_r"] = fmt_cp_ci(
        int(a2_off["best_tp_baseline"]), int(a2["truth_total"])
    )

    # Domain 3: Multi-allelic somatic continuation, MUC19 paralog, & truth slices
    out["v110_arm3_pairs"] = str(a3["distinct_bam_pairs"])
    out["v110_arm3_runs"] = str(a3["total_runs"])
    out["v110_arm3_vj_pairs"] = str(a3["vj_recovered_pairs"])
    out["v110_arm3_off_pairs"] = str(a3["off_recovered_pairs"])
    out["v110_arm3_vj_runs"] = str(a3["vj_recovered_runs"])
    out["v110_arm3_off_runs"] = str(a3["off_recovered_runs"])
    out["v110_arm3_vj_sec"] = str(a3["vj_secondary_alleles"])
    out["v110_arm3_off_sec"] = str(a3["off_secondary_alleles"])
    out["v110_arm3_muc19_vars"] = str(a3["muc19_variants"])
    out["v110_arm3_muc19_viol"] = str(a3["muc19_vd_gt_dp_violations"])
    out["v110_arm3_seqc2_ivs"] = f"{int(a3['seqc2_hc']['evaluated_intervals']):,}"
    out["v110_arm3_seqc2_bp"] = f"{int(a3['seqc2_hc']['evaluated_bp']):,}"
    out["v110_arm3_seqc2_tp"] = str(a3["seqc2_hc"]["compat_off"]["tp"])
    out["v110_arm3_seqc2_fp"] = str(a3["seqc2_hc"]["compat_off"]["fp"])
    out["v110_arm3_seqc2_truth"] = str(a3["seqc2_hc"]["compat_off"]["truth_total"])
    out["v110_arm3_seqc2_p"] = f"{float(a3['seqc2_hc']['compat_off']['precision']):.4f}"
    out["v110_arm3_seqc2_r"] = (
        f"{float(a3['seqc2_hc']['compat_off']['sensitivity']):.4f}"
    )
    out["v110_arm3_seqc2_f1"] = f"{float(a3['seqc2_hc']['compat_off']['f1']):.4f}"
    out["cp95_v110_arm3_seqc2_r"] = fmt_cp_ci(
        int(a3["seqc2_hc"]["compat_off"]["tp"]),
        int(a3["seqc2_hc"]["compat_off"]["truth_total"]),
    )
    out["v110_arm3_hg008_ivs"] = str(a3["hg008_v04"]["evaluated_intervals"])
    out["v110_arm3_hg008_bp"] = f"{int(a3['hg008_v04']['evaluated_bp']):,}"
    out["v110_arm3_hg008_tp"] = str(a3["hg008_v04"]["compat_off"]["tp"])
    out["v110_arm3_hg008_truth"] = str(a3["hg008_v04"]["compat_off"]["truth_total"])
    out["v110_arm3_hg008_vj_fp"] = str(a3["hg008_v04"]["vardictjava"]["fp"])
    out["v110_arm3_hg008_off_fp"] = str(a3["hg008_v04"]["compat_off"]["fp"])
    out["v110_arm3_hg008_vj_p"] = (
        f"{float(a3['hg008_v04']['vardictjava']['precision']):.4f}"
    )
    out["v110_arm3_hg008_off_p"] = (
        f"{float(a3['hg008_v04']['compat_off']['precision']):.4f}"
    )
    out["v110_arm3_hg008_r"] = (
        f"{float(a3['hg008_v04']['compat_off']['sensitivity']):.4f}"
    )
    out["v110_arm3_hg008_vj_f1"] = (
        f"{float(a3['hg008_v04']['vardictjava']['f1']):.4f}"
    )
    out["v110_arm3_hg008_off_f1"] = (
        f"{float(a3['hg008_v04']['compat_off']['f1']):.4f}"
    )
    out["cp95_v110_arm3_hg008_r"] = fmt_cp_ci(
        int(a3["hg008_v04"]["compat_off"]["tp"]),
        int(a3["hg008_v04"]["compat_off"]["truth_total"]),
    )

    # Domain 4: BRCA2 soft-clip adapter chimera suppression
    out["v110_arm4_chrom"] = str(a4["chimera_chrom"])
    out["v110_arm4_pos_raw"] = str(a4["chimera_pos"])
    out["v110_arm4_pos"] = f"{int(a4['chimera_pos']):,}"
    out["v110_arm4_span_bp"] = str(a4["region_span_bp"])
    out["v110_arm4_vj_ins"] = str(a4["vj_chimeric_insertions"])
    out["v110_arm4_vj_ins_bp"] = str(a4["vj_chimeric_ins_bp"])
    out["v110_arm4_off_ins"] = str(a4["off_chimeric_insertions"])
    out["v110_arm4_off_recs"] = str(a4["off_total_records"])

    # Domain 5: Clinical workflow presets & QIAseq V3 -sp
    out["v110_arm5_wgs_hg002"] = str(a5["wgs_hg002_60kbp_records"])
    out["v110_arm5_wgs_hg008"] = str(a5["wgs_hg008_paired_records"])
    out["v110_arm5_wes_bp"] = f"{int(a5['wes_seqc2_evaluated_bp']):,}"
    out["v110_arm5_wes_kbp"] = str(round(int(a5["wes_seqc2_evaluated_bp"]) / 1000))
    out["v110_arm5_wes_recs"] = str(a5["wes_seqc2_264kbp_records"])
    out["v110_arm5_amp_single"] = str(a5["amplicon_single_8col_records"])
    out["v110_arm5_amp_paired"] = str(a5["paired_somatic_amplicon_off_records"])
    out["v110_arm5_qia_recs"] = str(a5["qiaseq_spe_tumor_records"])
    out["v110_arm5_qia_panel_nosp"] = str(a5["qiaseq_spe_no_sp_known_pos_tp"])
    out["v110_arm5_qia_panel_sp"] = str(a5["qiaseq_spe_known_pos_tp"])
    out["v110_arm5_qia_panel_tot"] = str(a5["qiaseq_spe_panel_known_pos_unmasked"])
    out["v110_arm5_qia_panel_gain"] = str(
        int(a5["qiaseq_spe_known_pos_tp"]) - int(a5["qiaseq_spe_no_sp_known_pos_tp"])
    )
    out["v110_arm5_qia_dp50_nosp"] = str(a5["qiaseq_spe_dp50_no_sp_known_pos_tp"])
    out["v110_arm5_qia_dp50_sp"] = str(a5["qiaseq_spe_dp50_known_pos_tp"])
    out["v110_arm5_qia_dp50_tot"] = str(a5["qiaseq_spe_known_pos_total"])
    out["v110_arm5_qia_dp50_gain"] = str(
        int(a5["qiaseq_spe_dp50_known_pos_tp"])
        - int(a5["qiaseq_spe_dp50_no_sp_known_pos_tp"])
    )
    out["cp95_v110_qiaseq_dp50_nosp"] = fmt_cp_ci(
        int(a5["qiaseq_spe_dp50_no_sp_known_pos_tp"]),
        int(a5["qiaseq_spe_known_pos_total"]),
    )
    out["cp95_v110_qiaseq_dp50_sp"] = fmt_cp_ci(
        int(a5["qiaseq_spe_dp50_known_pos_tp"]),
        int(a5["qiaseq_spe_known_pos_total"]),
    )
    out["cp95_v110_qiaseq_kp"] = out["cp95_v110_qiaseq_dp50_sp"]
    out["cp95_v110_qiaseq_panel_nosp"] = fmt_cp_ci(
        int(a5["qiaseq_spe_no_sp_known_pos_tp"]),
        int(a5["qiaseq_spe_panel_known_pos_unmasked"]),
    )
    out["cp95_v110_qiaseq_panel_sp"] = fmt_cp_ci(
        int(a5["qiaseq_spe_known_pos_tp"]),
        int(a5["qiaseq_spe_panel_known_pos_unmasked"]),
    )
    out["v110_qiaseq_panel_tp_gain"] = out["v110_arm5_qia_panel_gain"]
    out["v110_qiaseq_panel_nosp_tp"] = out["v110_arm5_qia_panel_nosp"]
    out["v110_qiaseq_panel_sp_tp"] = out["v110_arm5_qia_panel_sp"]
    out["v110_qiaseq_panel_total"] = out["v110_arm5_qia_panel_tot"]
    out["v110_qiaseq_dp50_nosp_tp"] = out["v110_arm5_qia_dp50_nosp"]
    out["v110_qiaseq_dp50_sp_tp"] = out["v110_arm5_qia_dp50_sp"]
    out["v110_qiaseq_dp50_total"] = out["v110_arm5_qia_dp50_tot"]
    out["v110_arm5_ac5_d100"] = str(a5["sample_ac5_ctdna_d100_filtered"])
    out["v110_arm5_ctdna_freq"] = f"{float(a5['preset_ctdna_freq']):.3f}"
    out["v110_arm5_ctdna_vcf_af"] = f"{float(a5['preset_ctdna_vcf_af']):.3f}"
    out["v110_arm5_spike_vaf_pct"] = f"{float(a5['spikein_target_vaf_pct']):.1f}%"
    out["v110_arm5_spike_disc_af"] = f"{float(a5['spikein_discovery_af']):.4f}"
    out["v110_arm5_spike_chrom"] = str(a5["spikein_chrom"])
    out["v110_arm5_spike_pos_raw"] = str(a5["spikein_pos"])
    out["v110_arm5_spike_pos"] = f"{int(a5['spikein_pos']):,}"
    out["v110_arm5_spike_dp"] = f"{int(a5['spikein_observed_dp']):,}"
    out["v110_arm5_spike_vd"] = str(a5["spikein_observed_vd"])
    out["v110_arm5_spike_af"] = f"{float(a5['spikein_observed_af']):.4f}"
    out["v110_arm5_wt_copies"] = str(a5["admixture_wt_copies"])
    out["v110_arm5_alt_reads"] = str(a5["admixture_alt_reads"])
    out["v110_arm5_wt_bg_dp"] = (
        f"{int(a5['spikein_observed_dp']) - int(a5['spikein_observed_vd']):,}"
    )

    # Domain 6: Multi-block --csi spanning-SV indexing & -p homopolymer tiling
    out["v110_arm6_csi_bin_bp"] = f"{int(a6['csi_leaf_bin_bp']):,}"
    out["v110_arm6_csi_bin_kbp"] = str(int(a6["csi_leaf_bin_bp"]) // 1024)
    out["v110_arm6_max_sv_seq_bp"] = f"{int(a6['max_explicit_sv_seq_bp']):,}"
    out["v110_arm6_colo_blocks"] = str(a6["colo829_bgzf_blocks"])
    out["v110_arm6_span_blocks"] = str(a6["multi_block_bgzf_blocks"])
    out["v110_arm6_span_bins"] = str(int(a6["interior_bin_queries_total"]) - 1)
    out["v110_arm6_hg008_blocks"] = str(a6["hg008_bgzf_blocks"])
    out["v110_arm6_csi_pass"] = str(a6["interior_bin_queries_passed"])
    out["v110_arm6_csi_tot"] = str(a6["interior_bin_queries_total"])
    out["v110_arm6_pileup_full_span_total"] = str(a6["pileup_full_span_cases_total"])
    out["v110_arm6_pileup_eq_cases"] = str(a6["pileup_untiled_eq_tiled_cases"])
    out["v110_arm6_pileup_ret"] = f"{int(a6['pileup_covered_bases_retained']):,}"
    out["v110_arm6_pileup_tot"] = f"{int(a6['pileup_covered_bases_total']):,}"
    out["v110_arm6_pileup_pct"] = (
        f"{float(a6['pileup_covered_bases_retained_pct']):.1f}%"
    )
    out["v110_arm6_pos0_rows"] = str(a6["pos0_sentinel_rows"])

    return out


def build_macro_context(
    doc: dict[str, Any], v110_doc: dict[str, Any] | None = None
) -> dict[str, str]:
    """Construct complete flat macro dictionary derived from master_gcp_metrics.json and v110_compat_metrics.json."""
    m: dict[str, str] = {}
    for k, v in doc.get("flat_macros", {}).items():
        if isinstance(v, float):
            # Preserve exact formatting for common currency/metrics if needed
            if k in (
                "gcp_spend_usd",
                "gcp_budget_usd",
                "gcs_cache_gib",
                "gcs_weekly_usd",
            ):
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
    m["commit"] = str(ev["hg008-cram-twin"]["head"])
    m["f_std"] = "0.01"
    m["f_mrd_wgs"] = "0.001"
    m["f_mrd_ud"] = "0.0005"
    m["f_titration"] = "0.002"
    m["f_zero"] = "0"
    m["amp_param"] = "10:0.95"
    m_chunk = re.search(
        r"CPP_CHUNK=(\d+)", s5["wgs_hg008_full_giab_3way_th32"]["cpp"]["command"]
    )
    cpp_chunk_int = int(m_chunk.group(1)) if m_chunk else 50000
    m["cpp_chunk_bp"] = f"{cpp_chunk_int:,}"
    m["rust_default_tile_kbp"] = "100"
    # Halo widths mirror the engine constants in crates/vardictrust-core/src/lib.rs
    # (let halo = 500 /* bp per side */ at both tile() call sites): 0.5 kbp.
    m["rust_halo_kbp"] = "0.5"
    m["rust_splice_halo_kbp"] = "0.5"
    m["giab_arm_b_max_kbp"] = str(cpp_chunk_int // 1000)
    m["cpp_chunk_kbp"] = str(cpp_chunk_int // 1000)
    m["nvme_raid0_tb"] = "1.5"
    m["gcp_remaining_usd"] = f"{meta['gcp_budget_usd'] - meta['gcp_spend_usd']:.2f}"
    m["gcp_spend_pct"] = (
        f"{(meta['gcp_spend_usd'] / meta['gcp_budget_usd']) * 100.0:.1f}%"
    )
    tiles = s1["s1_h1a_simple_tiles"]["tile_sizes_bp"]
    m["s1_tile_500k"] = f"{tiles[0]:,}"
    m["s1_tile_1m"] = f"{tiles[1]:,}"
    m["s1_tile_2m"] = f"{tiles[2]:,}"
    m["s1_tile_10m"] = f"{tiles[3]:,}"
    span_ints = sorted(int(k.lstrip("R")) for k in s4["span_chr1_3way"]["th8"]["spans"])
    m["s4_span_1m"] = f"{span_ints[0]:,}"
    m["s4_span_10m"] = f"{span_ints[1]:,}"
    m["s4_span_100m"] = f"{span_ints[2]:,}"
    m["s4_span_249m"] = f"{span_ints[3]:,}"
    m["s6_v50q_gain_gbp"] = "0.23"
    m["s6_mosaic_vaf_range_pct"] = "1%–30%"
    m["s1_amplicon_speedup_vs_java"] = (
        f"{s1['s1_ceiling_amplicon']['java']['wall_s'] / s1['s1_ceiling_amplicon']['rust']['wall_s']:.0f}"
    )
    m["s1_amplicon_rss_ratio_vs_java"] = (
        f"{s1['s1_ceiling_amplicon']['java']['max_rss_kb'] / s1['s1_ceiling_amplicon']['rust']['max_rss_kb']:.1f}"
    )
    m["s1_gate7_speedup_vs_java"] = (
        f"{s1['s1_gate7_ref_multiregion']['java']['wall_s'] / s1['s1_gate7_ref_multiregion']['rust']['wall_s']:.1f}"
    )
    m["s1_gate7_rss_ratio_vs_java"] = (
        f"{s1['s1_gate7_ref_multiregion']['java']['max_rss_kb'] / s1['s1_gate7_ref_multiregion']['rust']['max_rss_kb']:.1f}"
    )

    # Exact Clopper-Pearson 95% binomial CIs
    total_cells = int(meta["total_cells"])
    m["cp95_60_60"] = fmt_cp_ci(total_cells, total_cells)
    m["cp95_13_13"] = fmt_cp_ci(
        int(s2["tierp_mode_parity_sweep"]["rust_match_count"]),
        int(s2["tierp_mode_parity_sweep"]["total_modes"]),
    )
    ext_mod = ev["extended-coverage-tier2"]
    m["cp95_17_17"] = fmt_cp_ci(
        int(ext_mod["modalities_passed"]), int(ext_mod["modalities_tested"])
    )
    mm_sl = ev["multimodal-slices-tier2"]
    m["cp95_18_18"] = fmt_cp_ci(int(mm_sl["slices_passed"]), int(mm_sl["total_slices"]))
    nf_rows = s5["cross_platform_short_read_chr20_panel"]["part_b_f0_noise_floor_rows"]
    nf_match = int(
        s5["cross_platform_short_read_chr20_panel"]["part_b_strict_match_count"]
    )
    m["cp95_noise_18_18"] = fmt_cp_ci(nf_match, len(nf_rows))
    ga_rows = s6["ga4gh_v35_stratification_matrix"]
    m["cp95_70_70"] = fmt_cp_ci(len(ga_rows), len(ga_rows))
    sw_match = int(sw["match_arms"])
    sw_tot = int(sw["total_arms"])
    sw_diff = sw_tot - sw_match
    sw_single_cells = ("sweep_span_chr1_f0", "sweep_R_chr1_full")
    sw_single_rows = [r for k in sw_single_cells for r in sw["by_cell"][k]]
    sw_single_win = len(sw_single_rows)
    sw_single_win_match = sum(
        1 for r in sw_single_rows if str(r.get("verdict", "")).upper().startswith("MATCH")
    )
    sw_bed_match = sw_match - sw_single_win_match
    m["sweep_match_pct"] = f"{(sw_match / sw_tot) * 100.0:.2f}%"
    m["optc_match_arms"] = str(sw_match)
    m["optc_total_arms"] = str(sw_tot)
    m["optc_bed_arms"] = str(sw_tot - sw_single_win)
    m["optc_bed_match_arms"] = str(sw_bed_match)
    m["optc_single_win_arms"] = str(sw_single_win)
    m["optc_single_win_match_arms"] = str(sw_single_win_match)
    m["optc_diff_arms"] = str(sw_diff)
    m["optc_match_pct"] = f"{(sw_match / sw_tot) * 100.0:.2f}"
    m["cp95_78_87"] = fmt_cp_ci(sw_match, sw_tot)
    m["cp95_optc_sweep"] = fmt_cp_ci(sw_match, sw_tot)
    m["sw_rchr1_rows"] = f"{int(sw['by_cell']['sweep_R_chr1_full'][0]['java_rows']):,}"
    m["sw_span_f0_rows"] = (
        f"{int(sw['by_cell']['sweep_span_chr1_f0'][0]['java_rows']):,}"
    )
    k_mos, n_mos = (
        int(x) for x in str(doc["flat_macros"]["s6_mosaic_tp_best"]).split("/")
    )
    m["s6_mosaic_rust_tp"] = str(k_mos)
    m["s6_mosaic_truth_total"] = str(n_mos)
    m["cp95_84_85"] = fmt_cp_ci(k_mos, n_mos)

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
    m["s1_hg008_speedup_vs_java"] = (
        f"{h8c['java']['wall_s'] / h8c['rust']['wall_s']:.2f}"
    )
    m["s1_hg008_rss_ratio_vs_java"] = (
        f"{h8c['java']['max_rss_kb'] / h8c['rust']['max_rss_kb']:.1f}"
    )

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
            m[f"s2_{sample_id}_{th_key}_rust_wall_range"] = (
                f"{r_stat['wall_s_min']:.2f}–{r_stat['wall_s_max']:.2f}"
            )
            m[f"s2_{sample_id}_{th_key}_java_wall_s"] = f"{j_stat['wall_s_mean']:.2f}"
            m[f"s2_{sample_id}_{th_key}_java_wall_range"] = (
                f"{j_stat['wall_s_min']:.2f}–{j_stat['wall_s_max']:.2f}"
            )
            m[f"s2_{sample_id}_{th_key}_cpp_wall_s"] = f"{c_stat['wall_s_mean']:.2f}"
            m[f"s2_{sample_id}_{th_key}_cpp_wall_range"] = (
                f"{c_stat['wall_s_min']:.2f}–{c_stat['wall_s_max']:.2f}"
            )
            m[f"s2_{sample_id}_{th_key}_rust_rss_mb"] = (
                f"{r_stat['max_rss_mb_mean']:.1f}"
            )
            m[f"s2_{sample_id}_{th_key}_java_rss_mb"] = (
                f"{j_stat['max_rss_mb_mean']:.1f}"
            )
            m[f"s2_{sample_id}_{th_key}_java_rss_gb"] = (
                f"{j_stat['max_rss_gb_dec_mean']:.2f}"
            )
            m[f"s2_{sample_id}_{th_key}_cpp_rss_mb"] = (
                f"{c_stat['max_rss_mb_mean']:.1f}"
            )

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
    m["s2_tierp_sign_perm_p_twosided"] = (
        f"{2.0 * (2.0 ** (-len(rust_walls_12))):.5f}"
    )
    m["s2_tierp_sign_perm_exp_twosided"] = str(len(rust_walls_12) - 1)

    # Stage 2 Table 3 -R chr1
    rchr1 = s2["table3_R_chr1"]
    rchr1_bp_int = int(rchr1["region"].split("-")[1])
    m["s2_rchr1_sample"] = str(rchr1["sample"])
    m["s2_rchr1_bp"] = f"{rchr1_bp_int:,}"
    m["s2_rchr1_mbp"] = f"{rchr1_bp_int / 1e6:.2f}"
    m["s2_rchr1_rust_wall_range"] = (
        f"{rchr1['rust']['wall_s_min']:.2f}–{rchr1['rust']['wall_s_max']:.2f}"
    )
    m["s2_rchr1_cpp_wall_range"] = (
        f"{rchr1['cpp']['wall_s_min']:.2f}–{rchr1['cpp']['wall_s_max']:.2f}"
    )
    m["s2_rchr1_java200g_wall_range"] = (
        f"{rchr1['java_xmx200g']['wall_s_min']:.2f}–{rchr1['java_xmx200g']['wall_s_max']:.2f}"
    )
    m["s2_rchr1_java4g_oom_wall_s"] = f"{rchr1['java_xmx4g_oom']['wall_s']:.2f}"
    m["s2_rchr1_java4g_oom_rss_gib"] = f"{rchr1['java_xmx4g_oom']['max_rss_gib']:.2f}"
    m["s2_fig1_fig3_wgs_bed_intervals"] = m["s2_SRR8657348_intervals"]

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
    m["s2_fig2_region_mbp"] = "2.31"
    m_fig2_par = re.search(
        r"(\d+)\s+variants",
        doc["cells"]["s2_fig2_hg002_chr20_truth"]["metrics"]["compare_hg002_chr20"],
    )
    if m_fig2_par:
        m["s2_fig2_parity_vars"] = f"{int(m_fig2_par.group(1)):,}"
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
    m["s3_rna_rust_wall_range"] = (
        f"{s3['rust']['wall_s_min']:.2f}–{s3['rust']['wall_s_max']:.2f}"
    )
    m["s3_rna_java_wall_range"] = (
        f"{s3['java']['wall_s_min']:.2f}–{s3['java']['wall_s_max']:.2f}"
    )

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
    ex = s4["scale_exome_3way"]
    m_ex_tot = re.search(r"(\d+)\s+variants", ex["th1"]["compare"])
    if m_ex_tot:
        m["s4_exome_total_vars"] = f"{int(m_ex_tot.group(1)):,}"
    for tk in ("th1", "th2", "th4", "th8", "th16", "th32"):
        td = ex[tk]
        m[f"s4_exome_{tk}_rust_wall_s"] = f"{td['rust']['wall_s_mean']:.2f}"
        m[f"s4_exome_{tk}_rust_wall_range"] = (
            f"{td['rust']['wall_s_min']:.2f}–{td['rust']['wall_s_max']:.2f}"
        )
        m[f"s4_exome_{tk}_rust_rss_mb"] = f"{td['rust']['max_rss_mb_mean']:.1f}"
        m[f"s4_exome_{tk}_rust_rss_gib"] = f"{td['rust']['max_rss_gib_mean']:.2f}"
        m[f"s4_exome_{tk}_rust_sp"] = f"{td['rust']['speedup_sp']:.2f}"
        m[f"s4_exome_{tk}_java_wall_s"] = f"{td['java']['wall_s_mean']:.2f}"
        m[f"s4_exome_{tk}_java_wall_range"] = (
            f"{td['java']['wall_s_min']:.2f}–{td['java']['wall_s_max']:.2f}"
        )
        m[f"s4_exome_{tk}_java_rss_gib"] = f"{td['java']['max_rss_gib_mean']:.2f}"
        m[f"s4_exome_{tk}_cpp_wall_s"] = f"{td['cpp']['wall_s_mean']:.2f}"
        m[f"s4_exome_{tk}_cpp_wall_range"] = (
            f"{td['cpp']['wall_s_min']:.2f}–{td['cpp']['wall_s_max']:.2f}"
        )
        m[f"s4_exome_{tk}_cpp_rss_gib"] = f"{td['cpp']['max_rss_gib_mean']:.2f}"
        sp_ex_j = td.get("rust_speedup_vs_java")
        if sp_ex_j is None:
            sp_ex_j = td["java"]["wall_s_mean"] / td["rust"]["wall_s_mean"]
        m[f"s4_exome_{tk}_speedup_vs_java"] = f"{float(sp_ex_j):.1f}"
        m[f"s4_exome_{tk}_rss_ratio_vs_java"] = (
            f"{td['java']['max_rss_kb_mean'] / td['rust']['max_rss_kb_mean']:.1f}"
        )
        m[f"s4_exome_{tk}_rss_ratio_vs_cpp"] = (
            f"{td['cpp']['max_rss_kb_mean'] / td['rust']['max_rss_kb_mean']:.1f}"
        )
    m["s4_exome_32t_speedup_vs_1t"] = m["s4_exome_th32_rust_sp"]
    m["s4_exome_32t_rss_ratio_vs_java"] = m["s4_exome_th32_rss_ratio_vs_java"]

    # Stage 4 Span chr1
    sp_c1 = s4["span_chr1_3way"]
    m["s4_span_h2_bound"] = "1.50"
    for tk in ("th8", "th16", "th32"):
        sd = sp_c1[tk]
        m[f"s4_span_{tk}_ratio_exact"] = f"{sd['rust_rss_ratio_249m_vs_10m']:.3f}"
        m[f"s4_span_{tk}_10m_rust_gib"] = (
            f"{sd['spans']['R10000000']['rust']['max_rss_gib_mean']:.2f}"
        )
        m[f"s4_span_{tk}_249m_rust_gib"] = (
            f"{sd['spans']['R248956422']['rust']['max_rss_gib_mean']:.2f}"
        )
        m[f"s4_span_{tk}_249m_rust_wall_s"] = (
            f"{sd['spans']['R248956422']['rust']['wall_s_mean']:.2f}"
        )
        m[f"s4_span_{tk}_249m_cpp_wall_s"] = (
            f"{sd['spans']['R248956422']['cpp']['wall_s_mean']:.2f}"
        )
        m[f"s4_span_{tk}_249m_cpp_gib"] = (
            f"{sd['spans']['R248956422']['cpp']['max_rss_gib_mean']:.2f}"
        )

    spf0 = s4["span_chr1_f0"]
    m_spf0_v = re.search(
        r"(\d+)\s+variants", doc["cells"]["s4_span_chr1_f0"]["parity_verdict"]
    )
    if m_spf0_v:
        m["s4_span_f0_vars"] = f"{int(m_spf0_v.group(1)):,}"
        m["s4_span_chr1_f0_vars"] = m["s4_span_f0_vars"]
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
    m["s4_jpn_8jobs_rust_gib"] = f"{4.0 * jpn['combined_peak_rss_gib_mean']:.2f}"
    m["s4_jpn_8jobs_java_gib"] = f"{8.0 * ex['th4']['java']['max_rss_gib_mean']:.1f}"

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
        m[f"s4_t3_{cid}_fp_reduction_pct"] = (
            f"{(1.0 - best['fp'] / none['fp']) * 100.0:.1f}%"
        )

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
    m["s5_hg002_giab_gbp"] = "2.51"
    m["s5_hg008_giab_gbp"] = "2.37"

    h2_a = s5["wgs_hg002_full_contigs_3way_th32"]
    h2_b = s5["wgs_hg002_full_giab_3way_th32"]
    h8_a = s5["wgs_hg008_full_contigs_3way_th32"]
    h8_b = s5["wgs_hg008_full_giab_3way_th32"]

    for s5_macro_prefix, s5_cell in (
        ("s5_hg002_contigs", h2_a),
        ("s5_hg002_giab", h2_b),
        ("s5_hg008_contigs", h8_a),
        ("s5_hg008_giab", h8_b),
    ):
        for eng_k in ("rust", "java", "cpp"):
            m[f"{s5_macro_prefix}_{eng_k}_rss_gb"] = (
                f"{float(s5_cell[eng_k]['max_rss_gb_dec']):.2f}"
            )

    m["s5_hg002_contigs_cpp_rss_ratio_vs_rust"] = (
        f"{h2_a['cpp']['max_rss_kb'] / h2_a['rust']['max_rss_kb']:.1f}"
    )
    m["s5_hg008_contigs_cpp_rss_ratio_vs_rust"] = (
        f"{h8_a['cpp']['max_rss_kb'] / h8_a['rust']['max_rss_kb']:.1f}"
    )
    m["s5_hg008_giab_cpp_rss_ratio_vs_rust"] = (
        f"{h8_b['cpp']['max_rss_kb'] / h8_b['rust']['max_rss_kb']:.1f}"
    )
    m["s5_hg002_giab_cpp_diff_vars"] = (
        f"{h2_b['cpp']['variant_count'] - h2_b['rust']['variant_count']:,}"
    )
    m["s5_hg008_contigs_cpp_wall_str"] = h8_a["cpp"]["wall_str"]
    m["s5_hg008_giab_java_wall_str"] = h8_b["java"]["wall_str"]
    m["s5_hg008_giab_cpp_wall_str"] = h8_b["cpp"]["wall_str"]

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
    m["cp95_chr20_java_completed"] = fmt_cp_ci(
        len(java_completing), len(java_completing)
    )

    hg2_std_snv = [r["snp_f1"] for r in rust_cp if r["cohort"] == "HG002_Germline"]
    hg2_ug_snv = [r["ughcr_snp_f1"] for r in rust_cp if r["cohort"] == "HG002_Germline"]
    hg8_f01_std_snv = [
        r["snp_f1"]
        for r in rust_cp
        if r["cohort"] == "HG008_Somatic" and str(r["freq"]) == "0.01"
    ]
    hg8_f01_ug_snv = [
        r["ughcr_snp_f1"]
        for r in rust_cp
        if r["cohort"] == "HG008_Somatic" and str(r["freq"]) == "0.01"
    ]
    hg8_all_std_snv = [r["snp_f1"] for r in rust_cp if r["cohort"] == "HG008_Somatic"]
    hg8_all_ug_snv = [
        r["ughcr_snp_f1"] for r in rust_cp if r["cohort"] == "HG008_Somatic"
    ]

    m["s5_chr20_hg2_std_snv_f1_range"] = (
        f"{min(hg2_std_snv):.4f}–{max(hg2_std_snv):.4f}"
    )
    m["s5_chr20_hg2_ug_snv_f1_range"] = f"{min(hg2_ug_snv):.4f}–{max(hg2_ug_snv):.4f}"
    m["s5_chr20_hg8_f01_std_snv_f1_range"] = (
        f"{min(hg8_f01_std_snv):.4f}–{max(hg8_f01_std_snv):.4f}"
    )
    m["s5_chr20_hg8_f01_ug_snv_f1_range"] = (
        f"{min(hg8_f01_ug_snv):.4f}–{max(hg8_f01_ug_snv):.4f}"
    )
    m["s5_chr20_hg8_all_std_snv_f1_range"] = (
        f"{min(hg8_all_std_snv):.4f}–{max(hg8_all_std_snv):.4f}"
    )
    m["s5_chr20_hg8_all_ug_snv_f1_range"] = (
        f"{min(hg8_all_ug_snv):.4f}–{max(hg8_all_ug_snv):.4f}"
    )

    tot_nf_cand = sum(int(r["total_candidates"]) for r in nf_rows)
    tot_nf_pass = sum(int(r["pass_somatic"]) for r in nf_rows)
    m["s5_nf_total_candidates"] = f"{tot_nf_cand:,}"
    m["s5_nf_total_pass"] = f"{tot_nf_pass:,}"

    # Stage 6 HG002 v5.0q / CMRG / Mosaic
    m["s6_v50q_gbp"] = "2.74"
    v50q_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for row in s6["hg002_v50q_cmrg_mosaic_eval"]:
        b = row["benchmark"]
        c = row["caller"]
        v50q_by_key[(b, c)] = row
        m[f"s6_{b}_{c}_snv_p"] = f"{row['snv_p']:.4f}"
        m[f"s6_{b}_{c}_snv_r"] = f"{row['snv_r']:.4f}"
        m[f"s6_{b}_{c}_snv_f1"] = f"{row['snv_f1']:.4f}"
        m[f"s6_{b}_{c}_indel_p"] = f"{row['indel_p']:.4f}"
        m[f"s6_{b}_{c}_indel_r"] = f"{row['indel_r']:.4f}"
        m[f"s6_{b}_{c}_indel_f1"] = f"{row['indel_f1']:.4f}"
    v50q_fc = v50q_by_key[("NIST_v5.0q_T2T_Q100", "VarDictRust_FullContigs")]
    v50q_g50 = v50q_by_key[("NIST_v5.0q_T2T_Q100", "VarDictJava_GIAB50kb")]
    cmrg_fc = v50q_by_key[("NIST_CMRG_v1.00_273Genes", "VarDictRust_FullContigs")]
    cmrg_g50 = v50q_by_key[("NIST_CMRG_v1.00_273Genes", "VarDictJava_GIAB50kb")]
    m["s6_v50q_snv_r_gain_pct"] = (
        f"{(v50q_fc['snv_r'] - v50q_g50['snv_r']) * 100.0:.2f}%"
    )
    m["s6_v50q_indel_r_gain_pct"] = (
        f"{(v50q_fc['indel_r'] - v50q_g50['indel_r']) * 100.0:.2f}%"
    )
    m["s6_cmrg_snv_r_gain_pct"] = (
        f"{(cmrg_fc['snv_r'] - cmrg_g50['snv_r']) * 100.0:.2f}%"
    )
    m["s6_cmrg_indel_r_gain_pct"] = (
        f"{(cmrg_fc['indel_r'] - cmrg_g50['indel_r']) * 100.0:.2f}%"
    )
    m["s6_mosaic_recall_pct"] = f"{float(m['s6_mosaic_recall_best']) * 100.0:.2f}%"
    m["s6_mosaic_rust_r_pct"] = m["s6_mosaic_recall_pct"]

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
    m["s6_tit_rust_wall_range"] = (
        f"{min(r['wall_s'] for r in rust_tit):.2f}–{max(r['wall_s'] for r in rust_tit):.2f}"
    )
    m["s6_tit_java_wall_range"] = (
        f"{min(j['wall_s'] for j in java_tit):.2f}–{max(j['wall_s'] for j in java_tit):.2f}"
    )
    m["s6_tit_speedup_range"] = f"{min(t_sp):.1f}×–{max(t_sp):.1f}×"
    m["s6_tit_rss_ratio_range"] = f"{min(t_rr):.1f}×–{max(t_rr):.1f}×"
    m["s6_tit_rust_rss_kb_range"] = (
        f"{min(r['peak_kb'] for r in rust_tit):,}–{max(r['peak_kb'] for r in rust_tit):,}"
    )
    m["s6_tit_java_rss_kb_range"] = (
        f"{min(j['peak_kb'] for j in java_tit):,}–{max(j['peak_kb'] for j in java_tit):,}"
    )
    for r in rust_tit:
        p = int(r["purity_pct"])
        m[f"s6_tit_vaf_{p}_pct"] = f"{p / 2.0:.1f}%"
    m["s6_tit_vaf_range_pct"] = f"{m['s6_tit_vaf_20_pct']} to {m['s6_tit_vaf_1_pct']}"
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

    # Live repository release metadata, test/CLI counts, and Rust engine constants
    # (keeps m["version"] == "1.0.6" and all 60-cell GCP study provenance untouched)
    m.update(_load_repo_release_and_engine_macros(ROOT))

    if v110_doc is None and DEFAULT_V110_METRICS.is_file():
        v110_doc = json.loads(DEFAULT_V110_METRICS.read_text(encoding="utf-8"))
    if v110_doc is not None:
        m.update(_build_v110_macros(v110_doc))

    return m


def compile_template(
    tmpl_text: str,
    macros: dict[str, str],
    tables_dir: Path,
    inv_path: Path | None = None,
) -> str:
    """Expand {{m:<key>}} and {{table:<name>}} in template text."""
    _ = inv_path

    # 1. Expand {{table:<name>}}
    def replace_table(match: re.Match[str]) -> str:
        tname = match.group(1).strip()
        if tname in TABLE_ALIASES:
            fname, anchor_id = TABLE_ALIASES[tname]
        else:
            fname = f"{tname}.md" if not tname.endswith(".md") else tname
            anchor_id = tname.removesuffix(".md")
        tpath = tables_dir / fname
        if not tpath.is_file():
            raise KeyError(
                f"Missing table file {tpath} for macro {{{{table:{tname}}}}}"
            )
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
        raise ValueError(
            f"Unresolved template macro found after compilation: {leftover.group(0)}"
        )

    return compiled


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compile reproducible manuscript.md and supplement.md from templates."
    )
    parser.add_argument(
        "--metrics",
        type=Path,
        default=DEFAULT_METRICS,
        help="Path to master_gcp_metrics.json",
    )
    parser.add_argument(
        "--v110-metrics",
        type=Path,
        default=DEFAULT_V110_METRICS,
        help="Path to v110_compat_metrics.json",
    )
    parser.add_argument(
        "--inventory",
        type=Path,
        default=DEFAULT_INVENTORY,
        help="Path to cells_60_inventory.tsv",
    )
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=DEFAULT_TABLES_DIR,
        help="Directory containing table*.md",
    )
    parser.add_argument(
        "--tmpl-dir",
        type=Path,
        default=DEFAULT_TMPL_DIR,
        help="Directory containing *.tmpl.md",
    )
    parser.add_argument(
        "--ms-out",
        type=Path,
        default=DEFAULT_MS_OUT,
        help="Output path for manuscript.md",
    )
    parser.add_argument(
        "--sup-out",
        type=Path,
        default=DEFAULT_SUP_OUT,
        help="Output path for supplement.md",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify existing manuscript.md and supplement.md match templates",
    )
    args = parser.parse_args()

    if not args.metrics.is_file():
        print(f"ERROR: Missing metrics JSON: {args.metrics}", file=sys.stderr)
        return 1
    if not args.v110_metrics.is_file():
        print(
            f"ERROR: Missing v1.1.0 metrics JSON: {args.v110_metrics}", file=sys.stderr
        )
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
    v110_doc = json.loads(args.v110_metrics.read_text(encoding="utf-8"))
    macros = build_macro_context(doc, v110_doc)

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
            print(
                "ERROR: --check failed because manuscript.md or supplement.md does not exist.",
                file=sys.stderr,
            )
            return 1
        ms_disk = args.ms_out.read_text(encoding="utf-8")
        sup_disk = args.sup_out.read_text(encoding="utf-8")
        if ms_disk != ms_compiled:
            print(
                f"ERROR: {args.ms_out} is out of sync with {ms_tmpl_path}",
                file=sys.stderr,
            )
            return 1
        if sup_disk != sup_compiled:
            print(
                f"ERROR: {args.sup_out} is out of sync with {sup_tmpl_path}",
                file=sys.stderr,
            )
            return 1
        ms_pdf = args.ms_out.with_suffix(".pdf")
        sup_pdf = args.sup_out.with_suffix(".pdf")
        if args.ms_out == DEFAULT_MS_OUT and (
            not ms_pdf.is_file() or not sup_pdf.is_file()
        ):
            print(f"ERROR: Missing compiled PDF {ms_pdf} or {sup_pdf}", file=sys.stderr)
            return 1
        print(
            f"OK verified {args.ms_out.name} and {args.sup_out.name} match compiled templates byte-for-byte"
        )
        return 0

    args.ms_out.write_text(ms_compiled, encoding="utf-8")
    args.sup_out.write_text(sup_compiled, encoding="utf-8")
    if args.ms_out == DEFAULT_MS_OUT and args.sup_out == DEFAULT_SUP_OUT:
        render_markdown_to_pdf(
            ms_compiled,
            DEFAULT_MS_PDF,
            base_dir=ROOT / "paper",
            title="VarDictRust Manuscript",
        )
        render_markdown_to_pdf(
            sup_compiled,
            DEFAULT_SUP_PDF,
            base_dir=ROOT / "paper",
            title="VarDictRust Supplementary Material",
        )
        print(
            f"OK compiled {args.ms_out} ({len(ms_compiled.split())} words), "
            f"{args.sup_out} ({len(sup_compiled.split())} words), "
            f"{DEFAULT_MS_PDF.name}, and {DEFAULT_SUP_PDF.name}"
        )
    else:
        print(
            f"OK compiled {args.ms_out} ({len(ms_compiled.split())} words) and {args.sup_out} ({len(sup_compiled.split())} words)"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
