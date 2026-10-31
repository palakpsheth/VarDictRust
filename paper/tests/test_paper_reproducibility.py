#!/usr/bin/env python3
"""End-to-end reproducibility and referee/biostatistics compliance tests for paper/manuscript.md and paper/supplement.md."""

from __future__ import annotations

import ast
import csv
import re
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD_SCRIPT = ROOT / "paper" / "build_manuscript.py"
MANUSCRIPT_TMPL = ROOT / "paper" / "templates" / "manuscript.tmpl.md"
SUPPLEMENT_TMPL = ROOT / "paper" / "templates" / "supplement.tmpl.md"
MANUSCRIPT = ROOT / "paper" / "manuscript.md"
SUPPLEMENT = ROOT / "paper" / "supplement.md"
INV_PATH = ROOT / "benchmarks" / "gcp-study-v106" / "cells_60_inventory.tsv"

FORBIDDEN_TEMPLATE_EMPIRICAL_STRINGS = {
    "0.9858",
    "0.8956",
    "0.9699",
    "0.9651",
    "0.6622",
    "0.6624",
    "0.9367",
    "0.8920",
    "0.9882",
    "812.75",
    "501.83",
    "339,698",
    "8,111,171",
    "9,643,769",
    "13,067,287",
    "10,195,118",
    "309,569",
    "213.70",
    "106.47",
    "228.90",
    "760.96",
    "12,668",
    "54,448",
}

EXPECTED_FIGURE_IDS = [
    "fig1_architecture_pipeline",
    "fig2_parity_and_tierp_replication",
    "fig3_truth_qual_and_titration",
    "fig4_memory_and_multicore_scaling",
    "fig5_full_wgs_arm_a_vs_arm_b",
    "fig6_cross_platform_and_ga4gh",
    "fig7_nextgen_truth_v50q_cmrg_mosaic",
    "fig8_modern_callers_and_low_vaf",
]

EXPECTED_TABLE_IDS = [
    "table1_parity_matrix",
    "table2_tierp_and_rchr1",
    "table3_truth_and_field_callers",
    "table4_memory_scaling_and_packing",
    "table5_full_wgs_arm_a_vs_b",
    "table6_cross_platform_chr20",
    "table7_v50q_cmrg_mosaic",
    "table8_hg008_modern_callers_and_titration",
]

EXPECTED_SUPP_TABLE_IDS = [
    "supp_table1_60cell_inventory",
    "supp_table2_option_c_87arms",
    "supp_table3_n3_technical_replicates",
    "supp_table4_qual_roc_and_noise_floor",
    "supp_table5_ga4gh_70cell_matrix",
    "supp_table6_sv_and_amplicon_splice",
]


class TestPaperReproducibility(unittest.TestCase):
    def test_manuscript_and_supplement_are_script_reproducible_and_complete(self) -> None:
        self.assertTrue(BUILD_SCRIPT.is_file(), f"Missing {BUILD_SCRIPT}")
        res = subprocess.run(
            [sys.executable, str(BUILD_SCRIPT), "--check"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=f"stderr={res.stderr}\nstdout={res.stdout}")

        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")
        combined = ms_text + "\n" + sup_text

        # 1. Zero unresolved template macros or NaN/TODO literals
        self.assertIsNone(re.search(r"\{\{[^}]+\}\}", combined), "Unresolved {{...}} macro found")
        self.assertIsNone(re.search(r"\b(nan|NaN|TODO|TBD|FIXME)\b", combined), "Forbidden literal found")

        # 2. All 60 directory entries appear in manuscript + supplement
        with INV_PATH.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 60)
        for r in rows:
            cid = r["cell_id"]
            self.assertIn(cid, combined, f"Cell {cid} missing from manuscript+supplement")

        # 3. All 8 Figures, 8 Main Tables, and 6 Supplementary Tables referenced
        for i in range(1, 9):
            self.assertIn(f"Figure {i}", ms_text)
            self.assertIn(f"Table {i}", ms_text)
        for i in range(1, 7):
            self.assertIn(f"Table S{i}", combined)

        # 4. Required Genome Biology / Bioinformatics sections & disclosures present
        for req in (
            "### Background",
            "### Results",
            "### Conclusions",
            "## Availability of data and materials",
            "5865a2a0b24a8c17",
            "c490c36317b8f25b",
            "--nosv",
            "Liss_lab",
        ):
            self.assertIn(req, ms_text, f"Missing required token/section: {req}")

    def test_templates_contain_no_hardcoded_empirical_numbers(self) -> None:
        self.assertTrue(MANUSCRIPT_TMPL.is_file(), f"Missing {MANUSCRIPT_TMPL}")
        self.assertTrue(SUPPLEMENT_TMPL.is_file(), f"Missing {SUPPLEMENT_TMPL}")

        for tmpl_path in (MANUSCRIPT_TMPL, SUPPLEMENT_TMPL):
            raw = tmpl_path.read_text(encoding="utf-8")
            # Strip {{m:...}}, {{table:...}}, and {{cells_60_ledger}} macros
            stripped = re.sub(r"\{\{[^}]+\}\}", "MACRO_TOKEN", raw)

            for bad_lit in FORBIDDEN_TEMPLATE_EMPIRICAL_STRINGS:
                self.assertNotIn(
                    bad_lit,
                    stripped,
                    f"Hardcoded empirical string {bad_lit!r} found outside macros in {tmpl_path.name}",
                )

            # Check for hardcoded empirical floats/ratios/percentages/comma-counts/timings/RSS outside macros
            patterns = [
                (r"\b\d+\.\d+[x×]\b", "hardcoded decimal speedup/ratio"),
                (r"\b\d+\.\d+%", "hardcoded decimal percentage"),
                (r"\b\d{1,3},\d{3}(?:,\d{3})?\b", "hardcoded comma-grouped count"),
                (r"\b\d+\.\d+s\b", "hardcoded decimal wall-time seconds"),
                (r"\b\d+\.\d+\s*(?:KB|MB|MiB|GB|GiB|TB)\b", "hardcoded decimal memory size"),
                (r"\b0\.\d{4}\b", "hardcoded 4-decimal accuracy metric"),
            ]
            for pat, desc in patterns:
                m = re.search(pat, stripped)
                self.assertIsNone(
                    m,
                    f"Found {desc} ({m.group(0) if m else ''!r}) outside {{m:...}} macros in {tmpl_path.name}",
                )

    def test_compiled_manuscript_and_supplement_have_zero_unresolved_macros(self) -> None:
        self.assertTrue(MANUSCRIPT.is_file(), f"Missing {MANUSCRIPT}")
        self.assertTrue(SUPPLEMENT.is_file(), f"Missing {SUPPLEMENT}")
        for doc_path in (MANUSCRIPT, SUPPLEMENT):
            txt = doc_path.read_text(encoding="utf-8")
            self.assertNotIn("{{", txt, f"Unresolved '{{{{' in {doc_path.name}")
            self.assertNotIn("}}", txt, f"Unresolved '}}}}' in {doc_path.name}")
            m = re.search(r"\b(NaN|nan|null|None|TODO|TBD|FIXME)\b", txt)
            self.assertIsNone(
                m,
                f"Forbidden token {m.group(0) if m else ''!r} found in {doc_path.name}",
            )

    def test_all_8_figures_and_14_tables_referenced(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")
        combined = ms_text + "\n" + sup_text

        for i in range(1, 9):
            self.assertIn(f"Figure {i}", ms_text)
            self.assertIn(f"figures/fig{i}.svg", ms_text)
            self.assertIn(f"Table {i}", ms_text)
        for i in range(1, 7):
            self.assertIn(f"Table S{i}", combined)

        for fid in EXPECTED_FIGURE_IDS:
            self.assertIn(fid, combined, f"Missing figure identifier {fid}")
        for tid in EXPECTED_TABLE_IDS:
            self.assertIn(tid, combined, f"Missing main table identifier {tid}")
        for stid in EXPECTED_SUPP_TABLE_IDS:
            self.assertIn(stid, combined, f"Missing supplementary table identifier {stid}")

    def test_all_60_cells_referenced_across_manuscript_and_supplement(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")
        combined = ms_text + "\n" + sup_text

        with INV_PATH.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 60)
        for r in rows:
            cid = r["cell_id"]
            self.assertIn(cid, combined, f"Cell {cid} missing from manuscript+supplement")
            self.assertIn(cid, sup_text, f"Cell {cid} missing from supplement.md")

    def test_referee_and_biostatistics_compliance(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")
        combined = ms_text + "\n" + sup_text

        required_phrases = [
            "Clopper-Pearson",
            "technical replicates",
            "Mann-Whitney",
            "--fix",
            "s2_fig1_fig3_field_callers",
            "rc=124",
            "rc=3",
            "--nosv",
            "Liss_lab",
            "in silico",
            "c3d-standard-32",
            "2a8e9a3",
            "```bash",
        ]
        for phrase in required_phrases:
            self.assertIn(phrase, combined, f"Missing required referee/biostatistics disclosure: {phrase!r}")

    def test_build_manuscript_script_has_no_hardcoded_empirical_strings(self) -> None:
        self.assertTrue(BUILD_SCRIPT.is_file(), f"Missing {BUILD_SCRIPT}")
        tree = ast.parse(BUILD_SCRIPT.read_text(encoding="utf-8"), filename=str(BUILD_SCRIPT))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                val = node.value
                self.assertIsNone(
                    re.search(r"\b\d{1,3},\d{3}(?:,\d{3})?\b", val),
                    f"Hardcoded comma-grouped count {val!r} at {BUILD_SCRIPT.name}:{getattr(node, 'lineno', '?')}",
                )
                self.assertIsNone(
                    re.search(r"\b\d+\.\d+%", val),
                    f"Hardcoded percentage string {val!r} at {BUILD_SCRIPT.name}:{getattr(node, 'lineno', '?')}",
                )

    def test_compiler_check_mode_is_idempotent(self) -> None:
        self.assertTrue(BUILD_SCRIPT.is_file(), f"Missing {BUILD_SCRIPT}")
        res1 = subprocess.run(
            [sys.executable, str(BUILD_SCRIPT), "--check"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        self.assertEqual(res1.returncode, 0, msg=f"stderr={res1.stderr}\nstdout={res1.stdout}")
        res2 = subprocess.run(
            [sys.executable, str(BUILD_SCRIPT), "--check"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        self.assertEqual(res2.returncode, 0, msg=f"stderr={res2.stderr}\nstdout={res2.stdout}")


if __name__ == "__main__":
    unittest.main()
