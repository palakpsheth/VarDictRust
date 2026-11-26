#!/usr/bin/env python3
"""Unit tests for paper/build_figures_and_tables.py (offline reproducibility & zero hardcoded metrics)."""

from __future__ import annotations

import ast
import csv
import json
import re
import subprocess
import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "paper" / "build_figures_and_tables.py"
FIG_DIR = ROOT / "paper" / "figures"
TAB_DIR = ROOT / "paper" / "tables"
INV_PATH = ROOT / "benchmarks" / "gcp-study-v106" / "cells_60_inventory.tsv"
V110_PATH = ROOT / "benchmarks" / "v110-compat-study" / "v110_compat_metrics.json"

ALLOWED_INT_CONSTANTS = (
    set(range(-20, 101))
    | {
        120,
        124,
        128,
        134,
        137,
        150,
        160,
        180,
        200,
        226,
        255,
        300,
        360,
        400,
        500,
        600,
        1000,
        1024,
        2000,
        3600,
        10000,
        50000,
        65536,
        100000,
        500000,
        1000000,
        1048576,
        2000000,
        10000000,
        100000000,
        248956422,
        1700000000,
    }
)

ALLOWED_FLOAT_MULTIPLIERS = {
    1e-9,
    1e-6,
    0.001,
    0.002,
    0.005,
    0.025,
    0.125,
    0.375,
    0.625,
    0.875,
    60.0,
    90.0,
    100.0,
    180.0,
    300.0,
    360.0,
    1000.0,
    1024.0,
    3600.0,
    1000000.0,
    1048576.0,
}

FORBIDDEN_DOMAIN_NUMBERS = {
    5090,
    8233,
    12668,
    18426,
    54448,
    120307,
    339698,
    8111171,
    9643769,
    9953338,
    10195118,
    13067287,
    0.9858,
    0.8956,
    0.9699,
    0.9651,
    0.6622,
    0.6624,
    0.9367,
    0.8920,
    0.9882,
    812.75,
    501.83,
    213.7,
    106.47,
    228.9,
    760.96,
    1985,
    9486,
    11612,
    13148,
    13152,
    455937,
    32683955,
    120999311,
    0.6389,
    0.6429,
    0.8293,
    0.5122,
    0.7293,
    0.3233,
}

FORBIDDEN_DOMAIN_STRINGS = {
    "0.9858",
    "0.8956",
    "0.9699",
    "0.9651",
    "812.75",
    "501.83",
    "339,698",
    "8,111,171",
    "9,643,769",
    "13,067,287",
    "10,195,118",
    "213.70",
    "106.47",
    "11,612",
    "9,486",
    "1,985",
    "0.6389",
    "0.6429",
}


class TestBuildFiguresAndTables(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if SCRIPT.is_file():
            res = subprocess.run(
                [sys.executable, str(SCRIPT), "--input", "benchmarks/gcp-study-v106/master_gcp_metrics.json"],
                cwd=str(ROOT),
                capture_output=True,
                text=True,
            )
            if res.returncode != 0:
                raise RuntimeError(f"build_figures_and_tables.py failed:\nstderr={res.stderr}\nstdout={res.stdout}")

    def test_build_and_verify_all_figures_and_tables(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        res = subprocess.run(
            [sys.executable, str(SCRIPT), "--verify"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=f"stderr={res.stderr}\nstdout={res.stdout}")

        for i in range(1, 10):
            svg = FIG_DIR / f"fig{i}.svg"
            pdf = FIG_DIR / f"fig{i}.pdf"
            self.assertTrue(svg.is_file() and svg.stat().st_size > 2000, f"Missing/empty {svg}")
            self.assertTrue(pdf.is_file() and pdf.stat().st_size > 2000, f"Missing/empty {pdf}")

        for i in range(1, 10):
            md = TAB_DIR / f"table{i}.md"
            tex = TAB_DIR / f"table{i}.tex"
            self.assertTrue(md.is_file() and md.stat().st_size > 200, f"Missing/empty {md}")
            self.assertTrue(tex.is_file() and tex.stat().st_size > 200, f"Missing/empty {tex}")
            self.assertIsNone(re.search(r"\b(nan|NaN)\b", md.read_text(encoding="utf-8")))

        for i in range(1, 8):
            smd = TAB_DIR / f"supp_table{i}.md"
            stex = TAB_DIR / f"supp_table{i}.tex"
            self.assertTrue(smd.is_file() and smd.stat().st_size > 200, f"Missing/empty {smd}")
            self.assertTrue(stex.is_file() and stex.stat().st_size > 200, f"Missing/empty {stex}")
            self.assertIsNone(re.search(r"\b(nan|NaN)\b", smd.read_text(encoding="utf-8")))

        # Ensure script never reads root benchmarks/evidence/ (single-host GCP isolation)
        src = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("benchmarks/evidence", src)
        ast.parse(src)

    def test_all_9_figures_generated_and_valid_svg_pdf(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        for i in range(1, 10):
            svg = FIG_DIR / f"fig{i}.svg"
            pdf = FIG_DIR / f"fig{i}.pdf"
            self.assertTrue(svg.is_file(), f"Missing {svg}")
            self.assertTrue(pdf.is_file(), f"Missing {pdf}")
            self.assertGreater(svg.stat().st_size, 2000, f"SVG too small: {svg}")
            self.assertGreater(pdf.stat().st_size, 2000, f"PDF too small: {pdf}")

            svg_text = svg.read_text(encoding="utf-8")
            # Validate XML well-formedness and searchable <text> elements (svg.fonttype = 'none')
            root = ET.fromstring(svg_text)
            self.assertTrue(root.tag.endswith("svg"), f"Invalid root tag in {svg}: {root.tag}")
            self.assertIn("<text", svg_text, f"Expected searchable <text> elements in {svg}")

            # Validate PDF header and EOF
            pdf_bytes = pdf.read_bytes()
            self.assertTrue(pdf_bytes.startswith(b"%PDF-"), f"Invalid PDF header in {pdf}")
            self.assertIn(b"%%EOF", pdf_bytes[-1024:], f"Missing %%EOF in {pdf}")

    def test_all_16_tables_generated_with_zero_nan(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        table_stems = [f"table{i}" for i in range(1, 10)] + [f"supp_table{i}" for i in range(1, 8)]
        self.assertEqual(len(table_stems), 16)
        for stem in table_stems:
            md = TAB_DIR / f"{stem}.md"
            tex = TAB_DIR / f"{stem}.tex"
            self.assertTrue(md.is_file() and md.stat().st_size > 200, f"Missing/empty {md}")
            self.assertTrue(tex.is_file() and tex.stat().st_size > 200, f"Missing/empty {tex}")
            md_txt = md.read_text(encoding="utf-8")
            tex_txt = tex.read_text(encoding="utf-8")
            self.assertIsNone(
                re.search(r"\b(nan|NaN|None|TODO|TBD)\b", md_txt),
                f"Forbidden token in {md}",
            )
            self.assertIsNone(
                re.search(r"\b(nan|NaN|None|TODO|TBD)\b", tex_txt),
                f"Forbidden token in {tex}",
            )
            self.assertIn("|", md_txt, f"Expected Markdown table pipes in {md}")
            self.assertIn("\\begin{tabular}", tex_txt, f"Expected LaTeX tabular in {tex}")
            self.assertIn(
                f"\\label{{tab:{stem.replace('_', '-')}}}",
                tex_txt,
                f"Expected clean hyphenated LaTeX label in {tex}",
            )

    def test_all_60_cells_represented_in_tables(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        with INV_PATH.open(encoding="utf-8") as f:
            inv_rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(inv_rows), 60)

        supp1_md = (TAB_DIR / "supp_table1.md").read_text(encoding="utf-8")
        supp1_tex = (TAB_DIR / "supp_table1.tex").read_text(encoding="utf-8")
        supp_tables_md = "\n".join(
            (TAB_DIR / f"supp_table{i}.md").read_text(encoding="utf-8") for i in range(1, 7)
        )

        for row in inv_rows:
            cid = row["cell_id"]
            self.assertIn(cid, supp1_md, f"Cell {cid} missing from supp_table1.md")
            self.assertIn(cid.replace("_", "\\_"), supp1_tex, f"Cell {cid} missing from supp_table1.tex")
            self.assertIn(cid, supp_tables_md, f"Cell {cid} missing from supp_table1..6.md")

        clean_table_stems = [f"table{i}" for i in range(1, 10)] + [f"supp_table{i}" for i in range(2, 8)]
        for stem in clean_table_stems:
            t_md_path = TAB_DIR / f"{stem}.md"
            t_md = t_md_path.read_text(encoding="utf-8")
            self.assertIsNone(
                re.search(
                    r"\b(s[1-6]_[a-z0-9_]+|sweep_[a-z0-9_]+|arm[1-6]_[a-z0-9_]*|evidence_certified)\b",
                    t_md,
                ),
                f"Raw internal cell_id or arm slug found in {t_md_path}",
            )
            for bad_tok in (
                "Option C",
                "Tier P",
                "Tier-P",
                "H1a",
                "Gate 7",
                "Stage 1",
                "Stage 2",
                "Stage 3",
                "Stage 4",
                "Stage 5",
                "Stage 6",
            ):
                self.assertNotIn(
                    bad_tok,
                    t_md,
                    f"Internal harness jargon {bad_tok!r} found in {t_md_path}",
                )

    def test_all_tables_have_citations_and_superscript_notes(self) -> None:
        table_stems = [f"table{i}" for i in range(1, 10)] + [f"supp_table{i}" for i in range(1, 8)]
        for stem in table_stems:
            md_path = TAB_DIR / f"{stem}.md"
            tex_path = TAB_DIR / f"{stem}.tex"
            md_txt = md_path.read_text(encoding="utf-8")
            tex_txt = tex_path.read_text(encoding="utf-8")

            self.assertRegex(md_txt, r"\[\d+(?:[–,-]\s*\d+)*\]", f"Missing bracketed citation [N] in {md_path}")
            self.assertIn("<sup>a</sup>", md_txt, f"Missing superscript note <sup>a</sup> in {md_path}")
            self.assertNotIn("`<sup>", md_txt, f"Backtick-escaped <sup> found in {md_path}")

            self.assertIn("$^{\\text{a}}$", tex_txt, f"Missing LaTeX superscript $^{{\\text{{a}}}}$ in {tex_path}")
            self.assertNotIn("<sup>", tex_txt, f"Unconverted raw HTML <sup> found in {tex_path}")
            self.assertIn("\\resizebox{\\textwidth}{!}{%", tex_txt, f"Missing \\resizebox in {tex_path}")
            self.assertIsNone(
                re.search(r"\\caption\{(Table|Supplementary Table)\s+", tex_txt),
                f"Duplicate leading Table prefix in LaTeX \\caption in {tex_path}",
            )

            for stale_phrase in (
                "Wall-Clock Timeout (>10 min)",
                "1 Threads",
                "1 sub-arms",
                "1 Sub-Arms",
                "17-Run",
                "25 evaluations",
            ):
                self.assertNotIn(stale_phrase, md_txt, f"Stale/inaccurate phrase {stale_phrase!r} in {md_path}")
                self.assertNotIn(stale_phrase, tex_txt, f"Stale/inaccurate phrase {stale_phrase!r} in {tex_path}")

    def test_no_hardcoded_numbers_in_builder_script(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        src = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("benchmarks/evidence", src)
        tree = ast.parse(src, filename=str(SCRIPT))

        violations: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant):
                val = node.value
                lineno = getattr(node, "lineno", "?")
                if isinstance(val, bool) or val is None:
                    continue
                if isinstance(val, int):
                    if val in FORBIDDEN_DOMAIN_NUMBERS or val not in ALLOWED_INT_CONSTANTS:
                        violations.append(f"line {lineno}: disallowed int literal {val}")
                elif isinstance(val, float):
                    if val in FORBIDDEN_DOMAIN_NUMBERS:
                        violations.append(f"line {lineno}: forbidden domain float literal {val}")
                    elif val not in ALLOWED_FLOAT_MULTIPLIERS:
                        # Allow styling/layout floats in [-25.0, 25.0] with at most 2 decimal places
                        if not (-25.0 <= val <= 25.0 and round(val, 2) == val):
                            violations.append(f"line {lineno}: disallowed float literal {val}")
                elif isinstance(val, str):
                    for bad_s in FORBIDDEN_DOMAIN_STRINGS:
                        if bad_s in val:
                            violations.append(f"line {lineno}: hardcoded metric string {bad_s!r} in {val!r}")
                    m_pat = re.search(r"\d+\.\d+x|\d+\.\d+%|\d{1,3},\d{3}", val)
                    if m_pat:
                        violations.append(
                            f"line {lineno}: hardcoded empirical ratio/percentage/count {m_pat.group(0)!r} in {val!r}"
                        )

        self.assertEqual(
            violations,
            [],
            "Hardcoded domain numbers found in paper/build_figures_and_tables.py:\n"
            + "\n".join(f"  - {v}" for v in violations),
        )

    def test_tables_free_of_false_cpp_mismatches_and_stale_sweep_diffs(self) -> None:
        supp1_md = (TAB_DIR / "supp_table1.md").read_text(encoding="utf-8")
        for line in supp1_md.splitlines():
            if "s5_wgs_hg002_full_giab_3way_th32" in line or "s6_low_vaf_titration_truth_curve" in line:
                self.assertNotIn("MISMATCH", line, f"False vardictcpp MISMATCH leaked into Supp Table S1: {line}")

        table1_md = (TAB_DIR / "table1.md").read_text(encoding="utf-8")
        self.assertIn("87/87 Sweep Sub-Arms Exact Match", table1_md)
        self.assertIn("24/24 Sub-Arms Exact Match", table1_md)
        self.assertIn("30/30 Sub-Arms Exact Match", table1_md)
        self.assertIn("1/1 Sub-Arm Exact Match", table1_md)
        self.assertIn("2/2 Sub-Arms Exact Match", table1_md)

    def test_v110_fig9_table9_and_supp_table7_hygiene_and_metrics(self) -> None:
        self.assertTrue(V110_PATH.is_file(), f"Missing {V110_PATH}")
        v110_doc = json.loads(V110_PATH.read_text(encoding="utf-8"))

        res = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--verify",
                "--v110-input",
                str(V110_PATH),
            ],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=f"stderr={res.stderr}\nstdout={res.stdout}")

        fig9_svg = (FIG_DIR / "fig9.svg").read_text(encoding="utf-8")
        ET.fromstring(fig9_svg)
        self.assertIn("SEQC2", fig9_svg)
        self.assertIn(">0<", re.sub(r"\s+", "", fig9_svg))

        t9_md = (TAB_DIR / "table9.md").read_text(encoding="utf-8")
        s7_md = (TAB_DIR / "supp_table7.md").read_text(encoding="utf-8")

        t9_data_rows = [
            ln for ln in t9_md.splitlines()
            if ln.strip().startswith("|") and not re.match(r"^\|\s*[:\-| ]+\|\s*$", ln.strip())
        ][1:]
        self.assertEqual(len(t9_data_rows), 6, f"Expected 6 empirical domain rows in table9.md, got {len(t9_data_rows)}")

        for sup_letter in ("a", "b", "c", "d", "e", "f"):
            self.assertIn(f"<sup>{sup_letter}</sup>", t9_md, f"Missing <sup>{sup_letter}</sup> in table9.md")
        self.assertIn(v110_doc["meta"]["study_fingerprint"], t9_md)
        for ref_num in ("1", "10", "12", "13", "19", "26", "27"):
            self.assertRegex(t9_md, rf"\[[^\]]*\b{ref_num}\b[^\]]*\]", f"Missing citation [{ref_num}] in table9.md")
        self.assertNotIn("43/91", t9_md, "Mismatched QIAseq -sp numerator/denominator (43/91) in table9.md")

        s7_data_rows = [
            ln for ln in s7_md.splitlines()
            if ln.strip().startswith("|") and not re.match(r"^\|\s*[:\-| ]+\|\s*$", ln.strip())
        ][1:]
        self.assertEqual(len(s7_data_rows), 35, f"Expected 35 switch_catalog rows in supp_table7.md, got {len(s7_data_rows)}")
        for entry in v110_doc["switch_catalog"]:
            self.assertIn(entry["name"], s7_md, f"Missing switch_catalog entry {entry['name']} in supp_table7.md")
        for sup_letter in ("a", "b", "c"):
            self.assertIn(f"<sup>{sup_letter}</sup>", s7_md, f"Missing <sup>{sup_letter}</sup> in supp_table7.md")
        for ref_num in ("1", "10", "12", "26", "27"):
            self.assertRegex(s7_md, rf"\[[^\]]*\b{ref_num}\b[^\]]*\]", f"Missing citation [{ref_num}] in supp_table7.md")

        bare_issue_re = re.compile(r"(?<![A-Za-z0-9_/-])#\d+\b")
        for label, md_txt in (("table9.md", t9_md), ("supp_table7.md", s7_md)):
            body_no_heading = "\n".join(
                ln for ln in md_txt.splitlines() if not ln.startswith("### ")
            )
            self.assertIsNone(
                bare_issue_re.search(body_no_heading),
                f"Bare #N issue reference found in {label}",
            )
            for code_span in re.findall(r"`([^`\n]+)`", md_txt):
                self.assertNotIn("|", code_span, f"Pipe character inside inline code span in {label}: {code_span!r}")
            stripped_code = re.sub(r"`[^`\n]+`", "", body_no_heading)
            stripped_sup = re.sub(r"</?(?:sup|sub)>", "", stripped_code)
            self.assertIsNone(
                re.search(r"<[A-Za-z0-9_\-±.+]+>", stripped_sup),
                f"Unbackticked <...> angle-bracket token in {label}",
            )



if __name__ == "__main__":
    unittest.main()
