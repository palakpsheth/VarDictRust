#!/usr/bin/env python3
"""Unit tests for paper/build_figures_and_tables.py (offline reproducibility & zero hardcoded metrics)."""

from __future__ import annotations

import ast
import csv
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

        for i in range(1, 9):
            svg = FIG_DIR / f"fig{i}.svg"
            pdf = FIG_DIR / f"fig{i}.pdf"
            self.assertTrue(svg.is_file() and svg.stat().st_size > 2000, f"Missing/empty {svg}")
            self.assertTrue(pdf.is_file() and pdf.stat().st_size > 2000, f"Missing/empty {pdf}")

        for i in range(1, 9):
            md = TAB_DIR / f"table{i}.md"
            tex = TAB_DIR / f"table{i}.tex"
            self.assertTrue(md.is_file() and md.stat().st_size > 200, f"Missing/empty {md}")
            self.assertTrue(tex.is_file() and tex.stat().st_size > 200, f"Missing/empty {tex}")
            self.assertIsNone(re.search(r"\b(nan|NaN)\b", md.read_text(encoding="utf-8")))

        for i in range(1, 7):
            smd = TAB_DIR / f"supp_table{i}.md"
            stex = TAB_DIR / f"supp_table{i}.tex"
            self.assertTrue(smd.is_file() and smd.stat().st_size > 200, f"Missing/empty {smd}")
            self.assertTrue(stex.is_file() and stex.stat().st_size > 200, f"Missing/empty {stex}")
            self.assertIsNone(re.search(r"\b(nan|NaN)\b", smd.read_text(encoding="utf-8")))

        # Ensure script never reads root benchmarks/evidence/ (single-host GCP isolation)
        src = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("benchmarks/evidence", src)
        ast.parse(src)

    def test_all_8_figures_generated_and_valid_svg_pdf(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        for i in range(1, 9):
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

    def test_all_14_tables_generated_with_zero_nan(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        table_stems = [f"table{i}" for i in range(1, 9)] + [f"supp_table{i}" for i in range(1, 7)]
        self.assertEqual(len(table_stems), 14)
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

    def test_all_60_cells_represented_in_tables(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        with INV_PATH.open(encoding="utf-8") as f:
            inv_rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(inv_rows), 60)

        supp1_md = (TAB_DIR / "supp_table1.md").read_text(encoding="utf-8")
        supp1_tex = (TAB_DIR / "supp_table1.tex").read_text(encoding="utf-8")
        all_tables_md = "\n".join(
            p.read_text(encoding="utf-8") for p in sorted(TAB_DIR.glob("*.md"))
        )

        for row in inv_rows:
            cid = row["cell_id"]
            self.assertIn(cid, supp1_md, f"Cell {cid} missing from supp_table1.md")
            self.assertIn(cid.replace("_", "\\_"), supp1_tex, f"Cell {cid} missing from supp_table1.tex")
            self.assertIn(cid, all_tables_md, f"Cell {cid} missing from generated markdown tables")

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

        self.assertEqual(
            violations,
            [],
            "Hardcoded domain numbers found in paper/build_figures_and_tables.py:\n"
            + "\n".join(f"  - {v}" for v in violations),
        )


if __name__ == "__main__":
    unittest.main()
