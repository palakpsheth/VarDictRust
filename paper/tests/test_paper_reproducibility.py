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
MANUSCRIPT_PDF = ROOT / "paper" / "manuscript.pdf"
SUPPLEMENT_PDF = ROOT / "paper" / "supplement.pdf"
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
    "218.83",
    "219.25",
    "5.05",
    "12.43",
    "106.47",
    "228.90",
    "760.96",
    "12,668",
    "54,448",
    "11,612",
    "13,148",
    "13,152",
    "455,937",
    "0.6389",
    "0.6429",
    "0.5610",
    "0.8780",
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
    "fig9_v110_compat_and_presets",
]

EXPECTED_TABLE_IDS = [
    "table1_parity_matrix",
    "table2_tierp_and_rchr1",
    "table3_truth_and_field_callers",
    "table4_memory_scaling_and_packing",
    "table5_full_wgs_arm_a_vs_b",
    "table6_cross_platform_chr20",
    "table7_v50q_cmrg_mosaic_and_hg008_v04",
    "table8_supply_chain_and_divergences",
    "table9_v110_compat_and_presets",
]

EXPECTED_SUPP_TABLE_IDS = [
    "supp_table1_60cell_inventory",
    "supp_table2_option_c_87arms",
    "supp_table3_n3_technical_replicates",
    "supp_table4_qual_roc_and_noise_floor",
    "supp_table5_ga4gh_70cell_matrix",
    "supp_table6_sv_and_amplicon_splice",
    "supp_table7_compat_32switches_catalog",
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

        # 3. All 9 Figures, 9 Main Tables, and 7 Supplementary Tables referenced
        for i in range(1, 10):
            self.assertIn(f"Figure {i}", ms_text)
            self.assertIn(f"Table {i}", ms_text)
        for i in range(1, 8):
            self.assertIn(f"Table S{i}", combined)

        # 4. Required Genome Biology / Bioinformatics sections & disclosures present
        for req in (
            "### Background",
            "### Results",
            "### Conclusions",
            "### Availability of data and materials",
            "### AI-assisted development and transparency disclosure",
            "Joachim Wolff",
            "10.1038/s41597-022-01359-6",
            "--nosv",
        ):
            self.assertIn(req, ms_text, f"Missing required token/section in manuscript.md: {req}")
        for sup_req in (
            "5865a2a0b24a8c17",
            "c490c36317b8f25b",
            "b18477ae0a7b07c1",
            "Liss_lab",
            "10.1038/s41597-022-01359-6",
        ):
            self.assertIn(sup_req, sup_text, f"Missing required provenance token in supplement.md: {sup_req}")

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

    def test_all_9_figures_and_16_tables_referenced(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")
        combined = ms_text + "\n" + sup_text

        for i in range(1, 10):
            self.assertIn(f"Figure {i}", ms_text)
            self.assertIn(f"figures/fig{i}.svg", ms_text)
            self.assertIn(f"Table {i}", ms_text)
        for i in range(1, 8):
            self.assertIn(f"Table S{i}", combined)

        for fid in EXPECTED_FIGURE_IDS:
            self.assertIn(fid, combined, f"Missing figure identifier {fid}")
        for tid in EXPECTED_TABLE_IDS:
            self.assertIn(tid, combined, f"Missing main table identifier {tid}")
        for stid in EXPECTED_SUPP_TABLE_IDS:
            self.assertIn(stid, combined, f"Missing supplementary table identifier {stid}")

    def test_structured_abstract_word_count_and_key_points(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        m_abs = re.search(
            r"^## Abstract\s*\n([\s\S]+?)\n\*\*Keywords:\*\*",
            ms_text,
            re.MULTILINE,
        )
        self.assertIsNotNone(m_abs, "Could not locate ## Abstract ... **Keywords:** block in manuscript.md")
        assert m_abs is not None
        abs_lines = [
            ln.strip()
            for ln in m_abs.group(1).splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")
        ]
        abs_words = " ".join(abs_lines).split()
        self.assertLessEqual(
            len(abs_words),
            350,
            f"Structured Abstract exceeds 350-word Genome Biology limit: {len(abs_words)} words",
        )
        self.assertGreaterEqual(
            len(abs_words),
            200,
            f"Structured Abstract unexpectedly short: {len(abs_words)} words",
        )

        m_kp = re.search(
            r"^>\s*\*\*Key Points\*\*\s*\n((?:>.*\n)+)",
            ms_text,
            re.MULTILINE,
        )
        self.assertIsNotNone(m_kp, "Missing > **Key Points** callout box in manuscript.md")
        assert m_kp is not None
        kp_bullets = [
            ln for ln in m_kp.group(1).splitlines() if re.match(r"^>\s*-\s+", ln)
        ]
        self.assertGreaterEqual(len(kp_bullets), 3, "Key Points box must contain >= 3 bullets")
        self.assertLessEqual(len(kp_bullets), 5, f"Key Points box exceeds 5 bullets: {len(kp_bullets)}")

    def test_all_60_cells_referenced_across_manuscript_and_supplement(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")

        with INV_PATH.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 60)
        for r in rows:
            cid = r["cell_id"]
            self.assertIn(cid, sup_text, f"Cell {cid} missing from supplement.md")
            self.assertNotIn(cid, ms_text, f"Internal cell_id {cid} leaked into main manuscript.md")

    def test_manuscript_free_of_internal_plan_jargon(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        forbidden_jargon = [
            "Option C",
            "Tier-P",
            "Tier P",
            "H1a",
            "Gate 7",
            "evidence_certified",
            "c3d-standard-32",
            "228.90",
            "Stage 1",
            "Stage 2",
            "Stage 3",
            "Stage 4",
            "Stage 5",
            "Stage 6",
            "s1_",
            "s2_",
            "s3_",
            "s4_",
            "s5_",
            "s6_",
            "arm1_",
            "arm2_",
            "arm3_",
            "arm4_",
            "arm5_",
            "arm6_",
        ]
        for term in forbidden_jargon:
            self.assertNotIn(
                term,
                ms_text,
                f"Internal plan jargon {term!r} found in main manuscript.md",
            )

    def test_prose_free_of_ai_cliches_and_throat_clearing(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")
        combined = ms_text + "\n" + sup_text
        banned_ai_phrases = [
            "central challenge",
            "cornerstone",
            "recurring tension",
            "underappreciated consequence",
            "merit explicit consideration",
            "architectural divide",
            "meticulously",
            "ballooned",
            "bloated",
            "delve",
            "tapestry",
            "testament",
            "paradigm shift",
            "underscores",
            "To ensure complete machine-traceable auditability",
            "claims of Wolff et al. [11], we replicated the Wolff et al. [11]",
            "palak@rewrites.bio",
            "github.com/rewrites-bio",
            "rewrites.bio) Initiative",
            "rewrites.bio) initiative",
            "Lukas Wolff",
            "10.1038/s41597-022-01359-2",
            "Landscape",
            "Clinical AML AML",
            "Scratch OCI",
            "matching `--preset wes`'s",
            "matching --preset wes's",
            "We thank",
            "Table 9,&",
        ]
        for phrase in banned_ai_phrases:
            self.assertNotIn(
                phrase,
                combined,
                f"Banned AI cliché or throat-clearing phrase {phrase!r} found in manuscript/supplement",
            )
        for req_author in (
            "https://www.linkedin.com/in/palak-p-sheth/",
            "palakpsheth@gmail.com",
            "https://github.com/palakpsheth/VarDictRust",
        ):
            self.assertIn(req_author, ms_text, f"Missing author/repo metadata {req_author!r} in manuscript.md")
        self.assertIn("https://www.linkedin.com/in/palak-p-sheth/", sup_text)
        self.assertIn("palakpsheth@gmail.com", sup_text)

    def test_referee_and_biostatistics_compliance(self) -> None:
        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")

        ms_required_phrases = [
            "Clopper-Pearson",
            "technical replicates",
            "Mann-Whitney",
            "--compat=off",
            "--nosv",
            "in silico",
            "```bash",
        ]
        for phrase in ms_required_phrases:
            self.assertIn(
                phrase,
                ms_text,
                f"Missing required scientific/biostatistics disclosure in manuscript.md: {phrase!r}",
            )

        sup_required_phrases = [
            "s2_fig1_fig3_field_callers",
            "rc=124",
            "rc=3",
            "Liss_lab",
            "c3d-standard-32",
            "2a8e9a3",
            "```bash",
        ]
        for phrase in sup_required_phrases:
            self.assertIn(
                phrase,
                sup_text,
                f"Missing required provenance disclosure in supplement.md: {phrase!r}",
            )

        bare_issue_re = re.compile(r"(?<![A-Za-z0-9_/-])#\d+\b")
        for tmpl_path in (MANUSCRIPT_TMPL, SUPPLEMENT_TMPL):
            raw_txt = tmpl_path.read_text(encoding="utf-8")
            no_fenced = re.sub(r"```[\s\S]*?```", "", raw_txt)
            body_lines = [ln for ln in no_fenced.splitlines() if not ln.lstrip().startswith("#")]
            m_bare = bare_issue_re.search("\n".join(body_lines))
            self.assertIsNone(
                m_bare,
                f"Bare #N issue reference {m_bare.group(0) if m_bare else ''!r} found in {tmpl_path.name}; repo-qualify as owner/repo#N",
            )

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

    def test_manuscript_and_supplement_pdfs_exist_and_valid(self) -> None:
        for pdf_path in (MANUSCRIPT_PDF, SUPPLEMENT_PDF):
            self.assertTrue(pdf_path.is_file(), f"Missing {pdf_path}")
            self.assertGreater(pdf_path.stat().st_size, 50_000, f"PDF unexpectedly small: {pdf_path}")
            raw = pdf_path.read_bytes()
            self.assertTrue(raw.startswith(b"%PDF-"), f"Invalid PDF header in {pdf_path}")
            self.assertIn(b"%%EOF", raw[-1024:], f"Missing %%EOF trailer in {pdf_path}")

    def test_paper_readme_cover_letter_matches_certified_macros(self) -> None:
        import json as _json

        sys.path.insert(0, str(ROOT / "paper"))
        import build_manuscript as bm  # type: ignore[import-not-found]

        gcp_doc = _json.loads(
            (ROOT / "benchmarks" / "gcp-study-v106" / "master_gcp_metrics.json").read_text(encoding="utf-8")
        )
        v110_doc = _json.loads(
            (ROOT / "benchmarks" / "v110-compat-study" / "v110_compat_metrics.json").read_text(encoding="utf-8")
        )
        macros = bm.build_macro_context(gcp_doc, v110_doc)
        readme_txt = (ROOT / "paper" / "README.md").read_text(encoding="utf-8")

        for macro_key in (
            "s5_h4_conda_pkgs",
            "s5_h4_musl_mb",
            "s2_rchr1_speedup_vs_java",
            "s2_rchr1_rss_ratio_vs_java",
            "s6_v50q_snv_r_gain_pct",
            "s6_v50q_indel_r_gain_pct",
            "s6_cmrg_snv_r_gain_pct",
            "s6_cmrg_indel_r_gain_pct",
            "v110_schema_defects_vj",
            "v110_arm2_f1_gain",
        ):
            self.assertIn(macro_key, macros, f"Missing macro {macro_key} in build_macro_context()")
            self.assertIn(
                macros[macro_key],
                readme_txt,
                f"paper/README.md missing certified macro {macro_key}={macros[macro_key]!r}",
            )

        self.assertIn("--reaggregate-gcp", readme_txt)
        self.assertIn("+31.7 pp", readme_txt)
        for stale_str in ("Landscape", "ecosystem", "121` Conda", "3.46 MB", "38.2×", "110.2×"):
            self.assertNotIn(stale_str, readme_txt, f"Stale/banned string {stale_str!r} in paper/README.md")

    def test_stage5_rss_sync_aggregate_verify_fail_closed_and_pdf_zwsp(self) -> None:
        import tempfile as _tempfile

        ms_text = MANUSCRIPT.read_text(encoding="utf-8")
        sup_text = SUPPLEMENT.read_text(encoding="utf-8")
        ms_plain = ms_text.replace("`", "")
        sup_plain = sup_text.replace("`", "")

        for expected_mem in ("5.05 GB", "12.43 GB", "218.83 GB", "32.68 GB", "94.06 GB"):
            self.assertIn(expected_mem, ms_plain, f"Expected Stage-5 RSS {expected_mem!r} in manuscript.md")
            self.assertIn(expected_mem, sup_plain, f"Expected Stage-5 RSS {expected_mem!r} in supplement.md")
        for stale_mem in ("4.93 GB", "12.14 GB", "213.7 GB", "214.1 GB"):
            self.assertNotIn(stale_mem, ms_plain, f"Stale GiB-labeled-as-GB {stale_mem!r} in manuscript.md")
            self.assertNotIn(stale_mem, sup_plain, f"Stale GiB-labeled-as-GB {stale_mem!r} in supplement.md")

        self.assertIn("29.4×", ms_plain)
        self.assertNotIn("29.3×", ms_plain)
        self.assertNotIn("29.3x", ms_plain)
        self.assertIn("0.00049", ms_plain)
        self.assertIn("0.00049", sup_plain)

        # aggregate_gcp_study.py --verify must fail closed when master_gcp_metrics.json is missing
        with _tempfile.TemporaryDirectory() as tmpdir:
            tmp_out = Path(tmpdir) / "out"
            tmp_raw = Path(tmpdir) / "raw"
            tmp_out.mkdir()
            tmp_raw.mkdir()
            res = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "paper" / "aggregate_gcp_study.py"),
                    "--verify",
                    "--out-dir",
                    str(tmp_out),
                    "--raw-dir",
                    str(tmp_raw),
                ],
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(res.returncode, 0)
            self.assertFalse((tmp_out / "master_gcp_metrics.json").exists())

        # render_pdf.py must not split 0/1 or 34/41 with &#8203; and must protect short code hyphens
        sys.path.insert(0, str(ROOT / "paper"))
        import render_pdf as rp  # type: ignore[import-not-found]

        rendered = rp._postprocess_body_html(
            "<td>0/1 and 1/0 and N/N and 34/41 and path_one/path_two</td><p><code>--compat=off</code></p>"
        )
        self.assertIn("0/1", rendered)
        self.assertIn("34/41", rendered)
        self.assertIn("path_&#8203;one/&#8203;path_&#8203;two", rendered)
        self.assertIn("&#8209;&#8209;compat=off", rendered)


if __name__ == "__main__":
    unittest.main()
