#!/usr/bin/env python3
"""Unit tests for paper/collect_v110_compat_metrics.py and benchmarks/v110-compat-study/."""

from __future__ import annotations

import ast
import csv
import importlib.util
import json
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import tomllib
import unittest
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "paper" / "collect_v110_compat_metrics.py"
TEST_FILE = Path(__file__).resolve()
OUT_DIR = ROOT / "benchmarks" / "v110-compat-study"
CONFIG_RS = ROOT / "crates" / "vardictrust-core" / "src" / "config.rs"
ROADMAP_MD = ROOT / "ROADMAP.md"
UPSTREAM_TOML = ROOT / "UPSTREAM.toml"
ECOSYSTEM_SURVEY_MD = ROOT / "docs" / "ecosystem-survey.md"


def _load_collector_module():
    spec = importlib.util.spec_from_file_location("collect_v110_compat_metrics", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestCollectV110CompatMetrics(unittest.TestCase):
    def test_verify_mode_succeeds_offline_without_tools_or_datasets(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing collector script: {SCRIPT}")
        self.assertTrue(OUT_DIR.is_dir(), f"Missing study bundle directory: {OUT_DIR}")
        self.assertEqual(
            list(OUT_DIR.rglob("*.md")),
            [],
            "benchmarks/v110-compat-study/ must contain zero .md files (Zensical not_in_nav rule)",
        )

        # Run --verify in a stripped environment with no bioinformatics tools or data paths
        stripped_env = {
            "PATH": "/usr/bin:/bin",
            "HOME": "/nonexistent_home",
            "G4_TOOLS_BIN": "/nonexistent_tools",
            "CONDA_ROOT": "/nonexistent_conda",
            "DATA_ROOT": "/nonexistent_data",
            "LC_ALL": "C.UTF-8",
        }
        res = subprocess.run(
            [sys.executable, str(SCRIPT), "--verify", "--out-dir", str(OUT_DIR)],
            cwd=str(ROOT),
            env=stripped_env,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(
            res.returncode,
            0,
            msg=f"Offline --verify failed:\nstderr={res.stderr}\nstdout={res.stdout}",
        )

    def test_tampered_bundle_files_rejected_by_verify(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing collector script: {SCRIPT}")
        self.assertTrue((OUT_DIR / "v110_compat_metrics.json").is_file())
        self.assertTrue((OUT_DIR / "v110_compat_arms.tsv").is_file())
        self.assertTrue((OUT_DIR / "SHA256SUMS.tsv").is_file())

        for target_name in ("v110_compat_metrics.json", "v110_compat_arms.tsv"):
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_bundle = Path(tmpdir) / "bundle"
                shutil.copytree(OUT_DIR, tmp_bundle)
                victim = tmp_bundle / target_name
                raw = bytearray(victim.read_bytes())
                raw[0] ^= 0x01
                victim.write_bytes(bytes(raw))

                res = subprocess.run(
                    [sys.executable, str(SCRIPT), "--verify", "--out-dir", str(tmp_bundle)],
                    cwd=str(ROOT),
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertNotEqual(
                    res.returncode,
                    0,
                    f"Expected --verify to reject 1-byte tamper in {target_name}",
                )

    def test_zero_path_leaks_and_sha256sums_format(self) -> None:
        metrics_path = OUT_DIR / "v110_compat_metrics.json"
        arms_path = OUT_DIR / "v110_compat_arms.tsv"
        sha_path = OUT_DIR / "SHA256SUMS.tsv"
        for p in (metrics_path, arms_path, sha_path):
            self.assertTrue(p.is_file(), f"Missing {p}")

        for p in (metrics_path, arms_path, sha_path):
            txt = p.read_text(encoding="utf-8")
            self.assertNotIn("/home/", txt, f"Local /home/ path leaked in {p.name}")
            self.assertNotIn(str(ROOT), txt, f"Workspace ROOT path leaked in {p.name}")
            self.assertNotIn("/tmp/", txt, f"Temporary /tmp/ path leaked in {p.name}")
            self.assertNotIn("vdr_v110_study_", txt, f"Tempdir prefix leaked in {p.name}")
            self.assertNotIn("miniconda3", txt, f"Conda path leaked in {p.name}")
            self.assertNotIn("miniforge3", txt, f"Conda path leaked in {p.name}")

        with sha_path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            self.assertEqual(reader.fieldnames, ["sha256", "path", "size_bytes"])
            rows = list(reader)

        self.assertEqual(
            [r["path"] for r in rows],
            ["v110_compat_arms.tsv", "v110_compat_metrics.json"],
        )
        mod = _load_collector_module()
        mod.verify_bundle(OUT_DIR)
        for r in rows:
            target = OUT_DIR / r["path"]
            self.assertEqual(mod.sha256_file(target), r["sha256"])
            self.assertEqual(str(target.stat().st_size), r["size_bytes"])

    def test_metrics_json_schema_and_six_arms_invariants(self) -> None:
        mod = _load_collector_module()
        mod.verify_bundle(OUT_DIR)
        metrics_path = OUT_DIR / "v110_compat_metrics.json"
        arms_path = OUT_DIR / "v110_compat_arms.tsv"
        doc = json.loads(metrics_path.read_text(encoding="utf-8"))
        a1 = doc["arms"]["arm1_vcf_schema_and_bcftools_norm"]
        for k in mod.DEFECT_KEYS:
            self.assertGreater(a1["vardictjava_defects"][k], 0, f"Expected VJ defect > 0 for {k}")
            self.assertEqual(a1["compat_off_defects"][k], 0, f"Expected compat_off defect == 0 for {k}")

        a6 = doc["arms"]["arm6_csi_spanning_sv_and_pileup_tiling"]
        self.assertEqual(a6["pileup_full_span_cases_total"], 2)
        self.assertEqual(a6["pileup_untiled_eq_tiled_cases"], a6["pileup_full_span_cases_total"])
        self.assertEqual(a6["pileup_covered_bases_retained"], a6["pileup_covered_bases_total"])

        # TSV rows
        with arms_path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            self.assertEqual(
                reader.fieldnames,
                [
                    "arm_id",
                    "arm_title",
                    "cohort_scope",
                    "provenance_type",
                    "vardictjava_result",
                    "compat_off_result",
                    "verdict",
                    "sanitized_cmd",
                ],
            )
            tsv_rows = list(reader)
        self.assertEqual(len(tsv_rows), 6)
        for row in tsv_rows:
            self.assertTrue(
                row["verdict"].startswith("PASS"),
                f"Expected PASS prefix in verdict, got {row['verdict']!r}",
            )
            tokens = shlex.split(row["sanitized_cmd"], posix=True)
            for idx_t, tok in enumerate(tokens):
                if idx_t > 0 and tokens[idx_t - 1] == "-b":
                    continue
                self.assertNotIn(
                    "|",
                    tok,
                    f"Unquoted shell pipe '|' in {row['arm_id']} sanitized_cmd: {row['sanitized_cmd']!r}",
                )
        self.assertIn('-b "$TUMOR_BAM|$NORMAL_BAM"', tsv_rows[0]["sanitized_cmd"])
        self.assertIn('-b "$TUMOR_BAM|$NORMAL_BAM"', tsv_rows[2]["sanitized_cmd"])
        self.assertNotIn("Strict VCFv4.2/v4.3 Schema", tsv_rows[0]["arm_title"])
        self.assertNotIn("strict schema-validator", doc["meta"]["scope_note"])

        # Arm 1 provenance clarity: 17 wet-lab BAM slices vs 1 synthetic SV BAM + IUPAC fixture
        self.assertIn("17", tsv_rows[0]["provenance_type"])
        self.assertIn("vartovcf_sv_contig1", tsv_rows[0]["provenance_type"])
        self.assertIn("IUPAC", tsv_rows[0]["provenance_type"])

        # Arm 2 full unthresholded + optimal-QUAL metrics and percentage-point (+31.7 pp) verdict
        for col_name in ("vardictjava_result", "compat_off_result"):
            self.assertIn("P=", tsv_rows[1][col_name])
            self.assertIn("R=", tsv_rows[1][col_name])
            self.assertIn("QUAL", tsv_rows[1][col_name])
        self.assertIn("+31.7 pp", tsv_rows[1]["verdict"])
        self.assertNotIn("+31.7%", tsv_rows[1]["verdict"])

        # Arm 5 full-panel /192 denominator alongside >=50x /91 denominator
        self.assertIn("24/192", tsv_rows[4]["vardictjava_result"])
        self.assertIn("14/91", tsv_rows[4]["vardictjava_result"])
        self.assertIn("43/192", tsv_rows[4]["compat_off_result"])
        self.assertIn("24/91", tsv_rows[4]["compat_off_result"])

    def test_explicit_minus_z_on_all_zero_based_bed_cli_calls_and_sanitize_order(self) -> None:
        src = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("import bisect", src, "Unused import bisect must be removed")
        tree = ast.parse(src, filename=str(SCRIPT))

        bed_col_lists_found = 0
        for node in ast.walk(tree):
            if isinstance(node, (ast.List, ast.Tuple)):
                str_elts = [
                    elt.value if isinstance(elt, ast.Constant) and isinstance(elt.value, str) else None
                    for elt in node.elts
                ]
                for i in range(len(str_elts) - 5):
                    if str_elts[i : i + 6] == ["-c", "1", "-S", "2", "-E", "3"]:
                        bed_col_lists_found += 1
                        prev_tok = str_elts[i - 1] if i > 0 else None
                        self.assertEqual(
                            prev_tok,
                            "-z",
                            f"Line {getattr(node, 'lineno', '?')}: explicit 0-based BED '-c 1 -S 2 -E 3' must have '-z' immediately before '-c'",
                        )
        self.assertGreaterEqual(bed_col_lists_found, 8)

        # Verify sanitize_text substitutes longer nested paths before shorter parent prefixes
        spec = importlib.util.spec_from_file_location("collect_v110_compat_metrics", SCRIPT)
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        # Worktree layout: data_root (/repo) is a parent prefix of root_dir (/repo/.worktrees/wt)
        sanitized = mod.sanitize_text(
            "cmd /repo/.worktrees/wt/sub/file.vcf /repo/.worktrees/wt/tmp_run/out.vcf",
            tmpdir="/repo/.worktrees/wt/tmp_run",
            data_root=Path("/repo"),
            root_dir=Path("/repo/.worktrees/wt"),
        )
        self.assertEqual(sanitized, "cmd $REPO_ROOT/sub/file.vcf $WORK_TMP/out.vcf")

    def test_wolff_2026_gigabyte_doi_fixed_in_upstream_bugs_md(self) -> None:
        for doc_path in (
            ROOT / "docs" / "upstream-bugs.md",
            ECOSYSTEM_SURVEY_MD,
            ROADMAP_MD,
            UPSTREAM_TOML,
        ):
            txt = doc_path.read_text(encoding="utf-8")
            self.assertIn("https://doi.org/10.46471/gigabyte.194", txt, f"Missing verified Wolff DOI in {doc_path}")
            self.assertNotIn("10.1093/gigabyte/giaf156", txt, f"Stale 404 Wolff DOI left in {doc_path}")

    def test_switch_catalog_matches_config_rs_in_order_and_hygiene(self) -> None:
        mod = _load_collector_module()
        expected_switches = mod.parse_compat_profile_switches(CONFIG_RS, ROADMAP_MD, require_exact_count=False)
        self.assertGreaterEqual(len(expected_switches), 32)

        metrics_path = OUT_DIR / "v110_compat_metrics.json"
        self.assertTrue(metrics_path.is_file(), f"Missing {metrics_path}")
        doc = json.loads(metrics_path.read_text(encoding="utf-8"))
        catalog = doc["switch_catalog"]
        self.assertEqual(len(catalog), 35)

        got_switches = [entry["name"] for entry in catalog[:32]]
        self.assertEqual(got_switches, expected_switches[:32])

        # Verify v1.1.1+ future growth (33rd switch in config.rs + ROADMAP.md) passes offline verify
        with tempfile.TemporaryDirectory() as tmpdir:
            tp = Path(tmpdir)
            cfg_copy = tp / "config.rs"
            rm_copy = tp / "ROADMAP.md"
            cfg_txt = CONFIG_RS.read_text(encoding="utf-8").replace(
                "pub canonical_vcf_gt: bool,",
                "pub canonical_vcf_gt: bool,\n    pub future_v111_switch: bool,",
            )
            cfg_copy.write_text(cfg_txt, encoding="utf-8")
            rm_copy.write_text(ROADMAP_MD.read_text(encoding="utf-8") + "\n`future_v111_switch`\n", encoding="utf-8")
            grown = mod.parse_compat_profile_switches(cfg_copy, rm_copy, require_exact_count=False)
            self.assertEqual(len(grown), len(expected_switches) + 1)
            with self.assertRaises(RuntimeError):
                mod.parse_compat_profile_switches(cfg_copy, rm_copy, require_exact_count=True)

        for entry in catalog:
            self.assertIn(entry["domain"], mod.ALLOWED_DOMAINS)

    def test_audit_vcf_schema_bgzf_and_intersect_beds_unit(self) -> None:
        mod = _load_collector_module()
        vcf_hdr = (
            "##fileformat=VCFv4.2\n"
            "##ALT=<ID=DEL,Description=\"Deletion\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        )
        recs = [
            # 1. Clean SNV + clean <DEL>
            "chr1\t10\t.\tA\tG\t50\tPASS\tTYPE=SNV;VD=5;HICNT=5;HICOV=20\tGT:DP:VD:AD:AF:ALD\t0/1:20:5:15,5:0.25:3,2",
            "chr1\t100\t.\tC\t<DEL>\t60\tPASS\tEND=200;SVTYPE=DEL;SVLEN=-100;TYPE=DEL\tGT:DP:VD:AD:AF:ALD\t0/1:30:8:22,8:0.2667:4,4",
            # 2. Defective rows covering all 5 defect classes
            "chr1\t20\t.\tA\tT\t40\tPASS\tEND=20;TYPE=SNV\tGT:DP:VD:AD:AF:ALD\t0/1:10:4:6,4:0.4:2,2",
            "chr1\t30\t.\tG\t.\t0\tPASS\t.\tGT:DP:VD:AD:AF:ALD\t0/0:10:0:10,0:0.0:0,0",
            "chr1\t40\t.\tC\tA\t30\tf0.30\tTYPE=SNV\tGT:DP:VD:AD:AF:ALD\t0/0:20:3:17,3:0.15:2,1",
            "chr1\t50\t.\tACGT\t<dup-2>\t40\tPASS\tEND=60;SVTYPE=DUP;SVLEN=10\tGT:DP:VD:AD:AF:ALD\t0/1:20:5:15,5:0.25:3,2",
            "chr1\t60\t.\tT\tR\t-5\tPASS\tTYPE=Complex;VD=6;HICNT=8;HICOV=5\tGT:DP:VD:AD:AF:ALD\t0/1:4:6:0,6:1.5:3,2",
        ]
        n_clean, def_clean = mod.audit_vcf_schema(vcf_hdr + "\n".join(recs[:2]) + "\n", check_ref_ok=True)
        self.assertEqual(n_clean, 2)
        self.assertEqual(sum(def_clean.values()), 0)

        n_all, def_all = mod.audit_vcf_schema(vcf_hdr + "\n".join(recs) + "\n", check_ref_ok=True)
        self.assertEqual(n_all, 7)
        for k in mod.DEFECT_KEYS:
            self.assertEqual(def_all[k], 1, f"Expected exactly 1 defect in {k}, got {def_all}")

        # Out-of-range genotype allele index (e.g., 1/2 on biallelic n_alt=1) must be flagged
        bad_gt_rec = "chr1\t70\t.\tA\tG\t50\tPASS\tTYPE=SNV\tGT:DP:VD:AD:AF:ALD\t1/2:20:5:15,5:0.25:3,2\n"
        _, def_bad_gt = mod.audit_vcf_schema(vcf_hdr + bad_gt_rec, check_ref_ok=True)
        self.assertEqual(def_bad_gt["noncanonical_gt"], 1)

        # Test intersect_beds and count_bgzf_blocks on synthetic inputs
        with tempfile.TemporaryDirectory() as tmpdir:
            tp = Path(tmpdir)
            ba = tp / "a.bed"
            bb = tp / "b.bed"
            bo = tp / "out.bed"
            ba.write_text("chr1\t100\t300\nchr1\t250\t400\n", encoding="utf-8")
            bb.write_text("chr1\t150\t350\n", encoding="utf-8")
            n_iv, tot_bp = mod.intersect_beds(ba, bb, bo)
            self.assertEqual((n_iv, tot_bp), (1, 200))
            self.assertEqual(bo.read_text(encoding="utf-8"), "chr1\t150\t350\n")

            def make_bgzf_block(payload: bytes) -> bytes:
                comp = zlib.compressobj(6, zlib.DEFLATED, -15)
                cdata = comp.compress(payload) + comp.flush()
                bsize = 18 + len(cdata) + 8
                hdr = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bsize - 1)
                trl = struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))
                return hdr + cdata + trl

            bgz = tp / "test.vcf.gz"
            bgz.write_bytes(make_bgzf_block(b"hello\n") + make_bgzf_block(b"world\n") + make_bgzf_block(b""))
            self.assertEqual(mod.count_bgzf_blocks(bgz), 2)

    def test_stdlib_purity_of_collector_and_test(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        stdlib = sys.stdlib_module_names
        for py_path in (SCRIPT, TEST_FILE):
            tree = ast.parse(py_path.read_text(encoding="utf-8"), filename=str(py_path))
            for node in ast.walk(tree):
                mods: list[str] = []
                if isinstance(node, ast.Import):
                    mods.extend(alias.name for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                    mods.append(node.module)
                for mod in mods:
                    top = mod.split(".")[0]
                    self.assertIn(
                        top,
                        stdlib,
                        f"Non-stdlib import {mod!r} in {py_path.relative_to(ROOT)}",
                    )

    def test_gong_2022_crossref_doi_fixed_in_roadmap_and_upstream_toml(self) -> None:
        roadmap_txt = ROADMAP_MD.read_text(encoding="utf-8")
        upstream_txt = UPSTREAM_TOML.read_text(encoding="utf-8")
        upstream_doc = tomllib.loads(upstream_txt)

        self.assertIn("https://doi.org/10.1038/s41597-022-01359-6", roadmap_txt)
        self.assertNotIn("10.1038/s41597-022-01359-2", roadmap_txt)

        dois = upstream_doc["compat_off_and_cli_lineage"]["clinical_benchmark_dois"]
        self.assertIn("https://doi.org/10.1038/s41597-022-01359-6", dois)
        self.assertNotIn("10.1038/s41597-022-01359-2", upstream_txt)


if __name__ == "__main__":
    unittest.main()
