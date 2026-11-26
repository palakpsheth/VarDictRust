#!/usr/bin/env python3
"""Unit tests for paper/collect_v110_compat_metrics.py and benchmarks/v110-compat-study/."""

from __future__ import annotations

import ast
import csv
import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "paper" / "collect_v110_compat_metrics.py"
TEST_FILE = Path(__file__).resolve()
OUT_DIR = ROOT / "benchmarks" / "v110-compat-study"
CONFIG_RS = ROOT / "crates" / "vardictrust-core" / "src" / "config.rs"
ROADMAP_MD = ROOT / "ROADMAP.md"
UPSTREAM_TOML = ROOT / "UPSTREAM.toml"

ALLOWED_DOMAINS = frozenset(
    {
        "VCF Schema & Genotypes",
        "Structural Variants & Indexing",
        "Somatic & Amplicon Calling",
        "CIGAR, Pileup & Aligner",
        "Depth, Rounding & Determinism",
    }
)

DEFECT_KEYS = (
    "nonsv_info_end",
    "ad_number_r_mismatch",
    "noncanonical_gt",
    "noncanonical_sv",
    "accounting_type_iupac_shift3",
)


def parse_compat_profile_fields(config_rs_path: Path) -> list[str]:
    text = config_rs_path.read_text(encoding="utf-8")
    m = re.search(r"pub struct CompatProfile\s*\{(.*?)\n\}", text, flags=re.DOTALL)
    if not m:
        raise AssertionError(f"Could not locate pub struct CompatProfile in {config_rs_path}")
    body = m.group(1)
    return re.findall(r"^\s*pub\s+([a-z0-9_]+)\s*:\s*bool\s*,", body, flags=re.MULTILINE)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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
        for r in rows:
            target = OUT_DIR / r["path"]
            self.assertEqual(sha256_file(target), r["sha256"])
            self.assertEqual(str(target.stat().st_size), r["size_bytes"])

    def test_metrics_json_schema_and_six_arms_invariants(self) -> None:
        metrics_path = OUT_DIR / "v110_compat_metrics.json"
        arms_path = OUT_DIR / "v110_compat_arms.tsv"
        self.assertTrue(metrics_path.is_file(), f"Missing {metrics_path}")
        self.assertTrue(arms_path.is_file(), f"Missing {arms_path}")

        doc = json.loads(metrics_path.read_text(encoding="utf-8"))
        meta = doc["meta"]
        self.assertEqual(meta["study_id"], "v110-compat-and-presets")
        self.assertEqual(meta["study_version"], "1.1.0")
        self.assertEqual(meta["study_fingerprint"], "b18477ae0a7b07c1")
        self.assertEqual(meta["study_perf_fingerprint"], "5c41be8047133842")
        self.assertEqual(meta["compat_switch_count"], 32)
        self.assertEqual(meta["preset_count"], 4)
        self.assertEqual(meta["var2vcf_flag_count"], 18)
        self.assertEqual(meta["total_arms"], 6)
        for tool_ver_key in ("rtg_version", "bcftools_version", "tabix_version", "samtools_version"):
            self.assertTrue(meta[tool_ver_key].strip(), f"Empty {tool_ver_key}")

        # Verify zero host-dependent wall-clock/RSS metrics anywhere in JSON
        serialized = json.dumps(doc)
        self.assertIsNone(
            re.search(r'"(wall_sec|elapsed_s|rss_kb|max_rss|runtime_s)"', serialized),
            "Host-dependent timing/RSS key leaked into v110_compat_metrics.json",
        )

        # Arm 1
        a1 = doc["arms"]["arm1_vcf_schema_and_bcftools_norm"]
        self.assertEqual(a1["total_callsets"], 19)
        self.assertGreaterEqual(a1["distinct_bams"], 14)
        for k in DEFECT_KEYS:
            self.assertGreater(a1["vardictjava_defects"][k], 0, f"Expected VJ defect > 0 for {k}")
            self.assertEqual(a1["compat_off_defects"][k], 0, f"Expected compat_off defect == 0 for {k}")
        self.assertEqual(a1["total_defects_vj"], sum(a1["vardictjava_defects"].values()))
        self.assertEqual(a1["total_defects_off"], 0)
        self.assertEqual(a1["bcftools_norm_check_ref_off_passed"], 19)
        self.assertEqual(a1["bcftools_norm_check_ref_off_total"], 19)
        self.assertLess(a1["bcftools_norm_check_ref_vj_passed"], 19)

        # Arm 2
        a2 = doc["arms"]["arm2_seqc2_diploid_vcfeval"]
        self.assertEqual(a2["discovery_af"], 0.05)
        self.assertEqual(a2["vcf_af"], 0.30)
        self.assertEqual(a2["gt_proj_tvaf_threshold"], 0.75)
        self.assertGreater(a2["vardictjava"]["gt00_carrier_somatic"], 0)
        self.assertEqual(a2["compat_off"]["gt00_carrier_somatic"], 0)
        self.assertGreater(a2["tp_baseline_gain"], 0)
        self.assertGreater(a2["f_measure_gain"], 0.0)
        self.assertGreater(a2["best_f_measure_gain"], 0.0)
        self.assertGreater(a2["compat_off"]["tp_baseline"], a2["vardictjava"]["tp_baseline"])
        self.assertGreater(a2["compat_off"]["f_measure"], a2["vardictjava"]["f_measure"])
        self.assertGreater(a2["compat_off"]["best_f_measure"], a2["vardictjava"]["best_f_measure"])

        # Arm 3
        a3 = doc["arms"]["arm3_somatic_multiallelic_and_truth"]
        self.assertEqual(a3["distinct_bam_pairs"], 5)
        self.assertEqual(a3["total_runs"], 6)
        self.assertEqual(a3["vj_recovered_pairs"], 0)
        self.assertEqual(a3["off_recovered_pairs"], 5)
        self.assertEqual(a3["vj_recovered_runs"], 0)
        self.assertEqual(a3["off_recovered_runs"], 6)
        self.assertGreaterEqual(a3["muc19_variants"], 20)
        self.assertEqual(a3["muc19_vd_gt_dp_violations"], 0)
        for cohort_key in ("seqc2_hc", "hg008_v04"):
            self.assertIn("vardictjava", a3[cohort_key])
            self.assertIn("compat_off", a3[cohort_key])
            self.assertGreater(a3[cohort_key]["compat_off"]["tp"], 0)

        # Arm 4
        a4 = doc["arms"]["arm4_brca2_adaptor_chimera"]
        self.assertEqual(a4["chimera_pos"], 32683955)
        self.assertGreaterEqual(a4["vj_chimeric_insertions"], 1)
        self.assertGreater(a4["vj_chimeric_ins_bp"], 30)
        self.assertEqual(a4["off_chimeric_insertions"], 0)

        # Arm 5
        a5 = doc["arms"]["arm5_clinical_presets_and_amplicon"]
        self.assertGreaterEqual(a5["wgs_hg002_60kbp_records"], 20)
        self.assertGreaterEqual(a5["wgs_hg008_paired_records"], 10)
        self.assertGreaterEqual(a5["wes_seqc2_264kbp_records"], 15)
        self.assertTrue(a5["wes_bed8_non_hijack_verified"])
        self.assertGreaterEqual(a5["amplicon_single_8col_records"], 1)
        self.assertFalse(a5["paired_somatic_amplicon_vj_supported"])
        self.assertGreaterEqual(a5["paired_somatic_amplicon_off_records"], 10)
        self.assertGreater(a5["qiaseq_spe_known_pos_tp"], a5["qiaseq_spe_no_sp_known_pos_tp"])
        self.assertGreaterEqual(a5["qiaseq_spe_known_pos_total"], a5["qiaseq_spe_known_pos_tp"])
        self.assertTrue(a5["spikein_detected"])
        self.assertEqual(a5["spikein_observed_vd"], 6)
        self.assertGreaterEqual(a5["spikein_observed_dp"], 1000)

        # Arm 6
        a6 = doc["arms"]["arm6_csi_spanning_sv_and_pileup_tiling"]
        self.assertEqual(a6["csi_leaf_bin_bp"], 16384)
        self.assertEqual(a6["max_explicit_sv_seq_bp"], 32768)
        self.assertGreaterEqual(a6["multi_block_bgzf_blocks"], 20)
        self.assertGreater(a6["interior_bin_queries_total"], 0)
        self.assertEqual(a6["interior_bin_queries_passed"], a6["interior_bin_queries_total"])
        self.assertEqual(a6["pos0_sentinel_rows"], 0)

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
        ub_md = (ROOT / "docs" / "upstream-bugs.md").read_text(encoding="utf-8")
        self.assertIn("https://doi.org/10.46471/gigabyte.194", ub_md)
        self.assertNotIn("10.1093/gigabyte/giaf156", ub_md)

    def test_switch_catalog_matches_config_rs_in_order_and_hygiene(self) -> None:
        expected_switches = parse_compat_profile_fields(CONFIG_RS)
        self.assertEqual(len(expected_switches), 32)

        metrics_path = OUT_DIR / "v110_compat_metrics.json"
        self.assertTrue(metrics_path.is_file(), f"Missing {metrics_path}")
        doc = json.loads(metrics_path.read_text(encoding="utf-8"))
        catalog = doc["switch_catalog"]
        self.assertEqual(len(catalog), 35)

        got_switches = [entry["name"] for entry in catalog[:32]]
        self.assertEqual(got_switches, expected_switches)

        cli_names = [entry["name"] for entry in catalog[32:]]
        self.assertEqual(
            cli_names,
            ["--vcf-all-alleles", "--amplicon-single-primer", "--preset"],
        )

        bare_issue_re = re.compile(r"(?<![A-Za-z0-9_/-])#\d+\b")
        forbidden_jargon_re = re.compile(
            r"\b(Gate\s*\d+|s[1-6]_[a-z0-9_]+|arm[1-6]_[a-z0-9_]+|Option C|Tier[- ]P|evidence_certified)\b"
        )
        for idx, entry in enumerate(catalog, start=1):
            self.assertEqual(entry["index"], idx)
            self.assertIn(entry["domain"], ALLOWED_DOMAINS)
            for field in (
                "name",
                "category",
                "domain",
                "upstream_refs",
                "vardictjava_behavior",
                "compat_off_behavior",
                "verification_arm",
            ):
                val = str(entry[field])
                self.assertTrue(val.strip(), f"Empty {field} in switch_catalog[{idx}]")
                self.assertIsNone(
                    bare_issue_re.search(val),
                    f"Bare #N citation in switch_catalog[{idx}].{field}: {val!r}",
                )
                self.assertIsNone(
                    forbidden_jargon_re.search(val),
                    f"Internal harness jargon in switch_catalog[{idx}].{field}: {val!r}",
                )
                # Verify no unescaped '<' or '|' outside inline backtick code spans
                outside_code = re.sub(r"`[^`]+`", "", val)
                self.assertNotIn(
                    "<",
                    outside_code,
                    f"Unbackticked '<' token in switch_catalog[{idx}].{field}: {val!r}",
                )
                self.assertNotIn(
                    "|",
                    outside_code,
                    f"Unescaped '|' pipe in switch_catalog[{idx}].{field}: {val!r}",
                )

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
