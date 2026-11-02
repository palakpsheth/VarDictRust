#!/usr/bin/env python3
"""
harness/tests/test_evidence_inputs.py

Verifies that:
1. `harness/record_evidence.py` records `inputs` (`path`, `size_bytes`, `sha256`)
   for explicit `--input` flags and file arguments in `cmd` (including `a.bam|b.bam`
   and companion `.fai`/`.bai` indices).
2. `bash harness/recertify.sh --dry-run` reports BED digests and exits 0 when
   staged inputs match committed evidence.
3. `bash harness/recertify.sh --dry-run` refuses to proceed (exits non-zero) when
   a staged input's SHA-256 differs from committed evidence unless
   `RECERTIFY_NEW_INPUTS=1` is set.
"""

from __future__ import annotations

import ast
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent
if str(ROOT_DIR / "harness") not in sys.path:
    sys.path.insert(0, str(ROOT_DIR / "harness"))

import perf_budget  # noqa: E402
import sync_evidence_docs  # noqa: E402


class TestEvidenceInputs(unittest.TestCase):
    def test_record_evidence_tracks_input_sha256_and_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            tpath = pathlib.Path(tmp)
            ref = tpath / "ref.fa"
            fai = tpath / "ref.fa.fai"
            bam1 = tpath / "tumor.bam"
            bai1 = tpath / "tumor.bam.bai"
            bam2 = tpath / "normal.bam"
            bed = tpath / "regions.bed"
            extra_in = tpath / "custom.txt"
            out_dir = tpath / "evidence"

            ref.write_bytes(b">chr1\nACGTACGT\n")
            fai.write_bytes(b"chr1\t8\t6\t8\t9\n")
            bam1.write_bytes(b"BAM1_PAYLOAD")
            bai1.write_bytes(b"BAI1_INDEX")
            bam2.write_bytes(b"BAM2_PAYLOAD")
            bed.write_bytes(b"chr1\t0\t8\treg1\n")
            extra_in.write_bytes(b"EXTRA_INPUT")

            cmd = [
                sys.executable,
                str(ROOT_DIR / "harness" / "record_evidence.py"),
                "--out-dir",
                str(out_dir),
                "--input",
                str(extra_in),
                "unit-input-test",
                "MATCH",
                "-",
                "--",
                "bash",
                "harness/compare.sh",
                str(ref),
                f"{bam1}|{bam2}",
                str(bed),
                "SAMPLE",
                "-th",
                "4",
            ]
            res = subprocess.run(cmd, cwd=ROOT_DIR, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, f"record_evidence.py failed: {res.stderr}")

            ev_file = out_dir / "unit-input-test.json"
            self.assertTrue(ev_file.is_file(), "evidence JSON was not created")
            data = json.loads(ev_file.read_text(encoding="utf-8"))
            self.assertIn("inputs", data, "record_evidence.py did not record 'inputs' list")

            by_path = {entry["path"]: entry for entry in data["inputs"]}
            for expected_file in (ref, fai, bam1, bai1, bam2, bed, extra_in):
                p_str = str(expected_file)
                self.assertIn(p_str, by_path, f"Missing {p_str} in recorded inputs: {by_path}")
                raw = expected_file.read_bytes()
                self.assertEqual(by_path[p_str]["size_bytes"], len(raw))
                self.assertEqual(
                    by_path[p_str]["sha256"],
                    hashlib.sha256(raw).hexdigest(),
                )

    def test_recertify_dry_run_reports_bed_digest_and_guards_mismatch(self):
        recertify_text = (ROOT_DIR / "harness" / "recertify.sh").read_text(encoding="utf-8")
        self.assertIn("--dry-run", recertify_text, "recertify.sh does not support --dry-run")
        self.assertIn("RECERTIFY_NEW_INPUTS", recertify_text, "recertify.sh does not check RECERTIFY_NEW_INPUTS")

        res = subprocess.run(
            ["bash", "harness/recertify.sh", "--dry-run"],
            cwd=ROOT_DIR,
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(
            res.returncode,
            0,
            f"recertify.sh --dry-run failed:\nstdout={res.stdout}\nstderr={res.stderr}",
        )
        self.assertRegex(
            res.stdout,
            r"df8cc38f4a269c0923120abf1d4f4429fc2fc92e55e6a5eda23468b5ec6e1bb1|[0-9a-f]{64}",
            "recertify.sh --dry-run did not report a 64-hex SHA-256 BED digest",
        )

        # Verify that a mismatched input digest in evidence is rejected unless RECERTIFY_NEW_INPUTS=1
        with tempfile.TemporaryDirectory() as tmp:
            tpath = pathlib.Path(tmp)
            staged_bed = tpath / "eval.bed"
            staged_bed.write_bytes(b"chr20\t100\t200\tamplicon_1\n")
            ev_dir = tpath / "evidence"
            ev_dir.mkdir()
            (ev_dir / "mock-tier2.json").write_text(
                json.dumps(
                    {
                        "name": "mock-tier2",
                        "result": "MATCH",
                        "fingerprint": "0000000000000000",
                        "inputs": [
                            {
                                "path": str(staged_bed),
                                "size_bytes": 24,
                                "sha256": "0" * 64,
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            env_fail = os.environ.copy()
            env_fail["RECERTIFY_EVIDENCE_DIR"] = str(ev_dir)
            env_fail.pop("RECERTIFY_NEW_INPUTS", None)
            res_fail = subprocess.run(
                ["bash", "harness/recertify.sh", "--dry-run"],
                cwd=ROOT_DIR,
                env=env_fail,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(
                res_fail.returncode,
                0,
                "recertify.sh --dry-run should fail when staged input sha256 mismatches committed evidence",
            )
            self.assertIn("RECERTIFY_NEW_INPUTS=1", res_fail.stderr + res_fail.stdout)

            env_ok = env_fail.copy()
            env_ok["RECERTIFY_NEW_INPUTS"] = "1"
            res_ok = subprocess.run(
                ["bash", "harness/recertify.sh", "--dry-run"],
                cwd=ROOT_DIR,
                env=env_ok,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                res_ok.returncode,
                0,
                f"recertify.sh --dry-run with RECERTIFY_NEW_INPUTS=1 should succeed: {res_ok.stderr}",
            )

            # Verify that needs_stamp invalidates receipts when RECERTIFY_NEW_INPUTS=1 or input digest changed
            (ev_dir / "amplicon-tier2.json").write_text(
                json.dumps(
                    {
                        "name": "amplicon-tier2",
                        "result": "MATCH",
                        "fingerprint": "same-fp",
                        "inputs": [
                            {
                                "path": str(staged_bed),
                                "size_bytes": len(staged_bed.read_bytes()),
                                "sha256": hashlib.sha256(staged_bed.read_bytes()).hexdigest(),
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            extract_needs_stamp = (
                'FP="same-fp" RECERTIFY_EVIDENCE_DIR="$1" bash -c '
                '"$(sed -n \'/^needs_stamp()/,/^}/p\' harness/recertify.sh); needs_stamp amplicon-tier2"'
            )
            res_reuse = subprocess.run(
                ["bash", "-c", extract_needs_stamp, "_", str(ev_dir)],
                cwd=ROOT_DIR,
                env={k: v for k, v in os.environ.items() if k != "RECERTIFY_NEW_INPUTS"},
            )
            self.assertEqual(res_reuse.returncode, 1, "unchanged input + matching FP should reuse receipt (exit 1)")

            res_override = subprocess.run(
                ["bash", "-c", extract_needs_stamp, "_", str(ev_dir)],
                cwd=ROOT_DIR,
                env={**os.environ, "RECERTIFY_NEW_INPUTS": "1"},
            )
            self.assertEqual(res_override.returncode, 0, "RECERTIFY_NEW_INPUTS=1 must force needs_stamp rerun (exit 0)")

            staged_bed.write_bytes(b"chr20\t100\t200\tamplicon_modified\n")
            res_changed = subprocess.run(
                ["bash", "-c", extract_needs_stamp, "_", str(ev_dir)],
                cwd=ROOT_DIR,
                env={k: v for k, v in os.environ.items() if k != "RECERTIFY_NEW_INPUTS"},
            )
            self.assertEqual(res_changed.returncode, 0, "changed input digest must force needs_stamp rerun (exit 0)")

    def test_sync_evidence_docs_updates_tier2_tables_and_summary_ranges(self):
        with tempfile.TemporaryDirectory() as tmp:
            tpath = pathlib.Path(tmp)
            rel_files = [
                "Cargo.toml",
                "README.md",
                "DIVERGENCES.md",
                "benchmarks/README.md",
                "docs/validation/omnibus-empirical-validation.md",
                "docs/testing/test-matrix-and-ci-audit.md",
                ".agents/skills/vardictrust-release/references/version-surface-map.md",
            ]
            for rel in rel_files:
                src = ROOT_DIR / rel
                dst = tpath / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

            ev_dst = tpath / "benchmarks" / "evidence"
            ev_dst.mkdir(parents=True, exist_ok=True)
            for jf in (ROOT_DIR / "benchmarks" / "evidence").glob("*.json"):
                shutil.copy2(jf, ev_dst / jf.name)

            sync_script = str(ROOT_DIR / "harness" / "sync_evidence_docs.py")
            res_init = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath), "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                res_init.returncode,
                0,
                f"Initial --check failed: {res_init.stderr}",
            )

            # Mutate Tier-2 and CRAM-twin timings/RSS in temp copy with deterministic synthetic values
            synthetic_tier2 = {
                "hg008-tier2.json": {
                    "java_time": {"wall": "1:47.00", "max_rss_kb": "3190024"},
                    "cpp_time": {"wall": "0:21.96", "max_rss_kb": "310908"},
                    "rust_time": {"wall": "0:38.20", "max_rss_kb": "531860"},
                },
                "seqc2-slices-tier2.json": {
                    "java_time": {"wall": "0:19.04", "max_rss_kb": "2545556"},
                    "cpp_time": {"wall": "0:14.79", "max_rss_kb": "304176"},
                    "rust_time": {"wall": "0:01.36", "max_rss_kb": "154168"},
                },
                "seqc2-exome-tier2.json": {
                    "java_time": {"wall": "17:27.99", "max_rss_kb": "8419392"},
                    "cpp_time": {"wall": "22:25.73", "max_rss_kb": "1762332"},
                    "rust_time": {"wall": "5:13.86", "max_rss_kb": "628712"},
                },
                "hg002-tier2.json": {
                    "java_time": {"wall": "2:43.25", "max_rss_kb": "7985080"},
                    "cpp_time": {"wall": "1:18.96", "max_rss_kb": "1333816"},
                    "rust_time": {"wall": "0:23.09", "max_rss_kb": "529500"},
                },
            }
            for fname, overrides in synthetic_tier2.items():
                fpath = ev_dst / fname
                fdata = json.loads(fpath.read_text(encoding="utf-8"))
                for k, v in overrides.items():
                    fdata[k].update(v)
                fpath.write_text(
                    json.dumps(fdata, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )

            cram_path = ev_dst / "hg008-cram-twin.json"
            cram_data = json.loads(cram_path.read_text(encoding="utf-8"))
            cram_data["rust_time"]["wall"] = "0:50.00"
            cram_data["rust_time"]["max_rss_kb"] = "256000"
            cram_data["rust_cram_time"]["wall"] = "1:15.00"
            cram_data["rust_cram_time"]["max_rss_kb"] = "1572864"
            cram_path.write_text(
                json.dumps(cram_data, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            scale_path = ev_dst / "scale.json"
            scale_data = json.loads(scale_path.read_text(encoding="utf-8"))
            scale_data["extra"]["wall_s"] = {
                "1": 450.0,
                "2": 225.0,
                "4": 150.0,
                "8": 138.46,
            }
            scale_data["extra"]["maxrss_kb"] = {
                "1000000": 524288,
                "10000000": 1258291,
                "100000000": 1310720,
                "248956422": 1447035,
            }
            scale_data["extra"]["speedup_8"] = round(450.0 / 138.46, 6)
            scale_data["extra"]["rss_ratio_chr1_vs_10mb"] = 1447035 / 1258291
            scale_data["extra"]["cpp_wall_s"] = {
                "1": 600.0,
                "2": 350.0,
                "4": 210.0,
                "8": 195.5,
            }
            scale_data["extra"]["cpp_maxrss_kb"] = {
                "1000000": 220000,
                "10000000": 2900000,
                "100000000": 22000000,
                "248956422": 26214400,
            }
            scale_data["extra"]["perf_fingerprint"] = "fedcba9876543210"
            scale_data["extra"]["cpus"] = "0-3,8-11"
            scale_data["extra"]["P"] = 4
            scale_path.write_text(
                json.dumps(scale_data, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            # --check must detect drift across all 6 surfaces
            res_drift = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath), "--check"],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res_drift.returncode, 0, "--check should fail when timings drift")
            for expected_surface in (
                "README.md",
                "DIVERGENCES.md",
                "benchmarks/README.md",
                "docs/validation/omnibus-empirical-validation.md",
                "docs/testing/test-matrix-and-ci-audit.md",
                ".agents/skills/vardictrust-release/references/version-surface-map.md",
            ):
                self.assertIn(expected_surface, res_drift.stderr)

            # Sync and verify updated values
            res_sync = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath)],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_sync.returncode, 0, f"Sync failed: {res_sync.stderr}")

            readme_text = (tpath / "README.md").read_text(encoding="utf-8")
            self.assertIn("**2.8×–14.0× faster on WGS/WES**", readme_text)
            self.assertIn("in **6.0×–16.5× less memory**", readme_text)
            self.assertIn("**2.8×–14.0× faster** on WGS/WES", readme_text)
            self.assertIn("**6.0×–16.5× less RAM** than VarDictJava on WGS/WES", readme_text)
            self.assertIn(
                "| **GIAB HG008 Somatic PDAC** | `4` (`-th 4`²) | `1:47.00` | `0:21.96` ([`CPP-9`](docs/upstream-bugs.md#cpp-9)) | **`0:38.20`** | **2.8× faster** | `3.04 GiB` | `303.6 MiB` | **`519.4 MiB`** | **6.0× less** | [`hg008-tier2.json`](benchmarks/evidence/hg008-tier2.json) |",
                readme_text,
            )
            self.assertIn(
                "| **GIAB HG008 CRAM-vs-BAM Twin** | `4` (`-th 4`²) | `0:50.00` (Rust BAM) | — | **`1:15.00` (Rust CRAM)** | **1.50× CRAM/BAM** | `250.0 MiB` (BAM) | — | **`1.50 GiB` (CRAM)** | **Bounded `512 MiB` cache** | [`hg008-cram-twin.json`](benchmarks/evidence/hg008-cram-twin.json) |",
                readme_text,
            )

            bench_text = (tpath / "benchmarks" / "README.md").read_text(encoding="utf-8")
            self.assertIn(
                "| `1:47.00` | `0:21.96` | **`0:38.20`** | **2.8× faster** | `3.04 GiB` | `303.6 MiB` | **`519.4 MiB`** | **6.0× less** | **[MISMATCH (`CPP-9`)](../docs/upstream-bugs.md#cpp-9)** | [`hg008-tier2.json`](evidence/hg008-tier2.json) |",
                bench_text,
            )
            self.assertIn(
                "| `0:50.00` (Rust BAM) | — | **`1:15.00` (Rust CRAM)** | **1.50× CRAM/BAM** | `250.0 MiB` (BAM) | — | **`1.50 GiB` (CRAM)** | **Bounded `512 MiB` cache** |",
                bench_text,
            )
            self.assertIn("$S_8 = 3.25\\times$, RSS 1.15×", bench_text)
            self.assertIn("pinned physical cores `0-3,8-11` for `scale.json`", bench_text)
            self.assertIn("(`perf_fingerprint: fedcba9876543210`)", bench_text)
            self.assertIn("| 1 core | 7m 30s (450.00 s) | 1.00× (baseline) |", bench_text)
            self.assertIn("| 2 cores | 3m 45s (225.00 s) | 2.00× |", bench_text)
            self.assertIn("| 4 cores | 2m 30s (150.00 s) | 3.00× |", bench_text)
            self.assertIn("| 8 threads | 2m 18s (138.46 s) | **3.25× parallel** |", bench_text)
            self.assertIn("| 1 Mbp | 512.0 MiB (524,288 KB) |", bench_text)
            self.assertIn("| 249 Mbp (Full chr1) | **1.38 GiB (1,447,035 KB)** |", bench_text)
            self.assertIn("**Observed RSS Ratio (249 Mbp vs 10 Mbp): 1.15× (`1.150×`)**", bench_text)
            self.assertIn(
                "at `-th 8` BED wall 195.50 s, full-span (249 Mbp) unchunked `-R` RSS 25.00 GiB (26,214,400 KB",
                bench_text,
            )

            omnibus_text = (tpath / "docs" / "validation" / "omnibus-empirical-validation.md").read_text(encoding="utf-8")
            self.assertIn("**3.25× parallel**", omnibus_text)
            self.assertIn("**1.15× flat RSS scaling**", omnibus_text)
            self.assertIn("**$1.15\\times$ (`1.150\\times`)**", omnibus_text)
            self.assertIn(
                "**3.25×** parallel speedup at 8 threads on chromosome 1 (3-run median) with flat memory scaling (**1.15×** (`1.150×`) RSS ratio: 1.20 GiB at 10 Mb vs 1.38 GiB at 249 Mb).",
                omnibus_text,
            )
            self.assertIn(
                "| **3.25× speedup, 1.15× flat RSS** | advisory 25.00 GiB unchunked `-R` RSS / 195.50 s BED wall ([CPP-5](../upstream-bugs.md#cpp-5)) |",
                omnibus_text,
            )
            self.assertIn(
                "| `1:47.00` | `0:21.96` | **`0:38.20`** | **2.8× faster** | `3.04 GiB` | `303.6 MiB` | **`519.4 MiB`** | **6.0× less** | [`hg008-tier2.json`](../../benchmarks/evidence/hg008-tier2.json) |",
                omnibus_text,
            )
            self.assertIn(
                "| `0:50.00` (Rust BAM) | — | **`1:15.00` (Rust CRAM)** | **1.50× CRAM/BAM** | `250.0 MiB` (BAM) | — | **`1.50 GiB` (CRAM)** | **Bounded `512 MiB` cache** | [`hg008-cram-twin.json`](../../benchmarks/evidence/hg008-cram-twin.json) |",
                omnibus_text,
            )

            matrix_text = (tpath / "docs" / "testing" / "test-matrix-and-ci-audit.md").read_text(encoding="utf-8")
            self.assertIn("(`cram_to_bam_wall_ratio <= 2.50`, achieved `1.50×`)", matrix_text)
            self.assertIn(
                "(`SCALE_CPUS=0-3,8-11` / `taskset -c 0-3,8-11`, $P = 4$ physical cores under pin).",
                matrix_text,
            )

            # --check must now exit 0
            res_after = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath), "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_after.returncode, 0, f"--check after sync failed: {res_after.stderr}")

    def test_python38_future_annotations_present(self):
        for rel in (
            "harness/sync_evidence_docs.py",
            "harness/perf_budget.py",
            "harness/tests/test_evidence_inputs.py",
        ):
            fpath = ROOT_DIR / rel
            tree = ast.parse(fpath.read_text(encoding="utf-8"), filename=str(fpath))
            imports = [
                node for node in tree.body if isinstance(node, (ast.Import, ast.ImportFrom))
            ]
            self.assertTrue(imports, f"No top-level imports found in {rel}")
            first = imports[0]
            self.assertIsInstance(first, ast.ImportFrom, f"First import in {rel} is not ImportFrom")
            self.assertEqual(first.module, "__future__", f"First import in {rel} is not from __future__")
            self.assertEqual(
                [alias.name for alias in first.names],
                ["annotations"],
                f"First import in {rel} is not 'from __future__ import annotations'",
            )

    def test_sync_evidence_docs_handles_null_cpp_metrics_and_roundtrips(self):
        with tempfile.TemporaryDirectory() as tmp:
            tpath = pathlib.Path(tmp)
            rel_files = [
                "Cargo.toml",
                "README.md",
                "DIVERGENCES.md",
                "benchmarks/README.md",
                "docs/validation/omnibus-empirical-validation.md",
                "docs/testing/test-matrix-and-ci-audit.md",
                ".agents/skills/vardictrust-release/references/version-surface-map.md",
            ]
            for rel in rel_files:
                src = ROOT_DIR / rel
                dst = tpath / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

            ev_dst = tpath / "benchmarks" / "evidence"
            ev_dst.mkdir(parents=True, exist_ok=True)
            for jf in (ROOT_DIR / "benchmarks" / "evidence").glob("*.json"):
                shutil.copy2(jf, ev_dst / jf.name)

            sync_script = str(ROOT_DIR / "harness" / "sync_evidence_docs.py")
            # Ensure initial state in tpath is synced with live evidence
            res_pre = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath)],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_pre.returncode, 0, f"Pre-sync failed: {res_pre.stderr}")

            md_files = [r for r in rel_files if r.endswith(".md")]
            initial_contents = {
                rel: (tpath / rel).read_text(encoding="utf-8") for rel in md_files
            }

            # Null out cpp_time on all 5 Tier-2 JSON files and cpp_wall_s/cpp_maxrss_kb on scale.json
            for t2 in (
                "amplicon-tier2.json",
                "hg002-tier2.json",
                "seqc2-exome-tier2.json",
                "seqc2-slices-tier2.json",
                "hg008-tier2.json",
            ):
                p = ev_dst / t2
                d = json.loads(p.read_text(encoding="utf-8"))
                d["cpp_time"] = None
                p.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")

            scale_p = ev_dst / "scale.json"
            scale_d = json.loads(scale_p.read_text(encoding="utf-8"))
            scale_d["extra"]["cpp_wall_s"] = {"1": None, "2": None, "4": None, "8": None}
            scale_d["extra"]["cpp_maxrss_kb"] = {
                "1000000": None,
                "10000000": None,
                "100000000": None,
                "248956422": None,
            }
            scale_p.write_text(json.dumps(scale_d, indent=2, sort_keys=True) + "\n", encoding="utf-8")

            res_null = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath)],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_null.returncode, 0, f"Null-C++ sync failed: {res_null.stderr}")
            res_null_check = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath), "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                res_null_check.returncode, 0, f"Null-C++ --check failed: {res_null_check.stderr}"
            )

            readme_null = (tpath / "README.md").read_text(encoding="utf-8")
            self.assertIn("| **GIAB HG008 Somatic PDAC** | `4` (`-th 4`²) | `1:47.00` | — | **`0:16.84`** |", readme_null)
            self.assertIn("and **— less RAM** than unchunked `vardictcpp` (`—`) on whole-chromosome 1.", readme_null)
            self.assertIn("| `7:15.96` (Rust `-th 1`³; Java `OOM`) | — ([`CPP-5`](docs/upstream-bugs.md#cpp-5)) | **`2:18.92`** |", readme_null)

            bench_null = (tpath / "benchmarks" / "README.md").read_text(encoding="utf-8")
            self.assertIn("| **GIAB HG008 Somatic PDAC** | SomaticMode (`tumor\\|normal`, `-th 4`²) | `1:47.00` | — | **`0:16.84`** |", bench_null)
            self.assertIn("at `-th 8` BED wall —, full-span (249 Mbp) unchunked `-R` RSS — (per [CPP-5](../docs/upstream-bugs.md#cpp-5)).", bench_null)

            omnibus_null = (tpath / "docs" / "validation" / "omnibus-empirical-validation.md").read_text(encoding="utf-8")
            self.assertIn("| advisory — ([CPP-5](../upstream-bugs.md#cpp-5)) |", omnibus_null)
            self.assertIn("| `7:15.96` (1 core Rust) | — | **`2:18.92` (8 cores)** |", omnibus_null)

            # Restore original JSON files and verify 100% lossless byte-identical round-trip
            for jf in (ROOT_DIR / "benchmarks" / "evidence").glob("*.json"):
                shutil.copy2(jf, ev_dst / jf.name)

            res_restore = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath)],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_restore.returncode, 0, f"Restore sync failed: {res_restore.stderr}")
            res_restore_check = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath), "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                res_restore_check.returncode, 0, f"Restore --check failed: {res_restore_check.stderr}"
            )

            for rel in md_files:
                restored = (tpath / rel).read_text(encoding="utf-8")
                self.assertEqual(
                    restored,
                    initial_contents[rel],
                    f"Lossless round-trip failed for {rel}",
                )

    def test_parse_wall_s_and_invalid_numeric_evidence_rejected(self):
        bad_walls = [
            "1:2:3:4",
            "0:00.00",
            "-0.0",
            "1.0:10.00",
            "1:60.00",
            "1:60:00",
            "-1:00",
            "nan",
            "inf",
            "",
            True,
        ]
        for bad in bad_walls:
            with self.subTest(func="parse_wall_s", val=bad):
                with self.assertRaises(ValueError):
                    sync_evidence_docs.parse_wall_s(bad)
            with self.subTest(func="perf_budget.secs", val=bad):
                with self.assertRaises(ValueError):
                    perf_budget.secs(bad)

        for bad_int in (3.7, True, 0, -5, "12.3"):
            with self.subTest(func="_require_pos_int", val=bad_int):
                with self.assertRaises(ValueError):
                    sync_evidence_docs._require_pos_int(bad_int, "test_int")

        self.assertEqual(sync_evidence_docs.format_wall_human_min_sec(59.9), "1m 00s")

        with tempfile.TemporaryDirectory() as tmp:
            tpath = pathlib.Path(tmp)
            rel_files = [
                "Cargo.toml",
                "README.md",
                "DIVERGENCES.md",
                "benchmarks/README.md",
                "docs/validation/omnibus-empirical-validation.md",
                "docs/testing/test-matrix-and-ci-audit.md",
                ".agents/skills/vardictrust-release/references/version-surface-map.md",
            ]
            for rel in rel_files:
                src = ROOT_DIR / rel
                dst = tpath / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

            ev_dst = tpath / "benchmarks" / "evidence"
            ev_dst.mkdir(parents=True, exist_ok=True)
            for jf in (ROOT_DIR / "benchmarks" / "evidence").glob("*.json"):
                shutil.copy2(jf, ev_dst / jf.name)

            sync_script = str(ROOT_DIR / "harness" / "sync_evidence_docs.py")

            # Case 1: scale.json extra.wall_s is missing, has "8": None, or has "8": 0.0
            scale_p = ev_dst / "scale.json"
            orig_scale = scale_p.read_text(encoding="utf-8")
            for mut_scale in (
                lambda d: d["extra"].pop("wall_s", None),
                lambda d: d["extra"]["wall_s"].__setitem__("8", None),
                lambda d: d["extra"]["wall_s"].__setitem__("8", 0.0),
            ):
                scale_d = json.loads(orig_scale)
                mut_scale(scale_d)
                scale_p.write_text(json.dumps(scale_d), encoding="utf-8")
                res1 = subprocess.run(
                    [sys.executable, sync_script, "--repo-root", str(tpath)],
                    capture_output=True,
                    text=True,
                )
                self.assertNotEqual(res1.returncode, 0)
                self.assertIn("ERROR:", res1.stderr)
            scale_p.write_text(orig_scale, encoding="utf-8")

            # Case 2: hg008-tier2.json rust_time.wall is "0:00.00" or rust_time.max_rss_kb is 0
            hg008_p = ev_dst / "hg008-tier2.json"
            orig_hg008 = hg008_p.read_text(encoding="utf-8")
            for field, bad_v in (("wall", "0:00.00"), ("max_rss_kb", 0)):
                hg008_d = json.loads(orig_hg008)
                hg008_d["rust_time"][field] = bad_v
                hg008_p.write_text(json.dumps(hg008_d), encoding="utf-8")
                res2 = subprocess.run(
                    [sys.executable, sync_script, "--repo-root", str(tpath)],
                    capture_output=True,
                    text=True,
                )
                self.assertNotEqual(res2.returncode, 0)
                self.assertIn("ERROR:", res2.stderr)
            hg008_p.write_text(orig_hg008, encoding="utf-8")

    def test_sync_evidence_docs_late_pattern_error_leaves_all_files_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            tpath = pathlib.Path(tmp)
            rel_files = [
                "Cargo.toml",
                "README.md",
                "DIVERGENCES.md",
                "benchmarks/README.md",
                "docs/validation/omnibus-empirical-validation.md",
                "docs/testing/test-matrix-and-ci-audit.md",
                ".agents/skills/vardictrust-release/references/version-surface-map.md",
            ]
            for rel in rel_files:
                src = ROOT_DIR / rel
                dst = tpath / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)

            ev_dst = tpath / "benchmarks" / "evidence"
            ev_dst.mkdir(parents=True, exist_ok=True)
            for jf in (ROOT_DIR / "benchmarks" / "evidence").glob("*.json"):
                shutil.copy2(jf, ev_dst / jf.name)

            # Mutate hg008-tier2.json so earlier markdown files would be staged for modification
            hg008_p = ev_dst / "hg008-tier2.json"
            hg008_d = json.loads(hg008_p.read_text(encoding="utf-8"))
            hg008_d["rust_time"]["wall"] = "0:38.20"
            hg008_p.write_text(
                json.dumps(hg008_d, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )

            # Corrupt the last surface (.agents/skills/vardictrust-release/references/version-surface-map.md)
            smap_p = tpath / ".agents/skills/vardictrust-release/references/version-surface-map.md"
            smap_p.write_text("# Corrupted surface map with no matching table rows\n", encoding="utf-8")

            md_files = [r for r in rel_files if r.endswith(".md")]
            digests_before = {
                rel: hashlib.sha256((tpath / rel).read_bytes()).hexdigest()
                for rel in md_files
            }

            sync_script = str(ROOT_DIR / "harness" / "sync_evidence_docs.py")
            res = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath)],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res.returncode, 0, "Expected non-zero exit on pattern error")
            self.assertIn("ERROR: required pattern not found", res.stderr)

            digests_after = {
                rel: hashlib.sha256((tpath / rel).read_bytes()).hexdigest()
                for rel in md_files
            }
            self.assertEqual(
                digests_before,
                digests_after,
                "Files were partially written despite a late pattern error",
            )


if __name__ == "__main__":
    unittest.main()
