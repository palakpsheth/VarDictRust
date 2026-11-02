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

            # Mutate hg008-tier2.json and hg008-cram-twin.json timings/RSS in temp copy
            hg008_path = ev_dst / "hg008-tier2.json"
            hg008_data = json.loads(hg008_path.read_text(encoding="utf-8"))
            hg008_data["rust_time"]["wall"] = "0:38.20"
            hg008_data["rust_time"]["max_rss_kb"] = "531860"
            hg008_path.write_text(
                json.dumps(hg008_data, indent=2, sort_keys=True) + "\n",
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
            scale_data["extra"]["speedup_8"] = 3.25
            scale_data["extra"]["rss_ratio_chr1_vs_10mb"] = 1.15
            scale_path.write_text(
                json.dumps(scale_data, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

            # --check must detect drift across all 4 surfaces
            res_drift = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath), "--check"],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res_drift.returncode, 0, "--check should fail when timings drift")
            for expected_surface in (
                "README.md",
                "benchmarks/README.md",
                "docs/validation/omnibus-empirical-validation.md",
                "docs/testing/test-matrix-and-ci-audit.md",
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

            omnibus_text = (tpath / "docs" / "validation" / "omnibus-empirical-validation.md").read_text(encoding="utf-8")
            self.assertIn("**3.25× parallel**", omnibus_text)
            self.assertIn("**1.15× flat RSS scaling**", omnibus_text)
            self.assertIn("**$1.15\\times$ (`1.106\\times`)**", omnibus_text)
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

            # --check must now exit 0
            res_after = subprocess.run(
                [sys.executable, sync_script, "--repo-root", str(tpath), "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_after.returncode, 0, f"--check after sync failed: {res_after.stderr}")


if __name__ == "__main__":
    unittest.main()
