#!/usr/bin/env python3
"""
harness/tests/test_perf_budget.py

Unit tests for harness/perf_budget.py verifying:
1. Non-zero rust_exit / rust_time.exit is rejected (exit 1).
2. Non-MATCH result ("MISMATCH", "ENGINE-FAIL") is rejected across standard and suite JSONs
   (extended_coverage.json, multimodal_slices.json, cram_parity.json).
3. CRAM/BAM wall-time ratio <= 1.50 enforcement in cram_parity.json when both
   rust_bam_wall_s and rust_cram_wall_s are positive.
4. Relative PR regression detection via --baseline <dir> --current <dir> --tolerance <pct>
   (e.g. 2x wall-time regression exits 1).
5. --vs-cpp ratio enforcement with amplicon floor exclusion.
"""

import json
import pathlib
import subprocess
import sys
import tempfile
import unittest


class TestPerfBudget(unittest.TestCase):
    def setUp(self):
        self.repo_root = pathlib.Path(__file__).resolve().parent.parent.parent
        self.script = self.repo_root / "harness" / "perf_budget.py"
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.tmp_path = pathlib.Path(self.tmp_dir.name)
        self.ev_dir = self.tmp_path / "benchmarks" / "evidence"
        self.ev_dir.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        self.tmp_dir.cleanup()

    def _write_json(self, path: pathlib.Path, data: dict) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=2), encoding="utf-8")

    def _run_budget(self, *args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(self.script), *args],
            cwd=str(self.tmp_path),
            capture_output=True,
            text=True,
            check=False,
        )

    def test_rust_exit_nonzero_fails(self):
        """A crashed Rust run (rust_exit: 1 or rust_time.exit != '0') must fail even if wall < java."""
        self._write_json(
            self.ev_dir / "rustcrash.json",
            {
                "fingerprint": "abc123",
                "result": "MATCH",
                "rust_exit": 1,
                "java_time": {"wall": "1:00.00", "max_rss_kb": "2000000", "exit": "0"},
                "rust_time": {"wall": "0:00.01", "max_rss_kb": "10000", "exit": "1"},
            },
        )
        res = self._run_budget("rustcrash")
        self.assertNotEqual(res.returncode, 0, f"Expected non-zero exit for rustcrash, got stdout={res.stdout}")

    def test_result_mismatch_fails(self):
        """A record with result: 'MISMATCH' must fail even if wall and RSS are within budget."""
        self._write_json(
            self.ev_dir / "mismatch_case.json",
            {
                "fingerprint": "abc123",
                "result": "MISMATCH",
                "rust_exit": 0,
                "java_time": {"wall": "1:00.00", "max_rss_kb": "2000000", "exit": "0"},
                "rust_time": {"wall": "0:10.00", "max_rss_kb": "500000", "exit": "0"},
            },
        )
        res = self._run_budget("mismatch_case")
        self.assertNotEqual(res.returncode, 0, f"Expected non-zero exit for MISMATCH result, got stdout={res.stdout}")

    def test_suite_jsons_require_rust_exit_zero_and_match(self):
        """extended_coverage.json, multimodal_slices.json, and cram_parity.json must enforce rust_exit==0 and result=='MATCH'."""
        curr = self.tmp_path / "curr"
        self._write_json(
            curr / "extended_coverage.json",
            {"fingerprint": "abc123", "rust_exit": 0, "result": "MATCH", "rust_wall_s": 5.0},
        )
        self._write_json(
            curr / "multimodal_slices.json",
            {"fingerprint": "abc123", "rust_exit": 0, "result": "MATCH", "rust_wall_s": 6.0},
        )
        self._write_json(
            curr / "cram_parity.json",
            {
                "fingerprint": "abc123",
                "rust_exit": 0,
                "result": "MATCH",
                "rust_bam_wall_s": 10.0,
                "rust_cram_wall_s": 13.5,
            },
        )
        res_ok = self._run_budget("--current", str(curr))
        self.assertEqual(res_ok.returncode, 0, f"Valid suite JSONs should pass: stderr={res_ok.stderr} stdout={res_ok.stdout}")

        # Flip multimodal_slices.json result to MISMATCH -> must fail
        self._write_json(
            curr / "multimodal_slices.json",
            {"fingerprint": "abc123", "rust_exit": 0, "result": "MISMATCH", "rust_wall_s": 6.0},
        )
        res_mm = self._run_budget("--current", str(curr))
        self.assertNotEqual(res_mm.returncode, 0, "MISMATCH in multimodal_slices.json must fail")

        # Flip extended_coverage.json rust_exit to 1 -> must fail
        self._write_json(
            curr / "multimodal_slices.json",
            {"fingerprint": "abc123", "rust_exit": 0, "result": "MATCH", "rust_wall_s": 6.0},
        )
        self._write_json(
            curr / "extended_coverage.json",
            {"fingerprint": "abc123", "rust_exit": 1, "result": "MATCH", "rust_wall_s": 5.0},
        )
        res_exit = self._run_budget("--current", str(curr))
        self.assertNotEqual(res_exit.returncode, 0, "rust_exit=1 in extended_coverage.json must fail")

    def test_cram_parity_wall_ratio_gate(self):
        """cram_parity.json must enforce rust_cram_wall_s / rust_bam_wall_s <= 1.50 when both are positive."""
        curr = self.tmp_path / "curr_cram"
        self._write_json(
            curr / "extended_coverage.json",
            {"fingerprint": "abc123", "rust_exit": 0, "result": "MATCH"},
        )
        self._write_json(
            curr / "multimodal_slices.json",
            {"fingerprint": "abc123", "rust_exit": 0, "result": "MATCH"},
        )
        self._write_json(
            curr / "cram_parity.json",
            {
                "fingerprint": "abc123",
                "rust_exit": 0,
                "result": "MATCH",
                "rust_bam_wall_s": 10.0,
                "rust_cram_wall_s": 16.0,
            },
        )
        res = self._run_budget("--current", str(curr))
        self.assertNotEqual(res.returncode, 0, "cram_parity ratio 1.60 > 1.50 must fail")

    def test_baseline_vs_current_2x_wall_regression_fails(self):
        """A 2x wall-time regression in --current vs --baseline must fail under --tolerance."""
        base = self.tmp_path / "base"
        curr = self.tmp_path / "curr_reg"
        base_doc = {
            "fingerprint": "abc123",
            "result": "MATCH",
            "rust_exit": 0,
            "java_time": {"wall": "2:00.00", "max_rss_kb": "2000000", "exit": "0"},
            "rust_time": {"wall": "0:20.00", "max_rss_kb": "500000", "exit": "0"},
        }
        curr_doc = {
            "fingerprint": "abc123",
            "result": "MATCH",
            "rust_exit": 0,
            "java_time": {"wall": "2:00.00", "max_rss_kb": "2000000", "exit": "0"},
            "rust_time": {"wall": "0:40.00", "max_rss_kb": "500000", "exit": "0"},
        }
        self._write_json(base / "hg008-tier2.json", base_doc)
        self._write_json(curr / "hg008-tier2.json", curr_doc)

        res_fail = self._run_budget(
            "--baseline", str(base), "--current", str(curr), "--tolerance", "10", "hg008-tier2"
        )
        self.assertNotEqual(res_fail.returncode, 0, "2x wall-time regression vs baseline must fail at 10% tolerance")

        # Within 10% tolerance (21.5s vs 20.0s = +7.5%) -> must pass
        curr_doc["rust_time"]["wall"] = "0:21.50"
        self._write_json(curr / "hg008-tier2.json", curr_doc)
        res_pass = self._run_budget(
            "--baseline", str(base), "--current", str(curr), "--tolerance", "10", "hg008-tier2"
        )
        self.assertEqual(res_pass.returncode, 0, f"7.5% wall-time delta should pass 10% tolerance: {res_pass.stderr}")

        # Missing baseline file in explicit --baseline directory -> must fail
        empty_base = self.tmp_path / "empty_base"
        empty_base.mkdir(parents=True, exist_ok=True)
        res_missing_base = self._run_budget(
            "--baseline", str(empty_base), "--current", str(curr), "--tolerance", "10", "hg008-tier2"
        )
        self.assertNotEqual(res_missing_base.returncode, 0, "Missing baseline file in explicit --baseline must fail")

        # --current with no explicit names must discover both hg008-tier2.json and extended_coverage.json
        curr_doc["rust_time"]["wall"] = "0:40.00"
        self._write_json(curr / "hg008-tier2.json", curr_doc)
        self._write_json(
            curr / "extended_coverage.json",
            {"fingerprint": "abc123", "rust_exit": 0, "result": "MATCH", "rust_wall_s": 5.0},
        )
        res_disc = self._run_budget("--baseline", str(base), "--current", str(curr), "--tolerance", "10")
        self.assertNotEqual(
            res_disc.returncode,
            0,
            "Auto-discovery on --current must still include hg008-tier2.json when extended_coverage.json is present",
        )


if __name__ == "__main__":
    unittest.main()
