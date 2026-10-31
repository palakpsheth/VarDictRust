#!/usr/bin/env python3
"""Unit tests for paper/aggregate_gcp_study.py (60-directory completeness & Class-0 schema)."""

import csv
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "paper" / "aggregate_gcp_study.py"
OUT_DIR = ROOT / "benchmarks" / "gcp-study-v106"


class TestAggregateGcpStudy(unittest.TestCase):
    def test_aggregate_and_verify_all_60_cells(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        res = subprocess.run(
            [sys.executable, str(SCRIPT), "--verify"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=f"stderr={res.stderr}\nstdout={res.stdout}")

        master_path = OUT_DIR / "master_gcp_metrics.json"
        inv_path = OUT_DIR / "cells_60_inventory.tsv"
        sha_path = OUT_DIR / "SHA256SUMS.tsv"
        self.assertTrue(master_path.is_file())
        self.assertTrue(inv_path.is_file())
        self.assertTrue(sha_path.is_file())

        doc = json.loads(master_path.read_text(encoding="utf-8"))
        self.assertEqual(doc["meta"]["version"], "1.0.6")
        self.assertEqual(doc["meta"]["fingerprint"], "5865a2a0b24a8c17")
        self.assertEqual(doc["meta"]["perf_fingerprint"], "c490c36317b8f25b")
        self.assertEqual(len(doc["cells"]), 60)

        with inv_path.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 60)
        for r in rows:
            self.assertEqual(r["exit_code"], "0", f"Cell {r['cell_id']} non-zero exit")
            self.assertTrue(r["paper_targets"].strip(), f"Orphan cell {r['cell_id']}")


if __name__ == "__main__":
    unittest.main()
