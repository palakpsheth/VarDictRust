#!/usr/bin/env python3
"""Self-test for harness/generate_linkedin_carousel.py (Class-0 integrity & output verification)."""

import json
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "harness" / "generate_linkedin_carousel.py"
EVIDENCE_DIR = ROOT / "benchmarks" / "evidence"


def read_png_dims(path: Path) -> tuple[int, int]:
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", f"{path} is not a valid PNG"
    assert data[12:16] == b"IHDR", f"{path} missing IHDR chunk"
    w, h = struct.unpack(">II", data[16:24])
    return w, h


class TestLinkedinCarousel(unittest.TestCase):
    def test_carousel_integrity_and_generation(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            # 1. Fail-closed test: mismatched fingerprint must exit non-zero
            bad_ev = tmp_path / "bad_fp_evidence"
            shutil.copytree(EVIDENCE_DIR, bad_ev)
            hg002_path = bad_ev / "hg002-tier2.json"
            doc = json.loads(hg002_path.read_text())
            doc["fingerprint"] = "0000000000000000"
            hg002_path.write_text(json.dumps(doc))
            res_bad_fp = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--evidence-dir",
                    str(bad_ev),
                    "--out-dir",
                    str(tmp_path / "out_bad_fp"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res_bad_fp.returncode, 0)
            self.assertIn("fingerprint", res_bad_fp.stderr.lower())

            # 2. Fail-closed test: non-passing result ("MISMATCH") must exit non-zero
            bad_res = tmp_path / "bad_res_evidence"
            shutil.copytree(EVIDENCE_DIR, bad_res)
            seqc2_path = bad_res / "seqc2-exome-tier2.json"
            doc2 = json.loads(seqc2_path.read_text())
            doc2["result"] = "MISMATCH"
            seqc2_path.write_text(json.dumps(doc2))
            res_bad_status = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--evidence-dir",
                    str(bad_res),
                    "--out-dir",
                    str(tmp_path / "out_bad_res"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res_bad_status.returncode, 0)

            # 3. Fail-closed test: missing required evidence file must exit non-zero
            missing_ev = tmp_path / "missing_evidence"
            shutil.copytree(EVIDENCE_DIR, missing_ev)
            (missing_ev / "scale.json").unlink()
            res_missing = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--evidence-dir",
                    str(missing_ev),
                    "--out-dir",
                    str(tmp_path / "out_missing"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res_missing.returncode, 0)

            # 4. Happy-path test on certified benchmarks/evidence/*.json
            out_dir = tmp_path / "out_ok"
            res_ok = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--evidence-dir",
                    str(EVIDENCE_DIR),
                    "--out-dir",
                    str(out_dir),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_ok.returncode, 0, msg=f"stderr={res_ok.stderr}")

            for idx in range(1, 11):
                png = out_dir / f"slide_{idx:02d}.png"
                svg = out_dir / f"slide_{idx:02d}.svg"
                self.assertTrue(png.is_file() and png.stat().st_size > 5000)
                self.assertTrue(svg.is_file() and svg.stat().st_size > 1000)
                self.assertEqual(read_png_dims(png), (1080, 1350))

            pdf = out_dir / "vardictrust_v1.0.3_carousel.pdf"
            self.assertTrue(pdf.is_file() and pdf.stat().st_size > 10000)

            post = out_dir / "post.md"
            self.assertTrue(post.is_file())
            post_text = post.read_text()
            for expected_token in ("27f3667b32ac14f5", "14.75s", "479 MB", "1.49 GB", "623,444", "227,562"):
                self.assertIn(expected_token, post_text)


if __name__ == "__main__":
    unittest.main()
