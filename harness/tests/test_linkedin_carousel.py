#!/usr/bin/env python3
"""Self-test for harness/generate_linkedin_carousel.py (Class-0 integrity & v1.0.6 60-cell GCP study verification)."""

import json
import re
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
GCP_METRICS = ROOT / "benchmarks" / "gcp-study-v106" / "master_gcp_metrics.json"


def read_png_dims(path: Path) -> tuple[int, int]:
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", f"{path} is not a valid PNG"
    assert data[12:16] == b"IHDR", f"{path} missing IHDR chunk"
    w, h = struct.unpack(">II", data[16:24])
    return w, h


class TestLinkedinCarousel(unittest.TestCase):
    def test_carousel_integrity_and_generation(self) -> None:
        self.assertTrue(SCRIPT.is_file(), f"Missing {SCRIPT}")
        self.assertTrue(GCP_METRICS.is_file(), f"Missing {GCP_METRICS}")

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

            # 4. Happy-path test on certified benchmarks/evidence/*.json + --gcp-metrics
            out_dir = tmp_path / "out_ok"
            res_ok = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--evidence-dir",
                    str(EVIDENCE_DIR),
                    "--gcp-metrics",
                    str(GCP_METRICS),
                    "--out-dir",
                    str(out_dir),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_ok.returncode, 0, msg=f"stderr={res_ok.stderr}")

            svg_texts: list[str] = []
            for idx in range(1, 11):
                png = out_dir / f"slide_{idx:02d}.png"
                svg = out_dir / f"slide_{idx:02d}.svg"
                self.assertTrue(png.is_file() and png.stat().st_size > 5000)
                self.assertTrue(svg.is_file() and svg.stat().st_size > 1000)
                self.assertEqual(read_png_dims(png), (1080, 1350))
                svg_texts.append(svg.read_text(encoding="utf-8"))

            pdf = out_dir / "vardictrust_v1.0.6_carousel.pdf"
            self.assertTrue(pdf.is_file() and pdf.stat().st_size > 10000)
            self.assertFalse((out_dir / "vardictrust_v1.0.3_carousel.pdf").exists())

            post = out_dir / "post.md"
            self.assertTrue(post.is_file())
            post_text = post.read_text(encoding="utf-8")
            all_rendered = "\n".join(svg_texts) + "\n" + post_text

            # Verify v1.0.6, fingerprint 5865a2a0b24a8c17, and 60-cell GCP study highlights
            for expected_token in (
                "v1.0.6",
                "5865a2a0b24a8c17",
                "623,444",
                "227,562",
                "8,111,171",
                "9,643,769",
                "13,067,287",
                "v5.0q",
                "CMRG",
                "84/85",
                "98.82%",
                "0.5%",
                "18/18",
            ):
                self.assertIn(expected_token, all_rendered)
                self.assertIn(expected_token, post_text)

            # Verify Slide 6+ specifically highlights Stage 5 & Stage 6 findings
            slides_6_to_10 = "\n".join(svg_texts[5:])
            for s6_token in ("v5.0q", "CMRG", "84/85", "0.5%", "18/18"):
                self.assertIn(s6_token, slides_6_to_10)

            # Verify zero stale versions/fingerprints or hardcoded local-Ryzen strings
            for forbidden in (
                "v1.0.3",
                "27f3667b32ac14f5",
                "5a4e93c556a2c709",
                "648.98s",
                "186.43s",
                "7:11.99",
                "7:46.34",
            ):
                self.assertNotIn(forbidden, all_rendered)

            # Verify zero negative "-%" strings in any slide SVG or post.md
            neg_pct = re.search(r"-\d+(?:\.\d+)?%", all_rendered)
            self.assertIsNone(
                neg_pct,
                f"Found negative percentage string in carousel output: {neg_pct.group(0) if neg_pct else ''!r}",
            )

            # 5. Fail-closed test: mismatched fingerprint in --gcp-metrics must exit non-zero
            bad_gcp = tmp_path / "bad_gcp_metrics.json"
            gcp_data = json.loads(GCP_METRICS.read_text(encoding="utf-8"))
            gcp_data["meta"]["fingerprint"] = "deadbeefdeadbeef"
            bad_gcp.write_text(json.dumps(gcp_data), encoding="utf-8")
            res_bad_gcp = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--evidence-dir",
                    str(EVIDENCE_DIR),
                    "--gcp-metrics",
                    str(bad_gcp),
                    "--out-dir",
                    str(tmp_path / "out_bad_gcp"),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res_bad_gcp.returncode, 0)
            self.assertIn("fingerprint", res_bad_gcp.stderr.lower())

            # 6. Evidence-only fallback test (when --gcp-metrics file is absent)
            out_fallback = tmp_path / "out_fallback"
            res_fb = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--evidence-dir",
                    str(EVIDENCE_DIR),
                    "--gcp-metrics",
                    str(tmp_path / "nonexistent_gcp_metrics.json"),
                    "--out-dir",
                    str(out_fallback),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_fb.returncode, 0, msg=f"stderr={res_fb.stderr}")
            fb_post = (out_fallback / "post.md").read_text(encoding="utf-8")
            fb_svgs = "\n".join((out_fallback / f"slide_{i:02d}.svg").read_text(encoding="utf-8") for i in range(1, 11))
            self.assertIsNone(re.search(r"-\d+(?:\.\d+)?%", fb_svgs + "\n" + fb_post))


if __name__ == "__main__":
    unittest.main()
