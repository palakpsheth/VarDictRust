#!/usr/bin/env python3
"""Self-test for harness/generate_linkedin_carousel.py (Class-0 integrity & v1.0.6 60-cell GCP study verification)."""

from __future__ import annotations

import importlib.util
import io
import json
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "harness" / "generate_linkedin_carousel.py"
EVIDENCE_DIR = ROOT / "benchmarks" / "evidence"
GCP_METRICS = ROOT / "benchmarks" / "gcp-study-v106" / "master_gcp_metrics.json"


def read_workspace_version(cargo_toml_path: Path = ROOT / "Cargo.toml") -> str:
    cargo_toml = cargo_toml_path.read_text(encoding="utf-8")
    ver_m = re.search(
        r'^\[workspace\.package\]\s*\n(?:(?!\[).*\n)*?version\s*=\s*"([^"]+)"',
        cargo_toml,
        re.MULTILINE,
    )
    assert ver_m is not None, f"Could not parse [workspace.package].version from {cargo_toml_path}"
    return f"v{ver_m.group(1)}"


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
        expected_ver = read_workspace_version()

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

            # 4. Fail-closed test: mismatched fingerprint in --gcp-metrics must exit non-zero
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

            if shutil.which("gs") is not None:
                # 5. Happy-path test on certified benchmarks/evidence/*.json + --gcp-metrics
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

                # Verify current Cargo.toml workspace version, evidence fingerprint, v1.0.6 GCP study fingerprint 5865a2a0b24a8c17, and 60-cell GCP study highlights
                expected_fp = json.loads((EVIDENCE_DIR / "hg002-tier2.json").read_text(encoding="utf-8"))[
                    "fingerprint"
                ]
                for expected_token in (
                    expected_ver,
                    "v1.0.6",
                    expected_fp,
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
                fb_svgs = "\n".join(
                    (out_fallback / f"slide_{i:02d}.svg").read_text(encoding="utf-8") for i in range(1, 11)
                )
                self.assertIn(expected_ver, fb_post)
                self.assertIn(expected_ver, fb_svgs)
                self.assertIsNone(re.search(r"-\d+(?:\.\d+)?%", fb_svgs + "\n" + fb_post))
            else:
                print(
                    "SKIP: Ghostscript ('gs') not found on PATH; skipping PDF/PNG happy-path rendering checks",
                    file=sys.stderr,
                )

    def test_dynamic_cargo_version_and_gs_less_fallback(self) -> None:
        # 1. Verify generate_linkedin_carousel.py reads [workspace.package].version dynamically from ROOT / "Cargo.toml"
        spec = importlib.util.spec_from_file_location("generate_linkedin_carousel", SCRIPT)
        self.assertIsNotNone(spec)
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        self.assertTrue(
            hasattr(mod, "ROOT"),
            "generate_linkedin_carousel.py must define ROOT = Path(__file__).resolve().parents[1]",
        )

        docs = mod.load_and_validate_evidence(EVIDENCE_DIR)
        gcp_doc = mod.load_gcp_metrics(GCP_METRICS, docs["hg002-tier2.json"]["fingerprint"])
        with tempfile.TemporaryDirectory() as tmp:
            fake_root = Path(tmp)
            (fake_root / "Cargo.toml").write_text(
                '[workspace]\nmembers = ["crates/vardictrust-core"]\n\n[workspace.package]\nversion = "1.0.8"\n',
                encoding="utf-8",
            )
            with mock.patch.object(mod, "ROOT", fake_root):
                slides, post_md = mod.build_slides(docs, "rustc 1.96.0", gcp_doc=gcp_doc)
            rendered_svgs = "\n".join(s.to_svg() for s in slides)
            self.assertIn("v1.0.8", post_md)
            self.assertIn("v1.0.8", rendered_svgs)
            self.assertNotIn("v1.0.7", post_md)
            self.assertNotIn("v1.0.7", rendered_svgs)

        # 2. Verify Ghostscript-less fallback executes fail-closed JSON checks and emits clear skip note
        orig_which = shutil.which
        orig_run = subprocess.run
        run_calls: list[list[str]] = []

        def tracked_run(cmd, *args, **kwargs):
            run_calls.append([str(c) for c in cmd])
            return orig_run(cmd, *args, **kwargs)

        stderr_buf = io.StringIO()
        with mock.patch(
            "shutil.which", side_effect=lambda name: None if name == "gs" else orig_which(name)
        ), mock.patch("subprocess.run", side_effect=tracked_run), mock.patch("sys.stderr", stderr_buf):
            self.test_carousel_integrity_and_generation()

        self.assertEqual(len(run_calls), 4, "All fail-closed JSON integrity checks must run when gs is absent")
        self.assertIn("SKIP: Ghostscript ('gs') not found on PATH", stderr_buf.getvalue())


if __name__ == "__main__":
    unittest.main()
