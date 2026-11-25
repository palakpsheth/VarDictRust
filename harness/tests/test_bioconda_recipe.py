#!/usr/bin/env python3
"""
harness/tests/test_bioconda_recipe.py

Automated verification of the canonical Bioconda recipe for VarDictRust.
Conforms strictly to the Python Standard Library purity constraint (zero PyPI packages).
"""

import importlib.util
import json
import os
import pathlib
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent
if str(ROOT_DIR / "harness") not in sys.path:
    sys.path.insert(0, str(ROOT_DIR / "harness"))

import toml_lite  # noqa: E402
RECIPE_DIR = ROOT_DIR / "packaging" / "bioconda" / "recipes" / "vardictrust"
META_PATH = RECIPE_DIR / "meta.yaml"
BUILD_PATH = RECIPE_DIR / "build.sh"
RUN_TEST_PATH = RECIPE_DIR / "run_test.sh"
NFCORE_DIR = ROOT_DIR / "packaging" / "nf-core"
NFCORE_MOD_DIR = NFCORE_DIR / "modules" / "vardictrust"
SKILLS_DIR = ROOT_DIR / ".agents" / "skills"
EXPECTED_SKILLS = (
    "vardictrust-issue",
    "vardictrust-pr",
    "vardictrust-profiling",
    "vardictrust-release",
)
ALLOWED_HTML_TAGS = frozenset([
    "a", "abbr", "b", "blockquote", "br", "code", "dd", "del", "details",
    "div", "dl", "dt", "em", "h1", "h2", "h3", "h4", "h5", "h6", "hr",
    "i", "img", "ins", "kbd", "li", "ol", "p", "pre", "q", "rp", "rt",
    "ruby", "s", "samp", "small", "span", "strong", "sub", "summary",
    "sup", "table", "tbody", "td", "tfoot", "th", "thead", "tr", "ul",
    "var", "wbr"
])


class TestBiocondaRecipe(unittest.TestCase):
    def setUp(self):
        self.assertTrue(RECIPE_DIR.is_dir(), f"Missing recipe directory: {RECIPE_DIR}")
        self.assertTrue(META_PATH.is_file(), f"Missing meta.yaml: {META_PATH}")
        self.assertTrue(BUILD_PATH.is_file(), f"Missing build.sh: {BUILD_PATH}")
        self.assertFalse(RUN_TEST_PATH.is_file(), f"run_test.sh should not exist in recipe: {RUN_TEST_PATH}")

    def test_recipe_files_executable(self):
        self.assertTrue(os.access(BUILD_PATH, os.X_OK), "build.sh must be executable")

    def test_build_script_content(self):
        content = BUILD_PATH.read_text(encoding="utf-8")
        self.assertIn("cargo-bundle-licenses", content, "build.sh must bundle licenses")
        self.assertIn("cargo install", content, "build.sh must run cargo install")
        self.assertIn("--locked", content, "build.sh must specify --locked")
        self.assertIn("--no-track", content, "build.sh must specify --no-track")
        self.assertIn("--root \"${PREFIX}\"", content, "build.sh must install into ${PREFIX}")
        self.assertIn("--path crates/vardictrust", content, "build.sh must target crates/vardictrust")
        self.assertIn("CARGO_HOME=", content, "build.sh must explicitly set CARGO_HOME for sandboxing")
        self.assertIn("rm -rf \"$CARGO_HOME\"", content, "build.sh must clean up CARGO_HOME cache")
        self.assertIn("trap 'rm -rf \"$CARGO_HOME\"' EXIT", content, "build.sh must set EXIT cleanup trap")
        self.assertIn(
            "CONDA_BUILD_CROSS_COMPILATION",
            content,
            "build.sh must guard native completion/manpage generation behind CONDA_BUILD_CROSS_COMPILATION",
        )
        self.assertIn("vardictrust.1", content, "build.sh must install manpage vardictrust.1")


    def test_cargo_manifest_license(self):
        cargo_cfg = toml_lite.load_path(ROOT_DIR / "Cargo.toml")
        ws_pkg = (cargo_cfg.get("workspace") or {}).get("package") or {}
        self.assertEqual(
            ws_pkg.get("license"),
            "MIT",
            "Root Cargo.toml must declare pure MIT SPDX license expression",
        )

    def test_cargo_workspace_msrv_inheritance(self):
        cargo_cfg = toml_lite.load_path(ROOT_DIR / "Cargo.toml")
        ws_pkg = (cargo_cfg.get("workspace") or {}).get("package") or {}
        self.assertEqual(
            ws_pkg.get("rust-version"),
            "1.91.0",
            'Workspace Cargo.toml must pin rust-version = "1.91.0" under [workspace.package]',
        )
        for member in ("vardictrust", "vardictrust-core"):
            member_path = ROOT_DIR / "crates" / member / "Cargo.toml"
            self.assertTrue(member_path.is_file(), f"Missing crate manifest: {member_path}")
            member_cfg = toml_lite.load_path(member_path)
            rust_ver = (member_cfg.get("package") or {}).get("rust-version")
            self.assertEqual(
                rust_ver,
                {"workspace": True},
                f"crates/{member}/Cargo.toml must inherit rust-version.workspace = true",
            )

    def test_fuzz_manifest_noodles_version_matches_workspace(self):
        root_cfg = toml_lite.load_path(ROOT_DIR / "Cargo.toml")
        fuzz_path = ROOT_DIR / "crates" / "vardictrust-core" / "fuzz" / "Cargo.toml"
        self.assertTrue(fuzz_path.is_file(), f"Missing fuzz manifest: {fuzz_path}")
        fuzz_cfg = toml_lite.load_path(fuzz_path)
        root_noodles = ((root_cfg.get("workspace") or {}).get("dependencies") or {}).get("noodles") or {}
        fuzz_noodles = (fuzz_cfg.get("dependencies") or {}).get("noodles") or {}
        self.assertTrue(root_noodles.get("version"), "Root Cargo.toml missing noodles dependency entry")
        self.assertTrue(
            fuzz_noodles.get("version"),
            "crates/vardictrust-core/fuzz/Cargo.toml missing noodles dependency entry",
        )
        self.assertEqual(
            root_noodles["version"],
            fuzz_noodles["version"],
            "crates/vardictrust-core/fuzz/Cargo.toml noodles version must match root Cargo.toml",
        )

        expected_fuzz_bins = {
            "fuzz_modify_cigar",
            "fuzz_realign",
            "fuzz_bed_parser",
            "fuzz_sam_aux",
            "fuzz_parse_and_realign",
            "fuzz_vcf_streamer",
            "fuzz_cram_and_fai",
            "fuzz_somatic_fisher_and_amplicon",
        }
        bins = fuzz_cfg.get("bin") or []
        self.assertEqual(len(bins), 8, f"Expected 8 [[bin]] targets in {fuzz_path}, got {len(bins)}")
        actual_fuzz_bins = {entry.get("name") for entry in bins if isinstance(entry, dict)}
        self.assertEqual(
            actual_fuzz_bins,
            expected_fuzz_bins,
            f"Fuzz Cargo.toml [[bin]] targets mismatch: got {sorted(a for a in actual_fuzz_bins if a)}",
        )
        fuzz_dir = fuzz_path.parent
        for entry in bins:
            rel_target = entry.get("path")
            self.assertTrue(rel_target, f"Missing path in fuzz [[bin]] entry: {entry}")
            target_file = fuzz_dir / rel_target
            self.assertTrue(target_file.is_file(), f"Missing fuzz target source on disk: {target_file}")

        runner_path = ROOT_DIR / "harness" / "run_deep_fuzz.sh"
        self.assertTrue(runner_path.is_file(), f"Missing deep fuzz runner: {runner_path}")
        self.assertTrue(os.access(runner_path, os.X_OK), "harness/run_deep_fuzz.sh must be executable")
        before_status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=str(ROOT_DIR),
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        with tempfile.TemporaryDirectory() as tmp_dir:
            probe_out_dir = pathlib.Path(tmp_dir) / "dry_run_out"
            env = dict(os.environ)
            env["FUZZ_OUT_DIR"] = str(probe_out_dir)
            dry = subprocess.run(
                ["bash", str(runner_path), "--dry-run"],
                cwd=str(ROOT_DIR),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                dry.returncode,
                0,
                f"run_deep_fuzz.sh --dry-run failed ({dry.returncode}):\n{dry.stderr}\n{dry.stdout}",
            )
            self.assertFalse(
                probe_out_dir.exists(),
                f"run_deep_fuzz.sh --dry-run must not create FUZZ_OUT_DIR: {probe_out_dir}",
            )
            for target_name in sorted(expected_fuzz_bins):
                self.assertIn(
                    target_name,
                    dry.stdout,
                    f"run_deep_fuzz.sh --dry-run schedule missing {target_name}",
                )

            # Verify artifact staging, stale-artifact cleanup across re-runs,
            # summary.json 1-to-1 alignment with $ART_OUT_DIR, and SIGTERM -> 130 exit.
            mock_bin = pathlib.Path(tmp_dir) / "mock_bin"
            mock_bin.mkdir(parents=True, exist_ok=True)
            mock_cargo = mock_bin / "cargo"
            mock_cargo.write_text(
                f"#!{sys.executable}\n"
                "import os, pathlib, sys, time\n"
                "args = sys.argv[1:]\n"
                "if args[:3] == ['+nightly', 'fuzz', '--version']:\n"
                "    print('cargo-fuzz 0.12.0')\n"
                "    sys.exit(0)\n"
                "if args[:3] == ['+nightly', 'fuzz', 'build']:\n"
                "    bmarker = os.environ.get('MOCK_BUILD_MARKER', '')\n"
                "    if bmarker:\n"
                "        pathlib.Path(bmarker).write_text('building', encoding='utf-8')\n"
                "        time.sleep(10)\n"
                "    sys.exit(0)\n"
                "if len(args) >= 4 and args[:3] == ['+nightly', 'fuzz', 'run']:\n"
                "    tgt = args[3]\n"
                "    shm_dir = os.environ.get('VDR_FUZZ_SHM_DIR', '')\n"
                "    if not shm_dir or not pathlib.Path(shm_dir).is_dir():\n"
                "        sys.stderr.write('missing VDR_FUZZ_SHM_DIR\\n')\n"
                "        sys.exit(99)\n"
                "    (pathlib.Path(shm_dir) / f'vdr-fuzz-{tgt}-leak').mkdir(exist_ok=True)\n"
                "    shm_rec = os.environ.get('MOCK_SHM_RECORD', '')\n"
                "    if shm_rec:\n"
                "        pathlib.Path(shm_rec).write_text(shm_dir, encoding='utf-8')\n"
                "    marker = os.environ.get('MOCK_WORKER_MARKER', '')\n"
                "    if marker:\n"
                "        pathlib.Path(marker).write_text(tgt, encoding='utf-8')\n"
                "        time.sleep(10)\n"
                "        sys.exit(0)\n"
                "    art_prefix = ''\n"
                "    for a in args:\n"
                "        if a.startswith('-artifact_prefix='):\n"
                "            art_prefix = a.split('=', 1)[1]\n"
                "    if os.environ.get('MOCK_CRASH_TARGET') == tgt and art_prefix:\n"
                "        p = pathlib.Path(art_prefix) / 'crash-mock123'\n"
                "        p.parent.mkdir(parents=True, exist_ok=True)\n"
                "        p.write_bytes(b'boom')\n"
                "        sys.exit(1)\n"
                "    print('#12\\tINITED cov: 10 ft: 15 corp: 1/1b exec/s: 12 rss: 32Mb')\n"
                "    print('Done 12 runs in 1 second(s)')\n"
                "    sys.exit(0)\n"
                "sys.exit(1)\n",
                encoding="utf-8",
            )
            mock_cargo.chmod(0o755)

            live_out_dir = pathlib.Path(tmp_dir) / "live_out"
            mock_env = dict(os.environ)
            mock_env["PATH"] = f"{mock_bin}:{mock_env.get('PATH', '')}"
            mock_env["FUZZ_OUT_DIR"] = str(live_out_dir)
            mock_env["MIN_AVAIL_GB"] = "0"
            mock_env["TMPDIR"] = tmp_dir

            # 1) Simulated crash run populates ART_OUT_DIR and fails summary.json 1-to-1
            crash_env = dict(mock_env, MOCK_CRASH_TARGET="fuzz_realign")
            crash_res = subprocess.run(
                ["bash", str(runner_path), "--fast-secs", "1", "--deep-secs", "1", "--jobs", "2"],
                cwd=str(ROOT_DIR),
                env=crash_env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(crash_res.returncode, 0, "Crash run must exit non-zero")
            staged_crash = live_out_dir / "artifacts" / "fuzz_realign__crash-mock123"
            self.assertTrue(staged_crash.is_file(), "Crash artifact must be staged in ART_OUT_DIR")
            crash_summary = json.loads((live_out_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(crash_summary["status"], "FAIL")
            self.assertEqual(crash_summary["total_artifacts"], 1)
            self.assertEqual(crash_summary["artifacts"], ["fuzz_realign__crash-mock123"])

            # 2) Clean re-run into same FUZZ_OUT_DIR clears stale artifacts, scopes /dev/shm cleanup
            # strictly to campaign VDR_FUZZ_SHM_DIR (preserving external /dev/shm/vdr-fuzz-*), and reports PASS
            shm_record_file = pathlib.Path(tmp_dir) / "recorded_shm_dir.txt"
            shm_base = pathlib.Path("/dev/shm") if os.access("/dev/shm", os.W_OK | os.X_OK) else pathlib.Path(tmp_dir)
            with tempfile.TemporaryDirectory(prefix="vdr-fuzz-external-", dir=str(shm_base)) as ext_shm:
                clean_env = dict(mock_env, MOCK_SHM_RECORD=str(shm_record_file))
                clean_res = subprocess.run(
                    ["bash", str(runner_path), "--fast-secs", "1", "--deep-secs", "1", "--jobs", "2"],
                    cwd=str(ROOT_DIR),
                    env=clean_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(
                    clean_res.returncode,
                    0,
                    f"Clean re-run failed ({clean_res.returncode}):\n{clean_res.stderr}\n{clean_res.stdout}",
                )
                self.assertTrue(
                    pathlib.Path(ext_shm).is_dir(),
                    "External same-user /dev/shm/vdr-fuzz-* directory must not be deleted by run_deep_fuzz.sh",
                )
                recorded_shm = pathlib.Path(shm_record_file.read_text(encoding="utf-8").strip())
                self.assertFalse(
                    recorded_shm.exists(),
                    f"Campaign-scoped VDR_FUZZ_SHM_DIR must be removed on exit: {recorded_shm}",
                )
            self.assertFalse(staged_crash.exists(), "Stale crash artifact must be removed on live re-run")
            clean_summary = json.loads((live_out_dir / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(clean_summary["status"], "PASS")
            self.assertEqual(clean_summary["total_artifacts"], 0)
            self.assertEqual(clean_summary["artifacts"], [])

            # 3) SIGTERM during worker wait exits 130 (not interrupted wait rc 143)
            marker_file = pathlib.Path(tmp_dir) / "worker_started"
            sig_env = dict(mock_env, MOCK_WORKER_MARKER=str(marker_file), FUZZ_OUT_DIR=str(pathlib.Path(tmp_dir) / "sig_out"))
            sig_proc = subprocess.Popen(
                ["bash", str(runner_path), "--fast-secs", "10", "--deep-secs", "10", "--jobs", "1"],
                cwd=str(ROOT_DIR),
                env=sig_env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            try:
                deadline = time.monotonic() + 5.0
                while not marker_file.exists() and time.monotonic() < deadline:
                    time.sleep(0.02)
                self.assertTrue(marker_file.exists(), "Mock fuzz worker did not start within 5s")
                time.sleep(0.05)
                sig_proc.send_signal(signal.SIGTERM)
                sig_rc = sig_proc.wait(timeout=10)
                self.assertEqual(sig_rc, 130, f"Expected exit 130 on SIGTERM, got {sig_rc}")
            finally:
                if sig_proc.poll() is None:
                    sig_proc.kill()
                    sig_proc.wait(timeout=5)
                if sig_proc.stdout:
                    sig_proc.stdout.close()
                if sig_proc.stderr:
                    sig_proc.stderr.close()

            # 4) SIGTERM during sequential pre-build clears stale summary.json/logs, interrupts (<4s), and exits 130
            bsig_out_dir = pathlib.Path(tmp_dir) / "bsig_out"
            bsig_out_dir.mkdir(parents=True, exist_ok=True)
            stale_summary = bsig_out_dir / "summary.json"
            stale_summary.write_text('{"status": "PASS"}', encoding="utf-8")
            bmarker_file = pathlib.Path(tmp_dir) / "build_started"
            bsig_env = dict(mock_env, MOCK_BUILD_MARKER=str(bmarker_file), FUZZ_OUT_DIR=str(bsig_out_dir))
            bsig_proc = subprocess.Popen(
                ["bash", str(runner_path), "--fast-secs", "10", "--deep-secs", "10", "--jobs", "1"],
                cwd=str(ROOT_DIR),
                env=bsig_env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            try:
                deadline = time.monotonic() + 5.0
                while not bmarker_file.exists() and time.monotonic() < deadline:
                    time.sleep(0.02)
                self.assertTrue(bmarker_file.exists(), "Mock fuzz pre-build did not start within 5s")
                time.sleep(0.05)
                t_sig = time.monotonic()
                bsig_proc.send_signal(signal.SIGTERM)
                bsig_rc = bsig_proc.wait(timeout=4)
                self.assertEqual(bsig_rc, 130, f"Expected exit 130 on SIGTERM during pre-build, got {bsig_rc}")
                self.assertLess(time.monotonic() - t_sig, 4.0, "Pre-build ignored SIGTERM instead of interrupting wait")
                self.assertFalse(stale_summary.exists(), "Stale PASS summary.json must be cleared before pre-build")
            finally:
                if bsig_proc.poll() is None:
                    bsig_proc.kill()
                    bsig_proc.wait(timeout=5)
                if bsig_proc.stdout:
                    bsig_proc.stdout.close()
                if bsig_proc.stderr:
                    bsig_proc.stderr.close()
        after_status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=str(ROOT_DIR),
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        self.assertEqual(
            before_status,
            after_status,
            f"run_deep_fuzz.sh tests dirtied git working tree:\n{after_status}",
        )

    def test_cargo_deny_license_allow_list(self):
        deny_path = ROOT_DIR / "deny.toml"
        content = deny_path.read_text(encoding="utf-8")
        self.assertNotIn(
            '"GPL-2.0-or-later"',
            content,
            "deny.toml must forbid GPL-2.0-or-later",
        )
        self.assertIn('"MIT"', content)
        self.assertIn("[graph]", content, "deny.toml must define [graph] section")
        self.assertIn("all-features = true", content, "deny.toml must set all-features = true")
        self.assertIn("[sources]", content, "deny.toml must define [sources] section")
        self.assertIn('unknown-registry = "deny"', content, "deny.toml must deny unknown registries")
        self.assertIn('unknown-git = "deny"', content, "deny.toml must deny unknown git sources")
        self.assertIn(
            'allow-registry = ["https://github.com/rust-lang/crates.io-index"]',
            content,
            "deny.toml must allow only crates.io-index",
        )

        audit_path = ROOT_DIR / ".cargo" / "audit.toml"
        self.assertTrue(audit_path.is_file(), f"Missing .cargo/audit.toml: {audit_path}")
        audit_content = audit_path.read_text(encoding="utf-8")
        self.assertIn('ignore = ["RUSTSEC-2024-0436"]', audit_content)
        self.assertIn('informational_warnings = ["unmaintained", "unsound", "notice"]', audit_content)
        self.assertIn('deny = ["warnings"]', audit_content)
        self.assertIn("quiet = false", audit_content)

        dependabot_path = ROOT_DIR / ".github" / "dependabot.yml"
        self.assertTrue(dependabot_path.is_file(), f"Missing .github/dependabot.yml: {dependabot_path}")
        dependabot_content = dependabot_path.read_text(encoding="utf-8")
        self.assertIn('package-ecosystem: "cargo"', dependabot_content)
        self.assertIn('package-ecosystem: "github-actions"', dependabot_content)
        self.assertIn('package-ecosystem: "pip"', dependabot_content)
        self.assertIn('directory: "/docs"', dependabot_content)
        self.assertIn("docs-dependencies:", dependabot_content)
        self.assertIn('dependency-type: "all"', dependabot_content)
        self.assertIn("transitive-dependencies:", dependabot_content)
        self.assertIn("exclude-patterns:", dependabot_content)
        expected_cargo_patterns = (
            "noodles*",
            "flate2*",
            "fishers_exact*",
            "rustc-hash*",
            "clap*",
            "rayon*",
            "anyhow*",
            "thiserror*",
            "log*",
            "env_logger*",
            "regex*",
            "tempfile*",
            "criterion*",
            "indexmap*",
            "tikv-jemallocator*",
            "pprof*",
            "inferno*",
            "dhat*",
            "jemalloc_pprof*",
            "libc*",
        )
        for pattern in expected_cargo_patterns:
            self.assertIn(
                f'"{pattern}"',
                dependabot_content,
                f"dependabot.yml missing direct Cargo pattern {pattern!r}",
            )

    def test_container_runner_security(self):
        runner_path = ROOT_DIR / "harness" / "test_bioconda_recipe.sh"
        content = runner_path.read_text(encoding="utf-8")
        self.assertIn(
            "BIOCONDA_BUILD_IMAGE=",
            content,
            "Runner must allow overriding the build image",
        )
        self.assertIn(
            "@sha256:",
            content,
            "Runner must pin build image by immutable digest",
        )
        self.assertIn(
            'WORKSPACE_MOUNT_MODE="${WORKSPACE_MOUNT_MODE:-ro}"',
            content,
            "Runner must default workspace mount to read-only",
        )

    def test_meta_yaml_required_sections(self):
        meta = META_PATH.read_text(encoding="utf-8")
        required_sections = [
            "package:",
            "source:",
            "build:",
            "requirements:",
            "test:",
            "about:",
            "extra:",
        ]
        for sec in required_sections:
            self.assertIn(sec, meta, f"meta.yaml missing required section: {sec}")

    def test_meta_yaml_metadata(self):
        meta = META_PATH.read_text(encoding="utf-8")
        self.assertIn('{% set name = "vardictrust" %}', meta, "Package name must be 'vardictrust'")
        ver_m = re.search(r'\{%\s*set\s+version\s*=\s*"([^"]+)"\s*%\}', meta)
        self.assertIsNotNone(ver_m, 'meta.yaml must define {% set version = "..." %}')
        recipe_version = ver_m.group(1)
        if recipe_version == "1.0.10":
            self.assertIn(
                "license: MIT AND GPL-2.0-or-later\n",
                meta,
                "Pinned 1.0.10 recipe must match v1.0.10 tarball license MIT AND GPL-2.0-or-later",
            )
            self.assertIn(
                "license_family: GPL2",
                meta,
                "Pinned 1.0.10 recipe must declare license_family: GPL2",
            )
        else:
            self.assertIn("license: MIT\n", meta, "License must be MIT")
            self.assertIn("license_family: MIT", meta, "License family must be MIT")
            self.assertNotIn("GPL", meta, "meta.yaml must not contain GPL")
        self.assertIn("LICENSE-3RD-PARTY", meta, "license_file must include LICENSE-3RD-PARTY")
        self.assertIn("THIRDPARTY.yml", meta, "license_file must include THIRDPARTY.yml")
        self.assertIn("recipe-maintainers:", meta, "extra must list recipe-maintainers")
        self.assertIn("- palakpsheth", meta, "Maintainer must include palakpsheth")
        self.assertIn("compiler('rust')", meta, "Requirements must require rust compiler")
        self.assertIn("compiler('c')", meta, "Requirements must require c compiler")
        self.assertIn("stdlib('c')", meta, "Requirements must require stdlib('c')")
        self.assertIn("- make", meta, "Requirements must include make for tikv-jemalloc-sys build.rs")
        self.assertIn("run_exports:", meta, "Bioconda linter requires run_exports pinning")
        self.assertIn("pin_subpackage", meta, "run_exports must use pin_subpackage")
        self.assertIn("binary_relocation: false", meta, "meta.yaml must declare binary_relocation: false")
        self.assertIn("osx-arm64", meta, "Must target osx-arm64")
        self.assertIn("linux-aarch64", meta, "Must target linux-aarch64")
        self.assertIn(
            "doc_url: https://vardictrust.readthedocs.io/",
            meta,
            "Bioconda meta.yaml doc_url must point to https://vardictrust.readthedocs.io/",
        )
        summary_m = re.search(r'^\s*summary:\s*"([^"]+)"\s*$', meta, re.MULTILINE)
        self.assertIsNotNone(summary_m, "meta.yaml must define double-quoted about.summary")
        summary = summary_m.group(1)
        self.assertLessEqual(
            len(summary),
            80,
            f"Bioconda about.summary must be <= 80 chars (got {len(summary)}): {summary!r}",
        )
        self.assertFalse(
            summary.endswith("."),
            f"Bioconda about.summary must not end with a period: {summary!r}",
        )
        self.assertTrue(summary.isascii(), f"Bioconda about.summary must be pure ASCII: {summary!r}")
        self.assertIn("VarDictJava", summary)
        self.assertNotIn("--compat", summary)
        self.assertIn("--compat=off", meta)
        self.assertIn("biotools:vardictrust", meta, "meta.yaml must list biotools:vardictrust")
        self.assertIsNone(
            re.search(r"\bbiotools:vardict\b", meta),
            "meta.yaml must replace biotools:vardict with biotools:vardictrust",
        )
        self.assertNotIn("ultra-fast", meta.lower(), "meta.yaml must not contain 'ultra-fast'")

    def test_readme_bioconda_package_name(self):
        readme_path = ROOT_DIR / "README.md"
        content = readme_path.read_text(encoding="utf-8")
        self.assertIn(
            "conda install -c bioconda -c conda-forge vardictrust",
            content,
            "README.md must document installation using official package name 'vardictrust'",
        )
        self.assertNotIn(
            "vardict-rust",
            content,
            "README.md must NOT reference erroneous package name 'vardict-rust'",
        )

    def test_conda_platform_identifiers(self):
        readme = (ROOT_DIR / "README.md").read_text(encoding="utf-8")
        changelog = (ROOT_DIR / "CHANGELOG.md").read_text(encoding="utf-8")
        for doc, name in [(readme, "README.md"), (changelog, "CHANGELOG.md")]:
            self.assertIn("osx-arm64", doc, f"{name} must use valid Conda platform identifier 'osx-arm64'")
            self.assertIn("osx-64", doc, f"{name} must use valid Conda platform identifier 'osx-64'")
            self.assertIn("linux-64", doc, f"{name} must use valid Conda platform identifier 'linux-64'")
            self.assertIn("linux-aarch64", doc, f"{name} must use valid Conda platform identifier 'linux-aarch64'")
            self.assertNotIn("macOS-arm64", doc, f"{name} must not use invalid Conda subdir 'macOS-arm64'")
            self.assertNotIn("macOS-64", doc, f"{name} must not use invalid Conda subdir 'macOS-64'")


class TestNfCoreModulePackaging(unittest.TestCase):
    """Verify in-tree Nextflow nf-core module structure, review fixes, and lockfiles."""

    def test_nfcore_module_files_exist(self):
        required = [
            NFCORE_DIR / "prepare_nfcore_pr.sh",
            NFCORE_MOD_DIR / "main.nf",
            NFCORE_MOD_DIR / "meta.yml",
            NFCORE_MOD_DIR / "environment.yml",
            NFCORE_MOD_DIR / "tests" / "main.nf.test",
            NFCORE_MOD_DIR / "tests" / "main.nf.test.snap",
            NFCORE_MOD_DIR / "tests" / "nextflow.config",
        ]
        for path in required:
            self.assertTrue(path.is_file(), f"Missing nf-core packaging file: {path}")
        self.assertTrue(
            os.access(NFCORE_DIR / "prepare_nfcore_pr.sh", os.X_OK),
            "prepare_nfcore_pr.sh must be executable",
        )
        script_content = (NFCORE_DIR / "prepare_nfcore_pr.sh").read_text(encoding="utf-8")
        self.assertIn(
            '${NFCORE_MODULES_DIR:?',
            script_content,
            "prepare_nfcore_pr.sh must require explicit NFCORE_MODULES_DIR",
        )
        self.assertNotIn(
            "/tmp/opencode/nf-core-modules",
            script_content,
            "prepare_nfcore_pr.sh must not default to world-writable /tmp path",
        )

    def test_nfcore_main_nf_and_meta_args_id(self):
        main_nf = (NFCORE_MOD_DIR / "main.nf").read_text(encoding="utf-8")
        meta_yml = (NFCORE_MOD_DIR / "meta.yml").read_text(encoding="utf-8")

        self.assertIn("def args = task.ext.args", main_nf)
        self.assertIn("def args2 = task.ext.args2", main_nf)
        self.assertNotIn("args3", main_nf, "main.nf must not define redundant args3")
        self.assertIn("| bgzip ${args2}", main_nf)

        self.assertIn('args_id: "$args"', meta_yml)
        self.assertIn('args_id: "$args2"', meta_yml)
        self.assertNotIn('args_id: "$args $args2"', meta_yml)
        self.assertNotIn('args_id: "$args3"', meta_yml)
        self.assertIn('documentation: "https://vardictrust.readthedocs.io/"', meta_yml)
        self.assertIn('pattern: "*.{bam,cram}"', meta_yml)
        self.assertIn("http://edamontology.org/format_3462", meta_yml)
        self.assertIn('pattern: "*.{bai,csi,crai}"', meta_yml)
        meta_lines = meta_yml.splitlines()
        self.assertTrue(
            meta_lines[1].startswith('description: "') and meta_lines[1].endswith('"'),
            f"nf-core meta.yml:2 description must be double-quoted: {meta_lines[1]!r}",
        )
        top_desc = meta_lines[1][len('description: "') : -1]
        self.assertLessEqual(
            len(top_desc),
            80,
            f"nf-core meta.yml:2 description must be <= 80 chars (got {len(top_desc)}): {top_desc!r}",
        )
        self.assertFalse(
            top_desc.endswith("."),
            f"nf-core meta.yml:2 description must not end with a period: {top_desc!r}",
        )
        self.assertIn("VarDictJava", top_desc)
        self.assertNotIn("--compat", top_desc)
        tool_block_m = re.search(
            r'^\s*-\s*"vardictrust":\s*\n\s*description:\s*("[^\n]+")',
            meta_yml,
            re.MULTILINE,
        )
        self.assertIsNotNone(
            tool_block_m,
            'nf-core meta.yml tools[0]."vardictrust".description must be double-quoted',
        )
        tool_desc = tool_block_m.group(1)
        self.assertIn("--compat=off", tool_desc)
        self.assertIn("--vcf-all-alleles", tool_desc)
        self.assertNotIn("--output", tool_desc)
        self.assertNotIn("--tbi", tool_desc)
        self.assertIn("description: BED with the regions of interest", meta_yml)

    def test_nfcore_tests_use_lowercase_channel_factory(self):
        nf_test = (NFCORE_MOD_DIR / "tests" / "main.nf.test").read_text(encoding="utf-8")
        self.assertNotIn("Channel.value(", nf_test, "main.nf.test must use lowercase channel.value(")
        self.assertEqual(
            nf_test.count("channel.value(["),
            4,
            "All 4 test blocks in main.nf.test must use channel.value([",
        )

    def test_nfcore_environment_and_conda_lockfiles(self):
        env_yml = (NFCORE_MOD_DIR / "environment.yml").read_text(encoding="utf-8")
        self.assertRegex(env_yml, r"bioconda::vardictrust=\d+\.\d+\.\d+")
        self.assertRegex(env_yml, r"bioconda::htslib=\d+\.\d+")

        meta_yml = (NFCORE_MOD_DIR / "meta.yml").read_text(encoding="utf-8")
        referenced_locks = re.findall(
            r"lock_file:\s*modules/nf-core/vardictrust/(\.conda-lock/[^\s]+)",
            meta_yml,
        )
        self.assertEqual(len(referenced_locks), 2, "meta.yml must reference amd64 and arm64 lock_files")
        expected_names = set()
        for rel_lock in referenced_locks:
            lock_path = NFCORE_MOD_DIR / rel_lock
            self.assertTrue(lock_path.is_file(), f"Missing .conda-lock file: {lock_path}")
            expected_names.add(lock_path.name)

        actual_names = {p.name for p in (NFCORE_MOD_DIR / ".conda-lock").glob("*.txt")}
        self.assertEqual(
            actual_names,
            expected_names,
            f"Orphaned .conda-lock files detected: {actual_names ^ expected_names}",
        )

    def test_readme_container_uris_match_nfcore_meta(self):
        readme = (ROOT_DIR / "README.md").read_text(encoding="utf-8")
        meta_yml = (NFCORE_MOD_DIR / "meta.yml").read_text(encoding="utf-8")
        self.assertNotIn(
            "ghcr.io/palakpsheth/vardictrust",
            readme,
            "README.md must not reference unpublished ghcr.io/palakpsheth/vardictrust images",
        )
        wave_uris = re.findall(
            r"(?:oras://)?community\.wave\.seqera\.io/library/htslib_vardictrust:[0-9a-f]{16}",
            readme,
        )
        self.assertGreaterEqual(
            len(wave_uris),
            2,
            "README.md must document verified multi-arch Seqera Wave container URIs",
        )
        for uri in wave_uris:
            self.assertIn(
                uri,
                meta_yml,
                f"README.md Wave container URI {uri!r} not found in packaging/nf-core/modules/vardictrust/meta.yml",
            )


class TestRepoAgentSkills(unittest.TestCase):
    """Verify .agents/skills/ SDO frontmatter, markdown hygiene, and release consistency validator."""

    def test_all_four_skills_structure_and_sdo_frontmatter(self):
        for skill_name in EXPECTED_SKILLS:
            skill_md = SKILLS_DIR / skill_name / "SKILL.md"
            self.assertTrue(skill_md.is_file(), f"Missing skill file: {skill_md}")
            text = skill_md.read_text(encoding="utf-8")
            lines = text.splitlines()
            self.assertLess(
                len(lines),
                350,
                f"{skill_md} has {len(lines)} lines (must stay < 350 lines)",
            )
            self.assertTrue(text.startswith("---\n"), f"{skill_md} must begin with YAML frontmatter")
            parts = text.split("---\n", 2)
            self.assertEqual(len(parts), 3, f"{skill_md} must have closing --- frontmatter delimiter")
            fm = parts[1]
            name_m = re.search(r"^name:\s*(\S+)\s*$", fm, re.MULTILINE)
            desc_m = re.search(r"^description:\s*(.+?)\s*$", fm, re.MULTILINE)
            self.assertIsNotNone(name_m, f"{skill_md} missing name field")
            self.assertEqual(name_m.group(1), skill_name)
            self.assertIsNotNone(desc_m, f"{skill_md} missing description field")
            desc = desc_m.group(1).strip()
            self.assertTrue(
                desc.startswith("Use when "),
                f"{skill_md} description must start with 'Use when ': {desc!r}",
            )
            self.assertLess(
                len(desc),
                250,
                f"{skill_md} description is {len(desc)} chars (must be < 250 chars)",
            )

    def test_skill_markdown_tables_and_angle_brackets(self):
        angle_pattern = re.compile(r"<([A-Za-z0-9_\-]+)>")
        for md_path in sorted(SKILLS_DIR.rglob("*.md")):
            rel_path = md_path.relative_to(ROOT_DIR)
            raw = md_path.read_text(encoding="utf-8")
            body = raw.split("---\n", 2)[2] if raw.startswith("---\n") else raw

            # 1. Check angle-bracket escaping outside code spans/blocks
            no_code_blocks = re.sub(r"```[\s\S]*?```", "", body)
            no_inline_code = re.sub(r"`[^`]*`", "", no_code_blocks)
            for match in angle_pattern.finditer(no_inline_code):
                tag = match.group(1).lower()
                if tag not in ALLOWED_HTML_TAGS and "@" not in tag and "." not in tag:
                    self.fail(
                        f"{rel_path}: unescaped angle bracket token '<{match.group(1)}>' must be backtick-escaped"
                    )

    def test_pr_review_governance_contract(self):
        surfaces = (
            ROOT_DIR / "AGENTS.md",
            ROOT_DIR / "CONTRIBUTING.md",
            ROOT_DIR / ".github" / "pull_request_template.md",
            SKILLS_DIR / "vardictrust-pr" / "SKILL.md",
        )
        required_tokens = (
            "/ponytail-review",
            "/rust-router",
            "/pragmatic-rust-guidelines",
            "copilot-pull-request-reviewer[bot]",
            "resolveReviewThread",
        )
        for surface in surfaces:
            self.assertTrue(surface.is_file(), f"Missing PR governance surface: {surface}")
            content = surface.read_text(encoding="utf-8")
            for token in required_tokens:
                self.assertIn(
                    token,
                    content,
                    f"{surface.relative_to(ROOT_DIR)} must document mandatory PR review requirement {token!r}",
                )
        skill_content = (SKILLS_DIR / "vardictrust-pr" / "SKILL.md").read_text(encoding="utf-8")
        contrib_content = (ROOT_DIR / "CONTRIBUTING.md").read_text(encoding="utf-8")
        self.assertIn("gh api --paginate repos/palakpsheth/VarDictRust/pulls/<PR>/reviews", skill_content)
        self.assertIn("gh api --paginate repos/palakpsheth/VarDictRust/pulls/<PR>/comments", skill_content)
        self.assertIn("gh api --paginate repos/palakpsheth/VarDictRust/commits/<HEAD_SHA>/check-runs", skill_content)
        self.assertIn("gh api graphql --paginate", skill_content)
        self.assertIn("pageInfo { hasNextPage endCursor }", skill_content)
        self.assertIn("--subject", skill_content)
        self.assertIn("--body-file", skill_content)
        self.assertIn("gh api --paginate repos/palakpsheth/VarDictRust/pulls/<PR>/reviews", contrib_content)
        pr_template_content = (ROOT_DIR / ".github" / "pull_request_template.md").read_text(encoding="utf-8")
        self.assertIn(
            "python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 hg008-cram-twin",
            pr_template_content,
        )

    def test_release_consistency_script_passes_and_catches_drift(self):
        script = SKILLS_DIR / "vardictrust-release" / "scripts" / "check_release_consistency.py"
        self.assertTrue(script.is_file(), f"Missing validator script: {script}")
        self.assertTrue(os.access(script, os.X_OK), "check_release_consistency.py must be executable")

        cargo_cfg = toml_lite.load_path(ROOT_DIR / "Cargo.toml")
        current_version = ((cargo_cfg.get("workspace") or {}).get("package") or {}).get("version")
        self.assertTrue(isinstance(current_version, str) and current_version, "Could not parse workspace version from Cargo.toml")

        # Positive check on current repo root
        res = subprocess.run(
            [sys.executable, str(script), "--version", current_version, "--phase", "pre-tag"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(
            res.returncode,
            0,
            f"check_release_consistency.py failed on {current_version}:\n{res.stderr}\n{res.stdout}",
        )

        # Negative check: asking for a mismatched version must fail non-zero
        res_bad = subprocess.run(
            [sys.executable, str(script), "--version", "9.9.9", "--phase", "pre-tag"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(res_bad.returncode, 0, "Validator must fail when version mismatches")
        self.assertIn("RELEASE CONSISTENCY FAILED", res_bad.stderr)

        # Negative check: missing `cargo` in PATH with --verify-cargo-list records actionable error (no traceback)
        res_nocargo = subprocess.run(
            [
                sys.executable,
                str(script),
                "--version",
                current_version,
                "--phase",
                "pre-tag",
                "--verify-cargo-list",
            ],
            env={"PATH": ""},
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(res_nocargo.returncode, 1)
        self.assertIn("RELEASE CONSISTENCY FAILED", res_nocargo.stderr)
        self.assertIn("could not execute 'cargo'", res_nocargo.stderr)
        self.assertNotIn("Traceback", res_nocargo.stderr)

        # Negative check: unregistered crates/*/tests/*.rs target is caught
        import importlib.util
        import tempfile

        spec = importlib.util.spec_from_file_location("check_release_consistency", script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_root = pathlib.Path(tmpdir)
            for pkg, targets in mod.REGISTERED_INTEGRATION_TESTS.items():
                tdir = fake_root / "crates" / pkg / "tests"
                tdir.mkdir(parents=True)
                for tgt in targets:
                    (tdir / f"{tgt}.rs").write_text("", encoding="utf-8")
            (fake_root / "crates" / "vardictrust-core" / "tests" / "sneaky_extra.rs").write_text(
                "", encoding="utf-8"
            )
            errs: list[str] = []
            self.assertFalse(mod.check_registered_integration_tests(fake_root, errs))
            self.assertTrue(any("sneaky_extra" in e for e in errs), f"Expected sneaky_extra error, got {errs}")

            # Certified release tree (no Current Workspace Test Suite block; counts in Verified Parity Baseline)
            # must pass both without --require-head-fingerprint (pre-tag CI on main / post-tag / post-bioconda / post-nfcore)
            # and with --require-head-fingerprint, while rejecting leftover dev blocks in post-tag/require_head_fp modes.
            tc = mod.EXPECTED_TEST_COUNTS
            ev_dir = fake_root / "benchmarks" / "evidence"
            ev_dir.mkdir(parents=True)
            (ev_dir / "hg008-tier2.json").write_text('{"fingerprint": "deadbeefcafe1234"}\n', encoding="utf-8")
            certified_div = (
                "# DIVERGENCES\n"
                "##source=VarDict_v1.8.2\n"
                f"Verified Parity Baseline (release/v1.0.8, commit abc1234):\n"
                f"- Integration Tests:  91 PASS / 1 XFAIL\n"
                f"- Core Library Tests: {tc['core_lib']} PASS / 0 FAIL\n"
                f"- Total Workspace:    {tc['total_workspace']} PASS / 0 FAIL\n"
                "- Final Acceptance:   8/8 PASS under release engine fingerprint `deadbeefcafe1234`\n"
            )
            (fake_root / "DIVERGENCES.md").write_text(certified_div, encoding="utf-8")
            for ph, req_fp in (("pre-tag", False), ("pre-tag", True), ("post-tag", False), ("post-bioconda", False), ("post-nfcore", False)):
                div_errs: list[str] = []
                mod.check_divergences_surface(fake_root, "1.0.8", div_errs, phase=ph, require_head_fp=req_fp)
                self.assertEqual(div_errs, [], f"Certified DIVERGENCES.md rejected for phase={ph}, require_head_fp={req_fp}: {div_errs}")

            unfolded_div = certified_div + (
                f"\nCurrent Workspace Test Suite (main / pre-v1.0.9):\n"
                f"- Core Library Tests: {tc['core_lib']} PASS / 0 FAIL\n"
                f"- Total Workspace:    {tc['total_workspace']} PASS / 0 FAIL\n"
            )
            (fake_root / "DIVERGENCES.md").write_text(unfolded_div, encoding="utf-8")
            post_tag_errs: list[str] = []
            mod.check_divergences_surface(fake_root, "1.0.8", post_tag_errs, phase="post-tag", require_head_fp=False)
            self.assertTrue(post_tag_errs, "post-tag must reject unfolded Current Workspace Test Suite block")

            sec_dir = fake_root / ".github"
            sec_dir.mkdir(parents=True, exist_ok=True)
            (sec_dir / "SECURITY.md").write_text(
                "| Version | Supported |\n| :--- | :--- |\n| **`1.1.x`** | :white_check_mark: |\n| `< 1.1.0` | :x: |\n",
                encoding="utf-8",
            )
            sec_ok_errs: list[str] = []
            mod.check_security_policy_surface(fake_root, "1.1.0", sec_ok_errs)
            self.assertEqual(sec_ok_errs, [], f"Valid 1.1.x SECURITY.md rejected: {sec_ok_errs}")

            sec_stale_errs: list[str] = []
            mod.check_security_policy_surface(fake_root, "1.2.0", sec_stale_errs)
            self.assertTrue(sec_stale_errs, "Stale 1.1.x SECURITY.md must fail for 1.2.0")

            (sec_dir / "SECURITY.md").write_text(
                "| Version | Supported |\n| :--- | :--- |\n| **`1.1.x`** | :white_check_mark: |\n| **`1.0.x`** | :white_check_mark: |\n",
                encoding="utf-8",
            )
            sec_dup_errs: list[str] = []
            mod.check_security_policy_surface(fake_root, "1.1.0", sec_dup_errs)
            self.assertTrue(sec_dup_errs, "Multiple :white_check_mark: rows in SECURITY.md must fail")

            (sec_dir / "SECURITY.md").write_text("", encoding="utf-8")
            sec_empty_errs: list[str] = []
            mod.check_security_policy_surface(fake_root, "1.1.0", sec_empty_errs)
            self.assertTrue(sec_empty_errs, "Empty 0-byte SECURITY.md must fail check_security_policy_surface")

            (sec_dir / "SECURITY.md").write_text(
                "| Version | Supported |\n| :--- | :--- |\n| **`1.1.x`** | :white_check_mark: |\n| `< 1.0.0` | :x: |\n",
                encoding="utf-8",
            )
            sec_stale_floor_errs: list[str] = []
            mod.check_security_policy_surface(fake_root, "1.1.0", sec_stale_floor_errs)
            self.assertTrue(
                sec_stale_floor_errs,
                "1.1.x active row paired with stale < 1.0.0 floor row must fail for 1.1.0",
            )


class TestLintRatchet(unittest.TestCase):
    """Verify 11-lint safety/hygiene + unsafe_code zero-allow ratchet across both crates and CI (AR-1)."""

    EXPECTED_11_LINTS = (
        "clippy::unwrap_used",
        "clippy::expect_used",
        "clippy::panic",
        "clippy::unreachable",
        "clippy::todo",
        "clippy::unimplemented",
        "clippy::dbg_macro",
        "clippy::print_stdout",
        "clippy::print_stderr",
        "clippy::string_slice",
        "clippy::indexing_slicing",
    )

    @staticmethod
    def _normalize_ws(text: str) -> str:
        return re.sub(r"\s+", "", text)

    def _populate_clean_workspace(self, root: pathlib.Path) -> None:
        shutil.rmtree(root / "crates", ignore_errors=True)
        core_src = root / "crates" / "vardictrust-core" / "src"
        cli_src = root / "crates" / "vardictrust" / "src"
        core_src.mkdir(parents=True, exist_ok=True)
        cli_src.mkdir(parents=True, exist_ok=True)

        lints_csv = ",\n        ".join(self.EXPECTED_11_LINTS)
        (core_src / "lib.rs").write_text(
            "#![forbid(unsafe_code)]\n"
            f"#![cfg_attr(\n    not(test),\n    deny(\n        {lints_csv},\n    )\n)]\n"
            "pub fn ok() {}\n",
            encoding="utf-8",
        )
        (cli_src / "cli.rs").write_text(
            "#![forbid(unsafe_code)]\n"
            "pub fn ok() {}\n",
            encoding="utf-8",
        )
        (cli_src / "profile.rs").write_text(
            '#![cfg_attr(not(feature = "profile-cpu"), forbid(unsafe_code))]\n'
            '#![cfg_attr(\n'
            '    feature = "profile-cpu",\n'
            '    deny(\n'
            '        unsafe_op_in_unsafe_fn,\n'
            '        clippy::undocumented_unsafe_blocks,\n'
            '        clippy::multiple_unsafe_ops_per_block\n'
            '    )\n'
            ')]\n'
            "pub fn ok() {}\n",
            encoding="utf-8",
        )
        (cli_src / "main.rs").write_text(
            '#![cfg_attr(\n'
            '    not(any(\n'
            '        feature = "profile-cpu",\n'
            '        feature = "jemalloc",\n'
            '        feature = "profile-heap"\n'
            '    )),\n'
            '    forbid(unsafe_code)\n'
            ')]\n'
            '#![cfg_attr(\n'
            '    all(\n'
            '        any(feature = "jemalloc", feature = "profile-heap"),\n'
            '        not(feature = "profile-cpu")\n'
            '    ),\n'
            '    deny(unsafe_code)\n'
            ')]\n'
            '#![cfg_attr(\n'
            '    feature = "profile-cpu",\n'
            '    deny(\n'
            '        unsafe_op_in_unsafe_fn,\n'
            '        clippy::undocumented_unsafe_blocks,\n'
            '        clippy::multiple_unsafe_ops_per_block\n'
            '    )\n'
            ')]\n'
            f"#![cfg_attr(\n    not(test),\n    deny(\n        {lints_csv},\n    )\n)]\n"
            "mod cli;\n"
            "mod profile;\n"
            '#[cfg(all(feature = "jemalloc", not(feature = "profile-heap"), target_os = "linux"))]\n'
            '#[cfg_attr(not(feature = "profile-cpu"), expect(unsafe_code, reason = "jemalloc C ABI symbol _rjem_malloc_conf requires #[export_name] on thin &[u8; N] static"))]\n'
            '#[export_name = "_rjem_malloc_conf"]\n'
            'pub static MALLOC_CONF: &[u8; 4] = b"ok!\\0";\n'
            "fn main() {}\n",
            encoding="utf-8",
        )

    def test_crate_headers_and_ci_wiring(self):
        lib_rs = ROOT_DIR / "crates" / "vardictrust-core" / "src" / "lib.rs"
        main_rs = ROOT_DIR / "crates" / "vardictrust" / "src" / "main.rs"
        cli_rs = ROOT_DIR / "crates" / "vardictrust" / "src" / "cli.rs"
        profile_rs = ROOT_DIR / "crates" / "vardictrust" / "src" / "profile.rs"
        pr_check = ROOT_DIR / "harness" / "pr_check.sh"
        ci_yml = ROOT_DIR / ".github" / "workflows" / "ci.yml"

        lib_content = lib_rs.read_text(encoding="utf-8")
        main_content = main_rs.read_text(encoding="utf-8")
        cli_content = cli_rs.read_text(encoding="utf-8")
        profile_content = profile_rs.read_text(encoding="utf-8")

        # 1) Both lib.rs and main.rs must deny all 11 safety/hygiene lints in non-test builds
        for target_name, content in (("vardictrust-core/src/lib.rs", lib_content), ("vardictrust/src/main.rs", main_content)):
            for lint in self.EXPECTED_11_LINTS:
                self.assertTrue(lint in content, f"{target_name} missing deny({lint})")

        # 2) Exact unsafe_code & unsafe discipline headers across lib.rs, cli.rs, profile.rs, and main.rs
        self.assertTrue("#![forbid(unsafe_code)]" in lib_content, "lib.rs missing #![forbid(unsafe_code)]")
        self.assertTrue("#![forbid(unsafe_code)]" in cli_content, "cli.rs missing #![forbid(unsafe_code)]")

        profile_norm = self._normalize_ws(profile_content)
        self.assertTrue(
            self._normalize_ws('#![cfg_attr(not(feature = "profile-cpu"), forbid(unsafe_code))]') in profile_norm,
            "profile.rs missing not(feature = \"profile-cpu\") forbid(unsafe_code)",
        )
        self.assertTrue(
            self._normalize_ws(
                '#![cfg_attr(feature = "profile-cpu", deny(unsafe_op_in_unsafe_fn, clippy::undocumented_unsafe_blocks, clippy::multiple_unsafe_ops_per_block))]'
            ) in profile_norm,
            "profile.rs missing feature = \"profile-cpu\" unsafe discipline deny header",
        )

        main_norm = self._normalize_ws(main_content)
        self.assertTrue(
            self._normalize_ws(
                '#![cfg_attr(not(any(feature = "profile-cpu", feature = "jemalloc", feature = "profile-heap")), forbid(unsafe_code))]'
            ) in main_norm,
            "main.rs missing forbid(unsafe_code) feature guard",
        )
        self.assertTrue(
            self._normalize_ws(
                '#![cfg_attr(all(any(feature = "jemalloc", feature = "profile-heap"), not(feature = "profile-cpu")), deny(unsafe_code))]'
            ) in main_norm,
            "main.rs missing deny(unsafe_code) feature guard",
        )
        # All 4 _rjem_malloc_conf statics carry cfg_attr(not(feature = "profile-cpu"), expect(unsafe_code, ...))
        self.assertEqual(
            main_content.count('#[export_name = "_rjem_malloc_conf"]'),
            4,
            "Expected 4 target/feature-split _rjem_malloc_conf statics in main.rs",
        )
        self.assertEqual(
            len(re.findall(r"\bexpect\s*\(\s*unsafe_code\b", main_content)),
            4,
            "expect(unsafe_code) must appear solely on the 4 _rjem_malloc_conf statics in main.rs",
        )

        # 3) pr_check.sh and .github/workflows/ci.yml Step 1 both run lint_ratchet.sh --scan-only
        #    plus all three Clippy feature configurations (--no-default-features, default, --all-features)
        pr_check_text = pr_check.read_text(encoding="utf-8")
        self.assertIn(
            "bash harness/lint_ratchet.sh --scan-only",
            pr_check_text,
        )
        ci_text = ci_yml.read_text(encoding="utf-8")
        step1_start = ci_text.find("1. Format & Lint Gate")
        step2_start = ci_text.find("2. Build Instrumented Workspace")
        self.assertNotEqual(step1_start, -1, "ci.yml missing Step 1")
        self.assertGreater(step2_start, step1_start)
        step1_block = ci_text[step1_start:step2_start]
        self.assertIn(
            "bash harness/lint_ratchet.sh --scan-only",
            step1_block,
            "ci.yml Step 1 must execute 'bash harness/lint_ratchet.sh --scan-only'",
        )
        for clippy_cmd in (
            "cargo clippy --workspace --all-targets --no-default-features -- -D warnings",
            "cargo clippy --workspace --all-targets -- -D warnings",
            "cargo clippy --all-targets --all-features -- -D warnings",
        ):
            self.assertIn(
                clippy_cmd,
                step1_block,
                f"ci.yml Step 1 must execute {clippy_cmd!r}",
            )
            self.assertIn(
                clippy_cmd,
                pr_check_text,
                f"pr_check.sh Stage 1 must execute {clippy_cmd!r}",
            )

    def test_lint_ratchet_scan_only_and_unauthorized_suppression_rejection(self):
        ratchet_script = ROOT_DIR / "harness" / "lint_ratchet.sh"
        self.assertTrue(ratchet_script.is_file(), f"Missing lint ratchet script: {ratchet_script}")
        self.assertTrue(os.access(ratchet_script, os.X_OK), "harness/lint_ratchet.sh must be executable")

        # 1) Live --scan-only check on workspace must exit 0 quickly
        t0 = time.monotonic()
        scan_res = subprocess.run(
            ["bash", str(ratchet_script), "--scan-only"],
            cwd=str(ROOT_DIR),
            capture_output=True,
            text=True,
            check=False,
        )
        elapsed = time.monotonic() - t0
        self.assertEqual(
            scan_res.returncode,
            0,
            f"lint_ratchet.sh --scan-only failed ({scan_res.returncode}):\n{scan_res.stderr}\n{scan_res.stdout}",
        )
        self.assertLess(elapsed, 2.0, f"lint_ratchet.sh --scan-only took {elapsed:.3f}s (expected < 0.5s)")

        # 2) Unknown flag must fail non-zero
        bad_flag = subprocess.run(
            ["bash", str(ratchet_script), "--unknown-flag"],
            cwd=str(ROOT_DIR),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(bad_flag.returncode, 0, "lint_ratchet.sh must reject unknown flags")

        # 3) Clean synthetic workspace passes, and injecting any unauthorized allow/expect in either crate fails
        unauthorized_attrs = (
            "#[allow(clippy::unwrap_used)]",
            "#[expect(unsafe_code)]",
            "#[allow(clippy::print_stderr)]",
            "#[allow(clippy::undocumented_unsafe_blocks)]",
            "#[allow(clippy::indexing_slicing)]",
        )
        injection_targets = (
            "crates/vardictrust-core/src/sneaky.rs",
            "crates/vardictrust/src/sneaky_cli.rs",
            "crates/vardictrust/src/main.rs",
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = pathlib.Path(tmp_dir)
            self._populate_clean_workspace(tmp_root)
            env = dict(os.environ, LINT_RATCHET_ROOT=tmp_dir)

            clean_res = subprocess.run(
                ["bash", str(ratchet_script), "--scan-only"],
                cwd=str(ROOT_DIR),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                clean_res.returncode,
                0,
                f"Clean synthetic workspace failed --scan-only:\n{clean_res.stderr}\n{clean_res.stdout}",
            )

            for rel_target in injection_targets:
                for bad_attr in unauthorized_attrs:
                    self._populate_clean_workspace(tmp_root)
                    target_path = tmp_root / rel_target
                    existing = target_path.read_text(encoding="utf-8") if target_path.exists() else ""
                    target_path.write_text(
                        existing
                        + f"\n{bad_attr}\npub fn unauthorized_bypass() {{}}\n"
                        + "mod tests {\n    #![cfg(test)]\n    #[allow(clippy::unwrap_used)]\n    fn ok_in_test() {}\n}\n",
                        encoding="utf-8",
                    )
                    viol_res = subprocess.run(
                        ["bash", str(ratchet_script), "--scan-only"],
                        cwd=str(ROOT_DIR),
                        env=env,
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    self.assertNotEqual(
                        viol_res.returncode,
                        0,
                        f"Unauthorized {bad_attr} in {rel_target} must fail lint_ratchet.sh --scan-only",
                    )
                    self.assertIn(rel_target, viol_res.stderr)

            # 4) Stripping required unsafe_code headers from cli.rs, profile.rs, or main.rs must also fail
            for broken_rel, broken_src in (
                ("crates/vardictrust/src/cli.rs", "pub fn missing_forbid() {}\n"),
                ("crates/vardictrust/src/profile.rs", "#![deny(unsafe_op_in_unsafe_fn)]\npub fn missing_forbid() {}\n"),
                (
                    "crates/vardictrust/src/main.rs",
                    "#![cfg_attr(not(test), deny(" + ", ".join(self.EXPECTED_11_LINTS) + "))]\nfn main() {}\n",
                ),
            ):
                self._populate_clean_workspace(tmp_root)
                (tmp_root / broken_rel).write_text(broken_src, encoding="utf-8")
                hdr_res = subprocess.run(
                    ["bash", str(ratchet_script), "--scan-only"],
                    cwd=str(ROOT_DIR),
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertNotEqual(
                    hdr_res.returncode,
                    0,
                    f"Missing unsafe_code header in {broken_rel} must fail lint_ratchet.sh --scan-only",
                )
                self.assertIn(broken_rel, hdr_res.stderr)


if __name__ == "__main__":
    unittest.main()
