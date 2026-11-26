#!/usr/bin/env python3
"""
harness/tests/test_release_workflows.py

Automated verification of VarDictRust CI/CD workflows, cargo-dist multi-stage
tag checkout, SLSA v1 provenance predicate/attestation verification, SHA-256
cargo-dist bootstrap, crates.io idempotent publishing, and GitHub repository
tag rulesets. Conforms strictly to the Python Standard Library purity constraint.
"""

import base64
import hashlib
import http.server
import importlib
import importlib.util
import io
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import textwrap
import threading
import unittest
import urllib.error

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent
if str(ROOT_DIR / "harness") not in sys.path:
    sys.path.insert(0, str(ROOT_DIR / "harness"))

import toml_lite  # noqa: E402


class TestReleaseAndCiWorkflows(unittest.TestCase):
    """Verify GitHub Actions CI/release workflows, SLSA attestation, and repo rulesets."""

    def test_ci_workflow_pipeline_and_concurrency_contract(self):
        workflows_dir = ROOT_DIR / ".github" / "workflows"
        for wf_name in ("ci.yml", "differential.yml", "docs.yml", "security.yml"):
            wf_path = workflows_dir / wf_name
            self.assertTrue(wf_path.is_file(), f"Missing workflow: {wf_path}")
            wf_text = wf_path.read_text(encoding="utf-8")
            self.assertIn("concurrency:", wf_text, f"{wf_name} missing concurrency block")
            self.assertIn(
                "cancel-in-progress: ${{ github.ref != 'refs/heads/main' }}",
                wf_text,
                f"{wf_name} missing cancel-in-progress guard",
            )
            if wf_name in ("ci.yml", "differential.yml", "docs.yml"):
                self.assertIn("workflow_dispatch:", wf_text, f"{wf_name} missing workflow_dispatch")

        ci_text = (workflows_dir / "ci.yml").read_text(encoding="utf-8")
        self.assertIn('RUST_TEST_THREADS: "2"', ci_text)
        self.assertIn('RAYON_NUM_THREADS: "2"', ci_text)
        self.assertIn('RUST_THREADS: "1"', ci_text)
        self.assertNotIn("apt-get install -y --no-install-recommends samtools", ci_text)
        self.assertIn("build-and-test:", ci_text)
        self.assertIn("name: Build, Unit, Integration, Parity & Coverage Gate", ci_text)

        expected_ci_steps = [
            "1. Format & Lint Gate (rustfmt & clippy)",
            "2. Build Instrumented Workspace & Test Targets",
            "3. Unit Tests (vardictrust-core & vardictrust)",
            "4. Integration & Hostile-Input Tests (Modes, Filters, CRAM & Bounded Memory)",
            "5. Corpus & Modality Parity Tests (110 Golden & Modality Cases)",
            "6. Code Coverage Gate (cargo-llvm-cov >= 90% Lines)",
            "7. Harness Diagnostic Self-Tests",
            "8. Release Consistency, Performance Budget & Cleanliness Gate",
            "9. Publish CI Telemetry Summary",
        ]
        last_pos = -1
        for step_label in expected_ci_steps:
            pos = ci_text.find(step_label)
            self.assertNotEqual(pos, -1, f"ci.yml missing step label: {step_label!r}")
            self.assertGreater(pos, last_pos, f"ci.yml step out of order: {step_label!r}")
            last_pos = pos

        self.assertIn("cargo llvm-cov clean --workspace", ci_text)
        self.assertIn("cargo llvm-cov --no-report --workspace --bins --tests -- --list >/dev/null", ci_text)
        self.assertIn("cargo llvm-cov --no-report --workspace --lib --bins -- --test-threads=2", ci_text)
        self.assertNotIn("--no-report --no-clean", ci_text)
        self.assertIn("-p vardictrust-core --test engine_behavior --test amplicon_stream_mem --test cram_stack", ci_text)
        self.assertIn("-p vardictrust-core --test corpus_parity -- --test-threads=2", ci_text)
        self.assertIn("cargo llvm-cov report --ignore-filename-regex '(^|/)tests/' --fail-under-lines 90", ci_text)
        self.assertIn("harness/tests/test_git_hooks.sh", ci_text)
        self.assertIn("harness/tests/test_release_workflows.py", ci_text)
        self.assertIn("--verify-cargo-list --cargo-llvm-cov", ci_text)

        diff_text = (workflows_dir / "differential.yml").read_text(encoding="utf-8")
        self.assertIn('RUST_TEST_THREADS: "2"', diff_text)
        self.assertIn('RAYON_NUM_THREADS: "2"', diff_text)
        self.assertIn("name: Dual-Engine Differential Parity (Tier-2)", diff_text)
        self.assertNotIn("paths:", diff_text, "differential.yml must not filter paths or required PR checks hang on docs PRs")
        expected_diff_steps = [
            "1. Build Release Binary & Verify Pinned Golden Store",
            "2. Tier-2 Micro Differential Parity (3-Way Java / Rust / C++)",
            "3. Tier-2 Standard Slices Differential Parity (10 Slices: s1a..s9)",
            "4. Tier-2 Extended Modality Differential Parity (17 Modes)",
            "5. Level 6 Oracle & BAM Edge-Case Self-Tests",
            "6. Publish Differential Telemetry Summary",
        ]
        last_diff_pos = -1
        for step_label in expected_diff_steps:
            pos = diff_text.find(step_label)
            self.assertNotEqual(pos, -1, f"differential.yml missing step label: {step_label!r}")
            self.assertGreater(pos, last_diff_pos, f"differential.yml step out of order: {step_label!r}")
            last_diff_pos = pos

        for wf_name in ("ci.yml", "differential.yml", "docs.yml", "security.yml", "release.yml", "wiki.yml"):
            wf_text = (workflows_dir / wf_name).read_text(encoding="utf-8")
            self.assertIn("persist-credentials: false", wf_text, f"{wf_name} missing persist-credentials: false")
            self.assertEqual(
                wf_text.count("persist-credentials: false"),
                wf_text.count("actions/checkout@"),
                f"{wf_name} must set persist-credentials: false on every actions/checkout step",
            )

        wiki_text = (workflows_dir / "wiki.yml").read_text(encoding="utf-8")
        self.assertIn("permissions:\n  contents: read", wiki_text)
        self.assertIn("    permissions:\n      contents: write", wiki_text)
        self.assertIn('echo "::add-mask::${AUTH}"', wiki_text)
        self.assertIn('trap \'rm -rf "${TMP}"\' EXIT', wiki_text)

        sec_text = (workflows_dir / "security.yml").read_text(encoding="utf-8")
        self.assertIn("dtolnay/rust-toolchain@", sec_text)
        self.assertIn("Swatinem/rust-cache@", sec_text)
        self.assertIn("set -euo pipefail", sec_text)
        self.assertIn("cargo deny --all-features check advisories licenses bans sources", sec_text)
        self.assertIn("cargo clippy --all-targets --all-features --message-format=json", sec_text)
        self.assertNotIn("needless_late_init", sec_text)
        self.assertNotIn("continue-on-error: true", sec_text)
        self.assertIn("steps.clippy.outcome == 'success'", sec_text)
        self.assertIn("steps.gitleaks.outcome", sec_text)
        self.assertIn("steps.deny.outcome", sec_text)
        self.assertIn("steps.audit.outcome", sec_text)
        self.assertIn("steps.clippy.outcome", sec_text)

        docs_text = (workflows_dir / "docs.yml").read_text(encoding="utf-8")
        self.assertIn('cache-dependency-path: "docs/requirements.txt"', docs_text)
        self.assertIn("python3 -m pip install -r docs/requirements.txt playwright==1.62.0", docs_text)
        docs_req_text = (ROOT_DIR / "docs" / "requirements.txt").read_text(encoding="utf-8")
        self.assertIn("zensical>=0.1.1,<0.2.0", docs_req_text)
        self.assertNotIn("playwright", docs_req_text)

        conda_action_text = (ROOT_DIR / ".github" / "actions" / "setup-conda-env" / "action.yml").read_text(
            encoding="utf-8"
        )
        self.assertIn("key: conda-env-v3-", conda_action_text)
        self.assertIn("'.github/actions/setup-conda-env/action.yml'", conda_action_text)
        self.assertIn("'harness/env/g4-tools.yml'", conda_action_text)
        self.assertNotIn("restore-keys:", conda_action_text)
        for pkg_pin in (
            "samtools=1.24",
            "bcftools=1.24",
            "htslib=1.24",
            "bedtools=2.31.1",
            "rtg-tools=3.13",
            "bwa=0.7.19",
        ):
            self.assertIn(pkg_pin, conda_action_text, f"setup-conda-env/action.yml missing {pkg_pin}")
        for forbidden_conda_pkg in ("ld_impl_linux-64", "libgcc-ng", "gitleaks"):
            self.assertNotIn(forbidden_conda_pkg, conda_action_text)

        rel_text = (workflows_dir / "release.yml").read_text(encoding="utf-8")
        self.assertIn('permissions:\n  "contents": "read"', rel_text)
        self.assertIn('    permissions:\n      "contents": "write"', rel_text)
        self.assertIn("workflow_dispatch:", rel_text)
        self.assertIn("DISPATCH_TAG: ${{ inputs.tag }}", rel_text)
        self.assertNotIn("${{ inputs.tag }}", rel_text.split("run: |", 1)[1] if "run: |" in rel_text else "")
        self.assertIn('"$WORKFLOW_REF" != "refs/heads/main"', rel_text)
        self.assertIn('RELEASE_COMMIT: "${{ needs.plan.outputs.release-commit }}"', rel_text)
        self.assertIn('RELEASE_TAG: "${{ needs.plan.outputs.tag }}"', rel_text)
        self.assertIn('gh release create "$RELEASE_TAG"', rel_text)
        self.assertIn('git merge-base --is-ancestor "$GITHUB_SHA" refs/remotes/origin/main', rel_text)
        self.assertIn('git merge-base --is-ancestor "$TAG_COMMIT" refs/remotes/origin/main', rel_text)
        self.assertIn(
            'DISALLOWED_DIFF="$(git diff --name-only "$TAG_COMMIT" "$GITHUB_SHA" | '
            "grep -Ev '^(dist-workspace\\.toml|\\.github/workflows/[^/]+\\.yml|harness/tests/test_[^/]+\\.(py|sh)|\\.agents/skills/.+)$' || true)\"",
            rel_text,
        )
        self.assertIn("checkout-ref: ${{ steps.guard.outputs.checkout_ref }}", rel_text)
        self.assertEqual(
            rel_text.count("ref: ${{ needs.plan.outputs.checkout-ref }}"),
            4,
            "All 4 downstream jobs in release.yml must check out needs.plan.outputs.checkout-ref",
        )
        self.assertNotIn("\n  announce:\n", rel_text, "Empty no-op announce job must be removed from release.yml")
        self.assertIn(
            "!fromJson(needs.host.outputs.val).announcement_is_prerelease",
            rel_text,
            "publish-crates if: condition must guard against prerelease tags",
        )
        pub_job_m = re.search(
            r"^  publish-crates:\n((?:  (?![A-Za-z0-9_-]+:).*\n|\n)+?)(?=^  [A-Za-z0-9_-]+:|\Z)",
            rel_text,
            re.MULTILINE,
        )
        self.assertIsNotNone(pub_job_m, "Missing publish-crates job in release.yml")
        pub_job_header = pub_job_m.group(1).split("    steps:\n", 1)[0]
        self.assertNotIn(
            "CARGO_REGISTRY_TOKEN",
            pub_job_header,
            "CARGO_REGISTRY_TOKEN must be step-scoped, not job-scoped on publish-crates",
        )
        self.assertEqual(
            rel_text.count("CARGO_REGISTRY_TOKEN: ${{ secrets.CARGO_REGISTRY_TOKEN }}"),
            2,
            "CARGO_REGISTRY_TOKEN must be scoped to ONLY the two cargo publish steps",
        )
        self.assertNotIn("ref: refs/tags/${{ needs.plan.outputs.tag }}", rel_text)
        self.assertEqual(
            rel_text.count("Verify Exact Tag Checkout Commit"),
            4,
            "build-local-artifacts, build-global-artifacts, host, and publish-crates must assert exact tag commit",
        )
        self.assertIn("name: Generate SLSA Provenance Predicate", rel_text)
        self.assertIn('predicate-type: "https://slsa.dev/provenance/v1"', rel_text)
        self.assertIn('predicate-path: "${{ runner.temp }}/slsa-predicate.json"', rel_text)
        self.assertIn("name: Verify Signed Attestation Bundle", rel_text)
        self.assertIn("plan-wt", rel_text)
        self.assertIn("build-wt", rel_text)
        self.assertIn("--artifacts=lies", rel_text)
        self.assertIn("retired GitHub Actions runner(s)", rel_text)

        dist_cfg = toml_lite.load_path(ROOT_DIR / "dist-workspace.toml")
        custom_runners = (dist_cfg.get("dist") or {}).get("github-custom-runners") or {}
        self.assertEqual(custom_runners.get("aarch64-apple-darwin"), "macos-15")
        dist_toml = (ROOT_DIR / "dist-workspace.toml").read_text(encoding="utf-8")
        for retired_runner in ("macos-12", "macos-13", "macos-14"):
            self.assertNotIn(retired_runner, dist_toml)

        deny_cfg = toml_lite.load_path(ROOT_DIR / "deny.toml")
        self.assertEqual((deny_cfg.get("advisories") or {}).get("yanked"), "deny")
        self.assertEqual((deny_cfg.get("bans") or {}).get("wildcards"), "deny")
        self.assertEqual((deny_cfg.get("sources") or {}).get("required-git-spec"), "rev")

        dependabot_text = (ROOT_DIR / ".github" / "dependabot.yml").read_text(encoding="utf-8")
        for eco in ('package-ecosystem: "cargo"', 'package-ecosystem: "github-actions"', 'package-ecosystem: "pip"'):
            self.assertIn(eco, dependabot_text)
        self.assertEqual(dependabot_text.count("cooldown:"), 3, "dependabot.yml must configure cooldown on all 3 ecosystems")
        self.assertEqual(dependabot_text.count("default-days: 7"), 3, "dependabot.yml must set default-days: 7 on all 3 ecosystems")

        nf_main_text = (ROOT_DIR / "packaging" / "nf-core" / "modules" / "vardictrust" / "main.nf").read_text(encoding="utf-8")
        self.assertIn(
            'def input = somatic ? \'-b "\' + bam_input[0] + \'|\' + bam_input[1] + \'"\' : "-b ${bam_input}"',
            nf_main_text,
            "main.nf must double-quote somatic -b tumor|normal so | is never parsed as a shell pipe",
        )
        self.assertIn('-N "${meta.id}"', nf_main_text)
        self.assertNotIn("def shQuote =", nf_main_text, "main.nf must follow standard nf-core module idiom without custom Groovy closures")

        cfg_text = (ROOT_DIR / "harness" / "configure_github_repo.py").read_text(encoding="utf-8")
        self.assertIn('"Accept": "application/vnd.github+json"', cfg_text)
        self.assertNotIn("application/vnd.github.v3+json", cfg_text)
        self.assertIn('"has_wiki": True', cfg_text)
        self.assertIn('"allow_merge_commit": False', cfg_text)
        self.assertIn('"delete_branch_on_merge": True', cfg_text)
        self.assertIn('"required_linear_history": True', cfg_text)
        self.assertIn('"secret_scanning_non_provider_patterns"', cfg_text)
        self.assertIn('"secret_scanning_validity_checks"', cfg_text)
        for ctx in (
            '"Build, Unit, Integration, Parity & Coverage Gate"',
            '"Dual-Engine Differential Parity (Tier-2)"',
            '"Security & Supply-Chain Scan"',
            '"MkDocs Strict Build Gate"',
            '"plan"',
        ):
            self.assertIn(ctx, cfg_text, f"configure_github_repo.py missing required check context {ctx}")

    def test_github_workflows_valid_yaml_and_pinned_shas(self):
        workflows_dir = ROOT_DIR / ".github" / "workflows"
        wf_files = sorted(workflows_dir.glob("*.yml"))
        self.assertGreaterEqual(len(wf_files), 6, "Expected at least 6 GitHub workflow YAML files")

        action_files = sorted((ROOT_DIR / ".github" / "actions").glob("*/action.yml"))
        self.assertGreaterEqual(len(action_files), 1, "Expected at least 1 local composite action YAML file")

        diff_text = (workflows_dir / "differential.yml").read_text(encoding="utf-8")
        self.assertIn("uses: ./.github/actions/setup-conda-env", diff_text)

        try:
            yaml_mod = importlib.import_module("yaml")
        except ImportError:
            yaml_mod = None

        uses_re = re.compile(r"^\s*(?:-\s*)?uses:\s*([^\s#]+)", re.MULTILINE)
        sha_re = re.compile(r"@[0-9a-f]{40}$")

        for wf_path in wf_files + action_files:
            wf_text = wf_path.read_text(encoding="utf-8")
            if yaml_mod is not None:
                parsed = yaml_mod.safe_load(wf_text)
                self.assertIsInstance(parsed, dict, f"{wf_path.name} did not parse to a YAML mapping")
            for match in uses_re.finditer(wf_text):
                action_ref = match.group(1).strip("\"'")
                if action_ref.startswith("./"):
                    target_action = ROOT_DIR / action_ref.removeprefix("./") / "action.yml"
                    self.assertTrue(target_action.is_file(), f"Missing local composite action: {target_action}")
                    continue
                self.assertRegex(
                    action_ref,
                    sha_re,
                    f"{wf_path.name}: action reference {action_ref!r} must be pinned to a 40-char lowercase hex commit SHA",
                )

    def test_cargo_dist_multi_stage_lies_tag_checkout_simulation(self):
        """Simulate release.yml 4-stage cargo-dist (--artifacts=lies) across plan and unified tag checkout."""
        dist_bin = shutil.which("dist")
        if not dist_bin:
            candidate = pathlib.Path.home() / ".cargo" / "bin" / "dist"
            if candidate.is_file() and os.access(candidate, os.X_OK):
                dist_bin = str(candidate)
        if not dist_bin:
            self.skipTest("cargo-dist ('dist') binary not installed on host; covered by release.yml PR job")

        cargo_cfg = toml_lite.load_path(ROOT_DIR / "Cargo.toml")
        ws_ver = ((cargo_cfg.get("workspace") or {}).get("package") or {}).get("version")
        self.assertTrue(isinstance(ws_ver, str) and ws_ver, "Could not parse workspace version from Cargo.toml")
        tag_arg = f"--tag=v{ws_ver}"

        env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
        cargo_bin_dir = str(pathlib.Path.home() / ".cargo" / "bin")
        env["PATH"] = f"{cargo_bin_dir}:{env.get('PATH', '')}"

        dist_paths = (
            "Cargo.toml",
            "Cargo.lock",
            "CHANGELOG.md",
            "dist-workspace.toml",
            "README.md",
            "LICENSE",
            "crates",
        )

        def init_isolated_git_repo(dest: pathlib.Path) -> str:
            subprocess.run(["git", "init", "-q"], cwd=str(dest), env=env, check=True)
            subprocess.run(["git", "config", "user.email", "ci@vardictrust.local"], cwd=str(dest), env=env, check=True)
            subprocess.run(["git", "config", "user.name", "VarDictRust CI"], cwd=str(dest), env=env, check=True)
            subprocess.run(["git", "add", "-A"], cwd=str(dest), env=env, check=True)
            subprocess.run(["git", "commit", "-qm", "dist-e2e"], cwd=str(dest), env=env, check=True)
            return subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=str(dest),
                env=env,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()

        def strip_custom_runners(toml_text: str) -> str:
            return re.sub(
                r"\n*\[dist\.github-custom-runners\]\n(?:(?!\[)[^\n]*\n?)*",
                "\n",
                toml_text,
            )

        def copy_working_tree_subset(dest: pathlib.Path, legacy_tag_fixture: bool = False) -> str:
            dest.mkdir(parents=True, exist_ok=True)
            for rel in dist_paths:
                src = ROOT_DIR / rel
                if not src.exists():
                    continue
                dst = dest / rel
                if src.is_dir():
                    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("target", "*.profraw"))
                else:
                    shutil.copy2(src, dst)
            if legacy_tag_fixture:
                dw = dest / "dist-workspace.toml"
                dw.write_text(strip_custom_runners(dw.read_text(encoding="utf-8")), encoding="utf-8")
            return init_isolated_git_repo(dest)

        tag_ref = f"refs/tags/v{ws_ver}"
        tag_check = subprocess.run(
            ["git", "-C", str(ROOT_DIR), "rev-parse", "--verify", f"{tag_ref}^{{commit}}"],
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            plan_wt = tmp_path / "plan-wt"
            build_wt = tmp_path / "build-wt"
            plan_sha = copy_working_tree_subset(plan_wt, legacy_tag_fixture=False)
            used_tag_archive = False
            if tag_check.returncode == 0:
                tag_dw = subprocess.run(
                    ["git", "-C", str(ROOT_DIR), "show", f"{tag_check.stdout.strip()}:dist-workspace.toml"],
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if tag_dw.returncode == 0 and "[dist.github-custom-runners]" not in tag_dw.stdout:
                    build_wt.mkdir(parents=True, exist_ok=True)
                    archive_bytes = subprocess.run(
                        ["git", "-C", str(ROOT_DIR), "archive", "--format=tar", tag_check.stdout.strip(), *dist_paths],
                        env=env,
                        check=True,
                        capture_output=True,
                    ).stdout
                    subprocess.run(["tar", "-xf", "-", "-C", str(build_wt)], input=archive_bytes, env=env, check=True)
                    build_sha = init_isolated_git_repo(build_wt)
                    used_tag_archive = True
            if not used_tag_archive:
                build_sha = copy_working_tree_subset(build_wt, legacy_tag_fixture=True)

            self.assertNotEqual(plan_sha, build_sha, "plan_wt and build_wt must be distinct git commits")
            self.assertIn(
                "[dist.github-custom-runners]",
                (plan_wt / "dist-workspace.toml").read_text(encoding="utf-8"),
            )
            self.assertNotIn(
                "[dist.github-custom-runners]",
                (build_wt / "dist-workspace.toml").read_text(encoding="utf-8"),
                "build_wt must exercise cross-revision compatibility against a prior tag without custom runners",
            )

            # Stage 1: plan (dist host --steps=create) on current working-tree config
            res_plan = subprocess.run(
                [dist_bin, "host", "--steps=create", tag_arg, "--output-format=json"],
                cwd=str(plan_wt),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                res_plan.returncode,
                0,
                f"dist host --steps=create failed:\n{res_plan.stderr}",
            )
            plan_manifest = json.loads(res_plan.stdout)
            include = (
                ((plan_manifest.get("ci") or {}).get("github") or {}).get("artifacts_matrix") or {}
            ).get("include") or []
            self.assertEqual(len(include), 5, f"Expected 5 target matrix entries, got {len(include)}")
            runners_by_target = {tuple(r.get("targets") or ()): r.get("runner") for r in include}
            self.assertEqual(runners_by_target.get(("aarch64-apple-darwin",)), "macos-15")
            retired = {"macos-11", "macos-12", "macos-13", "macos-14", "ubuntu-18.04", "ubuntu-20.04"}
            for targets, runner in runners_by_target.items():
                self.assertNotIn(runner, retired, f"Target {targets} uses retired runner {runner}")

            distrib_dir = build_wt / "target" / "distrib"
            distrib_dir.mkdir(parents=True, exist_ok=True)
            (distrib_dir / "plan-dist-manifest.json").write_text(res_plan.stdout, encoding="utf-8")

            # Stage 2: build-local-artifacts (--artifacts=lies) for all 5 targets on downstream tag tree
            for row in include:
                targets_slug = "-".join(row["targets"])
                dist_args = row["dist_args"].split()
                dist_args = ["--artifacts=lies" if a == "--artifacts=local" else a for a in dist_args]
                res_local = subprocess.run(
                    [dist_bin, "build", tag_arg, "--output-format=json", *dist_args],
                    cwd=str(build_wt),
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(
                    res_local.returncode,
                    0,
                    f"dist build --artifacts=lies failed for {targets_slug}:\n{res_local.stderr}",
                )
                (distrib_dir / f"{targets_slug}-dist-manifest.json").write_text(
                    res_local.stdout, encoding="utf-8"
                )

            # Stage 3: build-global-artifacts on downstream tag tree
            res_global = subprocess.run(
                [dist_bin, "build", tag_arg, "--output-format=json", "--artifacts=global"],
                cwd=str(build_wt),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                res_global.returncode,
                0,
                f"dist build --artifacts=global failed:\n{res_global.stderr}",
            )
            (distrib_dir / "global-dist-manifest.json").write_text(res_global.stdout, encoding="utf-8")

            # Stage 4: host (--steps=upload --steps=release) on downstream tag tree
            res_host = subprocess.run(
                [dist_bin, "host", tag_arg, "--steps=upload", "--steps=release", "--output-format=json"],
                cwd=str(build_wt),
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                res_host.returncode,
                0,
                f"dist host --steps=upload --steps=release failed:\n{res_host.stderr}",
            )
            host_manifest = json.loads(res_host.stdout)
            self.assertTrue(host_manifest.get("announcement_title"), "Missing announcement_title")
            self.assertTrue(host_manifest.get("announcement_github_body"), "Missing announcement_github_body")

            for expected_file in (
                "vardictrust-x86_64-unknown-linux-gnu.tar.xz",
                "vardictrust-x86_64-unknown-linux-musl.tar.xz",
                "vardictrust-aarch64-unknown-linux-musl.tar.xz",
                "vardictrust-aarch64-apple-darwin.tar.xz",
                "vardictrust-x86_64-apple-darwin.tar.xz",
                "source.tar.gz",
                "vardictrust-installer.sh",
                "sha256.sum",
            ):
                self.assertTrue(
                    (distrib_dir / expected_file).is_file(),
                    f"Expected staged artifact {expected_file} missing in {distrib_dir}",
                )

    def test_release_slsa_predicate_and_bundle_verification_logic(self):
        """Execute release.yml's inline SLSA predicate generator & bundle verifier for push and dispatch."""
        rel_text = (ROOT_DIR / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
        gen_m = re.search(
            r"- name: Generate SLSA Provenance Predicate\n.*?python3 - <<'PY'\n(.*?)\n\s+PY",
            rel_text,
            re.DOTALL,
        )
        ver_m = re.search(
            r"- name: Verify Signed Attestation Bundle\n.*?python3 - <<'PY'\n(.*?)\n\s+PY",
            rel_text,
            re.DOTALL,
        )
        self.assertIsNotNone(gen_m, "Could not extract Generate SLSA Provenance Predicate script from release.yml")
        self.assertIsNotNone(ver_m, "Could not extract Verify Signed Attestation Bundle script from release.yml")

        gen_code = textwrap.dedent(gen_m.group(1))
        ver_code = textwrap.dedent(ver_m.group(1))

        active_claims: dict[str, str] = {}

        class MockOIDCHandler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                if "audience=nobody" not in self.path or self.headers.get("Authorization") != "Bearer test-req-tok":
                    self.send_response(401)
                    self.end_headers()
                    return
                hdr = base64.urlsafe_b64encode(b'{"alg":"RS256","typ":"JWT"}').decode("ascii").rstrip("=")
                pld = base64.urlsafe_b64encode(json.dumps(active_claims).encode("utf-8")).decode("ascii").rstrip("=")
                jwt_tok = f"{hdr}.{pld}.sig"
                body = json.dumps({"value": jwt_tok}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, format, *args):
                pass

        server = http.server.HTTPServer(("127.0.0.1", 0), MockOIDCHandler)
        port = server.server_address[1]
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            tag_sha = "3637d507af4e140e3ed5b6a883b9792f4aed1ebe"
            wf_sha = "9999999999999999999999999999999999999999"
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = pathlib.Path(tmpdir)
                pred_file = tmp_path / "slsa-predicate.json"
                bundle_file = tmp_path / "attestation.json"

                # Case 1: workflow_dispatch on refs/heads/main targeting older refs/tags/v1.0.8
                active_claims.clear()
                active_claims.update(
                    {
                        "iss": "https://token.actions.githubusercontent.com",
                        "ref": "refs/heads/main",
                        "sha": wf_sha,
                        "repository": "palakpsheth/VarDictRust",
                        "event_name": "workflow_dispatch",
                        "job_workflow_ref": "palakpsheth/VarDictRust/.github/workflows/release.yml@refs/heads/main",
                        "workflow_ref": "palakpsheth/VarDictRust/.github/workflows/release.yml@refs/heads/main",
                        "repository_id": "1251726015",
                        "repository_owner_id": "13122115",
                        "runner_environment": "github-hosted",
                        "run_id": "44000000000",
                        "run_attempt": "1",
                    }
                )
                base_env = os.environ.copy()
                base_env.update(
                    {
                        "ACTIONS_ID_TOKEN_REQUEST_URL": f"http://127.0.0.1:{port}/token?api-version=2.0",
                        "ACTIONS_ID_TOKEN_REQUEST_TOKEN": "test-req-tok",
                        "GITHUB_SERVER_URL": "https://github.com",
                        "GITHUB_REPOSITORY": "palakpsheth/VarDictRust",
                        "GITHUB_REF": "refs/heads/main",
                        "GITHUB_SHA": wf_sha,
                        "EXPECTED_TAG": "v1.0.8",
                        "EXPECTED_REF": "refs/tags/v1.0.8",
                        "EXPECTED_COMMIT": tag_sha,
                        "PREDICATE_OUT": str(pred_file),
                        "BUNDLE_PATH": str(bundle_file),
                    }
                )
                res_gen = subprocess.run(
                    [sys.executable, "-c", gen_code],
                    env=base_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(res_gen.returncode, 0, f"Generator failed:\n{res_gen.stderr}")
                pred = json.loads(pred_file.read_text(encoding="utf-8"))
                self.assertEqual(
                    pred["buildDefinition"]["externalParameters"],
                    {
                        "workflow": {
                            "ref": "refs/heads/main",
                            "repository": "https://github.com/palakpsheth/VarDictRust",
                            "path": ".github/workflows/release.yml",
                        },
                        "inputs": {
                            "tag": "v1.0.8",
                        },
                    },
                )
                deps = pred["buildDefinition"]["resolvedDependencies"]
                self.assertEqual(len(deps), 2)
                self.assertEqual(
                    deps[0],
                    {
                        "uri": "git+https://github.com/palakpsheth/VarDictRust@refs/heads/main",
                        "digest": {"gitCommit": wf_sha},
                    },
                )
                self.assertEqual(
                    deps[1],
                    {
                        "uri": "git+https://github.com/palakpsheth/VarDictRust@refs/tags/v1.0.8",
                        "digest": {"gitCommit": tag_sha},
                    },
                )

                # Verify bundle verifier passes on valid signed statement and rejects wrong commit or swapped order
                stmt = {
                    "_type": "https://in-toto.io/Statement/v1",
                    "predicateType": "https://slsa.dev/provenance/v1",
                    "subject": [{"name": "vardictrust.tar.xz", "digest": {"sha256": "a" * 64}}],
                    "predicate": pred,
                }
                bundle_file.write_text(
                    json.dumps(
                        {
                            "dsseEnvelope": {
                                "payloadType": "application/vnd.in-toto+json",
                                "payload": base64.b64encode(json.dumps(stmt).encode("utf-8")).decode("ascii"),
                            }
                        }
                    )
                    + "\n",
                    encoding="utf-8",
                )
                res_ver = subprocess.run(
                    [sys.executable, "-c", ver_code],
                    env=base_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(res_ver.returncode, 0, f"Verifier failed:\n{res_ver.stderr}")

                bad_env = dict(base_env, EXPECTED_COMMIT="0" * 40)
                res_bad = subprocess.run(
                    [sys.executable, "-c", ver_code],
                    env=bad_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertNotEqual(res_bad.returncode, 0, "Verifier must reject mismatched release-commit")

                swapped_pred = json.loads(json.dumps(pred))
                swapped_pred["buildDefinition"]["resolvedDependencies"] = list(reversed(deps))
                swapped_stmt = dict(stmt, predicate=swapped_pred)
                bundle_file.write_text(
                    json.dumps(
                        {
                            "dsseEnvelope": {
                                "payloadType": "application/vnd.in-toto+json",
                                "payload": base64.b64encode(json.dumps(swapped_stmt).encode("utf-8")).decode("ascii"),
                            }
                        }
                    )
                    + "\n",
                    encoding="utf-8",
                )
                res_swapped = subprocess.run(
                    [sys.executable, "-c", ver_code],
                    env=base_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertNotEqual(
                    res_swapped.returncode,
                    0,
                    "Verifier must reject resolvedDependencies[0] != workflow ref/sha",
                )

                # Case 2: direct tag push (claims.ref == EXPECTED_REF and claims.sha == EXPECTED_COMMIT)
                active_claims.update(
                    {
                        "ref": "refs/tags/v1.0.8",
                        "sha": tag_sha,
                        "event_name": "push",
                        "job_workflow_ref": "palakpsheth/VarDictRust/.github/workflows/release.yml@refs/tags/v1.0.8",
                        "workflow_ref": "palakpsheth/VarDictRust/.github/workflows/release.yml@refs/tags/v1.0.8",
                    }
                )
                push_env = dict(base_env, GITHUB_REF="refs/tags/v1.0.8", GITHUB_SHA=tag_sha)
                res_push = subprocess.run(
                    [sys.executable, "-c", gen_code],
                    env=push_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(res_push.returncode, 0, f"Push generator failed:\n{res_push.stderr}")
                pred_push = json.loads(pred_file.read_text(encoding="utf-8"))
                self.assertNotIn("inputs", pred_push["buildDefinition"]["externalParameters"])
                self.assertEqual(
                    pred_push["buildDefinition"]["resolvedDependencies"],
                    [
                        {
                            "uri": "git+https://github.com/palakpsheth/VarDictRust@refs/tags/v1.0.8",
                            "digest": {"gitCommit": tag_sha},
                        }
                    ],
                )
        finally:
            server.shutdown()
            server.server_close()

    def test_release_workflow_sha256_bootstrap_and_harden_runner(self):
        """Verify SEC-05: 4-arch cargo-dist v0.33.0 SHA-256 verification and 5-job harden-runner."""
        rel_text = (ROOT_DIR / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")

        for forbidden in (
            "| sh",
            "cargo-dist-installer.sh",
            "sh.rustup.rs",
            "${{ matrix.install_dist.run }}",
            "dtolnay/rust-toolchain",
        ):
            self.assertNotIn(
                forbidden,
                rel_text,
                f"release.yml must not contain unpinned bootstrap construct {forbidden!r}",
            )
        self.assertIn(
            "unverified container build(s)",
            rel_text,
            "release.yml plan guard must reject container builds lacking preinstalled Rust toolchains",
        )

        # Harden-runner pinned and placed as first step in all 5 jobs
        harden_ref = "step-security/harden-runner@1d88384fd173af6ecb012dfd0264a306b711241a"
        expected_jobs = (
            "plan",
            "build-local-artifacts",
            "build-global-artifacts",
            "host",
            "publish-crates",
        )
        self.assertEqual(
            rel_text.count(f"uses: {harden_ref}"),
            len(expected_jobs),
            "harden-runner must appear once in each of the 5 release.yml jobs",
        )
        self.assertEqual(
            rel_text.count("egress-policy: audit"),
            len(expected_jobs),
            "Every harden-runner step in release.yml must set egress-policy: audit",
        )
        for job_name in expected_jobs:
            job_m = re.search(
                rf"^  {re.escape(job_name)}:\n((?:  (?![A-Za-z0-9_-]+:).*\n|\n)+?)(?=^  [A-Za-z0-9_-]+:|\Z)",
                rel_text,
                re.MULTILINE,
            )
            self.assertIsNotNone(job_m, f"Missing job {job_name!r} in release.yml")
            job_body = job_m.group(1)
            steps_m = re.search(r"^    steps:\n((?:      .*\n|\n)+)", job_body, re.MULTILINE)
            self.assertIsNotNone(steps_m, f"Job {job_name!r} missing steps block")
            steps_body = steps_m.group(1)
            first_step_lines = []
            for idx, ln in enumerate(steps_body.splitlines()):
                if idx > 0 and re.match(r"^      - ", ln):
                    break
                first_step_lines.append(ln)
            first_step = "\n".join(first_step_lines)
            self.assertIn(
                f"uses: {harden_ref}",
                first_step,
                f"First step in job {job_name!r} must be {harden_ref}",
            )
            self.assertIn(
                "egress-policy: audit",
                first_step,
                f"First step in job {job_name!r} must configure egress-policy: audit",
            )
            if job_name == "build-local-artifacts":
                self.assertIn(
                    "if: runner.os == 'Linux' && runner.arch == 'X64'",
                    first_step,
                    "Cross-OS job build-local-artifacts must gate harden-runner to Linux X64",
                )

        # Cache ~/.cache/cargo-dist and verify all 4 v0.33.0 archive SHA-256 digests in plan & build-local-artifacts
        self.assertEqual(
            rel_text.count("path: ~/.cache/cargo-dist"),
            2,
            "Both plan and build-local-artifacts must cache ~/.cache/cargo-dist",
        )
        self.assertEqual(
            rel_text.count("key: cargo-dist-archive-sha256-v0.33.0-${{ runner.os }}-${{ runner.arch }}"),
            2,
            "Both plan and build-local-artifacts must use cargo-dist-archive-sha256-v0.33.0 cache key",
        )

        expected_archives = {
            "cargo-dist-x86_64-unknown-linux-gnu.tar.xz": "4b3f0a5f0ebbdb798f6db649d01b32ba1518376b6f7a0502b7d92b75cc2c8293",
            "cargo-dist-aarch64-unknown-linux-gnu.tar.xz": "9c554ab21a58ad46eb9b6710f89633ccb26bc2577cf4b0c2d18a5b826a31913e",
            "cargo-dist-x86_64-apple-darwin.tar.xz": "6a49bfb61bd86770d79c27f3d2b40c6b2e71cde940d3d31a6ccaaffc124d7a29",
            "cargo-dist-aarch64-apple-darwin.tar.xz": "7b3cbe25511de01d74c0f5fcb7909edabd379bea9cfa284d93af5a3cdfa3247c",
        }
        for archive_name, digest in expected_archives.items():
            self.assertEqual(
                rel_text.count(archive_name),
                2,
                f"Archive {archive_name} must be pinned in both plan and build-local-artifacts",
            )
            self.assertEqual(
                rel_text.count(digest),
                2,
                f"SHA-256 {digest} for {archive_name} must be pinned in both plan and build-local-artifacts",
            )

        install_blocks = re.findall(
            r"      - name: Install dist\n        shell: bash\n        run: \|\n((?:          .*\n|\n)+)",
            rel_text,
        )
        self.assertEqual(len(install_blocks), 2, "Expected 2 'Install dist' bash blocks in release.yml")
        for block in install_blocks:
            script = textwrap.dedent(block)
            pos_dl_end = script.find("\nfi\n")
            pos_sha = script.find("hashlib.sha256")
            pos_tar = script.find("--strip-components=1")
            pos_path = script.find('"$GITHUB_PATH"')
            pos_ver = script.find('"$HOME/.cargo/bin/dist" --version')
            self.assertNotEqual(pos_dl_end, -1, "Install dist must guard download with if [ ! -f ... ]; then ... fi")
            self.assertGreater(pos_sha, pos_dl_end, "SHA-256 check must execute outside download block (covers cache hits)")
            self.assertGreater(pos_tar, pos_sha, "tar extraction (--strip-components=1) must occur AFTER SHA-256 verification")
            self.assertGreater(pos_path, pos_tar, "$GITHUB_PATH update must occur after extraction")
            self.assertGreater(pos_ver, pos_tar, "dist --version must occur AFTER SHA-256 verification and extraction")

        # Behavioral execution test on cache hit: tampered cached archive is rejected before tar runs
        script_to_test = textwrap.dedent(install_blocks[0])
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            fake_home = tmp_path / "home"
            cache_dir = fake_home / ".cache" / "cargo-dist"
            cache_dir.mkdir(parents=True)
            gh_path = tmp_path / "github_path"
            gh_path.write_text("", encoding="utf-8")

            bad_archive = cache_dir / "cargo-dist-x86_64-unknown-linux-gnu.tar.xz"
            bad_archive.write_bytes(b"tampered-bytes")
            fake_bin = tmp_path / "fakebin"
            fake_bin.mkdir()
            uname_sh = fake_bin / "uname"
            uname_sh.write_text(
                '#!/usr/bin/env bash\nif [ "${1:-}" = "-s" ]; then echo Linux; else echo x86_64; fi\n',
                encoding="utf-8",
            )
            uname_sh.chmod(0o755)
            test_env = dict(
                os.environ,
                HOME=str(fake_home),
                GITHUB_PATH=str(gh_path),
                PATH=f"{fake_bin}:{os.environ.get('PATH', '')}",
            )
            res_tampered = subprocess.run(
                ["bash", "-c", script_to_test],
                env=test_env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(res_tampered.returncode, 0, "Tampered cached cargo-dist archive must fail closed")
            self.assertIn("SHA-256 mismatch", res_tampered.stderr)
            self.assertFalse(
                (fake_home / ".cargo" / "bin" / "dist").exists(),
                "dist binary must never be extracted when SHA-256 verification fails",
            )

            nested_dir = tmp_path / "pkg" / "cargo-dist-x86_64-unknown-linux-gnu"
            nested_dir.mkdir(parents=True)
            mock_dist = nested_dir / "dist"
            mock_dist.write_text('#!/usr/bin/env bash\necho "cargo-dist 0.33.0"\n', encoding="utf-8")
            mock_dist.chmod(0o755)
            with tarfile.open(bad_archive, "w:xz") as tf:
                tf.add(mock_dist, arcname="cargo-dist-x86_64-unknown-linux-gnu/dist")
            synth_sha = hashlib.sha256(bad_archive.read_bytes()).hexdigest()
            patched_script = script_to_test.replace(
                "4b3f0a5f0ebbdb798f6db649d01b32ba1518376b6f7a0502b7d92b75cc2c8293",
                synth_sha,
            )
            res_valid = subprocess.run(
                ["bash", "-c", patched_script],
                env=test_env,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(res_valid.returncode, 0, f"Valid cached archive failed:\n{res_valid.stderr}")
            self.assertIn("cargo-dist 0.33.0", res_valid.stdout)
            self.assertIn(str(fake_home / ".cargo" / "bin"), gh_path.read_text(encoding="utf-8"))

    def test_release_publish_crates_idempotency_and_fast_fail_retry(self):
        """Verify AR-21: cargo publish steps treat already-uploaded as idempotent and only retry propagation lag."""
        rel_text = (ROOT_DIR / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
        core_m = re.search(
            r"      - name: Publish vardictrust-core to Crates\.io\n((?:        .*\n|\n)+?)(?=      - name:|^  [A-Za-z0-9_-]+:|\Z)",
            rel_text,
            re.MULTILINE,
        )
        bin_m = re.search(
            r"      - name: Publish vardictrust to Crates\.io\n((?:        .*\n|\n)+?)(?=      - name:|^  [A-Za-z0-9_-]+:|\Z)",
            rel_text,
            re.MULTILINE,
        )
        self.assertIsNotNone(core_m, "Missing Publish vardictrust-core to Crates.io step")
        self.assertIsNotNone(bin_m, "Missing Publish vardictrust to Crates.io step")
        for label, step_block in (("core", core_m.group(1)), ("bin", bin_m.group(1))):
            self.assertIn("shell: bash", step_block, f"{label} publish step missing shell: bash")
            self.assertIn("CARGO_REGISTRY_TOKEN: ${{ secrets.CARGO_REGISTRY_TOKEN }}", step_block)
            self.assertIn("set -euo pipefail", step_block, f"{label} publish step missing set -euo pipefail")
            self.assertIn("is already uploaded", step_block, f"{label} publish step missing already-uploaded guard")

        core_run_m = re.search(r"        run: \|\n((?:          .*\n|\n)+)", core_m.group(1))
        bin_run_m = re.search(r"        run: \|\n((?:          .*\n|\n)+)", bin_m.group(1))
        self.assertIsNotNone(core_run_m, "Publish vardictrust-core step must use multiline run: | script")
        self.assertIsNotNone(bin_run_m, "Publish vardictrust step must use multiline run: | script")
        core_script = textwrap.dedent(core_run_m.group(1))
        bin_script = textwrap.dedent(bin_run_m.group(1))

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            calls_file = tmp_path / "calls.txt"
            mode_file = tmp_path / "mode.txt"

            sleep_sh = fake_bin / "sleep"
            sleep_sh.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            sleep_sh.chmod(0o755)

            cargo_sh = fake_bin / "cargo"
            cargo_sh.write_text(
                f"""#!/usr/bin/env bash
count=$(( $(cat "{calls_file}" 2>/dev/null || echo 0) + 1 ))
echo "$count" > "{calls_file}"
mode="$(cat "{mode_file}")"
pkg="${{3:-vardictrust}}"
case "$mode" in
  already)
    echo "error: failed to publish to registry at https://crates.io: crate version \\`1.0.9\\` is already uploaded" >&2
    exit 101
    ;;
  already_index)
    echo "error: crate ${{pkg}}@1.0.9 already exists on crates.io index" >&2
    exit 101
    ;;
  lag_then_ok)
    if [ "$count" -eq 1 ]; then
      echo "error: no matching package named \\`vardictrust-core\\` found" >&2
      exit 101
    fi
    exit 0
    ;;
  auth_fail)
    echo "error: failed to authenticate to registry: invalid token" >&2
    exit 101
    ;;
esac
""",
                encoding="utf-8",
            )
            cargo_sh.chmod(0o755)
            test_env = dict(
                os.environ,
                PATH=f"{fake_bin}:{os.environ.get('PATH', '')}",
                RELEASE_TAG="v1.0.9",
                CARGO_REGISTRY_TOKEN="dummy",
            )

            # 1a. Already uploaded (upload error) -> exits 0 on attempt 1
            calls_file.write_text("0", encoding="utf-8")
            mode_file.write_text("already", encoding="utf-8")
            r1 = subprocess.run(["bash", "-c", bin_script], env=test_env, capture_output=True, text=True)
            self.assertEqual(r1.returncode, 0, r1.stderr)
            self.assertEqual(int(calls_file.read_text(encoding="utf-8").strip()), 1)

            # 1b. Sparse-index pre-check ("crate <pkg>@1.0.9 already exists on crates.io index") -> exits 0 on both steps
            for step_label, script in (("core", core_script), ("bin", bin_script)):
                calls_file.write_text("0", encoding="utf-8")
                mode_file.write_text("already_index", encoding="utf-8")
                r_idx = subprocess.run(["bash", "-c", script], env=test_env, capture_output=True, text=True)
                self.assertEqual(r_idx.returncode, 0, f"{step_label} failed on already_index:\n{r_idx.stderr}")
                self.assertEqual(int(calls_file.read_text(encoding="utf-8").strip()), 1)

            # 2. Index propagation lag -> retries and succeeds on attempt 2
            calls_file.write_text("0", encoding="utf-8")
            mode_file.write_text("lag_then_ok", encoding="utf-8")
            r2 = subprocess.run(["bash", "-c", bin_script], env=test_env, capture_output=True, text=True)
            self.assertEqual(r2.returncode, 0, r2.stderr)
            self.assertEqual(int(calls_file.read_text(encoding="utf-8").strip()), 2)

            # 3. Non-propagation failure (auth/compile) -> fails fast on attempt 1 without 30 retries
            calls_file.write_text("0", encoding="utf-8")
            mode_file.write_text("auth_fail", encoding="utf-8")
            r3 = subprocess.run(["bash", "-c", bin_script], env=test_env, capture_output=True, text=True)
            self.assertNotEqual(r3.returncode, 0)
            self.assertEqual(int(calls_file.read_text(encoding="utf-8").strip()), 1)

    def test_configure_github_repo_tag_ruleset_and_lockout_rollback(self):
        """Verify SEC-12b: configure_github_repo.py --dry-run, refs/tags/v* ruleset 24867690, and lockout rollback."""
        cfg_script = ROOT_DIR / "harness" / "configure_github_repo.py"
        cfg_text = cfg_script.read_text(encoding="utf-8")
        self.assertIn("def configure_tag_ruleset", cfg_text)
        self.assertIn("24867690", cfg_text)
        self.assertIn('"refs/tags/v*"', cfg_text)
        self.assertNotIn('"~refs/tags/v*"', cfg_text)
        self.assertIn("current_user_can_bypass", cfg_text)
        self.assertIn("urlopen(req, timeout=30)", cfg_text)
        self.assertIn("--dry-run", cfg_text)
        self.assertIn("--apply", cfg_text)

        clean_env = {"PATH": os.environ.get("PATH", "")}
        res_dry = subprocess.run(
            [sys.executable, str(cfg_script), "--dry-run"],
            env=clean_env,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(res_dry.returncode, 0, f"--dry-run failed:\n{res_dry.stderr}")
        self.assertIn("DRY-RUN", res_dry.stdout.upper())

        spec = importlib.util.spec_from_file_location("configure_github_repo", cfg_script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        self.assertTrue(hasattr(mod, "REPO_DESCRIPTION"), "configure_github_repo.py must define REPO_DESCRIPTION")
        self.assertLessEqual(len(mod.REPO_DESCRIPTION), 80)
        self.assertFalse(mod.REPO_DESCRIPTION.endswith("."))
        self.assertIn("VarDictJava", mod.REPO_DESCRIPTION)
        self.assertNotIn("--compat", mod.REPO_DESCRIPTION)
        site_desc = toml_lite.load_path(ROOT_DIR / "zensical.toml")["project"]["site_description"]
        self.assertLessEqual(len(site_desc), 80)
        self.assertFalse(site_desc.endswith("."))
        self.assertIn("VarDictJava", site_desc)
        self.assertNotIn("--compat", site_desc)
        self.assertTrue(hasattr(mod, "REPO_TOPICS"), "configure_github_repo.py must define REPO_TOPICS")
        self.assertEqual(len(mod.REPO_TOPICS), 14)
        self.assertEqual(mod.configure_metadata(dry_run=True)["description"], mod.REPO_DESCRIPTION)
        self.assertEqual(mod.configure_topics(dry_run=True)["names"], mod.REPO_TOPICS)

        # Case A: current_user_can_bypass == "always" -> keeps creation rule + RepositoryRole actor 5 (no rollback)
        calls_ok = []

        def mock_ok(endpoint, method="GET", data=None):
            calls_ok.append((endpoint, method, data))
            if method == "PUT":
                return {"id": 24867690}
            return {"id": 24867690, "current_user_can_bypass": "always"}

        mod.api_request = mock_ok
        mod.configure_tag_ruleset(dry_run=False)
        self.assertEqual(len(calls_ok), 2, f"Expected PUT + GET when bypass=='always', got {len(calls_ok)}")
        self.assertEqual(calls_ok[0][0], "/rulesets/24867690")
        self.assertEqual(calls_ok[0][1], "PUT")
        self.assertEqual(calls_ok[0][2]["conditions"]["ref_name"]["include"], ["refs/tags/v*"])
        self.assertIn({"type": "creation"}, calls_ok[0][2]["rules"])
        self.assertEqual(
            calls_ok[0][2]["bypass_actors"],
            [{"actor_id": 5, "actor_type": "RepositoryRole", "bypass_mode": "always"}],
        )
        self.assertEqual(calls_ok[1], ("/rulesets/24867690", "GET", None))

        # Case B: current_user_can_bypass != "always" -> triggers automatic rollback PUT without creation rule
        calls_deny = []

        def mock_deny(endpoint, method="GET", data=None):
            calls_deny.append((endpoint, method, data))
            if method == "PUT":
                return {"id": 24867690}
            return {"id": 24867690, "current_user_can_bypass": "never"}

        mod.api_request = mock_deny
        mod.configure_tag_ruleset(dry_run=False)
        self.assertEqual(len(calls_deny), 3, f"Expected PUT + GET + rollback PUT, got {len(calls_deny)}")
        rollback_payload = calls_deny[2][2]
        self.assertEqual(calls_deny[2][1], "PUT")
        self.assertEqual(rollback_payload["conditions"]["ref_name"]["include"], ["refs/tags/v*"])
        self.assertNotIn({"type": "creation"}, rollback_payload["rules"])
        self.assertEqual(
            rollback_payload["rules"],
            [{"type": "deletion"}, {"type": "update"}, {"type": "non_fast_forward"}],
        )
        self.assertEqual(rollback_payload["bypass_actors"], [])

        # Case C: verification GET raises Exception -> treats can_bypass as None and executes safe rollback PUT
        calls_err = []

        def mock_err(endpoint, method="GET", data=None):
            calls_err.append((endpoint, method, data))
            if method == "GET":
                raise RuntimeError("transient GitHub 502 on verification GET")
            return {"id": 24867690}

        mod.api_request = mock_err
        mod.configure_tag_ruleset(dry_run=False)
        self.assertEqual(len(calls_err), 3, f"Expected PUT + failed GET + rollback PUT, got {len(calls_err)}")
        self.assertEqual(calls_err[2][1], "PUT")
        self.assertNotIn({"type": "creation"}, calls_err[2][2]["rules"])
        self.assertEqual(calls_err[2][2]["bypass_actors"], [])

        # Case D: strict PUT raises HTTPError(422) -> falls back to PUT without bypass_actors/creation; non-422 re-raises
        calls_422 = []

        def mock_422(endpoint, method="GET", data=None):
            calls_422.append((endpoint, method, data))
            if method == "PUT" and data and data.get("bypass_actors"):
                raise urllib.error.HTTPError(
                    endpoint, 422, "Unprocessable Entity", hdrs=None, fp=io.BytesIO(b'{"message":"invalid bypass"}')
                )
            return {"id": 24867690}

        mod.api_request = mock_422
        mod.configure_tag_ruleset(dry_run=False)
        self.assertEqual(len(calls_422), 2, f"Expected failed strict PUT + fallback PUT, got {len(calls_422)}")
        self.assertEqual(calls_422[1][1], "PUT")
        self.assertNotIn({"type": "creation"}, calls_422[1][2]["rules"])
        self.assertEqual(calls_422[1][2]["bypass_actors"], [])

        def mock_500(endpoint, method="GET", data=None):
            raise urllib.error.HTTPError(
                endpoint, 500, "Internal Server Error", hdrs=None, fp=io.BytesIO(b'{"message":"boom"}')
            )

        mod.api_request = mock_500
        with self.assertRaises(urllib.error.HTTPError) as ctx_500:
            mod.configure_tag_ruleset(dry_run=False)
        self.assertEqual(ctx_500.exception.code, 500)

    def test_all_six_workflows_declare_explicit_job_timeouts(self):
        """Verify AR-5: every job across all 6 .github/workflows/*.yml files declares an integer timeout-minutes."""
        workflows_dir = ROOT_DIR / ".github" / "workflows"
        expected_timeouts = {
            "ci.yml": {"build-and-test": 45},
            "differential.yml": {"differential": 15},
            "docs.yml": {"build-docs": 20},
            "release.yml": {
                "plan": 15,
                "build-local-artifacts": 45,
                "build-global-artifacts": 20,
                "host": 20,
                "publish-crates": 20,
            },
            "security.yml": {"security-scan": 30},
            "wiki.yml": {"sync": 10},
        }
        actual_wf_names = {p.name for p in workflows_dir.glob("*.yml")}
        self.assertEqual(
            actual_wf_names,
            set(expected_timeouts.keys()),
            f"Unexpected .github/workflows/*.yml file set: {sorted(actual_wf_names)}",
        )

        for wf_name, expected_jobs in sorted(expected_timeouts.items()):
            wf_text = (workflows_dir / wf_name).read_text(encoding="utf-8")
            jobs_m = re.search(r"^jobs:\s*\n((?:[ \t]+.*\n|#.*\n|\n)+)", wf_text, re.MULTILINE)
            self.assertIsNotNone(jobs_m, f"{wf_name} missing top-level jobs: block")
            jobs_body = jobs_m.group(1)
            job_matches = list(
                re.finditer(
                    r"^  ([A-Za-z0-9_-]+):\s*(?:#.*)?\n((?:  (?![A-Za-z0-9_-]+:).*\n|\n)+?)(?=^  [A-Za-z0-9_-]+:|\Z)",
                    jobs_body,
                    re.MULTILINE,
                )
            )
            self.assertTrue(job_matches, f"{wf_name} parsed zero jobs under jobs:")
            parsed_timeouts: dict[str, int] = {}
            for jm in job_matches:
                job_id = jm.group(1)
                job_body = jm.group(2)
                job_header = job_body.split("    steps:\n", 1)[0]
                tm_m = re.search(r"^    timeout-minutes:\s*([0-9]+)\s*(?:#.*)?$", job_header, re.MULTILINE)
                self.assertIsNotNone(
                    tm_m,
                    f"{wf_name} job {job_id!r} missing explicit integer job-level 'timeout-minutes:'",
                )
                timeout_val = int(tm_m.group(1))
                self.assertGreater(timeout_val, 0, f"{wf_name} job {job_id!r} timeout-minutes must be > 0")
                parsed_timeouts[job_id] = timeout_val
            self.assertEqual(
                parsed_timeouts,
                expected_jobs,
                f"{wf_name} job timeout-minutes mismatch: got {parsed_timeouts}, expected {expected_jobs}",
            )

    def test_paper_reproduce_all_two_tier_gate_and_portability(self):
        """Verify AR-10b/AR-14: paper/reproduce_all.sh has zero hardcoded user-home literals, handles detached HEAD, and parses flags in any order."""
        repro_path = ROOT_DIR / "paper" / "reproduce_all.sh"
        self.assertTrue(repro_path.is_file(), f"Missing {repro_path}")
        repro_text = repro_path.read_text(encoding="utf-8")
        forbidden_home = "/" + "home/" + "palak"
        self.assertNotIn(forbidden_home, repro_text, "paper/reproduce_all.sh must not hardcode user home path")
        self.assertIn("$HOME/.local/bin/uv", repro_text)
        self.assertIn(
            "git symbolic-ref -q --short HEAD 2>/dev/null || true",
            repro_text,
            "paper/reproduce_all.sh must guard git symbolic-ref for detached-HEAD checkouts",
        )

        skill_path = ROOT_DIR / ".agents" / "skills" / "vardictrust-release" / "SKILL.md"
        skill_text = skill_path.read_text(encoding="utf-8")
        self.assertIn("python3 paper/build_manuscript.py", skill_text)
        self.assertIn("20 version-bearing files", skill_text)
        self.assertIn("18 pre-tag files", skill_text)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_root = pathlib.Path(tmpdir) / "repo"
            (tmp_root / "paper").mkdir(parents=True)
            (tmp_root / "harness").mkdir(parents=True)
            (tmp_root / ".venv" / "bin").mkdir(parents=True)
            (tmp_root / "mock_bin").mkdir(parents=True)
            raw_cells = tmp_root / "benchmarks" / "data" / "gcp-study-v106"
            for idx in range(60):
                (raw_cells / f"cell_{idx:02d}").mkdir(parents=True)

            copied_repro = tmp_root / "paper" / "reproduce_all.sh"
            shutil.copy2(repro_path, copied_repro)
            copied_repro.chmod(0o755)

            inv_log = pathlib.Path(tmpdir) / "invocations.log"
            stub_py = tmp_root / ".venv" / "bin" / "python"
            stub_py.write_text(
                f"#!{sys.executable}\n"
                "import pathlib, sys\n"
                f"log = pathlib.Path({str(inv_log)!r})\n"
                "with log.open('a', encoding='utf-8') as fh:\n"
                "    fh.write('PY ' + ' '.join(sys.argv[1:]) + '\\n')\n",
                encoding="utf-8",
            )
            stub_py.chmod(0o755)

            stub_pytest = tmp_root / ".venv" / "bin" / "pytest"
            stub_pytest.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            stub_pytest.chmod(0o755)

            stub_gcloud = tmp_root / "mock_bin" / "gcloud"
            stub_gcloud.write_text(
                "#!/usr/bin/env bash\n"
                f"printf 'GCLOUD %s\\n' \"$*\" >> {str(inv_log)!r}\n",
                encoding="utf-8",
            )
            stub_gcloud.chmod(0o755)

            stub_runner = tmp_root / "harness" / "gcp_paper_runner.sh"
            stub_runner.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            stub_runner.chmod(0o755)

            git_env = {
                k: v
                for k, v in os.environ.items()
                if not k.startswith("GIT_") and k not in ("GITHUB_HEAD_REF", "GITHUB_REF")
            }
            subprocess.run(["git", "init", "-q", "-b", "main"], cwd=str(tmp_root), env=git_env, check=True)
            subprocess.run(["git", "config", "user.email", "t@example.com"], cwd=str(tmp_root), env=git_env, check=True)
            subprocess.run(["git", "config", "user.name", "Test"], cwd=str(tmp_root), env=git_env, check=True)
            (tmp_root / "README.md").write_text("init\n", encoding="utf-8")
            subprocess.run(["git", "add", "README.md"], cwd=str(tmp_root), env=git_env, check=True)
            subprocess.run(["git", "commit", "-qm", "init"], cwd=str(tmp_root), env=git_env, check=True)
            subprocess.run(["git", "checkout", "-q", "--detach", "HEAD"], cwd=str(tmp_root), env=git_env, check=True)

            def run_repro(*flags: str, extra_env: dict[str, str] | None = None) -> tuple[subprocess.CompletedProcess[str], str]:
                inv_log.write_text("", encoding="utf-8")
                run_env = dict(git_env)
                run_env["PATH"] = f"{tmp_root / 'mock_bin'}:{run_env.get('PATH', '')}"
                run_env["PAPER_ALLOW_STALE_PDF"] = "1"
                if extra_env:
                    run_env.update(extra_env)
                res = subprocess.run(
                    ["bash", str(copied_repro), *flags],
                    cwd=str(tmp_root),
                    env=run_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                return res, inv_log.read_text(encoding="utf-8")

            # 1) Detached HEAD without release env/flags: exits 0 and omits --require-head-fingerprint
            r_det, log_det = run_repro()
            self.assertEqual(r_det.returncode, 0, f"Detached HEAD run failed:\n{r_det.stderr}\n{r_det.stdout}")
            self.assertNotIn("GCLOUD ", log_det)
            self.assertIn("check_release_consistency.py --phase pre-tag", log_det)
            self.assertNotIn("--require-head-fingerprint", log_det)

            # 2) Both flag orders (--sync-gcs + --require-head-fingerprint) work identically
            for flags in (
                ("--sync-gcs", "--require-head-fingerprint"),
                ("--require-head-fingerprint", "--sync-gcs"),
            ):
                r_both, log_both = run_repro(*flags)
                self.assertEqual(r_both.returncode, 0, f"Flags {flags} failed:\n{r_both.stderr}")
                self.assertIn("GCLOUD storage rsync", log_both)
                self.assertIn("check_release_consistency.py --phase pre-tag --require-head-fingerprint", log_both)

            # 3) Release env vars and release/* local branch automatically enable --require-head-fingerprint
            for env_case in (
                {"GITHUB_HEAD_REF": "release/v1.0.11"},
                {"GITHUB_REF": "refs/tags/v1.0.11"},
            ):
                r_env, log_env = run_repro(extra_env=env_case)
                self.assertEqual(r_env.returncode, 0, f"Env {env_case} failed:\n{r_env.stderr}")
                self.assertIn("check_release_consistency.py --phase pre-tag --require-head-fingerprint", log_env)

            subprocess.run(["git", "checkout", "-q", "-b", "release/v1.0.11"], cwd=str(tmp_root), env=git_env, check=True)
            r_rel_br, log_rel_br = run_repro()
            self.assertEqual(r_rel_br.returncode, 0, f"release/* branch failed:\n{r_rel_br.stderr}")
            self.assertIn("check_release_consistency.py --phase pre-tag --require-head-fingerprint", log_rel_br)

            # 4) Unknown flag fails fast with exit 2
            r_bad, _ = run_repro("--bogus-flag")
            self.assertEqual(r_bad.returncode, 2, f"Expected exit 2 on unknown flag, got {r_bad.returncode}")
            self.assertIn("ERROR: unknown flag '--bogus-flag'", r_bad.stderr)

    def test_check_release_consistency_keys_recipe_license_on_recipe_pin_version(self):
        """Verify --version 1.1.0 --phase pre-tag keys Bioconda/nf-core MIT checks on each recipe's own pinned version."""
        import importlib.util

        script_path = (
            ROOT_DIR
            / ".agents"
            / "skills"
            / "vardictrust-release"
            / "scripts"
            / "check_release_consistency.py"
        )
        spec = importlib.util.spec_from_file_location("check_release_consistency", script_path)
        self.assertIsNotNone(spec)
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        # 0) Live ROOT_DIR packaging state must always satisfy check_pre_tag_surfaces:
        cargo_ver_m = re.search(
            r'^\[workspace\.package\]\s*\nversion\s*=\s*"([^"]+)"',
            (ROOT_DIR / "Cargo.toml").read_text(encoding="utf-8"),
            re.MULTILINE,
        )
        self.assertIsNotNone(cargo_ver_m)
        live_ver = cargo_ver_m.group(1)
        errs_live: list[str] = []
        mod.check_pre_tag_surfaces(ROOT_DIR, live_ver, errs_live, phase="pre-tag")
        pkg_errs_live = [
            e
            for e in errs_live
            if "packaging/bioconda/recipes/vardictrust/meta.yaml" in e
            or "packaging/nf-core/modules/vardictrust/" in e
        ]
        self.assertEqual(
            pkg_errs_live,
            [],
            f"Live ROOT_DIR packaging failed check_pre_tag_surfaces: {pkg_errs_live}",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_root = pathlib.Path(tmpdir)
            for item in ROOT_DIR.iterdir():
                if item.name in {".git", ".worktrees", "target", "scratch", ".venv", ".venv_docs"}:
                    continue
                if item.name == "packaging":
                    shutil.copytree(item, tmp_root / "packaging")
                else:
                    os.symlink(item, tmp_root / item.name)

            bc_meta = tmp_root / "packaging" / "bioconda" / "recipes" / "vardictrust" / "meta.yaml"
            nf_env = tmp_root / "packaging" / "nf-core" / "modules" / "vardictrust" / "environment.yml"
            nf_meta = tmp_root / "packaging" / "nf-core" / "modules" / "vardictrust" / "meta.yml"

            # Seed deterministic legacy 1.0.10 GPL dual-licensed recipe state into tmp_root
            # so this state-machine test remains invariant before and after Phase 3A/3B:
            bc_1010 = re.sub(
                r'\{%\s*set\s+version\s*=\s*"[^"]*"\s*%\}',
                '{% set version = "1.0.10" %}',
                bc_meta.read_text(encoding="utf-8"),
                count=1,
            )
            bc_1010 = re.sub(
                r"(?m)^([ \t]*license:[ \t]*).+$",
                r"\g<1>MIT AND GPL-2.0-or-later",
                bc_1010,
                count=1,
            )
            bc_1010 = re.sub(
                r"(?m)^([ \t]*license_family:[ \t]*).+$",
                r"\g<1>GPL2",
                bc_1010,
                count=1,
            )
            bc_meta.write_text(bc_1010, encoding="utf-8")

            nf_env.write_text(
                re.sub(
                    r"bioconda::vardictrust=[0-9]+\.[0-9]+\.[0-9]+",
                    "bioconda::vardictrust=1.0.10",
                    nf_env.read_text(encoding="utf-8"),
                    count=1,
                ),
                encoding="utf-8",
            )
            nf_meta.write_text(
                re.sub(
                    r'(^\s*-\s*"vardictrust":\s*\n(?:(?!^\s*-\s*"[A-Za-z0-9_-]+":)[^\n]*\n)*?'
                    r"^(\s*)licence:\s*\n)"
                    r"(\s*-\s*\"[^\n]+\"\s*\n)+",
                    lambda m: f'{m.group(1)}{m.group(2)}  - "MIT"\n{m.group(2)}  - "GPL-2.0-or-later"\n',
                    nf_meta.read_text(encoding="utf-8"),
                    count=1,
                    flags=re.MULTILINE,
                ),
                encoding="utf-8",
            )

            # 1) Simulate Phase 1 pre-tag for 1.1.0 while Bioconda & nf-core recipes remain pinned to 1.0.10:
            errs_pretag: list[str] = []
            mod.check_pre_tag_surfaces(tmp_root, "1.1.0", errs_pretag, phase="pre-tag")
            pkg_errs_pretag = [
                e
                for e in errs_pretag
                if "packaging/bioconda/recipes/vardictrust/meta.yaml" in e
                or "packaging/nf-core/modules/vardictrust/" in e
            ]
            self.assertEqual(
                pkg_errs_pretag,
                [],
                f"pre-tag 1.1.0 must accept 1.0.10 recipe licenses before Phase 3A/3B: {pkg_errs_pretag}",
            )

            # 2) Once recipes bump their own pinned version to 1.1.0, stale GPL dual-licensing is rejected:
            bc_meta.write_text(
                bc_meta.read_text(encoding="utf-8").replace(
                    '{% set version = "1.0.10" %}', '{% set version = "1.1.0" %}'
                ),
                encoding="utf-8",
            )
            nf_env.write_text(
                nf_env.read_text(encoding="utf-8").replace(
                    "bioconda::vardictrust=1.0.10", "bioconda::vardictrust=1.1.0"
                ),
                encoding="utf-8",
            )
            errs_stale_gpl: list[str] = []
            mod.check_pre_tag_surfaces(tmp_root, "1.1.0", errs_stale_gpl, phase="post-nfcore")
            self.assertTrue(
                any("packaging/bioconda/recipes/vardictrust/meta.yaml" in e for e in errs_stale_gpl),
                f"Expected Bioconda pure-MIT error once meta.yaml version is 1.1.0, got: {errs_stale_gpl}",
            )
            self.assertTrue(
                any("packaging/nf-core/modules/vardictrust/meta.yml" in e for e in errs_stale_gpl),
                f"Expected nf-core pure-MIT error once environment.yml version is 1.1.0, got: {errs_stale_gpl}",
            )

            # 3) Switching both 1.1.0 recipes to pure MIT clears all recipe license errors:
            bc_meta.write_text(
                bc_meta.read_text(encoding="utf-8")
                .replace("license: MIT AND GPL-2.0-or-later", "license: MIT")
                .replace("license_family: GPL2", "license_family: MIT"),
                encoding="utf-8",
            )
            nf_meta.write_text(
                re.sub(
                    r'licence:\s*\n\s*-\s*"MIT"\s*\n\s*-\s*"GPL-2\.0-or-later"',
                    'licence:\n    - "MIT"',
                    nf_meta.read_text(encoding="utf-8"),
                ),
                encoding="utf-8",
            )
            errs_clean_mit: list[str] = []
            mod.check_pre_tag_surfaces(tmp_root, "1.1.0", errs_clean_mit, phase="post-nfcore")
            pkg_errs_clean = [
                e
                for e in errs_clean_mit
                if "packaging/bioconda/recipes/vardictrust/meta.yaml" in e
                or "packaging/nf-core/modules/vardictrust/" in e
            ]
            self.assertEqual(pkg_errs_clean, [], f"Expected 0 recipe license errors on pure MIT 1.1.0: {pkg_errs_clean}")


if __name__ == "__main__":
    unittest.main()
