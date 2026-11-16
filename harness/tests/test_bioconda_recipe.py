#!/usr/bin/env python3
"""
harness/tests/test_bioconda_recipe.py

Automated verification of the canonical Bioconda recipe for VarDictRust.
Conforms strictly to the Python Standard Library purity constraint (zero PyPI packages).
"""

import importlib
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent
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
        cargo_path = ROOT_DIR / "Cargo.toml"
        content = cargo_path.read_text(encoding="utf-8")
        self.assertIn(
            'license = "MIT AND GPL-2.0-or-later"',
            content,
            "Root Cargo.toml must match the recipe SPDX license expression",
        )

    def test_cargo_workspace_msrv_inheritance(self):
        root_cargo = (ROOT_DIR / "Cargo.toml").read_text(encoding="utf-8")
        self.assertRegex(
            root_cargo,
            re.compile(
                r'^\[workspace\.package\]\n(?:(?!\[)[^\n]*\n)*?rust-version\s*=\s*"1\.91\.0"',
                re.MULTILINE,
            ),
            'Workspace Cargo.toml must pin rust-version = "1.91.0" under [workspace.package]',
        )
        for member in ("vardictrust", "vardictrust-core"):
            member_path = ROOT_DIR / "crates" / member / "Cargo.toml"
            self.assertTrue(member_path.is_file(), f"Missing crate manifest: {member_path}")
            member_cargo = member_path.read_text(encoding="utf-8")
            self.assertRegex(
                member_cargo,
                re.compile(r"^rust-version\.workspace\s*=\s*true\b", re.MULTILINE),
                f"crates/{member}/Cargo.toml must inherit rust-version.workspace = true",
            )

    def test_fuzz_manifest_noodles_version_matches_workspace(self):
        noodles_re = re.compile(
            r'^noodles\s*=\s*\{[^\n]*\bversion\s*=\s*"([^"]+)"',
            re.MULTILINE,
        )
        root_cargo = (ROOT_DIR / "Cargo.toml").read_text(encoding="utf-8")
        fuzz_path = ROOT_DIR / "crates" / "vardictrust-core" / "fuzz" / "Cargo.toml"
        self.assertTrue(fuzz_path.is_file(), f"Missing fuzz manifest: {fuzz_path}")
        fuzz_cargo = fuzz_path.read_text(encoding="utf-8")
        root_m = noodles_re.search(root_cargo)
        fuzz_m = noodles_re.search(fuzz_cargo)
        self.assertIsNotNone(root_m, "Root Cargo.toml missing noodles dependency entry")
        self.assertIsNotNone(
            fuzz_m,
            "crates/vardictrust-core/fuzz/Cargo.toml missing noodles dependency entry",
        )
        self.assertEqual(
            root_m.group(1),
            fuzz_m.group(1),
            "crates/vardictrust-core/fuzz/Cargo.toml noodles version must match root Cargo.toml",
        )

    def test_cargo_deny_license_allow_list(self):
        deny_path = ROOT_DIR / "deny.toml"
        content = deny_path.read_text(encoding="utf-8")
        self.assertIn(
            '"GPL-2.0-or-later"',
            content,
            "deny.toml must allow GPL-2.0-or-later",
        )
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
        self.assertIn("license: MIT AND GPL-2.0-or-later", meta, "License must be MIT AND GPL-2.0-or-later")
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

            # 2. Check GFM table column balance
            in_code = False
            expected_cols = None
            for idx, line in enumerate(body.splitlines(), 1):
                stripped = line.strip()
                if stripped.startswith("```"):
                    in_code = not in_code
                    continue
                if in_code:
                    continue
                if stripped.startswith("|") and stripped.endswith("|"):
                    cells = [c.strip() for c in re.split(r"(?<!\\)\|", stripped)]
                    if cells and cells[0] == "":
                        cells.pop(0)
                    if cells and cells[-1] == "":
                        cells.pop(-1)
                    if expected_cols is None:
                        expected_cols = len(cells)
                    else:
                        self.assertEqual(
                            len(cells),
                            expected_cols,
                            f"{rel_path}: table row at body line {idx} has {len(cells)} cols, expected {expected_cols}",
                        )
                else:
                    expected_cols = None

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

        cargo_toml = (ROOT_DIR / "Cargo.toml").read_text(encoding="utf-8")
        ver_m = re.search(r'^\[workspace\.package\]\s*\nversion\s*=\s*"([^"]+)"', cargo_toml, re.MULTILINE)
        self.assertIsNotNone(ver_m, "Could not parse workspace version from Cargo.toml")
        current_version = ver_m.group(1)

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
            5,
            "All 5 downstream jobs in release.yml must check out needs.plan.outputs.checkout-ref",
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

        dist_toml = (ROOT_DIR / "dist-workspace.toml").read_text(encoding="utf-8")
        self.assertIn("[dist.github-custom-runners]", dist_toml)
        self.assertIn('aarch64-apple-darwin = "macos-15"', dist_toml)
        for retired_runner in ("macos-12", "macos-13", "macos-14"):
            self.assertNotIn(retired_runner, dist_toml)

        deny_text = (ROOT_DIR / "deny.toml").read_text(encoding="utf-8")
        self.assertIn('yanked = "deny"', deny_text)
        self.assertIn('wildcards = "deny"', deny_text)
        self.assertIn('required-git-spec = "rev"', deny_text)

        dependabot_text = (ROOT_DIR / ".github" / "dependabot.yml").read_text(encoding="utf-8")
        for eco in ('package-ecosystem: "cargo"', 'package-ecosystem: "github-actions"', 'package-ecosystem: "pip"'):
            self.assertIn(eco, dependabot_text)
        self.assertEqual(dependabot_text.count("cooldown:"), 3, "dependabot.yml must configure cooldown on all 3 ecosystems")
        self.assertEqual(dependabot_text.count("default-days: 7"), 3, "dependabot.yml must set default-days: 7 on all 3 ecosystems")

        nf_main_text = (ROOT_DIR / "packaging" / "nf-core" / "modules" / "vardictrust" / "main.nf").read_text(encoding="utf-8")
        self.assertEqual(
            nf_main_text.count("prefix ==~ /^[A-Za-z0-9._-]+$/"),
            2,
            "main.nf must validate task.ext.prefix in both script: and stub:",
        )
        self.assertEqual(
            nf_main_text.count('if (!("${meta.id}" ==~ /^[A-Za-z0-9._|-]+$/)) { error "Unsafe meta.id: ${meta.id}" }'),
            2,
            "main.nf must validate meta.id independently in both script: and stub:",
        )
        self.assertIn("def shQuote =", nf_main_text)
        self.assertIn("def fasta_arg = shQuote(fasta)", nf_main_text)
        self.assertIn("-G ${fasta_arg}", nf_main_text)
        self.assertEqual(
            nf_main_text.count('> "${prefix}.vcf.gz"'),
            2,
            'main.nf must quote > "${prefix}.vcf.gz" in both script: and stub:',
        )

        cfg_text = (ROOT_DIR / "harness" / "configure_github_repo.py").read_text(encoding="utf-8")
        self.assertIn('"has_wiki": True', cfg_text)
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

    def test_release_workflow_cargo_dist_pinned_shas_and_harden_runner(self):
        """Verify release.yml pins 4-arch cargo-dist SHA-256 digests, eliminates | sh, and hardens runners."""
        rel_path = ROOT_DIR / ".github" / "workflows" / "release.yml"
        self.assertTrue(rel_path.is_file(), "Missing .github/workflows/release.yml")
        rel_text = rel_path.read_text(encoding="utf-8")

        # Zero pipe-to-shell invocations
        self.assertNotIn("| sh", rel_text, "release.yml must not contain unpinned '| sh' pipe-to-shell invocations")
        self.assertNotIn("cargo-dist-installer.sh", rel_text, "release.yml must not use remote cargo-dist-installer.sh script")
        self.assertNotIn("sh.rustup.rs", rel_text, "release.yml must not use pipe-to-shell rustup installer")
        self.assertNotIn("${{ matrix.install_dist.run }}", rel_text, "release.yml must not execute unpinned matrix.install_dist.run")

        # Pinned 4-arch cargo-dist v0.33.0 SHA-256 digests
        expected_shas = {
            "Linux-x86_64": "4b3f0a5f0ebbdb798f6db649d01b32ba1518376b6f7a0502b7d92b75cc2c8293",
            "Linux-aarch64": "9c554ab21a58ad46eb9b6710f89633ccb26bc2577cf4b0c2d18a5b826a31913e",
            "Darwin-x86_64": "6a49bfb61bd86770d79c27f3d2b40c6b2e71cde940d3d31a6ccaaffc124d7a29",
            "Darwin-arm64": "7b3cbe25511de01d74c0f5fcb7909edabd379bea9cfa284d93af5a3cdfa3247c",
        }
        for arch_label, digest in expected_shas.items():
            self.assertIn(
                digest,
                rel_text,
                f"release.yml must pin {arch_label} cargo-dist v0.33.0 SHA-256 digest {digest}",
            )

        # Archive extraction flags, cargo/bin destination, path propagation, and cache salting
        self.assertIn("--strip-components=1", rel_text, "release.yml must extract cargo-dist tar.xz with --strip-components=1")
        self.assertIn("$HOME/.cargo/bin", rel_text, "release.yml must install cargo-dist binary to $HOME/.cargo/bin")
        self.assertIn("GITHUB_PATH", rel_text, "release.yml must append $HOME/.cargo/bin to GITHUB_PATH")
        self.assertIn("cargo-dist-archive-sha256-v0.33.0-", rel_text, "release.yml must salt cargo-dist cache key with cargo-dist-archive-sha256-v0.33.0-")
        self.assertIn("~/.cache/cargo-dist", rel_text, "release.yml must cache cargo-dist archive directory")

        # Verify dist --version never appears before verify_tar in plan or build-local-artifacts
        for job_name in ("plan", "build-local-artifacts"):
            job_split = rel_text.split(f"  {job_name}:")
            self.assertEqual(len(job_split), 2, f"release.yml missing job {job_name}")
            job_body = job_split[1].split("\n  build-")[0] if job_name == "plan" else job_split[1].split("\n  publish-")[0]
            self.assertIn("verify_tar", job_body, f"job {job_name} must define verify_tar")
            verify_idx = job_body.index("verify_tar")
            self.assertIn('"$HOME/.cargo/bin/dist" --version', job_body, f"job {job_name} must run dist --version")
            dist_ver_idx = job_body.index('"$HOME/.cargo/bin/dist" --version')
            self.assertGreater(dist_ver_idx, verify_idx, f"job {job_name} must not run dist --version before verify_tar")
            self.assertNotIn('exit 0', job_body[:verify_idx], f"job {job_name} must not have early exit shortcut before verify_tar")

        # harden-runner in all 6 jobs
        harden_action = "step-security/harden-runner@1d88384fd173af6ecb012dfd0264a306b711241a"
        self.assertEqual(
            rel_text.count(harden_action),
            6,
            f"release.yml must include harden-runner ({harden_action}) as first step in all 6 jobs",
        )
        self.assertIn("egress-policy: audit", rel_text, "harden-runner must configure egress-policy: audit")
        self.assertIn(
            "if: runner.os == 'Linux' && runner.arch == 'X64'",
            rel_text,
            "harden-runner on cross-OS jobs must be gated with 'if: runner.os == 'Linux' && runner.arch == 'X64''",
        )

    def test_configure_github_repo_tag_ruleset_and_security_settings(self):
        """Verify configure_github_repo.py API headers, merge commit restrictions, and lockout-safe tag ruleset."""
        cfg_path = ROOT_DIR / "harness" / "configure_github_repo.py"
        self.assertTrue(cfg_path.is_file(), "Missing harness/configure_github_repo.py")
        cfg_text = cfg_path.read_text(encoding="utf-8")

        self.assertIn(
            '"Accept": "application/vnd.github+json"',
            cfg_text,
            "configure_github_repo.py must use modern 'application/vnd.github+json' Accept header",
        )
        self.assertIn(
            '"allow_merge_commit": False',
            cfg_text,
            "configure_github_repo.py must set 'allow_merge_commit': False",
        )
        self.assertIn(
            '"MkDocs Strict Build Gate"',
            cfg_text,
            "configure_github_repo.py must retain 'MkDocs Strict Build Gate' in required status checks",
        )
        self.assertIn("def configure_tag_ruleset", cfg_text, "configure_github_repo.py must define configure_tag_ruleset()")
        self.assertIn("current_user_can_bypass", cfg_text, "configure_tag_ruleset must check current_user_can_bypass")
        self.assertIn("24867690", cfg_text, "configure_tag_ruleset must target existing ruleset ID 24867690")
        self.assertIn('"refs/tags/v*"', cfg_text, "configure_tag_ruleset must include 'refs/tags/v*'")
        self.assertNotIn('"~refs/tags/v*"', cfg_text, "configure_tag_ruleset must not use invalid '~refs/tags/v*'")

        help_run = subprocess.run(
            [sys.executable, str(cfg_path), "--help"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(help_run.returncode, 0, f"--help failed: {help_run.stderr}")
        self.assertIn("--dry-run", help_run.stdout)
        self.assertIn("--apply", help_run.stdout)

        dry_run = subprocess.run(
            [sys.executable, str(cfg_path), "--dry-run"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(dry_run.returncode, 0, f"--dry-run failed: {dry_run.stderr}")
        self.assertIn("Configured ruleset 24867690 for refs/tags/v*", dry_run.stdout)

        # Verify configure_tag_ruleset rolls back to safe payload if GET raises exception (Copilot 4580868305)
        import importlib.util
        import urllib.error
        from unittest.mock import patch

        spec = importlib.util.spec_from_file_location("configure_github_repo_mod", str(cfg_path))
        cfg_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cfg_mod)

        recorded_calls = []

        def mock_api_request(endpoint, method="GET", data=None, dry_run=False):
            recorded_calls.append((endpoint, method, data, dry_run))
            if method == "GET":
                raise urllib.error.URLError("Simulated connection timeout")
            return {"status": "ok"}

        with patch.object(cfg_mod, "api_request", side_effect=mock_api_request):
            cfg_mod.configure_tag_ruleset(dry_run=False)

        # Should have: 1. PUT strict, 2. GET (which raises), 3. PUT safe rollback
        self.assertEqual(len(recorded_calls), 3, f"Expected 3 API calls, got {len(recorded_calls)}")
        self.assertEqual(recorded_calls[0][1], "PUT")
        self.assertTrue(any(r.get("type") == "creation" for r in recorded_calls[0][2].get("rules", [])))
        self.assertEqual(recorded_calls[1][1], "GET")
        self.assertEqual(recorded_calls[2][1], "PUT")
        self.assertFalse(any(r.get("type") == "creation" for r in recorded_calls[2][2].get("rules", [])))

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
        import json

        dist_bin = shutil.which("dist")
        if not dist_bin:
            candidate = pathlib.Path.home() / ".cargo" / "bin" / "dist"
            if candidate.is_file() and os.access(candidate, os.X_OK):
                dist_bin = str(candidate)
        if not dist_bin:
            self.skipTest("cargo-dist ('dist') binary not installed on host; covered by release.yml PR job")

        cargo_toml = (ROOT_DIR / "Cargo.toml").read_text(encoding="utf-8")
        ver_m = re.search(
            r'^\[workspace\.package\]\n(?:(?!\[)[^\n]*\n)*?version\s*=\s*"([^"]+)"',
            cargo_toml,
            re.MULTILINE,
        )
        self.assertIsNotNone(ver_m, "Could not parse workspace version from Cargo.toml")
        ws_ver = ver_m.group(1)
        tag_arg = f"--tag=v{ws_ver}"

        # Strip GIT_* hook variables (GIT_DIR, GIT_INDEX_FILE, GIT_WORK_TREE) so subprocess git/dist
        # invocations work identically inside pre-commit/pre-push hooks and read-only mounts.
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
        import base64
        import http.server
        import json
        import threading

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
        import textwrap

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
                # GitHub's POST /repos/{owner}/{repo}/attestations validator enforces that
                # resolvedDependencies[0] matches externalParameters.workflow (Fulcio cert ref/sha).
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


if __name__ == "__main__":
    unittest.main()
