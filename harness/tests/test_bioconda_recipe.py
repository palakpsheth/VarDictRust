#!/usr/bin/env python3
"""
harness/tests/test_bioconda_recipe.py

Automated verification of the canonical Bioconda recipe for VarDictRust.
Conforms strictly to the Python Standard Library purity constraint (zero PyPI packages).
"""

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


    def test_cargo_manifest_license(self):
        cargo_path = ROOT_DIR / "Cargo.toml"
        content = cargo_path.read_text(encoding="utf-8")
        self.assertIn(
            'license = "MIT AND GPL-2.0-or-later"',
            content,
            "Root Cargo.toml must match the recipe SPDX license expression",
        )

    def test_cargo_deny_license_allow_list(self):
        deny_path = ROOT_DIR / "deny.toml"
        content = deny_path.read_text(encoding="utf-8")
        self.assertIn(
            '"GPL-2.0-or-later"',
            content,
            "deny.toml must allow GPL-2.0-or-later",
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
        self.assertIn("run_exports:", meta, "Bioconda linter requires run_exports pinning")
        self.assertIn("pin_subpackage", meta, "run_exports must use pin_subpackage")
        self.assertIn("binary_relocation: false", meta, "meta.yaml must declare binary_relocation: false")
        self.assertIn("osx-arm64", meta, "Must target osx-arm64")
        self.assertIn("linux-aarch64", meta, "Must target linux-aarch64")

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


if __name__ == "__main__":
    unittest.main()
