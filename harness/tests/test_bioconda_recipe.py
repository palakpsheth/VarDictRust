#!/usr/bin/env python3
"""
harness/tests/test_bioconda_recipe.py

Automated verification of the canonical Bioconda recipe for VarDictRust.
Conforms strictly to the Python Standard Library purity constraint (zero PyPI packages).
"""

import os
import pathlib
import unittest

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent
RECIPE_DIR = ROOT_DIR / "packaging" / "bioconda" / "recipes" / "vardictrust"
META_PATH = RECIPE_DIR / "meta.yaml"
BUILD_PATH = RECIPE_DIR / "build.sh"
RUN_TEST_PATH = RECIPE_DIR / "run_test.sh"


class TestBiocondaRecipe(unittest.TestCase):
    def setUp(self):
        self.assertTrue(RECIPE_DIR.is_dir(), f"Missing recipe directory: {RECIPE_DIR}")
        self.assertTrue(META_PATH.is_file(), f"Missing meta.yaml: {META_PATH}")
        self.assertTrue(BUILD_PATH.is_file(), f"Missing build.sh: {BUILD_PATH}")
        self.assertTrue(RUN_TEST_PATH.is_file(), f"Missing run_test.sh: {RUN_TEST_PATH}")

    def test_recipe_files_executable(self):
        self.assertTrue(os.access(BUILD_PATH, os.X_OK), "build.sh must be executable")
        self.assertTrue(os.access(RUN_TEST_PATH, os.X_OK), "run_test.sh must be executable")

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

    def test_run_test_negative_exit_contract(self):
        content = RUN_TEST_PATH.read_text(encoding="utf-8")
        self.assertIn(
            '[ "$status" -ne 2 ]',
            content,
            "run_test.sh must assert Clap usage error exit code 2 on invalid flags",
        )

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


if __name__ == "__main__":
    unittest.main()
