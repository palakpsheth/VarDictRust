#!/usr/bin/env python3
"""VarDictRust Release Consistency & Evidence Validator (stdlib-only).

Verifies that all version-bearing files, workspace test counts (619 total:
432 core lib + 110 corpus_parity + 29 engine_behavior + 2 amplicon_stream_mem
+ 3 cram_stack + 24 CLI bin + 19 CLI integration binary tests), engine
fingerprints (FP & PERF_FP), and 16 benchmarks/evidence/*.json SHA-256
manifests are synchronized for a target release version and phase.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

VALID_PHASES = ("pre-tag", "post-tag", "post-bioconda", "post-nfcore")

EXPECTED_TEST_COUNTS = {
    "total_workspace": 619,
    "core_lib": 432,
    "corpus_parity": 110,
    "engine_behavior": 29,
    "amplicon_stream_mem": 2,
    "cram_stack": 3,
    "cli_bin": 24,
    "cli_integration": 19,
    "harness_self_tests": 36,
}


def read_text(repo_root: Path, rel_path: str, errors: list[str]) -> str:
    path = repo_root / rel_path
    if not path.is_file():
        errors.append(f"MISSING FILE: {rel_path}")
        return ""
    return path.read_text(encoding="utf-8")


def require_regex(
    repo_root: Path,
    rel_path: str,
    pattern: str,
    description: str,
    errors: list[str],
    flags: int = re.MULTILINE,
) -> None:
    content = read_text(repo_root, rel_path, errors)
    if not content:
        return
    if not re.search(pattern, content, flags):
        errors.append(f"{rel_path}: failed check [{description}] (pattern: {pattern!r})")


def forbid_regex(
    repo_root: Path,
    rel_path: str,
    pattern: str,
    description: str,
    errors: list[str],
    flags: int = re.MULTILINE,
) -> None:
    content = read_text(repo_root, rel_path, errors)
    if not content:
        return
    match = re.search(pattern, content, flags)
    if match:
        errors.append(
            f"{rel_path}: found forbidden pattern [{description}]: {match.group(0)!r}"
        )


EXPECTED_CLI_HEADINGS = (
    "1. Core Inputs & Outputs",
    "2. VCF Post-Filtering (var2vcf)",
    "3. Coordinate & BED Column Mapping",
    "4. Variant Calling & Quality Thresholds",
    "5. Modality Flags: Somatic & Amplicon",
    "6. Structural Variants & Complex Indels",
    "7. Performance & Concurrency",
)
EXPECTED_CLI_FLAG_COUNT = 89


def check_cli_and_migration_surfaces(repo_root: Path, errors: list[str]) -> None:
    n_flags = EXPECTED_CLI_FLAG_COUNT
    n_groups = len(EXPECTED_CLI_HEADINGS)

    require_regex(
        repo_root,
        "README.md",
        rf"`{n_flags}` CLI flags across {n_groups} functional groups",
        f"README.md cites {n_flags} CLI flags across {n_groups} functional groups",
        errors,
    )
    forbid_regex(
        repo_root,
        "README.md",
        r"`87` CLI flags|6 functional groups",
        "README.md stale 87 CLI flags or 6 functional groups",
        errors,
    )

    readme_text = read_text(repo_root, "README.md", errors)
    if readme_text:
        sec_m = re.search(
            r"^## CLI Parameter Reference\n(.*?)(?=^## |\Z)",
            readme_text,
            re.MULTILINE | re.DOTALL,
        )
        if not sec_m:
            errors.append("README.md: missing '## CLI Parameter Reference' section")
        else:
            cli_sec = sec_m.group(1)
            found_headings = re.findall(r"^### (\d+\. [^\n]+)$", cli_sec, re.MULTILINE)
            if tuple(found_headings) != EXPECTED_CLI_HEADINGS:
                errors.append(
                    f"README.md: CLI Parameter Reference headings {found_headings} != {list(EXPECTED_CLI_HEADINGS)}"
                )
            data_rows = [
                ln
                for ln in cli_sec.splitlines()
                if ln.startswith("| ")
                and not ln.startswith("| Short | Long |")
                and not ln.startswith("| :---")
            ]
            if len(data_rows) != n_flags:
                errors.append(
                    f"README.md: CLI Parameter Reference tables contain {len(data_rows)} flag rows, expected {n_flags}"
                )

    roadmap_text = read_text(repo_root, "ROADMAP.md", errors)
    if roadmap_text:
        matches = re.findall(
            rf"`{n_flags}`-flag CLI surface across `{n_groups}` functional groups",
            roadmap_text,
        )
        if len(matches) < 2:
            errors.append(
                f"ROADMAP.md: expected >= 2 occurrences of '`{n_flags}`-flag CLI surface across `{n_groups}` functional groups', found {len(matches)}"
            )
    forbid_regex(
        repo_root,
        "ROADMAP.md",
        r"`87`-flag CLI surface|6 functional groups",
        "ROADMAP.md stale 87-flag CLI surface or 6 functional groups",
        errors,
    )

    forbid_regex(
        repo_root,
        "CHANGELOG.md",
        r"`87` total CLI flags across 6 functional groups",
        "CHANGELOG.md stale 87 total CLI flags across 6 functional groups",
        errors,
    )

    require_regex(
        repo_root,
        "docs/ecosystem-survey.md",
        rf"Explicit {n_flags}-flag Diátaxis reference categorized into {n_groups} functional groups",
        f"docs/ecosystem-survey.md §5 cites {n_flags}-flag reference across {n_groups} functional groups",
        errors,
    )
    require_regex(
        repo_root,
        "docs/ecosystem-survey.md",
        rf"`{n_flags}` CLI flags across `{n_groups}` functional groups",
        f"docs/ecosystem-survey.md §8.4 cites {n_flags} CLI flags across {n_groups} functional groups",
        errors,
    )
    forbid_regex(
        repo_root,
        "docs/ecosystem-survey.md",
        r"87-flag|`87` CLI flags|6 functional groups",
        "docs/ecosystem-survey.md stale 87-flag or 6 functional groups",
        errors,
    )

    require_regex(
        repo_root,
        "docs/testing/test-matrix-and-ci-audit.md",
        rf"across the {n_flags}-flag CLI surface \({n_groups} functional groups",
        f"test-matrix-and-ci-audit.md cites {n_flags}-flag CLI surface ({n_groups} functional groups)",
        errors,
    )
    forbid_regex(
        repo_root,
        "docs/testing/test-matrix-and-ci-audit.md",
        r"87-flag CLI surface|6 functional groups",
        "test-matrix-and-ci-audit.md stale 87-flag CLI surface or 6 functional groups",
        errors,
    )

    require_regex(
        repo_root,
        "wiki/Home.md",
        rf"Complete {n_flags}-flag reference across {n_groups} functional groups",
        f"wiki/Home.md cites {n_flags}-flag reference across {n_groups} functional groups",
        errors,
    )
    forbid_regex(
        repo_root,
        "wiki/Home.md",
        r"87-flag|6 functional groups",
        "wiki/Home.md stale 87-flag or 6 functional groups",
        errors,
    )

    migration_text = read_text(repo_root, "docs/migration.md", errors)
    if migration_text and "var2vcf_valid.pl" not in migration_text:
        errors.append("docs/migration.md: missing var2vcf_valid.pl migration table")
    require_regex(
        repo_root,
        "mkdocs.yml",
        r"^\s*-\s*Migration Guide:\s*docs/migration\.md\s*$",
        "mkdocs.yml nav includes Migration Guide: docs/migration.md",
        errors,
    )
    require_regex(
        repo_root,
        "wiki/_Sidebar.md",
        r"docs/migration\.md",
        "wiki/_Sidebar.md links to docs/migration.md",
        errors,
    )


def check_pre_tag_surfaces(
    repo_root: Path,
    version: str,
    errors: list[str],
    *,
    require_head_fp: bool = False,
) -> None:
    ver_esc = re.escape(version)

    # 1. Cargo.toml & MSRV sync
    require_regex(
        repo_root,
        "Cargo.toml",
        rf'^\[workspace\.package\]\s*\nversion\s*=\s*"{ver_esc}"',
        f"workspace.package version == {version}",
        errors,
    )
    cargo_toml_text = read_text(repo_root, "Cargo.toml", errors)
    if cargo_toml_text:
        msrv_m = re.search(
            r'^\[workspace\.package\](?:\n(?!\[)[^\n]*)*?\nrust-version\s*=\s*"([^"]+)"',
            cargo_toml_text,
            re.MULTILINE,
        )
        if not msrv_m:
            errors.append("Cargo.toml: missing rust-version under [workspace.package]")
        else:
            msrv = msrv_m.group(1)
            if msrv != "1.91.0":
                errors.append(
                    f'Cargo.toml: expected rust-version == "1.91.0", got {msrv!r}'
                )
            msrv_esc = re.escape(msrv)
            for member_toml in (
                "crates/vardictrust/Cargo.toml",
                "crates/vardictrust-core/Cargo.toml",
            ):
                require_regex(
                    repo_root,
                    member_toml,
                    r"^rust-version\.workspace\s*=\s*true\b",
                    f"{member_toml} inherits rust-version.workspace = true",
                    errors,
                )
            for doc_path in ("README.md", "CONTRIBUTING.md"):
                require_regex(
                    repo_root,
                    doc_path,
                    rf"rustc >= {msrv_esc}\b",
                    f"{doc_path} cites MSRV rustc >= {msrv}",
                    errors,
                )

    # 2. crates/vardictrust/Cargo.toml
    require_regex(
        repo_root,
        "crates/vardictrust/Cargo.toml",
        rf'vardictrust-core\s*=\s*\{{[^}}]*version\s*=\s*"{ver_esc}"',
        f"vardictrust-core dependency version == {version}",
        errors,
    )

    # 3. Cargo.lock
    for pkg in ("vardictrust", "vardictrust-core"):
        require_regex(
            repo_root,
            "Cargo.lock",
            rf'name\s*=\s*"{pkg}"\s*\nversion\s*=\s*"{ver_esc}"',
            f"Cargo.lock {pkg} version == {version}",
            errors,
        )

    # 4. CITATION.cff
    require_regex(
        repo_root,
        "CITATION.cff",
        rf"^version:\s*{ver_esc}\s*$",
        f"CITATION.cff version == {version}",
        errors,
    )

    # 5. CHANGELOG.md
    require_regex(
        repo_root,
        "CHANGELOG.md",
        r"^## \[Unreleased\]",
        "CHANGELOG.md has ## [Unreleased] header",
        errors,
    )
    require_regex(
        repo_root,
        "CHANGELOG.md",
        rf"^## \[{ver_esc}\]\s*-\s*\d{{4}}-\d{{2}}-\d{{2}}",
        f"CHANGELOG.md has ## [{version}] release section",
        errors,
    )
    require_regex(
        repo_root,
        "CHANGELOG.md",
        rf"^\[Unreleased\]:\s*https://github\.com/palakpsheth/VarDictRust/compare/v{ver_esc}\.\.\.HEAD",
        f"CHANGELOG.md [Unreleased] compares v{version}...HEAD",
        errors,
    )
    require_regex(
        repo_root,
        "CHANGELOG.md",
        rf"^\[{ver_esc}\]:\s*https://github\.com/palakpsheth/VarDictRust/compare/v\d+\.\d+\.\d+\.\.\.v{ver_esc}",
        f"CHANGELOG.md has [{version}] comparison link",
        errors,
    )
    require_regex(
        repo_root,
        "CHANGELOG.md",
        r"^\[1\.0\.2\]:\s*https://github\.com/palakpsheth/VarDictRust/compare/v1\.0\.1\.\.\.v1\.0\.2",
        "CHANGELOG.md has [1.0.2] comparison link",
        errors,
    )
    changelog_text = read_text(repo_root, "CHANGELOG.md", errors)
    if changelog_text:
        rel_headers = re.findall(
            r"^## \[(\d+\.\d+\.\d+)\]", changelog_text, re.MULTILINE
        )
        if not rel_headers:
            errors.append("CHANGELOG.md: found zero ## [X.Y.Z] release headers")
        for rel_ver in rel_headers:
            rel_esc = re.escape(rel_ver)
            if rel_ver == "1.0.0":
                footer_pat = r"^\[1\.0\.0\]:\s*https://github\.com/palakpsheth/VarDictRust/releases/tag/v1\.0\.0\s*$"
            else:
                footer_pat = rf"^\[{rel_esc}\]:\s*https://github\.com/palakpsheth/VarDictRust/compare/v\d+\.\d+\.\d+\.\.\.v{rel_esc}\s*$"
            if not re.search(footer_pat, changelog_text, re.MULTILINE):
                errors.append(
                    f"CHANGELOG.md: ## [{rel_ver}] header missing matching [{rel_ver}]: footer link"
                )
    forbid_regex(
        repo_root,
        "CHANGELOG.md",
        r"dirty_decay_ms:5000,\s*muzzy_decay_ms:5000|--max-ref-cache-mb|SingleRefSequenceProvider|round_half_even_f64_val|SmallVec<\[\(usize, usize\); 2\]>|dirty_decay_ms:0,narenas:4|harness/lib/oracle_platform\.sh|JavaHashMap<K:|m2f2",
        "stale MALLOC_CONF or superseded PR draft symbols in CHANGELOG.md",
        errors,
    )

    # 6. .github/ISSUE_TEMPLATE/bug_report.yml
    require_regex(
        repo_root,
        ".github/ISSUE_TEMPLATE/bug_report.yml",
        rf'placeholder:\s*"?vardictrust {ver_esc}\b',
        f"bug_report.yml placeholder version == {version}",
        errors,
    )

    # 7 & 8. Packaging helper scripts
    require_regex(
        repo_root,
        "packaging/bioconda/prepare_bioconda_pr.sh",
        rf'^VERSION="\$\{{1:-{ver_esc}\}}"',
        f"prepare_bioconda_pr.sh default VERSION == {version}",
        errors,
    )
    require_regex(
        repo_root,
        "packaging/nf-core/prepare_nfcore_pr.sh",
        rf'^VERSION="\$\{{1:-{ver_esc}\}}"',
        f"prepare_nfcore_pr.sh default VERSION == {version}",
        errors,
    )

    # 9. harness/run_self_tests.sh & harness/tests/
    tc = EXPECTED_TEST_COUNTS
    sh_cnt = len(list((repo_root / "harness" / "tests").glob("test_*.sh")))
    py_cnt = len(list((repo_root / "harness" / "tests").glob("test_*.py")))
    total_self_tests = sh_cnt + py_cnt
    if total_self_tests != tc["harness_self_tests"]:
        errors.append(
            f"harness/tests/ has {total_self_tests} scripts ({sh_cnt} .sh + {py_cnt} .py), "
            f"expected {tc['harness_self_tests']}"
        )
    require_regex(
        repo_root,
        "harness/run_self_tests.sh",
        rf"^# VarDictRust v{ver_esc}\b",
        f"harness/run_self_tests.sh header == v{version}",
        errors,
    )
    require_regex(
        repo_root,
        "harness/run_self_tests.sh",
        rf"all {tc['harness_self_tests']} harness diagnostic scripts \(25 shell tests\s*\n#\s*and 11 python test suites",
        f"harness/run_self_tests.sh cites {tc['harness_self_tests']} diagnostic scripts (25 shell + 11 python)",
        errors,
    )

    # 10. docs/release/vX.Y.Z-pr.md
    require_regex(
        repo_root,
        f"docs/release/v{version}-pr.md",
        rf"# Release v{ver_esc}:",
        f"docs/release/v{version}-pr.md title header",
        errors,
    )
    forbid_regex(
        repo_root,
        f"docs/release/v{version}-pr.md",
        r"--max-ref-cache-mb|SingleRefSequenceProvider|round_half_even_f64_val|SmallVec<\[\(usize, usize\); 2\]>|dirty_decay_ms:0,narenas:4|harness/lib/oracle_platform\.sh|JavaHashMap<K:|m2f2",
        f"docs/release/v{version}-pr.md superseded PR draft symbols",
        errors,
    )

    # 11. README.md
    require_regex(
        repo_root,
        "README.md",
        rf"--tag v{ver_esc}\b",
        f"README.md cargo install --tag v{version}",
        errors,
    )
    require_regex(
        repo_root,
        "README.md",
        rf"version\s*=\s*\{{{ver_esc}\}}",
        f"README.md BibTeX version == {version}",
        errors,
    )
    require_regex(
        repo_root,
        "README.md",
        (
            rf"{tc['total_workspace']} automated Cargo workspace tests "
            rf"\({tc['core_lib']} core library tests, "
            rf"{tc['corpus_parity']} corpus and modality parity tests, "
            rf"{tc['engine_behavior']} engine behavior tests, "
            rf"{tc['amplicon_stream_mem']} amplicon streaming memory tests, "
            rf"{tc['cram_stack']} CRAM stack tests, "
            rf"{tc['cli_bin']} CLI tests, and "
            rf"{tc['cli_integration']} CLI integration binary tests\)"
        ),
        (
            f"README.md exact {tc['total_workspace']} test breakdown "
            f"({tc['core_lib']} + {tc['corpus_parity']} + {tc['engine_behavior']} + "
            f"{tc['amplicon_stream_mem']} + {tc['cram_stack']} + {tc['cli_bin']} + "
            f"{tc['cli_integration']})"
        ),
        errors,
    )
    require_regex(
        repo_root,
        "README.md",
        rf"{tc['harness_self_tests']} harness diagnostic self-tests",
        f"README.md cites {tc['harness_self_tests']} harness diagnostic self-tests",
        errors,
    )
    forbid_regex(
        repo_root,
        "README.md",
        r"4\.3×–6\.9×|36\.4×|3\.48×|1\.42 GiB|3\.84×|2:23\.49|1:17:34|2:03\.11|13:04\.76|3:24\.57|25 harness diagnostic self-tests",
        "README.md stale v1.0.6 benchmark metrics or 25 self-tests",
        errors,
    )

    # 12. DIVERGENCES.md
    hg008_ev = repo_root / "benchmarks" / "evidence" / "hg008-tier2.json"
    baseline_fp = "33bc15964cc36b6d"
    if hg008_ev.is_file():
        try:
            baseline_fp = json.loads(hg008_ev.read_text(encoding="utf-8")).get(
                "fingerprint", baseline_fp
            )
        except Exception:
            pass
    require_regex(
        repo_root,
        "DIVERGENCES.md",
        rf"Verified Parity Baseline \(release/v{ver_esc}\b",
        f"DIVERGENCES.md cites Verified Parity Baseline (release/v{version}",
        errors,
    )
    require_regex(
        repo_root,
        "DIVERGENCES.md",
        r"##source=VarDict_v1\.8\.2",
        "DIVERGENCES.md documents ##source=VarDict_v1.8.2 byte-parity header",
        errors,
    )
    if require_head_fp:
        forbid_regex(
            repo_root,
            "DIVERGENCES.md",
            r"Current Workspace Test Suite \(main / pre-v",
            "DIVERGENCES.md pre-release Current Workspace Test Suite block must be folded into Verified Parity Baseline on release branches",
            errors,
        )
        require_regex(
            repo_root,
            "DIVERGENCES.md",
            rf"Verified Parity Baseline \(release/v{ver_esc}[^\n]*:(?:\n- [^\n]*)*?\n- Core Library Tests:\s+{tc['core_lib']} PASS[^\n]*\n- Total Workspace:\s+{tc['total_workspace']} PASS",
            f"DIVERGENCES.md Verified Parity Baseline (release/v{version}) cites {tc['core_lib']} core lib / {tc['total_workspace']} workspace tests",
            errors,
            flags=re.MULTILINE,
        )
    else:
        require_regex(
            repo_root,
            "DIVERGENCES.md",
            rf"Current Workspace Test Suite[^\n]*:\n- Core Library Tests:\s+{tc['core_lib']} PASS[^\n]*\n- Total Workspace:\s+{tc['total_workspace']} PASS",
            f"DIVERGENCES.md cites current {tc['core_lib']} core lib / {tc['total_workspace']} workspace tests",
            errors,
            flags=re.MULTILINE,
        )
    require_regex(
        repo_root,
        "DIVERGENCES.md",
        rf"{re.escape(baseline_fp)}",
        f"DIVERGENCES.md cites certified fingerprint {baseline_fp}",
        errors,
    )
    forbid_regex(
        repo_root,
        "DIVERGENCES.md",
        r"817d89c5fb34515d|profile_features|3\.58x speedup|1\.42 GiB|3\.48x speedup|affb71e",
        "DIVERGENCES.md stale fingerprint/test name/speedup/RSS/commit",
        errors,
    )

    # 13. CONTRIBUTING.md
    require_regex(
        repo_root,
        "CONTRIBUTING.md",
        rf"all {tc['total_workspace']} workspace unit and integration tests",
        f"CONTRIBUTING.md §5.4 cites {tc['total_workspace']} workspace tests",
        errors,
    )
    require_regex(
        repo_root,
        "CONTRIBUTING.md",
        r"aarch64-unknown-linux-musl",
        "CONTRIBUTING.md §9.1 cites aarch64-unknown-linux-musl",
        errors,
    )
    forbid_regex(
        repo_root,
        "CONTRIBUTING.md",
        r"aarch64-unknown-linux-gnu|all 233 workspace|all 277 workspace|all 289 workspace|all 293 workspace|all 294 workspace|all 296 workspace",
        "CONTRIBUTING.md stale target triple or 233/277/289/293/294/296 test count",
        errors,
    )

    # 14. benchmarks/README.md & harness/README.md
    require_regex(
        repo_root,
        "benchmarks/README.md",
        rf"v{ver_esc}",
        f"benchmarks/README.md cites v{version}",
        errors,
    )
    require_regex(
        repo_root,
        "harness/README.md",
        rf"{tc['harness_self_tests']} diagnostic test scripts",
        f"harness/README.md cites {tc['harness_self_tests']} diagnostic test scripts",
        errors,
    )
    forbid_regex(
        repo_root,
        "harness/README.md",
        r"\b427\b|25 diagnostic test scripts|25 test scripts|17 Shell Tests|8 Python Test Suites|4-target performance budget",
        "harness/README.md stale test counts or 4-target budget",
        errors,
    )

    # 15. docs/testing/test-matrix-and-ci-audit.md
    require_regex(
        repo_root,
        "docs/testing/test-matrix-and-ci-audit.md",
        rf"v{ver_esc}",
        f"test-matrix-and-ci-audit.md cites v{version}",
        errors,
    )
    require_regex(
        repo_root,
        "docs/testing/test-matrix-and-ci-audit.md",
        rf"\b{tc['total_workspace']}\b",
        f"test-matrix-and-ci-audit.md cites {tc['total_workspace']} workspace tests",
        errors,
    )
    require_regex(
        repo_root,
        "docs/testing/test-matrix-and-ci-audit.md",
        rf"{tc['harness_self_tests']} specialized scripts(?:\*\*)? \(25 shell tests and 11 Python unit test suites\)",
        f"test-matrix-and-ci-audit.md cites {tc['harness_self_tests']} specialized scripts (25 shell + 11 Python)",
        errors,
    )
    forbid_regex(
        repo_root,
        "docs/testing/test-matrix-and-ci-audit.md",
        r"a262f19b0a2eb8a8|\b233\b|\b192\b|\b277\b|\b289\b|\b293\b|\b294\b|\b296\b|\b243\b|\b263\b|\b269\b|\b159\b|\b139\b|3\.84×|3\.84x|1\.34 GiB|25 diagnostic self-tests|25 scripts|25 specialized scripts|ALL 25 HARNESS|--max-ref-cache-mb|4 continuous performance budgets|Four-Target Workload Budget|harness/lib/oracle_platform\.sh|harness/bench_profile\.sh|test_sync_sclip_after_conseq_preserves_unmatched_entries_for_later_sv_checks|test_get_reference_uses_config_extension_not_hardcoded_refext",
        "test-matrix-and-ci-audit.md stale fingerprint, test counts, scale metrics, or 4-target budget",
        errors,
    )

    # 16. docs/validation/omnibus-empirical-validation.md
    require_regex(
        repo_root,
        "docs/validation/omnibus-empirical-validation.md",
        rf"v{ver_esc}",
        f"omnibus-empirical-validation.md cites v{version}",
        errors,
    )

    # 17. mkdocs.yml, wiki/Home.md, wiki/_Sidebar.md
    require_regex(
        repo_root,
        "mkdocs.yml",
        rf"Release v{ver_esc} Certification:\s*docs/release/v{ver_esc}-pr\.md",
        f"mkdocs.yml nav includes Release v{version} Certification",
        errors,
    )
    require_regex(
        repo_root,
        "wiki/Home.md",
        rf"v{ver_esc}",
        f"wiki/Home.md cites v{version}",
        errors,
    )
    require_regex(
        repo_root,
        "wiki/Home.md",
        rf"\b{tc['total_workspace']}\b",
        f"wiki/Home.md cites {tc['total_workspace']} workspace tests",
        errors,
    )
    forbid_regex(
        repo_root,
        "wiki/Home.md",
        r"\b233\b|\b277\b|\b289\b|\b293\b|\b294\b|\b296\b",
        "wiki/Home.md stale 233/277/289/293/294/296 test count",
        errors,
    )
    require_regex(
        repo_root,
        "wiki/_Sidebar.md",
        rf"docs/release/v{ver_esc}-pr\.md",
        f"wiki/_Sidebar.md links to docs/release/v{version}-pr.md",
        errors,
    )

    # 18. CLI 89-flag / 7-group & Migration Guide surfaces
    check_cli_and_migration_surfaces(repo_root, errors)

    # 19. UPSTREAM.toml commit & comparator pins
    require_regex(
        repo_root,
        "UPSTREAM.toml",
        r'^\[vardictjava\](?:\n(?!\[)[^\n]*)*?\nversion\s*=\s*"1\.8\.3"(?:\n(?!\[)[^\n]*)*?\ncommit\s*=\s*"0886261"',
        "UPSTREAM.toml pins VarDictJava 1.8.3 commit 0886261",
        errors,
    )
    require_regex(
        repo_root,
        "UPSTREAM.toml",
        r'^\[vardict_scripts\](?:\n(?!\[)[^\n]*)*?\ncommit\s*=\s*"009e017"',
        "UPSTREAM.toml pins VarDict scripts commit 009e017",
        errors,
    )
    require_regex(
        repo_root,
        "UPSTREAM.toml",
        r'^\[advisory_comparators\.vardictcpp\](?:\n(?!\[)[^\n]*)*?\nversion\s*=\s*"3"(?:\n(?!\[)[^\n]*)*?\nbioconda_build\s*=\s*"hc52dbad_0"',
        "UPSTREAM.toml pins advisory vardictcpp version 3 build hc52dbad_0",
        errors,
    )

    # 20. .github/workflows/ci.yml test-count synchronization (:71 and :175-177)
    int_suite_total = (
        tc["engine_behavior"]
        + tc["amplicon_stream_mem"]
        + tc["cram_stack"]
        + tc["cli_integration"]
    )
    require_regex(
        repo_root,
        ".github/workflows/ci.yml",
        rf"\({tc['core_lib']} vardictrust-core --lib \+ {tc['cli_bin']} vardictrust --bins\)",
        f"ci.yml Step 3 log cites ({tc['core_lib']} vardictrust-core --lib + {tc['cli_bin']} vardictrust --bins)",
        errors,
    )
    require_regex(
        repo_root,
        ".github/workflows/ci.yml",
        rf"Running {tc['corpus_parity']} in-process golden corpus and modality parity tests",
        f"ci.yml Step 5 log cites {tc['corpus_parity']} corpus_parity cases",
        errors,
    )
    require_regex(
        repo_root,
        ".github/workflows/ci.yml",
        rf"Unit Tests \(\\`{tc['core_lib']} vardictrust-core\\` \+ \\`{tc['cli_bin']} vardictrust\\`\)",
        f"ci.yml Step 9 summary cites (`{tc['core_lib']} vardictrust-core` + `{tc['cli_bin']} vardictrust`)",
        errors,
    )
    require_regex(
        repo_root,
        ".github/workflows/ci.yml",
        rf"Integration & Hostile-Input Tests \(\\`{int_suite_total}\\` tests\)",
        f"ci.yml Step 9 summary cites (`{int_suite_total}` tests)",
        errors,
    )
    require_regex(
        repo_root,
        ".github/workflows/ci.yml",
        rf"Corpus & Modality Parity Tests \(\\`{tc['corpus_parity']}\\` cases\)",
        f"ci.yml Step 9 summary cites (`{tc['corpus_parity']}` cases)",
        errors,
    )


def check_evidence_and_fingerprints(
    repo_root: Path, require_head_fp: bool, errors: list[str]
) -> None:
    fp_script = repo_root / "harness" / "engine_fingerprint.sh"
    if not fp_script.is_file():
        errors.append("MISSING FILE: harness/engine_fingerprint.sh")
        return

    ev_dir = repo_root / "benchmarks" / "evidence"
    json_files = sorted(ev_dir.glob("*.json"))
    if len(json_files) != 16:
        errors.append(
            f"benchmarks/evidence/*.json count == {len(json_files)}, expected 16"
        )
        return

    hg008_data = json.loads((ev_dir / "hg008-tier2.json").read_text(encoding="utf-8"))
    baseline_fp = hg008_data.get("fingerprint", "")
    scale_data = json.loads((ev_dir / "scale.json").read_text(encoding="utf-8"))
    baseline_perf_fp = scale_data.get("extra", {}).get("perf_fingerprint", "")

    for jf in json_files:
        data = json.loads(jf.read_text(encoding="utf-8"))
        if jf.name == "scale.json":
            scale_fp = data.get("fingerprint")
            scale_perf = data.get("extra", {}).get("perf_fingerprint")
            if scale_fp != baseline_fp and scale_perf != baseline_perf_fp:
                errors.append(
                    f"{jf.name}: neither fingerprint {scale_fp!r} nor extra.perf_fingerprint {scale_perf!r} matches certified baseline ({baseline_fp!r} / {baseline_perf_fp!r})"
                )
        else:
            if data.get("fingerprint") != baseline_fp:
                errors.append(
                    f"{jf.name}: fingerprint {data.get('fingerprint')!r} != baseline {baseline_fp!r}"
                )

    # When cutting/certifying a release (--require-head-fingerprint), verify the
    # committed evidence matches the live git-index engine_fingerprint.sh.
    if require_head_fp and (repo_root / ".git").exists():
        res_fp = subprocess.run(
            ["bash", str(fp_script)],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
        res_perf = subprocess.run(
            ["bash", str(fp_script), "--perf"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
        if res_fp.returncode != 0:
            errors.append(
                f"harness/engine_fingerprint.sh failed (exit {res_fp.returncode}): {res_fp.stderr.strip()}"
            )
        if res_perf.returncode != 0:
            errors.append(
                f"harness/engine_fingerprint.sh --perf failed (exit {res_perf.returncode}): {res_perf.stderr.strip()}"
            )
        if res_fp.returncode == 0 and res_perf.returncode == 0:
            head_fp = res_fp.stdout.strip()
            head_perf_fp = res_perf.stdout.strip()
            if baseline_fp != head_fp:
                errors.append(
                    f"benchmarks/evidence/*.json baseline fingerprint {baseline_fp!r} != HEAD fingerprint {head_fp!r} (run bash harness/recertify.sh)"
                )
            scale_fp = scale_data.get("fingerprint")
            scale_perf = scale_data.get("extra", {}).get("perf_fingerprint")
            if scale_fp != head_fp and scale_perf != head_perf_fp:
                errors.append(
                    f"scale.json: neither fingerprint {scale_fp!r} != {head_fp!r} nor extra.perf_fingerprint {scale_perf!r} != {head_perf_fp!r}"
                )

    # Verify omnibus-empirical-validation.md §13 SHA-256 table matches on-disk JSONs
    omnibus_text = read_text(
        repo_root, "docs/validation/omnibus-empirical-validation.md", errors
    )
    if omnibus_text:
        pattern = re.compile(
            r"\|\s*\[`([a-z0-9_-]+\.json)`\]\([^)]+\)\s*\|[^|]+\|\s*`([0-9a-f]{64})`\s*\|"
        )
        matches = pattern.findall(omnibus_text)
        if len(matches) != 16:
            errors.append(
                f"omnibus-empirical-validation.md §13 SHA-256 table has {len(matches)} rows, expected 16"
            )
        for fname, expected_sha in matches:
            fpath = repo_root / "benchmarks" / "evidence" / fname
            if not fpath.is_file():
                errors.append(f"MISSING EVIDENCE FILE: benchmarks/evidence/{fname}")
                continue
            actual_sha = hashlib.sha256(fpath.read_bytes()).hexdigest()
            if actual_sha != expected_sha:
                errors.append(
                    f"benchmarks/evidence/{fname}: sha256 {actual_sha} != omnibus table {expected_sha}"
                )


def check_downstream_phases(
    repo_root: Path, version: str, phase: str, errors: list[str]
) -> None:
    ver_esc = re.escape(version)
    if phase in ("post-bioconda", "post-nfcore"):
        require_regex(
            repo_root,
            "packaging/bioconda/recipes/vardictrust/meta.yaml",
            rf'\{{%\s*set\s+version\s*=\s*"{ver_esc}"\s*%\}}',
            f"Bioconda meta.yaml version == {version}",
            errors,
        )
    if phase == "post-nfcore":
        require_regex(
            repo_root,
            "packaging/nf-core/modules/vardictrust/environment.yml",
            rf"bioconda::vardictrust={ver_esc}\b",
            f"nf-core environment.yml version == {version}",
            errors,
        )
        require_regex(
            repo_root,
            "packaging/nf-core/modules/vardictrust/tests/main.nf.test.snap",
            rf'"VARDICTRUST",\s*"vardictrust",\s*"{ver_esc}"',
            f"nf-core main.nf.test.snap version == {version}",
            errors,
            flags=re.MULTILINE | re.DOTALL,
        )


REGISTERED_INTEGRATION_TESTS: dict[str, tuple[str, ...]] = {
    "vardictrust-core": (
        "corpus_parity",
        "engine_behavior",
        "amplicon_stream_mem",
        "cram_stack",
    ),
    "vardictrust": ("profile_smoke", "robustness", "version"),
}


def check_registered_integration_tests(repo_root: Path, errors: list[str]) -> bool:
    ok = True
    for pkg, registered_tuple in REGISTERED_INTEGRATION_TESTS.items():
        registered = set(registered_tuple)
        tests_dir = repo_root / "crates" / pkg / "tests"
        on_disk = (
            {p.stem for p in tests_dir.glob("*.rs")} if tests_dir.is_dir() else set()
        )
        unregistered = sorted(on_disk - registered)
        missing = sorted(registered - on_disk)
        if unregistered:
            errors.append(
                f"crates/{pkg}/tests has unregistered integration test target(s): "
                f"{unregistered} (registered: {sorted(registered)})"
            )
            ok = False
        if missing:
            errors.append(
                f"crates/{pkg}/tests is missing registered integration test target(s): "
                f"{missing} (registered: {sorted(registered)})"
            )
            ok = False
    return ok


def _count_cargo_target_tests(
    repo_root: Path,
    cargo_args: list[str],
    errors: list[str],
    *,
    use_llvm_cov: bool = False,
) -> int | None:
    cmd = (
        ["cargo", "llvm-cov", "--no-report", *cargo_args, "--", "--list"]
        if use_llvm_cov
        else ["cargo", "test", *cargo_args, "--", "--list"]
    )
    cmd_label = " ".join(cmd)
    try:
        res = subprocess.run(
            cmd,
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        errors.append(f"{cmd_label} could not execute 'cargo': {exc}")
        return None
    if res.returncode != 0:
        errors.append(
            f"{cmd_label} failed (exit {res.returncode}): {res.stderr.strip()}"
        )
        return None
    return sum(
        1 for line in res.stdout.splitlines() if re.search(r": test$", line) is not None
    )


def check_cargo_test_counts(
    repo_root: Path, errors: list[str], *, use_llvm_cov: bool = False
) -> None:
    all_ok = check_registered_integration_tests(repo_root, errors)
    single_targets = (
        ("core_lib", ["-p", "vardictrust-core", "--lib"]),
        *(
            (tgt, ["-p", "vardictrust-core", "--test", tgt])
            for tgt in REGISTERED_INTEGRATION_TESTS["vardictrust-core"]
        ),
        ("cli_bin", ["-p", "vardictrust", "--bin", "vardictrust"]),
    )
    cli_integration_targets = REGISTERED_INTEGRATION_TESTS["vardictrust"]

    observed: dict[str, int] = {}
    for key, args in single_targets:
        cnt = _count_cargo_target_tests(
            repo_root, args, errors, use_llvm_cov=use_llvm_cov
        )
        if cnt is None:
            all_ok = False
        else:
            observed[key] = cnt

    cli_int_sum = 0
    cli_int_ok = True
    for tgt in cli_integration_targets:
        cnt = _count_cargo_target_tests(
            repo_root,
            ["-p", "vardictrust", "--test", tgt],
            errors,
            use_llvm_cov=use_llvm_cov,
        )
        if cnt is None:
            cli_int_ok = False
            all_ok = False
        else:
            cli_int_sum += cnt
    if cli_int_ok:
        observed["cli_integration"] = cli_int_sum

    if all_ok:
        observed["total_workspace"] = sum(observed.values())

    for key in (
        "core_lib",
        "corpus_parity",
        "engine_behavior",
        "amplicon_stream_mem",
        "cram_stack",
        "cli_bin",
        "cli_integration",
        "total_workspace",
    ):
        if key in observed and observed[key] != EXPECTED_TEST_COUNTS[key]:
            errors.append(
                f"Cargo test count mismatch for {key}: "
                f"expected {EXPECTED_TEST_COUNTS[key]}, "
                f"cargo --list reported {observed[key]}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate VarDictRust release version and evidence consistency."
    )
    parser.add_argument(
        "--version",
        default=None,
        help="Target SemVer (e.g. 1.0.7); defaults to [workspace.package].version in Cargo.toml",
    )
    parser.add_argument(
        "--phase",
        choices=VALID_PHASES,
        default="pre-tag",
        help="Release lifecycle phase to validate",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[4],
        help="Path to VarDictRust repository root",
    )
    parser.add_argument(
        "--require-head-fingerprint",
        action="store_true",
        help="Require benchmarks/evidence/*.json fingerprints to match live HEAD engine_fingerprint.sh (used during release certification)",
    )
    parser.add_argument(
        "--verify-cargo-list",
        action="store_true",
        help="Run cargo test -- --list across workspace library, binary, and integration targets (excluding benches and doctests) to verify EXPECTED_TEST_COUNTS",
    )
    parser.add_argument(
        "--cargo-llvm-cov",
        action="store_true",
        help="Invoke cargo llvm-cov --no-report ... -- --list with --verify-cargo-list to reuse warm instrumented CI binaries",
    )
    args = parser.parse_args()

    version = args.version
    if version is None:
        cargo_toml = (args.repo_root / "Cargo.toml").read_text(encoding="utf-8")
        ver_m = re.search(
            r'^\[workspace\.package\]\s*\nversion\s*=\s*"([^"]+)"',
            cargo_toml,
            re.MULTILINE,
        )
        if not ver_m:
            print("ERROR: could not read [workspace.package].version from Cargo.toml", file=sys.stderr)
            return 2
        version = ver_m.group(1)
    elif version.startswith("v") and re.match(r"^v\d+\.\d+\.\d+$", version):
        version = version[1:]

    if not re.match(r"^\d+\.\d+\.\d+$", version):
        print(f"ERROR: invalid SemVer version {version!r}", file=sys.stderr)
        return 2

    errors: list[str] = []
    check_pre_tag_surfaces(
        args.repo_root,
        version,
        errors,
        require_head_fp=args.require_head_fingerprint,
    )
    check_evidence_and_fingerprints(
        args.repo_root, args.require_head_fingerprint, errors
    )
    check_downstream_phases(args.repo_root, version, args.phase, errors)
    if args.verify_cargo_list:
        check_cargo_test_counts(
            args.repo_root, errors, use_llvm_cov=args.cargo_llvm_cov
        )
    else:
        check_registered_integration_tests(args.repo_root, errors)

    if errors:
        print(
            f"RELEASE CONSISTENCY FAILED (version={version}, phase={args.phase}, {len(errors)} error(s)):",
            file=sys.stderr,
        )
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    print(
        f"RELEASE CONSISTENCY OK: version={version}, phase={args.phase} (all surfaces and 16 evidence manifests verified)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
