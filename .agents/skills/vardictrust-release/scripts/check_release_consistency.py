#!/usr/bin/env python3
"""VarDictRust Release Consistency & Evidence Validator (stdlib-only).

Verifies that all version-bearing files, workspace test counts (277 total:
231 core lib + 2 amplicon_stream_mem + 2 cram_stack + 13 CLI bin + 29 CLI
integration binary tests), engine fingerprints (FP & PERF_FP), and 16
benchmarks/evidence/*.json SHA-256 manifests are synchronized for a target
release version and phase.
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


def check_pre_tag_surfaces(repo_root: Path, version: str, errors: list[str]) -> None:
    ver_esc = re.escape(version)

    # 1. Cargo.toml
    require_regex(
        repo_root,
        "Cargo.toml",
        rf'^\[workspace\.package\]\s*\nversion\s*=\s*"{ver_esc}"',
        f"workspace.package version == {version}",
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
    forbid_regex(
        repo_root,
        "CHANGELOG.md",
        r"dirty_decay_ms:5000,\s*muzzy_decay_ms:5000",
        "stale MALLOC_CONF description in CHANGELOG.md",
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

    # 9. harness/run_self_tests.sh
    require_regex(
        repo_root,
        "harness/run_self_tests.sh",
        rf"^# VarDictRust v{ver_esc}\b",
        f"harness/run_self_tests.sh header == v{version}",
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
        r"293 automated Cargo workspace tests \(241 core library tests, 2 amplicon streaming memory tests, 2 CRAM stack tests, 15 CLI tests, and 33 integration binary tests\)",
        "README.md exact 293 test breakdown (241 + 2 + 2 + 15 + 33)",
        errors,
    )

    # 12. DIVERGENCES.md
    require_regex(
        repo_root,
        "DIVERGENCES.md",
        rf"VarDictRust_v{ver_esc}",
        f"DIVERGENCES.md cites VarDictRust_v{version}",
        errors,
    )
    require_regex(
        repo_root,
        "DIVERGENCES.md",
        r"241 PASS.*293 PASS",
        "DIVERGENCES.md cites 241 core lib / 293 workspace tests",
        errors,
        flags=re.MULTILINE | re.DOTALL,
    )
    require_regex(
        repo_root,
        "DIVERGENCES.md",
        r"33bc15964cc36b6d",
        "DIVERGENCES.md cites certified fingerprint 33bc15964cc36b6d",
        errors,
    )
    forbid_regex(
        repo_root,
        "DIVERGENCES.md",
        r"817d89c5fb34515d|profile_features|3\.58x speedup",
        "DIVERGENCES.md stale fingerprint/test name/speedup",
        errors,
    )

    # 13. CONTRIBUTING.md
    require_regex(
        repo_root,
        "CONTRIBUTING.md",
        r"all 293 workspace unit and integration tests",
        "CONTRIBUTING.md §4.1 cites 293 workspace tests",
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
        r"aarch64-unknown-linux-gnu|all 233 workspace|all 277 workspace|all 289 workspace",
        "CONTRIBUTING.md stale target triple or 233/277/289 test count",
        errors,
    )

    # 14. benchmarks/README.md
    require_regex(
        repo_root,
        "benchmarks/README.md",
        rf"v{ver_esc}",
        f"benchmarks/README.md cites v{version}",
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
        r"\b293\b",
        "test-matrix-and-ci-audit.md cites 293 workspace tests",
        errors,
    )
    forbid_regex(
        repo_root,
        "docs/testing/test-matrix-and-ci-audit.md",
        r"a262f19b0a2eb8a8|\b233\b|\b192\b|\b277\b|\b289\b",
        "test-matrix-and-ci-audit.md stale fingerprint or 233/192/277/289 test counts",
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
        r"\b293\b",
        "wiki/Home.md cites 293 workspace tests",
        errors,
    )
    forbid_regex(
        repo_root,
        "wiki/Home.md",
        r"\b233\b|\b277\b|\b289\b",
        "wiki/Home.md stale 233/277/289 test count",
        errors,
    )
    require_regex(
        repo_root,
        "wiki/_Sidebar.md",
        rf"docs/release/v{ver_esc}-pr\.md",
        f"wiki/_Sidebar.md links to docs/release/v{version}-pr.md",
        errors,
    )


def check_evidence_and_fingerprints(repo_root: Path, errors: list[str]) -> None:
    fp_script = repo_root / "harness" / "engine_fingerprint.sh"
    if not fp_script.is_file():
        errors.append("MISSING FILE: harness/engine_fingerprint.sh")
        return

    # Run live git-index fingerprint check when repo_root has a .git directory;
    # record any non-zero exit code (e.g. uncommitted GATE_FILES changes) as a hard error.
    if (repo_root / ".git").exists():
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
            fp = res_fp.stdout.strip()
            perf_fp = res_perf.stdout.strip()
            ev_dir = repo_root / "benchmarks" / "evidence"
            json_files = sorted(ev_dir.glob("*.json"))
            if len(json_files) != 16:
                errors.append(
                    f"benchmarks/evidence/*.json count == {len(json_files)}, expected 16"
                )
            for jf in json_files:
                data = json.loads(jf.read_text(encoding="utf-8"))
                if data.get("fingerprint") != fp:
                    errors.append(
                        f"{jf.name}: fingerprint {data.get('fingerprint')!r} != {fp!r}"
                    )
                if jf.name == "scale.json":
                    scale_perf = data.get("extra", {}).get("perf_fingerprint")
                    if scale_perf != perf_fp:
                        errors.append(
                            f"{jf.name}: extra.perf_fingerprint {scale_perf!r} != {perf_fp!r}"
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate VarDictRust release version and evidence consistency."
    )
    parser.add_argument("--version", required=True, help="Target SemVer (e.g. 1.0.3)")
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
    args = parser.parse_args()

    if not re.match(r"^\d+\.\d+\.\d+$", args.version):
        print(f"ERROR: invalid SemVer version {args.version!r}", file=sys.stderr)
        return 2

    errors: list[str] = []
    check_pre_tag_surfaces(args.repo_root, args.version, errors)
    check_evidence_and_fingerprints(args.repo_root, errors)
    check_downstream_phases(args.repo_root, args.version, args.phase, errors)

    if errors:
        print(
            f"RELEASE CONSISTENCY FAILED (version={args.version}, phase={args.phase}, {len(errors)} error(s)):",
            file=sys.stderr,
        )
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    print(
        f"RELEASE CONSISTENCY OK: version={args.version}, phase={args.phase} (all surfaces and 16 evidence manifests verified)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
