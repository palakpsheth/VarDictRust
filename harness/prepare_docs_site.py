#!/usr/bin/env python3
"""
harness/prepare_docs_site.py — Stage repository markdown files for MkDocs site build.

Discovers all tracked markdown files in the repository and mirrors their directory
hierarchy inside `.docs_staged/`. By preserving the relative file paths, all
intra-document relative links and anchor hashes resolve identically on GitHub and ReadTheDocs.
"""

import pathlib
import shutil
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
STAGING_DIR = REPO_ROOT / ".docs_staged"

IGNORED_PREFIXES = (
    ".git/",
    ".github/",
    ".superpowers/",
    ".agents/",
    "target/",
    "scratch/",
    "VarDictJava/",
    "node_modules/",
    ".docs_staged/",
    "site/",
)


def stage_docs():
    """Mirror repository markdown files into .docs_staged."""
    if STAGING_DIR.exists():
        shutil.rmtree(STAGING_DIR)
    STAGING_DIR.mkdir(parents=True)

    # 1. Discover all tracked markdown documents outside ignored folders
    md_files = [
        p
        for p in REPO_ROOT.rglob("*.md")
        if not any(p.relative_to(REPO_ROOT).as_posix().startswith(ign) for ign in IGNORED_PREFIXES)
    ]

    for p in md_files:
        rel = p.relative_to(REPO_ROOT)
        target = STAGING_DIR / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)

    # 2. Stage .github/ documents into community/ (MkDocs ignores dot-directories by default)
    github_dir = REPO_ROOT / ".github"
    if github_dir.exists():
        community_dir = STAGING_DIR / "community"
        community_dir.mkdir(parents=True, exist_ok=True)
        for gh_md in github_dir.glob("*.md"):
            shutil.copy2(gh_md, community_dir / gh_md.name)

    print(f"prepare_docs_site: Staged {len(md_files)} markdown files into {STAGING_DIR.name}/")


if __name__ == "__main__":
    stage_docs()
