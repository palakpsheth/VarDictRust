#!/usr/bin/env python3
"""
harness/prepare_docs_site.py — Stage repository markdown files for MkDocs site build.

Discovers all tracked markdown files in the repository and mirrors their directory
hierarchy inside `.docs_staged/`. By preserving the relative file paths, all
intra-document relative links and anchor hashes resolve identically on GitHub and ReadTheDocs.

Additionally performs link normalization:
- Rewrites links pointing to non-markdown source/evidence files (.rs, .sh, .py, .json, .yml, LICENSE, etc.)
  to their canonical GitHub blob URLs (preventing 404s on static doc site).
- Rewrites links pointing to .github/*.md into their staged community/*.md equivalents.
- Converts raw HTML anchor tags (<a id="..."></a>) into attribute-list elements so MkDocs
  strict anchor validation can verify them.
"""

import os
import pathlib
import re
import shutil
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
STAGING_DIR = REPO_ROOT / ".docs_staged"
GITHUB_BLOB_BASE = "https://github.com/palakpsheth/VarDictRust/blob/main"

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
    ".worktrees/",
    "wiki/",
)

# Regex matching Markdown links [text](url) and badge links [![alt](img)](url)
LINK_RE = re.compile(
    r'(\[(?:!\[.*?\]\(.*?\)|[^\]])*\]\()([^\)\s]+)((?:\s+["\x27][^"\x27]*["\x27])?\))'
)

# Regex matching raw HTML anchors <a id="..." ...></a> or <a name="..." ...></a>
HTML_ANCHOR_RE = re.compile(
    r'<a\s+(?:[^>]*?\s+)?(?:id|name)=["\x27]([a-zA-Z0-9_\-]+)["\x27][^>]*>(?:</a>)?'
)


def rewrite_content(content: str, staged_file: pathlib.Path) -> str:
    """Rewrite links and anchors in a staged markdown file to prevent 404s."""
    rel_staged = staged_file.relative_to(STAGING_DIR)

    # Determine original source location in REPO_ROOT
    if rel_staged.parts[0] == "community":
        source_dir = REPO_ROOT / ".github"
    else:
        source_dir = (REPO_ROOT / rel_staged).parent

    def repl(match):
        prefix = match.group(1)
        target = match.group(2)
        suffix = match.group(3)

        # Do not rewrite web URLs, mailto, or in-page anchors
        if target.startswith(("http://", "https://", "mailto:", "ftp://", "#")):
            return match.group(0)

        url_path, sep, anchor = target.partition("#")
        if not url_path:
            return match.group(0)

        # Resolve target relative to original source file directory
        target_resolved = (source_dir / url_path).resolve()

        try:
            rel_to_repo = target_resolved.relative_to(REPO_ROOT).as_posix()
        except ValueError:
            # Outside repo root
            return match.group(0)

        # 1. Target is a .github/ markdown file
        if rel_to_repo.startswith(".github/") and rel_to_repo.endswith(".md"):
            dest_staged = STAGING_DIR / "community" / target_resolved.name
            new_rel = os.path.relpath(dest_staged, staged_file.parent)
            new_target = f"{new_rel}#{anchor}" if anchor else new_rel
            return f"{prefix}{new_target}{suffix}"

        # 2. Target is a non-markdown file or a directory
        if not rel_to_repo.endswith(".md") or target_resolved.is_dir():
            new_target = f"{GITHUB_BLOB_BASE}/{rel_to_repo}" + (f"#{anchor}" if anchor else "")
            return f"{prefix}{new_target}{suffix}"

        # 3. Target is a markdown file within the repository
        dest_staged = STAGING_DIR / rel_to_repo
        if dest_staged.exists() or (REPO_ROOT / rel_to_repo).exists():
            new_rel = os.path.relpath(dest_staged, staged_file.parent)
            new_target = f"{new_rel}#{anchor}" if anchor else new_rel
            return f"{prefix}{new_target}{suffix}"

        return match.group(0)

    # Step A: Rewrite markdown links
    content = LINK_RE.sub(repl, content)

    # Step B: Convert raw HTML anchors into Markdown attr_list elements for MkDocs validation
    content = HTML_ANCHOR_RE.sub(r"[](){: #\1 }", content)

    return content


def stage_docs():
    """Mirror repository markdown files into .docs_staged with link rewriting."""
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
    staged_gh_files = []
    if github_dir.exists():
        community_dir = STAGING_DIR / "community"
        community_dir.mkdir(parents=True, exist_ok=True)
        for gh_md in github_dir.glob("*.md"):
            dest = community_dir / gh_md.name
            shutil.copy2(gh_md, dest)
            staged_gh_files.append(dest)

    # 3. Rewrite relative links to non-md assets and .github docs across all staged files
    all_staged = list(STAGING_DIR.rglob("*.md"))
    for staged_md in all_staged:
        old_text = staged_md.read_text(encoding="utf-8")
        new_text = rewrite_content(old_text, staged_md)
        if new_text != old_text:
            staged_md.write_text(new_text, encoding="utf-8")

    print(
        f"prepare_docs_site: Staged {len(md_files) + len(staged_gh_files)} markdown files "
        f"into {STAGING_DIR.name}/ and normalized cross-references."
    )


if __name__ == "__main__":
    stage_docs()
