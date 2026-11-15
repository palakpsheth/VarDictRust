#!/usr/bin/env python3
"""
harness/tests/test_docs_ui.py

Automated UI and document integrity test suite for VarDictRust documentation.

1. Static Markdown Structure Validation (pure stdlib):
   - Table column consistency across header, delimiter, and rows.
   - Tag hygiene: non-HTML angle brackets (e.g. `<dup-N>`) must be backtick-escaped.
   - Markdown link format and anchor validity.

2. Browser UI Rendering with Playwright (dynamic import):
   - Loads GitHub-Flavored Markdown rendering of key documents (ROADMAP.md, README.md, etc.)
     into headless Chromium.
   - Verifies table DOM structure, row/column rendering, and absence of broken HTML nodes.
   - Verifies external link structures.
"""

import functools
import hashlib
import html
import http.server
import importlib
import json
import os
import pathlib
import re
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.parse
import urllib.request

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent

# Standard permitted HTML tags in Markdown
ALLOWED_HTML_TAGS = frozenset([
    "a", "abbr", "b", "blockquote", "br", "code", "dd", "del", "details",
    "div", "dl", "dt", "em", "h1", "h2", "h3", "h4", "h5", "h6", "hr",
    "i", "img", "ins", "kbd", "li", "ol", "p", "pre", "q", "rp", "rt",
    "ruby", "s", "samp", "small", "span", "strong", "sub", "summary",
    "sup", "table", "tbody", "td", "tfoot", "th", "thead", "tr", "ul",
    "var", "wbr"
])


def find_all_markdown_files():
    """Discover all tracked markdown files in the repository."""
    md_files = []
    for path in REPO_ROOT.rglob("*.md"):
        # Exclude git internals, cargo target, conda/python virtualenvs, scratch, and transient dev directories
        rel = str(path.relative_to(REPO_ROOT))
        if any(rel.startswith(ignore) for ignore in [
            ".git/", ".superpowers/", ".docs_staged/", "site/",
            "target/", "scratch/", "VarDictJava/", "node_modules/",
            ".venv/", "venv/", ".pytest_cache/", ".worktrees/",
            "benchmarks/data/", ".venv_docs/", ".cache/",
        ]):
            continue
        md_files.append(path)
    return sorted(md_files)


class TestMarkdownStaticHygiene(unittest.TestCase):
    """Test static formatting, table alignment, and escaping across all markdown documents."""

    def test_zensical_config_not_in_nav_and_validation(self):
        """Verify zensical.toml top-level not_in_nav, all 27 nav entries, and link validation."""
        toml_path = REPO_ROOT / "zensical.toml"
        self.assertTrue(toml_path.is_file(), "zensical.toml must exist at repository root")
        zensical_text = toml_path.read_text(encoding="utf-8")

        # 1. Verify not_in_nav multiline string
        self.assertRegex(
            zensical_text,
            re.compile(r'not_in_nav\s*=\s*"""', re.MULTILINE),
            "zensical.toml must declare 'not_in_nav = \"\"\"' for partial/historical docs",
        )
        for expected_excluded in (
            "benchmarks/gcp-prereg.md",
            "docs/validation/gcp-study-protocol.md",
            "docs/validation/scale-recalibration-proposal.md",
            "docs/validation/provenance/task-*.md",
            "paper/tables/*.md",
        ):
            self.assertIn(
                expected_excluded,
                zensical_text,
                f"zensical.toml not_in_nav must include {expected_excluded}",
            )

        # 2. Extract and verify all 27 nav entries
        nav_entries = re.findall(
            r'\{\s*"([^"]+)"\s*=\s*"([^"]+)"\s*\}',
            zensical_text,
        )
        self.assertGreaterEqual(
            len(nav_entries), 27,
            f"zensical.toml must define at least 27 navigation pages; found {len(nav_entries)}: {nav_entries}",
        )
        nav_dict = dict(nav_entries)
        self.assertEqual(nav_dict.get("Overview"), "README.md")
        self.assertEqual(nav_dict.get("Migration Guide"), "docs/migration.md")
        self.assertEqual(nav_dict.get("JVM & R Kernel Oracles"), "fixtures/jvm-oracles/README.md")
        self.assertEqual(nav_dict.get("Release v1.0.8 Certification"), "docs/release/v1.0.8-pr.md")

        # 3. Verify validation settings
        self.assertRegex(
            zensical_text,
            re.compile(r"\[project\.validation\]\s*\n(?:[^\n]*\n)*?invalid_links\s*=\s*true", re.MULTILINE),
            "zensical.toml must enforce invalid_links = true in [project.validation]",
        )
        self.assertRegex(
            zensical_text,
            re.compile(r"\[project\.validation\]\s*\n(?:[^\n]*\n)*?invalid_link_anchors\s*=\s*true", re.MULTILINE),
            "zensical.toml must enforce invalid_link_anchors = true in [project.validation]",
        )

        # 4. Verify prepare_docs_site.py staging exclusions
        prep_text = (REPO_ROOT / "harness" / "prepare_docs_site.py").read_text(encoding="utf-8")
        for excluded in ('".venv/"', '"venv/"', '".pytest_cache/"', '".cache/"', '".venv_docs/"'):
            self.assertIn(
                excluded,
                prep_text,
                f"prepare_docs_site.py IGNORED_PREFIXES must exclude {excluded}",
            )

    @staticmethod
    def _is_valid_zensical_pin(reqs_text):
        m = re.search(r"^zensical\s*([><=!~,\d.\s]+?)\s*(?:#.*)?$", reqs_text, re.MULTILINE)
        if not m:
            return False
        clauses = [c.strip() for c in m.group(1).split(",") if c.strip()]
        if len(clauses) != 2:
            return False
        has_lower = False
        has_upper = False
        for clause in clauses:
            low_m = re.fullmatch(r">=\s*0\.(\d+)(?:\.(\d+))?", clause)
            if low_m:
                minor = int(low_m.group(1))
                patch = int(low_m.group(2) or 0)
                if (0, minor, patch) >= (0, 1, 1):
                    has_lower = True
                continue
            if re.fullmatch(r"<\s*0\.2(?:\.0)?", clause):
                has_upper = True
        return has_lower and has_upper

    def test_zensical_requirements_and_dependabot_config(self):
        """Verify docs/requirements.txt pins zensical>=0.1.1,<0.2.0 and dependabot.yml ignores zensical >= 0.2.0."""
        reqs_text = (REPO_ROOT / "docs" / "requirements.txt").read_text(encoding="utf-8")
        self.assertTrue(
            self._is_valid_zensical_pin(reqs_text),
            "docs/requirements.txt must pin zensical>=0.1.1,<0.2.0",
        )
        self.assertNotIn("mkdocs", reqs_text, "docs/requirements.txt must not contain mkdocs")
        self.assertNotIn("mkdocs-material", reqs_text, "docs/requirements.txt must not contain mkdocs-material")
        self.assertNotIn("pymdown-extensions", reqs_text, "docs/requirements.txt must not contain pymdown-extensions")

        for compatible in (
            "zensical>=0.1.1,<0.2.0\n",
            "zensical>=0.1.2,<0.2.0\n",
            "zensical>=0.1.5,<0.2\n",
            "zensical<0.2.0,>=0.1.1\n",
        ):
            self.assertTrue(
                self._is_valid_zensical_pin(compatible),
                f"Validator must accept compatible Zensical update {compatible.strip()!r}",
            )
        for invalid in (
            "zensical>=0.1.0,<0.2.0\n",
            "zensical>=0.1.1\n",
            "zensical>=0.1.1,<0.3.0\n",
            "zensical>=0.2.0,<0.3.0\n",
        ):
            self.assertFalse(
                self._is_valid_zensical_pin(invalid),
                f"Validator must reject out-of-range or unpinned spec {invalid.strip()!r}",
            )

        dependabot_text = (REPO_ROOT / ".github" / "dependabot.yml").read_text(encoding="utf-8")
        pip_m = re.search(
            r'^  - package-ecosystem:\s*"pip"\n(?:(?!^  - package-ecosystem:)[^\n]*\n?)*',
            dependabot_text,
            re.MULTILINE,
        )
        self.assertIsNotNone(pip_m, ".github/dependabot.yml missing pip ecosystem block")
        pip_block = pip_m.group(0)
        self.assertRegex(
            pip_block,
            re.compile(r'^    directory:\s*"/docs"\s*$', re.MULTILINE),
            "pip ecosystem block in .github/dependabot.yml must target directory '/docs'",
        )
        ignore_m = re.search(
            r'^    ignore:\s*\n'
            r'      - dependency-name:\s*"zensical"\s*\n'
            r'        versions:\s*\[">= 0\.2\.0"\]\s*$',
            pip_block,
            re.MULTILINE,
        )
        self.assertIsNotNone(
            ignore_m,
            '.github/dependabot.yml /docs pip entry must declare 4-space ignore block for zensical >= 0.2.0',
        )
        groups_m = re.search(r"^    groups:\s*$", pip_block, re.MULTILINE)
        self.assertIsNotNone(groups_m, ".github/dependabot.yml /docs pip entry missing 4-space groups block")
        self.assertLess(
            ignore_m.start(),
            groups_m.start(),
            "ignore block in .github/dependabot.yml /docs pip entry must not be nested under groups",
        )

    def test_readthedocs_yaml_zensical_commands(self):
        """Verify .readthedocs.yaml uses ReadTheDocs v2 build.commands contract for Zensical."""
        rtd_path = REPO_ROOT / ".readthedocs.yaml"
        self.assertTrue(rtd_path.is_file(), ".readthedocs.yaml must exist at repository root")
        rtd_text = rtd_path.read_text(encoding="utf-8")

        # 1. Version 2
        self.assertRegex(rtd_text, r"(?m)^version:\s*2\b", ".readthedocs.yaml must define version: 2")

        # 2. Uses build.commands contract
        self.assertIn("commands:", rtd_text, ".readthedocs.yaml must define build.commands")

        # 3. Explicit virtualenv creation at ${READTHEDOCS_VIRTUALENV_PATH:-.venv_rtd}
        self.assertIn(
            'python3 -m venv "${READTHEDOCS_VIRTUALENV_PATH:-.venv_rtd}"',
            rtd_text,
            ".readthedocs.yaml must create an explicit virtualenv at ${READTHEDOCS_VIRTUALENV_PATH:-.venv_rtd}",
        )
        self.assertIn(
            '-r docs/requirements.txt',
            rtd_text,
            ".readthedocs.yaml must install docs/requirements.txt",
        )

        # 4. prepare_docs_site.py and zensical build --clean --strict
        self.assertIn(
            "harness/prepare_docs_site.py",
            rtd_text,
            ".readthedocs.yaml must run harness/prepare_docs_site.py",
        )
        self.assertIn(
            'zensical" build --clean --strict',
            rtd_text,
            ".readthedocs.yaml must run zensical build --clean --strict",
        )

        # 5. Copy built site to $READTHEDOCS_OUTPUT/html/
        self.assertIn(
            'cp -a site/. "$READTHEDOCS_OUTPUT/html/"',
            rtd_text,
            ".readthedocs.yaml must copy site/. to $READTHEDOCS_OUTPUT/html/",
        )

        # 6. Must NOT contain mutually exclusive mkdocs, sphinx, jobs, or top-level python blocks
        for forbidden in ("mkdocs:", "sphinx:", "jobs:", "\npython:"):
            self.assertNotIn(
                forbidden,
                rtd_text,
                f".readthedocs.yaml using build.commands must not contain {forbidden.strip()}",
            )

        # 7. Check pure stdlib YAML list item quoting syntax
        for line in rtd_text.splitlines():
            stripped = line.strip()
            if stripped.startswith('- "'):
                after_first = stripped[3:]
                close_quote_idx = after_first.find('"')
                if close_quote_idx != -1:
                    trailing = after_first[close_quote_idx + 1:].strip()
                    self.assertTrue(
                        trailing == "" or trailing.startswith("#"),
                        f"Found unquoted trailing scalar after double-quoted string in .readthedocs.yaml: {line!r}",
                    )

        # 8. Dynamic PyYAML validation when PyYAML is installed on host
        try:
            yaml_mod = importlib.import_module("yaml")
            parsed = yaml_mod.safe_load(rtd_text)
            self.assertIsInstance(parsed, dict, ".readthedocs.yaml must parse as a dictionary")
            self.assertEqual(parsed.get("version"), 2)
            self.assertIn("build", parsed)
            self.assertIn("commands", parsed["build"])
            commands = parsed["build"]["commands"]
            self.assertIsInstance(commands, list)
            self.assertEqual(
                len(commands), 7,
                f"Expected exactly 7 commands in .readthedocs.yaml build.commands; got {len(commands)}: {commands}",
            )
        except ImportError:
            pass

    def test_table_column_balance(self):
        """Verify that every markdown table has consistent column counts across all rows."""
        md_files = find_all_markdown_files()
        self.assertGreater(len(md_files), 10, "Should discover at least 10 markdown documents.")

        table_errors = []
        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            with open(md_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            in_code_block = False
            current_table_lines = []
            current_table_start = 0

            for idx, line in enumerate(lines, 1):
                stripped = line.strip()
                if stripped.startswith("```"):
                    in_code_block = not in_code_block
                    continue
                if in_code_block:
                    continue

                if stripped.startswith("|") and stripped.endswith("|"):
                    if not current_table_lines:
                        current_table_start = idx
                    current_table_lines.append((idx, stripped))
                else:
                    if current_table_lines:
                        # Validate completed table
                        self._validate_table(rel_path, current_table_start, current_table_lines, table_errors)
                        current_table_lines = []

            if current_table_lines:
                self._validate_table(rel_path, current_table_start, current_table_lines, table_errors)

        self.assertEqual(
            table_errors, [],
            "Table column count mismatches found:\n" + "\n".join(f"  - {e}" for e in table_errors)
        )

    def _validate_table(self, rel_path, start_line, table_lines, errors):
        """Helper to check column count equality within a single table."""
        if len(table_lines) < 2:
            return

        expected_cols = None
        for line_num, line_str in table_lines:
            # Split by | excluding escaped \|
            cells = [c.strip() for c in re.split(r"(?<!\\)\|", line_str)]
            # Strip empty leading and trailing cells from outer pipes
            if cells and cells[0] == "":
                cells.pop(0)
            if cells and cells[-1] == "":
                cells.pop(-1)

            col_count = len(cells)
            if expected_cols is None:
                expected_cols = col_count
            elif col_count != expected_cols:
                errors.append(
                    f"{rel_path}:{line_num}: table starting at line {start_line} has "
                    f"{col_count} columns; expected {expected_cols}"
                )

    def test_angle_bracket_escaping_in_tables(self):
        """Ensure non-HTML angle brackets (e.g. <dup-N>, <DUP>) are backtick-escaped."""
        md_files = find_all_markdown_files()
        escaping_errors = []

        angle_pattern = re.compile(r"<([A-Za-z0-9_\-]+)>")

        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            with open(md_path, "r", encoding="utf-8") as f:
                content = f.read()

            # Remove code blocks and inline code spans before inspecting raw text
            no_code_blocks = re.sub(r"```[\s\S]*?```", "", content)
            no_inline_code = re.sub(r"`[^`]*`", "", no_code_blocks)

            for match in angle_pattern.finditer(no_inline_code):
                tag = match.group(1).lower()
                if tag not in ALLOWED_HTML_TAGS:
                    # Check if it looks like an email or autolink
                    if "@" in tag or "." in tag:
                        continue
                    escaping_errors.append(
                        f"{rel_path}: unescaped angle bracket token '<{match.group(1)}>' "
                        "should be enclosed in backticks."
                    )

        self.assertEqual(
            escaping_errors, [],
            "Unescaped angle brackets found (can break HTML/GFM rendering):\n"
            + "\n".join(f"  - {e}" for e in escaping_errors)
        )

    def test_katex_math_hygiene(self):
        """Verify that markdown files do not contain unescaped underscores in \\text{} or malformed math."""
        md_files = find_all_markdown_files()
        math_errors = []

        # Pattern: \text{..._...} where _ is not escaped with backslash
        text_underscore_pattern = re.compile(r"\\text\{[^}]*(?<!\\)_[^}]*\}")

        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            with open(md_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            in_code_block = False
            for idx, line in enumerate(lines, 1):
                stripped = line.strip()
                if stripped.startswith("```"):
                    in_code_block = not in_code_block
                    continue
                if in_code_block:
                    continue

                for m in text_underscore_pattern.finditer(line):
                    math_errors.append(
                        f"{rel_path}:{idx}: unescaped underscore in '{m.group(0)}' "
                        "(causes KaTeX ParseError: '_' allowed only in math mode; use hyphens or \\_)"
                    )

        self.assertEqual(
            math_errors, [],
            "KaTeX math formatting errors found:\n" + "\n".join(f"  - {e}" for e in math_errors)
        )

    @staticmethod
    def _heading_to_slug(heading_line):
        """Convert a markdown heading line to its GitHub anchor slug.

        Follows GitHub Flavored Markdown (GFM) anchor generation rules:
        - Strip leading ATX hashes and subsequent whitespace.
        - Strip trailing ATX hashes and trailing whitespace.
        - Strip HTML tags (e.g. <code>...</code>, <a>...</a>).
        - Replace markdown link formatting with link text ([text](url) -> text).
        - Remove backticks and control characters.
        - Downcase according to ASCII.
        - Remove GFM punctuation characters (anything other than alphanumeric, whitespace, hyphen).
        - Replace spaces and tabs with hyphens without collapsing consecutive spaces
          (preserving double hyphens caused by removed punctuation such as ' & ').
        """
        text = re.sub(r"^#{1,6}[\t ]+", "", heading_line.strip())
        text = re.sub(r"[\t ]+#+[\t ]*$", "", text)
        text = re.sub(r"<[^>]+>", "", text)
        text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
        text = text.replace("`", "")
        text = re.sub(r"[\x00-\x1f]", "", text)
        text = text.lower()
        text = re.sub(r"[^\w\s-]", "", text)
        text = text.replace("\t", " ")
        slug = text.replace(" ", "-")
        return slug

    @staticmethod
    def _extract_all_anchors(content):
        """Extract all valid anchors (explicit HTML IDs/names and GFM heading slugs) from markdown content."""
        anchors = set()
        # Explicit HTML id/name anchors: <a id="...">, <span id="...">, <div id="...">, etc.
        for m in re.finditer(r'(?:id|name)=["\']([^"\']+)["\']', content, re.IGNORECASE):
            anchors.add(m.group(1))

        # GFM heading slugs: lines starting with 1 to 6 # followed by whitespace
        for line in content.splitlines():
            line_str = line.strip()
            if re.match(r"^#{1,6}\s+", line_str):
                slug = TestMarkdownStaticHygiene._heading_to_slug(line_str)
                if slug:
                    anchors.add(slug)
        return anchors

    def test_heading_to_slug_gfm_parity(self):
        """Verify _heading_to_slug produces exact GFM slugs for complex headings."""
        cases = [
            (
                "### 10.3 Release Packaging & SLSA Level 3 (`.github/workflows/release.yml`)",
                "103-release-packaging--slsa-level-3-githubworkflowsreleaseyml",
            ),
            (
                "## Level 1: Unit & Component Tests",
                "level-1-unit--component-tests",
            ),
            (
                "## Pillar 9: Multi-Modal Public Cohorts",
                "pillar-9-multi-modal-public-cohorts",
            ),
            (
                "### Heading with Trailing Hashes ###",
                "heading-with-trailing-hashes",
            ),
            (
                "### Section with <code>inline_code</code> and [link](https://example.com)",
                "section-with-inline_code-and-link",
            ),
            (
                "### Clean Title\r\n",
                "clean-title",
            ),
        ]
        for heading, expected_slug in cases:
            self.assertEqual(
                self._heading_to_slug(heading),
                expected_slug,
                f"Slug mismatch for '{heading}'",
            )

    def test_extract_all_anchors(self):
        """Verify _extract_all_anchors captures explicit HTML IDs and ATX heading slugs."""
        sample_doc = (
            "# Main Heading\n\n"
            '<a id="custom-anchor"></a>\n'
            "### Sub-Section & Details\n\n"
            '<div name="named-div">content</div>\n'
        )
        anchors = self._extract_all_anchors(sample_doc)
        self.assertIn("main-heading", anchors)
        self.assertIn("custom-anchor", anchors)
        self.assertIn("sub-section--details", anchors)
        self.assertIn("named-div", anchors)

    def test_wiki_links_and_anchors(self):
        """Verify wiki links and anchors are valid and compatible with GitHub Wiki."""
        wiki_dir = REPO_ROOT / "wiki"
        if not wiki_dir.exists():
            return

        wiki_page_files = list(wiki_dir.glob("*.md"))
        wiki_page_stems = {p.stem: p for p in wiki_page_files}

        link_pattern = re.compile(r"(?<!\!)\[([^\]]+)\]\(([^)\s]+)(?:\s+[\"'][^\"']*[\"'])?\)")
        errors = []

        for wiki_file in sorted(wiki_page_files):
            rel_path = wiki_file.relative_to(REPO_ROOT)
            with open(wiki_file, "r", encoding="utf-8") as f:
                content = f.read()

            no_code_blocks = re.sub(r"```[\s\S]*?```", "", content)
            no_inline_code = re.sub(r"`[^`\n]*`", "CODE", no_code_blocks)

            for match in link_pattern.finditer(no_inline_code):
                label = match.group(1)
                raw_url = match.group(2)
                url = urllib.parse.unquote(raw_url)

                # Rule 1: relative link starting with ../ will 404 on GitHub Wiki
                if raw_url.startswith("../") or url.startswith("../"):
                    errors.append(
                        f"{rel_path}: relative link '{raw_url}' for '{label}' will 404 on GitHub Wiki; "
                        "must use absolute GitHub URL or wiki page slug"
                    )
                    continue

                # Rule 2: in-page anchor (#anchor)
                if raw_url.startswith("#"):
                    anchor = raw_url[1:]
                    file_anchors = self._extract_all_anchors(content)
                    if anchor not in file_anchors and urllib.parse.unquote(anchor) not in file_anchors:
                        errors.append(
                            f"{rel_path}: in-page anchor '{raw_url}' for '{label}' not found in {rel_path}"
                        )
                    continue

                # Rule 3: internal wiki page link (extensionless or slug#anchor)
                if not (raw_url.startswith("http://") or raw_url.startswith("https://") or raw_url.startswith("mailto:")):
                    target_page, _, page_anchor = raw_url.partition("#")
                    target_page_clean = urllib.parse.unquote(target_page)

                    if target_page_clean.endswith(".md"):
                        errors.append(
                            f"{rel_path}: internal wiki link '{raw_url}' for '{label}' ends with '.md' "
                            "which triggers raw file download; use slug without .md"
                        )
                        continue

                    # Validate target wiki page existence
                    matched_file = None
                    if target_page_clean in wiki_page_stems:
                        matched_file = wiki_page_stems[target_page_clean]
                    else:
                        norm_target = re.sub(r"[\s_-]+", "-", target_page_clean.lower())
                        for stem, path in wiki_page_stems.items():
                            if re.sub(r"[\s_-]+", "-", stem.lower()) == norm_target:
                                matched_file = path
                                break

                    if not matched_file:
                        errors.append(
                            f"{rel_path}: internal wiki link target page '{target_page}' for '{label}' "
                            "does not exist in wiki/"
                        )
                    elif page_anchor:
                        with open(matched_file, "r", encoding="utf-8") as tf:
                            target_wiki_content = tf.read()
                        target_anchors = self._extract_all_anchors(target_wiki_content)
                        if page_anchor not in target_anchors and urllib.parse.unquote(page_anchor) not in target_anchors:
                            errors.append(
                                f"{rel_path}: anchor '#{page_anchor}' not found in wiki page "
                                f"'{matched_file.name}' (link '{raw_url}')"
                            )
                    continue

                # Rule 4: absolute links pointing to repo blob or tree
                if raw_url.startswith("https://github.com/palakpsheth/VarDictRust/blob/main/") or \
                   raw_url.startswith("https://github.com/palakpsheth/VarDictRust/tree/main/"):
                    prefix = "https://github.com/palakpsheth/VarDictRust/blob/main/" if \
                        raw_url.startswith("https://github.com/palakpsheth/VarDictRust/blob/main/") else \
                        "https://github.com/palakpsheth/VarDictRust/tree/main/"
                    target_ref = raw_url[len(prefix):]
                    if "#" in target_ref:
                        subpath, anchor = target_ref.split("#", 1)
                    else:
                        subpath, anchor = target_ref, None

                    subpath_unquoted = urllib.parse.unquote(subpath)
                    target_file = REPO_ROOT / subpath_unquoted
                    if not target_file.exists():
                        errors.append(
                            f"{rel_path}: target file does not exist: '{subpath}' (link '{raw_url}')"
                        )
                    elif anchor:
                        check_file = target_file if target_file.is_file() else (target_file / "README.md")
                        if not check_file.exists():
                            errors.append(
                                f"{rel_path}: target file for anchor does not exist: '{check_file}' (link '{raw_url}')"
                            )
                        else:
                            with open(check_file, "r", encoding="utf-8") as tf:
                                target_content = tf.read()
                            valid_anchors = self._extract_all_anchors(target_content)
                            if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                                errors.append(
                                    f"{rel_path}: anchor '#{anchor}' not found in '{check_file.relative_to(REPO_ROOT)}' (link '{raw_url}')"
                                )
                    continue

                # Rule 5: root repo link with anchor
                if raw_url.startswith("https://github.com/palakpsheth/VarDictRust") and "#" in raw_url:
                    anchor = raw_url.split("#", 1)[1]
                    readme_path = REPO_ROOT / "README.md"
                    if not readme_path.exists():
                        errors.append(
                            f"{rel_path}: README.md does not exist for anchor check (link '{raw_url}')"
                        )
                    else:
                        with open(readme_path, "r", encoding="utf-8") as rf:
                            readme_content = rf.read()
                        valid_anchors = self._extract_all_anchors(readme_content)
                        if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                            errors.append(
                                f"{rel_path}: anchor '#{anchor}' not found in README.md (link '{raw_url}')"
                            )
                    continue

        self.assertEqual(
            errors, [],
            "Wiki link / anchor hygiene errors found:\n" + "\n".join(f"  - {e}" for e in errors)
        )

    def test_readme_badges_camo_safe(self):
        """Verify README.md badges use Camo-safe endpoints and avoid raw readthedocs.org badge redirects."""
        readme_path = REPO_ROOT / "README.md"
        with open(readme_path, "r", encoding="utf-8") as f:
            content = f.read()
        self.assertNotIn(
            "readthedocs.org/projects/vardictrust/badge/",
            content,
            "README.md must not use raw readthedocs.org/projects/.../badge/ image URLs "
            "(Cloudflare 307/429 causes GitHub Camo HTTP 502); use img.shields.io/readthedocs/vardictrust.",
        )
        self.assertIn(
            "https://img.shields.io/readthedocs/vardictrust",
            content,
            "README.md should use img.shields.io/readthedocs/vardictrust for the Documentation Status badge.",
        )

    def test_repo_markdown_relative_links(self):
        """Verify relative file links and anchors across all tracked repository markdown files."""
        md_files = find_all_markdown_files()
        link_pattern = re.compile(r"(?<!\!)\[([^\]]+)\]\(([^)\s]+)(?:\s+[\"'][^\"']*[\"'])?\)")
        errors = []
        anchor_cache = {}
        tracked_files = set()
        tracked_dirs = set()
        try:
            tracked_out = subprocess.check_output(
                ["git", "-C", str(REPO_ROOT), "ls-files"],
                stderr=subprocess.DEVNULL,
                text=True,
            )
            for line in tracked_out.splitlines():
                if not line.strip():
                    continue
                tf = (REPO_ROOT / line).resolve()
                tracked_files.add(tf)
                parent = tf.parent
                while parent != REPO_ROOT and REPO_ROOT in parent.parents:
                    tracked_dirs.add(parent)
                    parent = parent.parent
            if tracked_files:
                tracked_dirs.add(REPO_ROOT)
        except (OSError, subprocess.CalledProcessError):
            pass

        for md_path in md_files:
            rel_path = md_path.relative_to(REPO_ROOT)
            if str(rel_path).startswith("wiki/"):
                continue

            with open(md_path, "r", encoding="utf-8") as f:
                content = f.read()

            no_code_blocks = re.sub(r"```[\s\S]*?```", "", content)
            no_inline_code = re.sub(r"`[^`\n]*`", "CODE", no_code_blocks)

            for match in link_pattern.finditer(no_inline_code):
                label = match.group(1)
                raw_url = match.group(2)

                if raw_url.startswith(("http://", "https://", "mailto:", "conversation://", "file://")):
                    continue

                if raw_url.startswith("#"):
                    anchor = raw_url[1:]
                    if not anchor or re.match(r"^L\d+(-L\d+)?$", anchor):
                        continue
                    if md_path not in anchor_cache:
                        anchor_cache[md_path] = self._extract_all_anchors(content)
                    valid_anchors = anchor_cache[md_path]
                    if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                        errors.append(
                            f"{rel_path}: broken in-page anchor '{raw_url}' for '{label}'"
                        )
                    continue

                target_part, _, anchor = raw_url.partition("#")
                target_clean = urllib.parse.unquote(target_part)
                if not target_clean:
                    continue

                resolved = (md_path.parent / target_clean).resolve()
                try:
                    resolved.relative_to(REPO_ROOT)
                except ValueError:
                    errors.append(
                        f"{rel_path}: relative link '{raw_url}' escapes repository root ('{resolved}')"
                    )
                    continue
                if not resolved.exists():
                    errors.append(
                        f"{rel_path}: broken relative link '{raw_url}' for '{label}' "
                        f"(target '{resolved}' does not exist)"
                    )
                    continue
                if tracked_files and resolved not in tracked_files and resolved not in tracked_dirs:
                    errors.append(
                        f"{rel_path}: relative link '{raw_url}' for '{label}' points to "
                        f"untracked path '{resolved.relative_to(REPO_ROOT)}'"
                    )
                    continue

                if anchor:
                    if re.match(r"^L\d+(-L\d+)?$", anchor):
                        continue
                    anchor_file = resolved if resolved.is_file() else (resolved / "README.md")
                    if anchor_file.is_file() and anchor_file.suffix.lower() == ".md":
                        if anchor_file not in anchor_cache:
                            with open(anchor_file, "r", encoding="utf-8") as tf:
                                anchor_cache[anchor_file] = self._extract_all_anchors(tf.read())
                        valid_anchors = anchor_cache[anchor_file]
                        if anchor not in valid_anchors and urllib.parse.unquote(anchor) not in valid_anchors:
                            errors.append(
                                f"{rel_path}: broken anchor '#{anchor}' in "
                                f"'{anchor_file.relative_to(REPO_ROOT)}' (link '{raw_url}')"
                            )

        self.assertEqual(
            errors, [],
            "Broken relative markdown links or anchors found:\n" + "\n".join(f"  - {e}" for e in errors),
        )


class TestZensicalBuiltSiteHtml(unittest.TestCase):
    """Verify built HTML output in site/ using pure-stdlib html.parser."""

    def setUp(self):
        self.site_dir = REPO_ROOT / "site"
        index_html = self.site_dir / "index.html"
        require_built = os.environ.get("VARDICTRUST_REQUIRE_BUILT_SITE") == "1"
        if not index_html.is_file():
            if require_built:
                self.fail(
                    "VARDICTRUST_REQUIRE_BUILT_SITE=1 is set, but site/index.html was not found. "
                    "Run 'python3 harness/prepare_docs_site.py && zensical build --clean --strict' first."
                )
            self.skipTest("site/index.html not found; run zensical build to test built HTML.")

    def test_site_html_structure_and_assets(self):
        """Verify built site contains valid navigation, anchor targets, mathjax, and mermaid."""
        from html.parser import HTMLParser

        class SiteParser(HTMLParser):
            def __init__(self):
                super().__init__()
                self.ids = set()
                self.links = []
                self.scripts = []
                self.has_mermaid = False

            def handle_starttag(self, tag, attrs):
                attr_dict = dict(attrs)
                if "id" in attr_dict:
                    self.ids.add(attr_dict["id"])
                if tag == "a" and "href" in attr_dict:
                    self.links.append(attr_dict["href"])
                if tag == "script" and "src" in attr_dict:
                    self.scripts.append(attr_dict["src"])
                classes = attr_dict.get("class", "").split()
                if "mermaid" in classes:
                    self.has_mermaid = True

        index_file = self.site_dir / "index.html"
        parser = SiteParser()
        parser.feed(index_file.read_text(encoding="utf-8"))

        # Verify mathjax script inclusion
        has_mathjax = any("mathjax" in s.lower() for s in parser.scripts)
        self.assertTrue(has_mathjax, "site/index.html must include MathJax script reference")

        # Verify navigation links exist
        self.assertGreater(len(parser.links), 5, "site/index.html must contain navigation links")

        # Check all HTML files in site/
        html_files = sorted(self.site_dir.rglob("*.html"))
        self.assertGreaterEqual(len(html_files), 10, "Expected at least 10 built HTML pages in site/")

        broken_anchors = []
        found_mermaid = False
        for html_file in html_files:
            p = SiteParser()
            try:
                p.feed(html_file.read_text(encoding="utf-8"))
            except Exception as e:
                broken_anchors.append(f"{html_file.relative_to(self.site_dir)}: parse error: {e}")
                continue
            if p.has_mermaid:
                found_mermaid = True
            for link in p.links:
                if link.startswith("#") and len(link) > 1:
                    anchor = urllib.parse.unquote(link[1:])
                    if anchor not in p.ids:
                        broken_anchors.append(
                            f"{html_file.relative_to(self.site_dir)}: anchor '{link}' not found in page IDs"
                        )
        self.assertTrue(found_mermaid, "Built site must contain Mermaid diagram blocks across pages")
        self.assertEqual(
            broken_anchors, [],
            "Broken in-page anchors detected across built HTML files:\n" + "\n".join(f"  - {e}" for e in broken_anchors)
        )


class TestPlaywrightBrowserRendering(unittest.TestCase):
    """Render documents with Playwright headless Chromium to verify real browser DOM layout."""

    sync_playwright_fn = None
    has_playwright = False

    @classmethod
    def setUpClass(cls):
        try:
            playwright_mod = importlib.import_module("playwright.sync_api")
            cls.sync_playwright_fn = staticmethod(playwright_mod.sync_playwright)
            cls.has_playwright = True
        except ImportError:
            cls.has_playwright = False

    def setUp(self):
        if not self.has_playwright:
            if os.environ.get("VARDICTRUST_REQUIRE_PLAYWRIGHT") == "1":
                self.fail("VARDICTRUST_REQUIRE_PLAYWRIGHT=1 is set, but playwright is not installed.")
            self.skipTest("playwright not installed on host; skipping browser rendering test.")

    def _setup_cdn_cache_routes(self, page):
        """Set up disk-caching route handlers for unpkg.com and jsdelivr.net CDN assets."""
        cache_dir = REPO_ROOT / ".cache" / "docs_ui_cdn"
        cache_dir.mkdir(parents=True, exist_ok=True)

        def handle_cdn_route(route):
            url = route.request.url
            url_key = hashlib.sha256(url.encode("utf-8")).hexdigest()
            body_file = cache_dir / f"{url_key}.body"
            meta_file = cache_dir / f"{url_key}.type"
            try:
                if body_file.is_file() and meta_file.is_file():
                    body = body_file.read_bytes()
                    content_type = meta_file.read_text(encoding="utf-8")
                else:
                    req = urllib.request.Request(
                        url,
                        headers={"User-Agent": "Mozilla/5.0 (Playwright Test)"},
                    )
                    with urllib.request.urlopen(req, timeout=15) as resp:
                        body = resp.read()
                        content_type = resp.headers.get_content_type() or "application/javascript"
                    body_file.write_bytes(body)
                    meta_file.write_text(content_type, encoding="utf-8")
                route.fulfill(status=200, content_type=content_type, body=body)
            except Exception:
                route.abort()

        page.route("**/*unpkg.com/**", handle_cdn_route)
        page.route("**/*jsdelivr.net/**", handle_cdn_route)

    @staticmethod
    def _fallback_markdown_to_html(markdown_text):
        """Pure stdlib fallback converting GitHub-Flavored Markdown to HTML when gh is offline/unauthenticated."""
        lines = markdown_text.splitlines()
        out = ['<div class="markdown-body">']
        i = 0
        n = len(lines)

        def render_inline(text):
            tokens = re.split(r"(`[^`]+`)", text)
            res = []
            for tok in tokens:
                if tok.startswith("`") and tok.endswith("`") and len(tok) >= 2:
                    res.append(f"<code>{html.escape(tok[1:-1])}</code>")
                else:
                    link_pattern = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")
                    last_idx = 0
                    for m in link_pattern.finditer(tok):
                        start, end = m.span()
                        res.append(html.escape(tok[last_idx:start]))
                        res.append(f'<a href="{html.escape(m.group(2))}">{html.escape(m.group(1))}</a>')
                        last_idx = end
                    res.append(html.escape(tok[last_idx:]))
            return "".join(res)

        def split_pipe_row(line):
            stripped = line.strip()
            if stripped.startswith("|"):
                stripped = stripped[1:]
            if stripped.endswith("|") and not stripped.endswith(r"\|"):
                stripped = stripped[:-1]
            raw_cells = re.split(r"(?<!\\)\|", stripped)
            return [c.replace(r"\|", "|").strip() for c in raw_cells]

        def is_delimiter_row(cells):
            if not cells:
                return False
            for c in cells:
                if not c or not re.fullmatch(r":?-+:?", c):
                    return False
            return True

        while i < n:
            line = lines[i]
            stripped = line.strip()

            if stripped.startswith("```"):
                code_lines = []
                i += 1
                while i < n and not lines[i].strip().startswith("```"):
                    code_lines.append(lines[i])
                    i += 1
                if i < n:
                    i += 1
                out.append(f"<pre><code>{html.escape(chr(10).join(code_lines))}</code></pre>")
                continue

            if not stripped:
                i += 1
                continue

            heading_match = re.match(r"^(#{1,6})\s+(.*)$", stripped)
            if heading_match:
                level = len(heading_match.group(1))
                heading_content = re.sub(r"\s+#+$", "", heading_match.group(2).strip())
                out.append(f"<h{level}>{render_inline(heading_content)}</h{level}>")
                i += 1
                continue

            if "|" in line and i + 1 < n and "|" in lines[i + 1]:
                header_cells = split_pipe_row(line)
                delim_cells = split_pipe_row(lines[i + 1])
                if is_delimiter_row(delim_cells):
                    table_html = ["<table><thead><tr>"]
                    for h in header_cells:
                        table_html.append(f"<th>{render_inline(h)}</th>")
                    table_html.append("</tr></thead><tbody>")
                    i += 2
                    while i < n and "|" in lines[i] and lines[i].strip():
                        body_cells = split_pipe_row(lines[i])
                        table_html.append("<tr>")
                        for idx in range(len(header_cells)):
                            cell_val = body_cells[idx] if idx < len(body_cells) else ""
                            table_html.append(f"<td>{render_inline(cell_val)}</td>")
                        table_html.append("</tr>")
                        i += 1
                    table_html.append("</tbody></table>")
                    out.append("".join(table_html))
                    continue

            out.append(f"<p>{render_inline(line)}</p>")
            i += 1

        out.append("</div>")
        return "\n".join(out)

    def _render_gfm_html(self, markdown_text):
        """Render markdown text to GitHub-Flavored HTML using gh api or pure-stdlib fallback."""
        try:
            res = subprocess.run(
                ["gh", "api", "/markdown", "-f", f"text={markdown_text}", "-f", "mode=gfm"],
                capture_output=True,
                text=True,
                check=True,
                timeout=10
            )
            if res.stdout and ("<h1" in res.stdout or "<table" in res.stdout or "<p" in res.stdout):
                return res.stdout
        except (subprocess.SubprocessError, FileNotFoundError):
            pass
        return self._fallback_markdown_to_html(markdown_text)

    def test_zensical_built_site_loopback_rendering(self):
        """Serve built site/ via loopback HTTP server and verify browser DOM across 4 key pages in headless Chromium."""
        site_dir = REPO_ROOT / "site"
        index_file = site_dir / "index.html"
        if not index_file.is_file():
            if os.environ.get("VARDICTRUST_REQUIRE_BUILT_SITE") == "1":
                self.fail("VARDICTRUST_REQUIRE_BUILT_SITE=1 is set, but site/index.html was not found.")
            self.skipTest("site/index.html not found; skipping loopback browser test.")

        class QuietHandler(http.server.SimpleHTTPRequestHandler):
            def log_message(self, format, *args):
                pass

        handler = functools.partial(QuietHandler, directory=str(site_dir))
        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        port = server.server_address[1]
        server_thread = threading.Thread(target=server.serve_forever, daemon=True)
        server_thread.start()

        try:
            with self.sync_playwright_fn() as p:
                browser = p.chromium.launch(headless=True)
                page = browser.new_page()

                # Disk-caching route handler for external CDNs (MathJax, fonts)
                self._setup_cdn_cache_routes(page)

                # Page 1: Homepage (/)
                page.goto(f"http://127.0.0.1:{port}/", wait_until="networkidle")
                title = page.title()
                self.assertIn("VarDictRust", title, f"Homepage title must mention VarDictRust: {title}")
                home_nav = page.query_selector_all("nav a, header a")
                self.assertGreater(len(home_nav), 0, "Homepage must render navigation elements")
                home_tables = page.query_selector_all("table")
                self.assertGreaterEqual(len(home_tables), 1, "Homepage must render benchmark table(s)")

                # Page 2: /ROADMAP/
                page.goto(f"http://127.0.0.1:{port}/ROADMAP/", wait_until="networkidle")
                roadmap_nav = page.query_selector_all("nav a, header a")
                self.assertGreater(len(roadmap_nav), 0, "ROADMAP page must render navigation elements")
                roadmap_tables = page.query_selector_all("table")
                self.assertGreaterEqual(len(roadmap_tables), 4, "ROADMAP page must render at least 4 tables")
                # Table 1: v1.1.0 candidates (4 columns, >= 10 rows)
                self.assertEqual(len(roadmap_tables[0].query_selector_all("th")), 4)
                self.assertGreaterEqual(len(roadmap_tables[0].query_selector_all("tbody tr")), 10)
                # Table 2: Landed as Default (4 columns, >= 2 rows)
                self.assertEqual(len(roadmap_tables[1].query_selector_all("th")), 4)
                self.assertGreaterEqual(len(roadmap_tables[1].query_selector_all("tbody tr")), 2)
                # Table 3: Needs a specification first (4 columns, >= 7 rows)
                self.assertEqual(len(roadmap_tables[2].query_selector_all("th")), 4)
                self.assertGreaterEqual(len(roadmap_tables[2].query_selector_all("tbody tr")), 7)
                # Table 4: Open upstream PRs reviewed (3 columns, >= 4 rows)
                self.assertEqual(len(roadmap_tables[3].query_selector_all("th")), 3)
                self.assertGreaterEqual(len(roadmap_tables[3].query_selector_all("tbody tr")), 4)

                # Page 3: /paper/manuscript/
                page.goto(f"http://127.0.0.1:{port}/paper/manuscript/", wait_until="networkidle")
                page.wait_for_selector("mjx-container", timeout=10000)
                paper_math = page.query_selector_all("mjx-container")
                self.assertGreater(len(paper_math), 5, "Manuscript must render >5 MathJax mjx-container equations")
                paper_math_errors = page.query_selector_all("mjx-merror")
                self.assertEqual(len(paper_math_errors), 0, "Manuscript must have 0 MathJax mjx-merror syntax errors")
                paper_tables = page.query_selector_all("table")
                self.assertGreaterEqual(len(paper_tables), 1, "Manuscript must render >= 1 tables")
                paper_svgs = page.query_selector_all('img[src$=".svg"]')
                self.assertGreaterEqual(len(paper_svgs), 1, "Manuscript must render SVG images")
                for svg_el in paper_svgs:
                    nw = svg_el.evaluate("el => el.naturalWidth")
                    self.assertGreater(nw, 0, f"SVG image must have naturalWidth > 0; got {nw} for {svg_el.get_attribute('src')}")

                # Page 4: /docs/validation/omnibus-empirical-validation/
                page.goto(f"http://127.0.0.1:{port}/docs/validation/omnibus-empirical-validation/", wait_until="networkidle")
                page.wait_for_selector("mjx-container", timeout=10000)
                val_math = page.query_selector_all("mjx-container")
                self.assertGreater(len(val_math), 5, "Omnibus validation must render >5 MathJax mjx-container equations")
                val_math_errors = page.query_selector_all("mjx-merror")
                self.assertEqual(len(val_math_errors), 0, "Omnibus validation must have 0 MathJax mjx-merror syntax errors")
                val_tables = page.query_selector_all("table")
                self.assertGreaterEqual(len(val_tables), 5, "Omnibus validation must render >= 5 tables")
                val_mermaids = page.query_selector_all(".mermaid, pre.mermaid")
                self.assertGreaterEqual(len(val_mermaids), 1, "Omnibus validation must contain Mermaid diagram element(s)")

                browser.close()
        finally:
            server.shutdown()
            server.server_close()

    def test_roadmap_browser_rendering(self):
        """Verify ROADMAP.md renders clean tables and links in headless Chromium."""
        roadmap_path = REPO_ROOT / "ROADMAP.md"
        with open(roadmap_path, "r", encoding="utf-8") as f:
            content = f.read()

        rendered_html = self._render_gfm_html(content)
        page_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>VarDictRust Roadmap</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; padding: 24px; }}
table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
th, td {{ border: 1px solid #d0d7de; padding: 6px 13px; text-align: left; }}
th {{ background-color: #f6f8fa; }}
code {{ font-family: monospace; background-color: rgba(175, 184, 193, 0.2); padding: 0.2em 0.4em; border-radius: 6px; }}
</style>
</head>
<body>
{rendered_html}
</body>
</html>
"""

        with self.sync_playwright_fn() as p:
            browser = p.chromium.launch(headless=True)
            page = browser.new_page()
            self._setup_cdn_cache_routes(page)
            page.set_content(page_content)

            # 1. Verify H1 header
            h1_text = page.inner_text("h1")
            self.assertIn("VarDictRust Roadmap", h1_text)

            # 2. Verify all 4 tables exist
            tables = page.query_selector_all("table")
            self.assertGreaterEqual(len(tables), 4, "ROADMAP.md should contain at least 4 tables.")

            # Table 1: v1.1.0 candidates (4 columns, >= 10 rows)
            table1 = tables[0]
            headers1 = [th.inner_text().strip() for th in table1.query_selector_all("th")]
            self.assertEqual(len(headers1), 4)
            rows1 = table1.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows1), 10)

            # Table 2: Landed as Default (4 columns, >= 2 rows)
            table2 = tables[1]
            headers2 = [th.inner_text().strip() for th in table2.query_selector_all("th")]
            self.assertEqual(len(headers2), 4)
            rows2 = table2.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows2), 2)

            # Table 3: Needs a specification first (4 columns, >= 7 rows)
            table3 = tables[2]
            headers3 = [th.inner_text().strip() for th in table3.query_selector_all("th")]
            self.assertEqual(len(headers3), 4)
            rows3 = table3.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows3), 7)

            # Table 4: Open upstream PRs reviewed (3 columns, >= 4 rows)
            table4 = tables[3]
            headers4 = [th.inner_text().strip() for th in table4.query_selector_all("th")]
            self.assertEqual(len(headers4), 3)
            rows4 = table4.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows4), 4)

            # 3. Verify links point to AstraZeneca-NGS issues/PRs or docs
            links = page.query_selector_all("a")
            self.assertGreater(len(links), 15, "Should have links for issues, PRs, and docs.")
            for link in links:
                href = link.get_attribute("href")
                self.assertIsNotNone(href)
                self.assertTrue(
                    href.startswith("https://github.com/AstraZeneca-NGS/VarDictJava/")
                    or href.startswith("docs/")
                    or href.startswith("#"),
                    f"Unexpected link URL: {href}"
                )

            # 4. Capture screenshot for verification
            with tempfile.TemporaryDirectory() as tmpdir:
                screenshot_path = os.path.join(tmpdir, "roadmap_browser_rendered.png")
                page.screenshot(path=screenshot_path, full_page=True)
                self.assertTrue(os.path.exists(screenshot_path))

            browser.close()

    def test_readme_browser_rendering(self):
        """Verify README.md renders cleanly with benchmark tables in headless Chromium."""
        readme_path = REPO_ROOT / "README.md"
        with open(readme_path, "r", encoding="utf-8") as f:
            content = f.read()

        rendered_html = self._render_gfm_html(content[:8000])  # First 8KB covers hook, badges, benchmarks
        page_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>VarDictRust README</title>
</head>
<body>
{rendered_html}
</body>
</html>
"""

        with self.sync_playwright_fn() as p:
            browser = p.chromium.launch(headless=True)
            page = browser.new_page()
            self._setup_cdn_cache_routes(page)
            page.set_content(page_content)

            h1_text = page.inner_text("h1")
            self.assertIn("VarDictRust", h1_text)

            tables = page.query_selector_all("table")
            self.assertGreaterEqual(len(tables), 1, "README.md should render empirical benchmark table.")

            browser.close()

    def test_omnibus_validation_rendering_and_katex(self):
        """Verify omnibus-empirical-validation.md renders cleanly with zero KaTeX errors."""
        omnibus_path = REPO_ROOT / "docs" / "validation" / "omnibus-empirical-validation.md"
        with open(omnibus_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Render first 25KB covering executive summary, master roster, NIST ROC math, and derivations
        sample_content = content[:25000]
        rendered_html = self._render_gfm_html(sample_content)

        page_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>VarDictRust Omnibus Validation</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.css">
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; padding: 24px; }}
table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
th, td {{ border: 1px solid #d0d7de; padding: 6px 13px; text-align: left; }}
th {{ background-color: #f6f8fa; }}
.katex-error {{ color: #cf222e; background-color: #ffebe9; padding: 2px 4px; border: 1px solid #cf222e; }}
</style>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.8/dist/contrib/auto-render.min.js"
        onload="renderMathInElement(document.body, {{delimiters: [{{left: '$$', right: '$$', display: true}}, {{left: '$', right: '$', display: false}}], throwOnError: false}});"></script>
</head>
<body>
{rendered_html}
</body>
</html>
"""

        with self.sync_playwright_fn() as p:
            browser = p.chromium.launch(headless=True)
            page = browser.new_page()
            self._setup_cdn_cache_routes(page)
            page.set_content(page_content)

            # Wait for KaTeX auto-render to process math
            try:
                page.wait_for_selector(".katex", timeout=5000)
            except Exception:
                pass

            # 1. Assert zero KaTeX parse error boxes
            katex_errors = page.query_selector_all(".katex-error")
            error_texts = [e.inner_text().strip() for e in katex_errors]
            self.assertEqual(
                error_texts, [],
                f"KaTeX parse error boxes detected in rendered DOM: {error_texts}"
            )

            # 2. Assert no KaTeX underscore error in page text
            content_html = page.content()
            self.assertNotIn("allowed only in math mode", content_html)

            # 3. Assert math elements were rendered
            katex_nodes = page.query_selector_all(".katex")
            self.assertGreater(len(katex_nodes), 5, "Should render mathematical expressions via KaTeX.")

            # 4. Save verification screenshot
            with tempfile.TemporaryDirectory() as tmpdir:
                screenshot_path = os.path.join(tmpdir, "omnibus_consolidated_browser_rendered.png")
                page.screenshot(path=screenshot_path, full_page=True)
                self.assertTrue(os.path.exists(screenshot_path))

            browser.close()


if __name__ == "__main__":
    unittest.main()
