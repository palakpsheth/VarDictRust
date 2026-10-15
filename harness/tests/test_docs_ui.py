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

import html
import importlib
import json
import os
import pathlib
import re
import subprocess
import sys
import unittest
import urllib.parse

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
        # Exclude git internals, cargo target, conda environments, scratch, and transient dev directories
        rel = str(path.relative_to(REPO_ROOT))
        if any(rel.startswith(ignore) for ignore in [
            ".git/", ".superpowers/", ".agents/", "target/", "scratch/", "VarDictJava/", "node_modules/", ".worktrees/"
        ]):
            continue
        md_files.append(path)
    return sorted(md_files)


class TestMarkdownStaticHygiene(unittest.TestCase):
    """Test static formatting, table alignment, and escaping across all markdown documents."""

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
            no_inline_code = re.sub(r"`[^`]*`", "", no_code_blocks)

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
            self.skipTest("playwright not installed on host; skipping browser rendering test.")

    def _render_gfm_html(self, markdown_text):
        """Render markdown text to GitHub-Flavored HTML using gh api or fallback."""
        try:
            res = subprocess.run(
                ["gh", "api", "/markdown", "-f", f"text={markdown_text}", "-f", "mode=gfm"],
                capture_output=True,
                text=True,
                check=True,
                timeout=10
            )
            return res.stdout
        except (subprocess.SubprocessError, FileNotFoundError):
            # Fallback wrapper for offline/isolated execution
            escaped = html.escape(markdown_text)
            return f"<pre>{escaped}</pre>"

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
            page.set_content(page_content)

            # 1. Verify H1 header
            h1_text = page.inner_text("h1")
            self.assertIn("VarDictRust Roadmap", h1_text)

            # 2. Verify all 3 tables exist
            tables = page.query_selector_all("table")
            self.assertGreaterEqual(len(tables), 3, "ROADMAP.md should contain at least 3 tables.")

            # Table 1: v1.1.0 candidates (4 columns, >= 10 rows)
            table1 = tables[0]
            headers1 = [th.inner_text().strip() for th in table1.query_selector_all("th")]
            self.assertEqual(len(headers1), 4)
            rows1 = table1.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows1), 10)

            # Table 2: Needs a specification first (4 columns, >= 7 rows)
            table2 = tables[1]
            headers2 = [th.inner_text().strip() for th in table2.query_selector_all("th")]
            self.assertEqual(len(headers2), 4)
            rows2 = table2.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows2), 7)

            # Table 3: Open upstream PRs reviewed (3 columns, >= 4 rows)
            table3 = tables[2]
            headers3 = [th.inner_text().strip() for th in table3.query_selector_all("th")]
            self.assertEqual(len(headers3), 3)
            rows3 = table3.query_selector_all("tbody tr")
            self.assertGreaterEqual(len(rows3), 4)

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
            screenshot_path = "/tmp/roadmap_browser_rendered.png"
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
            screenshot_path = "/tmp/omnibus_consolidated_browser_rendered.png"
            page.screenshot(path=screenshot_path, full_page=True)
            self.assertTrue(os.path.exists(screenshot_path))

            browser.close()


if __name__ == "__main__":
    unittest.main()
