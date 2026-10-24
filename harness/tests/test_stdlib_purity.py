#!/usr/bin/env python3
"""
harness/tests/test_stdlib_purity.py

Automated AST verification that all Python scripts in the harness and test suite
depend strictly on the Python Standard Library (zero PyPI packages).

Guarantees the VarDictRust zero-dependency harness contract: any standard Python
3.8+ interpreter can execute all tests and comparators immediately without
activating a virtual environment or installing pip wheels.
"""

import ast
import glob
import os
import pathlib
import sys
import unittest

# Standard library module names fallback for Python 3.8 / 3.9
FALLBACK_STDLIB = frozenset([
    "abc", "argparse", "ast", "asyncio", "base64", "bisect", "builtins",
    "bz2", "calendar", "cmath", "collections", "concurrent", "configparser",
    "contextlib", "copy", "csv", "ctypes", "dataclasses", "datetime",
    "decimal", "difflib", "dis", "distutils", "doctest", "email", "enum",
    "errno", "faulthandler", "filecmp", "fileinput", "fnmatch", "fractions",
    "functools", "gc", "getopt", "getpass", "gettext", "glob", "gzip",
    "hashlib", "heapq", "hmac", "html", "http", "idlelib", "imaplib",
    "imghdr", "importlib", "inspect", "io", "ipaddress", "itertools", "json",
    "keyword", "linecache", "locale", "logging", "lzma", "mailbox", "mailcap",
    "marshal", "math", "mimetypes", "mmap", "modulefinder", "multiprocessing",
    "netrc", "nntplib", "numbers", "operator", "optparse", "os", "pathlib",
    "pdb", "pickle", "pickletools", "pkgutil", "platform", "plistlib", "poplib",
    "posix", "posixpath", "pprint", "profile", "pstats", "pty", "pwd", "py_compile",
    "pyclbr", "pydoc", "queue", "quopri", "random", "re", "readline", "reprlib",
    "resource", "rlcompleter", "runpy", "sched", "secrets", "select", "selectors",
    "shelve", "shlex", "shutil", "signal", "site", "smtpd", "smtplib", "sndhdr",
    "socket", "socketserver", "spwd", "sqlite3", "ssl", "stat", "statistics",
    "string", "stringprep", "struct", "subprocess", "sunau", "symbol", "symtable",
    "sys", "sysconfig", "syslog", "tabnanny", "tarfile", "telnetlib", "tempfile",
    "termios", "test", "textwrap", "threading", "time", "timeit", "tkinter",
    "token", "tokenize", "trace", "traceback", "tracemalloc", "tty", "turtle",
    "turtledemo", "types", "typing", "unicodedata", "unittest", "urllib", "uu",
    "uuid", "venv", "warnings", "wave", "weakref", "webbrowser", "wsgiref",
    "xdrlib", "xml", "xmlrpc", "zipapp", "zipfile", "zipimport", "zlib"
])

ALLOWED_LOCAL_MODULES = frozenset([
    "select_slice_regions",
    "compare_vcf",
    "compare_var",
    "perf_budget",
    "record_evidence",
    "generate_validation_table",
    "convert_fastas",
])


class TestStdlibPurity(unittest.TestCase):
    """Verify that harness scripts have zero external PyPI dependencies."""

    def setUp(self):
        # Resolve repository root
        self.script_dir = pathlib.Path(__file__).resolve().parent
        self.harness_dir = self.script_dir.parent
        self.repo_root = self.harness_dir.parent

        if hasattr(sys, "stdlib_module_names"):
            self.stdlib_names = sys.stdlib_module_names
        else:
            self.stdlib_names = FALLBACK_STDLIB

    def test_all_harness_python_files_are_pure_stdlib(self):
        """Parse all Python files in harness/ and .agents/ and ensure every import is in stdlib."""
        py_files = sorted(
            glob.glob(str(self.harness_dir / "**" / "*.py"), recursive=True)
            + glob.glob(str(self.repo_root / ".agents" / "**" / "*.py"), recursive=True)
        )

        self.assertGreaterEqual(
            len(py_files), 10,
            f"Expected at least 10 Python files in {self.harness_dir}, found {len(py_files)}"
        )

        violations = []

        for file_path_str in py_files:
            file_path = pathlib.Path(file_path_str)
            rel_path = file_path.relative_to(self.repo_root)

            with open(file_path, "r", encoding="utf-8") as fh:
                try:
                    tree = ast.parse(fh.read(), filename=str(file_path))
                except SyntaxError as e:
                    self.fail(f"Syntax error in {rel_path}: {e}")

            for node in ast.walk(tree):
                imported_mods = []
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        imported_mods.append(alias.name)
                elif isinstance(node, ast.ImportFrom):
                    if node.module and node.level == 0:
                        imported_mods.append(node.module)

                for mod in imported_mods:
                    top_level = mod.split(".")[0]
                    # Check if top-level module is in stdlib or permitted local module
                    if top_level not in self.stdlib_names and top_level not in ALLOWED_LOCAL_MODULES:
                        violations.append(
                            f"{rel_path}: line {getattr(node, 'lineno', '?')} imports '{mod}' "
                            f"(top-level '{top_level}' is not in Python stdlib or local harness)"
                        )

        self.assertEqual(
            violations, [],
            "Detected unauthorized third-party PyPI dependencies in harness scripts:\n"
            + "\n".join(f"  - {v}" for v in violations)
            + "\n\nVarDictRust requires all harness scripts to use 100% Python standard library."
        )


if __name__ == "__main__":
    unittest.main()
