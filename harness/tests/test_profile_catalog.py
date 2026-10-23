#!/usr/bin/env python3
"""harness/tests/test_profile_catalog.py: harness/profile_catalog.py (CPU / dhat / heap catalog), stdlib only.

dhat fixtures are synthetic but use the exact dhat-0.3.3 frame-table shape captured from the micro fixture
(`0xADDR: symbol (file:line:col)`, `fs` innermost-first, ftbl[0] = "[root]"): a real capture is 1.5-3.5 MB, far
over the 20 KB fixture budget. The heap fixture (data/heap-micro.pprof.gz, 7.9 KB) is a REAL jemalloc_pprof
capture from fixtures/micro (VARDICT_PROFILE_HEAP, profile-heap build).
"""
import gzip
import json
import os
import pathlib
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
CAT = ROOT / "harness" / "profile_catalog.py"
HEAP_FIXTURE = ROOT / "harness" / "tests" / "data" / "heap-micro.pprof.gz"


def run_cat(d, extra_args=(), env=None):
    out = pathlib.Path(d) / "catalog.md"
    e = dict(os.environ)
    if env:
        e.update(env)
    subprocess.run([sys.executable, str(CAT), str(d), "-o", str(out), "--min-pct", "2", *extra_args],
                   check=True, env=e, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return out.read_text()


def dhat_doc():
    # Real-format frame strings; fs is innermost-first.
    return {"pps": [
        {"tb": 900, "tbk": 90, "mb": 10, "fs": [1, 2, 3]},   # leaf, std frame, vardict caller
        {"tb": 100, "tbk": 10, "mb": 5, "fs": [1, 4]}],
        "ftbl": ["[root]",
                 "0x5f1f1ddfddc8: alloc::raw_vec::RawVec<T,A>::grow_one (src/raw_vec/mod.rs:188:29)",
                 "0x5f1f1ddfd59a: alloc::vec::Vec<T,A>::push (src/vec/mod.rs:1:1)",
                 "0x5f1f1df1e98a: vardictrust_core::variations::get_variation (src/variations.rs:287:28)",
                 "0x5f1f1df1e98b: vardictrust_core::output::table::fmt_1dp (src/output/table.rs:9:9)"]}


def varint(n):
    out = bytearray()
    while True:
        b = n & 0x7F
        n >>= 7
        if n:
            out.append(b | 0x80)
        else:
            out.append(b)
            return bytes(out)


def fld(num, payload):
    """length-delimited field"""
    return varint(num << 3 | 2) + varint(len(payload)) + payload


def vfld(num, v):
    return varint(num << 3) + varint(v)


def pprof_bytes(samples, locations, functions, strings, sample_types=((1, 2),)):
    """Tiny profile.proto encoder (subset the catalog decodes).
    locations: {id: [function_id inner-first]}; functions: {id: string-index}; samples: [([loc ids leaf-first], [values])]"""
    b = b""
    for t, u in sample_types:
        b += fld(1, vfld(1, t) + vfld(2, u))
    for locs, vals in samples:
        b += fld(2, fld(1, b"".join(varint(x) for x in locs)) + fld(2, b"".join(varint(v) for v in vals)))
    for lid, fns in locations.items():
        b += fld(4, vfld(1, lid) + vfld(3, 0x1000 + lid) + b"".join(fld(4, vfld(1, f)) for f in fns))
    for fid, name in functions.items():
        b += fld(5, vfld(1, fid) + vfld(2, name))
    for s in strings:
        b += fld(6, s.encode())
    return b


class CatalogTest(unittest.TestCase):
    def tmp(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)   # no persistent artifacts
        return pathlib.Path(d)

    # ---- dhat -------------------------------------------------------------------------------
    def test_dhat_attributes_to_vardict_caller_not_alloc_leaf(self):
        d = self.tmp()
        (d / "dhat-heap.json").write_text(json.dumps(dhat_doc()))
        text = run_cat(d)
        self.assertLess(text.index("get_variation"), text.index("fmt_1dp"))
        self.assertNotIn("RawVec<T,A>::grow_one` |", text)   # allocator leaf must not be a ranked row
        self.assertNotIn("RawVec::grow_one` |", text)
        self.assertIn("90.0%", text)

    def test_dhat_real_format_frames_are_cleaned_and_columns_present(self):
        d = self.tmp()
        (d / "dhat-heap.json").write_text(json.dumps(dhat_doc()))
        text = run_cat(d)
        self.assertNotIn("0x5f1f", text.split("## ", 1)[1].split("Top stacks")[0])  # addresses stripped in rows
        for col in ("total_bytes", "total_blocks", "max_bytes"):
            self.assertIn(col, text)

    def test_dhat_rows_below_min_pct_fold_into_other(self):
        d = self.tmp()
        doc = dhat_doc()
        doc["pps"].append({"tb": 1, "tbk": 1, "mb": 1, "fs": [1, 5]})
        doc["ftbl"].append("0x1: vardictrust_core::tiny::rare (src/t.rs:1:1)")
        (d / "dhat-heap.json").write_text(json.dumps(doc))
        text = run_cat(d)
        self.assertNotIn("tiny::rare` |", text)
        self.assertIn("other", text)

    def test_dhat_collapsed_full_stacks_root_to_leaf(self):
        d = self.tmp()
        (d / "dhat-heap.json").write_text(json.dumps(dhat_doc()))
        run_cat(d)
        bytes_lines = (d / "dhat-bytes-collapsed.txt").read_text().splitlines()
        blocks_lines = (d / "dhat-blocks-collapsed.txt").read_text().splitlines()
        # fs [1,2,3] is innermost-first => root->leaf is get_variation;push;grow_one
        self.assertIn("vardictrust_core::variations::get_variation;alloc::vec::Vec<T,A>::push;"
                      "alloc::raw_vec::RawVec<T,A>::grow_one 900", bytes_lines)
        self.assertIn("vardictrust_core::variations::get_variation;alloc::vec::Vec<T,A>::push;"
                      "alloc::raw_vec::RawVec<T,A>::grow_one 90", blocks_lines)
        self.assertEqual(len(bytes_lines), 2)

    def test_dhat_keeps_at_most_top_three_stacks_per_row(self):
        d = self.tmp()
        doc = {"pps": [{"tb": 100 + i, "tbk": 1, "mb": 1, "fs": [1 + i, 7]} for i in range(6)],
               "ftbl": ["[root]"] + ["leaf%d" % i for i in range(6)] + ["vardictrust_core::row"]}
        (d / "dhat-heap.json").write_text(json.dumps(doc))
        text = run_cat(d, ("--min-pct", "0"))
        top = text.split("Top stacks", 1)[1]
        self.assertEqual(top.count(" B: "), 3)
        for i in (5, 4, 3):    # the three biggest survive
            self.assertIn("leaf%d" % i, top)
        self.assertNotIn("leaf0", top)

    # ---- cpu --------------------------------------------------------------------------------
    def test_cpu_collapsed_self_time_and_drop_fraction(self):
        d = self.tmp()
        (d / "cpu-collapsed.txt").write_text("main;run;get_variation 105\nmain;run;fmt_1dp 30\nmain;run;idle 15\n")
        (d / "cpu-meta.json").write_text(json.dumps({"ticks_total": 200, "samples_kept": 150}))
        text = run_cat(d)
        self.assertIn("70.0%", text)   # 105 of 150 kept samples
        # asymmetric on purpose: kept (75.0%) and dropped (25.0%) cannot be confused
        self.assertIn("(75.0% kept, 25.0% dropped;", text)
        self.assertIn("try_write", text)
        self.assertLess(text.index("get_variation"), text.index("fmt_1dp"))

    def test_cpu_attributes_to_innermost_vardict_frame(self):
        d = self.tmp()
        (d / "cpu-collapsed.txt").write_text(
            "t;vardictrust::main;vardictrust_core::run;vardictrust_core::parse;core::slice::sort 60\n"
            "t;vardictrust::main;vardictrust_core::run;libc::memcpy 40\n")
        text = run_cat(d)
        # 60 ticks sit under vardictrust_core::parse, 40 under vardictrust_core::run
        self.assertLess(text.index("vardictrust_core::parse"), text.index("vardictrust_core::run` |"))
        self.assertIn("60.0%", text)

    def test_cpu_hash_suffix_stripped_and_no_meta_ok(self):
        d = self.tmp()
        (d / "cpu-collapsed.txt").write_text("t;a::b::h0123456789abcdef 10\nt;a::b::h1111111111111111 10\n")
        text = run_cat(d)
        self.assertIn("100.0%", text)
        self.assertNotIn("h0123456789abcdef", text)

    # ---- heap -------------------------------------------------------------------------------
    def test_heap_real_pprof_fixture_decodes_and_attributes(self):
        d = self.tmp()
        (d / "heap.pb.gz").write_bytes(HEAP_FIXTURE.read_bytes())
        text = run_cat(d)
        self.assertIn("vardictrust_core::", text)
        lines = (d / "heap-collapsed.txt").read_text().splitlines()
        self.assertTrue(lines)
        # weights are the in-use bytes of the (single) sample type; every line is `stack weight`
        total = sum(int(ln.rsplit(" ", 1)[1]) for ln in lines)
        self.assertGreater(total, 0)
        # root->leaf: jemalloc's innermost frame (prof_backtrace_impl) is the LAST frame, never the first
        stacks = [ln.rsplit(" ", 1)[0].split(";") for ln in lines]
        self.assertTrue(all(st[0] != "prof_backtrace_impl" and st[-1] == "prof_backtrace_impl" for st in stacks), stacks)
        self.assertTrue(any(f.startswith("vardictrust_core::") for st in stacks for f in st))

    def test_heap_synthetic_inline_order_and_sample_type_choice(self):
        d = self.tmp()
        strings = ["", "alloc_objects", "count", "inuse_space", "bytes",
                   "malloc_leaf", "vardictrust_core::a::outer", "std::mid", "vardictrust_core::a::inlined"]
        # loc1: inlined pair [inner=fn1(malloc_leaf)]; loc2: [inner=fn3 inlined vardict, outer=fn2 vardict::outer]; loc3 std
        raw = pprof_bytes(
            samples=[([1, 3, 2], [7, 4096])],
            locations={1: [1], 2: [3, 2], 3: [4]},
            functions={1: 5, 2: 6, 3: 8, 4: 7},
            strings=strings,
            sample_types=((1, 2), (3, 4)))
        (d / "heap.pb.gz").write_bytes(gzip.compress(raw))
        text = run_cat(d)
        line = (d / "heap-collapsed.txt").read_text().strip()
        # leaf-first locs [1,3,2]; each loc inner-first => leaf-first frames:
        # malloc_leaf, std::mid, vardict::inlined, vardict::outer  => root->leaf reversed
        self.assertEqual(line, "vardictrust_core::a::outer;vardictrust_core::a::inlined;std::mid;malloc_leaf 4096")
        # innermost vardictrust_core frame is the inlined one
        self.assertIn("vardictrust_core::a::inlined", text.split("Heap", 1)[1])
        self.assertIn("4096", text)  # inuse_space (idx 1), not alloc_objects count 7

    def test_heap_without_function_names_reports_addresses_only(self):
        d = self.tmp()
        raw = pprof_bytes(samples=[([1], [64])], locations={1: []}, functions={}, strings=["", "inuse_space", "bytes"])
        (d / "heap.pb.gz").write_bytes(gzip.compress(raw))
        text = run_cat(d)
        self.assertIn("no function names", text)

    # ---- corrupt / truncated input: warn, keep the other files' catalog ---------------------------------
    def run_raw(self, d, *extra):
        out = pathlib.Path(d) / "catalog.md"
        r = subprocess.run([sys.executable, str(CAT), str(d), "-o", str(out), *extra],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        return r, out.read_text() if out.exists() else ""

    def test_truncated_gzip_heap_is_reported_other_files_kept(self):
        d = self.tmp()
        (d / "cpu-collapsed.txt").write_text("t;vardictrust_core::ok 5\n")
        (d / "heap.pb.gz").write_bytes(HEAP_FIXTURE.read_bytes()[:40])   # gzip stream cut short
        r, text = self.run_raw(d)
        self.assertEqual(r.returncode, 0)
        self.assertIn("vardictrust_core::ok", text)
        self.assertIn("heap.pb.gz` unreadable", text)
        self.assertIn("WARNING", r.stderr)
        self.assertIn("heap.pb.gz unreadable", r.stderr)

    def test_bad_collapsed_line_is_reported_other_files_kept(self):
        d = self.tmp()
        (d / "cpu-collapsed.txt").write_text("t;a 5\nthis line has no count\n")
        (d / "dhat-heap.json").write_text(json.dumps(dhat_doc()))
        r, text = self.run_raw(d)
        self.assertEqual(r.returncode, 0)
        self.assertIn("cpu-collapsed.txt` unreadable", text)
        self.assertIn("get_variation", text)   # dhat still catalogued

    def test_mid_message_truncated_pb_is_reported(self):
        d = self.tmp()
        raw = gzip.decompress(HEAP_FIXTURE.read_bytes())
        (d / "heap.pb.gz").write_bytes(gzip.compress(raw[:len(raw) // 2]))
        (d / "cpu-collapsed.txt").write_text("t;vardictrust_core::ok 5\n")
        r, text = self.run_raw(d)
        self.assertEqual(r.returncode, 0)
        self.assertIn("heap.pb.gz` unreadable", text)
        self.assertIn("truncated protobuf", text)
        self.assertIn("vardictrust_core::ok", text)

    def test_corrupt_dhat_json_reported_and_no_svg_attempt(self):
        d = self.tmp()
        (d / "dhat-heap.json").write_text("{not json")
        r, text = self.run_raw(d)
        self.assertEqual(r.returncode, 0)
        self.assertIn("dhat-heap.json` unreadable", text)
        self.assertFalse((d / "dhat-bytes-collapsed.txt").exists())

    def test_strict_flag_exits_nonzero_on_unreadable(self):
        d = self.tmp()
        (d / "cpu-collapsed.txt").write_text("garbage\n")
        r, text = self.run_raw(d, "--strict")
        self.assertEqual(r.returncode, 1)
        self.assertIn("unreadable", text)   # catalog still written

    # ---- partial dirs / svg links ------------------------------------------------------------
    def test_partial_dir_empty_and_multiple_dirs(self):
        d1, d2 = self.tmp(), self.tmp()
        (d2 / "cpu-collapsed.txt").write_text("t;x 5\n")
        out = d1 / "catalog.md"
        subprocess.run([sys.executable, str(CAT), str(d1), str(d2), "-o", str(out)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        text = out.read_text()
        self.assertIn(d1.name, text)
        self.assertIn(d2.name, text)
        self.assertIn("no profile data", text)

    def test_svgs_linked_relative_to_catalog(self):
        root = self.tmp()
        cell = root / "cellA"
        cell.mkdir()
        (cell / "cpu-collapsed.txt").write_text("t;x 5\n")
        (cell / "cpu-flamegraph.svg").write_text("<svg/>")
        (cell / "heap-flamegraph.svg").write_text("<svg/>")
        out = root / "catalog.md"
        subprocess.run([sys.executable, str(CAT), str(cell), "-o", str(out)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        text = out.read_text()
        self.assertIn("](cellA/cpu-flamegraph.svg)", text)
        self.assertIn("](cellA/heap-flamegraph.svg)", text)
        self.assertNotIn(str(root), text.replace("cellA", ""))  # relative links only, no absolute paths

    # ---- inferno-flamegraph hook ------------------------------------------------------------
    def stub_path(self, log, with_stub):
        sd = self.tmp()
        if with_stub:
            f = sd / "inferno-flamegraph"
            f.write_text(f'#!/usr/bin/env bash\necho "$@" >> {log}\ncat > /dev/null\necho "<svg stub/>"\n')
            f.chmod(f.stat().st_mode | stat.S_IEXEC)
        # keep python + coreutils reachable but make sure no real inferno-flamegraph shadows the stub-less case
        keep = [p for p in os.environ["PATH"].split(":") if not (pathlib.Path(p) / "inferno-flamegraph").exists()]
        return str(sd) + ":" + ":".join(keep)

    def test_svg_step_skipped_silently_when_inferno_absent(self):
        d = self.tmp()
        (d / "dhat-heap.json").write_text(json.dumps(dhat_doc()))
        out = d / "catalog.md"
        r = subprocess.run([sys.executable, str(CAT), str(d), "-o", str(out)], check=True,
                           env={**os.environ, "PATH": self.stub_path(None, False)},
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.assertFalse((d / "dhat-bytes.svg").exists())
        self.assertFalse((d / "dhat-blocks.svg").exists())
        self.assertTrue((d / "dhat-bytes-collapsed.txt").exists())
        self.assertEqual(r.returncode, 0)
        self.assertEqual(r.stderr.count("inferno-flamegraph"), 1)  # one hint line, not a failure

    def test_svg_step_invoked_via_path_stub_when_present(self):
        d = self.tmp()
        log = d / "inferno.log"
        (d / "dhat-heap.json").write_text(json.dumps(dhat_doc()))
        run_cat(d, env={"PATH": self.stub_path(log, True)})
        self.assertTrue((d / "dhat-bytes.svg").exists())
        self.assertTrue((d / "dhat-blocks.svg").exists())
        calls = log.read_text()
        self.assertEqual(calls.count("--title"), 2)
        text = (d / "catalog.md").read_text()
        self.assertIn("dhat-bytes.svg)", text)
        self.assertIn("dhat-blocks.svg)", text)


if __name__ == "__main__":
    unittest.main()
