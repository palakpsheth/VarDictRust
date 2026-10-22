#!/usr/bin/env python3
"""harness/tests/test_golden_store.py — golden baseline store (issue #31). stdlib only."""
import importlib.util
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

HARNESS = pathlib.Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("golden_store", HARNESS / "golden_store.py")
gs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gs)
CLI = [sys.executable, str(HARNESS / "golden_store.py")]
CLEAN_ENV = ("GOLDEN_MODE", "GOLDEN_TIME_OUT", "GOLDEN_RUN_ID", "AF_THR", "CPP_CHUNK", "GOLDEN_REHASH",
             "GOLDEN_ORACLE_CMD", "GOLDEN_CPP_CMD", "LC_ALL", "LANG", "STUB_EXIT")
HAS_TIME = os.path.exists(gs.TIME_BIN)


class StoreCase(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp(prefix="gs-test-"))
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.store = self.tmp / "store"
        self.env = {k: v for k, v in os.environ.items() if k not in CLEAN_ENV}
        self.env.update(GOLDEN_STORE_DIR=str(self.store), COUNT=str(self.tmp / "count"))
        self.ref = self.write("ref.fa", b">c\nACGT\n")
        self._saved = {k: os.environ.get(k) for k in ("GOLDEN_STORE_DIR", "COUNT") + CLEAN_ENV}
        for k in CLEAN_ENV:
            os.environ.pop(k, None)
        os.environ["GOLDEN_STORE_DIR"] = str(self.store)
        os.environ["COUNT"] = str(self.tmp / "count")
        self.addCleanup(self._restore_env)

    def _restore_env(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def write(self, name, data):
        p = self.tmp / name
        p.write_bytes(data)
        return p

    def runs(self):
        c = self.tmp / "count"
        return len(c.read_text().splitlines()) if c.exists() else 0

    def cli(self, *args, env=None):
        return subprocess.run(CLI + [str(a) for a in args], env={**self.env, **(env or {})},
                              capture_output=True, text=True)

    def raw(self, payload="payload", extra=(), env=None, kind="java", out=None, script=None):
        out = out or self.tmp / "out.txt"
        script = script or f'echo run >> "$COUNT"; echo {payload}; echo err-line >&2'
        return self.cli("raw", "--kind", kind, "--label", "t", "--out", out, "--err-out", self.tmp / "err.txt",
                        "--input", self.ref, *extra, "--", "sh", "-c", script, env=env)

    def events(self):
        log = self.store / "stats.log"
        return [l.split("\t")[2] for l in log.read_text().splitlines()] if log.exists() else []

    def manifests(self):
        return sorted((self.store / "cells").glob("*/manifest.json"))


class TestKey(StoreCase):
    def test_key_ignores_input_path_but_not_thread_count(self):
        other = self.write("elsewhere.fa", self.ref.read_bytes())
        k1, _ = gs.cell_key("java", [self.ref], ["-G", str(self.ref), "-th", "4", "-f", "0.1"])
        k2, _ = gs.cell_key("java", [other], ["-G", str(other), "-th", "4", "-f", "0.1"])
        k3, _ = gs.cell_key("java", [other], ["-G", str(other), "-th", "8", "-f", "0.1"])
        self.assertEqual(k1, k2)
        self.assertNotEqual(k1, k3, "-th stays in the key: a time file measured at -th 4 is not valid at -th 8")

    def test_literal_inputs_keep_their_path_spelling(self):
        other = self.write("elsewhere.fa", self.ref.read_bytes())
        k1, _ = gs.cell_key("java", [self.ref], [str(self.ref)], literal=[self.ref])
        k2, _ = gs.cell_key("java", [other], [str(other)], literal=[other])
        self.assertNotEqual(k1, k2, "##reference= embeds the -G path, so the spelling is part of the cell")

    def test_key_changes_with_content_flags_kind_wrapper_and_env(self):
        base, _ = gs.cell_key("java", [self.ref], ["-f", "0.1"])
        self.assertNotEqual(base, gs.cell_key("java", [self.ref], ["-f", "0.2"])[0])
        self.assertNotEqual(base, gs.cell_key("cpp", [self.ref], ["-f", "0.1"])[0])
        self.assertNotEqual(base, gs.cell_key("java", [self.ref], ["-f", "0.1"], wrapper="deadbeef")[0])
        for var, val in (("AF_THR", "0.5"), ("LC_ALL", "de_DE.UTF-8"), ("LANG", "de_DE.UTF-8")):
            os.environ[var] = val
            self.assertNotEqual(base, gs.cell_key("java", [self.ref], ["-f", "0.1"])[0], var)
            del os.environ[var]
        self.write("ref.fa", b">c\nACGTACGT\n")
        self.assertNotEqual(base, gs.cell_key("java", [self.ref], ["-f", "0.1"])[0])

    def test_sha_cache_invalidates_on_content_change_even_with_forged_mtime(self):
        before = gs.sha256_file(self.ref)
        st = self.ref.stat()
        time.sleep(0.05)  # ctime is kernel-coarse (jiffy); a real later edit always lands in a new tick
        self.ref.write_bytes(b">c\nTTTT\n")  # same size, different content
        os.utime(self.ref, ns=(st.st_atime_ns, st.st_mtime_ns))  # mtime forged back; ctime cannot be
        self.assertNotEqual(gs.sha256_file(self.ref), before, "ctime in the cache key must defeat a forged mtime")

    def test_golden_rehash_ignores_the_cache(self):
        real = gs.sha256_file(self.ref)
        cache = self.store / ".shacache.json"
        cache.write_text(json.dumps({k: "0" * 64 for k in json.loads(cache.read_text())}))
        self.assertEqual(gs.sha256_file(self.ref), "0" * 64)
        os.environ["GOLDEN_REHASH"] = "1"
        self.assertEqual(gs.sha256_file(self.ref), real)

    def test_missing_input_is_a_clear_error(self):
        r = self.cli("raw", "--kind", "java", "--label", "t", "--out", self.tmp / "o", "--input",
                     self.tmp / "nope.bam", "--", "true")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("input not found", r.stderr)


class TestRoundTrip(StoreCase):
    def test_miss_then_hit_skips_the_command_and_replays_stderr(self):
        r1 = self.raw()
        self.assertEqual(r1.returncode, 0, r1.stderr)
        self.assertEqual(self.runs(), 1)
        first = (self.tmp / "out.txt").read_text()
        (self.tmp / "out.txt").unlink()
        (self.tmp / "err.txt").unlink()
        r2 = self.raw()
        self.assertEqual(r2.returncode, 0, r2.stderr)
        self.assertEqual(self.runs(), 1, "hit must not execute the command")
        self.assertEqual((self.tmp / "out.txt").read_text(), first)
        self.assertEqual((self.tmp / "err.txt").read_text(), "err-line\n")
        self.assertEqual(self.events(), ["MISS", "HIT"])

    def test_stderr_is_forwarded_when_no_err_out(self):
        script = 'echo run >> "$COUNT"; echo p; echo cpp_chunk=50000 >&2'
        for _ in range(2):
            r = self.cli("raw", "--kind", "cpp", "--label", "t", "--out", self.tmp / "o", "--input", self.ref,
                         "--", "sh", "-c", script)
            self.assertEqual(r.stderr, "cpp_chunk=50000\n")
        self.assertEqual(self.runs(), 1)

    def test_large_artifacts_are_gzipped_small_are_not(self):
        big = 'echo run >> "$COUNT"; head -c 1200000 /dev/zero | tr "\\0" a'
        self.assertEqual(self.raw(script=big).returncode, 0)
        files = sorted(p.name for p in self.store.glob("cells/*/0/*"))
        self.assertIn("out.gz", files)
        self.assertIn("err", files)
        (self.tmp / "out.txt").unlink()
        self.assertEqual(self.raw(script=big).returncode, 0)
        self.assertEqual(self.runs(), 1)
        self.assertEqual((self.tmp / "out.txt").stat().st_size, 1200000)

    def test_mode_off_never_reads_or_writes_the_store(self):
        for _ in range(2):
            self.assertEqual(self.raw(env={"GOLDEN_MODE": "off"}).returncode, 0)
        self.assertEqual(self.runs(), 2)
        self.assertEqual(self.manifests(), [])
        self.assertEqual(self.events(), ["OFF", "OFF"])

    def test_mode_regen_reruns_and_appends_an_entry(self):
        self.raw()
        self.raw(env={"GOLDEN_MODE": "regen"})
        self.assertEqual(self.runs(), 2)
        m = json.loads(self.manifests()[0].read_text())
        self.assertEqual(len(m["entries"]), 2)
        self.assertEqual(self.events(), ["MISS", "REGEN"])

    def test_failed_command_is_never_pinned_and_exit_code_propagates(self):
        r = self.raw(script='echo run >> "$COUNT"; exit 7')
        self.assertEqual(r.returncode, 7)
        self.assertEqual(self.manifests(), [])
        self.assertEqual(self.events(), ["FAIL"])

    def test_bad_golden_mode_is_rejected(self):
        r = self.raw(env={"GOLDEN_MODE": "sometimes"})
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("GOLDEN_MODE", r.stderr)


class TestHardening(StoreCase):
    def test_pin_failure_does_not_fail_the_run(self):
        out = self.tmp / "o.txt"
        orig = gs.put

        def boom(*a, **k):
            raise OSError(28, "No space left on device")

        gs.put = boom
        self.addCleanup(setattr, gs, "put", orig)
        rc = gs.run_cached("java", "t", [self.ref], ["x"], ["sh", "-c", 'echo run >> "$COUNT"; echo hi'], str(out),
                           capture_stdout=True)
        self.assertEqual(rc, 0, "an hours-long oracle run must not be reported as failed because the pin failed")
        self.assertEqual(out.read_text(), "hi\n")
        self.assertEqual(self.events(), ["PINFAIL"])
        self.assertEqual(self.manifests(), [])

    def test_input_changed_during_the_run_is_not_pinned(self):
        r = self.raw(script='echo run >> "$COUNT"; printf changed-and-longer > "$IN"; echo payload', env={"IN": str(self.ref)})
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.events(), ["UNSTABLE"])
        self.assertEqual(self.manifests(), [])
        self.assertIn("changed during the run", r.stderr)

    def test_empty_or_implausible_output_is_refused(self):
        vcf = self.tmp / "x.vcf"
        for script, want in (("echo garbage", "REFUSED"), ("true", "REFUSED"), ("printf '#CHROM\\tPOS\\n'", "MISS")):
            r = self.cli("raw", "--kind", "java", "--label", "t", "--out", vcf, "--input", self.ref, "--", "sh", "-c", script)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(self.events()[-1], want, script)
        self.assertEqual(len(self.manifests()), 1, "only the plausible VCF is pinned")

    def test_corrupt_manifest_fails_loud_and_is_not_overwritten(self):
        self.raw()
        mf = self.manifests()[0]
        mf.write_text("{ not json")
        r = self.raw()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("is corrupt", r.stderr)
        self.assertEqual(mf.read_text(), "{ not json")
        self.assertEqual(self.runs(), 1)

    @unittest.skipUnless(HAS_TIME, "needs /usr/bin/time")
    def test_time_file_is_not_replayed_across_hosts(self):
        t1, t2 = self.tmp / "t1", self.tmp / "t2"
        self.raw(env={"GOLDEN_TIME_OUT": str(t1)})
        mf = self.manifests()[0]
        m = json.loads(mf.read_text())
        m["entries"][0]["host"]["nproc"] = 99999  # pinned on a different machine
        mf.write_text(json.dumps(m))
        self.raw(env={"GOLDEN_TIME_OUT": str(t2)})
        self.assertEqual(self.runs(), 2, "a timing consumer must re-measure on its own host")
        self.assertEqual(self.events(), ["MISS", "MISS"])
        self.raw()  # no timing consumer: the output alone replays fine
        self.assertEqual(self.runs(), 2)
        self.assertEqual(self.events()[-1], "HIT")

    @unittest.skipUnless(HAS_TIME and shutil.which("pgrep"), "needs /usr/bin/time and pgrep")
    def test_sigterm_kills_the_command_group_and_cleans_temp_files(self):
        tdir = self.tmp / "tmpd"
        tdir.mkdir()
        p = subprocess.Popen(CLI + ["raw", "--kind", "java", "--label", "t", "--out", str(self.tmp / "o"),
                                    "--input", str(self.ref), "--", "sh", "-c",
                                    'echo started >> "$COUNT"; exec sleep 31.7'],
                             env={**self.env, "TMPDIR": str(tdir)}, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(100):
            if self.runs() == 1:
                break
            time.sleep(0.05)
        p.terminate()
        self.assertEqual(p.wait(timeout=30), 143)
        self.assertEqual(list(tdir.iterdir()), [], "gs-err-*/gs-time-* temp files must not leak")
        self.assertNotEqual(subprocess.run(["pgrep", "-f", "sleep 31.7"], capture_output=True).returncode, 0,
                            "the oracle process group must not survive its parent")

    def test_concurrent_runs_of_one_cell_serialise_their_pins(self):
        procs = [subprocess.Popen(CLI + ["raw", "--kind", "java", "--label", "t", "--out", str(self.tmp / f"o{i}"),
                                         "--input", str(self.ref), "--", "sh", "-c",
                                         'echo run >> "$COUNT"; sleep 0.5; echo payload'],
                                  env=self.env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) for i in range(4)]
        self.assertEqual([p.wait(timeout=60) for p in procs], [0, 0, 0, 0])
        self.assertEqual(self.events().count("MISS"), 4)
        m = json.loads(self.manifests()[0].read_text())
        self.assertEqual(len(m["entries"]), 4, "no lost manifest update")
        cell = self.manifests()[0].parent
        self.assertEqual(sorted(p.name for p in cell.iterdir() if p.name.isdigit()), ["0", "1", "2", "3"])
        self.assertFalse((cell / ".lock").exists())


if __name__ == "__main__":
    unittest.main()
