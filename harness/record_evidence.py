#!/usr/bin/env python3
"""Usage: record_evidence.py <name> <result> <dir|-> -- <cmd...>
Writes benchmarks/evidence/<name>.json; reads optional <dir>/compare.txt, java.time, rust.time, extra.json,
plus the first 5 lines of <dir>/compare_cpp.txt into extra.cpp_head (advisory vardictcpp third arm),
plus <dir>/cpp.time into cpp_time (same tv() shape as java/rust for perf_budget.py)."""
import datetime, json, pathlib, re, subprocess, sys
name, result, d = sys.argv[1], sys.argv[2], sys.argv[3]
cmd = sys.argv[sys.argv.index("--") + 1:]
keep = None if d == "-" else pathlib.Path(d)
def tv(p):
    if keep is None or not (keep / p).exists(): return None
    t = (keep / p).read_text()
    g = lambda k: (re.search(rf"{re.escape(k)}: (.+)", t) or [None, None])[1]
    return {"wall": g("Elapsed (wall clock) time (h:mm:ss or m:ss)"),
            "max_rss_kb": g("Maximum resident set size (kbytes)"), "exit": g("Exit status")}
cmp_txt = (keep / "compare.txt").read_text() if keep and (keep / "compare.txt").exists() else ""
m = re.search(r"^MATCH \((\d+) variants\)", cmp_txt, re.M)
extra_in = json.loads((keep / "extra.json").read_text()) if keep and (keep / "extra.json").exists() else None
cpp_tv = tv("cpp.time")
PRE_EXEC_SKIP = ("disabled via CPP_SKIP", "conda not on PATH", "vardictcpp not available", "vardictcpp unavailable")
skip_r = str((extra_in or {}).get("cpp_skip_reason") or "") if isinstance(extra_in, dict) else ""
if cpp_tv and (cpp_tv.get("exit") != "0" or skip_r.startswith((*PRE_EXEC_SKIP, "engine failure"))):
    cpp_tv = None
rec = {"name": name, "result": result, "n_variants": int(m.group(1)) if m else None,
       "fingerprint": subprocess.check_output(["bash", "harness/engine_fingerprint.sh"], text=True).strip(),
       "head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
       "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
       "cmd": cmd, "java_time": tv("java.time"), "rust_time": tv("rust.time"),
       "cpp_time": cpp_tv,
       "compare_head": cmp_txt.splitlines()[:25]}
if extra_in is not None: rec["extra"] = extra_in
if skip_r.startswith(PRE_EXEC_SKIP):
    extra = rec.get("extra")
    if not isinstance(extra, dict):
        extra = rec["extra"] = {}
    extra["cpp_head"] = [f"CPP: SKIP: {skip_r}"]
elif keep is not None and (keep / "compare_cpp.txt").exists():
    cpp_head = (keep / "compare_cpp.txt").read_text().splitlines()[:5]
    extra = rec.get("extra")
    if not isinstance(extra, dict):
        extra = rec["extra"] = {}
    extra["cpp_head"] = cpp_head
out = pathlib.Path("benchmarks/evidence") / f"{name}.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n"); print(out)
