#!/usr/bin/env python3
"""Usage: record_evidence.py <name> <result> <dir|-> -- <cmd...>
Writes benchmarks/evidence/<name>.json; reads optional <dir>/compare.txt, java.time, rust.time, extra.json."""
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
rec = {"name": name, "result": result, "n_variants": int(m.group(1)) if m else None,
       "fingerprint": subprocess.check_output(["bash", "harness/engine_fingerprint.sh"], text=True).strip(),
       "head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
       "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
       "cmd": cmd, "java_time": tv("java.time"), "rust_time": tv("rust.time"),
       "compare_head": cmp_txt.splitlines()[:25]}
if keep and (keep / "extra.json").exists(): rec["extra"] = json.loads((keep / "extra.json").read_text())
out = pathlib.Path("benchmarks/evidence") / f"{name}.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n"); print(out)
