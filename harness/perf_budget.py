#!/usr/bin/env python3
"""Usage: perf_budget.py <evidence-name>...   Pre-registered (plan rev3, A9/OQ-13):
PASS iff for every named evidence: rust wall <= java wall AND rust peak RSS <= java peak RSS.
Inputs are the /usr/bin/time -v records compare.sh captured (both engines -th 8, same host, same courtesy).
Java side = vardict-java | R | var2vcf pipeline (the thing VarDictRust consolidates); RSS = largest
single child (getrusage RUSAGE_CHILDREN), i.e. the JVM vs vardictrust. Exit 0 = PASS.
Third arm (T3A-05): pinned vardictcpp 3 hc52dbad_0 is advisory only — per-evidence cpp wall/RSS
columns when that evidence carries cpp_time (same host/pins as the recorded evidence cmd,
LC_ALL=C time -v); a missing
cpp_time is reported as a loud SKIP and never affects the PASS verdict (gate stays java==rust)."""
import json, sys
def secs(w):
    p = [float(x) for x in w.split(":")]
    return p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) == 3 else 0)
rows, ok = {}, True
for n in sys.argv[1:]:
    d = json.load(open(f"benchmarks/evidence/{n}.json"))
    j, r = d.get("java_time"), d.get("rust_time")
    if not j or not r or None in (j["wall"], r["wall"], j["max_rss_kb"], r["max_rss_kb"]):
        sys.exit(f"perf_budget: {n} lacks java/rust timing")
    jw, rw, jm, rm = secs(j["wall"]), secs(r["wall"]), int(j["max_rss_kb"]), int(r["max_rss_kb"])
    p = rw <= jw and rm <= jm; ok = ok and p
    row = {"java_wall_s": jw, "rust_wall_s": rw, "java_rss_kb": jm, "rust_rss_kb": rm,
           "wall_ratio": round(rw / jw, 4), "rss_ratio": round(rm / jm, 4), "pass": p,
           "evidence_fingerprint": d["fingerprint"]}
    c = d.get("cpp_time")
    skip_r = str((d.get("extra") or {}).get("cpp_skip_reason") or "")
    non_exec = skip_r.startswith(("disabled via CPP_SKIP", "conda not on PATH",
                                  "vardictcpp not available", "vardictcpp unavailable",
                                  "engine failure"))
    if c and c.get("exit") == "0" and c.get("wall") is not None and c.get("max_rss_kb") is not None and not non_exec:
        cw, cm = secs(c["wall"]), int(c["max_rss_kb"])
        row.update({"cpp_wall_s": cw, "cpp_maxrss_kb": cm,
                    "cpp_wall_ratio": round(cw / jw, 4), "cpp_rss_ratio": round(cm / jm, 4)})
    else:
        row["cpp"] = "SKIP: no cpp_time in evidence"
    rows[n] = row
print(json.dumps({"budget": "rust<=java wall and peak RSS (A9)", "runs": rows, "pass": ok}, indent=2, sort_keys=True))
sys.exit(0 if ok else 1)
