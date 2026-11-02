#!/usr/bin/env python3
"""Usage: perf_budget.py [--baseline DIR] [--current DIR] [--tolerance PCT] [--vs-cpp SPEC] [evidence-name...]
Pre-registered (plan rev3, A9/OQ-13):
PASS iff for every named evidence:
  1. rust_exit == 0 (and rust_time.exit == "0", java_time.exit == "0" when present).
  2. result in ("MATCH", "PASS", "EQUAL") when present (and strictly "MATCH" for
     extended_coverage, multimodal_slices, and cram_parity suites).
  3. rust wall <= java wall AND rust peak RSS <= java peak RSS (when java_time & rust_time present).
  4. For cram_parity: when both rust_bam_wall_s and rust_cram_wall_s are positive,
     rust_cram_wall_s / rust_bam_wall_s <= 1.50.
  5. When --baseline DIR is provided, current wall and RSS must not exceed baseline by more than
     --tolerance PCT (default 10.0%).
  6. When --vs-cpp SPEC is provided, rust_wall_s / cpp_wall_s <= threshold (excluding amplicon floor).
Exit 0 = PASS, 1 = FAIL."""

import argparse
import json
import pathlib
import sys


SUITE_NAMES = frozenset([
    "extended_coverage",
    "extended-coverage-tier2",
    "multimodal_slices",
    "multimodal-slices-tier2",
    "cram_parity",
    "hg008-cram-twin",
])


def secs(w):
    if isinstance(w, (int, float)):
        return float(w)
    p = [float(x) for x in str(w).split(":")]
    return p[-1] + (60.0 * p[-2] if len(p) >= 2 else 0.0) + (3600.0 * p[-3] if len(p) == 3 else 0.0)


def extract_rust_wall(d):
    if d.get("rust_wall_s") is not None:
        return float(d["rust_wall_s"])
    r = d.get("rust_time")
    if isinstance(r, dict) and r.get("wall") is not None:
        return secs(r["wall"])
    if d.get("rust_cram_wall_s") is not None:
        return float(d["rust_cram_wall_s"])
    extra = d.get("extra") or {}
    if extra.get("rust_cram_wall_s") is not None:
        return float(extra["rust_cram_wall_s"])
    return None


def extract_rust_rss(d):
    if d.get("rust_rss_kb") is not None:
        return int(d["rust_rss_kb"])
    r = d.get("rust_time")
    if isinstance(r, dict) and r.get("max_rss_kb") is not None:
        return int(r["max_rss_kb"])
    return None


def parse_vs_cpp(spec):
    if not spec:
        return None, {}
    parts = spec.split(":")
    default_ratio = float(parts[0]) if parts[0] else None
    overrides = {}
    for item in parts[1:]:
        for kv in item.split(","):
            if "=" in kv:
                k, v = kv.split("=", 1)
                overrides[k.strip()] = float(v.strip())
    return default_ratio, overrides


def resolve_target_path(name, curr_dir):
    p = pathlib.Path(name)
    if p.is_file():
        return p.stem if p.suffix == ".json" else p.name, p
    cand1 = curr_dir / (name if name.endswith(".json") else f"{name}.json")
    if cand1.is_file():
        return pathlib.Path(name).stem if name.endswith(".json") else name, cand1
    cand2 = pathlib.Path("benchmarks/evidence") / (name if name.endswith(".json") else f"{name}.json")
    if cand2.is_file():
        return pathlib.Path(name).stem if name.endswith(".json") else name, cand2
    return name, cand1


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--baseline", type=str, default=None)
    parser.add_argument("--current", type=str, default=None)
    parser.add_argument("--tolerance", type=float, default=10.0)
    parser.add_argument("--vs-cpp", dest="vs_cpp", type=str, default=None)
    parser.add_argument("names", nargs="*")
    args = parser.parse_args()

    curr_dir = pathlib.Path(args.current) if args.current else pathlib.Path("benchmarks/evidence")
    base_dir = pathlib.Path(args.baseline) if args.baseline else None
    cpp_default, cpp_overrides = parse_vs_cpp(args.vs_cpp)

    target_list = list(args.names)
    if args.current:
        for auto_suite in ("extended_coverage", "multimodal_slices", "cram_parity"):
            if (curr_dir / f"{auto_suite}.json").is_file():
                if auto_suite not in target_list and f"{auto_suite}.json" not in target_list:
                    target_list.append(auto_suite)
    if not target_list:
        if curr_dir.is_dir():
            for f in sorted(curr_dir.glob("*.json")):
                if f.stem not in ("perf-budget", "perf_budget"):
                    target_list.append(f.stem)
        if not target_list:
            sys.exit("perf_budget: no evidence targets specified or found")

    rows = {}
    ok = True
    max_mult = 1.0 + (args.tolerance / 100.0)

    for raw_name in target_list:
        n, path = resolve_target_path(raw_name, curr_dir)
        if not path.is_file():
            sys.exit(f"perf_budget: {n} file not found at {path}")
        with open(path, "r", encoding="utf-8") as fh:
            d = json.load(fh)
        j, r = d.get("java_time"), d.get("rust_time")
        extra = d.get("extra") or {}

        is_suite_schema = (
            n in SUITE_NAMES
            or "rust_bam_wall_s" in d
            or "rust_cram_wall_s" in d
            or "rust_bam_wall_s" in extra
            or "rust_cram_wall_s" in extra
            or "modalities" in d
            or "slices" in d
            or (j is None and r is None and ("rust_wall_s" in d or ("result" in d and "rust_exit" in d)))
        )

        p = True

        # 1. Check exit status
        if d.get("rust_exit") is not None and int(d["rust_exit"]) != 0:
            p = False
        if isinstance(r, dict) and str(r.get("exit", "0")) != "0":
            p = False
        if isinstance(j, dict) and str(j.get("exit", "0")) != "0":
            p = False

        # 2. Check result status
        res_val = d.get("result")
        if n in SUITE_NAMES:
            if res_val != "MATCH":
                p = False
        elif res_val is not None:
            if res_val not in ("MATCH", "PASS", "EQUAL"):
                p = False

        if isinstance(d.get("modalities"), dict):
            if int(d.get("diff_count", 0)) != 0:
                p = False
            for mod_info in d["modalities"].values():
                if isinstance(mod_info, dict):
                    if mod_info.get("result") != "MATCH" or int(mod_info.get("rust_exit", 0)) != 0:
                        p = False
        if isinstance(d.get("slices"), dict):
            for sl_info in d["slices"].values():
                if isinstance(sl_info, dict):
                    if sl_info.get("result", "MATCH") != "MATCH" or int(sl_info.get("rust_exit", 0)) != 0:
                        p = False
        if d.get("total_slices") is not None and d.get("slices_passed") is not None:
            if int(d["slices_passed"]) != int(d["total_slices"]):
                p = False

        # 3. Check Java vs Rust wall & RSS when standard benchmark or java_time/rust_time present
        if not is_suite_schema or (j is not None and r is not None):
            if not j or not r or None in (j.get("wall"), r.get("wall"), j.get("max_rss_kb"), r.get("max_rss_kb")):
                sys.exit(f"perf_budget: {n} lacks java/rust timing")
            jw, rw, jm, rm = secs(j["wall"]), secs(r["wall"]), int(j["max_rss_kb"]), int(r["max_rss_kb"])
            p = p and (rw <= jw and rm <= jm)
            row = {
                "java_wall_s": jw,
                "rust_wall_s": rw,
                "java_rss_kb": jm,
                "rust_rss_kb": rm,
                "wall_ratio": round(rw / jw, 4),
                "rss_ratio": round(rm / jm, 4),
                "pass": p,
                "evidence_fingerprint": d.get("fingerprint", ""),
            }
            c = d.get("cpp_time")
            skip_r = str(extra.get("cpp_skip_reason") or "")
            non_exec = skip_r.startswith((
                "disabled via CPP_SKIP",
                "conda not on PATH",
                "vardictcpp not available",
                "vardictcpp unavailable",
                "engine failure",
            ))
            if (
                c
                and c.get("exit") == "0"
                and c.get("wall") is not None
                and c.get("max_rss_kb") is not None
                and not non_exec
            ):
                cw, cm = secs(c["wall"]), int(c["max_rss_kb"])
                row.update({
                    "cpp_wall_s": cw,
                    "cpp_maxrss_kb": cm,
                    "cpp_wall_ratio": round(cw / jw, 4),
                    "cpp_rss_ratio": round(cm / jm, 4),
                })
                cpp_lim = cpp_overrides.get(n, cpp_default)
                if cpp_lim is not None and not n.startswith("amplicon") and rw >= 0.05 and cw > 0:
                    row["rust_vs_cpp_wall_ratio"] = round(rw / cw, 4)
                    if (rw / cw) > cpp_lim:
                        p = False
            else:
                row["cpp"] = "SKIP: no cpp_time in evidence"
        else:
            rw = extract_rust_wall(d)
            rm = extract_rust_rss(d)
            row = {
                "evidence_fingerprint": d.get("fingerprint", ""),
                "result": res_val,
                "rust_exit": int(d.get("rust_exit", 0)) if d.get("rust_exit") is not None else 0,
                "pass": p,
            }
            if rw is not None:
                row["rust_wall_s"] = rw
            if rm is not None:
                row["rust_rss_kb"] = rm

        # 4. CRAM parity wall-time ratio check (rust_cram_wall_s / rust_bam_wall_s <= 1.50)
        bam_w = d.get("rust_bam_wall_s", extra.get("rust_bam_wall_s"))
        cram_w = d.get("rust_cram_wall_s", extra.get("rust_cram_wall_s"))
        if bam_w is not None and cram_w is not None:
            bam_wf, cram_wf = float(bam_w), float(cram_w)
            row["rust_bam_wall_s"] = bam_wf
            row["rust_cram_wall_s"] = cram_wf
            if bam_wf > 0.0 and cram_wf > 0.0:
                cram_ratio = cram_wf / bam_wf
                row["cram_bam_wall_ratio"] = round(cram_ratio, 4)
                if cram_ratio > 1.50:
                    p = False

        # 5. Baseline vs current relative regression check
        if base_dir is not None:
            _, bpath = resolve_target_path(n, base_dir)
            if bpath.is_file():
                with open(bpath, "r", encoding="utf-8") as bfh:
                    bd = json.load(bfh)
                base_rw = extract_rust_wall(bd)
                curr_rw = extract_rust_wall(d)
                if base_rw is not None and base_rw > 0.0 and curr_rw is not None:
                    rel_w = curr_rw / base_rw
                    row["baseline_wall_s"] = base_rw
                    row["baseline_wall_ratio"] = round(rel_w, 4)
                    if rel_w > max_mult:
                        p = False
                base_rm = extract_rust_rss(bd)
                curr_rm = extract_rust_rss(d)
                if base_rm is not None and base_rm > 0 and curr_rm is not None:
                    rel_m = curr_rm / base_rm
                    row["baseline_rss_kb"] = base_rm
                    row["baseline_rss_ratio"] = round(rel_m, 4)
                    if rel_m > max_mult:
                        p = False

        row["pass"] = p
        ok = ok and p
        rows[n] = row

    print(json.dumps({"budget": "rust<=java wall and peak RSS (A9)", "runs": rows, "pass": ok}, indent=2, sort_keys=True))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
