#!/usr/bin/env python3
"""Usage: record_evidence.py [--input <path>]... [--out-dir <dir>] <name> <result> <dir|-> -- <cmd...>
Writes benchmarks/evidence/<name>.json; reads optional <dir>/compare.txt, java.time, rust.time, extra.json,
plus the first 5 lines of <dir>/compare_cpp.txt into extra.cpp_head (advisory vardictcpp third arm),
plus <dir>/cpp.time into cpp_time (same tv() shape as java/rust for perf_budget.py),
and records SHA-256 + byte size of input files into inputs."""
import datetime
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys

args = sys.argv[1:]
sep_idx = args.index("--")
pre_args = args[:sep_idx]
cmd = args[sep_idx + 1 :]

explicit_inputs = []
out_dir = pathlib.Path("benchmarks/evidence")
pos_args = []
i = 0
while i < len(pre_args):
    a = pre_args[i]
    if a == "--input":
        explicit_inputs.append(pre_args[i + 1])
        i += 2
    elif a == "--out-dir":
        out_dir = pathlib.Path(pre_args[i + 1])
        i += 2
    else:
        pos_args.append(a)
        i += 1

name, result, d = pos_args[0], pos_args[1], pos_args[2]
keep = None if d == "-" else pathlib.Path(d)

INPUT_EXTS = (
    ".fa",
    ".fasta",
    ".fai",
    ".bam",
    ".bai",
    ".cram",
    ".crai",
    ".bed",
    ".vcf",
    ".vcf.gz",
    ".tbi",
    ".csi",
)


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def collect_inputs(cmd_tokens, extra_paths):
    candidates = list(extra_paths)
    for tok in cmd_tokens:
        if "=" in tok and not tok.startswith("-"):
            continue
        for part in tok.split("|"):
            part = part.strip()
            if not part or part.startswith("-"):
                continue
            if part.endswith(INPUT_EXTS):
                candidates.append(part)

    seen = set()
    recorded = []

    def add_path(p_str: str):
        p = pathlib.Path(p_str)
        if not p.is_file():
            return
        norm = str(p)
        if norm in seen:
            return
        seen.add(norm)
        st = p.stat()
        recorded.append(
            {
                "path": norm,
                "size_bytes": st.st_size,
                "sha256": sha256_file(p),
            }
        )

    for c in candidates:
        p = pathlib.Path(c)
        if not p.is_file():
            continue
        add_path(c)
        if c.endswith((".fa", ".fasta")):
            add_path(c + ".fai")
        elif c.endswith(".bam"):
            add_path(c + ".bai")
        elif c.endswith(".cram"):
            add_path(c + ".crai")
    return recorded


def tv(p):
    if keep is None or not (keep / p).exists():
        return None
    t = (keep / p).read_text()
    g = lambda k: (re.search(rf"{re.escape(k)}: (.+)", t) or [None, None])[1]
    return {
        "wall": g("Elapsed (wall clock) time (h:mm:ss or m:ss)"),
        "max_rss_kb": g("Maximum resident set size (kbytes)"),
        "exit": g("Exit status"),
    }


cmp_txt = (keep / "compare.txt").read_text() if keep and (keep / "compare.txt").exists() else ""
m = re.search(r"^MATCH \((\d+) variants\)", cmp_txt, re.M)
extra_in = (
    json.loads((keep / "extra.json").read_text())
    if keep and (keep / "extra.json").exists()
    else None
)
cpp_tv = tv("cpp.time")
PRE_EXEC_SKIP = (
    "disabled via CPP_SKIP",
    "conda not on PATH",
    "vardictcpp not available",
    "vardictcpp unavailable",
)
skip_r = str((extra_in or {}).get("cpp_skip_reason") or "") if isinstance(extra_in, dict) else ""
if cpp_tv and (cpp_tv.get("exit") != "0" or skip_r.startswith((*PRE_EXEC_SKIP, "engine failure"))):
    cpp_tv = None

fp_env = os.environ.get("RECORD_EVIDENCE_FP")
if fp_env:
    fp = fp_env.strip()
else:
    try:
        fp = subprocess.check_output(["bash", "harness/engine_fingerprint.sh"], text=True).strip()
    except subprocess.CalledProcessError:
        if out_dir != pathlib.Path("benchmarks/evidence"):
            fp = "0000000000000000"
        else:
            raise

inputs_list = collect_inputs(cmd, explicit_inputs)

rec = {
    "name": name,
    "result": result,
    "n_variants": int(m.group(1)) if m else None,
    "fingerprint": fp,
    "head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "cmd": cmd,
    "inputs": inputs_list,
    "java_time": tv("java.time"),
    "rust_time": tv("rust.time"),
    "cpp_time": cpp_tv,
    "compare_head": cmp_txt.splitlines()[:25],
}
if extra_in is not None:
    rec["extra"] = extra_in
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
out = out_dir / f"{name}.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
print(out)
