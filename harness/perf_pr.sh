#!/usr/bin/env bash
# harness/perf_pr.sh — Same-host interleaved Base vs Head PR Performance Gate
# Measures wall time, CPU time, and peak RSS across interleaved trials of BASE and HEAD
# binaries and fails if both median(head)/median(base) > threshold and min(head)/min(base) > threshold,
# or if any invocation exits non-zero, or if output bytes diverge (unless --allow-vcf-diff).
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

BASE_REF=""
BASE_BIN="${BASE_BIN:-}"
HEAD_BIN="${HEAD_BIN:-}"
THRESHOLD="1.10"
RSS_THRESHOLD="1.10"
RUNS="3"
SLICES=""
OUT_JSON=""
ALLOW_DIFF="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base) BASE_REF="$2"; shift 2 ;;
    --base-bin) BASE_BIN="$2"; shift 2 ;;
    --head-bin) HEAD_BIN="$2"; shift 2 ;;
    --threshold) THRESHOLD="$2"; shift 2 ;;
    --rss-threshold) RSS_THRESHOLD="$2"; shift 2 ;;
    --runs) RUNS="$2"; shift 2 ;;
    --slices) SLICES="$2"; shift 2 ;;
    --out-json) OUT_JSON="$2"; shift 2 ;;
    --allow-vcf-diff) ALLOW_DIFF="1"; shift ;;
    *) echo "perf_pr.sh: unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [ "${PERF_PR_SKIP:-0}" = "1" ]; then
  echo "perf_pr.sh: skipped via PERF_PR_SKIP=1"
  exit 0
fi

mkdir -p scratch/perf-pr
if [ -z "$OUT_JSON" ]; then
  OUT_JSON="scratch/perf-pr/perf_pr_latest.json"
fi

CLEANUP_WT=""
cleanup() {
  if [ -n "$CLEANUP_WT" ] && [ -d "$CLEANUP_WT" ]; then
    git worktree remove --force "$CLEANUP_WT" >/dev/null 2>&1 || rm -rf "$CLEANUP_WT"
  fi
}
trap cleanup EXIT

if [ -z "$HEAD_BIN" ]; then
  cargo build -q --release
  HEAD_BIN="$ROOT_DIR/target/release/vardictrust"
fi

if [ -z "$BASE_BIN" ]; then
  if [ -z "$BASE_REF" ]; then
    BASE_REF="$(git merge-base origin/main HEAD 2>/dev/null || git rev-parse HEAD)"
  fi
  CLEANUP_WT="$ROOT_DIR/scratch/perf-pr/base-wt-$$"
  git worktree add --detach "$CLEANUP_WT" "$BASE_REF" >/dev/null 2>&1
  CARGO_TARGET_DIR="$ROOT_DIR/scratch/perf-pr/base-target" \
    cargo build -q --release --manifest-path "$CLEANUP_WT/Cargo.toml"
  BASE_BIN="$ROOT_DIR/scratch/perf-pr/base-target/release/vardictrust"
fi

if [ ! -x "$BASE_BIN" ]; then
  echo "perf_pr.sh: base binary not executable: $BASE_BIN" >&2
  exit 1
fi
if [ ! -x "$HEAD_BIN" ]; then
  echo "perf_pr.sh: head binary not executable: $HEAD_BIN" >&2
  exit 1
fi

python3 - "$BASE_BIN" "$HEAD_BIN" "$THRESHOLD" "$RSS_THRESHOLD" "$RUNS" "$SLICES" "$OUT_JSON" "$ALLOW_DIFF" <<'PY'
import hashlib
import json
import os
import pathlib
import statistics
import subprocess
import sys
import tempfile
import time

base_bin, head_bin, thr_s, rss_thr_s, runs_s, slices_s, out_json, allow_diff_s = sys.argv[1:9]
threshold = float(thr_s)
rss_threshold = float(rss_thr_s)
runs = max(1, int(runs_s))
allow_diff = allow_diff_s == "1"

all_slices = {
    "micro": [
        "-G", "fixtures/micro/micro.ref.fa",
        "-b", "fixtures/micro/micro.bam",
        "-N", "SAMPLE",
        "-c", "1", "-S", "2", "-E", "3", "-g", "4",
        "-f", "0.01",
        "fixtures/micro/micro.bed",
    ],
    "hg002-sub": [
        "-G", "scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta",
        "-b", "scratch/g4/hg002/HG002.chr20.bam",
        "-N", "HG002",
        "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-th", "4",
        "-f", "0.01",
        "scratch/profiles/baseline/hg002-sub/sub.bed",
    ],
    "hg008-sub": [
        "-G", "scratch/g4/ref/giabv3/GRCh38_GIABv3.fa",
        "-b", "scratch/g4/hg008/HG008-T.slice.bam|scratch/g4/hg008/HG008-N.slice.bam",
        "-N", "HG008-T|HG008-N",
        "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-th", "4",
        "-f", "0.01",
        "scratch/g4/hg008/engine.bed",
    ],
}


def slice_inputs_exist(args):
    for i, a in enumerate(args):
        if a in ("-G", "-b") and i + 1 < len(args):
            for p in args[i + 1].split("|"):
                if not pathlib.Path(p).exists():
                    return False
    bed = args[-1]
    return pathlib.Path(bed).exists()


if slices_s:
    selected = [s.strip() for s in slices_s.split(",") if s.strip()]
else:
    selected = [k for k in ("hg002-sub", "hg008-sub") if slice_inputs_exist(all_slices[k])]
    if not selected:
        selected = ["micro"]

results = {}
overall_pass = True


def run_once(binary, args, capture=False):
    t0 = time.perf_counter()
    if capture:
        proc = subprocess.Popen([binary, *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        wall = time.perf_counter() - t0
        return proc.returncode, wall, 0, hashlib.sha256(out).hexdigest(), err.decode("utf-8", errors="replace")
    else:
        with open(os.devnull, "wb") as devnull, tempfile.TemporaryFile() as stderr_tmp:
            proc = subprocess.Popen([binary, *args], stdout=devnull, stderr=stderr_tmp)
            try:
                _, status, rusage = os.wait4(proc.pid, 0)
                rc = os.waitstatus_to_exitcode(status) if hasattr(os, "waitstatus_to_exitcode") else (status >> 8)
                max_rss = int(getattr(rusage, "ru_maxrss", 0))
            except AttributeError:
                rc = proc.wait()
                max_rss = 0
            stderr_tmp.seek(0)
            err = stderr_tmp.read().decode("utf-8", errors="replace")
        wall = time.perf_counter() - t0
        return rc, wall, max_rss, "", err


print("| Slice | Base Median (s) | Head Median (s) | Med Ratio | Min Ratio | Output SHA256 | Verdict |")
print("| :--- | ---: | ---: | ---: | ---: | :--- | :--- |")

for sl in selected:
    if sl not in all_slices:
        sys.exit(f"perf_pr.sh: unknown slice '{sl}'")
    args = all_slices[sl]
    if not slice_inputs_exist(args):
        sys.exit(f"perf_pr.sh: missing input files for slice '{sl}'")

    # Verification pass (captures output digest & checks exit code)
    rc_b, _, _, sha_b, err_b = run_once(base_bin, args, capture=True)
    if rc_b != 0:
        sys.exit(f"perf_pr.sh: base binary exited {rc_b} on {sl}: {err_b}")
    rc_h, _, _, sha_h, err_h = run_once(head_bin, args, capture=True)
    if rc_h != 0:
        sys.exit(f"perf_pr.sh: head binary exited {rc_h} on {sl}: {err_h}")

    sha_match = (sha_b == sha_h)
    if not sha_match and not allow_diff:
        print(f"| {sl} | - | - | - | - | MISMATCH | FAIL (output diff) |")
        sys.exit(f"perf_pr.sh: output digest mismatch on {sl}: base={sha_b[:12]} head={sha_h[:12]}")

    base_walls, head_walls = [], []
    base_rss, head_rss = [], []
    for _ in range(runs):
        rc_b, wb, rb, _, err_b = run_once(base_bin, args, capture=False)
        if rc_b != 0:
            sys.exit(f"perf_pr.sh: base binary exited {rc_b} during timing on {sl}: {err_b}")
        rc_h, wh, rh, _, err_h = run_once(head_bin, args, capture=False)
        if rc_h != 0:
            sys.exit(f"perf_pr.sh: head binary exited {rc_h} during timing on {sl}: {err_h}")
        base_walls.append(wb)
        head_walls.append(wh)
        if rb > 0:
            base_rss.append(rb)
        if rh > 0:
            head_rss.append(rh)

    med_b = statistics.median(base_walls)
    med_h = statistics.median(head_walls)
    min_b = min(base_walls)
    min_h = min(head_walls)

    med_ratio = (med_h / med_b) if med_b > 0 else 1.0
    min_ratio = (min_h / min_b) if min_b > 0 else 1.0

    # Fail if BOTH median and min exceed threshold
    wall_regressed = (med_ratio > threshold) and (min_ratio > threshold)
    rss_regressed = False
    rss_ratio = 1.0
    if base_rss and head_rss:
        med_rb = statistics.median(base_rss)
        med_rh = statistics.median(head_rss)
        if med_rb > 0:
            rss_ratio = med_rh / med_rb
            # Only enforce RSS ratio when baseline RSS is meaningful (> 32 MiB)
            if med_rb > 32768 and rss_ratio > rss_threshold and (min(head_rss) / max(1, min(base_rss))) > rss_threshold:
                rss_regressed = True

    sl_pass = (not wall_regressed) and (not rss_regressed) and (sha_match or allow_diff)
    overall_pass = overall_pass and sl_pass
    verdict = "PASS" if sl_pass else "FAIL"
    print(f"| {sl} | {med_b:.4f} | {med_h:.4f} | {med_ratio:.3f}x | {min_ratio:.3f}x | {sha_h[:12]} | {verdict} |")

    results[sl] = {
        "base_median_s": round(med_b, 6),
        "head_median_s": round(med_h, 6),
        "base_min_s": round(min_b, 6),
        "head_min_s": round(min_h, 6),
        "median_wall_ratio": round(med_ratio, 4),
        "min_wall_ratio": round(min_ratio, 4),
        "rss_ratio": round(rss_ratio, 4),
        "output_sha256_match": sha_match,
        "pass": sl_pass,
    }

summary = {
    "threshold": threshold,
    "runs": runs,
    "slices": results,
    "pass": overall_pass,
}
pathlib.Path(out_json).parent.mkdir(parents=True, exist_ok=True)
pathlib.Path(out_json).write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
sys.exit(0 if overall_pass else 1)
PY
