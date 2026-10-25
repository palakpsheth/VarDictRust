#!/usr/bin/env bash
# Usage: scale_bench.sh [--no-cpp] <ref> <bam> <bed> <sample> [outdir]   (outdir default: scratch/g4/scale)
# Certifying run (≈7 h with 3 runs; start between 05:15 and ~20:30, outside the 03:40-05:15 maintenance window):
#   EST_MIN=450 bash harness/courtesy.sh bash harness/scale_bench.sh --no-cpp <ref> <bam> <bed> <sample> scratch/g4/scale
# Pre-registered (plan rev2, recalibrated T3D-02 signed): PASS iff wall(th1)/wall(th8) >= bar(P) = 0.6*min(8,P) = 2.4
# (P=4 distinct physical cores under the taskset -c 0-7 pin, §2 lscpu count, recorded as P in extra.json;
# RSS leg host-independent, unchanged <= 1.5x) AND outputs identical across thread counts. Bar frozen after sign-off.
# Third arm (T3A-05): pinned vardictcpp 3 hc52dbad_0 runs beside each rust loop
# (advisory only; gate above unchanged). Probe facts (Task 1): THREAD -th N,
# CHUNK --chunk (cpp default; no explicit pin), VCF via stdout redirect
# (-o is numeric in cpp, not a file), -R takes only chr:start-end (never a BED path).
# Issue #31 protocol (signed §4): every cell is the MEDIAN of SCALE_RUNS (default 3) timed runs after one
# untimed warmup. The bench pins to SCALE_CPUS (default: first up-to-8 allowed CPUs, 0-7 on unpinned hosts;
# SCALE_CPUS= disables the pin); P is derived from the effective CPU set applied under taskset. extra.json["certifying"]
# needs: pinned (taskset enforced), niced by courtesy.sh, >=3 runs, P>=4 (the frozen bar is 0.6*P=2.4), and a clean
# engine tree (perf_fingerprint).
set -euo pipefail
NO_CPP=0
if [ "${1:-}" = "--no-cpp" ]; then NO_CPP=1; shift; fi
if [ "$NO_CPP" = 1 ]; then export CPP_SKIP=1; fi
REF=$1; BAM=$2; BED=$3; SM=$4; OUT=${5:-scratch/g4/scale}; BIN=${VARDICT_BIN:-target/release/vardictrust}; mkdir -p "$OUT"
rm -f "$OUT/extra.json"   # a failed run must not leave the previous run's verdict next to fresh partial logs
RUNS=${SCALE_RUNS:-3}
VARDICTCPP_BIN=${VARDICTCPP_BIN:-conda run -n vardict-oracle vardictcpp}
# courtesy gate: harness/courtesy.sh exports COURTESY_ACTIVE=1 and runs us at nice 19; both must hold to certify.
COURTESY=false
if [ "${COURTESY_ACTIVE:-0}" = "1" ] && [ "$(ps -o ni= -p $$ | tr -d ' ')" -ge 10 ]; then COURTESY=true; fi
if [ "$COURTESY" = false ] && [ "${SCALE_ALLOW_UNCOURTEOUS:-0}" != "1" ]; then
  echo "scale_bench: run via 'EST_MIN=<min> bash harness/courtesy.sh bash harness/scale_bench.sh ...' (SCALE_ALLOW_UNCOURTEOUS=1 allows a non-certifying advisory run)" >&2
  exit 2
fi
# CPU set (default: first up-to-8 permitted CPUs from current affinity, 0-7 on standard hosts). Plain lists/ranges
# only ("0-7", "0-3,8"); P counts distinct physical (core,socket) in the effective taskset affinity; pinned = an
# explicit taskset mask is enforced by the benchmark (SCALE_CPUS= opts out and sets pinned=false).
ALLOWED_AFF=$(taskset -cp $$ 2>/dev/null | sed 's/.*: //' || true)
[ -n "$ALLOWED_AFF" ] || ALLOWED_AFF="0-$(($(nproc --all) - 1))"
if [ -z "${SCALE_CPUS+x}" ]; then
  SCALE_CPUS=$(python3 - "$ALLOWED_AFF" <<'PY'
import sys
cpus = []
for part in sys.argv[1].split(","):
    a, _, b = part.partition("-")
    cpus.extend(range(int(a), int(b or a) + 1))
sel = sorted(set(cpus))[:8]
ranges = []
for c in sel:
    if ranges and c == ranges[-1][1] + 1:
        ranges[-1][1] = c
    else:
        ranges.append([c, c])
print(",".join(f"{a}-{b}" if a != b else str(a) for a, b in ranges))
PY
  )
fi
TASKSET=""
PINNED=false
if [ -n "$SCALE_CPUS" ]; then
  [[ "$SCALE_CPUS" =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]] || { echo "scale_bench: SCALE_CPUS must be a plain list/range like 0-7 or 0-3,8 (got '$SCALE_CPUS')" >&2; exit 2; }
  command -v taskset >/dev/null 2>&1 || { echo "scale_bench: taskset is required for SCALE_CPUS" >&2; exit 2; }
  TASKSET="taskset -c $SCALE_CPUS"
  EFF=$($TASKSET bash -c 'taskset -cp $$' 2>/dev/null | sed 's/.*: //' || true)
  [ -n "$EFF" ] || { echo "scale_bench: failed to apply taskset -c $SCALE_CPUS" >&2; exit 1; }
  PINNED=true
else
  EFF="$ALLOWED_AFF"
fi
command -v lscpu >/dev/null 2>&1 || { echo "scale_bench: lscpu is required to derive P" >&2; exit 2; }
read -r P_PHYS CPU_COUNT < <(python3 - "$EFF" <<'PY'
import subprocess, sys
cpus = set()
for part in sys.argv[1].split(","):
    a, _, b = part.partition("-")
    cpus.update(range(int(a), int(b or a) + 1))
rows = subprocess.check_output(["lscpu", "-p=CPU,CORE,SOCKET"], text=True).splitlines()
phys = {(r.split(",")[1], r.split(",")[2]) for r in rows if not r.startswith("#") and int(r.split(",")[0]) in cpus}
print(len(phys), len(cpus))
PY
)
[ -n "${P_PHYS:-}" ] || { echo "scale_bench: could not derive P from lscpu" >&2; exit 2; }
PERF_FP=$(bash harness/engine_fingerprint.sh --perf 2>/dev/null || echo dirty)   # fail fast, not after hours
GOVERNOR=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)
LOAD_START=$(cut -d' ' -f1-3 /proc/loadavg)
# Full-span size from the reference index (no hardcoded chr1 length).
CHR1=$(awk '$1=="chr1"||$1=="1"{print $1; exit}' "$REF.fai")
SPAN=$(awk '$1=="chr1"||$1=="1"{print $2; exit}' "$REF.fai")
[ -n "$SPAN" ] || { echo "scale_bench: no chr1/1 contig in $REF.fai" >&2; exit 2; }
SIZES=""
for s in 1000000 10000000 100000000; do if [ "$s" -lt "$SPAN" ]; then SIZES="$SIZES $s"; fi; done
SIZES="$SIZES $SPAN"
# NOTE(issue-6): P counts distinct (core,socket) under the effective CPU set — meaningful only when PINNED (unpinned P/bar advisory-only).
# cpp availability probe: a missing pin or CPP_SKIP=1 degrades to loud SKIP, never masks the gate.
CPP_OK="false"; CPP_SKIP_MSG="vardictcpp unavailable ($VARDICTCPP_BIN --version failed)"
if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
  CPP_SKIP_MSG="disabled via CPP_SKIP=1"
elif $VARDICTCPP_BIN --version >/dev/null 2>&1; then
  CPP_OK="true"
fi
# cell <name> <stdout-file> <cmd...>: RUNS timed executions; per-run `time -v` logs <OUT>/<name>.run<i>.time are
# preserved (also on failure); stdout of the last run is kept. Returns 1 naming the failed cell/run.
cell() {
  local name=$1 outf=$2 i rc; shift 2
  for ((i = 1; i <= RUNS; i++)); do
    rc=0
    LC_ALL=C /usr/bin/time -v -o "$OUT/$name.run$i.time" "$@" > "$outf" || rc=$?
    if [ "$rc" -ne 0 ]; then
      echo "scale_bench: CELL FAILED: $name run $i (exit $rc); partial time log: $OUT/$name.run$i.time" >&2
      return 1
    fi
  done
}
cpp_cell() {
  local name=$1 label=$2 i rc; shift 2
  for ((i = 1; i <= RUNS; i++)); do
    rc=0
    LC_ALL=C GOLDEN_TIME_OUT="$OUT/$name.run$i.time" $TASKSET \
      bash harness/golden_store.sh cpp "$label" "$REF" "$BAM" "$@" > /dev/null || rc=$?
    if [ "$rc" -ne 0 ]; then
      echo "scale_bench: CELL FAILED: $name run $i (exit $rc); partial time log: $OUT/$name.run$i.time" >&2
      return 1
    fi
  done
}
# warmup: one untimed pass over the smallest span (page cache, allocator, CPU frequency ramp).
$TASKSET "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -th 8 -R "$CHR1:1-1000000" > /dev/null \
  || { echo "scale_bench: warmup failed" >&2; exit 1; }
for th in 1 2 4 8; do
  cell "th$th" "$OUT/th$th.vcf" $TASKSET "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th "$th" "$BED" || exit 1
  # BED-capable path reuses run_cpp.sh via golden_store.sh (BED positional, java-mirror; pins
  # -c 1 -S 2 -E 3 -g 4 identical to the rust line, adds -vcf for stdout VCF).
  if [ "$CPP_OK" = "true" ]; then
    cpp_cell "cpp_th$th" "scale-th$th" "$BED" "$OUT/cpp_th$th.vcf" "$SM" -th "$th" -z \
      || { rm -f "$OUT"/cpp_th$th.run*.time "$OUT/cpp_th$th.vcf"; echo "CPP: SKIP: vardictcpp th$th failed" >&2; }
  else
    rm -f "$OUT"/cpp_th$th.run*.time "$OUT/cpp_th$th.vcf"
    echo "CPP: SKIP: $CPP_SKIP_MSG" >&2
  fi
done
for th in 2 4 8; do cmp <(grep -v '^##' "$OUT/th1.vcf") <(grep -v '^##' "$OUT/th$th.vcf"); done
for sz in $SIZES; do
  cell "R$sz" "$OUT/R$sz.vcf" $TASKSET "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -th 8 -R "$CHR1:1-$sz" || exit 1
  # cpp region arm: both cpp arms go through golden_store.sh cpp (run_cpp.sh); -R is mirrored as the equivalent 0-based BED [0,sz).
  if [ "$CPP_OK" = "true" ]; then
    printf '%s\t0\t%s\tR%s\n' "$CHR1" "$sz" "$sz" > "$OUT/R$sz.bed"
    cpp_cell "cpp_R$sz" "scale-R$sz" "$OUT/R$sz.bed" "$OUT/cpp_R$sz.vcf" "$SM" -th 8 -z \
      || { rm -f "$OUT"/cpp_R$sz.run*.time "$OUT/cpp_R$sz.vcf"; echo "CPP: SKIP: vardictcpp R$sz failed" >&2; }
  else
    rm -f "$OUT"/cpp_R$sz.run*.time "$OUT/cpp_R$sz.vcf"
    echo "CPP: SKIP: $CPP_SKIP_MSG" >&2
  fi
done
LOAD_END=$(cut -d' ' -f1-3 /proc/loadavg)
SIZES="$SIZES" SPAN="$SPAN" RUNS="$RUNS" COURTESY="$COURTESY" GOVERNOR="$GOVERNOR" LOAD_START="$LOAD_START" LOAD_END="$LOAD_END" EFF="$EFF" PERF_FP="$PERF_FP" \
python3 - "$OUT" "$PINNED" "$P_PHYS" > "$OUT/extra.json" <<'PY'
import json, os, pathlib, re, statistics, sys
o = pathlib.Path(sys.argv[1]); pinned = sys.argv[2] == "true"; P = int(sys.argv[3])
runs = int(os.environ["RUNS"]); sizes = [int(s) for s in os.environ["SIZES"].split()]; span = int(os.environ["SPAN"])
MIN_P = 4  # pre-registered: the frozen bar 2.4 = 0.6*P is defined at P=4; fewer physical cores cannot certify it
def get(f, k): return re.search(rf"{re.escape(k)}: (.+)", (o / f).read_text()).group(1)
def secs(w):
    p = [float(x) for x in w.split(":")]
    return p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) == 3 else 0)
WALL, RSS = "Elapsed (wall clock) time (h:mm:ss or m:ss)", "Maximum resident set size (kbytes)"
def med(prefix, key, conv, strict=True):
    vals = []
    for i in range(1, runs + 1):
        f = o / f"{prefix}.run{i}.time"
        if not f.exists(): return None
        txt = f.read_text()
        if not strict and not re.search(r"Exit status:\s+0\b", txt): return None
        vals.append(conv(get(f"{prefix}.run{i}.time", key)))
    return statistics.median(vals)
W = {t: round(med(f"th{t}", WALL, secs), 2) for t in (1, 2, 4, 8)}
M = {s: int(med(f"R{s}", RSS, int)) for s in sizes}
CW = {t: med(f"cpp_th{t}", WALL, secs, strict=False) for t in (1, 2, 4, 8)}
CM = {s: med(f"cpp_R{s}", RSS, int, strict=False) for s in sizes}
sp = W[1] / W[8] if W[8] > 0 else None
rr = M[span] / M[10000000] if 10000000 in M else None
bar = 0.6 * min(8, P)
perf = os.environ["PERF_FP"]
print(json.dumps({"wall_s": W, "maxrss_kb": M, "speedup_8": sp, "rss_ratio_chr1_vs_10mb": rr,
                  "P": P, "bar": bar,
                  "pass": sp is not None and sp >= bar and rr is not None and rr <= 1.5,
                  "cpp_wall_s": CW, "cpp_maxrss_kb": CM, "pinned": pinned,
                  "runs": runs, "cpus": os.environ["EFF"], "courtesy": os.environ["COURTESY"] == "true",
                  "governor": os.environ["GOVERNOR"], "loadavg_start": os.environ["LOAD_START"],
                  "loadavg_end": os.environ["LOAD_END"], "perf_fingerprint": perf,
                  "certifying": pinned and os.environ["COURTESY"] == "true" and runs >= 3 and P >= MIN_P and perf != "dirty",
                  "note": ">8 threads not measured (host policy A3); each cell = median of `runs` timed runs"}, indent=2))
PY
