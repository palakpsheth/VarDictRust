#!/usr/bin/env bash
# Usage: scale_bench.sh <ref> <bam> <bed> <sample> <outdir>
# Pre-registered (plan rev2): PASS iff wall(th1)/wall(th8) >= 4.8 AND maxRSS(chr1 full) <= 1.5 x maxRSS(10 Mb)
# AND outputs identical across thread counts. Thresholds must not change after seeing data.
# Third arm (T3A-05): pinned vardictcpp 3 hc52dbad_0 runs beside each rust loop
# (advisory only; gate above unchanged). Probe facts (Task 1): THREAD -th N,
# CHUNK --chunk (cpp default; no explicit pin), VCF via stdout redirect
# (-o is numeric in cpp, not a file), -R takes only chr:start-end (never a BED path).
set -euo pipefail
REF=$1; BAM=$2; BED=$3; SM=$4; OUT=$5; BIN=${VARDICT_BIN:-target/release/vardictrust}; mkdir -p "$OUT"
VARDICTCPP_BIN=${VARDICTCPP_BIN:-conda run -n vardict-oracle vardictcpp}
# taskset guard: pin all cells to 0-7 when available; recorded in extra.pinned.
TASKSET=""; PINNED="false"
if command -v taskset >/dev/null 2>&1; then TASKSET="taskset -c 0-7"; PINNED="true"; fi
# cpp availability probe: a missing pin degrades to loud SKIP, never masks the gate.
CPP_OK="false"
if $VARDICTCPP_BIN --version >/dev/null 2>&1; then CPP_OK="true"; fi
for th in 1 2 4 8; do
  LC_ALL=C /usr/bin/time -v -o "$OUT/th$th.time" $TASKSET "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th "$th" "$BED" > "$OUT/th$th.vcf"
  # BED-capable path reuses run_cpp.sh (BED positional, java-mirror; pins
  # -c 1 -S 2 -E 3 -g 4 identical to the rust line, adds -vcf for stdout VCF).
  if [ "$CPP_OK" = "true" ]; then
    LC_ALL=C /usr/bin/time -v -o "$OUT/cpp_th$th.time" $TASKSET bash harness/run_cpp.sh "$REF" "$BAM" "$BED" "$OUT/cpp_th$th.vcf" "$SM" -th "$th" -z
  else
    echo "CPP: SKIP: vardictcpp unavailable ($VARDICTCPP_BIN --version failed)" >&2
  fi
done
for th in 2 4 8; do cmp <(grep -v '^##' "$OUT/th1.vcf") <(grep -v '^##' "$OUT/th$th.vcf"); done
for sz in 1000000 10000000 100000000 248956422; do
  LC_ALL=C /usr/bin/time -v -o "$OUT/R$sz.time" $TASKSET "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -th 8 -R "chr1:1-$sz" > "$OUT/R$sz.vcf"
  # -R form mirrors the rust flags verbatim plus -vcf (stdout VCF, never -o file).
  if [ "$CPP_OK" = "true" ]; then
    LC_ALL=C /usr/bin/time -v -o "$OUT/cpp_R$sz.time" $TASKSET $VARDICTCPP_BIN -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -th 8 -R "chr1:1-$sz" -vcf > "$OUT/cpp_R$sz.vcf"
  else
    echo "CPP: SKIP: vardictcpp unavailable ($VARDICTCPP_BIN --version failed)" >&2
  fi
done
python3 - "$OUT" "$PINNED" > "$OUT/extra.json" <<'PY'
import json, re, sys, pathlib
o = pathlib.Path(sys.argv[1]); pinned = sys.argv[2] == "true"
def get(f, k): return re.search(rf"{re.escape(k)}: (.+)", (o / f).read_text()).group(1)
def opt(f, k):
    p = o / f
    if not p.exists(): return None
    m = re.search(rf"{re.escape(k)}: (.+)", p.read_text())
    return m.group(1) if m else None
def secs(w):
    p = [float(x) for x in w.split(":")]
    return p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) == 3 else 0)
W = {t: secs(get(f"th{t}.time", "Elapsed (wall clock) time (h:mm:ss or m:ss)")) for t in (1, 2, 4, 8)}
M = {s: int(get(f"R{s}.time", "Maximum resident set size (kbytes)")) for s in (1000000, 10000000, 100000000, 248956422)}
def opt_secs(f):
    v = opt(f, "Elapsed (wall clock) time (h:mm:ss or m:ss)")
    return secs(v) if v is not None else None
def opt_rss(f):
    v = opt(f, "Maximum resident set size (kbytes)")
    return int(v) if v is not None else None
CW = {t: opt_secs(f"cpp_th{t}.time") for t in (1, 2, 4, 8)}
CM = {s: opt_rss(f"cpp_R{s}.time") for s in (1000000, 10000000, 100000000, 248956422)}
sp = W[1] / W[8]; rr = M[248956422] / M[10000000]
print(json.dumps({"wall_s": W, "maxrss_kb": M, "speedup_8": sp, "rss_ratio_chr1_vs_10mb": rr,
                  "pass": sp >= 4.8 and rr <= 1.5,
                  "cpp_wall_s": CW, "cpp_maxrss_kb": CM, "pinned": pinned,
                  "note": ">8 threads not measured (host policy A3)"}, indent=2))
PY
