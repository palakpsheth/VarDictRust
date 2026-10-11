#!/usr/bin/env bash
# Usage: scale_bench.sh <ref> <bam> <bed> <sample> <outdir>
# Pre-registered (plan rev2): PASS iff wall(th1)/wall(th8) >= 4.8 AND maxRSS(chr1 full) <= 1.5 x maxRSS(10 Mb)
# AND outputs identical across thread counts. Thresholds must not change after seeing data.
set -euo pipefail
REF=$1; BAM=$2; BED=$3; SM=$4; OUT=$5; BIN=${VARDICT_BIN:-target/release/vardictrust}; mkdir -p "$OUT"
for th in 1 2 4 8; do
  /usr/bin/time -v -o "$OUT/th$th.time" "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th "$th" "$BED" > "$OUT/th$th.vcf"
done
for th in 2 4 8; do cmp <(grep -v '^##' "$OUT/th1.vcf") <(grep -v '^##' "$OUT/th$th.vcf"); done
for sz in 1000000 10000000 100000000 248956422; do
  /usr/bin/time -v -o "$OUT/R$sz.time" "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -th 8 -R "chr1:1-$sz" > "$OUT/R$sz.vcf"
done
python3 - "$OUT" > "$OUT/extra.json" <<'PY'
import json, re, sys, pathlib
o = pathlib.Path(sys.argv[1])
def get(f, k): return re.search(rf"{re.escape(k)}: (.+)", (o / f).read_text()).group(1)
def secs(w):
    p = [float(x) for x in w.split(":")]
    return p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) == 3 else 0)
W = {t: secs(get(f"th{t}.time", "Elapsed (wall clock) time (h:mm:ss or m:ss)")) for t in (1, 2, 4, 8)}
M = {s: int(get(f"R{s}.time", "Maximum resident set size (kbytes)")) for s in (1000000, 10000000, 100000000, 248956422)}
sp = W[1] / W[8]; rr = M[248956422] / M[10000000]
print(json.dumps({"wall_s": W, "maxrss_kb": M, "speedup_8": sp, "rss_ratio_chr1_vs_10mb": rr,
                  "pass": sp >= 4.8 and rr <= 1.5,
                  "note": ">8 threads not measured (host policy A3)"}, indent=2))
PY
