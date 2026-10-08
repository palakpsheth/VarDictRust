#!/usr/bin/env bash
set -euo pipefail
# Usage: compare.sh <ref> <bam|tumor|normal> <bed> [sample] [extra args...]
# Env: AF_THR (default 0.01) | KEEP_DIR (persist java.vcf rust.vcf compare.txt *.time)
#      MIN_VARIANTS (default 0) | CMP_TIMEOUT per engine (default 12h) | CMP_STRICT=1 (report-only, T11)
# Exit 0 = MATCH; non-zero = MISMATCH or engine failure. Decide on the EXIT CODE, never on text.
REF=$1; BAM=$2; BED=$3; SAMPLE=${4:-SAMPLE}; shift 4 || shift $#
if [ -n "${KEEP_DIR:-}" ]; then mkdir -p "$KEEP_DIR"; OUT=$KEEP_DIR
else OUT=$(mktemp -d); trap 'rm -rf "$OUT"' EXIT; fi
TO=${CMP_TIMEOUT:-12h}
# Mirror upstream IntegrationTest DEFAULT_ARGS ("-z -c 1 -S 2 -E 3 -g 4"): the
# corpus and VarDictJava's own suite always run zero-based. run_oracle.sh
# translates the flag-less -z into Java's explicit "-z 1" form.
/usr/bin/time -v -o "$OUT/java.time" timeout "$TO" \
  bash harness/run_oracle.sh "$REF" "$BAM" "$BED" "$OUT/java.vcf" "$SAMPLE" -z "$@" \
  || { echo "ENGINE-FAIL: oracle"; exit 2; }
/usr/bin/time -v -o "$OUT/rust.time" timeout "$TO" \
  cargo run -q --release -p vardict-cli -- -G "$REF" -f "${AF_THR:-0.01}" -N "$SAMPLE" \
  -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z "$@" "$BED" > "$OUT/rust.vcf" \
  || { echo "ENGINE-FAIL: vardictrust"; exit 2; }
rc=0
python3 harness/compare_vcf.py --min-variants "${MIN_VARIANTS:-0}" "$OUT/java.vcf" "$OUT/rust.vcf" \
  | tee "$OUT/compare.txt" || rc=$?
if [ "$rc" -eq 0 ]; then echo "MATCH"
else echo "MISMATCH"; diff "$OUT/java.vcf" "$OUT/rust.vcf" | head -50 || true; fi
exit "$rc"
