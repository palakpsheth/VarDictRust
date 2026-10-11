#!/usr/bin/env bash
set -euo pipefail
# Usage: compare.sh <ref> <bam|tumor|normal> <bed> [sample] [extra args...]
# Env: AF_THR (default 0.01) | KEEP_DIR (persist java.vcf rust.vcf compare.txt *.time +
#      cpp.vcf cpp.time compare_cpp.txt extra.json)
#      MIN_VARIANTS (default 0) | CMP_TIMEOUT per engine (default 12h) | CMP_STRICT=1 (report-only, T11)
#      CPP_SKIP=1 (skip advisory vardictcpp third arm; still prints `CPP: SKIP:`)
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
  cargo run -q --release -p vardictrust -- -G "$REF" -f "${AF_THR:-0.01}" -N "$SAMPLE" \
  -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z "$@" "$BED" > "$OUT/rust.vcf" \
  || { echo "ENGINE-FAIL: vardictrust"; exit 2; }
rc=0
STRICT_ARG=""
if [ "${CMP_STRICT:-0}" = "1" ] || [ "${CMP_STRICT:-0}" = "true" ]; then STRICT_ARG="--strict-bytes"; fi
python3 harness/compare_vcf.py $STRICT_ARG --min-variants "${MIN_VARIANTS:-0}" "$OUT/java.vcf" "$OUT/rust.vcf" \
  | tee "$OUT/compare.txt" || rc=$?
if [ "$rc" -eq 0 ]; then echo "MATCH"
else echo "MISMATCH"; diff "$OUT/java.vcf" "$OUT/rust.vcf" | head -50 || true; fi
# Third arm: vardictcpp (advisory only — stdout `CPP:` line + compare_cpp.txt;
# the gate verdict and exit code above are never touched by anything below).
# CPP_SKIP=1 escapes the arm; a missing conda/pin or an engine failure degrades
# to `CPP: SKIP: <reason>` with the java/rust verdict unchanged.
cpp_line=""
if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
  cpp_line="CPP: SKIP: disabled via CPP_SKIP=1"
elif ! command -v conda >/dev/null 2>&1; then
  cpp_line="CPP: SKIP: conda not on PATH"
elif LC_ALL=C /usr/bin/time -v -o "$OUT/cpp.time" timeout "$TO" \
    bash harness/run_cpp.sh "$REF" "$BAM" "$BED" "$OUT/cpp.vcf" "$SAMPLE" -z "$@" 2>"$OUT/cpp.err"; then
  # Empty-BED/variant-free runs leave a 0-byte file (run_cpp.sh), not a VCF:
  # fall back to a header-only VCF so compare_vcf.py sees a valid empty VCF.
  if [ ! -s "$OUT/cpp.vcf" ]; then
    grep '^#' "$OUT/java.vcf" > "$OUT/cpp.vcf" || printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$OUT/cpp.vcf"
  fi
  rc1=0; rc2=0
  python3 harness/compare_vcf.py $STRICT_ARG --min-variants "${MIN_VARIANTS:-0}" "$OUT/java.vcf" "$OUT/cpp.vcf" \
    > "$OUT/compare_cpp.txt" 2>&1 || rc1=$?
  echo "--- cpp vs rust ---" >> "$OUT/compare_cpp.txt"
  python3 harness/compare_vcf.py $STRICT_ARG --min-variants "${MIN_VARIANTS:-0}" "$OUT/cpp.vcf" "$OUT/rust.vcf" \
    >> "$OUT/compare_cpp.txt" 2>&1 || rc2=$?
  if [ "$rc1" -eq 0 ] && [ "$rc2" -eq 0 ]; then cpp_line="CPP: MATCH"; else cpp_line="CPP: MISMATCH"; fi
  echo "$cpp_line" >> "$OUT/compare_cpp.txt"
else
  cpp_exit=$?
  if [ "$cpp_exit" -eq 3 ]; then
    reason=$(head -1 "$OUT/cpp.err" 2>/dev/null || true)
    reason=${reason#SKIP: }
    cpp_line="CPP: SKIP: ${reason:-vardictcpp unavailable (run_cpp.sh exit 3)}"
  else
    echo "ENGINE-FAIL: cpp"
    cpp_line="CPP: SKIP: engine failure (exit $cpp_exit)"
  fi
  echo "$cpp_line" > "$OUT/compare_cpp.txt"
fi
if [ -n "$cpp_line" ]; then
  [ -f "$OUT/compare_cpp.txt" ] || echo "$cpp_line" > "$OUT/compare_cpp.txt"
  echo "$cpp_line"
fi
# Persist third-arm evidence keys for record_evidence.py (ev_cpp presence gate):
# success verdicts -> extra.cpp_version; any SKIP -> extra.cpp_skip_reason.
# Merge-only (never clobbers other extra.* keys) and never touches the gate exit.
if [ -n "${KEEP_DIR:-}" ]; then
  case "$cpp_line" in
    "CPP: SKIP: "*)
      cpp_key="cpp_skip_reason"; cpp_val="${cpp_line#CPP: SKIP: }" ;;
    *)
      cpp_key="cpp_version"
      cpp_val=$(conda list -n vardict-oracle vardictcpp 2>/dev/null | awk '/^vardictcpp[[:space:]]/{print "vardictcpp " $2 " " $3 " bioconda"}')
      [ -n "$cpp_val" ] || cpp_val=$(conda run -n vardict-oracle vardictcpp --version 2>/dev/null | head -n 1) ;;
  esac
  if [ -n "${cpp_val:-}" ]; then
    python3 - "$OUT/extra.json" "$cpp_key" "$cpp_val" <<'PYEOF' || true
import json, sys
p, k, v = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    with open(p) as f: d = json.load(f)
    if not isinstance(d, dict): d = {}
except Exception: d = {}
d[k] = v
with open(p, "w") as f: json.dump(d, f, indent=2, sort_keys=True)
PYEOF
  fi
fi
exit "$rc"
