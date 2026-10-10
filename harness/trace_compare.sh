#!/usr/bin/env bash
set -euo pipefail
# Usage: trace_compare.sh <ref.fa> <bam|tumor|normal.bam> <bed> [sample] [extra VarDict args...]
# Runs both VarDictJava and VarDictRust with verbose (-y) and diffs their traces.
REF=$1; BAM=$2; BED=$3; SAMPLE=${4:-SAMPLE}; shift 4 || shift $#
AF=${AF_THR:-0.01}

# VarDictJava CLI doesn't support concatenated options like -f0.001; split them
EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --*) EXTRA_ARGS+=("$arg") ;;
    -[a-zA-Z]*)
      opt="${arg:0:2}"
      val="${arg:2}"
      if [ -n "$val" ]; then
        EXTRA_ARGS+=("$opt" "$val")
      else
        EXTRA_ARGS+=("$arg")
      fi
      ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

echo "=== Running VarDictJava with verbose ==="
conda run -n vardict-oracle bash -c \
  'vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 -y "${@:8}" "$5" > "$6" 2> "$7"' \
  _ "$REF" "$AF" "$SAMPLE" "$BAM" "$BED" "$TMP/java.var" "$TMP/java_trace.txt" "${EXTRA_ARGS[@]}"

echo "=== Running VarDictRust with verbose ==="
cargo run -q --release -p vardictrust -- -z -G "$REF" -f "$AF" -N "$SAMPLE" \
  -b "$BAM" -c 1 -S 2 -E 3 -g 4 --emit-intermediate -y "${EXTRA_ARGS[@]}" "$BED" \
  > "$TMP/rust.var" 2> "$TMP/rust_trace.txt"

echo ""
echo "=== Output comparison ==="
if diff -q "$TMP/java.var" "$TMP/rust.var" > /dev/null 2>&1; then
  echo "OUTPUT: MATCH"
else
  echo "OUTPUT: MISMATCH"
  echo "--- Java-only lines (first 10):"
  comm -23 <(sort "$TMP/java.var") <(sort "$TMP/rust.var") | head -10
  echo "--- Rust-only lines (first 10):"
  comm -13 <(sort "$TMP/java.var") <(sort "$TMP/rust.var") | head -10
fi

echo ""
echo "=== Trace comparison (key events) ==="
echo "--- Java realigndel events:"
grep -c "Realigndel for:" "$TMP/java_trace.txt" 2>/dev/null || echo "0"
echo "--- Rust realigndel events:"
grep -c "Realigndel for:" "$TMP/rust_trace.txt" 2>/dev/null || echo "0"
echo ""
echo "--- Java realignlgdel events:"
grep -c "Working Realignlgdel\|Found Realignlgdel\|Found lgdel" "$TMP/java_trace.txt" 2>/dev/null || echo "0"
echo "--- Rust realignlgdel events:"
grep -c "Working Realignlgdel\|Found Realignlgdel\|Found lgdel" "$TMP/rust_trace.txt" 2>/dev/null || echo "0"
