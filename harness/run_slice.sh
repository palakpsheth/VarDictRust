#!/usr/bin/env bash
# Usage: run_slice.sh <id>  — one harness/slices.tsv row through compare.sh. Exit = compare.sh exit.
# (compare.sh also runs the advisory vardictcpp third arm: stdout carries a
# trailing `CPP: <MATCH|MISMATCH|SKIP ...>` line and KEEP_DIR gains
# cpp.vcf/cpp.time/compare_cpp.txt. The gate — and this script's exit — stays
# java==rust only; CPP_SKIP=1 disables the arm via environment. No behavior
# change needed at this layer.)
set -euo pipefail
ID=$1
row=$(awk -F'\t' -v id="$ID" '$1==id' harness/slices.tsv)
[ -n "$row" ] || { echo "run_slice: unknown id $ID" >&2; exit 2; }
IFS=$'\t' read -r _ REF BAM BED SAMPLE MINV EXTRA <<< "$row"
read -ra EX <<< "${EXTRA:-}"
KEEP_DIR="${SLICE_KEEP:-scratch/slices}/$ID" MIN_VARIANTS="$MINV" \
  bash harness/compare.sh "$REF" "$BAM" "$BED" "$SAMPLE" "${EX[@]}"
