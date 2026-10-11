#!/usr/bin/env bash
# Usage: run_slice.sh <id>  — one harness/slices.tsv row through compare.sh. Exit = compare.sh exit.
set -euo pipefail
ID=$1
row=$(awk -F'\t' -v id="$ID" '$1==id' harness/slices.tsv)
[ -n "$row" ] || { echo "run_slice: unknown id $ID" >&2; exit 2; }
IFS=$'\t' read -r _ REF BAM BED SAMPLE MINV EXTRA <<< "$row"
read -ra EX <<< "${EXTRA:-}"
KEEP_DIR="${SLICE_KEEP:-scratch/slices}/$ID" MIN_VARIANTS="$MINV" \
  bash harness/compare.sh "$REF" "$BAM" "$BED" "$SAMPLE" "${EX[@]}"
