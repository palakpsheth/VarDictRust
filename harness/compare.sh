#!/usr/bin/env bash
set -euo pipefail
# Usage: compare.sh <ref> <bam|tumor|normal> <bed> [sample] [extra args...]
REF=$1; BAM=$2; BED=$3; SAMPLE=${4:-SAMPLE}; shift 4 || shift $#
TMP=$(mktemp -d)
# Mirror upstream IntegrationTest DEFAULT_ARGS ("-z -c 1 -S 2 -E 3 -g 4"): the
# corpus and VarDictJava's own suite always run zero-based. run_oracle.sh
# translates the flag-less -z into Java's explicit "-z 1" form.
bash harness/run_oracle.sh "$REF" "$BAM" "$BED" "$TMP/java.vcf" "$SAMPLE" -z "$@"
cargo run -q --release -p vardict-cli -- -G "$REF" -f "${AF_THR:-0.01}" -N "$SAMPLE" \
  -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z "$@" "$BED" > "$TMP/rust.vcf"
python3 harness/compare_vcf.py "$TMP/java.vcf" "$TMP/rust.vcf"; rc=$?
[ $rc -eq 0 ] && echo "MATCH" || { echo "MISMATCH"; diff "$TMP/java.vcf" "$TMP/rust.vcf" | head -50; }
rm -rf "$TMP"; exit $rc
