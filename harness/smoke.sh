#!/usr/bin/env bash
# Release-artifact smoke gate: run ONLY the committed micro fixture through a
# built vardictrust binary and byte-compare the intermediate .var table against
# the dual-engine-validated golden/micro.expected.var.
#
# Usage: VARDICT_BIN=/path/to/vardictrust bash harness/smoke.sh
#   VARDICT_BIN defaults to target/release/vardictrust.
#
# Exit 0 + "SMOKE-PASS" on success; exit 1 + "SMOKE-FAIL" otherwise.
set -euo pipefail

BIN="${VARDICT_BIN:-target/release/vardictrust}"
FIXDIR="fixtures/micro"
EXPECTED="golden/micro.expected.var"

for f in "$BIN" "$FIXDIR/micro.ref.fa" "$FIXDIR/micro.ref.fa.fai" \
         "$FIXDIR/micro.bam" "$FIXDIR/micro.bam.bai" "$FIXDIR/micro.bed" "$EXPECTED"; do
    if [ ! -f "$f" ]; then
        echo "SMOKE-FAIL: missing $f (run from repo root?)"
        exit 1
    fi
done

TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT

"$BIN" -G "$FIXDIR/micro.ref.fa" -f 0.01 -N SAMPLE -b "$FIXDIR/micro.bam" \
    -c 1 -S 2 -E 3 -g 4 -z --emit-intermediate "$FIXDIR/micro.bed" > "$TMP"

if cmp -s "$EXPECTED" "$TMP"; then
    echo "SMOKE-PASS"
else
    echo "SMOKE-FAIL: output differs from $EXPECTED"
    diff "$EXPECTED" "$TMP" | head -20 || true
    exit 1
fi
