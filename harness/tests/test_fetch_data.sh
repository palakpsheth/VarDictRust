#!/usr/bin/env bash
set -euo pipefail

# 1. Validate harness/assets.tsv format
ASSETS="harness/assets.tsv"
[ -f "$ASSETS" ] || { echo "FAIL: $ASSETS not found" >&2; exit 1; }

header=$(head -n 1 "$ASSETS")
[ "$header" = "# group	name	url	md5	dest" ] || {
  echo "FAIL: Unexpected header in $ASSETS: $header" >&2; exit 1;
}

line_num=1
while IFS=$'\t' read -r group name url md5 dest extra; do
  line_num=$((line_num + 1))
  case "$group" in ''|\#*) continue ;; esac
  [ -n "$group" ] || { echo "FAIL line $line_num: missing group" >&2; exit 1; }
  [ -n "$name" ] || { echo "FAIL line $line_num: missing name" >&2; exit 1; }
  [ -n "$url" ] || { echo "FAIL line $line_num: missing url" >&2; exit 1; }
  [ -n "$md5" ] || { echo "FAIL line $line_num: missing md5" >&2; exit 1; }
  [ -n "$dest" ] || { echo "FAIL line $line_num: missing dest" >&2; exit 1; }
  [ -z "$extra" ] || { echo "FAIL line $line_num: too many columns" >&2; exit 1; }
  [[ "$url" =~ ^https?:// ]] || { echo "FAIL line $line_num: invalid url scheme $url" >&2; exit 1; }
done < "$ASSETS"

# Check required groups present
groups=$(awk '!/^#/ && NF {print $1}' "$ASSETS" | sort -u | tr '\n' ' ')
[[ "$groups" =~ "refs" ]] || { echo "FAIL: missing refs group" >&2; exit 1; }
[[ "$groups" =~ "seqc2" ]] || { echo "FAIL: missing seqc2 group" >&2; exit 1; }
[[ "$groups" =~ "hg008" ]] || { echo "FAIL: missing hg008 group" >&2; exit 1; }

# 2. Test harness/fetch_data.sh usage error on no arguments
set +e
bash harness/fetch_data.sh >/dev/null 2>&1
rc=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL: fetch_data.sh with no args should fail" >&2; exit 1; }

# 3. Test functional behavior with a temporary mock asset
TMPDIR=$(mktemp -d -t test_fetch_XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p "$TMPDIR/src" "$TMPDIR/data" "$TMPDIR/benchmarks"
echo "hello world" > "$TMPDIR/src/testfile.txt"
test_md5=$(md5sum "$TMPDIR/src/testfile.txt" | cut -d' ' -f1)

cat << EOF > "$TMPDIR/test_assets.tsv"
# group	name	url	md5	dest
testgrp	test1	file://$TMPDIR/src/testfile.txt	$test_md5	test/file1.txt
testgrp	test2	file://$TMPDIR/src/testfile.txt	-	test/file2.txt
badgrp	bad1	file://$TMPDIR/src/testfile.txt	00000000000000000000000000000000	test/bad.txt
EOF

# Test fetch matching group
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp

[ -f "$TMPDIR/data/test/file1.txt" ] || { echo "FAIL: file1.txt not fetched" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.done" ] || { echo "FAIL: file1.txt.done missing" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file2.txt" ] || { echo "FAIL: file2.txt not fetched" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file2.txt.done" ] || { echo "FAIL: file2.txt.done missing" >&2; exit 1; }

# Check MANIFEST
manifest_lines=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l)
[ "$manifest_lines" -eq 2 ] || { echo "FAIL: expected 2 lines in manifest, got $manifest_lines" >&2; exit 1; }

# Idempotency test: run again, should print 'have:' and not add new manifest lines
out_second=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp)
echo "$out_second" | grep -q "have: test1" || { echo "FAIL: second run did not report have: test1" >&2; exit 1; }
manifest_lines_after=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l)
[ "$manifest_lines_after" -eq 2 ] || { echo "FAIL: second run modified manifest" >&2; exit 1; }

# Bad md5 test: should fail
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh badgrp >/dev/null 2>&1
rc_bad=$?
set -e
[ "$rc_bad" -ne 0 ] || { echo "FAIL: bad md5 should have failed" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/bad.txt.done" ] || { echo "FAIL: bad md5 should not create .done file" >&2; exit 1; }

echo "PASS test_fetch_data"
