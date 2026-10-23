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
[[ "$groups" =~ "hg002" ]] || { echo "FAIL: missing hg002 group" >&2; exit 1; }


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
dashgrp	dash1	file://$TMPDIR/src/testfile.txt	-	test/dash.txt
badgrp	bad1	file://$TMPDIR/src/testfile.txt	00000000000000000000000000000000	test/bad.txt
EOF

# Test fetch matching group
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp

[ -f "$TMPDIR/data/test/file1.txt" ] || { echo "FAIL: file1.txt not fetched" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.done" ] || { echo "FAIL: file1.txt.done missing" >&2; exit 1; }

# Check MANIFEST
manifest_lines=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l)
[ "$manifest_lines" -eq 1 ] || { echo "FAIL: expected 1 line in manifest, got $manifest_lines" >&2; exit 1; }

# Idempotency test: run again, should print 'have:' and not add new manifest lines
out_second=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp)
echo "$out_second" | grep -q "have: test1" || { echo "FAIL: second run did not report have: test1" >&2; exit 1; }
manifest_lines_after=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l)
[ "$manifest_lines_after" -eq 1 ] || { echo "FAIL: second run modified manifest" >&2; exit 1; }

# Dash md5 test: should fail
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh dashgrp >/dev/null 2>&1
rc_dash=$?
set -e
[ "$rc_dash" -ne 0 ] || { echo "FAIL: dash md5 should have failed" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/dash.txt.done" ] || { echo "FAIL: dash md5 should not create .done file" >&2; exit 1; }

# Bad md5 test: should fail
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh badgrp >/dev/null 2>&1
rc_bad=$?
set -e
[ "$rc_bad" -ne 0 ] || { echo "FAIL: bad md5 should have failed" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/bad.txt.done" ] || { echo "FAIL: bad md5 should not create .done file" >&2; exit 1; }

# 4. Test stage_g4.sh usage & get_hg008_slice.sh --beds-only offline BED derivation
bash harness/stage_g4.sh --help | grep -q "stage_g4.sh" || { echo "FAIL: stage_g4.sh --help failed" >&2; exit 1; }
set +e
bash harness/stage_g4.sh unknown_target >/dev/null 2>&1
rc_unk=$?
set -e
[ "$rc_unk" -ne 0 ] || { echo "FAIL: stage_g4.sh unknown_target should fail" >&2; exit 1; }

mkdir -p "$TMPDIR/g4/hg008/truth"
cat << 'EOF' | gzip -c > "$TMPDIR/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10001	.	A	T	60	PASS	.
chr1	22001	.	C	G	60	PASS	.
chr1	35001	.	G	A	60	PASS	.
EOF
printf "chr1\t5000\t40000\n" > "$TMPDIR/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed"

DATA_DIR="$TMPDIR/g4" bash harness/get_hg008_slice.sh --beds-only >/dev/null
[ -s "$TMPDIR/g4/hg008/slices.bed" ] || { echo "FAIL: get_hg008_slice.sh did not create slices.bed" >&2; exit 1; }
[ -s "$TMPDIR/g4/hg008/engine.bed" ] || { echo "FAIL: get_hg008_slice.sh did not create engine.bed" >&2; exit 1; }
[ -s "$TMPDIR/g4/hg008/eval.bed" ]   || { echo "FAIL: get_hg008_slice.sh did not create eval.bed" >&2; exit 1; }
[ -s "$TMPDIR/g4/hg008/r1.bed" ]     || { echo "FAIL: get_hg008_slice.sh did not create r1.bed" >&2; exit 1; }
# Span is 9500..35501 (26001 bp > 20000 bp), so engine.bed must split into 2 tiles
[ "$(wc -l < "$TMPDIR/g4/hg008/engine.bed")" -eq 2 ] || { echo "FAIL: expected 2 tiles in engine.bed" >&2; exit 1; }

echo "PASS test_fetch_data"
