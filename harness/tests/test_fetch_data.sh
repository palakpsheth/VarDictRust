#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/lib/portable.sh"

# 1. Validate harness/assets.tsv format
ASSETS="harness/assets.tsv"
[ -f "$ASSETS" ] || { echo "FAIL: $ASSETS not found" >&2; exit 1; }

header=$(head -n 1 "$ASSETS")
[ "$header" = "# group	name	url	md5	sha256	dest" ] || {
  echo "FAIL: Unexpected header in $ASSETS: $header" >&2; exit 1;
}

line_num=1
while IFS=$'\t' read -r group name url md5 sha256 dest extra; do
  line_num=$((line_num + 1))
  case "$group" in ''|\#*) continue ;; esac
  [ -n "$group" ] || { echo "FAIL line $line_num: missing group" >&2; exit 1; }
  [ -n "$name" ] || { echo "FAIL line $line_num: missing name" >&2; exit 1; }
  [ -n "$url" ] || { echo "FAIL line $line_num: missing url" >&2; exit 1; }
  [ -n "$md5" ] || { echo "FAIL line $line_num: missing md5" >&2; exit 1; }
  [ -n "$sha256" ] || { echo "FAIL line $line_num: missing sha256" >&2; exit 1; }
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
test_md5=$(portable_md5 "$TMPDIR/src/testfile.txt")
test_sha256=$(portable_sha256 "$TMPDIR/src/testfile.txt")

cat << EOF > "$TMPDIR/test_assets.tsv"
# group	name	url	md5	sha256	dest
testgrp	test1	file://$TMPDIR/src/testfile.txt	$test_md5	$test_sha256	test/file1.txt
dashgrp	dash1	file://$TMPDIR/src/testfile.txt	$test_md5	-	test/dash.txt
tampergrp	tamper1	file://$TMPDIR/src/testfile.txt	$test_md5	0000000000000000000000000000000000000000000000000000000000000000	test/tamper.txt
badgrp	bad1	file://$TMPDIR/src/testfile.txt	00000000000000000000000000000000	-	test/bad.txt
prefixgrp	foo.vcf.gz	file://$TMPDIR/src/testfile.txt	$test_md5	$test_sha256	test/foo.vcf.gz
prefixgrp	foo.vcf.gz.tbi	file://$TMPDIR/src/testfile.txt	$test_md5	$test_sha256	test/foo.vcf.gz.tbi
EOF

# Test fetch matching group (happy-path SHA-256 + MANIFEST header & exact col-1 append)
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp

[ -f "$TMPDIR/data/test/file1.txt" ] || { echo "FAIL: file1.txt not fetched" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.done" ] || { echo "FAIL: file1.txt.done missing" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.sha256" ] || { echo "FAIL: file1.txt.sha256 marker missing" >&2; exit 1; }
[ "$(cat "$TMPDIR/data/test/file1.txt.sha256")" = "$test_sha256" ] || { echo "FAIL: file1.txt.sha256 content mismatch" >&2; exit 1; }

# Check MANIFEST format and entry
man_hdr=$(head -n 1 "$TMPDIR/benchmarks/MANIFEST")
[ "$man_hdr" = "# name	url	sha256	checksum_source	bytes	fetched" ] || {
  echo "FAIL: Unexpected MANIFEST header: $man_hdr" >&2; exit 1;
}
manifest_lines=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l)
[ "$manifest_lines" -eq 1 ] || { echo "FAIL: expected 1 line in manifest, got $manifest_lines" >&2; exit 1; }
man_entry=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST")
[[ "$man_entry" =~ ^test1$'\t'file://.*$'\t'"$test_sha256"$'\t'upstream$'\t'12$'\t' ]] || {
  echo "FAIL: MANIFEST entry did not record exact SHA-256 and metadata: $man_entry" >&2; exit 1;
}

# Idempotency test: run again, should print 'have:' and not add new manifest lines
out_second=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp)
echo "$out_second" | grep -q "have: test1" || { echo "FAIL: second run did not report have: test1" >&2; exit 1; }
manifest_lines_after=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l)
[ "$manifest_lines_after" -eq 1 ] || { echo "FAIL: second run modified manifest" >&2; exit 1; }

# Pre-migration promotion test: file exists with valid SHA-256 and legacy .done file, missing .sha256 marker
rm -f "$TMPDIR/data/test/file1.txt.sha256"
out_promote=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp)
[ -f "$TMPDIR/data/test/file1.txt.done" ] || { echo "FAIL: in-place promotion did not keep/create .done file" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.sha256" ] || { echo "FAIL: in-place promotion did not create .sha256 marker" >&2; exit 1; }
[ "$(cat "$TMPDIR/data/test/file1.txt.sha256")" = "$test_sha256" ] || { echo "FAIL: promoted .sha256 marker mismatch" >&2; exit 1; }
echo "$out_promote" | grep -Eq "(have|verified): test1" || { echo "FAIL: in-place promotion did not recognize existing file: $out_promote" >&2; exit 1; }

# Tampered payload with stale .done marker: must detect mismatch, re-fetch, and overwrite
rm -f "$TMPDIR/data/test/file1.txt.sha256"
echo "tampered content" > "$TMPDIR/data/test/file1.txt"
touch "$TMPDIR/data/test/file1.txt.done"
out_tamper_recover=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp)
[ "$(cat "$TMPDIR/data/test/file1.txt")" = "hello world" ] || { echo "FAIL: tampered file was not re-fetched and restored" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.sha256" ] || { echo "FAIL: re-fetch did not recreate .sha256 marker" >&2; exit 1; }
[ "$(cat "$TMPDIR/data/test/file1.txt.sha256")" = "$test_sha256" ] || { echo "FAIL: re-fetch .sha256 content mismatch" >&2; exit 1; }

# Tampered SHA-256 test: should fail non-zero and clean .part file
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh tampergrp >/dev/null 2>&1
rc_tamper=$?
set -e
[ "$rc_tamper" -ne 0 ] || { echo "FAIL: tampered SHA-256 should have failed" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/tamper.txt" ] || { echo "FAIL: tampered payload should not persist at dest" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/tamper.txt.part" ] || { echo "FAIL: tampered payload should clean .part file" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/tamper.txt.done" ] || { echo "FAIL: tampered payload should not create .done file" >&2; exit 1; }

# Missing SHA-256 (sha256 = "-"): should fail closed by default
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh dashgrp >/dev/null 2>&1
rc_dash=$?
set -e
[ "$rc_dash" -ne 0 ] || { echo "FAIL: missing sha256 should fail closed by default" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/dash.txt.done" ] || { echo "FAIL: missing sha256 should not create .done file" >&2; exit 1; }

# Opt-in --allow-legacy-md5: should succeed when MD5 matches
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh --allow-legacy-md5 dashgrp
[ -f "$TMPDIR/data/test/dash.txt" ] || { echo "FAIL: --allow-legacy-md5 did not fetch dash.txt" >&2; exit 1; }
[ -f "$TMPDIR/data/test/dash.txt.done" ] || { echo "FAIL: --allow-legacy-md5 did not create dash.txt.done" >&2; exit 1; }
[ -f "$TMPDIR/data/test/dash.txt.legacy_md5" ] || { echo "FAIL: --allow-legacy-md5 did not create dash.txt.legacy_md5" >&2; exit 1; }
[ "$(cat "$TMPDIR/data/test/dash.txt.legacy_md5")" = "$test_md5" ] || { echo "FAIL: dash.txt.legacy_md5 content mismatch" >&2; exit 1; }

# Subsequent run on legacy asset without --allow-legacy-md5 must fail closed even when dash.txt and markers exist
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh dashgrp >/dev/null 2>&1
rc_dash_subsequent=$?
set -e
[ "$rc_dash_subsequent" -ne 0 ] || { echo "FAIL: subsequent run without --allow-legacy-md5 should fail closed" >&2; exit 1; }

# Opt-in --allow-legacy-md5 with bad MD5: should fail
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh --allow-legacy-md5 badgrp >/dev/null 2>&1
rc_bad=$?
set -e
[ "$rc_bad" -ne 0 ] || { echo "FAIL: bad MD5 with --allow-legacy-md5 should fail" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/bad.txt.done" ] || { echo "FAIL: bad MD5 should not create .done file" >&2; exit 1; }

# Prefix collision test: foo.vcf.gz vs foo.vcf.gz.tbi
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh prefixgrp
[ -f "$TMPDIR/data/test/foo.vcf.gz.done" ] || { echo "FAIL: foo.vcf.gz.done missing" >&2; exit 1; }
[ -f "$TMPDIR/data/test/foo.vcf.gz.tbi.done" ] || { echo "FAIL: foo.vcf.gz.tbi.done missing" >&2; exit 1; }
foo_count=$(awk -F'\t' '$1 == "foo.vcf.gz" {count++} END{print count+0}' "$TMPDIR/benchmarks/MANIFEST")
footbi_count=$(awk -F'\t' '$1 == "foo.vcf.gz.tbi" {count++} END{print count+0}' "$TMPDIR/benchmarks/MANIFEST")
[ "$foo_count" -eq 1 ] || { echo "FAIL: expected exactly 1 foo.vcf.gz in MANIFEST, got $foo_count" >&2; exit 1; }
[ "$footbi_count" -eq 1 ] || { echo "FAIL: expected exactly 1 foo.vcf.gz.tbi in MANIFEST, got $footbi_count" >&2; exit 1; }

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

# Verify get_hg008_slice.sh preserves r1.bed mtime when slices.bed is unchanged
r1_mtime_before=$(python3 -c "import os, sys; print(os.stat(sys.argv[1]).st_mtime_ns)" "$TMPDIR/g4/hg008/r1.bed")
sleep 0.1
DATA_DIR="$TMPDIR/g4" bash harness/get_hg008_slice.sh --beds-only >/dev/null
r1_mtime_after=$(python3 -c "import os, sys; print(os.stat(sys.argv[1]).st_mtime_ns)" "$TMPDIR/g4/hg008/r1.bed")
[ "$r1_mtime_before" = "$r1_mtime_after" ] || { echo "FAIL: get_hg008_slice.sh should preserve r1.bed mtime when unchanged" >&2; exit 1; }

# Verify stage_g4.sh indexes GRCh38_GIABv3.fasta before linking .fa.fai and passes registered SRR7890879/SRR7890880 FASTQ paths
grep -q 'samtools faidx "\$gv/GRCh38_GIABv3.fasta"' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing samtools faidx on GRCh38_GIABv3.fasta" >&2; exit 1; }
grep -q 'seqc2_align.sh HCC1395 "\$sq/fastq/SRR7890879_1.fastq.gz" "\$sq/fastq/SRR7890879_2.fastq.gz"' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing SRR7890879 FASTQ paths" >&2; exit 1; }
grep -q 'seqc2_align.sh HCC1395BL "\$sq/fastq/SRR7890880_1.fastq.gz" "\$sq/fastq/SRR7890880_2.fastq.gz"' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing SRR7890880 FASTQ paths" >&2; exit 1; }
grep -q 'DATA_DIR="\$D" bash harness/get_hg002_slice.sh' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing DATA_DIR pass to get_hg002_slice.sh" >&2; exit 1; }
grep -q 'DATA_DIR="\$D" bash harness/get_hg008_slice.sh' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing DATA_DIR pass to get_hg008_slice.sh" >&2; exit 1; }

echo "PASS test_fetch_data"
