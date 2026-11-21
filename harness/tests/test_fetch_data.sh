#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/lib/portable.sh"

# 1. Validate 6-column harness/assets.tsv format (# group\tname\turl\tmd5\tsha256\tdest)
ASSETS="harness/assets.tsv"
[ -f "$ASSETS" ] || { echo "FAIL: $ASSETS not found" >&2; exit 1; }

header=$(head -n 1 "$ASSETS")
[ "$header" = "# group	name	url	md5	sha256	dest" ] || {
  echo "FAIL: Unexpected header in $ASSETS: $header" >&2; exit 1;
}

line_num=1
found_hg002_bai=0
while IFS=$'\t' read -r group name url md5 sha256 dest extra; do
  line_num=$((line_num + 1))
  case "$group" in ''|\#*) continue ;; esac
  [ -n "$group" ]  || { echo "FAIL line $line_num: missing group" >&2; exit 1; }
  [ -n "$name" ]   || { echo "FAIL line $line_num: missing name" >&2; exit 1; }
  [ -n "$url" ]    || { echo "FAIL line $line_num: missing url" >&2; exit 1; }
  [ -n "$md5" ]    || { echo "FAIL line $line_num: missing md5" >&2; exit 1; }
  [ -n "$sha256" ] || { echo "FAIL line $line_num: missing sha256" >&2; exit 1; }
  [ -n "$dest" ]   || { echo "FAIL line $line_num: missing dest" >&2; exit 1; }
  [ -z "$extra" ]  || { echo "FAIL line $line_num: too many columns" >&2; exit 1; }
  [[ "$url" =~ ^https?:// ]] || { echo "FAIL line $line_num: invalid url scheme $url" >&2; exit 1; }
  [[ "$md5" =~ ^([0-9a-f]{32}|-)$ ]] || { echo "FAIL line $line_num: invalid md5 '$md5'" >&2; exit 1; }
  [[ "$sha256" =~ ^([0-9a-f]{64}|-)$ ]] || { echo "FAIL line $line_num: invalid sha256 '$sha256'" >&2; exit 1; }
  case "$group" in
    refs|seqc2|hg008|hg002|stage6)
      [[ "$sha256" =~ ^[0-9a-f]{64}$ ]] || {
        echo "FAIL line $line_num: core group '$group' asset '$name' must have pinned 64-hex sha256" >&2
        exit 1
      }
      ;;
  esac
  if [ "$name" = "HG002.remote.bai" ]; then
    found_hg002_bai=1
    [ "$md5" = "a3c2b449df6509ca83fbd3fea22b9aee" ] || {
      echo "FAIL: HG002.remote.bai md5 not reconciled ($md5)" >&2; exit 1;
    }
    [ "$sha256" = "723f69de14b8adbbf1a60f23f9af007310011b51b13ee01d18a2a8daf7282f36" ] || {
      echo "FAIL: HG002.remote.bai sha256 mismatch ($sha256)" >&2; exit 1;
    }
  fi
done < "$ASSETS"
[ "$found_hg002_bai" -eq 1 ] || { echo "FAIL: HG002.remote.bai missing from $ASSETS" >&2; exit 1; }

# Check required groups present
groups=$(awk '!/^#/ && NF {print $1}' "$ASSETS" | sort -u | tr '\n' ' ')
[[ "$groups" =~ "refs" ]] || { echo "FAIL: missing refs group" >&2; exit 1; }
[[ "$groups" =~ "seqc2" ]] || { echo "FAIL: missing seqc2 group" >&2; exit 1; }
[[ "$groups" =~ "hg008" ]] || { echo "FAIL: missing hg008 group" >&2; exit 1; }
[[ "$groups" =~ "hg002" ]] || { echo "FAIL: missing hg002 group" >&2; exit 1; }
[[ "$groups" =~ "stage6" ]] || { echo "FAIL: missing stage6 group" >&2; exit 1; }

# 1b. Every stage6_fetch_verified target in harness/run_stage_6.sh must resolve to a 64-hex SHA-256 in harness/assets.tsv
python3 - "$ASSETS" "harness/run_stage_6.sh" <<'PY'
import pathlib, re, sys
assets_path, stage6_path = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
sha_by_name = {}
for line in assets_path.read_text(encoding="utf-8").splitlines():
    if not line or line.startswith("#"):
        continue
    cols = line.split("\t")
    if len(cols) == 6 and re.fullmatch(r"[0-9a-f]{64}", cols[4]):
        sha_by_name[cols[1]] = cols[4]
s6_text = stage6_path.read_text(encoding="utf-8")
direct_names = re.findall(r'stage6_fetch_verified\s+"([^"$]+)"', s6_text)
strat_names = [
    pathlib.PurePosixPath(m).name
    for m in re.findall(r'"[A-Za-z0-9_]+\|([A-Za-z0-9_./-]+\.bed\.gz)"', s6_text)
]
required = sorted(set(direct_names + strat_names))
assert len(required) >= 27, f"expected >=27 stage6_fetch_verified assets in run_stage_6.sh, parsed {len(required)}"
missing = [n for n in required if n not in sha_by_name]
if missing:
    sys.exit(f"FAIL: stage6_fetch_verified asset(s) missing 64-hex SHA-256 in {assets_path}: {missing}")
PY

# 2. Test harness/fetch_data.sh usage error on no arguments or flag-only invocation
set +e
bash harness/fetch_data.sh >/dev/null 2>&1
rc=$?
bash harness/fetch_data.sh --allow-legacy-md5 >/dev/null 2>&1
rc_flag_only=$?
set -e
[ "$rc" -ne 0 ] || { echo "FAIL: fetch_data.sh with no args should fail" >&2; exit 1; }
[ "$rc_flag_only" -ne 0 ] || { echo "FAIL: fetch_data.sh with only --allow-legacy-md5 should fail" >&2; exit 1; }

# 3. Test stage_verified_sha256 and fetch_data.sh with temporary mock assets
TMPDIR=$(mktemp -d -t test_fetch_XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p "$TMPDIR/src" "$TMPDIR/data" "$TMPDIR/benchmarks" "$TMPDIR/stage"
echo "hello world" > "$TMPDIR/src/testfile.txt"
echo "updated payload" > "$TMPDIR/src/updated.txt"
echo "tbi-bytes" > "$TMPDIR/src/foo.vcf.gz.tbi"
echo "vcf-bytes" > "$TMPDIR/src/foo.vcf.gz"

test_md5=$(portable_md5 "$TMPDIR/src/testfile.txt")
test_sha=$(portable_sha256 "$TMPDIR/src/testfile.txt")
upd_md5=$(portable_md5 "$TMPDIR/src/updated.txt")
upd_sha=$(portable_sha256 "$TMPDIR/src/updated.txt")
tbi_md5=$(portable_md5 "$TMPDIR/src/foo.vcf.gz.tbi")
tbi_sha=$(portable_sha256 "$TMPDIR/src/foo.vcf.gz.tbi")
vcf_md5=$(portable_md5 "$TMPDIR/src/foo.vcf.gz")
vcf_sha=$(portable_sha256 "$TMPDIR/src/foo.vcf.gz")

# 3a. Direct unit test of stage_verified_sha256
stage_verified_sha256 "file://$TMPDIR/src/testfile.txt" "$test_sha" "$TMPDIR/stage/good.txt"
[ -f "$TMPDIR/stage/good.txt" ] || { echo "FAIL: stage_verified_sha256 did not write good.txt" >&2; exit 1; }
[ ! -f "$TMPDIR/stage/good.txt.part" ] || { echo "FAIL: leftover .part file after stage_verified_sha256" >&2; exit 1; }

# Reuse in-place without network when existing file matches SHA-256
stage_verified_sha256 "file://$TMPDIR/src/does_not_exist.txt" "$test_sha" "$TMPDIR/stage/good.txt"

# Tampered SHA-256 rejected and cleans up .part
set +e
stage_verified_sha256 "file://$TMPDIR/src/testfile.txt" "0000000000000000000000000000000000000000000000000000000000000000" "$TMPDIR/stage/tampered.txt" >/dev/null 2>&1
rc_stage_bad=$?
set -e
[ "$rc_stage_bad" -ne 0 ] || { echo "FAIL: stage_verified_sha256 accepted wrong SHA-256" >&2; exit 1; }
[ ! -f "$TMPDIR/stage/tampered.txt" ] || { echo "FAIL: stage_verified_sha256 created destination on mismatch" >&2; exit 1; }
[ ! -f "$TMPDIR/stage/tampered.txt.part" ] || { echo "FAIL: stage_verified_sha256 left .part on mismatch" >&2; exit 1; }

cat << EOF > "$TMPDIR/test_assets.tsv"
# group	name	url	md5	sha256	dest
testgrp	test1	file://$TMPDIR/src/testfile.txt	$test_md5	$test_sha	test/file1.txt
badsha	bad1	file://$TMPDIR/src/testfile.txt	$test_md5	0000000000000000000000000000000000000000000000000000000000000000	test/bad_sha.txt
legacymd5	leg1	file://$TMPDIR/src/testfile.txt	$test_md5	-	test/legacy.txt
legacybad	legbad1	file://$TMPDIR/src/testfile.txt	00000000000000000000000000000000	-	test/legacy_bad.txt
dashboth	dash1	file://$TMPDIR/src/testfile.txt	-	-	test/dash.txt
promotegrp	promote1	file://$TMPDIR/src/never_download_me.txt	$test_md5	$test_sha	test/preexisting.txt
prefixgrp	foo.vcf.gz.tbi	file://$TMPDIR/src/foo.vcf.gz.tbi	$tbi_md5	$tbi_sha	test/foo.vcf.gz.tbi
prefixgrp	foo.vcf.gz	file://$TMPDIR/src/foo.vcf.gz	$vcf_md5	$vcf_sha	test/foo.vcf.gz
EOF

# 3b. Happy-path 6-column SHA-256 fetch + MANIFEST (name, url, digest, checksum_source) dedup
mkdir -p "$TMPDIR/benchmarks"
printf '# name\turl\tsha256\tchecksum_source\tbytes\tfetched\n' > "$TMPDIR/benchmarks/MANIFEST"
printf 'test1\tfile://%s/src/testfile.txt\t%s\tcomputed-locally\t12\t2026-10-08T00:00:00Z\n' "$TMPDIR" "$test_md5" >> "$TMPDIR/benchmarks/MANIFEST"

DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp

[ -f "$TMPDIR/data/test/file1.txt" ] || { echo "FAIL: file1.txt not fetched" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.done" ] || { echo "FAIL: file1.txt.done missing" >&2; exit 1; }
[ -f "$TMPDIR/data/test/file1.txt.sha256" ] || { echo "FAIL: file1.txt.sha256 missing" >&2; exit 1; }
[ "$(cat "$TMPDIR/data/test/file1.txt.sha256")" = "$test_sha" ] || { echo "FAIL: file1.txt.sha256 digest mismatch" >&2; exit 1; }

# Historical MD5 row must stay intact AND new SHA-256 receipt must be appended (2 rows for test1)
manifest_lines=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l | tr -d '[:space:]')
[ "$manifest_lines" -eq 2 ] || { echo "FAIL: expected historical MD5 + new SHA-256 rows (2), got $manifest_lines" >&2; exit 1; }
awk -F'\t' -v n="test1" -v d="$test_sha" -v s="upstream" '!/^#/ && $1 == n && $3 == d && $4 == s {f=1} END{exit !f}' "$TMPDIR/benchmarks/MANIFEST" || {
  echo "FAIL: new SHA-256 receipt was suppressed by historical MD5 row" >&2; exit 1;
}

#Also verify freshly initialized MANIFEST header when file does not exist yet
rm -f "$TMPDIR/fresh_MANIFEST"
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/fresh_MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp >/dev/null
man_hdr=$(head -n 1 "$TMPDIR/fresh_MANIFEST")
[ "$man_hdr" = "# name	url	sha256	checksum_source	bytes	fetched" ] || {
  echo "FAIL: unexpected MANIFEST header: $man_hdr" >&2; exit 1;
}

# Idempotency test: second run prints 'have:' and does not duplicate manifest line
out_second=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp)
grep -q "have: test1" <<< "$out_second" || { echo "FAIL: second run did not report have: test1" >&2; exit 1; }
manifest_lines_after=$(grep -v '^#' "$TMPDIR/benchmarks/MANIFEST" | wc -l | tr -d '[:space:]')
[ "$manifest_lines_after" -eq 2 ] || { echo "FAIL: second run modified manifest" >&2; exit 1; }

# Corrupting/truncating $out while .sha256 + .done remain intact MUST rehash, detect tampering, and re-fetch!
printf 'TAMPERED_AFTER_SHA256\n' > "$TMPDIR/data/test/file1.txt"
out_repair_sha=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh testgrp)
grep -q "fetch: test1" <<< "$out_repair_sha" || {
  echo "FAIL: modified cached SHA-256 payload was blindly accepted via .sha256 sidecar" >&2; exit 1;
}
[ "$(portable_sha256 "$TMPDIR/data/test/file1.txt")" = "$test_sha" ] || {
  echo "FAIL: modified cached SHA-256 payload was not repaired" >&2; exit 1;
}

# 3c. Tampered SHA-256 payload rejection
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh badsha >/dev/null 2>&1
rc_badsha=$?
set -e
[ "$rc_badsha" -ne 0 ] || { echo "FAIL: bad sha256 should have failed" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/bad_sha.txt" ] || { echo "FAIL: bad sha256 left bad_sha.txt" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/bad_sha.txt.part" ] || { echo "FAIL: bad sha256 left bad_sha.txt.part" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/bad_sha.txt.done" ] || { echo "FAIL: bad sha256 created .done" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/bad_sha.txt.sha256" ] || { echo "FAIL: bad sha256 created .sha256" >&2; exit 1; }

# 3d. Missing sha256='-' fails closed by default, succeeds ONLY with --allow-legacy-md5 + matching md5
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh legacymd5 >/dev/null 2>&1
rc_leg_default=$?
set -e
[ "$rc_leg_default" -ne 0 ] || { echo "FAIL: missing sha256='-' must fail closed without --allow-legacy-md5" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/legacy.txt.done" ] || { echo "FAIL: default run created legacy.txt.done" >&2; exit 1; }

DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh --allow-legacy-md5 legacymd5
[ -f "$TMPDIR/data/test/legacy.txt" ] || { echo "FAIL: --allow-legacy-md5 did not fetch legacy.txt" >&2; exit 1; }
[ -f "$TMPDIR/data/test/legacy.txt.legacy_md5" ] || { echo "FAIL: legacy.txt.legacy_md5 marker missing" >&2; exit 1; }
[ -f "$TMPDIR/data/test/legacy.txt.done" ] || { echo "FAIL: legacy.txt.done marker missing" >&2; exit 1; }

# Modifying/truncating $out while .legacy_md5 + .done remain intact MUST rehash MD5, detect tampering, and re-fetch!
printf 'TAMPERED_AFTER_LEGACY_MD5\n' > "$TMPDIR/data/test/legacy.txt"
out_repair_md5=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh --allow-legacy-md5 legacymd5)
grep -q "fetch: leg1" <<< "$out_repair_md5" || {
  echo "FAIL: modified legacy-MD5 payload was blindly accepted via .legacy_md5 sidecar" >&2; exit 1;
}
[ "$(portable_md5 "$TMPDIR/data/test/legacy.txt")" = "$test_md5" ] || {
  echo "FAIL: modified legacy-MD5 payload was not repaired" >&2; exit 1;
}

# Fail-closed check must still fire BEFORE inspecting .done / .legacy_md5 markers when --allow-legacy-md5 is omitted!
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh legacymd5 >/dev/null 2>&1
rc_leg_cached_noflag=$?
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh --allow-legacy-md5 legacybad >/dev/null 2>&1
rc_leg_badmd5=$?
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh --allow-legacy-md5 dashboth >/dev/null 2>&1
rc_dashboth=$?
set -e
[ "$rc_leg_cached_noflag" -ne 0 ] || { echo "FAIL: omitting --allow-legacy-md5 must fail closed even when .done/.legacy_md5 exist" >&2; exit 1; }
[ "$rc_leg_badmd5" -ne 0 ] || { echo "FAIL: --allow-legacy-md5 must still reject mismatched md5" >&2; exit 1; }
[ "$rc_dashboth" -ne 0 ] || { echo "FAIL: --allow-legacy-md5 must reject md5='-'" >&2; exit 1; }

# Upgrading a legacy-md5 row to a new pinned sha256 must NOT trust stale .legacy_md5/.done
cat << EOF > "$TMPDIR/upgrade_assets.tsv"
# group	name	url	md5	sha256	dest
legacymd5	leg1	file://$TMPDIR/src/updated.txt	$upd_md5	$upd_sha	test/legacy.txt
EOF
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/upgrade_assets.tsv" \
  bash harness/fetch_data.sh legacymd5
[ "$(portable_sha256 "$TMPDIR/data/test/legacy.txt")" = "$upd_sha" ] || {
  echo "FAIL: upgrading row from legacy_md5 to new sha256 did not refresh stale file" >&2; exit 1;
}
[ -f "$TMPDIR/data/test/legacy.txt.sha256" ] || { echo "FAIL: upgraded row missing .sha256 marker" >&2; exit 1; }
[ ! -f "$TMPDIR/data/test/legacy.txt.legacy_md5" ] || { echo "FAIL: upgraded row retained stale .legacy_md5" >&2; exit 1; }

# 3e. Existing local $out matching sha256 without $out.done promotes in-place (zero network download)
mkdir -p "$TMPDIR/data/test"
cp "$TMPDIR/src/testfile.txt" "$TMPDIR/data/test/preexisting.txt"
out_promote=$(DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh promotegrp)
grep -q "have: promote1" <<< "$out_promote" || { echo "FAIL: in-place promotion did not print have: promote1" >&2; exit 1; }
[ -f "$TMPDIR/data/test/preexisting.txt.done" ] || { echo "FAIL: in-place promotion did not touch .done" >&2; exit 1; }
[ -f "$TMPDIR/data/test/preexisting.txt.sha256" ] || { echo "FAIL: in-place promotion did not write .sha256" >&2; exit 1; }

# 3f. Prefix collision test: foo.vcf.gz.tbi fetched BEFORE foo.vcf.gz must record BOTH in MANIFEST
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/test_assets.tsv" \
  bash harness/fetch_data.sh prefixgrp
awk -F'\t' '!/^#/ && $1 == "foo.vcf.gz.tbi" {f=1} END{exit !f}' "$TMPDIR/benchmarks/MANIFEST" || {
  echo "FAIL: foo.vcf.gz.tbi missing from MANIFEST" >&2; exit 1;
}
awk -F'\t' '!/^#/ && $1 == "foo.vcf.gz" {f=1} END{exit !f}' "$TMPDIR/benchmarks/MANIFEST" || {
  echo "FAIL: foo.vcf.gz skipped in MANIFEST due to prefix collision with foo.vcf.gz.tbi" >&2; exit 1;
}

# 3g. Unverified dummy transformed refs MUST NOT bypass SHA-256 verification
mkdir -p "$TMPDIR/data/ref/d1vd1" "$TMPDIR/data/ref/giabv3"
printf ">chr1\nACGT\n" > "$TMPDIR/data/ref/d1vd1/GRCh38.d1.vd1.fa"
printf "chr1\t4\t6\t4\t5\n" > "$TMPDIR/data/ref/d1vd1/GRCh38.d1.vd1.fa.fai"
touch "$TMPDIR/data/ref/d1vd1/GRCh38.d1.vd1.fa.done"
printf ">chr1\nACGT\n" > "$TMPDIR/data/ref/giabv3/GRCh38_GIABv3.fasta"
touch "$TMPDIR/data/ref/giabv3/GRCh38_GIABv3.fasta.gz.done"
cat << EOF > "$TMPDIR/refs_assets.tsv"
# group	name	url	md5	sha256	dest
refs	GRCh38.d1.vd1.fa	file://$TMPDIR/src/unreachable_d1vd1.fa	7eb7c875fb968ab47186e01dbcbfe1b1	fbde121e898f00cc064eaafe5acc58940d7fe069ad14d34acd40c5d675e6b69b	ref/d1vd1/GRCh38.d1.vd1.fa
EOF
cat << EOF > "$TMPDIR/refs_giab_assets.tsv"
# group	name	url	md5	sha256	dest
refs	GRCh38_GIABv3.fasta.gz	file://$TMPDIR/src/unreachable_giabv3.fasta.gz	939ce19062d1462c09b88c55faca4d76	b37cb9a9c0f8716fdb13e2c5c60f08cb8be40de39c3088176f1f62d9e8eeeed8	ref/giabv3/GRCh38_GIABv3.fasta.gz
EOF
set +e
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/refs_assets.tsv" \
  bash harness/fetch_data.sh refs >/dev/null 2>&1
rc_fake_d1=$?
DATA_DIR="$TMPDIR/data" MANIFEST="$TMPDIR/benchmarks/MANIFEST" ASSETS_TSV="$TMPDIR/refs_giab_assets.tsv" \
  bash harness/fetch_data.sh refs >/dev/null 2>&1
rc_fake_gv=$?
set -e
[ "$rc_fake_d1" -ne 0 ] || { echo "FAIL: 4-base dummy GRCh38.d1.vd1.fa with .fai/.done bypassed SHA-256 check" >&2; exit 1; }
[ "$rc_fake_gv" -ne 0 ] || { echo "FAIL: 4-base dummy GRCh38_GIABv3.fasta bypassed transformed-ref SHA-256 check" >&2; exit 1; }

# 4. Test stage_g4.sh usage & get_hg008_slice.sh --beds-only SHA-256 validation + r1..r3 mtime/staleness
bash harness/stage_g4.sh --help | awk '/stage_g4\.sh/ { found=1 } END { exit !found }' || { echo "FAIL: stage_g4.sh --help failed" >&2; exit 1; }
set +e
bash harness/stage_g4.sh unknown_target >/dev/null 2>&1
rc_unk=$?
set -e
[ "$rc_unk" -ne 0 ] || { echo "FAIL: stage_g4.sh unknown_target should fail" >&2; exit 1; }

mkdir -p "$TMPDIR/hg008_src" "$TMPDIR/g4/hg008/truth"
cat << 'EOF' | gzip -c > "$TMPDIR/hg008_src/v03.vcf.gz"
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10001	.	A	T	60	PASS	.
chr1	22001	.	C	G	60	PASS	.
chr1	35001	.	G	A	60	PASS	.
EOF
printf "chr1\t5000\t40000\n" > "$TMPDIR/hg008_src/v03_all.bed"
cp "$TMPDIR/hg008_src/v03.vcf.gz" "$TMPDIR/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"
cp "$TMPDIR/hg008_src/v03_all.bed" "$TMPDIR/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed"

hg008_vcf_md5=$(portable_md5 "$TMPDIR/hg008_src/v03.vcf.gz")
hg008_vcf_sha=$(portable_sha256 "$TMPDIR/hg008_src/v03.vcf.gz")
hg008_bed_md5=$(portable_md5 "$TMPDIR/hg008_src/v03_all.bed")
hg008_bed_sha=$(portable_sha256 "$TMPDIR/hg008_src/v03_all.bed")

cat << EOF > "$TMPDIR/hg008_assets.tsv"
# group	name	url	md5	sha256	dest
hg008	HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz	file://$TMPDIR/hg008_src/v03.vcf.gz	$hg008_vcf_md5	$hg008_vcf_sha	hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz
hg008	HG008-T_somatic_smvar_benchmark_v0.3_all.bed	file://$TMPDIR/hg008_src/v03_all.bed	$hg008_bed_md5	$hg008_bed_sha	hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed
EOF

HG008_TRUTH_VCF_SHA256="$hg008_vcf_sha" HG008_TRUTH_BED_SHA256="$hg008_bed_sha" \
  DATA_DIR="$TMPDIR/g4" bash harness/get_hg008_slice.sh --beds-only >/dev/null
[ -s "$TMPDIR/g4/hg008/slices.bed" ] || { echo "FAIL: get_hg008_slice.sh did not create slices.bed" >&2; exit 1; }
[ -s "$TMPDIR/g4/hg008/engine.bed" ] || { echo "FAIL: get_hg008_slice.sh did not create engine.bed" >&2; exit 1; }
[ -s "$TMPDIR/g4/hg008/eval.bed" ]   || { echo "FAIL: get_hg008_slice.sh did not create eval.bed" >&2; exit 1; }
[ -s "$TMPDIR/g4/hg008/r1.bed" ]     || { echo "FAIL: get_hg008_slice.sh did not create r1.bed" >&2; exit 1; }
# Span is 9500..35501 (26001 bp > 20000 bp), so engine.bed must split into 2 tiles
[ "$(wc -l < "$TMPDIR/g4/hg008/engine.bed" | tr -d '[:space:]')" -eq 2 ] || { echo "FAIL: expected 2 tiles in engine.bed" >&2; exit 1; }

# Verify r1.bed mtime is preserved when slices.bed is unchanged, and updated when slices.bed changes
touch -t 202001020304 "$TMPDIR/g4/hg008/r1.bed"
mtime_before=$(python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$TMPDIR/g4/hg008/r1.bed")
HG008_TRUTH_VCF_SHA256="$hg008_vcf_sha" HG008_TRUTH_BED_SHA256="$hg008_bed_sha" \
  DATA_DIR="$TMPDIR/g4" bash harness/get_hg008_slice.sh --beds-only >/dev/null
mtime_same=$(python3 -c 'import os,sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$TMPDIR/g4/hg008/r1.bed")
[ "$mtime_before" = "$mtime_same" ] || { echo "FAIL: get_hg008_slice.sh touched r1.bed mtime when content was identical" >&2; exit 1; }

printf "chr2\t1000\t2000\tslice_override\n" > "$TMPDIR/g4/hg008/slices.bed"
HG008_TRUTH_VCF_SHA256="$hg008_vcf_sha" HG008_TRUTH_BED_SHA256="$hg008_bed_sha" \
  DATA_DIR="$TMPDIR/g4" bash harness/get_hg008_slice.sh --beds-only >/dev/null
grep -q '^chr2	1000	2000	slice_override$' "$TMPDIR/g4/hg008/r1.bed" || {
  echo "FAIL: get_hg008_slice.sh did not refresh stale r1.bed after slices.bed changed" >&2; exit 1;
}

# Non-empty modified cached HG008 truth BED (with stale .sha256 sidecar) MUST be detected, refreshed, and regenerate BEDs!
printf "chr9\t1\t10\n" > "$TMPDIR/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed"
printf '%s\n' "$hg008_bed_sha" > "$TMPDIR/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed.sha256"
HG008_TRUTH_VCF_SHA256="$hg008_vcf_sha" HG008_TRUTH_BED_SHA256="$hg008_bed_sha" \
  ASSETS_TSV="$TMPDIR/hg008_assets.tsv" MANIFEST="$TMPDIR/g4/MANIFEST" \
  DATA_DIR="$TMPDIR/g4" bash harness/get_hg008_slice.sh --beds-only >/dev/null
[ "$(portable_sha256 "$TMPDIR/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed")" = "$hg008_bed_sha" ] || {
  echo "FAIL: get_hg008_slice.sh accepted non-empty modified HG008 truth BED without refreshing" >&2; exit 1;
}
grep -q '^chr1	' "$TMPDIR/g4/hg008/r1.bed" || {
  echo "FAIL: get_hg008_slice.sh did not regenerate derived BEDs after repairing corrupted truth input" >&2; exit 1;
}

# 5. Verify DATA_DIR propagation, stage_verified_sha256 usage, and stage_g4.sh invariants
grep -q 'samtools faidx "\$gv/GRCh38_GIABv3.fasta"' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing samtools faidx on GRCh38_GIABv3.fasta" >&2; exit 1; }
grep -q 'seqc2_align.sh HCC1395 "\$sq/fastq/SRR7890879_1.fastq.gz" "\$sq/fastq/SRR7890879_2.fastq.gz"' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing SRR7890879 FASTQ paths" >&2; exit 1; }
grep -q 'seqc2_align.sh HCC1395BL "\$sq/fastq/SRR7890880_1.fastq.gz" "\$sq/fastq/SRR7890880_2.fastq.gz"' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh missing SRR7890880 FASTQ paths" >&2; exit 1; }
grep -q 'DATA_DIR="\$D" bash harness/get_hg002_slice.sh' harness/stage_g4.sh || { echo "FAIL: stage_g4.sh must pass DATA_DIR=\$D to get_hg002_slice.sh" >&2; exit 1; }
grep -q 'D="\${DATA_DIR:-scratch/g4}"' harness/get_hg002_slice.sh || { echo "FAIL: get_hg002_slice.sh must honour DATA_DIR" >&2; exit 1; }
grep -Eq 'stage_verified_sha256|stream_verified_bam_slice' harness/get_hg002_slice.sh || { echo "FAIL: get_hg002_slice.sh must use stream_verified_bam_slice" >&2; exit 1; }
grep -Eq 'stage_verified_sha256|stream_verified_bam_slice' harness/get_hg008_slice.sh || { echo "FAIL: get_hg008_slice.sh must use stream_verified_bam_slice" >&2; exit 1; }
grep -q 'lib/portable.sh' harness/get_multimodal_slices.sh || { echo "FAIL: get_multimodal_slices.sh must source portable.sh" >&2; exit 1; }
grep -Eq 'stage_verified_sha256|stream_verified_bam_slice' harness/get_multimodal_slices.sh || { echo "FAIL: get_multimodal_slices.sh must use stream_verified_bam_slice" >&2; exit 1; }

# 6. Regression test: get_hg002_slice.sh refreshes corrupted local truth BED/VCF/TBI even when stale .sha256/.done sidecars exist
mkdir -p "$TMPDIR/mockbin" "$TMPDIR/hg002_src" "$TMPDIR/g4_hg002/hg002/truth"
cat << 'EOF' > "$TMPDIR/mockbin/samtools"
#!/usr/bin/env bash
exit 0
EOF
cat << 'EOF' > "$TMPDIR/mockbin/bcftools"
#!/usr/bin/env bash
if [ "${1:-}" = "view" ] && [ "${2:-}" = "-H" ]; then
  printf 'chr20\t10000100\t.\tA\tG\t60\tPASS\t.\n'
  exit 0
fi
shift # drop 'view'
src=""
out=""
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    -r) shift 2 ;;
    -*) shift ;;
    *)  src="$1"; shift ;;
  esac
done
cp "$src" "$out"
EOF
cat << 'EOF' > "$TMPDIR/mockbin/tabix"
#!/usr/bin/env bash
touch "${!#}.tbi"
EOF
chmod +x "$TMPDIR/mockbin/samtools" "$TMPDIR/mockbin/bcftools" "$TMPDIR/mockbin/tabix"

# Valid cached chr20 BAM + BAI so Step 1 runs 100% offline
printf 'BAM' > "$TMPDIR/g4_hg002/hg002/HG002.chr20.bam"
printf 'BAI' > "$TMPDIR/g4_hg002/hg002/HG002.chr20.bam.bai"

# Valid upstream mock truth BED (50,000,000 bp span -> >1000 50kb tiles over chr20:10M-30M -> 400 tiles per interval, 3 intervals = 1200 tiles)
cat << 'EOF' > "$TMPDIR/hg002_src/bench.bed"
chr20	10000000	30000000
chr20	10000000	30000000
chr20	10000000	30000000
EOF
printf 'MOCK_VCF_GZ_PAYLOAD\n' > "$TMPDIR/hg002_src/bench.vcf.gz"
printf 'MOCK_TBI_PAYLOAD\n' > "$TMPDIR/hg002_src/bench.vcf.gz.tbi"

m_bed_md5=$(portable_md5 "$TMPDIR/hg002_src/bench.bed")
m_bed_sha=$(portable_sha256 "$TMPDIR/hg002_src/bench.bed")
m_vcf_md5=$(portable_md5 "$TMPDIR/hg002_src/bench.vcf.gz")
m_vcf_sha=$(portable_sha256 "$TMPDIR/hg002_src/bench.vcf.gz")
m_tbi_md5=$(portable_md5 "$TMPDIR/hg002_src/bench.vcf.gz.tbi")
m_tbi_sha=$(portable_sha256 "$TMPDIR/hg002_src/bench.vcf.gz.tbi")

cat << EOF > "$TMPDIR/hg002_assets.tsv"
# group	name	url	md5	sha256	dest
hg002	HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz	file://$TMPDIR/hg002_src/bench.vcf.gz	$m_vcf_md5	$m_vcf_sha	hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
hg002	HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi	file://$TMPDIR/hg002_src/bench.vcf.gz.tbi	$m_tbi_md5	$m_tbi_sha	hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
hg002	HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed	file://$TMPDIR/hg002_src/bench.bed	$m_bed_md5	$m_bed_sha	hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
EOF

# Seed CORRUPTED local truth files alongside seemingly valid .sha256 + .done sidecars!
corrupt_bed="$TMPDIR/g4_hg002/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
corrupt_vcf="$TMPDIR/g4_hg002/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
corrupt_tbi="${corrupt_vcf}.tbi"
printf 'CORRUPT_BED\n' > "$corrupt_bed"
printf '%s\n' "$m_bed_sha" > "${corrupt_bed}.sha256"
touch "${corrupt_bed}.done"
printf 'CORRUPT_VCF\n' > "$corrupt_vcf"
printf '%s\n' "$m_vcf_sha" > "${corrupt_vcf}.sha256"
touch "${corrupt_vcf}.done"
printf 'CORRUPT_TBI\n' > "$corrupt_tbi"
printf '%s\n' "$m_tbi_sha" > "${corrupt_tbi}.sha256"
touch "${corrupt_tbi}.done"

PATH="$TMPDIR/mockbin:$PATH" \
  DATA_DIR="$TMPDIR/g4_hg002" \
  ASSETS_TSV="$TMPDIR/hg002_assets.tsv" \
  MANIFEST="$TMPDIR/g4_hg002/MANIFEST" \
  HG002_BED_SHA256="$m_bed_sha" \
  HG002_VCF_SHA256="$m_vcf_sha" \
  HG002_TBI_SHA256="$m_tbi_sha" \
  bash harness/get_hg002_slice.sh >/dev/null

[ "$(portable_sha256 "$corrupt_bed")" = "$m_bed_sha" ] || { echo "FAIL: corrupted HG002 BED was not refreshed" >&2; exit 1; }
[ "$(portable_sha256 "$corrupt_vcf")" = "$m_vcf_sha" ] || { echo "FAIL: corrupted HG002 VCF was not refreshed" >&2; exit 1; }
[ "$(portable_sha256 "$corrupt_tbi")" = "$m_tbi_sha" ] || { echo "FAIL: corrupted HG002 TBI was not refreshed" >&2; exit 1; }
[ -s "$TMPDIR/g4_hg002/hg002/eval.bed" ] || { echo "FAIL: get_hg002_slice.sh did not produce eval.bed after repair" >&2; exit 1; }
[ -s "$TMPDIR/g4_hg002/hg002/truth/truth.vcf.gz" ] || { echo "FAIL: get_hg002_slice.sh did not produce truth.vcf.gz after repair" >&2; exit 1; }

echo "PASS test_fetch_data"
