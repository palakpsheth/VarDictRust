#!/usr/bin/env bash
set -euo pipefail
# Self-test for harness/concordance_sv.sh (SV Tier-3 gated spike).
T=$(mktemp -d -t test_concordance_sv_XXXXXX); trap 'rm -rf "$T"' EXIT; fail=0

# 1. Usage line: wrong arg count exits 2 (no truvari needed).
set +e
bash harness/concordance_sv.sh only-one-arg >/dev/null 2>&1; rc=$?
set -e
[ "$rc" -eq 2 ] || { echo "FAIL: usage exit=$rc, want 2"; fail=1; }

# 2. Plain-VCF contract: bgzipped+indexed required — loud rejection, no
# summary.json (validation precedes the truvari check, so this holds with
# or without truvari installed).
echo '##fileformat=VCFv4.2' > "$T/plain.vcf"
set +e
plain_out=$(bash harness/concordance_sv.sh "$T/plain.vcf" "$T/t.vcf.gz" "$T/r.bed" "$T/out_plain" 2>&1); rc=$?
set -e
[ "$rc" -eq 2 ] || { echo "FAIL: plain-vcf exit=$rc, want 2"; fail=1; }
[[ "$plain_out" == *"calls must be bgzipped + indexed VCF (samtools bgzip + tabix), got: "* ]] || { echo "FAIL: plain-vcf msg: $plain_out"; fail=1; }
[ ! -e "$T/out_plain/summary.json" ] || { echo "FAIL: summary.json written on plain-vcf rejection"; fail=1; }

# 3. Missing-truvari path (stub PATH hides any live truvari): SKIP + exit 0,
# outdir untouched (loud skip, never silent). Uses the valid .vcf.gz contract
# (dummy + index) so it reaches the SKIP branch.
S="$T/nostub"; mkdir -p "$S"
BASHBIN=$(command -v bash)
: > "$T/c.vcf.gz"; : > "$T/c.vcf.gz.tbi"
set +e
skip_out=$(PATH="$S" "$BASHBIN" harness/concordance_sv.sh "$T/c.vcf.gz" "$T/t.vcf.gz" "$T/r.bed" "$T/out_skip" 2>&1); rc=$?
set -e
[ "$rc" -eq 0 ] || { echo "FAIL: missing-truvari exit=$rc, want 0"; fail=1; }
[[ "$skip_out" == SKIP:*"truvari not installed"* ]] || { echo "FAIL: missing-truvari msg: $skip_out"; fail=1; }
[ ! -e "$T/out_skip" ] || { echo "FAIL: outdir touched on SKIP path"; fail=1; }

# 4. Ran path: only where truvari is actually installed; else SKIP cleanly.
# v5.x bench requires bgzip'd + indexed inputs, so compress first (read-only
# reuse of bgzip/tabix); SKIP if either half of the toolchain is missing.
if ! command -v truvari >/dev/null 2>&1; then
  echo "SKIP: truvari not installed (ran path untested)"
  [ $fail -eq 0 ] && echo "PASS test_concordance_sv"; exit $fail
fi
if ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1; then
  echo "SKIP: bgzip/tabix missing (ran path untested)"
  [ $fail -eq 0 ] && echo "PASS test_concordance_sv"; exit $fail
fi
cat > "$T/sv.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	101	.	N	<DEL>	.	PASS	SVTYPE=DEL;END=200;SVLEN=-100
EOF
printf 'chr1\t0\t1000\n' > "$T/r.bed"
bgzip -c "$T/sv.vcf" > "$T/sv.vcf.gz" && tabix -p vcf "$T/sv.vcf.gz" \
  || { echo "FAIL: bgzip/tabix fixture prep failed"; fail=1; }
bash harness/concordance_sv.sh "$T/sv.vcf.gz" "$T/sv.vcf.gz" "$T/r.bed" "$T/out_ran" >"$T/ran.stdout" \
  || { echo "FAIL: ran path exit $?"; fail=1; }
python3 - "$T/out_ran/summary.json" <<'PY' || { echo "FAIL: summary.json invalid"; fail=1; }
import json, sys
d = json.load(open(sys.argv[1]))
assert d["tool"] == "truvari", d
assert isinstance(d["tool_version"], str) and d["tool_version"].strip(), d
assert d["status"] == "ran", d
assert set(d) == {"tool", "tool_version", "status"}, d
PY

[ $fail -eq 0 ] && echo "PASS test_concordance_sv"; exit $fail
