#!/usr/bin/env bash
set -euo pipefail

if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
  exec conda run -n g4-tools bash "$0" "$@"
fi

TMPDIR=$(mktemp -d -t test_concordance_XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT
export GOLDEN_STORE_DIR="$TMPDIR/golden"

REF="fixtures/micro/micro.ref.fa"

# 1. Build SDF
rtg format -o "$TMPDIR/sdf" "$REF" >/dev/null

# 2. Get contig info from fai
CONTIG="chr1"
L=300

p1=$((L / 4))
p2=$((L / 2))
p3=$((3 * L / 4))

get_alt() {
  local ref_base=$1
  if [ "$ref_base" = "A" ]; then echo "C"; else echo "A"; fi
}

ref1=$(samtools faidx "$REF" "${CONTIG}:${p1}-${p1}" | tail -n +2 | tr -d '\r\n')
ref2=$(samtools faidx "$REF" "${CONTIG}:${p2}-${p2}" | tail -n +2 | tr -d '\r\n')
ref3=$(samtools faidx "$REF" "${CONTIG}:${p3}-${p3}" | tail -n +2 | tr -d '\r\n')

alt1=$(get_alt "$ref1")
alt2=$(get_alt "$ref2")
alt3=$(get_alt "$ref3")

# 3. Write truth VCF (sites-only)
cat << EOF > "$TMPDIR/truth.vcf"
##fileformat=VCFv4.2
##contig=<ID=chr1,length=300>
##contig=<ID=chr2,length=260>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
${CONTIG}	${p1}	.	${ref1}	${alt1}	.	PASS	.
${CONTIG}	${p2}	.	${ref2}	${alt2}	.	PASS	.
${CONTIG}	${p3}	.	${ref3}	${alt3}	.	PASS	.
EOF

bcftools view -Oz -o "$TMPDIR/truth.vcf.gz" "$TMPDIR/truth.vcf"
tabix -p vcf "$TMPDIR/truth.vcf.gz"

# 4. Write calls VCF
cat << EOF > "$TMPDIR/calls.vcf"
##fileformat=VCFv4.2
##contig=<ID=chr1,length=300>
##contig=<ID=chr2,length=260>
##INFO=<ID=STATUS,Number=1,Type=String,Description="Somatic status">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL
${CONTIG}	${p1}	.	${ref1}	${alt1}	100	PASS	STATUS=StrongSomatic	GT	0/1	0/0
${CONTIG}	${p2}	.	${ref2}	${alt2}	100	PASS	STATUS=StrongSomatic	GT	0/1	0/0
${CONTIG}	${p3}	.	${ref3}	${alt3}	100	PASS	STATUS=StrongSomatic	GT	0/1	0/0
EOF

# 5. Write eval BED
printf "%s\t0\t%d\n" "$CONTIG" "$L" > "$TMPDIR/eval.bed"

# 6. Run concordance on full set -> expect tp_baseline=3, fp=0, fn=0
bash harness/concordance.sh "$TMPDIR/calls.vcf" "$TMPDIR/truth.vcf.gz" "$TMPDIR/eval.bed" "$TMPDIR/sdf" "$TMPDIR/out" TUMOR >/dev/null

python3 - "$TMPDIR/out/summary.json" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
assert data["tp_baseline"] == 3, f"expected tp_baseline=3, got {data['tp_baseline']}"
assert data["fp"] == 0, f"expected fp=0, got {data['fp']}"
assert data["fn"] == 0, f"expected fn=0, got {data['fn']}"
PY

# 7. Drop one call -> expect fn=1
cat << EOF > "$TMPDIR/calls_drop1.vcf"
##fileformat=VCFv4.2
##contig=<ID=chr1,length=300>
##contig=<ID=chr2,length=260>
##INFO=<ID=STATUS,Number=1,Type=String,Description="Somatic status">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL
${CONTIG}	${p1}	.	${ref1}	${alt1}	100	PASS	STATUS=StrongSomatic	GT	0/1	0/0
${CONTIG}	${p2}	.	${ref2}	${alt2}	100	PASS	STATUS=StrongSomatic	GT	0/1	0/0
EOF

bash harness/concordance.sh "$TMPDIR/calls_drop1.vcf" "$TMPDIR/truth.vcf.gz" "$TMPDIR/eval.bed" "$TMPDIR/sdf" "$TMPDIR/out_drop1" TUMOR >/dev/null

python3 - "$TMPDIR/out_drop1/summary.json" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
assert data["fn"] == 1, f"expected fn=1, got {data['fn']}"
PY

# 8. Empty truth set -> expect non-zero exit
cat << EOF > "$TMPDIR/truth_empty.vcf"
##fileformat=VCFv4.2
##contig=<ID=chr1,length=300>
##contig=<ID=chr2,length=260>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

bcftools view -Oz -o "$TMPDIR/truth_empty.vcf.gz" "$TMPDIR/truth_empty.vcf"
tabix -p vcf "$TMPDIR/truth_empty.vcf.gz"

set +e
bash harness/concordance.sh "$TMPDIR/calls.vcf" "$TMPDIR/truth_empty.vcf.gz" "$TMPDIR/eval.bed" "$TMPDIR/sdf" "$TMPDIR/out_empty" TUMOR >/dev/null 2>&1
empty_rc=$?
set -e

[ "$empty_rc" -ne 0 ] || { echo "FAIL test_concordance: expected non-zero exit on empty truth" >&2; exit 1; }

echo "PASS test_concordance"
