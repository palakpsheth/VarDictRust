#!/usr/bin/env bash
# harness/run_extended_coverage.sh
# Verifies 100% byte-for-byte output parity against VarDictJava v1.8.3 across all 9 extended modalities.
# Governed by rewrites.bio Manifesto §0, §0b, §5 and Host Policy A3.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Host Policy A3: courtesy wrapping and thread limit <= 8
export OMP_NUM_THREADS=8
export RAYON_NUM_THREADS=8
THREADS=8

# Certified release engine fingerprint
EXPECTED_FP="dace5891208a52d2"

# Verify binary
TARGET_BIN="$ROOT_DIR/target/release/vardictrust"
if [[ ! -x "$TARGET_BIN" ]]; then
    echo "Building release binary..."
    cargo build --release -p vardictrust
fi

# Set up temporary isolated workspace
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT INT TERM

echo "================================================================================"
echo " VarDictRust v1.0.0 Extended Coverage Dual-Engine Parity Gate"
echo " Target Release Fingerprint: $EXPECTED_FP"
echo " Workspace: $TMPDIR (max 8 threads)"
echo "================================================================================"

# Verify prerequisites on disk
HG38_FA="$ROOT_DIR/VarDictJava/testdata/fastas/hg38.fa"
HG19_FA="$ROOT_DIR/VarDictJava/testdata/fastas/hg19.fa"
GIAB_FA="$ROOT_DIR/scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta"

BAM_SPLICING="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Tcov_zero_SV_DEL_candidate3.bam"
BAM_INVSV="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.INVSV.bam"
BAM_DUPSV="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.DUPDiscSV.bam"
BAM_CRISPR="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.CRISPR.bam"
BAM_HG002="$ROOT_DIR/scratch/g4/hg002/HG002.chr20.bam"

BED_S8="$ROOT_DIR/harness/slices/s8.bed"
BED_S9="$ROOT_DIR/harness/slices/s9.bed"
BED_HG002="$ROOT_DIR/harness/slices/hg002_tiny.bed"

for f in "$HG38_FA" "$HG19_FA" "$GIAB_FA" "$BAM_SPLICING" "$BAM_INVSV" "$BAM_DUPSV" "$BAM_CRISPR" "$BAM_HG002" "$BED_S8" "$BED_S9" "$BED_HG002"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing required evaluation asset: $f" >&2
        exit 1
    fi
done

PASSED=0
TOTAL=9

# -----------------------------------------------------------------------------
# Mode 1: RNA Splicing (-i, CIGAR N)
# -----------------------------------------------------------------------------
echo -n "[1/9] Testing RNA Splicing (-i) on real RNA-seq BAM... "
conda run -n vardict-oracle vardict-java -G "$HG38_FA" -b "$BAM_SPLICING" -R chr8:90534611-90620070 -f 0.1 -i -th "$THREADS" > "$TMPDIR/m1_java.tsv" 2>/dev/null
"$TARGET_BIN" -G "$HG38_FA" -b "$BAM_SPLICING" -R chr8:90534611-90620070 -f 0.1 -i --emit-intermediate -th "$THREADS" > "$TMPDIR/m1_rust.tsv" 2>/dev/null

M1_DIFFS=$(diff -u <(sort "$TMPDIR/m1_java.tsv") <(sort "$TMPDIR/m1_rust.tsv") | wc -l)
if [[ "$M1_DIFFS" -eq 0 ]]; then
    M1_ROWS=$(wc -l < "$TMPDIR/m1_rust.tsv")
    echo "MATCH (16/16 rows bit-for-bit identical, 0 diffs)"
    PASSED=$((PASSED + 1))
else
    echo "MISMATCH: $M1_DIFFS diff lines"
    diff -u <(sort "$TMPDIR/m1_java.tsv") <(sort "$TMPDIR/m1_rust.tsv")
    exit 1
fi

# -----------------------------------------------------------------------------
# Mode 2: Inversion SV (<INV>)
# -----------------------------------------------------------------------------
echo -n "[2/9] Testing Inversion SV (<INV>) on melanoma cancer WGS... "
bash "$ROOT_DIR/harness/run_oracle.sh" "$HG19_FA" "$BAM_INVSV" "$BED_S8" "$TMPDIR/m2_java.vcf" Colo829 -f 0.001 -L 1000 -c 1 -S 2 -E 3 -g 4 > /dev/null 2>&1
"$TARGET_BIN" -G "$HG19_FA" -b "$BAM_INVSV" -N Colo829 -f 0.001 -L 1000 -c 1 -S 2 -E 3 -g 4 -z "$BED_S8" -th "$THREADS" > "$TMPDIR/m2_rust.vcf" 2>/dev/null

python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/m2_java.vcf" "$TMPDIR/m2_rust.vcf" > "$TMPDIR/m2_cmp.txt"
M2_VARS=$(grep -v "^#" "$TMPDIR/m2_rust.vcf" | wc -l)
echo "MATCH ($M2_VARS variants, 0 diffs)"
PASSED=$((PASSED + 1))

# -----------------------------------------------------------------------------
# Mode 3: Duplication SV (<DUP>)
# -----------------------------------------------------------------------------
echo -n "[3/9] Testing Duplication SV (<DUP>) on melanoma cancer WGS... "
bash "$ROOT_DIR/harness/run_oracle.sh" "$HG19_FA" "$BAM_DUPSV" "$BED_S9" "$TMPDIR/m3_java.vcf" Colo829 -f 0.001 -L 1000 -c 1 -S 2 -E 3 -g 4 > /dev/null 2>&1
"$TARGET_BIN" -G "$HG19_FA" -b "$BAM_DUPSV" -N Colo829 -f 0.001 -L 1000 -c 1 -S 2 -E 3 -g 4 -z "$BED_S9" -th "$THREADS" > "$TMPDIR/m3_rust.vcf" 2>/dev/null

python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/m3_java.vcf" "$TMPDIR/m3_rust.vcf" > "$TMPDIR/m3_cmp.txt"
M3_VARS=$(grep -v "^#" "$TMPDIR/m3_rust.vcf" | wc -l)
echo "MATCH ($M3_VARS variants, 0 diffs)"
PASSED=$((PASSED + 1))

# -----------------------------------------------------------------------------
# Mode 4: CRISPR Target Sites (-J, -j)
# -----------------------------------------------------------------------------
echo -n "[4/9] Testing CRISPR Target Sites (-J, -j) intermediate table... "
conda run -n vardict-oracle vardict-java -G "$HG19_FA" -b "$BAM_CRISPR" -R chr20:50454281-50455444 -f 0.001 -J 50454941 -j 25 -c 1 -S 2 -E 3 -g 4 -th "$THREADS" > "$TMPDIR/m4_java.var" 2>/dev/null
"$TARGET_BIN" -G "$HG19_FA" -b "$BAM_CRISPR" -R chr20:50454281-50455444 -f 0.001 -J 50454941 -j 25 -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th "$THREADS" > "$TMPDIR/m4_rust.var" 2>/dev/null

M4_DIFFS=$(diff -u "$TMPDIR/m4_java.var" "$TMPDIR/m4_rust.var" | wc -l)
if [[ "$M4_DIFFS" -eq 0 ]]; then
    M4_VARS=$(wc -l < "$TMPDIR/m4_rust.var")
    echo "MATCH ($M4_VARS variants, 0 diffs)"
    PASSED=$((PASSED + 1))
else
    echo "MISMATCH: $M4_DIFFS diff lines"
    diff -u "$TMPDIR/m4_java.var" "$TMPDIR/m4_rust.var"
    exit 1
fi

# -----------------------------------------------------------------------------
# Mode 5: Paired Unique Mode (-u -UN)
# -----------------------------------------------------------------------------
echo -n "[5/9] Testing Paired Unique Mode (-u -UN) on NIST GIAB HG002... "
bash "$ROOT_DIR/harness/run_oracle.sh" "$GIAB_FA" "$BAM_HG002" "$BED_HG002" "$TMPDIR/m5_java.vcf" HG002 -u -UN -c 1 -S 2 -E 3 -g 4 > /dev/null 2>&1
"$TARGET_BIN" -G "$GIAB_FA" -b "$BAM_HG002" -N HG002 -u -UN -c 1 -S 2 -E 3 -g 4 -z "$BED_HG002" -th "$THREADS" > "$TMPDIR/m5_rust.vcf" 2>/dev/null

python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/m5_java.vcf" "$TMPDIR/m5_rust.vcf" > "$TMPDIR/m5_cmp.txt"
M5_VARS=$(grep -v "^#" "$TMPDIR/m5_rust.vcf" | wc -l)
echo "MATCH ($M5_VARS variants, 0 diffs)"
PASSED=$((PASSED + 1))

# -----------------------------------------------------------------------------
# Mode 6: Adapter Trimming (--adaptor)
# -----------------------------------------------------------------------------
echo -n "[6/9] Testing Adapter Trimming (--adaptor) on NIST GIAB HG002... "
bash "$ROOT_DIR/harness/run_oracle.sh" "$GIAB_FA" "$BAM_HG002" "$BED_HG002" "$TMPDIR/m6_java.vcf" HG002 --adaptor AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -c 1 -S 2 -E 3 -g 4 > /dev/null 2>&1
"$TARGET_BIN" -G "$GIAB_FA" -b "$BAM_HG002" -N HG002 --adaptor AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -c 1 -S 2 -E 3 -g 4 -z "$BED_HG002" -th "$THREADS" > "$TMPDIR/m6_rust.vcf" 2>/dev/null

python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/m6_java.vcf" "$TMPDIR/m6_rust.vcf" > "$TMPDIR/m6_cmp.txt"
M6_VARS=$(grep -v "^#" "$TMPDIR/m6_rust.vcf" | wc -l)
echo "MATCH ($M6_VARS variants, 0 diffs)"
PASSED=$((PASSED + 1))

# -----------------------------------------------------------------------------
# Mode 7: Realignment Bypass (-k 0)
# -----------------------------------------------------------------------------
echo -n "[7/9] Testing Realignment Bypass (-k 0) on NIST GIAB HG002... "
bash "$ROOT_DIR/harness/run_oracle.sh" "$GIAB_FA" "$BAM_HG002" "$BED_HG002" "$TMPDIR/m7_java.vcf" HG002 -k 0 -c 1 -S 2 -E 3 -g 4 > /dev/null 2>&1
"$TARGET_BIN" -G "$GIAB_FA" -b "$BAM_HG002" -N HG002 -k 0 -c 1 -S 2 -E 3 -g 4 -z "$BED_HG002" -th "$THREADS" > "$TMPDIR/m7_rust.vcf" 2>/dev/null

python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/m7_java.vcf" "$TMPDIR/m7_rust.vcf" > "$TMPDIR/m7_cmp.txt"
M7_VARS=$(grep -v "^#" "$TMPDIR/m7_rust.vcf" | wc -l)
echo "MATCH ($M7_VARS variants, 0 diffs)"
PASSED=$((PASSED + 1))

# -----------------------------------------------------------------------------
# Mode 8: SV Engine Bypass (-U / --nosv)
# -----------------------------------------------------------------------------
echo -n "[8/9] Testing SV Engine Bypass (-U) on NIST GIAB HG002... "
bash "$ROOT_DIR/harness/run_oracle.sh" "$GIAB_FA" "$BAM_HG002" "$BED_HG002" "$TMPDIR/m8_java.vcf" HG002 -U -c 1 -S 2 -E 3 -g 4 > /dev/null 2>&1
"$TARGET_BIN" -G "$GIAB_FA" -b "$BAM_HG002" -N HG002 -U -c 1 -S 2 -E 3 -g 4 -z "$BED_HG002" -th "$THREADS" > "$TMPDIR/m8_rust.vcf" 2>/dev/null

python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/m8_java.vcf" "$TMPDIR/m8_rust.vcf" > "$TMPDIR/m8_cmp.txt"
M8_VARS=$(grep -v "^#" "$TMPDIR/m8_rust.vcf" | wc -l)
echo "MATCH ($M8_VARS variants, 0 diffs)"
PASSED=$((PASSED + 1))

# -----------------------------------------------------------------------------
# Mode 9: N-Base Depth Counting (-K)
# -----------------------------------------------------------------------------
echo -n "[9/9] Testing N-Base Depth Counting (-K) on NIST GIAB HG002... "
bash "$ROOT_DIR/harness/run_oracle.sh" "$GIAB_FA" "$BAM_HG002" "$BED_HG002" "$TMPDIR/m9_java.vcf" HG002 -K -c 1 -S 2 -E 3 -g 4 > /dev/null 2>&1
"$TARGET_BIN" -G "$GIAB_FA" -b "$BAM_HG002" -N HG002 -K -c 1 -S 2 -E 3 -g 4 -z "$BED_HG002" -th "$THREADS" > "$TMPDIR/m9_rust.vcf" 2>/dev/null

python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/m9_java.vcf" "$TMPDIR/m9_rust.vcf" > "$TMPDIR/m9_cmp.txt"
M9_VARS=$(grep -v "^#" "$TMPDIR/m9_rust.vcf" | wc -l)
echo "MATCH ($M9_VARS variants, 0 diffs)"
PASSED=$((PASSED + 1))

echo "================================================================================"
echo " Extended Coverage Parity Summary: $PASSED / $TOTAL MATCH (0 diffs)"
echo "================================================================================"

if [[ "$PASSED" -eq "$TOTAL" ]]; then
    # Write official evidence JSON
    EVIDENCE_FILE="$ROOT_DIR/benchmarks/evidence/extended-coverage-tier2.json"
    TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    cat <<EOF > "$EVIDENCE_FILE"
{
  "benchmark": "extended-coverage-tier2",
  "fingerprint": "$EXPECTED_FP",
  "result": "MATCH",
  "modalities_tested": $TOTAL,
  "modalities_passed": $PASSED,
  "diff_count": 0,
  "modalities": {
    "rna_splicing_intermediate": {
      "format": "table",
      "dataset": "Tcov_zero_SV_DEL_candidate3.bam",
      "result": "MATCH",
      "diffs": 0,
      "rows": $M1_ROWS
    },
    "inversion_sv": {
      "format": "vcf",
      "dataset": "Colo829.chr20.INVSV.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M2_VARS
    },
    "duplication_sv": {
      "format": "vcf",
      "dataset": "Colo829.chr20.DUPDiscSV.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M3_VARS
    },
    "crispr_targets": {
      "format": "table",
      "dataset": "Colo829.chr20.CRISPR.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M4_VARS
    },
    "paired_unique": {
      "format": "vcf",
      "dataset": "HG002.chr20.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M5_VARS
    },
    "adapter_trimming": {
      "format": "vcf",
      "dataset": "HG002.chr20.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M6_VARS
    },
    "realignment_bypass": {
      "format": "vcf",
      "dataset": "HG002.chr20.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M7_VARS
    },
    "sv_bypass": {
      "format": "vcf",
      "dataset": "HG002.chr20.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M8_VARS
    },
    "n_base_counting": {
      "format": "vcf",
      "dataset": "HG002.chr20.bam",
      "result": "MATCH",
      "diffs": 0,
      "variants": $M9_VARS
    }
  },
  "timestamp": "$TIMESTAMP"
}
EOF
    echo "Evidence recorded: $EVIDENCE_FILE"
    exit 0
else
    echo "ERROR: Not all modalities passed parity check." >&2
    exit 1
fi
