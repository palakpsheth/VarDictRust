#!/usr/bin/env bash
# =============================================================================
# harness/run_multimodal_slices.sh
# Dual-engine validation runner across 18 non-overlapping biological slices
# from 6 diverse public sequencing datasets spanning 9 distinct chromosomes.
#
# Adheres to Host Policy A3:
# - Wraps execution with resource limits (max 8 threads)
# - Produces benchmarks/evidence/multimodal-slices-tier2.json
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$ROOT_DIR/scratch/g4/multimodal"
EVIDENCE_JSON="$ROOT_DIR/benchmarks/evidence/multimodal-slices-tier2.json"
TARGET_BIN="$ROOT_DIR/target/release/vardictrust"
THREADS=8

# Reference FASTAs
FA_B37="$ROOT_DIR/VarDictJava/testdata/fastas/GRCh37.fa"
FA_GIAB="$ROOT_DIR/scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta"
FA_SEQC2="$ROOT_DIR/scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa"
FA_HG38="$ROOT_DIR/VarDictJava/testdata/fastas/hg38.fa"

# Auto-materialize test FASTAs if missing
if [[ ! -f "$FA_B37" || ! -f "$FA_HG38" ]]; then
    python3 "$ROOT_DIR/harness/convert_fastas.py"
fi

CURRENT_FP=$(bash "$ROOT_DIR/harness/engine_fingerprint.sh")

# Advisory vardictcpp third arm (never gating): staple the conda pin presence
# into the evidence JSON for the ev_cpp presence gate (per-cell CPP lines in
# the log carry match/skip detail; this key only asserts the arm was wired).
# Canonical conda-pin probe (CPP_PIN_STR): keep pin format in sync with harness/run_extended_coverage.sh and harness/staple_cpp_extra.sh (cpp_val); hardcoded twin: harness/concordance_germline.sh + concordance_pair.sh.
if CPP_PIN_STR=$(conda list -n vardict-oracle vardictcpp 2>/dev/null | awk '/^vardictcpp[[:space:]]/{print "vardictcpp " $2 " " $3 " bioconda"}') && [ -n "$CPP_PIN_STR" ]; then
    CPP_JSON_LINE="'extra': {'cpp_version': '$CPP_PIN_STR'},"
else
    CPP_JSON_LINE="'extra': {'cpp_skip_reason': 'vardictcpp not available in conda env vardict-oracle'},"
fi

TMPDIR=$(mktemp -d /tmp/vardict_multimodal_XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

echo "================================================================================"
echo " VarDictRust v1.0.0 Multi-Modal Multi-Slice Dual-Engine Validation Suite"
echo " 18 Non-Overlapping Slices Across 6 Public Datasets & 9 Chromosomes"
echo " Release Engine Fingerprint: $CURRENT_FP"
echo " Temp Directory: $TMPDIR"
echo "================================================================================"


mkdir -p "$(dirname "$EVIDENCE_JSON")"

PASSED_COUNT=0
TOTAL_SLICES=18
TOTAL_VARIANTS=0
declare -a SLICE_RESULTS=()

export CMP_STRICT=1

function run_vcf_compare() {
    local slice_idx="$1"
    local slice_id="$2"
    local dataset="$3"
    local ref="$4"
    local bam_spec="$5"
    local bed_file="$6"
    local sample="$7"
    shift 7
    local extra_args=("$@")

    echo -n "[$slice_idx/$TOTAL_SLICES] Testing $dataset ($slice_id)... "

    local out_j="$TMPDIR/${slice_id}_java.vcf"
    local out_r="$TMPDIR/${slice_id}_rust.vcf"

    local start_t=$(date +%s%N)

    # 1. Oracle
    bash "$ROOT_DIR/harness/run_oracle.sh" "$ref" "$bam_spec" "$bed_file" "$out_j" "$sample" -z -th "$THREADS" "${extra_args[@]}" 2>/dev/null || {
        echo "FAIL (Oracle error)"
        exit 1
    }

    # 2. Rust
    "$TARGET_BIN" -G "$ref" -N "$sample" -b "$bam_spec" -c 1 -S 2 -E 3 -g 4 -z -th "$THREADS" "${extra_args[@]}" "$bed_file" > "$out_r" 2>/dev/null || {
        echo "FAIL (Rust error)"
        exit 1
    }

    local end_t=$(date +%s%N)
    local elapsed_ms=$(( (end_t - start_t) / 1000000 ))

    # 3. Comparator
    local cmp_out
    if cmp_out=$(python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$out_j" "$out_r" --eps 1e-6 2>&1); then
        local var_count=$(grep -v "^#" "$out_r" | wc -l)
        echo "MATCH ($var_count variants, 0 diffs, ${elapsed_ms}ms)"
        PASSED_COUNT=$((PASSED_COUNT + 1))
        TOTAL_VARIANTS=$((TOTAL_VARIANTS + var_count))
        local bed_str=$(cat "$bed_file" | tr '\t' ' ' | head -n 1)
        SLICE_RESULTS+=("{\"slice_id\": \"$slice_id\", \"dataset\": \"$dataset\", \"region\": \"$bed_str\", \"variants\": $var_count, \"status\": \"MATCH\", \"time_ms\": $elapsed_ms}")
    else
        echo "MISMATCH"
        echo "$cmp_out"
        diff -u "$out_j" "$out_r" | head -40 || true
        exit 1
    fi
}

function run_rna_slice() {
    local slice_idx="$1"
    local slice_id="$2"
    local dataset="$3"
    local ref="$4"
    local bam="$5"
    local region="$6"

    echo -n "[$slice_idx/$TOTAL_SLICES] Testing $dataset Splicing ($slice_id, $region)... "

    local out_j="$TMPDIR/${slice_id}_java.tsv"
    local out_r="$TMPDIR/${slice_id}_rust.tsv"

    local start_t=$(date +%s%N)

    # 1. Oracle Execution (4-column intron frequency table)
    conda run -n vardict-oracle vardict-java -G "$ref" -b "$bam" -R "$region" -f 0.1 -i -th "$THREADS" > "$out_j" 2>/dev/null || {
        echo "FAIL (Oracle RNA error)"
        exit 1
    }

    # 2. VarDictRust Execution
    "$TARGET_BIN" -G "$ref" -b "$bam" -R "$region" -f 0.1 -i --emit-intermediate -th "$THREADS" > "$out_r" 2>/dev/null || {
        echo "FAIL (VarDictRust RNA error)"
        exit 1
    }

    local end_t=$(date +%s%N)
    local elapsed_ms=$(( (end_t - start_t) / 1000000 ))

    local diff_lines=$(diff -u <(sort "$out_j") <(sort "$out_r") | wc -l)
    if [[ "$diff_lines" -eq 0 ]]; then
        local rows=$(wc -l < "$out_r")
        echo "MATCH ($rows splicing junctions, 0 diffs, ${elapsed_ms}ms)"
        PASSED_COUNT=$((PASSED_COUNT + 1))
        TOTAL_VARIANTS=$((TOTAL_VARIANTS + rows))
        SLICE_RESULTS+=("{\"slice_id\": \"$slice_id\", \"dataset\": \"$dataset\", \"region\": \"$region\", \"variants\": $rows, \"status\": \"MATCH\", \"time_ms\": $elapsed_ms}")
    else
        echo "MISMATCH: $diff_lines diff lines"
        diff -u <(sort "$out_j") <(sort "$out_r") | head -40
        exit 1
    fi

    # Advisory cpp third arm: vardictcpp ignores -i (no splice-table surface;
    # probed live T3A-04), so SKIP with reason, no engine run.
    echo "CPP: SKIP: $slice_id splicing (-i splice-table surface unsupported by vardictcpp)"
}

# -----------------------------------------------------------------------------
# 1. NIST GIAB HG005 (Han Chinese Diversity WGS 300x, GRCh38)
# -----------------------------------------------------------------------------
printf "chr3\t40000000\t40010000\n" > "$TMPDIR/hg005_1a.bed"
printf "chr6\t30000000\t30010000\n" > "$TMPDIR/hg005_1b.bed"
printf "chr21\t20000000\t20010000\n" > "$TMPDIR/hg005_1c.bed"

run_vcf_compare 1  "HG005_1A" "GIAB HG005 Diversity WGS" "$FA_GIAB" "$DATA_DIR/hg005_s4a.bam" "$TMPDIR/hg005_1a.bed" "HG005" -f 0.05
run_vcf_compare 2  "HG005_1B" "GIAB HG005 Diversity WGS" "$FA_GIAB" "$DATA_DIR/hg005_s4b.bam" "$TMPDIR/hg005_1b.bed" "HG005" -f 0.05
run_vcf_compare 3  "HG005_1C" "GIAB HG005 Diversity WGS" "$FA_GIAB" "$DATA_DIR/hg005_s4c.bam" "$TMPDIR/hg005_1c.bed" "HG005" -f 0.05

# -----------------------------------------------------------------------------
# 2. FDA SEQC2 Somatic WES (HCC1395 / HCC1395BL Breast Cancer, GRCh38)
# -----------------------------------------------------------------------------
printf "chr7\t140400000\t140420000\n" > "$TMPDIR/seqc2_2a.bed"
printf "chr11\t530000\t550000\n" > "$TMPDIR/seqc2_2b.bed"
printf "chr12\t25200000\t25250000\n" > "$TMPDIR/seqc2_2c.bed"

run_vcf_compare 4  "SEQC2_2A" "SEQC2 Somatic WES (BRAF)" "$FA_SEQC2" "$DATA_DIR/seqc2_s3a_tumor.bam|$DATA_DIR/seqc2_s3a_normal.bam" "$TMPDIR/seqc2_2a.bed" "TUMOR|NORMAL" -f 0.05 --fisher
run_vcf_compare 5  "SEQC2_2B" "SEQC2 Somatic WES (HRAS)" "$FA_SEQC2" "$DATA_DIR/seqc2_s3b_tumor.bam|$DATA_DIR/seqc2_s3b_normal.bam" "$TMPDIR/seqc2_2b.bed" "TUMOR|NORMAL" -f 0.05 --fisher
run_vcf_compare 6  "SEQC2_2C" "SEQC2 Somatic WES (KRAS)" "$FA_SEQC2" "$DATA_DIR/seqc2_s3c_tumor.bam|$DATA_DIR/seqc2_s3c_normal.bam" "$TMPDIR/seqc2_2c.bed" "TUMOR|NORMAL" -f 0.05 --fisher

# -----------------------------------------------------------------------------
# 3. NIST GIAB HG002 (Ashkenazi Germline WGS, GRCh38)
# -----------------------------------------------------------------------------
printf "chr20\t10000000\t10020000\n" > "$TMPDIR/hg002_3a.bed"
printf "chr20\t20000000\t20020000\n" > "$TMPDIR/hg002_3b.bed"
printf "chr20\t30000000\t30020000\n" > "$TMPDIR/hg002_3c.bed"

run_vcf_compare 7  "HG002_3A" "GIAB HG002 Germline WGS" "$FA_GIAB" "$ROOT_DIR/scratch/g4/hg002/HG002.chr20.bam" "$TMPDIR/hg002_3a.bed" "HG002"
run_vcf_compare 8  "HG002_3B" "GIAB HG002 Germline WGS" "$FA_GIAB" "$ROOT_DIR/scratch/g4/hg002/HG002.chr20.bam" "$TMPDIR/hg002_3b.bed" "HG002"
run_vcf_compare 9  "HG002_3C" "GIAB HG002 Germline WGS" "$FA_GIAB" "$ROOT_DIR/scratch/g4/hg002/HG002.chr20.bam" "$TMPDIR/hg002_3c.bed" "HG002"

# -----------------------------------------------------------------------------
# 4. NIST GIAB HG008 / HG008-T (Somatic Paired WGS, GRCh38)
# -----------------------------------------------------------------------------
HG008_BAMS="$ROOT_DIR/scratch/g4/hg008/HG008-T.slice.bam|$ROOT_DIR/scratch/g4/hg008/HG008-N.slice.bam"

run_vcf_compare 10 "HG008_4A" "GIAB HG008 Somatic WGS" "$FA_GIAB" "$HG008_BAMS" "$ROOT_DIR/scratch/g4/hg008/r1.bed" "HG008T|HG008N"
run_vcf_compare 11 "HG008_4B" "GIAB HG008 Somatic WGS" "$FA_GIAB" "$HG008_BAMS" "$ROOT_DIR/scratch/g4/hg008/r2.bed" "HG008T|HG008N"
run_vcf_compare 12 "HG008_4C" "GIAB HG008 Somatic WGS" "$FA_GIAB" "$HG008_BAMS" "$ROOT_DIR/scratch/g4/hg008/r3.bed" "HG008T|HG008N"

# -----------------------------------------------------------------------------
# 5. GIAB HG002 NA24385 Multi-Chromosome Deep WGS (GRCh37/b37)
# -----------------------------------------------------------------------------
printf "1\t246981051\t246981611\n" > "$TMPDIR/na24385_5a.bed"
printf "19\t58559581\t58567861\n" > "$TMPDIR/na24385_5b.bed"
printf "20\t36018823\t36039513\n" > "$TMPDIR/na24385_5c.bed"

BAM_NA24385_1="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_1_246981051_246981611.bam"
BAM_NA24385_19="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_19_58559581_58567861.bam"
BAM_NA24385_20="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_20_36018823_36039513.bam"

run_vcf_compare 13 "NA24385_5A" "NA24385 b37 WGS (chr1)" "$FA_B37" "$BAM_NA24385_1" "$TMPDIR/na24385_5a.bed" "NA24385" -f 0.001
run_vcf_compare 14 "NA24385_5B" "NA24385 b37 WGS (chr19)" "$FA_B37" "$BAM_NA24385_19" "$TMPDIR/na24385_5b.bed" "NA24385" -f 0.001
run_vcf_compare 15 "NA24385_5C" "NA24385 b37 WGS (chr20)" "$FA_B37" "$BAM_NA24385_20" "$TMPDIR/na24385_5c.bed" "NA24385" -f 0.001

# -----------------------------------------------------------------------------
# 6. RNA-Seq Transcriptome & Splicing Suite (hg38)
# -----------------------------------------------------------------------------
run_rna_slice 16 "RNA_6A" "RNA-Seq Splicing" "$FA_HG38" "$DATA_DIR/rna_s6a.bam" "chr8:90534611-90550000"
run_rna_slice 17 "RNA_6B" "RNA-Seq Splicing" "$FA_HG38" "$DATA_DIR/rna_s6b.bam" "chr8:90550000-90580000"
run_rna_slice 18 "RNA_6C" "RNA-Seq Splicing" "$FA_HG38" "$DATA_DIR/rna_s6c.bam" "chr8:90580000-90620070"

# -----------------------------------------------------------------------------
# Generate Machine-Readable Evidence JSON
# -----------------------------------------------------------------------------
python3 -c "
import json

data = {
    'suite': 'multimodal_multislice_parity_suite',
    'timestamp': '$(date -u +\"%Y-%m-%dT%H:%M:%SZ\")',
    'fingerprint': '$CURRENT_FP',
    'mode': 'tier-2-live-oracle-strict-bytes',
    'total_slices': $TOTAL_SLICES,
    'slices_passed': $PASSED_COUNT,
    'total_variants_evaluated': $TOTAL_VARIANTS,
    'float_tolerance_eps': 1e-6,
    $CPP_JSON_LINE
    'result': 'MATCH' if $PASSED_COUNT == $TOTAL_SLICES else 'MISMATCH',
    'slices': [$(IFS=,; echo "${SLICE_RESULTS[*]}")]
}

with open('$EVIDENCE_JSON', 'w') as f:
    json.dump(data, f, indent=2)
"

echo "================================================================================"
echo " Multi-Modal Parity Suite Result: $PASSED_COUNT/$TOTAL_SLICES MATCH (0 diffs)"
echo " Total Biological Variants Evaluated: $TOTAL_VARIANTS"
echo " Certified Evidence: $EVIDENCE_JSON"
echo "================================================================================"
