#!/usr/bin/env bash
# Edge-case parity test suite (D-51 / T05b)
# Tests 4 edge cases against both VarDictJava and VarDictRust:
# 1. empty.bam: valid SAM/BAM header, 0 alignments.
# 2. zero_cov: reads on chr1, bed region targeting a locus with zero coverage (chr2:1-100).
# 3. single_read.bam: exactly 1 alignment covering the target region.
# 4. contig_boundary.bam: reads at contig start (chr1:1-50) and contig end (chr1:251-300).
set -euo pipefail

DIR="${EDGE_CASES_DIR:-scratch/test_edge_cases}"
mkdir -p "$DIR"

REF="fixtures/micro/micro.ref.fa"
MICRO_BAM="fixtures/micro/micro.bam"
MICRO_BED="fixtures/micro/micro.bed"

if [ ! -f "$REF" ] || [ ! -f "$MICRO_BAM" ] || [ ! -f "$MICRO_BED" ]; then
  echo "FAIL: required micro fixture files not found" >&2
  exit 1
fi

cargo build -q --release -p vardictrust

fail=0

run_case() {
  local name="$1"
  local bam="$2"
  local bed="$3"
  local sample="$4"
  shift 4
  local out_dir="$DIR/out_${name}"
  rm -rf "$out_dir"
  mkdir -p "$out_dir"

  echo "=== Running edge case: $name ==="

  local out=""
  local rc=0
  out=$(KEEP_DIR="$out_dir" bash harness/compare.sh "$REF" "$bam" "$bed" "$sample" "$@" 2>&1) || rc=$?

  echo "$out"

  if [ "$rc" -ne 0 ]; then
    echo "FAIL: $name compare.sh exited with $rc" >&2
    fail=1
    return
  fi

  if ! grep -q "MATCH" <<<"$out"; then
    echo "FAIL: $name output does not contain MATCH" >&2
    fail=1
    return
  fi

  if grep -qi "panic" <<<"$out"; then
    echo "FAIL: $name detected panic in output" >&2
    fail=1
    return
  fi

  if [ ! -s "$out_dir/rust.vcf" ]; then
    echo "FAIL: $name rust.vcf is missing or empty" >&2
    fail=1
    return
  fi

  if ! grep -q "^##fileformat=VCF" "$out_dir/rust.vcf"; then
    echo "FAIL: $name rust.vcf missing ##fileformat header" >&2
    fail=1
    return
  fi

  if ! grep -q "^#CHROM" "$out_dir/rust.vcf"; then
    echo "FAIL: $name rust.vcf missing #CHROM line" >&2
    fail=1
    return
  fi

  echo "PASS: $name"
}

# --- 1. empty.bam: valid header, 0 alignments ---
echo "--- Preparing fixture: empty.bam ---"
samtools view -H "$MICRO_BAM" | samtools view -b -o "$DIR/empty.bam" -
samtools index "$DIR/empty.bam"
run_case "empty" "$DIR/empty.bam" "$MICRO_BED" "empty"

# --- 2. zero_cov: reads on chr1, query chr2:1-100 (zero coverage in window) ---
echo "--- Preparing fixture: zero_cov.bam ---"
samtools view -b "$MICRO_BAM" chr1 > "$DIR/zero_cov.bam"
samtools index "$DIR/zero_cov.bam"
printf "chr2\t0\t100\tzero_chr2\n" > "$DIR/zero_cov.bed"
run_case "zero_cov" "$DIR/zero_cov.bam" "$DIR/zero_cov.bed" "zero_cov" -R "chr2:1-100"

# --- 3. single_read.bam: exactly 1 alignment covering target region ---
echo "--- Preparing fixture: single_read.bam ---"
samtools view -H "$MICRO_BAM" > "$DIR/single_read.sam"
(set +o pipefail; samtools view "$MICRO_BAM" chr1 | head -n 1) >> "$DIR/single_read.sam"
samtools view -b "$DIR/single_read.sam" -o "$DIR/single_read.bam"
samtools index "$DIR/single_read.bam"
printf "chr1\t40\t100\tsingle_read\n" > "$DIR/single_read.bed"
run_case "single_read" "$DIR/single_read.bam" "$DIR/single_read.bed" "single_read"

# --- 4. contig_boundary.bam: reads at contig start (chr1:1-50) and contig end (chr1:251-300) ---
echo "--- Preparing fixture: contig_boundary.bam ---"
python3 -c '
with open("'"$REF"'") as f:
    lines = [l.strip() for l in f if l.strip() and not l.startswith(">")]
chr1 = "".join(lines[:5])
seq_start = chr1[0:50]
seq_end = chr1[250:300]
sam = f"""@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:300
@SQ\tSN:chr2\tLN:260
read_start\t0\tchr1\t1\t60\t50M\t*\t0\t0\t{seq_start}\t{"I"*50}
read_end\t0\tchr1\t251\t60\t50M\t*\t0\t0\t{seq_end}\t{"I"*50}
"""
with open("'"$DIR"'/contig_boundary.sam", "w") as f:
    f.write(sam)
'
samtools view -b "$DIR/contig_boundary.sam" | samtools sort -o "$DIR/contig_boundary.bam" -
samtools index "$DIR/contig_boundary.bam"
printf "chr1\t0\t50\tcontig_start\nchr1\t250\t300\tcontig_end\n" > "$DIR/contig_boundary.bed"
run_case "contig_boundary" "$DIR/contig_boundary.bam" "$DIR/contig_boundary.bed" "contig_boundary"

if [ "$fail" -eq 0 ]; then
  echo "PASS: test_edge_cases (all 4 edge cases matched)"
else
  echo "FAIL: test_edge_cases had failures" >&2
fi

exit "$fail"
