#!/usr/bin/env bash
# Usage: run_hg002.sh
# End-to-end first-class runner for NIST GIAB HG002 germline benchmark (chr20:10M-30M).
# Enforces:
#   1. Input asset and SDF preflight checks.
#   2. Dual-engine compare.sh parity (MATCH, 0 diffs).
#   3. NIST v4.2.1 truth concordance (EQUAL, dual thresholds).
set -euo pipefail

if ! command -v samtools &>/dev/null || ! command -v bcftools &>/dev/null || ! command -v tabix &>/dev/null || ! command -v rtg &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

DATA_DIR="scratch/g4"
OUTDIR="$DATA_DIR/hg002"
REF_FASTA="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
REF_SDF="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf"

# 1. Reference preflight
if [ ! -f "$REF_FASTA" ]; then
  if [ -f "$REF_FASTA.gz" ]; then
    echo "Decompressing reference genome $REF_FASTA.gz..."
    gzip -dc "$REF_FASTA.gz" > "$REF_FASTA"
    samtools faidx "$REF_FASTA"
  else
    echo "ERROR: Reference FASTA missing: $REF_FASTA" >&2
    echo "Run: bash harness/fetch_data.sh refs" >&2
    exit 1
  fi
fi

if [ ! -d "$REF_SDF" ]; then
  echo "Formatting RTG SDF for $REF_FASTA..."
  rtg format -o "$REF_SDF" "$REF_FASTA"
fi

# 2. HG002 data slice preflight
if [ ! -f "$OUTDIR/HG002.chr20.bam" ] || [ ! -f "$OUTDIR/HG002.chr20.bam.bai" ] || [ ! -f "$OUTDIR/eval.bed" ] || [ ! -f "$OUTDIR/truth/truth.vcf.gz" ]; then
  echo "Staging HG002 data slice via harness/get_hg002_slice.sh..."
  bash harness/get_hg002_slice.sh
fi

# 3. Dual-Engine Parity Check (Tier-2)
echo "=== Step 1: Running dual-engine parity on HG002 chr20:10M-30M ==="
EST_MIN=15 bash harness/courtesy.sh bash harness/compare.sh \
  "$REF_FASTA" \
  "$OUTDIR/HG002.chr20.bam" \
  "$OUTDIR/eval.bed" \
  HG002 -f 0.01 -th 8

echo "HG002 Tier-2 Parity: MATCH"

# 4. Truth Concordance Check (Dual-Threshold)
echo "=== Step 2: Running truth concordance against NIST v4.2.1 ==="
EST_MIN=10 bash harness/courtesy.sh bash harness/concordance_germline.sh \
  hg002-truth \
  "$OUTDIR/cmp/java.vcf" \
  "$OUTDIR/cmp/rust.vcf" \
  "$OUTDIR/truth/truth.vcf.gz" \
  "$OUTDIR/eval.bed" \
  "$REF_SDF" \
  "$OUTDIR/truth" \
  HG002

echo "HG002 Truth Concordance: EQUAL"
echo "HG002 first-class validation run completed successfully."
