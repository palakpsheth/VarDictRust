#!/usr/bin/env bash
# Stream HG002 chr20:10M-30M slice from AWS S3 mirror and prepare NIST v4.2.1 truth set.
set -euo pipefail

if ! command -v bcftools &>/dev/null || ! command -v tabix &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

OUTDIR="scratch/g4/hg002"
mkdir -p "$OUTDIR/truth" "$OUTDIR/cmp"

BAM_URL="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
TRUTH_VCF_URL="https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED_URL="https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

TARGET_BAM="$OUTDIR/HG002.chr20.bam"
TMP_BAM="$OUTDIR/.tmp.bam"
EVAL_BED="$OUTDIR/eval.bed"
TRUTH_VCF="$OUTDIR/truth/truth.vcf.gz"

REGION="chr20:10000000-30000000"

# 1. Download BAM index (.bai) if needed
LOCAL_BAI="$OUTDIR/HG002.remote.bai"
if [ ! -f "$LOCAL_BAI" ] || [ ! -s "$LOCAL_BAI" ]; then
  echo "Downloading remote BAI index (9.4MB) from ${BAM_URL}.bai..."
  curl -sSL "${BAM_URL}.bai" -o "$LOCAL_BAI"
fi

# 2. Stream BAM slice with 5-retry exponential backoff
if [ -f "$TARGET_BAM" ] && [ -f "$TARGET_BAM.bai" ] && samtools quickcheck "$TARGET_BAM"; then
  echo "HG002 chr20 BAM slice already exists and passes quickcheck: $TARGET_BAM"
else
  max_retries=5
  delay=2
  attempt=1
  success=0
  while [ $attempt -le $max_retries ]; do
    echo "Streaming $REGION from $BAM_URL using local index (attempt $attempt/$max_retries)..."
    rm -f "$TMP_BAM"
    if samtools view -@ 4 -b -X "$BAM_URL" "$LOCAL_BAI" "$REGION" > "$TMP_BAM"; then
      if samtools quickcheck "$TMP_BAM"; then
        mv "$TMP_BAM" "$TARGET_BAM"
        echo "Indexing $TARGET_BAM..."
        samtools index -@ 4 "$TARGET_BAM"
        success=1
        break
      else
        echo "samtools quickcheck failed on attempt $attempt" >&2
      fi
    else
      echo "samtools view failed on attempt $attempt" >&2
    fi
    attempt=$((attempt + 1))
    echo "Sleeping $delay seconds before retry..."
    sleep $delay
    delay=$((delay * 2))
  done
  if [ $success -ne 1 ]; then
    echo "Failed to stream BAM after $max_retries attempts" >&2
    exit 1
  fi
fi

# 2. Intersect and pre-tile NIST v4.2.1 BED into <=50kb chunks via awk
if [ -f "$EVAL_BED" ] && [ -s "$EVAL_BED" ]; then
  echo "Eval BED already exists: $EVAL_BED"
else
  echo "Fetching NIST v4.2.1 BED and pre-tiling to $EVAL_BED..."
  curl -sSL "$TRUTH_BED_URL" | awk '
    $1 == "chr20" {
      s = ($2 < 10000000 ? 10000000 : $2);
      e = ($3 > 30000000 ? 30000000 : $3);
      if (s < e) {
        step = 50000;
        for (cur = s; cur < e; cur += step) {
          cur_end = (cur + step < e ? cur + step : e);
          idx++;
          printf "%s\t%d\t%d\thg002_%d\n", $1, cur, cur_end, idx;
        }
      }
    }
  ' > "$EVAL_BED"
  echo "Created $EVAL_BED with $(wc -l < "$EVAL_BED") intervals"
fi

# 3. Fetch truth VCF for the region
if [ -f "$TRUTH_VCF" ] && [ -f "$TRUTH_VCF.tbi" ]; then
  echo "Truth VCF already exists: $TRUTH_VCF"
else
  echo "Extracting truth VCF for $REGION..."
  bcftools view -r "$REGION" "$TRUTH_VCF_URL" -Oz -o "$TRUTH_VCF"
  tabix -f -p vcf "$TRUTH_VCF"
  echo "Created $TRUTH_VCF"
fi

echo "HG002 slice preparation complete."
