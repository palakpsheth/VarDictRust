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
  curl -sSLf "${BAM_URL}.bai" -o "$LOCAL_BAI"
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

# 3. Intersect and pre-tile NIST v4.2.1 BED into <=50kb chunks via awk
LOCAL_BED="$OUTDIR/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
TMP_BED="$OUTDIR/.tmp.eval.bed"
if [ -f "$EVAL_BED" ] && [ -s "$EVAL_BED" ]; then
  echo "Eval BED already exists: $EVAL_BED"
else
  echo "Generating pre-tiled eval BED to $EVAL_BED..."
  awk_cmd='
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
  '
  if [ -f "$LOCAL_BED" ] && [ -s "$LOCAL_BED" ]; then
    awk "$awk_cmd" "$LOCAL_BED" > "$TMP_BED"
  else
    curl -sSLf "$TRUTH_BED_URL" | awk "$awk_cmd" > "$TMP_BED"
  fi
  num_intervals=$(wc -l < "$TMP_BED")
  [ "$num_intervals" -ge 1000 ] || { echo "FAIL: eval.bed interval count < 1000 (got $num_intervals)" >&2; exit 1; }
  mv "$TMP_BED" "$EVAL_BED"
  echo "Created $EVAL_BED with $num_intervals intervals"
fi

# 4. Fetch truth VCF for the region
LOCAL_VCF="$OUTDIR/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
if [ -f "$TRUTH_VCF" ] && [ -f "$TRUTH_VCF.tbi" ]; then
  echo "Truth VCF already exists: $TRUTH_VCF"
else
  echo "Extracting truth VCF for $REGION..."
  VCF_SOURCE="$TRUTH_VCF_URL"
  if [ -f "$LOCAL_VCF" ] && [ -f "$LOCAL_VCF.tbi" ]; then
    VCF_SOURCE="$LOCAL_VCF"
  fi
  bcftools view -r "$REGION" "$VCF_SOURCE" -Oz -o "$TRUTH_VCF"
  tabix -f -p vcf "$TRUTH_VCF"
  [ "$(bcftools view -H "$TRUTH_VCF" | head -n 1 | wc -l)" -gt 0 ] || { echo "FAIL: truth.vcf.gz is empty" >&2; exit 1; }
  echo "Created $TRUTH_VCF"
fi

echo "HG002 slice preparation complete."

