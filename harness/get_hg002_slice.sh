#!/usr/bin/env bash
# Stream HG002 chr20:10M-30M slice from AWS S3 mirror and prepare NIST v4.2.1 truth set.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/portable.sh"

if ! command -v bcftools &>/dev/null || ! command -v tabix &>/dev/null; then
  export PATH="$G4_TOOLS_BIN:$PATH"
fi

D="${DATA_DIR:-scratch/g4}"
OUTDIR="$D/hg002"
mkdir -p "$OUTDIR/truth" "$OUTDIR/cmp"

BAM_URL="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
BAM_BAI_SHA256="723f69de14b8adbbf1a60f23f9af007310011b51b13ee01d18a2a8daf7282f36"
TRUTH_VCF_URL="https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED_URL="https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

TARGET_BAM="$OUTDIR/HG002.chr20.bam"
TMP_BAM="$OUTDIR/.tmp.bam"
EVAL_BED="$OUTDIR/eval.bed"
TRUTH_VCF="$OUTDIR/truth/truth.vcf.gz"

REGION="chr20:10000000-30000000"

# 1. Check if BAM slice exists; if not, stage BAI and stream BAM
LOCAL_BAI="$OUTDIR/HG002.remote.bai"
if [ -s "$TARGET_BAM" ] && [ -s "${TARGET_BAM}.bai" ] && samtools quickcheck "$TARGET_BAM"; then
  echo "HG002 chr20 BAM slice already exists and passes quickcheck: $TARGET_BAM"
else
  stage_verified_sha256 "${BAM_URL}.bai" "$BAM_BAI_SHA256" "$LOCAL_BAI"
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
LOCAL_BED="$OUTDIR/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
TMP_BED="$OUTDIR/.tmp.eval.bed"
if [ -s "$EVAL_BED" ]; then
  echo "Eval BED already exists: $EVAL_BED"
else
  if [ ! -s "$LOCAL_BED" ]; then
    echo "Staging HG002 truth assets via fetch_data.sh..."
    DATA_DIR="$D" bash "$SCRIPT_DIR/fetch_data.sh" hg002
  fi
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
  awk "$awk_cmd" "$LOCAL_BED" > "$TMP_BED"
  num_intervals=$(wc -l < "$TMP_BED")
  [ "$num_intervals" -ge 1000 ] || { echo "FAIL: eval.bed interval count < 1000 (got $num_intervals)" >&2; exit 1; }
  mv "$TMP_BED" "$EVAL_BED"
  echo "Created $EVAL_BED with $num_intervals intervals"
fi

# 3. Fetch truth VCF for the region
LOCAL_VCF="$OUTDIR/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
if [ -s "$TRUTH_VCF" ] && [ -s "$TRUTH_VCF.tbi" ]; then
  echo "Truth VCF already exists: $TRUTH_VCF"
else
  if [ ! -s "$LOCAL_VCF" ] || [ ! -s "$LOCAL_VCF.tbi" ]; then
    echo "Staging HG002 truth assets via fetch_data.sh..."
    DATA_DIR="$D" bash "$SCRIPT_DIR/fetch_data.sh" hg002
  fi
  echo "Extracting truth VCF for $REGION..."
  bcftools view -r "$REGION" "$LOCAL_VCF" -Oz -o "$TRUTH_VCF"
  tabix -f -p vcf "$TRUTH_VCF"
  [ "$(bcftools view -H "$TRUTH_VCF" | head -n 1 | wc -l)" -gt 0 ] || { echo "FAIL: truth.vcf.gz is empty" >&2; exit 1; }
  echo "Created $TRUTH_VCF"
fi

echo "HG002 slice preparation complete."
