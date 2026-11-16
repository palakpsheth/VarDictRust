#!/usr/bin/env bash
# Stream HG002 chr20:10M-30M slice from AWS S3 mirror and prepare NIST v4.2.1 truth set.
# Usage: [DATA_DIR=scratch/g4] bash harness/get_hg002_slice.sh
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/portable.sh"

if ! command -v bcftools &>/dev/null || ! command -v tabix &>/dev/null || ! command -v samtools &>/dev/null; then
  ensure_g4_tools_path
fi

D="${DATA_DIR:-scratch/g4}"
OUTDIR="$D/hg002"
mkdir -p "$OUTDIR/truth" "$OUTDIR/cmp"

BAM_URL="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
HG002_BAI_SHA256="${HG002_BAI_SHA256:-723f69de14b8adbbf1a60f23f9af007310011b51b13ee01d18a2a8daf7282f36}"
HG002_BED_SHA256="${HG002_BED_SHA256:-fba9a57c36ec88e5d14ea3e259c8866c7935f4998d3ec0fa2d6c3962da5b5575}"
HG002_VCF_SHA256="${HG002_VCF_SHA256:-adb4d4a50048aa13353a06b84fcfcbca09a5d17525efaa4cea44f8822e81175c}"
HG002_TBI_SHA256="${HG002_TBI_SHA256:-23162712553413b10f88916c28bd1238deb703a7237be2910440773beea71894}"

TARGET_BAM="$OUTDIR/HG002.chr20.bam"
TMP_BAM="$OUTDIR/.tmp.bam"
EVAL_BED="$OUTDIR/eval.bed"
TRUTH_VCF="$OUTDIR/truth/truth.vcf.gz"
LOCAL_BAI="$OUTDIR/HG002.remote.bai"

REGION="chr20:10000000-30000000"

# 1. Stream BAM slice (check cached slice FIRST so offline runs never touch network or remote BAI)
if [ -s "$TARGET_BAM" ] && [ -s "${TARGET_BAM}.bai" ] && samtools quickcheck "$TARGET_BAM"; then
  echo "HG002 chr20 BAM slice already exists and passes quickcheck: $TARGET_BAM"
else
  echo "Staging SHA-256 verified remote BAI index from ${BAM_URL}.bai..."
  stage_verified_sha256 "${BAM_URL}.bai" "$HG002_BAI_SHA256" "$LOCAL_BAI"

  max_retries=5
  delay=2
  attempt=1
  success=0
  while [ $attempt -le $max_retries ]; do
    echo "Streaming $REGION from $BAM_URL using verified local index (attempt $attempt/$max_retries)..."
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
  echo "Generating pre-tiled eval BED to $EVAL_BED..."
  if [ ! -f "$LOCAL_BED" ] || [ "$(portable_sha256 "$LOCAL_BED")" != "$HG002_BED_SHA256" ]; then
    rm -f "$LOCAL_BED" "${LOCAL_BED}.sha256" "${LOCAL_BED}.legacy_md5" "${LOCAL_BED}.done"
    DATA_DIR="$D" bash harness/fetch_data.sh hg002
    [ "$(portable_sha256 "$LOCAL_BED")" = "$HG002_BED_SHA256" ] || {
      echo "FAIL: $LOCAL_BED SHA-256 mismatch after fetch_data.sh hg002" >&2
      exit 1
    }
  fi
  awk '
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
  ' "$LOCAL_BED" > "$TMP_BED"
  num_intervals=$(wc -l < "$TMP_BED" | tr -d '[:space:]')
  [ "$num_intervals" -ge 1000 ] || { echo "FAIL: eval.bed interval count < 1000 (got $num_intervals)" >&2; exit 1; }
  mv "$TMP_BED" "$EVAL_BED"
  echo "Created $EVAL_BED with $num_intervals intervals"
fi

# 3. Fetch truth VCF for the region from SHA-256 verified local benchmark VCF
LOCAL_VCF="$OUTDIR/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
if [ -s "$TRUTH_VCF" ] && [ -s "${TRUTH_VCF}.tbi" ]; then
  echo "Truth VCF already exists: $TRUTH_VCF"
else
  echo "Extracting truth VCF for $REGION..."
  if [ ! -f "$LOCAL_VCF" ] || [ "$(portable_sha256 "$LOCAL_VCF")" != "$HG002_VCF_SHA256" ]; then
    rm -f "$LOCAL_VCF" "${LOCAL_VCF}.sha256" "${LOCAL_VCF}.legacy_md5" "${LOCAL_VCF}.done"
  fi
  if [ ! -f "${LOCAL_VCF}.tbi" ] || [ "$(portable_sha256 "${LOCAL_VCF}.tbi")" != "$HG002_TBI_SHA256" ]; then
    rm -f "${LOCAL_VCF}.tbi" "${LOCAL_VCF}.tbi.sha256" "${LOCAL_VCF}.tbi.legacy_md5" "${LOCAL_VCF}.tbi.done"
  fi
  if [ ! -f "$LOCAL_VCF" ] || [ ! -f "${LOCAL_VCF}.tbi" ]; then
    DATA_DIR="$D" bash harness/fetch_data.sh hg002
  fi
  [ "$(portable_sha256 "$LOCAL_VCF")" = "$HG002_VCF_SHA256" ] || {
    echo "FAIL: $LOCAL_VCF SHA-256 mismatch after fetch_data.sh hg002" >&2
    exit 1
  }
  [ "$(portable_sha256 "${LOCAL_VCF}.tbi")" = "$HG002_TBI_SHA256" ] || {
    echo "FAIL: ${LOCAL_VCF}.tbi SHA-256 mismatch after fetch_data.sh hg002" >&2
    exit 1
  }
  bcftools view -r "$REGION" "$LOCAL_VCF" -Oz -o "${TRUTH_VCF}.tmp"
  mv "${TRUTH_VCF}.tmp" "$TRUTH_VCF"
  tabix -f -p vcf "$TRUTH_VCF"
  [ "$(bcftools view -H "$TRUTH_VCF" | head -n 1 | wc -l | tr -d '[:space:]')" -gt 0 ] || { echo "FAIL: truth.vcf.gz is empty" >&2; exit 1; }
  echo "Created $TRUTH_VCF"
fi

echo "HG002 slice preparation complete."
