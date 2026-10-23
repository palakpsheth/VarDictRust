#!/usr/bin/env bash
# Stream HG008 paired tumor/normal WGS slices from NIST GIAB FTP, derive engine/eval BEDs,
# and encode CRAM twins (HG008-T.slice.cram / HG008-N.slice.cram) for cram_twin.sh.
# Usage: [DATA_DIR=scratch/g4] bash harness/get_hg008_slice.sh [--beds-only]
set -euo pipefail

if ! command -v samtools &>/dev/null || ! command -v bedtools &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

BEDS_ONLY=0
if [ "${1:-}" = "--beds-only" ]; then
  BEDS_ONLY=1
  shift
fi

D="${DATA_DIR:-scratch/g4}"
OUTDIR="$D/hg008"
REF_FA="$D/ref/giabv3/GRCh38_GIABv3.fa"
mkdir -p "$OUTDIR/truth" "$OUTDIR/cmp" "$OUTDIR/cram"

TRUTH_VCF="$OUTDIR/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"
TRUTH_BED="$OUTDIR/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed"

if [ ! -f "$TRUTH_VCF" ] || [ ! -f "$TRUTH_BED" ]; then
  DATA_DIR="$D" bash harness/fetch_data.sh hg008
fi

# 1. Derive slices.bed (200 dense somatic truth clusters) via select_slice_regions.py
SLICES_BED="$OUTDIR/slices.bed"
if [ -f "$SLICES_BED" ] && [ -s "$SLICES_BED" ]; then
  echo "HG008 slices.bed already exists: $SLICES_BED"
else
  echo "Generating $SLICES_BED via select_slice_regions.py..."
  python3 harness/select_slice_regions.py \
    --truth "$TRUTH_VCF" \
    --eval  "$TRUTH_BED" \
    --chain 20000 --min-count 3 --pad 500 --max-span 100000 --max-regions 200 \
    --out   "$SLICES_BED.tmp"
  mv "$SLICES_BED.tmp" "$SLICES_BED"
fi

# 2. Derive engine.bed (<=20kb tiles of slices.bed) and r1.bed..r3.bed (multimodal slices 10-12)
ENGINE_BED="$OUTDIR/engine.bed"
if [ -f "$ENGINE_BED" ] && [ -s "$ENGINE_BED" ]; then
  echo "HG008 engine.bed already exists: $ENGINE_BED"
else
  echo "Tiling $SLICES_BED into <=20kb windows -> $ENGINE_BED..."
  awk -F'\t' -v OFS='\t' '{
    s = $2 + 0; e = $3 + 0;
    while (e - s > 20000) {
      print $1, s, s + 20000, $4;
      s += 20000;
    }
    print $1, s, e, $4;
  }' "$SLICES_BED" > "$ENGINE_BED.tmp"
  mv "$ENGINE_BED.tmp" "$ENGINE_BED"
fi

sed -n '1p' "$SLICES_BED" > "$OUTDIR/r1.bed"
sed -n '2p' "$SLICES_BED" > "$OUTDIR/r2.bed"
sed -n '3p' "$SLICES_BED" > "$OUTDIR/r3.bed"

# 3. Derive eval.bed (slices.bed intersected with NIST v0.3 benchmark all.bed)
EVAL_BED="$OUTDIR/eval.bed"
if [ -f "$EVAL_BED" ] && [ -s "$EVAL_BED" ]; then
  echo "HG008 eval.bed already exists: $EVAL_BED"
else
  echo "Intersecting $SLICES_BED with $TRUTH_BED -> $EVAL_BED..."
  bedtools intersect -a "$SLICES_BED" -b "$TRUTH_BED" \
    | cut -f1-3 | bedtools sort -i - | bedtools merge -i - > "$EVAL_BED.tmp"
  mv "$EVAL_BED.tmp" "$EVAL_BED"
fi

if [ "$BEDS_ONLY" -eq 1 ]; then
  echo "HG008 BED preparation complete (--beds-only)."
  exit 0
fi

# 4. Stream Tumor (195x) and Normal (169x) BAM slices from NIST GIAB FTP
BASE_URL="https://ftp.ncbi.nlm.nih.gov/giab/ftp/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313"
T_BAM_URL="$BASE_URL/HG008-T_Illumina_195x_GRCh38-GIABv3.bam"
N_BAM_URL="$BASE_URL/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam"

stream_hg008_bam() {
  local label="$1" bam_url="$2" local_bai="$3" out_bam="$4"
  if [ -f "$out_bam" ] && [ -f "$out_bam.bai" ] && samtools quickcheck "$out_bam"; then
    echo "HG008 $label slice BAM already exists and passes quickcheck: $out_bam"
    return 0
  fi
  if [ ! -f "$local_bai" ] || [ ! -s "$local_bai" ]; then
    echo "Downloading $label BAI index from ${bam_url}.bai..."
    curl -sSLf "${bam_url}.bai" -o "$local_bai"
  fi
  local max_retries=5 delay=2 attempt=1 tmp_unsorted="${out_bam%.slice.bam}.unsorted.bam"
  while [ $attempt -le $max_retries ]; do
    echo "Streaming $label slices from $bam_url (attempt $attempt/$max_retries)..."
    rm -f "$tmp_unsorted" "$out_bam.tmp"
    if samtools view -@4 -b -M -L "$SLICES_BED" -X -o "$tmp_unsorted" "$bam_url" "$local_bai" \
       && samtools sort -@4 -o "$out_bam.tmp" "$tmp_unsorted" \
       && samtools quickcheck "$out_bam.tmp"; then
      rm -f "$tmp_unsorted"
      mv "$out_bam.tmp" "$out_bam"
      samtools index -@4 "$out_bam"
      return 0
    fi
    attempt=$((attempt + 1))
    echo "Sleeping $delay seconds before retry..." >&2
    sleep $delay
    delay=$((delay * 2))
  done
  echo "Failed to stream $label BAM slice after $max_retries attempts" >&2
  return 1
}

stream_hg008_bam "Tumor (195x)"  "$T_BAM_URL" "$OUTDIR/HG008-T_Illumina_195x_GRCh38-GIABv3.bam.bai"   "$OUTDIR/HG008-T.slice.bam"
stream_hg008_bam "Normal (169x)" "$N_BAM_URL" "$OUTDIR/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam.bai" "$OUTDIR/HG008-N.slice.bam"

# 5. Encode CRAM twins (store_nm=1, store_md=1) for cram_twin.sh
[ -f "$REF_FA" ] || { echo "Missing reference $REF_FA (run: bash harness/stage_g4.sh refs)" >&2; exit 1; }

encode_cram_twin() {
  local in_bam="$1" out_cram="$2"
  if [ -f "$out_cram" ] && [ -f "$out_cram.crai" ] && samtools quickcheck "$out_cram"; then
    echo "CRAM twin already exists and passes quickcheck: $out_cram"
    return 0
  fi
  echo "Encoding CRAM twin $out_cram from $in_bam..."
  samtools view -@4 -C -T "$REF_FA" \
    --output-fmt-option store_nm=1 \
    --output-fmt-option store_md=1 \
    -o "$out_cram.tmp" "$in_bam"
  mv "$out_cram.tmp" "$out_cram"
  samtools index -@4 "$out_cram"
  samtools quickcheck "$out_cram"
}

encode_cram_twin "$OUTDIR/HG008-T.slice.bam" "$OUTDIR/HG008-T.slice.cram"
encode_cram_twin "$OUTDIR/HG008-N.slice.bam" "$OUTDIR/HG008-N.slice.cram"

echo "HG008 slice and CRAM twin preparation complete."
