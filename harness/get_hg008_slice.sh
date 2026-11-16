#!/usr/bin/env bash
# Stream HG008 paired tumor/normal WGS slices from NIST GIAB FTP, derive engine/eval BEDs,
# and encode CRAM twins (HG008-T.slice.cram / HG008-N.slice.cram) for cram_twin.sh.
# Usage: [DATA_DIR=scratch/g4] bash harness/get_hg008_slice.sh [--beds-only]
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/portable.sh"

if [ "${NO_G4_TOOLS:-0}" = 0 ] && [ -d "$G4_TOOLS_BIN" ] && { ! command -v samtools &>/dev/null || ! command -v bedtools &>/dev/null; }; then
  export PATH="$G4_TOOLS_BIN:$PATH"
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
TRUTH_VCF_SHA256="${HG008_TRUTH_VCF_SHA256:-6ff5b215e0dc634c85da6164cffbd5790233258fcc48c88a598331430b061cad}"
TRUTH_BED_SHA256="${HG008_TRUTH_BED_SHA256:-c6100933430322eb62c64e826a349c1861e26a029853502bd247f86d27e22726}"

truth_repaired=0
if [ ! -f "$TRUTH_VCF" ] || [ "$(portable_sha256 "$TRUTH_VCF")" != "$TRUTH_VCF_SHA256" ]; then
  rm -f "$TRUTH_VCF" "${TRUTH_VCF}.sha256" "${TRUTH_VCF}.legacy_md5" "${TRUTH_VCF}.done"
  truth_repaired=1
fi
if [ ! -f "$TRUTH_BED" ] || [ "$(portable_sha256 "$TRUTH_BED")" != "$TRUTH_BED_SHA256" ]; then
  rm -f "$TRUTH_BED" "${TRUTH_BED}.sha256" "${TRUTH_BED}.legacy_md5" "${TRUTH_BED}.done"
  truth_repaired=1
fi
if [ "$truth_repaired" -eq 1 ]; then
  DATA_DIR="$D" bash harness/fetch_data.sh hg008
  rm -f "$OUTDIR/slices.bed" "$OUTDIR/engine.bed" "$OUTDIR/eval.bed"
fi
[ "$(portable_sha256 "$TRUTH_VCF")" = "$TRUTH_VCF_SHA256" ] || {
  echo "FAIL: $TRUTH_VCF SHA-256 mismatch after fetch_data.sh hg008" >&2
  exit 1
}
[ "$(portable_sha256 "$TRUTH_BED")" = "$TRUTH_BED_SHA256" ] || {
  echo "FAIL: $TRUTH_BED SHA-256 mismatch after fetch_data.sh hg008" >&2
  exit 1
}

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

for idx in 1 2 3; do
  r_bed="$OUTDIR/r${idx}.bed"
  r_tmp="${r_bed}.tmp"
  sed -n "${idx}p" "$SLICES_BED" > "$r_tmp"
  if [ -f "$r_bed" ] && cmp -s "$r_tmp" "$r_bed"; then
    rm -f "$r_tmp"
  else
    mv -f "$r_tmp" "$r_bed"
  fi
done

# 3. Derive eval.bed (slices.bed intersected with NIST v0.3 benchmark all.bed)
EVAL_BED="$OUTDIR/eval.bed"
if [ -f "$EVAL_BED" ] && [ -s "$EVAL_BED" ]; then
  echo "HG008 eval.bed already exists: $EVAL_BED"
else
  echo "Intersecting $SLICES_BED with $TRUTH_BED -> $EVAL_BED..."
  if command -v bedtools >/dev/null 2>&1; then
    bedtools intersect -a "$SLICES_BED" -b "$TRUTH_BED" \
      | cut -f1-3 | bedtools sort -i - | bedtools merge -i - > "$EVAL_BED.tmp"
  else
    python3 - "$SLICES_BED" "$TRUTH_BED" > "$EVAL_BED.tmp" <<'PY'
import sys
def read_bed(p):
    out = {}
    with open(p) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            out.setdefault(parts[0], []).append((int(parts[1]), int(parts[2])))
    return out
a, b = read_bed(sys.argv[1]), read_bed(sys.argv[2])
for chrom in sorted(set(a) & set(b)):
    ints = sorted(
        (max(s1, s2), min(e1, e2))
        for s1, e1 in a[chrom]
        for s2, e2 in b[chrom]
        if max(s1, s2) < min(e1, e2)
    )
    merged = []
    for s, e in ints:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    for s, e in merged:
        print(f"{chrom}\t{s}\t{e}")
PY
  fi
  mv "$EVAL_BED.tmp" "$EVAL_BED"
fi

if [ "$BEDS_ONLY" -eq 1 ]; then
  echo "HG008 BED preparation complete (--beds-only)."
  exit 0
fi

# 4. Stream Tumor (195x) and Normal (169x) BAM slices from NIST GIAB FTP (verifies BAI via stage_verified_sha256)
BASE_URL="https://ftp.ncbi.nlm.nih.gov/giab/ftp/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313"
T_BAM_URL="$BASE_URL/HG008-T_Illumina_195x_GRCh38-GIABv3.bam"
N_BAM_URL="$BASE_URL/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam"
T_BAI_SHA256="1608f060248b3c88ddbf954fbfd3d0551a6bc4c27e7c4ae47f97722f7d828bb2"
N_BAI_SHA256="04f2be64201ece69523a7b9dd62a6f4d121c4f1759dfad9decaa5c04723afb2a"

stream_verified_bam_slice --bed "$SLICES_BED" --threads 4 \
  "Tumor (195x)" "$T_BAM_URL" "${T_BAM_URL}.bai" "$T_BAI_SHA256" \
  "$OUTDIR/HG008-T_Illumina_195x_GRCh38-GIABv3.bam.bai" "" "$OUTDIR/HG008-T.slice.bam"
stream_verified_bam_slice --bed "$SLICES_BED" --threads 4 \
  "Normal (169x)" "$N_BAM_URL" "${N_BAM_URL}.bai" "$N_BAI_SHA256" \
  "$OUTDIR/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam.bai" "" "$OUTDIR/HG008-N.slice.bam"

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
