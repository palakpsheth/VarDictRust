#!/usr/bin/env bash
# =============================================================================
# harness/stage_g4.sh
# Turnkey, idempotent data & reference staging orchestrator for VarDictRust.
#
# Bridges `harness/fetch_data.sh` (raw MD5-verified downloads) to ready-to-run
# `scratch/g4/` inputs (unpacked FASTAs, .fai indices, RTG SDF directories,
# BWA-aligned SEQC2 BAMs, derived engine/eval BEDs, remote BAM slices, and
# CRAM twins).
#
# Usage:
#   bash harness/stage_g4.sh [all|refs|amplicon|hg002|hg008|seqc2|seqc2-beds|multimodal]...
# =============================================================================
set -euo pipefail

if ! command -v samtools &>/dev/null || ! command -v bedtools &>/dev/null || ! command -v bcftools &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

D="${DATA_DIR:-scratch/g4}"

usage() {
  cat <<'EOF'
Usage: [DATA_DIR=scratch/g4] bash harness/stage_g4.sh <target>...

Targets:
  all         Stage all references and biological cohorts (refs, amplicon, hg002, hg008, seqc2, multimodal)
  refs        Fetch, unpack, index (.fai), and format RTG SDF (.sdf) for GRCh38.d1.vd1, GRCh38_GIABv3, and testdata FASTAs
  amplicon    Stage overlapping 8-column amplicon panel fixture (scratch/g4/amplicon/)
  hg002       Stream GIAB HG002 chr20:10M-30M BAM slice and NIST v4.2.1 truth (scratch/g4/hg002/)
  hg008       Stream GIAB HG008 paired PDAC BAM slices, derive BEDs, and encode CRAM twins (scratch/g4/hg008/)
  seqc2       Fetch SEQC2 WES FASTQs + truth, align HCC1395/HCC1395BL BAMs, and derive exome/slice/eval BEDs
  seqc2-beds  Derive SEQC2 exome.engine.bed, eval.bed, slices.bed, slices.engine.bed, and truth.vcf.gz from fetched truth/BED6
  multimodal  Stream 18 multi-modal slices across 6 cohorts (scratch/g4/multimodal/)
EOF
}

if [ $# -eq 0 ] || [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
  usage
  [ $# -eq 0 ] && exit 1 || exit 0
fi

stage_refs() {
  echo "=== [stage_g4: refs] ==="
  python3 harness/convert_fastas.py
  DATA_DIR="$D" bash harness/fetch_data.sh refs

  # 1. GRCh38.d1.vd1 (SEQC2 reference + BWA index + RTG SDF)
  local d1="$D/ref/d1vd1"
  if [ -f "$d1/GRCh38.d1.vd1.fa" ] && gzip -t "$d1/GRCh38.d1.vd1.fa" 2>/dev/null; then
    echo "Unpacking $d1/GRCh38.d1.vd1.fa tarball..."
    mv "$d1/GRCh38.d1.vd1.fa" "$d1/GRCh38.d1.vd1.fa.tar.gz"
    tar -xzf "$d1/GRCh38.d1.vd1.fa.tar.gz" -C "$d1"
    rm -f "$d1/GRCh38.d1.vd1.fa.tar.gz"
  fi
  if [ -f "$d1/GRCh38.d1.vd1_BWA.tar.gz" ] && [ ! -f "$d1/GRCh38.d1.vd1.fa.bwt" ]; then
    echo "Unpacking $d1/GRCh38.d1.vd1_BWA.tar.gz..."
    tar -xzf "$d1/GRCh38.d1.vd1_BWA.tar.gz" -C "$d1"
  fi
  if [ -f "$d1/GRCh38.d1.vd1.fa" ] && [ ! -f "$d1/GRCh38.d1.vd1.fa.fai" ]; then
    samtools faidx "$d1/GRCh38.d1.vd1.fa"
  fi
  if [ -f "$d1/GRCh38.d1.vd1.fa" ] && [ ! -d "$d1/GRCh38.d1.vd1.sdf" ]; then
    rtg format -o "$d1/GRCh38.d1.vd1.sdf" "$d1/GRCh38.d1.vd1.fa"
    ln -sfn GRCh38.d1.vd1.sdf "$d1/GRCh38.d1.vd1.fa.sdf"
  fi

  # 2. GRCh38_GIABv3 (HG002 & HG008 reference + RTG SDF)
  local gv="$D/ref/giabv3"
  if [ -f "$gv/GRCh38_GIABv3.fasta.gz" ] && [ ! -f "$gv/GRCh38_GIABv3.fasta" ]; then
    echo "Decompressing $gv/GRCh38_GIABv3.fasta.gz..."
    gzip -dc "$gv/GRCh38_GIABv3.fasta.gz" > "$gv/GRCh38_GIABv3.fasta.tmp"
    mv "$gv/GRCh38_GIABv3.fasta.tmp" "$gv/GRCh38_GIABv3.fasta"
  fi
  ln -sfn GRCh38_GIABv3.fasta "$gv/GRCh38_GIABv3.fa"
  ln -sfn GRCh38_GIABv3.fasta.fai "$gv/GRCh38_GIABv3.fa.fai"
  if [ -f "$gv/GRCh38_GIABv3.fasta" ] && [ ! -d "$gv/GRCh38_GIABv3.sdf" ]; then
    rtg format -o "$gv/GRCh38_GIABv3.sdf" "$gv/GRCh38_GIABv3.fasta"
    ln -sfn GRCh38_GIABv3.sdf "$gv/GRCh38_GIABv3.fa.sdf"
  fi
}

stage_seqc2_beds() {
  echo "=== [stage_g4: seqc2-beds] ==="
  local sq="$D/seqc2"
  mkdir -p "$sq/truth" "$sq/cmp-exome" "$sq/cmp-slices"
  if [ ! -f "$sq/exome.bed6.gz" ] || [ ! -f "$sq/truth/HC.bed" ]; then
    DATA_DIR="$D" bash harness/fetch_data.sh seqc2
  fi

  # 1. exome.engine.bed: 4-col merged primary-chrom intervals from SureSelect V6 exome.bed6.gz
  if [ ! -s "$sq/exome.engine.bed" ]; then
    echo "Deriving $sq/exome.engine.bed..."
    zcat "$sq/exome.bed6.gz" \
      | awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y|M)$/' \
      | cut -f1-4 \
      | sort -k1,1 -k2,2n -k3,3n \
      | bedtools merge -i - -c 4 -o first > "$sq/exome.engine.bed.tmp"
    mv "$sq/exome.engine.bed.tmp" "$sq/exome.engine.bed"
  fi

  # 2. eval.bed: High-confidence somatic truth regions intersected with exome targets
  if [ ! -s "$sq/eval.bed" ]; then
    echo "Deriving $sq/eval.bed..."
    bedtools intersect -a "$sq/truth/HC.bed" -b "$sq/exome.engine.bed" \
      | cut -f1-3 | bedtools sort -i - | bedtools merge -i - > "$sq/eval.bed.tmp"
    mv "$sq/eval.bed.tmp" "$sq/eval.bed"
  fi

  # 3. slices.bed & slices.engine.bed: 7 dense somatic truth clusters across chr7/9/11/12/17/22
  if [ ! -s "$sq/slices.bed" ]; then
    echo "Deriving $sq/slices.bed via select_slice_regions.py..."
    python3 harness/select_slice_regions.py \
      --truth "$sq/truth/sSNV.vcf.gz" \
      --truth "$sq/truth/sINDEL.vcf.gz" \
      --eval "$sq/eval.bed" \
      --exclude harness/data/grch38_exclude_ig_hla.bed \
      --chain 250000 --min-count 5 --pad 1000 --max-span 600000 --max-regions 20 \
      --out "$sq/slices.bed.tmp"
    mv "$sq/slices.bed.tmp" "$sq/slices.bed"
  fi
  if [ ! -s "$sq/slices.engine.bed" ]; then
    echo "Deriving $sq/slices.engine.bed..."
    bedtools intersect -a "$sq/exome.engine.bed" -b "$sq/slices.bed" \
      | cut -f1-4 > "$sq/slices.engine.bed.tmp"
    mv "$sq/slices.engine.bed.tmp" "$sq/slices.engine.bed"
  fi

  # 4. Combined somatic truth VCF (sSNV + sINDEL)
  if [ ! -s "$sq/truth/truth.vcf.gz" ] || [ ! -s "$sq/truth/truth.vcf.gz.tbi" ]; then
    echo "Building combined $sq/truth/truth.vcf.gz..."
    [ -f "$sq/truth/sSNV.vcf.gz.tbi" ]   || tabix -f -p vcf "$sq/truth/sSNV.vcf.gz"
    [ -f "$sq/truth/sINDEL.vcf.gz.tbi" ] || tabix -f -p vcf "$sq/truth/sINDEL.vcf.gz"
    bcftools concat -a -Oz -o "$sq/truth/truth.vcf.gz.tmp" \
      "$sq/truth/sSNV.vcf.gz" "$sq/truth/sINDEL.vcf.gz"
    mv "$sq/truth/truth.vcf.gz.tmp" "$sq/truth/truth.vcf.gz"
    tabix -f -p vcf "$sq/truth/truth.vcf.gz"
  fi
}

stage_seqc2() {
  stage_seqc2_beds
  local sq="$D/seqc2"
  local idx="$D/ref/d1vd1/GRCh38.d1.vd1.fa"
  if [ ! -f "$sq/HCC1395.bam" ] || ! samtools quickcheck "$sq/HCC1395.bam" 2>/dev/null; then
    echo "Aligning SEQC2 tumor BAM ($sq/HCC1395.bam)..."
    bash harness/seqc2_align.sh HCC1395 "$sq/fastq/T_1.fq.gz" "$sq/fastq/T_2.fq.gz" "$idx" "$sq/HCC1395.bam"
  fi
  if [ ! -f "$sq/HCC1395BL.bam" ] || ! samtools quickcheck "$sq/HCC1395BL.bam" 2>/dev/null; then
    echo "Aligning SEQC2 normal BAM ($sq/HCC1395BL.bam)..."
    bash harness/seqc2_align.sh HCC1395BL "$sq/fastq/N_1.fq.gz" "$sq/fastq/N_2.fq.gz" "$idx" "$sq/HCC1395BL.bam"
  fi
}

for target in "$@"; do
  case "$target" in
    all)
      stage_refs
      bash harness/get_amplicon_fixture.sh
      bash harness/get_hg002_slice.sh
      DATA_DIR="$D" bash harness/get_hg008_slice.sh
      stage_seqc2
      bash harness/get_multimodal_slices.sh
      ;;
    refs)       stage_refs ;;
    amplicon)   bash harness/get_amplicon_fixture.sh ;;
    hg002)      bash harness/get_hg002_slice.sh ;;
    hg008)      DATA_DIR="$D" bash harness/get_hg008_slice.sh ;;
    seqc2-beds) stage_seqc2_beds ;;
    seqc2)      stage_seqc2 ;;
    multimodal) bash harness/get_multimodal_slices.sh ;;
    *)
      echo "Unknown target: $target" >&2
      usage >&2
      exit 1
      ;;
  esac
done

echo "=== stage_g4 complete ($*) ==="
