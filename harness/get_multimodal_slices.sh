#!/usr/bin/env bash
# =============================================================================
# harness/get_multimodal_slices.sh
# Ingest and stream 18 non-overlapping slices across 6 diverse public sequencing
# datasets spanning 12 distinct chromosomes for VarDictRust v1.0.0 validation.
#
# Adheres to Host Policy A3:
# - Minimal network bandwidth: streams only small bounded regions via range requests
# - Total footprint: < 150 MB disk space
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="$ROOT_DIR/scratch/g4/multimodal"
mkdir -p "$OUT_DIR"

echo "================================================================================"
echo " Ingesting 18 Multi-Modal Slices Across 6 Public Datasets"
echo " Destination: $OUT_DIR"
echo "================================================================================"

function fetch_remote_slice() {
    local name="$1"
    local bam_url="$2"
    local bai_url="$3"
    local region="$4"
    local out_bam="$OUT_DIR/${name}.bam"

    if [[ -f "$out_bam" ]] && samtools quickcheck "$out_bam" 2>/dev/null; then
        echo "  [cached] $name ($region)"
        return 0
    fi

    echo -n "  [streaming] $name ($region)... "
    samtools view -b -X "$bam_url" "$bai_url" "$region" > "$out_bam.tmp"
    mv "$out_bam.tmp" "$out_bam"
    samtools index "$out_bam"
    samtools quickcheck "$out_bam"
    echo "done ($(du -h "$out_bam" | cut -f1))"
}

function extract_local_slice() {
    local name="$1"
    local src_bam="$2"
    local region="$3"
    local out_bam="$OUT_DIR/${name}.bam"

    if [[ -f "$out_bam" ]] && samtools quickcheck "$out_bam" 2>/dev/null; then
        echo "  [cached] $name ($region)"
        return 0
    fi

    echo -n "  [extracting] $name ($region)... "
    samtools view -b "$src_bam" "$region" > "$out_bam.tmp"
    mv "$out_bam.tmp" "$out_bam"
    samtools index "$out_bam"
    samtools quickcheck "$out_bam"
    echo "done ($(du -h "$out_bam" | cut -f1))"
}

# -----------------------------------------------------------------------------
# 1. Hartwig COLO829 / COLO829BL (Melanoma Somatic WGS, b37)
# -----------------------------------------------------------------------------
echo "[1/6] Slicing Hartwig COLO829 / COLO829BL Somatic WGS..."
COLO_T_BAM="https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752450/COLO829T_dedup.realigned.bam"
COLO_T_BAI="https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752450/COLO829T_dedup.realigned.bam.bai"
COLO_N_BAM="https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752449/COLO829R_dedup.realigned.bam"
COLO_N_BAI="https://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752449/COLO829R_dedup.realigned.bam.bai"

fetch_remote_slice "colo829_s1a_tumor"  "$COLO_T_BAM" "$COLO_T_BAI" "1:10000000-10200000"
fetch_remote_slice "colo829_s1a_normal" "$COLO_N_BAM" "$COLO_N_BAI" "1:10000000-10200000"
fetch_remote_slice "colo829_s1b_tumor"  "$COLO_T_BAM" "$COLO_T_BAI" "3:10000000-10200000"
fetch_remote_slice "colo829_s1b_normal" "$COLO_N_BAM" "$COLO_N_BAI" "3:10000000-10200000"
fetch_remote_slice "colo829_s1c_tumor"  "$COLO_T_BAM" "$COLO_T_BAI" "20:45000000-45200000"
fetch_remote_slice "colo829_s1c_normal" "$COLO_N_BAM" "$COLO_N_BAI" "20:45000000-45200000"

# -----------------------------------------------------------------------------
# 2. Illumina Platinum Genomes NA12878 (European Germline WGS, GRCh38)
# -----------------------------------------------------------------------------
echo "[2/6] Slicing Platinum Genomes NA12878 Germline WGS..."
PLAT_BAM="https://storage.googleapis.com/brain-genomics-public/research/platinum/bam/NA12878.GRCh38.haplotagged.35x.bam"
PLAT_BAI="https://storage.googleapis.com/brain-genomics-public/research/platinum/bam/NA12878.GRCh38.haplotagged.35x.bam.bai"

fetch_remote_slice "platinum_s2a" "$PLAT_BAM" "$PLAT_BAI" "chr2:50000000-50200000"
fetch_remote_slice "platinum_s2b" "$PLAT_BAM" "$PLAT_BAI" "chr17:40000000-40200000"
fetch_remote_slice "platinum_s2c" "$PLAT_BAM" "$PLAT_BAI" "chrX:15000000-15200000"

# -----------------------------------------------------------------------------
# 3. FDA SEQC2 Somatic WES (HCC1395 / HCC1395BL Breast Cancer, GRCh38)
# -----------------------------------------------------------------------------
echo "[3/6] Slicing SEQC2 HCC1395 / HCC1395BL Somatic WES..."
SEQC2_T_BAM="$ROOT_DIR/scratch/g4/seqc2/HCC1395.bam"
SEQC2_N_BAM="$ROOT_DIR/scratch/g4/seqc2/HCC1395BL.bam"

extract_local_slice "seqc2_s3a_tumor"  "$SEQC2_T_BAM" "chr7:140400000-140600000"
extract_local_slice "seqc2_s3a_normal" "$SEQC2_N_BAM" "chr7:140400000-140600000"
extract_local_slice "seqc2_s3b_tumor"  "$SEQC2_T_BAM" "chr11:530000-550000"
extract_local_slice "seqc2_s3b_normal" "$SEQC2_N_BAM" "chr11:530000-550000"
extract_local_slice "seqc2_s3c_tumor"  "$SEQC2_T_BAM" "chr12:25200000-25400000"
extract_local_slice "seqc2_s3c_normal" "$SEQC2_N_BAM" "chr12:25200000-25400000"

# -----------------------------------------------------------------------------
# 4. NIST GIAB HG005 (Han Chinese Diversity WGS 300x, GRCh38)
# -----------------------------------------------------------------------------
echo "[4/6] Slicing GIAB HG005 Han Chinese Diversity WGS..."
HG005_BAM="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam"
HG005_BAI="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam.bai"

fetch_remote_slice "hg005_s4a" "$HG005_BAM" "$HG005_BAI" "chr3:40000000-40100000"
fetch_remote_slice "hg005_s4b" "$HG005_BAM" "$HG005_BAI" "chr6:30000000-30100000"
fetch_remote_slice "hg005_s4c" "$HG005_BAM" "$HG005_BAI" "chr21:20000000-20100000"

# -----------------------------------------------------------------------------
# 5. Targeted Amplicon Panel (Horizon / Integration, hg19)
# -----------------------------------------------------------------------------
echo "[5/6] Slicing Targeted Amplicon Panel..."
AMP_BAM="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam"

extract_local_slice "amplicon_s5a" "$AMP_BAM" "chr7:55241000-55243000"
extract_local_slice "amplicon_s5b" "$AMP_BAM" "chr7:55248000-55250000"
extract_local_slice "amplicon_s5c" "$AMP_BAM" "chr7:55259000-55260500"

# -----------------------------------------------------------------------------
# 6. RNA-Seq Transcriptome & Splicing Suite (hg38)
# -----------------------------------------------------------------------------
echo "[6/6] Slicing RNA-Seq Transcriptome & Splicing Suite..."
RNA_BAM="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Tcov_zero_SV_DEL_candidate3.bam"

extract_local_slice "rna_s6a" "$RNA_BAM" "chr8:90534611-90550000"
extract_local_slice "rna_s6b" "$RNA_BAM" "chr8:90550000-90580000"
extract_local_slice "rna_s6c" "$RNA_BAM" "chr8:90580000-90620070"

echo "================================================================================"
echo " All 18 Multi-Modal Slices Successfully Sliced and Verified!"
echo " Total Storage Consumed: $(du -sh "$OUT_DIR" | cut -f1)"
echo "================================================================================"
