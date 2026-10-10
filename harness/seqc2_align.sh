#!/usr/bin/env bash
# Usage: seqc2_align.sh <SM> <R1.fq.gz> <R2.fq.gz> <bwa-index-prefix> <out.bam>  — <= 8 threads total (A3).
# Deliberately simpler than SEQC2's trimmomatic+Picard pipeline: affects absolute recall, never parity.
set -euo pipefail
SM=$1; R1=$2; R2=$3; IDX=$4; OUT=$5; tmp="$OUT.tmp"; mkdir -p "$tmp"
bwa mem -M -t 6 -R "@RG\tID:$SM\tSM:$SM\tPL:ILLUMINA\tLB:$SM" "$IDX" "$R1" "$R2" \
  | samtools fixmate -m - - | samtools sort -@1 -m 2G -T "$tmp/s" -o "$tmp/sorted.bam" -
samtools markdup -@2 "$tmp/sorted.bam" "$OUT"; samtools index "$OUT"; rm -rf "$tmp"
