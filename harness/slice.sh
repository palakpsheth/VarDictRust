#!/usr/bin/env bash
set -euo pipefail
# Usage: slice.sh <in.bam> <chr:start-end> <out.bam>   (requires samtools; subsets without loading whole file)
conda run -n vardict-oracle samtools view -b "$1" "$2" > "$3" && conda run -n vardict-oracle samtools index "$3"
