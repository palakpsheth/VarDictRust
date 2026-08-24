#!/usr/bin/env bash
set -euo pipefail
# Provision the VarDictJava 1.8.3 reference oracle via bioconda.
conda env list | grep -q vardict-oracle || \
  conda create -y -n vardict-oracle -c bioconda -c conda-forge \
    vardict-java=1.8.3 vardict=2019.06.04 samtools bcftools
echo "oracle ready: conda run -n vardict-oracle vardict-java --version"
