#!/usr/bin/env bash
set -euo pipefail
# Provision the VarDictJava 1.8.3 reference oracle via bioconda.
if ! conda env list | grep -q vardict-oracle; then
  solver_arg=""
  if conda create --help 2>&1 | grep -q -- '--solver'; then
    solver_arg="--solver=libmamba"
  fi
  conda create -y -n vardict-oracle $solver_arg -c bioconda -c conda-forge \
    vardict-java=1.8.3 vardict=2019.06.04 samtools bcftools
fi
echo "oracle ready: conda run -n vardict-oracle vardict-java --version"

