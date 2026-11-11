#!/usr/bin/env bash
set -euo pipefail
# Provision the VarDictJava 1.8.3 reference oracle via conda-forge + bioconda.
solver_arg=""
if conda create --help 2>&1 | grep -q -- '--solver'; then
  solver_arg="--solver=libmamba"
fi
if ! conda env list | grep -q vardict-oracle; then
  conda create -y -n vardict-oracle $solver_arg --override-channels -c conda-forge -c bioconda \
    vardict-java=1.8.3 vardict=2019.06.04 vardictcpp=3=hc52dbad_0 r-base=4.6.0 samtools bcftools
elif ! conda run -n vardict-oracle Rscript --version 2>&1 | grep -q '4\.6\.'; then
  conda install -y -n vardict-oracle $solver_arg --override-channels -c conda-forge -c bioconda \
    vardict-java=1.8.3 vardict=2019.06.04 vardictcpp=3=hc52dbad_0 r-base=4.6.0
fi
echo "oracle ready: conda run -n vardict-oracle vardict-java --version"

