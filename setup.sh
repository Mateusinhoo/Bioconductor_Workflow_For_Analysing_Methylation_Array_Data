# !/usr/bin/env bash
set -euo pipefail

# Ensure strict channel priority for reproducibility
conda config --set channel_priority strict >/dev/null 2>&1 || true
mamba -h >/dev/null 2>&1 || conda install -y -c conda-forge mamba

# Create env from binaries 
ENV_NAME="methylation-env"
[ "$(conda env list | awk '{print $1}' | grep -cx "$ENV_NAME" || true)" -gt 0 ] \
  && conda env remove -n "$ENV_NAME" -y

mamba env create -f environment.yml

# Activate and top up light R pkgs
conda activate "$ENV_NAME"
Rscript scripts/install_bioc.R

# Test
R -q -e 'source("scripts/check_env.R")'
echo " Environment ready: $ENV_NAME"
