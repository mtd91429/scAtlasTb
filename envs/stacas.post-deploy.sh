#!/usr/bin/env bash
# Install STACAS into the active `stacas` environment.
# STACAS is not available from conda-forge or bioconda, so the pinned GitHub release is
# installed from source (pure R, no compilation). Its dependencies come from conda
# (envs/stacas.yaml), so nothing else is installed or upgraded here.
# Snakemake runs this script automatically after creating the environment from
# envs/stacas.yaml (env_mode: from_yaml); envs/install_environment.sh runs it for
# locally installed environments.
set -euo pipefail

STACAS_VERSION=2.4.1

Rscript -e "remotes::install_github('carmonalab/STACAS@${STACAS_VERSION}', dependencies = FALSE, upgrade = 'never')"
Rscript -e "suppressPackageStartupMessages(library(STACAS)); stopifnot(packageVersion('STACAS') == '${STACAS_VERSION}'); cat('STACAS', as.character(packageVersion('STACAS')), 'installed\n')"
