#!/bin/bash
# Estimate heritability with BFMAP variance component analysis for a single trait
# Dependencies: BFMAP (https://github.com/jiang18/bfmap)
# Author: Junjian Wang
# Initiate date: 02/12/2026
# Current date: 03/31/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

BFMAP=~/bin/bfmap

TRAIT="MilkYield"                   # Trait name (must match column header in phenotype file)

PHENO_FILE="/path_to_the_file/${TRAIT}.csv"   # Phenotype file
GRM_FILE="/path_to_the_file/bfmap_grm"        # BFMAP GRM prefix (output from 02b)
OUT_PREFIX="/path_to_the_file/${TRAIT}"        # Output file prefix (BFMAP appends .varcomp.csv)

NUM_THREADS=10

# Optional: reliability/error weight column name in PHENO_FILE
# Set to "" to skip
ERROR_WEIGHT_NAME=""                 # e.g. "reliability" or "R-MilkYield"

# =============================================================================

echo "[$(date)] Estimating heritability for trait: ${TRAIT}"

ERROR_ARGS=""
if [[ -n "${ERROR_WEIGHT_NAME}" ]]; then
    ERROR_ARGS="--error_weight_name ${ERROR_WEIGHT_NAME}"
fi

"${BFMAP}" --varcomp \
    --phenotype       "${PHENO_FILE}" \
    --trait           "${TRAIT}" \
    --binary_grm_file "${GRM_FILE}" \
    --output          "${OUT_PREFIX}" \
    --num_threads     "${NUM_THREADS}" \
    ${ERROR_ARGS}

echo "[$(date)] Heritability estimation complete: ${OUT_PREFIX}.varcomp.csv"
