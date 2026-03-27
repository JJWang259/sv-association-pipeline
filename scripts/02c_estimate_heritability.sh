#!/bin/bash
# Estimate heritability with BFMAP variance component analysis for a single trait
# Dependencies: BFMAP (https://github.com/jiang18/bfmap)
# Author: Junjian Wang
# Initiate date: 02/12/2026
# Current date: 03/24/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

BFMAP=~/bin/bfmap                                 # Path to BFMAP executable

TRAIT="MilkYield"                                 # Trait name (must match column header in phenotype file)

PHENO_FILE="/path_to_the_file/${TRAIT}.csv"       # Phenotype file
GRM_FILE="/path_to_the_file/bfmap_grm"            # BFMAP GRM prefix (output from 02b)
OUT_PREFIX="/path_to_the_file/${TRAIT}"           # Output file prefix (BFMAP appends .varcomp.csv etc.)

NUM_THREADS=10

# =============================================================================

echo "[$(date)] Estimating heritability for trait: ${TRAIT}"

"${BFMAP}" --varcomp \
    --phenotype       "${PHENO_FILE}" \
    --trait           "${TRAIT}" \
    --binary_grm_file "${GRM_FILE}" \
    --output          "${OUT_PREFIX}" \
    --num_threads     "${NUM_THREADS}"

echo "[$(date)] Heritability estimation complete: ${OUT_PREFIX}.varcomp.csv"
