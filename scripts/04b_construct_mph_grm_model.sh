#!/bin/bash
# Construct GRM for MPH variant effect estimation
# Dependencies: MPH (https://github.com/jiang18/mph)
# Author: Junjian Wang
# Initiate date: 03/24/2026
# Current date: 03/24/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

MPH="$HOME/bin/mph"

GENO_FILE="/path_to_the_file/geno_model"   # PLINK binary prefix for model SNPs
SNP_INFO="/path_to_the_file/snp.info.csv"  # SNP info file
OUT_PREFIX="/path_to_the_file/mph_grm"     # Output GRM file prefix

NUM_THREADS=20

# =============================================================================

echo "[$(date)] Constructing MPH GRM"
mkdir -p "$(dirname "${OUT_PREFIX}")"

"${MPH}" --make_grm \
    --binary_genotype "${GENO_FILE}" \
    --snp_info        "${SNP_INFO}" \
    --num_threads     "${NUM_THREADS}" \
    --out             "${OUT_PREFIX}"

echo "[$(date)] GRM construction complete: ${OUT_PREFIX}"