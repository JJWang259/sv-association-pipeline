#!/bin/bash
# Construct genomic relationship matrix (GRM) with BFMAP
# Note: BFMAP GRM format is not compatible with MPH's GRM format
# Dependencies: BFMAP (https://github.com/jiang18/bfmap)
# Author: Junjian Wang
# Initiate date: 02/12/2026
# Current date: 03/24/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

BFMAP        =~/bin/bfmap                          # Path to BFMAP executable

GENO_FILE    ="/path_to_the_file/geno_model"          # PLINK binary prefix for model SNPs
SNP_INFO     "/path_to_the_file/snp.info.csv"  # SNP info file (shared with SLEMM and MPH)
OUT_PREFIX   "/path_to_the_file/bfmap_grm"         # Output GRM file prefix

GRM_TYPE     =2                                    # GRM type: 1 = centered, 2 = standardized
NUM_THREADS  =10

# =============================================================================

echo "[$(date)] Constructing BFMAP GRM"

"${BFMAP}" --compute_grm "${GRM_TYPE}" \
    --binary_genotype_file "${GENO_FILE}" \
    --snp_info_file        "${SNP_INFO}" \
    --output_file          "${OUT_PREFIX}" \
    --num_threads          "${NUM_THREADS}"

echo "[$(date)] GRM construction complete: ${OUT_PREFIX}"
