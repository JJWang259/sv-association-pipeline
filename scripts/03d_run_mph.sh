#!/bin/bash
# Polygenic variance partitioning with MPH MINQUE for a single trait
# Conditions on fine-mapped lead variants as fixed effects before
# partitioning residual genetic variance across functional annotation bins
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

TRAIT="MilkYield"

PHENO_FILE="/path_to_the_file/${TRAIT}.csv"       # Phenotype file
COV_FILE="/path_to_the_file/${TRAIT}.QTL.csv"     # Covariate file (output from 03b)
                                                   # Set to "" to run without covariates
                                                   # (e.g. if no QTL signals were detected)
GRM_LIST="/path_to_the_file/grms/grm_list.txt"    # GRM list (output from 03c)
                                                   # Each line is a GRM prefix; add or remove
                                                   # lines to include/exclude annotation partitions
OUT_PREFIX="/path_to_the_file/${TRAIT}.mph"        # Output file prefix

NUM_THREADS=20

# =============================================================================

mkdir -p "$(dirname "${OUT_PREFIX}")"
echo "[$(date)] Running MPH MINQUE for trait: ${TRAIT}"

# Build optional covariate argument
COV_ARGS=""
if [[ -n "${COV_FILE}" && -f "${COV_FILE}" ]]; then
    COV_ARGS="--covariate_file ${COV_FILE} --covariate_names all"
    echo "  Covariate file: ${COV_FILE}"
else
    echo "  No covariate file — running without fixed effects"
fi

/usr/bin/time -v "${MPH}" \
    --minque \
    --save_mem \
    --grm_list    "${GRM_LIST}" \
    --phenotype   "${PHENO_FILE}" \
    --trait       "${TRAIT}" \
    --num_threads "${NUM_THREADS}" \
    --out         "${OUT_PREFIX}" \
    ${COV_ARGS} \
    > "${OUT_PREFIX}.log" 2> "${OUT_PREFIX}.time"

echo "[$(date)] MPH MINQUE complete: ${OUT_PREFIX}.vc"