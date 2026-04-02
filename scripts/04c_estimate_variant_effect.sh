#!/bin/bash
# Estimate genotype-specific effects of a variant on a complex trait using MPH
# Uses AA and BB dummy variable covariates (AB as baseline) from 04a
# Dependencies: MPH (https://github.com/jiang18/mph)
# Author: Junjian Wang
# Initiate date: 03/24/2026
# Current date: 03/24/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

MPH=~/bin/mph

VARIANT="7_16507712_sv"                                   # Variant ID (must match 04a output)
TRAIT="MilkYield"                                         # Trait name

PHENO_FILE="/path_to_the_file/${TRAIT}.csv"               # Phenotype file
COV_FILE="/path_to_the_file/${VARIANT}.covariate.csv"     # Output from 04a
GRM_PREFIX="/path_to_the_file/mph_grm"                    # GRM prefix (output from 04b)
OUT_PREFIX="/path_to_the_file/${TRAIT}.${VARIANT}"        # Output file prefix

NUM_THREADS=20

# Optional: reliability/error weight column name in PHENO_FILE
# Set to "" to skip
ERROR_WEIGHT_NAME=""                                      # e.g. "reliability"

# =============================================================================

mkdir -p "$(dirname "${OUT_PREFIX}")"
echo "[$(date)] Estimating variant effect: ${VARIANT} | trait: ${TRAIT}"

# Generate temporary GRM list file and ensure cleanup on exit
GRM_LIST="$(dirname "${OUT_PREFIX}")/grm_list.temp"
echo "${GRM_PREFIX}" > "${GRM_LIST}"
trap 'rm -f "${GRM_LIST}"' EXIT

ERROR_ARGS=""
if [[ -n "${ERROR_WEIGHT_NAME}" ]]; then
    ERROR_ARGS="--error_weight ${ERROR_WEIGHT_NAME}"
fi

/usr/bin/time -v "${MPH}" \
    --minque \
    --save_mem \
    --grm_list        "${GRM_LIST}" \
    --phenotype       "${PHENO_FILE}" \
    --covariate_file  "${COV_FILE}" \
    --covariate_names all \
    --trait           "${TRAIT}" \
    --num_threads     "${NUM_THREADS}" \
    --out             "${OUT_PREFIX}" \
    ${ERROR_ARGS} \
    > "${OUT_PREFIX}.log" 2> "${OUT_PREFIX}.time"

echo "[$(date)] Complete: ${OUT_PREFIX}.mq.*"