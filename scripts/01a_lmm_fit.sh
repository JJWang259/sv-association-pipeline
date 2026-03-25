#!/bin/bash
# Fit linear mixed model (LMM) with SLEMM for a single trait
# Dependencies: SLEMM (https://github.com/jiang18/slemm)
# Author: Junjian Wang
# Initiate date: 02/12/2026
# Current date: 03/24/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

SLEMM="~/bin/slemm"               # Path to SLEMM executable
TRAIT="MilkYield"                        # Trait name (must match column header in phenotype file)
BFILE="/path_to_the_file/geno"           # Genotype file (PLINK binary prefix, SNPs only, QC-filtered)
PHENO_FILE="/path_to_the_file/${TRAIT}.csv"  # Phenotype file: must contain columns IID and <trait>
SNP_INFO="/path_to_the_file/snp.info.csv"    # SNP info file for SLEMM variance component estimation
OUT_PREFIX="path_to_the_file/${TRAIT}"  # Output file prefix (SLEMM appends .qs, .log, etc.)

# Compute settings
NUM_THREADS=20
MIN_MAF=0.01
MIN_HWE_PVAL=1e-9
NUM_QF_MARKERS=90

# Optional: covariate file
# Set COVAR_FILE to "" to skip; intercept is always included by default
# If provided, list covariate column names in COVAR_NAMES (space-separated)
COVAR_FILE="/path_to_the_file/covariates.csv"
COVAR_NAMES="intercept"                  # e.g. "intercept batch sex"

# Optional: error weight
# Set to "" to skip
# Column name in PHENO_FILE containing reliability or inverse-variance weights
ERROR_WEIGHT_NAME=""                     # e.g. "reliability"

# =============================================================================

echo "[$(date)] Fitting LMM for trait: ${TRAIT}"

# Build optional arguments
COVAR_ARGS=""
if [[ -n "${COVAR_FILE}" && -f "${COVAR_FILE}" ]]; then
    COVAR_ARGS="--covariate_file ${COVAR_FILE} --covariate_names ${COVAR_NAMES}"
fi

ERROR_ARGS=""
if [[ -n "${ERROR_WEIGHT_NAME}" ]]; then
    ERROR_ARGS="--error_weight_name ${ERROR_WEIGHT_NAME}"
fi

"${SLEMM}" --lmm \
    --phenotype_file  "${PHENO_FILE}" \
    --bfile           "${BFILE}" \
    --trait           "${TRAIT}" \
    --snp_info_file   "${SNP_INFO}" \
    --num_threads     "${NUM_THREADS}" \
    --min_maf         "${MIN_MAF}" \
    --min_hwe_pval    "${MIN_HWE_PVAL}" \
    --num_qf_markers  "${NUM_QF_MARKERS}" \
    --out             "${OUT_PREFIX}" \
    ${COVAR_ARGS} \
    ${ERROR_ARGS}

echo "[$(date)] LMM fitting complete: ${OUT_PREFIX}"
