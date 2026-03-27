#!/bin/bash
# Run genome-wide association analysis with SLEMM for a single trait
# Dependencies: SLEMM (https://github.com/jiang18/slemm)
# Author: Junjian Wang
# Initiate date: 02/12/2026
# Current date: 03/24/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

SLEMM_GWA=~/bin/slemm_gwa                # Path to slemm_gwa executable

TRAIT="MilkYield"                        # Trait name (must match 01a_lmm_fit.sh)

# Genotype file prefix for association testing
# For per-chromosome files, use {CHR} as a placeholder for the chromosome number
# Example: "/path/chr{CHR}.imputed" → /path/chr1.imputed, /path/chr2.imputed, ...
# For a single merged file, simply provide the path with no {CHR} placeholder
GENO_TEST_TEMPLATE="/path_to_the_file/chr{CHR}.imputed"

LMM_PREFIX="/path_to_the_file/${TRAIT}"  # Output prefix from 01a_lmm_fit.sh (without extension)
NUM_CHROMOSOMES=29                        # Number of autosomes (cattle=29, pig=18)
export OMP_NUM_THREADS=20

GWAS_OUT="/path_to_the_file/${TRAIT}_GWAS_All.txt" # Concatenated GWAS summary file

# =============================================================================

OUT_DIR="$(dirname "${GWAS_OUT}")"
mkdir -p "${OUT_DIR}"

echo "[$(date)] Running GWAS for trait: ${TRAIT}"

# Run per-chromosome association scan
for i in $(seq 1 "${NUM_CHROMOSOMES}"); do
    GENO_TEST="${GENO_TEST_TEMPLATE/\{CHR\}/${i}}"
    OUT_FILE="${OUT_DIR}/${TRAIT}_gwa.chr${i}.txt"

    echo "  [$(date)] Chromosome ${i} / ${NUM_CHROMOSOMES}"
    echo "    Using pfile: ${GENO_TEST}"

    "${SLEMM_GWA}" \
        --pfile  "${GENO_TEST}" \
        --slemm  "${LMM_PREFIX}" \
        --out    "${OUT_FILE}" \
        --chr    "${i}"
done

# Concatenate per-chromosome files into a single GWAS result file
head -1 "${OUT_DIR}/${TRAIT}_gwa.chr1.txt" > "${GWAS_OUT}"
for i in $(seq 1 "${NUM_CHROMOSOMES}"); do
    tail -n +2 "${OUT_DIR}/${TRAIT}_gwa.chr${i}.txt" >> "${GWAS_OUT}"
done

echo "[$(date)] GWAS complete: ${GWAS_OUT}"
