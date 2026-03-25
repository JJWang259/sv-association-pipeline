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

TRAIT="MilkYield"                         # Trait name (must match 01a_lmm_fit.sh)

BFILE="/path_to_the_file/geno_sv_snp"    # PLINK binary prefix for combined SNP+SV genotypes
LMM_PREFIX="/path_to_the_file/${TRAIT}"  # Output prefix from 01a_lmm_fit.sh (without extension)

NUM_CHROMOSOMES=29                        # Number of autosomes (cattle=29, pig=18)
export OMP_NUM_THREADS=10

# Output
OUT_DIR="/path_to_the_file/${TRAIT}"     # Directory for per-chromosome result files
GWAS_OUT="/path_to_the_file/${TRAIT}_GWAS_All.txt"  # Concatenated GWAS summary file

# =============================================================================

mkdir -p "${OUT_DIR}"

# Check that LMM model file exists
if [[ ! -f "${LMM_PREFIX}.qs" ]]; then
    echo "[ERROR] LMM model not found: ${LMM_PREFIX}.qs"
    echo "        Did you run 01a_lmm_fit.sh first?"
    exit 1
fi

echo "[$(date)] Running GWAS for trait: ${TRAIT}"

# Run per-chromosome association scan
for i in $(seq 1 "${NUM_CHROMOSOMES}"); do
    echo "  [$(date)] Chromosome ${i} / ${NUM_CHROMOSOMES}"
    "${SLEMM_GWA}" \
        --pfile  "${BFILE}" \
        --slemm  "${LMM_PREFIX}" \
        --out    "${OUT_DIR}/${TRAIT}_gwa.chr${i}.txt" \
        --chr    "${i}"
done

# Concatenate per-chromosome files into a single GWAS result file
head -1 "${OUT_DIR}/${TRAIT}_gwa.chr1.txt" > "${GWAS_OUT}"
for i in $(seq 1 "${NUM_CHROMOSOMES}"); do
    tail -n +2 "${OUT_DIR}/${TRAIT}_gwa.chr${i}.txt" >> "${GWAS_OUT}"
done

echo "[$(date)] GWAS complete: ${GWAS_OUT}"
