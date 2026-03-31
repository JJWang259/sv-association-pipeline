#!/bin/bash
# Bayesian fine-mapping with BFMAP forward selection for a single trait
# Dependencies: PLINK, BFMAP, awk
# Author: Junjian Wang
# Initiate date: 02/12/2026
# Current date: 03/31/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

BFMAP=~/bin/bfmap
PLINK=plink2

TRAIT="MilkYield"

REGION_FILE="/path_to_the_file/${TRAIT}_candidate_regions.csv"  # Output from 02a
PHENO_FILE="/path_to_the_file/${TRAIT}.csv"                      # Phenotype file
VARCOMP_FILE="/path_to_the_file/${TRAIT}.varcomp.csv"            # Output from 02c (h2 in row 2, col 2)
GRM_FILE="/path_to_the_file/bfmap_grm"                           # BFMAP GRM prefix (output from 02b)
OUT_DIR="/path_to_the_file/${TRAIT}_bfmap"                       # Output directory

# Genotype file(s) for fine-mapping
# For per-chromosome files, use {CHR} as a placeholder for the chromosome number
# For a single merged file, simply provide the path with no {CHR} placeholder
GENO_TEST="/path_to_the_file/chr{CHR}.imputed"

# Fine-mapping parameters
WINDOW=1000000    # Defines the candidate region as ±1 Mb around the top associated position
NUM_THREADS=10

# Optional: reliability/error weight column name in PHENO_FILE
# Set to "" to skip
ERROR_WEIGHT_NAME=""                 # e.g. "reliability" or "R-MilkYield"

# Optional: effective number of SNPs for forward selection stopping criterion
# If not set, BFMAP estimates it automatically using Li & Ji (2005)
# MEFF=100

# =============================================================================

mkdir -p "${OUT_DIR}"
echo "[$(date)] Fine-mapping trait: ${TRAIT}"

# Extract heritability from variance component file (row 2, col 2)
H2=$(awk -F',' 'NR==2 {print $2}' "${VARCOMP_FILE}")
echo "  Heritability (h2): ${H2}"

# Build optional arguments
ERROR_ARGS=""
if [[ -n "${ERROR_WEIGHT_NAME}" ]]; then
    ERROR_ARGS="--error_weight_name ${ERROR_WEIGHT_NAME}"
fi

MEFF_ARG=""
if [[ -n "${MEFF:-}" ]]; then
    MEFF_ARG="--meff ${MEFF}"
fi

# Write region index file header
RANGE_OUT="${OUT_DIR}/${TRAIT}.range.csv"
echo "region_id,chr,win_start,win_end,lead_variant,lead_pos" > "${RANGE_OUT}"

# ── Process each candidate region ────────────────────────────────────────────
# Skip header; fields: chr,region_start,region_end,region_range_bp,
#                      lead_variant,lead_pos,region_id,n_sig_variants,n_secondary_variants
tail -n +2 "${REGION_FILE}" | while IFS=',' read -r \
    chr region_start region_end region_range_bp \
    lead_variant lead_pos region_id n_sig n_sec; do

    echo "  [$(date)] Region: ${region_id} | chr${chr}:${region_start}-${region_end}"

    # Resolve genotype file for this chromosome
    GENO="${GENO_TEST/\{CHR\}/${chr}}"
    echo "    Using geno: ${GENO}"

    # Compute fine-mapping window
    # lead_pos may be ";"-separated if multiple variants share max CHISQ
    # win_start = pmin(min(lead_pos) - WINDOW, region_start), floor at 0
    # win_end   = pmax(max(lead_pos) + WINDOW, region_end)
    MIN_LEAD=$(echo "${lead_pos}" | tr ';' '\n' | sort -n | head -1)
    MAX_LEAD=$(echo "${lead_pos}" | tr ';' '\n' | sort -n | tail -1)

    WIN_START=$(( MIN_LEAD - WINDOW < region_start ? MIN_LEAD - WINDOW : region_start ))
    WIN_START=$(( WIN_START < 0 ? 0 : WIN_START ))
    WIN_END=$(( MAX_LEAD + WINDOW > region_end ? MAX_LEAD + WINDOW : region_end ))

    echo "    Window: chr${chr}:${WIN_START}-${WIN_END}"

    # Extract region genotypes
    "${PLINK}" \
        --pfile   "${GENO}" \
        --from-bp "${WIN_START}" \
        --to-bp   "${WIN_END}" \
        --chr     "${chr}" \
        --make-bed \
        --out     "${OUT_DIR}/${region_id}"

    BIM_FILE="${OUT_DIR}/${region_id}.bim"
    if [[ ! -f "${BIM_FILE}" ]]; then
        echo "  [WARNING] No genotypes extracted for region: ${region_id}"
        continue
    fi

    # Build variant list from extracted BIM file
    VARLIST_FILE="${OUT_DIR}/${region_id}.varlist"
    echo "SNP" > "${VARLIST_FILE}"
    awk '{print $2}' "${BIM_FILE}" >> "${VARLIST_FILE}"
    echo "    Variants in region: $(tail -n +2 "${VARLIST_FILE}" | wc -l)"

    # Forward selection fine-mapping
    "${BFMAP}" \
        --phenotype            "${PHENO_FILE}" \
        --trait                "${TRAIT}" \
        --snp_info_file        "${VARLIST_FILE}" \
        --binary_genotype_file "${OUT_DIR}/${region_id}" \
        --binary_grm           "${GRM_FILE}" \
        --heritability         "${H2}" \
        --output               "${OUT_DIR}/${region_id}" \
        --num_threads          "${NUM_THREADS}" \
        ${ERROR_ARGS} \
        ${MEFF_ARG}

    # Append region to range output
    echo "${region_id},${chr},${WIN_START},${WIN_END},${lead_variant},${lead_pos}" >> "${RANGE_OUT}"

done

echo "[$(date)] Fine-mapping complete. Range file: ${RANGE_OUT}"