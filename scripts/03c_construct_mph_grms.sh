#!/bin/bash
# Construct partitioned GRMs for MPH from a SNP info file
# Each column in the SNP info file defines one annotation partition
# Also generates grm_list.txt for use in 03d_run_mph.sh
# Dependencies: PLINK, MPH
# Author: Junjian Wang
# Initiate date: 02/12/2026
# Current date: 03/31/2026
# =============================================================================

set -euo pipefail

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

PLINK="plink"
MPH="$HOME/bin/mph"

GENO_FILE="/path_to_the_file/geno_test"            # PLINK binary prefix (same genotypes used in GWAS)
SNP_INFO="/path_to_the_file/snp.info.annot.csv"    # SNP info CSV: first column SNP ID,
                                                   # remaining columns define annotation partitions
OUT_DIR="/path_to_the_file/grms"                   # Output directory for GRM files

NUM_AUTOSOMES=29                                   # Number of autosomes (cattle=29, pig=18)
NUM_THREADS=20

# =============================================================================

mkdir -p "${OUT_DIR}"
GRM_LIST="${OUT_DIR}/grm_list.txt"
> "${GRM_LIST}"
SKIPPED=0

# Read partition names from header (skip first column)
IFS=',' read -r -a columns < "${SNP_INFO}"

echo "[$(date)] Constructing $(( ${#columns[@]} - 1 )) partitioned GRMs"

for partition in "${columns[@]:1}"; do
    echo "  [$(date)] Partition: ${partition}"

    EXTRACT="${OUT_DIR}/${partition}.extract.txt"
    OUT="${OUT_DIR}/${partition}"

    # Extract SNPs for this partition (column value == 1)
    awk -F',' -v col="${partition}" \
        'NR==1 { for (j=1;j<=NF;j++) if ($j==col) idx=j }
         NR>1  && $idx==1 { print $1 }' "${SNP_INFO}" > "${EXTRACT}"

    N=$(wc -l < "${EXTRACT}")
    if [[ "${N}" -eq 0 ]]; then
        echo "    [WARNING] No SNPs — skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    echo "    SNPs: ${N}"

    "${PLINK}" --bfile "${GENO_FILE}" --autosome-num "${NUM_AUTOSOMES}" \
        --extract "${EXTRACT}" --threads "${NUM_THREADS}" --make-bed --out "${OUT}"

    "${MPH}" --make_grm \
        --binary_geno "${OUT}" --snp_info "${SNP_INFO}" \
        --snp_weight "${partition}" --num_threads "${NUM_THREADS}" --out "${OUT}"

    # Remove intermediate PLINK files
    rm -f "${OUT}.bed" "${OUT}.bim" "${OUT}.fam" "${OUT}.log"

    echo "${OUT}" >> "${GRM_LIST}"
done

echo "[$(date)] Done. $(wc -l < "${GRM_LIST}") GRMs built, ${SKIPPED} partitions skipped -> ${GRM_LIST}"
