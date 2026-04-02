# Summarise genotype-specific effect estimates from MPH across multiple variants
# Reads BLUE estimates (AA and BB effects relative to AB baseline) and allele
# information from 04a and 04c outputs for all variants in the input list
# Dependencies: data.table
# Author: Junjian Wang
# Initiate date: 04/01/2026
# Current date: 04/01/2026
# =============================================================================

library(data.table)

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

# Variant list file: tab-delimited, columns: trait, variant
# Adjust col.names and sep if your file has different columns or format
VARIANT_LIST <- "/path_to_the_file/variant_list.txt"
vlist <- fread(VARIANT_LIST, header = FALSE, col.names = c("trait", "variant"))

# Path patterns — use {TRAIT} and {VARIANT} as placeholders
BLUE_FILE_PATTERN   <- "/path_to_the_file/{TRAIT}.{VARIANT}.mq.blue.csv"   # Output from 04c
ALLELE_FILE_PATTERN <- "/path_to_the_file/{VARIANT}.allele_coding.csv"     # Output from 04a

OUT_FILE <- "/path_to_the_file/variant_effects_summary.csv"

# =============================================================================

cat(sprintf("[%s] Summarising %d variant-trait combinations\n", Sys.time(), nrow(vlist)))

results <- rbindlist(lapply(seq_len(nrow(vlist)), function(i) {
    trait   <- vlist$trait[i]
    variant <- vlist$variant[i]

    blue_file   <- gsub("\\{TRAIT\\}",   trait,
                   gsub("\\{VARIANT\\}", variant, BLUE_FILE_PATTERN))
    allele_file <- gsub("\\{VARIANT\\}", variant, ALLELE_FILE_PATTERN)

    if (!file.exists(blue_file)) {
        cat(sprintf("  [WARNING] BLUE file not found: %s\n", blue_file))
        return(NULL)
    }
    if (!file.exists(allele_file)) {
        cat(sprintf("  [WARNING] Allele coding file not found: %s\n", allele_file))
        return(NULL)
    }

    blue   <- fread(blue_file)
    allele <- fread(allele_file)

    data.table(
        trait     = trait,
        variant   = variant,
        A_allele  = allele$A_allele[1],
        B_allele  = allele$B_allele[1],
        freq_A    = allele$freq_A[1],
        HWE_pval  = allele$HWE_pval[1],
        AA_effect = blue$blue[2],
        AA_se     = blue$se[2],
        BB_effect = blue$blue[3],
        BB_se     = blue$se[3]
    )
}), fill = TRUE)

fwrite(results, OUT_FILE, quote = FALSE, row.names = FALSE)

cat(sprintf("[%s] Done: %d results written -> %s\n",
            Sys.time(), nrow(results), OUT_FILE))