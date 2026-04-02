# Extract genotype dosage for a specific variant and recode as dummy variable
# covariates (AA = homozygous A_allele, BB = homozygous B_allele) with
# heterozygous (AB) as baseline, enabling genotype-specific effect estimation with MPH
# Note: A_allele is the minor allele (PLINK1) or reference allele (PLINK2);
#       B_allele is the major allele (PLINK1) or alternative allele (PLINK2)
# Dependencies: data.table
# Author: Junjian Wang
# Initiate date: 03/24/2026
# Current date: 03/24/2026
# =============================================================================

library(data.table)

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

VARIANT    <- "7_16507712_sv"                          # Variant ID to extract
GENO_FILE  <- "/path_to_the_file/geno_test"            # PLINK binary prefix
OUT_PREFIX <- paste0("/path_to_the_file/", VARIANT)    # Output file prefix

# PLINK format — choose one:
# PLINK2 (pgen/pvar/psam): PLINK = "plink2", PLINK_FLAG = "--pfile", RECODE_FLAG = "--export A"
# PLINK1 (bed/bim/fam):    PLINK = "plink",  PLINK_FLAG = "--bfile", RECODE_FLAG = "--recodeA"
PLINK       <- "plink2"
PLINK_FLAG  <- "--pfile"
RECODE_FLAG <- "--export A"

# =============================================================================

dir.create(dirname(OUT_PREFIX), recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[%s] Extracting genotype for variant: %s\n", Sys.time(), VARIANT))

# Extract variant genotype and compute HWE p-value via PLINK --hardy
system(sprintf("%s %s %s --snp %s %s --hardy --out %s",
               PLINK, PLINK_FLAG, GENO_FILE, VARIANT, RECODE_FLAG, OUT_PREFIX))

raw_file <- paste0(OUT_PREFIX, ".raw")
if (!file.exists(raw_file)) stop(sprintf("PLINK extraction failed: %s", raw_file))

raw       <- fread(raw_file)
geno_col  <- names(raw)[7]
genotypes <- raw[[geno_col]]

# Genotype counts
n_AA  <- sum(genotypes == 2, na.rm = TRUE)
n_AB  <- sum(genotypes == 1, na.rm = TRUE)
n_BB  <- sum(genotypes == 0, na.rm = TRUE)
n_tot <- n_AA + n_AB + n_BB

# A_allele frequency from genotype counts
freq_A <- (2 * n_AA + n_AB) / (2 * n_tot)

# HWE p-value from PLINK --hardy output (exact test, matches PLINK implementation)
# PLINK1 output: .hwe with column P
# PLINK2 output: .hardy with column P
hwe_pval <- tryCatch({
    hwe_file <- paste0(OUT_PREFIX, if (PLINK_FLAG == "--bfile") ".hwe" else ".hardy")
    hwe      <- fread(hwe_file)
    hwe$P[1]
}, error = function(e) NA_real_)

# Read allele information from BIM or PVAR
if (PLINK_FLAG == "--bfile") {
    allele_ref <- fread(paste0(GENO_FILE, ".bim"),
                        col.names = c("chr", "ID", "cm", "pos", "A_allele", "B_allele"))
} else {
    allele_ref <- fread(paste0(GENO_FILE, ".pvar"))
    setnames(allele_ref, c("#CHROM","POS","ID","REF","ALT"),
                         c("chr",   "pos","ID","A_allele","B_allele"), skip_absent = TRUE)
}
snp_info <- allele_ref[ID == VARIANT]
if (nrow(snp_info) > 0) {
    fwrite(data.table(
        variant  = VARIANT,
        A_allele = snp_info$A_allele[1],
        B_allele = snp_info$B_allele[1],
        freq_A   = round(freq_A, 6),
        HWE_pval = signif(hwe_pval, 4)
    ), paste0(OUT_PREFIX, ".allele_coding.csv"), quote = FALSE, row.names = FALSE)
}

# Create dummy variables: AA and BB with AB (heterozygous) as baseline
#   AA coefficient = effect of homozygous A_allele relative to heterozygous
#   BB coefficient = effect of homozygous B_allele relative to heterozygous
out <- data.table(
    IID       = raw$IID,
    intercept = 1L,
    AA        = as.integer(genotypes == 2),
    BB        = as.integer(genotypes == 0)
)

fwrite(out, paste0(OUT_PREFIX, ".covariate.csv"), quote = FALSE, row.names = FALSE)

# Remove intermediate PLINK files
for (ext in c(".raw", ".log", ".nosex", ".hwe", ".hardy")) {
    f <- paste0(OUT_PREFIX, ext)
    if (file.exists(f)) file.remove(f)
}

cat(sprintf("  AA (hom ref): %d | AB (het, baseline): %d | BB (hom alt): %d\n",
            n_AA, n_AB, n_BB))
cat(sprintf("  A_allele freq: %.4f | HWE p-value: %s\n", freq_A, signif(hwe_pval, 4)))
cat(sprintf("[%s] Done: %s.covariate.csv\n", Sys.time(), OUT_PREFIX))