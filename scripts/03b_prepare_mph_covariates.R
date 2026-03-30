# Prepare per-trait covariate files for MPH using fine-mapped lead variants
# Lead variants (SNPindex == 0) per signal are used as fixed-effect covariates
# to condition out large-effect QTL signals before polygenic variance partitioning
# Generates one covariate file per trait found in the fine-mapping summary table
# Dependencies: data.table
# Author: Junjian Wang
# Initiate date: 03/24/2026
# Current date: 03/24/2026
# =============================================================================

library(data.table)

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

FMAP_FILE <- "/path_to_the_file/fmap_all.csv"         # Output from 02e
GENO_FILE <- "/path_to_the_file/chr{CHR}.imputed"     # Use {CHR} for per-chromosome files
                                                       # or a plain path for a merged file
OUT_DIR   <- "./covariates"                            # Output directory for covariate files

# PLINK format — choose one:
# PLINK2 (pgen/pvar/psam): PLINK = "plink2", PLINK_FLAG = "--pfile", RECODE_FLAG = "--export A"
# PLINK1 (bed/bim/fam):    PLINK = "plink",  PLINK_FLAG = "--bfile", RECODE_FLAG = "--recodeA"
PLINK       <- "plink2"
PLINK_FLAG  <- "--pfile"
RECODE_FLAG <- "--export A"

# =============================================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

fmap       <- fread(FMAP_FILE)
trait_list <- unique(fmap$Trait)
cat(sprintf("[%s] Preparing MPH covariates for %d traits\n", Sys.time(), length(trait_list)))

for (trt in trait_list) {
    cat(sprintf("  [%s] Trait: %s\n", Sys.time(), trt))

    snp_dat  <- fmap[SNPindex == 0 & Trait == trt, .(SNPname, Chr)]
    out_file    <- file.path(OUT_DIR, paste0(trt, ".QTL.csv"))
    varlist_out <- file.path(OUT_DIR, paste0(trt, ".covariate.variants.csv"))
    tmp_dir     <- file.path(OUT_DIR, "tmp", trt)
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

    if (nrow(snp_dat) == 0) {
        cat("    No index variants found — writing intercept-only covariate file\n")
        fwrite(data.table(IID = character(), intercept = integer()), out_file, quote = FALSE, row.names = FALSE)
        fwrite(data.table(SNPname = character()), varlist_out, quote = FALSE, row.names = FALSE)
        next
    }

    cat(sprintf("    Index variants: %d across %d chromosome(s)\n",
                nrow(snp_dat), length(unique(snp_dat$Chr))))

    # Extract dosages per chromosome and merge by IID
    raw_list <- lapply(unique(snp_dat$Chr), function(chrom) {
        snps       <- snp_dat[Chr == chrom, SNPname]
        list_file  <- file.path(tmp_dir, sprintf("%s_chr%s.varlist", trt, chrom))
        raw_prefix <- file.path(tmp_dir, sprintf("%s_chr%s_QTL",    trt, chrom))
        fwrite(data.table(SNP = snps), list_file, col.names = TRUE, quote = FALSE)
        system(sprintf("%s %s %s --extract %s %s --out %s",
                       PLINK, PLINK_FLAG, gsub("\\{CHR\\}", chrom, GENO_FILE),
                       list_file, RECODE_FLAG, raw_prefix))
        raw_file <- paste0(raw_prefix, ".raw")
        if (!file.exists(raw_file)) { cat(sprintf("    [WARNING] PLINK failed for chr%s\n", chrom)); return(NULL) }
        raw <- fread(raw_file)
        raw[, c(1L, 3L:6L) := NULL]
        raw
    })

    out <- Reduce(function(a, b) merge(a, b, by = "IID", all = TRUE),
                  Filter(Negate(is.null), raw_list))
    out <- cbind(out[, .(IID)], intercept = 1L, out[, -"IID", with = FALSE])

    fwrite(out, out_file, quote = FALSE, row.names = FALSE)

    # Write variant list showing which SNPs were added as covariates
    fwrite(data.table(SNPname = snp_dat$SNPname), varlist_out, quote = FALSE, row.names = FALSE)

    cat(sprintf("    Written: %s (%d animals, %d covariates)\n",
                out_file, nrow(out), ncol(out) - 1L))
    cat(sprintf("    Covariate variants: %s\n", varlist_out))
}

cat(sprintf("[%s] Done.\n", Sys.time()))

# Remove temporary directory
unlink(file.path(OUT_DIR, "tmp"), recursive = TRUE)