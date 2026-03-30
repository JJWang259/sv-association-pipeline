# Functional enrichment analysis with GEMRICH for fine-mapped signals
# Estimates genetic effect enrichments across functional annotation categories
# using maximum likelihood estimation on BFMAP fine-mapping summary statistics
# Dependencies: gemrich, data.table
#   Install gemrich: devtools::install_github("jiang18/gemrich")
# Author: Junjian Wang
# Initiate date: 03/24/2026
# Current date: 03/24/2026
# =============================================================================

library(data.table)
library(gemrich)

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

FMAP_FILE     <- "/path_to_the_file/fmap_all.csv"    # Output from 02e
BED_FILE      <- "/path_to_the_file/annotation.bed"  # BED file: chr, start, end, category

# SNP universe for background proportion calculation
# CSV with two columns: chr, pos (no header needed if added via command below)
# Generate from PLINK1 BIM : awk '{print $1","$4}' geno.bim | sed '1s/^/chr,pos\n/' > snplist.csv
# Generate from PLINK2 PVAR: awk 'NR>1 {print $1","$2}' geno.pvar | sed '1s/^/chr,pos\n/' > snplist.csv
SNP_LIST_FILE <- "/path_to_the_file/snplist.csv"

PV_THRESHOLD <- 5e-5   # P-value threshold for enrichment analysis

# Trait group — edit trait names to match your fmap_all.csv
GROUP_NAME  <- "Production"
TRAIT_GROUP <- c("MilkYield", "FatYield", "ProteinYield")

# Category priority order for mutually exclusive annotation assignment
# Categories listed earlier take precedence when a SNP overlaps multiple categories
# Default follows standard genomic feature hierarchy: SV > CDS > UTR > promoter > intron > remaining
CATEGORY_LIST <- c("SV", "CDS", "UTR", "promoter", "intron")

# Output file prefix — default: <GROUP_NAME>.pv<PV_THRESHOLD> in current directory
# Change to a full path to write outputs elsewhere, e.g. "/path/gemrich/Production.pv5e-5"
OUT_PREFIX <- paste0(GROUP_NAME, ".pv", PV_THRESHOLD)

# =============================================================================

dir.create(dirname(OUT_PREFIX), recursive = TRUE, showWarnings = FALSE)

# Load data
cat(sprintf("[%s] Loading inputs\n", Sys.time()))
fmap    <- fread(FMAP_FILE)
bed     <- fread(BED_FILE, header = FALSE)
snplist <- fread(SNP_LIST_FILE)   # Must contain columns: chr, pos

# Map fine-mapped SNPs to annotation categories
# multi_cat = TRUE assigns each SNP to one mutually exclusive category by priority order
cat(sprintf("[%s] Mapping SNPs to annotations\n", Sys.time()))
snp2annot <- map_snp_annotation(fmap, bed,
                                category_list = CATEGORY_LIST,
                                multi_cat     = TRUE)

# Compute background category proportions from full SNP universe
cat(sprintf("[%s] Computing background category proportions\n", Sys.time()))
cat_prop <- calc_snp_category_prop(snplist, bed, category_list = CATEGORY_LIST)

# Run enrichment analysis
fmap_sub <- fmap[Trait %in% TRAIT_GROUP]
if (nrow(fmap_sub) == 0) stop(sprintf("No fine-mapped variants found for group: %s", GROUP_NAME))

cat(sprintf("[%s] Running enrichment: group=%s | pv=%s\n", Sys.time(), GROUP_NAME, PV_THRESHOLD))
mle_result <- estimate_category_enrichment(
    fmap_sub, snp2annot, cat_prop,
    pvalue_threshold = PV_THRESHOLD
)

# Save outputs
write.csv(mle_result$enrichment_mle,  paste0(OUT_PREFIX, ".enrichment.csv"), quote = FALSE, row.names = FALSE)
write.csv(mle_result$prob_mle,        paste0(OUT_PREFIX, ".mle.csv"),        quote = FALSE, row.names = FALSE)
write.csv(mle_result$prob_cov_matrix, paste0(OUT_PREFIX, ".var.csv"),        quote = FALSE, row.names = FALSE)
write.csv(mle_result$loglik,          paste0(OUT_PREFIX, ".logll.csv"),      quote = FALSE, row.names = FALSE)
write.csv(mle_result$counts,          paste0(OUT_PREFIX, ".counts.csv"),     quote = FALSE, row.names = FALSE)

cat(sprintf("[%s] Done -> %s.*\n", Sys.time(), OUT_PREFIX))