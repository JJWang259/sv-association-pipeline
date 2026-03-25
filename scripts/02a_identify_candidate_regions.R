# Identify fine-mapping candidate regions from GWAS summary statistics
# Dependencies: data.table
# Author: Junjian Wang
# Initiate date: 03/24/2026
# Current date: 03/24/2026
# =============================================================================

library(data.table)

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

TRAIT     <- "MilkYield"                               # Trait name
GWAS_FILE <- "/path_to_the_file/MilkYield_GWAS.txt"   # GWAS summary statistics
OUT_FILE  <- "/path_to_the_file/MilkYield_summary.csv" # Output summary table

# Column names in the GWAS file — adjust to match your software's output
# Note: PLINK2 uses "#CHROM" as the chromosome column name. Although "#" looks
#       like a comment character, fread() reads it correctly as a literal column
#       name. Set CHR_COL <- "#CHROM" and it will be matched automatically.
#
# SLEMM + PLINK2 input: CHR_COL = "#CHROM", POS_COL = "POS", P_COL = "P",  ID_COL = "ID",  CHISQ_COL = "CHISQ"
# SLEMM + PLINK1 input: CHR_COL = "CHR",    POS_COL = "BP",  P_COL = "P",  ID_COL = "SNP", CHISQ_COL = "CHISQ"
# GCTA output:          CHR_COL = "Chr",    POS_COL = "bp",  P_COL = "p",  ID_COL = "SNP", CHISQ_COL = "CHISQ" (pre-compute as (b/se)^2 before running)
CHR_COL   <- "#CHROM"   # default: SLEMM with PLINK2 input
POS_COL   <- "POS"
P_COL     <- "P"
ID_COL    <- "ID"
CHISQ_COL <- "CHISQ"    # column name for chi-square statistic

# Region definition parameters
p1                   <- 5e-7   # Primary significance threshold (lead variant cutoff)
p2                   <- 5e-6   # Secondary threshold (for counting supporting markers)
scan                 <- 5e6    # Minimum distance (bp) to define a new region
min_sig_variants     <- 1      # Minimum variants passing p1 to retain a region
min_secondary_variants <- 3    # Minimum variants passing p2 to retain a region

# =============================================================================

# -----------------------------------------------------------------------------
# Function: identify_candidate_regions
#
# Arguments:
#   gwas         - data.table of GWAS results (already loaded)
#   trait        - trait name (used for peak name labelling)
#   chr_col      - column name for chromosome
#   pos_col      - column name for base-pair position
#   p_col        - column name for p-value
#   id_col       - column name for variant ID
#   chisq_col    - column name for chi-square statistic (pre-compute as (b/se)^2 for GCTA)
#   p1           - primary significance threshold for lead variant
#   p2           - secondary threshold for counting supporting markers
#   scan                   - minimum gap (bp) between independent regions
#   min_sig_variants       - minimum variants passing p1 to retain a region (default: 1)
#   min_secondary_variants - minimum variants passing p2 to retain a region (default: 3)
#
# Returns:
#   data.table of candidate regions, or NULL if none found
#   Peak names are formatted as <trait>_<chr>:<lead_pos>
#   leadSNP and leadPOS columns contain all significant variants in the block
#   separated by ";" — used downstream to define fine-mapping boundaries
# -----------------------------------------------------------------------------
identify_candidate_regions <- function(gwas, trait,
                                       chr_col                = "#CHROM",
                                       pos_col                = "POS",
                                       p_col                  = "P",
                                       id_col                 = "ID",
                                       chisq_col              = "CHISQ",
                                       p1                     = 5e-7,
                                       p2                     = 5e-6,
                                       scan                   = 5e6,
                                       min_sig_variants       = 1,
                                       min_secondary_variants = 3) {

    # ── Validate columns ──────────────────────────────────────────────────────
    required_cols <- c(chr_col, pos_col, p_col, id_col, chisq_col)
    missing_cols <- setdiff(required_cols, names(gwas))
    if (length(missing_cols) > 0)
        stop(sprintf("Column(s) not found in GWAS file: %s\nAvailable columns: %s",
                     paste(missing_cols, collapse = ", "),
                     paste(names(gwas), collapse = ", ")))

    # ── Standardise column names ──────────────────────────────────────────────
    gwas <- copy(gwas)
    setnames(gwas, c(chr_col, pos_col, p_col, id_col, chisq_col),
                   c("CHROM",  "POS",   "P",   "ID",   "CHISQ"))
    setorder(gwas, CHROM, POS)

    # ── Define regions per chromosome ─────────────────────────────────────────
    # Helper: given a block of significant variants, build one region row
    make_row <- function(blk, full_chr, i, trait, p2) {
        # Lead variant: highest CHISQ; ";" separated if multiple variants share the maximum
        lead_idx  <- which(blk$CHISQ == max(blk$CHISQ))
        lead_pos  <- blk$POS[lead_idx]

        data.table(
            chr                  = i,
            region_start         = min(blk$POS),
            region_end           = max(blk$POS),
            region_range_bp      = max(blk$POS) - min(blk$POS),
            lead_variant         = paste(blk$ID[lead_idx],  collapse = ";"),
            lead_pos             = paste(lead_pos,           collapse = ";"),
            region_id            = paste0(trait, "_", i, ":", lead_pos[1]),
            n_sig_variants       = nrow(blk),
            n_secondary_variants = nrow(full_chr[P <= p2])
        )
    }

    rows <- list()

    for (i in sort(unique(gwas$CHROM))) {
        full_chr <- gwas[CHROM == i]
        gwa      <- full_chr[P <= p1]
        if (nrow(gwa) == 0) next

        gaps     <- c(0, diff(gwa$POS)) > scan
        block_id <- cumsum(gaps)

        for (b in unique(block_id)) {
            blk <- gwa[block_id == b]
            rows <- c(rows, list(make_row(blk, full_chr, i, trait, p2)))
        }
    }

    if (length(rows) == 0) return(NULL)
    result <- rbindlist(rows)
    result <- result[n_sig_variants       >= min_sig_variants &
                     n_secondary_variants >= min_secondary_variants]
    if (nrow(result) == 0) return(NULL)
    result
}

# =============================================================================
# Run
# =============================================================================

cat(sprintf("[%s] Processing trait: %s\n", Sys.time(), TRAIT))

# Read GWAS results
gwas <- fread(GWAS_FILE)
gwas <- gwas[complete.cases(gwas)]

# Identify candidate regions
result <- identify_candidate_regions(
    gwas                   = gwas,
    trait                  = TRAIT,
    chr_col                = CHR_COL,
    pos_col                = POS_COL,
    p_col                  = P_COL,
    id_col                 = ID_COL,
    chisq_col              = CHISQ_COL,
    p1                     = p1,
    p2                     = p2,
    scan                   = scan,
    min_sig_variants       = min_sig_variants,
    min_secondary_variants = min_secondary_variants
)

# Save output
if (is.null(result) || nrow(result) == 0) {
    cat(sprintf("  No significant regions found for trait: %s\n", TRAIT))
} else {
    dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
    fwrite(result, OUT_FILE, sep = ",", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Found %d candidate regions -> %s\n", nrow(result), OUT_FILE))
}

cat(sprintf("[%s] Done.\n", Sys.time()))
