# Aggregate BFMAP fine-mapping outputs into a deduplicated summary table
# Dependencies: data.table
# Author: Junjian Wang
# Initiate date: 03/24/2026
# Current date: 03/24/2026
# =============================================================================

library(data.table)

# =============================================================================
# User-defined inputs — edit this section before running
# =============================================================================

TRAIT_LIST <- c("MilkYield", "FatYield", "ProteinYield")

RANGE_FILE_PATTERN <- "/path_to_the_file/{TRAIT}_bfmap/{TRAIT}.range.csv"
BFMAP_FILE_PATTERN <- "/path_to_the_file/{TRAIT}_bfmap/{REGION_ID}.csv"

OUT_FILE   <- "/path_to_the_file/fmap_all.csv"

NUM_CORES  <- 4   # Parallel cores for loading BFMAP outputs across traits

# =============================================================================

library(parallel)

load_trait <- function(trt) {
    range_file <- gsub("\\{TRAIT\\}", trt, RANGE_FILE_PATTERN)
    if (!file.exists(range_file) || file.size(range_file) == 0) {
        cat(sprintf("[INFO] Skipping trait %s: no candidate regions found\n", trt))
        return(NULL)
    }
    rg <- fread(range_file)

    out <- rbindlist(lapply(seq_len(nrow(rg)), function(i) {
        fmap_file <- gsub("\\{REGION_ID\\}", rg$region_id[i],
                     gsub("\\{TRAIT\\}", trt, BFMAP_FILE_PATTERN))
        if (!file.exists(fmap_file)) {
            cat(sprintf("  [WARNING] Not found: %s\n", fmap_file))
            return(NULL)
        }
        fread(fmap_file)
    }), fill = TRUE)

    if (nrow(out) == 0) return(NULL)

    out[, signal := signal + cumsum(c(0, diff(signal) < 0))]
    out[, Trait  := trt]
    out
}

outall <- rbindlist(mclapply(TRAIT_LIST, load_trait, mc.cores = NUM_CORES), fill = TRUE)

if (nrow(outall) == 0)
    stop("No fine-mapping results found. Check that 02d ran successfully.")

outall[, signal := paste0(Trait, "_", signal)]

# Deduplicate: keep signal with lowest Pval at SNPindex == 0 per (Trait, SNPname)
keep_sigs  <- outall[SNPindex == 0][
    , .SD[which.min(Pval)], by = .(Trait, SNPname)
]$signal
outall     <- outall[signal %in% keep_sigs]

# Renumber signals sequentially within each trait
outall[, signal := paste0(Trait, "_",
    as.integer(factor(signal, levels = unique(signal))) - 1L),
    by = Trait]

fwrite(outall, OUT_FILE, quote = FALSE, row.names = FALSE)

cat(sprintf("[%s] Done: %d variants, %d signals -> %s\n",
            Sys.time(), nrow(outall), length(unique(outall$signal)), OUT_FILE))