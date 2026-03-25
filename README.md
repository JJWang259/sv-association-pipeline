# Genomics Association Pipeline

A collection of R and Bash scripts for large-scale genotype-phenotype association analysis, implementing the analytical framework used across multiple livestock genomics studies (see [Publications](#publications)).

The pipeline covers three steps: **genome-wide association studies (GWAS)**, **Bayesian fine-mapping**, and **functional enrichment analysis**. Originally developed for SV- and SNP-based GWAS in cattle and pigs, it is broadly applicable to any diploid population with genotype and phenotype data.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Input Data](#input-data)
- [Pipeline Steps](#pipeline-steps)
  - [Step 1: GWAS with SLEMM](#step-1-gwas-with-slemm)
  - [Step 2: Fine-Mapping with BFMAP](#step-2-fine-mapping-with-bfmap)
  - [Step 3: Functional Enrichment with MPH](#step-3-functional-enrichment-with-mph)
- [Output Files](#output-files)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

```
Genotypes (SNPs / SVs, PLINK binary)
         │
         ▼
┌──────────────────────┐
│  Step 1a             │  SLEMM --lmm
│  LMM model fitting   │  Fit variance components per trait
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 1b             │  slemm_gwa
│  Association scan    │  Genome-wide association statistics
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2a             │  R: identify_candidate_regions.R
│  Candidate regions   │  Define fine-mapping windows from GWAS hits
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2b             │  BFMAP (Bayesian fine-mapping)
│  Fine-mapping        │  Posterior inclusion probabilities (PIP) per variant
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2c             │  R: summarise_finemapping.R 
│  Summarise results   │  Aggregate PIPs across all traits and signals
└─────────┬────────────┘
          │
          ├─────────────────────────────────────┐
          ▼                                     ▼
┌──────────────────────┐           ┌──────────────────────┐
│  Step 3a + 3b        │           │  Step 3c             │
│  MPH MINQUE          │           │  GEMRICH             │
│  Polygenic variance  │           │  Large-effect MLE    │
│  partitioning across │           │  enrichment of fine- │
│  annotation bins     │           │  mapped signals      │
└──────────────────────┘           └──────────────────────┘
```

---

## Repository Structure

```
genomics-association-pipeline/
├── README.md
├── data/
│   ├── genotypes/                        # PLINK binary files (.bed/.bim/.fam)
│   ├── phenotypes/                       # Per-trait phenotype CSVs
│   ├── annotations/                      # Functional annotation matrix (snpinfo.csv)
│   └── covariates.csv                    # Shared covariate file
├── scripts/
│   ├── 01a_lmm_fit.sh                    # Step 1a: LMM model fitting (SLEMM)
│   ├── 01b_gwas_association.sh           # Step 1b: Chromosome-wise GWAS scan
│   ├── 02a_identify_candidate_regions.R  # Step 2a: Define candidate regions
│   ├── 02b_finemapping_bfmap.R           # Step 2b: BFMAP fine-mapping per trait
│   ├── 02c_summarise_finemapping.R       # Step 2c: Aggregate fine-mapping results
│   ├── 03a_prepare_mph_covariates.R      # Step 3a: Extract lead-SNP covariates
│   ├── 03b_run_mph.sh                    # Step 3b: MPH enrichment
│   └── 03c_gemrich_enrichment.R          # Step 3c: GEMRICH large-effect enrichment
└── results/
    ├── gwas/                             # Per-trait GWAS summary statistics
    ├── finemapping/                      # Per-region BFMAP output and PIP tables
    └── enrichment/
        ├── covariates/                   # Lead-SNP covariate files for MPH
        ├── mph/                          # MPH variance component estimates
        └── gemrich/                      # GEMRICH enrichment results
```

---

## Dependencies

### Software

| Tool | Version | Purpose | Link |
|------|---------|---------|------|
| SLEMM | any | GWAS | [jiang18/slemm](https://github.com/jiang18/slemm) |
| BFMAP | ≥ 0.65 | Bayesian fine-mapping | [jiang18/bfmap](https://github.com/jiang18/bfmap) |
| MPH | any | Polygenic variance partitioning | [jiang18/mph](https://github.com/jiang18/mph) |
| GEMRICH | any | Large-effect enrichment analysis | [jiang18/gemrich](https://github.com/jiang18/gemrich) |
| PLINK | any | Genotype data processing | [cog-genomics.org/plink](https://www.cog-genomics.org/plink/) |
| R | ≥ 4.1 | Data processing and scripting | [r-project.org](https://www.r-project.org) |

### R Packages

```r
install.packages(c("data.table", "parallel"))

# GEMRICH (install from GitHub)
devtools::install_github("jiang18/gemrich")
```
---

## Input Data

### Required files

| File | Format | Description |
|------|--------|-------------|
| `data/genotypes/geno_model.*` | PLINK binary | Model SNP genotypes for LMM fitting |
| `data/genotypes/geno_tests.*` | PLINK binary | Genotypes for GWAS and fine-mapping |
| `data/genotypes/grm.*` | Binary GRM | Genomic relationship matrix for BFMAP (computed from model SNPs) |
| `data/genotypes/grm_list.txt` | Text | List of partitioned GRM paths for MPH |
| `data/phenotypes/<trait>.csv` | CSV | Per-trait phenotype file (one animal per row) |
| `data/covariates.csv` | CSV | Shared covariates (e.g., intercept, batch effects) |

### Phenotype file format
```
IID,<trait>[,<reliability>]
HO123456,0.312,0.85
HO123457,1.120,0.91
```

The reliability column is optional. If included, specify its column name in
`ERROR_WEIGHT_NAME` in `01a_lmm_fit.sh`.


### Annotation file format (`data/annotations/snpinfo.csv`)

One row per SNP, one column per annotation category. Cell values are category labels.

```
SNP,V1,V2,V3,...
rs123,1,0,1,...
rs456,0,1,0,...
```

---

## Pipeline Steps

### Step 1: GWAS with SLEMM

GWAS is run in two substeps. First, a linear mixed model (LMM) is fit per trait to estimate variance components and compute the genomic relationship structure (`01a`). The fitted model is then used to run genome-wide association tests (`01b`), with results concatenated into a single per-trait file.

**Edit paths** at the top of each script (working directory, bfile prefix, phenotype directory, executable paths, number of threads, number of chromosomes).

```bash
# Step 1a — LMM fitting
bash scripts/01a_lmm_fit.sh

# Step 1b — Association scan (requires Step 1a output)
bash scripts/01b_gwas_association.sh
```

For a full description of SLEMM parameters, see the
[SLEMM documentation](https://github.com/jiang18/slemm).

**Output:** `results/gwas/<trait>_GWAS_All.txt` — genome-wide association statistics (columns: CHROM, POS, ID, REF, ALT, CHISQ, P, etc.)

---

### Step 2: Fine-Mapping with BFMAP

Candidate genomic regions are defined from GWAS hits, then submitted to BFMAP for Bayesian fine-mapping. Results are aggregated across all traits into a single summary table.

**Edit paths** at the top of each R script.

```bash
# Step 2a — Define candidate regions from GWAS (all traits)
Rscript scripts/02a_identify_candidate_regions.R

# Step 2b — Run BFMAP fine-mapping (one trait at a time, parallelise as needed)
Rscript scripts/02b_finemapping_bfmap.R MilkYield
Rscript scripts/02b_finemapping_bfmap.R FatYield
# ... repeat for each trait, or loop:
tail -n +2 config/traits.csv | cut -d, -f1 | \
  xargs -P 4 -I{} Rscript scripts/02b_finemapping_bfmap.R {}

# Step 2c — Aggregate fine-mapping results across all traits
Rscript scripts/02c_summarise_finemapping.R
```

**Key parameters in `02a_identify_candidate_regions.R`:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p1` | 5e-7 | Primary significance threshold for lead SNP |
| `p2` | 5e-6 | Secondary threshold for counting supporting markers |
| `scan` | 5,000,000 bp | Minimum gap to start a new region |
| `window` | 1,000,000 bp | Flanking window around lead variant |
| `min_markers2` | 3 | Minimum secondary markers to retain a region |

**Key parameters in `02b_finemapping_bfmap.R`:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pv_sig` | 0.05 | P-value cutoff for SNP list supplied to BFMAP |
| `min_markers2` | 3 | Minimum secondary markers to proceed with fine-mapping |
| `MEFF` | 100 | BFMAP maximum number of causal variants |
| `NUM_THREADS` | 15 | CPU threads for BFMAP |

**Output:**
- `results/finemapping/<trait>/<trait>.range` — candidate region coordinates
- `results/finemapping/<trait>/<peak>.csv` — per-region BFMAP results (PIP per variant)
- `results/finemapping/fmap_all_traits.csv` — aggregated fine-mapping table across all traits

---

### Step 3: Functional Enrichment with MPH

Fine-mapped lead variants are used as fixed-effect covariates to condition out large-effect QTL signals before partitioning residual genetic variance across functional annotation categories with MPH MINQUE.

```bash
# Step 3a — Prepare per-trait lead-SNP covariate files
Rscript scripts/03a_prepare_mph_covariates.R

# Step 3b — Run MPH MINQUE variance partitioning
bash scripts/03b_run_mph.sh
```

**Key parameters in `03b_run_mph.sh`:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_GRMS` | 21 | Number of GRM components in the partition |
| `NUM_THREADS` | 20 | CPU threads for MPH |

The GRM list (`data/genotypes/grm_list.txt`) should point to partitioned GRMs corresponding to functional annotation categories (e.g., CDS SNPs, UTR SNPs, intronic SNPs, intergenic SNPs, SVs). See [MPH documentation](https://github.com/jiang18/mph) for instructions on building partitioned GRMs.

**Output:**
- `results/enrichment/mph/<trait>.vc` — variance component estimates per annotation
- `results/enrichment/mph/logs/<trait>.log` — MPH log
- `results/enrichment/covariates/<trait>_QTL.csv` — lead-SNP covariate files

---

## Output Files

| File | Description |
|------|-------------|
| `results/gwas/<trait>_GWAS_All.txt` | Full GWAS summary statistics (all chromosomes) |
| `results/gwas/<trait>/<trait>_gwa.chr*.txt` | Per-chromosome association results |
| `results/finemapping/<trait>/<trait>_summary.txt` | Candidate region table from Step 2a |
| `results/finemapping/<trait>/<trait>.range` | Region coordinates used in BFMAP |
| `results/finemapping/<trait>/<peak>.csv` | Per-region BFMAP output (PIP per variant) |
| `results/finemapping/fmap_all_traits.csv` | Aggregated fine-mapping results |
| `results/enrichment/covariates/<trait>_QTL.csv` | Lead-SNP covariates for MPH |
| `results/enrichment/mph/<trait>.vc` | Variance component estimates per annotation |

---

## Citation

This pipeline was developed for and used in the following studies. If you find it helpful, please consider citing the relevant paper(s):

- **Yang L†, Wang J†, Kuhn K, Li W, Zanton G, Neupane M, Boschiero C, Cole JB, Li B, Li C, Baldwin RL VI, Van Tassell CP, Rosen BD, Smith TPL, Jiang J\*, Fang L\*, Ma L\*, Liu GE\*.** An enhanced pangenome reference panel enables accurate imputation of structural variation and large-scale genotype-phenotype analyses in dairy cattle. *(under review)*

- **[Authors].** Long-read sequencing reveals the impact of structural variation on gene expression and complex traits in pigs. *(citation)*

- **[Authors].** An integrated multi-tissue atlas of epigenomic landscapes and regulatory elements in the bovine genome. *(citation)*

Please also cite the underlying tools:

- **SLEMM:** Jiang J, et al. (2023). *Bioinformatics*, 39(4). https://doi.org/10.1093/bioinformatics/btad127
- **BFMAP:** Jiang J, Cole JB, Freebern E, et al. (2019). *Communications Biology*, 2, 212. https://doi.org/10.1038/s42003-019-0454-y
- **MPH:** Jiang J, et al. (2024). *Bioinformatics*, 40(5). https://doi.org/10.1093/bioinformatics/btae298
- **GEMRICH:** Jiang J, et al. https://github.com/jiang18/gemrich
- **PLINK:** Chang CC, et al. (2015). *GigaScience*, 4(1). https://doi.org/10.1186/s13742-015-0047-8

---

## Contact

**Junjian Wang**  
Department of Animal Science, North Carolina State University  
📧 jwang259@ncsu.edu  
🔗 https://github.com/JJWang259

For bug reports or questions, please open an [Issue](https://github.com/JJWang259/sv-association-pipeline/issues).
