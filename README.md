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
- [Publications](#publications)
- [Contact](#contact)

---

## Overview

```
Genotypes (SNPs / SVs, PLINK binary)
         │
         ▼
┌──────────────────────┐
│  Step 1a             │  01a_lmm_fit.sh
│  LMM model fitting   │  Fit variance components per trait with SLEMM --lmm
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 1b             │  01b_gwas_association.sh
│  Association scan    │  Genome-wide association statistics with slemm_gwa
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2a             │  02a_identify_candidate_regions.R
│  Candidate regions   │  Define candidate regions from GWAS hits
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2b             │  02b_construct_bfmap_grm.sh
│  Construct GRM       │  Build BFMAP GRM from model SNPs
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2c             │  02c_estimate_heritability.sh
│  Heritability        │  Estimate h2 per trait with BFMAP
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2d             │  02d_finemapping_bfmap.sh
│  Fine-mapping        │  Posterior inclusion probabilities (PIP) per variant
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2e             │  02e_summarise_finemapping.R
│  Aggregate results   │  Aggregate BFMAP outputs into a single summary table
└─────────┬────────────┘
          │
          ├─────────────────────────────────────┐
          ▼                                     ▼
┌──────────────────────┐           ┌──────────────────────┐
│  Step 3a             │           │  Step 3b + 3c        │
│  GEMRICH             │           │  MPH MINQUE          │
│  Large-effect MLE    │           │  Condition on fine-  │
│  enrichment of fine- │           │  mapped effects;     │
│  mapped signals      │           │  partition polygenic │
│                      │           │  variance across     │
│                      │           │  annotations         │
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
| `data/genotypes/geno_model.*` | PLINK binary | Model SNP genotypes for LMM fitting and GRM construction |
| `data/genotypes/geno_test.*` | PLINK binary | Genotypes for association testing and fine-mapping |
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

GWAS is run in two substeps. First, a linear mixed model (LMM) is fit per trait to estimate variance components (`01a`). The fitted model is then used to run genome-wide association tests (`01b`), with results concatenated into a single per-trait file.

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

Fine-mapping is run in five substeps. Candidate regions are first identified
from GWAS hits (2a), a GRM is constructed for BFMAP (2b), heritability is
estimated per trait (2c), fine-mapping is performed per region (2d), and
results are aggregated into a single summary table (2e).

Edit paths at the top of each script before running.

**2a: Identify candidate regions**
```bash
Rscript scripts/02a_identify_candidate_regions.R
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p1` | 5e-7 | Primary significance threshold for lead variant |
| `p2` | 5e-6 | Secondary threshold for counting supporting markers |
| `scan` | 5,000,000 bp | Minimum gap to define a new region |
| `min_sig_variants` | 1 | Minimum variants passing p1 to retain a region |
| `min_secondary_variants` | 3 | Minimum variants passing p2 to retain a region |

Output: `<trait>_candidate_regions.csv`

**2b: Construct BFMAP GRM**

Note: the BFMAP GRM format is not compatible with MPH's GRM format. This step only needs to be run once per genotype dataset.
```bash
bash scripts/02b_construct_bfmap_grm.sh
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `GRM_TYPE` | 2 | GRM type: 1 = centered, 2 = standardized |

Output: `bfmap_grm.*`

**2c: Estimate heritability**
```bash
bash scripts/02c_estimate_heritability.sh
```

Output: `<trait>.varcomp.csv`

**2d: Fine-mapping**
```bash
bash scripts/02d_finemapping_bfmap.sh
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `WINDOW` | 1,000,000 bp | Defines the candidate region as ±1 Mb around the top associated position for each significant GWAS peak |

Output: per-region BFMAP results (`<region_id>.csv`) and region index (`<trait>.range.csv`) in `<trait>_bfmap/`.

For a full description of BFMAP parameters, see the [BFMAP documentation](https://github.com/jiang18/bfmap).

**2e: Aggregate fine-mapping results**
```bash
Rscript scripts/02e_summarise_finemapping.R
```

Aggregates BFMAP outputs across one or more traits into a single deduplicated
summary table. Edit `TRAIT_LIST` to specify which traits to include.

### Step 3: Functional Enrichment

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

## Publications

This pipeline was developed for and used in the following studies. If you find it helpful, please consider citing the relevant paper(s):

- **Yang L†, Wang J†, Kuhn K, Li W, Zanton G, Neupane M, Boschiero C, Cole JB, Li B, Li C, Baldwin RL VI, Van Tassell CP, Rosen BD, Smith TPL, Jiang J\*, Fang L\*, Ma L\*, Liu GE\*.** An enhanced pangenome reference panel enables accurate imputation of structural variation and large-scale genotype-phenotype analyses in dairy cattle. *(under review)*

- **[Authors].** Long-read sequencing reveals the impact of structural variation on gene expression and complex traits in pigs. *(citation)*

- **[Authors].** An integrated multi-tissue atlas of epigenomic landscapes and regulatory elements in the bovine genome. *(citation)*

Please also cite the underlying tools:

- **SLEMM:** Cheng J, et al. (2023). *Bioinformatics*, 39(4). https://doi.org/10.1093/bioinformatics/btad127
- **BFMAP:** Jiang J, et al. (2019). *Communications Biology*, 2, 212. https://doi.org/10.1038/s42003-019-0454-y
- **MPH:** Jiang J. (2024). *Bioinformatics*, 40(5). https://doi.org/10.1093/bioinformatics/btae298
- **GEMRICH:** Jiang J. https://github.com/jiang18/gemrich
- **PLINK:** Chang CC, et al. (2015). *GigaScience*, 4(1). https://doi.org/10.1186/s13742-015-0047-8

---

## Contact

**Junjian Wang**  
Department of Animal Science, North Carolina State University  
📧 jwang259@ncsu.edu  
🔗 https://github.com/JJWang259

For bug reports or questions, please open an [Issue](https://github.com/JJWang259/sv-association-pipeline/issues).
