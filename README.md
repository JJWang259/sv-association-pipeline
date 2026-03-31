# Genomics Association Pipeline

A collection of R and Bash scripts for large-scale genotype-phenotype association analysis, implementing the analytical framework used across multiple livestock genomics studies (see [Publications](#publications)).

The pipeline covers three steps: **genome-wide association studies (GWAS)**, **Bayesian fine-mapping**, and **functional enrichment analysis**. Originally developed for genetic variant association analyses in farm animals, the
framework may also be applicable to other populations with individual-level genotype and phenotype data.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Input Data](#input-data)
- [Pipeline Steps](#pipeline-steps)
  - [Step 1: GWAS with SLEMM](#step-1-gwas-with-slemm)
  - [Step 2: Fine-Mapping with BFMAP](#step-2-fine-mapping-with-bfmap)
  - [Step 3: Functional Enrichment](#step-3-functional-enrichment)
- [Output Files](#output-files)
- [Publications](#publications)
- [Contact](#contact)

---

## Overview

```
Genotypes (SNPs / SVs, PLINK binary), phenotypes, and annotation files
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
│  Step 3a             │           │  Step 3b + 3c + 3d   │
│  GEMRICH             │           │  MPH                 │
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
└── scripts/
├── 01a_lmm_fit.sh                    # Step 1a: LMM model fitting (SLEMM)
├── 01b_gwas_association.sh           # Step 1b: GWAS (slemm_gwa)
├── 02a_identify_candidate_regions.R  # Step 2a: Define candidate regions
├── 02b_construct_bfmap_grm.sh        # Step 2b: Construct BFMAP GRM
├── 02c_estimate_heritability.sh      # Step 2c: Estimate heritability with BFMAP
├── 02d_finemapping_bfmap.sh          # Step 2d: BFMAP forward selection fine-mapping
├── 02e_summarise_finemapping.R       # Step 2e: Aggregate fine-mapping results
├── 03a_gemrich_enrichment.R          # Step 3a: GEMRICH large-effect enrichment
├── 03b_prepare_mph_covariates.R      # Step 3b: Extract lead-variant covariates
├── 03c_construct_mph_grms.sh         # Step 3c: Construct partitioned GRMs for MPH
└── 03d_run_mph.sh                    # Step 3d: MPH MINQUE variance partitioning
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
| `geno_model.*` | PLINK binary | Model SNP genotypes for LMM fitting and BFMAP GRM construction |
| `geno_test.*` | PLINK binary | Genotypes for association testing and fine-mapping |
| `<trait>.csv` | CSV | Per-trait phenotype file (columns: IID, trait[, reliability]) |
| `snp.info.csv` | CSV | SNP name list for SLEMM variance component estimation and BFMAP GRM construction |
| `snp.info.annot.csv` | CSV | Annotated SNP info for MPH partitioned GRM construction (first column SNP ID, remaining columns define annotation partitions) |
| `annotation.bed` | BED | Functional annotations for GEMRICH (columns: chr, start, end, category) |
| `snplist.csv` | CSV | SNP universe for GEMRICH background proportion calculation (columns: chr, pos) |

### Phenotype file format
```
IID,<trait>[,<reliability>]
HO123456,0.312,0.85
HO123457,1.120,0.91
```

The reliability column is optional. If included, specify its column name in `ERROR_WEIGHT_NAME` in `01a_lmm_fit.sh`, `02c_estimate_heritability.sh`, `02d_finemapping_bfmap.sh`, and `03d_run_mph.sh`.



### SNP info file format (`snp.info.csv`)

A single-column CSV listing SNP IDs used for variance component estimation
in SLEMM and GRM construction in BFMAP. No header required.
```
1_112_C
1_370_G
1_689_ACCTG
```

### Annotation file format (`snp.info.annot.csv`)

One row per SNP, one column per annotation partition. SNPs belonging to a partition have value `1`; all other cells are empty.

```
SNP,V1,V2,V3,...
rs123,,1,,...
rs456,1,,,...
```
### SNP list file for GEMRICH (`snplist.csv`)

```
chr,pos
1,112
1,370
```

Generate from PLINK1 BIM: `awk '{print $1","$4}' geno.bim | sed '1s/^/chr,pos\n/' > snplist.csv`

Generate from PLINK2 PVAR: `awk 'NR>1 {print $1","$2}' geno.pvar | sed '1s/^/chr,pos\n/' > snplist.csv`

---

## Pipeline Steps

### Step 1: GWAS with SLEMM

GWAS is run in two substeps. First, a linear mixed model (LMM) is fit per trait to estimate variance components (`01a`). The fitted model is then used to run genome-wide association tests (`01b`), with results concatenated into a single per-trait file.

Edit paths at the top of each script before running.

```bash
bash scripts/01a_lmm_fit.sh
bash scripts/01b_gwas_association.sh
```

`01b` supports both merged and per-chromosome genotype files via a `{CHR}` placeholder in `GENO_TEST_TEMPLATE`.

For a full description of SLEMM parameters, see the [SLEMM documentation](https://github.com/jiang18/slemm).

**Output:** `<trait>_GWAS_All.txt`

---

### Step 2: Fine-Mapping with BFMAP

Fine-mapping is run in five substeps. Candidate regions are first identified from significant GWAS peaks (2a), a GRM is constructed for BFMAP (2b), heritability is estimated per trait (2c), fine-mapping is performed per region (2d), and results are aggregated into a single summary table (2e). Edit paths at the top of each script before running.

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
Aggregates BFMAP outputs across one or more traits into a single summary table.

Output: `fmap_all.csv`

### Step 3: Functional Enrichment

Two complementary analyses are performed using the aggregated fine-mapping
results from Step 2e. GEMRICH estimates the enrichment of large-effect
fine-mapped signals across functional annotation categories. MPH MINQUE
estimates the enrichment of polygenic variance across annotation bins,
conditioning on fine-mapped signals as fixed effects to account for
large-effect QTL contributions.

---

**Step 3a: Large-effect enrichment with GEMRICH**

GEMRICH applies an MLE-based model to estimate enrichment of fine-mapped
large-effect signals across functional annotation categories.

```bash
Rscript scripts/03a_gemrich_enrichment.R
```

Output: `<group>.<pv>.enrichment.csv`, `.mle.csv`, `.var.csv`, `.logll.csv`, `.counts.csv`

For a full description of GEMRICH functions and usage, see the [GEMRICH documentation](https://github.com/jiang18/gemrich).

---

**Step 3b-d: Polygenic variance partitioning with MPH**

MPH partitions polygenic genetic variance across functional annotation
bins, conditioned on fine-mapped lead variants as fixed effects. Three
substeps are required.

**3b: Prepare covariates**

Extracts lead variants (`SNPindex == 0`) per signal from the fine-mapping
summary and recodes their genotypes as fixed-effect covariates. Generates
one covariate file per trait automatically.
```bash
Rscript scripts/03b_prepare_mph_covariates.R
```
Output: `<trait>.QTL.csv`

**3c: Construct partitioned GRMs**

Builds one GRM per annotation partition defined in the SNP info file and
generates `grm_list.txt` for use in 3d.
```bash
bash scripts/03c_construct_mph_grms.sh
```

Output: `<partition>.grm.*` files and `grm_list.txt`

**3d: Run MPH MINQUE**
```bash
bash scripts/03d_run_mph.sh
```

Output: `<trait>.mph.mq.*`

For a full description of MPH parameters and variance component output interpretation, see the [MPH documentation](https://jiang18.github.io/mph/).

---

## Output Files

| File | Description |
|------|-------------|
| `<trait>_GWAS_All.txt` | Full GWAS summary statistics from 01b |
| `<trait>_candidate_regions.csv` | Candidate region table from 02a |
| `<trait>_bfmap/<region_id>.csv` | Per-region BFMAP output (PIP per variant) from 02d |
| `fmap_all.csv` | Aggregated fine-mapping results from 02e |
| `<group>.<pv>.enrichment.csv` | GEMRICH enrichment estimates from 03a |
| `<trait>.mph.mq.vc.csv` | MPH variance component estimates per annotation from 03d|

---

## Publications

This pipeline was developed for and used in the following studies. If you find it helpful, please consider citing the relevant paper(s):

- **Yang L†, Wang J†, Kuhn K, Li W, Zanton G, Neupane M, Boschiero C, Cole JB, Li B, Li C, Baldwin RL VI, Van Tassell CP, Rosen BD, Smith TPL, Jiang J\*, Fang L\*, Ma L\*, Liu GE\*.** An enhanced pangenome reference panel enables accurate imputation of structural variation and large-scale genotype-phenotype analyses in dairy cattle. *(under review)*

- **Bao Q†, Wen H†, Zeng L†, Yang L, Wang J, He S, Yin H, Jiang Y, Qu X, Wang Z, Li X, Yang X, Teng J, Zhao P, Zhang D, Liu D, Cao G, Yu T, Ding R, Wang C, Zhang W, Tan Z, Zhu Y, Xia Y, Wang J, Zhao X, Tiezzi F, Gini C, Huang Y, See G, Schwab C, Xu R, Chen Z, Zhao Y, Xiang H, Zhou H, Ding X, Zhang Z, Tang Z, Li K, Maltecca C, Fang L\*, Jiang J\*, Yi G\*.** Long-read sequencing reveals the impact of structural variation on gene expression and complex traits in pigs. *(under review)*


Please also cite the underlying tools:

- **SLEMM:** Cheng J, et al. (2023). *Bioinformatics*, 39(4). https://doi.org/10.1093/bioinformatics/btad127
- **BFMAP:** Jiang J, et al. (2019). *Communications Biology*, 2, 212. https://doi.org/10.1038/s42003-019-0454-y; Wang J, et al. (2025). *Briefings in Bioinformatics*, 26(6). https://doi.org/10.1093/bib/bbaf614
- **MPH:** Jiang J. (2024). *Bioinformatics*, 40(5). https://doi.org/10.1093/bioinformatics/btae298
- **GEMRICH:** Jiang J. https://github.com/jiang18/gemrich
- **PLINK:** Chang CC, et al. (2015). *GigaScience*, 4(1). https://doi.org/10.1186/s13742-015-0047-8

---

## Contact

**Junjian Wang**  
Department of Animal Science, North Carolina State University  
📧 jwang259@ncsu.edu  
🔗 [jjwang259.github.io](https://jjwang259.github.io/)

**Jicai Jiang**  
Department of Animal Science, North Carolina State University  
📧 jicai_jiang@ncsu.edu  
🔗 [cals.ncsu.edu/animal-science/people/jicai-jiang](https://cals.ncsu.edu/animal-science/people/jicai-jiang/)

For bug reports or questions, please open an [Issue](https://github.com/JJWang259/sv-association-pipeline/issues).
