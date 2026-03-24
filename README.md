# Genomics Association Pipeline

A collection of R and Bash scripts for large-scale genotype-phenotype association analysis, implementing the analytical framework used across multiple livestock genomics studies (see [Publications](#publications)).

The pipeline covers three steps: **genome-wide association studies (GWAS)**, **Bayesian fine-mapping**, and **functional enrichment analysis**. Originally developed for SV- and SNP-based GWAS in cattle and pigs, it is broadly applicable to any diploid population with genotype and phenotype data.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Installation](#installation)
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
Imputed SVs + SNPs  (PLINK binary format)
         │
         ▼
┌──────────────────────┐
│  Step 1a             │  SLEMM --lmm
│  LMM model fitting   │  Fit variance components per trait
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 1b             │  slemm_gwa (per chromosome)
│  Association scan    │  Genome-wide SV-trait association statistics
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 2a             │  R: data.table
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
│  Step 2c             │  R: data.table
│  Summarise results   │  Aggregate PIPs across all traits and signals
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 3a             │  R + plink
│  Prepare covariates  │  Fine-mapped lead variants as fixed effects for MPH
└─────────┬────────────┘
          │
          ▼
┌──────────────────────┐
│  Step 3b             │  MPH --minque
│  Enrichment (MPH)    │  Variance partitioning across functional annotations
└──────────────────────┘
```

---

## Repository Structure

```
sv-association-pipeline/
├── README.md
├── config/
│   └── traits.csv                    # Trait list (trait, N, mean, sd)
├── data/
│   ├── genotypes/                    # PLINK binary files (.bed/.bim/.fam)
│   ├── phenotypes/                   # Per-trait phenotype CSVs
│   ├── covariates.csv                # Shared covariate file
│   └── varcomp/                      # Per-trait variance component estimates
├── scripts/
│   ├── 01a_lmm_fit.sh                # Step 1a: LMM model fitting (SLEMM)
│   ├── 01b_gwas_association.sh       # Step 1b: Chromosome-wise GWAS scan
│   ├── 02a_identify_candidate_regions.R  # Step 2a: Define candidate regions
│   ├── 02b_finemapping_bfmap.R       # Step 2b: BFMAP fine-mapping per trait
│   ├── 02c_summarise_finemapping.R   # Step 2c: Aggregate fine-mapping results
│   ├── 03a_prepare_mph_covariates.R  # Step 3a: Extract lead-SNP covariates
│   └── 03b_run_mph.sh                # Step 3b: MPH MINQUE enrichment analysis
└── results/
    ├── gwas/                         # Per-trait GWAS summary statistics
    ├── finemapping/                  # Per-region BFMAP output and PIP tables
    └── enrichment/                   # MPH variance components and enrichment
```

---

## Dependencies

### Software

| Tool | Version | Purpose | Link |
|------|---------|---------|------|
| SLEMM | ≥ 0.93 | LMM fitting and GWAS | [jiang18/slemm](https://github.com/jiang18/slemm) |
| BFMAP | ≥ 0.65 | Bayesian fine-mapping | [jiang18/bfmap](https://github.com/jiang18/bfmap) |
| MPH | latest | Variance partitioning and enrichment | [jiang18/mph](https://github.com/jiang18/mph) |
| PLINK 1.9 | ≥ 1.90 | Genotype extraction (`--recodeA`) | [plink 1.9](https://www.cog-genomics.org/plink/) |
| PLINK 2 | ≥ 2.0 | Genotype extraction (`--make-bed`) | [plink 2](https://www.cog-genomics.org/plink/2.0/) |
| R | ≥ 4.1 | Data processing and scripting | [r-project.org](https://www.r-project.org) |
| Bash | ≥ 4.0 | Shell scripting | — |

### R Packages

```r
install.packages("data.table")
```

---

## Installation

Clone this repository:

```bash
git clone https://github.com/JJWang259/sv-association-pipeline.git
cd sv-association-pipeline
```

Place (or symlink) the SLEMM, BFMAP, and MPH executables in `bin/`:

```bash
mkdir -p bin
cp /path/to/slemm   bin/slemm
cp /path/to/slemm_gwa bin/slemm_gwa
cp /path/to/bfmap   bin/bfmap
cp /path/to/mph     bin/mph
```

Alternatively, ensure they are available in your `$PATH` and update the `SLEMM`, `SLEMM_GWA`, `BFMAP`, and `MPH` variables at the top of each script.

Make shell scripts executable:

```bash
chmod +x scripts/*.sh
```

---

## Input Data

### Required files

| File | Format | Description |
|------|--------|-------------|
| `config/traits.csv` | CSV | Trait list with columns: `trait, N, mean, sd` |
| `data/genotypes/geno_qc.*` | PLINK binary | Quality-controlled SNP genotypes for LMM fitting |
| `data/genotypes/geno_sv_snp.*` | PLINK binary | Merged SNP + imputed SV genotypes for GWAS and fine-mapping |
| `data/genotypes/grm_10k.*` | Binary GRM | Genomic relationship matrix for BFMAP (computed from ~10k markers) |
| `data/genotypes/grm_list.txt` | Text | List of partitioned GRM paths for MPH |
| `data/phenotypes/<trait>.csv` | CSV | Per-trait phenotype file (one animal per row) |
| `data/covariates.csv` | CSV | Shared covariates (e.g., intercept, batch effects) |
| `data/varcomp/<trait>.varcomp.csv` | CSV | Variance component estimates for each trait |

### Phenotype file format

```
IID,<trait>
HO123456,0.312
HO123457,1.120
```

### Trait list format (`config/traits.csv`)

```
trait,N,mean,sd
MilkYield,50299,0.00,1.00
FatYield,50299,0.00,1.00
```

---

## Pipeline Steps

### Step 1: GWAS with SLEMM

GWAS is run in two substeps. First, a linear mixed model (LMM) is fit per trait to estimate variance components and compute the genomic relationship structure (`01a`). The fitted model is then used to run a chromosome-wise association scan across all variants (`01b`), with results concatenated into a single per-trait file.

**Edit paths** at the top of each script (working directory, bfile prefix, phenotype directory, executable paths, number of threads, number of chromosomes).

```bash
# Step 1a — LMM fitting
bash scripts/01a_lmm_fit.sh

# Step 1b — Association scan (requires Step 1a output)
bash scripts/01b_gwas_association.sh
```

**Key parameters in `01a_lmm_fit.sh`:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_THREADS` | 20 | CPU threads for SLEMM |
| `MIN_MAF` | 0.01 | Minimum minor allele frequency |
| `MIN_HWE_PVAL` | 1e-9 | Hardy-Weinberg equilibrium filter |
| `NUM_QF_MARKERS` | 90 | Number of QF markers for SLEMM |

**Key parameters in `01b_gwas_association.sh`:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_CHROMOSOMES` | 29 | Number of autosomes (29 cattle, 18 pig) |
| `OMP_NUM_THREADS` | 10 | OpenMP threads for slemm_gwa |

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

If you use this pipeline, please cite:

**Yang L†, Wang J†, Kuhn K, Li W, Zanton G, Neupane M, Boschiero C, Cole JB, Li B, Li C, Baldwin RL VI, Van Tassell CP, Rosen BD, Smith TPL, Jiang J\*, Fang L\*, Ma L\*, Liu GE\*.** An enhanced pangenome reference panel enables accurate imputation of structural variation and large-scale genotype-phenotype analyses in dairy cattle. *(under review)*

Please also cite the underlying tools:

- **SLEMM:** Jiang J, et al. (2023). *Genetics*, 224(1). https://doi.org/10.1093/genetics/iyad014
- **BFMAP:** Jiang J, et al. https://github.com/jiang18/bfmap
- **MPH / GEMRICH:** Jiang J, et al. https://github.com/jiang18/mph | https://github.com/jiang18/gemrich
- **PLINK 1.9/2:** Chang CC, et al. (2015). *GigaScience*, 4(1). https://doi.org/10.1186/s13742-015-0047-8

---

## Contact

**Junjian Wang**  
Department of Animal Science, North Carolina State University  
📧 jwang259@ncsu.edu  
🔗 https://github.com/JJWang259

For bug reports or questions, please open an [Issue](https://github.com/JJWang259/sv-association-pipeline/issues).
