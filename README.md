# sv-association-pipeline

An end-to-end, script-based pipeline for **structural variant (SV) association analysis**, designed to reproduce paper results and to serve as a general template for SV studies across species. The pipeline covers three major steps: genome-wide association studies (GWAS), statistical fine-mapping, and functional enrichment analysis of imputed SVs combined with SNPs.

## Table of Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input Data](#input-data)
- [Pipeline Steps](#pipeline-steps)
  - [Step 1: GWAS with SLEMM](#step-1-gwas-with-slemm)
  - [Step 2: Fine-Mapping with BFMAP](#step-2-fine-mapping-with-bfmap)
  - [Step 3: Functional Enrichment with GEMRICH](#step-3-functional-enrichment-with-gemrich)
- [Output Files](#output-files)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

This pipeline was developed to perform large-scale SV-trait association analyses in Holstein dairy cattle, integrating imputed SVs (from [HolPIP](https://github.com/JJWang259/sv-association-pipeline)) with high-density SNP data. It is designed to run on a local workstation or server and requires no job scheduler.

```
Imputed SVs + SNPs
        │
        ▼
  ┌─────────────┐
  │  Step 1     │  GWAS via SLEMM (linear mixed model)
  │  GWAS       │  → genome-wide SV-trait association statistics
  └──────┬──────┘
         │
         ▼
  ┌─────────────┐
  │  Step 2     │  Fine-mapping via BFMAP
  │  Fine-map   │  → posterior inclusion probabilities (PIP) per variant
  └──────┬──────┘
         │
         ▼
  ┌─────────────┐
  │  Step 3     │  Functional enrichment via GEMRICH
  │  Enrichment │  → enrichment of fine-mapped effects across annotations
  └─────────────┘
```

---

## Dependencies

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| [SLEMM](https://github.com/jiang18/slemm) | ≥ 0.93 | Linear mixed model GWAS |
| [BFMAP](https://github.com/jiang18/bfmap) | ≥ 0.65 | Bayesian fine-mapping |
| [GEMRICH](https://github.com/jiang18/gemrich) | latest | Functional enrichment analysis |
| R | ≥ 4.1 | Data processing and visualization |
| Bash | ≥ 4.0 | Script orchestration |


The pipeline includes:
1. **GWAS**: Association testing using SLEMM-GWA
2. **Fine-mapping**: Bayesian fine-mapping with BFMAP
3. **Enrichment Analysis**: 
   - MPH: Partitioned heritability analysis
   - GEMRICH: Functional enrichment of large-effect variants
---
## Installation

### Prerequisites
- R (>= 4.0)
- SLEMM/EMMAX
- BFMAP (v0.65)
- MPH (v0.52.1)
- GEMRICH

