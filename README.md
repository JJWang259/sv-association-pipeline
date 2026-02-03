# sv-association-pipeline

An end-to-end, script-based pipeline for **structural variant (SV) association analysis**, designed to reproduce paper results and to serve as a general template for SV studies across species. The workflow supports:

## Overview

This pipeline implements the statistical analysis workflow from:
> **Title of your paper**  
> Authors et al., Journal, Year

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

