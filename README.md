# Computational Analysis of Colorectal Cancer

## Author

**Dhiraj Kumar Biswal**  
GitHub: `codedemon404`

## Overview

This project analyzes gene-expression data from colorectal cancer tumor and matched normal tissue samples.

The analysis uses the GSE74602 dataset from the NCBI Gene Expression Omnibus (GEO).

## Objective

To identify genes showing significant differences in expression between colorectal cancer tumor tissue and matched normal tissue.

## Dataset

- GEO accession: GSE74602
- Platform: GPL6104
- Samples analyzed: 60
- Tumor samples: 30
- Normal samples: 30
- Matched tumor-normal pairs: 30
- Expression probes analyzed: 22,184

## Analysis

The workflow includes:

1. Loading the gene-expression matrix
2. Identifying matched tumor-normal samples
3. Mapping microarray probes to gene names
4. Differential expression analysis
5. Multiple-testing correction using FDR
6. Volcano plot visualization
7. Functional/pathway enrichment analysis

## Tools

- Python
- pandas
- NumPy
- SciPy
- Matplotlib
- g:Profiler

## Important Results

The analysis identified genes showing increased and decreased expression in colorectal tumor tissue compared with matched normal tissue.

Candidate differentially expressed genes included FOXQ1, UBE2C, KIAA1199, TOP2A, SCARA5 and ABCA8.

## Data Source

Gene-expression data were obtained from the NCBI Gene Expression Omnibus (GEO).

- GEO accession: GSE74602
- Platform: GPL6104

## AI Assistance

ChatGPT was used during development for code debugging, organization, and assistance with understanding the analysis workflow.

## Copyright

Copyright © 2026 Dhiraj Kumar Biswal.

This repository is provided for academic and research reference. No license is granted for reuse, redistribution, or submission of this work as another person's work without permission.
