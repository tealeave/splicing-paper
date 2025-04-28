# RNA Splicing Analysis Pipeline

**Author:** David Lin  
**Date:** April 28, 2025

## Overview

This repository contains a comprehensive RNA splicing analysis pipeline for studying differential gene expression and alternative splicing events in RNA-seq data. The pipeline integrates DESeq2 for differential expression analysis and rMATS/MASER for alternative splicing analysis.

## Data Description

The project analyzes RNA-seq data from two cell lines (MB468 and R8) at different time points:

- **MB468 cell line**: Time points 0, 30, 120, and 720 minutes
- **R8 cell line**: Time points 0 and 12 hours

### Input Data

- **Count Data** (`Counts/` directory):
  - `MBsimple_counts.txt`: Gene expression counts for MB468 cell line
  - `R8simple_counts.txt`: Gene expression counts for R8 cell line
  - `ALLsimple_counts.txt`: Combined gene expression counts

- **rMATS Results** (`RMATS_results/` directory):
  - Alternative splicing events detected by rMATS for different comparisons:
    - `468_0_v_468_12`: MB468 at 0h vs 12h
    - `468_0_v_R8_0`: MB468 vs R8 at 0h
    - `468_12_v_R8_12`: MB468 vs R8 at 12h
    - `R8_0_v_R8_12`: R8 at 0h vs 12h
  - Each comparison directory contains results for different splicing event types:
    - Skipped Exon (SE)
    - Retained Intron (RI)
    - Alternative 3' Splice Site (A3SS)
    - Alternative 5' Splice Site (A5SS)
    - Mutually Exclusive Exon (MXE)

## Pipeline Workflow

The analysis pipeline consists of three main components:

1. **Differential Expression Analysis** (DESeq2)
2. **Alternative Splicing Analysis** (rMATS/MASER)
3. **Integrated Reporting**

### Pipeline Structure

```
├── config/              # Configuration parameters
├── R/                   # Reusable R functions
├── scripts/             # Analysis scripts
├── data/                # Processed data files
├── output/              # Analysis results
│   ├── figures/         # Plots and visualizations
│   └── tables/          # Result tables
└── reports/             # R Markdown reports
```

### Running the Pipeline

The entire analysis can be executed with a single command:

```bash
Rscript scripts/run_pipeline.R
```

This will:
1. Run DESeq2 differential expression analysis
2. Run rMATS/MASER alternative splicing analysis
3. Generate an integrated HTML report
4. Log all steps to a timestamped log file in `output/logs/`

### Individual Analysis Steps

If you prefer to run individual components of the pipeline:

```bash
# Run only differential expression analysis
Rscript scripts/run_deseq2_analysis.R

# Run only alternative splicing analysis
Rscript scripts/run_rmats_analysis.R
```

## Analysis Components

### 1. Differential Expression Analysis (DESeq2)

The DESeq2 analysis:
- Normalizes RNA-seq count data
- Identifies differentially expressed genes between conditions
- Generates MA plots, PCA plots, and sample distance heatmaps
- Outputs tables of differentially expressed genes with statistics

Key outputs:
- Differential expression tables (`output/tables/`)
- Visualization plots (`output/figures/`)

### 2. Alternative Splicing Analysis (rMATS/MASER)

The rMATS/MASER analysis:
- Processes rMATS output files
- Identifies significant alternative splicing events
- Generates volcano plots for each splicing event type
- Outputs tables of significant splicing events

Key outputs:
- Splicing event tables (`output/tables/splicing/`)
- Volcano plots (`output/figures/splicing/`)

### 3. Integrated Analysis

The integrated analysis:
- Combines results from differential expression and alternative splicing analyses
- Identifies genes affected by both mechanisms
- Generates comprehensive HTML report with interactive tables and plots

## Dependencies

This pipeline uses `renv` for package management to ensure reproducibility. Required R packages include:

- DESeq2
- maser
- ggplot2
- pheatmap
- dplyr
- tidyr
- targets
- tarchetypes
- rmarkdown

To install dependencies:

```r
# Restore the project environment
renv::restore()
```

## Directory Structure

- `config/`: Configuration parameters
- `R/`: Reusable R functions
  - `deseq2_functions.R`: Functions for DESeq2 analysis
  - `rmats_functions.R`: Functions for rMATS/MASER analysis
- `scripts/`: Analysis scripts
  - `run_deseq2_analysis.R`: DESeq2 differential expression analysis
  - `run_rmats_analysis.R`: rMATS/MASER alternative splicing analysis
  - `run_pipeline.R`: Master script to run the entire pipeline
- `Counts/`: RNA-seq count data
- `RMATS_results/`: rMATS output files
- `output/`: Analysis results
  - `figures/`: Plots and visualizations
  - `tables/`: Result tables
  - `logs/`: Pipeline execution logs
- `reports/`: R Markdown reports
  - `integrated_report.Rmd`: Template for integrated analysis report

## License

This project is licensed under the terms of the LICENSE file included in this repository.