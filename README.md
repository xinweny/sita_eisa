# SITA Analysis Pipeline with EISA

A pipeline for the pairwise, comparative visualisation and analysis of stress-induced transcriptional attenuation (SITA) using exon-intron split analysis (EISA) of selected experimental conditions. Currently only supports studies deposited on the [NCBI GEO database](https://www.ncbi.nlm.nih.gov/geo/) involving Human or Mouse bulk RNA-seq data. This pipeline has a CLI and is scripted using Python and R, and is suitable for use with HPC systems.

## Table of Contents
TBD

## Background
SITA is a rapidly-induced process in cells upon stress induction, causing global transcriptional changes where the proportion of genes downregulated is far greater than those upregulated during this process. TBC

## Setup
### Scripts
TBD

### Directory structure
As the pipeline is designed to process all selected SRRs/GSMs for a single GSE with each run, it is currently necessary to follow the directory architecture below. Each GSE is given a separate folder, which contains the SraRunTable.txt metadata table. The subfolders fastq, bam, cache and processed are automatically generated upon running the first step of the pipeline.

TBD - make diagram

## Requirements
1. Python (≥ v3.x) and packages:
    - pandas (≥ v0.25.3)
    - bs4 (≥ v4.8.2)
2. R (≥ v4.x.x) and packages:
    - edgeR
    - eisaR
    - DESeq2
    - biomaRt
    - argparser
    - ggplot2
    - dplyr
    - glue
3. [Aspera Command Line Interface or Aspera Connect](https://www.ibm.com/products/aspera/downloads) added to PATH
4. [RSeQC](http://rseqc.sourceforge.net/) infer_experiment.py added to PATH

## Pipeline Workflow
TBD - make diagram

## Usage
### Input
This pipeline takes a single input - the SraRunTable.txt metadata table for a particular GSE obtained from the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/).

TBD - image showing how to download/select particular runs

### 1. Download fq.gz files
Two versions of scripts, aspera faster (choice #1) otherwise prefetch --> fastq_dump. TBC
Format and automation of file renaming to reduce sample mislabeling. TBC

### 2. Alignment and counting
Automation of... organism and RNA-seq strandedness detection. Need to modify file path to required .fa and .gtf files, and BED genome file for RSeQC. TBC

### 3. EISA and DESeq2 analysis
Pairwise comparison only. TBC
