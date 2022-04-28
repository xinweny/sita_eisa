# SITA Analysis Pipeline with EISA
A pipeline for the pairwise, comparative visualisation and analysis of stress-induced transcriptional attenuation (SITA) using exon-intron split analysis (EISA) of selected experimental conditions. This pipeline is written in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and has scripts in Python and R, and is suitable for use with HPC systems.

![Workflow of EISA RNA-seq pipeline](SITA EISA Pipeline Workflow.png?raw=true "Workflow of EISA RNA-seq pipeline")

## Background
SITA is a rapidly-induced process in cells upon stress induction, causing global transcriptional changes where the proportion of genes downregulated is far greater than those upregulated during this process.

<<<<<<< HEAD
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
=======
## Requirements
1. Snakemake
2. Python (≥ v3.x) and packages:
    - pandas (≥ v0.25.3)
    - bs4 (≥ v4.8.2)
3. R (≥ v4.x.x) and packages:
    - edgeR
    - eisaR
    - DESeq2
    - biomaRt
    - argparser
    - ggplot2
    - dplyr
    - glue
3. [Aspera Command Line Interface or Aspera Connect](https://www.ibm.com/products/aspera/downloads) added to PATH
4. [RSeQC](http://rseqc.sourceforge.net/) `infer_experiment.py` added to PATH

## Setup
### Directory structure
As the pipeline is designed to process all selected samples for a single GSE with each run, it is currently necessary to follow the directory architecture below. Each GSE is given a separate folder, which contains the `SraRunTable.txt` metadata table. The subfolder `fastq` is automatically generated if `01_download_fastqgz.py` is used, and `qc`, `bam`, `cache` and `processed` are automatically generated upon running the Snakemake pipeline.

TBD - make diagram

## Config files
The `config.yaml` stores the paths to the reference files of the respective organisms. By modifying the `config.yaml` file, this analysis can be done for any organism provided that the necessary reference files are available. The organism name (key) should match the scientific name as given in the `SraRunTable.txt` file under the `Organism` column.

### Reference files
1. FASTA reference genome
2. GTF annotation file
3. SQLITE TxDb file (easily generated from the GTF file using this tutorial [here](https://seandavi.github.io/ITR/transcriptdb.html))
4. BED genome file (downloadable from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables))

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
This pipeline takes a single input - the `SraRunTable.txt` metadata table for a particular GSE obtained from the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/).

TBD - image showing how to download/select particular runs

### Download fq.gz files
**(OPTIONAL)** The pipeline includes a script (`01_download_fastqgz.py`) for the download of compressed FASTQ files from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) for a particular GSE.

TBD: help

### Running the Snakemake pipeline
TBD snakemake command

#### 1. Alignment and counting

#### 2. Quality control
A `conda` environment containing the FastQC and MultiQC packages is provided as `qc.yml` in the `env` directory.

### EISA and DESeq2 analysis
Pairwise comparison only. TBC
