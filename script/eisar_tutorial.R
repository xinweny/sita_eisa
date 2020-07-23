# Tutorial: https://bioconductor.org/packages/release/bioc/vignettes/eisaR/inst/doc/eisaR.html
# Documentation: https://bioconductor.org/packages/release/bioc/manuals/eisaR/man/eisaR.pdf

#### Preparing the annotation ####
library(eisaR)

# Load TxDb object
txdbFile <- system.file("extdata", "hg19sub.sqlite", package = "eisaR")
txdb <- AnnotationDbi::loadDb(txdbFile)

# Extract filtered exonic and gene body regions
regS <- getRegionsFromTxDb(txdb = txdb, strandedData = TRUE)
regU <- getRegionsFromTxDb(txdb = txdb, strandedData = FALSE)

lengths(regS)
lengths(regU)

regS$exons

# Exporting to .gtf files
library(rtracklayer)

export(regS$exons, "hg19sub_exons_stranded.gtf")
export(regS$genebodies, "hg19sub_genebodies_stranded.gtf")

#### Quantify RNA-seq alignments in exons and introns
library(QuasR) # QuasR package for indexing and aligning short reads

# Copy sample data from package into current directory
file.copy(system.file(package = "QuasR", "extdata"), ".", recursive = TRUE)

# Align reads to a genome
proj <- qAlign(sampleFile = "extdata/samples_rna_single.txt",
               genome = "extdata/hg19sub.fa",
               aligner = "Rhisat2", splicedAlignment = TRUE)

alignmentStats(proj)