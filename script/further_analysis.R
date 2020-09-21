#### Packages ####
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(glue)

setwd("~/mrc/project/sita_eisa")

#### Functions ####
format_condition <- function (colnames) {
  new_cond <- gsub("_[0-9]*$", "", colnames)
  new_cond <- gsub("_rep[0-9]*$", "", new_cond)
  new_cond <- gsub("GSM[0-9]*_", "", new_cond)
  
  return(new_cond)
}

#### Load data ####
gse <- "GSE139592"

exon <- read.table(glue("processed/{gse}_ExonicCounts.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE) %>%
  as.matrix()
intron <- read.table(glue("processed/{gse}_IntronicCounts.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE) %>%
  as.matrix()

shared <- intersect(rownames(exon), rownames(intron))
exonsh <- exon[shared, ]
intronsh <- intron[shared, ]

# Format conditions
cond <- format_condition(colnames(exon))

# Make DESeq2 object
colData <- data.frame(row.names=colnames(exon),
                      condition=factor(cond))

exonDDS <- DESeq(DESeqDataSetFromMatrix(countData=exon,
                                      colData=colData,
                                      design=~ condition))
intronDDS <- DESeq(DESeqDataSetFromMatrix(countData=exon,
                                        colData=colData,
                                        design=~ condition))

exonRld <- vst(exonDDS, blind=TRUE)
intronRld <- vst(intronDDS, blind=TRUE)

# PCA plot
pcaExon <- plotPCA(exonRld) + ggtitle("Exon counts")
pcaIntron <- plotPCA(intronRld) + ggtitle("Intron counts")

grid.arrange(pcaExon, pcaIntron,
             nrow=1)

# Heatmap
par(mfrow=c(1, 2))
pheatmap(cor(assay(exonRld)),
        main="Exon counts")
pheatmap(cor(assay(intronRld)),
        main="Intron counts")
