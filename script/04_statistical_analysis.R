#### Packages ####
library(stringr)
library(glue)
library(dplyr)

#### Load data ####
setwd("/Users/Pomato/mrc/project/sita")

barbaraRes <- read.table("./processed/DEstats_barbara.txt", header=TRUE, sep="\t", row.names=1)
currentRes <- read.table("./processed/DEstats_xwy.txt", header=TRUE, sep="\t", row.names=1)

#### Checks ####
barbaraSigRes <- barbaraRes %>% filter(FDR < 0.05)
currentSigRes <- currentRes %>% filter(FDR < 0.05)

barbaraGenes <- rownames(barbaraSigRes)
currentGenes <- rownames(currentSigRes)

barbaraGenes <- gsub("\\.[0-9]+", "", barbaraGenes)
currentGenes <- gsub("\\.[0-9]+", "", currentGenes)

overlaps <- intersect(barbaraGenes, currentGenes)
nonoverlaps <- setdiff(barbaraGenes, currentGenes)

glue("No. of overlapping sig. genes: {length(overlaps)}")
glue("No. of non-overlapping sig. genes: {length(nonoverlaps)}")

#### Downstream analysis ####