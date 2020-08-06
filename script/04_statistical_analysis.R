#### Packages ####
library(stringr)
library(glue)

#### Load data ####
setwd("/Users/Pomato/mrc/project/sita")

barbaraRes <- read.table("./processed/DEstats_barbara.txt", header=TRUE, sep="\t", row.names=1)
currentRes <- read.table("./processed/DEstats_mainhg38.txt", header=TRUE, sep="\t", row.names=1) 

#### Checks ####
barbaraGenes <- rownames(barbaraRes)
currentGenes <- rownames(currentRes)

barbaraGenes <- gsub("\\.[0-9]+", "", barbaraGenes)
currentGenes <- gsub("\\.[0-9]+", "", currentGenes)

overlaps <- intersect(barbaraGenes, currentGenes)
nonoverlaps <- setdiff(barbaraGenes, currentGenes)

length(barbaraGenes[!(barbaraGenes %in% currentGenes)])

glue("No. of overlapping genes: {length(overlaps)}")
glue("No. of non-overlapping genes: {length(nonoverlaps)}")

#### Downstream analysis ####