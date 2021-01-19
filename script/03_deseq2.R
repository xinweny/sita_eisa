#### Packages ####
suppressPackageStartupMessages({
  library(DESeq2)
  library(glue)
  library(dplyr)
})

setwd("/Users/Pomato/mrc/project/sita/")

#### Load data ####
gse <- "GSE109834"
conditions <- c("PBS", "LPS")

countData <- read.table(glue("./processed/ExonicCounts_{gse}.txt"), header=TRUE, sep="\t", row.names=1)
selectedCountData <- countData %>% select(contains(conditions)) %>% as.matrix()

# Make column metadata
cond <- gsub("_rep[0-9]*.*", "", colnames(selectedCountData))
colData <- data.frame(row.names=colnames(selectedCountData),
                      condition=factor(cond, levels=conditions))

#### Run DESeq2 ####
dataset <- DESeqDataSetFromMatrix(countData=selectedCountData,
                              colData=colData,
                              design= ~ condition)

dds <- DESeq(dataset)
result <- results(dds, contrast=c('condition', conditions))
result <- result[complete.cases(result), ] # remove any rows with NA

#### Analyse results ####
# Top X DE genes

# MA plot
plotMA(result, main=glue("{conditions[1]} vs. {conditions[2]}"))