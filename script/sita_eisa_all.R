#### Packages ####
suppressPackageStartupMessages({
  library(eisaR)
  library(ggplot2)
  library(dplyr)
  library(glue)
  library(DESeq2)
  library(biomaRt)
})

setwd("/Users/Pomato/mrc/project/sita_eisa/processed/")

#### Functions ####
add_ensembl_symbol <- function (table) {
  genes <- row.names(table)
  
  if (grepl("ENSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "hsapiens_gene_ensembl"
    symbol <- "hgnc_symbol"
  } else if (grepl("ENSMUSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "mmusculus_gene_ensembl"
    symbol <- "mgi_symbol"
  }
  
  mart <- useDataset(ensemblDataset, useMart("ENSEMBL_MART_ENSEMBL", host="http://www.ensembl.org"))
  geneList <- getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id", symbol),
                    values=genes,
                    mart=mart)
  
  row.names(geneList) <- geneList[, 1]
  geneList[, 1] <- NULL
  
  table$geneSymbol <- geneList[, 1][match(rownames(table), rownames(geneList))]
  newTable <- table
  
  return(newTable)
}

format_condition <- function (colnames) {
  new_cond <- gsub("_[0-9]*$", "", colnames)
  new_cond <- gsub("_rep[0-9]*$", "", new_cond)
  new_cond <- gsub("GSM[0-9]*_", "", new_cond)
  
  return(new_cond)
}

run_eisa <- function (exon, intron) {
  # Format conditions
  cond <- format_condition(colnames(exon))
  
  conditions <- unique(cond)
  cond <- factor(cond, levels=rev(conditions))
  
  # Intersect exon and intron
  shared <- intersect(rownames(exon), rownames(intron))
  exonsh <- exon[shared, ]
  intronsh <- intron[shared, ]
  
  allsh <- exonsh + intronsh
  fracIn <- colSums(intronsh)/colSums(allsh)
  summary(fracIn)
  
  stopifnot(all(colnames(exonsh) == colnames(intronsh)))
  
  res_eisar <- runEISA(cntEx=exonsh, cntIn=intronsh,
                       cond=cond,
                       method=NULL,
                       modelSamples=TRUE,
                       geneSelection="filterByExpr",
                       statFramework="QLF",
                       effects="predFC",
                       pscnt=2,
                       recalcNormFactAfterFilt=TRUE,
                       recalcLibSizeAfterFilt=FALSE)
  
  message("No. of significant DE genes (FDR < 0.05): ", nrow(res_eisar$tab.ExIn %>% filter(FDR < 0.05)))
  
  # Add gene symbol
  res_eisar$tab.ExIn <- add_ensembl_symbol(res_eisar$tab.ExIn)

  # Ratio of up and down regulated significant genes
  nUp <- nrow(filter(res_eisar$tab.ExIn, FDR < 0.05 & logFC > 0))
  nDown <- nrow(filter(res_eisar$tab.ExIn, FDR < 0.05 & logFC < 0))
  ratioUpDown <- nUp / nDown
  message("Ratio of up:down regulated genes: ", ratioUpDown)
    
  # EISA MA plot
  MAplot <- ggplot(res_eisar$tab.ExIn, aes(x=logCPM, y=logFC)) +
    geom_point(color="lightgrey") +
    geom_point(data=res_eisar$tab.ExIn %>% filter(FDR < 0.05),
               color="red") +
    labs(title=glue("{gse}: {conditions[1]} vs. {conditions[2]}"),
         caption=glue("UP={nUp}, DOWN={nDown}")) +
    theme_bw() +
    theme(plot.title=element_text(size=15, face="bold"),
          plot.caption=element_text(size=15))
  
  # Save output
  png(glue("{gse}_eisaMAplot_{conditions[1]}.{conditions[2]}.png"))
  print(MAplot)
  dev.off()
  
  deGenes <- res_eisar$tab.ExIn %>% arrange(FDR, -logFC)
  write.table(deGenes, file=glue("{gse}_eisaDE_{conditions[1]}.{conditions[2]}.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
}

run_deseq2 <- function (exon, intron, paired=FALSE, alpha=0.05, lfcThresh=0) {
  # Make column metadata
  cond <- format_condition(colnames(exon))
  conditions <- unique(cond)
  
  if (paired) {
    colData <- data.frame(row.names=colnames(exon),
                          replicate=factor(gsub(".rep", "", colnames(exon)),
                                           levels=seq.int(1:(ncol(exon) / 2))),
                          condition=factor(cond,
                                           levels=conditions))
    
    dataset <- DESeqDataSetFromMatrix(countData=exon,
                                      colData=colData,
                                      design=~ replicate + condition)
  } else {
    colData <- data.frame(row.names=colnames(exon),
                          condition=factor(cond, levels=conditions))
    
    dataset <- DESeqDataSetFromMatrix(countData=exon,
                                      colData=colData,
                                      design=~ condition)
  }
  
  # Set reference level as control
  dataset$condition <- relevel(dataset$condition, ref=conditions[1])
  
  dataset <- dataset[rowSums(counts(dataset)) >= 10, ] # pre-filter genes with no reads
  
  ## Run DESeq
  dds <- DESeq(dataset)
  coef <- tail(resultsNames(dds), n=1)
  
  res <- results(dds, name=coef,
                 alpha=alpha,
                 lfcThreshold=lfcThresh,
                 altHypothesis='greaterAbs',
                 pAdjustMethod='BH')
  resLFC <- lfcShrink(dds, coef=coef, res=res)
  
  summary(res)
  message(glue("No. of significant DE genes (FDR < {alpha}): "), sum(resLFC$padj < alpha, na.rm=TRUE))
  
  # Add gene symbol
  resLFC <- add_ensembl_symbol(resLFC)
  
  # Save DESeq results table to output
  deGenes <- as.data.frame(resLFC) %>% arrange(padj, desc(log2FoldChange)) # order by adjusted p-value and FC
  write.table(deGenes,
              file=glue("{gse}_DESeq_{conditions[1]}.{conditions[2]}.txt"),
              row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
  
  # Ratio of up and down regulated significant genes
  nUp <- nrow(filter(deGenes, padj < alpha & log2FoldChange > 0))
  nDown <- nrow(filter(deGenes, padj < alpha & log2FoldChange < 0))
  ratioUpDown <- nUp / nDown
  message("Ratio of up:down regulated genes: ", ratioUpDown)
  
  ## Visualisation
  resLFC <- resLFC[order(-resLFC$padj),]
  
  # DESeq2 MA plot
  DESeq2::plotMA(resLFC, main=glue("{gse}: {conditions[1]} vs. {conditions[2]}
                                           UP={nUp}, DOWN={nDown}"))
  
  png(glue("{gse}_DESeqMAplot_{conditions[1]}.{conditions[2]}.png"))
  print(DESeq2::plotMA(resLFC, main=glue("{gse}: {conditions[1]} vs. {conditions[2]}
                                           UP={nUp}, DOWN={nDown}")))
  dev.off()
}
#### Load data ####
compList <- read.table("complist.txt",
                       header=TRUE, sep="\t", check.names=FALSE)

for (i in 1:nrow(compList)) {
  # Load GSE count table and select condition
  gse <- compList[i, "GSE"]
  selectConditions <- c(compList[i, "control"], compList[i, "condition"]) # control vs. treatment
  
  exon <- read.table(glue("{gse}_ExonicCounts.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  intron <- read.table(glue("{gse}_IntronicCounts.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  
  exon <- exon %>% dplyr::select(matches(selectConditions)) %>% as.matrix()
  intron <- intron %>% dplyr::select(matches(selectConditions)) %>% as.matrix()
  
  stopifnot(length(unique(format_condition(colnames(exon)))) == 2)
  
  # Run EISA
  # run_eisa(exon, intron)
  
  # Run DESeq2
  run_deseq2(exon, intron, alpha=0.01)
}