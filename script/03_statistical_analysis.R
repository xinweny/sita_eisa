#### Packages ####
library(edgeR)
library(eisaR)
library(ggplot2)
library(dplyr)

setwd("/Users/Pomato/mrc/project/sita")

#### Load count files ####
exonCount <- data.matrix(read.table("processed/ExonicCounts_THAP_HEK.txt",
                                    header=TRUE,
                                    sep=" ",
                                    row.names=1))
intronCount <- data.matrix(read.table("processed/IntronicCounts_THAP_HEK.txt",
                                      header=TRUE,
                                      sep=" ",
                                      row.names=1))

#### Normalisation ####
# Remove "width" column
exon <- exonCount[, colnames(exonCount) != "width"]
intron <- intronCount[, colnames(intronCount) != "width"]

# See fraction of introns
all <- exon + intron
fracInt <- colSums(intron) / colSums(all)
summary(fracInt)

# Scale counts to the mean library size, separately for exons and introns
Nex <- t(t(exon) / colSums(exon) * mean(colSums(exon)))
Nin <- t(t(intron) / colSums(intron) * mean(colSums(intron)))

# Log transform (add a pseudocount of 8)
NLex <- log2(Nex + 8)
NLin <- log2(Nin + 8)

#### Filtering ####
# Identifying and filtering genes with low counts in introns or exons
# length(rownames(exon))
# quantGenes <- rownames(exon)[rowMeans(NLex) > 5.0 & rowMeans(NLin) > 5.0]
# length(quantGenes)

#### Calculating changes in exon/intron counts between conditions ####
# Need to generalise here for sample name format:
Dex <- NLex[, c("THAP_RNAseq_rep1",
                "THAP_RNAseq_rep2",
                "THAP_RNAseq_rep3",
                "THAP_RNAseq_rep4")] - NLex[, c("DMSO_RNAseq_rep1",
                                                "DMSO_RNAseq_rep2",
                                                "DMSO_RNAseq_rep3",
                                                "DMSO_RNAseq_rep4")]
Din <- NLin[, c("THAP_RNAseq_rep1",
                "THAP_RNAseq_rep2",
                "THAP_RNAseq_rep3",
                "THAP_RNAseq_rep4")] - NLex[, c("DMSO_RNAseq_rep1",
                                                "DMSO_RNAseq_rep2",
                                                "DMSO_RNAseq_rep3",
                                                "DMSO_RNAseq_rep4")]
Dex.Din <- Dex - Din # changes in counts attributable to post-transcriptional regulation

# Apply filter on deltas
cor(Dex[quantGenes, 1], Dex[quantGenes, 2])
cor(Din[quantGenes, 1], Din[quantGenes, 2])
cor(Dex.Din[quantGenes,1], Dex.Din[quantGenes,2])

#### Statistical analysis with edgeR framework ####
# Convert count tables to DGEList
counts <- data.frame(Ex=exon, In=intron) # merge intron and exon count tables
y <- DGEList(counts=counts,
             genes=data.frame(ENTREZID=rownames(counts)))

# Filter and normalise
# y <- y[quantGenes, ]
y <- calcNormFactors(y) # normalises library sizes by finding a set of scaling factors for the library sizes that minimizes log-fold changes between samples for most genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) # recalcNormFactAfterFilt=TRUE

# Design matrix with interaction term
smpls <- rep(factor(c("001", "002", "003", "004", "005", "006", "007", "008")), 2)
c <- factor(c(rep("1.ex", 4), rep("2.ex", 4), rep("1.in", 4), rep("2.in", 4)),
            levels=c("1.in", "2.in", "1.ex", "2.ex"))
design <- model.matrix(~ smpls + c)
design <- design[, -9]
rownames(design) <- colnames(counts)

# Log transformation (with pseudocount 2)
# region <- factor(c(rep("ex", 8), rep("in", 8)), levels=c("in", "ex"))
# cond <- rep(factor(c(rep("DMSO", 4), rep("THAP", 4))), 2)
design2 <- model.matrix(~ c)
rownames(design2) <- colnames(counts)
logFC <- predFC(y, design2, prior.count=2) # use simpler model matrix

# Estimate model parameters and fit linear model
y <- estimateDisp(y, design) # with replicates
fit <- glmQLFit(y, design)

# Quasi-likelihood F test
qlf <- glmQLFTest(fit, coef=2)

# Results table
tt <- topTags(qlf, n=nrow(y), sort.by="none") # extracts the most differentially expressed genes from a test object
head(tt$table[order(tt$table$FDR, decreasing=FALSE), ])

#### Result Visualisation ####
sig <- tt$table$FDR < 0.05
sum(sig)

# MA plot
ggplot(data=tt$table, aes(x=logCPM, y=logFC)) +
  geom_point(color="lightgrey") +
  geom_point(data=tt$table %>% dplyr::filter(FDR < 0.05),
             color="red") +
  theme_bw()

# Volcano plot


# Delta I vs. Delta E

