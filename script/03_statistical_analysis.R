#### Packages ####
library(edgeR)

#### Load count files ####
# TBD
exonCount <- 
intronCount <-

#### Normalisation ####
# Remove "width" column
exon <- exonCount[, colnames(exonCount) != "width"]
intron <- intronCount[, colnames(intronCount) != "width"]

all <- exon + intron
fracInt <- colSums(intron) / colSums(all)
summary(fracInt)

# Scale counts to the mean library size, separately for exons and introns
Nex <- t(t(exon) / colSums(exon) * mean(colSums(exon)))
Nin <- t(t(intron) / colSums(intron) * mean(colSums(intron)))

# Log transformation (with pseudocount 8)
NLex <- log2(Nex + 8)
NLin <- log2(Nin + 8)

#### Filtering ####
# Identifying and filtering genes with low counts in introns or exons
length(rownames(exon))
quantGenes <- rownames(exon)[rowMeans(NLex) > 5.0 & rowMeans(NLin) > 5.0]
length(quantGenes)

#### Calculating changes in exon/intron counts between conditions ####
# Need to generalise here for sample name format:
Dex <- NLex[, c("DMSO_RNAseq_rep1")] - NLex[, c("THAP_RNAseq_rep4")]
Din <- NLin[, c("DMSO_RNAseq_rep1")] - NLex[, c("THAP_RNAseq_rep4")]
Dex.Din <- Dex - Din # changes in counts attributable to post-transcriptional regulation

# Apply filter on deltas
cor(Dex[quantGenes], Din[quantGenes])

#### Statistical analysis with edgeR framework ####
counts <- data.frame(Ex=exon, In=intron)
y <- DGEList(counts=counts,
             genes=data.frame(ENTREZID=rownames(counts)))

# Filter and normalise
y <- y[quantGenes, ]
y <- calcNormFactors(y) # normalises library sizes by finding a set of scaling factors for the library sizes that minimizes log-fold changes between samples for most genes

# Design matrix with interaction term
region <- factor(c("ex", "ex", "in", "in"),
                 levels=c("in", "ex")) # need to make factor generation flexible - from sample name
cond <- rep(factor(c("DMSO", "THAP")), 2) # need to make condition generation flexible - from sample name
design <- model.matrix(~ region * cond)
rownames(design) <- colnames(counts)

# Estimate model parameters
y <- estimateDisp(y, design) # with replicates
fit <- glmFit(y, design)