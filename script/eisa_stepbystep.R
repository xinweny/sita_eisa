#### Packages ####
library(edgeR)
library(eisaR)
library(ggplot2)
library(dplyr)
library(ExploreModelMatrix)
library(stringr)

setwd("/Users/Pomato/mrc/project/sita")

#### Load count files ####
exon <- read.table("processed/ExonicCounts_THAP_HEK_barbara.txt",
                   header=TRUE,
                   sep="\t",
                   row.names=1)
intron <- read.table("processed/IntronicCounts_THAP_HEK_barbara.txt",
                     header=TRUE,
                     sep="\t",
                     row.names=1)

# Filter for genes with at least 1 intron and exon
shared <- intersect(rownames(exon), rownames(intron))
exonsh <- exon[shared, ]
intronsh <- intron[shared, ]

# Remove extra columns
cntEx <- exonsh %>% select(contains("_rep")) %>% as.matrix()
cntIn <- intronsh %>% select(contains("_rep")) %>% as.matrix()

# Conditions
cond <- gsub("_RNAseq_rep[0-9]*.*", "", colnames(exon))
cond <- factor(cond, levels=unique(cond))

# No. of samples
nsmpls <- ncol(cntEx)

# See fraction of introns
all <- cntEx + cntIn
fracInt <- colSums(intron) / colSums(all)
summary(fracInt)

#### Normalisation ####
# Convert count tables to DGEList
cnt <- data.frame(Ex=cntEx, In=cntIn) # merge intron and exon count tables
y <- DGEList(counts=cnt,
             genes=data.frame(ENTREZID=rownames(cnt)))

# Normalise for library size
y$samples$norm.factors.exons <- rep(calcNormFactors(cntEx), 2) # normalises library sizes by finding a set of scaling factors for the library sizes that minimizes log-fold changes between samples for most genes
y$samples$norm.factors.introns <- rep(calcNormFactors(cntIn), 2)
y$samples$norm.factors.individual <- calcNormFactors(y$counts)

y$samples$lib.size.exons <- rep(colSums(cntEx), 2)
y$samples$lib.size.introns <- rep(colSums(cntIn), 2)
y$samples$lib.size.individual <- colSums(y$counts)

#### Filtering ####
# The filtering keeps genes that have count-per-million (CPM) above k in n samples, where k is determined by min.count and by the sample library sizes and n is determined by the design matrix.
quantGenes <- rownames(cntEx)[filterByExpr(y[, seq.int(nsmpls)],
                                           design=dsgn[seq.int(nsmpls), ]) &
                              filterByExpr(y[, nsmpls + seq.int(nsmpls)],
                                           design=dsgn[nsmpls + seq.int(nsmpls), ])]

pscnt <- 2

NLex <- cpm(y[, seq.int(nsmpls)], log=TRUE, prior.count=pscnt)
NLin <- cpm(y[, seq.int(nsmpls)], log=TRUE, prior.count=pscnt)

y <- y[quantGenes, ]

NLex <- NLex[quantGenes, ]
NLin <- NLin[quantGenes, ]

# Re-normalisation (recalcNormFactAfterFilt=TRUE)
y$samples$norm.factors.exons <- rep(calcNormFactors(y$counts[, seq.int(nsmpls)]), 2)
y$samples$norm.factors.introns <- rep(calcNormFactors(y$counts[, nsmpls + seq.int(nsmpls)]), 2)
y$sample$norm.factors.individual <- calcNormFactors(y$counts)

y$samples$norm.factors <- y$samples$norm.factors.exons # sizeFactor="exons"
y$samples$lib.size <- y$samples$lib.size.exons # sizeFactor="exons"

#### Design matrix and statistical model fitting ####
# Design matrix
cond2 <- rep(cond, 2L)
region <- factor(rep(c("ex", "in"), each=nsmpls),
                 levels=c("ex", "in"))
smpl <- factor(rep(sprintf("s%03d", seq.int(nsmpls)), 2)) # formatting strings

dsgn <- model.matrix(~ smpl)
c1.ex <- cond2 == levels(cond2)[1] & region == "ex"
c2.ex <- cond2 == levels(cond2)[2] & region == "ex"
dsgn <- cbind(dsgn, c1.ex, c2.ex)

rownames(dsgn) <- colnames(cnt)

# Fit statistical model
y <- estimateDisp(y, dsgn)

# Build contrast
contr <- (colnames(dsgn) == "c2.ex") - (colnames(dsgn) == "c1.ex")

# Quasi-likelihood F test
fit <- glmQLFit(y, dsgn)
tst.ExIn <- glmQLFTest(fit, contrast=contr)
tt.ExIn <- topTags(tst.ExIn, n=nrow(y), sort.by="none")

#### Calculate log-fold changes ####
lfc <- predFC(y, dsgn, prior.count=pscnt)
rownames(lfc) <- rownames(y)

# Changes in introns and/or exons between THAP vs. DMSO conditions
Dex <- rowSums(lfc[, c(3, 4)])
Din <- lfc[, 3]
Dex.Din <- lfc[, ncol(lfc)]

#### Visualisation ####
# Visualise design matrix
VisualizeDesign(
  sampleData=data.frame(
    treatment=str_extract(rownames(dsgn), "DMSO|THAP"),
    ctype=str_extract(rownames(dsgn), "Ex|In"),
    sample=gsub("Ex.|In.", "", rownames(dsgn)),
    row.names=rownames(dsgn)
  ),
  designFormula=NULL,
  designMatrix=dsgn
)

ggplot(tt.ExIn$table, aes(x=logCPM, y=logFC)) +
  geom_point(color="lightgrey") +
  geom_point(data=tt.ExIn$table %>% filter(FDR < 0.05), color="red") +
  theme_bw()

# Top 10 genes with most significant fold change
head(tt.ExIn$table %>% arrange(FDR), 10)
