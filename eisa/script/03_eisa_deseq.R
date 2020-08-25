#### Packages ####
suppressPackageStartupMessages({
  library(edgeR)
  library(eisaR)
  library(ggplot2)
  library(dplyr)
  library(argparser)
  library(glue)
  library(DESeq2)
})

setwd("/Users/Pomato/mrc/project/sita/")

#### Parser ####
# p <- arg_parser("Exon-Intron Split Analysis (EISA)")
# 
# p <- add_argument(p, "-e", help="path to exon counts file")
# p <- add_argument(p, "-i", help="path to intron counts file")
# 
# args <- parse_args(p)
# 
# # Check if arguments are provided
# stopifnot(!is.null(args$e),
#           !is.null(args$i))

#### Load counts files ####
gse <- "GSE103719"
conditions <- c("DMSO", "THAP") # control vs. treatment

exon <- read.table(glue("./processed/ExonicCounts_{gse}.txt"), header=TRUE, sep="\t", row.names=1) # args$e
intron <- read.table(glue("./processed/IntronicCounts_{gse}.txt"), header=TRUE, sep="\t", row.names=1) #args$i

exon <- exon %>% select(contains(conditions)) %>% as.matrix()
intron <- intron %>% select(contains(conditions)) %>% as.matrix()

#### EISA ####
# Filter for genes which have both exons and introns counted
shared <- intersect(rownames(exon), rownames(intron))
exonsh <- exon[shared, ]
intronsh <- intron[shared, ]

message("No. of genes with ≥ 1 exon and intron: ", nrow(exonsh))

# Checks
allsh <- exonsh + intronsh
fracIn <- colSums(intronsh)/colSums(allsh)
summary(fracIn)

# Format and select conditions for each sample
cond <- gsub("_rep[0-9]*.*", "", colnames(exonsh))
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

# Flip logFC
res_eisar$tab.ExIn$logFC <- -(res_eisar$tab.ExIn$logFC)

# EISA MA plot
MAplot <- ggplot(res_eisar$tab.ExIn, aes(x=logCPM, y=logFC)) +
  geom_point(color="lightgrey") +
  geom_point(data=res_eisar$tab.ExIn %>% filter(FDR < 0.05),
             color="red") +
  labs(title=glue("{conditions[1]} vs. {conditions[2]}")) +
  theme_bw() +
  theme(plot.title=element_text(size=15, face="bold"))

# Save output
png(glue("./processed/eisaMAplot_{gse}_{conditions[1]}.{conditions[2]}.png"))
print(MAplot)
dev.off()

#### DESeq2 ####
# Make column metadata
cond <- gsub("_rep[0-9]*.*", "", colnames(exon))
colData <- data.frame(row.names=colnames(exon),
                      condition=factor(cond, levels=conditions))

dataset <- DESeqDataSetFromMatrix(countData=exon,
                                  colData=colData,
                                  design= ~ condition)

# Run DESeq
dds <- DESeq(dataset)
result <- results(dds, contrast=c('condition', conditions), alpha=0.05)
result <- result[complete.cases(result), ] # remove any rows with NA
result <- result[order(result$padj), ] # order by adjusted p-value

message("No. of significant DE genes: ", sum(result$padj < 0.05, na.rm=TRUE))

# Save DESeq results table to output
write.table(result,
            file=glue("processed/DESeq_{gse}_{conditions[1]}.{conditions[2]}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t")

# DESeq2 MA plot
DESeq2::plotMA(result, main=glue("{conditions[1]} vs. {conditions[2]}"))
