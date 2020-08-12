#### Packages ####
suppressPackageStartupMessages({
  library(edgeR)
  library(eisaR)
  library(ggplot2)
  library(dplyr)
  library(argparser)
  library(glue)
})

setwd("/Users/Pomato/mrc/project/sita/");

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
gse <- "GSE135753"
conditions <- c("NT", "LPS") # control vs. treatment

exon <- read.table(glue("./processed/ExonicCounts_{gse}.txt"), header=TRUE, sep="\t", row.names=1) # args$e
intron <- read.table(glue("./processed/IntronicCounts_{gse}.txt"), header=TRUE, sep="\t", row.names=1) #args$i

exon <- exon %>% select(contains(conditions))
intron <- intron %>% select(contains(conditions))

exon <- exon %>% select(contains("_rep")) %>% as.matrix()
intron <- intron %>% select(contains("_rep")) %>% as.matrix()

#### Run EISA ####
# Filter for genes which have both exons and introns counted
shared <- intersect(rownames(exon), rownames(intron))
exonsh <- exon[shared, ]
intronsh <- intron[shared, ]

message("No. of genes with ≥ 1 exon and intron: ", nrow(exonsh))

# Checks
allsh <- exonsh + intronsh
fracIn <- colSums(intronsh)/colSums(allsh)
summary(fracIn)

# Format conditions for each sample
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

#### Visualisation ####
# MA plot
MAplot <- ggplot(res_eisar$tab.ExIn, aes(x=logCPM, y=logFC)) +
          geom_point(color="lightgrey") +
          geom_point(data=res_eisar$tab.ExIn %>% filter(FDR < 0.05),
                     color="red") +
          labs(title=glue("{conditions[1]} vs. {conditions[2]}")) +
          theme(plot.title=element_text(size=20)) +
          theme_bw()

ngenes <- nrow(res_eisar$tab.ExIn %>% filter(FDR < 0.05))
message("No. of significant DE genes (FDR < 0.05): ", ngenes)

# Save output
png(glue("./processed/eisaMAplot_{gse}.png"))
print(MAplot)
dev.off()

deGenes <- res_eisar$tab.ExIn %>% arrange(FDR)
write.table(deGenes, file=glue("./processed/DEstats_{gse}.txt"), sep="\t", row.names=TRUE, col.names=TRUE)
