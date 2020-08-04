#### Packages ####
suppressPackageStartupMessages({
  library(edgeR)
  library(eisaR)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(argparser)
  library(ExploreModelMatrix) 
})

#### Parser ####
p <- arg_parser("Exon-Intron Split Analysis (EISA)")

p <- add_argument(p, "-e", help="path to exon counts file")
p <- add_argument(p, "-i", help="path to intron counts file")

args <- parse_args(p)

# Check if arguments are provided
stopifnot(!is.null(args$e),
          !is.null(args$i))

#### Load counts files ####
exon <- read.table(args$e, header=TRUE, sep=" ", row.names=1)
intron <- read.table(args$i, header=TRUE, sep=" ", row.names=1)

exon <- exon %>% select(contains("_rep")) %>% as.matrix()
intron <- intron %>% select(contains("_rep")) %>% as.matrix()

#### Run EISA ####
# Filter for genes which have both exons and introns counted
shared <- intersect(rownames(exon), rownames(intron))
exonsh <- exon[shared, ]
intronsh <- intron[shared, ]

# Format conditions for each sample
cond <- gsub("_RNAseq_rep[0-9]*.*", "", colnames(exonsh))

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

#### Visualisation ####
# MA plot
MAplot <- ggplot(res_eisar$tab.ExIn, aes(x=logCPM, y=logFC)) +
          geom_point(color="lightgrey") +
          geom_point(data=res_eisar$tab.ExIn %>% filter(FDR < 0.05),
                     color="red") +
          theme_bw()

# Save output
png("eisaMAplot.png")
print(MAplot)
dev.off()

head(res_eisar$tab.ExIn %>% arrange(FDR))
