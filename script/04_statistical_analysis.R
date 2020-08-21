#### Packages ####
library(glue)
library(dplyr)

#### Load data ####
setwd("/Users/Pomato/mrc/project/sita")

gseCond <- list(GSE97744=c("non.polarized", "IFNg.LPS"),
                GSE103719=c("DMSO", "THAP"),
                GSE109834=c("PBS", "LPS"),
                GSE120807=c("unstimulated", "^LPS"),
                GSE135618=c("PBS_cells", "LPS_cells"),
                GSE135753=c("NT", "^LPS"))

gseSens <- list()

#### Checks ####
# Sensitivity analysis of EISA
for (gse in names(gseCond)) {
  message("Sensitivity analysis of ", gse)
  exon <- read.table(glue("processed/{gse}_ExonicCounts.txt"), header=TRUE, sep="\t", row.names=1)
  intron <- read.table(glue("processed/{gse}_IntronicCounts.txt"), header=TRUE, sep="\t", row.names=1)
  
  exon <- exon %>% dplyr::select(matches(gseCond[[gse]])) %>% as.matrix()
  intron <- intron %>% dplyr::select(matches(gseCond[[gse]])) %>% as.matrix()
  
  shared <- intersect(rownames(exon), rownames(intron))
  exonsh <- exon[shared, ]
  intronsh <- intron[shared, ]
  
  cond <- gsub("_rep[0-9]*.*", "", colnames(exonsh))
  conditions <- unique(cond)
  
  stopifnot(all(colnames(exonsh) == colnames(intronsh)))
  
  nGenes <- c()
  
  for (bool in c(FALSE, TRUE)) {
    res_eisar <- runEISA(cntEx=exonsh, cntIn=intronsh,
                         cond=cond,
                         method=NULL,
                         modelSamples=bool,
                         geneSelection="filterByExpr",
                         statFramework="QLF",
                         effects="predFC",
                         pscnt=2,
                         recalcNormFactAfterFilt=TRUE,
                         recalcLibSizeAfterFilt=FALSE)
    
    sigGenes <- nrow(res_eisar$tab.ExIn %>% filter(FDR < 0.05))
    nGenes <- append(nGenes, sigGenes)
  }
  
  diff <- nGenes[2] - nGenes[1]
  percent <- round(diff / nGenes[1] * 100, digits=2)
  statement <- glue("{diff} genes from {nGenes[1]} -> {nGenes[2]} (Δ {percent}%)")
  gseSens[[gse]] <- statement
}

#### Downstream analysis ####