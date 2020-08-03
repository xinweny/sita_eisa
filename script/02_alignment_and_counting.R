#### Packages ####
library(QuasR)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)

setwd("/rds/project/rs2099/rds-rs2099-toxgenomics/sita") # "/Users/Pomato/mrc/project/sita" "/rds/project/rs2099/rds-rs2099-toxgenomics/sita"

#### Alignment ####
proj <- qAlign(sampleFile="/home/xwy21/project/sita/list/sita_samples.txt", # "./list/sita_samples.txt" "./list/sita_sample.txt"
               genome="BSgenome.Hsapiens.UCSC.hg18",
               aligner="Rhisat2",
               splicedAlignment=TRUE,
               cacheDir="./cache")

#### Counting exons and introns ####
# TBD: make UI to choose stranded vs. unstranded
regStranded <- getRegionsFromTxDb(txdb=TxDb.Hsapiens.UCSC.hg18.knownGene, strandedData=TRUE)
# regUnstranded <- getRegionsFromTxDb(txdb=TxDb.Hsapiens.UCSC.hg18.knownGene, strandedData=FALSE)

exonCount <- qCount(proj, regStranded$exons, orientation="any", reportLevel="gene")
genebodyCount <- qCount(proj, regStranded$genebodies, orientation="any", reportLevel="gene")
intronCount <- genebodyCount - exonCount

# Convert coordinates to gene name - default given in Entrez ID
# TBD
# ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# entrezGenes <- rownames(exonCount)
# ensemblGenes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id", "entrezgene_id"), values=entrezGenes, mart=ensembl)

# Save counts to output
write.table(exonCount, file="processed/ExonicCounts_THAP_HEK.txt", row.names=TRUE, col.names=TRUE)
write.table(intronCount, file="processed/IntronicCounts_THAP_HEK.txt", row.names=TRUE, col.names=TRUE)
