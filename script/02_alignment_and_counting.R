#### Packages ####
suppressPackageStartupMessages({
  library(QuasR)
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(argparser)
  library(glue)
})

#### Parser ####
p <- arg_parser("QuasR alignment and counting for EISA")

p <- add_argument(p, "-i",
                  help="path to sample file")
p <- add_argument(p, "--stranded", flag=TRUE,
                  help="labels RNAseq data as stranded")

args <- parse_args(p)

# Check if arguments are provided
stopifnot(!is.null(args$i))

#### Alignment ####
proj <- qAlign(sampleFile=args$i,
               genome="BSgenome.Hsapiens.UCSC.hg18",
               aligner="Rhisat2",
               splicedAlignment=TRUE,
               cacheDir="./cache")

#### Counting exons and introns ####
if (args$stranded == TRUE) {
  regions <- getRegionsFromTxDb(txdb=TxDb.Hsapiens.UCSC.hg18.knownGene, strandedData=TRUE)
} else {
  regions <- getRegionsFromTxDb(txdb=TxDb.Hsapiens.UCSC.hg18.knownGene, strandedData=FALSE)
}

exonCount <- qCount(proj, regions$exons, orientation="any", reportLevel="gene")
genebodyCount <- qCount(proj, regions$genebodies, orientation="any", reportLevel="gene")
intronCount <- genebodyCount - exonCount

# Convert coordinates to gene name - default given in Entrez ID
# TBD
# ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# entrezGenes <- rownames(exonCount)
# ensemblGenes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id", "entrezgene_id"), values=entrezGenes, mart=ensembl)

# Save counts to output
date <- Sys.time()

write.table(exonCount, file=glue::glue("processed/{date}_ExonicCounts.txt"), row.names=TRUE, col.names=TRUE)
write.table(intronCount, file=glue::glue("processed/{date}_IntronicCounts.txt"), row.names=TRUE, col.names=TRUE)
