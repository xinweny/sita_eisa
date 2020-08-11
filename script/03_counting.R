#### Packages ####
suppressPackageStartupMessages({
  library(QuasR)
  library(eisaR)
  library(GenomicFeatures)
  library(argparser)
  library(AnnotationDbi)
})

#### Parser ####
p <- arg_parser("QuasR counting for EISA")

p <- add_argument(p, "-i",
                  help="path to sample file")
p <- add_argument(p, "--stranded", flag=TRUE,
                  help="labels RNAseq data as stranded")
p <- add_argument(p, "-g",
                  help="path to genome file (.fa, .fasta)")
p <- add_argument(p, "-a",
                  help="path to annotation file (txdb sqlite format)")

#### Load alignment ####
proj <- qAlign(sampleFile=args$i,
               genome=args$g,
               aligner="Rhisat2",
               splicedAlignment=TRUE,
               alignmentsDir="./bam",
               cacheDir="./cache")

#### Counting exons and introns ####
# Load TxDb
txdb <- loadDb(file=args$a)

# getRegionsFromTxDb also filters out genes with only 1 exon, have exons on > 1 chromosome/both strands and overlapping genes
if (args$stranded == TRUE) {
  regions <- getRegionsFromTxDb(txdb=txdb, strandedData=TRUE)
  orient <- "same"
} else {
  regions <- getRegionsFromTxDb(txdb=txdb, strandedData=FALSE)
  orient <- "any"
}

# Counting
exonCount <- qCount(proj, regions$exons, orientation=orient, reportLevel="gene")
genebodyCount <- qCount(proj, regions$genebodies, orientation=orient, reportLevel="gene")
intronCount <- genebodyCount - exonCount

# Save counts to output
write.table(exonCount, file="processed/ExonicCounts.txt", row.names=TRUE, col.names=TRUE, sep="\t")
write.table(intronCount, file="processed/IntronicCounts.txt", row.names=TRUE, col.names=TRUE, sep="\t")