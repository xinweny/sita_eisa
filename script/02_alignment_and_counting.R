#### Packages ####
suppressPackageStartupMessages({
  library(QuasR)
  library(eisaR)
  library(GenomicFeatures)
  library(argparser)
  library(AnnotationDbi)
})

#### Parser ####
p <- arg_parser("QuasR alignment and counting for EISA")

p <- add_argument(p, "-i",
                  help="path to sample file")
p <- add_argument(p, "--stranded", flag=TRUE,
                  help="labels RNAseq data as stranded")
p <- add_argument(p, "-g",
                  help="path to genome file (.fa, .fasta)")
p <- add_argument(p, "-a",
                  help="path to annotation file (txdb sqlite format)")

args <- parse_args(p)

# Check if arguments are provided
stopifnot(!is.null(args$i))

#### Alignment ####
proj <- qAlign(sampleFile=args$i,
               genome=args$g,
               aligner="Rhisat2",
               splicedAlignment=TRUE,
               alignmentsDir="./bam",
               cacheDir="./cache")

# Output alignment stats for each bam file
write.table(alignmentStats(proj), file="bam/AlignmentStats.txt", row.names=TRUE, col.names=TRUE)

#### Counting exons and introns ####
# Load TxDb
txdb <- loadDb(file=args$a)

# Select chromosomes
# chroms <- c(1:22, "X", "Y")
# chroms <- paste("chr", chroms, sep="") # add chr prefix

# seqlevels(txdb) <- chroms

if (args$stranded == TRUE) {
  regions <- getRegionsFromTxDb(txdb=txdb, strandedData=TRUE)
} else {
  regions <- getRegionsFromTxDb(txdb=txdb, strandedData=FALSE)
}

# Counting
exonCount <- qCount(proj, regions$exons, orientation="any", reportLevel="gene")
genebodyCount <- qCount(proj, regions$genebodies, orientation="any", reportLevel="gene")
intronCount <- genebodyCount - exonCount

# Save counts to output
write.table(exonCount, file="processed/ExonicCounts.txt", row.names=TRUE, col.names=TRUE, sep="\t")
write.table(intronCount, file="processed/IntronicCounts.txt", row.names=TRUE, col.names=TRUE, sep="\t")
