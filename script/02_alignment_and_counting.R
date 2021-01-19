#### Packages ####
suppressPackageStartupMessages({
  library(QuasR)
  library(eisaR)
  library(GenomicFeatures)
  library(argparser)
  library(AnnotationDbi)
  library(glue)
  library(dplyr)
})

#### Parser ####
p <- arg_parser("QuasR alignment and counting for EISA")

p <- add_argument(p, "-w",
                  help="work directory")
p <- add_argument(p, "-t",
                  help="no. of threads")
p <- add_argument(p, "-a",
                  help="spliced alignment (TRUE or FALSE)")
p <- add_argument(p, "-f",
                  help="path to .fa genome")
p <- add_argument(p, "-b",
                  help="path to .bed genome")
p <- add_argument(p, "-g",
                  help="path to .gtf file")
p <- add_argument(p, "-q",
                  help="path to QuasR sample file")

args <- parse_args(p)

# Set work directory
setwd(args$w)

#### Load data ####
gse <- tail(strsplit(getwd(), "/")[[1]], n=1)

genomeFile <- args$f
txdb <- makeTxDbFromGFF(args$g)
bedGenome <- args$b
sampleFile <- args$q

#### Make cluster ####
cl <- makeCluster(as.numeric(args$t))

#### Alignment ####
proj <- qAlign(sampleFile=sampleFile,
               genome=genomeFile,
               aligner="Rhisat2",
               splicedAlignment=as.logical(args$a),
               alignmentsDir="./bam",
               cacheDir="./cache",
               clObj=cl)

# Output alignment stats
write.table(alignmentStats(proj), file="bam/AlignmentStats.txt", row.names=TRUE, col.names=TRUE)

#### Counting exons and introns ####
# Check strandedness
bamFile <- Sys.glob("bam/*.bam")[1]
rseqcOutput <- system(glue("infer_experiment.py -r {bedGenome} -i '{bamFile}'"), intern=TRUE)
message(cat(rseqcOutput, sep="\n"))

libraryLayout <- strsplit(rseqcOutput[3], " ")[[1]][3]

strandThresh <- 0.9

percFwStrand <- as.numeric(tail(strsplit(rseqcOutput[5], " ")[[1]], n=1))
percRvStrand <- as.numeric(tail(strsplit(rseqcOutput[6], " ")[[1]], n=1))

stranded_fw <- percFwStrand > strandThresh
stranded_rv <- percRvStrand > strandThresh

message(glue("Stranded: {stranded_fw | stranded_rv}"))

# getRegionsFromTxDb also filters out genes with only 1 exon, have exons on > 1 chromosome/both strands and overlapping genes
if (stranded_fw == TRUE) {
  regions <- getRegionsFromTxDb(txdb=txdb, strandedData=TRUE)
  orient <- "same"
} else if (stranded_rv == TRUE) {
  regions <- getRegionsFromTxDb(txdb=txdb, strandedData=TRUE)
  orient <- "opposite"
} else {
  regions <- getRegionsFromTxDb(txdb=txdb, strandedData=FALSE)
  orient <- "any"
}

if (grepl("PAIRED", sampleFile, fixed=TRUE)) {
  read <- "first"
}  else {
  read <- "any"
}

# Counting
exonCount <- qCount(proj, regions$exons,
                    orientation=orient,
                    reportLevel="gene",
                    useRead=read,
                    clObj=cl)
genebodyCount <- qCount(proj, regions$genebodies,
                        orientation=orient,
                        reportLevel="gene",
                        useRead=read,
                        clObj=cl)

intronCount <- genebodyCount - exonCount

# Remove width column
exonCount <- exonCount[, -1]
intronCount <- intronCount[, -1]

# Remove version from genes
row.names(exonCount) <- gsub("\\.[0-9]*", "", rownames(exonCount))
row.names(intronCount) <- gsub("\\.[0-9]*", "", rownames(intronCount))

#### Save to output ####
write.table(exonCount, file=glue("processed/{gse}_ExonicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(intronCount, file=glue("processed/{gse}_IntronicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# Delete cache
system("rm -r cache/*")
