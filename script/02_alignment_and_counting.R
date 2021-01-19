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

<<<<<<< HEAD
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
=======
p <- add_argument(p, "-m",
                  help="path to SraRunTable.txt metadata file")
p <- add_argument(p, "-t",
                  help="no. of threads")

args <- parse_args(p)

#### Load data ####
gse <- tail(strsplit(getwd(), "/")[[1]], n=1)

# Make SampleFile
sampleFiles <- Sys.glob("SampleFile*.txt")

if (length(sampleFiles) == 0) {
  message("Creating QuasR sample file(s)...")
  
  system(glue("python3 /home/xwy21/project/sita_eisa/script/generate_samplefile.py -m {args$m} -e fq.gz"))
  sampleFiles <- Sys.glob("SampleFile*.txt")
} else {
  message("QuasR sample file(s) already created.")
}

metadata <- read.table(args$m, header=TRUE, sep=",")
organism <- metadata$Organism[1]
message(glue("Organism: {organism}"))

if (organism == "Homo sapiens") {
  genomeFile <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/human/GRCh38.primary_assembly.genome.fa"
  txdb <- loadDb(file="/rds/project/rs2099/rds-rs2099-toxgenomics/shared/human/txdb.hg38.gencode.v34.sqlite")
  bedGenome <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/human/hg38_Gencode_V24.bed"
} else if (organism == "Mus musculus") {
  genomeFile <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/mouse/GRCm38.primary_assembly.genome.fa"
  txdb <- loadDb(file="/rds/project/rs2099/rds-rs2099-toxgenomics/shared/mouse/txdb.gencode.vM25.annotation.sqlite")
  bedGenome <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/mouse/mm10_Gencode_VM18.bed"
}

#### Make cluster ####
cl <- makeCluster(as.numeric(args$t))

exonCountList <- list()
intronCountList <- list()

for (sampleFile in sampleFiles) {
  #### Alignment ####
  proj <- qAlign(sampleFile=sampleFile,
                 genome=genomeFile,
                 aligner="Rhisat2",
                 splicedAlignment=TRUE,
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

  percRvStrand <- as.numeric(tail(strsplit(rseqcOutput[6], " ")[[1]], n=1))
  percFwStrand <- as.numeric(tail(strsplit(rseqcOutput[5], " ")[[1]], n=1))
  
  stranded_fw <- percFwStrand > strandThresh
  stranded_rv <- percRvStrand > strandThresh
  
  # if (libraryLayout == "PairEnd") {
  #   percFwStrand <- as.numeric(tail(strsplit(rseqcOutput[5], " ")[[1]], n=1))
  #   # percRvStrand <- as.numeric(tail(strsplit(rseqcOutput[6], " ")[[1]], n=1))
  #   stranded <- percFwStrand > strandThresh # | percRvStrand > strandThresh
  # } else {
  #   percFwStrand <- as.numeric(tail(strsplit(rseqcOutput[5], " ")[[1]], n=1))
  #   percRvStrand <- as.numeric(tail(strsplit(rseqcOutput[6], " ")[[1]], n=1))
  #   stranded <- percFwStrand > strandThresh | percRvStrand > strandThresh
  # }
  
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
  
  if (length(sampleFiles) > 1) {
    exonCountList[[length(exonCountList) + 1]] <- exonCount
    intronCountList[[length(intronCountList) + 1]] <- intronCount
  }
}

#### Save to output ####
if (length(sampleFiles) > 1) {
  exonCount <- merge(exonCountList[[1]], exonCountList[[2]], by=0, all=TRUE)
  intronCount <- merge(intronCountList[[1]], intronCountList[[2]], by=0, all=TRUE)

  exonCount[is.na(exonCount)] <- 0
  intronCount[is.na(intronCount)] <- 0
}

>>>>>>> 659a09d652073412978158006cc0255ee77b6d5c
write.table(exonCount, file=glue("processed/{gse}_ExonicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(intronCount, file=glue("processed/{gse}_IntronicCounts.txt"), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# Delete cache
system("rm -r cache/*")
