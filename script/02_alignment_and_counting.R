#### Packages ####
suppressPackageStartupMessages({
  library(QuasR)
  library(eisaR)
  library(GenomicFeatures)
  library(argparser)
  library(AnnotationDbi)
  library(stringr)
})

#### Parser ####
p <- arg_parser("QuasR alignment and counting for EISA")

p <- add_argument(p, "-i",
                  help="path to sample file")
p <- add_argument(p, "-m",
                  help="model organism (human or mouse)")

args <- parse_args(p)

stopifnot(args$m == "human" | args$m == "mouse")

# Make sure SampleFile exists
if (!file.exists("SampleFile.txt")) {
  system("python3 /home/xwy21/project/sita/script/generate_samplefile.py -e fastq")
  args$i <- "SampleFile.txt"
}

if (args$m == "human") {
  genomeFile <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/human/GRCh38.primary_assembly.genome.fa"
  txdb <- loadDb(file="/rds/project/rs2099/rds-rs2099-toxgenomics/shared/human/txdb.hg38.gencode.v34.sqlite")
  bedGenome <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/human/hg38_Gencode_V24.bed"
} else {
  genomeFile <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/mouse/GRCm38.primary_assembly.genome.fa"
  txdb <- loadDb(file="/rds/project/rs2099/rds-rs2099-toxgenomics/shared/mouse/txdb.gencode.vM25.annotation.sqlite")
  bedGenome <- "/rds/project/rs2099/rds-rs2099-toxgenomics/shared/mouse/mm10_Gencode_VM18.bed"
}

#### Alignment ####
proj <- qAlign(sampleFile=args$i,
               genome=genomeFile,
               aligner="Rhisat2",
               splicedAlignment=TRUE,
               alignmentsDir="./bam",
               cacheDir="./cache")

# Output alignment stats for each bam file
write.table(alignmentStats(proj), file="bam/AlignmentStats.txt", row.names=TRUE, col.names=TRUE)

#### Counting exons and introns ####
# Check strandedness
bamFile <- Sys.glob("*.bam")[1]
rseqcOutput <- system(glue("infer_experiment.py -r {bedGenome} -i {bamFile}"), intern=TRUE)
percSameStrand <- as.numeric(tail(strsplit(rseqcOutput[5], " ")[[1]], n=1))
stranded <- percSameStrand > 0.9

# getRegionsFromTxDb also filters out genes with only 1 exon, have exons on > 1 chromosome/both strands and overlapping genes
if (stranded == TRUE) {
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
