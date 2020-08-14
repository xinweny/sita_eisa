#### Packages ####
suppressPackageStartupMessages({
  library(QuasR)
  library(eisaR)
  library(GenomicFeatures)
  library(argparser)
  library(AnnotationDbi)
  library(glue)
})

#### Parser ####
p <- arg_parser("QuasR alignment and counting for EISA")

p <- add_argument(p, "-m",
                  help="path to SraRunTable.txt metadata file")

args <- parse_args(p)

gse <- tail(strsplit(getwd(), "/")[[1]], n=1)

# Make SampleFile
system(glue("python3 /home/xwy21/project/sita/script/generate_samplefile.py -m {args$m} -e fastq"))

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

#### Alignment ####
proj <- qAlign(sampleFile="SampleFile.txt",
               genome=genomeFile,
               aligner="Rhisat2",
               splicedAlignment=TRUE,
               alignmentsDir="./bam",
               cacheDir="./cache")

# Output alignment stats for each bam file
write.table(alignmentStats(proj), file="bam/AlignmentStats.txt", row.names=TRUE, col.names=TRUE)

#### Counting exons and introns ####
# Check strandedness
bamFile <- Sys.glob("bam/*.bam")[1]
rseqcOutput <- system(glue("infer_experiment.py -r {bedGenome} -i {bamFile}"), intern=TRUE)
message(cat(rseqcOutput, sep="\n"))

percSameStrand <- as.numeric(tail(strsplit(rseqcOutput[5], " ")[[1]], n=1))
stranded <- percSameStrand > 0.9
message(glue("Stranded: {stranded}"))

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
write.table(exonCount, file=glue("processed/ExonicCounts_{gse}.txt"), row.names=TRUE, col.names=TRUE, sep="\t")
write.table(intronCount, file=glue("processed/IntronicCounts_{gse}.txt"), row.names=TRUE, col.names=TRUE, sep="\t")
