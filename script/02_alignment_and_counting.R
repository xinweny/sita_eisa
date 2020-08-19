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

p <- add_argument(p, "-m",
                  help="path to SraRunTable.txt metadata file")

args <- parse_args(p)

#### Load data ####
gse <- tail(strsplit(getwd(), "/")[[1]], n=1)

# Make SampleFile
sampleFiles <- Sys.glob("SampleFile*.txt")

if (length(sampleFiles) == 0) {
  message("Creating QuasR sample file(s)...")
  
  system(glue("python3 /home/xwy21/project/sita/script/generate_samplefile.py -m {args$m} -e fastq"))
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
cl <- makeCluster(8)

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
  exonCount <- qCount(proj, regions$exons, orientation=orient, reportLevel="gene", clObj=cl)
  genebodyCount <- qCount(proj, regions$genebodies, orientation=orient, reportLevel="gene",clObj=cl)
  intronCount <- genebodyCount - exonCount
  
  # Remove width column
  exonCount <- exonCount[, -1]
  intronCount <- intronCount[, -1]
  
  # Remove version from genes
  row.names(exonCount) <- gsub("\\.[0-9]*", "", rownames(exonCount))
  row.names(intronCount) <- gsub("\\.[0-9]*", "", rownames(intronCount))
  
  if (grepl("PAIRED", sampleFile, fixed=TRUE)) {
    counts <- c(exonCount, intronCount)
    
    countsList <- list()
    
    for (count in counts) {
      row.names(count) <- count[, 1]
      count[, 1] <- NULL
      
      xcols <- count %>% dplyr::select(contains(".x")) %>% colnames()
      ycols <- count %>% dplyr::select(contains(".y")) %>% colnames()
      newcols <- gsub("\\.x", "", xcols)
      
      for (i in seq_along(newcols)) {
        newcol <- as.data.frame(rowSums(count[, c(xcols[i], ycols[i])]))
        names(newcol) <- newcols[i]
        count <- cbind(count, newcol)
      }
      
      count <- subset(count, select=newcols)
      
      countsList[[length(countsList) + 1]] <- count
    }
    
    exonCount <- countsList[[1]]
    intronCount <- countsList[[2]]
  }
  
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

write.table(exonCount, file=glue("processed/ExonicCounts_{gse}.txt"), row.names=TRUE, col.names=TRUE, sep="\t")
write.table(intronCount, file=glue("processed/IntronicCounts_{gse}.txt"), row.names=TRUE, col.names=TRUE, sep="\t")
