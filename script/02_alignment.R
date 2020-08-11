#### Packages ####
suppressPackageStartupMessages({
  library(QuasR)
  library(GenomicFeatures)
  library(argparser)
})

#### Parser ####
p <- arg_parser("QuasR alignment for EISA")

p <- add_argument(p, "-i",
                  help="path to sample file")
p <- add_argument(p, "-g",
                  help="path to genome file (.fa, .fasta)")
p <- add_argument(p, "-a",
                  help="path to annotation file (txdb sqlite format)")

args <- parse_args(p)

# Make sure SampleFile exists
if (!file.exists("SampleFile.txt")) {
  system("python3 /home/xwy21/project/sita/script/generate_samplefile.py -e fastq")
}

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
write.table(alignmentStats(proj), file="bam/AlignmentStats.txt", row.names=TRUE, col.names=TRUE, sep="\t")
