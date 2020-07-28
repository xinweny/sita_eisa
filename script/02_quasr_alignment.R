#### Packages ####
library(QuasR)

setwd("/Users/Pomato/mrc/project/sita") # "/Users/Pomato/mrc/project/sita" "/rds/project/rs2099/rds-rs2099-toxgenomics/sita"

proj <- qAlign(sampleFile="./list/sita_sample.txt", # "./list/sita_samples.txt" "./list/sita_sample.txt"
               genome="BSgenome.Hsapiens.UCSC.hg18",
               aligner="Rhisat2",
               splicedAlignment=TRUE)