#R script to scale bedgraphs to 1
library(stringr)
library(readr)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
#read bedgraph
fileNames <- findInDirList(paste0(projPath, "/alignment/bedgraph"), ".bedgraph")
fileNames <- fileNames[!(grepl("2O|Undetermined", fileNames))]
fileNames
for (fileName in fileNames)
{
  #read bed
  #formula involve the max
  bed <- read.table(fileName, sep="\t", header=F, quote="")
  bed$V4 <- bed$V4/quantile(bed$V4, probs=0.999)
  write.table(bed, file=paste0(projPath,"/alignment/bedgraph/normalised1to0point999/",strsplit(basename(fileName), ".bedgraph")[[1]], ".normalised1to0point999.bedgraph"), sep="\t", col.names=F, row.names = F, quote = FALSE)
}
