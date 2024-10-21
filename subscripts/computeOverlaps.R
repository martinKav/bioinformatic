# Script for Computing Overlaps
library(GenomicRanges)
computeOverlaps <- function(fileTestName, fileControlName, mode= "default", quantile = 0)
{
  if(fileTestName == fileControlName){ return(c(1,1))}
  else{
    bed1 <- read.table(fileTestName, sep="\t")
    bed2 <- read.table(fileControlName, sep="\t")
    #if (mode=="default") will do later
    {
      gr1 <- GRanges(seqnames = bed1$V1, ranges = IRanges(start=bed1$V2, end=bed1$V3), score = bed1$V5)
      gr2 <- GRanges(seqnames = bed2$V1, ranges = IRanges(start=bed2$V2, end=bed2$V3), score = bed2$V5)
      gr1 <- gr1[score(gr1) >= quantile(score(gr1), probs=quantile)]
      gr2 <- gr2[score(gr2) >= quantile(score(gr2), probs=quantile)]
      testInControl <- sum(!is.na(findOverlaps(gr1, gr2, type="any", select="first")))/length(gr1)
      controlInTest <- sum(!is.na(findOverlaps(gr2, gr1, type="any", select="first")))/length(gr2)
      return(c(testInControl, controlInTest))
    }
  }
}

