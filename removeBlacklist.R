#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
blacklistFileName=args[1]
readFileName=args[2]
writeFileName=args[3]

library(stringr)
library(GenomicRanges)
library(here)

blacklist <- read.table(here(blacklistFileName), sep = "\t")
blackGR <- sortSeqlevels(GRanges(seqnames=blacklist$V1,
                   ranges=IRanges(start=blacklist$V2,
                   end=blacklist$V3)))
bed <- read.table(here(readFileName), header=F, sep="\t")
gr <- sortSeqlevels(GRanges(seqnames=bed$V1,
              ranges=IRanges(bed$V2,bed$V3),
              score=bed$V4))
wantedGR <-subsetByOverlaps(gr, blackGR, type="any", invert = TRUE)
grDF <- data.frame(chr=seqnames(wantedGR),
           start=start(wantedGR),
           end=end(wantedGR),
           score=score(wantedGR))
write.table(grDF, file=here(writeFileName), quote=F, row.names = F, col.names = F)


