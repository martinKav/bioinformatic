#source file for filter_peakfileName function
library(GenomicRanges)
library(stringr)
library(readr)
filter_peakFileName <- function(fileName, remove_non_seacr = "FALSE"){
  if (!(exists("blacklist.hg38"))){
    blackBed <- read.table(file = paste0(projPath, "/data/ncbi/hg38-blacklist.v2.bed" ), fill=T, header = F, sep = "\t")
    colnames(blackBed) <- c("chrom", "start", "end", "reason")
    blacklist.hg38 <- GRanges(ranges=IRanges(start=blackBed$start, end=blackBed$end),
                              seqnames=blackBed$chrom,
                              reason=blackBed$reason)
  }
  print(paste0("filtering " , fileName))
  if(!is.na(file.size(fileName)) && !(file.size(fileName) == 0)){
    bed <- read.table(file=fileName, sep= "\t")
    bed <-subset(bed, bed[,1] %in% c("chr1", "chr2", "chr3", "chr4","chr5", 
                                     "chr6", "chr7", "chr8", "chr9","chr10", 
                                     "chr11", "chr12", "chr13", "chr14","chr15",
                                     "chr16", "chr17", "chr18", "chr19","chr20",
                                     "chr21", "chr22", "chrX", "chrY"))
    if( !(any(grepl("seacr", strsplit(fileName, "\\.")[[1]]))) && (strsplit(fileName, "\\.")[[1]][length(strsplit(fileName, "\\.")[[1]])] == "narrowPeak" || any(grepl("MACS2", strsplit(fileName, "\\.")[[1]]))))
    {
      print("converting to seacr")
      gr <- GRanges(seqnames = bed$V1, ranges = IRanges(start= bed$V2, end= bed$V3),
                    score=bed$V5,TotalSignal=(bed$V7 * (((bed$V3)+1) - (bed$V2))),
                    MaxSignalRegion = paste0(bed$V1, ":", IRanges(start= bed$V2, end= bed$V3)))
    }else if(strsplit(fileName, "\\.")[[1]][length(strsplit(fileName, "\\.")[[1]])] != "bed") {print("ERROR. Encountered non bed, macs2 or narrowpeak file.")
    }else if(any(grepl("seacr", strsplit(fileName, "\\.")[[1]])))
    {
      gr <- GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2, end=bed$V3),
                    TotalSignal=bed$V4,score=bed$V5,
                    MaxSignalRegion=bed$V6)
    }
    gr
    if (length(gr) != 0){
      blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
      if (length(blacklisted) != 0){gr <- gr[-from(blacklisted)]}
      bed2 <- data.frame(chr=seqnames(gr), start= start(gr), end = end(gr), TotalSignal=gr$TotalSignal,score = score(gr), MaxSignalRegion=gr$MaxSignalRegion)
      if(length(bed2)!=0)
      {
        if (!(any(grepl("seacr", strsplit(fileName, "\\.")[[1]])))){
          write.table(bed2, file=paste0(str_sub(fileName, end = -(nchar(str_split(fileName, "\\.")[[1]][length(str_split(fileName, "\\.")[[1]])])+1)), "seacrised.bed"), sep="\t", row.names =F, col.names =F, quote=F, eol='\n')
          if(remove_non_seacr == "TRUE"){file.remove(fileName)}
          #on.exit(close(file(paste0(str_sub(fileName, end = -(nchar(str_split(fileName, "\\.")[[1]][length(str_split(fileName, "\\.")[[1]])])+1)), "seacrised.bed"),"wb")), add = TRUE)
          #on.exit(remove(bed, gr, bed2), add = TRUE)
        }else
          {write.table(bed2, file=fileName, sep="\t", row.names =F, col.names =F, quote=F, eol='\n')}
        print(paste0("Removed ", length(bed) - length(bed2), " peaks from ", fileName))
        #on.exit(close(fileName), add = TRUE)
        #on.exit(remove(bed, gr, bed2), add = TRUE)
      }else {print("Bed could not be written.")}
    }else {print("no GRanges made")}
  }
  else{ print(paste0(fileName, " is an empty file, or does not exist!"))}
}