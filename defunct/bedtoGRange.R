# Title: R code for converting bed files to Granges to make usable
# Author: Martin Kavanagh
# Date: 22/07/24
# ***************************************************************
# setup
rm(list=ls())
library(stringr)
library(readr)
library(ggplot2)
library(ChIPseeker)
library(ChIPQC)
library(rtracklayer)
library(txdbmaker)
setwd("C:/Users/mjmka/Documents/Bioinformatics")

# ***************************************************************
# The authors of all of the R packages for ChIP seq analysis were dumb
# instead of sticking with BED files, they use another file format, GRanges
# so I will convert my beds to GRanges

# then once they are in GRanges, we can visualise them in covplot and can compare them and convert them to txdb
# I would like to save them as GRanges

# so the SEdb data is actually very poor quality - lots of different cell types
# big aggregation - use kmeans to get similarity 

# setting dir
dir="C:/Users/mjmka/Documents/Bioinformatics"

# blacklist_hg38 <- prepare from bed file
blackBed <- read.table(file = paste0(dir, "/data/ncbi/hg38-blacklist.v2.bed" ), fill=T, header = F, sep = "\t")
colnames(blackBed) <- c("chrom", "start", "end", "reason")
blacklist.hg38 <- GRanges(ranges=IRanges(start=blackBed$start, end=blackBed$end),
                          seqnames=blackBed$chrom,
                          reason=blackBed$reason)

{
# get cell types of SEdb samples
typeDF <- read.csv(file = "R_beds/SEdb/SEdb 2.0-a comprehensive human and mouse Super-Enhancer database.csv", header= T)
head(typeDF)

#start off with the SEdb, eles first
sampleList = list.files(paste0(dir, "/R_beds/SEdb/beds/"))
sampleList <-sampleList[5:6]
for (sample in sampleList){
  subdir=sample
  fileList <- list.files(paste0(dir, "/R_beds/SEdb/beds/", subdir))
  gr_vector <- vector(mode="list", length=length(fileList))
  gr_vector_names <- c()
  i = 1
  if (str_sub(sample, -3, -1)=="ele")
  {
    for (file in fileList){
      bed <- read.table(paste0(dir, "/R_beds/SEdb/beds/", subdir, "/", file), header=T, fill = T)
      print(file)
      # bed <- bed[bed$ele_score >500,]
      gr <- GRanges(ranges=IRanges(start=bed$ele_start, end=bed$ele_end), 
                    seqnames=bed$ele_chr,
                    se_id=bed$se_id
                    ,cell_id=bed$cell_id,
                    cell_type = typeDF$Biosample.name[typeDF$Sample.ID == paste0("Sample_", str_sub(bed$cell_id[1], 4, 10))]
                    )
      blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
      if (length(blacklisted) != 0){
        gr <- gr[-from(blacklisted)]
      }
      print(file)
      gr_vector_names[i] <- gr$cell_id[1]
      gr_vector[[i]] <- gr
      i = i + 1
      if (!dir.exists(paste0(dir, "/R_beds/SEdb//gff3/", subdir))) {dir.create(paste0(dir, "/R_beds/SEdb//gff3/", subdir))}
      export(object = gr, con = paste0(dir, "/R_beds/SEdb//gff3/", subdir, "/", file,".gff3"))
      GH_PE_txdb<- txdbmaker::makeTxDbFromGFF(file =paste0(dir, "/R_beds/SEdb//gff3/", subdir, "/", file,".gff3"), dataSource="SEdb", organism="Homo sapiens")
      saveDb(GH_PE_txdb, file=paste0(dir, "/R_beds/SEdb//gff3/", subdir, "/", file,"_txdb.sqlite"))
      
    }
    print(file)
    names(gr_vector) <- gr_vector_names
    GR_list <- as(gr_vector, "GRangesList")
    GR_list
    saveRDS(object = GR_list,file = paste0(dir, "/R_beds/SEdb/GRangesLists/", subdir, "_filtered.GRangesList"))
    # p <- covplot(GR_list)
    # ggsave(filename = paste0(dir, "/R_beds/SEdb/Covplots/", subdir, ".covplot.svg"), plot = p)
    
    
  } else if (str_sub(sample, -3, -1)=="g38"){
    for (file in fileList){
      bed <- read.table(paste0(dir, "/R_beds/SEdb/beds/", subdir, "/", file), header=T, fill = T)
      # bed <- bed[bed$se_rank >500,]
      print(file)
      gr <- GRanges(ranges=IRanges(start=bed$se_start, end=bed$se_end),
                    seqnames=bed$se_chr,
                    se_id=bed$se_id,
                    cell_id=bed$cell_id
                    ,cell_type = typeDF$Biosample.name[typeDF$Sample.ID == paste0("Sample_", str_sub(bed$cell_id[1], 4, 10))]
                    ,se_rank=bed$se_rank)
      print(file)
      blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
      if (length(blacklisted) != 0){
        gr <- gr[-from(blacklisted)]
      }
      gr_vector_names[i] <- gr$cell_id[1]
      gr_vector[[i]] <- gr
      i = i + 1
      if (!dir.exists(paste0(dir, "/R_beds/SEdb//gff3/", subdir))) {dir.create(paste0(dir, "/R_beds/SEdb//gff3/", subdir))}
      export(object = gr, con = paste0(dir, "/R_beds/SEdb//gff3/", subdir, "/", file, ".gff3"))
      GH_PE_txdb<- txdbmaker::makeTxDbFromGFF(file =paste0(dir, "/R_beds/SEdb//gff3/", subdir, "/", file, ".gff3"), dataSource="SEdb", organism="Homo sapiens")
      saveDb(GH_PE_txdb, file=paste0(dir, "/R_beds/SEdb//gff3/", subdir, "/", file, "_txdb.sqlite"))
    }
    print(file)
    names(gr_vector) <- gr_vector_names
    GR_list <- as(gr_vector, "GRangesList")
    saveRDS(object = GR_list,file = paste0(dir, "/R_beds/SEdb/GRangesLists/", subdir, "_filtered.GRangesList"))
    # p <- covplot(GR_list)
    # ggsave(filename = paste0(dir, "/R_beds/SEdb/Covplots/", subdir, ".covplot.svg"), plot = p)
    
  } else {
    print("PROBLEMMM")
  }
}
# ***** this is for testing when the loop breaks

c(
  "~/Bioinformatics/R_beds/SEdb/beds/Human_Kidney_hg38/SE_02_0430_SE_hg38.bed",
  "~/Bioinformatics/R_beds/SEdb/beds/Human_Kidney_hg38/SE_02_0434_SE_hg38.bed"
)
sampleList = list.files(paste0(dir, "/R_beds/SEdb/beds"))
sampleList

# for (sample in sampleList){
sample = "Human_Kidney_hg38"
  subdir=sample
  fileList <- list.files(paste0(dir, "/R_beds/SEdb/beds/", subdir))
  fileList
  gr_vector <- vector(mode="list", length=length(fileList))
  gr_vector_names <- c()
  i = 1
  file ="SE_02_0430_SE_hg38.bed"
  # if (str_sub(sample, -3, -1)=="ele")
  # {
  #   for (file in fileList){
      bed <- read.table(paste0(dir, "/R_beds/SEdb/beds/", subdir, "/", file), header=T, fill = T)
      print(file)
      gr <- GRanges(ranges=IRanges(start=bed$ele_start, end=bed$ele_end),
                    seqnames=bed$ele_chr,
                    se_id=bed$se_id,
                    cell_id=bed$cell_id,
                    cell_type = typeDF$Biosample.name[typeDF$Sample.ID == paste0("Sample_", str_sub(bed$cell_id[1], 4, 10))])
      blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
      if (length(blacklisted) != 0){
        gr <- gr[-from(blacklisted)]
      }
      print(file)
      gr_vector_names[i] <- gr$cell_id[1]
      gr_vector[[i]] <- gr
      i = i + 1
    # }
    names(gr_vector) <- gr_vector_names
    print(file)
    GR_list <- as(gr_vector, "GRangesList")
    print(GR_list)
    saveRDS(object = GR_list,file = paste0(dir, "/R_beds/SEdb/GRangesLists/", subdir, "_filtered.GRangesList"))
    # p <- covplot(GR_list)
    # ggsave(filename = paste0(dir, "/R_beds/SEdb/Covplots/", subdir, ".covplot.svg"), plot = p)
  # } else if (str_sub(sample, -3, -1)=="g38"){
  #   for (file in fileList){
      bed <- read.table(paste0(dir, "/R_beds/SEdb/beds/", subdir, "/", file), header=T, fill = T)
      print(bed)
      gr <- GRanges(ranges=IRanges(start=bed$se_start, end=bed$se_end),
                    seqnames=bed$se_chr,
                    se_id=bed$se_id,
                    cell_id=bed$cell_id
                    ,cell_type = typeDF$Biosample.name[typeDF$Sample.ID == paste0("Sample_", str_sub(bed$cell_id[1], 4, 10))]
                    ,se_rank=bed$se_rank)
      print(gr)
      blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
      print(blacklisted)
      if (length(blacklisted) != 0){
      gr <- gr[-from(blacklisted)]
      }
      print(gr)
      gr_vector_names[i] <- gr$cell_id[1]
      gr_vector[[i]] <- gr
      i = i + 1
    # }
    names(gr_vector) <- gr_vector_names
    print(gr_vector)
    GR_list <- as(gr_vector, "GRangesList")
    print(file)
    saveRDS(object = GR_list,file = paste0(dir, "/R_beds/SEdb/GRangesLists/", subdir, "_filtered.GRangesList"))
    # p <- covplot(GR_list)
    # ggsave(filename = paste0(dir, "/R_beds/SEdb/Covplots/", subdir, ".covplot.svg"), plot = p)
  # } else {
  #   print("PROBLEMMM")
  # }
# }

# end testing fold
subdir=""
# reading the GRangesLists and printing a coverageplot, just change the subdir
GRL_read <- readRDS(file = paste0(dir, "/R_beds/SEdb/GRangesLists/", subdir, "_filtered.GRangesList"))
p <- covplot(GRL_read)
print(p)
}
{
# for genehancer we will have to change a lot, but there is only 1 loop so no reason to change anything
dir="C:/Users/mjmka/Documents/Bioinformatics"
# mmmmmm there is so much with everything included
bed <- read.table(paste0(dir, "/R_beds/genehancer/genehancer.txt"), header=F, fill = T)
head(bed)
colnames(bed) <- c("chrom","chromStart","chromEnd","name","score","strand",
                   "thickStart","thickEnd","colour","evidenceSources","elementType","eliteness")
BED <- subset(bed, elementType == Etype)
gr <- GRanges(ranges=IRanges(start=bed$chromStart, end=bed$chromEnd),
              seqnames=bed$chrom,
              score=bed$score,
              sitename=bed$name,
              colour=bed$colour)
gr
gr_vector <- vector(mode="list", length=3)
gr_vector_names <- c()
i = 1
for (Etype in unique(bed$elementType)){
  BED <- subset(bed, elementType == Etype)
  gr <- GRanges(ranges=IRanges(start=BED$chromStart, end=BED$chromEnd),
                seqnames=BED$chrom,
                score=BED$score,
                sitename=BED$name,
                colour=BED$colour)
  blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
  if (length(blacklisted) != 0){
    gr <- gr[-from(blacklisted)]
  }
  gr_vector[[i]] <- gr
  if(Etype == "Promoter/Enhancer"){ Etype = "PromEnha"}
  gr_vector_names[i] <- Etype
  i = i +1
  export(object = gr, con = paste0(dir, "/R_beds/genehancer/gff3/", Etype, ".gff3"))
  GH_PE_txdb<- txdbmaker::makeTxDbFromGFF(file =paste0(dir, "/R_beds/genehancer/gff3/", Etype, ".gff3"), dataSource="SEdb", organism="Homo sapiens")
  saveDb(GH_PE_txdb, file=paste0(dir, "/R_beds/genehancer/gff3/", Etype, "_txdb.sqlite"))
}
gr_vector
names(gr_vector) <- gr_vector_names
GR_list <- as(gr_vector, "GRangesList")
saveRDS(object = GR_list,file = paste0(dir, "/R_beds/genehancer/GRangesLists/genehancer_filtered.GRangesList"))

# doing the same for genehancer but raising the score threshold
bed <- bed[bed$score >500,]
bed
gr_vector <- vector(mode="list", length=3)
gr_vector_names <- c()
i = 1
for (Etype in unique(bed$elementType)){
  BED <- subset(bed, elementType == Etype)
  gr <- GRanges(ranges=IRanges(start=BED$chromStart, end=BED$chromEnd),
                seqnames=BED$chrom,
                score=BED$score,
                sitename=BED$name,
                colour=BED$colour)
  blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
  if (length(blacklisted) != 0){
    gr <- gr[-from(blacklisted)]
  }
  gr_vector[[i]] <- gr
  if(Etype == "Promoter/Enhancer"){ Etype = "PromEnha"}
  gr_vector_names[i] <- Etype
  i = i +1
  # export(object = gr, con = paste0(dir, "/R_beds/genehancer/gff3/", Etype, "_500plus.gff3"))
  # GH_PE_txdb<- txdbmaker::makeTxDbFromGFF(file =paste0(dir, "/R_beds/genehancer/gff3/", Etype, "_500plus.gff3"), dataSource="SEdb", organism="Homo sapiens")
  # saveDb(GH_PE_txdb, file=paste0(dir, "/R_beds/genehancer/gff3/", Etype, "_txdb_500plus.sqlite"))
}

gr_vector[1:2]
names(gr_vector) <- gr_vector_names
GR_list <- as(gr_vector[1:2], "GRangesList")
saveRDS(object = GR_list,file = paste0(dir, "/R_beds/genehancer/GRangesLists/genehancer_filtered_500plus.GRangesList"))
}
{
## now grange from peak bed, started with genehancer template
dir="C:/Users/mjmka/Documents/Bioinformatics"
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX", "IPSC_2O", "D7_2O", "Undetermined")
gr_vector <- vector(mode="list", length=length(sampleList))
gr_vector_names <- c()
i=1
for (hist in sampleList){
  bed <- read.table(paste0(dir, "/peakCalling/SEACR/", hist, "_seacr_top0.01_stringent.peaks.stringent.bed"), header=F, fill = T)
  colnames(bed) <- c("chr","start","end","TotalSignal","MaxSignal","MaxSignalRegion")
  gr <- GRanges(ranges=IRanges(start=bed$start, end=bed$end),
              seqnames=bed$chr,
              score=bed$MaxSignal,
              TotalSignal=bed$TotalSignal,
              MaxSignalRegion=bed$MaxSignalRegion)
  gr <-subset(gr, seqnames %in% c("chr1", "chr2", "chr3", "chr4",
                                  "chr5", "chr6", "chr7", "chr8", "chr9",
                                  "chr10", "chr11", "chr12", "chr13", "chr14",
                                  "chr15", "chr16", "chr17", "chr18", "chr19",
                                  "chr20", "chr21", "chr22", "chrX", "chrY", "chrM" ))
  
  seqlevels(gr) <- c("chr1", "chr2", "chr3", "chr4",
                            "chr5", "chr6", "chr7", "chr8", "chr9",
                            "chr10", "chr11", "chr12", "chr13", "chr14",
                            "chr15", "chr16", "chr17", "chr18", "chr19",
                            "chr20", "chr21", "chr22", "chrX", "chrY", "chrM" )
  blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
  if (length(blacklisted) != 0){
    gr <- gr[-from(blacklisted)]
  }
  gr_vector[[i]] <- gr
  gr_vector_names[i] <- hist
  i = i +1
  # export(object = gr, con = paste0(dir, "/R_beds/CUTandTAG/gff3/", hist, "_top0.01.peaks.stringent.gff3"))
  # GH_PE_txdb<- txdbmaker::makeTxDbFromGFF(file =paste0(dir, "/R_beds/CUTandTAG/gff3/", hist, "_top0.01.peaks.stringent.gff3"), dataSource="SEdb", organism="Homo sapiens")
  # saveDb(GH_PE_txdb, file=paste0(dir, "/R_beds/CUTandTAG/gff3/", hist, "_top0.01.peaks.stringent_txdb.sqlite"))
}
names(gr_vector) <- gr_vector_names
GR_list <- as(gr_vector, "GRangesList")
saveRDS(object = GR_list,file = paste0(dir, "/R_beds/CUTandTAG/sample_top0.01.peaks.stringent.GRangesList"))

# for relaxed peaks
dir="C:/Users/mjmka/Documents/Bioinformatics"
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX", "IPSC_2O", "D7_2O", "Undetermined")
gr_vector <- vector(mode="list", length=length(sampleList))
gr_vector_names <- c()
i=1
for (hist in sampleList){
  bed <- read.table(paste0(dir, "/peakCalling/SEACR/", hist, "_seacr_top0.01_relaxed.peaks.relaxed.bed"), header=F, fill = T)
  colnames(bed) <- c("chr","start","end","TotalSignal","MaxSignal","MaxSignalRegion")
  gr <- GRanges(ranges=IRanges(start=bed$start, end=bed$end),
                seqnames=bed$chr,
                score=bed$MaxSignal,
                TotalSignal=bed$TotalSignal,
                MaxSignalRegion=bed$MaxSignalRegion)
  gr <-subset(gr, seqnames %in% c("chr1", "chr2", "chr3", "chr4",
                                  "chr5", "chr6", "chr7", "chr8", "chr9",
                                  "chr10", "chr11", "chr12", "chr13", "chr14",
                                  "chr15", "chr16", "chr17", "chr18", "chr19",
                                  "chr20", "chr21", "chr22", "chrX", "chrY", "chrM" ))
  seqlevels(gr) <- c("chr1", "chr2", "chr3", "chr4",
                     "chr5", "chr6", "chr7", "chr8", "chr9",
                     "chr10", "chr11", "chr12", "chr13", "chr14",
                     "chr15", "chr16", "chr17", "chr18", "chr19",
                     "chr20", "chr21", "chr22", "chrX", "chrY", "chrM" )
  blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
  if (length(blacklisted) != 0){
    gr <- gr[-from(blacklisted)]
  }
  if (any(start(gr[seqnames(gr) == "chrM"]) <1)){
    start(gr[seqnames(gr) == "chrM"])[start(gr[seqnames(gr) == "chrM"])<1] <- 1
    print(paste0("adjusting chrM ", hist))
  }else{
    print(paste0("no chrM adjusting needed ", hist))
  }
  gr_vector[[i]] <- gr
  gr_vector_names[i] <- hist
  i = i +1
  # export(object = gr, con = paste0(dir, "/R_beds/CUTandTAG/gff3/", hist, "_top0.01.peaks.relaxed.gff3"))
  # GH_PE_txdb<- txdbmaker::makeTxDbFromGFF(file =paste0(dir, "/R_beds/CUTandTAG/gff3/", hist, "_top0.01.peaks.relaxed.gff3"), dataSource="SEdb", organism="Homo sapiens")
  # saveDb(GH_PE_txdb, file=paste0(dir, "/R_beds/CUTandTAG/gff3/", hist, "_top0.01.peaks.relaxed_txdb.sqlite"))
}
names(gr_vector) <- gr_vector_names
GR_list <- as(gr_vector, "GRangesList")
saveRDS(object = GR_list,file = paste0(dir, "/R_beds/CUTandTAG/sample_top0.01.peaks.relaxed.GRangesList"))
}

#for encode beds
projPath="C:/Users/mjmka/Documents/Bioinformatics"
dirList <-list.files(paste0(projPath, "/peakCalling/encode"))

GRlist_vector <- list()
for (Dir in dirList)
{
  subdirList <- list.files(paste0(projPath, "/peakCalling/encode/", Dir))
  gr_vector <- list()
    for (subdir in subdirList)
    {
      fileList <- list.files(paste0(projPath, "/peakCalling/encode/", Dir, "/", subdir))
      gr_vector <- vector(mode="list")
      for (file in fileList)
      {
       bed <- read.table(file = paste0(projPath, "/peakCalling/encode/", Dir, "/", subdir, "/", file), header = F, fill = T)[,1:5] 
       colnames(bed) <- c("chr","start","end", "name", "score")
       gr <- GRanges(ranges=IRanges(start=bed$start, end=bed$end),
                     seqnames=bed$chr,
                     score=bed$score,
                     name=bed$name)
       gr <-subset(gr, seqnames %in% c("chr1", "chr2", "chr3", "chr4",
                                       "chr5", "chr6", "chr7", "chr8", "chr9",
                                       "chr10", "chr11", "chr12", "chr13", "chr14",
                                       "chr15", "chr16", "chr17", "chr18", "chr19",
                                       "chr20", "chr21", "chr22", "chrX", "chrY", "chrM" ))
       seqlevels(gr) <- c("chr1", "chr2", "chr3", "chr4",
                          "chr5", "chr6", "chr7", "chr8", "chr9",
                          "chr10", "chr11", "chr12", "chr13", "chr14",
                          "chr15", "chr16", "chr17", "chr18", "chr19",
                          "chr20", "chr21", "chr22", "chrX", "chrY", "chrM" )
       blacklisted <- findOverlaps(gr, blacklist.hg38, type="within")
       if (length(blacklisted) != 0){
         gr <- gr[-from(blacklisted)]
       }
       gr_vector[[file]] <- gr
      }
    GRlist_vector[[subdir]] <- gr_vector
    }
}
GRlist <- as(GRlist_vector, "GRangesList")
saveRDS(GRlist_vector, file=paste0(projPath, "/R_beds/encode/ENCODE_GRanges.GRangesList"))

# for IPSC_yH2AX data
projPath="C:/Users/mjmka/Documents/Bioinformatics"
