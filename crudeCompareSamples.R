# Purpose: to filter our bed files and annotate them.
# Author: Martin Kavanagh
# Date: 01/08/24
# ***************************************************
rm(list=ls())
library(GenomicRanges)
library(stringr)
library(readr)
library(ggplot2)
library(VennDiagram)
setwd("C:/Users/mjmka/Documents/Bioinformatics")
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/filter_peakFileName.R")) # testing the filter but not needed anymore
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
source(paste0(projPath, "/scripts/subscripts/computeOverlaps.R"))
source(paste0(projPath, "/scripts/subscripts/annotateAndEnrich.R"))
# ***************************************************
# I want to find out how many of the peaks are shared between each sample
# this comparison will use score only
# blackBed <- read.table(file = paste0(projPath, "/data/ncbi/hg38-blacklist.v2.bed" ), fill=T, header = F, sep = "\t")
# colnames(blackBed) <- c("chrom", "start", "end", "reason")
# blacklist.hg38 <- GRanges(ranges=IRanges(start=blackBed$start, end=blackBed$end),
#                           seqnames=blackBed$chrom,
#                           reason=blackBed$reason)
# get hists we want
wantedHists=c("Homer", "encode", "CUTandTag")

fileNames <- vector(mode="character")
fileNames <- findInDirList(path = paste0(projPath, "/peakCalling"), wantedstrings = wantedHists)
fileNames <- fileNames[(!(grepl("stringent",fileNames[file.info(fileNames)$size != 0])))]
fileNames <- fileNames[(!(grepl("2O|Undetermined",fileNames)))]
fileNames
fileNames=fileNames[c(31, 47, 18, 11, 5, 24, 8)]
fileNames <- fileNames[-1]
fileNames
for(fileNme in fileNames){
  filter_peakFileName(fileNme, remove_non_seacr = TRUE)
}
# one thing I want to do though is check out findoverlaps

# overlapMatrix <- matrix(data = 1, nrow=length(fileNames), ncol=length(fileNames), dimnames =list(fileNames, fileNames))
# quantile=0.75
# i=0
# # computeOverlaps is a fancy function needed?
# for(fileTest in fileNames){
#   for(fileControl in fileNames){
#     temp_score_holder <-computeOverlaps(fileTest, fileControl, quantile)
#     overlapMatrix[fileTest, fileControl] <- temp_score_holder[1]
#     overlapMatrix[fileControl, fileTest] <- temp_score_holder[2]
#   }
#   fileNames <- fileNames[fileNames!= fileTest]
#   print(paste0("done", i))
# }
# 
  # for(fileName in fileNames)
  # {
  #   annotateAndEnrich(fileName)
  #   # score threshold,
  # 
  # }
### want to make venn diagram 
# import all annos into a list
GHtypeList=c("EnhancerOnly", "EnhancersAndPromoterEnhancers")
quantList=c( 0, 0.5, 0.75, 0.9, 0.99)
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX",
               "H1_H3K27ac","H1_H3K27me3", "H1_H3K4me1", "H9_H3K27ac","H9_H3K27me3", "H9_H3K4me1", "IPSC_ORLANDO_yH2AX")
myCol=c("red", "blue", "green")
GHList <- list()
wantedHists=c("orlando", "encode", "CUTandTag")

fileNames <- vector(mode="character")
fileNames <- findInDirList(path = paste0(projPath, "/annotationsAndEnrichment"), wantedstrings = wantedHists)
fileNames <- fileNames[(!(grepl("stringent",fileNames[file.info(fileNames)$size != 0])))]
fileNames <- fileNames[(!(grepl("2O|Undetermined",fileNames)))]
fileNames <- fileNames[((grepl("enhancerAnnotated",fileNames)))]
fileNames <- fileNames[(!(grepl("Entrez",fileNames)))]
fileNames

vennDiagramHists = c("IPSC_H3K27ac", "IPSC_yH2AX", "D7_yH2AX")

fileNames2 <- fileNames[grepl(paste0(vennDiagramHists, collapse="|"), fileNames)]
fileNames2
quant=0.99
for(quant in quantList)
{
  for (GHtype in GHtypeList)
  {
    if(GHtype=="EnhancerOnly"){fileNames3 <- fileNames2[(!(grepl("Promoter",fileNames2)))]}else{fileNames3 <- fileNames2}
    for( fileName in fileNames3)
    {
      histName <- sampleList[str_detect(fileName, sampleList)]
      bed <- read.table(file=fileName, header = T, sep="\t", quote="")
      bed <- bed[bed$mcols.score >= quantile(bed$mcols.score, probs=quant),]
      GHList[[histName]] <- unique(bed$GenehancerID[bed$distToGenehancer <200])
    }
    GHList[vennDiagramHists]
    venn.diagram(x = GHList[vennDiagramHists],
                 category.names = vennDiagramHists,
                 filename = paste0('plots/Annotation/Genehancer/VennDiagrams/test_hists_',quant,'_',GHtype  ,'_GH_venn_diagramm.png'),
                 output=FALSE,
                 
                 # Output features
                 imagetype="png" ,
                 height = 1080 , 
                 width = 1920 , 
                 resolution = 700,
                 compression = "lzw",
                 
                 # Circles
                 lwd = 2,
                 lty = 'blank',
                 fill = myCol,
                 
                 # Numbers
                 cex = .6,
                 fontface = "bold",
                 fontfamily = "sans",
                 
                 # Set names
                 cat.cex = 0.6,
                 cat.fontface = "bold",
                 cat.default.pos = "outer",
                 cat.pos = c(-27, 27, 135),
                 cat.dist = c(0.055, 0.055, 0.085),
                 cat.fontfamily = "sans",
                 rotation = 1)
  }
}



GHList
saveRDS(overlapMatrix, file = paste0(projPath, "/plots/Top0.25peaksOverlapCompare_K27ac_p300_H3k4me1_yH2AX.covMatrix"))

p <- prcomp(overlapMatrix)
summary(p)

p <-prcomp(overlapMatrix)
summary(p)
scoreThresholdList= c("top10percent", "top25percent")

for (scoreThreshold in scoreThresholdList)
{
  ourGR_list <- readRDS(file=paste0(projPath, "/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_", scoreThreshold ,"_controlgone_FIXED.peaks.GRangesList"))  
  for (wantedHist in c("H3K27ac", "H3K27me3"))
  {
    # init the matrix
    sampleList <- c(names(ourGR_list[grepl(wantedHist, names(ourGR_list), fixed = TRUE)]), names(encode_GRlist[grepl(wantedHist, names(encode_GRlist), fixed = TRUE)]))
    compDF <- matrix(nrow = length(sampleList), ncol = length(sampleList))
    colnames(compDF) <- sampleList
    rownames(compDF) <- sampleList
    newSampleList <-sampleList
    for (histName in newSampleList)
    {
      if(length(ourGR_list[[histName]]))
      {
        histTest <- ourGR_list[[histName]]
        for (histCompare in newSampleList)
        {
          if(length(ourGR_list[[histCompare]]))
          {
            histControl <- ourGR_list[[histCompare]] 
            compDF[histName, histCompare] <- length(subsetByOverlaps(histTest, histControl, maxgap = mxGap))/length(histTest)
            compDF[histCompare, histName] <- length(subsetByOverlaps(histControl, histTest, maxgap = mxGap))/length(histControl)
          }
          else if(length(encode_GRlist[[histCompare]]))
          {
            histControl <- encode_GRlist[[histCompare]]
            testInCompare=0
            compareInTest=0
            for (subHistCompare in names(histControl))
            {
              # length of the subset/length of the test 
              histSubControl <- histControl[[subHistCompare]]
              testInCompare = testInCompare + (length(subsetByOverlaps(histTest, histSubControl, maxgap = mxGap))/length(histTest))
              compareInTest = compareInTest + (length(subsetByOverlaps(histSubControl, histTest, maxgap = mxGap))/length(histSubControl))
            }
            compDF[histName, histCompare] <- testInCompare/length(histControl)
            compDF[histCompare, histName] <- compareInTest/length(histControl)
          }
          else {print("errror finding histControl")}
        }
      }
      else if(length(encode_GRlist[[histName]]))
      {#if its not ours, get it from encode
        histTest <- encode_GRlist[[histName]]
        
        for (histCompare in newSampleList)
        {
          if(length(ourGR_list[[histCompare]]))
          {
            histControl <- ourGR_list[[histCompare]] 
            testInCompare=0
            compareInTest=0
            for (subHistName in names(histTest))
            {
              histSubTest <- histTest[[subHistName]]
              testInCompare = testInCompare + (length(subsetByOverlaps(histSubTest, histControl, maxgap = mxGap))/length(histSubTest))
              compareInTest = compareInTest + (length(subsetByOverlaps(histControl, histSubTest, maxgap = mxGap))/length(histControl))
              
            }
            compDF[histName, histCompare] <- testInCompare/length(histTest)
            compDF[histCompare, histName] <- compareInTest/length(histTest)
          }
          else if(length(encode_GRlist[[histCompare]]))
          {
            histControl <- encode_GRlist[[histCompare]]
            testInCompare=0
            compareInTest=0
            for (subHistName in names(histTest))
            {
              histSubTest <- histTest[[subHistName]]
              for (subHistCompare in names(histControl))
              {
                histSubControl <- histControl[[subHistCompare]]
                testInCompare = testInCompare + (length(subsetByOverlaps(histSubTest, histSubControl, maxgap = mxGap))/length(histSubTest))
                compareInTest = compareInTest + (length(subsetByOverlaps(histSubControl, histSubTest, maxgap = mxGap))/length(histSubControl))
              }
              
            }
            compDF[histName, histCompare] <- testInCompare/(length(histTest)*length(histControl))
            compDF[histCompare, histName] <- compareInTest/(length(histTest)*length(histControl))
            
          }
          else {print("errror finding histControl")}
        }
      }
      else {print("errror finding histTest")}
      
      #remove histName from sampleList
      
      newSampleList <- sampleList[sampleList != histName]
    }
    write.table(compDF, file=paste0(projPath, "/R_beds/bulkData/", wantedHist, "_", scoreThreshold, "_crossComparison_", mxGap,"BPgap.txt"))
  }
}
bed <-read.table(paste0(projPath, "/peakCalling/encode", "/H1", "/H1_H3K27ac_ENCSR000ANP", "/ENCFF702RUG_MACS2_peaks.bed"), sep="\t")
bed
wantedHist="H3K27ac"
scoreThreshold="top25percent"
compDF <- read.table(file=paste0(projPath, "/R_beds/bulkData/", wantedHist, "_", scoreThreshold, "_crossComparison_", mxGap,"BPgap.txt"))
comp_pca <- prcomp(compDF, center=TRUE, scale=TRUE)
summary(comp_pca) 
comp_pcaDF <- data.frame(comp_pca$rotation)
compDF$names <- row.names(compDF)
compDF[grepl("H1", names(compDF), fixed = TRUE),"colors"]  <- "red"   
compDF[grepl("H9", names(compDF), fixed = TRUE),"colors"]  <- "green"   
compDF[grepl("Nephron", names(compDF), fixed = TRUE),"colors"]  <- "pink"   
compDF[grepl("NPC", names(compDF), fixed = TRUE),"colors"]  <- "pink"
compDF[grepl("IPSC", names(compDF), fixed = TRUE),"colors"]  <- "black"
compDF[grepl("D7", names(compDF), fixed = TRUE),"colors"]  <- "blue"   
comp_pcaDF


library(ggfortify) 
ggplot(data=comp_pca,mapping=aes(x=PC1, y=PC2))+
  geom_point(colour=compDF$colors, size=6)

biplot_comp <- biplot(comp_pca)    

library(threejs)
z <- comp_pca$rotation[,"PC3"]
x <- comp_pca$rotation[,"PC1"]
y <- comp_pca$rotation[,"PC2"]
scatterplot3js(x,y,z, color=compDF$colors)

library(plotly)
plot_ly(comp_pcaDF, x = comp_pcaDF[,"PC1"], y = comp_pcaDF[,"PC2"], 
        text = paste("Clarity: ", row.names(comp_pcaDF)),
        mode = "markers", color = compDF$colors, size = comp_pcaDF[,"PC3"])

# ok so this is clearly a sign that iPSCs are not good equivalents to ESCs
# lets try find some data online

"~/Bioinformatics/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_top10percent_controlgone_FIXED.peaks.GRangesList"
"~/Bioinformatics/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_top25percent_controlgone_FIXED.peaks.GRangesList"
"~/Bioinformatics/R_beds/encode/ENCODE_GRanges.GRangesList"

blackBed <- read.table(file = paste0(projPath, "/R_beds/hg38-blacklist.v2.bed" ), fill=T, header = F, sep = "\t")
colnames(blackBed) <- c("chrom", "start", "end", "reason")
blacklist.hg38 <- GRanges(ranges=IRanges(start=blackBed$start, end=blackBed$end),
                          seqnames=blackBed$chrom,
                          reason=blackBed$reason)



# is A in B the same as B in A? no!

# test granges just to confirm overlap behaviour
a <- data.frame(chr="chr1",
                start=c(1,10000,20000),
                end=c(50, 10020, 25000))
a
b <- data.frame(chr = "chr1",
                start=c(10,10010,20100),
                end=c(40,10200,21000))
grA <- GRanges( seqnames=a$chr,
                ranges=IRanges(start=a$start,
                               end=a$end))
grB <- GRanges( seqnames=b$chr,
                ranges=IRanges(start=b$start,
                               end=b$end))
subsetByOverlaps(grA, grB)
subsetByOverlaps(grA, grA, maxgap = 100000)

findOverlaps(grA, grB)
)


