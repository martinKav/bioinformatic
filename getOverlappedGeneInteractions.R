# this script gets overlapped gene interactions for genehancer
# Author: Martin Kavanagh
# Date: 14/08/24
rm(list=ls())
shell("cls")
library(plyr)
library(GenomicRanges)
library(stringr)
library(readr)
library(wrapr)
library(ggplot2)
library(VennDiagram)
library(R.utils)
filterByScoreGR <- function(GR, quant){ GR[score(GR) >= quantile(score(GR), probs=quant)]}
filterByTotalSignal <- function(GR, quant){ GR[GR$TotalSignal >= quantile(GR$TotalSignal, probs=quant)]}
mcolsRemove <- function(columnName){if(grepl("mcols", columnName)){newName <- strsplit(columnName, "mcols.")[[1]][2]; return(newName)}else{return(columnName)}}
setwd("C:/Users/mjmka/Documents/Bioinformatics")
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
source(paste0(projPath, "/scripts/subscripts/buildXwiseListFromVec.R"))

GHtypeList=c("EnhancerOnly", "EnhancersAndPromoterEnhancers")
quantList=c( 0, 0.5, 0.75, 0.9, 0.99)
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX",
               "H1_H3K27ac","H1_H3K27me3", "H1_H3K4me1", "H9_H3K27ac","H9_H3K27me3", "H9_H3K4me1", "IPSC_ORLANDO_yH2AX")
wantedCols <- c("seqnames", "start","end",  "TotalSignal", "score", "peakNumber", "GENENAME", "distanceToTSS", "GenehancerID","distToGenehancer", "GenehancerType", "interactedGeneName", "interactionValue")
histThird="ThirdLess"
elementsN=3
wantedHists=c( "CUTandTag" )
fileNames <- vector(mode="character")
fileNames <- findInDirList(path = paste0(projPath, "/annotationsAndEnrichment"), wantedstrings = wantedHists)
fileNames <- fileNames[(!(grepl("stringent",fileNames[file.info(fileNames)$size != 0])))]
fileNames <- fileNames[(!(grepl("2O|Undetermined",fileNames)))]
fileNames <- fileNames[((grepl("enhancerAnnotated",fileNames)))]
fileNames <- fileNames[(!(grepl("Entrez",fileNames)))]
fileNames <- fileNames[(!(grepl("EAaligned",fileNames)))]
for (GHtype in GHtypeList)
{
  if(GHtype == "EnhancerOnly"){fileNames2 <- fileNames[(!(grepl("Promoter",fileNames)))]}else{fileNames2 <- fileNames[((grepl("Promoter",fileNames)))]}
  
  wantedCols <- c("seqnames", "start","end",  "TotalSignal", "score", "peakNumber", "GENENAME", "distanceToTSS", "GenehancerID","distToGenehancer", "GenehancerType", "interactedGeneName", "interactionValue")
  # at this point we only really want a few columns: seqnames, ranges, GHid, score, GENENEAME, distanceToTSS, peakNumber, distToGenehancer, GenehancerType, interactedGeneName, interactionvalue
  GRlist <- list()
  print("grabbing beds")
  for( fileName in fileNames2)
  {
    histName <- sampleList[str_detect(fileName, sampleList)]
    bed <- read.table(file=fileName, header = T, sep="\t", quote="")
    if(colnames(bed)[1] == 1){print(paste0("The columns do not have names! for ", fileName))}
    colnames(bed) <- sapply(colnames(bed), mcolsRemove)
    bed <- bed[,match(wantedCols, colnames(bed))]
    bed <- bed[bed$distToGenehancer <200,]
    gr <- GRanges(seqnames=bed$seqnames, ranges=IRanges(start=bed$start, end=bed$end),
                  mcols=bed[,colnames(bed)[!(grepl("start|end|seqnames", colnames(bed)))]])
    colnames(mcols(gr)) <- sapply(colnames(mcols(gr)),mcolsRemove)
    GRlist[[histName]] <- gr
  }
  print("finished grabbing beds")
  vennDiagramHistList <- buildXwiseListFromVec(sampleList, elementsN = elementsN)
  for (vennDiagramHists in vennDiagramHistList)
  {
    if(!dir.exists(paste0(projPath, "/geneOverlaps/",paste0(vennDiagramHists,collapse = "-"))))
    {
      mkdirs(paste0(projPath, "/geneOverlaps/",paste0(vennDiagramHists,collapse = "-")))
      mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Scatterplots/",paste0(vennDiagramHists,collapse = "-")))
      mkdirs(paste0(projPath, '/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/',paste0(vennDiagramHists,collapse = "-")))
    } 
    # this loop will get us the overlap strictly only in two sets of data
    # I want a flexible function that will let me work with whatever number of elements
    # for overlaps, we need to make distinct from aUbUcUd vs an
    for(quant in quantList)
    {
      histList <- vennDiagramHists 
      # for each set of histones, find the pairwise comparisons - ANBNC, ANB!NC, ANB, ANC
      # could do the easy one and calculate ANBNC first, then subtract it from everything else ------- this is so big brain -> but when there are 4 or more elements mmmm
      newGRlist <- lapply(GRlist, filterByTotalSignal, quant)
      {
        for (histFirst in histList)
        {
          histSecond <- ""
          histList2 <- ""
          histList3 <- ""
          GR1 <- newGRlist[[histFirst]]
          if(length(histList) >1){
            histList2 <- histList[histList != histFirst]
            if(length(histList) >1){}
            for (histSecond in histList2)
            {
              histThird <- ""
              vennList <- list()
              vennList2 <- list()
              vennList2[[histFirst]] <- unique(newGRlist[[histFirst]]$GenehancerID)
              vennList2[[histSecond]] <- unique(newGRlist[[histSecond]]$GenehancerID)
              vennList[[histFirst]] <- unique(newGRlist[[histFirst]]$interactedGeneName)
              vennList[[histSecond]] <- unique(newGRlist[[histSecond]]$interactedGeneName)
              #make DF here
              GR2 <- newGRlist[[histSecond]]
              overlapNames <- unique(GR1$interactedGeneName[(GR1$interactedGeneName %in% GR2$interactedGeneName)])
              if(length(histList2) > 1)
              {
                histList3 <- histList2[histList2 != histSecond]
                for (histThird in histList3){
                  overlapNames <- overlapNames[!(overlapNames %in% newGRlist[[histThird]]$interactedGeneName)]
                  vennList[[histThird]] <- newGRlist[[histThird]]$interactedGeneName
                  vennList2[[histThird]] <- unique(newGRlist[[histThird]]$GenehancerID)}
              }
              geneDF <- data.frame(GENENAMES = overlapNames,
                                   histFirst = vaggregate(GR1[GR1$interactedGeneName %in% overlapNames]$TotalSignal, overlapNames, sum),
                                   histSecond = vaggregate(GR2[GR2$interactedGeneName %in% overlapNames]$TotalSignal, overlapNames, sum),
                                   Interaction1 = vaggregate(GR2[GR2$interactedGeneName %in% overlapNames]$interactionValue, overlapNames, sum),
                                   Interaction2 = vaggregate(GR1[GR1$interactedGeneName %in% overlapNames]$interactionValue, overlapNames, sum)
                                   )
              colnames(geneDF)[colnames(geneDF) == "histFirst"] <- paste0(histFirst, "_TotalSignal")
              colnames(geneDF)[colnames(geneDF) == "histSecond"] <- paste0(histSecond, "_TotalSignal")
              colnames(geneDF)[colnames(geneDF) == "Interaction1"] <- paste0(histFirst, "_Interaction")
              colnames(geneDF)[colnames(geneDF) == "Interaction2"] <- paste0(histSecond, "_Interaction")
              # write to file
             # write.table(geneDF, file = paste0(projPath, "/geneOverlaps/",paste0(vennDiagramHists,collapse = "-"),"/", histFirst, "In", histSecond, "_without", paste(histList3, collapse = "-"), "_overlap_", GHtype,"_quant", quant,".txt"), col.names = T, row.names=F,quote=FALSE, sep="\t")
              # scatterplot
              # scatterplot <- ggplot(data=geneDF, mapping=aes( x = (geneDF[,paste0(histFirst, "_TotalSignal")] * geneDF[,paste0(histFirst, "_Interaction")]),
              #                                                 y= geneDF[,paste0(histSecond, "_TotalSignal")] * geneDF[,paste0(histSecond, "_Interaction")]))+
              #   labs(title=paste0("TotalSignal Distribution of Interacted Genes between ", histFirst, " and ", histSecond, ". ", GHtype,". Quant=", quant),
              #        x=paste0(histFirst, " Log 10 Total Signal"), y= paste0(histSecond, " Log 10 Total Signal"))+
              #   geom_point(size=1)+
              #   scale_x_continuous(trans="log10")+
              #   scale_y_continuous(trans="log10")+
              #   theme_bw()
              # png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Scatterplots/",paste0(vennDiagramHists,collapse = "-"),"/", histFirst, "In", histSecond, "_without", paste(histList3, collapse = "-"), "_overlap_", GHtype,"_quant", quant,"_scatterplot.png"), width = 1920,
              #     height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
              # print(scatterplot + theme_bw(base_size = 15))
              # dev.off()
              
              if(length(vennList) < 3)
              {
                myColours <- c("red", "blue")
                myPos <- c(-27, 27)
                myDist<- c(0.02, 0.02)
              }
              else 
              {
                myColours <- c("red", "blue", "green")
                myPos <- c(-27, 27, 135)
                myDist <- c(0.055, 0.055, 0.085)
              }
              
              {
              #
              venn.diagram(x = vennList,
                           category.names = names(vennList),
                           filename = paste0(projPath, '/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/',paste0(vennDiagramHists,collapse = "-"),'/', histFirst,'-',histSecond,'-',histThird,'-',quant,'_',GHtype  ,'_SharedInteractedGenes_venn_diagramm.png'),
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
                           fill = myColours,

                           # Numbers
                           cex = .6,
                           fontface = "bold",
                           fontfamily = "sans",

                           # Set names
                           cat.cex = 0.6,
                           cat.fontface = "bold",
                           cat.default.pos = "outer",
                           cat.pos = myPos,
                           cat.dist = myDist,
                           cat.fontfamily = "sans")

              venn.diagram(x = vennList2,
                           category.names = names(vennList2),
                           filename = paste0(projPath, '/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/',paste0(vennDiagramHists,collapse = "-"),'/', histFirst,'-',histSecond,'-',histThird,'-',quant,'_',GHtype  ,'_SharedEnhancers_venn_diagramm.png'),
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
                           fill = myColours,

                           # Numbers
                           cex = .6,
                           fontface = "bold",
                           fontfamily = "sans",

                           # Set names
                           cat.cex = 0.6,
                           cat.fontface = "bold",
                           cat.default.pos = "outer",
                           cat.pos = myPos,
                           cat.dist = myDist,
                           cat.fontfamily = "sans")

              }
            }
            histList <- histList[histList != histFirst]
          }
          
        }
      }
      print(paste0("finished Quant of ", quant))
    }
    print("finished vennDiag List")
  }
  }
}
(195+123)/2
(28.887+29.683)/2
(801+353)/2
(34.409+35.561)/2
testDF <- data.frame(row.names = c("IPSC", "D7_DIFF", "D7_deltaCT"))
testDF$GADPH <- c("Undetermined", 25.159, 0)
testDF$LHX1 <- c("Undetermined", 29.285)
testDF$PAX2 <- c("Undetermined", 26.577)
testDF$SIX2 <- c("Undetermined", 34.985)
testDF$SIX2 <- c("Undetermined", )
29.285


library(venn)
library(venn)
venn(6, ilab=TRUE, zcolor = "style")

# write what we can, then do more analysis later