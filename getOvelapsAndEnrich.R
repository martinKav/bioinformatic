# this script will get the overlaps for each sample, graph them, run enrichment for them, findmotifs, and write the gene lists into a table
# Author: Martin Kavanagh
# Date: 18/08/24
# ********************************************************************
rm(list=ls())
shell("cls")
library(stringr)
library(readr)
library(plyr)
library(dplyr)
library(GenomicRanges)
library(venn)
library(UpSetR)
library(ggplot2)
library(ggpolypath)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
source(paste0(projPath, "/scripts/subscripts/buildXwiseListFromVec.R"))
source(paste0(projPath, "/scripts/subscripts/smallTools.R"))
source(paste0(projPath, "/scripts/subscripts/getGRListIntersections.R"))
source(paste0(projPath, "/scripts/subscripts/enrichGRange.R"))
# Lists
GHtypeList=c(
  #"EnhancerOnly", # only if you want, it doubles computational time however so not advised
  "EnhancersAndPromoterEnhancers")
quantList=c( 0.999, 0.9, 0.99, 0.75)
inQuantList=c(0.999, 0.99, 0.9, 0.75) # the quantile levels of our samples included in the intersections
outQuantList=c(0.5, 0.25, 0.1, 0) # the quantile levels of the samples excluded from the intersection 
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX")
wantedCols <- c("seqnames", "start","end",  "TotalSignal", "score", "peakNumber",  "distanceToTSS","distToGenehancer", "GenehancerID", "GenehancerType", "interactionValue","GENENAME", "interactedGeneName", "interactedGeneENSEMBL", "interactedGeneENTREZ")
# find files - we want the annotated peaks, grabs annotated peaks that end in txt
allFileNames <- findInDirList(paste0(projPath, "/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/"), ".txt")
file.remove(paste0(projPath, "/geneOverlaps/GeneOverlaps/EmptyOverlaps.txt")) # deletes the list of empty overlaps
GRlist <- list() # initialise an empty list to hold our GRanges
for (GHtype in GHtypeList){
  if (GHtype == "EnhancerOnly"){fileNames <- allFileNames[!(grepl("Promoter", allFileNames))]}else{fileNames <- allFileNames[(grepl("Promoter", allFileNames))]}
  for (fileName in fileNames){
    histName <- sampleList[str_detect(fileName, sampleList)] # grab the sample name from the file name
    bed <- read.table(file=fileName, sep="\t", header=T, quote="") # read the file
    if(colnames(bed)[1] == 1){print(paste0("The columns do not have names! for ", fileName))} # if the file has no column names, alert 
    colnames(bed) <- sapply(colnames(bed), mcolsRemove) # remove the mcols prefix that appears on the metadata columns of GRanges objects after construction
    bed <- bed[,match(wantedCols, colnames(bed))] # keeps the wanted columns of the bed data frame
    bed2 <- bed[bed$distToGenehancer <200 | bed$distanceToTSS <2000,] # keeps peaks that are 200bp away from an enhancer or 2000bp away from a TSS
    gr <- GRanges(seqnames=bed2$seqnames, ranges=IRanges(start=bed2$start, end=bed2$end),
                  mcols=bed2[,colnames(bed2)[!(grepl("start|end|seqnames", colnames(bed2)))]]) # makes a GRanges object from our bed data frame
    colnames(mcols(gr)) <- sapply(colnames(mcols(gr)),mcolsRemove) # remove the mcols prefix
    GRlist[[histName]] <- gr # add the GRanges to an array
  }
  print("finished grabbing beds")
  GRlist <- as(GRlist, "GRangesList") # turn our list of GRanges into a GRangesList
  for(inQuant in inQuantList) 
  {
    for(outQuant in outQuantList)
    {
      # get histCompareList, then filter by score, then get interactions
    for (N in 2:6)
    {
     histCompareList <- buildXwiseListFromVec(names(GRlist), elementsN = N) # build a list of combinations of comparisons given a list
     
     for (comparingGRs in histCompareList)
     {
         GRintersectionList <- getGRListIntersections(GRlist[comparingGRs], inQuant = inQuant, outQuant = outQuant) # get GRanges of these comparisons
         enrichGRange(GRlist = GRintersectionList, projPath = projPath, inQuant = inQuant, outQuant = outQuant, GHtype=GHtype, elementsN = N) # enrich these comparisons
         
         ###this code can generate a venn diagram and upset plot but the upset plot sucks - also the numbers aren't quite the same
         #newGRlist <- lapply(GRlist, filterByTotalSignal, quant)
        # for(GRintersection in GRintersectionList){enrichGRange(GRintersection, projPath, quant, GHtype)}
       # we can also just do a venn diagram without having to worry about the intersections - and UpsetR?
       # if(N != 1)
       # {
       #   mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/Venn/", N, "/GENENAME"))
       #   mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/Venn/", N, "/ENHANCERID"))
       #   vGeneNameList <- list()
       #   vEnhancerList <- list()
       #   for (histName in names(newGRlist[comparingGRs]))
       #   {
       #      vEnhancerList[[histName]] <- unique(c(newGRlist[[histName]]$GenehancerID[newGRlist[[histName]]$distToGenehancer<200]))
       #      vGeneNameList[[histName]] <- unique(c(newGRlist[[histName]]$GENENAME[newGRlist[[histName]]$distanceToTSS < 2000], newGRlist[[histName]]$interactedGeneName[newGRlist[[histName]]$distToGenehancer<200]))
       #   }
       #   geneVennPlot <- venn(vGeneNameList, ilabels = "counts", zcolor = "red,blue,green,yellow,purple,darkblue,darkgreen",
       #                        box=FALSE, ellipse = TRUE, ggplot=TRUE, size=0.5, ilcs=2, sncs=2, plotsize = 30)
       #   enhancerVennPlot <- venn(vEnhancerList, ilabels = "counts", zcolor = "red,blue,green,yellow,purple,darkblue,darkgreen",
       #                            box=FALSE, ellipse = TRUE, ggplot=TRUE, size=0.5, ilcs=2, sncs=2, plotsize =30)
       #   if(N !=6)
       #   {
       #   png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/Venn/", N,
       #                   "/GENENAME/", paste0(paste0(names(GRlist[comparingGRs]), collapse = "_AND_"), "_GENENAME_Quant",quant, "_GHtype_", GHtype, "_vennDiagram.png")),
       #                   width = 1920,height=1080,  pointsize = 100,antialias = "subpixel", type = "cairo-png")
       #   }else{png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/Venn/", N,
       #                         "/GENENAME/All_samples_GENENAME_Quant",quant, "_GHtype_", GHtype, "_vennDiagram.png"),
       #             width = 1920,height=1080,  pointsize = 100,antialias = "subpixel", type = "cairo-png")}
       #   
       #   print(geneVennPlot)
       #   dev.off()
       #   if(N !=6)
       #   {
       #   png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/Venn/", N,
       #                   "/ENHANCERID/", paste0(paste0(names(newGRlist[comparingGRs]), collapse = "_AND_"), "_ENHANCERID_Quant",quant, "_GHtype_", GHtype, "_vennDiagram.png")),
       #       width = 1920,height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
       #     }else{png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/VennDiagrams/Venn/", N,
       #                          "/ENHANCERID/All_samples_ENHANCERID_Quant",quant, "_GHtype_", GHtype, "_vennDiagram.png"),
       #              width = 1920,height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")}
       #   print(enhancerVennPlot)
       #   dev.off()
       #   
       #   if(N > 3)
       #   {
       #     mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/UpSet/", N, "/ENHANCERID"))
       #     mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/UpSet/", N, "/GENENAME"))
       #   geneUpset <- upset(fromList(vGeneNameList),mb.ratio = c(0.55, 0.45), order.by = "freq", number.angles = 2, point.size = 2, line.size = 1, 
       #         mainbar.y.label = "Number of Unique Genes Interacted With At Promoters and Enhancers", sets.x.label = "Genes Per Sample",
       #         text.scale = c(2.5, 3, 3, 3, 4, 3))
       #   enhancerUpset <- upset(fromList(vEnhancerList),mb.ratio = c(0.55, 0.45), order.by = "freq", number.angles = 2, point.size = 2, line.size = 1, 
       #                      mainbar.y.label = "Number of Unique Enhancers annotated to Each sample", sets.x.label = "Unique Enhancers Per Sample",
       #                      text.scale = c(2.5, 3, 3, 3, 4, 3))
       #   png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Upset/", N,
       #                   "/GENENAME/", paste0(paste0(names(GRlist[comparingGRs]), collapse = "_AND_"), "_GENENAME_Quant",quant, "_GHtype_", GHtype, "_Upset.png")),
       #       width = 1920,height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
       #   print(geneUpset)
       #   dev.off()
       #   png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Upset/", N,
       #                   "/ENHANCERID/", paste0(paste0(names(GRlist[comparingGRs]), collapse = "_AND_"), "_ENHANCERID_Quant",quant, "_GHtype_", GHtype, "_Upset.png")),
       #       width = 1920,height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
       #   print(enhancerUpset)
       #   dev.off()
       #   }
       }
     #}
    }
    }
  }
}

