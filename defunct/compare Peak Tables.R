library(dplyr)
library(plyr)
library(stringr)
library(readr)
library(GenomicRanges)
library(R.utils)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
mcolsRemove <- function(columnName){if(grepl("mcols", columnName)){newName <- strsplit(columnName, "mcols.")[[1]][2]; return(newName)}else{return(columnName)}}
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX",
               "H1_H3K27ac","H1_H3K27me3", "H1_H3K4me1", "H9_H3K27ac","H9_H3K27me3", "H9_H3K4me1", "IPSC_ORLANDO_yH2AX")
genesOfInterest = c("POU5F1", "SOX2","NANOG","SIX1", "SIX2","SIX4", "CITED1","CITED2", "OSR1", "EYA1", "LHX1", "PAX2", "PAX8", "HOXD11","SALL1", "WT1", 
                    "GDNF","GLI1", "GLI2", "GLI3", "SMAD2", "SMAD3","WNT1", "DVL1", "DVL2","DVL3", "FZD1", "FZD2", "FZD3", "NEUROG1", "NEUROG2", "TP53", "GAPDH", "ACTB",
                    "TBXT", "RUNX1", "PAX6", "ASCL1","BMP2","BMP4","BMP7", "FGFR1", "FGFR2", "KLF6", "MYF5", "MRF4", "MYOD", "EBP1", "PAX3", "PAX7")
wantedCols <- c("interactedGeneName", "score", "TotalSignal", "GENENAME", "distanceToTSS", "distToGenehancer", "interactionValue", "peakNumber" )
GHtypeList=c("EnhancerOnly","EnhancersAndPromoterEnhancers")
files <- c(
c(
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/IPSC_yH2AX_EnhancersAndPromoterEnhancers_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/IPSC_yH2AX_EnhancerOnly_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/IPSC_H3K27me3_EnhancersAndPromoterEnhancers_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/IPSC_H3K27me3_EnhancerOnly_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/IPSC_H3K27ac_EnhancersAndPromoterEnhancers_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/IPSC_H3K27ac_EnhancerOnly_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/D7_yH2AX_EnhancersAndPromoterEnhancers_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/D7_yH2AX_EnhancerOnly_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/D7_H3K27me3_EnhancersAndPromoterEnhancers_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/D7_H3K27me3_EnhancerOnly_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/D7_H3K27ac_EnhancersAndPromoterEnhancers_0_GHaligned_EnsemblEntrez.txt",
    "~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/D7_H3K27ac_EnhancerOnly_0_GHaligned_EnsemblEntrez.txt"
  ))
geneDFList <- list()
for (gene in genesOfInterest)
{
  for (GHtype in GHtypeList)
  {
    if(GHtype=="EnhancerOnly"){fileList <- files[!(grepl("Promoter", files))]} else{fileList <- files[(grepl("Promoter", files))]}
    # number interacting via promoter. number interacting via enhancer
    geneDF <- data.frame()
    for (file in fileList)
    {
      histName <- sampleList[str_detect(file, sampleList)]
      bed2 <- read.table(file = file, sep="\t", header=T, quote="")
      colnames(bed2) <- sapply(colnames(bed2), mcolsRemove)
      bed2 <- bed2[,match(wantedCols, colnames(bed2))]
      # subetting bed
      
      bed <- bed2[((bed2$interactedGeneName == gene)& bed2$distToGenehancer <200)|(bed2$GENENAME == gene& bed2$distanceToTSS < 3000),]
      # assigning cols
      if(length(unique(bed$peakNumber[(bed$GENENAME == gene& bed$distanceToTSS < 3000)])))
      {
        geneDF["PeakNumViaPromoter",histName] <- length(unique(bed$peakNumber[(bed$GENENAME == gene& bed$distanceToTSS < 3000)]))
        geneDF["SumPeakScorePromoter",histName] <- 
          sum(aggregate(bed[((bed$GENENAME == gene)& bed$distanceToTSS < 3000),
                            -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                        list(bed[(bed$GENENAME == gene)& bed$distanceToTSS < 3000,]$peakNumber), mean)$score)
        geneDF["AvgPeakScorePromoter",histName] <- 
          mean(aggregate(bed[((bed$GENENAME == gene)& bed$distanceToTSS < 3000),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[(bed$GENENAME == gene)& bed$distanceToTSS < 3000,]$peakNumber), mean)$score)
        geneDF["PropPeakScorePromoter",histName] <- 
          mean(aggregate(bed[((bed$GENENAME == gene)& bed$distanceToTSS < 3000),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[(bed$GENENAME == gene)& bed$distanceToTSS < 3000,]$peakNumber), mean)$score)/quantile(bed$score, probs=0.999)
        geneDF["SumPeakTotalSignalPromoter",histName] <- 
          sum(aggregate(bed[((bed$GENENAME == gene)& bed$distanceToTSS < 3000),
                            -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                        list(bed[(bed$GENENAME == gene)& bed$distanceToTSS < 3000,]$peakNumber), mean)$TotalSignal)
        geneDF["AvgPeakTotalSignalPromoter",histName] <- 
          mean(aggregate(bed[((bed$GENENAME == gene)& bed$distanceToTSS < 3000),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[(bed$GENENAME == gene)& bed$distanceToTSS < 3000,]$peakNumber), mean)$TotalSignal)
        geneDF["PropPeakTotalSignalPromoter",histName] <- 
          mean(aggregate(bed[((bed$GENENAME == gene)& bed$distanceToTSS < 3000),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[(bed$GENENAME == gene)& bed$distanceToTSS < 3000,]$peakNumber), mean)$TotalSignal)/quantile(bed$TotalSignal, probs=0.999)
      }
      if(length(unique(bed$peakNumber[(bed$interactedGeneName == gene) & bed$distToGenehancer <200])))
      {
        geneDF["PeakNumViaEnhancer",histName] <- length(unique(bed$peakNumber[(bed$interactedGeneName == gene) & bed$distToGenehancer <200]))
        geneDF["SumPeakScoreEnhancer",histName] <- 
          sum(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),
                            -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                        list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$score)
        geneDF["AvgPeakScoreEnhancer",histName] <- 
          mean(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$score)
        geneDF["PropPeakScoreEnhancer",histName] <- 
          mean(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$score)/quantile(bed$score, probs=0.999)
        geneDF["SumPeakTotalSignalEnhancer",histName]<- 
          sum(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),
                            -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                        list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$TotalSignal)
        geneDF["AvgPeakTotalSignalEnhancer",histName]<- 
          mean(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$TotalSignal)
        geneDF["PropPeakTotalSignalEnhancer",histName]<- 
          mean(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),
                             -which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))],
                         list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$TotalSignal)/quantile(bed$TotalSignal, probs=0.999)
        geneDF["SumPeakInteractionValue",histName] <-
          sum(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),-which(names(bed) %in% c("interactedGeneName",
                                                                                                                   "GENENAME", "peakNumber"))],list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$interactionValue)
        geneDF["AvgPeakInteractionValue",histName] <-
          mean(aggregate(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),-which(names(bed) %in% c("interactedGeneName",
                                                                                                                    "GENENAME", "peakNumber"))],list(bed[((bed$interactedGeneName == gene)& bed$distToGenehancer <200),]$peakNumber), mean)$interactionValue)
      }
      if((length(unique(bed$peakNumber[(bed$interactedGeneName == gene) & bed$distToGenehancer <200])))&(length(unique(bed$peakNumber[(bed$GENENAME == gene& bed$distanceToTSS < 3000)]))))
      {
        geneDF["SumPeakScoreBoth",histName] <-
          sum(aggregate(bed[,-which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))], list(bed$peakNumber), mean)$score)
        geneDF["AvgPeakScoreBoth",histName] <-
          mean(aggregate(bed[,-which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))], list(bed$peakNumber), mean)$score)
        geneDF["PropPeakScoreBoth",histName] <-
          mean(aggregate(bed[,-which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))], list(bed$peakNumber), mean)$score)/quantile(bed$score, probs=0.999)
        geneDF["SumPeakTotalSignalBoth",histName]<- 
          sum(aggregate(bed[,-which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))], list(bed$peakNumber), mean)$TotalSignal)
        geneDF["AvgPeakTotalSignalBoth",histName]<- 
          mean(aggregate(bed[,-which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))], list(bed$peakNumber), mean)$TotalSignal)
        geneDF["PropPeakTotalSignalBoth",histName]<- 
          mean(aggregate(bed[,-which(names(bed) %in% c("interactedGeneName", "GENENAME", "peakNumber"))], list(bed$peakNumber), mean)$TotalSignal)/quantile(bed$TotalSignal, probs=0.999)
      }
      
      
      geneDF
    }
    geneDFList[[gene]][[GHtype]] <- geneDF
    if(!dir.exists(paste0(projPath, "/curatedData/GHalignedGenesOfInterest/", gene, "/")))
    {
      mkdirs(paste0(projPath, "/curatedData/GHalignedGenesOfInterest/", gene, "/"))
    } 
    write.table(geneDFList[[gene]][[GHtype]], file=paste0(projPath, "/curatedData/GHalignedGenesOfInterest/", gene, "/" ,gene, "_sample_GHaligned_DataComparison_", GHtype , ".txt" ), sep = "\t", col.names=TRUE, row.names=TRUE, quote = FALSE)
  }
}
saveRDS(geneDFList, file= paste0(projPath, "/curatedData/GHalignedGenesOfInterest/totalcomparison.RDS"))
