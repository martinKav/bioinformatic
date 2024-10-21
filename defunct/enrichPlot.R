#R script dedicated to Enrichment
library(dplyr)
library(stringr)
library(readr)
library(clusterProfiler)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX", "H3K4me1")
GHtypeList=c("EnhancerOnly", "EnhancersAndPromoterEnhancers")
quantList=c(0,0.5,0.75,0.9,0.99)
# for each hist and quant level, we want to subset scores and search for entrez GO with ego and ekegg
# we need to tinker with graphs though so I will go 1 by 1 
fileNames <- findInDirList(paste0(projPath, "/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated"), ".txt")
for (fileName in fileNames)
{
  histName <- sampleList[str_detect(fileName, sampleList)]
  
  GHtype <- GHtypeList[str_detect(fileName, GHtypeList)]
  if(!dir.exists(paste0(projPath, "/plots/Enrichment/GHalign/", histName)))
  {
    mkdirs(paste0(projPath, "/plots/Enrichment/GHalign/", histName))
  } 
  DF <- read.table(file=fileName, header = T, sep="\t", quote="")
  for (quant in quantList)
  {
    ENTREZIDs <- DF$interactedGeneENTREZ[DF$mcols.score > quantile(DF$mcols.score, probs=quant)]
    ENTREZIDs <- ENTREZIDs[!(is.na(ENTREZIDs))]
    { # enrichGO dotplot
      png(file=paste0(projPath, "/plots/Enrichment/GHalign/", histName, "/", histName, "_Enrichment_Genehancer_", GHtype,"_", quant, "quantile_enrichGOdotplot.png"), width = 1920,
          height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
      ego <- enrichGO(gene = ENTREZIDs, 
                      keyType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)
      clusterResults <- data.frame(ego)
      dotplot <- dotplot(ego, showCategory=30, title=paste0("enrichGO Enrichment Analysis for ", histName, "(Quantile=", quant, ".GHtype= ", GHtype))
      print(dotplot)
      dev.off()
    } 
    { # KEGG dotplot
      
      ekegg <- enrichKEGG(gene = ENTREZIDs,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)
      if(any(ekegg@result$qvalue < 0.05))
      {
        png(file=paste0(projPath, "/plots/Enrichment/GHalign/", histName, "/", histName, "_Enrichment_Genehancer_", GHtype,"_", quant, "quantile_KEGGdotplot.png"), width = 1920,
            height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
      dotplot <- dotplot(ekegg, showCategory=30, title=paste0("KEGG Enrichment Analysis for ", histName, "(Quantile=", quant, ".GHtype= ", GHtype))
      print(dotplot)
      dev.off()
      }
    }
  }
}
