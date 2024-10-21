# setup
rm(list=ls())
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(rentrez)
library(rGREAT)
setwd("C:/Users/mjmka/Documents/Bioinformatics")
projPath="C:/Users/mjmka/Documents/Bioinformatics"
# ***********************************************
# REMOVING NAMES
scoreThresholdList = c("top10percent", "top25percent","top50percent")
for (scoreThreshold in scoreThresholdList)
{  
  GR_list <- readRDS(file = paste0(projPath, "/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_",scoreThreshold,"_controlgone_Filtered.peaks.GRangesList"))
  # fix names of each
  gr_vector=list()
  for( i in 1:length(GR_list)){
    name = names(GR_list[i])
    name
    gr <- unlist(GR_list[i])
    names(gr) <- NULL
    gr_vector[[name]] <- gr
  }
  GR_listFIX <- GRangesList(gr_vector)
  saveRDS(GR_listFIX, file=paste0(projPath, "/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_",scoreThreshold,"_controlgone_FIXED.peaks.GRangesList"))
}
# enriching and annotating

# annotate based on nearest peak
# read bed file
# then annotatePeak
# then give those IDs gene annotations
# then do a great analysis - GREAT requires GRanges, so we might as well convert to GRanges at the start
# then store in folder


scoreThreshold="top50percent"
#for (scoreThreshold in scoreThresholdList)
{
  GR_listFIX <- readRDS(file=paste0(projPath, "/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_",scoreThreshold,"_Filtered.peaks.GRangesList"))
  peakAnnoList <- lapply(GR_listFIX, annotatePeak, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                       tssRegion=c(-7000, 7000))
  # plotAnnoBar(peakAnnoList)
  # plotDistToTSS(peakAnnoList, title="Distribution of histone modifications \n relative to TSS")
  # 
  # genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
  # 
  # # Run KEGG analysis
  # compKEGG <- compareCluster(geneCluster = genes, 
  #                            fun = "enrichKEGG",
  #                            organism = "human",
  #                            pvalueCutoff  = 0.05, 
  #                            pAdjustMethod = "BH")
  # dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")
  # this is all cool but not very useful
  # lets just get enrichment and skidaddle
  for (histName in names(GR_listFIX))
  {
    peak <- data.frame(peakAnnoList[[histName]]@anno)
    entrez <- peak$geneId
    annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                             keys = entrez,
                                             columns = c("GENENAME"),
                                             keytype = "ENTREZID")
    annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
    
    
    annotatedPeaks <- left_join(peak, annotations_edb, by=c("geneId"="ENTREZID"), relationship = "many-to-many")
    unidentified_geneIds <- data.frame(geneID = unique(annotatedPeaks[is.na(annotatedPeaks$GENENAME),]$geneId),
                                       scrapedGeneName = NA)
    binSize=200
    for(start in 1:round(length(unidentified_geneIds$geneID)/binSize,digits=0))
    {
      if(start ==round(length(unidentified_geneIds$geneID)/binSize,digits=0)){last=((start-1)*binSize)+length(unidentified_geneIds$geneID)%%binSize}else{last = start*binSize}
      upload <- entrez_post(db="gene", id=c(unidentified_geneIds$geneID[(1 + binSize*(start-1)) : (last)]))
      bin <- entrez_summary(db="gene",id=c(unidentified_geneIds$geneID[(1 + binSize*(start-1)) : (last)]), use_history=TRUE)
      for (i in 1:length(bin))
      {
        unidentified_geneIds[unidentified_geneIds$geneID == bin[[i]]$uid,]$scrapedGeneName <- bin[[i]]$name
      }
      Sys.sleep(10)
    }
    colnames(unidentified_geneIds)  <- c("ENTREZID",  "GENENAME")   
    totalAnnos <- rbind(annotations_edb, unidentified_geneIds)
    betterannotatedPeaks <- left_join(peak, totalAnnos, by=c("geneId"="ENTREZID"), relationship = "many-to-many")
    
    write.table(betterannotatedPeaks, file=paste0(projPath, "/R_beds/CUTandTAG/annotation/", histName, "_seacr_relaxed_", scoreThreshold,"_controlgone_annotation.txt"), sep="\t", quote=F, row.names=F)
    
    # enrichment with rGREAT
    # gr <- GR_listFIX[[histName]]
    # gr2<-gr[score(gr) %in% sort(score(gr),decreasing = T)[1:5000]]
    # score(gr2) <- as.integer(score(gr2))
    # job <- submitGreatJob(gr2, genome = "hg38")
    # enrichmentTable <- getEnrichmentTables(job, ontology = "GO Biological Process", category = "GO",
    #                                        request_interval = 10, max_tries = 100, download_by = c("json", "tsv"),
    #                                        verbose = TRUE)
    # write.table(enrichmentTable$`GO Biological Process`, file = paste0(projPath, "/R_beds/CUTandTAG/enrichment/", histName, "_seacr_relaxed_controlgone_top5000enrichment.txt"), sep=" ")
    
    }
}



betterannotatedPeaks






peak <- GR_listFIX[[histName]]
plotPeakProf2(peak = peak, upstream = rel(1), downstream = rel(1),
              conf = 0.95, by = "gene", type = "body", nbin = 8000,
              TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, weightCol = "score",ignore_strand = F)

annotatedGR <- annotatePeak(peak, tssRegion = c(-3000,3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
plotAnnoBar(annotatedGR)
plotDistToTSS(annotatedGR, title="Distribution of IPSC_yH2AX peaks relative to TSS")
ego <- enrichGO(gene = entrez, 
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

clusterResults <- data.frame(ego)
clusterResults
dotplot(ego, showCategory=30)

ekegg <- enrichKEGG(gene = entrez,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

dotplot(ekegg)


