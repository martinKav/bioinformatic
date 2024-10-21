# source file for annotation and enrichment of a given peak file
library(dplyr)
library(stringr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(rentrez)
library(rGREAT)
library(R.utils)
# ***********************************************
annotateAndEnrich <- function(fileName, Quant = 0, greatEnrich = FALSE){
  if(!dir.exists(file.path(dirname(paste0(str_split(fileName, "peakCalling")[[1]][1],"annotationsAndEnrichment",str_split(str_split(fileName,"peakCalling")[[1]][2],".bed")[[1]][1], ".enrichment.txt")))))
  {
    mkdirs( file.path(dirname(paste0(str_split(fileName,
       "peakCalling")[[1]][1],"annotationsAndEnrichment",str_split(str_split( fileName, "peakCalling")[[1]][2],".bed")[[1]][1], ".enrichment.txt"))))
  }
  
  bed <- read.table(file = fileName, sep = "\t")
  bed
  gr <- GRanges(seqnames = bed$V1, ranges= IRanges(start=bed$V2, end=bed$V3),
                score=bed$V5, TotalSignal=bed$V4, MaxSignalRegion=bed$V6)
  gr$peakNumber <- c(paste0("Peak",1:length(gr)))
  gr <- gr[score(gr) >= quantile(score(gr), probs= Quant)]
  annotation <- annotatePeak(gr, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                             tssRegion = c(-1000, 1000))
  peakAnnotation <- as.data.frame(annotation)
  annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                           keys = peakAnnotation$geneId,
                                           columns = c("GENENAME"),
                                           keytype = "ENTREZID")
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  annotatedPeaks <- left_join(peakAnnotation, annotations_edb, by=c("geneId"="ENTREZID"), relationship = "many-to-many")
  if(any(is.na(annotatedPeaks$GENENAME)))
  {
    unidentified_geneIds <- data.frame(geneID = unique(annotatedPeaks[is.na(annotatedPeaks$GENENAME),]$geneId),
                                       scrapedGeneName = NA)
    
    binSize=200
    print(paste0("Will take ", ceiling(length(unidentified_geneIds$geneID)/binSize), " iterations"))
    for(start in 1:ceiling(length(unidentified_geneIds$geneID)/binSize))
    {
      print(paste0("iteration ", start))
      if(((length(unidentified_geneIds$geneID)/binSize)/start) <1){last=((start-1)*binSize)+length(unidentified_geneIds$geneID)%%binSize}else{last = start*binSize}
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
    betterannotatedPeaks <- left_join(peakAnnotation, totalAnnos, by=c("geneId"="ENTREZID"), relationship = "many-to-many")
    unique(betterannotatedPeaks$annotation)
    
    write.table(betterannotatedPeaks, 
                file=paste0(str_split(fileName, "peakCalling")[[1]][1], "annotationsAndEnrichment", str_split(str_split(fileName, "peakCalling")[[1]][2], ".bed")[[1]][1], Quant, ".annotation.txt"), sep="\t", quote=F, eol="\n") 
    
  }else{
  write.table(annotatedPeaks, 
              file=paste0(str_split(fileName, "peakCalling")[[1]][1], "annotationsAndEnrichment", str_split(str_split(fileName, "peakCalling")[[1]][2], ".bed")[[1]][1], Quant, ".annotation.txt"), sep="\t", quote=F, eol="\n")  
  }
  if(greatEnrich) # this is FALSE by default because we enrich the overlaps anyway
  {
    #enrichment with rGREAT
    gr2<-gr[score(gr) %in% sort(score(gr),decreasing = T)[1:5000]]
    score(gr2) <- as.integer(score(gr2))
    job <- submitGreatJob(gr2, genome = "hg38")
    enrichmentTable <- getEnrichmentTables(job, ontology = "GO Biological Process", category = "GO",
                                           request_interval = 10, max_tries = 100, download_by = c("json", "tsv"),
                                           verbose = TRUE)
    write.table(enrichmentTable$`GO Biological Process`, 
                file = paste0(str_split(fileName, "peakCalling")[[1]][1], "annotationsAndEnrichment", str_split(str_split(fileName, "peakCalling")[[1]][2], ".bed")[[1]][1], Quant, ".enrichment.txt"), sep="\t", quote=F, eol="\n") 
  }
}
