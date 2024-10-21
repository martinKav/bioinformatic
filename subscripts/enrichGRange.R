# this script will take a GRanges object (representing an overlap) and enrich it
library(stringr)
library(readr)
library(R.utils)
library(rGREAT)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library("org.Hs.eg.db")
library("reactome.db")
# it will need to extract the gene names for GO
# rGREAT can work off of the GRanges
# HOMER will need some setup
# BioCor clustering will also need some tinkering
getLower <- function(num1, num2){if(num1 >= num2){ return(num2)}else{return(num1)}}
enrichGRange <- function(GRlist, projPath, inQuant, outQuant, GHtype, elementsN)
{
  for (GRName in names(GRlist))
  {
    if(GHtype == "EnhancerOnly"){short <- "EO"}else{short <- "EPAP"}
    if (file.exists(paste0(projPath, "/geneOverlaps/GeneOverlaps/", N,"/GREAT/", GRName,"_",short,"_InQ",inQuant, "_OutQ", outQuant, "_GREAT.txt")) |
        file.exists(paste0(projPath, "/geneOverlaps/GeneOverlaps/", N,"/GREAT/", GRName,"_", GHtype,"_",short,"_InQ",inQuant, "_OutQ", outQuant, "_GREAT.txt"))){print("file exists!")}else{
    if(length(GRlist[[GRName]]) ==0){print(paste0(GRName,"_",short,"_InQ",inQuant, "_OutQ", outQuant, " is an empty list!"))
      cat(paste0(GRName,"_",short,"_InQ",inQuant, "_OutQ", outQuant, "\n"),file=paste0(projPath, "/geneOverlaps/GeneOverlaps/EmptyOverlaps.txt"),append=TRUE)}
    else
      {
    GR = GRlist[[GRName]]
    bed <- data.frame(GR)
    bed$association[bed$distanceToTSS < 2000] <- "Promoter"
    bed$association[bed$distToGenehancer < 200] <- "Enhancer"
    bed$association[bed$distanceToTSS < 2000 & bed$distToGenehancer < 200] <- "Both"
    bed$GenehancerID[bed$distToGenehancer > 200] <- NA
    bed$interactionValue[bed$distToGenehancer > 200] <- NA
    bed$interactedGeneName[bed$distToGenehancer > 200] <- NA
    bed$interactedGeneENSEMBL[bed$distToGenehancer > 200] <- NA
    bed$interactedGeneENTREZ[bed$distToGenehancer > 200] <- NA
    bed$GenehancerType[bed$distToGenehancer > 200] <- NA
    bed$GENENAME[bed$distanceToTSS > 2000] <- NA
    colnames(bed)[colnames(bed) == "GENENAME"] <- "PromoterAssociatedGene"
    colnames(bed)[colnames(bed) == "interactedGeneName"] <- "EnhancerAssociatedGene"
    colnames(bed)[colnames(bed) == "interactedGeneENSEMBL"] <- "EnhancerAssociatedGeneENSEMBL"
    colnames(bed)[colnames(bed) == "interactedGeneENTREZ"] <- "EnhancerAssociatedGeneENTREZ"
    colnames(bed)[colnames(bed) == "interactionValue"] <- "EnhancerGeneInteractionValue"
    colnames(bed)[colnames(bed) == "distToGenehancer"] <- "distanceToGenehancer"
    bed2 <- bed[,c("association", "PromoterAssociatedGene", "EnhancerAssociatedGene", "peakNumber","seqnames", "start", "end", "score", "TotalSignal",
                  "distanceToTSS", "distanceToGenehancer", "EnhancerGeneInteractionValue")]
    mkdirs(paste0(projPath, "/geneOverlaps/GeneOverlaps/",N, "/Condensed"))
    write.table(bed2, file=paste0(projPath, "/geneOverlaps/GeneOverlaps/", N,"/", "Condensed/", GRName,"_", short, "_InQ",inQuant, "_OutQ", outQuant, "_Condensed_Annotation.txt"), sep="\t", col.names=T, quote=FALSE, row.names=F)
    {
    print(paste0("Submitting GREAT"))
    GRe <-  GR[!duplicated(GR$peakNumber),]
    GRe<-GRe[GRe$TotalSignal %in% sort(GRe$TotalSignal,decreasing = T)[1:getLower(5000, length(GRe$TotalSignal))]]
    GRe$TotalSignal <- as.integer(GRe$TotalSignal)
    job <- submitGreatJob(GRe, genome = "hg38")
    enrichmentTable <- getEnrichmentTables(job, ontology = "GO Biological Process", category = "GO",
                                           request_interval = 10, max_tries = 100, download_by = c("json", "tsv"),
                                           verbose = TRUE)
    greatTable <- enrichmentTable$`GO Biological Process`[enrichmentTable$`GO Biological Process`$Hyper_Adjp_BH <0.05,][order(enrichmentTable$`GO Biological Process`[enrichmentTable$`GO Biological Process`$Hyper_Adjp_BH <0.05,]$Hyper_Fold_Enrichment, decreasing=T),]
    greatTable$name <- factor(greatTable$name, levels=rev(greatTable$name))
    
    greatDotPlot <- ggplot(greatTable[1:getLower(30,length(greatTable$name)),], aes(x=Hyper_Fold_Enrichment , y= name, colour= Hyper_Adjp_BH, size=Hyper_Gene_Set_Coverage))+
      geom_point()+
      scale_size_continuous(limits=c(0.001,1))+
      labs(title= paste0("GREAT Enrichment Analysis for \n", GRName, "\n(InQ=",inQuant, "OutQ=", outQuant, ".GHtype= ", GHtype ,")"))+
      theme_bw(base_size = 20)+
      theme(plot.title = element_text(size=18))
    mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Enrichment/", N, "/GREAT"))
    png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Enrichment/", N, "/GREAT/", GRName,"_",short, "_InQ",inQuant, "_OutQ", outQuant, "_GREAT.png"), width = 1920,
        height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
    print(greatDotPlot)
    while (!is.null(dev.list()))  dev.off()
    mkdirs(paste0(projPath, "/geneOverlaps/GeneOverlaps/",N, "/GREAT"))
    write.table(enrichmentTable$`GO Biological Process`, 
                file = paste0(projPath, "/geneOverlaps/GeneOverlaps/", N,"/GREAT/", GRName,"_",short, "_InQ",inQuant, "_OutQ", outQuant, "_GREAT.txt"), sep="\t", quote=F, eol="\n") 
    }
    
    #getting entrez
    GR$interactedGeneName[GR$distToGenehancer > 200] <- NA
    GR$GENENAME[GR$distanceToTSS > 2000] <- NA
    # getting straggler ENTREZ
    symbols <- unique(c(GR$GENENAME[!is.na(GR$GENENAME)],GR$interactedGeneName[!is.na(GR$interactedGeneName)]))
    if (length(symbols) <= 1){ print("only 1 gene or empty")}else{
    annotations_edb <- AnnotationDbi::select(org.Hs.eg.db, 
                             keys = symbols,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")
    
    
      ENTREZIDs <- as.character(annotations_edb$ENTREZID)
      ENTREZIDs <- ENTREZIDs[!(is.na(ENTREZIDs))]
      ENTREZIDs
    {
      { # enrichGO dotplot
        mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Enrichment/", N, "/enrichGO"))
        
        ego <- enrichGO(gene = ENTREZIDs, 
                        keyType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
        clusterResults <- data.frame(ego)
        if(nrow(clusterResults) != 0)
          {
          png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Enrichment/", N, "/enrichGO/", GRName,"_",short, "_InQ",inQuant, "_OutQ", outQuant, "_eGO.png"), width = 1920,
              height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
        dotplot <- dotplot(ego, showCategory=30, title=paste0("enrichGO Enrichment Analysis for ", GRName, "\n(InQ=",inQuant, "OutQ=", outQuant, ".GHtype= ", GHtype ,")"))
        print(dotplot)
        while (!is.null(dev.list()))  dev.off()
          }else{print(paste0(GRName, " empty eGO"))}
      } 
      { # KEGG dotplot
        mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Enrichment/", N, "/KEGG"))
        ekegg <- enrichKEGG(gene = ENTREZIDs,
                            organism = 'hsa',
                            pvalueCutoff = 0.05)
        if(is.null(ekegg)){print((paste0(GRName,"null ekegg")))}else
        if(any(ekegg@result$qvalue < 0.05) | (any(ekegg@result$p.adjust < 0.05) && anyNA(ekegg@result$qvalue)))
        {
          dotplot <- dotplot(ekegg, showCategory=30, title=paste0("KEGG Enrichment Analysis for ", GRName, "\n(InQ=",inQuant, "OutQ=", outQuant, ".GHtype= ", GHtype))
          if(nrow(dotplot$data))
            {
          png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Enrichment/", N, "/KEGG/", GRName,"_",short, "_InQ",inQuant, "_OutQ", outQuant, "_KEGG.png"), width = 1920,
              height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
          print(dotplot)
          while (!is.null(dev.list()))  dev.off()
            } else {print("Issue with kegg dotplot data")}
        }else{(print(paste0(GRName, "Empty Kegg")))}
      }
    }
    }
    }
    #clustering
    # genesKegg <- as.list(org.Hs.egPATH)
    # genesReact <- as.list(reactomeEXTID2PATHID)
    # genesReact
    # 
    # genesReact <- lapply(genesReact, function(x){
    #   unique(grep("R-HSA-", x, value = TRUE))
    # })
    # genesReact <- genesReact[lengths(genesReact) >= 1] 
    # GR[GR$TotalSignal %in%sort(GR$TotalSignal, decreasing = T)][1:getLower(5000, length(GR$TotalSignal))]$GENENAME
    # GRe[GRe$TotalSignal %in% sort(GRe$TotalSignal,decreasing = T)[1:getLower(5000, length(GRe$TotalSignal))]]
    # annotations_edb$SYMBOL %in% GR$GENENAME
    # annotations_edb[annotations_edb$SYMBOL %in% GR[GR$TotalSignal %in% sort(GR$TotalSignal, decreasing = T)][1:getLower(200, length(GR$TotalSignal))]$GENENAME,]
    # top500ENTREZ <- annotations_edb[annotations_edb$SYMBOL %in% GR[GR$TotalSignal %in% sort(GR$TotalSignal, decreasing = T)][1:getLower(100, length(GR$TotalSignal))]$GENENAME |
    #   annotations_edb$SYMBOL %in% GR[GR$TotalSignal %in% sort(GR$TotalSignal, decreasing = T)][1:getLower(100, length(GR$TotalSignal))]$interactedGeneName,]$ENTREZID
    # # clustering example - ok this is pretty portable
    # genes.id <- as.character(unique(top500ENTREZ))
    # genes.id <- mapIds(org.Hs.eg.db, keys = genes.id, keytype = "ENTREZID", 
    #                    column = "SYMBOL")
    # ## 'select()' returned 1:1 mapping between keys and columns
    # genes <- names(genes.id)
    # names(genes) <- genes.id
    # react <- mgeneSim(genes, genesReact)
    # kegg <- mgeneSim(genes, genesKegg)
    # ## We remove genes which are not in list (hence the warning):
    # nan <- genes %in% names(genesReact)
    # react <- react[nan, nan]
    # react
    # hc <- hclust(as.dist(1 - react))
    # hc
    #   
    #   
    # png(file=paste0(projPath, "/plots/Annotation/Genehancer/GeneOverlap/Enrichment/", N, "/REACT/Clustering/", GRName,"_",short,"_",quant, "_GREAT.png"), width = 1920,
    #     height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
    # plot(hc, main = "Similarities between genes - REACT") # shows relationship between genes
    # # cluster anal
    # dev.off()
    # mycl <- cutree(hc, h = 0.2)
    # mycl
    # clusters <- split(genes[nan], as.factor(mycl)) # puts them into a list
    # clusters
    # (clusters <- clusters[lengths(clusters) >= 2]) # removes clusters with only 1 gene
    # names(clusters) <- paste0("cluster", names(clusters)) #names them cluster1:8
    # clusters
    # sim_clus1 <- mclusterSim(clusters, genesReact) # then clusters the clusters
    # plot(hclust(as.dist(1 - sim_clus1)), 
    #      main = "Similarities between clusters by pathways") # shows similarities between clusters by pathways
    # sim_clus2 <- mclusterGeneSim(clusters, genesReact)
    # plot(hclust(as.dist(1 - sim_clus2)), 
    #      main ="Similarities between clusters by genes") #shows similarities between clusters by genes
    # 
    # kegg <- kegg[rowSums(is.na(kegg)) != ncol(kegg), rowSums(is.na(kegg)) != ncol(kegg)]
    # colnames(kegg)
    # hc <- hclust(as.dist(1-kegg))
    # plot(hc, main = "Similarlities between genes - KEGG")
    # mycl <- cutree(hc, h = 0.2)
    # mycl
    # clusters <- split(genes[colnames(kegg)], as.factor(mycl)) # puts them into a list
    # (clusters <- clusters[lengths(clusters) >= 2]) # removes clusters with only 1 gene
    # names(clusters) <- paste0("cluster", names(clusters)) #names them cluster1:8
    # sim_clus1 <- mclusterSim(clusters, genesKegg) # then clusters the clusters
    # plot(hclust(as.dist(1 - sim_clus1)), 
    #      main = "Similarities between clusters by pathways") # shows similarities between clusters by pathways
    # sim_clus2 <- mclusterGeneSim(clusters, genesKegg)
    # plot(hclust(as.dist(1 - sim_clus2)), 
    #      main ="Similarities between clusters by genes") #shows similarities between clusters by genes
    
  }
  #clustering
  # HOMER
}
}
  
  
  
  
