## This code use the GSEA method to check the CCMarker grab DGEA

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)
  
##### Version information ######
  # platform       x86_64-w64-mingw32          
  # arch           x86_64                      
  # os             mingw32                     
  # system         x86_64, mingw32             
  # status                                     
  # major          4                           
  # minor          1.2                         
  # year           2021                        
  # month          11                          
  # day            01                          
  # svn rev        81115                       
  # language       R                           
  # version.string R version 4.1.2 (2021-11-01)
  # nickname       Bird Hippie  
  
##### Current path and new folder setting ##### 
  PathName = setwd(getwd())
  RVersion = "20220114_GSEA_PBMC"
  dir.create(paste0(PathName,"/",RVersion))

# ##### Current path and new folder setting  ##### 
#   Date = "20220128"
#   SampleType = "PBMC"
#   Save.Path = paste0(getwd(),"/",Date,"_",SampleType)
#   
#   dir.create(Save.Path)
#   
#   Import.Path = ""
#   RawData.Path = ""
##### Load Packages ##### 
  library(Seurat)
  library(SeuratDisk)
  library(SeuratObject)
  library(magrittr)
  library(fgsea)

##### Function setting ##### 
  ## Call function
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_GSEA_LargeGeneSet.R")

##### Catch the Metadata #####
  PBMC.combined.Male <- PBMC.combined[ ,PBMC.combined@meta.data[["Sex"]] %in% c("Male")]
  PBMC.combined.Male.EO <- PBMC.combined[ ,PBMC.combined@meta.data[["Sex"]] %in% c("Male") & PBMC.combined@meta.data[["Cachexia"]] %in% c("EO")]
  # library(VGAM)
  PBMC.combined.Male.EO.df <- as.data.frame(PBMC.combined.Male.EO@assays[["RNA"]]@data)
  # data <- as(as.matrix(PBMC.combined.Male.EO@assays[["RNA"]]@data), 'sparseMatrix')


##### RNA-seq analysis in R #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(fgsea)
  
  # CellType.list <- as.character(unique(PBMC.combined@meta.data[["celltype"]]))
  CellType.list <- intersect_CellType

##### GSEAforCheck_M.Pos #####
  GSEAforCheck_M.Pos <- list()
  GSEAforCheck_M.Pos.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEAforCheck_M.Pos.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  
  Venn_Cachexia.Marker_Pos <- Venn_Cachexia.Marker_Pos[!unlist(lapply(Venn_Cachexia.Marker_Pos,is.null))]
  Venn_Cachexia.Marker_Neg <- Venn_Cachexia.Marker_Neg[!unlist(lapply(Venn_Cachexia.Marker_Neg,is.null))]
  
  for(i in 1:length(CellType.list)){
    Male.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_A"]]
    #  Male.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_A"]]
    #  Female.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_B"]]
    Female.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_B"]]
    Intersect.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Intersect_AB"]]
    #  Intersect.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Intersect_AB"]]
    
  gseaDat <- Cachexia.Marker.Male[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
  gseaDat <- data.frame(row.names(gseaDat),gseaDat)
  colnames(gseaDat)[[1]] <- c("Gene")
  ranks <- gseaDat$avg_log2FC
  names(ranks) <- gseaDat$Gene
  head(ranks)
  
  # barplot(sort(ranks, decreasing = T))
  pathwaysH <- list(union(Male.Gene.List.Pos,Intersect.Gene.List.Pos) %>% setdiff(.,Female.Gene.List.Neg))
  # pathwaysH <- list(Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_A"]])
  # pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Union_AB"]])
  names(pathwaysH) <- c("GeneforCheck")
  fgseaRes <- fgsea(pathwaysH, ranks, minSize=2, maxSize = 500)
  # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
  head(fgseaRes[order(padj, -abs(NES)), ], n=10)
  
  plotEnrichmen <- plotEnrichment(pathwaysH[["GeneforCheck"]], ranks)
  # plotEnrichmen
  Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichmen)
  names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichmen")
  GSEAforCheck_M.Pos[[i]] <- Sum
  names(GSEAforCheck_M.Pos)[[i]] <- paste0("Male.",CellType.list[i])
  
  fgseaRes2 <- data.frame("Pos",paste0("M.",CellType.list[i]),fgseaRes)
  colnames(fgseaRes2)[[1]] <- c("GeneType")
  colnames(fgseaRes2)[[2]] <- c("PhenoType")
  GSEAforCheck_M.Pos.df <- rbind(GSEAforCheck_M.Pos.df,fgseaRes2 )
  
  rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,
     Male.Gene.List.Pos,Male.Gene.List.Neg,Female.Gene.List.Pos,Female.Gene.List.Neg,Intersect.Gene.List.Pos,Intersect.Gene.List.Neg)
  }
  rm(i)

##### GSEAforCheck_M.Neg #####
  GSEAforCheck_M.Neg <- list()
  GSEAforCheck_M.Neg.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEAforCheck_M.Neg.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  
  for(i in 1:length(CellType.list)){
    #  Male.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_A"]]
    Male.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_A"]]
    Female.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_B"]]
    #  Female.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_B"]]
    #  Intersect.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Intersect_AB"]]
    Intersect.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Intersect_AB"]]
    
    gseaDat <- Cachexia.Marker.Male[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    head(ranks)
    
    # barplot(sort(ranks, decreasing = T))
    pathwaysH <- list(union(Male.Gene.List.Neg,Intersect.Gene.List.Neg) %>% setdiff(.,Female.Gene.List.Pos))
    # pathwaysH <- list(Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_A"]])
    # pathwaysHTTT <- list(Venn_Cachexia.Marker_Neg[["Venn_Cachexia.Marker.Mac_Neg"]][["Summary"]][["Union_AB"]])
    names(pathwaysH) <- c("GeneforCheck")
    fgseaRes <- fgsea(pathwaysH, ranks, minSize=2, maxSize = 500)
    # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
    head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    plotEnrichmen <- plotEnrichment(pathwaysH[["GeneforCheck"]], ranks)
    # plotEnrichmen
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichmen)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichmen")
    GSEAforCheck_M.Neg[[i]] <- Sum
    names(GSEAforCheck_M.Neg)[[i]] <- paste0("Male.",CellType.list[i])
    
    fgseaRes2 <- data.frame("Neg",paste0("M.",CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("GeneType")
    colnames(fgseaRes2)[[2]] <- c("PhenoType")
    GSEAforCheck_M.Neg.df <- rbind(GSEAforCheck_M.Neg.df,fgseaRes2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,
       Male.Gene.List.Pos,Male.Gene.List.Neg,Female.Gene.List.Pos,Female.Gene.List.Neg,Intersect.Gene.List.Pos,Intersect.Gene.List.Neg)
  }
  rm(i)

##### GSEAforCheck_F.Pos #####
  GSEAforCheck_F.Pos <- list()
  GSEAforCheck_F.Pos.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEAforCheck_F.Pos.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  for(i in 1:length(CellType.list)){
    #  Male.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_A"]]
    Male.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_A"]]
    Female.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_B"]]
    #  Female.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_B"]]
    Intersect.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Intersect_AB"]]
    #  Intersect.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Intersect_AB"]]
    
    gseaDat <- Cachexia.Marker.Female[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    head(ranks)
    
    # barplot(sort(ranks, decreasing = T))
    pathwaysH <- list(union(Female.Gene.List.Pos,Intersect.Gene.List.Pos) %>% setdiff(.,Male.Gene.List.Neg))
    # pathwaysH <- list(Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_A"]])
    # pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Union_AB"]])
    names(pathwaysH) <- c("GeneforCheck")
    fgseaRes <- fgsea(pathwaysH, ranks, minSize=2, maxSize = 500)
    # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
    head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    plotEnrichmen <- plotEnrichment(pathwaysH[["GeneforCheck"]], ranks)
    # plotEnrichmen
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichmen)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichmen")
    GSEAforCheck_F.Pos[[i]] <- Sum
    names(GSEAforCheck_F.Pos)[[i]] <- paste0("Female.",CellType.list[i])
    
    fgseaRes2 <- data.frame("Pos",paste0("F.",CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("GeneType")
    colnames(fgseaRes2)[[2]] <- c("PhenoType")
    GSEAforCheck_F.Pos.df <- rbind(GSEAforCheck_F.Pos.df,fgseaRes2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,
       Male.Gene.List.Pos,Male.Gene.List.Neg,Female.Gene.List.Pos,Female.Gene.List.Neg,Intersect.Gene.List.Pos,Intersect.Gene.List.Neg)
  }
  rm(i)


##### GSEAforCheck_F.Neg #####
  GSEAforCheck_F.Neg <- list()
  GSEAforCheck_F.Neg.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEAforCheck_F.Neg.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  for(i in 1:length(CellType.list)){
  try({
    Male.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_A"]]
    #  Male.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_A"]]
    #  Female.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Unique_B"]]
    Female.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_B"]]
    #  Intersect.Gene.List.Pos <- Venn_Cachexia.Marker_Pos[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Pos")]][["Summary"]][["Intersect_AB"]]
    Intersect.Gene.List.Neg <- Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Intersect_AB"]]
    
    gseaDat <- Cachexia.Marker.Female[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    head(ranks)
    
    # barplot(sort(ranks, decreasing = T))
    pathwaysH <- list(union(Female.Gene.List.Neg,Intersect.Gene.List.Neg) %>% setdiff(.,Male.Gene.List.Pos))
    # pathwaysH <- list(Venn_Cachexia.Marker_Neg[[paste0("Venn_Cachexia.Marker.",CellType.list[i],"_Neg")]][["Summary"]][["Unique_A"]])
    # pathwaysHTTT <- list(Venn_Cachexia.Marker_Neg[["Venn_Cachexia.Marker.Mac_Neg"]][["Summary"]][["Union_AB"]])
    names(pathwaysH) <- c("GeneforCheck")
    fgseaRes <- fgsea(pathwaysH, ranks, minSize=2, maxSize = 500)
    # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
    head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    plotEnrichmen <- plotEnrichment(pathwaysH[["GeneforCheck"]], ranks)
    # plotEnrichmen
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichmen)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichmen")
    GSEAforCheck_F.Neg[[i]] <- Sum
    names(GSEAforCheck_F.Neg)[[i]] <- paste0("Female.",CellType.list[i])
    
    fgseaRes2 <- data.frame("Neg",paste0("F.",CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("GeneType")
    colnames(fgseaRes2)[[2]] <- c("PhenoType")
    GSEAforCheck_F.Neg.df <- rbind(GSEAforCheck_F.Neg.df,fgseaRes2 )
      rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,
       Male.Gene.List.Pos,Male.Gene.List.Neg,Female.Gene.List.Pos,Female.Gene.List.Neg,Intersect.Gene.List.Pos,Intersect.Gene.List.Neg)
  })
    }
  rm(i)
  
##### GSEAforCheck.Sum #####
  GSEAforCheck.Sum <- rbind(GSEAforCheck_M.Pos.df,GSEAforCheck_M.Neg.df,GSEAforCheck_F.Pos.df,GSEAforCheck_F.Neg.df)
  GSEAforCheck.Sum <- GSEAforCheck.Sum[,!colnames(GSEAforCheck.Sum) %in% c("leadingEdge")]
  write.table(GSEAforCheck.Sum, file=paste0(PathName,"/",RVersion,"/Cachexia_PBMC_GSEAforCheck2.txt"),sep="\t",
              row.names=F, quote = FALSE)



## ------------------------------------------------------------------------------------------------------------------- ##

##### RNA-seq analysis in R #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(fgsea)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
  
  
  # Geneset from GSEA
  # H.all <- read.delim(paste0(PathName,"/h.all.v7.4.symbols.gmt"),header = F)
  H.all <- read.delim2(paste0(PathName,"/GSEA_Geneset_AllIndex.txt"),
                       col.names = 1:max(count.fields(paste0(PathName,"/GSEA_Geneset_AllIndex.txt"))),
                       header = F,sep = "\t")

##### GSEA.Large_Female #####
  GSEA.Large_Female <- list()
  GSEA.Large_Female.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA.Large_Female.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA.Large_Female.df.TOP <- GSEA.Large_Female.df
  
  for(i in 1:length(CellType.list)){
    
    gseaDat <- Cachexia.Marker.Female[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))
  
    
    GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,H.all,10)
    
    fgseaRes <- GSEA.Large.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA.Large.Output[["H.all.list"]] 
    
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
    
    topPathways <- GSEA.Large.Output[["topPathways"]]
    
    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.Female.",CellType.list[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.Female.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.Female.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA.Large_Female[[i]] <- Sum
    names(GSEA.Large_Female)[[i]] <- paste0("Female.",CellType.list[i])
    
    fgseaRes2 <- data.frame("Neg",paste0("F.",CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("GeneType")
    colnames(fgseaRes2)[[2]] <- c("PhenoType")
    GSEA.Large_Female.df <- rbind(GSEA.Large_Female.df,fgseaRes2 )
    
    topPathways2 <- data.frame("Neg",paste0("F.",CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("GeneType")
    colnames(topPathways2)[[2]] <- c("PhenoType")
    GSEA.Large_Female.df.TOP <- rbind(GSEA.Large_Female.df.TOP,topPathways2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
  
  }

##### GSEA.Large_Male #####
  library(fgsea)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
  
  GSEA.Large_Male <- list()
  GSEA.Large_Male.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA.Large_Male.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA.Large_Male.df.TOP <- GSEA.Large_Male.df
  
  for(i in 1:length(CellType.list)){
    
    gseaDat <- Cachexia.Marker.Male[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))
    
    
    GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,H.all,10)
    
    fgseaRes <- GSEA.Large.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA.Large.Output[["H.all.list"]] 
    
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
    
    topPathways <- GSEA.Large.Output[["topPathways"]]
    
    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.Male.",CellType.list[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.Male.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.Male.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA.Large_Male[[i]] <- Sum
    names(GSEA.Large_Male)[[i]] <- paste0("Male.",CellType.list[i])
    
    fgseaRes2 <- data.frame("Neg",paste0("M.",CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("GeneType")
    colnames(fgseaRes2)[[2]] <- c("PhenoType")
    GSEA.Large_Male.df <- rbind(GSEA.Large_Male.df,fgseaRes2 )
    
    topPathways2 <- data.frame("Neg",paste0("M.",CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("GeneType")
    colnames(topPathways2)[[2]] <- c("PhenoType")
    GSEA.Large_Male.df.TOP <- rbind(GSEA.Large_Male.df.TOP,topPathways2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
    
  }
  


##### GSEA.Large.Sum.TOP #####
  GSEA.Large.Sum.TOP <- rbind(GSEA.Large_Male.df.TOP,GSEA.Large_Female.df.TOP)
  GSEA.Large.Sum.TOP <- GSEA.Large.Sum.TOP[,!colnames(GSEA.Large.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA.Large.Sum.TOP, file=paste0(PathName,"/",RVersion,"/Cachexia_PBMC_GSEA_LargeTOP.txt"),sep="\t",
              row.names=F, quote = FALSE)


##### Bubble plot #####
  SEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP[abs(GSEA.Large.Sum.TOP$NES) > 1.7,]
  GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP[abs(GSEA.Large.Sum.TOP$padj) < 0.05,]
  # GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP[abs(GSEA.Large.Sum.TOP$padj) < 0.25,]
  # GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP.S[abs(GSEA.Large.Sum.TOP.S$pval) < 0.05,]
  library(ggplot2)
  library(scales)
  ggplot(GSEA.Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = padj, size = abs(NES))) + 
    geom_point() +
    scale_color_distiller()
  
  ggplot(GSEA.Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                           guide = "colourbar",midpoint = 0)
  
  BBPlot <- ggplot(GSEA.Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
            geom_point() +
            scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d", 
            guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
  BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal", xangle =90,
                            XtextSize=12,  YtextSize=8, AxisTitleSize=1)



##### Test Cachexia candidate markers #####
  # load(paste0(PathName,"/mouse_H_v5.RData"))
  # pathwaysH <- Mm.H
  # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500)
  
  pathwaysHTTT <- list(gseaDat$Gene[1:150])
  pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Unique_A"]])
  pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Union_AB"]])
  pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Unique_B"]])
  pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Intersect_AB"]])
  names(pathwaysHTTT) <- c("Test")
  fgseaRes <- fgsea(pathwaysHTTT, ranks, minSize=15, maxSize = 500)
  # fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
  head(fgseaRes[order(padj, -abs(NES)), ], n=10)
  
  plotEnrichment(pathwaysHTTT[["Test"]], ranks)

##### Using fgsea package #####
  # https://bioc.ism.ac.jp/packages/3.6/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html





# ## https://www.biostars.org/p/339934/
# library(msigdbr)
# library(fgsea)
# 
# #Retrieve human H (hallmark) gene set
# msigdbr_df <- msigdbr(species = "human", category = "H")
# head(msigdbr_df)
