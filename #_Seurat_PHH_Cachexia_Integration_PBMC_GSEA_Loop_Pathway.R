##### Reference #####
  ## RNA-seq analysis in R
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  # # R绘图实战|GSEA富集分析图 #https://zhuanlan.zhihu.com/p/358168557
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
  
##### Current path and new folder setting  #####
  Date = "20220128"
  SampleType = "PBMC"
  Save.Path = paste0(getwd(),"/",Date,"_",SampleType)
  dir.create(Save.Path)

  # Import.Path = ""
  # RawData.Path = ""

##### Load Packages #####   
  library(fgsea)
  library(Seurat)
  library(dplyr)
  library(ggplot2)

##### Function setting ##### 
  ## Call function    
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")

##### Load data #####    
  # Geneset from GSEA
  # Pathway.all <- read.delim(paste0(getwd(),"/h.all.v7.5.1.symbols.gmt"),header = F)
  Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset_Pathway_AllIndex.txt"),
                             col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset_Pathway_AllIndex.txt"))),
                             header = F,sep = "\t")
  #ok# Pathway.all <- read.delim(paste0(getwd(),"/GSEA_Geneset_Pathway_AllIndex.txt"),header = F)
  #ok# Pathway.all <- read.csv(paste0(getwd(),"/GSEA_Geneset_Pathway_AllIndex.txt"),header = F,sep = "\t")
  load("D:/Dropbox/##_GitHub/0-R/##_PHH_Lab/2021_Cachexia/Seurat_PHH_Cachexia_Integration_PBMC_FindCCMarkers.RData")
  
##### Converting the Human gene name to Mouse gene name ##### 
  #  Need to be optimized
  # (Method1) bind the different length of column (Cannot use rbind)
  # (Method2) Save the data as list first and than use do.call() to unlist to have dataframe

  ## (Ori method)
  Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*2))
  for (i in 1:nrow(Pathway.all)) {
    #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
    PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t() 
    colnames(PathwayN)="Temp"
    PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
    Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
  }

  Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
  colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
  
  rm(PathwayN)
  
  # assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)
  # assign(colnames(Pathway.all)[i],Pathway.all[,i])
  
  ## (Method1)
    # Refer # https://stackoverflow.com/questions/3699405/how-to-cbind-or-rbind-different-lengths-vectors-without-repeating-the-elements-o
    # How to cbind or rbind different lengths vectors without repeating the elements of the shorter vectors?
    ## Modify by Charlene: Can use in MultRow
    bind_diff <- function(x, y){
      if(ncol(x) > ncol(y)){
        len_diff <- ncol(x) - ncol(y)
        y <- data.frame(y, rep(NA, len_diff) %>% t() %>% as.data.frame())
        colnames(x) <- seq(1:ncol(x))
        colnames(y) <- seq(1:ncol(y))
      }else if(ncol(x) < ncol(y)){
        len_diff <- ncol(y) - ncol(x)
        x <- data.frame(x, rep(NA, len_diff) %>% t() %>% as.data.frame())
        colnames(x) <- seq(1:ncol(x))
        colnames(y) <- seq(1:ncol(y))
      }
      rbind(x, y) 
    }
  
    ## Converting
    for (i in 1:nrow(Pathway.all)) {
      PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()  %>% as.data.frame()
      colnames(PathwayN)="Temp"
      PathwayN <- HSsymbol2MMsymbol(PathwayN,"Temp")
      PathwayN <- PathwayN[!PathwayN$MM.symbol  == 0,]
      PathwayN <- PathwayN[!is.na(PathwayN$MM.symbol),]
      if(i==1){
        Pathway.all.MM <- unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame()
      }else{
        Pathway.all.MM <- bind_diff(Pathway.all.MM,unique(PathwayN$MM.symbol) %>% t()  %>% as.data.frame())
        # Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
      }
    }
    
    Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
    colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
    #Pathway.all.MM[Pathway.all.MM==0] <-NA
    
    rm(PathwayN)
    
    
    
##### GSEA.Large_Female #####
  library(fgsea)
  library(ggplot2)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
    
  GSEA.Large_Female <- list()
  GSEA.Large_Female.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA.Large_Female.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA.Large_Female.df.TOP <- GSEA.Large_Female.df
  
  CellType.set <- names(Cachexia.Marker.Female)
  for(i in 1:length(CellType.set)){
    
    # Extract Rank Data
    gseaDat <- Cachexia.Marker.Female[[paste0(CellType.set[i])]][["Cachexia.Marker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))

    GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
    # 共有 14 個警告 (用 warnings() 來顯示)
    # In serialize(data, node$con) : 載入時 'package:stats' 也許不能用
    # https://community.rstudio.com/t/error-when-running-parallelized-process-warning-in-serialize-package-stats-may-not-be-available-when-loading/110573/2
    
    fgseaRes <- GSEA.Large.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA.Large.Output[["Pathway.all.list"]] 
    
    # # plotEnrichment
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)

    # plotGseaTable    
    topPathways <- GSEA.Large.Output[["topPathways"]]
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.Female.",CellType.set[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.Female.",CellType.set[i],": ",as.character(topPathways[1,1])))
    # plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.Female.",CellType.set[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    # plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA.Large_Female[[i]] <- Sum
    names(GSEA.Large_Female)[[i]] <- paste0("Female.",CellType.set[i])
    
    fgseaRes2 <- data.frame(paste0("F.",CellType.set[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA.Large_Female.df <- rbind(GSEA.Large_Female.df,fgseaRes2 )
    
    topPathways2 <- data.frame(paste0("F.",CellType.set[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
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

for(i in 1:length(CellType.set)){
  
  gseaDat <- Cachexia.Marker.Male[[paste0(CellType.set[i])]][["Cachexia.Marker.All"]]
  gseaDat <- data.frame(row.names(gseaDat),gseaDat)
  colnames(gseaDat)[[1]] <- c("Gene")
  ranks <- gseaDat$avg_log2FC
  names(ranks) <- gseaDat$Gene
  # head(ranks)
  # barplot(sort(ranks, decreasing = T))
  
  
  #GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all,10)
  GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
  
  fgseaRes <- GSEA.Large.Output[["fgseaRes"]]
  # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
  
  pathwaysH <- GSEA.Large.Output[["Pathway.all.list"]] 
  
  # plot.new()
  # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
  
  topPathways <- GSEA.Large.Output[["topPathways"]]
  
  library(ggplot2)
  plot.new()
  plotGseaTable(pathwaysH[topPathways$pathway], 
                ranks, 
                fgseaRes, 
                gseaParam = 0.5) + title( paste0("PBMC.Male.",CellType.set[i]), adj = 0, line =3)
  
  plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.Male.",CellType.set[i],": ",as.character(topPathways[1,1])))
  #plotEnrichment_Pos1
  plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.Male.",CellType.set[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
  #plotEnrichment_Neg1
  
  Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
  names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
  GSEA.Large_Male[[i]] <- Sum
  names(GSEA.Large_Male)[[i]] <- paste0("Male.",CellType.set[i])
  
  fgseaRes2 <- data.frame(paste0("M.",CellType.set[i]),fgseaRes)
  colnames(fgseaRes2)[[1]] <- c("PhenoType")
  GSEA.Large_Male.df <- rbind(GSEA.Large_Male.df,fgseaRes2 )
  
  topPathways2 <- data.frame(paste0("M.",CellType.set[i]),topPathways)
  colnames(topPathways2)[[1]] <- c("PhenoType")
  GSEA.Large_Male.df.TOP <- rbind(GSEA.Large_Male.df.TOP,topPathways2 )
  
  rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
  
}



##### GSEA.Large.Sum.TOP #####
GSEA.Large.Sum.TOP <- rbind(GSEA.Large_Male.df.TOP,GSEA.Large_Female.df.TOP)
GSEA.Large.Sum.TOP <- GSEA.Large.Sum.TOP[,!colnames(GSEA.Large.Sum.TOP) %in% c("leadingEdge")]
write.table(GSEA.Large.Sum.TOP, file=paste0(PathName,"/",RVersion,"/Cachexia_PBMC_GSEA_Pathway_LargeTOP.txt"),sep="\t",
            row.names=F, quote = FALSE)


##### Bubble plot #####
GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP[abs(GSEA.Large.Sum.TOP$NES) > 1,]
GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP.S[abs(GSEA.Large.Sum.TOP.S$padj) < 0.05,]
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
BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                          XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)





##### Cell type annotation by 5 Index #####
library(ggplot2)
## M18229_GOBP_DNA_REPAIR
M18229_GOBP_DNA_REPAIR <- read.delim(paste0(PathName,"/M18229_GOBP_DNA_REPAIR.txt"),header=T, stringsAsFactors = FALSE)
M18229_GOBP_DNA_REPAIR <- data.frame(M18229_GOBP_DNA_REPAIR[-1,])
colnames(M18229_GOBP_DNA_REPAIR) <- c("M18229_GOBP_DNA_REPAIR")

M18229_GOBP_DNA_REPAIR_Mouse <- HSsymbol2MMsymbol(M18229_GOBP_DNA_REPAIR,"M18229_GOBP_DNA_REPAIR")
M18229_GOBP_DNA_REPAIR_Mouse2 <- unique(M18229_GOBP_DNA_REPAIR_Mouse[,2])
FeaturePlot(PBMC.combined, features = M18229_GOBP_DNA_REPAIR_Mouse2[1:19], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M18229_GOBP_DNA_REPAIR_Mouse2[20:34], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M18229_GOBP_DNA_REPAIR_Mouse2[35:52], min.cutoff = "q9")

FeaturePlot(PBMC.combined, features = c("Apex1","Actr2","Btg2","Bach1"), min.cutoff = "q9")


PBMC.combined <- AddModuleScore(object = PBMC.combined, features = as.data.frame(M18229_GOBP_DNA_REPAIR_Mouse2),
                                ctrl = 5, name = 'M18229_GOBP_DNA_REPAIR')
FeaturePlot(object = PBMC.combined, features = 'M18229_GOBP_DNA_REPAIR1')+
  scale_colour_gradient2(low = "#e8fff8", mid = "#00c480", high = "#19634a", 
                         guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

# FeaturePlot(object = PBMC.combined_M18229, features = 'M18229_GOBP_DNA_REPAIR1')+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

## M10934_GOBP_HISTONE_METHYLATION
M10934_GOBP_HISTONE_METHYLATION <- read.delim(paste0(PathName,"/M10934_GOBP_HISTONE_METHYLATION.txt"),header=T, stringsAsFactors = FALSE)
M10934_GOBP_HISTONE_METHYLATION <- data.frame(M10934_GOBP_HISTONE_METHYLATION[-1,])
colnames(M10934_GOBP_HISTONE_METHYLATION) <- c("M10934_GOBP_HISTONE_METHYLATION")

M10934_GOBP_HISTONE_METHYLATION_Mouse <- HSsymbol2MMsymbol(M10934_GOBP_HISTONE_METHYLATION,"M10934_GOBP_HISTONE_METHYLATION")
M10934_GOBP_HISTONE_METHYLATION_Mouse2 <- unique(M10934_GOBP_HISTONE_METHYLATION_Mouse[,2])
FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[1:18], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[18:39], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[40:56], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[57:73], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[74:92], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[93:108], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[109:126], min.cutoff = "q9")

FeaturePlot(PBMC.combined, features = M10934_GOBP_HISTONE_METHYLATION_Mouse2[1:4], min.cutoff = "q9")

PBMC.combined <- AddModuleScore(object = PBMC.combined, features = as.data.frame(M10934_GOBP_HISTONE_METHYLATION_Mouse2),
                                ctrl = 5, name = 'M10934_GOBP_HISTONE_METHYLATION')
FeaturePlot(object = PBMC.combined, features = 'M10934_GOBP_HISTONE_METHYLATION1')+
  scale_colour_gradient2(low = "#e8fff8", mid = "#00c480", high = "#19634a", 
                         guide = "colourbar",midpoint = 0, labs(fill ="Exp"))


## M10390_GOBP_RESPONSE_TO_ZINC_ION
M10390_GOBP_RESPONSE_TO_ZINC_ION <- read.delim(paste0(PathName,"/M10390_GOBP_RESPONSE_TO_ZINC_ION.txt"),header=T, stringsAsFactors = FALSE)
M10390_GOBP_RESPONSE_TO_ZINC_ION <- data.frame(M10390_GOBP_RESPONSE_TO_ZINC_ION[-1,])
colnames(M10390_GOBP_RESPONSE_TO_ZINC_ION) <- c("M10390_GOBP_RESPONSE_TO_ZINC_ION")

M10390_GOBP_RESPONSE_TO_ZINC_ION_Mouse <- HSsymbol2MMsymbol(M10390_GOBP_RESPONSE_TO_ZINC_ION,"M10390_GOBP_RESPONSE_TO_ZINC_ION")
M10390_GOBP_RESPONSE_TO_ZINC_ION_Mouse2 <- unique(M10390_GOBP_RESPONSE_TO_ZINC_ION_Mouse[,2])
FeaturePlot(PBMC.combined, features = M10390_GOBP_RESPONSE_TO_ZINC_ION_Mouse2[1:24], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M10390_GOBP_RESPONSE_TO_ZINC_ION_Mouse2[25:50], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

PBMC.combined <- AddModuleScore(object = PBMC.combined, features = as.data.frame(M10390_GOBP_RESPONSE_TO_ZINC_ION_Mouse2),
                                ctrl = 5, name = 'M10390_GOBP_RESPONSE_TO_ZINC_ION')
FeaturePlot(object = PBMC.combined, features = 'M10390_GOBP_RESPONSE_TO_ZINC_ION1')+
  scale_colour_gradient2(low = "#e8fff8", mid = "#00c480", high = "#19634a", 
                         guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

## M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- read.delim(paste0(PathName,"/M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt"),header=T, stringsAsFactors = FALSE)
M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- data.frame(M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION[-1,])
colnames(M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION) <- c("M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Mouse <- HSsymbol2MMsymbol(M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,"M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Mouse2 <- unique(M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Mouse[,2])
FeaturePlot(PBMC.combined, features = M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Mouse2[1:20], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Mouse2[24:40], min.cutoff = "q9")
FeaturePlot(PBMC.combined, features = M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Mouse2[41:57], min.cutoff = "q9")

FeaturePlot(PBMC.combined, features = c("Bgn","Areg","Cxcl2","Col1a2"), min.cutoff = "q9")


PBMC.combined <- AddModuleScore(object = PBMC.combined, features = as.data.frame(M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Mouse2),
                                ctrl = 5, name = 'M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')
FeaturePlot(object = PBMC.combined, features = 'M5930_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1')+
  scale_colour_gradient2(low = "#e8fff8", mid = "#00c480", high = "#19634a", 
                         guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

## Retrotransposition LINE-1
FeaturePlot(object = PBMC.combined, features = c("Lire1","Pol","Scn1a","Rap1a","Fut8","Arhgap15",
                                                 "Zfp352","Dennd1b","Gm9392","Mov10"))



