##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Version information ######
  # platform       x86_64-w64-mingw32          
  # arch           x86_64                      
  # os             mingw32                     
  # system         x86_64, mingw32             
  # status                                     
  # major          4                           
  # minor          1.1                         
  # year           2021                        
  # month          08                          
  # day            10                          
  # svn rev        80725                       
  # language       R                           
  # version.string R version 4.1.1 (2021-08-10)
  # nickname       Kick Things    

##### Current path and new folder setting #####
  Save.Path = paste0(getwd(),"/20220214_PBMC")
  # dir.create(Save.Path)
  SampleType = "PBMC"

##### Load library #####
  library(dplyr)
  library(stringr)
  library(scales)
  library(ggplot2)

##### Function setting  ##### 
  ## Call function  
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")

##### Abbreviated Notes #####
  # CCM: Cancer Cachexia Marker
  # SPA: Sex Pooled Analysis
  # SSA: Sex Separated Analysis
  # EO: Early Onset    
  # LO: Late Onset
  # CT: Cell Type
  
  # list:lt
  # dataframe: df

# ##### --------------------- Load RData --------------------- ##### 
#   load(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))
# 
# 
# ##### 09_0 GSEA Analysis (Geneset Prepare) #####
#   # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
#   # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
#   library(fgsea)
#   source("FUN_GSEA_LargeGeneSet.R")
#   source("FUN_HSsymbol2MMsymbol.R")
#   
#   # Geneset from GSEA
#   # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
#   Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"),
#                              col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"))),
#                              header = F,sep = "\t")
#   
#   # Convert Human gene to mouse
#   Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*1.5))
#   for (i in 1:nrow(Pathway.all)) {
#     #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
#     PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t() 
#     colnames(PathwayN)="Test"
#     PathwayN <- HSsymbol2MMsymbol(PathwayN,"Test")
#     Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
#     
#   }
#   
#   Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
#   colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))
#   
#   save.image(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))
# 
#   
##### --------------------- Load RData --------------------- ##### 
  load(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))
 
##### 09_1 GSEA Analysis (SPA) #####
  GSEA_Large <- list()
  GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large.df.TOP <- GSEA_Large.df
  
  
  pdf(file = paste0(Save.Path, "/PBMC_GSEA_CombineSex.pdf"),width = 15, height = 7 )
  
  for(i in 1:length(CellType.list)){
    
    gseaDat <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))
    
    
    #GSEA_Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all,10)
    GSEA_Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
    
    fgseaRes <- GSEA_Large.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA_Large.Output[["Pathway.all.list"]] 
    
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
    
    topPathways <- GSEA_Large.Output[["topPathways"]]
    
    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.",CellType.list[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large[[i]] <- Sum
    names(GSEA_Large)[[i]] <- paste0(CellType.list[i])
    
    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large.df <- rbind(GSEA_Large.df,fgseaRes2 )
    
    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large.df.TOP <- rbind(GSEA_Large.df.TOP,topPathways2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
    
  }
  
  
  ## GSEA_Large.Sum.TOP ##
  GSEA_Large.Sum.TOP <- rbind(GSEA_Large.df.TOP)
  GSEA_Large.Sum.TOP <- GSEA_Large.Sum.TOP[,!colnames(GSEA_Large.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large.Sum.TOP, file=paste0(Save.Path,"/Cachexia_PBMC_GSEA_Pathway_LargeTOP_CombineSex.txt"),sep="\t",
              row.names=F, quote = FALSE)
  
  
  ##### Bubble plot #####
    GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > 1,]
    GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < 0.05,]
    # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
    # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]
    library(ggplot2)
    library(scales)
    ggplot(GSEA_Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = padj, size = abs(NES))) + 
           geom_point() + scale_color_distiller()+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    ggplot(GSEA_Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
           geom_point() + scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d", 
                             guide = "colourbar",midpoint = 0)
    
    BBPlot <- ggplot(GSEA_Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
                     scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                     guide = "colourbar",midpoint = 0)+ geom_point() + 
                     theme(legend.position = "bottom") + theme_bw()+
                     theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    BBPlot <- BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                              XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
    BBPlot
    
    BBPlot2 <- BBPlot +theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
    BBPlot2
    print(BBPlot)
    print(BBPlot2)
  dev.off()
    
    
  ##### Clustering #####
    #GSEA.df <- GSEA_Large.df
    GSEA.df <- GSEA_Large.Sum.TOP.S
  
    
    # GSEA_Color.lt = list(red = "#ef476f",gray = "gray",blue = "#0077b6")
    GSEA_Color.lt = list(red = "#04873f",gray = "gray",blue = "#e3672d")
    
    GSEA.df$PhenoType <- factor(GSEA.df$PhenoType,
                                levels = sort(unique(as.character(GSEA.df$PhenoType))))
    
    BBPlot <- ggplot(GSEA.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ scale_x_discrete() 
    
    
    
    BBPlot
    
    BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                              XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
    
    
    df1 <- reshape2::dcast(GSEA.df,PhenoType~pathway,value.var = "NES")
    rownames(df1)<-df1$id
    
    df1.1<-df1[,2:ncol(df1)]
    df1.1[is.na(df1.1)] <- 0
    
    
    df1.1.clust.Pheno<-hclust(dist(df1.1))
    df1.1.clust.Pathway<-hclust(dist(t(df1.1)))
    
    
    library(ggtree)
    library(dplyr)
    PhenoType_Order <- data.frame(No=row.names(df1),PhenoType=df1[,1]) 
    PhenoType_Order$No <- as.numeric(PhenoType_Order$No)
    PhenoType_Treeclust_Order <- data.frame(No=df1.1.clust.Pheno[["order"]])
    PhenoType_Treeclust_Order <- left_join(PhenoType_Treeclust_Order,PhenoType_Order) 
    df1.1.clust.Pheno[["order"]] <- PhenoType_Treeclust_Order$PhenoType
    df1.1.clust.Pheno[["labels"]] <- PhenoType_Treeclust_Order$PhenoType
    
    p2 <- ggtree(df1.1.clust.Pheno)
    p2+
      geom_tiplab()+
      xlim(NA,7)
    p2.2<-p2+
      geom_tiplab()+
      xlim(NA,7)+
      #geom_tiplab(angle=90)+
      #theme_tree2()+
      layout_dendrogram()
    p2.2
    
    p3<-ggtree(df1.1.clust.Pathway)
    p3+
      geom_tiplab()+
      xlim(NA,7)
    ## Plot
    library(aplot)
    BBPlot_Cluster<- BBPlot %>%
      insert_left(p3,width = 0.2)
    
    BBPlot_Cluster
    
    BBPlotB <- BBPlot %>% 
      BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                     XtextSize=10,  YtextSize=7,AxisTitleSize=1, AspRat=2, XaThick=0.6, YaThick=0.6,
                     LegTextSize = 12,LegTitleSize=13)
    BBPlotB1 <- BBPlotB %>%
      insert_left(p3,width = 0.2)
    BBPlotB1
    
    # BBPlotB2 <-BBPlotB %>%
    #   insert_left(p3,width = 0.2)%>%
    #   insert_top(p2+layout_dendrogram(),height = 0.2)
    # BBPlotB2
    
    ## Export PDF 
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_PHH.pdf"),width = 17, height = 10 )
    BBPlotB1
    BBPlotB1
    #BBPlotB2
    dev.off()
    
    
  ##### Extract SubType #####
    
    ## T Cell
    # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
    GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]
    
    BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_size_area(max_size = 7)+
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    BBPlot_T

    BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                             XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
    # BBPlot_TB <- BBPlot_TB +theme(axis.title.y=element_blank(),
    #                  axis.text.y=element_blank(),
    #                  axis.ticks.y=element_blank())
    BBPlot_TB
    
    BBPlot_TB1 <- BBPlot_TB %>%
      insert_left(p3,width = 0.2)
    BBPlot_TB1
    
    
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SubType_T.pdf"),width = 17, height = 7 )
      BBPlot_TB
      BBPlot_TB1
    dev.off()
    
    
    ## Mac
    GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]
    
    BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_size_area(max_size = 5)+
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    BBPlot_Mac
    
    BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                 XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
    
    BBPlot_MacB1 <- BBPlot_MacB %>%
      insert_left(p3,width = 0.2)
    BBPlot_MacB1
    
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SubType_Mac.pdf"),width = 17, height = 20 )
      BBPlot_MacB
      BBPlot_MacB1
    dev.off()
    
    rm(GSEA_Large.Sum.TOP.S,p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
       df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)
    
  ##### save.image #####
    save.image(paste0(Save.Path,"/09_1_GSEA_Analysis_(SPA).RData"))    
  
##### 09_2 GSEA Analysis (SSA_MAle) #####
  GSEA_Large_Male <- list()
  GSEA_Large_Male.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large_Male.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large_Male.df.TOP <- GSEA_Large_Male.df
  
  
  pdf(file = paste0(Save.Path, "/PBMC_GSEA_SSA_Male.pdf"),width = 15, height = 7 )
  
  for(i in 1:length(CellType.list)){
    
    gseaDat <- CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))
    
    
    #GSEA_Large_Male.Output <- FUN_GSEA_Large_MaleGeneSet(ranks,Pathway.all,10)
    GSEA_Large_Male.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
    
    fgseaRes <- GSEA_Large_Male.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA_Large_Male.Output[["Pathway.all.list"]] 
    
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
    
    topPathways <- GSEA_Large_Male.Output[["topPathways"]]
    
    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.",CellType.list[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large_Male[[i]] <- Sum
    names(GSEA_Large_Male)[[i]] <- paste0(CellType.list[i])
    
    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large_Male.df <- rbind(GSEA_Large_Male.df,fgseaRes2 )
    
    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large_Male.df.TOP <- rbind(GSEA_Large_Male.df.TOP,topPathways2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
    
  }
  
  
  ## GSEA_Large_Male.Sum.TOP ##
  GSEA_Large_Male.Sum.TOP <- rbind(GSEA_Large_Male.df.TOP)
  GSEA_Large_Male.Sum.TOP <- GSEA_Large_Male.Sum.TOP[,!colnames(GSEA_Large_Male.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large_Male.Sum.TOP, file=paste0(Save.Path,"/Cachexia_PBMC_GSEA_Pathway_LargeTOP_SSA_Male.txt"),sep="\t",
              row.names=F, quote = FALSE)


  ##### Bubble plot #####
    GSEA_Large_Male.Sum.TOP.S <- GSEA_Large_Male.Sum.TOP[abs(GSEA_Large_Male.Sum.TOP$NES) > 1,]
    GSEA_Large_Male.Sum.TOP.S <- GSEA_Large_Male.Sum.TOP.S[abs(GSEA_Large_Male.Sum.TOP.S$padj) < 0.05,]
    # GSEA_Large_Male.Sum.TOP.S <- GSEA_Large_Male.Sum.TOP[abs(GSEA_Large_Male.Sum.TOP$padj) < 0.25,]
    # GSEA_Large_Male.Sum.TOP.S <- GSEA_Large_Male.Sum.TOP.S[abs(GSEA_Large_Male.Sum.TOP.S$pval) < 0.05,]
    library(ggplot2)
    library(scales)
    ggplot(GSEA_Large_Male.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = padj, size = abs(NES))) + 
      geom_point() +
      scale_color_distiller()
    
    ggplot(GSEA_Large_Male.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)
    
    BBPlot <- ggplot(GSEA_Large_Male.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
    BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                              XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
    
    
    print(BBPlot)
    dev.off()
    
    
  ##### Clustering #####
    #GSEA.df <- GSEA_Large.df
    GSEA.df <- GSEA_Large.Sum.TOP.S

    
    GSEA_Color.lt = list(red = "#ef476f",gray = "gray",blue = "#0077b6")
    
    GSEA.df$PhenoType <- factor(GSEA.df$PhenoType,
                                levels = sort(unique(as.character(GSEA.df$PhenoType))))
    
    BBPlot <- ggplot(GSEA.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ scale_x_discrete(limits = ) + theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    
    
    BBPlot
    
    BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                              XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
    
    
    df1 <- reshape2::dcast(GSEA.df,PhenoType~pathway,value.var = "NES")
    rownames(df1)<-df1$id
    
    df1.1<-df1[,2:ncol(df1)]
    df1.1[is.na(df1.1)] <- 0
    
    
    df1.1.clust.Pheno<-hclust(dist(df1.1))
    df1.1.clust.Pathway<-hclust(dist(t(df1.1)))
    
    
    library(ggtree)
    library(dplyr)
    PhenoType_Order <- data.frame(No=row.names(df1),PhenoType=df1[,1]) 
    PhenoType_Order$No <- as.numeric(PhenoType_Order$No)
    PhenoType_Treeclust_Order <- data.frame(No=df1.1.clust.Pheno[["order"]])
    PhenoType_Treeclust_Order <- left_join(PhenoType_Treeclust_Order,PhenoType_Order) 
    df1.1.clust.Pheno[["order"]] <- PhenoType_Treeclust_Order$PhenoType
    df1.1.clust.Pheno[["labels"]] <- PhenoType_Treeclust_Order$PhenoType
    
    p2 <- ggtree(df1.1.clust.Pheno)
    p2+
      geom_tiplab()+
      xlim(NA,7)
    p2.2<-p2+
      geom_tiplab()+
      xlim(NA,7)+
      #geom_tiplab(angle=90)+
      #theme_tree2()+
      layout_dendrogram()
    p2.2
    
    p3<-ggtree(df1.1.clust.Pathway)
    p3+
      geom_tiplab()+
      xlim(NA,7)
    ## Plot
    library(aplot)
    BBPlot_Cluster<- BBPlot %>%
      insert_left(p3,width = 0.2)
    
    BBPlot_Cluster
    
    BBPlotB <- BBPlot %>% 
      BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                     XtextSize=10,  YtextSize=7,AxisTitleSize=1, AspRat=2, XaThick=0.6, YaThick=0.6,
                     LegTextSize = 12,LegTitleSize=13)
    BBPlotB1 <-BBPlotB %>%
      insert_left(p3,width = 0.2)
    BBPlotB1
    
    BBPlotB2 <-BBPlotB %>%
      insert_left(p3,width = 0.2)%>%
      insert_top(p2+layout_dendrogram(),height = 0.2)
    BBPlotB2
    
    ## Export PDF 
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_PHH_Male.pdf"),width = 17, height = 10 )
    BBPlotB1
    BBPlotB1
    #BBPlotB2
    dev.off()
    
    
  ##### Extract SubType #####
    
    ## T Cell
    # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
    GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]
    
    BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_size_area(max_size = 7)+
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    BBPlot_T
    
    BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                             XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
    BBPlot_TB
    
    BBPlot_TB1 <- BBPlot_TB %>%
      insert_left(p3,width = 0.2)
    BBPlot_TB1
    
    
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SubType_T_Male.pdf"),width = 17, height = 7 )
      BBPlot_TB
      BBPlot_TB1
    dev.off()
    
    
    ## Mac
    GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]
    
    BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_size_area(max_size = 5)+
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    BBPlot_Mac
    
    BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                 XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
    BBPlot_MacB1 <- BBPlot_MacB %>%
      insert_left(p3,width = 0.2)
    BBPlot_MacB1
    
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SubType_Mac_Male.pdf"),width = 17, height = 20 )
      BBPlot_MacB
      BBPlot_MacB1
    dev.off()
    
    rm(GSEA_Large.Sum.TOP.S,p1,p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
       df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)
    save.image(paste0(Save.Path,"/09_2_GSEA_Analysis_(SSA_Male).RData"))      
    
##### 09_3 GSEA Analysis (SSA_Female) #####
  GSEA_Large_Female <- list()
  GSEA_Large_Female.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large_Female.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large_Female.df.TOP <- GSEA_Large_Female.df
  
  
  pdf(file = paste0(Save.Path, "/PBMC_GSEA_SSA_Female.pdf"),width = 15, height = 7 )
  
  for(i in 1:length(CellType.list)){
    
    gseaDat <- CCMarker_Female.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))
    
    
    #GSEA_Large_Female.Output <- FUN_GSEA_Large_FemaleGeneSet(ranks,Pathway.all,10)
    GSEA_Large_Female.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
    
    fgseaRes <- GSEA_Large_Female.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
    
    pathwaysH <- GSEA_Large_Female.Output[["Pathway.all.list"]] 
    
    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
    
    topPathways <- GSEA_Large_Female.Output[["topPathways"]]
    
    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway], 
                  ranks, 
                  fgseaRes, 
                  gseaParam = 0.5) + title( paste0("PBMC.",CellType.list[i]), adj = 0, line =3)
    
    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1
    
    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large_Female[[i]] <- Sum
    names(GSEA_Large_Female)[[i]] <- paste0(CellType.list[i])
    
    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large_Female.df <- rbind(GSEA_Large_Female.df,fgseaRes2 )
    
    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large_Female.df.TOP <- rbind(GSEA_Large_Female.df.TOP,topPathways2 )
    
    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
    
  }
  
  
  ## GSEA_Large_Female.Sum.TOP ##
  GSEA_Large_Female.Sum.TOP <- rbind(GSEA_Large_Female.df.TOP)
  GSEA_Large_Female.Sum.TOP <- GSEA_Large_Female.Sum.TOP[,!colnames(GSEA_Large_Female.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large_Female.Sum.TOP, file=paste0(Save.Path,"/Cachexia_PBMC_GSEA_Pathway_LargeTOP_SSA_Female.txt"),sep="\t",
              row.names=F, quote = FALSE)
  
  
  ##### Bubble plot #####
    GSEA_Large_Female.Sum.TOP.S <- GSEA_Large_Female.Sum.TOP[abs(GSEA_Large_Female.Sum.TOP$NES) > 1,]
    GSEA_Large_Female.Sum.TOP.S <- GSEA_Large_Female.Sum.TOP.S[abs(GSEA_Large_Female.Sum.TOP.S$padj) < 0.05,]
    # GSEA_Large_Female.Sum.TOP.S <- GSEA_Large_Female.Sum.TOP[abs(GSEA_Large_Female.Sum.TOP$padj) < 0.25,]
    # GSEA_Large_Female.Sum.TOP.S <- GSEA_Large_Female.Sum.TOP.S[abs(GSEA_Large_Female.Sum.TOP.S$pval) < 0.05,]
    library(ggplot2)
    library(scales)
    ggplot(GSEA_Large_Female.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = padj, size = abs(NES))) + 
      geom_point() +
      scale_color_distiller()
    
    ggplot(GSEA_Large_Female.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)
    
    BBPlot <- ggplot(GSEA_Large_Female.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                              XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
    
    
    print(BBPlot)
    dev.off()
    
    
  ##### Clustering #####
    #GSEA.df <- GSEA_Large.df
    GSEA.df <- GSEA_Large.Sum.TOP.S
    
    
    GSEA_Color.lt = list(red = "#ef476f",gray = "gray",blue = "#0077b6")
    
    GSEA.df$PhenoType <- factor(GSEA.df$PhenoType,
                                levels = sort(unique(as.character(GSEA.df$PhenoType))))
    
    BBPlot <- ggplot(GSEA.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ scale_x_discrete(limits = ) 
    
    
    
    BBPlot
    
    BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                              XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
    
    
    df1 <- reshape2::dcast(GSEA.df,PhenoType~pathway,value.var = "NES")
    rownames(df1)<-df1$id
    
    df1.1<-df1[,2:ncol(df1)]
    df1.1[is.na(df1.1)] <- 0
    
    
    df1.1.clust.Pheno<-hclust(dist(df1.1))
    df1.1.clust.Pathway<-hclust(dist(t(df1.1)))
    
    
    library(ggtree)
    library(dplyr)
    PhenoType_Order <- data.frame(No=row.names(df1),PhenoType=df1[,1]) 
    PhenoType_Order$No <- as.numeric(PhenoType_Order$No)
    PhenoType_Treeclust_Order <- data.frame(No=df1.1.clust.Pheno[["order"]])
    PhenoType_Treeclust_Order <- left_join(PhenoType_Treeclust_Order,PhenoType_Order) 
    df1.1.clust.Pheno[["order"]] <- PhenoType_Treeclust_Order$PhenoType
    df1.1.clust.Pheno[["labels"]] <- PhenoType_Treeclust_Order$PhenoType
    
    p2 <- ggtree(df1.1.clust.Pheno)
    p2+
      geom_tiplab()+
      xlim(NA,7)
    p2.2<-p2+
      geom_tiplab()+
      xlim(NA,7)+
      #geom_tiplab(angle=90)+
      #theme_tree2()+
      layout_dendrogram()
    p2.2
    
    p3<-ggtree(df1.1.clust.Pathway)
    p3+
      geom_tiplab()+
      xlim(NA,7)
    ## Plot
    library(aplot)
    BBPlot_Cluster<- BBPlot %>%
      insert_left(p3,width = 0.2)
    
    BBPlot_Cluster
    
    BBPlotB <- BBPlot %>% 
      BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                     XtextSize=10,  YtextSize=7,AxisTitleSize=1, AspRat=2, XaThick=0.6, YaThick=0.6,
                     LegTextSize = 12,LegTitleSize=13)
    BBPlotB1 <-BBPlotB %>%
      insert_left(p3,width = 0.2)
    BBPlotB1
    
    BBPlotB2 <-BBPlotB %>%
      insert_left(p3,width = 0.2)%>%
      insert_top(p2+layout_dendrogram(),height = 0.2)
    BBPlotB2
    
    ## Export PDF 
    pdf(file = paste0(PathName,"/",Version,"/PBMC_GSEA_Bubble_PHH_Female.pdf"),width = 17, height = 10 )
    BBPlotB1
    BBPlotB1
    #BBPlotB2
    dev.off()
    
    
  ##### Extract SubType #####
    
    ## T Cell
    # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
    GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]
    
    BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_size_area(max_size = 7)+
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    BBPlot_T
    
    BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                             XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
    BBPlot_TB
    
    BBPlot_TB1 <- BBPlot_TB %>%
      insert_left(p3,width = 0.2)
    BBPlot_TB1
    
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SubType_T_Female.pdf"),width = 17, height = 7 )
      BBPlot_TB
      BBPlot_TB1
    dev.off()
    
    
    ## Mac
    GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]
    
    BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_size_area(max_size = 5)+
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    BBPlot_Mac
    
    BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                 XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
    BBPlot_MacB1 <- BBPlot_MacB %>%
      insert_left(p3,width = 0.2)
    BBPlot_MacB1
    
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SubType_Mac_Female.pdf"),width = 17, height = 20 )
      BBPlot_MacB
      BBPlot_MacB1
    dev.off()
    
    rm(GSEA_Large.Sum.TOP.S,p1,p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
       df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)
    save.image(paste0(Save.Path,"/09_3_GSEA_Analysis_(SSA_Female).RData"))  
    
##### 09_4 GSEA Analysis (SSA) #####
  GSEA_Large.Sum.TOP <- rbind(GSEA_Large_Male.df.TOP,GSEA_Large_Female.df.TOP)
  GSEA_Large.Sum.TOP <- GSEA_Large.Sum.TOP[,!colnames(GSEA_Large.Sum.TOP) %in% c("leadingEdge")]
  write.table(GSEA_Large.Sum.TOP, file=paste0(Save.Path,"/PBMC_GSEA_Pathway_3Dataset_SSA.txt"),sep="\t",
              row.names=F, quote = FALSE)


  ##### Bubble plot #####
    GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > 1,]
    GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < 0.05,]
    # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
    # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]
    library(ggplot2)
    library(scales)
    
    pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SSA.pdf"),width = 17, height = 20 )
    
    ggplot(GSEA_Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = padj, size = abs(NES))) + 
      geom_point() +
      scale_color_distiller()
    
    ggplot(GSEA_Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f", 
                             guide = "colourbar",midpoint = 0)
    
    BBPlot <- ggplot(GSEA_Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
      geom_point() +
      scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d", 
                             guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
    BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                              XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)
    
    dev.off()
    
###########################################################################################
