##### Reference #####
# ggtree for Bubble plot # https://zhuanlan.zhihu.com/p/352143854
# Let's Plot 7: Clustered Dot Plots in the ggverse # https://davemcg.github.io/post/lets-plot-scrna-dotplots/

# ggtree annotation # http://www.randigriffin.com/2017/05/11/primate-phylogeny-ggtree.html 
# ggtree annotation # https://4va.github.io/biodatasci/r-ggtree.html
# clustered dotplot for single-cell RNAseq # https://divingintogeneticsandgenomics.rbind.io/post/clustered-dotplot-for-single-cell-rnaseq/

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Version information ######
  # _                           
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

##### Load libray #####
  library(dplyr)
  library(stringr)
  library(ggplot2)  
  library("magrittr")
  source("FUN_Beautify_ggplot.R")

##### Path and folder setting ##### 
  Save.Path = paste0(getwd(),"/20220119_GSEA_GO_PBMC")
  dir.create(Save.Path)
  SampleType = "PBMC"
  
  RawData.Path = ""
  
##### Load data ##### 
  GSEA.df <- read.csv(paste0(Save.Path,"/Cachexia_PBMC_GSEA_Pathway_LargeTOP_CombineSex.txt"),
                      sep = "\t")
  GSEA_No.df <- read.delim2(paste0(Save.Path,"/GO_Term.txt"),
                  col.names = 1:max(count.fields(paste0(Save.Path,"/GO_Term.txt"))),
                  header = F, sep = "\t")

##### Deal with data #####
  GSEA_No_In.df <- as.data.frame(GSEA_No.df[,1])
  colnames(GSEA_No_In.df) <- "Pathway"
  
  ## Test 
  #GSEA_No_In_GO.df <- str_match(GSEA_No_In.df$Pathway,"\\(GO:([0-9]{7})\\)")
  GSEA_No_In_GO.df <- data.frame(GO=str_extract(GSEA_No_In.df$Pathway,"\\(GO:([0-9]{7})\\)"))
  GSEA_No_In_Path.df <- data.frame(pathway=str_replace(GSEA_No_In.df$Pathway," \\(GO:([0-9]{7})\\)",""))
  GSEA_No_In.df <- cbind(GSEA_No_In.df,GSEA_No_In_Path.df,GSEA_No_In_GO.df)
  #NG# GSEA_No_In.df$GO <- str_replace(GSEA_No_In.df$GO,"[\\(\\)]","")
  GSEA_No_In.df$GO <- str_replace(GSEA_No_In.df$GO,"\\(","")
  GSEA_No_In.df$GO <- str_replace(GSEA_No_In.df$GO,"\\)","")
  GSEA_No_In.df$pathway <- tolower(GSEA_No_In.df$pathway)
  GSEA_No_In.df$pathway <- str_replace_all(GSEA_No_In.df$pathway,"[:punct:]"," ")
  GSEA_No_In.df$pathway <- str_replace_all(GSEA_No_In.df$pathway,"  "," ")
  # ## Try 
  # GSEATTT <- data.frame(Pathway=GSEA_No_In.df[grep("\\(GO:", GSEA_No_In.df$Pathway),])
  # GSEATTT2 <- data.frame(Pathway=grep("\\(GO:", GSEA_No_In.df$Pathway, value = TRUE))
  # 
  # GSEATTT3 <- regmatches(
  #     GSEA_No_In.df$Pathway, 
  #     regexec("\\(GO:([0-9]{7})\\)", GSEA_No_In.df$Pathway)
  #   )
  
  GSEA_M.df <- GSEA.df
  GSEA_M.df$pathway <- GSEA_M.df$pathway %>% str_replace(.,"GOBP_","") %>%
                       str_replace_all(.,"_"," ") %>%
                       tolower(.)
  
  
  GSEA_Fin.df <- left_join(GSEA_M.df,GSEA_No_In.df,by="pathway")
  # Test
  sum(is.na(GSEA_Fin.df$GO))
  
  # Remove the NA (Remove from new version of GO Term)
  # (20220125 Remove record: (1)covalent chromatin modification; (2)regulation of gene silencing)
  GSEA_Fin.df <- GSEA_Fin.df[!is.na(GSEA_Fin.df$GO),]
  
##### Export df for REVIGO #####  
  write.table( GSEA_Fin.df ,
             file = paste0(Save.Path,"/",SampleType, "_GOBP_forREVIGO.tsv"),
             sep = "\t",
             quote = F,
             row.names = F
  )
  
  GSEA_Fin_S.df <- GSEA_Fin.df %>% group_by(GO) %>% slice(which.min(padj))
  
  write.table( GSEA_Fin_S.df ,
               file = paste0(Save.Path,"/",SampleType, "_GOBP_forREVIGO_S.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )

##### Import REVIGO #####
# http://revigo.irb.hr/  
  RevigoTreeMap.df <- read.delim2(paste0(Save.Path,"/PBMC_GOBP_RevigoTreeMap.csv"),
              col.names = 1:max(count.fields(paste0(Save.Path,"/PBMC_GOBP_RevigoTreeMap.csv"))),
              header = F, sep = ",")
  
  
  RevigoTreeMap.df <- RevigoTreeMap.df[-1:-4,1:7]
  
  colnames(RevigoTreeMap.df) <- RevigoTreeMap.df[1,]
  colnames(RevigoTreeMap.df)[1] <- "GO"
  colnames(RevigoTreeMap.df)<- RevigoTreeMap.df %>% colnames() %>% gsub(" ","",.)
  
  GSEA_Fin.df <- left_join(GSEA_Fin.df,RevigoTreeMap.df,by="GO")
  length(unique(GSEA_Fin.df$Representative))
  
  
  GSEA_Fin.df2 <- GSEA_Fin.df
  GSEA_Fin.df2$Representative <- str_replace(GSEA_Fin.df2$Representative," ","")
  GSEA_Fin.df2$Name <- str_replace(GSEA_Fin.df2$Name," ","")
  
  # for (i in 1:nrow(GSEA_Fin.df2)) {
  #   if( GSEA_Fin.df2$Representative[!is.na(GSEA_Fin.df2$Representative[i]),] == "null"){
  #     GSEA_Fin.df2$Representative[i] <- GSEA_Fin.df2$pathway[i]
  #   }
  # }
  
  ## Replace null
  null.set <- GSEA_Fin.df2$Representative=="null" & !is.na(GSEA_Fin.df2$Representative) 
  GSEA_Fin.df2$Representative[null.set] <- GSEA_Fin.df2$Name[null.set]
  rm(null.set)
  
  length(unique(GSEA_Fin.df2$Representative))
  unique(GSEA_Fin.df2$Representative)
  
  ## Replace NA
  GSEA_Fin.df2[is.na(GSEA_Fin.df2$Representativ),]$Representative <- "other"
  length(unique(GSEA_Fin.df2$Representative))
  unique(GSEA_Fin.df2$Representative)
    
  GSEA_Fin.df <- GSEA_Fin.df2
  
##### Original Bubble Plot #####
  BBPlot <- ggplot(GSEA_Fin.df, aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
  
  BBPlot
  
  BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                            XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)

##
  BBPlot2 <- ggplot(GSEA_Fin.df,aes(x=GO , y = PhenoType, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
  
  BBPlot2
  
  BBPlot2 %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                            XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)



##### Clustering #####
  df1 <- reshape2::dcast(GSEA_Fin.df,PhenoType~GO,value.var = "NES")
  rownames(df1) <- df1$id
  
  df1.1 <- df1[,2:ncol(df1)]
  df1.1[is.na(df1.1)] <- 0
  
  
  df1.1.clust.Pheno <-hclust(dist(df1.1))
  df1.1.clust.GO <-hclust(dist(t(df1.1)))


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
  
  p3<-ggtree(df1.1.clust.GO)
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
                     XtextSize=10,  YtextSize=7,AxisTitleSize=1, AspRat=2, XaThick=0.6, YaThick=0.6)
  
  # Ok
  BBPlotB %>%
    insert_left(p3,width = 0.2)
  
  
  
  BBPlotB %>%
     insert_left(p3,width = 0.2)%>%
     insert_top(p2+layout_dendrogram(),height = 0.2)
  
  ### Try
  ha = HeatmapAnnotation(foo = 1:12)
  ha = HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE))
  # Heatmap(top_annotation = ha, right_annotation = row_ha)
  BBPlotB %>%
    insert_top(ha)
  
  
  BBPlotB + geom_bar(mapping = aes(x = PhenoType, y = Representative, fill = Representative))
  
  pdf(file = paste0(PathName,"/",RVersion,"/PBMC_GSEA_Bubble.pdf"),width = 17, height = 12 )
  BBPlotB %>%
    insert_left(p3,width = 0.2)
  
  BBPlotB %>%
    insert_left(p3,width = 0.2)%>%
    insert_top(p2+layout_dendrogram(),height = 0.2)
  dev.off()
  

##### Extract SubType #####

  ## T Cell
  # GSEA_T.df <- GSEA.df[grep("T",GSEA.df$PhenoType),]
  GSEA_T.df <- GSEA_Fin.df[GSEA_Fin.df$PhenoType %in% c("CD4+T","CD8+T","T"),]
  
  BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = GO, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_size_area(max_size = 7)+
    scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
  
  BBPlot_T
  
  BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                            XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  
  pdf(file = paste0(PathName,"/",RVersion,"/PBMC_GSEA_Bubble_SubType_T.pdf"),width = 17, height = 7 )
  BBPlot_TB %>%
    insert_left(p3,width = 0.2)
  dev.off()
  
  
  ## Mac
  GSEA_Mac.df <- GSEA.df[grep("Mac",GSEA.df$PhenoType),]
  
  BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
    geom_point() +
    scale_size_area(max_size = 5)+
    scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
  
  BBPlot_Mac
  
  BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                               XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  
  pdf(file = paste0(PathName,"/",RVersion,"/PBMC_GSEA_Bubble_SubType_Mac.pdf"),width = 17, height = 20 )
  BBPlot_MacB %>%
    insert_left(p3,width = 0.2)
  dev.off()

