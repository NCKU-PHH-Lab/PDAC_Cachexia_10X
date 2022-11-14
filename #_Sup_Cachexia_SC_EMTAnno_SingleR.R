## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)


##### Load Packages #####
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!require("ggpubr")) install.packages("ggpubr")

  library(ggpubr)
  library(tidyverse)
  library(Seurat)

  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("SingleR","scRNAseq","celldex","scran","scater","scuttle")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

#### Load data #####
  ## Load Ref
  load("D:/Dropbox/##_GitHub/#_scRNADataset/SeuratObject_CDS_PRJCA001063.RData")
  scRNA.SeuObj_Ref <- scRNA.SeuObj
  scRNA.SeuObj_Ref@meta.data$Type <- gsub("T","PT",scRNA.SeuObj_Ref@meta.data$Type)
  scRNA.SeuObj_Ref@meta.data$Type <- gsub("N","Ctrl",scRNA.SeuObj_Ref@meta.data$Type)

  # ## For small test
  # # CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref$CELL %in% sample(scRNA.SeuObj_Ref$CELL,1000)] ## For small test
  # CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref@meta.data[[1]] %in% sample(scRNA.SeuObj_Ref@meta.data[[1]],1000)] ## For small test
  # # ## For full data
  # # CTFeatures.SeuObj <- scRNA.SeuObj_Ref


  CTFeatures.SeuObj <- scRNA.SeuObj_Ref


  ## Load Tar

  # ## SeuObj_Tar
  # ## For small test
  # # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  # # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[[1]] %in% sample(scRNA.SeuObj@meta.data[[1]],1000)] ## For small test
  #
  # # load("D:/Dropbox/##_GitHub/#_scRNADataset/SeuratObject_GSE154778_PDAC_SC.RData")
  # load("D:/Dropbox/##_GitHub/##_Charlene/TumorHeterogeneity/2022-11-08_Com_PDAC/scRNA.SeuObj_CDS_PRJCA001063_Combine_ReBEC_Ori.RData")
  load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_SC_Main/06_Cell_type_annotation.RData")


  scRNA.SeuObj <- SC.combined

  # Clean up
  rm(list=setdiff(ls(), c("scRNA.SeuObj","CTFeatures.SeuObj","scRNA.SeuObj_Ref")))

  # # DefaultAssay(scRNA.SeuObj) <- "RNA"
  # # DefaultAssay(scRNA.SeuObj) <- "integrated"
  #
  # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
  # # Run the standard workflow for visualization and clustering
  # scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  # scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 160, verbose = FALSE)
  # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
  # scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:160)
  # scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)


##### Current path and new folder setting* #####
  ProjectName = paste0("EMTAnno_singleR_RefscRNA")
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }


##### Parameter setting* #####
  Remark = "PredbyscRNA" # c("PredbyCTDB","PredbyscRNA")
  RefType <- "BuiltIn_scRNA" # c("BuiltIn_celldex","BuiltIn_scRNA")
  celldexDatabase <- "HumanPrimaryCellAtlasData"
  # c("BlueprintEncodeData","DatabaseImmuneCellExpressionData","HumanPrimaryCellAtlasData","ImmGenData",
  #   "MonacoImmuneData","MouseRNAseqData","NovershternHematopoieticData")
  SingleR_DE_method <- "classic"


##### Set References #####
  if(RefType == "BuiltIn_celldex"){
    #### Database: Bulk reference setting for Cell type features ####
    library(celldex)

    if(celldexDatabase == "BlueprintEncodeData"){
      CTFeatures <- BlueprintEncodeData()
    }else if(celldexDatabase == "DatabaseImmuneCellExpressionData"){
      CTFeatures <- DatabaseImmuneCellExpressionData()
    }else if(celldexDatabase == "HumanPrimaryCellAtlasData"){
      CTFeatures <- HumanPrimaryCellAtlasData()
    }else if(celldexDatabase == "ImmGenData"){
      CTFeatures <- ImmGenData()
    }else if(celldexDatabase == "MonacoImmuneData"){
      CTFeatures <- MonacoImmuneData()
    }else if(celldexDatabase == "MouseRNAseqData"){
      CTFeatures <- MouseRNAseqData()
    }else if(celldexDatabase == "NovershternHematopoieticData"){
      CTFeatures <- NovershternHematopoieticData()
    }else{
      print("Error in database setting!")
    }


    #### Demo dataset ####
    # library(celldex)
    # hpca.se <- HumanPrimaryCellAtlasData()
    # hpca.se
  }else if(RefType =="BuiltIn_scRNA"){
    #### single-cell reference setting for Cell type features ####
    ## Prepossessing
    CTFeatures <- as.SingleCellExperiment(CTFeatures.SeuObj)
    CTFeatures$label <- CTFeatures@colData@listData[["Cell_type"]]
    CTFeatures <- CTFeatures[,!is.na(CTFeatures$label)]
    # CTFeatures <- logNormCounts(CTFeatures)
    #rm(CTFeatures.SeuObj)

    #### Demo dataset ####
    # library(scRNAseq)
    # sceM <- MuraroPancreasData()
    #
    # # One should normally do cell-based quality control at this point, but for
    # # brevity's sake, we will just remove the unlabelled libraries here.
    # sceM <- sceM[,!is.na(sceM$label)]
    #
    # # SingleR() expects reference datasets to be normalized and log-transformed.
    # library(scuttle)
    # sceM <- logNormCounts(sceM)

  }

  #########################################################################################################
  ##### Add scRNA-seq #####
  CTFeatures_Ori <- CTFeatures

  ##### Filter feature #####
  CTFeatures_Human <- HumanPrimaryCellAtlasData()
  CTFeatures_Human@colData@listData[["label.main"]] %>% unique()
  CTFeatures_Human_S <- CTFeatures_Human[,CTFeatures_Human@colData@listData[["label.main"]] %in%
                                           c("Endothelial_cells","Epithelial_cells","Fibroblasts",
                                             "Macrophage","T_cells","B_cell","Neuroepithelial_cell",
                                             "NK_cell","DC","Monocyte")]
  CTFeatures_Human_S2 <- CTFeatures_Human[,CTFeatures_Human@colData@listData[["label.main"]] %in%
                                           c("HSC_-G-CSF","HSC_CD34+","Myelocyte",
                                             "Pre-B_cell_CD34-","Pro-B_cell_CD34+","Pro-Myelocyte",
                                             "GMP")]

  CTFeatures_Human_sc <- as(CTFeatures_Human,"SingleCellExperiment")
  CTFeatures_Ori2 <- as(CTFeatures_Ori,"SummarizedExperiment")
  CTFeatures_Ori2 <- CTFeatures_Ori2[,1:ncol(CTFeatures_Human)]
  CTFeatures_Human_TTT <-  CTFeatures_Human[1:nrow(CTFeatures_Ori2),]
  # CTFeatures_Human_TTT <-  CTFeatures_Human[,1:ncol(CTFeatures_Ori2)]
  # colnames(CTFeatures_Human_TTT) <- colnames(CTFeatures_Ori2)
  CTFeatures <- cbind(CTFeatures_Human_TTT,CTFeatures_Ori2)


  #########################################################################################################



##### Set Target SeuObj #####
  ## Prepossessing
  scRNA <- as.SingleCellExperiment(scRNA.SeuObj)
  #### Demo dataset ####
  # library(scRNAseq)
  # hESCs <- LaMannoBrainData('human-es')
  # hESCs <- hESCs[,colSums(counts(hESCs)) > 0] # Remove libraries with no counts.
  # hESCs <- logNormCounts(hESCs)
  # hESCs <- hESCs[,1:100]

#### Run SingleR ####
  library(SingleR)
  if(RefType == "BuiltIn_celldex"){
    SingleR.lt <- SingleR(test = scRNA, ref = CTFeatures, assay.type.test=1,
                          labels = CTFeatures$label.main , de.method= SingleR_DE_method)#, de.method="wilcox") #  de.method = c("classic", "wilcox", "t")

  }else if(RefType =="BuiltIn_scRNA"){
    SingleR.lt <- SingleR(test = scRNA, ref = CTFeatures, assay.type.test=1,
                          labels = CTFeatures$label , de.method= SingleR_DE_method)#, de.method="wilcox") #  de.method = c("classic", "wilcox", "t")
  }

  SingleR.lt

  # Summarizing the distribution:
  CTCount_byCTDB.df <- table(SingleR.lt$labels) %>%
                       as.data.frame() %>%
                       dplyr::rename(Cell_Type = Var1, Count = Freq)

##### Annotation diagnostics #####
  p.ScoreHeatmap1 <- plotScoreHeatmap(SingleR.lt)
  p.ScoreHeatmap1
  p.DeltaDist1 <- plotDeltaDistribution(SingleR.lt, ncol = 3)
  p.DeltaDist1
  summary(is.na(SingleR.lt$pruned.labels))

  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_AnnoDiag.pdf"),
      width = 10,  height = 7
  )
    p.ScoreHeatmap1 %>% print()
    p.DeltaDist1 %>% print()
  dev.off()

  all.markers <- metadata(SingleR.lt)$de.genes
   scRNA@colData@listData[[paste0("labels_",SingleR_DE_method,"_",Remark)]] <- SingleR.lt$labels ## scRNA$labels <- SingleR.lt$labels

  # # Endothelial cell-related markers
  # library(scater)
  # plotHeatmap(scRNA, order_columns_by="labels",
  #             features = unique(unlist(all.markers[["Endothelial_cells"]])))



  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_HeatmapCTmarkers.pdf"),
      width = 12,  height = 7
  )
    for (i in 1:length(all.markers)) {
      plotHeatmap(scRNA, order_columns_by = paste0("labels_",SingleR_DE_method,"_",Remark),
                  features=unique(unlist(all.markers[[i]]))) %>% print()
    }
  dev.off()




  ## Plot UMAP
  scRNA.SeuObj@meta.data[[paste0("singleR_",SingleR_DE_method,"_",Remark)]]<- SingleR.lt$labels # scRNA.SeuObj$singleRPredbyCTDB <- SingleR.lt$labels
  p.CTPred1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by = paste0("singleR_",SingleR_DE_method,"_",Remark) ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CTPred1
  p.CT1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by ="Cell_type" ,label = TRUE, pt.size = 0.5) + NoLegend()
  p.CT1

  library(ggpubr)
  p.CTComp1 <- ggarrange(p.CT1, p.CTPred1, common.legend = TRUE, legend = "top")
  p.CTComp1

  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_CompareCTUMAP.pdf"),
      width = 12,  height = 7
  )
    p.CTComp1 %>% print()
  dev.off()



  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Remark,"_CompareCTUMAP.pdf"),
      width = 12,  height = 7
  )
    p.CTPred1 %>% BeautifyggPlot(XaThick=1,  YaThick=1, LegPos = c(1, 0.5)) %>%
                  print()
    ## Change color
    library(ggsci)
    p.CTPred1 + scale_color_npg() + coord_fixed(ratio = 1) + theme(legend.position = c(1, 0.5))

  dev.off()







##### Session information #####
  sessionInfo()
  ## Ref: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
  writeLines(capture.output(sessionInfo()), paste0(Save.Path,"/sessionInfo.txt"))



##### Save RData #####
  save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,".RData"))



##### ************************************************************************ #####
  scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA2 <- scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA  %>% as.character()
  TTT1 <- scRNA.SeuObj@meta.data
  TTT_Ref <- scRNA.SeuObj_Ref@meta.data
  TTT1[TTT1$CELL %in% TTT_Ref$CELL,]$singleR_classic_PredbyscRNA2 <- TTT_Ref$celltype %>% as.character()
  scRNA.SeuObj@meta.data <- TTT1

  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA2")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype")
##### ************************************************************************ #####


# ##### Verification (CellCheck) #####
#   #### Install ####
#   ## Check whether the installation of those packages is required
#   Package.set <- c("tidyverse","caret","cvms","DescTools","devtools","ggthemes")
#   for (i in 1:length(Package.set)) {
#     if (!requireNamespace(Package.set[i], quietly = TRUE)){
#       install.packages(Package.set[i])
#     }
#   }
#   ## Load Packages
#   # library(Seurat)
#   lapply(Package.set, library, character.only = TRUE)
#   rm(Package.set,i)
#
#   ## install CellCheck
#   # Install the CellCheck package
#   detach("package:CellCheck", unload = TRUE)
#   devtools::install_github("Charlene717/CellCheck")
#   # Load CellCheck
#   library(CellCheck)
#
#   #### Run CellCheck ####
#   ## Create check dataframe
#   CC.df <- scRNA.SeuObj@meta.data[,c("Cell_type","singleR_classic_PredbyscRNA", "singleR_classic_PredbyCTDB")]
#
#   CC.df <- data.frame(lapply(CC.df, as.character), stringsAsFactors=FALSE)
#
#   colnames(CC.df) <- c("Actual","Predict1","Predict2")
#   #CC.df$Actual <- as.character(CC.df$Actual)
#
#   CC.df$Predict2 <- gsub("_", " ", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("cells", "cell", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("Macrophage", "Macrophage cell", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("Fibroblasts", "Fibroblast cell", CC.df$Predict2)
#   CC.df$Predict2 <- gsub("Epithelial cell", "Ductal cell type 1", CC.df$Predict2)
#
#
#
#   CC.df[!CC.df$Predict2 %in% c(CC.df$Actual %>% unique()),]$Predict2 <- "Other"
#   # CC.df <- rbind(CC.df,"NotMatch")  #CC.df[nrow(CC.df)+1,1:ncol(CC.df)] <- "Other"
#
#   CC_Anno.df <- data.frame(TestID = c("Predict1","Predict2"),
#                            Tool = "singleR",
#                            Type = "PDAC",
#                            Set = c("singleRPredbyscRNA", "singleRPredbyCTDB"))
#
#   ## For one prediction
#   ## For one prediction
#   DisCMSet.lt = list(Mode = "One", Actual = "Actual", Predict = "Predict1" , FilterSet1 = "Tool", FilterSet2 = "singleR" , Remark = "") # Mode = c("One","Multiple")
#   BarChartSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
#   LinePlotSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
#   CCR_cm_DisMult.lt <- CellCheck_DisMult(CC.df, CC_Anno.df,
#                                          DisCMSet.lt = DisCMSet.lt,
#                                          BarChartSet.lt = BarChartSet.lt,
#                                          LinePlotSet.lt = LinePlotSet.lt,
#                                          Save.Path = Save.Path, ProjectName = paste0("CellCheck_",ProjectName))
#
#
#
# ##### Save RData #####
#   save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,".RData"))
#
#
#
#
