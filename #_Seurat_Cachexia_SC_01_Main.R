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
  # minor          1.2
  # year           2021
  # month          11
  # day            01
  # svn rev        81115
  # language       R
  # version.string R version 4.1.2 (2021-11-01)
  # nickname       Bird Hippie
  #--------------------------------#
  # R version 4.1.2
  # Seurat 4.0.2

##### Load Packages  #####
  library(Seurat)
  library(SeuratData)
  library(patchwork)
  library(ggplot2)
  library(ggpmisc)
  library(broom)

  library("stringr")
  library("magrittr")
  library("dplyr")

##### Function setting  #####
  ## Call function
  source("FUN_Cal_Mit.R")
  source("FUN_scRNAQC.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_Venn.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Beautify_UMAP.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")

##### Current path and new folder setting  #####
  Version = paste0(Sys.Date(),"_","SC_Main")
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path)

##### Load datasets  #####
  SC.data.TN138 <- Read10X(data.dir = paste0(getwd(),"/TN138/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN136 <- CreateSeuratObject(counts = SC.data.TN138, project = "EO.M", min.cells = 3, min.features = 200)
  SC.TN136
  SC.TN136@meta.data[["sample"]] <- rep(c("EO.M"), times=length(SC.TN136@meta.data[["orig.ident"]]))
  SC.TN136@meta.data[["ID"]] <- rep(c("TN138"), times=length(SC.TN136@meta.data[["orig.ident"]]))
  SC.TN136@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(SC.TN136@meta.data[["orig.ident"]]))  #EO: Early_Onset
  SC.TN136@meta.data[["Sex"]] <- rep(c("Male"), times=length(SC.TN136@meta.data[["orig.ident"]]))

  SC.data.TN139 <- Read10X(data.dir = paste0(getwd(),"/TN139/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN137 <- CreateSeuratObject(counts = SC.data.TN139, project = "LO.M", min.cells = 3, min.features = 200)
  SC.TN137
  SC.TN137@meta.data[["sample"]] <- rep(c("LO.M"), times=length(SC.TN137@meta.data[["orig.ident"]]))
  SC.TN137@meta.data[["ID"]] <- rep(c("TN139"), times=length(SC.TN137@meta.data[["orig.ident"]]))
  SC.TN137@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(SC.TN137@meta.data[["orig.ident"]]))  #LO: Late_Onset
  SC.TN137@meta.data[["Sex"]] <- rep(c("Male"), times=length(SC.TN137@meta.data[["orig.ident"]]))

  SC.data.TN146 <- Read10X(data.dir = paste0(getwd(),"/TN146/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN145 <- CreateSeuratObject(counts = SC.data.TN146, project = "LO.F", min.cells = 3, min.features = 200)
  SC.TN145
  SC.TN145@meta.data[["sample"]] <- rep(c("LO.F"), times=length(SC.TN145@meta.data[["orig.ident"]]))
  SC.TN145@meta.data[["ID"]] <- rep(c("TN146"), times=length(SC.TN145@meta.data[["orig.ident"]]))
  SC.TN145@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(SC.TN145@meta.data[["orig.ident"]]))
  SC.TN145@meta.data[["Sex"]] <- rep(c("Female"), times=length(SC.TN145@meta.data[["orig.ident"]]))

  SC.data.TN148 <- Read10X(data.dir = paste0(getwd(),"/TN148/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN147 <- CreateSeuratObject(counts = SC.data.TN148, project = "EO.F", min.cells = 3, min.features = 200)
  SC.TN147
  SC.TN147@meta.data[["sample"]] <- rep(c("EO.F"), times=length(SC.TN147@meta.data[["orig.ident"]]))
  SC.TN147@meta.data[["ID"]] <- rep(c("TN148"), times=length(SC.TN147@meta.data[["orig.ident"]]))
  SC.TN147@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(SC.TN147@meta.data[["orig.ident"]]))
  SC.TN147@meta.data[["Sex"]] <- rep(c("Female"), times=length(SC.TN147@meta.data[["orig.ident"]]))

##### Abbreviation Note #####
  # CCM: Cancer Cachexia Marker
  # SPA: Sex Pooled Analysis
  # SSA: Sex Separated Analysis
  # I:Intersection ; F: Female ; M: Male

  # EO: Early Onset
  # LO: Late Onset
  # CT: Cell Type

  # list:lt
  # dataframe: df

##### 01 Combine different datasets before QC  #####
  SC.list  <- c(SC.TN136,SC.TN137,SC.TN145,SC.TN147)
  rm(SC.TN136,SC.TN137,SC.TN145,SC.TN147, SC.data.TN138,SC.data.TN139,SC.data.TN146, SC.data.TN148)

  # normalize and identify variable features for each dataset independently
    set.seed(1) # Fix the seed
    SC.list <- lapply(X = SC.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

  # select features that are repeatedly variable across datasets for integration
    set.seed(1) # Fix the seed
    features <- SelectIntegrationFeatures(object.list = SC.list)

  ## Perform integration
    set.seed(1) # Fix the seed
    SC.anchors <- FindIntegrationAnchors(object.list = SC.list, anchor.features = features)
  # this command creates an 'integrated' data assay
    set.seed(1) # Fix the seed
    SC.combined <- IntegrateData(anchorset = SC.anchors)

    set.seed(1) # Fix the seed
    DefaultAssay(SC.combined) <- "integrated"

    #### Save RData ####
    save.image(paste0(Save.Path,"/01_Combine_different_datasets_before_QC.RData"))

##### 02 Quality Control  #####
  dir.create(paste0(Save.Path,"/SC_QC"))
  ## QC for all samples
  SC.combined_Ori <- SC.combined # Save the original obj
  #Test# SC.combined_Ori.list <- SplitObject(SC.combined_Ori, split.by = "ID")
  SC.combined_QCTry <- scRNAQC(SC.combined,FileName = paste0(Version,"/SC_QC/SC_QCTry"))

  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay

  rm(SC.anchors,SC.combined)

  ## QC for each sample for the new integration
  SC.TN136_QC <- scRNAQC(SC.list[[1]],FileName = paste0(Version,"/SC_QC/SC.TN136_QC"))
  SC.TN137_QC <- scRNAQC(SC.list[[2]],FileName = paste0(Version,"/SC_QC/SC.TN137_QC"))
  SC.TN145_QC <- scRNAQC(SC.list[[3]],FileName = paste0(Version,"/SC_QC/SC.TN145_QC"))
  SC.TN147_QC <- scRNAQC(SC.list[[4]],FileName = paste0(Version,"/SC_QC/SC.TN147_QC"))

  #### Save RData ####
  save.image(paste0(Save.Path,"/02_Quality_Control.RData"))


##### 03 Combine different data sets after QC  #####
  SC.list_QC  <- c(SC.TN136_QC ,SC.TN137_QC,SC.TN145_QC,SC.TN147_QC)
  rm(SC.TN136_QC ,SC.TN137_QC,SC.TN145_QC,SC.TN147_QC)

  # normalize and identify variable features for each dataset independently
    set.seed(1) # Fix the seed
    SC.list_QC <- lapply(X = SC.list_QC, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

  # select features that are repeatedly variable across datasets for integration
    set.seed(1) # Fix the seed
    features <- SelectIntegrationFeatures(object.list = SC.list_QC)

  ## Perform integration
    set.seed(1) # Fix the seed
    SC.anchors <- FindIntegrationAnchors(object.list = SC.list_QC, anchor.features = features)
    # this command creates an 'integrated' data assay
    set.seed(1) # Fix the seed
    SC.combined <- IntegrateData(anchorset = SC.anchors)


  ## Check QC
    scRNAQC(SC.combined,AddMitInf = "No",CheckOnly="Yes",FileName = paste0(Version,"/SC_QC/SC_QC_Check"))

  #### Save RData ####
  save.image(paste0(Save.Path,"/03_Combine_different_data_sets_after_QC.RData"))

##### 04 Perform an integrated analysis #####
  # Run the standard workflow for visualization and clustering

  # # # !!
  # set.seed(1) # Fix the seed
  # all.genes <- rownames(SC.combined)
  # SC.combined <- ScaleData(SC.combined, features = all.genes)

  ## Issues: re_clustering in seurat v3
  ## https://github.com/satijalab/seurat/issues/1528
  # DefaultAssay(SC.combined) <- "RNA"

  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(SC.combined) <- "integrated"

  ##?
  set.seed(1) # Fix the seed
  SC.combined <- ScaleData(SC.combined, verbose = FALSE)

  # ## Run if use filter
  # set.seed(1) # Fix the seed
  # SC.combined <- FindVariableFeatures(SC.combined)


  ### RunPCA
  # set.seed(1) # Fix the seed
  # SC.combined <- RunPCA(SC.combined, npcs = 30, verbose = FALSE)
  set.seed(1) # Fix the seed
  SC.combined <- RunPCA(SC.combined, features = VariableFeatures(object = SC.combined))

  print(SC.combined[["pca"]], dims = 1:5, nfeatures = 5)

  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/SC_PCA.pdf"),
    width = 10,  height = 8
  )
  VizDimLoadings(SC.combined, dims = 1:2, reduction = "pca")
  DimPlot(SC.combined, reduction = "pca")
  DimHeatmap(SC.combined, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(SC.combined, dims = 1:15, cells = 500, balanced = TRUE)
  DimHeatmap(SC.combined, dims = 16:30, cells = 500, balanced = TRUE)

  # # Determine the 'dimensionality' of the dataset
  # # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # # computation time
  # SC.combined <- JackStraw(SC.combined, num.replicate = 100)
  # SC.combined <- ScoreJackStraw(SC.combined, dims = 1:20)
  # JackStrawPlot(SC.combined, dims = 1:20)
  ElbowPlot(SC.combined, ndims = 50)
  dev.off()

  ElbowPlot(SC.combined, ndims = 50)

  ## Issues: RunUMAP causes R exit
  ## https://github.com/satijalab/seurat/issues/2259
  # The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
  # To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
  # This message will be shown once per session
  #### UMAP
  set.seed(1) # Fix the seed
  SC.combined <- RunUMAP(SC.combined, reduction = "pca", dims = 1:30)

  set.seed(1) # Fix the seed
  SC.combined <- FindNeighbors(SC.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  SC.combined <- FindClusters(SC.combined, resolution = 0.5)

  #### tSNE
  set.seed(1) # Fix the seed
  SC.combined <- RunTSNE(SC.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  SC.combined <- FindNeighbors(SC.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  SC.combined <- FindClusters(SC.combined, resolution = 0.5)


  ## Visualization
  DimPlot(SC.combined, reduction = "umap", group.by = "sample") %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)

        pdf(
          file = paste0(setwd(getwd()),"/",Version,"/SC_nlDR_Cluster.pdf"),
          width = 10,  height = 8
        )

        DimPlot(SC.combined, reduction = "umap", group.by = "sample") %>%
                BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)
        DimPlot(SC.combined, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14)

        DimPlot(SC.combined, reduction = "umap", ncol = 2,split.by = "sample", label = TRUE, label.size = 4) %>%
                BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                               SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9)

        DimPlot(SC.combined, reduction = "umap", ncol = 2,split.by = "Cachexia", label = TRUE, label.size = 4) %>%
                BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)
        DimPlot(SC.combined, reduction = "umap", ncol = 2,split.by = "Sex", label = TRUE, label.size = 4) %>%
          BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                         LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)

        ## tSNE
        DimPlot(SC.combined, reduction = "tsne", group.by = "sample") %>%
          BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)

        dev.off()
        # graphics.off()

  rm(SC.TN136_QC, SC.TN137_QC, SC.TN145_QC, SC.TN147_QC,SC.combined_QCTry)


  ##### Meta Table  #####

    ## Before QC
    Meta.df <- data.frame(matrix(nrow = 0,ncol = 3))
    colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
    Meta.df[1,1] <- c("EO.M")  # TN138
    Meta.df[1,2] <- ncol(SC.list[[1]]@assays[["RNA"]]@counts)
    Meta.df[1,3] <- nrow(SC.list[[1]]@assays[["RNA"]]@counts)

    Meta.df[2,1] <- c("LO.M")  # TN139
    Meta.df[2,2] <- ncol(SC.list[[2]]@assays[["RNA"]]@counts)
    Meta.df[2,3] <- nrow(SC.list[[2]]@assays[["RNA"]]@counts)

    Meta.df[3,1] <- c("LO.F")  # TN146
    Meta.df[3,2] <- ncol(SC.list[[3]]@assays[["RNA"]]@counts)
    Meta.df[3,3] <- nrow(SC.list[[3]]@assays[["RNA"]]@counts)

    Meta.df[4,1] <- c("EO.F")  # TN148
    Meta.df[4,2] <- ncol(SC.list[[4]]@assays[["RNA"]]@counts)
    Meta.df[4,3] <- nrow(SC.list[[4]]@assays[["RNA"]]@counts)

    # Summary to Meta table
    Meta.df[5,1] <- c("Summary")
    Meta.df[5,2] <- ncol(SC.combined_Ori@assays[["RNA"]]@counts)
    Meta.df[5,3] <- nrow(SC.combined_Ori@assays[["RNA"]]@counts)

    ## After QC
    colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
    Meta.df[6,1] <- c("EO.M.QC")  # TN138
    Meta.df[6,2] <- ncol(SC.list_QC[[1]]@assays[["RNA"]]@counts)
    Meta.df[6,3] <- nrow(SC.list_QC[[1]]@assays[["RNA"]]@counts)

    Meta.df[7,1] <- c("LO.M.QC")  # TN139
    Meta.df[7,2] <- ncol(SC.list_QC[[2]]@assays[["RNA"]]@counts)
    Meta.df[7,3] <- nrow(SC.list_QC[[2]]@assays[["RNA"]]@counts)

    Meta.df[8,1] <- c("LO.F.QC")  # TN146
    Meta.df[8,2] <- ncol(SC.list_QC[[3]]@assays[["RNA"]]@counts)
    Meta.df[8,3] <- nrow(SC.list_QC[[3]]@assays[["RNA"]]@counts)

    Meta.df[9,1] <- c("EO.F.QC")  # TN148
    Meta.df[9,2] <- ncol(SC.list_QC[[4]]@assays[["RNA"]]@counts)
    Meta.df[9,3] <- nrow(SC.list_QC[[4]]@assays[["RNA"]]@counts)

    # Summary to Meta table
    Meta.df[10,1] <- c("Summary")
    Meta.df[10,2] <- ncol(SC.combined@assays[["RNA"]]@counts)
    Meta.df[10,3] <- nrow(SC.combined@assays[["RNA"]]@counts)


    write.table( Meta.df ,
                 file = paste0(Save.Path,"/SC_CellCount_Meta.tsv"),
                 sep = "\t",
                 quote = F,
                 row.names = F
    )

  #### Save RData ####
  save.image(paste0(Save.Path,"/04_Perform_an_integrated_analysis.RData"))

##### 05 Identify conserved cell type markers  #####
  ## Identify conserved cell type markers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  set.seed(1) # Fix the seed
  SC.markers <- FindAllMarkers(SC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  library("magrittr")
  library("dplyr")

  # https://github.com/satijalab/seurat/issues/2960
  # Filter the top markers and plot the heatmap
  top_NSet = 7
  SC.markers %>%
    group_by(cluster) %>%
    top_n(n = top_NSet, wt = avg_log2FC) -> top_N
  SC.combined <- ScaleData(SC.combined, verbose = FALSE)
  DoHeatmap(SC.combined, features = top_N$gene) + NoLegend()
  write.table(top_N, file=paste0(Save.Path,"/SC_ClusterMarker_top",top_NSet,"Gene.txt"),sep="\t", row.names=T
              , quote = FALSE)
  write.table(SC.markers, file=paste0(Save.Path,"/SC_ClusterMarker_AllGene.txt"),sep="\t", row.names=T
              , quote = FALSE)

  pdf(
    file = paste0(Save.Path,"/SC_Heatmap_Cluster_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
  DoHeatmap(SC.combined, features = top_N$gene,size = 2,angle = 60) +
    scale_fill_gradient2(low="#5283ff",mid ="white", high ="#ff5c5c") +
    theme(axis.text.y = element_text(size  = 5)) +
    theme(legend.position = "bottom" )

  dev.off()


  # --------------- Check specific tissue marker --------------- #

  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/SC_nlDR_CTMarker.pdf"),
    width = 10,  height = 8
  )

  FeaturePlot(SC.combined, features = c("Krt19", "Prss1", "Chg8", "Cdh5", "Lum", "Rgs5", "Aif1",
                                        "Cd3d", "Ms4a1","Ms4a2"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  # Ductal cells
  FeaturePlot(SC.combined, features = c("Tff2", "Mmp7", "Tspan8", "Sox9", "Lcn2","Ceacam1","Ceacam5","Ceacam6","Krt19"),
              min.cutoff = "q9", ncol = 3, coord.fixed = 1)

  # Acinar cells
  # https://www.panglaodb.se/markers.html?cell_type=%27Acinar%20cells%27
  FeaturePlot(SC.combined, features = c("Cela3a","Prss1","Pnlip","Olfm4","Pnliprp1","Ctrc"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  # Pancreatic stellate cells
  # https://www.panglaodb.se/markers.html?cell_type=%27Pancreatic%20stellate%20cells%27
  FeaturePlot(SC.combined, features = c("Col6a1","Col6a2","Col6a3","Sfrp2","Thy1","Tnfaip6"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  # Fibroblast
  # https://www.novusbio.com/research-areas/cellular-markers/fibroblast-cell-markers.html
  FeaturePlot(SC.combined, features = c("Mas516","Fsp1","P4hb"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  # Fibroblast
  # https://cancerdiscovery.aacrjournals.org/content/9/8/1102
  FeaturePlot(SC.combined, features = c("Col1a1","Col3a1","Lum","Dcn"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  # Macrophages
  FeaturePlot(SC.combined, features = c("Mgl2","Clec4a2","Ccl12","Parp14","Fcgr1","Rbpj",
                                        "Naaa","March1",'Fgl2',"Cd68",'Tyrobp'), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  # Macrophages
  FeaturePlot(SC.combined, features = c("Chil3",'Cd86','Tnf',"Cxcl2"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  # M1
  FeaturePlot(SC.combined, features = c('Cd68','Cd80',"Tlr2","Tlr4"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  # M1
  FeaturePlot(SC.combined, features = c("Marco","Nos2","Tlr2","Cd80","Cd86","Csf2",
                                        "Tnf","Il1b","Il6","Tlr4","Cxcl2","Ifng","Il1r1"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)


  # M2
  FeaturePlot(SC.combined, features = c('IL-4','IL-13',"Tgm2","Mrc1","Marco"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  FeaturePlot(SC.combined, features = c("Cd4",'Cd8a','Chil3',"Il1b"), min.cutoff = "q9")
  FeaturePlot(SC.combined, features = c("Stat3",'Mtor','MyoD',"Erk","Nf","MuRF1","E3","Ifn","Sirt1","Lmcd1"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  FeaturePlot(SC.combined, features = c("Chil3","Csf1r","Mrc1","Pparg","Arg1","Cd163","Clec10a","Clec7a","Cd206",
                                        "Cd209","Ccl18","Fizz1"), min.cutoff = "q9", ncol = 3, coord.fixed = 1)


  # CAF
  # iCAF
  FeaturePlot(SC.combined, features = c('Clec3b','Col14a1',"Has1","IL6"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  FeaturePlot(SC.combined, features = c("Igf1","Ly6a","li6ra","Ly6c1","Il6","IL15ra"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  # myCAF
  FeaturePlot(SC.combined, features = c('Tagln','Thy1',"Col12a2","Thbs2"), min.cutoff = "q9")
  FeaturePlot(SC.combined, features = c("Cdkn2a", "Epha3", "Pdgfrb", "Myo1b", "Fap", "Heyl", "Inhba"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  # apCAF
  FeaturePlot(SC.combined, features = c('Cd74','Saa3',"Col12a2","Slpi"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)
  FeaturePlot(SC.combined, features = c("Nfe2l3","NKain4","Trem2","Cxadr","Cd74","Bcam","F11r","Ezr","Irf5","H2-Ab1"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  # PanCAF
  FeaturePlot(SC.combined, features = c('Col1a1','Col1a2',"Pdpn","Dcn"), min.cutoff = "q9",
              ncol = 3, coord.fixed = 1)

  dev.off()

  ## Summary
  markers.to.plot <- c("Sox9", "Csf1r","Lyz","Chil3","Col1a1","Col3a1","Dcn")

  FeaturePlot(SC.combined, features = markers.to.plot, min.cutoff = "q9", coord.fixed = 1 ,ncol =3)

  # Duc : "Sox9"; Macrophages: "Csf1r","Lyz","Chil3"; Fib: "Col1a1","Col3a1","Dcn"

  #### Save RData ####
  save.image(paste0(Save.Path,"/05_Identify_conserved_cell_type_markers.RData"))


##### 06 Cell type annotation  #####
  # SC.combined.copy <- SC.combined

  ## Duc: Ductal Cell;  Mac: Macrophages; Fib: Fibroblast
  ## Neu: Neutrophils; NK: NK Cell; Mast: Mast Cell; Ery: Erythrocytes;
  ## Thr: Thrombocytes
  SC.combined <- RenameIdents(SC.combined, `0` = "Duc1", `1` = "Duc2", `2` = "Duc3",
                              `3` = "Duc4", `4` = "Mac1", `5` = "Mac2", `6` = "Fib1", `7` = "Duc5",
                              `8` = "Mac3", `9` = "Duc6",
                              `10` = "Mac4", `11` = "Mac5", `12` = "Fib2",`13` = "Fib3",`14` = "Duc3")

  Cell_Type_Order.set <- c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" , "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
                           "Fib1", "Fib2", "Fib3")

  SC.combined$celltype <- Idents(SC.combined)
  # Idents(SC.combined) <- "celltype"

  ## Heatmap
  Heatmap_Color.lt <- list(low="#5283ff",mid ="#311669", high ="#ff5c5c")
  Heatmap_Color.lt <- list(low="#419e47",mid ="#311669", high ="#edff66")
  Heatmap_Color.lt <- list(low="#5283ff",mid ="white", high ="#ff5c5c")

  pdf(
    file = paste0(Save.Path,"/SC_Heatmap_CellType_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
  DoHeatmap(SC.combined, features = top_N$gene,size = 3,angle = 60) +
    scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                         mid = Heatmap_Color.lt[["mid"]],
                         high = Heatmap_Color.lt[["high"]]) +
    theme(axis.text.y = element_text(size  = 5)) +
    theme(legend.position = "bottom")+
    theme(aspect.ratio=1)
  dev.off()

  # ## Ch
  # # https://github.com/satijalab/seurat/issues/2960
  # set.seed(1) # Fix the seed
  # SC.combined <- ScaleData(SC.combined, verbose = FALSE)
  DoHeatmap(SC.combined, features = top_N$gene,size = 4,angle = 90) + NoLegend()
  DoHeatmap(SC.combined, features = top_N$gene,group.by = "celltype",size = 4,angle = 90) + NoLegend()

  # Color
  DoHeatmap(SC.combined, features = top_N$gene,group.by = "celltype",size = 3,angle = 60) +
    scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                         mid = Heatmap_Color.lt[["mid"]],
                         high = Heatmap_Color.lt[["high"]])


  DoHeatmap(SC.combined, features = top_N$gene,group.by = "seurat_clusters",size = 3,angle = 90) + NoLegend()

  pdf(
    file = paste0(Save.Path,"/SC_Heatmap_CellType_top",top_NSet,".pdf"),
    width = 10,  height = 8
  )
  DoHeatmap(SC.combined, features = top_N$gene,group.by = "celltype",size = 2,angle = 45) +
    scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                         mid = Heatmap_Color.lt[["mid"]],
                         high = Heatmap_Color.lt[["high"]])  +
    theme(axis.text.y = element_text(size = 5)) +
    theme(legend.position = "bottom" )

  dev.off()

  ## UMAP tSNE
  DimPlot(SC.combined, label = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5))

  DimPlot(SC.combined,group.by = "celltype",label.size = 7,label = TRUE,
          pt.size =2) %>% BeautifyUMAP(FileName = "_SC_nlDR_CellType")
  DimPlot(SC.combined,group.by = "sample",
          pt.size =0.5) %>% BeautifyUMAP(FileName = "_SC_nlDR_Sample")
  DimPlot(SC.combined,group.by = "seurat_clusters",label.size = 7, label = TRUE,
          pt.size =1) %>% BeautifyUMAP(FileName = "_SC_nlDR_Clusters")

  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/SC_nlDR_CellType_Sup.pdf"),
    width = 10,  height = 8
  )
  ##


  DimPlot(SC.combined, reduction = "umap", ncol = 2,split.by = "sample", label = TRUE) %>%
    BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                   SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9)

  DimPlot(SC.combined, reduction = "umap", ncol = 2,split.by = "Cachexia", label = TRUE, label.size = 4) %>%
    BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)

  DimPlot(SC.combined, reduction = "umap", ncol = 2,split.by = "Sex", label = TRUE, label.size = 4) %>%
    BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)

  ## tSNE
  DimPlot(SC.combined, reduction = "tsne", group.by = "sample") %>%
    BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)

  dev.off()


  ## DotPlot
  DotPlot_Color1.set <- c("#de3767", "#de3767", "#4169e1", "#4169e1")
  DotPlot_Color2.set <- c("#5b8e7d","#7b2cbf")
  DotPlot_Color3.set <- c("#de3767", "#4169e1")

  pdf(
    file = paste0(Save.Path,"/SC_DotPlot_CellType",".pdf"),
    width = 10,  height = 8
  )

  # https://satijalab.org/seurat/reference/dotplot
  DotPlot(SC.combined, features = markers.to.plot, cols = c("lightgrey", "blue"),
          dot.scale = 8) + RotatedAxis()%>%
    BeautifyggPlot(.,LegPos = "bottom",AxisTitleSize=1, TitleSize = 20, xangle =90,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1,XtextSize=12,  YtextSize=12)

  DotPlot(SC.combined, features = markers.to.plot, cols = DotPlot_Color1.set,
          dot.scale = 8, split.by = "sample") + RotatedAxis()


  # https://github.com/satijalab/seurat/issues/1541
  DotPlot(SC.combined, features = markers.to.plot, cols = DotPlot_Color2.set,
          dot.scale = 8, split.by = "Cachexia") + RotatedAxis()

  DotPlot(SC.combined, features = markers.to.plot, cols = DotPlot_Color3.set,
          dot.scale = 8, split.by = "Sex") + RotatedAxis()

  dev.off()

  rm(top_N, top_NSet)

  #### Save RData ####
  save.image(paste0(Save.Path,"/06_Cell_type_annotation.RData"))

  # ##### Export marker gene from specific cluster #####
  #   # For performing differential expression after integration, we switch back to the original data
  #   set.seed(1) # Fix the seed
  #   DefaultAssay(SC.combined) <- "RNA"
  #
  #   # nk.markers <- FindConservedMarkers(SC.combined, ident.1 = 6, grouping.var = "sample", verbose = FALSE)
  #   library(BiocManager)
  #   library(multtest)
  #   nk.markers <- FindConservedMarkers(SC.combined, ident.1 = 'NK', grouping.var = "sample", verbose = FALSE)
  #   head(nk.markers)
  #
  #   rm(nk.markers)

  # ##### Identify differential expressed genes across conditions  #####
  #   library(ggplot2)
  #   library(cowplot)
  #   theme_set(theme_cowplot())
  #   CD4T.cells <- subset(SC.combined, idents = "CD4+T")
  #   Idents(CD4T.cells) <- "Cachexia"
  #   avg.CD4T.cells <- as.data.frame(log1p(AverageExpression(CD4T.cells, verbose = FALSE)$RNA))
  #   avg.CD4T.cells$gene <- rownames(avg.CD4T.cells)
  #
  #   MacrophageM2 <- subset(SC.combined, idents = "Mac2")
  #   Idents(MacrophageM2) <- "Cachexia"
  #   avg.MacrophageM2 <- as.data.frame(log1p(AverageExpression(MacrophageM2, verbose = FALSE)$RNA))
  #   avg.MacrophageM2$gene <- rownames(avg.MacrophageM2)
  #
  #   genes.to.label = c("Sox17", "Mrpl15", "Lypla1", "Tcea1", "Rgs20", "Atp6v1h", "Rb1cc1", "4732440D04Rik", "St18")
  #   p1 <- ggplot(avg.CD4T.cells, aes(EO, LO)) + geom_point() + ggtitle("Cachexia T.cells")
  #   p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  #   p2 <- ggplot(avg.MacrophageM2, aes(EO, LO)) + geom_point() + ggtitle("Cachexia Macrophage")
  #   p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
  #   p1 + p2
  #   rm(p1 , p2 ,CD4T.cells, MacrophageM2, avg.CD4T.cells, avg.MacrophageM2)


##### 07 Count Cell number  #####
  Pheno.df <- data.frame(sample = SC.combined@meta.data[["sample"]],celltype = SC.combined@meta.data[["celltype"]],
                         Cachexia = SC.combined@meta.data[["Cachexia"]],Sex = SC.combined@meta.data[["Sex"]])

  Freq_sample.df <- table(Pheno.df$sample) %>% as.data.frame()
  Freq_CT.df <- table(Pheno.df$celltype) %>% as.data.frame()
  Freq_Cach.df <- table(Pheno.df$Cachexia) %>% as.data.frame()
  Freq_Sex.df <- table(Pheno.df$Sex) %>% as.data.frame()

  #
  Pheno_EO_M.df <- Pheno.df[Pheno.df$sample=="EO.M",]
  Pheno_LO_M.df <- Pheno.df[Pheno.df$sample=="LO.M",]
  Pheno_EO_F.df <- Pheno.df[Pheno.df$sample=="EO.F",]
  Pheno_LO_F.df <- Pheno.df[Pheno.df$sample=="LO.F",]

  # Count EO_M_CT
  Freq_EO_M_CT.df <- table(Pheno_EO_M.df$celltype) %>% as.data.frame()
  Freq_EO_M_CT.df <- data.frame(Type="EO.M",Freq_EO_M_CT.df)
  Freq_EO_M_CT.df$Percent <- Freq_EO_M_CT.df$Freq/sum(Freq_EO_M_CT.df$Freq)

  # Count LO_M_CT
  Freq_LO_M_CT.df <- table(Pheno_LO_M.df$celltype) %>% as.data.frame()
  Freq_LO_M_CT.df <- data.frame(Type="LO.M",Freq_LO_M_CT.df)
  Freq_LO_M_CT.df$Percent <- Freq_LO_M_CT.df$Freq/sum(Freq_LO_M_CT.df$Freq)

  # Count EO_F_CT
  Freq_EO_F_CT.df <- table(Pheno_EO_F.df$celltype) %>% as.data.frame()
  Freq_EO_F_CT.df <- data.frame(Type="EO.F",Freq_EO_F_CT.df)
  Freq_EO_F_CT.df$Percent <- Freq_EO_F_CT.df$Freq/sum(Freq_EO_F_CT.df$Freq)

  # Count LO_F_CT
  Freq_LO_F_CT.df <- table(Pheno_LO_F.df$celltype) %>% as.data.frame()
  Freq_LO_F_CT.df <- data.frame(Type="LO.F",Freq_LO_F_CT.df)
  Freq_LO_F_CT.df$Percent <- Freq_LO_F_CT.df$Freq/sum(Freq_LO_F_CT.df$Freq)

  # Combind all count of sample
  Freq_All.df <- rbind(Freq_EO_M_CT.df,Freq_LO_M_CT.df,
                       Freq_EO_F_CT.df,Freq_LO_F_CT.df)
  Freq_All.df <- data.frame(Index = row.names(Freq_All.df),Freq_All.df )
  colnames(Freq_All.df) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")
  # Freq_All.df$Index <- factor(Freq_All.df$Index,
  #                                              levels = Freq_All.df$Index)

  #### LinePlot ####
    # https://ithelp.ithome.com.tw/articles/10186047
    # Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
    #                                    levels = sort(unique(as.character(Freq_All.df$Cell_Type))))

    Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
                                    levels = Cell_Type_Order.set)

    CellNum_P1 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Number,
                                          colour = Pheno_Type, group = Pheno_Type)) +
      geom_line(linetype = "dashed",size=1.5) +
      geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

    CellNum_P1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                  XtextSize=15,  YtextSize=15, xangle = 90,
                                  LegTextSize = 15) +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P1
    CellNum_P1

    CellNum_P2 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Percent,
                                          colour = Pheno_Type, group = Pheno_Type)) +
      geom_line(linetype = "dashed",size=1.5) +
      geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

    CellNum_P2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=18, xangle = 90,
                                  LegTextSize = 15) +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P2
    CellNum_P2


  #### All type compare to Combine Sex ####
    ##
    Pheno_EO.df <- Pheno.df[Pheno.df$Cachexia=="EO",]
    Pheno_LO.df <- Pheno.df[Pheno.df$Cachexia=="LO",]

    # Count EO_CT
    Freq_EO_CT.df <- table(Pheno_EO.df$celltype) %>% as.data.frame()
    Freq_EO_CT.df <- data.frame(Type="EO",Freq_EO_CT.df)
    Freq_EO_CT.df$Percent <- Freq_EO_CT.df$Freq/sum(Freq_EO_CT.df$Freq)

    # Count LO_CT
    Freq_LO_CT.df <- table(Pheno_LO.df$celltype) %>% as.data.frame()
    Freq_LO_CT.df <- data.frame(Type="LO",Freq_LO_CT.df)
    Freq_LO_CT.df$Percent <- Freq_LO_CT.df$Freq/sum(Freq_LO_CT.df$Freq)

    # Combind all count of sample
    Freq_All_Ca.df <- rbind(Freq_EO_M_CT.df,Freq_LO_M_CT.df,
                            Freq_EO_F_CT.df,Freq_LO_F_CT.df,
                            Freq_EO_CT.df,Freq_LO_CT.df)

    Freq_All_Ca.df <- data.frame(Index = row.names(Freq_All_Ca.df),Freq_All_Ca.df )
    colnames(Freq_All_Ca.df) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")

    # Change the order
    # https://blog.csdn.net/weixin_48172266/article/details/117537465
    # CTOrder.set <- factor(Freq_All_Ca.df$Cell_Type,
    #                       levels = sort(unique(as.character(Freq_All_Ca.df$Cell_Type))))
    #
    # Freq_All_Ca.df <- Freq_All_Ca.df %>%
    #                   mutate(Cell_Type = CTOrder.set)

    Freq_All_Ca.df$Percent <- as.numeric(Freq_All_Ca.df$Percent)

    write.table( Freq_All_Ca.df ,
                 file = paste0(Save.Path,"/SC_CellCount_CT_Ca.tsv"),
                 sep = "\t",
                 quote = F,
                 row.names = F
    )

    # https://stackoverflow.com/questions/27350243/ggplot-line-graph-with-different-line-styles-and-markers/27350366
    # https://www.coder.work/article/6971741
    # https://stackoverflow.com/questions/11344561/controlling-line-color-and-line-type-in-ggplot-legend

  #### LinePlot ####
    Freq_All_Ca.df$Cell_Type <- factor(Freq_All_Ca.df$Cell_Type,
                                       levels = Cell_Type_Order.set)

    CellNum_P3 <- ggplot(Freq_All_Ca.df, aes(x = factor(Cell_Type), y = Number,
                                             colour = Pheno_Type,
                                             group = Pheno_Type,linetype=Pheno_Type
    )) +
      geom_line(size=1.5) +
      scale_linetype_manual(name="Pheno_Type",
                            values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                            labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")) +
      scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
      geom_point(shape = 12, size = 4, fill = "white") + theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


    CellNum_P3 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.4, 0.8),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=18,xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P3
    CellNum_P3

    library(eoffice)
    topptx(CellNum_P3,paste0(Save.Path,"/Temp.pptx"))

    CellNum_P4 <- ggplot(Freq_All_Ca.df, aes(x = factor(Cell_Type), y = Percent,
                                             colour = Pheno_Type,
                                             group = Pheno_Type,linetype=Pheno_Type
    )) +
      geom_line(size=1.5) +
      scale_linetype_manual(name="Pheno_Type",
                            values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                            labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")) +
      scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
      geom_point(shape = 12, size = 4, fill = "white") +
      theme(panel.border = element_blank(),panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
      #theme_set(theme_bw())+ # Remove the background
      theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

    CellNum_P4 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.15, 0.82),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P4
    CellNum_P4

  ##### Clean the dataframe #####
    CellPheno.lt <- list(Pheno.df, Pheno_EO.df, Pheno_LO.df,
                         Pheno_EO_M.df,Pheno_LO_M.df,Pheno_EO_F.df,Pheno_LO_F.df)
    names(CellPheno.lt) <- c("Pheno.df", "Pheno_EO.df", "Pheno_LO.df",
                             "Pheno_EO_M.df","Pheno_LO_M.df","Pheno_EO_F.df","Pheno_LO_F.df")
    rm(Pheno_EO.df, Pheno_LO.df, Pheno_EO_M.df,Pheno_LO_M.df,Pheno_EO_F.df,Pheno_LO_F.df)

    CellFreq.lt <- list(Freq_All_Ca.df, Freq_All.df, Freq_sample.df,
                        Freq_Sex.df,Freq_Cach.df, Freq_CT.df,
                        Freq_EO_CT.df, Freq_EO_F_CT.df, Freq_EO_M_CT.df,
                        Freq_LO_CT.df, Freq_LO_F_CT.df, Freq_LO_M_CT.df)
    names(CellFreq.lt) <- c("Freq_All_Ca.df", "Freq_All.df", "Freq_sample.df",
                            "Freq_Sex.df", "Freq_Cach.df", "Freq_CT.df",
                            "Freq_EO_CT.df", "Freq_EO_F_CT.df", "Freq_EO_M_CT.df",
                            "Freq_LO_CT.df", "Freq_LO_F_CT.df", "Freq_LO_M_CT.df")
    rm(Freq_sample.df,
       Freq_Sex.df,Freq_Cach.df, Freq_CT.df,
       Freq_EO_CT.df, Freq_EO_F_CT.df, Freq_EO_M_CT.df,
       Freq_LO_CT.df, Freq_LO_F_CT.df, Freq_LO_M_CT.df)

  ##### BarPlot #####
    # https://blog.gtwang.org/r/ggplot2-tutorial-layer-by-layer-plotting/3/
    colnames(Pheno.df) <- c("sample","Cell_Type","Cachexia","Sex")
    Pheno.df$Cell_Type <- factor(Pheno.df$Cell_Type,
                                 levels = sort(unique(as.character(Pheno.df$Cell_Type))))

    # sample
    BarPlot1_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=sample)) +
      geom_bar(position="dodge")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot1_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot1_1
    BarPlot1_1

    BarPlot1_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=sample)) +
      geom_bar(position="fill")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot1_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot1_2
    BarPlot1_2

    # Cachexia
    BarPlot2_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=Cachexia)) +
      geom_bar(position="dodge")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot2_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot2_1
    BarPlot2_1

    BarPlot2_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=Cachexia)) +
      geom_bar(position="fill")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot2_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot2_2
    BarPlot2_2

    # Sex
    BarPlot3_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=Sex)) +
      geom_bar(position="dodge")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot3_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot3_1
    BarPlot3_1
    BarPlot3_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=Sex)) +
      geom_bar(position="fill")+theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    BarPlot3_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                  LegTextSize = 15)  +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot3_2

    BarPlot3_2


  ##### Export PDF file #####
    pdf(file = paste0(Save.Path,"/SC_CellCount_LinePlot.pdf"),
        width = 7, height = 7 )
      CellNum_P4
      CellNum_P3
      CellNum_P1
      CellNum_P2
      BarPlot1_1
      BarPlot1_2
      BarPlot2_1
      BarPlot2_2
      BarPlot3_1
      BarPlot3_2
    dev.off() # graphics.off()

    rm(CellNum_P1, CellNum_P2, CellNum_P3, CellNum_P4, BarPlot1_1, BarPlot1_2,
       BarPlot2_1, BarPlot2_2, BarPlot3_1, BarPlot3_2)

    #### Save RData ####
    save.image(paste0(Save.Path,"/07_Count_Cell_number.RData"))


##### 08_1 Find CCmarker in different Cell type and VolcanoPlot (SSA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")
  SC.combined$celltype <- Idents(SC.combined)
  SC.combined$celltype.Cachexia <- paste(Idents(SC.combined), SC.combined$Cachexia, sep = "_")
  SC.combined$celltype.Cachexia.gender <- paste(Idents(SC.combined), SC.combined$Cachexia, SC.combined$Sex, sep = "_")
  Idents(SC.combined) <- "celltype.Cachexia.gender"

  SC.combined$Cachexia.gender <- paste(SC.combined$Cachexia, SC.combined$Sex, sep = "_")

  DefaultAssay(SC.combined) <- "RNA"

  ####-------------- Find Marker gene in Male --------------####
    CellType.list <- as.character(unique(SC.combined@meta.data[["celltype"]]))
    # CellType.list <- CellType.list[-9]

    dir.create(paste0(Save.Path,"/SC_SSA_Male_FindMarkers"))

    # About 15 mins
    CCMarker_Male.lt <- list()
    for(i in c(1:length(CellType.list))){
      try({
      CCMarker_Male.lt[[i]] <- Find_Markers(SC.combined,
                                        paste0(CellType.list[i],"_EO_Male"),
                                        paste0(CellType.list[i],"_LO_Male"),
                                        CellType.list[i],
                                        Path = Save.Path,
                                        ResultFolder = "SC_SSA_Male_FindMarkers")
      # names(CCMarker_Male.lt)[[i]] <- paste0("CCMarker_Male.lt.",CellType.list[i])
      names(CCMarker_Male.lt)[[i]] <- paste0(CellType.list[i])
      })
    }
    rm(i)

    CCMarker_Male.lt <- CCMarker_Male.lt[!unlist(lapply(CCMarker_Male.lt,is.null))]


    ## Generate pdf and tif file for Male VolcanoPlot
    dir.create(paste0(Save.Path,"/SC_SSA_Male_VolcanoPlot/"))

    pdf(file = paste0(Save.Path,"/SC_SSA_Male_VolcanoPlot/SC_SSA_Male_VolcanoPlot.pdf"),
        width = 7, height = 7 )
    for (i in 1:length(CellType.list)) {
      try({
      print(VolcanoPlot(CCMarker_Male.lt[[i]][["CCMarker.S"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+
            ggtitle(paste0("SC_Male_",CellType.list[i]))
            )
          })
    }
    dev.off() # graphics.off()
    rm(i)

    for (i in 1:length(CellType.list)) {
      try({
      tiff(file = paste0(Save.Path,"/SC_SSA_Male_VolcanoPlot/",CellType.list[i],".tif"),
           width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_Male.lt[[i]][["CCMarker.S"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Neg_List"]])+
              ggtitle(paste0("SC_Male_",CellType.list[i]))
            )

      graphics.off()
      })
    }
    rm(i)

  ####-------------- Find Marker gene in Female --------------####
    CellType.list <- as.character(unique(SC.combined@meta.data[["celltype"]]))
    # CellType.list <- CellType.list[-9] # Some cluster with cell lower than 3

    dir.create(paste0(Save.Path,"/SC_SSA_Female_FindMarkers"))

    # About 15 mins
    CCMarker_Female.lt <- list()
    for(i in c(1:length(CellType.list))){
      try({
        CCMarker_Female.lt[[i]] <- Find_Markers(SC.combined,
                                           paste0(CellType.list[i],"_EO_Female"),
                                           paste0(CellType.list[i],"_LO_Female"),
                                           CellType.list[i],
                                           Path = Save.Path,
                                           ResultFolder = "SC_SSA_Female_FindMarkers")
        # names(CCMarker_Female.lt)[[i]] <- paste0("CCMarker_Female.lt.",CellType.list[i])
        names(CCMarker_Female.lt)[[i]] <- paste0(CellType.list[i])
      })
    }
    rm(i)

    CCMarker_Female.lt <- CCMarker_Female.lt[!unlist(lapply(CCMarker_Female.lt,is.null))]

    ## Generate pdf and tif file for Female VolcanoPlot
    dir.create(paste0(Save.Path,"/SC_SSA_Female_VolcanoPlot/"))

    pdf(file = paste0(Save.Path,"/SC_SSA_Female_VolcanoPlot/SC_SSA_Female_VolcanoPlot.pdf"),
        width = 7, height = 7 )
    for (i in 1:length(CellType.list)) {
      try({
        print(VolcanoPlot(CCMarker_Female.lt[[i]][["CCMarker.S"]],
                          CCMarker_Female.lt[[i]][["CCMarker.S_Pos_List"]],
                          CCMarker_Female.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+
                ggtitle(paste0("SC_Female_",CellType.list[i]))
        )
      })
    }
    dev.off() # graphics.off()
    rm(i)

    for (i in 1:length(CellType.list)) {
      try({
        tiff(file = paste0(Save.Path,"/SC_SSA_Female_VolcanoPlot/",CellType.list[i],".tif"),
             width = 17, height = 17, units = "cm", res = 200)
        print(VolcanoPlot(CCMarker_Female.lt[[i]][["CCMarker.S"]],
                          CCMarker_Female.lt[[i]][["CCMarker.S_Pos_List"]],
                          CCMarker_Female.lt[[i]][["CCMarker.S_Neg_List"]])+
                ggtitle(paste0("SC_Female_",CellType.list[i]))
        )

        graphics.off()
      })
    }
    rm(i)

    #### Save RData ####
    save.image(paste0(Save.Path,"/08_1_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SSA).RData"))


##### 08_1 Find CCmarker in different Cell type and VennDiagrame (SSA_IntersectCT) ########
  ####-------------- Intersect_CellType --------------####
  #CCMarker_Male_Ori.lt <- CCMarker_Male.lt
  #CCMarker_Female_Ori.lt <- CCMarker_Female.lt
  #CellType_Ori.list <- CellType.list

  intersect_CellType <- intersect(names(CCMarker_Male.lt),names(CCMarker_Female.lt))

  CCMarker_Male.lt <- CCMarker_Male.lt[names(CCMarker_Male.lt) %in% intersect_CellType]
  CCMarker_Female.lt <- CCMarker_Female.lt[names(CCMarker_Female.lt) %in% intersect_CellType]

  CellType.list <- names(CCMarker_Male.lt)

  ####-------------- Venn Pos --------------####
  source("FUN_Venn.R")
  # pdf(file = paste0(Save.Path,"/SC_Female_VolcanoPlot.pdf"),width = 7, height = 7 )

  dir.create(paste0(Save.Path,"/SC_SSA_VennDiagrame"))
  Venn_CCMarker_Pos <- list()
  for(i in c(1:length(CellType.list))){
    try({
    Venn_CCMarker_Pos[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.S_Pos_List"]],
                                                          CCMarker_Female.lt[[paste0(CellType.list[i])]][["CCMarker.S_Pos_List"]],
                                                          CellType.list[i],"Pos","#9d0208","#f08080",SampleType="SC",
                                                          PathName = paste0(Save.Path,"/SC_SSA_VennDiagrame"))
    names(Venn_CCMarker_Pos)[[i]] <- paste0("Venn_CCMarker.",CellType.list[i],"_Pos")
    })
  }
  rm(i)

  ####-------------- Venn Neg --------------####
  Venn_CCMarker_Neg <- list()
  for(i in c(1:length(CellType.list))){
    try({
    Venn_CCMarker_Neg[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.S_Neg_List"]],
                                                          CCMarker_Female.lt[[paste0(CellType.list[i])]][["CCMarker.S_Neg_List"]],
                                                          CellType.list[i],"Neg","#00296b","#1368aa",SampleType="SC",
                                                          PathName = paste0(Save.Path,"/SC_SSA_VennDiagrame"))

    names(Venn_CCMarker_Neg)[[i]] <- paste0("Venn_CCMarker.",CellType.list[i],"_Neg")
    })
  }
  rm(i)

    #### Save RData ####
    save.image(paste0(Save.Path,"/08_1_Find_CCmarker_in_different_Cell_type_and_VennDiagrame(SSA_IntersectCT).RData"))


##### 08_2 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")

  Idents(SC.combined) <- "celltype.Cachexia"
  #CellType.list <- as.character(unique(SC.combined@meta.data[["celltype"]]))
  dir.create(paste0(Save.Path,"/SC_SPA_FindMarkers"))

    CCMarker_SPA.lt <- list()
    for(i in c(1:length(CellType.list))){
      try({
        CCMarker_SPA.lt[[i]] <- Find_Markers(SC.combined,
                                          paste0(CellType.list[i],"_EO"),
                                          paste0(CellType.list[i],"_LO"),
                                          CellType.list[i],
                                          Path = Save.Path,
                                          ResultFolder = "SC_SPA_FindMarkers")

        # names(CCMarker_SPA.lt)[[i]] <- paste0("CCMarker_SPA.lt.",CellType.list[i])
        names(CCMarker_SPA.lt)[[i]] <- paste0(CellType.list[i])
      })
    }
    rm(i)

    CCMarker_SPA.lt <- CCMarker_SPA.lt[!unlist(lapply(CCMarker_SPA.lt,is.null))]


    ## Generate pdf and tif file for VolcanoPlot
    dir.create(paste0(Save.Path,"/SC_SPA_VolcanoPlot/"))

    pdf(file = paste0(Save.Path,"/SC_SPA_VolcanoPlot/SC_SPA_VolcanoPlot.pdf"),width = 7, height = 7 )
      for (i in 1:length(CellType.list)) {
        try({
          print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                            CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                            CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+
                            ggtitle(paste0("SC_",CellType.list[i]))
          )
        })
      }
      # graphics.off()
    dev.off()
    rm(i)

    for (i in 1:length(CellType.list)) {
      try({
        tiff(file = paste0(Save.Path,"/SC_SPA_VolcanoPlot/SC_SPA_VolcanoPlot",CellType.list[i],".tif"), width = 17, height = 17, units = "cm", res = 200)
        print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                          CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                          CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]])+ ggtitle(paste0("SC_",CellType.list[i]))
        )

        graphics.off()
      })
    }
    rm(i)

    #### Save RData ####
    save.image(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))


#####------------------------------------------------------------------------------------------------------------#####

##### 09_0 GSEA Analysis (Geneset Prepare) #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(fgsea)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_GSEA_ggplot.R")

  # Geneset from GSEA
  # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
  Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"),
                             col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"))),
                             header = F,sep = "\t")

  # Convert Human gene to mouse
  Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*1.5))
  for (i in 1:nrow(Pathway.all)) {
    #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
    PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()
    colnames(PathwayN)="Test"
    PathwayN <- HSsymbol2MMsymbol(PathwayN,"Test")
    Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)

  }

  Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
  colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))

    #### Save RData ####
    save.image(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))

  ##### 09_1 GSEA Analysis (SPA) #####
    GSEA_Large <- list()
    GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
    colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
    GSEA_Large.df.TOP <- GSEA_Large.df

    dir.create(paste0(Save.Path,"/SC_GSEA"))


    pdf(file = paste0(Save.Path, "/SC_GSEA/SC_GSEA_SPA.pdf"),width = 15, height = 7 )

      for(i in 1:length(CellType.list)){

        gseaDat <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
        gseaDat <- data.frame(row.names(gseaDat),gseaDat)
        colnames(gseaDat)[[1]] <- c("Gene")
        ranks <- gseaDat$avg_log2FC
        names(ranks) <- gseaDat$Gene
        # head(ranks)
        # barplot(sort(ranks, decreasing = T))

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
                      gseaParam = 0.5) + title( paste0("SC.",CellType.list[i]), adj = 0, line =3)

        plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[1,1])))
        #plotEnrichment_Pos1
        plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
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
        GSEA_Large.df.TOP <- rbind(GSEA_Large.df.TOP, topPathways2)

        rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)

      }

    dev.off()

    ## GSEA_Large.Sum.TOP ##
    GSEA_Large.Sum.TOP <- rbind(GSEA_Large.df.TOP)
    GSEA_Large.Sum.TOP <- GSEA_Large.Sum.TOP[,!colnames(GSEA_Large.Sum.TOP) %in% c("leadingEdge")]
    write.table(GSEA_Large.Sum.TOP, file=paste0(Save.Path,"/SC_GSEA/SC_GSEA_Pathway_LargeTOP_SPA.txt"),sep="\t",
                row.names=F, quote = FALSE)

    ##### Bubble plot #####
      library(ggplot2)
      library(scales)
      GSEA_Color.lt = list(high = "#ef476f",mid = "white",low = "#0077b6")

      GSEA_Large.Sum.TOP$PhenoType <- factor(GSEA_Large.Sum.TOP$PhenoType,
                                        levels = Cell_Type_Order.set)

      GSEA_ggplot_SPA.lt <- GSEA_ggplot(GSEA_Large.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
      GSEA_Large.Sum.TOP.S <- GSEA_ggplot_SPA.lt[["GSEA_TOP.df"]]
      # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > 1,]
      # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < 0.05,]
      # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
      # GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]

        pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SPA.pdf"),width = 17, height = 12 )
          GSEA_ggplot_SPA.lt[["BBPlot_Ori"]]
          GSEA_ggplot_SPA.lt[["BBPlot"]]
          GSEA_ggplot_SPA.lt[["BBPlot2"]]
          GSEA_ggplot_SPA.lt[["BBPlotB1"]]
          GSEA_ggplot_SPA.lt[["BBPlotB1"]]
        dev.off()


    ##### Extract SubType #####
        ## T Cell
        # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
        GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]

        BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
          geom_point() +
          scale_size_area(max_size = 7)+
          scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                                 guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

        BBPlot_T

        BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                                 XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
        # BBPlot_TB <- BBPlot_TB +theme(axis.title.y=element_blank(),
        #                  axis.text.y=element_blank(),
        #                  axis.ticks.y=element_blank())
        BBPlot_TB

        BBPlot_TB1 <- BBPlot_TB %>%
          insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
        BBPlot_TB1


        pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SPA_SubType_T.pdf"),width = 17, height = 7 )
          BBPlot_TB
          BBPlot_TB1
        dev.off()


        ## Mac
        GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]

        BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
          geom_point() +
          scale_size_area(max_size = 5)+
          scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                                 guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

        BBPlot_Mac

        BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                                     XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)

        BBPlot_MacB1 <- BBPlot_MacB %>%
          insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
        BBPlot_MacB1

        pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SPA_SubType_Mac.pdf"),width = 17, height = 20 )
          BBPlot_MacB
          BBPlot_MacB1
        dev.off()

        rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
           df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

    ##### save.image #####
    save.image(paste0(Save.Path,"/09_1_GSEA_Analysis_(SPA).RData"))

  ##### 09_2 GSEA Analysis (SSA_MAle) #####
    GSEA_Large_Male <- list()
    GSEA_Large_Male.df <- as.data.frame(matrix(nrow=0,ncol=10))
    colnames(GSEA_Large_Male.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
    GSEA_Large_Male.df.TOP <- GSEA_Large_Male.df


    pdf(file = paste0(Save.Path, "/SC_GSEA/SC_GSEA_SSA_Male.pdf"),width = 15, height = 7 )

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
                      gseaParam = 0.5) + title( paste0("SC.",CellType.list[i]), adj = 0, line =3)

        plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[1,1])))
        #plotEnrichment_Pos1
        plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
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

    dev.off()

    ## GSEA_Large_Male.Sum.TOP ##
    GSEA_Large_Male.Sum.TOP <- rbind(GSEA_Large_Male.df.TOP)
    GSEA_Large_Male.Sum.TOP <- GSEA_Large_Male.Sum.TOP[,!colnames(GSEA_Large_Male.Sum.TOP) %in% c("leadingEdge")]
    write.table(GSEA_Large_Male.Sum.TOP, file=paste0(Save.Path,"/SC_GSEA/SC_GSEA_Pathway_LargeTOP_SSA_Male.txt"),sep="\t",
                row.names=F, quote = FALSE)



    ##### Bubble plot #####
      library(ggplot2)
      library(scales)

      GSEA_Large_Male.Sum.TOP$PhenoType <- factor(GSEA_Large_Male.Sum.TOP$PhenoType,
                                             levels = Cell_Type_Order.set)

      GSEA_ggplot_SSA_Male.lt <- GSEA_ggplot(GSEA_Large_Male.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
      GSEA_Large_Male.Sum.TOP.S <- GSEA_ggplot_SSA_Male.lt[["GSEA_TOP.df"]]

      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SSA_Male.pdf"),width = 17, height = 12 )
        GSEA_ggplot_SSA_Male.lt[["BBPlot_Ori"]]
        GSEA_ggplot_SSA_Male.lt[["BBPlot"]]
        GSEA_ggplot_SSA_Male.lt[["BBPlot2"]]
        GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
        GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
      dev.off()


    ##### Extract SubType #####

      ## T Cell
      # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
      GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]

      BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_T

      BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                               XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_TB

      BBPlot_TB1 <- BBPlot_TB %>%
        insert_left(GSEA_ggplot_SSA_Male.lt[["Y_Order"]],width = 0.2)
      BBPlot_TB1


      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SSA_Male_SubType_T.pdf"),width = 17, height = 7 )
      BBPlot_TB
      BBPlot_TB1
      dev.off()


      ## Mac
      GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]

      BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 5)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_Mac

      BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                   XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_MacB1 <- BBPlot_MacB %>%
        insert_left(GSEA_ggplot_SSA_Male.lt[["Y_Order"]],width = 0.2)
      BBPlot_MacB1

      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SSA_Male_Mac.pdf"),width = 17, height = 20 )
      BBPlot_MacB
      BBPlot_MacB1
      dev.off()

      rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
         df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

    ##### save.image #####
    save.image(paste0(Save.Path,"/09_2_GSEA_Analysis_(SSA_Male).RData"))

  ##### 09_3 GSEA Analysis (SSA_Female) #####
    GSEA_Large_Female <- list()
    GSEA_Large_Female.df <- as.data.frame(matrix(nrow=0,ncol=10))
    colnames(GSEA_Large_Female.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
    GSEA_Large_Female.df.TOP <- GSEA_Large_Female.df


    pdf(file = paste0(Save.Path, "/SC_GSEA/SC_GSEA_SSA_Female.pdf"),width = 15, height = 7 )

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
                      gseaParam = 0.5) + title( paste0("SC.",CellType.list[i]), adj = 0, line =3)

        plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[1,1])))
        #plotEnrichment_Pos1
        plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("SC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
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

    dev.off()

    ## GSEA_Large_Female.Sum.TOP ##
    GSEA_Large_Female.Sum.TOP <- rbind(GSEA_Large_Female.df.TOP)
    GSEA_Large_Female.Sum.TOP <- GSEA_Large_Female.Sum.TOP[,!colnames(GSEA_Large_Female.Sum.TOP) %in% c("leadingEdge")]
    write.table(GSEA_Large_Female.Sum.TOP, file=paste0(Save.Path,"/SC_GSEA/SC_GSEA_Pathway_LargeTOP_SSA_Female.txt"),sep="\t",
                row.names=F, quote = FALSE)

    ##### Bubble plot #####
      library(ggplot2)
      library(scales)

      GSEA_Large_Female.Sum.TOP$PhenoType <- factor(GSEA_Large_Female.Sum.TOP$PhenoType,
                                                  levels = Cell_Type_Order.set)

      GSEA_ggplot_SSA_Female.lt <- GSEA_ggplot(GSEA_Large_Female.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
      GSEA_Large_Female.Sum.TOP.S <- GSEA_ggplot_SSA_Female.lt[["GSEA_TOP.df"]]

      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SSA_Female.pdf"),width = 17, height = 12 )
        GSEA_ggplot_SSA_Female.lt[["BBPlot_Ori"]]
        GSEA_ggplot_SSA_Female.lt[["BBPlot"]]
        GSEA_ggplot_SSA_Female.lt[["BBPlot2"]]
        GSEA_ggplot_SSA_Female.lt[["BBPlotB1"]]
        GSEA_ggplot_SSA_Female.lt[["BBPlotB1"]]
      dev.off()


    ##### Extract SubType #####

      ## T Cell
      # GSEA_T.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
      GSEA_T.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]

      BBPlot_T <- ggplot(GSEA_T.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_T

      BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                               XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_TB

      BBPlot_TB1 <- BBPlot_TB %>%
        insert_left(GSEA_ggplot_SSA_Female.lt[["Y_Order"]],width = 0.2)
      BBPlot_TB1

      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SubType_T_Female.pdf"),width = 17, height = 7 )
      BBPlot_TB
      BBPlot_TB1
      dev.off()


      ## Mac
      GSEA_Mac.df <- GSEA_Large.Sum.TOP.S[grep("Mac",GSEA_Large.Sum.TOP.S$PhenoType),]

      BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 5)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_Mac

      BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                   XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_MacB1 <- BBPlot_MacB %>%
        insert_left(GSEA_ggplot_SSA_Female.lt[["Y_Order"]],width = 0.2)
      BBPlot_MacB1

      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SubType_Mac_Female.pdf"),width = 17, height = 20 )
      BBPlot_MacB
      BBPlot_MacB1
      dev.off()

      rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
         df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

    ##### save.image #####
    save.image(paste0(Save.Path,"/09_3_GSEA_Analysis_(SSA_Female).RData"))

  ##### 09_4 GSEA Analysis (SSA) #####
    GSEA_Large_Male.df.TOP2 <- GSEA_Large_Male.df.TOP
    GSEA_Large_Female.df.TOP2 <- GSEA_Large_Female.df.TOP
    GSEA_Large.df.TOP2 <- GSEA_Large.df.TOP

    GSEA_Large_Male.df.TOP2$PhenoType <- paste0("M_", GSEA_Large_Male.df.TOP2$PhenoType)
    GSEA_Large_Female.df.TOP2$PhenoType <- paste0("F_", GSEA_Large_Female.df.TOP2$PhenoType)
    GSEA_Large.df.TOP2$PhenoType <- paste0("SPA_", GSEA_Large.df.TOP2$PhenoType)


    GSEA_Large_SumTOP_Sex.df <- rbind(GSEA_Large.df.TOP2,GSEA_Large_Male.df.TOP2,GSEA_Large_Female.df.TOP2)
    GSEA_Large_SumTOP_Sex.df <- GSEA_Large_SumTOP_Sex.df[,!colnames(GSEA_Large_SumTOP_Sex.df) %in% c("leadingEdge")]

    GSEA_ggplot_SSA.lt <- GSEA_ggplot(GSEA_Large_SumTOP_Sex.df,NES_Th = 1.5, padj_Th = 0.01)
    GSEA_Large_SumTOP_Sex.df.S <- GSEA_ggplot_SSA.lt[["GSEA_TOP.df"]]

    write.table(GSEA_Large_SumTOP_Sex.df, file=paste0(Save.Path,"/SC_GSEA/SC_GSEA_Pathway_3Dataset_SSA.txt"),sep="\t",
                row.names=F, quote = FALSE)

    ##### Extract SubType (SC) #####
      ## SPA
      GSEA_Large_SumTOP_SPA.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("SPA",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

      BBPlot_SPA <- ggplot(GSEA_Large_SumTOP_SPA.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_SPA

      BBPlot_SPAB <- BBPlot_SPA %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                               XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_SPAB

      BBPlot_SPAB1 <- BBPlot_SPAB %>%
        insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
      BBPlot_SPAB1

      ## Female
      GSEA_Large_SumTOP_F.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("F_",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

      BBPlot_F <- ggplot(GSEA_Large_SumTOP_F.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_F

      BBPlot_F_B <- BBPlot_F %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                   XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_F_B

      BBPlot_F_B1 <- BBPlot_F_B %>%
        insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
      BBPlot_F_B1


      ## Male
      GSEA_Large_SumTOP_M.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("M_",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

      BBPlot_M <- ggplot(GSEA_Large_SumTOP_M.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_M

      BBPlot_M_B <- BBPlot_M %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_M_B

      BBPlot_M_B1 <- BBPlot_M_B %>%
        insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
      BBPlot_M_B1

      ## T Cell
      GSEA_Large_SumTOP_T.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("T$",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

      BBPlot_T <- ggplot(GSEA_Large_SumTOP_T.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_T

      BBPlot_T_B <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_T_B

      BBPlot_T_B1 <- BBPlot_T_B %>%
        insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
      BBPlot_T_B1

      ## Macrophage
      GSEA_Large_SumTOP_Mac.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("Mac",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

      BBPlot_Mac <- ggplot(GSEA_Large_SumTOP_Mac.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_Mac

      BBPlot_Mac_B <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_Mac_B

      BBPlot_Mac_B1 <- BBPlot_Mac_B %>%
        insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
      BBPlot_Mac_B1

      ## Mast
      GSEA_Large_SumTOP_Mast.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("Mast",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

      BBPlot_Mast <- ggplot(GSEA_Large_SumTOP_Mast.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
        geom_point() +
        scale_size_area(max_size = 7)+
        scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                               guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

      BBPlot_Mast

      BBPlot_Mast_B <- BBPlot_Mast %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                                    XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
      BBPlot_Mast_B

      BBPlot_Mast_B1 <- BBPlot_Mast_B %>%
        insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
      BBPlot_Mast_B1

      ##
      BBPlotB <- GSEA_ggplot_SSA.lt[["BBPlot"]] %>%
        BeautifyggPlot(LegPos = c(-0.5, -0.05),LegBox = "horizontal",LegDir="horizontal", xangle =90,
                       XtextSize=10,  YtextSize=7,AxisTitleSize=1, AspRat=0.8, XaThick=0.6, YaThick=0.6,
                       LegTextSize = 12,LegTitleSize=13)
      BBPlotB
      #BBPlotB <- GSEA_ggplot_SSA.lt[["BBPlot"]] + theme(aspect.ratio=1)

      BBPlotB1 <- BBPlotB %>%
        insert_left(GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.5)
      BBPlotB1


    ##### Export Bubble plot #####

      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_Sum.pdf"),width = 35, height = 17 )
        GSEA_ggplot_SSA.lt[["BBPlotB1"]]
        BBPlotB1
        BBPlot_SPAB
        BBPlot_SPAB1
        BBPlot_F_B
        BBPlot_F_B1
        BBPlot_M_B
        BBPlot_M_B1
        BBPlot_T_B
        BBPlot_T_B1
        BBPlot_Mac_B
        BBPlot_Mac_B1
        BBPlot_Mast_B
        BBPlot_Mast_B1
      dev.off()



    ##### Bubble plot #####
      GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df[abs(GSEA_Large_SumTOP_Sex.df$NES) > 1,]
      GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df.S[abs(GSEA_Large_SumTOP_Sex.df.S$padj) < 0.05,]
      # GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df[abs(GSEA_Large_SumTOP_Sex.df$padj) < 0.25,]
      # GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df.S[abs(GSEA_Large_SumTOP_Sex.df.S$pval) < 0.05,]
      library(ggplot2)
      library(scales)

      pdf(file = paste0(Save.Path,"/SC_GSEA/SC_GSEA_Bubble_SSA.pdf"),width = 17, height = 20 )


        ggplot(GSEA_Large_SumTOP_Sex.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
          geom_point() + scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f",
                                 guide = "colourbar",midpoint = 0)
        ggplot(GSEA_Large_SumTOP_Sex.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
          geom_point() + scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d",
                                                guide = "colourbar",midpoint = 0)


        BBPlot <- ggplot(GSEA_Large_SumTOP_Sex.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
          geom_point() + scale_colour_gradient2(low = "#ef476f", mid = "white", high = "#0077b6",
                                 guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")

        BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                  XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=1, XaThick=0.8, YaThick=0.8 ,OL_Thick = 1.5)

      dev.off()

    ##### save.image #####
    save.image(paste0(Save.Path,"/09_4_GSEA_Analysis_(SSA).RData"))


###########################################################################################

##### Interested genes #####

  FeaturePlot(SC.combined, features = c("Top2a", "Ptk2"), min.cutoff = "q9")
  Test <- intersect(row.names(SC.combined@assays[["RNA"]]@counts) ,  c("Top2a", "Ptk2"))
  Test2 <- intersect(row.names(SC.combined@assays[["integrated"]]@data) ,  c("Top2a", "Ptk2"))
  FeaturePlot(SC.combined, features = c("Kras", "Braf", "Myc","Egfr","Pik3ca"), min.cutoff = "q9",ncol = 3)

# -------------------------------------- #

###########################################################################################



# ## Export Seurat Object in 10X format
# # https://github.com/satijalab/seurat/issues/884
# library(DropletUtils)
# write10xCounts(SC.list_QC_Try[["TN138"]]@assays[["RNA"]]@counts, path = paste0(setwd(getwd()),"/Test"))
