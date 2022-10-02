##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  #####
  library(Seurat)
  library(SeuratData)
  library(patchwork)
  library(ggplot2)
  library(ggpmisc)
  library(broom)
  library(aplot)

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
  source("FUN_GSEA_Run_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")

##### Current path and new folder setting*  #####
  Version = paste0(Sys.Date(),"_","PBMC_Main")
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path)

##### Load datasets*  #####
  PBMC.data.TN138 <- Read10X(data.dir = paste0(getwd(),"/TN138/monocle/outs/filtered_gene_bc_matrices/mm10"))
  PBMC.TN138 <- CreateSeuratObject(counts = PBMC.data.TN138, project = "EO.M", min.cells = 3, min.features = 200)
  PBMC.TN138
  PBMC.TN138@meta.data[["sample"]] <- rep(c("EO.M"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))
  PBMC.TN138@meta.data[["ID"]] <- rep(c("TN138"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))
  PBMC.TN138@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))  #EO: Early_Onset
  PBMC.TN138@meta.data[["Sex"]] <- rep(c("Male"), times=length(PBMC.TN138@meta.data[["orig.ident"]]))

  PBMC.data.TN139 <- Read10X(data.dir = paste0(getwd(),"/TN139/monocle/outs/filtered_gene_bc_matrices/mm10"))
  PBMC.TN139 <- CreateSeuratObject(counts = PBMC.data.TN139, project = "LO.M", min.cells = 3, min.features = 200)
  PBMC.TN139
  PBMC.TN139@meta.data[["sample"]] <- rep(c("LO.M"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))
  PBMC.TN139@meta.data[["ID"]] <- rep(c("TN139"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))
  PBMC.TN139@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))  #LO: Late_Onset
  PBMC.TN139@meta.data[["Sex"]] <- rep(c("Male"), times=length(PBMC.TN139@meta.data[["orig.ident"]]))

  PBMC.data.TN146 <- Read10X(data.dir = paste0(getwd(),"/TN146/monocle/outs/filtered_gene_bc_matrices/mm10"))
  PBMC.TN146 <- CreateSeuratObject(counts = PBMC.data.TN146, project = "LO.F", min.cells = 3, min.features = 200)
  PBMC.TN146
  PBMC.TN146@meta.data[["sample"]] <- rep(c("LO.F"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))
  PBMC.TN146@meta.data[["ID"]] <- rep(c("TN146"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))
  PBMC.TN146@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))
  PBMC.TN146@meta.data[["Sex"]] <- rep(c("Female"), times=length(PBMC.TN146@meta.data[["orig.ident"]]))

  PBMC.data.TN148 <- Read10X(data.dir = paste0(getwd(),"/TN148/monocle/outs/filtered_gene_bc_matrices/mm10"))
  PBMC.TN148 <- CreateSeuratObject(counts = PBMC.data.TN148, project = "EO.F", min.cells = 3, min.features = 200)
  PBMC.TN148
  PBMC.TN148@meta.data[["sample"]] <- rep(c("EO.F"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))
  PBMC.TN148@meta.data[["ID"]] <- rep(c("TN148"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))
  PBMC.TN148@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))
  PBMC.TN148@meta.data[["Sex"]] <- rep(c("Female"), times=length(PBMC.TN148@meta.data[["orig.ident"]]))

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
  PBMC.list  <- c(PBMC.TN138,PBMC.TN139,PBMC.TN146,PBMC.TN148)
  rm(PBMC.TN138,PBMC.TN139,PBMC.TN146,PBMC.TN148, PBMC.data.TN138,PBMC.data.TN139,PBMC.data.TN146, PBMC.data.TN148)

  # normalize and identify variable features for each dataset independently
    set.seed(1) # Fix the seed
    PBMC.list <- lapply(X = PBMC.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

  # select features that are repeatedly variable across datasets for integration
    set.seed(1) # Fix the seed
    features <- SelectIntegrationFeatures(object.list = PBMC.list)

  ## Perform integration
    set.seed(1) # Fix the seed
    PBMC.anchors <- FindIntegrationAnchors(object.list = PBMC.list, anchor.features = features)
  # this command creates an 'integrated' data assay
    set.seed(1) # Fix the seed
    PBMC.combined <- IntegrateData(anchorset = PBMC.anchors)

    set.seed(1) # Fix the seed
    DefaultAssay(PBMC.combined) <- "integrated"

    #### Save RData ####
    save.image(paste0(Save.Path,"/01_Combine_different_datasets_before_QC.RData"))

##### 02 Quality Control  #####
  dir.create(paste0(Save.Path,"/PBMC_QC"))
  ## QC for all samples
  PBMC.combined_Ori <- PBMC.combined # Save the original obj
  #Test# PBMC.combined_Ori.list <- SplitObject(PBMC.combined_Ori, split.by = "ID")
  PBMC.combined_QCTry <- scRNAQC(PBMC.combined,FileName = paste0(Version,"/PBMC_QC/PBMC_QCTry"))

  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay

  rm(PBMC.anchors,PBMC.combined)

  ## QC for each sample for the new integration
  PBMC.TN138_QC <- scRNAQC(PBMC.list[[1]],FileName = paste0(Version,"/PBMC_QC/PBMC.TN138_QC"))
  PBMC.TN139_QC <- scRNAQC(PBMC.list[[2]],FileName = paste0(Version,"/PBMC_QC/PBMC.TN139_QC"))
  PBMC.TN146_QC <- scRNAQC(PBMC.list[[3]],FileName = paste0(Version,"/PBMC_QC/PBMC.TN146_QC"))
  PBMC.TN148_QC <- scRNAQC(PBMC.list[[4]],FileName = paste0(Version,"/PBMC_QC/PBMC.TN148_QC"))

  #### Save RData ####
  save.image(paste0(Save.Path,"/02_Quality_Control.RData"))


##### 03 Combine different data sets after QC  #####
  PBMC.list_QC  <- c(PBMC.TN138_QC ,PBMC.TN139_QC,PBMC.TN146_QC,PBMC.TN148_QC)
  rm(PBMC.TN138_QC ,PBMC.TN139_QC,PBMC.TN146_QC,PBMC.TN148_QC)

  # normalize and identify variable features for each dataset independently
    set.seed(1) # Fix the seed
    PBMC.list_QC <- lapply(X = PBMC.list_QC, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

  # select features that are repeatedly variable across datasets for integration
    set.seed(1) # Fix the seed
    features <- SelectIntegrationFeatures(object.list = PBMC.list_QC)

  ## Perform integration
    set.seed(1) # Fix the seed
    PBMC.anchors <- FindIntegrationAnchors(object.list = PBMC.list_QC, anchor.features = features)
    # this command creates an 'integrated' data assay
    set.seed(1) # Fix the seed
    PBMC.combined <- IntegrateData(anchorset = PBMC.anchors)


  ## Check QC
    scRNAQC(PBMC.combined,AddMitInf = "No",CheckOnly="Yes",FileName = paste0(Version,"/PBMC_QC/PBMC_QC_Check"))

  #### Save RData ####
  save.image(paste0(Save.Path,"/03_Combine_different_data_sets_after_QC.RData"))

##### 04 Perform an integrated analysis #####
  # Run the standard workflow for visualization and clustering

  # # # !!
  # set.seed(1) # Fix the seed
  # all.genes <- rownames(PBMC.combined)
  # PBMC.combined <- ScaleData(PBMC.combined, features = all.genes)

  ## Issues: re_clustering in seurat v3
  ## https://github.com/satijalab/seurat/issues/1528
  # DefaultAssay(PBMC.combined) <- "RNA"

  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(PBMC.combined) <- "integrated"

  ##?
  set.seed(1) # Fix the seed
  PBMC.combined <- ScaleData(PBMC.combined, verbose = FALSE)

  # ## Run if use filter
  # set.seed(1) # Fix the seed
  # PBMC.combined <- FindVariableFeatures(PBMC.combined)


  ### RunPCA
  # set.seed(1) # Fix the seed
  # PBMC.combined <- RunPCA(PBMC.combined, npcs = 30, verbose = FALSE)
  set.seed(1) # Fix the seed
  PBMC.combined <- RunPCA(PBMC.combined, features = VariableFeatures(object = PBMC.combined))

  print(PBMC.combined[["pca"]], dims = 1:5, nfeatures = 5)

  pdf(
    file = paste0(setwd(getwd()),"/",Version,"/PBMC_PCA.pdf"),
    width = 10,  height = 8
  )
  VizDimLoadings(PBMC.combined, dims = 1:2, reduction = "pca")
  DimPlot(PBMC.combined, reduction = "pca")
  DimHeatmap(PBMC.combined, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(PBMC.combined, dims = 1:15, cells = 500, balanced = TRUE)
  DimHeatmap(PBMC.combined, dims = 16:30, cells = 500, balanced = TRUE)

  # # Determine the 'dimensionality' of the dataset
  # # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # # computation time
  # PBMC.combined <- JackStraw(PBMC.combined, num.replicate = 100)
  # PBMC.combined <- ScoreJackStraw(PBMC.combined, dims = 1:20)
  # JackStrawPlot(PBMC.combined, dims = 1:20)
  ElbowPlot(PBMC.combined, ndims = 50)
  dev.off()

  ElbowPlot(PBMC.combined, ndims = 50)

  ## Issues: RunUMAP causes R exit
  ## https://github.com/satijalab/seurat/issues/2259
  # The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
  # To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
  # This message will be shown once per session
  #### UMAP
  set.seed(1) # Fix the seed
  PBMC.combined <- RunUMAP(PBMC.combined, reduction = "pca", dims = 1:30)

  set.seed(1) # Fix the seed
  PBMC.combined <- FindNeighbors(PBMC.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  PBMC.combined <- FindClusters(PBMC.combined, resolution = 0.5)

  #### tSNE
  set.seed(1) # Fix the seed
  PBMC.combined <- RunTSNE(PBMC.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  PBMC.combined <- FindNeighbors(PBMC.combined, reduction = "pca", dims = 1:30)
  set.seed(1) # Fix the seed
  PBMC.combined <- FindClusters(PBMC.combined, resolution = 0.5)


  ## Visualization
  DimPlot(PBMC.combined, reduction = "umap", group.by = "sample") %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)

        pdf(
          file = paste0(setwd(getwd()),"/",Version,"/PBMC_nlDR_Cluster.pdf"),
          width = 10,  height = 8
        )

        DimPlot(PBMC.combined, reduction = "umap", group.by = "sample") %>%
                BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)
        DimPlot(PBMC.combined, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14)

        DimPlot(PBMC.combined, reduction = "umap", ncol = 2,split.by = "sample", label = TRUE, label.size = 4) %>%
                BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                               SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9)

        DimPlot(PBMC.combined, reduction = "umap", ncol = 2,split.by = "Cachexia", label = TRUE, label.size = 4) %>%
                BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)
        DimPlot(PBMC.combined, reduction = "umap", ncol = 2,split.by = "Sex", label = TRUE, label.size = 4) %>%
          BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                         LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)

        ## tSNE
        DimPlot(PBMC.combined, reduction = "tsne", group.by = "sample") %>%
          BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)

        dev.off()
        # graphics.off()

  rm(PBMC.TN138_QC, PBMC.TN139_QC, PBMC.TN146_QC, PBMC.TN148_QC,PBMC.combined_QCTry)


  ##### Meta Table  #####

    ## Before QC
    Meta.df <- data.frame(matrix(nrow = 0,ncol = 3))
    colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
    Meta.df[1,1] <- c("EO.M")  # TN138
    Meta.df[1,2] <- ncol(PBMC.list[[1]]@assays[["RNA"]]@counts)
    Meta.df[1,3] <- nrow(PBMC.list[[1]]@assays[["RNA"]]@counts)

    Meta.df[2,1] <- c("LO.M")  # TN139
    Meta.df[2,2] <- ncol(PBMC.list[[2]]@assays[["RNA"]]@counts)
    Meta.df[2,3] <- nrow(PBMC.list[[2]]@assays[["RNA"]]@counts)

    Meta.df[3,1] <- c("LO.F")  # TN146
    Meta.df[3,2] <- ncol(PBMC.list[[3]]@assays[["RNA"]]@counts)
    Meta.df[3,3] <- nrow(PBMC.list[[3]]@assays[["RNA"]]@counts)

    Meta.df[4,1] <- c("EO.F")  # TN148
    Meta.df[4,2] <- ncol(PBMC.list[[4]]@assays[["RNA"]]@counts)
    Meta.df[4,3] <- nrow(PBMC.list[[4]]@assays[["RNA"]]@counts)

    # Summary to Meta table
    Meta.df[5,1] <- c("Summary")
    Meta.df[5,2] <- ncol(PBMC.combined_Ori@assays[["RNA"]]@counts)
    Meta.df[5,3] <- nrow(PBMC.combined_Ori@assays[["RNA"]]@counts)

    ## After QC
    colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
    Meta.df[6,1] <- c("EO.M.QC")  # TN138
    Meta.df[6,2] <- ncol(PBMC.list_QC[[1]]@assays[["RNA"]]@counts)
    Meta.df[6,3] <- nrow(PBMC.list_QC[[1]]@assays[["RNA"]]@counts)

    Meta.df[7,1] <- c("LO.M.QC")  # TN139
    Meta.df[7,2] <- ncol(PBMC.list_QC[[2]]@assays[["RNA"]]@counts)
    Meta.df[7,3] <- nrow(PBMC.list_QC[[2]]@assays[["RNA"]]@counts)

    Meta.df[8,1] <- c("LO.F.QC")  # TN146
    Meta.df[8,2] <- ncol(PBMC.list_QC[[3]]@assays[["RNA"]]@counts)
    Meta.df[8,3] <- nrow(PBMC.list_QC[[3]]@assays[["RNA"]]@counts)

    Meta.df[9,1] <- c("EO.F.QC")  # TN148
    Meta.df[9,2] <- ncol(PBMC.list_QC[[4]]@assays[["RNA"]]@counts)
    Meta.df[9,3] <- nrow(PBMC.list_QC[[4]]@assays[["RNA"]]@counts)

    # Summary to Meta table
    Meta.df[10,1] <- c("Summary")
    Meta.df[10,2] <- ncol(PBMC.combined@assays[["RNA"]]@counts)
    Meta.df[10,3] <- nrow(PBMC.combined@assays[["RNA"]]@counts)


    write.table( Meta.df ,
                 file = paste0(Save.Path,"/PBMC_CellCount_Meta.tsv"),
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
  PBMC.markers <- FindAllMarkers(PBMC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  library("magrittr")
  library("dplyr")

  # https://github.com/satijalab/seurat/issues/2960
  # Filter the top markers and plot the heatmap
  top_NSet = 7
  PBMC.markers %>%
    group_by(cluster) %>%
    top_n(n = top_NSet, wt = avg_log2FC) -> top_N
  PBMC.combined <- ScaleData(PBMC.combined, verbose = FALSE)
  DoHeatmap(PBMC.combined, features = top_N$gene) + NoLegend()
  write.table(top_N, file=paste0(Save.Path,"/PBMC_ClusterMarker_top",top_NSet,"Gene.txt"),sep="\t", row.names=T
              , quote = FALSE)
  write.table(PBMC.markers, file=paste0(Save.Path,"/PBMC_ClusterMarker_AllGene.txt"),sep="\t", row.names=T
              , quote = FALSE)

    pdf(
      file = paste0(Save.Path,"/PBMC_Heatmap_Cluster_top",top_NSet,".pdf"),
      width = 10,  height = 8
    )
      DoHeatmap(PBMC.combined, features = top_N$gene,size = 2,angle = 60) +
        scale_fill_gradient2(low="#5283ff",mid ="white", high ="#ff5c5c") +
        theme(axis.text.y = element_text(size  = 5)) +
        theme(legend.position = "bottom" )

    dev.off()


  # --------------- Check specific tissue marker --------------- #

    pdf(
      file = paste0(setwd(getwd()),"/",Version,"/PBMC_nlDR_CTMarker.pdf"),
      width = 10,  height = 8
    )

      # PMID: 31771616 #!!!!!!
      FeaturePlot(PBMC.combined, features = c("Cd3d", "Cd4", "Cd8a", "Csf1r", "Foxp3", "S100a9"), min.cutoff = "q9",
                  ncol = 3, coord.fixed = 1)
      # T Cell: Cd3d;  CD4+ T Cell: Cd4; CD8+ T Cell: Cd8a; Macrophages: Csf1r; regulatory T cells(Treg): Foxp3; Neutrophils: S100a9

      # PMID: 34296197 #!!!!!
      FeaturePlot(PBMC.combined, features = c("Cd3d", "Cd3e", "Lyz1", "Lyz2","Clu","Cd79a","Ms4a1","Nkg7","Gzmb"), min.cutoff = "q9",
                  ncol = 3, coord.fixed = 1)
      # T Cell: Cd3d,Cd3e;  Macrophages: Lyz; Mast Cell: Clu; B Cell: Cd79a,Ms4a1; NK Cell: Nkg7,Gzmb

      # http://biocc.hrbmu.edu.cn/CellMarker/
      # Mast cell
      FeaturePlot(PBMC.combined, features = c("Cd117", "Cd25","Cd203c","Slc18a2","Kit","Fcer1a","Cd9"), min.cutoff = "q9", coord.fixed = 1)
      # PMID: 30356731
      FeaturePlot(PBMC.combined, features = c("Cd9"), min.cutoff = "q9", coord.fixed = 1)
      # PMID: 34296197 #!!!!!
      FeaturePlot(PBMC.combined, features = c("Cpa3"), min.cutoff = "q9", coord.fixed = 1)
      # https://www.panglaodb.se/markers.html?cell_type=%27Mast%20cells%27

      ## Mac
      # Macrophage-Markers
      # https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
      ## M0
      FeaturePlot(PBMC.combined, features = c("Cd68", "Adgre1","Cd14","Csf1r","Ly6c1",
                                              "Cx3cr1","Fcgr1a","Itgam","Mertk"), min.cutoff = "q9", coord.fixed = 1)

      ## M1
      # M1 http://biocc.hrbmu.edu.cn/CellMarker/
      FeaturePlot(PBMC.combined, features = c("Cd16","Cd32","Cd64","Cd68","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1)
      # M1 https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
      FeaturePlot(PBMC.combined, features = c("Marco","Nos2","Tlr2","Cd80","Cd86","Csf2",
                                              "Tnf","Il1b","Il6","Tlr4","Cxcl2","Ifng","Il1r1"), min.cutoff = "q9", coord.fixed = 1)
      FeaturePlot(PBMC.combined, features = c("Il1a","Il1b","Il6","Nos2","Tlr2","Tlr4","Cd80","Cd86"), min.cutoff = "q9", coord.fixed = 1)


      ## M2
      # M2 http://biocc.hrbmu.edu.cn/CellMarker/
       FeaturePlot(PBMC.combined, features = c("Chil3","Csf1r","Mrc1","Pparg","Arg1","Cd163","Clec10a","Clec7a",
                                              "Cd206","Cd209","Ccl18","Fizz1"), min.cutoff = "q9", coord.fixed = 1)
      # https://www.biocompare.com/Editorial-Articles/566347-A-Guide-to-Macrophage-Markers/
       FeaturePlot(PBMC.combined, features = c("Cd115", "Cd206", "Pparg", "Arg1", "Cd163", "Cd301",
                                               "Dectin-1", "Pdcd1lg2", "Fizz1"), min.cutoff = "q9", coord.fixed = 1)

       FeaturePlot(PBMC.combined, features = c("Chil3"), min.cutoff = "q9", coord.fixed = 1)



      ## Tumor associated macrophage(TAM)
      FeaturePlot(PBMC.combined, features = c("Ccr2","Csf1r","Marco","Pdl2","Cd40","Ccl2","Csf1","Cd16"),
                  min.cutoff = "q9", coord.fixed = 1)

      # Erythrocytes
      FeaturePlot(PBMC.combined, features = c("Hbb-bs"), min.cutoff = "q9", coord.fixed = 1)
      # Platelet
      FeaturePlot(PBMC.combined, features = c("Ppbp"), min.cutoff = "q9", coord.fixed = 1)

      ## Summary
      markers.to.plot <- c("Cd3d","Cd3e", "Cd4","Cd8a", "Csf1r", "Lyz2","Chil3","Il1b", "S100a9","Nkg7",
                           "Gzmb", "Cd79a", "Ms4a1","Clu","Hbb-bs","Ppbp")

      FeaturePlot(PBMC.combined, features = markers.to.plot, min.cutoff = "q9", coord.fixed = 1)
      # T Cell: Cd3d,Cd3e;  CD4+ T Cell: Cd4; CD8+ T Cell: Cd8a; Macrophages: Csf1r,Lyz,Chil3;  Neutrophils: S100a9;
      # NK Cell: Nkg7,Gzmb; B Cell: Cd79a,Ms4a1; Mast Cell: Clu; Erythrocytes: Hbb-bs; Platelet: Ppbp


    dev.off()

    #### Save RData ####
    save.image(paste0(Save.Path,"/05_Identify_conserved_cell_type_markers.RData"))

##### 06 Cell type annotation  #####
  # PBMC.combined.copy <- PBMC.combined

  ## CD4+T: CD4+T Cell; CD8+T: CD8+T Cell; T: T Cell; B: B Cell; Mac: Macrophages;
  ## Neu: Neutrophils; NK: NK Cell; Mast: Mast Cell; Ery: Erythrocytes;
  ## Thr: Thrombocytes
  PBMC.combined <- RenameIdents(PBMC.combined, `0` = "CD4+T", `1` = "B", `2` = "Mac3",
                                `3` = "Neu", `4` = "CD8+T", `5` = "CD8+T", `6` = "Mac2", `7` = "CD4+T",
                                `8` = "NK", `9` = "Neu",`10` = "Mast", `11` = "T", `12` = "Ery", `13` = "Mac1",
                                `14` = "B", `15` = "B", `16` = "Other1", `17` = "Other2", `18` = "Neu")


  PBMC.combined$celltype <- Idents(PBMC.combined)
  PBMC.combined_OriCT <- PBMC.combined
  PBMC.combined <- PBMC.combined[,!grepl("Other",PBMC.combined$celltype)]
  PBMC.combined$celltype <- factor(PBMC.combined$celltype,levels = unique(PBMC.combined$celltype))


  # Cell_Type_Order.set <- c("T", "CD4+T", "CD8+T", "B" , "Mac1", "Mac2", "Mac3", "Mast", "NK", "Neu", "Ery", "Other1", "Other2")
  Cell_Type_Order.set <- c("Mac1", "Mac2", "Mac3","Neu", "T", "CD4+T", "CD8+T", "NK", "B" , "Mast",  "Ery")
  # PBMC.combined$celltype <- factor(PBMC.combined$celltype,
  #                                 levels = Cell_Type_Order.set)


  # Idents(PBMC.combined) <- "celltype"

    ## Heatmap

    Heatmap_Color.lt <- list(low="#5283ff",mid ="#311669", high ="#ff5c5c")
    Heatmap_Color.lt <- list(low="#419e47",mid ="#311669", high ="#edff66")
    Heatmap_Color.lt <- list(low="#5283ff",mid ="white", high ="#ff5c5c")

      pdf(
          file = paste0(Save.Path,"/PBMC_Heatmap_CellType_top",top_NSet,".pdf"),
          width = 10,  height = 8
      )
        DoHeatmap(PBMC.combined, features = top_N$gene,size = 3,angle = 60) +
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
    # PBMC.combined <- ScaleData(PBMC.combined, verbose = FALSE)
    DoHeatmap(PBMC.combined, features = top_N$gene,size = 4,angle = 90) + NoLegend()
    DoHeatmap(PBMC.combined, features = top_N$gene,group.by = "celltype",size = 4,angle = 90) + NoLegend()

    # Color
    DoHeatmap(PBMC.combined, features = top_N$gene,group.by = "celltype",size = 3,angle = 60) +
              scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                                   mid = Heatmap_Color.lt[["mid"]],
                                   high = Heatmap_Color.lt[["high"]])


    DoHeatmap(PBMC.combined, features = top_N$gene,group.by = "seurat_clusters",size = 3,angle = 90) + NoLegend()

    pdf(
        file = paste0(Save.Path,"/PBMC_Heatmap_CellType_top",top_NSet,".pdf"),
        width = 10,  height = 8
    )
      DoHeatmap(PBMC.combined, features = top_N$gene,group.by = "celltype",size = 2,angle = 45) +
                scale_fill_gradient2(low = Heatmap_Color.lt[["low"]],
                                     mid = Heatmap_Color.lt[["mid"]],
                                     high = Heatmap_Color.lt[["high"]])  +
                theme(axis.text.y = element_text(size = 5)) +
                theme(legend.position = "bottom" )

    dev.off()

  ## UMAP tSNE
    DimPlot(PBMC.combined, label = TRUE) %>% BeautifyggPlot(.,LegPos = c(1, 0.5))

    DimPlot(PBMC.combined,group.by = "celltype",label.size = 7,label = TRUE,
            pt.size =2) %>% BeautifyUMAP(FileName = "_PBMC_nlDR_CellType")
    DimPlot(PBMC.combined,group.by = "sample",
            pt.size =0.5) %>% BeautifyUMAP(FileName = "_PBMC_nlDR_Sample")
    DimPlot(PBMC.combined,group.by = "seurat_clusters",label.size = 7, label = TRUE,
            pt.size =1) %>% BeautifyUMAP(FileName = "_PBMC_nlDR_Clusters")

      pdf(
         file = paste0(setwd(getwd()),"/",Version,"/PBMC_nlDR_CellType_Sup.pdf"),
         width = 10,  height = 8
      )
        ##


        DimPlot(PBMC.combined, reduction = "umap", ncol = 2,split.by = "sample", label = TRUE) %>%
          BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                         SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9)

        DimPlot(PBMC.combined, reduction = "umap", ncol = 2,split.by = "Cachexia", label = TRUE, label.size = 4) %>%
          BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                         LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)

        DimPlot(PBMC.combined, reduction = "umap", ncol = 2,split.by = "Sex", label = TRUE, label.size = 4) %>%
          BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                         LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1)

        ## tSNE
        DimPlot(PBMC.combined, reduction = "tsne", group.by = "sample") %>%
          BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18)

      dev.off()


  ## DotPlot
    DotPlot_Color1.set <- c("#de3767", "#de3767", "#4169e1", "#4169e1")
    DotPlot_Color2.set <- c("#5b8e7d","#7b2cbf")
    DotPlot_Color3.set <- c("#de3767", "#4169e1")

    pdf(
      file = paste0(Save.Path,"/PBMC_DotPlot_CellType",".pdf"),
      width = 10,  height = 8
    )

      # https://satijalab.org/seurat/reference/dotplot
      DotPlot(PBMC.combined, features = markers.to.plot, cols = c("lightgrey", "blue"),
              dot.scale = 8) + RotatedAxis()%>%
        BeautifyggPlot(.,LegPos = "bottom",AxisTitleSize=1, TitleSize = 20, xangle =90,
                       LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1,XtextSize=12,  YtextSize=12)

      DotPlot(PBMC.combined, features = markers.to.plot, cols = DotPlot_Color1.set,
              dot.scale = 8, split.by = "sample") + RotatedAxis()


      # https://github.com/satijalab/seurat/issues/1541
      DotPlot(PBMC.combined, features = markers.to.plot, cols = DotPlot_Color2.set,
              dot.scale = 8, split.by = "Cachexia") + RotatedAxis()

      DotPlot(PBMC.combined, features = markers.to.plot, cols = DotPlot_Color3.set,
              dot.scale = 8, split.by = "Sex") + RotatedAxis()

    dev.off()

    rm(top_N, top_NSet)

    #### Save RData ####
    save.image(paste0(Save.Path,"/06_Cell_type_annotation.RData"))

  # ##### Export marker gene from specific cluster #####
  #   # For performing differential expression after integration, we switch back to the original data
  #   set.seed(1) # Fix the seed
  #   DefaultAssay(PBMC.combined) <- "RNA"
  #
  #   # nk.markers <- FindConservedMarkers(PBMC.combined, ident.1 = 6, grouping.var = "sample", verbose = FALSE)
  #   library(BiocManager)
  #   library(multtest)
  #   nk.markers <- FindConservedMarkers(PBMC.combined, ident.1 = 'NK', grouping.var = "sample", verbose = FALSE)
  #   head(nk.markers)
  #
  #   rm(nk.markers)

  # ##### Identify differential expressed genes across conditions  #####
  #   library(ggplot2)
  #   library(cowplot)
  #   theme_set(theme_cowplot())
  #   CD4T.cells <- subset(PBMC.combined, idents = "CD4+T")
  #   Idents(CD4T.cells) <- "Cachexia"
  #   avg.CD4T.cells <- as.data.frame(log1p(AverageExpression(CD4T.cells, verbose = FALSE)$RNA))
  #   avg.CD4T.cells$gene <- rownames(avg.CD4T.cells)
  #
  #   MacrophageM2 <- subset(PBMC.combined, idents = "Mac2")
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





###########################################################################################

##### Interested genes #####

  FeaturePlot(PBMC.combined, features = c("Top2a", "Ptk2"), min.cutoff = "q9")
  Test <- intersect(row.names(PBMC.combined@assays[["RNA"]]@counts) ,  c("Top2a", "Ptk2"))
  Test2 <- intersect(row.names(PBMC.combined@assays[["integrated"]]@data) ,  c("Top2a", "Ptk2"))
  FeaturePlot(PBMC.combined, features = c("Kras", "Braf", "Myc","Egfr","Pik3ca"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(PBMC.combined, features = c("Il4"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(PBMC.combined, features = c("Il4ra"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(PBMC.combined, features = c("Il13ra1"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(PBMC.combined, features = c("Stat6"), min.cutoff = "q9",ncol = 3)

# -------------------------------------- #

###########################################################################################



# ## Export Seurat Object in 10X format
# # https://github.com/satijalab/seurat/issues/884
# library(DropletUtils)
# write10xCounts(PBMC.list_QC_Try[["TN138"]]@assays[["RNA"]]@counts, path = paste0(setwd(getwd()),"/Test"))
