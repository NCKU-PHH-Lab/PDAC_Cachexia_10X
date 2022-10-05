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
  Version = paste0(Sys.Date(),"_","SC_Main")
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path)

##### Load datasets*  #####
  SC.data.TN136 <- Read10X(data.dir = paste0(getwd(),"/TN136/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN136 <- CreateSeuratObject(counts = SC.data.TN136, project = "EO.M", min.cells = 3, min.features = 200)
  SC.TN136
  SC.TN136@meta.data[["sample"]] <- rep(c("EO.M"), times=length(SC.TN136@meta.data[["orig.ident"]]))
  SC.TN136@meta.data[["ID"]] <- rep(c("TN136"), times=length(SC.TN136@meta.data[["orig.ident"]]))
  SC.TN136@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(SC.TN136@meta.data[["orig.ident"]]))  #EO: Early_Onset
  SC.TN136@meta.data[["Sex"]] <- rep(c("Male"), times=length(SC.TN136@meta.data[["orig.ident"]]))

  SC.data.TN137 <- Read10X(data.dir = paste0(getwd(),"/TN137/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN137 <- CreateSeuratObject(counts = SC.data.TN137, project = "LO.M", min.cells = 3, min.features = 200)
  SC.TN137
  SC.TN137@meta.data[["sample"]] <- rep(c("LO.M"), times=length(SC.TN137@meta.data[["orig.ident"]]))
  SC.TN137@meta.data[["ID"]] <- rep(c("TN137"), times=length(SC.TN137@meta.data[["orig.ident"]]))
  SC.TN137@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(SC.TN137@meta.data[["orig.ident"]]))  #LO: Late_Onset
  SC.TN137@meta.data[["Sex"]] <- rep(c("Male"), times=length(SC.TN137@meta.data[["orig.ident"]]))

  SC.data.TN145 <- Read10X(data.dir = paste0(getwd(),"/TN145/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN145 <- CreateSeuratObject(counts = SC.data.TN145, project = "LO.F", min.cells = 3, min.features = 200)
  SC.TN145
  SC.TN145@meta.data[["sample"]] <- rep(c("LO.F"), times=length(SC.TN145@meta.data[["orig.ident"]]))
  SC.TN145@meta.data[["ID"]] <- rep(c("TN145"), times=length(SC.TN145@meta.data[["orig.ident"]]))
  SC.TN145@meta.data[["Cachexia"]] <- rep(c("LO"), times=length(SC.TN145@meta.data[["orig.ident"]]))
  SC.TN145@meta.data[["Sex"]] <- rep(c("Female"), times=length(SC.TN145@meta.data[["orig.ident"]]))

  SC.data.TN147 <- Read10X(data.dir = paste0(getwd(),"/TN147/monocle/outs/filtered_gene_bc_matrices/mm10"))
  SC.TN147 <- CreateSeuratObject(counts = SC.data.TN147, project = "EO.F", min.cells = 3, min.features = 200)
  SC.TN147
  SC.TN147@meta.data[["sample"]] <- rep(c("EO.F"), times=length(SC.TN147@meta.data[["orig.ident"]]))
  SC.TN147@meta.data[["ID"]] <- rep(c("TN147"), times=length(SC.TN147@meta.data[["orig.ident"]]))
  SC.TN147@meta.data[["Cachexia"]] <- rep(c("EO"), times=length(SC.TN147@meta.data[["orig.ident"]]))
  SC.TN147@meta.data[["Sex"]] <- rep(c("Female"), times=length(SC.TN147@meta.data[["orig.ident"]]))


##### 01 Combine different datasets before QC  #####
  SC.list  <- c(SC.TN136,SC.TN137,SC.TN145,SC.TN147)
  rm(SC.TN136,SC.TN137,SC.TN145,SC.TN147, SC.data.TN136,SC.data.TN137,SC.data.TN145, SC.data.TN147)

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

  ## ScaleData
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
    #### Before QC ####
    Meta.df <- data.frame(matrix(nrow = 0,ncol = 3))
    colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
    Meta.df[1,1] <- c("EO.M")  # TN136
    Meta.df[1,2] <- ncol(SC.list[[1]]@assays[["RNA"]]@counts)
    Meta.df[1,3] <- nrow(SC.list[[1]]@assays[["RNA"]]@counts)

    Meta.df[2,1] <- c("LO.M")  # TN137
    Meta.df[2,2] <- ncol(SC.list[[2]]@assays[["RNA"]]@counts)
    Meta.df[2,3] <- nrow(SC.list[[2]]@assays[["RNA"]]@counts)

    Meta.df[3,1] <- c("LO.F")  # TN145
    Meta.df[3,2] <- ncol(SC.list[[3]]@assays[["RNA"]]@counts)
    Meta.df[3,3] <- nrow(SC.list[[3]]@assays[["RNA"]]@counts)

    Meta.df[4,1] <- c("EO.F")  # TN147
    Meta.df[4,2] <- ncol(SC.list[[4]]@assays[["RNA"]]@counts)
    Meta.df[4,3] <- nrow(SC.list[[4]]@assays[["RNA"]]@counts)

    # Summary to Meta table
    Meta.df[5,1] <- c("Summary")
    Meta.df[5,2] <- ncol(SC.combined_Ori@assays[["RNA"]]@counts)
    Meta.df[5,3] <- nrow(SC.combined_Ori@assays[["RNA"]]@counts)

    #### After QC ####
    colnames(Meta.df) <- c("NO.","Cell_Num","Gene_Num")
    Meta.df[6,1] <- c("EO.M.QC")  # TN136
    Meta.df[6,2] <- ncol(SC.list_QC[[1]]@assays[["RNA"]]@counts)
    Meta.df[6,3] <- nrow(SC.list_QC[[1]]@assays[["RNA"]]@counts)

    Meta.df[7,1] <- c("LO.M.QC")  # TN137
    Meta.df[7,2] <- ncol(SC.list_QC[[2]]@assays[["RNA"]]@counts)
    Meta.df[7,3] <- nrow(SC.list_QC[[2]]@assays[["RNA"]]@counts)

    Meta.df[8,1] <- c("LO.F.QC")  # TN145
    Meta.df[8,2] <- ncol(SC.list_QC[[3]]@assays[["RNA"]]@counts)
    Meta.df[8,3] <- nrow(SC.list_QC[[3]]@assays[["RNA"]]@counts)

    Meta.df[9,1] <- c("EO.F.QC")  # TN147
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
      ## Summary
      markers.to.plot <- c("Sox9", "Csf1r","Lyz","Chil3","Col1a1","Col3a1","Dcn")

      FeaturePlot(SC.combined, features = markers.to.plot, min.cutoff = "q9", coord.fixed = 1 ,ncol =3)

      # Duc : "Sox9"; Macrophages: "Csf1r","Lyz","Chil3"; Fib: "Col1a1","Col3a1","Dcn"


      FeaturePlot(SC.combined, features = c("Krt19", "Prss1", "Chg8", "Cdh5", "Lum", "Rgs5", "Aif1",
                                            "Cd3d", "Ms4a1","Ms4a2"), min.cutoff = "q9",
                  ncol = 3, coord.fixed = 1)
      # Ductal cells
      # https://pubmed.ncbi.nlm.nih.gov/24418153/
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

      ## Macrophages
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


      ## CAF
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
  # SC.combined$celltype <- factor(SC.combined$celltype,
  #                                  levels = Cell_Type_Order.set)


  SC.combined$celltype <- Idents(SC.combined)
  # Idents(SC.combined) <- "celltype"

    ## Heatmap
    # Heatmap_Color.lt <- list(low="#5283ff",mid ="#311669", high ="#ff5c5c")
    # Heatmap_Color.lt <- list(low="#419e47",mid ="#311669", high ="#edff66")
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



###########################################################################################
# -------------------------------------- #
##### Interested genes #####

  FeaturePlot(SC.combined, features = c("Top2a", "Ptk2"), min.cutoff = "q9")
  Test <- intersect(row.names(SC.combined@assays[["RNA"]]@counts) ,  c("Top2a", "Ptk2"))
  Test2 <- intersect(row.names(SC.combined@assays[["integrated"]]@data) ,  c("Top2a", "Ptk2"))
  FeaturePlot(SC.combined, features = c("Kras", "Braf", "Myc","Egfr","Pik3ca"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(SC.combined, features = c("Il4"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(SC.combined, features = c("Il4ra"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(SC.combined, features = c("Il13ra1"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(SC.combined, features = c("Stat6"), min.cutoff = "q9",ncol = 3)
  FeaturePlot(SC.combined, features = c("Il1b"), min.cutoff = "q9",ncol = 3)


  pdf(
    file = paste0(Save.Path,"/SC_UMAP_SubCTMarkers",".pdf"),
    width = 10,  height = 8
  )

    # https://pubmed.ncbi.nlm.nih.gov/31197017/
    # iCAF
    FeaturePlot(SC.combined, features = c("Clec3b","Col14a1","Has1","Il6"), min.cutoff = "q9",ncol = 2)
    # apCAF
    FeaturePlot(SC.combined, features = c("H2-Ab1","Cd74","Saa3","Slpi"), min.cutoff = "q9",ncol = 2)
    # myCAF
    FeaturePlot(SC.combined, features = c("Tagln","Thy1","Col12a1","Thbs2"), min.cutoff = "q9",ncol = 2)

    # EMT-like cells
    FeaturePlot(SC.combined, features = c("Cdkn2a","S100a6","Igfbp4","Sparc","Vim","Spp1"), min.cutoff = "q9",ncol = 2)
  dev.off()

  pdf(
    file = paste0(Save.Path,"/SC_UMAP_SubCTMarkers2",".pdf"),
    width = 10,  height = 14
  )

    # EMT-like cells
    FeaturePlot(SC.combined, features = c("Cdkn2a","S100a6","Igfbp4","Sparc","Vim","Spp1"), min.cutoff = "q9",ncol = 2)
  dev.off()

  rm(Test,Test2)
# -------------------------------------- #

###########################################################################################

# ## Export Seurat Object in 10X format
# # https://github.com/satijalab/seurat/issues/884
# library(DropletUtils)
# write10xCounts(SC.list_QC_Try[["TN136"]]@assays[["RNA"]]@counts, path = paste0(setwd(getwd()),"/Test"))
