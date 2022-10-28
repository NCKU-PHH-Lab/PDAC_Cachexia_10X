##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages  #####
library("tidyverse")
library(Seurat)


##### Load RData #####
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_PBMC_Main/06_Cell_type_annotation.RData")
# Clean up
rm(list=setdiff(ls(), c("PBMC.combined")))
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_SC_Main/06_Cell_type_annotation.RData")
rm(list=setdiff(ls(), c("PBMC.combined","SC.combined")))



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


PBMC.combined$celltype <- paste0("PBMC_",PBMC.combined$celltype)
PBMC.combined <- PBMC.combined[,grepl("Mac", PBMC.combined$celltype, ignore.case=TRUE)]


SC.combined$celltype <- paste0("SC_",SC.combined$celltype)
SC.combined <- SC.combined[,grepl("Mac", SC.combined$celltype, ignore.case=TRUE)]

# select features that are repeatedly variable across datasets for integration
set.seed(1) # Fix the seed
features <- SelectIntegrationFeatures(object.list = c(PBMC.combined,SC.combined))
## Perform integration
set.seed(1) # Fix the seed
SC.anchors <- FindIntegrationAnchors(object.list = c(PBMC.combined,SC.combined), anchor.features = features)
# this command creates an 'integrated' data assay
set.seed(1) # Fix the seed
scRNA.SeuObj <- IntegrateData(anchorset = SC.anchors)


##### 04 Perform an integrated analysis #####
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(scRNA.SeuObj) <- "integrated"

## ScaleData
set.seed(1) # Fix the seed
scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)

# ## Run if use filter
# set.seed(1) # Fix the seed
# scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)

### RunPCA
# set.seed(1) # Fix the seed
# scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 30, verbose = FALSE)
set.seed(1) # Fix the seed
# scRNA.SeuObj <- RunPCA(scRNA.SeuObj, features = VariableFeatures(object = scRNA.SeuObj))
scRNA.SeuObj <- RunPCA(scRNA.SeuObj)

print(scRNA.SeuObj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(scRNA.SeuObj, dims = 1:2, reduction = "pca")
DimPlot(scRNA.SeuObj, reduction = "pca")
DimHeatmap(scRNA.SeuObj, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(PBMC.combined, ndims = 50)

#### UMAP
set.seed(1) # Fix the seed
scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:30)

set.seed(1) # Fix the seed
scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:30)
set.seed(1) # Fix the seed
scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

DimPlot(scRNA.SeuObj, reduction = "umap", group.by = "sample") %>%
  BeautifyggPlot(.,TV= -5,TitleSize = 25,LegPos = c(0.85, 0.15),AxisTitleSize=1.2, LegTextSize = 18) %>% print()
DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, label.size = 7, repel = TRUE) %>%
  BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14) %>% print()

DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "sample", label = TRUE, label.size = 4) %>%
  BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, TitleSize = 20,
                 SubTitSize = 17, LegTextSize = 14, XaThick=0.9, YaThick=0.9) %>% print()

DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "Cachexia", label = TRUE, label.size = 4) %>%
  BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                 LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1) %>% print()
DimPlot(scRNA.SeuObj, reduction = "umap", ncol = 2,split.by = "Sex", label = TRUE, label.size = 4) %>%
  BeautifyggPlot(.,LegPos = "top",AxisTitleSize=1.2, TitleSize = 20,
                 LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1) %>% print()


DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype", label = TRUE, label.size = 7, repel = TRUE) %>%
  BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14) %>% print()




# scRNA_Sub.SeuObj <- scRNA.SeuObj[,str_subset(scRNA.SeuObj$celltype, pattern = "Mac")]
scRNA_Sub.SeuObj <- scRNA.SeuObj[,grepl("Mac", scRNA.SeuObj$celltype, ignore.case=TRUE)]
DimPlot(scRNA_Sub.SeuObj, reduction = "umap",group.by = "celltype", label = TRUE, label.size = 7, repel = TRUE) %>%
  BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14) %>% print()

library("scales")
library("ggsci")
show_col(pal_npg("nrc")(10))

DimPlot(scRNA_Sub.SeuObj, reduction = "umap",group.by = "celltype", pt.size =2,
        label = TRUE,label.size = 6,label.color="#1c1c19", repel = TRUE,cols =	pal_lancet("lanonc", alpha = 0.2)(10)) %>%
  BeautifyggPlot(.,LegPos = "bottom",AxisTitleSize=1.2, TitleSize = 20,
                 LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1) %>% print()


scRNA_Sub.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$celltype %in% c("SC_Mac1","SC_Mac2", "PBMC_Mac3")]
DimPlot(scRNA_Sub.SeuObj, reduction = "umap",group.by = "celltype", label = TRUE, label.size = 7, repel = TRUE) %>%
  BeautifyggPlot(.,LegPos = c(1, 0.5),AxisTitleSize=1.2, LegTextSize = 14) %>% print()



FeaturePlot(scRNA.SeuObj, features = c("Top2a", "Ptk2"), min.cutoff = "q9")



##### Count Metadata #####
Metadata.df <- scRNA.SeuObj@meta.data
