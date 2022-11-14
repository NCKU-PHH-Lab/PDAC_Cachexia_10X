

##### Load Data #####
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_SC_Main/06_Cell_type_annotation.RData")

## INTCHG: Interchangeable
## SubType Setting
if(SampleType == "PBMC"){
  # For PBMC
  scRNA.SeuObj <- PBMC.combined

  # Order the cell type
  CellType.Order = c("Mac1", "Mac2","Mac3","Neu","T","CD4+T","CD8+T","NK","B","Mast","Ery")
  scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
                                                 levels = CellType.Order)


}else if(SampleType == "SC"){
  # For SC
  scRNA.SeuObj <- SC.combined

  # Order the cell type
  CellType.Order = c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" , "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
                     "Fib1", "Fib2", "Fib3")
  scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
                                                 levels = CellType.Order)

}


rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Save.Path","CellType.Order")))

##### Load Package #####


library(Seurat)
library(tidyverse)
library(patchwork)

source("FUN_HSsymbol2MMsymbol.R")
source("FUN_Beautify_ggplot.R")
source("FUN_Beautify_UMAP.R")


##### Function (tsv)(New version of GSEA) #####
FUN_GeneSetAnno <- function(scRNA.SeuObj,
                            Set_GeneSet = "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.v2022.1.Mm.tsv",
                            Set_Path = paste0(getwd(),"/GSEA_Geneset/"))
{
  ## Load geneset
  GeneSet_Ori.set <- read.delim(paste0(Set_Path,Set_GeneSet),header=T, stringsAsFactors = FALSE)

  ## Clean up geneset
  GeneSet.df <- GeneSet_Ori.set[nrow(GeneSet_Ori.set),]
  GeneSetName <- colnames(GeneSet.df)[2]
  GeneSet.set <- GeneSet.df[,2]

  library(tidyr)
  GeneSet.df <- data.frame(strsplit(as.character(GeneSet.set), ","))
  colnames(GeneSet.df) <- GeneSetName


  ## Convert Gene name
  GeneSet_Species <- "Human" # c("Mouse","Human")

  if(GeneSet_Species == "Mouse"){
    GeneSet.df_Mouse <- GeneSet.df[,1]
  }else if(GeneSet_Species == "Human"){
    GeneSet.df_Mouse <- HSsymbol2MMsymbol(GeneSet.df,GeneSetName)
    GeneSet.df_Mouse <- unique(GeneSet.df_Mouse[,2])
  }

  # FeaturePlot(scRNA.SeuObj, features = GeneSet.df_Mouse[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")


  ## AddModuleScore
  # DefaultAssay(scRNA.SeuObj) <- "integrated"
  DefaultAssay(scRNA.SeuObj) <- "RNA"
  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(GeneSet.df_Mouse),
                                 ctrl = 100, name = GeneSetName)

  colnames(scRNA.SeuObj@meta.data)[ncol(scRNA.SeuObj@meta.data)] <- GeneSetName


  ## Plot UMAP
  plot.UMAP <- FeaturePlot(object = scRNA.SeuObj, features = paste0(GeneSetName), pt.size =1)+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           limits = c(-max(abs(scRNA.SeuObj@meta.data[,paste0(GeneSetName)])), max(abs(scRNA.SeuObj@meta.data[,paste0(GeneSetName)]))), # max(abs(scRNA.SeuObj@meta.data[,paste0(GeneSetName,"1")]))
                           # limits = c(-0.2, 0.2), # max(scRNA.SeuObj@meta.data[,paste0(GeneSetName)])
                           guide = "colourbar",midpoint = 0, labs(fill ="Score"))
  source("FUN_Beautify_UMAP.R")
  plot.UMAP %>% Beautify_UMAP(TitleSize = 15, LegTitleSize = 15, LegTextSize = 13) -> plot.UMAP
  plot.UMAP

  # FeaturePlot(object = scRNA.SeuObj, features = paste0(GeneSetName), split.by = "sample",
  #             cols = c("#0077b6", "grey",  "#ef476f"), pt.size =1)
  #
  # FeaturePlot(object = scRNA.SeuObj, features = paste0(GeneSetName), split.by = "Cachexia",
  #             cols = c("#0077b6", "grey",  "#ef476f"), pt.size =1)
  # FeaturePlot(object = scRNA.SeuObj, features = paste0(GeneSetName), split.by = "Cachexia")+
  #   scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
  #                          guide = "colourbar",midpoint = 0, labs(fill ="Exp"))


  ## Plot Violin

  plot.Violin <- VlnPlot(scRNA.SeuObj, features = paste0(GeneSetName), split.by = "Cachexia", group.by = "celltype",
                         pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  plot.Violin <- wrap_plots(plot.Violin = plot.Violin, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.12, -0.07),AxisTitleSize=1.7,
                                                                                     XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                                     LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  plot.Violin

  # plot.Violin1 <- VlnPlot(scRNA.SeuObj, features = paste0(GeneSetName), split.by = "sample", group.by = "celltype",
  #                   pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  # wrap_plots(plot.Violin1 = plot.Violin1, ncol = 1)

  ## Output setting
  Output.lt <- list(scRNA.SeuObj = scRNA.SeuObj,
                    UMAP = plot.UMAP,
                    Violin = plot.Violin)

  return(Output.lt)

}


## Example
AnnoResult.lt <-FUN_GeneSetAnno(scRNA.SeuObj,
                                Set_GeneSet = "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.v2022.1.Mm.tsv",
                                Set_Path = paste0(getwd(),"/GSEA_Geneset/"))

AnnoResult.lt <-FUN_GeneSetAnno(scRNA.SeuObj,
                                Set_GeneSet = "REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION.v2022.1.Mm.tsv",
                                Set_Path = paste0(getwd(),"/GSEA_Geneset/CATs/"))


scRNA.SeuObj <- AnnoResult.lt$scRNA.SeuObj
plot.UMAP <- AnnoResult.lt$UMAP
plot.UMAP

plot.Violin <- AnnoResult.lt$Violin
plot.Violin

##### Run Multiple GeneSet  #####
library(gtools)
# target.dir <- list.dirs(InputFolder)[-1]

GeneSetType <- "EMT"
InputFolder = paste0(getwd(),"/GSEA_Geneset/",GeneSetType,"/")
list.files <- list.files(InputFolder,full.names = F)
Nfiles = length(list.files)

for(i in 1:Nfiles)
{
  AnnoResult.lt <-FUN_GeneSetAnno(scRNA.SeuObj,
                                  Set_GeneSet = list.files[i],
                                  Set_Path = InputFolder)
  scRNA.SeuObj <- AnnoResult.lt$scRNA.SeuObj
}


rm(AnnoResult.lt,i)

MetaData <- scRNA.SeuObj@meta.data
EMT.df <- MetaData[,(ncol(MetaData)-Nfiles+1):ncol(MetaData)]



##### ComplexHeatmap #####
## https://jokergoo.github.io/ComplexHeatmap-reference/book/
library(ComplexHeatmap)
library(circlize)

Duc_Anno.df <- MetaData[grepl("Duc",MetaData$celltype), -(ncol(MetaData)-Nfiles+1):-ncol(MetaData)]
Duc_Anno.df <- data.frame(CellID = rownames(Duc_Anno.df), Duc_Anno.df)
Duc_Anno.df$celltype <- factor(Duc_Anno.df$celltype ,
                               levels = c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6"))


EMT.df <- data.frame(CellID = rownames(EMT.df), EMT.df)

Duc_EMT.df <- left_join(data.frame(CellID = Duc_Anno.df[,1]),EMT.df)



rownames(Duc_EMT.df) <- Duc_EMT.df[,1]
Duc_EMT.df <- Duc_EMT.df[,-1]


# ## Try Basic Heatmap
# Heatmap(Duc_EMT.df %>% t())




#### Set column annotation ####
sample = c("#2267a4", "#3d85c6", "#d5a6bd", "#c27ba0")
names(sample) <- Duc_Anno.df$sample %>% unique()

library(ggsci)
library(ggplot2)
# vignette( "ggsci") #Check the color setting
col3 = 	pal_d3("category20", alpha = 0.7)(length(CellType.Order))
colCT <- col3[1:length(CellType.Order)]
names(colCT) <- CellType.Order

ha_column_T = HeatmapAnnotation(
  Sample = Duc_Anno.df[,"sample"],  # anno_colum.df$sample_type,
  Cachexia = Duc_Anno.df[,"Cachexia"], # anno_colum.df$gender,
  Celltype = Duc_Anno.df[,"celltype"],
  col = list(Sample = sample,
             Celltype = colCT ,#pal_npg(), #colCT ,
             Cachexia = c("EOCX"="#5b517d", "PreCX"="#a095c7")), #,"Medium"="#b57545"
  show_legend = T,
  show_annotation_name = F
)




#### Plot Heatmap ####
# Set Heatmap color
# col_HMap <- c("#f0f4fc", "#6e8cc2", "#37558c")
col_HMap <- c( "#ffffff","#e697b3", "#c45e82")
col_HMap <- c( "#ffffff","#d44c5f", "#520e17")

### Plat GeneExpression Heatmap
## Reorder Heatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#row-and_column_orders
Heatmap(
  Duc_EMT.df %>% t(),
  cluster_rows = F, # Heatmap with/without clustering by rows
  cluster_columns = F, # Heatmap with/without clustering by columns
  column_order = order(Duc_Anno.df$celltype,Duc_Anno.df$Cachexia), ## Reorder Heatmap
  show_column_names = F,
  show_row_names = T,
  name = "Score",
  # # set color
  # col = colorRamp2(
  #   # c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
  #   c(min(GeneExp.df), max(GeneExp.df)*1/3 , max(GeneExp.df)),
  #
  #   col_HMap
  # ),
  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T,
  # right_annotation = ha_row
  # width =  length(CellType.Order)*unit(7, "mm"),
  # height = nrow( Duc_EMT.df)*unit(7, "mm"), # length(CellType.Order)*unit(15, "mm"),
  column_title = paste0(SampleType,"_",GeneSetType),
  column_title_gp = gpar(fontsize = 20, fontface = "bold")
) -> P.Heatmap_GeneExp

P.Heatmap_GeneExp %>% print







# ##### Old Version #####
# ##### ******************************************************************** #####
#
#
# ##### HALLMARK_TNFA_SIGNALING_VIA_NFKB #####
#   HALLMARK_TNFA_SIGNALING_VIA_NFKB <- read.delim(paste0(getwd(),"/GSEA_Geneset/HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt"),header=T, stringsAsFactors = FALSE)
#   HALLMARK_TNFA_SIGNALING_VIA_NFKB <- data.frame(HALLMARK_TNFA_SIGNALING_VIA_NFKB[-1,])
#   colnames(HALLMARK_TNFA_SIGNALING_VIA_NFKB) <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
#
#   HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse <- HSsymbol2MMsymbol(HALLMARK_TNFA_SIGNALING_VIA_NFKB,"HALLMARK_TNFA_SIGNALING_VIA_NFKB")
#   HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2 <- unique(HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2),
#                                   ctrl = 5, name = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')
#   FeaturePlot(object = scRNA.SeuObj, features = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB1', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB1', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots <- VlnPlot(scRNA.SeuObj, features = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB1"), split.by = "sample", group.by = "celltype",
#                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots = plots, ncol = 1)
#
#   plots <- VlnPlot(scRNA.SeuObj, features = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB1"), split.by = "Cachexia", group.by = "celltype",
#                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p1 <- wrap_plots(plots = plots, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                         XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                         LegTextSize = 15)  +
#                    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p1
#
#   ###
# ##### REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION #####
#   REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION.txt"),header=T, stringsAsFactors = FALSE)
#   REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION <- data.frame(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION[-1,])
#   colnames(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION) <- c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION")
#
#   REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse <- HSsymbol2MMsymbol(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION,"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION")
#   REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2 <- unique(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2),
#                                   ctrl = 5, name = 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION')
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots2 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1"), split.by = "sample", group.by = "celltype",
#                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots2 = plots2, ncol = 1)
#
#   plots2 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1"), split.by = "Cachexia", group.by = "celltype",
#                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p2 <- wrap_plots(plots2 = plots2, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                                 XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                                 LegTextSize = 15)  +
#     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p2
#
#
# ##### REACTOME_ANTIMICROBIAL_PEPTIDES #####
#   REACTOME_ANTIMICROBIAL_PEPTIDES <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_ANTIMICROBIAL_PEPTIDES.txt"),header=T, stringsAsFactors = FALSE)
#   REACTOME_ANTIMICROBIAL_PEPTIDES <- data.frame(REACTOME_ANTIMICROBIAL_PEPTIDES[-1,])
#   colnames(REACTOME_ANTIMICROBIAL_PEPTIDES) <- c("REACTOME_ANTIMICROBIAL_PEPTIDES")
#
#   REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse <- HSsymbol2MMsymbol(REACTOME_ANTIMICROBIAL_PEPTIDES,"REACTOME_ANTIMICROBIAL_PEPTIDES")
#   REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2 <- unique(REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2),
#                                   ctrl = 5, name = 'REACTOME_ANTIMICROBIAL_PEPTIDES')
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_ANTIMICROBIAL_PEPTIDES1', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_ANTIMICROBIAL_PEPTIDES1', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots3 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_ANTIMICROBIAL_PEPTIDES1"), split.by = "sample", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots3 = plots3, ncol = 1)
#
#   plots3 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_ANTIMICROBIAL_PEPTIDES1"), split.by = "Cachexia", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p3 <- wrap_plots(plots3 = plots3, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                                    XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                                    LegTextSize = 15)  +
#     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p3
#
#
# ##### REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION #####
#   REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION.txt"),header=T, stringsAsFactors = FALSE)
#   REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION <- data.frame(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION[-1,])
#   colnames(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION) <- c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION")
#
#   REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse <- HSsymbol2MMsymbol(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION,"REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION")
#   REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2 <- unique(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2),
#                                   ctrl = 5, name = 'REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION')
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots4 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1"), split.by = "sample", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots4 = plots4, ncol = 1)
#
#   plots4 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1"), split.by = "Cachexia", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p4 <- wrap_plots(plots4 = plots4, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                                    XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                                    LegTextSize = 15)  +
#     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p4
#
# ##### REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 #####
#   REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2.txt"),header=T, stringsAsFactors = FALSE)
#   REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 <- data.frame(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2[-1,])
#   colnames(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2) <- c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2")
#
#   REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse <- HSsymbol2MMsymbol(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2,"REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2")
#   REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2 <- unique(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2),
#                                   ctrl = 5, name = 'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2')
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots5 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21"), split.by = "sample", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots5 = plots5, ncol = 1)
#
#   plots5 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21"), split.by = "Cachexia", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p5 <- wrap_plots(plots5 = plots5, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                                    XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                                    LegTextSize = 15)  +
#     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p5
#
#
# ##### REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS #####
#   REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS.txt"),header=T, stringsAsFactors = FALSE)
#   REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS <- data.frame(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS[-1,])
#   colnames(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS) <- c("REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS")
#
#   REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse <- HSsymbol2MMsymbol(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS,"REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS")
#   REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2 <- unique(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2),
#                                   ctrl = 5, name = 'REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS')
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots6 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1"), split.by = "sample", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots6 = plots6, ncol = 1)
#
#   plots6 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1"), split.by = "Cachexia", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p6 <- wrap_plots(plots6 = plots6, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                                   XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                                   LegTextSize = 15)  +
#     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p6
#
# ##### KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION #####
#   KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION <- read.delim(paste0(getwd(),"/GSEA_Geneset/KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION.txt"),header=T, stringsAsFactors = FALSE)
#   KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION <- data.frame(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION[-1,])
#   colnames(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION) <- c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")
#
#   KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse <- HSsymbol2MMsymbol(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION,"KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")
#   KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2 <- unique(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2),
#                                   ctrl = 5, name = 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION')
#   FeaturePlot(object = scRNA.SeuObj, features = 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots7 <- VlnPlot(scRNA.SeuObj, features = c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1"), split.by = "sample", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots7 = plots7, ncol = 1)
#
#   plots7 <- VlnPlot(scRNA.SeuObj, features = c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1"), split.by = "Cachexia", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p7 <- wrap_plots(plots7 = plots7, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                                   XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                                   LegTextSize = 15)  +
#     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p7
#
#
# ##### REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX #####
#   REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX.txt"),header=T, stringsAsFactors = FALSE)
#   REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX <- data.frame(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX[-1,])
#   colnames(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX) <- c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")
#
#   REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse <- HSsymbol2MMsymbol(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX,"REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")
#   REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2 <- unique(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse[,2])
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2[1:24], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features = REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2[25:50], min.cutoff = "q9")
#   # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")
#
#   scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2),
#                                   ctrl = 5, name = 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX')
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1', split.by = "Cachexia")+
#     scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
#                            guide = "colourbar",midpoint = 0, labs(fill ="Exp"))
#
#   FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1', split.by = "sample",
#               cols = c("#0077b6", "white",  "#ef476f"))
#   plots8 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1"), split.by = "sample", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
#   wrap_plots(plots8 = plots8, ncol = 1)
#
#   plots8 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1"), split.by = "Cachexia", group.by = "celltype",
#                     pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
#   p8 <- wrap_plots(plots8 = plots8, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
#                                                                   XtextSize=18,  YtextSize=,18, xangle = 90,
#                                                                   LegTextSize = 15)  +
#     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
#
#   p8
#
#
# ##### Export PDF file #####
#   pdf(file = paste0(Save.Path,"/PBMC_GSEA_Violin.pdf"),
#       width = 7, height = 7 )
#     p1 %>% print()
#     p2 %>% print()
#     p3 %>% print()
#     p4 %>% print()
#     p5 %>% print()
#     p6 %>% print()
#     p7 %>% print()
#     p8 %>% print()
#   dev.off() # graphics.off()
#
#
# ##### Clean up the obj #####
#   rm(p1, p2, p3, p4, p5, p6, p7, p8)
#   rm(HALLMARK_TNFA_SIGNALING_VIA_NFKB,
#      HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse,
#      HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2,
#      REACTOME_ANTIMICROBIAL_PEPTIDES,
#      REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse,
#      REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2,
#      REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION,
#      REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse,
#      REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2,
#      REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION,
#      REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse,
#      REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2,
#      REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2,
#      REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse,
#      REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2,
#      REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS,
#      REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse,
#      REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2,
#      KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION,
#      KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse,
#      KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2,
#      REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX,
#      REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse,
#      REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2
#      )
#
#
