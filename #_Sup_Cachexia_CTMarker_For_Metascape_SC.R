##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load libray #####
library(tidyverse)
library(stringr)
library(ggplot2)
library(Seurat)

# ##### Load RData  #####
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_SC_Main/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData")


##### Current path and new folder setting #####
Version = paste0(Sys.Date(),"_","SC_CTMarker")
Save.Path = paste0(getwd(),"/",Version)
dir.create(Save.Path)
SampleType="SC"

PathName <- Save.Path

##### Cluster marker gene for Metascape #####
  Idents(SC.combined) <- "celltype"
  set.seed(1) # Fix the seed
  SC_markers_CT <- FindAllMarkers(SC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


  # Filter the top markers and plot the heatmap
  top_NSet = 25
  SC_markers_CT %>%
    group_by(cluster) %>%
    top_n(n = top_NSet, wt = avg_log2FC) -> top_N
  SC.combined <- ScaleData(SC.combined, verbose = FALSE)
  DoHeatmap(SC.combined, features = top_N$gene) + NoLegend()
  write.table(top_N, file=paste0(PathName,"/SC_CT_top",top_NSet,"Gene.txt"),sep="\t", row.names=T
              , quote = FALSE)

  Cluster_Marker.lst <- list()
  Cluster.list <- as.character(unique(top_N$cluster))

  for (i in 1:length(Cluster.list)) {
    top_N_Type <- top_N[top_N$cluster == Cluster.list[i],]
    Cluster_Marker.lst[[Cluster.list[i]]] <- top_N_Type[['gene']]

  }

  ## list to df
  Cluster_Marker.df <-""
  for(i in c(1:length(Cluster.list))){
    if(Cluster_Marker.df == ""){
      Cluster_Marker.df <- as.data.frame(Cluster_Marker.lst[[paste0(CellType.list[i])]])
    }
    else{
    Cluster_Marker.df <- cbind(Cluster_Marker.df,as.data.frame(Cluster_Marker.lst[[paste0(CellType.list[i])]]))
  }}
  colnames(Cluster_Marker.df) <- Cluster.list
  #Cluster_Marker.df <- as.data.frame(t(Cluster_Marker.df))
  write.table(Cluster_Marker.df,
              file=paste0(PathName,"/SC_CT_top",
                          top_NSet,"Gene_forMetascape.csv"),sep=",",
              row.names=F,col.names = T, quote = FALSE)
