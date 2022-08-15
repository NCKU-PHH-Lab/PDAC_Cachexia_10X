
##### Cluster marker gene for Metascape #####
  Idents(PBMC.combined) <- "celltype"
  set.seed(1) # Fix the seed
  PBMC_markers_CT <- FindAllMarkers(PBMC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Filter the top markers and plot the heatmap
  top_NSet = 25
  PBMC_markers_CT %>%
    group_by(cluster) %>%
    top_n(n = top_NSet, wt = avg_log2FC) -> top_N
  PBMC.combined <- ScaleData(PBMC.combined, verbose = FALSE)
  DoHeatmap(PBMC.combined, features = top_N$gene) + NoLegend()
  write.table(top_N, file=paste0(PathName,"/",RVersion,"/PBMC_CT_top",top_NSet,"Gene.txt"),sep="\t", row.names=T
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
              file=paste0(PathName,"/",RVersion,"/PBMC_CT_top",
                          top_NSet,"Gene_forMetascape.csv"),sep=",", 
              row.names=F,col.names = T, quote = FALSE)
  