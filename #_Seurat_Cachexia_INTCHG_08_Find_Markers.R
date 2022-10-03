# ## INTCHG: Interchangeable
#   scRNA.SeuObj <- PBMC.combined
#   # scRNA.SeuObj <- SC.combined

##### 08_1 Find CCmarker in different Cell type and VolcanoPlot (SSA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")
  scRNA.SeuObj$celltype <- Idents(scRNA.SeuObj)
  scRNA.SeuObj$celltype.Cachexia <- paste(Idents(scRNA.SeuObj), scRNA.SeuObj$Cachexia, sep = "_")
  scRNA.SeuObj$celltype.Cachexia.gender <- paste(Idents(scRNA.SeuObj), scRNA.SeuObj$Cachexia, scRNA.SeuObj$Sex, sep = "_")
  Idents(scRNA.SeuObj) <- "celltype.Cachexia.gender"

  scRNA.SeuObj$Cachexia.gender <- paste(scRNA.SeuObj$Cachexia, scRNA.SeuObj$Sex, sep = "_")

  DefaultAssay(scRNA.SeuObj) <- "RNA"

  ####-------------- Find Marker gene in Male --------------####
  CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[["celltype"]]))
  # CellType.list <- CellType.list[-9]

  dir.create(paste0(Save.Path,"/",SampleType,"_SSA_Male_FindMarkers"))

  # About 15 mins
  CCMarker_Male.lt <- list()
  for(i in c(1:length(CellType.list))){
    try({
      CCMarker_Male.lt[[i]] <- Find_Markers(scRNA.SeuObj,
                                            paste0(CellType.list[i],"_EO_Male"),
                                            paste0(CellType.list[i],"_LO_Male"),
                                            CellType.list[i],
                                            Path = Save.Path,
                                            ResultFolder = "",SampleType,"_SSA_Male_FindMarkers")
      # names(CCMarker_Male.lt)[[i]] <- paste0("CCMarker_Male.lt.",CellType.list[i])
      names(CCMarker_Male.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)

  CCMarker_Male.lt <- CCMarker_Male.lt[!unlist(lapply(CCMarker_Male.lt,is.null))]


  ## Generate pdf and tif file for Male VolcanoPlot
  dir.create(paste0(Save.Path,"/",SampleType,"_SSA_Male_VolcanoPlot/"))

  pdf(file = paste0(Save.Path,"/",SampleType,"_SSA_Male_VolcanoPlot/",SampleType,"_SSA_Male_VolcanoPlot.pdf"),
      width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_Male.lt[[i]][["CCMarker.S"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+
              ggtitle(paste0(SampleType,"_Male_",CellType.list[i]))
      )
    })
  }
  dev.off() # graphics.off()
  rm(i)

  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/",SampleType,"_SSA_Male_VolcanoPlot/",CellType.list[i],".tif"),
           width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_Male.lt[[i]][["CCMarker.S"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Male.lt[[i]][["CCMarker.S_Neg_List"]])+
              ggtitle(paste0(SampleType,"_Male_",CellType.list[i]))
      )

      graphics.off()
    })
  }
  rm(i)

  ####-------------- Find Marker gene in Female --------------####
  CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[["celltype"]]))
  # CellType.list <- CellType.list[-9] # Some cluster with cell lower than 3

  dir.create(paste0(Save.Path,"/",SampleType,"_SSA_Female_FindMarkers"))

  # About 15 mins
  CCMarker_Female.lt <- list()
  for(i in c(1:length(CellType.list))){
    try({
      CCMarker_Female.lt[[i]] <- Find_Markers(scRNA.SeuObj,
                                              paste0(CellType.list[i],"_EO_Female"),
                                              paste0(CellType.list[i],"_LO_Female"),
                                              CellType.list[i],
                                              Path = Save.Path,
                                              ResultFolder = paste0(SampleType,"_SSA_Female_FindMarkers"))
      # names(CCMarker_Female.lt)[[i]] <- paste0("CCMarker_Female.lt.",CellType.list[i])
      names(CCMarker_Female.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)

  CCMarker_Female.lt <- CCMarker_Female.lt[!unlist(lapply(CCMarker_Female.lt,is.null))]

  ## Generate pdf and tif file for Female VolcanoPlot
  dir.create(paste0(Save.Path,"/",SampleType,"_SSA_Female_VolcanoPlot/"))

  pdf(file = paste0(Save.Path,"/",SampleType,"_SSA_Female_VolcanoPlot/",SampleType,"_SSA_Female_VolcanoPlot.pdf"),
      width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_Female.lt[[i]][["CCMarker.S"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+
              ggtitle(paste0(SampleType,"_Female_",CellType.list[i]))
      )
    })
  }
  dev.off() # graphics.off()
  rm(i)

  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/",SampleType,"_SSA_Female_VolcanoPlot/",CellType.list[i],".tif"),
           width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_Female.lt[[i]][["CCMarker.S"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_Female.lt[[i]][["CCMarker.S_Neg_List"]])+
              ggtitle(paste0(SampleType,"_Female_",CellType.list[i]))
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
  # pdf(file = paste0(Save.Path,"/",SampleType,"_Female_VolcanoPlot.pdf"),width = 7, height = 7 )

  dir.create(paste0(Save.Path,"/",SampleType,"_SSA_VennDiagrame"))
  Venn_CCMarker_Pos <- list()
  for(i in c(1:length(CellType.list))){
    try({
      Venn_CCMarker_Pos[[i]] <- Venn_Intersect(CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.S_Pos_List"]],
                                               CCMarker_Female.lt[[paste0(CellType.list[i])]][["CCMarker.S_Pos_List"]],
                                               CellType.list[i],"Pos","#9d0208","#f08080", SampleType=SampleType,
                                               PathName = paste0(Save.Path,"/",SampleType,"_SSA_VennDiagrame"))
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
                                               CellType.list[i],"Neg","#00296b","#1368aa", SampleType=SampleType,
                                               PathName = paste0(Save.Path,"/",SampleType,"_SSA_VennDiagrame"))

      names(Venn_CCMarker_Neg)[[i]] <- paste0("Venn_CCMarker.",CellType.list[i],"_Neg")
    })
  }
  rm(i)

  #### Save RData ####
    save.image(paste0(Save.Path,"/08_1_Find_CCmarker_in_different_Cell_type_and_VennDiagrame(SSA_IntersectCT).RData"))


##### 08_2 Find CCmarker in different Cell type and VolcanoPlot (SPA) ########
  ### Define group by different phenotype ###
  source("FUN_Find_Markers.R")

  Idents(scRNA.SeuObj) <- "celltype.Cachexia"
  #CellType.list <- as.character(unique(scRNA.SeuObj@meta.data[["celltype"]]))
  dir.create(paste0(Save.Path,"/",SampleType,"_SPA_FindMarkers"))

  CCMarker_SPA.lt <- list()
  for(i in c(1:length(CellType.list))){
    try({
      CCMarker_SPA.lt[[i]] <- Find_Markers(scRNA.SeuObj,
                                           paste0(CellType.list[i],"_EO"),
                                           paste0(CellType.list[i],"_LO"),
                                           CellType.list[i],
                                           Path = Save.Path,
                                           ResultFolder = paste0(SampleType,"_SPA_FindMarkers"))

      # names(CCMarker_SPA.lt)[[i]] <- paste0("CCMarker_SPA.lt.",CellType.list[i])
      names(CCMarker_SPA.lt)[[i]] <- paste0(CellType.list[i])
    })
  }
  rm(i)

  CCMarker_SPA.lt <- CCMarker_SPA.lt[!unlist(lapply(CCMarker_SPA.lt,is.null))]


  ## Generate pdf and tif file for VolcanoPlot
  dir.create(paste0(Save.Path,"/",SampleType,"_SPA_VolcanoPlot/"))

  pdf(file = paste0(Save.Path,"/",SampleType,"_SPA_VolcanoPlot/",SampleType,"_SPA_VolcanoPlot.pdf"),width = 7, height = 7 )
  for (i in 1:length(CellType.list)) {
    try({
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]], ShowGeneNum = 6)+
              ggtitle(paste0(SampleType,"_",CellType.list[i]))
      )
    })
  }
  # graphics.off()
  dev.off()
  rm(i)

  for (i in 1:length(CellType.list)) {
    try({
      tiff(file = paste0(Save.Path,"/",SampleType,"_SPA_VolcanoPlot/",SampleType,"_SPA_VolcanoPlot",CellType.list[i],".tif"), width = 17, height = 17, units = "cm", res = 200)
      print(VolcanoPlot(CCMarker_SPA.lt[[i]][["CCMarker.S"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Pos_List"]],
                        CCMarker_SPA.lt[[i]][["CCMarker.S_Neg_List"]])+ ggtitle(paste0(SampleType,"_",CellType.list[i]))
      )

      graphics.off()
    })
  }
  rm(i)

  #### Save RData ####
    save.image(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))


#####------------------------------------------------------------------------------------------------------------#####

