## Ref: https://github.com/sqjin/CellChat
## Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

# SampleType = "PBMC"
#
# ## INTCHG: Interchangeable
# ## SubType Setting
#   if(SampleType == "PBMC"){
#     # For PBMC
#     scRNA.SeuObj <- PBMC.combined
#
#   }else if(SampleType == "SC"){
#     # For SC
#     scRNA.SeuObj <- SC.combined
#
#   }


# ##### Presetting ######
#   rm(list = ls()) # Clean variable
#   memory.limit(150000)

#### Installation and load the required libraries ####
  #### Basic installation ####
  ## Package.set
  Package.set <- c("tidyverse","CellChat","patchwork","NMF","ggalluvial","Seurat")
  ## Check whether the installation of those packages is required
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### BiocManager installation ####
  ## Package.set
  Package.set <- c("ComplexHeatmap")
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

##### Load CellChat object of each dataset and then merge together #####
  source("FUN_CellChatOne.R")

  # load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_PBMC_Main/09_4_GSEA_Analysis_(SSA).RData")
  # load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_Results_1stSubmission/2022-09-09_PBMC_Main//09_4_GSEA_Analysis_(SSA).RData")

  scRNA_EOCX.combined <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cachexia"]] %in% "EOCX"]
  scRNA_PreCX.combined <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["Cachexia"]] %in% "PreCX"]

  ## SubType Setting
  if(SampleType == "PBMC"){
    ## For PBMC
    Cell_Type_Order.set <- c("Mac1", "Mac2", "Mac3","Neu", "T", "CD4+T", "CD8+T",
                             "NK", "B" , "Mast",  "Ery")
  }else{
    ## For SC
    Cell_Type_Order.set <- c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" ,
                             "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
                             "Fib1", "Fib2", "Fib3")
  }


  scRNA_EOCX.combined$celltype <- factor(scRNA_EOCX.combined$celltype,
                                  levels = Cell_Type_Order.set)
  scRNA_PreCX.combined$celltype <- factor(scRNA_PreCX.combined$celltype,
                                      levels = Cell_Type_Order.set)


  ## ECM-Receptor
  CellChatOne(scRNA_EOCX.combined,
              signalingtype = "ECM-Receptor", projectName = "ECM_EOCX",
              save.path = paste0(Save.Path,"/",SampleType,"_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) ->   CellChat_ECM_EOCX.lt

  CellChatOne(scRNA_PreCX.combined,
              signalingtype = "ECM-Receptor", projectName = "ECM_PreCX",
              save.path = paste0(Save.Path,"/",SampleType,"_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) ->   CellChat_ECM_PreCX.lt

  ## Cell-Cell Contact
  CellChatOne(scRNA_EOCX.combined,
              signalingtype = "Cell-Cell Contact", projectName = "CC_EOCX",
              save.path = paste0(Save.Path,"/",SampleType,"_CellCell_Interaction"),
              groupby = "celltype",species =  "Mouse"
  ) -> CellChat_CC_EOCX.lt

  CellChatOne(scRNA_PreCX.combined,
              signalingtype = "Cell-Cell Contact", projectName = "CC_PreCX",
              save.path = paste0(Save.Path,"/",SampleType,"_CellCell_Interaction"),
              groupby = "celltype",species =  "Mouse"
  ) -> CellChat_CC_PreCX.lt


  ## Secreted Signaling
  CellChatOne(scRNA_EOCX.combined,
              signalingtype = "Secreted Signaling", projectName = "Secret_EOCX",
              save.path = paste0(Save.Path,"/",SampleType,"_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) -> CellChat_Secret_EOCX.lt

  CellChatOne(scRNA_PreCX.combined,
              signalingtype = "Secreted Signaling", projectName = "Secret_PreCX",
              save.path = paste0(Save.Path,"/",SampleType,"_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) -> CellChat_Secret_PreCX.lt

  ##### save.image #####
  save.image(paste0(Save.Path,"/010_Cell_Cell_Interaction.RData"))

##***************************************************************************##
