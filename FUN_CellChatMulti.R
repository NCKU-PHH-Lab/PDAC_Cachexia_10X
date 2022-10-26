## Ref: https://github.com/sqjin/CellChat
## Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

CellChatMulti <- function(seuratObject,
                        signalingtype = "ECM-Receptor", projectName = "ECM",
                        save.path = paste0(Save.Path,"/B04_CellCell_Interaction"),
                        groupby = "celltype", species = "Human" # species = c("Human","Mouse")
                        ){

# #### Load the required libraries ####
#   ## Check whether the installation of those packages is required from basic
#   Package.set <- c("tidyverse","CellChat","patchwork","reticulate","anndata","Seurat","NMF","ggalluvial")
#   for (i in 1:length(Package.set)) {
#     if (!requireNamespace(Package.set[i], quietly = TRUE)){
#       install.packages(Package.set[i])
#     }
#   }
#   ## Load Packages
#   lapply(Package.set, library, character.only = TRUE)
#   rm(Package.set,i)
#
#   ## Check whether the installation of those packages is required from BiocManager
#   if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   Package.set <- c("basilisk","zellkonverter","SeuratDisk")
#   for (i in 1:length(Package.set)) {
#     if (!requireNamespace(Package.set[i], quietly = TRUE)){
#       BiocManager::install(Package.set[i])
#     }
#   }
#   ## Load Packages
#   lapply(Package.set, library, character.only = TRUE)
#   rm(Package.set,i)
#
#   options(stringsAsFactors = FALSE)
#
# ##### Current path and new folder setting*  #####
#   SignalingType = signalingtype # Secreted Signaling, ECM-Receptor, Cell-Cell Contact
#   ProjectName = projectName # Secret, ECM, CC
#   Save_Path = save.path
#
#   ## Create new folder
#   if (!dir.exists(Save_Path)){
#     dir.create(Save_Path)
#   }
#
#
#   #### Identify and visualize incoming communication pattern of target cells ####
#
#   pdf(file = paste0(PathLR,"/",ProjectName,"_LRPair_GlobalPatterns_incoming.pdf"),
#       width = 12,  height = 8
#   )
#
#   try({
#   nPatterns = 4
#   cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
#   })
#   # river plot
#   try({
#   netAnalysis_river(cellchat, pattern = "incoming")
#   })
#   #> Please make sure you have load `library(ggalluvial)` when running this function
#
#   # dot plot
#   try({
#   netAnalysis_dot(cellchat, pattern = "incoming")
#
#   P.incoming
#   })
#   #graphics.off()
#   dev.off()
#
#
# ##### Part V: Save the CellChat object #####
#   saveRDS(cellchat, file = paste0(Save_Path,"/", ProjectName,"_CellChat.rds"))
#
#
# ##### Export Result #####
#     CellChat.lt <- list(CellChatObj = cellchat,
#                         DataBase_Use = CellChatDB.use,
#                         DB_Interact_Sig = DB_Interact_Sig.df
#     )
#
#     return(CellChat.lt)
}
