### Ref: Add P-values and Significance Levels to ggplots
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

## https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/plotting-models.html
## https://github.com/kassambara/ggpubr/issues/111

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

Save.Path <- c("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_SC_Main")
SampleType = "SC"

#### Installation and load the required libraries ####
#### Basic installation ####
## Package.set
Package.set <- c("tidyverse","CellChat","patchwork","NMF","ggalluvial","Seurat","ggpubr", "stringr", "Hmisc",
                 "circlize","ComplexHeatmap")
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

##### Function setting  #####
## Call function
source("FUN_CellChatOne.R")

##### Condition Setting ####
# ## INTCHG: Interchangeable
# SampleTypeSet = "PBMC"

## CellChat DB Set
CCDBType = "Secret" # c("ECM","CC","Secret")

# ##### Load Data* #####
# ## Load RData
# # load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_Results_1stSubmission/2022-09-09_PBMC_Main/06_Cell_type_annotation.RData")
# load(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))
#
# ## INTCHG: Interchangeable
# ## SubType Setting
# if(SampleType == "PBMC"){
#   ## For PBMC
#   scRNA.SeuObj <- PBMC.combined
#
#   # Order the cell type
#   CellType.Order = c("Mac1", "Mac2","Mac3","Neu","T","CD4+T","CD8+T","NK","B","Mast","Ery")
#   scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
#                                                  levels = CellType.Order)
#
#
# }else if(SampleType == "SC"){
#   ## For SC
#   scRNA.SeuObj <- SC.combined
#
#   # Order the cell type
#   CellType.Order = c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" , "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
#                      "Fib1", "Fib2", "Fib3")
#   scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
#                                                  levels = CellType.Order)
#
# }
#
#   # Clean up
#   rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Save.Path","CCDBType","CellType.Order",
#                           "CCMarker_Female.lt","CCMarker_Male.lt","CCMarker_SPA.lt")))

## Load CellChat rds
cellchat.EOCX <- readRDS(paste0(Save.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_EOCX_CellChat.rds"))
cellchat.PreCX <- readRDS(paste0(Save.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_PreCX_CellChat.rds"))

object.list <- list(PreCX = cellchat.PreCX, EOCX = cellchat.EOCX)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rm(object.list, cellchat.EOCX, cellchat.PreCX)


##### Current path and new folder setting  #####
Version = paste0(Sys.Date(),"_", SampleType, "_", CCDBType, "_CellChat_Bubble")
SaveCC.Path = paste0(Save.Path,"/",Version)
dir.create(SaveCC.Path)


##### Bubble Plot #####

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45, thresh = 0.01)

netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, thresh = 0.01)


sizeRange <- c(0,5)

# importing the ggplot2 library
library(ggplot2)


##### Export PDF #####

## Duc
pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Duc.pdf"),width = 30, height = 10 )
netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
  #scale_size(range = sizeRange, name="P-Value")+
  ggtitle(paste0("Allograft_Duc : ",CCDBType))+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))+
  theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0))

dev.off()

## Mac
pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Mac.pdf"),width = 20, height = 8 )
netVisual_bubble(cellchat, sources.use = c(7,8,9,10,11) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
  #scale_size(range = sizeRange, name="P-Value")+
  ggtitle(paste0("Allograft_Mac : ",CCDBType))+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))+
  theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0))

dev.off()

## Fib
pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Fib.pdf"),width = 20, height = 8 )
netVisual_bubble(cellchat, sources.use = c(12,13,14) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
  #scale_size(range = sizeRange, name="P-Value")+
  ggtitle(paste0("Allograft_Fib : ",CCDBType))+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))+
  theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0))

dev.off()

if(CCDBType == "ECM")
{
  ## Significant For ECM
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig.pdf"),width = 25, height = 20 )
  netVisual_bubble(cellchat, sources.use = c(1,2,3,4,7,8,9,12) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0)) -> p
  p %>% print()

  dev.off()

  ## Significant For ECM Duc
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig_Duc.pdf"),width = 20, height = 15 )
  netVisual_bubble(cellchat, sources.use = c(1,2,3,4) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft_Duc : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0))-> p
  p %>% print()

  dev.off()

  ## Significant For ECM Mac
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig_Mac.pdf"),width = 20, height = 10 )
  netVisual_bubble(cellchat, sources.use = c(7,8,9) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft_Mac : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=10, angle=0)) -> p
  p %>% print()

  dev.off()

  ## Significant For ECM Fib
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig_Fib.pdf"),width = 10, height = 13 )
  netVisual_bubble(cellchat, sources.use = c(12) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft_Fib : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=10, angle=0)) -> p
  p %>% print()

  dev.off()

  ## Significant For ECM Mac Fib
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig_Mac_Fib.pdf"),width = 20, height = 15 )
  netVisual_bubble(cellchat, sources.use = c(7,8,9,12) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft_Mac_Fib : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=10, angle=0)) -> p
  p %>% print()

  dev.off()
}

if(CCDBType == "CC")
{
  ## Significant For CC
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig.pdf"),width = 20, height = 4 )
  netVisual_bubble(cellchat, sources.use = c(1,2,3,4,7) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0)) -> p
  p %>% print()

  dev.off()
}

if(CCDBType == "Secret")
{
  ## Significant For Secret
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig.pdf"),width = 35, height = 15 )
  netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,7,8,9,12) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0)) -> p
  p %>% print()

  dev.off()

  ## Significant For Secret Duc
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig_Duc.pdf"),width = 30, height = 10 )
  netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft_Duc : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=12, angle=0)) -> p
  p %>% print()

  dev.off()

  ## Significant For ECM Mac Fib
  pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot_Sig_Mac_Fib.pdf"),width = 20, height = 10 )
  netVisual_bubble(cellchat, sources.use = c(7,8,9,12) , comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
    #scale_size(range = sizeRange, name="P-Value")+
    ggtitle(paste0("Allograft_Mac_Fib : ",CCDBType))+
    theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
    theme(axis.text.y = element_text(face="bold", color="black",size=10, angle=0)) -> p
  p %>% print()

  dev.off()
}

rm(p)

pdf(file = paste0(SaveCC.Path,"/",Version,"_LR_BubblePlot.pdf"),width = 40, height = 25 )
netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, thresh = 0.01)+
  #scale_size(range = sizeRange, name="P-Value")+
  ggtitle(paste0("Allograft : ",CCDBType))+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))+
  theme(axis.text.y = element_text(face="bold", color="black",size=9, angle=0))

dev.off()



