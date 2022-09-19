##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

source("FUN_netVisual_heatmap_Ch.R")

#### PBMC CC ####
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-16_PBMC_CellChat_Multi_CC/010_Cell_Cell_Interaction_Multi.RData")
cellchat_PBMC_CC <- cellchat
rm(cellchat)

load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-16_PBMC_CellChat_Multi_ECM/010_Cell_Cell_Interaction_Multi.RData")
cellchat_PBMC_ECM <- cellchat
rm(cellchat)

load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-16_PBMC_CellChat_Multi_Secret/010_Cell_Cell_Interaction_Multi.RData")
cellchat_PBMC_Secret <- cellchat
rm(cellchat)

#### SC CC ####
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-16_SC_CellChat_Multi_CC/010_Cell_Cell_Interaction_Multi.RData")
cellchat_SC_CC <- cellchat
rm(cellchat)

load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-16_SC_CellChat_Multi_ECM/010_Cell_Cell_Interaction_Multi.RData")
cellchat_SC_ECM <- cellchat
rm(cellchat)

load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-16_SC_CellChat_Multi_Secret/010_Cell_Cell_Interaction_Multi.RData")
cellchat_SC_Secret <- cellchat
rm(cellchat)


rm(list=setdiff(ls(), c("cellchat_PBMC_CC","cellchat_PBMC_ECM","cellchat_PBMC_Secret","netVisual_heatmap_Ch")))

##### Heatmap #####
#### Ori ####
ComHeatmap1 <- netVisual_heatmap(cellchat_PBMC_CC)
ComHeatmap1

ComHeatmap2 <- netVisual_heatmap(cellchat_PBMC_CC, measure = "weight",
                                 #color.heatmap = c("#2166ac", "#b2182b"),
)
ComHeatmap2

#> Do heatmap based on a merged object
ComHeatmap1 + ComHeatmap2

#### Charlene Modify ####
## PBMC_CC
ComHeatmapTest<- netVisual_heatmap(cellchat_PBMC_CC, measure = "weight")
ComHeatmapTest # -0.4, 0.2

ComHeatmap2_PBMC_CC_Ch <- netVisual_heatmap_Ch(cellchat_PBMC_CC, measure = "weight",
                                    #color.heatmap = c("#2166ac", "#b2182b"),
                                    MinSet = -0.4, MaxSet = 0.2)
ComHeatmap2_PBMC_CC_Ch

## PBMC_ECM
ComHeatmapTest <- netVisual_heatmap(cellchat_PBMC_ECM, measure = "weight")
ComHeatmapTest # -0.4, 0.4

ComHeatmap2_PBMC_ECM_Ch <- netVisual_heatmap_Ch(cellchat_PBMC_ECM, measure = "weight",
                                               #color.heatmap = c("#2166ac", "#b2182b"),
                                               MinSet = -0.4, MaxSet = 0.4)
ComHeatmap2_PBMC_ECM_Ch

## PBMC_Secret
ComHeatmapTest <- netVisual_heatmap(cellchat_PBMC_Secret, measure = "weight")
ComHeatmapTest # -0.2, 0.4

ComHeatmap2_PBMC_Secret_Ch <- netVisual_heatmap_Ch(cellchat_PBMC_Secret, measure = "weight",
                                                #color.heatmap = c("#2166ac", "#b2182b"),
                                                MinSet = -0.2, MaxSet = 0.4)
ComHeatmap2_PBMC_Secret_Ch




##### Current path and new folder setting*  #####
ProjectName = "PBMC_CellChat_Beautify_Plot" # Secret, ECM, CC
Version = paste0(Sys.Date(),"_",ProjectName)
Save.Path = paste0(getwd(),"/",Version)
## Create new folder
if (!dir.exists(Save.Path)){
  dir.create(Save.Path)
}

## Export PDF
pdf(file = paste0(Save.Path,"/",ProjectName,"_Compared_Heatmap.pdf"),
    width = 12,  height = 7
)

ComHeatmap2_PBMC_CC_Ch

dev.off()
