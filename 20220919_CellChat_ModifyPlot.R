source("FUN_netVisual_heatmap_Ch.R")

#### PBMC CC ####
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-16_PBMC_CellChat_Multi_CC/010_Cell_Cell_Interaction_Multi.RData")

cellchat_CC <- cellchat


### Heatmap ###
## Ori
ComHeatmap1 <- netVisual_heatmap(cellchat_CC)
ComHeatmap1

ComHeatmap2 <- netVisual_heatmap(cellchat_CC, measure = "weight",
                                 #color.heatmap = c("#2166ac", "#b2182b"),
)
ComHeatmap2

#> Do heatmap based on a merged object
ComHeatmap1 + ComHeatmap2

## Charlene Modify
ComHeatmap2_Ch <- netVisual_heatmap_Ch(cellchat_CC, measure = "weight",
                                    #color.heatmap = c("#2166ac", "#b2182b"),
                                    MinSet = -0.4, MaxSet = 0.2)
ComHeatmap2_Ch




