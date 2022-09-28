## Ref: https://github.com/sqjin/CellChat
## Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

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

  load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_SC_Main/09_4_GSEA_Analysis_(SSA).RData")
  SC_EO.combined <- SC.combined[,SC.combined@meta.data[["Cachexia"]] %in% "EO"]
  SC_LO.combined <- SC.combined[,SC.combined@meta.data[["Cachexia"]] %in% "LO"]

  Cell_Type_Order.set <- c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" , "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
                           "Fib1", "Fib2", "Fib3")
  SC_EO.combined$celltype <- factor(SC_EO.combined$celltype,
                                      levels = Cell_Type_Order.set)
  SC_LO.combined$celltype <- factor(SC_LO.combined$celltype,
                                      levels = Cell_Type_Order.set)


  ## ECM-Receptor
  CellChatOne(SC_EO.combined,
              signalingtype = "ECM-Receptor", projectName = "ECM_EO",
              save.path = paste0(Save.Path,"/SC_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) ->   CellChat_ECM_EO.lt

  CellChatOne(SC_LO.combined,
              signalingtype = "ECM-Receptor", projectName = "ECM_LO",
              save.path = paste0(Save.Path,"/SC_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) ->   CellChat_ECM_LO.lt

  ## Cell-Cell Contact
  CellChatOne(SC_EO.combined,
              signalingtype = "Cell-Cell Contact", projectName = "CC_EO",
              save.path = paste0(Save.Path,"/SC_CellCell_Interaction"),
              groupby = "celltype",species =  "Mouse"
  ) -> CellChat_CC_EO.lt

  CellChatOne(SC_LO.combined,
              signalingtype = "Cell-Cell Contact", projectName = "CC_LO",
              save.path = paste0(Save.Path,"/SC_CellCell_Interaction"),
              groupby = "celltype",species =  "Mouse"
  ) -> CellChat_CC_LO.lt


  ## Secreted Signaling
  CellChatOne(SC_EO.combined,
              signalingtype = "Secreted Signaling", projectName = "Secret_EO",
              save.path = paste0(Save.Path,"/SC_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) -> CellChat_Secret_EO.lt

  CellChatOne(SC_LO.combined,
              signalingtype = "Secreted Signaling", projectName = "Secret_LO",
              save.path = paste0(Save.Path,"/SC_CellCell_Interaction"),
              groupby = "celltype",species = "Mouse"
  ) -> CellChat_Secret_LO.lt

  # ##### save.image #####
  # save.image(paste0(Save.Path,"/010_Cell_Cell_Interaction.RData"))

##***************************************************************************##

##### Merge cellchat #####

  ##### Load rds #####
  cellchat.EO <- readRDS("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_SC_Main/SC_CellCell_Interaction/ECM_EO_CellChat.rds")
  cellchat.LO <- readRDS("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_SC_Main/SC_CellCell_Interaction/ECM_LO_CellChat.rds")

  object.list <- list(LO = cellchat.LO, EO = cellchat.EO)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))

  cellchat

  ##### Current path and new folder setting*  #####
  ProjectName = "SC_CellChat_Multi_ECM" # Secret, ECM, CC
  Version = paste0(Sys.Date(),"_",ProjectName)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }


#### Part I: Predict general principles of cell-cell communication #####
  CellType.set <- cellchat@meta[["celltype"]] %>% unique()
  Pathway.set <- c(cellchat@netP[["EO"]][["pathways"]],cellchat@netP[["LO"]][["pathways"]]) %>% unique()

  pdf(file = paste0(Save.Path,"/",ProjectName,"_01_Predict_general_principles.pdf"),
      width = 12,  height = 7
  )

  # Compare the total number of interactions and interaction strength
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2


## Compare the number of interactions and interaction strength among different cell populations
  # Differential number of interactions or interaction strength among different cell populations
  ## circle plot
  par(mfrow = c(1,2), xpd=TRUE)
  gg3 <- netVisual_diffInteraction(cellchat, weight.scale = T)
  gg4 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

  ## Heatmap
  gg5 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg6 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  gg5 + gg6


  weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  }

  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Weight of interactions - ", names(object.list)[i]))
  }

# ## Differential number of interactions or interaction strength among different cell types
#   group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
#   group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
#   object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
#   cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#
#   weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
#   par(mfrow = c(1,2), xpd=TRUE)
#   for (i in 1:length(object.list)) {
#     netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
#   }
#
#   par(mfrow = c(1,2), xpd=TRUE)
#   gg7 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
#   gg8 <-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

  ## Compare the major sources and targets in 2D space
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  }
  patchwork::wrap_plots(plots = gg)

  # gg9 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD14_Mono") # signaling.exclude = "MIF"
  # gg10 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "DC")
  # patchwork::wrap_plots(plots = list(gg9,gg10))

  gg9 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = CellType.set[1]) # signaling.exclude = "MIF"
  gg10 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = CellType.set[2])
  patchwork::wrap_plots(plots = list(gg9,gg10))

  dev.off()


  pdf(file = paste0(Save.Path,"/",ProjectName,"_01_Predict_general_principles_Compare_Sourcestargets_in2D.pdf"),
      width = 12,  height = 7
  )
    ## Compare the major sources and targets in 2D space
    num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
    weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
    gg <- list()
    for (i in 1:length(object.list)) {
      gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
    }
    patchwork::wrap_plots(plots = gg)

    for (i in 1:(length(CellType.set)/2)) {
      gg9 <- netAnalysis_signalingChanges_scatter(cellchat, CellType.set[i]) # signaling.exclude = "MIF"
      gg10 <- netAnalysis_signalingChanges_scatter(cellchat, CellType.set[i*2])
      patchwork::wrap_plots(plots = list(gg9,gg10)) %>% print()
    }
    rm(i)
  dev.off()


##### Part II: Identify the conserved and context-specific signaling pathways #####
  pdf(file = paste0(Save.Path,"/",ProjectName,"_02_Identify_the_conserved_and_context-specific_SigPath.pdf"),
      width = 7,  height = 7
  )

    ### Identify signaling groups based on their functional similarity
    #> Compute signaling network similarity for datasets 1 2
    cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
    #> Manifold learning of the signaling networks for datasets 1 2
    cellchat <- netEmbedding(cellchat, type = "functional")
    #> Classification learning of the signaling networks for datasets 1 2
    cellchat <- netClustering(cellchat, type = "functional")

    #> 2D visualization of signaling networks from datasets 1 2
    try({
      gg2_1 <- netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
      gg2_1 %>% print()
    })


    # ZoomIn
    try({
    gg2_2 <- netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
    gg2_2 %>% print()
    })

    ## Compute and visualize the pathway distance in the learned joint manifold
    #> Compute the distance of signaling networks between datasets 1 2
    gg2_3 <- rankSimilarity(cellchat, type = "functional")
    gg2_3 %>% print()

    ## Compute and visualize the pathway distance in the learned joint manifold
    #> Compute the distance of signaling networks between datasets 1 2
    gg2_4 <- rankSimilarity(cellchat, type = "functional")
    gg2_4 %>% print()


    ### Identify signaling groups based on structure similarity
    #> Compute signaling network similarity for datasets 1 2
    cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
    #> Manifold learning of the signaling networks for datasets 1 2
    cellchat <- netEmbedding(cellchat, type = "structural")
    #> Classification learning of the signaling networks for datasets 1 2
    cellchat <- netClustering(cellchat, type = "structural")
    #> 2D visualization of signaling networks from datasets 1 2
    gg2_5 <- netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
    gg2_5 %>% print()

    # ZoomIn
    gg2_6 <- netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
    gg2_6 %>% print()

    ## Compute and visualize the pathway distance in the learned joint manifold
    #> Compute the distance of signaling networks between datasets 1 2
    gg2_7 <- rankSimilarity(cellchat, type = "structural")
    gg2_7 %>% print()


    ## Identify and visualize the conserved and context-specific signaling pathways
    # Compare the overall information flow of each signaling pathway
    gg2_8 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2_9 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg2_8_9 <- gg2_8 + gg2_9
    gg2_8_9 %>% print()

    ## Compare outgoing (or incoming) signaling associated with each cell population
    library(ComplexHeatmap)

    i = 1
    ### combining all the identified signaling pathways from different datasets
    ## outgoing
    pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
    ht2_1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
    ht2_2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
    draw(ht2_1 + ht2_2, ht_gap = unit(0.5, "cm"))

    ## incoming
    ht2_3 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
    ht2_4 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
    draw(ht2_3 + ht2_4, ht_gap = unit(0.5, "cm"))

    ## all
    ht2_5 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
    ht2_6 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
    draw(ht2_5 + ht2_6, ht_gap = unit(0.5, "cm"))

  dev.off() # graphics.off()


##### Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs #####

  ## Identify dysfunctional signaling by using differential expression analysis ##
  # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
  pos.dataset = "EO"
  # define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  # perform differential expression analysis
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  #> Use the joint cell labels from the merged CellChat object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(cellchat, features.name = features.name)
  # extract the ligand-receptor pairs with upregulated ligands in EO
  net.up <- subsetCommunication(cellchat, net = net, datasets = "EO",ligand.logFC = 0.2, receptor.logFC = NULL)
  # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in LO, i.e.,downregulated in LS
  net.down <- subsetCommunication(cellchat, net = net, datasets = "LO",ligand.logFC = -0.1, receptor.logFC = -0.1)

  gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

  #### Bubble_Plot_Summarize ####
  pdf(file = paste0(Save.Path,"/",ProjectName,"_03_Identify_Up_down_signaling_ligand-receptor_pairs_Bubble_Sum.pdf"),
      width = 20,  height = 20
  )
    ## Identify dysfunctional signaling by comparing the communication probabities
    gg3_0 <- netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 45)
    print(gg3_0)

    gg3_1 <- netVisual_bubble(cellchat,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in EO", angle.x = 45, remove.isolate = T)
    gg3_2 <- netVisual_bubble(cellchat,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in EO", angle.x = 45, remove.isolate = T)

    print(gg3_1 + gg3_2)
    gg3_1 %>% print()
    gg3_2 %>% print()

    ## Identify dysfunctional signaling by using differential expression analysis ##
    pairLR.use.up = net.up[, "interaction_name", drop = F]
    gg3_3 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

    pairLR.use.down = net.down[, "interaction_name", drop = F]
    gg3_4 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    print(gg3_3 + gg3_4)
    gg3_3 %>% print()
    gg3_4 %>% print()
  dev.off()

  #### Bubble_Plot_CellTypeAll ####
  pdf(file = paste0(Save.Path,"/",ProjectName,"_03_Identify_Up_down_signaling_ligand-receptor_pairs_Bubble_AllCT.pdf"),
      width = 8,  height = 12
  )
    for (i in 1:length(CellType.set)) {
      try({
    ## Identify dysfunctional signaling by comparing the communication probabities
    p <- netVisual_bubble(cellchat, sources.use = i,  comparison = c(1, 2), angle.x = 45)
    print(p)

    p <- netVisual_bubble(cellchat, sources.use = i,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in EO", angle.x = 45, remove.isolate = T)
    print(p)
    p <- netVisual_bubble(cellchat, sources.use = i,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in EO", angle.x = 45, remove.isolate = T)
    print(p)
    p <- netVisual_bubble(cellchat, sources.use = i, pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    print(p)
    p <- netVisual_bubble(cellchat, sources.use = i, pairLR.use = pairLR.use.down,  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    print(p)

      })
  }
  dev.off()


  #### Chord diagram ####

  pdf(file = paste0(Save.Path,"/",ProjectName,"_03_Identify_Up_down_signaling_ligand-receptor_pairs.pdf"),
      width = 10,  height = 7
  )

    ## Identify dysfunctional signaling by comparing the communication probabities
    netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
    #> Comparing communications on a merged object

    gg3_1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LO", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg3_2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LO", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg3_1 + gg3_2




    pairLR.use.up = net.up[, "interaction_name", drop = F]
    gg3_3 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    #> Comparing communications on a merged object
    pairLR.use.down = net.down[, "interaction_name", drop = F]
    gg3_4 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    #> Comparing communications on a merged object
    gg3_3 + gg3_4

    # Chord diagram
    par(mfrow = c(1,2), xpd=TRUE)
    gg3_5 <-netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    gg3_6 <-netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

    ## Word Cloud
    try({
      # visualize the enriched ligands in the first condition
      computeEnrichmentScore(net.down, species = 'human')

      # visualize the enriched ligands in the second condition
      computeEnrichmentScore(net.up, species = 'human')
    })

    source("FUN_computeEnrichmentScore.R")
    par(mfrow = c(1,2), xpd=TRUE)
    # visualize the enriched ligands in the first condition
    computeEnrichmentScore(net.down, species = 'human')

    # visualize the enriched ligands in the second condition
    computeEnrichmentScore(net.up, species = 'human')


  dev.off()

##### Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram #####
  pdf(file = paste0(Save.Path,"/",ProjectName,"_04_Visually_compare_cell-cell_communication.pdf"),
      width = 12,  height = 7
  )
    # pathways.show <- c("CXCL")
    pathways.show.lt <- object.list[["EO"]]@netP[["pathways"]]

    for (j in 1:length(pathways.show.lt)) {
      pathways.show <- object.list[["EO"]]@netP[["pathways"]][j]



      weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
      par(mfrow = c(1,2), xpd=TRUE)
      for (i in 1:length(object.list)) {
        netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
      }

      #
      par(mfrow = c(1,2), xpd=TRUE)
      ht <- list()
      for (i in 1:length(object.list)) {
        ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
      }
      #> Do heatmap based on a single object
      #>
      #> Do heatmap based on a single object
      ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


      # Chord diagram
      par(mfrow = c(1,2), xpd=TRUE)
      for (i in 1:length(object.list)) {
        netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
      }


      # # Chord diagram
      # group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
      # names(group.cellType) <- levels(object.list[[1]]@idents)
      # par(mfrow = c(1,2), xpd=TRUE)
      # for (i in 1:length(object.list)) {
      #   netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
      # }
      # #> Plot the aggregated cell-cell communication network at the signaling pathway level
      # #> Plot the aggregated cell-cell communication network at the signaling pathway level
      # #>
      #
      # par(mfrow = c(1, 2), xpd=TRUE)
      # # compare all the interactions sending from Inflam.FIB to DC cells
      # for (i in 1:length(object.list)) {
      #   netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:8), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
      # }
      #
      #
      # # compare all the interactions sending from fibroblast to inflamatory immune cells
      # par(mfrow = c(1, 2), xpd=TRUE)
      # for (i in 1:length(object.list)) {
      #   netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
      # }
      #
      #
      # # show all the significant signaling pathways from fibroblast to immune cells
      # par(mfrow = c(1, 2), xpd=TRUE)
      # for (i in 1:length(object.list)) {
      #   netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
      # }


    }

  dev.off()


##### Part V: Compare the signaling gene expression distribution between different datasets #####
  pdf(file = paste0(Save.Path,"/",ProjectName,"_05_Compare_the_signaling_gene_expression_distribution.pdf"),
      width = 10,  height = 7
  )
    # pathways.show <- c("CXCL")
    pathways.show.lt <- object.list[["EO"]]@netP[["pathways"]]

    for (j in 1:length(pathways.show.lt)) {
      pathways.show <- object.list[["EO"]]@netP[["pathways"]][j]

      cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("EO", "LO")) # set factor level
      plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets",
                         colors.ggplot = T) + ggtitle(pathways.show) -> P
      P %>% print()
    }

  dev.off()
  #> The default behaviour of split.by has changed.
  #> Separate violin plots are now plotted side-by-side.
  #> To restore the old behaviour of a single split violin,
  #> set split.plot = TRUE.
  #>
  #> This message will be shown once per session.
  #> Scale for 'y' is already present. Adding another scale for 'y', which will
  #> replace the existing scale.
  #> Scale for 'y' is already present. Adding another scale for 'y', which will
  #> replace the existing scale.
  #> Scale for 'y' is already present. Adding another scale for 'y', which will
  #> replace the existing scale.

# ##### Save the merged CellChat object #####
  ## save rds
  saveRDS(cellchat, file = paste0(Save.Path,"/",ProjectName,"Cell_Cell_Interaction_Multi.rds"))

  ## save RData
  save.image(paste0(Save.Path,"/010_Cell_Cell_Interaction_Multi.RData"))

