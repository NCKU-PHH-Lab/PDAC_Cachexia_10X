##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("AUCell")) install.packages("AUCell"); library(AUCell)

if(!require("patchwork")) install.packages("patchwork"); library(patchwork)  # 使用 patchwork 组合多个图形

##### Set Para #####
## Set Para
Set_Dataset <- "PBMC" # "PBMC" # "SC"
Set_MarkerFile <- "KPCTAMs" # "KPCTAMs" # "HmMuConTAMs"

if(Set_Dataset == "SC"){
  Set_clusters_to_plot <- c("4", "5", "8", "10", "11")
}if(Set_Dataset == "PBMC"){
  Set_clusters_to_plot <- c("2", "6", "13")
}

## Set Export
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_Export <- paste0(Name_time_wo_micro,"_",Set_Dataset,"_",Set_MarkerFile,"_MarkerScore")
Name_ExportFolder <- paste0("Export")
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder


#### Load data ####
## load seuratObject
if(Set_Dataset == "PBMC"){
  load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-12_Results_FollowUpAnalysis/2022-10-17_PBMC_Main/06_Cell_type_annotation.RData")
  seuratObject <- PBMC.combined
}else if(Set_Dataset == "SC"){
  load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-12_Results_FollowUpAnalysis/2022-10-17_SC_Main/06_Cell_type_annotation.RData")
  seuratObject <- SC.combined
}

keep_lst <- c("Set_clusters_to_plot","Set_Dataset","Set_MarkerFile",
              "Name_time_wo_micro","Name_Export","Name_ExportFolder",
              "seuratObject")
rm(list = setdiff(ls(), keep_lst))

colnames(seuratObject@meta.data)[colnames(seuratObject@meta.data) == "celltype"] <- "Cell_Type"

## Preview
source("FUN_Beautify_ggplot.R")
source("FUN_Beautify_UMAP.R")
p1 <- DimPlot(seuratObject, reduction = "umap", group.by = "sample") # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
p2 <- DimPlot(seuratObject, reduction = "umap", group.by = "Cell_Type", label = TRUE) # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
p3 <- DimPlot(seuratObject, reduction = "umap", group.by = "seurat_clusters", label = TRUE) # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
p4 <- DimPlot(seuratObject, reduction = "umap", group.by = "Cachexia") # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
# p1 + p2 + p3 + p4
plot_combined_UMAP <- wrap_plots(list(p1,p2,p3,p4), ncol = 2)  # 根据需要调整列数
plot_combined_UMAP

#### load markers ####
if(Set_MarkerFile == "HmMuConTAMs"){
  file_path_Marker <- "D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/Input_Markers/PMID37914939/41586_2023_6685_MOESM4_ESM_HumanMouse_conserved_TAMsubsets.txt"
  data <- read.table(file_path_Marker, header=TRUE, sep="\t", skip = 1)  # 跳过第一行

}else if(Set_MarkerFile == "KPCTAMs"){
  file_path_Marker <- "D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/Input_Markers/PMID37914939/41586_2023_6685_MOESM4_ESM_scRNAseq_KPC_TAMs_Markers.txt"
  data <- read.table(file_path_Marker, header=TRUE, sep="\t")
  data <- data %>%
    dplyr::rename(Subset = Cell_population) %>%
    dplyr::rename(Gene_Mouse = Gene)
}


subset_genes <- data %>%
  dplyr::select(Subset, Gene_Mouse) %>%
  group_by(Subset) %>%
  summarize(Genes = list(Gene_Mouse))

# print(subset_genes)

# ## Turn to list
# ## Method1
# if(!require("purrr")) install.packages("purrr"); library(purrr)
# standard_list <- purrr::map(subset_genes$Genes, ~ .x)

# ## Method2
# standard_list <- lapply(subset_genes$Genes, function(x) x)

## Turn to list
## Method3
standard_list <- setNames(subset_genes$Genes, subset_genes$Subset)

#### AUCell ####
library(AUCell)

exprMatrix <- GetAssayData(seuratObject, assay = "RNA", slot = "data")
rankings <- AUCell::AUCell_buildRankings(exprMatrix)

expr_genes <- rownames(exprMatrix)

subset_scores <- list() # 初始化用于保存 AUCell 分数的列表
# 计算每个子集的 AUCell 分数
for(subset in names(standard_list)) {
  # 确保基因名称一致性
  gene_set <- standard_list[[subset]]
  valid_genes <- gene_set %in% expr_genes
  if (sum(valid_genes) > 0) {
    gene_set <- gene_set[valid_genes]
    try({
      auc_scores <- AUCell::AUCell_calcAUC(rankings, geneSet = gene_set)
      subset_scores[[subset]] = auc_scores
    })
  } else {
    subset_scores[[subset]] <- NA
  }
}

## Error ## print(subset_scores) # 检查计算结果

# 逐个检查每个子集的 AUCell 分数
for(subset in names(subset_scores)) {
  cat("Subset:", subset, "\n")
  print(head(subset_scores[[subset]]))  # 只打印每个子集的前几个分数
}

# ## Export txt
# write.table(subset_scores, file = "subset_scores.txt", sep = "\t", quote = FALSE)


# ## 将 AUCell 分数添加到 Seurat 对象的元数据中
# ## Method1
# for(subset in names(subset_scores)) {
#   # 假设 subset_scores[[subset]] 是一个向量，包含每个细胞的 AUCell 分数
#   # seuratObject <- AddMetaData(seuratObject, metadata = subset_scores[[subset]], col.name = paste0(subset, "_AUCell"))
#   seuratObject <- AddMetaData(seuratObject, metadata = as.numeric(subset_scores[[subset]]@assays@data@listData[["AUC"]]), col.name = paste0(subset, "_AUCell"))
#
# }

## 将 AUCell 分数添加到 Seurat 对象的元数据中
## Method2
for(subset in names(subset_scores)) {
  # 提取 AUCell 分数和细胞 ID
  auc_data <- subset_scores[[subset]]@assays@data@listData[["AUC"]]
  cell_ids <- colnames(auc_data)

  # 确保 Seurat 对象中的细胞 ID 与 AUCell 分数的细胞 ID 匹配
  if (all(colnames(seuratObject) %in% cell_ids)) {
    # 提取与 Seurat 对象细胞顺序相匹配的 AUCell 分数
    matched_auc_scores <- auc_data[, match(colnames(seuratObject), cell_ids)]

    # 添加匹配的 AUCell 分数到 Seurat 对象的元数据中
    metadata_col_name <- paste0(subset, "_AUCell")
    seuratObject[[metadata_col_name]] <- matched_auc_scores
  } else {
    warning(paste("Not all cells in Seurat object have matching AUCell scores in subset", subset))
  }
}

head(seuratObject@meta.data)

colnames(seuratObject@meta.data) <- gsub("\\.\\.", "+ ", colnames(seuratObject@meta.data))
colnames(seuratObject@meta.data) <- gsub("\\.", " ", colnames(seuratObject@meta.data))



##### Visualization #####
#### UMAP ####
plots <- list()
for(subset in names(subset_scores)) {
  try({
    p <- FeaturePlot(seuratObject, features = paste0(subset, "_AUCell"), reduction = "umap")
    plots[[subset]] <- p
  })
}

# 使用 patchwork 组合多个图形
if(!require("patchwork")) install.packages("patchwork"); library(patchwork)

plot_combined_UMAP_Score <- wrap_plots(plots, ncol = 2)  # 根据需要调整列数
plot_combined_UMAP_Score


#### VlnPlot ####
clusters_to_plot <- Set_clusters_to_plot # 定义要绘制的聚类

# 对于每个子集，绘制特定聚类的 AUCell 分数分布图
plots_2 <- list()
for(subset in names(subset_scores)) {
  feature_to_plot <- paste0(subset, "_AUCell")
  # 创建符合条件的子 Seurat 对象
  subset_seuratObject <- subset(seuratObject, subset = seurat_clusters %in% clusters_to_plot)
  p <- VlnPlot(subset_seuratObject, features = feature_to_plot, group.by = "seurat_clusters")
  plots_2[[subset]] <- p
}

# 使用 patchwork 组合多个图形
if(!require("patchwork")) install.packages("patchwork"); library(patchwork)

plot_combined_Vln_Score <- wrap_plots(plots_2, ncol = 2)  # 根据需要调整列数
plot_combined_Vln_Score

#### Export PDF ####
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult_AUCell.pdf"),
    width = 10, height = 10)
print(plot_combined_UMAP)
print(plot_combined_UMAP_Score)
print(plot_combined_Vln_Score)
dev.off()


####################################################################
#### Module Scores ####
# Add Module Scores to Seurat object
for(subset in names(standard_list)) {
  gene_set <- standard_list[[subset]]
  valid_genes <- gene_set %in% rownames(seuratObject)

  if (sum(valid_genes) >= 2) {  # Ensure there are at least two valid genes
    gene_set <- gene_set[valid_genes]

    # Further reduce the ctrl parameter value
    ctrl_size <- min(min(length(gene_set), 50), nrow(seuratObject) / 2)

    features_list <- list(gene_set)

    # Try to add module score and catch any errors
    tryCatch({
      seuratObject <- AddModuleScore(seuratObject, features = features_list, name = paste0(subset, "_ModuleScore"), ctrl = ctrl_size)
    }, error = function(e) {
      warning(paste("Error occurred while processing gene set", subset, ":", e$message))
    })
  } else {
    warning(paste("Gene set", subset, "has too few genes for Module Score calculation"))
  }
}

colnames(seuratObject@meta.data) <- gsub("\\.\\.", "+ ", colnames(seuratObject@meta.data))
colnames(seuratObject@meta.data) <- gsub("\\.", " ", colnames(seuratObject@meta.data))
colnames(seuratObject@meta.data) <- gsub("_ModuleScore1", "_ModuleScore", colnames(seuratObject@meta.data))

#### UMAP ####
# Example visualization of Module Scores on UMAP
plot_module_scores <- list()
for(subset in names(standard_list)) {
  try({
    feature_to_plot <- paste0(subset, "_ModuleScore") # The suffix "1" is added by Seurat to the Module Score names
    p <- FeaturePlot(seuratObject, features = feature_to_plot, reduction = "umap")
    plot_module_scores[[subset]] <- p
  })

}

if(!require("patchwork")) install.packages("patchwork"); library(patchwork)
plot_combined_UMAP_ModuleScore <- wrap_plots(plot_module_scores, ncol = 2) # Combine and display the plots
plot_combined_UMAP_ModuleScore


#### VlnPlot ####
clusters_to_plot <- Set_clusters_to_plot # 定义要绘制的聚类

# 对于每个子集，绘制特定聚类的 Module Score 分布图
plots_module_score <- list()
for(subset in names(standard_list)) {
  try({
    feature_to_plot <- paste0(subset, "_ModuleScore") # Module Score特征的名称
    # 创建符合条件的子 Seurat 对象
    subset_seuratObject <- subset(seuratObject, subset = seurat_clusters %in% clusters_to_plot)
    p <- VlnPlot(subset_seuratObject, features = feature_to_plot, group.by = "seurat_clusters")
    plots_module_score[[subset]] <- p
  })
}

plot_combined_Vln_ModuleScore <- wrap_plots(plots_module_score, ncol = 2)  # 根据需要调整列数
plot_combined_Vln_ModuleScore

#### Export PDF ####
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult_ModuleScore.pdf"),
    width = 12, height = 10)
print(plot_combined_UMAP)
print(plot_combined_UMAP_ModuleScore)
print(plot_combined_Vln_ModuleScore)
dev.off()


##### Export RData #####
# save.image(paste0(Name_ExportFolder,"/",Name_Export,".RData"))
save(seuratObject, file = paste0(Name_ExportFolder,"/",Name_Export,".RData"))

