##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("AUCell")) install.packages("AUCell"); library(AUCell)


#### Load data ####
## load seuratObject
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-12_Results_FollowUpAnalysis/2022-10-17_PBMC_Main/06_Cell_type_annotation.RData")
seuratObject <- PBMC.combined
rm(list = setdiff(ls(), "seuratObject"))

colnames(seuratObject@meta.data)[colnames(seuratObject@meta.data) == "celltype"] <- "Cell_Type"

## Preview
source("FUN_Beautify_ggplot.R")
source("FUN_Beautify_UMAP.R")
p1 <- DimPlot(seuratObject, reduction = "umap", group.by = "sample") # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
p2 <- DimPlot(seuratObject, reduction = "umap", group.by = "Cell_Type", label = TRUE) # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
p3 <- DimPlot(seuratObject, reduction = "umap", group.by = "seurat_clusters") # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
p4 <- DimPlot(seuratObject, reduction = "umap", group.by = "Cachexia") # %>% BeautifyggPlot(.,LegPos = c(0.85, 0.15),AxisTitleSize=1.1)
p1 + p2 + p3 + p4

## load markers
file_path <- "D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/Input_Markers/PMID37914939/41586_2023_6685_MOESM4_ESM_HumanMouse_conserved_TAMsubsets.txt"
data <- read.table(file_path, header=TRUE, sep="\t", skip = 1)  # 跳过第一行

file_path <- "D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/Input_Markers/PMID37914939/41586_2023_6685_MOESM4_ESM_scRNAseq_KPC_TAMs_Markers.txt"
data <- read.table(file_path, header=TRUE, sep="\t")  # 跳过第一行
data <- data %>%
  dplyr::rename(Subset = Cell_population) %>%
  dplyr::rename(Gene_Mouse = Gene)

subset_genes <- data %>%
  dplyr::select(Subset, Gene_Mouse) %>%
  group_by(Subset) %>%
  summarize(Genes = list(Gene_Mouse))

# print(subset_genes)

## Turn to list
# ## Method1
# if(!require("purrr")) install.packages("purrr"); library(purrr)
# standard_list <- purrr::map(subset_genes$Genes, ~ .x)

# ## Method2
# standard_list <- lapply(subset_genes$Genes, function(x) x)

## Method3
standard_list <- setNames(subset_genes$Genes, subset_genes$Subset)

#### AUCell ####
library(AUCell)

# 提取基因表达矩阵
exprMatrix <- GetAssayData(seuratObject, assay = "RNA", slot = "data")

# 创建基因排名
rankings <- AUCell::AUCell_buildRankings(exprMatrix)

expr_genes <- rownames(exprMatrix)
# 初始化用于保存 AUCell 分数的列表
subset_scores <- list()

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

# 检查计算结果
print(subset_scores)

# 将 AUCell 分数保存到文件中
write.table(subset_scores, file = "subset_scores.txt", sep = "\t", quote = FALSE)

# 逐个检查每个子集的 AUCell 分数
for(subset in names(subset_scores)) {
  cat("Subset:", subset, "\n")
  print(head(subset_scores[[subset]]))  # 只打印每个子集的前几个分数
}




# 将 AUCell 分数添加到 Seurat 对象的元数据中
for(subset in names(subset_scores)) {
  # 假设 subset_scores[[subset]] 是一个向量，包含每个细胞的 AUCell 分数
  # seuratObject <- AddMetaData(seuratObject, metadata = subset_scores[[subset]], col.name = paste0(subset, "_AUCell"))
  seuratObject <- AddMetaData(seuratObject, metadata = as.numeric(subset_scores[[subset]]@assays@data@listData[["AUC"]]), col.name = paste0(subset, "_AUCell"))

}

colnames(seuratObject@meta.data) <- gsub("\\.\\.", "+ ", colnames(seuratObject@meta.data))
colnames(seuratObject@meta.data) <- gsub("\\.", " ", colnames(seuratObject@meta.data))

# 绘制带有 AUCell 分数的 UMAP 组合图
plots <- list()
for(subset in names(subset_scores)) {
  try({
    p <- FeaturePlot(seuratObject, features = paste0(subset, "_AUCell"), reduction = "umap")
    plots[[subset]] <- p
  })
}

# 使用 patchwork 组合多个图形
if(!require("patchwork")) install.packages("patchwork"); library(patchwork)

combined_plot <- wrap_plots(plots, ncol = 3)  # 根据需要调整列数

# 显示组合图
combined_plot

combined_plot2 <- wrap_plots(list(p1,p2,p3,p4), ncol = 2)  # 根据需要调整列数
combined_plot2
