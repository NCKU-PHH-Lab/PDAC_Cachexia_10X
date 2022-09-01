### Ref: Add P-values and Significance Levels to ggplots
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if(!require("Seurat")) install.packages("Seurat")
library(Seurat)
if(!require("ggpubr")) install.packages("ggpubr")
library(ggpubr)

##### Load RData* #####
# load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-06-06_CC_PBMC/06_Cell_type_annotation.RData")
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-08-13_PBMC_Main/06_Cell_type_annotation.RData")

# Clean up
scRNA.SeuObj <- PBMC.combined
rm(list=setdiff(ls(), "scRNA.SeuObj"))


##### Current path and new folder setting  #####
TarGene <- "Chil3"
Version = paste0(Sys.Date(),"_","PBMC_Barplot_PVal_",TarGene)
Save.Path = paste0(getwd(),"/",Version)
dir.create(Save.Path)

##### Extract df #####
GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
Anno.df <- scRNA.SeuObj@meta.data
Anno.df <- data.frame(ID=row.names(Anno.df), Anno.df)

##### Data preprocessing #####
## Extract Target gene and combine to the annotation table
TarGene.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene,] %>% t() %>% as.data.frame()
TarGene.df <- data.frame(ID = row.names(TarGene.df), TarGene.df)
Anno.df <- left_join(Anno.df,TarGene.df)

## Clean up data
Anno.df <- Anno.df[!grepl("Other", Anno.df$celltype),]


##### Visualize the expression profile ####
ggboxplot(Anno.df, x = "celltype", y = TarGene, color = "celltype",
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(myeloma$DEPDC1), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova" )+ #, label.y = 1600)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                      # Pairwise comparison against all



# ## Example
# # Load myeloma data from GitHub
# myeloma <- read.delim("https://raw.githubusercontent.com/kassambara/data/master/myeloma.txt")
# # Perform the test
# compare_means(DEPDC1 ~ molecular_group,  data = myeloma,
#               ref.group = ".all.", method = "t.test")
#
# # Visualize the expression profile
# ggboxplot(myeloma, x = "molecular_group", y = "DEPDC1", color = "molecular_group",
#           add = "jitter", legend = "none") +
#   rotate_x_text(angle = 45)+
#   geom_hline(yintercept = mean(myeloma$DEPDC1), linetype = 2)+ # Add horizontal line at base mean
#   stat_compare_means(method = "anova", label.y = 1600)+        # Add global annova p-value
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = ".all.")                      # Pairwise comparison against all
#
