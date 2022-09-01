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
TarGene <- "Chil1"
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

#### EO LO ####
p <- ggboxplot(Anno.df, x = "Cachexia", y = TarGene,
               color = "Cachexia", palette = "jco",
               add = "jitter")
#  Add p-value
p1 <- p + stat_compare_means()
# Change method
p2 <- p + stat_compare_means(method = "t.test")

p1+p2
p1

#### Cell type & EO LO ####
# Visualize: Specify the comparisons you want
Anno.df$sample <- factor(Anno.df$sample ,levels =c("EO.M","EO.F", "LO.M","LO.F"))
my_comparisons <- list(  c("EO.M", "EO.F"), c("EO.F", "LO.M"), c("LO.M", "LO.F"), c("EO.F", "LO.F"), c("EO.M", "LO.M"), c("EO.M", "LO.F"))
ggboxplot(Anno.df, x = "sample", y = TarGene,
          color = "sample", palette = "jco",
          add = "jitter", legend = "none")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+ #, label.y = c(29, 35, 40))+
  stat_compare_means(label.x = 0.6, label.y = 25)



#### cell type ####
ggboxplot(Anno.df, x = "celltype", y = TarGene, color = "celltype",
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(Anno.df[,TarGene]), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 400)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.")                      # Pairwise comparison against all

#### Cell type & EO LO ####
# Box plot facetted by "celltype"
p <- ggboxplot(Anno.df, x = "Cachexia", y = TarGene,
               color = "Cachexia", palette = "jco",
               add = "jitter",
               facet.by = "celltype", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox.test")
# Or use significance symbol as label
p + stat_compare_means(label =  "p.signif", label.x = 1.5, method = "wilcox.test")


#### Cell type & EO.M EO.F LO.M LO.F ####
# Box plot facetted by "celltype"
my_comparisons <- list(  c("EO.M", "EO.F"), c("EO.F", "LO.M"), c("LO.M", "LO.F"), c("EO.F", "LO.F"), c("EO.M", "LO.M"), c("EO.M", "LO.F"))

p <- ggboxplot(Anno.df, x = "sample", y = TarGene,
               color = "sample", palette = "jco",
               add = "jitter",
               facet.by = "celltype", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")




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
