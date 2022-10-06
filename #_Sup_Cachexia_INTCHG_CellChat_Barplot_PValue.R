### Ref: Add P-values and Significance Levels to ggplots
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

## https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/plotting-models.html
## https://github.com/kassambara/ggpubr/issues/111

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
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_Results_1stSubmission/2022-09-09_PBMC_Main/06_Cell_type_annotation.RData")

## INTCHG: Interchangeable
  # For PBMC
  scRNA.SeuObj <- PBMC.combined
  SampleType = "PBMC"

  ## For SC
  # scRNA.SeuObj <- SC.combined
  # SampleType = "SC"

# Clean up
rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType")))

## Modify the Cachexia state name
scRNA.SeuObj@meta.data$Cachexia <-  gsub("EO", "EOCX", scRNA.SeuObj@meta.data$Cachexia)
scRNA.SeuObj@meta.data$Cachexia <-  gsub("LO", "PreCX", scRNA.SeuObj@meta.data$Cachexia)


## Order the cell type
  # For PBMC
  scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
                                                 levels =c("Mac1", "Mac2","Mac3","Neu","T","CD4+T","CD8+T",
                                                           "NK","B","Mast","Ery"))
  # # For SC
  # scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
  #                                                levels =c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" , "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
  #                                                          "Fib1", "Fib2", "Fib3"))


##### Current path and new folder setting  #####
TarGene <- c("Gp9")
Version = paste0(Sys.Date(),"_","PBMC_Barplot_PVal_",TarGene)
Save.Path = paste0(getwd(),"/",Version)
dir.create(Save.Path)

##### Extract df #####
## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

Anno.df <- scRNA.SeuObj@meta.data
Anno.df <- data.frame(ID=row.names(Anno.df), Anno.df)

## Save Ori
GeneExp_Ori.df <- GeneExp.df
Anno_Ori.df <- Anno.df
scRNA_Ori.SeuObj <- scRNA.SeuObj

##### Data preprocessing #####
## Extract Target gene and combine to the annotation table
TarGene.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene,] %>% t() %>% as.data.frame()
TarGene.df <- data.frame(ID = row.names(TarGene.df), TarGene.df)
Anno.df <- left_join(Anno.df,TarGene.df)

## Clean up data
Anno.df <- Anno.df[!grepl("Other", Anno.df$celltype),]
GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df$ID]
scRNA.SeuObj <- scRNA.SeuObj[,!grepl("Other", scRNA.SeuObj@meta.data[["celltype"]] )]


##### Visualize the expression profile ####
source("FUN_Beautify_ggplot.R")

## Set y position
LabelY <- max(Anno.df[,TarGene])

#
# ##### Export PDF #####
# pdf(file = paste0(Save.Path,"/",Version,"_BarplotMulti.pdf"),width = 13, height = 13 )
#   plt.ManyGroup2_Sum
#   plt.ManyGroup3_Sum
# dev.off()

#### Cell type & EO LO ####
# Box plot facetted by "celltype"
plt.ManyGroup2 <- ggboxplot(Anno.df, x = "celltype", y = TarGene,
                            color = "Cachexia", # palette = "jco",
                            add = "jitter", # short.panel.labs = T
                            ) +
  ylim(0, LabelY*1.2)


## Beautify ggPlot
plt.ManyGroup2 <- plt.ManyGroup2 +
    theme(axis.text.x = element_text(face="bold",  size = 17,angle = 90, hjust = 1, vjust = .5), # Change the size along the x axis
          axis.text.y = element_text(face="bold",size =  17), # Change the size along the y axis

          axis.line = element_line(colour = "black", size = 1.5, linetype = "solid"),
          # axis.title = element_text(size = rel(AxisTitleSize),face="bold"),
          # plot.title = element_text(color="black",
          #                           size=TitleSize,
          #                           face="bold.italic",
          #                           hjust = TH,vjust =TV), # margin = margin(t = 0.5, b = -7),
          # #     plot.background = element_rect(fill = 'chartreuse'),

          legend.title = element_text(size= 20, color = "black", face="bold"),
          legend.text = element_text(colour="black", size= 17,face="bold"),
          legend.background = element_rect(fill = alpha("white", 0.5)),
          legend.position = c(0.11, 0.96),legend.direction= "horizontal",
          #     plot.text = element_text(size = 20),
          #     aspect.ratio=AspRat
          )

plt.ManyGroup2

## Add PValue
plt.ManyGroup2 <- plt.ManyGroup2 +
                  stat_compare_means(aes(group = Cachexia),
                                     label =  "p.signif",label.x = 1.5,
                                     label.y = LabelY*1.2*0.9,
                                     method = "wilcox.test", size = 7)
plt.ManyGroup2

# # Use only p.format as label. Remove method name.
# p + stat_compare_means(label = "p.format", method = "wilcox.test", size = 7)
# # Or use significance symbol as label
# p + stat_compare_means(label =  "p.signif",label.x = 1.5, label.y = LabelY*0.9, method = "wilcox.test", size = 7)


##### 20221005 Try #####
plt.ManyGroup2_2_Sum <- ggboxplot(Anno.df, x = "celltype", y = TarGene,
                                  color = "Cachexia", #palette = "jco",
                                  add = "jitter",)
plt.ManyGroup2_2_Sum + stat_compare_means(aes(group = Cachexia), label = "p.format") -> TTT1
plt.ManyGroup2_2_Sum + stat_compare_means(aes(group = Cachexia), label = "p.signif") -> TTT2

TTT1/TTT2 -> TTT3
TTT3
TTT3/TTT2

# ## Violin
# ## http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/
# plt.ManyGroup2_2_Sum <- ggviolin(Anno_Sum.df, x = "celltype", y = TarGene,
#                                   color = "Cachexia", palette = "jco",
#                                   add = "jitter",)
# plt.ManyGroup2_2_Sum + stat_compare_means(aes(group = Cachexia), label = "p.format")
# plt.ManyGroup2_2_Sum + stat_compare_means(aes(group = Cachexia), label = "p.signif")


##### Save RData #####
save.image(paste0(Save.Path,"/",Version,".RData"))



##### Plot UMAP #####

DimPlot(scRNA.SeuObj, reduction = "umap", group.by = "celltype",label = T,label.size = 9) %>%
  BeautifyggPlot(.,LegPos = c(1.1, 0.5),AxisTitleSize=1.1) -> plt.UMAP1
plt.UMAP1

DimPlot(scRNA.SeuObj, reduction = "umap", group.by = "celltype",
        split.by = "sample",label = T,label.size = 5, ncol = 2) %>%
  BeautifyggPlot(.,LegPos = c(1.05, 0.5),AxisTitleSize=1.1 ,TV= 1,TH= 0.3) -> plt.UMAP2
plt.UMAP2


## https://github.com/satijalab/seurat/issues/2937
FeaturePlot(scRNA.SeuObj, features = TarGene, min.cutoff = "q9",
            split.by = "sample",ncol = 2, coord.fixed = 1)  & theme(legend.position = c(0.9,0.3)) -> plt.UMAP3


FeaturePlot(scRNA.SeuObj, features = TarGene, min.cutoff = "q9",
            split.by = "Cachexia",ncol = 2, coord.fixed = 1) & theme(legend.position = c(0.9,0.2)) -> plt.UMAP4

##### Export PDF #####
pdf(file = paste0(Save.Path,"/",Version,"_UMAP.pdf"),width = 15, height = 10 )
  plt.UMAP1
  plt.UMAP2
  plt.UMAP3
  plt.UMAP4
dev.off()

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

