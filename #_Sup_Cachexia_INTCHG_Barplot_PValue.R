### Ref: Add P-values and Significance Levels to ggplots
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

## https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/plotting-models.html
## https://github.com/kassambara/ggpubr/issues/111


# ##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

# Save.Path <- c("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_PBMC_Main")
Save.Path <- c("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-12_Results_FollowUpAnalysis/2022-10-17_PBMC_Main")


##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if(!require("Seurat")) install.packages("Seurat")
library(Seurat)
if(!require("ggpubr")) install.packages("ggpubr")
library(ggpubr)

##### Load RData* #####
# load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-04_PBMC_Main/06_Cell_type_annotation.RData")
load(paste0(Save.Path,"/06_Cell_type_annotation.RData"))

## INTCHG: Interchangeable
### SubType Setting
### Order the cell type
## INTCHG: Interchangeable
## SubType Setting
if(SampleType == "PBMC"){
  ## For PBMC
  scRNA.SeuObj <- PBMC.combined

  # Order the cell type
  CellType.Order = c("Mac1", "Mac2","Mac3","Neu","T","CD4+T","CD8+T","NK","B","Mast","Ery")
  scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
                                                 levels = CellType.Order)


}else if(SampleType == "SC"){
  ## For SC
  scRNA.SeuObj <- SC.combined

  # Order the cell type
  CellType.Order = c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" , "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
                     "Fib1", "Fib2", "Fib3")
  scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
                                                 levels = CellType.Order)

}

# Clean up
rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Save.Path", "CellType.Order")))


##### Current path and new folder setting  #####
TarGene <- "Il4"
Version = paste0(Sys.Date(),"_",SampleType,"_Barplot_PVal_",TarGene)
SaveSup.Path = paste0(getwd(),"/",Version)
dir.create(SaveSup.Path)

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

#### EOCX PreCX ####
plt.2Group <- ggboxplot(Anno.df, x = "Cachexia", y = TarGene,
                        color = "Cachexia", palette = "jco",
                        add = "jitter")
#  Add p-value
plt.2Group1 <- plt.2Group + stat_compare_means(label.x = 0.6, label.y = LabelY*1.05, size = 7)
# Change method
plt.2Group2 <- plt.2Group + stat_compare_means(method = "t.test")

plt.2Group1+plt.2Group2
plt.2Group1 %>% BeautifyggPlot(LegPos = c(1.08, 0.5),LegTitleSize=17 ,LegTextSize = 15,
                               XtextSize=25,  YtextSize=25,
                               AxisTitleSize=2) -> plt.2Group1
plt.2Group1

#### Cell type & EOCX PreCX ####
# Visualize: Specify the comparisons you want
## Ref(size): ## https://github.com/kassambara/ggpubr/issues/65
Anno.df$sample <- factor(Anno.df$sample ,levels =c("EOCX.M","EOCX.F", "PreCX.M","PreCX.F"))
my_comparisons <- list(  c("EOCX.M", "EOCX.F"), c("EOCX.F", "PreCX.M"), c("PreCX.M", "PreCX.F"), c("EOCX.M", "PreCX.M"), c("EOCX.F", "PreCX.F"), c("EOCX.M", "PreCX.F"))
plt.FewGroup1 <-  ggboxplot(Anno.df, x = "sample", y = TarGene,
                            color = "sample", palette = "jco",
                            add = "jitter", legend = "none")+
                    stat_compare_means(comparisons = my_comparisons, line.size=5,
                                       method = "wilcox.test",label = "p.signif",size = 7)+ #, label.y = c(29, 35, 40))+
                    stat_compare_means(label.x = 0.8, label.y = LabelY*1.8, size = 7)

plt.FewGroup1 %>% BeautifyggPlot(LegPos = c(1.08, 0.5),LegTitleSize=17 ,LegTextSize = 15,
                                XtextSize=25,  YtextSize=25,
                                AxisTitleSize=2) -> plt.FewGroup1
plt.FewGroup1

#### cell type ####
plt.ManyGroup1 <- ggboxplot(Anno.df, x = "celltype", y = TarGene, color = "celltype",
                            add = "jitter", legend = "none") +
                   rotate_x_text(angle = 45)+
                   geom_hline(yintercept = mean(Anno.df[,TarGene]), linetype = 2)+ # Add horizontal line at base mean
                   stat_compare_means(method = "anova", label.x = 1.5, label.y = LabelY*1.05, size = 7)+        # Add global annova p-value
                   stat_compare_means(label = "p.signif", method = "wilcox.test",
                                      ref.group = ".all.", size = 7)                      # Pairwise comparison against all

plt.ManyGroup1 %>% BeautifyggPlot(LegPos = c(1.1, 0.5),LegTitleSize=17 ,LegTextSize = 15,
                                 XtextSize=17,  YtextSize=17, xangle =90,
                                 AxisTitleSize=2) -> plt.ManyGroup1
plt.ManyGroup1




#### Cell type & EOCX PreCX ####
# Box plot facetted by "celltype"
plt.ManyGroup2 <- ggboxplot(Anno.df, x = "Cachexia", y = TarGene,
                            color = "Cachexia", palette = "jco",
                            add = "jitter",
                            facet.by = "celltype", short.panel.labs = T) +
  ylim(0, LabelY*1.2)+
  stat_compare_means(label.x = 1.1, label.y = LabelY*1.1, size = 4)

# # Use only p.format as label. Remove method name.
# p + stat_compare_means(label = "p.format", method = "wilcox.test", size = 7)
# # Or use significance symbol as label
# p + stat_compare_means(label =  "p.signif",label.x = 1.5, label.y = LabelY*0.9, method = "wilcox.test", size = 7)

plt.ManyGroup2 %>% BeautifyggPlot(LegPos = c(0.5, 1.1),LegTitleSize=17 ,LegTextSize = 15,
                                  LegBox = "horizontal",LegDir="horizontal",
                                  XtextSize=17,  YtextSize=17, xangle =0,
                                  XaThick=0,  YaThick=0,OL_Thick = 1.5,
                                  AxisTitleSize=2) -> plt.ManyGroup2
plt.ManyGroup2 <- plt.ManyGroup2 + stat_compare_means(label =  "p.signif",label.x = 1.5, label.y = LabelY*0.9, method = "wilcox.test", size = 7)

    ## Add all cell type plot to multiple plots
    Anno_Sum.df <- Anno.df
    Anno_Sum.df$celltype <- "Total"
    Anno_Sum.df <- rbind(Anno.df,Anno_Sum.df)

    Anno_Sum.df$sample <- factor(Anno_Sum.df$sample ,levels =c("EOCX.M","EOCX.F", "PreCX.M","PreCX.F"))
    Anno_Sum.df$celltype <- factor(Anno_Sum.df$celltype ,
                                   levels =c("Total",CellType.Order))

    # Box plot facetted by "celltype"
    plt.ManyGroup2_Sum <- ggboxplot(Anno_Sum.df, x = "Cachexia", y = TarGene,
                                color = "Cachexia", palette = "jco",
                                add = "jitter",
                                facet.by = "celltype", short.panel.labs = T) +
      ylim(0, LabelY*1.2)+
      stat_compare_means(label.x = 1.1, label.y = LabelY*1.1, size = 5)

    plt.ManyGroup2_Sum %>% BeautifyggPlot(LegPos = c(0.5, 1.1),LegTitleSize=17 ,LegTextSize = 15,
                                      LegBox = "horizontal",LegDir="horizontal",
                                      XtextSize=17,  YtextSize=17, xangle =0,
                                      XaThick=0,  YaThick=0,OL_Thick = 1.5,
                                      AxisTitleSize=2) -> plt.ManyGroup2_Sum
    plt.ManyGroup2_Sum <- plt.ManyGroup2_Sum + stat_compare_means(label =  "p.signif",label.x = 1.5, label.y = LabelY*0.9, method = "wilcox.test", size = 7)
    plt.ManyGroup2_Sum


#### Cell type & EOCX.M EOCX.F PreCX.M PreCX.F ####
# Box plot facetted by "celltype"
my_comparisons <- list(  c("EOCX.M", "EOCX.F"), c("EOCX.F", "PreCX.M"), c("PreCX.M", "PreCX.F"), c("EOCX.F", "PreCX.F"), c("EOCX.M", "PreCX.M"), c("EOCX.M", "PreCX.F"))

plt.ManyGroup3 <- ggboxplot(Anno.df, x = "sample", y = TarGene,
                            color = "sample", palette = "jco",
                            add = "jitter",
                            facet.by = "celltype", short.panel.labs = TRUE) +
  ylim(0, LabelY*2)+
  stat_compare_means(label.x = 1.25, label.y = LabelY*1.9, size = 4)

plt.ManyGroup3 %>% BeautifyggPlot(LegPos = c(0.5, 1.1),LegTitleSize=17 ,LegTextSize = 15,
                                  LegBox = "horizontal",LegDir="horizontal",
                                  XtextSize=15,  YtextSize=17, xangle =90,
                                  XaThick=0,  YaThick=0,OL_Thick = 1.5,
                                  AxisTitleSize=2) -> plt.ManyGroup3

# Use only p.format as label. Remove method name.
plt.ManyGroup3 <- plt.ManyGroup3 + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")
plt.ManyGroup3

    ## Add all cell type plot to multiple plots
    # Box plot facetted by "celltype"
    plt.ManyGroup3_Sum <- ggboxplot(Anno_Sum.df, x = "sample", y = TarGene,
                                color = "sample", palette = "jco",
                                add = "jitter",
                                facet.by = "celltype", short.panel.labs = TRUE) +
      ylim(0, LabelY*2)+
      stat_compare_means(label.x = 1.25, label.y = LabelY*1.9, size = 5)

    plt.ManyGroup3_Sum %>% BeautifyggPlot(LegPos = c(0.5, 1.1),LegTitleSize=17 ,LegTextSize = 15,
                                      LegBox = "horizontal",LegDir="horizontal",
                                      XtextSize=15,  YtextSize=17, xangle =90,
                                      XaThick=0,  YaThick=0,OL_Thick = 1.5,
                                      AxisTitleSize=2) -> plt.ManyGroup3_Sum

    # Use only p.format as label. Remove method name.
    plt.ManyGroup3_Sum <- plt.ManyGroup3_Sum + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")
    plt.ManyGroup3_Sum



##### Export PDF #####
pdf(file = paste0(SaveSup.Path,"/",Version,"_Barplot.pdf"),width = 10, height = 10 )
  plt.2Group1 %>% print()
  plt.FewGroup1 %>% print()
  plt.ManyGroup1 %>% print()
  plt.ManyGroup2 %>% print()
  plt.ManyGroup2_Sum %>% print()
  plt.ManyGroup3 %>% print()
  plt.ManyGroup3_Sum %>% print()
dev.off()

##### Export PDF #####
pdf(file = paste0(SaveSup.Path,"/",Version,"_BarplotMulti.pdf"),width = 13, height = 13 )
  plt.ManyGroup2_Sum %>% print()
  plt.ManyGroup3_Sum %>% print()
dev.off()



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
pdf(file = paste0(SaveSup.Path,"/",Version,"_UMAP.pdf"),width = 15, height = 10 )
  plt.UMAP1 %>% print()
  plt.UMAP2 %>% print()
  plt.UMAP3 %>% print()
  plt.UMAP4 %>% print()
dev.off()

# ##### Save RData #####
# save.image(paste0(Save.Path,"/",Version,".RData"))

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

