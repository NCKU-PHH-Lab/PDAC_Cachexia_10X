### Ref: Add P-values and Significance Levels to ggplots
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

## https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/plotting-models.html
## https://github.com/kassambara/ggpubr/issues/111

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


#### Installation and load the required libraries ####
#### Basic installation ####
## Package.set
Package.set <- c("tidyverse","CellChat","patchwork","NMF","ggalluvial","Seurat","ggpubr")
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


##### Load RData* #####
## Load Seurat.Obj
Load.Path <- "D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_Results_1stSubmission/2022-09-09_PBMC_Main"
load(paste0(Load.Path,"/06_Cell_type_annotation.RData"))
# load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_Results_1stSubmission/2022-09-09_PBMC_Main/06_Cell_type_annotation.RData")

  ## INTCHG: Interchangeable
  # For PBMC
  scRNA.SeuObj <- PBMC.combined
  SampleType = "PBMC"

  ## For SC
  # scRNA.SeuObj <- SC.combined
  # SampleType = "SC"

  # Clean up
  rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Load.Path")))

  ## Modify the Cachexia state name
  scRNA.SeuObj@meta.data$Cachexia <-  gsub("EO", "EOCX", scRNA.SeuObj@meta.data$Cachexia)
  scRNA.SeuObj@meta.data$Cachexia <-  gsub("LO", "PreCX", scRNA.SeuObj@meta.data$Cachexia)


  ## Order the cell type*
  # For PBMC
  CellType.Order = c("Mac1", "Mac2","Mac3","Neu","T","CD4+T","CD8+T","NK","B","Mast","Ery")
  scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
                                                 levels = CellType.Order)
  # # For SC
  # scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
  #                                                levels =c("Duc1", "Duc2", "Duc3", "Duc4", "Duc5", "Duc6" , "Mac1", "Mac2", "Mac3", "Mac4", "Mac5",
  #                                                          "Fib1", "Fib2", "Fib3"))


## Load CellChat rds
CCDBType = "ECM"
cellchat.EO <- readRDS(paste0(Load.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_EO_CellChat.rds"))
cellchat.LO <- readRDS(paste0(Load.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_LO_CellChat.rds"))

object.list <- list(LO = cellchat.LO, EO = cellchat.EO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rm(object.list, cellchat.EO, cellchat.LO)


##### Current path and new folder setting  #####
Version = paste0(Sys.Date(),"_", SampleType, "_CellChat_Barplot_PVal_",CCDBType)
Save.Path = paste0(getwd(),"/",Version)
dir.create(Save.Path)



##### Pathway and TarGene Setting  #####
pathways.show1 <-cellchat@netP[["EO"]][["pathways"]]
pathways.show2 <-cellchat@netP[["LO"]][["pathways"]]
pathways.show <- unique(pathways.show1,pathways.show2)
rm(pathways.show1,pathways.show2)

TarGene_All <- cellchat@data.signaling %>% rownames
LR.df <- rbind(cellchat@LR[["EO"]][["LRsig"]],cellchat@LR[["LO"]][["LRsig"]])

##### Summarize all signal #####
pdf(file = paste0(Save.Path,"/",Version,"_LR_BarplotMulti.pdf"),width = 15, height = 20 )
SummaryTable.df <-  as.data.frame(matrix(nrow=0,ncol=10))
colnames(SummaryTable.df) <- c( "celltype", ".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method","pathway_name")

for (j in 1:length(pathways.show)) {


  LR_Tar.df <- LR.df[LR.df$pathway_name == pathways.show[j],]

  # TarGene <- c("Vwf","Itga2b","Itgb3","Gp9")
  TarGene <- c(LR_Tar.df$ligand, LR_Tar.df$receptor) %>% unique()
  # TarGene <- c(LR_Tar.df$ligand[1], LR_Tar.df$receptor) %>% unique()
  TarGene <-intersect(TarGene,row.names(GeneExp.df))

  rm(LR_Tar.df)

  ##### Extract df #####
  ## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix
  ## Set y position
  LabelY <- max(GeneExp.df) %>% ceiling()

  Anno.df <- scRNA.SeuObj@meta.data
  Anno.df <- data.frame(ID=row.names(Anno.df), Anno.df)

  ## Save Ori
  GeneExp_Ori.df <- GeneExp.df
  Anno_Ori.df <- Anno.df
  scRNA_Ori.SeuObj <- scRNA.SeuObj

  ##### Data preprocessing #####
  ## Extract Target gene and combine to the annotation table
  TarGene.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene,] %>% t() %>% as.data.frame()
  max(TarGene.df)

  TarGene.df <- data.frame(ID = row.names(TarGene.df), TarGene.df)
  Anno.df <- left_join(Anno.df,TarGene.df)

  ## Clean up data
  Anno.df <- Anno.df[!grepl("Other", Anno.df$celltype),]
  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df$ID]
  scRNA.SeuObj <- scRNA.SeuObj[,!grepl("Other", scRNA.SeuObj@meta.data[["celltype"]] )]


  # ## Summary Statistic Table
  # #(Ori)# SummaryTable.df <- compare_means( Vwf ~ Cachexia, data = Anno.df, group.by = "celltype"	)
  #
  # # ## Error (Solved)
  # # TTT <- compare_means( Anno.df[,TarGene[1]] ~ Cachexia, data = Anno.df, group.by = "celltype"	)
  #
  # # https://stackoverflow.com/questions/44776446/compare-means-must-resolve-to-integer-column-positions-not-a-symbol-when-u
  # # convert string column name to name/symbol
  # f <- paste0(TarGene[1]," ~ Cachexia") # f <- "Vwf ~ Cachexia"
  # SummaryTable.df <- do.call("compare_means", list(as.formula(f), data=Anno.df, group.by = "celltype"))
  # rm(f)
  # SummaryTable.df$celltype <- factor(SummaryTable.df$celltype  ,levels = CellType.Order)
  # SummaryTable.df <- SummaryTable.df[order(SummaryTable.df$celltype), ]

  ##### Summary Statistic Table #####
  SummaryTable_Sub.df <-  as.data.frame(matrix(nrow=0,ncol=9))
  colnames(SummaryTable_Sub.df) <- c( "celltype", ".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method"  )
  for (i in 1:length(TarGene)) {
    # https://stackoverflow.com/questions/44776446/compare-means-must-resolve-to-integer-column-positions-not-a-symbol-when-u
    # convert string column name to name/symbol
    f <- paste0(TarGene[i]," ~ Cachexia") # f <- "Vwf ~ Cachexia"
    SummaryTable_Temp.df <- do.call("compare_means", list(as.formula(f), data=Anno.df, group.by = "celltype"))
    rm(f)
    SummaryTable_Temp.df$celltype <- factor(SummaryTable_Temp.df$celltype  ,levels = CellType.Order)
    SummaryTable_Temp.df <- SummaryTable_Temp.df[order(SummaryTable_Temp.df$celltype), ]

    if(c("****") %in% SummaryTable_Temp.df$p.signif || c("***") %in% SummaryTable_Temp.df$p.signif){
      SummaryTable_Sub.df <- rbind(SummaryTable_Sub.df,SummaryTable_Temp.df)
    }else{
      SummaryTable_Sub.df <- SummaryTable_Sub.df
    }

    rm(SummaryTable_Temp.df)



    # if(i==1){
    #   SummaryTable_Sub.df <- SummaryTable_Temp.df
    # }else{
    #   SummaryTable_Sub.df <- rbind(SummaryTable_Sub.df,SummaryTable_Temp.df)
    #   rm(SummaryTable_Temp.df)
    # }
  }
  SummaryTable_Sub.df$pathway_name <- pathways.show[j]
  SummaryTable.df <- rbind(SummaryTable.df, SummaryTable_Sub.df)

  TarGene <- SummaryTable_Sub.df$.y. %>% unique()

  ##### Visualize the expression profile #####
  source("FUN_Beautify_ggplot.R")

  # ##### Export PDF #####
  # pdf(file = paste0(Save.Path,"/",Version,"_BarplotMulti.pdf"),width = 13, height = 13 )
  #   plt.ManyGroup2_Sum
  #   plt.ManyGroup3_Sum
  # dev.off()

  #### Cell type & EO LO ####
  ## BarPlot
  for (i in 1:length(TarGene)) {
    ## Main ggPlot
    plt.ManyGroup2 <- ggboxplot(Anno.df, x = "celltype", y = TarGene[i],
                                color = "Cachexia", # palette = "jco",
                                add = "jitter", # short.panel.labs = T
    ) + ylim(0, LabelY*1.2)

    ## Beautify ggPlot
    plt.ManyGroup2 <- plt.ManyGroup2 +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(face="bold",size =  17), # Change the size along the y axis

            axis.line = element_line(colour = "black", size = 1.5, linetype = "solid"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 17,face="bold"),
            legend.position = "none",
      )

    ## Add PValue
    plt.ManyGroup2 <- plt.ManyGroup2 +
      stat_compare_means(aes(group = Cachexia),
                         label =  "p.signif",label.x = 1.5,
                         label.y = LabelY*1.2*0.9,
                         method = "wilcox.test", size = 7)

    if(i==1){
      plt.ManyGroupSum <- plt.ManyGroup2
      plt.ManyGroupSum <- plt.ManyGroupSum  + ggtitle(pathways.show[j]) + theme(
        plot.title = element_text(color="black", size=20, face="bold.italic"))+
        theme(legend.title = element_text(size= 17, color = "black", face="bold"),
              legend.text = element_text(colour="black", size= 17,face="bold"),
              legend.background = element_rect(fill = alpha("white", 0.5)),
              legend.position = c(0.11, 0.8), # legend.position = c(0.11, 0.96),
              legend.direction= "horizontal",
        )

    }else if(i>=1 && i<length(TarGene)){

      plt.ManyGroupSum <- plt.ManyGroupSum/plt.ManyGroup2

    }else{

      ## Main ggPlot
      plt.ManyGroup2 <- ggboxplot(Anno.df, x = "celltype", y = TarGene[i],
                                  color = "Cachexia", # palette = "jco",
                                  add = "jitter", # short.panel.labs = T
      ) + ylim(0, LabelY*1.2)

      ## Beautify ggPlot
      plt.ManyGroup2 <- plt.ManyGroup2 +
        theme(axis.text.x = element_text(face="bold",  size = 17,angle = 90, hjust = 1, vjust = .5), # Change the size along the x axis
              axis.text.y = element_text(face="bold",size =  17), # Change the size along the y axis

              axis.line = element_line(colour = "black", size = 1.5, linetype = "solid"),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = 17,face="bold"),

              # axis.title = element_text(size = rel(AxisTitleSize),face="bold"),
              # plot.title = element_text(color="black",
              #                           size=TitleSize,
              #                           face="bold.italic",
              #                           hjust = TH,vjust =TV), # margin = margin(t = 0.5, b = -7),
              # #     plot.background = element_rect(fill = 'chartreuse'),

              # legend.title = element_text(size= 17, color = "black", face="bold"),
              # legend.text = element_text(colour="black", size= 17,face="bold"),
              # legend.background = element_rect(fill = alpha("white", 0.5)),
              # legend.position = c(0.11, 0.8), # legend.position = c(0.11, 0.96),
              # legend.direction= "horizontal",
              #     plot.text = element_text(size = 20),
              #     aspect.ratio=AspRat
        )

      ## Add PValue
      plt.ManyGroup2 <- plt.ManyGroup2 +
        stat_compare_means(aes(group = Cachexia),
                           label =  "p.signif",label.x = 1.5,
                           label.y = LabelY*1.2*0.9,
                           method = "wilcox.test", size = 7)

      plt.ManyGroupSum <- plt.ManyGroupSum/plt.ManyGroup2
    }

  }

  plt.ManyGroupSum %>% print()


  # # Use only p.format as label. Remove method name.
  # p + stat_compare_means(label = "p.format", method = "wilcox.test", size = 7)
  # # Or use significance symbol as label
  # p + stat_compare_means(label =  "p.signif",label.x = 1.5, label.y = LabelY*0.9, method = "wilcox.test", size = 7)


  # ## Violin
  # ## http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/
  # plt.ManyGroup2_2_Sum <- ggviolin(Anno_Sum.df, x = "celltype", y = TarGene,
  #                                   color = "Cachexia", palette = "jco",
  #                                   add = "jitter",)
  # plt.ManyGroup2_2_Sum + stat_compare_means(aes(group = Cachexia), label = "p.format")
  # plt.ManyGroup2_2_Sum + stat_compare_means(aes(group = Cachexia), label = "p.signif")

}
dev.off()

##### Export TSV #####
write.table( SummaryTable.df ,
             file = paste0(Save.Path,"/",Version,"_LR_Stats.tsv"),
             sep = "\t",
             quote = F,
             row.names = F
)

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

