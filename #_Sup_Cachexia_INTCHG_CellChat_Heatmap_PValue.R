### Ref: Add P-values and Significance Levels to ggplots
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

## https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/plotting-models.html
## https://github.com/kassambara/ggpubr/issues/111

# ##### Presetting ######
# rm(list = ls()) # Clean variable
# memory.limit(150000)


#### Installation and load the required libraries ####
#### Basic installation ####
## Package.set
Package.set <- c("tidyverse","CellChat","patchwork","NMF","ggalluvial","Seurat","ggpubr", "stringr")
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

##### Condition Setting ####
# ## INTCHG: Interchangeable
# SampleTypeSet = "PBMC"

## CellChat DB Set
CCDBType = "ECM" # c("ECM","CC","Secret")

##### Load RData* #####

## Load Seurat.Obj
# load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_Results_1stSubmission/2022-09-09_PBMC_Main/06_Cell_type_annotation.RData")
load(paste0(Save.Path,"/06_Cell_type_annotation.RData"))

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
  rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Save.Path","CCDBType","CellType.Order")))

  # ## Modify the Cachexia state name
  # scRNA.SeuObj@meta.data$Cachexia <-  gsub("EO", "EOCX", scRNA.SeuObj@meta.data$Cachexia)
  # scRNA.SeuObj@meta.data$Cachexia <-  gsub("LO", "PreCX", scRNA.SeuObj@meta.data$Cachexia)
  #


## Load CellChat rds
cellchat.EOCX <- readRDS(paste0(Save.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_EOCX_CellChat.rds"))
cellchat.PreCX <- readRDS(paste0(Save.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_PreCX_CellChat.rds"))

object.list <- list(PreCX = cellchat.PreCX, EOCX = cellchat.EOCX)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rm(object.list, cellchat.EOCX, cellchat.PreCX)


##### Current path and new folder setting  #####
Version = paste0(Sys.Date(),"_", SampleType, "_", CCDBType, "_CellChat_PVal")
SaveCC.Path = paste0(Save.Path,"/",Version)
dir.create(SaveCC.Path)



##### Pathway and TarGene Setting  #####
pathways.show1 <-cellchat@netP[["EOCX"]][["pathways"]]
pathways.show2 <-cellchat@netP[["PreCX"]][["pathways"]]
pathways.show <- unique(pathways.show1,pathways.show2)
rm(pathways.show1,pathways.show2)

TarGene_All <- cellchat@data.signaling %>% rownames
LR.df <- rbind(cellchat@LR[["EOCX"]][["LRsig"]],cellchat@LR[["PreCX"]][["LRsig"]])

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


## Clean up data (Remove Other)
Anno.df <- Anno.df[!grepl("Other", Anno.df$celltype),]
GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df$ID]
scRNA.SeuObj <- scRNA.SeuObj[,!grepl("Other", scRNA.SeuObj@meta.data[["celltype"]] )]


##### Summarize all signal #####
SummaryTable.df <-  as.data.frame(matrix(nrow=0,ncol=10))
colnames(SummaryTable.df) <- c( "celltype", ".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method","pathway_name")


for (j in 1:length(pathways.show)) {
try({

  LR_Tar.df <- LR.df[LR.df$pathway_name == pathways.show[j],]

  # ## Method1
  # # TarGene <- c("Vwf","Itga2b","Itgb3","Gp9")
  # TarGene <- c(LR_Tar.df$ligand, LR_Tar.df$receptor) %>% unique()
  # # TarGene <- c(LR_Tar.df$ligand[1], LR_Tar.df$receptor) %>% unique()
  # TarGene <-intersect(TarGene,row.names(GeneExp.df))

  ## Method2
  library(stringr)
  library(Hmisc)
  TarGene <- LR_Tar.df$interaction_name %>%
             str_split(pattern = "_", n = Inf, simplify = FALSE) %>%
             unlist() %>%
             unique() %>% tolower() %>% capitalize()
  TarGeneH <- LR_Tar.df$interaction_name %>%
    str_split(pattern = "_", n = Inf, simplify = FALSE) %>%
    unlist() %>%
    unique()
  TarGeneH <- intersect(TarGeneH,row.names(GeneExp.df))
  TarGene <-intersect(TarGene,row.names(GeneExp.df))

  source("FUN_HSsymbol2MMsymbol.R")
  df <- TarGene %>% as.data.frame()
  colnames(df) <- "Gene"

  df1 <- HSsymbol2MMsymbol(df,"Gene")
  TarGeneM <- df1$MM.symbol
  TarGene <- c(TarGeneH,TarGeneM) %>% unique()
  rm(LR_Tar.df, df, df1)


  ##### Data preprocessing #####
  ## Extract Target gene and combine to the annotation table
  TarGene_Temp.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene,] %>% t() %>% as.data.frame()

  TarGene_Temp.df <- data.frame(ID = row.names(TarGene_Temp.df), TarGene_Temp.df)
  Anno_Temp.df <- left_join(Anno.df,TarGene_Temp.df)
  # rm(TarGene_Temp.df)


  # ## Summary Statistic Table
  # #(Ori)# SummaryTable.df <- compare_means( Vwf ~ Cachexia, data = Anno_Temp.df, group.by = "celltype"	)
  #
  # # ## Error (Solved)
  # # TTT <- compare_means( Anno_Temp.df[,TarGene[1]] ~ Cachexia, data = Anno_Temp.df, group.by = "celltype"	)
  #
  # # https://stackoverflow.com/questions/44776446/compare-means-must-resolve-to-integer-column-positions-not-a-symbol-when-u
  # # convert string column name to name/symbol
  # f <- paste0(TarGene[1]," ~ Cachexia") # f <- "Vwf ~ Cachexia"
  # SummaryTable.df <- do.call("compare_means", list(as.formula(f), data=Anno_Temp.df, group.by = "celltype"))
  # rm(f)
  # SummaryTable.df$celltype <- factor(SummaryTable.df$celltype  ,levels = CellType.Order)
  # SummaryTable.df <- SummaryTable.df[order(SummaryTable.df$celltype), ]

  ##### Summary Statistic Table #####

  SummaryTable_Sub.df <-  as.data.frame(matrix(nrow=0,ncol=9))
  colnames(SummaryTable_Sub.df) <- c( "celltype", ".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method"  )
  for (i in 1:length(TarGene)) {

    try({
      # https://stackoverflow.com/questions/44776446/compare-means-must-resolve-to-integer-column-positions-not-a-symbol-when-u
      # convert string column name to name/symbol
      f <- paste0(TarGene[i]," ~ Cachexia") # f <- "Vwf ~ Cachexia"
      SummaryTable_Temp.df <- do.call("compare_means", list(as.formula(f), data=Anno_Temp.df, group.by = "celltype"))
      rm(f)
      SummaryTable_Temp.df$celltype <- factor(SummaryTable_Temp.df$celltype  ,levels = CellType.Order)
      SummaryTable_Temp.df <- SummaryTable_Temp.df[order(SummaryTable_Temp.df$celltype), ]

      ## Filter
      if(c("****") %in% SummaryTable_Temp.df$p.signif || c("***") %in% SummaryTable_Temp.df$p.signif|| c("**") %in% SummaryTable_Temp.df$p.signif){
        SummaryTable_Sub.df <- rbind(SummaryTable_Sub.df,SummaryTable_Temp.df)
      }else{
        SummaryTable_Sub.df <- SummaryTable_Sub.df
      }

      # ## Withot Filter
      # SummaryTable_Sub.df <- rbind(SummaryTable_Sub.df,SummaryTable_Temp.df)


      rm(SummaryTable_Temp.df)



      # if(i==1){
      #   SummaryTable_Sub.df <- SummaryTable_Temp.df
      # }else{
      #   SummaryTable_Sub.df <- rbind(SummaryTable_Sub.df,SummaryTable_Temp.df)
      #   rm(SummaryTable_Temp.df)
      # }

    })
  }
  SummaryTable_Sub.df$pathway_name <- pathways.show[j]
  SummaryTable.df <- rbind(SummaryTable.df, SummaryTable_Sub.df)

  # TarGene <- SummaryTable_Sub.df$.y. %>% unique()
})

}
TarGene_Sum <- SummaryTable.df$.y. %>% unique()

## Extract Target gene and combine to the annotation table
TarGeneEXP.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene_Sum,] %>% t() %>% as.data.frame()

TarGeneEXP.df <- data.frame(ID = row.names(TarGeneEXP.df), TarGeneEXP.df)
Anno.df <- left_join(Anno.df,TarGeneEXP.df)

matrix.df <- Anno.df[,TarGene_Sum] %>% t()
colnames(matrix.df) <- Anno.df$ID

TarGeneAnno_Temp.df <- data.frame(gene=TarGene_Sum) %>% left_join(SummaryTable.df[,c("gene","pathway_name")]) %>% unique()


##### Export TSV #####
colnames(SummaryTable.df)[2] <- "gene"
SummaryTable.df <- relocate(SummaryTable.df,pathway_name,.before = gene)


write.table( SummaryTable.df ,
             file = paste0(SaveCC.Path,"/",Version,"_LR_Stats.tsv"),
             sep = "\t",
             quote = F,
             row.names = F
)

## ********************************************************************************************************************************* ##
##### ComplexHeatmap #####
## https://jokergoo.github.io/ComplexHeatmap-reference/book/

##### Load Packages #####
Package.set <- c("tidyverse","circlize","ComplexHeatmap","stringr")
## Check whether the installation of those packages is required from basic
for (i in 1:length(Package.set)) {
  if (!requireNamespace(Package.set[i], quietly = TRUE)){
    install.packages(Package.set[i])
  }
}
## Load Packages
lapply(Package.set, library, character.only = TRUE)
rm(Package.set,i)

##### Heatmap plotting #####
## Set column annotation
sample = c("#2267a4", "#8e7cc3", "#e06666", "#c27ba0")
names(sample) <- Anno.df$sample %>% unique()

library(ggsci)
library(ggplot2)
# vignette( "ggsci")
# col3 = pal_npg("nrc", alpha = 0.7)(length(CellType.Order))
col3 = 	pal_d3("category20", alpha = 0.7)(length(CellType.Order))
colCT <- col3[1:length(CellType.Order)]
names(colCT) <- CellType.Order

ha_column_T = HeatmapAnnotation(
  Sample = Anno.df[,"sample"],  # anno_colum.df$sample_type,
  Cachexia = Anno.df[,"Cachexia"], # anno_colum.df$gender,
  Celltype = Anno.df[,"celltype"],
  col = list(Sample = sample,
             Celltype = colCT ,#pal_npg(), #colCT ,
             Cachexia = c("EOCX"="#38761d", "PreCX"="#e69138")), #,"Medium"="#b57545"
             # Gender = c("Male"="#4382b5", "Female"="#c25988"), #,"Medium"="#b57545"
             # Celltype = c("High"="#db8051", "Low"="#c26334")), # #b6d4ca
  show_legend = T
)

ha_column_T2 = HeatmapAnnotation(
  Sample = Anno.df[,"sample"],  # anno_colum.df$sample_type,
  Cachexia = Anno.df[,"Cachexia"], # anno_colum.df$gender,
  Celltype = Anno.df[,"celltype"],
  col = list(Sample = sample,
             Celltype = colCT ,#pal_npg(), #colCT ,
             Cachexia = c("EOCX"="#38761d", "PreCX"="#e69138")), #,"Medium"="#b57545"
  show_legend = T,
  show_annotation_name = F
)

# # top_annotation = ha_column_T,
# top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4))),


# rm(sample)

# ## Set column annotation
# col_ha_column <- list(c("#9b6ab8", "#6e6970"),
#                       c("#4382b5", "#c25988"),
#                       c("#db8051", "#c26334"))
# names(col_ha_column) <- c(PhenoGroupType, PhenoGroupType2, TarGene_name)
#
# ha_column_Anno <- list(anno_colum.df[,PhenoGroupType],
#                        anno_colum.df[,PhenoGroupType2],
#                        anno_colum.df[,TarGene_name],
#                        col = col_ha_column %>% as.vector.factor(),
#                        show_legend = T)
# names(ha_column_Anno)[1:3] <- c(PhenoGroupType, PhenoGroupType2, TarGene_name)
#
# formals(HeatmapAnnotation)[names(ha_column_Anno)] <- ha_column_Anno
# formals(HeatmapAnnotation)[names(ha_column_Anno)] <- ha_column_Anno
# ha_column_T = HeatmapAnnotation()
#
# rm(col_ha_column, ha_column_Anno)

## Set row annotation
## Color setting
Pathway.set <- TarGeneAnno_Temp.df[,"pathway_name"] %>% unique()
col3 = pal_npg("nrc", alpha = 0.7)(length(Pathway.set))
colPT <- col3[1:length(Pathway.set)]
names(colPT) <- Pathway.set

TarGeneAnno_Temp2.df <- TarGeneAnno_Temp.df %>% group_by(gene) %>% count(pathway_name)
TarGeneAnno.df <- TarGeneAnno_Temp2.df %>% pivot_wider(names_from ="pathway_name" ,values_from = "n", values_fill =0) %>% t()
colnames(TarGeneAnno.df) <- TarGeneAnno.df[1,]
TarGeneAnno.df <- TarGeneAnno.df[-1,]
TarGeneAnno.df <- gsub("0","F",TarGeneAnno.df)
TarGeneAnno.df <- gsub("1","T",TarGeneAnno.df)

rownames(TarGeneAnno.df)[2]

colPT <- c("#45818e","#d0e0e3")
names(colPT) <- c("T","F")

ha_row = rowAnnotation(
  Pathway = TarGeneAnno.df[row.names(TarGeneAnno.df)[2],], # %>% as.character(),  # anno_colum.df$sample_type,
  col = list(#Sample = sample,
             Pathway = colPT
              ),
  show_legend = T,
  show_annotation_name = T
)



## Plot Heatmap

# Set Heatmap color
col_HMap <- c("#416db0", "#1a2938", "#bf627e")

# # Heatmap without clustering
# Heatmap(
#   matrix.df,
#   cluster_rows = F,
#   cluster_columns = F,
#   show_column_names = F,
#   show_row_names = F,
#   name = "GeneExp",
#   # set color
#   col = colorRamp2(
#     c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
#     col_HMap
#   ),
#   show_heatmap_legend = T,
#   use_raster = F,
#   top_annotation = ha_column_T,
#   # right_annotation = ha_row
# ) -> P.Heatmap1
#
# P.Heatmap1 %>% print

# # Heatmap with clustering
# Heatmap(
#   matrix.df,
#   # column_title = target_gene,
#   # column_title_side = "top",
#   cluster_rows = T,
#   cluster_columns = T,
#   show_column_names = F,
#   show_row_names = F,
#   name = "GeneExp",
#   # set color
#   col = colorRamp2(
#     c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
#     col_HMap
#   ),
#   show_heatmap_legend = T,
#   use_raster = F,
#   top_annotation = ha_column_T,
#   # right_annotation = ha_row
# ) -> P.Heatmap2
#
# P.Heatmap2 %>% print

# Reorder Heatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#row-and_column_orders
Heatmap(
  matrix.df,
  cluster_rows = F,
  cluster_columns = F,
  column_order = order(Anno.df$celltype,Anno.df$Cachexia),
  show_column_names = F,
  show_row_names = T,
  name = "GeneExp",
  # set color
  col = colorRamp2(
    c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
    col_HMap
  ),
  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T2,
  # right_annotation = ha_row
) -> P.Heatmap3

P.Heatmap3 %>% print


# ## Block annotation
# split = rep(1:le, each = 10)
#
# Heatmap(
#   matrix.df,
#   cluster_rows = T,
#   cluster_columns = F,
#
#   column_order = order(Anno.df$celltype,Anno.df$Cachexia),
#   show_column_names = F,
#   show_row_names = T,
#   name = "GeneExp",
#   # set color
#   col = colorRamp2(
#     c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
#     col_HMap
#   ),
#   show_heatmap_legend = T,
#   use_raster = F,
#   # top_annotation = ha_column_T,
#   top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4))),
#   # right_annotation = ha_row
#   column_split = Anno.df$celltype,
#
# ) -> P.Heatmap4
#
# P.Heatmap4 %>% print

# Test Combine Fig
P.Heatmap3+P.Heatmap3
# P.Heatmap3+P.Heatmap4


TarGeneAnnoMax.df <- TarGeneAnno.df %>% t() # %>% as.data.frame()
TarGeneAnnoMax.df <- TarGeneAnnoMax.df[row.names(matrix.df), ,drop=F]

colPT <- c("#45818e","#d0e0e3")
names(colPT) <- c("T","F")
col_HMap2 <- c("#45818e","#d0e0e3")
Heatmap(
  TarGeneAnnoMax.df ,
  cluster_rows = F,
  cluster_columns = F,
  # column_order = order(Anno.df$celltype,Anno.df$Cachexia),
  show_column_names = T,
  show_row_names = T,
  name = "Pathway",
  # set color
  col = col_HMap2,
  show_heatmap_legend = T,
  use_raster = F,
  # top_annotation = ha_column_T2,
  # right_annotation = ha_row
  width = ncol(TarGeneAnnoMax.df)*unit(5, "mm"),
  height = nrow(TarGeneAnnoMax.df)*unit(5, "mm")
) -> P.Heatmap3TTT

P.Heatmap3TTT %>% print

P.Heatmap3 + P.Heatmap3TTT

## ********************************************************************************************************************************* ##





## Set row annotation
## Color setting
col_exp <-  colorRamp2(
  c(min(anno_row.df$PValue), mean(anno_row.df$PValue), max(anno_row.df$PValue)),
  c("#3f705a", "#52bf8e","#b6d4ca")

)
col_exp2 <-  colorRamp2(
  c(min(anno_row.df$logFC), mean(anno_row.df$logFC), max(anno_row.df$logFC)),
  c("#488c67", "#333333","#edd493")
)

ha_row = rowAnnotation(
  p.value = anno_row.df$PValue,
  LogFC = anno_row.df$logFC,
  col = list(p.value = col_exp, LogFC = col_exp2),
  show_legend = T
)











##### Save RData #####
save.image(paste0(SaveCC.Path,"/",Version,"_LR_Stats_Heatmap.RData"))


