### ComplexHeatmap
### Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html

### Add P-values and Significance Levels to ggplots
### Ref: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

## https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/plotting-models.html
## https://github.com/kassambara/ggpubr/issues/111

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

#### Installation and load the required libraries ####
#### Basic installation ####
## Package.set
Package.set <- c("tidyverse","CellChat","patchwork","NMF","ggalluvial","Seurat","ggpubr", "stringr", "Hmisc",
                 "circlize","ComplexHeatmap")
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

## Set Path
Save.Path <- c("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_PBMC_Main")

## Set CellChat DB
CCDBType = "ECM" # c("ECM","CC","Secret")

##### Load Data* #####
## Load RData
# load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_Results_1stSubmission/2022-09-09_PBMC_Main/06_Cell_type_annotation.RData")
load(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))

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
  rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Save.Path","CCDBType","CellType.Order",
                          "CCMarker_Female.lt","CCMarker_Male.lt","CCMarker_SPA.lt")))

## Load CellChat rds
cellchat.EOCX <- readRDS(paste0(Save.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_EOCX_CellChat.rds"))
cellchat.PreCX <- readRDS(paste0(Save.Path,"/",SampleType,"_CellCell_Interaction/",CCDBType,"_PreCX_CellChat.rds"))

object.list <- list(PreCX = cellchat.PreCX, EOCX = cellchat.EOCX)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rm(object.list, cellchat.EOCX, cellchat.PreCX)


##### Current path and new folder setting*  #####
Version = paste0(Sys.Date(),"_", SampleType, "_", CCDBType, "_CellChat_Heatmap")
SaveCC.Path = paste0(Save.Path,"/",Version)
dir.create(SaveCC.Path)

##### Pathway and TarGene Setting  #####
Pathways_Sig.set <- unique(cellchat@netP[["EOCX"]][["pathways"]],
                           cellchat@netP[["PreCX"]][["pathways"]])

TarGene_Ori.set <- cellchat@data.signaling %>% rownames
LR.df <- rbind(cellchat@LR[["EOCX"]][["LRsig"]],cellchat@LR[["PreCX"]][["LRsig"]]) %>% unique()


##### Extract information #####
## Gene Expression
## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix
## Set y position
Colr_GE.set <- c(min(GeneExp.df), max(GeneExp.df) %>% ceiling())

## Annotation
# Cell annotation
Anno_Cell.df <- scRNA.SeuObj@meta.data
Anno_Cell.df <- data.frame(ID=row.names(Anno_Cell.df), Anno_Cell.df)

# Gene annotation
Ligand.set <- LR.df$ligand
Receptor.set <- LR.df$receptor %>%
                str_split(pattern = "_", n = Inf, simplify = FALSE) %>%
                unlist() %>% str_split(pattern = " ", n = Inf, simplify = FALSE) %>% unlist() %>%
                unique() %>% tolower() %>% capitalize()
# Receptor.set <- setdiff(Receptor.set, Ligand.set) %>% unlist()
Anno_TarGene.df <- data.frame(gene = TarGene_Ori.set, Ligand ="", Receptor ="")
for (i in 1:nrow(Anno_TarGene.df)) {
  if(Anno_TarGene.df[i,"gene"] %in% Ligand.set){
    Anno_TarGene.df[i,"Ligand"] = "T"
  }else{
    Anno_TarGene.df[i,"Ligand"] = "F"
  }
  if(Anno_TarGene.df[i,"gene"] %in% Receptor.set){
    Anno_TarGene.df[i,"Receptor"] = "T"
  }else{
    Anno_TarGene.df[i,"Receptor"] = "F"
  }
}
rm(i)

## Save Ori
OriSave.lt <- list(GeneExp.df = GeneExp.df,
               Anno_Cell.df = Anno_Cell.df,
               scRNA.SeuObj = scRNA.SeuObj,
               Anno_TarGene.df = Anno_TarGene.df)


##### Summarize all signal #####
SummaryTable.df <-  as.data.frame(matrix(nrow=0,ncol=10))
colnames(SummaryTable.df) <- c( "celltype", ".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method","pathway_name")

for (j in 1:length(Pathways_Sig.set)) {
try({

  #### Extract LR genes ####
  LR_Tar.df <- LR.df[LR.df$pathway_name == Pathways_Sig.set[j],]

  library(stringr)
  library(Hmisc)
  TarGene <- LR_Tar.df$interaction_name %>%
             str_split(pattern = "_", n = Inf, simplify = FALSE) %>%
             unlist() %>% unique() %>% tolower() %>% capitalize()
  ## Human Genes
  TarGeneH.set <- LR_Tar.df$interaction_name %>%
              str_split(pattern = "_", n = Inf, simplify = FALSE) %>%
              unlist() %>% unique()
  TarGeneH.set <- intersect(TarGeneH.set,row.names(GeneExp.df))
  TarGene <-intersect(TarGene,row.names(GeneExp.df))

  ## Mouse Genes
  source("FUN_HSsymbol2MMsymbol.R")
  df <- TarGene %>% as.data.frame()
  colnames(df) <- "Gene"

  df1 <- HSsymbol2MMsymbol(df,"Gene")
  TarGeneM.set <- df1$MM.symbol
  TarGene <- c(TarGeneH.set,TarGeneM.set) %>% unique()
  rm(LR_Tar.df, df, df1, TarGeneM.set, TarGeneH.set)


  ##### Data preprocessing #####
  ## Extract Target gene and combine to the annotation table
  TarGene_Temp.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene,] %>% t() %>% as.data.frame()

  TarGene_Temp.df <- data.frame(ID = row.names(TarGene_Temp.df), TarGene_Temp.df)
  Anno_Temp.df <- left_join(Anno_Cell.df,TarGene_Temp.df)

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
      # if(c("****") %in% SummaryTable_Temp.df$p.signif || c("***") %in% SummaryTable_Temp.df$p.signif|| c("**") %in% SummaryTable_Temp.df$p.signif){
      # if(c("****") %in% SummaryTable_Temp.df$p.signif){
      if(sum(SummaryTable_Temp.df$p.adj < 1.0e-7)>=1){
        SummaryTable_Sub.df <- rbind(SummaryTable_Sub.df,SummaryTable_Temp.df)
      }else{
        SummaryTable_Sub.df <- SummaryTable_Sub.df
      }

      # ## Withot Filter
      # SummaryTable_Sub.df <- rbind(SummaryTable_Sub.df,SummaryTable_Temp.df)


      rm(SummaryTable_Temp.df)

    })
  }
  SummaryTable_Sub.df$pathway_name <- Pathways_Sig.set[j]
  SummaryTable.df <- rbind(SummaryTable.df, SummaryTable_Sub.df)

  # TarGene <- SummaryTable_Sub.df$.y. %>% unique()
})

}

rm(TarGene_Temp.df, Anno_Temp.df,SummaryTable_Sub.df,i,j)



## Extract Target gene and combine to the annotation table
colnames(SummaryTable.df)[2] <- "gene"
SummaryTable.df <- relocate(SummaryTable.df,pathway_name,.before = gene)

TarGene_Sum.set <- SummaryTable.df$gene %>% unique()
TarGeneEXP.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene_Sum.set,] %>% t() %>% as.data.frame()

TarGeneEXP.df <- data.frame(ID = row.names(TarGeneEXP.df), TarGeneEXP.df)
Anno_Cell.df <- left_join(Anno_Cell.df,TarGeneEXP.df)

matrix.df <- Anno_Cell.df[,TarGene_Sum.set] %>% t()
colnames(matrix.df) <- Anno_Cell.df$ID

## TarGene with PathAnno Matrix
TarGeneAnno_Temp.df <- data.frame(gene=TarGene_Sum.set) %>% left_join(SummaryTable.df[,c("gene","pathway_name")]) %>% unique()
Pathway.set <- TarGeneAnno_Temp.df[,"pathway_name"] %>% unique()

# library(ggsci)
# col3 = pal_npg("nrc", alpha = 0.7)(length(Pathway.set))
# colPT <- col3[1:length(Pathway.set)]
# names(colPT) <- Pathway.set

TarGeneAnno_Temp2.df <- TarGeneAnno_Temp.df %>% group_by(gene) %>% count(pathway_name)
TarGeneAnno_Cell.df <- TarGeneAnno_Temp2.df %>% pivot_wider(names_from ="pathway_name" ,values_from = "n", values_fill =0) %>% t()
colnames(TarGeneAnno_Cell.df) <- TarGeneAnno_Cell.df[1,]
TarGeneAnno_Cell.df <- TarGeneAnno_Cell.df[-1,]
TarGeneAnno_Cell.df <- gsub("0","F",TarGeneAnno_Cell.df)
TarGeneAnno_Cell.df <- gsub("1","T",TarGeneAnno_Cell.df)
rm(TarGeneAnno_Temp.df, TarGeneAnno_Temp2.df)

TarGeneAnno.mtx <- TarGeneAnno_Cell.df %>% t() # %>% as.data.frame()
TarGeneAnno.mtx <- TarGeneAnno.mtx[row.names(matrix.df), ,drop=F]


##### Export TSV #####
write.table( SummaryTable.df ,
             file = paste0(SaveCC.Path,"/",Version,"_LR_Stats.tsv"),
             sep = "\t",
             quote = F,
             row.names = F
)



## ********************************************************************************************************************************* ##
##### ComplexHeatmap #####
## https://jokergoo.github.io/ComplexHeatmap-reference/book/
library(ComplexHeatmap)
library(circlize)

#### Set column annotation ####
sample = c("#2267a4", "#3d85c6", "#d5a6bd", "#c27ba0")
names(sample) <- Anno_Cell.df$sample %>% unique()

library(ggsci)
library(ggplot2)
# vignette( "ggsci") #Check the color setting
col3 = 	pal_d3("category20", alpha = 0.7)(length(CellType.Order))
colCT <- col3[1:length(CellType.Order)]
names(colCT) <- CellType.Order

ha_column_T = HeatmapAnnotation(
  Sample = Anno_Cell.df[,"sample"],  # anno_colum.df$sample_type,
  Cachexia = Anno_Cell.df[,"Cachexia"], # anno_colum.df$gender,
  Celltype = Anno_Cell.df[,"celltype"],
  col = list(Sample = sample,
             Celltype = colCT ,#pal_npg(), #colCT ,
             Cachexia = c("EOCX"="#5b517d", "PreCX"="#a095c7")), #,"Medium"="#b57545"
  show_legend = T,
  show_annotation_name = F
)




#### Plot Heatmap ####
# Set Heatmap color
col_HMap <- c("#f0f4fc", "#6e8cc2", "#37558c")
col_HMap <- c( "#ffffff","#e697b3", "#c45e82")


### Plat GeneExpression Heatmap
## Reorder Heatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#row-and_column_orders
Heatmap(
  matrix.df,
  cluster_rows = F, # Heatmap with/without clustering by rows
  cluster_columns = F, # Heatmap with/without clustering by columns
  column_order = order(Anno_Cell.df$celltype,Anno_Cell.df$Cachexia), ## Reorder Heatmap
  show_column_names = F,
  show_row_names = T,
  name = "GeneExp",
  # set color
  col = colorRamp2(
    # c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
    c(min(GeneExp.df), max(GeneExp.df)*2/3 , max(GeneExp.df)),

    col_HMap
  ),
  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T,
  # right_annotation = ha_row
  width = length(CellType.Order)*unit(7, "mm"),
  height = length(CellType.Order)*unit(17, "mm")
) -> P.Heatmap_GeneExp

P.Heatmap_GeneExp %>% print



### Plat Heatmap: Gene belong by which Pathways

#### Set row annotation ####
## Color setting
Anno_TarGene.df <- left_join(data.frame(gene=row.names(matrix.df)),Anno_TarGene.df )
colPT <- c("#646b6b","#b6baba")
names(colPT) <- c("T","F")

ha_row = rowAnnotation(
  Ligand = Anno_TarGene.df$Ligand,
  Receptor = Anno_TarGene.df$Receptor,
  col = list(Ligand = colPT,
             Receptor=colPT
             ),
  show_legend = T,
  show_annotation_name = T
)


col_HMap2 <- c("#d0e0e3","#45818e")
Heatmap(
  TarGeneAnno.mtx ,
  cluster_rows = F,
  cluster_columns = F,
  # column_order = order(Anno_Cell.df$celltype,Anno_Cell.df$Cachexia),
  show_column_names = T,
  show_row_names = T,
  name = "Pathway",
  # set color
  col = col_HMap2,
  show_heatmap_legend = T,
  use_raster = F,
  # top_annotation = ha_column_T,
  right_annotation = ha_row,
  width = ncol(TarGeneAnno.mtx)*unit(8, "mm"),
  height = length(CellType.Order)*unit(17, "mm") #nrow(TarGeneAnno.mtx)*unit(10, "mm")
) -> P.Heatmap_GenePath

P.Heatmap_GenePath %>% print

### Plat Heatmap: LogFC and FDR

## Create dataframe
for (i in 1:length(CCMarker_SPA.lt)) {
  statistics_Temp.df <- CCMarker_SPA.lt[[i]][["CCMarker.All"]]
  colnames(statistics_Temp.df)[1] <- "gene"
  statistics_Temp.df$celltype <- names(CCMarker_SPA.lt)[i]
  statistics_Temp.df <- relocate(statistics_Temp.df,celltype,.before = gene)

  if(i==1){
    statistics.df <- statistics_Temp.df
  }else{
    statistics.df <- rbind(statistics.df, statistics_Temp.df)
  }
  rm(statistics_Temp.df)
}


### FDR df Setting
statistics_FDR.df <- statistics.df[,c(1,2,7)] %>% pivot_wider(names_from = celltype, values_from = c("p_val_adj"))
statistics_FDR.mtx <- left_join(data.frame(gene=row.names(matrix.df)),statistics_FDR.df)
row.names(statistics_FDR.mtx) <- statistics_FDR.mtx[,1]
statistics_FDR.mtx <- statistics_FDR.mtx[,-1]
## Reoder the df
statistics_FDR.mtx <- statistics_FDR.mtx %>% select(CellType.Order) %>% as.matrix()

### LogFC df Setting
statistics_LogFC.df <- statistics.df[,c(1,2,4)] %>% pivot_wider(names_from = celltype, values_from = c("avg_log2FC"))
statistics_LogFC.mtx <- left_join(data.frame(gene=row.names(matrix.df)),statistics_LogFC.df)
row.names(statistics_LogFC.mtx) <- statistics_LogFC.mtx[,1]
statistics_LogFC.mtx <- statistics_LogFC.mtx[,-1]
## Reoder the df
statistics_LogFC.mtx <- statistics_LogFC.mtx %>% select(CellType.Order) %>% as.matrix()


### Plot FDR Heatmap
# Set column annotation
library(ggsci)
library(ggplot2)
# vignette( "ggsci") #Check the color setting
col3 = 	pal_d3("category20", alpha = 0.7)(length(CellType.Order))
colCT <- col3[1:length(CellType.Order)]
names(colCT) <- CellType.Order

ha_column_ST = HeatmapAnnotation(
  Celltype = colnames(statistics_FDR.mtx),
  col = list( Celltype = colCT),
  show_legend = F,
  show_annotation_name = F
)


# Set Heatmap color
col_HMapST <- c("#db9569","#f5d1ba","#f7ebe4")
col_HMapST <- c("#edffef","#a4ebab","#70b577")
col_HMapST <- c("#93db9b", "#b9f0bf","#edffef")
col_HMapST <- c("#6e8cc2", "#37558c","#f0f4fc")

col_HMapST <- c( "#ad8fd9", "#e8d9ff" , "#ffffff")

Heatmap(
  statistics_FDR.mtx,
  cluster_rows = F, # Heatmap with/without clustering by rows
  cluster_columns = F, # Heatmap with/without clustering by columns
  column_order = CellType.Order, ## Reorder Heatmap
  show_column_names = T,
  show_row_names = T,
  name = "FDR",
  # set color
  col = colorRamp2(
    # c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
    c(0, 0.0000001 , 1),
    col_HMapST
  ),
  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_ST,
  # right_annotation = ha_row,
  width = length(CellType.Order)*unit(7, "mm"),
  height = length(CellType.Order)*unit(17, "mm"),

  ## R plot pch symbols : The different point shapes available in R
  ## Ref: http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  # layer_fun = function(j, i, x, y, width, height, fill)
  # {
  #   v1 = pindex(statistics_FDR.mtx, i, j)
  #   l1 = v1 < 0.0000001
  #   grid.points(x[l1], y[l1], pch = 16, size = unit(4, "mm"))
  # }

  layer_fun = function(j, i, x, y, width, height, fill)
  {
    v1 = pindex(statistics_FDR.mtx, i, j)
    v2 = pindex(statistics_LogFC.mtx, i, j)

    l1 = v1 < 0.0000001
    grid.points(x[l1], y[l1], pch = 1, size = unit(4, "mm"))
    l3 = v1 < 0.0000001 & abs(v2) > 0.5
    grid.points(x[l3], y[l3], pch = 16, size = unit(4, "mm"))

  }


) -> P.Heatmap_FDR

P.Heatmap_FDR %>% print

### Plot LogFC Heatmap
# Set Heatmap color
col_HMapLog <- c("#4a5aa8","#ffffff","#c270b0")
Heatmap(
  statistics_LogFC.mtx,
  cluster_rows = F, # Heatmap with/without clustering by rows
  cluster_columns = F, # Heatmap with/without clustering by columns
  column_order = CellType.Order, ## Reorder Heatmap
  show_column_names = T,
  show_row_names = T,
  name = "LogFC",
  # set color
  col = colorRamp2(
    # c(min(matrix.df), matrix.df %>% unlist() %>% mean() , max(matrix.df)),
    c(min(statistics_LogFC.mtx), 0.01 , max(statistics_LogFC.mtx)),
    col_HMapLog
  ),
  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_ST,
  # right_annotation = ha_row,
  width = length(CellType.Order)*unit(7, "mm"),
  height = length(CellType.Order)*unit(17, "mm"),

  ## R plot pch symbols : The different point shapes available in R
  ## Ref: http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  layer_fun = function(j, i, x, y, width, height, fill)
  {
    v1 = pindex(statistics_FDR.mtx, i, j)
    v2 = pindex(statistics_LogFC.mtx, i, j)

    l2 = abs(v2) > 0.5
    grid.points(x[l2], y[l2], pch = 1, size = unit(4, "mm"))

    l3 = v1 < 0.0000001 & abs(v2) > 0.5
    grid.points(x[l3], y[l3], pch = 16, size = unit(4, "mm"))
  }
) -> P.Heatmap_LogFC

P.Heatmap_LogFC %>% print

#### Summary Plot ####
P.Heatmap_GeneExp + P.Heatmap_LogFC + P.Heatmap_FDR + P.Heatmap_GenePath


##### Export PDF #####
pdf(file = paste0(SaveCC.Path,"/",Version,"_LR.pdf"),width = 20, height = 20 )
P.Heatmap_GeneExp + P.Heatmap_LogFC + P.Heatmap_FDR + P.Heatmap_GenePath

dev.off()


# ##### Save RData #####
# save.image(paste0(SaveCC.Path,"/",Version,"_LR_Stats_Heatmap.RData"))


##******************************************************************************##
##### Record ####

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


#### Plot Heatmap ####

# ## Set row annotation
# ## Color setting
# col_exp <-  colorRamp2(
#   c(min(anno_row.df$PValue), mean(anno_row.df$PValue), max(anno_row.df$PValue)),
#   c("#3f705a", "#52bf8e","#b6d4ca")
#
# )
# col_exp2 <-  colorRamp2(
#   c(min(anno_row.df$logFC), mean(anno_row.df$logFC), max(anno_row.df$logFC)),
#   c("#488c67", "#333333","#edd493")
# )
#
# ha_row = rowAnnotation(
#   p.value = anno_row.df$PValue,
#   LogFC = anno_row.df$logFC,
#   col = list(p.value = col_exp, LogFC = col_exp2),
#   show_legend = T
# )
#


# ## Block annotation
# split = rep(1:le, each = 10)
#
# Heatmap(
#   matrix.df,
#   cluster_rows = T,
#   cluster_columns = F,
#
#   column_order = order(Anno_Cell.df$celltype,Anno_Cell.df$Cachexia),
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
#   column_split = Anno_Cell.df$celltype,
#
# ) -> P.Heatmap4
#
# P.Heatmap4 %>% print

# # Test Combine Fig
# P.Heatmap_GeneExp+P.Heatmap_GeneExp
# # P.Heatmap_GeneExp+P.Heatmap4

