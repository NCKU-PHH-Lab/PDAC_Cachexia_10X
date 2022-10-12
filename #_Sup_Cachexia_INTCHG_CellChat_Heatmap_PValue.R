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

if(SampleTypeSet == "PBMC"){
  ## For PBMC
  scRNA.SeuObj <- PBMC.combined
  SampleType = "PBMC"

  # Order the cell type
  CellType.Order = c("Mac1", "Mac2","Mac3","Neu","T","CD4+T","CD8+T","NK","B","Mast","Ery")
  scRNA.SeuObj@meta.data[["celltype"]] <- factor(scRNA.SeuObj@meta.data[["celltype"]] ,
                                                 levels = CellType.Order)


}else if(SampleTypeSet == "SC"){
  ## For SC
  scRNA.SeuObj <- SC.combined
  SampleType = "SC"

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
  TarGene.df <- GeneExp.df[row.names(GeneExp.df) %in% TarGene,] %>% t() %>% as.data.frame()

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

    try({
      # https://stackoverflow.com/questions/44776446/compare-means-must-resolve-to-integer-column-positions-not-a-symbol-when-u
      # convert string column name to name/symbol
      f <- paste0(TarGene[i]," ~ Cachexia") # f <- "Vwf ~ Cachexia"
      SummaryTable_Temp.df <- do.call("compare_means", list(as.formula(f), data=Anno.df, group.by = "celltype"))
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

  TarGene <- SummaryTable_Sub.df$.y. %>% unique()
  TarGene_Sum <- c(TarGene_Sum,TarGene)
})

}
TarGene_Sum <- SummaryTable.df$.y. %>% unique()

##### Export TSV #####
colnames(SummaryTable.df)[2] <- "gene"
SummaryTable.df <- relocate(SummaryTable.df,pathway_name,.before = gene)


write.table( SummaryTable.df ,
             file = paste0(SaveCC.Path,"/",Version,"_LR_Stats.tsv"),
             sep = "\t",
             quote = F,
             row.names = F
)

##### Save RData #####
save.image(paste0(SaveCC.Path,"/",Version,"_LR_Stats_Heatmap.RData"))


