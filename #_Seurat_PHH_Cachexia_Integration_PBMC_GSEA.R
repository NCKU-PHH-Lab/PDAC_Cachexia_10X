memory.limit(300000)
#####  Current path and new folder setting ##### 
PathName = setwd(getwd())
RVersion = "202109018_PBMC"
dir.create(paste0(PathName,"/",RVersion))

#####  Function setting ##### 

## Call function
filePath <- ""
# Import R files in the same folder
getFilePath <- function(fileName) {
  # Absolute path of project folder
  # path <- setwd("~")
  path <- setwd(getwd()) 
  # Combine strings without gaps
  # <<- Assigning values to global variable
  filePath <<- paste0(path ,"/" , fileName)  
  # Load file
  sourceObj <- source(filePath)
  return(sourceObj)
}


getFilePath("HSsymbol2MMsymbol.R")

#####
PBMC.combined.Male <- PBMC.combined[ ,PBMC.combined@meta.data[["Sex"]] %in% c("Male")]

library(Seurat)
library(SeuratDisk)
library(SeuratObject)


PBMC.combined.Male.EO <- PBMC.combined[ ,PBMC.combined@meta.data[["Sex"]] %in% c("Male") & PBMC.combined@meta.data[["Cachexia"]] %in% c("EO")]
# library(VGAM)
PBMC.combined.Male.EO.df <- as.data.frame(PBMC.combined.Male.EO@assays[["RNA"]]@data)
# data <- as(as.matrix(PBMC.combined.Male.EO@assays[["RNA"]]@data), 'sparseMatrix')



##### RNA-seq analysis in R #####
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
library(fgsea)
gseaDat <- Cachexia.Marker.Male[["Cachexia.Marker.Male.Mac"]][["Cachexia.Marker.All"]]
gseaDat <- data.frame(row.names(gseaDat),gseaDat)
colnames(gseaDat)[[1]] <- c("Gene")
ranks <- gseaDat$avg_log2FC
names(ranks) <- gseaDat$Gene
head(ranks)

barplot(sort(ranks, decreasing = T))

# load(paste0(PathName,"/mouse_H_v5.RData"))
# pathwaysH <- Mm.H
# fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500)

# Geneset from GSEA
H.all <- read.delim(paste0(PathName,"/h.all.v7.4.symbols.gmt"),header = F)
H.all.list <- list()
for (i in c(1:length(H.all[,1]))) {
  H.all.list.ori <- as.data.frame(t(H.all[i,3:length(H.all[i,])]))
  colnames(H.all.list.ori)[[1]] <- c("Gene")
  H.all.list.ori <- HSsymbol2MMsymbol(H.all.list.ori,"Gene")
  
  # Delete NA(or 0)
  H.all.list.ori <- H.all.list.ori[H.all.list.ori$MM.symbol!=0,]
  #Bug# H.all.list.ori <- H.all.list.ori[-which(H.all.list.ori$MM.symbol==0),]
  #Bug# H.all.list.ori <- H.all.list.ori[!which(H.all.list.ori$MM.symbol==0),]
  #NoUse# H.all.list[[i]] <- na.omit(H.all.list[[i]])
  #Error# H.all.list[[i]] <- H.all.list[!is.na(H.all.list)]
  
  H.all.list.ori <- unique(H.all.list.ori$MM.symbol)
  H.all.list[[i]] <- as.character(H.all.list.ori)
  # H.all.list[[i]] <- as.character(H.all[i,3:length(H.all[i,])])  

  rm(H.all.list.ori)
  names(H.all.list)[[i]] <- H.all[i,1]
}

# load(paste0(PathName,"/Full_annotation.RData"))
pathwaysH <- H.all.list 

fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500)

head(fgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)

library(magrittr)
library(dplyr)
topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)


##### Test Cachexia candidate markers #####
pathwaysHTTT <- list(gseaDat$Gene[1:150])
pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Unique_A"]])
pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Union_AB"]])
pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Unique_B"]])
pathwaysHTTT <- list(Venn_Cachexia.Marker_Pos[["Venn_Cachexia.Marker.Mac_Pos"]][["Summary"]][["Intersect_AB"]])
names(pathwaysHTTT) <- c("Test")
fgseaRes <- fgsea(pathwaysHTTT, ranks, minSize=15, maxSize = 500)
# fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=10)

plotEnrichment(pathwaysHTTT[["Test"]], ranks)
##### Using fgsea package #####
# https://bioc.ism.ac.jp/packages/3.6/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html





# ## https://www.biostars.org/p/339934/
# library(msigdbr)
# library(fgsea)
# 
# #Retrieve human H (hallmark) gene set
# msigdbr_df <- msigdbr(species = "human", category = "H")
# head(msigdbr_df)
