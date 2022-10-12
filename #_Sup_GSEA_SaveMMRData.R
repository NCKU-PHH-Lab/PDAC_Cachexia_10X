##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting #####
  Save.Path = paste0(getwd(),"/20220214_PBMC")
  # dir.create(Save.Path)
  # SampleType = "PBMC"

##### Load libray #####
  library(dplyr)
  library(stringr)
  library(scales)
  library(ggplot2)

##### Function setting  #####
  ## Call function
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")

##### Load RData  #####
  load(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))


##### 09_0 GSEA Analysis (Geneset Prepare) #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(fgsea)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")

  # Geneset from GSEA
  # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
  Pathway.all <- read.delim2(paste0(getwd(),"/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"),
                             col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"))),
                             header = F,sep = "\t")

  # Convert Human gene to mouse
  Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*1.5))
  for (i in 1:nrow(Pathway.all)) {
    #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
    PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t()
    colnames(PathwayN)="Test"
    PathwayN <- HSsymbol2MMsymbol(PathwayN,"Test")
    Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)

  }

  Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
  colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))

  save.image(paste0(Save.Path,"/09_0_GSEA_Analysis(Geneset_Prepare).RData"))

  rm(list=setdiff(ls(), "Pathway.all.MM"))
  save.image(paste0(Save.Path,"/GSEA_Analysis_Geneset.RData"))

