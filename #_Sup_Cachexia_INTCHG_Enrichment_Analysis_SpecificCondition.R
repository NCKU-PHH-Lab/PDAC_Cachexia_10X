## Tutorial:Gene Set Enrichment Analysis (fgsea)
## Ref: https://www.biostars.org/p/467197/
## RNA-seq analysis in R
## Ref: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html


##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

  Save.Path <- c("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_PBMC_Main")
  SampleType = "PBMC"

##### Load Packages #####
  # if(!require("tidyverse")) install.packages("tidyverse")
  # library(tidyverse)

  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


  #### BiocManager installation ####
  ## Set the desired organism
  # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db")   ##  c("org.Dm.eg.db")

  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("clusterProfiler","enrichplot","pathview") # c(organism,"fgsea","clusterProfiler","enrichplot","pathview")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  # options(stringsAsFactors = FALSE)

  # Sys.setlocale(category = "LC_ALL", locale = "UTF-8")

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")

##### Load RData* #####
  load(paste0(Save.Path,"/09_1_GSEA_Analysis_(SPA).RData"))

  ## INTCHG: Interchangeable
  ## SubType Setting
    if(SampleType == "PBMC"){
      # For PBMC
      scRNA.SeuObj <- PBMC.combined

    }else if(SampleType == "SC"){
      # For SC
      scRNA.SeuObj <- SC.combined

    }


  #### Clean up Object ####
  rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Save.Path",str_subset(objects(), pattern = "GSEA"))))
  # rm(list=setdiff(ls(), str_subset(objects(), pattern = "Venn")))

  ## Save Ori
  scRNA_Ori.SeuObj <- scRNA.SeuObj

  ## Clean up data (Delete other type)
  scRNA.SeuObj <- scRNA.SeuObj[,!grepl("Other", scRNA.SeuObj@meta.data[["celltype"]] )]

  ## Clean up data (Delete other type)

#####***************************************************************************#####
#####*  Plot by previous results *#####
##### Set condition #####
  SubType = "Neu"
  Set_FDR <- 0.05
  Set_NES <- 1

##### Extract df #####

  GSEA_Sub.df <- GSEA_Large.df.TOP[GSEA_Large.df.TOP$PhenoType %in% SubType,]
  # PBMC.combined$celltype <- factor(PBMC.combined$celltype,levels = unique(PBMC.combined$celltype))
  # GSEA_Sub.df$NES <- factor(GSEA_Sub.df$NES)

  ## Filter by FDR & NES
  GSEA_Sub.df <- GSEA_Sub.df[GSEA_Sub.df$padj <= Set_FDR & abs(GSEA_Sub.df$NES) > 1,]

##### Plot #####
  NumGenesetsPlt=15
  Barplot <- ggplot(GSEA_Sub.df, aes(NES, fct_reorder(pathway, NES), fill = padj), showCategory=(NumGenesetsPlt*2)) +
    geom_barh(stat='identity') +
    scale_fill_continuous(low = "#d45772", high = "#3b74bf", guide = guide_colorbar(reverse=TRUE)) +
    theme_minimal() + ylab(NULL)

  Barplot <- Barplot %>% BeautifyggPlot(LegPos = c(0.9, 0.15), AxisTitleSize=1.7, YtextSize=11,OL_Thick = 1.5)
  Barplot

#####***************************************************************************#####
#####*  Rerun GSEA *#####
##### Load Package #####
  library(DESeq2)
  library(org.Hs.eg.db)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(ggplot2)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_ggPlot_vline.R")
  source("FUN_GSEA_ANAL.R")
  source("FUN_DistrPlot.R")

##### Import setting and Import #####
  ## Import GSEA gene sets
  # InputGSEA <- "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"
  # InputGSEA <- "m5_go_bp_v0_3_symbols.gmt"  # InputGSEA <- "m2.all.v0.3.symbols.gmt"

  InputGSEA <- "m2.all.v0.3.symbols.gmt"
  InFOLName_GSEA <- "Input_Genesets"
  Pathway.all <- read.delim2(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA),
                             col.names = 1:max(count.fields(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA))),
                             header = F,sep = "\t")

##### Conditions setting* #####
  Group_Mode <- "GoupByPheno"   # c("GoupByPheno","GoupByGeneExp")
  TarGene_name <- "Chil3"
  PhenoGrp_name1 <- "Cachexia"
  PhenoGrp_name2 <- c("PreCX","EOCX")

  GeneExpSet.lt <- list(GeneExpMode = "Mean", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                        UpCutoff = 1, LowerCutoff = 1)

  if(Group_Mode == "GoupByGeneExp"){
    ## Group by GeneExp
    AnnoSet.lt <- list(GroupType = TarGene_name, GroupCompare = c("Low","High") )   ## DEG by GeneExp group
  }else{
    ## Group by Pheno
    AnnoSet.lt <- list(GroupType = PhenoGrp_name1, GroupCompare = PhenoGrp_name2 )
  }

  Thr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )

##### Current path and new folder setting* #####
  ProjectName = "CC10X"
  # SampleType = "PBMC"

  if(Group_Mode == "GoupByGeneExp"){
    ExportAnno = paste0(TarGene_name,"_",GeneExpSet.lt$GeneExpMode,"_",SubType)

  }else{
    ExportAnno = paste0(Group_Mode,"_",paste0(PhenoGrp_name2[1],PhenoGrp_name2[2],"_",SubType))

  }

  # ExportAnno = "Chil3Mean_PathM2"
  # ExportAnno = "Recur2Prim"

  ExportName = paste0(ProjectName,"_",SampleType,"_",ExportAnno)


  Version = paste0(Sys.Date(),"_",ProjectName,"_",SampleType,"_", ExportAnno)
  SaveSub.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(SaveSub.Path)){
    dir.create(SaveSub.Path)
  }

##### Update the genename ####
  ## Update the genename ##* Take very long time
  UpdateGene <- "No"  # UpdateGene <- c("Yes","No")
  if(UpdateGene == "Yes"){
    row.names(GeneExp.df) <- UpdateSymbolList(row.names(GeneExp.df))
  }

#************************************************************************************************************************#
##### Data preprocess setting #####
  ## Extract data from scRNA.SeuObj
  ## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

  Anno.df <- scRNA.SeuObj@meta.data
  Anno.df <- data.frame(ID=row.names(Anno.df), Anno.df)
  Anno.df <- left_join(data.frame("ID"=colnames(GeneExp.df)),
                       Anno.df)
  row.names(Anno.df) <- Anno.df[,1]

  ## Select Pheno column
  Anno_Ori.df <- Anno.df
  colnames(Anno.df)

  # PhenoColKeep.set <- c("X_INTEGRATION","X_PATIENT","histological_type","sample_type","gender")
  # Anno.df <- Anno.df[,c(PhenoColKeep.set)]
  # colnames(Anno.df)
  #
  # head(Anno.df)

  # ## Select Pheno row
  # PhenoRowKeep.set <- list(col="Cachexia" ,row=c("EOCX"))
  # Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]
  #
  # GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df[,1] ]
  # rm(PhenoRowKeep.set)

  ## Select Pheno row2
  PhenoRowKeep.set <- list(col="celltype" ,row=c(SubType))
  Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]

  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df[,1] ]
  rm(PhenoRowKeep.set)

  # ## Delete specific cell type
  # ## Clean up data
  # Anno.df <- Anno.df[!grepl("Other", Anno.df$celltype),]
  # GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df$ID]



#************************************************************************************************************************#
##### Visualization #####
  source("FUN_DistrPlot.R")
  ##### Group by gene expression 1: CutOff by total  #####
  Plot.DistrPlot <- FUN_DistrPlot(GeneExp.df,
                                  TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                  Save.Path = SaveSub.Path, ExportName = ExportName)
  Plot.DistrPlot_SD_Q <- Plot.DistrPlot[["TGeneDen_SD_Q.p"]]
  Plot.DistrPlot_SD_Q


#************************************************************************************************************************#
##### Grouping #####
  source("FUN_Group_GE.R")
  ##### Group by gene expression 1: CutOff by total  #####
  GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Anno.df,
                                    TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                    Save.Path = SaveSub.Path, ExportName = ExportName)
  Anno.df <- GeneExp_group.set[["AnnoNew.df"]]
  GeneExp_high.set <- GeneExp_group.set[["GeneExp_high.set"]]
  GeneExp_low.set <- GeneExp_group.set[["GeneExp_low.set"]]

  ##### Group by gene expression 2: CutOff by Comparison #####
  ## FUN Comparison (Visualization and value)

  ##### Group by phenotype #####


#************************************************************************************************************************#
##### Run Enrichment analysis in R #####
  #### Run DEG ####
  source("FUN_DEG_Analysis.R")
  DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                              GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
                              ThrSet = Thr.lt,
                              TarGeneName = TarGene_name, GroupMode = GeneExpSet.lt, SampleID = "ID",
                              Save.Path = SaveSub.Path, ExportName = ExportName, AnnoName = "AvB")
  DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]


  #### Run GSEA ####
  source("FUN_GSEA_ANAL.R")

  GSEA_Result.lt <- FUN_GSEA_ANAL(DE_Extract.df, CMGeneSet = Pathway.all,
                                  NumGenesetsPlt=15,
                                  TarGeneName = TarGene_name,
                                  ThrSet = Thr.lt, Species = "Homo sapiens", # Speices type can check by msigdbr_species()
                                  Save.Path = SaveSub.Path, ExportName = ExportName, AnnoName = "Path")





  #### Run ORA ####
  ## FUN ORA

#************************************************************************************************************************#
##### Build files for GSEA official input #####
  source("FUN_GSEA_ForOFFL.R")

  FUN_GSEA_ForOFFL(GeneExp.df, Group1 = GeneExp_high.set, Group2 = GeneExp_low.set,
                   GroupMode = Group_Mode,
                   TarGeneName = TarGene_name, GeneExpSet = GeneExpSet.lt,
                   Save.Path = SaveSub.Path, ExportName = ExportName,
                   AnnoName="Recur2Prim")

##### Build files for Metascape official input #####



# #### Save RData ####
#   save.image(paste0(Save.Path,"/GseaGo_",ExportName,".RData"))




