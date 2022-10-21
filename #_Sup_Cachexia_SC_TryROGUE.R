## Ref: https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)

  Save.Path <- c("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_SC_Main")

##### Load Packages #####
  # if(!require("tidyverse")) install.packages("tidyverse")
  # library(tidyverse)

  #### Basic installation ####
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc",
                   "stringr","magrittr","dplyr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("fgsea","AnnotationHub","ensembldb",
                   "SeuratDisk","monocle",
                   "SingleR","scRNAseq","celldex","scran")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


  #### GitHub installation ####
  # if (!require("devtools", quietly = TRUE))
  #   install.packages("devtools")
  # library(monocle)
  # devtools::install_github("cole-trapnell-lab/garnett")
  # devtools::install_github('cole-trapnell-lab/monocle3')
  # devtools::install_github("LTLA/SingleR")
  #
  # library(monocle3)
  # library(garnett)
  # # library(SingleR)

  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("PaulingLiu/ROGUE")
  library(ROGUE)


##### Function setting #####
  ## Call function
  source("FUN_Cal_Mit.R")
  source("FUN_CombineSeuObj.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Anno_SingleR.R")



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

##### Current path and new folder setting* #####
  ProjectName = "Com_ROGUE"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Extract data #####
  ## Gene GeneExp.dfession
  ## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

  ## Meta.data
  # Cell annotation
  Meta.data <- scRNA.SeuObj@meta.data
  Meta.data <- data.frame(ID=row.names(Meta.data), Meta.data)


##### Run ROGUE #####
  ## Filtering out low-abundance genes and low-quality cells
  GeneExp.df <- matr.filter(GeneExp.df, min.cells = 10, min.genes = 10)

  ## GeneExp.dfession entropy model
  ent.res <- SE_fun(GeneExp.df)
  head(ent.res)

  ## S-E plot
  SEplot(ent.res)

  ## ROGUE calculation
  rogue.value <- CalculateRogue(ent.res, platform = "UMI")
  rogue.value

  ## Calculate the ROGUE value of each putative cluster for each sample
  rogue.res <- rogue(GeneExp.df, labels = Meta.data$celltype, samples = Meta.data$sample, platform = "UMI", span = 0.6)
  rogue.res

  ## Visualize ROGUE values on a boxplot
  rogue.boxplot(rogue.res)




