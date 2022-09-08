##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Version information ######
  # platform       x86_64-w64-mingw32
  # arch           x86_64
  # os             mingw32
  # system         x86_64, mingw32
  # status
  # major          4
  # minor          1.1
  # year           2021
  # month          08
  # day            10
  # svn rev        80725
  # language       R
  # version.string R version 4.1.1 (2021-08-10)
  # nickname       Kick Things

##### Current path and new folder setting #####
  Save.Path = paste0(getwd(),"/2022-09-07_PBMC_Main")
  # dir.create(Save.Path)
  SampleType = "PBMC"

  RawData.Path = ""

##### Load libray #####
  library(dplyr)
  library(stringr)
  library(ggplot2)

##### Function setting  #####
  ## Call function
  source("FUN_CountGeneNum.R")
  source("FUN_Beautify_ggplot.R")

##### Abbreviated Notes #####
  # CCM: Cancer Cachexia Marker
  # SPA: Sex Pooled Analysis
  # SSA: Sex Separated Analysis
  # EO: Early Onset
  # LO: Late Onset
  # CT: Cell Type

  # list:lt
  # dataframe: df


##### Load data #####
  CCMList_SPA.df <- read.delim2(paste0(Save.Path,"/PBMC_CCMarker_SPA.tsv"))
  CCMList_SSA.df <- read.delim2(paste0(Save.Path,"/PBMC_CCMarker_SSA.tsv"))

  ## Modify the words
  library(stringr)
  CCMList_SPA.df$Type <- str_replace(CCMList_SPA.df$Type,"CCMarker.","")
  CCMList_SPA.df$Type <- str_replace(CCMList_SPA.df$Type,"Pos","EO")
  CCMList_SPA.df$Type <- str_replace(CCMList_SPA.df$Type,"Neg","LO")

  CCMList_SSA.df$Type <- str_replace(CCMList_SSA.df$Type,"Venn_CCMarker.","")
  CCMList_SSA.df$Type <- str_replace(CCMList_SSA.df$Type,"Pos","EO")
  CCMList_SSA.df$Type <- str_replace(CCMList_SSA.df$Type,"Neg","LO")

# ##### (Test) Specific Cell Type #####
#   CCM_SPA_Mac.df <- CCMList_SPA.df[grep("Mac", CCMList_SPA.df$Type), ]
#
#   CCM_SPA_Mac_EO.df <- CCM_SPA_Mac.df[grep("EO",CCM_SPA_Mac.df$Type), ]
#   CCM_SPA_Mac_LO.df <- CCM_SPA_Mac.df[grep("LO",CCM_SPA_Mac.df$Type), ]
#
#   ## Seperate genes ##
#   CCM_SPA_Mac_EO_Split.df <- strsplit(as.character(
#                                    CCM_SPA_Mac_EO.df$Marker), ", ")
#   rm( CCM_SPA_Mac.df,CCM_SPA_Mac_EO.df,CCM_SPA_Mac_LO.df, CCM_SPA_Mac_EO_Split.df)

##### Count gene number #####
  source("FUN_CountGeneNum.R")

  ##--------------------- SPA ---------------------##
  CCMTable = function( CCM_SPA_CT.set, CCM_SPA.set
  ){
    CCM_SPA_Split.lt <- strsplit(as.character(CCM_SPA.set), ", ")
    names(CCM_SPA_Split.lt) <- CCM_SPA_CT.set

    ## EO + LO
    CCM_SPA_Count.df <- CountGeneNum(CCM_SPA_Split.lt,mode=2)
    ## EO
    CCM_SPA_Split_EO.lt <- CCM_SPA_Split.lt[grep("EO", names(CCM_SPA_Split.lt))]
    CCM_SPA_Count_EO.df <- CountGeneNum(CCM_SPA_Split_EO.lt,mode=2)
    ## LO
    CCM_SPA_Split_LO.lt <- CCM_SPA_Split.lt[grep("LO", names(CCM_SPA_Split.lt))]
    CCM_SPA_Count_LO.df <- CountGeneNum(CCM_SPA_Split_LO.lt,mode=2)
    CCM_SPA_Count_LO.df$Freq <- -CCM_SPA_Count_LO.df$Freq

    # Combine
    CCM_SPA_Count_Summary.df <- full_join(CCM_SPA_Count.df, CCM_SPA_Count_EO.df, by="Var1")
    CCM_SPA_Count_Summary.df <- full_join(CCM_SPA_Count_Summary.df, CCM_SPA_Count_LO.df, by="Var1")
    colnames(CCM_SPA_Count_Summary.df) <- c("Genes","Total","EO","LO")
    CCM_SPA_Count_Summary.df[is.na(CCM_SPA_Count_Summary.df)] <- 0

    ## For each cell type
    CCM_SPA_Split_EO_CT.df <- CountGeneNum(CCM_SPA_Split_EO.lt,mode=1)
    CCM_SPA_Split_LO_CT.df <- CountGeneNum(CCM_SPA_Split_LO.lt,mode=1)
    CCM_SPA_Split_LO_CT.df[,2:ncol(CCM_SPA_Split_LO_CT.df)] <- -CCM_SPA_Split_LO_CT.df[,2:ncol(CCM_SPA_Split_LO_CT.df)]

    # Combine
    CCM_SPA_SumCT.df <- full_join(CCM_SPA_Count_Summary.df,
                                                  CCM_SPA_Split_EO_CT.df,by="Genes")
    CCM_SPA_SumCT.df <- full_join(CCM_SPA_SumCT.df,
                                                  CCM_SPA_Split_LO_CT.df,by="Genes")
    CCM_SPA_SumCT.df[is.na(CCM_SPA_SumCT.df)] <- 0

    ##
    rm(CCM_SPA_Count.df, CCM_SPA_Split_EO.lt, CCM_SPA_Split_LO.lt,
       CCM_SPA_Count_EO.df,CCM_SPA_Count_LO.df,CCM_SPA_Split_EO_CT.df,CCM_SPA_Split_LO_CT.df)
    return(CCM_SPA_SumCT.df)
    }


  CCM_SPA_CT.set <- CCMList_SPA.df[,1]
  CCM_SPA.set <- CCMList_SPA.df$Marker
  CCM_SPA_SumCT.df <- CCMTable( CCM_SPA_CT.set,CCM_SPA.set)

  ## Modify the colnames

  colnames(CCM_SPA_SumCT.df) <- str_replace(colnames(CCM_SPA_SumCT.df),
                                                      "EO","EO_SPA")
  colnames(CCM_SPA_SumCT.df) <- str_replace(colnames(CCM_SPA_SumCT.df),
                                                      "LO","LO_SPA")
  colnames(CCM_SPA_SumCT.df) <- str_replace(colnames(CCM_SPA_SumCT.df),
                                                       "^EO","All_EO")
  colnames(CCM_SPA_SumCT.df) <- str_replace(colnames(CCM_SPA_SumCT.df),
                                                       "^LO","All_LO")

  ## Check
  # CCM_SPA_SumCT.df2 <- CCMTable( CCM_SPA_CT.set,CCM_SPA.set)
  # sum(CCM_SPA_SumCT.df == CCM_SPA_SumCT.df2)

  ##--------------------- SSA ---------------------##
  ## SSA_Intersect_All
    CCM_SSA_CT.set <- CCMList_SSA.df[,1]
    CCMList_SSA.df[CCMList_SSA.df$Intersect=="",]$Intersect <- "NA"
    CCM_SSA_Inter.set <- CCMList_SSA.df$Intersect
    CCM_SSA_I_SumCT.df <- CCMTable(CCM_SSA_CT.set, CCM_SSA_Inter.set)
    CCM_SSA_I_SumCT.df <- CCM_SSA_I_SumCT.df[!CCM_SSA_I_SumCT.df$Genes=="NA",]

  ## SSA_Female_All
    CCM_SSA_F_CT.set <- CCMList_SSA.df[,1]
    #CCMList_SSA.df[CCMList_SSA.df$Female=="",]$Female <- 0 #"NA"
    CCM_SSA_F.set <- CCMList_SSA.df$Female
    CCM_SSA_F_SumCT.df <- CCMTable(CCM_SSA_F_CT.set, CCM_SSA_F.set)

  ## SSA_Male_All
    CCM_SSA_M_CT.set <- CCMList_SSA.df[,1]
    #CCMList_SSA.df[CCMList_SSA.df$Male=="",]$Male <- 0  #"NA"
    CCM_SSA_M.set <- CCMList_SSA.df$Male
    CCM_SSA_M_SumCT.df <- CCMTable(CCM_SSA_M_CT.set, CCM_SSA_M.set)

    CCM_SSA_SumCT.df <- full_join(CCM_SSA_F_SumCT.df,
                                            CCM_SSA_M_SumCT.df,by=c("Genes"))

    CCM_SSA_SumCT.df <- full_join(CCM_SSA_SumCT.df,
                                            CCM_SSA_I_SumCT.df,by=c("Genes"))

    colnames(CCM_SSA_SumCT.df) <- str_replace(colnames(CCM_SSA_SumCT.df),
                                                        "\\.x","_F")
    colnames(CCM_SSA_SumCT.df) <- str_replace(colnames(CCM_SSA_SumCT.df),
                                                        "\\.y","_M")
    colnames(CCM_SSA_SumCT.df) <- str_replace(colnames(CCM_SSA_SumCT.df),
                                                         "^EO","All_EO")
    colnames(CCM_SSA_SumCT.df) <- str_replace(colnames(CCM_SSA_SumCT.df),
                                                         "^LO","All_LO")

    colnames(CCM_SSA_SumCT.df) <- str_replace(colnames(CCM_SSA_SumCT.df),
                                                   "EO$","EO_I")
    colnames(CCM_SSA_SumCT.df) <- str_replace(colnames(CCM_SSA_SumCT.df),
                                                   "LO$","LO_I")

    # colnames(CCM_SSA_SumCT.df)[56:82] <- paste0(colnames(CCM_SSA_SumCT.df)[56:82],"_I")
    # head(CCM_SSA_SumCT.df)


    CCM_SSA_SumCT.df[is.na(CCM_SSA_SumCT.df)] <- 0
    head(CCM_SSA_SumCT.df)

##### Export the result #####
  # write.csv( CCM_SPA_SumCT.df[,-2] ,
  #            file = paste0(Save.Path,"/KLC_",SampleType, "_CCMCount_SPA_CT_Num.csv"),
  #            #sep = ",",
  #            quote = F,
  #            row.names = F
  # )

  write.csv( CCM_SPA_SumCT.df[,!(colnames(CCM_SPA_SumCT.df) %in%
                                        c("Total"))] ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMCount_SSA_CT_Num.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )


  write.csv( CCM_SSA_SumCT.df[,!(colnames(CCM_SSA_SumCT.df) %in%
                                            c("Total_F","Total_M","Total"))] ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMCount_SSA_CT_Num.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )


##### Combine two dataframe #####
  SPA <- CCM_SPA_SumCT.df[,!(colnames(CCM_SPA_SumCT.df) %in% c("Total"))]
  SSA <- CCM_SSA_SumCT.df[,!(colnames(CCM_SSA_SumCT.df) %in%
                               c("Total_F","Total_M","Total"))]
  CCM_All_SumCT.df <- data.frame(Genes=union(SPA$Genes,SSA$Genes))
  CCM_All_SumCT.df <- CCM_All_SumCT.df %>% left_join(.,SPA,by="Genes") %>%
                      left_join(.,SSA,by="Genes")
  CCM_All_SumCT.df[is.na(CCM_All_SumCT.df)] <- 0

  rm(SPA,SSA)

  # Strsplit # Refer: https://officeguide.cc/r-split-data-frame-string-column-into-multiple-columns-tutorial-examples/
  library(stringr)
  row.names(CCM_All_SumCT.df) <- CCM_All_SumCT.df[,1]
  CCM_All_SumCT.df <- CCM_All_SumCT.df[,-1]

  CCM_All_SumCT_Fin.df <- CCM_All_SumCT.df %>% t() %>% as.data.frame()
  CCM_All_SumCT_Fin.df <- data.frame(Pheno = row.names(CCM_All_SumCT_Fin.df),CCM_All_SumCT_Fin.df)

  All_PhenoSplit.df <- str_split_fixed(CCM_All_SumCT_Fin.df$Pheno, "_", 3)
  CCM_All_SumCT_Fin.df <- data.frame(All_PhenoSplit.df,CCM_All_SumCT_Fin.df)
  colnames(CCM_All_SumCT_Fin.df) <- c("Pheno1","Pheno2","Pheno3",
                              colnames(CCM_All_SumCT_Fin.df)[4:ncol(CCM_All_SumCT_Fin.df)])

  write.csv( CCM_All_SumCT_Fin.df ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMCount_SPA_SSA_CT_Num.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )

  rm(All_PhenoSplit.df, CCM_SPA_CT.set, CCM_SSA_CT.set, CCM_SPA.set, CCM_SSA_Inter.set,
     CCM_SSA_F_CT.set, CCM_SSA_M_CT.set,CCM_SSA_F.set, CCM_SSA_M.set)


#### Clean up the obj #####
  CCM_SumCT.lt <- list(CCM_All_SumCT_Fin.df, CCM_All_SumCT.df,
                       CCM_SPA_SumCT.df, CCM_SSA_SumCT.df,
                       CCM_SSA_I_SumCT.df, CCM_SSA_F_SumCT.df, CCM_SSA_M_SumCT.df)

  names(CCM_SumCT.lt) <- c("CCM_All_SumCT_Fin.df", "CCM_All_SumCT.df",
                           "CCM_SPA_SumCT.df", "CCM_SSA_SumCT.df",
                           "CCM_SSA_I_SumCT.df", "CCM_SSA_F_SumCT.df", "CCM_SSA_M_SumCT.df")

  rm(  CCM_All_SumCT.df,
       CCM_SPA_SumCT.df, CCM_SSA_SumCT.df,
       CCM_SSA_I_SumCT.df, CCM_SSA_F_SumCT.df, CCM_SSA_M_SumCT.df)

