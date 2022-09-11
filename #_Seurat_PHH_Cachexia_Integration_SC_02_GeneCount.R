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
  Version = paste0(Sys.Date(),"_","SC_Main")
  Save.Path = paste0(getwd(),"/",Version)
  dir.create(Save.Path)

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

  ##### Load RData  #####
  load(paste0(Save.Path,"/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData"))

##### SSA CCMarker genes( Exclude genes express opposite in Sex) #####
  ##-------------- Pos --------------##
  CCMarker_SSA_Pos.df = as.data.frame(matrix(nrow=length(names(Venn_CCMarker_Pos)),ncol= 4))
  row.names(CCMarker_SSA_Pos.df) <- names(Venn_CCMarker_Pos)
  colnames(CCMarker_SSA_Pos.df) <- c("Male","Female","Intersect","Union")

  for (i in 1:length(Venn_CCMarker_Pos)) {
    CCMarker_SSA_Pos.df[i,1] <- paste(c(Venn_CCMarker_Pos[[i]][["Summary"]][["Unique_A"]]), collapse = ", ")
    CCMarker_SSA_Pos.df[i,2] <- paste(c(Venn_CCMarker_Pos[[i]][["Summary"]][["Unique_B"]]), collapse = ", ")
    CCMarker_SSA_Pos.df[i,3] <- paste(c(Venn_CCMarker_Pos[[i]][["Summary"]][["Intersect_AB"]]), collapse = ", ")
    CCMarker_SSA_Pos.df[i,4] <- paste(c(Venn_CCMarker_Pos[[i]][["Summary"]][["Union_AB"]]), collapse = ", ")

  }
  rm(i)

  ##-------------- Neg --------------##
  CCMarker_SSA_Neg.df = as.data.frame(matrix(nrow=length(names(Venn_CCMarker_Neg)),ncol= 4))
  row.names(CCMarker_SSA_Neg.df) <- names(Venn_CCMarker_Neg)
  colnames(CCMarker_SSA_Neg.df) <- c("Male","Female","Intersect","Union")

  for (i in 1:length(Venn_CCMarker_Neg)) {
    CCMarker_SSA_Neg.df[i,1] <- paste(c(Venn_CCMarker_Neg[[i]][["Summary"]][["Unique_A"]]), collapse = ", ")
    CCMarker_SSA_Neg.df[i,2] <- paste(c(Venn_CCMarker_Neg[[i]][["Summary"]][["Unique_B"]]), collapse = ", ")
    CCMarker_SSA_Neg.df[i,3] <- paste(c(Venn_CCMarker_Neg[[i]][["Summary"]][["Intersect_AB"]]), collapse = ", ")
    CCMarker_SSA_Neg.df[i,4] <- paste(c(Venn_CCMarker_Neg[[i]][["Summary"]][["Union_AB"]]), collapse = ", ")

  }
  rm(i)

  ##-------------- Combine --------------##
  CCMList_SSA.df  <- rbind(CCMarker_SSA_Pos.df,CCMarker_SSA_Neg.df)
  CCMList_SSA.df  <- data.frame(Type = row.names(CCMList_SSA.df ),CCMList_SSA.df )

  rm(CCMarker_SSA_Pos.df,CCMarker_SSA_Neg.df)

  ## Export tsv
  write.table( CCMList_SSA.df  ,
               file = paste0(Save.Path,"/",SampleType,"_CCMarker_SSA.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )

##### SPA CCMarker genes #####
  ##-------------- Pos --------------##
    CCMarker_SPA_Pos.lt <- list()
    for(i in c(1:length(CellType.list))){
      try({
        CCMarker_SPA_Pos.lt[[i]] <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.S_Pos_List"]]
        names(CCMarker_SPA_Pos.lt)[[i]] <- paste0("CCMarker_",CellType.list[i],"_EO")
      })
    }
    rm(i)

    CCMarker_SPA_Pos.df = as.data.frame(matrix(nrow=length(names(CCMarker_SPA_Pos.lt)),ncol= 1))
    row.names(CCMarker_SPA_Pos.df) <- names(CCMarker_SPA_Pos.lt)
    colnames(CCMarker_SPA_Pos.df) <- c("Marker")

    for (i in 1:length(CCMarker_SPA_Pos.lt)) {
      CCMarker_SPA_Pos.df[i,1] <- paste(c(CCMarker_SPA_Pos.lt[[i]]), collapse = ", ")

    }
    rm(i)


  ##-------------- Neg --------------##
    CCMarker_SPA_Neg.lt <- list()
    for(i in c(1:length(CellType.list))){
      try({
        CCMarker_SPA_Neg.lt[[i]] <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.S_Neg_List"]]
        names(CCMarker_SPA_Neg.lt)[[i]] <- paste0("CCMarker_",CellType.list[i],"_LO")
      })
    }
    rm(i)

    CCMarker_SPA_Neg.df = as.data.frame(matrix(nrow=length(names(CCMarker_SPA_Neg.lt)),ncol= 1))
    row.names(CCMarker_SPA_Neg.df) <- names(CCMarker_SPA_Neg.lt)
    colnames(CCMarker_SPA_Neg.df) <- c("Marker")

    for (i in 1:length(CCMarker_SPA_Neg.lt)) {
      CCMarker_SPA_Neg.df[i,1] <- paste(c(CCMarker_SPA_Neg.lt[[i]]), collapse = ", ")
    }
    rm(i)

    CCMList_SPA.df  <- rbind(CCMarker_SPA_Pos.df,CCMarker_SPA_Neg.df)
    CCMList_SPA.df  <- data.frame(Type = row.names(CCMList_SPA.df ),CCMList_SPA.df )

    rm(CCMarker_SPA_Pos.df,CCMarker_SPA_Neg.df,CCMarker_SPA_Pos.lt,CCMarker_SPA_Neg.lt)

  ## Export tsv
    write.table( CCMList_SPA.df  ,
                 file = paste0(Save.Path,"/",SampleType,"_CCMarker_SPA.tsv"),
                 sep = "\t",
                 quote = F,
                 row.names = F
    )


  # ##### Load data #####
  #   CCMList_SPA.df <- read.delim2(paste0(Save.Path,"/SC_CCMarker_SPA.tsv"))
  #   CCMList_SSA.df <- read.delim2(paste0(Save.Path,"/SC_CCMarker_SSA.tsv"))

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
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMCount_SPA_CT_Num.csv"),
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
  colnames(CCM_All_SumCT_Fin.df) <- c("Cell_Type","State","Strategy",
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

  rm( CCM_SPA_SumCT.df, CCM_SSA_SumCT.df,
       CCM_SSA_I_SumCT.df, CCM_SSA_F_SumCT.df, CCM_SSA_M_SumCT.df)

##########################################################################


##### Visualization #####
  library(ggplot2)
  source("FUN_Beautify_ggplot.R")

  plot(x= seq(1:ncol(CCM_All_SumCT.df)), y= as.integer(CCM_All_SumCT.df[1,1:ncol(CCM_All_SumCT.df)]))

  ## EO
  # https://stackoverflow.com/questions/17423609/multiple-y-over-x-in-ggplot2
    #CCM_Count_Sum_EO.df <- CCM_All_SumCT.df[,c(1,3,6)] %>% arrange(desc(All_EO_SPA),desc(All_EO_SSA_I))
    # CCM_Count_Sum_EO.df <- CCM_All_SumCT.df[grep("All_EO", names(CCM_All_SumCT.df))] %>%
    #                      arrange(desc(All_EO_SPA),desc(All_EO_I))

    CCM_Count_Sum_EO.df <- CCM_All_SumCT.df[names(CCM_All_SumCT.df) %in% c("All_EO_SPA","All_EO_I")] %>%
                           arrange(desc(All_EO_SPA),desc(All_EO_I))
    CCM_Count_Sum_EO.df <- data.frame(Genes=row.names(CCM_Count_Sum_EO.df),CCM_Count_Sum_EO.df)
    colnames(CCM_Count_Sum_EO.df)[3] <-"All_EO_SSA_I"

    library(reshape2)
    df1 <- melt(CCM_Count_Sum_EO.df[1:500,], id = "Genes")
    df1$Genes <- factor(df1$Genes,levels = df1$Genes[1:500])

    source("FUN_Beautify_LinePlot_Ori.R")
    EO1.P <- ggplot(df1, aes(x = Genes, y = value, color = variable,
                             group = variable ,linetype = variable)) + theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

    EO1.P <- EO1.P %>% LinePlot_Ori(.,AspRat=1,LegPos = c(0.85, 0.9),AxisTitleSize=2,
                                    XtextSize=0,  YtextSize=25) +
      scale_color_manual(values = c('#f7598b','#9755ab'))
    EO1.P
    ##
    CCMarker_EO1 <- CCM_Count_Sum_EO.df[CCM_Count_Sum_EO.df$All_EO_SPA==1 &CCM_Count_Sum_EO.df$All_EO_SSA_I ==1,]
    CCMarker_EO2 <- CCM_Count_Sum_EO.df[CCM_Count_Sum_EO.df$All_EO_SPA==2 &CCM_Count_Sum_EO.df$All_EO_SSA_I ==2,]
    CCMarker_EO3 <- CCM_Count_Sum_EO.df[CCM_Count_Sum_EO.df$All_EO_SPA==3 &CCM_Count_Sum_EO.df$All_EO_SSA_I ==3,]

  ## LO
  # https://stackoverflow.com/questions/17423609/multiple-y-over-x-in-ggplot2
    #CCM_Count_Sum_LO.df <- CCM_All_SumCT.df[,c(1,4,7)] %>% arrange(desc(LO.SPA),desc(LO.SSA.I))
    CCM_Count_Sum_LO.df <- CCM_All_SumCT.df[names(CCM_All_SumCT.df) %in% c("All_LO_SPA","All_LO_I")] %>%
                           abs() %>%
                           arrange(desc(All_LO_SPA),desc(All_LO_I))
    CCM_Count_Sum_LO.df <- data.frame(Genes=row.names(CCM_Count_Sum_LO.df),CCM_Count_Sum_LO.df)
    colnames(CCM_Count_Sum_LO.df)[3] <-"All_LO_SSA_I"
    library(reshape2)
    df2 <- melt(CCM_Count_Sum_LO.df[1:500,], id = "Genes")
    df2$Genes <- factor(df2$Genes,levels = df2$Genes[1:500])
    LO1.P <- ggplot(df2, aes(x = Genes, y = value, color = variable,
                             group = variable ,linetype = variable))

    LO1.P <-LO1.P %>% LinePlot_Ori(.,AspRat=1,LegPos = c(0.85, 0.9),AxisTitleSize=2,
                                   XtextSize=0,  YtextSize=25) +
      scale_color_manual(values = c('#20ba87','#205eba'))
    LO1.P
    ##
    CCMarker_LO1 <- CCM_Count_Sum_LO.df[CCM_Count_Sum_LO.df$LO.SPA==1 &CCM_Count_Sum_LO.df$LO.SSA.I ==1,]
    CCMarker_LO2 <- CCM_Count_Sum_LO.df[CCM_Count_Sum_LO.df$LO.SPA==2 &CCM_Count_Sum_LO.df$LO.SSA.I ==2,]
    CCMarker_LO3 <- CCM_Count_Sum_LO.df[CCM_Count_Sum_LO.df$LO.SPA==3 &CCM_Count_Sum_LO.df$LO.SSA.I ==3,]


  #### Sex specific
  ## EO
  # CCM_Count_Sum_EO_SSA.df <- CCM_All_SumCT.df[,c(1,6,9,12)] %>% arrange(desc(All_EO_SSA_I),desc(EO.SSA.F),desc(EO.SSA.M))
    CCM_Count_Sum_EO_SSA.df <- CCM_All_SumCT.df[names(CCM_All_SumCT.df) %in% c("All_EO_I","All_EO_F","All_EO_M")] %>%
                            arrange(desc(All_EO_I),desc(All_EO_F),desc(All_EO_M))
    CCM_Count_Sum_EO_SSA.df <- data.frame(Genes=row.names(CCM_Count_Sum_EO_SSA.df),CCM_Count_Sum_EO_SSA.df)
    # colnames(CCM_Count_Sum_EO_SSA.df)[3] <-"All_EO_SSA_I"

    library(reshape2)
    df3 <- melt(CCM_Count_Sum_EO_SSA.df[1:500,], id = "Genes")
    df3$Genes <- factor(df3$Genes,levels = df3$Genes[1:500])
    EOSS.P <- ggplot(df3, aes(x = Genes, y = value, color = variable,
                              group = variable ,linetype = variable))

    EOSS.P <- EOSS.P %>% LinePlot_Ori(.,AspRat=1,LegPos = c(0.85, 0.89),AxisTitleSize=2,
                                      XtextSize=0,  YtextSize=25) +
      geom_point(shape = 12, size = 4, fill = "white") +
      scale_linetype_manual(#name="Pheno_Type",
        values=c("solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
        ) +  # labels=c("All_EO_I","All_EO_F","All_EO_M")
      scale_color_manual(values = c('#9755ab','#ff75da','#6db5fc'))

    CCM_Count_Sum_EO_SSA.df %>% arrange(desc(All_EO_F))
    EOSS.P


  ## LO
    #CCM_Count_Sum_LO_SSA.df <- CCM_All_SumCT.df[,c(1,7,10,13)] %>% arrange(desc(LO.SSA.I),desc(LO.SSA.F),desc(LO.SSA.M))
    CCM_Count_Sum_LO_SSA.df <- CCM_All_SumCT.df[names(CCM_All_SumCT.df) %in% c("All_LO_I","All_LO_F","All_LO_M")] %>% abs()%>%
      arrange(desc(All_LO_I),desc(All_LO_F),desc(All_LO_M))
    CCM_Count_Sum_LO_SSA.df <- data.frame(Genes=row.names(CCM_Count_Sum_LO_SSA.df),CCM_Count_Sum_LO_SSA.df)

    library(reshape2)
    df4 <- melt(CCM_Count_Sum_LO_SSA.df[1:500,], id = "Genes")
    df4$Genes <- factor(df4$Genes,levels = df4$Genes[1:500])
    LOSS.P <- ggplot(df4, aes(x = Genes, y = value, color = variable,
                              group = variable ,linetype = variable))

    LOSS.P <- LOSS.P %>% LinePlot_Ori(.,AspRat=1,LegPos = c(0.85, 0.85),AxisTitleSize=2,
                                      XtextSize=0,  YtextSize=25)+
      geom_point(shape = 12, size = 4, fill = "white") +
      scale_linetype_manual(#name="Pheno_Type",
        values=c("solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
        ) + # labels=c("All_LO_I","All_LO_F","All_LO_M")
      scale_color_manual(values = c('#9755ab','#ff75da','#6db5fc'))

    CCM_Count_Sum_LO_SSA.df %>% arrange(desc(All_LO_F))
    LOSS.P

##### Gene count vs Cluster #####
  # GF: Gene Frequency
    CCM_All_SumCT.df2 <- CCM_All_SumCT.df %>% abs()
    Base_GF <- data.frame(Var1=seq(from = 0, to = 9))
    Base_GF$Var1 <- as.factor(Base_GF$Var1)
    EO_SPA_GF <- data.frame(table(CCM_All_SumCT.df2$All_EO_SPA))
    colnames(EO_SPA_GF)[2] <-  'EO_SPA'
    EO_SPA_GF <- full_join(Base_GF, EO_SPA_GF)

    LO_SPA_GF <- data.frame(table(CCM_All_SumCT.df2$All_LO_SPA))
    colnames(LO_SPA_GF)[2] <-  'LO.SPA'

    Sum_GF <- full_join(EO_SPA_GF, LO_SPA_GF)

    EO_SSAI_GF <- data.frame(table(CCM_All_SumCT.df2$All_EO_I))
    colnames(EO_SSAI_GF)[2] <-  'EO.SSAI'

    Sum_GF <- full_join(Sum_GF, EO_SSAI_GF)

    LO_SSAI_GF <- data.frame(table(CCM_All_SumCT.df2$All_LO_I))
    colnames(LO_SSAI_GF)[2] <-  'LO.SSAI'

    Sum_GF <- full_join(Sum_GF, LO_SSAI_GF)

    EO_SSAF_GF <- data.frame(table(CCM_All_SumCT.df2$All_EO_F))
    colnames(EO_SSAF_GF)[2] <-  'EO.SSAF'

    Sum_GF <- full_join(Sum_GF, EO_SSAF_GF)

    LO_SSAF_GF <- data.frame(table(CCM_All_SumCT.df2$All_LO_F))
    colnames(LO_SSAF_GF)[2] <-  'LO.SSAF'

    Sum_GF <- full_join(Sum_GF, LO_SSAF_GF)

    EO_SSAM_GF <- data.frame(table(CCM_All_SumCT.df2$All_EO_M))
    colnames(EO_SSAM_GF)[2] <-  'EO.SSAM'

    Sum_GF <- full_join(Sum_GF, EO_SSAM_GF)

    LO_SSAM_GF <- data.frame(table(CCM_All_SumCT.df2$All_LO_M))
    colnames(LO_SSAM_GF)[2] <-  'LO.SSAM'

    Sum_GF <- full_join(Sum_GF, LO_SSAM_GF)
    colnames(Sum_GF)[1] <- "GeneCount"
    Sum_GF[is.na(Sum_GF)] <- 0
    rm(Base_GF,EO_SPA_GF, LO_SPA_GF, EO_SSAI_GF, LO_SSAI_GF, EO_SSAF_GF,
       LO_SSAF_GF, EO_SSAM_GF, LO_SSAM_GF)

  ## Plot
    library(ggplot2)
    library(reshape2)
    Sum_GF_Melt <- melt(Sum_GF[-1,], id = "GeneCount")
    #Sum_GF_Melt$GeneCount <- factor(Sum_GF_Melt$GeneCount,levels = Sum_GF_Melt$GeneCount)
    ##
    Sum_GF_Melt$GeneCount <- factor(Sum_GF_Melt$GeneCount)
    colnames(Sum_GF_Melt) <- c("Cluster_Count","Type","Gene_Count")

    Sum_GF_Melt$Cluster_Count <- as.numeric(Sum_GF_Melt$Cluster_Count)
    source("FUN_Beautify_LinePlot.R")
    Sum_GF.P <- ggplot(Sum_GF_Melt, aes(x = Cluster_Count, y = Gene_Count, color = Type,
                                        group = Type ,linetype = Type))

    Sum_GF.P <- Sum_GF.P %>% LinePlot(.,AspRat=1,LegPos = c(0.85, 0.8),AxisTitleSize=2,
                                      XtextSize=15,  YtextSize=15,
                                      scale_y_c =  seq(0, max(Sum_GF_Melt$Gene_Count)+30, 50))+
      scale_linetype_manual(#name="Pheno_Type",
        values=c("solid", "solid", "solid", "solid",
                 "dotted", "dotdash","dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
        labels=c("EO_SPA","LO.SPA","EO.SSAI","LO.SSAI",
                 "EO.SSAF","LO.SSAF","EO.SSAM","LO.SSAM")) +
      scale_color_manual(values = c('#ef476f','#0077b6','#e79bfa','#bc87c9',
                                    '#ff8fea','#c995c0','#3d84ff','#7099e0'))
    Sum_GF.P

  ## Log
  # Line Color Setting
    LTypeSet.lt <- list(
      values = c("solid", "solid", "solid", "solid", "dotted", "dotdash","dotted", "dotdash"), # values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
      labels = c("EO_SPA","LO.SPA","EO.SSAI","LO.SSAI","EO.SSAF","LO.SSAF","EO.SSAM","LO.SSAM")
    )
    LColorSet.lt <- list(
      values = c('#ef476f','#0077b6','#b372c2','#874299','#ff8fea','#bd559c','#7099e0','#151e75')
    )

    Sum_GF_Melt$`Log2(Gene_Count+1)` <- log2(Sum_GF_Melt$Gene_Count+1)
    Sum_GF_Melt$Cluster_Count <- as.numeric(Sum_GF_Melt$Cluster_Count)

    Sum_GF_log.P <- ggplot(Sum_GF_Melt, aes(x = Cluster_Count, y = `Log2(Gene_Count+1)`, color = Type,
                                            group = Type ,linetype = Type))

    Sum_GF_log.P <- Sum_GF_log.P %>% LinePlot(., AspRat=0.5, LegPos = c(0.85, 0.75), AxisTitleSize=1.5,
                                              XtextSize=20,  YtextSize=20,
                                              scale_x_c =  seq(0, max(Sum_GF_Melt$Cluster_Count)),
                                              scale_y_c =  seq(0, max(Sum_GF_Melt$Gene_Count)))+
      scale_linetype_manual(#name="Pheno_Type",
        values= LTypeSet.lt$values,
        labels= LTypeSet.lt$labels) +
      scale_color_manual(values= LColorSet.lt$values) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) # Remove grid

    Sum_GF_log.P


    Sum_GF_log.P2 <- Sum_GF_log.P %>% LinePlot(., AspRat=1, LegPos = c(0.85, 0.75), AxisTitleSize=1.5,
                                              XtextSize=20,  YtextSize=20,
                                              scale_x_c =  seq(0, max(Sum_GF_Melt$Cluster_Count)),
                                              scale_y_c =  seq(0, max(Sum_GF_Melt$Gene_Count)))+
      scale_linetype_manual(#name="Pheno_Type",
        values= LTypeSet.lt$values,
        labels= LTypeSet.lt$labels) +
      scale_color_manual(values= LColorSet.lt$values) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) # Remove grid

    Sum_GF_log.P2

  ##### Export PDF #####

  pdf(file = paste0(Save.Path,"/",SampleType,"_GeneCount.pdf"),width = 7, height = 7 )
    Sum_GF_log.P2
    Sum_GF.P
    EO1.P
    LO1.P
    EOSS.P
    LOSS.P
  dev.off()

  pdf(file = paste0(Save.Path,"/",SampleType,"_GeneCount_Main.pdf"),width = 14, height = 7 )
   Sum_GF_log.P
  dev.off()


##### Heatmap #####
  library(ComplexHeatmap)
  # H_CCM_SPA_CT.df <- CCM_SPA_Count_Sum_CT.df
  # H_CCM_SPA_CT.df <- left_join(H_CCM_SPA_CT.df,CCM_All_SumCT.df)
  #
  # row.names(H_CCM_SPA_CT.df) <- H_CCM_SPA_CT.df[,1]
  # H_CCM_SPA_CT.df <- H_CCM_SPA_CT.df[,-1]

  H_CCM_SPA_CT.df <- CCM_All_SumCT.df
  H_CCM_SPA_CT.df <- as.data.frame(t(H_CCM_SPA_CT.df))
  Heatmap(H_CCM_SPA_CT.df)

  CCM_TOP_EO.set <- c(as.character(CCMarker_EO1$Genes),as.character(CCMarker_EO2$Genes),
                      as.character(CCMarker_EO3$Genes),"Chil3","S100a8")
  CCM_TOP_LO.set <- c(as.character(CCMarker_LO1$Genes),as.character(CCMarker_LO2$Genes),
                      as.character(CCMarker_LO3$Genes))

  CCM_TOP.set <- c(CCM_TOP_EO.set,CCM_TOP_LO.set)

  CCM_Count_Sum_T_S.df <- H_CCM_SPA_CT.df[,colnames(H_CCM_SPA_CT.df) %in% CCM_TOP.set]

  library(circlize)
  col_fun = colorRamp2(c(0, 6, 12), c("white", "#B15BFF", "#5B00AE"))
  col_fun(seq(-3, 3))

  Heatmap(CCM_Count_Sum_T_S.df, name = "Num", col = col_fun)

##### Clean up the obj #####
  rm(CCMarker_EO1, CCMarker_EO2, CCMarker_EO3,
     CCMarker_LO1, CCMarker_LO2, CCMarker_LO3,
     df1, df2, df3, df4, EO1.P, EOSS.P,
     LO1.P, LOSS.P, CCM_All_SumCT.df2)

# ##### Export R.Data ######
#   save.image(paste0(Save.Path,"/S01_GeneCount.RData"))


##### SPA_SSA_Cell_Type_Match #####


  ##### Load data #####
  ## Way 1
  CCM_SPA.df <- read.delim2(paste0(Save.Path,"/KLC_",SampleType, "_CCMCount_SPA_CT_Num.csv"),
                            sep = ",")

  CCM_SSA.df <- read.delim2(paste0(Save.Path,"/KLC_",SampleType, "_CCMCount_SSA_CT_Num.csv"),
                            sep = ",")

  ## Way 2
  CCM_SPA_SumCT.df <- CCM_SumCT.lt[["CCM_SPA_SumCT.df"]]
  CCM_SPA.df <- CCM_SPA_SumCT.df[,!(colnames(CCM_SPA_SumCT.df) %in% c("Total"))]

  CCM_SSA_SumCT.df <- CCM_SumCT.lt[["CCM_SSA_SumCT.df"]]
  CCM_SSA.df <-CCM_SSA_SumCT.df[,!(colnames(CCM_SSA_SumCT.df) %in%
                                c("Total_F","Total_M","Total"))]
  ##### Prepossess #####
  row.names(CCM_SPA.df) <- CCM_SPA.df[,1]
  row.names(CCM_SSA.df) <- CCM_SSA.df[,1]

  CCM_All.df <- full_join(CCM_SPA.df,CCM_SSA.df,by="Genes")

  CCM_All.df[is.na(CCM_All.df)] <- 0

  ##### Plot Heatmap #####
  library(circlize)
  col_fun = colorRamp2(c(-12, -6, 0, 6, 12),
                       c("#003bbd","#5183f0","white", "#f051c2", "#bd0087"))
  #col_fun(seq(-3, 3))

  Heatmap(CCM_SPA.df[,-1], name = "Num", col = col_fun)

  ##### Count gene that Match in Cell Type by different Analysis #####

  ##
  #PhenoType_SPA.Set <- colnames(CCM_All.df)[4:27]
  PhenoType_SPA.Set <- grep("SPA",colnames(CCM_All.df),value=T)
  TempD <- grep("All",PhenoType_SPA.Set,value=T)
  PhenoType_SPA.Set <- PhenoType_SPA.Set[!PhenoType_SPA.Set %in% TempD]
  PhenoType_SPA.Set <- gsub("_SPA","",PhenoType_SPA.Set)
  ncol_Ori <- ncol(CCM_All.df)
  CCM_All_M.df <- CCM_All.df

  for (i in 1:length(PhenoType_SPA.Set)) {
    CCM_All_M.df <- CCM_All_M.df %>% mutate(Temp =
                                              abs(CCM_All_M.df[paste0(PhenoType_SPA.Set[i],"_SPA")]) == 1 &
                                              abs(CCM_All_M.df[paste0(PhenoType_SPA.Set[i],"_I")]) == 1)
    CCM_All_M.df$Temp <- as.numeric(CCM_All_M.df$Temp)
    colnames(CCM_All_M.df)[ncol_Ori+i] <- paste0(PhenoType_SPA.Set[i],"_Ma")
  }
  rm(i)

  # Check
  # t(CCM_All_M.df) %>% View()

  ## Summary
  CT_Ma.Set <- grep("_Ma",colnames(CCM_All_M.df),value=T)
  CCM_All_M.df$Match_Sum <- rowSums(CCM_All_M.df[,CT_Ma.Set])

  CT_Ma_EO.Set <- grep("EO",CT_Ma.Set,value=T)
  CCM_All_M.df$Match_EO_Sum <- rowSums(CCM_All_M.df[,CT_Ma_EO.Set])

  CT_Ma_LO.Set <- grep("LO",CT_Ma.Set,value=T)
  CCM_All_M.df$Match_LO_Sum <- rowSums(CCM_All_M.df[,CT_Ma_LO.Set])

  CCM_All_M_S.df <- CCM_All_M.df[,c("Genes",CT_Ma.Set,
                   "Match_Sum","Match_EO_Sum", "Match_LO_Sum")] %>%
                    arrange(desc(Match_Sum))
  CCM_All_M_S.df <- CCM_All_M_S.df[!CCM_All_M_S.df$Match_Sum == 0,]

  # CCM_All_M.df$Match_Sum <- rowSums(CCM_All_M.df[,107:130])
  # CCM_All_M.df$Match_EO_Sum <- rowSums(CCM_All_M.df[,107:118])
  # CCM_All_M.df$Match_LO_Sum <- rowSums(CCM_All_M.df[,119:130])
  # CCM_All_M_S.df <- CCM_All_M.df[,c(1,107:133)] %>% arrange(desc(Match_Sum))
  # CCM_All_M_S.df <- CCM_All_M_S.df[!CCM_All_M_S.df$Match_Sum == 0,]

  CCM_All_M_S_LO.df <- CCM_All_M_S.df[,grep("_LO",colnames(CCM_All_M_S.df))]
  CCM_All_M_S_LO.df <- -CCM_All_M_S_LO.df
  # CCM_All_M_S.df[,c(14:25,28)] <- CCM_All_M_S_LO.df
  # CCM_All_M_S.df[,grep("_LO",colnames(CCM_All_M_S.df))] <- CCM_All_M_S_LO.df
  rm(CCM_All_M_S_LO.df)

  CCM_All_M_S.df[str_detect(colnames(CCM_All_M_S.df),"LO")]  <- -CCM_All_M_S.df[str_detect(colnames(CCM_All_M_S.df),"LO")]

  ## Plot Heatmap
  minC <- min(CCM_All_M_S.df$Match_LO_Sum)
  maxC <- max(CCM_All_M_S.df$Match_EO_Sum)
  col_fun_Match = colorRamp2(c(minC,minC/2, 0,maxC/2,maxC),
                             c("#003bbd","#5183f0","white", "#f051c2", "#bd0087"))
  Heatmap(CCM_All_M_S.df[,-1], name = "Num", col = col_fun_Match)

  ##### Grab candidate gene by Level #####
  CCM_EO_1.set <- CCM_All_M_S.df[CCM_All_M_S.df$Match_EO_Sum > 1, 1]
  CCM_EO_2.set <- CCM_All_M_S.df[CCM_All_M_S.df$Match_EO_Sum <= 1 &
                                   CCM_All_M_S.df$Match_EO_Sum > 0  ,1]
  CCM_LO_1.set <- CCM_All_M_S.df[abs(CCM_All_M_S.df$Match_LO_Sum) > 1, 1]
  CCM_LO_2.set <- CCM_All_M_S.df[abs(CCM_All_M_S.df$Match_LO_Sum) <= 1 &
                                   abs(CCM_All_M_S.df$Match_LO_Sum) > 0,1]

  CCM_EO_1.str <- str_c(CCM_EO_1.set, collapse = ",")
  CCM_EO_2.str <- str_c(CCM_EO_2.set, collapse = ",")
  CCM_LO_1.str <- str_c(CCM_LO_1.set, collapse = ",")
  CCM_LO_2.str <- str_c(CCM_LO_2.set, collapse = ",")

  CCM_Level.df <- data.frame(EO1=CCM_EO_1.str,EO2=CCM_EO_2.str,
                             LO1=CCM_LO_1.str,LO2=CCM_LO_2.str) %>% t()
  colnames(CCM_Level.df) <- "CCMarker"

  CCM_Level_Num.df <- data.frame(EO1=length(CCM_EO_1.set),EO2=length(CCM_EO_2.set),
                                 LO1=length(CCM_LO_1.set),LO2=length(CCM_LO_2.set)) %>% t() %>% as.data.frame()
  colnames(CCM_Level_Num.df) <- "GeneNum"
  CCM_Level_Num.df$Type <- row.names(CCM_Level_Num.df)

  ## Bar plot
  library(ggplot2)
  CCM_Level_Num.p <- ggplot(data=CCM_Level_Num.df, aes(x=Type, y=GeneNum)) +
                            geom_bar(stat="identity")
  CCM_Level_Num.p

  rm(CCM_EO_1.set,CCM_EO_2.set,CCM_LO_1.set,CCM_LO_2.set,
     CCM_EO_1.str,CCM_EO_2.str,CCM_LO_1.str,CCM_LO_2.str)

  ##### Grab candidate gene by PhenoType #####
  ##
  PhenoType_Ma.Set <- colnames(CCM_All_M_S.df)
  PhenoType_Ma.Set <- PhenoType_Ma.Set[-c(1,(length(CCM_All_M_S.df)-2):length(CCM_All_M_S.df))]

  for (i in 1:length(PhenoType_Ma.Set)) {

    # CCM_MacM0_EO.df <- data.frame(CCM_All_M_S.df[CCM_All_M_S.df[PhenoType_Ma.Set[i]] == 1, 1])
    # colnames(CCM_MacM0_EO.df) <- PhenoType_Ma.Set[i]
    CCM_MacM0_EO.set <- CCM_All_M_S.df[abs(CCM_All_M_S.df[PhenoType_Ma.Set[i]]) == 1, 1]
    CCM_MacM0_EO.str <- str_c(CCM_MacM0_EO.set, collapse = ",")

    ## CCMarker set
    if(i==1){
      CCM_CT.df <- data.frame(Temp = CCM_MacM0_EO.str) %>% t()
      row.names(CCM_CT.df)[1] <- PhenoType_Ma.Set[i]
    }else{
      CCM_CT.df <- data.frame(Temp = CCM_MacM0_EO.str,CCM_CT.df%>% t()) %>% t()
      row.names(CCM_CT.df)[1] <- PhenoType_Ma.Set[i]
    }

    # CCM_CT_Num.df
    if(i==1){
      CCM_CT_Num.df <- data.frame(Temp=length(CCM_MacM0_EO.set)) %>%
        t() %>% as.data.frame()
      colnames(CCM_CT_Num.df)[1] <- "GeneNum"
      CCM_CT_Num.df$Type <- PhenoType_Ma.Set[i]
      row.names(CCM_CT_Num.df)[1] <-PhenoType_Ma.Set[i]
    }else{
      CCM_CT_Num_Temp.df <- data.frame(Temp=length(CCM_MacM0_EO.set)) %>%
        t() %>% as.data.frame()
      colnames(CCM_CT_Num_Temp.df)[1] <- "GeneNum"
      CCM_CT_Num_Temp.df$Type <- PhenoType_Ma.Set[i]
      row.names(CCM_CT_Num_Temp.df)[1] <-PhenoType_Ma.Set[i]

      CCM_CT_Num.df <- rbind(CCM_CT_Num_Temp.df, CCM_CT_Num.df)
    }
  }
  rm(i)

  CCM_CT_Num.df


  CCM_CT_Num_PSplit.df  <- str_split_fixed(CCM_CT_Num.df$Type, "_", 3)
  CCM_CT_Num.df <- data.frame(CCM_CT_Num_PSplit.df ,CCM_CT_Num.df)
  colnames(CCM_CT_Num.df) <- c("Cell_Type","State","Strategy",
                               colnames(CCM_CT_Num.df)[4:ncol(CCM_CT_Num.df)])
  CCM_CT_Num.df$Type2 <- CCM_CT_Num.df$Type
  CCM_CT_Num.df$Type2 <- str_replace(CCM_CT_Num.df$Type2,"_Ma","")



  ## Bar plot
  # Ref: https://officeguide.cc/r-ggplot2-bar-plot-tutorial-examples/
  require(ggplot2)

  # Plot1
    CCM_CT_Num.p <- ggplot(data=CCM_CT_Num.df, aes(x=Type2, y=GeneNum, fill=State)) +
      geom_bar(stat="identity")+
      #theme(axis.title.x = element_text(angle = 90,face="italic",colour = "darkred",size=14)) +
      theme(axis.text.x = element_text(face="bold",  size = 12,angle = 90,
                                       hjust = 1, vjust = .5)) + # Change the size along the x axis
      scale_fill_manual(values = c("#ef476f", "#0077b6"))

    CCM_CT_Num.p <- CCM_CT_Num.p %>% BeautifyggPlot(., AspRat=1, LegPos = c(0.9, 0.9), AxisTitleSize=1.5,
                              XtextSize= 15,  YtextSize= 20, xangle = 90,
                              LegTitleSize= 20 ,LegTextSize = 17) + labs(x="") +
                          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) # Remove grid

    CCM_CT_Num.p

    rm(CCM_CT_Num_Temp.df)

  # Plot2
    CCM_CT_Num.p2 <- ggplot(data=CCM_CT_Num.df, aes(x=Cell_Type, y=GeneNum, fill=State)) +
      geom_bar(stat="identity")+
      #theme(axis.title.x = element_text(angle = 90,face="italic",colour = "darkred",size=14)) +
      theme(axis.text.x = element_text(face="bold",  size = 12,angle = 90,
                                       hjust = 1, vjust = .5)) + # Change the size along the x axis
      scale_fill_manual(values = c("#ef476f", "#0077b6"))

    CCM_CT_Num.p2 <- CCM_CT_Num.p2 %>% BeautifyggPlot(., AspRat=1, LegPos = c(0.9, 0.9), AxisTitleSize=1.5,
                                                    XtextSize= 15,  YtextSize= 20, xangle = 90,
                                                    LegTitleSize= 20 ,LegTextSize = 17)+ labs(x="") +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) # Remove grid

    CCM_CT_Num.p2


  CCM.df <- CCM_CT.df
  CCM.df <- data.frame(Type = row.names(CCM.df), CCMarker = CCM.df)
  CCM.df$Type <- str_replace(CCM.df$Type,".T_","+T_")
  CCM.df <- full_join( CCM_CT_Num.df,CCM.df, by= "Type")
  # CCM_Level_Num.df <- data.frame(Type=row.names(CCM_Level_Num.df),CCM_Level_Num.df)

  CCM.df2 <- CCM_Level.df
  CCM.df2 <- data.frame(Type = row.names(CCM.df2), CCMarker = CCM.df2)
  CCM.df2 <- full_join( CCM_Level_Num.df,CCM.df2, by= "Type")
  CCM.df2 <- data.frame(Cell_Type="", State="", Strategy="",
                        GeneNum=CCM.df2[,1],Type=CCM.df2[,2],Type2="", CCMarker=CCM.df2[,3] )
  CCM.df <- rbind(CCM.df2,CCM.df)
  # ## Try
  # CCM_MacM0_EO.df <- data.frame(CCM_All_M_S.df[CCM_All_M_S.df$MacM0_EO_Ma == 1, 1])
  # colnames(CCM_MacM0_EO.df) <- "MacM0_EO_Ma"
  # CCM_MacM0_EO.set <- CCM_All_M_S.df[CCM_All_M_S.df$MacM0_EO_Ma == 1, 1]
  # CCM_MacM0_EO.str <- str_c(CCM_MacM0_EO.set, collapse = ",")
  #
  # CCM_CT_Num.df <- data.frame(EO1=CCM_MacM0_EO.str,CCM_CT_Num.df) %>% t()
  # row.names(CCM_Level.df)[1] <- "MacM0_EO_Ma"




  # ## Try
  # sum(CCM_All_M.df$MacM1_EO_Ma)
  # Genelist <- CCM_All_M.df[CCM_All_M.df$MacM1_EO_Ma == 1,1]
  # library(stringr)
  # Genelist2 <- str_c(Genelist, collapse = ", ")

  ## Total Table
  CCM_All_M2.df <- full_join(CCM_All.df, CCM_All_M_S.df, by="Genes")
  CCM_All_M2.df[is.na(CCM_All_M2.df)] <- 0

  ##### Export the result #####

  ##### Export PDF #####
    pdf(file = paste0(Save.Path,"/",SampleType,"_CCMark_CTCount.pdf"),width = 7, height = 7 )
      CCM_CT_Num.p
      CCM_CT_Num.p2
    dev.off()

  ##### Table #####
  ## Table For Heatmap
  write.csv( CCM_All_M_S.df ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMHeatmap_SPA_SSA_Match_S.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )

  write.csv( CCM_All_M2.df , # CCM_All_M2.df
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMHeatmap_SPA_SSA_Match_All.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )

  ## Table For Marker

  write.csv( CCM.df ,
             file = paste0(Save.Path,"/#_",SampleType, "_CCM_lt_SUM.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )

  write.table( CCM.df ,
               file = paste0(Save.Path,"/#_",SampleType, "_CCM_lt_SUM.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )



##### CCMarker bar plot #####
  Idents(SC.combined) <- SC.combined$celltype
  TTT <- CCM_All_M2.df %>% arrange(desc(Match_Sum))
  TTT <- TTT[TTT$Match_Sum>=2,]

  TTT$CCState <- ""
  TTT[TTT$Match_EO_Sum >0,]$CCState <- "EO"
  TTT[TTT$Match_EO_Sum <=0,]$CCState <- "LO"

  TTT$Genes <- factor(TTT$Genes,levels = TTT$Genes)
  TTT.p <- ggplot(data=TTT, aes(x=Genes, y=Match_Sum ,fill=CCState))  +
    geom_bar(stat="identity")+
    #theme(axis.title.x = element_text(angle = 90,face="italic",colour = "darkred",size=14)) +
    theme(axis.text.x = element_text(face="bold",  size = 12,angle = 90,
                                     hjust = 1, vjust = .5)) + # Change the size along the x axis
    scale_fill_manual(values = c("#ef476f", "#0077b6"))
  TTT.p
  TTT.p <- TTT.p %>% BeautifyggPlot(., AspRat=1, LegPos = c(0.9, 0.9), AxisTitleSize=1.5,
                                    XtextSize= 15,  YtextSize= 20, xangle = 90,
                                    LegTitleSize= 20 ,LegTextSize = 17) + labs(x="") +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) # Remove grid
  TTT.p

  TTT.DotPlot <- DotPlot(SC.combined, features = TTT$Genes, cols = c("#ef476f", "#0077b6"), dot.scale = 8, split.by = "Cachexia") + RotatedAxis()
  TTT.DotPlot <- TTT.DotPlot%>%
    BeautifyggPlot(.,LegPos = "bottom",AxisTitleSize=1, TitleSize = 20, xangle =90,
                   LegDir = "horizontal",SubTitSize = 17 , LegTextSize = 14, XaThick=1, YaThick=1,XtextSize=12,  YtextSize=12)
  TTT.DotPlot

  pdf(file = paste0(Save.Path,"/SC_CCM_Dot.pdf"),
      width = 10, height = 10 )
  TTT.p
  TTT.DotPlot
  dev.off() # graphics.off()

  # ##### Export R.Data ######
  #   save.image(paste0(Save.Path,"/S01_GeneCount.RData"))
