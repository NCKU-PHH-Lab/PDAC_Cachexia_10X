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

##### Path and folder setting #####
  Import.Path = paste0(getwd(),"/2022-09-07_PBMC_Main")
  Save.Path = paste0(getwd(),"/2022-09-07_PBMC_Main")
  dir.create(Save.Path)
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
  source("FUN_CountGeneNum.R")
  source("FUN_Beautify_LinePlot_Ori.R")

##### Abbreviated Notes #####
  # CCM: Cancer Cachexia Marker
  # SPA: Sex Pooled Analysis
  # SSA: Sex Separated Analysis
  # EO: Early Onset
  # LO: Late Onset
  # CT: Cell Type

  # list:lt
  # Sum:Sum

##### Load data #####
  CCM_SPA.df <- read.delim2(paste0(Import.Path,"/PBMC_CCMarker_SPA.tsv"))
  CCM_SSA.df <- read.delim2(paste0(Import.Path,"/PBMC_CCMarker_SSA.tsv"))

  ## Modify the words
  library(stringr)
  CCM_SPA.df$Type <- str_replace(CCM_SPA.df$Type,"Cachexia.Marker.","")
  CCM_SPA.df$Type <- str_replace(CCM_SPA.df$Type,"Pos","EO")
  CCM_SPA.df$Type <- str_replace(CCM_SPA.df$Type,"Neg","LO")

  CCM_SSA.df$Type <- str_replace(CCM_SSA.df$Type,"Venn_Cachexia.Marker.","")
  CCM_SSA.df$Type <- str_replace(CCM_SSA.df$Type,"Pos","EO")
  CCM_SSA.df$Type <- str_replace(CCM_SSA.df$Type,"Neg","LO")

##### (Test) Specific Cell Type #####
  CCM_SPA_Mac.df <- CCM_SPA.df[grep("Mac", CCM_SPA.df$Type), ]

  CCM_SPA_Mac_EO.df <- CCM_SPA_Mac.df[grep("EO",CCM_SPA_Mac.df$Type), ]
  CCM_SPA_Mac_LO.df <- CCM_SPA_Mac.df[grep("LO",CCM_SPA_Mac.df$Type), ]

  ## Seperate genes ##
  CCM_SPA_Mac_EO_Split.df <- strsplit(as.character(
                                   CCM_SPA_Mac_EO.df$Marker), ", ")
  rm( CCM_SPA_Mac.df,CCM_SPA_Mac_EO.df,CCM_SPA_Mac_LO.df, CCM_SPA_Mac_EO_Split.df)

##### Count gene number #####
  source("FUN_CountGeneNum.R")

  ##------- SPA --------------------------------------------------------------
    CCM_SPA_CT.set <- CCM_SPA.df[,1]
    CCM_SPA_Split.lt <- strsplit(as.character(CCM_SPA.df$Marker), ", ")
    names(CCM_SPA_Split.lt) <- CCM_SPA_CT.set

    ## All
    CCM_SPA_Count.df <- CountGeneNum(CCM_SPA_Split.lt,mode=2)
    ## EO
    CCM_SPA_Split_EO.lt <- CCM_SPA_Split.lt[grep("EO", names(CCM_SPA_Split.lt))]
    CCM_SPA_Count_EO.df <- CountGeneNum(CCM_SPA_Split_EO.lt,mode=2)
    ## LO
    CCM_SPA_Split_LO.lt <- CCM_SPA_Split.lt[grep("LO", names(CCM_SPA_Split.lt))]
    CCM_SPA_Count_LO.df <- CountGeneNum(CCM_SPA_Split_LO.lt,mode=2)

    CCM_SPA_Count_Sum.df <- full_join(CCM_SPA_Count.df, CCM_SPA_Count_EO.df, by="Var1")
    CCM_SPA_Count_Sum.df <- full_join(CCM_SPA_Count_Sum.df, CCM_SPA_Count_LO.df, by="Var1")
    colnames(CCM_SPA_Count_Sum.df) <- c("Genes","Total.SPA","EO.SPA","LO.SPA")
    CCM_SPA_Count_Sum.df[is.na(CCM_SPA_Count_Sum.df)] <- 0

    ## For each cell type
    CCM_SPA_Split_EO_CT.df <- CountGeneNum(CCM_SPA_Split_EO.lt,mode=1)
    CCM_SPA_Split_LO_CT.df <- CountGeneNum(CCM_SPA_Split_LO.lt,mode=1)
    CCM_SPA_Count_Sum_CT.df <- full_join(CCM_SPA_Count_Sum.df,
                                                  CCM_SPA_Split_EO_CT.df,by="Genes")
    CCM_SPA_Count_Sum_CT.df <- full_join(CCM_SPA_Count_Sum_CT.df,
                                                  CCM_SPA_Split_LO_CT.df,by="Genes")
    CCM_SPA_Count_Sum_CT.df[is.na(CCM_SPA_Count_Sum_CT.df)] <- 0

    ##
    rm(CCM_SPA_Count.df, CCM_SPA_Split_EO.lt, CCM_SPA_Split_LO.lt,
       CCM_SPA_Count_EO.df,CCM_SPA_Count_LO.df,CCM_SPA_Split_EO_CT.df,CCM_SPA_Split_LO_CT.df)

  ##------- SSA --------------------------------------------------------------
  ## SSA_Intersect_All
    CCM_SSA_CT.set <- CCM_SSA.df[,1]
    CCM_SSA_Split.lt <- strsplit(as.character(CCM_SSA.df$Intersect), ", ")
    names(CCM_SSA_Split.lt) <- CCM_SSA_CT.set

    ## All
    CCM_SSA_Count.df <- CountGeneNum(CCM_SSA_Split.lt,mode=2)

    ## EO
    CCM_SSA_Split_EO.lt <- CCM_SSA_Split.lt[grep("EO", names(CCM_SSA_Split.lt))]
    CCM_SSA_Count_EO.df <- CountGeneNum(CCM_SSA_Split_EO.lt,mode=2)

    ## LO
    CCM_SSA_Split_LO.lt <- CCM_SSA_Split.lt[grep("LO", names(CCM_SSA_Split.lt))]
    CCM_SSA_Count_LO.df <- CountGeneNum(CCM_SSA_Split_LO.lt,mode=2)


    CCM_SSA_Count_Sum.df <- full_join(CCM_SSA_Count.df, CCM_SSA_Count_EO.df, by="Var1")
    CCM_SSA_Count_Sum.df <- full_join(CCM_SSA_Count_Sum.df, CCM_SSA_Count_LO.df, by="Var1")
    colnames(CCM_SSA_Count_Sum.df) <- c("Genes","Total.SSA.I","EO.SSA.I","LO.SSA.I")
    CCM_SSA_Count_Sum.df[is.na(CCM_SSA_Count_Sum.df)] <- 0

    ##
    rm(CCM_SSA_Count.df, CCM_SSA_Split_EO.lt, CCM_SSA_Split_LO.lt,
       CCM_SSA_Count_EO.df,CCM_SSA_Count_LO.df)

  ## SSA_Female_All
    CCM_SSA_F_CT.set <- CCM_SSA.df[,1]
    CCM_SSA_F_Split.lt <- strsplit(as.character(CCM_SSA.df$Female), ", ")
    names(CCM_SSA_F_Split.lt) <- CCM_SSA_F_CT.set

    ## All
    CCM_SSA_F_Count.df <- CountGeneNum(CCM_SSA_F_Split.lt,mode=2)

    ## EO
    CCM_SSA_F_Split_EO.lt <- CCM_SSA_F_Split.lt[grep("EO", names(CCM_SSA_F_Split.lt))]
    CCM_SSA_F_Count_EO.df <- CountGeneNum(CCM_SSA_F_Split_EO.lt,mode=2)

    ## LO
    CCM_SSA_F_Split_LO.lt <- CCM_SSA_F_Split.lt[grep("LO", names(CCM_SSA_F_Split.lt))]
    CCM_SSA_F_Count_LO.df <- CountGeneNum(CCM_SSA_F_Split_LO.lt,mode=2)


    CCM_SSA_F_Count_Sum.df <- full_join(CCM_SSA_F_Count.df, CCM_SSA_F_Count_EO.df, by="Var1")
    CCM_SSA_F_Count_Sum.df <- full_join(CCM_SSA_F_Count_Sum.df, CCM_SSA_F_Count_LO.df, by="Var1")
    colnames(CCM_SSA_F_Count_Sum.df) <- c("Genes","Total.SSA.F","EO.SSA.F","LO.SSA.F")
    CCM_SSA_F_Count_Sum.df[is.na(CCM_SSA_F_Count_Sum.df)] <- 0

    ##
    rm(CCM_SSA_F_Count.df, CCM_SSA_F_Split_EO.lt, CCM_SSA_F_Split_LO.lt,
       CCM_SSA_F_Count_EO.df,CCM_SSA_F_Count_LO.df)

  ## SSA_Male_All
    CCM_SSA_M_CT.set <- CCM_SSA.df[,1]
    CCM_SSA_M_Split.lt <- strsplit(as.character(CCM_SSA.df$Male), ", ")
    names(CCM_SSA_M_Split.lt) <- CCM_SSA_M_CT.set

    ## All
    CCM_SSA_M_Count.df <- CountGeneNum(CCM_SSA_M_Split.lt,mode=2)

    ## EO
    CCM_SSA_M_Split_EO.lt <- CCM_SSA_M_Split.lt[grep("EO", names(CCM_SSA_M_Split.lt))]
    CCM_SSA_M_Count_EO.df <- CountGeneNum(CCM_SSA_M_Split_EO.lt,mode=2)

    ## LO
    CCM_SSA_M_Split_LO.lt <- CCM_SSA_M_Split.lt[grep("LO", names(CCM_SSA_M_Split.lt))]
    CCM_SSA_M_Count_LO.df <- CountGeneNum(CCM_SSA_M_Split_LO.lt,mode=2)


    CCM_SSA_M_Count_Sum.df <- full_join(CCM_SSA_M_Count.df, CCM_SSA_M_Count_EO.df, by="Var1")
    CCM_SSA_M_Count_Sum.df <- full_join(CCM_SSA_M_Count_Sum.df, CCM_SSA_M_Count_LO.df, by="Var1")
    colnames(CCM_SSA_M_Count_Sum.df) <- c("Genes","Total.SSA.M","EO.SSA.M","LO.SSA.M")
    CCM_SSA_M_Count_Sum.df[is.na(CCM_SSA_M_Count_Sum.df)] <- 0

    ##
    rm(CCM_SSA_M_Count.df, CCM_SSA_M_Split_EO.lt, CCM_SSA_M_Split_LO.lt,
       CCM_SSA_M_Count_EO.df,CCM_SSA_M_Count_LO.df)


  ## Compare SPA and SSA
    CCM_Count_Sum.df <- full_join( CCM_SPA_Count_Sum.df,CCM_SSA_Count_Sum.df,by="Genes")
    CCM_Count_Sum.df <- full_join( CCM_Count_Sum.df,CCM_SSA_F_Count_Sum.df,by="Genes")
    CCM_Count_Sum.df <- full_join( CCM_Count_Sum.df,CCM_SSA_M_Count_Sum.df,by="Genes")
    CCM_Count_Sum.df[is.na(CCM_Count_Sum.df)] <- 0


##### Export the result #####

  write.csv( CCM_SPA_Count_Sum_CT.df ,
             file = paste0(Save.Path,"/",SampleType, "_CCMCount_SPA_CT.df.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )


  write.csv( CCM_Count_Sum.df ,
             file = paste0(Save.Path,"/",SampleType, "_CCMCount_Sum.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )


##### Visualization #####
  library(ggplot2)
  source("FUN_Beautify_ggplot.R")

  plot(x= seq(2:ncol(CCM_Count_Sum.df)), y= as.integer(CCM_Count_Sum.df[1,2:ncol(CCM_Count_Sum.df)]))

  ## EO
  # https://stackoverflow.com/questions/17423609/multiple-y-over-x-in-ggplot2
  CCM_Count_Sum.df2 <- CCM_Count_Sum.df[,c(1,3,6)] %>% arrange(desc(EO.SPA),desc(EO.SSA.I))

  library(reshape2)
  df2 <- melt(CCM_Count_Sum.df2[1:500,], id = "Genes")
  df2$Genes <- factor(df2$Genes,levels = df2$Genes[1:500])

  source("FUN_Beautify_LinePlot_Ori.R")
  EO1.P <- ggplot(df2, aes(x = Genes, y = value, color = variable,
                        group = variable ,linetype = variable))

  EO1.P <- EO1.P %>% LinePlot_Ori(.,AspRat=1,LegPos = c(0.85, 0.9),AxisTitleSize=2,
                        XtextSize=0,  YtextSize=25) +
                        scale_color_manual(values = c('#f7598b','#9755ab'))
  EO1.P
  ##
  CCMarker_EO1 <- CCM_Count_Sum.df2[CCM_Count_Sum.df2$EO.SPA==1 &CCM_Count_Sum.df2$EO.SSA.I ==1,]
  CCMarker_EO2 <- CCM_Count_Sum.df2[CCM_Count_Sum.df2$EO.SPA==2 &CCM_Count_Sum.df2$EO.SSA.I ==2,]
  CCMarker_EO3 <- CCM_Count_Sum.df2[CCM_Count_Sum.df2$EO.SPA==3 &CCM_Count_Sum.df2$EO.SSA.I ==3,]

  ## LO
  # https://stackoverflow.com/questions/17423609/multiple-y-over-x-in-ggplot2
  CCM_Count_Sum_LO.df <- CCM_Count_Sum.df[,c(1,4,7)] %>% arrange(desc(LO.SPA),desc(LO.SSA.I))

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
  CCM_Count_Sum.df4 <- CCM_Count_Sum.df[,c(1,6,9,12)] %>% arrange(desc(EO.SSA.I),desc(EO.SSA.F),desc(EO.SSA.M))

  library(reshape2)
  df4 <- melt(CCM_Count_Sum.df4[1:500,], id = "Genes")
  df4$Genes <- factor(df4$Genes,levels = df4$Genes[1:500])
  EOSS.P <- ggplot(df4, aes(x = Genes, y = value, color = variable,
                        group = variable ,linetype = variable))

  EOSS.P <- EOSS.P %>% LinePlot_Ori(.,AspRat=1,LegPos = c(0.85, 0.89),AxisTitleSize=2,
                        XtextSize=0,  YtextSize=25) +
    geom_point(shape = 12, size = 4, fill = "white") +
    scale_linetype_manual(#name="Pheno_Type",
      values=c("solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
      labels=c("EO.SSA.I","EO.SSA.F","EO.SSA.M")) +
    scale_color_manual(values = c('#9755ab','#ff75da','#6db5fc'))

  CCM_Count_Sum.df4 %>% arrange(desc(EO.SSA.F))
  EOSS.P


  ## LO
  CCM_Count_Sum_LO_Sex.df <- CCM_Count_Sum.df[,c(1,7,10,13)] %>% arrange(desc(LO.SSA.I),desc(LO.SSA.F),desc(LO.SSA.M))

  library(reshape2)
  df4 <- melt(CCM_Count_Sum_LO_Sex.df[1:500,], id = "Genes")
  df4$Genes <- factor(df4$Genes,levels = df4$Genes[1:500])
  LOSS.P <- ggplot(df4, aes(x = Genes, y = value, color = variable,
                        group = variable ,linetype = variable))

  LOSS.P <- LOSS.P %>% LinePlot_Ori(.,AspRat=1,LegPos = c(0.85, 0.998),AxisTitleSize=2,
                        XtextSize=0,  YtextSize=25)+
    geom_point(shape = 12, size = 4, fill = "white") +
    scale_linetype_manual(#name="Pheno_Type",
      values=c("solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
      labels=c("EO.SSA.I","EO.SSA.F","EO.SSA.M")) +
    scale_color_manual(values = c('#9755ab','#ff75da','#6db5fc'))

  CCM_Count_Sum_LO_Sex.df %>% arrange(desc(LO.SSA.F))
  LOSS.P

##### Gene count vs Cluster #####
  # GF: Gene Frequency
  Base_GF <- data.frame(Var1=seq(from = 0, to = 9))
  Base_GF$Var1 <- as.factor(Base_GF$Var1)
  EO_SPA_GF <- data.frame(table(CCM_Count_Sum.df$EO.SPA))
  colnames(EO_SPA_GF)[2] <-  'EO.SPA'
  EO_SPA_GF <- full_join(Base_GF, EO_SPA_GF)

  LO_SPA_GF <- data.frame(table(CCM_Count_Sum.df$LO.SPA))
  colnames(LO_SPA_GF)[2] <-  'LO.SPA'

  Sum_GF <- full_join(EO_SPA_GF, LO_SPA_GF)

  EO_SSAI_GF <- data.frame(table(CCM_Count_Sum.df$EO.SSA.I))
  colnames(EO_SSAI_GF)[2] <-  'EO.SSAI'

  Sum_GF <- full_join(Sum_GF, EO_SSAI_GF)

  LO_SSAI_GF <- data.frame(table(CCM_Count_Sum.df$LO.SSA.I))
  colnames(LO_SSAI_GF)[2] <-  'LO.SSAI'

  Sum_GF <- full_join(Sum_GF, LO_SSAI_GF)

  EO_SSAF_GF <- data.frame(table(CCM_Count_Sum.df$EO.SSA.F))
  colnames(EO_SSAF_GF)[2] <-  'EO.SSAF'

  Sum_GF <- full_join(Sum_GF, EO_SSAF_GF)

  LO_SSAF_GF <- data.frame(table(CCM_Count_Sum.df$LO.SSA.F))
  colnames(LO_SSAF_GF)[2] <-  'LO.SSAF'

  Sum_GF <- full_join(Sum_GF, LO_SSAF_GF)

  EO_SSAM_GF <- data.frame(table(CCM_Count_Sum.df$EO.SSA.M))
  colnames(EO_SSAM_GF)[2] <-  'EO.SSAM'

  Sum_GF <- full_join(Sum_GF, EO_SSAM_GF)

  LO_SSAM_GF <- data.frame(table(CCM_Count_Sum.df$LO.SSA.M))
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

  Sum_GF.P <- Sum_GF.P %>% LinePlot(.,AspRat=1,LegPos = c(0.85, 0.9),AxisTitleSize=2,
                              XtextSize=15,  YtextSize=15,
                              scale_y_c =  seq(0, max(Sum_GF_Melt$Gene_Count)+30, 50))+
                              scale_linetype_manual(#name="Pheno_Type",
                                values=c("solid", "solid", "solid", "solid",
                                         "dotted", "dotdash","dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                                labels=c("EO.SPA","LO.SPA","EO.SSAI","LO.SSAI",
                                         "EO.SSAF","LO.SSAF","EO.SSAM","LO.SSAM")) +
                              scale_color_manual(values = c('#d736ff','#874299','#e79bfa','#bc87c9',
                                                            '#ff8fea','#c995c0','#3d84ff','#7099e0'))
  Sum_GF.P

  ## Log
  # Line Color Setting
    LTypeSet.lt <- list(
      values = c("solid", "solid", "solid", "solid", "dotted", "dotdash","dotted", "dotdash"), # values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
      labels = c("EO.SPA","LO.SPA","EO.SSAI","LO.SSAI","EO.SSAF","LO.SSAF","EO.SSAM","LO.SSAM")
    )
    LColorSet.lt <- list(
      values = c('#c79a6d','#784f24','#b372c2','#874299','#ff8fea','#bd559c','#7099e0','#151e75')
    )

    Sum_GF_Melt$`Log2(Gene_Count+1)` <- log2(Sum_GF_Melt$Gene_Count+1)
    Sum_GF_Melt$Cluster_Count <- as.numeric(Sum_GF_Melt$Cluster_Count)

    Sum_GF_log.P <- ggplot(Sum_GF_Melt, aes(x = Cluster_Count, y = `Log2(Gene_Count+1)`, color = Type,
                                        group = Type ,linetype = Type))

    Sum_GF_log.P <- Sum_GF_log.P %>% LinePlot(., AspRat=0.5, LegPos = c(0.9, 0.75), AxisTitleSize=1.5,
                                      XtextSize=20,  YtextSize=20,
                                      scale_x_c =  seq(0, max(Sum_GF_Melt$Cluster_Count)),
                                      scale_y_c =  seq(0, max(Sum_GF_Melt$Gene_Count)))+
      scale_linetype_manual(#name="Pheno_Type",
        values= LTypeSet.lt$values,
        labels= LTypeSet.lt$labels) +
      scale_color_manual(values= LColorSet.lt$values) +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) # Remove grid

    Sum_GF_log.P



##### Heatmap #####
  library(ComplexHeatmap)
  H_CCM_SPA_CT.df <- CCM_SPA_Count_Sum_CT.df
  H_CCM_SPA_CT.df <- left_join(H_CCM_SPA_CT.df,CCM_Count_Sum.df)
  row.names(H_CCM_SPA_CT.df) <- H_CCM_SPA_CT.df[,1]
  H_CCM_SPA_CT.df <- H_CCM_SPA_CT.df[,-1]

  H_CCM_SPA_CT.df <- as.data.frame(t(H_CCM_SPA_CT.df))
  Heatmap(H_CCM_SPA_CT.df)

  CCM_TOP_EO.set <- c(as.character(CCMarker_EO1$Genes),as.character(CCMarker_EO2$Genes),
                    as.character(CCMarker_EO3$Genes),"Chil3","S100a8")
  CCM_TOP_LO.set <- c(as.character(CCMarker_LO1$Genes),as.character(CCMarker_LO2$Genes),
                     as.character(CCMarker_LO3$Genes))

  CCM_TOP.set <- c(CCM_TOP_EO.set,CCM_TOP_LO.set)

  CCM_Count_Sum_T_S.df <- H_CCM_SPA_CT.df[row.names(H_CCM_SPA_CT.df) %in% CCM_TOP.set,]

  library(circlize)
  col_fun = colorRamp2(c(0, 6, 12), c("white", "#B15BFF", "#5B00AE"))
  col_fun(seq(-3, 3))

  Heatmap(CCM_Count_Sum_T_S.df[,-1:-3], name = "Num", col = col_fun)
  Heatmap(CCM_Count_Sum_T_S.df, name = "Num", col = col_fun)

##### Export PDF #####
  pdf(file = paste0(Save.Path,"/PBMC_GeneCount.pdf"),width = 7, height = 7 )
    EO1.P
    LO1.P
    EOSS.P
    LOSS.P
  dev.off()

  pdf(file = paste0(Save.Path,"/PBMC_GeneCount_Total.pdf"),width = 14, height = 7 )
    Sum_GF_log.P
  dev.off()

# ##### Redefine the Male & Female Count #####
#   CCM_Count_Sum_ReDe.df <- CCM_Count_Sum.df
#   CCM_Count_Sum_ReDe.df$EO.SSA.F <- CCM_Count_Sum_ReDe.df$EO.SSA.F+CCM_Count_Sum_ReDe.df$EO.SSA.I
#   CCM_Count_Sum_ReDe.df$EO.SSA.M <- CCM_Count_Sum_ReDe.df$EO.SSA.M+CCM_Count_Sum_ReDe.df$EO.SSA.I
#
#   CCM_Count_Sum_ReDe.df$LO.SSA.F <- CCM_Count_Sum_ReDe.df$LO.SSA.F+CCM_Count_Sum_ReDe.df$LO.SSA.I
#   CCM_Count_Sum_ReDe.df$LO.SSA.M <- CCM_Count_Sum_ReDe.df$LO.SSA.M+CCM_Count_Sum_ReDe.df$LO.SSA.I
#
#   CCM_Count_Sum_ReDe.df$Total.SSA.F <- CCM_Count_Sum_ReDe.df$Total.SSA.F+CCM_Count_Sum_ReDe.df$Total.SSA.I
#   CCM_Count_Sum_ReDe.df$Total.SSA.M <- CCM_Count_Sum_ReDe.df$Total.SSA.M+CCM_Count_Sum_ReDe.df$Total.SSA.I
#
#   # https://stackoverflow.com/questions/17423609/multiple-y-over-x-in-ggplot2
#   CCM_Count_Sum_ReDe.df3 <- CCM_Count_Sum_ReDe.df[,c(1,6,9,12)] %>% arrange(desc(EO.SSA.I),desc(EO.SSA.F),desc(EO.SSA.M))
#
#   library(reshape2)
#   df3 <- melt(CCM_Count_Sum_ReDe.df3[1:50,], id = "Genes")
#   df3$Genes <- factor(df3$Genes,levels = df3$Genes[1:50])
#   P1 <- ggplot(df3, aes(x = Genes, y = value, color = variable,
#                         group = variable ,linetype = variable)) +
#     geom_line(size = 1) + geom_point(size = 3) +
#     scale_color_brewer(palette = "Set1")+
#     geom_point(shape = 12, size = 4, fill = "white")+
#     theme_set(theme_bw())+ # Remove the background
#     theme(panel.grid.major=element_line(colour=NA))  # Remove the grid
#
#   P1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.95, 0.85),AxisTitleSize=1.5,
#                         XtextSize=15,  YtextSize=15) +
#     geom_point(shape = 12, size = 4, fill = "white") +
#     scale_linetype_manual(#name="Pheno_Type",
#                           values=c("solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
#                           labels=c("EO.SSA.I","EO.SSA.F","EO.SSA.M")) +
#     scale_color_manual(values = c('#7903ab','#ff75da','#6db5fc'))+
#     geom_point(shape = 12, size = 4, fill = "white")
#
#
#
#
#
#   CellNum_P3 <- ggplot(Pheno.mt.table_celltype3, aes(x = factor(Cell_Type), y = Number,
#                                                      colour = Pheno_Type,
#                                                      group = Pheno_Type,linetype=Pheno_Type
#   )) +
#     geom_line(size=1.5) +
#     scale_linetype_manual(name="Pheno_Type",
#                           values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
#                           labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")) +
#     scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
#     geom_point(shape = 12, size = 4, fill = "white")


