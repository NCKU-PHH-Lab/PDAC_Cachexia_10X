##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load libray #####
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(stringr)

##### Path and folder setting #####
  Import.Path = paste0(getwd(),"/2022-08-13_PBMC_Main")
  Save.Path = paste0(getwd(),"/2022-08-13_PBMC_Main")
  dir.create(Save.Path)
  SampleType = "PBMC"

  RawData.Path = ""

##### Load data #####
  CCM_SPA.df <- read.delim2(paste0(Import.Path,"/KLC_PBMC_CCMCount_SPA_CT_Num.csv"),
                            sep = ",")

  CCM_SSA.df <- read.delim2(paste0(Import.Path,"/KLC_PBMC_CCMCount_SSA_CT_Num.csv"),
                            sep = ",")


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

  # TTT <- CCM_All.df
  # ## Try
  # TTT2 <- TTT %>% mutate(T_EO_Match = CCM_All.df$T_EO ==1 & CCM_All.df$T_EO_I ==1)
  # sum(TTT2$T_EO_Match)
  # TTT2$T_EO_Match
  #
  # ## Try2
  # TTT2 <- TTT %>% mutate("WTF2" = CCM_All.df$"T_EO" ==1 & CCM_All.df$"T_EO_I" ==1)
  # TTT2$WTF2
  # sum(TTT2$WTF2)
  #
  # Genelist <- TTT2[TTT2$WTF2==T,1]
  #
  # library(stringr)
  # Genelist2 <- str_c(Genelist, collapse = ", ")
  # #Genelist %>% as.data.frame() -> Genelist

  ##
    PhenoType_SPA.Set <- colnames(CCM_All.df)[4:27]
    ncol_Ori <- ncol(CCM_All.df)
    CCM_All_M.df <- CCM_All.df

    for (i in 1:length(PhenoType_SPA.Set)) {
      CCM_All_M.df <- CCM_All_M.df %>% mutate(Temp =
                                      abs(CCM_All_M.df[PhenoType_SPA.Set[i]]) == 1 &
                                      abs(CCM_All_M.df[paste0(PhenoType_SPA.Set[i],"_I")]) == 1)
      CCM_All_M.df$Temp <- as.numeric(CCM_All_M.df$Temp)
      colnames(CCM_All_M.df)[ncol_Ori+i] <- paste0(PhenoType_SPA.Set[i],"_Ma")
    }
    rm(i)

    # Check
    # t(CCM_All_M.df) %>% View()

  ## Summary
    CCM_All_M.df$Match_Sum <- rowSums(CCM_All_M.df[,107:130])
    CCM_All_M.df$Match_EO_Sum <- rowSums(CCM_All_M.df[,107:118])
    CCM_All_M.df$Match_LO_Sum <- rowSums(CCM_All_M.df[,119:130])
    CCM_All_M_S.df <- CCM_All_M.df[,c(1,107:133)] %>% arrange(desc(Match_Sum))
    CCM_All_M_S.df <- CCM_All_M_S.df[!CCM_All_M_S.df$Match_Sum == 0,]


    CCM_All_M_S_LO.df <- CCM_All_M_S.df[,grep("_LO",colnames(CCM_All_M_S.df))]
    CCM_All_M_S_LO.df <- -CCM_All_M_S_LO.df
    CCM_All_M_S.df[,c(14:25,28)] <- CCM_All_M_S_LO.df
    rm(CCM_All_M_S_LO.df)

  ## Plot Heatmap
    minC <- min(CCM_All_M_S.df$Match_LO_Sum)
    maxC <- max(CCM_All_M_S.df$Match_EO_Sum)
    col_fun_Match = colorRamp2(c(minC,minC/2, 0,maxC/2,maxC),
                         c("#003bbd","#5183f0","white", "#f051c2", "#bd0087"))
    HeatmapTry <- Heatmap(CCM_All_M_S.df[,-1], name = "Num", col = col_fun_Match)
    HeatmapTry
    library(eoffice)
    topptx(HeatmapTry,"Temp.pptx")
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

    library(ggplot2)
    CCM_Level_Num.p <- ggplot(data=CCM_Level_Num.df, aes(x=Type, y=GeneNum)) +
         geom_bar(stat="identity")
    CCM_Level_Num.p

    rm(CCM_EO_1.set,CCM_EO_2.set,CCM_LO_1.set,CCM_LO_2.set,
       CCM_EO_1.str,CCM_EO_2.str,CCM_LO_1.str,CCM_LO_2.str)

    ##### Grab candidate gene by PhenoType #####
    ##
    PhenoType_Ma.Set <- colnames(CCM_All_M_S.df)
    PhenoType_Ma.Set <- PhenoType_Ma.Set[-c(1,26:28)]

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

    require(ggplot2)
    CCM_CT_Num.p <- ggplot(data=CCM_CT_Num.df, aes(x=Type, y=GeneNum)) +
      geom_bar(stat="identity")+
      #theme(axis.title.x = element_text(angle = 90,face="italic",colour = "darkred",size=14)) +
    theme(axis.text.x = element_text(face="bold",  size = 12,angle = 90, hjust = 1, vjust = .5)) # Change the size along the x axis

    CCM_CT_Num.p

    rm(CCM_CT_Num_Temp.df)

    CCM_Num.df <- rbind(CCM_Level_Num.df, CCM_CT_Num.df)
    CCM.df <- rbind(CCM_Level.df, CCM_CT.df)

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
  CCM_All_M2.df <- full_join(CCM_All.df,CCM_All_M_S.df,by="Genes")
  CCM_All_M2.df[is.na(CCM_All_M2.df)] <- 0

##### Export the result #####
  write.csv( CCM_All_M_S.df ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMHeatmap_SPA_SSA_Match_S.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )

  write.csv( CCM_All_M2.df ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCMHeatmap_SPA_SSA_Match_All.csv"),
             #sep = ",",
             quote = F,
             row.names = F
  )


  write.csv( CCM_Num.df ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCM_Num.csv"),
             #sep = ",",
             quote = F,
             row.names = T
  )

  write.csv( CCM.df ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCM.csv"),
             #sep = ",",
             quote = F,
             row.names = T
  )

  write.table( CCM.df ,
             file = paste0(Save.Path,"/KLC_",SampleType, "_CCM.tsv"),
             sep = "\t",
             quote = F,
             row.names = T
  )
