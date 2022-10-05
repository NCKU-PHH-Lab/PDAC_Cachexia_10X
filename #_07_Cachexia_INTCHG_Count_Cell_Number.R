## INTCHG: Interchangeable
  # For PBMC
  scRNA.SeuObj <- PBMC.combined
  SampleType = "PBMC"

  ## For SC
  # scRNA.SeuObj <- SC.combined
  # SampleType = "SC"

##### 07 Count Cell number  #####
  Pheno.df <- data.frame(sample = scRNA.SeuObj@meta.data[["sample"]],celltype = scRNA.SeuObj@meta.data[["celltype"]],
                         Cachexia = scRNA.SeuObj@meta.data[["Cachexia"]],Sex = scRNA.SeuObj@meta.data[["Sex"]])

  Freq_sample.df <- table(Pheno.df$sample) %>% as.data.frame()
  Freq_CT.df <- table(Pheno.df$celltype) %>% as.data.frame()
  Freq_Cach.df <- table(Pheno.df$Cachexia) %>% as.data.frame()
  Freq_Sex.df <- table(Pheno.df$Sex) %>% as.data.frame()

  #
  Pheno_EO_M.df <- Pheno.df[Pheno.df$sample=="EO.M",]
  Pheno_LO_M.df <- Pheno.df[Pheno.df$sample=="LO.M",]
  Pheno_EO_F.df <- Pheno.df[Pheno.df$sample=="EO.F",]
  Pheno_LO_F.df <- Pheno.df[Pheno.df$sample=="LO.F",]

  # Count EO_M_CT
  Freq_EO_M_CT.df <- table(Pheno_EO_M.df$celltype) %>% as.data.frame()
  Freq_EO_M_CT.df <- data.frame(Type="EO.M",Freq_EO_M_CT.df)
  Freq_EO_M_CT.df$Percent <- Freq_EO_M_CT.df$Freq/sum(Freq_EO_M_CT.df$Freq)

  # Count LO_M_CT
  Freq_LO_M_CT.df <- table(Pheno_LO_M.df$celltype) %>% as.data.frame()
  Freq_LO_M_CT.df <- data.frame(Type="LO.M",Freq_LO_M_CT.df)
  Freq_LO_M_CT.df$Percent <- Freq_LO_M_CT.df$Freq/sum(Freq_LO_M_CT.df$Freq)

  # Count EO_F_CT
  Freq_EO_F_CT.df <- table(Pheno_EO_F.df$celltype) %>% as.data.frame()
  Freq_EO_F_CT.df <- data.frame(Type="EO.F",Freq_EO_F_CT.df)
  Freq_EO_F_CT.df$Percent <- Freq_EO_F_CT.df$Freq/sum(Freq_EO_F_CT.df$Freq)

  # Count LO_F_CT
  Freq_LO_F_CT.df <- table(Pheno_LO_F.df$celltype) %>% as.data.frame()
  Freq_LO_F_CT.df <- data.frame(Type="LO.F",Freq_LO_F_CT.df)
  Freq_LO_F_CT.df$Percent <- Freq_LO_F_CT.df$Freq/sum(Freq_LO_F_CT.df$Freq)

  # Combind all count of sample
  Freq_All.df <- rbind(Freq_EO_M_CT.df,Freq_LO_M_CT.df,
                       Freq_EO_F_CT.df,Freq_LO_F_CT.df)
  Freq_All.df <- data.frame(Index = row.names(Freq_All.df),Freq_All.df )
  colnames(Freq_All.df) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")
  # Freq_All.df$Index <- factor(Freq_All.df$Index,
  #                                              levels = Freq_All.df$Index)

#### LinePlot ####
  # https://ithelp.ithome.com.tw/articles/10186047
  # Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
  #                                    levels = sort(unique(as.character(Freq_All.df$Cell_Type))))

  Freq_All.df$Cell_Type <- factor(Freq_All.df$Cell_Type,
                                  levels = Cell_Type_Order.set)

  CellNum_P1 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Number,
                                        colour = Pheno_Type, group = Pheno_Type)) +
    geom_line(linetype = "dashed",size=1.5) +
    geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  CellNum_P1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=15,  YtextSize=15, xangle = 90,
                                LegTextSize = 15) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P1
  CellNum_P1

  CellNum_P2 <- ggplot(Freq_All.df, aes(x = factor(Cell_Type), y = Percent,
                                        colour = Pheno_Type, group = Pheno_Type)) +
    geom_line(linetype = "dashed",size=1.5) +
    geom_point(shape = 12, size = 4, fill = "white")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  CellNum_P2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=18, xangle = 90,
                                LegTextSize = 15) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P2
  CellNum_P2


#### All type compare to Combine Sex ####
  ##
  Pheno_EO.df <- Pheno.df[Pheno.df$Cachexia=="EO",]
  Pheno_LO.df <- Pheno.df[Pheno.df$Cachexia=="LO",]

  # Count EO_CT
  Freq_EO_CT.df <- table(Pheno_EO.df$celltype) %>% as.data.frame()
  Freq_EO_CT.df <- data.frame(Type="EO",Freq_EO_CT.df)
  Freq_EO_CT.df$Percent <- Freq_EO_CT.df$Freq/sum(Freq_EO_CT.df$Freq)

  # Count LO_CT
  Freq_LO_CT.df <- table(Pheno_LO.df$celltype) %>% as.data.frame()
  Freq_LO_CT.df <- data.frame(Type="LO",Freq_LO_CT.df)
  Freq_LO_CT.df$Percent <- Freq_LO_CT.df$Freq/sum(Freq_LO_CT.df$Freq)

  # Combind all count of sample
  Freq_All_Ca.df <- rbind(Freq_EO_M_CT.df,Freq_LO_M_CT.df,
                          Freq_EO_F_CT.df,Freq_LO_F_CT.df,
                          Freq_EO_CT.df,Freq_LO_CT.df)

  Freq_All_Ca.df <- data.frame(Index = row.names(Freq_All_Ca.df),Freq_All_Ca.df )
  colnames(Freq_All_Ca.df) <- c("Index","Pheno_Type","Cell_Type","Number","Percent")

  # Change the order
  # https://blog.csdn.net/weixin_48172266/article/details/117537465
  # CTOrder.set <- factor(Freq_All_Ca.df$Cell_Type,
  #                       levels = sort(unique(as.character(Freq_All_Ca.df$Cell_Type))))
  #
  # Freq_All_Ca.df <- Freq_All_Ca.df %>%
  #                   mutate(Cell_Type = CTOrder.set)

  Freq_All_Ca.df$Percent <- as.numeric(Freq_All_Ca.df$Percent)

  write.table( Freq_All_Ca.df ,
               file = paste0(Save.Path,"/",SampleType,"_CellCount_CT_Ca.tsv"),
               sep = "\t",
               quote = F,
               row.names = F
  )

  # https://stackoverflow.com/questions/27350243/ggplot-line-graph-with-different-line-styles-and-markers/27350366
  # https://www.coder.work/article/6971741
  # https://stackoverflow.com/questions/11344561/controlling-line-color-and-line-type-in-ggplot-legend

  #### LinePlot ####
  Freq_All_Ca.df$Cell_Type <- factor(Freq_All_Ca.df$Cell_Type,
                                     levels = Cell_Type_Order.set)

  CellNum_P3 <- ggplot(Freq_All_Ca.df, aes(x = factor(Cell_Type), y = Number,
                                           colour = Pheno_Type,
                                           group = Pheno_Type,linetype=Pheno_Type
  )) +
    geom_line(size=1.5) +
    scale_linetype_manual(name="Pheno_Type",
                          values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                          labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")) +
    scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
    geom_point(shape = 12, size = 4, fill = "white") + theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


  CellNum_P3 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.4, 0.8),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=18,xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P3
  CellNum_P3

  library(eoffice)
  topptx(CellNum_P3,paste0(Save.Path,"/Temp.pptx"))

  CellNum_P4 <- ggplot(Freq_All_Ca.df, aes(x = factor(Cell_Type), y = Percent,
                                           colour = Pheno_Type,
                                           group = Pheno_Type,linetype=Pheno_Type
  )) +
    geom_line(size=1.5) +
    scale_linetype_manual(name="Pheno_Type",
                          values=c("solid",  "dotted","dotdash", "solid", "dotted", "dotdash"), #  values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
                          labels=c("EO","EO.F","EO.M","LO","LO.F","LO.M")) +
    scale_color_manual(values = c('#ba0449','#ff52bd','#f0679b','#3d3c99','#5292f2','#33aef5'))+
    geom_point(shape = 12, size = 4, fill = "white") +
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
    #theme_set(theme_bw())+ # Remove the background
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  CellNum_P4 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.15, 0.82),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> CellNum_P4
  CellNum_P4

##### Clean the dataframe #####
  CellPheno.lt <- list(Pheno.df, Pheno_EO.df, Pheno_LO.df,
                       Pheno_EO_M.df,Pheno_LO_M.df,Pheno_EO_F.df,Pheno_LO_F.df)
  names(CellPheno.lt) <- c("Pheno.df", "Pheno_EO.df", "Pheno_LO.df",
                           "Pheno_EO_M.df","Pheno_LO_M.df","Pheno_EO_F.df","Pheno_LO_F.df")
  rm(Pheno_EO.df, Pheno_LO.df, Pheno_EO_M.df,Pheno_LO_M.df,Pheno_EO_F.df,Pheno_LO_F.df)

  CellFreq.lt <- list(Freq_All_Ca.df, Freq_All.df, Freq_sample.df,
                      Freq_Sex.df,Freq_Cach.df, Freq_CT.df,
                      Freq_EO_CT.df, Freq_EO_F_CT.df, Freq_EO_M_CT.df,
                      Freq_LO_CT.df, Freq_LO_F_CT.df, Freq_LO_M_CT.df)
  names(CellFreq.lt) <- c("Freq_All_Ca.df", "Freq_All.df", "Freq_sample.df",
                          "Freq_Sex.df", "Freq_Cach.df", "Freq_CT.df",
                          "Freq_EO_CT.df", "Freq_EO_F_CT.df", "Freq_EO_M_CT.df",
                          "Freq_LO_CT.df", "Freq_LO_F_CT.df", "Freq_LO_M_CT.df")
  rm(Freq_sample.df,
     Freq_Sex.df,Freq_Cach.df, Freq_CT.df,
     Freq_EO_CT.df, Freq_EO_F_CT.df, Freq_EO_M_CT.df,
     Freq_LO_CT.df, Freq_LO_F_CT.df, Freq_LO_M_CT.df)

##### BarPlot #####
  # https://blog.gtwang.org/r/ggplot2-tutorial-layer-by-layer-plotting/3/
  colnames(Pheno.df) <- c("sample","Cell_Type","Cachexia","Sex")
  Pheno.df$Cell_Type <- factor(Pheno.df$Cell_Type,
                               levels = sort(unique(as.character(Pheno.df$Cell_Type))))

  # sample
  BarPlot1_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=sample)) +
    geom_bar(position="dodge")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot1_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot1_1
  BarPlot1_1

  BarPlot1_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=sample)) +
    geom_bar(position="fill")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot1_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot1_2
  BarPlot1_2

  # Cachexia
  BarPlot2_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=Cachexia)) +
    geom_bar(position="dodge")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot2_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot2_1
  BarPlot2_1

  BarPlot2_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=Cachexia)) +
    geom_bar(position="fill")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot2_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot2_2
  BarPlot2_2

  # Sex
  BarPlot3_1 <- ggplot(Pheno.df, aes(Cell_Type, fill=Sex)) +
    geom_bar(position="dodge")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot3_1 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot3_1
  BarPlot3_1
  BarPlot3_2 <- ggplot(Pheno.df, aes(Cell_Type, fill=Sex)) +
    geom_bar(position="fill")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  BarPlot3_2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) -> BarPlot3_2

  BarPlot3_2


##### Export PDF file #####
  pdf(file = paste0(Save.Path,"/",SampleType,"_CellCount_LinePlot.pdf"),
      width = 7, height = 7 )
    CellNum_P4
    CellNum_P3
    CellNum_P1
    CellNum_P2
    BarPlot1_1
    BarPlot1_2
    BarPlot2_1
    BarPlot2_2
    BarPlot3_1
    BarPlot3_2
  dev.off() # graphics.off()

  rm(CellNum_P1, CellNum_P2, CellNum_P3, CellNum_P4, BarPlot1_1, BarPlot1_2,
     BarPlot2_1, BarPlot2_2, BarPlot3_1, BarPlot3_2)

#### Save RData ####
  save.image(paste0(Save.Path,"/07_Count_Cell_number.RData"))
