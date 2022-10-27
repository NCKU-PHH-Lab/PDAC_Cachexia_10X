## INTCHG: Interchangeable
## SubType Setting
  if(SampleType == "PBMC"){
    # For PBMC
    scRNA.SeuObj <- PBMC.combined

  }else if(SampleType == "SC"){
    # For SC
    scRNA.SeuObj <- SC.combined

  }

  library(ggpubr)
##### 07 Count Cell number  #####
  Pheno.df <- data.frame(sample = scRNA.SeuObj@meta.data[["sample"]],celltype = scRNA.SeuObj@meta.data[["celltype"]],
                         Cachexia = scRNA.SeuObj@meta.data[["Cachexia"]],Sex = scRNA.SeuObj@meta.data[["Sex"]])

  Freq_sample.df <- table(Pheno.df$sample) %>% as.data.frame()
  Freq_CT.df <- table(Pheno.df$celltype) %>% as.data.frame()
  Freq_Cach.df <- table(Pheno.df$Cachexia) %>% as.data.frame()
  Freq_Sex.df <- table(Pheno.df$Sex) %>% as.data.frame()

  #
  Pheno_EOCX_M.df <- Pheno.df[Pheno.df$sample=="EOCX.M",]
  Pheno_PreCX_M.df <- Pheno.df[Pheno.df$sample=="PreCX.M",]
  Pheno_EOCX_F.df <- Pheno.df[Pheno.df$sample=="EOCX.F",]
  Pheno_PreCX_F.df <- Pheno.df[Pheno.df$sample=="PreCX.F",]

  # Count EOCX_M_CT
  Freq_EOCX_M_CT.df <- table(Pheno_EOCX_M.df$celltype) %>% as.data.frame()
  Freq_EOCX_M_CT.df <- data.frame(Type="EOCX.M",Freq_EOCX_M_CT.df)
  Freq_EOCX_M_CT.df$Percent <- Freq_EOCX_M_CT.df$Freq/sum(Freq_EOCX_M_CT.df$Freq)

  # Count PreCX_M_CT
  Freq_PreCX_M_CT.df <- table(Pheno_PreCX_M.df$celltype) %>% as.data.frame()
  Freq_PreCX_M_CT.df <- data.frame(Type="PreCX.M",Freq_PreCX_M_CT.df)
  Freq_PreCX_M_CT.df$Percent <- Freq_PreCX_M_CT.df$Freq/sum(Freq_PreCX_M_CT.df$Freq)

  # Count EOCX_F_CT
  Freq_EOCX_F_CT.df <- table(Pheno_EOCX_F.df$celltype) %>% as.data.frame()
  Freq_EOCX_F_CT.df <- data.frame(Type="EOCX.F",Freq_EOCX_F_CT.df)
  Freq_EOCX_F_CT.df$Percent <- Freq_EOCX_F_CT.df$Freq/sum(Freq_EOCX_F_CT.df$Freq)

  # Count PreCX_F_CT
  Freq_PreCX_F_CT.df <- table(Pheno_PreCX_F.df$celltype) %>% as.data.frame()
  Freq_PreCX_F_CT.df <- data.frame(Type="PreCX.F",Freq_PreCX_F_CT.df)
  Freq_PreCX_F_CT.df$Percent <- Freq_PreCX_F_CT.df$Freq/sum(Freq_PreCX_F_CT.df$Freq)

  # Combind all count of sample
  Freq_All.df <- rbind(Freq_EOCX_M_CT.df,Freq_PreCX_M_CT.df,
                       Freq_EOCX_F_CT.df,Freq_PreCX_F_CT.df)
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
  Pheno_EOCX.df <- Pheno.df[Pheno.df$Cachexia=="EOCX",]
  Pheno_PreCX.df <- Pheno.df[Pheno.df$Cachexia=="PreCX",]

  # Count EOCX_CT
  Freq_EOCX_CT.df <- table(Pheno_EOCX.df$celltype) %>% as.data.frame()
  Freq_EOCX_CT.df <- data.frame(Type="EOCX",Freq_EOCX_CT.df)
  Freq_EOCX_CT.df$Percent <- Freq_EOCX_CT.df$Freq/sum(Freq_EOCX_CT.df$Freq)

  # Count PreCX_CT
  Freq_PreCX_CT.df <- table(Pheno_PreCX.df$celltype) %>% as.data.frame()
  Freq_PreCX_CT.df <- data.frame(Type="PreCX",Freq_PreCX_CT.df)
  Freq_PreCX_CT.df$Percent <- Freq_PreCX_CT.df$Freq/sum(Freq_PreCX_CT.df$Freq)

  # Combind all count of sample
  Freq_All_Ca.df <- rbind(Freq_EOCX_M_CT.df,Freq_PreCX_M_CT.df,
                          Freq_EOCX_F_CT.df,Freq_PreCX_F_CT.df,
                          Freq_EOCX_CT.df,Freq_PreCX_CT.df)

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
                          labels=c("EOCX","EOCX.F","EOCX.M","PreCX","PreCX.F","PreCX.M")) +
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
                          labels=c("EOCX","EOCX.F","EOCX.M","PreCX","PreCX.F","PreCX.M")) +
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
  CellPheno.lt <- list(Pheno.df, Pheno_EOCX.df, Pheno_PreCX.df,
                       Pheno_EOCX_M.df,Pheno_PreCX_M.df,Pheno_EOCX_F.df,Pheno_PreCX_F.df)
  names(CellPheno.lt) <- c("Pheno.df", "Pheno_EOCX.df", "Pheno_PreCX.df",
                           "Pheno_EOCX_M.df","Pheno_PreCX_M.df","Pheno_EOCX_F.df","Pheno_PreCX_F.df")
  rm(Pheno_EOCX.df, Pheno_PreCX.df, Pheno_EOCX_M.df,Pheno_PreCX_M.df,Pheno_EOCX_F.df,Pheno_PreCX_F.df)

  CellFreq.lt <- list(Freq_All_Ca.df, Freq_All.df, Freq_sample.df,
                      Freq_Sex.df,Freq_Cach.df, Freq_CT.df,
                      Freq_EOCX_CT.df, Freq_EOCX_F_CT.df, Freq_EOCX_M_CT.df,
                      Freq_PreCX_CT.df, Freq_PreCX_F_CT.df, Freq_PreCX_M_CT.df)
  names(CellFreq.lt) <- c("Freq_All_Ca.df", "Freq_All.df", "Freq_sample.df",
                          "Freq_Sex.df", "Freq_Cach.df", "Freq_CT.df",
                          "Freq_EOCX_CT.df", "Freq_EOCX_F_CT.df", "Freq_EOCX_M_CT.df",
                          "Freq_PreCX_CT.df", "Freq_PreCX_F_CT.df", "Freq_PreCX_M_CT.df")
  rm(Freq_sample.df,
     Freq_Sex.df,Freq_Cach.df, Freq_CT.df,
     Freq_EOCX_CT.df, Freq_EOCX_F_CT.df, Freq_EOCX_M_CT.df,
     Freq_PreCX_CT.df, Freq_PreCX_F_CT.df, Freq_PreCX_M_CT.df)

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

#### New Barplot(2022/10/27) ####
  library(ggpubr)
  Freq_All2.df <- data.frame(Freq_All.df,str_split(Freq_All.df$Pheno_Type, pattern = "\\.",
                                                   n = Inf, simplify = TRUE))
  colnames(Freq_All2.df)[c((ncol(Freq_All2.df)-1),ncol(Freq_All2.df))] <- c("CCstate","Sex")

  ## Plot Number
  plt.ManyGroup <- ggboxplot(Freq_All2.df, x = "Cell_Type", y = "Number",
                             color = "CCstate",
                             fill = "CCstate",
                             lwd=0.8,
                             # palette = "jco",
                            add = "jitter", # short.panel.labs = T
  ) #+ ylim(0, LabelY*1.2)
  plt.ManyGroup %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                   XtextSize=18,  YtextSize=,18, xangle = 90,
                                   LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
    scale_color_manual(values = c('#ba0449','#3d3c99'))+
    scale_fill_manual(values =alpha( c('#ba0449','#3d3c99'),0.5)) -> plt.ManyGroup
  plt.ManyGroup

  ## Plot Percent
  plt.ManyGroup2 <- ggboxplot(Freq_All2.df, x = "Cell_Type", y = "Percent",
                             color = "CCstate",
                             fill = "CCstate",
                             lwd=0.8,
                             # palette = "jco",
                             add = "jitter", # short.panel.labs = T
  ) #+ ylim(0, LabelY*1.2)
  plt.ManyGroup2 %>% BeautifyggPlot(.,AspRat=1,LegPos = c(0.86, 0.85),AxisTitleSize=1.7,
                                   XtextSize=18,  YtextSize=,18, xangle = 90,
                                   LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
    scale_color_manual(values = c('#ba0449','#3d3c99'))+
    scale_fill_manual(values =alpha( c('#ba0449','#3d3c99'),0.5)) -> plt.ManyGroup2
  plt.ManyGroup2




##### Export PDF file #####
  pdf(file = paste0(Save.Path,"/",SampleType,"_CellCount_LinePlot.pdf"),
      width = 7, height = 7 )
    plt.ManyGroup2 %>% print()
    plt.ManyGroup %>% print()
    CellNum_P4 %>% print()
    CellNum_P3 %>% print()
    CellNum_P1 %>% print()
    CellNum_P2 %>% print()
    BarPlot1_1 %>% print()
    BarPlot1_2 %>% print()
    BarPlot2_1 %>% print()
    BarPlot2_2 %>% print()
    BarPlot3_1 %>% print()
    BarPlot3_2 %>% print()
  dev.off() # graphics.off()

  rm(CellNum_P1, CellNum_P2, CellNum_P3, CellNum_P4, BarPlot1_1, BarPlot1_2,
     BarPlot2_1, BarPlot2_2, BarPlot3_1, BarPlot3_2)

#### Save RData ####
  save.image(paste0(Save.Path,"/07_Count_Cell_number.RData"))
