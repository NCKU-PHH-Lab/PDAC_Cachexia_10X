PathName = setwd(getwd())
RVersion = "20211222"
dir.create(paste0(PathName,"/",RVersion))

SC.CCMar <- read.delim2(paste0(PathName,"/SC_CCMarker.CombineSex.tsv"))
##### Load libray #####
library(dplyr)


##### 3 Cell Type #####
SC.CCMar.Dut <- SC.CCMar[grep("Duc", SC.CCMar$Type), ]
SC.CCMar.Fib <- SC.CCMar[grep("Fib", SC.CCMar$Type), ]
SC.CCMar.Mac <- SC.CCMar[grep("Mac", SC.CCMar$Type), ]

SC.CCMar.Dut.Pos <- SC.CCMar.Dut[grep("Pos", SC.CCMar.Dut$Type), ]
SC.CCMar.Dut.Neg <- SC.CCMar.Dut[grep("Neg", SC.CCMar.Dut$Type), ]

SC.CCMar.Fib.Pos <- SC.CCMar.Fib[grep("Pos", SC.CCMar.Fib$Type), ]
SC.CCMar.Fib.Neg <- SC.CCMar.Fib[grep("Neg", SC.CCMar.Fib$Type), ]

SC.CCMar.Mac.Pos <- SC.CCMar.Mac[grep("Pos", SC.CCMar.Mac$Type), ]
SC.CCMar.Mac.Neg <- SC.CCMar.Mac[grep("Neg", SC.CCMar.Mac$Type), ]


#### Seperate genes #####
## Test
new.SC.CCMar.Dut.Pos <- strsplit(as.character(SC.CCMar.Dut.Pos$Marker), ", ")

## All
SC.CCMar.CT <- SC.CCMar[,1]
new.SC.CCMar <- strsplit(as.character(SC.CCMar$Marker), ", ")
names(new.SC.CCMar) <- SC.CCMar.CT 

new.SC.CCMar.df <- do.call(rbind.data.frame, new.SC.CCMar)
rownames(new.SC.CCMar.df) <- SC.CCMar.CT
colnames(new.SC.CCMar.df) <- seq(1:ncol(new.SC.CCMar.df))

# write.csv( new.SC.CCMar.df ,
#            file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_CellType.csv"),
#            #sep = ",",
#            quote = F,
#            row.names = T
# )

new.SC.CCMar.df.t <- t(new.SC.CCMar.df)

##### All Union Method1 #####
# https://stackoverflow.com/questions/13597091/stacking-columns-in-data-frame-into-one-column-in-r
library(reshape2)
new.SC.CCMar.df2 <- melt(t(new.SC.CCMar.df),  variable.name =  SC.CCMar.CT)
new.SC.CCMar.df3 <- data.frame(Genes=unique(new.SC.CCMar.df2$value))

write.csv( new.SC.CCMar.df3 ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Unique.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

source("FUN_HSsymbol2MMsymbol.R")
source("FUN_MMsymbol2HSsymbol.R")
library(org.Hs.eg.db)
library(DBI)

new.SC.CCMar.df3.Hm <- Alias2HSsymbol(new.SC.CCMar.df3,"Genes")
new.SC.CCMar.df3.Hm2 <- new.SC.CCMar.df3.Hm[!new.SC.CCMar.df3.Hm$HS.symbol == 0,]
write.csv( new.SC.CCMar.df3.Hm2 ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Unique.Hm.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

##### All Union Method2 #####
new.SC.CCMar.df2_2 <- do.call(rbind, lapply(new.SC.CCMar, as.data.frame)) 
colnames(new.SC.CCMar.df2_2)[1]="Genes"
new.SC.CCMar.df3_2 <- data.frame(Genes=unique(new.SC.CCMar.df2_2$Genes))

source("FUN_HSsymbol2MMsymbol.R")
source("FUN_MMsymbol2HSsymbol.R")
library(org.Hs.eg.db)
library(DBI)
new.SC.CCMar.df3_2.Hm <- Alias2HSsymbol(new.SC.CCMar.df3,"Genes")
new.SC.CCMar.df3_2.Hm2 <- new.SC.CCMar.df3_2.Hm[!new.SC.CCMar.df3_2.Hm$HS.symbol == 0,]
write.csv( new.SC.CCMar.df3_2.Hm2 ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Unique.Hm_Test.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)


##### Count gene number #####
## Test3 OK
df <- do.call(rbind, lapply(new.SC.CCMar, as.data.frame)) 
colnames(df )[1]="Genes"
df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))
 
##
CountGene<- function(new.SC.CCMar){
  df <- do.call(rbind, lapply(new.SC.CCMar, as.data.frame)) 
  colnames(df )[1]="Genes"
  df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))
  return(df_Count)
  }

new.SC.CCMar_Total_CountGene <- CountGene(new.SC.CCMar)

new.SC.CCMar_Pos <- new.SC.CCMar[grep("Pos", names(new.SC.CCMar))]
new.SC.CCMar_Pos_CountGene <- CountGene(new.SC.CCMar_Pos)

new.SC.CCMar_Neg <- new.SC.CCMar[grep("Neg", names(new.SC.CCMar))]
new.SC.CCMar_Neg_CountGene <- CountGene(new.SC.CCMar_Neg)

new.SC.CCMar_Total_CountGene_Summary <- full_join(new.SC.CCMar_Total_CountGene, new.SC.CCMar_Pos_CountGene, by="Var1")
new.SC.CCMar_Total_CountGene_Summary <- full_join(new.SC.CCMar_Total_CountGene_Summary, new.SC.CCMar_Neg_CountGene, by="Var1")  
colnames(new.SC.CCMar_Total_CountGene_Summary) <- c("Genes","Total","Malignant","Protective")
new.SC.CCMar_Total_CountGene_Summary[is.na(new.SC.CCMar_Total_CountGene_Summary)] <- 0
write.csv( new.SC.CCMar_Total_CountGene_Summary ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Count.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)



#NG# new.SC.CCMar_max <- matrix(unlist(new.SC.CCMar),ncol = length(new.SC.CCMar))

df <- new.SC.CCMar[[1]]
for (i in 2:length(new.SC.CCMar)) {
  df <- merge(df, new.SC.CCMar[[i]], by = "row.names", all=T)
  rownames(df) <- df$Row.names
  df <- df[ , !(names(df) %in% "Row.names")]  
}

colnames(df) <- names(new.SC.CCMar)
df %>% arrange(!!!syms(colnames(.))) -> df2

df3 <- arrange_all(df)

df0 <- df
df0[is.na(df0)] <- 0
df4 <- apply(df0, 2, sort,decreasing = T)
df4[df4==0] <- ""


write.csv( df4 ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_CellType.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

