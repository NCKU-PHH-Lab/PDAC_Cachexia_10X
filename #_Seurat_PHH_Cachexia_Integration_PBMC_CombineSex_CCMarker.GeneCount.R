##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

PathName = setwd(getwd())
RVersion = "2022-09-07_PBMC_Main"
dir.create(paste0(PathName,"/",RVersion))

PBMC.CCMar <- read.delim2(paste0(PathName,"/",RVersion,"/PBMC_CCMarker_SPA.tsv"))
##### Load libray #####
library(dplyr)


##### 3 Cell Type #####
PBMC.CCMar.Dut <- PBMC.CCMar[grep("Duc", PBMC.CCMar$Type), ]
PBMC.CCMar.Fib <- PBMC.CCMar[grep("Fib", PBMC.CCMar$Type), ]
PBMC.CCMar.Mac <- PBMC.CCMar[grep("Mac", PBMC.CCMar$Type), ]

PBMC.CCMar.Dut.Pos <- PBMC.CCMar.Dut[grep("Pos", PBMC.CCMar.Dut$Type), ]
PBMC.CCMar.Dut.Neg <- PBMC.CCMar.Dut[grep("Neg", PBMC.CCMar.Dut$Type), ]

PBMC.CCMar.Fib.Pos <- PBMC.CCMar.Fib[grep("Pos", PBMC.CCMar.Fib$Type), ]
PBMC.CCMar.Fib.Neg <- PBMC.CCMar.Fib[grep("Neg", PBMC.CCMar.Fib$Type), ]

PBMC.CCMar.Mac.Pos <- PBMC.CCMar.Mac[grep("Pos", PBMC.CCMar.Mac$Type), ]
PBMC.CCMar.Mac.Neg <- PBMC.CCMar.Mac[grep("Neg", PBMC.CCMar.Mac$Type), ]


#### Seperate genes #####
## Test
new.PBMC.CCMar.Dut.Pos <- strsplit(as.character(PBMC.CCMar.Dut.Pos$Marker), ", ")

## All
PBMC.CCMar.CT <- PBMC.CCMar[,1]
new.PBMC.CCMar <- strsplit(as.character(PBMC.CCMar$Marker), ", ")
names(new.PBMC.CCMar) <- PBMC.CCMar.CT

new.PBMC.CCMar.df <- do.call(rbind.data.frame, new.PBMC.CCMar)
rownames(new.PBMC.CCMar.df) <- PBMC.CCMar.CT
colnames(new.PBMC.CCMar.df) <- seq(1:ncol(new.PBMC.CCMar.df))

# write.csv( new.PBMC.CCMar.df ,
#            file = paste0(PathName,"/",RVersion,"/PBMC_CCMarkers_CellType.csv"),
#            #sep = ",",
#            quote = F,
#            row.names = T
# )

new.PBMC.CCMar.df.t <- t(new.PBMC.CCMar.df)

##### All Union Method1 #####
# https://stackoverflow.com/questions/13597091/stacking-columns-in-data-frame-into-one-column-in-r
library(reshape2)
new.PBMC.CCMar.df2 <- melt(t(new.PBMC.CCMar.df),  variable.name =  PBMC.CCMar.CT)
new.PBMC.CCMar.df3 <- data.frame(Genes=unique(new.PBMC.CCMar.df2$value))

write.csv( new.PBMC.CCMar.df3 ,
           file = paste0(PathName,"/",RVersion,"/PBMC_CCMarkers_Unique.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

source("FUN_HSsymbol2MMsymbol.R")
source("FUN_MMsymbol2HSsymbol.R")
library(org.Hs.eg.db)
library(DBI)

new.PBMC.CCMar.df3.Hm <- Alias2HSsymbol(new.PBMC.CCMar.df3,"Genes")
new.PBMC.CCMar.df3.Hm2 <- new.PBMC.CCMar.df3.Hm[!new.PBMC.CCMar.df3.Hm$HS.symbol == 0,]
write.csv( new.PBMC.CCMar.df3.Hm2 ,
           file = paste0(PathName,"/",RVersion,"/PBMC_CCMarkers_Unique.Hm.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

##### All Union Method2 #####
new.PBMC.CCMar.df2_2 <- do.call(rbind, lapply(new.PBMC.CCMar, as.data.frame))
colnames(new.PBMC.CCMar.df2_2)[1]="Genes"
new.PBMC.CCMar.df3_2 <- data.frame(Genes=unique(new.PBMC.CCMar.df2_2$Genes))

source("FUN_HSsymbol2MMsymbol.R")
source("FUN_MMsymbol2HSsymbol.R")
library(org.Hs.eg.db)
library(DBI)
new.PBMC.CCMar.df3_2.Hm <- Alias2HSsymbol(new.PBMC.CCMar.df3,"Genes")
new.PBMC.CCMar.df3_2.Hm2 <- new.PBMC.CCMar.df3_2.Hm[!new.PBMC.CCMar.df3_2.Hm$HS.symbol == 0,]
write.csv( new.PBMC.CCMar.df3_2.Hm2 ,
           file = paste0(PathName,"/",RVersion,"/PBMC_CCMarkers_Unique.Hm_Test.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)


##### Count gene number #####
## Test3 OK
df <- do.call(rbind, lapply(new.PBMC.CCMar, as.data.frame))
colnames(df )[1]="Genes"
df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
CountGene<- function(new.PBMC.CCMar){
  df <- do.call(rbind, lapply(new.PBMC.CCMar, as.data.frame))
  colnames(df )[1]="Genes"
  df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))
  return(df_Count)
  }

new.PBMC.CCMar_Total_CountGene <- CountGene(new.PBMC.CCMar)

new.PBMC.CCMar_Pos <- new.PBMC.CCMar[grep("EO", names(new.PBMC.CCMar))]
new.PBMC.CCMar_Pos_CountGene <- CountGene(new.PBMC.CCMar_Pos)

new.PBMC.CCMar_Neg <- new.PBMC.CCMar[grep("LO", names(new.PBMC.CCMar))]
new.PBMC.CCMar_Neg_CountGene <- CountGene(new.PBMC.CCMar_Neg)

new.PBMC.CCMar_Total_CountGene_Summary <- full_join(new.PBMC.CCMar_Total_CountGene, new.PBMC.CCMar_Pos_CountGene, by="Var1")
new.PBMC.CCMar_Total_CountGene_Summary <- full_join(new.PBMC.CCMar_Total_CountGene_Summary, new.PBMC.CCMar_Neg_CountGene, by="Var1")
colnames(new.PBMC.CCMar_Total_CountGene_Summary) <- c("Genes","Total","EO","LO")
new.PBMC.CCMar_Total_CountGene_Summary[is.na(new.PBMC.CCMar_Total_CountGene_Summary)] <- 0
write.csv( new.PBMC.CCMar_Total_CountGene_Summary ,
           file = paste0(PathName,"/",RVersion,"/PBMC_CCMarkers_Count.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)



#NG# new.PBMC.CCMar_max <- matrix(unlist(new.PBMC.CCMar),ncol = length(new.PBMC.CCMar))

df <- new.PBMC.CCMar[[1]]
for (i in 2:length(new.PBMC.CCMar)) {
  df <- merge(df, new.PBMC.CCMar[[i]], by = "row.names", all=T)
  rownames(df) <- df$Row.names
  df <- df[ , !(names(df) %in% "Row.names")]
}

colnames(df) <- names(new.PBMC.CCMar)
df %>% arrange(!!!syms(colnames(.))) -> df2

df3 <- arrange_all(df)

df0 <- df
df0[is.na(df0)] <- 0
df4 <- apply(df0, 2, sort,decreasing = T)
df4[df4==0] <- ""


write.csv( df4 ,
           file = paste0(PathName,"/",RVersion,"/PBMC_CCMarkers_CellType.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
##### Combine Sex Result #####

PBMC.CCMar.Sex <- read.delim2(paste0(PathName,"/",RVersion,"/PBMC_CCMarker_SSA.tsv"))

## All
PBMC.CCMar.Sex.CT <- PBMC.CCMar.Sex[,1]
new.PBMC.CCMar.Sex <- strsplit(as.character(PBMC.CCMar.Sex$Intersect), ", ")
names(new.PBMC.CCMar.Sex) <- PBMC.CCMar.Sex.CT


# df.Sex <- new.PBMC.CCMar.Sex[[1]]
# for (i in 2:length(new.PBMC.CCMar.Sex)) {
#   df.Sex <- merge(df.Sex, new.PBMC.CCMar.Sex[[i]], by = "row.names", all=T)
#   rownames(df.Sex) <- df.Sex$Row.names
#   df.Sex <- df.Sex[ , !(names(df.Sex) %in% "Row.names")]
# }
#
# colnames(df.Sex) <- names(new.PBMC.CCMar.Sex)
# df.Sex %>% arrange(!!!syms(colnames(.))) -> df.Sex2
#
# df.Sex3 <- arrange_all(df.Sex)
#
# df.Sex0 <- df.Sex
# df.Sex0[is.na(df.Sex0)] <- 0
# df.Sex4 <- apply(df.Sex0, 2, sort,decreasing = T)
# df.Sex4[df.Sex4==0] <- ""
#

new.PBMC.CCMar.Sex_Total_CountGene <- CountGene(new.PBMC.CCMar.Sex)

new.PBMC.CCMar.Sex_Pos <- new.PBMC.CCMar.Sex[grep("Pos", names(new.PBMC.CCMar.Sex))]
new.PBMC.CCMar.Sex_Pos_CountGene <- CountGene(new.PBMC.CCMar.Sex_Pos)

new.PBMC.CCMar.Sex_Neg <- new.PBMC.CCMar.Sex[grep("Neg", names(new.PBMC.CCMar.Sex))]
new.PBMC.CCMar.Sex_Neg_CountGene <- CountGene(new.PBMC.CCMar.Sex_Neg)

new.PBMC.CCMar.Sex_Total_CountGene_Summary <- full_join(new.PBMC.CCMar.Sex_Total_CountGene, new.PBMC.CCMar.Sex_Pos_CountGene, by="Var1")
new.PBMC.CCMar.Sex_Total_CountGene_Summary <- full_join(new.PBMC.CCMar.Sex_Total_CountGene_Summary, new.PBMC.CCMar.Sex_Neg_CountGene, by="Var1")
colnames(new.PBMC.CCMar.Sex_Total_CountGene_Summary) <- c("Genes","Total","EO","LO")
new.PBMC.CCMar.Sex_Total_CountGene_Summary[is.na(new.PBMC.CCMar.Sex_Total_CountGene_Summary)] <- 0


new.PBMC.CCMar.Sex_Total_CountGene_SUM <- full_join(new.PBMC.CCMar_Total_CountGene_Summary,new.PBMC.CCMar.Sex_Total_CountGene_Summary,by="Genes")
colnames(new.PBMC.CCMar.Sex_Total_CountGene_SUM) <- c("Genes","Total.Com","EO.Com","LO.Com","Total.SexInt","EO.SexInt","LO.SexInt")
new.PBMC.CCMar.Sex_Total_CountGene_SUM[is.na(new.PBMC.CCMar.Sex_Total_CountGene_SUM)] <- 0
write.csv( new.PBMC.CCMar.Sex_Total_CountGene_SUM ,
           file = paste0(PathName,"/",RVersion,"/PBMC_CCMarkers_SUM.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)
