
##### Path setting #####
PathName = setwd(getwd())
RVersion = "20211229"
dir.create(paste0(PathName,"/",RVersion))


##### Load Data #####
SC.Candidate <- read.delim2(paste0(PathName,"/SC_CCMarker.CombineSex.tsv"))

##### Split different cell Type #####
SC.Candidate.Dut <- SC.Candidate[grep("Duc", SC.Candidate$Type), ]
SC.Candidate.Fib <- SC.Candidate[grep("Fib", SC.Candidate$Type), ]
SC.Candidate.Mac <- SC.Candidate[grep("Mac", SC.Candidate$Type), ]

SC.Cand.Dut.Pos <- SC.Candidate.Dut[grep("Pos", SC.Candidate.Dut$Type), ]
SC.Cand.Dut.Neg <- SC.Candidate.Dut[grep("Neg", SC.Candidate.Dut$Type), ]

SC.Cand.Fib.Pos <- SC.Candidate.Fib[grep("Pos", SC.Candidate.Fib$Type), ]
SC.Cand.Fib.Neg <- SC.Candidate.Fib[grep("Neg", SC.Candidate.Fib$Type), ]

SC.Cand.Mac.Pos <- SC.Candidate.Mac[grep("Pos", SC.Candidate.Mac$Type), ]
SC.Cand.Mac.Neg <- SC.Candidate.Mac[grep("Neg", SC.Candidate.Mac$Type), ]



# ###
# library(tidyr)
# 
# new.SC.Candidate.Dut.Pos <- separate(SC.Candidate.Dut.Pos, Marker ,  ", ")


###
new.SC.Candidate.Dut.Pos <- strsplit(as.character(SC.Candidate.Dut.Pos$Marker), ", ")

### All
SC.Candidate.CT <- SC.Candidate[,1]
new.SC.Candidate <- strsplit(as.character(SC.Candidate$Marker), ", ")
names(new.SC.Candidate) <- SC.Candidate.CT 

new.SC.Candidate.df <- do.call(rbind.data.frame, new.SC.Candidate)
rownames(new.SC.Candidate.df) <- SC.Candidate.CT
colnames(new.SC.Candidate.df) <- seq(1:ncol(new.SC.Candidate.df))

write.csv( new.SC.Candidate.df ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_CellType.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

new.SC.Candidate.df.t <- t(new.SC.Candidate.df)

# https://stackoverflow.com/questions/13597091/stacking-columns-in-data-frame-into-one-column-in-r
library(reshape2)
new.SC.Candidate.df2 <- melt(t(new.SC.Candidate.df),  variable.name =  SC.Candidate.CT)
new.SC.Candidate.df3 <- data.frame(Genes=unique(new.SC.Candidate.df2$value))

write.csv( new.SC.Candidate.df3 ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Unique.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

source("FUN_HSsymbol2MMsymbol.R")
source("FUN_MMsymbol2HSsymbol.R")
library(org.Hs.eg.db)
library(DBI)

new.SC.Candidate.df3.Hm <- Alias2HSsymbol(new.SC.Candidate.df3,"Genes")
new.SC.Candidate.df3.Hm2 <- new.SC.Candidate.df3.Hm[!new.SC.Candidate.df3.Hm$HS.symbol == 0,]
write.csv( new.SC.Candidate.df3.Hm2 ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Unique.Hm.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)



##### Combine SubType to main 6 Type #####
library(reshape2)
CombineSubType <- function(SC.Cand.Dut.Pos, ColName = "SC.Cand.Dut.Pos"){
  # Dut
  SC.Cand.Dut.Pos.CT <- SC.Cand.Dut.Pos[,1]
  new.SC.Cand.Dut.Pos <- strsplit(as.character(SC.Cand.Dut.Pos$Marker), ", ")
  names(new.SC.Cand.Dut.Pos) <- SC.Cand.Dut.Pos.CT 
  
  new.SC.Cand.Dut.Pos.df <- do.call(rbind.data.frame, new.SC.Cand.Dut.Pos)
  rownames(new.SC.Cand.Dut.Pos.df) <- SC.Cand.Dut.Pos.CT
  colnames(new.SC.Cand.Dut.Pos.df) <- seq(1:ncol(new.SC.Cand.Dut.Pos.df))
  new.SC.Cand.Dut.Pos.df.t <- t(new.SC.Cand.Dut.Pos.df)
  
  # https://stackoverflow.com/questions/13597091/stacking-columns-in-data-frame-into-one-column-in-r
  # library(reshape2)
  new.SC.Cand.Dut.Pos.df2 <- melt(t(new.SC.Cand.Dut.Pos.df),  variable.name =  SC.Cand.Dut.Pos.CT)
  new.SC.Cand.Dut.Pos.df3 <- data.frame(Genes=unique(new.SC.Cand.Dut.Pos.df2$value))
  colnames(new.SC.Cand.Dut.Pos.df3) <- ColName
  return(new.SC.Cand.Dut.Pos.df3)
}

DupGene <- new.SC.Candidate.df2[duplicated(new.SC.Candidate.df2$value),]
DupGene2 <- DupGene[!duplicated(DupGene$value),]
new.SC.Cand.Dut.Pos <- CombineSubType(SC.Cand.Dut.Pos,"SC.Cand.Dut.Pos")
new.SC.Cand.Dut.Pos2 <- data.frame(SC.Cand.Dut.Pos=new.SC.Cand.Dut.Pos[!new.SC.Cand.Dut.Pos$SC.Cand.Dut.Pos %in%
                                               DupGene2$value,])
new.SC.Cand.Dut.Neg <- CombineSubType(SC.Cand.Dut.Neg,"SC.Cand.Dut.Neg")
new.SC.Cand.Dut.Neg2 <- data.frame(SC.Cand.Dut.Neg=new.SC.Cand.Dut.Neg[!new.SC.Cand.Dut.Neg$SC.Cand.Dut.Neg %in%
                                                                         DupGene2$value,])

new.SC.Cand.Fib.Pos <- CombineSubType(SC.Cand.Fib.Pos,"SC.Cand.Fib.Pos")
new.SC.Cand.Fib.Pos2 <- data.frame(SC.Cand.Fib.Pos=new.SC.Cand.Fib.Pos[!new.SC.Cand.Fib.Pos$SC.Cand.Fib.Pos %in%
                                                                         DupGene2$value,])
new.SC.Cand.Fib.Neg <- CombineSubType(SC.Cand.Fib.Neg,"SC.Cand.Fib.Neg") 
new.SC.Cand.Fib.Neg2 <- data.frame(SC.Cand.Fib.Neg=new.SC.Cand.Fib.Neg[!new.SC.Cand.Fib.Neg$SC.Cand.Fib.Neg %in%
                                                                         DupGene2$value,])


new.SC.Cand.Mac.Pos <- CombineSubType(SC.Cand.Mac.Pos,"SC.Cand.Mac.Pos")
new.SC.Cand.Mac.Pos2 <- data.frame(SC.Cand.Mac.Pos=new.SC.Cand.Mac.Pos[!new.SC.Cand.Mac.Pos$SC.Cand.Mac.Pos %in%
                                                                         DupGene2$value,])
new.SC.Cand.Mac.Neg <- CombineSubType(SC.Cand.Mac.Neg,"SC.Cand.Mac.Neg")
new.SC.Cand.Mac.Neg2 <- data.frame(SC.Cand.Mac.Neg=new.SC.Cand.Mac.Neg[!new.SC.Cand.Mac.Neg$SC.Cand.Mac.Neg %in%
                                                                         DupGene2$value,])

# https://statisticsglobe.com/create-data-frame-of-unequal-lengths-in-r


CCMarker_6Type.df <- data.frame(new.SC.Cand.Dut.Pos2,new.SC.Cand.Dut.Neg2)
max_length <- max(c(nrow(new.SC.Cand.Dut.Pos2), nrow(new.SC.Cand.Dut.Neg2),
                    nrow(new.SC.Cand.Fib.Pos2), nrow(new.SC.Cand.Fib.Neg2),
                    nrow(new.SC.Cand.Mac.Pos2), nrow(new.SC.Cand.Mac.Neg2)))    # Find out maximum length
max_length                                      # Print maximum length

CCMarker_6Type.df <- data.frame(Dut.Pos = c(new.SC.Cand.Dut.Pos2$SC.Cand.Dut.Pos,                 # Create data frame with unequal vectors
                                rep(NA, max_length - nrow(new.SC.Cand.Dut.Pos2))),
                                Dut.Neg = c(new.SC.Cand.Dut.Neg2$SC.Cand.Dut.Neg,
                                rep(NA, max_length - nrow(new.SC.Cand.Dut.Neg2))),
                                Fib.Pos = c(new.SC.Cand.Fib.Pos2$SC.Cand.Fib.Pos,
                                rep(NA, max_length - nrow(new.SC.Cand.Fib.Pos2))),
                                Fib.Neg = c(new.SC.Cand.Fib.Neg2$SC.Cand.Fib.Neg,
                                rep(NA, max_length - nrow(new.SC.Cand.Fib.Neg2))),
                                Mac.Pos = c(new.SC.Cand.Mac.Pos2$SC.Cand.Mac.Pos,
                                rep(NA, max_length - nrow(new.SC.Cand.Mac.Pos2))),
                                Mac.Neg = c(new.SC.Cand.Mac.Neg2$SC.Cand.Mac.Neg,
                                rep(NA, max_length - nrow(new.SC.Cand.Mac.Neg2)))
                                
                                )
CCMarker_6Type.df                                            # Print final data frame
CCMarker_6Type.df2 <- CCMarker_6Type.df
CCMarker_6Type.df2[is.na(CCMarker_6Type.df2)] <- ""
write.csv( CCMarker_6Type.df2 ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Cell6Type.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

##### Count gene number #####
## Test1
library(tibble)
df<-enframe(new.SC.Candidate)

## Test2
library(dplyr)
df<-bind_rows(new.SC.Candidate)
df<- tibble::as_tibble_row(new.SC.Candidate)

## Test3 OK
df <- do.call(rbind, lapply(new.SC.Candidate, as.data.frame)) 
colnames(df )[1]="Genes"
df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))
 
CountGene<- function(new.SC.Candidate){
  df <- do.call(rbind, lapply(new.SC.Candidate, as.data.frame)) 
  colnames(df )[1]="Genes"
  df_Count <- as.data.frame(table(df$Genes)) %>% arrange(.,desc(Freq))
  return(df_Count)
  }

new.SC.Candidate_Total_CountGene <- CountGene(new.SC.Candidate)

new.SC.Candidate_Pos <- new.SC.Candidate[grep("Pos", names(new.SC.Candidate))]
new.SC.Candidate_Pos_CountGene <- CountGene(new.SC.Candidate_Pos)

new.SC.Candidate_Neg <- new.SC.Candidate[grep("Neg", names(new.SC.Candidate))]
new.SC.Candidate_Neg_CountGene <- CountGene(new.SC.Candidate_Neg)

new.SC.Candidate_Total_CountGene_Summary <- full_join(new.SC.Candidate_Total_CountGene, new.SC.Candidate_Pos_CountGene, by="Var1")
new.SC.Candidate_Total_CountGene_Summary <- full_join(new.SC.Candidate_Total_CountGene_Summary, new.SC.Candidate_Neg_CountGene, by="Var1")  
colnames(new.SC.Candidate_Total_CountGene_Summary) <- c("Genes","Total","Malignant","Protective")
new.SC.Candidate_Total_CountGene_Summary[is.na(new.SC.Candidate_Total_CountGene_Summary)] <- 0
write.csv( new.SC.Candidate_Total_CountGene_Summary ,
           file = paste0(PathName,"/",RVersion,"/SC_CCMarkers_Count.csv"),
           #sep = ",",
           quote = F,
           row.names = F
)

## Test4
# new.SC.Candidate.df.t <- as.data.frame(new.SC.Candidate.df.t)
new.SC.Candidate.df.t2 <- matrix(nrow = nrow(new.SC.Candidate.df.t),ncol=ncol(new.SC.Candidate.df.t))
colnames(new.SC.Candidate.df.t2) <- colnames(new.SC.Candidate.df.t)
new.SC.Candidate.df.t3  <- ""
for(i in colnames(new.SC.Candidate.df.t)){
  TTT <- as.data.frame(unique(new.SC.Candidate.df.t[,i])) 
  colnames(TTT) <- colnames(new.SC.Candidate.df.t)[i]
#new.SC.Candidate.df.t2[1:nrow(TTT),i] <- TTT 
new.SC.Candidate.df.t3 <- cbind(new.SC.Candidate.df.t3,TTT )
}

names(new.SC.Candidate.df.t)[1]

