## Ref: https://www.jianshu.com/p/9e5b7ffcf80f

FeaturePlot(PBMC.combined, features = c("Marco"), min.cutoff = "q9")

scRNA.SeuObj <- PBMC.combined
GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

GeneExp_Marco.df <- GeneExp.df[row.names(GeneExp.df) == "Marco",]
gene_num <- ncol(GeneExp_Marco.df)
sum(GeneExp_Marco.df)

plot(density(GeneExp_Top2a.df %>% t() ), col=rainbow(gene_num)[1], lty=1,
     xlab = "Expression level", # ylim = c(0,1.5), main = ""
)


GeneExp_Top2a.df <- GeneExp.df[row.names(GeneExp.df) == "Top2a",]
gene_num <- ncol(GeneExp_Top2a.df)
sum(GeneExp_Top2a.df)

GeneExp_Sub.df <- GeneExp.df[row.names(GeneExp.df) %in% c("Top2a","Marco"),] %>% t() %>% as.data.frame()
gene_num <- ncol(GeneExp_Sub.df)
sum(GeneExp_Sub.df)


plot(density(GeneExp_Sub.df[,1]), col=rainbow(gene_num)[1], lty=1,
     xlab = "Expression level", main = names(GeneExp_Sub.df)[1])
polygon(density(GeneExp_Sub.df[,1]),col=rainbow(gene_num)[1])


# 绘制所有基因的密度分布图
plot(density(GeneExp_Sub.df[,1]), col=rainbow(gene_num)[1], lty=1,
     xlab = "Expression level", #ylim = c(0,1.5), main = ""
     )
# polygon(density(GeneExp_Sub.df[,1]),col=rainbow(gene_num)[1])
# 添加其他基因的密度曲线
for (i in seq(2,gene_num)){
  lines(density(GeneExp_Sub.df[,i]), col=rainbow(gene_num)[i], lty=i)
  #polygon(density(GeneExp_Sub.df[,i]),col=rainbow(gene_num)[i])
}
# 添加图例
legend("topright", inset = 0.02, title = "Gene", names(GeneExp_Sub.df),
       col = rainbow(gene_num), lty = seq(1,gene_num), bg = "gray")

