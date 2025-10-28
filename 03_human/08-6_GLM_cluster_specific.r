


library(tidyverse)
library(limma)
library(data.table)

rm(list = ls())
gc()


setwd("/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/")

pdata <- read.csv("Human_cellinfo.csv",header=TRUE,row.names=1)
head(pdata)

exprSet <- fread("Human_activity.csv",header=TRUE)
exprSet[1:4,1:4]

exprSet2 <- as.matrix(exprSet[,2:ncol(exprSet)])
rownames(exprSet2) <- exprSet$V1
exprSet2[1:4,1:4]





#===========整理分组数据为分组矩阵
pdata$preAnno_tmp <- "Other_neu"
pdata$preAnno_tmp[which(pdata$preAnno3=="Neu_7")] <- "Neu_7"
pdata$preAnno_tmp[which(pdata$preAnno3=="Non_neu")] <- "Non_neu"



table(pdata$preAnno_tmp)
group <- factor(pdata$preAnno_tmp,levels = unique(pdata$preAnno_tmp))



design <- model.matrix(~0 + group)

rownames(design) <- colnames(exprSet2)
colnames(design) <- levels(group)

# 差异比较矩阵
cont.matrix <- makeContrasts("Neu_7 vs Other_neu" = Neu_7-Other_neu, 
                             levels = colnames(design)
                             )
cont.matrix




#========== 定义阈值
logFCcutoff <- log2(0)
adjPvalueCutoff <- 0.3


#=========== 进行差异分析
fit <- lmFit(exprSet2, design)

#=========== 针对给定的对比计算估计系数和标准误差
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

jpeg("01-1_SA.jpg")
plotSA(fit2)
dev.off()



summary(decideTests(fit2, p.value = adjPvalueCutoff, lfc = logFCcutoff))


colnames(fit2)


tmp_vs_Other_neu <- topTreat(fit2, coef=1, n=Inf)
tmp_vs_Other_neu$t_abs <- abs(tmp_vs_Other_neu$t)


write.csv(tmp_vs_Other_neu,"./GLM/Neu_7_vs_Other_neu.csv",quote=FALSE)





Neu_6_vs_Other_neu <- read.csv("./GLM/Neu_6_vs_Other_neu.csv",header=TRUE,row.names=1)

#Neu_6 <- rownames(Neu_6_vs_Other_neu)[which(Neu_6_vs_Other_neu$t_abs > 40)]
Neu_6 <- rownames(arrange(Neu_6_vs_Other_neu,-t_abs))[1:20]

write.table(data.frame(TFs=Neu_6),"./GLM/Neu_Cluster_6_TFs.txt",quote=FALSE,row.names=FALSE)






Neu_1_vs_Other_neu <- read.csv("./GLM/Neu_1_vs_Other_neu.csv",header=TRUE,row.names=1)
Neu_2_vs_Other_neu <- read.csv("./GLM/Neu_2_vs_Other_neu.csv",header=TRUE,row.names=1)
Neu_3_vs_Other_neu <- read.csv("./GLM/Neu_3_vs_Other_neu.csv",header=TRUE,row.names=1)
Neu_4_vs_Other_neu <- read.csv("./GLM/Neu_4_vs_Other_neu.csv",header=TRUE,row.names=1)
Neu_5_vs_Other_neu <- read.csv("./GLM/Neu_5_vs_Other_neu.csv",header=TRUE,row.names=1)
Neu_6_vs_Other_neu <- read.csv("./GLM/Neu_6_vs_Other_neu.csv",header=TRUE,row.names=1)
Neu_7_vs_Other_neu <- read.csv("./GLM/Neu_7_vs_Other_neu.csv",header=TRUE,row.names=1)


'''
Neu_1 <- rownames(Neu_1_vs_Other_neu)[which(Neu_1_vs_Other_neu$t_abs > 40)]
Neu_2 <- rownames(Neu_2_vs_Other_neu)[which(Neu_2_vs_Other_neu$t_abs > 40)]
Neu_3 <- rownames(Neu_3_vs_Other_neu)[which(Neu_3_vs_Other_neu$t_abs > 40)]
Neu_4 <- rownames(Neu_4_vs_Other_neu)[which(Neu_4_vs_Other_neu$t_abs > 40)]
Neu_5 <- rownames(Neu_5_vs_Other_neu)[which(Neu_5_vs_Other_neu$t_abs > 40)]
Neu_6 <- rownames(Neu_6_vs_Other_neu)[which(Neu_6_vs_Other_neu$t_abs > 40)]
Neu_7 <- rownames(Neu_7_vs_Other_neu)[which(Neu_7_vs_Other_neu$t_abs > 40)]
'''


library(dplyr)

Neu_1 <- rownames(arrange(Neu_1_vs_Other_neu,-t_abs))[1:20]
Neu_2 <- rownames(arrange(Neu_2_vs_Other_neu,-t_abs))[1:20]
Neu_3 <- rownames(arrange(Neu_3_vs_Other_neu,-t_abs))[1:20]
Neu_4 <- rownames(arrange(Neu_4_vs_Other_neu,-t_abs))[1:20]
Neu_5 <- rownames(arrange(Neu_5_vs_Other_neu,-t_abs))[1:20]
Neu_6 <- rownames(arrange(Neu_6_vs_Other_neu,-t_abs))[1:20]
Neu_7 <- rownames(arrange(Neu_7_vs_Other_neu,-t_abs))[1:20]


TFs <- unique(c(Neu_1,Neu_2,Neu_3,Neu_4,Neu_5,Neu_6,Neu_7))

write.table(data.frame(TFs=TFs),"./GLM/Neu_Cluster_TFs.txt",quote=FALSE,row.names=FALSE)


















