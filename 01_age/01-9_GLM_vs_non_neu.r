


library(tidyverse)
library(limma)
library(data.table)

rm(list = ls())
gc()


setwd("/home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/")

pdata <- read.csv("Age_cellinfo.csv",header=TRUE,row.names=1)
head(pdata)

exprSet <- fread("Age_activity.csv",header=TRUE)
exprSet[1:4,1:4]

exprSet2 <- as.matrix(exprSet[,2:ncol(exprSet)])
rownames(exprSet2) <- exprSet$V1
exprSet2[1:4,1:4]


table(pdata$preAnno3)
group <- factor(pdata$preAnno3,levels = unique(pdata$preAnno3))



#===========整理分组数据为分组矩阵
design <- model.matrix(~0 + group)

rownames(design) <- colnames(exprSet2)
colnames(design) <- levels(group)

# 差异比较矩阵
cont.matrix <- makeContrasts("Neu_1 vs Non_neu" = Neu_1-Non_neu, 
                             "Neu_2 vs Non_neu" = Neu_2-Non_neu, 
                             "Neu_3 vs Non_neu" = Neu_3-Non_neu,
							 "Neu_4 vs Non_neu" = Neu_4-Non_neu, 
                             "Neu_5 vs Non_neu" = Neu_5-Non_neu, 
                             "Neu_6 vs Non_neu" = Neu_6-Non_neu,
							 "Neu_7 vs Non_neu" = Neu_7-Non_neu,
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

jpeg("01-1_Age_SA.jpg")
plotSA(fit2)
dev.off()



summary(decideTests(fit2, p.value = adjPvalueCutoff, lfc = logFCcutoff))
#       Neu_1 vs Non_neu Neu_2 vs Non_neu Neu_3 vs Non_neu Neu_4 vs Non_neu
#Down                320              238              251              247
#NotSig               25               31               28               29
#Up                  171              247              237              240
#       Neu_5 vs Non_neu Neu_6 vs Non_neu Neu_7 vs Non_neu
#Down                252              263              234
#NotSig               27               18               52
#Up                  237              235              230


colnames(fit2)


Neu_1_VS_Non_neu <- topTreat(fit2, coef=1, n=Inf)
Neu_2_VS_Non_neu <- topTreat(fit2, coef=2, n=Inf)
Neu_3_VS_Non_neu <- topTreat(fit2, coef=3, n=Inf)
Neu_4_VS_Non_neu <- topTreat(fit2, coef=4, n=Inf)
Neu_5_VS_Non_neu <- topTreat(fit2, coef=5, n=Inf)
Neu_6_VS_Non_neu <- topTreat(fit2, coef=6, n=Inf)
Neu_7_VS_Non_neu <- topTreat(fit2, coef=7, n=Inf)


Neu_1_VS_Non_neu$t_abs <- abs(Neu_1_VS_Non_neu$t)
Neu_2_VS_Non_neu$t_abs <- abs(Neu_2_VS_Non_neu$t)
Neu_3_VS_Non_neu$t_abs <- abs(Neu_3_VS_Non_neu$t)
Neu_4_VS_Non_neu$t_abs <- abs(Neu_4_VS_Non_neu$t)
Neu_5_VS_Non_neu$t_abs <- abs(Neu_5_VS_Non_neu$t)
Neu_6_VS_Non_neu$t_abs <- abs(Neu_6_VS_Non_neu$t)
Neu_7_VS_Non_neu$t_abs <- abs(Neu_7_VS_Non_neu$t)


write.csv(Neu_1_VS_Non_neu,"./Age_GLM/Neu_1_VS_Non_neu.csv",quote=FALSE)
write.csv(Neu_2_VS_Non_neu,"./Age_GLM/Neu_2_VS_Non_neu.csv",quote=FALSE)
write.csv(Neu_3_VS_Non_neu,"./Age_GLM/Neu_3_VS_Non_neu.csv",quote=FALSE)
write.csv(Neu_4_VS_Non_neu,"./Age_GLM/Neu_4_VS_Non_neu.csv",quote=FALSE)
write.csv(Neu_5_VS_Non_neu,"./Age_GLM/Neu_5_VS_Non_neu.csv",quote=FALSE)
write.csv(Neu_6_VS_Non_neu,"./Age_GLM/Neu_6_VS_Non_neu.csv",quote=FALSE)
write.csv(Neu_7_VS_Non_neu,"./Age_GLM/Neu_7_VS_Non_neu.csv",quote=FALSE)




Neu_1_VS_Non_neu <- read.csv("./Age_GLM/Neu_1_VS_Non_neu.csv",header=TRUE,row.names=1)
Neu_2_VS_Non_neu <- read.csv("./Age_GLM/Neu_2_VS_Non_neu.csv",header=TRUE,row.names=1)
Neu_3_VS_Non_neu <- read.csv("./Age_GLM/Neu_3_VS_Non_neu.csv",header=TRUE,row.names=1)
Neu_4_VS_Non_neu <- read.csv("./Age_GLM/Neu_4_VS_Non_neu.csv",header=TRUE,row.names=1)
Neu_5_VS_Non_neu <- read.csv("./Age_GLM/Neu_5_VS_Non_neu.csv",header=TRUE,row.names=1)
Neu_6_VS_Non_neu <- read.csv("./Age_GLM/Neu_6_VS_Non_neu.csv",header=TRUE,row.names=1)
Neu_7_VS_Non_neu <- read.csv("./Age_GLM/Neu_7_VS_Non_neu.csv",header=TRUE,row.names=1)


'''
Neu_1 <- rownames(Neu_1_VS_Non_neu)[which(Neu_1_VS_Non_neu$t_abs > 100)]
Neu_2 <- rownames(Neu_2_VS_Non_neu)[which(Neu_2_VS_Non_neu$t_abs > 100)]
Neu_3 <- rownames(Neu_3_VS_Non_neu)[which(Neu_3_VS_Non_neu$t_abs > 100)]
Neu_4 <- rownames(Neu_4_VS_Non_neu)[which(Neu_4_VS_Non_neu$t_abs > 100)]
Neu_5 <- rownames(Neu_5_VS_Non_neu)[which(Neu_5_VS_Non_neu$t_abs > 100)]
Neu_6 <- rownames(Neu_6_VS_Non_neu)[which(Neu_6_VS_Non_neu$t_abs > 100)]
Neu_7 <- rownames(Neu_7_VS_Non_neu)[which(Neu_7_VS_Non_neu$t_abs > 100)]
'''
library(dplyr)

Neu_1 <- rownames(arrange(Neu_1_VS_Non_neu,-t_abs))[1:20]
Neu_2 <- rownames(arrange(Neu_2_VS_Non_neu,-t_abs))[1:20]
Neu_3 <- rownames(arrange(Neu_3_VS_Non_neu,-t_abs))[1:20]
Neu_4 <- rownames(arrange(Neu_4_VS_Non_neu,-t_abs))[1:20]
Neu_5 <- rownames(arrange(Neu_5_VS_Non_neu,-t_abs))[1:20]
Neu_6 <- rownames(arrange(Neu_6_VS_Non_neu,-t_abs))[1:20]
Neu_7 <- rownames(arrange(Neu_7_VS_Non_neu,-t_abs))[1:20]


TFs <- unique(c(Neu_1,Neu_2,Neu_3,Neu_4,Neu_5,Neu_6,Neu_7))


write.table(data.frame(TFs=TFs),"./Age_GLM/Neu_VS_Non_neu_TFs.txt",quote=FALSE,row.names=FALSE)


























