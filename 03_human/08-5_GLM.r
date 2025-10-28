


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


table(pdata$preAnno3)
group <- factor(pdata$preAnno3,levels = unique(pdata$preAnno3))



#===========整理分组数据为分组矩阵
design <- model.matrix(~0 + group)

rownames(design) <- colnames(exprSet2)
colnames(design) <- levels(group)

# 差异比较矩阵
cont.matrix <- makeContrasts("Neu_2 vs Neu_1" = Neu_2-Neu_1, 
                             "Neu_3 vs Neu_2" = Neu_3-Neu_2,
							 "Neu_4 vs Neu_3" = Neu_4-Neu_3,
							 "Neu_5 vs Neu_4" = Neu_5-Neu_4,
							 "Neu_6 vs Neu_5" = Neu_6-Neu_5,
							 "Neu_7 vs Neu_6" = Neu_7-Neu_6,
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
#       Neu_1 vs Non_neu Neu_2 vs Non_neu Neu_3 vs Non_neu Neu_4 vs Non_neu
#Down                320              238              251              247
#NotSig               25               31               28               29
#Up                  171              247              237              240
#       Neu_5 vs Non_neu Neu_6 vs Non_neu Neu_7 vs Non_neu
#Down                252              263              234
#NotSig               27               18               52
#Up                  237              235              230


colnames(fit2)


Neu_2_VS_Neu_1 <- topTreat(fit2, coef=1, n=Inf)
Neu_3_VS_Neu_2 <- topTreat(fit2, coef=2, n=Inf)
Neu_4_VS_Neu_3 <- topTreat(fit2, coef=3, n=Inf)
Neu_5_VS_Neu_4 <- topTreat(fit2, coef=4, n=Inf)
Neu_6_VS_Neu_5 <- topTreat(fit2, coef=5, n=Inf)
Neu_7_VS_Neu_6 <- topTreat(fit2, coef=6, n=Inf)


Neu_2_VS_Neu_1$t_abs <- abs(Neu_2_VS_Neu_1$t)
Neu_3_VS_Neu_2$t_abs <- abs(Neu_3_VS_Neu_2$t)
Neu_4_VS_Neu_3$t_abs <- abs(Neu_4_VS_Neu_3$t)
Neu_5_VS_Neu_4$t_abs <- abs(Neu_5_VS_Neu_4$t)
Neu_6_VS_Neu_5$t_abs <- abs(Neu_6_VS_Neu_5$t)
Neu_7_VS_Neu_6$t_abs <- abs(Neu_7_VS_Neu_6$t)


write.csv(Neu_2_VS_Neu_1,"./GLM/Neu_2_VS_Neu_1.csv",quote=FALSE)
write.csv(Neu_3_VS_Neu_2,"./GLM/Neu_3_VS_Neu_2.csv",quote=FALSE)
write.csv(Neu_4_VS_Neu_3,"./GLM/Neu_4_VS_Neu_3.csv",quote=FALSE)
write.csv(Neu_5_VS_Neu_4,"./GLM/Neu_5_VS_Neu_4.csv",quote=FALSE)
write.csv(Neu_6_VS_Neu_5,"./GLM/Neu_6_VS_Neu_5.csv",quote=FALSE)
write.csv(Neu_7_VS_Neu_6,"./GLM/Neu_7_VS_Neu_6.csv",quote=FALSE)




Neu_2_VS_Neu_1 <- read.csv("./GLM/Neu_2_VS_Neu_1.csv",header=TRUE,row.names=1)
Neu_3_VS_Neu_2 <- read.csv("./GLM/Neu_3_VS_Neu_2.csv",header=TRUE,row.names=1)
Neu_4_VS_Neu_3 <- read.csv("./GLM/Neu_4_VS_Neu_3.csv",header=TRUE,row.names=1)
Neu_5_VS_Neu_4 <- read.csv("./GLM/Neu_5_VS_Neu_4.csv",header=TRUE,row.names=1)
Neu_6_VS_Neu_5 <- read.csv("./GLM/Neu_6_VS_Neu_5.csv",header=TRUE,row.names=1)
Neu_7_VS_Neu_6 <- read.csv("./GLM/Neu_7_VS_Neu_6.csv",header=TRUE,row.names=1)




Neu_2 <- rownames(Neu_2_VS_Neu_1)[which(Neu_2_VS_Neu_1$t_abs > 20)]
Neu_3 <- rownames(Neu_3_VS_Neu_2)[which(Neu_3_VS_Neu_2$t_abs > 20)]
Neu_4 <- rownames(Neu_4_VS_Neu_3)[which(Neu_4_VS_Neu_3$t_abs > 20)]
Neu_5 <- rownames(Neu_5_VS_Neu_4)[which(Neu_5_VS_Neu_4$t_abs > 20)]
Neu_6 <- rownames(Neu_6_VS_Neu_5)[which(Neu_6_VS_Neu_5$t_abs > 20)]
Neu_7 <- rownames(Neu_7_VS_Neu_6)[which(Neu_7_VS_Neu_6$t_abs > 20)]



TFs <- unique(c(Neu_2,Neu_3,Neu_4,Neu_5,Neu_6,Neu_7))

write.table(data.table(TFs=TFs),"./GLM/Neu_Cluster_TFs.txt",quote=FALSE,row.names=FALSE)


























