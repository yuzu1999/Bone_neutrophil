

#——————————————————————OVX————————————————————————————
setwd("/home/luchun/scRNA/Bone_neu/02_ovx2/")
rm(list=ls())
gc()


genes <- read.csv("OVX_Neu_preAnno4_deg.csv",header=TRUE,row.names=1)
head(genes)
dim(genes)
#[1] 63208     6


genes$pct_diff <- genes$pct_nz_group - genes$pct_nz_reference

write.csv(genes,"OVX_Neu_preAnno4_deg.csv",quote=FALSE,row.names = FALSE)



secreted <- read.table("/home/luchun/scRNA/Bone_neu/secreted_mice_1891.txt",header=FALSE)
head(secreted)


genes2 <- genes[which(genes$names %in% secreted$V1),]
head(genes2)
dim(genes2)
#[1] 2520    9


write.csv(genes2,"OVX_Neu_preAnno4_secreted.csv",quote=FALSE,row.names = FALSE)




#——————————————————————Age————————————————————————————
setwd("/home/luchun/scRNA/Bone_neu/01_age/")
rm(list=ls())
gc()


genes <- read.csv("Age_Neu_preAnno4_deg.csv",header=TRUE,row.names=1)
head(genes)
dim(genes)
#[1] 70528     8


genes$pct_diff <- genes$pct_nz_group - genes$pct_nz_reference

write.csv(genes,"Age_Neu_preAnno4_deg.csv",quote=FALSE,row.names = FALSE)



secreted <- read.table("/home/luchun/scRNA/Bone_neu/secreted_mice_1891.txt",header=FALSE)
head(secreted)


genes2 <- genes[which(genes$names %in% secreted$V1),]
head(genes2)
dim(genes2)
#[1] 3012    9


write.csv(genes2,"Age_Neu_preAnno4_deg_secreted.csv",quote=FALSE,row.names = FALSE)





#————————————————————————————————人————————————————————————————————
setwd("/home/luchun/scRNA/Bone_neu/04_cell/neu/")
rm(list=ls())
gc()


genes <- read.csv("Human_Neu_preAnno4_deg.csv",header=TRUE,row.names=1)
head(genes)
dim(genes)
#[1] 84036     8

genes$pct_diff <- genes$pct_nz_group - genes$pct_nz_reference

write.csv(genes,"Human_Neu_preAnno4_deg.csv",quote=FALSE,row.names = FALSE)


secreted <- read.table("/home/luchun/scRNA/Bone_neu/secreted_human_1891.txt",header=FALSE)
head(secreted)


genes2 <- genes[which(genes$names %in% secreted$V1),]
head(genes2)
dim(genes2)
#[1] 4692    9

write.csv(genes2,"Human_Neu_preAnno4_secreted.csv",quote=FALSE,row.names = FALSE)








