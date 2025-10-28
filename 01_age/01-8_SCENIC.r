


library(SCopeLoomR) 
#library(devtools)
#devtools::install_local("/home/luchun/Install/SCopeLoomR-master.zip")


library(AUCell)
library(SCENIC) 
#devtools::install_local("/home/luchun/Install/SCENIC-master.zip")


library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(dior)
library(SCENIC)


rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/")



sce_SCENIC <- open_loom("./data/Age_total_SCENIC.loom")

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)
ls(regulons)




len=c()
for(i in ls(regulons)){
	tmp <- length(regulons[[i]])
	len <- c(len,tmp)
}

len_names <- data.frame(tfs=names(regulons),len=paste(" (",len," genes)",sep=""))

len_names$new_name <- paste(substr(len_names$tfs, 1, nchar(len_names$tfs) - 3),len_names$len,sep="")

write.csv(len_names,"Age_len_names.csv",quote=FALSE,row.names=FALSE)





regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')


pydata <- dior::read_h5('Age_total_adata.h5')
table(pydata$preAnno)


pydata$preAnno2 <- "Neu"
pydata$preAnno2[which(pydata$preAnno == "non-neutrophil")] <- "Non_neu"

table(pydata$preAnno2)




pydata$preAnno3 <- "Non_neu"
pydata$preAnno3[which(pydata$preAnno == 'Elane+Prtn3+ progenitor')] <- "Neu_1"
pydata$preAnno3[which(pydata$preAnno == 'H2afz+Hmgb2+ proliferating')] <- "Neu_2"
pydata$preAnno3[which(pydata$preAnno == 'Ngp+Lcn2+ immature')] <- "Neu_3"
pydata$preAnno3[which(pydata$preAnno == 'Ifitm6+Ltf+ immature')] <- "Neu_4"
pydata$preAnno3[which(pydata$preAnno == 'Retnlg+Mmp8+ mature')] <- "Neu_5"
pydata$preAnno3[which(pydata$preAnno == 'Ccl6+Sell+ mature')] <- "Neu_6"
pydata$preAnno3[which(pydata$preAnno == 'Il1b+Srgn+ mature')] <- "Neu_7"
table(pydata$preAnno3)



cellinfo <- pydata@meta.data[,c('preAnno','preAnno2','preAnno3')] 
colnames(cellinfo) <- c('preAnno','preAnno2','preAnno3')

cellinfo$cell <- rownames(cellinfo)
dim(cellinfo)
#[1] 55337     4


head(cellinfo)

write.csv(cellinfo,"Age_cellinfo.csv",quote=FALSE)




#——————————————————计算细胞特意性的转录因子————————————————————————

#通过点图展示————————————————————————

cellTypes <-  as.data.frame(subset(cellinfo,select = 'preAnno'))
selectedResolution <- "preAnno"
sub_regulonAUC <- regulonAUC


identical(rownames(cellinfo),colnames(getAUC(sub_regulonAUC)))

AUC=getAUC(sub_regulonAUC)
dim(AUC)
#[1]   369 55337

write.csv(AUC,"Age_activity.csv",quote=FALSE)





#计算活性均值————————————————————————————
cellTypes <- data.frame(row.names = colnames(pydata), 
                        celltype = pydata$preAnno)


selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), cellTypes[,selectedResolution])
					   
					   
# 去除extened regulons
sub_regulonAUC <- regulonAUC

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
#[1]   369 55337


regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
		
dim(regulonActivity_byGroup)
#[1] 369   8

write.csv(regulonActivity_byGroup,"Age_activity_average.csv",quote=FALSE)




