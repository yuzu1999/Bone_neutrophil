
#————————————————————————————————————————————————————————————————————
library(dior)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
set.seed(101)

rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/01_age/")


scobj = dior::read_h5('05_neu_anno.h5')
scobj

scobj@assays$RNA@counts[1:10,1:40]

scobj$preAnno <- scobj$preAnno3

Idents(scobj) = 'preAnno'

head(scobj@meta.data)


scobj <- NormalizeData(scobj, normalization.method = "LogNormalize", scale.factor = 1e4)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj,features = rownames(scobj))




FeatureMatrix <- read.csv("Age_Neu_specific_markers.csv",header=TRUE)

head(FeatureMatrix)

colnames(FeatureMatrix) <- c("Neu_Pcna+Mpo+","Neu_Pcna+Mki67+","Neu_Ltf+Camp","Neu_Ltf+Mmp9+","Neu_Sell+S100a6+", "Neu_Sell+Csf3r+","Neu_Cxcr4+Il1b+")

FeatureMatrix2 <- FeatureMatrix[1:100,]




#画热图————————————————————————————————————
avgexp = AverageExpression(scobj, group.by = "preAnno", return.seurat = TRUE)
heatmat <- GetAssayData(object = avgexp, slot = 'data', assay = 'RNA')

cc <- colorRampPalette(c("white","white", "white", "#52A85F")) 


typelist = c("Neu_Pcna+Mpo+","Neu_Pcna+Mki67+","Neu_Ltf+Camp","Neu_Ltf+Mmp9+","Neu_Sell+S100a6+", "Neu_Sell+Csf3r+","Neu_Cxcr4+Il1b+")

col_ann <- data.frame(preAnno=typelist)
rownames(col_ann) <- typelist




cols1 <- c("#76afda","#5066a1","#b0d45d","#7fb961","#ffe788","#ffc556","#f06152")
names(cols1) <- typelist

cols2 <- c("#6a51a3","#D8D155","#2EBEBE")
names(cols2) <- c("Inflammation",'Oxidative stress','Apoptosis')

ann_colors <- list(preAnno = cols1,Gene=cols2)




your_colors <- colorRampPalette(rev(brewer.pal(n=11,name="RdBu")[c(1:3,5:7,9:11)]))(100)

critical_value <- 0.9
data_min <- -2
data_max <- 2

breaks_first_half <- seq(from = data_min, to = critical_value, length.out = 51)
breaks_second_half <- seq(from = critical_value, to = data_max, length.out = 51)[-1] 
custom_breaks <- c(breaks_first_half, breaks_second_half)


genes <- unique(as.character(as.matrix(FeatureMatrix2)))


genes <- c('Il1b','Il1r2','Cxcl2','Trem1','C5ar1','C5ar2','Cd14','Tnfrsf1a','Ccr1','Cxcr4','Nfe2l2','Txnip','Fth1','Ftl1','Srgn','Msrb1','Hmox1','Slc7a11','Osgin1','Slc40a1','Mcl1','Bcl2l11','Foxo3','Cdkn1b','Bcl3')

annotation_row <- data.frame(Gene=c(rep('Inflammation',10),rep("Oxidative stress",10),rep("Apoptosis",5)))
rownames(annotation_row) <- genes


jpeg("Heatmap_sen.jpg",height=600)
pheatmap(heatmat[genes, typelist],
           show_rownames = TRUE, show_colnames = FALSE,
		   fontsize_row = 12,
		   annotation_col=col_ann,
		   annotation_row = annotation_row,
			annotation_colors=ann_colors,
			color=your_colors,
			breaks = custom_breaks,
			border=FALSE,
           #color = cc(100),
		   gaps_row = c(10, 20),
		   angle_col = "45",
           cluster_rows = FALSE, cluster_cols = FALSE, scale = "row")
dev.off()


pdf("Heatmap_sen.pdf",height=8.75)
pheatmap(heatmat[genes, typelist],
           show_rownames = TRUE, show_colnames = FALSE,
		   fontsize_row = 12,
		   annotation_col=col_ann,
		   annotation_row = annotation_row,
			annotation_colors=ann_colors,
			color=your_colors,
			breaks = custom_breaks,
			border=FALSE,
           #color = cc(100),
		   gaps_row = c(10, 20),
		   angle_col = "45",
           cluster_rows = FALSE, cluster_cols = FALSE, scale = "row")
dev.off()




jpeg("Heatmap_top100.jpg",height=600)
pheatmap(heatmat[genes, typelist],
           show_rownames = FALSE, show_colnames = FALSE,
		   fontsize_row = 2,
		   annotation_col=col_ann,
			annotation_colors=ann_colors,
			color=cols,
			border=FALSE,
           #color = cc(100),
		   angle_col = "45",
           cluster_rows = FALSE, cluster_cols = FALSE, scale = "row")
dev.off()

pdf("Heatmap_top100.pdf",height=8.75)
pheatmap(heatmat[genes, typelist],
           show_rownames = FALSE, show_colnames = FALSE,
		   fontsize_row = 2,
		   annotation_col=col_ann,
			annotation_colors=ann_colors,
			color=cols,
			border=FALSE,
           #color = cc(100),
		   angle_col = "45",
           cluster_rows = FALSE, cluster_cols = FALSE, scale = "row")
dev.off()