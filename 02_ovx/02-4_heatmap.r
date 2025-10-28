

import omicverse as ov
#print(f"omicverse version: {ov.__version__}")
#ov.utils.ov_plot_set()
import os
import numpy as np
import scanpy as sc
#print(f"scanpy version: {sc.__version__}")
import anndata as ad
import pandas as pd


os.getcwd()  ##查看当前路径

os.chdir('/home/luchun/scRNA/Bone_neu/02_ovx/')

sc.settings.set_figure_params(dpi=50, facecolor="white")



adata=sc.read_h5ad("05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 19324 × 2000



adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 19324 × 15981


print(adata_counts.X)

ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 19324 × 15981

print(adata_counts.X)




# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['batch', 'dataset', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'celltype', 'preAnno', 'S_score', 'G2M_score', 'phase']],
    var=adata_counts.var[['gene_ids']]
)



adata=new_adata


import diopy
diopy.output.write_h5(adata,'05_neu_anno.h5')







#————————————————————————————————————————————————————————————————————
library(dior)
library(Seurat)
library(scProgram)
library(pheatmap)
library(RColorBrewer)
set.seed(101)

rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/02_ovx/")


scobj = dior::read_h5('05_neu_anno.h5')
scobj

scobj@assays$RNA@counts[1:10,1:40]


Idents(scobj) = 'preAnno'

head(scobj@meta.data)


scobj <- NormalizeData(scobj, normalization.method = "LogNormalize", scale.factor = 1e4)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj,features = rownames(scobj))




FeatureMatrix <- read.csv("OVX_Neu_specific_markers.csv",header=TRUE)

head(FeatureMatrix)

colnames(FeatureMatrix) <- c('Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature')

FeatureMatrix2 <- FeatureMatrix[1:100,]




#画热图————————————————————————————————————
avgexp = AverageExpression(scobj, group.by = "preAnno", return.seurat = T)
heatmat <- GetAssayData(object = avgexp, slot = 'data', assay = 'RNA')

cc <- colorRampPalette(c("white","white", "white", "#52A85F")) 


typelist = c('Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature')

col_ann <- data.frame(preAnno=typelist)
rownames(col_ann) <- typelist



ann_colors <- list(
  preAnno = c(`Elane+Prtn3+ progenitor` = "#76afda",`H2afz+Hmgb2+ proliferating` = "#8264CC",`Ngp+Lcn2+ immature` = "#009E73",`Ifitm6+Ltf+ immature` = "#F0E442",`Retnlg+Mmp8+ mature`="#f06152",`Ccl6+Sell+ mature`="#0072B2",`Il1b+Srgn+ mature`="#fc9272"))

cols <- colorRampPalette(rev(brewer.pal(n=11,name="RdBu")[c(1:3,5:7,9:11)]))(100)



genes <- unique(as.character(as.matrix(FeatureMatrix2)))

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

