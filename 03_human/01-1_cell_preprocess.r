
#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253355/suppl/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz



library(SeuratDisk)
library(Seurat)
rm(list=ls())
gc()

setwd("/home/luchun/scRNA/Bone_neu/data/")


sce_raw <- readRDS("GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds")
sce_raw
#An object of class Seurat
#29452 features across 82742 samples within 1 assay
#Active assay: RNA (29452 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, UMAP_dim30


sce_raw$cluster_anno_l2 <- as.character(sce_raw$cluster_anno_l2)
sce_raw$cluster_anno_coarse <- as.character(sce_raw$cluster_anno_coarse)
sce_raw$cluster_anno_l1 <- as.character(sce_raw$cluster_anno_l1)



sce_raw@assays$RNA@counts[1:4,1:4]

sce <- CreateSeuratObject(sce_raw@assays$RNA@counts,assay = "RNA",min.cells = 3, meta.data = sce_raw@meta.data)
mode(sce$cluster_anno_coarse)



SaveH5Seurat(sce, filename = "GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.h5Seurat")

Convert("GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.h5Seurat", dest = "h5ad")
