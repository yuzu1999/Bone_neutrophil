
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

os.chdir('/home/luchun/scRNA/Bone_neu/04_cell/mes/')

sc.settings.set_figure_params(dpi=50, facecolor="white")


adata=sc.read_h5ad("/home/luchun/scRNA/Bone_neu/data/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.h5ad")
adata
#AnnData object with n_obs × n_vars = 82742 × 29225

adata.raw = None


adata.obs['cluster_anno_l2'].unique()
adata.obs['cluster_anno_l2'].value_counts()


adata.var_names_make_unique()
adata.obs_names_make_unique()


adata.obs['orig.ident'].unique()
adata.obs['orig.ident'].value_counts()


adata.X=adata.X.astype(np.int64)



#————————————————————mes细胞保留————————————————————————————
adata_all=adata
adata_all.obs['cluster_anno_l2'].value_counts()


clutser_to_stay = ['Fibro-MSC','THY1+ MSC','Adipo-MSC','APOD+ MSC','Osteo-MSC','Osteoblast']
 
# 创建一个布尔掩码
mask = adata_all.obs['cluster_anno_l2'].isin(clutser_to_stay)
 
# 使用mask提取子集
adata_subset = adata_all[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 19868 × 29225

adata_subset.obs['cluster_anno_l2'].value_counts()


adata = adata_subset
sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 19868 × 27453

adata.obs['orig.ident'].unique()
adata.obs['orig.ident'].value_counts()





#——————————————————————Quality Control—————————————————————————————
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")


sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)


sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt","percent.mt"],
    jitter=0.4,size=0.3,
    groupby="orig.ident",
    multi_panel=True,
    show=False,
    save="_01-1.png"
)


sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",show=False,save="_01-2.png")





#—————————————————————————— 去除双细胞  —————————————————————————————
sc.pp.scrublet(adata, batch_key="orig.ident")
adata
#AnnData object with n_obs × n_vars = 19868 × 27453


adata.obs['predicted_doublet'].value_counts()
#False    19783
#True        85
#Name: predicted_doublet, dtype: int64


adata.write_h5ad('01_Raw.h5ad',compression='gzip')
#adata = sc.read_h5ad("01_Raw.h5ad")

 

adata = adata[adata.obs['predicted_doublet']== False, :]
adata
#View of AnnData object with n_obs × n_vars = 19783 × 27453



sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 19783 × 27446






#——————————————————————Normalization—————————————————————————————
#存储原始数据以便后续还原
ov.utils.store_layers(adata,layers='counts')
adata.layers['counts']=adata.X.copy()


sc.pp.normalize_total(adata)
sc.pp.log1p(adata)


adata.obs['orig.ident'].unique()
adata.obs['orig.ident'].value_counts()



#———————————————————Feature selection———————————————————————————
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="orig.ident")
sc.pl.highly_variable_genes(adata,show=False,save="_01-3.png")

'''
# 获取只有特异性基因的数据集
adata = adata[:, adata.var.highly_variable]
# 回归每个细胞的总计数和表达的线粒体基因的百分比的影响。
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# 将每个基因缩放到单位方差。阈值超过标准偏差 10。
sc.pp.scale(adata, max_value=10)
'''
'''
adata.var['highly_variable'].value_counts()

# 定义需要匹配的前缀
prefixes = ['MT-', 'RPS', 'RPL', 'HSP']

# 更新 'highly_variable' 列的值
for prefix in prefixes:
    adata.var.loc[adata.var.index.str.startswith(prefix), 'highly_variable'] = False


adata.var['highly_variable'].value_counts()

'''
 
 
adata = adata[:, adata.var.highly_variable]
adata
#View of AnnData object with n_obs × n_vars = 19783 × 2000


#sc.pp.regress_out(adata, ['percent.mt'])
sc.pp.scale(adata, max_value=10)


 
#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_01-4.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2", "VIM","CXCL12","LEPR","THY1","APOD","PPARG","LPL","RUNX2","SP7","IBSP"],size=8,ncols=4,cmap="Reds",show=False,save="_01-5.png")

sc.pl.umap(adata,color=['HBA1','HBA2',"pct_counts_mt"],cmap="Reds",show=False,save="_01-6.png")




#——————————————————————Clustering——————————————————————————————————
for res in [4]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
    
sc.pl.umap(adata, color=['leiden_res_4.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_01-7.png")


sc.pl.umap(adata, color=[ 'leiden_res_4.00'],show=False,save="_01-8.png")



#——————————————————Manual cell-type annotation——————————————————
marker_genes = {
    "Hematopoietic": ["PTPRC"],
    "HSPC": ['KIT','CD34','SPINK2', 'ADGRG1', 'CDK6'],
    "Lymp_prog": ["FLT3"],
    "Ery_prog": ["GATA1",'TFRC'],
    "Mega_prog": ["ITGA2B"],
    "Mye_prog": ['CSF1R','ELANE', 'MPO',"PRTN3"],
    "Eos/baso_prog": ["MS4A2"],
    "Neutrophil": ['LTF','CAMP','LCN2','CXCR4','CD101',
'MMP8','MMP9','CXCR2','S100A8', 'S100A9','CSF3R'],
    "Mast":["KIT","TPSAB1","HPGD","TPSB2",'MS4A2','GATA2'],
    #"Baso": ['CD200R3', 'MCPT8', 'PRSS34'],
    "Mono": ["CD14","LYZ","CSF1R","VCAN",'FN1', 'CCR2', 'F13A1'],
    "Macro": [ 'MRC1',"ITGAM","CD86",'C1QA', 'VCAM1',"ADGRE1"],
    "DC":['ITGAX', 'CD74','IRF8','BST2'],
    "T": ['CD3G','CD3D','CCL5',"CD4","FOXP3","CD8A"],
    "NK": ['KLRD1','NCAM1'],
    "Pre/pro_B":['VPREB1'],
    "B": ['CD79A','CD19','IGHM', 'EBF1'],
    "Plasma": ['JCHAIN', 'IGLC2', 'MZB1'],
    "Erythro":['HBA1','HBA2' ],
    
    "Mesenchymal":["CXCL12","PDGFRA"],
    "Osteolineage":['RUNX2',"SP7","IBSP","BGLAP","COL1A1"],
    "Adipo-lineage":["APOE","CEBPA","PPARG","LPL"],
    "Fibro":["DCN","PDPN","CSPG4"], 
    "Smooth_muscle":['ACTA2','RGS5','TAGLN'],
    "Endo":['CDH5','PECAM1','VWF'],
    "Chondrocyte":["SOX9"]   
    
}



sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_4.00', standard_scale="var",dendrogram=True,show=False,save="_01-9.png")


import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(18, 4))


sc.pl.violin(
    adata,
    ["n_genes_by_counts"],
    jitter=0.4,size=0.3,
    groupby="leiden_res_4.00",
    multi_panel=True,
    show=False,
    ax=ax,
    save="_01-10.png"
)






#—————————————————————— 想要剔除的细胞类型列表——————————————————
adata.obs['leiden_res_4.00'].value_counts()

adata
#AnnData object with n_obs × n_vars = 19783 × 2000

clutser_to_remove = ['27','22','50','37','49','51']
 
# 创建一个布尔掩码
mask = ~adata.obs['leiden_res_4.00'].isin(clutser_to_remove)
 
# 使用mask提取子集
adata_subset = adata[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 19081 × 2000


sc.pl.umap(adata_subset, color=[ 'leiden_res_4.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_test.png")
sc.pl.umap(adata_subset, color=['HBA1','HBA2',"pct_counts_mt"],cmap="Reds",show=False,save="_test2.png")



adata.write_h5ad('02_Raw_umap.h5ad',compression='gzip')

#adata = sc.read_h5ad("02_Raw_umap.h5ad")


adata_subset
#View of AnnData object with n_obs × n_vars = 19081 × 2000


adata_counts=adata_subset.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 19081 × 27446

print(adata_counts.X)

ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 19081 × 27446

print(adata_counts.X)




# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts','total_counts','pct_counts_mt']],
    var=adata_counts.var[['features']]
)



new_adata.write_h5ad('03_remove_contamin.h5ad',compression='gzip')








