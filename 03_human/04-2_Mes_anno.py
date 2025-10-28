
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

adata=sc.read_h5ad("03_remove_contamin.h5ad")
adata
#AnnData object with n_obs × n_vars = 19081 × 27446


adata.obs['cluster_anno_l2'].unique()
adata.obs['cluster_anno_l2'].value_counts()


adata.var_names_make_unique()
adata.obs_names_make_unique()


adata.obs['orig.ident'].unique()
adata.obs['orig.ident'].value_counts()


adata.X=adata.X.astype(np.int64)

sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 19081 × 27287





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
sc.pl.highly_variable_genes(adata,show=False,save="_02-1.png")

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
#View of AnnData object with n_obs × n_vars = 19081 × 2000


#sc.pp.regress_out(adata, ['percent.mt'])
sc.pp.scale(adata, max_value=10)





#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_02-2.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2", "VIM","CXCL12","LEPR","THY1","APOD","PPARG","LPL","RUNX2","SP7","IBSP"],size=8,ncols=4,show=False,save="_02-3.png")

sc.pl.umap(adata,color=['HBA1','HBA2',"pct_counts_mt"],cmap="Reds",show=False,save="_02-4.png")


#sc.pl.umap(adata,color=["MSRB1"],size=5,show=False,save="_test.png")




#——————————————————————Clustering——————————————————————————————————
for res in [0.5,1,1.5,2,2.5,3]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
  
  
sc.pl.umap(adata, color=['leiden_res_0.50', 'leiden_res_1.00','leiden_res_1.50','leiden_res_2.00','leiden_res_2.50','leiden_res_3.00'],ncols=3,legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_02-5.png")

sc.pl.umap(adata, color=['leiden_res_1.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_02-6.png")




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


sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_1.00', standard_scale="var",dendrogram=True,show=False,save="_02-7.png")




import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(18, 4))

sc.pl.violin(
    adata,
    ["percent.mt"],
    jitter=0.4,size=0.3,
    groupby="leiden_res_1.00",
    multi_panel=True,
    ax=ax,
    show=False,
    save="_02-8.png"
)



adata.write_h5ad('04_merge.h5ad',compression='gzip')


adata=sc.read_h5ad("04_merge.h5ad")
adata





#——————————————————————————bbknn——————————————————————————————————
sc.external.pp.bbknn(adata, batch_key="orig.ident") 
sc.tl.umap(adata)


adata.obs['cluster_anno_l2'] = pd.Categorical(adata.obs['cluster_anno_l2'],categories=['Fibro-MSC','THY1+ MSC','Adipo-MSC','APOD+ MSC','Osteo-MSC','Osteoblast'])    


sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2", "VIM","CXCL12","LEPR","THY1","APOD","PPARG","LPL","RUNX2","SP7","IBSP"],size=8,ncols=4,show=False,save="_02-9.png")

sc.pl.umap(adata,color=['HBA1','HBA2',"pct_counts_mt"],cmap="Reds",show=False,save="_02-10.png")




'''
#————————————————————————harmony————————————————————————————————
#https://support.parsebiosciences.com/hc/en-us/articles/7704577188500-How-to-analyze-a-1-million-cell-data-set-using-Scanpy-and-Harmony
import harmonypy
help(harmonypy.run_harmony)


import scanpy.external as sce

sce.pp.harmony_integrate(adata, key='orig.ident',theta=1,lamb=1)
#1,1
#1,2
#0,2
#0,1

#λ值越小，整合力度越大，Default lambda=1.调整范围一般在0.5-2之间
#θ值越大，聚类越多样，Default theta=2.

adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)


sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2", "VIM","CXCL12","LEPR","THY1","APOD","PPARG","LPL","RUNX2","SP7","IBSP"],size=8,ncols=4,show=False,save="_04-1.png")



adata.obs['cluster_anno_l2'] = pd.Categorical(adata.obs['cluster_anno_l2'],categories=["Fibro-MSC","APOD+ MSC","Osteoblast","Osteo-MSC","THY1+ MSC","Adipo-MSC"])     
adata.uns['cluster_anno_l2_color'] = ["#BC9DCC","#A65628","#54B0E4","#222F75","#1B9E77","#B2DF8A"]


sc.pl.umap(adata, color=['cluster_anno_l2'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_04-2.png")

'''






#——————————————————Manual cell-type annotation——————————————————

adata.obs['preAnno'] = adata.obs['cluster_anno_l2']

adata.obs['preAnno'] = pd.Categorical(adata.obs['preAnno'],categories=['Fibro-MSC','THY1+ MSC','Adipo-MSC','APOD+ MSC','Osteo-MSC','Osteoblast'])    
          

adata.uns['preAnno_colors'] = ["#428AC9","#129392","#FFCC4F","#F37E78","#883A96","#cca69c"] 
              
sc.pl.umap(adata, color='preAnno',legend_loc='on data',legend_fontsize=6, legend_fontoutline=2,size=10,show=False, save="_02-11.png")

sc.pl.umap(adata, color='preAnno',size=10,show=False, save="_02-12.png")


      
          
marker_genes = {
    "Mesenchymal":["CXCL12","PDGFRA"],
    "Osteolineage":['RUNX2',"SP7","IBSP","BGLAP","COL1A1"],
    "Adipo-lineage":["APOE","CEBPA","PPARG","LPL"],
    "Fibro":["DCN","PDPN","CSPG4"], 
    "Smooth_muscle":['ACTA2','RGS5','TAGLN'],
    "Endo":['CDH5','PECAM1','VWF'],
    "Chondrocyte":["SOX9"]   
}


sc.pl.dotplot(adata, marker_genes, groupby='preAnno', standard_scale="var",dendrogram=False,show=False,save="_02-13.png")






adata.write_h5ad('05_Anno.h5ad',compression='gzip')

#adata=sc.read_h5ad("05_Anno.h5ad")
#adata










