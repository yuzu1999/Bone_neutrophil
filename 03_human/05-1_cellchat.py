
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
os.chdir('/home/luchun/scRNA/Bone_neu/04_cell/cellchat2/')



#——————————————————————读入—————————————————————————————
neu = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/04_cell/neu/07_Anno_filter.h5ad")
neu
#AnnData object with n_obs × n_vars = 7048 × 2000


neu.obs['preAnno'].unique()
neu.obs['preAnno'].value_counts()


preAnno4 = dict()

for i in ['AZU1+KCNQ5+ progenitor', 'PDE4D+MCU+ proliferating', 'LTF+CAMP+ immature','S100A12+MMP9+ immature', 'S100A6+NCF1+ mature', 'NAMPT+IFITM2+ mature','NEAT1+FTH1+ mature']:
    preAnno4[str(i)] = "Undefined"


print(preAnno4)


####开始定义细胞类型
for i in ['AZU1+KCNQ5+ progenitor','PDE4D+MCU+ proliferating']:
    preAnno4[i] = "Neu_AZU1+"


for i in ['LTF+CAMP+ immature','S100A12+MMP9+ immature']:
    preAnno4[i] = "Neu_LTF+"


for i in ['S100A6+NCF1+ mature','NAMPT+IFITM2+ mature']:
    preAnno4[i] = "Neu_SELL+"
    

for i in ['NEAT1+FTH1+ mature']:
    preAnno4[i] = "Neu_IL1B+"



print(preAnno4)


neu.obs['preAnno4'] = neu.obs['preAnno'].map(preAnno4)

neu.obs['preAnno4'].unique()

neu.obs['preAnno4'] = pd.Categorical(neu.obs['preAnno4'],categories=["Neu_AZU1+","Neu_LTF+","Neu_SELL+","Neu_IL1B+"])  

neu.obs['preAnno'] = neu.obs['preAnno4']
neu.obs['preAnno'].unique()

neu2=neu.raw.to_adata().copy()
neu2
#AnnData object with n_obs × n_vars = 7048 × 21155


print(neu2.X)

ov.utils.retrieve_layers(neu2,layers='counts')
neu2
#AnnData object with n_obs × n_vars = 7048 × 21155

print(neu2.X)

# 使用counts的layer重新构建adata
new_neu = sc.AnnData(
    X=neu2.X,
    obs=neu2.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts','total_counts','pct_counts_mt','preAnno']],
    var=neu2.var[['features']]
)

new_neu






mes = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/04_cell/mes/05_Anno.h5ad")
mes
#AnnData object with n_obs × n_vars = 19081 × 2000

mes.obs['preAnno'].unique()
mes.obs['preAnno'].value_counts()


mes2=mes.raw.to_adata().copy()
mes2
#AnnData object with n_obs × n_vars = 19081 × 27287

print(mes2.X)

ov.utils.retrieve_layers(mes2,layers='counts')
mes2
#AnnData object with n_obs × n_vars = 19081 × 27287

print(mes2.X)


# 使用counts的layer重新构建adata
new_mes = sc.AnnData(
    X=mes2.X,
    obs=mes2.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts','total_counts','pct_counts_mt','preAnno']],
    var=mes2.var[['features']]
)

new_mes





adata=sc.concat([new_neu,new_mes],join='outer', merge='first')
adata
#AnnData object with n_obs × n_vars = 26129 × 27437

adata.var_names_make_unique()
adata.obs_names_make_unique()



adata.obs['preAnno'].unique()
adata.obs['preAnno'].value_counts()

adata.obs['preAnno'] = pd.Categorical(adata.obs['preAnno'],categories=["Neu_AZU1+","Neu_LTF+","Neu_SELL+","Neu_IL1B+",'Fibro-MSC','THY1+ MSC','Adipo-MSC','APOD+ MSC','Osteo-MSC','Osteoblast'])  
  
adata.obs['preAnno'].unique()







#——————————————————————Normalization—————————————————————————————
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)


adata.write_h5ad('./adata_cellchat.h5ad',compression='gzip')

#adata = sc.read_h5ad("/public/home/guest/scRNA/Bone_neu/adata_cellchat.h5ad")
#adata

#os.getcwd()  ##查看当前路径
#os.chdir('/public/home/guest/scRNA/Bone_neu/')




#————————————————————数据导出———————————————————————————————————————
mat=pd.DataFrame(data=adata.X.todense(),index=adata.obs_names,columns=adata.var_names)
mat.to_hdf("mat.h5","mat")


meta=pd.DataFrame(data=adata.obs)
meta.to_csv('metadata.tsv',sep="\t") 



















