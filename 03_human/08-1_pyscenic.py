
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

os.chdir('/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/')

sc.settings.set_figure_params(dpi=50, facecolor="white")



#——————————————————————其余细胞随机抽样——————————————————————————————
adata_all = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/04_cell/04_sce_anno.h5ad")
adata_all
#AnnData object with n_obs × n_vars = 80474 × 2000


adata_neu=sc.read_h5ad("/home/luchun/scRNA/Bone_neu/04_cell/neu/07_Anno_filter.h5ad")
adata_neu
#AnnData object with n_obs × n_vars = 7048 × 2000


neu_cell_ids = adata_neu.obs.index
adata_filtered = adata_all[~adata_all.obs.index.isin(neu_cell_ids)].copy()
adata_filtered
#AnnData object with n_obs × n_vars = 73430 × 2000



adata_filtered.obs['celltype'].unique()
adata_filtered.obs['celltype'].value_counts()


cell_types_to_include = [
    "HSPC","Mega_prog","Ery_prog","Erythroid","Baso","Mono_prog","Monocyte","Macrophage","pDC","T","NK","Pre/pro_B","B","Plasma","VSMC",'Fibro−MSC',"MSC","Endothelial"]


# 过滤出这些细胞类型
filtered_adata = adata_filtered[adata_filtered.obs['celltype'].isin(cell_types_to_include)]
filtered_adata 
#View of AnnData object with n_obs × n_vars = 69023 × 2000




adata1=filtered_adata.raw.to_adata().copy()
adata1
#AnnData object with n_obs × n_vars = 69023 × 29188
print(adata1.X)


ov.utils.retrieve_layers(adata1,layers='counts')
adata1
#AnnData object with n_obs × n_vars = 69023 × 29188
print(adata1.X)



adata1.obs['preAnno']="non-neutrophil"




# 使用counts的layer重新构建adata
adata1_2 = sc.AnnData(
    X=adata1.X,
    obs=adata1.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt','preAnno']],
    var=adata1.var[['features']]
)




adata2=adata_neu.raw.to_adata().copy()
adata2
#AnnData object with n_obs × n_vars = 7048 × 21155
print(adata2.X)

ov.utils.retrieve_layers(adata2,layers='counts')
adata2
#AnnData object with n_obs × n_vars = 7048 × 21155
print(adata2.X)

adata2.obs['preAnno'].unique()

# 使用counts的layer重新构建adata
adata2_2 = sc.AnnData(
    X=adata2.X,
    obs=adata2.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt','preAnno']],
    var=adata2.var[['features']]
)



#————————————————————————————loom——————————————————————————
adata=sc.concat([adata1_2,adata2_2],join='outer', merge='first')
adata
#AnnData object with n_obs × n_vars = 76071 × 29188

adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.obs['preAnno'].unique()


adata.write_h5ad('Human_total_adata.h5ad',compression='gzip')
#adata = sc.read_h5ad("Human_total_adata.h5ad")


import diopy
diopy.output.write_h5(adata,'Human_total_adata.h5')




import loompy as lp


'''
ranking_feather = pd.read_feather("/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
ranking_feather2 = pd.read_feather("/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

ranking_feather.columns
ranking_feather2.columns

gene_saved = list(set(ranking_feather.columns) & set(ranking_feather2.columns))
gene_saved2 = list(set(gene_saved) & set(adata.var_names))


len(gene_saved2)
#19162


# 使用布尔掩码过滤adata对象
adata2 = adata[:, gene_saved2]
adata2
#View of AnnData object with n_obs × n_vars = 76071 × 19162


adata=adata2
'''

print(adata.X)

row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names)
}


lp.create('/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data/Human_total_adata.loom', adata.X.transpose(), row_attrs, col_attrs)








#——————————————————————其余细胞随机抽样——————————————————————————————
adata_all = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/04_cell/04_sce_anno.h5ad")
adata_all
#AnnData object with n_obs × n_vars = 80474 × 2000


adata_neu=sc.read_h5ad("/home/luchun/scRNA/Bone_neu/04_cell/neu/07_Anno_filter.h5ad")
adata_neu
#AnnData object with n_obs × n_vars = 7048 × 2000


neu_cell_ids = adata_neu.obs.index
adata_filtered = adata_all[~adata_all.obs.index.isin(neu_cell_ids)].copy()
adata_filtered
#AnnData object with n_obs × n_vars = 73430 × 2000



adata_filtered.obs['celltype'].unique()
adata_filtered.obs['celltype'].value_counts()


cell_types_to_include = [
    "HSPC","Mega_prog","Ery_prog","Erythroid","Baso","Mono_prog","Monocyte","Macrophage","pDC","T","NK","Pre/pro_B","B","Plasma","VSMC",'Fibro−MSC',"MSC","Endothelial"]


# 过滤出这些细胞类型
filtered_adata = adata_filtered[adata_filtered.obs['celltype'].isin(cell_types_to_include)]
filtered_adata 
#View of AnnData object with n_obs × n_vars = 69023 × 2000


# 随机抽取1200个细胞
n_cells_to_sample = 1200
sampled_indices = filtered_adata.obs.sample(n=n_cells_to_sample, random_state=42).index
sampled_adata = filtered_adata[sampled_indices]

sampled_adata
#View of AnnData object with n_obs × n_vars = 1200 × 2000

sampled_adata.obs['celltype'].unique()
sampled_adata.obs['celltype'].value_counts()




adata_counts=sampled_adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 1200 × 29188
print(adata_counts.X)

ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 1200 × 29188
print(adata_counts.X)


adata_counts.obs['preAnno']="non-neutrophil"


# 使用counts的layer重新构建adata
adata1 = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt','preAnno']],
    var=adata_counts.var[['features']]
)





#——————————————————————————中性粒随机抽样——————————————————————————

adata_neu.obs['preAnno'].unique()
adata_neu.obs['preAnno'].value_counts()


n_cells_to_sample = 300
sampled_cells = []


# 遍历每种类型，随机抽取300个细胞
for cell_type in adata_neu.obs['preAnno'].unique():
    subset = adata_neu[adata_neu.obs['preAnno'] == cell_type]
    
    sampled_subset = subset.obs.sample(n=n_cells_to_sample, random_state=42).index
    
    sampled_cells.append(sampled_subset.tolist())



# 合并所有抽取的细胞索引
all_sampled_indices = pd.Index([cell for sublist in sampled_cells for cell in sublist])


# 创建新的adata对象
sampled_adata = adata_neu[all_sampled_indices]
sampled_adata
#View of AnnData object with n_obs × n_vars = 2100 × 2000

sampled_adata.obs['preAnno'].unique()
sampled_adata.obs['preAnno'].value_counts()




adata_counts=sampled_adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 2100 × 21155
print(adata_counts.X)


ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 2100 × 21155
print(adata_counts.X)



# 使用counts的layer重新构建adata
adata2 = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt','preAnno']],
    var=adata_counts.var[['features']]
)




#————————————————————————————loom——————————————————————————
adata=sc.concat([adata1,adata2],join='outer', merge='first')
adata
#AnnData object with n_obs × n_vars = 3300 × 29188



adata.write_h5ad('Human_sample_adata.h5ad',compression='gzip')

#adata = sc.read_h5ad("Human_sample_adata.h5ad")
#adata

import diopy
diopy.output.write_h5(adata,'Human_sample_adata.h5')



'''
ranking_feather = pd.read_feather("/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
ranking_feather2 = pd.read_feather("/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

ranking_feather.columns
ranking_feather2.columns

gene_saved = list(set(ranking_feather.columns) & set(ranking_feather2.columns))
gene_saved2 = list(set(gene_saved) & set(adata.var_names))

len(gene_saved2)
#19162

# 使用布尔掩码过滤adata对象
adata2 = adata[:, gene_saved2]
adata2
#View of AnnData object with n_obs × n_vars = 3300 × 19162


adata=adata2


'''

import loompy as lp


print(adata.X)

row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) 
}
lp.create('/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data/Human_sample_adata.loom', adata.X.transpose(), row_attrs, col_attrs)


























  
    