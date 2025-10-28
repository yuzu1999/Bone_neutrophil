
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

os.chdir('/home/luchun/scRNA/Bone_neu/01_age/')

sc.settings.set_figure_params(dpi=50, facecolor="white")



adata_all = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/03_age_ovx2/02_sce_anno.h5ad")
adata_all
#AnnData object with n_obs × n_vars = 102453 × 2000




#—————————————————————— 想要的细胞类型列表——————————————————
adata_all.obs['batch'].value_counts()


clutser_to_stay = ['4W', '12W', '24M']
 
# 创建一个布尔掩码
mask = adata_all.obs['batch'].isin(clutser_to_stay)
 
# 使用mask提取子集
adata_subset = adata_all[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 58698 × 2000



adata_counts=adata_subset.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 58698 × 21623
print(adata_counts.X)



ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 58698 × 21623
print(adata_counts.X)



# 使用counts的layer重新构建adata
adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['batch', 'dataset', 'n_genes_by_counts','total_counts','pct_counts_mt','celltype']],
    var=adata_counts.var[['gene_ids']]
)
adata
#AnnData object with n_obs × n_vars = 58698 × 21623



sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 58698 × 21216


adata.obs['batch'].value_counts()
adata.obs['dataset'].value_counts()

adata.obs['batch'] = pd.Categorical(adata.obs['batch'],categories=['4W', '12W', '24M'])






#——————————————————————Normalization—————————————————————————————
#存储原始数据以便后续还原
ov.utils.store_layers(adata,layers='counts')
adata.layers['counts']=adata.X.copy()


sc.pp.normalize_total(adata)
sc.pp.log1p(adata)




#———————————————————Feature selection———————————————————————————
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch")
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
prefixes = ['mt-', 'Rps', 'Rpl', 'Hsp']

# 更新 'highly_variable' 列的值
for prefix in prefixes:
    adata.var.loc[adata.var.index.str.startswith(prefix), 'highly_variable'] = False


adata.var['highly_variable'].value_counts()
#False    16187
#True      2957
 
'''



adata = adata[:, adata.var.highly_variable]
adata
#View of AnnData object with n_obs × n_vars = 58698 × 2000


sc.pp.scale(adata, max_value=10)



 
 
#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_02-2.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.pl.umap(adata,color=['batch',"celltype"],cmap="Reds",size=3,show=False,save="_02-3.png")

  
adata.obs['celltype'].unique()
adata.obs['celltype'].value_counts()                     
                     
adata.uns['celltype_colors'] = ["#abddff","#b20000","#167153","#b0d45d","#cca69c","#4a1486","#9e9ac8","#df65b0","#2873B3","#2EBEBE","#F0E442","#d73027","#fc9272","#e8743c","#ffc556", "#bf812d"]
                          
sc.pl.umap(adata, color='celltype',size=3,legend_loc='on data',legend_fontsize=5, legend_fontoutline=2,show=False, save="_02-4.png")
sc.pl.umap(adata, color='celltype',size=3,legend_loc='on data',legend_fontsize=5, legend_fontoutline=2,show=False, save="_02-4.pdf")

sc.pl.umap(adata, color='celltype',size=3,show=False, save="_02-5.png")
sc.pl.umap(adata, color='celltype',size=3,show=False, save="_02-5.pdf")



import pandas as pd
import numpy as np

adata.obs['celltype2'] = 'Other' 
adata.obs.loc[adata.obs['celltype'].isin(['Neutrophil', 'Neu_prog']), 'celltype2'] = 'Neutrophil'

adata.obs['celltype2'] = pd.Categorical(adata.obs['celltype2'],categories=['Neutrophil','Other'])   

adata.uns['celltype2_colors'] = ["#2873B3","lightgrey"]


sc.pl.umap(adata, color='celltype2',size=3,show=False, save="_02-5.png")
sc.pl.umap(adata, color='celltype2',size=3,show=False, save="_02-5.pdf")




adata.obs['batch'].unique()
adata.obs['batch'].value_counts()


sc.tl.embedding_density(adata, groupby="batch")
sc.pl.embedding_density(adata, groupby="batch",ncols=3,show=False, save="_02-5.png")



marker_genes = {
    "HSPC": ['Kit','Cd34','Adgrg1', 'Cdk6'],
    "Mega_prog": ["Itga2b"],
    "Ery_prog": ["Gata1",'Tfrc'],
    "Erythro":['Hbb-bt','Hbb-bs' , 'Hba-a1'],
    "Baso": ["Fcer1a",'Cd200r3', 'Mcpt8', 'Prss34'],
     
    "Mye_prog": ['Csf1r','Elane', 'Mpo',"Prtn3"],
    "Mono": ["Cd14","Vcan",'Fn1', 'Ccr2', 'F13a1'],
    "Macro": ['Mrc1',"Cd86",'C1qa', 'Vcam1',"Adgre1"],
    "Neutrophil": ['Csf3r','S100a8', 'S100a9'],

    "pDC":['Siglech','Bst2','Irf8','Irf7','Tcf4'],
    "T": ['Cd3g','Cd3d'],
    "NK": ['Klrd1','Klrk1','Klrb1c','Ncr1'],
    "Pre/pro_B":['Vpreb1'],
    "B": ['Cd79a','Cd19','Ighm'],
    "Plasma": ['Jchain', 'Iglc2', 'Mzb1']     
}

sc.pl.dotplot(adata, marker_genes, groupby='celltype', standard_scale="var",show=False,save="_02-6.png")

sc.pl.dotplot(adata, marker_genes, groupby='celltype', standard_scale="var",show=False,save="_02-6.pdf")





adata.write_h5ad('04_sce_anno.h5ad',compression='gzip')


#adata = sc.read_h5ad("04_sce_anno.h5ad")
#adata




#——————————————————————堆积条形图————————————————————————————

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


grp = "celltype"
color= 'celltype_colors'


# 提取数据
df = adata.obs[['batch', grp]].copy()

# 计算每个batch中每个celltype的占比
counts = df.groupby(['batch', grp]).size().unstack(fill_value=0)
proportions = counts.div(counts.sum(axis=1), axis=0)


colors = adata.uns[color]

# 绘制堆积条形图
ax = proportions.plot(kind='barh', stacked=True, figsize=(9,5), color=colors,width=0.8)
# 去掉背景网格线
ax.grid(False)
plt.xlabel('Proportion')
plt.ylabel('Sample')
plt.title('Proportions within each sample')
font_properties = FontProperties(size=12) 
plt.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left',prop=font_properties)
#plt.xticks(rotation=45)
plt.tight_layout()  # 自动调整布局，以适应标签

plt.savefig('./figures/proportion_stacked_02-6.pdf')




