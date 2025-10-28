
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

os.chdir('/home/luchun/scRNA/Bone_neu/02_ovx2/')

sc.settings.set_figure_params(dpi=50, facecolor="white")



adata = sc.read_h5ad("04_Neu_remove_contamin.h5ad")
adata
#AnnData object with n_obs × n_vars = 18784 × 15981


sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 18784 × 15938


adata.obs['batch'].value_counts()

adata.obs['batch'] = pd.Categorical(adata.obs['batch'],categories=["Sham","OVX"])




#——————————————————————Normalization—————————————————————————————
#存储原始数据以便后续还原
ov.utils.store_layers(adata,layers='counts')
adata.layers['counts']=adata.X.copy()


sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

#在评分之前，需对count数据进行log转化和scale
sc.pp.scale(adata)
print(adata.X)



#——————————————————————  细胞周期  ————————————————————————————————
# .strip()移除字符串头尾指定字符（默认为空格或换行符）或字符序列
cell_cycle_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/mice_cc_genes.txt')]

s_genes = cell_cycle_genes[:42] #列表的开头开始，一直到第 42个元素（但不包括第 42 个元素）的所有元素。
g2m_genes = cell_cycle_genes[42:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]


sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)


adata.obs['phase'].unique()
adata.obs['phase'].value_counts()
adata.obs['phase'] = pd.Categorical(adata.obs['phase'],categories=['G1','S', 'G2M'])
adata.uns['phase_colors'] = ["#167CB4","#FE9F2A","#EA4737"]



adata_cc_genes = adata[:, cell_cycle_genes]

sc.tl.pca(adata_cc_genes,use_highly_variable=False)
sc.pl.pca_scatter(adata_cc_genes, color='phase',show=False,save="_04-1.png")




#——————————————————————Normalization—————————————————————————————
#使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata.layers['counts'],
    obs=adata.obs,
    var=adata.var
)

adata = new_adata

ov.utils.store_layers(adata,layers='counts')
adata.layers['counts']=adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)




#———————————————————Feature selection———————————————————————————
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="batch")
sc.pl.highly_variable_genes(adata,show=False,save="_04-2.png")


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
#False    16128
#True      2957
'''
 
 

adata = adata[:, adata.var.highly_variable]
adata
#View of AnnData object with n_obs × n_vars = 18784 × 2000


sc.pp.scale(adata, max_value=10)



 
#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_04-3.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=15)
sc.tl.umap(adata)

sc.pl.umap(adata,color=["Mki67","Pcna","Mpo",'Elane','Ltf','Camp','Mmp9','Sell','S100a6','Csf3r','Cxcr4','Il1b'],cmap="Reds",size=8,ncols=5,show=False,save="_04-4.png")



adata.obs['batch'].unique()
adata.obs['batch'].value_counts()


sc.tl.embedding_density(adata, groupby="batch")
sc.pl.embedding_density(adata, groupby="batch",ncols=2,show=False, save="_04-5.png")




#——————————————————————Clustering——————————————————————————————————
for res in [0.5,1,1.5, 2,2.5, 3,3.5,4]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
    

sc.pl.umap(adata, color=['leiden_res_0.50','leiden_res_1.00','leiden_res_1.50','leiden_res_2.00','leiden_res_2.50','leiden_res_3.00','leiden_res_3.50','leiden_res_4.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=3,ncols=3,show=False,save="_04-6.png")

sc.pl.umap(adata, color=['leiden_res_1.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=3,show=False,save="_04-7.png")

sc.pl.umap(adata, color=[ 'leiden_res_4.00'],show=False,save="_04-8.png")

sc.pl.dotplot(adata, ['Cxcr4','Il1b'], groupby='leiden_res_4.00', standard_scale="var",dendrogram=True,show=False,save="_04-9.png")



#——————————————————Manual cell-type annotation——————————————————

marker_genes = {
    "Prog": ['Mpo', 'Elane', 'Prtn3', 'Ctsg', 'Ms4a3'],
    "Proliferating": ['Camp', 'Ptma', 'Stmn1','Mki67', 'Ppia', 'Hmgb2', 'Tuba1b', 'Pclaf', 'Anp32b', 'Top2a', 'Smc4', 'Tagln2', 'Birc5', 'Ncapd2'],
    "Immature": ['Ltf', 'Lcn2', 'Ngp', 'Cd177', 'S100a8'],
    "Mature": ['S100a11', 'Ccl6', 'Retnlg', 'S100a6', 'Sell', 'Slpi', 'Clec4d', 'Lcp1', 'Srgn', 'Fpr1', 'Tyrobp', 'Cd33', 'Anxa2', 'Selplg','Fos', 'C5ar1', 'Dusp1', 'Samhd1', 'Csf3r'],
    "Marker": ["Mki67","Pcna","Mpo",'Elane','Ltf','Camp','Mmp9','Sell','S100a6','Csf3r','Cxcr4','Il1b']
  
}

sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_1.00', standard_scale="var",dendrogram=True,show=False,save="_04-9.png")




import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(12, 4))


sc.pl.violin(
    adata,
    ['pct_counts_mt'],
    jitter=0.4,size=0.3,
    groupby="leiden_res_2.00",
    multi_panel=True,
    show=False,
    ax=ax,
    save="_04-11.png"
)

adata.write_h5ad('05_Neu_merge.h5ad',compression='gzip')




#——————————————————Manual cell-type annotation——————————————————

preAnno3 = dict()

for i in range(0,13):
    preAnno3[str(i)] = "Neutrophil"

print(preAnno3)


####开始定义细胞类型
for i in ['12']:
    preAnno3[i] = "Neu_Pcna+Mpo+"


for i in ['7','8']:
    preAnno3[i] = "Neu_Pcna+Mki67+"


for i in ['0']:
    preAnno3[i] = "Neu_Ltf+Camp"


for i in ['1']:
    preAnno3[i] = "Neu_Ltf+Mmp9+"


for i in ['5','6']:
    preAnno3[i] = "Neu_Sell+S100a6+"
    
for i in ['2','9','3','4','11']:
    preAnno3[i] = "Neu_Sell+Csf3r+"


for i in ['10']:
    preAnno3[i] = "Neu_Cxcr4+Il1b+"

    
print(preAnno3)


adata.obs['preAnno3'] = adata.obs['leiden_res_1.00'].map(preAnno3)

adata.obs.loc[adata.obs['leiden_res_4.00'].isin(['41', '45','46','55']), 'preAnno3'] = "Neu_Cxcr4+Il1b+"


adata.obs['preAnno3'] = pd.Categorical(adata.obs['preAnno3'],categories=["Neu_Pcna+Mpo+","Neu_Pcna+Mki67+","Neu_Ltf+Camp","Neu_Ltf+Mmp9+","Neu_Sell+S100a6+", "Neu_Sell+Csf3r+","Neu_Cxcr4+Il1b+"])                   
           
adata.uns['preAnno3_colors'] = ["#76afda","#5066a1","#b0d45d","#7fb961","#ffe788","#ffc556","#f06152"] 
 
sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.png")
sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.pdf")




adata.obs['preAnno3'].unique()

preAnno4 = dict()

for i in ["Neu_Pcna+Mpo+","Neu_Pcna+Mki67+","Neu_Ltf+Camp","Neu_Ltf+Mmp9+","Neu_Sell+S100a6+", "Neu_Sell+Csf3r+","Neu_Cxcr4+Il1b+"]:
    preAnno4[str(i)] = "Undefined"

print(preAnno4)


####开始定义细胞类型
for i in ["Neu_Pcna+Mpo+","Neu_Pcna+Mki67+"]:
    preAnno4[i] = "Neu_Pcna+"


for i in ["Neu_Ltf+Camp","Neu_Ltf+Mmp9+"]:
    preAnno4[i] = "Neu_Ltf+"


for i in ["Neu_Sell+S100a6+", "Neu_Sell+Csf3r+"]:
    preAnno4[i] = "Neu_Sell+"
    

for i in ["Neu_Cxcr4+Il1b+"]:
    preAnno4[i] = "Neu_Il1b+"



print(preAnno4)


adata.obs['preAnno4'] = adata.obs['preAnno3'].map(preAnno4)

adata.obs['preAnno4'].unique()

adata.obs['preAnno4'] = pd.Categorical(adata.obs['preAnno4'],categories=["Neu_Pcna+","Neu_Ltf+","Neu_Sell+","Neu_Il1b+"])  

adata.uns['preAnno4_colors'] = ["#5066a1","#7fb961","#ffc556","#f06152"]  

sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.png")
sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.pdf")



adata.write_h5ad('05_neu_anno.h5ad',compression='gzip')





#——————————————————————堆积条形图————————————————————————————

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


grp = "preAnno3"
color= 'preAnno3_colors'


# 提取数据
df = adata.obs[['batch', grp]].copy()

# 计算每个batch中每个celltype的占比
counts = df.groupby(['batch', grp]).size().unstack(fill_value=0)
proportions = counts.div(counts.sum(axis=1), axis=0)


colors = adata.uns[color]

# 绘制堆积条形图
ax = proportions.plot(kind='barh', stacked=True, figsize=(9,4), color=colors,width=0.8)
# 去掉背景网格线
ax.grid(False)
plt.xlabel('Proportion')
plt.ylabel('Sample')
plt.title('Proportions within each sample')
font_properties = FontProperties(size=15) 
plt.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left',prop=font_properties)
#plt.xticks(rotation=45)
plt.tight_layout()  # 自动调整布局，以适应标签

plt.savefig('./figures/proportion_stacked_04-15.pdf')





import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


grp = "preAnno4"
color= 'preAnno4_colors'


# 提取数据
df = adata.obs[['batch', grp]].copy()

# 计算每个batch中每个celltype的占比
counts = df.groupby(['batch', grp]).size().unstack(fill_value=0)
proportions = counts.div(counts.sum(axis=1), axis=0)


colors = adata.uns[color]

# 绘制堆积条形图
ax = proportions.plot(kind='barh', stacked=True, figsize=(9,4), color=colors,width=0.8)
# 去掉背景网格线
ax.grid(False)
plt.xlabel('Proportion')
plt.ylabel('Sample')
plt.title('Proportions within each sample')
font_properties = FontProperties(size=15) 
plt.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left',prop=font_properties)
#plt.xticks(rotation=45)
plt.tight_layout()  # 自动调整布局，以适应标签

plt.savefig('./figures/proportion_stacked_04-16.pdf')







#————————————————————————————————————————————————————————————————————



adata.uns['phase_colors'] = ["#167CB4","#FE9F2A","#EA4737"]

sc.pl.umap(adata, color='phase',size=3,show=False, save="_04-17.png")




import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['S_score', 'G2M_score'],multi_panel=True,jitter=0.4,size=0.3, groupby="preAnno",show=False)

color = ["#FE9F2A","#EA4737"]


fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.grid(False)
    
    # 遍历每个小提琴图的图形对象
    for collection in ax.collections:
        collection.set_facecolor(color[i])
        
        # 处理小提琴图的背景颜色
        if isinstance(collection, plt.Polygon):
            collection.set_facecolor(color[i])
        # 保持散点颜色为黑色
        elif isinstance(collection, PathCollection):
            collection.set_edgecolor('black')
            collection.set_facecolor('black')
        

plt.savefig('./figures/violin_04-18.pdf', bbox_inches='tight', pad_inches=1)





import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],multi_panel=True,jitter=0.4,size=0.3, groupby="preAnno",show=False)


fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-19.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')






import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['Srgn','Clec4d','Jaml','Csf3r','Tpd52','Zfp36','C5ar1','Tgfbi','Dazap2','Ptprc','Adam8','Slpi','Pla2g7','Cd44','Glipr1'],multi_panel=True,size=0,groupby="preAnno",show=False)


fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-20.png', bbox_inches='tight', pad_inches=1)

plt.close('all')




import matplotlib.pyplot as plt

# 创建一个新的图形
fig, axes = plt.subplots(nrows=3, ncols=5, figsize=(25, 15),subplot_kw=dict(aspect='equal'))  

genes = ['Srgn','Clec4d','Jaml','Csf3r','Tpd52','Zfp36','C5ar1','Tgfbi','Dazap2','Ptprc','Adam8','Slpi','Pla2g7','Cd44','Glipr1']

for i, gene in enumerate(genes):
    row = i // 5
    col = i % 5
    ax = axes[row, col]
    
    sc.pl.violin(adata,keys=[gene],groupby="preAnno",size=0,ax=ax, show=False)
    
    ax.set_title(f'Gene {gene}')
    
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
  
  
plt.tight_layout() 
plt.savefig('./figures/violin_04-20.png', bbox_inches='tight', pad_inches=1)

plt.close('all')







#———————————————————————— DEG  ————————————————————————————
adata=sc.read_h5ad("05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 18784 × 2000


adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 18784 × 15938

print(adata_counts.X)


ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 18784 × 15938

print(adata_counts.X)


# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['batch', 'dataset', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'celltype', 'preAnno','preAnno3','preAnno4', 'S_score', 'G2M_score', 'phase']],
    var=adata_counts.var[['gene_ids']]
)



new_adata.layers["counts"] = new_adata.X.copy()

sc.pp.normalize_total(new_adata)
sc.pp.log1p(new_adata)


adata=new_adata



#移除特定的基因————————————————————————————————————————————
gene_names = adata.var_names

# 定义要去除的前缀
prefixes = ['mt-', 'Rps', 'Rpl', 'Hsp']

keep_mask = [not any(gene_name.startswith(prefix) for prefix in prefixes) for gene_name in gene_names]


# 获取被过滤掉的基因名称
filtered_genes = [gene_name for gene_name, keep in zip(gene_names, keep_mask) if not keep]
len(filtered_genes)
#136


# 使用布尔掩码过滤adata对象
adata2 = adata[:, keep_mask]
adata2
#View of AnnData object with n_obs × n_vars = 18784 × 15802


adata = adata2



sc.tl.rank_genes_groups(adata, groupby="preAnno4", method="wilcoxon",pts=True)


celltype=adata.obs['preAnno4'].unique()
deg=sc.get.rank_genes_groups_df(adata,group=celltype) 
deg.to_csv('./OVX_Neu_preAnno4_deg.csv') #存储备份




adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X



import matplotlib.pyplot as plt

sc.pl.rank_genes_groups_heatmap(adata,groupby="preAnno",  layer="scaled",n_genes=15,show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 12), show=False)


fig = plt.gcf()

for ax in fig.axes:
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    

plt.savefig('./figures/heatmap_04-20.png', bbox_inches='tight', pad_inches=1)


plt.close('all')





#——————————————————SELL————————————————————————————————

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colors = ["#FFF5F0", "#FFF5F0","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"]

positions = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975,1]
#positions = [0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 1]

cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))


marker_genes = ['Cd63','Cd37','Cd33','Sell','Cd53']
marker_genes = ['Sell']

sc.pl.dotplot(adata, marker_genes, groupby='preAnno4', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-12.png")

sc.pl.dotplot(adata, marker_genes, groupby='preAnno4', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-12.pdf")

