
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



adata_all = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/03_age_ovx2/06_Neu_anno.h5ad")
adata_all
#AnnData object with n_obs × n_vars = 53904 × 2000




#—————————————————————— 想要的细胞类型列表——————————————————
adata_all.obs['batch'].value_counts()


clutser_to_stay = ['4W', '12W', '24M']
 
# 创建一个布尔掩码
mask = adata_all.obs['batch'].isin(clutser_to_stay)
 
# 使用mask提取子集
adata_subset = adata_all[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 34580 × 2000



adata_counts=adata_subset.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 34580 × 18544
print(adata_counts.X)



ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 34580 × 18544
print(adata_counts.X)



# 使用counts的layer重新构建adata
adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['batch', 'dataset', 'n_genes_by_counts','total_counts','pct_counts_mt','celltype','preAnno']],
    var=adata_counts.var[['gene_ids']]
)
adata
#AnnData object with n_obs × n_vars = 34580 × 18544



sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 34580 × 17769


adata.obs['batch'].value_counts()
adata.obs['dataset'].value_counts()

adata.obs['batch'] = pd.Categorical(adata.obs['batch'],categories=['4W', '12W', '24M'])






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
sc.pl.pca_scatter(adata_cc_genes, color='phase',show=False,save="_03-1.png")




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
sc.pl.highly_variable_genes(adata,show=False,save="_03-2.png")


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
#View of AnnData object with n_obs × n_vars = 34580 × 2000


sc.pp.scale(adata, max_value=10)



 
 
#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_03-3.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)

sc.pl.umap(adata,color=['batch','phase',"preAnno"],cmap="Reds",size=5,show=False,save="_03-4.png")

                       
adata.uns['preAnno_colors'] = ["#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272"]        


sc.pl.umap(adata, color='preAnno',size=5,show=False, save="_03-5.png")
sc.pl.umap(adata, color='preAnno',size=5,show=False, save="_03-5.pdf")


adata.obs['batch'].unique()
adata.obs['batch'].value_counts()


sc.tl.embedding_density(adata, groupby="batch")
sc.pl.embedding_density(adata, groupby="batch",ncols=3,show=False, save="_03-6.png")

sc.pl.umap(adata,color=['Elane','Prtn3','Mpo',"Mki67","Top2a",'H2afz','Hmgb2','Ltf','Camp','Ngp','Lcn2','Retnlg','Mmp8',"Sell","Cxcr2",'Ccl6','Srgn','Il1b',"Cxcr4"],ncols=5,cmap="Reds",size=8,show=False,save="_03-7.png")

sc.pl.umap(adata,color=['Ly6g'],ncols=5,cmap="Reds",size=8,show=False,save="_03-7.png")


adata.write_h5ad('05_neu_anno.h5ad',compression='gzip')


#adata = sc.read_h5ad("05_neu_anno.h5ad")
#adata




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

plt.savefig('./figures/proportion_stacked_03-7.pdf')







#———————————————————————— DEG  ————————————————————————————
adata=sc.read_h5ad("05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 34580 × 2000

adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 34580 × 17769
print(adata_counts.X)



ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 34580 × 17769
print(adata_counts.X)




# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['batch', 'dataset', 'n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'celltype', 'preAnno','preAnno2','preAnno3','preAnno4', 'S_score', 'G2M_score', 'phase']],
    var=adata_counts.var[['gene_ids']]
)



adata=new_adata



#移除特定的基因————————————————————————————————————————————
gene_names = adata.var_names

# 定义要去除的前缀
prefixes = ['mt-', 'Rps', 'Rpl', 'Hsp']

keep_mask = [not any(gene_name.startswith(prefix) for prefix in prefixes) for gene_name in gene_names]


# 获取被过滤掉的基因名称
filtered_genes = [gene_name for gene_name, keep in zip(gene_names, keep_mask) if not keep]
len(filtered_genes)
#137


# 使用布尔掩码过滤adata对象
adata2 = adata[:, keep_mask]
adata2
#View of AnnData object with n_obs × n_vars = 34580 × 17632


adata = adata2



ov.utils.store_layers(adata,layers='counts')
adata.layers['counts']=adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)




sc.tl.rank_genes_groups(adata, groupby="preAnno4", method="wilcoxon",pts=True)


celltype=adata.obs['preAnno4'].unique()
deg=sc.get.rank_genes_groups_df(adata,group=celltype) 
deg.to_csv('./Age_Neu_preAnno4_deg.csv') #存储备份




adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X



import matplotlib.pyplot as plt

sc.pl.rank_genes_groups_heatmap(adata,groupby="preAnno",  layer="scaled",n_genes=15,show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 12), show=False)


fig = plt.gcf()

for ax in fig.axes:
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    

plt.savefig('./figures/heatmap_03-8.png', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_03-8.pdf', bbox_inches='tight', pad_inches=1)


plt.close('all')




#——————————————————  热图  ————————————————————————————
adata=sc.read_h5ad("05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 34580 × 2000


adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 34580 × 17769

print(adata_counts.X)


adata=adata_counts


adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X


from collections import OrderedDict


df = pd.read_csv('Age_Neu_preAnno_deg.csv')

filtered_df = df[(df['logfoldchanges'] > 0) & (df['pvals_adj'] < 0.05)]

filtered_df.to_csv('Age_Neu_preAnno_filtered_deg.csv', index=False)




# 指定group的顺序
group_order = ['Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature']


group_genes_dict = OrderedDict()


grouped = df.groupby('group')

# 对每个组进行处理
for group in group_order:
    if group in grouped.groups:
        group_df = grouped.get_group(group)    
        top_genes = group_df.nlargest(20, 'scores')['names'].tolist()
        group_genes_dict[group] = top_genes
    else:
        group_genes_dict[group] = []  
        
        

print(group_genes_dict)


import matplotlib.pyplot as plt

sc.pl.heatmap(adata,group_genes_dict,groupby="preAnno",  layer="scaled",show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 12), show=False)


fig = plt.gcf()

for ax in fig.axes:
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    
plt.savefig('./figures/heatmap_03-8.svg',dpi=300, bbox_inches='tight',pad_inches=1)


plt.savefig('./figures/heatmap_03-8.png', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_03-8.pdf', bbox_inches='tight', pad_inches=1)







df = pd.read_csv('Age_Neu_preAnno_deg.csv')


group_order = [
    'Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature'
]



total_genes = list()


for group in group_order:
    group_df = df[df['group'] == group]
    top_genes= group_df.nlargest(20, 'scores')['names'].tolist()
    total_genes.extend(top_genes)
    
    
len(total_genes)
#140

total_genes

unique_genes=[]

for i in total_genes:
    if not i in unique_genes:
        unique_genes.append(i)

len(unique_genes)
#113


sc.settings.set_figure_params(dpi=10, facecolor="white")
sc.pl.umap(adata,color=unique_genes,ncols=10,cmap="Reds",size=5,show=False,save="_03-9.pdf")





#—————————————————   细胞周期 ———————————————————————————
adata=sc.read_h5ad("05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 34580 × 2000


adata.uns['phase_colors'] = ["#167CB4","#FE9F2A","#EA4737"]

sc.pl.umap(adata, color='phase',size=3,show=False, save="_03-10.png")
sc.pl.umap(adata, color='phase',size=3,show=False, save="_03-10.pdf")



import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['S_score', 'G2M_score'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

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
        

plt.savefig('./figures/violin_03-10.pdf', bbox_inches='tight', pad_inches=1)




#—————————————————   基因数量 ———————————————————————————
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_03-11.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')




#—————————————————— 重新命名 ————————————————————————

adata=sc.read_h5ad("05_neu_anno.h5ad")
adata


adata.obs['preAnno'].unique()

preAnno2 = dict()

for i in ['Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature']:
    preAnno2[str(i)] = "mNeu-1"


print(preAnno2)


####开始定义细胞类型
for i in ['H2afz+Hmgb2+ proliferating']:
    preAnno2[i] = "mNeu-2"


for i in ['Ngp+Lcn2+ immature']:
    preAnno2[i] = "mNeu-3"



for i in ['Ifitm6+Ltf+ immature']:
    preAnno2[i] = "mNeu-4"


for i in ['Retnlg+Mmp8+ mature']:
    preAnno2[i] = "mNeu-5"
    
for i in ['Ccl6+Sell+ mature']:
    preAnno2[i] = "mNeu-6"


for i in ['Il1b+Srgn+ mature']:
    preAnno2[i] = "mNeu-7"



print(preAnno2)


adata.obs['preAnno2'] = adata.obs['preAnno'].map(preAnno2)

adata.obs['preAnno2'].unique()
adata.uns['preAnno2_colors'] = ["#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272"]  
sc.pl.umap(adata, color='preAnno2',show=False, save="_test.png")






sc.pl.umap(adata,color=['Mpo','Elane','Mki67','Top2a','Pcna','Ltf','Camp','Isg15','Ifit1','Mmp8','Mmp9','Fcgr3','Sell','Ifitm2','Cxcr4','Il1b','Csf3r','Cxcr2','S100a6'],ncols=5,cmap="Reds",size=8,show=False,save="_test.png")




adata.obs['preAnno'].unique()

preAnno3 = dict()

for i in ['Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature']:
    preAnno3[str(i)] = "Undefined"


print(preAnno3)


####开始定义细胞类型
for i in ['Elane+Prtn3+ progenitor']:
    preAnno3[i] = "Neu_Pcna+Mpo+"


for i in ['H2afz+Hmgb2+ proliferating']:
    preAnno3[i] = "Neu_Pcna+Mki67+"


for i in ['Ngp+Lcn2+ immature']:
    preAnno3[i] = "Neu_Ltf+Camp"


for i in ['Ifitm6+Ltf+ immature']:
    preAnno3[i] = "Neu_Ltf+Mmp9+"


for i in ['Retnlg+Mmp8+ mature']:
    preAnno3[i] = "Neu_Sell+S100a6+"
    
for i in ['Ccl6+Sell+ mature']:
    preAnno3[i] = "Neu_Sell+Csf3r+"


for i in ['Il1b+Srgn+ mature']:
    preAnno3[i] = "Neu_Cxcr4+Il1b+"



print(preAnno3)


adata.obs['preAnno3'] = adata.obs['preAnno'].map(preAnno3)

adata.obs['preAnno3'].unique()

adata.obs['preAnno3'] = pd.Categorical(adata.obs['preAnno3'],categories=["Neu_Pcna+Mpo+","Neu_Pcna+Mki67+","Neu_Ltf+Camp","Neu_Ltf+Mmp9+","Neu_Sell+S100a6+", "Neu_Sell+Csf3r+","Neu_Cxcr4+Il1b+"])  

adata.uns['preAnno3_colors'] = ["#76afda","#5066a1","#b0d45d","#7fb961","#ffe788","#ffc556","#f06152"] 
 
sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.png")
sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.pdf")






adata.obs['preAnno'].unique()

preAnno4 = dict()

for i in ['Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature']:
    preAnno4[str(i)] = "Undefined"


print(preAnno4)


####开始定义细胞类型
for i in ['Elane+Prtn3+ progenitor','H2afz+Hmgb2+ proliferating']:
    preAnno4[i] = "Neu_Pcna+"


for i in ['Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature']:
    preAnno4[i] = "Neu_Ltf+"


for i in ['Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature']:
    preAnno4[i] = "Neu_Sell+"
    

for i in ['Il1b+Srgn+ mature']:
    preAnno4[i] = "Neu_Il1b+"



print(preAnno4)


adata.obs['preAnno4'] = adata.obs['preAnno'].map(preAnno4)

adata.obs['preAnno4'].unique()

adata.obs['preAnno4'] = pd.Categorical(adata.obs['preAnno4'],categories=["Neu_Pcna+","Neu_Ltf+","Neu_Sell+","Neu_Il1b+"])  

adata.uns['preAnno4_colors'] = ["#5066a1","#7fb961","#ffc556","#f06152"]  

sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.png")
sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.pdf")



adata.write_h5ad('05_neu_anno.h5ad',compression='gzip')



#—————————————————— 重新命名 ————————————————————————

adata=sc.read_h5ad("05_neu_anno.h5ad")
adata


adata.obs['preAnno'].unique()

preAnno2 = dict()

for i in ['Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature']:
    preAnno2[str(i)] = "mNeu-1"


print(preAnno2)


####开始定义细胞类型
for i in ['H2afz+Hmgb2+ proliferating']:
    preAnno2[i] = "mNeu-2"


for i in ['Ngp+Lcn2+ immature']:
    preAnno2[i] = "mNeu-3"



for i in ['Ifitm6+Ltf+ immature']:
    preAnno2[i] = "mNeu-4"


for i in ['Retnlg+Mmp8+ mature']:
    preAnno2[i] = "mNeu-5"
    
for i in ['Ccl6+Sell+ mature']:
    preAnno2[i] = "mNeu-6"


for i in ['Il1b+Srgn+ mature']:
    preAnno2[i] = "mNeu-7"



print(preAnno2)


adata.obs['preAnno2'] = adata.obs['preAnno'].map(preAnno2)

adata.obs['preAnno2'].unique()
adata.uns['preAnno2_colors'] = ["#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272"]  
sc.pl.umap(adata, color='preAnno2',show=False, save="_test.png")






sc.pl.umap(adata,color=['Mpo','Elane','Mki67','Top2a','Pcna','Ltf','Camp','Isg15','Ifit1','Mmp8','Mmp9','Fcgr3','Sell','Ifitm2','Cxcr4','Il1b','Csf3r','Cxcr2','S100a6','Fpr1'],ncols=5,cmap="Reds",size=8,show=False,save="_test.png")




adata.obs['preAnno'].unique()

preAnno3 = dict()

for i in ['Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature']:
    preAnno3[str(i)] = "Undefined"


print(preAnno3)


####开始定义细胞类型
for i in ['Elane+Prtn3+ progenitor']:
    preAnno3[i] = "Neu_Pcna+Mpo+"


for i in ['H2afz+Hmgb2+ proliferating']:
    preAnno3[i] = "Neu_Pcna+Mki67+"


for i in ['Ngp+Lcn2+ immature']:
    preAnno3[i] = "Neu_Ltf+Camp"


for i in ['Ifitm6+Ltf+ immature']:
    preAnno3[i] = "Neu_Ltf+Mmp9+"


for i in ['Retnlg+Mmp8+ mature']:
    preAnno3[i] = "Neu_Sell+Fpr1+"
    
for i in ['Ccl6+Sell+ mature']:
    preAnno3[i] = "Neu_Sell+Fcgr3+"


for i in ['Il1b+Srgn+ mature']:
    preAnno3[i] = "Neu_Cxcr4+Il1b+"



print(preAnno3)


adata.obs['preAnno3'] = adata.obs['preAnno'].map(preAnno3)

adata.obs['preAnno3'].unique()

adata.obs['preAnno3'] = pd.Categorical(adata.obs['preAnno3'],categories=["Neu_Pcna+Mpo+","Neu_Pcna+Mki67+","Neu_Ltf+Camp","Neu_Ltf+Mmp9+","Neu_Sell+Fpr1+", "Neu_Sell+Fcgr3+","Neu_Cxcr4+Il1b+"])  

adata.uns['preAnno3_colors'] = ["#76afda","#5066a1","#b0d45d","#7fb961","#ffe788","#ffc556","#f06152"] 
 
sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.png")
sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.pdf")






adata.obs['preAnno'].unique()

preAnno4 = dict()

for i in ['Elane+Prtn3+ progenitor', 'H2afz+Hmgb2+ proliferating', 'Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature', 'Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature','Il1b+Srgn+ mature']:
    preAnno4[str(i)] = "Undefined"


print(preAnno4)


####开始定义细胞类型
for i in ['Elane+Prtn3+ progenitor','H2afz+Hmgb2+ proliferating']:
    preAnno4[i] = "Neu_Pcna+"


for i in ['Ngp+Lcn2+ immature','Ifitm6+Ltf+ immature']:
    preAnno4[i] = "Neu_Ltf+"


for i in ['Retnlg+Mmp8+ mature', 'Ccl6+Sell+ mature']:
    preAnno4[i] = "Neu_Sell+"
    

for i in ['Il1b+Srgn+ mature']:
    preAnno4[i] = "Neu_Il1b+"



print(preAnno4)


adata.obs['preAnno4'] = adata.obs['preAnno'].map(preAnno4)

adata.obs['preAnno4'].unique()

adata.obs['preAnno4'] = pd.Categorical(adata.obs['preAnno4'],categories=["Neu_Pcna+","Neu_Ltf+","Neu_Sell+","Neu_Il1b+"])  

adata.uns['preAnno4_colors'] = ["#5066a1","#7fb961","#ffc556","#f06152"]  

sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.png")
sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.pdf")



adata.write_h5ad('05_neu_anno.h5ad',compression='gzip')





#————————————————————————————————————————————————————————————


clutser_to_stay = ["mNeu-3","mNeu-4","mNeu-5","mNeu-6","mNeu-7"]
 
# 创建一个布尔掩码
mask = adata.obs['preAnno2'].isin(clutser_to_stay)
 
# 使用mask提取子集
adata_subset = adata[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 30128 × 2000


import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection



sc.pl.violin(adata_subset, ['Cd63','Cd47','Cd177','Cd55','Cd37','Cd164'],multi_panel=True,jitter=0.4,size=0.3,groupby="preAnno2",show=False)


fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_03-12.png', bbox_inches='tight', pad_inches=1)





sc.pl.violin(adata_subset, ['Cd33','Cd52','Cd9','Sell','Il1r2','Cd53'],multi_panel=True,jitter=0.4,size=0.3,groupby="preAnno2",show=False)


fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_03-12-2.png', bbox_inches='tight', pad_inches=1)




adata=sc.read_h5ad("05_neu_anno.h5ad")
adata





import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colors = ["#FFF5F0", "#FFF5F0","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"]

positions = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975,1]
#positions = [0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 1]

cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))


marker_genes = ['Sell']

sc.pl.dotplot(adata, marker_genes, groupby='preAnno4', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-12.png")

sc.pl.dotplot(adata, marker_genes, groupby='preAnno4', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-12.pdf")



import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colors = ["#FFF5F0", "#FFF5F0","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"]

positions = [0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1]

cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))


marker_genes = ['Pcna','Mpo','Mki67','Ltf','Camp','Mmp9','Sell','Fpr1','S100a6','Fcgr3','Cxcr4','Il1b']

sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-13.png")
sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-13.pdf")





import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colors = ["#FFF5F0", "#FFF5F0","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"]

positions = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975,1]

cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))



marker_genes = {
   
    "Inflammation": ['Il1b','Il1r2','Cxcl2','Trem1','C5ar1','C5ar2','Cd14','Tnfrsf1a','Ccr1','Cxcr4'],
    "Oxidative stress": ['Nfe2l2','Txnip','Fth1','Ftl1','Srgn','Msrb1','Hmox1','Slc7a11','Osgin1','Slc40a1'],
    "Apoptosis": ['Mcl1','Bcl2l11','Foxo3','Cdkn1b','Bcl3']
}


sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-14.png")
sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-14.pdf")





#————————————————————导出——————————————————————————————

adata = sc.read_h5ad("05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 34580 × 2000


adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 34580 × 17769
print(adata_counts.X)


ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 34580 × 17769
print(adata_counts.X)

import diopy
diopy.output.write_h5(adata_counts,'05_neu_anno.h5')









