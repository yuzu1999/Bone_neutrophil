
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

os.chdir('/home/luchun/scRNA/Bone_neu/04_cell/neu/')

sc.settings.set_figure_params(dpi=50, facecolor="white")


adata=sc.read_h5ad("04_merge.h5ad")
adata
#AnnData object with n_obs × n_vars = 8997 × 2000





#—————————————————————— 想要剔除的细胞类型列表——————————————————
adata.obs['leiden_res_1.00'].value_counts()

adata
#AnnData object with n_obs × n_vars = 8997 × 2000


clutser_to_remove = ['15','4','1','10','0']
 
# 创建一个布尔掩码
mask = ~adata.obs['leiden_res_1.00'].isin(clutser_to_remove)
 
# 使用mask提取子集
adata_subset = adata[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 7048 × 2000


sc.pl.umap(adata_subset, color=['leiden_res_1.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_test.png")
sc.pl.umap(adata_subset, color=['HBA1','HBA2',"pct_counts_mt"],cmap="Reds",show=False,save="_test2.png")



adata_counts=adata_subset.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155


print(adata_counts.X)

ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155

print(adata_counts.X)




# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts','total_counts','pct_counts_mt']],
    var=adata_counts.var[['features']]
)



adata=new_adata


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
cell_cycle_genes = [x.strip() for x in open('/home/luchun/scRNA/LH/thymus/human/regev_lab_cell_cycle_genes.txt')]

s_genes = cell_cycle_genes[:43] #列表的开头开始，一直到第 42个元素（但不包括第 42 个元素）的所有元素。
g2m_genes = cell_cycle_genes[43:]

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

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="orig.ident")
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
prefixes = ['MT-', 'RPS', 'RPL', 'HSP']

# 更新 'highly_variable' 列的值
for prefix in prefixes:
    adata.var.loc[adata.var.index.str.startswith(prefix), 'highly_variable'] = False


adata.var['highly_variable'].value_counts()
'''


adata = adata[:, adata.var.highly_variable]
adata
#View of AnnData object with n_obs × n_vars = 7048 × 2000



#sc.pp.regress_out(adata, ['percent.mt'])
sc.pp.scale(adata, max_value=10)




 
#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_03-3.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)


sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2",'PTPRC','CSF1R',"VCAN","MPO",'ELANE',"CSF3R",'MKI67','PCNA','LTF','CAMP'],ncols=6,cmap="Reds",size=8,show=False,save="_03-4.png")

sc.pl.umap(adata,color=['HBA1','HBA2',"pct_counts_mt","phase"],cmap="Reds",show=False,save="_03-5.png")




#——————————————————————Clustering——————————————————————————————————
for res in [0.5,1,1.5,2,2.5,3]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
    
sc.pl.umap(adata, color=['leiden_res_0.50', 'leiden_res_1.00','leiden_res_1.50','leiden_res_2.00','leiden_res_2.50','leiden_res_3.00'],ncols=3,legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_03-6.png")

sc.pl.umap(adata, color=['leiden_res_1.50'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_03-7.png")




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


sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_0.50', standard_scale="var",dendrogram=True,show=False,save="_03-8.png")


marker_genes = {
   
    "Prog": ['ELANE', 'MPO',"PRTN3"],
    "Proliferating": ["MKI67","PCNA","TOP2A"],
    "Immature": ['LTF','CAMP','LCN2'],
    "Mature":["MME","SELL",'CXCR2','CSF3R'] ,
    "Age":['IL1B',"CXCL2","NFKBIA",'CXCR4']
}


sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_0.50', standard_scale="var",dendrogram=True,show=False,save="_03-9.png")



adata.write_h5ad('06_merge.h5ad',compression='gzip')


#adata=sc.read_h5ad("06_merge.h5ad")
#adata



#——————————————————————————bbknn——————————————————————————————————
sc.external.pp.bbknn(adata, batch_key="orig.ident") 
sc.tl.umap(adata)


sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2","MPO",'ELANE',"CSF3R",'MKI67','PCNA','LTF','CAMP','SELL','MME','CXCR4','IL1B','CEBPB'],ncols=5,cmap="Reds",size=8,show=False,save="_03-10.png")

sc.pl.umap(adata,color=['HBA1','HBA2',"pct_counts_mt","phase"],cmap="Reds",show=False,save="_03-11.png")





'''
#————————————————————————harmony————————————————————————————————
#https://support.parsebiosciences.com/hc/en-us/articles/7704577188500-How-to-analyze-a-1-million-cell-data-set-using-Scanpy-and-Harmony
import harmonypy
help(harmonypy.run_harmony)


import scanpy.external as sce

sce.pp.harmony_integrate(adata, key='orig.ident',theta=1,lamb=1)

adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)


sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2",'PTPRC','CSF1R',"VCAN","MPO",'ELANE',"CSF3R",'MKI67','PCNA','LTF','CAMP'],ncols=6,cmap="Reds",size=8,show=False,save="_02-5.png")

sc.pl.umap(adata,color=['HBA1','HBA2',"pct_counts_mt"],cmap="Reds",show=False,save="_02-6.png")

'''




#——————————————————————Clustering——————————————————————————————————
for res in [0.5,1,1.5,2,2.5,3]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
    
sc.pl.umap(adata, color=['leiden_res_0.50', 'leiden_res_1.00','leiden_res_1.50','leiden_res_2.00','leiden_res_2.50','leiden_res_3.00'],legend_loc='on data',ncols=3,legend_fontsize=8, legend_fontoutline=2,show=False,save="_03-12.png")

sc.pl.umap(adata, color=['leiden_res_1.50'],show=False,save="_03-13.png")




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


sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_2.00', standard_scale="var",dendrogram=True,show=False,save="_03-14.png")


marker_genes = {
   
    "Prog": ['ELANE', 'MPO',"PRTN3"],
    "Proliferating": ["MKI67","PCNA","TOP2A"],
    "Immature": ['LTF','CAMP','LCN2','CXCR4'],
    "Mature":["MME","SELL",'CXCR2','CSF3R'] ,
    "Age":['IL1B',"CXCL2","NFKBIA",'CEBPB']
}


sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_2.00', standard_scale="var",dendrogram=True,show=False,save="_03-15.png")





#——————————————————————定义发育阶段————————————————————————————
preAnno = dict()

for i in range(0,15):
    preAnno[str(i)] = "Undefined"

print(preAnno)


####开始定义细胞类型
for i in ['5','3','4','9']:
    preAnno[i] = "AZU1+KCNQ5+ progenitor"

for i in ['0']:
    preAnno[i] = "PDE4D+MCU+ proliferating"
    

for i in ['1','2']:
    preAnno[i] = "LTF+CAMP+ immature"

for i in ['10','6']:
    preAnno[i] = "S100A12+MMP9+ immature"
    
   
for i in ['12','11','14']:
    preAnno[i] = "S100A6+NCF1+ mature"

for i in ['8','13']:
    preAnno[i] = "NAMPT+IFITM2+ mature"
    
    
for i in ['7']:
    preAnno[i] = "NEAT1+FTH1+ mature"



  
print(preAnno)


adata.obs['preAnno'] = adata.obs['leiden_res_2.00'].map(preAnno)

adata.obs['preAnno'] = pd.Categorical(adata.obs['preAnno'],categories=["AZU1+KCNQ5+ progenitor","PDE4D+MCU+ proliferating","LTF+CAMP+ immature", "S100A12+MMP9+ immature","S100A6+NCF1+ mature","NAMPT+IFITM2+ mature","NEAT1+FTH1+ mature"])                   


adata.uns['preAnno_colors'] = ["#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272"] 
              
#sc.pl.umap(adata, color='preAnno',legend_loc='on data',legend_fontsize=4, legend_fontoutline=2,show=False, save="_03-16.png")

sc.pl.umap(adata, color='preAnno',show=False, save="_03-16.png")
sc.pl.umap(adata, color='preAnno',show=False, save="_03-16.pdf")


adata.uns['phase_colors'] = ["#167CB4","#FE9F2A","#EA4737"]

sc.pl.umap(adata, color='phase',show=False, save="_03-17.png")
sc.pl.umap(adata, color='phase',show=False, save="_03-17.pdf")


marker_genes = {
   
    "Prog": ['ELANE', 'MPO',"PRTN3"],
    "Proliferating": ["MKI67","PCNA","TOP2A"],
    "Immature": ['LTF','CAMP','LCN2','CXCR4'],
    "Mature":["MME",'CXCR2','CSF3R'] ,
    "Age":['IL1B',"CEBPB","NFKBIA",'S100A9','S100A6','NCF1'],
    "Gene":['CD63','LGALS3']
}


sc.pl.dotplot(adata, marker_genes, groupby='preAnno', standard_scale="var",dendrogram=False,show=False,save="_03-17.png")


sc.pl.umap(adata,color=['ELANE', 'MPO',"MKI67","TOP2A",'LTF','CAMP',"MME","SELL",'S100A12',"MMP9"],ncols=5,cmap="Reds",size=8,show=False,save="_03-18.png")

genes=['CD63','LGALS3','ELANE','MPO','AZU1','MKI67','TOP2A','PCNA','LTF','CAMP','LCN2','MMP8','CEBPE','CEACAM8','MMP9','ORM1','S100A12','IFIT1','CD177','MME','CXCR2','CSF3R','FCGR3B','ARG1','OLR1']
sc.pl.umap(adata,color=genes,ncols=5,cmap="Reds",size=8,show=False,save="_test.png")




#————————————细胞周期小提琴图——————————————————————————
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
        

plt.savefig('./figures/violin_03-19.pdf', bbox_inches='tight', pad_inches=1)




#————————————————基因数小提琴图————————————————————
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_03-20.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')



adata.write_h5ad('07_Anno_filter.h5ad',compression='gzip')

#adata=sc.read_h5ad("07_Anno_filter.h5ad")
#adata




#—————————————————————— DEG  ————————————————————————
adata=sc.read_h5ad("07_Anno_filter.h5ad")
adata


adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155


print(adata_counts.X)

ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155

print(adata_counts.X)




# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts','total_counts','pct_counts_mt','preAnno','preAnno2','preAnno3','preAnno4']],
    var=adata_counts.var[['features']]
)



adata=new_adata



#移除特定的基因————————————————————————————————————————————
gene_names = adata.var_names

# 定义要去除的前缀
prefixes = ['MT-', 'RPS', 'RPL', 'HSP']

keep_mask = [not any(gene_name.startswith(prefix) for prefix in prefixes) for gene_name in gene_names]


# 获取被过滤掉的基因名称
filtered_genes = [gene_name for gene_name, keep in zip(gene_names, keep_mask) if not keep]
len(filtered_genes)
#146


# 使用布尔掩码过滤adata对象
adata2 = adata[:, keep_mask]
adata2
#View of AnnData object with n_obs × n_vars = 7048 × 21009


adata = adata2



ov.utils.store_layers(adata,layers='counts')
adata.layers['counts']=adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)




sc.tl.rank_genes_groups(adata, groupby="preAnno4", method="wilcoxon",pts=True)


celltype=adata.obs['preAnno4'].unique()
deg=sc.get.rank_genes_groups_df(adata,group=celltype) 
deg.to_csv('./Human_Neu_preAnno4_deg.csv') #存储备份


adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X



#热图————————————————————————
import matplotlib.pyplot as plt

sc.pl.rank_genes_groups_heatmap(adata,groupby="preAnno",  layer="scaled",n_genes=15,show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 12), show=False)


fig = plt.gcf()

for ax in fig.axes:
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    

plt.savefig('./figures/heatmap_03-21.png', bbox_inches='tight', pad_inches=2)


plt.close('all')






#————————————————————————————————————————————————————————————————————
adata=sc.read_h5ad("07_Anno_filter.h5ad")
adata
#AnnData object with n_obs × n_vars = 7048 × 2000


adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155


print(adata_counts.X)


adata=adata_counts


adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X





#scprogram————————————————————————————
import csv
from collections import OrderedDict



df = pd.read_csv('Human_scprogram.csv')
group_genes_dict = OrderedDict((col, df[col].head(20).tolist()) for col in df.columns)



sc.pl.matrixplot(adata,group_genes_dict,groupby="preAnno", layer="scaled",cmap="RdBu_r",vmin=-2,vmax=2,show=False,save='_03-20.png')

sc.pl.dotplot(adata,group_genes_dict, groupby='preAnno', standard_scale="var",dendrogram=False,show=False,save="_03-20.pdf")




import matplotlib.pyplot as plt 
 
sc.pl.heatmap(adata,group_genes_dict,groupby="preAnno",  layer="scaled",show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 12), show=False)


fig = plt.gcf()

for ax in fig.axes:
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    
    
    
plt.savefig('./figures/heatmap_03-20.svg', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_03-20.png', bbox_inches='tight', pad_inches=1)







#DEG——————————————————————————————————

df = pd.read_csv('Neu_preAnno_deg.csv')

filtered_df = df[(df['logfoldchanges'] > 0) & (df['pvals_adj'] < 0.05)]

# 指定group的顺序
group_order = ["AZU1+KCNQ5+ progenitor","PDE4D+MCU+ proliferating","LTF+CAMP+ immature", "S100A12+MMP9+ immature","S100A6+NCF1+ mature","NAMPT+IFITM2+ mature","NEAT1+FTH1+ mature"]


group_genes_dict = OrderedDict()


grouped = filtered_df.groupby('group')

# 对每个组进行处理
for group in group_order:
    if group in grouped.groups:
        group_df = grouped.get_group(group)    
        top_genes = group_df.nlargest(20, 'logfoldchanges')['names'].tolist()
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

    
    
    
plt.savefig('./figures/heatmap_03-20.svg', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_03-20.png', bbox_inches='tight', pad_inches=1)





# specific DEG————————————————————————
df = pd.read_excel('Human_Neu_DEG.xlsx')

all_genes = set()

gene_dict = {}
for col in df.columns:
    genes = df[col].dropna().head(20).tolist()
    unique_genes = [gene for gene in genes if gene not in all_genes]
    all_genes.update(unique_genes)
    gene_dict[col] = unique_genes



print(gene_dict)



sc.pl.matrixplot(adata,gene_dict,groupby="preAnno", layer="scaled",cmap="RdBu_r",vmin=-2,vmax=2,show=False,save='_03-20.png')

sc.pl.dotplot(adata,gene_dict, groupby='preAnno', standard_scale="var",dendrogram=False,show=False,save="_03-20.pdf")




import matplotlib.pyplot as plt

sc.pl.heatmap(adata,gene_dict,groupby="preAnno",  layer="scaled",show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 12), show=False)


fig = plt.gcf()

for ax in fig.axes:
    ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    
    
    
plt.savefig('./figures/heatmap_03-20.svg', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_03-20.png', bbox_inches='tight', pad_inches=1)









#scprogram————————————————————————————

df = pd.read_csv('Human_scprogram.csv')

total_genes = []

# 遍历每一列，提取前50个基因并添加到列表中
for col in df.columns:
    total_genes.extend(df[col].head(20).tolist())


print(total_genes)


len(total_genes)
#140

total_genes

unique_genes=[]

for i in total_genes:
    if not i in unique_genes:
        unique_genes.append(i)


len(unique_genes)
#130

sc.settings.set_figure_params(dpi=10, facecolor="white")
sc.pl.umap(adata,color=unique_genes,ncols=10,cmap="Reds",size=5,show=False,save="_03-21.pdf")







#DEG——————————————————————————————————
df = pd.read_csv('Neu_preAnno_deg.csv')



group_order = ["AZU1+KCNQ5+ progenitor","PDE4D+MCU+ proliferating","LTF+CAMP+ immature", "S100A12+MMP9+ immature","S100A6+NCF1+ mature","NAMPT+IFITM2+ mature","NEAT1+FTH1+ mature"]



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
#110


sc.settings.set_figure_params(dpi=10, facecolor="white")
sc.pl.umap(adata,color=unique_genes,ncols=10,cmap="Reds",size=5,show=False,save="_03-21.pdf")





#————————————————————————————————————————————————————————————

import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['CXCL8','B2M','G0S2','C5AR1','LITAF','SRGN','MME','HLA-C','PLAUR','VNN2'],multi_panel=True,size=0,groupby="preAnno",show=False)


fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_03-22.png', bbox_inches='tight', pad_inches=1)

plt.savefig('./figures/violin_03-22.pdf', bbox_inches='tight', pad_inches=1)




sc.pl.violin(adata, ['MKI67', 'ANP32B', 'CCNE2', 'TOP2A', 'CENPF', 'SMC4'],multi_panel=True,size=0,groupby="preAnno",show=False)


fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_03-23.png', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/violin_03-23.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')


sc.pl.dotplot(adata, ['MKI67', 'ANP32B', 'CCNE2', 'TOP2A', 'CENPF', 'SMC4'], groupby='preAnno', standard_scale="var",dendrogram=False,show=False,save="_03-24.png")
sc.pl.dotplot(adata, ['MKI67', 'ANP32B', 'CCNE2', 'TOP2A', 'CENPF', 'SMC4'], groupby='preAnno', standard_scale="var",dendrogram=False,show=False,save="_03-24.pdf")






#—————————————————— 重新命名 ————————————————————————


adata=sc.read_h5ad("07_Anno_filter.h5ad")
adata


adata.obs['preAnno'].unique()

preAnno2 = dict()

for i in ['AZU1+KCNQ5+ progenitor', 'PDE4D+MCU+ proliferating', 'LTF+CAMP+ immature','S100A12+MMP9+ immature', 'S100A6+NCF1+ mature', 'NAMPT+IFITM2+ mature','NEAT1+FTH1+ mature']:
    preAnno2[str(i)] = "hNeu-1"


print(preAnno2)


####开始定义细胞类型
for i in ['PDE4D+MCU+ proliferating']:
    preAnno2[i] = "hNeu-2"


for i in ['LTF+CAMP+ immature']:
    preAnno2[i] = "hNeu-3"



for i in ['S100A12+MMP9+ immature']:
    preAnno2[i] = "hNeu-4"


for i in ['S100A6+NCF1+ mature']:
    preAnno2[i] = "hNeu-5"
    
for i in ['NAMPT+IFITM2+ mature']:
    preAnno2[i] = "hNeu-6"


for i in ['NEAT1+FTH1+ mature']:
    preAnno2[i] = "hNeu-7"



print(preAnno2)


adata.obs['preAnno2'] = adata.obs['preAnno'].map(preAnno2)

adata.obs['preAnno2'].unique()
adata.uns['preAnno2_colors'] = ["#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272"]  
sc.pl.umap(adata, color='preAnno2',show=False, save="_test.png")






sc.pl.umap(adata,color=['CD34','KIT','MPO','ELANE','AZU1','MKI67','TOP2A','PCNA','LTF','CAMP','ISG15','IFIT1','MMP8','MMP9','FCGR3B','SELL','IFITM2','CXCR4','IL1B','CSF3R','CXCR2','S100A6','S100A12','FPR1'],ncols=5,cmap="Reds",size=8,show=False,save="_test.png")




adata.obs['preAnno'].unique()

preAnno3 = dict()

for i in ['AZU1+KCNQ5+ progenitor', 'PDE4D+MCU+ proliferating', 'LTF+CAMP+ immature','S100A12+MMP9+ immature', 'S100A6+NCF1+ mature', 'NAMPT+IFITM2+ mature','NEAT1+FTH1+ mature']:
    preAnno3[str(i)] = "Undefined"


print(preAnno3)


####开始定义细胞类型
for i in ['AZU1+KCNQ5+ progenitor']:
    preAnno3[i] = "Neu_PCNA+MPO+"


for i in ['PDE4D+MCU+ proliferating']:
    preAnno3[i] = "Neu_PCNA+MKI67+"


for i in ['LTF+CAMP+ immature']:
    preAnno3[i] = "Neu_LTF+CAMP+"


for i in ['S100A12+MMP9+ immature']:
    preAnno3[i] = "Neu_LTF+MMP9+"


for i in ['S100A6+NCF1+ mature']:
    preAnno3[i] = "Neu_SELL+FPR1+"
    
for i in ['NAMPT+IFITM2+ mature']:
    preAnno3[i] = "Neu_SELL+FCGR3B+"


for i in ['NEAT1+FTH1+ mature']:
    preAnno3[i] = "Neu_CXCR4+IL1B+"



print(preAnno3)


adata.obs['preAnno3'] = adata.obs['preAnno'].map(preAnno3)

adata.obs['preAnno3'].unique()

adata.obs['preAnno3'] = pd.Categorical(adata.obs['preAnno3'],categories=["Neu_PCNA+MPO+","Neu_PCNA+MKI67+","Neu_LTF+CAMP+","Neu_LTF+MMP9+","Neu_SELL+FPR1+","Neu_SELL+FCGR3B+","Neu_CXCR4+IL1B+"])  

#adata.uns['preAnno3_colors'] = ["#5E4FA2","#3288BD","#ffe788","#ffc556","#F46D43","#D53E4F","#9E0142"]  

adata.uns['preAnno3_colors'] = ["#76afda","#5066a1","#b0d45d","#7fb961","#ffe788","#ffc556","#f06152"] 

sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.png")
sc.pl.umap(adata, color='preAnno3',show=False, save="_01-1.pdf")






adata.obs['preAnno'].unique()

preAnno4 = dict()

for i in ['AZU1+KCNQ5+ progenitor', 'PDE4D+MCU+ proliferating', 'LTF+CAMP+ immature','S100A12+MMP9+ immature', 'S100A6+NCF1+ mature', 'NAMPT+IFITM2+ mature','NEAT1+FTH1+ mature']:
    preAnno4[str(i)] = "Undefined"


print(preAnno4)


####开始定义细胞类型
for i in ['AZU1+KCNQ5+ progenitor','PDE4D+MCU+ proliferating']:
    preAnno4[i] = "Neu_PCNA+"


for i in ['LTF+CAMP+ immature','S100A12+MMP9+ immature']:
    preAnno4[i] = "Neu_LTF+"


for i in ['S100A6+NCF1+ mature','NAMPT+IFITM2+ mature']:
    preAnno4[i] = "Neu_SELL+"
    

for i in ['NEAT1+FTH1+ mature']:
    preAnno4[i] = "Neu_IL1B+"



print(preAnno4)


adata.obs['preAnno4'] = adata.obs['preAnno'].map(preAnno4)

adata.obs['preAnno4'].unique()

adata.obs['preAnno4'] = pd.Categorical(adata.obs['preAnno4'],categories=["Neu_PCNA+","Neu_LTF+","Neu_SELL+","Neu_IL1B+"])  

adata.uns['preAnno4_colors'] = ["#5066a1","#7fb961","#ffc556","#f06152"]  
sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.png")
sc.pl.umap(adata, color='preAnno4',show=False, save="_01-2.pdf")



adata.write_h5ad('07_Anno_filter.h5ad',compression='gzip')



adata = sc.read_h5ad("07_Anno_filter.h5ad")
adata




import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colors = ["#FFF5F0", "#FFF5F0","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"]

positions = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975,1]
#positions = [0, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 1]

cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))


marker_genes = ['SELL']

sc.pl.dotplot(adata, marker_genes, groupby='preAnno4', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-12.png")

sc.pl.dotplot(adata, marker_genes, groupby='preAnno4', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-12.pdf")





import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colors = ["#FFF5F0", "#FFF5F0","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"]

positions = [0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1]

cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))


marker_genes = ['PCNA','MPO','MKI67','LTF','CAMP','MMP9','SELL','FPR1','S100A6','FCGR3B','CXCR4','IL1B']

sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-13.png")
sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-13.pdf")



import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colors = ["#FFF5F0", "#FFF5F0","#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"]

positions = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975,1]

cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, colors)))



marker_genes = {
   
    "Inflammation": ['CXCL8', 'IL1B','IL1R1','IL1R2','CXCR4','TREM1','FCGR3B','TNFRSF1A','TNFRSF1B','CEBPB'],
    "Oxidative stress": ['NFE2L2','TXNIP','FTH1','OXSR1','SAT1','CYB5R4','RNASET2','SELENOK',"SLC6A6","FOXO3"],
    "Apoptosis": ['G0S2','FOXO3','BID','STAT3','BTG2']
}


sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-14.png")
sc.pl.dotplot(adata, marker_genes, groupby='preAnno3', standard_scale="var",dendrogram=False, cmap=cmap,show=False,swap_axes=True,save="_03-14.pdf")




#————————————————————导出——————————————————————————————

adata = sc.read_h5ad("07_Anno_filter.h5ad")
adata
#AnnData object with n_obs × n_vars = 7048 × 2000


adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155
print(adata_counts.X)


ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155
print(adata_counts.X)

import diopy
diopy.output.write_h5(adata_counts,'07_Anno_filter.h5')