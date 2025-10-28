
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

os.chdir('/home/luchun/scRNA/Bone_neu/04_cell/')

sc.settings.set_figure_params(dpi=50, facecolor="white")



adata = sc.read_h5ad("03_remove_contamin.h5ad")
adata
#AnnData object with n_obs × n_vars = 80474 × 29188


sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 80474 × 29188



#——————————————————————Normalization—————————————————————————————
#存储原始数据以便后续还原
ov.utils.store_layers(adata,layers='counts')
adata.layers['counts']=adata.X.copy()


sc.pp.normalize_total(adata)
sc.pp.log1p(adata)




#———————————————————Feature selection———————————————————————————
adata.raw = adata 

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="orig.ident")
sc.pl.highly_variable_genes(adata,show=False,save="_03-1.png")


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
#False    17796
#True      1991
'''


adata = adata[:, adata.var.highly_variable]
adata
#View of AnnData object with n_obs × n_vars = 80474 × 2000


#sc.pp.regress_out(adata, ['percent.mt'])
sc.pp.scale(adata, max_value=10)





#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_03-2.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.pl.umap(adata,color=['orig.ident',"cluster_anno_l2",'PTPRC','CSF1R',"CD34","MPO",'ELANE',"CSF3R"],ncols=4,cmap="Reds",size=8,show=False,save="_03-3.png")

sc.pl.umap(adata,color=["cluster_anno_l2"],legend_loc='on data',legend_fontsize=4, legend_fontoutline=2,show=False,save="_03-4.png")




#——————————————————————Clustering——————————————————————————————————
for res in [1.5, 2,2.5, 3,3.5,4]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
    
    
sc.pl.umap(adata, color=['leiden_res_1.50','leiden_res_2.00','leiden_res_2.50','leiden_res_3.00','leiden_res_3.50','leiden_res_4.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,ncols=3,show=False,save="_03-5.png")


sc.pl.umap(adata, color=[ 'leiden_res_2.50'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_03-6.png")

sc.pl.umap(adata, color=[ 'leiden_res_2.50'],show=False,save="_03-7.png")



#——————————————————Manual cell-type annotation——————————————————
marker_genes = {
    "Hematopoietic": ["PTPRC"],
    "HSPC": ['KIT','CD34','SPINK2', 'ADGRG1', 'CDK6'],
    "Lymp_prog": ["FLT3"],
    "Ery_prog": ["GATA1",'TFRC'],
    "Erythro":['HBA1','HBA2' ], 
    
    "Mega_prog": ["ITGA2B"],
    "Mye_prog": ['CSF1R','ELANE', 'MPO',"PRTN3"],
    "Eos/baso_prog": ["MS4A2"],
    "Neutrophil": ['LTF','CAMP','LCN2','CXCR4','CD101',
'MMP8','MMP9','CXCR2','S100A8', 'S100A9','CSF3R'],
    "Mast":["KIT","TPSAB1","HPGD","TPSB2",'MS4A2','GATA2'],
    #"Baso": ['CD200R3', 'MCPT8', 'PRSS34'],
    "Mono": ["CD14","LYZ","CSF1R","VCAN",'FN1', 'CCR2', 'F13A1'],
    "Macro": [ 'MRC1',"ITGAM","CD86",'C1QA', 'VCAM1',"ADGRE1"],
    "DC":['ITGAX', 'CD74','IRF8','BST2','IL3RA','LILRA4'],
    "cDC1":['CLEC9A','PTTG1','HIST1H4C'],
    "cDC2":['HLA-DQA1','HLA-DQB1','CD1C','PPA1','FCER1A'],
    "pDC":['IRF8','IRF7','LILRA4','TCF4','ITM2C'],
    "T": ['CD3G','CD3D','CCL5',"CD4","FOXP3","CD8A"],
    "NK": ['KLRD1','NCAM1'],
    "Pro_B":['DNTT','VPREB1','VPREB3','IGLL1'],
    "Immature_B":['TCL1A','VPREB3','CD79B','MZB1','IGLL5'],
    "B": ['CD79A','CD19','IGHM', 'EBF1'],
    "Plasma": ['JCHAIN', 'IGLC2', 'MZB1'],  
    "Mesenchymal":["CXCL12","PDGFRA"],
    "Osteolineage":['RUNX2',"SP7","IBSP","BGLAP","COL1A1"],
    "Adipo-lineage":["APOE","CEBPA","PPARG","LPL"],
    "Fibro":["DCN","PDPN","CSPG4"], 
    "Smooth_muscle":['ACTA2','RGS5','TAGLN'],
    "Endo":['CDH5','PECAM1','VWF'],
    "Chondrocyte":["SOX9"]   
    
}


sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_2.50', standard_scale="var",dendrogram=True,show=False,save="_03-8.png")


import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 4))


sc.pl.violin(
    adata,
    ["n_genes"],
    jitter=0.4,size=0.3,
    groupby="leiden_res_2.50",
    multi_panel=True,
    show=False,
    ax=ax,
    save="_03-9.png"
)



'''
#—————————————————————— 想要剔除的细胞类型列表——————————————————
adata.obs['leiden_res_2.50'].value_counts()

adata
#AnnData object with n_obs × n_vars = 80556 × 2000

clutser_to_remove = ['55','63']
 
# 创建一个布尔掩码
mask = ~adata.obs['leiden_res_2.50'].isin(clutser_to_remove)
 
# 使用mask提取子集
adata_subset = adata[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 80474 × 2000


sc.pl.umap(adata_subset, color=[ 'leiden_res_2.50'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_test.png")

sc.pl.umap(adata, color=['pct_counts_mt'],cmap="Reds",show=False,save="_test.png")



adata_counts=adata_subset.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 80474 × 29188
print(adata_counts.X)

ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 80474 × 29188
print(adata_counts.X)




# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'cluster_anno_l2', 'cluster_anno_coarse', 'cluster_anno_l1', 'Sex', 'Age', 'Ethnicity', 'n_genes_by_counts','total_counts','pct_counts_mt']],
    var=adata_counts.var[['features']]
)



new_adata.write_h5ad('03_remove_contamin.h5ad',compression='gzip')


'''






#————————————————————————Annotation————————————————————————————

preAnno = dict()

for i in range(0,61):
    preAnno[str(i)] = "HSPC"

print(preAnno)


####开始定义细胞类型

for i in ['60','25','26']:
    preAnno[i] = "Mega_prog"
    
    
for i in ['38']:
    preAnno[i] = "Ery_prog"


for i in ['43','16','21']:
    preAnno[i] = "Erythroid"
  

for i in ['48']:
    preAnno[i] = "Mono_prog"



for i in ['40','46','47']:
    preAnno[i] = "Neu_prog"
    
    
    
for i in ['23','22','44','15','50','34','45']:
    preAnno[i] = "Neutrophil"  
     
     
for i in ['51']:
    preAnno[i] = "Baso"



for i in ['49']:
    preAnno[i] = "Monocyte"
    
    
for i in ['39']:
    preAnno[i] = "Macrophage"
    
    
for i in ['52']:
    preAnno[i] = "pDC"
    

for i in ['35','37']:
    preAnno[i] = "T"
    
    
for i in ['36']:
    preAnno[i] = "NK"



for i in ['55','53','41','58','9']:
    preAnno[i] = "Pre/pro_B"
   
   
   
for i in ['42']:
    preAnno[i] = "B"

    
for i in ['29','3','20','6','0','1','5','4','7','2','8']:
    preAnno[i] = "Plasma"
    

   
for i in ['24','32']:
    preAnno[i] = "VSMC"


for i in ['56']:   
    preAnno[i] = 'Fibro−MSC'


for i in ['12','33']:   
    preAnno[i] = 'THY1+ MSC'
    
    
for i in ['11','10']:   
    preAnno[i] = 'Adipo-MSC'
    
    
for i in ['13']:   
    preAnno[i] = 'APOD+ MSC'
    
    
for i in ['19']:   
    preAnno[i] = 'Osteo-MSC'
    
    
for i in ['17','18']:   
    preAnno[i] = 'Osteoblast'


for i in ['14']:   
    preAnno[i] = "RNAlo MSC"
    
    

for i in ['30','31','27','28']:   
    preAnno[i] = "Endothelial"


    
    
print(preAnno)


adata.obs['celltype'] = adata.obs['leiden_res_2.50'].map(preAnno)

adata.obs['celltype'] = pd.Categorical(adata.obs['celltype'],categories=["HSPC","Mega_prog","Ery_prog","Erythroid","Baso","Mono_prog","Monocyte","Macrophage","Neu_prog","Neutrophil","pDC","T","NK","Pre/pro_B","B","Plasma","VSMC",'Fibro−MSC','THY1+ MSC','Adipo-MSC','APOD+ MSC','Osteo-MSC','Osteoblast',"RNAlo MSC","Endothelial"])                   
                              
adata.uns['celltype_colors'] = ["#abddff","#b20000","#167153","#b0d45d","#cca69c","#4a1486","#9e9ac8","#df65b0","#2873B3","#2EBEBE","#F0E442","#d73027","#fc9272","#e8743c","#ffc556", "#bf812d","#009E73", "#428AC9","#129392","#FFCC4F","#F37E78","#883A96","#cca69c","#999999","#CC79A7"]
           
           
sc.pl.umap(adata, color='celltype',size=3,legend_loc='on data',legend_fontsize=4, legend_fontoutline=2,show=False, save="_04-1.png")
sc.pl.umap(adata, color='celltype',size=3,legend_loc='on data',legend_fontsize=4, legend_fontoutline=2,show=False, save="_04-1.pdf")


sc.pl.umap(adata, color='celltype',size=3,show=False, save="_04-2.png")
sc.pl.umap(adata, color='celltype',size=3,show=False, save="_04-2.pdf")



marker_genes = { 
    "HSPC": ['SPINK2', 'ADGRG1', 'CDK6'],
    "Mega_prog": ["ITGA2B"],
    "Ery_prog": ["GATA1",'TFRC'],
    "Erythro":['HBA1','HBA2' ],
    #"Hematopoietic": ["PTPRC"],
    "Eos/baso_prog": ["MS4A2"],
    #"Baso": ['CD200R3', 'MCPT8', 'PRSS34'],
    "Mye_prog": ['CSF1R','ELANE', 'MPO',"PRTN3"],
    "Mono": ["CD14","LYZ","VCAN"],
    "Macro": [ 'MRC1','C1QA'],
    "Neutrophil": ['S100A8', 'S100A9','CSF3R'],
    "pDC":['IRF8','IRF7','LILRA4'],
    "T": ['CD3G','CD3D'],
    "NK": ['KLRD1','NCAM1'],
    "Pro_B":['DNTT','VPREB1','VPREB3','IGLL1'],
    "B": ['CD79A','IGHM'],
    "Plasma": ['JCHAIN', 'SDC1', 'MZB1'], 
    "Smooth_muscle":['ACTA2','RGS5','TAGLN'],
    "Fibro":["DCN","PDPN","CSPG4"],
    "Mesenchymal":["CXCL12","COL1A1","PDGFRA"],
    "THY1+ MSC":['THY1'],
    "Adipo-lineage":["APOE","CEBPA","PPARG","LPL"],
    "APOD+ MSC":['APOD'],
    "Osteolineage":['RUNX2',"SP7","IBSP","BGLAP"],
    "Endo":['CDH5','PECAM1','VWF']
}


sc.pl.dotplot(adata, marker_genes, groupby='celltype', standard_scale="var",dendrogram=False,show=False,save="_04-3.png")
sc.pl.dotplot(adata, marker_genes, groupby='celltype', standard_scale="var",dendrogram=False,show=False,save="_04-3.pdf")

adata.write_h5ad('04_sce_anno.h5ad',compression='gzip')

#adata = sc.read_h5ad("04_sce_anno.h5ad")
#adata





#——————————————————————堆积条形图————————————————————————————

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# 提取数据
df = adata.obs[['orig.ident', 'celltype']].copy()

# 计算每个batch中每个celltype的占比
counts = df.groupby(['orig.ident', 'celltype']).size().unstack(fill_value=0)
proportions = counts.div(counts.sum(axis=1), axis=0)


colors = adata.uns['celltype_colors']

# 绘制堆积条形图
ax = proportions.plot(kind='barh', stacked=True, figsize=(8,6), color=colors,width=0.8)
# 去掉背景网格线
ax.grid(False)
plt.xlabel('Proportion')
plt.ylabel('Sample')
plt.title('Proportions within each sample')
font_properties = FontProperties(size=10) 
plt.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left',prop=font_properties)
#plt.xticks(rotation=45)
plt.tight_layout()  # 自动调整布局，以适应标签

plt.savefig('./figures/proportion_stacked_04-3.pdf')








