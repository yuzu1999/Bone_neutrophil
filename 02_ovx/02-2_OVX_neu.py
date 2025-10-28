
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



adata_all = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/03_age_ovx2/06_Neu_anno.h5ad")
adata_all
#AnnData object with n_obs × n_vars = 53904 × 2000




#—————————————————————— 想要的细胞类型列表——————————————————
adata_all.obs['batch'].value_counts()


clutser_to_stay = ['Sham', 'OVX']
 
# 创建一个布尔掩码
mask = adata_all.obs['batch'].isin(clutser_to_stay)
 
# 使用mask提取子集
adata_subset = adata_all[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 19324 × 2000



adata_counts=adata_subset.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 19324 × 18544
print(adata_counts.X)



ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 19324 × 18544
print(adata_counts.X)



# 使用counts的layer重新构建adata
adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['batch', 'dataset', 'n_genes_by_counts','total_counts','pct_counts_mt','celltype','preAnno']],
    var=adata_counts.var[['gene_ids']]
)
adata
#AnnData object with n_obs × n_vars = 19324 × 18544



sc.pp.filter_genes(adata, min_cells=3)
adata
#AnnData object with n_obs × n_vars = 19324 × 15981


adata.obs['batch'].value_counts()
adata.obs['dataset'].value_counts()

adata.obs['batch'] = pd.Categorical(adata.obs['batch'],categories=['Sham', 'OVX'])






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
#View of AnnData object with n_obs × n_vars = 19324 × 2000


sc.pp.scale(adata, max_value=10)



 
 
#————————————————————Dimensionality Reduction——————————————————————

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False,save="_03-3.png")





#———————Nearest neighbor graph constuction and visualization——————
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=15)
sc.tl.umap(adata)

sc.pl.umap(adata,color=['batch','phase',"preAnno"],cmap="Reds",size=5,show=False,save="_03-4.png")

                       
adata.uns['preAnno_colors'] = ["#76afda","#8264CC","#009E73", "#F0E442", "#f06152","#0072B2","#fc9272"]        

sc.pl.umap(adata, color='preAnno',size=5,show=False, save="_03-5.png")
sc.pl.umap(adata, color='preAnno',size=5,show=False, save="_03-5.pdf")



#——————————————————————Clustering——————————————————————————————————
for res in [1.5, 2,2.5, 3,3.5,4]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(adata, color=['leiden_res_1.50','leiden_res_2.00','leiden_res_2.50','leiden_res_3.00','leiden_res_3.50','leiden_res_4.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=3,ncols=3,show=False,save="_03-6.png")

sc.pl.umap(adata, color=[ 'leiden_res_4.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=3,show=False,save="_03-7.png")


sc.pl.umap(adata, color=[ 'leiden_res_4.00'],show=False,save="_03-8.png")




#——————————————————Manual cell-type annotation——————————————————
marker_genes = {
    "Hematopoietic": ["Ptprc"],
    "HSPC": ['Kit','Cd34', 'Adgrg1', 'Cdk6'],
    "Lymp_prog": ["Flt3"],
    "Ery_prog": ["Gata1",'Tfrc'],
    "Mega_prog": ["Itga2b"],
    "Mye_prog": ['Csf1r','Elane', 'Mpo',"Prtn3"],
    "Eos/baso_prog": ["Ms4a2"],
    "Neutrophil": ['Ltf','Camp','Ngp','Lcn2','Cxcr4','Cd101',
'Mmp8','Mmp9','Cxcr2','Ly6g','S100a8', 'S100a9','Csf3r'],
    "Mast":["Kit","Fcer1a","Fcer2a"],
    "Baso": ['Cd200r3', 'Mcpt8', 'Prss34'],
    "Mono": ["Cd14","Vcan","Csf1r",'Fn1', 'Ccr2', 'F13a1'],
    "Macro": ["Itgam", 'Mrc1',"Cd86",'C1qa', 'Vcam1',"Adgre1"],
    "DC":['Siglech', 'Irf8','Bst2'],
    "T": ['Cd3g','Cd3d','Ccl5', 'Ms4a4b',"Cd4","Cd8a"],
    "NK": ['Klrd1','Ncam1'],
    "Pre/pro_B":['Vpreb1'],
    "B": ['Cd79a','Cd19','Ighm', 'Ebf1'],
    "Plasma": ['Jchain', 'Iglc2', 'Mzb1'],
    "Erythro":['Hbb-bt','Hbb-bs' , 'Hba-a1'],
    "Mesenchymal":["Cxcl12","Pdgfra","Vim","Lepr",'Spp1',"Mcam","Nt5e","Thy1","Eng"],
    "Osteolineage":['Runx2',"Ibsp","Bglap","Col1a1"],
    "Adipo-lineage":["Apoe","Cebpa","Pparg","Lpl"],
    "Fibro":["Dcn","Pdpn"], 
    "Smooth_muscle":['Myh11','Acta2','Tagln'],
    "Endo":["Pecam1",'Vwf'],
    "Chondrocyte":["Sox9"]   
    
}

sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_4.00', standard_scale="var",dendrogram=True,show=False,save="_03-9.png")



marker_genes = {
    "Prog": ['Mpo', 'Elane', 'Prtn3', 'Ctsg', 'Ms4a3'],
    "Proliferating": ['Camp', 'Ptma', 'Stmn1','Mki67', 'Ppia', 'Hmgb2', 'Tuba1b', 'Pclaf', 'Anp32b', 'Top2a', 'Smc4', 'Tagln2', 'Birc5', 'Ncapd2'],
    "Immature": ['Ltf', 'Lcn2', 'Ngp', 'Cd177', 'S100a8'],
    "Mature": ['S100a11', 'Ccl6', 'Retnlg', 'S100a6', 'Sell', 'Slpi', 'Clec4d', 'Lcp1', 'Srgn', 'Fpr1', 'Tyrobp', 'Cd33', 'Anxa2', 'Selplg','Fos', 'C5ar1', 'Dusp1', 'Samhd1', 'Csf3r'],
    "Senescent": ['Il1b', 'Zfp36', 'Pim1', 'Cxcl2', 'Dusp1', 'Nfkbia', 'Egr1', 'Btg2', 'Marcks', 'Mcl1', 'Cebpb', 'Nr4a1', 'Ptgs2']
  
}

sc.pl.dotplot(adata, marker_genes, groupby='leiden_res_4.00', standard_scale="var",dendrogram=True,show=False,save="_03-10.png")




import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(18, 4))


sc.pl.violin(
    adata,
    ['n_genes_by_counts'],
    jitter=0.4,size=0.3,
    groupby="leiden_res_4.00",
    multi_panel=True,
    show=False,
    ax=ax,
    save="_03-11.png"
)




#—————————————————————— 想要剔除的细胞类型列表——————————————————
adata.obs['leiden_res_4.00'].value_counts()


clutser_to_remove = ['49','51']
 
# 创建一个布尔掩码
mask = ~adata.obs['leiden_res_4.00'].isin(clutser_to_remove)
 
# 使用mask提取子集
adata_subset = adata[mask, :]
adata_subset
#View of AnnData object with n_obs × n_vars = 18784 × 2000


sc.pl.umap(adata_subset, color=[ 'leiden_res_4.00'],legend_loc='on data',legend_fontsize=8, legend_fontoutline=2,show=False,save="_test.png")


adata_counts=adata_subset.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 18784 × 15981

print(adata_counts.X)



ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 18784 × 15981
print(adata_counts.X)



# 使用counts的layer重新构建adata
new_adata = sc.AnnData(
    X=adata_counts.X,
    obs=adata_counts.obs[['batch', 'dataset', 'n_genes_by_counts','total_counts','pct_counts_mt','celltype','preAnno']],
    var=adata_counts.var[['gene_ids']]
)



adata.write_h5ad('03_Neu_raw.h5ad',compression='gzip')

new_adata.write_h5ad('04_Neu_remove_contamin.h5ad',compression='gzip')


#adata = sc.read_h5ad("03_Neu_raw.h5ad")
#adata






