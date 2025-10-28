

import sctour as sct
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

os.chdir('/home/luchun/scRNA/Bone_neu/01_age/scTour/')

sc.settings.set_figure_params(dpi=50, facecolor="white")




adata=sc.read_h5ad("/home/luchun/scRNA/Bone_neu/01_age/05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 19324 × 2000


adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 19324 × 15981
print(adata_counts.X)




ov.utils.retrieve_layers(adata_counts,layers='counts')
adata_counts
#AnnData object with n_obs × n_vars = 19324 × 15981
print(adata_counts.X)


dense_array = adata_counts.X.toarray()
float_array = dense_array.astype(np.float32)
print(float_array.dtype)

adata_counts.X = float_array


print(adata_counts.X)



adata=adata_counts

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000, subset=True)





#——————————————————————Train the scTour model————————————————————————
tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
tnode.train()



#——————————————————————Infer cellular dynamics——————————————————————
adata.obs['ptime'] = tnode.get_time()

sc.pl.umap(adata, color='ptime',  show=False,save="_01-1.png")


mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs


adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])


#adata = adata[np.argsort(adata.obs['ptime'].values), :]


sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='preAnno',frameon=False, show=False, legend_loc='none', size=5, alpha=1,save="./figures/stream_01-2.png")

sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='preAnno',frameon=False, show=False, legend_loc='none', size=5, alpha=1,save="./figures/stream_01-2.pdf")




adata.write_h5ad('Age_scTour.h5ad',compression='gzip')

#adata=sc.read_h5ad("Age_scTour.h5ad")
#adata









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
 


sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='preAnno3',frameon=False, show=False, legend_loc='none', size=5, alpha=1,save="./figures/stream_01-3.png")

sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='preAnno3',frameon=False, show=False, legend_loc='none', size=5, alpha=1,save="./figures/stream_01-3.pdf")






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

sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='preAnno4',frameon=False, show=False, legend_loc='none', size=5, alpha=1,save="./figures/stream_01-4.png")

sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', color='preAnno4',frameon=False, show=False, legend_loc='none', size=5, alpha=1,save="./figures/stream_01-4.pdf")