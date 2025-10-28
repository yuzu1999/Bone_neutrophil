
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

os.chdir('/home/luchun/scRNA/Bone_neu/02_ovx/')

sc.settings.set_figure_params(dpi=50, facecolor="white")


adata=sc.read_h5ad("05_neu_anno.h5ad")
adata
#AnnData object with n_obs × n_vars = 19324 × 2000



#—————————————————————————————— 热图————————————————————————————
adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 19324 × 15981

print(adata_counts.X)

adata=adata_counts


#sc.pp.scale(adata)
#adata.layers["scaled"] = adata.X

adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X

print(adata.X)




proliferating_genes = ['Cks2','Mki67','Ranbp1','Spc24','Anp32b','H2afx','Ccne2','Cdc20','Smc2','Pmf1','Top2a','Esco2','Cks1b','Cenpf','Cdkn3','Tuba1c','Kif22','Stmn1','Smc4']

import matplotlib.pyplot as plt

sc.pl.heatmap(adata,proliferating_genes,groupby="preAnno",  layer="scaled",show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 9), show=False)



fig = plt.gcf()

for ax in fig.axes:
    #ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    

plt.savefig('./figures/heatmap_04-1.png', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_04-1.pdf', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_04-1.svg', bbox_inches='tight', pad_inches=1)


plt.close('all')







glycolysis_genes = ['Gpi1','Aldoa','Gapdh','Pkm','Eno1','Pgam1','Tpi1','Pgk1','Hk1','Hk2','Pfkl','Pfkm','Pklr']
#'Hkdc1','Pgam2','Eno2',
import matplotlib.pyplot as plt

sc.pl.heatmap(adata,glycolysis_genes,groupby="preAnno",  layer="scaled",show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 7), show=False)



fig = plt.gcf()

for ax in fig.axes:
    #ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    

plt.savefig('./figures/heatmap_04-2.png', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_04-2.pdf', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_04-2.svg', bbox_inches='tight', pad_inches=1)

plt.close('all')







glucose_genes = ['Slc2a8','Slc5a11','Slc5a3','Slc2a1','Slc2a10','Slc2a3','Slc2a4','Slc2a9','Gm5134']

#'Slc2a2', 'Slc45a1', 'Slc5a1', 'Slc5a10', 'Slc5a2', 'Slc5a4b', 'Slc5a9''Slc2a5',

import matplotlib.pyplot as plt

sc.pl.heatmap(adata,glucose_genes,groupby="preAnno",  layer="scaled",show_gene_labels=True,vmin=-2,vmax=2,cmap="RdBu_r",dendrogram=False,swap_axes=True,figsize=(8, 7), show=False)



fig = plt.gcf()

for ax in fig.axes:
    #ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)

    

plt.savefig('./figures/heatmap_04-3.png', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_04-3.pdf', bbox_inches='tight', pad_inches=1)
plt.savefig('./figures/heatmap_04-3.svg', bbox_inches='tight', pad_inches=1)


plt.close('all')







#—————————————————————————————score————————————————————————————————

maturation_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/maturation_list.txt')]

sc.tl.score_genes(adata,maturation_genes,score_name='Maturation score')

sc.pl.umap(adata,color=['Maturation score'],cmap="Reds",size=5,show=False,save="_04-4.png")
sc.pl.umap(adata,color=['Maturation score'],cmap="Reds",size=5,show=False,save="_04-4.pdf")




import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['Maturation score'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-5.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')











Azurophil_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Azurophil_lis.txt')]
sc.tl.score_genes(adata,Azurophil_genes,score_name='Azurophil score')



Specific_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Specific_lis.txt')]
sc.tl.score_genes(adata,Specific_genes,score_name='Specific score')



Gelatinase_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Gelatinase_lis.txt')]
sc.tl.score_genes(adata,Gelatinase_genes,score_name='Gelatinase score')



Secretory_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Secretory_lis.txt')]
sc.tl.score_genes(adata,Secretory_genes,score_name='Secretory score')



sc.pl.umap(adata,color=['Azurophil score','Specific score','Gelatinase score','Secretory score'],ncols=4,cmap="Reds",size=5,show=False,save="_04-6.png")
sc.pl.umap(adata,color=['Azurophil score','Specific score','Gelatinase score','Secretory score'],ncols=4,cmap="Reds",size=5,show=False,save="_04-6.pdf")



import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['Azurophil score','Specific score','Gelatinase score','Secretory score'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-7.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')














Phagocytosis_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Phagocytosis_lis.txt')]
sc.tl.score_genes(adata,Phagocytosis_genes,score_name='Phagocytosis score')



Chemotaxis_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Chemotaxis_lis.txt')]
sc.tl.score_genes(adata,Chemotaxis_genes,score_name='Chemotaxis score')



Activation_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Activation_lis.txt')]
sc.tl.score_genes(adata,Activation_genes,score_name='Activation score')



sc.pl.umap(adata,color=['Phagocytosis score','Chemotaxis score','Activation score'],ncols=3,cmap="Reds",size=5,show=False,save="_04-8.png")

sc.pl.umap(adata,color=['Phagocytosis score','Chemotaxis score','Activation score'],ncols=3,cmap="Reds",size=5,show=False,save="_04-8.pdf")




import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['Phagocytosis score','Chemotaxis score','Activation score'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-9.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')














Glycolysis_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Glycolysis_lis.txt')]
sc.tl.score_genes(adata,Glycolysis_genes,score_name='Glycolysis score')



Oxidative_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Oxidative.txt')]
sc.tl.score_genes(adata,Oxidative_genes,score_name='Oxidative phosphorylation score')



Electron_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Electron.txt')]
sc.tl.score_genes(adata,Electron_genes,score_name='Electron transport chain score')



Tricarboxylic_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/Tricarboxylic.txt')]
sc.tl.score_genes(adata,Tricarboxylic_genes,score_name='Tricarboxylic acid cycle score')



sc.pl.umap(adata,color=['Glycolysis score','Oxidative phosphorylation score','Electron transport chain score','Tricarboxylic acid cycle score'],ncols=4,cmap="Reds",size=5,show=False,save="_04-10.png")
sc.pl.umap(adata,color=['Glycolysis score','Oxidative phosphorylation score','Electron transport chain score','Tricarboxylic acid cycle score'],ncols=4,cmap="Reds",size=5,show=False,save="_04-10.pdf")



import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['Glycolysis score','Oxidative phosphorylation score','Electron transport chain score','Tricarboxylic acid cycle score'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-11.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')


