
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


adata=sc.read_h5ad("07_Anno_filter.h5ad")
adata
#AnnData object with n_obs × n_vars = 7048 × 2000



#—————————————————————————————— 热图————————————————————————————
adata_counts=adata.raw.to_adata().copy()
adata_counts
#AnnData object with n_obs × n_vars = 7048 × 21155


print(adata_counts.X)

adata=adata_counts


#sc.pp.scale(adata)
#adata.layers["scaled"] = adata.X

adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X

print(adata.X)




proliferating_genes = ['CKS2','MKI67','RANBP1','SPC24','ANP32B','H2AFX','CCNE2','CDC20','SMC2','PMF1','TOP2A','ESCO2','CKS1B','CENPF','CDKN3','TUBA1C','KIF22','STMN1','SMC4']

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







glycolysis_genes = ['GPI','ALDOA','GAPDH','PKM','ENO1','PGAM1','TPI1','PGK1','HK1','HK2','PFKL','PFKM','PGAM2','ENO2','PKLR']
#'HKDC1',
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







glucose_genes = ['SLC2A8','SLC5A2','SLC2A5','SLC5A11','SLC5A3','SLC2A1','SLC2A10','SLC2A3','SLC2A4','SLC2A9','SLC5A10','SLC5A9']

#'GM5134', 'SLC2A2', 'SLC45A1', 'SLC5A1', 'SLC5A4A', 'SLC5A4B'


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

maturation_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/DIFFERENTIATION.txt')]

sc.tl.score_genes(adata,maturation_genes,score_name='Differentiation score')

sc.pl.umap(adata,color=['Differentiation score'],cmap="Reds",show=False,save="_04-4.png")
sc.pl.umap(adata,color=['Differentiation score'],cmap="Reds",show=False,save="_04-4.pdf")




import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['Differentiation score'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-5.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')











Azurophil_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/AZUROPHIL.txt')]
sc.tl.score_genes(adata,Azurophil_genes,score_name='Azurophil score')



Specific_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/SPECIFIC.txt')]
sc.tl.score_genes(adata,Specific_genes,score_name='Specific score')



Tertiary_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/TERTIARY.txt')]
sc.tl.score_genes(adata,Tertiary_genes,score_name='Tertiary score')



Secretory_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/SECRETORY.txt')]
sc.tl.score_genes(adata,Secretory_genes,score_name='Secretory score')



sc.pl.umap(adata,color=['Azurophil score','Specific score','Tertiary score','Secretory score'],ncols=4,cmap="Reds",show=False,save="_04-6.png")
sc.pl.umap(adata,color=['Azurophil score','Specific score','Tertiary score','Secretory score'],ncols=4,cmap="Reds",show=False,save="_04-6.pdf")



import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

sc.pl.violin(adata, ['Azurophil score','Specific score','Tertiary score','Secretory score'],multi_panel=True,size=0, groupby="preAnno",show=False)
#jitter=0.4,size=0.3,

fig = plt.gcf()

for i, ax in enumerate(fig.axes): 
            
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.grid(False)
    
    

plt.savefig('./figures/violin_04-7.pdf', bbox_inches='tight', pad_inches=1)

plt.close('all')














Phagocytosis_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/PHAGOCYTOSIS.txt')]
sc.tl.score_genes(adata,Phagocytosis_genes,score_name='Phagocytosis score')



Chemotaxis_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/CHEMOTAXIS.txt')]
sc.tl.score_genes(adata,Chemotaxis_genes,score_name='Chemotaxis score')



Activation_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/ACTIVATION.txt')]
sc.tl.score_genes(adata,Activation_genes,score_name='Activation score')



sc.pl.umap(adata,color=['Phagocytosis score','Chemotaxis score','Activation score'],ncols=3,cmap="Reds",show=False,save="_04-8.png")

sc.pl.umap(adata,color=['Phagocytosis score','Chemotaxis score','Activation score'],ncols=3,cmap="Reds",show=False,save="_04-8.pdf")




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














Glycolysis_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/GLYCOLYSIS.txt')]
sc.tl.score_genes(adata,Glycolysis_genes,score_name='Glycolysis score')



Oxidative_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/OXIDATIVE.txt')]
sc.tl.score_genes(adata,Oxidative_genes,score_name='Oxidative phosphorylation score')



Electron_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/ELECTRON.txt')]
sc.tl.score_genes(adata,Electron_genes,score_name='Electron transport chain score')



Tricarboxylic_genes = [x.strip() for x in open('/home/luchun/scRNA/Bone_neu/04_cell/neu/TRICARBOXYLIC.txt')]
sc.tl.score_genes(adata,Tricarboxylic_genes,score_name='Tricarboxylic acid cycle score')



sc.pl.umap(adata,color=['Glycolysis score','Oxidative phosphorylation score','Electron transport chain score','Tricarboxylic acid cycle score'],ncols=4,cmap="Reds",show=False,save="_04-10.png")
sc.pl.umap(adata,color=['Glycolysis score','Oxidative phosphorylation score','Electron transport chain score','Tricarboxylic acid cycle score'],ncols=4,cmap="Reds",show=False,save="_04-10.pdf")



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


