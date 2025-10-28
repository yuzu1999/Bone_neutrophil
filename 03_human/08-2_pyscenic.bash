
#——————————————————————  pyscenic  ————————————————————————————————
docker pull aertslab/pyscenic:0.12.1

cd /home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data


docker run -it --rm \
    -v /home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data:/data \
    aertslab/pyscenic:0.12.1 pyscenic grn \
        --num_workers 15 \
		--method grnboost2 \
        --output /data/Human_adjacencies.tsv \
        /data/Human_sample_adata.loom \
        /data/allTFs_hg38.txt





docker run -it --rm \
    -v /home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data:/data \
    aertslab/pyscenic:0.12.1 pyscenic ctx \
        /data/Human_adjacencies.tsv \
		
		/data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        /data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        
		--annotations_fname /data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
		
		--auc_threshold 0.05 \
		
        --expression_mtx_fname /data/Human_total_adata.loom \
        --mode "custom_multiprocessing" \
        --output /data/Human_total_regulons.csv \
        --num_workers 10 \
		--mask_dropouts 
		#--rho_mask_dropouts
		#--no_pruning




docker run -it --rm \
    -v /home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data:/data \
    aertslab/pyscenic:0.12.1 pyscenic aucell \
        /data/Human_total_adata.loom \
        /data/Human_total_regulons.csv \
        -o /data/Human_total_SCENIC.loom \
        --num_workers 5
		
		
		



'''
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

os.chdir('/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/data/')

sc.settings.set_figure_params(dpi=50, facecolor="white")

adata = sc.read_h5ad("/home/luchun/scRNA/Bone_neu/04_cell/pyscenic/Human_sample_adata.h5ad")
adata

# ex_matrix is my expression matrix
ranking_feather = pd.read_feather("hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
ranking_feather2 = pd.read_feather("hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

ranking_feather.columns
ranking_feather2.columns

gene_saved = list(set(ranking_feather.columns) & set(ranking_feather2.columns))

overlap_values = list(set(adata.var_names) & set(ranking_feather.columns))
len(overlap_values)

ex_matrix = ex_matrix.loc[overlap_values, :]


'''








