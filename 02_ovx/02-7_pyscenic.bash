
#——————————————————————  pyscenic  ————————————————————————————————
docker pull aertslab/pyscenic:0.12.1

cd /home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/data/


docker run -it --rm \
    -v /home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/data:/data \
    aertslab/pyscenic:0.12.1 pyscenic grn \
        --num_workers 6 \
		--method grnboost2 \
        --output /data/OVX_adjacencies.tsv \
        /data/OVX_sample_adata.loom \
        /data/allTFs_mm.txt





docker run -it --rm \
    -v /home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/data:/data \
    aertslab/pyscenic:0.12.1 pyscenic ctx \
        /data/OVX_adjacencies.tsv \
		
		/data/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        /data/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        
		--annotations_fname /data/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
		
		--auc_threshold 0.05 \
		
        --expression_mtx_fname /data/OVX_sample_adata.loom \
        --mode "custom_multiprocessing" \
        --output /data/OVX_regulons.csv \
        --num_workers 6 \
		--mask_dropouts 
		#--rho_mask_dropouts
		#--no_pruning





docker run -it --rm \
    -v /home/luchun/scRNA/Bone_neu/03_age_ovx2/pyscenic/data:/data \
    aertslab/pyscenic:0.12.1 pyscenic aucell \
        /data/OVX_total_adata.loom \
        /data/OVX_total_regulons.csv \
        -o /data/OVX_total_SCENIC.loom \
        --num_workers 5
		
		
		






