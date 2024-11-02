# remove double let by scrublet
# https://cellgeni.github.io/notebooks/html/new-doublets-scrublet.html
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

adata = sc.read_10x_mtx('/mnt/DATA/RNAseq/SCRNA/cellranger_output/20241020_6.31_All/Female_Ctrl_B/outs/raw_feature_bc_matrix', cache=False)

# scrub = scr.Scrublet(adata.X, expected_doublet_rate = 0.076)
scrub = scr.Scrublet(adata.X)

adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)

scrub.plot_histogram()

scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
adata.obs
adata.obs['predicted_doublets'].value_counts()
# pd.DataFrame(adata.obs).to_csv("scrublet_calls.tsv",sep = '\t',header = False)