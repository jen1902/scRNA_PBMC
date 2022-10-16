### Evaluation by Ingest
## Reference: 1. https://uppsala.instructure.com/courses/52011 2.https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html 
import scanpy as sc
import scanorama
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_versions()
sc.settings.set_figure_params(dpi=120, facecolor='white')

save_file = 'data/COVID/results/clustered_covid_2h2c_v2.h5ad'
adata = sc.read_h5ad(save_file)
print(adata)
sc.pl.umap(adata, color=['louvain_0.2'], wspace=0.5, legend_loc='on data', legend_fontsize=12)

for sample in ['healthy_1', 'healthy_2', 'covid_1', 'covid_2']:
    sc.pl.umap(adata, color='sample', groups=[sample], size=20, title=sample, legend_loc='on data', legend_fontsize=12)

sc.pl.umap(adata[adata.obs["sample"] == "healthy_1",:], color=['louvain_0.2'], size=20, title='healthy_1', legend_fontsize=15, wspace=1.0)
sc.pl.umap(adata[adata.obs["sample"] == "healthy_2",:], color=['louvain_0.2'], size=20, title='healthy_2', legend_fontsize=15, wspace=1.0)
sc.pl.umap(adata[adata.obs["sample"] == "covid_1",:], color=['louvain_0.2'], size=20, title='covid_1', legend_fontsize=15, wspace=1.0)
sc.pl.umap(adata[adata.obs["sample"] == "covid_2",:], color=['louvain_0.2'], size=20, title='covid_2', legend_fontsize=15, wspace=1.0)


## Load reference PBMC data
adata_ref = sc.datasets.pbmc3k_processed()  # this is an earlier version of the dataset from the pbmc3k tutorial

## Extract the intersection of variable genes
var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

sc.tl.pca(adata_ref, svd_solver='arpack')
sc.pp.neighbors(adata_ref, n_pcs = 25, n_neighbors = 20, random_state=0)
sc.tl.umap(adata_ref, random_state=0)

## UMAP of ref data
sc.pl.umap(adata_ref, color='louvain', size=20,  title='ref data', legend_loc='on data', legend_fontsize=8, wspace=0.5)

print(adata)
sc.tl.ingest(adata, adata_ref, obs='louvain')
print(adata)
adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix colors
sc.pl.umap(adata, color=['louvain', 'louvain_0.2'], wspace=1.0)     

