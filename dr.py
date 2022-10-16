### Dimension Reduction
## Reference: 1. https://uppsala.instructure.com/courses/52011 2.https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html 
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


adata = sc.read_h5ad('data/COVID/filtered_covid_2h2c_v2.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


adata.raw = adata 

### Dimension reduction - PCA
# Compute highly variable genes  (features/genes are important in our dataset to distinguish cell types)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
print("Highly variable genes: %d"%sum(adata.var.highly_variable))

#plot variable genes
sc.pl.highly_variable_genes(adata)

# subset for variable genes in the dataset
adata = adata[:, adata.var['highly_variable']]

# regress out unwanted variables
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# scale data, clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

### PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='sample')
print(adata)

#plot variance ratio
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 40)

### Visualization tSNE
sc.tl.tsne(adata, n_pcs = 25, random_state=0)
sc.pl.tsne(adata, color='sample')

### Visualization UMAP
sc.pp.neighbors(adata, n_pcs = 25, n_neighbors = 20, random_state=0)
sc.pl.umap(adata, color='sample', edges=True)

adata.write_h5ad('data/COVID/results/dr_covid_2h2c_v2.h5ad')
