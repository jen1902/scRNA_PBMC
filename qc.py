###Merge Data & Quality Control
## Reference: 1. https://uppsala.instructure.com/courses/52011 2.https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html 
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


### Load data
data_healthy1 = sc.read_10x_mtx(
    'data/COVID/PBMC_healthy [22CBD6-00]/',  # the directory with the `.mtx` file
    var_names='gene_symbols',              
    cache=True)                            
data_healthy1.var_names_make_unique()  
print(data_healthy1)

data_healthy2 = sc.read_10x_mtx(
    'data/COVID/PBMC_healthy [52BA23-00]/',  
    var_names='gene_symbols',              
    cache=True)                             
data_healthy2.var_names_make_unique()  
print(data_healthy2)

data_cov1 = sc.read_10x_mtx(
    'data/COVID/PBMC_hospitalized [2DB5F9-00]/', 
    var_names='gene_symbols',  
    cache=True)                             
data_cov1.var_names_make_unique()  
print(data_cov1)

data_cov2 = sc.read_10x_mtx(
    'data/COVID/PBMC_hospitalized [460F8D-00]/', 
    var_names='gene_symbols',               
    cache=True)                              
data_cov2.var_names_make_unique()  
print(data_cov2)



### Create one merged object
# add metadata
data_healthy1.obs['type']="Healthy"
data_healthy1.obs['sample']="healthy_1"
data_healthy2.obs['type']="Healthy"
data_healthy2.obs['sample']="healthy_2"
data_cov1.obs['type']="Covid"
data_cov1.obs['sample']="covid_1"
data_cov2.obs['type']="Covid"
data_cov2.obs['sample']="covid_2"


# merge into one object
adata = data_healthy1.concatenate(data_healthy2, data_cov1, data_cov2)
print(adata.obs['sample'].value_counts())
print(adata)

# QC metrics: mitochondiral gene
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


### QC control
## Visualization before filtering
# Histogram plots
plt.hist(adata.obs['n_genes_by_counts'], bins=100)
plt.title('Genes by counts before filtering')       #need to remove cells with deted genes less than500, by sc.pp.filter_cells(adata, min_genes=500)
plt.xlabel('N genes')
plt.ylabel('N cells')
plt.show()

plt.hist(adata.obs['total_counts'], bins=100)
plt.title('Total counts before filtering')      #need to delete genes
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.show()

# Violin plot
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], groupby='sample', rotation=45, title='Violin plot before filtering by sample')
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, title='Violin plot before filtering')

# Gene fraction visualization
sc.pl.highest_expr_genes(adata, n_top=20)


## Basic filtering
sc.pp.filter_cells(adata, min_genes=500)    #cells with fewer than 500 genes are removed
sc.pp.filter_genes(adata, min_cells=10)     #genes expressed by fewer than 10 cells are removed

## Mito filtering
adata = adata[adata.obs['pct_counts_mt'] < 10, :]
print("Remaining cells %d"%adata.n_obs)


## Visualization after filtering
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# Histogram plots
plt.hist(adata.obs['n_genes_by_counts'], bins=100)
plt.title('Genes by counts after filtering')
plt.xlabel('N genes')
plt.ylabel('N cells')
plt.show()

plt.hist(adata.obs['total_counts'], bins=100)
plt.title('Total counts after filtering')
plt.xlabel('Total counts')
plt.ylabel('N cells')
plt.show()

# Violin plot
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], groupby='sample', rotation=45, title='Violin plot after filtering by sample')
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, title='Violin plot after filtering')


## Gene filtering
malat1 = adata.var_names.str.startswith('MALAT1')
mito_genes = adata.var_names.str.startswith('MT-')      # we need to redefine the mito_genes since they were first calculated on the full object before removing low expressed genes.
hb_genes = adata.var_names.str.contains('^HB[^(P)]')

remove = np.add(mito_genes, malat1)
remove = np.add(remove, hb_genes)
keep = np.invert(remove)

adata = adata[:,keep]
print(adata)
print(adata.n_obs, adata.n_vars)
sc.pl.highest_expr_genes(adata, n_top=20)


adata = adata[adata.obs.pct_counts_mt < 10, :]   #keep cells with % mt genes<10



## Total count normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, title='after normalization')
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], groupby='sample', rotation=45, title='after normalization')

# save normalized counts in raw slot.
adata.raw = adata

## Doublet prediction
scrub = scr.Scrublet(adata.raw.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
scrub.plot_histogram()

sum(adata.obs['predicted_doublets'])

adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
sc.pl.violin(adata, 'n_genes_by_counts',
             jitter=0.4, groupby = 'doublet_info', rotation=45)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.pl.umap(adata, color=['doublet_scores','doublet_info','sample'])


adata = adata.raw.to_adata()

adata = adata[adata.obs['doublet_info'] == 'False',:]
print(adata.shape)

adata.write('data/COVID/filtered_covid_2h2c_v2.h5ad')