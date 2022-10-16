### Clustering
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from matplotlib.pyplot import rc_context

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=120)

save_file = 'data/COVID/results/dr_covid_2h2c_v2.h5ad'
adata = sc.read_h5ad(save_file)
print(adata)



### Louvain clustering
sc.tl.louvain(adata, resolution = 0.2, key_added = "louvain_0.2", random_state=0)
sc.tl.louvain(adata, resolution = 0.35, key_added = "louvain_0.35", random_state=0)
sc.tl.louvain(adata, resolution = 0.5, key_added = "louvain_0.5", random_state=0)

sc.pl.umap(adata, color=['louvain_0.2', 'louvain_0.35', 'louvain_0.5'])

### rename the cluster
new_cluster_names = ['CD4 T cell', 'FCGR3A+ monocyte', 'CD14+ monocyte', 'CD8 T cell', 'NK cell', 'B cell', 'Convntional DCs', 'Plasmacytoid DCs']
adata.rename_categories('louvain_0.2', new_cluster_names)

### plot results
for sample in ['healthy_1', 'healthy_2', 'covid_1', 'covid_2']:
    sc.pl.umap(adata, color='sample', groups=[sample], size=20, title=sample, legend_loc='on data', legend_fontsize=12, wspace=2.0)

#sc.tl.dendrogram(adata, groupby = "louvain_0.2")
#sc.pl.dendrogram(adata, groupby = "louvain_0.2")


#fig, axs = plt.subplots(2, 2, constrained_layout=True)    #加入wspace axis
sc.pl.umap(adata[adata.obs["sample"] == "healthy_1",:], color=['louvain_0.2'], size=20, title='healthy_1', legend_fontsize=8, wspace=1.0, legend_loc='on data')
sc.pl.umap(adata[adata.obs["sample"] == "healthy_2",:], color=['louvain_0.2'], size=20, title='healthy_2', legend_fontsize=8, wspace=1.0, legend_loc='on data')
sc.pl.umap(adata[adata.obs["sample"] == "covid_1",:], color=['louvain_0.2'], size=20, title='covid_1', legend_fontsize=8, wspace=1.0, legend_loc='on data')
sc.pl.umap(adata[adata.obs["sample"] == "covid_2",:], color=['louvain_0.2'], size=20, title='covid_2', legend_fontsize=8, wspace=1.0, legend_loc='on data')
#plt.show()

genes = ['CD3E', 'CD3D', 'LYZ','IL7R', 'CD14', 'MS4A7', 'FCGR3A', 'CD8A', 'NKG7', 'GNLY', 'CD79A', 'MS4A1', 'FCER1A', 'CST3', 'IL3RA', 'ITM2C', 'PPBP']


sc.pl.dotplot(adata, genes, groupby='louvain_0.2', figsize=(10, 8))

sc.tl.rank_genes_groups(adata, 'louvain_0.2', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="wilcoxon")

sc.tl.rank_genes_groups(adata, 'louvain_0.2', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


### Save clustered data
adata.write_h5ad('data/COVID/results/clustered_covid_2h2c_v2.h5ad')
print(adata)



#### Additional clustering methods
### K means clustering
# extract pca coordinates
X_pca = adata.obsm['X_pca'] ###need to change

# kmeans with k=8
kmeans = KMeans(n_clusters=8, random_state=0).fit(X_pca)
adata.obs['kmeans8'] = kmeans.labels_.astype(str)

# kmeans with k=10
kmeans = KMeans(n_clusters=10, random_state=0).fit(X_pca)
adata.obs['kmeans10'] = kmeans.labels_.astype(str)

sc.pl.umap(adata, color=['kmeans8', 'kmeans10'])

sc.tl.dendrogram(adata, groupby = "kmeans10")
sc.pl.dendrogram(adata, groupby = "kmeans10")


### Leiden clustering
sc.tl.leiden(adata, key_added = "leiden_1.0") # default resolution in 1.0
sc.tl.leiden(adata, resolution = 0.3, key_added = "leiden_0.3")
sc.tl.leiden(adata, resolution = 0.2, key_added = "leiden_0.2")
sc.tl.leiden(adata, resolution = 0.25, key_added = "leiden_0.25")

sc.pl.umap(adata, color=['leiden_0.2', 'leiden_0.3', 'leiden_1.0','leiden_0.25'])