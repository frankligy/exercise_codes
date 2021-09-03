import scanpy as sc
import umap
import numpy as np
from sklearn.decomposition import PCA

tnc2_adata = sc.read('/Users/ligk2e/Downloads/TNC2/adata_after_scanpy_recipe_rna_0.5_1_2_3_umap_True.h5ad')
tnc1_adata = sc.read('/Users/ligk2e/Downloads/TNC/adata_after_scanpy_recipe_rna_0.5_1_2_3_umap_True.h5ad')
common_genes = list(set(tnc2_adata.var_names).intersection(set(tnc1_adata.var_names)))
tnc2_adata = tnc2_adata[:,common_genes]
tnc1_adata = tnc1_adata[:,common_genes]
X2 = tnc2_adata.X.toarray()
X1 = tnc1_adata.X.toarray()

# train on TNC1
tnc1_pca_model = PCA(n_components=50).fit(X1)
tnc1_pca_score = tnc1_pca_model.transform(X1)
tnc1_umap_model = umap.UMAP().fit(tnc1_pca_score)
tnc1_umap_coord = tnc1_umap_model.embedding_

# transform on TNC2
tnc2_pca_score = tnc1_pca_model.transform(X2)
tnc2_umap_coord = tnc1_umap_model.transform(tnc2_pca_score)

# visualize
tnc1_adata.obsm['X_umap'] = tnc1_umap_coord
tnc2_adata.obsm['X_umap'] = tnc2_umap_coord
sc.pl.umap(tnc1_adata,color='sctri_rna_leiden_0.5')
sc.pl.umap(tnc2_adata,color='sctri_rna_leiden_0.5')


















