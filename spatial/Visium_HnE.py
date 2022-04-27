import scanpy as sc
import anndata as ad
import squidpy as sq
import numpy as np
import pandas as pd

# load
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()

# plot
sc.pl.spatial(adata, color="cluster")


# extract image summary feature
for scale in [1.0, 2.0]:
    feature_name = f"features_summary_scale{scale}"
    sq.im.calculate_image_features(adata,img.compute(),features="summary",key_added=feature_name,n_jobs=4,scale=scale,)
adata.obsm["features"] = pd.concat([adata.obsm[f] for f in adata.obsm.keys() if "features_summary" in f], axis="columns")
adata.obsm["features"].columns = ad.utils.make_index_unique(adata.obsm["features"].columns)

# cluster on summary image feature
def cluster_features(features: pd.DataFrame, like=None) -> pd.Series:
    if like is not None:
        features = features.filter(like=like)
    adata = ad.AnnData(features)
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    return adata.obs["leiden"]
adata.obs["features_cluster"] = cluster_features(adata.obsm["features"], like="summary")
sc.pl.spatial(adata, color=["features_cluster", "cluster"])

# neighbor enrichment
sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")

# co-occurence
sq.gr.co_occurrence(adata, cluster_key="cluster")
sq.pl.co_occurrence(adata,cluster_key="cluster",clusters="Hippocampus",figsize=(8, 4))
