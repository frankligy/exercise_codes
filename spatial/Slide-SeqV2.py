import scanpy as sc
import squidpy as sq


# load
adata = sq.datasets.slideseqv2()

# plot
sc.pl.spatial(adata, color="cluster", spot_size=30)

# neighborhood enrichment analysis
sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster", method="single", cmap="inferno", vmin=-50, vmax=100)

# ripley L, spatial homogeneity
mode = "L"
sq.gr.ripley(adata, cluster_key="cluster", mode=mode, max_dist=500)
sq.pl.ripley(adata, cluster_key="cluster", mode=mode)

# spatial variable gene, moranI
sq.gr.spatial_autocorr(adata, mode="moran")
adata.uns["moranI"].head(10)
sc.pl.spatial(adata,color=["Ttr", "Plp1", "Mbp", "Hpca", "Enpp2"],spot_size=30)


