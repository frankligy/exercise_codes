import scanpy as sc
import squidpy as sq

adata = sq.datasets.imc()
sc.pl.spatial(adata,color='cell type',spot_size=10)

# co-occurence
sq.gr.co_occurrence(adata,cluster_key='cell type')
sq.pl.co_occurrence(adata,cluster_key='cell type',clusters=['basal CK tumor cell','T cells'])

# spatial neighbor graph
sq.gr.spatial_neighbors(adata)

# neighbor enrichment
sq.gr.nhood_enrichment(adata,cluster_key='cell type')
sq.pl.nhood_enrichment(adata,cluster_key='cell type')

# interaction matrix
sq.gr.interaction_matrix(adata,cluster_key='cell type')
sq.pl.interaction_matrix(adata, cluster_key="cell type")

# centrality
sq.gr.centrality_scores(adata,cluster_key="cell type")
sq.pl.centrality_scores(adata, cluster_key="cell type")