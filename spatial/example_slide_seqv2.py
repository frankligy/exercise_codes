import numpy as np
import pandas as pd
import squidpy as sq
import scanpy as sc
import anndata as ad
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

adata = sq.datasets.slideseqv2()

# derive spatial cluster based on spatial coord
sq.gr.spatial_neighbors(adata,coord_type='generic',n_neighs=6)
sc.tl.leiden(adata,key_added='leiden_spatial_coord_knn',neighbors_key='spatial_neighbors',resolution=0.5)
sc.pl.spatial(adata,color='leiden_spatial_coord_knn',spot_size=20)
plt.savefig('spatial_coord_knn.pdf',bbox_inches='tight');plt.close()

# rna cluster
sc.pl.spatial(adata,color='cluster',spot_size=20)
plt.savefig('rna_cluster.pdf',bbox_inches='tight');plt.close()

# derive spatial cluster based on graph node importance
spatial_cluster_node_importance(adata,neighbor_key='spatial_distances',sparse=True)
sc.pl.spatial(adata,color='leiden_spatial_importance',spot_size=20)
plt.savefig('spatial_importance_knn.pdf',bbox_inches='tight');plt.close()

# derive spatial cluster from tissue image (squidpy existing pipeline)

# cluster-level spatial stability
cluster_level_spatial_stability(adata,'cluster',method='centrality')
sc.pl.spatial(adata,color='degree_centrality',spot_size=20)
plt.savefig('degree_centrality.pdf',bbox_inches='tight');plt.close()

cluster_level_spatial_stability(adata,'cluster',method='spread')
sc.pl.spatial(adata,color='spread',spot_size=20)
plt.savefig('spread.pdf',bbox_inches='tight');plt.close()

cluster_level_spatial_stability(adata,'cluster',method='assortativity',neighbor_key='spatial_distances',sparse=True)
sc.pl.spatial(adata,color='assortativity',spot_size=20)
plt.savefig('assortativity.pdf',bbox_inches='tight');plt.close()

cluster_level_spatial_stability(adata,'cluster',method='number_connected_components',neighbor_key='spatial_distances',sparse=True)
sc.pl.spatial(adata,color='number_connected_components',spot_size=20)
plt.savefig('number_connected_components.pdf',bbox_inches='tight');plt.close()