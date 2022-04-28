import numpy as np
import pandas as pd
import squidpy as sq
import scanpy as sc
import anndata as ad
import networkx as nx


# node importance
def spatial_feature_graph_node_importance_clustering(adata,neighbor_key,sparse,exclude=['vitality','betweeness'],key_added='leiden_spatial_importance'):
    adjacency_matrix = adata.obsp[neighbor_key]
    if sparse:
        G = nx.from_scipy_sparse_matrix(adjacency_matrix)
    else:
        G = nx.from_numpy_matrix(adjacency_matrix)
    values_dict = {}
    degree_dict = dict(G.degree(weight='weight')); values_dict['degree'] = degree_dict
    cc_dict = nx.clustering(G, weight='weight'); values_dict['clustering_coeffient'] = cc_dict
    pr_dict = nx.pagerank(G, weight='weight'); values_dict['pagerank_score'] = pr_dict
    if 'vitality' not in exclude:
        vt_dict = nx.closeness_vitality(G, weight='weight'); values_dict['vitality'] = vt_dict
    if 'betweeness' not in exclude:
        btc_dict = nx.betweenness_centrality(G, weight='weight'); values_dict['betweeness_centrality'] = btc_dict
    for k,v in values_dict.items():
        adata.obs[k] = list(v.values())
    adata_new = ad.AnnData(adata.obs[list(values_dict.keys())])
    sc.pp.scale(adata_new)
    sc.pp.neighbors(adata_new)
    sc.tl.leiden(adata_new)
    mapping = adata_new.obs['leiden'].to_dict()
    adata.obs[key_added] = adata.obs_names.map(mapping).values
    return adata


# stability
def cluster_level_spatial_stability(adata,key,method):
    if method == 'centrality':
        sq.gr.centrality_scores(adata, cluster_key='cluster')
        df = adata.uns['cluster_centrality_scores']
        for key in df.columns:
            adata.obs[key] = adata.obs['cluster'].map(df[key].to_dict()).values
            adata.obs[key] = adata.obs[key].astype('float')
    elif method =='assortativity':









