import numpy as np
import pandas as pd
import squidpy as sq
import scanpy as sc
import anndata as ad
import networkx as nx


# node importance
def spatial_cluster_node_importance(adata,neighbor_key,sparse,exclude=['vitality','betweeness'],key_added='leiden_spatial_importance'):
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
    sc.tl.leiden(adata_new,resolution=0.5)
    mapping = adata_new.obs['leiden'].to_dict()
    adata.obs[key_added] = adata.obs_names.map(mapping).values
    return adata


# stability
def cluster_level_spatial_stability(adata,key,method,neighbor_key=None,sparse=None):
    if method == 'centrality':
        sq.gr.centrality_scores(adata, cluster_key='cluster')
        df = adata.uns['cluster_centrality_scores']
        for key in df.columns:
            adata.obs[key] = adata.obs['cluster'].map(df[key].to_dict()).values
            adata.obs[key] = adata.obs[key].astype('float')
    elif method =='spread':
        sq.gr.ripley(adata, cluster_key='cluster', mode='L')
        df = adata.uns['cluster_ripley_L']['L_stat']
        q = 0.5
        mapping = {}
        for c, subdf in df.groupby(by='cluster'):
            t = subdf.shape[0]
            s = round(t * q)
            v = subdf.iloc[s, :]['stats']
            mapping[c] = v
        adata.obs['spread'] = adata.obs['cluster'].map(mapping).values
        adata.obs['spread'] = adata.obs['spread'].astype('float')
    elif method == 'assortativity':
        adjacency_matrix = adata.obsp[neighbor_key]
        if sparse:
            G = nx.from_scipy_sparse_matrix(adjacency_matrix)
        else:
            G = nx.from_numpy_matrix(adjacency_matrix)
        index2cluster = pd.Series(index=np.arange(adata.shape[0]), data=adata.obs[key].values).to_dict()
        all_cluster = adata.obs[key].cat.categories.tolist()
        all_index = np.arange(len(all_cluster))
        cluster2order = pd.Series(index=all_cluster, data=all_index).to_dict()
        nx.set_node_attributes(G, index2cluster, 'cluster')
        mix_mat = nx.attribute_mixing_matrix(G, attribute='cluster', mapping=cluster2order)
        mapping = {}
        for c, o in cluster2order.items():
            self_mix = mix_mat[o, o] / mix_mat[o, :].sum()
            mapping[c] = o
        adata.obs['assortativity'] = adata.obs['cluster'].map(mapping).values
        adata.obs['assortativity'] = adata.obs['assortativity'].astype('float')
    elif method == 'number_connected_components':
        adjacency_matrix = adata.obsp[neighbor_key]
        if sparse:
            G = nx.from_scipy_sparse_matrix(adjacency_matrix)
        else:
            G = nx.from_numpy_matrix(adjacency_matrix)
        adata.obs['index'] = np.arange(adata.shape[0])
        cluster2index = adata.obs.loc[:, [key, 'index']].set_index(keys='index').reset_index().groupby(by=key)[
            'index'].apply(
            lambda x: x.tolist()).to_dict()
        mapping = {}
        for c, indices in cluster2index.items():
            subgraph = G.subgraph(nodes=indices)
            n = nx.number_connected_components(subgraph)
            mapping[c] = n
        adata.obs['number_connected_components'] = adata.obs['cluster'].map(mapping).values
        adata.obs['number_connected_components'] = adata.obs['number_connected_components'].astype('float')










