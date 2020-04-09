#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 17:37:22 2020

@author: ligk2e
"""

import scanpy as sc
import pandas as pd
import seaborn as sns

sc.settings.verbosity = 1
sc.logging.print_versions()
sc.settings.set_figure_params(dpi = 80,frameon=False,figsize=(3,3))

'''
Before running the command below, make sure you pip install ipywidges,
then make sure you have Jupyter Notebook installed in this environment.
'''
adata_ref = sc.datasets.pbmc3k_processed()  # reference dataset
adata = sc.datasets.pbmc68k_reduced() # dataset that you want to query labels and embeddings


# Only consider the common genes between two dataset, results in 208 genes
var_names = adata_ref.var_names.intersection(adata.var_names) # type: pandas.core.indexes.base.Index
adata_ref = adata_ref[:,var_names]
adata = adata[:,var_names]

# pca, neighbors, UMAP for adata_ref dataset, finally will get a UMAP plot
sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref,color = 'louvain')

'''
A little digress toward understanding adata.obsm['X_pca']:
    it will be an array with length of the number of observations(cells),
    each sub-element is also an array with leagth of the number of PCs, 
    each value means the projection of each cell on each PC.
'''

# map labels and embeddings from adata_ref to adata based on chosen representaion, here is adata_ref.obsm['X_pca']
sc.tl.ingest(adata,adata_ref,obs='louvain')
adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors']
sc.pl.umap(adata, color=['louvain', 'bulk_labels'], wspace=0.5)

'''
Intuitive understanding: anchor the cells in query dataset that shows strong similarity with cells in reference dataset,
then based upon that, we could pass the labels in reference dataset to query dataset. The refernce dataset should contain
enough biological variabtion to meaningfully accommate query data.

'''
# concatenate the reference dataset and query dataset 
adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
adata_concat.obs.louvain = adata_concat.obs.louvain.astype('category')
adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories,inplace = True) # fix category ordering
adata_concat.uns['louvain_colors'] = adata_ref.uns['louvain_colors'] # fix category colors
sc.pl.umap(adata_concat,color=['batch','louvain'])

# what if using bbknn? 
sc.tl.pca(adata_concat)
%%time    # magic syntax, calculate the processing time of excecuting following command, they must be shown in the same Ipython chunk
sc.external.pp.bbknn(adata_concat,batch_key = 'batch')
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat,color=['batch','louvain'])

##################################################################################################
###############         Pancreas Dataset     #################################################
######################################################################################################


# load, inspect, remove 5 minority classes
adata_all = sc.read('data/pancreas.h5ad', backup_url='https://www.dropbox.com/s/qj1jlm9w10wmt0u/pancreas.h5ad?dl=1')
adata_all.shape
counts = adata_all.obs.celltype.value_counts()
counts   # it is a pandas.core.series.Series.

# remove the last 5 minority class
minority_classes = counts.index[-5:].tolist()        # get the minority classes
adata_all = adata_all[~adata_all.obs.celltype.isin(minority_classes)]
adata_all.obs.celltype.cat.reorder_categories(counts.index[:-5].tolist(), inplace=True)

sc.pp.pca(adata_all)
sc.pp.neighbors(adata_all)
sc.tl.umap(adata_all)


'''
Above line could generate a very decent example of what batch effect is and what our goals are. by using BBKNN, we speficy 
the batch, then it will recalculate the value, eliminate the batch effect.
'''
sc.pl.umap(adata_all,color=['batch','celltype'], palette = sc.pl.palettes.vega_20_scanpy)

# using BBKNN
sc.external.pp.bbknn(adata_all,batch_key = 'batch')
sc.tl.umap(adata_all)
sc.pl.umap(adata_all,color = ['batch','celltype'])

# using ingest
adata_ref = adata_all[adata_all.obs.batch == '0']
sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref,color = 'celltype')

adatas = [adata_all[adata_all.obs.batch == i].copy() for i in ['1','2','3']]
sc.settings.verbosity = 2  # a bit more logging
for iadata, adata in enumerate(adatas):  
    print(f'... integrating batch {iadata+1}')  # formatted stirng literal
    adata.obs['celltype_orig'] = adata.obs.celltype  # save the original cell type
    sc.tl.ingest(adata, adata_ref, obs='celltype')
    
adata_concat = adata_ref.concatenate(adatas)

adata_concat.obs.celltype = adata_concat.obs.celltype.astype('category')
adata_concat.obs.celltype.cat.reorder_categories(adata_ref.obs.celltype.cat.categories, inplace=True)  # fix category ordering
adata_concat.uns['celltype_colors'] = adata_ref.uns['celltype_colors']  # fix category coloring
sc.pl.umap(adata_concat, color=['batch', 'celltype'])

# evalute ingest performance by comparing the concatenated one with cell origal labelling
adata_query = adata_concat[adata_concat.obs.batch.isin(['1', '2', '3'])]
sc.pl.umap(adata_query, color=['batch', 'celltype', 'celltype_orig'], wspace=0.4)

# consistency chack by pd.crosstab function
obs_query = adata_query.obs
conserved_categories = obs_query.celltype.cat.categories.intersection(obs_query.celltype_orig.cat.categories)  # intersected categories
obs_query_conserved = obs_query.loc[obs_query.celltype.isin(conserved_categories) & obs_query.celltype_orig.isin(conserved_categories)]  # intersect categories
obs_query_conserved.celltype.cat.remove_unused_categories(inplace=True)  # remove unused categoriyes
obs_query_conserved.celltype_orig.cat.remove_unused_categories(inplace=True)  # remove unused categoriyes
obs_query_conserved.celltype_orig.cat.reorder_categories(obs_query_conserved.celltype.cat.categories, inplace=True)  # fix category ordering
pd.crosstab(obs_query_conserved.celltype, obs_query_conserved.celltype_orig)
pd.crosstab(adata_query.obs.celltype, adata_query.obs.celltype_orig)

# easily visualize each batch in the concatenated umap plot
sc.tl.embedding_density(adata_concat, groupby='batch')
sc.pl.embedding_density(adata_concat, groupby='batch')
for batch in ['1', '2', '3']:
    sc.pl.umap(adata_concat, color='batch', groups=[batch])








































