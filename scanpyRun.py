#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:35:51 2020

@author: ligk2e
"""

'''
conda create -n scanpyEnv python=3.7
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda scanpy
'''


'''
Download pbmc3k_filtered_gene_bc_matrices.tar.gz 
from http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
Using curl, wget or directly downloading and move the tar.gz to /scanpy/data
tar -xzf pbmc_filtered_gene_bc_matrices.tar.gz
'''

import numpy as np
import pandas as pd
import scanpy as sc

 
sc.settings.verbosity =3  # errors(0), warnings(1), info(2),hints(3)
sc.logging.print_versions()  # print all dependency and their version
sc.settings.set_figure_params(dpi=80)

import os
os.chdir('/Users/ligk2e/Desktop/scanpy/')

results_file = './write/pbmc3k.h5ad'
help(sc.read_10x_mtx)  # can always access document by typing "help"
adata = sc.read_10x_mtx('./data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols',cache=True)
adata.var_names_make_unique()

# should dissect and understand the structure of adata object, X, obs, vars, uns components, refer to gmail bookmark, scanpy tutorial

sc.pl.highest_expr_genes(adata,n_top=20,)

sc.pp.filter_cells(adata,min_genes=200)
sc.pp.filter_genes(adata,min_cells=3)   # sc.pp.filter_genes: here '.' is indicative of sc is a directory, pp is a subdirectory, they are called module

mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:,mito_genes].X, axis =1).A1 / np.sum(adata.X, axis =1).A1  # here I guess adata is inheriting 
# some property of its parental class, so adata[:,mito_genes] could work, otherwise, it doesn't possess this kind of method
# np.sum and matrixName.sum are the same meaning, A1 is a np method, returning a flattened darray.
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

sc.pl.violin(adata,['n_genes','n_counts','percent_mito'],jitter=0.4, multi_panel=True)
sc.pl.scatter(adata,x='n_counts',y='percent_mito')
sc.pl.scatter(adata,x='n_counts',y='n_genes')

adata  # easily to check the X, obs, var components ,they are all dataframe

adata = adata[adata.obs.n_genes < 2500, :]
adata = adata[adata.obs.percent_mito < 0.05, :]

sc.pp.normalize_total(adata,target_sum = 1e4) # normalize counts per cell, per 10,000 reads

sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=3,min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata = adata[:,adata.var.highly_variable]

adata = adata.copy()  # different from official tutorials
sc.pp.regress_out(adata, ['n_counts', 'percent_mito']) # regress out effects of 'n_counts' and 'percent_mito'
sc.pp.scale(adata,max_value=10) # scale to unit variance 

'''                      
After preprocessing, we are goona to PCA, computing neighbor, clustering, marker gene finding
'''

sc.tl.pca(adata,svd_solver='arpack')
sc.pl.pca(adata,color='CST3')
sc.pl.pca_variance_ratio(adata,log=True)  # img is stored somewhere in the object, so this command don't return anything, results NoneType
adata.write(results_file)   # save the whole file

# adata.uns is a collection.orderedDict, {'pca':{'variance':[],'variance_ratio':[]}}

sc.pp.neighbors(adata,n_neighbors=10,n_pcs=40)  # compute neighbor as part of generating UMAP?
sc.tl.umap(adata)  #UMAP could be indicative of global connectivity compared to tSNE
sc.pl.umap(adata,color=['CST3','NKG7','PPBP'])
sc.pl.umap(adata,color=['CST3','NKG7','PPBP'],use_raw=False)
import leidenalg
sc.tl.leiden(adata)
sc.pl.umap(adata,color=['leiden','CST3','NKG7'])
adata.write(results_file)

sc.tl.rank_genes_groups(adata,'leiden',method='t-test') # other methods like Wilconxon rank-sum, MAST, limma, DESeq2, diffxpy
sc.pl.rank_genes_groups(adata,n_genes=25,sharey=False)
sc.settings.verbosity = 2
sc.tl.rank_genes_groups(adata,'leiden',method='wilcoxon')
sc.pl.rank_genes_groups(adata,n_genes=25, sharey=False)
adata.write(results_file)

sc.tl.rank_genes_groups(adata,'leiden',method='logreg') # logistic regression to do comparison, it is multi-variate, all before are uni-variate
sc.pl.rank_genes_groups(adata,n_genes=25, sharey = False)

marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']

adata = sc.read(results_file)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

sc.tl.rank_genes_groups(adata,'leiden',groups=['0'],reference='1',method='wilcoxon')
sc.pl.rank_genes_groups(adata,groups=['0'],n_genes=20)
sc.pl.rank_genes_groups_violin(adata,groups='0',n_genes=8)
adata = sc.read(results_file)
sc.pl.rank_genes_groups_violin(adata,groups='0',n_genes=8)
sc.pl.violin(adata,['CST3','NKG7','PPBP'],groupby='leiden')

# mark the cell types
new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes','add_one'] # here we got 9 clusters, I don't know why
adata.rename_categories('leiden', new_cluster_names)

sc.pl.umap(adata,color='leiden',legend_loc='on data',title='',frameon = False, save='.pdf')
ax = sc.pl.dotplot(adata,marker_genes,groupby='leiden')
ax = sc.pl.stacked_violin(adata,marker_genes,groupby='leiden',rotation=90)
adata
adata.write(results_file,compression='gzip')

# adata.obs[['n_counts','leiden']].to_csv(path)


























