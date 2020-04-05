#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:57:31 2020

@author: ligk2e
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc  
import fa2

'''
1. the most irritating part of that: dependencies confliction, please strictly obey what the tutorial ask you to do, 
in this script, just install scanpy using pip, then it will automatically have numpy, pandas and matplotlib, 
then pip install python-igraph, not igraph, no need to improt that. Then pip install fa2.

2. trajectory inference, batch effect elimination, spatial transcriptomics are more specific topics, clustering is the most
basic one, then we need to know how to label your cell, it is cubersome and needs a lot of knowledge. Plus, it depends on
what questions you are going to solve.

3. scanpy is faster and has its unique function, Seurat is better than plotting and deploying other related R packages

4. useful single cell database or repository
    4.1 GEO
    4.2 Human Cell Atlas
    4.3 single cell portal, Broad institure, nice cell labelling
    4.4 UCSC cell browser, just for visulization 


'''

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/paul15.h5ad'
sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures

adata = sc.datasets.paul15() # have install and import ipywidgets
adata.X = adata.X.astype('float64')

sc.pp.recipe_zheng17(adata)
sc.tl.pca(adata,svd_solver='arpack')
sc.pp.neighbors(adata,n_neighbors=4,n_pcs=20)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata,color='paul15_clusters',legend_loc='on data')

# optional: denoise the graph, computing the distances within a few diffusion components amounts to denoising the graph
sc.tl.diffmap(adata)
sc.pp.neighbors(adata,n_neighbors=10,use_rep='X_diffmap')
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata,color='paul15_clusters',legend_loc='on data')

import louvain
sc.tl.louvain(adata,resolution=1.0)

sc.tl.paga(adata,groups='louvain')
sc.pl.paga(adata,color=['louvain','Hba-a2','Elane','Irf8'])
sc.pl.paga(adata,color=['louvain','Itga2b','Prss34','Cma1'])
adata.obs['louvain'].cat.categories
adata.obs['louvain_anno'] = adata.obs['louvain']
adata.obs['louvain_anno'].cat.categories = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10/Ery', '11', '12', \
       '13', '14', '15', '16/Stem', '17', '18', '19/Neu', '20/Mk', '21', '22/Baso', '23', '24/Mo']
# now use annotated clusters for PAGA
sc.tl.paga(adata,groups='louvain_anno')
sc.pl.paga(adata,threshold=0.03,show=False)    ########## undirect graph, reflecting mutual relationshop 

# Recomputing the embedding using PAGA-initialization
sc.tl.draw_graph(adata,init_pos='paga')
sc.pl.draw_graph(adata,color=['louvain_anno','Itga2b','Prss34','Cma1'],legend_loc='on data')

# choose colors of clusters a bit more consistently
pl.figure(figsize=(8,2))  # (width,height)
for i in range(28):
    pl.scatter(i,1,c=sc.pl.palettes.zeileis_28[i],s=200)   #usage of a palettes, an array contain different colors
pl.show()
zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['louvain_anno_colors'])

new_colors[[16]] = zeileis_colors[[12]]  # Stem colors / green
new_colors[[10, 17, 5, 3, 15, 6, 18, 13, 7, 12]] = zeileis_colors[[5, 5, 5, 5, 11, 11, 10, 9, 21, 21]]  # Ery colors / red
new_colors[[20, 8]] = zeileis_colors[[17, 16]]  # Mk early Ery colors / yellow
new_colors[[4, 0]] = zeileis_colors[[2, 8]]  # lymph progenitors / grey
new_colors[[22]] = zeileis_colors[[18]]  # Baso / turquoise
new_colors[[19, 14, 2]] = zeileis_colors[[6, 6, 6]]  # Neu / light blue
new_colors[[24, 9, 1, 11]] = zeileis_colors[[0, 0, 0, 0]]  # Mo / dark blue
new_colors[[21, 23]] = zeileis_colors[[25, 25]]  # outliers / grey

adata.uns['louvain_anno_colors'] = new_colors

sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)  ######## embedding plot




# specifying a root cell, reconstrucing PAGA path(trajectory)

adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain_anno']  == '16/Stem')[0]
sc.tl.dpt(adata)                                           ############# trajectory plot
# gene name will be the raw in PAGA path
gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
              'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
              'Irf8', 'Csf1r', 'Ctsg']                     # monocyte
# introduct raw data for visualization of pseudotime
adata_raw = sc.datasets.paul15()
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw
sc.pl.draw_graph(adata, color=['louvain_anno', 'dpt_pseudotime'], legend_loc='on data')

# formally construct PAGA path
paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
         ('neutrophils', [16, 0, 4, 2, 14, 19]),
         ('monocytes', [16, 0, 4, 11, 1, 9, 24])]
adata.obs['distance'] = adata.obs['dpt_pseudotime']
adata.obs['clusters'] = adata.obs['louvain_anno']  # just a cosmetic change
adata.uns['clusters_colors'] = adata.uns['louvain_anno_colors']


_, axs = pl.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12}) # they don't need totalplot, so just use _
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata, path, gene_names,
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=50,
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Greys',
        groups_key='clusters',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('./write/paga_path_{}.csv'.format(descr))
pl.savefig('./figures/paga_path_paul15.pdf')
pl.show()    ###### PAGA graph




















