import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import anndata as ad

# tissue image
img_dapi = sq.im.ImageContainer(img='/Users/ligk2e/Desktop/tmp/spatial/D1-2/32753-Slide3_D1-2_DAPI.tiff')
img_bf = sq.im.ImageContainer(img='/Users/ligk2e/Desktop/tmp/spatial/D1-2/32753-Slide3_D1-2_brightfield.tiff')
img_raw = sq.im.ImageContainer(img='/Users/ligk2e/Desktop/tmp/spatial/D1-2/32753-Slide3_D1-2_raw.tiff')
img_data_list = []
for img in [img_dapi,img_bf,img_raw]:
    img.compute()
    img_data_list.append(img['image'].data.squeeze())
merged_img_channels = np.stack(img_data_list,axis=2)
img = sq.im.ImageContainer(img=merged_img_channels,layer='image',dims=('y','x','channels'),scale=1.0)
img.show(layer='image',channelwise=True)


# rna
df = pd.read_csv('/Users/ligk2e/Desktop/tmp/spatial/D1-2/32753-Slide3_D1-2_results.txt',sep='\t',header=None)
df.drop(columns=4,inplace=True)
df.columns = ['y','x','count','gene']
spot_barcode = ['y'+'_'+str(y)+'_'+'x'+'_'+str(x) for y,x in zip(df['y'],df['x'])]
df['spot_barcode'] = spot_barcode
df.groupby(by=['spot_barcode','gene'])
new_df = pd.pivot_table(data=df,index='spot_barcode',columns='gene',values='count',aggfunc=np.mean,fill_value=0)
coord_df = df.loc[:,['spot_barcode','y','x']].set_index(keys='spot_barcode').loc[new_df.index,:].drop_duplicates()
adata = ad.AnnData(new_df)
adata.obsm['spatial'] = coord_df.values


# cluster rna
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.leiden(adata,resolution=0.1)
sc.pl.spatial(adata,color='leiden',spot_size=1)
adata.write('tmp.h5ad')
adata = sc.read('tmp.h5ad')
import sys
sys.path.insert(0,'.')
from preprocessing import umap_color_exceed_102
umap_color_exceed_102(adata,key='leiden')
