#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import os
import sys
import scanpy as sc
import anndata as ad
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
import harmonypy as hm


'''load data'''
# cite adt 42794 × 52
cite_adt_adata = small_txt_to_adata(int_file='../CITE-Seq_ND251/ADT-TotalVI/totalvi_denoised_adt_values_Xuan_citeseq2.txt',gene_is_index=False)

# cite rna 41520 × 32738
large_txt_to_mtx(int_file='/data/salomonis2/Grimes/RNA/scRNA-Seq/10x-Genomics/LGCHMC35_Xuan/211202_Grimes_GSL-PY-2514_TEA-Seq/cellRanger-ADT/SoupX/CITE-Seq/MergedFiles-clean2.txt',out_folder='./data_cite_rna',gene_is_index=True,type_convert_to='float32')
cite_rna_adata = mtx_to_adata(int_folder='./data_cite_rna',gene_is_index=True,feature='genes')


# tea adt  36771 × 52
tea_adt_adata = small_txt_to_adata(int_file='../TEA-Seq_ND251/ADT-TotalVI-corrected-log2/TotalVI-clean-log2.txt',gene_is_index=True)
tea_adt_adata.var_names = [item.replace('__','--')for item in tea_adt_adata.var_names]
tea_adt_adata.var.rename(index={'Fc_RI--_ADT':'Fc_RI_--ADT'},inplace=True)


# tea rna 36771 × 32738
large_txt_to_mtx(int_file='../TEA-Seq_ND251/RNA/RNA_SoupX-corrected-0.5/MergedFiles-QC.txt',out_folder='./data_tea_rna',gene_is_index=True,type_convert_to='float32')
tea_rna_adata = mtx_to_adata('./data_tea_rna',gene_is_index=True,feature='genes')

# tea atac 51919 × 533227
tea_atac_adata = mtx_to_adata('../TEA-Seq_ND251/ATAC',gene_is_index=False,feature='features')

'''process data'''
common_cite = list(set(cite_adt_adata.obs_names).intersection(set(cite_rna_adata.obs_names)))  # 41520
common_tea = list(set(tea_adt_adata.obs_names).intersection(set(tea_rna_adata.obs_names)).intersection(set(tea_atac_adata.obs_names))) # 36771
cite_adt_adata = cite_adt_adata[common_cite,:]
cite_rna_adata = cite_rna_adata[common_cite,:]
tea_adt_adata = tea_adt_adata[common_tea,cite_adt_adata.var_names]
tea_rna_adata = tea_rna_adata[common_tea,cite_rna_adata.var_names]
tea_atac_adata = tea_atac_adata[common_tea,:]
for i,adata in enumerate([cite_adt_adata,cite_rna_adata,tea_adt_adata,tea_rna_adata,tea_atac_adata]):
    adata.obs.index.name = None
    adata.var.index.name = None
    adata.write('./processed/adata_{}.h5ad'.format(i))

#View of AnnData object with n_obs × n_vars = 41520 × 52
#View of AnnData object with n_obs × n_vars = 41520 × 32738
#View of AnnData object with n_obs × n_vars = 36771 × 52
#View of AnnData object with n_obs × n_vars = 36771 × 32738
#View of AnnData object with n_obs × n_vars = 36771 × 533227


'''get pca for rna and run harmony'''
cite_rna_adata = sc.read('processed/cite_rna_adata.h5ad')
tea_rna_adata = sc.read('processed/tea_rna_adata.h5ad')
merge_rna_adata = ad.concat([cite_rna_adata,tea_rna_adata],axis=0,join='outer',merge='first',label='batch',keys=['cite','tea'])
identity = []
for item in merge_rna_adata.obs_names:
    if 'HSC' in item:
        identity.append('HSC')
    elif 'LMPP' in item:
        identity.append('LMPP')
    elif 'MPP' in item:
        identity.append('MPP')
    elif '34' in item:
        identity.append('34') 
assert len(identity) == merge_rna_adata.shape[0]
merge_rna_adata.obs['identity'] = identity
merge_rna_adata.obs['concat'] = [b+'_'+i for b,i in zip(merge_rna_adata.obs['batch'],merge_rna_adata.obs['identity'])]

def get_pca_for_rna(adata):
    from scipy.sparse import csr_matrix,issparse
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=3000)
    adata.raw = adata
    adata = adata[:,adata.var['highly_variable']]
    sc.pp.regress_out(adata,['total_counts','pct_counts_mt'])
    sc.pp.scale(adata,max_value=10)
    sc.tl.pca(adata,n_comps=50)
    adata = adata.raw.to_adata()
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    return adata

    
merge_rna_adata = get_pca_for_rna(merge_rna_adata)
ho = hm.run_harmony(merge_rna_adata.obsm['X_pca'],merge_rna_adata.obs,['batch'])
merge_rna_adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(merge_rna_adata,use_rep='X_pca_harmony')
sc.tl.umap(merge_rna_adata)
umap_dual_view_save(merge_rna_adata,cols=['batch','identity','concat'])
merge_rna_adata.write('merge_rna_adata.h5ad')

'''additional analysis on the merge_rna_adata'''
merge_rna_adata = sc.read('merge_rna_adata.h5ad')
# pd.DataFrame(data=merge_rna_adata.obsm['X_umap'],index=merge_rna_adata.obs_names,columns=['umap_x','umap_y']).to_csv('merge_rna_adata_coord.txt',sep='\t')
add_annotations(merge_rna_adata,inputs='../TEA-Seq_ND251/scTriangulate/two-TF-IDF/groups.scTriangulate-R2.txt',cols_input=['label'],index_col=0,cols_output=['sctri_tea'],kind='disk')
merge_rna_adata_tea = merge_rna_adata[merge_rna_adata.obs['batch']=='tea',:]
# umap_dual_view_save(merge_rna_adata_tea,cols=['batch','identity','sctri_tea'])
merge_rna_adata_cite = merge_rna_adata[merge_rna_adata.obs['batch']=='cite',:]
knn_X_train = merge_rna_adata_tea.obsm['X_pca_harmony']
knn_X_test = merge_rna_adata_cite.obsm['X_pca_harmony']
knn_Y_train = merge_rna_adata_tea.obs['sctri_tea'].values
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder
le = LabelEncoder().fit(knn_Y_train)
knn_Y_train_encoded = le.transform(knn_Y_train)
model = KNeighborsClassifier(n_neighbors=15,weights='distance')
model.fit(knn_X_train,knn_Y_train_encoded)
knn_Y_test = le.inverse_transform(model.predict(knn_X_test))
merge_rna_adata_cite.obs['sctri_tea'] = knn_Y_test
# umap_dual_view_save(merge_rna_adata_cite,cols=['batch','identity','sctri_tea'])
merge_rna_adata_after_transfer = ad.concat([merge_rna_adata_cite,merge_rna_adata_tea],axis=0,join='outer',merge='first')
umap_dual_view_save(merge_rna_adata_after_transfer,cols=['batch','identity','sctri_tea'])
print(merge_rna_adata_after_transfer)
merge_rna_adata_after_transfer.write('merge_rna_adata_after_transfer.h5ad')
merge_rna_adata_after_transfer.obs['sctri_tea'].to_csv('after_transfer_sctri_tea_label.txt',sep='\t')






