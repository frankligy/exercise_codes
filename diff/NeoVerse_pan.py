#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from tqdm import tqdm
import mygene
import subprocess
import re
import argparse
import pickle
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import bisect
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from collections import Counter
import anndata as ad
from scipy.sparse import csr_matrix
from collections import Counter
from ast import literal_eval
import math
import argparse
import json
from copy import deepcopy

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'




def integrate_hla():
    data = []
    for s in samples:
        optitype_tmp = os.path.join(ROOTDIR,s,'optitype_tmp')
        f = subprocess.run('find {} -type f -name "*.tsv"'.format(optitype_tmp),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0]
        df = pd.read_csv(f,sep='\t',index_col=0)
        result = pd.DataFrame(index=df.columns[:6],data={'hla_allele':df.iloc[0,:6].values,'total_read_mapped_HLA_locus':np.full(shape=6,fill_value=df.iat[0,6])})
        df = result
        data.append([s] + df['hla_allele'].tolist() + [df.iat[0,1]])
    final = pd.DataFrame.from_records(data=data,columns=['sample','A1','A2','B1','B2','C1','C2','total_read_mapped_HLA_locus']).set_index(keys='sample')
    final.to_csv(os.path.join(OUTDIR,'hla_types.txt'),sep='\t')
    return final



def get_canonical_transcript():
    gtf = pd.read_csv(VARIANT_ENSEMBL_GTF,sep='\t',skiprows=5,header=None)
    gtf.columns = ['chrom','source','feature','start','end','score','strand','phase','attribute']
    gtf_gene = gtf.loc[gtf['feature']=='gene',:]
    pat1 = re.compile(r'gene_id "(ENSG\d+)";')
    pat2 = re.compile(r'gene_name "(.+?)";')
    dic1 = {}
    dic_strand = {}
    for row in gtf_gene.itertuples():    
        strand = row.strand
        ensg = re.search(pat1,row.attribute).group(1)
        try:
            symbol = re.search(pat2,row.attribute).group(1)
        except:   # lncRNA, no gene_name
            symbol = 'unknown'
        dic1[ensg] = symbol
        dic_strand[ensg] = strand
    
    gtf_transcript = gtf.loc[gtf['feature']=='transcript',:]
    pat1 = re.compile(r'gene_id "(ENSG\d+)";')
    pat2 = re.compile(r'transcript_id "(ENST\d+)";')
    dic2 = {}
    for row in gtf_transcript.itertuples():
        if 'Ensembl_canonical' in row.attribute:
            ensg = re.search(pat1,row.attribute).group(1)
            enst = re.search(pat2,row.attribute).group(1)
            dic2[enst] = ensg
    
    with open(os.path.join(OUTDIR,'canonical.txt'),'w') as f:
        f.write('ENST\tENSG\tsymbol\tstrand\n')
        for k,v in dic2.items():
            f.write('{}\t{}\t{}\t{}\n'.format(k,v,dic1[v],dic_strand[v]))

    return dic1,dic2



def integrate_gene(mode,ensg=None,symbol=None,tpm=0,plot_type='boxplot+boxplot',cat_dic=None,gene_path=None):
    if mode == 'compute':
        # process gtex
        gtex = pd.read_csv(GTEX_GENE,sep='\t',skiprows=2,index_col=0)
        gtex = gtex.loc[~gtex.index.str.contains('_PAR'),:]
        gtex.index = [item.split('.')[0] for item in gtex.index]
        gtex.drop(labels=['Cells - Cultured fibroblasts', 'Cells - EBV-transformed lymphocytes'],axis=1,inplace=True)
        gtex.drop(labels=['Ovary','Prostate','Testis','Vagina','Adrenal Gland','Cervix - Endocervix','Cervix - Ectocervix','Fallopian Tube','Pituitary'],axis=1,inplace=True)
        ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
        max_median = gtex.iloc[:,1:].values.max(axis=1)  # n_gene
        ensg2value = pd.Series(index=gtex.index.tolist(),data=max_median).to_dict()

        # process biotype and membrane type
        df = pd.read_csv(GTF_GENE,sep='\t',skiprows=5,header=None)
        df = df.loc[df.iloc[:,2]=='gene',:]
        attr = df.iloc[:,-1].tolist()
        pat1 = re.compile(r'gene_id "(ENSG\d+)\..+"')
        pat2 = re.compile(r'gene_type "(.+?)"')
        ensg2biotype = {}
        for item in attr:
            ensg = re.search(pat1,item).group(1)
            biotype = re.search(pat2,item).group(1)
            ensg2biotype[ensg] = biotype
        
        df = pd.read_csv(MEMBRANE_GENE,sep='\t')
        membrane_ensg = set(df['Ens'].tolist())

        # load the gene
        final = pd.read_csv(gene_path,sep='\t',index_col=0)
        median = np.median(final.values,axis=1)  # n_gene
        result = pd.DataFrame(data={'median_tumor':median},index=final.index)
        result['max_median_gtex'] = result.index.map(ensg2value).fillna(-1).values
        result = result.loc[result['max_median_gtex']!=-1,:]   # only focus on ensg that has gtex entry
        result['median_tumor'] = [1e-5 if item==0 else item for item in result['median_tumor']]
        result['max_median_gtex'] = [1e-5 if item==0 else item for item in result['max_median_gtex']]
        result['logFC'] = np.log2(result['median_tumor'].values / result['max_median_gtex'].values)
        result['gene_symbol'] = result.index.map(ensg2symbol).values

        nu = pd.read_csv(NUORF,sep='\t',index_col=0)
        t_lnc_pseudo = nu.loc[(nu['plotType'].isin(['Pseudogene','lncRNA'])) & (nu['geneType']!='protein_coding'),:]['geneId'].tolist()
        data = []
        for item in t_lnc_pseudo:
            if isinstance(item,str) and item.startswith('ENSG'):
                ensg = item.split('.')[0]
                data.append(ensg)
        data = set(data)   # 7310

        result['translatable_lnc_pseudo'] = [True if item in data else False for item in result.index]
        result['biotype'] = result.index.map(ensg2biotype).fillna('unknown').values
        result['is_membrane'] = result.index.isin(membrane_ensg)

        result.to_csv(os.path.join(atlas_path,'gene_lfc.txt'),sep='\t')
    
    elif mode == 'plot':
        # process gtex
        gtex = pd.read_csv(GTEX_GENE,sep='\t',skiprows=2,index_col=0)
        cond = ~gtex.columns.isin(['Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts'])
        gtex = gtex.loc[:,cond]
        gtex.index = [item.split('.')[0] for item in gtex.index]
        ensg2symbol = pd.Series(index=gtex.index.tolist(),data=gtex['Description'].tolist()).to_dict()
        series = gtex.loc[ensg,:].iloc[1:]

        # process tumor cohort
        final = pd.read_csv(gene_path,sep='\t',index_col=0)
        tumor_expr = final.loc[ensg,:].values

        # plot
        t_pt, n_pt = plot_type.split('+')
        fig = plt.figure(figsize=(15,6))
        gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.1,0.9),wspace=0.2)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1],sharey=ax1)
        if t_pt == 'swarmplot':
            sns.swarmplot(data=tumor_expr,color='red',ax=ax1)
        elif t_pt == 'boxplot':
            sns.boxplot(data=tumor_expr,color='red',ax=ax1)
        elif t_pt == 'swarmplot_cat':
            sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
            from colors import pick_n_colors
            select_colors = pick_n_colors(len(cat_dic))
            for i,(cat,samples) in enumerate(cat_dic.items()):
                samples = [item + ',TPM' for item in samples]
                tumor_expr = final.loc[ensg,samples].values
                sns.swarmplot(data=tumor_expr,color=select_colors[i],ax=ax1)
            import matplotlib.lines as mlines
            ax1.legend(handles=[mlines.Line2D([],[],marker='o',linestyle='',color=i) for i in select_colors],labels=list(cat_dic.keys()),bbox_to_anchor=(0,1),loc='upper right') 
        elif t_pt == 'swarmplot_cat_box':
            sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
            from colors import pick_n_colors
            select_colors = pick_n_colors(len(cat_dic))
            for i,(cat,samples) in enumerate(cat_dic.items()):
                samples = [item + ',TPM' for item in samples]
                tumor_expr = final.loc[ensg,samples].values
                sns.swarmplot(data=tumor_expr,color=select_colors[i],ax=ax1)
            import matplotlib.lines as mlines
            ax1.legend(handles=[mlines.Line2D([],[],marker='o',linestyle='',color=i) for i in select_colors],labels=list(cat_dic.keys()),bbox_to_anchor=(0,1),loc='upper right')
            tumor_expr = final.loc[ensg,:].values
            bp = ax1.boxplot(x=[tumor_expr],positions=[0],patch_artist=False)
            

        ax1.set_ylabel('TPM')
        ax1.set_xlabel('tumor')

        if n_pt == 'barplot':
            ax2.bar(x=np.arange(len(series)),height=series.values,width=0.8,color='green')
        elif n_pt == 'boxplot':
            # # build a h5ad
            # tmp_df = pd.read_csv(GTEX_GENE_ALL,sep='\t',index_col=0,skiprows=2).iloc[:,1:]
            # adata = ad.AnnData(X=csr_matrix(tmp_df.values),obs=pd.DataFrame(index=tmp_df.index),var=pd.DataFrame(index=tmp_df.columns))
            # meta = pd.read_csv(GTEX_META,sep='\t',index_col=0)
            # common_samples = list(set(adata.var_names).intersection(set(meta.index)))
            # meta = meta.loc[common_samples,:]
            # mapping = meta['SMTSD'].to_dict()
            # adata.var['tissue'] = [mapping.get(item,'unknown') for item in adata.var_names]
            # adata.write('gtex_gene_all.h5ad')

            # use the h5ad
            adata = ad.read_h5ad(GTEX_GENE_ALL_H5AD)  # 56200 × 17382
            adata.obs_names = [item.split('.')[0] for item in adata.obs_names]
            adata.obs_names_make_unique()
            adata_gene = adata[[ensg],:]
            normal_expr_list = []
            for t in series.index.tolist():
                values = adata_gene[:,adata_gene.var['tissue']==t].X.toarray().reshape(-1)
                normal_expr_list.append(values)
            bp = ax2.boxplot(x=normal_expr_list,positions=np.arange(len(normal_expr_list)),patch_artist=True)
            for flier in bp['fliers']:
                flier.set_markersize(1)
                flier.set_marker('o')
            for box in bp['boxes']:
                box.set_facecolor('green')
                box.set_edgecolor('black')
                box.set_linewidth(1)

        ax2.set_xticks(np.arange(len(series)))
        ax2.set_xticklabels(series.index.tolist(),fontsize=10,rotation=90)
        ax2.set_xlabel('GTEx Normal')
        
        fig.suptitle('{},{}'.format(ensg,symbol))
        plt.savefig(os.path.join(atlas_path,'{}_{}_expr_{}.{}'.format(ensg,symbol,plot_type,image_format)),bbox_inches='tight')
        plt.close()

    elif mode == 'fasta':
        mat = pd.read_csv(os.path.join(OUTDIR,'genes.txt'),sep='\t',index_col=0)
        meta = pd.read_csv(os.path.join(OUTDIR,'genes_lfc.txt'),sep='\t',index_col=0)
        ensg_protein_coding = set(meta.loc[meta['biotype']=='protein_coding',:].index)
        ensg2symbol,enst2ensg = get_canonical_transcript()
        ensg2enst = {v:k for k,v in enst2ensg.items()}
        canonical = set(enst2ensg.keys())
        ensg2protein = {}
        with open(PROTEIN,'r') as in_handle:
            for title,seq in SimpleFastaParser(in_handle):
                pat1 = r'gene:(ENSG\d+)\.\d+'
                pat2 = r'transcript:(ENST\d+)\.\d+'
                ensg = re.search(pat1,title).group(1)
                enst = re.search(pat2,title).group(1)
                if enst in canonical and ensg in ensg_protein_coding:
                    ensg2protein[ensg] = seq
        with open(os.path.join(OUTDIR,'ensembl_protein.fasta'),'w') as f:   # 19116 entries
            for ensg,seq in ensg2protein.items():
                f.write('>{}|{}|{}\n{}\n'.format(ensg,ensg2enst[ensg],ensg2symbol[ensg],seq))
        for c in mat.columns:
            sample = c.rstrip(',TPM')
            series = mat[c]
            expressed_ensg = set(series.loc[series>tpm].index)
            ensg2tpm = series.to_dict()
            valid_ensg = expressed_ensg.intersection(set(ensg2protein.keys()))
            print('# valid for {}: {}'.format(sample,len(valid_ensg)))
            with open(os.path.join(OUTDIR,db_fasta_dir,sample,'tmp_gene.fasta'),'w') as f:
                for ensg in valid_ensg:
                    f.write('>{}|{}|{}|tpm:{}\n'.format(ensg,ensg2enst[ensg],ensg2symbol[ensg],ensg2tpm[ensg]))
                    f.write('{}\n'.format(ensg2protein[ensg]))



def integrate_fusion(mode):
    if mode == 'compute':
        whole_data = []
        for s in samples:
            data = []
            star_fusion_tmp = os.path.join(ROOTDIR,s,'star_fusion_tmp')
            df = pd.read_csv(os.path.join(star_fusion_tmp,'StarFusionOut','star-fusion.fusion_predictions.abridged.coding_effect.tsv'),sep='\t')
            df.rename(columns={'#FusionName':'FusionName'},inplace=True)
            for row in df.itertuples():
                jrc = row.JunctionReadCount
                src = row.SpanningFragCount
                st = row.SpliceType
                lg = row.LeftGene
                lb = row.LeftBreakpoint
                rg = row.RightGene
                rb = row.RightBreakpoint
                ldas = row.LargeAnchorSupport
                ffpm = row.FFPM
                anno = row.annots
                coding_effect = row.FUSION_CDS
                cds_left = row.CDS_LEFT_ID
                cds_right = row.CDS_RIGHT_ID
                data.append((','.join([lg,lb,rg,rb]),jrc,src,st,ldas,ffpm,anno,coding_effect,cds_left,cds_right))
            result = pd.DataFrame.from_records(data,columns=['feature','jrc','src','st','ldas','ffpm','anno','coding_effect','cds_left','cds_right']).set_index(keys='feature')
            df = result.reset_index()
            whole_data.append(df)
        final = pd.concat(whole_data,axis=0,keys=samples).reset_index(level=-2)
        final.rename(columns={'level_0':'sample'},inplace=True)
        col = []
        for item in final['feature']:
            lg,lb,rg,rb = item.split(',')
            lgs = lg.split('^')[0]  # even when only ENSG is available, this will still return the ensg, so good
            rgs = rg.split('^')[0]
            fusion = '-'.join([lgs,rgs])
            col.append(fusion)
        final['fusion'] = col
        final.reset_index(drop=True,inplace=True)
        final.to_csv(os.path.join(OUTDIR,'fusion.txt'),sep='\t')

    elif mode == 'fasta':
        # remove existing tmp_fusion.fasta
        for s in samples:
            tmp_file = os.path.join(OUTDIR,db_fasta_dir,s,'tmp_fusion.fasta')
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        # now generate fasta
        df = pd.read_csv(os.path.join(OUTDIR,'fusion.txt'),sep='\t',index_col=0)
        red_herrings = set(['GTEx_recurrent_StarF2019','BodyMap','DGD_PARALOGS','HGNC_GENEFAM','Greger_Normal','Babiceanu_Normal','ConjoinG'])
        for row in df.itertuples():
            cds = row.coding_effect
            anno = row.anno
            fusion = row.fusion
            cds_left = row.cds_left
            cds_right = row.cds_right
            sample = row.sample
            ffpm = row.ffpm
            cond = True
            for item in anno:
                if item.replace('"','') in red_herrings:
                    cond = False
            if cds != '.' and cond:  
                for i in range(len(cds)-1):
                    if cds[i].islower() and cds[i+1].isupper():
                        bp_index = i  # the python index for the last nt in the first section
                        first_section = cds[:i+1]
                        second_section = cds[i+1:]
                        break
                phase = (bp_index + 1) % 3  # exon phase explanation: https://www.biostars.org/p/63864/, lend {phase} nt
                '''
                this phase is not the same as how exon phase is defined, phase = 2 means actual phase = 1
                so if phase = 2, lend 1 nt to the last codon, so you only need to retrieve 3 * {MAX_PEP_LEN - 1} + 2 nt in the first section, and 3 * {MAX_PEP_LEN - 1} + 1 nt in the second section
                so if phase = 1, lend 2 nt to the last codon, so you only need to retrieve 3 * {MAX_PEP_LEN - 1} + 1 nt in the first section, and 3 * {MAX_PEP_LEN - 1} + 2 nt in the second section
                so if phase = 0, lend 0(3) nt to the last codon, so you only need dto retrive 3 * {MAX_PEP_LEN - 1} in the first section, and 3 * {MAX_PEP_LEN - 1} nt in the second section
                '''
                if phase == 2:
                    n_need_1 = 3 * (MAX_PEP_LEN - 1) + 2
                    n_need_2 = 3 * (MAX_PEP_LEN - 1) + 1
                    first_retrieve = first_section[-n_need_1:]
                    second_retrieve = second_section[:n_need_2]
                    retrieve = first_retrieve + second_retrieve
                elif phase == 1:
                    n_need_1 = 3 * (MAX_PEP_LEN - 1) + 1
                    n_need_2 = 3 * (MAX_PEP_LEN - 1) + 2
                    first_retrieve = first_section[-n_need_1:]
                    second_retrieve = second_section[:n_need_2]
                    retrieve = first_retrieve + second_retrieve
                elif phase == 0:
                    n_need_1 = 3 * (MAX_PEP_LEN - 1) 
                    n_need_2 = 3 * (MAX_PEP_LEN - 1) 
                    first_retrieve = first_section[-n_need_1:]
                    second_retrieve = second_section[:n_need_2]
                    retrieve = first_retrieve + second_retrieve
                try:
                    assert len(retrieve) % 3 == 0
                except AssertionError:
                    print('starfusion error, cds generated not a multiple of three, skipped')
                    print(sample,fusion,bp_index,retrieve)
                    continue
                else:
                    mini_protein = str(Seq(retrieve).translate(to_stop=True))  # because the starfusion cds include the stop codon
                    with open(os.path.join(OUTDIR,db_fasta_dir,sample,'tmp_fusion.fasta'),'a') as f:
                        f.write('>{}|{}|{}|ffpm:{}\n{}\n'.format(fusion,cds_left,cds_right,ffpm,mini_protein))





def integrate_pathogen(mode,n_cutoff=0,f_cutoff=0,included_strains=None,pathogen_path=None,total_sample=None,pathogen_rec_path=None,plot_strain=None,ylim=None):
    if mode == 'compute':
        df = pd.read_csv(pathogen_path,sep='\t')
        df = df.loc[df['taxonomy'].str.contains('|s__',regex=False),:]
        df = df.loc[df['count']>=n_cutoff,:]
        df = df.loc[~df['taxonomy'].str.contains('|s__Homo sapiens',regex=False),:]

        final = df
        vc = final['taxonomy'].value_counts()
        f_cutoff = f_cutoff * total_sample
        recurrent_taxonomy = set(vc.loc[vc > f_cutoff].index.tolist())
        final = final.loc[final['taxonomy'].isin(recurrent_taxonomy),:]
        final['strain'] = [item.split('|s__')[1] for item in final['taxonomy']]  # frequent taxonomy only and extracted strain

        if final.shape[0] > 0:
            order = {s:np.median(sub_df['count'].values) for s,sub_df in final.groupby(by='strain')}
            order = dict(sorted(order.items(), key=lambda item: item[1],reverse=True))
            order = list(order.keys())
            fig,ax = plt.subplots(figsize=(6.4,12))
            sns.stripplot(data=final,y='strain',x='count',size=8,order=order,ax=ax)
            ax.tick_params(axis='y', labelsize=5)
            ax.tick_params(axis='x', bottom=False,top=True,labelbottom=False,labeltop=True)
            plt.savefig(os.path.join(atlas_path,'pathogen_frequency.{}'.format(image_format)),bbox_inches='tight')
            plt.close()
        
        # # pathogen blacklist
        # meta = pd.read_csv(T2T_NORMAL_META,sep='\t',index_col=0,header=None)
        # meta.columns = ['srr']
        # tissue2srr = meta['srr'].to_dict()
        # normal_pathogen_list = []
        # all_srrs = []
        # all_tissues = []
        # for t,ss in tissue2srr.items():
        #     for s in ss.split(','):
        #         kraken2_result_path = os.path.join(T2T_NORMAL_DIR,s,'test_report.txt')
        #         kraken2_df = pd.read_csv(kraken2_result_path,sep='\t',header=None)
        #         kraken2_df = kraken2_df.loc[kraken2_df[0].str.contains('|s__',regex=False),:]
        #         normal_pathogen_list.append(kraken2_df)
        #         all_srrs.append(s)
        #         all_tissues.append(t)
        # all_uids = [s+','+t for s,t in zip(all_srrs,all_tissues)]
        # normal_pathogen = pd.concat(normal_pathogen_list,axis=0,keys=all_uids).reset_index(level=-2)
        # normal_pathogen.columns = ['uids','strain','count']
        # normal_pathogen['sample'] = [item.split(',')[0] for item in normal_pathogen['uids']]
        # normal_pathogen['tissue'] = [item.split(',')[1] for item in normal_pathogen['uids']]
        # normal_pathogen.to_csv(NORMAL_PATHOGEN_PATH,sep='\t',index=None)

        pathogen_rec = vc.to_frame()
        pathogen_rec.columns =['count']
        pathogen_rec['strain'] = [item.split('|s__')[1] for item in pathogen_rec.index]
        normal_pathogen = pd.read_csv(NORMAL_PATHOGEN_PATH,sep='\t')
        pathogen_dict = {}
        for strain,sub_df in normal_pathogen.groupby(by='strain'):
            pathogen_dict[strain] = {}
            for tissue,sub_df2 in sub_df.groupby(by='tissue'):
                per_tissue_strain = ','.join([str(item) for item in sub_df2['count'].values])
                pathogen_dict[strain][tissue] = per_tissue_strain
        pathogen_rec['normal'] = pathogen_rec.index.map(pathogen_dict).values
        pathogen_rec.to_csv(os.path.join(atlas_path,'pathogen_rec.txt'),sep='\t')

    elif mode == 'plot':
        normal_pathogen = pd.read_csv(NORMAL_PATHOGEN_PATH,sep='\t')
        normal_pathogen['strain'] = [item.split('|s__')[1] for item in normal_pathogen['strain']]
        sample2uq = {}
        for sample,sub_df in normal_pathogen.groupby(by='sample'):
            all_counts = sub_df['count'].values
            all_counts = all_counts[all_counts!=0]
            upper_quantile = np.quantile(all_counts,0.75)
            sample2uq[sample] = upper_quantile
        
        subset = normal_pathogen.loc[normal_pathogen['strain']==plot_strain,:]
        subset.set_index(keys='sample',inplace=True)
        subset['uq'] = subset.index.map(sample2uq).values
        subset['normalized_value'] = np.log2(subset['count'].values / subset['uq'].values)
        normal_normalized_all_values = {t:sub_df['normalized_value'].values.tolist() for t,sub_df in subset.groupby(by='tissue')}

        pathogen_path=os.path.join(atlas_path,'pathogen_all.txt')
        df = pd.read_csv(pathogen_path,sep='\t')
        df = df.loc[df['taxonomy'].str.contains('|s__',regex=False),:]
        df = df.loc[~df['taxonomy'].str.contains('|s__Homo sapiens',regex=False),:]
        df['strain'] = [item.split('|s__')[1] for item in df['taxonomy']] 
        sample2uq = {}
        for sample,sub_df in df.groupby(by='sample'):
            all_counts = sub_df['count'].values
            all_counts = all_counts[all_counts!=0]
            upper_quantile = np.quantile(all_counts,0.75)
            sample2uq[sample] = upper_quantile
        subset = df.loc[df['strain']==plot_strain,:]
        subset.set_index(keys='sample',inplace=True)
        subset['uq'] = subset.index.map(sample2uq).values
        subset['normalized_value'] = np.log2(subset['count'].values / subset['uq'].values)
        tumor_normalized_all_values = subset['normalized_value'].values.tolist()

        
        # plot
        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.1,0.9),wspace=0.2)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1],sharey=ax1)  
        sns.boxplot(data=tumor_normalized_all_values,color='red',ax=ax1)    
        ax1.set_ylabel('value')
        ax1.set_xlabel('tumor')
        if ylim is not None:
            ax1.set_ylim(ylim)
        if len(normal_normalized_all_values) > 0:
            bp = ax2.boxplot(x=list(normal_normalized_all_values.values()),positions=np.arange(len(normal_normalized_all_values)),patch_artist=True)
            for flier in bp['fliers']:
                flier.set_markersize(1)
                flier.set_marker('o')
            for box in bp['boxes']:
                box.set_facecolor('green')
                box.set_edgecolor('black')
                box.set_linewidth(1)
            ax2.set_xticks(np.arange(len(normal_normalized_all_values)))
            ax2.set_xticklabels(list(normal_normalized_all_values.keys()),fontsize=10,rotation=90)
            ax2.set_xlabel('GTEx Normal')
        fig.suptitle(plot_strain)
        plt.savefig(os.path.join(atlas_path,'{}.{}'.format(plot_strain,image_format)),bbox_inches='tight')
        plt.close()
        



    elif mode == 'fasta':
        OUTDIR = atlas_path
        final = pd.read_csv(os.path.join(OUTDIR,'pathogen_rec.txt'),sep='\t')
        species = included_strains
        ref_proteome_table = pd.read_csv(REF_PROTEOME,sep='\t',index_col=0)
        other_proteome_table = pd.read_csv(OTHER_PROTEOME,sep='\t',index_col=0)
        proteome_table = pd.concat([ref_proteome_table,other_proteome_table],axis=0)
        proteome_table = proteome_table.loc[proteome_table['Organism'].notna(),:]
        valid_dict = {}
        for s in species:
            hits = proteome_table.loc[proteome_table['Organism'].str.contains(s,regex=False),:]
            if hits.shape[0] == 1:
                upid = hits.iloc[0,:].name
                valid_dict[s] = upid
            elif hits.shape[0] == 0:
                continue
            else:
                custom_order = ['Standard','Close to standard (high value)','Close to standard (low value)','Outlier (high value)','Outlier (low value)','Unknown']
                hits['CPD'] = pd.Categorical(hits['CPD'], categories=custom_order, ordered=True)
                col = []
                for item in hits['BUSCO']:
                    if isinstance(item,str):
                        completeness = float(item.split('[')[0].lstrip('C:').rstrip('%'))
                    else:
                        completeness = 0
                    col.append(completeness)
                hits['completness'] = col
                hits.sort_values(by=['CPD','completness'],ascending=[True,False],inplace=True)
                upid = hits.iloc[0,:].name
                valid_dict[s] = upid
        for s,upid in valid_dict.items():
            url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A{}%29%29'.format(upid)
            download_file_name = os.path.join(OUTDIR,'{}_{}.fasta.gz'.format(s.replace(' ','_'),upid))
            subprocess.run(['curl','-o',download_file_name,url])
            subprocess.run(['gunzip',download_file_name])

        valid_species = list(set(species).intersection(set(valid_dict.keys())))
        defacto_final = final.loc[final['strain'].isin(valid_species)]
        for row in defacto_final.itertuples():
            strain = row.strain
            fasta_file_name = os.path.join(OUTDIR,'{}_{}.fasta'.format(strain.replace(' ','_'),valid_dict[strain]))
            subprocess.run(['cp',fasta_file_name,os.path.join(OUTDIR,db_fasta_dir)])
        for s,upid in valid_dict.items():
            subprocess.run("rm {}".format(os.path.join(OUTDIR,'{}_{}.fasta'.format(s.replace(' ','_'),upid))),shell=True)




def check_orientation(pileupread,strand,library):
    read = pileupread.alignment
    if library == 'first_strand':
        strand = '+' if strand == '+' else '-'
    elif library == 'second_strand':
        strand = '-' if strand == '+' else '+'


    if strand == '+':
        # F2R1
        if read.is_read1 and read.is_reverse:   
            cond = True
        elif read.is_read2 and not read.is_reverse:
            cond = True
        else:
            cond = False
    elif strand == '-':
        # F1R2
        if read.is_read1 and not read.is_reverse:   
            cond = True
        elif read.is_read2 and read.is_reverse:
            cond = True
        else:
            cond = False    

    return cond 

def count_n_correct_reads(sample,strand,library,variant,chrom,pos):
    samfile = pysam.AlignmentFile(os.path.join(ROOTDIR,sample,'star_tmp','{}_secondAligned.sortedByCoord.out.bam'.format(sample)),'rb') 
    counter = 0
    pos = int(pos)
    for pileupcolumn in samfile.pileup(chrom,pos-1,pos,truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip and check_orientation(pileupread,strand,library):
                if pileupread.alignment.query_sequence[pileupread.query_position] == variant:
                    counter += 1
    return counter

def integrate_variant(mode,library='first_strand',gene_path=None,tpm=0,mutation_rec_path=None,protein_path=None,fs_dict=None):
    if mode == 'fasta':
        variants = pd.read_csv(mutation_rec_path,sep='\t')
        variants = variants.loc[(variants['ensg'].str.contains(r'^ENSG')) &                  
                                (~variants['mutation'].str.contains(r'^HLA')) &
                                (~variants['mutation'].str.contains(r'^IG')) &
                                (~variants['mutation'].str.contains(r'^TRAV')) & 
                                (~variants['mutation'].str.contains(r'^TRAJ')) &
                                (~variants['mutation'].str.contains(r'^TRBV')) & 
                                (~variants['mutation'].str.contains(r'^TRBD')) & 
                                (~variants['mutation'].str.contains(r'^TRBJ')) &
                                ((variants['n_samples']>1) | (variants['is_driver'])), :]   

        # load the hg38 seq for frameshift
        dict_fa = {}  # {chrom:seq}
        with open(HG38_SEQ,'r') as in_handle:
            for title,seq in SimpleFastaParser(in_handle):
                dict_fa[title] = seq

        # start to generate peptide
        ori_gene_fasta = {}
        variants_fasta = {}
        with open(protein_path,'r') as in_handle:
            for title,seq in SimpleFastaParser(in_handle):
                ensg,enst,symbol = title.split('|')
                ori_gene_fasta[ensg] = (title,seq)
        for row in variants.itertuples():
            if row.type == 'missense_variant':
                try:
                    info = ori_gene_fasta[row.ensg]
                except:
                    continue
                else:
                    aa_title,aa_seq = ori_gene_fasta[row.ensg]
                    pat = r'@p\.([ARNDCQEGHILKMFPSTWYV])(\d+)([ARNDCQEGHILKMFPSTWYV])'
                    match = re.search(pat,row.mutation)
                    ref_aa = match.group(1)
                    pos_aa = int(match.group(2))
                    alt_aa = match.group(3)
                    try:
                        documented_ref_aa = aa_seq[pos_aa-1]    # maybe the position is not based on canonical isoform
                    except:
                        continue
                    try:
                        assert documented_ref_aa == ref_aa    # correct
                    except AssertionError:
                        if documented_ref_aa == alt_aa:   # genome version difference
                            continue
                        else:           # weird edge case, do not consider to improve precision
                            continue
                    if MAX_PEP_LEN > pos_aa:
                        actual_needed = pos_aa
                    else:
                        actual_needed = MAX_PEP_LEN
                    variant_seq = aa_seq[pos_aa-1-(actual_needed-1):pos_aa-1] + alt_aa + aa_seq[pos_aa:pos_aa+(MAX_PEP_LEN-1)]
                    variants_fasta[row.mutation] = (row.n_samples,row.median_dna_vaf,row.ensg,row.genetic,row.type,variant_seq)

            elif row.type == 'inframe_insertion':  # all very low
                try:
                    info = ori_gene_fasta[row.ensg]
                except:
                    continue
                else:
                    aa_title,aa_seq = ori_gene_fasta[row.ensg]
                    if row.mutation.endswith('dup'):
                        gene,ce = row.mutation.split('@')
                        ce = ce.split('dup')[0].split('p.')[1]
                        if '_' in ce:
                            first,second = ce.split('_')
                            first = int(first[1:])
                            second = int(second[1:])
                        else:
                            first = int(ce[1:])
                            second = first
                        preceding = aa_seq[:first-1]
                        following = aa_seq[second:]
                        impacted = aa_seq[first-1:second]
                        updated = impacted + impacted
                        n_aa_needed = MAX_PEP_LEN - len(updated)
                        variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]
                        variants_fasta[row.mutation] = (row.n_samples,row.median_dna_vaf,row.ensg,row.genetic,row.type,variant_seq)
                    elif 'ins' in row.mutation:
                        gene,ce = row.mutation.split('@')
                        ce,insert = ce.split('ins')
                        ce = ce.split('p.')[1]
                        first,second = ce.split('_')
                        first = int(first[1:])
                        second = int(second[1:])
                        preceding = aa_seq[:first-1]
                        following = aa_seq[second:]
                        updated = aa_seq[first-1] + insert + aa_seq[second-1]
                        n_aa_needed = MAX_PEP_LEN - len(updated)
                        variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]
                        variants_fasta[row.mutation] = (row.n_samples,row.median_dna_vaf,row.ensg,row.genetic,row.type,variant_seq)
            elif row.type == 'inframe_deletion':
                try:
                    info = ori_gene_fasta[row.ensg]
                except:
                    continue
                else:
                    aa_title,aa_seq = ori_gene_fasta[row.ensg]
                    gene,ce = row.mutation.split('@')
                    if 'delins' in ce:
                        ce,delins = ce.split('delins')
                        ce = ce.split('p.')[1]
                        first,second = ce.split('_')
                        first = int(first[1:])
                        second = int(second[1:])
                        preceding = aa_seq[:first-1]
                        following = aa_seq[second:]
                        updated = delins
                        n_aa_needed = MAX_PEP_LEN - len(updated)
                        variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]
                        variants_fasta[row.mutation] = (row.n_samples,row.median_dna_vaf,row.ensg,row.genetic,row.type,variant_seq)
                    elif ce.endswith('del'):
                        ce = ce.split('del')[0]
                        ce = ce.split('p.')[1]
                        if '_' in ce:
                            first,second = ce.split('_')
                            first = int(first[1:])
                            second = int(second[1:])
                        else:      
                            first = int(ce[1:])
                            second = int(ce[1:])
                        preceding = aa_seq[:first-1]
                        following = aa_seq[second:]
                        updated = ''
                        n_aa_needed = MAX_PEP_LEN - 1
                        variant_seq = preceding[-n_aa_needed:] + updated + following[:n_aa_needed]
                        variants_fasta[row.mutation] = (row.n_samples,row.median_dna_vaf,row.ensg,row.genetic,row.type,variant_seq)
            elif row.type == 'frameshift_variant':
                try:
                    info = ori_gene_fasta[row.ensg]
                except:
                    continue
                else:
                    aa_title,aa_seq = ori_gene_fasta[row.ensg]
                    gene,ce = row.mutation.split('@')
                    if gene not in fs_dict.keys():
                        continue
                    coord,replace = row.genetic.split(';') 
                    full_coord = coord
                    before, after = replace.split('/')
                    chrom,coord = coord.split(':')
                    start,end = coord.split('-')
                    start,end = int(start),int(end)
                    if full_coord not in fs_dict[gene].keys():
                        continue
                    strand, mode, n_codon, p_codon = fs_dict[gene][full_coord]  

                    '''
                    example fs_dict, idea is the following second stretch will be derived from dna, and the first nt will start to translate
                    fs_dict = {
                        'NPM1':{'chr5:171410539-171410540':['+','mode1',287,1]},
                        'DNMT3A':{'chr2:25247662-25247662':['-','mode2',314,3]},
                        'TP53':{'chr17:7670685-7670685':['-','mode2',341,3]},
                    }
                    '''

                    if mode == 'mode2':
                        start = start - 1
                        end = end + 1
                        after = ''
                    preceding = aa_seq[:n_codon-1]
                    n_aa_needed = MAX_PEP_LEN - 1
                    first_part = preceding[-n_aa_needed:]
                    if strand == '+':
                        stretch = dict_fa[chrom][start-1-(p_codon-1):start] + after + dict_fa[chrom][end-1:end-1+1000]
                    elif strand == '-':
                        stretch = dict_fa[chrom][start-1-1000:start] + after + dict_fa[chrom][end-1:end+(p_codon-1)]
                        stretch = str(Seq(stretch).reverse_complement())
                    second_part = str(Seq(stretch).translate(to_stop=False)).split('*')[0]
                    if len(second_part) > 0:
                        variant_seq = first_part + second_part
                        variants_fasta[row.mutation] = (row.n_samples,row.median_dna_vaf,row.ensg,row.genetic,row.type,variant_seq)



        # final write out
        with open(os.path.join(atlas_path,db_fasta_dir,'mutation.fasta'),'w') as f:
            for k,v in variants_fasta.items():
                n_samples, median_dna_vaf, ensg, genetic, typ, seq = v
                f.write('>{}|{}|{}|{}|{}|{}\n{}\n'.format(k.replace('@','|'),str(n_samples),str(round(median_dna_vaf,2)),ensg,genetic.replace(';','|'),typ,seq))





        


    elif mode == 'compute':
        # merge the vcf files, make sure htslib/1.3, samtools are loaded
        tmp_dir = os.path.join(OUTDIR,'NeoVerse_tmp_variant_concat')
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        
        for s in samples:
            file_path = os.path.join(ROOTDIR,s,'variant_tmp','variants.vcf')
            new_path = os.path.join(tmp_dir,'{}.vcf'.format(s))
            cmd = 'bcftools filter -i \'FILTER="PASS"\' -o {} {}'.format(new_path,file_path)  # choose precision over sensitivity per Opossum paper
            subprocess.run(cmd,shell=True)
            subprocess.run("bgzip {}".format(new_path),shell=True)
            subprocess.run("tabix -p vcf {}".format(new_path+'.gz'),shell=True)

        samples_string = ''
        for s in samples:
            samples_string += '{}.vcf.gz '.format(s)
        old_dir = os.getcwd()
        os.chdir(tmp_dir)
        subprocess.run("bcftools merge --force-samples {}>merged.vcf".format(samples_string),shell=True)

        with open('samples.txt','w') as f:
            for i,s in enumerate(samples):
                if i == 0:
                    f.write('output {}\n'.format(s))
                else:
                    index = i+1
                    f.write('{}:output {}\n'.format(index,s))
        subprocess.run("bcftools reheader -s samples.txt -o merged_rename.vcf merged.vcf",shell=True)
        subprocess.run("rm merged.vcf",shell=True)
        subprocess.run("rm *.vcf.gz*",shell=True)
        os.chdir(old_dir)

        # liftover to hg38 and then run VEP to get coding effect mutation
        with open(os.path.join(tmp_dir,'merged_rename.vcf'),'r') as f:
            count = 0
            for line in f:
                if line.startswith('##'):
                    count += 1
                elif line.startswith('#CHROM'):
                    break
        skiprows = count

        vcf = pd.read_csv(os.path.join(tmp_dir,'merged_rename.vcf'),sep='\t',skiprows=skiprows)
        vcf.rename(columns={'#CHROM':'CHROM'},inplace=True)
        with open(os.path.join(tmp_dir,'prelift.bed'),'w') as f:
            for row in vcf.itertuples():
                uid = 'var{}'.format(row.Index+1)
                chrom = row.CHROM
                pos = row.POS
                f.write('{}\t{}\t{}\t{}\n'.format(chrom,int(pos)-1,pos,uid))  # bed file is 0-based and half-open
        subprocess.run("{} {} {} {} {}".format(LIFTOVER_PROGRAM,os.path.join(tmp_dir,'prelift.bed'),LIFTOVER_CHAIN['hs1_hg38'],os.path.join(tmp_dir,'postlift.bed'),os.path.join(tmp_dir,'unmapped.txt')),shell=True)
        mapping = pd.read_csv(os.path.join(tmp_dir,'postlift.bed'),sep='\t',header=None,index_col=3)
        mapping.columns = ['CHROM','POS_1','POS']
        valid_vcf = vcf.iloc[[int(item.split('var')[1])-1 for item in mapping.index],:]
        valid_vcf['CHROM'] = mapping['CHROM'].values
        valid_vcf['POS'] = mapping['POS'].values
        valid_chrom = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
        valid_vcf = valid_vcf.loc[valid_vcf['CHROM'].isin(set(valid_chrom)),:]
        valid_vcf.rename(columns={'CHROM':'#CHROM'},inplace=True)
        valid_vcf.to_csv(os.path.join(tmp_dir,'tmp_vcf.vcf'),sep='\t',index=None)
        subprocess.run("head -n {} {} > {}".format(skiprows,os.path.join(tmp_dir,'merged_rename.vcf'),os.path.join(tmp_dir,'tmp.header')),shell=True)
        subprocess.run("cat {} {} > {}".format(os.path.join(tmp_dir,'tmp.header'),os.path.join(tmp_dir,'tmp_vcf.vcf'),os.path.join(tmp_dir,'to_vep.vcf')),shell=True)
        subprocess.run("rm {}".format(os.path.join(tmp_dir,'tmp*')),shell=True)

        '''
        VEP is very tricky, what ended up working is
        module load singularity/3.9.8  # must be this one, not the 3.1 version
        singularity pull --name vep.sif docker://ensemblorg/ensembl-vep
        singularity exec vep.sif INSTALL.pl -c $HOME/vep_data -a cf -s homo_sapiens -y GRCh38  # you must install that to $HOME/vep_data, I keep a copy in the working dir
        '''
        subprocess.run("bcftools sort -o {} {}".format(os.path.join(HOME_DIR,'to_vep_sorted.vcf'),os.path.join(tmp_dir,'to_vep.vcf')),shell=True)
        subprocess.run("singularity exec {} vep --dir {} -i {} --cache --offline".format(VEP_SIF,os.path.join(HOME_DIR,'vep_data'),os.path.join(HOME_DIR,'to_vep_sorted.vcf')),shell=True)
        if not os.path.exists(os.path.join(HOME_DIR,'vep_output')):
            subprocess.run("mkdir {}".format(os.path.join(HOME_DIR,'vep_output')),shell=True)
        subprocess.run("cp {} {}".format(os.path.join(HOME_DIR,'variant_effect_output*'),os.path.join(HOME_DIR,'vep_output')),shell=True)
        subprocess.run("rm {}".format(os.path.join(HOME_DIR,'to_vep_sorted.vcf')),shell=True)
        subprocess.run("rm {}".format(os.path.join(HOME_DIR,'variant_effect_output*')),shell=True)
        subprocess.run("mv {} {}".format(os.path.join(HOME_DIR,'vep_output'),tmp_dir),shell=True)

        # get missense, deletion and insertion
        cmd = 'awk \'$7=="missense_variant" {{print}}\' {} > {}'.format(os.path.join(tmp_dir,'vep_output','variant_effect_output.txt'),os.path.join(tmp_dir,'vep_output','variant_effect_output_missense_variant.txt'))
        subprocess.run(cmd,shell=True)
        cmd = 'awk \'$7=="inframe_deletion" {{print}}\' {} > {}'.format(os.path.join(tmp_dir,'vep_output','variant_effect_output.txt'),os.path.join(tmp_dir,'vep_output','variant_effect_output_inframe_deletion.txt'))
        subprocess.run(cmd,shell=True)
        cmd = 'awk \'$7=="inframe_insertion" {{print}}\' {} > {}'.format(os.path.join(tmp_dir,'vep_output','variant_effect_output.txt'),os.path.join(tmp_dir,'vep_output','variant_effect_output_inframe_insertion.txt'))
        subprocess.run(cmd,shell=True)

        _,_ = get_canonical_transcript()
        canonical = pd.read_csv(os.path.join(OUTDIR,'canonical.txt'),sep='\t')

        for category in ['missense_variant','inframe_deletion','inframe_insertion']:
            df = pd.read_csv(os.path.join(tmp_dir,'vep_output','variant_effect_output_{}.txt'.format(category)),sep='\t',header=None)
            df.columns = ['Uploaded_variation','Location','Allele','Gene','Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','Extra']
            df = df.loc[df['Feature'].isin(set(canonical['ENST'])),:]   # 10250 out of 10792 genes have effect on canonical isoform, so its' ok to drop non-canonical isoform
            not_multi_allele = []
            for item in df['Uploaded_variation']:
                if len(item.split('/')) > 2:
                    not_multi_allele.append(False)
                else:
                    not_multi_allele.append(True)
            df = df.loc[not_multi_allele,:]
            df.to_csv(os.path.join(tmp_dir,'vep_output','variant_effect_output_{}.txt'.format(category)),sep='\t',index=None)

        '''
        for insertion
            ** chr1_234609200_-/GCT means 234609199 A->AGCT
            ** protein 16-17 means put the new one between 16 and 17, or if there's number, then just find that

        for deletion
            ** chr1_1295458_GCCCCTGCG/- means chr1	1295457	AGCCCCTGCG	A
            ** protein just means where it was deleted

        For missense
            ** easy
        '''
        dic = {}  # {chr1:942176:ENSG00000187634_Gcc/Tcc_A467S}
        for category in ['missense_variant','inframe_insertion','inframe_deletion']:
            df = pd.read_csv(os.path.join(tmp_dir,'vep_output','variant_effect_output_{}.txt'.format(category)),sep='\t')
            for row in df.itertuples():
                ref,alt = row.Amino_acids.split('/')
                pos = str(row.Protein_position)
                assembly = ref + pos + alt
                if category == 'missense_variant':
                    key = row.Location
                elif category == 'inframe_insertion':
                    key = row.Location.split('-')[0]
                elif category == 'inframe_deletion':
                    tmp = row.Location.split('-')[0]
                    tmp1,tmp2 = tmp.split(':')
                    correct_tmp2 = str(int(tmp2) - 1)
                    key = ':'.join([tmp1,correct_tmp2])
                dic[key] = '_'.join([row.Gene,row.Codons,assembly]) 

        mapping['uid'] = [':'.join([item1,str(item2)]) for item1,item2 in zip(mapping['CHROM'],mapping['POS'])]
        mapping['coding_effect'] = mapping['uid'].map(dic).fillna('invalid').values
        mapping = mapping.loc[mapping['coding_effect']!='invalid',:]

        valid_vcf = vcf.iloc[[int(item.split('var')[1])-1 for item in mapping.index],:]
        valid_vcf['coding_effect'] = mapping['coding_effect'].values
        valid_vcf['INFO'] = [';'.join([item1,'CE={}'.format(item2)]) for item1,item2 in zip(valid_vcf['INFO'],valid_vcf['coding_effect'])]
        valid_vcf.rename(columns={'CHROM':'#CHROM'},inplace=True)
        valid_vcf.drop(columns=['coding_effect'],inplace=True)
        valid_vcf.to_csv(os.path.join(tmp_dir,'tmp_vcf.vcf'),sep='\t',index=None)
        subprocess.run("head -n {} {} > {}".format(skiprows,os.path.join(tmp_dir,'merged_rename.vcf'),os.path.join(tmp_dir,'tmp.header')),shell=True)
        subprocess.run("sed -i '{}i\{}' {}".format(skiprows,'##INFO=<ID=CE,Number=1,Type=String,Description="Coding Effect for this variant">',os.path.join(tmp_dir,'tmp.header')),shell=True)
        subprocess.run("cat {} {} > {}".format(os.path.join(tmp_dir,'tmp.header'),os.path.join(tmp_dir,'tmp_vcf.vcf'),os.path.join(tmp_dir,'merged_rename_ce.vcf')),shell=True)
        subprocess.run("rm {}".format(os.path.join(tmp_dir,'tmp*')),shell=True)


        # add dbSNP
        ## subprocess.run("bcftools query -f '%CHROM\t%POS\t%{}\n' {} | bgzip -c > {}".format('FREQ',VARIANT_DBSNP,os.path.join(tmp_dir,'annot.txt.gz')),shell=True)
        ## subprocess.run("tabix -s1 -b2 -e2 {}".format(os.path.join(tmp_dir,'annot.txt.gz')),shell=True)
        subprocess.run("echo -e '##INFO=<ID={},Number=.,Type=String,Description=\"Population Frequency\">' > {}".format('FREQ',os.path.join(tmp_dir,'hdr.txt')),shell=True)
        subprocess.run("bcftools annotate -a {} -h {} -c CHROM,POS,INFO/FREQ -o {} {}".format(VARIANT_DBSNP_PROCESSED,os.path.join(tmp_dir,'hdr.txt'),os.path.join(tmp_dir,'merged_rename_ce_dbsnp.vcf'),os.path.join(tmp_dir,'merged_rename_ce.vcf')),shell=True)


        # add REDIT, first liftover REDIT, then annotate
        ## with open(VARIANT_REDIT,'r') as f1, open(os.path.join(tmp_dir,'redit_prelift.bed'),'w') as f2:
        ##     next(f1)
        ##     i = 1
        ##     for line in tqdm(f1):
        ##         line = line.rstrip('\n')
        ##         contig,pos,ref,alt,strand = line.split('\t')[:5]
        ##         col1 = contig  # chrom
        ##         col2 = int(pos) - 1 # start, 0 based
        ##         col3 = int(pos)
        ##         col4 = col1 + ':' + str(col3)   # original hg38 pos for cross-reference later
        ##         f2.write('{}\t{}\t{}\t{}\n'.format(col1,col2,col3,col4))
        ##         i += 1

        ## subprocess.run("{} {} {} {} {}".format(LIFTOVER_PROGRAM,os.path.join(tmp_dir,'redit_prelift.bed'),LIFTOVER_CHAIN['hg38_hs1'],os.path.join(tmp_dir,'redit_postlift.bed'),os.path.join(tmp_dir,'unmapped.txt')),shell=True)
        ## subprocess.run("sort -k1V -k2n -k3n {} | bgzip -c > {}.gz".format(os.path.join(tmp_dir,'redit_postlift.bed'),os.path.join(tmp_dir,'redit_postlift.bed')),shell=True)
        ## subprocess.run("tabix -s 1 -b 2 -e 3 {}.gz".format(os.path.join(tmp_dir,'redit_postlift.bed')),shell=True)
        subprocess.run("echo -e '##INFO=<ID=REDI,Number=1,Type=String,Description=\"siteID from REDIportal\">' > {}".format(os.path.join(tmp_dir,'tmp.header')),shell=True)
        subprocess.run("bcftools annotate -a {} -h {} -c CHROM,FROM,TO,REDI -o {} {}".format(VARIANT_REDIT_PROCESSED,os.path.join(tmp_dir,'tmp.header'),os.path.join(tmp_dir,'merged_rename_ce_dbsnp_redit.vcf'),os.path.join(tmp_dir,'merged_rename_ce_dbsnp.vcf')),shell=True)

        # annotate by cosmic
        # cosmic = pd.read_csv(VARIANT_COSMIC,sep='\t')
        # cosmic = cosmic.loc[cosmic['MUTATION_DESCRIPTION']=='missense_variant',:]
        # cosmic = cosmic.loc[cosmic['MUTATION_SOMATIC_STATUS'].isin(['Confirmed somatic variant','Reported in another cancer sample as somatic']),:]
        # with open('cosmic_prelift.bed','w') as f:
        #     for item1,item2 in zip(cosmic['CHROMOSOME'],cosmic['GENOME_START']):
        #             col1 = 'chr' + str(item1)  # chrom
        #             col2 = int(item2) - 1 # start, 0 based
        #             col3 = int(item2)
        #             col4 = col1 + ':' + str(col3)   # original hg38 pos for cross-reference later
        #             f.write('{}\t{}\t{}\t{}\n'.format(col1,col2,col3,col4))
        # subprocess.run("{} {} {} {} {}".format(LIFTOVER_PROGRAM,'cosmic_prelift.bed',LIFTOVER_CHAIN['hg38_hs1'],'cosmic_postlift.bed','unmapped.txt'),shell=True)
        # subprocess.run("sort -k1V -k2n -k3n {} | bgzip -c > {}.gz".format('cosmic_postlift.bed','cosmic_postlift.bed'),shell=True)
        # subprocess.run("tabix -s 1 -b 2 -e 3 {}.gz".format('cosmic_postlift.bed')),shell=True)
        subprocess.run("echo -e '##INFO=<ID=COSMIC,Number=1,Type=String,Description=\"siteID from COSMIC\">' > {}".format(os.path.join(tmp_dir,'tmp.header')),shell=True)
        subprocess.run("bcftools annotate -a {} -h {} -c CHROM,FROM,TO,COSMIC -o {} {}".format(VARIANT_COSMIC_PROCESSED,os.path.join(tmp_dir,'tmp.header'),os.path.join(tmp_dir,'merged_rename_ce_dbsnp_redit_cosmic.vcf'),os.path.join(tmp_dir,'merged_rename_ce_dbsnp_redit.vcf')),shell=True)
        subprocess.run("rm {}".format(os.path.join(tmp_dir,'tmp.header')),shell=True)
        subprocess.run("rm {}".format(os.path.join(tmp_dir,'merged_rename_ce.vcf')),shell=True)
        subprocess.run("rm {}".format(os.path.join(tmp_dir,'merged_rename_ce_dbsnp.vcf')),shell=True)
        subprocess.run("rm {}".format(os.path.join(tmp_dir,'merged_rename_ce_dbsnp_redit.vcf')),shell=True)

        # analyze the data
        with open(os.path.join(tmp_dir,'merged_rename_ce_dbsnp_redit_cosmic.vcf'),'r') as f:
            count = 0
            for line in f:
                if line.startswith('##'):
                    count += 1
                elif line.startswith('#CHROM'):
                    break
        skiprows = count

        df = pd.read_csv(os.path.join(tmp_dir,'merged_rename_ce_dbsnp_redit_cosmic.vcf'),sep='\t',skiprows=skiprows)
        col = [] # dbsnp pass loose
        col1 = [] # dbsnp pass stringent
        for item in df['INFO']:
            dic = dict(pair.split('=') for pair in item.split(';'))
            try:
                freq = dic['FREQ']
            except KeyError:
                col.append(True)
                col1.append(True)
            else:
                col1.append(False)
                freq_dic = dict(pair.split(':') for pair in freq.split('|'))
                ref_freq_list = []
                for k,v in freq_dic.items():
                    ref_freq = float(v.split(',')[0])
                    ref_freq_list.append(ref_freq)
                if np.all(np.array(ref_freq_list) > 0.999):  # somatic mutation in dbsnp usually MAF < 0.001
                    col.append(True)
                else:
                    col.append(False)
        df['pass_dbsnp_stringent'] = col1
        df['pass_dbsnp_loose'] = col

        col = [] # redit pass
        for item in df['INFO']:
            dic = dict(pair.split('=') for pair in item.split(';'))
            try:
                redit = dic['REDI']
            except KeyError:
                col.append(True)
            else:
                col.append(False)
        df['pass_redit'] = col

        col = [] # recurrent
        sub_df = df.loc[:,samples]
        empty = './.:.:.:.:.:.'
        df['n_sample'] = sub_df.apply(lambda x:len(x.loc[x!=empty]),axis=1,result_type='reduce').values

        col = [] # cosmic pass
        for item in df['INFO']:
            dic = dict(pair.split('=') for pair in item.split(';'))
            try:
                redit = dic['COSMIC']
            except KeyError:
                col.append(False)
            else:
                col.append(True)
        df['cosmic_check'] = col        


        # most covered sample name
        def get_max_nr_sample(x):
            empty = './.:.:.:.:.:.'
            x = x.loc[x!=empty]
            sub_samples = x.index.tolist()
            lis = []
            for item in x:
                try:
                    lis.append(int(item.split(':')[-2]))
                except:
                    lis.append(int(item.split(':')[-2].split(',')[0]))  # two variant allele, like 5,5 or 2,2, for now just take the first one
            max_sample = sub_samples[np.argmax(lis)]
            return max_sample
        col = [] # sample name
        sub_df = df.loc[:,samples]
        df['max_sample'] = sub_df.apply(lambda x:get_max_nr_sample(x),axis=1,result_type='reduce').values


        def get_nr(x):
            empty = './.:.:.:.:.:.'
            x = x.loc[x!=empty]
            lis = []
            for item in x:
                try:
                    lis.append(int(item.split(':')[-2]))
                except:
                    lis.append(int(item.split(':')[-2].split(',')[0]))  # two variant allele, like 5,5 or 2,2, for now just take the first one
            lis = [str(i) for i in lis]
            return ','.join(lis)

        col = [] #  NR/TC
        sub_df = df.loc[:,samples]
        df['nr_tc'] = sub_df.apply(lambda x:get_nr(x),axis=1,result_type='reduce').values

        def get_nv(x):
            empty = './.:.:.:.:.:.'
            x = x.loc[x!=empty]
            lis = []
            for item in x:
                try:
                    lis.append(int(item.split(':')[-1]))
                except:
                    lis.append(int(item.split(':')[-1].split(',')[0]))  # two variant allele, like 5,5 or 2,2, for now just take the first one
            lis = [str(i) for i in lis]
            return ','.join(lis)

        col = [] #  NV/TR
        sub_df = df.loc[:,samples]
        df['nv_tr'] = sub_df.apply(lambda x:get_nv(x),axis=1,result_type='reduce').values


        def get_clonality(x):
            empty = './.:.:.:.:.:.'
            x = x.loc[x!=empty]
            lis_nr = []
            for item in x:
                try:
                    lis_nr.append(int(item.split(':')[-2]))
                except:
                    lis_nr.append(int(item.split(':')[-2].split(',')[0]))  # two variant allele, like 5,5 or 2,2, for now just take the first one

            lis_nv = []
            for item in x:
                try:
                    lis_nv.append(int(item.split(':')[-1]))
                except:
                    lis_nv.append(int(item.split(':')[-1].split(',')[0]))  # two variant allele, like 5,5 or 2,2, for now just take the first one

            clonality = (np.array(lis_nv) / np.array(lis_nr)).tolist()
            clonality = [str(round(item,2)) for item in clonality]

            return ','.join(clonality)

        col = [] #  clonality
        sub_df = df.loc[:,samples]
        df['clonality'] = sub_df.apply(lambda x:get_clonality(x),axis=1,result_type='reduce').values

        col1 = [] # max clonality
        col2 = [] # median clonality
        for item in df['clonality']:
            cs = np.array(item.split(','))
            cs = [float(item) for item in cs]
            col1.append(np.max(cs))
            col2.append(np.median(cs))
        df['max_clonality'] = col1
        df['median_clonality'] = col2
            
        # annotate based on gene
        gene = pd.read_csv(gene_path,sep='\t',index_col=0)
        ensg2median = gene['median_tumor'].to_dict()
        ensg2symbol = gene['gene_symbol'].to_dict()
        col = [] # gene
        for item in df['INFO']:
            dic = dict(pair.split('=') for pair in item.split(';'))
            col.append(dic['CE'].split('_')[0])
        df['ensg'] = col
        df['median'] = df['ensg'].map(ensg2median).fillna(0).values
        df['symbol'] = df['ensg'].map(ensg2symbol).fillna('unknown').values

        # annotate based on strand
        canonical = pd.read_csv(os.path.join(OUTDIR,'canonical.txt'),sep='\t')
        ensg2strand = pd.Series(index=canonical['ENSG'].values,data=canonical['strand'].values).to_dict()
        df['strand'] = df['ensg'].map(ensg2strand).fillna('unknown').values

        df.to_csv(os.path.join(OUTDIR,'merged_rename_ce_dbsnp_redit_cosmic_annotate.txt'),sep='\t',index=None)

        # check strandedness
        df = pd.read_csv(os.path.join(OUTDIR,'merged_rename_ce_dbsnp_redit_cosmic_annotate.txt'),sep='\t')
        df.rename(columns={'#CHROM':'CHROM'},inplace=True)
        if library == 'first_strand' or library == 'second_strand':
            col = []
            for row in tqdm(df.itertuples(),total=df.shape[0]):
                counter = count_n_correct_reads(row.max_sample,row.strand,library,row.ALT,row.CHROM,row.POS)
                col.append(counter)
            df['strandedness_library_check'] = col
        elif library == 'no':
            df['strandedness_library_check'] = np.full(df.shape[0],fill_value=50)
        # further get rid of multi-alleles
        col = [False if ',' in item else True for item in df['ALT']]
        df = df.loc[col,:]
        df.to_csv(os.path.join(OUTDIR,'variants.txt'),sep='\t',index=None)

        # filter
        cond = df['cosmic_check'] | (df['pass_dbsnp_loose'] & df['pass_dbsnp_stringent'] & df['pass_redit'] &(df['max_clonality'] < 0.95) & (df['median'] > 1))
        df = df.loc[cond,:].sort_values(by='n_sample',ascending=False)
        df.to_csv(os.path.join(OUTDIR,'variants_possible_somatic.txt'),sep='\t',index=None)







def annotate_splicing_by_te(vc_file,te_gtf):
    # build the dic for TE, stranded and unstranded version
    te_df = pd.read_csv(TE_GTF_HG38,sep='\t',header=None)
    te_df.columns = ['chrom','source','feature','start','end','score','strand','phase','attribute']

    chrom_mapping = {
        1:'chr1',
        2:'chr2',
        3:'chr3',
        4:'chr4',
        5:'chr5',
        6:'chr6',
        7:'chr7',
        8:'chr8',
        '9':'chr9',
        10:'chr10',
        11:'chr11',
        12:'chr12',
        13:'chr13',
        14:'chr14',
        15:'chr15',
        16:'chr16',
        17:'chr17',
        18:'chr18',
        19:'chr19',
        20:'chr20',
        21:'chr21',
        22:'chr22',
        'X':'chrX',
        'Y':'chrY'
    }

    te_df['chrom'] = te_df['chrom'].map(chrom_mapping).fillna('unknown').values
    te_df = te_df.loc[te_df['chrom']!='unknown',:]


    big_dict = {}
    for chrom,sub_df in te_df.groupby(by='chrom'):
        big_dict[chrom] = {}
        for strand,sub_df2 in sub_df.groupby(by='strand'):
            big_dict[chrom][strand] = [[],[]]   # first element stores the start and end, second element stores the info
            for row in sub_df2.itertuples():
                exons = (int(row.start),int(row.end))
                gid,tid,fid,cid,_ = row.attribute.split(';')
                info = [gid.split(' "')[1].rstrip('"'),tid.split(' "')[1].rstrip('"'),fid.split(' "')[1].rstrip('"'),cid.split(' "')[1].rstrip('"')]
                big_dict[chrom][strand][0].extend([exons[0],exons[1]])
                big_dict[chrom][strand][1].append(info)
    big_dict['chrM'] = {'+':[[],[]],'-':[[],[]]}


    # now start to annotate 
    col = []
    for row in tqdm(vc_file.itertuples(),total=vc_file.shape[0]):
        if row.has_known_ss and row.inferred_strand != 0:
            uid = row.location
            strand = row.inferred_strand
            chrom,coord = uid.split(':')
            start,end = coord.split('-')
            start,end = int(start),int(end)
            ss1, ss2 = row.gene_symbol.split(',')
            if strand == '+' and ss1 == 'None':
                infos = []
                for te_stuff, te_strand in zip([big_dict[chrom]['+'],big_dict[chrom]['-']],['+','-']):
                    lm_pos = bisect.bisect_left(te_stuff[0],start) # if odd, good, 3 means fall into second regions, left most if already present
                    if lm_pos % 2 == 1:
                        i_region = (lm_pos + 1) // 2 - 1  # minus 1 because python 0 based
                        info = deepcopy(te_stuff[1][i_region])
                        info.append('donor')
                        info.append(te_strand)
                        block = [te_stuff[0][lm_pos-1],start]
                        block = [str(item) for item in block]
                        info.extend(block)
                        info = ','.join(info)
                        infos.append(info)
                col.append(';'.join(infos))


            elif strand == '+' and ss2 == 'None':
                infos = []
                for te_stuff, te_strand in zip([big_dict[chrom]['+']],['+']):  # no anti-sense because it is acceptor
                    rm_pos = bisect.bisect_right(te_stuff[0],end)   # right most
                    if rm_pos % 2 == 1:
                        i_region = (rm_pos + 1) // 2 - 1  
                        info = deepcopy(te_stuff[1][i_region])
                        info.append('acceptor')   # acceptor
                        info.append(te_strand)
                        block = [end,te_stuff[0][rm_pos]]   # not pos+1 any more, just pos
                        block = [str(item) for item in block]
                        info.extend(block)
                        info = ','.join(info)
                        infos.append(info)
                col.append(';'.join(infos))

            elif strand == '-' and ss1 == 'None':
                infos = []
                for te_stuff,te_strand in zip([big_dict[chrom]['-']],['-']):  # same as acceptor, but - strand
                    lm_pos = bisect.bisect_left(te_stuff[0],start)   # back to left most
                    if lm_pos % 2 == 1:
                        i_region = (lm_pos + 1) // 2 - 1
                        info = deepcopy(te_stuff[1][i_region])
                        info.append('acceptor')  
                        info.append(te_strand)
                        block = [te_stuff[0][lm_pos-1],start]   # back to left most
                        block = [str(item) for item in block]
                        info.extend(block)
                        info = ','.join(info)
                        infos.append(info)
                col.append(';'.join(infos))

            elif strand == '-' and ss2 == 'None':
                infos = []
                for te_stuff,te_strand in zip([big_dict[chrom]['-'],big_dict[chrom]['+']],['-','+']):  # need to consider anti-sense again
                    rm_pos = bisect.bisect_right(te_stuff[0],end)   # back to right most
                    if rm_pos % 2 == 1:
                        i_region = (rm_pos + 1) // 2 - 1
                        info = deepcopy(te_stuff[1][i_region])
                        info.append('donor')  
                        info.append(te_strand)
                        block = [end,te_stuff[0][rm_pos]]   # back to right most
                        block = [str(item) for item in block]
                        info.extend(block)
                        info = ','.join(info)
                        infos.append(info)
                col.append(';'.join(infos))

            else:   # known ss, do a broader round of te annotation, no
                col.append('')
            

        else:
            col.append('')


    vc_file['te_info'] = col
    return vc_file


def integrate_splicing(mode,events=None,splicing_path=None,total_sample=None):

    if mode == 'fasta':
        vc = pd.read_csv(splicing_path,sep='\t',index_col=0)
        cond = (vc['n_sample'] > total_sample * 0.2) & (vc['ave_count_normal'] < 1) & (vc['ave_count_tumor'] > 10) & (vc['logFC'] > 4) & (vc['has_known_ss'])
        vc = vc.loc[cond,:]

        dict_fa = {}  # {chrom:seq}
        with open(HG38_SEQ,'r') as in_handle:
            for title,seq in SimpleFastaParser(in_handle):
                dict_fa[title] = seq

        j2idb = pd.read_csv(junction2peptide_db,sep='\t',index_col=0)
        j2idb = j2idb.loc[j2idb['is_valid'],:]
        j2i = j2idb['protein'].to_dict()
        j2enst = j2idb['enst'].to_dict()
        all_j = set(j2idb.index)

        # now start to get sequence
        junction2peptide = {}
        for item,strand,te_info,g_symbol in tqdm(zip(vc.index,vc['inferred_strand'],vc['te_info'],vc['gene_symbol']),total=vc.shape[0]):
            info = []
            chrom,coord = item.split(':')
            start,end = coord.split('-')
            start,end = int(start),int(end)
            overhang = 100
            if isinstance(te_info,float):   # for some reason, still coerced to float(nan)

                if strand == '+':
                    first_section = dict_fa[chrom][start-overhang:start]
                    second_section = dict_fa[chrom][end-1:end-1+overhang]
                elif strand == '-':
                    old_first_section = dict_fa[chrom][start-overhang:start]
                    old_second_section = dict_fa[chrom][end-1:end-1+overhang]
                    # reverse complement process
                    first_section = Seq(old_second_section).reverse_complement()
                    second_section = Seq(old_first_section).reverse_complement()
            else:
                if ';' in te_info:  # must be single one, no ;
                    for te_element in te_info.split(';'):
                        if te_element.split(',')[5] == strand:
                            te_info = te_element
                            break

                gid,tid,fid,cid,typ,t_strand,t_start,t_end = te_info.split(',')
                t_start,t_end = int(t_start),int(t_end)

                if strand == '+':
                    if typ == 'donor':
                        first_section = dict_fa[chrom][t_start-1:start]
                        second_section = dict_fa[chrom][end-1:end-1+overhang]
                    elif typ == 'acceptor':
                        first_section = dict_fa[chrom][start-overhang:start]
                        second_section = dict_fa[chrom][end-1:t_end]
                elif strand == '-':
                    if typ == 'donor':
                        old_first_section = dict_fa[chrom][start-overhang:start]
                        old_second_section = dict_fa[chrom][end-1:t_end]
                    elif typ == 'acceptor':
                        old_first_section = dict_fa[chrom][t_start-1:start]
                        old_second_section = dict_fa[chrom][end-1:end-1+overhang]  
                    first_section = Seq(old_second_section).reverse_complement()
                    second_section = Seq(old_first_section).reverse_complement()

            
            # decide whether and how to correct
            if not isinstance(te_info,float):
                if typ == 'donor' and len(first_section) > 3 * MAX_PEP_LEN_SPLICING:
                    correct_first = True
                    correct_second = False
                elif typ == 'acceptor' and len(second_section) > 3 * MAX_PEP_LEN_SPLICING:
                    correct_first = False
                    correct_second = True
            else:
                correct_first, correct_second = False,False


            # consider phase = 0
            if not correct_first:
                defacto_first = first_section[-(3*(MAX_PEP_LEN_SPLICING-1)):]
            else:
                CORRECT_MAX_PEP_LEN_SPLICING = len(first_section) // 3 + 1
                defacto_first = first_section[-(3*(CORRECT_MAX_PEP_LEN_SPLICING-1)):]
            if not correct_second:
                defacto_second = second_section[:(3*(MAX_PEP_LEN_SPLICING-1))]
            else:
                CORRECT_MAX_PEP_LEN_SPLICING = len(second_section) // 3 + 1
                defacto_second = second_section[:(3*(CORRECT_MAX_PEP_LEN_SPLICING-1))]
            defacto_junction = defacto_first + defacto_second
            defacto_peptide = str(Seq(defacto_junction).translate(to_stop=False))
            for i,small_pep in enumerate(defacto_peptide.split('*')):
                if len(small_pep) >= MIN_PEP_LEN_SPLICING:
                    info.append((strand,0,small_pep,te_info,g_symbol,i+1))

            # consider phase = 1
            if not correct_first:
                defacto_first = first_section[-(3*(MAX_PEP_LEN_SPLICING-1)+2):]
            else:
                CORRECT_MAX_PEP_LEN_SPLICING = (len(first_section)-2) // 3 + 1
                defacto_first = first_section[-(3*(CORRECT_MAX_PEP_LEN_SPLICING-1)+2):]
            if not correct_second:
                defacto_second = second_section[:(3*(MAX_PEP_LEN_SPLICING-1)+1)]
            else:
                CORRECT_MAX_PEP_LEN_SPLICING = (len(second_section)-1) // 3 + 1
                defacto_second = second_section[:(3*(CORRECT_MAX_PEP_LEN_SPLICING-1)+1)]
            defacto_junction = defacto_first + defacto_second
            defacto_peptide = str(Seq(defacto_junction).translate(to_stop=False))
            for i,small_pep in enumerate(defacto_peptide.split('*')):
                if len(small_pep) >= MIN_PEP_LEN_SPLICING:
                    info.append((strand,1,small_pep,te_info,g_symbol,i+1)) 

            # consider phase = 2
            if not correct_first:
                defacto_first = first_section[-(3*(MAX_PEP_LEN_SPLICING-1)+1):]
            else:
                CORRECT_MAX_PEP_LEN_SPLICING = (len(first_section)-1) // 3 + 1
                defacto_first = first_section[-(3*(CORRECT_MAX_PEP_LEN_SPLICING-1)+1):]
            if not correct_second:
                defacto_second = second_section[:(3*(MAX_PEP_LEN_SPLICING-1)+2)]
            else:
                CORRECT_MAX_PEP_LEN_SPLICING = (len(second_section)-2) // 3 + 1
                defacto_second = second_section[:(3*(CORRECT_MAX_PEP_LEN_SPLICING-1)+2)]
            defacto_junction = defacto_first + defacto_second
            defacto_peptide = str(Seq(defacto_junction).translate(to_stop=False))
            for i,small_pep in enumerate(defacto_peptide.split('*')):
                if len(small_pep) >= MIN_PEP_LEN_SPLICING:
                    info.append((strand,2,small_pep,te_info,g_symbol,i+1)) 
            junction2peptide[item] = info  



        # write to fasta
        with open(os.path.join(atlas_path,db_fasta_dir,'splicing.fasta'),'w') as f:
            for row in vc.itertuples():
                location = row.Index
                junction = '|'.join([location,str(row.n_sample),str(round(row.ave_count_tumor,2)),str(round(row.ave_count_normal,2)),str(round(row.logFC,2))])
                info = junction2peptide[location]
                for item in info:
                    strand = item[0]
                    phase = item[1]
                    seq = item[2]
                    about_te = item[3]
                    g_symbol = item[4]
                    num_id = item[5]
                    if isinstance(about_te,float):
                        if location not in all_j:   # treat as a normal novel splicing junction
                            f.write('>{}|{}|{}|{}|{}\n{}\n'.format(junction,strand,str(phase),g_symbol,str(num_id),seq))
                        else:
                            seq = j2i[location]
                            enst = j2enst[location]
                            f.write('>{}|{}|{}|{}\n{}\n'.format(junction,strand,enst,g_symbol,seq))
                            break   # for this, you only need to do once fo each location
                    else:
                        f.write('>{}|{}|{}|TE_info:{}|{}|{}\n{}\n'.format(junction,strand,str(phase),about_te,g_symbol,str(num_id),seq))
                        

    elif mode == 'plot':
        adata_gtex = ad.read_h5ad(SPLICING_GTEX_HG38)
        adata_tcga = ad.read_h5ad(SPLICING_TCGA_HG38)
        # build gtex_map
        gtex_map = pd.read_csv(SPLICING_GTEX_MAPPING_HG38,sep='\t',header=None)
        gtex_map.columns = ['chrom','start_1','end','uid']
        gtex_mapping = {}
        for row in gtex_map.itertuples():
            gtex_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid
        # build tcga_map
        tcga_map = pd.read_csv(SPLICING_TCGA_MAPPING_HG38,sep='\t',header=None)
        tcga_map.columns = ['chrom','start_1','end','uid']
        tcga_mapping = {}
        for row in tcga_map.itertuples():
            tcga_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid
        # start to plot
        for event in events:
            uid = gtex_mapping.get(event,None)
            if uid is not None:
                adata_gtex_single = adata_gtex[uid,:]
            else:
                adata_gtex_single = ad.AnnData(X=csr_matrix(np.full((1,adata_gtex.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0],'std':[0]},index=[uid]),var=adata_gtex.var)
                print('imputing gtex for {}'.format(event))

            uid = tcga_mapping.get(event,None)
            if uid is not None:
                adata_tcga_single = adata_tcga[uid,:]
            else:
                adata_tcga_single = ad.AnnData(X=csr_matrix(np.full((1,adata_tcga.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0],'std':[0]},index=[uid]),var=adata_tcga.var)
                adata_tcga_single.obs_names = [uid]
                print('imputing tcga for {}'.format(event))

            adata_single = ad.concat([adata_gtex_single,adata_tcga_single],axis=1,join='outer',merge='first')
            fig = plt.figure(figsize=(18,6))
            gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.1,0.9),wspace=0.2)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1],sharey=ax1)
            # for tumor
            all_splicing = pd.read_csv(splicing_path,sep='\t')
            tumor_expr = all_splicing.loc[all_splicing['Unnamed: 0']==event,'count'].values.tolist()
            sns.boxplot(data=tumor_expr,color='red',ax=ax1)
            # for normal
            normal_expr_list = []
            all_tissues = adata_single.var['tissue'].unique()
            for t in all_tissues:
                values = adata_single[:,adata_single.var['tissue']==t].X.toarray().reshape(-1)
                normal_expr_list.append(values)
            bp = ax2.boxplot(x=normal_expr_list,positions=np.arange(len(normal_expr_list)),patch_artist=True)
            for flier in bp['fliers']:
                flier.set_markersize(1)
                flier.set_marker('o')
            for box in bp['boxes']:
                box.set_facecolor('green')
                box.set_edgecolor('black')
                box.set_linewidth(1)
            ax2.set_xticks(np.arange(len(all_tissues)))
            ax2.set_xticklabels(all_tissues,fontsize=10,rotation=90)
            ax2.set_xlabel('Normal')
            fig.suptitle('{}'.format(event))
            plt.savefig(os.path.join(atlas_path,'{}_splicing.{}'.format(event,image_format)),bbox_inches='tight')
            plt.close()
            


            





    elif mode == 'compute':
        # read in two h5ad and the mapping file
        adata_gtex = ad.read_h5ad(SPLICING_GTEX_HG38)
        adata_tcga = ad.read_h5ad(SPLICING_TCGA_HG38)
        # build gtex_map
        gtex_map = pd.read_csv(SPLICING_GTEX_MAPPING_HG38,sep='\t',header=None)
        gtex_map.columns = ['chrom','start_1','end','uid']
        gtex_mapping = {}
        for row in gtex_map.itertuples():
            gtex_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid
        # build tcga_map
        tcga_map = pd.read_csv(SPLICING_TCGA_MAPPING_HG38,sep='\t',header=None)
        tcga_map.columns = ['chrom','start_1','end','uid']
        tcga_mapping = {}
        for row in tcga_map.itertuples():
            tcga_mapping['{}:{}-{}'.format(row.chrom,str(int(row.start_1+1)),row.end)] = row.uid

        # process to only get mean count
        uid2mean_gtex = adata_gtex.obs['mean'].to_dict()
        uid2mean_tcga = adata_tcga.obs['mean'].to_dict()

        # read in the splicing_all
        final = pd.read_csv(splicing_path,sep='\t')
        final.columns = ['location','count','sample']

        # get splicing_rec
        data = []
        for event,sub_df in tqdm(final.groupby(by='location')):
            data.append((event,sub_df.shape[0],sub_df['count'].values.mean()))
        splicing_rec = pd.DataFrame.from_records(data,columns=['location','n_sample','ave_count_tumor'])

        # add normal tissue
        col = []
        for event in tqdm(splicing_rec['location']):
            # retrieve gtex
            uid = gtex_mapping.get(event,None)
            if uid is not None:
                gtex_expr = uid2mean_gtex[uid]
            else:
                gtex_expr = 0
            # retrieve tcga
            uid = tcga_mapping.get(event,None)
            if uid is not None:
                tcga_expr = uid2mean_tcga[uid]
            else:
                tcga_expr = 0
            # combine
            gtex_sf = 2629 / (2629 + 701)
            tcga_sf = 701 / (2629 + 701)
            expr_normal = gtex_expr * gtex_sf + tcga_expr * tcga_sf
            col.append(expr_normal)
        splicing_rec['ave_count_normal'] = col

        # calculate lfc
        splicing_rec['ave_count_tumor'] = [1e-5 if item==0 else item for item in splicing_rec['ave_count_tumor']]
        splicing_rec['ave_count_normal'] = [1e-5 if item==0 else item for item in splicing_rec['ave_count_normal']]
        splicing_rec['logFC'] = np.log2(splicing_rec['ave_count_tumor'].values / splicing_rec['ave_count_normal'].values)

        # infer strand
        gtf = pd.read_csv(VARIANT_ENSEMBL_GTF,sep='\t',skiprows=5,header=None)
        gtf.columns = ['chrom','source','feature','start','end','score','strand','phase','attribute']

        chrom_mapping = {
            1:'chr1',
            2:'chr2',
            3:'chr3',
            4:'chr4',
            5:'chr5',
            6:'chr6',
            7:'chr7',
            8:'chr8',
            9:'chr9',
            10:'chr10',
            11:'chr11',
            12:'chr12',
            13:'chr13',
            14:'chr14',
            15:'chr15',
            16:'chr16',
            17:'chr17',
            18:'chr18',
            19:'chr19',
            20:'chr20',
            21:'chr21',
            22:'chr22',
            'X':'chrX',
            'Y':'chrY',
            'M':'chrM',
        }

        gtf['chrom'] = gtf['chrom'].map(chrom_mapping).fillna('unknown').values
        gtf = gtf.loc[gtf['chrom']!='unknown',:]

        gtf = gtf.loc[gtf['feature']=='exon',:]
        pat = re.compile(r'gene_name "(.+?)"')
        cond = []
        symbol = []
        for item in gtf.attribute:
            match = re.search(pat,item)
            if not match:
                cond.append(False)
                symbol.append(None)
            else:
                gs = match.group(1)
                if 'MSTRG' in gs or gs.startswith('AC') or gs.startswith('AL') or gs.startswith('AP') or '-AS' in gs:
                    cond.append(False)
                    symbol.append(None)
                else:
                    cond.append(True)
                    symbol.append(gs)
        gtf['symbol'] = symbol
        gtf = gtf.loc[cond,:]

        known_splicing_site = {
            'chr1':{},
            'chr2':{},
            'chr3':{},
            'chr4':{},
            'chr5':{},
            'chr6':{},
            'chr7':{},
            'chr8':{},
            'chr9':{},
            'chr10':{},
            'chr11':{},
            'chr12':{},
            'chr13':{},
            'chr14':{},
            'chr15':{},
            'chr16':{},
            'chr17':{},
            'chr18':{},
            'chr19':{},
            'chr20':{},
            'chr21':{},
            'chr22':{},
            'chrX':{},
            'chrY':{},
            'chrM':{},
        }
        for chrom,sub_df in gtf.groupby(by='chrom'):
            for strand,sub_df2 in sub_df.groupby(by='strand'):
                all_start = sub_df2['start'].tolist()
                all_end = sub_df2['end'].tolist()
                known_sites = set(all_start + all_end)
                known_splicing_site[chrom][strand] = known_sites
        # rescue chrM
        known_splicing_site['chrM']['+'] = set()
        known_splicing_site['chrM']['-'] = set()

        # build another map mapping each location to the gene based on GTF
        location2gene = {}
        for chrom,sub_df in tqdm(gtf.groupby(by='chrom')):
            for strand,sub_df2 in sub_df.groupby(by='strand'):
                for symbol,sub_df3 in sub_df2.groupby(by='symbol'):
                    all_start = sub_df3['start'].tolist()
                    all_end = sub_df3['end'].tolist()
                    known_sites = set(all_start + all_end)
                    for site in known_sites:
                        location_id = '_'.join([chrom,strand,str(site)])
                        location2gene[location_id] = symbol

        col = []
        col_strand = []
        for item in splicing_rec['location']:
            chrom,coord = item.split(':')
            start,end = coord.split('-')
            if int(start) in known_splicing_site[chrom]['+'] or int(end) in known_splicing_site[chrom]['+']:
                col.append(True)
                col_strand.append('+')
            elif int(start) in known_splicing_site[chrom]['-'] or int(end) in known_splicing_site[chrom]['-']:
                col.append(True)
                col_strand.append('-')
            else:
                col.append(False)
                col_strand.append('0')
        splicing_rec['has_known_ss'] = col
        splicing_rec['inferred_strand'] = col_strand

        # add gene symbol
        col = []
        for row in splicing_rec.itertuples():
            chrom,coord = row.location.split(':')
            start,end = coord.split('-')
            if row.has_known_ss:
                start_location_id = '_'.join([chrom,row.inferred_strand,str(start)])
                end_location_id = '_'.join([chrom,row.inferred_strand,str(end)])
                start_symbol = location2gene.get(start_location_id,None)
                end_symbol = location2gene.get(end_location_id,None)
                symbol = ','.join([str(start_symbol),str(end_symbol)])
            else:
                symbol = None
            col.append(symbol)

        splicing_rec['gene_symbol'] = col

        # add te_chimeric
        splicing_rec = annotate_splicing_by_te(splicing_rec,TE_GTF_HG38)

        # write out
        splicing_rec.to_csv(os.path.join(atlas_path,'splicing_rec.txt'),sep='\t',index=None)

            






def integrate_erv(mode,ervs=None,mean_normal_cutoff=0.5,median_tumor_cutoff=1,frac_sample=0.15,min_lfc_cutoff=0.58,te_dir_path=None,total_number=None):
    if mode == 'compute':
        # you still manifest
        manifest = pd.read_csv(manifest_path,sep='\t')
        manifest = manifest.loc[manifest['sample_id'].notna(),:]
        tid = [int(item.split('-')[-1][:-1]) for item in manifest['sample_id']]
        cond = [False if item == 10 or item ==11 else True for item in tid]  # 01 primary, 02 recurrent, 03 blood primary, 06 metastatic, 07 additional metastatic,10 blood normal, 11 solid normal
        manifest = manifest.loc[cond,:]
        print('valid number: {}'.format(manifest.shape[0]))
        sample_names = [item.split('.bam')[0] for item in manifest['name']]
        sample_ids = [item for item in manifest['sample_id']]
        mapping_n2i = {n:i for n,i in zip(sample_names,sample_ids)}
        mapping_i2s = {item:'-'.join(item.split('-')[:-1]) for item in sample_ids}

        # combined all tumor ERV
        old_dir = os.getcwd()
        os.chdir(te_dir_path)
        all_files = subprocess.run("for f in *_TElocal_out.cntTable; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
        os.chdir(old_dir)
        tumor_erv_index = None   # contain feature
        tumor_erv_list = {}    # {s:array,}
        aux_list = []
        valid_files = set([item + '_TElocal_out.cntTable' for item in mapping_n2i.keys()])
        reduced_files = []
        for s in tqdm(all_files):
            if s in valid_files:
                reduced_files.append(s)
                telocal_df = pd.read_csv(os.path.join(te_dir_path,s),sep='\t',skiprows=1,header=None)
                telocal_df.columns = ['feature','count']
                total_count = telocal_df['count'].sum()
                tmp = telocal_df['count'].values
                tmp = tmp[tmp!=0]
                upper_quantile = np.quantile(tmp,0.75)
                telocal_erv_df = telocal_df.set_index(keys='feature')
                tumor_erv_index = telocal_erv_df.index.to_list()
                tumor_erv_list[s] = telocal_erv_df['count'].values
                aux_list.append((total_count,upper_quantile))
        tumor_erv = pd.DataFrame(data=tumor_erv_list,index=tumor_erv_index)

        tumor_erv = ad.AnnData(X=csr_matrix(tumor_erv.values),obs=pd.DataFrame(index=tumor_erv.index),var=pd.DataFrame(index=tumor_erv.columns))
        aux_df = pd.DataFrame.from_records(data=aux_list,columns=['total_count','upper_quantile'])
        aux_df.index = reduced_files

        tumor_erv.layers['cpm'] = csr_matrix(tumor_erv.X / aux_df['total_count'].values.reshape(1,-1) * 1e6)
        tumor_erv.write(os.path.join(atlas_path,'tumor_erv.h5ad'))
        aux_df.to_csv(os.path.join(atlas_path,'tumor_erv_aux_df.txt'),sep='\t')

        # # combined all normal ERV
        # normal_erv_list = []
        # meta = pd.read_csv(HG38_NORMAL_META,sep='\t',index_col=0,header=None)
        # meta.columns = ['srr']
        # tissue2srr = meta['srr'].to_dict()
        # all_srrs = []
        # all_tissues = []
        # all_aux = []
        # for t,ss in tqdm(tissue2srr.items()):
        #     for s in ss.split(','):
        #         telocal_result_path = os.path.join(HG38_NORMAL_DIR,s,'TElocal_out.cntTable')
        #         telocal_df = pd.read_csv(telocal_result_path,sep='\t',skiprows=1,header=None)
        #         telocal_df.columns = ['feature','count']
        #         tmp = telocal_df['count'].values
        #         tmp = tmp[tmp!=0]
        #         upper_quantile = np.quantile(tmp,0.75)
        #         total_count = telocal_df['count'].sum()

        #         telocal_erv_df = telocal_df.set_index(keys='feature')
        #         normal_erv_list.append(telocal_erv_df)
        #         all_srrs.append(s)
        #         all_tissues.append(t)
        #         all_aux.append((total_count,upper_quantile))
        # all_uids = [s+','+t for s,t in zip(all_srrs,all_tissues)]
        # normal_erv = pd.concat(normal_erv_list,axis=1,keys=all_uids,join='outer').fillna(value=0)
        # normal_erv.columns = [l1 for l1,l2 in normal_erv.columns.tolist()]
        # normal_erv.to_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv.txt'),sep='\t')
        # normal_aux_df = pd.DataFrame.from_records(data=all_aux,columns=['total_count','upper_quantile'])
        # normal_aux_df.index = all_uids
        # normal_aux_df.to_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv_aux_df.txt'),sep='\t')


        # compare
        tumor_erv = ad.read_h5ad(os.path.join(atlas_path,'tumor_erv.h5ad'))
        tumor_erv_aux = pd.read_csv(os.path.join(atlas_path,'tumor_erv_aux_df.txt'),sep='\t',index_col=0)
        normal_erv = pd.read_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv.txt'),sep='\t',index_col=0)
        normal_erv_aux = pd.read_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv_aux_df.txt'),sep='\t',index_col=0)

        # normalize using aux
        t_data = csr_matrix(tumor_erv.layers['cpm'])
        n_data = normal_erv.values / normal_erv_aux['total_count'].values.reshape(1,-1) * 1e6

        # now start to diff
        median_tumor = np.zeros(t_data.shape[0])
        n_sample_tumor = np.zeros(t_data.shape[0])
        for i in tqdm(range(t_data.shape[0])):
            row_data = t_data[i].toarray().flatten()  
            median_tumor[i] = np.median(row_data)
            n_sample_tumor[i] = np.count_nonzero(row_data)

        mean_normal = np.mean(n_data,axis=1)
        tumor_dict = {}
        normal_dict = {}
        n_sample_dict = {}
        for erv,value in zip(tumor_erv.obs_names,median_tumor):
            tumor_dict[erv] = value if value !=0 else 1e-5
        for erv,value in zip(normal_erv.index,mean_normal):
            normal_dict[erv] = value if value !=0 else 1e-5
        for erv,value in zip(tumor_erv.obs_names,n_sample_tumor):
            n_sample_dict[erv] = value
        col1 = list(set(tumor_dict.keys()).union(set(normal_dict.keys())))
        col2 = [tumor_dict.get(item,1e-5) for item in col1]
        col3 = [normal_dict.get(item,1e-5) for item in col1]
        col4 = np.log2(np.array(col2)/np.array(col3))
        col5 = [n_sample_dict.get(item,0) for item in col1]
        result = pd.DataFrame(data={'erv':col1,'median_tumor':col2,'mean_normal':col3,'logfc':col4,'n_sample':col5}).set_index(keys='erv')
        result = result.sort_values(by='n_sample',ascending=False)
        result.to_csv(os.path.join(atlas_path,'ERV.txt'),sep='\t')

    elif mode == 'plot':
        tumor_erv = ad.read_h5ad(os.path.join(atlas_path,'tumor_erv.h5ad'))
        tumor_erv_aux = pd.read_csv(os.path.join(atlas_path,'tumor_erv_aux_df.txt'),sep='\t',index_col=0)
        normal_erv = pd.read_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv.txt'),sep='\t',index_col=0)
        normal_erv_aux = pd.read_csv(os.path.join(HG38_NORMAL_DIR,'normal_erv_aux_df.txt'),sep='\t',index_col=0)

        # normalize using aux
        n_data = normal_erv.values / normal_erv_aux['total_count'].values.reshape(1,-1) * 1e6

        # reorg
        tumor_erv = ad.read_h5ad(os.path.join(atlas_path,'tumor_erv.h5ad'))
        tumor_erv.X = tumor_erv.layers['cpm']
        del tumor_erv.layers['cpm']
        normal_erv = pd.DataFrame(data=n_data,index=normal_erv.index,columns=normal_erv.columns)

        # use simple name
        tumor_erv.obs_names = [index.split(':')[0] for index in tumor_erv.obs_names]
        normal_erv.index = [index.split(':')[0] for index in normal_erv.index]

        for erv in ervs:
            tumor_expr = tumor_erv[erv,:].X.toarray().flatten()
            try:
                series = normal_erv.loc[erv,:]
            except:
                series = pd.Series(index=normal_erv.columns,data=np.full(normal_erv.shape[1],0),name='erv')
            normal_df = series.to_frame()

            fig = plt.figure()
            gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.1,0.9),wspace=0.2)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1],sharey=ax1)
            sns.boxplot(data=tumor_expr,color='red',ax=ax1)
            ax1.set_ylabel('value')
            ax1.set_xlabel('tumor')
            normal_df['tissue'] = [item.split(',')[1] for item in normal_df.index]
            normal_tissue_expr_dict = {t:sub_df[erv].values for t,sub_df in normal_df.groupby(by='tissue')}
            bp = ax2.boxplot(x=list(normal_tissue_expr_dict.values()),positions=np.arange(len(normal_tissue_expr_dict)),patch_artist=True)
            for flier in bp['fliers']:
                flier.set_markersize(1)
                flier.set_marker('o')
            for box in bp['boxes']:
                box.set_facecolor('green')
                box.set_edgecolor('black')
                box.set_linewidth(1)
            ax2.set_xticks(np.arange(len(normal_tissue_expr_dict)))
            ax2.set_xticklabels(list(normal_tissue_expr_dict.keys()),fontsize=10,rotation=90)
            ax2.set_xlabel('GTEx Normal')
            
            fig.suptitle('{}'.format(erv))
            plt.savefig(os.path.join(atlas_path,'{}_expr.{}'.format(erv.replace('/','_'),image_format)),bbox_inches='tight')
            plt.close()


    elif mode == 'fasta':
        erv = pd.read_csv(os.path.join(atlas_path,'ERV.txt'),sep='\t',index_col=0)
        cond = [True if ':' in item else False for item in erv.index]
        erv = erv.loc[cond,:]

        cond_erv = [True if item.split(':')[-1]=='LTR' else False for item in erv.index]
        cond_l1 = [True if (item.split(':')[-1]=='LINE') and (item.split(':')[-2]=='L1') else False for item in erv.index]
        cond_sva = [True if (item.split(':')[-1]=='Retroposon') and (item.split(':')[-2]=='SVA') else False for item in erv.index]
        cond_alu = [True if (item.split(':')[-1]=='SINE') and (item.split(':')[-2]=='Alu') else False for item in erv.index]
        cond_dna_transposon = [True if item.split(':')[-1]=='DNA' else False for item in erv.index]

        cond = np.any([cond_erv,cond_l1,cond_sva,cond_alu,cond_dna_transposon],axis=0).tolist()
        erv = erv.loc[cond,:]

        erv = erv.loc[(erv['logfc']>min_lfc_cutoff) & (erv['n_sample'] > total_number * frac_sample) & (erv['mean_normal']<mean_normal_cutoff) & (erv['median_tumor']>median_tumor_cutoff),:]
        erv.to_csv(os.path.join(atlas_path,'good_erv.txt'),sep='\t')
        print('total erv is {}'.format(erv.shape[0]))

        tid2cid = {}
        for item in erv.index:
            tid,gid,fid,cid = item.split(':')
            tid2cid[tid] = cid
        erv.index = [item.split(':')[0] for item in erv.index]
        erv2logfc = erv['logfc'].to_dict()
        erv2median_tumor = erv['median_tumor'].to_dict()
        erv2mean_normal = erv['mean_normal'].to_dict()
        erv2n_sample = erv['n_sample'].to_dict()
        erv2coord = pd.read_csv(TE_ANNO_HG38,sep='\t',index_col=0)['chromsome:start-stop:strand'].to_dict()  # 'LTR19-int_dup96': 'chrY:9995133-9995460:+'

        # hg38 dict_fa
        dict_fa = {}  
        with open(HG38_SEQ,'r') as in_handle:
            for title,seq in SimpleFastaParser(in_handle):
                dict_fa[title] = seq

        # load HERV splice
        db = pd.read_csv(HERV_splice_db,sep='\t',index_col=0)
        dup2coords = {}
        for item,herv_strand in zip(db.index,db['strand']):
            coords = [herv_strand]   # [-,chr1:777-888,chr1:888-899]
            for dup in item.split(','):
                coords.append(':'.join(erv2coord[dup].split(':')[:-1]))
            for dup in item.split(','):
                dup2coords[dup] = coords

        
        erv2orf = {}  
        for element in tqdm(erv.index,total=erv.shape[0]):
            final_pep_list = []
            coord = erv2coord[element]
            typ = tid2cid[element]
            chrom,coord,strand = coord.split(':')
            start,end = coord.split('-')
            dnas = []
            if chrom in set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']):
                # how you normally get dna
                start = int(start)
                end = int(end)
                dna = dict_fa[chrom][start-1:end]
                if strand == '-':
                    dna = str(Seq(dna).reverse_complement())

                # classify to consider special case
                if typ in ['SINE','LINE','DNA']:
                    dnas.append(dna)
                
                elif typ in ['LTR','Retroposon']:
                    
                    if element in dup2coords.keys():   # use spliced full-length to overwrite the old dna, special in special
                        splice_coords = dup2coords[element]
                        herv_strand = splice_coords[0]
                        herv_coords = splice_coords[1:]
                        dna = ''
                        for coord in herv_coords:
                            chrom,coord = coord.split(':')
                            start,end = coord.split('-')
                            start,end = int(start),int(end)
                            p_dna = dict_fa[chrom][start-1:end]
                            dna += p_dna
                        if herv_strand == '-':
                            dna = str(Seq(dna).reverse_complement())
                    

                    dnas.append(dna)
                    dna_anti = str(Seq(dna).reverse_complement())
                    dnas.append(dna_anti)

                
                # for each dna, translate from first, second and third position
                if len(dnas) == 1:
                    raw_pep_1 = Seq(dna).translate(to_stop=False)
                    raw_pep_2 = Seq(dna[1:]).translate(to_stop=False)
                    raw_pep_3 = Seq(dna[2:]).translate(to_stop=False)
                    for i,raw_pep in enumerate([raw_pep_1,raw_pep_2,raw_pep_3]):
                        short_pep_list = raw_pep.split('*')
                        short_pep_list = [(i+1,pep,'sense') for pep in short_pep_list if len(pep) >= MIN_PEP_LEN]
                        final_pep_list.extend(short_pep_list)
                elif len(dnas) == 2:
                    # first do sense:
                    dna = dnas[0]
                    raw_pep_1 = Seq(dna).translate(to_stop=False)
                    raw_pep_2 = Seq(dna[1:]).translate(to_stop=False)
                    raw_pep_3 = Seq(dna[2:]).translate(to_stop=False)
                    for i,raw_pep in enumerate([raw_pep_1,raw_pep_2,raw_pep_3]):
                        short_pep_list = raw_pep.split('*')
                        short_pep_list = [(i+1,pep,'sense') for pep in short_pep_list if len(pep) >= MIN_PEP_LEN]
                        final_pep_list.extend(short_pep_list)

                    # then do anti-sense
                    dna = dnas[1]
                    raw_pep_1 = Seq(dna).translate(to_stop=False)
                    raw_pep_2 = Seq(dna[1:]).translate(to_stop=False)
                    raw_pep_3 = Seq(dna[2:]).translate(to_stop=False)
                    for i,raw_pep in enumerate([raw_pep_1,raw_pep_2,raw_pep_3]):
                        short_pep_list = raw_pep.split('*')
                        short_pep_list = [(i+1,pep,'anti-sense') for pep in short_pep_list if len(pep) >= MIN_PEP_LEN]
                        final_pep_list.extend(short_pep_list)

            erv2orf[element] = final_pep_list

        with open(os.path.join(atlas_path,db_fasta_dir,'TE_self_translate.fasta'),'w') as f:
            for erv,orfs in erv2orf.items():
                coord = erv2coord[erv]
                for i,orf in enumerate(orfs):
                    phase,orf,sense = orf
                    f.write('>{}|{}|{}|{}|{}|{}|{}|{}|{}\n{}\n'.format(erv,coord,erv2median_tumor[erv],erv2mean_normal[erv],round(erv2logfc[erv],2),erv2n_sample[erv],str(i+1),phase,sense,orf))

            

    



def grab_intron_retention(star_tmp,gene_model_path):
    df = pd.read_csv(gene_model_path,sep='\t',header=None)
    df.columns = ['gene','chromsome','strand','exon_intron','start','end']
    df = df.loc[df['exon_intron'].str.startswith('I'),:]
    tbx_denovo = pysam.TabixFile('{}/denovo_sorted.gtf.gz'.format(star_tmp))
    tbx_guided = pysam.TabixFile('{}/guided_sorted.gtf.gz'.format(star_tmp))
    cond = []
    expr = []
    transcript = []
    exon_number_list = []
    pat_cov = re.compile(r'; cov "(.+?)";')
    pat_tid = re.compile(r'; transcript_id "(STRG.+?)";')
    pat_exon = re.compile(r'; exon_number "(\d+)";')
    pat_ref_gene = re.compile(r'; ref_gene_id "(.+?)";')
    counter_ambiguous = 0
    for row in tqdm(df.itertuples(),total=df.shape[0]):
        gene = row.gene
        start = row.start
        end = row.end
        chromsome = row.chromsome
        strand = row.strand
        flag = False
        cov = None
        tid = None
        exon_number = None

        # first fetch guided to see if this intron is overlapping with another genes' exon, make sure not because then the reads are ambiguous and result in false positive
        proceed = True
        for row in tbx_guided.fetch(chromsome,start-1,end):
            try:
                feature = row.split('\t')[2]
            except:  # means not overlapping with any features, good!
                break
            else:   
                if feature == 'transcript':
                    attr = row.split('\t')[-1]
                    try:
                        ref_gene = re.search(pat_ref_gene,attr).group(1)  # CHM13_G0000101
                    except:
                        # means it doesn't associated with known gene, it will be displayed as gene_id "STRG.3", good
                        continue
                    else:
                        if ref_gene not in gene: # CHM13_G0000101,gene_symbol, not in means this is another gene transcript, then no good
                            proceed = False    
                            counter_ambiguous += 1
                            break
            
        # based on whether to proceed  or not, we can now see if they can be overlapping with any denovo exons
        if proceed:
            for row in tbx_denovo.fetch(chromsome,start-1,end):
                try:
                    feature = row.split('\t')[2]
                except:
                    break
                else:
                    if feature == 'exon':
                        e_start = int(row.split('\t')[3])
                        e_end = int(row.split('\t')[4])
                        if e_start < start and e_end > end:
                            attribute = row.split('\t')[8]
                            cov = float(re.search(pat_cov,attribute).group(1))
                            tid = re.search(pat_tid,attribute).group(1)
                            exon_number = re.search(pat_exon,attribute).group(1)
                            flag = True
                            break
        # no matter proceed or not, you should append them
        cond.append(flag)
        expr.append(cov)
        transcript.append(tid)
        exon_number_list.append(exon_number)
    df['cond'] = cond
    df['expr'] = expr
    df['transcript'] = transcript
    df['exon_number'] = exon_number_list
    df = df.loc[df['cond'],:]
    print(star_tmp,counter_ambiguous)

    col_feature = []
    col_expr = []
    col_transcript = []
    col_exon_number = []
    for c1,c2,c3,c4,c5,c6,c7,c8,c9 in zip(df['gene'],df['chromsome'],df['strand'],df['exon_intron'],df['start'],df['end'],df['expr'],df['transcript'],df['exon_number']):
        col_feature.append(','.join([c1,c2,c3,c4,str(c5),str(c6)]))
        col_expr.append(c7)
        col_transcript.append(c8)
        col_exon_number.append(c9)
    result = pd.DataFrame(data={'feature':col_feature,'coverage':col_expr,'transcript':col_transcript,'exon_number':col_exon_number}).set_index(keys='feature')
    return result

'''ORF finder from SNAF'''
def score_GC(sequence):
    length_seq = len(sequence)
    counter = Counter(sequence)
    GC_percent = (counter.get('G',0) + counter.get('C',0)) / length_seq
    return GC_percent

def score_coding_potential(sequence):
    # coding frequency table is from GenScript webpage
    usage_dict = {'TTT':16.9,'TTC':20.4,'TTA':7.2,'TTG':12.6,'TAT':12.0,'TAC':15.6,'TAA':0.7,'TAG':0.5,
                  'CTT':12.8,'CTC':19.4,'CTA':6.9,'CTG':40.3,'CAT':10.4,'CAC':14.9,'CAA':11.8,'CAG':34.6,
                  'ATT':15.7,'ATC':21.4,'ATA':7.1,'ATG':22.3,'AAT':16.7,'AAC':19.5,'AAA':24.0,'AAG':32.9,
                  'GTT':10.9,'GTC':14.6,'GTA':7.0,'GTG':28.9,'GAT':22.3,'GAC':26.0,'GAA':29.0,'GAG':40.8,
                  'TCT':14.6,'TCC':17.4,'TCA':11.7,'TCG':4.5,'TGT':9.9,'TGC':12.2,'TGA':1.3,'TGG':12.8,
                  'CCT':17.3,'CCC':20.0,'CCA':16.7,'CCG':7.0,'CGT':4.7,'CGC':10.9,'CGA':6.3,'CGG':11.9,
                  'ACT':12.8,'ACC':19.2,'ACA':14.8,'ACG':6.2,'AGT':11.9,'AGC':19.4,'AGA':11.5,'AGG':11.4,
                  'GCT':18.6,'GCC':28.5,'GCA':16.0,'GCG':7.6,'GGT':10.8,'GGC':22.8,'GGA':16.3,'GGG':16.4} 
    # do a normaliztion for each triplet, then for all the triplet's sum, divided by the number of triplet
    min_freq = min(list(usage_dict.values()))
    max_freq = max(list(usage_dict.values()))
    norm_usage_dict = {k:(v-min_freq)/(max_freq-min_freq) for k,v in usage_dict.items()}      
    length_seq = len(sequence)
    num_triplet = length_seq / 3
    i = 0   
    score = 0
    while i + 2 < length_seq:
        triplet = sequence[i:i+3:1]
        score_tri = norm_usage_dict[triplet]
        score += score_tri
        i += 3
    score_potential = score/num_triplet # scale by the number of triplet in the sequence
    return score_potential

def orf2pep(orf):
    assert len(orf) % 3 == 0
    pep = str(Seq(orf).translate(to_stop=False))
    assert '*' not in pep
    return pep

def transcript2orf(cdna):
    candidate_orfs = []
    p_start = re.compile(r'ATG')
    p_end = re.compile(r'(TAA|TGA|TAG)')
    ms = list(re.finditer(p_start,cdna))   
    if len(ms) > 0:
        for s in ms:
            s_pos = int(s.start())
            me = list(re.finditer(p_end,cdna))
            if len(me) > 0:
                for e in me:
                    e_pos = int(e.start())
                    if s_pos < e_pos and s_pos % 3 == e_pos % 3:
                        orf = cdna[s_pos:e_pos]
                        valid = True
                        mo = list(re.finditer(p_end,orf))   # still need to check whether there is stop codon in the orf
                        if len(mo) > 0:
                            for o in mo:
                                o_pos = int(o.start())
                                if o_pos % 3 == 0:
                                    valid = not valid
                                    break
                        if valid:
                            candidate_orfs.append(orf)
            else:  # maybe the stop codon not in the cdna, the last reading codon is the last three bases
                orf = cdna[s_pos:]
                valid = True
                if len(orf) % 3 != 0:
                    valid = not valid
                mo = list(re.finditer(p_end,orf))
                if len(mo) > 0:
                    for o in mo:
                        o_pos = int(o.start())
                        if o_pos % 3 == 0:
                            valid = not valid
                            break
                if valid:
                    candidate_orfs.append(orf)
    return candidate_orfs



def prioritize_orf(candidate_orfs,min_len=30*3,tol_len=8*3):
    max_orf = ''
    max_length = 0
    max_score = 0
    for orf in candidate_orfs:
        if len(orf) < min_len:
            continue
        score = score_GC(orf) + score_coding_potential(orf)
        if len(orf) - max_length > tol_len:
            max_length = len(orf)
            max_score = score
            max_orf = orf
        else:
            if len(orf) > max_length and score > max_score:
                max_length = len(orf)
                max_score = score
                max_orf = orf  
    return max_orf         

def run_orffinder(cdna,method):
    if method == 'snaf':
        candidate_orfs = transcript2orf(cdna)
        max_orf = prioritize_orf(candidate_orfs)
    elif method == 'ncbi':
        # first need to download https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/
        # then conda install libuv: conda install conda-forge::libuv, it is in ./lib
        # add to LD_LIBRARY_PATH: export LD_LIBRARY_PATH=/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/lib:$LD_LIBRARY_PATH
        # ./ORFfinder -in test.fasta -s 0 -n true -strand plus -out test.out.fasta -outfmt 1
        # it returns cds including stop codon
        input_file_path = os.path.join(OUTDIR,'tmp.ncbi.orffinder.input.fasta')
        output_file_path = os.path.join(OUTDIR,'tmp.ncbi.orffinder.output.fasta')
        with open(input_file_path,'w') as f:
            f.write('>{}\n{}\n'.format('cdna',cdna))
        subprocess.run("{} -in {} -s 0 -n true -strand plus -out {} -outfmt 1".format(NCBI_ORFFINDER_PATH,input_file_path,output_file_path),shell=True)
        with open(output_file_path,'r') as in_handle:
            candidate_orfs = [seq[:-3] for title,seq in SimpleFastaParser(in_handle)]
            max_orf = prioritize_orf(candidate_orfs)
        subprocess.run("rm {}".format(input_file_path),shell=True)
        subprocess.run("rm {}".format(output_file_path),shell=True)
    return max_orf
    


def integrate_intron_retention(mode,normal_recurrency=0,intron_path=None,total_number=None,intron_peptide_path=None,plot_intron=None):

    if mode == 'compute':
        snaf_db = pd.read_csv(SPLICING_SNAF,sep='\t',index_col=0)
        snaf_db['strand'] = [item.split('(')[1].rstrip(')') for item in snaf_db['coordinate']]
        uid2strand = snaf_db['strand'].to_dict()

        blacklist = pd.read_csv(SPLICING_BLACKLIST,sep='\t')
        blacklist.columns = ['chrom','start','end','uid']
        blacklist_start_pos = []   # remember, for reverse strand, this start means forward strand start
        blacklist_end_pos = []
        pat1 = r'^ENSG\d+:E.+?-I.+?$'
        pat2 = r'^ENSG\d+:I.+?-E.+?$'
        for row in blacklist.itertuples():
            actual_start = row.start + 1
            actual_end = row.end
            if actual_end == actual_start + 1:   # only care about intron retention
                strand = uid2strand[row.uid]
                chrom = row.chrom
                if re.search(pat1,row.uid) and strand == '+' or (re.search(pat2,row.uid) and strand == '-'):   # E1.4-I1.1
                    blacklist_start_pos.append('{}:{}'.format(chrom,str(actual_end)))
                elif (re.search(pat2,row.uid) and strand == '+') or (re.search(pat1,row.uid) and strand == '-'): # I1.1-E2.1
                    blacklist_end_pos.append('{}:{}'.format(chrom,str(actual_start)))
        blacklist_start_pos = set(blacklist_start_pos)
        blacklist_end_pos = set(blacklist_end_pos)

        # readin intron
        intron = pd.read_csv(intron_path,sep='\t')
        data = []
        for event,sub_df in tqdm(intron.groupby(by='feature')):
            data.append((event,sub_df.shape[0],round(sub_df['coverage'].values.mean(),2)))
        intron_rec = pd.DataFrame.from_records(data,columns=['intron','n_sample','mean_coverage'])

        # whether in snaf normal db
        cond = []
        for event in intron_rec['intron']:
            gid,gs,chrom,strand,intron,start,end = event.split(',')
            if 'MSTRG' in gs or gs.startswith('AC') or gs.startswith('AL') or gs.startswith('AP') or '-AS' in gs:
                cond.append(False)
                continue
            start = '{}:{}'.format(chrom,start)
            end = '{}:{}'.format(chrom,end)
            if (not start in blacklist_start_pos) and (not end in blacklist_end_pos):
                cond.append(True)
            else:
                cond.append(False)
        intron_rec['cond'] = cond


        # intron blacklist
        # meta = pd.read_csv(HG38_NORMAL_META,sep='\t',index_col=0,header=None)
        # meta.columns = ['srr']
        # tissue2srr = meta['srr'].to_dict()
        # detected_introns = []
        # for t,ss in tissue2srr.items():
        #     for s in ss.split(','):
        #         intron_result_path = os.path.join(HG38_NORMAL_DIR,s,'intron.txt')
        #         intron_df = pd.read_csv(intron_result_path,sep='\t',index_col=0)
        #         detected_introns.extend([(intron,t,s,cov) for intron,cov in zip(intron_df.index,intron_df['coverage'])])
        # normal_intron = pd.DataFrame.from_records(data=detected_introns,columns=['intron','tissue','sample','coverage'])
        # normal_intron.to_csv(NORMAL_INTRON_PATH,sep='\t',index=None)
        
        normal_intron = pd.read_csv(NORMAL_INTRON_PATH,sep='\t')
        normal_intron_dict = {}
        for intron,sub_df in tqdm(normal_intron.groupby(by='intron')):
            normal_intron_dict[intron] = {}
            for tissue,sub_df2 in sub_df.groupby(by='tissue'):
                normal_intron_dict[intron][tissue] = ','.join([str(item) for item in sub_df2['coverage'].round(2)])
        intron_rec['normal'] = intron_rec['intron'].map(normal_intron_dict)
        intron_rec.to_csv(os.path.join(atlas_path,'intron_rec.txt'),sep='\t',index=None)

    elif mode == 'plot':
        # tumor
        intron = pd.read_csv(intron_path,sep='\t')
        tumor_expr = intron.loc[intron['feature']==plot_intron,'coverage'].values
        # normal
        normal_intron = pd.read_csv(NORMAL_INTRON_PATH,sep='\t')
        all_tissues = normal_intron['tissue'].unique()
        subset = normal_intron.loc[normal_intron['intron']==plot_intron,:]
        normal_expr_dict = {}
        for tissue,sub_df in subset.groupby(by='tissue'):
            normal_expr_dict[tissue] = sub_df['coverage'].values
        # fill up
        missing_tissues = set(all_tissues).difference(set(normal_expr_dict.keys()))
        missing_tissues_expr_dict = {t:np.full(5,fill_value=0) for t in missing_tissues}
        normal_expr_dict.update(missing_tissues_expr_dict)

        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.1,0.9),wspace=0.2)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1],sharey=ax1)
        sns.boxplot(data=tumor_expr,color='red',ax=ax1)
        ax1.set_ylabel('value')
        ax1.set_xlabel('tumor')

        bp = ax2.boxplot(x=list(normal_expr_dict.values()),positions=np.arange(len(normal_expr_dict)),patch_artist=True)
        for flier in bp['fliers']:
            flier.set_markersize(1)
            flier.set_marker('o')
        for box in bp['boxes']:
            box.set_facecolor('green')
            box.set_edgecolor('black')
            box.set_linewidth(1)
        ax2.set_xticks(np.arange(len(normal_expr_dict)))
        ax2.set_xticklabels(list(normal_expr_dict.keys()),fontsize=0.5,rotation=90)
        ax2.set_xlabel('GTEx Normal')
        
        fig.suptitle('{}'.format(plot_intron))
        plt.savefig(os.path.join(atlas_path,'{}_expr.{}'.format(plot_intron,image_format)))
        plt.close() 



    elif mode == 'fasta':
        final = pd.read_csv(intron_path,sep='\t',index_col=0)
        final = final.loc[final['cond'],:]
        final['normal'] = [literal_eval(item) if isinstance(item,str) else {} for item in final['normal']]
        final['normal_recurrency'] = [len(item) for item in final['normal']]
        final = final.loc[final['normal_recurrency'] <= normal_recurrency, :]
        final = final.loc[final['n_sample'] > 0.15*total_number,:]

        # get intron2peptide
        final_peptide = pd.read_csv(os.path.join(intron_peptide_path),sep='\t')
        intron2peptide = {}
        for uid,sub_df in tqdm(final_peptide.groupby(by='uid')):
            intron2peptide.setdefault(uid,[]).extend(list(set(sub_df['intron_pepseq'].values.tolist())))

        # start to write
        with open(os.path.join(atlas_path,db_fasta_dir,'intron.fasta'),'w') as f:
            for row in final.itertuples():
                peps = intron2peptide.get(row.Index,None)
                if peps is not None:
                    for i,pep in enumerate(peps):
                        f.write('>{}|{}|{}|{}|intron\n{}\n'.format(row.Index,row.n_sample,str(round(row.mean_coverage,2)),str(i+1),pep))




def integrate_circ_rna(mode):
    if mode == 'compute':
        tumor_circ_list = []
        for s in tqdm(samples):
            circ_te_tmp = os.path.join(ROOTDIR,s,'circ_te_tmp')
            circ_df = pd.read_csv(os.path.join(circ_te_tmp,'circularRNA_known.txt'),sep='\t',header=None)
            circ_df = circ_df.iloc[:,[0,1,2,3,5]]
            circ_df.columns = ['chrom','start','end','name_read','strand']
            # grab reads
            col = []
            for item in circ_df['name_read']:
                _, read = item.split('/')
                col.append(read)
            circ_df['count'] = col
            circ_df.drop(columns=['name_read'],inplace=True)
            tumor_circ_list.append(circ_df)
        tumor_circ = pd.concat(tumor_circ_list,axis=0,keys=samples).reset_index(level=-2)
        tumor_circ.rename(columns={'level_0':'sample'},inplace=True)
        tumor_circ['uid'] = [item1 + ':' + str(item2) + '-' + str(item3) for item1,item2,item3 in zip(tumor_circ['chrom'],tumor_circ['start'],tumor_circ['end'])]

        normal_circ = pd.read_csv(CIRC_RNA_NORMAL,sep=',',header=None)
        blacklist = set(normal_circ.iloc[:,4].tolist())
        tumor_circ['cond'] = ~tumor_circ['uid'].isin(blacklist)
        tumor_circ.to_csv(os.path.join(OUTDIR,'circ_rna.txt'),sep='\t',index=None)

    elif mode == 'fasta':
        tumor_circ = pd.read_csv(os.path.join(OUTDIR,'circ_rna.txt'),sep='\t')
        # filter and fix the +1 issue (https://github.com/YangLab/CIRCexplorer2/issues/75), it seems that if the splicing site is 3', it doesn't need to +1
        tumor_circ = tumor_circ.loc[tumor_circ['cond'],:]
        tumor_circ['start'] = [int(item) + 1 for item in tumor_circ['start']] 
        tumor_circ['uid'] = [item1 + ':' + str(item2) + '-' + str(item3) for item1,item2,item3 in zip(tumor_circ['chrom'],tumor_circ['start'],tumor_circ['end'])]
        dict_fa = {}  
        with open(HG38_SEQ,'r') as in_handle:
            for title,seq in SimpleFastaParser(in_handle):
                dict_fa[title] = seq
        junction2peptide = {}
        for item,strand in tqdm(zip(tumor_circ['uid'],tumor_circ['strand']),total=tumor_circ.shape[0]):
            info = []
            chrom,coord = item.split(':')
            start,end = coord.split('-')
            start,end = int(start),int(end)
            overhang = 50
            # below is the unique to circ, because it's the reverse regarding what is first section, what is second section
            # see chr1:400-500, strand +, first section is 450-500, second section is 400-450
            if strand == '+':
                first_section = dict_fa[chrom][end-overhang:end]
                second_section = dict_fa[chrom][start-1:start-1+overhang]
            elif strand == '-':
                old_first_section = dict_fa[chrom][start-1:start-1+overhang]
                old_second_section = dict_fa[chrom][end-overhang:end]
                # reverse complement process
                first_section = Seq(old_second_section).reverse_complement()
                second_section = Seq(old_first_section).reverse_complement()
            # consider phase = 0
            defacto_first = first_section[-(3*(MAX_PEP_LEN-1)):]
            defacto_second = second_section[:(3*(MAX_PEP_LEN-1))]
            defacto_junction = defacto_first + defacto_second
            defacto_peptide = str(Seq(defacto_junction).translate(to_stop=False))
            if '*' not in defacto_peptide:
                info.append((strand,0,defacto_peptide))
            # consider phase = 1
            defacto_first = first_section[-(3*(MAX_PEP_LEN-1)+2):]
            defacto_second = second_section[:(3*(MAX_PEP_LEN-1)+1)]
            defacto_junction = defacto_first + defacto_second
            defacto_peptide = str(Seq(defacto_junction).translate(to_stop=False))
            if '*' not in defacto_peptide:
                info.append((strand,1,defacto_peptide)) 
            # consider phase = 2
            defacto_first = first_section[-(3*(MAX_PEP_LEN-1)+1):]
            defacto_second = second_section[:(3*(MAX_PEP_LEN-1)+2)]
            defacto_junction = defacto_first + defacto_second
            defacto_peptide = str(Seq(defacto_junction).translate(to_stop=False))
            if '*' not in defacto_peptide:
                info.append((strand,2,defacto_peptide)) 
            junction2peptide[item] = info   
        # distribute
        for sample,sub_df in tumor_circ.groupby(by='sample'):
            with open(os.path.join(OUTDIR,db_fasta_dir,sample,'tmp_circ.fasta'),'w') as f:
                for junction,count in zip(sub_df['uid'],sub_df['count']):
                    info = junction2peptide[junction]
                    for item in info:
                        strand = item[0]
                        phase = item[1]
                        seq = item[2]
                        f.write('>{}|{}|{}|{}|{}\n{}\n'.format(junction,strand,str(phase),count,'circRNA',seq))




def integrate_nuorf(mode,modality=None,class2=False,obj=None,t_dic=None,n_dic=None,immuno_dir=None,technology=None,col=None,metadata_path=None,final_path=None):
    if mode == 'fasta':
        if class2:
            with open(NUORF_FASTA,'r') as in_handle, open(os.path.join(atlas_path,db_fasta_dir,'nuorf.fasta'),'w') as out_handle:
                for title,seq in SimpleFastaParser(in_handle):
                    if len(seq) >= MIN_PEP_LEN:
                        out_handle.write('>{}\n{}\n'.format(title,seq))
        else:
            subprocess.run(['cp',NUORF_FASTA,os.path.join(atlas_path,db_fasta_dir,'nuorf.fasta')])
    elif mode == 'plot':
        if modality == 'qualitative':
            normal_hits = t_dic # {tissue:value}
            normal_hits.update(n_dic)
            array_data = np.array(list(normal_hits.values())).reshape((1,-1))
            tissues = list(normal_hits.keys())
            fig,ax = plt.subplots()
            cax = ax.imshow(array_data,cmap='bwr')
            ax.set_xticks(np.arange(len(normal_hits)))
            ax.set_xticks(np.arange(len(normal_hits))+0.5,minor=True)
            ax.set_xticklabels(tissues,fontsize=1,rotation=60)
            ax.grid(which='minor', axis='x', color='black', linestyle='-', linewidth=0.5)
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
            ax.tick_params(axis='x', which='minor', bottom=False, labelbottom=False)
            fig.colorbar(cax, ax=ax)
            plt.savefig(os.path.join(atlas_path,'nuorf_{}_{}.{}'.format(modality,obj,image_format)),bbox_inches='tight')
            plt.close()
        elif modality == 'quantitative':
            normal = pd.read_csv(HLA_LIGAND_ATLAS_INTENSITY,sep='\t')
            peptide = obj
            if technology == 'bruker':
                file_name = 'accumulatedMsmsScans.txt'
            elif technology == 'orbitrap':
                file_name = 'msmsScans_new.txt'
            all_txts = subprocess.run('find {} -type f -name "{}"'.format(immuno_dir,file_name),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
            tumor_dfs = []
            tumor_raws = []
            tumor_raw_dict = {}
            # read in metadata
            metadata = pd.read_csv(metadata_path,sep='\t')
            raw2bio = {}
            for row in metadata.itertuples():
                raw2bio[row.sample.split('.raw')[0]] = row.biology
            # read in final
            hla_dict = literal_eval(pd.read_csv(final_path,sep='\t',index_col=0).at[peptide,'presented_by_each_sample_hla'])
            # now process
            for txt in all_txts:
                msms = pd.read_csv(txt,sep='\t')
                msms = msms.loc[msms['Reverse'] != '+',:]
                msms = msms.loc[msms['Identified']=='+',:]
                for raw,sub_df in msms.groupby(by='Raw file'):
                    each_raw_data = []
                    sub_df['Precursor intensity'] = sub_df['Precursor intensity'].fillna(value=0)
                    sub_df = sub_df.sort_values(by='Precursor intensity',ascending=True)
                    sub_df['percentile'] = [(i+1)/sub_df.shape[0] for i in range(sub_df.shape[0])]
                    for p,sub_df2 in sub_df.groupby(by='Sequence'):
                        intensity = sub_df2['Precursor intensity'].values.max()
                        percentile = sub_df2['percentile'].values.max()
                        each_raw_data.append((p,intensity,percentile,raw2bio[raw]))
                    each_raw_df = pd.DataFrame.from_records(data=each_raw_data,columns=['peptide','intensity','percentile','bio'])
                    each_raw_df = each_raw_df.loc[each_raw_df['intensity']>0,:]
                    tumor_raw_dict[raw] = each_raw_df
                    upper = np.quantile(each_raw_df['intensity'].values,0.75)
                    each_raw_df['norm'] = np.log2(each_raw_df['intensity'].values/upper)
                    each_raw_df['norm'] = np.log2(each_raw_df['intensity'].values/upper)
                    tumor_dfs.append(each_raw_df)
                    tumor_raws.append(raw)
            final = pd.concat(tumor_dfs,axis=0,keys=tumor_raws).reset_index(level=-2).rename(columns={'level_0':'raw_file'})
            final['log_intensity'] = np.log2(final['intensity'].values)

            # you can first draw the rank abundance plot
            final_pep = final.loc[final['peptide']==peptide,:]
            n = final_pep.shape[0]
            n_col = 3
            n_row = n // n_col + 1
            fig,axes = plt.subplots(nrows=n_row,ncols=n_col,gridspec_kw={'wspace':0.5,'hspace':0.5})
            axes = axes.flatten()
            for i,row in enumerate(final_pep.itertuples()):
                df = tumor_raw_dict[row.raw_file]
                df['log_intensity'] = np.log2(df['intensity'].values)
                df = df.sort_values(by='log_intensity')
                x = np.arange(df.shape[0])
                y = df['log_intensity'].values.tolist()
                s = [10 if item == peptide else 0.5 for item in df['peptide']]
                c = ['r' if item == peptide else 'k' for item in df['peptide']]
                axes[i].scatter(x,y,s,c)
            plt.savefig(os.path.join(atlas_path,'{}_rank_abundance.{}'.format(peptide,image_format)))
            plt.close()

            normal_normalized_all_values = {}
            for tissue,sub_df in normal.loc[normal['peptide']==peptide,:].groupby(by='tissue'):
                normal_normalized_all_values[tissue] = sub_df[col].values

            fig = plt.figure()
            gs = mpl.gridspec.GridSpec(nrows=1,ncols=2,width_ratios=(0.1,0.9),wspace=0.2)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1],sharey=ax1)  

            peptide_tumor = final.loc[final['peptide']==peptide,:]


            col1 = []
            for item in peptide_tumor['bio']:
                hla_list = [element[0].split('HLA-')[1] if element[0].startswith('HLA-') else 'None' for element in hla_dict[item]]
                hla_list = list(set(hla_list))
                hla_list.sort()
                col1.append(';'.join(hla_list))
            peptide_tumor['hla'] = col1


            unique_hla = peptide_tumor['hla'].unique()
            from colors import pick_n_colors
            select_colors = pick_n_colors(len(unique_hla))
            for i,hla in enumerate(unique_hla):
                tumor_expr = peptide_tumor.loc[peptide_tumor['hla']==hla,col].values
                sns.swarmplot(data=tumor_expr,color=select_colors[i],ax=ax1)
                import matplotlib.lines as mlines
            ax2.legend(handles=[mlines.Line2D([],[],marker='o',linestyle='',color=i) for i in select_colors],labels=unique_hla.tolist(),bbox_to_anchor=(1,1),loc='upper left') 
            ax1.set_ylabel('value')
            ax1.set_xlabel('tumor')
            if len(normal_normalized_all_values) > 0:
                bp = ax2.boxplot(x=list(normal_normalized_all_values.values()),positions=np.arange(len(normal_normalized_all_values)),patch_artist=True)
                for flier in bp['fliers']:
                    flier.set_markersize(1)
                    flier.set_marker('o')
                for box in bp['boxes']:
                    box.set_facecolor('green')
                    box.set_edgecolor('black')
                    box.set_linewidth(1)
                ax2.set_xticks(np.arange(len(normal_normalized_all_values)))
                ax2.set_xticklabels(list(normal_normalized_all_values.keys()),fontsize=10,rotation=90)
            fig.suptitle(peptide)
            plt.savefig(os.path.join(atlas_path,'{}_{}.{}'.format(peptide,col,image_format)),bbox_inches='tight')
            plt.close()




if __name__ == '__main__':

    '''
    This script is an adaptation for pan-cancer analysis to reuse some of the function from NeoVerse
    '''

    parser = argparse.ArgumentParser(description='NeoVerse-py')
    parser.add_argument('--config',type=str,default='',help='path to the config json')
    parser.add_argument('--running_mode',type=str,default=None,help='non_interactive or interactive')
    parser.add_argument('--sub_mode',type=str,default=None,help='if interactive, which sub_mode')
    args = parser.parse_args()

    config_path = args.config
    running_mode = args.running_mode
    sub_mode = args.sub_mode

    with open(config_path,'r') as f:
        config = json.load(f)['python_script']
    tunable = config['tunable_parameters']
    fixed = config['fixed_parameters']


    MIN_PEP_LEN = tunable['min_pep_len']
    MAX_PEP_LEN = tunable['max_pep_len']
    image_format = tunable['image_format']
    MAX_PEP_LEN_SPLICING = MAX_PEP_LEN
    MIN_PEP_LEN_SPLICING = MIN_PEP_LEN
    atlas_path = tunable['atlas_path']
    GTEX_GENE = fixed['gtex_median']
    GTEX_GENE_ALL = fixed['gtex_all']
    GTEX_META = fixed['gtex_meta']
    GTEX_GENE_ALL_H5AD = fixed['gtex_all_h5ad']
    GTF_GENE = fixed['hg38_gtf']
    MEMBRANE_GENE = fixed['membrane_gene']
    NUORF = fixed['nuorf']
    PROTEIN = fixed['protein']
    # I manually switched MUC19 based on original uniprot file as there's stop codon, also GPATCH4 canonical one
    ENSEMBL_PROTEIN_PATH = fixed['ensembl_protein_path']


    VARIANT_ENSEMBL_GTF = fixed['variant_ensembl_gtf']

    SPLICING_SNAF = fixed['splicing_snaf']
    SPLICING_BLACKLIST = fixed['splicing_blacklist']  # in-house version is postlift.bed
    SPLICING_GTEX_HG38 = fixed['splicing_gtex_hg38']
    SPLICING_TCGA_HG38 = fixed['splicing_tcga_hg38']
    SPLICING_GTEX_MAPPING_HG38 = fixed['splicing_gtex_mapping_hg38']
    SPLICING_TCGA_MAPPING_HG38 = fixed['splicing_tcga_mapping_hg38']

    REF_PROTEOME = fixed['ref_proteome']
    OTHER_PROTEOME = fixed['other_proteome']

    TE_GTF_HG38 = fixed['te_gtf_hg38']  # in-house version is t2t
    TE_ANNO_HG38 = fixed['te_anno_hg38']
    HG38_NORMAL_DIR = fixed['normal_hg38_dir']
    HG38_NORMAL_META = fixed['normal_hg38_meta']
    T2T_NORMAL_META = fixed['normal_t2t_meta']
    T2T_NORMAL_DIR = fixed['normal_t2t_dir']
    ERV_ANNO = fixed['geve_hg38_anno']
    ERV_SEQ = fixed['geve_hg38_seq']
    ERV_BED = fixed['geve_hg38_bed']
    NORMAL_PATHOGEN_PATH = fixed['combined_normal_pathogen_t2t']
    NORMAL_INTRON_PATH = fixed['combined_normal_intron_hg38']   # differ



    HG38_SEQ = fixed['hg38_fa']

    NUORF_FASTA = fixed['nuorf_fa']
    IEATLAS_NORMAL = fixed['ieatlas_normal']
    IEATLAS_CANCER = fixed['ieatlas_cancer']
    junction2peptide_db = fixed['junction2isoformDB_path']
    HERV_splice_db = fixed['HERV_splice_db_path']
    manifest_path = tunable['manifest_path']

    db_fasta_dir = tunable['db_fasta_dir']
    total = tunable['total']
    class2 = tunable['class2']
    te_dir_path = tunable['te_dir_path']

    HLA_LIGAND_ATLAS_INTENSITY = fixed['hla_ligand_atlas_intensity']



    if running_mode == 'non_interactive':
        integrate_gene(mode='compute',gene_path=os.path.join(atlas_path,'gene_tpm.txt'))
        integrate_splicing(mode='compute',splicing_path=os.path.join(atlas_path,'splicing_all.txt'))
        integrate_splicing(mode='fasta',splicing_path=os.path.join(atlas_path,'splicing_rec.txt'),total_sample=total)
        integrate_intron_retention(mode='compute',intron_path=os.path.join(atlas_path,'intron_all.txt'))
        integrate_intron_retention(mode='fasta',intron_path=os.path.join(atlas_path,'intron_rec.txt'),total_number=total,intron_peptide_path=os.path.join(atlas_path,'intron_peptide_all.txt'))
        integrate_variant(mode='fasta',mutation_rec_path=os.path.join(atlas_path,'mutation_rec.txt'),protein_path=ENSEMBL_PROTEIN_PATH,fs_dict={})
        integrate_nuorf(mode='fasta',class2=class2)
        integrate_pathogen(mode='compute',pathogen_path=os.path.join(atlas_path,'pathogen_all.txt'),total_sample=total)
        integrate_erv(mode='compute',te_dir_path=te_dir_path)
        integrate_erv(mode='fasta',total_number=total)


    elif running_mode == 'interactive':

        with open(config_path,'r') as f:
            config = json.load(f)['python_interact']
        
        image_format = config['image_format']

        if sub_mode == 'plot_gene':
            for symbol,ensg in config['gene']['markers'].items():
                integrate_gene(mode='plot',ensg=ensg,symbol=symbol,plot_type='boxplot+boxplot',cat_dic=config['gene']['cat_dic'],gene_path=os.path.join(atlas_path,'gene_tpm.txt'))

        elif sub_mode == 'plot_pathogen':

            for strain in config['pathogen']['plot_strains']:
                integrate_pathogen(mode='plot',plot_strain=strain)
        
        elif sub_mode == 'fasta_pathogen':
            integrate_pathogen(mode='fasta',included_strains=config['pathogen']['included_strains']) 

        elif sub_mode == 'fs_variant':
            integrate_variant(mode='fasta',mutation_rec_path=os.path.join(atlas_path,'mutation_rec.txt'),protein_path=ENSEMBL_PROTEIN_PATH,fs_dict=config['variant']['fs_dict'])

        elif sub_mode == 'plot_splicing':
            events = list(config['splicing']['events'])
            integrate_splicing(mode='plot',events=config['splicing']['events'],splicing_path=os.path.join(atlas_path,'splicing_all.txt'))

        elif sub_mode == 'plot_erv':
            ervs = config['erv']['markers']
            integrate_erv(mode='plot',ervs=ervs)

        elif sub_mode == 'plot_ir':
            for ir in config['intron_retention']['introns']:
                integrate_intron_retention(mode='plot',plot_intron=ir)
        
        elif sub_mode == 'plot_peptide_quali':
            integrate_nuorf(mode='plot',modality='qualitative',obj=config['nuorf']['obj'],t_dic=config['nuorf']['t_dic'],n_dic=config['nuorf']['n_dic'])

        elif sub_mode == 'plot_peptide_quanti':
            integrate_nuorf(mode='plot',modality='quantitative',obj=config['nuorf']['obj'],immuno_dir=config['nuorf']['immuno_dir'],technology=config['nuorf']['technology'],col=config['nuorf']['col'],metadata_path=config['nuorf']['metadata_path'],final_path=config['nuorf']['final_path'])













 
