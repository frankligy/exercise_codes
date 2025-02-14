#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
import subprocess
from tqdm import tqdm
import json
import multiprocessing as mp
from io import StringIO
import bisect
import argparse

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


def classify_source(item):
    if item.startswith('REV_'):
        identity = 'reverse'
    elif item.startswith('CON_'):
        identity = 'contaminant'
    elif 'missense_variant' in item or 'inframe' in item or 'frameshift_variant' in item:
        identity = 'variant'
    elif item.startswith('tr') or item.startswith('sp'):
        identity = 'pathogen'
    elif 'fusion' in item:
        identity = 'fusion' 
    elif 'intron' in item:
        identity = 'intron_retention'
    elif 'nuORF' in item:
        identity = 'nuORF'
    elif item.startswith('chr'):
        if 'TE_info' in item:
            identity = 'TE_chimeric_transcript'
        else:
            identity = 'splicing'
    elif '_dup' in item or item.endswith('sense'):
        identity = 'ERV'
    elif item.startswith('ENSG'):
        identity = 'self_gene'
    elif item.startswith('nc'):
        identity = 'nc_isoform'
    else:
        identity = 'unknown'
    return identity

def rederive_fdr(msms,fdr):
    df = msms
    df = df.sort_values(by='PEP',ascending=True)
    cond = [1 if isinstance(item,str) else 0 for item in df['Reverse']]
    total_decoy = np.cumsum(cond)
    prop_decoy = total_decoy / (np.arange(len(total_decoy)) + 1)
    col1 = np.where(prop_decoy<fdr,'+',None)
    col2 = [None if item2 == '+' else item1 for item1,item2 in zip(col1,df['Reverse'])]
    df['Identified'] = col2
    return df

def get_stats_spectra():
    # how many submitted and identified (no FDR) MS2
    data = []
    for k,v in result_dict.items():
        technology = technology_dict[k]
        if technology == 'orbitrap':
            file_name = 'msmsScans_new.txt'
        elif technology == 'bruker':
            file_name = 'accumulatedMsmsScans_new.txt'
        msms_path = os.path.join(v,'combined','txt',file_name)
        msms = pd.read_csv(msms_path,sep='\t')
        n_submitted = msms.shape[0]
        msms = msms.loc[msms['Reverse'] != '+',:]   # no matter how, get rid of reverse map first
        msms_identified = msms.loc[msms['Identified']=='+',:]
        n_identified = msms_identified.shape[0]
        msms_decent = msms.loc[msms['Score']>=40,:]
        n_decent = msms_decent.shape[0]
        if technology == 'orbitrap':
            msms_all = msms.loc[msms['Sequence'].notna(),:]
            n_all = msms_all.shape[0]
        elif technology == 'bruker':
            cond = [True if len(item) > 1 else False for item in msms['Sequence']]
            msms_all = msms.loc[cond,:]
            n_all = msms_all.shape[0]
        data.append((k,n_submitted,n_all,n_decent,n_identified))
    if technology == 'orbitrap':
        df = pd.DataFrame.from_records(data=data,columns=['sample','Total_MS2_spectra','n_all','n_decent(score>40)','n_identified (FDR)']).set_index(keys='sample')
    elif technology == 'bruker':
        df = pd.DataFrame.from_records(data=data,columns=['sample','Total_MS2_spectra','n_all','n_decent(score>40)','n_identified (FDR)']).set_index(keys='sample')
    df.sort_values(by='Total_MS2_spectra',ascending=False,inplace=True)
    df.to_csv(os.path.join(outdir,'spectra.txt'),sep='\t')
    fig,ax = plt.subplots()
    ax = df.plot.bar()
    plt.savefig(os.path.join(outdir,'submitted_ms2.pdf'),bbox_inches='tight')
    plt.close()

def get_stats_source():
    # summarize source
    data = []
    data_peptide = []
    for k,v in result_dict.items():
        technology = technology_dict[k]
        if technology == 'orbitrap':
            file_name = 'msmsScans_new.txt'
        elif technology == 'bruker':
            file_name = 'accumulatedMsmsScans_new.txt'
        msms_path = os.path.join(v,'combined','txt',file_name)
        msms = pd.read_csv(msms_path,sep='\t')
        msms = msms.loc[msms['Reverse'] != '+',:]   # no matter how, get rid of reverse map first
        if mode == 'fdr':
            msms = msms.loc[msms['Identified']=='+',:]
        elif mode == 'score':
            msms = msms.loc[msms['Score']>=40,:]
        elif mode == 'all':
            msms = msms.loc[msms['Sequence'].notna(),:]
        final_identity_col = []
        for item,pep in zip(msms['Proteins'],msms['Sequence']):
            lis = item.split(';')
            identities = list(set([classify_source(element) for element in lis]))
            if 'contaminant' in identities:
                final_identity = 'contaminant'
            elif 'reverse' in identities:
                final_identity = 'reverse'
            elif 'self_gene' in identities:
                final_identity = 'self_gene'
            else:
                if len(identities) == 1:
                    final_identity = identities[0]
                else:
                    if 'nc_isoform' in identities:
                        identities.remove('nc_isoform')
                        if len(identities) == 1:
                            final_identity = identities[0]
                        else:
                            final_identity = 'non_canonical_ambiguous'   
                    else:

                        if 'nuORF' in identities:
                            if (pep not in bl) and (pep in wl):
                                identities.remove('nuORF')
                                if len(identities) == 1:
                                    final_identity = identities[0]
                                else:
                                    final_identity = 'non_canonical_ambiguous'
                            else:
                                final_identity = 'non_canonical_ambiguous'   
                        else:
                            identities.sort()
                            if identities == ['ERV','TE_chimeric_transcript']:
                                final_identity = 'ERV'
                            else:
                                final_identity = 'non_canonical_ambiguous'
                    

            final_identity_col.append(final_identity)
        msms['final_identity'] = final_identity_col
        n_self_gene = msms.loc[msms['final_identity']=='self_gene',:].shape[0]
        n_splicing = msms.loc[msms['final_identity']=='splicing',:].shape[0]
        n_rna_edit = msms.loc[msms['final_identity']=='rna_edit',:].shape[0]
        n_erv = msms.loc[msms['final_identity']=='ERV',:].shape[0]
        n_circ = msms.loc[msms['final_identity']=='circRNA',:].shape[0]
        n_te_chimeric = msms.loc[msms['final_identity']=='TE_chimeric_transcript',:].shape[0]
        n_intron_retention = msms.loc[msms['final_identity']=='intron_retention',:].shape[0]
        n_pathogen = msms.loc[msms['final_identity']=='pathogen',:].shape[0]
        n_variant = msms.loc[msms['final_identity']=='variant',:].shape[0]
        n_fusion = msms.loc[msms['final_identity']=='fusion',:].shape[0]
        n_contaminant = msms.loc[msms['final_identity']=='contaminant',:].shape[0]
        n_nuorf = msms.loc[msms['final_identity']=='nuORF',:].shape[0]
        n_ambiguous = msms.loc[msms['final_identity']=='non_canonical_ambiguous',:].shape[0]
        n_nc_isoform = msms.loc[msms['final_identity']=='nc_isoform',:].shape[0]
        data.append((k,n_self_gene,n_splicing,n_rna_edit,n_erv,n_circ,n_te_chimeric,n_intron_retention,n_pathogen,n_variant,n_fusion,n_contaminant,n_nuorf,n_ambiguous,n_nc_isoform))

        n_self_gene = len(msms.loc[msms['final_identity']=='self_gene',:]['Sequence'].unique())
        n_splicing = len(msms.loc[msms['final_identity']=='splicing',:]['Sequence'].unique())
        n_rna_edit = len(msms.loc[msms['final_identity']=='rna_edit',:]['Sequence'].unique())
        n_erv = len(msms.loc[msms['final_identity']=='ERV',:]['Sequence'].unique())
        n_circ = len(msms.loc[msms['final_identity']=='circRNA',:]['Sequence'].unique())
        n_te_chimeric = len(msms.loc[msms['final_identity']=='TE_chimeric_transcript',:]['Sequence'].unique())
        n_intron_retention = len(msms.loc[msms['final_identity']=='intron_retention',:]['Sequence'].unique())
        n_pathogen = len(msms.loc[msms['final_identity']=='pathogen',:]['Sequence'].unique())
        n_variant = len(msms.loc[msms['final_identity']=='variant',:]['Sequence'].unique())
        n_fusion = len(msms.loc[msms['final_identity']=='fusion',:]['Sequence'].unique())
        n_contaminant = len(msms.loc[msms['final_identity']=='contaminant',:]['Sequence'].unique())
        n_nuorf = len(msms.loc[msms['final_identity']=='nuORF',:]['Sequence'].unique())
        n_ambiguous = len(msms.loc[msms['final_identity']=='non_canonical_ambiguous',:]['Sequence'].unique())
        n_nc_isoform = len(msms.loc[msms['final_identity']=='nc_isoform',:]['Sequence'].unique())
        data_peptide.append((k,n_self_gene,n_splicing,n_rna_edit,n_erv,n_circ,n_te_chimeric,n_intron_retention,n_pathogen,n_variant,n_fusion,n_contaminant,n_nuorf,n_ambiguous,n_nc_isoform))

    df = pd.DataFrame.from_records(data=data,columns=['sample','n_self_gene','n_splicing','n_rna_edit','n_erv','n_circ','n_te_chimeric','n_intron_retention','n_pathogen','n_variant','n_fusion','n_contaminant','n_nuorf','non_canonical_ambiguous','nc_isoform'])
    df.to_csv(os.path.join(outdir,'source_spectrum.txt'),sep='\t',index=None)

    df = pd.DataFrame.from_records(data=data_peptide,columns=['sample','n_self_gene','n_splicing','n_rna_edit','n_erv','n_circ','n_te_chimeric','n_intron_retention','n_pathogen','n_variant','n_fusion','n_contaminant','n_nuorf','non_canonical_ambiguous','nc_isoform'])
    df.to_csv(os.path.join(outdir,'source_peptide.txt'),sep='\t',index=None)

def get_each_category():
    # get each category
    for abe in ['self_gene','splicing','rna_edit','ERV','circRNA','TE_chimeric_transcript','intron_retention','pathogen','variant','fusion','non_canonical_ambiguous','nuORF','nc_isoform','unknown']:
        data_abe = []
        for k,v in result_dict.items():
            technology = technology_dict[k]
            if technology == 'orbitrap':
                file_name = 'msmsScans_new.txt'
            elif technology == 'bruker':
                file_name = 'accumulatedMsmsScans_new.txt'
            msms_path = os.path.join(v,'combined','txt',file_name)
            msms = pd.read_csv(msms_path,sep='\t')
            msms = msms.loc[msms['Reverse'] != '+',:]   # no matter how, get rid of reverse map first
            if mode == 'fdr':
                msms = msms.loc[msms['Identified']=='+',:]
            elif mode == 'score':
                msms = msms.loc[msms['Score']>=40,:]
            elif mode == 'all':
                msms = msms.loc[msms['Sequence'].notna(),:]
            final_identity_col = []
            for item,pep in zip(msms['Proteins'],msms['Sequence']):
                lis = item.split(';')
                identities = list(set([classify_source(element) for element in lis]))

                if 'contaminant' in identities:
                    final_identity = 'contaminant'
                elif 'reverse' in identities:
                    final_identity = 'reverse'
                elif 'self_gene' in identities:
                    final_identity = 'self_gene'
                else:
                    if len(identities) == 1:
                        final_identity = identities[0]
                    else:
                        if 'nc_isoform' in identities:
                            identities.remove('nc_isoform')
                            if len(identities) == 1:
                                final_identity = identities[0]
                            else:
                                final_identity = 'non_canonical_ambiguous'   
                        else:
                            if 'nuORF' in identities:
                                if (pep not in bl) and (pep in wl):
                                    identities.remove('nuORF')
                                    if len(identities) == 1:
                                        final_identity = identities[0]
                                    else:
                                        final_identity = 'non_canonical_ambiguous'
                                else:
                                    final_identity = 'non_canonical_ambiguous'   
                            else:
                                identities.sort()
                                if identities == ['ERV','TE_chimeric_transcript']:
                                    final_identity = 'ERV'
                                else:
                                    final_identity = 'non_canonical_ambiguous'
                final_identity_col.append(final_identity)
            msms['final_identity'] = final_identity_col
            msms = msms.loc[msms['final_identity']==abe,:]
            data_abe.append(msms)
        result_abe = pd.concat(data_abe,axis=0,keys=list(result_dict.keys())).reset_index(level=-2)
        result_abe.to_csv(os.path.join(outdir,'{}_neoantigen.txt'.format(abe)),sep='\t',index=None)


    # for nuorf, filter by epitope in normal
    nuorf = pd.read_csv(os.path.join(outdir,'nuORF_neoantigen.txt'),sep='\t')
    nuorf['in_normal'] = [True if item in bl else False for item in nuorf['Sequence']]
    nuorf['in_cancer'] = [True if item in wl else False for item in nuorf['Sequence']]
    nuorf.to_csv(os.path.join(outdir,'nuORF_neoantigen.txt'),sep='\t',index=None)

def run_netMHCpan(software_path,peptides,hlas,length,cmd_num=1,tmp_dir=None,tmp_name=None):
    '''
    on top of SNAF function, I added -BA and extract different column from the result
    also make sure to create scratch folder
    right now, only consider cmd1
    '''
    # set the default
    if tmp_dir is None:
        tmp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')
    if tmp_name is None:
        tmp_name = 'input_{}.pep'.format(os.getpid())
    # reformat/create to the strings that we need
    peptides_path = os.path.join(tmp_dir,tmp_name)
    with open(peptides_path,'w') as f:
        for pep in peptides:
            f.write('{}\n'.format(pep))
    # netMHCpan has a hard limit for 1024 characters for -a, so to be safe, 90 HLA as max for one goal, each input format 11 character 'HLA-B57:01,'
    if len(hlas) <= 90:
        hla_strings = ','.join(hlas)
        if cmd_num == 1:
            reconstruct = ' || '.join(['$3 == ' + '"{}"'.format(pep) for pep in peptides])
            cmd = '{} -p {} -a {} -BA -l {} | awk \'BEGIN {{OFS = "\\t"}} {{if ({}) {{print $3, length($3),$2,$13,$16,$18}}}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
            '''
            ../external/netMHCpan-4.1/netMHCpan -p ./test.pep -a HLA-A01:01,HLA-A02:01 -BA -l 9 | awk 'BEGIN {OFS = "\t"} {if ($3 == "AAAWYLWEV" || $3 == "AAGLQDCTM" || $3 == "AARNIVRRA") {print $3, length($3),$2,$13,$16,$18}}'
            '''
        try:
            df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
            df.columns = ['peptide','mer','hla','score','binding_nM','identity']
        except:   # no stdout, just no candidates
            df = pd.DataFrame(columns=['peptide','mer','hla','score','binding_nM','identity'])


    else:
        total = len(hlas)   # 91,137,180
        df_store = []
        i = 0 
        while i < total:
            if i + 90 <= total:
                batch_hlas = hlas[i:i+90]
            else:
                batch_hlas = hlas[i:]
            hla_strings = ','.join(batch_hlas)
            if cmd_num == 1:
                reconstruct = ' || '.join(['$3 == ' + '"{}"'.format(pep) for pep in peptides])
                cmd = '{} -p {} -a {} -BA -l {} | awk \'BEGIN {{OFS = "\\t"}} {{if ({}) {{print $3, length($3),$2,$13,$16,$18}}}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
                '''
                ../external/netMHCpan-4.1/netMHCpan -p ./test.pep -a HLA-A01:01,HLA-A02:01 -BA -l 9 | awk 'BEGIN {OFS = "\t"} {if ($3 == "AAAWYLWEV" || $3 == "AAGLQDCTM" || $3 == "AARNIVRRA") {print $3, length($3),$2,$13,$16,$18}}'
                '''
            try:
                df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
                df.columns = ['peptide','mer','hla','score','binding_nM','identity']
            except:   # no stdout, just no candidates
                df = pd.DataFrame(columns=['peptide','mer','hla','score','binding_nM','identity'])
            df_store.append(df)
            i += 90  
        df = pd.concat(df_store,axis=0)      

    # remove the scratch_pid folder
    os.remove(peptides_path)

    return df

def run_netMHCpan_wrapper(peptide,hlas):
    netMHCpan_path = fixed_config['netMHCpan_path']
    df = run_netMHCpan(software_path=netMHCpan_path,
                       peptides=[peptide],
                       hlas=hlas,
                       length=len(peptide),
                       cmd_num=1,
                       tmp_dir=None,
                       tmp_name=None)
    return df

def run_netMHCpanII(peptide,hlas):
    '''
    df = run_netMHCpanII(peptides=['ASQKRPSQRHGSKASQKRPSQRHGS'],hlas=['DRB1_0101','DRB1_0102','HLA-DPA10103-DPB10201','HLA-DQA10201-DQB10201'])
    '''
    
    software_path = fixed_config['netMHCIIpan_path']
    length = len(peptide[0])
    tmp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')
    tmp_name = 'input_{}.fas'.format(os.getpid())
    peptides_path = os.path.join(tmp_dir,tmp_name)
    with open(peptides_path,'w') as f:
        for i,pep in enumerate(peptide):
            f.write('>seq{}\n{}\n'.format(str(i+1),pep))
    '''
    ../netMHCIIpan -inptype 0 -f example.fsa -a DRB1_0101,DRB1_0102 -length 15 -BA | awk 'BEGIN {OFS = "\t"} { if ($3 == "ASQKRPSQRHGSKYR" || $3 == "ASQKRPSQRHGSKYR") {print $3,length($3),$2,$10,$11,$14}}'
    '''

    # netMHCpan has a hard limit for 1024 characters for -a, so to be safe, 40 HLA as max for one goal HLA-DQA10201-DQB10201, DRB1_0101
    if len(hlas) <= 40:
        hla_strings = ','.join(hlas)
        reconstruct = ' || '.join(['$3 == ' + '"{}"'.format(pep) for pep in peptide])
        cmd = '{} -inptype 0 -f {} -a {} -BA -length {} | awk \'BEGIN {{OFS = "\\t"}} {{if ({}) {{print $3,length($3),$2,$10,$14,$16}}}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
        try:
            df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
            df.columns = ['peptide','mer','hla','score','binding_nM','identity']
        except:   # no stdout, just no candidates
            df = pd.DataFrame(columns=['peptide','mer','hla','score','binding_nM','identity'])


    else:
        total = len(hlas)   # 91,137,180
        df_store = []
        i = 0 
        while i < total:
            if i + 40 <= total:
                batch_hlas = hlas[i:i+40]
            else:
                batch_hlas = hlas[i:]
            hla_strings = ','.join(batch_hlas)
            reconstruct = ' || '.join(['$3 == ' + '"{}"'.format(pep) for pep in peptide])
            cmd = '{} -inptype 0 -f {} -a {} -BA -length {} | awk \'BEGIN {{OFS = "\\t"}} {{if ({}) {{print $3,length($3),$2,$10,$14,$16}}}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
            try:
                df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
                df.columns = ['peptide','mer','hla','score','binding_nM','identity']
            except:   # no stdout, just no candidates
                df = pd.DataFrame(columns=['peptide','mer','hla','score','binding_nM','identity'])
            df_store.append(df)
            i += 40 
        df = pd.concat(df_store,axis=0)      

    # remove the scratch_pid folder
    os.remove(peptides_path)
    return df




def hla_formatting(pre,pre_type,post_type):
    if pre_type == 'netMHCpan_output' and post_type == 'netMHCpan_input':  # HLA-A*01:01 to HLA-A01:01
        post = [hla.replace('*','') for hla in pre]
    elif pre_type == 'netMHCpan_input' and post_type == 'netMHCpan_output':  # HLA-A01:01 to HLA-A*01:01
        post = [hla[:5] + '*' + hla[5:] for hla in pre]
    elif pre_type == 'netMHCpan_output' and post_type == 'deepimmuno':  # HLA-A*01:01 to HLA-A*0101
        post = [hla.replace(':','') for hla in pre]
    elif pre_type == 'deepimmuno' and post_type == 'netMHCpan_output': # HLA-A*0101 to HLA-A*01:01
        post = []
        for hla in pre:
            first = hla[:6] # HLA-A*
            second = hla[6:8] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[8:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    elif pre_type == 'deepimmuno_nostar' and post_type == 'netMHCpan_input':   # HLA-A0101 to HLA-A01:01
        post = []
        for hla in pre:
            first = hla[:5] # HLA-A
            second = hla[5:7] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[7:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    return post

def modify_df(df):
    # change level_0
    col = []
    for study,sample in zip(df['level_0'],df['Raw file']):
        uid = study + ',' + sample
        bio = uid2bio[uid]
        col.append(bio)
    df['level_0'] = col
    return df


def inquiry(tup,filter=True,add_hlas=None,mode='i'):
    # get n_psm and highest_abundance
    pep,typ = tup
    df = modified_df_neo_dict[typ]
    df = df.loc[df['Sequence']==pep,:]
    n_psm = df.shape[0]
    is_only_rescore = False
    if 'Identified_rescore' in df.columns and (not np.any(df['Identified_vanilla'].values)):
        is_only_rescore = True

    try:
        highest_abundance = round(max([float(item) for item in df['Precursor intensity']]))
    except:  # maxquant weird edge case
        highest_abundance = None
    try:
        highest_score = round(max([float(item) for item in df['Score']]))
    except:  # maxquant weird edge case
        highest_score = None
    try:
        best_pep = min([float(item) for item in df['PEP']])
    except:  # maxquant weird edge case
        best_pep = None
    samples = ','.join([str(item) for item in df['level_0'].unique().tolist()])
    # whether can be presented by sample-specific alleles
    info_dic = {}
    for s in samples.split(','):
        info = []
        raw_hla = bio2hla[s]
        if isinstance(raw_hla,str):
            if mode == 'i':
                hlas = ['HLA-' + item.replace('*','') for item in raw_hla.split(',')]  # HLA-A02:01
                result = run_netMHCpan_wrapper(pep,hlas)
            elif mode == 'ii':
                hlas = [item for item in raw_hla.split(',')]
                result = run_netMHCpanII([pep],hlas)
            if filter:
                result = result.loc[result['identity'].isin(['SB','WB']),['hla','score','binding_nM','identity']]
            else:
                result = result.loc[:,['hla','score','binding_nM','identity']]
            for row in result.itertuples():
                info.append((row.hla,row.score,row.binding_nM,row.identity))
        else:
            info.append((None,None,None,None))
        info_dic[s] = info
    # additional hla alleles to query
    if add_hlas is not None:
        info = []
        if mode == 'i':
            result = run_netMHCpan_wrapper(pep,add_hlas)
        elif mode == 'ii':
            result = run_netMHCpanII([pep],add_hlas)
        if filter:
            result = result.loc[result['identity'].isin(['SB','WB']),['hla','score','binding_nM','identity']]
        else:
            result = result.loc[:,['hla','score','binding_nM','identity']]
        for row in result.itertuples():
            info.append((row.hla,row.score,row.binding_nM,row.identity))
    return pep,typ,n_psm,is_only_rescore,highest_abundance,highest_score,best_pep,samples,info_dic,info




def split_array_to_chunks(array,cores=None):
    if not isinstance(array,list):
        raise Exception('split_array_to_chunks function works for list, not ndarray')
    array_index = np.arange(len(array))
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(array_index,cores)
    sub_arrays = []
    for sub_index in sub_indices:
        item_in_group = []
        for i in sub_index:
            item_in_group.append(array[i])
        sub_arrays.append(item_in_group)
    return sub_arrays

def modify_df(df):
    # change level_0
    col = []
    for study,sample in zip(df['level_0'],df['Raw file']):
        uid = study + ',' + sample
        bio = uid2bio[uid]
        col.append(bio)
    df['level_0'] = col
    return df

def create_full_list():
    full_list = []
    for abe in categories:
        if abe == 'self_gene':
            lis = []
            gene_lfc = pd.read_csv(os.path.join(atlas_dir,'gene_lfc.txt'),sep='\t',index_col=0)
            gene_lfc = gene_lfc.loc[(gene_lfc['median_tumor']>20) & (gene_lfc['median_tumor'] > gene_lfc['max_median_gtex']),:]
            real_common = list(set(gene_lfc.index).intersection(set(common)))
            gene_lfc = gene_lfc.loc[real_common,:]
            ts_ensgs = set(gene_lfc.index.tolist())
            pat = re.compile(r'(ENSG\d+)\|ENST')
            df = pd.read_csv(os.path.join(outdir,'self_gene_neoantigen.txt'),sep='\t')
            for item1,item2 in zip(df['Sequence'],df['Proteins']):
                ensgs = []
                for actual_hit in item2.split(';'):
                    match = re.search(pat,actual_hit)
                    if match:
                        ensg = match.group(1)
                        if ensg in ts_ensgs:
                            lis.append((item1,abe))
                            break
            full_list.extend(list(set(lis)))
            modified_df_neo_dict[abe] = modify_df(df)
                        
        elif abe == 'nuORF':
            lis = []
            df = pd.read_csv(os.path.join(outdir,'nuORF_neoantigen.txt'),sep='\t')
            df = df.loc[(~df['in_normal']) & df['in_cancer'],:]
            unique = df['Sequence'].unique().tolist()
            for item in unique:
                lis.append((item,abe))
            full_list.extend(lis)
            modified_df_neo_dict[abe] = modify_df(df)
        else:
            lis = []
            df = pd.read_csv(os.path.join(outdir,'{}_neoantigen.txt'.format(abe)),sep='\t')
            unique = df['Sequence'].unique().tolist()
            for item in unique:
                lis.append((item,abe))
            full_list.extend(lis)
            modified_df_neo_dict[abe] = modify_df(df)


    full_list = full_list + canonical_spike_in
    return full_list

def each_chunk_func(sub_full_list):
    data = []
    for tup in tqdm(sub_full_list):
        data.append(inquiry(tup,filter,final_hla,inquiry_mode))   
    final = pd.DataFrame.from_records(data,columns=['pep','typ','n_psm','is_only_rescore','highest_abundance','highest_score','best_pep','samples','presented_by_each_sample_hla','additional_query'])
    return final

def get_main_binding_result():
    pool = mp.Pool(processes=cores)
    list_of_sub_full_list = split_array_to_chunks(full_list,cores)
    r = [pool.apply_async(func=each_chunk_func,args=(sub_full_list,)) for sub_full_list in list_of_sub_full_list]
    pool.close()
    pool.join()
    results = []
    for collect in r:
        result = collect.get()
        results.append(result)
    final = pd.concat(results,axis=0)
    return final

def added_relative_abundance():
    # add relative abundance
    ori_dir = os.getcwd()
    os.chdir(outdir)
    antigen_files = subprocess.run("for f in *_neoantigen.txt; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    os.chdir(outdir)
    total = []
    for antigen_file in antigen_files:
        df = pd.read_csv(os.path.join(outdir,antigen_file),sep='\t')
        total.append(df)
    total = pd.concat(total,axis=0)
    max_rows = []
    for peptide,sub_df in total.groupby(by='Sequence'):
        sub_df['Precursor intensity'] = sub_df['Precursor intensity'].fillna(value=0)
        sub_df.sort_values(by='Precursor intensity',ascending=False,inplace=True)
        max_rows.append(sub_df.iloc[[0],:])
    max_total = pd.concat(max_rows,axis=0).sort_values(by='Precursor intensity',ascending=False)
    max_total['percentile'] = [1-(i+1)/max_total.shape[0] for i in range(max_total.shape[0])]
    seq2pert = pd.Series(index=max_total['Sequence'].values.tolist(),data=max_total['percentile'].values.tolist()).to_dict()
    final['relative_abundance'] = final['pep'].map(seq2pert).values
    return final

def added_detailed_abundance():
    ori_dir = os.getcwd()
    os.chdir(outdir)
    antigen_files = subprocess.run("for f in *_neoantigen.txt; do echo $f; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    os.chdir(ori_dir)
    total = []
    for antigen_file in antigen_files:
        df = pd.read_csv(os.path.join(outdir,antigen_file),sep='\t')
        total.append(df)
    total = pd.concat(total,axis=0)
    raw_file_intensity = {}
    for raw_file, sub_df in total.groupby(by='Raw file'):
        sub_df['Precursor intensity'] = sub_df['Precursor intensity'].fillna(value=0)
        raw_file_intensity[raw_file] = np.sort(sub_df['Precursor intensity'].values)  # ascending order
    seq2pep_norm_intensity = {}
    for peptide,sub_df in total.groupby(by='Sequence'):
        if peptide in all_need_peptides:
            sub_df['Precursor intensity'] = sub_df['Precursor intensity'].fillna(value=0)
            peptide_norm_intensity = []
            for raw_file, sub_df2 in sub_df.groupby(by='Raw file'):
                all_intensity = raw_file_intensity[raw_file]
                max_intensity = sub_df2['Precursor intensity'].values.max()
                # normalized_intensity = np.interp(x=max_intensity,xp=[all_intensity.min(),all_intensity.max()],fp=[0,1])
                normalized_intensity = bisect.bisect(all_intensity,max_intensity) / len(all_intensity)
                peptide_norm_intensity.append(normalized_intensity)
            seq2pep_norm_intensity[peptide] = peptide_norm_intensity
    final['detailed_intensity'] = final['pep'].map(seq2pep_norm_intensity).values
    return final

def added_source():
    col = []
    # first load in all neoantigen file
    all_df_dict = {}
    for abe in inquiry_categories:
        path = '{}_neoantigen.txt'.format(abe)
        antigen_df = pd.read_csv(os.path.join(outdir,path),sep='\t')
        all_df_dict[abe] = antigen_df
    # now start to do that
    for row in tqdm(final.itertuples()):
        typ = row.typ
        pep = row.pep
        antigen_df = all_df_dict[typ]
        pep2protein = {k:v for k,v in zip(antigen_df['Sequence'],antigen_df['Proteins'])}
        col.append(pep2protein[pep])
    final['source'] = col
    return final

def added_canonical_info():
    # add canonical gene information
    gene_lfc_path = os.path.join(atlas_dir,'gene_lfc.txt')
    self_gene_path = os.path.join(outdir,'self_gene_neoantigen.txt')
    depmap_path = fixed_config['depmap_path']
    gprofile_convert_path = fixed_config['depmap_gprofile_convert_path']
    hla_ligand_atlas_path = fixed_config['hla_ligand_atlas_path']

    gene_lfc = pd.read_csv(gene_lfc_path,sep='\t',index_col=0)
    self_gene = pd.read_csv(self_gene_path,sep='\t',index_col=0)
    depmap = pd.read_csv(depmap_path,sep=',',index_col=0)
    hla_ligand_atlas = pd.read_csv(hla_ligand_atlas_path,sep='\t',index_col=0)

    ensg2mt = gene_lfc['median_tumor'].to_dict()
    ensg2mmg = gene_lfc['max_median_gtex'].to_dict()
    ensg2symbol = gene_lfc['gene_symbol'].to_dict()

    # depmap
    convert = pd.read_csv(gprofile_convert_path,sep=',',index_col=0)['converted_alias'].to_dict()
    ensg2depmap = {}
    for symbol in depmap.columns:
        ensg = convert.get(symbol,None)
        if ensg is not None:
            ensg2depmap[ensg] = np.nanmedian(depmap[symbol].values)

    # hla ligand
    if inquiry_mode == 'i':
        hla_ligand_atlas = hla_ligand_atlas.loc[hla_ligand_atlas['hla_class']!='HLA-II',:]
    elif inquiry_mode == 'ii':
        hla_ligand_atlas = hla_ligand_atlas.loc[hla_ligand_atlas['hla_class']!='HLA-I',:]
    peptide2tissues = pd.Series(index=hla_ligand_atlas['peptide_sequence'].values,data=hla_ligand_atlas['tissues'].values).to_dict()


    peptide2ensg = {}
    pat = re.compile(r'(ENSG\d+)\|ENST')
    for item1,item2 in zip(self_gene['Sequence'],self_gene['Proteins']):
        for actual_hit in item2.split(';'):
            match = re.search(pat,actual_hit)
            if match:
                ensg = match.group(1)
                peptide2ensg.setdefault(item1,[]).append(ensg)

    col1 = []   # ensgs
    col2 = []   # mm
    col3 = []   # mmg
    col4 = []   # symbol
    col5 = []   # unique
    col6 = []   # depmap
    col7 = []   # hla ligand atlas
 
    for item,typ in zip(final['pep'],final['typ']):
        if typ == 'self_gene':
            ensgs = peptide2ensg[item]
            ensgs = list(set(ensgs))
            if len(ensgs) == 1:
                unique_ensg = ensgs[0]
                col1.append(unique_ensg)
                col2.append(ensg2mt[unique_ensg])
                col3.append(ensg2mmg[unique_ensg])
                col4.append(ensg2symbol[unique_ensg])
                col5.append(True)
                col6.append(ensg2depmap.get(unique_ensg,None))
                col7.append(peptide2tissues.get(item,None))

            else:
                col1.append(ensgs)
                col2.append([ensg2mt.get(ensg,None) for ensg in ensgs])
                col3.append([ensg2mmg.get(ensg,None) for ensg in ensgs])
                col4.append([ensg2symbol.get(ensg,None) for ensg in ensgs])
                col5.append(False)
                col6.append([ensg2depmap.get(ensg,None) for ensg in ensgs])
                col7.append([peptide2tissues.get(ensg,None) for ensg in ensgs])
        else:
            col1.append(None)
            col2.append(None)
            col3.append(None)
            col4.append(None)
            col5.append(None)
            col6.append(None)
            col7.append(None)

    final['ensgs'] = col1
    final['median_tumor'] = col2
    final['max_median_gtex'] = col3
    final['gene_symbol'] = col4
    final['unique'] = col5
    final['depmap_median'] = col6
    final['hla_ligand_atlas'] = col7
    return final

def get_common_and_hla():
    # get the canonical gene cand
    bayests_cutoff = 0.3
    bayests_xy_path = fixed_config['bayests_xy_path']
    bayests = pd.read_csv(bayests_xy_path,sep='\t',index_col=0)
    ts_gene = bayests.loc[bayests['BayesTS']<bayests_cutoff,:].index.tolist()
    common = ts_gene


    # if it's target data, use deg 
    try:
        deg = pd.read_csv(os.path.join(atlas_dir,'deg.txt'),sep='\t',index_col=1)
    except:
        pass
    else:
        deg = deg.loc[(deg['adjp']<0.05) & (deg['Log2(Fold Change)']>0.58),:].index.tolist()
        common = list(set(ts_gene).intersection(set(deg)))

    added = [] if fixed_config['added_genes'] is None else fixed_config['added_genes']
    common = common + added

    if inquiry_mode == 'ii':
        gene_lfc = pd.read_csv(os.path.join(atlas_dir,'gene_lfc.txt'),sep='\t',index_col=0)
        common = gene_lfc.index.tolist()

    # update this
    canonical_spike_in = []

    cutoff = 0.05
    us_hla = pd.read_csv(fixed_config['us_hla_path'],sep=',',index_col=0)
    frequent_hla = hla_formatting(us_hla.loc[us_hla['Percent US population']>cutoff,:].index.tolist(),'deepimmuno_nostar','netMHCpan_input') 

    # update this
    additional_hla = []
    final_hla = list(set(frequent_hla).union(set(additional_hla)))

    if inquiry_mode == 'ii':
        final_hla = ['HLA-DPA10103-DPB10201','HLA-DQA10102-DQB10501']

    # ### overwrite here for hla-ii
    # final_hla = ['HLA-DPA10103-DPB10201','HLA-DQA10102-DQB10501']

    return common,final_hla,canonical_spike_in


def added_non_canonical_unique():
    col = []
    for item1,item2,item3 in zip(final['typ'],final['source'],final['unique']):
        if item1 != 'self_gene':
            if len(item2.split(';')) == 1:
                col.append(True)
            elif len(item2.split(';')) >= 2:
                if item1 != 'nuORF' and 'nuORF' in item2:
                    col.append(True)
                elif item1 == 'ERV' and 'TE_info' in item2:
                    col.append(True)
                else:
                    col.append(False)
            else:
                col.append(False)

        else:
            col.append(item3)

    final['unique'] = col
    return final

def added_nuorf_type():
    col = []
    df = pd.read_csv(nuorf_db,sep='\t',index_col=0)
    uid2typ = df['plotType'].to_dict()
    for item1,item2 in zip(final['typ'],final['source']):
        if item1 == 'nuORF':
            uid = item2.split('|')[0]
            typ = uid2typ.get(uid,'unknown')
            col.append(typ)
        else:
            col.append(None)
    final['nuorf_type'] = col
    return final



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Configure MaxQuant Python')
    parser.add_argument('--config',type=str,default='',help='path to the config json')
    args = parser.parse_args()
    config_path = args.config
    with open(config_path,'r') as f:
        config = json.load(f)['maxquant']

    tunable_config = config['tunable_parameters']
    fixed_config = config['fixed_parameters']

    root_dir = tunable_config['immunopeptidome_dir']
    atlas_dir = tunable_config['outdir']
    result_dict = {item:os.path.join(root_dir,item) for item in tunable_config['result_dict'].values()}  # CC1: path/to/immuno/CC1, when immuno file name is not interpretable
    mode = tunable_config['mode']
    technology_dict = tunable_config['technology']
    cores = tunable_config['cores']
    inquiry_mode = tunable_config['inquiry_mode']
    filter = True
    inquiry_categories = fixed_config['inquiry_categories']
    nuorf_db = fixed_config['nuorf']

    metadata_path = tunable_config['metadata_path']
    metadata = pd.read_csv(metadata_path,sep='\t')
    uid2bio = {}
    bio2hla = {}
    for row in metadata.itertuples():
        if row.sample.endswith('.raw'):
            ext = '.raw'
        elif row.sample.endswith('.d'):
            ext = '.d'
        if isinstance(row.batch,str):
            uid = row.study + '_{}'.format(row.batch) + ',' + row.sample.split(ext)[0]
        else:
            uid = row.study + ',' + row.sample.split(ext)[0]
        bio = row.biology
        hla = row.HLA
        uid2bio[uid] = bio
        bio2hla[bio] = hla

    modified_df_neo_dict = {}  # basically, change the level_0 to biology, also the batch folder needs to be prefixed with study name


    # outdir
    if mode == 'fdr':
        folder = 'fdr'
    elif mode == 'score':
        folder = 'score'
    elif mode == 'all':
        folder = 'all'
    outdir = os.path.join(atlas_dir,fixed_config['antigen_dir'],folder)  
    
    categories = inquiry_categories
    normal = pd.read_csv(fixed_config['ieatlas_normal'],sep='\t')
    cancer = pd.read_csv(fixed_config['ieatlas_cancer'],sep='\t')
    bl = set(normal['Sequence'].tolist())
    wl = set(cancer['Sequence'].tolist())

    if not os.path.exists(outdir):
        os.makedirs(outdir)



    # summarization
    get_stats_spectra()
    get_stats_source()
    get_each_category()

    # generate results
    if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')):
        os.mkdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch'))
    common, final_hla, canonical_spike_in = get_common_and_hla()

    full_list = create_full_list()
    final = get_main_binding_result()

    all_need_peptides = set(final['pep'].values.tolist())
    final = added_relative_abundance()
    final = added_detailed_abundance()
    final = added_source()
    final = added_canonical_info()
    final = added_non_canonical_unique()
    final = added_nuorf_type()
    final.to_csv(os.path.join(outdir,'final_enhanced.txt'),sep='\t',index=None)









