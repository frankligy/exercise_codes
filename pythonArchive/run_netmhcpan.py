#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import os,sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
from io import StringIO
import subprocess


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
    # HLA-A02:01
    df = run_netMHCpan(software_path=netMHCpan_path,
                       peptides=[peptide],
                       hlas=hlas,
                       length=len(peptide),
                       cmd_num=1,
                       tmp_dir=None,
                       tmp_name=None)
    return df

def each_chunk_func(sub_list):
    list_df = []
    for pep in tqdm(sub_list):
        peptide = pep
        hlas = pep2hla[peptide]
        df = run_netMHCpan_wrapper(peptide,hlas)
        list_df.append(df)
    final = pd.concat(list_df,axis=0)
    return final


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


# step1: create pep2hla from the test.out.fasta
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq



# df = pd.read_csv('final.txt',sep='\t')
# hlas = ['HLA-A02:01','HLA-A24:02']
# pep2hla = {}
# for row in df.itertuples():
#     pep2hla[row.pep] = hlas
#     pep2hla[row.mutated_pep] = hlas


# netMHCpan_path = '/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/netMHCpan-4.1/netMHCpan'
# cores = 1

# full_list = list(pep2hla.keys())
# pool = mp.Pool(processes=cores)
# list_of_sub_full_list = split_array_to_chunks(full_list,cores)
# r = [pool.apply_async(func=each_chunk_func,args=(sub_full_list,)) for sub_full_list in list_of_sub_full_list]
# pool.close()
# pool.join()
# results = []
# for collect in r:
#     result = collect.get()
#     results.append(result)
# final = pd.concat(results,axis=0)
# final.to_csv('final_results.txt',sep='\t')

# then return back
df = pd.read_csv('final.txt',sep='\t')
result = pd.read_csv('final_results.txt',sep='\t')
dic = {}
for hla,sub_df in result.groupby(by='hla'):
    d = pd.Series(index=sub_df['peptide'].values,data=sub_df['binding_nM'].values).to_dict()
    dic[hla] = d

df['pep_A0201'] = [dic['HLA-A*02:01'][item] for item in df['pep']]
df['mpep_A0201'] = [dic['HLA-A*02:01'][item] for item in df['mutated_pep']]
df['pep_A2402'] = [dic['HLA-A*24:02'][item] for item in df['pep']]
df['mpep_A2402'] = [dic['HLA-A*24:02'][item] for item in df['mutated_pep']]

df.to_csv('to_plot.txt',sep='\t',index=None)

