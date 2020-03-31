#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 09:38:53 2020

@author: ligk2e
"""

import os
os.chdir('/Users/ligk2e/Desktop')
import pandas as pd
import pickle
import pick_peptide as pp
import ast
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def extract_EnsID(df):
    
    UID = list(df['UID'])
    EnsID_array = []
    for item in UID:
        EnsID = item.split('|')[0].split(':')[1]
        EnsID_array.append(EnsID)
    return EnsID_array

def write_list_to_file(list):
    with open('EnsID4query.txt','w') as f1:
        f1.writelines('%s\n' % EnsID_i for EnsID_i in list)
    return None

def find_membrane_EnsID(conversion_table,uniprot_info):
    dict_EnsID_uni = {}
    for i in range(conversion_table.shape[0]):
        dict_EnsID_uni[conversion_table['EnsID'][i]] = conversion_table['uniprot_entry'][i]
    extracellular = []
    for key in dict_EnsID_uni.keys():
        if dict_EnsID_uni[key] in list(uniprot_info['Entry']):  # remember generator with list
            extracellular.append(key)
    return extracellular, dict_EnsID_uni
            

def representative_tran_and_whole_tran(df):
    # get the narrowed(3745-318) whole transcript from saved pickle file
    with open('col1_i.p','rb') as f1:
        exam_col = pickle.load(f1)
    truth_table = [True if item in extracellular_gene else False for item in EnsID]
    narrow_whole_tran = []
    for i in range(len(truth_table)):
        if truth_table[i]:
            narrow_whole_tran.append(exam_col[i][2])
    # one on one to find correponding representative and their whole transcript        
    representative = []
    whole_tran = []
    for i in range(df.shape[0]):
        event = list(df['exam_match_tran'])[i]
        event = ast.literal_eval(event)
        # initiate the condition to find the most likely representative
        flag = -1
        max_seq = ''
        max_length = 0
        max_item_score = 0
        position = 0    # how to use flag to record the position info and when to use index/item to loop
        for tran in event:   
            if len(tran) == 0:
                flag += 1
            if not len(tran)==0: 
                flag += 1
                temp1 = len(tran)
                add_score = pp.score_GC(tran) + pp.score_coding_bias(tran)
                if temp1 > max_length:
                    max_length = temp1
                    max_item_score = add_score
                    max_seq = tran
                    position = flag
                elif temp1 == max_length:
                    if add_score > max_item_score:
                        max_length = temp1
                        max_item_score = add_score
                        max_seq = tran
                        position =flag
        try:
            representative.append(max_seq) 
            whole_tran.append(narrow_whole_tran[i][position])
        except:
            representative.append('')
            whole_tran.append('')
    
    return representative,whole_tran
        #df['representative_tran'] = representative  


                        
def check_if_good_representative(df):
    condition_array = []
    for i in range(df.shape[0]):
        repre = list(df['representative_tran'])[i]
        if repre:
            junction = list(df['exam_seq'])[i]
            whole = list(df['whole_tran'])[i]
            start_repre = whole.find(repre)
            end_repre = whole.find(repre) + len(repre)
            start_junction = whole.find(junction)
            end_junction = whole.find(junction) + len(junction)
            if start_junction <= end_repre and end_junction >= start_repre:
                condition_array.append(True)
            else:
                condition_array.append(False)
        else:
            condition_array.append(False)
    df['good_repre'] = condition_array
    df_filtered = df[df['good_repre']==True]    # how to access all the column name and how to drop one column
    return df_filtered
            
def alignment_to_uniprot(df,dict_fa,dict_EnsID_uni):
    notebook = []
    result_array = []
    repre_aa_array = []
    for i in range(df.shape[0]):
        EnsID = list(df['UID'])[i].split('|')[0].split(':')[1]
        target_aa = dict_fa[dict_EnsID_uni[EnsID]]
        repre_aa = str(Seq(list(df['representative_tran'])[i],generic_dna).translate(to_stop=False))
    # shotgun and align
        bucket = chop_sequence(repre_aa,10)
        notes = []
        for j in range(len(bucket)):
            frag = bucket[j]
            if frag in target_aa:
                notes.append(True)
            else:
                notes.append(False)
        result = neoantigen_iden(notes)
        notebook.append(notes)
        result_array.append(result)
        repre_aa_array.append(repre_aa)
    return notebook,result_array,repre_aa_array
        
def neoantigen_iden(notes):
    if sum(notes) > 0 and sum(notes) < len(notes):
        result = 'partially aligned'
    elif sum(notes) == 0:
        result = 'no align'
    elif sum(notes) == len(notes):
        result = 'totally align'
    return result
    
          

    
def read_uniprot_seq(path):
    dict_fa = {}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            uniID = title.split('|')[1]
            dict_fa[uniID] = seq
    return dict_fa       

def chop_sequence(seq,kmer):   # how to splice sequence, elegant way to use range
    frag_bucket = []
    for i in range(0,len(seq),kmer):
        try:
            frag_bucket.append(seq[i:i+kmer])
        except:
            frag_bucket.append(seq[i:])
    return frag_bucket

if __name__ == "__main__":
    
    ### narrow down the list

    
    # get the EnsID and batch query, get the one-to-one conversion table between EnsID and uniprot iD
    df_ori_version2 = pd.read_csv('/Users/ligk2e/Desktop/add_exon_itself.txt',sep='\t')
    EnsID = extract_EnsID(df_ori_version2)
    write_list_to_file(EnsID)
    # load in the query rsult downloading from uniprot
    EnsID_to_uniprot = pd.read_csv('/Users/ligk2e/Desktop/query_result.tab',sep='\t',header=None,
   names=['EnsID','isoforms','uniprot_entry','Entry_name','protein','length','topology','gene_name'],
   skiprows=1)
    # load all the extracellular protein in human from uniprot
    uniprot_info = pd.read_csv('/Users/ligk2e/Desktop/uniprot_info.tab',sep='\t')
    # find extracellular increasingly spliced EnsID
    extracellular_gene,dict_EnsID_uni = find_membrane_EnsID(EnsID_to_uniprot,uniprot_info)
    # finally narrow down the dataframe to 221 rows
  
    df_ori_version2['condition'] = [True if item in extracellular_gene else False for item in EnsID]
    # this could also be achived by apply function or lambda function
    df_ori_narrow = df_ori_version2[df_ori_version2['condition'] == True]
   
    
    
    #### get the representative and check if it is a good representative

    representative,whole_tran = representative_tran_and_whole_tran(df_ori_narrow)
    df_ori_narrow['representative_tran'] = representative
    df_ori_narrow['whole_tran'] = whole_tran# don't to_csv then read_csv
    df_ori_narrow_good_repre = check_if_good_representative(df_ori_narrow)
    
    ### final alignment
    dict_uni_fa = read_uniprot_seq('/Users/ligk2e/Desktop/uniprot_canonical.fasta')
    comment,alignment,repre_aa = alignment_to_uniprot(df_ori_narrow_good_repre,dict_uni_fa,dict_EnsID_uni)   
    df_ori_narrow_good_repre['comment'] = comment
    df_ori_narrow_good_repre['alignment'] = alignment
    
   
            
    # mannaully check the case
    df_ori_narrow_good_repre['repre_aa'] = repre_aa  

     # export for presentation
    df_ori_narrow_good_repre[['UID','alignment','repre_aa']].to_csv('tenta_result.txt',sep='\t',header=True,index=False)      
    df_ori_narrow_good_repre.to_csv('final_result.txt',sep='\t',header=True,index=False)  
    
    b = df_ori_narrow_good_repre
   
    
    # ENSG00000282228   - strand use case
    # ENSG00000110514   +
    # ENSG00000243646

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    