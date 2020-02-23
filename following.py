#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 17:29:27 2020

@author: ligk2e
"""
import os
os.chdir('/Users/ligk2e/Desktop')
import branch2 as br2
import task1_mod as mich
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser




def retrieve_junction_site(df_ori):
    exam_seq,back_seq = [],[]
    for i in range(df_ori.shape[0]):
        temp = mich.UID(df_ori,i)
        EnsID = list(temp.keys())[0].split(':')[1]
        exam_site = list(temp.values())[0][0]
        back_site = list(temp.values())[0][1]
        exam_site_1 = subexon_tran(exam_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
        exam_site_2 = subexon_tran(exam_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
        exam_seq.append(exam_site_1 + exam_site_2)
        back_site_1 = subexon_tran(back_site.split('-')[0],EnsID,dict_exonCoords,dict_fa,'site1')
        back_site_2 = subexon_tran(back_site.split('-')[1],EnsID,dict_exonCoords,dict_fa,'site2')
        back_seq.append(back_site_1 + back_site_2)
        
    df_ori['exam_seq'] = exam_seq
    df_ori['back_seq'] = back_seq
    return df_ori
        
        
def subexon_tran(subexon,EnsID,dict_exonCoords,dict_fa,flag): # flag means if it is site_1 or site_2
    try:
        attrs = dict_exonCoords[EnsID][subexon]
        exon_seq = br2.query_from_dict_fa(dict_fa,attrs[2],attrs[3],EnsID,attrs[1])
    except:
#        print(subexon,'t')
        try:
            suffix = subexon.split('_')[1]
        except:
            exon_seq = '*'   # '*' means possible gene fusion
            print(subexon,'possible fusion gene event',EnsID)
        else:
            subexon = subexon.split('_')[0]
            try:
                attrs = dict_exonCoords[EnsID][subexon]
            except:
                exon_seq = '#'   # '# means 'U0.1_32130925'
                print(subexon,'splicing occurs in UTR',EnsID)
            else:
                if flag == 'site1':           
                    exon_seq = br2.query_from_dict_fa(dict_fa,suffix,attrs[3],EnsID,attrs[1])
                elif flag == 'site2':
                    exon_seq = br2.query_from_dict_fa(dict_fa,attrs[2],suffix,EnsID,attrs[1])
    return exon_seq
    
        
if __name__ == "__main__":
    df_ori_v2 = pd.read_csv('/Users/ligk2e/Desktop/match_i_v1_Feb20_2.txt',sep='\t')
    dict_exonCoords = br2.exonCoords_to_dict('/Users/ligk2e/Desktop/project/Hs_Ensembl_exon.txt','\t')
    dict_fa = br2.fasta_to_dict('/Users/ligk2e/Desktop/project/Hs_gene-seq-2000_flank.fa')
    new_df = retrieve_junction_site(df_ori_v2)
    new_df.to_csv('/Users/ligk2e/Desktop/add_exon_itself.txt',sep='\t',header=True,index=False)
    
    
    
    
    
    
    
    