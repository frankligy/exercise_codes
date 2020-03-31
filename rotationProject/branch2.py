#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 14:05:41 2020

@author: ligk2e
"""
import os
os.chdir('/Users/ligk2e/Desktop/')
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
import task1_mod as mich
from decimal import Decimal as D
import pickle
import debug


def match_with_exonlist(df_ori,df_exonlist,dict_exonCoords):
#    sum = 0
    col1 = []
    col2 = []
    
    for i in range(df_ori.shape[0]):
        temp=mich.UID(df_ori,i)
        EnsID=list(temp.keys())[0].split(':')[1]
#        Exons_examined = exon_update(temp,0,EnsID)
#        Exons_back = exon_update(temp,1,EnsID)
        Exons_examined = exon_extract(temp,0,EnsID)
        Exons_back = exon_extract(temp,1,EnsID)
        col1.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_examined))
        col2.append(core_match(df_exonlist,dict_exonCoords,EnsID,Exons_back))
        
    return col1,col2


def exon_update(temp,pos,EnsID):
    Exons_former = list(temp.values())[0][pos].split('-')[0]
    #print(Exons_former)
    if len(Exons_former) > 7:  # non-canonical/novel splicing sites
        print(Exons_former + ' in ' + EnsID + ' is a non-canonocal case\n')
        Exons_former_update = Exons_former
    elif Exons_former.startswith('E'):
        Exons_former_nume = Exons_former.lstrip('E')
        Exons_former_nume_update = str(D(Exons_former_nume) + D('0.1'))
        Exons_former_update = 'E' + Exons_former_nume_update
    elif Exons_former.startswith('I'):
        Exons_former_nume = Exons_former.lstrip('I')
        Exons_former_nume_update = str(D(Exons_former_nume) + D('0.1'))
        Exons_former_update = 'I' + Exons_former_nume_update
    Exons_latter = list(temp.values())[0][pos].split('-')[1]
    Exons = Exons_former_update + '|' + Exons_latter
    print(Exons)
    return Exons

def exon_extract(temp,pos,EnsID):
    Exons = list(temp.values())[0][pos].split('-')[0] + '|' + list(temp.values())[0][pos].split('-')[1]
    return Exons

def core_match(df_exonlist,dict_exonCoords,EnsID,Exons):
   
    try:
        df_certain = df_exonlist[df_exonlist['EnsGID'] == EnsID]
    except:
        final_fullAA = []
        peek_pep = []
        full_transcript_store = []
        
    final_fullAA = []
    peek_pep = []
    full_transcript_store = []
    for item in list(df_certain['Exons']):
        full_transcript=''
#        splicing_transcript=''
        if Exons in item:
            Exonlist = item.split('|')
            for j in range(len(Exonlist)):
                coords = dict_exonCoords[EnsID][Exonlist[j]]
                frag = frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                full_transcript += frag
#                if check_exonlist(Exonlist,j):   
#                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
#                    
#                else:
#                    frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1]) 
#                full_transcript += frag
            full_transcript = full_transcript.replace('\n','')
#            splicing_transcript = full_transcript.replace('\n','')
            full_transcript_store.append(full_transcript)   
            pot_fullAA=mich.translate(full_transcript)
            peek_pep.append(pot_fullAA)
            max_fullAA=find_longest_AA(list(pot_fullAA.values()))
            final_fullAA.append(max_fullAA)
        else:
            peek_pep.append('')
            final_fullAA.append('')
            full_transcript_store.append('')
    result = [final_fullAA,peek_pep,full_transcript_store]
    return result


def check_exonlist_stra1(exonlist,index):
    dict = {}
    for subexon in exonlist:
        exon_num = subexon.split('.')[0]
        if exon_num in dict:
            dict[exon_num].append(subexon)
        else:
            dict[exon_num] = []
            dict[exon_num].append(subexon)
    # check
    query = exonlist[index]
    if len(dict[query.split('.')[0]]) == 1:
        return True
    else:
        if query == dict[query.split('.')[0]][0]:
            return True
        else:
            return False
        

        
            
def fasta_to_dict(path):
    dict_fa = {}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            temp_list = []
            EnsID = title.split('|')[0]
            chro = title.split('|')[1]
            start = title.split('|')[2]
            end = title.split('|')[3]
            temp_list=[chro,start,end,seq]
            dict_fa[EnsID] = temp_list
    return dict_fa
        
def query_from_dict_fa(dict_fa,abs_start,abs_end,EnsID,strand):
    #print(strand)
    if strand == '+':
        
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        #print(type(start_index))
        end_index = int(abs_end) - start + 1 + 2000
        exon_seq = seq[start_index:end_index]
    
    elif strand == '-':
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        #print(type(start_index))
        end_index = int(abs_end) - start + 1 + 2000 # python range/slice doesn't include end point
        exon_seq_1 = seq[start_index:end_index]
        s = Seq(exon_seq_1,generic_dna)
        exon_seq = str(s.reverse_complement())
    return exon_seq
    
    
           


def exonCoords_to_dict(path,delimiter):
    coords=[]
    dict_exonCoords={}
    with open(path,'r') as file:
        next(file)
        for line in file:
            dict_temp={}
            items = line.split('\t')
            coords=(items[2],items[3],items[4],items[5])
            if items[0] in dict_exonCoords:
                dict_exonCoords[items[0]][items[1]] = coords
            else:
                dict_exonCoords[items[0]] = {}
                dict_exonCoords[items[0]][items[1]] = coords
    # final structure {'EnsID':{E1:[chr,strand,start,end],E2:[chr,strand,start,end]}}
    return dict_exonCoords

            
def find_longest_AA(listAA):
    max=0
    for item in listAA:
        try:
            stop_pos = item.index('*') # return only first occurence
            length = len(item[:stop_pos])
        except ValueError:
            length=len(item)
        if int(length) > max:
            max = int(length)
            max_item = item
    return max_item
    
    
if __name__ == "__main__":
    df_ori = pd.read_csv('/Users/ligk2e/Desktop/df_increase.txt',sep='\t')
    df_exonlist = pd.read_csv('/Users/ligk2e/Desktop/project/mRNA-ExonIDs.txt',sep='\t',
                              header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    dict_exonCoords = debug.debug_boundary_issue('/Users/ligk2e/Desktop/project/Hs_Ensembl_exon.txt')
    dict_fa = fasta_to_dict('/Users/ligk2e/Desktop/project/Hs_gene-seq-2000_flank.fa')
    col1,col2 = match_with_exonlist(df_ori,df_exonlist,dict_exonCoords)
    with open('col1_i.p','wb') as col1_file:
        pickle.dump(col1,col1_file)
    with open('col2_i.p','wb') as col2_file:
        pickle.dump(col2,col2_file)
        
    # load pickle:
    # with open('col1.p','rb') as col1_file:
    #     col1 = pickle.load(col1_file)
    #df_ori.to_csv('/Users/ligk2e/Desktop/result_branch2.txt',sep='\t',header=True,index=False)

# col1[0][0][1]
# col1[0][1][1]    
