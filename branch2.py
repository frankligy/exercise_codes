#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 14:05:41 2020

@author: ligk2e
"""
import pandas as pd
import task1_mod as mich
import os

def match_with_exonlist(df_ori,df_exonlist,dict_exonCoords):
    
    for i in range(df_ori.shape[0]):
        amino_acid_seq=[]
        temp=mich.UID(df_ori,i)
        EnsID=list(temp.keys())[0].split(':')[1]
        Exons=list(temp.values())[0][1].split('-')[0] + '|' + list(temp.values())[0][1].split('-')[1]
        try:
            df_certain = df_exonlist[df_exonlist['EnsGID'] == EnsID]
        except:
            #print('{} can not match up with anything'.format(EnsID))
            continue
        for item in list(df_certain['Exons']):
            full_transcript=''
            fullAA=[]
            if Exons in item:
                Exonlist = item.split('|')
                for i in range(len(Exonlist)):
                    coords = dict_exonCoords[Exonlist[i]]
                    path='http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=' + coords[0] + ':' + coords[1] + ',' + coords[2]
                    frag=mich.web_scraping(path)
                    full_transcript += frag
                pot_fullAA=mich.translate(full_transcript)
                max_fullAA=find_longest_AA(list(pot_fullAA.values()))
                fullAA.append(max_fullAA)
        if fullAA:
            final_fullAA = find_longest_AA(fullAA)
            amino_acid_seq.append(final_fullAA)
        else:
            amino_acid_seq.append(fullAA)
    df_ori['fullAA'] = amino_acid_seq
    return df_ori      
            

def exonCoords_to_dict(path,delimiter):
    coords=[]
    dict_exonCoords={}
    with open(path,'r') as file:
        for line in file:
            items = line.split('\t')
            coords=(items[2],items[4],items[5])
            dict_exonCoords[items[1]]=coords
    return dict_exonCoords
            
def find_longest_AA(listAA):
    max=0
    for item in listAA:
        try:
            stop_pos = item.index('*') # return only first occurence
            length = item[:stop_pos-1]
        except ValueError:
            length=len(item)
            if length > max:
                max = length
                max_item = item
    return max_item
    
    
if __name__ == "__main__":
    df_ori = pd.read_csv('df_decrease.txt',sep='\t')
    df_exonlist = pd.read_csv('project/mRNA-ExonIDs.txt',sep='\t')
    dict_exonCoords = exonCoords_to_dict('project/Hs_Ensembl_exon.txt','\t')
    result = match_with_exonlist(df_ori,df_exonlist,dict_exonCoords)
    result.to_csv('ban2.txt',sep='\t',header=True,index=False)


    
