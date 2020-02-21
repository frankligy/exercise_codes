#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 20:27:50 2020

@author: ligk2e
"""

import os
os.chdir("/Users/ligk2e/Desktop")
import pickle
import pandas as pd



    

def score_GC(sequence):
    GC_content = 0
    length_seq = len(sequence)
    for nt in sequence:
        if nt == 'G' or nt == 'C':
            GC_content += 1
    GC_percent = GC_content / length_seq
    return GC_percent
            
def score_coding_bias(sequence):
    # coding frequency table is from GenScript webpage
    usage_dict = {'TTT':16.9,'TTC':20.4,'TTA':7.2,'TTG':12.6,'TAT':12.0,'TAC':15.6,'TAA':0.7,'TAG':0.5,
                  'CTT':12.8,'CTC':19.4,'CTA':6.9,'CTG':40.3,'CAT':10.4,'CAC':14.9,'CAA':11.8,'CAG':34.6,
                  'ATT':15.7,'ATC':21.4,'ATA':7.1,'ATG':22.3,'AAT':16.7,'AAC':19.5,'AAA':24.0,'AAG':32.9,
                  'GTT':10.9,'GTC':14.6,'GTA':7.0,'GTG':28.9,'GAT':22.3,'GAC':26.0,'GAA':29.0,'GAG':40.8,
                  'TCT':14.6,'TCC':17.4,'TCA':11.7,'TCG':4.5,'TGT':9.9,'TGC':12.2,'TGA':1.3,'TGG':12.8,
                  'CCT':17.3,'CCC':20.0,'CCA':16.7,'CCG':7.0,'CGT':4.7,'CGC':10.9,'CGA':6.3,'CGG':11.9,
                  'ACT':12.8,'ACC':19.2,'ACA':14.8,'ACG':6.2,'AGT':11.9,'AGC':19.4,'AGA':11.5,'AGG':11.4,
                  'GCT':18.6,'GCC':28.5,'GCA':16.0,'GCG':7.6,'GGT':10.8,'GGC':22.8,'GGA':16.3,'GGG':16.4} 
    length_seq = len(sequence)
    i = 0   
    score = 0
    while i < length_seq - 2:
        triplet = sequence[i:i+3:1]
        score_tri = usage_dict[triplet]
        score += score_tri
        i += 3
    score_bias = 1/score * 100 # take reciprocal and scale them by multipling 100
    return score_bias
    


def readingframe2peptide(frame_dict):
    frag_comp_array = []
    for pep_seq in frame_dict.values():
        frag_array1 = pep_seq.split('*')
        for frag in frag_array1:
            if 'M' not in frag or len(frag) == 0:
                continue
            else:
                index_M = frag.index('M')
                frag_comp = frag[index_M:]
                frag_comp_array.append(frag_comp)
    #print(frag_comp_array)
    #pick most likely one: length, GC content, coding frequency
    max_seq = ''
    max_length = 0
    max_item_score = 0
    #global count
    for item in frag_comp_array:
        temp1 = len(item)
        add_score = score_GC(item)
        if temp1 > max_length:
            max_length = temp1
            max_item_score = add_score
            max_seq = item
        elif temp1 == max_length:
            if add_score > max_item_score:
                max_length = temp1
                max_item_score = add_score
                max_seq = item
            elif add_score == max_item_score:
                #count += 1
                print('Even considering GC and coding frequency are not able to differentiate them')            
    return max_seq

def final_conversion(col_pickle_file):
    with open(col_pickle_file,'rb') as col_file:
        col = pickle.load(col_file)
    output_array = []
    for event in col:
        temp_array = []
        for transcript in event[1]:
            if transcript == '':
                temp_array.append(transcript)
            else:
                max_pep = readingframe2peptide(transcript) 
                temp_array.append(max_pep)
        output_array.append(temp_array)
    return output_array
               
    

if __name__ == "__main__":
   #count =0
   df_ori = pd.read_csv('/Users/ligk2e/Desktop/df_increase.txt',sep='\t')
   output_exam = final_conversion('/Users/ligk2e/Desktop/col1_i.p')
   output_back = final_conversion('/Users/ligk2e/Desktop/col2_i.p')
   df_ori['exam_match'] = output_exam
   df_ori['back_match'] = output_back
   #print(count)
   df_ori.to_csv('/Users/ligk2e/Desktop/match_i_v1_Feb20.txt',sep='\t',header=True,index=False)
   output_decrease = [output_exam,output_back]
   with open('increase.p','wb') as f1:
       pickle.dump(output_decrease,f1)
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    