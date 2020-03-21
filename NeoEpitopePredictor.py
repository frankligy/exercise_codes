#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 19:43:43 2020

@author: ligk2e
"""

import os
os.chdir('/Users/ligk2e/Desktop/project_breast/RBM47/')
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
from decimal import Decimal as D
#import pickle
import regex
import re
from time import process_time
import collections
#import ast
import matplotlib.pyplot as plt
import math

############################################################################################
# part1: branch2.py, match to existing peptides
#############################################################################################

def RatingFunction(x):   # dPSI is float64 type
    if x > 0:
        return True  # anything like good/bad, because you still need to use a == good to slice the datafrome
    else:
        return False # or use lambda function

def GetIncreasedPart(df):
    df['sign'] = df['dPSI'].apply(lambda x: True if x>0 else False) # x = lambda a,b:a+b; x(5,6)
    # Method2: apply(RatingFunction) 
    # Method3: list comprehension, see main function example when got extracellular instances    
    df_ori = df[df['sign']==True]
    df_ori = df_ori.drop(columns=['sign'])  # how to drop a column
    return df_ori




def UID(df, i):
    uid = list(df['UID'])[i]
    #uid = uid.split(':')[1]
    gene = uid.split(':')[0]
    dict = {}
    gene = gene + ':' + uid.split('|')[1].split(':')[0]
    x = uid.split('|')
    dict[gene] = [x[0].split(':')[2]]
    dict[gene].append(x[1].split(':')[1])
    #print ((uid.split('|')[0],gene + ':' + uid.split('|')[1]))
    #return ((uid.split('|')[0],gene + ':' + uid.split('|')[1]))
    #print dict
    #{'gene:ENSid':[E22-33,E34-56]}
    return dict

def match_with_exonlist(df_ori,df_exonlist,dict_exonCoords):
#    sum = 0
    col1 = []
    col2 = []
    
    for i in range(df_ori.shape[0]):
        temp=UID(df_ori,i)
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
#        final_fullAA = []
#        peek_pep = []
        full_transcript_store = []
        
#    final_fullAA = []
#    peek_pep = []
    full_transcript_store = []
    for item in list(df_certain['Exons']):
        full_transcript=''
#        splicing_transcript=''
        if Exons in item:
            Exonlist = item.split('|')
            for j in range(len(Exonlist)):
                coords = dict_exonCoords[EnsID][Exonlist[j]]
                strand = coords[1]
                judge = check_exonlist_general(Exonlist,j,strand)
                if strand == '+' and judge:   
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1]) # corresponds to abs_start, abs_end, strand
                elif strand == '+' and not judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],int(coords[3])-1,EnsID,coords[1]) 
                elif strand == '-' and judge:
                    frag = query_from_dict_fa(dict_fa,coords[2],coords[3],EnsID,coords[1])
                elif strand == '-' and not judge:
                    frag = query_from_dict_fa(dict_fa,int(coords[2])+1,coords[3],EnsID,coords[1])  # because of the weird
                    # expression of minus strand, need to draw an illustrator to visulize that.
                full_transcript += frag
            full_transcript = full_transcript.replace('\n','')
#            splicing_transcript = full_transcript.replace('\n','')
            full_transcript_store.append(full_transcript)   
#            pot_fullAA=mich.translate(full_transcript)
#            peek_pep.append(pot_fullAA)
#            max_fullAA=find_longest_AA(list(pot_fullAA.values()))
#            final_fullAA.append(max_fullAA)
        else:
#            peek_pep.append('')
#            final_fullAA.append('')
            full_transcript_store.append('')
#    result = [final_fullAA,peek_pep,full_transcript_store]
    return full_transcript_store


def check_exonlist_general(exonlist,index,strand):
    dict = {}
    for subexon in exonlist:
        exon_num = subexon.split('.')[0]
        subexon_num = subexon.split('.')[1]
        if exon_num in dict:
            dict[exon_num].append(subexon_num)
        else:
            dict[exon_num] = []
            dict[exon_num].append(subexon_num)  # E14 > 1,2,4,5
    # check
    query = exonlist[index]
    query_exon_num = query.split('.')[0]   #E14.1
    query_subexon_num = int(query.split('.')[1])   #2, it is a int
    if strand == '+':
        if str(query_subexon_num + 1) in dict[query_exon_num]:
            return False
        else:
            return True
    else:
        if str(query_subexon_num + 1) in dict[query_exon_num]:
            return False
        else:
            return True
        
        
        
            
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
        seq_reverse = dict_fa[EnsID][3]
        seq_forward = str(Seq(seq_reverse,generic_dna).reverse_complement())  # Hs_gene.fa restore the reverse strand info
        start_index = int(abs_start) - start + 2000
        #print(type(start_index))
        end_index = int(abs_end) - start + 1 + 2000 # python range/slice doesn't include end point
        exon_seq_1 = seq_forward[start_index:end_index]
        s = Seq(exon_seq_1,generic_dna)
        exon_seq = str(s.reverse_complement())
    return exon_seq

##############################################################################################    
# part2: pick_peptide.py, find the most likely ORF

###############################################################################################    

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
    # do a normaliztion for each triplet, then for all the triplet's sum, divided by the number of triplet
    min_freq = 4.5
    max_freq = 40.8
    norm_usage_dict = {}
    for codon,freq in usage_dict.items():
        norm_usage_dict[codon] = float((D(freq) - D(min_freq)) / (D(max_freq) - D(min_freq)))        
    length_seq = len(sequence)
    num_triplet = length_seq/3
    i = 0   
    score = 0
    while i < length_seq - 2:
        triplet = sequence[i:i+3:1]
        score_tri = norm_usage_dict[triplet]
        score += score_tri
        i += 3
    score_bias = score/num_triplet # take reciprocal and scale them by multipling 100
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

def transcript2peptide(cdna_sequence):
    # TAA,TGA,TAG
    # tips: find/index and split function && their corresponding re version for multiple stirng
    reading_manners = []
    reading_manners.append(cdna_sequence[0:])
    reading_manners.append(cdna_sequence[1:])
    reading_manners.append(cdna_sequence[2:])
    frag_comp_array = []
    for manner in reading_manners:       
#        frag_array = re.split(r'TAA|TGA|TAG',manner), it is for multiple condition
        pos = []
        for m in re.finditer(r'TAA|TGA|TAG',manner):   # for multiple instances
            if m.start() % 3 == 0:
                pos.append(m.start())
        frag_array = pos_to_frags(pos,manner)
        for frag in frag_array:
            if 'ATG' not in frag or len(frag) == 0:
                continue
            else:
                for n in re.finditer('ATG',frag):
                    if (len(frag) - n.start()) % 3 == 0:
                        frag_comp = frag[n.start():]
                        frag_comp_array.append(frag_comp)
                        break
                    else:
                        continue
                    
    #######################  # We think if you only has longer length(0-7) but add_score is not higher than original one, you are FAlSE
    max_seq = ''
    max_length = 0
    max_item_score = 0
#    global count
    for item in frag_comp_array:
        temp1 = len(item)
        add_score = score_GC(item) + score_coding_bias(item)
        if (temp1 - max_length) >= 8:
            max_length = temp1
            max_item_score = add_score
            max_seq = item
        elif (temp1 - max_length) >= 0 and (temp1 - max_length) < 8:
            if add_score >= max_item_score:
                max_length = temp1
                max_item_score = add_score
                max_seq = item
            else:
#                count += 1
                print('equal length but less likely to be a true ORF or longer length but less likely to be a true ORF',
                      add_score,max_item_score) 
    max_seq_tran = max_seq
    max_seq_aa = str(Seq(max_seq,generic_dna).translate(to_stop=False))
    max_seq_final = [max_seq_tran,max_seq_aa,frag_comp_array]
    return max_seq_final     

def pos_to_frags(pos,sequence):
    frag_array = []
    if pos: #tips: this elegant way to define criteria        
        frag_array.append(sequence[0:pos[0]])
        i = 0
        while i < len(pos)-1:
            frag_array.append(sequence[pos[i]+3:pos[i+1]])
            i += 1
        last_seq = sequence[pos[-1]+3:]
        if not any(codon in last_seq for codon in ['TAA','TAG','TGA']):
            frag_array.append(sequence[pos[-1]+3:])
#    else:
#        pass
    return frag_array
        
    
def final_conversion(col):
#    with open(col_pickle_file,'rb') as col_file:
#        col = pickle.load(col_file)
    output_array_aa = []
    output_array_tran = []
#    output_peek = []
    for event in col:
        temp_array_aa = []
        temp_array_tran = []
#        temp_peek = []
        for transcript in event:
            if transcript == '':
                temp_array_aa.append(transcript)
                temp_array_tran.append(transcript)
#                temp_peek.append(transcript)
            else:
                max_pep = transcript2peptide(transcript)[1]
                max_tran = transcript2peptide(transcript)[0]
#                peek = transcript2peptide(transcript)[2]
                temp_array_aa.append(max_pep)
                temp_array_tran.append(max_tran)
#                temp_peek.append(peek)
        output_array_aa.append(temp_array_aa)
        output_array_tran.append(temp_array_tran)
#        output_peek.append(temp_peek)
        output_array = [output_array_aa,output_array_tran]
    return output_array           


def exonCoords_to_dict(path,delimiter):
    coords=[]
    dict_exonCoords={}
    with open(path,'r') as file:
        next(file)
        for line in file:
#            dict_temp={}
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

#################################################################################################
# part3: following.py   find the junction sites' sequence.   
#################################################################################################

def retrieve_junction_site(df_ori):
    exam_seq,back_seq = [],[]
    for i in range(df_ori.shape[0]):
        temp = UID(df_ori,i)
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
        exon_seq = query_from_dict_fa(dict_fa,attrs[2],attrs[3],EnsID,attrs[1])  
    except:
#        print(subexon,'t')
        try:
            suffix = subexon.split('_')[1]
        except:
            exon_seq = '*'   # '*' means possible gene fusion
            print(subexon,'possible fusion gene event',EnsID)   # normal subexon without cognate coordinate in dict will fall here
        else:
            subexon = subexon.split('_')[0]
            try:
                attrs = dict_exonCoords[EnsID][subexon]
            except:
                exon_seq = '#'   # '# means 'U0.1_32130925'
                print(subexon,'splicing occurs in UTR',EnsID)
            else:
                if flag == 'site2':           
                    exon_seq = query_from_dict_fa(dict_fa,suffix,attrs[3],EnsID,attrs[1])  # chr,strand, start,end
                elif flag == 'site1':
                    exon_seq = query_from_dict_fa(dict_fa,attrs[2],suffix,EnsID,attrs[1])
    return exon_seq


##################################################################################################
# part4: narrow down to extracellular instances and check if good representative and do seeding alignment   
###################################################################################################
    
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
   
    exam_col = col1
    truth_table = [True if item in extracellular_gene else False for item in EnsID]
    narrow_whole_tran = []
    for i in range(len(truth_table)):
        if truth_table[i]:
            narrow_whole_tran.append(exam_col[i])
    # one on one to find correponding representative and their whole transcript        
    representative = []
    whole_tran = []
    position_array = []
    for i in range(df.shape[0]):
        event = list(df['exam_match_tran'])[i]
#        event = ast.literal_eval(event)
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
                add_score = score_GC(tran) + score_coding_bias(tran)
                if (temp1 - max_length) >= 8:
                    max_length = temp1
                    max_item_score = add_score
                    max_seq = tran
                    position = flag
                elif (temp1 - max_length) >= 0 and (temp1 - max_length) < 8:
                    if add_score >= max_item_score:
                        max_length = temp1
                        max_item_score = add_score
                        max_seq = tran
                        position =flag
        try:
            representative.append(max_seq) 
            whole_tran.append(narrow_whole_tran[i][position])
            position_array.append(position)
        except:
            representative.append('')
            whole_tran.append('')
            position_array.append(-1)
    
    return representative,whole_tran,position_array
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
            
    # we have to apply fuzzy matching, because junction consists of former part and latter part(i.e. E6.3|E8.1)
    # In my way to handle overlapping 1nt, the former one will always stay constant, but latter one in the whole
    # transcript, it might get trimmed but here we don't trim it, so there might be 1 overhang in jucntion seq.
    
            pattern = regex.compile('(%s){d<=1}' % junction) 
            start_junction = pattern.search(whole).span()[0]
            end_junction = pattern.search(whole).span()[1] - 1
            
#            start_junction = whole.find(junction)
#            end_junction = whole.find(junction) + len(junction)
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
        result = 'not aligned'
    elif sum(notes) == len(notes):
        result = 'totally aligned'
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

##############################################################################################
# part5: interrogating chromosome stataistics
#############################################################################################
    

def ChroDistribution(df):
    chro_array = []
    for i in range(df.shape[0]):
        ensid = list(df['UID'])[i].split('|')[0].split(':')[1]
        chro = dict_fa[ensid][0]
        chro_array.append(chro)
    freq = collections.Counter(chro_array)
    return freq

'''courtesy by user on stackoverflow'''
def Round2Precision(value,precision:int=0,mode:str=''): # default argument with specified type
    assert precision >= 0 # if true, continue, otherwise raise assertError, using for self-check
    value *= 10 ** precision # if you wanna round by precision, have to do that
    method = round   # round will base on >.5 or <.5
    if mode.lower() == 'up': 
        method = math.ceil     # always round up
    elif mode.lower() == 'down':
        method = math.floor   # always round down
    answer = '{0:.{1}f}'.format(method(value)/10**precision,precision)   
    return float(answer) 


def PlotChroScarse(chro_dict,path):
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.75,0.75])  #[left,bottom, width, height]
    scarse = {}
    for chro,attr in chro_dict.items():
        chro_s = re.split(r'chr',chro)[-1]
        scarse[chro_s] = Round2Precision(attr[0]/attr[1],2)  # genes per 1 Millon bp
    x_axis = list(scarse.keys())
    y_axis = list(scarse.values())
    ax.bar(x_axis,y_axis)
    ax.set(xlabel='chromosome',ylabel='genes per 1 Million bp',title='crowdness of human chromosome')
    
    #ax.legend()    
    #plt.show()
    fig.savefig(path)
    plt.close(fig)
    
        
    

if __name__ == "__main__":
    start_time = process_time()
    # get increased part
    df = pd.read_csv('PSI.RBM47Deletion_vs_noRBM47Deletion.txt',sep='\t')
    
    # load the files, return matched result and cognate whole transcript sequence
    df_ori = GetIncreasedPart(df)
    df_exonlist = pd.read_csv('mRNA-ExonIDs.txt',sep='\t',
                              header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
    dict_exonCoords = exonCoords_to_dict('Hs_Ensembl_exon.txt','\t')
    dict_fa = fasta_to_dict('Hs_gene-seq-2000_flank.fa')
    col1,col2 = match_with_exonlist(df_ori,df_exonlist,dict_exonCoords)
    
    # derive the most likely ORF for each whole transcript sequence
    output_exam = final_conversion(col1)
    output_back = final_conversion(col2)
    output_exam_aa,output_exam_tran = output_exam[0],output_exam[1]
    output_back_aa,output_back_tran = output_back[0],output_back[1]

    df_ori['exam_match_aa'] = output_exam_aa
    df_ori['back_match_aa'] = output_back_aa
    df_ori['exam_match_tran'] = output_exam_tran
    df_ori['back_match_tran'] = output_back_tran
    
    # derive the junction site sequence and add two columns to df_ori
    new_df = retrieve_junction_site(df_ori)
    
    # get extracellur, part4
    EnsID = extract_EnsID(df_ori)
    write_list_to_file(EnsID)  # write all EnsID to a list then batch query on uniprot to get corresponding uniprot ID
    
    
    
    # Pause here to get EnsID query list to upload to Uniprot to get query_result
    '''
    go to uniprot, retrieve/mapping, upload the file we just got, choose from Ensembl to uniprotKB
    '''
    #
    #########################################################################################
    ########################################################################################
    
    
    EnsID_to_uniprot = pd.read_csv('query_result.tab',sep='\t',header=None,
       names=['EnsID','isoforms','uniprot_entry','Entry_name','protein','length','topology','gene_name'], #change the column name if needed
       skiprows=1)     # load in the query result, EnsID to Uniprot ID 
    # aboved step could use API, https://www.uniprot.org/help/uploadlists
    
    # based on EnsID to Uniprot ID relationship, narrow the df_ori to df_ori_narrow which only has extracellular one
    uniprot_info = pd.read_csv('uniprot_info.tab',sep='\t')  # load all uniprot human membrane protein info
    extracellular_gene,dict_EnsID_uni = find_membrane_EnsID(EnsID_to_uniprot,uniprot_info)
    df_ori['condition'] = [True if item in extracellular_gene else False for item in EnsID]
    # this could also be achived by apply function or lambda function
    df_ori_narrow = df_ori[df_ori['condition'] == True]
    # aboved step could use API, https://www.uniprot.org/help/uploadlists
    
   #### get the representative and check if it is a good representative

    representative,whole_tran,position_array = representative_tran_and_whole_tran(df_ori_narrow)
    df_ori_narrow['representative_tran'] = representative
    df_ori_narrow['whole_tran'] = whole_tran# don't to_csv then read_csv
    df_ori_narrow['postion'] = position_array
    df_ori_narrow_good_repre = check_if_good_representative(df_ori_narrow)
    
    ### final alignment
    dict_uni_fa = read_uniprot_seq('uniprot_canonical.fasta')
    comment,alignment,repre_aa = alignment_to_uniprot(df_ori_narrow_good_repre,dict_uni_fa,dict_EnsID_uni)   
    df_ori_narrow_good_repre['comment'] = comment
    df_ori_narrow_good_repre['alignment'] = alignment
    
               
    # mannaully check the case
    df_ori_narrow_good_repre['repre_aa'] = repre_aa  # this operation might catch a caveat, .loc[rowInd,colInd] = value
    
    #write them out
    df_ori_narrow_good_repre.to_csv('final_result.txt',sep='\t',header=True,index=False)      
    df_ori_narrow.to_csv('extracellur_narrow_down.txt',sep='\t',header=True,index=False)
   
    
    # ENSG00000282228   - strand use case
    # ENSG00000110514   +
    # ENSG00000243646

# focus on all novel splicing events
    novel_splicing = df_ori_narrow[df_ori_narrow['whole_tran']=='']
    novel_splicing.to_csv('novel_splicing.txt',sep='\t',header=True,index=False)
    
    
# seperate totally aligned, not aligned and partially aligned instances
    totally = df_ori_narrow_good_repre[df_ori_narrow_good_repre['alignment'] == 'totally aligned']
    partially = df_ori_narrow_good_repre[df_ori_narrow_good_repre['alignment'] == 'partially aligned'] 
    not_aligned = df_ori_narrow_good_repre[df_ori_narrow_good_repre['alignment'] == 'not aligned'] 
    
    partially.to_csv('partially.txt',sep='\t',header=True,index=False)
    not_aligned.to_csv('not_aligned.txt',sep='\t',header=True,index=False)
    
    # summarize the distribution of splicing event in df_all and df_increased
    freq_all = ChroDistribution(df)
    freq_increased = ChroDistribution(df_ori)
#    print(freq_all,freq_increased)
    
    # continue exploit on chromosomes
    chro_dict = {
            'chr1': [1961,248,'Metacentric'],    #[1961genes,248 or so million bp, type of centromere]
            'chr2': [1194,242,'Submetacentric'],
            'chr3': [1024,198,'Metacentric'],
            'chr4': [727,190,'Submetacentric'],
            'chr5': [839,181,'Submetacentric'],
            'chr6': [996,170,'Submetacentric'],
            'chr7': [862,159,'Submetacentric'],
            'chr8': [646,145,'Submetacentric'],
            'chr9': [739,138,'Submetacentric'],
            'chr10': [706,133,'Submetacentric'],
            'chr11': [1224,135,'Submetacentric'],
            'chr12': [988,133,'Submetacentric'],
            'chr13': [308,114,'Acrocentric'],
            'chr14': [583,107,'Acrocentric'],
            'chr15': [561,101,'Acrocentric'],
            'chr16': [795,90,'Metacentric'],
            'chr17': [1124,83,'Submetacentric'],
            'chr18': [261,80,'Submetacentric'],
            'chr19': [1357,58,'Metacentric'],
            'chr20': [516,64,'Metacentric'],
            'chr21': [215,46,'Acrocentric'],
            'chr22': [417,50,'Acrocentric'],
            'chrX': [804,156,'Submetacentric'],
            'chrY': [63,57,'Acrocentric']}  
    PlotChroScarse(chro_dict,'human chromosome genes distribution.pdf')
    
    end_time = process_time()
    print("from {0} to {1},consuming {2} seconds".format(start_time,end_time,end_time - start_time))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    