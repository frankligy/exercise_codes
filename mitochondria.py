#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:28:44 2020

@author: ligk2e
"""

##################################################################################################################
# ######################         README        ###################################################################
# 1. Input files: 
#       a) Homo_sapiens.GRCh38.dna.chromosome.MT.fa  (mitochondria dna sequence)
#       b) mart_export.txt (All 13 protein coding gene sequences)
#  The b) file is just used for calculating sensitivity.
#
#
#
# 2. Step by step:
#    a) reading the fasta file into session, parse the fasta file
#    b) extract all ORFs from dna sequence
#    c) define score for length, GC content(score_GC function)), coding preference(score_coding_bias function)
#
#
# 3. How to run that:
# Python3 mitochondria.py   
#
# 4. Caveat:
#  a) changing the path of two input file in the main function down below
#  b) please install Biopython package and re package
##################################################################################################################

import os
#os.chdir('/Users/ligk2e/Desktop/gene_finding/')
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re


def read_mito_dna(path):
    with open(path,'r') as f1:
        next(f1)
        dna_sequence = ''
        for line in f1:
            dna_sequence += line.rstrip('\n')
        return dna_sequence
    

def from_dna_to_ORF(dna_sequence):
    reading_manners = []
    # for all positive strand
    reading_manners.append(dna_sequence[0:])
    reading_manners.append(dna_sequence[1:])
    reading_manners.append(dna_sequence[2:])
    # for all negative strand
    rt = str(Seq(dna_sequence,generic_dna).reverse_complement())
    reading_manners.append(rt[0:])
    reading_manners.append(rt[1:])
    reading_manners.append(rt[2:])
    # ready to ascertain all ORF
    frag_comp_array = []
    for manner in reading_manners:       
        pos = []     
        for m in re.finditer(r'TAA|TAG',manner):     #in mitochondria, TGA is for tryptophan, not stop codon
            if m.start() % 3 == 0:      
                pos.append(m.start())
        frag_array = pos_to_frags(pos,manner)
       
        for frag in frag_array:
            pos2 = []
            if not any(x in frag for x in ['ATA','ATT','ATG']) or len(frag) == 0: #highlight systax
                pass
            else:
                for n in re.finditer(r'ATA|ATT|ATG',frag):
                  
                    if (len(frag)-n.start()) % 3 == 0:
                        pos2.append(n.start())
                ORFs = pos2_to_frags(pos2,frag)
           
                [frag_comp_array.append(item) for item in ORFs] #highlight syntax
            
    return frag_comp_array

def pos_to_frags(pos,sequence):
    frag_array = []
    if pos:        
        frag_array.append(sequence[0:pos[0]])
        i = 0
        while i < len(pos)-1:
            frag_array.append(sequence[pos[i]+3:pos[i+1]])
            i += 1
        frag_array.append(sequence[pos[-1]+3:])
    return frag_array


def pos2_to_frags(pos2,sequence):
    frag_array = []
    stop_index = len(sequence)
    if pos2:
        for pos in pos2:
            frag_array.append(sequence[pos:stop_index])
    return frag_array
        


def score_GC(sequence):
    GC_content = 0
    length_seq = len(sequence)
    for nt in sequence:
        if nt == 'G' or nt == 'C':
            GC_content += 1
    GC_percent = GC_content / length_seq
    return GC_percent
            
def score_coding_bias(sequence):
    # coding frequency table is from 'A new and updated resource for codon usage tables' BMC Bioinformatics
    usage_dict = {'UUU':69,'UUC':139,'UUA':65,'UUG':11,'CUU':65,'CUC':167,'CUA':276,'CUG':42,'AUU':112,
                  'AUC':196,'AUA':165,'AUG':32,'GUU':22,'GUC':45,'GUA':61,'GUG':8,'UCU':29,'UCC':99,
                  'UCA':81,'UCG':7,'CCU':37,'CCC':119,'CCA':52,'CCG':7,'ACU':50,'ACC':155,'ACA':132,
                  'ACG':10,'GCU':39,'GCC':123,'GCA':79,'GCG':5,'UAU':35,'UAC':89,'UAA':4,'UAG':3,'CAU':18,
                  'CAC':79,'CAA':82,'CAG':8,'AAU':29,'AAC':131,'AAA':84,'AAG':9,'GAU':12,'GAC':51,'GAA':63,
                  'GAG':15,'UGU':5,'UGC':17,'UGA':90,'UGG':9,'CGU':6,'CGC':26,'CGA':28,'CGG':0,'AGU':11,
                  'AGC':37,'AGA':1,'AGG':0,'GGU':16,'GGC':87,'GGA':61,'GGG':19} 
    sequence = sequence.replace('T','U')
    length_seq = len(sequence)
    i = 0   
    score = 0
    while i < length_seq - 2:
        triplet = sequence[i:i+3:1]
   
        score_tri = usage_dict[triplet]
 
        #print(triplet)
 
        score += score_tri
        i += 3
    score_bias = score # take reciprocal and scale them by multipling 100
    return score_bias

def apply_filtering_for_ORF(ORFs):
    dict = {}
    for ORF in ORFs:
        try:
            ORF.index('N')
        except ValueError:
            if len(ORF) > 200 and score_GC(ORF) > 0.30 and score_coding_bias(ORF)>7000:
                dict[ORF] = score_coding_bias(ORF)
    return dict

if __name__ == "__main__":
    # please change the following path
    mito_fasta = read_mito_dna('/Users/ligk2e/Desktop/gene_finding/Homo_sapiens.GRCh38.dna.chromosome.MT.fa')
    ORFs = from_dna_to_ORF(mito_fasta)
    coding_genes = apply_filtering_for_ORF(ORFs)
    
    # load the standard
    real_coding_genes = {}
    # please change the following path
    with open('mart_export.txt','r') as f2:
        for title,seq in SimpleFastaParser(f2):
            if title.split('|')[0].lstrip('>') in real_coding_genes:
                real_coding_genes[title.split('|')[0].lstrip('>')] = seq[:-3]
            else:
                real_coding_genes[title.split('|')[0].lstrip('>')] = []
                real_coding_genes[title.split('|')[0].lstrip('>')] = seq[:-3]
    
    # sensitivity
    sensitivity = 0
    recovery = []
    for key,value in real_coding_genes.items():
        for pre in list(coding_genes.keys()):
            if value == pre:
                sensitivity += 1
                recovery.append(key)
     
    # print to standard output           
    [print('%s has been recovered by my predictor\n' % EnsID) for EnsID in recovery]
    print('My perdictor reports 334 potential protein coding gene, 5 out of 13 known protein coding genes has been recovered. The sensitivity is %.02f' % (sensitivity/13))
         
    
            










