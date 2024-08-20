#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import os,sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
from io import StringIO
import subprocess
import bisect
from copy import deepcopy
import random 


from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq


# record all mutant pos regarding the dna sequence
normal = {}
with open('Non_mutant_L1Hs_consensus.fa','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        normal[title] = seq
normal = normal['L1PA1']

# run orffinder, make sure configure
# export LD_LIBRARY_PATH=/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/lib:$LD_LIBRARY_PATH
# NCBI_ORFFINDER_PATH = '/gpfs/data/yarmarkovichlab/chordoma/NeoVerse_analysis/ORFfinder'
# subprocess.run("{} -in {} -s 0 -n true -strand plus -out {} -outfmt 1".format(NCBI_ORFFINDER_PATH,'Non_mutant_L1Hs_consensus.fa','Non_mutant_L1Hs_consensus_orf.fa'),shell=True)

# in-silico translation and record whether a peptide is group 1 or group 2
normal_orf = {}
with open('Non_mutant_L1Hs_consensus_orf.fa','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        normal_orf[title] = seq

mer = 9
span = 3 * mer

def mutate_codon(codon):
    dic = {
        'TCA':('TTA',0.523),
        'TCT':('TTT',0.30),
        'TCC':('TTC',0.10),
        'TCG':('TTG',0.05)
    }

    subs = set(dic.keys())
    if codon not in subs:
        return codon
    else:
        replace, prob = dic[codon]
        if random.random() < prob:
            return replace
        else:
            return codon


group1_pos = set([1,3,4,5,6,7,8])
group2_pos = set([2,9])
data = []
for k,v in normal_orf.items():
    start = int(k.split('lcl|L1PA1:')[1].split('-')[0]) - 1  
    effective_length = len(v) - 3   # remove the included stop codon
    for i in range(0,effective_length - span + 1,3):
        cds_start = start + i   
        cds_end = start + i + span - 1   
        cds = normal[cds_start:cds_end+1]

        mutated_cds = ''
        mutated_pos = []
        for i_, j in enumerate(range(0,len(cds)-3+1,3)):
            codon = cds[j:j+3]
            mutated_codon = mutate_codon(codon)
            mutated_cds += mutated_codon
            if mutated_codon != codon:
                mutated_pos.append(i_ + 1)
            
        pep = str(Seq(cds).translate(to_stop=False))
        mutated_pep = str(Seq(mutated_cds).translate(to_stop=False))
        assert '*' not in pep
        assert '*' not in mutated_pep


        if len(mutated_pos) > 0:
            group1 = False
            group2 = False
            if len(group1_pos.intersection(set(mutated_pos))) > 0:
                group1 = True
            if len(group2_pos.intersection(set(mutated_pos))) > 0:
                group2 = True
            if group1 and group2:
                identity = 'both'
            elif group1 and not group2:
                identity = 'group1'
            elif group2 and not group1:
                identity = 'group2'
            data.append((pep,mutated_pep,identity))
        

final = pd.DataFrame.from_records(data,columns=['pep','mutated_pep','identity'])
final.to_csv('final.txt',sep='\t',index=None)


# use run_netmhcpan.py to get binding results for A0201 and A2402


        
        









