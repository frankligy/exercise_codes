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


from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq


# record all mutant pos regarding the dna sequence
normal = {}
with open('Non_mutant_L1Hs_consensus.fa','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        normal[title] = seq
normal = normal['L1PA1']

mutant = {}
with open('L1Hs_Mutated_Correct_Frequencies.fa','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        mutant[title] = seq
mutant = mutant['L1PA1']

mutant_pos = []    # 0-based
for i in range(len(mutant)):
    if normal[i] == 'C' and mutant[i] == 'T':
        mutant_pos.append(i)


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


data = []
for k,v in normal_orf.items():
    start = int(k.split('lcl|L1PA1:')[1].split('-')[0]) - 1  
    effective_length = len(v) - 3   # remove the included stop codon
    for i in range(0,effective_length - span + 1,3):
        cds_start = start + i   
        cds_end = start + i + span - 1   
        cds = normal[cds_start:cds_end+1]
        pep = str(Seq(cds).translate(to_stop=False))
        assert '*' not in pep
        # grab which mutant position fall into each peptide
        start_pos = bisect.bisect_left(mutant_pos,cds_start)    
        end_pos = bisect.bisect_right(mutant_pos,cds_end)
        if start_pos < end_pos:
            included_mutant_pos = mutant_pos[start_pos:end_pos]
        else:
            included_mutant_pos = []
        # determine which position of the peptide will be mutated
        for pos in included_mutant_pos:
            p = pos - cds_start
            pep_p = p // 3 + 1
            mutated_cds = deepcopy(cds)
            mutated_cds[p] == 'T'
            mutated_pep = str(Seq(mutated_cds).translate(to_stop=False))
            if pep == 'DCGGVGGGG':
                print(i,cds_start,cds_end,cds,start_pos,end_pos,included_mutant_pos,p,pep_p);sys.exit('stop')
            if not '*' in mutated_pep:
                if pep_p == 2 or pep_p == 9:
                    identity = 'group2'
                    data.append((pep,mutated_pep,identity))
                else:
                    identity = 'group1'
                    data.append((pep,mutated_pep,identity))

final = pd.DataFrame.from_records(data,columns=['pep','mutated_pep','identity'])
final['differ'] = [False if item1 == item2 else True for item1, item2 in zip(final['pep'],final['mutated_pep'])]
final = final.loc[final['differ'],:]

'''
           pep mutated_pep identity
0    FDELREEGF   FDELREEGF   group2
1    DELREEGFR   DELREEGFR   group1
2    ELREEGFRR   ELREEGFRR   group1
3    LREEGFRRS   LREEGFRRS   group1
4    REEGFRRSN   REEGFRRSN   group1
..         ...         ...      ...
806  DCGGVGGGG   DCGGVGGGG   group1
807  CGGVGGGGR   CGGVGGGGR   group1
808  GGVGGGGRD   GGVGGGGRD   group1
809  GVGGGGRDS   GVGGGGRDS   group2
810  VGGGGRDSI   VGGGGRDSI   group1
'''


        
        









