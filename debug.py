#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:52:07 2020

@author: ligk2e
"""

import os
os.chdir('/Users/ligk2e/Desktop')
import re
import pandas as pd



def debug_boundary_issue(path):
    ### construct the special dictionary
    dict = {}
    with open(path,'r') as f1:
        next(f1)
        for line in f1:
            items = line.split('\t')
            exon_num = items[1].split('.')[0]    # E1.1 will become E1, I2.1 will be I2
            attr = [items[2],items[3],items[4],items[5]]   # they are chr, strand, start, end respectively
            if items[0] in dict:    # items[0] is EnsID
                if exon_num in dict[items[0]]:
                    if not items[1] in dict[items[0]][exon_num]:
                        dict[items[0]][exon_num][items[1]] = attr                    
                else:
                    dict[items[0]][exon_num] = {}
                    if not items[1] in dict[items[0]][exon_num]:
                        dict[items[0]][exon_num][items[1]] = attr
                           
            else:
                dict[items[0]] = {}            
                if exon_num in dict[items[0]]:
                    if not items[1] in dict[items[0]][exon_num]:
                        dict[items[0]][exon_num][items[1]] = attr
                else:
                    dict[items[0]][exon_num] = {}
                    if not items[1] in dict[items[0]][exon_num]:
                        dict[items[0]][exon_num][items[1]] = attr
                        
        # this will be a three layer lattice structure:
        # EnsID > E1 > E1.1 -- [attr], three keys.
                        
    ### make the change to starting position and return the dictionary as original format, exactly sneak into without any trace
    bigger_new_dict = {}
    for gene in dict.keys():
        new_dict = {}   # as we made the change, we are simutaneously constructing the new dictionary
        for exon in dict[gene].keys():        
            sub_exon_array = list(dict[gene][exon].keys())        
            for sub_exon in dict[gene][exon].keys():
                attr = dict[gene][exon][sub_exon]
                strand = attr[1]
                if strand == '+':
                    if sub_exon == sub_exon_array[-1]:
                        new_dict[sub_exon] = dict[gene][exon][sub_exon]
                    else:
                        end_point = dict[gene][exon][sub_exon][3]
                        new_end_point = int(end_point) - 1 
                        dict[gene][exon][sub_exon][3] = str(new_end_point)
                        new_dict[sub_exon] = dict[gene][exon][sub_exon]
                else: # strand == '-'
                    if sub_exon == sub_exon_array[0]:
                        new_dict[sub_exon] = dict[gene][exon][sub_exon]
                    else:
                        end_point = dict[gene][exon][sub_exon][3]
                        new_end_point = int(end_point) - 1 
                        dict[gene][exon][sub_exon][3] = str(new_end_point)
                        new_dict[sub_exon] = dict[gene][exon][sub_exon]
                    
        bigger_new_dict[gene] = new_dict
    return bigger_new_dict

if __name__ == "__main__":
    dict_exonCoords = debug_boundary_issue('/Users/ligk2e/Desktop/project/Hs_Ensembl_exon.txt')
    

                
                
            
            
            