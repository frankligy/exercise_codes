#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 11:43:54 2019

@author: LiGk2e
"""

import os

os.chdir('/Users/ligk2e/Downloads/')

ensembl_genes = []
with open('genelist.txt','r') as f1:
    for line in f1:
        ensembl_genes.append(line.strip('\n'))

print(len(ensembl_genes)) 
       
dict = {}   
max_tx = {}


     
f2 = open('en_export.txt').readlines()[1:]
for gene in ensembl_genes:
    for line in f2:
        items = line.strip('\n').split('\t')
        if items[0] == gene:
            if gene in dict:
                dict[gene].append(items[-2])
            else:
                dict[gene]=[]
                dict[gene].append(items[-2])
        else:
            length=len(dict[gene])
            del f2[0:length]
            dict[gene] = list(map(int,dict[gene]))
            max_tx[gene] = max(dict[gene])
            
            break
    
                
with open('max_tx.txt','a+') as f3:
    for key,value in max_tx.items():
        f3.write('{}\t{}\n'.format(key,value))
                 
        
        
    

   

        
    
                
                
            
        