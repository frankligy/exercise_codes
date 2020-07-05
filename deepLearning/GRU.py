#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 11:34:51 2020

@author: ligk2e
"""

import pandas as pd
import xlrd
import numpy as np
import json
import collections
import os
from sklearn.model_selection import train_test_split,cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.svm import SVC
import pickle

# load HLA-A*0201 sequence
seq1 = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA'       
seq2 = 'GSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQ'
seq3 = 'RTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'

# load pos and neg dataset
data = pd.read_excel('/Users/ligk2e/Desktop/immunogenecity/A0201.xlsx',index=None)

# pre-process the data
col = []
for i in range(data.shape[0]):
    peptide = data['peptide'].iloc[i]
    peptide_new = peptide.rstrip(' ')
    peptide_new2 = peptide.replace(' ','')
    col.append(peptide_new2)
new_df = pd.DataFrame({'peptide':col})
data.update(new_df)

col = []
for i in range(data.shape[0]):
    peptide = data['peptide'].iloc[i]
    if len(peptide) in set([8,9,10,11]): col.append(True)
    else: col.append(False)
data['col'] = col
data = data[data['col']]
data = data.drop(columns=['col'])
data = data.set_index(pd.Index(np.arange(data.shape[0])))
        



# pre-process a0201.json file
with open('/Users/ligk2e/Desktop/immunogenecity/scrapy/imgt/imgt/a0201.json') as f:
    content = json.load(f)   # using load instead of loads, then content will be a big list

dic1,dic2 = {},{}    
for dic in content:
    for key,value in dic.items():
        if key.startswith('alpha1_interaction') and value != "This ligand is not from 8-11 mer":
            for i,j in value.items():
                i = int(i)
                try:
                    dic1[i].extend(j)
                except KeyError:
                    dic1[i] = []
                    dic1[i].extend(j)
        elif key.startswith('alpha2_interaction') and value != "This ligand is not from 8-11 mer":
            for i,j in value.items():
                i = int(i)
                try:
                    dic2[i].extend(j)
                except KeyError:
                    dic2[i] = []
                    dic2[i].extend(j)

            
for key,value in dic1.items():
    dic1[key] = collections.Counter(value) 

for key,value in dic2.items():
    dic2[key] = collections.Counter(value) 


# load 40 AACP
index = np.zeros([40,20,20])
annotation = []
with open('/Users/ligk2e/Desktop/immunogenecity/file.txt') as f:
    NR = -1
    for line in f:
        NR += 1 
        name = line.rstrip('\n')
        annotation.append(name.rstrip('_.txt'))
        index[NR,:,:] = np.loadtxt(os.path.join('/Users/ligk2e/Desktop/immunogenecity/AAindex_norm',name), delimiter='\t')                                       
        



# Encoding X,Y
transform = {
        'A':0,
        'R':1,
        'N':2,
        'D':3,
        'C':4,
        'Q':5,
        'E':6,
        'G':7,
        'H':8,
        'I':9,
        'L':10,
        'K':11,
        'M':12,
        'F':13,
        'P':14,
        'S':15,
        'T':16,
        'W':17,
        'Y':18,
        'V':19,
        }
pyramid = {
        8: {1:(1),2:(2),3:(3),4:(4),5:(5),6:(6),7:(7),8:(8)},
        9: {1:(1,2),2:(2,3),3:(3,4),4:(4,5),5:(5,6),6:(6,7),7:(7,8),8:(8,9)},
        10: {1:(1,2,3),2:(2,3,4),3:(3,4,5),4:(4,5,6),5:(5,6,7),6:(6,7,8),7:(7,8,9),8:(8,9,10)},
        11: {1:(1,2,3,4),2:(2,3,4,5),3:(3,4,5,6),4:(4,5,6,7),5:(5,6,7,8),6:(6,7,8,9),7:(7,8,9,10),8:(8,9,10,11)}
        }


def calculation(aa,m,k,dic1,dic2,seq1,seq2):  # kth type of interaction, mth relative position
    goal = 0
    total = 0
    inter = index[k,:,:]
    length = len(aa)

    positions = pyramid[length][m]  # the absolute positions that we will consider for each relative position  
    for i in positions:
        query1 = transform[aa[i-1]]
        for pos,counter in dic1[m].items():
            total += counter
            query2 = transform[seq1[pos-1]]
            value = inter[query1,query2]
            weighted_value = value * counter
            goal += weighted_value
        for pos,counter in dic2[m].items():

            try:
                query2 = transform[seq2[pos-1]]
            except:
                continue
            else:
                total += counter
                value = inter[query1,query2]
                weighted_value = value * counter
                goal += weighted_value 
    final = goal/(total*len(positions))
    return final
            
def calculation1(aa,m,k,dic1,dic2,seq1,seq2):  # don't take weight for frequent occured hla position
    goal = 0
    total = 0
    inter = index[k,:,:]
    length = len(aa)

    positions = pyramid[length][m]  # the absolute positions that we will consider for each relative position  
    for i in positions:
        query1 = transform[aa[i-1]]
        for pos,counter in dic1[m].items():
            total += counter
            query2 = transform[seq1[pos-1]]
            value = inter[query1,query2]
            weighted_value = value * counter
            goal += weighted_value
        for pos,counter in dic2[m].items():

            try:
                query2 = transform[seq2[pos-1]]
            except:
                continue
            else:
                total += counter
                value = inter[query1,query2]
                #weighted_value = value * counter
                goal += value 
    final = goal/(len(positions)*(len(dic1[m])+len(dic2[m])))
    return final       

def calculation2(aa,m,k,dic1,dic2,seq1,seq2):  # don't take weight for frequent occured hla position
                                                # and give P2, P3, P4 extra weight
    goal = 0
    total = 0
    inter = index[k,:,:]
    length = len(aa)
    extra_dic = {1:1,2:1.4,3:1.2,4:1.6,5:1,6:1,7:1,8:1}
    positions = pyramid[length][m]  # the absolute positions that we will consider for each relative position  
    for i in positions:
        query1 = transform[aa[i-1]]
        for pos,counter in dic1[m].items():
            total += counter
            query2 = transform[seq1[pos-1]]
            value = inter[query1,query2]
            weighted_value = value * counter
            goal += weighted_value
        for pos,counter in dic2[m].items():

            try:
                query2 = transform[seq2[pos-1]]
            except:
                continue
            else:
                total += counter
                value = inter[query1,query2]
                #weighted_value = value * counter
                goal += value 
    final = extra_dic[m]*goal/(len(positions)*(len(dic1[m])+len(dic2[m])))
    return final  

X = np.zeros((data.shape[0],8*40))
for i in range(X.shape[0]):
    peptide = data.iloc[i]['peptide']
    for j in range(X.shape[1]):   
        X[i,j]=calculation2(peptide,j%8+1,j//8,dic1,dic2,seq1,seq2)

Y = data['cond'].values





# model part1: one-hotting encoding, three layer GRU + attention

# concatenate peptide with seq1 and seq2

from torch.utils.data import Dataset, DataLoader, random_split
import torch
from torch.nn.utils.rnn import pad_sequence

torch.manual_seed(0)

whole = []
for i in range(data.shape[0]):
    peptide = data.iloc[i]['peptide']
    cat = seq1 + peptide + seq2
    whole.append(cat)
data_new = data.join(pd.Series(whole,name='whole'))
    

class my_dataset(Dataset):
    def __init__(self,original):
        self.original = original
        self.new_data = torch.transpose(self.padded_to_max_len(),1,2)   # it is a torch
        self.label = self.get_Y()   # torch [1,209]

        
    def __len__(self):
        return self.original.shape[0]
    
    def __getitem__(self,idx):

        
        return (self.new_data[idx],self.label[0,idx])
    
    def padded_to_max_len(self):
        tmp = []
        for i in range(self.original.shape[0]):
            whole = self.original.iloc[i]['whole']
            #print(whole)
            onehot = torch.from_numpy(my_dataset_X.onehot_peptide(whole)) # length*20
            #print(onehot,onehot.size(),onehot[0,:])
            tmp.append(onehot)
        #print(tmp[0][0,:])
        stacked_equal_length = pad_sequence(tmp,batch_first=True)  #[number of item,length,20]
        
        return stacked_equal_length
    
    def get_Y(self):
        return torch.from_numpy(self.original['cond'].values).view(1,-1).float()
            
        
    
    @staticmethod
    def onehot_peptide(pep):
        result = np.zeros([20,len(pep)])
        for i in range(result.shape[1]):
            result[:,i] = my_dataset_X.onehot_aa(pep[i])
        return np.transpose(result)
    
    @staticmethod
    def onehot_aa(a):
        transform = {
        'A':0,
        'R':1,
        'N':2,
        'D':3,
        'C':4,
        'Q':5,
        'E':6,
        'G':7,
        'H':8,
        'I':9,
        'L':10,
        'K':11,
        'M':12,
        'F':13,
        'P':14,
        'S':15,
        'T':16,
        'W':17,
        'Y':18,
        'V':19,
        }
        mat = torch.eye(len(transform))
        return mat[transform[a],:]
        
        
data_torch = my_dataset(data_new)


#training_set, validation_set, testing_set = random_split(data_torch,(200,5,4))
training_set, testing_set = random_split(data_torch,(170,39))
        

data_torch_loader = DataLoader(data_torch,batch_size=10,shuffle=True,drop_last=True)
for i,batch in enumerate(data_torch_loader):
    print(i,batch)
    break


# let's build the model





































