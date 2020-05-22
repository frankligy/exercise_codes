#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 19:33:04 2020

@author: ligk2e
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
import csv
import numpy as np
import random

base='ATCG'
baseRNA = 'AUCG'


device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')  # device object

def seqtopad(sequence,motlen):
    rows=len(sequence)+2*motlen-2
    S=np.empty([rows,4])
    base= baseRNA
    for i in range(rows):
        for j in range(4):
            if i-motlen+1<len(sequence) and sequence[i-motlen+1]=='N' or i<motlen-1 or i>len(sequence)+motlen-2:
                S[i,j]=np.float32(0.25)
            elif sequence[i-motlen+1]==base[j]:
                S[i,j]=np.float32(1)
            else:
                S[i,j]=np.float32(0)
    return np.transpose(S)


class Spliter():   
    def __init__(self,good,random_=False):   # good is seperable based on some index
        if random_==False: 
            self.good = good
        else: 
            random.shuffle(good)
            self.good=good            
        self.queue = []  # to store splited data
    
    def split_df(self,n):
        dim  = self.good.shape[0]
        lis = [i for i in range(dim)]
        size = len(lis)//n + 1
        for j in range(0,len(lis),size):
            part_index = lis[j:j+size]
            part_df = self.good.iloc[part_index]
            
            self.queue.append(part_df)
            
    def split_ndarray(self,n):
        dim = len(self.good)
        lis = [i for i in range(dim)]
        size = len(lis)//n + 1
        for j in range(0,len(lis),size):
            part_ndarray = self.good[j:j+size]
            self.queue.append(part_ndarray)
            
        
class clipseq_dataset(Dataset):

    def __init__(self,xy=None):
        self.x_data=np.asarray([el[0] for el in xy],dtype=np.float32)
        self.y_data =np.asarray([el[1] for el in xy ],dtype=np.float32)
        self.x_data = torch.from_numpy(self.x_data)
        self.y_data = torch.from_numpy(self.y_data)
        self.len=len(self.x_data)
      

    def __getitem__(self, index):
        return self.x_data[index], self.y_data[index]

    def __len__(self):
        return self.len




# read file
train_set = []
with open('/Users/ligk2e/Downloads/clip_test','r') as f1:
    next(f1)
    data=csv.reader(f1,delimiter='\t')
    for row in data:
        train_set.append([seqtopad(row[0],24),row[1]])
SpliterRAM = Spliter(train_set,True)
SpliterRAM.split_ndarray(3)

firstvalid=SpliterRAM.queue[0]
secondvalid=SpliterRAM.queue[1]
thirdvalid=SpliterRAM.queue[2]
firsttrain=secondvalid+thirdvalid
secondtrain=firstvalid+thirdvalid
thirdtrain=firstvalid+secondvalid

train1_dataset=clipseq_dataset(firsttrain)
train2_dataset=clipseq_dataset(secondtrain)
train3_dataset=clipseq_dataset(thirdtrain)
valid1_dataset=clipseq_dataset(firstvalid)
valid2_dataset=clipseq_dataset(secondvalid)
valid3_dataset=clipseq_dataset(thirdvalid)

train_loader1 = DataLoader(dataset=train1_dataset,batch_size=128,shuffle=True)
train_loader2 = DataLoader(dataset=train2_dataset,batch_size=128,shuffle=True)
train_loader3 = DataLoader(dataset=train3_dataset,batch_size=128,shuffle=True)
valid1_loader = DataLoader(dataset=valid1_dataset,batch_size=128,shuffle=True)
valid2_loader = DataLoader(dataset=valid2_dataset,batch_size=128,shuffle=True)
valid3_loader = DataLoader(dataset=valid3_dataset,batch_size=128,shuffle=True)

for i, (data, target) in enumerate(train_loader1):
    print(type(data))  # tensor
    print(data.size())  # [128,4,147]
    b = data
    a = data.to(device)
    print(type(a))
    break

for i, (data, target) in enumerate(train_loader1):
    print(i)


'''
data is the entity that will feed into the model;

calibration: 40 different combination of hyperparameters, then choose optimal learning step(how many batch as a learning phase), 
for each condition combination, 3-cross validation, get average ROC score. Best combination will be used for training setting.

When training, also training 5 times, take best model as final model.

'''




















        