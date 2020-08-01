#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 10:01:40 2020

@author: ligk2e
"""

import pandas as pd
import xlrd
import numpy as np
import json
import collections
import os
import time
import random
import itertools
from Bio.SubsMat import MatrixInfo

from torch.utils.data import Dataset, DataLoader, random_split,Subset
import torch
from torch.nn.utils.rnn import pad_sequence
import torch.nn as nn
import torch.nn.functional as F

from prettytable import PrettyTable

def clean_series(series):  # give a pandas series
    
    if series.dtype == object:  # pandas will store str as object since string has variable length, you can use astype('|S')
        clean = []
        for item in series:
            item = item.lstrip(' ')   # remove leading whitespace
            item = item.rstrip(' ')   # remove trailing whitespace
            item = item.replace(' ','')  # replace all whitespace in the middle
            clean.append(item)
    else:
        clean = series
        
    
        
    return pd.Series(clean)


def clean_data_frame(data):  # give a pandas dataFrame
    
    peptide_clean = clean_series(data['peptide'])
    hla_clean = clean_series(data['HLA'])
    immunogenecity_clean = clean_series(data['immunogenecity'])
    
    data_clean = pd.concat([peptide_clean,hla_clean,immunogenecity_clean],axis=1)
    data_clean.columns = ['peptide','HLA','immunogenecity']
    
    return data_clean
    


    
def clean_hla(data_clean):
    cond = []
    hla = data_clean['HLA']
    for i in hla:
        if not ':' in i: cond.append(False)
        else: cond.append(True)
    data_clean = data_clean.join(pd.Series(cond,name='cond'))
    data_clean2 = data_clean.loc[data_clean['cond']]
    data_clean2 = data_clean2.drop(columns=['cond'])
    data_clean2 = data_clean2.set_index(pd.Index(np.arange(data_clean2.shape[0])))
    return data_clean2

def convert_hla(data_clean2):
    hla = data_clean2['HLA']
    new = []
    for i in hla:
        if not '*' in i: 
            new.append(i[0:5]+'*'+i[5:])
        else:
            new.append(i)
    data_clean2.update({'HLA':new})
    return data_clean2
            

        
def read_hla(path):
    dic = {}
    with open(path,'r') as f:
        flag = False
        hla = ''
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                hla = line.split(',')[0].lstrip('>')     
                if dic.get(hla) == None:
                    dic[hla] = []
                    flag = True
                else: 
                    flag = True
                
            else:
                if flag == True:
                    dic[hla].append(line)
                    flag = False
    return dic
    #{'hla':[seq1,seq2],}
    
    
    
def get_imgt_contact_site(path,hla,hla_seq):
    hla_short = ''.join(hla[4:].split(':')).replace('*','')
    # pre-process json file
    with open(os.path.join(path,'{0}.json'.format(hla_short))) as f:
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
        
        
    # let's extract all possible contact site
    seq1 = hla_seq[hla][0]
    seq2 = hla_seq[hla][1]
    
    seq1_index = []
    seq2_index = []
    seq1_contact = []
    seq2_contact = []
    for counter in dic1.values():
        counter = sort_by_value(counter,reverse=True)
        for k,v in counter.items():
            if v > 0: seq1_index.append(k)
    
    for i in sorted(list(set(seq1_index))):
        try:
            seq1_contact.append(seq1[i])
        except:
            continue
    
    for counter in dic2.values():
        counter = sort_by_value(counter,reverse=True)
        for k,v in counter.items():
            if v > 0: seq2_index.append(k)  
   
    for i in sorted(list(set(seq2_index))):
        try:
            seq2_contact.append(seq2[i])
        except:
            continue
            
    return seq1_contact,seq2_contact             
              


def sort_by_value(dic,reverse,num=None):
    zipped= zip(dic.keys(),dic.values())
    zipped_list=[i for i in zipped]
    sorted_list = sorted(zipped_list,key=lambda x:x[1],reverse=reverse)
    
    if num == None: num = len(dic)

    
    sorted_list_shortlist = sorted_list[:num]

    back_to_dic = {i[0]:i[1] for i in sorted_list_shortlist}
    return back_to_dic
    
    
def get_whole_seq(path,data_clean):
    whole = []
    for i in range(data_clean.shape[0]):
        print(i)
        peptide = data_clean.iloc[i]['peptide']
        hla = data_clean.iloc[i]['HLA']
        try:
            seq1_contact,seq2_contact= get_imgt_contact_site(path,hla,hla_seq)
        except:
            hla = rescue_unknown_hla(hla,dic_inventory)
            seq1_contact,seq2_contact= get_imgt_contact_site(path,hla,hla_seq)
            
        cat = ''.join(seq1_contact) + peptide + ''.join(seq2_contact)

        whole.append(cat)
    data_whole = data_clean.join(pd.Series(whole,name='whole'))
    return data_whole

        

def dict_inventory(inventory):
    dicA,dicB,dicC = {},{},{}
    dic = {'A':dicA,'B':dicB,'C':dicC}
    
    for hla in inventory:
        type_ = hla[0]  # A,B,C
        first2 = hla[1:3] # 01
        last2 = hla[3:5]  # 01
        try:
            dic[type_][first2].append(last2)
        except KeyError:
            dic[type_][first2] = []
            dic[type_][first2].append(last2)
            
    return dic

def rescue_unknown_hla(hla,dic_inventory):
    print(hla)
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[9:11]
    big_category = dic_inventory[type_]
    if not big_category.get(first2) == None:
        small_category = big_category.get(first2)
        distance = [abs(int(last2)-int(i)) for i in small_category]
        optimal = min(zip(small_category,distance),key=lambda x:x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(first2) + ':' + str(optimal)
    else:
        small_category = list(big_category.keys())
        distance = [abs(int(first2)-int(i)) for i in small_category]   
        optimal = min(zip(small_category,distance),key=lambda x:x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(optimal) + ':' + str(big_category[optimal][0])
    

class MyDataSet(Dataset):
    def __init__(self,original):
        self.original = original
        self.check_validity()
        self.new_data = self.encode_peptide()
        self.y = self.get_y()
        
    def __len__(self):
        return self.y.shape[1]
    
    def __getitem__(self,idx):
        X = self.new_data[idx,:,:]
        y = self.y[0,idx]
        return (X,y)
    
    def encode_peptide(self):
        whole = self.original['whole']
        for_padding = []
        for i in range(whole.shape[0]):
            print(i)
            length = len(whole[i])
            result = np.zeros([length,20])
            for j in range(len(whole[i])):   # iterate each amino acid
                result[j,:] = MyDataSet.blosum50_aa(whole[i][j])
            for_padding.append(torch.from_numpy(result).float())
        # Here, we should have a list for_padding that [torch(length1*20),torch(length2*20),...]

        
        
        # padded = pad_sequence(for_padding,batch_first=True,padding_value= -1.0)
        padded = self.padding(for_padding)
        
        
        
        # padded should be torch(40000,max_length,20)
        return padded
    
    
    def get_y(self):
        y = torch.from_numpy(self.original['immunogenecity'].values).view(1,-1).float()  # torch(1,40000)
        return y
    

    def padding(self,lis):   # accept for_padding
        
        # find max_length
        max_length=max([torch.shape[0] for torch in lis])
        
        # padding
        bucket = []
        for item in lis:
            length = item.shape[0]
            gap = max_length - length
            if gap % 2 == 0:  # even number
                gapped_left, gapped_right = gap // 2, gap //2  # will be an int
            else:  # odd number
                if np.random.uniform() < 0.5:  # randomly decide which side will have one more padded value
                    gapped_left = gap // 2
                    gapped_right = gap - gapped_left
                else:
                    gapped_right = gap // 2
                    gapped_left = gap - gapped_right
                    
            padding_left = torch.empty([gapped_left,20]).fill_(-1.0)
            padding_right = torch.empty([gapped_right,20]).fill_(-1.0)
            final = torch.cat([padding_left,item,padding_right],dim=0)
            bucket.append(final)  
        
        # here we should have a list bucket that [torch(max_length*20),torch(max_length*20),...]
        padded_version = torch.stack(bucket,dim=0)
        
        self.max_length = max_length
        return padded_version
            
                
        
    def check_validity(self):
        cond = []
        for i in range(self.original.shape[0]):
            whole = self.original['whole'].iloc[i]
            if 'X' in set(list(whole)): cond.append(False)
            else: cond.append(True)
        self.original=self.original.join(pd.Series(cond,name='cond'))
        self.original=self.original.loc[self.original['cond']]
        self.original=self.original.drop(columns=['cond'])
        self.original=self.original.set_index(pd.Index(np.arange(self.original.shape[0])))    
    
    @staticmethod
    def blosum50_aa(a):
        a = a.upper()
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
        transform_rev = {
        0:'A',
        1:'R',
        2:'N',
        3:'D',
        4:'C',
        5:'Q',
        6:'E',
        7:'G',
        8:'H',
        9:'I',
        10:'L',
        11:'K',
        12:'M',
        13:'F',
        14:'P',
        15:'S',
        16:'T',
        17:'W',
        18:'Y',
        19:'V',
        }
        dic = MatrixInfo.blosum50

        matrix = np.zeros([20,20])
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                try:
                    matrix[i,j] = dic[(transform_rev[i],transform_rev[j])] 
                except KeyError:
                    matrix[i,j] = dic[(transform_rev[j],transform_rev[i])]
                
                
        return matrix[:,transform[a]]
        

def calculate_conv2d_dimension(dim_in,padding,dilation,kernel,stride):
    dim_out = (dim_in + 2*padding - dilation*(kernel-1) - 1) / stride + 1
    return int(dim_out)  
                

class dilated_CNN(nn.Module):
    def __init__(self,length,channel=1,filters=16,kernel=(10,20),stride=(1,1),padding=(0,0),dilation=(2,1),hidden=100):
        super(dilated_CNN,self).__init__()
        self.length = length
        self.channel = channel
        self.filters = filters
        self.kernel = kernel
        self.stride = stride
        self.padding = padding
        self.dilation = dilation
        self.hidden = hidden
        
        self.layer1 = nn.Sequential(
            nn.Conv2d(self.channel, self.filters, self.kernel, dilation=self.dilation),
            nn.BatchNorm2d(self.filters),
            nn.ReLU(),
            nn.MaxPool2d((3,1),stride=self.stride)
            )

        
        layer1_conv2d_output_H = calculate_conv2d_dimension(self.length, self.padding[0], self.dilation[0], self.kernel[0], 
                                                            self.stride[0])
        layer1_conv2d_output_W = calculate_conv2d_dimension(20, self.padding[1], self.dilation[1], self.kernel[1], 
                                                            self.stride[1])    
        
        layer1_maxpool2d_output_H = calculate_conv2d_dimension(layer1_conv2d_output_H,padding=0,dilation=1,kernel=3,stride=1)
        layer1_maxpool2d_output_W = calculate_conv2d_dimension(layer1_conv2d_output_W,padding=0,dilation=1,kernel=1,stride=1)      
        
        print((layer1_conv2d_output_H,layer1_conv2d_output_W),(layer1_maxpool2d_output_H,layer1_maxpool2d_output_W))
        
        self.layer2 = nn.Sequential(
            nn.Conv2d(self.filters,self.filters, kernel_size=(10,1),dilation=self.dilation),
            nn.BatchNorm2d(self.filters),
            nn.ReLU(),
            nn.MaxPool2d((3,1),stride=self.stride)
            )
        
        layer2_conv2d_output_H = calculate_conv2d_dimension(layer1_maxpool2d_output_H, self.padding[0], 
                                                            self.dilation[0], 10, self.stride[0])
        layer2_conv2d_output_W = calculate_conv2d_dimension(layer1_maxpool2d_output_W, self.padding[1], 
                                                            self.dilation[1], 1, self.stride[1])
        
        layer2_maxpool2d_output_H = calculate_conv2d_dimension(layer2_conv2d_output_H,padding=0,dilation=1,kernel=3,
                                                               stride=1)
        layer2_maxpool2d_output_W = calculate_conv2d_dimension(layer2_conv2d_output_W,padding=0,dilation=1,kernel=1,
                                                               stride=1)
        
        print((layer2_conv2d_output_H,layer2_conv2d_output_W),(layer2_maxpool2d_output_H,layer2_maxpool2d_output_W))
        
        self.layer3 = nn.Sequential(
            nn.Dropout(p=0.25),
            nn.Linear(self.filters * layer2_maxpool2d_output_H * layer2_maxpool2d_output_W, self.hidden),
            nn.Dropout(p=0.25),
            nn.ReLU(),
            nn.Linear(self.hidden,2)
        )
        
        self.Sigmoid = nn.Sigmoid()
        
    def forward(self,x):   # x: [batch,channel,h=max_length,w=20]
        out = self.layer1(x)

        out = self.layer2(out)
        out = out.view(out.shape[0],-1)

        out = self.layer3(out)
        out = self.Sigmoid(out)
        return out
        

def balancedBinaryLoader(dataset,batch_size):
    
    dic = {'0':[],'1':[]}
    for i in range(len(dataset)):
        X = dataset[i][0]
        y = dataset[i][1]
        if y == 1: dic['1'].append(i)
        elif y == 0: dic['0'].append(i)
        
    #print(dic)
    

    
    sample_size = batch_size // 2  # will be an int, make sure batch_size is an even number
    
    negative = Subset(dataset,dic['0']) # dataset.Subset object
    positive = Subset(dataset,dic['1'])
    # print(len(positive),type(positive)) 
    
    negative_loader = DataLoader(negative,batch_size=sample_size,shuffle=True,drop_last=True)
    positive_loader = DataLoader(positive,batch_size=sample_size,shuffle=True,drop_last=True)
    



    neg_chunks_X = []
    neg_chunks_y = []
    for idx,(X,y) in enumerate(negative_loader):
        neg_chunks_X.append(X)
        neg_chunks_y.append(y)

    
    pos_chunks_X = []
    pos_chunks_y = []
    for idx,(X,y) in enumerate(positive_loader):
        pos_chunks_X.append(X)
        pos_chunks_y.append(y)

    
    pos_chunks_X_cycle = pos_chunks_X * 10
    pos_chunks_y_cycle = pos_chunks_y * 10


    chunks_X_list = []
    chunks_y_list = []    
    for i in range(len(neg_chunks_X)):
        chunks_X = torch.cat([neg_chunks_X[i],pos_chunks_X_cycle[i]],dim=0)
        chunks_y = torch.cat([neg_chunks_y[i],pos_chunks_y_cycle[i]],dim=0)
        chunks_X_list.append(chunks_X)
        chunks_y_list.append(chunks_y)
        
    
        

    
    loader = list(zip(chunks_X_list,chunks_y_list)) # zip can only be iterated once
    return loader    

if __name__ == '__main__':
    
    start_time = time.time()
    
    data = pd.read_excel('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data.xlsx')
    data_clean = clean_data_frame(data)
    data_clean = clean_hla(data_clean)
    data_clean = convert_hla(data_clean)
    
    hla_seq = read_hla('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_seq/hla_seq.txt')
    inventory = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/imgt_contact/inventory.txt',sep='\n',header=None)[0].tolist() 
    dic_inventory = dict_inventory(inventory)   
    data_whole = get_whole_seq('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/imgt_contact',data_clean)
    
    
    data_torch = MyDataSet(data_whole)

    #training_set, validation_set, testing_set = random_split(data_torch,(200,5,4))
    training_set, testing_set = random_split(data_torch,(40000,3519))
            
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Starting to work
    # data_torch_training_loader = DataLoader(training_set,batch_size=512,shuffle=True,drop_last=True)
    # data_torch_testing_loader = DataLoader(testing_set,batch_size=512,shuffle=True,drop_last=True)
    
    data_torch_training_loader = balancedBinaryLoader(training_set,batch_size=64)
    
    
    model = dilated_CNN(data_torch.max_length)
    
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
    scheduler = scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, factor=0.1, patience=2, verbose=True)
    # if it observe a non-decreasing loss, give you another (patience-1) more chances, if still not decrease, will reduce 
    # learning rate to factor*learning rate
    loss_f=nn.CrossEntropyLoss()

    
    num_epochs = 5
    for epoch in range(num_epochs):
        loss_list = []
        acc_list = []
        for i in data_torch_training_loader:


            X = i[0].unsqueeze(1).float().to(device)
            y = i[1].long().to(device)
            #print(X.size(),y,size())
            optimizer.zero_grad()
            
            y_pred = model(X)
            loss = loss_f(y_pred,y)
            loss.backward()
            optimizer.step()
            loss_list.append(loss.item())
            
            num_correct = 0
            num_samples = 0
            _,predictions = y_pred.max(1)

            num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it

            num_samples  += predictions.size(0)

            acc_list.append(float(num_correct)/float(num_samples)*100)
            
        loss,acc = sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    
        scheduler.step(loss)
        print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}%'.format(epoch+1,num_epochs,loss,acc))
        
        
                
        
                    
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    